//! What a small request may cost: the bounds that hold when the work is in the records and not in
//! their number.
//!
//! The limits of the server are written in terms of variants, and a variant is not a unit of cost.
//! A structural variant makes a row for every transcript it overlaps, so a body of 7 KB, a hundred
//! deletions of a chromosome arm, held 2.4 GB on a genome of 6,000 transcripts before the first byte
//! of the answer went out; a record with 200,000 alternate alleles took 108 s and was answered
//! after the time limit; a position of 2^64 - 1 overflowed. An independent review found each of them,
//! and these tests keep them found.

use axum::body::Body;
use axum::http::{Request, StatusCode};
use fastvep_web::context::AnnotationContext;
use fastvep_web::handlers::{AppState, ServerConfig, SharedState};
use fastvep_web::limits::Limits;
use serde_json::{json, Value};
use std::sync::Arc;
use std::time::{Duration, Instant};
use tower::ServiceExt;

/// `genes` one-transcript genes of 3 kb on chromosome 17, 10 kb apart.
fn synthetic_gff3(genes: usize) -> String {
    let mut s = String::from("##gff-version 3\n");
    for i in 0..genes {
        let start = 10_000 + i as u64 * 10_000;
        let end = start + 3_000;
        let (id, tx) = (format!("G{i}"), format!("T{i}"));
        s.push_str(&format!(
            "17\ttest\tgene\t{start}\t{end}\t.\t+\t.\tID=gene:{id};Name={id};biotype=protein_coding;gene_id={id}\n\
             17\ttest\tmRNA\t{start}\t{end}\t.\t+\t.\tID=transcript:{tx};Parent=gene:{id};Name={id}-201;biotype=protein_coding;tag=Ensembl_canonical;transcript_id={tx}\n\
             17\ttest\texon\t{start}\t{e1}\t.\t+\t.\tID=exon:{tx}-E1;Parent=transcript:{tx};rank=1\n\
             17\ttest\tCDS\t{start}\t{e1}\t.\t+\t0\tID=CDS:{tx};Parent=transcript:{tx};protein_id=P{i}\n\
             17\ttest\texon\t{s2}\t{end}\t.\t+\t.\tID=exon:{tx}-E2;Parent=transcript:{tx};rank=2\n\
             17\ttest\tCDS\t{s2}\t{end}\t.\t+\t0\tID=CDS:{tx};Parent=transcript:{tx};protein_id=P{i}\n",
            e1 = start + 200,
            s2 = end - 200,
        ));
    }
    s
}

fn state(genes: usize, limits: Limits) -> AppState {
    let mut ctx = AnnotationContext::new(None, None, None, 5000).expect("context");
    ctx.update_gff3_text(&synthetic_gff3(genes)).expect("gff3");
    Arc::new(SharedState::new(
        ctx,
        "synthetic",
        ServerConfig {
            data_dir: None,
            sa_dir: None,
            stats_file: None,
            distance: 5000,
            allow_model_replacement: false,
            acmg_allow_partial: false,
            limits,
        },
    ))
}

async fn post(state: &AppState, body: Value) -> (StatusCode, Value) {
    let router = fastvep_web::build_router(Arc::clone(state), 10_485_760, 8);
    let response = router
        .oneshot(
            Request::builder()
                .method("POST")
                .uri("/api/annotate")
                .header("content-type", "application/json")
                .body(Body::from(body.to_string()))
                .unwrap(),
        )
        .await
        .unwrap();
    let status = response.status();
    let bytes = axum::body::to_bytes(response.into_body(), usize::MAX)
        .await
        .unwrap();
    (
        status,
        serde_json::from_slice(&bytes).unwrap_or(Value::Null),
    )
}

/// A deletion from the first base of the chromosome to `end`, which overlaps every transcript.
fn deletion(i: usize, end: u64) -> String {
    format!(
        "17\t{}\tdel{i}\tN\t<DEL>\t50\tPASS\tSVTYPE=DEL;END={end};SVLEN=-{}",
        1 + i as u64,
        end - 1
    )
}

fn error_of(body: &Value) -> String {
    body["error"].as_str().unwrap_or_default().to_string()
}

#[tokio::test]
async fn variants_that_overlap_many_transcripts_are_annotated_a_few_at_a_time() {
    // 100 transcripts, so a deletion over the chromosome is 100 rows, and a limit under what two of
    // them make: the request is refused after the first batch that passes it, and a batch is one
    // variant. With the hundred variants of the old batches it would have been annotated whole,
    // all ten, and refused after that.
    let end = 10_000 + 100 * 10_000;
    let st = state(
        100,
        Limits {
            batch_rows: 150,
            max_response_bytes: 150_000,
            ..Limits::default()
        },
    );
    let vcf: Vec<String> = (0..10).map(|i| deletion(i, end)).collect();
    let (status, body) = post(&st, json!({ "vcf": vcf.join("\n") })).await;
    assert_eq!(status, StatusCode::PAYLOAD_TOO_LARGE, "{body}");
    let message = error_of(&body);
    assert!(message.contains("of 10 variants were done"), "{message}");
    assert!(!message.contains("(10 of 10"), "{message}");
    let done: usize = message
        .split("(")
        .nth(1)
        .and_then(|t| t.split(' ').next())
        .and_then(|n| n.parse().ok())
        .unwrap_or(99);
    // A variant of this genome is about 40 KB of answer and the limit is 147 KiB, so the fourth
    // variant is the one that passes it. All ten is the old batch.
    assert!(
        (2..=6).contains(&done),
        "{done} variants were annotated before the refusal: {message}"
    );
}

#[tokio::test]
async fn a_variant_whose_answer_alone_cannot_fit_is_refused_before_it_is_annotated() {
    let end = 10_000 + 400 * 10_000;
    let st = state(
        400,
        Limits {
            max_response_bytes: 100_000,
            ..Limits::default()
        },
    );
    let started = Instant::now();
    let (status, body) = post(&st, json!({ "vcf": deletion(0, end) })).await;
    assert_eq!(status, StatusCode::PAYLOAD_TOO_LARGE, "{body}");
    let message = error_of(&body);
    assert!(message.contains("VCF line 1 overlaps 4"), "{message}");
    assert!(message.contains("`pick`"), "{message}");
    assert!(
        started.elapsed() < Duration::from_secs(2),
        "refused after {:?}: it was annotated",
        started.elapsed()
    );
    // The answer the client is told to ask for fits, and is made.
    let (status, body) = post(&st, json!({ "vcf": deletion(0, end), "pick": true })).await;
    assert_eq!(status, StatusCode::OK, "{body}");
    assert_eq!(body["count"], 1);
}

#[tokio::test]
async fn a_record_with_more_alternate_alleles_than_the_server_annotates_is_refused_with_its_line() {
    let st = state(2, Limits::default());
    let alts = |n: usize| {
        (0..n)
            .map(|i| ["A", "C", "G", "T"][i % 4].repeat(1 + i / 4))
            .collect::<Vec<_>>()
            .join(",")
    };
    // The reference base of 17:10100 is not known (no FASTA), so any ALT is a different allele; the
    // alleles of the list are distinct.
    let line = |n: usize| format!("17\t10100\t.\tC\t{}\t.\t.\t.", alts(n).replace("C,", "CC,"));
    let (status, body) = post(&st, json!({ "vcf": format!("#header\n{}", line(257)) })).await;
    assert_eq!(status, StatusCode::BAD_REQUEST, "{body}");
    let message = error_of(&body);
    assert!(
        message.contains("VCF line 2 has 257 alternate alleles") && message.contains("256"),
        "{message}"
    );
    let (status, body) = post(&st, json!({ "vcf": line(256) })).await;
    assert_eq!(status, StatusCode::OK, "{body}");
}

#[tokio::test]
async fn a_position_beyond_the_range_of_a_coordinate_is_refused_not_wrapped() {
    let st = state(2, Limits::default());
    for line in [
        "17\t18446744073709551615\t.\tA\tG\t.\t.\t.".to_string(),
        "17\t4294967296\t.\tA\tG\t.\t.\t.".to_string(),
        "17\t20000\t.\tN\t<DEL>\t.\t.\tSVTYPE=DEL;END=18446744073709551615".to_string(),
    ] {
        let (status, body) = post(&st, json!({ "vcf": line })).await;
        assert_eq!(status, StatusCode::BAD_REQUEST, "{line}: {body}");
        assert!(error_of(&body).contains("VCF line 1 names a position or END beyond 4294967295"));
    }
    // The largest coordinate there is room for is annotated.
    let (status, body) = post(&st, json!({ "vcf": "17\t4294967295\t.\tA\tG\t.\t.\t." })).await;
    assert_eq!(status, StatusCode::OK, "{body}");
}

#[tokio::test]
async fn the_time_an_annotation_is_allowed_runs_from_its_work_and_not_from_the_wait_for_a_place() {
    // One place, taken for 700 ms, and a request that is allowed 500 ms of work: it waits, is given
    // the place, and does a variant's work in far less than 500 ms. Counted from the request, as it
    // was, the 700 ms of waiting were spent before the first batch and the answer was a 413.
    let st = state(
        2,
        Limits {
            max_parallel: 1,
            queue_timeout: Duration::from_secs(5),
            max_request_time: Duration::from_millis(500),
            ..Limits::default()
        },
    );
    let place = Arc::clone(&st.compute).acquire_owned().await.unwrap();
    let release = tokio::spawn(async move {
        tokio::time::sleep(Duration::from_millis(700)).await;
        drop(place);
    });
    let (status, body) = post(&st, json!({ "vcf": "17\t10100\t.\tA\tG\t.\t.\t." })).await;
    release.await.unwrap();
    assert_eq!(status, StatusCode::OK, "{body}");
    assert_eq!(body["count"], 1);
}
