//! Contract tests for the JSON API documented in `docs/API.md`.
//!
//! These drive the real router through `oneshot`, so they cover routing, the
//! middleware stack, and serialization the way a client meets them. The point
//! is that the request and response shapes people write clients against are
//! asserted somewhere, not just described in a document that can drift.

use axum::body::Body;
use axum::http::{Request, StatusCode};
use fastvep_web::context::AnnotationContext;
use fastvep_web::handlers::{AppState, SharedState};
use serde_json::{json, Value};
use std::sync::atomic::AtomicU64;
use std::sync::{Arc, RwLock};
use tower::ServiceExt;

const MINI_GFF3: &str = include_str!("../fixtures/mini.gff3");

/// One transcript, no FASTA, no supplementary annotations. Deliberately the
/// minimum a caller can stand up, so a failure here is the API's fault rather
/// than the fixture's.
fn test_state() -> AppState {
    let mut ctx = AnnotationContext::new(None, None, None, 5000).expect("build context");
    ctx.update_gff3_text(MINI_GFF3).expect("load mini gff3");
    Arc::new(SharedState {
        ctx: RwLock::new(ctx),
        data_dir: None,
        sa_dir: None,
        // No stats file: these tests must not write to the working directory.
        stats_file: None,
        total_variants: AtomicU64::new(0),
        total_genomes: AtomicU64::new(0),
    })
}

fn router(state: &AppState) -> axum::Router {
    fastvep_web::build_router(Arc::clone(state), 10_485_760, 8)
}

async fn read_json(resp: axum::response::Response) -> (StatusCode, Value) {
    let status = resp.status();
    let bytes = axum::body::to_bytes(resp.into_body(), usize::MAX)
        .await
        .expect("read body");
    let value = serde_json::from_slice(&bytes).expect("body is JSON");
    (status, value)
}

async fn get(state: &AppState, uri: &str) -> (StatusCode, Value) {
    let resp = router(state)
        .oneshot(Request::builder().uri(uri).body(Body::empty()).unwrap())
        .await
        .unwrap();
    read_json(resp).await
}

async fn post_json(state: &AppState, uri: &str, body: Value) -> (StatusCode, Value) {
    let resp = router(state)
        .oneshot(
            Request::builder()
                .method("POST")
                .uri(uri)
                .header("content-type", "application/json")
                .body(Body::from(body.to_string()))
                .unwrap(),
        )
        .await
        .unwrap();
    read_json(resp).await
}

fn vcf_line(pos: u32, id: &str) -> String {
    format!("17\t{}\t{}\tG\tA\t50\tPASS\t.", pos, id)
}

#[tokio::test]
async fn status_reports_what_is_loaded() {
    let (status, body) = get(&test_state(), "/api/status").await;

    assert_eq!(status, StatusCode::OK);
    assert_eq!(body["status"], "ok");
    assert_eq!(body["backend"], true);
    assert_eq!(body["version"], env!("CARGO_PKG_VERSION"));
    assert_eq!(body["transcripts"], 1);
    // Callers are told to gate on these two before trusting an annotation, so
    // they have to be present and honest even when nothing is loaded.
    assert_eq!(body["has_fasta"], false);
    assert_eq!(body["sa_sources"], json!([]));
}

#[tokio::test]
async fn annotate_returns_one_result_per_record() {
    let vcf = [
        vcf_line(1100, "."), // exon 1, inside the CDS
        vcf_line(2000, "."), // between the two exons
        vcf_line(500, "."),  // before the gene, within --distance
    ]
    .join("\n");

    let (status, body) = post_json(&test_state(), "/api/annotate", json!({ "vcf": vcf })).await;

    assert_eq!(status, StatusCode::OK);
    assert_eq!(body["count"], 3);
    let results = body["results"].as_array().unwrap();
    assert_eq!(results.len(), 3);

    let terms: Vec<&str> = results
        .iter()
        .map(|r| r["most_severe_consequence"].as_str().unwrap())
        .collect();
    assert_eq!(
        terms,
        [
            "coding_sequence_variant",
            "intron_variant",
            "upstream_gene_variant"
        ]
    );

    // Record order is the caller's order: clients index into `results`
    // positionally against what they sent.
    let starts: Vec<u64> = results
        .iter()
        .map(|r| r["start"].as_u64().unwrap())
        .collect();
    assert_eq!(starts, [1100, 2000, 500]);

    let first = &results[0]["transcript_consequences"][0];
    assert_eq!(first["transcript_id"], "TXA");
    assert_eq!(first["gene_symbol"], "GENEA");
    assert_eq!(first["biotype"], "protein_coding");
    assert_eq!(
        first["consequence_terms"],
        json!(["coding_sequence_variant"])
    );
}

#[tokio::test]
async fn annotate_accepts_a_vcf_header_and_carries_the_id_through() {
    let vcf = format!(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n{}",
        vcf_line(1100, "rs123456")
    );

    let (status, body) = post_json(&test_state(), "/api/annotate", json!({ "vcf": vcf })).await;

    assert_eq!(status, StatusCode::OK);
    assert_eq!(body["count"], 1);
    // Header lines are skipped rather than parsed as records, and the ID
    // column round-trips so callers can match responses to their input.
    assert_eq!(body["results"][0]["id"], "rs123456");
}

#[tokio::test]
async fn coding_variant_without_a_fasta_reports_the_generic_term() {
    // The hazard documented under "Always pass --fasta": with no reference
    // sequence the server cannot read the codon, so a coding change comes back
    // as `coding_sequence_variant` with MODIFIER impact instead of the specific
    // term and its real impact. It does not error, and the response looks
    // perfectly well-formed, which is what makes it worth pinning here.
    let (status, body) = post_json(
        &test_state(),
        "/api/annotate",
        json!({ "vcf": vcf_line(1100, ".") }),
    )
    .await;

    assert_eq!(status, StatusCode::OK);
    assert_eq!(
        body["results"][0]["most_severe_consequence"],
        "coding_sequence_variant"
    );
    assert_eq!(
        body["results"][0]["transcript_consequences"][0]["impact"],
        "MODIFIER"
    );
}

#[tokio::test]
async fn acmg_is_attached_only_when_requested() {
    let state = test_state();
    let vcf = vcf_line(1100, ".");

    let (status, off) = post_json(&state, "/api/annotate", json!({ "vcf": vcf })).await;
    assert_eq!(status, StatusCode::OK);
    let plain = &off["results"][0]["transcript_consequences"][0];
    // Assert the consequence itself is really there first: `get` on a missing
    // path also returns None, which would make the check below pass for the
    // wrong reason if the response shape ever changed underneath it.
    assert!(plain.is_object(), "expected a transcript consequence");
    assert!(
        plain.get("acmg").is_none(),
        "acmg must be absent by default so clients do not pay for it unasked"
    );

    let (status, on) =
        post_json(&state, "/api/annotate", json!({ "vcf": vcf, "acmg": true })).await;
    assert_eq!(status, StatusCode::OK);
    let acmg = &on["results"][0]["transcript_consequences"][0]["acmg"];
    assert!(acmg.is_object(), "acmg object missing when requested");
    assert!(
        acmg["classification"].is_string(),
        "acmg block must carry a classification"
    );
    assert!(
        acmg["criteria"].is_array(),
        "acmg block must carry per-criterion verdicts"
    );
}

#[tokio::test]
async fn annotate_rejects_an_empty_vcf() {
    let state = test_state();

    for body in [json!({ "vcf": "" }), json!({})] {
        let (status, body) = post_json(&state, "/api/annotate", body).await;
        assert_eq!(status, StatusCode::BAD_REQUEST);
        // 400 messages are caller-facing and documented; 500 messages are not.
        assert_eq!(body["error"], "No VCF data provided");
    }
}

#[tokio::test]
async fn oversized_body_is_rejected_rather_than_annotated() {
    // The body cap is the boundary between "small-N queries" and "use the
    // CLI", so it has to actually stop a large request.
    let app = fastvep_web::build_router(test_state(), 64, 8);
    let body = json!({ "vcf": vcf_line(1100, ".").repeat(50) }).to_string();

    let resp = app
        .oneshot(
            Request::builder()
                .method("POST")
                .uri("/api/annotate")
                .header("content-type", "application/json")
                .body(Body::from(body))
                .unwrap(),
        )
        .await
        .unwrap();

    assert_eq!(resp.status(), StatusCode::PAYLOAD_TOO_LARGE);
}

#[tokio::test]
async fn genome_endpoints_are_inert_without_a_data_dir() {
    let state = test_state();

    let (status, body) = get(&state, "/api/genomes").await;
    assert_eq!(status, StatusCode::OK);
    assert_eq!(body["genomes"], json!([]));

    let (status, body) = post_json(&state, "/api/load-genome", json!({ "name": "human" })).await;
    assert_eq!(status, StatusCode::BAD_REQUEST);
    assert_eq!(body["error"], "No data directory configured");
}

#[tokio::test]
async fn upload_gff3_replaces_the_active_gene_model() {
    let state = test_state();

    let resp = router(&state)
        .oneshot(
            Request::builder()
                .method("POST")
                .uri("/api/upload-gff3")
                .body(Body::from(MINI_GFF3))
                .unwrap(),
        )
        .await
        .unwrap();
    let (status, body) = read_json(resp).await;

    assert_eq!(status, StatusCode::OK);
    assert_eq!(body["genes"], 1);
    assert_eq!(body["transcripts"], 1);

    // The swap is visible to the next caller, not just to this request.
    let (_, status_body) = get(&state, "/api/status").await;
    assert_eq!(status_body["transcripts"], 1);
}
