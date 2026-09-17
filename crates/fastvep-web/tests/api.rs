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
const PICK_GFF3: &str = include_str!("../fixtures/pick.gff3");

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
        // The default a server ships with: the gene model is whatever it was
        // started with, and no request can move it (#113).
        allow_model_replacement: false,
        gene_model_generation: AtomicU64::new(0),
    })
}

/// The same, started with `--allow-model-replacement` - the single-user
/// desktop configuration the browser UI's genome switching needs.
fn replaceable_state() -> AppState {
    let state = test_state();
    let mut ctx = AnnotationContext::new(None, None, None, 5000).expect("build context");
    ctx.update_gff3_text(MINI_GFF3).expect("load mini gff3");
    Arc::new(SharedState {
        ctx: RwLock::new(ctx),
        data_dir: state.data_dir.clone(),
        sa_dir: state.sa_dir.clone(),
        stats_file: None,
        total_variants: AtomicU64::new(0),
        total_genomes: AtomicU64::new(0),
        allow_model_replacement: true,
        gene_model_generation: AtomicU64::new(0),
    })
}

/// Two overlapping transcripts of one gene, the first-seen one non-canonical.
/// Kept separate from [`test_state`] so the single-transcript counts the other
/// tests assert stay meaningful.
fn pick_state() -> AppState {
    let mut ctx = AnnotationContext::new(None, None, None, 5000).expect("build context");
    ctx.update_gff3_text(PICK_GFF3).expect("load pick gff3");
    Arc::new(SharedState {
        ctx: RwLock::new(ctx),
        data_dir: None,
        sa_dir: None,
        stats_file: None,
        total_variants: AtomicU64::new(0),
        total_genomes: AtomicU64::new(0),
        allow_model_replacement: false,
        gene_model_generation: AtomicU64::new(0),
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
async fn pick_returns_exactly_one_transcript_and_it_is_the_canonical_one() {
    // Regression test for the drift between the two annotation drivers. The
    // CLI ran VEP's `--pick_order` hierarchy; this path kept a transcript when
    // it was canonical *or* the first one seen, which returned both rows here
    // with the non-canonical TXB first. A client reading
    // `transcript_consequences[0]` therefore got the wrong transcript from a
    // request that had explicitly asked for one answer.
    let state = pick_state();
    let vcf = vcf_line(1100, "pick");

    let (status, all) = post_json(&state, "/api/annotate", json!({ "vcf": &vcf })).await;
    assert_eq!(status, StatusCode::OK);
    let unpicked = all["results"][0]["transcript_consequences"]
        .as_array()
        .expect("transcript consequences");
    // Both transcripts really do overlap the variant, so the pick below is
    // making a choice rather than describing a set that only had one member.
    let ids: Vec<&str> = unpicked
        .iter()
        .map(|t| t["transcript_id"].as_str().unwrap_or_default())
        .collect();
    assert!(ids.contains(&"TXA") && ids.contains(&"TXB"), "got {ids:?}");

    let (status, picked) = post_json(
        &state,
        "/api/annotate",
        json!({ "vcf": &vcf, "pick": true }),
    )
    .await;
    assert_eq!(status, StatusCode::OK);
    let kept = picked["results"][0]["transcript_consequences"]
        .as_array()
        .expect("transcript consequences");
    assert_eq!(kept.len(), 1, "pick must reduce to one row, got {kept:?}");
    assert_eq!(
        kept[0]["transcript_id"], "TXA",
        "pick must keep the canonical transcript, not the first one seen"
    );
}

#[tokio::test]
async fn flag_pick_keeps_every_transcript_and_marks_the_one_pick_would_have_kept() {
    // `pick` and `flag_pick` must agree about the winner, or the flag means
    // something other than the option it is named after. Worth asserting over
    // HTTP and not only in `pick::`: the six switches reach this handler
    // through `#[serde(flatten)]`, which is a runtime contract - a name that
    // does not deserialize fails silently as "no pick requested".
    let state = pick_state();
    let vcf = vcf_line(1100, "flagpick");

    let (status, flagged) = post_json(
        &state,
        "/api/annotate",
        json!({ "vcf": &vcf, "flag_pick": true }),
    )
    .await;
    assert_eq!(status, StatusCode::OK);
    let all = flagged["results"][0]["transcript_consequences"]
        .as_array()
        .expect("transcript consequences");
    assert_eq!(all.len(), 2, "flag_pick retains every transcript: {all:?}");

    let marked: Vec<&str> = all
        .iter()
        .filter(|t| t["pick"].as_u64() == Some(1))
        .map(|t| t["transcript_id"].as_str().unwrap_or_default())
        .collect();
    assert_eq!(
        marked,
        vec!["TXA"],
        "flag_pick must mark the transcript `pick` keeps, in {all:?}"
    );
    // Absent, not 0, on the rest - the same convention as `canonical`.
    assert!(
        all.iter()
            .filter(|t| t["transcript_id"] != "TXA")
            .all(|t| t.get("pick").is_none()),
        "unpicked entries should not carry the key: {all:?}"
    );

    // And a request that asks for nothing says nothing about picking.
    let (status, plain) = post_json(&state, "/api/annotate", json!({ "vcf": &vcf })).await;
    assert_eq!(status, StatusCode::OK);
    assert!(
        plain["results"][0]["transcript_consequences"]
            .as_array()
            .expect("transcript consequences")
            .iter()
            .all(|t| t.get("pick").is_none()),
        "a plain request should not mention pick"
    );
}

#[tokio::test]
async fn every_pick_switch_is_accepted_under_its_documented_name() {
    // docs/API.md lists these six names. A typo in one of them deserializes as
    // `false` and the request is answered without the pick it asked for, which
    // is a wrong answer that looks like a right one.
    let state = pick_state();
    let vcf = vcf_line(1100, "names");

    for (name, reduces) in [
        ("pick", true),
        ("pick_allele", true),
        ("pick_allele_gene", true),
        ("flag_pick", false),
        ("flag_pick_allele", false),
        ("flag_pick_allele_gene", false),
    ] {
        let (status, body) =
            post_json(&state, "/api/annotate", json!({ "vcf": &vcf, name: true })).await;
        assert_eq!(status, StatusCode::OK, "{name}");
        let rows = body["results"][0]["transcript_consequences"]
            .as_array()
            .unwrap_or_else(|| panic!("{name}: no transcript consequences"));
        if reduces {
            assert_eq!(rows.len(), 1, "{name} should reduce to one row: {rows:?}");
        } else {
            assert_eq!(rows.len(), 2, "{name} should retain both rows: {rows:?}");
            assert_eq!(
                rows.iter()
                    .filter(|t| t["pick"].as_u64() == Some(1))
                    .count(),
                1,
                "{name} should mark exactly one row: {rows:?}"
            );
        }
    }
}

#[tokio::test]
async fn pick_does_not_drop_alleles_at_a_site_with_no_transcripts() {
    // An intergenic site is scaffolded one row per *alt allele*, not one row
    // per transcript, so running the pick hierarchy over those rows keeps one
    // and silently loses the other alt. There is no transcript to pick here,
    // so `pick` must leave the site alone.
    let state = pick_state();
    let vcf = "17\t900000\t.\tG\tA,T\t50\tPASS\t.";

    for pick in [false, true] {
        let (status, body) =
            post_json(&state, "/api/annotate", json!({ "vcf": vcf, "pick": pick })).await;
        assert_eq!(status, StatusCode::OK);
        let alleles: Vec<&str> = body["results"][0]["transcript_consequences"]
            .as_array()
            .expect("transcript consequences")
            .iter()
            .map(|t| t["variant_allele"].as_str().unwrap_or_default())
            .collect();
        assert_eq!(alleles, ["A", "T"], "pick={pick} lost an alt allele");
    }
}

#[tokio::test]
async fn pick_keeps_every_allele_of_the_transcript_it_picks() {
    // `pick` reduces to one transcript, not to one row: the JSON carries one
    // entry per (transcript, allele), so a biallelic site keeps two.
    let state = pick_state();
    let vcf = "17\t1100\t.\tG\tA,T\t50\tPASS\t.";

    let (status, body) =
        post_json(&state, "/api/annotate", json!({ "vcf": vcf, "pick": true })).await;
    assert_eq!(status, StatusCode::OK);
    let kept: Vec<(&str, &str)> = body["results"][0]["transcript_consequences"]
        .as_array()
        .expect("transcript consequences")
        .iter()
        .map(|t| {
            (
                t["transcript_id"].as_str().unwrap_or_default(),
                t["variant_allele"].as_str().unwrap_or_default(),
            )
        })
        .collect();
    assert_eq!(kept, [("TXA", "A"), ("TXA", "T")]);
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
    let state = replaceable_state();

    let (status, body) = get(&state, "/api/genomes").await;
    assert_eq!(status, StatusCode::OK);
    assert_eq!(body["genomes"], json!([]));

    let (status, body) = post_json(&state, "/api/load-genome", json!({ "name": "human" })).await;
    assert_eq!(status, StatusCode::BAD_REQUEST);
    assert_eq!(body["error"], "No data directory configured");
}

#[tokio::test]
async fn upload_gff3_replaces_the_active_gene_model() {
    let state = replaceable_state();

    let resp = post_gff3(&state, MINI_GFF3).await;
    let (status, body) = read_json(resp).await;

    assert_eq!(status, StatusCode::OK);
    assert_eq!(body["genes"], 1);
    assert_eq!(body["transcripts"], 1);
    // The replacement is counted, so a client can tell one model from the next
    // even though every upload reports the same `gff3_source`.
    assert_eq!(body["generation"], 1);

    // The swap is visible to the next caller, not just to this request.
    let (_, status_body) = get(&state, "/api/status").await;
    assert_eq!(status_body["transcripts"], 1);
    assert_eq!(status_body["gene_model_generation"], 1);
}

// ---- #113: the gene model a request was answered against ----

async fn post_gff3(state: &AppState, gff3: &str) -> axum::response::Response {
    router(state)
        .oneshot(
            Request::builder()
                .method("POST")
                .uri("/api/upload-gff3")
                .body(Body::from(gff3.to_string()))
                .unwrap(),
        )
        .await
        .unwrap()
}

#[tokio::test]
async fn replacing_the_gene_model_is_refused_unless_the_server_allows_it() {
    // The documented deployment is a lab server several colleagues share.
    // There, one client posting a GFF3 changed what every other client was
    // answered against - both authorised, neither doing anything wrong, so a
    // firewall does not help. Off by default makes the gene model a property
    // of how the server was started.
    let state = test_state();

    let (status, body) = read_json(post_gff3(&state, MINI_GFF3).await).await;
    assert_eq!(status, StatusCode::FORBIDDEN);
    let message = body["error"].as_str().unwrap();
    assert!(
        message.contains("--allow-model-replacement"),
        "the refusal has to name the flag that lifts it: {message}"
    );

    let (status, _) = post_json(&state, "/api/load-genome", json!({ "name": "human" })).await;
    assert_eq!(
        status,
        StatusCode::FORBIDDEN,
        "/api/load-genome replaces the model the same way /api/upload-gff3 does"
    );

    // Refused, and nothing moved.
    let (_, after) = get(&state, "/api/status").await;
    assert_eq!(after["transcripts"], 1);
    assert_eq!(after["gene_model_generation"], 0);
}

#[tokio::test]
async fn status_advertises_whether_the_model_can_be_replaced() {
    // So a UI can offer the model-replacing controls only when they will
    // work, rather than finding out with a 403 the user cannot interpret.
    let (_, locked) = get(&test_state(), "/api/status").await;
    assert_eq!(locked["allow_model_replacement"], false);

    let (_, open) = get(&replaceable_state(), "/api/status").await;
    assert_eq!(open["allow_model_replacement"], true);
}

#[tokio::test]
async fn annotate_reports_the_gene_model_that_answered() {
    // Before this, an `/api/annotate` response carried no indication of what
    // it was computed against: a client whose gene model had been replaced
    // mid-batch got a well-formed HTTP 200 answer against the wrong model,
    // with nothing to assert on. `/api/status` showed the swap, but only to a
    // client that polled it between every request.
    let state = replaceable_state();

    let (status, before) = post_json(
        &state,
        "/api/annotate",
        json!({ "vcf": vcf_line(43_100_000, "v1") }),
    )
    .await;
    assert_eq!(status, StatusCode::OK);
    assert_eq!(before["gene_model"]["transcripts"], 1);
    assert_eq!(before["gene_model"]["generation"], 0);

    // Another client replaces the model. The request below is byte-identical
    // to the one above.
    let (status, _) = read_json(post_gff3(&state, PICK_GFF3).await).await;
    assert_eq!(status, StatusCode::OK);

    let (status, after) = post_json(
        &state,
        "/api/annotate",
        json!({ "vcf": vcf_line(43_100_000, "v1") }),
    )
    .await;
    assert_eq!(status, StatusCode::OK);
    assert_eq!(
        after["gene_model"]["generation"], 1,
        "the same request answered against a different model has to say so"
    );
    assert_eq!(after["gene_model"]["gff3_source"], "user-upload");
    assert_ne!(
        before["gene_model"]["generation"], after["gene_model"]["generation"],
        "an undetectable wrong answer is what this field exists to prevent"
    );
}
