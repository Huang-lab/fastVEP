//! Contract tests for the JSON API documented in `docs/API.md`.
//!
//! These drive the real router through `oneshot`, so they cover routing, the
//! middleware stack, and serialization the way a client meets them. The point
//! is that the request and response shapes people write clients against are
//! asserted somewhere, not just described in a document that can drift.

use axum::body::Body;
use axum::http::{Request, StatusCode};
use fastvep_web::context::AnnotationContext;
use fastvep_web::genomes::default_genome_name;
use fastvep_web::handlers::{AppState, ServerConfig, SharedState};
use fastvep_web::limits::Limits;
use serde_json::{json, Value};
use std::path::PathBuf;
use std::sync::Arc;
use tower::ServiceExt;

const MINI_GFF3: &str = include_str!("../fixtures/mini.gff3");
const PICK_GFF3: &str = include_str!("../fixtures/pick.gff3");

/// What these tests start a server with. No stats file: they must not write to
/// the working directory.
fn config(data_dir: Option<PathBuf>, allow_model_replacement: bool) -> ServerConfig {
    ServerConfig {
        data_dir,
        sa_dir: None,
        stats_file: None,
        distance: 5000,
        allow_model_replacement,
        acmg_allow_partial: false,
        limits: Limits::default(),
    }
}

/// One transcript, no FASTA, no supplementary annotations. Deliberately the
/// minimum a caller can stand up, so a failure here is the API's fault rather
/// than the fixture's.
fn test_state() -> AppState {
    let mut ctx = AnnotationContext::new(None, None, None, 5000).expect("build context");
    ctx.update_gff3_text(MINI_GFF3).expect("load mini gff3");
    // The default a server ships with: the gene model is whatever it was
    // started with, and no request can move it (#113).
    Arc::new(SharedState::new(ctx, "mini", config(None, false)))
}

/// The same, started with `--allow-model-replacement` - the single-user
/// desktop configuration the browser UI's genome switching needs.
fn replaceable_state() -> AppState {
    let mut ctx = AnnotationContext::new(None, None, None, 5000).expect("build context");
    ctx.update_gff3_text(MINI_GFF3).expect("load mini gff3");
    Arc::new(SharedState::new(ctx, "mini", config(None, true)))
}

/// Two overlapping transcripts of one gene, the first-seen one non-canonical.
/// Kept separate from [`test_state`] so the single-transcript counts the other
/// tests assert stay meaningful.
fn pick_state() -> AppState {
    let mut ctx = AnnotationContext::new(None, None, None, 5000).expect("build context");
    ctx.update_gff3_text(PICK_GFF3).expect("load pick gff3");
    Arc::new(SharedState::new(ctx, "pick", config(None, false)))
}

/// A data directory of two genomes that answer differently at 17:1100: `alpha`
/// holds one transcript there and `beta` holds two. The server starts with
/// `alpha`, the way the hosted one starts with hg38, so a request that gets
/// the wrong genome shows in how many transcripts come back.
fn genome_dir(label: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!(
        "fastvep-web-genomes-{}-{}",
        label,
        std::process::id()
    ));
    let _ = std::fs::remove_dir_all(&dir);
    for (name, file, gff3) in [
        ("alpha", "mini.gff3", MINI_GFF3),
        ("beta", "pick.gff3", PICK_GFF3),
    ] {
        std::fs::create_dir_all(dir.join(name)).unwrap();
        std::fs::write(dir.join(name).join(file), gff3).unwrap();
    }
    dir
}

fn multi_genome_state(label: &str, allow_model_replacement: bool) -> AppState {
    let dir = genome_dir(label);
    let gff3 = dir.join("alpha").join("mini.gff3");
    let ctx = AnnotationContext::new(gff3.to_str(), None, None, 5000).expect("build context");
    let name = default_genome_name(gff3.to_str(), Some(dir.as_path()));
    Arc::new(SharedState::new(
        ctx,
        name,
        config(Some(dir), allow_model_replacement),
    ))
}

/// [`multi_genome_state`] started with `--acmg-allow-partial`, for the tests that look at what
/// ACMG returns and not at whether the genome may be classified.
fn multi_genome_state_allowing_partial_acmg(label: &str) -> AppState {
    let dir = genome_dir(label);
    let gff3 = dir.join("alpha").join("mini.gff3");
    let ctx = AnnotationContext::new(gff3.to_str(), None, None, 5000).expect("build context");
    let name = default_genome_name(gff3.to_str(), Some(dir.as_path()));
    Arc::new(SharedState::new(
        ctx,
        name,
        ServerConfig {
            acmg_allow_partial: true,
            ..config(Some(dir), false)
        },
    ))
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

/// The same one-transcript server, started with `--acmg-allow-partial`: the fixture has no
/// supplementary annotations, so without the switch it would refuse every request for ACMG.
fn partial_acmg_state() -> AppState {
    let mut ctx = AnnotationContext::new(None, None, None, 5000).expect("build context");
    ctx.update_gff3_text(MINI_GFF3).expect("load mini gff3");
    Arc::new(SharedState::new(
        ctx,
        "mini",
        ServerConfig {
            acmg_allow_partial: true,
            ..config(None, false)
        },
    ))
}

#[tokio::test]
async fn acmg_is_attached_only_when_requested() {
    let state = partial_acmg_state();
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
async fn without_a_data_dir_the_only_genome_is_the_one_the_server_started_with() {
    let state = replaceable_state();

    // The list says so, rather than coming back empty and leaving a client to
    // wonder whether the server has a genome at all.
    let (status, body) = get(&state, "/api/genomes").await;
    assert_eq!(status, StatusCode::OK);
    assert_eq!(body["default"], "mini");
    let genomes = body["genomes"].as_array().expect("genomes");
    assert_eq!(genomes.len(), 1, "got {genomes:?}");
    assert_eq!(genomes[0]["name"], "mini");
    assert_eq!(genomes[0]["default"], true);
    assert_eq!(genomes[0]["loaded"], true);
    assert_eq!(genomes[0]["transcripts"], 1);

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

// --- choosing the genome per request ---------------------------------------
//
// The hosted site serves several genomes from one process, with
// `--allow-model-replacement` off. A request names its genome and is answered
// from that genome's own model, so nothing is replaced and no client can move
// another's answers. What these assert is the part that must not go wrong
// quietly: the answer comes from the genome that was asked for, the response
// says which, and a name that is not a genome is never answered from the
// default instead.

/// The number of transcripts a one-variant request at 17:1100 returns: `alpha`
/// has one there and `beta` has two.
async fn transcripts_at_1100(state: &AppState, genome: Option<&str>) -> (usize, Value) {
    let mut body = json!({ "vcf": vcf_line(1100, "v") });
    if let Some(genome) = genome {
        body["genome"] = json!(genome);
    }
    let (status, response) = post_json(state, "/api/annotate", body).await;
    assert_eq!(status, StatusCode::OK, "{response}");
    let n = response["results"][0]["transcript_consequences"]
        .as_array()
        .expect("transcript consequences")
        .len();
    (n, response["gene_model"].clone())
}

#[tokio::test]
async fn status_names_the_default_genome() {
    let state = multi_genome_state("status", false);
    let (_, body) = get(&state, "/api/status").await;
    // This is how a client learns that the server answers per request: a
    // server that predates it does not send the field.
    assert_eq!(body["default_genome"], "alpha");
}

#[tokio::test]
async fn a_request_without_a_genome_is_answered_from_the_default_and_says_so() {
    let state = multi_genome_state("default", false);
    let (n, gene_model) = transcripts_at_1100(&state, None).await;
    assert_eq!(n, 1);
    assert_eq!(gene_model["genome"], "alpha");
}

#[tokio::test]
async fn a_request_is_answered_from_the_genome_it_names() {
    let state = multi_genome_state("named", false);

    let (n, gene_model) = transcripts_at_1100(&state, Some("beta")).await;
    assert_eq!(n, 2, "beta has two transcripts at this position");
    assert_eq!(gene_model["genome"], "beta");
    assert_eq!(gene_model["gff3_source"], "pick.gff3");
    assert_eq!(gene_model["transcripts"], 2);

    // Naming the genome the server started with is the same as naming none.
    let (n, gene_model) = transcripts_at_1100(&state, Some("alpha")).await;
    assert_eq!(n, 1);
    assert_eq!(gene_model["genome"], "alpha");

    // Choosing moved nothing: this needed no `--allow-model-replacement`, and
    // the default answers exactly as before.
    let (_, status) = get(&state, "/api/status").await;
    assert_eq!(status["default_genome"], "alpha");
    assert_eq!(status["transcripts"], 1);
    assert_eq!(status["gene_model_generation"], 0);
    let (n, _) = transcripts_at_1100(&state, None).await;
    assert_eq!(n, 1);
}

#[tokio::test]
async fn the_other_annotate_options_apply_to_a_named_genome() {
    let state = multi_genome_state_allowing_partial_acmg("options");
    let (status, body) = post_json(
        &state,
        "/api/annotate",
        json!({ "vcf": vcf_line(1100, "v"), "genome": "beta", "pick": true }),
    )
    .await;
    assert_eq!(status, StatusCode::OK);
    let kept = body["results"][0]["transcript_consequences"]
        .as_array()
        .expect("transcript consequences");
    assert_eq!(kept.len(), 1, "pick reduces a named genome's answer too");
    assert_eq!(kept[0]["transcript_id"], "TXA");
    assert_eq!(body["gene_model"]["genome"], "beta");

    let (status, body) = post_json(
        &state,
        "/api/annotate",
        json!({ "vcf": vcf_line(1100, "v"), "genome": "beta", "acmg": true }),
    )
    .await;
    assert_eq!(status, StatusCode::OK);
    assert!(
        body["results"][0]["transcript_consequences"][0]
            .get("acmg")
            .is_some(),
        "acmg was requested: {body}"
    );
}

#[tokio::test]
async fn an_unknown_genome_is_refused_and_never_answered_from_the_default() {
    let state = multi_genome_state("unknown", false);
    // A misspelling that quietly returned the default's consequences would
    // look exactly like a correct answer, which is the failure that matters.
    for name in [
        "gamma",
        "Alpha",
        "alph",
        "../alpha",
        "alpha/../beta",
        "a\\b",
    ] {
        let (status, body) = post_json(
            &state,
            "/api/annotate",
            json!({ "vcf": vcf_line(1100, "v"), "genome": name }),
        )
        .await;
        assert_eq!(status, StatusCode::BAD_REQUEST, "{name}: {body}");
        assert!(body.get("results").is_none(), "{name}: {body}");
        assert!(
            body["error"].as_str().unwrap().contains("not found"),
            "{name}: {body}"
        );
    }
    // None of those left anything behind.
    assert_eq!(state.genomes.loads_started(), 0);
}

#[tokio::test]
async fn an_empty_genome_name_is_refused_rather_than_read_as_no_genome() {
    let state = multi_genome_state("empty", false);
    let (status, body) = post_json(
        &state,
        "/api/annotate",
        json!({ "vcf": vcf_line(1100, "v"), "genome": "" }),
    )
    .await;
    assert_eq!(status, StatusCode::BAD_REQUEST, "{body}");
}

#[tokio::test]
async fn a_server_without_a_data_dir_serves_only_the_genome_it_started_with() {
    let state = test_state();
    let (status, body) = post_json(
        &state,
        "/api/annotate",
        json!({ "vcf": vcf_line(1100, "v"), "genome": "mouse" }),
    )
    .await;
    assert_eq!(status, StatusCode::BAD_REQUEST, "{body}");
    let message = body["error"].as_str().unwrap();
    assert!(
        message.contains("--data-dir") && message.contains("'mini'"),
        "the refusal says what the server does serve: {message}"
    );

    let (n, gene_model) = transcripts_at_1100(&state, Some("mini")).await;
    assert_eq!(n, 1);
    assert_eq!(gene_model["genome"], "mini");
}

#[tokio::test]
async fn the_genome_list_says_what_is_loaded_and_which_is_the_default() {
    let state = multi_genome_state("list", false);

    let (status, before) = get(&state, "/api/genomes").await;
    assert_eq!(status, StatusCode::OK);
    assert_eq!(before["default"], "alpha");
    let genomes = before["genomes"].as_array().expect("genomes");
    let by_name = |name: &str| genomes.iter().find(|g| g["name"] == name).expect(name);
    assert_eq!(by_name("alpha")["default"], true);
    assert_eq!(by_name("alpha")["loaded"], true);
    assert_eq!(by_name("alpha")["transcripts"], 1);
    assert_eq!(by_name("alpha")["gff3"], "mini.gff3");
    assert_eq!(by_name("beta")["default"], false);
    assert_eq!(by_name("beta")["loaded"], false);
    assert_eq!(by_name("beta")["transcripts"], Value::Null);
    assert_eq!(by_name("beta")["gff3"], "pick.gff3");

    transcripts_at_1100(&state, Some("beta")).await;

    let (_, after) = get(&state, "/api/genomes").await;
    let genomes = after["genomes"].as_array().expect("genomes");
    let beta = genomes.iter().find(|g| g["name"] == "beta").unwrap();
    assert_eq!(beta["loaded"], true);
    assert_eq!(beta["transcripts"], 2);
}

#[tokio::test(flavor = "multi_thread", worker_threads = 4)]
async fn concurrent_first_requests_for_a_genome_share_one_load() {
    let state = multi_genome_state("concurrent", false);
    let requests = (0..8).map(|_| {
        let state = Arc::clone(&state);
        async move { transcripts_at_1100(&state, Some("beta")).await }
    });
    let answers = futures_join_all(requests).await;
    for (n, gene_model) in &answers {
        assert_eq!(*n, 2);
        assert_eq!(gene_model["genome"], "beta");
    }
    assert_eq!(
        state.genomes.loads_started(),
        1,
        "eight simultaneous first requests are one load, not eight"
    );
}

/// `futures::future::join_all` without the dependency: the requests are
/// independent, so polling them from separate tasks is the same thing.
async fn futures_join_all<F>(futures: impl IntoIterator<Item = F>) -> Vec<F::Output>
where
    F: std::future::Future + Send + 'static,
    F::Output: Send + 'static,
{
    let handles: Vec<_> = futures.into_iter().map(tokio::spawn).collect();
    let mut out = Vec::with_capacity(handles.len());
    for handle in handles {
        out.push(handle.await.expect("task"));
    }
    out
}

#[tokio::test]
async fn preloading_loads_every_other_genome_and_skips_the_default() {
    let state = multi_genome_state("preload", false);
    assert_eq!(state.genomes.loads_started(), 0);

    let loaded = Arc::clone(&state).preload_genomes().await;
    assert_eq!(loaded, 1, "beta; alpha is the one the server started with");
    assert_eq!(state.genomes.loads_started(), 1);

    let (_, listing) = get(&state, "/api/genomes").await;
    let genomes = listing["genomes"].as_array().expect("genomes");
    assert!(
        genomes.iter().all(|g| g["loaded"] == true),
        "everything is loaded: {genomes:?}"
    );
}

#[tokio::test]
async fn named_runs_are_the_ones_counted_as_genome_sessions() {
    let state = multi_genome_state("counter", false);
    transcripts_at_1100(&state, None).await;
    let (_, after_default) = get(&state, "/api/status").await;
    assert_eq!(after_default["total_variants"], 1);
    assert_eq!(after_default["total_genomes"], 0);

    transcripts_at_1100(&state, Some("beta")).await;
    transcripts_at_1100(&state, Some("alpha")).await;
    let (_, after_named) = get(&state, "/api/status").await;
    assert_eq!(after_named["total_variants"], 3);
    assert_eq!(after_named["total_genomes"], 2);
}

#[tokio::test]
async fn replacing_the_default_model_renames_it_so_the_old_name_is_not_answered_from_the_new_model()
{
    // The trap in routing by name: `alpha` is the default, a client replaces
    // the default with `beta`, and a later request for `alpha` is answered
    // from the model that is now in the default slot, which is beta's. Both
    // answers are well-formed, and only one of them is about alpha.
    let state = multi_genome_state("rename", true);

    let (status, _) = post_json(&state, "/api/load-genome", json!({ "name": "beta" })).await;
    assert_eq!(status, StatusCode::OK);
    let (_, now) = get(&state, "/api/status").await;
    assert_eq!(now["default_genome"], "beta");

    let (n, gene_model) = transcripts_at_1100(&state, Some("alpha")).await;
    assert_eq!(
        n, 1,
        "alpha has one transcript here; beta, now the default, has two"
    );
    assert_eq!(gene_model["genome"], "alpha");
    assert_eq!(
        gene_model["generation"], 0,
        "alpha was loaded, not replaced"
    );

    let (n, gene_model) = transcripts_at_1100(&state, Some("beta")).await;
    assert_eq!(n, 2);
    assert_eq!(gene_model["genome"], "beta");
    assert_eq!(
        gene_model["generation"], 1,
        "beta answered from the default slot, which has been replaced once"
    );

    // An upload is not any genome of the data directory either.
    let (status, _) = read_json(post_gff3(&state, PICK_GFF3).await).await;
    assert_eq!(status, StatusCode::OK);
    let (_, now) = get(&state, "/api/status").await;
    assert_eq!(now["default_genome"], "user-upload");
    let (n, gene_model) = transcripts_at_1100(&state, Some("alpha")).await;
    assert_eq!(n, 1, "the upload has two transcripts here; alpha has one");
    assert_eq!(gene_model["genome"], "alpha");
}

// ---------------------------------------------------------------------------
// ACMG needs the reference stack
// ---------------------------------------------------------------------------

/// Write a database directory that holds all nine sources of the reference stack, each with
/// one record, so that the files open and carry the keys the classifier reads them by.
fn write_reference_stack(sa_dir: &std::path::Path) {
    use fastvep_sa::common::{GeneRecord, IntervalRecord};
    use fastvep_sa::fields::{Field, FieldType};
    use fastvep_sa::gene::{GeneHeader, GeneIndex};
    use fastvep_sa::interval::{IntervalHeader, IntervalIndex};
    use fastvep_sa::writer_v2::{Osa2Metadata, Osa2Record, Osa2Writer};

    std::fs::create_dir_all(sa_dir).unwrap();
    let field = Field {
        field: "AF".into(),
        alias: "af".into(),
        ftype: FieldType::Float,
        multiplier: 1_000_000,
        zigzag: false,
        missing_value: u32::MAX,
        missing_string: ".".into(),
        description: "allele frequency".into(),
    };
    for (stem, key) in [
        ("gnomad", "gnomad"),
        ("clinvar", "clinvar"),
        ("revel", "revel"),
        ("spliceai", "spliceAI"),
        ("phylop", "phylop"),
    ] {
        let metadata = Osa2Metadata {
            format_version: 2,
            name: key.into(),
            version: "test".into(),
            assembly: "GRCh38".into(),
            json_key: key.into(),
            match_by_allele: true,
            is_array: false,
            is_positional: false,
            chunk_bits: 20,
            description: String::new(),
        };
        // One record on chromosome 1 and one on 17, the chromosome the gene models of these tests
        // have their transcripts on: a gnomAD that held chromosome 1 alone would not be a complete
        // stack for them.
        let records: Vec<Osa2Record> = ["chr1", "chr17"]
            .iter()
            .map(|chrom| Osa2Record {
                chrom: (*chrom).into(),
                position: 1000,
                ref_allele: b"A".to_vec(),
                alt_allele: b"G".to_vec(),
                values: vec![field.encode_float(0.25)],
                json_blob: None,
            })
            .collect();
        let file = std::fs::File::create(sa_dir.join(format!("{stem}.osa2"))).unwrap();
        Osa2Writer::new(metadata, vec![field.clone()])
            .write_all(std::io::BufWriter::new(file), &records)
            .unwrap();
    }
    for key in ["gnomad_genes", "omim", "clinvar_protein"] {
        let mut index = GeneIndex::new(GeneHeader {
            schema_version: fastvep_sa::common::SCHEMA_VERSION,
            json_key: key.into(),
            name: key.into(),
            version: "test".into(),
            assembly: "GRCh38".into(),
        });
        index.add(GeneRecord {
            gene_symbol: "NOSUCHGENE".into(),
            json: "{}".into(),
        });
        let mut file = std::fs::File::create(sa_dir.join(format!("{key}.oga"))).unwrap();
        index.write_to(&mut file).unwrap();
    }
    let mut repeats = IntervalIndex::new(IntervalHeader {
        schema_version: fastvep_sa::common::SCHEMA_VERSION,
        json_key: "repeatmasker".into(),
        name: "RepeatMasker".into(),
        version: "test".into(),
        assembly: "GRCh38".into(),
    });
    repeats.add(IntervalRecord {
        chrom: "chr1".into(),
        start: 10,
        end: 20,
        json: "{}".into(),
    });
    repeats.sort();
    let mut file = std::fs::File::create(sa_dir.join("repeatmasker.osi")).unwrap();
    repeats.write_to(&mut file).unwrap();
}

/// Three genomes: `alpha` (the default, with the whole stack in its `sa/`), `beta` with no
/// supplementary annotations, and `gamma` with the whole stack, which is not loaded until asked
/// for.
fn genomes_with_and_without_the_stack(label: &str, allow_partial: bool) -> AppState {
    let dir = genome_dir(label);
    write_reference_stack(&dir.join("alpha").join("sa"));
    std::fs::create_dir_all(dir.join("gamma")).unwrap();
    std::fs::write(dir.join("gamma").join("mini.gff3"), MINI_GFF3).unwrap();
    write_reference_stack(&dir.join("gamma").join("sa"));
    let gff3 = dir.join("alpha").join("mini.gff3");
    let sa = dir.join("alpha").join("sa");
    let ctx =
        AnnotationContext::new(gff3.to_str(), None, sa.to_str(), 5000).expect("build context");
    let name = default_genome_name(gff3.to_str(), Some(dir.as_path()));
    Arc::new(SharedState::new(
        ctx,
        name,
        ServerConfig {
            acmg_allow_partial: allow_partial,
            ..config(Some(dir), false)
        },
    ))
}

#[tokio::test]
async fn status_says_whether_acmg_will_be_answered_and_what_stands_in_the_way() {
    let (status, body) = get(&test_state(), "/api/status").await;
    assert_eq!(status, StatusCode::OK);
    let acmg = &body["acmg"];
    assert_eq!(acmg["ready"], false, "{acmg}");
    assert_eq!(acmg["available"], false, "{acmg}");
    assert_eq!(acmg["loaded"], json!([]));
    let missing = acmg["missing"]
        .as_array()
        .expect("a list of missing sources");
    assert_eq!(missing.len(), 9);
    let gnomad = missing
        .iter()
        .find(|m| m["id"] == "gnomad")
        .expect("gnomAD is listed");
    assert_eq!(gnomad["drives"], json!(["BA1", "BS1", "BS2", "PM2"]));
    assert!(gnomad["without"]
        .as_str()
        .unwrap()
        .contains("population frequency"));

    let (_, body) = get(
        &genomes_with_and_without_the_stack("status-full", false),
        "/api/status",
    )
    .await;
    assert_eq!(body["acmg"]["ready"], true, "{}", body["acmg"]);
    assert_eq!(body["acmg"]["available"], true);
    assert_eq!(body["acmg"]["missing"], json!([]));
    assert_eq!(body["acmg"]["loaded"].as_array().unwrap().len(), 9);
}

#[tokio::test]
async fn a_request_for_acmg_on_a_genome_without_the_stack_is_refused_and_says_why() {
    let state = test_state();
    let vcf = vcf_line(1100, ".");

    let (status, body) =
        post_json(&state, "/api/annotate", json!({ "vcf": vcf, "acmg": true })).await;
    assert_eq!(status, StatusCode::UNPROCESSABLE_ENTITY, "{body}");
    assert!(
        body.get("results").is_none(),
        "nothing is annotated: {body}"
    );
    let message = body["error"].as_str().unwrap();
    assert!(
        message.contains("not available for genome 'mini'"),
        "{message}"
    );
    assert!(message.contains("9 of the 9 sources"), "{message}");
    assert_eq!(body["acmg"]["available"], false);
    assert_eq!(body["acmg"]["missing"].as_array().unwrap().len(), 9);

    // The same variant without ACMG is answered: only the classification is refused.
    let (status, body) = post_json(&state, "/api/annotate", json!({ "vcf": vcf })).await;
    assert_eq!(status, StatusCode::OK, "{body}");
    assert!(body.get("acmg_stack").is_none(), "{body}");
    // A refused request is not a run, and is not counted as one.
    assert_eq!(
        state
            .total_variants
            .load(std::sync::atomic::Ordering::Relaxed),
        1
    );
}

#[tokio::test]
async fn a_partial_stack_is_answered_when_the_operator_allowed_it_and_the_response_says_so() {
    let (status, body) = post_json(
        &partial_acmg_state(),
        "/api/annotate",
        json!({ "vcf": vcf_line(1100, "."), "acmg": true }),
    )
    .await;
    assert_eq!(status, StatusCode::OK, "{body}");
    assert_eq!(body["acmg_stack"]["complete"], false);
    assert_eq!(body["acmg_stack"]["missing"].as_array().unwrap().len(), 9);
    assert!(body["acmg_stack"]["missing"]
        .as_array()
        .unwrap()
        .contains(&json!("gnomad")));
    // And PM2 does not read the missing database as absence from controls.
    let criteria = body["results"][0]["transcript_consequences"][0]["acmg"]["criteria"]
        .as_array()
        .unwrap();
    let pm2 = criteria
        .iter()
        .find(|c| c["code"].as_str().unwrap().starts_with("PM2"))
        .expect("PM2 is reported");
    assert_eq!(pm2["evaluated"], false, "{pm2}");
    assert_eq!(pm2["met"], false, "{pm2}");
}

#[tokio::test]
async fn a_genome_with_the_whole_stack_is_classified_and_the_response_says_the_stack_was_complete()
{
    let state = genomes_with_and_without_the_stack("full", false);
    // The default genome (alpha) and one loaded on request (gamma) both hold it.
    for genome in [None, Some("gamma")] {
        let mut request = json!({ "vcf": vcf_line(1100, "."), "acmg": true });
        if let Some(g) = genome {
            request["genome"] = json!(g);
        }
        let (status, body) = post_json(&state, "/api/annotate", request).await;
        assert_eq!(status, StatusCode::OK, "{genome:?}: {body}");
        assert_eq!(
            body["acmg_stack"],
            json!({ "complete": true, "missing": [] }),
            "{genome:?}"
        );
        let acmg = &body["results"][0]["transcript_consequences"][0]["acmg"];
        assert!(acmg["classification"].is_string(), "{genome:?}: {body}");
    }
    // A genome of the same server that lacks the stack is still refused.
    let (status, body) = post_json(
        &state,
        "/api/annotate",
        json!({ "vcf": vcf_line(1100, "."), "acmg": true, "genome": "beta" }),
    )
    .await;
    assert_eq!(status, StatusCode::UNPROCESSABLE_ENTITY, "{body}");
    assert!(body["error"].as_str().unwrap().contains("genome 'beta'"));
}

/// `alpha` with the whole stack except that its gnomAD holds chromosome 1 alone, while the gene
/// model has its transcripts on chromosome 17: what a gnomAD of one chromosome, built by a script
/// that gave every chromosome the same name, looks like to the classifier.
fn state_with_a_gnomad_of_another_chromosome(label: &str, allow_partial: bool) -> AppState {
    let dir = genome_dir(label);
    let sa = dir.join("alpha").join("sa");
    write_reference_stack(&sa);
    write_osa2(&sa, "gnomad", "gnomad", "test", 2);
    let gff3 = dir.join("alpha").join("mini.gff3");
    let ctx =
        AnnotationContext::new(gff3.to_str(), None, sa.to_str(), 5000).expect("build context");
    let name = default_genome_name(gff3.to_str(), Some(dir.as_path()));
    Arc::new(SharedState::new(
        ctx,
        name,
        ServerConfig {
            acmg_allow_partial: allow_partial,
            ..config(Some(dir), false)
        },
    ))
}

#[tokio::test]
async fn a_gnomad_of_other_chromosomes_than_the_gene_model_is_not_a_complete_stack() {
    let state = state_with_a_gnomad_of_another_chromosome("gnomad-elsewhere", false);
    let (_, status) = get(&state, "/api/status").await;
    let acmg = &status["acmg"];
    assert_eq!(acmg["ready"], false, "{acmg}");
    assert_eq!(acmg["available"], false, "{acmg}");
    assert!(
        !acmg["loaded"]
            .as_array()
            .unwrap()
            .contains(&json!("gnomad")),
        "{acmg}"
    );
    let first = &acmg["missing"][0];
    assert_eq!(first["id"], "gnomad", "{acmg}");
    assert!(
        first["label"]
            .as_str()
            .unwrap()
            .contains("none for chromosomes 17"),
        "{first}"
    );
    assert!(
        first["without"]
            .as_str()
            .unwrap()
            .contains("PM2 would fire"),
        "{first}"
    );
    // The other eight are in place: only gnomAD is named.
    assert_eq!(acmg["missing"].as_array().unwrap().len(), 1, "{acmg}");

    let (code, body) = post_json(
        &state,
        "/api/annotate",
        json!({ "vcf": vcf_line(1100, "."), "acmg": true }),
    )
    .await;
    assert_eq!(code, StatusCode::UNPROCESSABLE_ENTITY, "{body}");
    assert!(
        body["error"]
            .as_str()
            .unwrap()
            .contains("none for chromosomes 17"),
        "{body}"
    );
}

#[tokio::test]
async fn with_partial_acmg_allowed_a_variant_on_a_chromosome_gnomad_lacks_is_not_called_absent() {
    let state = state_with_a_gnomad_of_another_chromosome("gnomad-elsewhere-allowed", true);
    let (code, body) = post_json(
        &state,
        "/api/annotate",
        json!({ "vcf": vcf_line(1100, "."), "acmg": true }),
    )
    .await;
    assert_eq!(code, StatusCode::OK, "{body}");
    assert_eq!(body["acmg_stack"]["complete"], false, "{body}");
    assert_eq!(body["acmg_stack"]["missing"], json!(["gnomad"]), "{body}");
    let criteria = body["results"][0]["transcript_consequences"][0]["acmg"]["criteria"]
        .as_array()
        .unwrap();
    // Reported as PM2 or, at the strength the default gives it, PM2_Supporting.
    let pm2 = criteria
        .iter()
        .find(|c| {
            c["code"]
                .as_str()
                .is_some_and(|code| code.starts_with("PM2"))
        })
        .expect("PM2 is reported");
    assert_eq!(pm2["met"], false, "{pm2}");
    assert_eq!(pm2["evaluated"], false, "{pm2}");
    assert!(
        pm2["summary"]
            .as_str()
            .unwrap()
            .contains("holds no record of chromosome 17"),
        "{pm2}"
    );
}

/// Corrupt the payload of one ZIP entry in place, leaving every offset and length untouched, so
/// that the reader still finds the entry (the database opens, and lists the chromosome) and fails
/// only when it has to inflate it.
fn corrupt_entry_payload(bytes: &mut [u8], name: &str) -> bool {
    let name_bytes = name.as_bytes();
    let mut i = 0;
    while i + 30 + name_bytes.len() <= bytes.len() {
        if &bytes[i..i + 4] == b"PK\x03\x04" {
            let name_len = u16::from_le_bytes([bytes[i + 26], bytes[i + 27]]) as usize;
            let extra_len = u16::from_le_bytes([bytes[i + 28], bytes[i + 29]]) as usize;
            let name_start = i + 30;
            if name_len == name_bytes.len()
                && &bytes[name_start..name_start + name_len] == name_bytes
            {
                let data = name_start + name_len + extra_len;
                for b in bytes.iter_mut().skip(data).take(16) {
                    *b ^= 0xFF;
                }
                return true;
            }
        }
        i += 1;
    }
    false
}

#[tokio::test]
async fn a_database_that_cannot_be_read_at_query_time_says_so_in_the_answer_and_in_the_status() {
    use fastvep_sa::fields::{Field, FieldType};
    use fastvep_sa::writer_v2::{Osa2Metadata, Osa2Record, Osa2Writer};

    let dir = genome_dir("unreadable-chunk");
    let sa = dir.join("alpha").join("sa");
    write_reference_stack(&sa);
    // A gnomAD that opens, holds chromosome 17, and cannot inflate the chunk the variant is in: the
    // answer for that variant is "no record", and nothing in it said that this was not "none".
    let field = Field {
        field: "AF".into(),
        alias: "af".into(),
        ftype: FieldType::Float,
        multiplier: 1_000_000,
        zigzag: false,
        missing_value: u32::MAX,
        missing_string: ".".into(),
        description: String::new(),
    };
    let metadata = Osa2Metadata {
        format_version: 2,
        name: "gnomad".into(),
        version: "test".into(),
        assembly: "GRCh38".into(),
        json_key: "gnomad".into(),
        match_by_allele: true,
        is_array: false,
        is_positional: false,
        chunk_bits: 20,
        description: String::new(),
    };
    let record = Osa2Record {
        chrom: "chr17".into(),
        position: 1100,
        ref_allele: b"G".to_vec(),
        alt_allele: b"A".to_vec(),
        values: vec![field.encode_float(0.25)],
        json_blob: None,
    };
    let path = sa.join("gnomad.osa2");
    let file = std::fs::File::create(&path).unwrap();
    Osa2Writer::new(metadata, vec![field])
        .write_all(std::io::BufWriter::new(file), &[record])
        .unwrap();
    let mut bytes = std::fs::read(&path).unwrap();
    assert!(corrupt_entry_payload(
        &mut bytes,
        "fastsa/chr17/0/var32.bin"
    ));
    std::fs::write(&path, bytes).unwrap();

    let gff3 = dir.join("alpha").join("mini.gff3");
    let ctx =
        AnnotationContext::new(gff3.to_str(), None, sa.to_str(), 5000).expect("build context");
    let name = default_genome_name(gff3.to_str(), Some(dir.as_path()));
    let state = Arc::new(SharedState::new(ctx, name, config(Some(dir), false)));

    let (code, body) = post_json(
        &state,
        "/api/annotate",
        json!({ "vcf": vcf_line(1100, ".") }),
    )
    .await;
    assert_eq!(code, StatusCode::OK, "{body}");
    let warnings = body["warnings"]
        .as_array()
        .expect("the answer carries warnings");
    assert!(
        warnings.iter().any(|w| w
            .as_str()
            .unwrap()
            .contains("gnomad could not be read for 1 variant lookup")),
        "{warnings:?}"
    );
    assert!(
        warnings[0]
            .as_str()
            .unwrap()
            .contains("not the same as having none"),
        "{warnings:?}"
    );

    // The tally is the server's, so the status says it too, for the monitor that never sends a VCF.
    let (_, status) = get(&state, "/api/status").await;
    assert!(
        status["sa_lookup_errors"]["gnomad"].as_u64().unwrap() >= 1,
        "{status}"
    );
}

#[tokio::test]
async fn the_genome_list_says_for_each_genome_whether_it_can_be_classified_before_it_is_loaded() {
    let state = genomes_with_and_without_the_stack("listing", false);
    let (_, body) = get(&state, "/api/genomes").await;
    let by_name = |name: &str| -> Value {
        body["genomes"]
            .as_array()
            .unwrap()
            .iter()
            .find(|g| g["name"] == name)
            .unwrap_or_else(|| panic!("{name} is listed: {body}"))
            .clone()
    };
    // The default one: loaded, with the stack.
    assert_eq!(by_name("alpha")["acmg"]["ready"], true);
    // gamma has not been asked for, and its directory says it can be classified.
    let gamma = by_name("gamma");
    assert_eq!(gamma["loaded"], false);
    assert_eq!(gamma["acmg"]["ready"], true, "{gamma}");
    assert_eq!(gamma["acmg"]["loaded"].as_array().unwrap().len(), 9);
    // beta has no sa/ directory at all.
    let beta = by_name("beta");
    assert_eq!(beta["acmg"]["ready"], false);
    assert_eq!(beta["acmg"]["missing"].as_array().unwrap().len(), 9);

    // Loading gamma changes nothing about what the list says of it.
    let (status, _) = post_json(
        &state,
        "/api/annotate",
        json!({ "vcf": vcf_line(1100, "."), "genome": "gamma" }),
    )
    .await;
    assert_eq!(status, StatusCode::OK);
    let (_, body) = get(&state, "/api/genomes").await;
    let gamma = body["genomes"]
        .as_array()
        .unwrap()
        .iter()
        .find(|g| g["name"] == "gamma")
        .unwrap();
    assert_eq!(gamma["loaded"], true);
    assert_eq!(gamma["acmg"]["ready"], true);
}

#[tokio::test]
async fn with_partial_acmg_allowed_the_list_says_available_but_not_ready() {
    let state = genomes_with_and_without_the_stack("listing-partial", true);
    let (_, body) = get(&state, "/api/genomes").await;
    let beta = body["genomes"]
        .as_array()
        .unwrap()
        .iter()
        .find(|g| g["name"] == "beta")
        .unwrap();
    assert_eq!(beta["acmg"]["ready"], false);
    assert_eq!(beta["acmg"]["available"], true);
}

// ---------------------------------------------------------------------------
// What one request may cost
// ---------------------------------------------------------------------------

/// The one-transcript server started with these limits.
fn state_with_limits(limits: Limits) -> AppState {
    let mut ctx = AnnotationContext::new(None, None, None, 5000).expect("build context");
    ctx.update_gff3_text(MINI_GFF3).expect("load mini gff3");
    Arc::new(SharedState::new(
        ctx,
        "mini",
        ServerConfig {
            limits,
            ..config(None, false)
        },
    ))
}

fn many_records(n: u32) -> String {
    (0..n)
        .map(|i| vcf_line(1000 + i, &format!("v{i}")))
        .collect::<Vec<_>>()
        .join("\n")
}

#[tokio::test]
async fn a_line_that_cannot_be_read_is_a_400_that_names_it_and_nothing_is_annotated() {
    let state = test_state();
    for (vcf, line) in [
        ("hello world", 1),
        ("17\t0\t.\tG\tA\t50\tPASS\t.", 1),
        (
            "17\t1100\t.\tG\tA\t50\tPASS\t.\n17\tfoo\t.\tG\tA\t50\tPASS\t.",
            2,
        ),
        ("##fileformat=VCFv4.2\n#CHROM\tPOS\n17\t1100\t.\tG", 3),
    ] {
        let (status, body) = post_json(&state, "/api/annotate", json!({ "vcf": vcf })).await;
        assert_eq!(status, StatusCode::BAD_REQUEST, "{vcf:?}: {body}");
        let message = body["error"].as_str().unwrap();
        assert!(
            message.contains(&format!("VCF line {line} cannot be read")),
            "{vcf:?}: {message}"
        );
        assert!(body.get("results").is_none(), "{body}");
    }
    // Blank lines and comments among the records are not records, as the parser reads them.
    let (status, body) = post_json(
        &state,
        "/api/annotate",
        json!({ "vcf": format!("\n# a comment\n{}\n\n", vcf_line(1100, ".")) }),
    )
    .await;
    assert_eq!(status, StatusCode::OK, "{body}");
    assert_eq!(body["count"], 1);
    // A refused request is not a run.
    assert_eq!(
        state
            .total_variants
            .load(std::sync::atomic::Ordering::Relaxed),
        1
    );
}

#[tokio::test]
async fn annotating_in_batches_gives_the_answer_of_annotating_all_at_once() {
    let state = state_with_limits(Limits {
        batch_variants: 3,
        ..Limits::default()
    });
    let vcf = format!(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n{}",
        many_records(10)
    );
    let (status, body) = post_json(&state, "/api/annotate", json!({ "vcf": vcf })).await;
    assert_eq!(status, StatusCode::OK, "{body}");
    assert_eq!(body["count"], 10);

    let whole = {
        let guard = state.ctx.read().unwrap();
        guard.annotate_vcf_text_with_acmg(&vcf, None, None).unwrap()
    };
    assert_eq!(
        body["results"],
        Value::Array(whole),
        "four batches of three, three, three and one"
    );
    // The fields that follow the results are all there, whichever order they come in.
    for field in ["count", "time_ms", "gene_model"] {
        assert!(body.get(field).is_some(), "{field}: {body}");
    }
}

#[tokio::test]
async fn an_answer_that_passes_the_response_limit_is_refused_and_the_server_goes_on() {
    let state = state_with_limits(Limits {
        max_response_bytes: 2_048,
        batch_variants: 2,
        ..Limits::default()
    });
    let (status, body) =
        post_json(&state, "/api/annotate", json!({ "vcf": many_records(10) })).await;
    assert_eq!(status, StatusCode::PAYLOAD_TOO_LARGE, "{body}");
    let message = body["error"].as_str().unwrap();
    assert!(message.contains("larger than the 2 KiB"), "{message}");
    assert!(message.contains("of 10 variants were done"), "{message}");
    assert!(message.contains("`pick`"), "{message}");
    assert!(body.get("results").is_none());
    assert_eq!(
        state
            .total_variants
            .load(std::sync::atomic::Ordering::Relaxed),
        0
    );

    // The server is as it was: a request that fits is answered.
    let (status, body) = post_json(
        &state,
        "/api/annotate",
        json!({ "vcf": vcf_line(1100, "."), "pick": true }),
    )
    .await;
    assert_eq!(status, StatusCode::OK, "{body}");
}

#[tokio::test]
async fn a_request_that_runs_past_the_time_limit_is_refused_with_what_was_done() {
    let state = state_with_limits(Limits {
        max_request_time: std::time::Duration::ZERO,
        batch_variants: 1,
        ..Limits::default()
    });
    let (status, body) =
        post_json(&state, "/api/annotate", json!({ "vcf": many_records(3) })).await;
    assert_eq!(status, StatusCode::PAYLOAD_TOO_LARGE, "{body}");
    assert!(
        body["error"]
            .as_str()
            .unwrap()
            .contains("longer than the 0 seconds"),
        "{body}"
    );
    assert!(
        body["error"]
            .as_str()
            .unwrap()
            .contains("of 3 variants were done"),
        "{body}"
    );
}

#[tokio::test]
async fn a_request_that_finds_every_place_taken_is_turned_away_and_told_when_to_come_back() {
    let state = state_with_limits(Limits {
        max_parallel: 1,
        queue_timeout: std::time::Duration::from_millis(50),
        ..Limits::default()
    });
    let held = Arc::clone(&state.compute).acquire_owned().await.unwrap();
    let resp = router(&state)
        .oneshot(
            Request::builder()
                .method("POST")
                .uri("/api/annotate")
                .header("content-type", "application/json")
                .body(Body::from(
                    json!({ "vcf": vcf_line(1100, ".") }).to_string(),
                ))
                .unwrap(),
        )
        .await
        .unwrap();
    assert_eq!(resp.status(), StatusCode::SERVICE_UNAVAILABLE);
    assert_eq!(resp.headers().get("retry-after").unwrap(), "5");
    let (_, body) = read_json(resp).await;
    assert!(body["error"].as_str().unwrap().contains("busy"), "{body}");

    // Reads do not need a place, and the place comes back.
    let (status, _) = get(&state, "/api/status").await;
    assert_eq!(status, StatusCode::OK);
    drop(held);
    let (status, body) = post_json(
        &state,
        "/api/annotate",
        json!({ "vcf": vcf_line(1100, ".") }),
    )
    .await;
    assert_eq!(status, StatusCode::OK, "{body}");
}

#[tokio::test]
async fn a_client_that_goes_away_stops_the_annotation_and_frees_its_place() {
    let state = state_with_limits(Limits {
        max_parallel: 1,
        batch_variants: 1,
        ..Limits::default()
    });
    let work = post_json(
        &state,
        "/api/annotate",
        json!({ "vcf": many_records(30_000) }),
    );
    // The client gives up almost at once; axum drops the handler's future.
    let _ = tokio::time::timeout(std::time::Duration::from_millis(30), work).await;
    // The annotation sees the flag between batches and stops: it never counts as a run, and its
    // place comes back long before 30,000 batches could have run.
    let started = std::time::Instant::now();
    while state.compute.available_permits() == 0
        && started.elapsed() < std::time::Duration::from_secs(5)
    {
        tokio::time::sleep(std::time::Duration::from_millis(20)).await;
    }
    assert_eq!(state.compute.available_permits(), 1, "the place came back");
    assert!(
        started.elapsed() < std::time::Duration::from_secs(4),
        "{:?}",
        started.elapsed()
    );
    assert_eq!(
        state
            .total_variants
            .load(std::sync::atomic::Ordering::Relaxed),
        0
    );
}

#[tokio::test]
async fn status_reports_what_one_request_may_cost() {
    let (_, body) = get(&test_state(), "/api/status").await;
    let limits = &body["limits"];
    assert_eq!(limits["max_body_bytes"], 10 * 1024 * 1024);
    assert_eq!(limits["max_response_bytes"], 128 * 1024 * 1024);
    assert_eq!(limits["max_request_seconds"], 90);
    assert!(limits["max_parallel"].as_u64().unwrap() >= 1);
}

// ---------------------------------------------------------------------------
// What each source is, and what was not loaded
// ---------------------------------------------------------------------------

fn write_osa2(dir: &std::path::Path, stem: &str, key: &str, version: &str, records: usize) {
    use fastvep_sa::fields::{Field, FieldType};
    use fastvep_sa::writer_v2::{Osa2Metadata, Osa2Record, Osa2Writer};
    let field = Field {
        field: "AF".into(),
        alias: "af".into(),
        ftype: FieldType::Float,
        multiplier: 1_000_000,
        zigzag: false,
        missing_value: u32::MAX,
        missing_string: ".".into(),
        description: "allele frequency".into(),
    };
    let metadata = Osa2Metadata {
        format_version: 2,
        name: key.into(),
        version: version.into(),
        assembly: "GRCh38".into(),
        json_key: key.into(),
        match_by_allele: true,
        is_array: false,
        is_positional: false,
        chunk_bits: 20,
        description: String::new(),
    };
    let rows: Vec<Osa2Record> = (0..records)
        .map(|i| Osa2Record {
            chrom: "chr1".into(),
            position: 1000 + i as u32,
            ref_allele: b"A".to_vec(),
            alt_allele: b"G".to_vec(),
            values: vec![field.encode_float(0.25)],
            json_blob: None,
        })
        .collect();
    let file = std::fs::File::create(dir.join(format!("{stem}.osa2"))).unwrap();
    Osa2Writer::new(metadata, vec![field])
        .write_all(std::io::BufWriter::new(file), &rows)
        .unwrap();
}

#[tokio::test]
async fn status_says_what_each_source_is_and_which_database_files_were_not_loaded() {
    let dir = genome_dir("skipped");
    let sa = dir.join("alpha").join("sa");
    std::fs::create_dir_all(&sa).unwrap();
    // A source in two shards, dated; an empty database, which is a build that found nothing; and a
    // truncated copy.
    write_osa2(&sa, "gnomad_chr1", "gnomad", "2026-11-02", 2);
    write_osa2(&sa, "gnomad_chr2", "gnomad", "2026-11-02", 2);
    write_osa2(&sa, "dbsnp", "dbsnp", "latest", 0);
    std::fs::write(
        sa.join("clinvar.osa2"),
        b"PK truncated in the middle of a copy",
    )
    .unwrap();
    let gff3 = dir.join("alpha").join("mini.gff3");
    let ctx =
        AnnotationContext::new(gff3.to_str(), None, sa.to_str(), 5000).expect("build context");
    let state = Arc::new(SharedState::new(
        ctx,
        default_genome_name(gff3.to_str(), Some(dir.as_path())),
        config(Some(dir), false),
    ));

    let (_, body) = get(&state, "/api/status").await;
    // One source, held in two files, named once and dated.
    assert_eq!(body["sa_sources"], json!(["gnomad"]));
    let info = body["sa_info"].as_array().unwrap();
    assert_eq!(info.len(), 1, "{info:?}");
    assert_eq!(info[0]["key"], "gnomad");
    assert_eq!(info[0]["version"], "2026-11-02");
    assert_eq!(info[0]["assembly"], "GRCh38");
    assert_eq!(info[0]["files"], 2);
    // The two files that did not load, by name and by class of reason, and nothing of the server's paths.
    let skipped = body["sa_skipped"].as_array().unwrap();
    let mut by_file: Vec<(String, String)> = skipped
        .iter()
        .map(|s| {
            (
                s["file"].as_str().unwrap().into(),
                s["reason"].as_str().unwrap().into(),
            )
        })
        .collect();
    by_file.sort();
    assert_eq!(
        by_file,
        [
            ("clinvar.osa2".to_string(), "unopenable".to_string()),
            ("dbsnp.osa2".to_string(), "empty".to_string()),
        ]
    );
    assert!(
        !body.to_string().contains(sa.to_str().unwrap()),
        "no path of the server is sent"
    );

    // The genome list says the same of the default, and of a genome that is not loaded.
    let (_, listing) = get(&state, "/api/genomes").await;
    let alpha = listing["genomes"]
        .as_array()
        .unwrap()
        .iter()
        .find(|g| g["name"] == "alpha")
        .unwrap();
    assert_eq!(alpha["sa_skipped"].as_array().unwrap().len(), 2);
    let beta = listing["genomes"]
        .as_array()
        .unwrap()
        .iter()
        .find(|g| g["name"] == "beta")
        .unwrap();
    assert_eq!(beta["sa_skipped"], json!([]));
}

#[tokio::test]
async fn the_genome_list_reports_the_skipped_files_of_a_genome_that_is_not_loaded() {
    let dir = genome_dir("skipped-unloaded");
    let sa = dir.join("beta").join("sa");
    std::fs::create_dir_all(&sa).unwrap();
    write_osa2(&sa, "dbsnp", "dbsnp", "latest", 0);
    let gff3 = dir.join("alpha").join("mini.gff3");
    let ctx = AnnotationContext::new(gff3.to_str(), None, None, 5000).expect("build context");
    let state = Arc::new(SharedState::new(
        ctx,
        default_genome_name(gff3.to_str(), Some(dir.as_path())),
        config(Some(dir), false),
    ));
    let (_, listing) = get(&state, "/api/genomes").await;
    let beta = listing["genomes"]
        .as_array()
        .unwrap()
        .iter()
        .find(|g| g["name"] == "beta")
        .unwrap();
    assert_eq!(beta["loaded"], false);
    assert_eq!(
        beta["sa_skipped"],
        json!([{ "file": "dbsnp.osa2", "reason": "empty" }])
    );
    // The empty database is not a source of the stack and is not counted as one.
    assert_eq!(beta["acmg"]["loaded"], json!([]));
    // And it is not an annotation: a directory whose only database is empty is listed as having
    // none, where it used to be listed as annotated and the page said so.
    assert_eq!(beta["has_sa"], false);
}

#[tokio::test]
async fn a_genome_is_listed_as_annotated_when_a_database_opened_whether_or_not_it_is_in_the_stack()
{
    let dir = genome_dir("has-sa-opened");
    let sa = dir.join("beta").join("sa");
    std::fs::create_dir_all(&sa).unwrap();
    // dbSNP is a database and not one of the nine sources the classifier reads.
    write_osa2(&sa, "dbsnp", "dbsnp", "latest", 2);
    let gff3 = dir.join("alpha").join("mini.gff3");
    let ctx = AnnotationContext::new(gff3.to_str(), None, None, 5000).expect("build context");
    let state = Arc::new(SharedState::new(
        ctx,
        default_genome_name(gff3.to_str(), Some(dir.as_path())),
        config(Some(dir), false),
    ));
    let (_, listing) = get(&state, "/api/genomes").await;
    let beta = listing["genomes"]
        .as_array()
        .unwrap()
        .iter()
        .find(|g| g["name"] == "beta")
        .unwrap();
    assert_eq!(beta["has_sa"], true);
    assert_eq!(beta["sa_skipped"], json!([]));
    assert_eq!(beta["acmg"]["loaded"], json!([]));
}
