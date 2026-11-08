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
    let state = multi_genome_state("options", false);
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
