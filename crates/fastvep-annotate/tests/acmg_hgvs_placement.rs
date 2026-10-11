//! The library and web path's copy of the placement regressions in
//! `fastvep-cli/tests/acmg_hgvs_placement.rs`.
//!
//! `AnnotationContext::annotate_vcf_text_with_acmg` carries its own per-variant
//! loop, separate from the CLI batch pipeline's, and both feed one classifier.
//! Nothing structurally keeps the two in step, and the placement is worked out
//! *after* `--pick` has deleted rows, by looking each survivor's transcript up
//! again - a step that exists once per loop. So each loop needs its own test of
//! the whole chain: parse the variant, predict, pick, place, classify.
//!
//! The fixture is the CLI test's: three forward-strand exons, exon 1 ending
//! `GTAT` and intron 1 opening `GTATGTAT`, so a `GTAT` duplication there is a
//! repeat that straddles the boundary.

use fastvep_annotate::pick::PickFlags;
use fastvep_annotate::AnnotationContext;
use fastvep_classification::AcmgConfig;
use serde_json::Value;
use std::fs::File;
use std::io::Write;
use std::path::Path;

fn reference() -> String {
    let exon1 = format!("ATG{}GCGTAT", "GCC".repeat(17));
    let intron1 = format!("GTATGTATCC{}CTAG", "C".repeat(86));
    let exon2 = format!("GAGGCC{}", "GCC".repeat(18));
    let intron2 = format!("GT{}AG", "C".repeat(96));
    let exon3 = format!("CTG{}TAA", "GCC".repeat(18));
    let seq = format!("{exon1}{intron1}{exon2}{intron2}{exon3}");
    assert_eq!(seq.len(), 380);
    seq
}

const GFF3: &str = "##gff-version 3\n\
1\ttest\tgene\t1\t380\t.\t+\t.\tID=gene:ENSG_PLC;Name=PLCTEST;gene_name=PLCTEST;biotype=protein_coding\n\
1\ttest\tmRNA\t1\t380\t.\t+\t.\tID=transcript:ENST_PLC;Parent=gene:ENSG_PLC;biotype=protein_coding\n\
1\ttest\texon\t1\t60\t.\t+\t.\tID=exon:E1;Parent=transcript:ENST_PLC;rank=1\n\
1\ttest\texon\t161\t220\t.\t+\t.\tID=exon:E2;Parent=transcript:ENST_PLC;rank=2\n\
1\ttest\texon\t321\t380\t.\t+\t.\tID=exon:E3;Parent=transcript:ENST_PLC;rank=3\n\
1\ttest\tCDS\t1\t60\t.\t+\t0\tID=CDS:ENSP_PLC;Parent=transcript:ENST_PLC\n\
1\ttest\tCDS\t161\t220\t.\t+\t0\tID=CDS:ENSP_PLC;Parent=transcript:ENST_PLC\n\
1\ttest\tCDS\t321\t380\t.\t+\t0\tID=CDS:ENSP_PLC;Parent=transcript:ENST_PLC\n";

/// A curated loss-of-function mechanism is what lets PVS1 run without any
/// constraint or disease-gene database on disk.
const ACMG_CONFIG: &str = "[gene_overrides.PLCTEST]\nmechanism = \"LOF\"\n";

fn write(path: &Path, contents: &str) {
    File::create(path)
        .unwrap()
        .write_all(contents.as_bytes())
        .unwrap();
}

/// Annotate `records` (`POS REF ALT` rows) through the library entry point.
fn annotate(records: &[(u64, String, String)]) -> Vec<Value> {
    annotate_over(
        GFF3,
        ACMG_CONFIG,
        PickFlags::default(),
        &vcf_text(records, None),
    )
}

/// The VCF the fixture's tests send. With `proband_gt` the records carry a FORMAT
/// column and one sample, `P`, holding that genotype.
fn vcf_text(records: &[(u64, String, String)], proband_gt: Option<&str>) -> String {
    let mut vcf =
        String::from("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
    if proband_gt.is_some() {
        vcf.push_str("\tFORMAT\tP");
    }
    vcf.push('\n');
    for (pos, r, a) in records {
        vcf.push_str(&format!("1\t{pos}\t.\t{r}\t{a}\t.\tPASS\t."));
        if let Some(gt) = proband_gt {
            vcf.push_str(&format!("\tGT:DP:GQ\t{gt}:30:99"));
        }
        vcf.push('\n');
    }
    vcf
}

/// Annotate a VCF over `gff3` with the ACMG config `acmg_toml` and the pick `pick`.
fn annotate_over(gff3_text: &str, acmg_toml: &str, pick: PickFlags, vcf: &str) -> Vec<Value> {
    let tmp = tempfile::tempdir().unwrap();
    let gff3 = tmp.path().join("test.gff3");
    let fasta = tmp.path().join("test.fa");
    let config = tmp.path().join("acmg.toml");
    write(&gff3, gff3_text);
    write(&fasta, &format!(">1\n{}\n", reference()));
    write(&config, acmg_toml);

    let mut ctx = AnnotationContext::new(
        Some(gff3.to_str().unwrap()),
        Some(fasta.to_str().unwrap()),
        None,
        5000,
    )
    .expect("the fixture should load");
    ctx.hgvs = true;
    let acmg = AcmgConfig::from_toml_file(config.to_str().unwrap()).unwrap();

    ctx.annotate_vcf_text_with_acmg(vcf, pick.resolve(), Some(&acmg))
        .expect("annotation should succeed")
}

fn transcript(variant: &Value) -> &Value {
    let all = variant["transcript_consequences"].as_array().unwrap();
    assert_eq!(
        all.len(),
        1,
        "one transcript and allele expected: {variant}"
    );
    &all[0]
}

fn terms(variant: &Value) -> Vec<String> {
    terms_of(transcript(variant))
}

fn terms_of(consequence: &Value) -> Vec<String> {
    consequence["consequence_terms"]
        .as_array()
        .unwrap()
        .iter()
        .map(|t| t.as_str().unwrap().to_string())
        .collect()
}

fn criterion<'a>(variant: &'a Value, code: &str) -> &'a Value {
    transcript(variant)["acmg"]["criteria"]
        .as_array()
        .unwrap()
        .iter()
        .find(|c| c["code"] == code || c["code"].as_str().is_some_and(|s| s.starts_with(code)))
        .unwrap_or_else(|| panic!("no {code} criterion in {variant}"))
}

fn met(variant: &Value, code: &str) -> bool {
    criterion(variant, code)["met"].as_bool().unwrap()
}

/// MSH6 `c.3801+3_3801+6dup` in miniature: the record anchors on the exon's last
/// base, the term is `frameshift_variant`, and HGVS places four intronic bases.
#[test]
fn lib_pvs1_stands_down_for_a_frameshift_term_that_hgvs_places_in_the_intron() {
    let seq = reference();
    let anchor = &seq[59..60];
    let out = annotate(&[(60, anchor.to_string(), format!("{anchor}GTAT"))]);
    let v = &out[0];
    assert!(
        terms(v).contains(&"frameshift_variant".to_string()),
        "{:?}",
        terms(v)
    );
    assert!(
        transcript(v)["hgvsc"]
            .as_str()
            .unwrap()
            .contains("c.60+5_60+8dup"),
        "{}",
        transcript(v)["hgvsc"]
    );
    assert!(!met(v, "PVS1"), "{}", criterion(v, "PVS1"));
}

/// NBEAL2 shape: the same duplication anchored inside the exon.
#[test]
fn lib_pvs1_stands_down_for_an_exon_anchored_duplication_hgvs_places_in_the_intron() {
    let seq = reference();
    let anchor = &seq[55..56];
    let out = annotate(&[(56, anchor.to_string(), format!("{anchor}GTAT"))]);
    let v = &out[0];
    assert!(
        terms(v).contains(&"frameshift_variant".to_string()),
        "{:?}",
        terms(v)
    );
    assert!(!met(v, "PVS1"), "{}", criterion(v, "PVS1"));
    assert_eq!(
        criterion(v, "PVS1")["details"]["hgvs_placement"],
        serde_json::json!({"kind": "intronic", "offset": 8})
    );
}

#[test]
fn lib_pvs1_still_fires_for_a_frameshift_that_is_in_the_exon() {
    let seq = reference();
    let out = annotate(&[(189, seq[188..190].to_string(), seq[188..189].to_string())]);
    let v = &out[0];
    assert!(
        terms(v).contains(&"frameshift_variant".to_string()),
        "{:?}",
        terms(v)
    );
    assert!(met(v, "PVS1"), "{}", criterion(v, "PVS1"));
}

/// The `TAGG>T` of the benchmark's MLH1, RYR1, ATM, SORD and RYR2 rows:
/// `splice_acceptor_variant` by anchor, an in-frame exonic deletion once HGVS has
/// placed it.
#[test]
fn lib_pvs1_stands_down_for_an_acceptor_deletion_hgvs_places_in_frame_in_the_exon() {
    let seq = reference();
    let out = annotate(&[(158, seq[157..161].to_string(), seq[157..158].to_string())]);
    let v = &out[0];
    assert!(
        terms(v).contains(&"splice_acceptor_variant".to_string()),
        "{:?}",
        terms(v)
    );
    assert!(!met(v, "PVS1"), "{}", criterion(v, "PVS1"));
}

#[test]
fn lib_pvs1_still_fires_when_the_acceptor_dinucleotide_is_really_destroyed() {
    let seq = reference();
    let out = annotate(&[(158, seq[157..160].to_string(), seq[157..158].to_string())]);
    let v = &out[0];
    assert!(
        terms(v).contains(&"splice_acceptor_variant".to_string()),
        "{:?}",
        terms(v)
    );
    assert!(met(v, "PVS1"), "{}", criterion(v, "PVS1"));
}

/// Walker 2023's BP7 exon-edge exclusion reaches the criterion on this path too.
#[test]
fn lib_bp7_is_excluded_at_the_first_base_and_last_three_bases_of_an_exon() {
    let seq = reference();
    let rec = |pos: usize, alt: &str| (pos as u64, seq[pos - 1..pos].to_string(), alt.to_string());
    let out = annotate(&[rec(190, "A"), rec(217, "A"), rec(220, "A"), rec(321, "T")]);
    let edge = |i: usize| criterion(&out[i], "BP7")["details"]["at_exon_edge"].clone();
    assert_eq!(edge(0), Value::Bool(false), "mid-exon");
    assert_eq!(edge(1), Value::Bool(false), "fourth base from the end");
    assert_eq!(edge(2), Value::Bool(true), "last base of an exon");
    assert_eq!(edge(3), Value::Bool(true), "first base of an exon");
}

/// Two models of one gene over the same reference. `ENST_A` is `GFF3`'s, whose
/// exon 1 ends at base 60. `ENST_B`'s exon 1 runs on to base 81, so the `GTAT`
/// repeat that straddles A's exon 1 boundary lies wholly inside B's exon: one
/// record is placed in the intron by A and in the exon, as a frameshift, by B.
/// `canonical` is the model carrying `Ensembl_canonical`, which `--pick` keeps.
fn two_model_gff3(canonical: &str) -> String {
    let model = |id: &str, exon1_end: u32| {
        let tag = if id == canonical {
            ";tag=Ensembl_canonical"
        } else {
            ""
        };
        format!(
            "1\ttest\tmRNA\t1\t380\t.\t+\t.\tID=transcript:ENST_{id};Parent=gene:ENSG_PLC;biotype=protein_coding{tag}\n\
             1\ttest\texon\t1\t{exon1_end}\t.\t+\t.\tID=exon:{id}1;Parent=transcript:ENST_{id};rank=1\n\
             1\ttest\texon\t161\t220\t.\t+\t.\tID=exon:{id}2;Parent=transcript:ENST_{id};rank=2\n\
             1\ttest\texon\t321\t380\t.\t+\t.\tID=exon:{id}3;Parent=transcript:ENST_{id};rank=3\n\
             1\ttest\tCDS\t1\t{exon1_end}\t.\t+\t0\tID=CDS:ENSP_{id};Parent=transcript:ENST_{id}\n\
             1\ttest\tCDS\t161\t220\t.\t+\t0\tID=CDS:ENSP_{id};Parent=transcript:ENST_{id}\n\
             1\ttest\tCDS\t321\t380\t.\t+\t0\tID=CDS:ENSP_{id};Parent=transcript:ENST_{id}\n"
        )
    };
    format!(
        "##gff-version 3\n\
         1\ttest\tgene\t1\t380\t.\t+\t.\tID=gene:ENSG_PLC;Name=PLCTEST;gene_name=PLCTEST;biotype=protein_coding\n\
         {}{}",
        model("A", 60),
        model("B", 81)
    )
}

/// The NBEAL2-shape record the two-model tests share: `GTAT` inserted before the
/// first base of the repeat. Model A places it in the intron, past the donor, so
/// PVS1 stands down; model B's exon 1 runs on past the whole repeat, so B places
/// it in the exon, out of frame, and PVS1 stands.
fn repeat_insertion() -> (u64, String, String) {
    let seq = reference();
    let anchor = &seq[55..56];
    (56, anchor.to_string(), format!("{anchor}GTAT"))
}

/// The consequence a variant reports for one model, whichever other rows it has.
fn model<'a>(variant: &'a Value, id: &str) -> &'a Value {
    variant["transcript_consequences"]
        .as_array()
        .unwrap()
        .iter()
        .find(|tc| {
            tc["transcript_id"]
                .as_str()
                .is_some_and(|t| t.starts_with(&format!("ENST_{id}")))
        })
        .unwrap_or_else(|| panic!("no consequence for ENST_{id} in {variant}"))
}

fn pvs1_of(consequence: &Value) -> &Value {
    consequence["acmg"]["criteria"]
        .as_array()
        .unwrap()
        .iter()
        .find(|c| c["code"].as_str().is_some_and(|s| s.starts_with("PVS1")))
        .unwrap_or_else(|| panic!("no PVS1 criterion in {consequence}"))
}

/// Placement is worked out after `--pick` has deleted rows, by looking each
/// survivor's own transcript up again. A survivor that was handed another
/// model's exons would still pass every test over a one-transcript fixture, so
/// this keeps the *second* of two models and asks which one the criteria read.
#[test]
fn lib_pick_keeping_the_second_model_places_the_change_against_that_model() {
    let pick = PickFlags {
        pick: true,
        ..Default::default()
    };
    let vcf = vcf_text(&[repeat_insertion()], None);
    let out = annotate_over(&two_model_gff3("B"), ACMG_CONFIG, pick, &vcf);
    let rows = out[0]["transcript_consequences"].as_array().unwrap();
    assert_eq!(rows.len(), 1, "--pick keeps one row: {}", out[0]);
    let kept = model(&out[0], "B");

    assert!(
        terms_of(kept).contains(&"frameshift_variant".to_string()),
        "{:?}",
        terms_of(kept)
    );
    let pvs1 = pvs1_of(kept);
    assert_eq!(
        pvs1["details"]["hgvs_placement"],
        serde_json::json!({"kind": "exonic", "in_frame": false}),
        "B's exon 1 holds the whole repeat, so B places the change in its exon: {pvs1}"
    );
    assert!(pvs1["met"].as_bool().unwrap(), "{pvs1}");
}

/// The mirror, with the first model kept: its own exon 1 ends inside the repeat.
#[test]
fn lib_pick_keeping_the_first_model_places_the_change_against_that_model() {
    let pick = PickFlags {
        pick: true,
        ..Default::default()
    };
    let vcf = vcf_text(&[repeat_insertion()], None);
    let out = annotate_over(&two_model_gff3("A"), ACMG_CONFIG, pick, &vcf);
    let rows = out[0]["transcript_consequences"].as_array().unwrap();
    assert_eq!(rows.len(), 1, "--pick keeps one row: {}", out[0]);
    let pvs1 = pvs1_of(model(&out[0], "A"));
    assert_eq!(
        pvs1["details"]["hgvs_placement"],
        serde_json::json!({"kind": "intronic", "offset": 8})
    );
    assert!(!pvs1["met"].as_bool().unwrap(), "{pvs1}");
}

/// `--flag-pick` keeps every row and marks one, so each row is classified in
/// place and has to be paired with its own model by position.
#[test]
fn lib_flag_pick_places_each_row_against_its_own_model() {
    let pick = PickFlags {
        flag_pick: true,
        ..Default::default()
    };
    let vcf = vcf_text(&[repeat_insertion()], None);
    let out = annotate_over(&two_model_gff3("B"), ACMG_CONFIG, pick, &vcf);
    let rows = out[0]["transcript_consequences"].as_array().unwrap();
    assert_eq!(rows.len(), 2, "--flag-pick keeps both rows: {}", out[0]);

    let (a, b) = (model(&out[0], "A"), model(&out[0], "B"));
    assert_eq!(b["pick"], 1, "B is the canonical model: {b}");
    assert!(a.get("pick").is_none_or(|p| p != 1), "{a}");
    assert_eq!(
        pvs1_of(a)["details"]["hgvs_placement"],
        serde_json::json!({"kind": "intronic", "offset": 8})
    );
    assert!(!pvs1_of(a)["met"].as_bool().unwrap(), "{}", pvs1_of(a));
    assert_eq!(
        pvs1_of(b)["details"]["hgvs_placement"],
        serde_json::json!({"kind": "exonic", "in_frame": false})
    );
    assert!(pvs1_of(b)["met"].as_bool().unwrap(), "{}", pvs1_of(b));
}

/// A deletion that removes a whole exon with both of its splice sites is not
/// placed anywhere that spares the exon, so there is nothing to stand down for:
/// PVS1 is decided by the record's terms, as it was before there was a
/// placement. The span helper was tested alone; nothing drove this through the
/// pipeline.
#[test]
fn lib_pvs1_still_fires_for_a_deletion_that_removes_a_whole_exon() {
    let seq = reference();
    // POS 158 deletes 159..222: exon 2 (161..220) with the acceptor `AG` before
    // it and the donor `GT` after it.
    let out = annotate(&[(158, seq[157..222].to_string(), seq[157..158].to_string())]);
    let v = &out[0];
    let pvs1 = criterion(v, "PVS1");
    assert_eq!(
        pvs1["details"]["hgvs_placement"],
        serde_json::json!({"kind": "spans"}),
        "{pvs1}"
    );
    assert!(
        !pvs1["summary"]
            .as_str()
            .unwrap()
            .contains("HGVS places the change"),
        "PVS1 stood down for a deletion that spans an exon: {pvs1}"
    );
    assert!(met(v, "PVS1"), "{pvs1}");
}

/// The compound-heterozygote pass classifies a variant a second time once its
/// companions are known, from the splice context stored on the allele. A pass
/// that rebuilt the input without it would hand the placed variant back its
/// PVS1, so this runs two heterozygous variants of one gene and checks both that
/// the second pass happened and that the placement survived it.
#[test]
fn lib_compound_het_second_pass_keeps_the_placement() {
    let seq = reference();
    // Exon 2 holds the second variant, a missense at its fifth codon.
    let missense = (170u64, seq[169..170].to_string(), "A".to_string());
    let vcf = vcf_text(&[repeat_insertion(), missense], Some("0/1"));
    let config = format!("{ACMG_CONFIG}[trio]\nproband = \"P\"\n");
    let out = annotate_over(GFF3, &config, PickFlags::default(), &vcf);

    let dup = &out[0];
    let bp2 = criterion(dup, "BP2");
    assert!(
        !bp2["summary"]
            .as_str()
            .unwrap()
            .contains("No companion variants available"),
        "the companion pass did not run for this variant: {bp2}"
    );
    let pvs1 = criterion(dup, "PVS1");
    assert!(!met(dup, "PVS1"), "{pvs1}");
    assert_eq!(
        pvs1["details"]["hgvs_placement"],
        serde_json::json!({"kind": "intronic", "offset": 8})
    );
}
