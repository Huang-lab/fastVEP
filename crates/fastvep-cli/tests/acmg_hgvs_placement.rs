//! End-to-end regressions for criteria that read *where the change is*, as
//! opposed to *how the record was written*.
//!
//! A consequence term is decided from the position the VCF record anchors the
//! change at, and an indel in a repeat can be written at several positions that
//! all produce the same edited sequence. HGVS numbers it at the most 3' one.
//! The round-6 medical-genetics review found three places where the two
//! disagreed and the ACMG criteria followed the record instead of the sequence:
//!
//! * MSH6 `c.3801+3_3801+6dup`: the record anchors on the last base of exon 8,
//!   so the term is `frameshift_variant`, and the change is four intronic bases.
//! * An acceptor-straddling deletion, the `TAGG>T` of the benchmark's MLH1, RYR1,
//!   ATM, SORD and RYR2 rows (FOXP2's `CCAG>C`, one unit of a CAG repeat that
//!   runs across the acceptor, reaches the same place): the record deletes the
//!   `AG` and the first exon base, so the term is `splice_acceptor_variant`, and
//!   HGVS writes an in-frame deletion wholly inside the exon.
//! * BP7 on a synonymous change at the first base or last three bases of an
//!   exon, which Walker 2023 excludes and the pipeline never told the criterion
//!   about.
//!
//! These drive the real `run_annotate` entry point over a GFF3 + FASTA + VCF on
//! disk, so they fail if any link between parsing the variant and writing the
//! classification regresses. The CLI pipeline carries its own copy of the
//! per-allele loop, separate from the one in `fastvep-annotate`; nothing
//! structurally keeps the two in step, so each needs its own test.

use fastvep_cli::pipeline::{run_annotate, AnnotateConfig, PickFlags};
use serde_json::Value;
use std::fs::File;
use std::io::Write;
use std::path::Path;

/// Three exons on contig `1`, forward strand, the whole of each coding.
///
/// ```text
///   exon 1   1..60    ATG GCC x17 GCG TAT          (60 bases, ends GTAT)
///   intron 1 61..160  GTATGTATCC C... CTAG          (opens with GTAT x2, ends in AG)
///   exon 2   161..220 GAG GCC GCC x18               (60 bases, opens GAGGCC)
///   intron 2 221..320 GT C... AG
///   exon 3   321..380 CTG GCC x18 TAA               (60 bases, ends in the stop codon)
/// ```
///
/// Exon 1 ends `GTAT` and intron 1 opens `GTATGTAT`, so a `GTAT` duplication
/// there is a repeat that straddles the boundary and can be anchored in the exon
/// or on its last base. Every other codon is an alanine bar the first of each exon
/// and the stop, so no change below lands on an accidental stop.
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

/// Annotate `records` (`POS REF ALT` rows) and return the parsed JSON array.
fn annotate(records: &[(u64, String, String)]) -> Vec<Value> {
    annotate_over(
        GFF3,
        ACMG_CONFIG,
        PickFlags::default(),
        None,
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

/// Annotate a VCF over `gff3` with the ACMG config `acmg_toml`, the pick `pick`
/// and, when given, `proband` as the trio's proband.
fn annotate_over(
    gff3_text: &str,
    acmg_toml: &str,
    pick: PickFlags,
    proband: Option<&str>,
    vcf: &str,
) -> Vec<Value> {
    let tmp = tempfile::tempdir().unwrap();
    let gff3 = tmp.path().join("test.gff3");
    let fasta = tmp.path().join("test.fa");
    let vcf_in = tmp.path().join("in.vcf");
    let out = tmp.path().join("out.json");
    let config = tmp.path().join("acmg.toml");
    write(&gff3, gff3_text);
    write(&fasta, &format!(">1\n{}\n", reference()));
    write(&config, acmg_toml);
    write(&vcf_in, vcf);

    run_annotate(AnnotateConfig {
        input: vcf_in.to_string_lossy().into(),
        output: out.to_string_lossy().into(),
        gff3: vec![gff3.to_string_lossy().into()],
        fasta: Some(fasta.to_string_lossy().into()),
        output_format: "json".into(),
        pick,
        hgvs: true,
        distance: 5000,
        cache_dir: None,
        transcript_cache: None,
        sa_dir: None,
        sa_only: false,
        acmg: true,
        acmg_config: Some(config.to_string_lossy().into()),
        pick_order: None,
        functional_evidence: None,
        proband: proband.map(str::to_string),
        mother: None,
        father: None,
        gene_list: None,
        explicit_alleles: false,
        qc_rules: None,
        show_progress: false,
        compat: Default::default(),
        columns: Default::default(),
    })
    .expect("annotation should succeed");

    let text = std::fs::read_to_string(&out).unwrap();
    serde_json::from_str(&text).expect("the output should be a JSON array")
}

/// The transcript consequence for the only transcript and allele, as JSON.
fn transcript(variant: &Value) -> &Value {
    let all = variant["transcript_consequences"].as_array().unwrap();
    assert_eq!(
        all.len(),
        1,
        "one transcript and allele expected: {variant}"
    );
    &all[0]
}

/// The consequence of one ALT of a multi-allelic record, as JSON.
fn allele<'a>(variant: &'a Value, alt: &str) -> &'a Value {
    variant["transcript_consequences"]
        .as_array()
        .unwrap()
        .iter()
        .find(|tc| tc["variant_allele"] == alt)
        .unwrap_or_else(|| panic!("no consequence for {alt} in {variant}"))
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

fn criterion_of<'a>(consequence: &'a Value, code: &str) -> &'a Value {
    consequence["acmg"]["criteria"]
        .as_array()
        .unwrap()
        .iter()
        .find(|c| c["code"] == code || c["code"].as_str().is_some_and(|s| s.starts_with(code)))
        .unwrap_or_else(|| panic!("no {code} criterion in {consequence}"))
}

fn criterion<'a>(variant: &'a Value, code: &str) -> &'a Value {
    criterion_of(transcript(variant), code)
}

fn met(variant: &Value, code: &str) -> bool {
    criterion(variant, code)["met"].as_bool().unwrap()
}

/// MSH6 `c.3801+3_3801+6dup`, in miniature. The record inserts `GTAT` after the
/// last base of exon 1, where the intron already opens `GTAT`; the term is
/// `frameshift_variant` and HGVS writes a duplication of intronic bases, so the
/// exon is untouched and the donor `GT` is intact. PVS1 must not fire.
#[test]
fn cli_pvs1_stands_down_for_a_frameshift_term_that_hgvs_places_in_the_intron() {
    let seq = reference();
    let anchor = &seq[59..60];
    let out = annotate(&[(60, anchor.to_string(), format!("{anchor}GTAT"))]);
    let v = &out[0];

    // Guard the premise: both halves of the disagreement must be present, or
    // the assertions below pass for the wrong reason.
    assert!(
        terms(v).contains(&"frameshift_variant".to_string()),
        "the record should read as a frameshift, got {:?}",
        terms(v)
    );
    let hgvsc = transcript(v)["hgvsc"].as_str().unwrap();
    assert!(
        hgvsc.contains("c.60+5_60+8dup"),
        "HGVS should number it as a duplication four bases into the intron, got {hgvsc}"
    );

    assert!(
        !met(v, "PVS1"),
        "PVS1 fired on a frameshift term that HGVS places wholly in the intron: {}",
        criterion(v, "PVS1")
    );
}

/// NBEAL2 `c.6560+4_6560+20dup` arrives anchored *inside* the exon: the caller
/// left-aligned a repeat that straddles the boundary, so the record reaches no
/// intronic base the way it is written. Nothing about the record says where HGVS
/// will put it, which is why the placement walks the reference instead of reading
/// the position.
#[test]
fn cli_pvs1_stands_down_for_an_exon_anchored_duplication_hgvs_places_in_the_intron() {
    let seq = reference();
    // Exon 1 ends GTAT at 57..60: insert GTAT before its first base.
    let anchor = &seq[55..56];
    let out = annotate(&[(56, anchor.to_string(), format!("{anchor}GTAT"))]);
    let v = &out[0];

    assert!(
        terms(v).contains(&"frameshift_variant".to_string()),
        "the record should read as a frameshift, got {:?}",
        terms(v)
    );
    assert!(
        terms(v).iter().all(|t| !t.contains("intron")),
        "the record reaches no intronic base as written: {:?}",
        terms(v)
    );
    let hgvsc = transcript(v)["hgvsc"].as_str().unwrap();
    assert!(hgvsc.contains("c.60+5_60+8dup"), "{hgvsc}");

    assert!(!met(v, "PVS1"), "{}", criterion(v, "PVS1"));
    assert_eq!(
        criterion(v, "PVS1")["details"]["hgvs_placement"],
        serde_json::json!({"kind": "intronic", "offset": 8})
    );
}

/// A multi-allelic record is classified allele by allele. The first ALT is the
/// MSH6-shape duplication; the second is a stop-gain at the same base, which is
/// as null as a variant gets. One must not be decided by the other.
#[test]
fn cli_placement_is_per_allele_on_a_multi_allelic_record() {
    let seq = reference();
    assert_eq!(&seq[59..60], "T");
    let out = annotate(&[(60, "T".into(), "TGTAT,A".into())]);
    let v = &out[0];

    let dup = allele(v, "TGTAT");
    assert!(terms_of(dup).contains(&"frameshift_variant".to_string()));
    assert!(!criterion_of(dup, "PVS1")["met"].as_bool().unwrap());

    let stop = allele(v, "A");
    assert!(
        terms_of(stop).contains(&"stop_gained".to_string()),
        "{:?}",
        terms_of(stop)
    );
    assert!(criterion_of(stop, "PVS1")["met"].as_bool().unwrap());
}

/// A deletion that reaches both an exon and its intron and has nowhere to slide
/// is not placed anywhere that spares either, so there is nothing to stand down
/// for: PVS1 is decided by the record's terms, as it was before there was a
/// placement.
#[test]
fn cli_pvs1_still_fires_for_a_deletion_spanning_the_exon_and_its_donor() {
    let seq = reference();
    // Exon 2 ends ...GCC and intron 2 opens GT. Delete 219..221 (CCG): the base
    // after the deletion, T, is not the C it starts with, so it cannot slide.
    let out = annotate(&[(218, seq[217..221].to_string(), seq[217..218].to_string())]);
    let v = &out[0];
    assert!(
        terms(v).contains(&"splice_donor_variant".to_string()),
        "{:?}",
        terms(v)
    );
    assert!(met(v, "PVS1"), "{}", criterion(v, "PVS1"));
    assert_eq!(
        criterion(v, "PVS1")["details"]["hgvs_placement"],
        serde_json::json!({"kind": "spans"})
    );
}

/// The converse, so the gate cannot pass by never firing: a one-base deletion
/// in the middle of an exon is a frameshift wherever it is written, and PVS1
/// stays.
#[test]
fn cli_pvs1_still_fires_for_a_frameshift_that_is_in_the_exon() {
    let seq = reference();
    // Delete base 190 of exon 2, anchored on the base before it: REF is two
    // bases, ALT is one.
    let out = annotate(&[(189, seq[188..190].to_string(), seq[188..189].to_string())]);
    let v = &out[0];
    assert!(
        terms(v).contains(&"frameshift_variant".to_string()),
        "{:?}",
        terms(v)
    );
    assert!(met(v, "PVS1"), "{}", criterion(v, "PVS1"));
}

/// The `TAGG>T` shape at an acceptor, as in the benchmark's MLH1, RYR1, ATM, SORD
/// and RYR2 rows. The record deletes the intron's final `AG` and the exon's first
/// base, so the term is `splice_acceptor_variant`; but exon 2 opens `GAGG`, so
/// the same sequence is the exon's `AGG` deleted, in frame, with the acceptor
/// `AG` untouched, and HGVS writes it that way. Not a null variant. FOXP2's
/// `CCAG>C` is the same class by another route: one unit of a CAG repeat that
/// runs across the acceptor, which HGVS writes `c.618_620del` and the reviewer
/// read as "Q209del".
#[test]
fn cli_pvs1_stands_down_for_an_acceptor_deletion_hgvs_places_in_frame_in_the_exon() {
    let seq = reference();
    // POS 158 REF TAGG ALT T deletes 159..161.
    let out = annotate(&[(158, seq[157..161].to_string(), seq[157..158].to_string())]);
    let v = &out[0];

    assert!(
        terms(v).contains(&"splice_acceptor_variant".to_string()),
        "the record should read as a splice acceptor variant, got {:?}",
        terms(v)
    );
    let hgvsc = transcript(v)["hgvsc"].as_str().unwrap();
    assert!(
        hgvsc.ends_with("c.62_64del"),
        "HGVS should number it as an in-frame deletion inside exon 2, got {hgvsc}"
    );

    assert!(
        !met(v, "PVS1"),
        "PVS1 fired on an in-frame exonic deletion that only the record's anchoring calls a splice variant: {}",
        criterion(v, "PVS1")
    );
}

/// The same acceptor, destroyed: deleting the intron's final `AG` where nothing
/// repeats it cannot be rewritten anywhere that spares the dinucleotide, so it
/// is the null variant PVS1 is for.
#[test]
fn cli_pvs1_still_fires_when_the_acceptor_dinucleotide_is_really_destroyed() {
    let seq = reference();
    // POS 158 REF TAG ALT T deletes 159..160, the AG itself. The base after it
    // is the exon's G, not an A, so the deletion has nowhere to slide.
    let out = annotate(&[(158, seq[157..160].to_string(), seq[157..158].to_string())]);
    let v = &out[0];
    assert!(
        terms(v).contains(&"splice_acceptor_variant".to_string()),
        "{:?}",
        terms(v)
    );
    let hgvsc = transcript(v)["hgvsc"].as_str().unwrap();
    assert!(hgvsc.contains("c.61-2_61-1del"), "{hgvsc}");
    assert!(met(v, "PVS1"), "{}", criterion(v, "PVS1"));
}

/// Walker 2023: BP7 must not be applied to a synonymous change at the first
/// base or the last three bases of an exon. The pipeline never populated the
/// field the criterion reads, so the exclusion that `docs/ACMG.md` describes did
/// not exist: CDH23 p.Thr2454= (ClinVar Pathogenic) and COL7A1 p.Pro1699=
/// (Likely pathogenic) both collected BP7 at an exon's last bases.
#[test]
fn cli_bp7_is_excluded_at_the_first_base_and_last_three_bases_of_an_exon() {
    let seq = reference();
    let rec = |pos: usize, alt: &str| (pos as u64, seq[pos - 1..pos].to_string(), alt.to_string());
    // 220 is the last base of exon 2 (third base of an alanine codon); 321 is
    // the first base of exon 3 (CTG -> TTG, leucine); 217 is the fourth base
    // from the end of exon 2 and so is not at the edge; 190 is mid-exon.
    let out = annotate(&[rec(190, "A"), rec(217, "A"), rec(220, "A"), rec(321, "T")]);
    let edge = |i: usize| {
        let c = criterion(&out[i], "BP7");
        c["details"]["at_exon_edge"].clone()
    };
    for (i, label) in [
        "mid-exon",
        "fourth from the end",
        "last base",
        "first base of exon 3",
    ]
    .iter()
    .enumerate()
    {
        assert!(
            terms(&out[i]).contains(&"synonymous_variant".to_string()),
            "{label}: expected a synonymous change, got {:?}",
            terms(&out[i])
        );
    }
    assert_eq!(edge(0), Value::Bool(false), "mid-exon is not an exon edge");
    assert_eq!(
        edge(1),
        Value::Bool(false),
        "the fourth base from the end is not an exon edge"
    );
    assert_eq!(
        edge(2),
        Value::Bool(true),
        "the last base of an exon is an exon edge"
    );
    assert_eq!(
        edge(3),
        Value::Bool(true),
        "the first base of an exon is an exon edge"
    );
    assert!(!met(&out[2], "BP7") && !met(&out[3], "BP7"));
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

/// Placement is worked out after `--pick` has deleted rows, by looking each
/// survivor's own transcript up again, in a loop of the CLI's own. A survivor
/// that was handed another model's exons would still pass every test over a
/// one-transcript fixture, so this keeps the *second* of two models and asks
/// which one the criteria read.
#[test]
fn cli_pick_keeping_the_second_model_places_the_change_against_that_model() {
    let pick = PickFlags {
        pick: true,
        ..Default::default()
    };
    let vcf = vcf_text(&[repeat_insertion()], None);
    let out = annotate_over(&two_model_gff3("B"), ACMG_CONFIG, pick, None, &vcf);
    let rows = out[0]["transcript_consequences"].as_array().unwrap();
    assert_eq!(rows.len(), 1, "--pick keeps one row: {}", out[0]);
    let kept = model(&out[0], "B");

    assert!(
        terms_of(kept).contains(&"frameshift_variant".to_string()),
        "{:?}",
        terms_of(kept)
    );
    let pvs1 = criterion_of(kept, "PVS1");
    assert_eq!(
        pvs1["details"]["hgvs_placement"],
        serde_json::json!({"kind": "exonic", "in_frame": false}),
        "B's exon 1 holds the whole repeat, so B places the change in its exon: {pvs1}"
    );
    assert!(pvs1["met"].as_bool().unwrap(), "{pvs1}");
}

/// The mirror, with the first model kept: its own exon 1 ends inside the repeat.
#[test]
fn cli_pick_keeping_the_first_model_places_the_change_against_that_model() {
    let pick = PickFlags {
        pick: true,
        ..Default::default()
    };
    let vcf = vcf_text(&[repeat_insertion()], None);
    let out = annotate_over(&two_model_gff3("A"), ACMG_CONFIG, pick, None, &vcf);
    let rows = out[0]["transcript_consequences"].as_array().unwrap();
    assert_eq!(rows.len(), 1, "--pick keeps one row: {}", out[0]);
    let pvs1 = criterion_of(model(&out[0], "A"), "PVS1");
    assert_eq!(
        pvs1["details"]["hgvs_placement"],
        serde_json::json!({"kind": "intronic", "offset": 8})
    );
    assert!(!pvs1["met"].as_bool().unwrap(), "{pvs1}");
}

/// `--flag-pick` keeps every row and marks one, so each row is classified in
/// place and has to be paired with its own model by position.
#[test]
fn cli_flag_pick_places_each_row_against_its_own_model() {
    let pick = PickFlags {
        flag_pick: true,
        ..Default::default()
    };
    let vcf = vcf_text(&[repeat_insertion()], None);
    let out = annotate_over(&two_model_gff3("B"), ACMG_CONFIG, pick, None, &vcf);
    let rows = out[0]["transcript_consequences"].as_array().unwrap();
    assert_eq!(rows.len(), 2, "--flag-pick keeps both rows: {}", out[0]);

    let (a, b) = (model(&out[0], "A"), model(&out[0], "B"));
    assert_eq!(b["pick"], 1, "B is the canonical model: {b}");
    assert!(a.get("pick").is_none_or(|p| p != 1), "{a}");
    let (pvs1_a, pvs1_b) = (criterion_of(a, "PVS1"), criterion_of(b, "PVS1"));
    assert_eq!(
        pvs1_a["details"]["hgvs_placement"],
        serde_json::json!({"kind": "intronic", "offset": 8})
    );
    assert!(!pvs1_a["met"].as_bool().unwrap(), "{pvs1_a}");
    assert_eq!(
        pvs1_b["details"]["hgvs_placement"],
        serde_json::json!({"kind": "exonic", "in_frame": false})
    );
    assert!(pvs1_b["met"].as_bool().unwrap(), "{pvs1_b}");
}

/// A deletion that removes a whole exon with both of its splice sites is not
/// placed anywhere that spares the exon, so there is nothing to stand down for:
/// PVS1 is decided by the record's terms, as it was before there was a
/// placement. The span helper was tested alone; nothing drove this through the
/// pipeline.
#[test]
fn cli_pvs1_still_fires_for_a_deletion_that_removes_a_whole_exon() {
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
/// companions are known, from the splice context stored on the allele, in a loop
/// of the CLI's own. A pass that rebuilt the input without it would hand the
/// placed variant back its PVS1, so this runs two heterozygous variants of one
/// gene and checks both that the second pass happened and that the placement
/// survived it.
#[test]
fn cli_compound_het_second_pass_keeps_the_placement() {
    let seq = reference();
    // Exon 2 holds the second variant, a missense at its fifth codon.
    let missense = (170u64, seq[169..170].to_string(), "A".to_string());
    let vcf = vcf_text(&[repeat_insertion(), missense], Some("0/1"));
    let out = annotate_over(GFF3, ACMG_CONFIG, PickFlags::default(), Some("P"), &vcf);

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
