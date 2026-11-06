//! End-to-end regressions for genes with more than one mode of inheritance, and for
//! reading the inheritance ClinGen curates.
//!
//! The phenotypes of a ClinGen build end in a tag the converter writes from the
//! curation's own columns, `(ClinGen Definitive/AD, MONDO:0007254)`. Until the tag
//! was read, the classifier found inheritance only where a disease's *name* said
//! "dominant" or "recessive", so a gene ClinGen curates as autosomal recessive and
//! calls "phenylketonuria" was read as having no inheritance at all, and PM2 held
//! it to the dominant bar. These drive the real `run_annotate` over a GFF3, a FASTA
//! and databases built by `run_sa_build` from files laid out the way
//! `clingen_gdv_to_oga.py` lays them out, so they fail if any link from the
//! curated field to the written classification regresses.
//!
//! The gene has an autosomal dominant disease and an autosomal recessive one, and
//! the variant sits at an allele frequency that is above the dominant PM2 bar
//! (4e-5) and below the recessive one (7e-5). The call for the gene as a whole
//! takes the strict bar; each mode of inheritance is classified on its own beside
//! it.

use fastvep_cli::pipeline::{run_annotate, run_sa_build, AnnotateConfig, PickFlags};
use serde_json::Value;
use std::fs;
use std::path::Path;

/// One gene on the forward strand with one three-exon transcript, and a
/// missense change at 1:165 (`GCC` to `GAC`, Ala to Asp).
const GFF3: &str = "##gff-version 3\n\
1\ttest\tgene\t1\t380\t.\t+\t.\tID=gene:ENSG_DUAL;Name=DUALTEST;gene_name=DUALTEST;biotype=protein_coding\n\
1\ttest\tmRNA\t1\t380\t.\t+\t.\tID=transcript:ENST_DUAL;Parent=gene:ENSG_DUAL;biotype=protein_coding\n\
1\ttest\texon\t1\t60\t.\t+\t.\tID=exon:E1;Parent=transcript:ENST_DUAL;rank=1\n\
1\ttest\texon\t161\t220\t.\t+\t.\tID=exon:E2;Parent=transcript:ENST_DUAL;rank=2\n\
1\ttest\texon\t321\t380\t.\t+\t.\tID=exon:E3;Parent=transcript:ENST_DUAL;rank=3\n\
1\ttest\tCDS\t1\t60\t.\t+\t0\tID=CDS:ENSP_DUAL;Parent=transcript:ENST_DUAL\n\
1\ttest\tCDS\t161\t220\t.\t+\t0\tID=CDS:ENSP_DUAL;Parent=transcript:ENST_DUAL\n\
1\ttest\tCDS\t321\t380\t.\t+\t0\tID=CDS:ENSP_DUAL;Parent=transcript:ENST_DUAL\n";

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

const AD_DISEASE: &str = "a dominant cancer predisposition (ClinGen Definitive/AD, MONDO:0000001)";
const AR_DISEASE: &str = "a recessive syndrome (ClinGen Definitive/AR, MONDO:0000002)";

/// The file `clingen_gdv_to_oga.py` writes: a comment header, then one row per
/// gene with 13 tab-separated columns, the gene symbol in column 8 and the
/// phenotypes joined with `; ` in column 12.
fn genemap(gene: &str, phenotypes: &[&str]) -> String {
    let mut cols = vec![""; 13];
    cols[5] = "0";
    cols[8] = gene;
    let joined = phenotypes.join("; ");
    cols[12] = &joined;
    format!(
        "# Synthetic genemap2 from ClinGen Gene-Disease Validity\n\
         # col5=mim_proxy (always 0), col8=gene_symbol, col12=phenotypes\n{}\n",
        cols.join("\t")
    )
}

/// A gnomAD record for the variant at `AF`, with `AC` of its alleles among 1.4 M.
fn gnomad_vcf(af: f64, ac: u64) -> String {
    format!(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
         1\t165\t.\tC\tA\t.\tPASS\tAF={af};AC={ac};AN=1400000;nhomalt=0\n"
    )
}

const INPUT_VCF: &str = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n1\t165\t.\tC\tA\t.\tPASS\t.\n";

/// Everything a run reads, on disk: the gene models, the reference, the input and
/// a `sa` directory holding the `.oga` and the gnomAD database built by
/// `run_sa_build` from `phenotypes` and from gnomAD at `af`.
fn world(phenotypes: &[&str], af: f64, ac: u64) -> tempfile::TempDir {
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();
    let sa = dir.join("sa");
    fs::create_dir(&sa).unwrap();

    fs::write(dir.join("genemap.txt"), genemap("DUALTEST", phenotypes)).unwrap();
    run_sa_build(
        "omim",
        dir.join("genemap.txt").to_str().unwrap(),
        sa.join("omim").to_str().unwrap(),
        "GRCh38",
        None,
        &[],
        false,
    )
    .unwrap();
    fs::write(dir.join("gnomad.vcf"), gnomad_vcf(af, ac)).unwrap();
    run_sa_build(
        "gnomad",
        dir.join("gnomad.vcf").to_str().unwrap(),
        sa.join("gnomad").to_str().unwrap(),
        "GRCh38",
        None,
        &[],
        false,
    )
    .unwrap();

    fs::write(dir.join("genes.gff3"), GFF3).unwrap();
    fs::write(dir.join("genome.fa"), format!(">1\n{}\n", reference())).unwrap();
    fs::write(dir.join("in.vcf"), INPUT_VCF).unwrap();
    tmp
}

/// Annotate the variant over a gene with `phenotypes` and gnomAD at `af`, with
/// `acmg` on or off, and return what the CLI wrote.
fn annotate(phenotypes: &[&str], af: f64, ac: u64, format: &str, acmg: bool) -> String {
    let tmp = world(phenotypes, af, ac);
    let dir = tmp.path();
    let out = dir.join(format!("out.{format}"));

    run_annotate(AnnotateConfig {
        input: dir.join("in.vcf").to_string_lossy().into(),
        output: out.to_string_lossy().into(),
        gff3: vec![dir.join("genes.gff3").to_string_lossy().into()],
        fasta: Some(dir.join("genome.fa").to_string_lossy().into()),
        output_format: format.into(),
        pick: PickFlags::default(),
        hgvs: true,
        distance: 5000,
        cache_dir: None,
        transcript_cache: None,
        sa_dir: Some(dir.join("sa").to_string_lossy().into()),
        sa_only: false,
        acmg,
        acmg_config: None,
        pick_order: None,
        functional_evidence: None,
        proband: None,
        mother: None,
        father: None,
        gene_list: None,
        explicit_alleles: false,
        qc_rules: None,
        show_progress: false,
        compat: Default::default(),
        columns: Default::default(),
        exclude: Default::default(),
    })
    .expect("annotation should succeed");
    fs::read_to_string(Path::new(&out)).unwrap()
}

/// The same world through the library's own per-variant loop, which is the one
/// the web server runs.
fn annotate_with_the_library(phenotypes: &[&str], af: f64, ac: u64) -> String {
    let tmp = world(phenotypes, af, ac);
    let dir = tmp.path();
    let mut ctx = fastvep_annotate::AnnotationContext::new(
        Some(dir.join("genes.gff3").to_str().unwrap()),
        Some(dir.join("genome.fa").to_str().unwrap()),
        Some(dir.join("sa").to_str().unwrap()),
        5000,
    )
    .expect("the world should load");
    ctx.hgvs = true;
    ctx.acmg_config = Some(fastvep_classification::AcmgConfig::default());
    let vcf = fs::read_to_string(dir.join("in.vcf")).unwrap();
    let variants = ctx
        .annotate_vcf_text(&vcf, None)
        .expect("annotation should succeed");
    serde_json::to_string(&variants).unwrap()
}

/// The `acmg` object of the one transcript row of the one variant.
fn acmg_of(json_text: &str) -> Value {
    let variants: Vec<Value> = serde_json::from_str(json_text).unwrap();
    let rows = variants[0]["transcript_consequences"].as_array().unwrap();
    assert_eq!(rows.len(), 1, "one transcript expected: {}", variants[0]);
    rows[0]["acmg"].clone()
}

fn met(result: &Value, prefix: &str) -> bool {
    result["criteria"]
        .as_array()
        .unwrap()
        .iter()
        .any(|c| c["met"] == true && c["code"].as_str().unwrap().starts_with(prefix))
}

fn criterion<'a>(result: &'a Value, prefix: &str) -> &'a Value {
    result["criteria"]
        .as_array()
        .unwrap()
        .iter()
        .find(|c| c["code"].as_str().unwrap().starts_with(prefix))
        .unwrap_or_else(|| panic!("no {prefix} in {result}"))
}

#[test]
fn a_gene_with_a_dominant_and_a_recessive_disease_is_classified_against_each() {
    let acmg = acmg_of(&annotate(&[AD_DISEASE, AR_DISEASE], 5e-5, 70, "json", true));

    // The call for the gene as a whole reads it as both, which is the strict bar.
    assert!(!met(&acmg, "PM2"), "{acmg}");
    assert_eq!(
        criterion(&acmg, "PM2")["details"]["inheritance_basis"],
        "AD_or_unknown"
    );

    let per = acmg["per_disorder"].as_array().expect("per_disorder");
    let codes: Vec<&str> = per
        .iter()
        .map(|d| d["inheritance"].as_str().unwrap())
        .collect();
    assert_eq!(codes, ["AD", "AR"]);

    let (ad, ar) = (&per[0], &per[1]);
    assert!(!met(ad, "PM2"), "5e-5 is above the dominant bar: {ad}");
    assert!(met(ar, "PM2"), "5e-5 is below the recessive bar: {ar}");
    assert_eq!(criterion(ar, "PM2")["details"]["inheritance_basis"], "AR");

    // Each entry names the diseases it was made against, as ClinGen curates them.
    assert_eq!(ad["disorders"][0]["mondo"], "MONDO:0000001");
    assert_eq!(ad["disorders"][0]["validity"], "Definitive");
    assert_eq!(
        ad["disorders"][0]["label"],
        "a dominant cancer predisposition"
    );
    assert_eq!(ar["disorders"][0]["mondo"], "MONDO:0000002");
    // And reads like a classification: the call, its shorthand and its criteria.
    for entry in per {
        assert!(entry["classification"].is_string());
        assert!(entry["shorthand"].is_string());
        assert!(entry["criteria"].is_array());
    }
}

#[test]
fn a_gene_curated_as_recessive_takes_the_recessive_bar_whatever_its_disease_is_called() {
    // The label says nothing about inheritance. Before the tag was read, this
    // gene was read as having none and PM2 was held to the dominant bar.
    let acmg = acmg_of(&annotate(
        &["phenylketonuria (ClinGen Definitive/AR, MONDO:0009861)"],
        5e-5,
        70,
        "json",
        true,
    ));
    assert!(met(&acmg, "PM2"), "{acmg}");
    assert_eq!(
        criterion(&acmg, "PM2")["details"]["inheritance_basis"],
        "AR"
    );
    assert!(
        acmg.get("per_disorder").is_none(),
        "a gene with one mode of inheritance has nothing to separate: {acmg}"
    );
}

#[test]
fn a_gene_curated_as_dominant_takes_the_strict_bar() {
    let acmg = acmg_of(&annotate(
        &["Kabuki syndrome 1 (ClinGen Definitive/AD, MONDO:0008624)"],
        5e-5,
        70,
        "json",
        true,
    ));
    assert!(!met(&acmg, "PM2"), "{acmg}");
    assert_eq!(criterion(&acmg, "PM2")["details"]["is_dominant"], true);
}

/// The CSQ of a VCF run: the names of its fields, from the header, and the
/// values of the one record's one transcript entry.
fn csq(vcf: &str) -> (Vec<String>, Vec<String>) {
    let header = vcf
        .lines()
        .find(|l| l.starts_with("##INFO=<ID=CSQ,"))
        .expect("a CSQ header");
    let format = header
        .split("Format: ")
        .nth(1)
        .unwrap()
        .trim_end_matches("\">")
        .to_string();
    let names: Vec<String> = format.split('|').map(str::to_string).collect();
    let record = vcf.lines().find(|l| !l.starts_with('#')).unwrap();
    let info = record.split('\t').nth(7).unwrap();
    let csq = info
        .split(';')
        .find_map(|kv| kv.strip_prefix("CSQ="))
        .expect("a CSQ value");
    let entries: Vec<&str> = csq.split(',').collect();
    assert_eq!(entries.len(), 1, "one transcript expected: {csq}");
    (names, entries[0].split('|').map(str::to_string).collect())
}

#[test]
fn the_per_disorder_calls_are_in_the_csq_of_an_acmg_run() {
    let (names, values) = csq(&annotate(&[AD_DISEASE, AR_DISEASE], 5e-5, 70, "vcf", true));
    assert_eq!(names.last().map(String::as_str), Some("ACMG_DISORDERS"));
    assert_eq!(names.len(), values.len());
    let cell = values.last().unwrap();
    let entries: Vec<&str> = cell.split('&').collect();
    // `INHERITANCE:CLASS:CRITERIA`, the criteria joined with `+`: PM2 is met in the
    // recessive context only.
    assert_eq!(entries.len(), 2, "{cell}");
    assert!(entries[0].starts_with("AD:"), "{cell}");
    assert!(entries[1].starts_with("AR:"), "{cell}");
    assert!(!entries[0].contains("PM2"), "{cell}");
    assert!(entries[1].contains("PM2_Supporting"), "{cell}");
    // It sits after the columns every reader of the CSQ already parses.
    let acmg_criteria = names.iter().position(|n| n == "ACMG_CRITERIA").unwrap();
    assert_eq!(names[acmg_criteria + 1], "ACMG_DISORDERS");
}

#[test]
fn the_column_is_empty_for_a_gene_with_one_mode_of_inheritance() {
    let (names, values) = csq(&annotate(
        &["phenylketonuria (ClinGen Definitive/AR, MONDO:0009861)"],
        5e-5,
        70,
        "vcf",
        true,
    ));
    assert_eq!(names.last().map(String::as_str), Some("ACMG_DISORDERS"));
    assert_eq!(values.last().map(String::as_str), Some(""));
}

#[test]
fn a_run_without_acmg_has_the_columns_it_always_had() {
    let (names, _) = csq(&annotate(&[AD_DISEASE, AR_DISEASE], 5e-5, 70, "vcf", false));
    assert!(!names.iter().any(|n| n == "ACMG_DISORDERS"), "{names:?}");
    assert_eq!(names.len(), 49, "the 49 fields every reader already parses");
}

/// The CLI and the library each carry their own per-allele loop, and nothing
/// structurally keeps them in step: the same gene and variant have to come out
/// with the same call for the gene and the same call for each of its modes.
#[test]
fn the_library_loop_classifies_each_disorder_as_the_cli_does() {
    let phenotypes = [AD_DISEASE, AR_DISEASE];
    let cli = acmg_of(&annotate(&phenotypes, 5e-5, 70, "json", true));
    let library = acmg_of(&annotate_with_the_library(&phenotypes, 5e-5, 70));

    assert_eq!(cli["shorthand"], library["shorthand"]);
    let per = |v: &Value| -> Vec<(String, String, bool)> {
        v["per_disorder"]
            .as_array()
            .expect("per_disorder")
            .iter()
            .map(|d| {
                (
                    d["inheritance"].as_str().unwrap().to_string(),
                    d["shorthand"].as_str().unwrap().to_string(),
                    met(d, "PM2"),
                )
            })
            .collect()
    };
    assert_eq!(per(&cli), per(&library));
    assert_eq!(
        per(&library),
        [
            ("AD".to_string(), "VUS".to_string(), false),
            ("AR".to_string(), "VUS".to_string(), true)
        ]
    );
}
