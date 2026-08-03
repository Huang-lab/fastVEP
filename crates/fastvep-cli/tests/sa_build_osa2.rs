//! Integration tests for the v2 (`.osa2`) gnomAD builder wired into
//! `fastvep sa-build --format osa2`.
//!
//! Verifies the full pipeline: parse a gnomAD VCF, stream it into a `.osa2`
//! file, then read it back through the same `Osa2Reader` the annotation
//! pipeline uses and confirm the annotations round-trip. Also guards the
//! crash-safety and source-gating contracts.

use fastvep_cache::annotation::{AnnotationProvider, AnnotationValue};
use fastvep_cli::pipeline::{
    run_sa_build, run_sa_build_format, run_sa_build_v2, source_supports_osa2,
};
use fastvep_sa::reader::SaReader;
use fastvep_sa::reader_v2::Osa2Reader;
use std::fs;

const GNOMAD_VCF: &str = "\
##fileformat=VCFv4.2
##INFO=<ID=AF,Number=A,Type=Float,Description=\"af\">
##INFO=<ID=AN,Number=1,Type=Integer,Description=\"an\">
##INFO=<ID=AC,Number=A,Type=Integer,Description=\"ac\">
##INFO=<ID=nhomalt,Number=A,Type=Integer,Description=\"hc\">
##INFO=<ID=AF_nfe,Number=A,Type=Float,Description=\"nfe af\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\t.\tA\tG\t.\tPASS\tAF=0.01;AN=1000;AC=10;nhomalt=1;AF_nfe=0.02
chr1\t200\t.\tGATTACA\tG\t.\tPASS\tAF=0.005;AN=1000;AC=5;nhomalt=0
chr2\t150\t.\tG\tA\t.\tPASS\tAF=0.3;AN=1000;AC=300;nhomalt=40
";

fn af_of(json: &str) -> Option<f64> {
    // Values look like {"allAf":1.000000e-2,...}
    let start = json.find("\"allAf\":")? + "\"allAf\":".len();
    let rest = &json[start..];
    let end = rest
        .find(',')
        .unwrap_or_else(|| rest.find('}').unwrap_or(rest.len()));
    rest[..end].parse().ok()
}

fn json_at(reader: &Osa2Reader, chrom: &str, pos: u64, r: &str, a: &str) -> Option<String> {
    match reader.annotate_position(chrom, pos, r, a).unwrap()? {
        AnnotationValue::Json(j) | AnnotationValue::Positional(j) => Some(j),
        other => panic!("unexpected value: {:?}", other),
    }
}

#[test]
fn gnomad_osa2_build_round_trips() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("gnomad.vcf");
    let output = dir.path().join("gnomad"); // builder appends .osa2
    fs::write(&input, GNOMAD_VCF).unwrap();

    run_sa_build_v2(
        "gnomad",
        input.to_str().unwrap(),
        output.to_str().unwrap(),
        "GRCh38",
        false,
    )
    .unwrap();

    let osa2 = output.with_extension("osa2");
    assert!(osa2.exists(), ".osa2 file should be created");

    let reader = Osa2Reader::open(&osa2).unwrap();
    assert_eq!(reader.json_key(), "gnomad");

    // SNV: AF encodes and decodes ~0.01.
    let j = json_at(&reader, "chr1", 100, "A", "G").expect("chr1:100 A>G annotated");
    let af = af_of(&j).expect("allAf present");
    assert!((af - 0.01).abs() < 1e-5, "AF round-trip off: {j}");
    assert!(j.contains("\"allAc\":10"), "AC missing: {j}");
    assert!(j.contains("\"nfeAf\":"), "per-pop nfe AF missing: {j}");

    // Indel (long variant, ref+alt > 4 bases) must round-trip its own values.
    let j = json_at(&reader, "chr1", 200, "GATTACA", "G").expect("chr1:200 indel annotated");
    let af = af_of(&j).expect("indel allAf present");
    assert!((af - 0.005).abs() < 1e-5, "indel AF round-trip off: {j}");
    assert!(j.contains("\"allAc\":5"), "indel AC missing: {j}");

    // Different chromosome.
    let j = json_at(&reader, "chr2", 150, "G", "A").expect("chr2:150 annotated");
    assert!(
        af_of(&j).map(|af| (af - 0.3).abs() < 1e-5).unwrap_or(false),
        "chr2 AF: {j}"
    );

    // A definite miss.
    assert!(json_at(&reader, "chr1", 999, "A", "G").is_none());
    // Wrong allele misses.
    assert!(json_at(&reader, "chr1", 100, "A", "T").is_none());
}

const ONEKG_VCF: &str = "\
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t10001\t.\tA\tG\t.\t.\tAF=0.15;AFR_AF=0.20;EUR_AF=0.10
chr1\t20000\t.\tC\tT\t.\t.\tAF=0.4;EAS_AF=0.5;SAS_AF=0.33
chr2\t500\t.\tGATTACA\tG\t.\t.\tAF=0.02;AMR_AF=0.03
";

const TOPMED_VCF: &str = "\
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\t.\tA\tG\t.\t.\tAF=0.05;AC=500;AN=10000
chr1\t200\t.\tGATTACA\tG\t.\t.\tAF=0.001;AC=10;AN=10000
chr2\t300\t.\tT\tC\t.\t.\tAF=0.9;AC=9000;AN=10000
";

/// Parse the two JSON payloads and assert every numeric field matches within a
/// small tolerance — the definitive "v2 output equals v1 output" check.
fn assert_json_equivalent(v1: &str, v2: &str) {
    let a: serde_json::Value = serde_json::from_str(v1).expect("v1 json");
    let b: serde_json::Value = serde_json::from_str(v2).expect("v2 json");
    let a = a.as_object().expect("v1 object");
    let b = b.as_object().expect("v2 object");
    assert_eq!(
        a.len(),
        b.len(),
        "different key counts:\n v1={v1}\n v2={v2}"
    );
    for (k, av) in a {
        let bv = b
            .get(k)
            .unwrap_or_else(|| panic!("v2 missing key {k}: {v2}"));
        match (av.as_f64(), bv.as_f64()) {
            (Some(x), Some(y)) => assert!(
                (x - y).abs() <= 1e-6 * x.abs().max(1.0),
                "field {k} differs: v1={x} v2={y}"
            ),
            _ => assert_eq!(av, bv, "field {k} differs (non-numeric)"),
        }
    }
}

/// Build the same VCF as v1 `.osa` and v2 `.osa2`, then confirm both readers
/// return equivalent annotations at every site — including an indel.
fn assert_v1_v2_parity(source: &str, vcf: &str, sites: &[(&str, u64, &str, &str)]) {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("in.vcf");
    fs::write(&input, vcf).unwrap();

    let v1_base = dir.path().join("v1");
    let v2_base = dir.path().join("v2");
    run_sa_build(
        source,
        input.to_str().unwrap(),
        v1_base.to_str().unwrap(),
        "GRCh38",
        None,
        &[],
        false,
    )
    .unwrap();
    run_sa_build_v2(
        source,
        input.to_str().unwrap(),
        v2_base.to_str().unwrap(),
        "GRCh38",
        false,
    )
    .unwrap();

    let v1 = SaReader::open(&v1_base.with_extension("osa")).unwrap();
    let v2 = Osa2Reader::open(&v2_base.with_extension("osa2")).unwrap();
    assert_eq!(
        v1.json_key(),
        v2.json_key(),
        "json_key must match between formats"
    );

    for &(chrom, pos, r, a) in sites {
        let j1 = match v1.annotate_position(chrom, pos, r, a).unwrap() {
            Some(AnnotationValue::Json(j)) | Some(AnnotationValue::Positional(j)) => j,
            other => panic!("v1 missing {chrom}:{pos} {r}>{a}: {other:?}"),
        };
        let j2 = match v2.annotate_position(chrom, pos, r, a).unwrap() {
            Some(AnnotationValue::Json(j)) | Some(AnnotationValue::Positional(j)) => j,
            other => panic!("v2 missing {chrom}:{pos} {r}>{a}: {other:?}"),
        };
        assert_json_equivalent(&j1, &j2);
    }
}

#[test]
fn onekg_v1_v2_output_parity() {
    assert_v1_v2_parity(
        "onekg",
        ONEKG_VCF,
        &[
            ("chr1", 10001, "A", "G"),
            ("chr1", 20000, "C", "T"),
            ("chr2", 500, "GATTACA", "G"), // indel
        ],
    );
}

#[test]
fn topmed_v1_v2_output_parity() {
    assert_v1_v2_parity(
        "topmed",
        TOPMED_VCF,
        &[
            ("chr1", 100, "A", "G"),
            ("chr1", 200, "GATTACA", "G"), // indel
            ("chr2", 300, "T", "C"),
        ],
    );
}

// AlphaMissense TSV: leading comment/license lines, a header, then data rows.
// Column order matches the real Zenodo release. Includes a multi-base ref/alt
// row to exercise the long-variant (kmer16) path alongside the categorical
// class column.
const ALPHAMISSENSE_TSV: &str = "\
# Copyright 2023 DeepMind Technologies Limited
#
#CHROM\tPOS\tREF\tALT\tgenome\tuniprot_id\ttranscript_id\tprotein_variant\tam_pathogenicity\tam_class
chr1\t100\tA\tG\thg38\tQ1\tENST1\tV2L\t0.2937\tlikely_benign
chr1\t100\tA\tT\thg38\tQ1\tENST1\tV2M\t0.9800\tlikely_pathogenic
chr1\t200\tGATTACA\tG\thg38\tQ1\tENST1\tX9Y\t0.4500\tambiguous
chr2\t150\tC\tG\thg38\tQ2\tENST2\tR3Q\t0.0125\tlikely_benign
";

#[test]
fn alphamissense_v1_v2_output_parity() {
    // Every field at every site must match between v1 .osa and v2 .osa2 —
    // including the categorical amClass and the long (indel-shaped) variant.
    assert_v1_v2_parity(
        "alphamissense",
        ALPHAMISSENSE_TSV,
        &[
            ("chr1", 100, "A", "G"),
            ("chr1", 100, "A", "T"),
            ("chr2", 150, "C", "G"),
            ("chr1", 200, "GATTACA", "G"), // long variant + categorical
        ],
    );
}

#[test]
fn alphamissense_osa2_round_trips_score_and_class() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("am.tsv");
    let output = dir.path().join("am"); // builder appends .osa2
    fs::write(&input, ALPHAMISSENSE_TSV).unwrap();

    run_sa_build_v2(
        "alphamissense",
        input.to_str().unwrap(),
        output.to_str().unwrap(),
        "GRCh38",
        false,
    )
    .unwrap();

    let reader = Osa2Reader::open(&output.with_extension("osa2")).unwrap();
    assert_eq!(reader.json_key(), "alphaMissense");

    // SNV: score + benign class round-trip.
    let j = json_at(&reader, "chr1", 100, "A", "G").expect("chr1:100 A>G");
    assert!(j.contains("\"amClass\":\"likely_benign\""), "class: {j}");
    let score = {
        let key = "\"amPathogenicity\":";
        let start = j.find(key).expect("score key") + key.len();
        let rest = &j[start..];
        let end = rest.find(|c| c == ',' || c == '}').unwrap();
        rest[..end].parse::<f64>().expect("score float")
    };
    assert!((score - 0.2937).abs() < 1e-5, "score round-trip: {j}");

    // Same position, different allele resolves to its own class.
    let j = json_at(&reader, "chr1", 100, "A", "T").expect("chr1:100 A>T");
    assert!(
        j.contains("\"amClass\":\"likely_pathogenic\""),
        "class: {j}"
    );

    // Long variant carries the ambiguous class.
    let j = json_at(&reader, "chr1", 200, "GATTACA", "G").expect("chr1:200 indel");
    assert!(j.contains("\"amClass\":\"ambiguous\""), "class: {j}");

    // A definite miss.
    assert!(json_at(&reader, "chr1", 999, "A", "G").is_none());
}

// dbSNP VCF: RS IDs + CAF (per-ALT allele frequencies), a bare SNV with no
// CAF, a multi-allelic site, and an indel (long variant).
const DBSNP_VCF: &str = "\
##fileformat=VCFv4.0
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\trs101\tA\tG\t.\t.\tRS=101;CAF=0.99,0.01
chr1\t200\trs202\tC\tT\t.\t.\tRS=202
chr1\t300\trs303\tA\tC,T\t.\t.\tRS=303;CAF=0.90,0.08,0.02
chr2\t150\trs404\tGATTACA\tG\t.\t.\tRS=404;CAF=0.995,0.005
";

#[test]
fn dbsnp_v1_v2_output_parity() {
    // Whole-record-blob v2 must return byte-equivalent JSON to v1 at every
    // site, including the no-MAF record, both alleles of the multi-allelic
    // site, and the indel (long variant).
    assert_v1_v2_parity(
        "dbsnp",
        DBSNP_VCF,
        &[
            ("chr1", 100, "A", "G"),
            ("chr1", 200, "C", "T"),
            ("chr1", 300, "A", "C"),
            ("chr1", 300, "A", "T"),
            ("chr2", 150, "GATTACA", "G"), // long variant
        ],
    );
}

#[test]
fn dbsnp_osa2_stores_id_and_maf() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("dbsnp.vcf");
    let output = dir.path().join("dbsnp");
    fs::write(&input, DBSNP_VCF).unwrap();

    run_sa_build_v2(
        "dbsnp",
        input.to_str().unwrap(),
        output.to_str().unwrap(),
        "GRCh38",
        false,
    )
    .unwrap();

    let reader = Osa2Reader::open(&output.with_extension("osa2")).unwrap();
    assert_eq!(reader.json_key(), "dbsnp");

    let j = json_at(&reader, "chr1", 100, "A", "G").expect("chr1:100");
    assert!(j.contains("\"id\":\"rs101\""), "id: {j}");
    assert!(j.contains("\"globalMaf\":"), "maf: {j}");

    // No-CAF record carries the id but no globalMaf.
    let j = json_at(&reader, "chr1", 200, "C", "T").expect("chr1:200");
    assert!(j.contains("\"id\":\"rs202\""), "id: {j}");
    assert!(!j.contains("globalMaf"), "should have no MAF: {j}");

    // Indel round-trips its blob.
    let j = json_at(&reader, "chr2", 150, "GATTACA", "G").expect("chr2:150 indel");
    assert!(j.contains("\"id\":\"rs404\""), "indel id: {j}");

    assert!(json_at(&reader, "chr1", 999, "A", "G").is_none());
}

// COSMIC coding-mutations VCF: COSV id + GENE + CNT, plus a bare record
// (no gene/count) and an indel.
const COSMIC_VCF: &str = "\
##fileformat=VCFv4.1
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\tCOSV1\tA\tG\t.\t.\tGENE=TP53;CNT=15
chr1\t250\tCOSV2\tC\tT\t.\t.\tGENE=BRAF
chr2\t150\tCOSV3\tGATTACA\tG\t.\t.\tGENE=EGFR;CNT=3
";

#[test]
fn cosmic_v1_v2_output_parity() {
    assert_v1_v2_parity(
        "cosmic",
        COSMIC_VCF,
        &[
            ("chr1", 100, "A", "G"),
            ("chr1", 250, "C", "T"),
            ("chr2", 150, "GATTACA", "G"), // long variant
        ],
    );
}

#[test]
fn cosmic_osa2_stores_id_gene_count() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("cosmic.vcf");
    let output = dir.path().join("cosmic");
    fs::write(&input, COSMIC_VCF).unwrap();

    run_sa_build_v2(
        "cosmic",
        input.to_str().unwrap(),
        output.to_str().unwrap(),
        "GRCh38",
        false,
    )
    .unwrap();

    let reader = Osa2Reader::open(&output.with_extension("osa2")).unwrap();
    assert_eq!(reader.json_key(), "cosmic");

    let j = json_at(&reader, "chr1", 100, "A", "G").expect("chr1:100");
    assert!(j.contains("\"id\":\"COSV1\""), "id: {j}");
    assert!(j.contains("\"gene\":\"TP53\""), "gene: {j}");
    assert!(j.contains("\"count\":15"), "count: {j}");

    // No CNT → no count key.
    let j = json_at(&reader, "chr1", 250, "C", "T").expect("chr1:250");
    assert!(j.contains("\"gene\":\"BRAF\""), "gene: {j}");
    assert!(!j.contains("count"), "should have no count: {j}");
}

// ClinVar VCF: significance/phenotype arrays, review status, population AFs
// (a nested-object payload), plus an indel.
const CLINVAR_VCF: &str = "\
##fileformat=VCFv4.1
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\trs1\tA\tG\t.\t.\tCLNSIG=Pathogenic|Likely_pathogenic;CLNREVSTAT=criteria_provided,_multiple_submitters;CLNDN=Breast_cancer|Ovarian_cancer;CLNVC=single_nucleotide_variant;AF_EXAC=0.00054
chr1\t250\trs2\tC\tT\t.\t.\tCLNSIG=Benign;CLNREVSTAT=criteria_provided,_single_submitter;CLNDN=not_provided
chr2\t150\trs3\tGATTACA\tG\t.\t.\tCLNSIG=Uncertain_significance;CLNDN=Cardiomyopathy
";

#[test]
fn clinvar_v1_v2_output_parity() {
    // ClinVar's nested arrays (significance, phenotypes) must survive the
    // whole-record-blob round-trip byte-for-byte, at every site including the
    // indel.
    assert_v1_v2_parity(
        "clinvar",
        CLINVAR_VCF,
        &[
            ("chr1", 100, "A", "G"),
            ("chr1", 250, "C", "T"),
            ("chr2", 150, "GATTACA", "G"), // long variant
        ],
    );
}

#[test]
fn clinvar_osa2_preserves_arrays_and_is_array_flag() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("clinvar.vcf");
    let output = dir.path().join("clinvar");
    fs::write(&input, CLINVAR_VCF).unwrap();

    run_sa_build_v2(
        "clinvar",
        input.to_str().unwrap(),
        output.to_str().unwrap(),
        "GRCh38",
        false,
    )
    .unwrap();

    let reader = Osa2Reader::open(&output.with_extension("osa2")).unwrap();
    assert_eq!(reader.json_key(), "clinvar");
    // is_array metadata is preserved so downstream consumers behave as with v1.
    assert!(
        reader.metadata().is_array,
        "clinvar .osa2 must keep is_array=true"
    );

    let j = json_at(&reader, "chr1", 100, "A", "G").expect("chr1:100");
    // Nested significance array with two entries round-trips.
    assert!(
        j.contains("\"significance\":[\"Pathogenic\",\"Likely_pathogenic\"]"),
        "significance array: {j}"
    );
    assert!(
        j.contains("\"phenotypes\":[\"Breast_cancer\",\"Ovarian_cancer\"]"),
        "phenotypes: {j}"
    );
    assert!(j.contains("\"afExac\":0.00054"), "population AF: {j}");
    // The reconstructed JSON must still parse as a valid object.
    let v: serde_json::Value = serde_json::from_str(&j).expect("valid JSON");
    assert!(v.get("significance").and_then(|s| s.as_array()).is_some());
}

// REVEL CSV: chr,hg19_pos,grch38_pos,ref,alt,aaref,aaalt,REVEL (position from
// column 2, the GRCh38 coordinate).
const REVEL_CSV: &str = "\
chr,hg19_pos,grch38_pos,ref,alt,aaref,aaalt,REVEL
1,35142,100,G,A,T,M,0.027
1,35142,100,G,C,T,S,0.842
2,50000,200,C,T,R,Q,0.5
";

#[test]
fn revel_v1_v2_output_parity() {
    assert_v1_v2_parity(
        "revel",
        REVEL_CSV,
        &[
            ("chr1", 100, "G", "A"),
            ("chr1", 100, "G", "C"),
            ("chr2", 200, "C", "T"),
        ],
    );
}

#[test]
fn revel_osa2_stores_score() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("revel.csv");
    let output = dir.path().join("revel");
    fs::write(&input, REVEL_CSV).unwrap();

    run_sa_build_v2(
        "revel",
        input.to_str().unwrap(),
        output.to_str().unwrap(),
        "GRCh38",
        false,
    )
    .unwrap();

    let reader = Osa2Reader::open(&output.with_extension("osa2")).unwrap();
    assert_eq!(reader.json_key(), "revel");
    let j = json_at(&reader, "chr1", 100, "G", "C").expect("chr1:100 G>C");
    assert!(j.contains("\"score\":0.842"), "score: {j}");
    // Same position, other allele resolves to its own score.
    let j = json_at(&reader, "chr1", 100, "G", "A").expect("chr1:100 G>A");
    assert!(j.contains("\"score\":0.027"), "score: {j}");
}

// PrimateAI TSV: chr, pos, ref, alt, score (tab-separated).
const PRIMATEAI_TSV: &str = "\
chr\tpos\tref\talt\tscore
1\t100\tA\tG\t0.1234
1\t100\tA\tT\t0.9010
2\t200\tC\tG\t0.5500
";

#[test]
fn primateai_v1_v2_output_parity() {
    assert_v1_v2_parity(
        "primateai",
        PRIMATEAI_TSV,
        &[
            ("chr1", 100, "A", "G"),
            ("chr1", 100, "A", "T"),
            ("chr2", 200, "C", "G"),
        ],
    );
}

#[test]
fn primateai_osa2_stores_score() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("primateai.tsv");
    let output = dir.path().join("primateai");
    fs::write(&input, PRIMATEAI_TSV).unwrap();

    run_sa_build_v2(
        "primateai",
        input.to_str().unwrap(),
        output.to_str().unwrap(),
        "GRCh38",
        false,
    )
    .unwrap();

    let reader = Osa2Reader::open(&output.with_extension("osa2")).unwrap();
    assert_eq!(reader.json_key(), "primateAI");
    let j = json_at(&reader, "chr1", 100, "A", "T").expect("chr1:100 A>T");
    assert!(j.contains("\"score\":0.901"), "score: {j}");
}

// dbNSFP TSV: header-driven column detection, composite SIFT/PolyPhen
// prediction strings, plus a record missing SIFT.
const DBNSFP_TSV: &str = "\
#chr\tpos(1-based)\tref\talt\tSIFT_score\tSIFT_pred\tPolyphen2_HDIV_score\tPolyphen2_HDIV_pred
1\t100\tA\tG\t0.012\tD\t0.998\tD
1\t150\tC\tT\t0.512\tT\t0.045\tB
2\t200\tG\tA\t.\t.\t0.900\tP
";

#[test]
fn dbnsfp_v1_v2_output_parity() {
    assert_v1_v2_parity(
        "dbnsfp",
        DBNSFP_TSV,
        &[
            ("chr1", 100, "A", "G"),
            ("chr1", 150, "C", "T"),
            ("chr2", 200, "G", "A"), // SIFT absent, PolyPhen present
        ],
    );
}

#[test]
fn dbnsfp_osa2_stores_predictions() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("dbnsfp.tsv");
    let output = dir.path().join("dbnsfp");
    fs::write(&input, DBNSFP_TSV).unwrap();

    run_sa_build_v2(
        "dbnsfp",
        input.to_str().unwrap(),
        output.to_str().unwrap(),
        "GRCh38",
        false,
    )
    .unwrap();

    let reader = Osa2Reader::open(&output.with_extension("osa2")).unwrap();
    assert_eq!(reader.json_key(), "dbnsfp");
    let j = json_at(&reader, "chr1", 100, "A", "G").expect("chr1:100");
    assert!(j.contains("\"sift\":"), "sift: {j}");
    assert!(j.contains("\"polyphen\":"), "polyphen: {j}");
    // The record with no SIFT keeps its PolyPhen only.
    let j = json_at(&reader, "chr2", 200, "G", "A").expect("chr2:200");
    assert!(!j.contains("\"sift\":"), "should have no sift: {j}");
    assert!(j.contains("\"polyphen\":"), "polyphen: {j}");
}

// Positional per-base score TSV (chrom, 1-based pos, score) — the shape
// PhyloP (auto-detected) / GERP / DANN consume.
const SCORE_TSV: &str = "\
chr1\t100\t2.345
chr1\t101\t-1.5
chr1\t500000\t0.001
chr2\t50\t3.14
";

/// Positional sources return a bare-number `Positional` payload, not a JSON
/// object, so compare the raw payloads (and confirm both formats agree a miss
/// is a miss). Alleles are irrelevant to positional lookup.
fn assert_v1_v2_positional_parity(source: &str, data: &str, sites: &[(&str, u64, bool)]) {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("scores.tsv");
    fs::write(&input, data).unwrap();
    let v1_base = dir.path().join("v1");
    let v2_base = dir.path().join("v2");
    run_sa_build(
        source,
        input.to_str().unwrap(),
        v1_base.to_str().unwrap(),
        "GRCh38",
        None,
        &[],
        false,
    )
    .unwrap();
    run_sa_build_v2(
        source,
        input.to_str().unwrap(),
        v2_base.to_str().unwrap(),
        "GRCh38",
        false,
    )
    .unwrap();

    let v1 = SaReader::open(&v1_base.with_extension("osa")).unwrap();
    let v2 = Osa2Reader::open(&v2_base.with_extension("osa2")).unwrap();
    assert_eq!(v1.json_key(), v2.json_key());
    assert!(
        v1.metadata().is_positional && v2.metadata().is_positional,
        "both must be positional"
    );

    for &(chrom, pos, expect_hit) in sites {
        let a = v1.annotate_position(chrom, pos, "A", "G").unwrap();
        let b = v2.annotate_position(chrom, pos, "A", "G").unwrap();
        match (a, b) {
            (Some(AnnotationValue::Positional(x)), Some(AnnotationValue::Positional(y))) => {
                assert!(expect_hit, "{chrom}:{pos} unexpectedly hit");
                assert_eq!(x, y, "positional payload differs at {chrom}:{pos}");
            }
            (None, None) => assert!(!expect_hit, "{chrom}:{pos} unexpectedly missed"),
            (x, y) => panic!("v1/v2 disagree at {chrom}:{pos}: {x:?} vs {y:?}"),
        }
    }
}

#[test]
fn phylop_v1_v2_positional_parity() {
    assert_v1_v2_positional_parity(
        "phylop",
        SCORE_TSV,
        &[
            ("chr1", 100, true),
            ("chr1", 101, true),
            ("chr1", 500000, true),
            ("chr2", 50, true),
            ("chr1", 999, false), // miss
        ],
    );
}

#[test]
fn gerp_v1_v2_positional_parity() {
    // GERP uses the plain score-TSV path (no wig auto-detect).
    assert_v1_v2_positional_parity(
        "gerp",
        SCORE_TSV,
        &[
            ("chr1", 100, true),
            ("chr2", 50, true),
            ("chr2", 999, false),
        ],
    );
}

#[test]
fn positional_lookup_ignores_alleles() {
    // A positional source matches by coordinate alone: any query allele at a
    // covered position returns the same score.
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("scores.tsv");
    fs::write(&input, SCORE_TSV).unwrap();
    let out = dir.path().join("phylop");
    run_sa_build_v2(
        "phylop",
        input.to_str().unwrap(),
        out.to_str().unwrap(),
        "GRCh38",
        false,
    )
    .unwrap();

    let reader = Osa2Reader::open(&out.with_extension("osa2")).unwrap();
    assert_eq!(reader.json_key(), "phylop");
    let a = json_at(&reader, "chr1", 100, "A", "G").expect("A>G");
    let b = json_at(&reader, "chr1", 100, "C", "T").expect("C>T");
    let c = json_at(&reader, "chr1", 100, "GATTACA", "G").expect("indel-shaped query");
    assert_eq!(a, "2.345");
    assert_eq!(
        a, b,
        "different alleles must return the same positional score"
    );
    assert_eq!(a, c, "even an indel-shaped query resolves by position only");
}

#[test]
fn osa2_rejects_unsupported_source() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("in.vcf");
    fs::write(&input, GNOMAD_VCF).unwrap();
    let output = dir.path().join("out");

    // mitomap has no v2 encoder — forcing --format osa2 must error, not
    // silently build v1.
    let err = run_sa_build_v2(
        "mitomap",
        input.to_str().unwrap(),
        output.to_str().unwrap(),
        "GRCh38",
        false,
    )
    .unwrap_err();
    assert!(
        err.to_string()
            .contains("currently supported for --source gnomad"),
        "expected source-gating error, got: {err}"
    );
    assert!(
        !output.with_extension("osa2").exists(),
        "no file on rejected build"
    );
}

#[test]
fn source_supports_osa2_membership() {
    // Numeric, blob, and positional sources are v2-capable; gene-level,
    // MitoMap, SpliceAI, and custom sources are not.
    for s in [
        "gnomad",
        "onekg",
        "1000g",
        "topmed",
        "alphamissense",
        "dbsnp",
        "cosmic",
        "clinvar",
        "revel",
        "primateai",
        "dbnsfp",
        "phylop",
        "gerp",
        "dann",
    ] {
        assert!(source_supports_osa2(s), "{s} should support osa2");
    }
    for s in ["mitomap", "spliceai", "omim", "gnomad_genes", "custom_vcf"] {
        assert!(!source_supports_osa2(s), "{s} should NOT support osa2");
    }
}

#[test]
fn format_auto_picks_v2_for_supported_and_v1_for_others() {
    let dir = tempfile::tempdir().unwrap();

    // gnomAD supports v2 → auto must produce .osa2 (and not .osa).
    let gnomad_in = dir.path().join("gnomad.vcf");
    fs::write(&gnomad_in, GNOMAD_VCF).unwrap();
    let gnomad_out = dir.path().join("gnomad");
    run_sa_build_format(
        "auto",
        "gnomad",
        gnomad_in.to_str().unwrap(),
        gnomad_out.to_str().unwrap(),
        "GRCh38",
        None,
        &[],
        false,
    )
    .unwrap();
    assert!(
        gnomad_out.with_extension("osa2").exists(),
        "auto+gnomad should build .osa2"
    );
    assert!(
        !gnomad_out.with_extension("osa").exists(),
        "auto+gnomad should not build .osa"
    );

    // MitoMap has no v2 encoder → auto must fall back to v1 .osa.
    let mito_in = dir.path().join("mitomap.vcf");
    fs::write(
        &mito_in,
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
         chrMT\t100\t.\tA\tG\t.\t.\tDisease=Test;DiseaseStatus=Confirmed\n",
    )
    .unwrap();
    let mito_out = dir.path().join("mitomap");
    run_sa_build_format(
        "auto",
        "mitomap",
        mito_in.to_str().unwrap(),
        mito_out.to_str().unwrap(),
        "GRCh38",
        None,
        &[],
        false,
    )
    .unwrap();
    assert!(
        mito_out.with_extension("osa").exists(),
        "auto+mitomap should build .osa"
    );
    assert!(
        !mito_out.with_extension("osa2").exists(),
        "auto+mitomap should not build .osa2"
    );
}
