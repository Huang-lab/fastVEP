//! `sa-build --extra-info`: INFO keys stored in addition to a source's own columns.
//!
//! The columns take their type from the input's `##INFO` header, so these tests hold that a
//! v1 `.osa` and a v2 `.osa2` built with the same entries answer with the same values, that a
//! String key reads back as its text and not as the index the `.osa2` stores, and that a build
//! which cannot store what was asked for stops instead of writing a database without it.

use fastvep_cache::annotation::{AnnotationProvider, AnnotationValue};
use fastvep_cli::pipeline::{run_sa_build, run_sa_build_v2, with_extra_info};
use fastvep_sa::reader::SaReader;
use fastvep_sa::reader_v2::Osa2Reader;
use std::fs;
use std::path::Path;

const GNOMAD_VCF: &str = "\
##fileformat=VCFv4.2
##INFO=<ID=AF,Number=A,Type=Float,Description=\"af\">
##INFO=<ID=AN,Number=1,Type=Integer,Description=\"an\">
##INFO=<ID=AC,Number=A,Type=Integer,Description=\"ac\">
##INFO=<ID=AC_nfe,Number=A,Type=Integer,Description=\"Alternate allele count, nfe\">
##INFO=<ID=nhomalt_nfe,Number=A,Type=Integer,Description=\"Homozygote count, nfe\">
##INFO=<ID=fafmax_faf99_max,Number=A,Type=Float,Description=\"Max FAF, 99% CI\">
##INFO=<ID=fafmax_faf95_max_gen_anc,Number=A,Type=String,Description=\"Group of the max FAF\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\t.\tA\tG\t.\tPASS\tAF=0.01;AN=1000;AC=10;AC_nfe=6;nhomalt_nfe=1;fafmax_faf99_max=0.0003;fafmax_faf95_max_gen_anc=nfe
chr1\t200\t.\tGATTACA\tG\t.\tPASS\tAF=0.005;AN=1000;AC=5;AC_nfe=2;fafmax_faf95_max_gen_anc=afr
";

const ENTRIES: &[&str] = &[
    "AC_nfe=nfeAc",
    "nhomalt_nfe=nfeHc",
    "fafmax_faf99_max:af=faf99Max",
    "fafmax_faf95_max_gen_anc=faf95MaxGenAnc",
];

fn entries(list: &[&str]) -> Vec<String> {
    list.iter().map(|s| s.to_string()).collect()
}

fn json_of(value: Option<AnnotationValue>) -> serde_json::Value {
    match value {
        Some(AnnotationValue::Json(j)) => serde_json::from_str(&j).unwrap(),
        other => panic!("expected a JSON annotation, got {other:?}"),
    }
}

fn build_both(dir: &Path, list: &[&str]) -> (SaReader, Osa2Reader) {
    let input = dir.join("gnomad.vcf");
    fs::write(&input, GNOMAD_VCF).unwrap();
    let (v1, v2) = (dir.join("v1"), dir.join("v2"));
    with_extra_info(entries(list), || {
        run_sa_build(
            "gnomad",
            input.to_str().unwrap(),
            v1.to_str().unwrap(),
            "GRCh38",
            None,
            &[],
            false,
        )?;
        run_sa_build_v2(
            "gnomad",
            input.to_str().unwrap(),
            v2.to_str().unwrap(),
            "GRCh38",
            None,
            &[],
            false,
        )
    })
    .unwrap();
    (
        SaReader::open(&v1.with_extension("osa")).unwrap(),
        Osa2Reader::open(&v2.with_extension("osa2")).unwrap(),
    )
}

#[test]
fn v1_and_v2_answer_with_the_same_extra_values() {
    let dir = tempfile::tempdir().unwrap();
    let (v1, v2) = build_both(dir.path(), ENTRIES);

    for (pos, r, a) in [(100, "A", "G"), (200, "GATTACA", "G")] {
        let one = json_of(v1.annotate_position("chr1", pos, r, a).unwrap());
        let two = json_of(v2.annotate_position("chr1", pos, r, a).unwrap());
        for key in ["nfeAc", "nfeHc", "faf95MaxGenAnc"] {
            assert_eq!(
                one.get(key),
                two.get(key),
                "{key} at {pos}: v1={one} v2={two}"
            );
        }
        let faf = |v: &serde_json::Value| v.get("faf99Max").and_then(|x| x.as_f64());
        match (faf(&one), faf(&two)) {
            (Some(x), Some(y)) => assert!((x - y).abs() < 1e-6, "faf99Max at {pos}: {x} vs {y}"),
            (None, None) => {}
            (x, y) => panic!("faf99Max at {pos}: v1={x:?} v2={y:?}"),
        }
    }

    let first = json_of(v2.annotate_position("chr1", 100, "A", "G").unwrap());
    assert_eq!(first["nfeAc"], 6);
    assert_eq!(first["nfeHc"], 1);
    assert_eq!(first["faf95MaxGenAnc"], "nfe");
    let second = json_of(v2.annotate_position("chr1", 200, "GATTACA", "G").unwrap());
    assert_eq!(second["faf95MaxGenAnc"], "afr");
}

#[test]
fn a_build_without_extra_info_has_none_of_the_keys() {
    let dir = tempfile::tempdir().unwrap();
    let (_, v2) = build_both(dir.path(), &[]);
    let json = json_of(v2.annotate_position("chr1", 100, "A", "G").unwrap());
    for key in ["AC_nfe", "nfeAc", "faf99Max", "faf95MaxGenAnc"] {
        assert!(json.get(key).is_none(), "{key} in {json}");
    }
}

#[test]
fn a_build_that_cannot_store_what_was_asked_for_stops() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("in.vcf");
    fs::write(&input, GNOMAD_VCF).unwrap();
    let run = |source: &str, list: &[&str]| {
        with_extra_info(entries(list), || {
            run_sa_build_v2(
                source,
                input.to_str().unwrap(),
                dir.path().join("out").to_str().unwrap(),
                "GRCh38",
                None,
                &[],
                false,
            )
        })
        .unwrap_err()
        .to_string()
    };

    assert!(run("gnomad", &["fafmax_faf99_max"]).contains("needs an encoding"));
    assert!(run("gnomad", &["AC_joint_nfe"]).contains("declares no ##INFO=<ID=AC_joint_nfe>"));
    assert!(run("gnomad", &["AC_nfe=allAc"]).contains("already writes the JSON key 'allAc'"));
    assert!(run("topmed", &["AC_nfe"]).contains("--extra-info is supported for --source gnomad"));
    assert!(!dir.path().join("out.osa2").exists());
}
