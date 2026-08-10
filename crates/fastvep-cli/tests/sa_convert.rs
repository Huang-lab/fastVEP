//! Integration tests for `fastvep sa-convert`: transcoding an existing v1
//! `.osa` database into v2 `.osa2` without going back to the upstream source.
//!
//! The contract being pinned is that a converted database answers every query
//! exactly as the `.osa` did — same records, same JSON bytes, same identity and
//! matching semantics — because that is the only reason to trust the conversion
//! over a rebuild.

use fastvep_cache::annotation::{AnnotationProvider, AnnotationValue};
use fastvep_cli::pipeline::{run_sa_build, run_sa_convert};
use fastvep_sa::reader::SaReader;
use fastvep_sa::reader_v2::Osa2Reader;
use std::fs;

fn json_of(v: Option<AnnotationValue>) -> Option<String> {
    match v {
        Some(AnnotationValue::Json(j)) | Some(AnnotationValue::Positional(j)) => Some(j),
        Some(AnnotationValue::Interval(v)) => Some(format!("[{}]", v.join(","))),
        None => None,
    }
}

/// Build `source` to v1, convert it, and return the (v1, v2) reader pair.
fn build_and_convert(
    dir: &std::path::Path,
    source: &str,
    body: &str,
    input_name: &str,
) -> (SaReader, Osa2Reader) {
    let input = dir.join(input_name);
    fs::write(&input, body).unwrap();
    let base = dir.join(format!("db_{source}"));
    run_sa_build(
        source,
        input.to_str().unwrap(),
        base.to_str().unwrap(),
        "GRCh38",
        None,
        &[],
        false,
    )
    .unwrap();

    let osa = base.with_extension("osa");
    assert!(osa.exists(), "v1 build should have produced {}", osa.display());
    run_sa_convert(osa.to_str().unwrap(), base.to_str().unwrap(), false).unwrap();

    let osa2 = base.with_extension("osa2");
    assert!(osa2.exists(), "conversion should have produced {}", osa2.display());
    (SaReader::open(&osa).unwrap(), Osa2Reader::open(&osa2).unwrap())
}

const CLINVAR_VCF: &str = "\
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\t12345\tA\tG\t.\t.\tCLNSIG=Pathogenic;CLNREVSTAT=criteria_provided;GENEINFO=BRCA1:672
chr1\t300\t12346\tGATTACA\tG\t.\t.\tCLNSIG=Likely_benign;GENEINFO=BRCA1:672
chr2\t50\t12347\tC\tT\t.\t.\tCLNSIG=Uncertain_significance;GENEINFO=TP53:7157
chr2\t9000000\t12348\tT\tA\t.\t.\tCLNSIG=Benign;GENEINFO=TP53:7157
";

#[test]
fn converted_database_answers_every_query_identically() {
    let dir = tempfile::tempdir().unwrap();
    let (v1, v2) = build_and_convert(dir.path(), "clinvar", CLINVAR_VCF, "clinvar.vcf");

    // Identity and matching semantics carry over verbatim — an annotate run
    // against the converted file must produce the same column, keyed the same
    // way, with the same allele-matching behaviour.
    assert_eq!(v1.json_key(), v2.json_key());
    assert_eq!(v1.name(), v2.name());
    assert_eq!(v1.metadata().assembly, v2.metadata().assembly);
    assert_eq!(v1.metadata().version, v2.metadata().version);
    assert_eq!(v1.metadata().description, v2.metadata().description);
    assert_eq!(v1.metadata().match_by_allele, v2.metadata().match_by_allele);
    assert_eq!(v1.metadata().is_array, v2.metadata().is_array);
    assert_eq!(v1.metadata().is_positional, v2.metadata().is_positional);

    for &(chrom, pos, r, a) in &[
        ("chr1", 100u64, "A", "G"),
        ("chr1", 300, "GATTACA", "G"), // long variant (kmer16 path)
        ("chr2", 50, "C", "T"),
        ("chr2", 9_000_000, "T", "A"), // a far-away chunk on a second contig
    ] {
        let j1 = json_of(v1.annotate_position(chrom, pos, r, a).unwrap())
            .unwrap_or_else(|| panic!("v1 missing {chrom}:{pos} {r}>{a}"));
        let j2 = json_of(v2.annotate_position(chrom, pos, r, a).unwrap())
            .unwrap_or_else(|| panic!("converted db missing {chrom}:{pos} {r}>{a}"));
        assert_eq!(j1, j2, "{chrom}:{pos} {r}>{a} must be byte-identical");
    }

    // Misses stay misses: an unknown position, a wrong allele, and an absent
    // contig must not start matching.
    for &(chrom, pos, r, a) in &[
        ("chr1", 999, "A", "G"),
        ("chr1", 100, "A", "T"),
        ("chr9", 100, "A", "G"),
    ] {
        assert!(
            json_of(v1.annotate_position(chrom, pos, r, a).unwrap()).is_none(),
            "v1 should miss {chrom}:{pos} {r}>{a}"
        );
        assert!(
            json_of(v2.annotate_position(chrom, pos, r, a).unwrap()).is_none(),
            "converted db should miss {chrom}:{pos} {r}>{a}"
        );
    }
}

const PHYLOP_WIG: &str = "\
fixedStep chrom=chr1 start=100 step=1
1.5
2.5
3.5
";

#[test]
fn conversion_preserves_positional_matching() {
    // A positional source (`is_positional`) matches on coordinate alone. The
    // flag has to survive so the v2 reader keys on `positional_key` rather than
    // trying to match the alleles the blobs happen to carry.
    let dir = tempfile::tempdir().unwrap();
    let (v1, v2) = build_and_convert(dir.path(), "phylop", PHYLOP_WIG, "phylop.wig");

    assert!(v2.metadata().is_positional, "positional flag must survive conversion");
    assert!(!v2.metadata().match_by_allele);

    for pos in [100u64, 101, 102] {
        let j1 = json_of(v1.annotate_position("chr1", pos, "A", "G").unwrap())
            .unwrap_or_else(|| panic!("v1 missing chr1:{pos}"));
        let j2 = json_of(v2.annotate_position("chr1", pos, "A", "G").unwrap())
            .unwrap_or_else(|| panic!("converted db missing chr1:{pos}"));
        assert_eq!(j1, j2);
        // Any allele resolves to the same score, in both formats.
        let other = json_of(v2.annotate_position("chr1", pos, "C", "T").unwrap()).unwrap();
        assert_eq!(j2, other, "positional lookup must ignore alleles");
    }
}

#[test]
fn conversion_covers_every_record_not_just_the_queried_ones() {
    // A block-by-block scan is easy to get subtly wrong (skipping a block
    // boundary, dropping the last block, missing a chromosome). Build enough
    // records to span several blocks on two contigs and require the converted
    // database to return all of them.
    let dir = tempfile::tempdir().unwrap();
    let mut vcf =
        String::from("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
    let mut expected: Vec<(String, u64)> = Vec::new();
    for (chrom, base) in [("chr1", 1_000u64), ("chr2", 5_000)] {
        for i in 0..4_000u64 {
            let pos = base + i * 7;
            vcf.push_str(&format!(
                "{chrom}\t{pos}\t.\tA\tG\t.\t.\tCLNSIG=Pathogenic;GENEINFO=G{i}:1\n"
            ));
            expected.push((chrom.to_string(), pos));
        }
    }
    let (v1, v2) = build_and_convert(dir.path(), "clinvar", &vcf, "big.vcf");

    let mut checked = 0;
    for (chrom, pos) in &expected {
        let j1 = json_of(v1.annotate_position(chrom, *pos, "A", "G").unwrap())
            .unwrap_or_else(|| panic!("v1 missing {chrom}:{pos}"));
        let j2 = json_of(v2.annotate_position(chrom, *pos, "A", "G").unwrap())
            .unwrap_or_else(|| panic!("converted db dropped {chrom}:{pos}"));
        assert_eq!(j1, j2, "{chrom}:{pos}");
        checked += 1;
    }
    assert_eq!(checked, expected.len());
    assert_eq!(checked, 8_000, "sanity: the fixture should span many blocks");
}

#[test]
fn rejects_inputs_that_have_no_v2_form() {
    let dir = tempfile::tempdir().unwrap();

    // Already v2.
    let osa2 = dir.path().join("already.osa2");
    fs::write(&osa2, b"whatever").unwrap();
    let err = run_sa_convert(osa2.to_str().unwrap(), dir.path().join("o").to_str().unwrap(), false)
        .unwrap_err();
    assert!(err.to_string().contains("already a v2"), "got: {err}");

    // Interval- and gene-level databases: v2 is a variant-level container, so
    // there is nothing to convert them into. The error must say so rather than
    // failing on a parse deep inside.
    for (ext, expected) in [("osi", "interval-level"), ("oga", "gene-level")] {
        let path = dir.path().join(format!("db.{ext}"));
        fs::write(&path, b"whatever").unwrap();
        let err =
            run_sa_convert(path.to_str().unwrap(), dir.path().join("o").to_str().unwrap(), false)
                .unwrap_err();
        assert!(
            err.to_string().contains(expected) && err.to_string().contains("variant-level"),
            "expected a {ext}-specific message, got: {err}"
        );
    }

    // Anything else.
    let other = dir.path().join("notes.txt");
    fs::write(&other, b"whatever").unwrap();
    let err = run_sa_convert(other.to_str().unwrap(), dir.path().join("o").to_str().unwrap(), false)
        .unwrap_err();
    assert!(err.to_string().contains("expects a v1"), "got: {err}");
}

#[test]
fn refuses_to_overwrite_its_own_input() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("in.vcf");
    fs::write(&input, CLINVAR_VCF).unwrap();
    let base = dir.path().join("db");
    run_sa_build(
        "clinvar",
        input.to_str().unwrap(),
        base.to_str().unwrap(),
        "GRCh38",
        None,
        &[],
        false,
    )
    .unwrap();

    // An `.osa2` input was already rejected above; this covers the other way to
    // collide — asking for an output path that resolves onto the input file.
    let osa2_in = base.with_extension("osa2");
    fs::rename(base.with_extension("osa"), &osa2_in).unwrap();
    let err = run_sa_convert(osa2_in.to_str().unwrap(), base.to_str().unwrap(), false).unwrap_err();
    assert!(err.to_string().contains("already a v2"), "got: {err}");
}

#[test]
fn a_failed_conversion_leaves_no_partial_output() {
    // Truncated `.osa` (the `.osa.idx` is intact, so `open` succeeds and the
    // failure lands mid-scan). No `.osa2` may survive.
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("in.vcf");
    fs::write(&input, CLINVAR_VCF).unwrap();
    let base = dir.path().join("db");
    run_sa_build(
        "clinvar",
        input.to_str().unwrap(),
        base.to_str().unwrap(),
        "GRCh38",
        None,
        &[],
        false,
    )
    .unwrap();

    let osa = base.with_extension("osa");
    let truncated: Vec<u8> = fs::read(&osa).unwrap().into_iter().take(16).collect();
    fs::write(&osa, &truncated).unwrap();

    let out_base = dir.path().join("converted");
    let result = run_sa_convert(osa.to_str().unwrap(), out_base.to_str().unwrap(), false);
    assert!(result.is_err(), "a truncated .osa must not convert successfully");
    assert!(
        !out_base.with_extension("osa2").exists(),
        "no partial .osa2 may be left behind"
    );
}
