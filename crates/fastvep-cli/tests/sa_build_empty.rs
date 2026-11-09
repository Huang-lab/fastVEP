//! A database that holds nothing, and the sources that were never loaded.
//!
//! The dbSNP database that was meant to serve GRCh37 on fastvep.org was ten bytes of header. It had
//! been built from a release that names its contigs by RefSeq accession, none of them matched the
//! assembly, and `sa-build` wrote an empty database and exited with success. A server then opened it
//! without a word and answered every query with "no record", which for a frequency source reads as
//! "absent from the population". Two things were wrong and both are tested here: the builder
//! must refuse to write a database of nothing, and a loader must refuse to serve one that already
//! exists (and say so, in a way a CLI user sees).

use fastvep_annotate::{load_gene_providers_reporting, load_sa_providers_reporting};
use fastvep_cli::pipeline::{
    run_sa_build, run_sa_build_format, run_sa_build_v2, run_sa_convert, with_allow_empty,
};
use std::fs;
use std::path::Path;

/// A dbSNP-shaped VCF whose contigs are RefSeq accessions of no assembly the builder knows, which
/// is what an input for the wrong build looks like to the parser.
const DBSNP_FOR_ANOTHER_BUILD: &str = "##fileformat=VCFv4.2\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
NC_999999.1\t100\trs1\tA\tG\t.\t.\t.\n\
NC_999999.1\t200\trs2\tC\tT\t.\t.\t.\n";

const DBSNP: &str = "##fileformat=VCFv4.2\n\
##fileDate=2026-11-02\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
1\t100\trs1\tA\tG\t.\t.\t.\n\
1\t200\trs2\tC\tT\t.\t.\t.\n";

fn write(dir: &Path, name: &str, text: &str) -> String {
    let path = dir.join(name);
    fs::write(&path, text).unwrap();
    path.to_string_lossy().into_owned()
}

fn out(dir: &Path, stem: &str) -> String {
    dir.join(stem).to_string_lossy().into_owned()
}

#[test]
fn a_build_that_finds_no_record_is_an_error_and_leaves_nothing_that_loads_v2() {
    let tmp = tempfile::tempdir().unwrap();
    let input = write(tmp.path(), "dbsnp.vcf", DBSNP_FOR_ANOTHER_BUILD);
    let err = run_sa_build_v2(
        "dbsnp",
        &input,
        &out(tmp.path(), "dbsnp"),
        "GRCh37",
        None,
        &[],
        false,
    )
    .unwrap_err();
    let message = format!("{err:#}");
    assert!(message.contains("no record was written"), "{message}");
    assert!(message.contains("RefSeq accession"), "{message}");
    assert!(message.contains("--allow-empty"), "{message}");
    assert!(!tmp.path().join("dbsnp.osa2").exists());
}

#[test]
fn a_build_that_finds_no_record_is_an_error_and_leaves_nothing_that_loads_v1() {
    let tmp = tempfile::tempdir().unwrap();
    let input = write(tmp.path(), "dbsnp.vcf", DBSNP_FOR_ANOTHER_BUILD);
    // The streaming v1 builders (dbSNP among them) write the file as they go, so the refusal has to
    // take it back.
    let err = run_sa_build_format(
        "osa",
        "dbsnp",
        &input,
        &out(tmp.path(), "dbsnp"),
        "GRCh37",
        None,
        &[],
        false,
    )
    .unwrap_err();
    assert!(format!("{err:#}").contains("no record was written"));
    assert!(!tmp.path().join("dbsnp.osa").exists());
    assert!(!tmp.path().join("dbsnp.osa.idx").exists());

    // A source that is buffered before it is written refuses before it writes.
    let clinvar = write(
        tmp.path(),
        "clinvar.vcf",
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nNC_999999.1\t100\t1\tA\tG\t.\t.\tCLNSIG=Pathogenic\n",
    );
    let err = run_sa_build(
        "clinvar",
        &clinvar,
        &out(tmp.path(), "clinvar"),
        "GRCh38",
        None,
        &[],
        false,
    )
    .unwrap_err();
    assert!(format!("{err:#}").contains("no record was written"));
    assert!(!tmp.path().join("clinvar.osa").exists());
}

#[test]
fn an_empty_database_is_written_only_when_it_is_asked_for_and_the_loader_does_not_serve_it() {
    let tmp = tempfile::tempdir().unwrap();
    let sa = tmp.path().join("sa");
    fs::create_dir(&sa).unwrap();
    let input = write(tmp.path(), "dbsnp.vcf", DBSNP_FOR_ANOTHER_BUILD);
    with_allow_empty(true, || {
        run_sa_build_v2(
            "dbsnp",
            &input,
            &out(&sa, "dbsnp"),
            "GRCh37",
            None,
            &[],
            false,
        )
    })
    .expect("an empty database can be asked for");
    assert!(sa.join("dbsnp.osa2").exists());
    // The flag is for that call alone.
    assert!(run_sa_build_v2(
        "dbsnp",
        &input,
        &out(&sa, "again"),
        "GRCh37",
        None,
        &[],
        false
    )
    .is_err());

    // And a real one beside it.
    let good = write(tmp.path(), "good.vcf", DBSNP);
    run_sa_build_v2(
        "dbsnp",
        &good,
        &out(&sa, "dbsnp_real"),
        "GRCh38",
        None,
        &[],
        false,
    )
    .unwrap();

    let (providers, skipped) = load_sa_providers_reporting(&sa).unwrap();
    assert_eq!(providers.len(), 1, "{skipped:?}");
    assert_eq!(providers[0].json_key(), "dbsnp");
    assert!(!providers[0].is_empty());
    assert_eq!(skipped.len(), 1);
    assert!(skipped[0].path.ends_with("dbsnp.osa2"));
    assert!(
        skipped[0].reason.contains("holds no records"),
        "{}",
        skipped[0].reason
    );
    assert!(skipped[0].to_string().contains("was not loaded"));
}

#[test]
fn an_empty_v1_database_that_already_exists_is_not_served_either() {
    // What the hosted site had: written by an older builder, so the only way to get one now is to ask.
    let tmp = tempfile::tempdir().unwrap();
    let sa = tmp.path().join("sa");
    fs::create_dir(&sa).unwrap();
    let input = write(tmp.path(), "dbsnp.vcf", DBSNP_FOR_ANOTHER_BUILD);
    with_allow_empty(true, || {
        run_sa_build_format(
            "osa",
            "dbsnp",
            &input,
            &out(&sa, "dbsnp"),
            "GRCh37",
            None,
            &[],
            false,
        )
    })
    .unwrap();
    assert!(sa.join("dbsnp.osa").exists());

    let (providers, skipped) = load_sa_providers_reporting(&sa).unwrap();
    assert!(providers.is_empty());
    assert_eq!(skipped.len(), 1, "{skipped:?}");
    assert!(skipped[0].reason.contains("holds no records"));

    // Converting it to v2 is the same refusal.
    let err = run_sa_convert(
        sa.join("dbsnp.osa").to_str().unwrap(),
        &out(tmp.path(), "converted"),
        20,
        false,
    )
    .unwrap_err();
    assert!(format!("{err:#}").contains("no record was written"));
}

#[test]
fn a_gene_database_with_no_gene_is_refused_too() {
    let tmp = tempfile::tempdir().unwrap();
    let input = write(
        tmp.path(),
        "genemap.txt",
        "# Generated by test\n# nothing else\n",
    );
    let err = run_sa_build(
        "omim",
        &input,
        &out(tmp.path(), "omim"),
        "GRCh38",
        None,
        &[],
        false,
    )
    .unwrap_err();
    assert!(format!("{err:#}").contains("no record was written"));
    assert!(!tmp.path().join("omim.oga").exists());

    let sa = tmp.path().join("sa");
    fs::create_dir(&sa).unwrap();
    with_allow_empty(true, || {
        run_sa_build(
            "omim",
            &input,
            &out(&sa, "omim"),
            "GRCh38",
            None,
            &[],
            false,
        )
    })
    .unwrap();
    let (genes, skipped) = load_gene_providers_reporting(&sa).unwrap();
    assert!(genes.is_empty());
    assert_eq!(skipped.len(), 1);
    assert!(skipped[0].reason.contains("holds no records"));
}

#[test]
fn a_file_that_cannot_be_opened_is_reported_with_the_reason_and_the_others_are_still_loaded() {
    let tmp = tempfile::tempdir().unwrap();
    let sa = tmp.path().join("sa");
    fs::create_dir(&sa).unwrap();
    let good = write(tmp.path(), "good.vcf", DBSNP);
    run_sa_build_v2(
        "dbsnp",
        &good,
        &out(&sa, "a_good"),
        "GRCh38",
        None,
        &[],
        false,
    )
    .unwrap();
    fs::write(sa.join("b_truncated.osa2"), b"PK\x03\x04 cut off mid-copy").unwrap();
    fs::write(sa.join("notes.txt"), "not a database").unwrap();

    let (providers, skipped) = load_sa_providers_reporting(&sa).unwrap();
    assert_eq!(providers.len(), 1);
    assert_eq!(skipped.len(), 1);
    assert!(skipped[0].path.ends_with("b_truncated.osa2"));
    assert!(!skipped[0].reason.is_empty());
    assert!(
        !skipped[0].reason.contains("holds no records"),
        "{}",
        skipped[0].reason
    );
}

#[test]
fn a_release_that_dates_itself_is_stamped_with_its_date() {
    let tmp = tempfile::tempdir().unwrap();
    let sa = tmp.path().join("sa");
    fs::create_dir(&sa).unwrap();
    let input = write(tmp.path(), "dbsnp.vcf", DBSNP);
    run_sa_build_v2("dbsnp", &input, &out(&sa, "v2"), "GRCh38", None, &[], false).unwrap();
    run_sa_build_format(
        "osa",
        "dbsnp",
        &input,
        &out(&sa, "v1"),
        "GRCh38",
        None,
        &[],
        false,
    )
    .unwrap();
    let (providers, skipped) = load_sa_providers_reporting(&sa).unwrap();
    assert!(skipped.is_empty(), "{skipped:?}");
    assert_eq!(providers.len(), 2);
    for p in &providers {
        assert_eq!(p.metadata().version, "2026-11-02", "{}", p.name());
    }

    // One that does not stays `latest`, as before: nothing is invented.
    let undated = write(
        tmp.path(),
        "undated.vcf",
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n1\t100\trs1\tA\tG\t.\t.\t.\n",
    );
    let other = tmp.path().join("other");
    fs::create_dir(&other).unwrap();
    run_sa_build_v2(
        "dbsnp",
        &undated,
        &out(&other, "undated"),
        "GRCh38",
        None,
        &[],
        false,
    )
    .unwrap();
    let (providers, _) = load_sa_providers_reporting(&other).unwrap();
    assert_eq!(providers[0].metadata().version, "latest");
}
