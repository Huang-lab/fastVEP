//! `--cache-dir`: reading a VEP variation cache the way VEP converts it, and
//! answering the questions issue #139 raised about it.
//!
//! The fixtures under `fixtures/vep_cache/` are made by `build.sh` with `bgzip`
//! and the same `tabix` command `convert_cache.pl` runs, so the `.csi` the first
//! tests read is the real thing and not one written by the library under test.
//! Where an expectation about *which known variants match an input* appears
//! below, the allele logic is VEP's own `get_matched_variant_alleles` (the
//! differential test is `vep_match_golden.rs`); the candidate window, the
//! failed-variant filter, the allele-less record rule and the order of the
//! output are VEP's as read in ensembl-vep release/115, and each is named where
//! it is checked.

use fastvep_cache::info::CacheInfo;
use fastvep_cache::providers::{MatchedVariant, TabixVariationProvider, VariationProvider};
use fastvep_cache::variation::{VariationCacheRecord, VariationTabixReader};
use std::path::{Path, PathBuf};

fn fixture(kind: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("fixtures/vep_cache")
        .join(kind)
}

fn info(dir: &Path) -> CacheInfo {
    CacheInfo::from_file(&dir.join("info.txt")).expect("the fixture has an info.txt")
}

fn open_dir(dir: &Path) -> anyhow::Result<VariationTabixReader> {
    VariationTabixReader::new(dir, &info(dir).variation_cols)
}

fn open(kind: &str) -> VariationTabixReader {
    open_dir(&fixture(kind)).unwrap_or_else(|e| panic!("opening fixture {kind}: {e:#}"))
}

fn open_error(dir: &Path) -> String {
    match open_dir(dir) {
        Ok(_) => panic!("{} should have been refused", dir.display()),
        Err(e) => format!("{e:#}"),
    }
}

fn names(records: &[VariationCacheRecord]) -> Vec<&str> {
    records.iter().map(|r| r.variation_name.as_str()).collect()
}

fn copy_dir(from: &Path, to: &Path) {
    std::fs::create_dir_all(to).unwrap();
    for entry in std::fs::read_dir(from).unwrap() {
        let entry = entry.unwrap();
        let target = to.join(entry.file_name());
        if entry.file_type().unwrap().is_dir() {
            copy_dir(&entry.path(), &target);
        } else {
            std::fs::copy(entry.path(), target).unwrap();
        }
    }
}

// -- the index ---------------------------------------------------------------

/// Issue #139, claim 4. `convert_cache.pl` indexes with `tabix -C`, which
/// writes a `.csi`; a cache VEP converted has no `.tbi`. The reader used to look
/// for a `.tbi` only, found none, and answered every lookup with "nothing known".
#[test]
fn a_cache_with_only_a_csi_index_is_searched() {
    let dir = fixture("csi");
    assert!(dir.join("10/all_vars.gz.csi").is_file());
    assert!(
        !dir.join("10/all_vars.gz.tbi").exists(),
        "the premise: no .tbi to fall back on"
    );

    let reader = open("csi");
    let found = reader.query("10", 10400, 10401).unwrap();
    assert_eq!(names(&found), ["rsFIX_DEL_MIN", "rsFIX_DEL_ANCH"]);
}

#[test]
fn a_cache_with_a_tbi_index_still_works_and_answers_the_same() {
    let csi = open("csi");
    let tbi = open("tbi");
    for (lo, hi) in [
        (10400, 10401),
        (10408, 10409),
        (10414, 10414),
        (10419, 10421),
        (1, 20000),
    ] {
        let from_csi = csi.query("10", lo, hi).unwrap();
        let from_tbi = tbi.query("10", lo, hi).unwrap();
        assert_eq!(names(&from_csi), names(&from_tbi), "window {lo}..={hi}");
        assert!(
            !from_csi.is_empty(),
            "window {lo}..={hi} should not be empty"
        );
    }
}

/// A directory that holds both indexes is read through the `.csi`, the one the
/// `tabix` command prefers (checked with tabix 1.23.1). In `both/` the `.csi` is
/// right and the `.tbi` is built on the end column, so reading the `.tbi` would
/// be refused.
#[test]
fn a_directory_with_both_indexes_is_read_through_the_csi() {
    let dir = fixture("both");
    assert!(dir.join("10/all_vars.gz.csi").is_file());
    assert!(dir.join("10/all_vars.gz.tbi").is_file());

    let reader = open("both");
    let found = reader.query("10", 10400, 10401).unwrap();
    assert_eq!(names(&found), ["rsFIX_DEL_MIN", "rsFIX_DEL_ANCH"]);
}

/// A variation file that cannot be searched has to be refused when the cache is
/// opened, with the file named and the command that fixes it. Reading it as "no
/// known variants" is the failure the whole issue is about.
#[test]
fn a_variation_file_without_an_index_is_refused_when_the_cache_is_opened() {
    let message = open_error(&fixture("noindex"));
    assert!(message.contains("all_vars.gz"), "{message}");
    assert!(message.contains("no .csi or .tbi index"), "{message}");
    assert!(
        message.contains("tabix -C -s 1 -b 5 -e 5"),
        "the message gives the command convert_cache.pl runs: {message}"
    );
}

/// An index built on the end column is a valid tabix index and answers queries,
/// with the wrong records. info.txt says where `start` is; the index says what
/// it was built on; they have to agree.
#[test]
fn an_index_built_on_another_column_is_refused() {
    let message = open_error(&fixture("endcol"));
    assert!(message.contains("indexes columns"), "{message}");
    assert!(message.contains("tabix -C -s 1 -b 5 -e 5"), "{message}");
}

#[test]
fn a_directory_without_variation_files_is_refused() {
    let dir = tempfile::tempdir().unwrap();
    std::fs::create_dir(dir.path().join("1")).unwrap();
    std::fs::copy(fixture("csi").join("info.txt"), dir.path().join("info.txt")).unwrap();
    let message = open_error(dir.path());
    assert!(
        message.contains("holds no `{chr}/all_vars.gz`"),
        "{message}"
    );
    assert!(
        message.contains("unconverted"),
        "it says what an unconverted cache is: {message}"
    );
}

/// The chromosomes are the directories, as in VEP. An info.txt that has a
/// `valid_chromosomes` line, however it is written, does not change which
/// directories are read or make a good cache unreadable.
#[test]
fn a_valid_chromosomes_line_in_info_txt_is_not_consulted() {
    let dir = tempfile::tempdir().unwrap();
    copy_dir(&fixture("csi"), dir.path());
    let info_path = dir.path().join("info.txt");
    let mut text = std::fs::read_to_string(&info_path).unwrap();
    // tab-separated, which info.rs would read as a single name
    text.push_str("valid_chromosomes\t10\tMT\n");
    std::fs::write(&info_path, text).unwrap();

    let reader = open_dir(dir.path()).expect("the directories decide, not the line");
    assert_eq!(reader.chromosomes(), ["10", "MT"]);
    assert_eq!(
        names(&reader.query("chrM", 100, 100).unwrap()),
        ["rsFIX_MT"]
    );
}

#[test]
fn a_corrupt_index_is_refused_when_the_cache_is_opened() {
    let dir = tempfile::tempdir().unwrap();
    copy_dir(&fixture("csi"), dir.path());
    std::fs::write(dir.path().join("10/all_vars.gz.csi"), b"not an index").unwrap();
    let message = open_error(dir.path());
    assert!(message.contains("all_vars.gz.csi"), "{message}");
}

// -- the window --------------------------------------------------------------

/// `convert_cache.pl` indexes the start column alone, and VEP asks for
/// `start - 1 ..= end + 1` of the input, so a record is found by where it
/// starts. Its own end plays no part.
#[test]
fn a_record_is_found_by_its_start_and_not_by_overlap() {
    let reader = open("csi");

    // An insertion is stored with end below start. It was lost to an overlap
    // test (`record.end < input_start`).
    let found = reader.query("10", 10408, 10409).unwrap();
    assert_eq!(names(&found), ["rsFIX_INS_MIN"]);
    assert_eq!((found[0].start, found[0].end), (10409, 10408));

    // The anchored deletion spans 10401..=10406 but starts at 10401: a window
    // inside its span and past its start does not contain it.
    assert!(reader.query("10", 10405, 10408).unwrap().is_empty());

    // Records that share a start come back in file order, and the reader does
    // not filter: that is the provider's job.
    let at_10414 = reader.query("10", 10414, 10414).unwrap();
    assert_eq!(
        names(&at_10414),
        [
            "COSVFIX1",
            "rsFIX_SNV",
            "rsFIX_SNV_FAILED",
            "rsFIX_SNV_MULTI"
        ]
    );

    assert!(reader.query("10", 5000, 5100).unwrap().is_empty());
    // an inverted window is empty, not an error
    assert!(reader.query("10", 10414, 10413).unwrap().is_empty());
}

#[test]
fn records_carry_what_the_cache_line_says() {
    let reader = open("csi");
    let at_10414 = reader.query("10", 10414, 10414).unwrap();
    let by_name = |n: &str| at_10414.iter().find(|r| r.variation_name == n).unwrap();

    let snv = by_name("rsFIX_SNV");
    assert_eq!(snv.allele_string, "C/A");
    assert_eq!(snv.strand, 1);
    assert_eq!(snv.minor_allele.as_deref(), Some("A"));
    assert_eq!(snv.minor_allele_freq, Some(0.002));
    assert_eq!(snv.clin_sig.as_deref(), Some("benign"));
    assert_eq!(
        snv.frequencies.get("AFR").map(String::as_str),
        Some("A:0.0100")
    );
    assert!(!snv.frequencies.contains_key("AMR"), "`.` is no value");
    assert!(!snv.failed && !snv.somatic);

    assert!(by_name("rsFIX_SNV_FAILED").failed);
    assert!(by_name("COSVFIX1").somatic);
    assert_eq!(by_name("COSVFIX1").allele_string, "COSMIC_MUTATION");
}

// -- contigs -----------------------------------------------------------------

/// The input says `chrM`; the cache directory is `MT`. VEP loads a
/// `chr_synonyms.txt` found in the cache directory; the mechanical `chr` and
/// mitochondrial rules cover the rest.
#[test]
fn contig_names_resolve_to_the_cache_directories() {
    let with_file = open("csi");
    for (query, want) in [
        ("10", Some("10")),
        ("chr10", Some("10")),
        ("NC_000010.11", Some("10")),
        ("MT", Some("MT")),
        ("chrM", Some("MT")),
        ("M", Some("MT")),
        ("NC_012920.1", Some("MT")),
        ("chr11", None),
        ("chrUn_KI270742v1", None),
    ] {
        assert_eq!(with_file.normalize_chrom(query), want, "query {query}");
    }
    assert_eq!(with_file.chromosomes(), ["10", "MT"]);

    // No synonyms file: the mechanical rules still reach MT from chrM, and an
    // accession, which has no rule, is simply unknown.
    let without_file = open("tbi");
    assert_eq!(without_file.normalize_chrom("chrM"), Some("MT"));
    assert_eq!(without_file.normalize_chrom("NC_012920.1"), None);

    // And a query by any spelling reads the right file.
    let by_alias = with_file.query("chrM", 100, 100).unwrap();
    assert_eq!(names(&by_alias), ["rsFIX_MT"]);
    assert!(with_file
        .query("chrUn_KI270742v1", 1, 1_000_000)
        .unwrap()
        .is_empty());
}

// -- the provider ------------------------------------------------------------

fn provider(kind: &str) -> TabixVariationProvider {
    let dir = fixture(kind);
    TabixVariationProvider::new(&dir, &info(&dir)).expect("the fixture cache opens")
}

fn alts(a: &[&str]) -> Vec<String> {
    a.iter().map(ToString::to_string).collect()
}

fn lookup(
    p: &TabixVariationProvider,
    chrom: &str,
    (start, end): (u64, u64),
    reference: &str,
    alt: &[&str],
) -> Vec<Vec<String>> {
    p.get_matched_variants(chrom, start, end, reference, &alts(alt))
        .unwrap_or_else(|e| panic!("lookup at {chrom}:{start}: {e:#}"))
        .iter()
        .map(|per_alt| per_alt.iter().map(|m| m.name.clone()).collect())
        .collect()
}

fn strs(v: &[&[&str]]) -> Vec<Vec<String>> {
    v.iter().map(|names| alts(names)).collect()
}

/// The shapes issue #139 reports, as a site and its known variants. Each
/// answer is what VEP's `--check_existing` lists, in the order it lists them.
#[test]
fn each_shape_the_issue_reports_finds_its_known_variants() {
    let p = provider("csi");

    // `AACCCT>A` at 10400, minimised to ACCCT/- at 10401; the cache has the
    // deletion anchored (ACCCTA/A). The minimal deletion that is a different
    // sequence at 10400 is not a match for it, nor it for that.
    assert_eq!(
        lookup(&p, "chr10", (10401, 10405), "ACCCT", &["-"]),
        strs(&[&["rsFIX_DEL_ANCH"]])
    );
    assert_eq!(
        lookup(&p, "chr10", (10400, 10404), "AACCC", &["-"]),
        strs(&[&["rsFIX_DEL_MIN"]])
    );

    // An insertion stored with end < start, and one stored anchored.
    assert_eq!(
        lookup(&p, "chr10", (10409, 10408), "-", &["T"]),
        strs(&[&["rsFIX_INS_MIN"]])
    );
    assert_eq!(
        lookup(&p, "chr10", (10421, 10420), "-", &["A"]),
        strs(&[&["rsFIX_INS_ANCH"]])
    );

    // The SNV site. Four records start there: COSVFIX1 (somatic, no alleles),
    // rsFIX_SNV (C/A), rsFIX_SNV_FAILED (C/A, failed) and rsFIX_SNV_MULTI
    // (C/A/G). The failed one is never listed. The record without alleles
    // matches any allele at the same start and end. Germline comes before
    // somatic whatever the file order is.
    assert_eq!(
        lookup(&p, "chr10", (10414, 10414), "C", &["A"]),
        strs(&[&["rsFIX_SNV", "rsFIX_SNV_MULTI", "COSVFIX1"]])
    );
    assert_eq!(
        lookup(&p, "chr10", (10414, 10414), "C", &["G"]),
        strs(&[&["rsFIX_SNV_MULTI", "COSVFIX1"]])
    );
    assert_eq!(
        lookup(&p, "chr10", (10414, 10414), "C", &["T"]),
        strs(&[&["COSVFIX1"]])
    );

    // One record, two alleles: each allele gets its own list.
    assert_eq!(
        lookup(&p, "chr10", (10414, 10414), "C", &["A", "G"]),
        strs(&[
            &["rsFIX_SNV", "rsFIX_SNV_MULTI", "COSVFIX1"],
            &["rsFIX_SNV_MULTI", "COSVFIX1"]
        ])
    );

    // The mitochondrion, by the input's own spelling.
    assert_eq!(
        lookup(&p, "chrM", (100, 100), "A", &["G"]),
        strs(&[&["rsFIX_MT"]])
    );

    // Nothing there.
    assert_eq!(
        lookup(&p, "chr10", (10500, 10500), "A", &["C"]),
        strs(&[&[]])
    );
}

#[test]
fn the_tbi_form_of_the_cache_answers_identically() {
    let (csi, tbi) = (provider("csi"), provider("tbi"));
    for (chrom, span, reference, alt) in [
        ("chr10", (10401, 10405), "ACCCT", vec!["-"]),
        ("chr10", (10409, 10408), "-", vec!["T"]),
        ("chr10", (10414, 10414), "C", vec!["A", "G"]),
        ("chrM", (100, 100), "A", vec!["G"]),
    ] {
        assert_eq!(
            lookup(&csi, chrom, span, reference, &alt),
            lookup(&tbi, chrom, span, reference, &alt)
        );
    }
}

#[test]
fn a_matched_variant_carries_its_record_and_the_frequencies_of_the_matched_allele() {
    let p = provider("csi");
    let found: Vec<MatchedVariant> = p
        .get_matched_variants("chr10", 10414, 10414, "C", &alts(&["A"]))
        .unwrap()
        .remove(0);

    let snv = &found[0];
    assert_eq!(snv.name, "rsFIX_SNV");
    assert_eq!(snv.matched_allele, "A");
    assert_eq!(snv.minor_allele.as_deref(), Some("A"));
    assert_eq!(snv.minor_allele_freq, Some(0.002));
    assert_eq!(snv.clin_sig.as_deref(), Some("benign"));
    assert_eq!(snv.frequencies.get("AFR"), Some(&0.01));
    assert_eq!(snv.frequencies.get("minor_allele_freq"), Some(&0.002));

    // For an allele the record has no frequency for, there is none to report.
    let g = p
        .get_matched_variants("chr10", 10414, 10414, "C", &alts(&["G"]))
        .unwrap()
        .remove(0);
    let multi = g.iter().find(|m| m.name == "rsFIX_SNV_MULTI").unwrap();
    assert_eq!(multi.matched_allele, "G");
    assert!(multi.frequencies.is_empty());

    // The record without alleles has no matched allele to key a frequency by.
    let cosmic = found.iter().find(|m| m.name == "COSVFIX1").unwrap();
    assert!(cosmic.somatic);
    assert_eq!(cosmic.matched_allele, "");
}

// -- what the run tells the user ---------------------------------------------

#[test]
fn a_run_that_found_no_contig_in_the_cache_says_so() {
    let p = provider("csi");
    assert_eq!(p.finish(), None, "nothing looked up, nothing to say");

    lookup(&p, "chr10", (10414, 10414), "C", &["A"]);
    assert_eq!(p.finish(), None, "every contig resolved");

    // Two sites on a contig the cache does not hold.
    lookup(&p, "chrUn_KI270742v1", (500, 500), "A", &["C"]);
    lookup(&p, "chrUn_KI270742v1", (900, 900), "A", &["C"]);
    let partial = p.finish().expect("some sites were on unknown contigs");
    assert!(
        partial.starts_with("warning: 2 of 3 variants are on contigs"),
        "{partial}"
    );
    assert!(partial.contains("chrUn_KI270742v1 (2)"), "{partial}");
    assert!(
        partial.contains("Existing_variation is empty for those"),
        "{partial}"
    );
    assert!(
        partial.contains("The cache has 10, MT"),
        "the cache's contigs: {partial}"
    );
    assert!(
        partial.contains("chr_synonyms.txt"),
        "the way out: {partial}"
    );
}

#[test]
fn a_run_with_no_resolvable_contig_at_all_says_that_instead() {
    let p = provider("csi");
    lookup(&p, "NC_000010.99", (10414, 10414), "C", &["A"]);
    let all = p.finish().unwrap();
    assert!(
        all.contains("none of the input's variants is on a contig"),
        "{all}"
    );
    assert!(all.contains("empty for all of them"), "{all}");
}

// -- errors are errors -------------------------------------------------------

/// The lookup used to swallow its errors (`.unwrap_or_default()` at the call
/// site), so a cache that could not be read looked like a cache with nothing in
/// it. A lookup that cannot be done is an error, and only the chromosome it
/// concerns.
#[test]
fn a_lookup_that_cannot_read_the_cache_is_an_error() {
    let dir = tempfile::tempdir().unwrap();
    copy_dir(&fixture("csi"), dir.path());
    let p = TabixVariationProvider::new(dir.path(), &info(dir.path())).unwrap();

    // The index goes bad after the cache was validated and opened.
    std::fs::write(dir.path().join("10/all_vars.gz.csi"), b"damaged").unwrap();

    let err = p
        .get_matched_variants("chr10", 10414, 10414, "C", &alts(&["A"]))
        .unwrap_err();
    let message = format!("{err:#}");
    assert!(message.contains("all_vars.gz.csi"), "{message}");

    // Another chromosome's file is unaffected.
    assert_eq!(
        lookup(&p, "chrM", (100, 100), "A", &["G"]),
        strs(&[&["rsFIX_MT"]])
    );
}

#[test]
fn opening_a_provider_on_a_bad_cache_names_the_directory_and_the_info() {
    let err = match TabixVariationProvider::new(&fixture("noindex"), &info(&fixture("noindex"))) {
        Ok(_) => panic!("a cache without an index is not usable"),
        Err(e) => format!("{e:#}"),
    };
    assert!(err.contains("--cache-dir"), "{err}");
    assert!(err.contains("var_type tabix"), "{err}");
    assert!(err.contains("no .csi or .tbi index"), "{err}");
}
