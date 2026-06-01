//! Regression tests for issue #37: `chr*` vs bare contig naming.
//!
//! Before the fix, `standard_chrom_map()` declared the canonical contig
//! list as `["1","2",…,"22","X","Y","MT"]`, so every `.osa.idx` was keyed
//! by the bare name even when the source VCF used `chr*`. Reader queries
//! coming from a modern (chr-prefixed) input VCF then missed the index
//! and returned `null` for every gnomAD annotation. The fix has two
//! halves: (a) the writer now canonicalizes to `chr*`, and (b) the
//! reader walks aliases so existing bare-keyed `.osa` files still serve
//! `chr*` queries without a rebuild.

use fastvep_cache::annotation::AnnotationProvider;
use fastvep_sa::common::AnnotationRecord;
use fastvep_sa::index::{IndexHeader, SaIndex};
use fastvep_sa::reader::SaReader;
use fastvep_sa::writer::SaWriter;

fn gnomad_header() -> IndexHeader {
    IndexHeader {
        schema_version: fastvep_sa::common::SCHEMA_VERSION,
        json_key: "gnomad".into(),
        name: "gnomAD".into(),
        version: "v4.1".into(),
        description: "Test".into(),
        assembly: "GRCh38".into(),
        match_by_allele: true,
        is_array: false,
        is_positional: false,
    }
}

/// A `chr*`-style query must hit an index built with the bare-style key
/// (the broken legacy state) AND a bare-style query must continue to
/// resolve, so existing on-disk databases keep serving traffic without a
/// rebuild.
#[test]
fn reader_tolerates_chr_vs_bare_naming() {
    let records = vec![AnnotationRecord {
        chrom_idx: 0,
        position: 10266,
        ref_allele: "A".into(),
        alt_allele: "G".into(),
        json: r#"{"allAf":1.0e-5}"#.into(),
    }];

    // Case A: index stores bare-style key ("10") — the pre-fix state.
    {
        let dir = tempfile::tempdir().unwrap();
        let base = dir.path().join("legacy_bare");
        let mut writer = SaWriter::new(gnomad_header());
        writer
            .write_to_files(&base, records.clone().into_iter(), &["10".to_string()])
            .unwrap();
        let reader = SaReader::open(&base.with_extension("osa")).unwrap();

        assert!(
            reader
                .annotate_position("chr10", 10266, "A", "G")
                .unwrap()
                .is_some(),
            "chr10 query must hit a bare-keyed index"
        );
        assert!(
            reader
                .annotate_position("10", 10266, "A", "G")
                .unwrap()
                .is_some(),
            "bare query must hit a bare-keyed index"
        );

        // preload() uses the same alias logic.
        reader.preload("chr10", &[10266]).unwrap();
    }

    // Case B: index stores canonical chr*-style key — the new default
    // for fresh builds.
    {
        let dir = tempfile::tempdir().unwrap();
        let base = dir.path().join("canonical_chr");
        let mut writer = SaWriter::new(gnomad_header());
        writer
            .write_to_files(&base, records.into_iter(), &["chr10".to_string()])
            .unwrap();
        let reader = SaReader::open(&base.with_extension("osa")).unwrap();

        assert!(
            reader
                .annotate_position("chr10", 10266, "A", "G")
                .unwrap()
                .is_some(),
            "chr10 query must hit a chr-keyed index"
        );
        assert!(
            reader
                .annotate_position("10", 10266, "A", "G")
                .unwrap()
                .is_some(),
            "bare query must hit a chr-keyed index"
        );
    }
}

/// Mitochondrial contigs collect four aliases in the wild (UCSC `chrM`,
/// NCBI `MT`, plus the bare `M` and the non-standard `chrMT`). All four
/// must round-trip in both directions so MT variants don't drop
/// depending on which convention the source and input happen to use.
#[test]
fn mitochondrial_aliases_round_trip() {
    let records = vec![AnnotationRecord {
        chrom_idx: 0,
        position: 3243,
        ref_allele: "A".into(),
        alt_allele: "G".into(),
        json: r#"{"allAf":1.0e-4}"#.into(),
    }];

    for stored_as in ["chrM", "MT"] {
        let dir = tempfile::tempdir().unwrap();
        let base = dir.path().join(format!("mito_{}", stored_as));
        let mut writer = SaWriter::new(gnomad_header());
        writer
            .write_to_files(
                &base,
                records.clone().into_iter(),
                &[stored_as.to_string()],
            )
            .unwrap();
        let reader = SaReader::open(&base.with_extension("osa")).unwrap();

        for queried_as in ["chrM", "M", "MT", "chrMT"] {
            assert!(
                reader
                    .annotate_position(queried_as, 3243, "A", "G")
                    .unwrap()
                    .is_some(),
                "query {} must hit index keyed by {}",
                queried_as,
                stored_as
            );
        }
    }
}

/// `find_blocks_range` shares the new alias-resolving lookup with
/// `find_blocks`. The range path is used by interval-style providers,
/// so a regression on its side wouldn't be caught by the
/// annotate-position tests above.
#[test]
fn find_blocks_range_uses_alias_resolution() {
    let dir = tempfile::tempdir().unwrap();
    let base = dir.path().join("range");
    let records: Vec<AnnotationRecord> = (0..10)
        .map(|i| AnnotationRecord {
            chrom_idx: 0,
            position: 1_000 + i * 50,
            ref_allele: "A".into(),
            alt_allele: "G".into(),
            json: format!(r#"{{"i":{}}}"#, i),
        })
        .collect();
    let mut writer = SaWriter::new(gnomad_header());
    writer
        .write_to_files(&base, records.into_iter(), &["chr5".to_string()])
        .unwrap();

    // Re-open the index directly so we can probe the range lookup.
    let idx_path = base.with_extension("osa.idx");
    let mut f = std::fs::File::open(&idx_path).unwrap();
    let index = SaIndex::read_from(&mut f).unwrap();

    let by_chr = index.find_blocks_range("chr5", 1_000, 2_000);
    let by_bare = index.find_blocks_range("5", 1_000, 2_000);
    assert!(!by_chr.is_empty(), "chr5 range query must return blocks");
    assert_eq!(
        by_chr.len(),
        by_bare.len(),
        "chr5 and 5 must return the same range blocks"
    );

    // A truly unknown contig must still return empty.
    assert!(index.find_blocks_range("chrZZ", 1_000, 2_000).is_empty());
}
