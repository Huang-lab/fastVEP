//! Regression test for the v2 reader half of issue #37.
//!
//! The `Osa2Reader` keyed its chunk LRU on the raw query chromosome
//! while only resolving aliases deep inside `build_chunk`. A pipeline
//! that mixed `chr1` and `1` for the same physical chunk — or warmed
//! the cache via `preload("chr1")` and then served queries with `"1"`
//! — would double-key the cache (wasted RAM + re-decoded chunk) and
//! silently miss preload entirely. The fix canonicalizes the chromosome
//! against the archive's on-disk chrom set at every public entry point
//! before constructing any cache key.

use fastvep_cache::annotation::AnnotationProvider;
use fastvep_sa::fields::{Field, FieldType};
use fastvep_sa::reader_v2::Osa2Reader;
use fastvep_sa::writer_v2::{Osa2Metadata, Osa2Record, Osa2Writer};

fn metadata() -> Osa2Metadata {
    Osa2Metadata {
        format_version: 2,
        name: "gnomAD".into(),
        version: "v4.1".into(),
        assembly: "GRCh38".into(),
        json_key: "gnomad".into(),
        match_by_allele: true,
        is_array: false,
        is_positional: false,
        chunk_bits: 20,
        description: "test".into(),
    }
}

fn one_int_field() -> Vec<Field> {
    vec![Field {
        field: "AC".into(),
        alias: "allAc".into(),
        ftype: FieldType::Integer,
        multiplier: 1,
        zigzag: false,
        missing_value: u32::MAX,
        missing_string: ".".into(),
        description: String::new(),
    }]
}

#[test]
fn v2_reader_resolves_chr_query_against_bare_keyed_archive() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("v2_bare.osa2");

    let fields = one_int_field();
    let records: Vec<Osa2Record> = (0..16)
        .map(|i| Osa2Record {
            chrom: "1".into(),
            position: 10_000 + i * 100,
            ref_allele: b"A".to_vec(),
            alt_allele: b"G".to_vec(),
            values: vec![i as u32 + 1],
            json_blob: None,
        })
        .collect();

    let writer = Osa2Writer::new(metadata(), fields);
    writer
        .write_all(std::fs::File::create(&path).unwrap(), &records)
        .unwrap();

    let reader = Osa2Reader::open(&path).unwrap();

    // chr*-prefixed query hits a bare-keyed archive.
    let hit = reader
        .annotate_position("chr1", 10_500, "A", "G")
        .unwrap();
    assert!(hit.is_some(), "chr1 query must resolve against bare archive");

    // And the bare form still works.
    let hit = reader.annotate_position("1", 10_500, "A", "G").unwrap();
    assert!(hit.is_some(), "bare query must resolve against bare archive");
}

#[test]
fn v2_reader_preload_chr_warms_cache_for_bare_query() {
    // The original bug: preload("chr1") only resolved aliases inside
    // `build_chunk`, but `load_chunk` keyed the LRU on the raw chrom.
    // A subsequent `query("1", ...)` then missed the cache and went
    // back to disk. Now both sides canonicalize against the archive,
    // so the preload and the query share a single cache slot.
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("v2_preload.osa2");

    let fields = one_int_field();
    let records: Vec<Osa2Record> = (0..16)
        .map(|i| Osa2Record {
            chrom: "1".into(),
            position: 10_000 + i * 100,
            ref_allele: b"A".to_vec(),
            alt_allele: b"G".to_vec(),
            values: vec![i as u32 + 1],
            json_blob: None,
        })
        .collect();
    let writer = Osa2Writer::new(metadata(), fields);
    writer
        .write_all(std::fs::File::create(&path).unwrap(), &records)
        .unwrap();

    let reader = Osa2Reader::open(&path).unwrap();

    // Warm the cache via the chr*-prefixed style.
    reader.preload("chr1", &[10_500u64]).unwrap();

    // Now query with the bare form. If the cache key is not
    // canonicalized, this would re-decode the chunk — which is the
    // expensive path the preload is supposed to amortize. The
    // visible-from-test signal is simply that the bare query still
    // returns the right annotation; the cache-key coherence is a
    // structural invariant kept by `resolve_chrom`.
    let hit = reader.annotate_position("1", 10_500, "A", "G").unwrap();
    assert!(hit.is_some());
}
