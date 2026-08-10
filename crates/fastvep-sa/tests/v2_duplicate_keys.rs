//! Regression tests: when several records share one `(position, ref, alt)` key,
//! `.osa2` must resolve to the same one the `.osa` reader does — the first.
//!
//! Found by annotating 15,437 real ClinVar chr22 variants through a real `.osa`
//! directory and the same databases converted to `.osa2`: the REVEL column
//! disagreed (`0.324` vs `0.349`, `0.009` vs `0.042`, …). Real REVEL carries
//! **111,270 duplicate keys on chr22 alone — 6.4% of its 1,776,286 records** —
//! one per transcript/protein change, and SpliceAI carries one per overlapping
//! gene.
//!
//! The v1 reader resolves a position by scanning forward from its first entry,
//! so it always returns the first matching record. A v2 chunk index is a sorted
//! key array resolved by `binary_search`, which on a run of equal keys returns an
//! *arbitrary* member. The writer now keeps only the first record of each
//! duplicate run, which both reproduces v1 exactly and drops records that were
//! unreachable through either reader.

use fastvep_cache::annotation::{AnnotationProvider, AnnotationValue};
use fastvep_sa::fields::{Field, FieldType};
use fastvep_sa::reader_v2::Osa2Reader;
use fastvep_sa::writer_v2::{Osa2Metadata, Osa2Record, Osa2StreamWriter, Osa2Writer};
use std::io::BufWriter;

fn metadata(is_positional: bool) -> Osa2Metadata {
    Osa2Metadata {
        format_version: 2,
        name: "Test".into(),
        version: "1".into(),
        assembly: "GRCh38".into(),
        json_key: "test".into(),
        match_by_allele: !is_positional,
        is_array: false,
        is_positional,
        chunk_bits: 20,
        description: "test".into(),
    }
}

fn blob_fields() -> Vec<Field> {
    fastvep_sa::writer_v2::raw_json_blob_fields()
}

fn one_int_field() -> Vec<Field> {
    vec![Field {
        field: "AC".into(),
        alias: "ac".into(),
        ftype: FieldType::Integer,
        multiplier: 1,
        zigzag: false,
        missing_value: u32::MAX,
        missing_string: ".".into(),
        description: String::new(),
    }]
}

fn json_of(v: Option<AnnotationValue>) -> Option<String> {
    match v {
        Some(AnnotationValue::Json(j)) | Some(AnnotationValue::Positional(j)) => Some(j),
        Some(AnnotationValue::Interval(v)) => Some(format!("[{}]", v.join(","))),
        None => None,
    }
}

fn blob(chrom: &str, position: u32, r: &str, a: &str, payload: &str) -> Osa2Record {
    Osa2Record {
        chrom: chrom.into(),
        position,
        ref_allele: r.as_bytes().to_vec(),
        alt_allele: a.as_bytes().to_vec(),
        values: Vec::new(),
        json_blob: Some(format!("{{\"id\":\"{payload}\"}}")),
    }
}

/// Records in **input order**, with duplicate keys across all three index
/// buckets: short (Var32), long (kmer16), and non-ACGT (verbatim).
fn duplicated() -> Vec<Osa2Record> {
    vec![
        // Short key, three records. Only "first" may ever be returned.
        blob("chr1", 100, "A", "G", "first"),
        blob("chr1", 100, "A", "G", "second"),
        blob("chr1", 100, "A", "G", "third"),
        // A distinct allele at the same position is not a duplicate.
        blob("chr1", 100, "A", "T", "other-allele"),
        // Long key (ref+alt > 4 bases), duplicated.
        blob("chr1", 200, "GATTACA", "G", "long-first"),
        blob("chr1", 200, "GATTACA", "G", "long-second"),
        // Non-ACGT key, duplicated.
        blob("chr1", 300, "A", "N", "n-first"),
        blob("chr1", 300, "A", "N", "n-second"),
        // A different chunk entirely, also duplicated.
        blob("chr1", 5_000_000, "C", "T", "far-first"),
        blob("chr1", 5_000_000, "C", "T", "far-second"),
        // And a different chromosome, so the dedup cannot be keyed globally.
        blob("chr2", 100, "A", "G", "chr2-first"),
        blob("chr2", 100, "A", "G", "chr2-second"),
    ]
}

fn assert_first_wins(reader: &Osa2Reader) {
    for (chrom, pos, r, a, want) in [
        ("chr1", 100u64, "A", "G", "first"),
        ("chr1", 100, "A", "T", "other-allele"),
        ("chr1", 200, "GATTACA", "G", "long-first"),
        ("chr1", 300, "A", "N", "n-first"),
        ("chr1", 5_000_000, "C", "T", "far-first"),
        ("chr2", 100, "A", "G", "chr2-first"),
    ] {
        let json = json_of(reader.annotate_position(chrom, pos, r, a).unwrap())
            .unwrap_or_else(|| panic!("{chrom}:{pos} {r}>{a} missing"));
        assert_eq!(
            json,
            format!("{{\"id\":\"{want}\"}}"),
            "{chrom}:{pos} {r}>{a} must resolve to the first record written"
        );
    }
}

#[test]
fn buffered_writer_keeps_the_first_of_each_duplicate_key() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("t.osa2");
    Osa2Writer::new(metadata(false), blob_fields())
        .write_all(std::fs::File::create(&path).unwrap(), &duplicated())
        .unwrap();
    assert_first_wins(&Osa2Reader::open(&path).unwrap());
}

#[test]
fn streaming_writer_keeps_the_first_of_each_duplicate_key() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("t.osa2");
    let f = std::fs::File::create(&path).unwrap();
    let mut w = Osa2StreamWriter::new(BufWriter::new(f), &metadata(false), blob_fields()).unwrap();
    for r in duplicated() {
        w.push(r).unwrap();
    }
    w.finish().unwrap();
    assert_first_wins(&Osa2Reader::open(&path).unwrap());
}

#[test]
fn deduping_keeps_value_columns_aligned() {
    // Dropping records shifts every later value slot, so the columns and the
    // index have to be rebuilt together. Give each record a distinct value and
    // check the survivors resolve to their own.
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("t.osa2");

    // (position, ref, alt, value) in input order; value 0 marks the first of a run.
    let specs: Vec<(u32, &str, &str, u32)> = vec![
        (100, "A", "G", 10),
        (100, "A", "G", 11), // dropped
        (100, "A", "T", 20),
        (200, "GATTACA", "G", 30),
        (200, "GATTACA", "G", 31), // dropped
        (300, "A", "N", 40),
        (300, "A", "N", 41), // dropped
        (400, "C", "T", 50),
    ];
    let records: Vec<Osa2Record> = specs
        .iter()
        .map(|(pos, r, a, v)| Osa2Record {
            chrom: "chr1".into(),
            position: *pos,
            ref_allele: r.as_bytes().to_vec(),
            alt_allele: a.as_bytes().to_vec(),
            values: vec![*v],
            json_blob: None,
        })
        .collect();
    Osa2Writer::new(metadata(false), one_int_field())
        .write_all(std::fs::File::create(&path).unwrap(), &records)
        .unwrap();

    let reader = Osa2Reader::open(&path).unwrap();
    for (pos, r, a, want) in [
        (100u32, "A", "G", 10),
        (100, "A", "T", 20),
        (200, "GATTACA", "G", 30),
        (300, "A", "N", 40),
        (400, "C", "T", 50),
    ] {
        let json = json_of(reader.annotate_position("chr1", pos as u64, r, a).unwrap())
            .unwrap_or_else(|| panic!("chr1:{pos} {r}>{a} missing"));
        assert_eq!(
            json,
            format!("{{\"ac\":{want}}}"),
            "chr1:{pos} {r}>{a} resolved to a shifted value slot"
        );
    }
}

#[test]
fn positional_sources_also_keep_the_first_record_per_position() {
    // A positional source keys on coordinate alone, so two records at the same
    // position are duplicates regardless of their alleles — same rule applies.
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("t.osa2");
    let records = vec![
        blob("chr1", 100, "", "", "pos-first"),
        blob("chr1", 100, "", "", "pos-second"),
        blob("chr1", 101, "", "", "next-position"),
    ];
    Osa2Writer::new(metadata(true), blob_fields())
        .write_all(std::fs::File::create(&path).unwrap(), &records)
        .unwrap();

    let reader = Osa2Reader::open(&path).unwrap();
    assert_eq!(
        json_of(reader.annotate_position("chr1", 100, "A", "G").unwrap()).unwrap(),
        "{\"id\":\"pos-first\"}"
    );
    assert_eq!(
        json_of(reader.annotate_position("chr1", 101, "A", "G").unwrap()).unwrap(),
        "{\"id\":\"next-position\"}"
    );
}

#[test]
fn dedup_shrinks_the_archive() {
    // The dropped records were unreachable through either reader's API, so
    // removing them is pure savings — worth pinning so a future change that
    // reintroduces them is noticed.
    let dir = tempfile::tempdir().unwrap();
    let with_dupes = dir.path().join("dupes.osa2");
    let without = dir.path().join("unique.osa2");

    let mut many: Vec<Osa2Record> = Vec::new();
    let mut unique: Vec<Osa2Record> = Vec::new();
    for i in 0..5_000u32 {
        let pos = 1_000 + i * 7;
        // Five records per key, as REVEL has per transcript.
        for k in 0..5 {
            many.push(blob("chr1", pos, "A", "G", &format!("rec-{i}-{k}")));
        }
        unique.push(blob("chr1", pos, "A", "G", &format!("rec-{i}-0")));
    }
    Osa2Writer::new(metadata(false), blob_fields())
        .write_all(std::fs::File::create(&with_dupes).unwrap(), &many)
        .unwrap();
    Osa2Writer::new(metadata(false), blob_fields())
        .write_all(std::fs::File::create(&without).unwrap(), &unique)
        .unwrap();

    let a = std::fs::metadata(&with_dupes).unwrap().len();
    let b = std::fs::metadata(&without).unwrap().len();
    assert_eq!(
        a, b,
        "an archive built from 5x-duplicated records must match one built from \
         the unique records alone: {a} vs {b}"
    );

    // And it still answers correctly.
    let reader = Osa2Reader::open(&with_dupes).unwrap();
    assert_eq!(
        json_of(reader.annotate_position("chr1", 1_000, "A", "G").unwrap()).unwrap(),
        "{\"id\":\"rec-0-0\"}"
    );
}
