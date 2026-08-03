//! Regression tests for long-variant (indel) handling in the .osa2 format.
//!
//! Variants whose combined ref+alt length exceeds 4 bases cannot be Var32
//! encoded and are stored in a separate `too-long.enc` block. These tests
//! verify that such variants round-trip their *values* correctly when they
//! share a chunk with short variants — and when a chunk contains only long
//! variants. gnomAD (the first v2 source) is full of indels, so this path
//! must be correct.

use fastvep_cache::annotation::{AnnotationProvider, AnnotationValue};
use fastvep_sa::fields::{Field, FieldType};
use fastvep_sa::reader_v2::Osa2Reader;
use fastvep_sa::writer_v2::{Osa2Metadata, Osa2Record, Osa2Writer};

fn metadata() -> Osa2Metadata {
    Osa2Metadata {
        format_version: 2,
        name: "test".into(),
        version: "1.0".into(),
        assembly: "GRCh38".into(),
        json_key: "test".into(),
        match_by_allele: true,
        is_array: false,
        is_positional: false,
        chunk_bits: 20,
        description: String::new(),
    }
}

fn ac_field() -> Vec<Field> {
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

fn ac_of(reader: &Osa2Reader, pos: u64, r: &str, a: &str) -> Option<i64> {
    let val = reader.annotate_position("chr1", pos, r, a).unwrap()?;
    let json = match val {
        AnnotationValue::Json(j) | AnnotationValue::Positional(j) => j,
        other => panic!("unexpected value: {:?}", other),
    };
    // JSON looks like {"allAc":<n>}
    let needle = "\"allAc\":";
    let start = json.find(needle)? + needle.len();
    let rest = &json[start..];
    let end = rest
        .find(|c: char| !c.is_ascii_digit() && c != '-')
        .unwrap_or(rest.len());
    rest[..end].parse().ok()
}

/// A chunk holding both short SNVs and long indels: every variant must return
/// its own value. Before the fix, long variants indexed the short-only,
/// var32-sorted value array with their input-order position and returned the
/// wrong AC (or a miss).
#[test]
fn long_and_short_in_same_chunk_return_own_values() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("mixed.osa2");
    let fields = ac_field();

    // All positions share chunk 0 (position >> 20 == 0). Interleave short and
    // long variants; give each a distinct AC so a mis-index is detectable.
    let records = vec![
        Osa2Record {
            chrom: "chr1".into(),
            position: 100,
            ref_allele: b"A".to_vec(),
            alt_allele: b"G".to_vec(),
            values: vec![11],
            json_blob: None,
        },
        Osa2Record {
            chrom: "chr1".into(),
            position: 200,
            ref_allele: b"ACGTACGT".to_vec(),
            alt_allele: b"A".to_vec(),
            values: vec![22],
            json_blob: None,
        },
        Osa2Record {
            chrom: "chr1".into(),
            position: 300,
            ref_allele: b"C".to_vec(),
            alt_allele: b"T".to_vec(),
            values: vec![33],
            json_blob: None,
        },
        Osa2Record {
            chrom: "chr1".into(),
            position: 400,
            ref_allele: b"G".to_vec(),
            alt_allele: b"GAAAA".to_vec(),
            values: vec![44],
            json_blob: None,
        },
    ];

    let writer = Osa2Writer::new(metadata(), fields);
    writer
        .write_all(std::fs::File::create(&path).unwrap(), &records)
        .unwrap();

    let reader = Osa2Reader::open(&path).unwrap();
    assert_eq!(ac_of(&reader, 100, "A", "G"), Some(11), "short SNV #1");
    assert_eq!(
        ac_of(&reader, 200, "ACGTACGT", "A"),
        Some(22),
        "long deletion"
    );
    assert_eq!(ac_of(&reader, 300, "C", "T"), Some(33), "short SNV #2");
    assert_eq!(
        ac_of(&reader, 400, "G", "GAAAA"),
        Some(44),
        "long insertion"
    );
}

/// A chunk containing only long variants must still be queryable. Before the
/// fix, `Chunk::is_empty()` looked only at `var32s`, so a long-only chunk was
/// treated as empty and every lookup missed.
#[test]
fn long_only_chunk_is_queryable() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("longonly.osa2");
    let fields = ac_field();

    let records = vec![
        Osa2Record {
            chrom: "chr1".into(),
            position: 500,
            ref_allele: b"ACGTACGT".to_vec(),
            alt_allele: b"A".to_vec(),
            values: vec![77],
            json_blob: None,
        },
        Osa2Record {
            chrom: "chr1".into(),
            position: 600,
            ref_allele: b"G".to_vec(),
            alt_allele: b"GTTTT".to_vec(),
            values: vec![88],
            json_blob: None,
        },
    ];

    let writer = Osa2Writer::new(metadata(), fields);
    writer
        .write_all(std::fs::File::create(&path).unwrap(), &records)
        .unwrap();

    let reader = Osa2Reader::open(&path).unwrap();
    assert_eq!(
        ac_of(&reader, 500, "ACGTACGT", "A"),
        Some(77),
        "long deletion in long-only chunk"
    );
    assert_eq!(
        ac_of(&reader, 600, "G", "GTTTT"),
        Some(88),
        "long insertion in long-only chunk"
    );
}
