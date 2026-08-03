//! Round-trip and contract tests for `Osa2StreamWriter`, the bounded-memory
//! streaming builder used for genome-scale `.osa2` sources.

use fastvep_cache::annotation::{AnnotationProvider, AnnotationValue};
use fastvep_sa::fields::{Field, FieldType};
use fastvep_sa::reader_v2::Osa2Reader;
use fastvep_sa::writer_v2::{Osa2Metadata, Osa2Record, Osa2StreamWriter};

fn metadata() -> Osa2Metadata {
    Osa2Metadata {
        format_version: 2,
        name: "stream".into(),
        version: "1.0".into(),
        assembly: "GRCh38".into(),
        json_key: "stream".into(),
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
    let needle = "\"allAc\":";
    let start = json.find(needle)? + needle.len();
    let rest = &json[start..];
    let end = rest
        .find(|c: char| !c.is_ascii_digit() && c != '-')
        .unwrap_or(rest.len());
    rest[..end].parse().ok()
}

fn rec(pos: u32, r: &[u8], a: &[u8], ac: u32) -> Osa2Record {
    Osa2Record {
        chrom: "chr1".into(),
        position: pos,
        ref_allele: r.to_vec(),
        alt_allele: a.to_vec(),
        values: vec![ac],
        json_blob: None,
    }
}

/// Streamed records spanning several 1 MB chunks (plus one indel) must round
/// trip exactly, matching what the buffered writer would produce.
#[test]
fn stream_writer_round_trips_across_chunks() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("stream.osa2");

    // chunk_bits = 20 => chunk size ~1,048,576. Place records in chunks 0, 1, 2.
    let records = vec![
        rec(100, b"A", b"G", 1),
        rec(500, b"ACGTACGT", b"A", 2), // long variant in chunk 0
        rec(1_100_000, b"C", b"T", 3),  // chunk 1
        rec(1_100_500, b"G", b"C", 4),  // chunk 1
        rec(2_200_000, b"T", b"A", 5),  // chunk 2
    ];

    let file = std::fs::File::create(&path).unwrap();
    let mut w =
        Osa2StreamWriter::new(std::io::BufWriter::new(file), &metadata(), ac_field()).unwrap();
    for r in &records {
        w.push(r.clone()).unwrap();
    }
    w.finish().unwrap();

    let reader = Osa2Reader::open(&path).unwrap();
    assert_eq!(ac_of(&reader, 100, "A", "G"), Some(1));
    assert_eq!(
        ac_of(&reader, 500, "ACGTACGT", "A"),
        Some(2),
        "long variant"
    );
    assert_eq!(ac_of(&reader, 1_100_000, "C", "T"), Some(3));
    assert_eq!(ac_of(&reader, 1_100_500, "G", "C"), Some(4));
    assert_eq!(ac_of(&reader, 2_200_000, "T", "A"), Some(5));
    // A definite miss.
    assert_eq!(ac_of(&reader, 999, "A", "G"), None);
}

/// Reopening an already-flushed (chrom, chunk) — i.e. unsorted input — must be
/// rejected rather than silently emitting a duplicate ZIP entry.
#[test]
fn stream_writer_rejects_unsorted_input() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("bad.osa2");
    let file = std::fs::File::create(&path).unwrap();
    let mut w = Osa2StreamWriter::new(file, &metadata(), ac_field()).unwrap();

    w.push(rec(100, b"A", b"G", 1)).unwrap();
    w.push(rec(2_200_000, b"C", b"T", 2)).unwrap(); // chunk 2
                                                    // Back to chunk 0 — out of order; chunk 0 was already flushed.
    let err = w.push(rec(200, b"A", b"C", 3)).unwrap_err();
    assert!(
        err.to_string().contains("not sorted"),
        "expected unsorted-input error, got: {err}"
    );
}
