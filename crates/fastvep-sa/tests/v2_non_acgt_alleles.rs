//! Regression tests: `.osa2` must not drop variants whose alleles contain
//! bytes outside `ACGT`.
//!
//! Found by running the full real ClinVar release (4,438,232 records) through
//! both formats and querying every v1 record against the v2 database: **668
//! records were retrievable from the `.osa` and missing from the `.osa2`.**
//! Their alleles all carry `N` runs (ClinVar's `Microsatellite` and long
//! `Insertion` entries), and both v2 key encodings — Var32 and kmer16 — are two
//! bits per base, so neither can represent them. The writer dropped them with a
//! `log::warn!` that the CLI installs no logger to display, making it silent
//! data loss on a source that already built to v2 by default.
//!
//! Such variants now go in a third per-chunk bucket that stores the allele bytes
//! verbatim and is binary-searched on `(position, ref, alt)`.

use fastvep_cache::annotation::{AnnotationProvider, AnnotationValue};
use fastvep_sa::fields::{Field, FieldType};
use fastvep_sa::reader_v2::Osa2Reader;
use fastvep_sa::writer_v2::{Osa2Metadata, Osa2Record, Osa2StreamWriter, Osa2Writer};
use std::io::BufWriter;

fn metadata() -> Osa2Metadata {
    Osa2Metadata {
        format_version: 2,
        name: "Test".into(),
        version: "1".into(),
        assembly: "GRCh38".into(),
        json_key: "test".into(),
        match_by_allele: true,
        is_array: false,
        is_positional: false,
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

/// The allele shapes that used to be dropped, alongside ordinary ones that must
/// keep working. Ordering here is deliberately not sorted, and positions are
/// spread so short/long/non-ACGT variants share chunks.
fn mixed_records() -> Vec<(u32, &'static str, &'static str, &'static str)> {
    vec![
        // (position, ref, alt, payload)
        (100, "A", "G", "snv-short"),
        // Short but non-ACGT: fails Var32 encoding.
        (100, "A", "N", "short-alt-N"),
        (101, "N", "A", "short-ref-N"),
        // IUPAC ambiguity codes, not just N.
        (102, "A", "R", "iupac-R"),
        (103, "Y", "A", "iupac-Y"),
        // Long (ref+alt > 4) and non-ACGT: fails kmer16 encoding. This is the
        // exact shape of the real ClinVar Microsatellite entries.
        (
            200,
            "C",
            "CAAAAGAACCATTTNNNNNNNNNNAAAAAAAAAA",
            "long-alt-with-N-run",
        ),
        (201, "GATTACANNNN", "G", "long-ref-with-N-run"),
        // Ordinary long variant in the same neighbourhood.
        (202, "GATTACA", "G", "long-acgt"),
        // Lowercase must still be treated as ACGT, not shunted to the bucket.
        (300, "a", "g", "lowercase-snv"),
        (301, "gattaca", "g", "lowercase-long"),
        // A far chunk holding ONLY a non-ACGT variant.
        (5_000_000, "A", "N", "lonely-N-in-its-own-chunk"),
    ]
}

fn blob_records() -> Vec<Osa2Record> {
    mixed_records()
        .into_iter()
        .map(|(pos, r, a, payload)| Osa2Record {
            chrom: "chr1".into(),
            position: pos,
            ref_allele: r.as_bytes().to_vec(),
            alt_allele: a.as_bytes().to_vec(),
            values: Vec::new(),
            json_blob: Some(format!("{{\"id\":\"{payload}\"}}")),
        })
        .collect()
}

/// Assert every record written is retrievable with its own alleles, carrying its
/// own payload — i.e. no drops and no cross-wired value slots.
fn assert_all_retrievable(reader: &Osa2Reader) {
    for (pos, r, a, payload) in mixed_records() {
        let json = json_of(reader.annotate_position("chr1", pos as u64, r, a).unwrap())
            .unwrap_or_else(|| panic!("chr1:{pos} {r}>{a} ({payload}) was dropped"));
        assert!(
            json.contains(payload),
            "chr1:{pos} {r}>{a} returned the wrong record's payload: {json}"
        );
    }
}

#[test]
fn buffered_writer_keeps_non_acgt_alleles() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("t.osa2");
    let mut records = blob_records();
    records.sort_by_key(|r| r.position);

    Osa2Writer::new(metadata(), blob_fields())
        .write_all(std::fs::File::create(&path).unwrap(), &records)
        .unwrap();

    assert_all_retrievable(&Osa2Reader::open(&path).unwrap());
}

#[test]
fn streaming_writer_keeps_non_acgt_alleles() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("t.osa2");
    let mut records = blob_records();
    records.sort_by_key(|r| r.position);

    let f = std::fs::File::create(&path).unwrap();
    let mut w = Osa2StreamWriter::new(BufWriter::new(f), &metadata(), blob_fields()).unwrap();
    for r in records {
        w.push(r).unwrap();
    }
    w.finish().unwrap();

    assert_all_retrievable(&Osa2Reader::open(&path).unwrap());
}

#[test]
fn non_acgt_alleles_resolve_to_their_own_numeric_columns() {
    // The bucket's `idx` must point at the right slot in the parallel value
    // arrays, which is exactly what the original long-variant `idx` bug got
    // wrong. Give each record a distinct value and check every one.
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("t.osa2");

    let specs = mixed_records();
    let mut records: Vec<Osa2Record> = specs
        .iter()
        .enumerate()
        .map(|(i, (pos, r, a, _))| Osa2Record {
            chrom: "chr1".into(),
            position: *pos,
            ref_allele: r.as_bytes().to_vec(),
            alt_allele: a.as_bytes().to_vec(),
            values: vec![i as u32 + 1],
            json_blob: None,
        })
        .collect();
    records.sort_by_key(|r| r.position);

    Osa2Writer::new(metadata(), one_int_field())
        .write_all(std::fs::File::create(&path).unwrap(), &records)
        .unwrap();

    let reader = Osa2Reader::open(&path).unwrap();
    for (i, (pos, r, a, label)) in specs.iter().enumerate() {
        let json = json_of(reader.annotate_position("chr1", *pos as u64, r, a).unwrap())
            .unwrap_or_else(|| panic!("chr1:{pos} {r}>{a} ({label}) was dropped"));
        assert_eq!(
            json,
            format!("{{\"ac\":{}}}", i + 1),
            "chr1:{pos} {r}>{a} ({label}) resolved to the wrong value slot"
        );
    }
}

#[test]
fn a_wrong_non_acgt_allele_still_misses() {
    // The bucket must match exactly, not merely by position — otherwise
    // recovering these records would trade data loss for false positives.
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("t.osa2");
    let mut records = blob_records();
    records.sort_by_key(|r| r.position);
    Osa2Writer::new(metadata(), blob_fields())
        .write_all(std::fs::File::create(&path).unwrap(), &records)
        .unwrap();
    let reader = Osa2Reader::open(&path).unwrap();

    // Right position, different non-ACGT allele.
    assert!(json_of(reader.annotate_position("chr1", 100, "A", "NN").unwrap()).is_none());
    assert!(json_of(reader.annotate_position("chr1", 100, "N", "N").unwrap()).is_none());
    // Right allele shape, wrong position.
    assert!(json_of(reader.annotate_position("chr1", 999, "A", "N").unwrap()).is_none());
    // A non-ACGT query against a position holding only ACGT records.
    assert!(json_of(
        reader
            .annotate_position("chr1", 202, "GATTACN", "G")
            .unwrap()
    )
    .is_none());
    // And the ACGT record at a position that also holds a non-ACGT one is
    // unaffected by the bucket's presence.
    assert!(json_of(reader.annotate_position("chr1", 100, "A", "G").unwrap()).is_some());
}

#[test]
fn archives_without_non_acgt_alleles_carry_no_extra_entry() {
    // The `other.enc` sub-entry is only emitted when a chunk actually needs it,
    // so ordinary sources produce the same archive they did before the bucket
    // existed — and a reader that predates it never encounters the entry.
    let dir = tempfile::tempdir().unwrap();
    let clean = dir.path().join("clean.osa2");
    let dirty = dir.path().join("dirty.osa2");

    let acgt_only: Vec<Osa2Record> = [(100u32, "A", "G"), (200, "GATTACA", "G")]
        .iter()
        .map(|(pos, r, a)| Osa2Record {
            chrom: "chr1".into(),
            position: *pos,
            ref_allele: r.as_bytes().to_vec(),
            alt_allele: a.as_bytes().to_vec(),
            values: Vec::new(),
            json_blob: Some("{\"id\":\"x\"}".into()),
        })
        .collect();
    Osa2Writer::new(metadata(), blob_fields())
        .write_all(std::fs::File::create(&clean).unwrap(), &acgt_only)
        .unwrap();

    let mut with_n = acgt_only.clone();
    with_n.push(Osa2Record {
        chrom: "chr1".into(),
        position: 300,
        ref_allele: b"A".to_vec(),
        alt_allele: b"N".to_vec(),
        values: Vec::new(),
        json_blob: Some("{\"id\":\"n\"}".into()),
    });
    Osa2Writer::new(metadata(), blob_fields())
        .write_all(std::fs::File::create(&dirty).unwrap(), &with_n)
        .unwrap();

    let entry_names = |p: &std::path::Path| -> Vec<String> {
        let f = std::fs::File::open(p).unwrap();
        let mut zip = zip::ZipArchive::new(f).unwrap();
        (0..zip.len())
            .map(|i| zip.by_index(i).unwrap().name().to_string())
            .collect()
    };

    let clean_entries = entry_names(&clean);
    assert!(
        !clean_entries.iter().any(|n| n.ends_with("other.enc")),
        "an ACGT-only archive must not carry other.enc: {clean_entries:?}"
    );
    let dirty_entries = entry_names(&dirty);
    assert!(
        dirty_entries.iter().any(|n| n.ends_with("other.enc")),
        "an archive with an N allele must carry other.enc: {dirty_entries:?}"
    );
}
