//! Chunk width is a storage decision, never an answer-changing one.
//!
//! `chunk_bits` sets how much of the genome one chunk covers, and it is what
//! makes a blob-encoded dense source usable: the default 1 Mb puts so many
//! records in one chunk that a query has to decompress hundreds of megabytes of
//! JSON to read one of them, and past the limit every lookup in that megabase
//! fails outright (#101). The remedy is to convert such a source with narrower
//! chunks, which is only safe advice if narrowing them cannot change what a
//! query returns.
//!
//! It reaches further than the blob path. The width also decides the position
//! bits of every Var32 key, which is what a short-variant lookup binary-searches
//! on, so writer and reader have to agree about it through the file's metadata
//! rather than by assumption.

use fastvep_cache::annotation::{AnnotationProvider, AnnotationValue};
use fastvep_sa::fields::{Field, FieldType};
use fastvep_sa::reader_v2::Osa2Reader;
use fastvep_sa::writer_v2::{Osa2Metadata, Osa2Record, Osa2Writer};

fn metadata(chunk_bits: u32) -> Osa2Metadata {
    Osa2Metadata {
        format_version: 2,
        name: "Widths".into(),
        version: "1".into(),
        assembly: "GRCh38".into(),
        json_key: "widths".into(),
        match_by_allele: true,
        is_array: false,
        is_positional: false,
        chunk_bits,
        description: "test".into(),
    }
}

/// One whole-record JSON blob, the encoding `sa-convert` produces.
fn blob_field() -> Vec<Field> {
    vec![Field {
        field: String::new(),
        alias: String::new(),
        ftype: FieldType::JsonBlob,
        multiplier: 1,
        zigzag: false,
        missing_value: u32::MAX,
        missing_string: ".".into(),
        description: String::new(),
    }]
}

/// Records spread over several megabases and deliberately sitting on chunk
/// boundaries, with the allele shapes that take each of the reader's three
/// lookup paths: 2-bit packable short variants, long ones, and non-ACGT ones.
fn records() -> Vec<Osa2Record> {
    let mut out = Vec::new();
    let mut push = |chrom: &str, pos: u32, r: &[u8], a: &[u8], tag: &str| {
        out.push(Osa2Record {
            chrom: chrom.into(),
            position: pos,
            ref_allele: r.to_vec(),
            alt_allele: a.to_vec(),
            values: Vec::new(),
            json_blob: Some(format!("{{\"tag\":\"{tag}\",\"pos\":{pos}}}")),
        });
    };
    // Every power-of-two boundary between 2^14 and 2^20 is a chunk edge for one
    // of the widths under test, so include the bases either side of each.
    let mut positions: Vec<u32> = Vec::new();
    for bits in 14..=20u32 {
        let w = 1u32 << bits;
        for mult in 1..=3u32 {
            let edge = w * mult;
            positions.extend([edge.saturating_sub(1), edge, edge + 1]);
        }
    }
    positions.extend([1, 2, 1_000, 500_000, 3_000_000, 4_194_303]);
    positions.sort_unstable();
    positions.dedup();

    for (i, &p) in positions.iter().enumerate() {
        for chrom in ["chr1", "chr2"] {
            match i % 3 {
                0 => push(chrom, p, b"A", b"G", "short"),
                1 => push(chrom, p, b"ACGTACGT", b"T", "long"),
                _ => push(chrom, p, b"N", b"A", "nonacgt"),
            }
        }
    }
    out.sort_by(|a, b| a.chrom.cmp(&b.chrom).then(a.position.cmp(&b.position)));
    out
}

#[test]
fn every_chunk_width_answers_identically() {
    let dir = tempfile::tempdir().unwrap();
    let recs = records();
    assert!(recs.len() > 100, "need a meaningful record set");

    // Ask about every stored record, plus positions that are absent, so a miss
    // has to stay a miss at every width too.
    let mut queries: Vec<(String, u64, Vec<u8>, Vec<u8>)> = recs
        .iter()
        .map(|r| {
            (
                r.chrom.clone(),
                r.position as u64,
                r.ref_allele.clone(),
                r.alt_allele.clone(),
            )
        })
        .collect();
    for r in &recs {
        queries.push((
            r.chrom.clone(),
            r.position as u64 + 7,
            b"A".to_vec(),
            b"C".to_vec(),
        ));
        queries.push((
            r.chrom.clone(),
            r.position as u64,
            b"A".to_vec(),
            b"TTTT".to_vec(),
        ));
    }

    let mut baseline: Option<Vec<Option<String>>> = None;
    for chunk_bits in [20u32, 18, 16, 14, 12, 8, 1] {
        let path = dir.path().join(format!("w{chunk_bits}.osa2"));
        Osa2Writer::new(metadata(chunk_bits), blob_field())
            .write_all(std::fs::File::create(&path).unwrap(), &recs)
            .unwrap();
        let reader = Osa2Reader::open(&path).unwrap();

        let answers: Vec<Option<String>> = queries
            .iter()
            .map(|(c, p, r, a)| {
                let got = reader
                    .annotate_position(
                        c,
                        *p,
                        std::str::from_utf8(r).unwrap(),
                        std::str::from_utf8(a).unwrap(),
                    )
                    .unwrap_or_else(|e| panic!("chunk_bits {chunk_bits}: lookup failed: {e}"));
                got.map(|v| match v {
                    AnnotationValue::Json(j) | AnnotationValue::Positional(j) => j,
                    AnnotationValue::Interval(v) => v.join(","),
                })
            })
            .collect();

        // The record set must actually be found, or the comparison is vacuous.
        let hits = answers.iter().filter(|a| a.is_some()).count();
        assert!(
            hits >= recs.len(),
            "chunk_bits {chunk_bits}: only {hits} hits for {} records",
            recs.len()
        );

        match &baseline {
            None => baseline = Some(answers),
            Some(want) => {
                assert_eq!(
                    &answers, want,
                    "chunk_bits {chunk_bits} answered differently from the 1 Mb layout"
                );
            }
        }
    }
}
