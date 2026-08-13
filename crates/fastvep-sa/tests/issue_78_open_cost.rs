//! Issue #78: opening an `.osa2` must not scale with the number of ZIP entries
//! in random I/O.
//!
//! The reported symptom was a 25-30 minute silent stall before annotation
//! started after adding per-chromosome gnomAD shards to `--sa-dir`. The cause
//! was `open` resolving every entry's data offset from its *local* file header,
//! which is one random read per entry scattered across a multi-GB file. These
//! tests pin the two halves of the contract: `open` touches a constant number
//! of local headers regardless of archive size, and the deferred resolution
//! still yields byte-identical lookups.

use fastvep_cache::annotation::{AnnotationProvider, AnnotationValue};
use fastvep_sa::fields::{Field, FieldType};
use fastvep_sa::reader_v2::Osa2Reader;
use fastvep_sa::writer_v2::{Osa2Metadata, Osa2Record, Osa2StreamWriter};
use std::path::{Path, PathBuf};
use tempfile::TempDir;

/// Multiple value columns so each chunk becomes several ZIP entries, the way a
/// real gnomAD shard does.
fn fields() -> Vec<Field> {
    (0..4)
        .map(|i| Field {
            field: format!("F{i}"),
            alias: format!("f{i}"),
            ftype: FieldType::Integer,
            multiplier: 1,
            zigzag: false,
            missing_value: u32::MAX,
            missing_string: ".".into(),
            description: format!("field {i}"),
        })
        .collect()
}

fn metadata(chunk_bits: u32) -> Osa2Metadata {
    Osa2Metadata {
        format_version: 2,
        name: "issue78".into(),
        version: "test".into(),
        assembly: "GRCh38".into(),
        json_key: "issue78".into(),
        match_by_allele: true,
        is_array: false,
        is_positional: false,
        chunk_bits,
        description: String::new(),
    }
}

/// One record every `1 << chunk_bits` bases, so `n_chunks` chunks are emitted
/// and the entry count grows linearly with `n_chunks`.
fn build(dir: &Path, name: &str, n_chunks: u32, chunk_bits: u32) -> PathBuf {
    let flds = fields();
    let path = dir.join(name);
    let file = std::io::BufWriter::new(std::fs::File::create(&path).unwrap());
    let mut writer = Osa2StreamWriter::new(file, &metadata(chunk_bits), flds.clone()).unwrap();
    for c in 0..n_chunks {
        for k in 0..4u32 {
            let position = (c << chunk_bits) + 10 + k;
            writer
                .push(Osa2Record {
                    chrom: "chr1".into(),
                    position,
                    ref_allele: b"A".to_vec(),
                    alt_allele: b"G".to_vec(),
                    values: (0..flds.len()).map(|f| position + f as u32).collect(),
                    json_blob: None,
                })
                .unwrap();
        }
    }
    writer.finish().unwrap();
    path
}

/// The core regression: a 16x bigger archive must cost the *same* number of
/// local-header reads at open. Before the fix this was one per entry.
#[test]
fn open_cost_is_independent_of_entry_count() {
    let dir = TempDir::new().unwrap();
    let small = Osa2Reader::open(&build(dir.path(), "small.osa2", 8, 12)).unwrap();
    let big = Osa2Reader::open(&build(dir.path(), "big.osa2", 128, 12)).unwrap();

    assert!(
        big.entry_count() > small.entry_count() * 8,
        "fixture should scale: {} vs {}",
        small.entry_count(),
        big.entry_count()
    );
    assert_eq!(
        small.header_read_count(),
        big.header_read_count(),
        "open must read a constant number of local headers, not one per entry \
         (small={} entries, big={} entries)",
        small.entry_count(),
        big.entry_count()
    );
    // Only the metadata + config prelude is read at open.
    assert!(
        big.header_read_count() <= 4,
        "open read {} local headers; expected only the prelude entries",
        big.header_read_count()
    );
}

/// Deferring the offset resolution must not change what a lookup returns, and
/// a repeated lookup must not re-parse the header it already memoized.
#[test]
fn deferred_offsets_still_resolve_correctly() {
    let dir = TempDir::new().unwrap();
    let reader = Osa2Reader::open(&build(dir.path(), "q.osa2", 64, 12)).unwrap();

    let mut seen = 0;
    for c in 0..64u32 {
        for k in 0..4u32 {
            let position = ((c << 12) + 10 + k) as u64;
            let value = reader
                .annotate_position("chr1", position, "A", "G")
                .unwrap()
                .unwrap_or_else(|| panic!("missing annotation at chr1:{position}"));
            let AnnotationValue::Json(json) = value else {
                panic!("expected an allele-matched JSON value");
            };
            for f in 0..4u32 {
                assert!(
                    json.contains(&format!("\"f{f}\":{}", position as u32 + f)),
                    "value column f{f} wrong at chr1:{position}: {json}"
                );
            }
            seen += 1;
        }
    }
    assert_eq!(seen, 256);

    // A miss must stay a miss, not become an error, once offsets are lazy.
    assert!(reader
        .annotate_position("chr1", 999_999_999, "A", "G")
        .unwrap()
        .is_none());

    let after_first_pass = reader.header_read_count();
    for c in 0..64u32 {
        reader
            .annotate_position("chr1", ((c << 12) + 10) as u64, "A", "G")
            .unwrap();
    }
    assert_eq!(
        after_first_pass,
        reader.header_read_count(),
        "a second pass must reuse the memoized data offsets"
    );
}

/// Concurrent first-touch of the same entries must not deadlock, double-count
/// incorrectly, or produce wrong values - the `OnceLock` fill is racy by design.
#[test]
fn concurrent_first_touch_is_consistent() {
    use rayon::prelude::*;

    let dir = TempDir::new().unwrap();
    let reader = Osa2Reader::open(&build(dir.path(), "par.osa2", 64, 12)).unwrap();

    let hits: usize = (0..64u32)
        .into_par_iter()
        .flat_map(|c| {
            (0..4u32)
                .into_par_iter()
                .map(move |k| ((c << 12) + 10 + k) as u64)
        })
        .map(|position| {
            reader
                .annotate_position("chr1", position, "A", "G")
                .unwrap()
                .is_some() as usize
        })
        .sum();
    assert_eq!(hits, 256);
}
