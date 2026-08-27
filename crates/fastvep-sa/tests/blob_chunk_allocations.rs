//! Pins the cost of loading a JSON-blob chunk to the chunk, not to its rows.
//!
//! A chunk holds every record in its genomic width, and reading one of them
//! used to build a `String` per row - one heap allocation and one copy each,
//! paid on the way to returning a single record. For the blob-encoded form of a
//! source that scores nearly every base that is hundreds of thousands of
//! allocations per lookup (#101).
//!
//! This file installs a counting global allocator, so it deliberately holds
//! exactly ONE test: `cargo test` runs the tests in a binary concurrently, and a
//! second test allocating on another thread would be counted here.

use fastvep_cache::annotation::AnnotationProvider;
use fastvep_sa::fields::{Field, FieldType};
use fastvep_sa::reader_v2::Osa2Reader;
use fastvep_sa::writer_v2::{Osa2Metadata, Osa2Record, Osa2Writer};
use std::alloc::{GlobalAlloc, Layout, System};
use std::sync::atomic::{AtomicUsize, Ordering};

static ALLOCS: AtomicUsize = AtomicUsize::new(0);

struct Counting;

unsafe impl GlobalAlloc for Counting {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        ALLOCS.fetch_add(1, Ordering::Relaxed);
        unsafe { System.alloc(layout) }
    }
    unsafe fn dealloc(&self, ptr: *mut u8, layout: Layout) {
        unsafe { System.dealloc(ptr, layout) }
    }
    unsafe fn realloc(&self, ptr: *mut u8, layout: Layout, new_size: usize) -> *mut u8 {
        ALLOCS.fetch_add(1, Ordering::Relaxed);
        unsafe { System.realloc(ptr, layout, new_size) }
    }
}

#[global_allocator]
static ALLOCATOR: Counting = Counting;

fn allocations_during<F: FnOnce()>(f: F) -> usize {
    let before = ALLOCS.load(Ordering::Relaxed);
    f();
    ALLOCS.load(Ordering::Relaxed) - before
}

/// A whole-record blob field: empty alias, so `reconstruct_json` returns the
/// stored object verbatim. This is what `sa-convert` writes.
fn whole_record_blob_field() -> Vec<Field> {
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

const ROWS: u32 = 20_000;

#[test]
fn loading_a_blob_chunk_does_not_allocate_per_row() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("blobs.osa2");

    // Every record lands in chunk 0, so one lookup loads all ROWS of them.
    let records: Vec<Osa2Record> = (0..ROWS)
        .map(|i| Osa2Record {
            chrom: "chr1".into(),
            position: 1_000 + i,
            ref_allele: b"A".to_vec(),
            alt_allele: b"G".to_vec(),
            values: Vec::new(),
            json_blob: Some(format!("{{\"score\":{i},\"symbol\":\"GENE{i}\"}}")),
        })
        .collect();
    Osa2Writer::new(
        Osa2Metadata {
            format_version: 2,
            name: "Blobs".into(),
            version: "1".into(),
            assembly: "GRCh38".into(),
            json_key: "blobs".into(),
            match_by_allele: true,
            is_array: false,
            is_positional: false,
            chunk_bits: 20,
            description: "test".into(),
        },
        whole_record_blob_field(),
    )
    .write_all(std::fs::File::create(&path).unwrap(), &records)
    .unwrap();

    let reader = Osa2Reader::open(&path).unwrap();

    // The very first lookup is the one that builds the chunk, so measure it
    // before anything has warmed the cache.
    let mut answer = None;
    let allocs = allocations_during(|| {
        answer = reader.annotate_position("chr1", 1_777, "A", "G").unwrap();
    });
    assert!(answer.is_some(), "the record should have been found");

    // The chunk is built from a fixed number of buffers - the ZIP entries, the
    // decompressed blob, and its offsets - so the count must not scale with
    // ROWS. Measured at 77; the budget is loose because the point of the
    // assertion is the shape of the growth, not the constant. A `String` per row
    // would put this above ROWS.
    assert!(
        allocs < ROWS as usize / 10,
        "building a {ROWS}-row blob chunk took {allocs} allocations; \
         that scales with the row count rather than the chunk"
    );

    // And every row still reads back byte-for-byte what was stored - indexing
    // the buffer must delimit the rows exactly as splitting it did.
    for probe in [0u32, 1, ROWS / 2, ROWS - 1] {
        let got = reader
            .annotate_position("chr1", (1_000 + probe) as u64, "A", "G")
            .unwrap()
            .expect("row should be present");
        let json = match got {
            fastvep_cache::annotation::AnnotationValue::Json(j) => j,
            other => panic!("expected a JSON value for row {probe}, got {other:?}"),
        };
        assert_eq!(
            json,
            format!("{{\"score\":{probe},\"symbol\":\"GENE{probe}\"}}"),
            "row {probe} did not round-trip"
        );
    }
}
