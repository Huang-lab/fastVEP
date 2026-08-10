//! Pins the allocation behaviour of the SA query path.
//!
//! Chromosome resolution used to call `chrom_aliases`, which builds a
//! `Vec<String>` (two or three heap allocations) on every query. With N
//! providers in a `--sa-dir` and millions of variants that is tens of millions
//! of allocations to resolve a name that is constant per (reader, chromosome).
//! Both readers now invert the alias rules once at open into a map, so a query
//! resolves its chromosome with one hash lookup and no allocation.
//!
//! This file installs a counting global allocator, so it deliberately holds
//! exactly ONE test: `cargo test` runs the tests in a binary concurrently, and a
//! second test allocating on another thread would be counted here.

use fastvep_cache::annotation::AnnotationProvider;
use fastvep_sa::common::AnnotationRecord;
use fastvep_sa::index::IndexHeader;
use fastvep_sa::fields::{Field, FieldType};
use fastvep_sa::interval::{IntervalHeader, IntervalIndex};
use fastvep_sa::reader::SaReader;
use fastvep_sa::reader_v2::Osa2Reader;
use fastvep_sa::writer::SaWriter;
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

/// Allocations performed while running `f`.
fn allocations_during<F: FnOnce()>(f: F) -> usize {
    let before = ALLOCS.load(Ordering::Relaxed);
    f();
    ALLOCS.load(Ordering::Relaxed) - before
}

fn v2_metadata() -> Osa2Metadata {
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

fn v1_header() -> IndexHeader {
    IndexHeader {
        schema_version: fastvep_sa::common::SCHEMA_VERSION,
        json_key: "test".into(),
        name: "Test".into(),
        version: "1".into(),
        description: "test".into(),
        assembly: "GRCh38".into(),
        match_by_allele: true,
        is_array: false,
        is_positional: false,
    }
}

#[test]
fn resolving_a_chromosome_allocates_nothing_on_any_reader() {
    let dir = tempfile::tempdir().unwrap();

    // ---- v2 `.osa2` ----
    let v2_path = dir.path().join("t.osa2");
    let records: Vec<Osa2Record> = (0..64)
        .map(|i| Osa2Record {
            chrom: "chr1".into(),
            position: 10_000 + i * 10,
            ref_allele: b"A".to_vec(),
            alt_allele: b"G".to_vec(),
            values: vec![i + 1],
            json_blob: None,
        })
        .collect();
    Osa2Writer::new(v2_metadata(), one_int_field())
        .write_all(std::fs::File::create(&v2_path).unwrap(), &records)
        .unwrap();
    let v2 = Osa2Reader::open(&v2_path).unwrap();

    // ---- v1 `.osa` ----
    let v1_base = dir.path().join("t");
    let v1_records: Vec<AnnotationRecord> = (0..64u32)
        .map(|i| AnnotationRecord {
            chrom_idx: 0,
            position: 10_000 + i * 10,
            ref_allele: "A".into(),
            alt_allele: "G".into(),
            json: "{\"ac\":1}".into(),
        })
        .collect();
    let mut writer = SaWriter::new(v1_header());
    writer
        .write_to_files(&v1_base, v1_records.into_iter(), &["chr1".to_string()])
        .unwrap();
    let v1 = SaReader::open(&v1_base.with_extension("osa")).unwrap();

    // ---- `.osi` ----
    let mut osi = IntervalIndex::new(IntervalHeader {
        schema_version: fastvep_sa::common::SCHEMA_VERSION,
        json_key: "iv".into(),
        name: "IV".into(),
        version: "1".into(),
        assembly: "GRCh38".into(),
    });
    for i in 0..64u32 {
        osi.add(fastvep_sa::common::IntervalRecord {
            chrom: "chr1".into(),
            start: 10_000 + i * 10,
            end: 10_005 + i * 10,
            json: "{}".into(),
        });
    }
    osi.sort();
    let osi_path = dir.path().join("t.osi");
    osi.write_to(&mut std::fs::File::create(&osi_path).unwrap())
        .unwrap();
    let osi = fastvep_sa::interval::OsiReader::open(&osi_path).unwrap();

    let readers: [(&str, &dyn AnnotationProvider); 3] =
        [("osa2", &v2), ("osa", &v1), ("osi", &osi)];

    // Warm every cache and every lazily-resolved offset first, so the measured
    // window covers only steady-state lookups.
    for (_, r) in &readers {
        for _ in 0..3 {
            r.annotate_position("chr1", 10_000, "A", "G").unwrap();
            r.annotate_position("chr2", 10_000, "A", "G").unwrap();
        }
    }

    for (name, reader) in &readers {
        // A query for a contig the database does not carry resolves to a miss
        // and returns. That is pure chromosome resolution with nothing else in
        // the way, so it must allocate exactly nothing.
        let allocs = allocations_during(|| {
            for _ in 0..1_000 {
                assert!(reader
                    .annotate_position("chr2", 10_000, "A", "G")
                    .unwrap()
                    .is_none());
            }
        });
        assert_eq!(
            allocs, 0,
            "{name}: 1000 absent-contig queries allocated {allocs} times; \
             chromosome resolution must be allocation-free"
        );

        // The same, via an alias of an absent contig — the path that used to
        // build the alias Vec before discovering the miss.
        let allocs = allocations_during(|| {
            for _ in 0..1_000 {
                assert!(reader.annotate_position("2", 10_000, "A", "G").unwrap().is_none());
            }
        });
        assert_eq!(allocs, 0, "{name}: aliased absent-contig queries allocated {allocs} times");

        // A hit necessarily allocates (it returns an owned JSON String), but
        // resolving via an alias must cost the same as the canonical spelling.
        let canonical = allocations_during(|| {
            for _ in 0..1_000 {
                assert!(reader
                    .annotate_position("chr1", 10_000, "A", "G")
                    .unwrap()
                    .is_some());
            }
        });
        let aliased = allocations_during(|| {
            for _ in 0..1_000 {
                assert!(reader.annotate_position("1", 10_000, "A", "G").unwrap().is_some());
            }
        });
        assert_eq!(
            canonical, aliased,
            "{name}: alias resolution cost {aliased} allocations vs {canonical} for the \
             canonical name; resolving a spelling must be free"
        );
    }
}
