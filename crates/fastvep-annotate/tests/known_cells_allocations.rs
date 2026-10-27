//! Pins what working out VEP's known-variant columns costs, per variant.
//!
//! [`known_cells`] runs once for each variant that the variation cache attached a
//! known variant to, in the parallel phase of `run_annotate`, and its answer is
//! shared by every transcript row of the variant. Three things are held:
//!
//! - a variant with no known variant, which is most of a callset, costs nothing;
//! - a run that did not ask for the columns costs nothing either, whatever the
//!   variant has;
//! - one that did pays for the text of the cells it writes and a few scratch
//!   lists, and the count does not grow with the number of populations a record
//!   lists beyond one small string for each.
//!
//! This file installs a counting global allocator, so it deliberately holds
//! exactly ONE test: `cargo test` runs the tests of a binary concurrently, and a
//! second one allocating on another thread would be counted here.

use fastvep_annotate::known_cells;
use fastvep_cache::providers::{KnownRecord, KnownSite};
use fastvep_cache::variation::{AlleleMatch, Matched, VariationCacheRecord};
use fastvep_core::Allele;
use fastvep_io::known::KnownFlags;
use std::alloc::{GlobalAlloc, Layout, System};
use std::collections::HashMap;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;

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

/// A dbSNP-like record: `AF`, five 1000 Genomes populations, gnomAD exomes and
/// genomes with several sub-populations, a clinical significance and a phenotype.
fn site() -> KnownSite {
    let mut frequencies = HashMap::new();
    for (key, value) in [
        ("AF", "A:0.0123"),
        ("AFR", "A:0.0345"),
        ("AMR", "A:0.0067"),
        ("EAS", "A:0"),
        ("EUR", "A:0.0102"),
        ("SAS", "A:0.0003"),
        ("gnomADe", "A:8.291e-06"),
        ("gnomADe_AFR", "A:2e-05"),
        ("gnomADe_NFE", "A:1.1e-05"),
        ("gnomADg", "A:0.0004"),
        ("gnomADg_AFR", "A:0.0012"),
        ("gnomADg_NFE", "A:0.0002"),
    ] {
        frequencies.insert(key.to_string(), value.to_string());
    }
    let record = VariationCacheRecord {
        variation_name: "rs1".into(),
        failed: false,
        somatic: false,
        start: 100,
        end: 100,
        allele_string: "C/A".into(),
        strand: 1,
        minor_allele: None,
        minor_allele_freq: None,
        clin_sig: Some("likely_benign,benign".into()),
        clin_sig_allele: Some("A:likely_benign".into()),
        clin_sig_ref_allele: Some("C".into()),
        phenotype_or_disease: true,
        pubmed: vec!["123456".into(), "234567".into()],
        frequencies,
    };
    KnownSite {
        records: vec![KnownRecord {
            record: Arc::new(record),
            matched: Matched::Alleles(vec![AlleleMatch {
                alt_index: 0,
                record_allele: "A".into(),
            }]),
        }],
    }
}

#[test]
fn working_out_the_known_cells_costs_nothing_unless_there_is_something_to_write() {
    let alleles = vec![Allele::Sequence(b"A".to_vec())];
    let site = site();
    let everything = KnownFlags::everything();

    // warm anything lazily initialised
    let _ = known_cells(&site, "C", &alleles, everything);

    let mut cells = Vec::new();
    let nothing_known = allocations_during(|| {
        cells = known_cells(&KnownSite::default(), "C", &alleles, everything);
    });
    assert!(cells.is_empty());
    assert_eq!(
        nothing_known, 0,
        "a variant with no known variant must cost nothing"
    );

    let not_asked_for = allocations_during(|| {
        cells = known_cells(&site, "C", &alleles, KnownFlags::default());
    });
    assert!(cells.is_empty());
    assert_eq!(
        not_asked_for, 0,
        "a run without the flags must cost nothing"
    );

    let asked_for = allocations_during(|| {
        cells = known_cells(&site, "C", &alleles, everything);
    });
    assert_eq!(cells.len(), 1);
    let columns = cells[0].cells.iter().count();
    assert!(
        columns >= 14,
        "the record has 14 cells to write, found {columns}"
    );

    // Measured: 147 when every population read built its own scratch lists and a `String`
    // for each value, 79 once those were kept across populations and borrowed, and 67 with
    // a list of one moved rather than joined. What is left is the text of each cell, its list
    // and the lookups of one allele; the bound is above that, and the rows only borrow it.
    eprintln!("known_cells, 14 cells from 12 populations: {asked_for} allocations");
    assert!(
        asked_for <= 80,
        "working out 14 cells took {asked_for} allocations, over the 80 budget: a population read \
         must not build scratch lists or a `String` of its own"
    );
}
