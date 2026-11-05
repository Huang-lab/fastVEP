//! Pins the allocation behaviour of the classifier's reading of a gene's disease annotation.
//!
//! PM2, BS2, PM3, BP2 and the validity gate ask what a gene's diseases say about inheritance and
//! validity, once for every allele and transcript of every variant, so a question that built a
//! `String` or a `Vec` would be built millions of times in a run. `crate::disease` reads the stored
//! phenotype strings in place and says so in its documentation; this test holds it to that. The
//! calls that do allocate, `disorder_contexts` and `restricted_annotation`, are asked only of a gene
//! with more than one mode of inheritance, after `has_several_disorder_contexts` has said so, which
//! is itself answered without an allocation.
//!
//! The counting allocator counts the allocations of the thread that is measured and of no other: the
//! test harness has a thread of its own that allocates while a test runs, and a count of the whole
//! process took those into the window (see `fastvep-sa/tests/query_path_allocations.rs`).

use fastvep_classification::sa_extract::OmimData;
use std::alloc::{GlobalAlloc, Layout, System};
use std::cell::Cell;
use std::hint::black_box;

thread_local! {
    static ALLOCS: Cell<usize> = const { Cell::new(0) };
}

struct Counting;

unsafe impl GlobalAlloc for Counting {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        ALLOCS.with(|count| count.set(count.get() + 1));
        unsafe { System.alloc(layout) }
    }
    unsafe fn dealloc(&self, ptr: *mut u8, layout: Layout) {
        unsafe { System.dealloc(ptr, layout) }
    }
    unsafe fn realloc(&self, ptr: *mut u8, layout: Layout, new_size: usize) -> *mut u8 {
        ALLOCS.with(|count| count.set(count.get() + 1));
        unsafe { System.realloc(ptr, layout, new_size) }
    }
}

#[global_allocator]
static ALLOCATOR: Counting = Counting;

/// Allocations performed while running `f`.
fn allocations_during<F: FnOnce()>(f: F) -> usize {
    let before = ALLOCS.with(Cell::get);
    f();
    ALLOCS.with(Cell::get) - before
}

fn gene(phenotypes: &[&str]) -> OmimData {
    OmimData {
        mim_number: None,
        phenotypes: Some(phenotypes.iter().map(|s| s.to_string()).collect()),
    }
}

/// A gene as the ClinGen build writes one: an autosomal dominant and an autosomal recessive disease
/// that are established, an X-linked one, a weak one, and the untagged piece that the builder makes of
/// a label with a `;`.
fn clingen_gene_with_several_modes() -> OmimData {
    gene(&[
        "Bohring-Opitz syndrome (ClinGen Definitive/AD, MONDO:0011510)",
        "a recessive disease of the same gene (ClinGen Strong/AR, MONDO:0000002)",
        "an X-linked disorder, said to be dominant (ClinGen Moderate/XL, MONDO:0000003)",
        "a disputed relationship (ClinGen Disputed/AD, MONDO:0000004)",
        "lung disease, immunodeficiency, and chromosome breakage syndrome",
        " (ClinGen Limited/AR, MONDO:0014984)",
    ])
}

/// A gene as an OMIM `genemap2.txt` build writes one: free text, no tag.
fn omim_gene() -> OmimData {
    gene(&[
        "Charcot-Marie-Tooth disease, axonal, type 2N, 613287 (3), Autosomal recessive",
        "Spastic paraplegia 7, 607259 (3), Autosomal dominant",
        "Deafness, autosomal dominant 3A, 601544 (3)",
    ])
}

#[test]
fn reading_what_a_gene_says_about_inheritance_and_validity_allocates_nothing() {
    for (name, omim) in [
        ("ClinGen", clingen_gene_with_several_modes()),
        ("OMIM", omim_gene()),
    ] {
        // Anything the first call sets up once is not what this measures.
        black_box(omim.inheritance());
        let allocations = allocations_during(|| {
            for _ in 0..1_000 {
                black_box(omim.inheritance());
                black_box(omim.has_dominant_inheritance());
                black_box(omim.has_recessive_inheritance());
                black_box(omim.has_established_relationship());
                black_box(omim.has_no_valid_relationship());
                black_box(omim.has_several_disorder_contexts());
                black_box(omim.entries().count());
            }
        });
        assert_eq!(
            allocations, 0,
            "{name}: reading the annotation allocated {allocations} times in 1,000 rounds of seven questions; \
             it is asked once for each allele and transcript of every variant"
        );
    }
}

/// The calls that do allocate are the ones a gene with one mode of inheritance never makes, because
/// `has_several_disorder_contexts` is false for it and allocates nothing in saying so.
#[test]
fn a_gene_with_one_mode_of_inheritance_is_known_to_be_one_without_an_allocation() {
    let single = gene(&[
        "disease one (ClinGen Definitive/AD, MONDO:0000001)",
        "disease two (ClinGen Moderate/AD, MONDO:0000002)",
        "disease three (ClinGen Limited/AR, MONDO:0000003)",
    ]);
    black_box(single.has_several_disorder_contexts());
    let allocations = allocations_during(|| {
        for _ in 0..1_000 {
            assert!(!black_box(&single).has_several_disorder_contexts());
        }
    });
    assert_eq!(allocations, 0);

    let several = clingen_gene_with_several_modes();
    assert!(several.has_several_disorder_contexts());
    // The counter counts: a call that builds the contexts is seen to allocate, so a zero above is a
    // measurement and not a broken allocator.
    let mut contexts = Vec::new();
    let allocations = allocations_during(|| contexts = several.disorder_contexts());
    assert!(allocations > 0);
    assert_eq!(
        contexts.len(),
        3,
        "AD, AR and XL: the weak ones make no context"
    );
}
