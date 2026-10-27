//! Pins what asking the in-memory FASTA reader for a contig's length costs in heap allocations, for
//! a name the FASTA has: none.
//!
//! The sequence provider asks for the length of the contig before every fetch (to wrap a region
//! over the origin of a circular one), so a run asks once for each record that overlaps a
//! transcript and writes an HGVS description. The aliases of a name (`chr1` for `1`, `chrM` for
//! `MT`) are a `Vec` of `String`s, which the lookup built for every call before it tried the name
//! it was given, which is the one that is there almost every time.
//!
//! This file installs a counting global allocator, so it deliberately holds exactly ONE test:
//! `cargo test` runs the tests of a binary concurrently, and a second test allocating on another
//! thread would be counted here.

use fastvep_cache::fasta::FastaReader;
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

#[test]
fn the_length_of_a_contig_the_fasta_has_is_asked_for_without_an_allocation() {
    let fasta = ">1\nACGTACGTAC\n>chrX\nACGT\n>MT\nACGTACGTACGTACGTACGT\n";
    let reader = FastaReader::from_reader(fasta.as_bytes()).unwrap();

    // the names the file has, and the spellings that reach them through the aliases
    assert_eq!(reader.sequence_length("1"), Some(10));
    assert_eq!(reader.sequence_length("chrX"), Some(4));
    assert_eq!(reader.sequence_length("MT"), Some(20));
    assert_eq!(reader.sequence_length("chr1"), Some(10));
    assert_eq!(reader.sequence_length("X"), Some(4));
    assert_eq!(reader.sequence_length("chrM"), Some(20));
    assert_eq!(reader.sequence_length("2"), None);

    let before = ALLOCS.load(Ordering::Relaxed);
    let mut total = 0u64;
    for _ in 0..10_000 {
        for name in ["1", "chrX", "MT"] {
            total += reader.sequence_length(name).unwrap();
        }
    }
    let allocations = ALLOCS.load(Ordering::Relaxed) - before;
    assert_eq!(total, 340_000);
    assert_eq!(
        allocations, 0,
        "30,000 lookups of names the FASTA has allocated {allocations} times"
    );
}
