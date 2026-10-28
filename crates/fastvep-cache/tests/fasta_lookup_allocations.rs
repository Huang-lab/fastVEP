//! Pins what asking the in-memory FASTA reader for a contig's length costs in heap allocations, for
//! a name the FASTA has: none.
//!
//! The sequence provider asks for the length of the contig before every fetch (to wrap a region
//! over the origin of a circular one), so a run asks once for each record that overlaps a
//! transcript and writes an HGVS description. The aliases of a name (`chr1` for `1`, `chrM` for
//! `MT`) are a `Vec` of `String`s, which the lookup built for every call before it tried the name
//! it was given, which is the one that is there almost every time.
//!
//! This file installs a counting global allocator that counts the allocations of the thread that is
//! measured and of no other: the test harness has a thread of its own that allocates while a test
//! runs, and a count of the whole process took those into the window, which failed a test of this
//! kind on a CI runner.

use fastvep_cache::fasta::FastaReader;
use std::alloc::{GlobalAlloc, Layout, System};
use std::cell::Cell;

thread_local! {
    // The allocations of the calling thread only. The test harness has a thread of its own that allocates
    // while a test runs, and a count of the whole process took those into the window that is measured.
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

    let before = ALLOCS.with(Cell::get);
    let mut total = 0u64;
    for _ in 0..10_000 {
        for name in ["1", "chrX", "MT"] {
            total += reader.sequence_length(name).unwrap();
        }
    }
    let allocations = ALLOCS.with(Cell::get) - before;
    assert_eq!(total, 340_000);
    assert_eq!(
        allocations, 0,
        "30,000 lookups of names the FASTA has allocated {allocations} times"
    );
}
