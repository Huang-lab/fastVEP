//! Reading a record that is not a structural variant costs what its alleles cost, and not what its
//! INFO column does.
//!
//! A record whose alleles are not bases alone (a `*` beside the alternate of a spanning deletion,
//! `N`, a lower-case or an IUPAC base) has to be told from a structural variant, and a first reading
//! of that test built a `Vec` for each of the INFO column's pieces and the joined ALTs, for every
//! such record: 21 allocations for an INFO of a few keys became 173 for one of 150, and the time
//! went from 0.45 to 5.4 microseconds. The test of the INFO is now one search for `SVTYPE`, which a
//! record that is not a structural variant does not have, and the ALTs are joined only when one of
//! them could be symbolic.
//!
//! The allocations of the calling thread only are counted: the test harness has a thread of its own
//! that allocates while a test runs.

use fastvep_io::vcf::parse_vcf_line;
use std::alloc::{GlobalAlloc, Layout, System};
use std::cell::Cell;

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

fn allocations_of(line: &str) -> usize {
    let before = ALLOCS.with(Cell::get);
    let parsed = parse_vcf_line(line).expect("the record should parse");
    let after = ALLOCS.with(Cell::get);
    std::hint::black_box(parsed);
    after - before
}

fn record(alleles: (&str, &str), keys: usize) -> String {
    let info: Vec<String> = (0..keys).map(|k| format!("KEY{k}=value{k}")).collect();
    format!(
        "chr1\t100\t.\t{}\t{}\t.\t.\t{}",
        alleles.0,
        alleles.1,
        info.join(";")
    )
}

#[test]
fn a_record_that_is_not_a_structural_variant_does_not_pay_for_its_info_columns_width() {
    // warm the first-use costs (a lazily built table, a thread-local) before anything is counted
    allocations_of(&record(("A", "C"), 2));
    for alleles in [
        ("A", "C"),
        ("AT", "A,ATT"),
        ("A", "*"),
        ("A", "G,*"),
        ("N", "C"),
        ("a", "c"),
        ("A", "R"),
    ] {
        let narrow = allocations_of(&record(alleles, 2));
        let wide = allocations_of(&record(alleles, 150));
        assert!(
            wide <= narrow,
            "{alleles:?}: {narrow} allocations for an INFO of 2 keys, {wide} for one of 150"
        );
    }
}
