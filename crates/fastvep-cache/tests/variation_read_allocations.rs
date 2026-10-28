//! Pins what the read path of the variation cache costs in heap allocations, per variant.
//!
//! `VariationTabixReader::query` runs once for each variant of a run, on one thread, before the
//! parallel annotation starts, so everything it allocates is paid for by the whole callset in
//! sequence. The sliding window of the reader (a sorted input asks for a window that overlaps the one
//! before it, and its lines are read once) left it the lines that are new to the window, the
//! records a window returns and the chunks the index finds. This holds that, over a sorted sweep
//! of the dense part of the fixture with the windows VEP asks for (a base either side of the
//! variant), and over windows that no record is in reach of, which cost what finding that there is
//! no chunk costs.
//!
//! Measured on the fixture, which has 20 columns and few of them populated: 12.4 allocations a
//! lookup over the dense part (1,283 records returned for 1,500 lookups) and 1 over the gap. The
//! budget is that with a quarter on top, for the allocator and the compiler of the day. It is
//! tight enough that parsing one record more than a lookup returns (the line past the window, say)
//! costs 18.4 and fails it. A real cache has more columns filled, so a lookup costs more there.
//!
//! What it does not hold is how often a line is read: a line read again is not allocated again
//! (the buffer is reused, and only the lines of the window are kept), so a reader that seeks back
//! to the index's chunk start for every variant costs the same allocations as one that carries on.
//! That is held by the unit test of the reader (`a_sorted_sweep_reads_each_line_about_once`), which
//! counts the lines read and saw 11,362,930 of them for 6,667 windows with the carry-on disabled,
//! against a limit of 12,000.
//!
//! This file installs a counting global allocator that counts the allocations of the thread that is
//! measured and of no other: the test harness has a thread of its own that allocates while a test
//! runs, and a count of the whole process took those into the window, which failed a test of this
//! kind on a CI runner.

use fastvep_cache::info::CacheInfo;
use fastvep_cache::variation::VariationTabixReader;
use std::alloc::{GlobalAlloc, Layout, System};
use std::cell::Cell;
use std::path::Path;

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

fn allocations_during<F: FnOnce()>(f: F) -> usize {
    let before = ALLOCS.with(Cell::get);
    f();
    ALLOCS.with(Cell::get) - before
}

#[test]
fn a_sorted_sweep_of_lookups_costs_what_the_lines_and_records_it_meets_cost() {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("fixtures/vep_cache_window");
    let info = CacheInfo::from_file(&dir.join("info.txt")).unwrap();
    let reader = VariationTabixReader::new(&dir, &info.variation_cols).unwrap();

    // The dense part of the fixture: a record every one to seven bases from 100,000, which a
    // sweep of one variant every three bases meets once each. Warm first: the file is opened and
    // its index read on the first query, which is the run's and not a variant's.
    let first = reader.query("7", 100_000, 100_002).unwrap();
    assert!(!first.is_empty());
    let sweep = 1_500u64;
    let mut returned = 0usize;
    let dense = allocations_during(|| {
        for i in 1..=sweep {
            let at = 100_000 + 3 * i;
            returned += reader.query("7", at - 1, at + 1).unwrap().len();
        }
    });
    assert!(returned > 500, "the sweep found {returned} records");
    let per_query = dense as f64 / sweep as f64;
    assert!(
        per_query < BUDGET_DENSE,
        "{per_query:.1} allocations a lookup over the dense part ({dense} in {sweep}, {returned} \
         records returned), where {BUDGET_DENSE} was the budget"
    );

    // A window with nothing in reach (the gap of a megabase between the dense regions): what an
    // index query costs, which is a vector of chunks or none at all.
    let empty = allocations_during(|| {
        for i in 0..sweep {
            let at = 600_000 + 400 * i;
            assert!(reader.query("7", at - 1, at + 1).unwrap().is_empty());
        }
    });
    let per_empty = empty as f64 / sweep as f64;
    assert!(
        per_empty < BUDGET_EMPTY,
        "{per_empty:.1} allocations a lookup over the gap ({empty} in {sweep}), where {BUDGET_EMPTY} was the budget"
    );
}

/// Allocations a lookup may cost on average over the dense sweep. See the module documentation.
const BUDGET_DENSE: f64 = 16.0;
/// ... and over a window in which there is nothing.
const BUDGET_EMPTY: f64 = 2.0;
