//! Pins the allocation behaviour of the intronic 3'-shift.
//!
//! The shift walked the reference one base at a time, calling `fetch_sequence`
//! per step - and twice per step for a deletion, once for each end - with every
//! call returning an owned `Vec`. A variant sliding 60 bases down a repeat cost
//! 60 heap allocations, or 120 as a deletion, and the walk runs once per
//! (variant x transcript x allele): a gene-dense locus with twenty overlapping
//! transcripts pays it twenty times for the same variant.
//!
//! It now reads a block at a time, so the count is set by how many blocks the
//! walk crosses rather than by how far it travels.
//!
//! This file installs a counting global allocator, so it deliberately holds
//! exactly ONE test: `cargo test` runs the tests in a binary concurrently, and a
//! second test allocating on another thread would be counted here.

use anyhow::{anyhow, Result};
use fastvep_annotate::three_prime_shift_intronic;
use fastvep_cache::providers::SequenceProvider;
use fastvep_core::{Allele, Strand};
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

/// Owns its sequence, so a fetch allocates exactly once - the same shape the
/// real FASTA readers have, which return an owned `Vec` per call.
struct MemRef(Vec<u8>);

impl SequenceProvider for MemRef {
    fn fetch_sequence(&self, _chrom: &str, start: u64, end: u64) -> Result<Vec<u8>> {
        if start < 1 || end < start {
            return Err(anyhow!("bad range"));
        }
        let s0 = (start - 1) as usize;
        if s0 >= self.0.len() {
            return Err(anyhow!("past contig end"));
        }
        Ok(self.0[s0..(end as usize).min(self.0.len())].to_vec())
    }
}

#[test]
fn the_shift_reads_the_reference_in_blocks_not_a_base_at_a_time() {
    // 4,000 bases of `TG`, an intron a variant can slide a long way down.
    let reference = MemRef(b"TG".repeat(2_000).to_vec());
    let intron = (1u64, 4_000u64);

    // An insertion travelling the full repeat: 3,999 steps.
    let mut landed = (0, 0);
    let insertion = allocations_during(|| {
        landed = three_prime_shift_intronic(
            &reference,
            "1",
            1,
            0,
            &Allele::Deletion,
            &Allele::Sequence(b"TG".to_vec()),
            Strand::Forward,
            intron.0,
            intron.1,
        );
    });
    // The insertion point ends one past the last intronic base it slid over.
    assert_eq!(landed.0, 4_001, "the walk should cross the whole repeat");

    // A deletion travelling the same distance, reading both of its ends.
    let mut deleted = (0, 0);
    let deletion = allocations_during(|| {
        deleted = three_prime_shift_intronic(
            &reference,
            "1",
            1,
            2,
            &Allele::Sequence(b"TG".to_vec()),
            &Allele::Deletion,
            Strand::Forward,
            intron.0,
            intron.1,
        );
    });
    assert_eq!(deleted.1, 4_000, "the walk should cross the whole repeat");

    // The window opens at 16 bases and doubles to 1,024, so a 4,000-base walk
    // pays a handful of small refills and then a few large ones; the deletion
    // pays that once per end. Per-base reads cost ~4,000 and ~8,000. The budget
    // is loose enough to survive a change to the growth curve and tight enough
    // that a return to per-base reads fails it.
    assert!(
        insertion <= 64,
        "insertion shift allocated {insertion} times crossing 4,000 bases"
    );
    assert!(
        deletion <= 128,
        "deletion shift allocated {deletion} times crossing 4,000 bases"
    );
}
