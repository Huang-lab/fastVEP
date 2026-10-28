//! Pins the allocation behaviour of `splice_context`.
//!
//! It runs once per reported (variant x transcript x allele) under `--acmg`, so
//! a `String` or `Vec` built here is built millions of times. Two things are
//! held:
//!
//! * the part every row pays - the intronic offset, and the exon edge for a
//!   synonymous change - allocates nothing, because it is read straight off the
//!   transcript's exon list rather than through a sorted copy of it;
//! * the placement, which walks the reference, is taken only for an indel whose
//!   terms ask for it, and costs a handful of allocations (the clipped alleles
//!   and one block of reference) that do not grow with how far the change slides.
//!
//! This file installs a counting global allocator that counts the allocations of the thread that is
//! measured and of no other: the test harness has a thread of its own that allocates while a test
//! runs, and a count of the whole process took those into the window, which failed a test of this
//! kind on a CI runner.

use anyhow::{anyhow, Result};
use fastvep_annotate::splice_context;
use fastvep_cache::providers::SequenceProvider;
use fastvep_core::{Allele, Consequence, HgvsPlacement, Strand};
use fastvep_genome::{Exon, Gene, Transcript};
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

fn allocations_during<F: FnOnce()>(f: F) -> usize {
    let before = ALLOCS.with(Cell::get);
    f();
    ALLOCS.with(Cell::get) - before
}

/// Owns its sequence, so a fetch allocates exactly once - the shape the real
/// FASTA readers have.
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

/// Exon 1 is 1..=20, the intron 21..=2080 and exon 2 2081..=2100. The intron
/// opens with 2,000 bases of `AGGT` repeat, so an insertion of `AGGT` slides the
/// whole way down it and stops where the repeat does.
fn fixture() -> (Transcript, MemRef) {
    let reference = MemRef(
        format!(
            "{}{}{}",
            "A".repeat(20),
            "AGGT".repeat(500) + &"C".repeat(60),
            "C".repeat(20)
        )
        .as_bytes()
        .to_vec(),
    );
    let strand = Strand::Forward;
    let exon = |start: u64, end: u64, rank: u32| Exon {
        stable_id: format!("ENSE{rank}"),
        start,
        end,
        strand,
        phase: 0,
        end_phase: 0,
        rank,
    };
    let transcript = Transcript {
        stable_id: "ENST00000000001".into(),
        version: Some(1),
        gene: Gene {
            stable_id: "ENSG00000000001".into(),
            symbol: Some("TEST".into()),
            symbol_source: None,
            hgnc_id: None,
            biotype: "protein_coding".into(),
            chromosome: "1".into(),
            start: 1,
            end: 2100,
            strand,
        },
        biotype: "protein_coding".into(),
        chromosome: "1".into(),
        start: 1,
        end: 2100,
        strand,
        exons: vec![exon(1, 20, 1), exon(2081, 2100, 2)],
        translation: None,
        cdna_coding_start: Some(1),
        cdna_coding_end: Some(40),
        coding_region_start: None,
        coding_region_end: None,
        spliced_seq: None,
        translateable_seq: None,
        peptide: None,
        canonical: true,
        mane_select: None,
        mane_plus_clinical: None,
        tsl: None,
        appris: None,
        ccds: None,
        protein_id: None,
        protein_version: None,
        swissprot: vec![],
        trembl: vec![],
        uniparc: vec![],
        refseq_id: None,
        source: None,
        gencode_primary: false,
        flags: vec![],
        codon_table_start_phase: 0,
        initiator_met: false,
        cds_flags_from_sequence: false,
        readthrough: false,
    };
    (transcript, reference)
}

#[test]
fn the_context_a_row_pays_for_does_not_allocate_and_the_placement_does_not_scale_with_the_slide() {
    let (transcript, reference) = fixture();
    let a = Allele::Sequence(b"A".to_vec());
    let g = Allele::Sequence(b"G".to_vec());

    // Every row pays for the offset. Warm everything first so one-time costs
    // (lazy statics, the test harness) are not counted against it.
    let _ = splice_context(Some(&reference), "1", &transcript, 25, 25, &a, &g, &[]);

    let mut intronic = None;
    let offset_only = allocations_during(|| {
        intronic = Some(splice_context(
            Some(&reference),
            "1",
            &transcript,
            25,
            25,
            &a,
            &g,
            &[Consequence::IntronVariant],
        ));
    });
    assert_eq!(intronic.unwrap().intron_offset, Some(5));
    assert_eq!(
        offset_only, 0,
        "the offset is read off the exon list in place"
    );

    // A synonymous change adds the exon edge, which is the same walk.
    let mut edge = None;
    let with_exon_edge = allocations_during(|| {
        edge = Some(splice_context(
            Some(&reference),
            "1",
            &transcript,
            20,
            20,
            &a,
            &g,
            &[Consequence::SynonymousVariant],
        ));
    });
    assert_eq!(edge.unwrap().at_exon_edge, Some(true));
    assert_eq!(with_exon_edge, 0);

    // A substitution with a null term is not walked for a placement.
    let mut substitution = None;
    let no_walk = allocations_during(|| {
        substitution = Some(splice_context(
            Some(&reference),
            "1",
            &transcript,
            21,
            21,
            &a,
            &g,
            &[Consequence::SpliceDonorVariant],
        ));
    });
    assert_eq!(substitution.unwrap().hgvs_placement, None);
    assert_eq!(
        no_walk, 0,
        "a substitution has nowhere to be placed but where it is"
    );

    // An insertion that slides 2,000 bases down the repeat, and one that slides
    // four, and both come to rest where the repeat ends: 60 bases from the
    // acceptor, so the placement is the nearer of the two splice sites. The
    // count is set by the blocks of reference read, not by the distance, so the
    // long slide costs a few blocks more and not 500 times as much.
    let unit = Allele::Sequence(b"AGGT".to_vec());
    let slide = |anchor_end: u64| {
        let mut placed = None;
        let n = allocations_during(|| {
            placed = Some(splice_context(
                Some(&reference),
                "1",
                &transcript,
                anchor_end + 1,
                anchor_end,
                &Allele::Deletion,
                &unit,
                &[Consequence::FrameshiftVariant],
            ));
        });
        (placed.unwrap(), n)
    };
    let (short, short_allocations) = slide(2016);
    let (long, long_allocations) = slide(20);
    assert!(
        long.hgvs_placement == Some(HgvsPlacement::Intronic { offset: -60 }),
        "{:?}",
        long.hgvs_placement
    );
    assert!(
        short.hgvs_placement == Some(HgvsPlacement::Intronic { offset: -60 }),
        "{:?}",
        short.hgvs_placement
    );
    assert!(
        short_allocations <= 6,
        "a short slide cost {short_allocations} allocations"
    );
    assert!(
        long_allocations <= short_allocations + 12,
        "a 2,000-base slide cost {long_allocations} allocations against {short_allocations} for a short one"
    );
}
