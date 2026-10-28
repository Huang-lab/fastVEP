//! Pins what asking a transcript where a base is costs in heap allocations: none.
//!
//! The consequence code asks which exon or intron a position is in, what its cDNA position is and
//! whether a stop there escapes decay, several times for every (variant x overlapping transcript x
//! allele). The answers are read off the exons the transcript holds, which the GFF3 reader keeps in
//! transcript order, so nothing here has a reason to copy and sort them. It did: an allocation and a
//! sort of the exons for every one of those questions, which was 11,700 of the 133,328 allocations
//! of annotating 2,000 records of the world fixture.
//!
//! This file installs a counting global allocator that counts the allocations of the thread that is
//! measured and of no other: the test harness has a thread of its own that allocates while a test
//! runs, and a count of the whole process took those into the window, which failed a test of this
//! kind on a CI runner.

use fastvep_core::Strand;
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

/// A transcript of four exons (and so three introns) on `strand`, with the exons stored in the
/// order of the transcript: ascending on the forward strand, descending on the reverse.
fn transcript(strand: Strand) -> Transcript {
    let mut spans = [(100u64, 199u64), (300, 399), (500, 599), (700, 799)];
    if strand == Strand::Reverse {
        spans.reverse();
    }
    let exons = spans
        .iter()
        .enumerate()
        .map(|(i, &(start, end))| Exon {
            stable_id: format!("ENSE_{i}"),
            start,
            end,
            strand,
            phase: -1,
            end_phase: -1,
            rank: i as u32 + 1,
        })
        .collect();
    Transcript {
        stable_id: "ENST_EXONS".into(),
        version: Some(1),
        gene: Gene {
            stable_id: "ENSG_EXONS".into(),
            symbol: None,
            symbol_source: None,
            hgnc_id: None,
            biotype: "protein_coding".into(),
            chromosome: "1".into(),
            start: 100,
            end: 799,
            strand,
        },
        biotype: "protein_coding".into(),
        chromosome: "1".into(),
        start: 100,
        end: 799,
        strand,
        exons,
        translation: None,
        cdna_coding_start: Some(21),
        cdna_coding_end: Some(380),
        coding_region_start: Some(120),
        coding_region_end: Some(780),
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
        initiator_met: false,
        cds_flags_from_sequence: false,
        readthrough: false,
        codon_table_start_phase: 0,
    }
}

#[test]
fn asking_a_transcript_where_a_base_is_allocates_nothing() {
    for strand in [Strand::Forward, Strand::Reverse] {
        let t = transcript(strand);
        // in each exon and each intron, on an edge, before the first exon and after the last
        let positions = [
            50u64, 100, 150, 199, 250, 300, 399, 450, 500, 599, 650, 700, 799, 900,
        ];
        // the answers come first, so that nothing the loop below allocates is counted
        let answers: Vec<_> = positions
            .iter()
            .map(|&p| {
                (
                    t.genomic_to_cdna(p),
                    t.exon_at(p),
                    t.exon_overlapping(p, p + 3),
                    t.intron_at(p),
                    t.intron_bounds_at(p),
                    t.genomic_to_intronic_cdna(p),
                )
            })
            .collect();
        assert!(answers.iter().any(|a| a.0.is_some()), "{strand:?}");
        assert!(answers.iter().any(|a| a.3.is_some()), "{strand:?}");
        assert_eq!(t.escapes_nmd(10), Some(false), "{strand:?}");

        let before = ALLOCS.with(Cell::get);
        let mut seen = 0usize;
        for _ in 0..10_000 {
            for &p in &positions {
                seen += usize::from(t.genomic_to_cdna(p).is_some());
                seen += usize::from(t.exon_at(p).is_some());
                seen += usize::from(t.exon_overlapping(p, p + 3).is_some());
                seen += usize::from(t.intron_at(p).is_some());
                seen += usize::from(t.intron_bounds_at(p).is_some());
                seen += usize::from(t.genomic_to_intronic_cdna(p).is_some());
                seen += usize::from(t.escapes_nmd(p).is_some());
            }
        }
        let allocations = ALLOCS.with(Cell::get) - before;
        assert!(seen > 0);
        assert_eq!(
            allocations, 0,
            "{strand:?}: asking where a base is allocated {allocations} times over 980,000 questions"
        );
    }
}
