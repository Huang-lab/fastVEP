//! Pins what `--total_length` costs in heap allocations: none.
//!
//! The per-variant loop gives a row its totals once per (variant x overlapping transcript x allele),
//! and a run with the flag does that for every row it writes. The totals are three lengths the
//! transcript already has (the exons added up, the translateable sequence, the peptide) and three
//! tests of extent, so nothing here has a reason to build a `String` or a `Vec`.
//!
//! This file installs a counting global allocator, so it deliberately holds exactly ONE test:
//! `cargo test` runs the tests of a binary concurrently, and a second test allocating on another
//! thread would be counted here.

use fastvep_core::Strand;
use fastvep_genome::{Exon, Gene, Transcript, Translation};
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

/// Two exons of 17 and 15 bases with a 20-base intron between them, and a CDS of
/// `ATG GCT AAA GCT TGG TAA` (M A K A W *) from the 10th base of the first to the 10th of the second.
const EXON_1: &str = "ACGTACGTAATGGCTAA";
const EXON_2: &str = "AGCTTGGTAAACGTA";

fn transcript() -> Transcript {
    let (e1, e2) = ((1u64, 17u64), (38u64, 52u64));
    let exon = |n: u32, (start, end): (u64, u64)| Exon {
        stable_id: format!("ENSE_{n}"),
        start,
        end,
        strand: Strand::Forward,
        phase: -1,
        end_phase: -1,
        rank: n,
    };
    let mut t = Transcript {
        stable_id: "ENST_ALLOC".into(),
        version: Some(1),
        gene: Gene {
            stable_id: "ENSG_ALLOC".into(),
            symbol: None,
            symbol_source: None,
            hgnc_id: None,
            biotype: "protein_coding".into(),
            chromosome: "1".into(),
            start: 1,
            end: 52,
            strand: Strand::Forward,
        },
        biotype: "protein_coding".into(),
        chromosome: "1".into(),
        start: 1,
        end: 52,
        strand: Strand::Forward,
        exons: vec![exon(1, e1), exon(2, e2)],
        translation: Some(Translation {
            stable_id: "ENSP_ALLOC".into(),
            genomic_start: 10,
            genomic_end: 47,
            start_exon_rank: 1,
            start_exon_offset: 9,
            end_exon_rank: 2,
            end_exon_offset: 9,
        }),
        cdna_coding_start: Some(10),
        cdna_coding_end: Some(27),
        coding_region_start: Some(10),
        coding_region_end: Some(47),
        spliced_seq: None,
        translateable_seq: None,
        peptide: None,
        canonical: true,
        mane_select: None,
        mane_plus_clinical: None,
        tsl: None,
        appris: None,
        ccds: None,
        protein_id: Some("ENSP_ALLOC".into()),
        protein_version: Some(1),
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
    };
    let genome = format!("{EXON_1}{}{EXON_2}", "N".repeat(20));
    t.build_sequences(|_, start, end| {
        Ok(genome.as_bytes()[(start - 1) as usize..end as usize].to_vec())
    })
    .unwrap();
    t
}

#[test]
fn giving_a_row_its_totals_allocates_nothing() {
    let t = transcript();
    assert_eq!(t.cds_length(), Some(18));
    assert_eq!(t.peptide_length(), Some(5));

    // a change in the 5' UTR, in the first codon, across the intron, in the intron, in the 3' UTR,
    // upstream of the transcript, and an insertion before the initiator
    let changes: [(u64, u64); 7] = [
        (5, 5),
        (11, 11),
        (15, 40),
        (25, 25),
        (50, 50),
        (0, 0),
        (10, 9),
    ];
    let locations: Vec<_> = changes.iter().map(|&(s, e)| t.location(s, e)).collect();

    let before = ALLOCS.load(Ordering::Relaxed);
    let mut totals = 0usize;
    for _ in 0..10_000 {
        for (&(s, e), location) in changes.iter().zip(&locations) {
            let mut location = *location;
            t.apply_total_length(&mut location, s, e);
            totals += usize::from(location.cdna.total.is_some())
                + usize::from(location.cds.total.is_some())
                + usize::from(location.protein.total.is_some());
        }
    }
    let allocations = ALLOCS.load(Ordering::Relaxed) - before;
    assert_eq!(allocations, 0, "{allocations} allocations");
    // the loop did work: five of the seven changes touch an exon, and three of them the CDS
    // too (the two UTR ones are `exon` and not `coding`, the intron one and the one upstream
    // are neither)
    assert_eq!(totals, 10_000 * (5 + 3 + 3), "{totals}");
}
