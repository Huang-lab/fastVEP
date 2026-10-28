//! Pins what the HGVSp of an allele costs in heap allocations when there is a genome to walk over.
//!
//! `hgvsp_allocations.rs` does it with no sequence provider, where the 3'-rule has nothing to read.
//! With one, an insertion or a deletion is carried over the repeat it sits in by a walk over the
//! genome (`shift_length`), and when the walk crosses an intron the same change is placed again on
//! the transcript (`crossed_an_intron`, `walk_the_transcript`). The loop asks once per (variant x
//! overlapping transcript x allele), so what a request allocates is paid millions of times, and the
//! walk is the part whose cost could grow with the length of the repeat. This holds the number for
//! each kind of request, and holds that a repeat nearly six times as long costs about what the
//! short one does.
//!
//! This file installs a counting global allocator, so it deliberately holds exactly ONE test:
//! `cargo test` runs the tests of a binary concurrently, and a second test allocating on another
//! thread would be counted here.

use fastvep_annotate::{hgvsp_for_allele, HgvspRequest};
use fastvep_cache::fasta::FastaReader;
use fastvep_cache::providers::{FastaSequenceProvider, SequenceProvider};
use fastvep_consequence::ConsequencePredictor;
use fastvep_core::{Allele, GenomicPosition, Strand};
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

fn allocations_during<T>(f: impl FnOnce() -> T) -> (T, usize) {
    let before = ALLOCS.load(Ordering::Relaxed);
    let out = f();
    (out, ALLOCS.load(Ordering::Relaxed) - before)
}

/// Nine bases of 5' UTR, `ATG GCT AAG ACA CAC` and the first exon's end, then a one-base intron
/// whose `A` carries on the `AC` repeat the exon ends in, and a second exon that starts with the
/// rest of it.
const EXON1: &str = "ACGTACGTAATGGCTAAGACACAC";
const INTRON: &str = "A";
/// the repeat, `CACAC` and `AC` times `units`, then a stop codon in frame and a 3' UTR
fn exon2(units: usize) -> String {
    // 15 bases of CDS in the first exon, 5 + 2 * units here, and `x` more to make the stop in frame
    let x = (3 - (15 + 5 + 2 * units) % 3) % 3;
    format!("CACAC{}{}TAAACGTAC", "AC".repeat(units), "G".repeat(x))
}

struct World {
    transcript: Transcript,
    provider: FastaSequenceProvider,
}

fn world(units: usize) -> World {
    let (exon2, genome) = {
        let e2 = exon2(units);
        (e2.clone(), format!("{EXON1}{INTRON}{e2}"))
    };
    let (len1, len2) = (EXON1.len() as u64, exon2.len() as u64);
    let end = len1 + INTRON.len() as u64 + len2;
    let cds_end = (len1 + len2) as usize - "ACGTAC".len();
    let place = |i: usize| -> u64 {
        if (i as u64) < len1 {
            i as u64 + 1
        } else {
            i as u64 + INTRON.len() as u64 + 1
        }
    };
    let (coding_start, coding_end) = (place(9), place(cds_end - 1));
    let mut transcript = Transcript {
        stable_id: "ENST_WALK".into(),
        version: Some(1),
        gene: Gene {
            stable_id: "ENSG_WALK".into(),
            symbol: None,
            symbol_source: None,
            hgnc_id: None,
            biotype: "protein_coding".into(),
            chromosome: "1".into(),
            start: 1,
            end,
            strand: Strand::Forward,
        },
        biotype: "protein_coding".into(),
        chromosome: "1".into(),
        start: 1,
        end,
        strand: Strand::Forward,
        exons: vec![
            Exon {
                stable_id: "E0".into(),
                start: 1,
                end: len1,
                strand: Strand::Forward,
                phase: -1,
                end_phase: -1,
                rank: 1,
            },
            Exon {
                stable_id: "E1".into(),
                start: len1 + INTRON.len() as u64 + 1,
                end,
                strand: Strand::Forward,
                phase: -1,
                end_phase: -1,
                rank: 2,
            },
        ],
        translation: Some(Translation {
            stable_id: "ENSP_WALK".into(),
            genomic_start: coding_start,
            genomic_end: coding_end,
            start_exon_rank: 1,
            start_exon_offset: 9,
            end_exon_rank: 1,
            end_exon_offset: cds_end as u64,
        }),
        cdna_coding_start: Some(10),
        cdna_coding_end: Some(cds_end as u64),
        coding_region_start: Some(coding_start),
        coding_region_end: Some(coding_end),
        spliced_seq: None,
        translateable_seq: None,
        peptide: None,
        canonical: true,
        mane_select: None,
        mane_plus_clinical: None,
        tsl: None,
        appris: None,
        ccds: None,
        protein_id: Some("ENSP_WALK".into()),
        protein_version: Some(1),
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
    let fasta = format!(">1\n{genome}\n");
    let provider = FastaSequenceProvider::new(FastaReader::from_reader(fasta.as_bytes()).unwrap());
    transcript
        .build_sequences(|chrom, start, end| {
            provider
                .fetch_sequence(chrom, start, end)
                .map_err(|e| e.to_string())
        })
        .unwrap();
    World {
        transcript,
        provider,
    }
}

/// The notation of one change with the genome to walk over, and what it cost.
fn hgvsp(
    w: &World,
    start: u64,
    end: u64,
    ref_allele: Allele,
    alt_allele: Allele,
) -> (Option<String>, usize) {
    let predictor = ConsequencePredictor::new(5000, 5000);
    let position = GenomicPosition::new("1", start, end, Strand::Forward);
    let result = predictor.predict(
        &position,
        &ref_allele,
        std::slice::from_ref(&alt_allele),
        &[&w.transcript],
        None,
    );
    let ac = &result.transcript_consequences[0].allele_consequences[0];
    allocations_during(|| {
        hgvsp_for_allele(&HgvspRequest {
            predictor: &predictor,
            transcript: &w.transcript,
            provider: Some(&w.provider as &dyn SequenceProvider),
            chrom: "1",
            start,
            end,
            ref_allele: &ref_allele,
            alt_allele: &alt_allele,
            location: &ac.location,
            consequences: &ac.consequences,
            amino_acids: ac
                .amino_acids
                .as_ref()
                .map(|(r, a)| (r.as_str(), a.as_str())),
            coding_sequence_kept: ac.coding_sequence_kept,
        })
    })
}

#[test]
fn a_walk_over_the_genome_costs_what_it_has_to_and_not_what_the_repeat_is_long() {
    let short = world(2);
    let long = world(40);
    assert!(short
        .transcript
        .peptide
        .as_deref()
        .is_some_and(|p| p.ends_with('*')));
    assert!(long
        .transcript
        .peptide
        .as_deref()
        .is_some_and(|p| p.ends_with('*')));
    let seq = |s: &str| Allele::Sequence(s.as_bytes().to_vec());

    // warm: the first request of a transcript is the run's, not a variant's
    let _ = hgvsp(&short, 3, 3, seq("G"), seq("T"));

    // a change in the UTR touches no residue: nothing to walk, nothing allocated
    let (none, utr) = hgvsp(&short, 3, 3, seq("G"), seq("T"));
    assert_eq!(none, None);
    assert_eq!(utr, 0, "a UTR change allocated {utr} times");

    // a missense is not carried anywhere: the name of the protein and the notation
    let (missense, snv) = hgvsp(&short, 13, 13, seq("G"), seq("A"));
    assert!(missense.is_some());
    assert!(snv <= BUDGET_SNV, "a missense allocated {snv} times");

    // `AC` inserted in the repeat inside the second exon, which has no intron behind it to cross
    let (inside, within) = hgvsp(&short, 28, 27, Allele::Deletion, seq("AC"));
    assert!(inside.is_some());
    assert!(
        within <= BUDGET_WITHIN,
        "an insertion walked over the genome allocated {within} times"
    );

    // `AC` inserted at the start of the repeat the first exon ends in, which the walk over the genome
    // carries across the intron and into the second exon: the change is placed again on the transcript
    let (crossed, across) = hgvsp(&short, 19, 18, Allele::Deletion, seq("AC"));
    assert!(crossed.is_some());
    assert!(
        across <= BUDGET_ACROSS,
        "an insertion carried across an intron allocated {across} times"
    );

    // ... and over a repeat of 92 bases and not 16: the walk reads the bases it walks over from one
    // read of the genome, so the repeat being longer costs the vector that holds them growing
    // and not an allocation for each base
    let (crossed_long, across_long) = hgvsp(&long, 19, 18, Allele::Deletion, seq("AC"));
    assert!(crossed_long.is_some());
    let (inside_long, within_long) = hgvsp(&long, 28, 27, Allele::Deletion, seq("AC"));
    assert!(inside_long.is_some());
    assert!(
        within_long <= BUDGET_WITHIN && across_long <= BUDGET_ACROSS + SLACK_FOR_A_LONG_REPEAT,
        "a repeat of 92 bases: {within_long} and {across_long} allocations, against {within} and {across} for one of 16"
    );
}

// Measured on this transcript (a debug build, as the others of this directory are held): 0 for a
// change in the UTR, 2 for a missense, 31 for an insertion walked over the genome within an exon
// (the same for a repeat of 16 bases and one of 92), 47 for one carried across the intron (57 for
// the longer repeat). The budgets are those with a quarter on top, for the allocator and the
// compiler of the day: a walk that allocated for each base it passes would add 76 to the last two.
const BUDGET_SNV: usize = 4;
const BUDGET_WITHIN: usize = 40;
const BUDGET_ACROSS: usize = 60;
const SLACK_FOR_A_LONG_REPEAT: usize = 15;
