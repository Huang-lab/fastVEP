//! Pins what the HGVSp of an allele costs in heap allocations.
//!
//! The per-variant loop asks for it once per (variant x overlapping transcript x allele), and
//! nearly every one of those is a change that has no protein notation at all (an intron, a UTR,
//! a flank) or a missense in a gene that the rest of the run is already busy with. The work is
//! done in the order that lets the common ones leave early: the checks that cost nothing first,
//! the name of the protein, which is a `String`, last.
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

/// Nine bases of UTR, a CDS of `ATG GCT AAA GCT TGG TAA` (M A K A W *), and five more of UTR,
/// on one forward-strand exon.
const GENOME: &str = "ACGTACGTAATGGCTAAAGCTTGGTAAACGTA";

fn transcript() -> Transcript {
    let end = GENOME.len() as u64;
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
            end,
            strand: Strand::Forward,
        },
        biotype: "protein_coding".into(),
        chromosome: "1".into(),
        start: 1,
        end,
        strand: Strand::Forward,
        exons: vec![Exon {
            stable_id: "ENSE_ALLOC".into(),
            start: 1,
            end,
            strand: Strand::Forward,
            phase: -1,
            end_phase: -1,
            rank: 1,
        }],
        translation: Some(Translation {
            stable_id: "ENSP_ALLOC".into(),
            genomic_start: 10,
            genomic_end: 27,
            start_exon_rank: 1,
            start_exon_offset: 9,
            end_exon_rank: 1,
            end_exon_offset: 26,
        }),
        cdna_coding_start: Some(10),
        cdna_coding_end: Some(27),
        coding_region_start: Some(10),
        coding_region_end: Some(27),
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
        codon_table_start_phase: 0,
        initiator_met: false,
        cds_flags_from_sequence: false,
        readthrough: false,
    };
    let fasta = format!(">1\n{GENOME}\n");
    let provider = FastaSequenceProvider::new(FastaReader::from_reader(fasta.as_bytes()).unwrap());
    t.build_sequences(|chrom, start, end| {
        provider
            .fetch_sequence(chrom, start, end)
            .map_err(|e| e.to_string())
    })
    .unwrap();
    t
}

/// The notation of one change, and what it cost.
fn hgvsp(
    predictor: &ConsequencePredictor,
    transcript: &Transcript,
    start: u64,
    end: u64,
    ref_allele: Allele,
    alt_allele: Allele,
) -> (Option<String>, usize) {
    let position = GenomicPosition::new("1", start, end, Strand::Forward);
    let result = predictor.predict(
        &position,
        &ref_allele,
        std::slice::from_ref(&alt_allele),
        &[transcript],
        None,
    );
    let ac = &result.transcript_consequences[0].allele_consequences[0];
    allocations_during(|| {
        hgvsp_for_allele(&HgvspRequest {
            predictor,
            transcript,
            provider: None,
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
        })
    })
}

#[test]
fn a_notation_costs_what_it_has_to_and_a_change_without_one_costs_nothing() {
    let predictor = ConsequencePredictor::new(5000, 5000);
    let transcript = transcript();
    assert_eq!(transcript.peptide.as_deref(), Some("MAKAW*"));

    // In the UTR, and past the end of the transcript: no protein, nothing allocated
    let (none, utr) = hgvsp(
        &predictor,
        &transcript,
        3,
        3,
        Allele::Sequence(b"G".to_vec()),
        Allele::Sequence(b"T".to_vec()),
    );
    assert_eq!(none, None);
    assert_eq!(utr, 0, "a UTR change allocated {utr} times");
    let (none, flank) = hgvsp(
        &predictor,
        &transcript,
        40,
        40,
        Allele::Sequence(b"A".to_vec()),
        Allele::Sequence(b"T".to_vec()),
    );
    assert_eq!(none, None);
    assert_eq!(
        flank, 0,
        "a change past the transcript allocated {flank} times"
    );

    // A missense (G>A at 13, GCT to ACT): the protein's name and the notation
    let (missense, cost) = hgvsp(
        &predictor,
        &transcript,
        13,
        13,
        Allele::Sequence(b"G".to_vec()),
        Allele::Sequence(b"A".to_vec()),
    );
    assert_eq!(missense.as_deref(), Some("ENSP_ALLOC.1:p.Ala2Thr"));
    assert!(cost <= 2, "a missense allocated {cost} times");

    // An in-frame deletion of the Lys (16-18): the general path, which builds strings, but
    // not one per base of the transcript
    let (deletion, cost) = hgvsp(
        &predictor,
        &transcript,
        16,
        18,
        Allele::Sequence(b"AAA".to_vec()),
        Allele::Deletion,
    );
    assert_eq!(deletion.as_deref(), Some("ENSP_ALLOC.1:p.Lys3del"));
    assert!(cost <= 40, "an in-frame deletion allocated {cost} times");
}
