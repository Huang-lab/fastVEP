//! Regression tests for the backward transcript scan (issue #126).
//!
//! `IndexedTranscriptProvider` prunes its backward scan with a per-chromosome
//! max-end array. Getting that array's direction wrong silently drops every
//! transcript that is still open past the last transcript to *begin* on the
//! chromosome — a shape that exists in every real gene model but in none of the
//! hand-built fixtures, so it survived the whole unit suite.
//!
//! These tests are hermetic: the first pins the real Ensembl shape by
//! coordinate, the second searches for the shape at random.

use fastvep_cache::providers::{
    IndexedTranscriptProvider, MemoryTranscriptProvider, TranscriptProvider,
};
use fastvep_core::Strand;
use fastvep_genome::{Exon, Gene, Transcript};
use std::collections::BTreeSet;
use std::sync::Arc;

fn transcript(id: &str, chrom: &str, start: u64, end: u64) -> Transcript {
    Transcript {
        stable_id: Arc::from(id),
        version: None,
        gene: Gene {
            stable_id: "ENSG_1".into(),
            symbol: None,
            symbol_source: None,
            hgnc_id: None,
            biotype: "lncRNA".into(),
            chromosome: chrom.into(),
            start,
            end,
            strand: Strand::Forward,
        },
        biotype: "lncRNA".into(),
        chromosome: chrom.into(),
        start,
        end,
        strand: Strand::Forward,
        exons: vec![Exon {
            stable_id: "ENSE_1".into(),
            start,
            end,
            strand: Strand::Forward,
            phase: 0,
            end_phase: 0,
            rank: 1,
        }],
        translation: None,
        cdna_coding_start: None,
        cdna_coding_end: None,
        coding_region_start: None,
        coding_region_end: None,
        spliced_seq: None,
        translateable_seq: None,
        peptide: None,
        canonical: false,
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
    }
}

fn ids(rows: Vec<&Transcript>) -> BTreeSet<String> {
    rows.into_iter().map(|t| t.stable_id.to_string()).collect()
}

/// The 29 transcripts on Ensembl 116 chr22 that stay open past the end of
/// ENST00000427528, the last transcript on the chromosome to begin. A query
/// anywhere in 50,799,124-50,801,309 overlaps some of them, and the whole
/// 2,186 bp window was reported as intergenic before this was fixed.
const CHR22_TAIL: &[(&str, u64, u64)] = &[
    ("ENST00000467796", 50756831, 50799540),
    ("ENST00000826711", 50756876, 50799436),
    ("ENST00000826712", 50757134, 50799126),
    ("ENST00000462238", 50783676, 50799667),
    ("ENST00000687913", 50783702, 50799539),
    ("ENST00000651346", 50783719, 50799546),
    ("ENST00000480246", 50783740, 50801309),
    ("ENST00000826714", 50783752, 50799434),
    ("ENST00000690684", 50783754, 50799523),
    ("ENST00000826720", 50783755, 50799178),
    ("ENST00000687360", 50783763, 50799436),
    ("ENST00000826722", 50783764, 50799178),
    ("ENST00000826715", 50783764, 50799434),
    ("ENST00000701498", 50783766, 50799213),
    ("ENST00000826718", 50783766, 50799335),
    ("ENST00000701520", 50783769, 50799360),
    ("ENST00000700830", 50783773, 50799481),
    ("ENST00000826721", 50783775, 50799193),
    ("ENST00000826723", 50783785, 50799178),
    ("ENST00000826716", 50783785, 50799434),
    ("ENST00000826713", 50783785, 50799539),
    ("ENST00000826717", 50783799, 50799434),
    ("ENST00000701035", 50783809, 50799178),
    ("ENST00000826724", 50783817, 50799193),
    ("ENST00000826726", 50783832, 50799190),
    ("ENST00000826725", 50783832, 50799193),
    ("ENST00000826728", 50783850, 50799125),
    ("ENST00000826727", 50783875, 50799181),
    ("ENST00000826719", 50783878, 50799430),
    ("ENST00000427528", 50798655, 50799123),
];

#[test]
fn chr22_tail_window_is_not_reported_as_intergenic() {
    let transcripts: Vec<Transcript> = CHR22_TAIL
        .iter()
        .map(|&(id, s, e)| transcript(id, "22", s, e))
        .collect();
    let indexed = IndexedTranscriptProvider::new(transcripts.clone());
    let linear = MemoryTranscriptProvider::new(transcripts);

    // Every base of the affected window, plus its two boundaries.
    for pos in [
        50_799_123u64,
        50_799_124,
        50_800_000,
        50_801_309,
        50_801_310,
    ] {
        let expected = ids(linear.get_transcripts("22", pos, pos).unwrap());
        let actual = ids(indexed.get_transcripts("22", pos, pos).unwrap());
        assert_eq!(actual, expected, "22:{pos}");
    }

    // The window's interior must never come back empty.
    for pos in (50_799_124..=50_801_309).step_by(7) {
        let actual = indexed.get_transcripts("22", pos, pos).unwrap();
        assert!(
            !actual.is_empty(),
            "22:{pos} lies inside 29 transcripts but the index returned none"
        );
    }
}

/// Random nested interval sets, checked against a plain filter. The generator
/// deliberately produces long transcripts that outlive every later start, which
/// is the shape the hand-built fixtures never had.
#[test]
fn indexed_provider_matches_plain_filter_on_random_nested_sets() {
    // A fixed-seed LCG keeps this deterministic without a proptest dependency.
    let mut state = 0x2545_F491_4F6C_DD1Du64;
    let mut next = move |bound: u64| {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state % bound
    };

    for case in 0..500 {
        let n = 1 + next(40) as usize;
        let mut transcripts = Vec::with_capacity(n);
        for i in 0..n {
            let start = 1 + next(10_000);
            // Mostly short, occasionally long enough to span many later starts.
            let len = if next(4) == 0 { next(9_000) } else { next(200) };
            transcripts.push(transcript(&format!("T{i}"), "chr1", start, start + len));
        }
        let indexed = IndexedTranscriptProvider::new(transcripts.clone());
        let linear = MemoryTranscriptProvider::new(transcripts);

        for _ in 0..200 {
            let start = 1 + next(11_000);
            let end = start + next(300);
            let expected = ids(linear.get_transcripts("chr1", start, end).unwrap());
            let actual = ids(indexed.get_transcripts("chr1", start, end).unwrap());
            assert_eq!(actual, expected, "case {case}, query {start}-{end}");
        }
    }
}

/// An empty chromosome and a single-transcript chromosome must not panic on the
/// max-end array's boundary.
#[test]
fn degenerate_chromosomes_are_safe() {
    let empty = IndexedTranscriptProvider::new(vec![]);
    assert!(empty.get_transcripts("chr1", 1, 10).unwrap().is_empty());

    let single = IndexedTranscriptProvider::new(vec![transcript("T", "chr1", 100, 200)]);
    assert!(single.get_transcripts("chr1", 1, 99).unwrap().is_empty());
    assert_eq!(single.get_transcripts("chr1", 100, 100).unwrap().len(), 1);
    assert_eq!(single.get_transcripts("chr1", 200, 200).unwrap().len(), 1);
    assert!(single.get_transcripts("chr1", 201, 300).unwrap().is_empty());
}
