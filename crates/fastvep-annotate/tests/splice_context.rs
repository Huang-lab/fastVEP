//! What the classification pass is told about where an allele sits against its
//! transcript: `splice_context`, and the lookup that finds each surviving row's
//! transcript after `--pick` has deleted the others.

use anyhow::{anyhow, Result};
use fastvep_annotate::{splice_context, RowTranscripts};
use fastvep_cache::providers::SequenceProvider;
use fastvep_core::{Allele, Consequence, HgvsPlacement, SpliceContext, Strand};
use fastvep_genome::{Exon, Gene, Transcript};
use std::sync::Arc;

/// Reference for one contig, 1-based inclusive, `Err` past the end.
struct MemRef(String);

impl SequenceProvider for MemRef {
    fn fetch_sequence(&self, _chrom: &str, start: u64, end: u64) -> Result<Vec<u8>> {
        if start < 1 || end < start {
            return Err(anyhow!("bad range"));
        }
        let b = self.0.as_bytes();
        let s0 = (start - 1) as usize;
        if s0 >= b.len() {
            return Err(anyhow!("past contig end"));
        }
        Ok(b[s0..(end as usize).min(b.len())].to_vec())
    }
}

/// Exon 1 is 1..=20, the intron 21..=80 (opening `AGGTAGGT`), exon 2 81..=100.
fn reference() -> MemRef {
    MemRef(format!(
        "{}AGGTAGGT{}{}",
        "A".repeat(20),
        "C".repeat(52),
        "G".repeat(20)
    ))
}

fn transcript(id: &str) -> Transcript {
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
    Transcript {
        stable_id: id.into(),
        version: Some(1),
        gene: Gene {
            stable_id: "ENSG00000000001".into(),
            symbol: Some("TEST".into()),
            symbol_source: None,
            hgnc_id: None,
            biotype: "protein_coding".into(),
            chromosome: "1".into(),
            start: 1,
            end: 100,
            strand,
        },
        biotype: "protein_coding".into(),
        chromosome: "1".into(),
        start: 1,
        end: 100,
        strand,
        exons: vec![exon(1, 20, 1), exon(81, 100, 2)],
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
    }
}

fn seq(bases: &str) -> Allele {
    Allele::Sequence(bases.as_bytes().to_vec())
}

fn context(
    start: u64,
    end: u64,
    ref_allele: &Allele,
    alt_allele: &Allele,
    consequences: &[Consequence],
) -> SpliceContext {
    let reference = reference();
    splice_context(
        Some(&reference),
        "1",
        &transcript("ENST00000000001"),
        start,
        end,
        ref_allele,
        alt_allele,
        consequences,
    )
}

/// The offset is cheap and several criteria read it, so every row gets one. The
/// placement walks the reference and the exon edge is BP7's alone, so each is
/// worked out only for an allele whose terms ask for it.
#[test]
fn only_the_parts_a_criterion_reads_are_worked_out() {
    // An intronic substitution: +5, and nothing asks for the rest.
    let intronic = context(25, 25, &seq("C"), &seq("T"), &[Consequence::IntronVariant]);
    assert_eq!(
        intronic,
        SpliceContext {
            intron_offset: Some(5),
            hgvs_placement: None,
            at_exon_edge: None,
        }
    );

    // A synonymous change at an exon's last base asks for the exon edge.
    let edge = context(
        20,
        20,
        &seq("A"),
        &seq("G"),
        &[Consequence::SynonymousVariant],
    );
    assert_eq!(edge.at_exon_edge, Some(true));
    assert_eq!(edge.hgvs_placement, None);
    let mid = context(
        10,
        10,
        &seq("A"),
        &seq("G"),
        &[Consequence::SynonymousVariant],
    );
    assert_eq!(mid.at_exon_edge, Some(false));

    // A deletion with a splice term asks for the placement, and the offset is
    // where the record sits while the placement is where HGVS puts it.
    let deletion = context(
        21,
        24,
        &seq("AGGT"),
        &Allele::Deletion,
        &[Consequence::SpliceDonorVariant],
    );
    assert_eq!(deletion.intron_offset, Some(1));
    assert_eq!(
        deletion.hgvs_placement,
        Some(HgvsPlacement::Intronic { offset: 5 })
    );
    assert_eq!(deletion.at_exon_edge, None);

    // An insertion with a frameshift term, on the exon's last base.
    let insertion = context(
        21,
        20,
        &Allele::Deletion,
        &seq("AGGT"),
        &[Consequence::FrameshiftVariant],
    );
    assert_eq!(
        insertion.hgvs_placement,
        Some(HgvsPlacement::Intronic { offset: 8 })
    );
}

/// A substitution has one way of being written, so it is never walked for a
/// placement however null its term: the placement would only repeat where it is.
#[test]
fn a_substitution_is_not_given_a_placement() {
    for terms in [
        vec![Consequence::StopGained],
        vec![Consequence::SpliceDonorVariant],
        vec![Consequence::StartLost],
    ] {
        let sub = context(21, 21, &seq("A"), &seq("T"), &terms);
        assert_eq!(sub.hgvs_placement, None, "{terms:?}");
    }
    // Equal-length replacements included.
    let mnv = context(
        21,
        22,
        &seq("AG"),
        &seq("TT"),
        &[Consequence::SpliceDonorVariant],
    );
    assert_eq!(mnv.hgvs_placement, None);
}

/// Without a reference there is no shift to take, so an indel is placed where it
/// sits - the same fallback the HGVSc takes.
#[test]
fn no_reference_places_the_change_where_it_sits() {
    let ctx = splice_context(
        None,
        "1",
        &transcript("ENST00000000001"),
        21,
        24,
        &seq("AGGT"),
        &Allele::Deletion,
        &[Consequence::SpliceDonorVariant],
    );
    assert_eq!(
        ctx.hgvs_placement,
        Some(HgvsPlacement::Intronic { offset: 1 })
    );
}

fn ids(names: &[&str]) -> Vec<Arc<str>> {
    names.iter().map(|n| Arc::<str>::from(*n)).collect()
}

/// A row is found by the identity of the id it carries, and the pick only ever
/// deletes rows, so survivors arrive in prediction order.
#[test]
fn row_transcripts_follow_the_survivors_of_a_pick() {
    let models: Vec<Transcript> = ["T1", "T2", "T3", "T4"]
        .iter()
        .map(|id| transcript(id))
        .collect();
    // Make each transcript's id the same `Arc` a predicted row would carry.
    let row_ids: Vec<Arc<str>> = models.iter().map(|t| t.stable_id.clone()).collect();
    let overlapping: Vec<&Transcript> = models.iter().collect();

    // A pick kept the second and the fourth.
    let mut rows = RowTranscripts::new(&overlapping);
    assert_eq!(&*rows.of(&row_ids[1]).unwrap().stable_id, "T2");
    assert_eq!(&*rows.of(&row_ids[3]).unwrap().stable_id, "T4");

    // No pick: every row, in order.
    let mut rows = RowTranscripts::new(&overlapping);
    for (i, id) in row_ids.iter().enumerate() {
        assert!(std::ptr::eq(rows.of(id).unwrap(), overlapping[i]));
    }

    // A row with no transcript, such as a `--sa-only` scaffold.
    assert!(RowTranscripts::new(&overlapping)
        .of(&Arc::<str>::from("-"))
        .is_none());
}

/// Two models that share an id - one GFF3 loaded twice under different labels -
/// are told apart by the identity of the id, where comparing its text would hand
/// the second row the first model.
#[test]
fn row_transcripts_tell_apart_two_models_sharing_an_id() {
    let first = transcript("ENST_DUP");
    let mut second = transcript("ENST_DUP");
    second.exons.pop(); // a different exon structure under the same id
    let overlapping: Vec<&Transcript> = vec![&first, &second];

    // Only the second survived the pick.
    let mut rows = RowTranscripts::new(&overlapping);
    let found = rows.of(&second.stable_id).unwrap();
    assert!(std::ptr::eq(found, &second));
    assert_eq!(found.exons.len(), 1);
}

/// Rows out of order, or carrying an id the predictor minted itself, still
/// resolve: the in-order walk is the fast path and not the contract.
#[test]
fn row_transcripts_fall_back_when_rows_are_out_of_order_or_not_the_same_arc() {
    let models: Vec<Transcript> = ["T1", "T2", "T3"].iter().map(|id| transcript(id)).collect();
    let overlapping: Vec<&Transcript> = models.iter().collect();
    let row_ids: Vec<Arc<str>> = models.iter().map(|t| t.stable_id.clone()).collect();

    let mut rows = RowTranscripts::new(&overlapping);
    assert_eq!(&*rows.of(&row_ids[2]).unwrap().stable_id, "T3");
    // Behind the cursor now.
    assert_eq!(&*rows.of(&row_ids[0]).unwrap().stable_id, "T1");

    // Equal text, a different `Arc`.
    let minted = ids(&["T2"]);
    assert!(std::ptr::eq(
        RowTranscripts::new(&overlapping).of(&minted[0]).unwrap(),
        overlapping[1]
    ));
}
