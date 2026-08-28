//! The HGVS 3'-rule holds for intronic indels, on both strands.
//!
//! An insertion inside a repeat can be written at any position the repeat allows
//! and describes the same variant every time, so every spelling must produce one
//! HGVSc. fastVEP did not: the dup anchor was derived by walking the *unshifted*
//! insertion point one base at a time while the next base repeated the current
//! one, which slides through a homopolymer and nowhere else. A `TG` insertion in
//! a `TGTG…` repeat therefore never moved.
//!
//! Measured against real Ensembl VEP 115.1 (`ensemblorg/ensembl-vep:release_115.1`,
//! `--gff` mode, Ensembl 115 GRCh38 GFF3 + FASTA) over six genome-wide variants
//! written both ways: VEP returned one string on all 168 transcript rows,
//! fastVEP returned two on all 168. Genome-wide over a 20,241-variant HG002
//! sample that was 1,390 of 1,538 disagreeing HGVSc rows.
//!
//! These tests go through `hgvsc_intronic_shifted`, which is the intronic path
//! both annotation loops call, so a regression in either one fails here.

use anyhow::{anyhow, Result};
use fastvep_annotate::hgvsc_intronic_shifted;
use fastvep_cache::providers::SequenceProvider;
use fastvep_core::{Allele, Strand};
use fastvep_genome::{Exon, Gene, Transcript};

/// Reference for one contig, 1-based inclusive, `Err` past the end - the same
/// contract the FASTA readers keep.
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

/// Exon 1 at 1..=20, intron at 21..=80, exon 2 at 81..=100. The intron is a
/// `TG` repeat, so an insertion of `TG` slides its whole length.
fn reference() -> MemRef {
    MemRef(format!(
        "{}{}{}",
        "A".repeat(20),
        "TG".repeat(30),
        "C".repeat(20)
    ))
}

fn transcript(strand: Strand) -> Transcript {
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

fn revcomp(bases: &[u8]) -> Vec<u8> {
    bases
        .iter()
        .rev()
        .map(|b| match b.to_ascii_uppercase() {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            other => other,
        })
        .collect()
}

/// Annotate one spelling of an intronic insertion, the way both per-variant
/// loops do: the shift reads the genomic alleles, the notation is written from
/// the transcript-oriented pair.
fn hgvsc_for_insertion(strand: Strand, vcf_pos: u64, ins: &[u8]) -> Option<String> {
    let reference = reference();
    let tr = transcript(strand);
    let genomic_alt = Allele::Sequence(ins.to_vec());
    let hgvs_alt = match strand {
        Strand::Forward => genomic_alt.clone(),
        Strand::Reverse => Allele::Sequence(revcomp(ins)),
    };
    hgvsc_intronic_shifted(
        Some(&reference),
        "1",
        &tr,
        "ENST00000000001.1",
        // Ensembl's convention for an insertion: start is one past end.
        vcf_pos + 1,
        vcf_pos,
        &Allele::Deletion,
        &genomic_alt,
        &Allele::Deletion,
        &hgvs_alt,
        Some(1),
        Some(40),
    )
}

/// Every position of the repeat is the same variant, so every spelling must
/// produce the same HGVSc. This is the property that failed on 168 of 168 rows.
#[test]
fn every_spelling_of_an_intronic_insertion_agrees_on_the_forward_strand() {
    // The repeat runs 21..=80. Inserting after position `p` carries the rotation
    // that position implies: `TG` at an even offset, `GT` at an odd one.
    let spellings: Vec<String> = (20..80)
        .map(|p| {
            let ins: &[u8] = if (p - 20) % 2 == 0 { b"TG" } else { b"GT" };
            hgvsc_for_insertion(Strand::Forward, p, ins)
                .unwrap_or_else(|| panic!("no HGVSc at position {p}"))
        })
        .collect();
    let first = &spellings[0];
    assert!(
        spellings.iter().all(|s| s == first),
        "60 spellings of one variant produced {} distinct strings: {:?}",
        spellings
            .iter()
            .collect::<std::collections::HashSet<_>>()
            .len(),
        spellings.iter().collect::<std::collections::HashSet<_>>()
    );
    // The 3'-most placement on the forward strand is the end of the repeat, two
    // bases before exon 2, so the block is written from the downstream anchor.
    assert_eq!(first, "ENST00000000001.1:c.21-2_21-1dup");
}

/// 3' runs towards lower coordinates on the reverse strand, so the same repeat
/// normalises to its other end - and the tool must not simply keep whichever end
/// the VCF happened to name.
#[test]
fn every_spelling_of_an_intronic_insertion_agrees_on_the_reverse_strand() {
    let spellings: Vec<String> = (20..80)
        .map(|p| {
            let ins: &[u8] = if (p - 20) % 2 == 0 { b"TG" } else { b"GT" };
            hgvsc_for_insertion(Strand::Reverse, p, ins)
                .unwrap_or_else(|| panic!("no HGVSc at position {p}"))
        })
        .collect();
    let first = &spellings[0];
    assert!(
        spellings.iter().all(|s| s == first),
        "60 spellings of one variant produced {:?}",
        spellings.iter().collect::<std::collections::HashSet<_>>()
    );
    assert!(
        first.ends_with("dup"),
        "an insertion that duplicates the repeat it sits in is a dup, got {first}"
    );
}

/// The two strands must not land on the same genomic block: the rule is
/// direction-sensitive, and a shift that ignored strand would pick the same end
/// of the repeat for both.
///
/// The two `c.` strings coincide, which is not a coincidence and not agreement:
/// `c.21-1` counts one base 5' of exon 2's first base in *transcript* order, so
/// it is genomic 80 on the forward strand and genomic 21 on the reverse. The
/// same name denotes mirrored blocks, which is exactly why the genomic span is
/// the thing worth asserting.
#[test]
fn the_two_strands_normalise_to_opposite_ends_of_the_repeat() {
    let reference = reference();
    let fwd = fastvep_annotate::intronic_dup_span(&reference, "1", 81, b"TG", 60, Strand::Forward);
    let rev = fastvep_annotate::intronic_dup_span(&reference, "1", 21, b"TG", 0, Strand::Reverse);
    assert_eq!(fwd, Some((79, 80)), "forward 3' is the end nearer exon 2");
    assert_eq!(rev, Some((21, 22)), "reverse 3' is the end nearer exon 1");
}

/// An insertion is written over both bases it sits between, and across the
/// middle of an intron those two anchor to different exons: `+n` counts from the
/// exon before, `-m` from the one after. Inferring the second coordinate as
/// `offset + 1` produced `c.20+30_20+31ins…`, naming a base past the half the
/// donor-side offsets reach.
#[test]
fn an_insertion_at_the_intron_midpoint_is_written_from_both_exons() {
    // Intron 21..=80, so genomic 50 is c.20+30 and genomic 51 is c.21-30.
    // `AAA` matches neither side of the `TG` repeat, so nothing shifts and the
    // coordinates are the ones under test.
    let out = hgvsc_for_insertion(Strand::Forward, 50, b"AAA");
    assert_eq!(
        out.as_deref(),
        Some("ENST00000000001.1:c.20+30_21-30insAAA")
    );
}

/// A deletion follows the same rule, and shares the same shift.
#[test]
fn every_spelling_of_an_intronic_deletion_agrees() {
    let reference = reference();
    let tr = transcript(Strand::Forward);
    let spellings: Vec<String> = (21..79)
        .step_by(2)
        .map(|start| {
            let deleted: &[u8] = b"TG";
            hgvsc_intronic_shifted(
                Some(&reference),
                "1",
                &tr,
                "ENST00000000001.1",
                start,
                start + 1,
                &Allele::Sequence(deleted.to_vec()),
                &Allele::Deletion,
                &Allele::Sequence(deleted.to_vec()),
                &Allele::Deletion,
                Some(1),
                Some(40),
            )
            .unwrap_or_else(|| panic!("no HGVSc at position {start}"))
        })
        .collect();
    let first = &spellings[0];
    assert!(
        spellings.iter().all(|s| s == first),
        "spellings of one deletion disagreed: {:?}",
        spellings.iter().collect::<std::collections::HashSet<_>>()
    );
    assert_eq!(first, "ENST00000000001.1:c.21-2_21-1del");
}

/// Without a reference there is nothing to shift against, and the variant is
/// still described - at the position the VCF gave it.
#[test]
fn no_reference_still_produces_a_description() {
    let tr = transcript(Strand::Forward);
    let out = hgvsc_intronic_shifted(
        None,
        "1",
        &tr,
        "ENST00000000001.1",
        41,
        40,
        &Allele::Deletion,
        &Allele::Sequence(b"TG".to_vec()),
        &Allele::Deletion,
        &Allele::Sequence(b"TG".to_vec()),
        Some(1),
        Some(40),
    );
    assert_eq!(out.as_deref(), Some("ENST00000000001.1:c.20+20_20+21insTG"));
}

/// The inserted string rotates as the variant travels, and `hgvs_alt` is already
/// in transcript orientation - where the shift always runs 3' whichever way the
/// transcript does, so the rotation is always leftward. Rotating right on the
/// reverse strand named the right position with the wrong bases:
/// `c.21-21_21-20insGTC` where the same transcript calls the same variant
/// `insCGT` when the VCF spells it one base over.
///
/// This only shows on a *partial* shift. A shift of a whole period, or one that
/// ends on a duplication, hides it - which is why 716 of 717 reverse-strand
/// intronic `ins` rows in the HG002 sample agreed with VEP regardless.
#[test]
fn a_partial_shift_rotates_the_insert_the_same_way_on_both_strands() {
    // Deliberately imperfect repeats, so the shift stops mid-period.
    let patterns = [
        "ACGACGACTTTGCA",
        "AACAACAAGGTTCC",
        "GTGTGTTTACCAGG",
        "TTATTATTGCCAAG",
    ];
    let inserts: [&[u8]; 5] = [b"ACG", b"AAC", b"GT", b"TTA", b"ACGA"];
    let mut disagreed = Vec::new();
    for pattern in patterns {
        let seq = format!("{}{}{}", "N".repeat(20), pattern.repeat(10), "N".repeat(20));
        let reference = MemRef(seq.clone());
        for ins in inserts {
            for strand in [Strand::Forward, Strand::Reverse] {
                for base in [40u64, 61, 77] {
                    let forms = equivalent_spellings(&seq, base, ins);
                    if forms.len() < 2 {
                        continue;
                    }
                    // Every spelling must produce a string. Dropping the ones
                    // that do not would let a regression that stopped naming
                    // them pass with a set of size one.
                    let answers: std::collections::BTreeSet<String> = forms
                        .iter()
                        .map(|(pos, genomic)| {
                            let hgvs_alt = match strand {
                                Strand::Forward => Allele::Sequence(genomic.clone()),
                                Strand::Reverse => Allele::Sequence(revcomp(genomic)),
                            };
                            hgvsc_intronic_shifted(
                                Some(&reference),
                                "1",
                                &wide_transcript(strand),
                                "ENST00000000001.1",
                                pos + 1,
                                *pos,
                                &Allele::Deletion,
                                &Allele::Sequence(genomic.clone()),
                                &Allele::Deletion,
                                &hgvs_alt,
                                Some(1),
                                Some(40),
                            )
                            .unwrap_or_else(|| {
                                panic!(
                                    "no HGVSc for {pattern} ins={} {strand:?} at {pos}",
                                    String::from_utf8_lossy(genomic)
                                )
                            })
                        })
                        .collect();
                    if answers.len() > 1 {
                        disagreed.push(format!(
                            "{pattern} ins={} {strand:?} at {base}: {answers:?}",
                            String::from_utf8_lossy(ins)
                        ));
                    }
                }
            }
        }
    }
    assert!(
        disagreed.is_empty(),
        "{} configurations gave more than one answer:\n  {}",
        disagreed.len(),
        disagreed.join("\n  ")
    );
}

/// Every way the reference lets one insertion be written: slide it 3' along the
/// genome one base at a time, rotating the inserted string as it goes, for as
/// long as the sequence is unchanged.
fn equivalent_spellings(reference: &str, pos: u64, ins: &[u8]) -> Vec<(u64, Vec<u8>)> {
    let bases = reference.as_bytes();
    let len = ins.len();
    let mut out = vec![(pos, ins.to_vec())];
    let mut k = 0usize;
    while (pos as usize) + k < bases.len()
        && bases[(pos as usize) + k].eq_ignore_ascii_case(&ins[k % len])
    {
        k += 1;
        let mut rotated = ins.to_vec();
        rotated.rotate_left(k % len);
        out.push((pos + k as u64, rotated));
    }
    out
}

/// Exons far enough apart that the patterns above sit well inside the intron.
fn wide_transcript(strand: Strand) -> Transcript {
    let mut t = transcript(strand);
    t.end = 200;
    t.gene.end = 200;
    t.exons[1].start = 181;
    t.exons[1].end = 200;
    t
}
