//! The HGVSp of one allele: what Ensembl's `hgvs_protein` is given, worked out from a
//! transcript and a change.
//!
//! Both per-variant loops call [`hgvsp_for_allele`], so there is one place that decides what the
//! notation is built from. They used to carry a routing chain each, which drifted: a rule added
//! to one was missing from the other until a probe found it.
//!
//! What `hgvs_protein` is built from is not the change as predicted. Ensembl moves an insertion
//! or a deletion to its 3'-most position first (`_return_3prime`), maps *that* through the
//! transcript, and reads the codon window of the moved change: a deletion the 3'-rule carries
//! over an exon boundary has no protein position and so no notation, however clearly the VCF's
//! own position is in the CDS, and an insertion carried to the far side of a codon boundary is
//! described against the residues there. The peptide-level shift that follows (HGVS counts from
//! the first residue that differs) is a second, separate step in `fastvep_hgvs::protein`.
//!
//! Three things here are not Ensembl's, because what Ensembl writes in them is not the protein
//! change of the variant (`docs/VEP_DIVERGENCE.md`):
//!
//! - **The walk is over the transcript where it left the exon.** A change in one exon that the
//!   3'-rule's walk over the genome carries out of it by way of an intron is a different change of
//!   the transcript where it ends, or none: Ensembl writes no notation, or one for another residue.
//!   [`walk_the_transcript`] walks the spliced sequence instead.
//! - **A change that loses the initiator is `p.Met1?`**, as HGVS has it, where Ensembl writes a
//!   range, two residue names, or nothing ([`start_loss`]).
//! - **A row with no term for the coding sequence has no protein notation** ([`is_only_a_utr_variant`]),
//!   and neither has a change that leaves the whole coding sequence behind it, which has no protein
//!   change to name (`coding_sequence_kept`).

use std::borrow::Cow;
use std::cell::OnceCell;

use fastvep_cache::providers::SequenceProvider;
use fastvep_consequence::ConsequencePredictor;
use fastvep_core::{Allele, Consequence, Strand, TranscriptLocation, TranscriptSpan};
use fastvep_genome::Transcript;
use fastvep_hgvs::{hgvs_protein, EditedTranslation, ProteinChange};

use crate::hgvs_normalize::{is_shiftable_indel, three_prime_shift_genomic, SHIFT_WINDOW};

/// One allele against one transcript, as the per-variant loops hold it.
pub struct HgvspRequest<'a> {
    pub predictor: &'a ConsequencePredictor,
    pub transcript: &'a Transcript,
    /// The reference the 3'-rule reads. Without it nothing is shifted, which is the answer
    /// for every change that has no repeat to shift over, and the right one when there is no
    /// reference to ask.
    pub provider: Option<&'a dyn SequenceProvider>,
    pub chrom: &'a str,
    /// Ensembl's own `start` and `end` of the change: the first and last base it touches, and
    /// for an insertion `start = end + 1`, the two bases it sits between.
    pub start: u64,
    pub end: u64,
    pub ref_allele: &'a Allele,
    pub alt_allele: &'a Allele,
    /// Where the change as it stands sits in the transcript.
    pub location: &'a TranscriptLocation,
    /// The consequences of the change, all of them, before any is dropped for output.
    pub consequences: &'a [Consequence],
    /// The residue windows of the change as it stands, which are the `Amino_acids` pair
    /// (`-` for none). `None` when it has no window.
    pub amino_acids: Option<(&'a str, &'a str)>,
    /// The change leaves the whole coding sequence behind it, so the protein is the protein
    /// (`AlleleConsequenceResult::coding_sequence_kept`).
    pub coding_sequence_kept: bool,
}

/// The edited translation, built the first time a notation asks for it.
struct Edited<'a> {
    predictor: &'a ConsequencePredictor,
    transcript: &'a Transcript,
    alt: &'a Allele,
    cds: fastvep_core::TranscriptSpan,
    translation: OnceCell<Option<String>>,
}

impl EditedTranslation for Edited<'_> {
    fn alt_translation(&self) -> Option<&str> {
        self.translation
            .get_or_init(|| {
                self.predictor
                    .edited_translation(self.transcript, self.alt, self.cds)
            })
            .as_deref()
    }
}

/// How far Ensembl's `_genomic_shift` moves the change along the transcript, in bases.
///
/// Only a pure insertion or a pure deletion is shifted, and only by the length of the repeat it
/// sits in, read from a window of the chromosome a thousand bases either side
/// (`SHIFT_WINDOW`) whatever the transcript is.
fn shift_length(r: &HgvspRequest<'_>) -> u64 {
    let Some(provider) = r.provider else { return 0 };
    if !is_shiftable_indel(r.ref_allele, r.alt_allele) {
        return 0;
    }
    let (s, _) = three_prime_shift_genomic(
        provider,
        r.chrom,
        r.start,
        r.end,
        r.ref_allele,
        r.alt_allele,
        r.transcript.strand,
        r.start.saturating_sub(SHIFT_WINDOW),
        r.end.saturating_add(SHIFT_WINDOW),
    );
    match r.transcript.strand {
        Strand::Forward => s.saturating_sub(r.start),
        Strand::Reverse => r.start.saturating_sub(s),
    }
}

/// The index of the exon of `transcript` that holds the genomic base `pos`.
fn exon_holding(transcript: &Transcript, pos: u64) -> Option<usize> {
    transcript
        .exons
        .iter()
        .position(|e| e.start <= pos && pos <= e.end)
}

/// Whether the 3'-rule's walk over the genome, `shift` bases from the change, went out of the exon
/// the change is in by way of an intron.
///
/// The genome and the transcript are the same sequence for as long as the walk stays in the exon,
/// and they are not once it has crossed an intron: the bases it read in the intron are not in the
/// transcript, the base after it that did not repeat may not be the one that follows in the
/// transcript, and the string it rotated was turned once for each base of the genome it passed and
/// not of the transcript. What the walk ends at is then a different change of the transcript, and
/// a different protein, from the one the VCF has (a 4-base insertion at the end of an exon that
/// the intron repeats is `c.24+1_24+2ins`, and Ensembl writes no HGVSp for it; one that crosses the
/// whole of a short intron is described against the residue after the one it changes, and with
/// `Ter?` where the stop is there to be counted).
///
/// A walk that went past the end of the transcript is not one that crossed an intron: that is
/// the other half of what Ensembl writes nothing for, and stays as it is.
fn crossed_an_intron(r: &HgvspRequest<'_>, shift: u64) -> bool {
    let transcript = r.transcript;
    let Some(home) = home_exon(r) else {
        return false;
    };
    // the base the walk stopped on, the first that did not repeat
    let next = match transcript.strand {
        Strand::Forward => r.end + shift + 1,
        Strand::Reverse => r.start.saturating_sub(shift + 1),
    };
    if next < transcript.start || next > transcript.end {
        return false;
    }
    exon_holding(transcript, next) != Some(home)
}

/// The exon the change as the VCF has it is in: the one that holds all of it, both bases of an
/// insertion's pair included, or for an insertion at the edge of an exon the one it is beside, with
/// the exon's own last base (or first) for one flank and the intron's for the other. Ensembl calls
/// such an insertion coding (`c.1617+1dup` is `frameshift_variant`) and the transcript has a place
/// for it: between the exon's last base and the next exon's first.
fn home_exon(r: &HgvspRequest<'_>) -> Option<usize> {
    let transcript = r.transcript;
    let insertion = r.start == r.end + 1;
    let (lo, hi) = (r.start.min(r.end), r.start.max(r.end));
    match (exon_holding(transcript, lo), exon_holding(transcript, hi)) {
        (Some(a), Some(b)) if a == b => Some(a),
        (Some(a), None) if insertion && lo == transcript.exons[a].end => Some(a),
        (None, Some(b)) if insertion && hi == transcript.exons[b].start => Some(b),
        _ => None,
    }
}

/// The 3'-rule over the sequence of the transcript instead of the genome: where in the transcript
/// the change ends up, and how many of its bases it passed to get there. `None` where the
/// transcript has no spliced sequence or the change has no place in it.
///
/// Only for a change [`crossed_an_intron`] says the genome walk is not the transcript's.
fn walk_the_transcript(r: &HgvspRequest<'_>) -> Option<(TranscriptLocation, u64)> {
    let transcript = r.transcript;
    let spliced = transcript.spliced_seq.as_deref()?.as_bytes();
    let (lo, hi) = (r.start.min(r.end), r.start.max(r.end));
    let forward = transcript.strand == Strand::Forward;
    if r.start == r.end + 1 {
        // The insertion sits after the base `c` of the transcript: after its 5' flank, or, where
        // that one is the intron's, in front of its 3' flank (`home_exon`)
        let (five, three) = if forward { (lo, hi) } else { (hi, lo) };
        let c = match (
            transcript.genomic_to_cdna(five),
            transcript.genomic_to_cdna(three),
        ) {
            (Some(c), _) => c,
            (None, Some(c)) => c.checked_sub(1)?,
            (None, None) => return None,
        };
        let inserted = if forward {
            r.alt_allele.as_bytes().to_vec()
        } else {
            fastvep_genome::codon::reverse_complement(r.alt_allele.as_bytes())
        };
        if inserted.is_empty() {
            return None;
        }
        // slides over a base that is the next base of the inserted string, which turns it
        let mut passed = 0u64;
        while spliced
            .get((c + passed) as usize)
            .is_some_and(|b| *b == inserted[passed as usize % inserted.len()].to_ascii_uppercase())
        {
            passed += 1;
        }
        let at = c + passed;
        Some((
            transcript.location_in_cdna(TranscriptSpan::new(Some(at + 1), Some(at)), true),
            passed,
        ))
    } else {
        // the ends in the order of the transcript, the 5' one first
        let (c_lo, c_hi) = (
            transcript.genomic_to_cdna(lo)?,
            transcript.genomic_to_cdna(hi)?,
        );
        let (first, last) = if forward { (c_lo, c_hi) } else { (c_hi, c_lo) };
        // slides over a base that is the first base of what is deleted
        let mut passed = 0u64;
        while spliced
            .get((last + passed) as usize)
            .zip(spliced.get((first + passed - 1) as usize))
            .is_some_and(|(next, front)| next == front)
        {
            passed += 1;
        }
        Some((
            transcript.location_in_cdna(
                TranscriptSpan::new(Some(first + passed), Some(last + passed)),
                false,
            ),
            passed,
        ))
    }
}

/// An inserted string as it reads after the change moved `shift` bases 3': one rotation per
/// base, in the direction of travel, which is left for a transcript on the forward strand and
/// right for one on the reverse.
fn rotate(alt: &Allele, shift: u64, strand: Strand) -> Allele {
    match alt {
        Allele::Sequence(bases) if !bases.is_empty() => {
            let mut rotated = bases.clone();
            let k = (shift % bases.len() as u64) as usize;
            match strand {
                Strand::Forward => rotated.rotate_left(k),
                Strand::Reverse => rotated.rotate_right(k),
            }
            Allele::Sequence(rotated)
        }
        other => other.clone(),
    }
}

/// The HGVSp of one allele, `None` where Ensembl writes none.
///
/// Except for a change that loses the initiator, which has one whatever Ensembl does: HGVS names
/// it `p.Met1?`, and Ensembl writes nothing for the ones it cannot place in the protein (a
/// deletion that begins in the 5' UTR and runs into the start codon has no protein start, `?-1`).
/// That is 138 of the 2,762 start-loss rows of the ClinVar 2-star+ sample, and 11,217 of the 96,350
/// of complete transcripts in VEP 116.3's Ensembl cache run over 45,916 changes at the start and
/// stop codons and exon edges of chr21 and chr22.
pub fn hgvsp_for_allele(r: &HgvspRequest<'_>) -> Option<String> {
    hgvsp_as_ensembl_writes_it(r).or_else(|| start_loss(r))
}

/// `p.Met1?` for a change the consequence code calls `start_lost`, in the protein it is a change of.
fn start_loss(r: &HgvspRequest<'_>) -> Option<String> {
    if !r.consequences.contains(&Consequence::StartLost)
        || !(r.alt_allele.is_plain_bases() && r.ref_allele.is_plain_bases())
    {
        return None;
    }
    let first = *r.transcript.peptide.as_deref()?.as_bytes().first()?;
    Some(format!(
        "{}:p.{}1?",
        r.transcript.versioned_protein_id()?,
        fastvep_genome::codon::aa_one_to_three(first)
    ))
}

/// Whether the change is in the UTR as it was written and has no term about the coding sequence:
/// it touches no residue, however a 3'-rule moves it.
fn is_only_a_utr_variant(consequences: &[Consequence]) -> bool {
    use Consequence::*;
    let in_utr = consequences
        .iter()
        .any(|c| matches!(c, FivePrimeUtrVariant | ThreePrimeUtrVariant));
    let in_cds = consequences.iter().any(|c| {
        matches!(
            c,
            StopGained
                | FrameshiftVariant
                | StopLost
                | StartLost
                | InframeInsertion
                | InframeDeletion
                | MissenseVariant
                | ProteinAlteringVariant
                | IncompleteTerminalCodonVariant
                | StartRetainedVariant
                | StopRetainedVariant
                | SynonymousVariant
                | CodingSequenceVariant
        )
    });
    in_utr && !in_cds
}

fn hgvsp_as_ensembl_writes_it(r: &HgvspRequest<'_>) -> Option<String> {
    let transcript = r.transcript;
    // `return undef if variation_feature_seq =~ m/[^ACGT\-]/` for the allele, and the same
    // for the reference, whose peptide is the other half of the notation: with an `N` for
    // the reference base the reference allele has no peptide and `hgvs_protein` returns
    // undef at `return undef unless $hgvs_notation->{ref}`.
    if !(r.alt_allele.is_plain_bases() && r.ref_allele.is_plain_bases()) {
        return None;
    }
    if !transcript.pre_coding(r.start, r.end) {
        return None;
    }
    // Ensembl asks `pre_coding` of the two bases an insertion sits between and moves the change
    // with the 3'-rule afterwards, so an insertion between the last base of the 5' UTR and the
    // initiator that repeats the codon beside it has a protein change written for it
    // (`p.Met1AsnfsTer?`, `p.Met1_Ala2insMet`, `p.Met1dup`) on a row whose consequence is
    // `5_prime_UTR_variant` and nothing else: the consequence is of the change where the VCF puts
    // it, the notation of the change where the 3'-rule moved it, and the row says both. A row with
    // no term for the coding sequence gets no protein notation. 3,422 of the 350,819 UTR-only rows
    // of VEP 116.3's Ensembl cache run (the complete transcripts of chr21 and chr22, 45,916 changes
    // at the start and stop codons and exon edges) and 386 of the 52,001 of its RefSeq cache run;
    // none of the 196,229 of the ClinVar 2-star+ sample or of the 1,500 of HG002.
    if is_only_a_utr_variant(r.consequences) {
        return None;
    }
    // The edited sequence still has the whole coding sequence in it and no ATG where the initiator
    // was, so the first ATG is the one the transcript had. Ensembl takes the notation from the
    // window of codons at the old coding start and writes `p.Met1?` or a frameshift for a protein
    // that has no change (the consequence step leaves `frameshift_variant` and `start_lost` off
    // such a change, `EditTests::cds_kept`).
    if r.coding_sequence_kept {
        return None;
    }
    let reference_translation = transcript.peptide.as_deref()?;

    // The change is moved first, and what follows is about the change in its new place.
    let mut shift = shift_length(r);
    let (location, alt): (TranscriptLocation, Cow<'_, Allele>) = if shift == 0 {
        (*r.location, Cow::Borrowed(r.alt_allele))
    } else {
        let walked = if crossed_an_intron(r, shift) {
            Some(walk_the_transcript(r)?)
        } else {
            None
        };
        match walked {
            Some((location, passed)) => {
                shift = passed;
                // not moved at all, in the transcript, is not moved
                let alt = if passed == 0 {
                    Cow::Borrowed(r.alt_allele)
                } else {
                    Cow::Owned(rotate(r.alt_allele, passed, transcript.strand))
                };
                (if passed == 0 { *r.location } else { location }, alt)
            }
            None => {
                let offset = |p: u64| match transcript.strand {
                    Strand::Forward => p + shift,
                    Strand::Reverse => p.saturating_sub(shift),
                };
                (
                    transcript.location(offset(r.start), offset(r.end)),
                    Cow::Owned(rotate(r.alt_allele, shift, transcript.strand)),
                )
            }
        }
    };
    // `translation_start && translation_end`: both ends of the moved change have to be protein
    let (start, end) = (location.protein.start?, location.protein.end?);
    if start == 0 || end == 0 {
        return None;
    }

    // A change that was not moved has the window the consequence code already built
    let (ref_peptide, alt_peptide): (Cow<'_, str>, Cow<'_, str>) = match r.amino_acids {
        Some((ref_aas, alt_aas)) if shift == 0 => (Cow::Borrowed(ref_aas), Cow::Borrowed(alt_aas)),
        _ => {
            let window = r.predictor.protein_window(transcript, &alt, &location)?;
            (
                Cow::Owned(window.ref_peptide),
                Cow::Owned(window.alt_peptide),
            )
        }
    };
    let has = |c: Consequence| r.consequences.contains(&c);
    // `frameshift` is the one predicate `hgvs_protein` evaluates afresh, and it reads the CDS
    // span of the change *after* it was moved, where `start_lost` and `stop_lost` are cached
    // from the first evaluation. So a single-base insertion in front of the initiator is not a
    // `frameshift_variant` of its own, and once the 3'-rule has carried it onto the first codon
    // it is one, with a notation to match: `p.Met1AsnfsTer?`.
    let cds_len = (location.cds.end? as i64 - location.cds.start? as i64 + 1).max(0) as usize;
    let moves_frame =
        !has(Consequence::IncompleteTerminalCodonVariant) && alt.len().abs_diff(cds_len) % 3 != 0;
    // It is not Ensembl's predicate, though. That one answers no for a change that starts on the
    // terminator or keeps it (release 116 settles which *term* such a change earns, and the frame
    // is not what it is settling), and the notation is then built from the window alone. That is
    // sound where the window reaches a stop, since nothing past it is protein, and where it does
    // not it ends part-way through a codon whose `X` Ensembl reads as a terminator:
    // `p.Ter583=` beside `stop_lost`, `p.Ter201delinsTyrTer` for a change whose next residue
    // nobody has looked at. The frame did move there, and the edited translation says what comes
    // of it.
    // The terminator is asked about where the change was before it moved: `_get_peptide_alleles`
    // is cached from the first evaluation.
    let exempt = has(Consequence::StopRetainedVariant)
        || r.amino_acids
            .is_some_and(|(ref_aas, _)| ref_aas.starts_with('*'));
    let frameshift = moves_frame && (!exempt || !alt_peptide.contains('*'));

    let edited = Edited {
        predictor: r.predictor,
        transcript,
        alt: &alt,
        cds: location.cds,
        translation: OnceCell::new(),
    };
    // Last, because it allocates and most (allele, transcript) pairs never get this far
    let ref_name = transcript.versioned_protein_id()?;
    hgvs_protein(
        &ProteinChange {
            ref_name: &ref_name,
            translation_start: start,
            translation_end: end,
            ref_peptide: &ref_peptide,
            alt_peptide: &alt_peptide,
            reference_translation,
            frameshift,
            start_lost: has(Consequence::StartLost),
            stop_lost: has(Consequence::StopLost),
        },
        &edited,
    )
}
