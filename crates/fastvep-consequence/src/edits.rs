//! The questions Ensembl VEP 116 asks of the *edited* transcript sequence, answered.
//!
//! Most of `VariationEffect.pm` compares two peptides, but five predicates do not trust a
//! peptide and edit the transcript instead: `_ins_del_start_altered`, `_inv_start_altered`,
//! `_snp_start_altered`, `_ins_del_stop_altered` and its "consider insertion length" twin
//! `_ins_del_stop_altered_cil`, with the two `_overlaps_stop_codon` tests that guard them. They
//! splice the allele into the sequence at the variant's cDNA or CDS position, as Ensembl's
//! mapper reports it (so an insertion is the pair `(c + 1, c)` and replaces nothing), and read one
//! codon of the result. Release 116 changed what the stop-side ones are asked and when
//! (`stop_retained` consults the `_cil` pair whenever the alternate peptide holds an `X`), so
//! the answers are kept together, here, and only computed for a change that overlaps the
//! initiator or the terminator: that is a handful of variants in a million, which is why the
//! coordinate tests come first and nothing is allocated before one of them holds.
//!
//! Every test below is a transcription, with the guards in the order Ensembl has them, and the
//! guards that look pointless are the point: `_ins_del_stop_altered_cil` returns "not altered"
//! when the variant has no CDS coordinate, so a deletion that starts in the 3' UTR and covers
//! the terminator is `stop_retained_variant` in 116.3, not `stop_lost`.

use std::borrow::Cow;

use fastvep_core::{Strand, TranscriptLocation};
use fastvep_genome::{CodonTable, Transcript};

/// A change, as the edit tests need it.
pub(crate) struct Edit<'a> {
    pub transcript: &'a Transcript,
    pub table: &'a CodonTable,
    /// Ensembl's cDNA and CDS positions of the change (`TranscriptVariation::cdna_start` ..):
    /// `?` where an end is not in an exon or not coding, the pair `(c + 1, c)` for an insertion.
    pub location: &'a TranscriptLocation,
    /// The genomic span as Ensembl holds it: `start = end + 1` for an insertion.
    pub start: u64,
    pub end: u64,
    /// The alternate allele as the VCF gives it, in genomic orientation; empty for `-`.
    pub alt: &'a [u8],
    /// Bases the reference allele covers and the alternate supplies (`ref_length`, `alt_length`).
    pub ref_len: usize,
    pub alt_len: usize,
}

/// What the edited sequence says. Every field is false for a change that does not reach the
/// codon in question, and the whole struct is the default when the transcript has no sequence.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub(crate) struct EditTests {
    /// `_overlaps_start_codon`
    pub overlaps_initiator: bool,
    pub ins_del_start_altered: bool,
    pub inv_start_altered: bool,
    pub snp_start_altered: bool,
    /// `_overlaps_stop_codon`, by cDNA position
    pub overlaps_stop: bool,
    /// `_overlaps_stop_codon_cil`, by genomic position with an insertion's length added
    pub overlaps_stop_cil: bool,
    pub ins_del_stop_altered: bool,
    pub ins_del_stop_altered_cil: bool,
}

/// The complement of a base, ambiguity codes included, as Ensembl's `reverse_comp` has it
/// (`tr/acgtrymkswhbvdnACGTRYMKSWHBVDN/tgcayrkmswdvbhnTGCAYRKMSWDVBHN/`): `R` is `Y` on the other
/// strand, and an `N` stays one.
pub(crate) fn complement(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        b'a' => b't',
        b't' => b'a',
        b'c' => b'g',
        b'g' => b'c',
        b'R' => b'Y',
        b'Y' => b'R',
        b'K' => b'M',
        b'M' => b'K',
        b'B' => b'V',
        b'V' => b'B',
        b'D' => b'H',
        b'H' => b'D',
        b'r' => b'y',
        b'y' => b'r',
        b'k' => b'm',
        b'm' => b'k',
        b'b' => b'v',
        b'v' => b'b',
        b'd' => b'h',
        b'h' => b'd',
        other => other,
    }
}

/// Alternate bases in the transcript's orientation (`feature_seq`). On the reverse strand the VCF
/// gives them in genomic order, so they are reversed *as well as* complemented: complementing in
/// place reported HPS4 c.1060_1061delTCinsAG as `stop_gained` when it is synonymous. Borrowed on
/// the forward strand, where they are the allele as given.
pub(crate) fn oriented(bases: &[u8], strand: Strand) -> Cow<'_, [u8]> {
    match strand {
        Strand::Forward => Cow::Borrowed(bases),
        Strand::Reverse => Cow::Owned(bases.iter().rev().map(|&b| complement(b)).collect()),
    }
}

fn overlap(a_start: i64, a_end: i64, b_start: i64, b_end: i64) -> bool {
    a_end >= b_start && a_start <= b_end
}

/// `substr($s, $offset, $length) = $insert`, on a copy.
fn splice(seq: &[u8], offset: usize, remove: usize, insert: &[u8]) -> Vec<u8> {
    let offset = offset.min(seq.len());
    let end = (offset + remove).min(seq.len());
    let mut out = Vec::with_capacity(seq.len() + insert.len());
    out.extend_from_slice(&seq[..offset]);
    out.extend_from_slice(insert);
    out.extend_from_slice(&seq[end..]);
    out
}

impl EditTests {
    pub(crate) fn compute(e: &Edit<'_>) -> EditTests {
        let t = e.transcript;
        let (Some(cdna_cs), Some(cdna_ce)) = (t.cdna_coding_start, t.cdna_coding_end) else {
            return EditTests::default();
        };
        let changes_length = e.ref_len != e.alt_len;
        let cdna = e.location.cdna;
        let (cdna_start, cdna_end) = match (cdna.start, cdna.end) {
            (Some(s), Some(en)) if s > 0 && en > 0 => (Some(s as i64), Some(en as i64)),
            _ => (None, None),
        };

        // The coordinate half of the three `_overlaps_*` tests, which is all most changes need:
        // the flags of the transcript are looked at only once one of them holds.
        //
        // `_overlaps_start_codon`: the cDNA span against the three bases at the coding start.
        let at_initiator = matches!((cdna_start, cdna_end), (Some(s), Some(en))
            if overlap(s, en, cdna_cs as i64, cdna_cs as i64 + 2));
        // `_overlaps_stop_codon`, by cDNA position
        let at_terminator = matches!((cdna_start, cdna_end), (Some(s), Some(en))
            if overlap(s, en, cdna_ce as i64 - 2, cdna_ce as i64));
        // `_overlaps_stop_codon_cil`: by genomic position, with the length of an inserted
        // sequence added to the start in the direction the transcript runs, so an insertion
        // counts as overlapping only if the whole of it still fits in the codon
        let at_terminator_cil = {
            let (mut v_start, v_end) = (e.start as i64, e.end as i64);
            if v_end < v_start
                && !e.alt.is_empty()
                && e.alt
                    .iter()
                    .all(|b| matches!(b, b'A' | b'C' | b'G' | b'T' | b'N'))
            {
                let len = e.alt.len() as i64;
                v_start += if t.strand == Strand::Forward {
                    len
                } else {
                    -len
                };
            }
            match (t.coding_region_start, t.coding_region_end, t.strand) {
                (Some(_), Some(cre), Strand::Forward) => {
                    overlap(v_start, v_end, cre as i64 - 2, cre as i64)
                }
                (Some(crs), Some(_), Strand::Reverse) => {
                    overlap(v_start, v_end, crs as i64, crs as i64 + 2)
                }
                _ => false,
            }
        };
        if !(at_initiator || at_terminator || at_terminator_cil) {
            return EditTests::default();
        }

        let flagged = |flag: &str| t.flags.iter().any(|f| f == flag);
        // A `cds_start_NF` transcript, or one whose CDS begins part-way through a codon, does not
        // claim to carry an initiator. The phase is the same fact as the flag in a different
        // field, and without it every length-changing variant reaching that end would be
        // `start_lost`, because the three bases at the coding start are then the tail of a codon
        // that starts before them and never read as ATG.
        let overlaps_initiator =
            at_initiator && !flagged("cds_start_NF") && t.codon_table_start_phase == 0;
        let cds_end_known = !flagged("cds_end_NF");
        let overlaps_stop = at_terminator && cds_end_known;
        let overlaps_stop_cil = at_terminator_cil && cds_end_known;
        let mut tests = EditTests {
            overlaps_initiator,
            overlaps_stop,
            overlaps_stop_cil,
            ..EditTests::default()
        };
        if !(overlaps_initiator || overlaps_stop || overlaps_stop_cil) {
            return tests;
        }

        let Some(spliced) = t.spliced_seq.as_deref().map(str::as_bytes) else {
            return tests;
        };
        let Some(translateable) = t.translateable_seq.as_deref().map(str::as_bytes) else {
            return tests;
        };
        let (Some(cdna_start), Some(cdna_end)) = (cdna_start, cdna_end) else {
            return tests;
        };
        let replaced = (cdna_end - cdna_start + 1).max(0) as usize;
        let alt = oriented(e.alt, t.strand);
        let alt = &*alt;
        let atg = |seq: &[u8], at: usize| seq.get(at..at + 3) == Some(&b"ATG"[..]);

        if overlaps_initiator {
            // the 5' UTR and the coding sequence, as Ensembl builds them: `_five_prime_utr` then
            // `_translateable_seq`
            let utr_len = cdna_cs.saturating_sub(1) as usize;
            let utr_and_translateable = &spliced[..(cdna_ce as usize).min(spliced.len())];
            if utr_and_translateable.len() >= utr_len {
                let edited = splice(
                    utr_and_translateable,
                    (cdna_start - 1) as usize,
                    replaced,
                    alt,
                );
                let tail_equals_translateable = edited.len() >= translateable.len()
                    && edited[edited.len() - translateable.len()..] == *translateable;

                // `_ins_del_start_altered`: an edit that leaves the 5' UTR alone and puts an ATG
                // at the old coding start is not one, and neither is one that merely brings the
                // whole coding sequence along behind inserted bases (atg -> aCGAtg).
                //
                // Ensembl asks the first of these only `if ($utr)`, which is a guard against a
                // missing object and has the effect that a transcript with no 5' UTR is never
                // told that its ATG is still there: an insertion inside the start codon that
                // keeps it (`atg` -> `atGg`) is `start_lost` on the transcripts of DYRK1A that
                // begin at the initiator and `start_retained_variant` on the ones that have a UTR,
                // for one edit and one codon. An empty UTR is unchanged like any other.
                if changes_length {
                    tests.ins_del_start_altered = if edited.len() >= utr_len
                        && edited[..utr_len] == utr_and_translateable[..utr_len]
                        && atg(&edited, utr_len)
                    {
                        false
                    } else if edited.len() < translateable.len() {
                        true
                    } else {
                        !tail_equals_translateable
                    };
                }

                // `_inv_start_altered`, which despite the name judges every kind of change: is the
                // codon that now sits where the initiator did an ATG? It needs a 5' UTR and gives
                // up when the edit reaches past the end of the sequence it edits.
                tests.inv_start_altered = utr_len > 0
                    && (cdna_end as usize) <= utr_and_translateable.len()
                    && !atg(&edited, utr_len);

                // `_snp_start_altered`, for a change that keeps the length
                if !changes_length {
                    let from = edited.len().saturating_sub(translateable.len());
                    tests.snp_start_altered = !atg(&edited, from);
                }
            }
        }

        // `_ins_del_stop_altered` and its `_cil` twin edit the coding sequence and the 3' UTR at
        // the CDS position, and read the codon that now sits where the terminator did. They give
        // up, and say "not altered", for a variant with no CDS start.
        let Some(cds_start) = e.location.cds.start.filter(|s| *s > 0) else {
            return tests;
        };
        if changes_length && translateable.len() >= 3 {
            let edited_stop = |overlaps: bool| -> bool {
                if !overlaps {
                    return false;
                }
                let mut cds_and_utr = translateable.to_vec();
                cds_and_utr.extend_from_slice(&spliced[(cdna_ce as usize).min(spliced.len())..]);
                let edited = splice(&cds_and_utr, (cds_start - 1) as usize, replaced, alt);
                if edited.len() < translateable.len() {
                    return true;
                }
                let at = translateable.len() - 3;
                let codon = [edited[at], edited[at + 1], edited[at + 2]];
                e.table.translate(&codon.map(|b| b.to_ascii_uppercase())) != b'*'
            };
            tests.ins_del_stop_altered = edited_stop(overlaps_stop);
            tests.ins_del_stop_altered_cil = edited_stop(overlaps_stop_cil);
        }
        tests
    }
}
