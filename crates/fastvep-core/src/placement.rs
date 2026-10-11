//! Where HGVS places an indel, against the exons and introns of a transcript.

/// The position HGVS gives a change, read against one transcript.
///
/// An indel in a repeat can be written at several positions, and every one of
/// them produces the same edited sequence. HGVS names the most 3'. A VCF record
/// anchors the change wherever its caller left-aligned it, and a consequence
/// term is decided from that, so the two can disagree about whether the change
/// touches an exon at all: MSH6 `2:47806358 G>GGTAT` is `frameshift_variant`
/// because its record anchors on the exon's last base, and `c.3801+3_3801+6dup`,
/// four intronic bases, in HGVS.
///
/// Criteria that ask "does this alter coding sequence, or a splice
/// dinucleotide?" are asking about the sequence, not about how the record was
/// written, and this is how they find out. The consequence terms themselves stay
/// where they are, because they are what Ensembl VEP reports and the output is
/// compared against it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HgvsPlacement {
    /// Wholly inside one intron. `offset` is the signed distance from the exon
    /// boundary of the nearest base the change covers: `+N` counted from the
    /// donor, `-N` from the acceptor. An insertion is the two bases it sits
    /// between, so one past `+2` leaves the dinucleotide on either side of it
    /// untouched.
    Intronic { offset: i64 },
    /// Wholly inside one exon. `in_frame` says the net change in length is a
    /// multiple of three; whether the change also creates a stop codon is what
    /// the consequence terms are for, and this does not say.
    Exonic { in_frame: bool },
    /// Anything else: it reaches both an exon and an intron, crosses an exon
    /// between two introns, covers a whole intron, or does not touch the
    /// transcript. Nothing can be concluded from it, which is a different thing
    /// from a placement not having been worked out.
    Spans,
}

impl HgvsPlacement {
    /// Whether the change lies in an intron beyond the canonical `±1`/`±2`
    /// dinucleotides, so that it alters no exon base and no canonical splice
    /// dinucleotide whichever way the record happened to be written.
    pub fn is_intronic_beyond_canonical(&self) -> bool {
        matches!(self, Self::Intronic { offset } if offset.abs() > 2)
    }
}

/// Where one allele sits against the exons and introns of one transcript, as the
/// ACMG criteria need it.
///
/// Worked out once per reported (transcript, allele) row, after `--pick` has
/// chosen which rows those are, and by nothing but a run that classifies:
/// finding any of it walks the transcript's exons, and the placement walks the
/// reference as well. `Default` is "not worked out", which is a different thing
/// from every field's own answer to it.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct SpliceContext {
    /// Signed distance in bp from the nearest exon boundary to the nearest
    /// intronic base the change covers: `+N` counted from the donor, `-N` from
    /// the acceptor, `None` for a change that reaches no intronic base.
    ///
    /// Read off the transcript at the position the record sits, which is where
    /// the consequence terms are decided, so the two agree: a change this calls
    /// `+8` does not arrive carrying `splice_donor_region_variant`. It is not
    /// read off the HGVSc, which is a *display* form and 3'-shifted: a deletion
    /// of an acceptor's own `G` is written `c.1686-1del` where the shift has
    /// nowhere to go and `c.1686del` where the first exonic base repeats it, the
    /// same variant with one spelling carrying no offset at all.
    pub intron_offset: Option<i64>,
    /// Where HGVS places the change, for the alleles whose consequence terms a
    /// criterion would read it for. `None` is "not asked", and is not
    /// [`HgvsPlacement::Spans`], which is an answer.
    pub hgvs_placement: Option<HgvsPlacement>,
    /// Whether the change touches the first base of an exon or one of its last
    /// three, the positions Walker 2023 excludes from BP7. `None` is "not
    /// asked": only a synonymous change is.
    pub at_exon_edge: Option<bool>,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn only_an_intronic_placement_past_two_clears_the_canonical_dinucleotide() {
        for offset in [-9, -3, 3, 8] {
            assert!(HgvsPlacement::Intronic { offset }.is_intronic_beyond_canonical());
        }
        for offset in [-2, -1, 1, 2] {
            assert!(!HgvsPlacement::Intronic { offset }.is_intronic_beyond_canonical());
        }
        assert!(!HgvsPlacement::Exonic { in_frame: true }.is_intronic_beyond_canonical());
        assert!(!HgvsPlacement::Spans.is_intronic_beyond_canonical());
    }
}
