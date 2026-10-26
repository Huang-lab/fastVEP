//! Where a change sits in a transcript, in the form Ensembl VEP writes it.
//!
//! `cDNA_position`, `CDS_position`, `Protein_position`, `EXON` and `INTRON` are not
//! properties of a base but of a span, and the span may begin or end where the transcript
//! has no coordinate at all. A deletion that starts in an intron and ends in the next exon
//! is `?-7` in the cDNA, one that runs past the end of the coding sequence is `601-?`, and an
//! insertion is the pair of positions it sits between. Collapsing such a span to the one
//! coordinate it has, as an earlier version of this crate did, names a change that removed
//! sixteen bases as the single base at the edge.

use std::fmt::Write;

/// One coordinate system's two ends for a change (cDNA, CDS or protein), in transcript
/// order, the way Ensembl's `format_coords` is given them.
///
/// `None` is an end that falls where the transcript has no such coordinate: an intronic
/// or flanking base for the cDNA, and those or the UTRs for the CDS and the protein.
/// For an insertion `start` is one more than `end` (it sits between them).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct TranscriptSpan {
    pub start: Option<u64>,
    pub end: Option<u64>,
}

impl TranscriptSpan {
    pub const NONE: Self = Self {
        start: None,
        end: None,
    };

    pub const fn new(start: Option<u64>, end: Option<u64>) -> Self {
        Self { start, end }
    }

    /// Neither end has a coordinate: VEP's cell is `-`.
    pub fn is_none(&self) -> bool {
        self.start.is_none() && self.end.is_none()
    }

    /// Ensembl's `format_coords` (`Bio::EnsEMBL::VEP::Utils`): `a-b` low end first,
    /// `a` when the ends agree, `a-?` or `?-b` when one is unknown. Writes nothing when
    /// neither is known, so the caller chooses between VEP's `-` and the empty CSQ field.
    pub fn write_to(&self, buf: &mut String) {
        match (self.start, self.end) {
            (Some(s), Some(e)) if s == e => {
                let _ = write!(buf, "{s}");
            }
            (Some(s), Some(e)) => {
                let _ = write!(buf, "{}-{}", s.min(e), s.max(e));
            }
            (Some(s), None) => {
                let _ = write!(buf, "{s}-?");
            }
            (None, Some(e)) => {
                let _ = write!(buf, "?-{e}");
            }
            (None, None) => {}
        }
    }

    /// The two ends as VEP's JSON output has them. It builds `*_start` and `*_end` by
    /// splitting the formatted cell on `-`, and drops a `?`, so an end with no
    /// coordinate is simply absent, and a span of one position has both.
    pub fn json_bounds(&self) -> (Option<u64>, Option<u64>) {
        match (self.start, self.end) {
            (Some(s), Some(e)) => (Some(s.min(e)), Some(s.max(e))),
            (Some(s), None) => (Some(s), None),
            (None, Some(e)) => (None, Some(e)),
            (None, None) => (None, None),
        }
    }

    /// The lowest coordinate the span has, for a caller that wants a position and not a
    /// cell.
    pub fn lowest(&self) -> Option<u64> {
        match (self.start, self.end) {
            (Some(s), Some(e)) => Some(s.min(e)),
            (Some(s), None) => Some(s),
            (None, e) => e,
        }
    }
}

/// The exons (or introns) a change overlaps, numbered in transcript order, and how many
/// the transcript has: VEP's `EXON` and `INTRON` columns, `3/10` or, for a change that
/// overlaps several, `2-3/10`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Ordinals {
    pub first: u32,
    pub last: u32,
    pub total: u32,
}

impl Ordinals {
    pub fn write_to(&self, buf: &mut String) {
        if self.first == self.last {
            let _ = write!(buf, "{}/{}", self.first, self.total);
        } else {
            let _ = write!(buf, "{}-{}/{}", self.first, self.last, self.total);
        }
    }
}

/// The five location cells of a transcript row.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct TranscriptLocation {
    pub cdna: TranscriptSpan,
    pub cds: TranscriptSpan,
    pub protein: TranscriptSpan,
    pub exon: Option<Ordinals>,
    pub intron: Option<Ordinals>,
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cell(start: Option<u64>, end: Option<u64>) -> String {
        let mut s = String::new();
        TranscriptSpan::new(start, end).write_to(&mut s);
        s
    }

    /// `format_coords` of Ensembl VEP 116 (`Utils.pm`), case by case.
    #[test]
    fn spans_are_written_as_format_coords_writes_them() {
        assert_eq!(cell(Some(5), Some(5)), "5");
        assert_eq!(cell(Some(5), Some(9)), "5-9");
        // an insertion's pair arrives start > end and is written low to high
        assert_eq!(cell(Some(602), Some(601)), "601-602");
        assert_eq!(cell(Some(7), None), "7-?");
        assert_eq!(cell(None, Some(7)), "?-7");
        assert_eq!(cell(None, None), "");
    }

    #[test]
    fn json_drops_the_unknown_end_and_sorts_the_known_ones() {
        let b = |s, e| TranscriptSpan::new(s, e).json_bounds();
        assert_eq!(b(Some(602), Some(601)), (Some(601), Some(602)));
        assert_eq!(b(Some(7), Some(7)), (Some(7), Some(7)));
        assert_eq!(b(Some(7), None), (Some(7), None));
        assert_eq!(b(None, Some(7)), (None, Some(7)));
    }

    #[test]
    fn ordinals_name_a_range_only_when_there_is_one() {
        let mut s = String::new();
        Ordinals {
            first: 3,
            last: 3,
            total: 10,
        }
        .write_to(&mut s);
        assert_eq!(s, "3/10");
        s.clear();
        Ordinals {
            first: 2,
            last: 4,
            total: 10,
        }
        .write_to(&mut s);
        assert_eq!(s, "2-4/10");
    }
}
