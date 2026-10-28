use serde::{Deserialize, Serialize};

/// Classification of variant type for dispatch between small-variant and SV pipelines.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum VariantType {
    /// Single nucleotide variant.
    Snv,
    /// Insertion of one or more bases.
    Insertion,
    /// Deletion of one or more bases.
    Deletion,
    /// Combined insertion and deletion.
    Indel,
    /// Multi-nucleotide variant (substitution of >1 base).
    Mnv,
    // -- Structural variant types --
    /// Copy number variation (unspecified direction).
    CopyNumberVariation,
    /// Copy number loss.
    CopyNumberLoss,
    /// Copy number gain.
    CopyNumberGain,
    /// Tandem duplication.
    TandemDuplication,
    /// Inversion.
    Inversion,
    /// Translocation breakend.
    TranslocationBreakend,
    /// Short tandem repeat expansion/contraction.
    ShortTandemRepeatVariation,
    /// Unknown or unclassified variant.
    Unknown,
}

impl VariantType {
    /// Returns true for structural variant types that need the SV annotation pipeline.
    pub fn is_structural(self) -> bool {
        matches!(
            self,
            Self::CopyNumberVariation
                | Self::CopyNumberLoss
                | Self::CopyNumberGain
                | Self::TandemDuplication
                | Self::Inversion
                | Self::TranslocationBreakend
                | Self::ShortTandemRepeatVariation
        )
    }

    /// Returns true for small variant types handled by the standard pipeline.
    pub fn is_small(self) -> bool {
        matches!(
            self,
            Self::Snv | Self::Insertion | Self::Deletion | Self::Indel | Self::Mnv
        )
    }
}

/// Strand orientation of a genomic feature.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum Strand {
    Forward,
    Reverse,
}

impl Strand {
    pub fn from_int(val: i8) -> Self {
        if val >= 0 {
            Strand::Forward
        } else {
            Strand::Reverse
        }
    }

    pub fn as_int(self) -> i8 {
        match self {
            Strand::Forward => 1,
            Strand::Reverse => -1,
        }
    }
}

/// A position on a reference genome.
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct GenomicPosition {
    pub chromosome: String,
    /// 1-based start coordinate (Ensembl convention).
    pub start: u64,
    /// 1-based end coordinate, inclusive.
    pub end: u64,
    pub strand: Strand,
}

impl GenomicPosition {
    pub fn new(chromosome: impl Into<String>, start: u64, end: u64, strand: Strand) -> Self {
        Self {
            chromosome: chromosome.into(),
            start,
            end,
            strand,
        }
    }

    /// Length of the spanned region in bases.
    pub fn length(&self) -> u64 {
        self.end.saturating_sub(self.start) + 1
    }

    /// Check whether two positions overlap on the same chromosome.
    pub fn overlaps(&self, other: &GenomicPosition) -> bool {
        self.chromosome == other.chromosome && self.start <= other.end && other.start <= self.end
    }
}

/// Representation of a variant allele.
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum Allele {
    /// A sequence of nucleotide bases (A, C, G, T, or N).
    Sequence(Vec<u8>),
    /// A deletion represented as "-".
    Deletion,
    /// A missing/unspecified allele "*".
    Missing,
    /// A symbolic allele such as `<DEL>`, `<DUP>`, `<INV>`, `<CNV>`, `<INS>`, `<BND>`.
    Symbolic(String),
}

impl Allele {
    /// Create an allele from a string representation.
    // Inherent constructor, deliberately infallible: an unrecognised input
    // maps to a catch-all variant rather than an error, so the fallible
    // `std::str::FromStr` contract does not fit.
    #[allow(clippy::should_implement_trait)]
    pub fn from_str(s: &str) -> Self {
        match s {
            "-" => Allele::Deletion,
            "*" => Allele::Missing,
            _ if s.starts_with('<') && s.ends_with('>') => Allele::Symbolic(s.to_string()),
            _ => Allele::Sequence(s.as_bytes().to_vec()),
        }
    }

    /// Get the length of the allele in bases (0 for deletion/missing/symbolic).
    pub fn len(&self) -> usize {
        match self {
            Allele::Sequence(bases) => bases.len(),
            Allele::Deletion | Allele::Missing | Allele::Symbolic(_) => 0,
        }
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Return the sequence bytes, or an empty slice for deletion/missing/symbolic.
    pub fn as_bytes(&self) -> &[u8] {
        match self {
            Allele::Sequence(bases) => bases,
            Allele::Deletion | Allele::Missing | Allele::Symbolic(_) => &[],
        }
    }

    /// Returns true if this is a symbolic allele (e.g., `<DEL>`, `<DUP>`).
    pub fn is_symbolic(&self) -> bool {
        matches!(self, Allele::Symbolic(_))
    }

    /// Whether this allele is something Ensembl can translate: a run of `A`, `C`, `G` and `T`
    /// in either case, or the `-` of a deletion (`seq_is_unambiguous_dna`, `^[ACGT-]+$`).
    ///
    /// `N` and the other ambiguity codes are not, and neither are `*` and a symbolic allele.
    /// For an allele that fails this, Ensembl has no peptide, so no `Amino_acids`, no
    /// `HGVSp`, no `HGVSc`, and no consequence term that is about the protein: a change
    /// that moves the frame is a `frameshift_variant` and anything else in the coding
    /// sequence is `coding_sequence_variant`.
    pub fn is_plain_bases(&self) -> bool {
        match self {
            Allele::Deletion => true,
            Allele::Sequence(bases) => bases
                .iter()
                .all(|b| matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T')),
            Allele::Missing | Allele::Symbolic(_) => false,
        }
    }

    /// The length Ensembl gives this allele (`seq_length`): the number of bases, ambiguity
    /// codes included (`^[ACGTUMRWSYKVHDBXN-]+$`), 0 for `-`, and none for `*` and a
    /// symbolic allele, whose length nobody can say.
    pub fn known_length(&self) -> Option<usize> {
        match self {
            Allele::Deletion => Some(0),
            Allele::Sequence(bases)
                if bases.iter().all(|b| {
                    matches!(
                        b.to_ascii_uppercase(),
                        b'A' | b'C'
                            | b'G'
                            | b'T'
                            | b'U'
                            | b'M'
                            | b'R'
                            | b'W'
                            | b'S'
                            | b'Y'
                            | b'K'
                            | b'V'
                            | b'H'
                            | b'D'
                            | b'B'
                            | b'X'
                            | b'N'
                    )
                }) =>
            {
                Some(bases.len())
            }
            _ => None,
        }
    }
}

impl std::fmt::Display for Allele {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Allele::Sequence(bases) => {
                write!(f, "{}", std::str::from_utf8(bases).unwrap_or("?"))
            }
            Allele::Deletion => write!(f, "-"),
            Allele::Missing => write!(f, "*"),
            Allele::Symbolic(s) => write!(f, "{}", s),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_genomic_position_overlap() {
        let a = GenomicPosition::new("chr1", 100, 200, Strand::Forward);
        let b = GenomicPosition::new("chr1", 150, 250, Strand::Forward);
        let c = GenomicPosition::new("chr1", 300, 400, Strand::Forward);
        let d = GenomicPosition::new("chr2", 100, 200, Strand::Forward);

        assert!(a.overlaps(&b));
        assert!(b.overlaps(&a));
        assert!(!a.overlaps(&c));
        assert!(!a.overlaps(&d));
    }

    #[test]
    fn test_allele_from_str() {
        assert_eq!(Allele::from_str("-"), Allele::Deletion);
        assert_eq!(Allele::from_str("*"), Allele::Missing);
        assert_eq!(Allele::from_str("ACGT"), Allele::Sequence(b"ACGT".to_vec()));
    }

    #[test]
    fn only_a_run_of_acgt_or_a_deletion_is_plain_bases() {
        let is = |s: &str| Allele::from_str(s).is_plain_bases();
        assert!(is("ACGT") && is("acgt") && is("-") && is("T"));
        assert!(!is("N") && !is("ACGN") && !is("R") && !is("*") && !is("<DEL>"));
    }

    #[test]
    fn a_length_is_known_for_bases_with_ambiguity_codes_and_not_for_star_or_symbolic() {
        let len = |s: &str| Allele::from_str(s).known_length();
        assert_eq!(len("ACGTN"), Some(5));
        assert_eq!(len("-"), Some(0));
        assert_eq!(len("NN"), Some(2));
        assert_eq!(len("*"), None);
        assert_eq!(len("<*>"), None);
        assert_eq!(len("<INS>"), None);
    }

    #[test]
    fn test_strand_round_trip() {
        assert_eq!(Strand::from_int(1), Strand::Forward);
        assert_eq!(Strand::from_int(-1), Strand::Reverse);
        assert_eq!(Strand::Forward.as_int(), 1);
        assert_eq!(Strand::Reverse.as_int(), -1);
    }
}
