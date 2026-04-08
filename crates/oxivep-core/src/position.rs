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
    pub fn from_str(s: &str) -> Self {
        match s {
            "-" => Allele::Deletion,
            "*" => Allele::Missing,
            _ if s.starts_with('<') && s.ends_with('>') => {
                Allele::Symbolic(s.to_string())
            }
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
        assert_eq!(
            Allele::from_str("ACGT"),
            Allele::Sequence(b"ACGT".to_vec())
        );
    }

    #[test]
    fn test_strand_round_trip() {
        assert_eq!(Strand::from_int(1), Strand::Forward);
        assert_eq!(Strand::from_int(-1), Strand::Reverse);
        assert_eq!(Strand::Forward.as_int(), 1);
        assert_eq!(Strand::Reverse.as_int(), -1);
    }
}
