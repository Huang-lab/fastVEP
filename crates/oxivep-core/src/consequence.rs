use serde::{Deserialize, Serialize};

/// Subjective impact classification of a consequence type.
/// Ordered from most to least severe.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub enum Impact {
    /// Disruptive: truncating, frameshifting, or splice-site-destroying.
    High,
    /// Non-disruptive change to protein: missense, in-frame indel.
    Moderate,
    /// Mostly harmless: synonymous, splice region.
    Low,
    /// Non-coding or distant: intronic, upstream, downstream, intergenic.
    Modifier,
}

impl Impact {
    /// Return the impact as an uppercase static string, avoiding format!("{:?}").
    pub fn as_str(self) -> &'static str {
        match self {
            Impact::High => "HIGH",
            Impact::Moderate => "MODERATE",
            Impact::Low => "LOW",
            Impact::Modifier => "MODIFIER",
        }
    }
}

/// Sequence Ontology consequence terms used by Ensembl VEP.
///
/// Variants are ordered by severity rank (lower rank = more severe).
/// Ranks and SO terms match Ensembl VEP release 115.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum Consequence {
    TranscriptAblation,
    SpliceAcceptorVariant,
    SpliceDonorVariant,
    StopGained,
    FrameshiftVariant,
    StopLost,
    StartLost,
    TranscriptAmplification,
    FeatureElongation,
    FeatureTruncation,
    InframeInsertion,
    InframeDeletion,
    MissenseVariant,
    ProteinAlteringVariant,
    SpliceRegionVariant,
    SpliceDonorFifthBaseVariant,
    SpliceDonorRegionVariant,
    SplicePolypyrimidineTractVariant,
    IncompleteTerminalCodonVariant,
    StartRetainedVariant,
    StopRetainedVariant,
    SynonymousVariant,
    CodingSequenceVariant,
    MatureMirnaVariant,
    FivePrimeUtrVariant,
    ThreePrimeUtrVariant,
    NonCodingTranscriptExonVariant,
    IntronVariant,
    NmdTranscriptVariant,
    NonCodingTranscriptVariant,
    CodingTranscriptVariant,
    UpstreamGeneVariant,
    DownstreamGeneVariant,
    TfbsAblation,
    TfbsAmplification,
    TfBindingSiteVariant,
    RegulatoryRegionAblation,
    RegulatoryRegionAmplification,
    RegulatoryRegionVariant,
    IntergenicVariant,
    SequenceVariant,
    // -- Structural variant consequence terms (Nirvana parity) --
    CopyNumberChange,
    CopyNumberIncrease,
    CopyNumberDecrease,
    ShortTandemRepeatChange,
    ShortTandemRepeatExpansion,
    ShortTandemRepeatContraction,
    /// Unidirectional gene fusion (from breakend analysis).
    UnidirectionalGeneFusion,
    /// Transcript variant (generic SV overlap with transcript).
    TranscriptVariant,
}

impl Consequence {
    /// Severity rank (lower = more severe). Matches Ensembl VEP ordering.
    pub fn rank(self) -> u32 {
        match self {
            Self::TranscriptAblation => 1,
            Self::SpliceAcceptorVariant => 2,
            Self::SpliceDonorVariant => 3,
            Self::StopGained => 4,
            Self::FrameshiftVariant => 5,
            Self::StopLost => 6,
            Self::StartLost => 7,
            Self::TranscriptAmplification => 8,
            Self::FeatureElongation => 9,
            Self::FeatureTruncation => 10,
            Self::InframeInsertion => 11,
            Self::InframeDeletion => 12,
            Self::MissenseVariant => 13,
            Self::ProteinAlteringVariant => 14,
            Self::SpliceRegionVariant => 15,
            Self::SpliceDonorFifthBaseVariant => 16,
            Self::SpliceDonorRegionVariant => 17,
            Self::SplicePolypyrimidineTractVariant => 18,
            Self::IncompleteTerminalCodonVariant => 19,
            Self::StartRetainedVariant => 20,
            Self::StopRetainedVariant => 21,
            Self::SynonymousVariant => 22,
            Self::CodingSequenceVariant => 23,
            Self::MatureMirnaVariant => 24,
            Self::FivePrimeUtrVariant => 25,
            Self::ThreePrimeUtrVariant => 26,
            Self::NonCodingTranscriptExonVariant => 27,
            Self::IntronVariant => 28,
            Self::NmdTranscriptVariant => 29,
            Self::NonCodingTranscriptVariant => 30,
            Self::CodingTranscriptVariant => 31,
            Self::UpstreamGeneVariant => 32,
            Self::DownstreamGeneVariant => 33,
            Self::TfbsAblation => 34,
            Self::TfbsAmplification => 35,
            Self::TfBindingSiteVariant => 36,
            Self::RegulatoryRegionAblation => 37,
            Self::RegulatoryRegionAmplification => 38,
            Self::RegulatoryRegionVariant => 39,
            Self::IntergenicVariant => 40,
            Self::SequenceVariant => 41,
            Self::CopyNumberChange => 42,
            Self::CopyNumberIncrease => 43,
            Self::CopyNumberDecrease => 44,
            Self::ShortTandemRepeatChange => 45,
            Self::ShortTandemRepeatExpansion => 46,
            Self::ShortTandemRepeatContraction => 47,
            Self::UnidirectionalGeneFusion => 48,
            Self::TranscriptVariant => 49,
        }
    }

    /// The impact classification for this consequence.
    pub fn impact(self) -> Impact {
        match self {
            Self::TranscriptAblation
            | Self::SpliceAcceptorVariant
            | Self::SpliceDonorVariant
            | Self::StopGained
            | Self::FrameshiftVariant
            | Self::StopLost
            | Self::StartLost
            | Self::TranscriptAmplification
            | Self::TfbsAblation
            | Self::RegulatoryRegionAblation => Impact::High,

            Self::InframeInsertion
            | Self::InframeDeletion
            | Self::MissenseVariant
            | Self::ProteinAlteringVariant
            | Self::RegulatoryRegionAmplification
            | Self::TfbsAmplification => Impact::Moderate,

            Self::SpliceRegionVariant
            | Self::SpliceDonorFifthBaseVariant
            | Self::SpliceDonorRegionVariant
            | Self::SplicePolypyrimidineTractVariant
            | Self::IncompleteTerminalCodonVariant
            | Self::StartRetainedVariant
            | Self::StopRetainedVariant
            | Self::SynonymousVariant => Impact::Low,

            Self::CodingSequenceVariant
            | Self::MatureMirnaVariant
            | Self::FivePrimeUtrVariant
            | Self::ThreePrimeUtrVariant
            | Self::NonCodingTranscriptExonVariant
            | Self::IntronVariant
            | Self::NmdTranscriptVariant
            | Self::NonCodingTranscriptVariant
            | Self::CodingTranscriptVariant
            | Self::UpstreamGeneVariant
            | Self::DownstreamGeneVariant
            | Self::TfBindingSiteVariant
            | Self::RegulatoryRegionVariant
            | Self::IntergenicVariant
            | Self::SequenceVariant
            | Self::FeatureElongation
            | Self::FeatureTruncation
            | Self::CopyNumberChange
            | Self::CopyNumberIncrease
            | Self::CopyNumberDecrease
            | Self::ShortTandemRepeatChange
            | Self::ShortTandemRepeatExpansion
            | Self::ShortTandemRepeatContraction
            | Self::UnidirectionalGeneFusion
            | Self::TranscriptVariant => Impact::Modifier,
        }
    }

    /// The SO term as used in VEP output (snake_case).
    pub fn so_term(self) -> &'static str {
        match self {
            Self::TranscriptAblation => "transcript_ablation",
            Self::SpliceAcceptorVariant => "splice_acceptor_variant",
            Self::SpliceDonorVariant => "splice_donor_variant",
            Self::StopGained => "stop_gained",
            Self::FrameshiftVariant => "frameshift_variant",
            Self::StopLost => "stop_lost",
            Self::StartLost => "start_lost",
            Self::TranscriptAmplification => "transcript_amplification",
            Self::FeatureElongation => "feature_elongation",
            Self::FeatureTruncation => "feature_truncation",
            Self::InframeInsertion => "inframe_insertion",
            Self::InframeDeletion => "inframe_deletion",
            Self::MissenseVariant => "missense_variant",
            Self::ProteinAlteringVariant => "protein_altering_variant",
            Self::SpliceRegionVariant => "splice_region_variant",
            Self::SpliceDonorFifthBaseVariant => "splice_donor_5th_base_variant",
            Self::SpliceDonorRegionVariant => "splice_donor_region_variant",
            Self::SplicePolypyrimidineTractVariant => "splice_polypyrimidine_tract_variant",
            Self::IncompleteTerminalCodonVariant => "incomplete_terminal_codon_variant",
            Self::StartRetainedVariant => "start_retained_variant",
            Self::StopRetainedVariant => "stop_retained_variant",
            Self::SynonymousVariant => "synonymous_variant",
            Self::CodingSequenceVariant => "coding_sequence_variant",
            Self::MatureMirnaVariant => "mature_miRNA_variant",
            Self::FivePrimeUtrVariant => "5_prime_UTR_variant",
            Self::ThreePrimeUtrVariant => "3_prime_UTR_variant",
            Self::NonCodingTranscriptExonVariant => "non_coding_transcript_exon_variant",
            Self::IntronVariant => "intron_variant",
            Self::NmdTranscriptVariant => "NMD_transcript_variant",
            Self::NonCodingTranscriptVariant => "non_coding_transcript_variant",
            Self::CodingTranscriptVariant => "coding_transcript_variant",
            Self::UpstreamGeneVariant => "upstream_gene_variant",
            Self::DownstreamGeneVariant => "downstream_gene_variant",
            Self::TfbsAblation => "TFBS_ablation",
            Self::TfbsAmplification => "TFBS_amplification",
            Self::TfBindingSiteVariant => "TF_binding_site_variant",
            Self::RegulatoryRegionAblation => "regulatory_region_ablation",
            Self::RegulatoryRegionAmplification => "regulatory_region_amplification",
            Self::RegulatoryRegionVariant => "regulatory_region_variant",
            Self::IntergenicVariant => "intergenic_variant",
            Self::SequenceVariant => "sequence_variant",
            Self::CopyNumberChange => "copy_number_change",
            Self::CopyNumberIncrease => "copy_number_increase",
            Self::CopyNumberDecrease => "copy_number_decrease",
            Self::ShortTandemRepeatChange => "short_tandem_repeat_change",
            Self::ShortTandemRepeatExpansion => "short_tandem_repeat_expansion",
            Self::ShortTandemRepeatContraction => "short_tandem_repeat_contraction",
            Self::UnidirectionalGeneFusion => "unidirectional_gene_fusion",
            Self::TranscriptVariant => "transcript_variant",
        }
    }

    /// Parse a consequence from its SO term string.
    pub fn from_so_term(term: &str) -> Option<Self> {
        match term {
            "transcript_ablation" => Some(Self::TranscriptAblation),
            "splice_acceptor_variant" => Some(Self::SpliceAcceptorVariant),
            "splice_donor_variant" => Some(Self::SpliceDonorVariant),
            "stop_gained" => Some(Self::StopGained),
            "frameshift_variant" => Some(Self::FrameshiftVariant),
            "stop_lost" => Some(Self::StopLost),
            "start_lost" => Some(Self::StartLost),
            "transcript_amplification" => Some(Self::TranscriptAmplification),
            "feature_elongation" => Some(Self::FeatureElongation),
            "feature_truncation" => Some(Self::FeatureTruncation),
            "inframe_insertion" => Some(Self::InframeInsertion),
            "inframe_deletion" => Some(Self::InframeDeletion),
            "missense_variant" => Some(Self::MissenseVariant),
            "protein_altering_variant" => Some(Self::ProteinAlteringVariant),
            "splice_region_variant" => Some(Self::SpliceRegionVariant),
            "splice_donor_5th_base_variant" => Some(Self::SpliceDonorFifthBaseVariant),
            "splice_donor_region_variant" => Some(Self::SpliceDonorRegionVariant),
            "splice_polypyrimidine_tract_variant" => {
                Some(Self::SplicePolypyrimidineTractVariant)
            }
            "incomplete_terminal_codon_variant" => Some(Self::IncompleteTerminalCodonVariant),
            "start_retained_variant" => Some(Self::StartRetainedVariant),
            "stop_retained_variant" => Some(Self::StopRetainedVariant),
            "synonymous_variant" => Some(Self::SynonymousVariant),
            "coding_sequence_variant" => Some(Self::CodingSequenceVariant),
            "mature_miRNA_variant" => Some(Self::MatureMirnaVariant),
            "5_prime_UTR_variant" => Some(Self::FivePrimeUtrVariant),
            "3_prime_UTR_variant" => Some(Self::ThreePrimeUtrVariant),
            "non_coding_transcript_exon_variant" => Some(Self::NonCodingTranscriptExonVariant),
            "intron_variant" => Some(Self::IntronVariant),
            "NMD_transcript_variant" => Some(Self::NmdTranscriptVariant),
            "non_coding_transcript_variant" => Some(Self::NonCodingTranscriptVariant),
            "coding_transcript_variant" => Some(Self::CodingTranscriptVariant),
            "upstream_gene_variant" => Some(Self::UpstreamGeneVariant),
            "downstream_gene_variant" => Some(Self::DownstreamGeneVariant),
            "TFBS_ablation" => Some(Self::TfbsAblation),
            "TFBS_amplification" => Some(Self::TfbsAmplification),
            "TF_binding_site_variant" => Some(Self::TfBindingSiteVariant),
            "regulatory_region_ablation" => Some(Self::RegulatoryRegionAblation),
            "regulatory_region_amplification" => Some(Self::RegulatoryRegionAmplification),
            "regulatory_region_variant" => Some(Self::RegulatoryRegionVariant),
            "intergenic_variant" => Some(Self::IntergenicVariant),
            "sequence_variant" => Some(Self::SequenceVariant),
            "copy_number_change" => Some(Self::CopyNumberChange),
            "copy_number_increase" => Some(Self::CopyNumberIncrease),
            "copy_number_decrease" => Some(Self::CopyNumberDecrease),
            "short_tandem_repeat_change" => Some(Self::ShortTandemRepeatChange),
            "short_tandem_repeat_expansion" => Some(Self::ShortTandemRepeatExpansion),
            "short_tandem_repeat_contraction" => Some(Self::ShortTandemRepeatContraction),
            "unidirectional_gene_fusion" => Some(Self::UnidirectionalGeneFusion),
            "transcript_variant" => Some(Self::TranscriptVariant),
            _ => None,
        }
    }

    /// Return the most severe consequence from a list.
    pub fn most_severe(consequences: &[Consequence]) -> Option<Consequence> {
        consequences.iter().copied().min_by_key(|c| c.rank())
    }

    /// Return the most severe impact from a list of consequences.
    pub fn worst_impact(consequences: &[Consequence]) -> Option<Impact> {
        consequences.iter().map(|c| c.impact()).min()
    }
}

impl PartialOrd for Consequence {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Consequence {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.rank().cmp(&other.rank())
    }
}

impl std::fmt::Display for Consequence {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.so_term())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_consequence_ranking() {
        assert!(Consequence::TranscriptAblation < Consequence::IntergenicVariant);
        assert!(Consequence::StopGained < Consequence::MissenseVariant);
        assert!(Consequence::MissenseVariant < Consequence::SynonymousVariant);
        assert!(Consequence::SynonymousVariant < Consequence::IntronVariant);
    }

    #[test]
    fn test_impact_ordering() {
        assert!(Impact::High < Impact::Moderate);
        assert!(Impact::Moderate < Impact::Low);
        assert!(Impact::Low < Impact::Modifier);
    }

    #[test]
    fn test_consequence_impact() {
        assert_eq!(Consequence::StopGained.impact(), Impact::High);
        assert_eq!(Consequence::MissenseVariant.impact(), Impact::Moderate);
        assert_eq!(Consequence::SynonymousVariant.impact(), Impact::Low);
        assert_eq!(Consequence::IntronVariant.impact(), Impact::Modifier);
    }

    #[test]
    fn test_so_term_round_trip() {
        let consequences = [
            Consequence::TranscriptAblation,
            Consequence::MissenseVariant,
            Consequence::SpliceRegionVariant,
            Consequence::IntergenicVariant,
            Consequence::FivePrimeUtrVariant,
            Consequence::NmdTranscriptVariant,
            Consequence::TfbsAblation,
        ];
        for c in consequences {
            let term = c.so_term();
            let parsed = Consequence::from_so_term(term).unwrap();
            assert_eq!(c, parsed, "round-trip failed for {term}");
        }
    }

    #[test]
    fn test_most_severe() {
        let cs = vec![
            Consequence::IntronVariant,
            Consequence::MissenseVariant,
            Consequence::SynonymousVariant,
        ];
        assert_eq!(Consequence::most_severe(&cs), Some(Consequence::MissenseVariant));
    }
}
