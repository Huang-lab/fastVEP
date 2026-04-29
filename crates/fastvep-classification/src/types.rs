use serde::{Deserialize, Serialize};

/// Strength level of an ACMG-AMP evidence criterion.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize, Deserialize)]
pub enum EvidenceStrength {
    Supporting,
    Moderate,
    Strong,
    VeryStrong,
    Standalone,
}

impl EvidenceStrength {
    pub fn as_str(&self) -> &'static str {
        match self {
            Self::Supporting => "Supporting",
            Self::Moderate => "Moderate",
            Self::Strong => "Strong",
            Self::VeryStrong => "Very_Strong",
            Self::Standalone => "Stand_Alone",
        }
    }
}

/// Direction of evidence: pathogenic or benign.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum EvidenceDirection {
    Pathogenic,
    Benign,
}

/// A single ACMG-AMP evidence criterion with its evaluation result.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EvidenceCriterion {
    /// Criterion code, e.g. "PVS1", "PM2_Supporting", "PP3"
    pub code: String,
    /// Whether this is pathogenic or benign evidence
    pub direction: EvidenceDirection,
    /// Strength level (possibly modified from default per ClinGen SVI)
    pub strength: EvidenceStrength,
    /// Default strength per Richards et al. 2015
    pub default_strength: EvidenceStrength,
    /// Whether this criterion was triggered (evidence met)
    pub met: bool,
    /// Whether this criterion could be evaluated from available data
    pub evaluated: bool,
    /// Human-readable explanation of why met/not met
    pub summary: String,
    /// Key data points used in evaluation (for transparency)
    #[serde(skip_serializing_if = "is_null_value")]
    pub details: serde_json::Value,
}

fn is_null_value(v: &serde_json::Value) -> bool {
    v.is_null()
}

/// Five-tier ACMG-AMP variant classification.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum AcmgClassification {
    Pathogenic,
    LikelyPathogenic,
    UncertainSignificance,
    LikelyBenign,
    Benign,
}

impl AcmgClassification {
    pub fn as_str(&self) -> &'static str {
        match self {
            Self::Pathogenic => "Pathogenic",
            Self::LikelyPathogenic => "Likely_pathogenic",
            Self::UncertainSignificance => "Uncertain_significance",
            Self::LikelyBenign => "Likely_benign",
            Self::Benign => "Benign",
        }
    }

    pub fn shorthand(&self) -> &'static str {
        match self {
            Self::Pathogenic => "P",
            Self::LikelyPathogenic => "LP",
            Self::UncertainSignificance => "VUS",
            Self::LikelyBenign => "LB",
            Self::Benign => "B",
        }
    }
}

/// Full ACMG-AMP classification result for one allele-transcript pair.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AcmgResult {
    /// Final 5-tier classification
    pub classification: AcmgClassification,
    /// Shorthand code (P, LP, VUS, LB, B)
    pub shorthand: String,
    /// All criteria that were evaluated
    pub criteria: Vec<EvidenceCriterion>,
    /// Which combination rule was triggered (e.g. "PVS1+PM")
    #[serde(skip_serializing_if = "Option::is_none")]
    pub triggered_rule: Option<String>,
    /// Summary counts of met criteria by direction/strength
    pub counts: EvidenceCounts,
}

/// Counts of met evidence criteria by direction and strength.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct EvidenceCounts {
    pub pathogenic_very_strong: u8,
    pub pathogenic_strong: u8,
    pub pathogenic_moderate: u8,
    pub pathogenic_supporting: u8,
    pub benign_standalone: u8,
    pub benign_strong: u8,
    pub benign_supporting: u8,
}

impl EvidenceCounts {
    /// Count met criteria from a list of evaluated criteria.
    pub fn from_criteria(criteria: &[EvidenceCriterion]) -> Self {
        let mut counts = Self::default();
        for c in criteria {
            if !c.met {
                continue;
            }
            match (c.direction, c.strength) {
                (EvidenceDirection::Pathogenic, EvidenceStrength::VeryStrong) => {
                    counts.pathogenic_very_strong += 1
                }
                (EvidenceDirection::Pathogenic, EvidenceStrength::Strong) => {
                    counts.pathogenic_strong += 1
                }
                (EvidenceDirection::Pathogenic, EvidenceStrength::Moderate) => {
                    counts.pathogenic_moderate += 1
                }
                (EvidenceDirection::Pathogenic, EvidenceStrength::Supporting) => {
                    counts.pathogenic_supporting += 1
                }
                (EvidenceDirection::Benign, EvidenceStrength::Standalone) => {
                    counts.benign_standalone += 1
                }
                (EvidenceDirection::Benign, EvidenceStrength::Strong) => {
                    counts.benign_strong += 1
                }
                (EvidenceDirection::Benign, EvidenceStrength::Supporting) => {
                    counts.benign_supporting += 1
                }
                // Benign Moderate is not standard in ACMG but we count it as strong
                (EvidenceDirection::Benign, EvidenceStrength::Moderate) => {
                    counts.benign_strong += 1
                }
                // Benign Very Strong (e.g. BP4_Very_Strong from REVEL ≤ 0.003 per
                // Pejaver 2022) reaches the Benign classification on its own under
                // the Tavtigian Bayesian framework. Counting it as 2 BS triggers
                // the existing ≥2 BS → Benign rule without introducing a new slot.
                (EvidenceDirection::Benign, EvidenceStrength::VeryStrong) => {
                    counts.benign_strong += 2
                }
                // Pathogenic + Standalone is not defined in ACMG/AMP.
                (EvidenceDirection::Pathogenic, EvidenceStrength::Standalone) => {}
            }
        }
        counts
    }

    pub fn has_pathogenic(&self) -> bool {
        self.pathogenic_very_strong > 0
            || self.pathogenic_strong > 0
            || self.pathogenic_moderate > 0
            || self.pathogenic_supporting > 0
    }

    pub fn has_benign(&self) -> bool {
        self.benign_standalone > 0 || self.benign_strong > 0 || self.benign_supporting > 0
    }
}
