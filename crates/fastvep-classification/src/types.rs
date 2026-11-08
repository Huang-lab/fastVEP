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
    ///
    /// Not written when null, so it has to be allowed to be missing when read
    /// back: without the default, a result with a criterion that has no details
    /// serialised to JSON that `AcmgResult` itself could not deserialise.
    #[serde(default, skip_serializing_if = "is_null_value")]
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
    /// The same variant classified once for each mode of inheritance the gene's
    /// established diseases fall under, when there is more than one.
    ///
    /// Empty for the great majority of genes, which have one. A gene with an
    /// autosomal dominant disease and an autosomal recessive one - 124 genes in
    /// ClinGen's export of 11 October 2026 - has criteria that depend on the
    /// inheritance (PM2's frequency threshold, which arm of BS2 applies, PM3, BP2),
    /// and the call above reads its diseases together: PM2's threshold, BS2's arm and
    /// BP2's trans arm take the reading for a gene that is neither dominant only
    /// nor recessive only, which is the cautious one for each, and PM3 asks only
    /// that one of the diseases be recessive. ClinGen's guidance on genes with
    /// several disorders (January 2024) is to answer each against its own
    /// disease, and this is that answer. The call above is unchanged by it: no
    /// disorder is picked on the user's behalf.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub per_disorder: Vec<DisorderClassification>,
}

/// A classification made against one mode of inheritance of a multi-disorder
/// gene: the diseases that share it, and everything an [`AcmgResult`] carries.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DisorderClassification {
    /// ClinGen's code for the mode the diseases share (`AD`, `AR`, `XL`, `SD`,
    /// `MT`, `UD`), or `unspecified` for the entries that carry no curated mode
    /// (a curation whose code ClinGen left blank or this crate does not know).
    pub inheritance: String,
    /// The established diseases the call was made against.
    pub disorders: Vec<DisorderRef>,
    /// The classification itself, flattened so that an entry reads like the call
    /// above it.
    #[serde(flatten)]
    pub result: AcmgResult,
}

/// A disease as ClinGen curates it, named in a [`DisorderClassification`].
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DisorderRef {
    /// The disease label.
    pub label: String,
    /// The MONDO identifier, when the curation carries one.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub mondo: Option<String>,
    /// ClinGen's gene-disease validity classification of the curation.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub validity: Option<String>,
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
                (EvidenceDirection::Benign, EvidenceStrength::Strong) => counts.benign_strong += 1,
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
