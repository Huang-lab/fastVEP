use crate::config::AcmgConfig;
use crate::sa_extract::ClassificationInput;
use crate::types::{EvidenceCriterion, EvidenceDirection, EvidenceStrength};

/// Evaluate all pathogenic strong criteria: PS1, PS2, PS3, PS4.
pub fn evaluate_all(
    input: &ClassificationInput,
    config: &AcmgConfig,
) -> Vec<EvidenceCriterion> {
    vec![
        evaluate_ps1(input, config),
        evaluate_ps2(input, config),
        evaluate_ps3(input, config),
        evaluate_ps4(input, config),
    ]
}

/// PS1: Same amino acid change as a previously established pathogenic variant.
///
/// True PS1 requires knowing that a different nucleotide change at the same position
/// with the same amino acid change is pathogenic. This needs a protein-position ClinVar
/// index not currently available, so we mark as not evaluated.
fn evaluate_ps1(
    _input: &ClassificationInput,
    _config: &AcmgConfig,
) -> EvidenceCriterion {
    EvidenceCriterion {
        code: "PS1".to_string(),
        direction: EvidenceDirection::Pathogenic,
        strength: EvidenceStrength::Strong,
        default_strength: EvidenceStrength::Strong,
        met: false,
        evaluated: false,
        summary: "Not evaluated: requires protein-position ClinVar index to check if a different nucleotide change at the same position with the same AA change is pathogenic".to_string(),
        details: serde_json::Value::Null,
    }
}

/// PS2: De novo (both maternity and paternity confirmed) in a patient with the disease.
fn evaluate_ps2(
    _input: &ClassificationInput,
    _config: &AcmgConfig,
) -> EvidenceCriterion {
    EvidenceCriterion {
        code: "PS2".to_string(),
        direction: EvidenceDirection::Pathogenic,
        strength: EvidenceStrength::Strong,
        default_strength: EvidenceStrength::Strong,
        met: false,
        evaluated: false,
        summary: "Not evaluated: requires de novo status from trio/family analysis".to_string(),
        details: serde_json::Value::Null,
    }
}

/// PS3: Well-established in vitro or in vivo functional studies.
fn evaluate_ps3(
    _input: &ClassificationInput,
    _config: &AcmgConfig,
) -> EvidenceCriterion {
    EvidenceCriterion {
        code: "PS3".to_string(),
        direction: EvidenceDirection::Pathogenic,
        strength: EvidenceStrength::Strong,
        default_strength: EvidenceStrength::Strong,
        met: false,
        evaluated: false,
        summary: "Not evaluated: requires functional study data".to_string(),
        details: serde_json::Value::Null,
    }
}

/// PS4: Prevalence of the variant in affected individuals is significantly increased
/// compared with controls.
///
/// Approximated using ClinVar pathogenic with high review confidence.
fn evaluate_ps4(
    input: &ClassificationInput,
    _config: &AcmgConfig,
) -> EvidenceCriterion {
    let mut details = serde_json::Map::new();

    let (met, summary) = if let Some(ref clinvar) = input.clinvar {
        let stars = clinvar.review_stars();
        let is_pathogenic = clinvar.has_pathogenic();
        details.insert("clinvar_pathogenic".into(), serde_json::json!(is_pathogenic));
        details.insert("review_stars".into(), serde_json::json!(stars));

        if is_pathogenic && stars >= 3 {
            (
                true,
                format!(
                    "ClinVar pathogenic with {}-star review (expert panel or practice guideline)",
                    stars
                ),
            )
        } else {
            (
                false,
                format!(
                    "ClinVar significance: {:?}, review: {}-star (needs >=3 stars for PS4)",
                    clinvar.significance.as_deref().unwrap_or(&[]),
                    stars
                ),
            )
        }
    } else {
        (false, "No ClinVar data available".to_string())
    };

    EvidenceCriterion {
        code: "PS4".to_string(),
        direction: EvidenceDirection::Pathogenic,
        strength: EvidenceStrength::Strong,
        default_strength: EvidenceStrength::Strong,
        met,
        evaluated: input.clinvar.is_some(),
        summary,
        details: serde_json::Value::Object(details),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sa_extract::ClinvarData;
    use fastvep_core::{Consequence, Impact};

    fn make_input(clinvar: Option<ClinvarData>) -> ClassificationInput {
        ClassificationInput {
            consequences: vec![Consequence::MissenseVariant],
            impact: Impact::Moderate,
            gene_symbol: Some("TP53".to_string()),
            is_canonical: true,
            amino_acids: Some(("R".to_string(), "H".to_string())),
            protein_position: Some(175),
            gnomad: None,
            clinvar,
            revel: None,
            splice_ai: None,
            dbnsfp: None,
            phylop: None,
            gerp: None,
            gene_constraints: None,
            omim: None,
        }
    }

    #[test]
    fn test_ps4_expert_panel_pathogenic() {
        let input = make_input(Some(ClinvarData {
            significance: Some(vec!["Pathogenic".to_string()]),
            review_status: Some("reviewed_by_expert_panel".to_string()),
            ..Default::default()
        }));
        let result = evaluate_ps4(&input, &AcmgConfig::default());
        assert!(result.met);
        assert!(result.evaluated);
    }

    #[test]
    fn test_ps4_single_submitter_not_enough() {
        let input = make_input(Some(ClinvarData {
            significance: Some(vec!["Pathogenic".to_string()]),
            review_status: Some("criteria_provided,_single_submitter".to_string()),
            ..Default::default()
        }));
        let result = evaluate_ps4(&input, &AcmgConfig::default());
        assert!(!result.met);
    }

    #[test]
    fn test_ps1_not_evaluated() {
        let input = make_input(None);
        let result = evaluate_ps1(&input, &AcmgConfig::default());
        assert!(!result.evaluated);
        assert!(!result.met);
    }
}
