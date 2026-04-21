use fastvep_core::Consequence;

use crate::config::AcmgConfig;
use crate::sa_extract::ClassificationInput;
use crate::types::{EvidenceCriterion, EvidenceDirection, EvidenceStrength};

/// Evaluate all pathogenic moderate criteria: PM1, PM2, PM3, PM4, PM5, PM6.
pub fn evaluate_all(
    input: &ClassificationInput,
    config: &AcmgConfig,
) -> Vec<EvidenceCriterion> {
    vec![
        evaluate_pm1(input, config),
        evaluate_pm2(input, config),
        evaluate_pm3(input, config),
        evaluate_pm4(input, config),
        evaluate_pm5(input, config),
        evaluate_pm6(input, config),
    ]
}

/// PM1: Located in a mutational hot spot and/or critical functional domain.
fn evaluate_pm1(
    _input: &ClassificationInput,
    _config: &AcmgConfig,
) -> EvidenceCriterion {
    EvidenceCriterion {
        code: "PM1".to_string(),
        direction: EvidenceDirection::Pathogenic,
        strength: EvidenceStrength::Moderate,
        default_strength: EvidenceStrength::Moderate,
        met: false,
        evaluated: false,
        summary: "Not evaluated: requires UniProt/InterPro functional domain data".to_string(),
        details: serde_json::Value::Null,
    }
}

/// PM2: Absent from controls (or at extremely low frequency if recessive).
///
/// Per ClinGen SVI recommendation, PM2 is downgraded to Supporting strength by default.
fn evaluate_pm2(
    input: &ClassificationInput,
    config: &AcmgConfig,
) -> EvidenceCriterion {
    let threshold = config.effective_pm2_threshold(input.gene_symbol.as_deref());
    let strength = if config.pm2_downgrade_to_supporting {
        EvidenceStrength::Supporting
    } else {
        EvidenceStrength::Moderate
    };

    let code = if config.pm2_downgrade_to_supporting {
        "PM2_Supporting".to_string()
    } else {
        "PM2".to_string()
    };

    let mut details = serde_json::Map::new();
    details.insert("af_threshold".into(), serde_json::json!(threshold));

    let (met, summary) = if let Some(ref gnomad) = input.gnomad {
        let af = gnomad.all_af.unwrap_or(0.0);
        details.insert("gnomad_allAf".into(), serde_json::json!(af));

        if af <= threshold {
            (
                true,
                format!(
                    "Rare in gnomAD (AF={:.6}, threshold={:.6})",
                    af, threshold
                ),
            )
        } else {
            (
                false,
                format!(
                    "Not rare in gnomAD (AF={:.6}, threshold={:.6})",
                    af, threshold
                ),
            )
        }
    } else {
        // Absent from gnomAD entirely
        details.insert("gnomad_allAf".into(), serde_json::Value::Null);
        (true, "Absent from gnomAD".to_string())
    };

    EvidenceCriterion {
        code,
        direction: EvidenceDirection::Pathogenic,
        strength,
        default_strength: EvidenceStrength::Moderate,
        met,
        evaluated: true,
        summary,
        details: serde_json::Value::Object(details),
    }
}

/// PM3: For recessive disorders, detected in trans with a pathogenic variant.
fn evaluate_pm3(
    _input: &ClassificationInput,
    _config: &AcmgConfig,
) -> EvidenceCriterion {
    EvidenceCriterion {
        code: "PM3".to_string(),
        direction: EvidenceDirection::Pathogenic,
        strength: EvidenceStrength::Moderate,
        default_strength: EvidenceStrength::Moderate,
        met: false,
        evaluated: false,
        summary: "Not evaluated: requires phasing data from trio/family analysis".to_string(),
        details: serde_json::Value::Null,
    }
}

/// PM4: Protein length changes due to in-frame deletions/insertions in non-repeat region,
/// or stop-loss variants.
fn evaluate_pm4(
    input: &ClassificationInput,
    _config: &AcmgConfig,
) -> EvidenceCriterion {
    let is_length_change = input.consequences.iter().any(|c| {
        matches!(
            c,
            Consequence::InframeInsertion | Consequence::InframeDeletion | Consequence::StopLost
        )
    });

    let mut details = serde_json::Map::new();
    if is_length_change {
        let types: Vec<&str> = input
            .consequences
            .iter()
            .filter(|c| {
                matches!(
                    c,
                    Consequence::InframeInsertion
                        | Consequence::InframeDeletion
                        | Consequence::StopLost
                )
            })
            .map(|c| c.so_term())
            .collect();
        details.insert("consequence_types".into(), serde_json::json!(types));
    }

    let summary = if is_length_change {
        "Protein length-changing variant (in-frame indel or stop-loss)".to_string()
    } else {
        "Not a protein length-changing variant".to_string()
    };

    EvidenceCriterion {
        code: "PM4".to_string(),
        direction: EvidenceDirection::Pathogenic,
        strength: EvidenceStrength::Moderate,
        default_strength: EvidenceStrength::Moderate,
        met: is_length_change,
        evaluated: true,
        summary,
        details: serde_json::Value::Object(details),
    }
}

/// PM5: Novel missense change at an amino acid residue where a different pathogenic
/// missense change has been seen before.
fn evaluate_pm5(
    _input: &ClassificationInput,
    _config: &AcmgConfig,
) -> EvidenceCriterion {
    EvidenceCriterion {
        code: "PM5".to_string(),
        direction: EvidenceDirection::Pathogenic,
        strength: EvidenceStrength::Moderate,
        default_strength: EvidenceStrength::Moderate,
        met: false,
        evaluated: false,
        summary: "Not evaluated: requires protein-position ClinVar index to check for different pathogenic missense at same residue".to_string(),
        details: serde_json::Value::Null,
    }
}

/// PM6: Assumed de novo, but without confirmation of paternity and maternity.
fn evaluate_pm6(
    _input: &ClassificationInput,
    _config: &AcmgConfig,
) -> EvidenceCriterion {
    EvidenceCriterion {
        code: "PM6".to_string(),
        direction: EvidenceDirection::Pathogenic,
        strength: EvidenceStrength::Moderate,
        default_strength: EvidenceStrength::Moderate,
        met: false,
        evaluated: false,
        summary: "Not evaluated: requires family/trio analysis data".to_string(),
        details: serde_json::Value::Null,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sa_extract::GnomadData;
    use fastvep_core::Impact;

    fn make_input(
        consequences: Vec<Consequence>,
        gnomad: Option<GnomadData>,
    ) -> ClassificationInput {
        ClassificationInput {
            consequences,
            impact: Impact::Moderate,
            gene_symbol: Some("TEST".to_string()),
            is_canonical: true,
            amino_acids: None,
            protein_position: None,
            gnomad,
            clinvar: None,
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
    fn test_pm2_absent_from_gnomad() {
        let input = make_input(vec![Consequence::MissenseVariant], None);
        let result = evaluate_pm2(&input, &AcmgConfig::default());
        assert!(result.met);
        assert_eq!(result.strength, EvidenceStrength::Supporting); // Downgraded per SVI
        assert_eq!(result.code, "PM2_Supporting");
    }

    #[test]
    fn test_pm2_rare_in_gnomad() {
        let input = make_input(
            vec![Consequence::MissenseVariant],
            Some(GnomadData {
                all_af: Some(0.00005),
                ..Default::default()
            }),
        );
        let result = evaluate_pm2(&input, &AcmgConfig::default());
        assert!(result.met);
    }

    #[test]
    fn test_pm2_common_in_gnomad() {
        let input = make_input(
            vec![Consequence::MissenseVariant],
            Some(GnomadData {
                all_af: Some(0.01),
                ..Default::default()
            }),
        );
        let result = evaluate_pm2(&input, &AcmgConfig::default());
        assert!(!result.met);
    }

    #[test]
    fn test_pm2_not_downgraded() {
        let mut config = AcmgConfig::default();
        config.pm2_downgrade_to_supporting = false;
        let input = make_input(vec![Consequence::MissenseVariant], None);
        let result = evaluate_pm2(&input, &config);
        assert!(result.met);
        assert_eq!(result.strength, EvidenceStrength::Moderate);
        assert_eq!(result.code, "PM2");
    }

    #[test]
    fn test_pm4_inframe_deletion() {
        let input = make_input(vec![Consequence::InframeDeletion], None);
        let result = evaluate_pm4(&input, &AcmgConfig::default());
        assert!(result.met);
    }

    #[test]
    fn test_pm4_stop_lost() {
        let input = make_input(vec![Consequence::StopLost], None);
        let result = evaluate_pm4(&input, &AcmgConfig::default());
        assert!(result.met);
    }

    #[test]
    fn test_pm4_missense_not_met() {
        let input = make_input(vec![Consequence::MissenseVariant], None);
        let result = evaluate_pm4(&input, &AcmgConfig::default());
        assert!(!result.met);
    }
}
