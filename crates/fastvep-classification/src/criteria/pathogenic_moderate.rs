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
///
/// Approximated using ClinVar pathogenic variant density as a hotspot proxy:
/// if >=N pathogenic variants exist within ±W amino acid positions, the region
/// is considered a hotspot.
fn evaluate_pm1(
    input: &ClassificationInput,
    config: &AcmgConfig,
) -> EvidenceCriterion {
    let mut details = serde_json::Map::new();
    let window = config.pm1_hotspot_window;
    let threshold = config.pm1_hotspot_min_pathogenic;
    details.insert("hotspot_window".into(), serde_json::json!(window));
    details.insert("hotspot_threshold".into(), serde_json::json!(threshold));

    let prot_pos = match input.protein_position {
        Some(pos) => pos,
        None => {
            return EvidenceCriterion {
                code: "PM1".to_string(),
                direction: EvidenceDirection::Pathogenic,
                strength: EvidenceStrength::Moderate,
                default_strength: EvidenceStrength::Moderate,
                met: false,
                evaluated: false,
                summary: "Protein position not available".to_string(),
                details: serde_json::Value::Object(details),
            };
        }
    };

    details.insert("protein_position".into(), serde_json::json!(prot_pos));

    if let Some(ref cpd) = input.clinvar_protein {
        let low = prot_pos.saturating_sub(window);
        let high = prot_pos + window;
        let nearby_pathogenic: usize = cpd
            .protein_variants
            .iter()
            .filter(|v| v.pos >= low && v.pos <= high && v.sig.to_lowercase().contains("pathogenic"))
            .count();

        details.insert("nearby_pathogenic_count".into(), serde_json::json!(nearby_pathogenic));

        let met = nearby_pathogenic >= threshold as usize;
        let summary = if met {
            format!(
                "Mutational hotspot: {} pathogenic variants within ±{} AA of position {} (threshold: {})",
                nearby_pathogenic, window, prot_pos, threshold
            )
        } else {
            format!(
                "Not a hotspot: {} pathogenic variants within ±{} AA of position {} (threshold: {})",
                nearby_pathogenic, window, prot_pos, threshold
            )
        };

        EvidenceCriterion {
            code: "PM1".to_string(),
            direction: EvidenceDirection::Pathogenic,
            strength: EvidenceStrength::Moderate,
            default_strength: EvidenceStrength::Moderate,
            met,
            evaluated: true,
            summary,
            details: serde_json::Value::Object(details),
        }
    } else {
        EvidenceCriterion {
            code: "PM1".to_string(),
            direction: EvidenceDirection::Pathogenic,
            strength: EvidenceStrength::Moderate,
            default_strength: EvidenceStrength::Moderate,
            met: false,
            evaluated: false,
            summary: "ClinVar protein-position index not available for hotspot analysis".to_string(),
            details: serde_json::Value::Object(details),
        }
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
///
/// Can be evaluated with phased VCF data and multi-variant gene context.
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
        summary: "Requires phased VCF with compound heterozygote analysis to detect in-trans with pathogenic variant".to_string(),
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
///
/// Uses the ClinVar protein-position index to check if pathogenic variants
/// with a DIFFERENT amino acid change exist at the same protein position.
fn evaluate_pm5(
    input: &ClassificationInput,
    _config: &AcmgConfig,
) -> EvidenceCriterion {
    let is_missense = input
        .consequences
        .iter()
        .any(|c| matches!(c, Consequence::MissenseVariant));

    let mut details = serde_json::Map::new();

    if !is_missense {
        return EvidenceCriterion {
            code: "PM5".to_string(),
            direction: EvidenceDirection::Pathogenic,
            strength: EvidenceStrength::Moderate,
            default_strength: EvidenceStrength::Moderate,
            met: false,
            evaluated: true,
            summary: "Not a missense variant".to_string(),
            details: serde_json::Value::Object(details),
        };
    }

    let (prot_pos, _ref_aa, alt_aa) = match (&input.protein_position, &input.amino_acids) {
        (Some(pos), Some((r, a))) => (*pos, r.as_str(), a.as_str()),
        _ => {
            return EvidenceCriterion {
                code: "PM5".to_string(),
                direction: EvidenceDirection::Pathogenic,
                strength: EvidenceStrength::Moderate,
                default_strength: EvidenceStrength::Moderate,
                met: false,
                evaluated: false,
                summary: "Protein position or amino acid change not available".to_string(),
                details: serde_json::Value::Object(details),
            };
        }
    };

    details.insert("protein_position".into(), serde_json::json!(prot_pos));
    details.insert("alt_aa".into(), serde_json::json!(alt_aa));

    if let Some(ref cpd) = input.clinvar_protein {
        // Find pathogenic variants at same position with DIFFERENT amino acid change
        let different_aa_matches: Vec<&crate::sa_extract::ClinvarProteinVariant> = cpd
            .protein_variants
            .iter()
            .filter(|v| {
                v.pos == prot_pos
                    && v.alt_aa != alt_aa
                    && v.sig.to_lowercase().contains("pathogenic")
            })
            .collect();

        details.insert(
            "different_aa_pathogenic_count".into(),
            serde_json::json!(different_aa_matches.len()),
        );

        if !different_aa_matches.is_empty() {
            let other_aas: Vec<&str> = different_aa_matches.iter().map(|v| v.alt_aa.as_str()).collect();
            details.insert("other_pathogenic_aas".into(), serde_json::json!(other_aas));

            return EvidenceCriterion {
                code: "PM5".to_string(),
                direction: EvidenceDirection::Pathogenic,
                strength: EvidenceStrength::Moderate,
                default_strength: EvidenceStrength::Moderate,
                met: true,
                evaluated: true,
                summary: format!(
                    "Different pathogenic missense at same residue {} (other AA changes: {})",
                    prot_pos,
                    other_aas.join(", ")
                ),
                details: serde_json::Value::Object(details),
            };
        }

        EvidenceCriterion {
            code: "PM5".to_string(),
            direction: EvidenceDirection::Pathogenic,
            strength: EvidenceStrength::Moderate,
            default_strength: EvidenceStrength::Moderate,
            met: false,
            evaluated: true,
            summary: format!("No different pathogenic missense at position {}", prot_pos),
            details: serde_json::Value::Object(details),
        }
    } else {
        EvidenceCriterion {
            code: "PM5".to_string(),
            direction: EvidenceDirection::Pathogenic,
            strength: EvidenceStrength::Moderate,
            default_strength: EvidenceStrength::Moderate,
            met: false,
            evaluated: false,
            summary: "ClinVar protein-position index not available".to_string(),
            details: serde_json::Value::Object(details),
        }
    }
}

/// PM6: Assumed de novo, but without confirmation of paternity and maternity.
///
/// Can be evaluated with partial trio data (one parent only).
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
        summary: "Requires trio VCF with at least one parent to assess assumed de novo status".to_string(),
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
