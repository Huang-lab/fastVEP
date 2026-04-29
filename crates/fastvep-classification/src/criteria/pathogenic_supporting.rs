use fastvep_core::Consequence;

use crate::config::AcmgConfig;
use crate::sa_extract::ClassificationInput;
use crate::types::{EvidenceCriterion, EvidenceDirection, EvidenceStrength};

/// Evaluate all pathogenic supporting criteria: PP1, PP2, PP3, PP4, PP5.
pub fn evaluate_all(
    input: &ClassificationInput,
    config: &AcmgConfig,
) -> Vec<EvidenceCriterion> {
    let mut criteria = vec![
        evaluate_pp1(input, config),
        evaluate_pp2(input, config),
        evaluate_pp3(input, config),
        evaluate_pp4(input, config),
    ];
    if config.use_pp5_bp6 {
        criteria.push(evaluate_pp5(input, config));
    }
    criteria
}

/// PP1: Co-segregation with disease in multiple affected family members.
fn evaluate_pp1(
    _input: &ClassificationInput,
    _config: &AcmgConfig,
) -> EvidenceCriterion {
    EvidenceCriterion {
        code: "PP1".to_string(),
        direction: EvidenceDirection::Pathogenic,
        strength: EvidenceStrength::Supporting,
        default_strength: EvidenceStrength::Supporting,
        met: false,
        evaluated: false,
        summary: "Requires multi-generation pedigree with affection status to assess co-segregation".to_string(),
        details: serde_json::Value::Null,
    }
}

/// PP2: Missense variant in a gene that has a low rate of benign missense variation
/// and in which missense variants are a common mechanism of disease.
fn evaluate_pp2(
    input: &ClassificationInput,
    config: &AcmgConfig,
) -> EvidenceCriterion {
    let is_missense = input
        .consequences
        .iter()
        .any(|c| matches!(c, Consequence::MissenseVariant));

    let mut details = serde_json::Map::new();
    details.insert("is_missense".into(), serde_json::json!(is_missense));

    let (met, evaluated, summary) = if !is_missense {
        (false, true, "Not a missense variant".to_string())
    } else if let Some(ref gc) = input.gene_constraints {
        if let Some(mis_z) = gc.mis_z {
            details.insert("misZ".into(), serde_json::json!(mis_z));
            details.insert(
                "threshold".into(),
                serde_json::json!(config.pp2_misz_threshold),
            );
            if mis_z > config.pp2_misz_threshold {
                (
                    true,
                    true,
                    format!(
                        "Missense in gene with high missense constraint (misZ={:.2}, threshold={:.2})",
                        mis_z, config.pp2_misz_threshold
                    ),
                )
            } else {
                (
                    false,
                    true,
                    format!(
                        "Gene missense constraint not significant (misZ={:.2}, threshold={:.2})",
                        mis_z, config.pp2_misz_threshold
                    ),
                )
            }
        } else {
            (
                false,
                false,
                "Missense constraint Z-score not available".to_string(),
            )
        }
    } else {
        (
            false,
            false,
            "No gene constraint data available".to_string(),
        )
    };

    EvidenceCriterion {
        code: "PP2".to_string(),
        direction: EvidenceDirection::Pathogenic,
        strength: EvidenceStrength::Supporting,
        default_strength: EvidenceStrength::Supporting,
        met,
        evaluated,
        summary,
        details: serde_json::Value::Object(details),
    }
}

/// PP3: Multiple lines of computational evidence support a deleterious effect.
///
/// Uses ClinGen SVI calibrated REVEL thresholds. PP3 can be elevated to
/// Moderate or Strong strength based on REVEL score.
fn evaluate_pp3(
    input: &ClassificationInput,
    config: &AcmgConfig,
) -> EvidenceCriterion {
    let mut details = serde_json::Map::new();
    let mut evidence_lines: Vec<String> = Vec::new();

    // Primary: REVEL score (ClinGen SVI calibrated thresholds)
    let (revel_met, revel_strength) = if let Some(ref revel) = input.revel {
        if let Some(score) = revel.score {
            details.insert("revel_score".into(), serde_json::json!(score));
            if score >= config.pp3_revel_strong {
                evidence_lines.push(format!("REVEL={:.3} (Strong)", score));
                (true, Some(EvidenceStrength::Strong))
            } else if score >= config.pp3_revel_moderate {
                evidence_lines.push(format!("REVEL={:.3} (Moderate)", score));
                (true, Some(EvidenceStrength::Moderate))
            } else if score >= config.pp3_revel_supporting {
                evidence_lines.push(format!("REVEL={:.3} (Supporting)", score));
                (true, Some(EvidenceStrength::Supporting))
            } else {
                evidence_lines.push(format!("REVEL={:.3} (below threshold)", score));
                (false, None)
            }
        } else {
            (false, None)
        }
    } else {
        (false, None)
    };

    // Secondary: SIFT/PolyPhen/PhyloP/GERP consensus (when REVEL unavailable)
    let mut computational_votes = 0u8;
    let mut computational_total = 0u8;

    if let Some(ref dbnsfp) = input.dbnsfp {
        if let Some(sift) = dbnsfp.parse_sift() {
            details.insert("sift".into(), serde_json::json!(sift.prediction));
            computational_total += 1;
            if sift.prediction.contains("deleterious") && !sift.prediction.contains("tolerated") {
                computational_votes += 1;
                evidence_lines.push(format!("SIFT={}", sift.prediction));
            }
        }
        if let Some(pp) = dbnsfp.parse_polyphen() {
            details.insert("polyphen".into(), serde_json::json!(pp.prediction));
            computational_total += 1;
            if pp.prediction.contains("damaging") {
                computational_votes += 1;
                evidence_lines.push(format!("PolyPhen={}", pp.prediction));
            }
        }
    }

    if let Some(phylop) = input.phylop {
        details.insert("phylop".into(), serde_json::json!(phylop));
        computational_total += 1;
        if phylop > config.phylop_conserved {
            computational_votes += 1;
            evidence_lines.push(format!("PhyloP={:.2} (conserved)", phylop));
        }
    }

    if let Some(gerp) = input.gerp {
        details.insert("gerp".into(), serde_json::json!(gerp));
        computational_total += 1;
        if gerp > config.gerp_conserved {
            computational_votes += 1;
            evidence_lines.push(format!("GERP={:.2} (constrained)", gerp));
        }
    }

    // SpliceAI evidence (for variants near splice sites)
    if let Some(ref splice) = input.splice_ai {
        if let Some(max_ds) = splice.max_delta_score() {
            details.insert("spliceai_max_ds".into(), serde_json::json!(max_ds));
            if max_ds >= config.spliceai_pathogenic {
                evidence_lines.push(format!("SpliceAI max_ds={:.2}", max_ds));
            }
        }
    }

    // Determine final result
    let (met, strength) = if revel_met {
        (true, revel_strength.unwrap_or(EvidenceStrength::Supporting))
    } else if computational_total >= 3 && computational_votes >= 3 {
        // Consensus of at least 3 out of 4 computational tools
        (true, EvidenceStrength::Supporting)
    } else if let Some(ref splice) = input.splice_ai {
        if splice
            .max_delta_score()
            .map_or(false, |ds| ds >= config.spliceai_strong)
        {
            (true, EvidenceStrength::Strong)
        } else if splice
            .max_delta_score()
            .map_or(false, |ds| ds >= config.spliceai_pathogenic)
        {
            (true, EvidenceStrength::Supporting)
        } else {
            (false, EvidenceStrength::Supporting)
        }
    } else {
        (false, EvidenceStrength::Supporting)
    };

    let evaluated = input.revel.is_some()
        || computational_total >= 2
        || input.splice_ai.is_some();

    details.insert("evidence_lines".into(), serde_json::json!(evidence_lines));

    let summary = if met {
        format!(
            "Computational evidence supports deleterious effect ({}): {}",
            strength.as_str(),
            evidence_lines.join("; ")
        )
    } else if evaluated {
        "Computational evidence does not support deleterious effect".to_string()
    } else {
        "Insufficient computational prediction data available".to_string()
    };

    let code = if met && strength != EvidenceStrength::Supporting {
        format!("PP3_{}", strength.as_str())
    } else {
        "PP3".to_string()
    };

    EvidenceCriterion {
        code,
        direction: EvidenceDirection::Pathogenic,
        strength,
        default_strength: EvidenceStrength::Supporting,
        met,
        evaluated,
        summary,
        details: serde_json::Value::Object(details),
    }
}

/// PP4: Patient's phenotype or family history is highly specific for a disease
/// with a single genetic etiology.
fn evaluate_pp4(
    _input: &ClassificationInput,
    _config: &AcmgConfig,
) -> EvidenceCriterion {
    EvidenceCriterion {
        code: "PP4".to_string(),
        direction: EvidenceDirection::Pathogenic,
        strength: EvidenceStrength::Supporting,
        default_strength: EvidenceStrength::Supporting,
        met: false,
        evaluated: false,
        summary: "Requires patient HPO phenotype terms matched to disease-gene associations".to_string(),
        details: serde_json::Value::Null,
    }
}

/// PP5: Reputable source recently reports variant as pathogenic, but evidence not available.
///
/// Note: ClinGen SVI recommends against using PP5 without reviewing underlying evidence.
/// Disabled by default (config.use_pp5_bp6 = false).
fn evaluate_pp5(
    input: &ClassificationInput,
    _config: &AcmgConfig,
) -> EvidenceCriterion {
    let mut details = serde_json::Map::new();
    details.insert(
        "svi_note".into(),
        serde_json::json!("ClinGen SVI recommends against using PP5 without reviewing underlying evidence"),
    );

    let (met, evaluated, summary) = if let Some(ref clinvar) = input.clinvar {
        let stars = clinvar.review_stars();
        let is_pathogenic = clinvar.has_pathogenic();
        details.insert("clinvar_pathogenic".into(), serde_json::json!(is_pathogenic));
        details.insert("review_stars".into(), serde_json::json!(stars));

        if is_pathogenic && stars >= 2 {
            (
                true,
                true,
                format!(
                    "ClinVar pathogenic with {}-star review (use with caution per SVI)",
                    stars
                ),
            )
        } else {
            (
                false,
                true,
                format!("ClinVar not pathogenic or insufficient review ({} stars)", stars),
            )
        }
    } else {
        (false, false, "No ClinVar data available".to_string())
    };

    EvidenceCriterion {
        code: "PP5".to_string(),
        direction: EvidenceDirection::Pathogenic,
        strength: EvidenceStrength::Supporting,
        default_strength: EvidenceStrength::Supporting,
        met,
        evaluated,
        summary,
        details: serde_json::Value::Object(details),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sa_extract::{GnomadGeneData, RevelData, SpliceAiData};
    use fastvep_core::Impact;

    fn make_input_with_revel(score: f64) -> ClassificationInput {
        ClassificationInput {
            consequences: vec![Consequence::MissenseVariant],
            impact: Impact::Moderate,
            gene_symbol: Some("TEST".to_string()),
            is_canonical: true,
            amino_acids: None,
            protein_position: None,
            gnomad: None,
            clinvar: None,
            revel: Some(RevelData {
                score: Some(score),
            }),
            splice_ai: None,
            dbnsfp: None,
            phylop: None,
            gerp: None,
            gene_constraints: None,
            omim: None,
            clinvar_protein: None,
            in_repeat_region: None,
            proband_genotype: None,
            mother_genotype: None,
            father_genotype: None,
            companion_variants: vec![],
        }
    }

    #[test]
    fn test_pp3_revel_strong() {
        let input = make_input_with_revel(0.95);
        let result = evaluate_pp3(&input, &AcmgConfig::default());
        assert!(result.met);
        assert_eq!(result.strength, EvidenceStrength::Strong);
        assert!(result.code.contains("Strong"));
    }

    #[test]
    fn test_pp3_revel_moderate() {
        let input = make_input_with_revel(0.80);
        let result = evaluate_pp3(&input, &AcmgConfig::default());
        assert!(result.met);
        assert_eq!(result.strength, EvidenceStrength::Moderate);
    }

    #[test]
    fn test_pp3_revel_supporting() {
        let input = make_input_with_revel(0.65);
        let result = evaluate_pp3(&input, &AcmgConfig::default());
        assert!(result.met);
        assert_eq!(result.strength, EvidenceStrength::Supporting);
    }

    #[test]
    fn test_pp3_revel_below_threshold() {
        let input = make_input_with_revel(0.50);
        let result = evaluate_pp3(&input, &AcmgConfig::default());
        assert!(!result.met);
    }

    #[test]
    fn test_pp3_spliceai_strong() {
        let input = ClassificationInput {
            consequences: vec![Consequence::SpliceRegionVariant],
            impact: Impact::Low,
            gene_symbol: Some("TEST".to_string()),
            is_canonical: true,
            amino_acids: None,
            protein_position: None,
            gnomad: None,
            clinvar: None,
            revel: None,
            splice_ai: Some(SpliceAiData {
                ds_ag: Some(0.01),
                ds_al: Some(0.95),
                ds_dg: Some(0.02),
                ds_dl: Some(0.01),
                ..Default::default()
            }),
            dbnsfp: None,
            phylop: None,
            gerp: None,
            gene_constraints: None,
            omim: None,
            clinvar_protein: None,
            in_repeat_region: None,
            proband_genotype: None,
            mother_genotype: None,
            father_genotype: None,
            companion_variants: vec![],
        };
        let result = evaluate_pp3(&input, &AcmgConfig::default());
        assert!(result.met);
        assert_eq!(result.strength, EvidenceStrength::Strong);
    }

    #[test]
    fn test_pp2_high_misz() {
        let input = ClassificationInput {
            consequences: vec![Consequence::MissenseVariant],
            impact: Impact::Moderate,
            gene_symbol: Some("TEST".to_string()),
            is_canonical: true,
            amino_acids: None,
            protein_position: None,
            gnomad: None,
            clinvar: None,
            revel: None,
            splice_ai: None,
            dbnsfp: None,
            phylop: None,
            gerp: None,
            gene_constraints: Some(GnomadGeneData {
                mis_z: Some(4.5),
                ..Default::default()
            }),
            omim: None,
            clinvar_protein: None,
            in_repeat_region: None,
            proband_genotype: None,
            mother_genotype: None,
            father_genotype: None,
            companion_variants: vec![],
        };
        let result = evaluate_pp2(&input, &AcmgConfig::default());
        assert!(result.met);
    }

    #[test]
    fn test_pp2_low_misz() {
        let input = ClassificationInput {
            consequences: vec![Consequence::MissenseVariant],
            impact: Impact::Moderate,
            gene_symbol: Some("TEST".to_string()),
            is_canonical: true,
            amino_acids: None,
            protein_position: None,
            gnomad: None,
            clinvar: None,
            revel: None,
            splice_ai: None,
            dbnsfp: None,
            phylop: None,
            gerp: None,
            gene_constraints: Some(GnomadGeneData {
                mis_z: Some(1.5),
                ..Default::default()
            }),
            omim: None,
            clinvar_protein: None,
            in_repeat_region: None,
            proband_genotype: None,
            mother_genotype: None,
            father_genotype: None,
            companion_variants: vec![],
        };
        let result = evaluate_pp2(&input, &AcmgConfig::default());
        assert!(!result.met);
    }
}
