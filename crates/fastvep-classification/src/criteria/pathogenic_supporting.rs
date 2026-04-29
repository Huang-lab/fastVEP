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
/// Per ClinGen SVI calibration (Pejaver et al. 2022, AJHG): REVEL is applied
/// only to **missense** variants and uses calibrated bands for Supporting /
/// Moderate / Strong. Pejaver explicitly recommends a single calibrated tool
/// rather than ad-hoc consensus across SIFT/PolyPhen/PhyloP/GERP, so the
/// previous ≥3-of-4 consensus path has been removed; those scores are still
/// captured in `details` for transparency.
///
/// Per Walker 2023 (ClinGen SVI Splicing Subgroup): SpliceAI ≥ 0.2 yields
/// PP3 at *Supporting* strength only. SpliceAI alone does not reach Strong;
/// experimental RNA evidence (PVS1_RNA / PS3) is required for Strong splicing
/// claims.
fn evaluate_pp3(
    input: &ClassificationInput,
    config: &AcmgConfig,
) -> EvidenceCriterion {
    let mut details = serde_json::Map::new();
    let mut evidence_lines: Vec<String> = Vec::new();

    let is_missense = input
        .consequences
        .iter()
        .any(|c| matches!(c, Consequence::MissenseVariant));
    details.insert("is_missense".into(), serde_json::json!(is_missense));

    // Primary missense path: REVEL with calibrated bands (Pejaver 2022).
    // Only applied to missense variants — REVEL is undefined / uncalibrated for
    // other consequence types.
    let (revel_met, revel_strength) = if is_missense {
        if let Some(ref revel) = input.revel {
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
        }
    } else {
        if let Some(ref revel) = input.revel {
            if let Some(score) = revel.score {
                details.insert("revel_score".into(), serde_json::json!(score));
                details.insert(
                    "revel_skipped_reason".into(),
                    serde_json::json!(
                        "REVEL is calibrated for missense only (Pejaver 2022); not applied to non-missense consequences"
                    ),
                );
            }
        }
        (false, None)
    };

    // Capture transparency-only secondary scores (no longer used for PP3 firing
    // per Pejaver 2022 — single calibrated tool only).
    if let Some(ref dbnsfp) = input.dbnsfp {
        if let Some(sift) = dbnsfp.parse_sift() {
            details.insert("sift".into(), serde_json::json!(sift.prediction));
        }
        if let Some(pp) = dbnsfp.parse_polyphen() {
            details.insert("polyphen".into(), serde_json::json!(pp.prediction));
        }
    }
    if let Some(phylop) = input.phylop {
        details.insert("phylop".into(), serde_json::json!(phylop));
    }
    if let Some(gerp) = input.gerp {
        details.insert("gerp".into(), serde_json::json!(gerp));
    }

    // Splicing path: SpliceAI ≥ spliceai_pathogenic → PP3 Supporting.
    // Per Walker 2023, SpliceAI alone does not reach Strong even at very high
    // delta scores; that requires experimental RNA evidence (PVS1_RNA / PS3).
    let splice_supporting = if let Some(ref splice) = input.splice_ai {
        if let Some(max_ds) = splice.max_delta_score() {
            details.insert("spliceai_max_ds".into(), serde_json::json!(max_ds));
            if max_ds >= config.spliceai_pathogenic {
                evidence_lines.push(format!(
                    "SpliceAI max_ds={:.2} (Supporting per Walker 2023)",
                    max_ds
                ));
                true
            } else {
                false
            }
        } else {
            false
        }
    } else {
        false
    };

    // Resolve final strength.
    // - REVEL (missense) takes precedence and is already strength-graded.
    // - Otherwise, SpliceAI Supporting may fire for any consequence (canonical
    //   splice variants typically also fire PVS1; anti-double-counting between
    //   PP3 and PVS1 for splice is handled in `criteria::reconcile_evidence`).
    let (met, strength, source) = if revel_met {
        (
            true,
            revel_strength.unwrap_or(EvidenceStrength::Supporting),
            "revel_missense",
        )
    } else if splice_supporting {
        (true, EvidenceStrength::Supporting, "spliceai")
    } else {
        (false, EvidenceStrength::Supporting, "none")
    };

    let evaluated = (is_missense && input.revel.is_some()) || input.splice_ai.is_some();

    details.insert("pp3_source".into(), serde_json::json!(source));
    details.insert("evidence_lines".into(), serde_json::json!(evidence_lines));

    let summary = if met {
        format!(
            "Computational evidence supports deleterious effect ({}): {}",
            strength.as_str(),
            evidence_lines.join("; ")
        )
    } else if evaluated {
        "Computational evidence does not support deleterious effect".to_string()
    } else if !is_missense && input.splice_ai.is_none() {
        "PP3 requires REVEL (missense) or SpliceAI; neither available".to_string()
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
            hgvs_c: None,
            predicted_nmd: None,
            protein_truncation_pct: None,
            is_last_exon: None,
            in_critical_region: None,
            alt_start_codon_distance: None,
            same_splice_position_pathogenic: None,
            in_repeat_region: None,
            at_exon_edge: None,
            intronic_offset: None,
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
    fn test_pp3_spliceai_caps_at_supporting() {
        // Walker 2023 (ClinGen SVI Splicing Subgroup): SpliceAI alone tops
        // out at PP3 Supporting, even at very high delta scores. Strong
        // splicing evidence requires experimental RNA assay (PVS1_RNA / PS3).
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
            hgvs_c: None,
            predicted_nmd: None,
            protein_truncation_pct: None,
            is_last_exon: None,
            in_critical_region: None,
            alt_start_codon_distance: None,
            same_splice_position_pathogenic: None,
            in_repeat_region: None,
            at_exon_edge: None,
            intronic_offset: None,
            proband_genotype: None,
            mother_genotype: None,
            father_genotype: None,
            companion_variants: vec![],
        };
        let result = evaluate_pp3(&input, &AcmgConfig::default());
        assert!(result.met);
        assert_eq!(result.strength, EvidenceStrength::Supporting);
        assert_eq!(result.code, "PP3");
    }

    #[test]
    fn test_pp3_revel_ignored_for_non_missense() {
        // Pejaver 2022: REVEL is calibrated for missense only. A high REVEL
        // score on a frameshift / synonymous / splice variant must NOT fire PP3.
        let mut input = make_input_with_revel(0.95);
        input.consequences = vec![Consequence::FrameshiftVariant];
        let result = evaluate_pp3(&input, &AcmgConfig::default());
        assert!(!result.met);
        // REVEL score is still recorded in details for transparency.
        let details = result.details.as_object().unwrap();
        assert!(details.contains_key("revel_score"));
        assert!(details.contains_key("revel_skipped_reason"));
    }

    #[test]
    fn test_pp3_consensus_dropped() {
        // Pre-PR1 behavior allowed ≥3-of-4 SIFT/PolyPhen/PhyloP/GERP consensus
        // to fire PP3 Supporting in the absence of REVEL. Pejaver 2022 explicitly
        // recommends a single calibrated tool, not ad-hoc consensus, so this
        // path has been removed. Without REVEL or SpliceAI, PP3 must NOT fire
        // even when SIFT/PolyPhen/PhyloP/GERP all signal pathogenic.
        let mut input = ClassificationInput {
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
            phylop: Some(5.0),
            gerp: Some(5.0),
            gene_constraints: None,
            omim: None,
            clinvar_protein: None,
            in_repeat_region: None,
            proband_genotype: None,
            mother_genotype: None,
            father_genotype: None,
            companion_variants: vec![],
            at_exon_edge: None,
            intronic_offset: None,
            hgvs_c: None,
            predicted_nmd: None,
            protein_truncation_pct: None,
            is_last_exon: None,
            in_critical_region: None,
            alt_start_codon_distance: None,
            same_splice_position_pathogenic: None,
        };
        // Synthesize a dbNSFP entry with deleterious SIFT + damaging PolyPhen
        // by going through the same JSON path the evaluator uses.
        input.dbnsfp = Some(crate::sa_extract::DbNsfpData {
            sift: Some("deleterious(0.000)".to_string()),
            polyphen: Some("probably_damaging(0.998)".to_string()),
        });
        let result = evaluate_pp3(&input, &AcmgConfig::default());
        assert!(!result.met);
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
            hgvs_c: None,
            predicted_nmd: None,
            protein_truncation_pct: None,
            is_last_exon: None,
            in_critical_region: None,
            alt_start_codon_distance: None,
            same_splice_position_pathogenic: None,
            in_repeat_region: None,
            at_exon_edge: None,
            intronic_offset: None,
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
            hgvs_c: None,
            predicted_nmd: None,
            protein_truncation_pct: None,
            is_last_exon: None,
            in_critical_region: None,
            alt_start_codon_distance: None,
            same_splice_position_pathogenic: None,
            in_repeat_region: None,
            at_exon_edge: None,
            intronic_offset: None,
            proband_genotype: None,
            mother_genotype: None,
            father_genotype: None,
            companion_variants: vec![],
        };
        let result = evaluate_pp2(&input, &AcmgConfig::default());
        assert!(!result.met);
    }
}
