use fastvep_core::Consequence;

use crate::config::AcmgConfig;
use crate::sa_extract::ClassificationInput;
use crate::types::{EvidenceCriterion, EvidenceDirection, EvidenceStrength};

/// Evaluate all benign supporting criteria: BP1, BP2, BP3, BP4, BP5, BP6, BP7.
pub fn evaluate_all(
    input: &ClassificationInput,
    config: &AcmgConfig,
) -> Vec<EvidenceCriterion> {
    let mut criteria = vec![
        evaluate_bp1(input, config),
        evaluate_bp2(input, config),
        evaluate_bp3(input, config),
        evaluate_bp4(input, config),
        evaluate_bp5(input, config),
        evaluate_bp7(input, config),
    ];
    if config.use_pp5_bp6 {
        criteria.push(evaluate_bp6(input, config));
    }
    criteria
}

/// BP1: Missense variant in a gene for which primarily truncating variants are known
/// to cause disease.
///
/// Approximated: gene has high pLI (LOF-intolerant) but low missense Z (missense-tolerant).
fn evaluate_bp1(
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
        let pli_high = gc
            .pli
            .map_or(false, |p| p >= config.pli_lof_intolerant);
        let misz_low = gc.mis_z.map_or(false, |z| z < 2.0);

        if let Some(pli) = gc.pli {
            details.insert("pLI".into(), serde_json::json!(pli));
        }
        if let Some(misz) = gc.mis_z {
            details.insert("misZ".into(), serde_json::json!(misz));
        }

        if pli_high && misz_low {
            (
                true,
                true,
                format!(
                    "Missense in gene where primarily truncating variants cause disease (pLI={:.2}, misZ={:.2})",
                    gc.pli.unwrap_or(0.0),
                    gc.mis_z.unwrap_or(0.0)
                ),
            )
        } else {
            (
                false,
                true,
                format!(
                    "Gene does not meet BP1 criteria (pLI={:.2}, misZ={:.2}; needs pLI>={:.1} and misZ<2.0)",
                    gc.pli.unwrap_or(0.0),
                    gc.mis_z.unwrap_or(0.0),
                    config.pli_lof_intolerant
                ),
            )
        }
    } else {
        (false, false, "No gene constraint data available".to_string())
    };

    EvidenceCriterion {
        code: "BP1".to_string(),
        direction: EvidenceDirection::Benign,
        strength: EvidenceStrength::Supporting,
        default_strength: EvidenceStrength::Supporting,
        met,
        evaluated,
        summary,
        details: serde_json::Value::Object(details),
    }
}

/// BP2: Observed in trans with a pathogenic variant for a fully penetrant dominant
/// gene/disorder or observed in cis with a pathogenic variant in any inheritance pattern.
///
/// Two triggers:
/// 1. Dominant disorders: variant is in trans with a ClinVar pathogenic variant (phased).
/// 2. Any inheritance: variant is in cis with a ClinVar pathogenic variant (phased).
fn evaluate_bp2(
    input: &ClassificationInput,
    _config: &AcmgConfig,
) -> EvidenceCriterion {
    let mut details = serde_json::Map::new();

    if input.companion_variants.is_empty() {
        return EvidenceCriterion {
            code: "BP2".to_string(),
            direction: EvidenceDirection::Benign,
            strength: EvidenceStrength::Supporting,
            default_strength: EvidenceStrength::Supporting,
            met: false,
            evaluated: false,
            summary: "No companion variants available for in-trans/in-cis analysis".to_string(),
            details: serde_json::Value::Null,
        };
    }

    let is_dominant = input
        .omim
        .as_ref()
        .map_or(false, |o| o.has_dominant_inheritance());
    details.insert("is_dominant_gene".into(), serde_json::json!(is_dominant));

    // Check for in-cis with pathogenic (any inheritance pattern)
    let in_cis_pathogenic: Vec<&crate::sa_extract::CompanionVariant> = input
        .companion_variants
        .iter()
        .filter(|cv| cv.is_clinvar_pathogenic && cv.is_in_trans == Some(false))
        .collect();

    // Check for in-trans with pathogenic in dominant gene
    let in_trans_pathogenic: Vec<&crate::sa_extract::CompanionVariant> = input
        .companion_variants
        .iter()
        .filter(|cv| cv.is_clinvar_pathogenic && cv.is_in_trans == Some(true))
        .collect();

    details.insert("in_cis_pathogenic_count".into(), serde_json::json!(in_cis_pathogenic.len()));
    details.insert("in_trans_pathogenic_count".into(), serde_json::json!(in_trans_pathogenic.len()));

    // Collect HGVSc for reporting
    let cis_ids: Vec<String> = in_cis_pathogenic
        .iter()
        .filter_map(|cv| cv.hgvsc.clone())
        .collect();
    let trans_ids: Vec<String> = in_trans_pathogenic
        .iter()
        .filter_map(|cv| cv.hgvsc.clone())
        .collect();
    if !cis_ids.is_empty() {
        details.insert("in_cis_hgvsc".into(), serde_json::json!(cis_ids));
    }
    if !trans_ids.is_empty() {
        details.insert("in_trans_hgvsc".into(), serde_json::json!(trans_ids));
    }

    // Trigger 1: In cis with pathogenic (any inheritance)
    if !in_cis_pathogenic.is_empty() {
        return EvidenceCriterion {
            code: "BP2".to_string(),
            direction: EvidenceDirection::Benign,
            strength: EvidenceStrength::Supporting,
            default_strength: EvidenceStrength::Supporting,
            met: true,
            evaluated: true,
            summary: format!(
                "Observed in cis (same haplotype) with ClinVar pathogenic variant ({} companion(s))",
                in_cis_pathogenic.len()
            ),
            details: serde_json::Value::Object(details),
        };
    }

    // Trigger 2: In trans with pathogenic in dominant gene
    if is_dominant && !in_trans_pathogenic.is_empty() {
        return EvidenceCriterion {
            code: "BP2".to_string(),
            direction: EvidenceDirection::Benign,
            strength: EvidenceStrength::Supporting,
            default_strength: EvidenceStrength::Supporting,
            met: true,
            evaluated: true,
            summary: format!(
                "Observed in trans with ClinVar pathogenic variant in dominant gene ({} companion(s))",
                in_trans_pathogenic.len()
            ),
            details: serde_json::Value::Object(details),
        };
    }

    // Check if we have any phased data at all (to determine if evaluated)
    let has_phased_data = input
        .companion_variants
        .iter()
        .any(|cv| cv.is_in_trans.is_some());

    let summary = if has_phased_data {
        "Phased companion variants do not meet BP2 criteria".to_string()
    } else {
        "Companion variants are unphased; BP2 requires phased data to determine cis/trans".to_string()
    };

    EvidenceCriterion {
        code: "BP2".to_string(),
        direction: EvidenceDirection::Benign,
        strength: EvidenceStrength::Supporting,
        default_strength: EvidenceStrength::Supporting,
        met: false,
        evaluated: has_phased_data,
        summary,
        details: serde_json::Value::Object(details),
    }
}

/// BP3: In-frame deletions/insertions in a repetitive region without a known function.
///
/// Uses RepeatMasker interval annotations when available.
fn evaluate_bp3(
    input: &ClassificationInput,
    _config: &AcmgConfig,
) -> EvidenceCriterion {
    let is_inframe_indel = input.consequences.iter().any(|c| {
        matches!(
            c,
            Consequence::InframeInsertion | Consequence::InframeDeletion
        )
    });

    let mut details = serde_json::Map::new();
    details.insert(
        "is_inframe_indel".into(),
        serde_json::json!(is_inframe_indel),
    );

    if let Some(in_repeat) = input.in_repeat_region {
        details.insert("in_repeat_region".into(), serde_json::json!(in_repeat));
        let met = is_inframe_indel && in_repeat;
        let summary = if met {
            "In-frame indel in a repetitive region".to_string()
        } else if is_inframe_indel && !in_repeat {
            "In-frame indel but not in a repetitive region".to_string()
        } else {
            "Not an in-frame insertion or deletion".to_string()
        };

        EvidenceCriterion {
            code: "BP3".to_string(),
            direction: EvidenceDirection::Benign,
            strength: EvidenceStrength::Supporting,
            default_strength: EvidenceStrength::Supporting,
            met,
            evaluated: true,
            summary,
            details: serde_json::Value::Object(details),
        }
    } else {
        let summary = if is_inframe_indel {
            "In-frame indel detected, but repeat region data not available (load RepeatMasker .osi)".to_string()
        } else {
            "Not an in-frame insertion or deletion".to_string()
        };

        EvidenceCriterion {
            code: "BP3".to_string(),
            direction: EvidenceDirection::Benign,
            strength: EvidenceStrength::Supporting,
            default_strength: EvidenceStrength::Supporting,
            met: false,
            evaluated: !is_inframe_indel, // evaluated if not applicable; not evaluated if it's relevant but no data
            summary,
            details: serde_json::Value::Object(details),
        }
    }
}

/// BP4: Multiple lines of computational evidence suggest no impact on gene or gene product.
///
/// Per ClinGen SVI calibration (Pejaver et al. 2022, AJHG): REVEL is applied
/// only to **missense** variants and uses calibrated bands for Supporting /
/// Moderate / Strong / Very Strong. The Very Strong band (REVEL ≤ 0.003) is
/// reached only by REVEL among the 13 calibrated tools. Pejaver explicitly
/// recommends a single calibrated tool over ad-hoc consensus, so the previous
/// SIFT/PolyPhen/PhyloP ≥2-of-3 consensus path has been removed; those scores
/// remain in `details` for transparency.
///
/// Per Walker 2023 (ClinGen SVI Splicing Subgroup): SpliceAI ≤ 0.1 yields BP4
/// at Supporting strength when a SpliceAI score is available (scores between
/// 0.1 and 0.2 are uninformative).
fn evaluate_bp4(
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

    // Primary missense path: REVEL with calibrated bands (Pejaver 2022),
    // including the Very Strong band (REVEL ≤ 0.003) which only REVEL reaches.
    let (revel_met, revel_strength) = if is_missense {
        if let Some(ref revel) = input.revel {
            if let Some(score) = revel.score {
                details.insert("revel_score".into(), serde_json::json!(score));
                if score <= config.bp4_revel_very_strong {
                    evidence_lines.push(format!("REVEL={:.3} (Very Strong benign)", score));
                    (true, Some(EvidenceStrength::VeryStrong))
                } else if score <= config.bp4_revel_strong {
                    evidence_lines.push(format!("REVEL={:.3} (Strong benign)", score));
                    (true, Some(EvidenceStrength::Strong))
                } else if score <= config.bp4_revel_moderate {
                    evidence_lines.push(format!("REVEL={:.3} (Moderate benign)", score));
                    (true, Some(EvidenceStrength::Moderate))
                } else if score <= config.bp4_revel_supporting {
                    evidence_lines.push(format!("REVEL={:.3} (Supporting benign)", score));
                    (true, Some(EvidenceStrength::Supporting))
                } else {
                    evidence_lines.push(format!("REVEL={:.3} (above benign threshold)", score));
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

    // Capture transparency-only secondary scores (no longer used for BP4 firing
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

    // Splicing path: SpliceAI ≤ spliceai_benign → BP4 Supporting (Walker 2023).
    let splice_supporting = if let Some(ref splice) = input.splice_ai {
        if let Some(max_ds) = splice.max_delta_score() {
            details.insert("spliceai_max_ds".into(), serde_json::json!(max_ds));
            if max_ds <= config.spliceai_benign {
                evidence_lines.push(format!(
                    "SpliceAI max_ds={:.2} (Supporting benign per Walker 2023)",
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

    // Resolve final strength. REVEL (missense, already strength-graded) takes
    // precedence over splice. If both fire, prefer the stronger.
    let (met, strength) = match (revel_met, splice_supporting) {
        (true, _) => (
            true,
            revel_strength.unwrap_or(EvidenceStrength::Supporting),
        ),
        (false, true) => (true, EvidenceStrength::Supporting),
        (false, false) => (false, EvidenceStrength::Supporting),
    };

    let evaluated = (is_missense && input.revel.is_some()) || input.splice_ai.is_some();
    details.insert("evidence_lines".into(), serde_json::json!(evidence_lines));

    let summary = if met {
        format!(
            "Computational evidence supports benign ({}): {}",
            strength.as_str(),
            evidence_lines.join("; ")
        )
    } else if evaluated {
        "Computational evidence does not support benign classification".to_string()
    } else if !is_missense && input.splice_ai.is_none() {
        "BP4 requires REVEL (missense) or SpliceAI; neither available".to_string()
    } else {
        "Insufficient computational prediction data".to_string()
    };

    let code = if met && strength != EvidenceStrength::Supporting {
        format!("BP4_{}", strength.as_str())
    } else {
        "BP4".to_string()
    };

    EvidenceCriterion {
        code,
        direction: EvidenceDirection::Benign,
        strength,
        default_strength: EvidenceStrength::Supporting,
        met,
        evaluated,
        summary,
        details: serde_json::Value::Object(details),
    }
}

/// BP5: Variant found in a case with an alternate molecular basis for disease.
fn evaluate_bp5(
    _input: &ClassificationInput,
    _config: &AcmgConfig,
) -> EvidenceCriterion {
    EvidenceCriterion {
        code: "BP5".to_string(),
        direction: EvidenceDirection::Benign,
        strength: EvidenceStrength::Supporting,
        default_strength: EvidenceStrength::Supporting,
        met: false,
        evaluated: false,
        summary: "Requires case-level knowledge of other confirmed pathogenic variants explaining the phenotype".to_string(),
        details: serde_json::Value::Null,
    }
}

/// BP6: Reputable source recently reports variant as benign, but evidence not available.
///
/// Note: ClinGen SVI recommends against using BP6 without reviewing underlying evidence.
/// Disabled by default (config.use_pp5_bp6 = false).
fn evaluate_bp6(
    input: &ClassificationInput,
    _config: &AcmgConfig,
) -> EvidenceCriterion {
    let mut details = serde_json::Map::new();
    details.insert(
        "svi_note".into(),
        serde_json::json!("ClinGen SVI recommends against using BP6 without reviewing underlying evidence"),
    );

    let (met, evaluated, summary) = if let Some(ref clinvar) = input.clinvar {
        let stars = clinvar.review_stars();
        let is_benign = clinvar.has_benign();
        details.insert("clinvar_benign".into(), serde_json::json!(is_benign));
        details.insert("review_stars".into(), serde_json::json!(stars));

        if is_benign && stars >= 2 {
            (
                true,
                true,
                format!(
                    "ClinVar benign with {}-star review (use with caution per SVI)",
                    stars
                ),
            )
        } else {
            (
                false,
                true,
                format!(
                    "ClinVar not benign or insufficient review ({} stars)",
                    stars
                ),
            )
        }
    } else {
        (false, false, "No ClinVar data available".to_string())
    };

    EvidenceCriterion {
        code: "BP6".to_string(),
        direction: EvidenceDirection::Benign,
        strength: EvidenceStrength::Supporting,
        default_strength: EvidenceStrength::Supporting,
        met,
        evaluated,
        summary,
        details: serde_json::Value::Object(details),
    }
}

/// BP7: A synonymous (silent) variant for which splicing prediction algorithms
/// predict no impact to the splice consensus sequence AND the nucleotide is not
/// highly conserved.
fn evaluate_bp7(
    input: &ClassificationInput,
    config: &AcmgConfig,
) -> EvidenceCriterion {
    let is_synonymous = input
        .consequences
        .iter()
        .any(|c| matches!(c, Consequence::SynonymousVariant));

    let mut details = serde_json::Map::new();
    details.insert("is_synonymous".into(), serde_json::json!(is_synonymous));

    if !is_synonymous {
        return EvidenceCriterion {
            code: "BP7".to_string(),
            direction: EvidenceDirection::Benign,
            strength: EvidenceStrength::Supporting,
            default_strength: EvidenceStrength::Supporting,
            met: false,
            evaluated: true,
            summary: "Not a synonymous variant".to_string(),
            details: serde_json::Value::Object(details),
        };
    }

    // Check SpliceAI: no predicted splice impact
    let no_splice_impact = if let Some(ref splice) = input.splice_ai {
        let max_ds = splice.max_delta_score().unwrap_or(0.0);
        details.insert("spliceai_max_ds".into(), serde_json::json!(max_ds));
        max_ds < config.spliceai_pathogenic
    } else {
        // If no SpliceAI data, we conservatively assume no impact
        details.insert("spliceai_max_ds".into(), serde_json::Value::Null);
        true
    };

    // Check PhyloP: not highly conserved
    let not_conserved = if let Some(phylop) = input.phylop {
        details.insert("phylop".into(), serde_json::json!(phylop));
        phylop < config.phylop_conserved
    } else {
        // Conservative: if no data, don't assume not conserved
        details.insert("phylop".into(), serde_json::Value::Null);
        false
    };

    let splice_ai_available = input.splice_ai.is_some();
    let phylop_available = input.phylop.is_some();
    let evaluated = splice_ai_available || phylop_available;
    let met = is_synonymous && no_splice_impact && not_conserved;

    let summary = if met {
        format!(
            "Synonymous variant with no predicted splice impact (SpliceAI max_ds={:.2}) and not conserved (PhyloP={:.2})",
            input.splice_ai.as_ref().and_then(|s| s.max_delta_score()).unwrap_or(0.0),
            input.phylop.unwrap_or(0.0)
        )
    } else if is_synonymous && !no_splice_impact {
        "Synonymous but predicted to affect splicing".to_string()
    } else if is_synonymous && !not_conserved {
        "Synonymous but position is highly conserved".to_string()
    } else {
        "Synonymous variant but insufficient data to confirm BP7".to_string()
    };

    EvidenceCriterion {
        code: "BP7".to_string(),
        direction: EvidenceDirection::Benign,
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

    fn make_input(
        consequences: Vec<Consequence>,
        revel_score: Option<f64>,
        splice_ai_max_ds: Option<f64>,
        phylop: Option<f64>,
        gene_constraints: Option<GnomadGeneData>,
    ) -> ClassificationInput {
        ClassificationInput {
            consequences,
            impact: Impact::Moderate,
            gene_symbol: Some("TEST".to_string()),
            is_canonical: true,
            amino_acids: None,
            protein_position: None,
            gnomad: None,
            clinvar: None,
            revel: revel_score.map(|s| RevelData { score: Some(s) }),
            splice_ai: splice_ai_max_ds.map(|ds| SpliceAiData {
                ds_al: Some(ds),
                ..Default::default()
            }),
            dbnsfp: None,
            phylop,
            gerp: None,
            gene_constraints,
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
    fn test_bp4_revel_strong_benign() {
        let input = make_input(
            vec![Consequence::MissenseVariant],
            Some(0.01),
            None,
            None,
            None,
        );
        let result = evaluate_bp4(&input, &AcmgConfig::default());
        assert!(result.met);
        assert_eq!(result.strength, EvidenceStrength::Strong);
    }

    #[test]
    fn test_bp4_revel_very_strong_benign() {
        // Pejaver 2022: REVEL ≤ 0.003 → Very Strong benign. Only REVEL reaches
        // this band among the 13 calibrated tools.
        let input = make_input(
            vec![Consequence::MissenseVariant],
            Some(0.002),
            None,
            None,
            None,
        );
        let result = evaluate_bp4(&input, &AcmgConfig::default());
        assert!(result.met);
        assert_eq!(result.strength, EvidenceStrength::VeryStrong);
        assert_eq!(result.code, "BP4_Very_Strong");
    }

    #[test]
    fn test_bp4_revel_ignored_for_non_missense() {
        // REVEL is calibrated for missense only. A low REVEL score on a
        // synonymous variant must NOT fire BP4 (BP7 is the synonymous code).
        let input = make_input(
            vec![Consequence::SynonymousVariant],
            Some(0.001),
            None,
            None,
            None,
        );
        let result = evaluate_bp4(&input, &AcmgConfig::default());
        assert!(!result.met);
    }

    #[test]
    fn test_bp4_spliceai_low_supporting() {
        // Walker 2023: SpliceAI ≤ 0.1 → BP4 Supporting for non-canonical
        // splice context. Non-missense variant with low SpliceAI alone fires.
        let input = make_input(
            vec![Consequence::IntronVariant],
            None,
            Some(0.05),
            None,
            None,
        );
        let result = evaluate_bp4(&input, &AcmgConfig::default());
        assert!(result.met);
        assert_eq!(result.strength, EvidenceStrength::Supporting);
    }

    #[test]
    fn test_bp4_spliceai_uninformative_zone() {
        // Walker 2023: SpliceAI scores 0.10 < ds < 0.20 are uninformative —
        // BP4 must not fire and PP3 (splice) must not fire.
        let input = make_input(
            vec![Consequence::IntronVariant],
            None,
            Some(0.15),
            None,
            None,
        );
        let result = evaluate_bp4(&input, &AcmgConfig::default());
        assert!(!result.met);
    }

    #[test]
    fn test_bp4_revel_supporting_benign() {
        let input = make_input(
            vec![Consequence::MissenseVariant],
            Some(0.25),
            None,
            None,
            None,
        );
        let result = evaluate_bp4(&input, &AcmgConfig::default());
        assert!(result.met);
        assert_eq!(result.strength, EvidenceStrength::Supporting);
    }

    #[test]
    fn test_bp4_revel_above_benign() {
        let input = make_input(
            vec![Consequence::MissenseVariant],
            Some(0.50),
            None,
            None,
            None,
        );
        let result = evaluate_bp4(&input, &AcmgConfig::default());
        assert!(!result.met);
    }

    #[test]
    fn test_bp7_synonymous_no_splice_not_conserved() {
        let input = make_input(
            vec![Consequence::SynonymousVariant],
            None,
            Some(0.05),
            Some(0.5),
            None,
        );
        let result = evaluate_bp7(&input, &AcmgConfig::default());
        assert!(result.met);
    }

    #[test]
    fn test_bp7_synonymous_splice_impact() {
        let input = make_input(
            vec![Consequence::SynonymousVariant],
            None,
            Some(0.50),
            Some(0.5),
            None,
        );
        let result = evaluate_bp7(&input, &AcmgConfig::default());
        assert!(!result.met);
    }

    #[test]
    fn test_bp7_synonymous_conserved() {
        let input = make_input(
            vec![Consequence::SynonymousVariant],
            None,
            Some(0.05),
            Some(5.0),
            None,
        );
        let result = evaluate_bp7(&input, &AcmgConfig::default());
        assert!(!result.met);
    }

    #[test]
    fn test_bp1_lof_intolerant_missense_tolerant() {
        let input = make_input(
            vec![Consequence::MissenseVariant],
            None,
            None,
            None,
            Some(GnomadGeneData {
                pli: Some(0.99),
                mis_z: Some(1.0),
                ..Default::default()
            }),
        );
        let result = evaluate_bp1(&input, &AcmgConfig::default());
        assert!(result.met);
    }

    #[test]
    fn test_bp1_not_lof_intolerant() {
        let input = make_input(
            vec![Consequence::MissenseVariant],
            None,
            None,
            None,
            Some(GnomadGeneData {
                pli: Some(0.5),
                mis_z: Some(1.0),
                ..Default::default()
            }),
        );
        let result = evaluate_bp1(&input, &AcmgConfig::default());
        assert!(!result.met);
    }
}
