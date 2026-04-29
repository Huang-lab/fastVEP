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
/// Uses the ClinVar protein-position index (.oga) to check if pathogenic variants
/// with the same amino acid change exist at the same protein position.
fn evaluate_ps1(
    input: &ClassificationInput,
    _config: &AcmgConfig,
) -> EvidenceCriterion {
    let is_missense = input
        .consequences
        .iter()
        .any(|c| matches!(c, fastvep_core::Consequence::MissenseVariant));

    let mut details = serde_json::Map::new();
    details.insert("is_missense".into(), serde_json::json!(is_missense));

    if !is_missense {
        return EvidenceCriterion {
            code: "PS1".to_string(),
            direction: EvidenceDirection::Pathogenic,
            strength: EvidenceStrength::Strong,
            default_strength: EvidenceStrength::Strong,
            met: false,
            evaluated: true,
            summary: "Not a missense variant".to_string(),
            details: serde_json::Value::Object(details),
        };
    }

    let (prot_pos, ref_aa, alt_aa) = match (&input.protein_position, &input.amino_acids) {
        (Some(pos), Some((r, a))) => (*pos, r.as_str(), a.as_str()),
        _ => {
            return EvidenceCriterion {
                code: "PS1".to_string(),
                direction: EvidenceDirection::Pathogenic,
                strength: EvidenceStrength::Strong,
                default_strength: EvidenceStrength::Strong,
                met: false,
                evaluated: false,
                summary: "Protein position or amino acid change not available".to_string(),
                details: serde_json::Value::Object(details),
            };
        }
    };

    details.insert("protein_position".into(), serde_json::json!(prot_pos));
    details.insert("ref_aa".into(), serde_json::json!(ref_aa));
    details.insert("alt_aa".into(), serde_json::json!(alt_aa));

    if let Some(ref cpd) = input.clinvar_protein {
        let matches: Vec<&crate::sa_extract::ClinvarProteinVariant> = cpd
            .protein_variants
            .iter()
            .filter(|v| v.pos == prot_pos && v.alt_aa == alt_aa && v.sig.to_lowercase().contains("pathogenic"))
            .collect();

        details.insert("matching_pathogenic_count".into(), serde_json::json!(matches.len()));

        if !matches.is_empty() {
            return EvidenceCriterion {
                code: "PS1".to_string(),
                direction: EvidenceDirection::Pathogenic,
                strength: EvidenceStrength::Strong,
                default_strength: EvidenceStrength::Strong,
                met: true,
                evaluated: true,
                summary: format!(
                    "Same amino acid change (p.{}{}{}>) is pathogenic in ClinVar ({} entries at protein position {})",
                    ref_aa, prot_pos, alt_aa, matches.len(), prot_pos
                ),
                details: serde_json::Value::Object(details),
            };
        }

        EvidenceCriterion {
            code: "PS1".to_string(),
            direction: EvidenceDirection::Pathogenic,
            strength: EvidenceStrength::Strong,
            default_strength: EvidenceStrength::Strong,
            met: false,
            evaluated: true,
            summary: format!(
                "No pathogenic ClinVar variant with same AA change at position {}",
                prot_pos
            ),
            details: serde_json::Value::Object(details),
        }
    } else {
        EvidenceCriterion {
            code: "PS1".to_string(),
            direction: EvidenceDirection::Pathogenic,
            strength: EvidenceStrength::Strong,
            default_strength: EvidenceStrength::Strong,
            met: false,
            evaluated: false,
            summary: "ClinVar protein-position index not available".to_string(),
            details: serde_json::Value::Object(details),
        }
    }
}

/// PS2: De novo (both maternity and paternity confirmed) in a patient with the disease.
///
/// Requires trio VCF with both parents. Proband carries the variant, both parents are
/// homozygous reference, and all three pass depth/quality thresholds.
fn evaluate_ps2(
    input: &ClassificationInput,
    config: &AcmgConfig,
) -> EvidenceCriterion {
    let mut details = serde_json::Map::new();

    let trio = match &config.trio {
        Some(t) => t,
        None => {
            return EvidenceCriterion {
                code: "PS2".to_string(),
                direction: EvidenceDirection::Pathogenic,
                strength: EvidenceStrength::Strong,
                default_strength: EvidenceStrength::Strong,
                met: false,
                evaluated: false,
                summary: "Requires trio VCF with --proband/--mother/--father sample names to assess de novo status".to_string(),
                details: serde_json::Value::Null,
            };
        }
    };

    // PS2 requires BOTH parents
    if trio.mother.is_none() || trio.father.is_none() {
        return EvidenceCriterion {
            code: "PS2".to_string(),
            direction: EvidenceDirection::Pathogenic,
            strength: EvidenceStrength::Strong,
            default_strength: EvidenceStrength::Strong,
            met: false,
            evaluated: false,
            summary: "PS2 requires both mother and father in trio configuration; use PM6 for partial trio".to_string(),
            details: serde_json::Value::Null,
        };
    }

    let proband_gt = match &input.proband_genotype {
        Some(gt) => gt,
        None => {
            return EvidenceCriterion {
                code: "PS2".to_string(),
                direction: EvidenceDirection::Pathogenic,
                strength: EvidenceStrength::Strong,
                default_strength: EvidenceStrength::Strong,
                met: false,
                evaluated: false,
                summary: "Proband genotype not available for this variant".to_string(),
                details: serde_json::Value::Null,
            };
        }
    };

    let mother_gt = match &input.mother_genotype {
        Some(gt) => gt,
        None => {
            return EvidenceCriterion {
                code: "PS2".to_string(),
                direction: EvidenceDirection::Pathogenic,
                strength: EvidenceStrength::Strong,
                default_strength: EvidenceStrength::Strong,
                met: false,
                evaluated: false,
                summary: "Mother genotype not available for this variant".to_string(),
                details: serde_json::Value::Null,
            };
        }
    };

    let father_gt = match &input.father_genotype {
        Some(gt) => gt,
        None => {
            return EvidenceCriterion {
                code: "PS2".to_string(),
                direction: EvidenceDirection::Pathogenic,
                strength: EvidenceStrength::Strong,
                default_strength: EvidenceStrength::Strong,
                met: false,
                evaluated: false,
                summary: "Father genotype not available for this variant".to_string(),
                details: serde_json::Value::Null,
            };
        }
    };

    let min_dp = trio.min_depth;
    let min_gq = trio.min_gq;
    details.insert("min_depth".into(), serde_json::json!(min_dp));
    details.insert("min_gq".into(), serde_json::json!(min_gq));
    details.insert("proband_carries_variant".into(), serde_json::json!(proband_gt.carries_variant()));
    details.insert("mother_hom_ref".into(), serde_json::json!(mother_gt.is_hom_ref));
    details.insert("father_hom_ref".into(), serde_json::json!(father_gt.is_hom_ref));
    details.insert("proband_depth".into(), serde_json::json!(proband_gt.depth));
    details.insert("mother_depth".into(), serde_json::json!(mother_gt.depth));
    details.insert("father_depth".into(), serde_json::json!(father_gt.depth));
    details.insert("proband_gq".into(), serde_json::json!(proband_gt.quality));
    details.insert("mother_gq".into(), serde_json::json!(mother_gt.quality));
    details.insert("father_gq".into(), serde_json::json!(father_gt.quality));

    let proband_carries = proband_gt.carries_variant();
    let mother_ref = mother_gt.is_hom_ref;
    let father_ref = father_gt.is_hom_ref;
    let proband_qc = proband_gt.passes_quality(min_dp, min_gq);
    let mother_qc = mother_gt.passes_quality(min_dp, min_gq);
    let father_qc = father_gt.passes_quality(min_dp, min_gq);

    if !proband_carries {
        return EvidenceCriterion {
            code: "PS2".to_string(),
            direction: EvidenceDirection::Pathogenic,
            strength: EvidenceStrength::Strong,
            default_strength: EvidenceStrength::Strong,
            met: false,
            evaluated: true,
            summary: "Proband does not carry the variant allele".to_string(),
            details: serde_json::Value::Object(details),
        };
    }

    if !proband_qc || !mother_qc || !father_qc {
        let mut fail_reasons = Vec::new();
        if !proband_qc { fail_reasons.push("proband"); }
        if !mother_qc { fail_reasons.push("mother"); }
        if !father_qc { fail_reasons.push("father"); }
        return EvidenceCriterion {
            code: "PS2".to_string(),
            direction: EvidenceDirection::Pathogenic,
            strength: EvidenceStrength::Strong,
            default_strength: EvidenceStrength::Strong,
            met: false,
            evaluated: true,
            summary: format!(
                "Genotype quality insufficient for {}: requires DP>={} and GQ>={}",
                fail_reasons.join(", "), min_dp, min_gq
            ),
            details: serde_json::Value::Object(details),
        };
    }

    let is_de_novo = proband_carries && mother_ref && father_ref;
    details.insert("is_de_novo".into(), serde_json::json!(is_de_novo));

    let summary = if is_de_novo {
        "De novo variant: proband carries variant, both parents homozygous reference, all pass quality thresholds".to_string()
    } else {
        let mut reasons = Vec::new();
        if !mother_ref { reasons.push("mother is not hom_ref"); }
        if !father_ref { reasons.push("father is not hom_ref"); }
        format!("Not de novo: {}", reasons.join(", "))
    };

    EvidenceCriterion {
        code: "PS2".to_string(),
        direction: EvidenceDirection::Pathogenic,
        strength: EvidenceStrength::Strong,
        default_strength: EvidenceStrength::Strong,
        met: is_de_novo,
        evaluated: true,
        summary,
        details: serde_json::Value::Object(details),
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
        summary: "Requires curated functional study evidence (in vitro/in vivo assays) — not automatable from variant data".to_string(),
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
            clinvar_protein: None,
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
