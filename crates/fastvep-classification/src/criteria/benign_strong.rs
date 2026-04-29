use crate::config::AcmgConfig;
use crate::sa_extract::ClassificationInput;
use crate::types::{EvidenceCriterion, EvidenceDirection, EvidenceStrength};

/// Evaluate all benign strong criteria: BS1, BS2, BS3, BS4.
pub fn evaluate_all(
    input: &ClassificationInput,
    config: &AcmgConfig,
) -> Vec<EvidenceCriterion> {
    vec![
        evaluate_bs1(input, config),
        evaluate_bs2(input, config),
        evaluate_bs3(input, config),
        evaluate_bs4(input, config),
    ]
}

/// BS1: Allele frequency is greater than expected for disorder.
fn evaluate_bs1(
    input: &ClassificationInput,
    config: &AcmgConfig,
) -> EvidenceCriterion {
    let threshold = config.effective_bs1_threshold(input.gene_symbol.as_deref());

    let mut details = serde_json::Map::new();
    details.insert("af_threshold".into(), serde_json::json!(threshold));

    let (met, summary) = if let Some(ref gnomad) = input.gnomad {
        // ClinGen SVI gnomAD v4 guidance (March 2024): require minimum AN
        // before BS1 fires, same as BA1. Treat missing AN as NotEvaluated
        // (the SVI guidance is a requirement, not an opt-in).
        match gnomad.all_an {
            Some(an) if an >= config.min_an_for_frequency_criteria => {}
            other => {
                details.insert(
                    "min_an_for_frequency_criteria".into(),
                    serde_json::json!(config.min_an_for_frequency_criteria),
                );
                let summary = match other {
                    Some(an) => {
                        details.insert("an_below_minimum".into(), serde_json::json!(an));
                        format!(
                            "BS1 not evaluated: gnomAD AN={} below minimum {} (gnomAD v4 guidance)",
                            an, config.min_an_for_frequency_criteria
                        )
                    }
                    None => {
                        details.insert("an_missing".into(), serde_json::json!(true));
                        format!(
                            "BS1 not evaluated: gnomAD AN unavailable; minimum {} required (gnomAD v4 guidance)",
                            config.min_an_for_frequency_criteria
                        )
                    }
                };
                return EvidenceCriterion {
                    code: "BS1".to_string(),
                    direction: EvidenceDirection::Benign,
                    strength: EvidenceStrength::Strong,
                    default_strength: EvidenceStrength::Strong,
                    met: false,
                    evaluated: false,
                    summary,
                    details: serde_json::Value::Object(details),
                };
            }
        }
        let af = gnomad.all_af.unwrap_or(0.0);
        details.insert("gnomad_allAf".into(), serde_json::json!(af));

        // BS1 should not fire if BA1 would fire (BA1 takes precedence)
        let max_pop_af = gnomad.max_pop_af().unwrap_or(0.0);
        if max_pop_af > config.ba1_af_threshold {
            (
                false,
                format!(
                    "BA1 takes precedence (max pop AF={:.4} > BA1 threshold {:.2})",
                    max_pop_af, config.ba1_af_threshold
                ),
            )
        } else if af > threshold {
            (
                true,
                format!(
                    "Allele frequency ({:.6}) exceeds expected for disorder (threshold={:.4})",
                    af, threshold
                ),
            )
        } else {
            (
                false,
                format!(
                    "Allele frequency ({:.6}) within expected range (threshold={:.4})",
                    af, threshold
                ),
            )
        }
    } else {
        (false, "No gnomAD data available".to_string())
    };

    EvidenceCriterion {
        code: "BS1".to_string(),
        direction: EvidenceDirection::Benign,
        strength: EvidenceStrength::Strong,
        default_strength: EvidenceStrength::Strong,
        met,
        evaluated: input.gnomad.is_some(),
        summary,
        details: serde_json::Value::Object(details),
    }
}

/// BS2: Observed in a healthy adult individual for a recessive (homozygous),
/// dominant (heterozygous), or X-linked (hemizygous) disorder with full penetrance
/// expected at an early age.
///
/// For dominant disorders: presence in gnomAD at any frequency with adequate sample size.
/// For recessive disorders: homozygotes in gnomAD.
/// Inheritance pattern inferred from OMIM phenotype annotations when available.
fn evaluate_bs2(
    input: &ClassificationInput,
    _config: &AcmgConfig,
) -> EvidenceCriterion {
    let mut details = serde_json::Map::new();

    let is_dominant = input
        .omim
        .as_ref()
        .map_or(false, |o| o.has_dominant_inheritance());
    let is_recessive = input
        .omim
        .as_ref()
        .map_or(false, |o| o.has_recessive_inheritance());
    details.insert("omim_dominant".into(), serde_json::json!(is_dominant));
    details.insert("omim_recessive".into(), serde_json::json!(is_recessive));

    let (met, evaluated, summary) = if let Some(ref gnomad) = input.gnomad {
        let hc = gnomad.all_hc.unwrap_or(0);
        let an = gnomad.all_an.unwrap_or(0);
        let af = gnomad.all_af.unwrap_or(0.0);
        details.insert("gnomad_allHc".into(), serde_json::json!(hc));
        details.insert("gnomad_allAn".into(), serde_json::json!(an));

        if is_dominant && af > 0.0 && an > 10000 {
            // For dominant fully-penetrant disorders, any carrier in a large population
            // suggests benign (the individual would be affected if pathogenic)
            (
                true,
                true,
                format!(
                    "Present in gnomAD (AF={:.6}, AN={}) for autosomal dominant disorder, suggesting tolerated",
                    af, an
                ),
            )
        } else if hc > 0 {
            // For recessive or unknown inheritance: homozygotes suggest benign
            (
                true,
                true,
                format!(
                    "Observed as homozygous in gnomAD ({} homozygotes), suggesting tolerated in healthy adults",
                    hc
                ),
            )
        } else {
            (
                false,
                true,
                "No homozygotes observed in gnomAD".to_string(),
            )
        }
    } else {
        (false, false, "No gnomAD data available".to_string())
    };

    EvidenceCriterion {
        code: "BS2".to_string(),
        direction: EvidenceDirection::Benign,
        strength: EvidenceStrength::Strong,
        default_strength: EvidenceStrength::Strong,
        met,
        evaluated,
        summary,
        details: serde_json::Value::Object(details),
    }
}

/// BS3: Well-established in vitro or in vivo functional studies show no damaging effect.
fn evaluate_bs3(
    _input: &ClassificationInput,
    _config: &AcmgConfig,
) -> EvidenceCriterion {
    EvidenceCriterion {
        code: "BS3".to_string(),
        direction: EvidenceDirection::Benign,
        strength: EvidenceStrength::Strong,
        default_strength: EvidenceStrength::Strong,
        met: false,
        evaluated: false,
        summary: "Requires curated functional study evidence showing no damaging effect — not automatable from variant data".to_string(),
        details: serde_json::Value::Null,
    }
}

/// BS4: Lack of segregation in affected members of a family.
fn evaluate_bs4(
    _input: &ClassificationInput,
    _config: &AcmgConfig,
) -> EvidenceCriterion {
    EvidenceCriterion {
        code: "BS4".to_string(),
        direction: EvidenceDirection::Benign,
        strength: EvidenceStrength::Strong,
        default_strength: EvidenceStrength::Strong,
        met: false,
        evaluated: false,
        summary: "Requires multi-generation pedigree with affection status to assess lack of segregation".to_string(),
        details: serde_json::Value::Null,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sa_extract::GnomadData;
    use fastvep_core::Impact;

    fn make_input(gnomad: Option<GnomadData>) -> ClassificationInput {
        ClassificationInput {
            consequences: vec![],
            impact: Impact::Modifier,
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
    fn test_bs1_above_threshold() {
        let input = make_input(Some(GnomadData {
            all_af: Some(0.02),
            all_an: Some(100_000),
            ..Default::default()
        }));
        let result = evaluate_bs1(&input, &AcmgConfig::default());
        assert!(result.met);
    }

    #[test]
    fn test_bs1_below_threshold() {
        let input = make_input(Some(GnomadData {
            all_af: Some(0.001),
            all_an: Some(100_000),
            ..Default::default()
        }));
        let result = evaluate_bs1(&input, &AcmgConfig::default());
        assert!(!result.met);
    }

    #[test]
    fn test_bs1_ba1_takes_precedence() {
        let input = make_input(Some(GnomadData {
            all_af: Some(0.10),
            afr_af: Some(0.10),
            all_an: Some(100_000),
            ..Default::default()
        }));
        let result = evaluate_bs1(&input, &AcmgConfig::default());
        assert!(!result.met); // BA1 would fire, so BS1 should not
    }

    #[test]
    fn test_bs1_low_an_not_evaluated() {
        // gnomAD v4 guidance: AN below 2000 → NotEvaluated, even at high AF.
        let input = make_input(Some(GnomadData {
            all_af: Some(0.02),
            all_an: Some(500),
            ..Default::default()
        }));
        let result = evaluate_bs1(&input, &AcmgConfig::default());
        assert!(!result.met);
        assert!(!result.evaluated);
        assert!(result.summary.contains("below minimum"));
    }

    #[test]
    fn test_bs1_missing_an_not_evaluated() {
        // gnomAD record present but AN is None → NotEvaluated, never fires.
        let input = make_input(Some(GnomadData {
            all_af: Some(0.02),
            all_an: None,
            ..Default::default()
        }));
        let result = evaluate_bs1(&input, &AcmgConfig::default());
        assert!(!result.met);
        assert!(!result.evaluated);
        assert!(result.summary.contains("AN unavailable"));
    }

    #[test]
    fn test_bs2_homozygotes_present() {
        let input = make_input(Some(GnomadData {
            all_hc: Some(5),
            ..Default::default()
        }));
        let result = evaluate_bs2(&input, &AcmgConfig::default());
        assert!(result.met);
    }

    #[test]
    fn test_bs2_no_homozygotes() {
        let input = make_input(Some(GnomadData {
            all_hc: Some(0),
            ..Default::default()
        }));
        let result = evaluate_bs2(&input, &AcmgConfig::default());
        assert!(!result.met);
    }
}
