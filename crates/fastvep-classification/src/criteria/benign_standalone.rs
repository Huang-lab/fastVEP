use crate::config::{AcmgConfig, Ba1Exception};
use crate::sa_extract::ClassificationInput;
use crate::types::{EvidenceCriterion, EvidenceDirection, EvidenceStrength};

/// BA1: Allele frequency is >5% in any general continental population dataset.
///
/// Per ClinGen SVI updated recommendation (Ghosh et al. 2018, Hum Mutat),
/// BA1 must NOT be applied to a defined exception list of well-known
/// high-AF variants whose pathogenicity is established (e.g. HFE c.845G>A
/// p.Cys282Tyr — hereditary hemochromatosis; F5 / GJB2 founder alleles).
/// The exception list is configurable via `config.ba1_exceptions` so VCEPs
/// can extend it.
pub fn evaluate_ba1(
    input: &ClassificationInput,
    config: &AcmgConfig,
) -> EvidenceCriterion {
    let mut details = serde_json::Map::new();
    details.insert(
        "af_threshold".into(),
        serde_json::json!(config.ba1_af_threshold),
    );

    // Frequencies from a homology-confounded gene (or for a variant ClinVar
    // labels low-penetrance) cannot support a standalone benign call either.
    if let Some(reason) = super::frequency_gate::benign_blocker(input, config) {
        details.insert("frequency_blocked".into(), serde_json::json!(reason.clone()));
        return EvidenceCriterion {
            code: "BA1".to_string(),
            direction: EvidenceDirection::Benign,
            strength: EvidenceStrength::Standalone,
            default_strength: EvidenceStrength::Standalone,
            met: false,
            evaluated: false,
            summary: format!("BA1 not evaluated: {}", reason),
            details: serde_json::Value::Object(details),
        };
    }

    // Detect whether this allele is on the BA1 exception list. We can only
    // match when both gene_symbol and hgvs_c are populated.
    let exception_match: Option<&Ba1Exception> =
        match (input.gene_symbol.as_deref(), input.hgvs_c.as_deref()) {
            (Some(g), Some(h)) => config
                .ba1_exceptions
                .iter()
                .find(|e| e.gene.eq_ignore_ascii_case(g) && e.hgvs_c.eq_ignore_ascii_case(h)),
            _ => None,
        };

    if let Some(exc) = exception_match {
        details.insert("ba1_exception_match".into(), serde_json::json!(true));
        details.insert("ba1_exception_gene".into(), serde_json::json!(exc.gene));
        details.insert("ba1_exception_hgvs_c".into(), serde_json::json!(exc.hgvs_c));
        if let Some(reason) = &exc.reason {
            details.insert("ba1_exception_reason".into(), serde_json::json!(reason));
        }
        return EvidenceCriterion {
            code: "BA1".to_string(),
            direction: EvidenceDirection::Benign,
            strength: EvidenceStrength::Standalone,
            default_strength: EvidenceStrength::Standalone,
            met: false,
            evaluated: true,
            summary: format!(
                "{} {} is on the ClinGen BA1 exception list — BA1 cannot fire ({})",
                exc.gene,
                exc.hgvs_c,
                exc.reason.as_deref().unwrap_or("Ghosh 2018")
            ),
            details: serde_json::Value::Object(details),
        };
    }

    let (met, summary) = if let Some(ref gnomad) = input.gnomad {
        // ClinGen SVI gnomAD v4 guidance (March 2024): require minimum AN
        // before frequency-based criteria fire — protects against noisy AF
        // estimates in poorly-covered populations. Treat missing AN the
        // same as below-minimum: the SVI guidance is a *requirement*, not
        // an opt-in check, so we cannot fire BA1 without confirming AN.
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
                            "BA1 not evaluated: gnomAD AN={} below minimum {} (gnomAD v4 guidance)",
                            an, config.min_an_for_frequency_criteria
                        )
                    }
                    None => {
                        details.insert("an_missing".into(), serde_json::json!(true));
                        format!(
                            "BA1 not evaluated: gnomAD AN unavailable; minimum {} required (gnomAD v4 guidance)",
                            config.min_an_for_frequency_criteria
                        )
                    }
                };
                return EvidenceCriterion {
                    code: "BA1".to_string(),
                    direction: EvidenceDirection::Benign,
                    strength: EvidenceStrength::Standalone,
                    default_strength: EvidenceStrength::Standalone,
                    met: false,
                    evaluated: false,
                    summary,
                    details: serde_json::Value::Object(details),
                };
            }
        }
        let test_af = super::frequency_gate::benign_test_af(gnomad, config);
        if let Some((af, af_label)) = test_af {
            details.insert("max_pop_af".into(), serde_json::json!(af));
            details.insert("af_statistic".into(), serde_json::json!(af_label));
            if let Some(point) = gnomad.max_pop_af() {
                details.insert("max_pop_point_af".into(), serde_json::json!(point));
            }

            // Add per-population breakdown for transparency
            let mut pop_afs = serde_json::Map::new();
            if let Some(v) = gnomad.all_af { pop_afs.insert("all".into(), serde_json::json!(v)); }
            if let Some(v) = gnomad.afr_af { pop_afs.insert("afr".into(), serde_json::json!(v)); }
            if let Some(v) = gnomad.nfe_af { pop_afs.insert("nfe".into(), serde_json::json!(v)); }
            if let Some(v) = gnomad.eas_af { pop_afs.insert("eas".into(), serde_json::json!(v)); }
            if let Some(v) = gnomad.amr_af { pop_afs.insert("amr".into(), serde_json::json!(v)); }
            if let Some(v) = gnomad.asj_af { pop_afs.insert("asj".into(), serde_json::json!(v)); }
            if let Some(v) = gnomad.fin_af { pop_afs.insert("fin".into(), serde_json::json!(v)); }
            if let Some(v) = gnomad.sas_af { pop_afs.insert("sas".into(), serde_json::json!(v)); }
            details.insert("population_afs".into(), serde_json::Value::Object(pop_afs));

            if af > config.ba1_af_threshold {
                (
                    true,
                    format!(
                        "Common variant: {}={:.4} exceeds {:.2} threshold",
                        af_label, af, config.ba1_af_threshold
                    ),
                )
            } else {
                (
                    false,
                    format!(
                        "{}={:.6} does not exceed {:.2} threshold",
                        af_label, af, config.ba1_af_threshold
                    ),
                )
            }
        } else {
            (false, "gnomAD data present but no allele frequency values".to_string())
        }
    } else {
        (false, "No gnomAD population frequency data available".to_string())
    };

    EvidenceCriterion {
        code: "BA1".to_string(),
        direction: EvidenceDirection::Benign,
        strength: EvidenceStrength::Standalone,
        default_strength: EvidenceStrength::Standalone,
        met,
        evaluated: input.gnomad.is_some(),
        summary,
        details: serde_json::Value::Object(details),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_support::minimal_input;
    use crate::sa_extract::GnomadData;

    #[test]
    fn test_ba1_common_variant() {
        let input = ClassificationInput {
            impact: fastvep_core::Impact::Modifier,
            gene_symbol: None,
            is_canonical: false,
            gnomad: Some(GnomadData {
                all_af: Some(0.10),
                afr_af: Some(0.15),
                all_an: Some(100_000),
                ..Default::default()
            }),
            ..minimal_input()
        };
        let result = evaluate_ba1(&input, &AcmgConfig::default());
        assert!(result.met);
        assert_eq!(result.strength, EvidenceStrength::Standalone);
    }

    #[test]
    fn test_ba1_rare_variant() {
        let input = ClassificationInput {
            impact: fastvep_core::Impact::Modifier,
            gene_symbol: None,
            is_canonical: false,
            gnomad: Some(GnomadData {
                all_af: Some(0.001),
                ..Default::default()
            }),
            ..minimal_input()
        };
        let result = evaluate_ba1(&input, &AcmgConfig::default());
        assert!(!result.met);
    }

    #[test]
    fn test_ba1_exception_list_blocks_known_pathogenic_high_af() {
        // HFE c.845G>A (p.Cys282Tyr) — hereditary hemochromatosis. ~10% AF in
        // European populations but pathogenic. Per Ghosh 2018, BA1 must NOT fire.
        let input = ClassificationInput {
            impact: fastvep_core::Impact::Modifier,
            gene_symbol: Some("HFE".to_string()),
            gnomad: Some(GnomadData {
                all_af: Some(0.06),
                nfe_af: Some(0.10),
                ..Default::default()
            }),
            hgvs_c: Some("c.845G>A".to_string()),
            ..minimal_input()
        };
        let result = evaluate_ba1(&input, &AcmgConfig::default());
        assert!(!result.met, "BA1 must not fire for HFE c.845G>A (Ghosh 2018 exception)");
        assert!(result.evaluated);
        assert!(result.summary.contains("exception"));
    }

    #[test]
    fn test_ba1_exception_match_is_case_insensitive() {
        // Pipeline may emit "C.845G>a" or other casing — match must still work.
        let input = ClassificationInput {
            impact: fastvep_core::Impact::Modifier,
            gene_symbol: Some("hfe".to_string()),
            gnomad: Some(GnomadData { all_af: Some(0.10), ..Default::default() }),
            hgvs_c: Some("C.845G>A".to_string()),
            ..minimal_input()
        };
        let result = evaluate_ba1(&input, &AcmgConfig::default());
        assert!(!result.met);
    }

    #[test]
    fn test_ba1_high_af_non_exception_still_fires() {
        // Same gene (HFE) but a different c. notation NOT on the exception
        // list → BA1 still fires at high AF.
        let input = ClassificationInput {
            impact: fastvep_core::Impact::Modifier,
            gene_symbol: Some("HFE".to_string()),
            gnomad: Some(GnomadData { all_af: Some(0.10), all_an: Some(100_000), ..Default::default() }),
            hgvs_c: Some("c.999A>T".to_string()),
            ..minimal_input()
        };
        let result = evaluate_ba1(&input, &AcmgConfig::default());
        assert!(result.met);
    }

    #[test]
    fn test_ba1_low_an_not_evaluated() {
        // PR10 (gnomAD v4 guidance): AN below 2000 → NotEvaluated, not Benign,
        // even at high AF — the frequency estimate is unreliable.
        let input = ClassificationInput {
            impact: fastvep_core::Impact::Modifier,
            gene_symbol: None,
            is_canonical: false,
            gnomad: Some(GnomadData {
                all_af: Some(0.10),
                all_an: Some(500), // way below 2000
                ..Default::default()
            }),
            ..minimal_input()
        };
        let result = evaluate_ba1(&input, &AcmgConfig::default());
        assert!(!result.met);
        assert!(!result.evaluated);
        assert!(result.summary.contains("below minimum"));
    }

    #[test]
    fn test_ba1_one_pop_above_threshold() {
        let input = ClassificationInput {
            impact: fastvep_core::Impact::Modifier,
            gene_symbol: None,
            is_canonical: false,
            gnomad: Some(GnomadData {
                all_af: Some(0.02),
                eas_af: Some(0.06), // Only EAS above 5%
                all_an: Some(100_000),
                ..Default::default()
            }),
            ..minimal_input()
        };
        let result = evaluate_ba1(&input, &AcmgConfig::default());
        assert!(result.met);
    }
}
