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

/// A reason why frequency-based benign evidence cannot be assessed for this
/// variant at all. Shared by BS1 and BS2 (and mirrored by BA1 / PM2) so that a
/// gene or variant whose population data is untrustworthy reports
/// NotEvaluated rather than firing on the untrustworthy number.
pub(crate) fn frequency_evidence_blocker(
    input: &ClassificationInput,
    config: &AcmgConfig,
) -> Option<String> {
    if config.is_homology_unreliable(input.gene_symbol.as_deref()) {
        return Some(format!(
            "gene {} has paralogue/pseudogene homology that makes population frequencies unreliable (Mandelker 2016, PMID 27228465)",
            input.gene_symbol.as_deref().unwrap_or("?")
        ));
    }

    if config.clinvar_low_penetrance_blocks_benign_frequency {
        if let Some(ref cv) = input.clinvar {
            if cv.review_stars() >= 2 {
                if let Some(term) = cv.low_penetrance_term() {
                    return Some(format!(
                        "ClinVar ({} stars) reports this variant as \"{}\"; a low-penetrance or risk allele is expected to be frequent and is outside BS2's full-penetrance precondition",
                        cv.review_stars(),
                        term
                    ));
                }
            }
        }
    }

    None
}

/// 95 % lower confidence bound on a Poisson rate that produced `observed`
/// events, via the chi-square relation `lower = 0.5 * chi2_{0.05, 2k}`.
/// Tabulated for the small counts that matter here and approximated with the
/// Wilson-Hilferty expansion above the table.
fn poisson_lower_95(observed: u64) -> f64 {
    const TABLE: [f64; 21] = [
        0.0, 0.0513, 0.3554, 0.8177, 1.3663, 1.9702, 2.6130, 3.2853, 3.9808,
        4.6952, 5.4254, 6.1690, 6.9242, 7.6896, 8.4639, 9.2463, 10.0360,
        10.8324, 11.6350, 12.4432, 13.2547,
    ];
    if (observed as usize) < TABLE.len() {
        return TABLE[observed as usize];
    }
    let k = observed as f64;
    // Wilson-Hilferty: lower ≈ k * (1 - 1/(9k) - 1.645/(3*sqrt(k)))^3
    let term = 1.0 - 1.0 / (9.0 * k) - 1.645 / (3.0 * k.sqrt());
    (k * term * term * term).max(0.0)
}

/// BS1: Allele frequency is greater than expected for disorder.
fn evaluate_bs1(
    input: &ClassificationInput,
    config: &AcmgConfig,
) -> EvidenceCriterion {
    let threshold = config.effective_bs1_threshold(input.gene_symbol.as_deref());

    let mut details = serde_json::Map::new();
    details.insert("af_threshold".into(), serde_json::json!(threshold));

    if let Some(reason) = frequency_evidence_blocker(input, config) {
        details.insert("frequency_blocked".into(), serde_json::json!(reason.clone()));
        return EvidenceCriterion {
            code: "BS1".to_string(),
            direction: EvidenceDirection::Benign,
            strength: EvidenceStrength::Strong,
            default_strength: EvidenceStrength::Strong,
            met: false,
            evaluated: false,
            summary: format!("BS1 not evaluated: {}", reason),
            details: serde_json::Value::Object(details),
        };
    }

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
        // ClinGen SVI guidance applies BS1 against the **max-population
        // AF** (mirroring BA1), not the cohort-wide allAf. Using the
        // cohort AF would let a 5% variant in a single subpopulation slip
        // under a 1 % BS1 threshold whenever the global cohort happens
        // to dilute it. METHODS.md aligns: "max population AF" for both
        // BA1 and BS1.
        let max_pop_af = gnomad.max_pop_af().unwrap_or(0.0);
        let cohort_af = gnomad.all_af.unwrap_or(0.0);
        details.insert("gnomad_allAf".into(), serde_json::json!(cohort_af));
        details.insert("gnomad_max_pop_af".into(), serde_json::json!(max_pop_af));

        // BS1 should not fire if BA1 would fire (BA1 takes precedence)
        if max_pop_af > config.ba1_af_threshold {
            (
                false,
                format!(
                    "BA1 takes precedence (max pop AF={:.4} > BA1 threshold {:.2})",
                    max_pop_af, config.ba1_af_threshold
                ),
            )
        } else if max_pop_af > threshold {
            (
                true,
                format!(
                    "Max-pop AF ({:.6}) exceeds expected for disorder (threshold={:.4})",
                    max_pop_af, threshold
                ),
            )
        } else {
            (
                false,
                format!(
                    "Max-pop AF ({:.6}) within expected range (threshold={:.4})",
                    max_pop_af, threshold
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
/// dominant (heterozygous), or X-linked (hemizygous) disorder with full
/// penetrance expected at an early age (Richards 2015).
///
/// - **Recessive** (or unknown inheritance): require ≥1 homozygote in gnomAD.
/// - **Dominant**: require AC ≥ `bs2_ad_min_ac` (default 5) — Richards 2015
///   says "observed in unaffected adult", which is not the same as a single
///   carrier of a novel allele in a 100K cohort. Singletons / doubletons
///   are sequencing-noise plausibility, not evidence the variant is
///   tolerated. ClinGen VCEPs commonly use AC ≥ 5 (Hereditary Cancer
///   VCEP, Lynch Syndrome curation guide).
///
/// Inheritance is inferred from the disease-gene `.oga` (ClinGen GDV
/// preferred, OMIM accepted as legacy).
fn evaluate_bs2(
    input: &ClassificationInput,
    config: &AcmgConfig,
) -> EvidenceCriterion {
    let mut details = serde_json::Map::new();

    if let Some(reason) = frequency_evidence_blocker(input, config) {
        details.insert("frequency_blocked".into(), serde_json::json!(reason.clone()));
        return EvidenceCriterion {
            code: "BS2".to_string(),
            direction: EvidenceDirection::Benign,
            strength: EvidenceStrength::Strong,
            default_strength: EvidenceStrength::Strong,
            met: false,
            evaluated: false,
            summary: format!("BS2 not evaluated: {}", reason),
            details: serde_json::Value::Object(details),
        };
    }

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
        let ac = gnomad.all_ac.unwrap_or(0);
        details.insert("gnomad_allHc".into(), serde_json::json!(hc));
        details.insert("gnomad_allAn".into(), serde_json::json!(an));
        details.insert("gnomad_allAc".into(), serde_json::json!(ac));

        if is_dominant && !is_recessive && ac >= config.bs2_ad_min_ac {
            // For AD-only genes: ≥`bs2_ad_min_ac` heterozygote observations
            // in healthy adults (gnomAD). Singletons / doubletons of a
            // novel allele are not BS2 evidence.
            (
                true,
                true,
                format!(
                    "Observed in gnomAD as ≥{} unaffected heterozygotes (AC={}) for autosomal-dominant disorder",
                    config.bs2_ad_min_ac, ac
                ),
            )
        } else if hc > 0 {
            // Recessive / X-linked / unknown inheritance. Richards 2015 asks
            // for observation "in a healthy adult ... with full penetrance
            // expected at an early age", so the question is not "is there a
            // homozygote?" but "are there more homozygotes than the disorder
            // itself could account for?". Compare the 95 % lower bound on the
            // homozygote frequency against the maximum credible disease
            // prevalence: that makes the criterion scale with the size of the
            // cohort behind the observation, which is what a single
            // homozygote in gnomAD v4 (730 K individuals) fails.
            let individuals = (an / 2).max(1);
            let hom_freq_lower_95 = poisson_lower_95(hc) / individuals as f64;
            details.insert("gnomad_individuals".into(), serde_json::json!(individuals));
            details.insert(
                "hom_freq_lower_95".into(),
                serde_json::json!(hom_freq_lower_95),
            );
            details.insert(
                "prevalence_threshold".into(),
                serde_json::json!(config.bs2_hom_prevalence_threshold),
            );
            details.insert("bs2_ar_min_hom".into(), serde_json::json!(config.bs2_ar_min_hom));

            if hc < config.bs2_ar_min_hom {
                (
                    false,
                    true,
                    format!(
                        "{} homozygote(s) in gnomAD is below the floor of {}; too few to establish tolerance in healthy adults",
                        hc, config.bs2_ar_min_hom
                    ),
                )
            } else if hom_freq_lower_95 > config.bs2_hom_prevalence_threshold {
                (
                    true,
                    true,
                    format!(
                        "{} homozygotes among {} gnomAD individuals (95% lower bound on homozygote frequency {:.2e} exceeds the {:.0e} prevalence bar), so homozygosity is tolerated in healthy adults",
                        hc, individuals, hom_freq_lower_95, config.bs2_hom_prevalence_threshold
                    ),
                )
            } else {
                (
                    false,
                    true,
                    format!(
                        "{} homozygotes among {} gnomAD individuals gives a 95% lower bound of {:.2e}, which does not exceed the {:.0e} prevalence bar; consistent with a late-onset, reduced-penetrance or variably expressive disorder rather than tolerance",
                        hc, individuals, hom_freq_lower_95, config.bs2_hom_prevalence_threshold
                    ),
                )
            }
        } else if is_dominant && !is_recessive {
            (
                false,
                true,
                format!(
                    "AC={} below BS2 threshold ({} required for AD)",
                    ac, config.bs2_ad_min_ac
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
            is_pure_insertion: None,
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

    fn make_input_omim(
        gnomad: Option<GnomadData>,
        omim: Option<crate::sa_extract::OmimData>,
    ) -> ClassificationInput {
        let mut i = make_input(gnomad);
        i.omim = omim;
        i
    }

    #[test]
    fn test_bs2_ad_gene_singleton_does_not_fire() {
        // AD gene + AC=1: a single heterozygote is not "observed in
        // healthy adult" per Richards 2015 BS2; default threshold AC≥5.
        use crate::sa_extract::OmimData;
        let input = make_input_omim(
            Some(GnomadData {
                all_ac: Some(1),
                all_hc: Some(0),
                all_af: Some(0.000001),
                all_an: Some(1_000_000),
                ..Default::default()
            }),
            Some(OmimData {
                mim_number: None,
                phenotypes: Some(vec!["dominant disorder".into()]),
            }),
        );
        let r = evaluate_bs2(&input, &AcmgConfig::default());
        assert!(!r.met, "AC=1 should not fire AD BS2; got {}", r.summary);
        assert!(r.summary.contains("below BS2 threshold"));
    }

    #[test]
    fn test_bs2_ad_gene_meets_threshold_fires() {
        // AD gene + AC≥5 (default `bs2_ad_min_ac`) → BS2 fires.
        use crate::sa_extract::OmimData;
        let input = make_input_omim(
            Some(GnomadData {
                all_ac: Some(7),
                all_hc: Some(0),
                all_af: Some(7e-6),
                all_an: Some(1_000_000),
                ..Default::default()
            }),
            Some(OmimData {
                mim_number: None,
                phenotypes: Some(vec!["dominant disorder".into()]),
            }),
        );
        let r = evaluate_bs2(&input, &AcmgConfig::default());
        assert!(r.met, "AC=7 ≥ default 5 should fire AD BS2");
        assert!(r.summary.contains("autosomal-dominant"));
    }

    #[test]
    fn test_bs2_ad_gene_min_ac_configurable() {
        // Config knob lets a stricter VCEP raise the threshold.
        use crate::sa_extract::OmimData;
        let mut cfg = AcmgConfig::default();
        cfg.bs2_ad_min_ac = 20;
        let input = make_input_omim(
            Some(GnomadData {
                all_ac: Some(7),
                all_hc: Some(0),
                ..Default::default()
            }),
            Some(OmimData {
                mim_number: None,
                phenotypes: Some(vec!["dominant disorder".into()]),
            }),
        );
        let r = evaluate_bs2(&input, &cfg);
        assert!(!r.met, "AC=7 < raised threshold 20 should not fire");
    }

    fn recessive_input(hc: u64, an: u64) -> ClassificationInput {
        use crate::sa_extract::OmimData;
        make_input_omim(
            Some(GnomadData {
                all_ac: Some(hc * 2 + 10),
                all_hc: Some(hc),
                all_an: Some(an),
                ..Default::default()
            }),
            Some(OmimData {
                mim_number: None,
                phenotypes: Some(vec!["recessive disorder".into()]),
            }),
        )
    }

    #[test]
    fn test_bs2_ar_single_homozygote_in_large_cohort_does_not_fire() {
        // gnomAD v4 scale (730 K individuals). One homozygote is what you
        // expect for a late-onset or reduced-penetrance recessive disorder, so
        // it must not be read as tolerance. This was the single largest source
        // of false-benign calls in the round-2 medical-genetics review.
        let r = evaluate_bs2(&recessive_input(1, 1_460_000), &AcmgConfig::default());
        assert!(!r.met, "1 homozygote among 730 K individuals must not fire BS2");
    }

    #[test]
    fn test_bs2_ar_many_homozygotes_in_large_cohort_fires() {
        // 30 homozygotes among 730 K individuals puts the 95 % lower bound on
        // the homozygote frequency above the 1e-5 prevalence bar.
        let r = evaluate_bs2(&recessive_input(30, 1_460_000), &AcmgConfig::default());
        assert!(r.met, "30 homozygotes among 730 K individuals should fire BS2");
        assert!(r.summary.contains("tolerated in healthy adults"));
    }

    #[test]
    fn test_bs2_ar_scales_with_cohort_size() {
        // The same raw count means different things at different cohort sizes.
        // 3 homozygotes among 2 K individuals is a homozygote frequency far
        // above any Mendelian prevalence; 3 among 730 K is not.
        let small = evaluate_bs2(&recessive_input(3, 4_000), &AcmgConfig::default());
        let large = evaluate_bs2(&recessive_input(3, 1_460_000), &AcmgConfig::default());
        assert!(small.met, "3 homozygotes among 2 K individuals should fire BS2");
        assert!(
            !large.met,
            "3 homozygotes among 730 K individuals should not fire BS2"
        );
    }

    #[test]
    fn test_bs2_homology_gene_not_evaluated() {
        // CYP21A2 frequencies are confounded by CYP21A1P (Mandelker 2016).
        let mut input = recessive_input(30, 1_460_000);
        input.gene_symbol = Some("CYP21A2".to_string());
        let r = evaluate_bs2(&input, &AcmgConfig::default());
        assert!(!r.met);
        assert!(!r.evaluated, "BS2 should be NotEvaluated, not a negative call");
        assert!(r.summary.contains("homology"));
    }

    #[test]
    fn test_bs1_low_penetrance_clinvar_term_not_evaluated() {
        use crate::sa_extract::ClinvarData;
        let mut input = make_input(Some(GnomadData {
            all_af: Some(0.037),
            all_an: Some(1_460_000),
            ..Default::default()
        }));
        input.clinvar = Some(ClinvarData {
            significance: Some(vec!["Pathogenic,_low_penetrance".into()]),
            review_status: Some("criteria_provided,_multiple_submitters,_no_conflicts".into()),
            ..Default::default()
        });
        let r = evaluate_bs1(&input, &AcmgConfig::default());
        assert!(!r.met);
        assert!(!r.evaluated);
        assert!(r.summary.contains("low_penetrance") || r.summary.contains("low penetrance"));
    }

    // ── BS1 (max-pop AF) ──

    #[test]
    fn test_bs1_uses_max_pop_af_not_cohort_af() {
        // Cohort AF below threshold but max-pop AF above: BS1 fires.
        // (Pre-fix this would have used `all_af` and missed the
        // single-population enrichment.)
        let input = make_input(Some(GnomadData {
            all_af: Some(0.001),
            // 5 % in EAS — well above default BS1 threshold of 1 %.
            eas_af: Some(0.05),
            all_an: Some(2_000_000),
            ..Default::default()
        }));
        let r = evaluate_bs1(&input, &AcmgConfig::default());
        assert!(r.met, "max-pop AF should drive BS1, not cohort AF");
        assert!(r.summary.contains("Max-pop"));
    }
}
