use super::frequency_gate::benign_blocker;
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

    if let Some(reason) = benign_blocker(input, config) {
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
        // ClinGen SVI guidance applies BS1 across populations rather than to
        // the cohort-wide allAf: using the cohort AF would let a 5 % variant in
        // a single subpopulation slip under a 1 % BS1 threshold whenever the
        // global cohort happens to dilute it. `benign_test_af` supplies that
        // cross-population number - the filtering AF where the database
        // provides one, otherwise the population maximum - and BA1 reads the
        // same helper, so the two cannot disagree about how frequent the
        // variant is.
        let (test_af, af_label) = super::frequency_gate::benign_test_af(gnomad, config)
            .unwrap_or((0.0, "max population AF"));
        let cohort_af = gnomad.all_af.unwrap_or(0.0);
        details.insert("gnomad_allAf".into(), serde_json::json!(cohort_af));
        details.insert("gnomad_max_pop_af".into(), serde_json::json!(test_af));
        details.insert("af_statistic".into(), serde_json::json!(af_label));

        // BS1 should not fire if BA1 would fire (BA1 takes precedence)
        if test_af > config.ba1_af_threshold {
            (
                false,
                format!(
                    "BA1 takes precedence ({}={:.4} > BA1 threshold {:.2})",
                    af_label, test_af, config.ba1_af_threshold
                ),
            )
        } else if test_af > threshold {
            (
                true,
                format!(
                    "{} ({:.6}) exceeds expected for disorder (threshold={:.4})",
                    af_label, test_af, threshold
                ),
            )
        } else {
            (
                false,
                format!(
                    "{} ({:.6}) within expected range (threshold={:.4})",
                    af_label, test_af, threshold
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
/// - **Recessive / X-linked** (or unknown inheritance): count the individuals
///   with no functional copy - homozygotes plus, on a non-PAR sex-chromosome
///   site, hemizygotes - and require both an absolute floor of
///   `bs2_ar_min_hom` and a 95 % lower bound on their frequency above
///   `bs2_hom_prevalence_threshold`. The criterion therefore scales with the
///   size of the cohort behind the observation: one homozygote among gnomAD's
///   730 K individuals is what a late-onset or reduced-penetrance disorder
///   looks like, not evidence of tolerance.
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

    if let Some(reason) = benign_blocker(input, config) {
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
        .is_some_and(|o| o.has_dominant_inheritance());
    let is_recessive = input
        .omim
        .as_ref()
        .is_some_and(|o| o.has_recessive_inheritance());
    details.insert("omim_dominant".into(), serde_json::json!(is_dominant));
    details.insert("omim_recessive".into(), serde_json::json!(is_recessive));

    let (met, evaluated, summary) = if let Some(ref gnomad) = input.gnomad {
        let hc = gnomad.all_hc.unwrap_or(0);
        let an = gnomad.all_an.unwrap_or(0);
        let ac = gnomad.all_ac.unwrap_or(0);
        details.insert("gnomad_allHc".into(), serde_json::json!(hc));
        details.insert("gnomad_allAn".into(), serde_json::json!(an));
        details.insert("gnomad_allAc".into(), serde_json::json!(ac));

        // Hemizygous observations, available only for a non-PAR sex-chromosome
        // site in a database built with the XY columns. Everywhere else these
        // are 0 and every expression below collapses to the autosomal form.
        let hemizygotes = gnomad.hemizygote_count().unwrap_or(0);
        let hemizygote_an = gnomad.hemizygote_individuals().unwrap_or(0);
        if hemizygotes > 0 || hemizygote_an > 0 {
            details.insert("gnomad_hemizygotes".into(), serde_json::json!(hemizygotes));
            details.insert("gnomad_xy_individuals".into(), serde_json::json!(hemizygote_an));
        }

        // Individuals surveyed. XX samples contribute two alleles each and XY
        // samples on a non-PAR site contribute one, so the cohort size is
        // `(AN - AN_XY) / 2 + AN_XY`. With no XY data this is exactly `AN / 2`,
        // the autosomal count used before hemizygotes were extracted.
        let individuals = (an.saturating_sub(hemizygote_an) / 2 + hemizygote_an).max(1);

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
        } else if hc + hemizygotes > 0 {
            // Recessive / X-linked / unknown inheritance. Richards 2015 asks
            // for observation "in a healthy adult ... with full penetrance
            // expected at an early age", so the question is not "is there a
            // homozygote?" but "are there more individuals lacking a working
            // copy than the disorder itself could account for?". Compare the
            // 95 % lower bound on their frequency against the maximum credible
            // disease prevalence: that makes the criterion scale with the size
            // of the cohort behind the observation, which is what a single
            // homozygote in gnomAD v4 (730 K individuals) fails.
            //
            // Homozygotes and hemizygotes are counted together because they are
            // the same observation genetically: an individual with no
            // functional allele. Richards 2015 lists both. gnomAD calls XY
            // samples haploid outside the pseudoautosomal regions, so an
            // X-linked hemizygote appears in `AC_XY` and in *neither*
            // `nhomalt` nor `nhomalt_XY` - which is why X-linked genes
            // (FMR1, IDS in the round-2 review) previously looked as though
            // gnomAD had never seen a null individual.
            let null_individuals = hc + hemizygotes;
            let freq_lower_95 = poisson_lower_95(null_individuals) / individuals as f64;
            details.insert("gnomad_individuals".into(), serde_json::json!(individuals));
            details.insert("null_individuals".into(), serde_json::json!(null_individuals));
            details.insert(
                "hom_freq_lower_95".into(),
                serde_json::json!(freq_lower_95),
            );
            details.insert(
                "prevalence_threshold".into(),
                serde_json::json!(config.bs2_hom_prevalence_threshold),
            );
            details.insert("bs2_ar_min_hom".into(), serde_json::json!(config.bs2_ar_min_hom));

            // Describe what was actually counted, so a reviewer can tell an
            // X-linked hemizygous observation from an autosomal homozygous one.
            let observed = if hemizygotes > 0 && hc > 0 {
                format!("{} homozygote(s) and {} hemizygote(s)", hc, hemizygotes)
            } else if hemizygotes > 0 {
                format!("{} hemizygote(s)", hemizygotes)
            } else {
                format!("{} homozygote(s)", hc)
            };

            if null_individuals < config.bs2_ar_min_hom {
                (
                    false,
                    true,
                    format!(
                        "{} in gnomAD is below the floor of {}; too few to establish tolerance in healthy adults",
                        observed, config.bs2_ar_min_hom
                    ),
                )
            } else if freq_lower_95 > config.bs2_hom_prevalence_threshold {
                (
                    true,
                    true,
                    format!(
                        "{} among {} gnomAD individuals (95% lower bound on their frequency {:.2e} exceeds the {:.0e} prevalence bar), so complete loss of this allele is tolerated in healthy adults",
                        observed, individuals, freq_lower_95, config.bs2_hom_prevalence_threshold
                    ),
                )
            } else {
                (
                    false,
                    true,
                    format!(
                        "{} among {} gnomAD individuals gives a 95% lower bound of {:.2e}, which does not exceed the {:.0e} prevalence bar; consistent with a late-onset, reduced-penetrance or variably expressive disorder rather than tolerance",
                        observed, individuals, freq_lower_95, config.bs2_hom_prevalence_threshold
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

/// BS3: Well-established in vitro or in vivo functional studies show no
/// damaging effect.
///
/// The benign-direction twin of PS3, read from the same curated
/// `--functional-evidence` file. See [`crate::functional`].
fn evaluate_bs3(
    input: &ClassificationInput,
    _config: &AcmgConfig,
) -> EvidenceCriterion {
    super::functional_criterion(
        input,
        crate::functional::FunctionalCriterion::Bs3,
        EvidenceDirection::Benign,
        "Requires curated functional study evidence showing no damaging effect — supply one with --functional-evidence",
    )
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
    use crate::test_support::minimal_input;
    use crate::sa_extract::GnomadData;

    fn make_input(gnomad: Option<GnomadData>) -> ClassificationInput {
        ClassificationInput {
            gnomad,
            ..minimal_input()
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
        let cfg = AcmgConfig {
            bs2_ad_min_ac: 20,
            ..Default::default()
        };
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
    fn test_bs2_prevalence_bar_default_is_the_measured_one() {
        // The bar was chosen from a sweep of the full ClinVar 2-star+
        // benchmark, not from convention, and the value is quoted in
        // docs/ACMG.md and the config doc comment. Changing it silently would
        // leave those claims wrong, so pin it here: move it deliberately or
        // not at all.
        assert_eq!(AcmgConfig::default().bs2_hom_prevalence_threshold, 1e-3);
    }

    #[test]
    fn test_bs2_common_mendelian_disorder_is_not_called_tolerant() {
        // The case the 1e-3 bar exists for. A recessive disorder with a
        // prevalence near 1 in 1,000 - hearing loss, alpha-1 antitrypsin
        // deficiency, familial Mediterranean fever in a high-prevalence
        // population - will genuinely have affected homozygotes in gnomAD.
        // At 500 homozygotes among 730 K the observation is still consistent
        // with disease rather than with tolerance, and a 1e-5 bar would have
        // called it benign evidence.
        let input = recessive_input(500, 1_460_000);
        assert!(!evaluate_bs2(&input, &AcmgConfig::default()).met);

        let permissive = AcmgConfig {
            bs2_hom_prevalence_threshold: 1e-5,
            ..Default::default()
        };
        assert!(
            evaluate_bs2(&input, &permissive).met,
            "the old bar is what made this look like tolerance"
        );
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
        // 1,000 homozygotes among 730 K individuals puts the 95 % lower bound
        // on their frequency above the 1e-3 prevalence bar, i.e. above the
        // prevalence of the most common Mendelian disorders.
        let r = evaluate_bs2(&recessive_input(1_000, 1_460_000), &AcmgConfig::default());
        assert!(r.met, "1,000 homozygotes among 730 K individuals should fire BS2");
        assert!(r.summary.contains("tolerated in healthy adults"));
    }

    #[test]
    fn test_bs2_ar_scales_with_cohort_size() {
        // The same raw count means different things at different cohort sizes.
        // 3 homozygotes among 500 individuals is a homozygote frequency far
        // above any Mendelian prevalence; 3 among 730 K is not.
        let small = evaluate_bs2(&recessive_input(3, 1_000), &AcmgConfig::default());
        let large = evaluate_bs2(&recessive_input(3, 1_460_000), &AcmgConfig::default());
        assert!(small.met, "3 homozygotes among 500 individuals should fire BS2");
        assert!(
            !large.met,
            "3 homozygotes among 730 K individuals should not fire BS2"
        );
    }

    /// A non-PAR chrX site: `hemizygotes` XY carriers among `an_xy` XY samples,
    /// with no XX homozygotes. Modelled on the real gnomAD v4.1 IDS-region
    /// records, where AN_XY runs about a third of AN.
    fn x_linked_input(hemizygotes: u64, an: u64, an_xy: u64) -> ClassificationInput {
        use crate::sa_extract::OmimData;
        make_input_omim(
            Some(GnomadData {
                all_ac: Some(hemizygotes + 10),
                all_hc: Some(0),
                all_an: Some(an),
                non_par: true,
                all_ac_xy: Some(hemizygotes),
                all_an_xy: Some(an_xy),
                ..Default::default()
            }),
            Some(OmimData {
                mim_number: None,
                phenotypes: Some(vec!["X-linked recessive disorder".into()]),
            }),
        )
    }

    #[test]
    fn test_bs2_counts_hemizygotes_on_x() {
        // gnomAD calls XY samples haploid outside the PAR, so a hemizygote is
        // recorded in AC_XY and in neither nhomalt nor nhomalt_XY. Before this,
        // an X-linked gene with thousands of hemizygous carriers looked to BS2
        // as though gnomAD had never seen a null individual, which is why FMR1
        // and IDS were flagged in the round-2 medical-genetics review.
        let r = evaluate_bs2(&x_linked_input(1275, 1_039_980, 327_798), &AcmgConfig::default());
        assert!(r.met, "1275 hemizygotes among 327 K XY samples should fire BS2");
        assert!(r.summary.contains("hemizygote"), "got: {}", r.summary);
    }

    #[test]
    fn test_bs2_single_hemizygote_does_not_fire() {
        // The cohort-scaling test applies to hemizygotes exactly as it does to
        // homozygotes: one observation in a large cohort is what a late-onset
        // or reduced-penetrance disorder looks like.
        let r = evaluate_bs2(&x_linked_input(1, 1_039_980, 327_798), &AcmgConfig::default());
        assert!(!r.met);
    }

    #[test]
    fn test_bs2_par_site_does_not_read_xy_counts_as_hemizygotes() {
        // Inside a pseudoautosomal region XY samples are diploid, so AC_XY is
        // an allele count and not a count of null individuals. Without the
        // non_par flag those carriers must not be counted.
        let mut input = x_linked_input(1275, 1_039_980, 327_798);
        input.gnomad.as_mut().unwrap().non_par = false;
        let r = evaluate_bs2(&input, &AcmgConfig::default());
        assert!(!r.met, "AC_XY must only count as hemizygotes outside the PAR");
    }

    #[test]
    fn test_bs2_autosomal_result_unchanged_by_hemizygote_support() {
        // An autosomal record carries no XY columns, so the individual count
        // must stay AN/2 and the outcome must match the pre-hemizygote
        // behaviour exactly.
        let with_hom = evaluate_bs2(&recessive_input(1_000, 1_460_000), &AcmgConfig::default());
        assert!(with_hom.met);
        assert!(with_hom.summary.contains("homozygote"));
        assert!(!with_hom.summary.contains("hemizygote"));
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
        assert!(r.summary.contains("max population AF"), "got: {}", r.summary);
    }

    // ── BS1 against the filtering allele frequency (Whiffin 2017) ──

    #[test]
    fn test_bs1_prefers_filtering_af_over_point_estimate() {
        // A 3 % point estimate in one group crosses the 1 % BS1 bar, but a
        // filtering AF of 0.2 % says that estimate rests on very few alleles.
        // Testing the point estimate is what turns thin sampling in one
        // population into a benign call.
        let input = make_input(Some(GnomadData {
            all_af: Some(0.0004),
            mid_af: Some(0.03),
            all_an: Some(2_000_000),
            faf95_max: Some(0.002),
            ..Default::default()
        }));
        let r = evaluate_bs1(&input, &AcmgConfig::default());
        assert!(!r.met, "BS1 must not fire on a point estimate the FAF contradicts");
        assert!(r.summary.contains("filtering AF"), "got: {}", r.summary);
    }

    #[test]
    fn test_bs1_fires_when_filtering_af_is_itself_high() {
        // A genuinely common allele: the CI lower bound is still above the bar.
        let input = make_input(Some(GnomadData {
            all_af: Some(0.03),
            mid_af: Some(0.03),
            all_an: Some(2_000_000),
            faf95_max: Some(0.028),
            ..Default::default()
        }));
        let r = evaluate_bs1(&input, &AcmgConfig::default());
        assert!(r.met);
    }

    #[test]
    fn test_bs1_falls_back_to_point_af_without_faf_columns() {
        // An annotation database built before the FAF columns existed must
        // behave exactly as it did before.
        let input = make_input(Some(GnomadData {
            all_af: Some(0.001),
            eas_af: Some(0.05),
            all_an: Some(2_000_000),
            ..Default::default()
        }));
        let r = evaluate_bs1(&input, &AcmgConfig::default());
        assert!(r.met);
        assert!(r.summary.contains("max population AF"));
    }

    #[test]
    fn test_filtering_af_can_be_disabled() {
        let input = make_input(Some(GnomadData {
            all_af: Some(0.0004),
            mid_af: Some(0.03),
            all_an: Some(2_000_000),
            faf95_max: Some(0.002),
            ..Default::default()
        }));
        let config = AcmgConfig {
            use_filtering_af: false,
            ..Default::default()
        };
        let r = evaluate_bs1(&input, &config);
        assert!(r.met, "with the FAF disabled BS1 reverts to the point estimate");
    }
}
