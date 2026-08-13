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
    // Per-gene where a VCEP has published one, global otherwise.
    let threshold = config.effective_ba1_threshold(input.gene_symbol.as_deref());

    let mut details = serde_json::Map::new();
    details.insert("af_threshold".into(), serde_json::json!(threshold));

    // Detect whether this allele is on the BA1 exception list. Coordinate
    // first, HGVS second - see `Ba1Exception::matches` for why the order is
    // load-bearing rather than a preference.
    //
    // Checked ahead of the generic frequency gate below, which would otherwise
    // answer first for these variants and give a vaguer reason. Both block
    // BA1, but "on the ClinGen exception list" is the citable one, and it is
    // an evaluated conclusion rather than an inability to evaluate. HFE
    // c.845G>A reaches both: ClinVar labels it low-penetrance *and* Ghosh 2018
    // lists it.
    let exception_match: Option<&Ba1Exception> =
        super::frequency_gate::curated_exception(input, config, "BA1");

    if let Some(exc) = exception_match {
        // Same wording for all three frequency criteria. It says "curated"
        // rather than "ClinGen" because the list is no longer only Ghosh
        // 2018's nine: the hypomorphic entries are curated from the disease
        // mechanism and have no published table behind them.
        let summary = super::frequency_gate::record_exception("BA1", exc, &mut details);
        return EvidenceCriterion {
            code: "BA1".to_string(),
            direction: EvidenceDirection::Benign,
            strength: EvidenceStrength::Standalone,
            default_strength: EvidenceStrength::Standalone,
            met: false,
            evaluated: true,
            summary,
            details: serde_json::Value::Object(details),
        };
    }

    // Frequencies from a homology-confounded gene (or for a variant ClinVar
    // labels low-penetrance) cannot support a standalone benign call either.
    if let Some(blocker) = super::frequency_gate::benign_blocker(input, config) {
        blocker.record(&mut details);
        super::frequency_gate::note_withheld_benign_frequency(
            input, config, threshold, &mut details,
        );
        return EvidenceCriterion {
            code: "BA1".to_string(),
            direction: EvidenceDirection::Benign,
            strength: EvidenceStrength::Standalone,
            default_strength: EvidenceStrength::Standalone,
            met: false,
            evaluated: false,
            summary: format!("BA1 not evaluated: {}", blocker.reason),
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

            if af > threshold {
                (
                    true,
                    format!(
                        "Common variant: {}={:.3e} exceeds the {:.2e} threshold",
                        af_label, af, threshold
                    ),
                )
            } else {
                (
                    false,
                    format!(
                        "{}={:.3e} does not exceed the {:.2e} threshold",
                        af_label, af, threshold
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

    /// A gnomAD record at `af`, in `gene`.
    fn at_frequency(gene: &str, af: f64) -> ClassificationInput {
        ClassificationInput {
            impact: fastvep_core::Impact::Modifier,
            gene_symbol: Some(gene.to_string()),
            gnomad: Some(GnomadData {
                all_af: Some(af),
                all_an: Some(100_000),
                ..Default::default()
            }),
            ..minimal_input()
        }
    }

    fn with_gene_ba1(gene: &str, threshold: f64) -> AcmgConfig {
        let mut cfg = AcmgConfig::default();
        cfg.gene_overrides.insert(
            gene.to_string(),
            crate::config::GeneOverride {
                mechanism: None,
                ba1_af_threshold: Some(threshold),
                bs1_af_threshold: None,
                pm2_af_threshold: None,
                disabled_criteria: vec![],
                strength_overrides: std::collections::HashMap::new(),
                disorders: std::collections::HashMap::new(),
            },
        );
        cfg
    }

    #[test]
    fn test_ba1_uses_a_published_per_gene_bar_in_both_directions() {
        // CDKL5's published bar is 8.3e-5, six hundred times tighter than the
        // 0.05 default: a variant at 1e-4 is common for that disorder and not
        // remotely common in general.
        let cfg = with_gene_ba1("CDKL5", 8.3e-5);
        assert!(evaluate_ba1(&at_frequency("CDKL5", 1e-4), &cfg).met);
        assert!(!evaluate_ba1(&at_frequency("CDKL5", 5e-5), &cfg).met);
        // A gene the table says nothing about keeps the global bar.
        assert!(!evaluate_ba1(&at_frequency("TP53", 1e-4), &cfg).met);

        // ABCA4's runs the other way: 0.163, three times looser than default,
        // so 10 % is not standalone-benign evidence there.
        let cfg = with_gene_ba1("ABCA4", 0.163);
        assert!(!evaluate_ba1(&at_frequency("ABCA4", 0.10), &cfg).met);
        assert!(evaluate_ba1(&at_frequency("ABCA4", 0.20), &cfg).met);
    }

    /// Every variant on the Ghosh 2018 list, in the form the pipeline actually
    /// emits: `ENST00000357618.10:c.845G>A`, not a bare `c.845G>A`.
    ///
    /// This is the test that was missing. The exception matcher compared the
    /// whole HGVS string, the unit tests fed it bare `c.` tokens the pipeline
    /// never produces, and the list shipped inert - BA1 fired on ACAD9, HFE
    /// p.His63Asp, MEFV p.Pro369Ser and PIBF1 regardless, calling all four
    /// Benign. Verified end to end afterwards, not just here.
    #[test]
    fn test_every_ghosh_2018_exception_blocks_ba1_in_pipeline_form() {
        // (gene, transcript-prefixed HGVS as fastVEP emits it, chrom, pos, ref, alt)
        let cases = [
            ("ACAD9", "ENST00000308982.12:c.-45_-44insTAAG", "3", 128_879_647, "C", "CTAAG"),
            ("GJB2", "ENST00000382848.5:c.109G>A", "13", 20_189_473, "C", "T"),
            ("HFE", "ENST00000357618.10:c.187C>G", "6", 26_090_951, "C", "G"),
            ("HFE", "ENST00000357618.10:c.845G>A", "6", 26_092_913, "G", "A"),
            ("MEFV", "ENST00000219596.6:c.1105C>T", "16", 3_249_586, "G", "A"),
            ("MEFV", "ENST00000219596.6:c.1223G>A", "16", 3_249_468, "C", "T"),
            ("PIBF1", "ENST00000326291.11:c.1214G>A", "13", 72_835_359, "G", "A"),
            ("ACADS", "ENST00000242592.9:c.511C>T", "12", 120_737_875, "C", "T"),
            // BTD is the one no HGVS match can reach: Ghosh lists c.1330G>C on
            // NM_000060.4, fastVEP reports c.1270G>C on the Ensembl canonical.
            ("BTD", "ENST00000643237.3:c.1270G>C", "3", 15_645_186, "G", "C"),
        ];
        let cfg = AcmgConfig::default();
        for (gene, hgvs, chrom, pos, ref_allele, alt) in cases {
            let input = ClassificationInput {
                impact: fastvep_core::Impact::Modifier,
                gene_symbol: Some(gene.to_string()),
                hgvs_c: Some(hgvs.to_string()),
                variant_coordinates: Some((
                    chrom.to_string(),
                    pos,
                    ref_allele.to_string(),
                    alt.to_string(),
                )),
                gnomad: Some(GnomadData {
                    all_af: Some(0.20),
                    all_an: Some(100_000),
                    faf95_max: Some(0.20),
                    ..Default::default()
                }),
                ..minimal_input()
            };
            let result = evaluate_ba1(&input, &cfg);
            assert!(!result.met, "{gene} {hgvs}: BA1 fired despite Ghosh 2018");
            assert!(
                result.summary.contains("exception list"),
                "{gene} {hgvs}: blocked, but not by the exception list: {}",
                result.summary
            );
        }
    }

    #[test]
    fn test_ba1_exception_matches_on_coordinate_when_hgvs_cannot() {
        // BTD again, with no HGVS at all. The coordinate is what identifies
        // the variant; a c. token identifies it only relative to a transcript.
        let input = ClassificationInput {
            impact: fastvep_core::Impact::Modifier,
            gene_symbol: Some("BTD".to_string()),
            hgvs_c: None,
            variant_coordinates: Some(("chr3".to_string(), 15_645_186, "G".to_string(), "C".to_string())),
            gnomad: Some(GnomadData {
                all_af: Some(0.20),
                all_an: Some(100_000),
                faf95_max: Some(0.20),
                ..Default::default()
            }),
            ..minimal_input()
        };
        let result = evaluate_ba1(&input, &AcmgConfig::default());
        assert!(!result.met, "{}", result.summary);
        assert!(result.summary.contains("exception list"), "{}", result.summary);
    }

    #[test]
    fn test_ba1_exception_does_not_match_a_neighbouring_variant() {
        // One base away, same gene: not the listed variant, and BA1 must fire.
        let input = ClassificationInput {
            impact: fastvep_core::Impact::Modifier,
            gene_symbol: Some("BTD".to_string()),
            hgvs_c: Some("ENST00000643237.3:c.1271G>C".to_string()),
            variant_coordinates: Some(("3".to_string(), 15_645_187, "G".to_string(), "C".to_string())),
            gnomad: Some(GnomadData {
                all_af: Some(0.20),
                all_an: Some(100_000),
                faf95_max: Some(0.20),
                ..Default::default()
            }),
            ..minimal_input()
        };
        let result = evaluate_ba1(&input, &AcmgConfig::default());
        assert!(result.met, "{}", result.summary);
    }

    #[test]
    fn test_ba1_reports_tight_thresholds_readably() {
        // The summary used to format both numbers with `{:.2}`, which prints a
        // bar of 8.3e-5 as "0.00" - fine while every bar was 0.05, actively
        // misleading once per-gene bars arrived.
        let cfg = with_gene_ba1("CDKL5", 8.3e-5);
        let summary = evaluate_ba1(&at_frequency("CDKL5", 1e-4), &cfg).summary;
        assert!(summary.contains("8.30e-5"), "{summary}");
        assert!(!summary.contains("0.00 "), "{summary}");
    }
}
