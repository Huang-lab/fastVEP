use super::frequency_gate::benign_blocker;
use crate::config::AcmgConfig;
use crate::sa_extract::ClassificationInput;
use crate::types::{EvidenceCriterion, EvidenceDirection, EvidenceStrength};

/// Evaluate all benign strong criteria: BS1, BS2, BS3, BS4.
pub fn evaluate_all(input: &ClassificationInput, config: &AcmgConfig) -> Vec<EvidenceCriterion> {
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
        0.0, 0.0513, 0.3554, 0.8177, 1.3663, 1.9702, 2.6130, 3.2853, 3.9808, 4.6952, 5.4254,
        6.1690, 6.9242, 7.6896, 8.4639, 9.2463, 10.0360, 10.8324, 11.6350, 12.4432, 13.2547,
    ];
    if (observed as usize) < TABLE.len() {
        return TABLE[observed as usize];
    }
    let k = observed as f64;
    // Wilson-Hilferty: lower ≈ k * (1 - 1/(9k) - 1.645/(3*sqrt(k)))^3
    let term = 1.0 - 1.0 / (9.0 * k) - 1.645 / (3.0 * k.sqrt());
    (k * term * term * term).max(0.0)
}

/// 95 % lower confidence bound on the frequency of `count` events among
/// `individuals`: the quantity both arms of BS2 hold against a prevalence.
fn lower_bound_frequency(count: u64, individuals: u64) -> f64 {
    poisson_lower_95(count) / individuals as f64
}

/// BS1: Allele frequency is greater than expected for disorder.
fn evaluate_bs1(input: &ClassificationInput, config: &AcmgConfig) -> EvidenceCriterion {
    let threshold = config.effective_bs1_threshold(input.gene_symbol.as_deref());

    let mut details = serde_json::Map::new();
    details.insert("af_threshold".into(), serde_json::json!(threshold));

    // Checked ahead of the generic frequency gate, which would otherwise answer
    // first and give a vaguer reason. "On the curated exception list" is an
    // evaluated conclusion with a citation, not an inability to evaluate.
    if let Some(exc) = super::frequency_gate::curated_exception(input, config, "BS1") {
        let summary = super::frequency_gate::record_exception("BS1", exc, &mut details);
        return EvidenceCriterion {
            code: "BS1".to_string(),
            direction: EvidenceDirection::Benign,
            strength: EvidenceStrength::Strong,
            default_strength: EvidenceStrength::Strong,
            met: false,
            evaluated: true,
            summary,
            details: serde_json::Value::Object(details),
        };
    }

    if let Some(blocker) = benign_blocker(input, config) {
        blocker.record(&mut details);
        super::frequency_gate::note_withheld_benign_frequency(
            input,
            config,
            threshold,
            &mut details,
        );
        return EvidenceCriterion {
            code: "BS1".to_string(),
            direction: EvidenceDirection::Benign,
            strength: EvidenceStrength::Strong,
            default_strength: EvidenceStrength::Strong,
            met: false,
            evaluated: false,
            summary: format!("BS1 not evaluated: {}", blocker.reason),
            details: serde_json::Value::Object(details),
        };
    }

    let (met, summary) = if let Some(ref gnomad) = input.gnomad {
        // Coverage first, the same requirement as BA1: see `insufficient_an`.
        if let Some(summary) =
            super::frequency_gate::insufficient_an("BS1", gnomad.all_an, config, &mut details)
        {
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

        // BS1 does not fire where BA1 already covers the frequency. The bar
        // read here is the gene's own, not the global default: a gene whose
        // VCEP publishes a looser BA1 - ABCA4's is 0.163 - would otherwise
        // have BS1 stand down at 0.06 while BA1 declines too, leaving a common
        // variant with no frequency evidence in either direction.
        //
        // A curated BA1 exception does not reopen BS1. Ghosh 2018's nine are
        // removed from BA1 because they are frequent *and* not benign, and
        // collecting Strong benign evidence on HFE p.Cys282Tyr instead of
        // Standalone benign evidence would miss the point of the list. An
        // entry that means to leave BS1 available says so through `blocks`.
        let ba1_bar = config.effective_ba1_threshold(input.gene_symbol.as_deref());
        if test_af > ba1_bar {
            (
                false,
                format!(
                    "BA1 takes precedence ({}={:.4} > BA1 threshold {:.2e})",
                    af_label, test_af, ba1_bar
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
/// Both arms ask the same question - are there more healthy individuals with the
/// allele than the disorder could account for? - of the individuals the disorder
/// would affect, so both compare a 95 % Poisson lower bound on their frequency
/// against a maximum credible prevalence. That makes the criterion scale with
/// the size of the cohort behind the observation: one homozygote among gnomAD's
/// 730 K individuals is what a late-onset or reduced-penetrance disorder looks
/// like, not evidence of tolerance, and nor is a handful of heterozygous
/// carriers of a dominant one.
///
/// - **Recessive / X-linked** (or unknown inheritance): count the individuals
///   with no functional copy - homozygotes plus, on a non-PAR sex-chromosome
///   site, hemizygotes - and require both an absolute floor of
///   `bs2_ar_min_hom` and a lower bound on their frequency above
///   `bs2_hom_prevalence_threshold`.
/// - **Dominant**: count the individuals carrying the allele at all - `AC - nhomalt`,
///   so a homozygote counts once - and require both an absolute floor of
///   `bs2_ad_min_ac` alleles and a lower bound on the carrier frequency above
///   `bs2_ad_prevalence_threshold`. The floor is the old rule and is kept:
///   Richards 2015 says "observed in unaffected adult", and a singleton or
///   doubleton is sequencing-noise plausibility, not evidence the variant is
///   tolerated; ClinGen VCEPs commonly use AC ≥ 5 (Hereditary Cancer VCEP,
///   Lynch Syndrome curation guide). The prevalence bar is what the floor alone
///   could not do, which is tell nine carriers among 1.4 M alleles from nine
///   hundred.
///
/// A gene is dominant here when its disease annotation says so in words: an
/// OMIM-style inheritance ("autosomal dominant") or "dominant" in the disease
/// name. ClinGen's own `/AD` tag, which the `.oga` carries on most of its
/// entries, is not read, so a gene known to ClinGen only by that tag takes the
/// recessive arm. Reading it for BS2 was measured and not done: the criterion
/// asks for full penetrance at an early age, which is a property of the
/// disorder and not of its inheritance, and six ClinVar Pathogenic or Likely
/// pathogenic variants moved toward benign (limitation 20 in `docs/ACMG.md`
/// names them). The predicate is also shared with PM2, PM3 and BP2.
///
/// BS2 is not evaluated for a record with no AN, since there is then no cohort
/// size to test a count against. It does not take BA1's and BS1's minimum AN:
/// we read that as SVI guidance about the noise of a frequency estimate, and
/// the lower bound used here already accounts for a small cohort.
///
/// Inheritance is inferred from the disease-gene `.oga` (ClinGen GDV
/// preferred, OMIM accepted as legacy).
fn evaluate_bs2(input: &ClassificationInput, config: &AcmgConfig) -> EvidenceCriterion {
    let mut details = serde_json::Map::new();

    // A hypomorphic allele is expected to be seen homozygous in unaffected
    // people - that is the disease model, not evidence against it.
    if let Some(exc) = super::frequency_gate::curated_exception(input, config, "BS2") {
        let summary = super::frequency_gate::record_exception("BS2", exc, &mut details);
        return EvidenceCriterion {
            code: "BS2".to_string(),
            direction: EvidenceDirection::Benign,
            strength: EvidenceStrength::Strong,
            default_strength: EvidenceStrength::Strong,
            met: false,
            evaluated: true,
            summary,
            details: serde_json::Value::Object(details),
        };
    }

    if let Some(blocker) = benign_blocker(input, config) {
        blocker.record(&mut details);
        return EvidenceCriterion {
            code: "BS2".to_string(),
            direction: EvidenceDirection::Benign,
            strength: EvidenceStrength::Strong,
            default_strength: EvidenceStrength::Strong,
            met: false,
            evaluated: false,
            summary: format!("BS2 not evaluated: {}", blocker.reason),
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
        // The cohort size below is read from AN, and a record with none has no
        // cohort to test a count against: `.max(1)` used to turn that into one
        // individual and make any bar vacuous. A small AN is not that problem.
        // The lower bound already accounts for how few individuals there were,
        // and we read BA1's and BS1's minimum AN as SVI guidance about the noise
        // of a frequency *estimate*, which a count test does not make. Applying
        // it here, with the carrier bar off, moved 693 ClinVar 2-star+ calls,
        // 690 of them away from the ClinVar classification and 3 toward it.
        if matches!(gnomad.all_an, None | Some(0)) {
            details.insert("an_missing".into(), serde_json::json!(true));
            return EvidenceCriterion {
                code: "BS2".to_string(),
                direction: EvidenceDirection::Benign,
                strength: EvidenceStrength::Strong,
                default_strength: EvidenceStrength::Strong,
                met: false,
                evaluated: false,
                summary: "BS2 not evaluated: gnomAD AN unavailable, so there is no cohort size to test the count against".to_string(),
                details: serde_json::Value::Object(details),
            };
        }
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
            details.insert(
                "gnomad_xy_individuals".into(),
                serde_json::json!(hemizygote_an),
            );
        }

        // Individuals surveyed. XX samples contribute two alleles each and XY
        // samples on a non-PAR site contribute one, so the cohort size is
        // `(AN - AN_XY) / 2 + AN_XY`. With no XY data this is exactly `AN / 2`,
        // the autosomal count used before hemizygotes were extracted.
        let individuals = (an.saturating_sub(hemizygote_an) / 2 + hemizygote_an).max(1);
        details.insert("gnomad_individuals".into(), serde_json::json!(individuals));

        if is_dominant && !is_recessive {
            // A dominant-only gene is decided here and nowhere else: a variant
            // that clears the floor but not the bar must not fall into the
            // recessive arm below and be explained, or worse fired, by its
            // homozygotes.
            let carriers = ac.saturating_sub(hc);
            let carrier_freq_lower_95 = lower_bound_frequency(carriers, individuals);
            let bar = config.bs2_ad_prevalence_threshold;
            details.insert("gnomad_carriers".into(), serde_json::json!(carriers));
            details.insert(
                "carrier_freq_lower_95".into(),
                serde_json::json!(carrier_freq_lower_95),
            );
            details.insert("ad_prevalence_threshold".into(), serde_json::json!(bar));
            details.insert(
                "bs2_ad_min_ac".into(),
                serde_json::json!(config.bs2_ad_min_ac),
            );

            let homozygous = if hc > 0 {
                format!(", {hc} of them homozygous")
            } else {
                String::new()
            };
            if ac < config.bs2_ad_min_ac {
                (
                    false,
                    true,
                    format!(
                        "AC={} below BS2 threshold ({} required for a dominant disorder)",
                        ac, config.bs2_ad_min_ac
                    ),
                )
            } else if bar <= 0.0 {
                // The bar is switched off, so the floor is the whole test.
                (
                    true,
                    true,
                    format!(
                        "Observed in gnomAD as ≥{} alleles (AC={}) for a dominant disorder; no prevalence bar is configured",
                        config.bs2_ad_min_ac, ac
                    ),
                )
            } else if carrier_freq_lower_95 > bar {
                (
                    true,
                    true,
                    format!(
                        "{} of {} gnomAD individuals carry this allele{} (AC={}); the 95% lower bound on the carrier frequency, {:.2e}, exceeds the {:e} prevalence bar for a dominant disorder, so more healthy adults carry it than the disorder could account for",
                        carriers, individuals, homozygous, ac, carrier_freq_lower_95, bar
                    ),
                )
            } else {
                (
                    false,
                    true,
                    format!(
                        "{} of {} gnomAD individuals carry this allele{} (AC={}), a 95% lower bound on the carrier frequency of {:.2e} that does not exceed the {:e} prevalence bar for a dominant disorder; this few carriers is what a rare or reduced-penetrance dominant disorder looks like, not evidence of tolerance",
                        carriers, individuals, homozygous, ac, carrier_freq_lower_95, bar
                    ),
                )
            }
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
            let freq_lower_95 = lower_bound_frequency(null_individuals, individuals);
            details.insert(
                "null_individuals".into(),
                serde_json::json!(null_individuals),
            );
            details.insert("hom_freq_lower_95".into(), serde_json::json!(freq_lower_95));
            details.insert(
                "prevalence_threshold".into(),
                serde_json::json!(config.bs2_hom_prevalence_threshold),
            );
            details.insert(
                "bs2_ar_min_hom".into(),
                serde_json::json!(config.bs2_ar_min_hom),
            );

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
                        "{} among {} gnomAD individuals (95% lower bound on their frequency {:.2e} exceeds the {:e} prevalence bar), so complete loss of this allele is tolerated in healthy adults",
                        observed, individuals, freq_lower_95, config.bs2_hom_prevalence_threshold
                    ),
                )
            } else {
                (
                    false,
                    true,
                    format!(
                        "{} among {} gnomAD individuals gives a 95% lower bound of {:.2e}, which does not exceed the {:e} prevalence bar; consistent with a late-onset, reduced-penetrance or variably expressive disorder rather than tolerance",
                        observed, individuals, freq_lower_95, config.bs2_hom_prevalence_threshold
                    ),
                )
            }
        } else {
            (false, true, "No homozygotes observed in gnomAD".to_string())
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
fn evaluate_bs3(input: &ClassificationInput, _config: &AcmgConfig) -> EvidenceCriterion {
    super::functional_criterion(
        input,
        crate::functional::FunctionalCriterion::Bs3,
        EvidenceDirection::Benign,
        "Requires curated functional study evidence showing no damaging effect — supply one with --functional-evidence",
    )
}

/// BS4: Lack of segregation in affected members of a family.
fn evaluate_bs4(_input: &ClassificationInput, _config: &AcmgConfig) -> EvidenceCriterion {
    EvidenceCriterion {
        code: "BS4".to_string(),
        direction: EvidenceDirection::Benign,
        strength: EvidenceStrength::Strong,
        default_strength: EvidenceStrength::Strong,
        met: false,
        evaluated: false,
        summary:
            "Requires multi-generation pedigree with affection status to assess lack of segregation"
                .to_string(),
        details: serde_json::Value::Null,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sa_extract::GnomadData;
    use crate::test_support::minimal_input;

    fn make_input(gnomad: Option<GnomadData>) -> ClassificationInput {
        ClassificationInput {
            gnomad,
            ..minimal_input()
        }
    }

    /// The three shipped hypomorphic alleles, in the form the pipeline
    /// actually produces: a transcript-prefixed HGVS and a coordinate. Each is
    /// frequent enough and homozygous often enough that BS1 and BS2 would both
    /// fire without the exception, and each is pathogenic in trans with a null
    /// allele.
    #[test]
    fn test_hypomorphic_alleles_block_bs1_and_bs2_in_pipeline_form() {
        let cases = [
            (
                "GAA",
                "ENST00000302262.8:c.-32-13T>G",
                "17",
                80_104_542u64,
                "T",
                "G",
                5.4e-3,
                23u64,
            ),
            (
                "CFTR",
                "ENST00000003084.11:c.1210-11T>G",
                "7",
                117_548_630,
                "T",
                "G",
                9.8e-3,
                27,
            ),
            (
                "SPTA1",
                "ENST00000643759.2:c.4339-99C>T",
                "1",
                158_643_524,
                "G",
                "A",
                6.6e-3,
                40,
            ),
        ];
        let config = AcmgConfig::default();
        for (gene, hgvs, chrom, pos, r, a, af, hc) in cases {
            let input = ClassificationInput {
                gene_symbol: Some(gene.to_string()),
                hgvs_c: Some(hgvs.to_string()),
                variant_coordinates: Some((chrom.to_string(), pos, r.to_string(), a.to_string())),
                gnomad: Some(GnomadData {
                    all_af: Some(af),
                    all_an: Some(1_400_000),
                    all_hc: Some(hc),
                    faf95_max: Some(af),
                    ..Default::default()
                }),
                ..minimal_input()
            };
            let bs1 = evaluate_bs1(&input, &config);
            assert!(!bs1.met, "{} BS1 must not fire: {}", gene, bs1.summary);
            assert!(bs1.evaluated, "{} BS1 is a conclusion, not a gap", gene);
            assert!(bs1.summary.contains("exception list"), "{}", bs1.summary);

            let bs2 = evaluate_bs2(&input, &config);
            assert!(!bs2.met, "{} BS2 must not fire: {}", gene, bs2.summary);
            assert!(bs2.summary.contains("exception list"), "{}", bs2.summary);
        }
    }

    /// The three alleles the round-6 review named, in the form the pipeline
    /// produces them. A wrong coordinate, allele or chromosome spelling, or a
    /// wrong HGVS token, would leave one route to the entry inert, and the
    /// benchmark that would show it is not part of CI. Either route alone is
    /// enough for a match, so with both supplied a typo in one is hidden by the
    /// other: each is exercised on its own.
    #[test]
    fn test_reviewer_named_alleles_block_bs1_and_bs2_in_pipeline_form() {
        // (gene, HGVS as the pipeline writes it, chrom, pos, ref, alt, gnomAD
        // AF, population-maximum AF, homozygotes, AN)
        let cases = [
            (
                "G6PD",
                "ENST00000393562.10:c.563C>T",
                "X",
                154_534_419u64,
                "G",
                "A",
                1.5e-3,
                4.4e-2,
                21u64,
                1_097_898u64,
            ),
            (
                "G6PD",
                "ENST00000393562.10:c.202G>A",
                "X",
                154_536_002,
                "C",
                "T",
                3.7e-3,
                0.131,
                173,
                1_098_048,
            ),
            (
                "SERPINC1",
                "ENST00000367698.4:c.439A>G",
                "1",
                173_911_984,
                "T",
                "C",
                1.8e-4,
                7.0e-3,
                3,
                1_461_832,
            ),
        ];
        let config = AcmgConfig::default();
        for (gene, hgvs, chrom, pos, r, a, af, pop_af, hc, an) in cases {
            for (route, hgvs_c, coordinates) in [
                (
                    "by coordinate",
                    None,
                    Some((chrom.to_string(), pos, r.to_string(), a.to_string())),
                ),
                ("by HGVS", Some(hgvs.to_string()), None),
            ] {
                let input = ClassificationInput {
                    gene_symbol: Some(gene.to_string()),
                    hgvs_c,
                    variant_coordinates: coordinates,
                    gnomad: Some(GnomadData {
                        all_af: Some(af),
                        all_an: Some(an),
                        all_hc: Some(hc),
                        faf95_max: Some(pop_af),
                        ..Default::default()
                    }),
                    ..minimal_input()
                };
                let bs1 = evaluate_bs1(&input, &config);
                assert!(
                    !bs1.met,
                    "{gene} {hgvs} {route}: BS1 must not fire: {}",
                    bs1.summary
                );
                assert!(bs1.summary.contains("exception list"), "{}", bs1.summary);
                let bs2 = evaluate_bs2(&input, &config);
                assert!(
                    !bs2.met,
                    "{gene} {hgvs} {route}: BS2 must not fire: {}",
                    bs2.summary
                );
                assert!(bs2.summary.contains("exception list"), "{}", bs2.summary);
                assert_eq!(
                    bs2.details["exception_gene"],
                    serde_json::json!(gene),
                    "{gene} {hgvs} {route}"
                );
            }
        }
    }

    #[test]
    fn test_a_ba1_only_exception_does_not_reach_bs1() {
        // Ghosh's entries are scoped to BA1 by `blocks`, so BS1 evaluates
        // normally on them. Below the BA1 bar, where BA1's precedence rule is
        // not in the way, that means BS1 fires.
        let input = ClassificationInput {
            gene_symbol: Some("GJB2".to_string()),
            hgvs_c: Some("ENST00000382844.2:c.109G>A".to_string()),
            variant_coordinates: Some((
                "13".to_string(),
                20_189_473,
                "C".to_string(),
                "T".to_string(),
            )),
            gnomad: Some(GnomadData {
                all_af: Some(0.02),
                all_an: Some(1_400_000),
                faf95_max: Some(0.02),
                ..Default::default()
            }),
            ..minimal_input()
        };
        let bs1 = evaluate_bs1(&input, &AcmgConfig::default());
        assert!(bs1.met, "{}", bs1.summary);
    }

    #[test]
    fn test_bs1_defers_to_the_genes_own_ba1_bar() {
        // ABCA4's published BA1 is 0.163, three times looser than the global
        // 0.05. At 0.06 the old global comparison stood BS1 down while BA1
        // also declined, so a 6 % allele collected no frequency evidence at
        // all. Against the gene's own bar, BS1 fires.
        let mut config = AcmgConfig::default();
        config.gene_overrides.insert(
            "ABCA4".to_string(),
            crate::config::GeneOverride {
                mechanism: None,
                ba1_af_threshold: Some(0.163),
                bs1_af_threshold: Some(0.0163),
                pm2_af_threshold: None,
                disabled_criteria: vec![],
                strength_overrides: std::collections::HashMap::new(),
                disorders: std::collections::HashMap::new(),
            },
        );
        let input = ClassificationInput {
            gene_symbol: Some("ABCA4".to_string()),
            gnomad: Some(GnomadData {
                all_af: Some(0.06),
                all_an: Some(1_400_000),
                faf95_max: Some(0.06),
                ..Default::default()
            }),
            ..minimal_input()
        };
        let bs1 = evaluate_bs1(&input, &config);
        assert!(bs1.met, "{}", bs1.summary);
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
        // Five homozygotes among the thousand individuals a 2,000-allele site
        // surveys: a lower bound of 2e-3, above the 1e-3 bar. The same five in
        // gnomAD's 730 K individuals would not be (see the AR tests below).
        let input = make_input(Some(GnomadData {
            all_hc: Some(5),
            all_an: Some(2_000),
            ..Default::default()
        }));
        let result = evaluate_bs2(&input, &AcmgConfig::default());
        assert!(result.met, "{}", result.summary);
    }

    #[test]
    fn test_bs2_no_homozygotes() {
        let input = make_input(Some(GnomadData {
            all_hc: Some(0),
            all_an: Some(1_460_000),
            ..Default::default()
        }));
        let result = evaluate_bs2(&input, &AcmgConfig::default());
        assert!(!result.met);
        assert!(result.evaluated);
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

    /// A dominant gene's record: `ac` alleles, `hc` of them in homozygotes, among
    /// `an` called chromosomes.
    fn dominant_input(ac: u64, hc: u64, an: Option<u64>) -> ClassificationInput {
        use crate::sa_extract::OmimData;
        make_input_omim(
            Some(GnomadData {
                all_ac: Some(ac),
                all_hc: Some(hc),
                all_an: an,
                ..Default::default()
            }),
            Some(OmimData {
                mim_number: None,
                phenotypes: Some(vec!["dominant disorder".into()]),
            }),
        )
    }

    #[test]
    fn test_bs2_ad_gene_many_carriers_fires() {
        // 500 carriers among the 500 K individuals of a million-allele site: a
        // 95% lower bound of about 9e-4, above the 3e-4 bar, so more healthy
        // adults carry this allele than a dominant disorder could account for.
        let r = evaluate_bs2(
            &dominant_input(500, 0, Some(1_000_000)),
            &AcmgConfig::default(),
        );
        assert!(r.met, "{}", r.summary);
        assert!(r.summary.contains("carry this allele"), "{}", r.summary);
        assert_eq!(r.details["gnomad_carriers"], serde_json::json!(500));
        assert_eq!(r.details["gnomad_individuals"], serde_json::json!(500_000));
    }

    /// The case the floor alone got wrong. SNRNP200 p.Ser1087Leu has nine
    /// carriers among 1.4 M alleles - ClinVar Pathogenic - and AC ≥ 5 called it
    /// tolerated.
    #[test]
    fn test_bs2_ad_gene_few_carriers_in_a_large_cohort_does_not_fire() {
        let r = evaluate_bs2(
            &dominant_input(9, 0, Some(1_461_868)),
            &AcmgConfig::default(),
        );
        assert!(!r.met, "{}", r.summary);
        assert!(r.evaluated);
        assert!(r.summary.contains("does not exceed"), "{}", r.summary);
        assert!(r.summary.contains("3e-4"), "{}", r.summary);
    }

    /// `AC - nhomalt` counts individuals: a homozygote contributes two alleles
    /// and is one person. 480 alleles of which 240 are homozygotes is 240
    /// carriers, below the bar, where 480 heterozygotes would clear it.
    #[test]
    fn test_bs2_ad_gene_counts_a_homozygote_once() {
        let an = Some(1_460_000);
        let heterozygous = evaluate_bs2(&dominant_input(480, 0, an), &AcmgConfig::default());
        assert!(heterozygous.met, "{}", heterozygous.summary);
        let homozygous = evaluate_bs2(&dominant_input(480, 240, an), &AcmgConfig::default());
        assert!(!homozygous.met, "{}", homozygous.summary);
        assert_eq!(
            homozygous.details["gnomad_carriers"],
            serde_json::json!(240)
        );
    }

    /// A dominant-only gene is decided in its own arm. A variant that clears the
    /// floor but not the bar, and has homozygotes, must not fall into the
    /// recessive arm and be explained - or fired - by them.
    #[test]
    fn test_bs2_ad_gene_that_fails_the_bar_does_not_fall_into_the_recessive_arm() {
        let r = evaluate_bs2(
            &dominant_input(20, 3, Some(1_461_868)),
            &AcmgConfig::default(),
        );
        assert!(!r.met);
        assert!(r.summary.contains("dominant disorder"), "{}", r.summary);
        assert!(!r.summary.contains("late-onset"), "{}", r.summary);
        assert!(r.details.get("hom_freq_lower_95").is_none());
        assert!(r.details.get("null_individuals").is_none());

        // And with the recessive bar lowered to where the homozygotes alone
        // would have fired it, the dominant verdict still stands.
        let permissive_recessive = AcmgConfig {
            bs2_hom_prevalence_threshold: 1e-9,
            bs2_ar_min_hom: 1,
            ..Default::default()
        };
        assert!(
            !evaluate_bs2(
                &dominant_input(20, 3, Some(1_461_868)),
                &permissive_recessive
            )
            .met
        );
    }

    #[test]
    fn test_bs2_ad_gene_singleton_floor_and_bar_apply_together() {
        // 7 carriers in a 20,000-allele site clears the bar (lower bound 3.3e-4)
        // and the default floor of 5, and so fires; raising the floor to 20 is
        // enough to stop it. Neither condition replaces the other.
        let input = dominant_input(7, 0, Some(20_000));
        assert!(evaluate_bs2(&input, &AcmgConfig::default()).met);
        let strict = AcmgConfig {
            bs2_ad_min_ac: 20,
            ..Default::default()
        };
        let r = evaluate_bs2(&input, &strict);
        assert!(!r.met, "AC=7 < raised floor 20 should not fire");
        assert!(r.summary.contains("below BS2 threshold"), "{}", r.summary);
    }

    #[test]
    fn test_bs2_ad_prevalence_bar_default_is_the_measured_one() {
        // Chosen from the sweep in docs/ACMG.md and quoted in the config doc
        // comment: move it deliberately or not at all.
        assert_eq!(AcmgConfig::default().bs2_ad_prevalence_threshold, 3e-4);
    }

    #[test]
    fn test_bs2_ad_bar_can_be_switched_off_and_is_reported_as_configured() {
        let off = AcmgConfig {
            bs2_ad_prevalence_threshold: 0.0,
            ..Default::default()
        };
        // With the bar off the floor is the whole test, as it was.
        let r = evaluate_bs2(&dominant_input(9, 0, Some(1_461_868)), &off);
        assert!(r.met, "{}", r.summary);
        assert!(r.summary.contains("no prevalence bar"), "{}", r.summary);

        // A bar with more than one significant digit is printed whole.
        let odd = AcmgConfig {
            bs2_ad_prevalence_threshold: 1.5e-3,
            ..Default::default()
        };
        let r = evaluate_bs2(&dominant_input(9, 0, Some(1_461_868)), &odd);
        assert!(r.summary.contains("1.5e-3"), "{}", r.summary);
    }

    /// On a non-PAR sex-chromosome site an XY sample contributes one allele to
    /// AN, so the individuals surveyed are `(AN - AN_XY) / 2 + AN_XY` and not
    /// `AN / 2`. 200 carriers among 300,000 XX and 400,000 XY individuals have a
    /// lower bound of 2.5e-4, under the 3e-4 bar; read as 500,000 individuals
    /// they would have 3.5e-4 and clear it.
    #[test]
    fn test_bs2_ad_gene_counts_xy_individuals_once_on_a_non_par_site() {
        let mut input = dominant_input(200, 0, Some(1_000_000));
        {
            let gnomad = input.gnomad.as_mut().unwrap();
            gnomad.non_par = true;
            gnomad.all_an_xy = Some(400_000);
            gnomad.all_ac_xy = Some(0);
        }
        let r = evaluate_bs2(&input, &AcmgConfig::default());
        assert_eq!(r.details["gnomad_individuals"], serde_json::json!(700_000));
        assert!(!r.met, "{}", r.summary);

        // The same counts on an autosome survey 500,000 individuals.
        input.gnomad.as_mut().unwrap().non_par = false;
        let r = evaluate_bs2(&input, &AcmgConfig::default());
        assert_eq!(r.details["gnomad_individuals"], serde_json::json!(500_000));
        assert!(r.met, "{}", r.summary);
    }

    /// Every BS2 verdict rests on this bound: a table up to 20 events and the
    /// Wilson-Hilferty approximation above it, standing in for
    /// `0.5 * chi2_{0.05, 2k}`. The values are that quantile, from scipy.
    #[test]
    fn test_poisson_lower_95_agrees_with_the_exact_bound() {
        const EXACT: [f64; 21] = [
            0.0, 0.051293, 0.355362, 0.817691, 1.366319, 1.970150, 2.613013, 3.285316, 3.980823,
            4.695228, 5.425405, 6.169007, 6.924246, 7.689641, 8.463899, 9.246259, 10.035990,
            10.832100, 11.634300, 12.442000, 13.254700,
        ];
        for (k, exact) in EXACT.iter().enumerate() {
            let got = poisson_lower_95(k as u64);
            assert!((got - exact).abs() < 2e-3, "k={k}: {got} against {exact}");
        }
        // Above the table the approximation stays within 0.02% of the exact bound.
        for (k, exact) in [
            (21u64, 14.0720),
            (50, 38.9647),
            (100, 84.1393),
            (300, 272.0900),
            (719, 675.4704),
            (7409, 7267.9891),
        ] {
            let got = poisson_lower_95(k);
            assert!(
                (got - exact).abs() / exact < 2e-4,
                "k={k}: {got} against {exact}"
            );
        }
    }

    /// Both arms read the cohort size from AN, and without one a single
    /// individual's worth of cohort made any bar vacuous.
    #[test]
    fn test_bs2_is_not_evaluated_without_a_cohort_size() {
        for (label, input) in [
            ("dominant, AN missing", dominant_input(9, 0, None)),
            ("dominant, AN zero", dominant_input(9, 0, Some(0))),
            ("recessive, AN missing", {
                let mut i = recessive_input(5, 1_460_000);
                i.gnomad.as_mut().unwrap().all_an = None;
                i
            }),
        ] {
            let r = evaluate_bs2(&input, &AcmgConfig::default());
            assert!(!r.met, "{label}");
            assert!(!r.evaluated, "{label}: {}", r.summary);
            assert!(
                r.summary.starts_with("BS2 not evaluated"),
                "{label}: {}",
                r.summary
            );
            assert_eq!(r.details["an_missing"], serde_json::json!(true), "{label}");
        }
    }

    /// A small AN is evaluated, not withheld: BA1's and BS1's minimum AN guards
    /// a frequency estimate against noise, and a count test over a small cohort
    /// is already as strict as the cohort is small. Nine carriers among 250
    /// individuals is a carrier frequency far above the bar, and a few
    /// homozygotes among as many are above the recessive one.
    #[test]
    fn test_bs2_reads_a_small_cohort_as_the_small_cohort_it_is() {
        let dominant = evaluate_bs2(&dominant_input(9, 0, Some(500)), &AcmgConfig::default());
        assert!(dominant.met, "{}", dominant.summary);
        assert_eq!(
            dominant.details["gnomad_individuals"],
            serde_json::json!(250)
        );
        let recessive = evaluate_bs2(&recessive_input(5, 500), &AcmgConfig::default());
        assert!(recessive.met, "{}", recessive.summary);
        // And the same counts in a large cohort are not evidence of anything.
        assert!(
            !evaluate_bs2(
                &dominant_input(9, 0, Some(1_461_868)),
                &AcmgConfig::default()
            )
            .met
        );
        assert!(!evaluate_bs2(&recessive_input(5, 1_461_868), &AcmgConfig::default()).met);
    }

    /// A gene annotated both ways takes the recessive arm, as it always has.
    #[test]
    fn test_bs2_gene_that_is_dominant_and_recessive_takes_the_recessive_arm() {
        use crate::sa_extract::OmimData;
        let mut input = recessive_input(1_000, 1_460_000);
        input.omim = Some(OmimData {
            mim_number: None,
            phenotypes: Some(vec![
                "dominant disorder".into(),
                "recessive disorder".into(),
            ]),
        });
        let r = evaluate_bs2(&input, &AcmgConfig::default());
        assert!(r.met, "{}", r.summary);
        assert!(r.summary.contains("homozygote"), "{}", r.summary);
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
        assert!(
            !r.met,
            "1 homozygote among 730 K individuals must not fire BS2"
        );
    }

    #[test]
    fn test_bs2_ar_many_homozygotes_in_large_cohort_fires() {
        // 1,000 homozygotes among 730 K individuals puts the 95 % lower bound
        // on their frequency above the 1e-3 prevalence bar, i.e. above the
        // prevalence of the most common Mendelian disorders.
        let r = evaluate_bs2(&recessive_input(1_000, 1_460_000), &AcmgConfig::default());
        assert!(
            r.met,
            "1,000 homozygotes among 730 K individuals should fire BS2"
        );
        assert!(r.summary.contains("tolerated in healthy adults"));
    }

    #[test]
    fn test_bs2_ar_scales_with_cohort_size() {
        // The same raw count means different things at different cohort sizes.
        // 4 homozygotes among 1,000 individuals - the smallest cohort BS2 will
        // read, 2,000 alleles - is a homozygote frequency above any Mendelian
        // prevalence; 4 among 730 K is not.
        let small = evaluate_bs2(&recessive_input(4, 2_000), &AcmgConfig::default());
        let large = evaluate_bs2(&recessive_input(4, 1_460_000), &AcmgConfig::default());
        assert!(
            small.met,
            "4 homozygotes among 1,000 individuals should fire BS2: {}",
            small.summary
        );
        assert!(
            !large.met,
            "4 homozygotes among 730 K individuals should not fire BS2"
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
        let r = evaluate_bs2(
            &x_linked_input(1275, 1_039_980, 327_798),
            &AcmgConfig::default(),
        );
        assert!(
            r.met,
            "1275 hemizygotes among 327 K XY samples should fire BS2"
        );
        assert!(r.summary.contains("hemizygote"), "got: {}", r.summary);
    }

    #[test]
    fn test_bs2_single_hemizygote_does_not_fire() {
        // The cohort-scaling test applies to hemizygotes exactly as it does to
        // homozygotes: one observation in a large cohort is what a late-onset
        // or reduced-penetrance disorder looks like.
        let r = evaluate_bs2(
            &x_linked_input(1, 1_039_980, 327_798),
            &AcmgConfig::default(),
        );
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
        assert!(
            !r.met,
            "AC_XY must only count as hemizygotes outside the PAR"
        );
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
        assert!(
            !r.evaluated,
            "BS2 should be NotEvaluated, not a negative call"
        );
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
        assert!(
            r.summary.contains("max population AF"),
            "got: {}",
            r.summary
        );
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
        assert!(
            !r.met,
            "BS1 must not fire on a point estimate the FAF contradicts"
        );
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
        assert!(
            r.met,
            "with the FAF disabled BS1 reverts to the point estimate"
        );
    }
}
