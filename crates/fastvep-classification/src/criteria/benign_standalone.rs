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
        // estimates in poorly-covered populations.
        let an_ok = gnomad
            .all_an
            .map_or(true, |an| an >= config.min_an_for_frequency_criteria);
        if !an_ok {
            details.insert(
                "an_below_minimum".into(),
                serde_json::json!(gnomad.all_an),
            );
            details.insert(
                "min_an_for_frequency_criteria".into(),
                serde_json::json!(config.min_an_for_frequency_criteria),
            );
            return EvidenceCriterion {
                code: "BA1".to_string(),
                direction: EvidenceDirection::Benign,
                strength: EvidenceStrength::Standalone,
                default_strength: EvidenceStrength::Standalone,
                met: false,
                evaluated: false,
                summary: format!(
                    "BA1 not evaluated: gnomAD AN={} below minimum {} (gnomAD v4 guidance)",
                    gnomad.all_an.unwrap_or(0),
                    config.min_an_for_frequency_criteria
                ),
                details: serde_json::Value::Object(details),
            };
        }
        let max_af = gnomad.max_pop_af();
        if let Some(af) = max_af {
            details.insert("max_pop_af".into(), serde_json::json!(af));

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
                        "Common variant: max population AF={:.4} exceeds {:.2} threshold",
                        af, config.ba1_af_threshold
                    ),
                )
            } else {
                (
                    false,
                    format!(
                        "Max population AF={:.6} does not exceed {:.2} threshold",
                        af, config.ba1_af_threshold
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
    use crate::sa_extract::GnomadData;

    #[test]
    fn test_ba1_common_variant() {
        let input = ClassificationInput {
            consequences: vec![],
            impact: fastvep_core::Impact::Modifier,
            gene_symbol: None,
            is_canonical: false,
            amino_acids: None,
            protein_position: None,
            gnomad: Some(GnomadData {
                all_af: Some(0.10),
                afr_af: Some(0.15),
                ..Default::default()
            }),
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
        };
        let result = evaluate_ba1(&input, &AcmgConfig::default());
        assert!(result.met);
        assert_eq!(result.strength, EvidenceStrength::Standalone);
    }

    #[test]
    fn test_ba1_rare_variant() {
        let input = ClassificationInput {
            consequences: vec![],
            impact: fastvep_core::Impact::Modifier,
            gene_symbol: None,
            is_canonical: false,
            amino_acids: None,
            protein_position: None,
            gnomad: Some(GnomadData {
                all_af: Some(0.001),
                ..Default::default()
            }),
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
        };
        let result = evaluate_ba1(&input, &AcmgConfig::default());
        assert!(!result.met);
    }

    #[test]
    fn test_ba1_exception_list_blocks_known_pathogenic_high_af() {
        // HFE c.845G>A (p.Cys282Tyr) — hereditary hemochromatosis. ~10% AF in
        // European populations but pathogenic. Per Ghosh 2018, BA1 must NOT fire.
        let input = ClassificationInput {
            consequences: vec![],
            impact: fastvep_core::Impact::Modifier,
            gene_symbol: Some("HFE".to_string()),
            is_canonical: true,
            amino_acids: None,
            protein_position: None,
            gnomad: Some(GnomadData {
                all_af: Some(0.06),
                nfe_af: Some(0.10),
                ..Default::default()
            }),
            clinvar: None,
            revel: None,
            splice_ai: None,
            dbnsfp: None,
            phylop: None,
            gerp: None,
            gene_constraints: None,
            omim: None,
            clinvar_protein: None,
            hgvs_c: Some("c.845G>A".to_string()),
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
        let result = evaluate_ba1(&input, &AcmgConfig::default());
        assert!(!result.met, "BA1 must not fire for HFE c.845G>A (Ghosh 2018 exception)");
        assert!(result.evaluated);
        assert!(result.summary.contains("exception"));
    }

    #[test]
    fn test_ba1_exception_match_is_case_insensitive() {
        // Pipeline may emit "C.845G>a" or other casing — match must still work.
        let input = ClassificationInput {
            consequences: vec![],
            impact: fastvep_core::Impact::Modifier,
            gene_symbol: Some("hfe".to_string()),
            is_canonical: true,
            amino_acids: None,
            protein_position: None,
            gnomad: Some(GnomadData { all_af: Some(0.10), ..Default::default() }),
            clinvar: None,
            revel: None,
            splice_ai: None,
            dbnsfp: None,
            phylop: None,
            gerp: None,
            gene_constraints: None,
            omim: None,
            clinvar_protein: None,
            hgvs_c: Some("C.845G>A".to_string()),
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
        let result = evaluate_ba1(&input, &AcmgConfig::default());
        assert!(!result.met);
    }

    #[test]
    fn test_ba1_high_af_non_exception_still_fires() {
        // Same gene (HFE) but a different c. notation NOT on the exception
        // list → BA1 still fires at high AF.
        let input = ClassificationInput {
            consequences: vec![],
            impact: fastvep_core::Impact::Modifier,
            gene_symbol: Some("HFE".to_string()),
            is_canonical: true,
            amino_acids: None,
            protein_position: None,
            gnomad: Some(GnomadData { all_af: Some(0.10), ..Default::default() }),
            clinvar: None,
            revel: None,
            splice_ai: None,
            dbnsfp: None,
            phylop: None,
            gerp: None,
            gene_constraints: None,
            omim: None,
            clinvar_protein: None,
            hgvs_c: Some("c.999A>T".to_string()),
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
        let result = evaluate_ba1(&input, &AcmgConfig::default());
        assert!(result.met);
    }

    #[test]
    fn test_ba1_low_an_not_evaluated() {
        // PR10 (gnomAD v4 guidance): AN below 2000 → NotEvaluated, not Benign,
        // even at high AF — the frequency estimate is unreliable.
        let input = ClassificationInput {
            consequences: vec![],
            impact: fastvep_core::Impact::Modifier,
            gene_symbol: None,
            is_canonical: false,
            amino_acids: None,
            protein_position: None,
            gnomad: Some(GnomadData {
                all_af: Some(0.10),
                all_an: Some(500), // way below 2000
                ..Default::default()
            }),
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
        };
        let result = evaluate_ba1(&input, &AcmgConfig::default());
        assert!(!result.met);
        assert!(!result.evaluated);
        assert!(result.summary.contains("below minimum"));
    }

    #[test]
    fn test_ba1_one_pop_above_threshold() {
        let input = ClassificationInput {
            consequences: vec![],
            impact: fastvep_core::Impact::Modifier,
            gene_symbol: None,
            is_canonical: false,
            amino_acids: None,
            protein_position: None,
            gnomad: Some(GnomadData {
                all_af: Some(0.02),
                eas_af: Some(0.06), // Only EAS above 5%
                ..Default::default()
            }),
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
        };
        let result = evaluate_ba1(&input, &AcmgConfig::default());
        assert!(result.met);
    }
}
