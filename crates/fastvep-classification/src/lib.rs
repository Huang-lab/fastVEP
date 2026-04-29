//! ACMG-AMP variant classification engine for fastVEP.
//!
//! Implements the 28 evidence criteria from Richards et al. 2015
//! (Standards and guidelines for the interpretation of sequence variants)
//! and produces a 5-tier classification: Pathogenic, Likely Pathogenic,
//! Uncertain Significance (VUS), Likely Benign, Benign.
//!
//! Incorporates ClinGen SVI (Sequence Variant Interpretation) working group
//! recommendations including calibrated REVEL thresholds for PP3/BP4 and
//! PM2 downgrade to Supporting strength.

pub mod combiner;
pub mod config;
pub mod criteria;
pub mod sa_extract;
pub mod types;

pub use config::{AcmgConfig, TrioConfig};
pub use sa_extract::{
    extract_classification_input, ClassificationInput, CompanionVariant, GenotypeInfo,
};
pub use types::{AcmgClassification, AcmgResult, EvidenceCounts, EvidenceCriterion};

/// Classify a variant using the ACMG-AMP framework.
///
/// Takes a `ClassificationInput` (extracted from pipeline annotation data)
/// and an `AcmgConfig` (with thresholds and gene-specific overrides).
///
/// Returns an `AcmgResult` containing the classification, all evaluated
/// criteria, triggered rule, and evidence counts.
pub fn classify(input: &ClassificationInput, config: &AcmgConfig) -> AcmgResult {
    // Evaluate all 28 criteria
    let criteria = criteria::evaluate_all_criteria(input, config);

    // Count met criteria by direction/strength
    let counts = EvidenceCounts::from_criteria(&criteria);

    // Apply combination rules
    let (classification, triggered_rule) = combiner::combine(&criteria);

    AcmgResult {
        shorthand: classification.shorthand().to_string(),
        classification,
        criteria,
        triggered_rule,
        counts,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sa_extract::*;
    use fastvep_core::{Consequence, Impact};

    /// Helper to build a ClassificationInput with common defaults.
    fn make_input(
        consequences: Vec<Consequence>,
        impact: Impact,
        gene_symbol: &str,
    ) -> ClassificationInput {
        ClassificationInput {
            consequences,
            impact,
            gene_symbol: Some(gene_symbol.to_string()),
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
    fn test_classify_common_variant_benign() {
        let mut input = make_input(
            vec![Consequence::MissenseVariant],
            Impact::Moderate,
            "TEST",
        );
        input.gnomad = Some(GnomadData {
            all_af: Some(0.10),
            afr_af: Some(0.15),
            ..Default::default()
        });
        let result = classify(&input, &AcmgConfig::default());
        assert_eq!(result.classification, AcmgClassification::Benign);
        assert_eq!(result.shorthand, "B");
        assert!(result.triggered_rule.as_deref() == Some("BA1"));
    }

    #[test]
    fn test_classify_frameshift_lof_gene_likely_pathogenic() {
        let mut input = make_input(
            vec![Consequence::FrameshiftVariant],
            Impact::High,
            "BRCA1",
        );
        // Absent from gnomAD → PM2_Supporting
        // High pLI → PVS1
        input.gene_constraints = Some(GnomadGeneData {
            pli: Some(1.0),
            loeuf: Some(0.03),
            mis_z: Some(2.5),
            syn_z: Some(0.5),
        });
        let config = AcmgConfig::default();
        let result = classify(&input, &config);

        // PVS1 (VeryStrong) + PM2_Supporting (Supporting) = PVS=1, PP=1
        // Per ClinGen SVI (Sept 2020): PVS + >=1 PP → Likely Pathogenic
        // Bayesian Post_P = 0.988, within LP range (0.90–0.99)
        assert!(result.counts.pathogenic_very_strong >= 1); // PVS1
        assert!(result.counts.pathogenic_supporting >= 1); // PM2_Supporting
        assert_eq!(result.classification, AcmgClassification::LikelyPathogenic);
    }

    #[test]
    fn test_classify_frameshift_revel_does_not_drive_pathogenic() {
        // Pre-PR1, a high REVEL score on a frameshift incorrectly fired
        // PP3_Strong, pushing PVS1 + REVEL → Pathogenic via the PVS+PS rule.
        // Per Pejaver 2022, REVEL is calibrated for missense only and must
        // not contribute to non-missense classification. Without other
        // pathogenic-Strong evidence, PVS1 + PM2_Supporting tops out at LP
        // via the ClinGen SVI PVS+PP rule.
        let mut input = make_input(
            vec![Consequence::FrameshiftVariant],
            Impact::High,
            "BRCA1",
        );
        input.gene_constraints = Some(GnomadGeneData {
            pli: Some(1.0),
            loeuf: Some(0.03),
            mis_z: Some(2.5),
            syn_z: Some(0.5),
        });
        input.revel = Some(RevelData { score: Some(0.95) });

        let config = AcmgConfig::default();
        let result = classify(&input, &config);
        assert_eq!(result.classification, AcmgClassification::LikelyPathogenic);
    }

    #[test]
    fn test_classify_synonymous_no_splice_not_conserved() {
        let mut input = make_input(
            vec![Consequence::SynonymousVariant],
            Impact::Low,
            "TEST",
        );
        input.splice_ai = Some(SpliceAiData {
            ds_ag: Some(0.01),
            ds_al: Some(0.02),
            ds_dg: Some(0.01),
            ds_dl: Some(0.01),
            ..Default::default()
        });
        input.phylop = Some(0.3);
        input.revel = Some(RevelData { score: Some(0.10) });
        // Provide gnomAD data to prevent PM2_Supporting from firing (which would cause conflict)
        input.gnomad = Some(GnomadData {
            all_af: Some(0.005), // Above PM2 threshold (0.0001) but below BS1 (0.01)
            ..Default::default()
        });

        let result = classify(&input, &AcmgConfig::default());
        // BP7 (synonymous + no splice + not conserved): benign supporting
        // BP4_Moderate (REVEL=0.10 <= 0.183): counts as benign strong
        // BS + BP → Likely Benign
        assert_eq!(result.classification, AcmgClassification::LikelyBenign);
    }

    #[test]
    fn test_classify_vus_no_data() {
        let input = make_input(
            vec![Consequence::MissenseVariant],
            Impact::Moderate,
            "UNKNOWN_GENE",
        );
        let result = classify(&input, &AcmgConfig::default());
        // No SA data at all → most criteria not evaluable
        // PM2_Supporting fires (absent from gnomAD)
        // But just 1 supporting isn't enough for any LP rule
        assert_eq!(
            result.classification,
            AcmgClassification::UncertainSignificance
        );
    }

    #[test]
    fn test_classify_conflicting_evidence() {
        // When both pathogenic AND benign criteria are met but benign
        // doesn't reach standalone/>=2 Strong thresholds, it's VUS
        let mut input = make_input(
            vec![Consequence::FrameshiftVariant],
            Impact::High,
            "GENE",
        );
        input.gene_constraints = Some(GnomadGeneData {
            pli: Some(1.0),
            loeuf: Some(0.03),
            ..Default::default()
        });
        // PVS1 should fire (frameshift in LOF-intolerant gene)
        // BS1 fires (AF > 0.01) but BS2 does NOT (no homozygotes)
        // So we have 1 BS (not 2), plus some BP → conflicting with PVS1
        input.gnomad = Some(GnomadData {
            all_af: Some(0.02),
            all_hc: Some(0), // No homozygotes → BS2 does not fire
            ..Default::default()
        });

        let result = classify(&input, &AcmgConfig::default());
        // PVS1 (pathogenic) + BS1 (benign) = conflicting → VUS
        assert_eq!(
            result.classification,
            AcmgClassification::UncertainSignificance
        );
        assert!(result
            .triggered_rule
            .as_deref()
            .unwrap_or("")
            .contains("Conflicting"));
    }

    #[test]
    fn test_acmg_result_serialization() {
        let input = make_input(
            vec![Consequence::MissenseVariant],
            Impact::Moderate,
            "TEST",
        );
        let result = classify(&input, &AcmgConfig::default());
        let json = serde_json::to_value(&result).unwrap();

        assert!(json.get("classification").is_some());
        assert!(json.get("shorthand").is_some());
        assert!(json.get("criteria").is_some());
        assert!(json.get("counts").is_some());

        let criteria = json.get("criteria").unwrap().as_array().unwrap();
        assert!(!criteria.is_empty());

        // Each criterion should have required fields
        for c in criteria {
            assert!(c.get("code").is_some());
            assert!(c.get("direction").is_some());
            assert!(c.get("strength").is_some());
            assert!(c.get("met").is_some());
            assert!(c.get("evaluated").is_some());
            assert!(c.get("summary").is_some());
        }
    }
}
