use crate::config::AcmgConfig;
use crate::sa_extract::ClassificationInput;
use crate::types::{EvidenceCriterion, EvidenceDirection, EvidenceStrength};

/// BA1: Allele frequency is >5% in any general continental population dataset.
pub fn evaluate_ba1(
    input: &ClassificationInput,
    config: &AcmgConfig,
) -> EvidenceCriterion {
    let mut details = serde_json::Map::new();
    details.insert(
        "af_threshold".into(),
        serde_json::json!(config.ba1_af_threshold),
    );

    let (met, summary) = if let Some(ref gnomad) = input.gnomad {
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
        };
        let result = evaluate_ba1(&input, &AcmgConfig::default());
        assert!(!result.met);
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
        };
        let result = evaluate_ba1(&input, &AcmgConfig::default());
        assert!(result.met);
    }
}
