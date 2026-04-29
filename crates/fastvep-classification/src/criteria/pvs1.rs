use fastvep_core::Consequence;

use crate::config::AcmgConfig;
use crate::sa_extract::ClassificationInput;
use crate::types::{EvidenceCriterion, EvidenceDirection, EvidenceStrength};

/// PVS1: Null variant (nonsense, frameshift, canonical ±1/2 splice, initiation codon,
/// single/multi-exon deletion) in a gene where loss-of-function is a known mechanism of disease.
pub fn evaluate_pvs1(input: &ClassificationInput, config: &AcmgConfig) -> EvidenceCriterion {
    let is_null_variant = input.consequences.iter().any(|c| {
        matches!(
            c,
            Consequence::StopGained
                | Consequence::FrameshiftVariant
                | Consequence::SpliceAcceptorVariant
                | Consequence::SpliceDonorVariant
                | Consequence::StartLost
                | Consequence::TranscriptAblation
        )
    });

    let is_lof_gene = is_lof_intolerant_gene(input, config);

    let met = is_null_variant && is_lof_gene;
    let evaluated = true;

    let mut details = serde_json::Map::new();
    details.insert(
        "is_null_variant".into(),
        serde_json::Value::Bool(is_null_variant),
    );
    details.insert(
        "is_lof_gene".into(),
        serde_json::Value::Bool(is_lof_gene),
    );
    if let Some(ref gc) = input.gene_constraints {
        if let Some(pli) = gc.pli {
            details.insert("pLI".into(), serde_json::json!(pli));
        }
        if let Some(loeuf) = gc.loeuf {
            details.insert("LOEUF".into(), serde_json::json!(loeuf));
        }
    }
    if is_null_variant {
        let null_types: Vec<&str> = input
            .consequences
            .iter()
            .filter(|c| {
                matches!(
                    c,
                    Consequence::StopGained
                        | Consequence::FrameshiftVariant
                        | Consequence::SpliceAcceptorVariant
                        | Consequence::SpliceDonorVariant
                        | Consequence::StartLost
                        | Consequence::TranscriptAblation
                )
            })
            .map(|c| c.so_term())
            .collect();
        details.insert("null_consequence_types".into(), serde_json::json!(null_types));
    }

    let gene_name = input.gene_symbol.as_deref().unwrap_or("unknown");
    let summary = if met {
        format!(
            "Null variant in LOF-intolerant gene {}{}",
            gene_name,
            format_gene_constraint_summary(input, config)
        )
    } else if is_null_variant && !is_lof_gene {
        format!(
            "Null variant but gene {} is not established as LOF-intolerant",
            gene_name
        )
    } else {
        "Not a null variant (nonsense/frameshift/splice/start-lost)".to_string()
    };

    EvidenceCriterion {
        code: "PVS1".to_string(),
        direction: EvidenceDirection::Pathogenic,
        strength: EvidenceStrength::VeryStrong,
        default_strength: EvidenceStrength::VeryStrong,
        met,
        evaluated,
        summary,
        details: serde_json::Value::Object(details),
    }
}

/// Determine if a gene is LOF-intolerant using available constraint data.
fn is_lof_intolerant_gene(input: &ClassificationInput, config: &AcmgConfig) -> bool {
    // Check gene constraint scores
    if let Some(ref gc) = input.gene_constraints {
        if gc.pli.map_or(false, |p| p >= config.pli_lof_intolerant) {
            return true;
        }
        if gc.loeuf.map_or(false, |l| l <= config.loeuf_lof_intolerant) {
            return true;
        }
    }

    // Check gene-specific override for LOF mechanism
    if let Some(gene) = input.gene_symbol.as_deref() {
        if let Some(override_cfg) = config.gene_override(gene) {
            if let Some(ref mechanism) = override_cfg.mechanism {
                if mechanism.contains("LOF") {
                    return true;
                }
            }
        }
    }

    // Check OMIM for disease associations (proxy for disease gene)
    if let Some(ref omim) = input.omim {
        if omim
            .phenotypes
            .as_ref()
            .map_or(false, |p| !p.is_empty())
        {
            return true;
        }
    }

    false
}

fn format_gene_constraint_summary(input: &ClassificationInput, _config: &AcmgConfig) -> String {
    if let Some(ref gc) = input.gene_constraints {
        let mut parts = Vec::new();
        if let Some(pli) = gc.pli {
            parts.push(format!("pLI={:.2}", pli));
        }
        if let Some(loeuf) = gc.loeuf {
            parts.push(format!("LOEUF={:.2}", loeuf));
        }
        if !parts.is_empty() {
            return format!(" ({})", parts.join(", "));
        }
    }
    String::new()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sa_extract::{GnomadGeneData, OmimData};
    use fastvep_core::Impact;

    fn make_input(consequences: Vec<Consequence>, gene_constraints: Option<GnomadGeneData>, omim: Option<OmimData>) -> ClassificationInput {
        ClassificationInput {
            consequences,
            impact: Impact::High,
            gene_symbol: Some("BRCA1".to_string()),
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
            gene_constraints,
            omim,
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
    fn test_pvs1_frameshift_lof_gene() {
        let input = make_input(
            vec![Consequence::FrameshiftVariant],
            Some(GnomadGeneData { pli: Some(1.0), loeuf: Some(0.03), ..Default::default() }),
            None,
        );
        let result = evaluate_pvs1(&input, &AcmgConfig::default());
        assert!(result.met);
        assert!(result.evaluated);
        assert_eq!(result.strength, EvidenceStrength::VeryStrong);
    }

    #[test]
    fn test_pvs1_missense_not_null() {
        let input = make_input(
            vec![Consequence::MissenseVariant],
            Some(GnomadGeneData { pli: Some(1.0), loeuf: Some(0.03), ..Default::default() }),
            None,
        );
        let result = evaluate_pvs1(&input, &AcmgConfig::default());
        assert!(!result.met);
    }

    #[test]
    fn test_pvs1_null_variant_no_constraint_data() {
        let input = make_input(vec![Consequence::StopGained], None, None);
        let result = evaluate_pvs1(&input, &AcmgConfig::default());
        assert!(!result.met); // No gene constraint data = not LOF-intolerant
    }

    #[test]
    fn test_pvs1_null_variant_omim_disease_gene() {
        let input = make_input(
            vec![Consequence::StopGained],
            None,
            Some(OmimData { mim_number: Some(113705), phenotypes: Some(vec!["Breast cancer".to_string()]) }),
        );
        let result = evaluate_pvs1(&input, &AcmgConfig::default());
        assert!(result.met); // OMIM disease association is proxy for LOF gene
    }
}
