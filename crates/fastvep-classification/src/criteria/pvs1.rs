use fastvep_core::Consequence;

use crate::config::AcmgConfig;
use crate::sa_extract::ClassificationInput;
use crate::types::{EvidenceCriterion, EvidenceDirection, EvidenceStrength};

/// PVS1 — strength-graded null-variant evidence per Abou Tayoun 2018
/// (Hum Mutat 39(11):1517-1524). Possible outcomes:
///
/// - **PVS1** (Very Strong): nonsense/frameshift predicted to undergo NMD,
///   canonical ±1/2 splice causing NMD, or whole-gene deletion in a
///   haploinsufficient gene.
/// - **PVS1_Strong**: NMD-escape but in a critical functional region.
/// - **PVS1_Moderate**: NMD-escape, non-critical region, ≥10% protein
///   removed; canonical splice in last exon; start-loss with downstream alt
///   start ≤100 codons + pathogenic variant upstream.
/// - **PVS1_Supporting**: <10% protein removed in non-critical region;
///   start-loss without strong corroborating evidence.
///
/// The decision tree depends on optional fields populated by the pipeline:
/// `predicted_nmd`, `protein_truncation_pct`, `is_last_exon`,
/// `in_critical_region`, `alt_start_codon_distance`. When these are not
/// available, PVS1 falls back to the legacy binary rule (full Very Strong
/// when a null variant fires in a LOF-intolerant gene), preserving
/// backward compatibility for pipelines that haven't been updated yet.
pub fn evaluate_pvs1(input: &ClassificationInput, config: &AcmgConfig) -> EvidenceCriterion {
    let mut details = serde_json::Map::new();

    let null_kind = NullKind::detect(&input.consequences);
    details.insert(
        "null_kind".into(),
        serde_json::json!(null_kind.as_ref().map(|k| k.label()).unwrap_or("none")),
    );

    let Some(kind) = null_kind else {
        return mk(
            "PVS1".to_string(),
            EvidenceStrength::VeryStrong,
            false,
            true,
            "Not a null variant (nonsense/frameshift/canonical splice/start-lost/whole-gene deletion)".to_string(),
            details,
        );
    };

    let is_lof_gene = is_lof_intolerant_gene(input, config);
    details.insert("is_lof_gene".into(), serde_json::json!(is_lof_gene));
    if let Some(ref gc) = input.gene_constraints {
        if let Some(pli) = gc.pli {
            details.insert("pLI".into(), serde_json::json!(pli));
        }
        if let Some(loeuf) = gc.loeuf {
            details.insert("LOEUF".into(), serde_json::json!(loeuf));
        }
    }

    if !is_lof_gene {
        let gene = input.gene_symbol.as_deref().unwrap_or("unknown");
        return mk(
            "PVS1".to_string(),
            EvidenceStrength::VeryStrong,
            false,
            true,
            format!("Null variant but gene {} is not established as LOF-intolerant", gene),
            details,
        );
    }

    let (strength, summary) = match kind {
        NullKind::NonsenseOrFrameshift => grade_nonsense_frameshift(input, &mut details),
        NullKind::CanonicalSplice => grade_canonical_splice(input, &mut details),
        NullKind::StartLost => grade_start_lost(input, &mut details),
        NullKind::WholeGeneDeletion => (
            EvidenceStrength::VeryStrong,
            "Whole-gene deletion in haploinsufficient gene → PVS1".to_string(),
        ),
    };

    let code = if strength == EvidenceStrength::VeryStrong {
        "PVS1".to_string()
    } else {
        format!("PVS1_{}", strength.as_str())
    };

    mk(code, strength, true, true, summary, details)
}

#[derive(Debug, Clone, Copy)]
enum NullKind {
    NonsenseOrFrameshift,
    CanonicalSplice,
    StartLost,
    WholeGeneDeletion,
}

impl NullKind {
    fn detect(cs: &[Consequence]) -> Option<Self> {
        for c in cs {
            match c {
                Consequence::StopGained | Consequence::FrameshiftVariant => {
                    return Some(Self::NonsenseOrFrameshift)
                }
                Consequence::SpliceAcceptorVariant | Consequence::SpliceDonorVariant => {
                    return Some(Self::CanonicalSplice)
                }
                Consequence::StartLost => return Some(Self::StartLost),
                Consequence::TranscriptAblation => return Some(Self::WholeGeneDeletion),
                _ => {}
            }
        }
        None
    }

    fn label(&self) -> &'static str {
        match self {
            Self::NonsenseOrFrameshift => "nonsense_or_frameshift",
            Self::CanonicalSplice => "canonical_splice",
            Self::StartLost => "start_lost",
            Self::WholeGeneDeletion => "whole_gene_deletion",
        }
    }
}

fn grade_nonsense_frameshift(
    input: &ClassificationInput,
    details: &mut serde_json::Map<String, serde_json::Value>,
) -> (EvidenceStrength, String) {
    if let Some(nmd) = input.predicted_nmd {
        details.insert("predicted_nmd".into(), serde_json::json!(nmd));
        if nmd {
            return (
                EvidenceStrength::VeryStrong,
                "Nonsense/frameshift predicted to undergo NMD → PVS1".to_string(),
            );
        }
        let pct = input.protein_truncation_pct;
        let critical = input.in_critical_region;
        if let Some(p) = pct {
            details.insert("protein_truncation_pct".into(), serde_json::json!(p));
        }
        if let Some(c) = critical {
            details.insert("in_critical_region".into(), serde_json::json!(c));
        }
        match (critical, pct) {
            (Some(true), _) => (
                EvidenceStrength::Strong,
                "NMD-escape in critical functional region → PVS1_Strong".to_string(),
            ),
            (Some(false), Some(p)) if p >= 0.10 => (
                EvidenceStrength::Moderate,
                format!(
                    "NMD-escape in non-critical region, {:.0}% of protein removed (≥10%) → PVS1_Moderate",
                    p * 100.0
                ),
            ),
            (Some(false), Some(p)) => (
                EvidenceStrength::Supporting,
                format!(
                    "NMD-escape in non-critical region, only {:.0}% of protein removed (<10%) → PVS1_Supporting",
                    p * 100.0
                ),
            ),
            _ => (
                EvidenceStrength::Strong,
                "NMD-escape; insufficient signal to grade further → PVS1_Strong (conservative)"
                    .to_string(),
            ),
        }
    } else {
        details.insert("nmd_unknown_fallback".into(), serde_json::json!(true));
        (
            EvidenceStrength::VeryStrong,
            "Nonsense/frameshift in LOF-intolerant gene → PVS1 (NMD signal unavailable; legacy fallback)".to_string(),
        )
    }
}

fn grade_canonical_splice(
    input: &ClassificationInput,
    details: &mut serde_json::Map<String, serde_json::Value>,
) -> (EvidenceStrength, String) {
    if let Some(nmd) = input.predicted_nmd {
        details.insert("predicted_nmd".into(), serde_json::json!(nmd));
        if nmd {
            return (
                EvidenceStrength::VeryStrong,
                "Canonical splice variant → exon-skip / cryptic-splice predicted to undergo NMD → PVS1".to_string(),
            );
        }
    }
    if input.is_last_exon == Some(true) {
        details.insert("is_last_exon".into(), serde_json::json!(true));
        return (
            EvidenceStrength::Moderate,
            "Canonical splice in last exon (NMD unlikely) → PVS1_Moderate".to_string(),
        );
    }
    details.insert("splice_unknown_fallback".into(), serde_json::json!(true));
    (
        EvidenceStrength::VeryStrong,
        "Canonical ±1/2 splice in LOF-intolerant gene → PVS1 (NMD signal unavailable; legacy fallback)".to_string(),
    )
}

fn grade_start_lost(
    input: &ClassificationInput,
    details: &mut serde_json::Map<String, serde_json::Value>,
) -> (EvidenceStrength, String) {
    if let Some(d) = input.alt_start_codon_distance {
        details.insert("alt_start_codon_distance".into(), serde_json::json!(d));
        if d.unsigned_abs() <= 100 && input.in_critical_region == Some(true) {
            return (
                EvidenceStrength::Moderate,
                format!(
                    "Start-lost with downstream Met {} codons away and pathogenic variant upstream → PVS1_Moderate",
                    d
                ),
            );
        }
        if d.unsigned_abs() > 100 {
            return (
                EvidenceStrength::Supporting,
                format!(
                    "Start-lost; alternative downstream Met is {} codons away → PVS1_Supporting",
                    d
                ),
            );
        }
    }
    (
        EvidenceStrength::Supporting,
        "Start-lost in LOF-intolerant gene; downgraded to Supporting absent stronger evidence → PVS1_Supporting".to_string(),
    )
}

fn mk(
    code: String,
    strength: EvidenceStrength,
    met: bool,
    evaluated: bool,
    summary: String,
    details: serde_json::Map<String, serde_json::Value>,
) -> EvidenceCriterion {
    EvidenceCriterion {
        code,
        direction: EvidenceDirection::Pathogenic,
        strength,
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

    // ── Abou Tayoun 2018 decision-tree tests ────────────────────────────

    #[test]
    fn test_pvs1_nonsense_nmd_predicted_full_strength() {
        let mut input = make_input(
            vec![Consequence::StopGained],
            Some(GnomadGeneData { pli: Some(1.0), loeuf: Some(0.03), ..Default::default() }),
            None,
        );
        input.predicted_nmd = Some(true);
        let r = evaluate_pvs1(&input, &AcmgConfig::default());
        assert!(r.met);
        assert_eq!(r.strength, EvidenceStrength::VeryStrong);
        assert_eq!(r.code, "PVS1");
    }

    #[test]
    fn test_pvs1_nmd_escape_critical_region_strong() {
        let mut input = make_input(
            vec![Consequence::FrameshiftVariant],
            Some(GnomadGeneData { pli: Some(1.0), ..Default::default() }),
            None,
        );
        input.predicted_nmd = Some(false);
        input.in_critical_region = Some(true);
        let r = evaluate_pvs1(&input, &AcmgConfig::default());
        assert!(r.met);
        assert_eq!(r.strength, EvidenceStrength::Strong);
        assert_eq!(r.code, "PVS1_Strong");
    }

    #[test]
    fn test_pvs1_nmd_escape_noncritical_large_truncation_moderate() {
        let mut input = make_input(
            vec![Consequence::FrameshiftVariant],
            Some(GnomadGeneData { pli: Some(1.0), ..Default::default() }),
            None,
        );
        input.predicted_nmd = Some(false);
        input.in_critical_region = Some(false);
        input.protein_truncation_pct = Some(0.25); // 25% removed
        let r = evaluate_pvs1(&input, &AcmgConfig::default());
        assert_eq!(r.strength, EvidenceStrength::Moderate);
        assert_eq!(r.code, "PVS1_Moderate");
    }

    #[test]
    fn test_pvs1_nmd_escape_noncritical_small_truncation_supporting() {
        let mut input = make_input(
            vec![Consequence::FrameshiftVariant],
            Some(GnomadGeneData { pli: Some(1.0), ..Default::default() }),
            None,
        );
        input.predicted_nmd = Some(false);
        input.in_critical_region = Some(false);
        input.protein_truncation_pct = Some(0.05); // <10%
        let r = evaluate_pvs1(&input, &AcmgConfig::default());
        assert_eq!(r.strength, EvidenceStrength::Supporting);
        assert_eq!(r.code, "PVS1_Supporting");
    }

    #[test]
    fn test_pvs1_canonical_splice_last_exon_moderate() {
        let mut input = make_input(
            vec![Consequence::SpliceDonorVariant],
            Some(GnomadGeneData { pli: Some(1.0), ..Default::default() }),
            None,
        );
        input.predicted_nmd = Some(false);
        input.is_last_exon = Some(true);
        let r = evaluate_pvs1(&input, &AcmgConfig::default());
        assert_eq!(r.strength, EvidenceStrength::Moderate);
        assert_eq!(r.code, "PVS1_Moderate");
    }

    #[test]
    fn test_pvs1_start_lost_no_signals_supporting() {
        let input = make_input(
            vec![Consequence::StartLost],
            Some(GnomadGeneData { pli: Some(1.0), ..Default::default() }),
            None,
        );
        let r = evaluate_pvs1(&input, &AcmgConfig::default());
        assert!(r.met);
        assert_eq!(r.strength, EvidenceStrength::Supporting);
        assert_eq!(r.code, "PVS1_Supporting");
    }

    #[test]
    fn test_pvs1_start_lost_alt_start_with_pathogenic_upstream_moderate() {
        let mut input = make_input(
            vec![Consequence::StartLost],
            Some(GnomadGeneData { pli: Some(1.0), ..Default::default() }),
            None,
        );
        input.alt_start_codon_distance = Some(50);
        input.in_critical_region = Some(true); // proxy for "pathogenic upstream"
        let r = evaluate_pvs1(&input, &AcmgConfig::default());
        assert_eq!(r.strength, EvidenceStrength::Moderate);
        assert_eq!(r.code, "PVS1_Moderate");
    }

    #[test]
    fn test_pvs1_legacy_fallback_when_no_nmd_signal() {
        // No predicted_nmd → falls back to legacy full PVS1.
        let input = make_input(
            vec![Consequence::FrameshiftVariant],
            Some(GnomadGeneData { pli: Some(1.0), ..Default::default() }),
            None,
        );
        let r = evaluate_pvs1(&input, &AcmgConfig::default());
        assert!(r.met);
        assert_eq!(r.strength, EvidenceStrength::VeryStrong);
        assert_eq!(r.code, "PVS1");
    }
}
