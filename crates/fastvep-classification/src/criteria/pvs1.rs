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
/// The decision tree runs on optional fields the pipeline derives from the
/// transcript: `predicted_nmd` / `nmd_escape_50nt`, `protein_truncation_pct`,
/// `is_last_exon`, `in_critical_region`, `alt_start_codon_distance`. When a
/// field is absent the tree falls back to the legacy binary rule (full Very
/// Strong for a null variant in a LOF-intolerant gene), so a pipeline that
/// does not populate them behaves exactly as before.
///
/// Two of the tree's branches are still unreachable: `alt_start_codon_distance`
/// (start-loss) and `same_splice_position_pathogenic` (PS1's splice track) have
/// no plumbing yet, so start-loss stays at PVS1_Supporting.
///
/// Which NMD prediction feeds the tree is a configuration choice - see
/// [`AcmgConfig::pvs1_nmd_50nt_rule`], which documents the measured trade
/// between guideline faithfulness and ClinVar concordance.
pub fn evaluate_pvs1(input: &ClassificationInput, config: &AcmgConfig) -> EvidenceCriterion {
    let mut details = serde_json::Map::new();

    let null_kind = NullKind::detect(&input.consequences);
    details.insert(
        "null_kind".into(),
        serde_json::json!(null_kind.as_ref().map(|k| k.label()).unwrap_or("none")),
    );

    let Some(mut kind) = null_kind else {
        return mk(
            "PVS1".to_string(),
            EvidenceStrength::VeryStrong,
            false,
            true,
            "Not a null variant (nonsense/frameshift/canonical splice/start-lost/whole-gene deletion)".to_string(),
            details,
        );
    };

    // Two gene-level preconditions, both ahead of the constraint test, because
    // constraint cannot answer either of them. The first asks whether the gene
    // causes disease at all (B7); the second whether it causes disease by
    // losing function (B6). A gene can be highly constrained and fail either.
    if let Some(reason) = super::gene_disease::validity_blocker(input, config) {
        details.insert("gene_disease_validity".into(), serde_json::json!(false));
        return mk(
            "PVS1".to_string(),
            EvidenceStrength::VeryStrong,
            false,
            false,
            reason,
            details,
        );
    }
    if let Some(reason) = super::gene_disease::lof_mechanism_blocker(input, config) {
        details.insert(
            "mechanism".into(),
            serde_json::json!(config.effective_mechanism(input.gene_symbol.as_deref())),
        );
        return mk(
            "PVS1".to_string(),
            EvidenceStrength::VeryStrong,
            false,
            true,
            reason,
            details,
        );
    }

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

    // Canonical-splice offset gate (Abou Tayoun 2018 / Walker 2023): PVS1's
    // splice track applies only to the canonical ±1/±2 dinucleotide. When the
    // pipeline provides an intronic offset that falls outside ±2, the variant
    // is not a canonical null splice variant — it was labeled splice_donor /
    // splice_acceptor by an indel spanning the region or by a deep-intronic
    // call. In that case PVS1's splice track does not apply. But if the same
    // indel ALSO deletes coding sequence (frameshift / nonsense) or hits the
    // start codon, PVS1 still applies via that null track — re-grade as such
    // rather than discarding genuine LOF evidence. Only when there is no
    // coding-null consequence is PVS1 dropped entirely (defer to SpliceAI/PP3).
    if matches!(kind, NullKind::CanonicalSplice) {
        if let Some(offset) = input.intronic_offset {
            details.insert("intronic_offset".into(), serde_json::json!(offset));
            if offset.abs() > 2 {
                match NullKind::detect_non_splice(&input.consequences) {
                    Some(coding_null) => {
                        details.insert(
                            "splice_offset_regraded_to".into(),
                            serde_json::json!(coding_null.label()),
                        );
                        kind = coding_null;
                    }
                    None => {
                        return mk(
                            "PVS1".to_string(),
                            EvidenceStrength::VeryStrong,
                            false,
                            true,
                            format!(
                                "Splice consequence at intronic offset {:+} is outside the canonical ±1/±2 dinucleotide and no coding-null consequence is present → PVS1 not applicable (defer to SpliceAI/PP3)",
                                offset
                            ),
                            details,
                        );
                    }
                }
            }
        }

        // Splice-prediction consistency gate (Walker 2023, ClinGen SVI
        // Splicing Subgroup). PVS1's canonical track assumes the ±1/±2
        // dinucleotide is destroyed. Two situations break that assumption and
        // both showed up in the medical-genetics review:
        //
        //  1. SpliceAI is confidently benign (≤ `spliceai_benign`). A genuine
        //     canonical ±1/±2 change is the easiest call SpliceAI makes, so a
        //     score of ~0 means the positional call is wrong - typically an
        //     indel that merely overlaps the region, or a repeat-context
        //     deletion with an ambiguous alignment. ATM c.?  `GTAATC>G`
        //     (SpliceAI 0.00) and KMT2C `C>CT` (0.05) are both ClinVar
        //     benign/likely-benign and both collected PVS1 this way.
        //  2. The variant is a pure insertion or duplication. Inserting bases
        //     beside, or even inside, the dinucleotide does not necessarily
        //     destroy it: `PTEN c.802-2dupA` and `BRIP1 c.2258-2dup` add the
        //     base the acceptor already carries, so the intron still ends AG.
        //     Without positive splice evidence (SpliceAI ≥ `spliceai_pathogenic`)
        //     PVS1 must not fire on those.
        //
        // In both cases a coding-null consequence on the same allele still
        // carries PVS1 through its own track, so re-grade rather than discard.
        if matches!(kind, NullKind::CanonicalSplice) {
            let spliceai_max = input.splice_ai.as_ref().and_then(|s| s.max_delta_score());
            if let Some(ds) = spliceai_max {
                details.insert("spliceai_max_ds".into(), serde_json::json!(ds));
            }
            let contradicted_by_spliceai =
                spliceai_max.is_some_and(|ds| ds <= config.spliceai_benign);
            let unsupported_insertion = input.is_pure_insertion == Some(true)
                && !spliceai_max.is_some_and(|ds| ds >= config.spliceai_pathogenic);

            if contradicted_by_spliceai || unsupported_insertion {
                let reason = if contradicted_by_spliceai {
                    format!(
                        "SpliceAI max_ds={:.2} ≤ {:.2} contradicts loss of the canonical ±1/±2 site",
                        spliceai_max.unwrap_or(0.0),
                        config.spliceai_benign
                    )
                } else {
                    "pure insertion/duplication at the canonical site may leave the ±1/±2 dinucleotide intact, and no positive SpliceAI support is available".to_string()
                };
                match NullKind::detect_non_splice(&input.consequences) {
                    Some(coding_null) => {
                        details.insert(
                            "splice_evidence_regraded_to".into(),
                            serde_json::json!(coding_null.label()),
                        );
                        kind = coding_null;
                    }
                    None => {
                        details.insert(
                            "splice_evidence_conflict".into(),
                            serde_json::json!(reason.clone()),
                        );
                        return mk(
                            "PVS1".to_string(),
                            EvidenceStrength::VeryStrong,
                            false,
                            true,
                            format!("PVS1 splice track not applicable: {} (defer to PP3/BP4)", reason),
                            details,
                        );
                    }
                }
            }
        }
    }

    // Which NMD prediction the tree runs on. The exact 50-nt rule is the one
    // Abou Tayoun 2018 states; the last-exon proxy is what fastVEP has always
    // used and stays the default because the exact rule trades ClinVar
    // pathogenic recall for faithfulness - see `pvs1_nmd_50nt_rule`.
    let predicted_nmd = if config.pvs1_nmd_50nt_rule {
        input
            .nmd_escape_50nt
            .map(|escapes| !escapes)
            .or(input.predicted_nmd)
    } else {
        input.predicted_nmd
    };

    let (strength, summary) = match kind {
        NullKind::NonsenseOrFrameshift => {
            grade_nonsense_frameshift(input, predicted_nmd, &mut details)
        }
        NullKind::CanonicalSplice => grade_canonical_splice(input, predicted_nmd, &mut details),
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
        // Scan all consequences and pick the most severe null kind so the
        // result doesn't depend on input ordering (e.g. when both
        // splice_donor_variant and frameshift_variant appear, splice wins
        // deterministically). Severity rank below matches Ensembl VEP's
        // null-variant ordering.
        let mut best: Option<Self> = None;
        for c in cs {
            let kind = match c {
                Consequence::TranscriptAblation => Some(Self::WholeGeneDeletion),
                Consequence::SpliceAcceptorVariant | Consequence::SpliceDonorVariant => {
                    Some(Self::CanonicalSplice)
                }
                Consequence::StopGained | Consequence::FrameshiftVariant => {
                    Some(Self::NonsenseOrFrameshift)
                }
                Consequence::StartLost => Some(Self::StartLost),
                _ => None,
            };
            if let Some(k) = kind {
                best = Some(match best {
                    None => k,
                    Some(prev) if k.severity_rank() > prev.severity_rank() => k,
                    Some(prev) => prev,
                });
            }
        }
        best
    }

    /// Detect the most severe NON-splice null kind (whole-gene deletion,
    /// nonsense/frameshift, start-lost). Used when a non-canonical splice
    /// offset disqualifies the splice track but a co-occurring coding-null
    /// consequence still justifies PVS1.
    fn detect_non_splice(cs: &[Consequence]) -> Option<Self> {
        let mut best: Option<Self> = None;
        for c in cs {
            let kind = match c {
                Consequence::TranscriptAblation => Some(Self::WholeGeneDeletion),
                Consequence::StopGained | Consequence::FrameshiftVariant => {
                    Some(Self::NonsenseOrFrameshift)
                }
                Consequence::StartLost => Some(Self::StartLost),
                _ => None,
            };
            if let Some(k) = kind {
                best = Some(match best {
                    None => k,
                    Some(prev) if k.severity_rank() > prev.severity_rank() => k,
                    Some(prev) => prev,
                });
            }
        }
        best
    }

    fn severity_rank(&self) -> u8 {
        match self {
            Self::WholeGeneDeletion => 4,
            Self::CanonicalSplice => 3,
            Self::NonsenseOrFrameshift => 2,
            Self::StartLost => 1,
        }
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
    predicted_nmd: Option<bool>,
    details: &mut serde_json::Map<String, serde_json::Value>,
) -> (EvidenceStrength, String) {
    if let Some(nmd) = predicted_nmd {
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
            _ if input.is_last_exon == Some(true) => (
                EvidenceStrength::Moderate,
                "NMD-escape PTC in last exon; truncation magnitude unknown → PVS1_Moderate (conservative)".to_string(),
            ),
            _ => (
                EvidenceStrength::VeryStrong,
                "NMD-escape; grading signals incomplete → PVS1 (legacy fallback)".to_string(),
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
    predicted_nmd: Option<bool>,
    details: &mut serde_json::Map<String, serde_json::Value>,
) -> (EvidenceStrength, String) {
    if let Some(nmd) = predicted_nmd {
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
        // alt_start_codon_distance is the downstream distance in codons; the
        // pipeline produces a non-negative value. A negative value is
        // out-of-contract — treat it as "no usable signal" rather than
        // silently abs() it (which would let upstream / invalid distances
        // qualify the variant for PVS1_Moderate or _Supporting).
        if d < 0 {
            details.insert("alt_start_codon_distance_invalid".into(), serde_json::json!(true));
        } else {
            let d_codons = d as u64;
            // Note the field is being borrowed for a different question here.
            // The pipeline derives `in_critical_region` as "ClinVar has a
            // pathogenic variant *downstream* of the truncation", which is what
            // the nonsense/frameshift track needs; the start-loss track wants
            // "pathogenic variant upstream of the alternative start". The two
            // only coincide because a start-loss truncates from residue 1. This
            // branch is unreachable today (`alt_start_codon_distance` has no
            // plumbing) and the signal must be split before it becomes live.
            if d_codons <= 100 && input.in_critical_region == Some(true) {
                return (
                    EvidenceStrength::Moderate,
                    format!(
                        "Start-lost with downstream Met {} codons away and pathogenic variant upstream → PVS1_Moderate",
                        d_codons
                    ),
                );
            }
            if d_codons > 100 {
                return (
                    EvidenceStrength::Supporting,
                    format!(
                        "Start-lost; alternative downstream Met is {} codons away → PVS1_Supporting",
                        d_codons
                    ),
                );
            }
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
        if gc.pli.is_some_and(|p| p >= config.pli_lof_intolerant) {
            return true;
        }
        if gc.loeuf.is_some_and(|l| l <= config.loeuf_lof_intolerant) {
            return true;
        }
    }

    // A curated LOF mechanism enables PVS1 even where constraint does not
    // reach the threshold. Resolved through `effective_mechanism` so the
    // shipped table and a user's `gene_overrides` are read the same way here
    // as in the gain-of-function gate above; otherwise a gene could be
    // LOF-enabled by one map and GOF-blocked by the other.
    if let Some(mechanism) = config.effective_mechanism(input.gene_symbol.as_deref()) {
        if mechanism.to_ascii_uppercase().contains("LOF") {
            return true;
        }
    }

    // Disease-gene fallback per ClinGen SVI / Abou Tayoun 2018: when
    // gnomAD constraints don't reach the LOF threshold, accept a
    // curated disease-gene association as evidence the gene is
    // LOF-relevant. The .oga is populated from ClinGen Gene-Disease
    // Validity (preferred — Strong/Definitive/Moderate only) or OMIM
    // `genemap2.txt` (legacy). Both share the `omim` json_key.
    if let Some(ref omim) = input.omim {
        if omim
            .phenotypes
            .as_ref()
            .is_some_and(|p| !p.is_empty())
        {
            return true;
        }
    }

    false
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_support::minimal_input;
    use crate::sa_extract::{GnomadGeneData, OmimData, SpliceAiData};
    use fastvep_core::Impact;

    fn make_input(consequences: Vec<Consequence>, gene_constraints: Option<GnomadGeneData>, omim: Option<OmimData>) -> ClassificationInput {
        ClassificationInput {
            consequences,
            impact: Impact::High,
            gene_symbol: Some("BRCA1".to_string()),
            gene_constraints,
            omim,
            ..minimal_input()
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

    #[test]
    fn test_pvs1_canonical_splice_deep_offset_not_applicable() {
        // splice_donor/acceptor consequence but the HGVS offset is beyond ±2
        // (e.g. an indel spanning into the intron, c.4001+12_4001+15del). PVS1
        // must not fire — splicing signal is left to SpliceAI/PP3.
        let mut input = make_input(
            vec![Consequence::SpliceDonorVariant],
            Some(GnomadGeneData { pli: Some(1.0), loeuf: Some(0.03), ..Default::default() }),
            None,
        );
        input.intronic_offset = Some(12);
        let r = evaluate_pvs1(&input, &AcmgConfig::default());
        assert!(!r.met);
        assert!(r.summary.contains("canonical ±1/±2"));
    }

    #[test]
    fn test_pvs1_canonical_splice_within_2_still_fires() {
        // Canonical ±1/2 splice (offset = -2) keeps PVS1.
        let mut input = make_input(
            vec![Consequence::SpliceAcceptorVariant],
            Some(GnomadGeneData { pli: Some(1.0), loeuf: Some(0.03), ..Default::default() }),
            None,
        );
        input.intronic_offset = Some(-2);
        let r = evaluate_pvs1(&input, &AcmgConfig::default());
        assert!(r.met);
        assert_eq!(r.strength, EvidenceStrength::VeryStrong);
    }

    #[test]
    fn test_pvs1_deep_offset_with_coding_frameshift_still_fires() {
        // An indel called BOTH splice_donor (deep offset >2, e.g. an exon→intron
        // deletion) AND frameshift must keep PVS1 via the frameshift track — the
        // genuine coding LOF must not be discarded by the splice offset gate.
        let mut input = make_input(
            vec![Consequence::SpliceDonorVariant, Consequence::FrameshiftVariant],
            Some(GnomadGeneData { pli: Some(1.0), loeuf: Some(0.03), ..Default::default() }),
            None,
        );
        input.intronic_offset = Some(7); // outside canonical ±1/±2
        let r = evaluate_pvs1(&input, &AcmgConfig::default());
        assert!(r.met, "frameshift LOF should keep PVS1 despite deep splice offset");
        // No NMD/grading signals → frameshift legacy fallback → Very Strong.
        assert_eq!(r.strength, EvidenceStrength::VeryStrong);
    }

    #[test]
    fn test_pvs1_last_exon_nonsense_downgraded_to_moderate() {
        // A PTC in the last exon escapes NMD; with no finer truncation signals
        // it downgrades to PVS1_Moderate instead of the legacy Very Strong.
        let mut input = make_input(
            vec![Consequence::StopGained],
            Some(GnomadGeneData { pli: Some(1.0), loeuf: Some(0.03), ..Default::default() }),
            None,
        );
        input.predicted_nmd = Some(false);
        input.is_last_exon = Some(true);
        let r = evaluate_pvs1(&input, &AcmgConfig::default());
        assert!(r.met);
        assert_eq!(r.strength, EvidenceStrength::Moderate);
        assert_eq!(r.code, "PVS1_Moderate");
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

    fn lof_gene() -> Option<GnomadGeneData> {
        Some(GnomadGeneData { pli: Some(1.0), loeuf: Some(0.03), ..Default::default() })
    }

    #[test]
    fn test_pvs1_splice_dropped_when_spliceai_contradicts() {
        // A genuine canonical +-1/+-2 change is the easiest call SpliceAI
        // makes. A score of ~0 means the positional call is wrong -- typically
        // an indel overlapping the region, or a repeat-context deletion with an
        // ambiguous alignment. ATM GTAATC>G (SpliceAI 0.00) is ClinVar
        // likely-benign and was collecting PVS1 this way.
        let mut input = make_input(vec![Consequence::SpliceDonorVariant], lof_gene(), None);
        input.splice_ai = Some(SpliceAiData { ds_dl: Some(0.0), ..Default::default() });
        let r = evaluate_pvs1(&input, &AcmgConfig::default());
        assert!(!r.met);
        assert!(r.summary.contains("contradicts"));
    }

    #[test]
    fn test_pvs1_splice_kept_when_spliceai_supports() {
        let mut input = make_input(vec![Consequence::SpliceDonorVariant], lof_gene(), None);
        input.splice_ai = Some(SpliceAiData { ds_dl: Some(0.95), ..Default::default() });
        let r = evaluate_pvs1(&input, &AcmgConfig::default());
        assert!(r.met);
    }

    #[test]
    fn test_pvs1_splice_dropped_for_unsupported_insertion() {
        // PTEN c.802-2dupA and BRIP1 c.2258-2dup add the base the acceptor
        // already carries, so the intron still ends AG. Without positive
        // SpliceAI support PVS1 must not fire.
        let mut input = make_input(vec![Consequence::SpliceAcceptorVariant], lof_gene(), None);
        input.is_pure_insertion = Some(true);
        let r = evaluate_pvs1(&input, &AcmgConfig::default());
        assert!(!r.met);
        assert!(r.summary.contains("insertion"));
    }

    #[test]
    fn test_pvs1_splice_insertion_kept_with_spliceai_support() {
        let mut input = make_input(vec![Consequence::SpliceAcceptorVariant], lof_gene(), None);
        input.is_pure_insertion = Some(true);
        input.splice_ai = Some(SpliceAiData { ds_al: Some(0.99), ..Default::default() });
        let r = evaluate_pvs1(&input, &AcmgConfig::default());
        assert!(r.met);
    }

    #[test]
    fn test_pvs1_splice_conflict_regrades_to_coding_null() {
        // An indel that both overlaps the splice site and deletes coding
        // sequence still carries PVS1 through the frameshift track.
        let mut input = make_input(
            vec![Consequence::SpliceDonorVariant, Consequence::FrameshiftVariant],
            lof_gene(),
            None,
        );
        input.splice_ai = Some(SpliceAiData { ds_dl: Some(0.0), ..Default::default() });
        let r = evaluate_pvs1(&input, &AcmgConfig::default());
        assert!(r.met, "coding-null track should still carry PVS1");
    }

    // ── B7: gene-disease validity ────────────────────────────────────────

    /// A frameshift in a maximally constrained gene: everything PVS1 wants,
    /// so whatever stops it here is the gate under test and nothing else.
    fn constrained_frameshift(gene: &str) -> ClassificationInput {
        let mut input = make_input(vec![Consequence::FrameshiftVariant], lof_gene(), None);
        input.gene_symbol = Some(gene.to_string());
        input
    }

    #[test]
    fn test_pvs1_blocked_when_clingen_curated_the_gene_as_invalid() {
        let mut input = constrained_frameshift("RYK");
        input.omim = Some(OmimData {
            mim_number: Some(0),
            phenotypes: Some(vec![
                "some proposed disease (ClinGen Disputed/AD, MONDO:0000001)".into(),
            ]),
        }); // ClinGen curated it and found nothing
        let r = evaluate_pvs1(&input, &AcmgConfig::default());
        assert!(!r.met);
        assert!(!r.evaluated, "unknown gene-disease validity is not an assessment");
        assert!(
            r.summary.contains("no_valid_gene_disease_relationship"),
            "got: {}",
            r.summary
        );
    }

    #[test]
    fn test_pvs1_unaffected_when_no_gene_disease_source_is_loaded() {
        // The back-compat case that makes the gate safe to default on: without
        // an .oga, pLI alone still carries PVS1 exactly as before.
        let input = constrained_frameshift("RYK");
        assert!(input.omim.is_none());
        assert!(evaluate_pvs1(&input, &AcmgConfig::default()).met);
    }

    #[test]
    fn test_pvs1_fires_for_a_gene_the_source_lists() {
        let mut input = constrained_frameshift("BRCA1");
        input.omim = Some(OmimData {
            mim_number: Some(0),
            phenotypes: Some(vec!["hereditary breast cancer (ClinGen Definitive/AD)".into()]),
        });
        assert!(evaluate_pvs1(&input, &AcmgConfig::default()).met);
    }


    #[test]
    fn test_pvs1_survives_for_a_gene_clingen_has_not_curated() {
        // The v10 regression: SPAST, ABCB11, FLG and LAMB3 all cause disease
        // and are all absent from ClinGen GDV. Blocking on absence cost 1,497
        // truth-pathogenic PVS1 firings.
        let mut input = constrained_frameshift("SPAST");
        input.omim = None;
        assert!(evaluate_pvs1(&input, &AcmgConfig::default()).met);
    }

    // ── B6: disease mechanism ────────────────────────────────────────────

    #[test]
    fn test_pvs1_blocked_for_a_gain_of_function_gene() {
        // PCSK9 is the case that shows constraint cannot answer this: it is a
        // constrained gene where the null allele *lowers* LDL. Its pLI would
        // carry PVS1 straight through without the mechanism gate.
        let input = constrained_frameshift("PCSK9");
        let r = evaluate_pvs1(&input, &AcmgConfig::default());
        assert!(!r.met);
        assert!(r.evaluated, "we did assess it; PVS1 is not applicable");
        assert!(
            r.summary.contains("mechanism_not_loss_of_function"),
            "got: {}",
            r.summary
        );
    }

    #[test]
    fn test_pvs1_survives_for_a_gene_with_both_mechanisms() {
        // RYR1: malignant hyperthermia is gain of function, but the congenital
        // myopathies are loss of function, so a null allele is still
        // pathogenic for one of the two diseases.
        assert!(evaluate_pvs1(&constrained_frameshift("RYR1"), &AcmgConfig::default()).met);
    }

    #[test]
    fn test_curated_lof_mechanism_enables_pvs1_without_constraint_data() {
        // The mechanism table has to work in both directions, or adding it
        // would quietly remove the pre-existing "curated LOF enables PVS1"
        // path that read `gene_overrides` alone.
        let mut input = make_input(vec![Consequence::FrameshiftVariant], None, None);
        input.gene_symbol = Some("MYGENE".to_string());
        let mut config = AcmgConfig::default();
        config
            .gene_mechanisms
            .insert("MYGENE".to_string(), "LOF".to_string());
        assert!(evaluate_pvs1(&input, &config).met);
    }

    // ── The 50-nt NMD rule, and what switching it on costs ───────────────

    /// A PTC in the last 50 nt of the penultimate exon: the one place the two
    /// NMD signals disagree. MSH6 `c.3978dup` is the real case - exon 9 of 10,
    /// eight bases from the junction, 2.6 % of the protein removed, and a
    /// region ClinVar has pathogenic variants in.
    fn penultimate_exon_escape() -> ClassificationInput {
        let mut input = make_input(vec![Consequence::FrameshiftVariant], lof_gene(), None);
        input.is_last_exon = Some(false);
        input.predicted_nmd = Some(true); // last-exon proxy: "will decay"
        input.nmd_escape_50nt = Some(true); // measured: escapes
        input.in_critical_region = Some(true);
        input.protein_truncation_pct = Some(0.026);
        input
    }

    #[test]
    fn test_50nt_rule_off_by_default_keeps_very_strong() {
        let r = evaluate_pvs1(&penultimate_exon_escape(), &AcmgConfig::default());
        assert_eq!(r.code, "PVS1");
        assert_eq!(r.strength, EvidenceStrength::VeryStrong);
    }

    #[test]
    fn test_50nt_rule_on_grades_the_escape_down() {
        let config = AcmgConfig { pvs1_nmd_50nt_rule: true, ..Default::default() };
        let r = evaluate_pvs1(&penultimate_exon_escape(), &config);
        assert_eq!(r.code, "PVS1_Strong", "Abou Tayoun: NMD escape in a critical region");
        assert_eq!(r.strength, EvidenceStrength::Strong);
    }

    #[test]
    fn test_50nt_rule_falls_back_to_the_proxy_when_unmeasured() {
        // Intronic variants have no cDNA coordinate, so the measurement is
        // absent and the proxy must still answer.
        let mut input = penultimate_exon_escape();
        input.nmd_escape_50nt = None;
        let config = AcmgConfig { pvs1_nmd_50nt_rule: true, ..Default::default() };
        assert_eq!(evaluate_pvs1(&input, &config).strength, EvidenceStrength::VeryStrong);
    }

    #[test]
    fn test_50nt_rule_does_not_disturb_a_mid_gene_ptc() {
        // The overwhelming majority of null variants: both signals agree that
        // the message decays, so the flag changes nothing either way.
        let mut input = make_input(vec![Consequence::StopGained], lof_gene(), None);
        input.is_last_exon = Some(false);
        input.predicted_nmd = Some(true);
        input.nmd_escape_50nt = Some(false);
        for on in [false, true] {
            let config = AcmgConfig { pvs1_nmd_50nt_rule: on, ..Default::default() };
            let r = evaluate_pvs1(&input, &config);
            assert_eq!(r.code, "PVS1", "flag={on}");
        }
    }

    #[test]
    fn test_mechanism_gate_can_be_switched_off() {
        let config = AcmgConfig { mechanism_gates_pvs1: false, ..Default::default() };
        assert!(evaluate_pvs1(&constrained_frameshift("PCSK9"), &config).met);
    }
}
