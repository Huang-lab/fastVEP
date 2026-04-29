pub mod pvs1;
pub mod pathogenic_strong;
pub mod pathogenic_moderate;
pub mod pathogenic_supporting;
pub mod benign_standalone;
pub mod benign_strong;
pub mod benign_supporting;

use crate::config::AcmgConfig;
use crate::sa_extract::ClassificationInput;
use crate::types::{EvidenceCriterion, EvidenceStrength};

/// Evaluate all 28 ACMG-AMP criteria and return the full list.
///
/// After per-criterion evaluation, a reconciliation pass suppresses
/// computational evidence (PP3/BP4) that would double-count with criteria
/// covering the same molecular signal — see `reconcile_evidence`.
pub fn evaluate_all_criteria(
    input: &ClassificationInput,
    config: &AcmgConfig,
) -> Vec<EvidenceCriterion> {
    let gene = input.gene_symbol.as_deref();

    let mut criteria = Vec::with_capacity(28);

    // Pathogenic Very Strong
    let pvs1 = pvs1::evaluate_pvs1(input, config);
    if !is_disabled(gene, &pvs1.code, config) {
        criteria.push(pvs1);
    }

    // Pathogenic Strong
    for c in pathogenic_strong::evaluate_all(input, config) {
        if !is_disabled(gene, &c.code, config) {
            criteria.push(c);
        }
    }

    // Pathogenic Moderate
    for c in pathogenic_moderate::evaluate_all(input, config) {
        if !is_disabled(gene, &c.code, config) {
            criteria.push(c);
        }
    }

    // Pathogenic Supporting
    for c in pathogenic_supporting::evaluate_all(input, config) {
        if !is_disabled(gene, &c.code, config) {
            criteria.push(c);
        }
    }

    // Benign Standalone
    let ba1 = benign_standalone::evaluate_ba1(input, config);
    if !is_disabled(gene, &ba1.code, config) {
        criteria.push(ba1);
    }

    // Benign Strong
    for c in benign_strong::evaluate_all(input, config) {
        if !is_disabled(gene, &c.code, config) {
            criteria.push(c);
        }
    }

    // Benign Supporting
    for c in benign_supporting::evaluate_all(input, config) {
        if !is_disabled(gene, &c.code, config) {
            criteria.push(c);
        }
    }

    reconcile_evidence(&mut criteria);

    criteria
}

fn is_disabled(gene: Option<&str>, code: &str, config: &AcmgConfig) -> bool {
    gene.map_or(false, |g| config.is_criterion_disabled(g, code))
}

/// Suppress double-counted evidence per ClinGen SVI guidance.
///
/// Why this exists: PP3/BP4 are *computational* surrogates for molecular
/// effects that other criteria capture more directly. Letting both fire
/// inflates the evidence weight for the same biological signal. Pejaver 2022
/// (PP3/BP4 calibration) and Walker 2023 (splicing) explicitly call this out.
///
/// Rules applied (each leaves a trail in the affected criterion's `summary`
/// and `details.suppressed_by` so downstream consumers see why a code that
/// could have fired didn't):
///
/// 1. **PVS1 + PP3 (splice)** — if PVS1 fires for a canonical splice variant
///    and PP3 was met from SpliceAI, the splicing signal is already counted.
///    PP3 is suppressed (Walker 2023).
///
/// 2. **PS1 + PP3 (REVEL)** — PS1 means a known pathogenic missense exists at
///    the same residue + same alt AA; that's stronger residue-level evidence
///    than REVEL. PP3 is suppressed when both fire on missense (Pejaver 2022).
///
/// 3. **PM5 + PP3 (REVEL)** — same logic as PS1: PM5 covers a different alt
///    AA at a known pathogenic residue. PP3 is suppressed (Pejaver 2022).
///
/// 4. **PP3 + PM1 strength cap** — Pejaver 2022 recommends capping the
///    combined evidence weight at Strong (4 Tavtigian points). When PP3 is
///    elevated to Strong (the only case where the combined points exceed 4),
///    PM1 is suppressed. PP3 at Moderate or Supporting + PM1 stays within the
///    cap and both fire.
fn reconcile_evidence(criteria: &mut [EvidenceCriterion]) {
    // Collect the firing state of each criterion of interest before we mutate
    // anything. Using indices avoids borrow issues with mutable iteration.
    let mut pvs1_met = false;
    let mut ps1_met = false;
    let mut pm5_met = false;
    let mut pm1_met = false;
    let mut pp3_idx: Option<usize> = None;
    let mut pm1_idx: Option<usize> = None;

    for (i, c) in criteria.iter().enumerate() {
        if !c.met {
            continue;
        }
        // Match by code prefix so any graded variants (e.g. PVS1_Strong from
        // the Abou Tayoun decision tree) participate in reconciliation. An
        // exact match would silently miss graded codes once the pipeline
        // populates the PVS1 grading signals.
        let code = c.code.as_str();
        if code == "PVS1" || code.starts_with("PVS1_") {
            pvs1_met = true;
        } else if code == "PS1" || code.starts_with("PS1_") {
            ps1_met = true;
        } else if code == "PM5" || code.starts_with("PM5_") {
            pm5_met = true;
        } else if code == "PM1" || code.starts_with("PM1_") {
            pm1_met = true;
            pm1_idx = Some(i);
        } else if code.starts_with("PP3") {
            pp3_idx = Some(i);
        }
    }

    let Some(pp3_i) = pp3_idx else {
        // No PP3 firing — nothing to reconcile on the pathogenic computational side.
        return;
    };

    // Inspect what evidence actually drove PP3 firing. The PP3 evaluator
    // stamps `details.pp3_source` with one of "revel_missense" / "spliceai" /
    // "none" so we don't need to second-guess from raw scores.
    let pp3_source = criteria[pp3_i]
        .details
        .get("pp3_source")
        .and_then(|v| v.as_str())
        .unwrap_or("");
    let pp3_from_revel = pp3_source == "revel_missense";
    let pp3_from_splice = pp3_source == "spliceai";
    let pp3_strength = criteria[pp3_i].strength;

    // Rule 1: PVS1 + PP3(splice) — drop PP3.
    if pvs1_met && pp3_from_splice {
        suppress(
            &mut criteria[pp3_i],
            "Suppressed: PVS1 already counts the splicing signal that drove PP3 (Walker 2023).",
        );
        return;
    }

    // Rules 2 & 3: PS1 / PM5 + PP3(REVEL) — drop PP3.
    if pp3_from_revel && (ps1_met || pm5_met) {
        let reason = if ps1_met && pm5_met {
            "Suppressed: PS1 and PM5 already count the residue-level pathogenicity that REVEL captures (Pejaver 2022)."
        } else if ps1_met {
            "Suppressed: PS1 already counts the residue-level pathogenicity that REVEL captures (Pejaver 2022)."
        } else {
            "Suppressed: PM5 already counts the residue-level pathogenicity that REVEL captures (Pejaver 2022)."
        };
        suppress(&mut criteria[pp3_i], reason);
        return;
    }

    // Rule 4: PP3 + PM1 cap — when PP3 is elevated to Strong, drop PM1 to keep
    // the combined evidence ≤ Strong (4 Tavtigian points). Lower PP3 strengths
    // already sit within the cap and both can stand.
    if pm1_met && pp3_strength == EvidenceStrength::Strong {
        if let Some(pm1_i) = pm1_idx {
            suppress(
                &mut criteria[pm1_i],
                "Suppressed: PP3_Strong + PM1 would exceed the Strong combined cap (Pejaver 2022); keeping the stronger code.",
            );
        }
    }
}

fn suppress(c: &mut EvidenceCriterion, reason: &str) {
    c.met = false;
    c.summary = format!("{} (was {})", reason, c.summary);
    if let serde_json::Value::Object(ref mut map) = c.details {
        map.insert("suppressed_by_reconcile".into(), serde_json::json!(reason));
    } else {
        let mut map = serde_json::Map::new();
        map.insert("suppressed_by_reconcile".into(), serde_json::json!(reason));
        c.details = serde_json::Value::Object(map);
    }
}

#[cfg(test)]
mod reconcile_tests {
    use super::*;
    use crate::types::{EvidenceDirection, EvidenceStrength};

    fn met(code: &str, dir: EvidenceDirection, strength: EvidenceStrength) -> EvidenceCriterion {
        EvidenceCriterion {
            code: code.to_string(),
            direction: dir,
            strength,
            default_strength: strength,
            met: true,
            evaluated: true,
            summary: String::new(),
            details: serde_json::Value::Object(serde_json::Map::new()),
        }
    }

    fn pp3_with_source(strength: EvidenceStrength, source: &str) -> EvidenceCriterion {
        let code = if strength == EvidenceStrength::Supporting {
            "PP3".to_string()
        } else {
            format!("PP3_{}", strength.as_str())
        };
        let mut details = serde_json::Map::new();
        details.insert("pp3_source".into(), serde_json::json!(source));
        EvidenceCriterion {
            code,
            direction: EvidenceDirection::Pathogenic,
            strength,
            default_strength: EvidenceStrength::Supporting,
            met: true,
            evaluated: true,
            summary: "test".to_string(),
            details: serde_json::Value::Object(details),
        }
    }

    fn find<'a>(criteria: &'a [EvidenceCriterion], code_prefix: &str) -> &'a EvidenceCriterion {
        criteria
            .iter()
            .find(|c| c.code.starts_with(code_prefix))
            .expect("criterion missing")
    }

    #[test]
    fn pvs1_plus_pp3_splice_suppresses_pp3() {
        // Walker 2023: PVS1 already counts the splicing signal; PP3 from
        // SpliceAI must not double-count.
        let mut criteria = vec![
            met("PVS1", EvidenceDirection::Pathogenic, EvidenceStrength::VeryStrong),
            pp3_with_source(EvidenceStrength::Supporting, "spliceai"),
        ];
        reconcile_evidence(&mut criteria);
        assert!(find(&criteria, "PVS1").met, "PVS1 should remain met");
        assert!(!find(&criteria, "PP3").met, "PP3(splice) should be suppressed");
    }

    #[test]
    fn pvs1_does_not_suppress_pp3_revel() {
        // PP3 driven by REVEL (missense) is unrelated to splicing — PVS1 firing
        // for some other null variant in the same call must not suppress it.
        let mut criteria = vec![
            met("PVS1", EvidenceDirection::Pathogenic, EvidenceStrength::VeryStrong),
            pp3_with_source(EvidenceStrength::Strong, "revel_missense"),
        ];
        reconcile_evidence(&mut criteria);
        assert!(find(&criteria, "PP3").met, "PP3(REVEL) should not be suppressed by PVS1");
    }

    #[test]
    fn ps1_plus_pp3_revel_suppresses_pp3() {
        // Pejaver 2022: PS1 already covers the residue-level pathogenicity
        // that REVEL captures.
        let mut criteria = vec![
            met("PS1", EvidenceDirection::Pathogenic, EvidenceStrength::Strong),
            pp3_with_source(EvidenceStrength::Strong, "revel_missense"),
        ];
        reconcile_evidence(&mut criteria);
        assert!(find(&criteria, "PS1").met);
        assert!(!find(&criteria, "PP3").met, "PP3(REVEL) should be suppressed by PS1");
    }

    #[test]
    fn pm5_plus_pp3_revel_suppresses_pp3() {
        let mut criteria = vec![
            met("PM5", EvidenceDirection::Pathogenic, EvidenceStrength::Moderate),
            pp3_with_source(EvidenceStrength::Moderate, "revel_missense"),
        ];
        reconcile_evidence(&mut criteria);
        assert!(find(&criteria, "PM5").met);
        assert!(!find(&criteria, "PP3").met, "PP3(REVEL) should be suppressed by PM5");
    }

    #[test]
    fn ps1_does_not_suppress_pp3_splice() {
        // PS1 is missense; if PP3 came from SpliceAI it's a different signal.
        let mut criteria = vec![
            met("PS1", EvidenceDirection::Pathogenic, EvidenceStrength::Strong),
            pp3_with_source(EvidenceStrength::Supporting, "spliceai"),
        ];
        reconcile_evidence(&mut criteria);
        assert!(find(&criteria, "PP3").met, "PP3(splice) should not be suppressed by PS1");
    }

    #[test]
    fn pp3_strong_plus_pm1_caps_pm1() {
        // Pejaver 2022: limit PP3+PM1 combined to Strong (4 points).
        // PP3_Strong (4) + PM1 (2) = 6 > 4 → drop PM1.
        let mut criteria = vec![
            pp3_with_source(EvidenceStrength::Strong, "revel_missense"),
            met("PM1", EvidenceDirection::Pathogenic, EvidenceStrength::Moderate),
        ];
        reconcile_evidence(&mut criteria);
        assert!(find(&criteria, "PP3").met, "PP3_Strong should remain met");
        assert!(!find(&criteria, "PM1").met, "PM1 should be suppressed under PP3_Strong cap");
    }

    #[test]
    fn pp3_moderate_plus_pm1_keeps_both() {
        // PP3_Moderate (2) + PM1 (2) = 4 → at the Strong cap, both can stand.
        let mut criteria = vec![
            pp3_with_source(EvidenceStrength::Moderate, "revel_missense"),
            met("PM1", EvidenceDirection::Pathogenic, EvidenceStrength::Moderate),
        ];
        reconcile_evidence(&mut criteria);
        assert!(find(&criteria, "PP3").met);
        assert!(find(&criteria, "PM1").met);
    }

    #[test]
    fn pp3_supporting_plus_pm1_keeps_both() {
        // PP3 (1) + PM1 (2) = 3 < 4 → both stand.
        let mut criteria = vec![
            pp3_with_source(EvidenceStrength::Supporting, "revel_missense"),
            met("PM1", EvidenceDirection::Pathogenic, EvidenceStrength::Moderate),
        ];
        reconcile_evidence(&mut criteria);
        assert!(find(&criteria, "PP3").met);
        assert!(find(&criteria, "PM1").met);
    }

    #[test]
    fn graded_pvs1_codes_still_suppress_pp3_splice() {
        // Regression guard: an exact-match `"PVS1"` check would silently miss
        // graded codes from the Abou Tayoun decision tree (PVS1_Strong /
        // PVS1_Moderate / PVS1_Supporting). The reconcile pass uses prefix
        // matching so PVS1_* + PP3(splice) still drops PP3 (Walker 2023).
        for graded_code in ["PVS1_Strong", "PVS1_Moderate", "PVS1_Supporting"] {
            let mut criteria = vec![
                met(graded_code, EvidenceDirection::Pathogenic, EvidenceStrength::Strong),
                pp3_with_source(EvidenceStrength::Supporting, "spliceai"),
            ];
            reconcile_evidence(&mut criteria);
            assert!(
                !find(&criteria, "PP3").met,
                "PP3(splice) should be suppressed when {} fires",
                graded_code
            );
        }
    }

    #[test]
    fn graded_pm1_code_still_capped_by_pp3_strong() {
        // Same regression guard for PM1 graded codes (if ever introduced).
        let mut criteria = vec![
            pp3_with_source(EvidenceStrength::Strong, "revel_missense"),
            met("PM1_Strong", EvidenceDirection::Pathogenic, EvidenceStrength::Strong),
        ];
        reconcile_evidence(&mut criteria);
        assert!(!find(&criteria, "PM1").met, "graded PM1 should be capped under PP3_Strong");
    }
}
