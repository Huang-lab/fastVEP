pub(crate) mod frequency_gate;
pub(crate) use frequency_gate::{
    FREQUENCY_BLOCKED_SITE_LEVEL, FREQUENCY_BLOCKED_WOULD_BE_BENIGN,
};
pub(crate) mod gene_disease;
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
    gene.is_some_and(|g| config.is_criterion_disabled(g, code))
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
///
/// It also applies one *upgrade*, which is not suppression but needs the same
/// cross-criterion view:
///
/// 5. **PS1 (splice) graded against PVS1** - Walker 2023 Table 3 sets PS1's
///    strength for a canonical-dinucleotide variant from the PVS1 code the
///    variant already carries: `PS1_Supporting` under a full-strength PVS1, and
///    full `PS1` under a downgraded one, so that a variant whose null evidence
///    was reduced is not also denied the borrowed clinical evidence.
fn reconcile_evidence(criteria: &mut [EvidenceCriterion]) {
    grade_ps1_splice_against_pvs1(criteria);

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

    // Rule 0: PVS1 + PM1 - a null variant already carries the strongest
    // possible statement about the region it destroys, so adding residue-level
    // hotspot evidence on top double-counts. The round-2 medical-genetics
    // review raised this directly on CBS, MSH6 and RYR1: "PM1 is called with
    // PVS1? Where is the evidence for PM1 coming?". Applied before the PP3
    // rules because it does not depend on PP3 firing.
    if pvs1_met && pm1_met {
        if let Some(pm1_i) = pm1_idx {
            suppress(
                &mut criteria[pm1_i],
                "Suppressed: PVS1 already counts the loss of this region; PM1 residue-level hotspot evidence would double-count (Abou Tayoun 2018).",
            );
            pm1_met = false;
        }
    }

    // Rule 0b: curated functional evidence outranks every computational
    // prediction. The ClinGen SVI ordering is explicit - a well-established
    // assay is stronger evidence than a predictor, and PP3/BP4/BP7 must not be
    // used to argue against sound experimental data. The round-2 review's OCA2
    // c.1503A>G is the case: a synonymous variant with published splice-defect
    // data, where fastVEP fired BP7 and BP4 and called it benign.
    //
    // Both directions of prediction are suppressed by evidence in either
    // direction. That is deliberate. Where an assay says "damaging" and REVEL
    // says "benign", the predictor is not a dissenting vote to be counted, it
    // is superseded; and the reverse holds just as strongly.
    if criteria.iter().any(|c| c.met && (c.code == "PS3" || c.code == "BS3")) {
        let evidence = criteria
            .iter()
            .find(|c| c.met && (c.code == "PS3" || c.code == "BS3"))
            .map(|c| c.code.clone())
            .unwrap_or_default();
        let reason = format!(
            "Superseded by functional evidence: {} rests on a curated assay, and ClinGen SVI ranks experimental evidence above computational prediction.",
            evidence
        );
        let mut suppressed_pp3 = false;
        for c in criteria.iter_mut() {
            if !c.met {
                continue;
            }
            if c.code.starts_with("PP3") || c.code.starts_with("BP4") || c.code.starts_with("BP7") {
                suppress(c, &reason);
                suppressed_pp3 |= c.code.starts_with("PP3");
            }
        }
        // PP3 is gone, so every rule below it has nothing left to reconcile.
        if suppressed_pp3 {
            return;
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

/// Build PS3 or BS3 from the run's curated functional evidence, if it names
/// this variant.
///
/// Shared by both because they are the same lookup in opposite directions, and
/// because the strength has to come from the file in both cases: Brnich 2020
/// treats assay strength as a judgement about the assay, not a constant.
fn functional_criterion(
    input: &ClassificationInput,
    want: crate::functional::FunctionalCriterion,
    direction: crate::types::EvidenceDirection,
    absent_summary: &str,
) -> EvidenceCriterion {
    let code = want.code().to_string();
    let entry = input
        .functional_evidence
        .as_ref()
        .filter(|e| e.criterion == want);

    let Some(entry) = entry else {
        return EvidenceCriterion {
            code,
            direction,
            strength: EvidenceStrength::Strong,
            default_strength: EvidenceStrength::Strong,
            met: false,
            evaluated: false,
            summary: absent_summary.to_string(),
            details: serde_json::Value::Null,
        };
    };

    let mut details = serde_json::Map::new();
    details.insert("source".into(), serde_json::json!("curated functional evidence"));
    details.insert("strength".into(), serde_json::json!(entry.strength.as_str()));
    if let Some(ref pmid) = entry.pmid {
        details.insert("pmid".into(), serde_json::json!(pmid));
    }
    if let Some(ref note) = entry.note {
        details.insert("note".into(), serde_json::json!(note));
    }

    // Cite the source in the summary itself. A curator reading a report should
    // not have to open the details blob to find out which paper this rests on.
    let mut summary = format!(
        "Curated functional evidence: {} at {} strength",
        if want == crate::functional::FunctionalCriterion::Ps3 {
            "assay shows a damaging effect"
        } else {
            "assay shows no damaging effect"
        },
        entry.strength.as_str(),
    );
    if let Some(ref pmid) = entry.pmid {
        summary.push_str(&format!(" (PMID {})", pmid));
    }
    if let Some(ref note) = entry.note {
        summary.push_str(&format!(": {}", note));
    }

    EvidenceCriterion {
        code,
        direction,
        strength: entry.strength,
        default_strength: EvidenceStrength::Strong,
        met: true,
        evaluated: true,
        summary,
        details: serde_json::Value::Object(details),
    }
}

/// Apply Walker 2023 Table 3's PVS1-dependent grading of PS1's splice path.
///
/// The splice evaluator emits the conservative `PS1_Supporting` row, which is
/// correct whenever PVS1 stands at full Very Strong. Table 3's fifth row raises
/// it to full `PS1` when the variant's PVS1 was itself downgraded - the
/// reduction exists only to stop the variant under assessment out-scoring the
/// pathogenic variant it borrows evidence from, and a downgraded PVS1 leaves
/// room for it.
///
/// A variant with no PVS1 at all (a gene with no established LOF mechanism, say)
/// falls under no row of Table 3, and keeps the Supporting default.
fn grade_ps1_splice_against_pvs1(criteria: &mut [EvidenceCriterion]) {
    let pvs1_strength = criteria
        .iter()
        .find(|c| c.met && (c.code == "PVS1" || c.code.starts_with("PVS1_")))
        .map(|c| c.strength);
    let Some(pvs1_strength) = pvs1_strength else {
        return;
    };
    if pvs1_strength >= EvidenceStrength::VeryStrong {
        return;
    }
    for c in criteria.iter_mut() {
        let is_splice_ps1 = c.met
            && c.code.starts_with("PS1")
            && c.details.get("ps1_path").and_then(|v| v.as_str()) == Some("splice_rna_match");
        if !is_splice_ps1 {
            continue;
        }
        c.code = "PS1".to_string();
        c.strength = EvidenceStrength::Strong;
        c.summary = format!(
            "{}; raised to Strong because PVS1 is downgraded to {} on this variant (Walker 2023 Table 3)",
            c.summary,
            pvs1_strength.as_str()
        );
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
    use crate::test_support::minimal_input;
    use crate::types::{EvidenceDirection, EvidenceStrength};
    use fastvep_core::Consequence;

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

    fn not_met(code: &str, dir: EvidenceDirection, strength: EvidenceStrength) -> EvidenceCriterion {
        EvidenceCriterion { met: false, ..met(code, dir, strength) }
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

    // ── B8: curated functional evidence ──────────────────────────────────

    #[test]
    fn functional_evidence_supersedes_every_predictor() {
        // The SVI ordering: PP3/BP4/BP7 must not be used to argue against
        // sound experimental data. OCA2 c.1503A>G is the case - a synonymous
        // variant with published splice-defect data, where fastVEP was firing
        // BP7 and BP4 and calling it benign.
        let mut criteria = vec![
            met("PS3", EvidenceDirection::Pathogenic, EvidenceStrength::Strong),
            pp3_with_source(EvidenceStrength::Strong, "revel_missense"),
            met("BP4", EvidenceDirection::Benign, EvidenceStrength::Supporting),
            met("BP7", EvidenceDirection::Benign, EvidenceStrength::Supporting),
            met("PM2_Supporting", EvidenceDirection::Pathogenic, EvidenceStrength::Supporting),
        ];
        reconcile_evidence(&mut criteria);

        assert!(find(&criteria, "PS3").met, "the evidence itself must survive");
        for code in ["PP3", "BP4", "BP7"] {
            let c = find(&criteria, code);
            assert!(!c.met, "{code} should be superseded");
            assert!(c.summary.contains("Superseded by functional evidence"), "got: {}", c.summary);
        }
        assert!(
            find(&criteria, "PM2_Supporting").met,
            "only the computational predictors are superseded, not the frequency evidence",
        );
    }

    #[test]
    fn benign_functional_evidence_also_supersedes_pathogenic_prediction() {
        // Symmetry matters: an assay showing no damaging effect outranks a
        // high REVEL exactly as a damaging assay outranks a low one.
        let mut criteria = vec![
            met("BS3", EvidenceDirection::Benign, EvidenceStrength::Strong),
            pp3_with_source(EvidenceStrength::Strong, "revel_missense"),
        ];
        reconcile_evidence(&mut criteria);
        assert!(find(&criteria, "BS3").met);
        assert!(!find(&criteria, "PP3").met);
    }

    #[test]
    fn predictors_stand_when_there_is_no_functional_evidence() {
        // A PS3 that was evaluated but not met must not suppress anything.
        let mut criteria = vec![
            not_met("PS3", EvidenceDirection::Pathogenic, EvidenceStrength::Strong),
            met("BP4", EvidenceDirection::Benign, EvidenceStrength::Supporting),
        ];
        reconcile_evidence(&mut criteria);
        assert!(find(&criteria, "BP4").met);
    }

    /// The shared implementation behind PS3 and BS3. `pathogenic_strong` and
    /// `benign_strong` each wrap it in a one-liner, so exercising it here
    /// covers both.
    fn ps3(input: &ClassificationInput) -> EvidenceCriterion {
        functional_criterion(
            input,
            crate::functional::FunctionalCriterion::Ps3,
            EvidenceDirection::Pathogenic,
            "Requires curated functional study evidence (in vitro/in vivo assays) — supply one with --functional-evidence",
        )
    }

    fn bs3(input: &ClassificationInput) -> EvidenceCriterion {
        functional_criterion(
            input,
            crate::functional::FunctionalCriterion::Bs3,
            EvidenceDirection::Benign,
            "Requires curated functional study evidence showing no damaging effect — supply one with --functional-evidence",
        )
    }

    fn with_functional(
        criterion: crate::functional::FunctionalCriterion,
        strength: EvidenceStrength,
    ) -> ClassificationInput {
        let mut input = minimal_input();
        input.consequences = vec![Consequence::MissenseVariant];
        input.functional_evidence = Some(crate::functional::FunctionalEvidence {
            criterion,
            strength,
            pmid: Some("29625052".to_string()),
            note: Some("minigene shows exon skipping".to_string()),
        });
        input
    }

    #[test]
    fn test_functional_criteria_are_not_evaluated_without_a_file() {
        // The honest default: no curated evidence means no assertion either
        // way, not an absence of effect.
        let input = minimal_input();
        for c in [
            ps3(&input),
            bs3(&input),
        ] {
            assert!(!c.met);
            assert!(!c.evaluated);
            assert!(c.summary.contains("--functional-evidence"), "got: {}", c.summary);
        }
    }

    #[test]
    fn test_ps3_fires_from_curated_evidence_and_cites_it() {
        let input = with_functional(
            crate::functional::FunctionalCriterion::Ps3,
            EvidenceStrength::Strong,
        );
        let c = ps3(&input);
        assert!(c.met && c.evaluated);
        assert_eq!(c.strength, EvidenceStrength::Strong);
        assert!(c.summary.contains("PMID 29625052"), "got: {}", c.summary);
        assert!(c.summary.contains("exon skipping"), "got: {}", c.summary);
    }

    #[test]
    fn test_curator_assigned_strength_is_honoured() {
        // Brnich 2020: a weak assay supports weak evidence. The file says so
        // and the criterion has to carry it, or the strength column is a lie.
        let input = with_functional(
            crate::functional::FunctionalCriterion::Ps3,
            EvidenceStrength::Supporting,
        );
        let c = ps3(&input);
        assert!(c.met);
        assert_eq!(c.strength, EvidenceStrength::Supporting);
        assert_eq!(c.default_strength, EvidenceStrength::Strong);
    }

    #[test]
    fn test_bs3_evidence_does_not_fire_ps3() {
        // The two criteria read the same field and must not be confused.
        let input = with_functional(
            crate::functional::FunctionalCriterion::Bs3,
            EvidenceStrength::Strong,
        );
        assert!(!ps3(&input).met);
        assert!(bs3(&input).met);
    }

    // ── Walker 2023 Table 3: PS1's splice path graded against PVS1 ───────

    fn splice_ps1_supporting() -> EvidenceCriterion {
        let mut c = met("PS1_Supporting", EvidenceDirection::Pathogenic, EvidenceStrength::Supporting);
        c.details = serde_json::json!({"ps1_path": "splice_rna_match"});
        c
    }

    #[test]
    fn test_splice_ps1_stays_supporting_under_a_full_strength_pvs1() {
        let mut criteria = vec![
            met("PVS1", EvidenceDirection::Pathogenic, EvidenceStrength::VeryStrong),
            splice_ps1_supporting(),
        ];
        grade_ps1_splice_against_pvs1(&mut criteria);
        assert_eq!(criteria[1].code, "PS1_Supporting");
        assert_eq!(criteria[1].strength, EvidenceStrength::Supporting);
    }

    #[test]
    fn test_splice_ps1_rises_to_strong_under_a_downgraded_pvs1() {
        for downgraded in [
            EvidenceStrength::Strong,
            EvidenceStrength::Moderate,
            EvidenceStrength::Supporting,
        ] {
            let mut criteria = vec![
                met(
                    &format!("PVS1_{}", downgraded.as_str()),
                    EvidenceDirection::Pathogenic,
                    downgraded,
                ),
                splice_ps1_supporting(),
            ];
            grade_ps1_splice_against_pvs1(&mut criteria);
            assert_eq!(criteria[1].code, "PS1", "PVS1_{:?}", downgraded);
            assert_eq!(criteria[1].strength, EvidenceStrength::Strong, "PVS1_{:?}", downgraded);
        }
    }

    #[test]
    fn test_splice_ps1_stays_supporting_with_no_pvs1_at_all() {
        // No row of Table 3 covers a canonical splice variant carrying no PVS1
        // (a gene with no established LOF mechanism, say), so the conservative
        // strength stands.
        let mut criteria = vec![
            not_met("PVS1", EvidenceDirection::Pathogenic, EvidenceStrength::VeryStrong),
            splice_ps1_supporting(),
        ];
        grade_ps1_splice_against_pvs1(&mut criteria);
        assert_eq!(criteria[1].code, "PS1_Supporting");
    }

    #[test]
    fn test_missense_ps1_is_untouched_by_the_splice_grading() {
        // The missense path is Richards 2015 PS1 at full Strong and has
        // nothing to do with Table 3.
        let mut ps1 = met("PS1", EvidenceDirection::Pathogenic, EvidenceStrength::Strong);
        ps1.details = serde_json::json!({"ps1_path": "missense_aa_match"});
        let mut criteria = vec![
            met("PVS1_Moderate", EvidenceDirection::Pathogenic, EvidenceStrength::Moderate),
            ps1,
        ];
        grade_ps1_splice_against_pvs1(&mut criteria);
        assert_eq!(criteria[1].strength, EvidenceStrength::Strong);
        assert!(criteria[1].summary.is_empty(), "{}", criteria[1].summary);
    }
}
