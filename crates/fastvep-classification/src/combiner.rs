use crate::types::{AcmgClassification, EvidenceCounts, EvidenceCriterion};

/// Apply ACMG-AMP combination rules to determine final classification.
///
/// Returns the classification and the name of the triggered rule.
///
/// Rules: pathogenic (8 combinations) and benign (BA1 / >=2 BS / BS+BP /
/// >=2 BP) are evaluated independently. If both directions reach a
/// definite call (P/LP and B/LB), the result is VUS (Conflicting).
/// Otherwise the directional call wins.
///
/// Includes the ClinGen SVI novel combination rule (Sept 2020) that
/// PVS + >=1 PP → Likely Pathogenic, added to compensate for the PM2
/// downgrade to Supporting (Bayesian Post_P = 0.988, within LP range).
pub fn combine(criteria: &[EvidenceCriterion]) -> (AcmgClassification, Option<String>) {
    let counts = EvidenceCounts::from_criteria(criteria);
    let pvs = counts.pathogenic_very_strong;
    let ps = counts.pathogenic_strong;
    let pm = counts.pathogenic_moderate;
    let pp = counts.pathogenic_supporting;
    let ba = counts.benign_standalone;
    let bs = counts.benign_strong;
    let bp = counts.benign_supporting;

    let pathogenic_call = compute_pathogenic_call(pvs, ps, pm, pp);
    let benign_call = compute_benign_call(ba, bs, bp);

    match (pathogenic_call, benign_call) {
        // Both directions reach a definite call → conflict.
        (Some((p_cls, p_rule)), Some((b_cls, b_rule)))
            if is_definite(p_cls) && is_definite(b_cls) =>
        {
            (
                AcmgClassification::UncertainSignificance,
                Some(format!(
                    "Conflicting evidence: pathogenic rules → {} ({}), benign rules → {} ({})",
                    p_cls.shorthand(),
                    p_rule,
                    b_cls.shorthand(),
                    b_rule
                )),
            )
        }
        // Only pathogenic side reaches a call.
        (Some((cls, rule)), _) => (cls, Some(rule)),
        // Only benign side reaches a call.
        (None, Some((cls, rule))) => (cls, Some(rule)),
        // Neither side fired — VUS.
        (None, None) => (AcmgClassification::UncertainSignificance, None),
    }
}

/// Returns the strongest pathogenic call (or None) along with the rule name.
fn compute_pathogenic_call(
    pvs: u8,
    ps: u8,
    pm: u8,
    pp: u8,
) -> Option<(AcmgClassification, String)> {
    // ── Pathogenic (8 combinations, most specific first) ──
    if pvs >= 1 && ps >= 1 {
        return Some((AcmgClassification::Pathogenic, "PVS + >=1 PS".to_string()));
    }
    if pvs >= 1 && pm >= 2 {
        return Some((AcmgClassification::Pathogenic, "PVS + >=2 PM".to_string()));
    }
    if pvs >= 1 && pm >= 1 && pp >= 1 {
        return Some((AcmgClassification::Pathogenic, "PVS + PM + PP".to_string()));
    }
    if pvs >= 1 && pp >= 2 {
        return Some((AcmgClassification::Pathogenic, "PVS + >=2 PP".to_string()));
    }
    if ps >= 2 {
        return Some((AcmgClassification::Pathogenic, ">=2 PS".to_string()));
    }
    if ps >= 1 && pm >= 3 {
        return Some((AcmgClassification::Pathogenic, "PS + >=3 PM".to_string()));
    }
    if ps >= 1 && pm >= 2 && pp >= 2 {
        return Some((
            AcmgClassification::Pathogenic,
            "PS + 2 PM + >=2 PP".to_string(),
        ));
    }
    if ps >= 1 && pm >= 1 && pp >= 4 {
        return Some((
            AcmgClassification::Pathogenic,
            "PS + PM + >=4 PP".to_string(),
        ));
    }

    // ── Likely Pathogenic (7 combinations, includes ClinGen SVI PVS+PP) ──
    if pvs >= 1 && pm >= 1 {
        return Some((AcmgClassification::LikelyPathogenic, "PVS + PM".to_string()));
    }
    // ClinGen SVI Sept 2020: PVS + ≥1 PP → LP (compensates PM2 downgrade).
    // Bayesian Post_P = 0.988, within LP range (0.90–0.99).
    if pvs >= 1 && pp >= 1 {
        return Some((
            AcmgClassification::LikelyPathogenic,
            "PVS + >=1 PP (SVI)".to_string(),
        ));
    }
    if ps >= 1 && (1..=2).contains(&pm) {
        return Some((
            AcmgClassification::LikelyPathogenic,
            "PS + 1-2 PM".to_string(),
        ));
    }
    if ps >= 1 && pp >= 2 {
        return Some((
            AcmgClassification::LikelyPathogenic,
            "PS + >=2 PP".to_string(),
        ));
    }
    if pm >= 3 {
        return Some((AcmgClassification::LikelyPathogenic, ">=3 PM".to_string()));
    }
    if pm >= 2 && pp >= 2 {
        return Some((
            AcmgClassification::LikelyPathogenic,
            "2 PM + >=2 PP".to_string(),
        ));
    }
    if pm >= 1 && pp >= 4 {
        return Some((
            AcmgClassification::LikelyPathogenic,
            "PM + >=4 PP".to_string(),
        ));
    }

    None
}

/// Returns the strongest benign call (or None) along with the rule name.
fn compute_benign_call(ba: u8, bs: u8, bp: u8) -> Option<(AcmgClassification, String)> {
    if ba >= 1 {
        return Some((AcmgClassification::Benign, "BA1".to_string()));
    }
    if bs >= 2 {
        return Some((AcmgClassification::Benign, ">=2 BS".to_string()));
    }
    if bs >= 1 && bp >= 1 {
        return Some((AcmgClassification::LikelyBenign, "BS + BP".to_string()));
    }
    if bp >= 2 {
        return Some((AcmgClassification::LikelyBenign, ">=2 BP".to_string()));
    }
    None
}

/// A "definite" call is one that reaches P/LP or B/LB — anything strong
/// enough that mixing it with the opposite direction warrants a Conflicting
/// label rather than letting the stronger side win.
fn is_definite(cls: AcmgClassification) -> bool {
    matches!(
        cls,
        AcmgClassification::Pathogenic
            | AcmgClassification::LikelyPathogenic
            | AcmgClassification::Benign
            | AcmgClassification::LikelyBenign
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::{EvidenceDirection, EvidenceStrength};

    fn make_criterion(
        code: &str,
        direction: EvidenceDirection,
        strength: EvidenceStrength,
        met: bool,
    ) -> EvidenceCriterion {
        EvidenceCriterion {
            code: code.to_string(),
            direction,
            strength,
            default_strength: strength,
            met,
            evaluated: true,
            summary: String::new(),
            details: serde_json::Value::Null,
        }
    }

    fn met(code: &str, dir: EvidenceDirection, strength: EvidenceStrength) -> EvidenceCriterion {
        make_criterion(code, dir, strength, true)
    }

    fn not_met(
        code: &str,
        dir: EvidenceDirection,
        strength: EvidenceStrength,
    ) -> EvidenceCriterion {
        make_criterion(code, dir, strength, false)
    }

    use EvidenceDirection::*;
    use EvidenceStrength::*;

    // ── Benign Rules ──

    #[test]
    fn test_ba1_standalone_benign() {
        let criteria = vec![met("BA1", Benign, Standalone)];
        let (cls, rule) = combine(&criteria);
        assert_eq!(cls, AcmgClassification::Benign);
        assert_eq!(rule.unwrap(), "BA1");
    }

    #[test]
    fn test_two_bs_benign() {
        let criteria = vec![
            met("BS1", Benign, Strong),
            met("BS2", Benign, Strong),
        ];
        let (cls, rule) = combine(&criteria);
        assert_eq!(cls, AcmgClassification::Benign);
        assert_eq!(rule.unwrap(), ">=2 BS");
    }

    // ── Pathogenic Rules ──

    #[test]
    fn test_pvs_plus_ps_pathogenic() {
        let criteria = vec![
            met("PVS1", Pathogenic, VeryStrong),
            met("PS4", Pathogenic, Strong),
        ];
        let (cls, rule) = combine(&criteria);
        assert_eq!(cls, AcmgClassification::Pathogenic);
        assert_eq!(rule.unwrap(), "PVS + >=1 PS");
    }

    #[test]
    fn test_pvs_plus_2pm_pathogenic() {
        let criteria = vec![
            met("PVS1", Pathogenic, VeryStrong),
            met("PM2", Pathogenic, Moderate),
            met("PM4", Pathogenic, Moderate),
        ];
        let (cls, rule) = combine(&criteria);
        assert_eq!(cls, AcmgClassification::Pathogenic);
        assert_eq!(rule.unwrap(), "PVS + >=2 PM");
    }

    #[test]
    fn test_pvs_plus_pm_plus_pp_pathogenic() {
        let criteria = vec![
            met("PVS1", Pathogenic, VeryStrong),
            met("PM2", Pathogenic, Moderate),
            met("PP3", Pathogenic, Supporting),
        ];
        let (cls, rule) = combine(&criteria);
        assert_eq!(cls, AcmgClassification::Pathogenic);
        assert_eq!(rule.unwrap(), "PVS + PM + PP");
    }

    #[test]
    fn test_pvs_plus_2pp_pathogenic() {
        let criteria = vec![
            met("PVS1", Pathogenic, VeryStrong),
            met("PP2", Pathogenic, Supporting),
            met("PP3", Pathogenic, Supporting),
        ];
        let (cls, rule) = combine(&criteria);
        assert_eq!(cls, AcmgClassification::Pathogenic);
        assert_eq!(rule.unwrap(), "PVS + >=2 PP");
    }

    #[test]
    fn test_two_ps_pathogenic() {
        let criteria = vec![
            met("PS1", Pathogenic, Strong),
            met("PS4", Pathogenic, Strong),
        ];
        let (cls, rule) = combine(&criteria);
        assert_eq!(cls, AcmgClassification::Pathogenic);
        assert_eq!(rule.unwrap(), ">=2 PS");
    }

    #[test]
    fn test_ps_plus_3pm_pathogenic() {
        let criteria = vec![
            met("PS4", Pathogenic, Strong),
            met("PM1", Pathogenic, Moderate),
            met("PM2", Pathogenic, Moderate),
            met("PM4", Pathogenic, Moderate),
        ];
        let (cls, rule) = combine(&criteria);
        assert_eq!(cls, AcmgClassification::Pathogenic);
        assert_eq!(rule.unwrap(), "PS + >=3 PM");
    }

    // ── Likely Pathogenic Rules ──

    #[test]
    fn test_pvs_plus_pm_likely_pathogenic() {
        let criteria = vec![
            met("PVS1", Pathogenic, VeryStrong),
            met("PM2", Pathogenic, Moderate),
        ];
        let (cls, rule) = combine(&criteria);
        assert_eq!(cls, AcmgClassification::LikelyPathogenic);
        assert_eq!(rule.unwrap(), "PVS + PM");
    }

    #[test]
    fn test_ps_plus_pm_likely_pathogenic() {
        let criteria = vec![
            met("PS4", Pathogenic, Strong),
            met("PM2", Pathogenic, Moderate),
        ];
        let (cls, rule) = combine(&criteria);
        assert_eq!(cls, AcmgClassification::LikelyPathogenic);
        assert_eq!(rule.unwrap(), "PS + 1-2 PM");
    }

    #[test]
    fn test_ps_plus_2pp_likely_pathogenic() {
        let criteria = vec![
            met("PS4", Pathogenic, Strong),
            met("PP2", Pathogenic, Supporting),
            met("PP3", Pathogenic, Supporting),
        ];
        let (cls, rule) = combine(&criteria);
        assert_eq!(cls, AcmgClassification::LikelyPathogenic);
        assert_eq!(rule.unwrap(), "PS + >=2 PP");
    }

    #[test]
    fn test_3pm_likely_pathogenic() {
        let criteria = vec![
            met("PM1", Pathogenic, Moderate),
            met("PM2", Pathogenic, Moderate),
            met("PM4", Pathogenic, Moderate),
        ];
        let (cls, rule) = combine(&criteria);
        assert_eq!(cls, AcmgClassification::LikelyPathogenic);
        assert_eq!(rule.unwrap(), ">=3 PM");
    }

    #[test]
    fn test_2pm_2pp_likely_pathogenic() {
        let criteria = vec![
            met("PM2", Pathogenic, Moderate),
            met("PM4", Pathogenic, Moderate),
            met("PP2", Pathogenic, Supporting),
            met("PP3", Pathogenic, Supporting),
        ];
        let (cls, rule) = combine(&criteria);
        assert_eq!(cls, AcmgClassification::LikelyPathogenic);
        assert_eq!(rule.unwrap(), "2 PM + >=2 PP");
    }

    // ── ClinGen SVI PVS + PP Rule ──

    #[test]
    fn test_pvs_plus_pp_likely_pathogenic_svi() {
        // ClinGen SVI (Sept 2020): PVS + >=1 PP → LP
        // This is the key rule for PVS1 + PM2_Supporting
        let criteria = vec![
            met("PVS1", Pathogenic, VeryStrong),
            met("PM2_Supporting", Pathogenic, Supporting),
        ];
        let (cls, rule) = combine(&criteria);
        assert_eq!(cls, AcmgClassification::LikelyPathogenic);
        assert_eq!(rule.unwrap(), "PVS + >=1 PP (SVI)");
    }

    // ── Likely Benign Rules ──

    #[test]
    fn test_bs_plus_bp_likely_benign() {
        let criteria = vec![
            met("BS1", Benign, Strong),
            met("BP7", Benign, Supporting),
        ];
        let (cls, rule) = combine(&criteria);
        assert_eq!(cls, AcmgClassification::LikelyBenign);
        assert_eq!(rule.unwrap(), "BS + BP");
    }

    #[test]
    fn test_2bp_likely_benign() {
        let criteria = vec![
            met("BP4", Benign, Supporting),
            met("BP7", Benign, Supporting),
        ];
        let (cls, rule) = combine(&criteria);
        assert_eq!(cls, AcmgClassification::LikelyBenign);
        assert_eq!(rule.unwrap(), ">=2 BP");
    }

    // ── Conflicting Evidence ──

    #[test]
    fn test_conflicting_evidence_definite_both_directions() {
        // Pathogenic-direction definite (PVS+PM → LP) + Benign-direction
        // definite (≥2 BS → Benign) ⇒ Conflicting → VUS.
        let criteria = vec![
            met("PVS1", Pathogenic, VeryStrong),
            met("PM2", Pathogenic, Moderate),
            met("BS1", Benign, Strong),
            met("BS2", Benign, Strong),
        ];
        let (cls, rule) = combine(&criteria);
        assert_eq!(cls, AcmgClassification::UncertainSignificance);
        assert!(rule.unwrap().contains("Conflicting"));
    }

    #[test]
    fn test_pvs1_with_lone_bs_does_not_auto_conflict() {
        // PR9 fix: PVS1 alone + BS1 alone reach NEITHER directional call
        // (need PVS+PS/PM/PP for pathogenic, ≥2 BS for benign), so result
        // should be plain VUS, not "Conflicting". Pre-PR9 short-circuit
        // would have flagged it as Conflicting.
        let criteria = vec![
            met("PVS1", Pathogenic, VeryStrong),
            met("BS1", Benign, Strong),
        ];
        let (cls, rule) = combine(&criteria);
        assert_eq!(cls, AcmgClassification::UncertainSignificance);
        assert!(rule.is_none(), "expected no rule, got {:?}", rule);
    }

    #[test]
    fn test_pvs_plus_bp4_supporting_does_not_auto_conflict() {
        // PVS1 + PM2_Supporting + BP4_Supporting under PR9: pathogenic side
        // fires PVS+>=1 PP → LP; benign side has only 1 BP, sub-threshold for
        // any benign rule. Result should be LP, not auto-Conflicting.
        let criteria = vec![
            met("PVS1", Pathogenic, VeryStrong),
            met("PM2_Supporting", Pathogenic, Supporting),
            met("BP4", Benign, Supporting),
        ];
        let (cls, _) = combine(&criteria);
        assert_eq!(cls, AcmgClassification::LikelyPathogenic);
    }

    // ── Default VUS ──

    #[test]
    fn test_insufficient_evidence_vus() {
        let criteria = vec![
            met("PM2", Pathogenic, Supporting), // Note: Supporting due to SVI downgrade
        ];
        let (cls, rule) = combine(&criteria);
        assert_eq!(cls, AcmgClassification::UncertainSignificance);
        assert!(rule.is_none());
    }

    #[test]
    fn test_no_criteria_met_vus() {
        let criteria = vec![
            not_met("PVS1", Pathogenic, VeryStrong),
            not_met("BA1", Benign, Standalone),
        ];
        let (cls, rule) = combine(&criteria);
        assert_eq!(cls, AcmgClassification::UncertainSignificance);
        assert!(rule.is_none());
    }
}
