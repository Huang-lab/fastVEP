use crate::types::{
    AcmgClassification, EvidenceCounts, EvidenceCriterion, EvidenceDirection, EvidenceStrength,
};

/// Apply ACMG-AMP combination rules to determine final classification.
///
/// Returns the classification and the name of the triggered rule.
///
/// Rules: pathogenic (8 combinations) and benign (BA1, `>=2` BS, BS+BP,
/// `>=2` BP) are evaluated independently. If both directions reach a
/// definite call (P/LP and B/LB), the result is VUS (Conflicting).
/// Otherwise the directional call wins.
///
/// Includes the ClinGen SVI novel combination rule (Sept 2020) that
/// PVS + >=1 PP → Likely Pathogenic, added to compensate for the PM2
/// downgrade to Supporting (Bayesian Post_P = 0.988, within LP range).
pub fn combine(
    criteria: &[EvidenceCriterion],
    config: &crate::config::AcmgConfig,
) -> (AcmgClassification, Option<String>) {
    combine_inner(criteria, config.use_point_system)
}

/// The ClinGen SVI point system (Tavtigian et al. 2020, Genet Med 22:1735).
///
/// The 2015 combining table and this point system are two encodings of the same
/// Bayesian model, and where they disagree the SVI's own analysis is that the
/// table has gaps. The clearest is a lone PVS1: 8 points, squarely inside the
/// Likely Pathogenic band, but no row of Table 5 matches it, so the table gives
/// Uncertain Significance. Run v10 had **2,319 truth-pathogenic variants land in
/// VUS on `PVS1` and nothing else** - a nonsense or frameshift variant in a
/// haploinsufficient disease gene, which is close to the least uncertain thing a
/// classifier sees.
///
/// Points are `1 / 2 / 4 / 8` for Supporting / Moderate / Strong / Very Strong,
/// negated for benign criteria, with a standalone benign criterion at `-8`.
/// Bands: `>= 10` Pathogenic, `6..=9` Likely Pathogenic, `-6..=-2` Likely
/// Benign, `<= -7` Benign, everything else Uncertain.
///
/// The Likely Benign floor is the one deliberate deviation from Tavtigian's
/// published table, which puts it at `-1`. At `-1` a single BP4 - one
/// supporting benign criterion, on its own - is enough for a Likely Benign
/// call, where Richards 2015 requires two benign supporting criteria or a
/// strong plus a supporting. Measured on the benchmark sample, the `-1` floor
/// produced 36 false-benign calls against 1, and **22 of those 36 were a lone
/// BP4**. Taking the stricter of the two schemes on the benign side costs
/// almost nothing and removes the whole class: a missed diagnosis and a variant
/// left uncertain are not comparable errors.
///
/// Conflict handling falls out of the arithmetic rather than needing a separate
/// rule: evidence in both directions cancels toward zero, which is the VUS band.
/// Pathogenic-side call from the ClinGen SVI point system, or `None` when the
/// pathogenic evidence does not reach Likely Pathogenic.
///
/// Points are `1 / 2 / 4 / 8` for Supporting / Moderate / Strong / Very Strong.
/// Bands: `>= 10` Pathogenic, `6..=9` Likely Pathogenic.
///
/// **Only the pathogenic direction is scored this way**, and the benign
/// direction keeps the Richards 2015 rules. That split is deliberate and was
/// measured, not assumed. Tavtigian's Likely Benign band opens at `-1`, so a
/// single BP4 - one supporting benign criterion on its own - is enough for a
/// Likely Benign call, where Richards requires two benign supporting criteria
/// or a strong plus a supporting. On the benchmark sample the full point system
/// scored 36 false-benign calls against the table's 1, and **22 of the 36 were
/// a lone BP4**. Tightening the band instead (floor `-2`) removed those but
/// took benign recall from 56.3 % down to 45.6 %, because a lone BP4 carries a
/// great many correct benign calls too.
///
/// So the two schemes are used where each is the safer one. The point system's
/// documented gap is on the pathogenic side: a lone PVS1 is 8 points, inside
/// Likely Pathogenic, but matches no row of Table 5 and so returns Uncertain -
/// which put 2,319 truth-pathogenic variants in VUS in run v10. Its divergence
/// on the benign side runs the other way, toward calling variants benign on
/// thinner evidence, and a missed diagnosis is not the same kind of error as a
/// variant left uncertain.
fn pathogenic_call_by_points(
    criteria: &[EvidenceCriterion],
) -> Option<(AcmgClassification, String)> {
    let points: i32 = criteria
        .iter()
        .filter(|c| c.met && c.direction == EvidenceDirection::Pathogenic)
        .map(|c| match c.strength {
            EvidenceStrength::Supporting => 1,
            EvidenceStrength::Moderate => 2,
            EvidenceStrength::Strong => 4,
            EvidenceStrength::VeryStrong | EvidenceStrength::Standalone => 8,
        })
        .sum();

    let cls = if points >= 10 {
        AcmgClassification::Pathogenic
    } else if points >= 6 {
        AcmgClassification::LikelyPathogenic
    } else {
        return None;
    };
    Some((cls, format!("{} pathogenic points (Tavtigian 2020)", points)))
}

/// The Richards 2015 combining table on both sides. The tests use it directly
/// to pin table behaviour independently of the configured default.
#[cfg(test)]
fn combine_by_table(criteria: &[EvidenceCriterion]) -> (AcmgClassification, Option<String>) {
    combine_inner(criteria, false)
}

fn combine_inner(
    criteria: &[EvidenceCriterion],
    use_points: bool,
) -> (AcmgClassification, Option<String>) {
    let counts = EvidenceCounts::from_criteria(criteria);
    let pvs = counts.pathogenic_very_strong;
    let ps = counts.pathogenic_strong;
    let pm = counts.pathogenic_moderate;
    let pp = counts.pathogenic_supporting;
    let ba = counts.benign_standalone;
    let bs = counts.benign_strong;
    let bp = counts.benign_supporting;

    let pathogenic_call = if use_points {
        pathogenic_call_by_points(criteria)
    } else {
        compute_pathogenic_call(pvs, ps, pm, pp)
    };
    let benign_call = compute_benign_call(ba, bs, bp);

    // Richards 2015: "If a variant does not fulfil criteria using either of
    // these sets, or the evidence for benign and pathogenic is conflicting,
    // the variant defaults to uncertain significance." The rule below catches
    // the case where one direction reaches a definite call while the other
    // holds evidence at Strong or above without itself reaching a call.
    //
    // That asymmetry produced most of the false-benign calls in the round-2
    // medical-genetics review: a met PS1 (Strong) plus two benign supporting
    // criteria was reported Likely Benign, with the Strong pathogenic evidence
    // silently discarded. 36 of the 78 such rows carried a met PS1 or PVS1.
    let strong_pathogenic_present = pvs >= 1 || ps >= 1;
    let strong_benign_present = ba >= 1 || bs >= 1;

    if let Some((cls, rule)) = &benign_call {
        if is_definite(*cls) && strong_pathogenic_present && pathogenic_call.is_none() {
            return (
                AcmgClassification::UncertainSignificance,
                Some(format!(
                    "Conflicting evidence: benign rules → {} ({}), but {} pathogenic criteri{} at Strong or above also met",
                    cls.shorthand(),
                    rule,
                    pvs + ps,
                    if pvs + ps == 1 { "on is" } else { "a are" }
                )),
            );
        }
    }
    if let Some((cls, rule)) = &pathogenic_call {
        if is_definite(*cls) && strong_benign_present && benign_call.is_none() {
            return (
                AcmgClassification::UncertainSignificance,
                Some(format!(
                    "Conflicting evidence: pathogenic rules → {} ({}), but {} benign criteri{} at Strong or above also met",
                    cls.shorthand(),
                    rule,
                    ba + bs,
                    if ba + bs == 1 { "on is" } else { "a are" }
                )),
            );
        }
    }

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
        let (cls, rule) = combine_by_table(&criteria);
        assert_eq!(cls, AcmgClassification::Benign);
        assert_eq!(rule.unwrap(), "BA1");
    }

    #[test]
    fn test_two_bs_benign() {
        let criteria = vec![
            met("BS1", Benign, Strong),
            met("BS2", Benign, Strong),
        ];
        let (cls, rule) = combine_by_table(&criteria);
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
        let (cls, rule) = combine_by_table(&criteria);
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
        let (cls, rule) = combine_by_table(&criteria);
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
        let (cls, rule) = combine_by_table(&criteria);
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
        let (cls, rule) = combine_by_table(&criteria);
        assert_eq!(cls, AcmgClassification::Pathogenic);
        assert_eq!(rule.unwrap(), "PVS + >=2 PP");
    }

    #[test]
    fn test_two_ps_pathogenic() {
        let criteria = vec![
            met("PS1", Pathogenic, Strong),
            met("PS4", Pathogenic, Strong),
        ];
        let (cls, rule) = combine_by_table(&criteria);
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
        let (cls, rule) = combine_by_table(&criteria);
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
        let (cls, rule) = combine_by_table(&criteria);
        assert_eq!(cls, AcmgClassification::LikelyPathogenic);
        assert_eq!(rule.unwrap(), "PVS + PM");
    }

    #[test]
    fn test_ps_plus_pm_likely_pathogenic() {
        let criteria = vec![
            met("PS4", Pathogenic, Strong),
            met("PM2", Pathogenic, Moderate),
        ];
        let (cls, rule) = combine_by_table(&criteria);
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
        let (cls, rule) = combine_by_table(&criteria);
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
        let (cls, rule) = combine_by_table(&criteria);
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
        let (cls, rule) = combine_by_table(&criteria);
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
        let (cls, rule) = combine_by_table(&criteria);
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
        let (cls, rule) = combine_by_table(&criteria);
        assert_eq!(cls, AcmgClassification::LikelyBenign);
        assert_eq!(rule.unwrap(), "BS + BP");
    }

    #[test]
    fn test_2bp_likely_benign() {
        let criteria = vec![
            met("BP4", Benign, Supporting),
            met("BP7", Benign, Supporting),
        ];
        let (cls, rule) = combine_by_table(&criteria);
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
        let (cls, rule) = combine_by_table(&criteria);
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
        let (cls, rule) = combine_by_table(&criteria);
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
        let (cls, _) = combine_by_table(&criteria);
        assert_eq!(cls, AcmgClassification::LikelyPathogenic);
    }

    // ── Default VUS ──

    #[test]
    fn test_insufficient_evidence_vus() {
        let criteria = vec![
            met("PM2", Pathogenic, Supporting), // Note: Supporting due to SVI downgrade
        ];
        let (cls, rule) = combine_by_table(&criteria);
        assert_eq!(cls, AcmgClassification::UncertainSignificance);
        assert!(rule.is_none());
    }

    #[test]
    fn test_no_criteria_met_vus() {
        let criteria = vec![
            not_met("PVS1", Pathogenic, VeryStrong),
            not_met("BA1", Benign, Standalone),
        ];
        let (cls, rule) = combine_by_table(&criteria);
        assert_eq!(cls, AcmgClassification::UncertainSignificance);
        assert!(rule.is_none());
    }

    // ── Asymmetric conflict guard (Richards 2015) ──

    #[test]
    fn test_strong_pathogenic_blocks_likely_benign_call() {
        // PS1 (Strong) alongside two BP criteria used to be reported Likely
        // Benign with the Strong pathogenic evidence silently discarded. 36 of
        // the 78 false-benign rows in the round-2 medical-genetics review had
        // exactly this shape.
        let criteria = vec![
            met("PS1", Pathogenic, Strong),
            met("BP1", Benign, Supporting),
            met("BP4", Benign, Supporting),
        ];
        let (cls, rule) = combine_by_table(&criteria);
        assert_eq!(cls, AcmgClassification::UncertainSignificance);
        assert!(rule.unwrap().contains("Conflicting evidence"));
    }

    #[test]
    fn test_very_strong_pathogenic_blocks_benign_call() {
        let criteria = vec![
            met("PVS1", Pathogenic, VeryStrong),
            met("BS1", Benign, Strong),
            met("BS2", Benign, Strong),
        ];
        let (cls, _) = combine_by_table(&criteria);
        assert_eq!(cls, AcmgClassification::UncertainSignificance);
    }

    #[test]
    fn test_strong_benign_blocks_likely_pathogenic_call() {
        // Mirror direction: a met BS1 alongside PM+PP evidence that reaches LP.
        let criteria = vec![
            met("PM1", Pathogenic, Moderate),
            met("PM2", Pathogenic, Moderate),
            met("PP2", Pathogenic, Supporting),
            met("PP3", Pathogenic, Supporting),
            met("BS1", Benign, Strong),
        ];
        let (cls, rule) = combine_by_table(&criteria);
        assert_eq!(cls, AcmgClassification::UncertainSignificance);
        assert!(rule.unwrap().contains("Conflicting evidence"));
    }

    #[test]
    fn test_supporting_benign_alone_does_not_block_pathogenic_call() {
        // The guard is deliberately limited to Strong-or-above opposing
        // evidence. A lone BP against a PVS1 is ordinary noise, not a conflict,
        // and must not turn every pathogenic call into VUS.
        let criteria = vec![
            met("PVS1", Pathogenic, VeryStrong),
            met("PM2", Pathogenic, Moderate),
            met("BP3", Benign, Supporting),
        ];
        let (cls, _) = combine_by_table(&criteria);
        assert_eq!(cls, AcmgClassification::LikelyPathogenic);
    }

    // ── ClinGen SVI point system (pathogenic side) ───────────────────────

    fn points_of(criteria: &[EvidenceCriterion]) -> (AcmgClassification, Option<String>) {
        combine_inner(criteria, true)
    }

    #[test]
    fn lone_pvs1_is_likely_pathogenic_under_points_and_vus_under_the_table() {
        // The disagreement that motivated the change. 2,319 truth-pathogenic
        // variants sat in VUS on exactly this signature in run v10.
        let criteria = vec![met("PVS1", EvidenceDirection::Pathogenic, EvidenceStrength::VeryStrong)];
        assert_eq!(points_of(&criteria).0, AcmgClassification::LikelyPathogenic);
        assert_eq!(
            combine_by_table(&criteria).0,
            AcmgClassification::UncertainSignificance
        );
    }

    #[test]
    fn ten_points_reaches_pathogenic() {
        // PVS1 (8) + PM2_Supporting (1) = 9 → LP; add a Moderate → 11 → P.
        let mut criteria = vec![
            met("PVS1", EvidenceDirection::Pathogenic, EvidenceStrength::VeryStrong),
            met("PM2_Supporting", EvidenceDirection::Pathogenic, EvidenceStrength::Supporting),
        ];
        assert_eq!(points_of(&criteria).0, AcmgClassification::LikelyPathogenic);
        criteria.push(met("PM1", EvidenceDirection::Pathogenic, EvidenceStrength::Moderate));
        assert_eq!(points_of(&criteria).0, AcmgClassification::Pathogenic);
    }

    #[test]
    fn five_pathogenic_points_stay_uncertain() {
        // PM1 (2) + PM5 (2) + PM2_Supporting (1) = 5, one short of the band.
        let criteria = vec![
            met("PM1", EvidenceDirection::Pathogenic, EvidenceStrength::Moderate),
            met("PM5", EvidenceDirection::Pathogenic, EvidenceStrength::Moderate),
            met("PM2_Supporting", EvidenceDirection::Pathogenic, EvidenceStrength::Supporting),
        ];
        assert_eq!(points_of(&criteria).0, AcmgClassification::UncertainSignificance);
    }

    #[test]
    fn the_benign_side_keeps_the_2015_table_under_points() {
        // The measured reason for the split: Tavtigian's Likely Benign band
        // opens at -1, so a lone BP4 would be a benign call. It must not be.
        let criteria = vec![met("BP4", EvidenceDirection::Benign, EvidenceStrength::Supporting)];
        assert_eq!(points_of(&criteria).0, AcmgClassification::UncertainSignificance);
        // Two supporting benign criteria do reach Likely Benign, as Richards says.
        let criteria = vec![
            met("BP4", EvidenceDirection::Benign, EvidenceStrength::Supporting),
            met("BP7", EvidenceDirection::Benign, EvidenceStrength::Supporting),
        ];
        assert_eq!(points_of(&criteria).0, AcmgClassification::LikelyBenign);
    }

    #[test]
    fn point_scoring_does_not_disable_the_conflict_guard() {
        // A lone PVS1 now reaches LP on the pathogenic side, so the guard that
        // sends definite-vs-definite to VUS has to keep working over it.
        let criteria = vec![
            met("PVS1", EvidenceDirection::Pathogenic, EvidenceStrength::VeryStrong),
            met("BA1", EvidenceDirection::Benign, EvidenceStrength::Standalone),
        ];
        let (cls, rule) = points_of(&criteria);
        assert_eq!(cls, AcmgClassification::UncertainSignificance);
        assert!(rule.unwrap_or_default().contains("Conflicting"));
    }

    #[test]
    fn graded_subcodes_score_at_their_effective_strength() {
        // PVS1_Moderate is 2 points, not 8 — the grading has to reach the sum.
        let criteria = vec![
            met("PVS1_Moderate", EvidenceDirection::Pathogenic, EvidenceStrength::Moderate),
            met("PM2_Supporting", EvidenceDirection::Pathogenic, EvidenceStrength::Supporting),
        ];
        assert_eq!(points_of(&criteria).0, AcmgClassification::UncertainSignificance);
    }
}
