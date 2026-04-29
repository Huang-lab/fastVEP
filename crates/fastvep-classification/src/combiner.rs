use crate::types::{AcmgClassification, EvidenceCounts, EvidenceCriterion};

/// Apply ACMG-AMP combination rules to determine final classification.
///
/// Returns the classification and the name of the triggered rule.
///
/// Rules are applied in this order:
/// 1. Benign (BA1 standalone, or >=2 BS)
/// 2. Check for conflicting evidence (pathogenic + benign criteria both met -> VUS)
/// 3. Pathogenic (8 combinations)
/// 4. Likely Pathogenic (7 combinations, includes ClinGen SVI PVS+PP rule)
/// 5. Likely Benign (BS+BP, or >=2 BP)
/// 6. Default: VUS
///
/// Note: Includes the ClinGen SVI novel combination rule (Sept 2020) that
/// PVS + >=1 PP → Likely Pathogenic. This rule was added to compensate for
/// the PM2 downgrade to Supporting, ensuring that PVS1 + PM2_Supporting
/// still reaches LP (Bayesian Post_P = 0.988, within LP range 0.90–0.99).
pub fn combine(criteria: &[EvidenceCriterion]) -> (AcmgClassification, Option<String>) {
    let counts = EvidenceCounts::from_criteria(criteria);

    let pvs = counts.pathogenic_very_strong;
    let ps = counts.pathogenic_strong;
    let pm = counts.pathogenic_moderate;
    let pp = counts.pathogenic_supporting;
    let ba = counts.benign_standalone;
    let bs = counts.benign_strong;
    let bp = counts.benign_supporting;

    // ── Benign (BA1 is standalone) ──
    if ba >= 1 {
        return (AcmgClassification::Benign, Some("BA1".to_string()));
    }
    if bs >= 2 {
        return (AcmgClassification::Benign, Some(">=2 BS".to_string()));
    }

    // ── Conflicting evidence → VUS ──
    // If both pathogenic and benign criteria are met, classify as VUS
    if counts.has_pathogenic() && counts.has_benign() {
        return (
            AcmgClassification::UncertainSignificance,
            Some("Conflicting pathogenic and benign evidence".to_string()),
        );
    }

    // ── Pathogenic (8 combinations, checked in order of specificity) ──
    if pvs >= 1 && ps >= 1 {
        return (
            AcmgClassification::Pathogenic,
            Some("PVS + >=1 PS".to_string()),
        );
    }
    if pvs >= 1 && pm >= 2 {
        return (
            AcmgClassification::Pathogenic,
            Some("PVS + >=2 PM".to_string()),
        );
    }
    if pvs >= 1 && pm >= 1 && pp >= 1 {
        return (
            AcmgClassification::Pathogenic,
            Some("PVS + PM + PP".to_string()),
        );
    }
    if pvs >= 1 && pp >= 2 {
        return (
            AcmgClassification::Pathogenic,
            Some("PVS + >=2 PP".to_string()),
        );
    }
    if ps >= 2 {
        return (
            AcmgClassification::Pathogenic,
            Some(">=2 PS".to_string()),
        );
    }
    if ps >= 1 && pm >= 3 {
        return (
            AcmgClassification::Pathogenic,
            Some("PS + >=3 PM".to_string()),
        );
    }
    if ps >= 1 && pm >= 2 && pp >= 2 {
        return (
            AcmgClassification::Pathogenic,
            Some("PS + 2 PM + >=2 PP".to_string()),
        );
    }
    if ps >= 1 && pm >= 1 && pp >= 4 {
        return (
            AcmgClassification::Pathogenic,
            Some("PS + PM + >=4 PP".to_string()),
        );
    }

    // ── Likely Pathogenic (7 combinations) ──
    if pvs >= 1 && pm >= 1 {
        return (
            AcmgClassification::LikelyPathogenic,
            Some("PVS + PM".to_string()),
        );
    }
    // ClinGen SVI novel rule (Sept 2020): PVS + >=1 PP → LP
    // Bayesian Post_P = 0.988, within LP range (0.90–0.99).
    // Compensates for PM2 downgrade to Supporting, so PVS1 + PM2_Supporting → LP.
    if pvs >= 1 && pp >= 1 {
        return (
            AcmgClassification::LikelyPathogenic,
            Some("PVS + >=1 PP (SVI)".to_string()),
        );
    }
    if ps >= 1 && (1..=2).contains(&pm) {
        return (
            AcmgClassification::LikelyPathogenic,
            Some("PS + 1-2 PM".to_string()),
        );
    }
    if ps >= 1 && pp >= 2 {
        return (
            AcmgClassification::LikelyPathogenic,
            Some("PS + >=2 PP".to_string()),
        );
    }
    if pm >= 3 {
        return (
            AcmgClassification::LikelyPathogenic,
            Some(">=3 PM".to_string()),
        );
    }
    if pm >= 2 && pp >= 2 {
        return (
            AcmgClassification::LikelyPathogenic,
            Some("2 PM + >=2 PP".to_string()),
        );
    }
    if pm >= 1 && pp >= 4 {
        return (
            AcmgClassification::LikelyPathogenic,
            Some("PM + >=4 PP".to_string()),
        );
    }

    // ── Likely Benign ──
    if bs >= 1 && bp >= 1 {
        return (
            AcmgClassification::LikelyBenign,
            Some("BS + BP".to_string()),
        );
    }
    if bp >= 2 {
        return (
            AcmgClassification::LikelyBenign,
            Some(">=2 BP".to_string()),
        );
    }

    // ── Default: Uncertain Significance ──
    (AcmgClassification::UncertainSignificance, None)
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
    fn test_conflicting_evidence_vus() {
        let criteria = vec![
            met("PVS1", Pathogenic, VeryStrong),
            met("BS1", Benign, Strong),
        ];
        let (cls, rule) = combine(&criteria);
        assert_eq!(cls, AcmgClassification::UncertainSignificance);
        assert!(rule.unwrap().contains("Conflicting"));
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
