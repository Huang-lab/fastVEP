use fastvep_core::Consequence;

use crate::config::AcmgConfig;
use crate::sa_extract::ClassificationInput;
use crate::types::{EvidenceCriterion, EvidenceDirection, EvidenceStrength};

/// Evaluate all pathogenic moderate criteria: PM1, PM2, PM3, PM4, PM5, PM6.
pub fn evaluate_all(
    input: &ClassificationInput,
    config: &AcmgConfig,
) -> Vec<EvidenceCriterion> {
    vec![
        evaluate_pm1(input, config),
        evaluate_pm2(input, config),
        evaluate_pm3(input, config),
        evaluate_pm4(input, config),
        evaluate_pm5(input, config),
        evaluate_pm6(input, config),
    ]
}

/// PM1: Located in a mutational hot spot and/or critical functional domain.
///
/// Approximated using ClinVar pathogenic variant density as a hotspot proxy:
/// if >=N pathogenic variants exist within ±W amino acid positions, the region
/// is considered a hotspot.
fn evaluate_pm1(
    input: &ClassificationInput,
    config: &AcmgConfig,
) -> EvidenceCriterion {
    let mut details = serde_json::Map::new();
    let window = config.pm1_hotspot_window;
    let threshold = config.pm1_hotspot_min_pathogenic;
    details.insert("hotspot_window".into(), serde_json::json!(window));
    details.insert("hotspot_threshold".into(), serde_json::json!(threshold));

    // PM1 is residue-level evidence: a hotspot or critical domain argues that
    // *this* amino-acid substitution matters. It is defined for missense and
    // in-frame changes. Applied to a frameshift or nonsense variant it adds
    // nothing PVS1 does not already carry, and the round-2 review flagged
    // exactly that stacking on CBS, MSH6 and RYR1 ("PM1 is called with PVS1?
    // Where is the evidence for PM1 coming?"). The PVS1 co-occurrence case is
    // additionally suppressed in the reconciliation pass.
    let pm1_eligible_consequence = input.consequences.iter().any(|c| {
        matches!(
            c,
            Consequence::MissenseVariant
                | Consequence::InframeInsertion
                | Consequence::InframeDeletion
                | Consequence::ProteinAlteringVariant
                | Consequence::StopLost
                | Consequence::StartLost
        )
    });
    if !pm1_eligible_consequence {
        return EvidenceCriterion {
            code: "PM1".to_string(),
            direction: EvidenceDirection::Pathogenic,
            strength: EvidenceStrength::Moderate,
            default_strength: EvidenceStrength::Moderate,
            met: false,
            evaluated: true,
            summary: "PM1 applies to missense and in-frame changes; residue-level hotspot evidence does not apply to this consequence".to_string(),
            details: serde_json::Value::Object(details),
        };
    }

    // PM1 reads a cluster of pathogenic variants as marking a critical region.
    // Where the gene has no established disease relationship there is no
    // disease for that cluster to be critical to, and the ClinVar entries
    // behind it are exactly the assertions the validity curation declined to
    // accept. Checked after the consequence test so that a frameshift still
    // gets the more specific "wrong consequence for PM1" answer.
    if let Some(reason) = super::gene_disease::validity_blocker(input, config) {
        return EvidenceCriterion {
            code: "PM1".to_string(),
            direction: EvidenceDirection::Pathogenic,
            strength: EvidenceStrength::Moderate,
            default_strength: EvidenceStrength::Moderate,
            met: false,
            evaluated: false,
            summary: reason,
            details: serde_json::Value::Object(details),
        };
    }

    let prot_pos = match input.protein_position {
        Some(pos) => pos,
        None => {
            return EvidenceCriterion {
                code: "PM1".to_string(),
                direction: EvidenceDirection::Pathogenic,
                strength: EvidenceStrength::Moderate,
                default_strength: EvidenceStrength::Moderate,
                met: false,
                evaluated: false,
                summary: "Protein position not available".to_string(),
                details: serde_json::Value::Object(details),
            };
        }
    };

    details.insert("protein_position".into(), serde_json::json!(prot_pos));

    if let Some(ref cpd) = input.clinvar_protein {
        let low = prot_pos.saturating_sub(window);
        let high = prot_pos + window;
        // Substring match `"pathogenic"` is intentional: it counts both
        // ClinVar `Pathogenic` and `Likely_pathogenic` records as
        // hotspot evidence. Most VCEPs (and ClinGen PM1 calibration)
        // include LP variants in residue-density counts because LP
        // variants are still ≥90 % posterior probability of
        // pathogenicity per the Tavtigian Bayesian framework — the
        // hotspot signal is robust to LP/P aggregation.
        let raw_nearby: usize = cpd
            .protein_variants
            .iter()
            .filter(|v| v.pos >= low && v.pos <= high && v.sig.to_lowercase().contains("pathogenic"))
            .count();

        // The index is built from ClinVar, so when the variant being
        // classified is itself ClinVar pathogenic its own record is one of the
        // neighbours counted here. Counting it makes PM1 partly self-derived
        // (and inflates any ClinVar-based benchmark), so discount it.
        let self_contribution = if config.exclude_self_from_clinvar_evidence
            && input.clinvar.as_ref().is_some_and(|c| c.has_pathogenic())
        {
            1
        } else {
            0
        };
        let nearby_pathogenic = raw_nearby.saturating_sub(self_contribution);

        details.insert("nearby_pathogenic_count".into(), serde_json::json!(nearby_pathogenic));
        if self_contribution > 0 {
            details.insert("self_excluded_from_count".into(), serde_json::json!(true));
            details.insert("nearby_pathogenic_raw".into(), serde_json::json!(raw_nearby));
        }

        // The other half of PM1's definition. Richards 2015 asks for a hotspot
        // or critical domain "**without benign variation**", and fastVEP counted
        // only the pathogenic half for as long as the index carried only that
        // half. A window that also holds benign missense is not a region where
        // any substitution is damaging - it is a region ClinVar has looked at
        // from both sides.
        //
        // Measured on the reviewer's rows: this resolves MSH2 p.Gly315Val (3
        // pathogenic and 5 benign neighbours; the MSH2 VCEP specification
        // cspec GN137 excludes PM1 for the gene outright) and leaves TP53
        // p.Arg248 alone (23 pathogenic, 0 benign). It does *not* resolve her
        // CHD7 and PTCH1 rows, where the objection was that the gene has no
        // hotspots at all rather than that this window has benign variation;
        // both still show 3 pathogenic and 0 benign neighbours here.
        //
        // Only applied when the index actually carries benign assertions.
        // Against an older `.oga` the count is structurally zero and testing it
        // would silently pass every window.
        //
        // Leave-one-out applies here for the same reason it applies to the
        // pathogenic count: the index is ClinVar, so a variant ClinVar calls
        // benign is one of its own benign neighbours, and letting it veto PM1
        // is that variant's ClinVar label deciding its own classification.
        let benign_tested = cpd.benign_indexed;
        let nearby_benign: usize = if benign_tested {
            let raw = cpd
                .protein_variants
                .iter()
                .filter(|v| v.pos >= low && v.pos <= high && v.sig.to_lowercase().contains("benign"))
                .count();
            let self_benign = if config.exclude_self_from_clinvar_evidence
                && input.clinvar.as_ref().is_some_and(|c| c.has_benign())
            {
                1
            } else {
                0
            };
            raw.saturating_sub(self_benign)
        } else {
            0
        };
        if benign_tested {
            details.insert("nearby_benign_count".into(), serde_json::json!(nearby_benign));
            details.insert(
                "max_benign_in_window".into(),
                serde_json::json!(config.pm1_max_benign_in_window),
            );
        }
        let benign_variation = benign_tested
            && nearby_benign > config.pm1_max_benign_in_window as usize;

        let met = nearby_pathogenic >= threshold as usize && !benign_variation;
        let summary = if benign_variation {
            format!(
                "Not a hotspot: {} pathogenic variants within ±{} AA of position {}, but the same window carries {} benign/likely-benign missense variants in ClinVar; PM1 requires a region without benign variation (Richards 2015)",
                nearby_pathogenic, window, prot_pos, nearby_benign
            )
        } else if met {
            format!(
                "Mutational hotspot: {} pathogenic variants within ±{} AA of position {} (threshold: {})",
                nearby_pathogenic, window, prot_pos, threshold
            )
        } else {
            format!(
                "Not a hotspot: {} pathogenic variants within ±{} AA of position {} (threshold: {})",
                nearby_pathogenic, window, prot_pos, threshold
            )
        };

        EvidenceCriterion {
            code: "PM1".to_string(),
            direction: EvidenceDirection::Pathogenic,
            strength: EvidenceStrength::Moderate,
            default_strength: EvidenceStrength::Moderate,
            met,
            evaluated: true,
            summary,
            details: serde_json::Value::Object(details),
        }
    } else {
        EvidenceCriterion {
            code: "PM1".to_string(),
            direction: EvidenceDirection::Pathogenic,
            strength: EvidenceStrength::Moderate,
            default_strength: EvidenceStrength::Moderate,
            met: false,
            evaluated: false,
            summary: "ClinVar protein-position index not available for hotspot analysis".to_string(),
            details: serde_json::Value::Object(details),
        }
    }
}

/// PM2: Absent from controls (or at extremely low frequency if recessive).
///
/// Per ClinGen SVI v1.0 (Sept 2020):
/// - Strength is Supporting (downgraded from Moderate by default).
/// - Use raw gnomAD allele frequency (NOT filtering allele frequency / FAF).
/// - Threshold depends on inheritance:
///     * AD / unknown: strict absence (AC = 0 or AF = 0).
///     * AR: AF ≤ 0.00007 (0.007%).
///
/// Inheritance is inferred from OMIM phenotypes (`OmimData::has_recessive_inheritance` /
/// `has_dominant_inheritance`). When a per-gene `pm2_af_threshold` override is
/// configured, that value wins regardless of inheritance.
fn evaluate_pm2(
    input: &ClassificationInput,
    config: &AcmgConfig,
) -> EvidenceCriterion {
    // "Absent from population databases" is only evidence when the database
    // could have seen the variant. Where reads pile up on a paralogue, or
    // gnomAD itself rejected the site, both a frequency and the absence of a
    // frequency are artefacts. This is the same gate BA1/BS1/BS2 apply, so a
    // site is never trusted for one frequency criterion and distrusted for
    // another - only the benign-direction preconditions differ.
    if let Some(blocker) = super::frequency_gate::data_blocker(input, config) {
        let mut details = serde_json::Map::new();
        blocker.record(&mut details);
        return EvidenceCriterion {
            code: if config.pm2_downgrade_to_supporting { "PM2_Supporting".to_string() } else { "PM2".to_string() },
            direction: EvidenceDirection::Pathogenic,
            strength: if config.pm2_downgrade_to_supporting { EvidenceStrength::Supporting } else { EvidenceStrength::Moderate },
            default_strength: EvidenceStrength::Moderate,
            met: false,
            evaluated: false,
            summary: format!("PM2 not evaluated: {}", blocker.reason),
            details: serde_json::Value::Object(details),
        };
    }

    let strength = if config.pm2_downgrade_to_supporting {
        EvidenceStrength::Supporting
    } else {
        EvidenceStrength::Moderate
    };
    let code = if config.pm2_downgrade_to_supporting {
        "PM2_Supporting".to_string()
    } else {
        "PM2".to_string()
    };

    // Determine the effective threshold and which inheritance rule applied.
    // Order of precedence:
    //   1. Per-gene override (config.gene_overrides[GENE].pm2_af_threshold)
    //   2. AR inheritance (from OMIM) → config.pm2_ar_af_threshold
    //   3. AD or unknown → config.pm2_ad_af_threshold (default 0.0 = strict absence)
    //
    // The legacy single-threshold field `config.pm2_af_threshold` is kept as a
    // fallback when neither AD nor AR-specific knobs are configured (back-compat).
    let gene = input.gene_symbol.as_deref();
    let gene_specific_threshold = gene.and_then(|g| {
        config
            .gene_overrides
            .get(g)
            .and_then(|o| o.pm2_af_threshold)
    });

    let is_recessive = input
        .omim
        .as_ref()
        .is_some_and(|o| o.has_recessive_inheritance());
    let is_dominant = input
        .omim
        .as_ref()
        .is_some_and(|o| o.has_dominant_inheritance());

    let (threshold, inheritance_basis): (f64, &'static str) = if let Some(t) = gene_specific_threshold {
        (t, "gene_override")
    } else if is_recessive && !is_dominant {
        (config.pm2_ar_af_threshold, "AR")
    } else {
        (config.pm2_ad_af_threshold, "AD_or_unknown")
    };

    let mut details = serde_json::Map::new();
    details.insert("af_threshold".into(), serde_json::json!(threshold));
    details.insert("inheritance_basis".into(), serde_json::json!(inheritance_basis));
    details.insert("is_recessive".into(), serde_json::json!(is_recessive));
    details.insert("is_dominant".into(), serde_json::json!(is_dominant));

    let (met, evaluated, summary) = if let Some(ref gnomad) = input.gnomad {
        details.insert("gnomad_allAf".into(), serde_json::json!(gnomad.all_af));
        details.insert("gnomad_allAc".into(), serde_json::json!(gnomad.all_ac));

        // For strict absence (threshold = 0.0), require AC and AF to both be
        // PRESENT and equal to zero — treating missing fields as zero would
        // call PM2 on incomplete gnomAD records. For non-zero thresholds
        // (e.g. AR 0.00007), require AF present and ≤ threshold.
        if threshold == 0.0 {
            match (gnomad.all_ac, gnomad.all_af) {
                (Some(0), Some(0.0)) => (
                    true,
                    true,
                    format!(
                        "Absent in gnomAD (AC=0, AF=0, inheritance={})",
                        inheritance_basis
                    ),
                ),
                (Some(ac), Some(af)) => (
                    false,
                    true,
                    format!(
                        "Not absent in gnomAD (AF={:.6}, AC={}, inheritance={})",
                        af, ac, inheritance_basis
                    ),
                ),
                _ => (
                    false,
                    false,
                    format!(
                        "PM2 not evaluated: gnomAD record present but AC/AF missing (inheritance={})",
                        inheritance_basis
                    ),
                ),
            }
        } else {
            match gnomad.all_af {
                Some(af) if af <= threshold => (
                    true,
                    true,
                    format!(
                        "Rare in gnomAD (AF={:.6}, threshold={:.6}, inheritance={})",
                        af, threshold, inheritance_basis
                    ),
                ),
                Some(af) => (
                    false,
                    true,
                    format!(
                        "Not rare enough in gnomAD (AF={:.6}, threshold={:.6}, inheritance={})",
                        af, threshold, inheritance_basis
                    ),
                ),
                None => (
                    false,
                    false,
                    format!(
                        "PM2 not evaluated: gnomAD record present but AF missing (inheritance={})",
                        inheritance_basis
                    ),
                ),
            }
        }
    } else if config.pm2_absent_when_no_record {
        // No gnomAD record at the variant. ClinGen SVI v1.0: "absent or
        // extremely rare in population databases". The natural reading of
        // a missing record is that gnomAD has never observed the variant,
        // i.e. it IS absent. Fire PM2 (downgraded to Supporting per SVI).
        //
        // Configurable via `pm2_absent_when_no_record` (default true).
        // Set false for partial-coverage runs where gnomAD .osa covers
        // only some input regions and you want PM2 silenced outside that
        // coverage. Note that when `--sa-dir` includes no gnomAD .osa at
        // all, `input.gnomad` is None for every variant — disable this
        // flag in that case (or load gnomAD) to avoid firing PM2 globally.
        details.insert("gnomad_allAf".into(), serde_json::Value::Null);
        details.insert("pm2_absent_when_no_record".into(), serde_json::json!(true));

        // Frequency backstop: even without a gnomAD record, ClinVar ships
        // ExAC / 1000G / ESP allele frequencies. When those show the variant
        // is not extremely rare, the "absent" assumption is wrong — do not
        // fire PM2. This catches common variants that silently fail to match
        // gnomAD (e.g. indels: GIGYF2, NOTCH2; or coverage gaps: ASPM,
        // TOR1AIP1) in the ClinVar 2★+ benchmark. The bar is the PM2 rarity
        // threshold, floored at the AR "extremely rare" cutoff so the AD
        // strict-absence case (threshold 0.0) doesn't reject on noise.
        let crosscheck_cutoff = threshold.max(config.pm2_ar_af_threshold);
        let clinvar_pop_af = input.clinvar.as_ref().and_then(|c| c.max_pop_af());
        if let Some(af) = clinvar_pop_af {
            details.insert("clinvar_max_pop_af".into(), serde_json::json!(af));
        }
        match clinvar_pop_af {
            Some(af) if af > crosscheck_cutoff => (
                false,
                true,
                format!(
                    "No gnomAD record, but ClinVar population AF={:.6} exceeds the PM2 rarity bar ({:.6}); not treated as absent (inheritance={})",
                    af, crosscheck_cutoff, inheritance_basis
                ),
            ),
            _ => (
                true,
                true,
                format!(
                    "Absent in gnomAD: no record at this position/allele (inheritance={}, treating as absent per pm2_absent_when_no_record)",
                    inheritance_basis
                ),
            ),
        }
    } else {
        // Strict-coverage stance: cannot distinguish "absent from gnomAD"
        // from "gnomAD .osa not loaded for this region". Mark NotEvaluated.
        details.insert("gnomad_allAf".into(), serde_json::Value::Null);
        (
            false,
            false,
            "PM2 not evaluated: no gnomAD annotation present (pm2_absent_when_no_record disabled).".to_string(),
        )
    };

    EvidenceCriterion {
        code,
        direction: EvidenceDirection::Pathogenic,
        strength,
        default_strength: EvidenceStrength::Moderate,
        met,
        evaluated,
        summary,
        details: serde_json::Value::Object(details),
    }
}

/// PM3: For recessive disorders, detected in trans with a pathogenic variant.
///
/// Implements the **ClinGen SVI PM3 v1.0** points-based scoring framework.
/// Each qualifying companion / homozygous occurrence contributes a point
/// value depending on phasing × variant classification:
///
/// | Scenario | Points |
/// |----------|--------|
/// | Confirmed in-trans + co-occurring **Pathogenic** companion | 1.0 |
/// | Confirmed in-trans + co-occurring **Likely Pathogenic** | 0.5 |
/// | Phase unknown + co-occurring Pathogenic | 0.5 |
/// | Phase unknown + co-occurring Likely Pathogenic | 0.25 |
/// | Homozygous occurrence (proband hom-alt) | 0.5 each, capped at 1.0 |
///
/// The total point value maps to PM3 strength:
///
/// | Total | Strength |
/// |-------|----------|
/// | < 0.5 | not met |
/// | ≥ 0.5 | PM3_Supporting |
/// | ≥ 1.0 | PM3 (Moderate) |
/// | ≥ 2.0 | PM3_Strong |
/// | ≥ 4.0 | PM3_VeryStrong |
///
/// Companions in cis with a pathogenic variant are excluded (those count
/// toward BP2 instead). Requires AR inheritance from OMIM.
fn evaluate_pm3(
    input: &ClassificationInput,
    _config: &AcmgConfig,
) -> EvidenceCriterion {
    let mut details = serde_json::Map::new();

    // Recessive inheritance gate.
    let is_recessive = input
        .omim
        .as_ref()
        .is_some_and(|o| o.has_recessive_inheritance());
    details.insert("is_recessive_gene".into(), serde_json::json!(is_recessive));

    if !is_recessive {
        return mk_pm3(
            "PM3".to_string(),
            EvidenceStrength::Moderate,
            false,
            true,
            "Gene does not have autosomal recessive inheritance (PM3 requires recessive disorder)".to_string(),
            details,
        );
    }

    let proband = input.proband_genotype.as_ref();
    let proband_het = proband.is_some_and(|g| g.is_het);
    let proband_hom_alt = proband.is_some_and(|g| g.is_hom_alt);
    details.insert("proband_het".into(), serde_json::json!(proband_het));
    details.insert("proband_hom_alt".into(), serde_json::json!(proband_hom_alt));

    if !proband_het && !proband_hom_alt {
        return mk_pm3(
            "PM3".to_string(),
            EvidenceStrength::Moderate,
            false,
            proband.is_some(),
            if proband.is_some() {
                "Proband is neither het nor hom-alt for this variant (PM3 requires presence)".to_string()
            } else {
                "Proband genotype not available; PM3 requires trio VCF for compound-het analysis".to_string()
            },
            details,
        );
    }

    // Score each contributing observation.
    let mut total: f64 = 0.0;
    let mut breakdown: Vec<String> = Vec::new();

    // Homozygous occurrence (proband hom-alt) earns 0.5 pt, capped at 1.0
    // total across all hom contributions. We model this single proband as one
    // hom occurrence (so the 1.0 cap can never bind yet); a full pedigree
    // workflow would aggregate multiple hom occurrences against that cap.
    if proband_hom_alt {
        let pts: f64 = 0.5;
        total += pts;
        breakdown.push(format!("homozygous_proband:+{:.2}", pts));
    }

    // Compound-het companions in trans / phase-unknown.
    for cv in &input.companion_variants {
        if !cv.proband_het {
            continue;
        }
        // In-cis companions go to BP2, not PM3.
        if cv.is_in_trans == Some(false) {
            continue;
        }
        let confirmed_trans = cv.is_in_trans == Some(true);
        let pts = match (confirmed_trans, cv.is_clinvar_pathogenic, cv.is_clinvar_likely_pathogenic) {
            (true, true, _) => 1.0,
            (true, _, true) => 0.5,
            (false, true, _) => 0.5,
            (false, _, true) => 0.25,
            _ => 0.0,
        };
        if pts == 0.0 {
            continue;
        }
        let label = match (confirmed_trans, cv.is_clinvar_pathogenic, cv.is_clinvar_likely_pathogenic) {
            (true, true, _) => "trans+P",
            (true, _, true) => "trans+LP",
            (false, true, _) => "unphased+P",
            (false, _, true) => "unphased+LP",
            _ => "skipped",
        };
        let label = if let Some(ref hgvs) = cv.hgvsc {
            format!("{}({}):+{:.2}", label, hgvs, pts)
        } else {
            format!("{}:+{:.2}", label, pts)
        };
        breakdown.push(label);
        total += pts;
    }

    details.insert("total_points".into(), serde_json::json!(total));
    details.insert("breakdown".into(), serde_json::json!(breakdown));

    let (strength, code) = if total >= 4.0 {
        (EvidenceStrength::VeryStrong, "PM3_Very_Strong".to_string())
    } else if total >= 2.0 {
        (EvidenceStrength::Strong, "PM3_Strong".to_string())
    } else if total >= 1.0 {
        (EvidenceStrength::Moderate, "PM3".to_string())
    } else if total >= 0.5 {
        (EvidenceStrength::Supporting, "PM3_Supporting".to_string())
    } else {
        return mk_pm3(
            "PM3".to_string(),
            EvidenceStrength::Moderate,
            false,
            true,
            "PM3 points = 0; no qualifying compound-het / homozygous observation".to_string(),
            details,
        );
    };

    let summary = format!(
        "PM3 v1.0 points = {:.2} → {} ({})",
        total,
        strength.as_str(),
        breakdown.join(", ")
    );
    mk_pm3(code, strength, true, true, summary, details)
}

fn mk_pm3(
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
        default_strength: EvidenceStrength::Moderate,
        met,
        evaluated,
        summary,
        details: serde_json::Value::Object(details),
    }
}

/// PM4: Protein length changes due to in-frame deletions/insertions in non-repeat region,
/// or stop-loss variants.
fn evaluate_pm4(
    input: &ClassificationInput,
    _config: &AcmgConfig,
) -> EvidenceCriterion {
    let is_length_change = input.consequences.iter().any(|c| {
        matches!(
            c,
            Consequence::InframeInsertion | Consequence::InframeDeletion | Consequence::StopLost
        )
    });

    let mut details = serde_json::Map::new();
    if is_length_change {
        let types: Vec<&str> = input
            .consequences
            .iter()
            .filter(|c| {
                matches!(
                    c,
                    Consequence::InframeInsertion
                        | Consequence::InframeDeletion
                        | Consequence::StopLost
                )
            })
            .map(|c| c.so_term())
            .collect();
        details.insert("consequence_types".into(), serde_json::json!(types));
    }

    let summary = if is_length_change {
        "Protein length-changing variant (in-frame indel or stop-loss)".to_string()
    } else {
        "Not a protein length-changing variant".to_string()
    };

    EvidenceCriterion {
        code: "PM4".to_string(),
        direction: EvidenceDirection::Pathogenic,
        strength: EvidenceStrength::Moderate,
        default_strength: EvidenceStrength::Moderate,
        met: is_length_change,
        evaluated: true,
        summary,
        details: serde_json::Value::Object(details),
    }
}

/// PM5: Novel missense change at an amino acid residue where a different pathogenic
/// missense change has been seen before.
///
/// Uses the ClinVar protein-position index to check if pathogenic variants
/// with a DIFFERENT amino acid change exist at the same protein position.
fn evaluate_pm5(
    input: &ClassificationInput,
    _config: &AcmgConfig,
) -> EvidenceCriterion {
    let is_missense = input
        .consequences
        .iter()
        .any(|c| matches!(c, Consequence::MissenseVariant));

    let mut details = serde_json::Map::new();

    if !is_missense {
        return EvidenceCriterion {
            code: "PM5".to_string(),
            direction: EvidenceDirection::Pathogenic,
            strength: EvidenceStrength::Moderate,
            default_strength: EvidenceStrength::Moderate,
            met: false,
            evaluated: true,
            summary: "Not a missense variant".to_string(),
            details: serde_json::Value::Object(details),
        };
    }

    let (prot_pos, _ref_aa, alt_aa) = match (&input.protein_position, &input.amino_acids) {
        (Some(pos), Some((r, a))) => (*pos, r.as_str(), a.as_str()),
        _ => {
            return EvidenceCriterion {
                code: "PM5".to_string(),
                direction: EvidenceDirection::Pathogenic,
                strength: EvidenceStrength::Moderate,
                default_strength: EvidenceStrength::Moderate,
                met: false,
                evaluated: false,
                summary: "Protein position or amino acid change not available".to_string(),
                details: serde_json::Value::Object(details),
            };
        }
    };

    details.insert("protein_position".into(), serde_json::json!(prot_pos));
    details.insert("alt_aa".into(), serde_json::json!(alt_aa));

    if let Some(ref cpd) = input.clinvar_protein {
        // Find pathogenic variants at same position with DIFFERENT amino acid change.
        // Substring match `"pathogenic"` includes both `Pathogenic` and
        // `Likely_pathogenic` ClinVar significance values; PM5 fires when ANY
        // P/LP variant exists at the same residue with a different alt AA
        // (LP at the residue is sufficient evidence of a same-codon hotspot).
        let different_aa_matches: Vec<&crate::sa_extract::ClinvarProteinVariant> = cpd
            .protein_variants
            .iter()
            .filter(|v| {
                v.pos == prot_pos
                    && v.alt_aa != alt_aa
                    && v.sig.to_lowercase().contains("pathogenic")
            })
            .collect();

        details.insert(
            "different_aa_pathogenic_count".into(),
            serde_json::json!(different_aa_matches.len()),
        );

        if !different_aa_matches.is_empty() {
            let other_aas: Vec<&str> = different_aa_matches.iter().map(|v| v.alt_aa.as_str()).collect();
            details.insert("other_pathogenic_aas".into(), serde_json::json!(other_aas));

            return EvidenceCriterion {
                code: "PM5".to_string(),
                direction: EvidenceDirection::Pathogenic,
                strength: EvidenceStrength::Moderate,
                default_strength: EvidenceStrength::Moderate,
                met: true,
                evaluated: true,
                summary: format!(
                    "Different pathogenic missense at same residue {} (other AA changes: {})",
                    prot_pos,
                    other_aas.join(", ")
                ),
                details: serde_json::Value::Object(details),
            };
        }

        EvidenceCriterion {
            code: "PM5".to_string(),
            direction: EvidenceDirection::Pathogenic,
            strength: EvidenceStrength::Moderate,
            default_strength: EvidenceStrength::Moderate,
            met: false,
            evaluated: true,
            summary: format!("No different pathogenic missense at position {}", prot_pos),
            details: serde_json::Value::Object(details),
        }
    } else {
        EvidenceCriterion {
            code: "PM5".to_string(),
            direction: EvidenceDirection::Pathogenic,
            strength: EvidenceStrength::Moderate,
            default_strength: EvidenceStrength::Moderate,
            met: false,
            evaluated: false,
            summary: "ClinVar protein-position index not available".to_string(),
            details: serde_json::Value::Object(details),
        }
    }
}

/// PM6: Assumed de novo, but without confirmation of paternity and maternity.
///
/// Fires when the proband carries the variant and only partial parental data is available
/// (one parent specified or one parent fails quality), and the available parent(s) are hom_ref.
/// PS2 and PM6 are mutually exclusive: if full trio data passes quality, PS2 takes priority.
fn evaluate_pm6(
    input: &ClassificationInput,
    config: &AcmgConfig,
) -> EvidenceCriterion {
    let mut details = serde_json::Map::new();

    let trio = match &config.trio {
        Some(t) => t,
        None => {
            return EvidenceCriterion {
                code: "PM6".to_string(),
                direction: EvidenceDirection::Pathogenic,
                strength: EvidenceStrength::Moderate,
                default_strength: EvidenceStrength::Moderate,
                met: false,
                evaluated: false,
                summary: "Requires trio VCF with at least one parent to assess assumed de novo status".to_string(),
                details: serde_json::Value::Null,
            };
        }
    };

    // If both parents are configured, both genotypes are present, and all pass quality,
    // then PS2 should fire instead. PM6 should NOT fire.
    let both_parents_configured = trio.mother.is_some() && trio.father.is_some();
    let both_parents_present = input.mother_genotype.is_some() && input.father_genotype.is_some();
    let min_dp = trio.min_depth;
    let min_gq = trio.min_gq;

    if both_parents_configured && both_parents_present {
        let mother_qc = input.mother_genotype.as_ref().unwrap().passes_quality(min_dp, min_gq);
        let father_qc = input.father_genotype.as_ref().unwrap().passes_quality(min_dp, min_gq);
        let proband_qc = input.proband_genotype.as_ref().is_some_and(|g| g.passes_quality(min_dp, min_gq));
        if mother_qc && father_qc && proband_qc {
            // Full trio with good quality: PS2 applies instead
            return EvidenceCriterion {
                code: "PM6".to_string(),
                direction: EvidenceDirection::Pathogenic,
                strength: EvidenceStrength::Moderate,
                default_strength: EvidenceStrength::Moderate,
                met: false,
                evaluated: true,
                summary: "Both parents available with sufficient quality; PS2 applies instead of PM6".to_string(),
                details: serde_json::Value::Null,
            };
        }
    }

    let proband_gt = match &input.proband_genotype {
        Some(gt) => gt,
        None => {
            return EvidenceCriterion {
                code: "PM6".to_string(),
                direction: EvidenceDirection::Pathogenic,
                strength: EvidenceStrength::Moderate,
                default_strength: EvidenceStrength::Moderate,
                met: false,
                evaluated: false,
                summary: "Proband genotype not available for this variant".to_string(),
                details: serde_json::Value::Null,
            };
        }
    };

    if !proband_gt.carries_variant() {
        return EvidenceCriterion {
            code: "PM6".to_string(),
            direction: EvidenceDirection::Pathogenic,
            strength: EvidenceStrength::Moderate,
            default_strength: EvidenceStrength::Moderate,
            met: false,
            evaluated: true,
            summary: "Proband does not carry the variant allele".to_string(),
            details: serde_json::Value::Null,
        };
    }

    details.insert("proband_carries_variant".into(), serde_json::json!(true));

    // Check available parent(s) -- at least one must be hom_ref and pass quality
    let mut available_parents_ref = 0u32;
    let mut available_parents_count = 0u32;

    if let Some(ref mother_gt) = input.mother_genotype {
        if mother_gt.passes_quality(min_dp, min_gq) {
            available_parents_count += 1;
            details.insert("mother_hom_ref".into(), serde_json::json!(mother_gt.is_hom_ref));
            if mother_gt.is_hom_ref {
                available_parents_ref += 1;
            }
        } else {
            details.insert("mother_quality_fail".into(), serde_json::json!(true));
        }
    }

    if let Some(ref father_gt) = input.father_genotype {
        if father_gt.passes_quality(min_dp, min_gq) {
            available_parents_count += 1;
            details.insert("father_hom_ref".into(), serde_json::json!(father_gt.is_hom_ref));
            if father_gt.is_hom_ref {
                available_parents_ref += 1;
            }
        } else {
            details.insert("father_quality_fail".into(), serde_json::json!(true));
        }
    }

    details.insert("available_parents_passing_qc".into(), serde_json::json!(available_parents_count));
    details.insert("available_parents_hom_ref".into(), serde_json::json!(available_parents_ref));

    if available_parents_count == 0 {
        return EvidenceCriterion {
            code: "PM6".to_string(),
            direction: EvidenceDirection::Pathogenic,
            strength: EvidenceStrength::Moderate,
            default_strength: EvidenceStrength::Moderate,
            met: false,
            evaluated: false,
            summary: "No parent genotype data passing quality thresholds available".to_string(),
            details: serde_json::Value::Object(details),
        };
    }

    let met = available_parents_ref > 0 && available_parents_ref == available_parents_count;
    let summary = if met {
        format!(
            "Assumed de novo: proband carries variant, {} of {} available parent(s) are hom_ref (partial trio confirmation)",
            available_parents_ref, available_parents_count
        )
    } else {
        format!(
            "Not assumed de novo: {} of {} available parent(s) carry the variant",
            available_parents_count - available_parents_ref, available_parents_count
        )
    };

    EvidenceCriterion {
        code: "PM6".to_string(),
        direction: EvidenceDirection::Pathogenic,
        strength: EvidenceStrength::Moderate,
        default_strength: EvidenceStrength::Moderate,
        met,
        evaluated: true,
        summary,
        details: serde_json::Value::Object(details),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_support::minimal_input;
    use crate::sa_extract::{ClinvarData, GnomadData, OmimData};
    use fastvep_core::Impact;

    fn make_input(
        consequences: Vec<Consequence>,
        gnomad: Option<GnomadData>,
    ) -> ClassificationInput {
        ClassificationInput {
            consequences,
            impact: Impact::Moderate,
            gnomad,
            ..minimal_input()
        }
    }

    #[test]
    fn test_pm2_no_gnomad_record_treated_as_absent_by_default() {
        // Default config has `pm2_absent_when_no_record = true`: a missing
        // gnomAD record is interpreted as "variant is absent from gnomAD"
        // per ClinGen SVI v1.0 (PM2 = absent or extremely rare). PM2 fires.
        let input = make_input(vec![Consequence::MissenseVariant], None);
        let result = evaluate_pm2(&input, &AcmgConfig::default());
        assert!(result.met);
        assert!(result.evaluated);
        assert!(result.summary.contains("Absent in gnomAD"));
    }

    #[test]
    fn test_pm2_no_gnomad_record_suppressed_by_common_clinvar_af() {
        // Backstop: no gnomAD record, but ClinVar ships a population AF showing
        // the variant is common (e.g. an indel that failed to match gnomAD, or
        // a coverage gap like ASPM/TOR1AIP1). PM2 must NOT fire.
        let mut input = make_input(vec![Consequence::MissenseVariant], None);
        input.clinvar = Some(ClinvarData {
            af_tgp: Some(0.0113),
            ..Default::default()
        });
        let result = evaluate_pm2(&input, &AcmgConfig::default());
        assert!(!result.met);
        assert!(result.summary.contains("ClinVar population AF"));
    }

    #[test]
    fn test_pm2_no_gnomad_record_still_fires_when_clinvar_af_rare() {
        // A ClinVar AF below the rarity bar does not block PM2.
        let mut input = make_input(vec![Consequence::MissenseVariant], None);
        input.clinvar = Some(ClinvarData {
            af_tgp: Some(0.00001),
            ..Default::default()
        });
        let result = evaluate_pm2(&input, &AcmgConfig::default());
        assert!(result.met);
    }

    #[test]
    fn test_pm2_no_gnomad_record_not_evaluated_when_strict() {
        // Setting `pm2_absent_when_no_record = false` reverts to the
        // strict-coverage stance: PM2 NotEvaluated when no record present.
        // Use this when the loaded gnomAD .osa covers only some input
        // regions and you want PM2 silenced outside that coverage.
        let cfg = AcmgConfig {
            pm2_absent_when_no_record: false,
            ..Default::default()
        };
        let input = make_input(vec![Consequence::MissenseVariant], None);
        let result = evaluate_pm2(&input, &cfg);
        assert!(!result.met);
        assert!(!result.evaluated);
        assert!(result.summary.contains("not evaluated"));
    }

    #[test]
    fn test_pm2_truly_absent_with_gnomad_record_fires() {
        // Real "absent from gnomAD" means a gnomAD record exists with AC=0
        // and AF=0 (the variant was tested for and found absent). PM2 fires.
        let input = make_input(
            vec![Consequence::MissenseVariant],
            Some(GnomadData {
                all_ac: Some(0),
                all_af: Some(0.0),
                ..Default::default()
            }),
        );
        let result = evaluate_pm2(&input, &AcmgConfig::default());
        assert!(result.met);
        assert_eq!(result.strength, EvidenceStrength::Supporting);
        assert_eq!(result.code, "PM2_Supporting");
    }

    #[test]
    fn test_pm2_unknown_inheritance_requires_strict_absence() {
        // Per ClinGen SVI v1.0: AD/unknown-inheritance defaults to strict
        // absence (AC=0). A variant with AF=0.00005 but a record in gnomAD is
        // NOT absent and must NOT fire PM2 in this configuration.
        let input = make_input(
            vec![Consequence::MissenseVariant],
            Some(GnomadData {
                all_af: Some(0.00005),
                all_ac: Some(1),
                ..Default::default()
            }),
        );
        let result = evaluate_pm2(&input, &AcmgConfig::default());
        assert!(!result.met);
    }

    #[test]
    fn test_pm2_ar_gene_under_threshold_fires() {
        // AR gene (OMIM phenotype contains "autosomal recessive") + AF below
        // 0.00007 → PM2_Supporting fires per ClinGen SVI v1.0.
        let mut input = make_input(
            vec![Consequence::MissenseVariant],
            Some(GnomadData {
                all_af: Some(0.00006),
                all_ac: Some(2),
                ..Default::default()
            }),
        );
        input.omim = Some(OmimData {
            mim_number: None,
            phenotypes: Some(vec!["Cystic fibrosis, autosomal recessive".to_string()]),
        });
        let result = evaluate_pm2(&input, &AcmgConfig::default());
        assert!(result.met);
        assert!(result.summary.contains("AR"));
    }

    #[test]
    fn test_pm2_ar_gene_above_threshold_does_not_fire() {
        // AR gene with AF > 0.00007 → PM2 must not fire.
        let mut input = make_input(
            vec![Consequence::MissenseVariant],
            Some(GnomadData {
                all_af: Some(0.0001),
                all_ac: Some(5),
                ..Default::default()
            }),
        );
        input.omim = Some(OmimData {
            mim_number: None,
            phenotypes: Some(vec!["Some disease, autosomal recessive".to_string()]),
        });
        let result = evaluate_pm2(&input, &AcmgConfig::default());
        assert!(!result.met);
    }

    /// A dominant-inheritance variant at the given gnomAD allele frequency.
    fn dominant_at(af: f64, ac: u64) -> ClassificationInput {
        let mut input = make_input(
            vec![Consequence::MissenseVariant],
            Some(GnomadData {
                all_af: Some(af),
                all_ac: Some(ac),
                ..Default::default()
            }),
        );
        input.omim = Some(OmimData {
            mim_number: None,
            phenotypes: Some(vec!["Some disease, autosomal dominant".to_string()]),
        });
        input
    }

    #[test]
    fn test_pm2_ad_gene_singleton_is_extremely_rare_not_present() {
        // A single allele among gnomAD v4's ~800,000 individuals is what PM2
        // was asking about when it said "absent from controls" of a cohort an
        // order of magnitude smaller. Strict absence used to reject this.
        let result = evaluate_pm2(&dominant_at(0.000005, 1), &AcmgConfig::default());
        assert!(result.met, "got: {}", result.summary);
    }

    #[test]
    fn test_pm2_ad_gene_above_the_bar_still_does_not_fire() {
        // The bar has to keep rejecting something, or it is not a bar.
        let result = evaluate_pm2(&dominant_at(0.0005, 400), &AcmgConfig::default());
        assert!(!result.met, "got: {}", result.summary);
    }

    #[test]
    fn test_pm2_ad_default_bar_is_the_measured_one() {
        // Pinned so the default cannot drift away from the figure the config
        // doc comment and docs/ACMG.md both quote.
        assert!((AcmgConfig::default().pm2_ad_af_threshold - 0.00004).abs() < 1e-12);
    }

    #[test]
    fn test_pm2_strict_absence_remains_available() {
        // The literal Richards 2015 reading stays one config line away, for a
        // lab that wants it.
        let strict = AcmgConfig { pm2_ad_af_threshold: 0.0, ..Default::default() };
        assert!(!evaluate_pm2(&dominant_at(0.000005, 1), &strict).met);
        assert!(evaluate_pm2(&dominant_at(0.0, 0), &strict).met);
    }

    #[test]
    fn test_pm2_common_in_gnomad() {
        let input = make_input(
            vec![Consequence::MissenseVariant],
            Some(GnomadData {
                all_af: Some(0.01),
                ..Default::default()
            }),
        );
        let result = evaluate_pm2(&input, &AcmgConfig::default());
        assert!(!result.met);
    }

    #[test]
    fn test_pm2_gene_override_takes_precedence_over_inheritance() {
        // A per-gene pm2_af_threshold override should win even when OMIM says AR.
        let mut config = AcmgConfig::default();
        config.gene_overrides.insert(
            "TEST".to_string(),
            crate::config::GeneOverride {
                mechanism: None,
                ba1_af_threshold: None,
                bs1_af_threshold: None,
                pm2_af_threshold: Some(0.001),
                disabled_criteria: vec![],
                strength_overrides: Default::default(),
                disorders: Default::default(),
            },
        );
        let mut input = make_input(
            vec![Consequence::MissenseVariant],
            Some(GnomadData {
                all_af: Some(0.0005),
                all_ac: Some(20),
                ..Default::default()
            }),
        );
        input.omim = Some(OmimData {
            mim_number: None,
            phenotypes: Some(vec!["Test, autosomal recessive".to_string()]),
        });
        // Override threshold = 0.001; AF = 0.0005 → PM2 fires under override.
        let result = evaluate_pm2(&input, &config);
        assert!(result.met);
        assert!(result.summary.contains("gene_override"));
    }

    #[test]
    fn test_pm2_not_downgraded() {
        // When the SVI downgrade is disabled, PM2 fires at Moderate strength
        // — but still requires real gnomAD data confirming absence (AC=0).
        let config = AcmgConfig {
            pm2_downgrade_to_supporting: false,
            ..Default::default()
        };
        let input = make_input(
            vec![Consequence::MissenseVariant],
            Some(GnomadData {
                all_ac: Some(0),
                all_af: Some(0.0),
                ..Default::default()
            }),
        );
        let result = evaluate_pm2(&input, &config);
        assert!(result.met);
        assert_eq!(result.strength, EvidenceStrength::Moderate);
        assert_eq!(result.code, "PM2");
    }

    #[test]
    fn test_pm4_inframe_deletion() {
        let input = make_input(vec![Consequence::InframeDeletion], None);
        let result = evaluate_pm4(&input, &AcmgConfig::default());
        assert!(result.met);
    }

    #[test]
    fn test_pm4_stop_lost() {
        let input = make_input(vec![Consequence::StopLost], None);
        let result = evaluate_pm4(&input, &AcmgConfig::default());
        assert!(result.met);
    }

    #[test]
    fn test_pm4_missense_not_met() {
        let input = make_input(vec![Consequence::MissenseVariant], None);
        let result = evaluate_pm4(&input, &AcmgConfig::default());
        assert!(!result.met);
    }

    // ── PM3 v1.0 points scoring ────────────────────────────────────────

    use crate::sa_extract::{CompanionVariant, GenotypeInfo};

    fn ar_input_with_proband(het: bool, hom_alt: bool) -> ClassificationInput {
        let mut input = make_input(vec![Consequence::MissenseVariant], None);
        input.omim = Some(OmimData {
            mim_number: None,
            phenotypes: Some(vec!["Cystic fibrosis, autosomal recessive".to_string()]),
        });
        input.proband_genotype = Some(GenotypeInfo {
            is_het: het,
            is_hom_alt: hom_alt,
            is_hom_ref: !het && !hom_alt,
            is_missing: false,
            is_phased: false,
            depth: Some(30),
            quality: Some(50),
            alt_allele_index: if het || hom_alt { Some(1) } else { None },
        });
        input
    }

    fn cv(p: bool, lp: bool, in_trans: Option<bool>, het: bool) -> CompanionVariant {
        CompanionVariant {
            is_clinvar_pathogenic: p,
            is_clinvar_likely_pathogenic: lp,
            is_in_trans: in_trans,
            proband_het: het,
            hgvsc: None,
        }
    }

    #[test]
    fn test_pm3_not_recessive_gene_does_not_fire() {
        let input = make_input(vec![Consequence::MissenseVariant], None);
        let r = evaluate_pm3(&input, &AcmgConfig::default());
        assert!(!r.met);
    }

    #[test]
    fn test_pm3_in_trans_pathogenic_moderate() {
        // 1 confirmed in-trans + Pathogenic = 1.0 pt → PM3 (Moderate)
        let mut input = ar_input_with_proband(true, false);
        input.companion_variants = vec![cv(true, false, Some(true), true)];
        let r = evaluate_pm3(&input, &AcmgConfig::default());
        assert!(r.met);
        assert_eq!(r.strength, EvidenceStrength::Moderate);
        assert_eq!(r.code, "PM3");
    }

    #[test]
    fn test_pm3_in_trans_lp_supporting() {
        // 1 confirmed in-trans + LP = 0.5 pt → PM3_Supporting
        let mut input = ar_input_with_proband(true, false);
        input.companion_variants = vec![cv(false, true, Some(true), true)];
        let r = evaluate_pm3(&input, &AcmgConfig::default());
        assert_eq!(r.strength, EvidenceStrength::Supporting);
        assert_eq!(r.code, "PM3_Supporting");
    }

    #[test]
    fn test_pm3_unphased_pathogenic_supporting() {
        // 1 phase-unknown + Pathogenic = 0.5 pt → PM3_Supporting
        let mut input = ar_input_with_proband(true, false);
        input.companion_variants = vec![cv(true, false, None, true)];
        let r = evaluate_pm3(&input, &AcmgConfig::default());
        assert_eq!(r.strength, EvidenceStrength::Supporting);
    }

    #[test]
    fn test_pm3_two_in_trans_pathogenic_strong() {
        // 2 confirmed in-trans + P = 2.0 pt → PM3_Strong
        let mut input = ar_input_with_proband(true, false);
        input.companion_variants = vec![
            cv(true, false, Some(true), true),
            cv(true, false, Some(true), true),
        ];
        let r = evaluate_pm3(&input, &AcmgConfig::default());
        assert_eq!(r.strength, EvidenceStrength::Strong);
        assert_eq!(r.code, "PM3_Strong");
    }

    #[test]
    fn test_pm3_in_cis_does_not_score() {
        // In-cis companion is excluded from PM3 (it's a BP2 case).
        let mut input = ar_input_with_proband(true, false);
        input.companion_variants = vec![cv(true, false, Some(false), true)];
        let r = evaluate_pm3(&input, &AcmgConfig::default());
        assert!(!r.met);
    }

    #[test]
    fn test_pm3_homozygous_proband_capped_at_one_point() {
        // Hom-alt proband alone earns 0.5 pt → PM3_Supporting.
        let input = ar_input_with_proband(false, true);
        let r = evaluate_pm3(&input, &AcmgConfig::default());
        assert_eq!(r.strength, EvidenceStrength::Supporting);
    }

    #[test]
    fn test_pm3_homozygous_plus_in_trans_p_combines() {
        // Hom-alt (0.5) + in-trans P (1.0) = 1.5 pt → PM3 (Moderate)
        let mut input = ar_input_with_proband(false, true);
        input.companion_variants = vec![cv(true, false, Some(true), true)];
        let r = evaluate_pm3(&input, &AcmgConfig::default());
        assert_eq!(r.strength, EvidenceStrength::Moderate);
    }

    #[test]
    fn test_pm3_four_in_trans_p_very_strong() {
        let mut input = ar_input_with_proband(true, false);
        input.companion_variants = vec![
            cv(true, false, Some(true), true),
            cv(true, false, Some(true), true),
            cv(true, false, Some(true), true),
            cv(true, false, Some(true), true),
        ];
        let r = evaluate_pm3(&input, &AcmgConfig::default());
        assert_eq!(r.strength, EvidenceStrength::VeryStrong);
        assert_eq!(r.code, "PM3_Very_Strong");
    }

    // ── B7: gene-disease validity gates PM1 ──────────────────────────────

    /// A missense variant sitting in a dense cluster of ClinVar pathogenic
    /// variants: PM1's hotspot condition is satisfied outright.
    fn hotspot_missense(gene: &str) -> ClassificationInput {
        use crate::sa_extract::{ClinvarProteinData, ClinvarProteinVariant};
        let neighbours = (1..=4)
            .map(|i| ClinvarProteinVariant {
                pos: 100 + i,
                ref_aa: "A".into(),
                alt_aa: "V".into(),
                sig: "Pathogenic".into(),
                n: 2,
            })
            .collect();
        ClassificationInput {
            consequences: vec![Consequence::MissenseVariant],
            impact: Impact::Moderate,
            gene_symbol: Some(gene.to_string()),
            protein_position: Some(100),
            amino_acids: Some(("A".into(), "T".into())),
            clinvar_protein: Some(ClinvarProteinData {
                protein_variants: neighbours,
                benign_indexed: true,
            ..Default::default()
            }),
            ..minimal_input()
        }
    }

    #[test]
    fn test_pm1_fires_on_a_hotspot_without_a_gene_disease_source() {
        assert!(evaluate_pm1(&hotspot_missense("ARMC9"), &AcmgConfig::default()).met);
    }

    #[test]
    fn test_pm1_blocked_when_clingen_curated_the_gene_as_invalid() {
        // The cluster PM1 reads is made of ClinVar assertions in a gene whose
        // disease relationship curation declined to accept. Counting them as
        // evidence of a critical region assumes the conclusion.
        let mut input = hotspot_missense("ARMC9");
        input.omim = Some(OmimData {
            mim_number: Some(0),
            phenotypes: Some(vec![
                "some proposed disease (ClinGen Disputed/AD, MONDO:0000001)".into(),
            ]),
        });
        let r = evaluate_pm1(&input, &AcmgConfig::default());
        assert!(!r.met);
        assert!(!r.evaluated);
        assert!(
            r.summary.contains("no_valid_gene_disease_relationship"),
            "got: {}",
            r.summary
        );
    }

    #[test]
    fn test_pm1_fires_for_a_gene_the_source_lists() {
        let mut input = hotspot_missense("BRCA1");
        input.omim = Some(OmimData {
            mim_number: Some(0),
            phenotypes: Some(vec!["hereditary breast cancer (ClinGen Definitive/AD)".into()]),
        });
        assert!(evaluate_pm1(&input, &AcmgConfig::default()).met);
    }

    #[test]
    fn test_pm1_survives_for_a_gene_clingen_has_not_curated() {
        let mut input = hotspot_missense("SPAST");
        input.omim = None;
        assert!(evaluate_pm1(&input, &AcmgConfig::default()).met);
    }

    // ── PM1's "without benign variation" half ────────────────────────────

    /// Add `n` benign missense entries inside PM1's window.
    fn with_benign_neighbours(mut input: ClassificationInput, n: u64) -> ClassificationInput {
        use crate::sa_extract::ClinvarProteinVariant;
        let cpd = input.clinvar_protein.as_mut().expect("helper builds one");
        for i in 0..n {
            cpd.protein_variants.push(ClinvarProteinVariant {
                pos: 100 + i,
                ref_aa: "A".into(),
                alt_aa: "G".into(),
                sig: if i % 2 == 0 { "Benign".into() } else { "Likely_benign".into() },
                n: 1,
            });
        }
        input
    }

    #[test]
    fn test_pm1_does_not_fire_where_clinvar_has_benign_variation() {
        // Richards 2015 asks for a hot spot "without benign variation". The
        // pathogenic cluster is untouched; what changes is that the same window
        // has been looked at from the other side too.
        let input = with_benign_neighbours(hotspot_missense("CHD7"), 1);
        let r = evaluate_pm1(&input, &AcmgConfig::default());
        assert!(!r.met);
        assert!(r.evaluated, "declined on the evidence, not for lack of it");
        assert!(r.summary.contains("without benign variation"), "got: {}", r.summary);
    }

    #[test]
    fn test_pm1_benign_tolerance_is_configurable() {
        let input = with_benign_neighbours(hotspot_missense("CHD7"), 2);
        let tolerant = AcmgConfig {
            pm1_max_benign_in_window: 2,
            ..AcmgConfig::default()
        };
        assert!(evaluate_pm1(&input, &tolerant).met);
        let strict = AcmgConfig {
            pm1_max_benign_in_window: 1,
            ..AcmgConfig::default()
        };
        assert!(!evaluate_pm1(&input, &strict).met);
    }

    #[test]
    fn test_pm1_benign_test_is_skipped_for_an_index_without_benign_entries() {
        // An older `.oga` has no benign entries because the builder never wrote
        // them, not because the gene has none. Testing anyway would pass every
        // window and read as evidence that was never gathered.
        let mut input = hotspot_missense("CHD7");
        input.clinvar_protein.as_mut().unwrap().benign_indexed = false;
        let r = evaluate_pm1(&input, &AcmgConfig::default());
        assert!(r.met);
        assert!(
            r.details.get("nearby_benign_count").is_none(),
            "must not report a count it could not measure"
        );
    }

    #[test]
    fn test_pm1_benign_test_excludes_the_variants_own_clinvar_record() {
        // The index is ClinVar. A variant ClinVar calls benign is one of its
        // own benign neighbours, and letting that veto PM1 would be the
        // variant's ClinVar label deciding its own classification - the same
        // circularity the self-match exclusion removed from PS1.
        use crate::sa_extract::ClinvarData;
        let mut input = with_benign_neighbours(hotspot_missense("CHD7"), 1);
        input.clinvar = Some(ClinvarData {
            significance: Some(vec!["Benign".into()]),
            ..Default::default()
        });
        assert!(evaluate_pm1(&input, &AcmgConfig::default()).met);

        let circular = AcmgConfig {
            exclude_self_from_clinvar_evidence: false,
            ..AcmgConfig::default()
        };
        assert!(!evaluate_pm1(&input, &circular).met);
    }

    #[test]
    fn test_pm1_ignores_benign_variation_outside_the_window() {
        use crate::sa_extract::ClinvarProteinVariant;
        let mut input = hotspot_missense("CHD7");
        let window = AcmgConfig::default().pm1_hotspot_window;
        input
            .clinvar_protein
            .as_mut()
            .unwrap()
            .protein_variants
            .push(ClinvarProteinVariant {
                pos: 100 + window + 1,
                ref_aa: "A".into(),
                alt_aa: "G".into(),
                sig: "Benign".into(),
                n: 1,
            });
        assert!(evaluate_pm1(&input, &AcmgConfig::default()).met);
    }
}
