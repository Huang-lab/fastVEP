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
        let nearby_pathogenic: usize = cpd
            .protein_variants
            .iter()
            .filter(|v| v.pos >= low && v.pos <= high && v.sig.to_lowercase().contains("pathogenic"))
            .count();

        details.insert("nearby_pathogenic_count".into(), serde_json::json!(nearby_pathogenic));

        let met = nearby_pathogenic >= threshold as usize;
        let summary = if met {
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
        .map_or(false, |o| o.has_recessive_inheritance());
    let is_dominant = input
        .omim
        .as_ref()
        .map_or(false, |o| o.has_dominant_inheritance());

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

    let (met, summary) = if let Some(ref gnomad) = input.gnomad {
        let af = gnomad.all_af.unwrap_or(0.0);
        let ac = gnomad.all_ac.unwrap_or(0);
        details.insert("gnomad_allAf".into(), serde_json::json!(af));
        details.insert("gnomad_allAc".into(), serde_json::json!(ac));

        // For strict absence (threshold = 0.0), require AC = 0 AND AF = 0.
        // For non-zero thresholds (e.g. AR 0.00007), allow AF ≤ threshold.
        let met = if threshold == 0.0 {
            ac == 0 && af == 0.0
        } else {
            af <= threshold
        };
        let summary = if met {
            format!(
                "{} in gnomAD (AF={:.6}, AC={}, threshold={:.6}, inheritance={})",
                if threshold == 0.0 { "Absent" } else { "Rare" },
                af,
                ac,
                threshold,
                inheritance_basis
            )
        } else {
            format!(
                "Not rare enough in gnomAD (AF={:.6}, AC={}, threshold={:.6}, inheritance={})",
                af, ac, threshold, inheritance_basis
            )
        };
        (met, summary)
    } else {
        // Absent from gnomAD entirely (no record at all).
        details.insert("gnomad_allAf".into(), serde_json::Value::Null);
        (
            true,
            format!(
                "Absent from gnomAD (no record) — meets PM2 under {} rule",
                inheritance_basis
            ),
        )
    };

    EvidenceCriterion {
        code,
        direction: EvidenceDirection::Pathogenic,
        strength,
        default_strength: EvidenceStrength::Moderate,
        met,
        evaluated: true,
        summary,
        details: serde_json::Value::Object(details),
    }
}

/// PM3: For recessive disorders, detected in trans with a pathogenic variant.
///
/// Requires: gene has recessive inheritance (from OMIM), proband is het for this variant,
/// and a companion variant in the same gene is ClinVar pathogenic and proband is het for it.
/// If phased, requires is_in_trans=true. If unphased, PM3 still applies with a limitation note.
fn evaluate_pm3(
    input: &ClassificationInput,
    _config: &AcmgConfig,
) -> EvidenceCriterion {
    let mut details = serde_json::Map::new();

    // Check if gene has recessive inheritance
    let is_recessive = input
        .omim
        .as_ref()
        .map_or(false, |o| o.has_recessive_inheritance());
    details.insert("is_recessive_gene".into(), serde_json::json!(is_recessive));

    if !is_recessive {
        return EvidenceCriterion {
            code: "PM3".to_string(),
            direction: EvidenceDirection::Pathogenic,
            strength: EvidenceStrength::Moderate,
            default_strength: EvidenceStrength::Moderate,
            met: false,
            evaluated: true,
            summary: "Gene does not have autosomal recessive inheritance (PM3 requires recessive disorder)".to_string(),
            details: serde_json::Value::Object(details),
        };
    }

    // Check if proband is het for this variant
    let proband_het = input
        .proband_genotype
        .as_ref()
        .map_or(false, |g| g.is_het);
    details.insert("proband_het".into(), serde_json::json!(proband_het));

    if !proband_het {
        return EvidenceCriterion {
            code: "PM3".to_string(),
            direction: EvidenceDirection::Pathogenic,
            strength: EvidenceStrength::Moderate,
            default_strength: EvidenceStrength::Moderate,
            met: false,
            evaluated: input.proband_genotype.is_some(),
            summary: if input.proband_genotype.is_some() {
                "Proband is not heterozygous for this variant (PM3 requires het in compound-het context)".to_string()
            } else {
                "Proband genotype not available; requires trio VCF for compound-het analysis".to_string()
            },
            details: serde_json::Value::Object(details),
        };
    }

    // Check companion variants: need at least one that is ClinVar pathogenic and proband is het
    let qualifying_companions: Vec<&crate::sa_extract::CompanionVariant> = input
        .companion_variants
        .iter()
        .filter(|cv| cv.is_clinvar_pathogenic && cv.proband_het)
        .collect();

    details.insert(
        "qualifying_companion_count".into(),
        serde_json::json!(qualifying_companions.len()),
    );

    if qualifying_companions.is_empty() {
        let companion_total = input.companion_variants.len();
        return EvidenceCriterion {
            code: "PM3".to_string(),
            direction: EvidenceDirection::Pathogenic,
            strength: EvidenceStrength::Moderate,
            default_strength: EvidenceStrength::Moderate,
            met: false,
            evaluated: true,
            summary: if companion_total == 0 {
                "No companion variants in the same gene for compound-het analysis".to_string()
            } else {
                format!(
                    "No qualifying companion: {} companion variant(s) in gene, but none are ClinVar pathogenic with proband het",
                    companion_total
                )
            },
            details: serde_json::Value::Object(details),
        };
    }

    // Check phase information
    let has_phased_trans = qualifying_companions
        .iter()
        .any(|cv| cv.is_in_trans == Some(true));
    let has_phased_cis = qualifying_companions
        .iter()
        .any(|cv| cv.is_in_trans == Some(false));
    let all_unphased = qualifying_companions
        .iter()
        .all(|cv| cv.is_in_trans.is_none());

    details.insert("has_phased_trans".into(), serde_json::json!(has_phased_trans));
    details.insert("has_phased_cis".into(), serde_json::json!(has_phased_cis));
    details.insert("all_unphased".into(), serde_json::json!(all_unphased));

    // Collect companion HGVSc for reporting
    let companion_ids: Vec<String> = qualifying_companions
        .iter()
        .filter_map(|cv| cv.hgvsc.clone())
        .collect();
    if !companion_ids.is_empty() {
        details.insert("companion_hgvsc".into(), serde_json::json!(companion_ids));
    }

    if has_phased_cis && !has_phased_trans {
        // All qualifying companions are in cis -- not PM3 (but may be BP2)
        return EvidenceCriterion {
            code: "PM3".to_string(),
            direction: EvidenceDirection::Pathogenic,
            strength: EvidenceStrength::Moderate,
            default_strength: EvidenceStrength::Moderate,
            met: false,
            evaluated: true,
            summary: "Phased data shows companion pathogenic variant is in cis (same haplotype), not trans".to_string(),
            details: serde_json::Value::Object(details),
        };
    }

    let met = has_phased_trans || all_unphased;
    let summary = if has_phased_trans {
        format!(
            "Compound heterozygote: in trans with ClinVar pathogenic variant in recessive gene (phased, {} companion(s))",
            qualifying_companions.len()
        )
    } else {
        format!(
            "Compound heterozygote: co-occurs with ClinVar pathogenic variant in recessive gene (unphased data, {} companion(s) -- phase not confirmed)",
            qualifying_companions.len()
        )
    };

    EvidenceCriterion {
        code: "PM3".to_string(),
        direction: EvidenceDirection::Pathogenic,
        strength: EvidenceStrength::Moderate,
        default_strength: EvidenceStrength::Moderate,
        met,
        evaluated: true,
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
        // Find pathogenic variants at same position with DIFFERENT amino acid change
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
        let proband_qc = input.proband_genotype.as_ref().map_or(false, |g| g.passes_quality(min_dp, min_gq));
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
    use crate::sa_extract::{GnomadData, OmimData};
    use fastvep_core::Impact;

    fn make_input(
        consequences: Vec<Consequence>,
        gnomad: Option<GnomadData>,
    ) -> ClassificationInput {
        ClassificationInput {
            consequences,
            impact: Impact::Moderate,
            gene_symbol: Some("TEST".to_string()),
            is_canonical: true,
            amino_acids: None,
            protein_position: None,
            gnomad,
            clinvar: None,
            revel: None,
            splice_ai: None,
            dbnsfp: None,
            phylop: None,
            gerp: None,
            gene_constraints: None,
            omim: None,
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
    fn test_pm2_absent_from_gnomad() {
        let input = make_input(vec![Consequence::MissenseVariant], None);
        let result = evaluate_pm2(&input, &AcmgConfig::default());
        assert!(result.met);
        assert_eq!(result.strength, EvidenceStrength::Supporting); // Downgraded per SVI
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

    #[test]
    fn test_pm2_ad_gene_with_one_allele_does_not_fire() {
        // AD gene + any AC > 0 → not absent → PM2 must not fire under SVI v1.0.
        let mut input = make_input(
            vec![Consequence::MissenseVariant],
            Some(GnomadData {
                all_af: Some(0.000005),
                all_ac: Some(1),
                ..Default::default()
            }),
        );
        input.omim = Some(OmimData {
            mim_number: None,
            phenotypes: Some(vec!["Some disease, autosomal dominant".to_string()]),
        });
        let result = evaluate_pm2(&input, &AcmgConfig::default());
        assert!(!result.met);
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
                bs1_af_threshold: None,
                pm2_af_threshold: Some(0.001),
                disabled_criteria: vec![],
                strength_overrides: Default::default(),
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
        let mut config = AcmgConfig::default();
        config.pm2_downgrade_to_supporting = false;
        let input = make_input(vec![Consequence::MissenseVariant], None);
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
}
