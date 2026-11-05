//! ACMG-AMP variant classification engine for fastVEP.
//!
//! Implements the 28 evidence criteria from Richards et al. 2015 (PP5 and BP6
//! are off by default)
//! (Standards and guidelines for the interpretation of sequence variants)
//! and produces a 5-tier classification: Pathogenic, Likely Pathogenic,
//! Uncertain Significance (VUS), Likely Benign, Benign.
//!
//! Incorporates ClinGen SVI (Sequence Variant Interpretation) working group
//! recommendations including calibrated REVEL thresholds for PP3/BP4 and
//! PM2 downgrade to Supporting strength.

pub mod combiner;
pub mod config;
pub mod criteria;
pub mod disease;
pub mod functional;
pub mod sa_extract;
#[cfg(test)]
mod test_support;
pub mod types;

pub use config::{AcmgConfig, TrioConfig};
pub use disease::{ClingenTag, DiseaseEntry, DisorderContext, Inheritance, ModeOfInheritance};
pub use functional::{FunctionalCriterion, FunctionalEvidence, FunctionalEvidenceIndex};
pub use sa_extract::{
    extract_classification_input, is_pure_insertion, reads_exon_edge, reads_hgvs_placement,
    same_splice_position_pathogenic, ClassificationInput, CompanionVariant, GenotypeInfo,
};
pub use types::{
    AcmgClassification, AcmgResult, DisorderClassification, DisorderRef, EvidenceCounts,
    EvidenceCriterion,
};

/// Classify a variant using the ACMG-AMP framework.
///
/// Takes a `ClassificationInput` (extracted from pipeline annotation data)
/// and an `AcmgConfig` (with thresholds and gene-specific overrides).
///
/// Returns an `AcmgResult` containing the classification, all evaluated
/// criteria, triggered rule, and evidence counts.
///
/// For a gene whose established diseases fall under more than one mode of
/// inheritance, the result also carries [`AcmgResult::per_disorder`]: the same
/// variant classified once per mode, against only the diseases that have it.
/// The headline call is the same either way.
pub fn classify(input: &ClassificationInput, config: &AcmgConfig) -> AcmgResult {
    let mut result = classify_against_gene(input, config);
    result.per_disorder = classify_each_disorder(input, config);
    result
}

/// One classification of the variant against everything `input` says about its
/// gene.
fn classify_against_gene(input: &ClassificationInput, config: &AcmgConfig) -> AcmgResult {
    // Evaluate every criterion (26 by default: PP5 and BP6 are off)
    let criteria = criteria::evaluate_all_criteria(input, config);

    // Count met criteria by direction/strength
    let counts = EvidenceCounts::from_criteria(&criteria);

    // Apply combination rules
    let (classification, triggered_rule) = combiner::combine(&criteria, config);

    AcmgResult {
        shorthand: classification.shorthand().to_string(),
        classification,
        criteria,
        triggered_rule,
        counts,
        per_disorder: Vec::new(),
    }
}

/// The variant classified against each mode of inheritance of its gene, or
/// nothing when the gene has one.
///
/// A disorder context is the gene's disease annotation cut down to the
/// established diseases that share a mode of inheritance, and every criterion
/// that reads inheritance or gene-disease validity reads it from that
/// annotation, so classifying against a context is classifying the same input
/// with the annotation replaced. Nothing else about the input changes: the
/// variant, its frequencies and the gene's constraint are not disorder-specific.
fn classify_each_disorder(
    input: &ClassificationInput,
    config: &AcmgConfig,
) -> Vec<DisorderClassification> {
    let Some(annotation) = input.omim.as_ref() else {
        return Vec::new();
    };
    // Asked of every variant of every annotated gene, so it is answered before
    // anything is built.
    if !annotation.has_several_disorder_contexts() {
        return Vec::new();
    }
    annotation
        .disorder_contexts()
        .iter()
        .map(|context| {
            let mut contextual = input.clone();
            contextual.omim = Some(context.restricted_annotation());
            DisorderClassification {
                inheritance: context.inheritance_code().to_string(),
                disorders: context
                    .entries
                    .iter()
                    .map(|entry| DisorderRef {
                        label: entry.label.to_string(),
                        mondo: entry.clingen.and_then(|tag| tag.mondo).map(str::to_string),
                        validity: entry.clingen.map(|tag| tag.classification.to_string()),
                    })
                    .collect(),
                result: classify_against_gene(&contextual, config),
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sa_extract::*;
    use crate::test_support::minimal_input;
    use fastvep_core::{Consequence, Impact};

    /// Helper to build a ClassificationInput with common defaults.
    fn make_input(
        consequences: Vec<Consequence>,
        impact: Impact,
        gene_symbol: &str,
    ) -> ClassificationInput {
        ClassificationInput {
            consequences,
            impact,
            gene_symbol: Some(gene_symbol.to_string()),
            ..minimal_input()
        }
    }

    #[test]
    fn test_classify_common_variant_benign() {
        let mut input = make_input(vec![Consequence::MissenseVariant], Impact::Moderate, "TEST");
        input.gnomad = Some(GnomadData {
            all_af: Some(0.10),
            afr_af: Some(0.15),
            all_an: Some(100_000),
            ..Default::default()
        });
        let result = classify(&input, &AcmgConfig::default());
        assert_eq!(result.classification, AcmgClassification::Benign);
        assert_eq!(result.shorthand, "B");
        // BA1 is standalone benign: -8 points, inside the Benign band on its
        // own. The rule string names the arithmetic under the point system and
        // the rule id under the table, so assert on the evidence instead.
        assert!(result.criteria.iter().any(|c| c.code == "BA1" && c.met));
        let table = AcmgConfig {
            use_point_system: false,
            ..Default::default()
        };
        assert_eq!(
            classify(&input, &table).triggered_rule.as_deref(),
            Some("BA1")
        );
    }

    #[test]
    fn test_classify_frameshift_lof_gene_likely_pathogenic() {
        let mut input = make_input(vec![Consequence::FrameshiftVariant], Impact::High, "BRCA1");
        // High pLI → PVS1
        input.gene_constraints = Some(GnomadGeneData {
            pli: Some(1.0),
            loeuf: Some(0.03),
            mis_z: Some(2.5),
            syn_z: Some(0.5),
        });
        // Confirmed-absent gnomAD record (AC=0, AF=0) → PM2_Supporting fires.
        // Pre-fix the test relied on PM2 firing when gnomAD was None, which
        // was the buggy "no-data == absent" behavior fixed in PR23.
        input.gnomad = Some(GnomadData {
            all_ac: Some(0),
            all_af: Some(0.0),
            ..Default::default()
        });
        let config = AcmgConfig::default();
        let result = classify(&input, &config);

        // PVS1 (VeryStrong) + PM2_Supporting (Supporting) = PVS=1, PP=1
        // Per ClinGen SVI (Sept 2020): PVS + >=1 PP → Likely Pathogenic
        assert!(result.counts.pathogenic_very_strong >= 1); // PVS1
        assert!(result.counts.pathogenic_supporting >= 1); // PM2_Supporting
        assert_eq!(result.classification, AcmgClassification::LikelyPathogenic);
    }

    #[test]
    fn test_classify_frameshift_revel_does_not_drive_pathogenic() {
        // Pre-PR1, a high REVEL score on a frameshift incorrectly fired
        // PP3_Strong, pushing PVS1 + REVEL → Pathogenic via the PVS+PS rule.
        // Per Pejaver 2022, REVEL is calibrated for missense only and must
        // not contribute to non-missense classification. Without other
        // pathogenic-Strong evidence, PVS1 + PM2_Supporting tops out at LP
        // via the ClinGen SVI PVS+PP rule.
        let mut input = make_input(vec![Consequence::FrameshiftVariant], Impact::High, "BRCA1");
        input.gene_constraints = Some(GnomadGeneData {
            pli: Some(1.0),
            loeuf: Some(0.03),
            mis_z: Some(2.5),
            syn_z: Some(0.5),
        });
        input.revel = Some(RevelData { score: Some(0.95) });
        // Confirmed-absent gnomAD record so PM2_Supporting fires.
        input.gnomad = Some(GnomadData {
            all_ac: Some(0),
            all_af: Some(0.0),
            ..Default::default()
        });

        let config = AcmgConfig::default();
        let result = classify(&input, &config);
        assert_eq!(result.classification, AcmgClassification::LikelyPathogenic);
    }

    #[test]
    fn test_classify_synonymous_no_splice_not_conserved() {
        let mut input = make_input(vec![Consequence::SynonymousVariant], Impact::Low, "TEST");
        input.splice_ai = Some(SpliceAiData {
            ds_ag: Some(0.01),
            ds_al: Some(0.02),
            ds_dg: Some(0.01),
            ds_dl: Some(0.01),
            ..Default::default()
        });
        input.phylop = Some(0.3);
        input.revel = Some(RevelData { score: Some(0.10) });
        // Provide gnomAD data to prevent PM2_Supporting from firing (which would cause conflict)
        input.gnomad = Some(GnomadData {
            all_af: Some(0.005), // Above PM2 threshold (0.0001) but below BS1 (0.01)
            ..Default::default()
        });

        let result = classify(&input, &AcmgConfig::default());
        // BP7 (synonymous + no splice + not conserved): benign supporting
        // BP4_Moderate (REVEL=0.10 <= 0.183): counts as benign strong
        // BS + BP → Likely Benign
        assert_eq!(result.classification, AcmgClassification::LikelyBenign);
    }

    #[test]
    fn test_classify_vus_no_data() {
        let input = make_input(
            vec![Consequence::MissenseVariant],
            Impact::Moderate,
            "UNKNOWN_GENE",
        );
        let result = classify(&input, &AcmgConfig::default());
        // No SA data at all → most criteria not evaluable
        // PM2_Supporting fires (absent from gnomAD)
        // But just 1 supporting isn't enough for any LP rule
        assert_eq!(
            result.classification,
            AcmgClassification::UncertainSignificance
        );
    }

    #[test]
    fn test_classify_conflicting_evidence() {
        // PR9: combiner conflict-gating fix. PVS1 alone + BS1 alone do NOT
        // reach a definite call on either side (PVS1 needs PS/PM/PP to fire
        // a pathogenic rule, and BS1 alone is sub-threshold for Benign).
        // Result is plain VUS without a "Conflicting" label.
        let mut input = make_input(vec![Consequence::FrameshiftVariant], Impact::High, "GENE");
        input.gene_constraints = Some(GnomadGeneData {
            pli: Some(1.0),
            loeuf: Some(0.03),
            ..Default::default()
        });
        input.gnomad = Some(GnomadData {
            all_af: Some(0.02),
            all_hc: Some(0),
            ..Default::default()
        });

        // AF = 0.02 is above the PM2 bar, so PM2 does not fire; BS1 is not
        // evaluated because the record carries no AN. That leaves PVS1 alone,
        // and the two combining schemes genuinely disagree about it.
        //
        // Point system (the default): 8 points, inside the Likely Pathogenic
        // band. A frameshift in a haploinsufficient disease gene is not an
        // uncertain finding, and run v10 had 2,319 ClinVar-Pathogenic variants
        // sitting in VUS on exactly this signature.
        assert_eq!(
            classify(&input, &AcmgConfig::default()).classification,
            AcmgClassification::LikelyPathogenic
        );

        // Richards 2015 Table 5: no row matches one Very Strong alone, so the
        // table returns Uncertain Significance.
        let table = AcmgConfig {
            use_point_system: false,
            ..Default::default()
        };
        let result = classify(&input, &table);
        assert_eq!(
            result.classification,
            AcmgClassification::UncertainSignificance
        );
        assert!(
            result.triggered_rule.is_none()
                || !result
                    .triggered_rule
                    .as_deref()
                    .unwrap_or("")
                    .contains("Conflicting"),
            "PR9 expects no Conflicting label here; got {:?}",
            result.triggered_rule
        );
    }

    #[test]
    fn test_acmg_result_serialization() {
        let input = make_input(vec![Consequence::MissenseVariant], Impact::Moderate, "TEST");
        let result = classify(&input, &AcmgConfig::default());
        let json = serde_json::to_value(&result).unwrap();

        assert!(json.get("classification").is_some());
        assert!(json.get("shorthand").is_some());
        assert!(json.get("criteria").is_some());
        assert!(json.get("counts").is_some());

        let criteria = json.get("criteria").unwrap().as_array().unwrap();
        assert!(!criteria.is_empty());

        // Each criterion should have required fields
        for c in criteria {
            assert!(c.get("code").is_some());
            assert!(c.get("direction").is_some());
            assert!(c.get("strength").is_some());
            assert!(c.get("met").is_some());
            assert!(c.get("evaluated").is_some());
            assert!(c.get("summary").is_some());
        }
    }

    /// The whole chain on one variant: a missense in a splice region that REVEL
    /// calls damaging and SpliceAI calls splice-neutral. PP3 fires from REVEL;
    /// BP4 fires from SpliceAI, because the splice region keeps the splice path
    /// open; and the reconciliation pass drops the BP4, so the benign side does
    /// not collect a criterion that says nothing about the protein.
    #[test]
    fn test_classify_damaging_missense_in_a_splice_region_does_not_collect_a_splice_bp4() {
        let mut input = make_input(
            vec![
                Consequence::MissenseVariant,
                Consequence::SpliceRegionVariant,
            ],
            Impact::Moderate,
            "TEST",
        );
        input.revel = Some(RevelData { score: Some(0.95) });
        input.splice_ai = Some(SpliceAiData {
            ds_al: Some(0.0),
            ..Default::default()
        });
        let result = classify(&input, &AcmgConfig::default());
        let by_code = |prefix: &str| {
            result
                .criteria
                .iter()
                .find(|c| c.code.starts_with(prefix))
                .unwrap_or_else(|| panic!("no {prefix}"))
        };
        assert!(by_code("PP3").met, "{}", by_code("PP3").summary);
        let bp4 = by_code("BP4");
        assert!(!bp4.met, "{}", bp4.summary);
        assert_eq!(bp4.details["bp4_source"], serde_json::json!("spliceai"));
        assert!(bp4.details.get("suppressed_by_reconcile").is_some());
    }

    // ── Genes with more than one mode of inheritance ──

    /// A gene with an autosomal dominant disease and an autosomal recessive one,
    /// as ATM, ABCC8 and BRCA1 are in ClinGen's export.
    fn dual_inheritance_input(gnomad: GnomadData) -> ClassificationInput {
        let mut input = make_input(vec![Consequence::MissenseVariant], Impact::Moderate, "DUAL");
        input.gnomad = Some(gnomad);
        input.omim = Some(OmimData {
            mim_number: None,
            phenotypes: Some(vec![
                "a dominant cancer predisposition (ClinGen Definitive/AD, MONDO:0000001)".into(),
                "a recessive syndrome (ClinGen Definitive/AR, MONDO:0000002)".into(),
                "a disputed one (ClinGen Disputed/XL, MONDO:0000003)".into(),
            ]),
        });
        input
    }

    fn met(result: &AcmgResult, prefix: &str) -> bool {
        result
            .criteria
            .iter()
            .any(|c| c.met && c.code.starts_with(prefix))
    }

    #[test]
    fn test_pm2_takes_the_recessive_bar_in_the_recessive_context_only() {
        // 5e-5 is above the dominant bar (4e-5) and below the recessive one
        // (7e-5). The gene is both, so the headline reads it as unknown and
        // applies the strict bar; each context applies its own.
        let input = dual_inheritance_input(GnomadData {
            all_ac: Some(70),
            all_af: Some(5e-5),
            all_an: Some(1_400_000),
            all_hc: Some(0),
            ..Default::default()
        });
        let result = classify(&input, &AcmgConfig::default());
        assert!(!met(&result, "PM2"), "the headline stays on the strict bar");

        let codes: Vec<_> = result
            .per_disorder
            .iter()
            .map(|d| d.inheritance.as_str())
            .collect();
        assert_eq!(
            codes,
            ["AD", "AR"],
            "the disputed X-linked curation is not a context"
        );
        let by = |code: &str| {
            &result
                .per_disorder
                .iter()
                .find(|d| d.inheritance == code)
                .expect("context")
                .result
        };
        assert!(!met(by("AD"), "PM2"), "5e-5 is above the dominant bar");
        assert!(met(by("AR"), "PM2"), "5e-5 is below the recessive bar");
    }

    #[test]
    fn test_bs2_counts_carriers_in_the_dominant_context_and_homozygotes_in_the_recessive() {
        // 500 carriers among 500 K individuals: more than a dominant disorder
        // can account for, and not a single homozygote.
        let input = dual_inheritance_input(GnomadData {
            all_ac: Some(500),
            all_af: Some(5e-4),
            all_an: Some(1_000_000),
            all_hc: Some(0),
            ..Default::default()
        });
        let result = classify(&input, &AcmgConfig::default());
        assert!(
            !met(&result, "BS2"),
            "the headline takes the arm for a gene that is both"
        );
        let ad = result
            .per_disorder
            .iter()
            .find(|d| d.inheritance == "AD")
            .unwrap();
        let ar = result
            .per_disorder
            .iter()
            .find(|d| d.inheritance == "AR")
            .unwrap();
        assert!(met(&ad.result, "BS2"), "{:#?}", ad.result.criteria);
        assert!(
            !met(&ar.result, "BS2"),
            "no homozygotes, so nothing for the recessive arm"
        );
    }

    #[test]
    fn test_a_context_names_the_diseases_it_was_made_against() {
        let input = dual_inheritance_input(GnomadData::default());
        let result = classify(&input, &AcmgConfig::default());
        let ad = &result.per_disorder[0];
        assert_eq!(ad.disorders.len(), 1);
        assert_eq!(ad.disorders[0].label, "a dominant cancer predisposition");
        assert_eq!(ad.disorders[0].mondo.as_deref(), Some("MONDO:0000001"));
        assert_eq!(ad.disorders[0].validity.as_deref(), Some("Definitive"));
    }

    #[test]
    fn test_a_gene_with_one_mode_of_inheritance_has_no_per_disorder_results() {
        let mut input = dual_inheritance_input(GnomadData::default());
        input.omim = Some(OmimData {
            mim_number: None,
            phenotypes: Some(vec![
                "a (ClinGen Definitive/AD, MONDO:0000001)".into(),
                "b (ClinGen Strong/AD, MONDO:0000002)".into(),
            ]),
        });
        let result = classify(&input, &AcmgConfig::default());
        assert!(result.per_disorder.is_empty());
        let json = serde_json::to_value(&result).unwrap();
        assert!(
            json.get("per_disorder").is_none(),
            "an empty list is not written: a reader of the old shape sees the old shape"
        );
        // No annotation at all is the same.
        input.omim = None;
        assert!(classify(&input, &AcmgConfig::default())
            .per_disorder
            .is_empty());
    }

    #[test]
    fn test_per_disorder_results_do_not_change_the_headline_call() {
        // The headline is what it would be without the feature: classifying the
        // same input with the per-disorder pass skipped gives the same answer.
        let input = dual_inheritance_input(GnomadData {
            all_ac: Some(500),
            all_af: Some(5e-4),
            all_an: Some(1_000_000),
            all_hc: Some(0),
            ..Default::default()
        });
        let config = AcmgConfig::default();
        let with = classify(&input, &config);
        let without = classify_against_gene(&input, &config);
        assert_eq!(with.classification, without.classification);
        assert_eq!(with.shorthand, without.shorthand);
        assert_eq!(with.triggered_rule, without.triggered_rule);
        let codes = |r: &AcmgResult| {
            r.criteria
                .iter()
                .map(|c| (c.code.clone(), c.met))
                .collect::<Vec<_>>()
        };
        assert_eq!(codes(&with), codes(&without));

        // The comparison above classifies one input twice, so a change that made the headline
        // one mode's answer would pass it. This pins what the headline is: the gene's diseases
        // read together, which is neither mode's reading.
        let reading = |r: &AcmgResult| {
            let bs2 = r
                .criteria
                .iter()
                .find(|c| c.code.starts_with("BS2"))
                .expect("BS2 is evaluated");
            (
                bs2.details["omim_dominant"].as_bool(),
                bs2.details["omim_recessive"].as_bool(),
            )
        };
        assert_eq!(reading(&with), (Some(true), Some(true)));
        assert_eq!(
            reading(&with.per_disorder[0].result),
            (Some(true), Some(false))
        );
        assert_eq!(
            reading(&with.per_disorder[1].result),
            (Some(false), Some(true))
        );
    }

    /// One companion in trans with a ClinVar pathogenic variant is PM3's observation for a recessive
    /// disease and BP2's argument against a variant in a fully penetrant dominant one, and a gene that has
    /// both does not get both for it: the headline takes PM3, which asks only that one disease be
    /// recessive, and not BP2, which needs the diseases to be dominant and none recessive; each mode
    /// takes the criterion that is its own.
    #[test]
    fn test_one_companion_in_trans_is_pm3_or_bp2_and_not_both_in_a_gene_with_both_modes() {
        use crate::sa_extract::{CompanionVariant, GenotypeInfo};
        let mut input = dual_inheritance_input(GnomadData::default());
        input.proband_genotype = Some(GenotypeInfo {
            is_het: true,
            is_hom_alt: false,
            is_hom_ref: false,
            is_missing: false,
            is_phased: true,
            depth: Some(30),
            quality: Some(50),
            alt_allele_index: Some(1),
        });
        input.companion_variants = vec![CompanionVariant {
            is_clinvar_pathogenic: true,
            is_clinvar_likely_pathogenic: false,
            is_in_trans: Some(true),
            proband_het: true,
            hgvsc: None,
        }];
        let result = classify(&input, &AcmgConfig::default());
        assert!(met(&result, "PM3"), "the gene has a recessive disease");
        assert!(
            !met(&result, "BP2"),
            "and a recessive disease explains the companion"
        );

        let (ad, ar) = (&result.per_disorder[0], &result.per_disorder[1]);
        assert_eq!(
            (ad.inheritance.as_str(), ar.inheritance.as_str()),
            ("AD", "AR")
        );
        assert!(
            met(&ad.result, "BP2") && !met(&ad.result, "PM3"),
            "the dominant disease: BP2"
        );
        assert!(
            met(&ar.result, "PM3") && !met(&ar.result, "BP2"),
            "the recessive disease: PM3"
        );
    }

    #[test]
    fn test_per_disorder_results_serialise_and_read_back() {
        let input = dual_inheritance_input(GnomadData::default());
        let result = classify(&input, &AcmgConfig::default());
        let json = serde_json::to_value(&result).unwrap();
        let entries = json["per_disorder"].as_array().expect("an array");
        assert_eq!(entries.len(), 2);
        // Flattened: an entry reads like the headline above it.
        for entry in entries {
            assert!(entry.get("classification").is_some());
            assert!(entry.get("shorthand").is_some());
            assert!(entry["criteria"].as_array().is_some_and(|c| !c.is_empty()));
            assert!(entry.get("inheritance").is_some());
            assert!(entry["disorders"].as_array().is_some_and(|d| !d.is_empty()));
            assert!(entry.get("per_disorder").is_none(), "nothing nests");
        }
        let back: AcmgResult = serde_json::from_value(json).unwrap();
        assert_eq!(back.per_disorder.len(), 2);
        assert_eq!(back.per_disorder[1].inheritance, "AR");
    }

    /// The two calls the tag first got wrong, TTR p.Val142Ile (722 gnomAD alleles)
    /// and VWF p.Arg854Gln (7,432), each a record that fires BS1 and, in a gene
    /// read as dominant, BS2's carrier arm too. Two Strong benign criteria are
    /// Benign under the Richards rules, so counting the frequency twice called a
    /// ClinVar-Pathogenic variant Benign. The record below is a synthetic one of
    /// that kind, not either variant's own counts; the variant is a VUS, because
    /// BS1 stands once.
    #[test]
    fn test_a_frequency_counted_by_bs1_is_not_counted_again_by_bs2() {
        let mut input = make_input(vec![Consequence::MissenseVariant], Impact::Moderate, "TTR");
        input.gnomad = Some(GnomadData {
            all_ac: Some(15_000),
            all_hc: Some(50),
            all_an: Some(1_460_000),
            all_af: Some(0.0103),
            afr_af: Some(0.012),
            ..Default::default()
        });
        input.omim = Some(OmimData {
            mim_number: None,
            phenotypes: Some(vec![
                "hereditary amyloidosis (ClinGen Definitive/AD, MONDO:0017132)".into(),
            ]),
        });
        let result = classify(&input, &AcmgConfig::default());
        assert!(met(&result, "BS1"), "{:#?}", result.criteria);
        let bs2 = result.criteria.iter().find(|c| c.code == "BS2").unwrap();
        assert!(!bs2.met);
        assert!(bs2.details.get("suppressed_by_reconcile").is_some());
        assert_eq!(
            result.classification,
            AcmgClassification::UncertainSignificance,
            "{:?}",
            result.triggered_rule
        );
    }
}
