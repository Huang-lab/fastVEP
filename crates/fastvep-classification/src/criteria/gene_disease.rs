//! Gene-level disease knowledge shared by the criteria that presuppose it.
//!
//! Three of the pathogenic criteria are not statements about a variant on its
//! own. PVS1 says "this null variant destroys a gene whose loss causes
//! disease". PP2 says "missense is how this gene causes disease". PM1 says
//! "this residue is where the disease-causing changes cluster". Each rests on
//! a claim about the gene that the variant-level evidence cannot supply, and
//! each becomes unfalsifiable when that claim is untested.
//!
//! Two gates live here, both keyed on the gene rather than the variant:
//!
//! * [`validity_blocker`] - is there an established gene-disease relationship
//!   at all? (B7; ClinGen Gene-Disease Validity, Strehlow 2024.)
//! * [`lof_mechanism_blocker`] - is loss of function the mechanism, or does
//!   this gene cause disease by gain of function? (B6; Abou Tayoun 2018.)
//!
//! Both degrade to "no blocker" when the knowledge they read is absent, so a
//! run without a gene-disease source, or a gene with no curated mechanism,
//! behaves exactly as it did before these gates existed. That matters more
//! here than for the frequency gates: absence of evidence about a gene is the
//! normal state for most of the genome, and treating it as evidence of absence
//! would suppress PVS1 across every gene no curator has reached yet.

use crate::config::AcmgConfig;
use crate::sa_extract::ClassificationInput;

/// A reason the gene's disease relationship has been curated and found wanting,
/// and so cannot support PVS1, PP2 or PM1.
///
/// **Blocks on positive evidence, never on absence.** ClinGen has curated
/// roughly 2,400 genes to Definitive/Strong/Moderate, so "not in the file" is
/// overwhelmingly "nobody has got to it yet" - a statement about ClinGen's
/// progress, not about the gene. Only a gene ClinGen curated and classified as
/// Limited, Disputed, Refuted or No Known Disease Relationship carries evidence
/// that the relationship does not hold.
///
/// The first version of this gate blocked on absence, and run v10 measured what
/// that costs: **1,497 truth-pathogenic null variants lost PVS1**, in genes like
/// SPAST, ABCB11, FLG and LAMB3 that cause disease and simply are not curated
/// yet. The genes it was meant to catch (RYK, GIGYF2) are absent from ClinGen
/// too, so absence-blocking was not even selecting for the right thing - it was
/// selecting for curation coverage.
pub(crate) fn validity_blocker(input: &ClassificationInput, config: &AcmgConfig) -> Option<String> {
    if !config.require_gene_disease_validity {
        return None;
    }
    let omim = input.omim.as_ref()?;
    if !omim.has_no_valid_relationship() {
        return None;
    }
    Some(format!(
        "no_valid_gene_disease_relationship: every ClinGen curation of {} is Limited, Disputed, Refuted or No Known Disease Relationship, so there is no established disease for this variant to be pathogenic for",
        input.gene_symbol.as_deref().unwrap_or("the gene"),
    ))
}

/// A reason loss of function is not how this gene causes disease, and so
/// PVS1's whole premise does not hold.
///
/// Abou Tayoun 2018 makes this the first question in the PVS1 decision tree:
/// "is LOF a known mechanism of disease for this gene?" A gene that causes
/// disease only through gain of function - a constitutively active receptor, a
/// channel that stays open, a protein that acquires a new toxic interaction -
/// is one where knocking the allele out is the *benign* outcome. Haploinsufficiency
/// intolerance measured by pLI does not distinguish the two: a GoF gene under
/// strong purifying selection scores as constrained exactly like a
/// haploinsufficient one, which is why constraint alone must not be allowed to
/// carry PVS1.
///
/// Consults the mechanism resolved by
/// [`AcmgConfig::effective_mechanism`] - a user's `gene_overrides` entry if
/// there is one, otherwise the curated `gene_mechanisms` table.
///
/// Blocks only on a *recognised* non-LOF mechanism. An unrecognised string is
/// treated as "not curated" rather than as gain of function, so a typo in a
/// user's TOML cannot silently disable PVS1 for a gene: the failure mode of
/// the gate is to do nothing, never to suppress evidence by accident.
pub(crate) fn lof_mechanism_blocker(
    input: &ClassificationInput,
    config: &AcmgConfig,
) -> Option<String> {
    if !config.mechanism_gates_pvs1 {
        return None;
    }
    let gene = input.gene_symbol.as_deref().unwrap_or("the gene");
    let mechanism = config.effective_mechanism(input.gene_symbol.as_deref())?;

    // "LOF_and_GOF" names a gene with a loss-of-function disease as well as a
    // gain-of-function one - RYR1 being the standing example - and a null
    // allele is still pathogenic for the first. Testing for LOF before GOF is
    // what keeps those genes evaluable.
    let m = mechanism.to_ascii_uppercase();
    if m.contains("LOF") || m.contains("LOSS_OF_FUNCTION") {
        return None;
    }
    let non_lof = ["GOF", "GAIN_OF_FUNCTION", "DOMINANT_NEGATIVE", "DN"];
    if !non_lof.iter().any(|t| m == *t) {
        return None;
    }

    Some(format!(
        "mechanism_not_loss_of_function: {} causes disease by a mechanism that is not loss of function ({}), so a null allele is not the disease mechanism and PVS1 does not apply (Abou Tayoun 2018)",
        gene, mechanism,
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::GeneOverride;
    use crate::sa_extract::OmimData;
    use crate::test_support::minimal_input;
    use std::collections::HashMap;

    fn gene_override(mechanism: &str) -> GeneOverride {
        GeneOverride {
            mechanism: Some(mechanism.to_string()),
            ba1_af_threshold: None,
            bs1_af_threshold: None,
            pm2_af_threshold: None,
            disabled_criteria: vec![],
            strength_overrides: HashMap::new(),
            disorders: HashMap::new(),
        }
    }

    fn with_mechanism(gene: &str, mechanism: &str) -> (ClassificationInput, AcmgConfig) {
        let mut input = minimal_input();
        input.gene_symbol = Some(gene.to_string());
        let mut config = AcmgConfig::default();
        config
            .gene_overrides
            .insert(gene.to_string(), gene_override(mechanism));
        (input, config)
    }

    fn with_phenotypes(entries: &[&str]) -> ClassificationInput {
        let mut input = minimal_input();
        input.omim = Some(OmimData {
            mim_number: Some(0),
            phenotypes: Some(entries.iter().map(|s| s.to_string()).collect()),
        });
        input
    }

    #[test]
    fn test_gene_absent_from_the_source_never_blocks() {
        // The rule that matters most. ClinGen has curated ~2,400 genes, so a
        // gene it has not reached says nothing about the gene. Blocking on
        // absence cost 1,497 truth-pathogenic PVS1 firings in run v10.
        let input = minimal_input();
        assert!(input.omim.is_none());
        assert!(validity_blocker(&input, &AcmgConfig::default()).is_none());
    }

    #[test]
    fn test_established_relationship_does_not_block() {
        for class in ["Definitive", "Strong", "Moderate"] {
            let input = with_phenotypes(&[&format!(
                "a real disease (ClinGen {class}/AD, MONDO:0000002)"
            )]);
            assert!(
                validity_blocker(&input, &AcmgConfig::default()).is_none(),
                "{class} meets the SVI evidence threshold",
            );
        }
    }

    #[test]
    fn test_curated_as_weak_or_refuted_blocks() {
        for class in [
            "Limited",
            "Disputed",
            "Refuted",
            "No Known Disease Relationship",
        ] {
            let input = with_phenotypes(&[&format!(
                "a proposed disease (ClinGen {class}/AD, MONDO:0000001)"
            )]);
            let reason = validity_blocker(&input, &AcmgConfig::default())
                .unwrap_or_else(|| panic!("{class} must block"));
            assert!(
                reason.contains("no_valid_gene_disease_relationship"),
                "got: {reason}"
            );
            assert!(
                reason.contains("TEST"),
                "the reason must name the gene: {reason}"
            );
        }
    }

    #[test]
    fn test_one_established_disease_rescues_a_gene_with_disputed_ones() {
        // Genes routinely carry several proposed relationships of differing
        // quality. One that holds is enough for the criteria to apply.
        let input = with_phenotypes(&[
            "a proposed disease (ClinGen Disputed/AD, MONDO:0000001)",
            "a real disease (ClinGen Definitive/AR, MONDO:0000002)",
        ]);
        assert!(validity_blocker(&input, &AcmgConfig::default()).is_none());
    }

    #[test]
    fn test_omim_entries_without_a_clingen_tag_count_as_established() {
        // The legacy genemap2 source has no classification tier, so its
        // phenotypes must not all read as unclassified-and-therefore-weak.
        let input = with_phenotypes(&["Breast cancer, 114480 (3), Autosomal dominant"]);
        assert!(validity_blocker(&input, &AcmgConfig::default()).is_none());
    }

    #[test]
    fn test_empty_phenotype_list_does_not_block() {
        // A record carrying only a mimNumber is not a curation finding.
        let mut input = minimal_input();
        for phenotypes in [None, Some(vec![])] {
            input.omim = Some(OmimData {
                mim_number: Some(0),
                phenotypes,
            });
            assert!(validity_blocker(&input, &AcmgConfig::default()).is_none());
        }
    }

    #[test]
    fn test_validity_gate_is_switchable() {
        let input = with_phenotypes(&["x (ClinGen Refuted/AD, MONDO:0000001)"]);
        let config = AcmgConfig {
            require_gene_disease_validity: false,
            ..Default::default()
        };
        assert!(validity_blocker(&input, &config).is_none());
    }

    #[test]
    fn test_gof_only_mechanism_blocks_pvs1() {
        let (input, config) = with_mechanism("PCSK9", "GOF");
        let reason = lof_mechanism_blocker(&input, &config).expect("GOF-only must block");
        assert!(
            reason.contains("mechanism_not_loss_of_function"),
            "got: {reason}"
        );
        assert!(reason.contains("PCSK9"));
    }

    #[test]
    fn test_lof_and_gof_gene_still_has_a_lof_disease() {
        // RYR1 is the case: malignant hyperthermia is gain of function, the
        // congenital myopathies are loss of function. A null allele is still
        // pathogenic for one of the two, so PVS1 must survive.
        for mechanism in ["LOF_and_GOF", "GOF_and_LOF", "lof_and_gof"] {
            let (input, config) = with_mechanism("RYR1", mechanism);
            assert!(
                lof_mechanism_blocker(&input, &config).is_none(),
                "{mechanism} must not block PVS1",
            );
        }
    }

    #[test]
    fn test_unknown_mechanism_does_not_block() {
        // The common case by far. Absence of a curated mechanism is not
        // evidence of gain of function.
        let input = minimal_input();
        assert!(lof_mechanism_blocker(&input, &AcmgConfig::default()).is_none());
        let (input, config) = with_mechanism("BRCA1", "LOF");
        assert!(lof_mechanism_blocker(&input, &config).is_none());
    }

    #[test]
    fn test_unrecognised_mechanism_string_does_not_block() {
        // A misspelling must fail safe. Blocking on anything that merely fails
        // to say "LOF" would let a typo silently delete PVS1 for a gene.
        for mechanism in ["gian_of_function", "unknown", "?", ""] {
            let (input, config) = with_mechanism("SOMEGENE", mechanism);
            assert!(
                lof_mechanism_blocker(&input, &config).is_none(),
                "{mechanism:?} is not a recognised non-LOF mechanism and must not block",
            );
        }
    }

    #[test]
    fn test_dominant_negative_blocks_like_gain_of_function() {
        // The question PVS1 asks is whether LOF is the mechanism, not whether
        // GOF is. MYH7 is dominant-negative and the ClinGen specification says
        // PVS1 does not apply.
        let (input, config) = with_mechanism("MYH7", "DOMINANT_NEGATIVE");
        assert!(lof_mechanism_blocker(&input, &config).is_some());
    }

    #[test]
    fn test_curated_table_applies_without_a_user_override() {
        // The shipped table has to work on its own; requiring a TOML entry
        // would leave the default config with a gate that never fires.
        let mut input = minimal_input();
        input.gene_symbol = Some("PCSK9".to_string());
        let reason = lof_mechanism_blocker(&input, &AcmgConfig::default())
            .expect("the curated GOF table must apply by default");
        assert!(reason.contains("PCSK9"), "got: {reason}");

        input.gene_symbol = Some("RYR1".to_string());
        assert!(
            lof_mechanism_blocker(&input, &AcmgConfig::default()).is_none(),
            "RYR1 is curated LOF_and_GOF and must stay evaluable",
        );
    }

    #[test]
    fn test_user_override_beats_the_curated_table() {
        let (input, config) = with_mechanism("PCSK9", "LOF");
        assert!(
            lof_mechanism_blocker(&input, &config).is_none(),
            "an explicit gene_overrides mechanism must win over the shipped table",
        );
    }

    #[test]
    fn test_setting_one_gene_override_does_not_wipe_the_curated_table() {
        // The reason the curated mechanisms live in their own map: a TOML that
        // sets a single [gene_overrides.X] block replaces `gene_overrides`
        // wholesale, and would take the shipped table with it.
        let config: AcmgConfig = toml::from_str("[gene_overrides.BRCA1]\nmechanism = \"LOF\"\n")
            .expect("config must parse");
        let mut input = minimal_input();
        input.gene_symbol = Some("PCSK9".to_string());
        assert!(lof_mechanism_blocker(&input, &config).is_some());
    }

    #[test]
    fn test_mechanism_gate_is_switchable() {
        let (input, mut config) = with_mechanism("PCSK9", "GOF");
        config.mechanism_gates_pvs1 = false;
        assert!(lof_mechanism_blocker(&input, &config).is_none());
    }
}
