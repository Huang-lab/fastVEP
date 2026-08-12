//! Preconditions shared by every criterion that reasons from population
//! allele frequency: BA1, BS1, BS2 and PM2.
//!
//! All four ask a question of the same underlying gnomAD record, so they have
//! to agree about when that record can be believed. Keeping the checks here,
//! rather than repeating them per criterion, is what stops a gene being
//! trusted for one frequency criterion and distrusted for another.
//!
//! Two layers:
//!
//! * [`data_blocker`] - reasons the population data itself is not usable.
//!   Direction-neutral, so it gates the pathogenic-direction PM2 ("absent from
//!   controls") exactly as it gates the benign-direction BA1/BS1/BS2. When the
//!   frequency cannot be believed, neither its presence nor its absence is
//!   evidence.
//! * [`benign_blocker`] - `data_blocker` plus the reasons that apply only to
//!   benign-direction criteria, where a high frequency is expected and so is
//!   not evidence against pathogenicity.
//!
//! Every check degrades to "no blocker" when the annotation database predates
//! the field it reads, so an older database behaves exactly as it did before
//! these gates existed.

use crate::config::AcmgConfig;
use crate::sa_extract::ClassificationInput;

/// A reason the population-frequency record for this variant cannot support
/// any inference, in either direction.
pub(crate) fn data_blocker(
    input: &ClassificationInput,
    config: &AcmgConfig,
) -> Option<String> {
    // Curated gene-level homology. Reads from a paralogue or pseudogene mismap
    // onto the gene of interest, so both the frequency and the absence of a
    // frequency are artefacts (Mandelker 2016, PMID 27228465).
    if config.is_homology_unreliable(input.gene_symbol.as_deref()) {
        return Some(format!(
            "gene {} has paralogue/pseudogene homology that makes population frequencies unreliable (Mandelker 2016, PMID 27228465)",
            input.gene_symbol.as_deref().unwrap_or("?")
        ));
    }

    let gnomad = input.gnomad.as_ref()?;

    // gnomAD's own QC verdict. A non-PASS record is not a measurement of a
    // population frequency: AC0 means every genotype behind the call was
    // filtered out, AS_VQSR that the site failed the variant-quality model,
    // and InbreedingCoeff that the genotypes are distributed unlike real
    // population data.
    if let Some(filter) = gnomad.failed_filter() {
        return Some(format!(
            "gnomAD FILTER={}; the record failed gnomAD's own quality control and is not a measurement of a population frequency",
            filter
        ));
    }

    // Per-site homology, the same concern as the curated gene list but
    // resolved at the site rather than the gene. Gated separately because it
    // is the more aggressive of the two: it fires on individual sites inside
    // otherwise well-behaved genes.
    if config.gnomad_region_flags_block_frequency {
        if let Some(region) = gnomad.unreliable_region() {
            return Some(format!(
                "gnomAD flags this site as falling in a {}, where short-read allele frequencies are systematically unreliable",
                region
            ));
        }
    }

    None
}

/// A reason frequency-based *benign* evidence cannot be assessed. Adds the
/// benign-only preconditions on top of [`data_blocker`].
pub(crate) fn benign_blocker(
    input: &ClassificationInput,
    config: &AcmgConfig,
) -> Option<String> {
    if let Some(reason) = data_blocker(input, config) {
        return Some(reason);
    }

    // Richards 2015 conditions BS1/BS2 on a disorder with full penetrance
    // expected at an early age. ClinVar's controlled vocabulary marks the
    // variants that are known exceptions, and a variant that is *expected* to
    // be frequent cannot have its frequency read as evidence against
    // pathogenicity.
    if config.clinvar_low_penetrance_blocks_benign_frequency {
        if let Some(ref cv) = input.clinvar {
            if cv.review_stars() >= 2 {
                if let Some(term) = cv.low_penetrance_term() {
                    return Some(format!(
                        "ClinVar ({} stars) reports this variant as \"{}\"; a low-penetrance or risk allele is expected to be frequent and is outside BS2's full-penetrance precondition",
                        cv.review_stars(),
                        term
                    ));
                }
            }
        }
    }

    None
}

/// The allele frequency that BA1 and BS1 should test against, with a label
/// naming the statistic so the criterion summary can report which was used.
///
/// Prefers gnomAD's filtering allele frequency - the lower bound of the 95 %
/// confidence interval on the frequency, maximised over genetic-ancestry
/// groups (Whiffin 2017, PMID 28518168). Testing a point estimate instead
/// treats "AF 0.06 measured from 200 alleles in one small group" as though it
/// were as solid as "AF 0.06 measured from 700,000", and the first is the
/// shape a spuriously high subpopulation frequency usually takes.
///
/// Note what this does and does not address in the round-2 review's founder
/// point. It removes the small-sample half: a frequency inflated by thin
/// sampling in one group no longer reaches BA1/BS1. It does not remove the
/// other half, a pathogenic founder allele that really is common in a
/// well-sampled population - its FAF is high too, and only the curated BA1
/// exception list (Ghosh 2018) keeps BA1 off those.
///
/// Falls back to the population-maximum point estimate when the database
/// predates the FAF columns, so behaviour is unchanged against an old
/// database.
pub(crate) fn benign_test_af(
    gnomad: &crate::sa_extract::GnomadData,
    config: &AcmgConfig,
) -> Option<(f64, &'static str)> {
    if config.use_filtering_af {
        if let Some(faf) = gnomad.filtering_af() {
            return Some((faf, "filtering AF (95% CI lower bound, max across ancestry groups)"));
        }
    }
    gnomad.max_pop_af().map(|af| (af, "max population AF"))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sa_extract::GnomadData;
    use crate::test_support::minimal_input;

    fn input_with(gnomad: GnomadData) -> ClassificationInput {
        let mut input = minimal_input();
        input.gnomad = Some(gnomad);
        input
    }

    #[test]
    fn test_pass_record_is_not_blocked() {
        let input = input_with(GnomadData {
            all_af: Some(0.01),
            ..Default::default()
        });
        assert!(data_blocker(&input, &AcmgConfig::default()).is_none());
    }

    #[test]
    fn test_failed_filter_blocks_in_both_directions() {
        for gnomad in [
            GnomadData { filter_ac0: true, ..Default::default() },
            GnomadData { filter_vqsr: true, ..Default::default() },
            GnomadData { filter_inbreeding: true, ..Default::default() },
        ] {
            let input = input_with(gnomad);
            let reason = data_blocker(&input, &AcmgConfig::default());
            assert!(reason.is_some(), "a non-PASS gnomAD record must block");
            assert!(reason.unwrap().contains("FILTER="));
        }
    }

    #[test]
    fn test_segdup_and_lcr_block_when_enabled() {
        for (gnomad, expected) in [
            (
                GnomadData { segdup: true, ..Default::default() },
                "segmental duplication",
            ),
            (
                GnomadData { lcr: true, ..Default::default() },
                "low-complexity region",
            ),
        ] {
            let input = input_with(gnomad);
            let reason = data_blocker(&input, &AcmgConfig::default()).expect("should block");
            assert!(reason.contains(expected), "got: {reason}");
        }
    }

    #[test]
    fn test_region_flags_can_be_disabled() {
        let input = input_with(GnomadData { segdup: true, ..Default::default() });
        let config = AcmgConfig {
            gnomad_region_flags_block_frequency: false,
            ..Default::default()
        };
        assert!(data_blocker(&input, &config).is_none());
        // The FILTER gate is not covered by that switch: a record gnomAD
        // itself rejected is never usable.
        let input = input_with(GnomadData { filter_ac0: true, ..Default::default() });
        assert!(data_blocker(&input, &config).is_some());
    }

    #[test]
    fn test_database_without_the_new_fields_is_not_blocked() {
        // Everything a pre-v4-QC database can say about a clean common variant.
        // The gates must be invisible to it, or upgrading fastVEP without
        // rebuilding the annotation database would silently change calls.
        let input = input_with(GnomadData {
            all_af: Some(0.2),
            all_an: Some(1_460_000),
            all_ac: Some(292_000),
            all_hc: Some(30_000),
            ..Default::default()
        });
        assert!(data_blocker(&input, &AcmgConfig::default()).is_none());
        assert!(benign_blocker(&input, &AcmgConfig::default()).is_none());
    }

    #[test]
    fn test_benign_blocker_adds_low_penetrance_on_top_of_data_blocker() {
        use crate::sa_extract::ClinvarData;
        let mut input = input_with(GnomadData {
            all_af: Some(0.037),
            ..Default::default()
        });
        input.clinvar = Some(ClinvarData {
            significance: Some(vec!["Pathogenic,_low_penetrance".into()]),
            review_status: Some("criteria_provided,_multiple_submitters,_no_conflicts".into()),
            ..Default::default()
        });
        // Direction-neutral gate says nothing; the benign-only gate blocks.
        assert!(data_blocker(&input, &AcmgConfig::default()).is_none());
        assert!(benign_blocker(&input, &AcmgConfig::default()).is_some());
    }
}
