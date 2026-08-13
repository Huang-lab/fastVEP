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

/// Why a frequency criterion could not be assessed, and whether the reason is a
/// statement about *this site's short-read data* or about something else.
///
/// The distinction is what [`crate::combiner`] needs to decide whether a
/// pathogenic call resting on read-derived evidence is still safe. A gnomAD
/// FILTER, a low-complexity tract and a segmental duplication all say the same
/// thing: alignments here are unreliable. That undermines the frameshift or
/// canonical-splice call PVS1 rests on just as much as it undermines the allele
/// frequency, because both are read from the same pileup. A curated
/// homology-gene entry or a ClinVar low-penetrance label says something
/// narrower - that *the frequency* is confounded - and leaves the consequence
/// call alone.
pub(crate) struct FrequencyBlocker {
    pub reason: String,
    /// `true` for the read-level vetoes (gnomAD FILTER, low-complexity or
    /// segmental-duplication flags); `false` for gene-level homology curation
    /// and the ClinVar low-penetrance precondition.
    pub site_level: bool,
}

impl FrequencyBlocker {
    fn gene_level(reason: String) -> Self {
        Self { reason, site_level: false }
    }

    fn site_level(reason: String) -> Self {
        Self { reason, site_level: true }
    }

    /// Record this blocker in a criterion's `details`, under the two keys the
    /// combiner and the review tooling read.
    pub(crate) fn record(&self, details: &mut serde_json::Map<String, serde_json::Value>) {
        details.insert("frequency_blocked".into(), serde_json::json!(self.reason));
        details.insert(
            FREQUENCY_BLOCKED_SITE_LEVEL.into(),
            serde_json::json!(self.site_level),
        );
    }
}

/// Detail key carrying [`FrequencyBlocker::site_level`]. Named here so the
/// combiner reads the same string the criteria write.
pub(crate) const FREQUENCY_BLOCKED_SITE_LEVEL: &str = "frequency_blocked_site_level";

/// Detail key set by [`note_withheld_benign_frequency`].
pub(crate) const FREQUENCY_BLOCKED_WOULD_BE_BENIGN: &str = "frequency_blocked_would_be_benign";

/// Record whether the frequency the classifier just refused to believe would
/// have carried benign evidence had it been believed.
///
/// Refusing to read a frequency is only *materially* one-sided when there was
/// something in it to read. A site that fails gnomAD's FILTER with two carriers
/// costs the benign side nothing when it is withheld; the same veto over an
/// allele seen in a third of the population withholds BA1. The combiner needs
/// to tell those apart, because demoting every call at every flagged site would
/// cost thousands of correct pathogenic calls to avoid a few dozen wrong ones -
/// measured on the ClinVar 2-star+ set, 2,630 against 18.
pub(crate) fn note_withheld_benign_frequency(
    input: &ClassificationInput,
    config: &AcmgConfig,
    threshold: f64,
    details: &mut serde_json::Map<String, serde_json::Value>,
) {
    let Some((af, _)) = input.gnomad.as_ref().and_then(|g| benign_test_af(g, config)) else {
        return;
    };
    if af > threshold {
        details.insert(FREQUENCY_BLOCKED_WOULD_BE_BENIGN.into(), serde_json::json!(true));
        details.insert("frequency_blocked_af".into(), serde_json::json!(af));
    }
}

/// The curated exception forbidding `code` on this variant, if any.
///
/// Coordinate first, HGVS second - see [`crate::config::Ba1Exception::matches`]
/// for why that order is load-bearing rather than a preference.
pub(crate) fn curated_exception<'a>(
    input: &ClassificationInput,
    config: &'a AcmgConfig,
    code: &str,
) -> Option<&'a crate::config::Ba1Exception> {
    let coordinates = input
        .variant_coordinates
        .as_ref()
        .map(|(chrom, pos, r, a)| (chrom.as_str(), *pos, r.as_str(), a.as_str()));
    config.frequency_exception(
        code,
        input.gene_symbol.as_deref(),
        input.hgvs_c.as_deref(),
        coordinates,
    )
}

/// Record a curated exception in a criterion's `details` and return the
/// summary line for it.
pub(crate) fn record_exception(
    code: &str,
    exc: &crate::config::Ba1Exception,
    details: &mut serde_json::Map<String, serde_json::Value>,
) -> String {
    details.insert("exception_match".into(), serde_json::json!(true));
    details.insert("exception_gene".into(), serde_json::json!(exc.gene));
    details.insert("exception_hgvs_c".into(), serde_json::json!(exc.hgvs_c));
    if let Some(reason) = &exc.reason {
        details.insert("exception_reason".into(), serde_json::json!(reason));
    }
    format!(
        "{} {} is on the curated frequency-exception list, so {} cannot fire ({})",
        exc.gene,
        exc.hgvs_c,
        code,
        exc.reason.as_deref().unwrap_or("no reason recorded"),
    )
}

/// A reason the population-frequency record for this variant cannot support
/// any inference, in either direction.
pub(crate) fn data_blocker(
    input: &ClassificationInput,
    config: &AcmgConfig,
) -> Option<FrequencyBlocker> {
    // Curated gene-level homology. Reads from a paralogue or pseudogene mismap
    // onto the gene of interest, so both the frequency and the absence of a
    // frequency are artefacts (Mandelker 2016, PMID 27228465).
    if config.is_homology_unreliable(input.gene_symbol.as_deref()) {
        return Some(FrequencyBlocker::gene_level(format!(
            "gene {} has paralogue/pseudogene homology that makes population frequencies unreliable (Mandelker 2016, PMID 27228465)",
            input.gene_symbol.as_deref().unwrap_or("?")
        )));
    }

    let gnomad = input.gnomad.as_ref()?;

    // gnomAD's own QC verdict. A non-PASS record is not a measurement of a
    // population frequency: AC0 means every genotype behind the call was
    // filtered out, AS_VQSR that the site failed the variant-quality model,
    // and InbreedingCoeff that the genotypes are distributed unlike real
    // population data.
    if let Some(filter) = gnomad.failed_filter() {
        return Some(FrequencyBlocker::site_level(format!(
            "gnomAD FILTER={}; the record failed gnomAD's own quality control and is not a measurement of a population frequency",
            filter
        )));
    }

    // Per-site homology, the same concern as the curated gene list but
    // resolved at the site rather than the gene. Gated separately because it
    // is the more aggressive of the two: it fires on individual sites inside
    // otherwise well-behaved genes.
    if config.gnomad_region_flags_block_frequency {
        if let Some(region) = gnomad.unreliable_region() {
            return Some(FrequencyBlocker::site_level(format!(
                "gnomAD flags this site as falling in a {}, where short-read allele frequencies are systematically unreliable",
                region
            )));
        }
    }

    None
}

/// A reason frequency-based *benign* evidence cannot be assessed. Adds the
/// benign-only preconditions on top of [`data_blocker`].
pub(crate) fn benign_blocker(
    input: &ClassificationInput,
    config: &AcmgConfig,
) -> Option<FrequencyBlocker> {
    if let Some(blocker) = data_blocker(input, config) {
        return Some(blocker);
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
                    return Some(FrequencyBlocker::gene_level(format!(
                        "ClinVar ({} stars) reports this variant as \"{}\"; a low-penetrance or risk allele is expected to be frequent and is outside BS2's full-penetrance precondition",
                        cv.review_stars(),
                        term
                    )));
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
            let blocker = data_blocker(&input, &AcmgConfig::default())
                .expect("a non-PASS gnomAD record must block");
            assert!(blocker.reason.contains("FILTER="));
            assert!(
                blocker.site_level,
                "gnomAD's own FILTER is a verdict on this site's reads"
            );
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
            let blocker = data_blocker(&input, &AcmgConfig::default()).expect("should block");
            assert!(blocker.reason.contains(expected), "got: {}", blocker.reason);
            assert!(blocker.site_level, "a region flag is a verdict on this site's reads");
        }
    }

    #[test]
    fn test_gene_level_vetoes_are_not_site_level() {
        // Both of these say the *frequency* is confounded. Neither says the
        // alignment at this position is wrong, so neither may be used to
        // discount a consequence call.
        let mut input = minimal_input();
        input.gene_symbol = Some("PMS2".to_string());
        let config = AcmgConfig::default();
        if let Some(blocker) = data_blocker(&input, &config) {
            assert!(!blocker.site_level, "curated homology is a gene-level claim");
        }

        let mut input = input_with(GnomadData { all_af: Some(0.01), ..Default::default() });
        input.clinvar = Some(crate::sa_extract::ClinvarData {
            significance: Some(vec!["Pathogenic,_low_penetrance".to_string()]),
            review_status: Some("criteria_provided,_multiple_submitters,_no_conflicts".to_string()),
            ..Default::default()
        });
        let blocker = benign_blocker(&input, &config).expect("low penetrance must block");
        assert!(
            !blocker.site_level,
            "a low-penetrance label says the frequency is expected, not that the reads are wrong"
        );
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
