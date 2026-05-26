//! LOFTEE (Loss-Of-Function Transcript Effect Estimator) for fastVEP.
//!
//! Classifies putative loss-of-function variants (stop-gained, frameshift,
//! splice-site disrupting) as High-Confidence (HC), Low-Confidence (LC),
//! or Other-Splice (OS) based on LOFTEE filtering criteria.

mod filters;
pub mod conservation;
pub mod maxentscan;
pub mod motifs;
pub mod splice;
pub mod svm;

pub use filters::*;

use fastvep_cache::providers::SequenceProvider;
use fastvep_core::Consequence;
use fastvep_genome::Transcript;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// LOFTEE confidence levels.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum Confidence {
    /// High Confidence: canonical LoF with no failing filters.
    HC,
    /// Low Confidence: canonical LoF that failed one or more filters.
    LC,
    /// Other Splice: LOFTEE-predicted splice disruption not called by VEP.
    OS,
}

impl std::fmt::Display for Confidence {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Confidence::HC => write!(f, "HC"),
            Confidence::LC => write!(f, "LC"),
            Confidence::OS => write!(f, "OS"),
        }
    }
}

/// Result of LOFTEE evaluation for a single allele annotation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LofteeResult {
    pub confidence: String,
    pub filters: Vec<String>,
    pub flags: Vec<String>,
    pub info: HashMap<String, String>,
}

/// Configuration for LOFTEE evaluation.
pub struct LofteeConfig {
    pub min_intron_size: u64,
    pub ancestral_seq_provider: Option<Box<dyn SequenceProvider + Send + Sync>>,
    /// Splice prediction data (Phase 3). None = splice predictions disabled.
    pub splice_data: Option<splice::SpliceData>,
    /// GERP conservation provider (Phase 4). None = use 50bp fallback for END_TRUNC.
    pub gerp_provider: Option<Box<dyn conservation::ConservationProvider>>,
    /// PhyloCSF provider (Phase 4). None = skip PhyloCSF flags.
    pub phylocsf_provider: Option<Box<dyn conservation::PhyloCsfProvider>>,
}

impl Default for LofteeConfig {
    fn default() -> Self {
        Self {
            min_intron_size: 15,
            ancestral_seq_provider: None,
            splice_data: None,
            gerp_provider: None,
            phylocsf_provider: None,
        }
    }
}

/// Input data for LOFTEE evaluation.
pub struct LofteeInput<'a> {
    pub transcript: &'a Transcript,
    pub consequences: &'a [Consequence],
    pub exon: Option<(u32, u32)>,
    pub intron: Option<(u32, u32)>,
    pub cds_end: Option<u64>,
    pub variant_start: u64,
    pub variant_end: u64,
    pub allele: &'a str,
    pub ref_allele: &'a str,
    pub flags: &'a [String],
}

/// Check if the consequence set contains a putative LoF consequence.
pub fn has_lof_consequence(consequences: &[Consequence]) -> bool {
    consequences.iter().any(|c| matches!(c,
        Consequence::StopGained
        | Consequence::FrameshiftVariant
        | Consequence::SpliceAcceptorVariant
        | Consequence::SpliceDonorVariant
    ))
}

/// Evaluate LOFTEE filters for a single allele annotation.
pub fn evaluate(
    input: &LofteeInput,
    config: &LofteeConfig,
    seq_provider: Option<&(dyn SequenceProvider + Send + Sync)>,
) -> LofteeResult {
    let mut filters: Vec<String> = Vec::new();
    let mut flags: Vec<String> = Vec::new();
    let mut info: HashMap<String, String> = HashMap::new();

    let consequences = input.consequences;
    let tr = input.transcript;

    let is_stop_gained = consequences.contains(&Consequence::StopGained);
    let is_frameshift = consequences.contains(&Consequence::FrameshiftVariant);
    let is_splice_donor = consequences.contains(&Consequence::SpliceDonorVariant);
    let is_splice_acceptor = consequences.contains(&Consequence::SpliceAcceptorVariant);
    let is_vep_splice_lof = is_splice_donor || is_splice_acceptor;
    let is_other_lof = is_stop_gained || is_frameshift;

    let is_upstream = consequences.contains(&Consequence::UpstreamGeneVariant);
    let is_downstream = consequences.contains(&Consequence::DownstreamGeneVariant);
    let is_genic = !is_upstream && !is_downstream;
    let is_utr = filters::check_5utr_splice(tr, input.variant_start, input.variant_end)
        || filters::check_3utr_splice(tr, input.variant_start, input.variant_end);

    // Splice predictions (Phase 3): run on genic, non-UTR, non-stop/frameshift variants
    let mut loftee_splice_lof = false;
    let mut _lof_position: Option<u64> = None;

    if let Some(ref splice_data) = config.splice_data {
        if is_genic && !is_utr && !is_other_lof {
            // Extended splice: check if variant disrupts an annotated splice site
            if let Some(intron) = input.intron {
                let splice_result = splice::get_effect_on_splice(
                    tr,
                    intron.0 as usize,
                    input.variant_start,
                    input.variant_end,
                    input.ref_allele,
                    input.allele,
                    is_vep_splice_lof,
                    splice_data,
                    seq_provider,
                );

                if let Some(result) = splice_result {
                    for (k, v) in &result.features {
                        info.insert(k.clone(), v.clone());
                    }

                    if result.disrupting {
                        loftee_splice_lof = true;
                        info.insert(
                            format!("{}_DISRUPTING", result.splice_type),
                            "true".into(),
                        );
                    } else if is_vep_splice_lof {
                        filters.push(format!("NON_{}_DISRUPTING", result.splice_type));
                    }
                }
            }

            // De novo donor: check if variant creates a new strong donor
            if !loftee_splice_lof {
                if let Some(exon) = input.exon {
                    let denovo = splice::check_for_denovo_donor(
                        tr,
                        exon.0 as usize,
                        input.variant_start,
                        input.variant_end,
                        input.ref_allele,
                        input.allele,
                        splice_data,
                        seq_provider,
                    );

                    if let Some(result) = denovo {
                        if result.probability > splice_data.denovo_donor_cutoff {
                            info.insert("DE_NOVO_DONOR".into(), "true".into());
                            _lof_position = Some(result.lof_pos);
                            loftee_splice_lof = result.lof;
                        }
                        for (k, v) in &result.features {
                            info.insert(k.clone(), v.clone());
                        }
                    }
                }
            }
        }
    }

    // Determine initial confidence
    let mut confidence = if is_vep_splice_lof || is_other_lof {
        Confidence::HC
    } else if loftee_splice_lof {
        Confidence::OS
    } else {
        return LofteeResult {
            confidence: String::new(),
            filters: vec![],
            flags: vec![],
            info,
        };
    };

    // === END_TRUNC filter (stop-gained / frameshift in exonic context) ===
    if is_other_lof {
        if let Some(cds_end) = input.cds_end {
            let cds_length = filters::get_cds_length_fast(tr);
            if cds_length > 0 {
                let percentile = cds_end as f64 / cds_length as f64;
                info.insert("PERCENTILE".into(), format!("{:.4}", percentile));
            }

            if let Some((_, _total_exons)) = input.exon {
                let last_exon_coding_len =
                    filters::get_last_exon_coding_length(tr);
                let bp_dist = filters::get_bp_dist_to_stop(tr, input.variant_start);
                let d = bp_dist as i64 - last_exon_coding_len as i64;
                info.insert("BP_DIST".into(), bp_dist.to_string());
                info.insert("DIST_FROM_LAST_EXON".into(), d.to_string());
                info.insert(
                    "50_BP_RULE".into(),
                    if d <= 50 { "FAIL" } else { "PASS" }.into(),
                );

                // END_TRUNC: if GERP available, require d <= 50 AND gerp_dist <= 180
                // Otherwise fall back to just d <= 50 (conservative over-filtering)
                if d <= 50 {
                    if let Some(ref gerp) = config.gerp_provider {
                        let (gerp_dist, _) = conservation::get_gerp_weighted_dist(
                            tr,
                            input.variant_start,
                            gerp.as_ref(),
                        );
                        info.insert("GERP_DIST".into(), format!("{:.2}", gerp_dist));
                        if gerp_dist <= 180.0 {
                            filters.push("END_TRUNC".into());
                        }
                    } else {
                        // No GERP data: conservative fallback
                        filters.push("END_TRUNC".into());
                    }
                }
            } else {
                flags.push("NO_EXON_NUMBER".into());
            }
        }
    }

    // === Exonic filters (stop-gained / frameshift) ===
    if is_other_lof && input.exon.is_some() {
        let (exon_num, total_exons) = input.exon.unwrap();

        if filters::check_exon_intron_undef(input.exon) {
            filters.push("EXON_INTRON_UNDEF".into());
        } else if filters::check_single_exon(total_exons) {
            flags.push("SINGLE_EXON".into());
        } else {
            if filters::check_incomplete_cds(&input.flags) {
                filters.push("INCOMPLETE_CDS".into());
            }
        }

        // PhyloCSF conservation check (Phase 4)
        if let Some(ref phylocsf) = config.phylocsf_provider {
            if let Some((flag, phylo_info)) =
                conservation::check_phylocsf(tr, exon_num, phylocsf.as_ref())
            {
                flags.push(flag);
                for (k, v) in phylo_info {
                    info.insert(k, v);
                }
            } else {
                info.insert("PHYLOCSF_TOO_SHORT".into(), "true".into());
            }
        }
    }

    // === Intronic filters (splice variants) ===
    if input.intron.is_some() {
        let (intron_idx, _total_introns) = input.intron.unwrap();

        if filters::check_intron_annotation_errors(input.intron) {
            filters.push("EXON_INTRON_UNDEF".into());
        } else {
            // Intron size filter
            let intron_size = filters::get_intron_size(tr, intron_idx as usize);
            info.insert("INTRON_SIZE".into(), intron_size.to_string());
            if intron_size < config.min_intron_size {
                filters.push("SMALL_INTRON".into());
            }

            if is_vep_splice_lof {
                // GC_TO_GT_DONOR: splice donor GC→GT mutation
                if is_splice_donor {
                    if filters::check_gc_to_gt_donor(
                        tr,
                        intron_idx as usize,
                        input.ref_allele,
                        input.allele,
                        seq_provider,
                    ) {
                        filters.push("GC_TO_GT_DONOR".into());
                    }
                }

                // NON_CAN_SPLICE: non-canonical splice site motif
                if filters::check_non_canonical_splice(
                    tr,
                    intron_idx as usize,
                    seq_provider,
                ) {
                    flags.push("NON_CAN_SPLICE".into());
                }

                // UTR splice filters
                if filters::check_5utr_splice(
                    tr,
                    input.variant_start,
                    input.variant_end,
                ) {
                    filters.push("5UTR_SPLICE".into());
                }
                if filters::check_3utr_splice(
                    tr,
                    input.variant_start,
                    input.variant_end,
                ) {
                    filters.push("3UTR_SPLICE".into());
                }
            }

            // NAGNAG_SITE: splice acceptor AG.AG motif
            if is_splice_acceptor {
                if filters::check_nagnag(tr, intron_idx as usize, seq_provider) {
                    flags.push("NAGNAG_SITE".into());
                }
            }
        }
    }

    // === ANC_ALLELE filter (Phase 2) ===
    if let Some(ref anc_provider) = config.ancestral_seq_provider {
        if filters::check_ancestral_allele(
            &tr.chromosome,
            input.variant_start,
            input.variant_end,
            input.allele,
            input.ref_allele,
            anc_provider.as_ref(),
        ) {
            filters.push("ANC_ALLELE".into());
        }
    }

    // === Confidence downgrade ===
    // If there are any filters, downgrade from HC → LC (or OS → LC)
    if !filters.is_empty() {
        confidence = Confidence::LC;
    }

    LofteeResult {
        confidence: confidence.to_string(),
        filters,
        flags,
        info,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fastvep_core::Strand;
    use fastvep_genome::{Exon, Gene, Transcript, Translation};

    pub(crate) fn make_multi_exon_transcript() -> Transcript {
        Transcript {
            stable_id: "ENST00000001".into(),
            version: None,
            gene: Gene {
                stable_id: "ENSG00000001".into(),
                symbol: Some("TEST".into()),
                symbol_source: None,
                hgnc_id: None,
                biotype: "protein_coding".into(),
                chromosome: "chr1".into(),
                start: 1000,
                end: 5000,
                strand: Strand::Forward,
            },
            biotype: "protein_coding".into(),
            chromosome: "chr1".into(),
            start: 1000,
            end: 5000,
            strand: Strand::Forward,
            exons: vec![
                Exon {
                    stable_id: "E1".into(),
                    start: 1000,
                    end: 1200,
                    strand: Strand::Forward,
                    phase: 0,
                    end_phase: 0,
                    rank: 1,
                },
                Exon {
                    stable_id: "E2".into(),
                    start: 2000,
                    end: 2300,
                    strand: Strand::Forward,
                    phase: 0,
                    end_phase: 0,
                    rank: 2,
                },
                Exon {
                    stable_id: "E3".into(),
                    start: 4000,
                    end: 5000,
                    strand: Strand::Forward,
                    phase: 0,
                    end_phase: 0,
                    rank: 3,
                },
            ],
            translation: Some(Translation {
                stable_id: "ENSP1".into(),
                genomic_start: 1050,
                genomic_end: 4500,
                start_exon_rank: 1,
                start_exon_offset: 50,
                end_exon_rank: 3,
                end_exon_offset: 500,
            }),
            cdna_coding_start: Some(51),
            cdna_coding_end: Some(952),
            coding_region_start: Some(1050),
            coding_region_end: Some(4500),
            spliced_seq: None,
            translateable_seq: None,
            peptide: None,
            canonical: true,
            mane_select: None,
            mane_plus_clinical: None,
            tsl: None,
            appris: None,
            ccds: None,
            protein_id: None,
            protein_version: None,
            swissprot: vec![],
            trembl: vec![],
            uniparc: vec![],
            refseq_id: None,
            source: None,
            gencode_primary: false,
            flags: vec![],
            codon_table_start_phase: 0,
        }
    }

    fn make_single_exon_transcript() -> Transcript {
        let mut tr = make_multi_exon_transcript();
        tr.exons = vec![Exon {
            stable_id: "E1".into(),
            start: 1000,
            end: 2000,
            strand: Strand::Forward,
            phase: 0,
            end_phase: 0,
            rank: 1,
        }];
        tr
    }

    #[test]
    fn test_has_lof_consequence() {
        assert!(has_lof_consequence(&[Consequence::StopGained]));
        assert!(has_lof_consequence(&[Consequence::FrameshiftVariant]));
        assert!(has_lof_consequence(&[Consequence::SpliceAcceptorVariant]));
        assert!(has_lof_consequence(&[Consequence::SpliceDonorVariant]));
        assert!(!has_lof_consequence(&[Consequence::MissenseVariant]));
        assert!(!has_lof_consequence(&[Consequence::SynonymousVariant]));
        assert!(!has_lof_consequence(&[]));
    }

    #[test]
    fn test_evaluate_stop_gained_hc() {
        let tr = make_multi_exon_transcript();
        let config = LofteeConfig::default();
        let input = LofteeInput {
            transcript: &tr,
            consequences: &[Consequence::StopGained],
            exon: Some((1, 3)),
            intron: None,
            cds_end: Some(300), // early in CDS, should pass END_TRUNC
            variant_start: 2100,
            variant_end: 2100,
            allele: "T",
            ref_allele: "C",
            flags: &[],
        };
        let result = evaluate(&input, &config, None);
        assert_eq!(result.confidence, "HC");
        assert!(result.filters.is_empty());
    }

    #[test]
    fn test_evaluate_stop_gained_end_trunc() {
        let tr = make_multi_exon_transcript();
        let config = LofteeConfig::default();
        // Variant near the end of CDS in the last exon
        let input = LofteeInput {
            transcript: &tr,
            consequences: &[Consequence::StopGained],
            exon: Some((2, 3)),
            intron: None,
            cds_end: Some(890), // near end of CDS
            variant_start: 4480, // close to coding_region_end (4500)
            variant_end: 4480,
            allele: "T",
            ref_allele: "C",
            flags: &[],
        };
        let result = evaluate(&input, &config, None);
        assert_eq!(result.confidence, "LC");
        assert!(result.filters.contains(&"END_TRUNC".into()));
    }

    #[test]
    fn test_evaluate_single_exon_flag() {
        let tr = make_single_exon_transcript();
        let config = LofteeConfig::default();
        let input = LofteeInput {
            transcript: &tr,
            consequences: &[Consequence::StopGained],
            exon: Some((0, 1)),
            intron: None,
            cds_end: Some(300),
            variant_start: 1300,
            variant_end: 1300,
            allele: "T",
            ref_allele: "C",
            flags: &[],
        };
        let result = evaluate(&input, &config, None);
        // SINGLE_EXON is a flag, not a filter, so confidence stays HC
        assert_eq!(result.confidence, "HC");
        assert!(result.flags.contains(&"SINGLE_EXON".into()));
    }

    #[test]
    fn test_evaluate_incomplete_cds() {
        let tr = make_multi_exon_transcript();
        let config = LofteeConfig::default();
        let input = LofteeInput {
            transcript: &tr,
            consequences: &[Consequence::StopGained],
            exon: Some((1, 3)),
            intron: None,
            cds_end: Some(300),
            variant_start: 2100,
            variant_end: 2100,
            allele: "T",
            ref_allele: "C",
            flags: &["cds_start_NF".to_string()],
        };
        let result = evaluate(&input, &config, None);
        assert_eq!(result.confidence, "LC");
        assert!(result.filters.contains(&"INCOMPLETE_CDS".into()));
    }

    #[test]
    fn test_evaluate_small_intron() {
        // Build a transcript with a very small intron (10bp)
        let mut tr = make_multi_exon_transcript();
        tr.exons = vec![
            Exon {
                stable_id: "E1".into(),
                start: 1000,
                end: 1200,
                strand: Strand::Forward,
                phase: 0,
                end_phase: 0,
                rank: 1,
            },
            Exon {
                stable_id: "E2".into(),
                start: 1211, // intron is only 10bp (1201-1210)
                end: 1400,
                strand: Strand::Forward,
                phase: 0,
                end_phase: 0,
                rank: 2,
            },
        ];
        let config = LofteeConfig::default(); // min_intron_size = 15
        let input = LofteeInput {
            transcript: &tr,
            consequences: &[Consequence::SpliceDonorVariant],
            exon: None,
            intron: Some((0, 1)),
            cds_end: None,
            variant_start: 1201,
            variant_end: 1202,
            allele: "T",
            ref_allele: "G",
            flags: &[],
        };
        let result = evaluate(&input, &config, None);
        assert_eq!(result.confidence, "LC");
        assert!(result.filters.contains(&"SMALL_INTRON".into()));
    }

    #[test]
    fn test_evaluate_splice_donor_hc() {
        let tr = make_multi_exon_transcript();
        let config = LofteeConfig::default();
        let input = LofteeInput {
            transcript: &tr,
            consequences: &[Consequence::SpliceDonorVariant],
            exon: None,
            intron: Some((0, 2)),
            cds_end: None,
            variant_start: 1201,
            variant_end: 1202,
            allele: "T",
            ref_allele: "G",
            flags: &[],
        };
        let result = evaluate(&input, &config, None);
        assert_eq!(result.confidence, "HC");
        assert!(result.filters.is_empty());
    }

    #[test]
    fn test_evaluate_non_lof_returns_empty() {
        let tr = make_multi_exon_transcript();
        let config = LofteeConfig::default();
        let input = LofteeInput {
            transcript: &tr,
            consequences: &[Consequence::MissenseVariant],
            exon: Some((1, 3)),
            intron: None,
            cds_end: Some(300),
            variant_start: 2100,
            variant_end: 2100,
            allele: "T",
            ref_allele: "C",
            flags: &[],
        };
        let result = evaluate(&input, &config, None);
        assert!(result.confidence.is_empty());
    }
}
