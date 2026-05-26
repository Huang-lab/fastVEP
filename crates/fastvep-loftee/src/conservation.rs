//! Conservation data readers for GERP and PhyloCSF.
//!
//! Phase 4: Implements GERP-weighted distance for END_TRUNC and
//! PhyloCSF-based PHYLOCSF_WEAK / PHYLOCSF_UNLIKELY_ORF flags.

use fastvep_core::Strand;
use fastvep_genome::Transcript;
use std::collections::HashMap;

/// Trait for looking up per-base conservation scores.
pub trait ConservationProvider: Send + Sync {
    /// Fetch the GERP score for a single position.
    fn gerp_score(&self, chrom: &str, pos: u64) -> Option<f64>;

    /// Sum GERP scores across a range [start, end] inclusive.
    fn gerp_sum(&self, chrom: &str, start: u64, end: u64) -> f64 {
        let mut sum = 0.0;
        for pos in start..=end {
            if let Some(score) = self.gerp_score(chrom, pos) {
                sum += score;
            }
        }
        sum
    }
}

/// Trait for looking up PhyloCSF exon-level scores.
pub trait PhyloCsfProvider: Send + Sync {
    /// Look up PhyloCSF scores for a given transcript + exon.
    /// Returns (corresponding_orf_score, max_score) or None if too short.
    fn phylocsf_scores(
        &self,
        transcript_id: &str,
        exon_number: u32,
    ) -> Option<PhyloCsfResult>;
}

#[derive(Debug, Clone)]
pub struct PhyloCsfResult {
    pub corresponding_orf_score: f64,
    pub max_score: f64,
}

/// Compute GERP-weighted distance from a variant to the stop codon.
///
/// Iterates through exons from the variant position to the stop codon,
/// weighting each base by its GERP score. This is the Phase 4 upgrade
/// to the simple base-pair distance used in Phase 1.
///
/// Returns (gerp_weighted_dist, bp_dist).
pub fn get_gerp_weighted_dist(
    tr: &Transcript,
    variant_pos: u64,
    gerp_provider: &dyn ConservationProvider,
) -> (f64, u64) {
    let sorted = sorted_exons(tr);

    let stop_pos = match tr.strand {
        Strand::Forward => tr.coding_region_end,
        Strand::Reverse => tr.coding_region_start,
    };
    let stop_pos = match stop_pos {
        Some(p) => p,
        None => return (0.0, 0),
    };

    let mut gerp_dist = 0.0;
    let mut bp_dist = 0u64;

    for exon in &sorted {
        // Determine the coding portion of this exon
        let (exon_coding_start, exon_coding_end) = match tr.strand {
            Strand::Forward => {
                let cs = exon.start.max(tr.coding_region_start.unwrap_or(exon.start));
                let ce = exon.end.min(stop_pos);
                if cs > ce || exon.end < tr.coding_region_start.unwrap_or(0) {
                    continue;
                }
                (cs, ce)
            }
            Strand::Reverse => {
                let cs = exon.start.max(stop_pos);
                let ce = exon.end.min(tr.coding_region_end.unwrap_or(exon.end));
                if cs > ce || exon.start > tr.coding_region_end.unwrap_or(u64::MAX) {
                    continue;
                }
                (cs, ce)
            }
        };

        // Determine the region in this exon that's after the variant
        let region_start = match tr.strand {
            Strand::Forward => variant_pos.max(exon_coding_start),
            Strand::Reverse => exon_coding_start,
        };
        let region_end = match tr.strand {
            Strand::Forward => exon_coding_end,
            Strand::Reverse => variant_pos.min(exon_coding_end),
        };

        if region_start > region_end {
            continue;
        }

        let region_bp = region_end - region_start + 1;
        bp_dist += region_bp;
        gerp_dist += gerp_provider.gerp_sum(&tr.chromosome, region_start, region_end);
    }

    (gerp_dist, bp_dist)
}

/// Check PhyloCSF conservation for an exonic variant.
/// Returns (flag_name, info_entries) or None if no flag triggered.
pub fn check_phylocsf(
    tr: &Transcript,
    exon_number: u32,
    phylocsf_provider: &dyn PhyloCsfProvider,
) -> Option<(String, HashMap<String, String>)> {
    let result = phylocsf_provider.phylocsf_scores(&tr.stable_id, exon_number)?;

    let mut info = HashMap::new();
    info.insert(
        "ANN_ORF".to_string(),
        format!("{:.4}", result.corresponding_orf_score),
    );
    info.insert("MAX_ORF".to_string(), format!("{:.4}", result.max_score));

    if result.corresponding_orf_score < 0.0 {
        let flag = if result.max_score > 0.0 {
            "PHYLOCSF_UNLIKELY_ORF"
        } else {
            "PHYLOCSF_WEAK"
        };
        Some((flag.to_string(), info))
    } else {
        None
    }
}

fn sorted_exons(tr: &Transcript) -> Vec<&fastvep_genome::Exon> {
    let mut exons: Vec<&fastvep_genome::Exon> = tr.exons.iter().collect();
    match tr.strand {
        Strand::Forward => exons.sort_by_key(|e| e.start),
        Strand::Reverse => exons.sort_by(|a, b| b.start.cmp(&a.start)),
    }
    exons
}

#[cfg(test)]
mod tests {
    use super::*;

    struct MockGerpProvider {
        scores: HashMap<u64, f64>,
    }

    impl ConservationProvider for MockGerpProvider {
        fn gerp_score(&self, _chrom: &str, pos: u64) -> Option<f64> {
            self.scores.get(&pos).copied()
        }
    }

    struct MockPhyloCsfProvider;

    impl PhyloCsfProvider for MockPhyloCsfProvider {
        fn phylocsf_scores(&self, _tid: &str, exon: u32) -> Option<PhyloCsfResult> {
            match exon {
                1 => Some(PhyloCsfResult {
                    corresponding_orf_score: -5.0,
                    max_score: 10.0,
                }),
                2 => Some(PhyloCsfResult {
                    corresponding_orf_score: -5.0,
                    max_score: -3.0,
                }),
                3 => Some(PhyloCsfResult {
                    corresponding_orf_score: 5.0,
                    max_score: 10.0,
                }),
                _ => None,
            }
        }
    }

    #[test]
    fn test_check_phylocsf_unlikely_orf() {
        let tr = crate::tests::make_multi_exon_transcript();
        let provider = MockPhyloCsfProvider;
        // Exon 1: corr_orf < 0, max > 0 → PHYLOCSF_UNLIKELY_ORF
        let result = check_phylocsf(&tr, 1, &provider);
        assert!(result.is_some());
        let (flag, _) = result.unwrap();
        assert_eq!(flag, "PHYLOCSF_UNLIKELY_ORF");
    }

    #[test]
    fn test_check_phylocsf_weak() {
        let tr = crate::tests::make_multi_exon_transcript();
        let provider = MockPhyloCsfProvider;
        // Exon 2: corr_orf < 0, max < 0 → PHYLOCSF_WEAK
        let result = check_phylocsf(&tr, 2, &provider);
        assert!(result.is_some());
        let (flag, _) = result.unwrap();
        assert_eq!(flag, "PHYLOCSF_WEAK");
    }

    #[test]
    fn test_check_phylocsf_no_flag() {
        let tr = crate::tests::make_multi_exon_transcript();
        let provider = MockPhyloCsfProvider;
        // Exon 3: corr_orf > 0 → no flag
        let result = check_phylocsf(&tr, 3, &provider);
        assert!(result.is_none());
    }

    #[test]
    fn test_gerp_weighted_dist() {
        let tr = crate::tests::make_multi_exon_transcript();
        let mut scores = HashMap::new();
        // Give some positions high GERP scores
        for pos in 4400..=4500 {
            scores.insert(pos, 2.0);
        }
        let provider = MockGerpProvider { scores };
        let (gerp_dist, bp_dist) = get_gerp_weighted_dist(&tr, 4400, &provider);
        // 4400 to 4500 = 101 positions, each with GERP 2.0
        assert_eq!(bp_dist, 101);
        assert!((gerp_dist - 202.0).abs() < 0.001);
    }
}
