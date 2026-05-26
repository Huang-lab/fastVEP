//! Extended splice site prediction and de novo donor detection.
//!
//! Ports the `extended_splice.pl` and `de_novo_donor.pl` LOFTEE algorithms.

use crate::maxentscan::MaxEntScanData;
use crate::motifs::{MotifSet, MotifType};
use crate::svm::{self, SvmModel};
use fastvep_cache::providers::SequenceProvider;
use fastvep_core::Strand;
use fastvep_genome::Transcript;
use std::collections::HashMap;
use std::path::Path;

/// Pre-loaded splice prediction data.
pub struct SpliceData {
    pub mes: MaxEntScanData,
    pub donor_motifs: MotifSet,
    pub acceptor_motifs: MotifSet,
    pub donor_model: HashMap<String, f64>,
    pub acceptor_model: HashMap<String, f64>,
    pub donor_svm: SvmModel,
    // Thresholds
    pub donor_disruption_cutoff: f64,
    pub acceptor_disruption_cutoff: f64,
    pub donor_disruption_mes_cutoff: f64,
    pub acceptor_disruption_mes_cutoff: f64,
    pub donor_rescue_cutoff: f64,
    pub acceptor_rescue_cutoff: f64,
    pub weak_donor_cutoff: f64,
    pub max_scan_distance: usize,
    pub max_denovo_donor_distance: usize,
    pub denovo_donor_cutoff: f64,
    pub sre_flanksize: usize,
    pub exonic_denovo_only: bool,
}

impl SpliceData {
    /// Load all splice prediction data from a LOFTEE data directory.
    pub fn load(loftee_dir: &Path) -> Result<Self, String> {
        let mes = MaxEntScanData::load(loftee_dir)?;
        let splice_data_dir = loftee_dir.join("splice_data");
        let donor_motifs = MotifSet::load(&splice_data_dir.join("donor_motifs"))?;
        let acceptor_motifs = MotifSet::load(&splice_data_dir.join("acceptor_motifs"))?;
        let donor_model = svm::load_logreg_model(&splice_data_dir.join("donor_model.txt"))?;
        let acceptor_model = svm::load_logreg_model(&splice_data_dir.join("acceptor_model.txt"))?;
        let donor_svm = SvmModel::load(&splice_data_dir.join("de_novo_donor_SVM"))?;

        Ok(Self {
            mes,
            donor_motifs,
            acceptor_motifs,
            donor_model,
            acceptor_model,
            donor_svm,
            donor_disruption_cutoff: 0.98,
            acceptor_disruption_cutoff: 0.99,
            donor_disruption_mes_cutoff: 6.0,
            acceptor_disruption_mes_cutoff: 7.0,
            donor_rescue_cutoff: 8.5,
            acceptor_rescue_cutoff: 8.5,
            weak_donor_cutoff: -4.0,
            max_scan_distance: 15,
            max_denovo_donor_distance: 200,
            denovo_donor_cutoff: 0.995,
            sre_flanksize: 100,
            exonic_denovo_only: true,
        })
    }
}

/// Result of extended splice analysis.
#[derive(Debug, Clone)]
pub struct SpliceResult {
    /// Whether the variant disrupts the splice site.
    pub disrupting: bool,
    /// Type of splice site affected ("DONOR" or "ACCEPTOR").
    pub splice_type: String,
    /// Features computed during analysis.
    pub features: HashMap<String, String>,
    /// Info strings for LoF_info field.
    pub info: Vec<String>,
}

/// Result of de novo donor analysis.
#[derive(Debug, Clone)]
pub struct DeNovoDonorResult {
    /// SVM probability that a de novo donor is created.
    pub probability: f64,
    /// Features computed during analysis.
    pub features: HashMap<String, String>,
    /// Whether this is a loss of function (frameshift or stop introduced).
    pub lof: bool,
    /// Position of the de novo junction relative to the original.
    pub lof_pos: u64,
}

/// Check whether a variant disrupts an existing splice site.
/// Returns None if the variant doesn't overlap the extended splice region.
pub fn get_effect_on_splice(
    tr: &Transcript,
    intron_idx: usize,
    variant_start: u64,
    variant_end: u64,
    ref_allele: &str,
    alt_allele: &str,
    is_vep_splice_lof: bool,
    splice_data: &SpliceData,
    seq_provider: Option<&(dyn SequenceProvider + Send + Sync)>,
) -> Option<SpliceResult> {
    let sp = seq_provider?;
    let sorted = sorted_exons(tr);
    if intron_idx >= sorted.len().saturating_sub(1) {
        return None;
    }

    // Determine intron boundaries
    let (intron_start, intron_end) = match tr.strand {
        Strand::Forward => (
            sorted[intron_idx].end + 1,
            sorted[intron_idx + 1].start - 1,
        ),
        Strand::Reverse => (
            sorted[intron_idx + 1].end + 1,
            sorted[intron_idx].start - 1,
        ),
    };

    // Check overlap with donor region (9bp around intron start: -3 to +5)
    let donor_start = intron_start.saturating_sub(3);
    let donor_end = intron_start + 5;
    let overlaps_donor = variant_start <= donor_end && variant_end >= donor_start;

    // Check overlap with acceptor region (23bp around intron end: -19 to +3)
    let acceptor_start = intron_end.saturating_sub(19);
    let acceptor_end = intron_end + 3;
    let overlaps_acceptor = variant_start <= acceptor_end && variant_end >= acceptor_start;

    if !overlaps_donor && !overlaps_acceptor {
        return None;
    }

    let mut features = HashMap::new();
    let mut info = Vec::new();

    let splice_type = if overlaps_donor { "DONOR" } else { "ACCEPTOR" };

    // Get reference and variant splice site sequences
    let (ref_seq_start, ref_seq_end, _seq_len) = if overlaps_donor {
        (intron_start.saturating_sub(3), intron_start + 5, 9usize)
    } else {
        (intron_end.saturating_sub(19), intron_end + 3, 23usize)
    };

    let ref_bytes = sp
        .fetch_sequence(&tr.chromosome, ref_seq_start, ref_seq_end)
        .ok()?;

    // Handle strand
    let ref_seq = if tr.strand == Strand::Reverse {
        reverse_complement(&ref_bytes)
    } else {
        ref_bytes.iter().map(|b| b.to_ascii_uppercase()).collect()
    };

    // Score reference
    let ref_mes = if overlaps_donor {
        splice_data.mes.score_donor(&ref_seq)?
    } else {
        splice_data.mes.score_acceptor(&ref_seq)?
    };

    // Create variant sequence by mutating the reference
    let var_seq = mutate_splice_seq(
        &ref_seq,
        variant_start,
        variant_end,
        ref_allele,
        alt_allele,
        ref_seq_start,
        tr.strand,
    );

    let var_mes = if overlaps_donor {
        splice_data.mes.score_donor(&var_seq).unwrap_or(f64::NEG_INFINITY)
    } else {
        splice_data.mes.score_acceptor(&var_seq).unwrap_or(f64::NEG_INFINITY)
    };

    let mes_diff = ref_mes - var_mes;

    features.insert(
        format!("{}_MES_DIFF", splice_type),
        format!("{:.4}", mes_diff),
    );
    features.insert(
        format!("MUTANT_{}_MES", splice_type),
        format!("{:.4}", var_mes),
    );

    // Compute SRE features
    let motifs = if overlaps_donor {
        &splice_data.donor_motifs
    } else {
        &splice_data.acceptor_motifs
    };

    let flank_start = ref_seq_start.saturating_sub(splice_data.sre_flanksize as u64);
    let flank_end = ref_seq_end + splice_data.sre_flanksize as u64;
    if let Ok(flank_bytes) = sp.fetch_sequence(&tr.chromosome, flank_start, flank_end) {
        let flank = if tr.strand == Strand::Reverse {
            reverse_complement(&flank_bytes)
        } else {
            flank_bytes.iter().map(|b| b.to_ascii_uppercase()).collect()
        };

        let ese_count = motifs.scan(&flank, MotifType::Ese);
        let ess_count = motifs.scan(&flank, MotifType::Ess);
        let ise_count = motifs.scan(&flank, MotifType::Ise);
        let iss_count = motifs.scan(&flank, MotifType::Iss);

        features.insert(format!("{}_ESE", splice_type), ese_count.to_string());
        features.insert(format!("{}_ESS", splice_type), ess_count.to_string());
        features.insert(format!("{}_ISE", splice_type), ise_count.to_string());
        features.insert(format!("{}_ISS", splice_type), iss_count.to_string());
    }

    // Determine disruption using MES-only mode (no GERP in Phase 3)
    let cutoff = if overlaps_donor {
        splice_data.donor_disruption_mes_cutoff
    } else {
        splice_data.acceptor_disruption_mes_cutoff
    };
    let disrupting = mes_diff > cutoff;

    features.insert(
        format!("{}_DISRUPTION_PROB", splice_type),
        format!("{:.4}", mes_diff),
    );

    if disrupting {
        info.push(format!("{}_DISRUPTING", splice_type));
    } else if is_vep_splice_lof {
        // VEP called it splice_donor/acceptor_variant but MES says not disrupted
        info.push(format!("NON_{}_DISRUPTING", splice_type));
    }

    Some(SpliceResult {
        disrupting,
        splice_type: splice_type.to_string(),
        features,
        info,
    })
}

/// Check for de novo donor splice site creation.
pub fn check_for_denovo_donor(
    tr: &Transcript,
    exon_idx: usize,
    variant_start: u64,
    variant_end: u64,
    ref_allele: &str,
    alt_allele: &str,
    splice_data: &SpliceData,
    seq_provider: Option<&(dyn SequenceProvider + Send + Sync)>,
) -> Option<DeNovoDonorResult> {
    let sp = seq_provider?;
    let sorted = sorted_exons(tr);
    let n_exons = sorted.len();

    // Skip if in last exon (no downstream intron)
    if exon_idx >= n_exons.saturating_sub(1) {
        return None;
    }

    let exon = sorted[exon_idx];
    let exon_length = (exon.end - exon.start + 1) as usize;

    // Check distance from exon boundary
    let dist_from_boundary = match tr.strand {
        Strand::Forward => {
            if variant_start > exon.end {
                (variant_start - exon.end) as usize
            } else {
                0
            }
        }
        Strand::Reverse => {
            if variant_end < exon.start {
                (exon.start - variant_end) as usize
            } else {
                0
            }
        }
    };

    if splice_data.exonic_denovo_only && dist_from_boundary > 6 {
        return None;
    }
    if dist_from_boundary > splice_data.max_denovo_donor_distance {
        return None;
    }

    // Get sequence context: 200bp flanking
    let flank = 200u64;
    let (intron_start, intron_end) = match tr.strand {
        Strand::Forward => (exon.end + 1, sorted[exon_idx + 1].start - 1),
        Strand::Reverse => (sorted[exon_idx + 1].end + 1, exon.start - 1),
    };

    let seq_start = match tr.strand {
        Strand::Forward => exon.start.saturating_sub(flank),
        Strand::Reverse => intron_start.saturating_sub(flank),
    };
    let seq_end = match tr.strand {
        Strand::Forward => intron_end + flank,
        Strand::Reverse => exon.end + flank,
    };

    let raw_seq = sp.fetch_sequence(&tr.chromosome, seq_start, seq_end).ok()?;
    let full_seq = if tr.strand == Strand::Reverse {
        reverse_complement(&raw_seq)
    } else {
        raw_seq.iter().map(|b| b.to_ascii_uppercase()).collect()
    };

    // Reference donor position
    let ref_junc = flank as usize + exon_length;
    if ref_junc + 6 > full_seq.len() || ref_junc < 3 {
        return None;
    }

    let ref_donor = &full_seq[ref_junc - 3..ref_junc + 6];
    let ref_mes = splice_data.mes.score_donor(ref_donor)?;

    if ref_mes < splice_data.weak_donor_cutoff {
        return None;
    }
    // Validate GT consensus
    if ref_donor[3] != b'G' || ref_donor[4] != b'T' {
        return None;
    }

    // Compute variant position relative to sequence start
    let var_pos_in_seq = match tr.strand {
        Strand::Forward => (variant_start - seq_start) as usize,
        Strand::Reverse => (seq_end - variant_end) as usize,
    };

    // Mutate sequence
    let (var_seq, nt_delta) = mutate_full_seq(
        &full_seq,
        var_pos_in_seq,
        ref_allele,
        alt_allele,
        tr.strand,
    );

    let effective_exon_length = (exon_length as i64 + nt_delta) as usize;
    let new_ref_junc = flank as usize + effective_exon_length;

    // Scan for de novo donor sites
    let adj = if nt_delta > 0 { nt_delta as usize } else { 0 };
    let scan_start = var_pos_in_seq.saturating_sub(8);
    let scan_end = (var_pos_in_seq + adj + 1).min(var_seq.len().saturating_sub(9));

    let mut best_prob = 0.0;
    let mut best_features = HashMap::new();
    let mut best_lof = false;
    let mut best_lof_pos = 0u64;

    for pos in scan_start..=scan_end {
        if pos + 9 > var_seq.len() {
            continue;
        }
        let new_junc = pos + 3;
        if new_junc == new_ref_junc {
            continue;
        }

        let current_exon_len = new_junc.saturating_sub(flank as usize);
        if current_exon_len == 0 {
            continue;
        }

        // Check minimum intron length
        let intron_len = var_seq.len().saturating_sub(flank as usize + new_junc);
        if intron_len < 70 {
            continue;
        }

        let candidate = &var_seq[pos..pos + 9];
        let alt_mes = match splice_data.mes.score_donor(candidate) {
            Some(s) => s,
            None => continue,
        };

        // Filter: must not be dramatically worse than reference
        if alt_mes - ref_mes < -15.0 {
            continue;
        }

        // Compute SRE features
        let flanksize = splice_data.sre_flanksize;
        let ese_delta = compute_sre_delta(
            &var_seq,
            new_junc,
            new_ref_junc,
            flanksize,
            &splice_data.donor_motifs,
            MotifType::Ese,
            true,
        );
        let ess_delta = compute_sre_delta(
            &var_seq,
            new_junc,
            new_ref_junc,
            flanksize,
            &splice_data.donor_motifs,
            MotifType::Ess,
            true,
        );
        let ise_delta = compute_sre_delta(
            &var_seq,
            new_junc,
            new_ref_junc,
            flanksize,
            &splice_data.donor_motifs,
            MotifType::Ise,
            false,
        );
        let iss_delta = compute_sre_delta(
            &var_seq,
            new_junc,
            new_ref_junc,
            flanksize,
            &splice_data.donor_motifs,
            MotifType::Iss,
            false,
        );

        let mes_diff = if new_junc > new_ref_junc {
            alt_mes - ref_mes
        } else {
            ref_mes - alt_mes
        };

        let mut svm_features = HashMap::new();
        svm_features.insert("ese".to_string(), ese_delta);
        svm_features.insert("ess".to_string(), ess_delta);
        svm_features.insert("ise".to_string(), ise_delta);
        svm_features.insert("iss".to_string(), iss_delta);
        svm_features.insert("MESdiff".to_string(), mes_diff);

        let mut pr = splice_data.donor_svm.predict(&svm_features);
        // Flip probability for exon truncation
        if new_junc < new_ref_junc {
            pr = 1.0 - pr;
        }

        if pr > best_prob {
            best_prob = pr;
            let delta = new_junc as i64 - new_ref_junc as i64;
            let lof = if delta.unsigned_abs() as usize % 3 != 0 {
                true // frameshift
            } else {
                check_stop_introduced(&var_seq, new_junc, new_ref_junc)
            };
            best_lof = lof;
            best_lof_pos = match tr.strand {
                Strand::Forward => exon.start + new_junc as u64 - flank,
                Strand::Reverse => exon.end - new_junc as u64 + flank,
            };
            best_features.clear();
            best_features.insert("ese".to_string(), format!("{:.4}", ese_delta));
            best_features.insert("ess".to_string(), format!("{:.4}", ess_delta));
            best_features.insert("ise".to_string(), format!("{:.4}", ise_delta));
            best_features.insert("iss".to_string(), format!("{:.4}", iss_delta));
            best_features.insert("MESdiff".to_string(), format!("{:.4}", mes_diff));
            best_features.insert("de_novo_donor_prob".to_string(), format!("{:.6}", pr));
        }
    }

    if best_prob > 0.0 {
        Some(DeNovoDonorResult {
            probability: best_prob,
            features: best_features,
            lof: best_lof,
            lof_pos: best_lof_pos,
        })
    } else {
        None
    }
}

// --- Helper functions ---

fn sorted_exons(tr: &Transcript) -> Vec<&fastvep_genome::Exon> {
    let mut exons: Vec<&fastvep_genome::Exon> = tr.exons.iter().collect();
    match tr.strand {
        Strand::Forward => exons.sort_by_key(|e| e.start),
        Strand::Reverse => exons.sort_by(|a, b| b.start.cmp(&a.start)),
    }
    exons
}

fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' | b'a' => b'T',
            b'T' | b't' => b'A',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            other => other,
        })
        .collect()
}

/// Mutate a splice site sequence (9bp or 23bp) at the variant position.
fn mutate_splice_seq(
    ref_seq: &[u8],
    variant_start: u64,
    _variant_end: u64,
    _ref_allele: &str,
    alt_allele: &str,
    seq_genomic_start: u64,
    strand: Strand,
) -> Vec<u8> {
    let mut result = ref_seq.to_vec();
    let offset = match strand {
        Strand::Forward => (variant_start - seq_genomic_start) as usize,
        Strand::Reverse => (ref_seq.len() - 1) - (variant_start - seq_genomic_start) as usize,
    };
    if offset < result.len() && alt_allele.len() == 1 {
        let alt_byte = if strand == Strand::Reverse {
            complement(alt_allele.as_bytes()[0])
        } else {
            alt_allele.as_bytes()[0].to_ascii_uppercase()
        };
        result[offset] = alt_byte;
    }
    result
}

/// Mutate a full-length sequence for de novo donor analysis.
fn mutate_full_seq(
    seq: &[u8],
    var_pos: usize,
    ref_allele: &str,
    alt_allele: &str,
    strand: Strand,
) -> (Vec<u8>, i64) {
    let mut result = seq.to_vec();
    let alt_bases: Vec<u8> = if strand == Strand::Reverse {
        alt_allele
            .as_bytes()
            .iter()
            .rev()
            .map(|&b| complement(b))
            .collect()
    } else {
        alt_allele.as_bytes().iter().map(|b| b.to_ascii_uppercase()).collect()
    };

    if ref_allele == "-" {
        // Insertion
        for (j, &b) in alt_bases.iter().enumerate() {
            result.insert(var_pos + j, b);
        }
        (result, alt_bases.len() as i64)
    } else if alt_allele == "-" {
        // Deletion
        let del_len = ref_allele.len().min(result.len() - var_pos);
        result.drain(var_pos..var_pos + del_len);
        (result, -(del_len as i64))
    } else {
        // SNP or MNP
        for (j, &b) in alt_bases.iter().enumerate() {
            if var_pos + j < result.len() {
                result[var_pos + j] = b;
            }
        }
        (result, alt_bases.len() as i64 - ref_allele.len() as i64)
    }
}

fn complement(b: u8) -> u8 {
    match b {
        b'A' | b'a' => b'T',
        b'T' | b't' => b'A',
        b'C' | b'c' => b'G',
        b'G' | b'g' => b'C',
        other => other,
    }
}

/// Compute SRE motif count delta between new and reference junctions.
fn compute_sre_delta(
    seq: &[u8],
    new_junc: usize,
    ref_junc: usize,
    flanksize: usize,
    motifs: &MotifSet,
    motif_type: MotifType,
    upstream: bool,
) -> f64 {
    let len = seq.len();
    let (new_start, new_end) = if upstream {
        (new_junc.saturating_sub(flanksize), new_junc.min(len))
    } else {
        (new_junc.min(len), (new_junc + flanksize).min(len))
    };
    let (ref_start, ref_end) = if upstream {
        (ref_junc.saturating_sub(flanksize), ref_junc.min(len))
    } else {
        (ref_junc.min(len), (ref_junc + flanksize).min(len))
    };

    let new_count = if new_start < new_end {
        motifs.scan(&seq[new_start..new_end], motif_type)
    } else {
        0
    };
    let ref_count = if ref_start < ref_end {
        motifs.scan(&seq[ref_start..ref_end], motif_type)
    } else {
        0
    };

    new_count as f64 - ref_count as f64
}

/// Check if exon extension/truncation introduces a stop codon.
fn check_stop_introduced(seq: &[u8], new_junc: usize, ref_junc: usize) -> bool {
    let (start, end) = if new_junc > ref_junc {
        (ref_junc, new_junc)
    } else {
        (new_junc, ref_junc)
    };

    // Check codons in the extended/truncated region
    let region = &seq[start..end.min(seq.len())];
    for i in (0..region.len()).step_by(3) {
        if i + 3 <= region.len() {
            let codon = &region[i..i + 3];
            let codon_upper: Vec<u8> = codon.iter().map(|b| b.to_ascii_uppercase()).collect();
            if codon_upper == b"TAG" || codon_upper == b"TAA" || codon_upper == b"TGA" {
                return true;
            }
        }
    }
    false
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reverse_complement() {
        assert_eq!(reverse_complement(b"ACGT"), b"ACGT");
        assert_eq!(reverse_complement(b"AAAA"), b"TTTT");
        assert_eq!(reverse_complement(b"GCAT"), b"ATGC");
    }

    #[test]
    fn test_mutate_full_seq_snp() {
        let seq = b"ACGTACGT";
        let (result, delta) = mutate_full_seq(seq, 2, "G", "A", Strand::Forward);
        assert_eq!(result, b"ACATACGT");
        assert_eq!(delta, 0);
    }

    #[test]
    fn test_mutate_full_seq_insertion() {
        let seq = b"ACGTACGT";
        let (result, delta) = mutate_full_seq(seq, 2, "-", "AA", Strand::Forward);
        assert_eq!(result, b"ACAAGTACGT");
        assert_eq!(delta, 2);
    }

    #[test]
    fn test_mutate_full_seq_deletion() {
        let seq = b"ACGTACGT";
        let (result, delta) = mutate_full_seq(seq, 2, "GT", "-", Strand::Forward);
        assert_eq!(result, b"ACACGT");
        assert_eq!(delta, -2);
    }

    #[test]
    fn test_check_stop_introduced() {
        // Region contains TAG
        assert!(check_stop_introduced(b"AAATAGGGG", 0, 6));
        // No stop codon
        assert!(!check_stop_introduced(b"AAAGCTGGG", 0, 6));
        // TAA
        assert!(check_stop_introduced(b"TAAGGGGGG", 0, 3));
        // TGA
        assert!(check_stop_introduced(b"TGAGGGGGG", 0, 3));
    }
}
