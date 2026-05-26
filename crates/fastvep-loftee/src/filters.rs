//! LOFTEE filter implementations.
//!
//! Each filter returns a boolean indicating whether the filter condition is met
//! (true = variant fails the filter and should be downgraded).

use fastvep_cache::providers::SequenceProvider;
use fastvep_core::Strand;
use fastvep_genome::Transcript;

// --- Structural filters (no external data needed) ---

/// Check if the transcript has undefined exon/intron boundaries.
/// For exonic variants: exon number should be defined.
/// For intronic variants: intron number should be defined.
pub fn check_exon_intron_undef(exon: Option<(u32, u32)>) -> bool {
    // If we're told the variant is exonic but exon is None, it's undefined
    // The caller passes the exon/intron info from AlleleAnnotation
    exon.is_none()
}

/// Check if the intron annotation is missing/undefined.
pub fn check_intron_annotation_errors(intron: Option<(u32, u32)>) -> bool {
    intron.is_none()
}

/// Check if the transcript has only a single exon.
pub fn check_single_exon(total_exons: u32) -> bool {
    total_exons == 1
}

/// Check if the transcript has an incomplete CDS (missing start or stop codon).
pub fn check_incomplete_cds(flags: &[String]) -> bool {
    flags.iter().any(|f| f == "cds_start_NF" || f == "cds_end_NF")
}

/// Get CDS length (fast approximation using cDNA coordinates).
pub fn get_cds_length_fast(tr: &Transcript) -> u64 {
    match (tr.cdna_coding_start, tr.cdna_coding_end) {
        (Some(start), Some(end)) => end - start + 1,
        _ => 0,
    }
}

/// Get the intron size for a given intron index (0-based).
pub fn get_intron_size(tr: &Transcript, intron_idx: usize) -> u64 {
    let sorted = sorted_exons(tr);
    if intron_idx >= sorted.len().saturating_sub(1) {
        return 0;
    }
    match tr.strand {
        Strand::Forward => {
            let intron_start = sorted[intron_idx].end + 1;
            let intron_end = sorted[intron_idx + 1].start - 1;
            if intron_end >= intron_start {
                intron_end - intron_start + 1
            } else {
                0
            }
        }
        Strand::Reverse => {
            // sorted in descending start order for reverse strand
            let intron_start = sorted[intron_idx + 1].end + 1;
            let intron_end = sorted[intron_idx].start - 1;
            if intron_end >= intron_start {
                intron_end - intron_start + 1
            } else {
                0
            }
        }
    }
}

/// Get the coding length of the last exon (distance from last exon start to stop codon).
pub fn get_last_exon_coding_length(tr: &Transcript) -> u64 {
    let stop_pos = match tr.strand {
        Strand::Forward => tr.coding_region_end,
        Strand::Reverse => tr.coding_region_start,
    };
    let stop_pos = match stop_pos {
        Some(p) => p,
        None => return 0,
    };

    let sorted = sorted_exons(tr);
    // Find the exon containing the stop codon
    for exon in sorted.iter().rev() {
        match tr.strand {
            Strand::Forward => {
                if exon.start <= stop_pos && exon.end >= stop_pos {
                    return stop_pos - exon.start;
                }
            }
            Strand::Reverse => {
                if exon.start <= stop_pos && exon.end >= stop_pos {
                    return exon.end - stop_pos;
                }
            }
        }
    }
    0
}

/// Get the base-pair distance from the variant to the stop codon (in CDS space).
pub fn get_bp_dist_to_stop(tr: &Transcript, variant_genomic_pos: u64) -> u64 {
    // Convert variant genomic position to cDNA, then compute distance to coding end
    let variant_cdna = tr.genomic_to_cdna(variant_genomic_pos);
    let coding_end = tr.cdna_coding_end;

    match (variant_cdna, coding_end) {
        (Some(vpos), Some(cend)) => {
            if cend > vpos {
                cend - vpos
            } else {
                0
            }
        }
        _ => 0,
    }
}

// --- Sequence-based filters ---

/// Check if a splice donor variant is a GC→GT mutation.
/// This means the intron starts with GC (non-canonical) and the variant changes C→T,
/// converting it to the canonical GT donor.
pub fn check_gc_to_gt_donor(
    tr: &Transcript,
    intron_idx: usize,
    ref_allele: &str,
    alt_allele: &str,
    seq_provider: Option<&(dyn SequenceProvider + Send + Sync)>,
) -> bool {
    let sp = match seq_provider {
        Some(sp) => sp,
        None => return false,
    };

    let intron_start_seq = get_intron_start_dinuc(tr, intron_idx, sp);
    match intron_start_seq {
        Some(dinuc) => {
            // Intron starts with GC, and ref is C, alt is T → GC→GT
            dinuc == "GC" && ref_allele == "C" && alt_allele == "T"
        }
        None => false,
    }
}

/// Check if the intron has a non-canonical splice motif (not GT..AG).
pub fn check_non_canonical_splice(
    tr: &Transcript,
    intron_idx: usize,
    seq_provider: Option<&(dyn SequenceProvider + Send + Sync)>,
) -> bool {
    let sp = match seq_provider {
        Some(sp) => sp,
        None => return false,
    };

    let start_dinuc = get_intron_start_dinuc(tr, intron_idx, sp);
    let end_dinuc = get_intron_end_dinuc(tr, intron_idx, sp);

    let non_canonical_start = start_dinuc.as_deref().map_or(false, |d| d != "GT");
    let non_canonical_end = end_dinuc.as_deref().map_or(false, |d| d != "AG");

    non_canonical_start || non_canonical_end
}

/// Check for NAGNAG splice acceptor site (AG.AG motif around the splice site).
pub fn check_nagnag(
    tr: &Transcript,
    intron_idx: usize,
    seq_provider: Option<&(dyn SequenceProvider + Send + Sync)>,
) -> bool {
    let sp = match seq_provider {
        Some(sp) => sp,
        None => return false,
    };

    // Get 9bp context around splice acceptor (4bp before + 1bp variant + 4bp after)
    let sorted = sorted_exons(tr);
    if intron_idx >= sorted.len().saturating_sub(1) {
        return false;
    }

    // The acceptor site is at the start of the downstream exon
    let acceptor_pos = match tr.strand {
        Strand::Forward => sorted[intron_idx + 1].start,
        Strand::Reverse => sorted[intron_idx + 1].end,
    };

    // Fetch 9 bases centered on acceptor boundary (4 bases either side)
    let (fetch_start, fetch_end) = (acceptor_pos.saturating_sub(4), acceptor_pos + 4);
    let seq = sp
        .fetch_sequence(&tr.chromosome, fetch_start, fetch_end)
        .ok();

    match seq {
        Some(bytes) => {
            let seq_str = if tr.strand == Strand::Reverse {
                let rc: Vec<u8> = bytes.iter().rev().map(|&b| complement(b)).collect();
                String::from_utf8_lossy(&rc).to_uppercase()
            } else {
                String::from_utf8_lossy(&bytes).to_uppercase()
            };
            // Look for AG.AG pattern (. is any single nucleotide)
            if seq_str.len() >= 5 {
                for i in 0..seq_str.len().saturating_sub(4) {
                    let window = &seq_str[i..i + 5];
                    if window.starts_with("AG") && window.ends_with("AG") {
                        return true;
                    }
                }
            }
            false
        }
        None => false,
    }
}

/// Check if the variant is in the 5' UTR region (before coding start).
pub fn check_5utr_splice(tr: &Transcript, variant_start: u64, variant_end: u64) -> bool {
    match tr.strand {
        Strand::Forward => {
            if let Some(crs) = tr.coding_region_start {
                variant_end < crs
            } else {
                false
            }
        }
        Strand::Reverse => {
            if let Some(cre) = tr.coding_region_end {
                variant_start > cre
            } else {
                false
            }
        }
    }
}

/// Check if the variant is in the 3' UTR region (after coding end).
pub fn check_3utr_splice(tr: &Transcript, variant_start: u64, variant_end: u64) -> bool {
    match tr.strand {
        Strand::Forward => {
            if let Some(cre) = tr.coding_region_end {
                variant_start > cre
            } else {
                false
            }
        }
        Strand::Reverse => {
            if let Some(crs) = tr.coding_region_start {
                variant_end < crs
            } else {
                false
            }
        }
    }
}

/// Check if the ALT allele matches the ancestral allele (ANC_ALLELE filter).
pub fn check_ancestral_allele(
    chrom: &str,
    variant_start: u64,
    variant_end: u64,
    alt_allele: &str,
    ref_allele: &str,
    ancestral_provider: &(dyn SequenceProvider + Send + Sync),
) -> bool {
    // Only apply to SNPs for now (matching Perl: ignoring indels)
    if ref_allele.contains('-') || alt_allele.contains('-') || ref_allele.len() != 1 || alt_allele.len() != 1 {
        return false;
    }

    match ancestral_provider.fetch_sequence(chrom, variant_start, variant_end) {
        Ok(seq) => {
            let ancestral = String::from_utf8_lossy(&seq).to_uppercase();
            ancestral == alt_allele.to_uppercase()
        }
        Err(_) => false,
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

/// Get the first 2 bases of an intron (donor site dinucleotide).
fn get_intron_start_dinuc(
    tr: &Transcript,
    intron_idx: usize,
    sp: &(dyn SequenceProvider + Send + Sync),
) -> Option<String> {
    let sorted = sorted_exons(tr);
    if intron_idx >= sorted.len().saturating_sub(1) {
        return None;
    }
    let (start, end) = match tr.strand {
        Strand::Forward => (sorted[intron_idx].end + 1, sorted[intron_idx].end + 2),
        Strand::Reverse => (sorted[intron_idx].start - 2, sorted[intron_idx].start - 1),
    };
    let seq = sp.fetch_sequence(&tr.chromosome, start, end).ok()?;
    let s = if tr.strand == Strand::Reverse {
        let rc: Vec<u8> = seq.iter().rev().map(|&b| complement(b)).collect();
        String::from_utf8_lossy(&rc).to_uppercase()
    } else {
        String::from_utf8_lossy(&seq).to_uppercase()
    };
    Some(s)
}

/// Get the last 2 bases of an intron (acceptor site dinucleotide).
fn get_intron_end_dinuc(
    tr: &Transcript,
    intron_idx: usize,
    sp: &(dyn SequenceProvider + Send + Sync),
) -> Option<String> {
    let sorted = sorted_exons(tr);
    if intron_idx >= sorted.len().saturating_sub(1) {
        return None;
    }
    let (start, end) = match tr.strand {
        Strand::Forward => (sorted[intron_idx + 1].start - 2, sorted[intron_idx + 1].start - 1),
        Strand::Reverse => (sorted[intron_idx + 1].end + 1, sorted[intron_idx + 1].end + 2),
    };
    let seq = sp.fetch_sequence(&tr.chromosome, start, end).ok()?;
    let s = if tr.strand == Strand::Reverse {
        let rc: Vec<u8> = seq.iter().rev().map(|&b| complement(b)).collect();
        String::from_utf8_lossy(&rc).to_uppercase()
    } else {
        String::from_utf8_lossy(&seq).to_uppercase()
    };
    Some(s)
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

#[cfg(test)]
mod tests {
    use super::*;
    use fastvep_core::Strand;
    use fastvep_genome::Exon;

    fn make_transcript() -> Transcript {
        crate::tests::make_multi_exon_transcript()
    }

    #[test]
    fn test_check_single_exon() {
        assert!(check_single_exon(1));
        assert!(!check_single_exon(2));
        assert!(!check_single_exon(3));
    }

    #[test]
    fn test_check_incomplete_cds() {
        assert!(!check_incomplete_cds(&[]));
        assert!(check_incomplete_cds(&["cds_start_NF".to_string()]));
        assert!(check_incomplete_cds(&["cds_end_NF".to_string()]));
        assert!(check_incomplete_cds(&[
            "cds_start_NF".to_string(),
            "cds_end_NF".to_string(),
        ]));
        assert!(!check_incomplete_cds(&["some_other_flag".to_string()]));
    }

    #[test]
    fn test_check_exon_intron_undef() {
        assert!(check_exon_intron_undef(None));
        assert!(!check_exon_intron_undef(Some((1, 3))));
    }

    #[test]
    fn test_get_cds_length_fast() {
        let tr = make_transcript();
        // cdna_coding_end(952) - cdna_coding_start(51) + 1 = 902
        assert_eq!(get_cds_length_fast(&tr), 902);
    }

    #[test]
    fn test_get_intron_size() {
        let tr = make_transcript();
        // Intron 0: between exon1 (end=1200) and exon2 (start=2000)
        // Size = 2000 - 1200 - 1 = 799
        assert_eq!(get_intron_size(&tr, 0), 799);
        // Intron 1: between exon2 (end=2300) and exon3 (start=4000)
        // Size = 4000 - 2300 - 1 = 1699
        assert_eq!(get_intron_size(&tr, 1), 1699);
    }

    #[test]
    fn test_get_intron_size_small() {
        let mut tr = make_transcript();
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
                start: 1211, // intron is 10bp: 1201-1210
                end: 1400,
                strand: Strand::Forward,
                phase: 0,
                end_phase: 0,
                rank: 2,
            },
        ];
        assert_eq!(get_intron_size(&tr, 0), 10);
    }

    #[test]
    fn test_get_last_exon_coding_length() {
        let tr = make_transcript();
        // coding_region_end = 4500, last exon (E3) starts at 4000
        // last_exon_coding_length = 4500 - 4000 = 500
        assert_eq!(get_last_exon_coding_length(&tr), 500);
    }

    #[test]
    fn test_get_bp_dist_to_stop() {
        let tr = make_transcript();
        // Variant at genomic 2100 in exon 2
        // cdna pos of 2100 = 201 (exon1 len) + (2100 - 2000) + 1 = 302
        // cdna_coding_end = 952
        // dist = 952 - 302 = 650
        let dist = get_bp_dist_to_stop(&tr, 2100);
        assert_eq!(dist, 650);
    }

    #[test]
    fn test_5utr_splice_forward() {
        let tr = make_transcript();
        // coding_region_start = 1050
        // Variant at 1020-1020 is before CDS start → 5'UTR
        assert!(check_5utr_splice(&tr, 1020, 1020));
        // Variant at 1100-1100 is after CDS start → not 5'UTR
        assert!(!check_5utr_splice(&tr, 1100, 1100));
    }

    #[test]
    fn test_3utr_splice_forward() {
        let tr = make_transcript();
        // coding_region_end = 4500
        // Variant at 4600-4600 is after CDS end → 3'UTR
        assert!(check_3utr_splice(&tr, 4600, 4600));
        // Variant at 4400-4400 is before CDS end → not 3'UTR
        assert!(!check_3utr_splice(&tr, 4400, 4400));
    }

    #[test]
    fn test_5utr_3utr_reverse_strand() {
        let mut tr = make_transcript();
        tr.strand = Strand::Reverse;
        // For reverse strand:
        // 5'UTR is genomically AFTER coding_region_end
        // 3'UTR is genomically BEFORE coding_region_start
        assert!(check_5utr_splice(&tr, 4600, 4600)); // after CRE (4500)
        assert!(!check_5utr_splice(&tr, 1020, 1020));
        assert!(check_3utr_splice(&tr, 1020, 1020)); // before CRS (1050)
        assert!(!check_3utr_splice(&tr, 4600, 4600));
    }

    #[test]
    fn test_ancestral_allele_snp() {
        // Mock ancestral sequence provider
        struct MockAncProvider;
        impl SequenceProvider for MockAncProvider {
            fn fetch_sequence(
                &self,
                _chrom: &str,
                _start: u64,
                _end: u64,
            ) -> anyhow::Result<Vec<u8>> {
                Ok(b"T".to_vec())
            }
        }
        let provider = MockAncProvider;

        // ALT matches ancestral → ANC_ALLELE should trigger
        assert!(check_ancestral_allele("chr1", 100, 100, "T", "C", &provider));
        // ALT does not match ancestral
        assert!(!check_ancestral_allele("chr1", 100, 100, "A", "C", &provider));
        // Indel → skip
        assert!(!check_ancestral_allele("chr1", 100, 100, "-", "C", &provider));
    }
}
