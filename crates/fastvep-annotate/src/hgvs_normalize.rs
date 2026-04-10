//! HGVS normalization utilities: insertion-to-duplication conversion and 3' shifting.
//!
//! These functions post-process HGVSc strings to convert insertion notation (`ins`)
//! to duplication notation (`dup`) when the inserted bases match the adjacent reference,
//! and to 3'-shift intronic variants per HGVS conventions.

use fastvep_cache::providers::SequenceProvider;

/// Convert an intronic insertion HGVSc to duplication notation (coding transcript).
pub fn convert_ins_to_dup(
    hgvsc: &str,
    intron_offset: i64,
    ins_len: u64,
    nearest_exon_cdna_pos: u64,
    coding_start: u64,
    coding_end: Option<u64>,
) -> Option<String> {
    let prefix_end = hgvsc
        .find(":c.")
        .map(|i| i + 3)
        .or_else(|| hgvsc.find(":n.").map(|i| i + 3))?;
    let prefix = &hgvsc[..prefix_end];

    let build_pos = |cdna: u64, off: i64| -> String {
        let raw = cdna as i64 - coding_start as i64 + 1;
        let cp = if raw <= 0 { raw - 1 } else { raw };
        if cp < 0 {
            if off > 0 {
                format!("{}+{}", cp, off)
            } else {
                format!("{}{}", cp, off)
            }
        } else if coding_end.is_some() && cdna > coding_end.unwrap() {
            let u = cdna - coding_end.unwrap();
            if off > 0 {
                format!("*{}+{}", u, off)
            } else {
                format!("*{}{}", u, off)
            }
        } else if off > 0 {
            format!("{}+{}", cp, off)
        } else {
            format!("{}{}", cp, off)
        }
    };

    if ins_len == 1 {
        let pos = build_pos(nearest_exon_cdna_pos, intron_offset);
        Some(format!("{}{}dup", prefix, pos))
    } else {
        let start_offset = intron_offset - ins_len as i64 + 1;
        let start_pos = build_pos(nearest_exon_cdna_pos, start_offset);
        let end_pos = build_pos(nearest_exon_cdna_pos, intron_offset);
        Some(format!("{}{}_{}dup", prefix, start_pos, end_pos))
    }
}

/// Convert intronic insertion to dup notation with explicit start/end offsets (coding).
pub fn convert_ins_to_dup_range(
    hgvsc: &str,
    start_offset: i64,
    end_offset: i64,
    nearest_exon_cdna_pos: u64,
    coding_start: u64,
    coding_end: Option<u64>,
) -> Option<String> {
    let prefix_end = hgvsc
        .find(":c.")
        .map(|i| i + 3)
        .or_else(|| hgvsc.find(":n.").map(|i| i + 3))?;
    let prefix = &hgvsc[..prefix_end];

    let build_pos = |cdna: u64, off: i64| -> String {
        let raw = cdna as i64 - coding_start as i64 + 1;
        let cp = if raw <= 0 { raw - 1 } else { raw };
        if cp < 0 {
            if off > 0 {
                format!("{}+{}", cp, off)
            } else {
                format!("{}{}", cp, off)
            }
        } else if coding_end.is_some() && cdna > coding_end.unwrap() {
            let u = cdna - coding_end.unwrap();
            if off > 0 {
                format!("*{}+{}", u, off)
            } else {
                format!("*{}{}", u, off)
            }
        } else if off > 0 {
            format!("{}+{}", cp, off)
        } else {
            format!("{}{}", cp, off)
        }
    };

    if start_offset == end_offset {
        let pos = build_pos(nearest_exon_cdna_pos, start_offset);
        Some(format!("{}{}dup", prefix, pos))
    } else {
        let start_pos = build_pos(nearest_exon_cdna_pos, start_offset);
        let end_pos = build_pos(nearest_exon_cdna_pos, end_offset);
        Some(format!("{}{}_{}dup", prefix, start_pos, end_pos))
    }
}

/// Convert intronic insertion to dup notation with explicit start/end offsets (non-coding).
pub fn convert_ins_to_dup_range_noncoding(
    hgvsc: &str,
    start_offset: i64,
    end_offset: i64,
    nearest_exon_cdna_pos: u64,
) -> Option<String> {
    let prefix_end = hgvsc
        .find(":n.")
        .map(|i| i + 3)
        .or_else(|| hgvsc.find(":c.").map(|i| i + 3))?;
    let prefix = &hgvsc[..prefix_end];

    let build_pos = |off: i64| -> String {
        if off > 0 {
            format!("{}+{}", nearest_exon_cdna_pos, off)
        } else {
            format!("{}{}", nearest_exon_cdna_pos, off)
        }
    };

    if start_offset == end_offset {
        let pos = build_pos(start_offset);
        Some(format!("{}{}dup", prefix, pos))
    } else {
        let start_pos = build_pos(start_offset);
        let end_pos = build_pos(end_offset);
        Some(format!("{}{}_{}dup", prefix, start_pos, end_pos))
    }
}

/// Convert an intronic insertion HGVSc to duplication notation (non-coding transcript).
pub fn convert_ins_to_dup_noncoding(
    hgvsc: &str,
    intron_offset: i64,
    ins_len: u64,
    nearest_exon_cdna_pos: u64,
) -> Option<String> {
    let prefix_end = hgvsc
        .find(":n.")
        .map(|i| i + 3)
        .or_else(|| hgvsc.find(":c.").map(|i| i + 3))?;
    let prefix = &hgvsc[..prefix_end];

    let build_pos = |off: i64| -> String {
        if off > 0 {
            format!("{}+{}", nearest_exon_cdna_pos, off)
        } else {
            format!("{}{}", nearest_exon_cdna_pos, off)
        }
    };

    if ins_len == 1 {
        let pos = build_pos(intron_offset);
        Some(format!("{}{}dup", prefix, pos))
    } else {
        let start_offset = intron_offset - ins_len as i64 + 1;
        let start_pos = build_pos(start_offset);
        let end_pos = build_pos(intron_offset);
        Some(format!("{}{}_{}dup", prefix, start_pos, end_pos))
    }
}

/// 3' shift an intronic indel along the transcript direction.
///
/// HGVS requires variants to be described at the most 3' position.
/// For intronic deletions and insertions/dups in repetitive regions,
/// the position must be shifted toward the 3' end of the transcript.
///
/// Returns the shifted genomic start and end positions.
pub fn three_prime_shift_intronic(
    seq_provider: &dyn SequenceProvider,
    chrom: &str,
    start: u64,
    end: u64,
    ref_allele: &fastvep_core::Allele,
    alt_allele: &fastvep_core::Allele,
    strand: fastvep_core::Strand,
    intron_genomic_start: u64,
    intron_genomic_end: u64,
) -> (u64, u64) {
    use fastvep_core::Allele;

    match (ref_allele, alt_allele) {
        // Deletion: shift the deleted bases toward 3' end
        (Allele::Sequence(ref_bases), Allele::Deletion) if !ref_bases.is_empty() => {
            let mut s = start;
            let mut e = end;

            match strand {
                fastvep_core::Strand::Forward => {
                    loop {
                        let next_pos = e + 1;
                        if next_pos > intron_genomic_end {
                            break;
                        }
                        let next_base = match seq_provider.fetch_sequence(chrom, next_pos, next_pos)
                        {
                            Ok(seq) if seq.len() == 1 => seq[0].to_ascii_uppercase(),
                            _ => break,
                        };
                        let first_base = match seq_provider.fetch_sequence(chrom, s, s) {
                            Ok(seq) if seq.len() == 1 => seq[0].to_ascii_uppercase(),
                            _ => break,
                        };
                        if next_base == first_base {
                            s += 1;
                            e += 1;
                        } else {
                            break;
                        }
                    }
                }
                fastvep_core::Strand::Reverse => {
                    loop {
                        if s == 0 || s - 1 < intron_genomic_start {
                            break;
                        }
                        let prev_pos = s - 1;
                        let prev_base =
                            match seq_provider.fetch_sequence(chrom, prev_pos, prev_pos) {
                                Ok(seq) if seq.len() == 1 => seq[0].to_ascii_uppercase(),
                                _ => break,
                            };
                        let last_base = match seq_provider.fetch_sequence(chrom, e, e) {
                            Ok(seq) if seq.len() == 1 => seq[0].to_ascii_uppercase(),
                            _ => break,
                        };
                        if prev_base == last_base {
                            s -= 1;
                            e -= 1;
                        } else {
                            break;
                        }
                    }
                }
            }
            (s, e)
        }
        // Insertion/dup: shift toward 3' end using the actual inserted bases
        (Allele::Deletion, Allele::Sequence(ins_bases)) if !ins_bases.is_empty() => {
            let ins_len = ins_bases.len();
            let mut pos = start;
            let genomic_ins: Vec<u8> =
                ins_bases.iter().map(|b| b.to_ascii_uppercase()).collect();

            match strand {
                fastvep_core::Strand::Forward => {
                    let mut shift_count = 0u64;
                    loop {
                        if pos > intron_genomic_end {
                            break;
                        }
                        let check_base =
                            match seq_provider.fetch_sequence(chrom, pos, pos) {
                                Ok(seq) if seq.len() == 1 => seq[0].to_ascii_uppercase(),
                                _ => break,
                            };
                        let idx = (shift_count as usize) % ins_len;
                        if check_base == genomic_ins[idx] {
                            pos += 1;
                            shift_count += 1;
                        } else {
                            break;
                        }
                    }
                }
                fastvep_core::Strand::Reverse => {
                    let mut shift_count = 0u64;
                    loop {
                        if pos == 0 || pos - 1 < intron_genomic_start {
                            break;
                        }
                        let check_pos = pos - 1;
                        let check_base =
                            match seq_provider.fetch_sequence(chrom, check_pos, check_pos) {
                                Ok(seq) if seq.len() == 1 => seq[0].to_ascii_uppercase(),
                                _ => break,
                            };
                        let idx = ins_len - 1 - (shift_count as usize % ins_len);
                        if check_base == genomic_ins[idx] {
                            pos -= 1;
                            shift_count += 1;
                        } else {
                            break;
                        }
                    }
                }
            }
            (pos, pos.saturating_sub(1))
        }
        _ => (start, end),
    }
}
