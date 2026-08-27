use fastvep_core::Allele;

/// Format one cDNA position as an HGVS `c.` coordinate.
///
/// Three numbering schemes meet on a coding transcript: `-N` counting back from
/// the initiator, `N` inside the CDS, `*N` counting on from the terminator.
/// Which one applies is a property of the *position*, not of the variant, and
/// deciding it per variant is what let a span that crosses a boundary be
/// numbered in the wrong scheme: a delins over FOXL2's stop codon came out
/// `c.1127_1135`, naming four bases past the end of an 1131-base CDS, where
/// VEP writes `c.1127_*4` (#100).
///
/// `start_phase` is Ensembl's offset for a CDS whose first codon is incomplete,
/// and applies only inside the CDS - there is no position 0, so the 5' UTR
/// numbering has neither the `+ 1` nor the phase.
fn cds_coord(
    cdna_pos: u64,
    coding_start: u64,
    coding_end: Option<u64>,
    start_phase: u64,
) -> String {
    if cdna_pos < coding_start {
        return format!("{}", cdna_pos as i64 - coding_start as i64);
    }
    if let Some(coding_end) = coding_end {
        if cdna_pos > coding_end {
            return format!("*{}", cdna_pos - coding_end);
        }
    }
    format!(
        "{}",
        cdna_pos as i64 - coding_start as i64 + 1 + start_phase as i64
    )
}

/// Format a cDNA span as an HGVS `c.` coordinate or coordinate range,
/// collapsing a single-position span to one coordinate.
fn cds_span(
    cdna_start: u64,
    cdna_end: u64,
    coding_start: u64,
    coding_end: Option<u64>,
    start_phase: u64,
) -> String {
    let start = cds_coord(cdna_start, coding_start, coding_end, start_phase);
    let end = cds_coord(cdna_end, coding_start, coding_end, start_phase);
    if start == end {
        start
    } else {
        format!("{}_{}", start, end)
    }
}

/// Generate HGVSc (coding DNA) notation.
///
/// Uses the transcript ID as the reference and cDNA position numbering.
/// Format: ENST00000001.1:c.123A>G
pub fn hgvsc(
    transcript_id: &str,
    cdna_start: u64,
    cdna_end: u64,
    ref_allele: &Allele,
    alt_allele: &Allele,
    coding_start: u64,
    coding_end: Option<u64>,
) -> Option<String> {
    hgvsc_with_seq(
        transcript_id,
        cdna_start,
        cdna_end,
        ref_allele,
        alt_allele,
        coding_start,
        coding_end,
        None,
        0,
    )
}

/// Generate HGVSc with optional 3' shifting for deletions/insertions.
///
/// When `spliced_seq` is provided, deletions in repetitive regions are shifted
/// to the most 3' position per HGVS nomenclature standard.
// Each argument is an independent coordinate, allele or flag with no
// natural grouping; bundling them into a struct would only move the
// argument list to the call site.
#[allow(clippy::too_many_arguments)]
pub fn hgvsc_with_seq(
    transcript_id: &str,
    cdna_start: u64,
    cdna_end: u64,
    ref_allele: &Allele,
    alt_allele: &Allele,
    coding_start: u64,
    coding_end: Option<u64>,
    spliced_seq: Option<&str>,
    start_phase: u64,
) -> Option<String> {
    let prefix = format!("{}:c.", transcript_id);

    let pos_str = cds_span(cdna_start, cdna_end, coding_start, coding_end, start_phase);

    let notation = match (ref_allele, alt_allele) {
        // SNV
        (Allele::Sequence(ref_bases), Allele::Sequence(alt_bases))
            if ref_bases.len() == 1 && alt_bases.len() == 1 =>
        {
            format!(
                "{}{}{}>{}",
                prefix, pos_str, ref_bases[0] as char, alt_bases[0] as char
            )
        }
        // Deletion — apply HGVS 3' shifting in repetitive regions
        (Allele::Sequence(_), Allele::Deletion) => {
            // Try 3' shifting if we have the transcript sequence
            let (shifted_start, shifted_end) = if let Some(seq) = spliced_seq {
                let seq_bytes = seq.as_bytes();
                let mut s = (cdna_start - 1) as usize; // 0-based start
                let mut e = (cdna_end - 1) as usize; // 0-based end (inclusive)
                                                     // Shift right while the base at position end+1 matches base at start
                while e + 1 < seq_bytes.len()
                    && seq_bytes[e + 1].eq_ignore_ascii_case(&seq_bytes[s])
                {
                    s += 1;
                    e += 1;
                }
                (s as u64 + 1, e as u64 + 1) // back to 1-based
            } else {
                (cdna_start, cdna_end)
            };
            let shifted_pos = cds_span(
                shifted_start,
                shifted_end,
                coding_start,
                coding_end,
                start_phase,
            );
            format!("{}{}del", prefix, shifted_pos)
        }
        // Insertion — normalize coordinates: ensure ins_before < ins_after
        (Allele::Deletion, Allele::Sequence(alt_bases)) => {
            let ins_before_cdna = cdna_start.min(cdna_end); // base before insertion
            let ins_after_cdna = cdna_start.max(cdna_end); // base after insertion

            // Check for duplication: if inserted bases match the preceding OR following sequence
            let is_dup = if let Some(seq) = spliced_seq {
                let seq_bytes = seq.as_bytes();
                let ins_len = alt_bases.len();
                let before_pos = (ins_before_cdna - 1) as usize;
                // Check preceding sequence
                let dup_before = if before_pos + 1 >= ins_len && before_pos < seq_bytes.len() {
                    let preceding = &seq_bytes[before_pos + 1 - ins_len..=before_pos];
                    preceding
                        .iter()
                        .zip(alt_bases.iter())
                        .all(|(a, b)| a.eq_ignore_ascii_case(b))
                } else {
                    false
                };
                // Check following sequence
                let dup_after = if !dup_before {
                    let after_pos = ins_after_cdna as usize; // 0-based index of base after insertion
                    if after_pos + ins_len <= seq_bytes.len() {
                        let following = &seq_bytes[after_pos..after_pos + ins_len];
                        following
                            .iter()
                            .zip(alt_bases.iter())
                            .all(|(a, b)| a.eq_ignore_ascii_case(b))
                    } else {
                        false
                    }
                } else {
                    false
                };
                dup_before || dup_after
            } else {
                false
            };

            // A duplication is written over the duplicated bases, which end at
            // the base before the insertion point. A plain insertion is written
            // over the two bases it sits between, and those two can straddle
            // the end of the CDS, so both spans go through the same numbering.
            let ins_pos_str = if is_dup {
                let ins_len = alt_bases.len() as u64;
                cds_span(
                    ins_before_cdna.saturating_sub(ins_len.saturating_sub(1)),
                    ins_before_cdna,
                    coding_start,
                    coding_end,
                    start_phase,
                )
            } else {
                cds_span(
                    ins_before_cdna,
                    ins_after_cdna,
                    coding_start,
                    coding_end,
                    start_phase,
                )
            };

            if is_dup {
                format!("{}{}dup", prefix, ins_pos_str)
            } else {
                format!(
                    "{}{}ins{}",
                    prefix,
                    ins_pos_str,
                    std::str::from_utf8(alt_bases).unwrap_or("?")
                )
            }
        }
        // MNV or complex
        (Allele::Sequence(_), Allele::Sequence(alt_bases)) => {
            format!(
                "{}{}delins{}",
                prefix,
                pos_str,
                std::str::from_utf8(alt_bases).unwrap_or("?")
            )
        }
        _ => return None,
    };

    Some(notation)
}

/// Generate HGVSc notation for an intronic variant.
///
/// Uses offset notation from the nearest exon boundary:
/// - Positive offset: c.151+5A>G (donor side)
/// - Negative offset: c.152-3A>G (acceptor side)
pub fn hgvsc_intronic(
    transcript_id: &str,
    nearest_exon_cdna_pos: u64,
    intron_offset: i64,
    ref_allele: &Allele,
    alt_allele: &Allele,
    coding_start: u64,
    coding_end: Option<u64>,
) -> Option<String> {
    hgvsc_intronic_range(
        transcript_id,
        nearest_exon_cdna_pos,
        intron_offset,
        None, // end_cdna_pos
        None, // end_offset
        ref_allele,
        alt_allele,
        coding_start,
        coding_end,
    )
}

/// Generate HGVSc intronic notation with optional end position for multi-base variants.
// Each argument is an independent coordinate, allele or flag with no
// natural grouping; bundling them into a struct would only move the
// argument list to the call site.
#[allow(clippy::too_many_arguments)]
pub fn hgvsc_intronic_range(
    transcript_id: &str,
    nearest_exon_cdna_pos: u64,
    intron_offset: i64,
    end_cdna_pos: Option<u64>,
    end_intron_offset: Option<i64>,
    ref_allele: &Allele,
    alt_allele: &Allele,
    coding_start: u64,
    coding_end: Option<u64>,
) -> Option<String> {
    let prefix = format!("{}:c.", transcript_id);

    // Convert cDNA pos to CDS-relative position.
    // For CDS positions: cds_pos = cdna - coding_start + 1 (so position 1 = first CDS base)
    // For 5'UTR positions: cds_pos = cdna - coding_start (no +1, since there's no position 0;
    //   position -1 = last 5'UTR base at cdna = coding_start - 1)
    let raw_cds_pos = nearest_exon_cdna_pos as i64 - coding_start as i64 + 1;
    let cds_pos = if raw_cds_pos <= 0 {
        raw_cds_pos - 1
    } else {
        raw_cds_pos
    };

    // Build the position string with offset
    let pos_str = if cds_pos < 0 {
        // 5' UTR
        if intron_offset > 0 {
            format!("{}+{}", cds_pos, intron_offset)
        } else {
            format!("{}{}", cds_pos, intron_offset) // offset is already negative
        }
    } else if coding_end.is_some_and(|ce| nearest_exon_cdna_pos > ce) {
        // 3' UTR
        let utr_offset = nearest_exon_cdna_pos - coding_end.unwrap();
        if intron_offset > 0 {
            format!("*{}+{}", utr_offset, intron_offset)
        } else {
            format!("*{}{}", utr_offset, intron_offset)
        }
    } else if intron_offset > 0 {
        format!("{}+{}", cds_pos, intron_offset)
    } else {
        format!("{}{}", cds_pos, intron_offset) // offset is already negative
    };

    // Format the variant
    let notation = match (ref_allele, alt_allele) {
        (Allele::Sequence(ref_bases), Allele::Sequence(alt_bases))
            if ref_bases.len() == 1 && alt_bases.len() == 1 =>
        {
            format!(
                "{}{}{}>{}",
                prefix, pos_str, ref_bases[0] as char, alt_bases[0] as char
            )
        }
        (Allele::Sequence(ref_bases), Allele::Deletion) => {
            // A single-base deletion, or one with no end offset from the
            // caller, is written from the start position alone; `filter` binds
            // the offset only in the case that actually uses it.
            if let Some(e_offset) = end_intron_offset.filter(|_| ref_bases.len() > 1) {
                // Multi-base intronic deletion: use end position from caller
                let e_cdna = end_cdna_pos.unwrap_or(nearest_exon_cdna_pos);
                let e_raw = e_cdna as i64 - coding_start as i64 + 1;
                let e_cds = if e_raw <= 0 { e_raw - 1 } else { e_raw };
                let end_pos_str = if e_cds < 0 {
                    if e_offset > 0 {
                        format!("{}+{}", e_cds, e_offset)
                    } else {
                        format!("{}{}", e_cds, e_offset)
                    }
                } else if coding_end.is_some_and(|ce| e_cdna > ce) {
                    let utr = e_cdna - coding_end.unwrap();
                    if e_offset > 0 {
                        format!("*{}+{}", utr, e_offset)
                    } else {
                        format!("*{}{}", utr, e_offset)
                    }
                } else if e_offset > 0 {
                    format!("{}+{}", e_cds, e_offset)
                } else {
                    format!("{}{}", e_cds, e_offset)
                };

                // For same-sign offsets: smaller absolute value first (closer to exon)
                // For positive: +4035 before +4038
                // For negative: -2625 before -2616
                let same_sign =
                    (intron_offset > 0 && e_offset > 0) || (intron_offset < 0 && e_offset < 0);
                let (start_str, end_str) = if same_sign && intron_offset >= e_offset {
                    (end_pos_str, pos_str.clone())
                } else {
                    (pos_str.clone(), end_pos_str)
                };
                format!("{}{}_{}del", prefix, start_str, end_str)
            } else {
                format!("{}{}del", prefix, pos_str)
            }
        }
        // Insertion in intron — show range: c.X+N_X+N+1insB
        (Allele::Deletion, Allele::Sequence(alt_bases)) => {
            let ins_str = std::str::from_utf8(alt_bases).unwrap_or("?");
            let next_offset = intron_offset + 1;
            let build_pos = |cdna: u64, off: i64| -> String {
                let raw = cdna as i64 - coding_start as i64 + 1;
                let cp = if raw <= 0 { raw - 1 } else { raw }; // skip position 0 for 5'UTR
                if cp < 0 {
                    if off > 0 {
                        format!("{}+{}", cp, off)
                    } else {
                        format!("{}{}", cp, off)
                    }
                } else if coding_end.is_some_and(|ce| cdna > ce) {
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
            let end_str = build_pos(nearest_exon_cdna_pos, next_offset);
            format!("{}{}_{}ins{}", prefix, pos_str, end_str, ins_str)
        }
        _ => return None,
    };

    Some(notation)
}

/// Generate HGVSc notation for non-coding transcripts using `n.` prefix.
///
/// Non-coding transcripts (lncRNA, retained_intron, etc.) use cDNA position
/// directly with `n.` prefix, e.g. `ENST00000472807.5:n.1234A>G`.
pub fn hgvsc_noncoding(
    transcript_id: &str,
    cdna_start: u64,
    cdna_end: u64,
    ref_allele: &Allele,
    alt_allele: &Allele,
) -> Option<String> {
    let prefix = format!("{}:n.", transcript_id);
    let (pos_min, pos_max) = if cdna_start <= cdna_end {
        (cdna_start, cdna_end)
    } else {
        (cdna_end, cdna_start)
    };
    let pos_str = if pos_min == pos_max {
        format!("{}", pos_min)
    } else {
        format!("{}_{}", pos_min, pos_max)
    };

    let notation = match (ref_allele, alt_allele) {
        (Allele::Sequence(ref_bases), Allele::Sequence(alt_bases))
            if ref_bases.len() == 1 && alt_bases.len() == 1 =>
        {
            format!(
                "{}{}{}>{}",
                prefix, pos_str, ref_bases[0] as char, alt_bases[0] as char
            )
        }
        (Allele::Sequence(_), Allele::Deletion) => {
            format!("{}{}del", prefix, pos_str)
        }
        (Allele::Deletion, Allele::Sequence(alt_bases)) => {
            let ins_pos = format!("{}_{}", pos_max, pos_max + 1);
            format!(
                "{}{}ins{}",
                prefix,
                ins_pos,
                std::str::from_utf8(alt_bases).unwrap_or("?")
            )
        }
        (Allele::Sequence(_), Allele::Sequence(alt_bases)) => {
            format!(
                "{}{}delins{}",
                prefix,
                pos_str,
                std::str::from_utf8(alt_bases).unwrap_or("?")
            )
        }
        _ => return None,
    };
    Some(notation)
}

/// Generate HGVSc intronic notation for non-coding transcripts using `n.` prefix.
pub fn hgvsc_noncoding_intronic(
    transcript_id: &str,
    nearest_exon_cdna_pos: u64,
    intron_offset: i64,
    ref_allele: &Allele,
    alt_allele: &Allele,
) -> Option<String> {
    hgvsc_noncoding_intronic_range(
        transcript_id,
        nearest_exon_cdna_pos,
        intron_offset,
        None,
        None,
        ref_allele,
        alt_allele,
    )
}

/// Generate HGVSc intronic notation for non-coding transcripts with optional end position.
pub fn hgvsc_noncoding_intronic_range(
    transcript_id: &str,
    nearest_exon_cdna_pos: u64,
    intron_offset: i64,
    end_cdna_pos: Option<u64>,
    end_intron_offset: Option<i64>,
    ref_allele: &Allele,
    alt_allele: &Allele,
) -> Option<String> {
    let prefix = format!("{}:n.", transcript_id);
    let pos_str = if intron_offset > 0 {
        format!("{}+{}", nearest_exon_cdna_pos, intron_offset)
    } else {
        format!("{}{}", nearest_exon_cdna_pos, intron_offset)
    };

    let notation = match (ref_allele, alt_allele) {
        (Allele::Sequence(ref_bases), Allele::Sequence(alt_bases))
            if ref_bases.len() == 1 && alt_bases.len() == 1 =>
        {
            format!(
                "{}{}{}>{}",
                prefix, pos_str, ref_bases[0] as char, alt_bases[0] as char
            )
        }
        (Allele::Sequence(ref_bases), Allele::Deletion) => {
            // A single-base deletion, or one with no end offset from the
            // caller, is written from the start position alone; `filter` binds
            // the offset only in the case that actually uses it.
            if let Some(e_offset) = end_intron_offset.filter(|_| ref_bases.len() > 1) {
                let e_cdna = end_cdna_pos.unwrap_or(nearest_exon_cdna_pos);
                let end_pos_str = if e_offset > 0 {
                    format!("{}+{}", e_cdna, e_offset)
                } else {
                    format!("{}{}", e_cdna, e_offset)
                };
                // Order: smaller offset first
                let (s, e) = if intron_offset < e_offset {
                    (pos_str.clone(), end_pos_str)
                } else {
                    (end_pos_str, pos_str.clone())
                };
                format!("{}{}_{}del", prefix, s, e)
            } else {
                format!("{}{}del", prefix, pos_str)
            }
        }
        // Insertion in intron — show range: n.X+N_X+N+1insB
        (Allele::Deletion, Allele::Sequence(alt_bases)) => {
            let ins_str = std::str::from_utf8(alt_bases).unwrap_or("?");
            let next_offset = intron_offset + 1;
            let end_pos_str = if next_offset > 0 {
                format!("{}+{}", nearest_exon_cdna_pos, next_offset)
            } else {
                format!("{}{}", nearest_exon_cdna_pos, next_offset)
            };
            format!("{}{}_{}ins{}", prefix, pos_str, end_pos_str, ins_str)
        }
        _ => return None,
    };
    Some(notation)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hgvsc_snv() {
        let result = hgvsc(
            "ENST00000001",
            151,
            151,
            &Allele::Sequence(b"A".to_vec()),
            &Allele::Sequence(b"G".to_vec()),
            51,
            None,
        );
        assert_eq!(result, Some("ENST00000001:c.101A>G".to_string()));
    }

    #[test]
    fn test_hgvsc_deletion() {
        let result = hgvsc(
            "ENST00000001",
            54,
            56,
            &Allele::Sequence(b"ACG".to_vec()),
            &Allele::Deletion,
            51,
            None,
        );
        assert_eq!(result, Some("ENST00000001:c.4_6del".to_string()));
    }

    #[test]
    fn test_hgvsc_insertion() {
        let result = hgvsc(
            "ENST00000001",
            54,
            53,
            &Allele::Deletion,
            &Allele::Sequence(b"TTT".to_vec()),
            51,
            None,
        );
        assert_eq!(result, Some("ENST00000001:c.3_4insTTT".to_string()));
    }

    #[test]
    fn test_hgvsc_5_utr() {
        let result = hgvsc(
            "ENST00000001",
            10,
            10,
            &Allele::Sequence(b"A".to_vec()),
            &Allele::Sequence(b"G".to_vec()),
            51,
            None,
        );
        // 5'UTR: position = cdna - coding_start = 10 - 51 = -41
        assert_eq!(result, Some("ENST00000001:c.-41A>G".to_string()));
    }

    #[test]
    fn test_hgvsc_3_utr() {
        // coding_end = 1003 in cDNA, variant at cDNA 1021
        let result = hgvsc(
            "ENST00000001",
            1021,
            1021,
            &Allele::Sequence(b"G".to_vec()),
            &Allele::Sequence(b"A".to_vec()),
            51,
            Some(1003),
        );
        // utr_offset = 1021 - 1003 = 18
        assert_eq!(result, Some("ENST00000001:c.*18G>A".to_string()));
    }

    #[test]
    fn test_hgvsc_intronic_donor() {
        // Variant 5 bases into intron, near donor (upstream exon)
        // nearest exon cDNA pos = 201 (end of exon 1), offset = +5
        // CDS pos = 201 - 51 + 1 = 151
        let result = hgvsc_intronic(
            "ENST00000001",
            201,
            5,
            &Allele::Sequence(b"G".to_vec()),
            &Allele::Sequence(b"A".to_vec()),
            51,
            None,
        );
        assert_eq!(result, Some("ENST00000001:c.151+5G>A".to_string()));
    }

    #[test]
    fn test_hgvsc_intronic_acceptor() {
        // Variant 3 bases before exon 2, near acceptor
        // nearest exon cDNA pos = 202 (start of exon 2), offset = -3
        // CDS pos = 202 - 51 + 1 = 152
        let result = hgvsc_intronic(
            "ENST00000001",
            202,
            -3,
            &Allele::Sequence(b"A".to_vec()),
            &Allele::Sequence(b"G".to_vec()),
            51,
            None,
        );
        assert_eq!(result, Some("ENST00000001:c.152-3A>G".to_string()));
    }

    #[test]
    fn test_hgvsc_deletion_3prime_shift() {
        // Sequence: AACTTTTGA at CDS positions 1-9 (coding_start=1, cdna positions 1-9)
        // Deletion of T at cdna position 4 (cds pos 4) should shift to position 7
        // because TTTT is repetitive — shift to the most 3' T
        let seq = "AACTTTTGA";
        let result = hgvsc_with_seq(
            "ENST00000001",
            4,
            4, // cdna_start=4, cdna_end=4 (the first T in TTTT)
            &Allele::Sequence(b"T".to_vec()),
            &Allele::Deletion,
            1,
            None,
            Some(seq),
            0,
        );
        // Should shift from pos 4 to pos 7 (last T in the TTTT run)
        assert_eq!(result, Some("ENST00000001:c.7del".to_string()));
    }

    #[test]
    fn test_hgvsc_deletion_no_shift() {
        // Non-repetitive: deletion at position 1 (A) — no shifting possible
        let seq = "ACGTACGT";
        let result = hgvsc_with_seq(
            "ENST00000001",
            1,
            1,
            &Allele::Sequence(b"A".to_vec()),
            &Allele::Deletion,
            1,
            None,
            Some(seq),
            0,
        );
        assert_eq!(result, Some("ENST00000001:c.1del".to_string()));
    }

    #[test]
    fn test_hgvsc_noncoding_snv() {
        let result = hgvsc_noncoding(
            "ENST00000472807.5",
            100,
            100,
            &Allele::Sequence(b"A".to_vec()),
            &Allele::Sequence(b"G".to_vec()),
        );
        assert_eq!(result, Some("ENST00000472807.5:n.100A>G".to_string()));
    }

    // ---- spans that cross a CDS boundary (#100) ----

    /// A span that starts inside the CDS and ends past the terminator numbers
    /// each end in its own scheme.
    ///
    /// It used to number both ends as CDS positions, so FOXL2's
    /// `c.1127_*4delinsCG` came out `c.1127_1135` - four bases past the end of
    /// an 1131-base CDS, a coordinate the transcript does not have.
    #[test]
    fn a_delins_reaching_past_the_terminator_numbers_its_end_from_the_stop() {
        // CDS is cDNA 51..1181 (1131 bases). cDNA 1177..1185 spans c.1127
        // to four bases into the 3' UTR.
        let result = hgvsc(
            "ENST00000648323.1",
            1177,
            1185,
            &Allele::Sequence(b"GCTCTCAGA".to_vec()),
            &Allele::Sequence(b"CG".to_vec()),
            51,
            Some(1181),
        );
        assert_eq!(
            result,
            Some("ENST00000648323.1:c.1127_*4delinsCG".to_string())
        );
    }

    /// The mirror at the other end: a span starting in the 5' UTR keeps its end
    /// coordinate, which used to be dropped entirely (`c.-4delinsTT`).
    #[test]
    fn a_delins_reaching_out_of_the_five_prime_utr_keeps_both_coordinates() {
        // CDS starts at cDNA 149. cDNA 145..153 spans c.-4 to c.5.
        let result = hgvsc(
            "ENST00000383165.4",
            145,
            153,
            &Allele::Sequence(b"GGAGATGAC".to_vec()),
            &Allele::Sequence(b"TT".to_vec()),
            149,
            Some(676),
        );
        assert_eq!(result, Some("ENST00000383165.4:c.-4_5delinsTT".to_string()));
    }

    /// A deletion whose end falls inside the CDS numbers that end as a CDS
    /// position. The 5' UTR branch used to reuse its own formula for both ends,
    /// leaving the end one too low - `c.-9_1del` for a span ending at c.2.
    #[test]
    fn a_deletion_from_the_five_prime_utr_into_the_cds_numbers_its_end_in_the_cds() {
        // CDS starts at cDNA 11; cDNA 2..12 spans c.-9 to c.2.
        let result = hgvsc(
            "ENST00000001",
            2,
            12,
            &Allele::Sequence(b"CCCATACCTCG".to_vec()),
            &Allele::Deletion,
            11,
            Some(100),
        );
        assert_eq!(result, Some("ENST00000001:c.-9_2del".to_string()));
    }

    /// There is no `c.0`, so a span reaching back across the initiator has to
    /// skip it. `c.-1_1` names two positions where three bases were changed.
    #[test]
    fn a_span_across_the_initiator_skips_the_nonexistent_position_zero() {
        // CDS starts at cDNA 11; cDNA 9..11 spans c.-2, c.-1, c.1.
        let result = hgvsc(
            "ENST00000001",
            9,
            11,
            &Allele::Sequence(b"CCA".to_vec()),
            &Allele::Deletion,
            11,
            Some(100),
        );
        assert_eq!(result, Some("ENST00000001:c.-2_1del".to_string()));
    }

    /// The last base of the CDS is a CDS position, not `*0`. The 3' UTR branch
    /// used to reach that coordinate by arithmetic and emit `c.*0_*1`.
    #[test]
    fn the_last_cds_base_is_never_numbered_star_zero() {
        // CDS is cDNA 11..40; cDNA 40..41 straddles the terminator's last base
        // and the first 3' UTR base.
        let result = hgvsc(
            "ENST00000001",
            40,
            41,
            &Allele::Deletion,
            &Allele::Sequence(b"AC".to_vec()),
            11,
            Some(40),
        );
        assert_eq!(result, Some("ENST00000001:c.30_*1insAC".to_string()));
    }

    /// A CDS position carries Ensembl's phase offset for an incomplete first
    /// codon. Deletions used to skip it while every other path applied it, so
    /// the same base was `c.716` as a substitution and `c.715` as a deletion.
    #[test]
    fn a_deletion_carries_the_same_phase_offset_as_a_substitution() {
        let args = (100u64, 100u64, 51u64, Some(1000u64), 2u64);
        let (cdna_s, cdna_e, coding_start, coding_end, phase) = args;
        let sub = hgvsc_with_seq(
            "ENST00000001",
            cdna_s,
            cdna_e,
            &Allele::Sequence(b"G".to_vec()),
            &Allele::Sequence(b"A".to_vec()),
            coding_start,
            coding_end,
            None,
            phase,
        );
        let del = hgvsc_with_seq(
            "ENST00000001",
            cdna_s,
            cdna_e,
            &Allele::Sequence(b"G".to_vec()),
            &Allele::Deletion,
            coding_start,
            coding_end,
            None,
            phase,
        );
        assert_eq!(sub, Some("ENST00000001:c.52G>A".to_string()));
        assert_eq!(del, Some("ENST00000001:c.52del".to_string()));
    }
}
