//! The intronic HGVSc path: 3'-shifting, and naming a shifted insertion as the
//! duplication it is.
//!
//! [`hgvsc_intronic_shifted`] is the entry point, and both annotation loops call
//! it. They used to carry a copy each - the CLI's shifted, the library's did
//! not - so the same intronic duplication came out normalised from
//! `fastvep annotate` and unnormalised from the server.

use fastvep_cache::providers::SequenceProvider;

/// Convert intronic insertion to dup notation with explicit start/end positions
/// (coding).
///
/// Each end carries its own `(exon anchor, offset)` pair. A span deep inside one
/// intron shares an anchor, but one running past the intron's midpoint is
/// written from the exon on either side - `c.5044+27_5045-47dup` - so the two
/// ends are not interchangeable.
pub fn convert_ins_to_dup_range(
    hgvsc: &str,
    start: (u64, i64),
    end: (u64, i64),
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
        // An offset of 0 is an exonic end, written as the anchor alone. Folding
        // it into the negative arm renders `c.21` as `c.210`.
        let anchor = match coding_end.filter(|&ce| cp >= 0 && cdna > ce) {
            Some(ce) => format!("*{}", cdna - ce),
            None => format!("{}", cp),
        };
        match off.cmp(&0) {
            std::cmp::Ordering::Greater => format!("{}+{}", anchor, off),
            std::cmp::Ordering::Less => format!("{}{}", anchor, off),
            std::cmp::Ordering::Equal => anchor,
        }
    };

    if start == end {
        Some(format!("{}{}dup", prefix, build_pos(start.0, start.1)))
    } else {
        Some(format!(
            "{}{}_{}dup",
            prefix,
            build_pos(start.0, start.1),
            build_pos(end.0, end.1)
        ))
    }
}

/// Convert intronic insertion to dup notation with explicit start/end positions
/// (non-coding). See [`convert_ins_to_dup_range`] for why each end carries its
/// own anchor.
pub fn convert_ins_to_dup_range_noncoding(
    hgvsc: &str,
    start: (u64, i64),
    end: (u64, i64),
) -> Option<String> {
    let prefix_end = hgvsc
        .find(":n.")
        .map(|i| i + 3)
        .or_else(|| hgvsc.find(":c.").map(|i| i + 3))?;
    let prefix = &hgvsc[..prefix_end];

    let build_pos = |cdna: u64, off: i64| -> String {
        match off.cmp(&0) {
            std::cmp::Ordering::Greater => format!("{}+{}", cdna, off),
            std::cmp::Ordering::Less => format!("{}{}", cdna, off),
            std::cmp::Ordering::Equal => format!("{}", cdna),
        }
    };

    if start == end {
        Some(format!("{}{}dup", prefix, build_pos(start.0, start.1)))
    } else {
        Some(format!(
            "{}{}_{}dup",
            prefix,
            build_pos(start.0, start.1),
            build_pos(end.0, end.1)
        ))
    }
}

/// 3' shift an intronic indel along the transcript direction.
///
/// HGVS requires variants to be described at the most 3' position.
/// For intronic deletions and insertions/dups in repetitive regions,
/// the position must be shifted toward the 3' end of the transcript.
///
/// Returns the shifted genomic start and end positions.
// Each argument is an independent coordinate, allele or flag with no
// natural grouping; bundling them into a struct would only move the
// argument list to the call site.
#[allow(clippy::too_many_arguments)]
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
                fastvep_core::Strand::Forward => loop {
                    let next_pos = e + 1;
                    if next_pos > intron_genomic_end {
                        break;
                    }
                    let next_base = match seq_provider.fetch_sequence(chrom, next_pos, next_pos) {
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
                },
                fastvep_core::Strand::Reverse => loop {
                    if s == 0 || s - 1 < intron_genomic_start {
                        break;
                    }
                    let prev_pos = s - 1;
                    let prev_base = match seq_provider.fetch_sequence(chrom, prev_pos, prev_pos) {
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
                },
            }
            (s, e)
        }
        // Insertion/dup: shift toward 3' end using the actual inserted bases
        (Allele::Deletion, Allele::Sequence(ins_bases)) if !ins_bases.is_empty() => {
            let ins_len = ins_bases.len();
            let mut pos = start;
            let genomic_ins: Vec<u8> = ins_bases.iter().map(|b| b.to_ascii_uppercase()).collect();

            match strand {
                fastvep_core::Strand::Forward => {
                    let mut shift_count = 0u64;
                    loop {
                        if pos > intron_genomic_end {
                            break;
                        }
                        let check_base = match seq_provider.fetch_sequence(chrom, pos, pos) {
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

/// The block a 3'-shifted intronic insertion duplicates, in genomic coordinates.
///
/// After a maximal 3'-shift the duplicated copy can only sit immediately 5' of
/// the insertion point *in transcript orientation*; anything further 3' would
/// have been shifted over. The inserted string rotates one base per position
/// shifted, so the block is compared against that rotated form and not against
/// the string the VCF carried.
///
/// `shifted_start` follows the same convention [`three_prime_shift_intronic`]
/// returns: the insertion sits between genomic `shifted_start - 1` and
/// `shifted_start`.
///
/// The dup anchor used to be derived by re-shifting the *unshifted* insertion
/// point through `three_prime_shift_intronic` over a single position, which
/// walks while the next base repeats the current one - a homopolymer test. A
/// `TG` insertion in a `TGTGTG…` repeat therefore never moved at all and named
/// the copy 16 bases 5' of the one HGVS asks for. That accounted for 1,390 of
/// the 1,538 HGVSc rows disagreeing with real Ensembl VEP 115.1 over a
/// genome-wide HG002 sample, and it is why the two tools split by strand: this
/// walk stays where the VCF put the variant, so fastVEP's anchor was always the
/// genomically-left one whichever way the transcript ran.
pub fn intronic_dup_span(
    seq_provider: &dyn SequenceProvider,
    chrom: &str,
    shifted_start: u64,
    ins_bases: &[u8],
    shift: u64,
    strand: fastvep_core::Strand,
) -> Option<(u64, u64)> {
    let len = ins_bases.len();
    if len == 0 {
        return None;
    }
    // One base of rotation per position travelled, in the direction of travel:
    // shifting 3' on the forward strand walks the first base to the end, and on
    // the reverse strand - where 3' runs towards lower coordinates - the last
    // base walks to the front.
    let rot = (shift % len as u64) as usize;
    let mut rotated = ins_bases.to_vec();
    match strand {
        fastvep_core::Strand::Forward => rotated.rotate_left(rot),
        fastvep_core::Strand::Reverse => rotated.rotate_right(rot),
    }

    // 5' of the insertion point is genomically before it on the forward strand
    // and after it on the reverse.
    let (lo, hi) = match strand {
        fastvep_core::Strand::Forward => {
            let hi = shifted_start.checked_sub(1)?;
            (hi.checked_sub(len as u64 - 1)?, hi)
        }
        fastvep_core::Strand::Reverse => (shifted_start, shifted_start + len as u64 - 1),
    };
    if lo == 0 {
        return None;
    }

    let block = seq_provider.fetch_sequence_slice(chrom, lo, hi).ok()?;
    if block.len() != len {
        return None;
    }
    block
        .iter()
        .zip(rotated.iter())
        .all(|(a, b)| a.eq_ignore_ascii_case(b))
        .then_some((lo, hi))
}

/// Rewrite a 3'-shifted intronic insertion as a duplication, when it is one.
///
/// `hgvsc` is the insertion notation already built for the shifted position, and
/// is returned rewritten. `None` means the insertion does not duplicate an
/// adjacent block, or the block leaves the intron it started in - a range
/// written across that boundary names bases in the wrong intron, so the
/// insertion notation it already carries is the correct description.
///
/// `coding_start` is `None` for a transcript numbered from its first base, which
/// selects `n.` numbering.
// The arguments are independent coordinates, alleles and transcript state with
// no natural grouping; a struct would only move the list to the call site.
#[allow(clippy::too_many_arguments)]
pub fn intronic_ins_as_dup(
    seq_provider: &dyn SequenceProvider,
    chrom: &str,
    transcript: &fastvep_genome::Transcript,
    hgvsc: &str,
    shifted_start: u64,
    ins_bases: &[u8],
    shift: u64,
    coding_start: Option<u64>,
    coding_end: Option<u64>,
) -> Option<String> {
    let (lo, hi) = intronic_dup_span(
        seq_provider,
        chrom,
        shifted_start,
        ins_bases,
        shift,
        transcript.strand,
    )?;
    // Transcript order, not genomic: on the reverse strand the block's 5' end is
    // its higher coordinate.
    let (first, last) = match transcript.strand {
        fastvep_core::Strand::Forward => (lo, hi),
        fastvep_core::Strand::Reverse => (hi, lo),
    };
    let start = transcript.genomic_to_intronic_cdna(first)?;
    let end = transcript.genomic_to_intronic_cdna(last)?;
    // Both ends must lie in the *same* intron. Offsets count from their own
    // exon and do not run through zero, so a block reaching into the next
    // intron - or into an exon - cannot be written as one range: stepping back
    // from `+1` does not arrive at `-2`, it names bases in the intron before.
    // Writing the crossing range anyway put PVS1's offset gate at `-2` for a
    // duplication sitting on the donor and called two ClinVar-benign MSH6 and
    // DSP variants likely pathogenic.
    //
    // A span running past the intron's own midpoint is fine, and Ensembl writes
    // it from the exon on either side: `c.5044+27_5045-47dup`. Requiring one
    // shared anchor instead of one shared intron left 160 such rows as `ins`.
    if transcript.intron_bounds_at(first) != transcript.intron_bounds_at(last) {
        return None;
    }
    match coding_start {
        Some(cs) => convert_ins_to_dup_range(hgvsc, start, end, cs, coding_end),
        None => convert_ins_to_dup_range_noncoding(hgvsc, start, end),
    }
}

/// Build the HGVSc for a variant reaching into an intron: 3'-shifted, and
/// written as a duplication where the shifted insertion sits against the block
/// it copies.
///
/// Both annotation loops call this. They had drifted: the CLI's copy shifted and
/// converted to `dup`, the library's did neither, so the same intronic
/// duplication came out normalised from `fastvep annotate` and unnormalised from
/// the server and the web UI.
///
/// `genomic_ref` and `genomic_alt` are the alleles as the VCF carried them;
/// `hgvs_ref` and `hgvs_alt` are the same pair in transcript orientation. The
/// shift reads the reference, so it needs the genomic pair; the notation is
/// written from the transcript pair.
///
/// `coding_start` is `None` for a transcript numbered from its first base, which
/// selects `n.` numbering.
// The arguments are independent coordinates, alleles and transcript state with
// no natural grouping; a struct would only move the list to the call site.
#[allow(clippy::too_many_arguments)]
pub fn hgvsc_intronic_shifted(
    seq_provider: Option<&dyn SequenceProvider>,
    chrom: &str,
    transcript: &fastvep_genome::Transcript,
    versioned_tid: &str,
    var_start: u64,
    var_end: u64,
    genomic_ref: &fastvep_core::Allele,
    genomic_alt: &fastvep_core::Allele,
    hgvs_ref: &fastvep_core::Allele,
    hgvs_alt: &fastvep_core::Allele,
    coding_start: Option<u64>,
    coding_end: Option<u64>,
) -> Option<String> {
    use fastvep_core::{Allele, Strand};

    let is_insertion = matches!(
        (hgvs_ref, hgvs_alt),
        (Allele::Deletion, Allele::Sequence(_))
    );
    let is_indel = is_insertion
        || matches!(
            (hgvs_ref, hgvs_alt),
            (Allele::Sequence(_), Allele::Deletion)
        );

    // The walk is bounded by the intron the variant sits in, and an insertion
    // sits *between* two bases, so `var_start` alone does not always find it:
    // one written against the last intronic base has `var_start` in the exon
    // beyond it, and a reverse-strand transcript travels 3' back down into that
    // intron. Reading `var_end` when `var_start` lands outside places it.
    let bounds = transcript
        .intron_bounds_at(var_start)
        .or_else(|| transcript.intron_bounds_at(var_end));
    let (shifted_start, shifted_end) = match seq_provider.filter(|_| is_indel).zip(bounds) {
        Some((sp, (intron_start, intron_end))) => three_prime_shift_intronic(
            sp,
            chrom,
            var_start,
            var_end,
            genomic_ref,
            genomic_alt,
            transcript.strand,
            intron_start,
            intron_end,
        ),
        None => (var_start, var_end),
    };
    // Distance travelled, which is also how far the inserted string rotated.
    let shift = match transcript.strand {
        Strand::Forward => shifted_start.saturating_sub(var_start),
        Strand::Reverse => var_start.saturating_sub(shifted_start),
    };
    // `hgvs_alt` is already in transcript orientation, where the shift always
    // travels 3' whichever way the transcript runs, so the string always rotates
    // left - the strand is spent before this point, on reverse-complementing it.
    //
    // Rotating right on the reverse strand instead, which is what this did,
    // named the right position with the wrong bases: `c.21-21_21-20insGTC` for a
    // variant the same transcript calls `insCGT` when the VCF spells it one base
    // over. Genomically that rotation is correct and [`intronic_dup_span`] keeps
    // it, because that function reads the reference rather than the transcript.
    let shifted_alt = match hgvs_alt {
        Allele::Sequence(ins) if is_insertion && shift > 0 && !ins.is_empty() => {
            let mut rotated = ins.clone();
            let k = (shift % rotated.len() as u64) as usize;
            rotated.rotate_left(k);
            Allele::Sequence(rotated)
        }
        other => other.clone(),
    };

    // An insertion is written over the two bases it sits between, and both are
    // mapped: `shifted_end` and `shifted_start` are the pair, in transcript
    // order on the forward strand and reversed on the reverse.
    //
    // Letting the renderer infer the second coordinate as `offset + 1` instead
    // breaks across the middle of an intron, where `+n` counts from one exon and
    // `-m` from the next: an insertion between `c.20+30` and `c.21-30` came out
    // `c.20+30_20+31ins…`, naming a base past the half the donor-side offsets
    // reach.
    let (anchor_pos, second_pos) = if is_insertion {
        match transcript.strand {
            Strand::Forward => (shifted_end, Some(shifted_end + 1)),
            Strand::Reverse => (shifted_end + 1, Some(shifted_end)),
        }
    } else {
        (
            shifted_start,
            (shifted_start != shifted_end).then_some(shifted_end),
        )
    };
    let (cdna_pos, offset) = crate::intronic_or_exonic_cdna(transcript, anchor_pos)?;
    let (end_cdna, end_offset) = second_pos
        .and_then(|p| crate::intronic_or_exonic_cdna(transcript, p))
        .map(|(c, o)| (Some(c), Some(o)))
        .unwrap_or((None, None));

    let hgvsc = match coding_start {
        Some(cs) => fastvep_hgvs::hgvsc_intronic_range(
            versioned_tid,
            cdna_pos,
            offset,
            end_cdna,
            end_offset,
            hgvs_ref,
            &shifted_alt,
            cs,
            coding_end,
        ),
        None => fastvep_hgvs::hgvsc_noncoding_intronic_range(
            versioned_tid,
            cdna_pos,
            offset,
            end_cdna,
            end_offset,
            hgvs_ref,
            &shifted_alt,
        ),
    }?;

    if let (true, Allele::Sequence(ins), Some(sp)) = (is_insertion, genomic_alt, seq_provider) {
        if hgvsc.contains("ins") && !ins.is_empty() {
            if let Some(dup) = intronic_ins_as_dup(
                sp,
                chrom,
                transcript,
                &hgvsc,
                shifted_start,
                ins,
                shift,
                coding_start,
                coding_end,
            ) {
                return Some(dup);
            }
        }
    }
    Some(hgvsc)
}

#[cfg(test)]
mod tests {
    use super::*;
    use anyhow::{anyhow, Result};
    use fastvep_core::Strand;
    use fastvep_genome::{Exon, Gene, Transcript};

    /// Minimal `SequenceProvider` over a 1-based reference string for one contig,
    /// mirroring the real readers' contract: 1-based inclusive, `Err` past the end.
    struct StrRef(&'static str);
    impl SequenceProvider for StrRef {
        fn fetch_sequence(&self, _chrom: &str, start: u64, end: u64) -> Result<Vec<u8>> {
            if start < 1 || end < start {
                return Err(anyhow!("bad range"));
            }
            let b = self.0.as_bytes();
            let s0 = (start - 1) as usize;
            if s0 >= b.len() {
                return Err(anyhow!("past contig end"));
            }
            Ok(b[s0..(end as usize).min(b.len())].to_vec())
        }
    }

    /// Two exons on `strand` with one intron between them, so an intronic
    /// position has an anchor on either side. Exon 1 is 1..=20, exon 2 is
    /// 81..=100, and the intron is 21..=80.
    fn transcript(strand: Strand) -> Transcript {
        let exon = |start: u64, end: u64, rank: u32| Exon {
            stable_id: format!("ENSE{}", rank),
            start,
            end,
            strand,
            phase: 0,
            end_phase: 0,
            rank,
        };
        Transcript {
            stable_id: "ENST00000000001".into(),
            version: Some(1),
            gene: Gene {
                stable_id: "ENSG00000000001".into(),
                symbol: Some("TEST".into()),
                symbol_source: None,
                hgnc_id: None,
                biotype: "protein_coding".into(),
                chromosome: "1".into(),
                start: 1,
                end: 100,
                strand,
            },
            biotype: "protein_coding".into(),
            chromosome: "1".into(),
            start: 1,
            end: 100,
            strand,
            exons: vec![exon(1, 20, 1), exon(81, 100, 2)],
            translation: None,
            cdna_coding_start: Some(1),
            cdna_coding_end: Some(40),
            coding_region_start: None,
            coding_region_end: None,
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

    /// 20 exonic bases, then a `TG` repeat filling the intron, then exon 2. An
    /// insertion of `TG` anywhere in that repeat is the same variant.
    const TG_REPEAT: &str = "AAAAAAAAAAAAAAAAAAAA\
                             TGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTG\
                             CCCCCCCCCCCCCCCCCCCC";

    /// The duplicated block sits immediately 5' of the *shifted* insertion
    /// point, so a `TG` insertion that travelled the whole repeat names the last
    /// two bases of it, not the two it was written against.
    #[test]
    fn dup_span_follows_the_shifted_insertion_not_the_vcf_position() {
        let r = StrRef(TG_REPEAT);
        // Insertion point right after exon 1; the repeat runs 21..=80, so the
        // maximal 3' shift on the forward strand travels its full 60 bases.
        let span = intronic_dup_span(&r, "1", 81, b"TG", 60, Strand::Forward);
        assert_eq!(span, Some((79, 80)));
    }

    /// The inserted string rotates one base per position travelled. After an odd
    /// shift of a two-base insert the block is `GT`, and reading the unrotated
    /// `TG` against the reference would reject a duplication that is real.
    #[test]
    fn dup_span_compares_against_the_rotated_insert() {
        let r = StrRef(TG_REPEAT);
        // An odd shift lands the insertion point between a G and a T.
        let span = intronic_dup_span(&r, "1", 80, b"TG", 59, Strand::Forward);
        assert_eq!(span, Some((78, 79)));
        let block = r.fetch_sequence("1", 78, 79).unwrap();
        assert_eq!(&block, b"GT", "the block is the rotated form, not `TG`");
    }

    /// 3' on the reverse strand runs towards lower coordinates, so the block a
    /// reverse-strand insertion duplicates lies *after* the insertion point.
    #[test]
    fn dup_span_reads_the_other_side_on_the_reverse_strand() {
        let r = StrRef(TG_REPEAT);
        let span = intronic_dup_span(&r, "1", 21, b"TG", 0, Strand::Reverse);
        assert_eq!(span, Some((21, 22)));
    }

    /// A one-base walk through a homopolymer is not the repeat test: this is the
    /// shape the old dup anchor used, and it stopped after one base of `TGTG…`.
    #[test]
    fn dup_span_is_none_when_the_block_does_not_repeat() {
        let r = StrRef(TG_REPEAT);
        // Inside the exon's poly-A run, a `TG` insert duplicates nothing.
        assert_eq!(
            intronic_dup_span(&r, "1", 15, b"TG", 0, Strand::Forward),
            None
        );
    }

    #[test]
    fn dup_span_declines_rather_than_reading_past_the_contig() {
        let r = StrRef(TG_REPEAT);
        assert_eq!(
            intronic_dup_span(&r, "1", 1, b"TG", 0, Strand::Forward),
            None
        );
        assert_eq!(
            intronic_dup_span(&r, "1", 0, b"TG", 0, Strand::Forward),
            None
        );
        assert_eq!(
            intronic_dup_span(&r, "1", 81, b"", 0, Strand::Forward),
            None
        );
    }

    /// The whole rewrite: a shifted intronic insertion becomes a `dup` naming the
    /// block at its 3' end, written from the anchor its own side of the intron.
    #[test]
    fn intronic_insertion_becomes_a_dup_at_its_shifted_position() {
        let r = StrRef(TG_REPEAT);
        let tr = transcript(Strand::Forward);
        let out = intronic_ins_as_dup(
            &r,
            "1",
            &tr,
            "ENST00000000001.1:c.20+59_20+60insTG",
            81,
            b"TG",
            60,
            Some(1),
            Some(40),
        );
        // 79 and 80 are the last two intronic bases, one and two before exon 2,
        // so they are written from the *downstream* anchor: c.21-2_21-1.
        assert_eq!(out.as_deref(), Some("ENST00000000001.1:c.21-2_21-1dup"));
    }

    /// A block reaching out of the intron cannot be written as one offset range -
    /// offsets count from their own exon and do not run through zero - so the
    /// insertion notation it already carries is kept.
    #[test]
    fn a_block_leaving_the_intron_keeps_its_insertion_notation() {
        // Exon 2 begins with `C`s, so a `CC` insert at the intron's 3' edge
        // duplicates a block that starts inside the exon.
        let r = StrRef(TG_REPEAT);
        let tr = transcript(Strand::Reverse);
        let out = intronic_ins_as_dup(
            &r,
            "1",
            &tr,
            "ENST00000000001.1:c.21-1_21insCC",
            81,
            b"CC",
            0,
            Some(1),
            Some(40),
        );
        assert_eq!(out, None);
    }

    /// The property the whole shift exists to provide, and the one that broke:
    /// every spelling of the same insertion inside a repeat must name the same
    /// block. Walking every position of the `TG` repeat, each with the rotation
    /// that spelling carries, must land on one answer.
    ///
    /// Against real Ensembl VEP 115.1 this was 0 of 168 transcript rows agreeing
    /// across the two spellings of six genome-wide variants; VEP agreed on all
    /// 168.
    #[test]
    fn every_spelling_of_one_insertion_names_the_same_block() {
        let r = StrRef(TG_REPEAT);
        // The repeat runs 21..=80. An insertion written at offset `k` into it
        // carries the insert rotated by `k`, and has `60 - k` left to travel.
        let answers: std::collections::HashSet<_> = (0..60)
            .map(|k| {
                let mut ins = b"TG".to_vec();
                ins.rotate_left(k % 2);
                intronic_dup_span(&r, "1", 81, &ins, (60 - k) as u64, Strand::Forward)
            })
            .collect();
        assert_eq!(
            answers,
            [Some((79, 80))].into_iter().collect(),
            "every spelling must name the last two bases of the repeat"
        );
    }

    /// A span running past the intron's own midpoint is still one intron, and is
    /// written from the exon on either side. Requiring one shared anchor instead
    /// left 160 such rows as `ins` where Ensembl writes `c.5044+27_5045-47dup`.
    #[test]
    fn a_span_crossing_the_intron_midpoint_is_still_one_dup() {
        // A 60-base insert filling the whole intron duplicates all of it.
        let r = StrRef(TG_REPEAT);
        let tr = transcript(Strand::Forward);
        let ins: Vec<u8> = TG_REPEAT.as_bytes()[20..80].to_vec();
        let out = intronic_ins_as_dup(
            &r,
            "1",
            &tr,
            "ENST00000000001.1:c.20+60_21-0ins…",
            81,
            &ins,
            0,
            Some(1),
            Some(40),
        );
        assert_eq!(out.as_deref(), Some("ENST00000000001.1:c.20+1_21-1dup"));
    }
}
