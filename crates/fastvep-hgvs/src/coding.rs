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

/// Whether `alt` is the reverse complement of `ref_bases`, which is what makes
/// an equal-length replacement an inversion rather than a delins.
fn is_reverse_complement(ref_bases: &[u8], alt: &[u8]) -> bool {
    ref_bases.len() == alt.len()
        && ref_bases
            .iter()
            .rev()
            .zip(alt.iter())
            .all(|(&r, &a)| complement(r) == a.to_ascii_uppercase())
}

fn complement(base: u8) -> u8 {
    match base.to_ascii_uppercase() {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        other => other,
    }
}

/// One HGVS change on a spliced transcript: the shape, and the cDNA span it
/// applies to.
///
/// Coding and non-coding transcripts differ only in how the coordinates are
/// *printed* - `c.` numbering counts from the initiator and the terminator,
/// `n.` numbering counts from the transcript's first base - so everything that
/// decides which shape a variant has belongs here rather than in either
/// renderer. It used to live in the coding one alone, which is why a non-coding
/// insertion was written `n.218_219insA` where VEP writes `n.217_218insA`, and
/// never collapsed to `dup`: 4,688 rows over a 6,600-variant ClinVar sample.
enum Change {
    Substitution {
        pos: u64,
        from: u8,
        to: u8,
    },
    Deletion {
        start: u64,
        end: u64,
    },
    Duplication {
        start: u64,
        end: u64,
    },
    Insertion {
        before: u64,
        after: u64,
        bases: String,
    },
    Inversion {
        start: u64,
        end: u64,
    },
    Delins {
        start: u64,
        end: u64,
        bases: String,
    },
}

/// Read a variant on a spliced transcript as an HGVS change.
///
/// `spliced_seq` is what makes the 3'-rule available. Without it a deletion is
/// written where the VCF put it and an insertion never collapses to `dup`,
/// which is a valid description but not the one HGVS asks for.
fn describe(
    cdna_start: u64,
    cdna_end: u64,
    ref_allele: &Allele,
    alt_allele: &Allele,
    spliced_seq: Option<&str>,
) -> Option<Change> {
    // Every span below is read low-to-high. The pair arrives strand-ordered, so
    // a reverse-strand deletion has its ends the other way round, and the 3'
    // shift walks off in the wrong direction from a start that is really the
    // end. Only an insertion carries meaning in the order, and it is read
    // through `min`/`max` too.
    let (cdna_lo, cdna_hi) = (cdna_start.min(cdna_end), cdna_start.max(cdna_end));
    // cDNA numbering is 1-based, so a zero here is a coordinate that was never
    // mapped. Everything below indexes the spliced sequence from it.
    if cdna_lo == 0 {
        return None;
    }
    match (ref_allele, alt_allele) {
        (Allele::Sequence(ref_bases), Allele::Sequence(alt_bases))
            if ref_bases.len() == 1 && alt_bases.len() == 1 =>
        {
            Some(Change::Substitution {
                pos: cdna_lo,
                from: ref_bases[0],
                to: alt_bases[0],
            })
        }
        // A deletion slides 3' while the base past its end repeats the base at
        // its start, which leaves the deleted sequence unchanged.
        (Allele::Sequence(_), Allele::Deletion) => {
            let (mut start, mut end) = (cdna_lo, cdna_hi);
            if let Some(seq) = spliced_seq {
                let seq_bytes = seq.as_bytes();
                let (mut s, mut e) = ((start - 1) as usize, (end - 1) as usize);
                while e + 1 < seq_bytes.len()
                    && seq_bytes[e + 1].eq_ignore_ascii_case(&seq_bytes[s])
                {
                    s += 1;
                    e += 1;
                }
                (start, end) = (s as u64 + 1, e as u64 + 1);
            }
            Some(Change::Deletion { start, end })
        }
        // HGVS 3'-rule, then duplication. The order is the whole point: an
        // insertion inside a repeat can be written at any of several positions,
        // HGVS picks the most 3' of them, and only *there* does the inserted
        // sequence sit against the copy it duplicates. Testing for a duplication
        // at the unshifted position finds one only when the variant was already
        // written 3'-most, which is why we wrote `c.4172_4173insCACCAG` where VEP
        // writes `c.4177_4182dup` - 1,109 of 1,808 in-frame insertion rows over a
        // 6,600-variant ClinVar sample.
        //
        // The shift is a rotation: stepping the insertion point over a base that
        // repeats the next base of the inserted sequence leaves the transcript
        // unchanged, so the inserted string rotates by one each step. Ensembl
        // writes the rotated form - `c.3263_3264insCGATAGCAG` for an insertion of
        // `GATAGCAGC` shifted by eight.
        (Allele::Deletion, Allele::Sequence(alt_bases)) => {
            let (mut before, mut after) = (cdna_lo, cdna_hi);
            let len = alt_bases.len();
            let mut rotation = 0usize;
            if let Some(seq) = spliced_seq {
                let seq_bytes = seq.as_bytes();
                let mut shift = 0usize;
                while let Some(&next) = seq_bytes.get(before as usize) {
                    if !next.eq_ignore_ascii_case(&alt_bases[shift % len]) {
                        break;
                    }
                    shift += 1;
                    before += 1;
                    after += 1;
                }
                rotation = shift % len;
                // After a maximal shift the only place a duplicated copy can be
                // is immediately before the insertion point; anything after it
                // would have been shifted over.
                let end = (before - 1) as usize;
                if end + 1 >= len
                    && end < seq_bytes.len()
                    && seq_bytes[end + 1 - len..=end]
                        .iter()
                        .zip(alt_bases.iter().cycle().skip(rotation))
                        .all(|(a, b)| a.eq_ignore_ascii_case(b))
                {
                    // `before` is the last duplicated base, so the span it
                    // closes is `len` bases long ending there.
                    return Some(Change::Duplication {
                        start: before + 1 - len as u64,
                        end: before,
                    });
                }
            }
            Some(Change::Insertion {
                before,
                after,
                bases: alt_bases
                    .iter()
                    .cycle()
                    .skip(rotation)
                    .take(len)
                    .map(|&b| b as char)
                    .collect(),
            })
        }
        // An equal-length replacement by the reverse complement is an inversion,
        // and HGVS has a notation for it. Checked against real VEP 115.1 over a
        // 6,600-variant ClinVar sample: all 121 same-length multi-base variants
        // whose alternate is the reverse complement of the reference are written
        // `inv`, and no `delins` row is one.
        (Allele::Sequence(ref_bases), Allele::Sequence(alt_bases)) => {
            let (start, end) = (cdna_lo, cdna_hi);
            if ref_bases.len() > 1 && is_reverse_complement(ref_bases, alt_bases) {
                Some(Change::Inversion { start, end })
            } else {
                Some(Change::Delins {
                    start,
                    end,
                    bases: std::str::from_utf8(alt_bases).ok()?.to_string(),
                })
            }
        }
        _ => None,
    }
}

/// Generate HGVSc with 3' shifting, duplication and inversion detection.
///
/// When `spliced_seq` is provided, indels in repetitive regions are shifted to
/// the most 3' position per HGVS nomenclature standard.
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

    // Ensembl applies the CDS phase offset to HGVSc for a single-base
    // substitution and to nothing else. `hgvs_transcript`
    // (`TranscriptVariationAllele.pm` release/115) takes the position from
    // `$tv->cds_start` when `$vf->var_class eq 'SNP'`, and `cds_start` carries
    // the offset - it is why a phase-1 transcript reports `CDS_position` 2286
    // where `cDNA_position` is 2285. Every other class of change goes through
    // `_get_cDNA_position`, which computes `cdna + 1 - cdna_coding_start` and
    // has no phase term at all.
    //
    // So the asymmetry is Ensembl's. #102 unified the three copies of this
    // numbering on the grounds that they had drifted apart; they had, but this
    // particular difference was not drift, and unifying it moved 5,216 ClinVar
    // indel coordinates off VEP's by one or two bases. Only a run against real
    // VEP tells the two conventions apart - both read as reasonable from the
    // code alone, and neither produces a malformed coordinate.
    const NO_PHASE: u64 = 0;
    let span = |a: u64, b: u64, phase: u64| cds_span(a, b, coding_start, coding_end, phase);

    Some(
        match describe(cdna_start, cdna_end, ref_allele, alt_allele, spliced_seq)? {
            Change::Substitution { pos, from, to } => format!(
                "{}{}{}>{}",
                prefix,
                span(pos, pos, start_phase),
                from as char,
                to as char
            ),
            Change::Deletion { start, end } => {
                format!("{}{}del", prefix, span(start, end, NO_PHASE))
            }
            Change::Duplication { start, end } => {
                format!("{}{}dup", prefix, span(start, end, NO_PHASE))
            }
            // An insertion is written over the two bases it sits between, and
            // those two can straddle the end of the CDS, so both go through the
            // same numbering.
            Change::Insertion {
                before,
                after,
                bases,
            } => format!("{}{}ins{}", prefix, span(before, after, NO_PHASE), bases),
            Change::Inversion { start, end } => {
                format!("{}{}inv", prefix, span(start, end, NO_PHASE))
            }
            Change::Delins { start, end, bases } => {
                format!("{}{}delins{}", prefix, span(start, end, NO_PHASE), bases)
            }
        },
    )
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

/// Format one intronic position as an HGVS coordinate: an exonic anchor plus a
/// signed offset, `c.151+5` or `c.*12-3`.
///
/// The anchor's own numbering is the same three-scheme problem `cds_coord`
/// solves, minus the phase - Ensembl applies that to substitutions in the CDS
/// and nowhere else. This used to be written out four times inside
/// `hgvsc_intronic_range`, once per allele shape, which is how two of the four
/// came to disagree about the 5' UTR.
fn intronic_coord(
    cdna_pos: u64,
    offset: i64,
    coding_start: u64,
    coding_end: Option<u64>,
) -> String {
    // Offset 0 is an exonic endpoint, which is written as the anchor alone.
    // A span can have one of each - `c.764_771+9delins…` runs from the last
    // coding bases of an exon into the intron - and writing only its intronic
    // end says the change happens nine bases into the intron when it starts in
    // the exon. PVS1 reads that offset to decide whether a splice consequence
    // reached the canonical dinucleotide, and stood itself down on 15
    // ClinVar-pathogenic donor and acceptor deletions on the strength of it.
    let anchor = if let Some(ce) = coding_end.filter(|&ce| cdna_pos > ce) {
        format!("*{}", cdna_pos - ce)
    } else {
        let raw = cdna_pos as i64 - coding_start as i64 + 1;
        // There is no position 0: the base before the initiator is -1.
        format!("{}", if raw <= 0 { raw - 1 } else { raw })
    };
    match offset.cmp(&0) {
        std::cmp::Ordering::Greater => format!("{}+{}", anchor, offset),
        std::cmp::Ordering::Less => format!("{}{}", anchor, offset),
        std::cmp::Ordering::Equal => anchor,
    }
}

/// Generate HGVSc for a variant reaching into an intron, in `c.` numbering.
///
/// Each end is a `(cDNA anchor, offset)` pair, and an offset of 0 marks an
/// exonic end - the two are not interchangeable, so a span with one of each is
/// written `c.764_771+9`. The end pair is optional: a single-base change, and
/// any change the caller could not map at both ends, is written from the start
/// alone.
// Each argument is an independent coordinate or allele; grouping them into a
// struct would only move the argument list to the call site.
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
    let start = (nearest_exon_cdna_pos, intron_offset);
    let end = end_cdna_pos.zip(end_intron_offset);
    let show = |(c, o): (u64, i64)| intronic_coord(c, o, coding_start, coding_end);

    // The caller's pair is strand-ordered - it comes from the variant's genomic
    // ends - so on the reverse strand the two arrive back to front. Transcript
    // order is `(anchor, offset)` ascending: a donor-side offset is positive and
    // grows away from the exon, an acceptor-side one is negative and grows
    // towards it, and the anchor separates the two sides of the same intron.
    let span = |ref_len: usize| -> String {
        match end.filter(|_| ref_len > 1) {
            Some(e) => {
                let (lo, hi) = if start <= e { (start, e) } else { (e, start) };
                format!("{}_{}", show(lo), show(hi))
            }
            None => show(start),
        }
    };

    let notation = match (ref_allele, alt_allele) {
        (Allele::Sequence(ref_bases), Allele::Sequence(alt_bases))
            if ref_bases.len() == 1 && alt_bases.len() == 1 =>
        {
            format!(
                "{}{}{}>{}",
                prefix,
                show(start),
                ref_bases[0] as char,
                alt_bases[0] as char
            )
        }
        (Allele::Sequence(ref_bases), Allele::Deletion) => {
            format!("{}{}del", prefix, span(ref_bases.len()))
        }
        // An insertion sits between two positions and is written over both, so
        // it needs the pair whatever its length.
        (Allele::Deletion, Allele::Sequence(alt_bases)) => {
            let between = match end {
                Some(e) => {
                    let (lo, hi) = if start <= e { (start, e) } else { (e, start) };
                    format!("{}_{}", show(lo), show(hi))
                }
                // With only one end mapped, the other is the next position along.
                None => format!(
                    "{}_{}",
                    show(start),
                    show((nearest_exon_cdna_pos, intron_offset + 1))
                ),
            };
            format!(
                "{}{}ins{}",
                prefix,
                between,
                std::str::from_utf8(alt_bases).unwrap_or("?")
            )
        }
        // A replacement by the reverse complement is an inversion; anything else
        // that changes the sequence in place is a delins. Neither had an arm
        // here at all, so 5,343 intronic rows over a 6,600-variant ClinVar sample
        // carried no HGVSc where VEP writes one.
        (Allele::Sequence(ref_bases), Allele::Sequence(alt_bases)) => {
            if ref_bases.len() > 1 && is_reverse_complement(ref_bases, alt_bases) {
                format!("{}{}inv", prefix, span(ref_bases.len()))
            } else {
                format!(
                    "{}{}delins{}",
                    prefix,
                    span(ref_bases.len()),
                    std::str::from_utf8(alt_bases).unwrap_or("?")
                )
            }
        }
        _ => return None,
    };
    Some(notation)
}

/// Generate HGVSc notation for non-coding transcripts using `n.` prefix.
///
/// Non-coding transcripts (lncRNA, retained_intron, etc.) number from the
/// transcript's first base rather than from the initiator, but the change they
/// describe is the same one - see [`describe`], which both renderers share.
pub fn hgvsc_noncoding(
    transcript_id: &str,
    cdna_start: u64,
    cdna_end: u64,
    ref_allele: &Allele,
    alt_allele: &Allele,
    spliced_seq: Option<&str>,
) -> Option<String> {
    let prefix = format!("{}:n.", transcript_id);
    let span = |a: u64, b: u64| {
        if a == b {
            format!("{}", a)
        } else {
            format!("{}_{}", a.min(b), a.max(b))
        }
    };
    Some(
        match describe(cdna_start, cdna_end, ref_allele, alt_allele, spliced_seq)? {
            Change::Substitution { pos, from, to } => {
                format!("{}{}{}>{}", prefix, pos, from as char, to as char)
            }
            Change::Deletion { start, end } => format!("{}{}del", prefix, span(start, end)),
            Change::Duplication { start, end } => format!("{}{}dup", prefix, span(start, end)),
            Change::Insertion {
                before,
                after,
                bases,
            } => format!("{}{}ins{}", prefix, span(before, after), bases),
            Change::Inversion { start, end } => format!("{}{}inv", prefix, span(start, end)),
            Change::Delins { start, end, bases } => {
                format!("{}{}delins{}", prefix, span(start, end), bases)
            }
        },
    )
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

/// Generate HGVSc for a variant spanning intronic positions of a non-coding
/// transcript, in `n.` numbering.
///
/// The same shapes as [`hgvsc_intronic_range`], with the anchor printed as a
/// plain cDNA position - there is no initiator to count from.
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
    let start = (nearest_exon_cdna_pos, intron_offset);
    let end = end_cdna_pos.zip(end_intron_offset);
    let show = |(c, o): (u64, i64)| match o.cmp(&0) {
        // Offset 0 is an exonic endpoint; see `intronic_coord`.
        std::cmp::Ordering::Greater => format!("{}+{}", c, o),
        std::cmp::Ordering::Less => format!("{}{}", c, o),
        std::cmp::Ordering::Equal => format!("{}", c),
    };
    // Transcript order is `(anchor, offset)` ascending; the caller's pair is
    // strand-ordered, so on the reverse strand it arrives the other way round.
    let span = |ref_len: usize| -> String {
        match end.filter(|_| ref_len > 1) {
            Some(e) => {
                let (lo, hi) = if start <= e { (start, e) } else { (e, start) };
                format!("{}_{}", show(lo), show(hi))
            }
            None => show(start),
        }
    };

    let notation = match (ref_allele, alt_allele) {
        (Allele::Sequence(ref_bases), Allele::Sequence(alt_bases))
            if ref_bases.len() == 1 && alt_bases.len() == 1 =>
        {
            format!(
                "{}{}{}>{}",
                prefix,
                show(start),
                ref_bases[0] as char,
                alt_bases[0] as char
            )
        }
        (Allele::Sequence(ref_bases), Allele::Deletion) => {
            format!("{}{}del", prefix, span(ref_bases.len()))
        }
        (Allele::Deletion, Allele::Sequence(alt_bases)) => {
            let between = match end {
                Some(e) => {
                    let (lo, hi) = if start <= e { (start, e) } else { (e, start) };
                    format!("{}_{}", show(lo), show(hi))
                }
                None => format!(
                    "{}_{}",
                    show(start),
                    show((nearest_exon_cdna_pos, intron_offset + 1))
                ),
            };
            format!(
                "{}{}ins{}",
                prefix,
                between,
                std::str::from_utf8(alt_bases).unwrap_or("?")
            )
        }
        (Allele::Sequence(ref_bases), Allele::Sequence(alt_bases)) => {
            if ref_bases.len() > 1 && is_reverse_complement(ref_bases, alt_bases) {
                format!("{}{}inv", prefix, span(ref_bases.len()))
            } else {
                format!(
                    "{}{}delins{}",
                    prefix,
                    span(ref_bases.len()),
                    std::str::from_utf8(alt_bases).unwrap_or("?")
                )
            }
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
            None,
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

    /// Ensembl's phase offset for an incomplete first codon reaches HGVSc for a
    /// single-base substitution and for nothing else, so the same base is
    /// `c.52G>A` as a substitution and `c.50del` as a deletion.
    ///
    /// That asymmetry looks like drift and is not. `hgvs_transcript` reads the
    /// position from `$tv->cds_start` - which carries the offset - only when
    /// `$vf->var_class eq 'SNP'`, and from `_get_cDNA_position` otherwise, where
    /// the arithmetic is `cdna + 1 - cdna_coding_start` with no phase term.
    /// #102 removed the difference on the reasoning that three copies of this
    /// numbering had drifted apart, and moved 5,216 ClinVar indel coordinates
    /// off VEP's by one or two bases. Real VEP 115.1 is the arbiter: nothing in
    /// the shape of either answer says which is right.
    #[test]
    fn only_a_substitution_carries_the_cds_phase_offset() {
        let (cdna_s, cdna_e, coding_start, coding_end, phase) =
            (100u64, 100u64, 51u64, Some(1000u64), 2u64);
        let at = |reference: Allele, alternate: Allele| {
            hgvsc_with_seq(
                "ENST00000001",
                cdna_s,
                cdna_e,
                &reference,
                &alternate,
                coding_start,
                coding_end,
                None,
                phase,
            )
        };
        let g = || Allele::Sequence(b"G".to_vec());

        assert_eq!(
            at(g(), Allele::Sequence(b"A".to_vec())),
            Some("ENST00000001:c.52G>A".to_string()),
            "a substitution is numbered from cds_start, which includes the phase"
        );
        for (label, alternate) in [
            ("deletion", Allele::Deletion),
            ("delins", Allele::Sequence(b"AC".to_vec())),
        ] {
            let got = at(g(), alternate);
            assert!(
                got.as_deref().is_some_and(|d| d.contains("c.50")),
                "a {label} is numbered without the phase, expected c.50, got {got:?}"
            );
        }
        // An insertion is written over the two bases it sits between, and both
        // of those are numbered the same way. It arrives as a coordinate pair
        // rather than a single position, so it does not go through `at`.
        assert_eq!(
            hgvsc_with_seq(
                "ENST00000001",
                cdna_s,
                cdna_e + 1,
                &Allele::Deletion,
                &Allele::Sequence(b"AC".to_vec()),
                coding_start,
                coding_end,
                None,
                phase,
            ),
            Some("ENST00000001:c.50_51insAC".to_string())
        );
    }
}

#[cfg(test)]
mod shape_tests {
    use super::*;

    /// An insertion inside a repeat is written at its most 3' position, and only
    /// there does it read as a duplication.
    ///
    /// Checked against real VEP 115.1: `c.4172_4173insCACCAG` is `c.4177_4182dup`
    /// to Ensembl, and testing for the duplication before shifting found one on
    /// 699 of 1,808 in-frame insertion rows instead of 1,808.
    #[test]
    fn an_insertion_in_a_repeat_shifts_before_it_can_be_a_duplication() {
        // cDNA 1-based: 1 2 3 4 5 6 7 8 9 10
        //               A C G T T T T C A  G
        let seq = "ACGTTTTCAG";
        let ins = |before: u64, bases: &str| {
            hgvsc_with_seq(
                "T",
                before + 1,
                before,
                &Allele::Deletion,
                &Allele::Sequence(bases.as_bytes().to_vec()),
                1,
                Some(10),
                Some(seq),
                0,
            )
            .unwrap()
        };
        // Inserting a T before the run of four slides to its end, where the base
        // in front of it is the T it duplicates.
        assert_eq!(ins(4, "T"), "T:c.7dup");
        // Inserting a base that does not repeat what follows stays put.
        assert_eq!(ins(4, "A"), "T:c.4_5insA");
        // The shift rotates the inserted string, and the rotated form is what is
        // written: inserting `TC` before the run of T's is `CT` one base into
        // it, which is as far as it goes.
        assert_eq!(ins(4, "TC"), "T:c.5_6insCT");
        // A pair that does repeat the run collapses to a two-base duplication.
        assert_eq!(ins(4, "TT"), "T:c.6_7dup");
    }

    /// An equal-length replacement by the reverse complement is an inversion.
    ///
    /// Verified against real VEP 115.1 over a 6,600-variant ClinVar sample: all
    /// 121 same-length multi-base variants whose alternate is the reverse
    /// complement of the reference are written `inv`, and none of the 20,781
    /// `delins` rows is one.
    #[test]
    fn a_replacement_by_the_reverse_complement_is_an_inversion() {
        let call = |r: &str, a: &str| {
            hgvsc_with_seq(
                "T",
                3,
                4,
                &Allele::Sequence(r.as_bytes().to_vec()),
                &Allele::Sequence(a.as_bytes().to_vec()),
                1,
                Some(10),
                None,
                0,
            )
            .unwrap()
        };
        assert_eq!(call("CA", "TG"), "T:c.3_4inv");
        assert_eq!(call("CA", "GT"), "T:c.3_4delinsGT");
        // A single base is never an inversion, however it complements.
        assert_eq!(
            hgvsc_with_seq(
                "T",
                3,
                3,
                &Allele::Sequence(b"C".to_vec()),
                &Allele::Sequence(b"G".to_vec()),
                1,
                Some(10),
                None,
                0,
            )
            .unwrap(),
            "T:c.3C>G"
        );
    }

    /// A multi-base change inside an intron had no arm at all, so it was
    /// annotated with no HGVSc where VEP writes one - 5,343 rows over a
    /// 6,600-variant ClinVar sample.
    #[test]
    fn a_multi_base_intronic_change_is_written_as_a_range() {
        let call = |r: &str, a: Option<&str>| {
            hgvsc_intronic_range(
                "T",
                96,
                21764,
                Some(96),
                Some(21765),
                &Allele::Sequence(r.as_bytes().to_vec()),
                &match a {
                    Some(a) => Allele::Sequence(a.as_bytes().to_vec()),
                    None => Allele::Deletion,
                },
                1,
                Some(300),
            )
            .unwrap()
        };
        assert_eq!(call("TC", Some("CG")), "T:c.96+21764_96+21765delinsCG");
        assert_eq!(call("TC", Some("GA")), "T:c.96+21764_96+21765inv");
        assert_eq!(call("TC", None), "T:c.96+21764_96+21765del");
    }

    /// The two ends of an intronic span arrive strand-ordered, so on the reverse
    /// strand they arrive back to front. Transcript order is `(anchor, offset)`
    /// ascending on both sides of the intron.
    #[test]
    fn an_intronic_span_is_written_in_transcript_order() {
        let call = |(c1, o1): (u64, i64), (c2, o2): (u64, i64)| {
            hgvsc_intronic_range(
                "T",
                c1,
                o1,
                Some(c2),
                Some(o2),
                &Allele::Sequence(b"TC".to_vec()),
                &Allele::Deletion,
                1,
                Some(300),
            )
            .unwrap()
        };
        assert_eq!(call((96, 5), (96, 4)), "T:c.96+4_96+5del");
        assert_eq!(call((97, -5), (97, -6)), "T:c.97-6_97-5del");
        // Opposite sides of one intron: the anchor separates them.
        assert_eq!(call((97, -6), (96, 5)), "T:c.96+5_97-6del");
    }
}
