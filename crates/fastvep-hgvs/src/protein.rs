use fastvep_core::Strand;
use fastvep_genome::codon::{aa_one_to_three, CodonTable};

/// Generate HGVSp (protein) notation.
///
/// Format: ENSP00000001:p.Arg41Lys (missense)
///         ENSP00000001:p.Arg41Ter (stop gained)
///         ENSP00000001:p.Arg41= (synonymous)
///         ENSP00000001:p.Arg41fs (frameshift)
pub fn hgvsp(
    protein_id: &str,
    protein_pos: u64,
    ref_aa: u8,
    alt_aa: u8,
    is_frameshift: bool,
) -> Option<String> {
    let prefix = format!("{}:p.", protein_id);
    let ref_aa3 = aa_one_to_three(ref_aa);

    if is_frameshift {
        return Some(format!("{}{}{}fs", prefix, ref_aa3, protein_pos));
    }

    if ref_aa == alt_aa {
        // Synonymous
        return Some(format!("{}{}{}=", prefix, ref_aa3, protein_pos));
    }

    let alt_aa3 = aa_one_to_three(alt_aa);

    if alt_aa == b'*' {
        // Stop gained
        return Some(format!("{}{}{}{}", prefix, ref_aa3, protein_pos, alt_aa3));
    }

    if ref_aa == b'*' {
        // Stop lost - extension
        return Some(format!("{}{}{}ext*?", prefix, alt_aa3, protein_pos));
    }

    // Missense
    Some(format!("{}{}{}{}", prefix, ref_aa3, protein_pos, alt_aa3))
}

/// Render a 1-based inclusive residue range, e.g. `Gly41` or `Asn587_Asp600`.
fn residue_span(first_pos: u64, residues: &[u8]) -> String {
    let first = aa_one_to_three(residues[0]);
    if residues.len() == 1 {
        format!("{}{}", first, first_pos)
    } else {
        let last = aa_one_to_three(residues[residues.len() - 1]);
        format!(
            "{}{}_{}{}",
            first,
            first_pos,
            last,
            first_pos + residues.len() as u64 - 1
        )
    }
}

/// Render a change that takes out the initiation codon, whose downstream
/// consequence is unknown: `Glu1?` for a single residue, `MetAsnIle1_?3` for a
/// run. The residues are still named - the annotation does say what was
/// affected - and the `?` marks the part that cannot be resolved from sequence.
///
/// `start` is the position of the first residue named. It is not always 1: a
/// change that begins at the initiator can share its leading residue with the
/// replacement, and Ensembl trims that off first. `XP/XAKSTVGA` on
/// ENSP00000500558 is `p.Pro2?` in VEP 115.1, not `p.XaaPro1_?2`.
///
/// The run's second coordinate is written as the last residue named, which is
/// what `p.Met1_?N` reads as. Not verified: every multi-residue case in the
/// ClinVar sample starts at residue 1, where the last position and the residue
/// count are the same number, so nothing here distinguishes them.
fn uncertain_from_initiator(start: u64, residues: &[u8]) -> String {
    if residues.len() == 1 {
        format!("{}{}?", aa_one_to_three(residues[0]), start)
    } else {
        format!(
            "{}{}_?{}",
            three_letter(residues),
            start,
            start + residues.len() as u64 - 1
        )
    }
}

fn three_letter(residues: &[u8]) -> String {
    residues.iter().map(|&b| aa_one_to_three(b)).collect()
}

/// Describe the change using only the residues the caller supplied, without
/// consulting the peptide. Always available when there is at least one
/// reference residue, and the fallback whenever the peptide is missing or
/// cannot be trusted — a valid, unshifted description beats emitting nothing.
fn unshifted_description(
    prefix: &str,
    start: u64,
    reference: &[u8],
    alternate: &[u8],
) -> Option<String> {
    if reference.is_empty() {
        return None;
    }
    let range = residue_span(start, reference);
    if alternate.is_empty() {
        Some(format!("{}{}del", prefix, range))
    } else {
        Some(format!(
            "{}{}delins{}",
            prefix,
            range,
            three_letter(alternate)
        ))
    }
}

/// Whether `peptide` actually carries `reference` at `protein_start`.
///
/// A transcript whose sequence disagrees with its own coordinates would
/// otherwise produce a confident, well-formed, wrong description — worse than
/// no normalisation at all.
fn peptide_carries(peptide: &[u8], protein_start: u64, reference: &[u8]) -> bool {
    let Some(lo) = protein_start.checked_sub(1).map(|v| v as usize) else {
        return false;
    };
    if reference.is_empty() {
        // Pure insertion: only the flanking position is read.
        return lo <= peptide.len();
    }
    peptide.get(lo..lo + reference.len()) == Some(reference)
}

/// Whether `protein_start` names the *end* of the affected span rather than its
/// start.
///
/// `protein_start` is derived from the genomic left edge. On a transcript that
/// runs left to right that edge is the first affected residue; on one that runs
/// right to left it is the last. The residues themselves are built from the
/// lower of the two CDS coordinates for a shrinking change and from `cds_start`
/// otherwise (see `predict_coding_consequence`), so it is exactly the shrinking
/// change on the reverse strand whose span arrives anchored at its end. An
/// insertion is anchored at `protein_start` on either strand.
///
/// This is a property of the coordinate convention, not of the sequence, which
/// is why it can be decided here rather than guessed from the peptide.
fn anchored_at_span_end(strand: Strand, reference_len: usize, alternate_len: usize) -> bool {
    strand == Strand::Reverse && alternate_len < reference_len
}

/// The positions at which `reference` may sit, given the caller's anchor, most
/// likely first.
///
/// `protein_start` is the first affected residue on the forward strand, but for
/// a shrinking change it can be the *last*: the field is derived from the
/// genomic left edge, which is the end of the affected range when the transcript
/// runs right to left (see #89). Both readings describe the same span, so a
/// reference of length n is anchored either at `protein_start` or at
/// `protein_start - (n - 1)`.
///
/// Those two and no others. The pair is one span read from either end, so any
/// further candidate would be a scan for a coincidental repeat rather than a
/// consequence of the coordinate convention. An empty or single-residue
/// reference has one candidate, which leaves insertions - where `protein_start`
/// is already the far end of the pair - on exactly the path they had before.
///
/// `end_first` puts them in the order [`anchored_at_span_end`] determines, which
/// is what makes the choice between them a reading of the convention rather than
/// a coincidence. Taking the first *corroborated* candidate is not enough on its
/// own: where the reference is periodic with period n-1 both ends corroborate,
/// and before #96 the wrong one won whenever it was tried first. `MEGEGEA` with
/// `EGE` arriving at `protein_start` 4 is the smallest case - residues 2-4 and
/// 4-6 both read `EGE`, and deleting them leaves different proteins (`MGEA` and
/// `MEGA`), so the 3'-rule cannot reconcile the two descriptions afterwards.
///
/// The other end stays as a fallback rather than being dropped. The residues are
/// not guaranteed to sit at either scalar - over 37,122 in-frame ClinVar rows
/// they sat at `protein_start` in 71.1% of cases, one residue earlier in 11.9%,
/// and at neither in 17.0% - so ordering is all the strand can honestly buy. The
/// peptide still decides.
fn anchor_candidates(
    protein_start: u64,
    reference_len: usize,
    end_first: bool,
) -> [Option<u64>; 2] {
    let from_end = reference_len
        .checked_sub(1)
        .filter(|&back| back > 0)
        .and_then(|back| protein_start.checked_sub(back as u64))
        // Residues are numbered from 1, so 0 is not a position to try.
        .filter(|&anchor| anchor > 0);
    match from_end {
        Some(other) if end_first => [Some(other), Some(protein_start)],
        _ => [Some(protein_start), from_end],
    }
}

/// Generate HGVSp notation for an in-frame indel — deletion, delins, insertion
/// or duplication.
///
/// `ref_aas` are the affected reference residues (one-letter, in order);
/// `alt_aas` is the replacement ("-" or empty for a pure deletion, longer than
/// `ref_aas` for an insertion).
///
///   ENSP0:p.Phe157del            single-residue deletion
///   ENSP0:p.Tyr43_Gln45del       multi-residue deletion
///   ENSP0:p.Asn2173_Leu2174delinsLys   delins
///   ENSP0:p.Ser92_Ser93insGly    insertion
///   ENSP0:p.Asn587_Asp600dup     insertion that repeats what it follows
///
/// This is the only correct rendering for an in-frame indel: `hgvsp` compares
/// just the first residue of each side, so it would describe these as a
/// substitution — synonymous when that residue is unchanged, missense when it
/// differs, and neither is true of a variant that changes the protein's length.
///
/// `ref_peptide` is the reference protein in one-letter codes, residue 1 at
/// index 0 — pass `Transcript::peptide`, which is translated in frame from
/// `codon_table_start_phase` and with the right codon table for the contig. Any
/// terminator it carries is trimmed here. Given one, insertions and deletions
/// are normalised per the HGVS 3'-rule - shifted as far C-terminal as they can
/// go, with duplications collapsed to `dup`. `delins` is not shifted.
///
/// This agrees with Ensembl VEP except within the last `n - 1` residues of the
/// protein, for a change of `n` residues. VEP's `_shift_3prime`
/// (`Bio::EnsEMBL::Variation::TranscriptVariationAllele`) bounds its scan at
/// `length(post_seq) - n` while comparing a single residue per step, so it
/// halts `n - 1` residues short of the terminus. We shift to the terminus,
/// which is what the 3'-rule specifies and what our HGVSc for the same variant
/// already describes. See issue #94.
///
/// `strand` says which end of `ref_aas` the caller's `protein_start` names,
/// which the coordinate convention fixes rather than leaves open: see
/// [`anchored_at_span_end`]. Pass the strand of the transcript the residues were
/// derived from - not the strand the caller finds convenient - because a wrong
/// one here reorders the anchor candidates and, where the reference is periodic,
/// selects a span the variant does not touch.
///
/// Every peptide-dependent step degrades to the unshifted description rather
/// than failing: a peptide that is absent, too short, or inconsistent with
/// `protein_start` still yields valid HGVS, just unnormalised.
pub fn hgvsp_inframe_indel(
    protein_id: &str,
    protein_start: u64,
    ref_aas: &str,
    alt_aas: &str,
    ref_peptide: Option<&[u8]>,
    strand: Strand,
) -> Option<String> {
    let strip = |s: &str| -> Vec<u8> {
        if s == "-" {
            Vec::new()
        } else {
            s.bytes().collect()
        }
    };
    let original_ref = strip(ref_aas);
    // Nothing past a terminator the change introduces is translated, so those
    // residues are not part of the protein and must not be named. Ensembl's
    // `Amino_acids` column keeps the whole translated window - `SL/MEP*S` - but
    // its HGVSp for the same row is `p.Ser269_Leu270delinsMetGluProTer`, and
    // naming the `S` after the `Ter` would describe a residue that does not
    // exist. The terminator itself is kept: it is the last thing the protein has.
    let original_alt = {
        let all = strip(alt_aas);
        match all.iter().position(|&b| b == b'*') {
            Some(terminator) => all[..=terminator].to_vec(),
            None => all,
        }
    };
    let prefix = format!("{}:p.", protein_id);
    let fallback = || unshifted_description(&prefix, protein_start, &original_ref, &original_alt);

    // Transcript::peptide runs to cdna_coding_end, so it carries the terminator
    // as a final `*` — and anything translated past an internal one on a
    // mis-annotated CDS. Those residues are not part of the protein, so bound
    // the peptide at the first terminator before shifting against it.
    let ref_peptide = ref_peptide.map(|p| match p.iter().position(|&b| b == b'*') {
        Some(terminator) => &p[..terminator],
        None => p,
    });
    // Everything below is positioned relative to an anchor, so none of it is
    // safe unless the peptide corroborates that the caller's residues really sit
    // there. They do not always: for a shrinking change like `FF/F` the call
    // sites pass the end of the affected range rather than its start, so try
    // reading the span from its other end before giving up.
    //
    // Trusting only `protein_start` was not merely a missed normalisation. The
    // un-normalised description takes its residue letters from `ref_aas` and its
    // numbers from the anchor, so a wrong anchor emits a position the protein
    // does not have - `p.Ser1092_Phe1094delinsPhe` on ENSP00000247087 where
    // residues 1092-1094 are FQP and the SSF is at 1090-1092. The letters and
    // the numbers contradict each other. Over 47,013 ClinVar indels that was
    // 17,654 descriptions.
    let Some((peptide, protein_start)) = ref_peptide.and_then(|p| {
        let end_first = anchored_at_span_end(strand, original_ref.len(), original_alt.len());
        anchor_candidates(protein_start, original_ref.len(), end_first)
            .into_iter()
            .flatten()
            .find(|&anchor| peptide_carries(p, anchor, &original_ref))
            .map(|anchor| (p, anchor))
    }) else {
        return fallback();
    };

    // Ensembl's `start_lost` (`VariationEffect.pm` release/115, l. 851): the
    // replacement begins at the initiator and preserves the reference residues
    // at neither end of itself. What the ribosome does then - start at a
    // downstream ATG, or not at all - is not something the sequence tells us, so
    // HGVS marks the consequence unknown with `?` rather than naming residues of
    // a protein that may never be made. Evaluated on the untrimmed residues,
    // because that is what the predicate reads: for `XP/XAKSTVGA` the shared
    // leading `X` trims away and the description starts at residue 2, and VEP
    // still writes `p.Pro2?`.
    //
    // The pure-deletion arm below has its own, older test for this and keeps it:
    // it decides after 3'-shifting, which can move the deletion off residue 1.
    let start_lost = protein_start == 1
        && !original_alt.starts_with(&original_ref[..])
        && !original_alt.ends_with(&original_ref[..]);

    // Reduce to the minimal changed region: residues shared at either end are
    // not part of the description. Trimming the prefix moves the start right.
    let mut reference = original_ref.clone();
    let mut alternate = original_alt.clone();
    let mut start = protein_start;
    while !reference.is_empty() && !alternate.is_empty() && reference[0] == alternate[0] {
        reference.remove(0);
        alternate.remove(0);
        start += 1;
    }
    while !reference.is_empty() && !alternate.is_empty() && reference.last() == alternate.last() {
        reference.pop();
        alternate.pop();
    }

    if reference.is_empty() && alternate.is_empty() {
        return None;
    }

    if reference.is_empty() {
        // Pure insertion, sitting between residues (start - 1) and start.
        let Some(mut at) = start.checked_sub(1).map(|v| v as usize) else {
            return fallback();
        };
        let mut inserted = alternate;
        // 3'-rule: slide right while the residue the insertion sits in front of
        // is the one it would place there.
        while at < peptide.len() && peptide[at] == inserted[0] {
            inserted.rotate_left(1);
            at += 1;
        }
        // A duplication is an insertion whose residues repeat those immediately
        // before it.
        let preceding = at
            .checked_sub(inserted.len())
            .and_then(|lo| peptide.get(lo..at));
        if preceding == Some(inserted.as_slice()) {
            let dup_start = (at - inserted.len() + 1) as u64;
            return Some(format!(
                "{}{}dup",
                prefix,
                residue_span(dup_start, &inserted)
            ));
        }
        // Otherwise name the flanking pair. At a terminus there is no pair, so
        // fall back rather than drop the annotation.
        match (
            at.checked_sub(1).and_then(|i| peptide.get(i)),
            peptide.get(at),
        ) {
            (Some(&before), Some(&after)) => Some(format!(
                "{}{}{}_{}{}ins{}",
                prefix,
                aa_one_to_three(before),
                at,
                aa_one_to_three(after),
                at + 1,
                three_letter(&inserted)
            )),
            _ => fallback(),
        }
    } else if alternate.is_empty() {
        // Pure deletion of `reference` starting at `start`.
        let mut at = start;
        let mut residues = reference;
        {
            let len = residues.len();
            // 3'-rule: slide right while the residue following the deleted block
            // repeats the first deleted residue.
            loop {
                let first = at.checked_sub(1).and_then(|i| peptide.get(i as usize));
                let following = peptide.get(at as usize + len - 1);
                match (first, following) {
                    (Some(a), Some(b)) if a == b => at += 1,
                    _ => break,
                }
            }
            match at
                .checked_sub(1)
                .and_then(|lo| peptide.get(lo as usize..lo as usize + len))
            {
                Some(block) => residues = block.to_vec(),
                None => return fallback(),
            }
        }
        if at == 1 {
            // The deletion removes the initiation codon. What the ribosome does
            // instead - start downstream, or not at all - is not something the
            // sequence tells us, so HGVS marks the consequence unknown with `?`
            // rather than describing a protein that may never be made. This is
            // the form Ensembl VEP emits.
            return Some(format!(
                "{}{}",
                prefix,
                uncertain_from_initiator(at, &residues)
            ));
        }
        Some(format!("{}{}del", prefix, residue_span(at, &residues)))
    } else if start_lost {
        Some(format!(
            "{}{}",
            prefix,
            uncertain_from_initiator(start, &reference)
        ))
    } else {
        // Replacement of one residue run by another.
        Some(format!(
            "{}{}delins{}",
            prefix,
            residue_span(start, &reference),
            three_letter(&alternate)
        ))
    }
}

/// Generate HGVSp notation for a frameshift variant.
///
/// Scans the frameshifted sequence to find the first changed amino acid and
/// the position of the new stop codon.
///
/// Format: ENSP00000001:p.Ala498ProfsTer28
///   - Ala498 = first amino acid that changes (ref)
///   - Pro = new amino acid at that position
///   - Ter28 = new stop codon 28 positions downstream
///
/// `codon_table` lets the caller select the genetic code to translate with —
/// pass the vertebrate mitochondrial table (NCBI table 2) for MT transcripts
/// so AGA/AGG/ATA/TGA are read correctly instead of with the standard code.
pub fn hgvsp_frameshift(
    protein_id: &str,
    ref_translateable: &[u8],
    alt_translateable: &[u8],
    affected_codon_start: usize, // 0-based codon index where the frameshift starts
    codon_table: &CodonTable,
) -> Option<String> {
    let prefix = format!("{}:p.", protein_id);

    // Translate both sequences from the affected codon onwards
    let ref_start = affected_codon_start * 3;
    if ref_start + 3 > ref_translateable.len() {
        return None;
    }
    if ref_start > alt_translateable.len() {
        return None;
    }

    let ref_peptide: Vec<u8> = ref_translateable[ref_start..]
        .chunks(3)
        .filter(|c| c.len() == 3)
        .map(|c| codon_table.translate(&[c[0], c[1], c[2]]))
        .collect();

    let alt_peptide: Vec<u8> = alt_translateable[ref_start..]
        .chunks(3)
        .filter(|c| c.len() == 3)
        .map(|c| codon_table.translate(&[c[0], c[1], c[2]]))
        .collect();

    // Find the first position where amino acids differ
    let mut first_changed_offset = 0;
    for i in 0..ref_peptide.len().min(alt_peptide.len()) {
        if ref_peptide[i] != alt_peptide[i] {
            first_changed_offset = i;
            break;
        }
        // If we reach a stop codon in ref before finding a change,
        // the change starts at this position
        if ref_peptide[i] == b'*' {
            first_changed_offset = i;
            break;
        }
        first_changed_offset = i + 1;
    }

    if first_changed_offset >= ref_peptide.len() && first_changed_offset >= alt_peptide.len() {
        return None;
    }

    let first_changed_pos = affected_codon_start + first_changed_offset + 1; // 1-based
    let ref_aa = if first_changed_offset < ref_peptide.len() {
        ref_peptide[first_changed_offset]
    } else {
        b'X'
    };
    let alt_aa = if first_changed_offset < alt_peptide.len() {
        alt_peptide[first_changed_offset]
    } else {
        b'X'
    };

    let ref_aa3 = aa_one_to_three(ref_aa);
    let alt_aa3 = aa_one_to_three(alt_aa);

    // Find the new stop codon position in the alt sequence.
    // If the sequence contains unresolved (X) amino acids, use Ter? to indicate uncertainty.
    let mut stop_dist = None;
    let mut hit_unresolved = false;
    let unresolved_count = alt_peptide[first_changed_offset..]
        .iter()
        .take(10)
        .filter(|&&aa| aa == b'X')
        .count();
    let mostly_unresolved = unresolved_count > 5;
    if !mostly_unresolved {
        for (offset, &aa) in alt_peptide.iter().enumerate().skip(first_changed_offset) {
            if aa == b'*' {
                stop_dist = Some(offset - first_changed_offset + 1);
                break;
            }
            if aa == b'X' {
                hit_unresolved = true;
            }
        }
    } else {
        hit_unresolved = true;
    }

    if let Some(d) = stop_dist {
        Some(format!(
            "{}{}{}{}fsTer{}",
            prefix, ref_aa3, first_changed_pos, alt_aa3, d
        ))
    } else if hit_unresolved || mostly_unresolved {
        // Sequence has unresolved regions - can't determine stop position
        Some(format!(
            "{}{}{}{}fsTer?",
            prefix, ref_aa3, first_changed_pos, alt_aa3
        ))
    } else {
        // No stop found and sequence is clean - true extension
        Some(format!(
            "{}{}{}{}fsTer?",
            prefix, ref_aa3, first_changed_pos, alt_aa3
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fastvep_genome::mitochondrial_codon_table;

    #[test]
    fn test_hgvsp_frameshift_mitochondrial_table_differs() {
        // Same ref/alt translateable sequences, only the codon table differs.
        // Codon 0 changes (Arg CGT -> Pro CCC, same under both tables), so
        // the frameshift starts there regardless of table. Codon 1 is TGA:
        // a stop under the standard table but Trp under the vertebrate
        // mitochondrial table (NCBI table 2), so the two tables must find
        // the new stop codon (Ter) at different downstream distances.
        let ref_translateable = b"CGTCGTCGTCGT"; // Arg Arg Arg Arg
        let alt_translateable = b"CCCTGAAAATAA"; // Pro TGA(*/W) Lys TAA(*)

        let standard = CodonTable::standard();
        let mitochondrial = mitochondrial_codon_table();

        let standard_result =
            hgvsp_frameshift("ENSP1", ref_translateable, alt_translateable, 0, &standard);
        let mito_result = hgvsp_frameshift(
            "ENSP1",
            ref_translateable,
            alt_translateable,
            0,
            &mitochondrial,
        );

        // Standard table: TGA is a stop, so the new terminator is 2 codons in.
        assert_eq!(standard_result, Some("ENSP1:p.Arg1ProfsTer2".to_string()));
        // Mitochondrial table: TGA reads as Trp, so translation continues
        // past it to the real stop (TAA) 4 codons in.
        assert_eq!(mito_result, Some("ENSP1:p.Arg1ProfsTer4".to_string()));
        assert_ne!(standard_result, mito_result);
    }

    #[test]
    fn test_hgvsp_frameshift_short_alt_translateable_returns_none() {
        // Regression: there's a bounds check guarding `ref_translateable`
        // (`ref_start + 3 > ref_translateable.len()`) but nothing equivalent
        // guarded `alt_translateable[ref_start..]` on the next line. If the
        // alt sequence is shorter than `ref_start`, that slice must not
        // panic ("start index out of range") -- it should return None, same
        // as the existing ref-side guard.
        let ref_translateable = b"CGTCGTCGTCGT"; // 12 bases, ref_start=3 is in-bounds
        let alt_translateable = b"CC"; // only 2 bases -- shorter than ref_start (3)

        let standard = CodonTable::standard();
        let result = hgvsp_frameshift("ENSP1", ref_translateable, alt_translateable, 1, &standard);
        assert_eq!(result, None);
    }

    #[test]
    fn test_hgvsp_missense() {
        let result = hgvsp("ENSP00000001", 41, b'R', b'K', false);
        assert_eq!(result, Some("ENSP00000001:p.Arg41Lys".to_string()));
    }

    #[test]
    fn test_hgvsp_synonymous() {
        let result = hgvsp("ENSP00000001", 41, b'R', b'R', false);
        assert_eq!(result, Some("ENSP00000001:p.Arg41=".to_string()));
    }

    #[test]
    fn test_hgvsp_inframe_deletion_single() {
        // single-residue in-frame deletion
        let r = hgvsp_inframe_indel("ENSP00000001", 157, "F", "-", None, Strand::Forward);
        assert_eq!(r, Some("ENSP00000001:p.Phe157del".to_string()));
    }

    #[test]
    fn test_hgvsp_inframe_deletion_range() {
        // multi-residue in-frame deletion (regression for the p.Tyr43??? bug)
        let r = hgvsp_inframe_indel("ENSP00000001", 43, "YXQ", "-", None, Strand::Forward);
        assert_eq!(r, Some("ENSP00000001:p.Tyr43_Gln45del".to_string()));
    }

    #[test]
    fn test_hgvsp_inframe_delins() {
        // in-frame deletion-insertion
        let r = hgvsp_inframe_indel("ENSP00000001", 2173, "NL", "K", None, Strand::Forward);
        assert_eq!(
            r,
            Some("ENSP00000001:p.Asn2173_Leu2174delinsLys".to_string())
        );
    }

    // In-frame insertions previously fell through to the substitution branch,
    // which compares only the first residue of each side. Every case below is
    // real fastVEP output from a clinical panel, checked against Ensembl VEP.

    /// Build a reference peptide with `residues` placed at 1-based `at`, padded
    /// with a filler residue that cannot be confused with the payload.
    fn peptide_with(at: u64, residues: &str, length: usize) -> Vec<u8> {
        let mut pep = vec![b'M'; length];
        for (i, b) in residues.bytes().enumerate() {
            pep[at as usize - 1 + i] = b;
        }
        pep
    }

    #[test]
    fn test_hgvsp_inframe_insertion_collapses_to_duplication() {
        // C8A c.553_554insGGA, amino_acids "W/WR" at 185 — previously p.Trp185=.
        // Residue 186 is already Arg, so inserting Arg duplicates it.
        let pep = peptide_with(185, "WRQ", 200);
        let r = hgvsp_inframe_indel("ENSP00000001", 185, "W", "WR", Some(&pep), Strand::Forward);
        assert_eq!(r, Some("ENSP00000001:p.Arg186dup".to_string()));
    }

    #[test]
    fn test_hgvsp_inframe_insertion_multi_residue_duplication() {
        // FLT3 c.1759_1800dup, a 14-codon ITD — previously p.Asp600=. The delins
        // spelling runs to ~55 characters; the duplication form is 18.
        let dup = "NEYFYVDFREYEYD";
        let pep = peptide_with(587, &format!("{}K", dup), 700);
        let alt = format!("D{}", dup);
        let r = hgvsp_inframe_indel("ENSP00000001", 600, "D", &alt, Some(&pep), Strand::Forward);
        assert_eq!(r, Some("ENSP00000001:p.Asn587_Asp600dup".to_string()));
    }

    #[test]
    fn test_hgvsp_inframe_insertion_true_insertion_uses_ins_form() {
        // ITPKB c.275_276insGGT, amino_acids "S/SG" at 92 — previously p.Ser92=.
        // Gly does not repeat the preceding residues, so it stays an insertion.
        let pep = peptide_with(92, "SSK", 200);
        let r = hgvsp_inframe_indel("ENSP00000001", 92, "S", "SG", Some(&pep), Strand::Forward);
        assert_eq!(r, Some("ENSP00000001:p.Ser92_Ser93insGly".to_string()));
    }

    #[test]
    fn test_hgvsp_inframe_deletion_is_three_prime_shifted() {
        // A deletion inside a homopolymer run is reported at the most C-terminal
        // position it can occupy: deleting one Ala from AAA at 2..4 is p.Ala4del.
        let pep: Vec<u8> = "MAAAGK".bytes().collect();
        let r = hgvsp_inframe_indel("ENSP00000001", 2, "A", "-", Some(&pep), Strand::Forward);
        assert_eq!(r, Some("ENSP00000001:p.Ala4del".to_string()));
    }

    #[test]
    fn test_hgvsp_inframe_deletion_shifts_onto_the_final_residue() {
        // The 3'-rule runs to the end of the protein, not to one residue short
        // of it. Ensembl VEP stops `n - 1` residues early for an n-residue
        // change, so this pins our answer against a future attempt to reproduce
        // that off-by-one in the name of compatibility. Real case: NT5C2
        // c.1674_1679del is p.Glu560_Glu561del on a 561-residue protein, which
        // VEP reports as p.Glu559_Glu560del.
        let pep: Vec<u8> = "MKEEEEE*".bytes().collect();
        let r = hgvsp_inframe_indel("ENSP00000001", 3, "EE", "-", Some(&pep), Strand::Forward);
        assert_eq!(r, Some("ENSP00000001:p.Glu6_Glu7del".to_string()));
    }

    #[test]
    fn test_hgvsp_inframe_indel_without_peptide_stays_valid() {
        // No sequence context: emit an unshifted but well-formed description
        // rather than nothing, and never a substitution shape.
        let r = hgvsp_inframe_indel("ENSP00000001", 185, "W", "WR", None, Strand::Forward);
        assert_eq!(r, Some("ENSP00000001:p.Trp185delinsTrpArg".to_string()));
        let d = hgvsp_inframe_indel("ENSP00000001", 157, "F", "-", None, Strand::Forward);
        assert_eq!(d, Some("ENSP00000001:p.Phe157del".to_string()));
    }

    // The peptide is caller-derived, so every index into it has to survive a
    // transcript whose sequence disagrees with its own coordinates. Each case
    // below aborted the process before the bounds checks were added.

    /// `(case name, protein_start, ref_aas, alt_aas, peptide, expected)`.
    type UnusablePeptideCase<'a> = (&'a str, u64, &'a str, &'a str, &'a [u8], Option<&'a str>);

    #[test]
    fn test_hgvsp_inframe_indel_survives_unusable_peptides() {
        let short: Vec<u8> = "MAAAGK".bytes().collect();
        let cases: Vec<UnusablePeptideCase<'_>> = vec![
            // Deleted block overruns the peptide end.
            (
                "block overruns end",
                6,
                "KX",
                "-",
                &short,
                Some("ENSP00000001:p.Lys6_Xaa7del"),
            ),
            // protein_start past the peptide entirely.
            (
                "start beyond peptide",
                100,
                "AK",
                "-",
                &short,
                Some("ENSP00000001:p.Ala100_Lys101del"),
            ),
            // Insertion anchored past the peptide.
            ("insertion beyond peptide", 500, "-", "R", &short, None),
            // Empty peptide (truncated transcript).
            (
                "empty peptide",
                1,
                "F",
                "-",
                &[],
                Some("ENSP00000001:p.Phe1del"),
            ),
            // Insertion one residue past the C-terminus.
            (
                "insertion at C-terminus",
                7,
                "K",
                "KG",
                &short,
                Some("ENSP00000001:p.Lys7delinsLysGly"),
            ),
            // protein_start of 0 would underflow a 1-based conversion.
            (
                "zero protein_start",
                0,
                "F",
                "-",
                &short,
                Some("ENSP00000001:p.Phe0del"),
            ),
        ];
        for (name, start, reference, alternate, pep, expected) in cases {
            let got = hgvsp_inframe_indel(
                "ENSP00000001",
                start,
                reference,
                alternate,
                Some(pep),
                Strand::Forward,
            );
            assert_eq!(got.as_deref(), expected, "case: {name}");
        }
    }

    #[test]
    fn test_hgvsp_inframe_indel_does_not_shift_onto_the_terminator() {
        // Transcript::peptide ends with `*`. A change abutting the stop must not
        // shift onto it or name it as a flanking residue: Ter is not a residue
        // of the protein, and a position at or past it is not a real position.
        let pep: Vec<u8> = "MKKG*".bytes().collect();

        // Deleting one Lys from the KK run shifts to the 3'-most Lys (3), not
        // onto Gly4 or the terminator.
        let deletion =
            hgvsp_inframe_indel("ENSP00000001", 2, "K", "-", Some(&pep), Strand::Forward);
        assert_eq!(deletion, Some("ENSP00000001:p.Lys3del".to_string()));

        // An insertion immediately before the terminator has no residue on its
        // 3' side once the stop is excluded, so it falls back rather than
        // emitting a Ter-flanked range.
        let insertion =
            hgvsp_inframe_indel("ENSP00000001", 4, "G", "GS", Some(&pep), Strand::Forward);
        assert_eq!(
            insertion,
            Some("ENSP00000001:p.Gly4delinsGlySer".to_string())
        );
        for out in [deletion, insertion] {
            let out = out.unwrap();
            assert!(!out.contains("Ter"), "terminator named as a residue: {out}");
        }
    }

    #[test]
    fn test_hgvsp_inframe_indel_ignores_a_peptide_that_disagrees() {
        // The peptide says Gln at 2; the caller says Phe. Trusting the peptide
        // would emit a confident, well-formed, wrong description (p.Gln6del).
        // Fall back to the caller's residues instead.
        let pep: Vec<u8> = "MQQQQQK".bytes().collect();
        let r = hgvsp_inframe_indel("ENSP00000001", 2, "F", "-", Some(&pep), Strand::Forward);
        assert_eq!(r, Some("ENSP00000001:p.Phe2del".to_string()));
    }

    #[test]
    fn test_hgvsp_inframe_indel_reads_the_span_from_either_end() {
        // For a shrinking change the call sites do not always pass the start of
        // `ref_aas`: real ClinVar output has `Protein_position` 328-329 with
        // `Amino_acids` FF/F arriving here as protein_start 329, the end of the
        // pair rather than its start. Both numbers name the same span, so the
        // description must not depend on which end the caller happened to send.
        //
        // The peptide is M R I F F A S M, so the F pair is at 4-5. Which end
        // the caller sends follows from the strand (#96), so the two readings
        // are paired with the strand that produces them - but the answer has to
        // be the same either way, which is what this test is for.
        let pep: Vec<u8> = "MRIFFASM".bytes().collect();
        for (anchor, strand) in [(4u64, Strand::Forward), (5, Strand::Reverse)] {
            assert_eq!(
                hgvsp_inframe_indel("ENSP00000001", anchor, "FF", "F", Some(&pep), strand),
                Some("ENSP00000001:p.Phe5del".to_string()),
                "anchor {anchor} on {strand:?} should describe the same deletion"
            );
        }

        // And the pairing is a preference, not a requirement: a reverse-strand
        // caller whose residues really do sit at `protein_start` still
        // normalises, because the other end stays as a fallback.
        assert_eq!(
            hgvsp_inframe_indel("ENSP00000001", 4, "FF", "F", Some(&pep), Strand::Reverse),
            Some("ENSP00000001:p.Phe5del".to_string()),
            "the unpreferred end must still be tried"
        );

        // This case used to take the un-normalised path and emit
        // `p.Phe5_Phe6delinsPhe`, naming residue 6 as Phe when the peptide has
        // Ala there. Reading the span from its other end is what fixes it: the
        // letters and the numbers now agree with the protein.
    }

    #[test]
    fn a_deletion_of_the_initiation_codon_is_described_as_unresolvable() {
        // Removing the start codon leaves the protein's fate undetermined: the
        // ribosome may initiate downstream, or not at all, and the sequence does
        // not say which. Describing it as an ordinary deletion asserts a protein
        // that begins where nothing says it begins.
        //
        // Both shapes come from real ClinVar rows checked against Ensembl VEP:
        // KCNA2 `MNII/I` at residues 1-4, and POLE `EA/A` at 1-2.
        let kcna2: Vec<u8> = "MNIIDIVAIIPY".bytes().collect();
        assert_eq!(
            hgvsp_inframe_indel(
                "ENSP00000491354",
                1,
                "MNII",
                "I",
                Some(&kcna2),
                Strand::Forward
            ),
            Some("ENSP00000491354:p.MetAsnIle1_?3".to_string())
        );

        let pole: Vec<u8> = "EAKRQ".bytes().collect();
        assert_eq!(
            hgvsp_inframe_indel(
                "ENSP00000500921",
                1,
                "EA",
                "A",
                Some(&pole),
                Strand::Forward
            ),
            Some("ENSP00000500921:p.Glu1?".to_string())
        );
    }

    #[test]
    fn a_change_that_loses_the_initiator_is_unresolvable() {
        // The marker is specific to the initiation codon. A deletion anywhere
        // else is an ordinary `del`, including one that starts at residue 2.
        let pep: Vec<u8> = "MKFFASM".bytes().collect();
        assert_eq!(
            hgvsp_inframe_indel("ENSP00000001", 2, "K", "-", Some(&pep), Strand::Forward),
            Some("ENSP00000001:p.Lys2del".to_string())
        );

        // A delins that takes out the initiator is *not* an ordinary delins.
        // Replacing `MK` with `W` removes the ATG, so where translation begins
        // is exactly what is no longer known, and Ensembl marks it `?` for every
        // shape of change that loses the start - `p.Met1?` for `M/T`,
        // `p.MetAla1_?2` for `MA/IS`, `p.Pro2?` for `XP/XAKSTVGA`, all from real
        // VEP 115.1. Naming the replacement instead described a protein that may
        // never be made.
        assert_eq!(
            hgvsp_inframe_indel("ENSP00000001", 1, "MK", "W", Some(&pep), Strand::Forward),
            Some("ENSP00000001:p.MetLys1_?2".to_string())
        );

        // A replacement that keeps the reference residues at one end of itself
        // has not lost the start. `MK` -> `MWK` is an insertion between them,
        // and it is described as one.
        assert_eq!(
            hgvsp_inframe_indel("ENSP00000001", 1, "MK", "MWK", Some(&pep), Strand::Forward),
            Some("ENSP00000001:p.Met1_Lys2insTrp".to_string())
        );
    }

    #[test]
    fn test_hgvsp_inframe_indel_still_declines_when_neither_end_corroborates() {
        // Reading from the other end is a consequence of the coordinate
        // convention, not a search. When the peptide carries the reference at
        // neither end of the span, there is no evidence for any anchor and the
        // un-normalised description is still the honest answer.
        let pep: Vec<u8> = "MRIFFASM".bytes().collect();
        for strand in [Strand::Forward, Strand::Reverse] {
            let r = hgvsp_inframe_indel("ENSP00000001", 5, "KK", "K", Some(&pep), strand);
            assert_eq!(
                r,
                Some("ENSP00000001:p.Lys5_Lys6delinsLys".to_string()),
                "on {strand:?} an uncorroborated span must stay unshifted"
            );
        }

        // A single residue has only one candidate: there is no other end to try,
        // so an uncorroborated anchor cannot be rescued and must not be guessed.
        let r = hgvsp_inframe_indel("ENSP00000001", 2, "W", "-", Some(&pep), Strand::Forward);
        assert_eq!(r, Some("ENSP00000001:p.Trp2del".to_string()));
    }

    #[test]
    fn anchor_candidates_offers_the_other_end_only_when_there_is_one() {
        // An empty or single reference has one reading, which keeps insertions -
        // where protein_start is already the far end of the pair - unchanged.
        assert_eq!(anchor_candidates(10, 0, false), [Some(10), None]);
        assert_eq!(anchor_candidates(10, 1, false), [Some(10), None]);
        assert_eq!(anchor_candidates(10, 2, false), [Some(10), Some(9)]);
        assert_eq!(anchor_candidates(10, 5, false), [Some(10), Some(6)]);
        // Never underflows past residue 1.
        assert_eq!(anchor_candidates(2, 5, false), [Some(2), None]);
        assert_eq!(anchor_candidates(1, 2, false), [Some(1), None]);
    }

    #[test]
    fn anchor_candidates_puts_the_determined_end_first() {
        // Same pair, opposite order: which one is tried first is what decides a
        // periodic reference, where both corroborate (#96).
        assert_eq!(anchor_candidates(10, 5, true), [Some(6), Some(10)]);
        assert_eq!(anchor_candidates(10, 5, false), [Some(10), Some(6)]);

        // With only one candidate there is nothing to order, and the request for
        // the other end must not invent one or lose the one there is.
        assert_eq!(anchor_candidates(10, 1, true), [Some(10), None]);
        assert_eq!(anchor_candidates(10, 0, true), [Some(10), None]);
        assert_eq!(anchor_candidates(2, 5, true), [Some(2), None]);
        assert_eq!(anchor_candidates(1, 2, true), [Some(1), None]);
    }

    #[test]
    fn only_a_shrinking_change_on_the_reverse_strand_is_anchored_at_its_end() {
        // The four combinations, because the rule is a reading of how the
        // coordinates were built (see `predict_coding_consequence`): the
        // residues come from the lower CDS coordinate for a shrinking change and
        // from `cds_start` otherwise, and only on the reverse strand are those
        // two different ends of the span.
        assert!(anchored_at_span_end(Strand::Reverse, 3, 0), "deletion");
        assert!(
            anchored_at_span_end(Strand::Reverse, 3, 1),
            "shrinking delins"
        );
        assert!(
            !anchored_at_span_end(Strand::Forward, 3, 0),
            "forward strand"
        );
        assert!(!anchored_at_span_end(Strand::Reverse, 1, 3), "insertion");
        assert!(
            !anchored_at_span_end(Strand::Reverse, 3, 3),
            "equal-length replacement does not shrink"
        );
    }

    #[test]
    fn a_periodic_reference_on_the_reverse_strand_names_the_span_the_caller_meant() {
        // Issue #96. `EGE` sits at residues 2-4 *and* at 4-6, so both anchors
        // are corroborated and taking whichever came first picked the wrong one.
        // These are not two spellings of one variant: deleting 2-4 leaves MGEA
        // and deleting 4-6 leaves MEGA, so the 3'-rule cannot merge them
        // afterwards - the wrong answer named residues the variant never touched.
        let pep: Vec<u8> = "MEGEGEA".bytes().collect();
        assert_eq!(
            hgvsp_inframe_indel("ENSP00000001", 4, "EGE", "-", Some(&pep), Strand::Reverse),
            Some("ENSP00000001:p.Glu2_Glu4del".to_string())
        );

        // Period 3 at n = 4, the same shape one residue longer.
        let pep: Vec<u8> = "MABCABCA".bytes().collect();
        assert_eq!(
            hgvsp_inframe_indel("ENSP00000001", 5, "ABCA", "-", Some(&pep), Strand::Reverse),
            Some("ENSP00000001:p.Ala2_Ala5del".to_string())
        );
    }

    #[test]
    fn a_periodic_reference_on_the_forward_strand_is_read_from_the_start() {
        // The mirror of the case above, and the reason the fix is an ordering
        // rather than a preference for the earlier residue: on the forward strand
        // `protein_start` *is* the first affected residue, so the same peptide
        // and the same reference must resolve to the other span.
        let pep: Vec<u8> = "MEGEGEA".bytes().collect();
        assert_eq!(
            hgvsp_inframe_indel("ENSP00000001", 4, "EGE", "-", Some(&pep), Strand::Forward),
            Some("ENSP00000001:p.Glu4_Glu6del".to_string())
        );
        // And the caller that means residues 2-4 on the forward strand says so.
        assert_eq!(
            hgvsp_inframe_indel("ENSP00000001", 2, "EGE", "-", Some(&pep), Strand::Forward),
            Some("ENSP00000001:p.Glu2_Glu4del".to_string())
        );
    }

    /// Deterministic pseudo-random source. A seeded LCG rather than a `rand`
    /// dependency: the sweep below has to fail the same way twice, or a failure
    /// cannot be investigated.
    struct Lcg(u64);

    impl Lcg {
        fn below(&mut self, n: usize) -> usize {
            self.0 = self
                .0
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            ((self.0 >> 33) as usize) % n
        }
    }

    /// Apply an emitted HGVSp description back to the reference peptide, and
    /// check on the way that every residue it names is the residue the peptide
    /// has at that position.
    ///
    /// This is the inverse of the description, which is what makes it a check
    /// worth having: a span that is well-formed, corroborated and *wrong* still
    /// reconstructs the wrong protein.
    fn apply_to_peptide(description: &str, peptide: &str) -> String {
        let one = |aa3: &str| -> char {
            match aa3 {
                "Ala" => 'A',
                "Glu" => 'E',
                "Gly" => 'G',
                "Lys" => 'K',
                "Met" => 'M',
                other => panic!("unexpected residue {other} in {description}"),
            }
        };
        let body = description
            .split(":p.")
            .nth(1)
            .unwrap_or_else(|| panic!("no :p. in {description}"));
        let (span, replacement) = match body.split_once("delins") {
            Some((span, inserted)) => (
                span,
                inserted
                    .as_bytes()
                    .chunks(3)
                    .map(|c| one(std::str::from_utf8(c).unwrap()))
                    .collect::<String>(),
            ),
            None => (
                body.strip_suffix("del")
                    .unwrap_or_else(|| panic!("neither del nor delins: {description}")),
                String::new(),
            ),
        };

        // `Glu2` or `Glu2_Glu4`, and both ends are checked against the peptide.
        let mut bounds = Vec::new();
        for part in span.split('_') {
            let (aa3, digits) = part.split_at(3);
            let pos: usize = digits
                .parse()
                .unwrap_or_else(|_| panic!("unparsable position in {description}"));
            assert_eq!(
                peptide.chars().nth(pos - 1),
                Some(one(aa3)),
                "{description} names {aa3} at {pos}, peptide has {:?} ({peptide})",
                peptide.chars().nth(pos - 1)
            );
            bounds.push(pos);
        }
        let lo = bounds[0];
        let hi = *bounds.last().unwrap();
        assert!(lo <= hi, "inverted span in {description}");
        format!("{}{}{}", &peptide[..lo - 1], replacement, &peptide[hi..])
    }

    #[test]
    fn every_description_reconstructs_the_protein_the_variant_produces() {
        // The property that matters, over 4,000 shrinking in-frame changes on a
        // four-residue alphabet chosen to make periodic references common: apply
        // the description back to the reference peptide and you must get the
        // protein the variant actually produces. A description can be
        // well-formed, name real residues, and still fail this - that is exactly
        // what #96 was, and what #91 was before it.
        //
        // Both strands, each with the anchor its coordinate convention produces:
        // the start of the span on the forward strand, the end of it on the
        // reverse. Spans start at residue 2 or later, because a deletion reaching
        // the initiation codon is deliberately described as unresolvable (`?`)
        // rather than as a reconstructible event.
        const ALPHABET: &[u8] = b"AEGK";
        let mut rng = Lcg(0x5EED_1234_9ABC_DEF0);
        let mut checked = 0;

        for _ in 0..4_000 {
            let len = 6 + rng.below(14);
            let peptide: String = (0..len)
                .map(|_| ALPHABET[rng.below(ALPHABET.len())] as char)
                .collect();
            let n = 1 + rng.below(4);
            if len < n + 2 {
                continue;
            }
            // 1-based, and never residue 1.
            let lo = 2 + rng.below(len - n);
            let hi = lo + n - 1;
            let reference = &peptide[lo - 1..hi];
            // Shrinking: a pure deletion, or a replacement by fewer residues.
            let m = rng.below(n);
            let replacement: String = (0..m)
                .map(|_| ALPHABET[rng.below(ALPHABET.len())] as char)
                .collect();
            let alt_aas = if m == 0 { "-".to_string() } else { replacement };
            let expected = format!(
                "{}{}{}",
                &peptide[..lo - 1],
                alt_aas.replace('-', ""),
                &peptide[hi..]
            );

            for (anchor, strand) in [(lo as u64, Strand::Forward), (hi as u64, Strand::Reverse)] {
                let got = hgvsp_inframe_indel(
                    "P",
                    anchor,
                    reference,
                    &alt_aas,
                    Some(peptide.as_bytes()),
                    strand,
                )
                .unwrap_or_else(|| {
                    panic!("no description for {reference}/{alt_aas} at {lo}-{hi} in {peptide}")
                });
                assert_eq!(
                    apply_to_peptide(&got, &peptide),
                    expected,
                    "{got} does not reconstruct {expected} from {peptide} \
                     ({reference}/{alt_aas} at {lo}-{hi}, anchor {anchor} on {strand:?})"
                );
                checked += 1;
            }
        }
        assert!(checked > 5_000, "sweep covered only {checked} cases");
    }

    #[test]
    fn a_two_residue_homopolymer_reads_the_same_from_either_end() {
        // Where the reference is a homopolymer inside a longer run both anchors
        // are corroborated at n = 2, but the 3'-shift converges on one answer, so
        // the strand cannot change it. Worth pinning: it bounds the exposure the
        // ordering fix was needed for to n >= 3.
        let pep: Vec<u8> = "MAKKKA".bytes().collect();
        for (anchor, strand) in [
            (3u64, Strand::Forward),
            (4, Strand::Forward),
            (3, Strand::Reverse),
            (4, Strand::Reverse),
        ] {
            assert_eq!(
                hgvsp_inframe_indel("ENSP00000001", anchor, "KK", "-", Some(&pep), strand),
                Some("ENSP00000001:p.Lys4_Lys5del".to_string()),
                "anchor {anchor} on {strand:?}"
            );
        }
    }

    #[test]
    fn test_hgvsp_inframe_indel_never_drops_a_terminal_insertion() {
        // An insertion with no flanking pair still has to produce something —
        // returning None would be a regression against the pre-normalisation
        // behaviour, which always emitted a (wrong) substitution.
        let pep: Vec<u8> = "MKKRSTV".bytes().collect();
        for start in [1u64, 8] {
            // Both strands: an insertion is anchored at `protein_start` whichever
            // way the transcript runs, so the strand must make no difference here.
            for strand in [Strand::Forward, Strand::Reverse] {
                let got = hgvsp_inframe_indel("ENSP00000001", start, "G", "GG", Some(&pep), strand);
                assert!(
                    got.is_some(),
                    "dropped annotation at protein_start {start} on {strand:?}"
                );
            }
        }
    }

    #[test]
    fn test_hgvsp_inframe_indel_never_emits_substitution_shape() {
        // Every in-frame insertion observed in a clinical panel run, as
        // (protein_start, amino_acids ref, amino_acids alt, prior output).
        // All fifteen previously rendered in a substitution shape: eight as
        // synonymous, seven as a plausible missense. Two (FLT3 at 598 and 596)
        // are ITDs, where a missense reading is clinically misleading.
        let insertions: &[(u64, &str, &str, &str)] = &[
            (41, "G", "AG", "p.Gly41Ala"),
            (185, "W", "WR", "p.Trp185="),
            (92, "S", "SG", "p.Ser92="),
            (2927, "E", "DE", "p.Glu2927Asp"),
            (375, "Q", "PLGPAKPPAQQ", "p.Gln375Pro"),
            (1829, "G", "GSSG", "p.Gly1829="),
            (510, "P", "QP", "p.Pro510Gln"),
            (510, "P", "QQP", "p.Pro510Gln"),
            (498, "Q", "QQ", "p.Gln498="),
            (600, "D", "DNEYFYVDFREYEYD", "p.Asp600="),
            (598, "E", "DVDFREYE", "p.Glu598Asp"),
            (596, "E", "VPSDNEYFYVDFRE", "p.Glu596Val"),
            (439, "I", "IKKK", "p.Ile439="),
            (
                1688,
                "S",
                "CSKDLEAFNPESKELLDLVEFTNEIQTLLGSS",
                "p.Ser1688Cys",
            ),
            (188, "S", "SD", "p.Ser188="),
        ];
        let deletions: &[(u64, &str, &str)] = &[(157, "F", "-"), (43, "YXQ", "-")];

        // Run each case twice: with no peptide, and with a peptide that really
        // carries the stated reference residues at the stated position.
        for &(start, ref_aas, alt_aas, prior) in insertions {
            let pep = peptide_with(start, ref_aas, start as usize + ref_aas.len() + 64);
            for context in [None, Some(pep.as_slice())] {
                // Both strands. An insertion is anchored at `protein_start`
                // regardless (#96), so every assertion below has to hold either
                // way; a strand-dependent answer here would be a bug.
                for strand in [Strand::Forward, Strand::Reverse] {
                    let out = hgvsp_inframe_indel(
                        "ENSP00000001",
                        start,
                        ref_aas,
                        alt_aas,
                        context,
                        strand,
                    )
                    .expect("in-frame indel must produce a protein description");
                    let change = out.split(":p.").nth(1).unwrap();
                    // "delins" contains both "del" and "ins", so require one of the
                    // whole forms rather than a substring of another.
                    assert!(
                        change.ends_with("del")
                            || change.ends_with("dup")
                            || change.contains("delins")
                            || change.contains("ins"),
                        "{ref_aas}/{alt_aas} at {start} (was {prior}) rendered \
                     without an indel form: {out}"
                    );
                    assert!(!out.contains('?'), "placeholder residue in {out}");
                    assert!(!out.ends_with('='), "synonymous shape: {out} (was {prior})");
                    assert_ne!(
                        out.split(':').nth(1).unwrap(),
                        prior,
                        "still emitting {prior}"
                    );
                }
            }
        }

        for &(start, ref_aas, alt_aas) in deletions {
            let pep = peptide_with(start, ref_aas, start as usize + ref_aas.len() + 64);
            for context in [None, Some(pep.as_slice())] {
                for strand in [Strand::Forward, Strand::Reverse] {
                    let out = hgvsp_inframe_indel(
                        "ENSP00000001",
                        start,
                        ref_aas,
                        alt_aas,
                        context,
                        strand,
                    )
                    .expect("deletion must produce a protein description");
                    assert!(
                        out.ends_with("del"),
                        "deletion regressed on {strand:?}: {out}"
                    );
                }
            }
        }
    }

    #[test]
    fn test_hgvsp_stop_gained() {
        let result = hgvsp("ENSP00000001", 41, b'R', b'*', false);
        assert_eq!(result, Some("ENSP00000001:p.Arg41Ter".to_string()));
    }

    #[test]
    fn test_hgvsp_frameshift() {
        let result = hgvsp("ENSP00000001", 41, b'R', b'X', true);
        assert_eq!(result, Some("ENSP00000001:p.Arg41fs".to_string()));
    }

    #[test]
    fn test_hgvsp_stop_lost() {
        let result = hgvsp("ENSP00000001", 100, b'*', b'R', false);
        assert_eq!(result, Some("ENSP00000001:p.Arg100ext*?".to_string()));
    }
}
