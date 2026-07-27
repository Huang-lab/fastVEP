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
        return Some(format!("{}{}{}{}",prefix, ref_aa3, protein_pos, alt_aa3));
    }

    if ref_aa == b'*' {
        // Stop lost - extension
        return Some(format!("{}{}{}ext*?",prefix, alt_aa3, protein_pos));
    }

    // Missense
    Some(format!("{}{}{}{}", prefix, ref_aa3, protein_pos, alt_aa3))
}

/// Generate HGVSp notation for an in-frame deletion or delins.
///
/// `ref_aas` are the affected reference residues (one-letter, in order) starting
/// at `protein_start`; `alt_aas` is the replacement ("-" or empty for a pure
/// deletion). Produces standard HGVS:
///   ENSP0:p.Phe157del                 (single-residue deletion)
///   ENSP0:p.Tyr43_Gln45del            (multi-residue deletion)
///   ENSP0:p.Asn2173_Leu2174delinsLys  (in-frame delins)
///
/// This replaces the incorrect missense-style output (e.g. `p.Tyr43???`, where
/// the deletion marker `-` had no three-letter code) for in-frame deletions.
pub fn hgvsp_inframe_deletion(
    protein_id: &str,
    protein_start: u64,
    ref_aas: &str,
    alt_aas: &str,
) -> Option<String> {
    let ref_bytes = ref_aas.as_bytes();
    if ref_bytes.is_empty() {
        return None;
    }
    let prefix = format!("{}:p.", protein_id);
    let first = aa_one_to_three(ref_bytes[0]);
    let range = if ref_bytes.len() == 1 {
        format!("{}{}", first, protein_start)
    } else {
        let last = aa_one_to_three(ref_bytes[ref_bytes.len() - 1]);
        let end = protein_start + ref_bytes.len() as u64 - 1;
        format!("{}{}_{}{}", first, protein_start, last, end)
    };
    if alt_aas.is_empty() || alt_aas == "-" {
        Some(format!("{}{}del", prefix, range))
    } else {
        let ins: String = alt_aas.bytes().map(aa_one_to_three).collect();
        Some(format!("{}{}delins{}", prefix, range, ins))
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
    let unresolved_count = alt_peptide[first_changed_offset..].iter()
        .take(10)
        .filter(|&&aa| aa == b'X')
        .count();
    let mostly_unresolved = unresolved_count > 5;
    if !mostly_unresolved {
        for i in first_changed_offset..alt_peptide.len() {
            if alt_peptide[i] == b'*' {
                stop_dist = Some(i - first_changed_offset + 1);
                break;
            }
            if alt_peptide[i] == b'X' {
                hit_unresolved = true;
            }
        }
    } else {
        hit_unresolved = true;
    }

    if let Some(d) = stop_dist {
        Some(format!("{}{}{}{}fsTer{}", prefix, ref_aa3, first_changed_pos, alt_aa3, d))
    } else if hit_unresolved || mostly_unresolved {
        // Sequence has unresolved regions - can't determine stop position
        Some(format!("{}{}{}{}fsTer?", prefix, ref_aa3, first_changed_pos, alt_aa3))
    } else {
        // No stop found and sequence is clean - true extension
        Some(format!("{}{}{}{}fsTer?", prefix, ref_aa3, first_changed_pos, alt_aa3))
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
        let mito_result =
            hgvsp_frameshift("ENSP1", ref_translateable, alt_translateable, 0, &mitochondrial);

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
        let result =
            hgvsp_frameshift("ENSP1", ref_translateable, alt_translateable, 1, &standard);
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
        let r = hgvsp_inframe_deletion("ENSP00000001", 157, "F", "-");
        assert_eq!(r, Some("ENSP00000001:p.Phe157del".to_string()));
    }

    #[test]
    fn test_hgvsp_inframe_deletion_range() {
        // multi-residue in-frame deletion (regression for the p.Tyr43??? bug)
        let r = hgvsp_inframe_deletion("ENSP00000001", 43, "YXQ", "-");
        assert_eq!(r, Some("ENSP00000001:p.Tyr43_Gln45del".to_string()));
    }

    #[test]
    fn test_hgvsp_inframe_delins() {
        // in-frame deletion-insertion
        let r = hgvsp_inframe_deletion("ENSP00000001", 2173, "NL", "K");
        assert_eq!(r, Some("ENSP00000001:p.Asn2173_Leu2174delinsLys".to_string()));
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
