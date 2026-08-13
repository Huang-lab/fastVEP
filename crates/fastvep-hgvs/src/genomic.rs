use fastvep_core::Allele;

/// Generate HGVSg (genomic) notation.
///
/// Format: chromosome:g.positionREF>ALT (for SNVs)
///         chromosome:g.start_enddelREF (for deletions)
///         chromosome:g.pos_pos+1insALT (for insertions)
pub fn hgvsg(
    chrom: &str,
    start: u64,
    end: u64,
    ref_allele: &Allele,
    alt_allele: &Allele,
) -> String {
    let prefix = format!("{}:g.", chrom);

    match (ref_allele, alt_allele) {
        // SNV
        (Allele::Sequence(ref_bases), Allele::Sequence(alt_bases))
            if ref_bases.len() == 1 && alt_bases.len() == 1 =>
        {
            format!(
                "{}{}{}>{}",
                prefix, start, ref_bases[0] as char, alt_bases[0] as char
            )
        }
        // Deletion
        (Allele::Sequence(ref_bases), Allele::Deletion) => {
            if ref_bases.len() == 1 {
                format!("{}{}del", prefix, start)
            } else {
                format!("{}{}_{}del", prefix, start, end)
            }
        }
        // Insertion (ref is deletion marker, meaning zero-length ref interval)
        (Allele::Deletion, Allele::Sequence(alt_bases)) => {
            // In Ensembl coords, insertion is between end and end+1 (since start > end)
            format!(
                "{}{}_{}ins{}",
                prefix,
                end,
                end + 1,
                std::str::from_utf8(alt_bases).unwrap_or("?")
            )
        }
        // MNV (substitution of multiple bases)
        (Allele::Sequence(ref_bases), Allele::Sequence(alt_bases))
            if ref_bases.len() == alt_bases.len() && ref_bases.len() > 1 =>
        {
            format!(
                "{}{}_{}delins{}",
                prefix,
                start,
                end,
                std::str::from_utf8(alt_bases).unwrap_or("?")
            )
        }
        // Complex indel (different lengths)
        (Allele::Sequence(_), Allele::Sequence(alt_bases)) => {
            if start == end {
                format!(
                    "{}{}delins{}",
                    prefix,
                    start,
                    std::str::from_utf8(alt_bases).unwrap_or("?")
                )
            } else {
                format!(
                    "{}{}_{}delins{}",
                    prefix,
                    start,
                    end,
                    std::str::from_utf8(alt_bases).unwrap_or("?")
                )
            }
        }
        _ => format!("{}{}?", prefix, start),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hgvsg_snv() {
        let result = hgvsg(
            "chr1",
            12345,
            12345,
            &Allele::Sequence(b"A".to_vec()),
            &Allele::Sequence(b"G".to_vec()),
        );
        assert_eq!(result, "chr1:g.12345A>G");
    }

    #[test]
    fn test_hgvsg_deletion_single() {
        let result = hgvsg(
            "chr1",
            100,
            100,
            &Allele::Sequence(b"A".to_vec()),
            &Allele::Deletion,
        );
        assert_eq!(result, "chr1:g.100del");
    }

    #[test]
    fn test_hgvsg_deletion_multi() {
        let result = hgvsg(
            "chr1",
            100,
            102,
            &Allele::Sequence(b"ACG".to_vec()),
            &Allele::Deletion,
        );
        assert_eq!(result, "chr1:g.100_102del");
    }

    #[test]
    fn test_hgvsg_insertion() {
        let result = hgvsg(
            "chr1",
            101,
            100, // insertion: start > end in Ensembl coords
            &Allele::Deletion,
            &Allele::Sequence(b"TCG".to_vec()),
        );
        assert_eq!(result, "chr1:g.100_101insTCG");
    }

    #[test]
    fn test_hgvsg_mnv() {
        let result = hgvsg(
            "chr1",
            100,
            101,
            &Allele::Sequence(b"AC".to_vec()),
            &Allele::Sequence(b"GT".to_vec()),
        );
        assert_eq!(result, "chr1:g.100_101delinsGT");
    }
}
