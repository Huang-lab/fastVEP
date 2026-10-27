//! Mitochondrial genome-specific handling.
//!
//! Handles circular coordinate wrapping and the mitochondrial codon table.

use crate::CodonTable;

/// Length of the human mitochondrial genome (rCRS reference).
pub const MT_LENGTH: u64 = 16569;

/// Returns true if the chromosome name indicates mitochondrial DNA.
///
/// Asked of every FASTA fetch and every codon table lookup, so it compares without making a
/// lower-case copy of the name: the four names are ASCII, and no other string lower-cases to one.
pub fn is_mitochondrial(chrom: &str) -> bool {
    ["mt", "chrm", "chrmt", "m"]
        .iter()
        .any(|name| chrom.eq_ignore_ascii_case(name))
}

/// Wrap a position around a circular genome of the given `length`.
///
/// Generalization of [`wrap_position`] for circular contigs whose length
/// isn't the human rCRS (e.g. a non-human mitochondrial genome loaded from a
/// FASTA with a different `MT`/`chrM` record length). Positions beyond
/// `length` wrap back to the beginning; `pos == 0` or `length == 0` are
/// returned unchanged (there's no sensible wrap for a zero-length contig or
/// the sentinel 0 position).
pub fn wrap_position_for(pos: u64, length: u64) -> u64 {
    if pos == 0 || length == 0 {
        return pos;
    }
    ((pos - 1) % length) + 1
}

/// Wrap a position around the circular mitochondrial genome (human rCRS length).
/// Positions > MT_LENGTH wrap to the beginning.
pub fn wrap_position(pos: u64) -> u64 {
    wrap_position_for(pos, MT_LENGTH)
}

/// The vertebrate mitochondrial codon table (NCBI translation table 2).
///
/// Differences from standard table:
/// - AGA, AGG = Stop (not Arg)
/// - ATA = Met (not Ile)
/// - TGA = Trp (not Stop)
pub fn mitochondrial_codon_table() -> CodonTable {
    CodonTable::from_ncbi_table(2)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_is_mitochondrial() {
        assert!(is_mitochondrial("MT"));
        assert!(is_mitochondrial("chrM"));
        assert!(is_mitochondrial("chrMT"));
        assert!(is_mitochondrial("M"));
        assert!(!is_mitochondrial("chr1"));
        assert!(!is_mitochondrial("chrX"));
    }

    /// What it answered when it lower-cased a copy of the name, for every spelling of the four
    /// names in either case and for names that are not one.
    #[test]
    fn a_name_is_mitochondrial_as_it_was_when_a_lower_case_copy_was_compared() {
        let before = |chrom: &str| {
            let c = chrom.to_lowercase();
            c == "mt" || c == "chrm" || c == "chrmt" || c == "m"
        };
        let mut names: Vec<String> = Vec::new();
        for base in ["mt", "chrm", "chrmt", "m"] {
            // every case of every letter
            for mask in 0..(1u32 << base.len()) {
                names.push(
                    base.chars()
                        .enumerate()
                        .map(|(i, c)| {
                            if mask >> i & 1 == 1 {
                                c.to_ascii_uppercase()
                            } else {
                                c
                            }
                        })
                        .collect(),
                );
            }
        }
        for other in [
            "",
            "chr1",
            "1",
            "X",
            "chrX",
            "MTx",
            "xMT",
            "chr_M",
            "M ",
            " M",
            "chrUn_KI270442v1",
            "NC_012920.1",
            "\u{ff2d}\u{ff34}",
            "\u{212a}",
            "m\u{0307}",
            "chrm\u{0131}",
        ] {
            names.push(other.to_string());
        }
        for name in &names {
            assert_eq!(is_mitochondrial(name), before(name), "{name:?}");
        }
        assert!(names.iter().filter(|n| is_mitochondrial(n)).count() > 30);
    }

    #[test]
    fn test_wrap_position() {
        assert_eq!(wrap_position(1), 1);
        assert_eq!(wrap_position(16569), 16569);
        assert_eq!(wrap_position(16570), 1);
        assert_eq!(wrap_position(16571), 2);
        assert_eq!(wrap_position(0), 0);
    }

    #[test]
    fn test_wrap_position_for_arbitrary_length() {
        // Same shape as wrap_position, but parameterized for organisms whose
        // mitochondrial contig isn't the human rCRS length (e.g. mouse mtDNA,
        // ~16299 bp).
        assert_eq!(wrap_position_for(1, 100), 1);
        assert_eq!(wrap_position_for(100, 100), 100);
        assert_eq!(wrap_position_for(101, 100), 1);
        assert_eq!(wrap_position_for(105, 100), 5);
        assert_eq!(wrap_position_for(0, 100), 0);
        assert_eq!(wrap_position_for(5, 0), 5);
        // Delegation: wrap_position(pos) == wrap_position_for(pos, MT_LENGTH).
        assert_eq!(wrap_position_for(16570, MT_LENGTH), wrap_position(16570));
    }

    #[test]
    fn test_mt_codon_table() {
        let table = mitochondrial_codon_table();
        // ATA = Met in MT (Ile in standard)
        assert_eq!(table.translate(b"ATA"), b'M');
        // TGA = Trp in MT (Stop in standard)
        assert_eq!(table.translate(b"TGA"), b'W');
        // AGA = Stop in MT (Arg in standard)
        assert_eq!(table.translate(b"AGA"), b'*');
        // Normal codons unchanged
        assert_eq!(table.translate(b"ATG"), b'M');
        assert_eq!(table.translate(b"TAA"), b'*');
    }
}
