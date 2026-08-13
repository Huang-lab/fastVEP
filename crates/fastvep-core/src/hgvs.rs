//! Small parsers over HGVS `c.` notation that more than one crate needs.
//!
//! These live here rather than beside their first caller because both the
//! classifier (which reads an offset off the transcript annotation it is
//! handed) and the supplementary-annotation builders (which read one off a
//! ClinVar `Name` column) have to agree exactly on what `c.376-2` means. Two
//! copies of a parser this fiddly drift.

/// Signed offset in bp from the nearest exon boundary, parsed from an HGVS
/// `c.` expression. Positive after a donor (`c.123+5` → `+5`), negative before
/// an acceptor (`c.123-15` → `-15`).
///
/// The offset is the HGVS `+N` / `-N` token that *follows* the CDS position
/// number, so it is always immediately preceded by a digit. That distinguishes
/// it from the leading sign of a UTR position (`c.-23…`, where the `-` is
/// preceded by `.`) and from transcript-version dots (`ENST….7:c.…`). For range
/// variants (e.g. `c.4001+12_4001+15del`) the endpoint nearest the boundary
/// (smallest |offset|) is returned. Returns `None` for purely exonic variants
/// (no such token) - e.g. `c.5098G>C`, `c.*1411T>A`.
pub fn parse_intronic_offset(hgvs_c: &str) -> Option<i64> {
    let bytes = hgvs_c.as_bytes();
    let mut best: Option<i64> = None;
    let mut i = 0;
    while i < bytes.len() {
        let b = bytes[i];
        if (b == b'+' || b == b'-') && i > 0 && bytes[i - 1].is_ascii_digit() {
            let sign = if b == b'-' { -1i64 } else { 1i64 };
            let mut j = i + 1;
            let mut val: i64 = 0;
            let mut any = false;
            while j < bytes.len() && bytes[j].is_ascii_digit() {
                val = val.saturating_mul(10).saturating_add((bytes[j] - b'0') as i64);
                any = true;
                j += 1;
            }
            if any {
                let off = sign * val;
                best = Some(match best {
                    Some(prev) if prev.abs() <= off.abs() => prev,
                    _ => off,
                });
            }
            i = j;
        } else {
            i += 1;
        }
    }
    best
}

/// Whether an intronic offset falls on one of the two canonical splice
/// dinucleotide bases, i.e. `+1`, `+2`, `-1` or `-2`.
pub fn is_canonical_dinucleotide_offset(offset: i64) -> bool {
    (1..=2).contains(&offset.unsigned_abs())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_intronic_offset() {
        // Exonic: no `+N`/`-N` token following the CDS number.
        assert_eq!(parse_intronic_offset("ENST00000272371.7:c.5098G>C"), None);
        assert_eq!(parse_intronic_offset("c.*1411T>A"), None);
        // Canonical splice positions.
        assert_eq!(parse_intronic_offset("c.964+1G>A"), Some(1));
        assert_eq!(parse_intronic_offset("ENST00000378156.9:c.2818-2A>."), Some(-2));
        // Deep intronic.
        assert_eq!(parse_intronic_offset("c.4001+12_4001+15del"), Some(12));
        assert_eq!(parse_intronic_offset("n.162-24414C>T"), Some(-24414));
        // Ranges return the endpoint nearest the boundary.
        assert_eq!(parse_intronic_offset("c.366-1_366+2dup"), Some(-1));
        assert_eq!(parse_intronic_offset("c.541-30_541-2dup"), Some(-2));
        // A leading UTR sign is not an offset; the trailing one is.
        assert_eq!(parse_intronic_offset("c.-23+1G>A"), Some(1));
        assert_eq!(parse_intronic_offset("c.-23-1G>A"), Some(-1));
    }

    #[test]
    fn test_is_canonical_dinucleotide_offset() {
        for off in [1, 2, -1, -2] {
            assert!(is_canonical_dinucleotide_offset(off), "{off}");
        }
        for off in [0, 3, -3, 7, -21] {
            assert!(!is_canonical_dinucleotide_offset(off), "{off}");
        }
    }
}
