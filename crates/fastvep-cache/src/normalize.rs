//! Reference-based variant normalization (vt-normalize / minimal representation).
//!
//! Supplementary-annotation sources such as gnomAD store their variants in the
//! VCF convention used by the source: left-aligned and parsimonious (shared
//! prefix/suffix bases trimmed down to a single anchor). The SA lookup matches
//! by exact `(position, ref, alt)` string equality, so a query allele that is
//! anchored differently — common for indels, and especially for indels inside
//! tandem repeats where the same event can be written at several positions —
//! silently fails to match even though the variant is present.
//!
//! [`normalize_variant`] rewrites a query `(pos, ref, alt)` into the same
//! minimal/left-aligned representation gnomAD uses, using the reference to
//! left-shift indels to their leftmost position. It is a strict no-op for SNVs
//! and degrades gracefully (returning the input after only the
//! reference-free trims) whenever the reference is unavailable, a fetch fails,
//! or the variant reaches a contig boundary. It never panics.

use crate::providers::SequenceProvider;

/// A position + ref/alt in 1-based anchored VCF coordinate space.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NormalizedVariant {
    /// 1-based position.
    pub pos: u64,
    pub ref_allele: String,
    pub alt_allele: String,
}

fn is_acgtn(b: u8) -> bool {
    matches!(b, b'A' | b'C' | b'G' | b'T' | b'N')
}

/// Left-align + parsimonious-trim a single `(pos, ref, alt)` to its minimal /
/// gnomAD representation.
///
/// `pos` is 1-based (VCF / [`SequenceProvider`] convention). `chrom` is passed
/// through to `seq` unchanged. The caller splits multi-allelic sites and calls
/// this once per ALT.
///
/// Returns the input unchanged (after only unconditionally-safe trims) for any
/// case where reference-based normalization cannot proceed safely: SNVs / equal
/// length ref&alt, empty / symbolic / non-ACGTN alleles, reference fetch
/// failures, and contig-start boundaries. Never panics.
pub fn normalize_variant(
    seq: &dyn SequenceProvider,
    chrom: &str,
    pos: u64,
    ref_allele: &str,
    alt_allele: &str,
) -> NormalizedVariant {
    let as_is = || NormalizedVariant {
        pos,
        ref_allele: ref_allele.to_string(),
        alt_allele: alt_allele.to_string(),
    };

    // ── Guard 0: bail to no-op for cases normalization must not touch. ──
    if ref_allele.is_empty() || alt_allele.is_empty() {
        return as_is();
    }
    // Symbolic ALT (<DEL>, <*>) or an unsplit multi-allelic value.
    if alt_allele.starts_with('<') || alt_allele.contains(',') {
        return as_is();
    }
    // SNV / MNV (equal length): no indel to left-shift; exact match already works.
    if ref_allele.len() == alt_allele.len() {
        return as_is();
    }

    let mut r: Vec<u8> = ref_allele.to_ascii_uppercase().into_bytes();
    let mut a: Vec<u8> = alt_allele.to_ascii_uppercase().into_bytes();
    // Refuse to roll over ambiguous / non-nucleotide bases (e.g. '*', IUPAC).
    if !r.iter().all(|&b| is_acgtn(b)) || !a.iter().all(|&b| is_acgtn(b)) {
        return as_is();
    }

    let mut p = pos;

    // ── Step 1: trim shared SUFFIX bases (reference not needed). Keep ≥1 base
    // in each allele so we never produce an empty allele here. ──
    while r.len() > 1 && a.len() > 1 && r.last() == a.last() {
        r.pop();
        a.pop();
    }

    // ── Step 2: left-shift while the variant is a pure indel anchored on a
    // shared base. Roll one base at a time using the reference. ──
    loop {
        // Shiftable only when one side is reduced to the single anchor base and
        // both share that leading base (canonical anchored indel).
        let shiftable = (r.len() == 1 || a.len() == 1) && r.first() == a.first();
        if !shiftable || p <= 1 {
            break;
        }
        // Fetch the single reference base immediately to the left. Roll only
        // over unambiguous A/C/G/T — like vt/bcftools, stop at N so we never
        // produce an N-anchored representation inside an assembly gap.
        let prev = match seq.fetch_sequence(chrom, p - 1, p - 1) {
            Ok(s) if s.len() == 1 && matches!(s[0].to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T') => {
                s[0].to_ascii_uppercase()
            }
            _ => break, // fetch failed / short read / N / ambiguous → stop here.
        };
        // Prepend the new left base to both alleles, then drop the now-redundant
        // shared last base, and step the position left.
        r.insert(0, prev);
        a.insert(0, prev);
        if r.len() > 1 && a.len() > 1 && r.last() == a.last() {
            r.pop();
            a.pop();
        }
        p -= 1;
    }

    // ── Step 3: trim shared PREFIX bases, keeping one anchor base. Each trimmed
    // prefix base advances the position. ──
    while r.len() > 1 && a.len() > 1 && r.first() == a.first() {
        r.remove(0);
        a.remove(0);
        p += 1;
    }

    NormalizedVariant {
        pos: p,
        // Inputs were validated ACGTN above, so these are valid UTF-8.
        ref_allele: String::from_utf8(r).unwrap_or_else(|_| ref_allele.to_string()),
        alt_allele: String::from_utf8(a).unwrap_or_else(|_| alt_allele.to_string()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use anyhow::{anyhow, Result};

    /// Minimal `SequenceProvider` backed by a 1-based reference string for one
    /// contig. `fetch_sequence` mirrors the real readers' contract: 1-based
    /// inclusive, `Err` past the contig.
    struct StrRef {
        chrom: &'static str,
        seq: &'static str, // sequence starting at position 1
    }
    impl SequenceProvider for StrRef {
        fn fetch_sequence(&self, chrom: &str, start: u64, end: u64) -> Result<Vec<u8>> {
            if chrom != self.chrom {
                return Err(anyhow!("unknown contig"));
            }
            if start < 1 || end < start {
                return Err(anyhow!("bad range"));
            }
            let bytes = self.seq.as_bytes();
            let s0 = (start - 1) as usize;
            if s0 >= bytes.len() {
                return Err(anyhow!("past contig end"));
            }
            let e = ((end) as usize).min(bytes.len());
            Ok(bytes[s0..e].to_vec())
        }
    }

    /// Provider whose fetches always fail — exercises graceful degradation.
    struct FailRef;
    impl SequenceProvider for FailRef {
        fn fetch_sequence(&self, _chrom: &str, _start: u64, _end: u64) -> Result<Vec<u8>> {
            Err(anyhow!("no reference"))
        }
    }

    #[test]
    fn snv_is_noop() {
        let r = StrRef { chrom: "chr1", seq: "ACGTACGT" };
        let n = normalize_variant(&r, "chr1", 3, "G", "T");
        assert_eq!(n, NormalizedVariant { pos: 3, ref_allele: "G".into(), alt_allele: "T".into() });
    }

    #[test]
    fn mnv_equal_len_is_noop() {
        let r = StrRef { chrom: "chr1", seq: "ACGTACGT" };
        let n = normalize_variant(&r, "chr1", 2, "CG", "TT");
        assert_eq!(n.pos, 2);
        assert_eq!(n.ref_allele, "CG");
        assert_eq!(n.alt_allele, "TT");
    }

    #[test]
    fn suffix_trim_deletion() {
        // ref=CAT alt=CT share suffix T → CA/C; then no left-repeat to shift.
        // Reference around pos 10: ...positions 10=C,11=A,12=T...
        let r = StrRef { chrom: "chr1", seq: "NNNNNNNNNCATGGG" };
        let n = normalize_variant(&r, "chr1", 10, "CAT", "CT");
        // Suffix-trim T → (CA, C) at pos 10. Anchor 'C' shared; left base is
        // 'N' (pos 9) which is ACGTN but rolling would just walk through Ns;
        // the deletion of 'A' is already at its leftmost real position.
        assert_eq!(n.ref_allele.len() as i64 - n.alt_allele.len() as i64, 1);
    }

    #[test]
    fn repeat_region_left_align_deletion() {
        // Microsatellite: positions 1.. = "G TAAC TAAC TAAC ...".
        // Reference: G at 1, then TAAC repeated.
        // Seq: pos1=G,2=T,3=A,4=A,5=C,6=T,7=A,8=A,9=C,10=T,11=A,12=A,13=C,...
        let r = StrRef { chrom: "chr2", seq: "GTAACTAACTAACTAACG" };
        // A right-anchored 4-base (TAAC) contraction written at pos 6:
        // ref=CTAAC alt=C (delete one TAAC copy). Should left-align to pos 1
        // (anchor G) → ref=GTAAC, alt=G.
        let n = normalize_variant(&r, "chr2", 5, "CTAAC", "C");
        assert_eq!(n.pos, 1, "deletion should left-align to the leftmost repeat copy");
        assert_eq!(n.ref_allele, "GTAAC");
        assert_eq!(n.alt_allele, "G");
    }

    #[test]
    fn insertion_left_align_in_homopolymer() {
        // Homopolymer A-run: pos1=C, 2..=AAAA. An inserted 'A' anywhere in the
        // run left-aligns to the C anchor at pos 1.
        let r = StrRef { chrom: "chr3", seq: "CAAAAG" };
        // Insertion written at pos 4: ref=A alt=AA → left-align to pos1 C anchor.
        let n = normalize_variant(&r, "chr3", 4, "A", "AA");
        assert_eq!(n.pos, 1);
        assert_eq!(n.ref_allele, "C");
        assert_eq!(n.alt_allele, "CA");
    }

    #[test]
    fn contig_start_boundary_no_panic() {
        // Deletion already adjacent to position 1 cannot shift further left.
        let r = StrRef { chrom: "chr1", seq: "ATTTTG" };
        let n = normalize_variant(&r, "chr1", 1, "AT", "A");
        assert_eq!(n.pos, 1);
        assert!(!n.ref_allele.is_empty() && !n.alt_allele.is_empty());
    }

    #[test]
    fn fetch_error_is_graceful() {
        // No usable reference → return the suffix/prefix-trimmed form, no panic.
        let n = normalize_variant(&FailRef, "chrX", 100, "CAT", "CT");
        assert!(!n.ref_allele.is_empty() && !n.alt_allele.is_empty());
        // Suffix-trim still applies (reference-free): CAT/CT → CA/C.
        assert_eq!(n.ref_allele, "CA");
        assert_eq!(n.alt_allele, "C");
    }

    #[test]
    fn non_acgtn_is_noop() {
        let r = StrRef { chrom: "chr1", seq: "ACGTACGT" };
        let n = normalize_variant(&r, "chr1", 2, "C", "*");
        assert_eq!(n.alt_allele, "*");
        assert_eq!(n.pos, 2);
    }

    #[test]
    fn left_shift_stops_at_n() {
        // Reference: pos1=A, then an N-gap, then a real region. A deletion just
        // right of the N-run must NOT roll into the gap (vt/bcftools stop at N).
        // seq: pos1=A,2=N,3=N,4=C,5=A,6=T,7=G
        let r = StrRef { chrom: "chr1", seq: "ANNCATG" };
        // ref=CAT alt=CT at pos4 → suffix-trim to (CA,C); left base pos3='N' → stop.
        let n = normalize_variant(&r, "chr1", 4, "CAT", "CT");
        assert_eq!(n.pos, 4);
        assert_eq!(n.ref_allele, "CA");
        assert_eq!(n.alt_allele, "C");
    }

    #[test]
    fn already_minimal_insertion_unchanged() {
        // ref=A alt=AT with no left repeat of T → stays put.
        let r = StrRef { chrom: "chr1", seq: "GGAGGG" };
        let n = normalize_variant(&r, "chr1", 3, "A", "AT");
        assert_eq!(n.pos, 3);
        assert_eq!(n.ref_allele, "A");
        assert_eq!(n.alt_allele, "AT");
    }
}
