//! Long variant encoding using 2-bit-per-base packing into u32 words.
//!
//! For variants where ref+alt exceeds 4 bases (the Var32 limit), we pack
//! the allele sequences into a vector of u32 words at 16 bases per word.

use anyhow::Result;
use serde::{Deserialize, Serialize};

/// A variant too long for Var32 encoding.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LongVariant {
    /// Genomic position (full, not within-chunk).
    pub position: u32,
    /// Index into the parallel value arrays for this chunk.
    pub idx: u32,
    /// Encoded allele sequence: [ref_len, alt_len, packed_bases...]
    pub sequence: Vec<u32>,
}

impl PartialEq for LongVariant {
    fn eq(&self, other: &Self) -> bool {
        self.position == other.position && self.sequence == other.sequence
    }
}

impl Eq for LongVariant {}

impl PartialOrd for LongVariant {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for LongVariant {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.position
            .cmp(&other.position)
            .then_with(|| self.sequence.cmp(&other.sequence))
    }
}

/// A variant whose alleles cannot be 2-bit packed at all, because they contain
/// a byte outside {A,C,G,T} — `N` runs and IUPAC ambiguity codes, which real
/// ClinVar carries on ~0.015% of records.
///
/// Neither Var32 nor [`LongVariant`] can represent these: both encodings are
/// 2 bits per base by construction. Rather than dropping them (which is what
/// the writer used to do, with a `log::warn!` that no CLI logger was installed
/// to show), a chunk keeps them in a third, small bucket that stores the allele
/// bytes verbatim and is binary-searched on `(position, ref, alt)`.
///
/// Ordering is `(position, ref_allele, alt_allele)` so the writer's sort and the
/// reader's `binary_search` agree; `idx` is excluded from comparison for the
/// same reason it is in `LongVariant` — the query side does not know it.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OtherVariant {
    /// Genomic position (full, not within-chunk).
    pub position: u32,
    /// Index into the parallel value arrays for this chunk.
    pub idx: u32,
    /// REF allele bytes, verbatim.
    pub ref_allele: Vec<u8>,
    /// ALT allele bytes, verbatim.
    pub alt_allele: Vec<u8>,
}

impl PartialEq for OtherVariant {
    fn eq(&self, other: &Self) -> bool {
        self.position == other.position
            && self.ref_allele == other.ref_allele
            && self.alt_allele == other.alt_allele
    }
}

impl Eq for OtherVariant {}

impl PartialOrd for OtherVariant {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for OtherVariant {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.position
            .cmp(&other.position)
            .then_with(|| self.ref_allele.cmp(&other.ref_allele))
            .then_with(|| self.alt_allele.cmp(&other.alt_allele))
    }
}

/// Whether both alleles are 2-bit packable, i.e. contain only `ACGTacgt`.
///
/// The single predicate both the writer and the reader use to decide whether a
/// variant belongs in the Var32/kmer16 index or in the verbatim
/// [`OtherVariant`] bucket. Keeping it in one place is what stops the two sides
/// from disagreeing about where a variant was filed.
#[inline]
pub fn is_acgt_only(ref_allele: &[u8], alt_allele: &[u8]) -> bool {
    ref_allele
        .iter()
        .chain(alt_allele.iter())
        .all(|b| base_to_bits(*b).is_some())
}

/// DNA base to 2-bit encoding. Returns `None` for non-ACGT bytes so callers
/// can surface the problem instead of silently encoding `N`/IUPAC as 'T'.
#[inline]
fn base_to_bits(b: u8) -> Option<u32> {
    match b {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

/// Encode ref and alt alleles into a kmer16 sequence vector.
///
/// Format: `[ref_len as u32, alt_len as u32, packed_bases...]`
/// Each u32 word holds 16 bases at 2 bits each.
///
/// Returns `None` if any base is not in {A,C,G,T,a,c,g,t}. Earlier revisions
/// silently mapped non-ACGT bytes to 'T', which caused encode/decode round
/// trips to diverge from the input sequence.
pub fn encode_var(ref_allele: &[u8], alt_allele: &[u8]) -> Option<Vec<u32>> {
    let total_bases = ref_allele.len() + alt_allele.len();
    let num_words = total_bases.div_ceil(16);

    let mut result = Vec::with_capacity(2 + num_words);
    result.push(ref_allele.len() as u32);
    result.push(alt_allele.len() as u32);

    let mut word: u32 = 0;
    let mut bit_pos: u32 = 0;

    for &b in ref_allele.iter().chain(alt_allele.iter()) {
        let bits = base_to_bits(b)?;
        word |= bits << bit_pos;
        bit_pos += 2;
        if bit_pos >= 32 {
            result.push(word);
            word = 0;
            bit_pos = 0;
        }
    }
    if bit_pos > 0 {
        result.push(word);
    }

    Some(result)
}

/// Decode a kmer16 sequence vector back to (ref_allele, alt_allele).
///
/// Returns `Err` if the sequence is malformed or claims a length that
/// exceeds the bases packed into its trailing words. Earlier revisions
/// silently returned empty vectors here, which produced false-negative
/// lookups against malformed chunks instead of surfacing the corruption.
pub fn decode_var(sequence: &[u32]) -> Result<(Vec<u8>, Vec<u8>)> {
    if sequence.len() < 2 {
        anyhow::bail!("kmer16 sequence too short: need ref_len/alt_len header");
    }
    let ref_len = sequence[0] as usize;
    let alt_len = sequence[1] as usize;
    let total = ref_len
        .checked_add(alt_len)
        .ok_or_else(|| anyhow::anyhow!("kmer16 ref_len + alt_len overflows usize"))?;

    // Each trailing word holds 16 bases; verify the encoded sequence has
    // capacity for the claimed total before decoding.
    let max_bases = sequence.len().saturating_sub(2).saturating_mul(16);
    if total > max_bases {
        anyhow::bail!(
            "kmer16 sequence claims {} bases but only {} fit in payload",
            total,
            max_bases
        );
    }

    let bases_decode = *b"ACGT";
    let mut all_bases = Vec::with_capacity(total);

    let mut count = 0;
    'outer: for &word in &sequence[2..] {
        for shift in (0..32).step_by(2) {
            if count >= total {
                break 'outer;
            }
            let idx = ((word >> shift) & 0x3) as usize;
            all_bases.push(bases_decode[idx]);
            count += 1;
        }
    }

    if all_bases.len() < total {
        anyhow::bail!(
            "kmer16 decoded only {} of {} expected bases",
            all_bases.len(),
            total
        );
    }

    let ref_allele = all_bases[..ref_len].to_vec();
    let alt_allele = all_bases[ref_len..ref_len + alt_len].to_vec();
    Ok((ref_allele, alt_allele))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_encode_decode_short() {
        let seq = encode_var(b"ACGT", b"TGCA").unwrap();
        let (r, a) = decode_var(&seq).unwrap();
        assert_eq!(r, b"ACGT");
        assert_eq!(a, b"TGCA");
    }

    #[test]
    fn test_encode_decode_long() {
        let ref_allele = b"ACGTACGTACGTACGT"; // 16 bases
        let alt_allele = b"TGCATGCATGCATGCA"; // 16 bases
        let seq = encode_var(ref_allele, alt_allele).unwrap();
        let (r, a) = decode_var(&seq).unwrap();
        assert_eq!(r, ref_allele);
        assert_eq!(a, alt_allele);
    }

    #[test]
    fn test_long_variant_ordering() {
        let a = LongVariant {
            position: 100,
            idx: 0,
            sequence: encode_var(b"ACGTAC", b"T").unwrap(),
        };
        let b = LongVariant {
            position: 200,
            idx: 1,
            sequence: encode_var(b"ACGTAC", b"T").unwrap(),
        };
        let c = LongVariant {
            position: 100,
            idx: 2,
            sequence: encode_var(b"TTTTT", b"A").unwrap(),
        };
        assert!(a < b); // position 100 < 200
        assert!(a != c); // same position, different sequence
    }

    #[test]
    fn test_equality_ignores_idx() {
        let a = LongVariant { position: 100, idx: 0, sequence: encode_var(b"ACGTAC", b"T").unwrap() };
        let b = LongVariant { position: 100, idx: 999, sequence: encode_var(b"ACGTAC", b"T").unwrap() };
        assert_eq!(a, b); // idx is ignored in PartialEq
    }

    #[test]
    fn test_encode_var_rejects_non_acgt() {
        assert!(encode_var(b"ACGTN", b"A").is_none());
        assert!(encode_var(b"A", b"ACGTN").is_none());
        assert!(encode_var(b"R", b"Y").is_none()); // IUPAC
        assert!(encode_var(b"acgt", b"A").is_some()); // lowercase ACGT ok
    }

    #[test]
    fn test_decode_var_errors_on_undersized_payload() {
        // Claims 32 bases but provides only one trailing word (16 bases).
        let bogus = vec![16u32, 16u32, 0u32];
        let err = decode_var(&bogus).unwrap_err();
        assert!(err.to_string().contains("only 16 fit"));
    }

    #[test]
    fn test_decode_var_errors_on_truncated_header() {
        assert!(decode_var(&[]).is_err());
        assert!(decode_var(&[5u32]).is_err());
    }
}
