//! Long variant encoding using 2-bit-per-base packing into u32 words.
//!
//! For variants where ref+alt exceeds 4 bases (the Var32 limit), we pack
//! the allele sequences into a vector of u32 words at 16 bases per word.

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

/// DNA base to 2-bit encoding.
#[inline]
fn base_to_bits(b: u8) -> u32 {
    match b {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        _ => 3, // T and anything else
    }
}

/// Encode ref and alt alleles into a kmer16 sequence vector.
///
/// Format: `[ref_len as u32, alt_len as u32, packed_bases...]`
/// Each u32 word holds 16 bases at 2 bits each.
pub fn encode_var(ref_allele: &[u8], alt_allele: &[u8]) -> Vec<u32> {
    let total_bases = ref_allele.len() + alt_allele.len();
    let num_words = (total_bases + 15) / 16;

    let mut result = Vec::with_capacity(2 + num_words);
    result.push(ref_allele.len() as u32);
    result.push(alt_allele.len() as u32);

    let mut word: u32 = 0;
    let mut bit_pos: u32 = 0;

    for &b in ref_allele.iter().chain(alt_allele.iter()) {
        word |= base_to_bits(b) << bit_pos;
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

    result
}

/// Decode a kmer16 sequence vector back to (ref_allele, alt_allele).
///
/// Returns empty vectors if the sequence is malformed or claims a length that
/// exceeds the bases packed into its trailing words.
pub fn decode_var(sequence: &[u32]) -> (Vec<u8>, Vec<u8>) {
    if sequence.len() < 2 {
        return (Vec::new(), Vec::new());
    }
    let ref_len = sequence[0] as usize;
    let alt_len = sequence[1] as usize;
    let total = match ref_len.checked_add(alt_len) {
        Some(t) => t,
        None => return (Vec::new(), Vec::new()),
    };

    // Each trailing word holds 16 bases; verify the encoded sequence has
    // capacity for the claimed total before decoding.
    let max_bases = sequence.len().saturating_sub(2).saturating_mul(16);
    if total > max_bases {
        return (Vec::new(), Vec::new());
    }

    let bases_decode = [b'A', b'C', b'G', b'T'];
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
        return (Vec::new(), Vec::new());
    }

    let ref_allele = all_bases[..ref_len].to_vec();
    let alt_allele = all_bases[ref_len..ref_len + alt_len].to_vec();
    (ref_allele, alt_allele)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_encode_decode_short() {
        let seq = encode_var(b"ACGT", b"TGCA");
        let (r, a) = decode_var(&seq);
        assert_eq!(r, b"ACGT");
        assert_eq!(a, b"TGCA");
    }

    #[test]
    fn test_encode_decode_long() {
        let ref_allele = b"ACGTACGTACGTACGT"; // 16 bases
        let alt_allele = b"TGCATGCATGCATGCA"; // 16 bases
        let seq = encode_var(ref_allele, alt_allele);
        let (r, a) = decode_var(&seq);
        assert_eq!(r, ref_allele);
        assert_eq!(a, alt_allele);
    }

    #[test]
    fn test_long_variant_ordering() {
        let a = LongVariant {
            position: 100,
            idx: 0,
            sequence: encode_var(b"ACGTAC", b"T"),
        };
        let b = LongVariant {
            position: 200,
            idx: 1,
            sequence: encode_var(b"ACGTAC", b"T"),
        };
        let c = LongVariant {
            position: 100,
            idx: 2,
            sequence: encode_var(b"TTTTT", b"A"),
        };
        assert!(a < b); // position 100 < 200
        assert!(a != c); // same position, different sequence
    }

    #[test]
    fn test_equality_ignores_idx() {
        let a = LongVariant { position: 100, idx: 0, sequence: encode_var(b"ACGTAC", b"T") };
        let b = LongVariant { position: 100, idx: 999, sequence: encode_var(b"ACGTAC", b"T") };
        assert_eq!(a, b); // idx is ignored in PartialEq
    }
}
