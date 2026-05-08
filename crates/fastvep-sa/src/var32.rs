//! Compact 32-bit variant encoding for fast binary search.
//!
//! Packs a variant's within-chunk position, reference allele length, alternate
//! allele length, and encoded DNA bases into a single `u32`. This enables
//! binary search on sorted arrays of 4-byte keys instead of variable-length
//! string comparisons.
//!
//! Layout (32 bits):
//! ```text
//! [position: 20 bits][rlen: 2 bits][alen: 2 bits][enc: 8 bits]
//! ```
//!
//! - `position`: Within-chunk offset (0 to 2^20 - 1 = 1,048,575)
//! - `rlen`: Reference allele length minus 1 (0-3 → 1-4 bases)
//! - `alen`: Alternate allele length minus 1 (0-3 → 1-4 bases)
//! - `enc`: Up to 4 DNA bases packed at 2 bits each (A=0, C=1, G=2, T=3)
//!
//! Variants with combined ref+alt length > 4 bases are "long" and stored
//! separately via [`super::kmer16`].

/// Number of bits used for the within-chunk position.
pub const CHUNK_BITS: u32 = 20;

/// Maximum combined ref+alt bases that fit in a Var32.
pub const MAX_COMBINED_LEN: usize = 4;

/// DNA base to 2-bit encoding lookup (indexed by ASCII byte).
const BASE_ENC: [u8; 128] = {
    let mut table = [3u8; 128]; // default T
    table[b'A' as usize] = 0;
    table[b'a' as usize] = 0;
    table[b'C' as usize] = 1;
    table[b'c' as usize] = 1;
    table[b'G' as usize] = 2;
    table[b'g' as usize] = 2;
    table[b'T' as usize] = 3;
    table[b't' as usize] = 3;
    table
};

/// 2-bit encoding back to DNA base.
const BASE_DEC: [u8; 4] = [b'A', b'C', b'G', b'T'];

/// Returns true if the variant is too long for Var32 encoding.
#[inline]
pub fn is_long(ref_len: usize, alt_len: usize) -> bool {
    ref_len + alt_len > MAX_COMBINED_LEN || ref_len == 0 || alt_len == 0 || ref_len > 4 || alt_len > 4
}

/// Encode a variant into a compact 32-bit representation.
///
/// `pos` is the within-chunk position (i.e., `genomic_pos & ((1 << CHUNK_BITS) - 1)`).
/// Returns `None` if the variant is too long for Var32.
#[inline]
pub fn encode(pos: u32, ref_allele: &[u8], alt_allele: &[u8]) -> Option<u32> {
    let rlen = ref_allele.len();
    let alen = alt_allele.len();
    if is_long(rlen, alen) {
        return None;
    }

    let mut v: u32 = 0;
    // Position in bits 12..31
    v |= (pos & 0xF_FFFF) << 12;
    // rlen-1 in bits 10..11
    v |= ((rlen as u32 - 1) & 0x3) << 10;
    // alen-1 in bits 8..9
    v |= ((alen as u32 - 1) & 0x3) << 8;

    // Encode bases: ref then alt, 2 bits each, packed into bits 0..7
    let mut enc: u8 = 0;
    let mut bit_pos = 6; // Start from high bits of the 8-bit field
    for &b in ref_allele.iter().chain(alt_allele.iter()) {
        let idx = (b as usize).min(127);
        enc |= BASE_ENC[idx] << bit_pos;
        if bit_pos >= 2 {
            bit_pos -= 2;
        }
    }
    v |= enc as u32;

    Some(v)
}

/// Decode a Var32 back to (within-chunk position, ref_allele, alt_allele).
///
/// `is_long` enforces `rlen + alen <= MAX_COMBINED_LEN` (4) at encode time,
/// so decoding consumes at most 4 base slots from the 8-bit packed field.
/// Shifts here are statically bounded to 0, 2, 4, 6.
#[inline]
pub fn decode(v: u32) -> (u32, Vec<u8>, Vec<u8>) {
    let pos = (v >> 12) & 0xF_FFFF;
    let rlen = (((v >> 10) & 0x3) + 1) as usize;
    let alen = (((v >> 8) & 0x3) + 1) as usize;
    let enc = (v & 0xFF) as u8;

    let mut ref_allele = Vec::with_capacity(rlen);
    let mut alt_allele = Vec::with_capacity(alen);

    // bit_pos goes 6, 4, 2, 0 across at most 4 bases.
    let mut shift: u32 = 6;
    let total = rlen + alen;
    let total = total.min(MAX_COMBINED_LEN);
    for i in 0..total {
        let idx = ((enc >> shift) & 0x3) as usize;
        if i < rlen {
            ref_allele.push(BASE_DEC[idx]);
        } else {
            alt_allele.push(BASE_DEC[idx]);
        }
        shift = shift.saturating_sub(2);
    }

    (pos, ref_allele, alt_allele)
}

/// Extract the within-chunk position from a genomic position.
#[inline]
pub fn chunk_position(genomic_pos: u32) -> u32 {
    genomic_pos & ((1 << CHUNK_BITS) - 1)
}

/// Extract the chunk ID from a genomic position.
#[inline]
pub fn chunk_id(genomic_pos: u32) -> u32 {
    genomic_pos >> CHUNK_BITS
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_snv_round_trip() {
        // A>G at position 12345
        let pos = chunk_position(12345);
        let encoded = encode(pos, b"A", b"G").unwrap();
        let (dec_pos, dec_ref, dec_alt) = decode(encoded);
        assert_eq!(dec_pos, pos);
        assert_eq!(dec_ref, b"A");
        assert_eq!(dec_alt, b"G");
    }

    #[test]
    fn test_all_single_base_combinations() {
        let bases = [b'A', b'C', b'G', b'T'];
        for &r in &bases {
            for &a in &bases {
                let pos = 100;
                let encoded = encode(pos, &[r], &[a]).unwrap();
                let (dp, dr, da) = decode(encoded);
                assert_eq!(dp, pos);
                assert_eq!(dr, &[r]);
                assert_eq!(da, &[a]);
            }
        }
    }

    #[test]
    fn test_multi_base() {
        // 2-base ref, 2-base alt (total 4 = max)
        let pos = 500;
        let encoded = encode(pos, b"AC", b"GT").unwrap();
        let (dp, dr, da) = decode(encoded);
        assert_eq!(dp, pos);
        assert_eq!(dr, b"AC");
        assert_eq!(da, b"GT");
    }

    #[test]
    fn test_three_plus_one() {
        // 3-base ref, 1-base alt
        let pos = 999;
        let encoded = encode(pos, b"ACG", b"T").unwrap();
        let (dp, dr, da) = decode(encoded);
        assert_eq!(dp, pos);
        assert_eq!(dr, b"ACG");
        assert_eq!(da, b"T");
    }

    #[test]
    fn test_too_long_returns_none() {
        // 3+2 = 5 bases, exceeds MAX_COMBINED_LEN
        assert!(encode(0, b"ACG", b"TT").is_none());
        // Empty alleles
        assert!(encode(0, b"", b"A").is_none());
        assert!(encode(0, b"A", b"").is_none());
    }

    #[test]
    fn test_sorting_order() {
        // Var32 values should sort by position first (position is in high bits)
        let v1 = encode(100, b"A", b"G").unwrap();
        let v2 = encode(200, b"A", b"G").unwrap();
        let v3 = encode(100, b"A", b"T").unwrap();
        assert!(v1 < v2); // position 100 < 200
        // Same position: sorted by allele encoding
        assert!(v1 != v3);
    }

    #[test]
    fn test_round_trip_all_combined_lengths() {
        // Exhaustively round-trip every (rlen, alen) shape allowed by the
        // encoder using a representative base for each slot. Catches the
        // historic bug where decoding the 4th base used a stale shift.
        for rlen in 1..=4usize {
            for alen in 1..=4usize {
                if rlen + alen > MAX_COMBINED_LEN { continue; }
                let r: Vec<u8> = (0..rlen).map(|i| BASE_DEC[i % 4]).collect();
                let a: Vec<u8> = (0..alen).map(|i| BASE_DEC[(i + 1) % 4]).collect();
                let encoded = encode(7, &r, &a)
                    .unwrap_or_else(|| panic!("encode failed for rlen={}, alen={}", rlen, alen));
                let (dp, dr, da) = decode(encoded);
                assert_eq!(dp, 7);
                assert_eq!(dr, r, "ref mismatch for rlen={} alen={}", rlen, alen);
                assert_eq!(da, a, "alt mismatch for rlen={} alen={}", rlen, alen);
            }
        }
    }

    #[test]
    fn test_chunk_id_and_position() {
        let genomic = 1_500_000u32;
        let cid = chunk_id(genomic);
        let cpos = chunk_position(genomic);
        assert_eq!(cid, 1); // 1_500_000 >> 20 = 1
        assert_eq!(cpos, 1_500_000 & 0xF_FFFF);
        // Reconstruct
        assert_eq!((cid << CHUNK_BITS) | cpos, genomic);
    }
}
