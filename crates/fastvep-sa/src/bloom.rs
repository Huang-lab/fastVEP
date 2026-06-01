//! Simple Bloom filter for fast negative lookups on genomic positions.
//!
//! Before loading and searching a chunk, check the Bloom filter to skip
//! decompression entirely when a position is definitely not present.
//! False positives are acceptable (they just cause an unnecessary chunk load);
//! false negatives never occur.

/// A compact Bloom filter backed by a bit array.
pub struct BloomFilter {
    bits: Vec<u64>,
    num_bits: usize,
    num_hashes: u32,
}

impl BloomFilter {
    /// Create a Bloom filter sized for `expected_items` with target false positive rate.
    pub fn new(expected_items: usize, fp_rate: f64) -> Self {
        let num_bits = optimal_num_bits(expected_items, fp_rate).max(64);
        let num_hashes = optimal_num_hashes(expected_items, num_bits);
        let num_words = (num_bits + 63) / 64;
        Self {
            bits: vec![0u64; num_words],
            num_bits,
            num_hashes,
        }
    }

    /// Insert a position into the filter.
    pub fn insert(&mut self, position: u32) {
        for i in 0..self.num_hashes {
            let bit = self.hash(position, i) % self.num_bits;
            self.bits[bit / 64] |= 1u64 << (bit % 64);
        }
    }

    /// Check if a position might be present (true = maybe, false = definitely not).
    #[inline]
    pub fn might_contain(&self, position: u32) -> bool {
        for i in 0..self.num_hashes {
            let bit = self.hash(position, i) % self.num_bits;
            if self.bits[bit / 64] & (1u64 << (bit % 64)) == 0 {
                return false;
            }
        }
        true
    }

    /// Simple double-hashing using mixing functions.
    #[inline]
    fn hash(&self, position: u32, i: u32) -> usize {
        let h1 = fmix32(position) as usize;
        let h2 = fmix32(position.wrapping_add(0x9e3779b9)) as usize;
        h1.wrapping_add(h2.wrapping_mul(i as usize))
    }
}

/// Murmur3 finalization mix for u32.
#[inline]
fn fmix32(mut h: u32) -> u32 {
    h ^= h >> 16;
    h = h.wrapping_mul(0x85ebca6b);
    h ^= h >> 13;
    h = h.wrapping_mul(0xc2b2ae35);
    h ^= h >> 16;
    h
}

/// Maximum number of hash functions used. Capped to keep `might_contain`
/// fast and bounded.
const MAX_HASHES: u32 = 16;

/// Upper bound on bits per element. The optimal formula
/// `m = -n * ln(p) / (ln 2)^2` runs through `as f64` and loses integer
/// precision once `n` exceeds `2^53` (the f64 mantissa). For tiny `p` or
/// huge `n` the unbounded result can also be 0 or +∞. Capping the
/// bits/element ratio keeps the result in a sane range no matter what the
/// caller asks for.
///
/// 64 bits/element is generous (well below `f64`'s precision boundary and
/// far past the ~14 bits/element that already gives FPR ≈ 1e-4 at the
/// optimal hash count). A reasonable genomic deployment will never hit it;
/// the cap exists solely as a numerical-safety belt.
const MAX_BITS_PER_ELEMENT: usize = 64;

/// Headroom we hold back from `usize::MAX` so that downstream word-count
/// math like `(num_bits + 63) / 64` cannot overflow even in debug builds.
/// Choosing 1024 leaves room for any conceivable round-up by the caller.
const NUM_BITS_HEADROOM: usize = 1024;

fn optimal_num_bits(n: usize, p: f64) -> usize {
    if n == 0 {
        return 64;
    }
    // Clamp false-positive rate to a sane open interval so `p.ln()` is finite
    // and negative.
    let p = p.clamp(1e-12, 0.5);
    // `saturating_mul` can return `usize::MAX` for astronomically large `n`.
    // Subtract a small headroom so downstream callers (e.g.
    // `BloomFilter::new`, which computes `(num_bits + 63) / 64`) cannot
    // overflow on the round-up to whole 64-bit words.
    let hard_cap = n
        .saturating_mul(MAX_BITS_PER_ELEMENT)
        .min(usize::MAX - NUM_BITS_HEADROOM);
    let ln2_sq = std::f64::consts::LN_2 * std::f64::consts::LN_2;
    let bits = (-(n as f64) * p.ln() / ln2_sq).ceil();
    if !bits.is_finite() || bits <= 0.0 {
        return hard_cap.max(64);
    }
    // Compare bits as f64 against the hard cap as f64; convert back via min
    // to avoid an out-of-range f64-to-usize cast on astronomically large n.
    let capped = bits.min(hard_cap as f64);
    if capped <= 0.0 {
        64
    } else {
        capped as usize
    }
}

fn optimal_num_hashes(n: usize, m: usize) -> u32 {
    if n == 0 {
        return 1;
    }
    let k = (m as f64 / n as f64) * std::f64::consts::LN_2;
    if !k.is_finite() {
        return 1;
    }
    (k.ceil() as u32).max(1).min(MAX_HASHES)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bloom_basic() {
        let mut bloom = BloomFilter::new(1000, 0.01);

        // Insert some positions
        for i in 0..100 {
            bloom.insert(i * 100);
        }

        // All inserted positions should be found
        for i in 0..100 {
            assert!(bloom.might_contain(i * 100), "Missing position {}", i * 100);
        }

        // Most non-inserted positions should NOT be found (allow ~1% FP)
        let mut false_positives = 0;
        for i in 0..1000 {
            let pos = i * 100 + 50; // Offsets that were never inserted
            if bloom.might_contain(pos) {
                false_positives += 1;
            }
        }
        // With 1% FP rate and 1000 queries, expect ~10 false positives
        assert!(false_positives < 50, "Too many false positives: {}", false_positives);
    }

    #[test]
    fn test_bloom_empty() {
        let bloom = BloomFilter::new(100, 0.01);
        assert!(!bloom.might_contain(12345));
    }

    #[test]
    fn test_optimal_num_bits_huge_n_stays_finite() {
        // At 2^53 elements (the f64 integer-precision limit) the previous
        // formula could lose precision and yield zero. The hard cap on
        // bits/element keeps the result finite and bounded.
        let bits = optimal_num_bits(1usize << 53, 1e-12);
        assert!(bits >= 64);
        assert!(bits <= (1usize << 53).saturating_mul(MAX_BITS_PER_ELEMENT));
    }

    #[test]
    fn test_optimal_num_bits_does_not_saturate_to_usize_max() {
        // For an `n` so large that `n * MAX_BITS_PER_ELEMENT` saturates
        // `usize::MAX`, the helper must still leave room for the caller's
        // `(num_bits + 63) / 64` round-up. Earlier the saturating multiply
        // could return `usize::MAX` directly, which then overflowed in
        // `BloomFilter::new`.
        let bits = optimal_num_bits(usize::MAX, 1e-12);
        assert!(bits <= usize::MAX - NUM_BITS_HEADROOM);
        // And the headroom is sufficient for the downstream word-count math.
        let _ = bits.checked_add(63).expect("(num_bits + 63) must not overflow");
    }

    #[test]
    fn test_bloom_new_does_not_panic_for_huge_n() {
        // Smoke test: construction with absurd `n` must not panic. We don't
        // actually allocate the bit vector at this size — we just want to
        // prove the arithmetic guards work.
        let n = (1usize << 40) - 1;
        let bits = optimal_num_bits(n, 1e-9);
        let words = bits
            .checked_add(63)
            .map(|x| x / 64)
            .expect("word-count math must not overflow");
        assert!(words >= 1);
    }

    #[test]
    fn test_optimal_num_bits_basic_shape() {
        // Sanity: small n with reasonable p falls in a sensible range.
        let bits = optimal_num_bits(1000, 0.01);
        assert!(bits >= 1000 && bits <= 1000 * MAX_BITS_PER_ELEMENT);
    }
}
