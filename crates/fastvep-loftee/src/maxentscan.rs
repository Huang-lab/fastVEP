//! MaxEntScan splice site scoring.
//!
//! Ports the MaxEntScan algorithm (Yeo & Burge 2004) used by LOFTEE
//! for donor (5' splice site, 9bp) and acceptor (3' splice site, 23bp) scoring.

use std::io::{BufRead, BufReader};
use std::path::Path;

/// Pre-loaded MaxEntScan scoring tables.
pub struct MaxEntScanData {
    /// Donor context scores: 16384 entries indexed by hashseq of 7-mer (positions [0,1,2,5,6,7,8]).
    pub me2x5: Vec<f64>,
    /// Acceptor MaxEnt lookup tables (9 tables).
    pub metables: Vec<Vec<f64>>,
}

// Background nucleotide frequencies
const BGD_A: f64 = 0.27;
const BGD_C: f64 = 0.23;
const BGD_G: f64 = 0.23;
const BGD_T: f64 = 0.27;

// Donor consensus position 3 (first base of GT)
const CONS1_A: f64 = 0.004;
const CONS1_C: f64 = 0.0032;
const CONS1_G: f64 = 0.9896;
const CONS1_T: f64 = 0.0032;

// Donor consensus position 4 (second base of GT)
const CONS2_A: f64 = 0.0034;
const CONS2_C: f64 = 0.0039;
const CONS2_G: f64 = 0.0042;
const CONS2_T: f64 = 0.9884;

// Acceptor consensus position 18 (first base of AG)
const ACONS1_A: f64 = 0.9903;
const ACONS1_C: f64 = 0.0032;
const ACONS1_G: f64 = 0.0034;
const ACONS1_T: f64 = 0.0030;

// Acceptor consensus position 19 (second base of AG)
const ACONS2_A: f64 = 0.0027;
const ACONS2_C: f64 = 0.0037;
const ACONS2_G: f64 = 0.9905;
const ACONS2_T: f64 = 0.0030;

impl MaxEntScanData {
    /// Load MaxEntScan data from a LOFTEE data directory.
    /// Expects: <dir>/maxEntScan/me2x5 and <dir>/maxEntScan/splicemodels/me2x3acc{1-9}
    pub fn load(loftee_dir: &Path) -> Result<Self, String> {
        let mes_dir = loftee_dir.join("maxEntScan");

        // Load me2x5 (donor context scores)
        let me2x5_path = mes_dir.join("me2x5");
        let me2x5 = load_score_table(&me2x5_path)?;
        if me2x5.len() != 16384 {
            return Err(format!(
                "me2x5 should have 16384 entries, got {}",
                me2x5.len()
            ));
        }

        // Load 9 acceptor tables
        let models_dir = mes_dir.join("splicemodels");
        let mut metables = Vec::with_capacity(9);
        for i in 1..=9 {
            let path = models_dir.join(format!("me2x3acc{}", i));
            let table = load_score_table(&path)?;
            metables.push(table);
        }

        Ok(Self { me2x5, metables })
    }

    /// Score a 9-base donor splice site sequence.
    /// Returns the MaxEntScan donor score (log2 scale).
    pub fn score_donor(&self, seq: &[u8]) -> Option<f64> {
        if seq.len() != 9 {
            return None;
        }
        self.score_donor_impl(seq)
    }

    fn score_donor_impl(&self, seq: &[u8]) -> Option<f64> {
        let cons_score = score_consensus_donor(seq)?;
        // Get context bases [0,1,2,5,6,7,8] (skip positions 3,4)
        let rest: Vec<u8> = [seq[0], seq[1], seq[2], seq[5], seq[6], seq[7], seq[8]].to_vec();
        let idx = hashseq(&rest)?;
        if idx >= self.me2x5.len() {
            return None;
        }
        let context_score = self.me2x5[idx];
        let raw = cons_score * context_score;
        if raw <= 0.0 {
            return None;
        }
        Some(raw.log2())
    }

    /// Score a 23-base acceptor splice site sequence.
    /// Returns the MaxEntScan acceptor score (log2 scale).
    pub fn score_acceptor(&self, seq: &[u8]) -> Option<f64> {
        if seq.len() != 23 {
            return None;
        }
        self.score_acceptor_impl(seq)
    }

    fn score_acceptor_impl(&self, seq: &[u8]) -> Option<f64> {
        let cons_score = score_consensus_acceptor(seq)?;
        // Get rest: positions [0..18] + [20..23] (skip positions 18,19 = AG consensus)
        let mut rest = Vec::with_capacity(21);
        rest.extend_from_slice(&seq[0..18]);
        rest.extend_from_slice(&seq[20..23]);

        let me_score = max_ent_score(&rest, &self.metables)?;
        let raw = cons_score * me_score;
        if raw <= 0.0 {
            return None;
        }
        Some(raw.log2())
    }
}

/// Donor consensus scoring (positions 3 and 4 of the 9-base sequence).
fn score_consensus_donor(seq: &[u8]) -> Option<f64> {
    let c1 = match seq[3] {
        b'A' | b'a' => CONS1_A,
        b'C' | b'c' => CONS1_C,
        b'G' | b'g' => CONS1_G,
        b'T' | b't' => CONS1_T,
        _ => return None,
    };
    let c2 = match seq[4] {
        b'A' | b'a' => CONS2_A,
        b'C' | b'c' => CONS2_C,
        b'G' | b'g' => CONS2_G,
        b'T' | b't' => CONS2_T,
        _ => return None,
    };
    let b1 = bgd(seq[3])?;
    let b2 = bgd(seq[4])?;
    Some((c1 * c2) / (b1 * b2))
}

/// Acceptor consensus scoring (positions 18 and 19 of the 23-base sequence).
fn score_consensus_acceptor(seq: &[u8]) -> Option<f64> {
    let c1 = match seq[18] {
        b'A' | b'a' => ACONS1_A,
        b'C' | b'c' => ACONS1_C,
        b'G' | b'g' => ACONS1_G,
        b'T' | b't' => ACONS1_T,
        _ => return None,
    };
    let c2 = match seq[19] {
        b'A' | b'a' => ACONS2_A,
        b'C' | b'c' => ACONS2_C,
        b'G' | b'g' => ACONS2_G,
        b'T' | b't' => ACONS2_T,
        _ => return None,
    };
    let b1 = bgd(seq[18])?;
    let b2 = bgd(seq[19])?;
    Some((c1 * c2) / (b1 * b2))
}

/// MaxEnt score for the 21-base acceptor context (after removing AG consensus).
fn max_ent_score(seq: &[u8], metables: &[Vec<f64>]) -> Option<f64> {
    if seq.len() != 21 || metables.len() != 9 {
        return None;
    }

    // 9 overlapping windows into the 21-base sequence
    let indices = [
        (0, 7),   // table 0: positions 0-6
        (7, 14),  // table 1: positions 7-13
        (14, 21), // table 2: positions 14-20
        (4, 11),  // table 3: positions 4-10
        (11, 18), // table 4: positions 11-17
        (4, 7),   // table 5: positions 4-6 (3-mer)
        (7, 11),  // table 6: positions 7-10 (4-mer)
        (11, 14), // table 7: positions 11-13 (3-mer)
        (14, 18), // table 8: positions 14-17 (4-mer)
    ];

    let mut scores = [0.0f64; 9];
    for (i, &(start, end)) in indices.iter().enumerate() {
        let subseq = &seq[start..end];
        let idx = hashseq(subseq)?;
        if idx >= metables[i].len() {
            return None;
        }
        scores[i] = metables[i][idx];
    }

    // Final score: (sc0 * sc1 * sc2 * sc3 * sc4) / (sc5 * sc6 * sc7 * sc8)
    let numerator = scores[0] * scores[1] * scores[2] * scores[3] * scores[4];
    let denominator = scores[5] * scores[6] * scores[7] * scores[8];
    if denominator == 0.0 {
        return None;
    }
    Some(numerator / denominator)
}

/// Convert a DNA sequence to a base-4 index.
fn hashseq(seq: &[u8]) -> Option<usize> {
    let mut result = 0usize;
    for (i, &base) in seq.iter().enumerate() {
        let val = match base {
            b'A' | b'a' => 0,
            b'C' | b'c' => 1,
            b'G' | b'g' => 2,
            b'T' | b't' => 3,
            _ => return None,
        };
        result += val * 4usize.pow((seq.len() - i - 1) as u32);
    }
    Some(result)
}

fn bgd(base: u8) -> Option<f64> {
    match base {
        b'A' | b'a' => Some(BGD_A),
        b'C' | b'c' => Some(BGD_C),
        b'G' | b'g' => Some(BGD_G),
        b'T' | b't' => Some(BGD_T),
        _ => None,
    }
}

/// Load a score table from a file (one float per line).
fn load_score_table(path: &Path) -> Result<Vec<f64>, String> {
    let file = std::fs::File::open(path)
        .map_err(|e| format!("Cannot open {}: {}", path.display(), e))?;
    let reader = BufReader::new(file);
    let mut scores = Vec::new();
    for line in reader.lines() {
        let line = line.map_err(|e| format!("Read error in {}: {}", path.display(), e))?;
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        let score: f64 = trimmed
            .parse()
            .map_err(|e| format!("Parse error in {}: {} (line: '{}')", path.display(), e, trimmed))?;
        scores.push(score);
    }
    Ok(scores)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hashseq() {
        assert_eq!(hashseq(b"A"), Some(0));
        assert_eq!(hashseq(b"C"), Some(1));
        assert_eq!(hashseq(b"G"), Some(2));
        assert_eq!(hashseq(b"T"), Some(3));
        assert_eq!(hashseq(b"AA"), Some(0));
        assert_eq!(hashseq(b"AC"), Some(1));
        assert_eq!(hashseq(b"CA"), Some(4));
        assert_eq!(hashseq(b"TT"), Some(15));
        // 7-mer AAAAAAA = 0
        assert_eq!(hashseq(b"AAAAAAA"), Some(0));
        // 7-mer TTTTTTT = 4^7 - 1 = 16383
        assert_eq!(hashseq(b"TTTTTTT"), Some(16383));
        // Invalid base
        assert_eq!(hashseq(b"N"), None);
    }

    #[test]
    fn test_score_consensus_donor() {
        // GT at positions 3-4 should give high score
        let gt_seq = b"AAAGTAAAA";
        let score = score_consensus_donor(gt_seq).unwrap();
        // Expected: (0.9896 * 0.9884) / (0.23 * 0.27) ≈ 15.74
        assert!(score > 10.0);

        // Non-GT should give low score
        let aa_seq = b"AAAAAAAA\x00"; // 9 bytes
        // Use proper 9-byte slice
        let low_seq = b"AAAAAAAAA";
        let low_score = score_consensus_donor(low_seq).unwrap();
        assert!(low_score < 1.0);
    }

    #[test]
    fn test_score_consensus_acceptor() {
        // AG at positions 18-19 should give high score
        let mut seq = [b'A'; 23];
        seq[18] = b'A';
        seq[19] = b'G';
        let score = score_consensus_acceptor(&seq).unwrap();
        // Expected: (0.9903 * 0.9905) / (0.27 * 0.23) ≈ 15.79
        assert!(score > 10.0);
    }
}
