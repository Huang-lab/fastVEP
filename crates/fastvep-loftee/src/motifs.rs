//! Splice regulatory element (SRE) motif scanning.
//!
//! Scans sequences for exonic/intronic splice enhancers and silencers (ESE, ESS, ISE, ISS).

use std::collections::HashSet;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// Pre-loaded motif sets for donor or acceptor splice contexts.
#[derive(Debug, Clone)]
pub struct MotifSet {
    pub ese: HashSet<String>,
    pub ess: HashSet<String>,
    pub ise: HashSet<String>,
    pub iss: HashSet<String>,
}

impl MotifSet {
    /// Load motif files from a directory (expects ese.txt, ess.txt, ise.txt, iss.txt).
    pub fn load(dir: &Path) -> Result<Self, String> {
        Ok(Self {
            ese: load_motif_file(&dir.join("ese.txt"))?,
            ess: load_motif_file(&dir.join("ess.txt"))?,
            ise: load_motif_file(&dir.join("ise.txt"))?,
            iss: load_motif_file(&dir.join("iss.txt"))?,
        })
    }

    /// Scan a sequence for motif hits of the given type.
    pub fn scan(&self, seq: &[u8], motif_type: MotifType) -> usize {
        let motifs = match motif_type {
            MotifType::Ese => &self.ese,
            MotifType::Ess => &self.ess,
            MotifType::Ise => &self.ise,
            MotifType::Iss => &self.iss,
        };
        scan_seq(seq, motifs)
    }
}

#[derive(Debug, Clone, Copy)]
pub enum MotifType {
    Ese,
    Ess,
    Ise,
    Iss,
}

/// Scan a sequence for k-mer hits in a motif set (k = 6, 7, 8).
fn scan_seq(seq: &[u8], motifs: &HashSet<String>) -> usize {
    let seq_upper: Vec<u8> = seq.iter().map(|&b| b.to_ascii_uppercase()).collect();
    let seq_str = String::from_utf8_lossy(&seq_upper);
    let len = seq_str.len();
    let mut hits = 0;

    for i in 0..len {
        for k in [6, 7, 8] {
            if i + k <= len {
                let kmer = &seq_str[i..i + k];
                if motifs.contains(kmer) {
                    hits += 1;
                }
            }
        }
    }
    hits
}

fn load_motif_file(path: &Path) -> Result<HashSet<String>, String> {
    let file = std::fs::File::open(path)
        .map_err(|e| format!("Cannot open {}: {}", path.display(), e))?;
    let reader = BufReader::new(file);
    let mut motifs = HashSet::new();
    for line in reader.lines() {
        let line = line.map_err(|e| format!("Read error in {}: {}", path.display(), e))?;
        let trimmed = line.trim().to_uppercase();
        if !trimmed.is_empty() {
            motifs.insert(trimmed);
        }
    }
    Ok(motifs)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_scan_seq_basic() {
        let mut motifs = HashSet::new();
        motifs.insert("AAAAAA".to_string()); // 6-mer
        // Sequence with two overlapping hits
        assert_eq!(scan_seq(b"AAAAAAA", &motifs), 2); // positions 0-5 and 1-6
    }

    #[test]
    fn test_scan_seq_no_hits() {
        let mut motifs = HashSet::new();
        motifs.insert("GGGGGG".to_string());
        assert_eq!(scan_seq(b"AAAAAA", &motifs), 0);
    }

    #[test]
    fn test_scan_seq_case_insensitive() {
        let mut motifs = HashSet::new();
        motifs.insert("AAAAAA".to_string());
        assert_eq!(scan_seq(b"aaaaaa", &motifs), 1);
    }
}
