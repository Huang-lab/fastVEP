//! Generic position-level score parsers (PhyloP, GERP, DANN, etc.).
//!
//! These parsers produce positional AnnotationRecords where the JSON is
//! just the numeric score value as a string. The SaWriter stores them
//! as positional annotations (match_by_allele=false, is_positional=true).
//!
//! PhyloP/GERP are per-base, genome-wide sources (~3 billion positions for
//! hg38) — denser than any VCF-based source. Use `iter_score_tsv` /
//! `iter_wigfix` for streaming builds; `parse_score_tsv` / `parse_wigfix`
//! buffer everything in memory and are retained for tests and small inputs.

use crate::common::AnnotationRecord;
use anyhow::{Context, Result};
use std::collections::HashMap;
use std::io::BufRead;

/// Stream a tab-separated score file (chrom, pos, score) as `AnnotationRecord`s
/// without buffering the whole file in memory.
///
/// Supports formats like:
/// - BED-like: `chr1\t12345\t12346\t2.345`  (4 columns: chrom, start, end, score)
/// - Simple:   `chr1\t12345\t2.345`          (3 columns: chrom, pos, score)
///
/// Positions are 0-based in BED format (converted to 1-based internally)
/// or 1-based in simple format. Set `zero_based` accordingly. The input must
/// already be sorted by chromosome and position.
pub fn iter_score_tsv<R: BufRead>(
    reader: R,
    chrom_to_idx: &HashMap<String, u16>,
    zero_based: bool,
) -> ScoreTsvIter<'_, R> {
    ScoreTsvIter { lines: reader.lines(), chrom_to_idx, zero_based }
}

pub struct ScoreTsvIter<'a, R: BufRead> {
    lines: std::io::Lines<R>,
    chrom_to_idx: &'a HashMap<String, u16>,
    zero_based: bool,
}

impl<R: BufRead> Iterator for ScoreTsvIter<'_, R> {
    type Item = Result<AnnotationRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            let line = match self.lines.next()? {
                Ok(l) => l,
                Err(e) => return Some(Err(e).context("Reading score TSV line")),
            };
            if line.starts_with('#') || line.starts_with("track") || line.is_empty() {
                continue;
            }

            let fields: Vec<&str> = line.split('\t').collect();

            let (chrom_str, pos_str, score_str) = match fields.len() {
                // BED 4-column: chrom, start (0-based), end, score
                4.. => (fields[0], fields[1], fields[3]),
                // Simple 3-column: chrom, pos, score
                3 => (fields[0], fields[1], fields[2]),
                _ => continue,
            };

            let chrom = normalize_chrom(chrom_str);
            let chrom_idx = match self.chrom_to_idx.get(&chrom) {
                Some(&idx) => idx,
                None => continue,
            };

            let pos: u32 = match pos_str.parse::<u32>() {
                Ok(p) => {
                    if self.zero_based { p + 1 } else { p }
                }
                Err(_) => continue,
            };

            let score: f64 = match score_str.trim().parse() {
                Ok(s) => s,
                Err(_) => continue,
            };

            return Some(Ok(AnnotationRecord {
                chrom_idx,
                position: pos,
                ref_allele: String::new(),
                alt_allele: String::new(),
                json: format_score(score),
            }));
        }
    }
}

/// Parse a tab-separated score file into sorted `AnnotationRecord`s.
///
/// Loads all records into memory — suitable for tests and small inputs.
/// For genome-wide sources use `iter_score_tsv` via the pipeline instead.
pub fn parse_score_tsv<R: BufRead>(
    reader: R,
    chrom_to_idx: &HashMap<String, u16>,
    zero_based: bool,
) -> Result<Vec<AnnotationRecord>> {
    let mut records: Vec<_> =
        iter_score_tsv(reader, chrom_to_idx, zero_based).collect::<Result<_>>()?;
    records.sort_by(|a, b| a.chrom_idx.cmp(&b.chrom_idx).then(a.position.cmp(&b.position)));
    Ok(records)
}

/// Stream a UCSC wiggle fixed-step (wigFix) file as `AnnotationRecord`s
/// without buffering the whole file in memory.
///
/// Format:
/// ```text
/// fixedStep chrom=chr1 start=10001 step=1
/// 0.064
/// -0.002
/// ...
/// ```
///
/// The input must already be sorted by chromosome (all standard UCSC
/// per-chromosome releases are; concatenating them in chromosome order
/// preserves this).
pub fn iter_wigfix<R: BufRead>(reader: R, chrom_to_idx: &HashMap<String, u16>) -> WigFixIter<'_, R> {
    WigFixIter { lines: reader.lines(), chrom_to_idx, current_chrom_idx: None, current_pos: 0, step: 1 }
}

pub struct WigFixIter<'a, R: BufRead> {
    lines: std::io::Lines<R>,
    chrom_to_idx: &'a HashMap<String, u16>,
    current_chrom_idx: Option<u16>,
    current_pos: u32,
    step: u32,
}

impl<R: BufRead> Iterator for WigFixIter<'_, R> {
    type Item = Result<AnnotationRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            let line = match self.lines.next()? {
                Ok(l) => l,
                Err(e) => return Some(Err(e).context("Reading wigFix line")),
            };

            if line.starts_with("fixedStep") {
                // Parse header: fixedStep chrom=chr1 start=10001 step=1
                let mut chrom = None;
                let mut start = None;
                self.step = 1;

                for part in line.split_whitespace().skip(1) {
                    if let Some((key, val)) = part.split_once('=') {
                        match key {
                            "chrom" => chrom = Some(normalize_chrom(val)),
                            "start" => start = val.parse().ok(),
                            "step" => self.step = val.parse().unwrap_or(1),
                            _ => {}
                        }
                    }
                }

                self.current_chrom_idx = chrom.as_ref().and_then(|c| self.chrom_to_idx.get(c)).copied();
                self.current_pos = start.unwrap_or(1);
                continue;
            }

            if line.starts_with('#') || line.starts_with("track") || line.is_empty() {
                continue;
            }

            if let Some(chrom_idx) = self.current_chrom_idx {
                let pos = self.current_pos;
                self.current_pos += self.step;
                if let Ok(score) = line.trim().parse::<f64>() {
                    return Some(Ok(AnnotationRecord {
                        chrom_idx,
                        position: pos,
                        ref_allele: String::new(),
                        alt_allele: String::new(),
                        json: format_score(score),
                    }));
                }
                continue;
            }
        }
    }
}

/// Parse a UCSC wiggle fixed-step (wigFix) file into sorted AnnotationRecords.
///
/// Loads all records into memory — suitable for tests and small inputs.
/// For genome-wide builds use `iter_wigfix` via the pipeline instead.
pub fn parse_wigfix<R: BufRead>(
    reader: R,
    chrom_to_idx: &HashMap<String, u16>,
) -> Result<Vec<AnnotationRecord>> {
    let mut records: Vec<_> = iter_wigfix(reader, chrom_to_idx).collect::<Result<_>>()?;
    records.sort_by(|a, b| a.chrom_idx.cmp(&b.chrom_idx).then(a.position.cmp(&b.position)));
    Ok(records)
}

/// Format a score compactly: drop trailing zeros, max 4 decimal places.
fn format_score(score: f64) -> String {
    if score == 0.0 {
        return "0".into();
    }
    // Use up to 4 decimal places, strip trailing zeros
    let s = format!("{:.4}", score);
    let s = s.trim_end_matches('0');
    let s = s.trim_end_matches('.');
    s.to_string()
}

fn normalize_chrom(chrom: &str) -> String {
    if chrom.starts_with("chr") {
        chrom.to_string()
    } else {
        format!("chr{}", chrom)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn test_chrom_map() -> HashMap<String, u16> {
        let mut m = HashMap::new();
        m.insert("chr1".into(), 0);
        m.insert("chr2".into(), 1);
        m
    }

    #[test]
    fn test_parse_score_tsv_bed4() {
        let data = "\
chr1\t99\t100\t2.345
chr1\t199\t200\t-1.5
chr2\t49\t50\t0.001
";
        let records = parse_score_tsv(data.as_bytes(), &test_chrom_map(), true).unwrap();
        assert_eq!(records.len(), 3);
        assert_eq!(records[0].position, 100); // 0-based 99 -> 1-based 100
        assert_eq!(records[0].json, "2.345");
        assert_eq!(records[1].json, "-1.5");
        assert_eq!(records[2].chrom_idx, 1);
        assert_eq!(records[2].json, "0.001");
    }

    #[test]
    fn test_parse_score_tsv_simple() {
        let data = "chr1\t100\t3.14\nchr1\t200\t-0.5\n";
        let records = parse_score_tsv(data.as_bytes(), &test_chrom_map(), false).unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].position, 100);
        assert_eq!(records[0].json, "3.14");
    }

    #[test]
    fn test_parse_wigfix() {
        let data = "\
fixedStep chrom=chr1 start=10001 step=1
0.064
-0.002
1.500
fixedStep chrom=chr2 start=5000 step=1
-3.210
";
        let records = parse_wigfix(data.as_bytes(), &test_chrom_map()).unwrap();
        assert_eq!(records.len(), 4);
        // chr1 records
        assert_eq!(records[0].chrom_idx, 0);
        assert_eq!(records[0].position, 10001);
        assert_eq!(records[0].json, "0.064");
        assert_eq!(records[1].position, 10002);
        assert_eq!(records[1].json, "-0.002");
        assert_eq!(records[2].position, 10003);
        assert_eq!(records[2].json, "1.5");
        // chr2 record
        assert_eq!(records[3].chrom_idx, 1);
        assert_eq!(records[3].position, 5000);
        assert_eq!(records[3].json, "-3.21");
    }

    #[test]
    fn test_format_score() {
        assert_eq!(format_score(2.3450), "2.345");
        assert_eq!(format_score(0.0), "0");
        assert_eq!(format_score(-1.0), "-1");
        assert_eq!(format_score(0.12345), "0.1235"); // rounds to 4dp
        assert_eq!(format_score(3.0001), "3.0001");
    }
}
