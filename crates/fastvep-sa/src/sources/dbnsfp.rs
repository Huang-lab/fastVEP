//! dbNSFP parser for building .osa annotation files.
//!
//! dbNSFP provides pre-computed functional predictions (SIFT, PolyPhen,
//! REVEL, CADD, etc.) for all possible missense variants.
//!
//! This parser extracts SIFT, PolyPhen, AlphaMissense, and BayesDel predictions.

use crate::common::AnnotationRecord;
use anyhow::{Context, Result};
use std::collections::HashMap;
use std::io::BufRead;

/// Parse a dbNSFP TSV file to extract SIFT and PolyPhen predictions.
///
/// Expected header columns include:
/// `#chr`, `pos(1-based)`, `ref`, `alt`, `SIFT_score`, `SIFT_pred`,
/// `Polyphen2_HDIV_score`, `Polyphen2_HDIV_pred`
///
/// Column indices are auto-detected from the header row.
pub fn parse_dbnsfp<R: BufRead>(
    reader: R,
    chrom_to_idx: &HashMap<String, u16>,
) -> Result<Vec<AnnotationRecord>> {
    let mut records = Vec::new();
    let mut col_indices: Option<DbNsfpColumns> = None;

    for line in reader.lines() {
        let line = line.context("Reading dbNSFP line")?;

        if line.starts_with('#') || line.starts_with("chr\t") {
            // Parse header to find column indices
            let header = line.trim_start_matches('#');
            col_indices = Some(DbNsfpColumns::from_header(header)?);
            continue;
        }

        if line.is_empty() {
            continue;
        }

        let cols = match &col_indices {
            Some(c) => c,
            None => continue,
        };

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() <= cols.max_idx() {
            continue;
        }

        let chrom = normalize_chrom(fields[cols.chr]);
        let chrom_idx = match chrom_to_idx.get(&chrom) {
            Some(&idx) => idx,
            None => continue,
        };

        let pos: u32 = match fields[cols.pos].parse() {
            Ok(p) => p,
            Err(_) => continue,
        };

        let ref_allele = fields[cols.ref_col].to_string();
        let alt_allele = fields[cols.alt].to_string();

        let mut parts = Vec::new();

        // SIFT
        if let Some(idx) = cols.sift_score {
            let score_str = fields[idx];
            if score_str != "." {
                // May have multiple scores separated by ";" — take the first
                let score = score_str.split(';').next().unwrap_or(".");
                if let Ok(s) = score.parse::<f64>() {
                    let pred = cols.sift_pred.and_then(|i| {
                        let p = fields[i].split(';').next()?;
                        if p == "." { None } else { Some(p) }
                    });
                    let pred_str = match pred {
                        Some("D") => "deleterious",
                        Some("T") => "tolerated",
                        _ => "",
                    };
                    if !pred_str.is_empty() {
                        parts.push(format!("\"sift\":\"{}({:.3})\"", pred_str, s));
                    } else {
                        parts.push(format!("\"sift\":\"{:.3}\"", s));
                    }
                }
            }
        }

        // PolyPhen
        if let Some(idx) = cols.polyphen_score {
            let score_str = fields[idx];
            if score_str != "." {
                let score = score_str.split(';').next().unwrap_or(".");
                if let Ok(s) = score.parse::<f64>() {
                    let pred = cols.polyphen_pred.and_then(|i| {
                        let p = fields[i].split(';').next()?;
                        if p == "." { None } else { Some(p) }
                    });
                    let pred_str = match pred {
                        Some("D") => "probably_damaging",
                        Some("P") => "possibly_damaging",
                        Some("B") => "benign",
                        _ => "",
                    };
                    if !pred_str.is_empty() {
                        parts.push(format!("\"polyphen\":\"{}({:.3})\"", pred_str, s));
                    } else {
                        parts.push(format!("\"polyphen\":\"{:.3}\"", s));
                    }
                }
            }
        }

        // AlphaMissense — calibrated pathogenicity score (0-1). Emitted numeric,
        // like REVEL, so the classifier can threshold it directly.
        if let Some(idx) = cols.alphamissense_score {
            if let Some(v) = first_score(fields[idx]) {
                parts.push(format!("\"alphaMissense\":{}", v));
            }
        }

        // BayesDel (no-AF) — frequency-independent deleteriousness score.
        if let Some(idx) = cols.bayesdel_score {
            if let Some(v) = first_score(fields[idx]) {
                parts.push(format!("\"bayesDel\":{}", v));
            }
        }

        if parts.is_empty() {
            continue;
        }

        records.push(AnnotationRecord {
            chrom_idx,
            position: pos,
            ref_allele,
            alt_allele,
            json: format!("{{{}}}", parts.join(",")),
        });
    }

    records.sort_by(|a, b| a.chrom_idx.cmp(&b.chrom_idx).then(a.position.cmp(&b.position)));
    Ok(records)
}

struct DbNsfpColumns {
    chr: usize,
    pos: usize,
    ref_col: usize,
    alt: usize,
    sift_score: Option<usize>,
    sift_pred: Option<usize>,
    polyphen_score: Option<usize>,
    polyphen_pred: Option<usize>,
    alphamissense_score: Option<usize>,
    bayesdel_score: Option<usize>,
}

impl DbNsfpColumns {
    fn from_header(header: &str) -> Result<Self> {
        let fields: Vec<&str> = header.split('\t').collect();
        let find = |names: &[&str]| -> Option<usize> {
            fields.iter().position(|f| {
                let fl = f.to_lowercase();
                names.iter().any(|n| fl == *n)
            })
        };

        Ok(Self {
            chr: find(&["chr", "#chr"]).unwrap_or(0),
            pos: find(&["pos(1-based)", "pos", "hg38_pos"]).unwrap_or(1),
            ref_col: find(&["ref", "ref_allele"]).unwrap_or(2),
            alt: find(&["alt", "alt_allele"]).unwrap_or(3),
            sift_score: find(&["sift_score"]),
            sift_pred: find(&["sift_pred"]),
            polyphen_score: find(&["polyphen2_hdiv_score"]),
            polyphen_pred: find(&["polyphen2_hdiv_pred"]),
            // Calibrated missense predictors: AlphaMissense and BayesDel. The no-AF
            // BayesDel variant is frequency-independent — the ACMG-appropriate choice
            // for PP3/BP4 (avoids circularity with the frequency criteria).
            alphamissense_score: find(&["alphamissense_score"]),
            bayesdel_score: find(&["bayesdel_noaf_score"]),
        })
    }

    fn max_idx(&self) -> usize {
        let mut m = self.alt;
        for opt in [
            self.sift_score,
            self.sift_pred,
            self.polyphen_score,
            self.polyphen_pred,
            self.alphamissense_score,
            self.bayesdel_score,
        ] {
            if let Some(i) = opt {
                m = m.max(i);
            }
        }
        m
    }
}

/// Take the first non-missing (";"-separated) value from a dbNSFP cell and parse
/// it as a float. dbNSFP uses "." for missing and may list one value per transcript.
fn first_score(cell: &str) -> Option<f64> {
    cell.split(';')
        .find(|v| *v != "." && !v.is_empty())?
        .parse::<f64>()
        .ok()
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

    #[test]
    fn test_parse_dbnsfp() {
        let data = "\
#chr\tpos(1-based)\tref\talt\tSIFT_score\tSIFT_pred\tPolyphen2_HDIV_score\tPolyphen2_HDIV_pred\tAlphaMissense_score\tBayesDel_noAF_score
1\t10001\tA\tG\t0.032\tD\t0.998\tD\t0.9876\t0.512
1\t10001\tA\tC\t0.450\tT\t0.100\tB\t.\t-0.231
1\t10002\tC\tT\t.\t.\t.\t.\t.\t.
";
        let mut chrom_map = HashMap::new();
        chrom_map.insert("chr1".into(), 0u16);

        let records = parse_dbnsfp(data.as_bytes(), &chrom_map).unwrap();
        // Third line has all dots, should be skipped
        assert_eq!(records.len(), 2);

        assert!(records[0].json.contains("deleterious(0.032)"));
        assert!(records[0].json.contains("probably_damaging(0.998)"));

        assert!(records[0].json.contains("\"alphaMissense\":0.9876"));
        assert!(records[0].json.contains("\"bayesDel\":0.512"));

        assert!(records[1].json.contains("tolerated(0.450)"));
        assert!(records[1].json.contains("benign(0.100)"));
        // AlphaMissense missing (".") -> key omitted; BayesDel present (negative).
        assert!(!records[1].json.contains("alphaMissense"));
        assert!(records[1].json.contains("\"bayesDel\":-0.231"));
    }
}
