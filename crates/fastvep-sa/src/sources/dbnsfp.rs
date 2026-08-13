//! dbNSFP parser for building .osa annotation files.
//!
//! dbNSFP provides pre-computed functional predictions (SIFT, PolyPhen,
//! REVEL, CADD, etc.) for all possible missense variants.
//!
//! This parser extracts SIFT and PolyPhen predictions specifically.

use crate::common::AnnotationRecord;
use crate::writer_v2::Osa2Metadata;
use anyhow::{Context, Result};
use std::collections::HashMap;
use std::io::BufRead;

/// Standard dbNSFP `.osa2` metadata. dbNSFP's payload is composite prediction
/// strings (`{"sift":"D(0.012)","polyphen":..}`) that don't decompose into
/// numeric u32 columns, so it is stored as a whole-record JSON blob (see
/// [`crate::writer_v2::raw_json_blob_fields`]): byte-identical to v1, with v2's
/// chunk-level zstd shrinking the database.
pub fn dbnsfp_osa2_metadata(assembly: &str) -> Osa2Metadata {
    Osa2Metadata {
        format_version: 2,
        name: "dbNSFP".into(),
        version: "latest".into(),
        assembly: assembly.into(),
        json_key: "dbnsfp".into(),
        match_by_allele: true,
        is_array: false,
        is_positional: false,
        chunk_bits: 20,
        description: format!("dbNSFP SIFT/PolyPhen predictions for {assembly}"),
    }
}

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
                        if p == "." {
                            None
                        } else {
                            Some(p)
                        }
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
                        if p == "." {
                            None
                        } else {
                            Some(p)
                        }
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

    records.sort_by(|a, b| {
        a.chrom_idx
            .cmp(&b.chrom_idx)
            .then(a.position.cmp(&b.position))
    });
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
        })
    }

    fn max_idx(&self) -> usize {
        let mut m = self.alt;
        for i in [
            self.sift_score,
            self.sift_pred,
            self.polyphen_score,
            self.polyphen_pred,
        ]
        .into_iter()
        .flatten()
        {
            m = m.max(i);
        }
        m
    }
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
#chr\tpos(1-based)\tref\talt\tSIFT_score\tSIFT_pred\tPolyphen2_HDIV_score\tPolyphen2_HDIV_pred
1\t10001\tA\tG\t0.032\tD\t0.998\tD
1\t10001\tA\tC\t0.450\tT\t0.100\tB
1\t10002\tC\tT\t.\t.\t.\t.
";
        let mut chrom_map = HashMap::new();
        chrom_map.insert("chr1".into(), 0u16);

        let records = parse_dbnsfp(data.as_bytes(), &chrom_map).unwrap();
        // Third line has all dots, should be skipped
        assert_eq!(records.len(), 2);

        assert!(records[0].json.contains("deleterious(0.032)"));
        assert!(records[0].json.contains("probably_damaging(0.998)"));

        assert!(records[1].json.contains("tolerated(0.450)"));
        assert!(records[1].json.contains("benign(0.100)"));
    }
}
