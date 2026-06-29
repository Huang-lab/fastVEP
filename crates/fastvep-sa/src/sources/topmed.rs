//! TOPMed population frequency parser.
//!
//! The full TOPMed freeze VCF contains ~450 million records. Use
//! `iter_topmed_vcf` for streaming builds; `parse_topmed_vcf` is retained for
//! tests and small inputs.

use crate::common::AnnotationRecord;
use anyhow::{Context, Result};
use std::collections::HashMap;
use std::collections::VecDeque;
use std::io::BufRead;

/// Stream a coordinate-sorted TOPMed VCF as `AnnotationRecord`s without
/// buffering the whole file in memory.
///
/// The input must already be sorted by chromosome and position (all standard
/// TOPMed freeze releases are).
pub fn iter_topmed_vcf<'a, R: BufRead>(
    reader: R,
    chrom_to_idx: &'a HashMap<String, u16>,
) -> TopmedRecordIter<'a, R> {
    TopmedRecordIter {
        lines: reader.lines(),
        chrom_to_idx,
        pending: VecDeque::new(),
    }
}

pub struct TopmedRecordIter<'a, R: BufRead> {
    lines: std::io::Lines<R>,
    chrom_to_idx: &'a HashMap<String, u16>,
    pending: VecDeque<AnnotationRecord>,
}

impl<R: BufRead> Iterator for TopmedRecordIter<'_, R> {
    type Item = Result<AnnotationRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            if let Some(record) = self.pending.pop_front() {
                return Some(Ok(record));
            }

            let line = match self.lines.next()? {
                Ok(l) => l,
                Err(e) => return Some(Err(e).context("Reading TOPMed VCF line")),
            };

            if line.starts_with('#') {
                continue;
            }

            let fields: Vec<&str> = line.splitn(9, '\t').collect();
            if fields.len() < 8 {
                continue;
            }

            let chrom = normalize_chrom(fields[0]);
            let chrom_idx = match self.chrom_to_idx.get(&chrom) {
                Some(&i) => i,
                None => continue,
            };
            let pos: u32 = match fields[1].parse() {
                Ok(p) => p,
                Err(_) => continue,
            };
            let ref_allele = fields[3].to_string();
            let alt_field = fields[4];
            let info = fields[7];
            let info_map = parse_info(info);

            let alts: Vec<&str> = alt_field.split(',').collect();
            let all_afs = split_vals(info_map.get("AF").map(|s| s.as_str()));
            let all_acs = split_vals(info_map.get("AC").map(|s| s.as_str()));

            for (i, alt) in alts.iter().enumerate() {
                if *alt == "." || *alt == "*" {
                    continue;
                }
                let mut parts = Vec::new();
                if let Some(af) = all_afs.get(i).and_then(|s| s.parse::<f64>().ok()) {
                    parts.push(format!("\"allAf\":{:.6e}", af));
                }
                if let Some(ac) = all_acs.get(i) {
                    parts.push(format!("\"allAc\":{}", ac));
                }
                if let Some(an) = info_map.get("AN") {
                    parts.push(format!("\"allAn\":{}", an));
                }
                if parts.is_empty() {
                    continue;
                }
                self.pending.push_back(AnnotationRecord {
                    chrom_idx,
                    position: pos,
                    ref_allele: ref_allele.clone(),
                    alt_allele: alt.to_string(),
                    json: format!("{{{}}}", parts.join(",")),
                });
            }
        }
    }
}

/// Parse a TOPMed freeze VCF into sorted `AnnotationRecord`s.
///
/// Loads all records into memory — suitable for tests and small inputs.
/// For the full TOPMed release use `iter_topmed_vcf` via the pipeline instead.
pub fn parse_topmed_vcf<R: BufRead>(
    reader: R,
    chrom_to_idx: &HashMap<String, u16>,
) -> Result<Vec<AnnotationRecord>> {
    let mut records: Vec<_> = iter_topmed_vcf(reader, chrom_to_idx).collect::<Result<_>>()?;
    records.sort_by(|a, b| a.chrom_idx.cmp(&b.chrom_idx).then(a.position.cmp(&b.position)));
    Ok(records)
}

fn parse_info(info: &str) -> HashMap<String, String> {
    let mut m = HashMap::new();
    for p in info.split(';') {
        if let Some((k, v)) = p.split_once('=') {
            m.insert(k.into(), v.into());
        }
    }
    m
}

fn split_vals(v: Option<&str>) -> Vec<String> {
    v.map(|s| s.split(',').map(|x| x.to_string()).collect())
        .unwrap_or_default()
}

fn normalize_chrom(c: &str) -> String {
    if c.starts_with("chr") {
        c.to_string()
    } else {
        format!("chr{}", c)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_topmed() {
        let vcf = "#h\nchr1\t100\t.\tA\tG\t.\t.\tAF=0.05;AC=500;AN=10000\n";
        let mut m = HashMap::new();
        m.insert("chr1".into(), 0u16);
        let recs = parse_topmed_vcf(vcf.as_bytes(), &m).unwrap();
        assert_eq!(recs.len(), 1);
        assert!(recs[0].json.contains("\"allAf\":"));
    }
}
