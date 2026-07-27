//! TOPMed population frequency parser.
//!
//! The full TOPMed freeze VCF contains ~450 million records. Use
//! `iter_topmed_vcf` for streaming builds; `parse_topmed_vcf` is retained for
//! tests and small inputs.

use crate::common::AnnotationRecord;
use crate::fields::{Field, FieldType};
use crate::writer_v2::Osa2Metadata;
use anyhow::{Context, Result};
use std::collections::HashMap;
use std::collections::VecDeque;
use std::io::BufRead;

/// Canonical TOPMed `.osa2` field schema. Aliases match the JSON keys the v1
/// builder emits (`allAf`, `allAc`, `allAn`) and AF renders in the same
/// scientific notation. AC/AN are stored exactly; AF is quantized into a u32
/// column at a fixed absolute resolution (`1 / 2_000_000`), so v1 and v2 agree
/// for every TOPMed frequency (all above that floor) but AF is not guaranteed
/// byte-identical for arbitrarily fine values. See gnomAD's `AF_MULTIPLIER`.
pub fn topmed_osa2_fields() -> Vec<Field> {
    vec![
        Field {
            field: "allAf".into(), alias: "allAf".into(), ftype: FieldType::Float,
            multiplier: 2_000_000, zigzag: false, missing_value: u32::MAX,
            missing_string: ".".into(), description: "Global allele frequency".into(),
        },
        Field {
            field: "allAc".into(), alias: "allAc".into(), ftype: FieldType::Integer,
            multiplier: 1, zigzag: false, missing_value: u32::MAX,
            missing_string: ".".into(), description: "Allele count".into(),
        },
        Field {
            field: "allAn".into(), alias: "allAn".into(), ftype: FieldType::Integer,
            multiplier: 1, zigzag: false, missing_value: u32::MAX,
            missing_string: ".".into(), description: "Total allele number".into(),
        },
    ]
}

/// Standard TOPMed `.osa2` metadata (mirrors the v1 header:
/// `json_key = "topmed"`, allele-matched, non-positional).
pub fn topmed_osa2_metadata(assembly: &str) -> Osa2Metadata {
    Osa2Metadata {
        format_version: 2,
        name: "TOPMed".into(),
        version: "latest".into(),
        assembly: assembly.into(),
        json_key: "topmed".into(),
        match_by_allele: true,
        is_array: false,
        is_positional: false,
        chunk_bits: 20,
        description: format!("TOPMed population frequencies for {assembly}"),
    }
}

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
                // AC/AN are written unquoted, so each must be a validated
                // non-negative integer — raw INFO-field text would otherwise
                // land in the JSON as a bare, unquoted token and break every
                // downstream serde_json::from_str on this record. Mirrors
                // the CNT validation in sources/cosmic.rs.
                if let Some(ac) = all_acs.get(i).and_then(|s| s.parse::<u64>().ok()) {
                    parts.push(format!("\"allAc\":{}", ac));
                }
                if let Some(an) = info_map.get("AN").and_then(|s| s.parse::<u64>().ok()) {
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

    #[test]
    fn test_garbage_ac_an_omitted_not_emitted_unescaped() {
        // AC/AN are spliced into the JSON unquoted, so malformed upstream
        // text must be dropped rather than emitted as a bare token that
        // would break serde_json::from_str on the whole record.
        let vcf = "#h\nchr1\t100\t.\tA\tG\t.\t.\tAF=0.05;AC=garbage;AN=not_a_number\n";
        let mut m = HashMap::new();
        m.insert("chr1".into(), 0u16);
        let recs = parse_topmed_vcf(vcf.as_bytes(), &m).unwrap();
        assert_eq!(recs.len(), 1);
        assert!(!recs[0].json.contains("allAc"));
        assert!(!recs[0].json.contains("allAn"));
        let v: serde_json::Value = serde_json::from_str(&recs[0].json).unwrap();
        assert!(v.get("allAf").is_some());
    }

    #[test]
    fn test_topmed_osa2_bridge_encodes_values() {
        let vcf = "#h\nchr1\t100\t.\tA\tG\t.\t.\tAF=0.05;AC=500;AN=10000\n";
        let mut m = HashMap::new();
        m.insert("chr1".into(), 0u16);
        let recs = parse_topmed_vcf(vcf.as_bytes(), &m).unwrap();
        let fields = topmed_osa2_fields();
        let o = crate::writer_v2::osa2_record_from_v1(&recs[0], "chr1".into(), &fields).unwrap();
        assert_eq!(o.values[0], fields[0].encode_float(0.05)); // allAf
        assert_eq!(o.values[1], 500); // allAc
        assert_eq!(o.values[2], 10000); // allAn
    }
}
