//! dbSNP VCF parser for building .osa annotation files.
//!
//! Parses dbSNP's VCF release to extract RS IDs and global MAF.
//!
//! The full NCBI dbSNP VCF contains ~800 million records. Use
//! `iter_dbsnp_vcf` for streaming builds; `parse_dbsnp_vcf` is retained for
//! tests and small inputs.

use crate::common::AnnotationRecord;
use anyhow::{Context, Result};
use std::collections::HashMap;
use std::collections::VecDeque;
use std::io::BufRead;

/// Stream a coordinate-sorted dbSNP VCF as `AnnotationRecord`s without
/// buffering the whole file in memory.
///
/// The input must already be sorted by chromosome and position (all standard
/// NCBI dbSNP releases are).
pub fn iter_dbsnp_vcf<'a, R: BufRead>(
    reader: R,
    chrom_to_idx: &'a HashMap<String, u16>,
) -> DbsnpRecordIter<'a, R> {
    DbsnpRecordIter {
        lines: reader.lines(),
        chrom_to_idx,
        pending: VecDeque::new(),
    }
}

pub struct DbsnpRecordIter<'a, R: BufRead> {
    lines: std::io::Lines<R>,
    chrom_to_idx: &'a HashMap<String, u16>,
    pending: VecDeque<AnnotationRecord>,
}

impl<R: BufRead> Iterator for DbsnpRecordIter<'_, R> {
    type Item = Result<AnnotationRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            if let Some(record) = self.pending.pop_front() {
                return Some(Ok(record));
            }

            let line = match self.lines.next()? {
                Ok(l) => l,
                Err(e) => return Some(Err(e).context("Reading dbSNP VCF line")),
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
                Some(&idx) => idx,
                None => continue,
            };

            let pos: u32 = match fields[1].parse() {
                Ok(p) => p,
                Err(_) => continue,
            };

            let id = fields[2];
            let ref_allele = fields[3].to_string();
            let alt_field = fields[4];
            let info = fields[7];

            let rs_id = if id.starts_with("rs") {
                id.to_string()
            } else {
                let info_map = parse_info(info);
                match info_map.get("RS") {
                    Some(rs) => format!("rs{}", rs),
                    None => continue,
                }
            };

            let info_map = parse_info(info);
            // dbSNP's CAF is `ref_freq,alt1_freq,alt2_freq,...` in ALT order
            // (i.e. index i+1 is the frequency for the i-th ALT); index it
            // per-alt below rather than once, or every ALT past the first in
            // a multi-allelic record gets the first ALT's frequency.
            let caf_parts: Option<Vec<&str>> =
                info_map.get("CAF").map(|caf| caf.split(',').collect());

            for (i, alt) in alt_field.split(',').enumerate() {
                if alt == "." || alt == "*" {
                    continue;
                }
                let freq = caf_parts
                    .as_ref()
                    .and_then(|parts| parts.get(i + 1))
                    .and_then(|s| s.parse::<f64>().ok());
                let mut parts = vec![format!("\"id\":\"{}\"", rs_id)];
                if let Some(f) = freq {
                    parts.push(format!("\"globalMaf\":{:.6e}", f));
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

/// Parse a dbSNP VCF and produce sorted `AnnotationRecord`s.
///
/// Loads all records into memory — suitable for tests and small inputs.
/// For the full NCBI release use `iter_dbsnp_vcf` via the pipeline instead.
pub fn parse_dbsnp_vcf<R: BufRead>(
    reader: R,
    chrom_to_idx: &HashMap<String, u16>,
) -> Result<Vec<AnnotationRecord>> {
    let mut records: Vec<_> = iter_dbsnp_vcf(reader, chrom_to_idx).collect::<Result<_>>()?;
    records.sort_by(|a, b| a.chrom_idx.cmp(&b.chrom_idx).then(a.position.cmp(&b.position)));
    Ok(records)
}

fn parse_info(info: &str) -> HashMap<String, String> {
    let mut map = HashMap::new();
    for pair in info.split(';') {
        if let Some((key, value)) = pair.split_once('=') {
            map.insert(key.to_string(), value.to_string());
        }
    }
    map
}

fn normalize_chrom(chrom: &str) -> String {
    // NCBI's dbSNP VCF names contigs by RefSeq accession (`NC_000001.11`).
    // Leave those untouched so the lookup hits the accession keys the builder
    // seeds into the chromosome map; mangling them to `chrNC_000001.11` is what
    // produced "0 records parsed" (issue #51).
    if chrom.starts_with("chr") || fastvep_core::looks_like_refseq_accession(chrom) {
        chrom.to_string()
    } else {
        format!("chr{}", chrom)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_dbsnp_vcf() {
        let vcf = "\
##fileformat=VCFv4.0
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
1\t10019\trs775809821\tTA\tT\t.\t.\tRS=775809821;CAF=0.9998,0.0002
1\t10039\trs978760828\tA\tC\t.\t.\tRS=978760828
";

        let mut chrom_map = HashMap::new();
        chrom_map.insert("chr1".to_string(), 0u16);

        let records = parse_dbsnp_vcf(vcf.as_bytes(), &chrom_map).unwrap();
        assert_eq!(records.len(), 2);

        assert_eq!(records[0].position, 10019);
        assert!(records[0].json.contains("rs775809821"));
        assert!(records[0].json.contains("globalMaf"));

        assert_eq!(records[1].position, 10039);
        assert!(records[1].json.contains("rs978760828"));
        assert!(!records[1].json.contains("globalMaf")); // No CAF
    }

    #[test]
    fn test_parse_dbsnp_multiallelic_caf_indexed_per_alt() {
        // CAF=ref_freq,alt1_freq,alt2_freq — each ALT must get its own
        // frequency, not the first ALT's frequency repeated for every ALT.
        let vcf = "\
##fileformat=VCFv4.0
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
1\t10019\trs1\tA\tC,T\t.\t.\tRS=1;CAF=0.90,0.08,0.02
";
        let mut chrom_map = HashMap::new();
        chrom_map.insert("chr1".to_string(), 0u16);

        let records = parse_dbsnp_vcf(vcf.as_bytes(), &chrom_map).unwrap();
        assert_eq!(records.len(), 2);

        let c = records.iter().find(|r| r.alt_allele == "C").unwrap();
        assert!(c.json.contains("8.000000e-2"), "C should get CAF index 1 (0.08): {}", c.json);

        let t = records.iter().find(|r| r.alt_allele == "T").unwrap();
        assert!(t.json.contains("2.000000e-2"), "T should get CAF index 2 (0.02), not C's frequency: {}", t.json);
    }

    #[test]
    fn test_parse_dbsnp_refseq_accessions() {
        // The real NCBI dbSNP release (GCF_000001405.40) names contigs by
        // RefSeq accession, not `1`/`chr1`. Regression for issue #51: these
        // must resolve when the chromosome map carries the accession key.
        let vcf = "\
##fileformat=VCFv4.0
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
NC_000001.11\t10019\trs775809821\tTA\tT\t.\t.\tRS=775809821;CAF=0.9998,0.0002
NC_000023.11\t100\trs1\tA\tG\t.\t.\tRS=1
";
        let mut chrom_map = HashMap::new();
        chrom_map.insert("NC_000001.11".to_string(), 0u16);
        chrom_map.insert("NC_000023.11".to_string(), 22u16);

        let records = parse_dbsnp_vcf(vcf.as_bytes(), &chrom_map).unwrap();
        assert_eq!(records.len(), 2, "RefSeq-accession contigs were skipped");
        assert_eq!(records[0].chrom_idx, 0);
        assert!(records[0].json.contains("rs775809821"));
        assert_eq!(records[1].chrom_idx, 22);
        assert!(records[1].json.contains("rs1"));
    }
}
