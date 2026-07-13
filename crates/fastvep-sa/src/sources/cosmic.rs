//! COSMIC VCF/TSV parser for building .osa annotation files.
//!
//! Extracts somatic mutation data from COSMIC's coding mutations file.

use crate::common::{escape_json, AnnotationRecord};
use anyhow::{Context, Result};
use std::collections::HashMap;
use std::io::BufRead;

/// Parse a COSMIC VCF (CosmicCodingMuts.vcf) into sorted AnnotationRecords.
pub fn parse_cosmic_vcf<R: BufRead>(
    reader: R,
    chrom_to_idx: &HashMap<String, u16>,
) -> Result<Vec<AnnotationRecord>> {
    let mut records = Vec::new();

    for line in reader.lines() {
        let line = line.context("Reading COSMIC VCF")?;
        if line.starts_with('#') {
            continue;
        }

        let fields: Vec<&str> = line.splitn(9, '\t').collect();
        if fields.len() < 8 {
            continue;
        }

        let chrom = normalize_chrom(fields[0]);
        let chrom_idx = match chrom_to_idx.get(&chrom) {
            Some(&idx) => idx,
            None => continue,
        };

        let pos: u32 = match fields[1].parse() {
            Ok(p) => p,
            Err(_) => continue,
        };

        let id = fields[2]; // COSV ID
        let ref_allele = fields[3].to_string();
        let alt_field = fields[4];
        let info = fields[7];

        let info_map = parse_info(info);
        let gene = info_map.get("GENE").cloned().unwrap_or_default();
        let cnt = info_map.get("CNT").cloned().unwrap_or_default();

        for alt in alt_field.split(',') {
            if alt == "." || alt == "*" {
                continue;
            }

            let mut parts = Vec::new();
            if !id.is_empty() && id != "." {
                parts.push(format!("\"id\":\"{}\"", escape_json(id)));
            }
            if !gene.is_empty() && gene != "." {
                parts.push(format!("\"gene\":\"{}\"", escape_json(&gene)));
            }
            // CNT is written unquoted, so it must be a validated integer —
            // the VCF missing-value sentinel "." (or any other garbage)
            // would otherwise land in the JSON as a bare, unquoted token
            // and break every downstream serde_json::from_str on this record.
            if let Ok(count) = cnt.parse::<u64>() {
                parts.push(format!("\"count\":{}", count));
            }

            records.push(AnnotationRecord {
                chrom_idx,
                position: pos,
                ref_allele: ref_allele.clone(),
                alt_allele: alt.to_string(),
                json: format!("{{{}}}", parts.join(",")),
            });
        }
    }

    records.sort_by(|a, b| a.chrom_idx.cmp(&b.chrom_idx).then(a.position.cmp(&b.position)));
    Ok(records)
}

fn parse_info(info: &str) -> HashMap<String, String> {
    let mut map = HashMap::new();
    for pair in info.split(';') {
        if let Some((k, v)) = pair.split_once('=') {
            map.insert(k.to_string(), v.to_string());
        }
    }
    map
}

fn normalize_chrom(chrom: &str) -> String {
    if chrom.starts_with("chr") { chrom.to_string() } else { format!("chr{}", chrom) }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_cosmic() {
        let vcf = "#header\n1\t10001\tCOSV123\tA\tG\t.\t.\tGENE=TP53;CNT=15\n";
        let mut m = HashMap::new();
        m.insert("chr1".into(), 0u16);
        let records = parse_cosmic_vcf(vcf.as_bytes(), &m).unwrap();
        assert_eq!(records.len(), 1);
        assert!(records[0].json.contains("COSV123"));
        assert!(records[0].json.contains("TP53"));
    }

    #[test]
    fn test_parse_cosmic_escapes_gene_field_for_valid_json() {
        // A GENE value containing a double quote or backslash must not
        // produce invalid JSON in the .osa record.
        let vcf = "#header\n1\t10001\tCOSV123\tA\tG\t.\t.\tGENE=WEIRD\"GENE\\NAME;CNT=1\n";
        let mut m = HashMap::new();
        m.insert("chr1".into(), 0u16);
        let records = parse_cosmic_vcf(vcf.as_bytes(), &m).unwrap();
        assert_eq!(records.len(), 1);
        let parsed: serde_json::Value = serde_json::from_str(&records[0].json)
            .expect("COSMIC record must be valid JSON even with quotes/backslashes in GENE");
        assert_eq!(parsed["gene"], "WEIRD\"GENE\\NAME");
    }

    #[test]
    fn test_parse_cosmic_missing_value_sentinel_is_omitted_not_emitted_raw() {
        // GENE="." (unmapped/intergenic) and CNT="." (missing count) are
        // VCF's missing-value sentinel, not real values. GENE="." must not
        // be stored as a literal gene symbol, and CNT="." must not be
        // written unquoted into JSON (which would make it invalid).
        let vcf = "#header\n1\t10001\tCOSV123\tA\tG\t.\t.\tGENE=.;CNT=.\n";
        let mut m = HashMap::new();
        m.insert("chr1".into(), 0u16);
        let records = parse_cosmic_vcf(vcf.as_bytes(), &m).unwrap();
        assert_eq!(records.len(), 1);
        let json = &records[0].json;
        let _parsed: serde_json::Value =
            serde_json::from_str(json).expect("record must still be valid JSON: {json}");
        assert!(!json.contains("gene"), "GENE=. must be omitted: {json}");
        assert!(!json.contains("count"), "CNT=. must be omitted: {json}");
    }
}
