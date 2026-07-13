//! MitoMap mitochondrial variant parser.

use crate::common::{escape_json, AnnotationRecord};
use anyhow::{Context, Result};
use std::collections::HashMap;
use std::io::BufRead;

/// Parse a MitoMap variants TSV (pos, ref, alt, disease, status).
pub fn parse_mitomap<R: BufRead>(
    reader: R,
    chrom_to_idx: &HashMap<String, u16>,
) -> Result<Vec<AnnotationRecord>> {
    let mut records = Vec::new();
    let mt_idx = chrom_to_idx.get("chrM").or_else(|| chrom_to_idx.get("chrMT"));
    let chrom_idx = match mt_idx { Some(&i) => i, None => return Ok(records) };

    for line in reader.lines() {
        let line = line.context("Reading MitoMap")?;
        if line.starts_with('#') || line.starts_with("Position") || line.is_empty() { continue; }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 4 { continue; }

        let pos: u32 = match fields[0].trim().parse() { Ok(p) => p, Err(_) => continue };
        let ref_allele = fields[1].trim().to_string();
        let alt_allele = fields[2].trim().to_string();
        let disease = fields.get(3).unwrap_or(&"").trim();

        let mut parts = Vec::new();
        if !disease.is_empty() && disease != "." {
            parts.push(format!("\"disease\":\"{}\"", escape_json(disease)));
        }
        if let Some(status) = fields.get(4) {
            let status = status.trim();
            if !status.is_empty() && status != "." {
                parts.push(format!("\"status\":\"{}\"", escape_json(status)));
            }
        }

        records.push(AnnotationRecord {
            chrom_idx, position: pos,
            ref_allele, alt_allele,
            json: format!("{{{}}}", parts.join(",")),
        });
    }
    records.sort_by_key(|r| r.position);
    Ok(records)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_parse_mitomap() {
        let data = "Position\tRef\tAlt\tDisease\tStatus\n3243\tA\tG\tMELAS\tConfirmed\n";
        let mut m = HashMap::new();
        m.insert("chrM".into(), 24u16);
        let recs = parse_mitomap(data.as_bytes(), &m).unwrap();
        assert_eq!(recs.len(), 1);
        assert_eq!(recs[0].position, 3243);
        assert!(recs[0].json.contains("MELAS"));
    }
}
