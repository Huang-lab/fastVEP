//! MitoMap mitochondrial variant parser.

use crate::common::{escape_json, AnnotationRecord};
use crate::writer_v2::Osa2Metadata;
use anyhow::{Context, Result};
use std::collections::HashMap;
use std::io::BufRead;

/// Standard MitoMap `.osa2` metadata. The payload is free-text disease and
/// review-status strings, so it is stored as a whole-record JSON blob per
/// variant (see [`crate::writer_v2::raw_json_blob_fields`]) rather than split
/// into numeric columns.
///
/// MitoMap is a few thousand records, so v2 buys no measurable speed or space
/// here. It is wired up so that a `--sa-dir` can be entirely v2: with every
/// variant-level source building to `.osa2`, an operator never has to reason
/// about which of their sources landed in which format, and the v1 reader is
/// needed only to read databases built before the migration.
pub fn mitomap_osa2_metadata(assembly: &str) -> Osa2Metadata {
    Osa2Metadata {
        format_version: 2,
        name: "MitoMap".into(),
        version: "latest".into(),
        assembly: assembly.into(),
        json_key: "mitomap".into(),
        match_by_allele: true,
        is_array: false,
        is_positional: false,
        chunk_bits: 20,
        description: format!("MitoMap mitochondrial variants for {assembly}"),
    }
}

/// Parse a MitoMap variants TSV (pos, ref, alt, disease, status).
pub fn parse_mitomap<R: BufRead>(
    reader: R,
    chrom_to_idx: &HashMap<String, u16>,
) -> Result<Vec<AnnotationRecord>> {
    let mut records = Vec::new();
    let mt_idx = chrom_to_idx
        .get("chrM")
        .or_else(|| chrom_to_idx.get("chrMT"));
    let chrom_idx = match mt_idx {
        Some(&i) => i,
        None => return Ok(records),
    };

    for line in reader.lines() {
        let line = line.context("Reading MitoMap")?;
        if line.starts_with('#') || line.starts_with("Position") || line.is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 4 {
            continue;
        }

        let pos: u32 = match fields[0].trim().parse() {
            Ok(p) => p,
            Err(_) => continue,
        };
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
            chrom_idx,
            position: pos,
            ref_allele,
            alt_allele,
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
