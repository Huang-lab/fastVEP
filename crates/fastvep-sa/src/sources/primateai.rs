//! PrimateAI score parser for building .osa annotation files.
//!
//! PrimateAI provides pathogenicity predictions trained on primate variation.
//! Input: TSV with columns chr, pos, ref, alt, primateDL_score.

use crate::common::AnnotationRecord;
use crate::writer_v2::Osa2Metadata;
use anyhow::{Context, Result};
use std::collections::HashMap;
use std::io::BufRead;

/// Standard PrimateAI `.osa2` metadata. Stored as a whole-record JSON blob
/// (see [`crate::writer_v2::raw_json_blob_fields`]): the `{"score":..}` payload
/// rides through byte-for-byte while v2's chunk-level zstd shrinks the database.
pub fn primateai_osa2_metadata(assembly: &str) -> Osa2Metadata {
    Osa2Metadata {
        format_version: 2,
        name: "PrimateAI".into(),
        version: "latest".into(),
        assembly: assembly.into(),
        json_key: "primateAI".into(),
        match_by_allele: true,
        is_array: false,
        is_positional: false,
        chunk_bits: 20,
        description: format!("PrimateAI pathogenicity predictions for {assembly}"),
    }
}

/// Parse a PrimateAI TSV file into sorted AnnotationRecords.
///
/// Expected columns: chr, pos, ref, alt, score (tab-separated).
pub fn parse_primateai<R: BufRead>(
    reader: R,
    chrom_to_idx: &HashMap<String, u16>,
) -> Result<Vec<AnnotationRecord>> {
    let mut records = Vec::new();

    for line in reader.lines() {
        let line = line.context("Reading PrimateAI line")?;
        if line.starts_with('#') || line.starts_with("chr\t") || line.is_empty() {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 5 {
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

        let ref_allele = fields[2].to_string();
        let alt_allele = fields[3].to_string();

        let score: f64 = match fields[4].trim().parse() {
            Ok(s) => s,
            Err(_) => continue,
        };

        let json = format!(r#"{{"score":{:.4}}}"#, score);

        records.push(AnnotationRecord {
            chrom_idx,
            position: pos,
            ref_allele,
            alt_allele,
            json,
        });
    }

    records.sort_by(|a, b| a.chrom_idx.cmp(&b.chrom_idx).then(a.position.cmp(&b.position)));
    Ok(records)
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
    fn test_parse_primateai() {
        let data = "\
chr\tpos\tref\talt\tscore
chr1\t10001\tA\tG\t0.432
chr1\t10001\tA\tC\t0.876
chr1\t10002\tC\tT\t0.123
";
        let mut chrom_map = HashMap::new();
        chrom_map.insert("chr1".into(), 0u16);

        let records = parse_primateai(data.as_bytes(), &chrom_map).unwrap();
        assert_eq!(records.len(), 3);
        assert_eq!(records[0].position, 10001);
        assert_eq!(records[0].alt_allele, "G");
        assert!(records[0].json.contains("0.432"));
    }
}
