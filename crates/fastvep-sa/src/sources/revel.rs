//! REVEL score parser for building .osa annotation files.
//!
//! REVEL provides missense pathogenicity predictions as allele-specific scores.
//! Input format: CSV with columns chr, hg19_pos, grch38_pos, ref, alt, REVEL.

use crate::common::AnnotationRecord;
use crate::writer_v2::Osa2Metadata;
use anyhow::{Context, Result};
use std::collections::HashMap;
use std::io::BufRead;

/// Standard REVEL `.osa2` metadata. REVEL's payload is a single `{"score":..}`
/// object per allele; it is stored as a whole-record JSON blob (see
/// [`crate::writer_v2::raw_json_blob_fields`]) so v2 output is byte-identical
/// to v1 — the score's fixed-decimal text rides through untouched — while v2's
/// chunk-level zstd of the blob column shrinks the database.
pub fn revel_osa2_metadata(assembly: &str) -> Osa2Metadata {
    Osa2Metadata {
        format_version: 2,
        name: "REVEL".into(),
        version: "latest".into(),
        assembly: assembly.into(),
        json_key: "revel".into(),
        match_by_allele: true,
        is_array: false,
        is_positional: false,
        chunk_bits: 20,
        description: format!("REVEL missense pathogenicity scores for {assembly}"),
    }
}

/// Column of the REVEL CSV that holds the coordinate for `assembly`.
///
/// One file carries both builds side by side, so indexing the wrong column
/// does not fail: a GRCh37 database built from `grch38_pos` loads, answers,
/// and returns the score of a different base (#133). Anything that is not
/// GRCh37/hg19 reads `grch38_pos`, which is what every earlier release did.
pub fn position_column(assembly: &str) -> usize {
    match assembly.to_ascii_lowercase().as_str() {
        "grch37" | "hg19" => 1,
        _ => 2,
    }
}

/// Parse a REVEL score file (CSV) into sorted AnnotationRecords.
///
/// REVEL distributes scores as CSV: chr, hg19_pos, grch38_pos, ref, alt, aaref, aaalt, REVEL
/// `pos_column` picks which of the two coordinate columns is indexed; use
/// [`position_column`] to derive it from the build's `--assembly`.
pub fn parse_revel<R: BufRead>(
    reader: R,
    chrom_to_idx: &HashMap<String, u16>,
    pos_column: usize,
) -> Result<Vec<AnnotationRecord>> {
    let mut records = Vec::new();

    for line in reader.lines() {
        let line = line.context("Reading REVEL line")?;
        if line.starts_with("chr,") || line.starts_with('#') || line.is_empty() {
            continue; // Skip header
        }

        let fields: Vec<&str> = line.split(',').collect();
        if fields.len() < 8 {
            continue;
        }

        let chrom = normalize_chrom(fields[0]);
        let chrom_idx = match chrom_to_idx.get(&chrom) {
            Some(&idx) => idx,
            None => continue,
        };

        let pos: u32 = match fields[pos_column].parse() {
            Ok(p) => p,
            Err(_) => continue,
        };

        let ref_allele = fields[3].to_string();
        let alt_allele = fields[4].to_string();

        let score: f64 = match fields[7].trim().parse() {
            Ok(s) => s,
            Err(_) => continue,
        };

        let json = format!(r#"{{"score":{:.4}}}"#, score).replace(".0000}", ".0}");

        records.push(AnnotationRecord {
            chrom_idx,
            position: pos,
            ref_allele,
            alt_allele,
            json,
        });
    }

    records.sort_by(|a, b| {
        a.chrom_idx
            .cmp(&b.chrom_idx)
            .then(a.position.cmp(&b.position))
    });
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
    fn position_column_follows_the_assembly() {
        assert_eq!(position_column("GRCh38"), 2);
        assert_eq!(position_column("GRCh37"), 1);
        assert_eq!(position_column("hg19"), 1);
    }

    #[test]
    fn test_parse_revel() {
        let data = "\
chr,hg19_pos,grch38_pos,ref,alt,aaref,aaalt,REVEL
1,35142,35142,G,A,T,M,0.027
1,35142,35142,G,C,T,S,0.035
1,35143,35143,C,A,T,N,0.842
";
        let mut chrom_map = HashMap::new();
        chrom_map.insert("chr1".into(), 0u16);

        let records = parse_revel(data.as_bytes(), &chrom_map, 2).unwrap();
        assert_eq!(records.len(), 3);
        assert_eq!(records[0].position, 35142);
        assert_eq!(records[0].ref_allele, "G");
        assert_eq!(records[0].alt_allele, "A");
        assert!(records[0].json.contains("0.027"));
        assert_eq!(records[2].position, 35143);
        assert!(records[2].json.contains("0.842"));
    }
}
