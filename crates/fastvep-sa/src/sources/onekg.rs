//! 1000 Genomes population frequency parser.

use crate::common::AnnotationRecord;
use crate::fields::{Field, FieldType};
use crate::writer_v2::Osa2Metadata;
use anyhow::{Context, Result};
use std::collections::HashMap;
use std::io::BufRead;

const POPS: &[&str] = &[
    "AFR", "AMR", "EAS", "EUR", "SAS",
];

fn af_field(alias: &str, description: &str) -> Field {
    Field {
        field: alias.into(),
        alias: alias.into(),
        ftype: FieldType::Float,
        multiplier: 2_000_000,
        zigzag: false,
        missing_value: u32::MAX,
        missing_string: ".".into(),
        description: description.into(),
    }
}

/// Canonical 1000 Genomes `.osa2` field schema. Aliases match the JSON keys
/// the v1 builder emits (`allAf`, `<pop>Af`) so v1 and v2 databases produce
/// identical output. Frequencies render in scientific notation, matching v1.
pub fn onekg_osa2_fields() -> Vec<Field> {
    let mut fields = vec![af_field("allAf", "Global allele frequency")];
    for pop in POPS {
        let lower = pop.to_lowercase();
        fields.push(af_field(
            &format!("{lower}Af"),
            &format!("{pop} allele frequency"),
        ));
    }
    fields
}

/// Standard 1000 Genomes `.osa2` metadata (mirrors the v1 header:
/// `json_key = "oneKg"`, allele-matched, non-positional).
pub fn onekg_osa2_metadata(assembly: &str) -> Osa2Metadata {
    Osa2Metadata {
        format_version: 2,
        name: "1000 Genomes".into(),
        version: "latest".into(),
        assembly: assembly.into(),
        json_key: "oneKg".into(),
        match_by_allele: true,
        is_array: false,
        is_positional: false,
        chunk_bits: 20,
        description: format!("1000 Genomes population frequencies for {assembly}"),
    }
}

/// Parse a 1000 Genomes sites-only VCF into sorted AnnotationRecords.
pub fn parse_onekg_vcf<R: BufRead>(
    reader: R,
    chrom_to_idx: &HashMap<String, u16>,
) -> Result<Vec<AnnotationRecord>> {
    let mut records = Vec::new();

    for line in reader.lines() {
        let line = line.context("Reading 1000G VCF")?;
        if line.starts_with('#') { continue; }

        let fields: Vec<&str> = line.splitn(9, '\t').collect();
        if fields.len() < 8 { continue; }

        let chrom = normalize_chrom(fields[0]);
        let chrom_idx = match chrom_to_idx.get(&chrom) { Some(&i) => i, None => continue };
        let pos: u32 = match fields[1].parse() { Ok(p) => p, Err(_) => continue };
        let ref_allele = fields[3].to_string();
        let alt_field = fields[4];
        let info = fields[7];
        let info_map = parse_info(info);

        let alts: Vec<&str> = alt_field.split(',').collect();
        let all_afs = split_vals(info_map.get("AF").map(|s| s.as_str()));

        for (i, alt) in alts.iter().enumerate() {
            if *alt == "." || *alt == "*" { continue; }

            let mut parts = Vec::new();
            if let Some(af) = all_afs.get(i).and_then(|s| s.parse::<f64>().ok()) {
                parts.push(format!("\"allAf\":{:.6e}", af));
            }
            for pop in POPS {
                let key = format!("{}_AF", pop);
                if let Some(val) = info_map.get(&key) {
                    let vals = split_vals(Some(val.as_str()));
                    if let Some(f) = vals.get(i).and_then(|s| s.parse::<f64>().ok()) {
                        parts.push(format!("\"{}Af\":{:.6e}", pop.to_lowercase(), f));
                    }
                }
            }
            if parts.is_empty() { continue; }
            records.push(AnnotationRecord {
                chrom_idx, position: pos,
                ref_allele: ref_allele.clone(), alt_allele: alt.to_string(),
                json: format!("{{{}}}", parts.join(",")),
            });
        }
    }
    records.sort_by(|a, b| a.chrom_idx.cmp(&b.chrom_idx).then(a.position.cmp(&b.position)));
    Ok(records)
}

fn parse_info(info: &str) -> HashMap<String, String> {
    let mut m = HashMap::new();
    for p in info.split(';') { if let Some((k, v)) = p.split_once('=') { m.insert(k.into(), v.into()); } }
    m
}
fn split_vals(v: Option<&str>) -> Vec<String> {
    v.map(|s| s.split(',').map(|x| x.to_string()).collect()).unwrap_or_default()
}
fn normalize_chrom(c: &str) -> String {
    if c.starts_with("chr") { c.to_string() } else { format!("chr{}", c) }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_parse_onekg() {
        let vcf = "#h\nchr1\t10001\t.\tA\tG\t.\t.\tAF=0.15;AFR_AF=0.20;EUR_AF=0.10\n";
        let mut m = HashMap::new();
        m.insert("chr1".into(), 0u16);
        let recs = parse_onekg_vcf(vcf.as_bytes(), &m).unwrap();
        assert_eq!(recs.len(), 1);
        assert!(recs[0].json.contains("\"allAf\":"));
        assert!(recs[0].json.contains("\"afrAf\":"));
    }

    #[test]
    fn test_onekg_osa2_bridge_encodes_values() {
        let vcf = "#h\nchr1\t10001\t.\tA\tG\t.\t.\tAF=0.15;AFR_AF=0.20;EUR_AF=0.10\n";
        let mut m = HashMap::new();
        m.insert("chr1".into(), 0u16);
        let recs = parse_onekg_vcf(vcf.as_bytes(), &m).unwrap();
        let fields = onekg_osa2_fields();
        let o = crate::writer_v2::osa2_record_from_v1(&recs[0], "chr1".into(), &fields).unwrap();

        let idx = |alias: &str| fields.iter().position(|f| f.alias == alias).unwrap();
        assert_eq!(o.values[idx("allAf")], fields[idx("allAf")].encode_float(0.15));
        assert_eq!(o.values[idx("afrAf")], fields[idx("afrAf")].encode_float(0.20));
        assert_eq!(o.values[idx("eurAf")], fields[idx("eurAf")].encode_float(0.10));
        // Populations absent from the record encode as the missing sentinel.
        assert_eq!(o.values[idx("amrAf")], u32::MAX);
        assert_eq!(o.values[idx("sasAf")], u32::MAX);
    }
}
