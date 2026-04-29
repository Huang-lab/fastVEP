//! ClinVar protein-position index builder for .oga files.
//!
//! Parses ClinVar VCF to extract protein-level data for pathogenic/likely-pathogenic
//! missense variants, enabling PS1, PM5, and PM1 (hotspot) ACMG criteria evaluation.

use crate::common::GeneRecord;
use anyhow::{Context, Result};
use std::collections::HashMap;
use std::io::BufRead;

/// A pathogenic missense variant at a specific protein position.
#[derive(Debug, Clone)]
struct ProteinVariant {
    pos: u64,
    ref_aa: String,
    alt_aa: String,
    sig: String,
}

/// Parse a ClinVar VCF file and produce gene-level records containing
/// pathogenic/likely-pathogenic missense variants indexed by protein position.
///
/// Output: Vec<GeneRecord> where each record's JSON has the structure:
/// `{"proteinVariants":[{"pos":175,"refAa":"R","altAa":"H","sig":"Pathogenic"}, ...]}`
pub fn parse_clinvar_protein_vcf<R: BufRead>(reader: R) -> Result<Vec<GeneRecord>> {
    // Collect pathogenic missense variants per gene
    let mut gene_variants: HashMap<String, Vec<ProteinVariant>> = HashMap::new();

    for line in reader.lines() {
        let line = line.context("Reading ClinVar VCF line")?;
        if line.starts_with('#') {
            continue;
        }

        let fields: Vec<&str> = line.splitn(9, '\t').collect();
        if fields.len() < 8 {
            continue;
        }

        let info = fields[7];
        let info_map = parse_info(info);

        // Only process pathogenic/likely-pathogenic variants
        let clnsig = match info_map.get("CLNSIG") {
            Some(sig) => sig.to_lowercase(),
            None => continue,
        };
        if !clnsig.contains("pathogenic") || clnsig.contains("conflicting") {
            continue;
        }

        // Check molecular consequence for missense
        let mc = info_map.get("MC").map(|s| s.as_str()).unwrap_or("");
        if !mc.contains("missense_variant") {
            continue;
        }

        // Extract gene symbol from GENEINFO=GENE:ID
        let gene_symbol = match info_map.get("GENEINFO") {
            Some(gi) => {
                if let Some(colon_pos) = gi.find(':') {
                    gi[..colon_pos].to_string()
                } else {
                    gi.to_string()
                }
            }
            None => continue,
        };

        // Try to extract protein change from MC field or CLNHGVS
        // MC format: SO:0001583|missense_variant  (no protein info)
        // CLNHGVS format: NC_000017.11:g.7676154G>A (genomic, not protein)
        // We need to parse protein change from alternative fields

        // Try CLNHGVS for protein notation (some entries have multiple, pipe-separated)
        let mut protein_variant = None;

        // Parse from MC field which sometimes includes protein change info
        // Format can be: "SO:0001583|missense_variant|NP_000537.3:p.Arg175His"
        for part in mc.split(',') {
            if let Some(prot) = extract_protein_from_mc(part) {
                protein_variant = Some(prot);
                break;
            }
        }

        // Also try parsing from CLNHGVS (may contain p. notation)
        if protein_variant.is_none() {
            if let Some(hgvs) = info_map.get("CLNHGVS") {
                for part in hgvs.split(',') {
                    if let Some(prot) = parse_protein_hgvs(part) {
                        protein_variant = Some(prot);
                        break;
                    }
                }
            }
        }

        if let Some(pv) = protein_variant {
            let sig_clean = if clnsig.contains("likely") {
                "Likely_pathogenic".to_string()
            } else {
                "Pathogenic".to_string()
            };

            gene_variants
                .entry(gene_symbol)
                .or_default()
                .push(ProteinVariant {
                    pos: pv.0,
                    ref_aa: pv.1,
                    alt_aa: pv.2,
                    sig: sig_clean,
                });
        }
    }

    // Convert to GeneRecords
    let mut records: Vec<GeneRecord> = gene_variants
        .into_iter()
        .map(|(gene, variants)| {
            // Deduplicate by (pos, ref_aa, alt_aa)
            let mut unique: HashMap<(u64, String, String), String> = HashMap::new();
            for v in &variants {
                unique
                    .entry((v.pos, v.ref_aa.clone(), v.alt_aa.clone()))
                    .or_insert_with(|| v.sig.clone());
            }

            let variant_jsons: Vec<String> = unique
                .iter()
                .map(|((pos, ref_aa, alt_aa), sig)| {
                    format!(
                        r#"{{"pos":{},"refAa":"{}","altAa":"{}","sig":"{}"}}"#,
                        pos, ref_aa, alt_aa, sig
                    )
                })
                .collect();

            let json = format!(r#"{{"proteinVariants":[{}]}}"#, variant_jsons.join(","));

            GeneRecord {
                gene_symbol: gene,
                json,
            }
        })
        .collect();

    records.sort_by(|a, b| a.gene_symbol.cmp(&b.gene_symbol));
    Ok(records)
}

/// Extract protein position and amino acid change from MC field component.
/// MC can contain entries like: "SO:0001583|missense_variant" (no protein info usually)
/// but some ClinVar entries have extended format.
fn extract_protein_from_mc(mc_part: &str) -> Option<(u64, String, String)> {
    // Look for p. notation in the MC field
    if let Some(p_idx) = mc_part.find(":p.") {
        return parse_protein_hgvs(&mc_part[p_idx + 1..]);
    }
    None
}

/// Parse a protein HGVS expression like "p.Arg175His" or "p.R175H"
/// Returns (position, ref_aa, alt_aa) using single-letter codes.
fn parse_protein_hgvs(hgvs: &str) -> Option<(u64, String, String)> {
    let p_str = if let Some(idx) = hgvs.find("p.") {
        &hgvs[idx + 2..]
    } else {
        return None;
    };

    // Try three-letter codes first: "Arg175His"
    if let Some(result) = parse_three_letter_protein(p_str) {
        return Some(result);
    }

    // Try single-letter codes: "R175H"
    if p_str.len() >= 3 {
        let first = p_str.chars().next()?;
        if first.is_ascii_uppercase() {
            // Extract digits
            let digits: String = p_str[1..].chars().take_while(|c| c.is_ascii_digit()).collect();
            if let Ok(pos) = digits.parse::<u64>() {
                let rest = &p_str[1 + digits.len()..];
                if let Some(alt_aa) = rest.chars().next() {
                    if alt_aa.is_ascii_uppercase() && alt_aa != '*' {
                        return Some((pos, first.to_string(), alt_aa.to_string()));
                    }
                }
            }
        }
    }

    None
}

/// Parse three-letter amino acid protein change like "Arg175His"
fn parse_three_letter_protein(s: &str) -> Option<(u64, String, String)> {
    let aa_map: HashMap<&str, &str> = [
        ("Ala", "A"), ("Arg", "R"), ("Asn", "N"), ("Asp", "D"), ("Cys", "C"),
        ("Gln", "Q"), ("Glu", "E"), ("Gly", "G"), ("His", "H"), ("Ile", "I"),
        ("Leu", "L"), ("Lys", "K"), ("Met", "M"), ("Phe", "F"), ("Pro", "P"),
        ("Ser", "S"), ("Thr", "T"), ("Trp", "W"), ("Tyr", "Y"), ("Val", "V"),
        ("Sec", "U"), ("Pyl", "O"), ("Ter", "*"),
    ]
    .iter()
    .copied()
    .collect();

    // Find ref AA (first 3 chars)
    if s.len() < 4 {
        return None;
    }
    let ref_three = &s[..3];
    let ref_aa = aa_map.get(ref_three)?;

    // Extract position digits
    let rest = &s[3..];
    let digits: String = rest.chars().take_while(|c| c.is_ascii_digit()).collect();
    let pos = digits.parse::<u64>().ok()?;

    // Find alt AA
    let after_digits = &rest[digits.len()..];
    if after_digits.len() < 3 {
        return None;
    }
    let alt_three = &after_digits[..3];
    let alt_aa = aa_map.get(alt_three)?;

    // Skip stop codon/terminator variants (not missense)
    if *alt_aa == "*" {
        return None;
    }

    Some((pos, ref_aa.to_string(), alt_aa.to_string()))
}

fn parse_info(info: &str) -> HashMap<String, String> {
    let mut map = HashMap::new();
    for part in info.split(';') {
        if let Some(eq_pos) = part.find('=') {
            let key = &part[..eq_pos];
            let val = &part[eq_pos + 1..];
            map.insert(key.to_string(), val.to_string());
        }
    }
    map
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_protein_hgvs_three_letter() {
        let result = parse_protein_hgvs("p.Arg175His").unwrap();
        assert_eq!(result, (175, "R".to_string(), "H".to_string()));
    }

    #[test]
    fn test_parse_protein_hgvs_single_letter() {
        let result = parse_protein_hgvs("p.R175H").unwrap();
        assert_eq!(result, (175, "R".to_string(), "H".to_string()));
    }

    #[test]
    fn test_parse_protein_hgvs_with_prefix() {
        let result = parse_protein_hgvs("NP_000537.3:p.Arg175His").unwrap();
        assert_eq!(result, (175, "R".to_string(), "H".to_string()));
    }

    #[test]
    fn test_parse_protein_hgvs_stop_codon_rejected() {
        assert!(parse_protein_hgvs("p.Arg175Ter").is_none());
    }

    #[test]
    fn test_parse_three_letter_protein() {
        let result = parse_three_letter_protein("Cys315Met").unwrap();
        assert_eq!(result, (315, "C".to_string(), "M".to_string()));
    }
}
