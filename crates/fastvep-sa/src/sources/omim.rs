//! Disease-gene annotation parser for building .oga files.
//!
//! Reads OMIM `genemap2.txt` layout (13-column TSV): gene symbol,
//! identifier, and a semicolon-separated phenotype list. Two real-world
//! inputs use this layout:
//!
//! - **ClinGen Gene-Disease Validity (GDV)** — the SVI-preferred source
//!   (Abou Tayoun 2018; Strehlow et al. 2024). A multi-curator scored
//!   rubric with explicit Definitive/Strong/Moderate classifications.
//!   Convert from the public CSV with `clingen_gdv_to_oga.py`.
//! - **OMIM `genemap2.txt`** — the legacy source (registration-gated at
//!   omim.org). Supported for back-compat; ClinGen GDV is preferred.
//!
//! Both populate the same `omim` json_key. The fastvep PVS1 evaluator
//! consults this data via `OmimData::has_recessive_inheritance()` /
//! `has_dominant_inheritance()` and as a disease-gene fallback when
//! gnomAD constraints don't cross the LOF threshold.

use crate::common::{escape_json, GeneRecord};
use anyhow::{Context, Result};
use std::io::BufRead;

/// Parse OMIM genemap2.txt into GeneRecords.
///
/// Format: tab-separated (0-indexed columns) with columns including:
/// - Column 5: MIM Number
/// - Column 8: Approved Gene Symbol (may be blank for phenotype/QTL-only rows)
/// - Column 12: Phenotypes (semicolon-separated, each with MIM and inheritance)
pub fn parse_omim_genemap<R: BufRead>(reader: R) -> Result<Vec<GeneRecord>> {
    let mut records = Vec::new();

    for line in reader.lines() {
        let line = line.context("Reading OMIM line")?;
        if line.starts_with('#') || line.is_empty() {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 13 {
            continue;
        }

        let mim_number = fields[5]; // MIM Number
        let gene_symbols = fields[8]; // Approved Gene Symbol

        // Use the primary (first) gene symbol
        let gene_symbol = gene_symbols.split(',').next().unwrap_or("").trim();
        if gene_symbol.is_empty() {
            continue;
        }

        let phenotypes_raw = if fields.len() > 12 { fields[12] } else { "" };

        let mut parts = Vec::new();

        if !mim_number.is_empty() && mim_number != "." {
            parts.push(format!("\"mimNumber\":{}", mim_number));
        }

        if !phenotypes_raw.is_empty() && phenotypes_raw != "." {
            let phenotypes: Vec<String> = phenotypes_raw
                .split(';')
                .filter(|p| !p.trim().is_empty())
                .map(|p| {
                    let p = p.trim();
                    format!("\"{}\"", escape_json(p))
                })
                .collect();
            if !phenotypes.is_empty() {
                parts.push(format!("\"phenotypes\":[{}]", phenotypes.join(",")));
            }
        }

        if parts.is_empty() {
            continue;
        }

        records.push(GeneRecord {
            gene_symbol: gene_symbol.to_string(),
            json: format!("{{{}}}", parts.join(",")),
        });
    }

    Ok(records)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_omim() {
        // Real genemap2.txt column layout (0-indexed):
        // 0 Chromosome | 1 Genomic Position Start | 2 Genomic Position End |
        // 3 Cyto Location | 4 Computed Cyto Location | 5 MIM Number |
        // 6 Gene/Locus And Other Related Symbols | 7 Gene Name |
        // 8 Approved Gene Symbol | 9 Entrez Gene ID | 10 Ensembl Gene ID |
        // 11 Comments | 12 Phenotypes | 13 Mouse Gene Symbol/ID
        let data = "# Generated\n\
                     # Copyright OMIM\n\
                     chr17\t43044295\t43125483\t17q21.31\t\t113705\tBRCA1,RNF53\tBreast cancer 1, early onset\tBRCA1\t672\tENSG00000012048\t\tBreast cancer, 114480 (3), Autosomal dominant; Ovarian cancer, 167000 (3)\t\n";

        let records = parse_omim_genemap(data.as_bytes()).unwrap();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].gene_symbol, "BRCA1");
        assert!(records[0].json.contains("113705"));
        assert!(records[0].json.contains("Breast cancer"));
    }

    #[test]
    fn test_parse_omim_skips_rows_without_approved_symbol() {
        // QTL/phenotype-only rows have no Approved Gene Symbol (column 8)
        // but do have MIM Number and aliases; they should be skipped rather
        // than mis-keyed under an alias or MIM number.
        let data = "chr1\t1\t27600000\t1p36\t\t612367\tALPQTL2\tAlkaline phosphatase, plasma level of, QTL 2\t\t100196914\t\t\t\t\n";

        let records = parse_omim_genemap(data.as_bytes()).unwrap();
        assert!(records.is_empty());
    }
}
