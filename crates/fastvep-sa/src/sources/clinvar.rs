//! ClinVar VCF parser for building .osa annotation files.
//!
//! Parses ClinVar's VCF release to extract clinical significance,
//! review status, disease names, and accession numbers.

use crate::common::AnnotationRecord;
use anyhow::{Context, Result};
use std::collections::HashMap;
use std::io::BufRead;

/// Parse a ClinVar VCF file and produce sorted `AnnotationRecord`s.
///
/// The VCF must be from NCBI's ClinVar release (clinvar.vcf.gz).
/// Records are sorted by (chrom_idx, position) for `SaWriter`.
pub fn parse_clinvar_vcf<R: BufRead>(
    reader: R,
    chrom_to_idx: &HashMap<String, u16>,
) -> Result<Vec<AnnotationRecord>> {
    let mut records = Vec::new();

    for line in reader.lines() {
        let line = line.context("Reading ClinVar VCF line")?;
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
            None => continue, // Skip unrecognized chromosomes
        };

        let pos: u32 = match fields[1].parse() {
            Ok(p) => p,
            Err(_) => continue,
        };

        // ClinVar's VCF release carries the integer VariationID in the ID column
        // (column 3); rsIDs, when present, live in the RS INFO field instead.
        let variation_id = fields[2];
        let ref_allele = fields[3].to_string();
        let alt_field = fields[4];

        // Parse INFO field
        let info = fields[7];
        let info_map = parse_info(info);

        let clnsig = info_map.get("CLNSIG").cloned().unwrap_or_default();
        let clnrevstat = info_map.get("CLNREVSTAT").cloned().unwrap_or_default();
        let clndn = info_map.get("CLNDN").cloned().unwrap_or_default();
        let clnacc = info_map.get("CLNVC").cloned(); // variant class
        let clnid = info_map.get("CLNVCSO").cloned(); // SO accession
        let clnsigconf = info_map.get("CLNSIGCONF").cloned(); // conflicting-significance breakdown
        let clndisdb = info_map.get("CLNDISDB").cloned(); // disease-database cross-references
        let mc = info_map.get("MC").cloned(); // molecular consequence (SO term | label)
        // ClinVar-distributed population allele frequencies (ExAC / 1000G / ESP).
        // Used as a frequency backstop by PM2 when gnomAD has no record.
        let af_exac = info_map.get("AF_EXAC").and_then(|s| s.parse::<f64>().ok());
        let af_tgp = info_map.get("AF_TGP").and_then(|s| s.parse::<f64>().ok());
        let af_esp = info_map.get("AF_ESP").and_then(|s| s.parse::<f64>().ok());

        // Handle multi-allelic: each ALT gets its own record
        for alt in alt_field.split(',') {
            if alt == "." || alt == "*" {
                continue;
            }

            let json = build_clinvar_json(
                variation_id,
                &clnsig,
                &clnrevstat,
                &clndn,
                clnacc.as_deref(),
                clnid.as_deref(),
                clnsigconf.as_deref(),
                clndisdb.as_deref(),
                mc.as_deref(),
                af_exac,
                af_tgp,
                af_esp,
            );

            records.push(AnnotationRecord {
                chrom_idx,
                position: pos,
                ref_allele: ref_allele.clone(),
                alt_allele: alt.to_string(),
                json,
            });
        }
    }

    // Sort by (chrom_idx, position) as required by SaWriter
    records.sort_by(|a, b| a.chrom_idx.cmp(&b.chrom_idx).then(a.position.cmp(&b.position)));

    Ok(records)
}

#[allow(clippy::too_many_arguments)]
fn build_clinvar_json(
    variation_id: &str,
    clnsig: &str,
    clnrevstat: &str,
    clndn: &str,
    clnvc: Option<&str>,
    clnvcso: Option<&str>,
    clnsigconf: Option<&str>,
    clndisdb: Option<&str>,
    mc: Option<&str>,
    af_exac: Option<f64>,
    af_tgp: Option<f64>,
    af_esp: Option<f64>,
) -> String {
    let mut parts = Vec::new();

    if !clnsig.is_empty() {
        let sigs: Vec<String> = clnsig
            .split('|')
            .map(|s| format!("\"{}\"", escape_json(s)))
            .collect();
        parts.push(format!("\"significance\":[{}]", sigs.join(",")));
    }

    if !clnrevstat.is_empty() {
        parts.push(format!(
            "\"reviewStatus\":\"{}\"",
            escape_json(clnrevstat)
        ));
    }

    if !clndn.is_empty() && clndn != "not_provided" {
        let diseases: Vec<String> = clndn
            .split('|')
            .filter(|s| *s != "not_provided")
            .map(|s| format!("\"{}\"", escape_json(s)))
            .collect();
        if !diseases.is_empty() {
            parts.push(format!("\"phenotypes\":[{}]", diseases.join(",")));
        }
    }

    if let Some(vc) = clnvc {
        parts.push(format!("\"variantClass\":\"{}\"", escape_json(vc)));
    }

    if let Some(vcso) = clnvcso {
        parts.push(format!("\"soAccession\":\"{}\"", escape_json(vcso)));
    }

    // Population AFs are finite floats (parsed via f64::from_str); Display emits
    // plain decimal JSON numbers (never NaN/inf), so this stays valid JSON.
    if let Some(af) = af_exac {
        parts.push(format!("\"afExac\":{}", af));
    }
    if let Some(af) = af_tgp {
        parts.push(format!("\"afTgp\":{}", af));
    }
    if let Some(af) = af_esp {
        parts.push(format!("\"afEsp\":{}", af));
    }

    // VariationID (the VCF ID column). Emit only when present (non-empty, not ".").
    if !variation_id.is_empty() && variation_id != "." {
        parts.push(format!("\"variationId\":\"{}\"", escape_json(variation_id)));
    }

    // Review confidence as 0-4 "gold stars" derived from CLNREVSTAT. 0 is a
    // meaningful "no assertion" signal, so it is emitted rather than suppressed.
    if !clnrevstat.is_empty() {
        parts.push(format!(
            "\"goldStars\":{}",
            clnrevstat_to_gold_stars(clnrevstat)
        ));
    }

    // Auxiliary, curator-facing detail — JSON-only (not projected to the pipe).
    if let Some(c) = clnsigconf {
        parts.push(format!("\"clnSigConf\":\"{}\"", escape_json(c)));
    }
    if let Some(d) = clndisdb {
        if d != "not_provided" {
            parts.push(format!("\"clnDisDb\":\"{}\"", escape_json(d)));
        }
    }
    if let Some(m) = mc {
        parts.push(format!("\"molecularConsequence\":\"{}\"", escape_json(m)));
    }

    format!("{{{}}}", parts.join(","))
}

/// Map ClinVar `CLNREVSTAT` to the 0-4 "gold star" review-confidence scale.
///
/// Checked most-confident first; the 2-star "multiple_submitters,_no_conflicts"
/// case is matched before the 1-star conflicting / single-submitter cases.
fn clnrevstat_to_gold_stars(clnrevstat: &str) -> u8 {
    if clnrevstat.contains("practice_guideline") {
        4
    } else if clnrevstat.contains("reviewed_by_expert_panel") {
        3
    } else if clnrevstat.contains("criteria_provided,_multiple_submitters,_no_conflicts") {
        2
    } else if clnrevstat.contains("criteria_provided,_conflicting")
        || clnrevstat.contains("criteria_provided,_single_submitter")
    {
        1
    } else {
        0
    }
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
    if chrom.starts_with("chr") {
        chrom.to_string()
    } else {
        format!("chr{}", chrom)
    }
}

fn escape_json(s: &str) -> String {
    s.replace('\\', "\\\\")
        .replace('"', "\\\"")
        .replace('\n', "\\n")
        .replace('\r', "\\r")
        .replace('\t', "\\t")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_clinvar_vcf() {
        let vcf = "\
##fileformat=VCFv4.1
##INFO=<ID=CLNSIG,Number=.,Type=String>
##INFO=<ID=CLNREVSTAT,Number=.,Type=String>
##INFO=<ID=CLNDN,Number=.,Type=String>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
1\t12345\t12500\tA\tG\t.\t.\tCLNSIG=Pathogenic;CLNREVSTAT=reviewed_by_expert_panel;CLNDN=Breast_cancer;CLNSIGCONF=Pathogenic(3);CLNDISDB=MedGen:C0006142,OMIM:114480;MC=SO:0001583|missense_variant
1\t67890\t67900\tC\tT\t.\t.\tCLNSIG=Benign;CLNREVSTAT=criteria_provided,_single_submitter;CLNDN=not_provided
";

        let mut chrom_map = HashMap::new();
        chrom_map.insert("chr1".to_string(), 0u16);

        let records = parse_clinvar_vcf(vcf.as_bytes(), &chrom_map).unwrap();
        assert_eq!(records.len(), 2);

        assert_eq!(records[0].position, 12345);
        assert!(records[0].json.contains("Pathogenic"));
        assert!(records[0].json.contains("Breast_cancer"));
        // VariationID is the integer VCF ID column (not an rsID).
        assert!(records[0].json.contains("\"variationId\":\"12500\""));
        // reviewed_by_expert_panel -> 3 gold stars.
        assert!(records[0].json.contains("\"goldStars\":3"));
        assert!(records[0].json.contains("\"clnSigConf\":\"Pathogenic(3)\""));
        assert!(records[0].json.contains("OMIM:114480"));
        assert!(records[0].json.contains("missense_variant"));

        assert_eq!(records[1].position, 67890);
        assert!(records[1].json.contains("Benign"));
        assert!(records[1].json.contains("\"variationId\":\"67900\""));
        // single_submitter -> 1 gold star.
        assert!(records[1].json.contains("\"goldStars\":1"));
        // "not_provided" should be filtered out
        assert!(!records[1].json.contains("phenotypes"));
    }

    #[test]
    fn test_parse_clinvar_population_af() {
        // ExAC / 1000G / ESP allele frequencies are emitted as numeric JSON
        // (afExac / afTgp / afEsp) for the PM2 frequency backstop. Variants
        // without these INFO keys must not emit the keys at all.
        let vcf = "\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
1\t12345\trs1\tA\tG\t.\t.\tCLNSIG=Pathogenic;AF_EXAC=0.00054;AF_TGP=0.0016
1\t67890\trs2\tC\tT\t.\t.\tCLNSIG=Pathogenic
";
        let mut chrom_map = HashMap::new();
        chrom_map.insert("chr1".to_string(), 0u16);

        let records = parse_clinvar_vcf(vcf.as_bytes(), &chrom_map).unwrap();
        assert_eq!(records.len(), 2);
        assert!(records[0].json.contains("\"afExac\":0.00054"));
        assert!(records[0].json.contains("\"afTgp\":0.0016"));
        assert!(!records[0].json.contains("afEsp"));
        // No AF INFO at all → no AF keys emitted.
        assert!(!records[1].json.contains("afExac"));
        assert!(!records[1].json.contains("afTgp"));
        assert!(!records[1].json.contains("afEsp"));
    }

    #[test]
    fn test_clnrevstat_to_gold_stars() {
        assert_eq!(clnrevstat_to_gold_stars("practice_guideline"), 4);
        assert_eq!(clnrevstat_to_gold_stars("reviewed_by_expert_panel"), 3);
        assert_eq!(
            clnrevstat_to_gold_stars("criteria_provided,_multiple_submitters,_no_conflicts"),
            2
        );
        assert_eq!(
            clnrevstat_to_gold_stars("criteria_provided,_conflicting_classifications"),
            1
        );
        assert_eq!(
            clnrevstat_to_gold_stars("criteria_provided,_single_submitter"),
            1
        );
        assert_eq!(clnrevstat_to_gold_stars("no_assertion_criteria_provided"), 0);
        assert_eq!(clnrevstat_to_gold_stars(""), 0);
    }
}
