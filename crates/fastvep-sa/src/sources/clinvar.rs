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
        // ClinVar-distributed population allele frequencies (ExAC / 1000G / ESP).
        // Used as a frequency backstop by PM2 when gnomAD has no record. Reject
        // non-finite values: `f64::from_str` accepts "inf"/"nan", whose Display
        // form is invalid JSON and would poison the entire record's JSON string
        // (silently dropping its ClinVar annotation). Real ClinVar AFs are
        // always finite decimals, so this only guards against malformed input.
        let parse_af = |k: &str| {
            info_map
                .get(k)
                .and_then(|s| s.parse::<f64>().ok())
                .filter(|v| v.is_finite())
        };
        let af_exac = parse_af("AF_EXAC");
        let af_tgp = parse_af("AF_TGP");
        let af_esp = parse_af("AF_ESP");

        // Handle multi-allelic: each ALT gets its own record
        for alt in alt_field.split(',') {
            if alt == "." || alt == "*" {
                continue;
            }

            let json = build_clinvar_json(
                &clnsig,
                &clnrevstat,
                &clndn,
                clnacc.as_deref(),
                clnid.as_deref(),
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
    clnsig: &str,
    clnrevstat: &str,
    clndn: &str,
    clnvc: Option<&str>,
    clnvcso: Option<&str>,
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

    format!("{{{}}}", parts.join(","))
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
1\t12345\trs123\tA\tG\t.\t.\tCLNSIG=Pathogenic;CLNREVSTAT=criteria_provided,_multiple_submitters,_no_conflicts;CLNDN=Breast_cancer
1\t67890\trs456\tC\tT\t.\t.\tCLNSIG=Benign;CLNREVSTAT=criteria_provided,_single_submitter;CLNDN=not_provided
";

        let mut chrom_map = HashMap::new();
        chrom_map.insert("chr1".to_string(), 0u16);

        let records = parse_clinvar_vcf(vcf.as_bytes(), &chrom_map).unwrap();
        assert_eq!(records.len(), 2);

        assert_eq!(records[0].position, 12345);
        assert!(records[0].json.contains("Pathogenic"));
        assert!(records[0].json.contains("Breast_cancer"));

        assert_eq!(records[1].position, 67890);
        assert!(records[1].json.contains("Benign"));
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
    fn test_parse_clinvar_non_finite_af_dropped() {
        // `f64::from_str` accepts "inf"/"nan"; their Display form is invalid
        // JSON and would poison the whole record. A malformed AF must be
        // dropped (treated as absent), leaving the rest of the record intact
        // and the JSON parseable.
        let vcf = "\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
1\t12345\trs1\tA\tG\t.\t.\tCLNSIG=Pathogenic;AF_EXAC=inf;AF_TGP=nan;AF_ESP=0.0016
";
        let mut chrom_map = HashMap::new();
        chrom_map.insert("chr1".to_string(), 0u16);

        let records = parse_clinvar_vcf(vcf.as_bytes(), &chrom_map).unwrap();
        assert_eq!(records.len(), 1);
        // Non-finite AFs dropped; the finite one survives.
        assert!(!records[0].json.contains("afExac"));
        assert!(!records[0].json.contains("afTgp"));
        assert!(records[0].json.contains("\"afEsp\":0.0016"));
        // The emitted JSON must still be valid.
        let v: serde_json::Value = serde_json::from_str(&records[0].json).unwrap();
        assert!(v.get("afExac").is_none());
        assert_eq!(v.get("afEsp").and_then(|x| x.as_f64()), Some(0.0016));
    }
}
