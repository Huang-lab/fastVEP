//! gnomAD VCF parser for building .osa annotation files.
//!
//! Parses gnomAD's sites-only VCF to extract allele frequencies
//! per population, allele counts, and homozygote counts.
//!
//! Auto-detects the INFO field-naming convention used by the input:
//!
//! | Release | Top-level | Per-population |
//! |---------|-----------|----------------|
//! | v2.x (exomes/genomes) | `AF` / `AC` / `AN` / `nhomalt` | `AF_afr` (underscore) |
//! | v3.x (genomes)        | `AF` / `AC` / `AN` / `nhomalt` | none at top level (only subset-qualified, e.g. `AF-non_neuro-afr`) |
//! | v4.0/v4.1 (exomes/genomes) | `AF` / `AC` / `AN` / `nhomalt` | `AF_afr` (underscore) |
//! | v4.1 joint            | `AF_joint` / `AC_joint` / `AN_joint` / `nhomalt_joint` | `AF_joint_afr` |
//!
//! Some locally-processed files expose simple per-population fields with a
//! hyphen separator (`AF-afr`). The parser detects and handles this case.
//!
//! We pick the scheme by scanning `##INFO=<ID=...>` header lines.

use crate::common::AnnotationRecord;
use anyhow::{Context, Result};
use std::collections::{HashMap, HashSet, VecDeque};
use std::io::BufRead;

/// Population keys to extract from gnomAD VCF INFO field. Covers both
/// gnomAD v2.1 codes (`oth`) and the v4.1 codes (`mid`, `remaining`); a
/// missing key is silently skipped per VCF, so listing all is harmless.
const POPULATIONS: &[&str] = &[
    "afr", "amr", "asj", "eas", "fin", "mid", "nfe", "oth", "remaining", "sas",
];

/// INFO field names for a particular gnomAD release flavor.
///
/// Built from the VCF header so we use whatever names the upstream file
/// actually exposes. See module docs for the v4.1-joint vs. v4.0 split.
#[derive(Debug, Clone)]
struct FieldNames {
    af: String,
    an: String,
    ac: String,
    nhomalt: String,
    /// Format string for per-population AF, with `{}` substituted for the
    /// population code (e.g., `"AF_{}"` or `"AF_joint_{}"`).
    af_pop_template: String,
}

impl FieldNames {
    fn standard() -> Self {
        Self {
            af: "AF".into(),
            an: "AN".into(),
            ac: "AC".into(),
            nhomalt: "nhomalt".into(),
            af_pop_template: "AF_{}".into(),
        }
    }

    fn joint() -> Self {
        Self {
            af: "AF_joint".into(),
            an: "AN_joint".into(),
            ac: "AC_joint".into(),
            nhomalt: "nhomalt_joint".into(),
            af_pop_template: "AF_joint_{}".into(),
        }
    }

    fn pop_key(&self, pop: &str) -> String {
        self.af_pop_template.replace("{}", pop)
    }
}

/// Pick a field-naming scheme based on which INFO IDs the VCF header
/// declares. Prefer standard names when present; fall back to the joint
/// names if the standard `AF` is absent but `AF_joint` is declared.
///
/// Some locally-processed files may use a hyphen separator for per-population
/// fields (`AF-afr`) rather than the standard underscore (`AF_afr`). We
/// detect this by checking whether any `AF-<pop>` field appears in the header.
fn detect_field_names(info_ids: &HashSet<String>) -> FieldNames {
    if info_ids.contains("AF") {
        let mut names = FieldNames::standard();
        if POPULATIONS
            .iter()
            .any(|p| info_ids.contains(&format!("AF-{p}")))
        {
            names.af_pop_template = "AF-{}".into();
        }
        names
    } else if info_ids.contains("AF_joint") {
        FieldNames::joint()
    } else {
        FieldNames::standard()
    }
}

/// Extract the `ID=` value from a `##INFO=<ID=...,Number=...,...>` header
/// line. Returns `None` for non-INFO lines or malformed entries.
///
/// Handles the two real-world quirks of VCF headers:
/// - The trailing `>` when `ID` is the last attribute
///   (`##INFO=<...,ID=AF>` must yield `"AF"`, not `"AF>"`).
/// - Commas inside quoted `Description="foo, bar"` values, which a naive
///   `split(',')` would split on. We walk the body tracking quote state
///   so attributes can appear in any order.
fn parse_info_id(line: &str) -> Option<&str> {
    let body = line.strip_prefix("##INFO=<")?.strip_suffix('>')?;
    let bytes = body.as_bytes();
    let mut i = 0;
    while i < bytes.len() {
        // Skip leading whitespace between attributes.
        while i < bytes.len() && bytes[i] == b' ' {
            i += 1;
        }
        let key_start = i;
        while i < bytes.len() && bytes[i] != b'=' && bytes[i] != b',' {
            i += 1;
        }
        let key = &body[key_start..i];
        // Bare attribute with no value: skip the separator and continue.
        if i >= bytes.len() || bytes[i] == b',' {
            if i < bytes.len() {
                i += 1;
            }
            continue;
        }
        // Consume '=' and read the value, respecting double-quoted strings.
        i += 1;
        let value_start = i;
        let mut in_quotes = false;
        while i < bytes.len() {
            match bytes[i] {
                b'"' => in_quotes = !in_quotes,
                b',' if !in_quotes => break,
                _ => {}
            }
            i += 1;
        }
        if key == "ID" {
            let mut value = &body[value_start..i];
            if value.len() >= 2 && value.starts_with('"') && value.ends_with('"') {
                value = &value[1..value.len() - 1];
            }
            return Some(value);
        }
        if i < bytes.len() {
            i += 1; // skip the ',' separator
        }
    }
    None
}

/// Parse a gnomAD sites-only VCF and produce sorted `AnnotationRecord`s.
///
/// Collects all records into memory and sorts. For genome-scale VCFs use
/// `iter_gnomad_vcf` instead, which streams one block at a time.
pub fn parse_gnomad_vcf<R: BufRead>(
    reader: R,
    chrom_to_idx: &HashMap<String, u16>,
) -> Result<Vec<AnnotationRecord>> {
    let mut records: Vec<_> = iter_gnomad_vcf(reader, chrom_to_idx).collect::<Result<_>>()?;
    records.sort_by(|a, b| a.chrom_idx.cmp(&b.chrom_idx).then(a.position.cmp(&b.position)));
    Ok(records)
}

/// Stream a coordinate-sorted gnomAD VCF as `AnnotationRecord`s without
/// buffering the whole file in memory.
///
/// The input must already be sorted by chromosome and position (all standard
/// gnomAD releases are). The writer will detect and error on out-of-order
/// records.
pub fn iter_gnomad_vcf<'a, R: BufRead>(
    reader: R,
    chrom_to_idx: &'a HashMap<String, u16>,
) -> GnomadRecordIter<'a, R> {
    GnomadRecordIter {
        lines: reader.lines(),
        chrom_to_idx,
        pending: VecDeque::new(),
        info_ids: HashSet::new(),
        field_names: None,
    }
}

pub struct GnomadRecordIter<'a, R: BufRead> {
    lines: std::io::Lines<R>,
    chrom_to_idx: &'a HashMap<String, u16>,
    pending: VecDeque<AnnotationRecord>,
    info_ids: HashSet<String>,
    field_names: Option<FieldNames>,
}

impl<R: BufRead> Iterator for GnomadRecordIter<'_, R> {
    type Item = Result<AnnotationRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            if let Some(record) = self.pending.pop_front() {
                return Some(Ok(record));
            }

            let line = match self.lines.next()? {
                Ok(l) => l,
                Err(e) => return Some(Err(e).context("Reading gnomAD VCF line")),
            };

            if line.starts_with('#') {
                if let Some(id) = parse_info_id(&line) {
                    self.info_ids.insert(id.to_string());
                }
                continue;
            }

            if self.field_names.is_none() {
                self.field_names = Some(detect_field_names(&self.info_ids));
            }
            let field_names = self.field_names.as_ref().unwrap();

            let fields: Vec<&str> = line.splitn(9, '\t').collect();
            if fields.len() < 8 {
                continue;
            }

            let chrom = normalize_chrom(fields[0]);
            let chrom_idx = match self.chrom_to_idx.get(&chrom) {
                Some(&idx) => idx,
                None => continue,
            };

            let pos: u32 = match fields[1].parse() {
                Ok(p) => p,
                Err(_) => continue,
            };

            let ref_allele = fields[3].to_string();
            let alt_field = fields[4];
            let info = fields[7];

            let info_map = parse_info(info);
            let alts: Vec<&str> = alt_field.split(',').collect();
            let all_afs = split_info_values(info_map.get(&field_names.af).map(|s| s.as_str()));
            let all_ans = split_info_values(info_map.get(&field_names.an).map(|s| s.as_str()));
            let all_acs = split_info_values(info_map.get(&field_names.ac).map(|s| s.as_str()));
            let all_nhomalt =
                split_info_values(info_map.get(&field_names.nhomalt).map(|s| s.as_str()));

            for (i, alt) in alts.iter().enumerate() {
                if *alt == "." || *alt == "*" {
                    continue;
                }
                let json = build_gnomad_json(
                    all_afs.get(i).map(|s| s.as_str()),
                    all_ans.first().map(|s| s.as_str()),
                    all_acs.get(i).map(|s| s.as_str()),
                    all_nhomalt.get(i).map(|s| s.as_str()),
                    &info_map,
                    i,
                    field_names,
                );
                self.pending.push_back(AnnotationRecord {
                    chrom_idx,
                    position: pos,
                    ref_allele: ref_allele.clone(),
                    alt_allele: alt.to_string(),
                    json,
                });
            }
        }
    }
}

fn build_gnomad_json(
    af: Option<&str>,
    an: Option<&str>,
    ac: Option<&str>,
    nhomalt: Option<&str>,
    info_map: &HashMap<String, String>,
    allele_idx: usize,
    field_names: &FieldNames,
) -> String {
    let mut parts = Vec::new();

    if let Some(af_str) = af {
        if let Ok(f) = af_str.parse::<f64>() {
            parts.push(format!("\"allAf\":{:.6e}", f));
        }
    }

    if let Some(an_str) = an {
        parts.push(format!("\"allAn\":{}", an_str));
    }

    if let Some(ac_str) = ac {
        parts.push(format!("\"allAc\":{}", ac_str));
    }

    if let Some(nh) = nhomalt {
        parts.push(format!("\"allHc\":{}", nh));
    }

    // Per-population AFs
    for pop in POPULATIONS {
        let key = field_names.pop_key(pop);
        if let Some(val) = info_map.get(&key) {
            let vals = split_info_values(Some(val.as_str()));
            if let Some(af_str) = vals.get(allele_idx) {
                if let Ok(f) = af_str.parse::<f64>() {
                    parts.push(format!("\"{}Af\":{:.6e}", pop, f));
                }
            }
        }
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

fn split_info_values(value: Option<&str>) -> Vec<String> {
    match value {
        Some(v) => v.split(',').map(|s| s.to_string()).collect(),
        None => Vec::new(),
    }
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
    fn test_parse_gnomad_vcf() {
        let vcf = "\
##fileformat=VCFv4.2
##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">
##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total allele number\">
##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t10001\t.\tA\tG\t.\tPASS\tAF=0.001;AN=150000;AC=150;nhomalt=2;AF_afr=0.002;AF_nfe=0.0005
chr1\t20000\t.\tC\tT,A\t.\tPASS\tAF=0.01,0.005;AN=140000;AC=1400,700;nhomalt=10,3;AF_eas=0.02,0.01
";

        let mut chrom_map = HashMap::new();
        chrom_map.insert("chr1".to_string(), 0u16);

        let records = parse_gnomad_vcf(vcf.as_bytes(), &chrom_map).unwrap();
        assert_eq!(records.len(), 3); // 1 SNV + 2 from multi-allelic

        // First record
        assert_eq!(records[0].position, 10001);
        assert!(records[0].json.contains("\"allAf\":"));
        assert!(records[0].json.contains("\"afrAf\":"));
        assert!(records[0].json.contains("\"nfeAf\":"));

        // Multi-allelic: second alt
        assert_eq!(records[2].position, 20000);
        assert_eq!(records[2].alt_allele, "A");
    }

    #[test]
    fn test_parse_gnomad_v41_joint_vcf() {
        // Regression for issue #39: v4.1 joint release uses *_joint suffixes
        // and FV_GNOMAD came out empty because the parser only looked for
        // the bare AF/AC/AN names.
        let vcf = "\
##fileformat=VCFv4.2
##INFO=<ID=AF_joint,Number=A,Type=Float,Description=\"Joint allele frequency\">
##INFO=<ID=AN_joint,Number=1,Type=Integer,Description=\"Joint total allele number\">
##INFO=<ID=AC_joint,Number=A,Type=Integer,Description=\"Joint allele count\">
##INFO=<ID=nhomalt_joint,Number=A,Type=Integer,Description=\"Joint homozygote count\">
##INFO=<ID=AF_joint_afr,Number=A,Type=Float,Description=\"Joint AF AFR\">
##INFO=<ID=AF_joint_nfe,Number=A,Type=Float,Description=\"Joint AF NFE\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t10001\t.\tA\tG\t.\tPASS\tAF_joint=0.001;AN_joint=150000;AC_joint=150;nhomalt_joint=2;AF_joint_afr=0.002;AF_joint_nfe=0.0005
chr1\t20000\t.\tC\tT,A\t.\tPASS\tAF_joint=0.01,0.005;AN_joint=140000;AC_joint=1400,700;nhomalt_joint=10,3;AF_joint_eas=0.02,0.01
";

        let mut chrom_map = HashMap::new();
        chrom_map.insert("chr1".to_string(), 0u16);

        let records = parse_gnomad_vcf(vcf.as_bytes(), &chrom_map).unwrap();
        assert_eq!(records.len(), 3);

        // First record: all frequency fields populated, not empty.
        assert_eq!(records[0].position, 10001);
        assert!(
            records[0].json.contains("\"allAf\":"),
            "missing allAf in: {}",
            records[0].json
        );
        assert!(
            records[0].json.contains("\"allAn\":150000"),
            "missing allAn in: {}",
            records[0].json
        );
        assert!(
            records[0].json.contains("\"allAc\":150"),
            "missing allAc in: {}",
            records[0].json
        );
        assert!(
            records[0].json.contains("\"allHc\":2"),
            "missing allHc in: {}",
            records[0].json
        );
        assert!(
            records[0].json.contains("\"afrAf\":"),
            "missing afrAf in: {}",
            records[0].json
        );
        assert!(
            records[0].json.contains("\"nfeAf\":"),
            "missing nfeAf in: {}",
            records[0].json
        );

        // Multi-allelic: second alt also gets per-allele values.
        assert_eq!(records[2].position, 20000);
        assert_eq!(records[2].alt_allele, "A");
        assert!(
            records[2].json.contains("\"allAc\":700"),
            "second alt should pick second AC value: {}",
            records[2].json
        );
    }

    #[test]
    fn test_detect_field_names_prefers_standard() {
        let mut ids = HashSet::new();
        ids.insert("AF".into());
        ids.insert("AF_joint".into());
        let names = detect_field_names(&ids);
        assert_eq!(names.af, "AF");
    }

    #[test]
    fn test_detect_field_names_hyphen_separator() {
        // Some locally-processed files use AF-afr (hyphen) instead of the standard AF_afr.
        let mut ids = HashSet::new();
        ids.insert("AF".into());
        ids.insert("AN".into());
        ids.insert("AF-afr".into());
        ids.insert("AF-nfe".into());
        let names = detect_field_names(&ids);
        assert_eq!(names.af, "AF");
        assert_eq!(names.pop_key("afr"), "AF-afr");
        assert_eq!(names.pop_key("nfe"), "AF-nfe");
    }

    #[test]
    fn test_detect_field_names_falls_back_to_joint() {
        let mut ids = HashSet::new();
        ids.insert("AF_joint".into());
        ids.insert("AC_joint".into());
        let names = detect_field_names(&ids);
        assert_eq!(names.af, "AF_joint");
        assert_eq!(names.pop_key("nfe"), "AF_joint_nfe");
    }

    #[test]
    fn test_parse_info_id_standard() {
        let line = "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">";
        assert_eq!(parse_info_id(line), Some("AF"));
    }

    #[test]
    fn test_parse_info_id_trailing_angle_bracket() {
        // ID is the last attribute — the closing '>' must not be captured.
        let line = "##INFO=<Number=A,Type=Float,ID=AF>";
        assert_eq!(parse_info_id(line), Some("AF"));
    }

    #[test]
    fn test_parse_info_id_reordered_with_quoted_comma() {
        // Description quoted string contains commas — must not split inside it.
        let line =
            "##INFO=<Number=A,Type=Float,Description=\"AF, joint, multi-pop\",ID=AF_joint>";
        assert_eq!(parse_info_id(line), Some("AF_joint"));
    }

    #[test]
    fn test_parse_info_id_quoted_id_value() {
        let line = "##INFO=<ID=\"weird_id\",Number=A,Type=Float,Description=\"x\">";
        assert_eq!(parse_info_id(line), Some("weird_id"));
    }

    #[test]
    fn test_parse_info_id_non_info_line() {
        assert_eq!(parse_info_id("##fileformat=VCFv4.2"), None);
        assert_eq!(parse_info_id("#CHROM\tPOS\tID"), None);
    }

    #[test]
    fn test_parse_info_id_malformed_no_closing_bracket() {
        // Missing trailing '>' — refuse to parse rather than guess.
        let line = "##INFO=<ID=AF,Number=A";
        assert_eq!(parse_info_id(line), None);
    }
}
