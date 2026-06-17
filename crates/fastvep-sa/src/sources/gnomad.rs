//! gnomAD VCF parser for building .osa annotation files.
//!
//! Parses gnomAD's sites-only VCF to extract allele frequencies
//! per population, allele counts, and homozygote counts.
//!
//! Auto-detects the INFO field-naming convention used by the input. gnomAD
//! v2.1 / v3 / v4.0 (exomes, genomes) use bare `AF` / `AC` / `AN` /
//! `nhomalt` with `AF_<pop>` for population subsets. The v4.1 *joint*
//! release (`gnomad.joint.v4.1.sites.*.vcf.bgz`) instead exposes
//! `AF_joint` / `AC_joint` / `AN_joint` / `nhomalt_joint` and
//! `AF_joint_<pop>` for per-population frequencies. We pick the scheme by
//! scanning `##INFO=<ID=...>` header lines.

use crate::common::AnnotationRecord;
use anyhow::{Context, Result};
use std::collections::{HashMap, HashSet};
use std::io::BufRead;

/// Population keys to extract from gnomAD VCF INFO field. Covers both
/// gnomAD v2.1 codes (`oth`) and the v4.1 codes (`mid`, `remaining`); a
/// missing key is silently skipped per VCF, so listing all is harmless.
const POPULATIONS: &[&str] = &[
    "afr", "ami", "amr", "asj", "eas", "fin", "mid", "nfe", "oth", "remaining", "sas",
];

/// Joint-VCF region flags (valueless INFO flags). No bare equivalent — only the
/// v4.1 joint dataset declares these; emitted as booleans when present.
const FLAGS: &[(&str, &str)] = &[
    ("fail_interval_qc", "failIntervalQc"),
    ("outside_broad_capture_region", "outsideBroadCaptureRegion"),
    ("outside_ukb_capture_region", "outsideUkbCaptureRegion"),
    ("not_called_in_exomes", "notCalledInExomes"),
    ("not_called_in_genomes", "notCalledInGenomes"),
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
    // --- ACMG-relevant fields, per-scheme like the above ---
    /// Group-max AF and the ancestry group that produced it.
    grpmax_af: String,
    grpmax_group: String,
    /// Filtering allele frequency (grpmax-based, Poisson 95% / 99% CI).
    faf95: String,
    faf99: String,
    /// Per-population homozygote count template (`{}` = population code).
    nhomalt_pop_template: String,
    /// Sex-stratified AF.
    xx_af: String,
    xy_af: String,
}

impl FieldNames {
    fn standard() -> Self {
        Self {
            af: "AF".into(),
            an: "AN".into(),
            ac: "AC".into(),
            nhomalt: "nhomalt".into(),
            af_pop_template: "AF_{}".into(),
            grpmax_af: "AF_grpmax".into(),
            grpmax_group: "grpmax".into(),
            faf95: "fafmax_faf95_max".into(),
            faf99: "fafmax_faf99_max".into(),
            nhomalt_pop_template: "nhomalt_{}".into(),
            xx_af: "AF_XX".into(),
            xy_af: "AF_XY".into(),
        }
    }

    fn joint() -> Self {
        Self {
            af: "AF_joint".into(),
            an: "AN_joint".into(),
            ac: "AC_joint".into(),
            nhomalt: "nhomalt_joint".into(),
            af_pop_template: "AF_joint_{}".into(),
            grpmax_af: "AF_grpmax_joint".into(),
            grpmax_group: "grpmax_joint".into(),
            faf95: "fafmax_faf95_max_joint".into(),
            faf99: "fafmax_faf99_max_joint".into(),
            nhomalt_pop_template: "nhomalt_joint_{}".into(),
            xx_af: "AF_joint_XX".into(),
            xy_af: "AF_joint_XY".into(),
        }
    }

    fn pop_key(&self, pop: &str) -> String {
        self.af_pop_template.replace("{}", pop)
    }

    fn nhomalt_pop_key(&self, pop: &str) -> String {
        self.nhomalt_pop_template.replace("{}", pop)
    }
}

/// Pick a field-naming scheme based on which INFO IDs the VCF header
/// declares. Prefer standard names when present; fall back to the joint
/// names if the standard `AF` is absent but `AF_joint` is declared.
fn detect_field_names(info_ids: &HashSet<String>) -> FieldNames {
    if info_ids.contains("AF") {
        FieldNames::standard()
    } else if info_ids.contains("AF_joint") {
        FieldNames::joint()
    } else {
        // Nothing declared — default to standard so behavior is unchanged
        // for malformed or header-less inputs.
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
pub fn parse_gnomad_vcf<R: BufRead>(
    reader: R,
    chrom_to_idx: &HashMap<String, u16>,
) -> Result<Vec<AnnotationRecord>> {
    let mut records = Vec::new();
    let mut info_ids: HashSet<String> = HashSet::new();
    let mut field_names = FieldNames::standard();
    let mut header_done = false;

    for line in reader.lines() {
        let line = line.context("Reading gnomAD VCF line")?;
        if line.starts_with('#') {
            if let Some(id) = parse_info_id(&line) {
                info_ids.insert(id.to_string());
            }
            continue;
        }

        // First data line: finalize the field-naming choice.
        if !header_done {
            field_names = detect_field_names(&info_ids);
            header_done = true;
        }

        let fields: Vec<&str> = line.splitn(9, '\t').collect();
        if fields.len() < 8 {
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

        let ref_allele = fields[3].to_string();
        let alt_field = fields[4];
        let info = fields[7];
        // FILTER column (index 6): recorded verbatim, no gating. PASS / "." -> absent.
        let filter = fields
            .get(6)
            .copied()
            .map(str::trim)
            .filter(|f| !matches!(*f, "PASS" | "." | ""));

        let info_map = parse_info(info);

        // Handle multi-allelic: split allele-specific fields by comma
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
                all_ans.first().map(|s| s.as_str()), // AN is typically single-valued
                all_acs.get(i).map(|s| s.as_str()),
                all_nhomalt.get(i).map(|s| s.as_str()),
                &info_map,
                i,
                &field_names,
                filter,
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

    records.sort_by(|a, b| a.chrom_idx.cmp(&b.chrom_idx).then(a.position.cmp(&b.position)));
    Ok(records)
}

#[allow(clippy::too_many_arguments)]
fn build_gnomad_json(
    af: Option<&str>,
    an: Option<&str>,
    ac: Option<&str>,
    nhomalt: Option<&str>,
    info_map: &HashMap<String, String>,
    allele_idx: usize,
    field_names: &FieldNames,
    filter: Option<&str>,
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

    // --- Group-max AF, filtering AF (faf95/faf99), FILTER, sex AF, region flags. ---
    let getf = |key: &str| -> Option<f64> {
        allele_value(info_map, key, allele_idx).and_then(|v| v.parse::<f64>().ok())
    };
    // Group-max AF + the ancestry group that produced it.
    if let Some(f) = getf(&field_names.grpmax_af) {
        parts.push(format!("\"grpmaxAf\":{:.6e}", f));
    }
    if let Some(g) = allele_value(info_map, &field_names.grpmax_group, allele_idx) {
        parts.push(format!("\"grpmaxGroup\":\"{}\"", escape_json(g)));
    }
    // Filtering allele frequency (grpmax-based; ClinGen-recommended for BA1/BS1).
    if let Some(f) = getf(&field_names.faf95) {
        parts.push(format!("\"faf95\":{:.6e}", f));
    }
    if let Some(f) = getf(&field_names.faf99) {
        parts.push(format!("\"faf99\":{:.6e}", f));
    }
    // FILTER status — recorded only when not PASS (absence == passed). Verbatim,
    // no gating; policy is left to the classifier.
    if let Some(f) = filter {
        parts.push(format!("\"filter\":\"{}\"", escape_json(f)));
    }
    // Sex-stratified AF (JSON-only).
    if let Some(f) = getf(&field_names.xx_af) {
        parts.push(format!("\"xxAf\":{:.6e}", f));
    }
    if let Some(f) = getf(&field_names.xy_af) {
        parts.push(format!("\"xyAf\":{:.6e}", f));
    }
    // Per-population homozygote count (JSON-only, for recessive / population BS2).
    for pop in POPULATIONS {
        if let Some(n) = allele_value(info_map, &field_names.nhomalt_pop_key(pop), allele_idx)
            .and_then(|v| v.parse::<i64>().ok())
        {
            parts.push(format!("\"{}Hc\":{}", pop, n));
        }
    }
    // Region flags (joint-only; emitted true when present).
    for (flag, alias) in FLAGS {
        if info_map.contains_key(*flag) {
            parts.push(format!("\"{}\":true", alias));
        }
    }

    format!("{{{}}}", parts.join(","))
}

/// Look up an allele-indexed INFO value: split by comma, take `allele_idx`
/// (falling back to the first element for single-valued fields such as `AN`),
/// and treat "." / empty as absent.
fn allele_value<'a>(
    info_map: &'a HashMap<String, String>,
    key: &str,
    allele_idx: usize,
) -> Option<&'a str> {
    let raw = info_map.get(key)?;
    let parts: Vec<&str> = raw.split(',').collect();
    let v = parts.get(allele_idx).copied().or_else(|| parts.first().copied())?;
    if v == "." || v.is_empty() {
        None
    } else {
        Some(v)
    }
}

fn escape_json(s: &str) -> String {
    s.replace('\\', "\\\\").replace('"', "\\\"")
}

fn parse_info(info: &str) -> HashMap<String, String> {
    let mut map = HashMap::new();
    for pair in info.split(';') {
        match pair.split_once('=') {
            Some((key, value)) => {
                map.insert(key.to_string(), value.to_string());
            }
            // Valueless INFO flag (e.g. the joint region flags) — record
            // presence with an empty value so `contains_key` works.
            None if !pair.is_empty() => {
                map.insert(pair.to_string(), String::new());
            }
            None => {}
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
    fn test_gnomad_acmg_fields_joint() {
        // ACMG fields on the joint scheme: grpmax,
        // faf95/faf99, verbatim FILTER, per-pop nhomalt, sex AF, region flags.
        // Header declares AF_joint (not AF) so the joint scheme is selected.
        let vcf = "\
##fileformat=VCFv4.2
##INFO=<ID=AF_joint,Number=A,Type=Float,Description=\"joint AF\">
##INFO=<ID=AF_grpmax_joint,Number=A,Type=Float,Description=\"grpmax AF\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\t.\tA\tG\t.\tGENOMES_FILTERED\tAF_joint=1.5e-4;AN_joint=1000000;AC_joint=150;nhomalt_joint=3;AF_joint_afr=2.0e-4;nhomalt_joint_afr=1;AF_grpmax_joint=2.0e-4;grpmax_joint=afr;fafmax_faf95_max_joint=1.0e-4;fafmax_faf99_max_joint=5.0e-5;AF_joint_XX=1.6e-4;not_called_in_exomes
";
        let mut chrom_map = HashMap::new();
        chrom_map.insert("chr1".to_string(), 0u16);

        let records = parse_gnomad_vcf(vcf.as_bytes(), &chrom_map).unwrap();
        assert_eq!(records.len(), 1);
        let j = &records[0].json;

        assert!(j.contains("\"allAf\":1.500000e-4"), "{j}");
        assert!(j.contains("\"allHc\":3"), "{j}");
        assert!(j.contains("\"grpmaxAf\":2.000000e-4"), "{j}");
        assert!(j.contains("\"grpmaxGroup\":\"afr\""), "{j}");
        assert!(j.contains("\"faf95\":1.000000e-4"), "{j}");
        assert!(j.contains("\"faf99\":5.000000e-5"), "{j}");
        assert!(j.contains("\"afrAf\":2.000000e-4"), "{j}");
        assert!(j.contains("\"afrHc\":1"), "{j}");
        assert!(j.contains("\"xxAf\":1.600000e-4"), "{j}");
        assert!(j.contains("\"notCalledInExomes\":true"), "{j}");
        // Filtered site recorded verbatim (no boolean, no gating).
        assert!(j.contains("\"filter\":\"GENOMES_FILTERED\""), "{j}");
        // Per-pop AC/AN are intentionally not emitted (JSON-only AF + nhomalt).
        assert!(!j.contains("\"afrAc\""), "{j}");
        assert!(!j.contains("\"afrAn\""), "{j}");
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
