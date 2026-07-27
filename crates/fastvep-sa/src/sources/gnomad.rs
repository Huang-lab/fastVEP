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
use crate::fields::{Field, FieldType};
use crate::writer_v2::{Osa2Metadata, Osa2Record};
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

    // AN/AC/nhomalt are written unquoted, so each must be a validated
    // non-negative integer — raw INFO-field text (e.g. the "." missing-value
    // sentinel or other garbage) would otherwise land in the JSON as a bare,
    // unquoted token and break every downstream serde_json::from_str on this
    // record. Mirrors the CNT validation in sources/cosmic.rs.
    if let Some(an_str) = an {
        if let Ok(n) = an_str.parse::<u64>() {
            parts.push(format!("\"allAn\":{}", n));
        }
    }

    if let Some(ac_str) = ac {
        if let Ok(n) = ac_str.parse::<u64>() {
            parts.push(format!("\"allAc\":{}", n));
        }
    }

    if let Some(nh) = nhomalt {
        if let Ok(n) = nh.parse::<u64>() {
            parts.push(format!("\"allHc\":{}", n));
        }
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

// =============================================================================
// v2 (.osa2) encoding
// =============================================================================

/// Multiplier for allele-frequency floats stored as u32 (`stored = af * M`).
/// 2e6 keeps ~6-7 significant digits across the 0..1 AF range while staying
/// well under `u32::MAX`.
const AF_MULTIPLIER: u32 = 2_000_000;

fn af_field(field: &str, alias: &str, description: &str) -> Field {
    Field {
        field: field.into(),
        alias: alias.into(),
        ftype: FieldType::Float,
        multiplier: AF_MULTIPLIER,
        zigzag: false,
        missing_value: u32::MAX,
        missing_string: ".".into(),
        description: description.into(),
    }
}

fn count_field(field: &str, alias: &str, description: &str) -> Field {
    Field {
        field: field.into(),
        alias: alias.into(),
        ftype: FieldType::Integer,
        multiplier: 1,
        zigzag: false,
        missing_value: u32::MAX,
        missing_string: ".".into(),
        description: description.into(),
    }
}

/// Canonical gnomAD v2 (`.osa2`) field schema. The value vector produced by
/// [`iter_gnomad_osa2`] is parallel to this list, in this exact order: global
/// AF / AN / AC / nhomalt followed by per-population AF for each entry in
/// [`POPULATIONS`]. The aliases match the JSON keys the v1 builder emits
/// (`allAf`, `allAn`, `allAc`, `allHc`, `<pop>Af`) so downstream consumers see
/// identical output regardless of which format backs the database.
pub fn gnomad_osa2_fields() -> Vec<Field> {
    let mut fields = vec![
        af_field("AF", "allAf", "Global allele frequency"),
        count_field("AN", "allAn", "Total allele number"),
        count_field("AC", "allAc", "Allele count"),
        count_field("nhomalt", "allHc", "Homozygote count"),
    ];
    for pop in POPULATIONS {
        fields.push(af_field(
            &format!("AF_{pop}"),
            &format!("{pop}Af"),
            &format!("{} allele frequency", pop.to_uppercase()),
        ));
    }
    fields
}

/// Standard gnomAD `.osa2` metadata. Mirrors the v1 header (`json_key =
/// "gnomad"`, allele-matched, non-positional) so the annotation pipeline
/// treats a v1 and v2 gnomAD database identically.
pub fn gnomad_osa2_metadata(assembly: &str) -> Osa2Metadata {
    Osa2Metadata {
        format_version: 2,
        name: "gnomAD".into(),
        version: "latest".into(),
        assembly: assembly.into(),
        json_key: "gnomad".into(),
        match_by_allele: true,
        is_array: false,
        is_positional: false,
        chunk_bits: 20,
        description: format!("gnomAD population frequencies for {assembly}"),
    }
}

/// Stream a coordinate-sorted gnomAD VCF as `.osa2` records with parallel u32
/// value arrays (parallel to [`gnomad_osa2_fields`]), without buffering the
/// whole file. Reuses the same INFO field-name detection as the v1 path, so
/// v2/v3/v4/joint releases and hyphen-separated local files are all handled.
pub fn iter_gnomad_osa2<'a, R: BufRead>(
    reader: R,
    chrom_to_idx: &'a HashMap<String, u16>,
) -> GnomadOsa2Iter<'a, R> {
    GnomadOsa2Iter {
        lines: reader.lines(),
        chrom_to_idx,
        fields: gnomad_osa2_fields(),
        pending: VecDeque::new(),
        info_ids: HashSet::new(),
        field_names: None,
    }
}

pub struct GnomadOsa2Iter<'a, R: BufRead> {
    lines: std::io::Lines<R>,
    chrom_to_idx: &'a HashMap<String, u16>,
    fields: Vec<Field>,
    pending: VecDeque<Osa2Record>,
    info_ids: HashSet<String>,
    field_names: Option<FieldNames>,
}

impl<R: BufRead> Iterator for GnomadOsa2Iter<'_, R> {
    type Item = Result<Osa2Record>;

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

            let cols: Vec<&str> = line.splitn(9, '\t').collect();
            if cols.len() < 8 {
                continue;
            }

            let chrom = normalize_chrom(cols[0]);
            if !self.chrom_to_idx.contains_key(&chrom) {
                continue;
            }
            let pos: u32 = match cols[1].parse() {
                Ok(p) => p,
                Err(_) => continue,
            };
            let ref_allele = cols[3].as_bytes().to_vec();
            let alt_field = cols[4];
            let info_map = parse_info(cols[7]);

            let all_afs = split_info_values(info_map.get(&field_names.af).map(|s| s.as_str()));
            // AN is a per-site (Number=1) value; AF/AC/nhomalt are per-allele.
            let an = all_from(&info_map, &field_names.an);
            let all_acs = split_info_values(info_map.get(&field_names.ac).map(|s| s.as_str()));
            let all_nh = split_info_values(info_map.get(&field_names.nhomalt).map(|s| s.as_str()));

            for (ai, alt) in alt_field.split(',').enumerate() {
                if alt == "." || alt == "*" {
                    continue;
                }
                // Value order MUST match `gnomad_osa2_fields()`.
                let mut values = Vec::with_capacity(self.fields.len());
                values.push(enc_float(&self.fields[0], all_afs.get(ai).map(|s| s.as_str())));
                values.push(enc_int(&self.fields[1], an.first().map(|s| s.as_str())));
                values.push(enc_int(&self.fields[2], all_acs.get(ai).map(|s| s.as_str())));
                values.push(enc_int(&self.fields[3], all_nh.get(ai).map(|s| s.as_str())));
                for (pi, pop) in POPULATIONS.iter().enumerate() {
                    let key = field_names.pop_key(pop);
                    let vals = split_info_values(info_map.get(&key).map(|s| s.as_str()));
                    values.push(enc_float(&self.fields[4 + pi], vals.get(ai).map(|s| s.as_str())));
                }

                self.pending.push_back(Osa2Record {
                    chrom: chrom.clone(),
                    position: pos,
                    ref_allele: ref_allele.clone(),
                    alt_allele: alt.as_bytes().to_vec(),
                    values,
                    json_blob: None,
                });
            }
        }
    }
}

fn all_from(info: &HashMap<String, String>, key: &str) -> Vec<String> {
    split_info_values(info.get(key).map(|s| s.as_str()))
}

fn enc_float(field: &Field, raw: Option<&str>) -> u32 {
    match raw.and_then(|s| s.parse::<f64>().ok()) {
        Some(f) => field.encode_float(f),
        None => field.missing_value,
    }
}

fn enc_int(field: &Field, raw: Option<&str>) -> u32 {
    match raw.and_then(|s| s.parse::<i64>().ok()) {
        Some(i) => field.encode_int(i),
        None => field.missing_value,
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
    fn test_garbage_an_ac_nhomalt_omitted_not_emitted_unescaped() {
        // AN/AC/nhomalt are spliced into the JSON unquoted, so malformed
        // upstream text (e.g. the "." missing-value sentinel or other
        // garbage) must be dropped rather than emitted as a bare token that
        // would break serde_json::from_str on the whole record.
        let vcf = "\
##fileformat=VCFv4.2
##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">
##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total allele number\">
##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t10001\t.\tA\tG\t.\tPASS\tAF=0.001;AN=not_a_number;AC=garbage;nhomalt=.
";
        let mut chrom_map = HashMap::new();
        chrom_map.insert("chr1".to_string(), 0u16);

        let records = parse_gnomad_vcf(vcf.as_bytes(), &chrom_map).unwrap();
        assert_eq!(records.len(), 1);
        // Malformed fields must be entirely absent, not emitted as garbage.
        assert!(!records[0].json.contains("allAn"));
        assert!(!records[0].json.contains("allAc"));
        assert!(!records[0].json.contains("allHc"));
        // The emitted JSON must still be valid.
        let v: serde_json::Value = serde_json::from_str(&records[0].json).unwrap();
        assert!(v.get("allAf").is_some());
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

    // ---- v2 (.osa2) encoder ----

    fn chr1_map() -> HashMap<String, u16> {
        let mut m = HashMap::new();
        m.insert("chr1".to_string(), 0u16);
        m
    }

    #[test]
    fn test_osa2_fields_order_matches_value_layout() {
        // The iterator indexes `fields[0..4]` for the global stats and
        // `fields[4 + pi]` per population, so this ordering is load-bearing.
        let fields = gnomad_osa2_fields();
        assert_eq!(fields.len(), 4 + POPULATIONS.len());
        assert_eq!(fields[0].alias, "allAf");
        assert_eq!(fields[1].alias, "allAn");
        assert_eq!(fields[2].alias, "allAc");
        assert_eq!(fields[3].alias, "allHc");
        for (pi, pop) in POPULATIONS.iter().enumerate() {
            assert_eq!(fields[4 + pi].alias, format!("{pop}Af"));
        }
    }

    #[test]
    fn test_iter_gnomad_osa2_encodes_values() {
        let vcf = "\
##fileformat=VCFv4.2
##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">
##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total allele number\">
##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t10001\t.\tA\tG\t.\tPASS\tAF=0.001;AN=150000;AC=150;nhomalt=2;AF_afr=0.002;AF_nfe=0.0005
chr1\t20000\t.\tC\tT,A\t.\tPASS\tAF=0.01,0.005;AN=140000;AC=1400,700;nhomalt=10,3;AF_eas=0.02,0.01
";
        let map = chr1_map();
        let recs: Vec<Osa2Record> =
            iter_gnomad_osa2(vcf.as_bytes(), &map).collect::<Result<_>>().unwrap();
        assert_eq!(recs.len(), 3); // 1 SNV + 2 from the multi-allelic site

        let fields = gnomad_osa2_fields();
        // First record: AF=0.001 -> 0.001 * 2_000_000 = 2000.
        assert_eq!(recs[0].position, 10001);
        assert_eq!(recs[0].values[0], 2000);
        assert_eq!(recs[0].values[1], 150000); // AN
        assert_eq!(recs[0].values[2], 150); // AC
        assert_eq!(recs[0].values[3], 2); // nhomalt
        // afr AF present, sas AF missing.
        let afr_idx = 4 + POPULATIONS.iter().position(|p| *p == "afr").unwrap();
        let sas_idx = 4 + POPULATIONS.iter().position(|p| *p == "sas").unwrap();
        assert_eq!(recs[0].values[afr_idx], fields[afr_idx].encode_float(0.002));
        assert_eq!(recs[0].values[sas_idx], u32::MAX); // missing sentinel

        // Multi-allelic second alt picks the second per-allele values.
        assert_eq!(recs[2].position, 20000);
        assert_eq!(recs[2].alt_allele, b"A");
        assert_eq!(recs[2].values[2], 700); // AC second alt
    }

    #[test]
    fn test_iter_gnomad_osa2_handles_joint_release() {
        // v4.1 joint naming must be detected here too (regression parity with
        // the v1 parser).
        let vcf = "\
##fileformat=VCFv4.2
##INFO=<ID=AF_joint,Number=A,Type=Float,Description=\"Joint AF\">
##INFO=<ID=AN_joint,Number=1,Type=Integer,Description=\"Joint AN\">
##INFO=<ID=AC_joint,Number=A,Type=Integer,Description=\"Joint AC\">
##INFO=<ID=AF_joint_nfe,Number=A,Type=Float,Description=\"Joint AF NFE\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t10001\t.\tA\tG\t.\tPASS\tAF_joint=0.001;AN_joint=150000;AC_joint=150;AF_joint_nfe=0.0005
";
        let map = chr1_map();
        let recs: Vec<Osa2Record> =
            iter_gnomad_osa2(vcf.as_bytes(), &map).collect::<Result<_>>().unwrap();
        assert_eq!(recs.len(), 1);
        assert_eq!(recs[0].values[0], 2000); // AF_joint 0.001 encoded
        assert_eq!(recs[0].values[1], 150000); // AN_joint
        let nfe_idx = 4 + POPULATIONS.iter().position(|p| *p == "nfe").unwrap();
        assert_ne!(recs[0].values[nfe_idx], u32::MAX); // nfe populated from joint
    }

    #[test]
    fn test_iter_gnomad_osa2_skips_unknown_contigs() {
        let vcf = "\
##fileformat=VCFv4.2
##INFO=<ID=AF,Number=A,Type=Float,Description=\"AF\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr2\t100\t.\tA\tG\t.\tPASS\tAF=0.5
chr1\t200\t.\tA\tG\t.\tPASS\tAF=0.5
";
        let map = chr1_map(); // only chr1 is valid
        let recs: Vec<Osa2Record> =
            iter_gnomad_osa2(vcf.as_bytes(), &map).collect::<Result<_>>().unwrap();
        assert_eq!(recs.len(), 1);
        assert_eq!(recs[0].chrom, "chr1");
        assert_eq!(recs[0].position, 200);
    }
}
