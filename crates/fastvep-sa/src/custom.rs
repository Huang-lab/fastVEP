//! Custom annotation providers.
//!
//! Allows users to provide their own annotation files (VCF, BED, TSV)
//! that get indexed and queried at annotation time.

use crate::common::AnnotationRecord;
use anyhow::{Context, Result};
use std::collections::{BTreeMap, HashMap};
use std::io::BufRead;

/// Parse a custom VCF annotation file.
///
/// Extracts specified INFO fields as JSON annotations. If `info_fields` is
/// empty, every INFO field present on a record is included.
///
/// Multi-allelic handling: for each ALT, the parser emits its own
/// `AnnotationRecord`, and INFO values that look like per-allele lists
/// (i.e. comma-separated with `n_alts` elements, matching VCF
/// `Number=A`, or `n_alts+1` matching `Number=R`) are split so each ALT
/// gets only its own slice. Values whose comma-count doesn't match are
/// kept whole — this is the conservative thing to do for ad-hoc INFO
/// fields whose `Number` we don't know without parsing the header.
pub fn parse_custom_vcf<R: BufRead>(
    reader: R,
    chrom_to_idx: &HashMap<String, u16>,
    name: &str,
    info_fields: &[String],
) -> Result<Vec<AnnotationRecord>> {
    let mut records = Vec::new();

    for line in reader.lines() {
        let line = line.context("Reading custom VCF")?;
        // Strip CRLF from Windows-produced VCFs. `BufRead::lines` only
        // strips `\n`, so without this the last INFO value carries a
        // trailing `\r` that breaks downstream serde_json parsing.
        let line = line.trim_end_matches('\r');
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

        // Pre-split ALTs once so we can both iterate to emit records and
        // know `n_alts` for per-allele INFO splitting.
        let alts: Vec<&str> = alt_field
            .split(',')
            .filter(|a| *a != "." && *a != "*")
            .collect();
        if alts.is_empty() {
            continue;
        }
        let n_alts = alts.len();

        for (alt_idx, alt) in alts.iter().enumerate() {
            // Build the JSON object for *this specific ALT* using serde_json
            // so escaping is correct for control chars, tabs, embedded
            // quotes, etc. Per-allele INFO arrays (Number=A / Number=R)
            // get the right slice; everything else is shared verbatim
            // across alts.
            let mut obj = serde_json::Map::new();
            let mut push = |key: &str, val: &str| {
                obj.insert(key.to_string(), serde_json::Value::String(val.to_string()));
            };
            if info_fields.is_empty() {
                for (key, val) in &info_map {
                    let v = pick_per_allele(val, alt_idx, n_alts);
                    push(key, v);
                }
            } else {
                for field in info_fields {
                    if let Some(val) = info_map.get(field.as_str()) {
                        let v = pick_per_allele(val, alt_idx, n_alts);
                        push(field, v);
                    }
                }
            }
            if obj.is_empty() {
                obj.insert(
                    "source".to_string(),
                    serde_json::Value::String(name.to_string()),
                );
            }
            // `serde_json::to_string` on a `Map<String, Value>` produces
            // a well-formed JSON object string.
            let json = serde_json::to_string(&obj)
                .unwrap_or_else(|_| "{}".to_string());

            records.push(AnnotationRecord {
                chrom_idx,
                position: pos,
                ref_allele: ref_allele.clone(),
                alt_allele: (*alt).to_string(),
                json,
            });
        }
    }

    records.sort_by(|a, b| a.chrom_idx.cmp(&b.chrom_idx).then(a.position.cmp(&b.position)));
    Ok(records)
}

/// Pick the per-allele value out of a possibly-arrayed INFO field.
///
/// Auto-detection only fires for multi-allelic records (`n_alts > 1`) —
/// for a bi-allelic record (`n_alts == 1`) a 2-element value like
/// `"AC_male,AC_female"` is indistinguishable from a Number=R array and
/// silently mis-splitting it would corrupt categorical fields. Bi-allelic
/// custom VCFs are the typical user input (the standard workflow runs
/// `bcftools norm -m -any` upstream of sa-build); they keep the value
/// verbatim, which is the safe behaviour.
///
/// For `n_alts > 1`:
/// - If `val` has exactly `n_alts` comma-separated elements, treat as Number=A.
/// - If `val` has exactly `n_alts + 1` elements, treat as Number=R (skip REF).
/// - Otherwise, return `val` unchanged (shared across all alts).
fn pick_per_allele(val: &str, alt_idx: usize, n_alts: usize) -> &str {
    // Fast path: no comma → can't be a per-allele list.
    if !val.contains(',') {
        return val;
    }
    // Single-ALT records are too ambiguous to auto-split safely; see the
    // doc above.
    if n_alts < 2 {
        return val;
    }
    let pieces: Vec<&str> = val.split(',').collect();
    if pieces.len() == n_alts {
        return pieces[alt_idx];
    }
    if pieces.len() == n_alts + 1 {
        // Number=R: first slot is REF.
        return pieces[alt_idx + 1];
    }
    val
}

/// Parse a custom BED annotation file into interval records.
///
/// Format: chrom, start (0-based), end, name, [score], [additional fields]
pub fn parse_custom_bed<R: BufRead>(
    reader: R,
    chrom_to_idx: &HashMap<String, u16>,
) -> Result<Vec<crate::common::IntervalRecord>> {
    let mut records = Vec::new();

    for line in reader.lines() {
        let line = line.context("Reading custom BED")?;
        // Strip trailing CR (Windows CRLF) — without this the `end` field
        // parse fails for every record on a CRLF-terminated file and the
        // whole BED silently produces 0 intervals.
        let line = line.trim_end_matches('\r');
        if line.starts_with('#') || line.starts_with("track") || line.is_empty() {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 3 {
            continue;
        }

        let chrom = normalize_chrom(fields[0]);
        if !chrom_to_idx.contains_key(&chrom) {
            continue;
        }

        let start_bed: u32 = match fields[1].parse::<u32>() {
            Ok(s) => s,
            Err(_) => continue,
        };
        let end: u32 = match fields[2].parse::<u32>() {
            Ok(e) => e,
            Err(_) => continue,
        };
        // BED is 0-based half-open. Convert to fastVEP's 1-based closed.
        // Guard against `start_bed == u32::MAX` (saturating add) and
        // against malformed `end <= start` (zero-/negative-width intervals
        // that would otherwise survive into the index and silently miss
        // every query).
        let start = start_bed.saturating_add(1);
        if end < start {
            continue;
        }

        let name = fields.get(3).unwrap_or(&".").to_string();
        let score = fields.get(4).unwrap_or(&".").to_string();

        // Build the JSON via serde_json so embedded control chars / tabs
        // / quotes in the name field don't break downstream consumers.
        let mut obj = serde_json::Map::new();
        if name != "." {
            obj.insert("name".to_string(), serde_json::Value::String(name.clone()));
        }
        if score != "." {
            if let Ok(s) = score.parse::<f64>() {
                if let Some(n) = serde_json::Number::from_f64(s) {
                    obj.insert("score".to_string(), serde_json::Value::Number(n));
                }
            }
        }
        let json = serde_json::to_string(&obj).unwrap_or_else(|_| "{}".to_string());

        records.push(crate::common::IntervalRecord {
            chrom,
            start,
            end,
            json,
        });
    }

    Ok(records)
}

/// Parse a VCF INFO field into a deterministic `(key → value)` map.
///
/// `BTreeMap` (not `HashMap`) so iteration order is stable across runs —
/// this directly affects the byte layout of the resulting `.osa` and is
/// what makes the build content-hash-reproducible.
fn parse_info(info: &str) -> BTreeMap<String, String> {
    let mut map = BTreeMap::new();
    for pair in info.split(';') {
        let pair = pair.trim();
        if pair.is_empty() {
            continue;
        }
        match pair.split_once('=') {
            Some((k, v)) => {
                map.insert(k.to_string(), v.to_string());
            }
            // Flag-style INFO entry (e.g. `SOMATIC`, `H3`). VCF spec calls
            // these `Number=0,Type=Flag`. Store them as the JSON true-ish
            // string so user filters that look up the key by name find a
            // truthy value rather than missing it entirely.
            None => {
                map.insert(pair.to_string(), "true".to_string());
            }
        }
    }
    map
}

fn normalize_chrom(chrom: &str) -> String {
    if chrom.starts_with("chr") { chrom.to_string() } else { format!("chr{}", chrom) }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_custom_vcf() {
        let vcf = "#h\nchr1\t100\t.\tA\tG\t.\t.\tMY_SCORE=0.95;MY_FLAG=true\n";
        let mut m = HashMap::new();
        m.insert("chr1".into(), 0u16);

        // All fields
        let recs = parse_custom_vcf(vcf.as_bytes(), &m, "test", &[]).unwrap();
        assert_eq!(recs.len(), 1);
        assert!(recs[0].json.contains("MY_SCORE"));

        // Specific fields
        let recs = parse_custom_vcf(
            vcf.as_bytes(), &m, "test",
            &["MY_SCORE".to_string()],
        ).unwrap();
        assert!(recs[0].json.contains("MY_SCORE"));
        assert!(!recs[0].json.contains("MY_FLAG"));
    }

    #[test]
    fn test_parse_custom_bed() {
        let bed = "chr1\t99\t200\tregion1\t0.5\nchr1\t499\t600\tregion2\n";
        let mut m = HashMap::new();
        m.insert("chr1".into(), 0u16);
        let recs = parse_custom_bed(bed.as_bytes(), &m).unwrap();
        assert_eq!(recs.len(), 2);
        assert_eq!(recs[0].start, 100); // 0-based -> 1-based
        assert_eq!(recs[0].end, 200);
        assert!(recs[0].json.contains("region1"));
        assert!(recs[0].json.contains("0.5"));
    }

    #[test]
    fn test_parse_custom_vcf_multiallelic_splits_per_allele_info() {
        // AF is `Number=A` (one per ALT). Without splitting, the user-
        // facing JSON would say AF=0.1,0.9 for both alts, which is wrong.
        // FLAG (with no `=`) is a Number=0 / Type=Flag-style entry — we
        // store it as `"true"` so name-based lookups don't miss it.
        let vcf = "#h\nchr1\t100\t.\tA\tG,T\t.\t.\tAF=0.1,0.9;DP=50;FLAG\n";
        let mut m = HashMap::new();
        m.insert("chr1".into(), 0u16);
        let recs = parse_custom_vcf(vcf.as_bytes(), &m, "test", &[]).unwrap();
        assert_eq!(recs.len(), 2);
        // Records are sorted by chrom+pos but the per-ALT emission order
        // within a position is preserved, so [0]=G, [1]=T.
        let g_rec = recs.iter().find(|r| r.alt_allele == "G").unwrap();
        let t_rec = recs.iter().find(|r| r.alt_allele == "T").unwrap();
        assert!(g_rec.json.contains(r#""AF":"0.1""#), "{}", g_rec.json);
        assert!(t_rec.json.contains(r#""AF":"0.9""#), "{}", t_rec.json);
        // DP (Number=1) is shared across ALTs unchanged.
        assert!(g_rec.json.contains(r#""DP":"50""#));
        assert!(t_rec.json.contains(r#""DP":"50""#));
        // Flag-only INFO field stored as `true`.
        assert!(g_rec.json.contains(r#""FLAG":"true""#));
    }

    #[test]
    fn test_parse_custom_vcf_number_R_per_allele() {
        // AD is `Number=R` (REF + each ALT). With 2 ALTs we have 3 values
        // and each ALT should pick its own slot (skipping the REF slot).
        let vcf = "#h\nchr1\t100\t.\tA\tG,T\t.\t.\tAD=80,12,8\n";
        let mut m = HashMap::new();
        m.insert("chr1".into(), 0u16);
        let recs = parse_custom_vcf(vcf.as_bytes(), &m, "test", &[]).unwrap();
        let g = recs.iter().find(|r| r.alt_allele == "G").unwrap();
        let t = recs.iter().find(|r| r.alt_allele == "T").unwrap();
        assert!(g.json.contains(r#""AD":"12""#), "{}", g.json);
        assert!(t.json.contains(r#""AD":"8""#), "{}", t.json);
    }

    #[test]
    fn test_parse_custom_vcf_handles_crlf_line_endings() {
        // BufRead::lines strips only `\n`. Without an explicit CRLF trim,
        // the last INFO value carries a trailing `\r` and the resulting
        // JSON contains a literal CR, which fails downstream serde_json
        // parsing. Verify the round-trip JSON is parseable.
        let vcf = "##h\r\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\r\nchr1\t100\t.\tA\tG\t.\t.\tCLIN=hot\r\n";
        let mut m = HashMap::new();
        m.insert("chr1".into(), 0u16);
        let recs = parse_custom_vcf(vcf.as_bytes(), &m, "t", &[]).unwrap();
        assert_eq!(recs.len(), 1);
        let val: serde_json::Value = serde_json::from_str(&recs[0].json)
            .expect("JSON must be parseable; CRLF likely leaked into the value");
        assert_eq!(val["CLIN"], "hot");
    }

    #[test]
    fn test_parse_custom_vcf_escapes_quotes_and_backslashes() {
        // INFO values containing literal `"` or `\` — the old hand-rolled
        // `escape_json` got these right too, but exercising the path
        // through serde_json confirms we haven't regressed.
        let vcf = "#h\nchr1\t100\t.\tA\tG\t.\t.\tDESC=has\"quote\\and\\\\bs\n";
        let mut m = HashMap::new();
        m.insert("chr1".into(), 0u16);
        let recs = parse_custom_vcf(vcf.as_bytes(), &m, "t", &[]).unwrap();
        assert_eq!(recs.len(), 1);
        // The serialised JSON must round-trip through serde_json. The old
        // hand-rolled escaper produced output that *some* strict parsers
        // rejected; we now use serde_json end-to-end.
        let val: serde_json::Value = serde_json::from_str(&recs[0].json).expect("parseable JSON");
        assert_eq!(val["DESC"], "has\"quote\\and\\\\bs");
    }

    #[test]
    fn test_parse_custom_vcf_info_order_is_deterministic() {
        // BTreeMap iteration yields sorted keys, so the byte layout of
        // the resulting JSON is stable across runs and identical inputs
        // produce identical .osa contents (content-hash reproducible).
        let vcf = "#h\nchr1\t100\t.\tA\tG\t.\t.\tZED=1;APPLE=2;MID=3\n";
        let mut m = HashMap::new();
        m.insert("chr1".into(), 0u16);
        let r1 = parse_custom_vcf(vcf.as_bytes(), &m, "t", &[]).unwrap();
        let r2 = parse_custom_vcf(vcf.as_bytes(), &m, "t", &[]).unwrap();
        assert_eq!(r1[0].json, r2[0].json);
        // Sorted: APPLE, MID, ZED.
        assert!(r1[0].json.find("APPLE").unwrap() < r1[0].json.find("MID").unwrap());
        assert!(r1[0].json.find("MID").unwrap() < r1[0].json.find("ZED").unwrap());
    }

    #[test]
    fn test_pick_per_allele_skips_split_for_biallelic_records() {
        // Bi-allelic record with a 2-value INFO that *looks* like Number=R
        // could be a genuine Number=2 categorical — the older code would
        // mis-split. With n_alts < 2 we now keep the value whole.
        let vcf = "#h\nchr1\t100\t.\tA\tG\t.\t.\tCAT=foo,bar\n";
        let mut m = HashMap::new();
        m.insert("chr1".into(), 0u16);
        let recs = parse_custom_vcf(vcf.as_bytes(), &m, "t", &[]).unwrap();
        assert_eq!(recs.len(), 1);
        let val: serde_json::Value = serde_json::from_str(&recs[0].json).unwrap();
        assert_eq!(val["CAT"], "foo,bar", "biallelic 2-value field must stay intact");
    }

    #[test]
    fn test_parse_custom_bed_handles_crlf() {
        let bed = "chr1\t99\t200\tregion1\t0.5\r\nchr1\t499\t600\tregion2\r\n";
        let mut m = HashMap::new();
        m.insert("chr1".into(), 0u16);
        let recs = parse_custom_bed(bed.as_bytes(), &m).unwrap();
        assert_eq!(recs.len(), 2, "CRLF must not break end-field parsing");
        let val: serde_json::Value = serde_json::from_str(&recs[0].json).unwrap();
        assert_eq!(val["name"], "region1");
    }

    #[test]
    fn test_parse_custom_bed_handles_pathological_inputs() {
        let mut m = HashMap::new();
        m.insert("chr1".into(), 0u16);

        // start == u32::MAX would overflow the +1 conversion; we use
        // saturating_add so the resulting interval is still emitted (with
        // start == u32::MAX) instead of panicking in debug builds.
        let near_max = format!("chr1\t{}\t{}\tx\n", u32::MAX, u32::MAX);
        let recs = parse_custom_bed(near_max.as_bytes(), &m).unwrap();
        assert_eq!(recs.len(), 1);
        assert_eq!(recs[0].start, u32::MAX);

        // end < start (after 0→1 conversion) is malformed BED — must be
        // skipped, not silently stored as a phantom interval.
        let bad = "chr1\t100\t50\trev\n";
        let recs = parse_custom_bed(bad.as_bytes(), &m).unwrap();
        assert!(recs.is_empty());

        // Empty-after-comments file is valid; just produces zero records.
        let comments = "# header\ntrack name=x\n# more\n";
        let recs = parse_custom_bed(comments.as_bytes(), &m).unwrap();
        assert!(recs.is_empty());
    }
}
