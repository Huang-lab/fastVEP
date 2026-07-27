//! AlphaMissense pathogenicity parser.
//!
//! AlphaMissense (Cheng et al., Science 2023; Zenodo record 8208688) predicts
//! the pathogenicity of every possible single-amino-acid missense variant. The
//! genome-coordinate releases (`AlphaMissense_hg38.tsv.gz`,
//! `AlphaMissense_hg19.tsv.gz`) are coordinate-sorted, tab-separated files with
//! a leading comment/license block and a header line:
//!
//! ```text
//! #CHROM  POS  REF  ALT  genome  uniprot_id  transcript_id  protein_variant  am_pathogenicity  am_class
//! chr1    69094   G    T    hg38    Q8NH21   ENST00000335137.4  V2L   0.2937  likely_benign
//! ```
//!
//! Each row carries an allele-specific pathogenicity score in `[0, 1]` and a
//! three-level class (`likely_benign` / `ambiguous` / `likely_pathogenic`).
//! The canonical `AlphaMissense_hg38.tsv.gz` holds one row per genomic variant
//! (canonical transcript), so it maps cleanly onto allele-matched supplementary
//! annotation.
//!
//! This is a numeric-plus-small-categorical payload — a sweet spot for the v2
//! `.osa2` format (one u32 score column plus a u32 class-index column against a
//! 3-entry string table), so AlphaMissense builds natively into v2 rather than
//! going through the generic v1→v2 JSON bridge (which cannot encode
//! categoricals). The v1 `.osa` path builds each record's JSON through the very
//! same `Field`/`format_value` code the v2 reader reconstructs with, so the two
//! formats emit byte-identical annotations.

use crate::common::AnnotationRecord;
use crate::fields::{format_value, Field, FieldType};
use crate::writer_v2::{Osa2Metadata, Osa2Record};
use anyhow::{Context, Result};
use std::collections::HashMap;
use std::io::BufRead;

/// The three AlphaMissense classes, in the order that defines their stored
/// u32 index (0 = `likely_benign`, 1 = `ambiguous`, 2 = `likely_pathogenic`).
/// This is the categorical string table written into the `.osa2` archive.
pub const AM_CLASSES: &[&str] = &["likely_benign", "ambiguous", "likely_pathogenic"];

/// Field index of the pathogenicity score within [`alphamissense_osa2_fields`].
const SCORE_FIELD: usize = 0;
/// Field index of the class within [`alphamissense_osa2_fields`].
const CLASS_FIELD: usize = 1;

/// Canonical AlphaMissense `.osa2` field schema: a float pathogenicity score
/// plus a categorical class. Aliases are the JSON keys the annotation output
/// carries (`amPathogenicity`, `amClass`).
///
/// The score uses a 1e6 multiplier — AlphaMissense distributes four-decimal
/// scores, so six decimals of storage precision are lossless.
pub fn alphamissense_osa2_fields() -> Vec<Field> {
    vec![
        Field {
            field: "am_pathogenicity".into(),
            alias: "amPathogenicity".into(),
            ftype: FieldType::Float,
            multiplier: 1_000_000,
            zigzag: false,
            missing_value: u32::MAX,
            missing_string: ".".into(),
            description: "AlphaMissense pathogenicity score".into(),
        },
        Field {
            field: "am_class".into(),
            alias: "amClass".into(),
            ftype: FieldType::Categorical,
            multiplier: 1,
            zigzag: false,
            missing_value: u32::MAX,
            missing_string: ".".into(),
            description: "AlphaMissense class (likely_benign/ambiguous/likely_pathogenic)".into(),
        },
    ]
}

/// The categorical string table for the `amClass` field, paired with its field
/// index. Threaded into the v2 streaming writer via `set_string_table`.
pub fn alphamissense_string_tables() -> Vec<(usize, Vec<String>)> {
    vec![(
        CLASS_FIELD,
        AM_CLASSES.iter().map(|s| s.to_string()).collect(),
    )]
}

/// Standard AlphaMissense `.osa2` metadata (`json_key = "alphaMissense"`,
/// allele-matched, non-positional).
pub fn alphamissense_osa2_metadata(assembly: &str) -> Osa2Metadata {
    Osa2Metadata {
        format_version: 2,
        name: "AlphaMissense".into(),
        version: "latest".into(),
        assembly: assembly.into(),
        json_key: "alphaMissense".into(),
        match_by_allele: true,
        is_array: false,
        is_positional: false,
        chunk_bits: 20,
        description: format!("AlphaMissense pathogenicity predictions for {assembly}"),
    }
}

/// One parsed AlphaMissense row, reduced to the fields we store.
struct AmRow {
    chrom_idx: u16,
    position: u32,
    ref_allele: String,
    alt_allele: String,
    score: f64,
    /// Index into [`AM_CLASSES`], or `None` when the class is absent/unknown.
    class_idx: Option<u32>,
}

/// Parse a single AlphaMissense data line. Returns `None` for comment/header
/// lines, records whose chromosome is not in `chrom_to_idx`, and rows that are
/// malformed (too few columns, unparseable position/score).
fn parse_row(line: &str, chrom_to_idx: &HashMap<String, u16>) -> Option<AmRow> {
    if line.is_empty() || line.starts_with('#') {
        return None;
    }
    let fields: Vec<&str> = line.split('\t').collect();
    // CHROM POS REF ALT genome uniprot_id transcript_id protein_variant
    // am_pathogenicity am_class  => 10 columns.
    if fields.len() < 10 {
        return None;
    }

    let chrom = normalize_chrom(fields[0]);
    let chrom_idx = *chrom_to_idx.get(&chrom)?;
    let position: u32 = fields[1].parse().ok()?;
    let ref_allele = fields[2].to_string();
    let alt_allele = fields[3].to_string();
    let score: f64 = fields[8].trim().parse().ok()?;
    let class_idx = AM_CLASSES
        .iter()
        .position(|c| *c == fields[9].trim())
        .map(|i| i as u32);

    Some(AmRow {
        chrom_idx,
        position,
        ref_allele,
        alt_allele,
        score,
        class_idx,
    })
}

fn normalize_chrom(c: &str) -> String {
    if c.starts_with("chr") {
        c.to_string()
    } else {
        format!("chr{}", c)
    }
}

/// Build the flat-object JSON for one row using the very same `Field` encode →
/// `format_value` path the v2 reader reconstructs with, so v1 `.osa` and v2
/// `.osa2` emit byte-identical annotations. Fields are emitted in config order
/// (`amPathogenicity`, then `amClass`); a missing value is omitted, matching
/// both the reader's `reconstruct_json` and the v1 convention.
fn row_json(row: &AmRow, fields: &[Field], class_table: &[String]) -> String {
    let mut parts: Vec<String> = Vec::with_capacity(fields.len());

    let score_field = &fields[SCORE_FIELD];
    let stored = score_field.encode_float(row.score);
    if stored != score_field.missing_value {
        let v = format_value(score_field, stored, None);
        if v != "null" {
            parts.push(format!("\"{}\":{}", score_field.alias, v));
        }
    }

    let class_field = &fields[CLASS_FIELD];
    if let Some(idx) = row.class_idx {
        let v = format_value(class_field, idx, Some(class_table));
        if v != "null" {
            parts.push(format!("\"{}\":{}", class_field.alias, v));
        }
    }

    format!("{{{}}}", parts.join(","))
}

/// Stream an AlphaMissense TSV as v1 `AnnotationRecord`s without buffering the
/// whole file. The full hg38 release is ~71M rows, so streaming is required.
/// The input must already be coordinate-sorted (all AlphaMissense releases are).
pub fn iter_alphamissense_tsv<'a, R: BufRead>(
    reader: R,
    chrom_to_idx: &'a HashMap<String, u16>,
) -> AlphaMissenseV1Iter<'a, R> {
    AlphaMissenseV1Iter {
        lines: reader.lines(),
        chrom_to_idx,
        fields: alphamissense_osa2_fields(),
        class_table: AM_CLASSES.iter().map(|s| s.to_string()).collect(),
    }
}

pub struct AlphaMissenseV1Iter<'a, R: BufRead> {
    lines: std::io::Lines<R>,
    chrom_to_idx: &'a HashMap<String, u16>,
    fields: Vec<Field>,
    class_table: Vec<String>,
}

impl<R: BufRead> Iterator for AlphaMissenseV1Iter<'_, R> {
    type Item = Result<AnnotationRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            let line = match self.lines.next()? {
                Ok(l) => l,
                Err(e) => return Some(Err(e).context("Reading AlphaMissense line")),
            };
            let Some(row) = parse_row(&line, self.chrom_to_idx) else {
                continue;
            };
            let json = row_json(&row, &self.fields, &self.class_table);
            return Some(Ok(AnnotationRecord {
                chrom_idx: row.chrom_idx,
                position: row.position,
                ref_allele: row.ref_allele,
                alt_allele: row.alt_allele,
                json,
            }));
        }
    }
}

/// Stream an AlphaMissense TSV directly as v2 `Osa2Record`s (score column +
/// class-index column). Requires `chrom_list` to map the parser's chrom index
/// back to the on-disk chromosome name the `.osa2` layout keys on.
pub fn iter_alphamissense_osa2<'a, R: BufRead>(
    reader: R,
    chrom_to_idx: &'a HashMap<String, u16>,
    chrom_list: &'a [String],
) -> AlphaMissenseV2Iter<'a, R> {
    AlphaMissenseV2Iter {
        lines: reader.lines(),
        chrom_to_idx,
        chrom_list,
        fields: alphamissense_osa2_fields(),
    }
}

pub struct AlphaMissenseV2Iter<'a, R: BufRead> {
    lines: std::io::Lines<R>,
    chrom_to_idx: &'a HashMap<String, u16>,
    chrom_list: &'a [String],
    fields: Vec<Field>,
}

impl<R: BufRead> Iterator for AlphaMissenseV2Iter<'_, R> {
    type Item = Result<Osa2Record>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            let line = match self.lines.next()? {
                Ok(l) => l,
                Err(e) => return Some(Err(e).context("Reading AlphaMissense line")),
            };
            let Some(row) = parse_row(&line, self.chrom_to_idx) else {
                continue;
            };
            let chrom = match self.chrom_list.get(row.chrom_idx as usize) {
                Some(c) => c.clone(),
                None => {
                    return Some(Err(anyhow::anyhow!(
                        "chrom_idx {} out of range",
                        row.chrom_idx
                    )))
                }
            };
            let values = vec![
                self.fields[SCORE_FIELD].encode_float(row.score),
                row.class_idx.unwrap_or(self.fields[CLASS_FIELD].missing_value),
            ];
            return Some(Ok(Osa2Record {
                chrom,
                position: row.position,
                ref_allele: row.ref_allele.into_bytes(),
                alt_allele: row.alt_allele.into_bytes(),
                values,
                json_blob: None,
            }));
        }
    }
}

/// Parse an AlphaMissense TSV into sorted v1 `AnnotationRecord`s.
///
/// Loads everything into memory — for tests and small inputs. The full release
/// should stream through `iter_alphamissense_tsv` via the pipeline instead.
pub fn parse_alphamissense<R: BufRead>(
    reader: R,
    chrom_to_idx: &HashMap<String, u16>,
) -> Result<Vec<AnnotationRecord>> {
    let mut records: Vec<_> =
        iter_alphamissense_tsv(reader, chrom_to_idx).collect::<Result<_>>()?;
    records.sort_by(|a, b| a.chrom_idx.cmp(&b.chrom_idx).then(a.position.cmp(&b.position)));
    Ok(records)
}

#[cfg(test)]
mod tests {
    use super::*;

    const SAMPLE: &str = "\
# Copyright 2023 DeepMind Technologies Limited
#
#CHROM\tPOS\tREF\tALT\tgenome\tuniprot_id\ttranscript_id\tprotein_variant\tam_pathogenicity\tam_class
chr1\t69094\tG\tT\thg38\tQ8NH21\tENST00000335137.4\tV2L\t0.2937\tlikely_benign
chr1\t69095\tT\tA\thg38\tQ8NH21\tENST00000335137.4\tF2Y\t0.9021\tlikely_pathogenic
chr1\t69095\tT\tC\thg38\tQ8NH21\tENST00000335137.4\tF2L\t0.4500\tambiguous
";

    fn chrom_map() -> HashMap<String, u16> {
        let mut m = HashMap::new();
        m.insert("chr1".into(), 0u16);
        m
    }

    #[test]
    fn parses_rows_skipping_header_and_comments() {
        let recs = parse_alphamissense(SAMPLE.as_bytes(), &chrom_map()).unwrap();
        assert_eq!(recs.len(), 3);
        assert_eq!(recs[0].position, 69094);
        assert_eq!(recs[0].ref_allele, "G");
        assert_eq!(recs[0].alt_allele, "T");
    }

    #[test]
    fn v1_json_carries_score_and_class() {
        let recs = parse_alphamissense(SAMPLE.as_bytes(), &chrom_map()).unwrap();
        // Score is emitted in scientific notation (shared float formatting).
        assert!(recs[0].json.contains("\"amPathogenicity\":"), "{}", recs[0].json);
        assert!(recs[0].json.contains("\"amClass\":\"likely_benign\""), "{}", recs[0].json);
        assert!(recs[1].json.contains("\"amClass\":\"likely_pathogenic\""), "{}", recs[1].json);
        assert!(recs[2].json.contains("\"amClass\":\"ambiguous\""), "{}", recs[2].json);
    }

    #[test]
    fn v2_records_encode_score_and_class_index() {
        let chrom_list = vec!["chr1".to_string()];
        let recs: Vec<_> = iter_alphamissense_osa2(SAMPLE.as_bytes(), &chrom_map(), &chrom_list)
            .collect::<Result<_>>()
            .unwrap();
        assert_eq!(recs.len(), 3);
        let fields = alphamissense_osa2_fields();
        assert_eq!(recs[0].values[SCORE_FIELD], fields[SCORE_FIELD].encode_float(0.2937));
        assert_eq!(recs[0].values[CLASS_FIELD], 0); // likely_benign
        assert_eq!(recs[1].values[CLASS_FIELD], 2); // likely_pathogenic
        assert_eq!(recs[2].values[CLASS_FIELD], 1); // ambiguous
    }

    #[test]
    fn v1_and_v2_agree_on_json_shape() {
        // The v1 JSON is built through the same Field/format_value path the v2
        // reader reconstructs with, so a row's v1 JSON must match a manual
        // reconstruction from the v2-encoded values.
        let chrom_list = vec!["chr1".to_string()];
        let v1 = parse_alphamissense(SAMPLE.as_bytes(), &chrom_map()).unwrap();
        let fields = alphamissense_osa2_fields();
        let table: Vec<String> = AM_CLASSES.iter().map(|s| s.to_string()).collect();
        let v2: Vec<_> = iter_alphamissense_osa2(SAMPLE.as_bytes(), &chrom_map(), &chrom_list)
            .collect::<Result<_>>()
            .unwrap();
        for (r1, r2) in v1.iter().zip(v2.iter()) {
            let score = format_value(&fields[SCORE_FIELD], r2.values[SCORE_FIELD], None);
            let class = format_value(&fields[CLASS_FIELD], r2.values[CLASS_FIELD], Some(&table));
            let expected = format!(
                "{{\"amPathogenicity\":{},\"amClass\":{}}}",
                score, class
            );
            assert_eq!(r1.json, expected);
        }
    }

    #[test]
    fn unknown_class_encodes_missing() {
        let line = "chr1\t100\tA\tG\thg38\tU\tENST1\tX1Y\t0.5\tsomething_else";
        let row = parse_row(line, &chrom_map()).unwrap();
        assert_eq!(row.class_idx, None);
    }
}
