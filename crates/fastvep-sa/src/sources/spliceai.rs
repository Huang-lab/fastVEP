//! SpliceAI VCF parser for building .osa (v1) and .osa2 (v2) annotation files.
//!
//! SpliceAI provides splice site effect predictions with delta scores
//! for acceptor/donor gain/loss and their positions.
//!
//! SpliceAI is the densest supplementary source fastVEP ships against - the
//! precomputed whole-genome SNV release is billions of records - and its payload
//! is eight numbers plus a gene symbol. That makes it the best fit of any source
//! for the v2 `.osa2` layout: eight parallel u32 columns plus one categorical
//! index column against a gene string table, instead of a ~120-byte JSON object
//! repeated per record. So SpliceAI builds natively into v2 (see
//! [`iter_spliceai_osa2`]) rather than going through the whole-record JSON blob
//! bridge.
//!
//! As with AlphaMissense, the v1 `.osa` path builds each record's JSON through
//! the very same [`Field`]/[`format_value`] code the v2 reader reconstructs
//! with, so the two formats emit byte-identical annotations.

use crate::common::AnnotationRecord;
use crate::fields::{format_value, Field, FieldType};
use crate::writer_v2::{Osa2Metadata, Osa2Record};
use anyhow::{Context, Result};
use std::cell::RefCell;
use std::collections::HashMap;
use std::collections::VecDeque;
use std::io::BufRead;
use std::rc::Rc;

/// Field index of the gene symbol within [`spliceai_osa2_fields`]. The gene is
/// the only categorical field, so this is the index its string table is filed
/// under.
pub const GENE_FIELD: usize = 8;

/// Upper bound on distinct gene symbols admitted into the categorical string
/// table. Real SpliceAI releases carry a few tens of thousands of symbols; a
/// table running into the millions means the `gene` column of the input is not
/// a gene symbol at all, and silently accumulating it would grow unbounded in
/// memory during a multi-hour genome-scale build.
const MAX_GENE_TABLE: usize = 1_000_000;

/// Canonical SpliceAI `.osa2` field schema.
///
/// Field order is the JSON key order both formats emit. It is deliberately
/// alphabetical by alias (`dpAg`..`dsDl`, then `gene`) because that is the order
/// the historical v1 builder produced - it built the object with `serde_json`,
/// whose default `Map` is a `BTreeMap` - so routing v1 through `format_value`
/// changes only the *numeric formatting* of the payload, not its key order.
///
/// * delta scores (`DS_*`) are floats in `[0, 1]`; the released files carry two
///   decimals, so a 1e6 multiplier is comfortably lossless.
/// * delta positions (`DP_*`) are signed offsets in bases, hence zigzag
///   integers.
pub fn spliceai_osa2_fields() -> Vec<Field> {
    let dp = |field: &str, alias: &str, what: &str| Field {
        field: field.into(),
        alias: alias.into(),
        ftype: FieldType::Integer,
        multiplier: 1,
        zigzag: true,
        missing_value: u32::MAX,
        missing_string: ".".into(),
        description: format!("SpliceAI delta position, {what}"),
    };
    let ds = |field: &str, alias: &str, what: &str| Field {
        field: field.into(),
        alias: alias.into(),
        ftype: FieldType::Float,
        multiplier: 1_000_000,
        zigzag: false,
        missing_value: u32::MAX,
        missing_string: ".".into(),
        description: format!("SpliceAI delta score, {what}"),
    };
    vec![
        dp("DP_AG", "dpAg", "acceptor gain"),
        dp("DP_AL", "dpAl", "acceptor loss"),
        dp("DP_DG", "dpDg", "donor gain"),
        dp("DP_DL", "dpDl", "donor loss"),
        ds("DS_AG", "dsAg", "acceptor gain"),
        ds("DS_AL", "dsAl", "acceptor loss"),
        ds("DS_DG", "dsDg", "donor gain"),
        ds("DS_DL", "dsDl", "donor loss"),
        Field {
            field: "SYMBOL".into(),
            alias: "gene".into(),
            ftype: FieldType::Categorical,
            multiplier: 1,
            zigzag: false,
            missing_value: u32::MAX,
            missing_string: ".".into(),
            description: "Gene symbol the prediction is relative to".into(),
        },
    ]
}

/// Standard SpliceAI `.osa2` metadata (`json_key = "spliceAI"`, allele-matched,
/// non-positional) - the same view the v1 `.osa` header carries.
pub fn spliceai_osa2_metadata(assembly: &str) -> Osa2Metadata {
    Osa2Metadata {
        format_version: 2,
        name: "SpliceAI".into(),
        version: "latest".into(),
        assembly: assembly.into(),
        json_key: "spliceAI".into(),
        match_by_allele: true,
        is_array: false,
        is_positional: false,
        chunk_bits: 20,
        description: format!("SpliceAI splice site predictions for {assembly}"),
    }
}

/// Accumulates the distinct gene symbols seen while a SpliceAI input streams
/// past, assigning each a stable u32 index for the `gene` column.
///
/// The `.osa2` categorical string table is global to the archive, but the
/// symbols are only discovered as records go by, so the table cannot be known
/// up front the way AlphaMissense's fixed three-class table can. A cheap shared
/// handle lets the record iterator intern as it goes while the build driver
/// keeps a clone to hand the finished table to the writer just before the
/// archive is finalized. Single-threaded by construction: the v2 build loop
/// pushes records on one thread.
#[derive(Clone, Default)]
pub struct GeneInterner(Rc<RefCell<GeneTable>>);

#[derive(Default)]
struct GeneTable {
    index: HashMap<String, u32>,
    order: Vec<String>,
}

impl GeneInterner {
    pub fn new() -> Self {
        Self::default()
    }

    /// Index for `gene`, assigning a fresh one on first sight.
    ///
    /// Errors past [`MAX_GENE_TABLE`] distinct symbols rather than growing
    /// without bound - see that constant for why.
    pub fn intern(&self, gene: &str) -> Result<u32> {
        let mut table = self.0.borrow_mut();
        if let Some(&idx) = table.index.get(gene) {
            return Ok(idx);
        }
        if table.order.len() >= MAX_GENE_TABLE {
            anyhow::bail!(
                "SpliceAI gene column has more than {} distinct values (at '{}'); \
                 this does not look like a SpliceAI VCF",
                MAX_GENE_TABLE,
                gene
            );
        }
        let idx = table.order.len() as u32;
        table.order.push(gene.to_string());
        table.index.insert(gene.to_string(), idx);
        Ok(idx)
    }

    /// The string table in index order, ready for `set_string_table`.
    pub fn table(&self) -> Vec<String> {
        self.0.borrow().order.clone()
    }

    /// The `(field_idx, table)` pairing the v2 writer takes.
    pub fn string_tables(&self) -> Vec<(usize, Vec<String>)> {
        vec![(GENE_FIELD, self.table())]
    }
}

/// One parsed `SpliceAI=` entry: the variant key plus the nine stored values.
struct SpliceAiRow {
    chrom_idx: u16,
    position: u32,
    ref_allele: String,
    alt_allele: String,
    gene: String,
    /// `[DP_AG, DP_AL, DP_DG, DP_DL]`, in field-config order.
    dp: [i32; 4],
    /// `[DS_AG, DS_AL, DS_DG, DS_DL]`, in field-config order.
    ds: [f64; 4],
}

/// Parse a SpliceAI VCF and produce sorted AnnotationRecords.
///
/// SpliceAI INFO field format:
/// `SpliceAI=A|GENE|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL`
pub fn parse_spliceai_vcf<R: BufRead>(
    reader: R,
    chrom_to_idx: &HashMap<String, u16>,
) -> Result<Vec<AnnotationRecord>> {
    let mut records: Vec<_> = iter_spliceai_vcf(reader, chrom_to_idx).collect::<Result<_>>()?;
    records.sort_by(|a, b| a.chrom_idx.cmp(&b.chrom_idx).then(a.position.cmp(&b.position)));
    Ok(records)
}

/// Stream a coordinate-sorted SpliceAI VCF as AnnotationRecords.
///
/// This avoids retaining the whole genome-wide SpliceAI VCF in memory before
/// writing fastSA. The input must already be sorted by chromosome and position.
pub fn iter_spliceai_vcf<'a, R: BufRead>(
    reader: R,
    chrom_to_idx: &'a HashMap<String, u16>,
) -> SpliceAiRecordIter<'a, R> {
    SpliceAiRecordIter {
        rows: SpliceAiRowIter::new(reader, chrom_to_idx),
        fields: spliceai_osa2_fields(),
    }
}

pub struct SpliceAiRecordIter<'a, R: BufRead> {
    rows: SpliceAiRowIter<'a, R>,
    fields: Vec<Field>,
}

impl<R: BufRead> Iterator for SpliceAiRecordIter<'_, R> {
    type Item = Result<AnnotationRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        let row = match self.rows.next()? {
            Ok(row) => row,
            Err(e) => return Some(Err(e)),
        };
        let json = row_json(&row, &self.fields);
        Some(Ok(AnnotationRecord {
            chrom_idx: row.chrom_idx,
            position: row.position,
            ref_allele: row.ref_allele,
            alt_allele: row.alt_allele,
            json,
        }))
    }
}

/// Stream a coordinate-sorted SpliceAI VCF directly as v2 `Osa2Record`s:
/// eight numeric columns plus a `gene` index interned into `genes`.
///
/// Every parsed entry is emitted, including the several a variant in overlapping
/// genes produces for one (position, ref, alt) key. Collapsing those to the
/// first - which is the only one either reader can return - is the writer's job
/// (see `writer_v2::write_chunk_entries`), so that every source gets the same
/// treatment rather than each parser reimplementing it.
pub fn iter_spliceai_osa2<'a, R: BufRead>(
    reader: R,
    chrom_to_idx: &'a HashMap<String, u16>,
    chrom_list: &'a [String],
    genes: GeneInterner,
) -> SpliceAiOsa2Iter<'a, R> {
    SpliceAiOsa2Iter {
        rows: SpliceAiRowIter::new(reader, chrom_to_idx),
        chrom_list,
        fields: spliceai_osa2_fields(),
        genes,
    }
}

pub struct SpliceAiOsa2Iter<'a, R: BufRead> {
    rows: SpliceAiRowIter<'a, R>,
    chrom_list: &'a [String],
    fields: Vec<Field>,
    genes: GeneInterner,
}

impl<R: BufRead> Iterator for SpliceAiOsa2Iter<'_, R> {
    type Item = Result<Osa2Record>;

    fn next(&mut self) -> Option<Self::Item> {
        let row = match self.rows.next()? {
            Ok(row) => row,
            Err(e) => return Some(Err(e)),
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
        let gene_idx = match self.genes.intern(&row.gene) {
            Ok(idx) => idx,
            Err(e) => return Some(Err(e)),
        };
        let values = encode_values(&row, &self.fields, gene_idx);
        Some(Ok(Osa2Record {
            chrom,
            position: row.position,
            ref_allele: row.ref_allele.into_bytes(),
            alt_allele: row.alt_allele.into_bytes(),
            values,
            json_blob: None,
        }))
    }
}

/// Encode one row into the parallel u32 value vector, in field-config order.
fn encode_values(row: &SpliceAiRow, fields: &[Field], gene_idx: u32) -> Vec<u32> {
    let mut values = Vec::with_capacity(fields.len());
    for (i, dp) in row.dp.iter().enumerate() {
        values.push(fields[i].encode_int(*dp as i64));
    }
    for (i, ds) in row.ds.iter().enumerate() {
        values.push(fields[4 + i].encode_float(*ds));
    }
    values.push(gene_idx);
    debug_assert_eq!(values.len(), fields.len());
    values
}

/// Build the flat-object JSON for one row through the same `Field` encode →
/// `format_value` path the v2 reader reconstructs with, so v1 `.osa` and v2
/// `.osa2` emit byte-identical annotations. Fields are emitted in config order
/// and a missing value is omitted, matching `Chunk::reconstruct_json`.
fn row_json(row: &SpliceAiRow, fields: &[Field]) -> String {
    // The gene's stored index is 0 against a one-entry table: the v2 archive
    // uses a global table, but the rendered value depends only on the string it
    // resolves to, so a local table gives the identical output here.
    let gene_table = [row.gene.clone()];
    let values = encode_values(row, fields, 0);

    let mut parts: Vec<String> = Vec::with_capacity(fields.len());
    for (i, field) in fields.iter().enumerate() {
        let stored = values[i];
        if stored == field.missing_value {
            continue;
        }
        let strings = if field.ftype == FieldType::Categorical {
            Some(&gene_table[..])
        } else {
            None
        };
        let rendered = format_value(field, stored, strings);
        if rendered != "null" {
            parts.push(format!("\"{}\":{}", field.alias, rendered));
        }
    }
    format!("{{{}}}", parts.join(","))
}

/// The shared front end: turns VCF lines into [`SpliceAiRow`]s. Both the v1 and
/// v2 iterators wrap this so the two formats can never drift on parsing,
/// chromosome normalization, or which records are skipped.
struct SpliceAiRowIter<'a, R: BufRead> {
    lines: std::io::Lines<R>,
    chrom_to_idx: &'a HashMap<String, u16>,
    pending: VecDeque<SpliceAiRow>,
}

impl<'a, R: BufRead> SpliceAiRowIter<'a, R> {
    fn new(reader: R, chrom_to_idx: &'a HashMap<String, u16>) -> Self {
        Self {
            lines: reader.lines(),
            chrom_to_idx,
            pending: VecDeque::new(),
        }
    }
}

impl<R: BufRead> Iterator for SpliceAiRowIter<'_, R> {
    type Item = Result<SpliceAiRow>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            if let Some(row) = self.pending.pop_front() {
                return Some(Ok(row));
            }

            let line = match self.lines.next()? {
                Ok(line) => line,
                Err(err) => return Some(Err(err).context("Reading SpliceAI VCF line")),
            };
            if line.starts_with('#') {
                continue;
            }

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

            let ref_allele = fields[3];
            let info = fields[7];

            for pair in info.split(';') {
                if let Some(val) = pair.strip_prefix("SpliceAI=") {
                    for entry in val.split(',') {
                        if let Some(row) = parse_spliceai_entry(chrom_idx, pos, ref_allele, entry) {
                            self.pending.push_back(row);
                        }
                    }
                }
            }
        }
    }
}

fn parse_spliceai_entry(
    chrom_idx: u16,
    position: u32,
    ref_allele: &str,
    entry: &str,
) -> Option<SpliceAiRow> {
    let parts: Vec<&str> = entry.split('|').collect();
    if parts.len() < 10 {
        return None;
    }

    Some(SpliceAiRow {
        chrom_idx,
        position,
        ref_allele: ref_allele.to_string(),
        alt_allele: parts[0].to_string(),
        gene: parts[1].to_string(),
        ds: [
            parts[2].parse().ok()?,
            parts[3].parse().ok()?,
            parts[4].parse().ok()?,
            parts[5].parse().ok()?,
        ],
        dp: [
            parts[6].parse().ok()?,
            parts[7].parse().ok()?,
            parts[8].parse().ok()?,
            parts[9].parse().ok()?,
        ],
    })
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

    const SAMPLE: &str = "\
##fileformat=VCFv4.0
##INFO=<ID=SpliceAI,Number=.,Type=String>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
1\t25000\t.\tA\tG\t.\t.\tSpliceAI=G|GENE1|0.01|0.00|0.85|0.00|5|-28|2|-13
1\t30000\t.\tC\tT,A\t.\t.\tSpliceAI=T|GENE2|0.00|0.10|0.00|0.92|3|-5|10|-2,A|GENE2|0.50|0.00|0.00|0.00|7|0|0|0
";

    fn chrom_map() -> HashMap<String, u16> {
        let mut m = HashMap::new();
        m.insert("chr1".into(), 0u16);
        m
    }

    fn chrom_list() -> Vec<String> {
        vec!["chr1".to_string()]
    }

    /// Field indices, for readability in the assertions below.
    const DP_AG: usize = 0;
    const DS_AG: usize = 4;
    const DS_DG: usize = 6;
    const DS_DL: usize = 7;

    #[test]
    fn test_parse_spliceai_vcf() {
        let records = parse_spliceai_vcf(SAMPLE.as_bytes(), &chrom_map()).unwrap();
        assert_eq!(records.len(), 3);

        assert_eq!(records[0].position, 25000);
        assert_eq!(records[0].alt_allele, "G");
        let first: serde_json::Value = serde_json::from_str(&records[0].json).unwrap();
        assert_eq!(first["gene"], "GENE1");
        assert!((first["dsDg"].as_f64().unwrap() - 0.85).abs() < 1e-9);

        assert_eq!(records[1].position, 30000);
        assert_eq!(records[1].alt_allele, "T");
        let second: serde_json::Value = serde_json::from_str(&records[1].json).unwrap();
        assert!((second["dsDl"].as_f64().unwrap() - 0.92).abs() < 1e-9);

        assert_eq!(records[2].position, 30000);
        assert_eq!(records[2].alt_allele, "A");
        let third: serde_json::Value = serde_json::from_str(&records[2].json).unwrap();
        assert!((third["dsAg"].as_f64().unwrap() - 0.50).abs() < 1e-9);
    }

    #[test]
    fn test_iter_spliceai_vcf_streams_records() {
        let records: Vec<_> = iter_spliceai_vcf(SAMPLE.as_bytes(), &chrom_map())
            .collect::<Result<Vec<_>>>()
            .unwrap();

        assert_eq!(records.len(), 3);
        assert_eq!(records[0].position, 25000);
        assert_eq!(records[0].alt_allele, "G");
        assert_eq!(records[1].position, 30000);
        assert_eq!(records[1].alt_allele, "T");
        assert_eq!(records[2].position, 30000);
        assert_eq!(records[2].alt_allele, "A");
    }

    #[test]
    fn test_parse_spliceai_vcf_escapes_gene_for_json() {
        let vcf = "\
##fileformat=VCFv4.0
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
1\t25000\t.\tA\tG\t.\t.\tSpliceAI=G|GENE\"1|0.01|0.00|0.85|0.00|5|-28|2|-13
";
        let records = parse_spliceai_vcf(vcf.as_bytes(), &chrom_map()).unwrap();
        assert_eq!(records.len(), 1);
        let value: serde_json::Value = serde_json::from_str(&records[0].json).unwrap();
        assert_eq!(value["gene"], "GENE\"1");
    }

    #[test]
    fn v1_json_keys_stay_in_the_historical_alphabetical_order() {
        // The pre-v2 builder emitted the object through serde_json, whose
        // default Map is a BTreeMap, so keys came out alphabetically. Routing
        // v1 through `format_value` must not reshuffle them.
        let records = parse_spliceai_vcf(SAMPLE.as_bytes(), &chrom_map()).unwrap();
        let keys: Vec<&str> = records[0]
            .json
            .trim_matches(|c| c == '{' || c == '}')
            .split(',')
            .map(|kv| kv.split(':').next().unwrap().trim_matches('"'))
            .collect();
        assert_eq!(
            keys,
            vec!["dpAg", "dpAl", "dpDg", "dpDl", "dsAg", "dsAl", "dsDg", "dsDl", "gene"]
        );
    }

    #[test]
    fn v2_records_encode_scores_positions_and_gene_index() {
        let genes = GeneInterner::new();
        let recs: Vec<_> =
            iter_spliceai_osa2(SAMPLE.as_bytes(), &chrom_map(), &chrom_list(), genes.clone())
                .collect::<Result<_>>()
                .unwrap();
        assert_eq!(recs.len(), 3);

        let fields = spliceai_osa2_fields();
        assert_eq!(recs[0].chrom, "chr1");
        assert_eq!(recs[0].position, 25000);
        assert_eq!(recs[0].values[DS_DG], fields[DS_DG].encode_float(0.85));
        assert_eq!(recs[0].values[DP_AG], fields[DP_AG].encode_int(5));
        // Negative delta positions survive via zigzag.
        assert_eq!(fields[1].decode_int(recs[0].values[1]), -28);

        assert_eq!(recs[1].values[DS_DL], fields[DS_DL].encode_float(0.92));
        assert_eq!(recs[2].values[DS_AG], fields[DS_AG].encode_float(0.50));

        // Two distinct symbols across the sample, interned in first-seen order.
        assert_eq!(genes.table(), vec!["GENE1".to_string(), "GENE2".to_string()]);
        assert_eq!(recs[0].values[GENE_FIELD], 0);
        assert_eq!(recs[1].values[GENE_FIELD], 1);
        assert_eq!(recs[2].values[GENE_FIELD], 1);
        assert_eq!(genes.string_tables(), vec![(GENE_FIELD, genes.table())]);
    }

    #[test]
    fn v1_and_v2_reconstruct_identical_json() {
        // The whole point of routing both formats through Field/format_value:
        // a record's v1 JSON must equal what the v2 reader rebuilds from the
        // encoded columns and the global gene table.
        let genes = GeneInterner::new();
        let v1 = parse_spliceai_vcf(SAMPLE.as_bytes(), &chrom_map()).unwrap();
        let v2: Vec<_> =
            iter_spliceai_osa2(SAMPLE.as_bytes(), &chrom_map(), &chrom_list(), genes.clone())
                .collect::<Result<_>>()
                .unwrap();
        let fields = spliceai_osa2_fields();
        let table = genes.table();

        assert_eq!(v1.len(), v2.len());
        for (r1, r2) in v1.iter().zip(v2.iter()) {
            let mut parts = Vec::new();
            for (i, field) in fields.iter().enumerate() {
                let stored = r2.values[i];
                if stored == field.missing_value {
                    continue;
                }
                let strings = if field.ftype == FieldType::Categorical {
                    Some(&table[..])
                } else {
                    None
                };
                let rendered = format_value(field, stored, strings);
                if rendered != "null" {
                    parts.push(format!("\"{}\":{}", field.alias, rendered));
                }
            }
            assert_eq!(r1.json, format!("{{{}}}", parts.join(",")));
        }
    }

    #[test]
    fn multi_gene_entries_are_all_emitted_in_input_order() {
        // A variant in two overlapping genes yields two entries for one
        // (position, ref, alt) key. Both iterators emit both, in input order —
        // collapsing to the first is the writer's job, and it relies on that
        // order to know which one "first" is.
        let vcf = "\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
1\t500\t.\tA\tG\t.\t.\tSpliceAI=G|FIRST|0.10|0.00|0.00|0.00|1|0|0|0,G|SECOND|0.90|0.00|0.00|0.00|2|0|0|0
1\t500\t.\tA\tT\t.\t.\tSpliceAI=T|THIRD|0.30|0.00|0.00|0.00|3|0|0|0
";
        let v1 = parse_spliceai_vcf(vcf.as_bytes(), &chrom_map()).unwrap();
        assert_eq!(v1.len(), 3);

        let genes = GeneInterner::new();
        let v2: Vec<_> =
            iter_spliceai_osa2(vcf.as_bytes(), &chrom_map(), &chrom_list(), genes.clone())
                .collect::<Result<_>>()
                .unwrap();
        assert_eq!(v2.len(), 3, "the iterator must not drop the second gene");

        let fields = spliceai_osa2_fields();
        let gene_of = |r: &Osa2Record| genes.table()[r.values[GENE_FIELD] as usize].clone();
        assert_eq!(gene_of(&v2[0]), "FIRST");
        assert_eq!(v2[0].values[DS_AG], fields[DS_AG].encode_float(0.10));
        assert_eq!(gene_of(&v2[1]), "SECOND");
        assert_eq!(v2[1].values[DS_AG], fields[DS_AG].encode_float(0.90));
        assert_eq!(gene_of(&v2[2]), "THIRD");
        assert_eq!(v2[2].alt_allele, b"T");
    }

    #[test]
    fn gene_interner_assigns_stable_indices_and_dedups() {
        let genes = GeneInterner::new();
        assert_eq!(genes.intern("A").unwrap(), 0);
        assert_eq!(genes.intern("B").unwrap(), 1);
        assert_eq!(genes.intern("A").unwrap(), 0);
        assert_eq!(genes.table(), vec!["A".to_string(), "B".to_string()]);
        // Clones share the same table.
        let other = genes.clone();
        assert_eq!(other.intern("C").unwrap(), 2);
        assert_eq!(genes.table().len(), 3);
    }

    #[test]
    fn malformed_entries_are_skipped_by_both_formats() {
        // Too few pipe-separated parts, and an unparseable score.
        let vcf = "\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
1\t100\t.\tA\tG\t.\t.\tSpliceAI=G|GENE1|0.01
1\t200\t.\tA\tG\t.\t.\tSpliceAI=G|GENE1|abc|0.00|0.00|0.00|1|2|3|4
1\t300\t.\tA\tG\t.\t.\tSpliceAI=G|GENE1|0.01|0.00|0.00|0.00|1|2|3|4
";
        let v1 = parse_spliceai_vcf(vcf.as_bytes(), &chrom_map()).unwrap();
        let v2: Vec<_> =
            iter_spliceai_osa2(vcf.as_bytes(), &chrom_map(), &chrom_list(), GeneInterner::new())
                .collect::<Result<_>>()
                .unwrap();
        assert_eq!(v1.len(), 1);
        assert_eq!(v2.len(), 1);
        assert_eq!(v1[0].position, 300);
        assert_eq!(v2[0].position, 300);
    }
}
