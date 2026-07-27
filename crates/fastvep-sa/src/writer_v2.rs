//! Writer for .osa2 format (ZIP-based chunked annotation files).
//!
//! Organizes annotations into ~1MB genomic chunks with parallel u32 value
//! arrays, sorted Var32 keys, and delta encoding for efficient compression.

use crate::chunk::delta_encode;
use crate::common::AnnotationRecord;
use crate::fields::{Field, FieldType};
use crate::kmer16::{self, LongVariant};
use crate::var32;
use anyhow::{Context, Result};
use std::io::{Seek, Write};

/// Metadata for the .osa2 file.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct Osa2Metadata {
    pub format_version: u32,
    pub name: String,
    pub version: String,
    pub assembly: String,
    pub json_key: String,
    pub match_by_allele: bool,
    pub is_array: bool,
    pub is_positional: bool,
    pub chunk_bits: u32,
    pub description: String,
}

/// A record to write into the .osa2 file.
#[derive(Debug, Clone)]
pub struct Osa2Record {
    pub chrom: String,
    pub position: u32,
    pub ref_allele: Vec<u8>,
    pub alt_allele: Vec<u8>,
    /// Parallel field values (same order as the Field config).
    pub values: Vec<u32>,
    /// Optional JSON blob for JsonBlob fields.
    pub json_blob: Option<String>,
}

/// Builds an .osa2 file from sorted records.
pub struct Osa2Writer {
    metadata: Osa2Metadata,
    fields: Vec<Field>,
    /// Categorical string tables: field_idx -> list of unique strings.
    string_tables: Vec<Vec<String>>,
}

impl Osa2Writer {
    pub fn new(metadata: Osa2Metadata, fields: Vec<Field>) -> Self {
        let string_tables = fields.iter().map(|_| Vec::new()).collect();
        Self { metadata, fields, string_tables }
    }

    /// Set the string table for a categorical field.
    pub fn set_string_table(&mut self, field_idx: usize, strings: Vec<String>) {
        if field_idx < self.string_tables.len() {
            self.string_tables[field_idx] = strings;
        }
    }

    /// Write all records to a .osa2 ZIP file.
    ///
    /// Records MUST be sorted by (chrom, position) so that all records for a
    /// given (chrom, chunk) are contiguous.
    pub fn write_all<W: Write + Seek>(
        &self,
        writer: W,
        records: &[Osa2Record],
    ) -> Result<()> {
        let mut zip = zip::ZipWriter::new(writer);
        let options = default_options();

        write_prelude(&mut zip, options, &self.metadata, &self.fields)?;
        write_string_tables(&mut zip, options, &self.fields, &self.string_tables)?;

        // Walk contiguous (chrom, chunk_id) runs. Because `records` is sorted,
        // every record for a chunk is adjacent, so a single pass slices them
        // without an intermediate index map.
        let chunk_bits = self.metadata.chunk_bits;
        let mut start = 0;
        while start < records.len() {
            let chrom = records[start].chrom.as_str();
            let cid = records[start].position >> chunk_bits;
            let mut end = start + 1;
            while end < records.len()
                && records[end].chrom == chrom
                && (records[end].position >> chunk_bits) == cid
            {
                end += 1;
            }
            write_chunk_entries(
                &mut zip,
                options,
                &records[start..end],
                &self.fields,
                chrom,
                cid,
                chunk_bits,
                self.metadata.is_positional,
            )?;
            start = end;
        }

        // `ZipWriter::finish` returns the inner writer without flushing it.
        // When that writer is a `BufWriter`, dropping it flushes but *discards*
        // any IO error (e.g. a disk that fills on the final write), which would
        // leave a truncated archive reported as success. Flush explicitly so
        // the error propagates and the caller's cleanup removes the partial file.
        let mut inner = zip.finish()?;
        inner.flush()?;
        Ok(())
    }
}

/// Default ZIP entry options (Deflate compression) shared by both writers.
fn default_options() -> zip::write::SimpleFileOptions {
    zip::write::SimpleFileOptions::default().compression_method(zip::CompressionMethod::Deflated)
}

/// Write the metadata + field-config entries that every .osa2 archive begins
/// with.
fn write_prelude<W: Write + Seek>(
    zip: &mut zip::ZipWriter<W>,
    options: zip::write::SimpleFileOptions,
    metadata: &Osa2Metadata,
    fields: &[Field],
) -> Result<()> {
    zip.start_file("fastsa/metadata.json", options)?;
    serde_json::to_writer_pretty(&mut *zip, metadata)?;
    zip.start_file("fastsa/config.json", options)?;
    serde_json::to_writer_pretty(&mut *zip, fields)?;
    Ok(())
}

/// Write categorical string tables. The reader looks these up by name, so
/// their position within the archive is irrelevant — the streaming writer
/// emits them last, the buffered writer emits them up front.
fn write_string_tables<W: Write + Seek>(
    zip: &mut zip::ZipWriter<W>,
    options: zip::write::SimpleFileOptions,
    fields: &[Field],
    string_tables: &[Vec<String>],
) -> Result<()> {
    for (i, field) in fields.iter().enumerate() {
        if field.ftype == FieldType::Categorical && !string_tables[i].is_empty() {
            zip.start_file(format!("fastsa/strings/{}.txt", field.alias), options)?;
            for s in &string_tables[i] {
                writeln!(zip, "{}", s)?;
            }
        }
    }
    Ok(())
}

/// Encode one chunk's records (already restricted to a single (chrom, chunk))
/// into the archive.
///
/// The value/JSON-blob columns are laid out as
/// `[short variants in Var32-sorted order] ++ [long variants in sorted order]`,
/// and each `LongVariant.idx` is set to the slot that variant occupies in that
/// combined layout. This keeps a single value slot per variant that both the
/// short (`binary_search` position) and long (`LongVariant.idx`) lookup paths
/// resolve to consistently. An earlier revision set `idx` to the record's
/// input-order position and only stored short variants in the value columns,
/// so any chunk mixing short and long variants returned the wrong values for
/// its indels (and long-only chunks were unreadable).
fn write_chunk_entries<W: Write + Seek>(
    zip: &mut zip::ZipWriter<W>,
    options: zip::write::SimpleFileOptions,
    chunk_records: &[Osa2Record],
    fields: &[Field],
    chrom: &str,
    chunk_id: u32,
    chunk_bits: u32,
    is_positional: bool,
) -> Result<()> {
    let chunk_mask = (1u32 << chunk_bits) - 1;

    // Partition into short (Var32) and long (kmer16) entries, remembering each
    // one's position within `chunk_records` so the value columns can be built
    // in the combined short-then-long order.
    let mut short_entries: Vec<(u32, usize)> = Vec::new(); // (var32_key, local_idx)
    let mut long_entries: Vec<(LongVariant, usize)> = Vec::new(); // (variant, local_idx)

    for (local_idx, record) in chunk_records.iter().enumerate() {
        let within_chunk_pos = record.position & chunk_mask;

        if is_positional {
            // Positional sources match by coordinate alone; key on position
            // only (alleles are empty) and never take the long-variant path.
            short_entries.push((var32::positional_key(within_chunk_pos), local_idx));
        } else if var32::is_long(record.ref_allele.len(), record.alt_allele.len()) {
            // Skip variants whose alleles contain non-ACGT bases. Earlier
            // revisions silently encoded them as runs of 'T', producing index
            // entries that could never be retrieved with their original allele
            // string.
            let Some(sequence) = kmer16::encode_var(&record.ref_allele, &record.alt_allele) else {
                log::warn!(
                    "Skipping long variant at {}:{} with non-ACGT allele (ref={:?} alt={:?})",
                    chrom,
                    record.position,
                    String::from_utf8_lossy(&record.ref_allele),
                    String::from_utf8_lossy(&record.alt_allele),
                );
                continue;
            };
            long_entries.push((
                LongVariant { position: record.position, idx: 0, sequence },
                local_idx,
            ));
        } else if let Some(key) =
            var32::encode(within_chunk_pos, &record.ref_allele, &record.alt_allele)
        {
            short_entries.push((key, local_idx));
        } else {
            // Short variant that fails Var32 encoding only when it contains a
            // non-ACGT base; same skip-with-warning policy.
            log::warn!(
                "Skipping short variant at {}:{} with non-ACGT allele (ref={:?} alt={:?})",
                chrom,
                record.position,
                String::from_utf8_lossy(&record.ref_allele),
                String::from_utf8_lossy(&record.alt_allele),
            );
        }
    }

    short_entries.sort_by_key(|(key, _)| *key);
    long_entries.sort_by(|(a, _), (b, _)| a.cmp(b));

    // Value slot layout: short variants first (parallel to the sorted var32
    // key array, so `binary_search` positions map directly), then long
    // variants. Assign each long variant its slot index.
    let short_count = short_entries.len();
    for (rank, (lv, _)) in long_entries.iter_mut().enumerate() {
        lv.idx = (short_count + rank) as u32;
    }
    let value_order: Vec<usize> = short_entries
        .iter()
        .map(|(_, li)| *li)
        .chain(long_entries.iter().map(|(_, li)| *li))
        .collect();

    let var32s: Vec<u32> = short_entries.iter().map(|(k, _)| *k).collect();
    let delta_var32s = delta_encode(&var32s);

    let prefix = format!("fastsa/{}/{}/", chrom, chunk_id);
    zip.start_file(format!("{}var32.bin", prefix), options)?;
    write_u32_array(zip, &delta_var32s)?;

    if !long_entries.is_empty() {
        let longs: Vec<&LongVariant> = long_entries.iter().map(|(lv, _)| lv).collect();
        zip.start_file(format!("{}too-long.enc", prefix), options)?;
        let data = bincode::serialize(&longs)?;
        zip.write_all(&data)?;
    }

    // Parallel value arrays, one per non-JsonBlob field, in the combined
    // short-then-long value-slot order.
    for (fi, field) in fields.iter().enumerate() {
        if field.ftype == FieldType::JsonBlob {
            continue;
        }
        let values: Vec<u32> = value_order
            .iter()
            .map(|&li| chunk_records[li].values.get(fi).copied().unwrap_or(field.missing_value))
            .collect();
        zip.start_file(format!("{}{}.bin", prefix, field.alias), options)?;
        write_u32_array(zip, &values)?;
    }

    // JSON blobs, if any field uses them.
    if fields.iter().any(|f| f.ftype == FieldType::JsonBlob) {
        let blobs: Vec<&str> = value_order
            .iter()
            .map(|&li| chunk_records[li].json_blob.as_deref().unwrap_or(""))
            .collect();
        if blobs.iter().any(|b| !b.is_empty()) {
            zip.start_file(format!("{}json_blobs.zst", prefix), options)?;
            let joined = blobs.join("\n");
            let compressed = zstd::encode_all(joined.as_bytes(), 3)?;
            zip.write_all(&compressed)?;
        }
    }

    Ok(())
}

/// Streaming writer for `.osa2` files. Consumes records sorted by
/// (chrom, position) and flushes one chunk at a time, so peak memory is
/// bounded by the densest ~1 MB genomic window rather than the whole file.
/// This is the path used to build genome-scale sources (e.g. gnomAD).
///
/// Contract: records MUST arrive sorted by (chrom, position). A record that
/// reopens an already-flushed (chrom, chunk) is rejected, since it would emit
/// a duplicate ZIP entry and signals unsorted input.
pub struct Osa2StreamWriter<W: Write + Seek> {
    zip: zip::ZipWriter<W>,
    options: zip::write::SimpleFileOptions,
    fields: Vec<Field>,
    string_tables: Vec<Vec<String>>,
    chunk_bits: u32,
    /// Positional sources key on coordinate alone (see `write_chunk_entries`).
    is_positional: bool,
    /// Records accumulated for the chunk currently being filled.
    current: Option<(String, u32, Vec<Osa2Record>)>,
    /// (chrom, chunk_id) runs already flushed, to detect unsorted input.
    seen: std::collections::HashSet<(String, u32)>,
}

impl<W: Write + Seek> Osa2StreamWriter<W> {
    /// Create a streaming writer, emitting the metadata + config prelude
    /// immediately.
    pub fn new(writer: W, metadata: &Osa2Metadata, fields: Vec<Field>) -> Result<Self> {
        let mut zip = zip::ZipWriter::new(writer);
        let options = default_options();
        write_prelude(&mut zip, options, metadata, &fields)?;
        let string_tables = fields.iter().map(|_| Vec::new()).collect();
        Ok(Self {
            zip,
            options,
            chunk_bits: metadata.chunk_bits,
            is_positional: metadata.is_positional,
            fields,
            string_tables,
            current: None,
            seen: std::collections::HashSet::new(),
        })
    }

    /// Set the string table for a categorical field. Must be called before
    /// `finish`; tables are written when the archive is finalized.
    pub fn set_string_table(&mut self, field_idx: usize, strings: Vec<String>) {
        if field_idx < self.string_tables.len() {
            self.string_tables[field_idx] = strings;
        }
    }

    /// Add one record. Flushes the previous chunk when the (chrom, chunk_id)
    /// boundary is crossed.
    pub fn push(&mut self, record: Osa2Record) -> Result<()> {
        let cid = record.position >> self.chunk_bits;
        let same_chunk = matches!(
            &self.current,
            Some((chrom, ccid, _)) if *chrom == record.chrom && *ccid == cid
        );
        if same_chunk {
            self.current.as_mut().unwrap().2.push(record);
        } else {
            self.flush_current()?;
            if !self.seen.insert((record.chrom.clone(), cid)) {
                anyhow::bail!(
                    "records not sorted by (chrom, position): chunk {}/{} reopened after being written",
                    record.chrom,
                    cid
                );
            }
            self.current = Some((record.chrom.clone(), cid, vec![record]));
        }
        Ok(())
    }

    fn flush_current(&mut self) -> Result<()> {
        if let Some((chrom, cid, buf)) = self.current.take() {
            write_chunk_entries(
                &mut self.zip,
                self.options,
                &buf,
                &self.fields,
                &chrom,
                cid,
                self.chunk_bits,
                self.is_positional,
            )?;
        }
        Ok(())
    }

    /// Flush the final chunk, write string tables, and finalize the archive.
    pub fn finish(mut self) -> Result<()> {
        self.flush_current()?;
        write_string_tables(&mut self.zip, self.options, &self.fields, &self.string_tables)?;
        // Flush the inner writer explicitly — `ZipWriter::finish` hands it back
        // unflushed, and a `BufWriter`'s Drop swallows flush errors, which would
        // silently leave a truncated `.osa2` on a full disk. See `write_all`.
        let mut inner = self.zip.finish()?;
        inner.flush()?;
        Ok(())
    }
}

/// Field schema for a source stored as one opaque whole-record JSON blob per
/// variant, with no decomposed value columns.
///
/// Used for the string/array payloads (ClinVar, dbSNP, COSMIC) that don't fit
/// the parallel u32 layout — high-cardinality ID strings and nested arrays —
/// but still benefit substantially from v2's chunk-level zstd of the blob
/// column (it compresses a whole chunk's JSON records together, exploiting
/// cross-record redundancy the v1 per-block scheme can't). The single field's
/// **empty alias** is the signal to [`crate::chunk::Chunk::reconstruct_json`]
/// that the stored blob is the complete record object and must be emitted
/// verbatim rather than nested under a key — so v2 output is byte-identical to
/// the v1 `.osa` this record came from.
pub fn raw_json_blob_fields() -> Vec<Field> {
    vec![Field {
        field: String::new(),
        alias: String::new(),
        ftype: FieldType::JsonBlob,
        multiplier: 1,
        zigzag: false,
        missing_value: u32::MAX,
        missing_string: ".".into(),
        description: "Whole-record JSON blob".into(),
    }]
}

/// Bridge a v1 [`AnnotationRecord`] into a whole-record-blob [`Osa2Record`]:
/// the entire v1 JSON string is stored as the blob (paired with
/// [`raw_json_blob_fields`]), so the v2 reader returns exactly the bytes the v1
/// reader would have. Reuses the existing, well-tested v1 parsers wholesale for
/// sources whose payload is opaque to the numeric u32 encoding.
pub fn osa2_raw_blob_from_v1(record: &AnnotationRecord, chrom: String) -> Osa2Record {
    Osa2Record {
        chrom,
        position: record.position,
        ref_allele: record.ref_allele.as_bytes().to_vec(),
        alt_allele: record.alt_allele.as_bytes().to_vec(),
        values: Vec::new(),
        json_blob: Some(record.json.clone()),
    }
}

/// Build an [`Osa2Record`] from a v1 [`AnnotationRecord`] whose `json` is a
/// flat object of scalar values keyed by field alias.
///
/// This bridges the existing, well-tested v1 source parsers to the v2 format:
/// the parser emits `AnnotationRecord`s exactly as it does for the `.osa`
/// builder, and this function re-encodes each field's value into the parallel
/// u32 layout by reading it back out of the JSON object by alias. A field
/// absent from the JSON is stored as its `missing_value` (and thus omitted on
/// read), matching how the v1 output would have simply not carried that key.
///
/// Only `Float`, `Integer`, and `Flag` fields are supported — the shape
/// produced by frequency/score sources. `Categorical`/`JsonBlob` fields are
/// rejected, since those need a source-specific encoder (string tables / blob
/// handling) rather than this generic bridge.
pub fn osa2_record_from_v1(
    record: &AnnotationRecord,
    chrom: String,
    fields: &[Field],
) -> Result<Osa2Record> {
    let parsed: serde_json::Value = serde_json::from_str(&record.json)
        .with_context(|| format!("Parsing SA JSON at {}:{}", chrom, record.position))?;
    let obj = parsed.as_object();

    let mut values = Vec::with_capacity(fields.len());
    for field in fields {
        let raw = obj.and_then(|m| m.get(&field.alias));
        let encoded = match field.ftype {
            FieldType::Float => match raw.and_then(|v| v.as_f64()) {
                Some(f) => field.encode_float(f),
                None => field.missing_value,
            },
            FieldType::Integer => match raw {
                Some(v) if v.as_i64().is_some() => field.encode_int(v.as_i64().unwrap()),
                // JSON may carry an integer-valued float (e.g. 42.0); accept it.
                Some(v) if v.as_f64().is_some() => field.encode_int(v.as_f64().unwrap() as i64),
                _ => field.missing_value,
            },
            FieldType::Flag => match raw {
                Some(serde_json::Value::Bool(b)) => u32::from(*b),
                Some(v) if v.as_i64() == Some(1) => 1,
                _ => 0,
            },
            FieldType::Categorical | FieldType::JsonBlob => {
                anyhow::bail!(
                    "osa2_record_from_v1 does not support {:?} fields (field '{}'); \
                     use a source-specific encoder",
                    field.ftype,
                    field.alias
                );
            }
        };
        values.push(encoded);
    }

    Ok(Osa2Record {
        chrom,
        position: record.position,
        ref_allele: record.ref_allele.as_bytes().to_vec(),
        alt_allele: record.alt_allele.as_bytes().to_vec(),
        values,
        json_blob: None,
    })
}

/// Write a u32 array as [4B count][4B * count values].
fn write_u32_array<W: Write>(writer: &mut W, values: &[u32]) -> Result<()> {
    writer.write_all(&(values.len() as u32).to_le_bytes())?;
    for &v in values {
        writer.write_all(&v.to_le_bytes())?;
    }
    Ok(())
}

/// Read a u32 array from [4B count][4B * count values].
pub fn read_u32_array(data: &[u8]) -> Result<Vec<u32>> {
    if data.len() < 4 {
        anyhow::bail!("u32 array too short");
    }
    let count = u32::from_le_bytes(data[0..4].try_into()?) as usize;
    // `count` comes from an untrusted .osa2 chunk sub-file. Bound it against
    // how many 4-byte elements the remaining bytes could actually hold before
    // allocating, so a corrupted file claiming `count = u32::MAX` can't force
    // a multi-GB allocation attempt (mirrors block.rs's `decompress` guard).
    let remaining = data.len() - 4;
    if count > remaining / 4 {
        anyhow::bail!("u32 array claims {} elements, exceeds data size", count);
    }
    let mut values = Vec::with_capacity(count);
    let mut offset = 4;
    for _ in 0..count {
        if offset + 4 > data.len() {
            anyhow::bail!("u32 array truncated");
        }
        values.push(u32::from_le_bytes(data[offset..offset + 4].try_into()?));
        offset += 4;
    }
    Ok(values)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_u32_array_round_trip() {
        let values = vec![1, 2, 3, 100, 200];
        let mut buf = Vec::new();
        write_u32_array(&mut buf, &values).unwrap();
        let decoded = read_u32_array(&buf).unwrap();
        assert_eq!(decoded, values);
    }

    #[test]
    fn test_read_u32_array_rejects_oversized_claimed_count() {
        // A corrupted/hostile .osa2 chunk sub-file claiming `count = u32::MAX`
        // with only a few actual trailing bytes must be rejected up front,
        // not passed to Vec::with_capacity(count) (which would attempt a
        // multi-GB allocation and abort the process).
        let mut buf = Vec::new();
        buf.extend_from_slice(&u32::MAX.to_le_bytes());
        buf.extend_from_slice(&1u32.to_le_bytes()); // a few trailing bytes only
        let err = read_u32_array(&buf).unwrap_err();
        assert!(err.to_string().contains("exceeds data size"));
    }
}
