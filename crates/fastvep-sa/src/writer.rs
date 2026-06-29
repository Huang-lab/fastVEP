//! Writer for .osa position/allele-level annotation files.
//!
//! Records must be added in chromosome-sorted, position-sorted order.
//! The writer accumulates entries into blocks, compresses them, and writes
//! to the data file while building the index.

use crate::block::{BlockEntry, SaBlock};
use crate::common::{AnnotationRecord, DEFAULT_BLOCK_SIZE, OSA_MAGIC, SCHEMA_VERSION};
use crate::index::{BlockRef, IndexHeader, SaIndex};
use anyhow::{Context, Result};
use std::io::{BufWriter, Write};
use std::path::Path;

/// Builds an .osa data file and its .osa.idx index file.
pub struct SaWriter {
    index: SaIndex,
    block: SaBlock,
    current_chrom: Option<String>,
    last_key: Option<(u16, u32)>,
    /// Chromosome name -> numeric index mapping.
    chrom_names: Vec<String>,
    data_offset: u64,
}

impl SaWriter {
    pub fn new(header: IndexHeader) -> Self {
        Self {
            index: SaIndex::new(header),
            block: SaBlock::new(DEFAULT_BLOCK_SIZE),
            current_chrom: None,
            last_key: None,
            chrom_names: Vec::new(),
            data_offset: 0,
        }
    }

    /// Build .osa and .osa.idx from an iterator of sorted annotation records.
    ///
    /// Records MUST be sorted by (chrom_idx, position).
    /// `chrom_map` maps chrom_idx -> chromosome name string.
    pub fn write_all<W: Write>(
        &mut self,
        data_writer: &mut W,
        records: impl Iterator<Item = AnnotationRecord>,
        chrom_map: &[String],
    ) -> Result<()> {
        self.chrom_names = chrom_map.to_vec();

        // Write data file header
        data_writer.write_all(OSA_MAGIC)?;
        data_writer.write_all(&SCHEMA_VERSION.to_le_bytes())?;
        self.data_offset = (OSA_MAGIC.len() + 2) as u64;

        for record in records {
            self.write_record(data_writer, record, chrom_map)?;
        }

        // Flush remaining
        self.flush_block(data_writer)?;
        Ok(())
    }

    /// Build .osa and .osa.idx from an iterator that can surface parse errors.
    ///
    /// Records MUST be sorted by (chrom_idx, position).
    pub fn write_all_results<W: Write>(
        &mut self,
        data_writer: &mut W,
        records: impl Iterator<Item = Result<AnnotationRecord>>,
        chrom_map: &[String],
    ) -> Result<()> {
        self.chrom_names = chrom_map.to_vec();

        data_writer.write_all(OSA_MAGIC)?;
        data_writer.write_all(&SCHEMA_VERSION.to_le_bytes())?;
        self.data_offset = (OSA_MAGIC.len() + 2) as u64;

        for record in records {
            self.write_record(data_writer, record?, chrom_map)?;
        }

        self.flush_block(data_writer)?;
        Ok(())
    }

    fn write_record<W: Write>(
        &mut self,
        data_writer: &mut W,
        record: AnnotationRecord,
        chrom_map: &[String],
    ) -> Result<()> {
        if let Some((last_chrom, last_pos)) = self.last_key {
            if (record.chrom_idx, record.position) < (last_chrom, last_pos) {
                anyhow::bail!(
                    "SA records are not sorted: previous chrom_idx={}, position={}; current chrom_idx={}, position={}. \
                     The streaming .osa builder requires input sorted by chromosome (chr1..chr22,X,Y,M) then position \
                     — sort the source file (e.g. `bcftools sort` / `sort -k1,1 -k2,2n`) and rebuild.",
                    last_chrom,
                    last_pos,
                    record.chrom_idx,
                    record.position
                );
            }
        }
        self.last_key = Some((record.chrom_idx, record.position));

        let chrom_name = &chrom_map[record.chrom_idx as usize];

        // If we've moved to a new chromosome, flush the current block
        if self.current_chrom.as_ref() != Some(chrom_name) {
            self.flush_block(data_writer)?;
            self.current_chrom = Some(chrom_name.clone());
        }

        let entry = BlockEntry {
            position: record.position,
            ref_allele: record.ref_allele,
            alt_allele: record.alt_allele,
            json: record.json,
        };

        if !self.block.add(entry.clone()) {
            // Block is full, flush and retry
            self.flush_block(data_writer)?;
            assert!(self.block.add(entry), "Single entry exceeds block size");
        }

        Ok(())
    }

    fn flush_block<W: Write>(&mut self, writer: &mut W) -> Result<()> {
        if self.block.is_empty() {
            return Ok(());
        }

        let chrom = self.current_chrom.as_ref().unwrap().clone();
        let start_pos = self.block.start_position().unwrap();
        let end_pos = self.block.end_position().unwrap();

        let compressed = self.block.compress()?;
        let compressed_len = compressed.len() as u32;

        // Write compressed block length prefix + data
        writer.write_all(&compressed_len.to_le_bytes())?;
        writer.write_all(&compressed)?;

        self.index.add_block(
            &chrom,
            BlockRef {
                start_pos,
                end_pos,
                file_offset: self.data_offset,
                compressed_len,
            },
        );

        self.data_offset += 4 + compressed_len as u64;
        self.block.clear();
        Ok(())
    }

    /// Write the index file.
    pub fn write_index<W: Write>(&self, writer: &mut W) -> Result<()> {
        self.index.write_to(writer)
    }

    /// Convenience: write .osa and .osa.idx to files at the given base path.
    pub fn write_to_files(
        &mut self,
        base_path: &Path,
        records: impl Iterator<Item = AnnotationRecord>,
        chrom_map: &[String],
    ) -> Result<()> {
        let data_path = base_path.with_extension("osa");
        let idx_path = base_path.with_extension("osa.idx");

        let data_file = std::fs::File::create(&data_path)
            .with_context(|| format!("Creating output file {} (does the output directory exist?)", data_path.display()))?;
        let mut data_writer = BufWriter::new(data_file);
        self.write_all(&mut data_writer, records, chrom_map)?;
        data_writer.flush()?;

        let idx_file = std::fs::File::create(&idx_path)
            .with_context(|| format!("Creating index file {}", idx_path.display()))?;
        let mut idx_writer = BufWriter::new(idx_file);
        self.write_index(&mut idx_writer)?;
        idx_writer.flush()?;

        Ok(())
    }

    /// Convenience: write .osa and .osa.idx to files from fallible records.
    pub fn write_results_to_files(
        &mut self,
        base_path: &Path,
        records: impl Iterator<Item = Result<AnnotationRecord>>,
        chrom_map: &[String],
    ) -> Result<()> {
        let data_path = base_path.with_extension("osa");
        let idx_path = base_path.with_extension("osa.idx");

        let data_file = std::fs::File::create(&data_path)
            .with_context(|| format!("Creating output file {} (does the output directory exist?)", data_path.display()))?;
        let mut data_writer = BufWriter::new(data_file);
        self.write_all_results(&mut data_writer, records, chrom_map)?;
        data_writer.flush()?;

        let idx_file = std::fs::File::create(&idx_path)
            .with_context(|| format!("Creating index file {}", idx_path.display()))?;
        let mut idx_writer = BufWriter::new(idx_file);
        self.write_index(&mut idx_writer)?;
        idx_writer.flush()?;

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn header() -> IndexHeader {
        IndexHeader {
            schema_version: SCHEMA_VERSION,
            json_key: "test".into(),
            name: "Test".into(),
            version: "test".into(),
            description: "test".into(),
            assembly: "GRCh38".into(),
            match_by_allele: true,
            is_array: false,
            is_positional: false,
        }
    }

    fn record(chrom_idx: u16, position: u32) -> AnnotationRecord {
        AnnotationRecord {
            chrom_idx,
            position,
            ref_allele: "A".into(),
            alt_allele: "G".into(),
            json: "{}".into(),
        }
    }

    #[test]
    fn write_all_rejects_unsorted_records() {
        let mut writer = SaWriter::new(header());
        let mut out = Vec::new();
        let err = writer
            .write_all(
                &mut out,
                vec![record(0, 20), record(0, 10)].into_iter(),
                &["1".into()],
            )
            .unwrap_err();

        assert!(err.to_string().contains("SA records are not sorted"));
    }
}
