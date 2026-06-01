//! Index structures for fastSA files.
//!
//! The index maps (chromosome, position) -> file offset so the reader can
//! seek directly to the relevant compressed block.

use crate::common::{chrom_aliases, MAX_INDEX_PAYLOAD, OSA_MAGIC, SCHEMA_VERSION};
use anyhow::Result;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::io::{Read, Write};

/// A reference to a single compressed block within the data file.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BlockRef {
    /// First position covered by this block.
    pub start_pos: u32,
    /// Last position covered by this block.
    pub end_pos: u32,
    /// Byte offset in the data file where the compressed block starts.
    pub file_offset: u64,
    /// Length of the compressed block in bytes.
    pub compressed_len: u32,
}

/// Metadata stored in the index header.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct IndexHeader {
    pub schema_version: u16,
    pub json_key: String,
    pub name: String,
    pub version: String,
    pub description: String,
    pub assembly: String,
    pub match_by_allele: bool,
    pub is_array: bool,
    pub is_positional: bool,
}

/// The complete index for an OSA file.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SaIndex {
    pub header: IndexHeader,
    /// Chromosome name -> list of block references, sorted by start_pos.
    pub chromosomes: HashMap<String, Vec<BlockRef>>,
}

impl SaIndex {
    /// Create a new empty index with the given header.
    pub fn new(header: IndexHeader) -> Self {
        Self {
            header,
            chromosomes: HashMap::new(),
        }
    }

    /// Add a block reference for a chromosome.
    pub fn add_block(&mut self, chrom: &str, block_ref: BlockRef) {
        self.chromosomes
            .entry(chrom.to_string())
            .or_default()
            .push(block_ref);
    }

    /// Resolve a chromosome name to the on-disk key, tolerating `chr*` vs
    /// bare naming and mitochondrial aliases. Returns `None` if no alias is
    /// present in the index.
    fn blocks_for_chrom(&self, chrom: &str) -> Option<&Vec<BlockRef>> {
        for alias in chrom_aliases(chrom) {
            if let Some(blocks) = self.chromosomes.get(&alias) {
                return Some(blocks);
            }
        }
        None
    }

    /// Find the block(s) that may contain the given position on a chromosome.
    /// Returns the file offset and compressed length of each matching block.
    pub fn find_blocks(&self, chrom: &str, position: u32) -> Vec<&BlockRef> {
        let blocks = match self.blocks_for_chrom(chrom) {
            Some(b) => b,
            None => return Vec::new(),
        };

        // Binary search for the first block whose end_pos >= position
        let start_idx = blocks.partition_point(|b| b.end_pos < position);

        // Collect blocks that overlap [position, position]
        let mut result = Vec::new();
        for block in &blocks[start_idx..] {
            if block.start_pos > position {
                break;
            }
            result.push(block);
        }
        result
    }

    /// Find all blocks that overlap the given range [start, end].
    pub fn find_blocks_range(&self, chrom: &str, start: u32, end: u32) -> Vec<&BlockRef> {
        let blocks = match self.blocks_for_chrom(chrom) {
            Some(b) => b,
            None => return Vec::new(),
        };

        let start_idx = blocks.partition_point(|b| b.end_pos < start);

        let mut result = Vec::new();
        for block in &blocks[start_idx..] {
            if block.start_pos > end {
                break;
            }
            result.push(block);
        }
        result
    }

    /// Serialize the index to a writer.
    pub fn write_to<W: Write>(&self, writer: &mut W) -> Result<()> {
        // Magic + schema version
        writer.write_all(OSA_MAGIC)?;
        writer.write_all(&SCHEMA_VERSION.to_le_bytes())?;

        // Serialize the rest with bincode
        let data = bincode::serialize(self)?;
        writer.write_all(&(data.len() as u64).to_le_bytes())?;
        writer.write_all(&data)?;
        Ok(())
    }

    /// Deserialize the index from a reader.
    pub fn read_from<R: Read>(reader: &mut R) -> Result<Self> {
        // Verify magic
        let mut magic = [0u8; 8];
        reader.read_exact(&mut magic)?;
        if &magic != OSA_MAGIC {
            anyhow::bail!(
                "Invalid OSA index magic: expected {:?}, got {:?}",
                OSA_MAGIC,
                magic
            );
        }

        // Verify schema version
        let mut version_bytes = [0u8; 2];
        reader.read_exact(&mut version_bytes)?;
        let version = u16::from_le_bytes(version_bytes);
        if version != SCHEMA_VERSION {
            anyhow::bail!(
                "Unsupported schema version: expected {}, got {}",
                SCHEMA_VERSION,
                version
            );
        }

        // Read bincode data
        let mut len_bytes = [0u8; 8];
        reader.read_exact(&mut len_bytes)?;
        let len_u64 = u64::from_le_bytes(len_bytes);
        if len_u64 > MAX_INDEX_PAYLOAD {
            anyhow::bail!(
                "Index payload size {} exceeds limit {}",
                len_u64,
                MAX_INDEX_PAYLOAD
            );
        }
        let len: usize = len_u64
            .try_into()
            .map_err(|_| anyhow::anyhow!("Index payload size {} exceeds usize", len_u64))?;

        let mut data = vec![0u8; len];
        reader.read_exact(&mut data)?;

        let index: SaIndex = bincode::deserialize(&data)?;
        Ok(index)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_index_round_trip() {
        let header = IndexHeader {
            schema_version: SCHEMA_VERSION,
            json_key: "test".into(),
            name: "Test Source".into(),
            version: "1.0".into(),
            description: "Test".into(),
            assembly: "GRCh38".into(),
            match_by_allele: true,
            is_array: false,
            is_positional: false,
        };

        let mut index = SaIndex::new(header);
        index.add_block("chr1", BlockRef {
            start_pos: 100,
            end_pos: 500,
            file_offset: 0,
            compressed_len: 1024,
        });
        index.add_block("chr1", BlockRef {
            start_pos: 501,
            end_pos: 1000,
            file_offset: 1024,
            compressed_len: 2048,
        });
        index.add_block("chr2", BlockRef {
            start_pos: 1,
            end_pos: 300,
            file_offset: 3072,
            compressed_len: 512,
        });

        // Serialize and deserialize
        let mut buf = Vec::new();
        index.write_to(&mut buf).unwrap();

        let mut cursor = std::io::Cursor::new(buf);
        let loaded = SaIndex::read_from(&mut cursor).unwrap();

        assert_eq!(loaded.header.json_key, "test");
        assert_eq!(loaded.chromosomes.len(), 2);
        assert_eq!(loaded.chromosomes["chr1"].len(), 2);
        assert_eq!(loaded.chromosomes["chr2"].len(), 1);
    }

    #[test]
    fn test_find_blocks() {
        let header = IndexHeader {
            schema_version: SCHEMA_VERSION,
            json_key: "test".into(),
            name: "Test".into(),
            version: "1.0".into(),
            description: "".into(),
            assembly: "GRCh38".into(),
            match_by_allele: true,
            is_array: false,
            is_positional: false,
        };

        let mut index = SaIndex::new(header);
        index.add_block("chr1", BlockRef {
            start_pos: 100,
            end_pos: 500,
            file_offset: 0,
            compressed_len: 100,
        });
        index.add_block("chr1", BlockRef {
            start_pos: 501,
            end_pos: 1000,
            file_offset: 100,
            compressed_len: 200,
        });

        // Position in first block
        let blocks = index.find_blocks("chr1", 250);
        assert_eq!(blocks.len(), 1);
        assert_eq!(blocks[0].start_pos, 100);

        // Position in second block
        let blocks = index.find_blocks("chr1", 750);
        assert_eq!(blocks.len(), 1);
        assert_eq!(blocks[0].start_pos, 501);

        // Position before any block
        let blocks = index.find_blocks("chr1", 50);
        assert_eq!(blocks.len(), 0);

        // Position after all blocks
        let blocks = index.find_blocks("chr1", 1500);
        assert_eq!(blocks.len(), 0);

        // Missing chromosome
        let blocks = index.find_blocks("chr99", 100);
        assert_eq!(blocks.len(), 0);
    }
}
