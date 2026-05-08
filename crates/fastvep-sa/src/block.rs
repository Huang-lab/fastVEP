//! Block-level compression for SA data.
//!
//! Each block stores annotations for a contiguous range of positions within a
//! chromosome. Items are serialized as length-prefixed entries, then the entire
//! block is zstd-compressed.

use crate::common::ZSTD_LEVEL;
use anyhow::Result;
use std::io::Read;

/// Hard cap on block decompressed size (256 MiB). Defends against zstd bombs
/// in maliciously crafted .osa files. Real blocks are typically <= 8 MiB.
const MAX_BLOCK_DECOMPRESSED: usize = 256 * 1024 * 1024;

/// A single entry within a block: position + ref + alt + json.
#[derive(Debug, Clone)]
pub struct BlockEntry {
    pub position: u32,
    pub ref_allele: String,
    pub alt_allele: String,
    pub json: String,
}

/// An in-memory block that accumulates entries and compresses them.
pub struct SaBlock {
    entries: Vec<BlockEntry>,
    uncompressed_size: usize,
    max_size: usize,
}

impl SaBlock {
    /// Create a new block with the given maximum uncompressed size.
    pub fn new(max_size: usize) -> Self {
        Self {
            entries: Vec::new(),
            uncompressed_size: 0,
            max_size,
        }
    }

    /// Try to add an entry. Returns `false` if the block is full (entry not added).
    pub fn add(&mut self, entry: BlockEntry) -> bool {
        let entry_size = 4 + 2 + entry.ref_allele.len() + 2 + entry.alt_allele.len() + 4 + entry.json.len();
        if !self.entries.is_empty() && self.uncompressed_size + entry_size > self.max_size {
            return false;
        }
        self.uncompressed_size += entry_size;
        self.entries.push(entry);
        true
    }

    /// Returns true if the block has no entries.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// The first position in this block.
    pub fn start_position(&self) -> Option<u32> {
        self.entries.first().map(|e| e.position)
    }

    /// The last position in this block.
    pub fn end_position(&self) -> Option<u32> {
        self.entries.last().map(|e| e.position)
    }

    /// Number of entries.
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// Serialize and compress the block. Entries are sorted by position before
    /// compression to enable binary search on decompressed data.
    pub fn compress(&self) -> Result<Vec<u8>> {
        // Sort entries by position for binary search support
        let mut sorted: Vec<usize> = (0..self.entries.len()).collect();
        sorted.sort_by_key(|&i| self.entries[i].position);

        let mut raw = Vec::with_capacity(self.uncompressed_size);

        // Write number of entries
        raw.extend_from_slice(&(self.entries.len() as u32).to_le_bytes());

        for &idx in &sorted {
            let entry = &self.entries[idx];
            // Position (4 bytes)
            raw.extend_from_slice(&entry.position.to_le_bytes());
            // Ref allele (2-byte length + data)
            raw.extend_from_slice(&(entry.ref_allele.len() as u16).to_le_bytes());
            raw.extend_from_slice(entry.ref_allele.as_bytes());
            // Alt allele (2-byte length + data)
            raw.extend_from_slice(&(entry.alt_allele.len() as u16).to_le_bytes());
            raw.extend_from_slice(entry.alt_allele.as_bytes());
            // JSON (4-byte length + data)
            raw.extend_from_slice(&(entry.json.len() as u32).to_le_bytes());
            raw.extend_from_slice(entry.json.as_bytes());
        }

        let compressed = zstd::encode_all(raw.as_slice(), ZSTD_LEVEL)?;
        Ok(compressed)
    }

    /// Decompress and deserialize a block from compressed bytes.
    pub fn decompress(data: &[u8]) -> Result<Vec<BlockEntry>> {
        // Streaming-decompress with a hard cap so a zstd bomb can never force
        // an oversized allocation. `decode_all` would have to materialize the
        // entire output before we could measure it.
        let mut decoder = zstd::stream::Decoder::new(data)?;
        let mut raw = Vec::new();
        (&mut decoder)
            .take(MAX_BLOCK_DECOMPRESSED as u64 + 1)
            .read_to_end(&mut raw)?;
        if raw.len() > MAX_BLOCK_DECOMPRESSED {
            anyhow::bail!(
                "Decompressed block exceeds limit ({} bytes)",
                MAX_BLOCK_DECOMPRESSED
            );
        }
        let mut cursor: usize = 0;

        // Helper: ensure `n` more bytes are available starting at cursor.
        let need = |cursor: usize, n: usize, raw: &[u8]| -> Result<()> {
            let end = cursor
                .checked_add(n)
                .ok_or_else(|| anyhow::anyhow!("Block cursor overflow"))?;
            if end > raw.len() {
                anyhow::bail!("Unexpected end of block data");
            }
            Ok(())
        };

        need(cursor, 4, &raw)?;
        let count = u32::from_le_bytes(raw[cursor..cursor + 4].try_into()?) as usize;
        cursor += 4;

        // Each entry requires at least 12 bytes after the count field:
        // 4 position + 2 ref_len + 0 ref + 2 alt_len + 0 alt + 4 json_len + 0 json.
        // Validate against the bytes remaining in the block, and do not use the
        // untrusted count for large upfront allocation.
        let remaining = raw.len() - cursor;
        if count > remaining / 12 {
            anyhow::bail!("Block claims {} entries, exceeds data size", count);
        }

        let mut entries = Vec::new();
        for _ in 0..count {
            // Position
            need(cursor, 4, &raw)?;
            let position = u32::from_le_bytes(raw[cursor..cursor + 4].try_into()?);
            cursor += 4;

            // Ref allele
            need(cursor, 2, &raw)?;
            let ref_len = u16::from_le_bytes(raw[cursor..cursor + 2].try_into()?) as usize;
            cursor += 2;
            need(cursor, ref_len, &raw)?;
            let ref_allele = std::str::from_utf8(&raw[cursor..cursor + ref_len])?.to_string();
            cursor += ref_len;

            // Alt allele
            need(cursor, 2, &raw)?;
            let alt_len = u16::from_le_bytes(raw[cursor..cursor + 2].try_into()?) as usize;
            cursor += 2;
            need(cursor, alt_len, &raw)?;
            let alt_allele = std::str::from_utf8(&raw[cursor..cursor + alt_len])?.to_string();
            cursor += alt_len;

            // JSON
            need(cursor, 4, &raw)?;
            let json_len = u32::from_le_bytes(raw[cursor..cursor + 4].try_into()?) as usize;
            cursor += 4;
            need(cursor, json_len, &raw)?;
            let json = std::str::from_utf8(&raw[cursor..cursor + json_len])?.to_string();
            cursor += json_len;

            entries.push(BlockEntry {
                position,
                ref_allele,
                alt_allele,
                json,
            });
        }

        Ok(entries)
    }

    /// Binary search for a variant in sorted decompressed entries.
    ///
    /// Uses Var32 encoding to find the entry index, then verifies the full
    /// allele match (handles long variants that can't be Var32-encoded).
    pub fn find_by_position(
        entries: &[BlockEntry],
        position: u32,
        ref_allele: &str,
        alt_allele: &str,
        is_positional: bool,
    ) -> Option<usize> {
        if is_positional {
            // Positional: binary search by position only
            let idx = entries.partition_point(|e| e.position < position);
            if idx < entries.len() && entries[idx].position == position {
                return Some(idx);
            }
            return None;
        }

        // Find the first entry at this position via binary search
        let start = entries.partition_point(|e| e.position < position);

        // Linear scan among entries at the same position (usually just 1-3)
        for i in start..entries.len() {
            if entries[i].position != position {
                break;
            }
            if entries[i].ref_allele == ref_allele && entries[i].alt_allele == alt_allele {
                return Some(i);
            }
        }
        None
    }

    /// Reset the block for reuse.
    pub fn clear(&mut self) {
        self.entries.clear();
        self.uncompressed_size = 0;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_block_round_trip() {
        let mut block = SaBlock::new(1024 * 1024);
        for i in 0..100 {
            assert!(block.add(BlockEntry {
                position: 1000 + i,
                ref_allele: "A".into(),
                alt_allele: "G".into(),
                json: format!(r#"{{"score":{}}}"#, i),
            }));
        }

        assert_eq!(block.len(), 100);
        assert_eq!(block.start_position(), Some(1000));
        assert_eq!(block.end_position(), Some(1099));

        let compressed = block.compress().unwrap();
        let entries = SaBlock::decompress(&compressed).unwrap();
        assert_eq!(entries.len(), 100);
        assert_eq!(entries[0].position, 1000);
        assert_eq!(entries[0].ref_allele, "A");
        assert_eq!(entries[0].alt_allele, "G");
        assert_eq!(entries[0].json, r#"{"score":0}"#);
        assert_eq!(entries[99].position, 1099);
        assert_eq!(entries[99].json, r#"{"score":99}"#);
    }

    #[test]
    fn test_block_full() {
        // Each entry is roughly: 4 + 2+1 + 2+1 + 4+7 = 21 bytes.
        // Set max_size to 25 so first fits, second doesn't.
        let mut block = SaBlock::new(25);
        // First entry always accepted (even if > max_size)
        assert!(block.add(BlockEntry {
            position: 1,
            ref_allele: "A".into(),
            alt_allele: "G".into(),
            json: r#"{"x":1}"#.into(),
        }));
        // Second should be rejected (block full: 21 + 21 = 42 > 25)
        assert!(!block.add(BlockEntry {
            position: 2,
            ref_allele: "A".into(),
            alt_allele: "G".into(),
            json: r#"{"x":2}"#.into(),
        }));
    }

    #[test]
    fn test_decompress_truncated_returns_error() {
        // Build a valid block, then truncate the compressed payload to feed
        // the decompressor short input. It must error rather than panic.
        let mut block = SaBlock::new(1024);
        block.add(BlockEntry {
            position: 1,
            ref_allele: "ACGT".into(),
            alt_allele: "T".into(),
            json: r#"{"x":1}"#.into(),
        });
        let compressed = block.compress().unwrap();

        // Decompress the full block once to capture the raw layout, then
        // construct a truncated raw buffer and re-compress it so the
        // bounds-check path inside decompress() is exercised on each missing
        // length field.
        let raw = zstd::decode_all(compressed.as_slice()).unwrap();
        for cut in 0..raw.len() {
            let truncated = zstd::encode_all(&raw[..cut], 3).unwrap();
            // Either Ok with empty entries, or Err — never a panic.
            let _ = SaBlock::decompress(&truncated);
        }
    }

    #[test]
    fn test_decompress_lying_count_rejected() {
        // Hand-craft a small zstd-compressed buffer that claims a huge entry
        // count but supplies no entry data. Must error.
        let mut raw = Vec::new();
        raw.extend_from_slice(&u32::MAX.to_le_bytes());
        let compressed = zstd::encode_all(raw.as_slice(), 3).unwrap();
        let result = SaBlock::decompress(&compressed);
        assert!(result.is_err(), "expected error for absurd entry count");
    }

    #[test]
    fn test_block_clear() {
        let mut block = SaBlock::new(1024);
        block.add(BlockEntry {
            position: 1,
            ref_allele: "A".into(),
            alt_allele: "G".into(),
            json: "{}".into(),
        });
        assert!(!block.is_empty());
        block.clear();
        assert!(block.is_empty());
        assert_eq!(block.len(), 0);
    }
}
