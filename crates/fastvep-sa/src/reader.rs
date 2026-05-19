//! Reader for .osa position/allele-level annotation files.
//!
//! Uses memory-mapped I/O for the data file and binary search on the index
//! for O(log n) block lookups. Decompressed blocks are held in a thread-safe
//! LRU cache shared across batches and across queries on the same block.

use crate::block::{BlockEntry, SaBlock};
use crate::common::OSA_MAGIC;
use crate::index::{BlockRef, SaIndex};
use anyhow::Result;
use fastvep_cache::annotation::{AnnotationProvider, AnnotationValue, SaMetadata};
use lru::LruCache;
use memmap2::Mmap;
use std::fs::File;
use std::num::NonZeroUsize;
use std::path::Path;
use std::sync::{Arc, Mutex};

/// Number of decompressed blocks retained in the LRU cache.
///
/// VCF inputs are sorted by (chrom, pos) and a batch of 1024 variants typically
/// spans only a handful of blocks, but the per-transcript/per-allele query
/// pattern hits each block many times. 64 blocks is enough to keep the
/// working set hot across sorted-input batches without bounding RSS too high.
const BLOCK_CACHE_CAPACITY: usize = 64;

/// Reader for .osa annotation files.
///
/// Thread-safety: the reader is `Send + Sync`. Block decompression results are
/// cached in a `Mutex<LruCache>`; the lock is held only for the brief lookup
/// or insert and is released across the decompression and find_match steps.
pub struct SaReader {
    mmap: Mmap,
    index: SaIndex,
    metadata: SaMetadata,
    /// LRU cache of decompressed blocks, keyed by file offset.
    ///
    /// `Arc` lets a worker clone a reference to the block payload and drop the
    /// mutex before searching, so other workers can hit the cache concurrently.
    block_cache: Mutex<LruCache<u64, Arc<Vec<BlockEntry>>>>,
}

impl SaReader {
    /// Open an .osa + .osa.idx file pair.
    pub fn open(data_path: &Path) -> Result<Self> {
        let idx_path = data_path.with_extension("osa.idx");

        let mut idx_file = File::open(&idx_path)?;
        let index = SaIndex::read_from(&mut idx_file)?;

        let data_file = File::open(data_path)?;
        let mmap = unsafe { Mmap::map(&data_file)? };

        if mmap.len() < 10 || &mmap[..8] != OSA_MAGIC {
            anyhow::bail!("Invalid OSA data file: bad magic");
        }

        let metadata = SaMetadata {
            name: index.header.name.clone(),
            version: index.header.version.clone(),
            description: index.header.description.clone(),
            assembly: index.header.assembly.clone(),
            json_key: index.header.json_key.clone(),
            match_by_allele: index.header.match_by_allele,
            is_array: index.header.is_array,
            is_positional: index.header.is_positional,
        };

        let capacity = NonZeroUsize::new(BLOCK_CACHE_CAPACITY)
            .expect("BLOCK_CACHE_CAPACITY is a non-zero compile-time constant");

        Ok(Self {
            mmap,
            index,
            metadata,
            block_cache: Mutex::new(LruCache::new(capacity)),
        })
    }

    /// Decompress a block straight from the mmap. Pure: touches no cache state.
    fn decompress_block(&self, file_offset: u64, compressed_len: u32) -> Result<Vec<BlockEntry>> {
        let offset: usize = file_offset
            .try_into()
            .map_err(|_| anyhow::anyhow!("Block offset {} too large for usize", file_offset))?;
        // Data file layout per block: [4-byte compressed_len] [compressed_data]
        let data_start = offset
            .checked_add(4)
            .ok_or_else(|| anyhow::anyhow!("Block offset overflow"))?;
        let data_end = data_start
            .checked_add(compressed_len as usize)
            .ok_or_else(|| anyhow::anyhow!("Block end offset overflow"))?;

        if data_end > self.mmap.len() {
            anyhow::bail!("Block extends beyond data file");
        }

        SaBlock::decompress(&self.mmap[data_start..data_end])
    }

    /// Return the decompressed block at the given file offset, hitting or
    /// populating the LRU cache as needed.
    fn get_block(&self, block_ref: &BlockRef) -> Result<Arc<Vec<BlockEntry>>> {
        // Fast path: cache hit.
        {
            let mut cache = self
                .block_cache
                .lock()
                .map_err(|_| anyhow::anyhow!("SA block cache mutex poisoned"))?;
            if let Some(arc) = cache.get(&block_ref.file_offset) {
                return Ok(Arc::clone(arc));
            }
        }

        // Slow path: decompress without holding the lock so other workers can
        // serve their own queries from the cache concurrently. If two threads
        // race on the same missing block they each decompress once; the second
        // `put` simply replaces an identical entry — acceptable for an LRU.
        let entries = self.decompress_block(block_ref.file_offset, block_ref.compressed_len)?;
        let arc = Arc::new(entries);

        let mut cache = self
            .block_cache
            .lock()
            .map_err(|_| anyhow::anyhow!("SA block cache mutex poisoned"))?;
        cache.put(block_ref.file_offset, Arc::clone(&arc));
        Ok(arc)
    }

    /// Query annotations for a specific position and allele.
    fn query(
        &self,
        chrom: &str,
        position: u32,
        ref_allele: &str,
        alt_allele: &str,
    ) -> Result<Option<String>> {
        let block_refs = self.index.find_blocks(chrom, position);
        for block_ref in block_refs {
            let entries = self.get_block(block_ref)?;
            if let Some(json) = self.find_match(&entries, position, ref_allele, alt_allele) {
                return Ok(Some(json));
            }
        }
        Ok(None)
    }

    fn find_match(
        &self,
        entries: &[BlockEntry],
        position: u32,
        ref_allele: &str,
        alt_allele: &str,
    ) -> Option<String> {
        let allele_ref = if self.metadata.match_by_allele { ref_allele } else { "" };
        let allele_alt = if self.metadata.match_by_allele { alt_allele } else { "" };

        SaBlock::find_by_position(
            entries,
            position,
            allele_ref,
            allele_alt,
            self.metadata.is_positional,
        )
        .map(|idx| entries[idx].json.clone())
    }
}

impl AnnotationProvider for SaReader {
    fn name(&self) -> &str {
        &self.metadata.name
    }

    fn json_key(&self) -> &str {
        &self.metadata.json_key
    }

    fn metadata(&self) -> &SaMetadata {
        &self.metadata
    }

    fn annotate_position(
        &self,
        chrom: &str,
        pos: u64,
        ref_allele: &str,
        alt_allele: &str,
    ) -> Result<Option<AnnotationValue>> {
        let position: u32 = pos
            .try_into()
            .map_err(|_| anyhow::anyhow!("Position {} exceeds u32::MAX", pos))?;
        match self.query(chrom, position, ref_allele, alt_allele)? {
            Some(json) => {
                if self.metadata.is_positional {
                    Ok(Some(AnnotationValue::Positional(json)))
                } else {
                    Ok(Some(AnnotationValue::Json(json)))
                }
            }
            None => Ok(None),
        }
    }

    /// Decompress (and cache) the blocks containing each requested position.
    ///
    /// Unlike a range-based preload, this only touches blocks that actually
    /// hold at least one queried position, so a batch that straddles a wide
    /// region but lands in only a few blocks does not pay for everything in
    /// between. Already-cached blocks are no-ops.
    fn preload(&self, chrom: &str, positions: &[u64]) -> Result<()> {
        if positions.is_empty() {
            return Ok(());
        }

        let blocks = match self.index.chromosomes.get(chrom) {
            Some(b) => b.as_slice(),
            None => return Ok(()),
        };
        if blocks.is_empty() {
            return Ok(());
        }

        // Sort + dedup positions so the sweep across blocks is monotonic.
        let max_u32 = u32::MAX as u64;
        let mut positions_u32: Vec<u32> = Vec::with_capacity(positions.len());
        for &p in positions {
            if p > max_u32 {
                anyhow::bail!("Position {} exceeds u32::MAX", p);
            }
            positions_u32.push(p as u32);
        }
        positions_u32.sort_unstable();
        positions_u32.dedup();

        // Single forward pass: for each position, advance to the first block
        // whose end >= pos; if that block also starts <= pos, decompress it
        // (once per offset). Blocks are sorted by start_pos.
        let mut block_idx = 0usize;
        let mut last_loaded: Option<u64> = None;
        for &pos in &positions_u32 {
            while block_idx < blocks.len() && blocks[block_idx].end_pos < pos {
                block_idx += 1;
            }
            if block_idx >= blocks.len() {
                break;
            }
            let block_ref = &blocks[block_idx];
            if block_ref.start_pos > pos {
                continue; // position falls in a gap between blocks
            }
            if last_loaded == Some(block_ref.file_offset) {
                continue; // multiple positions inside the same block
            }
            self.get_block(block_ref)?;
            last_loaded = Some(block_ref.file_offset);
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::{AnnotationRecord, SCHEMA_VERSION};
    use crate::index::IndexHeader;
    use crate::writer::SaWriter;
    use tempfile::TempDir;

    fn header(match_by_allele: bool, is_positional: bool) -> IndexHeader {
        IndexHeader {
            schema_version: SCHEMA_VERSION,
            json_key: "test".into(),
            name: "Test".into(),
            version: "1.0".into(),
            description: "".into(),
            assembly: "GRCh38".into(),
            match_by_allele,
            is_array: false,
            is_positional,
        }
    }

    fn write_fixture(path: &Path, records: Vec<AnnotationRecord>) {
        let chrom_map = vec!["chr1".to_string()];
        let mut writer = SaWriter::new(header(true, false));
        writer
            .write_to_files(path, records.into_iter(), &chrom_map)
            .unwrap();
    }

    #[test]
    fn query_roundtrip_via_block_cache() {
        let dir = TempDir::new().unwrap();
        let base = dir.path().join("test");
        write_fixture(
            &base,
            (0..100)
                .map(|i| AnnotationRecord {
                    chrom_idx: 0,
                    position: 1000 + i,
                    ref_allele: "A".into(),
                    alt_allele: "G".into(),
                    json: format!(r#"{{"i":{}}}"#, i),
                })
                .collect(),
        );

        let reader = SaReader::open(&base.with_extension("osa")).unwrap();
        let ann = reader
            .annotate_position("chr1", 1042, "A", "G")
            .unwrap()
            .unwrap();
        match ann {
            AnnotationValue::Json(j) => assert!(j.contains(r#""i":42"#)),
            other => panic!("expected JSON value, got {:?}", other),
        }

        // Cache hit on second query of same block — exercises the fast path.
        let again = reader
            .annotate_position("chr1", 1043, "A", "G")
            .unwrap()
            .unwrap();
        match again {
            AnnotationValue::Json(j) => assert!(j.contains(r#""i":43"#)),
            _ => unreachable!(),
        }
    }

    #[test]
    fn preload_only_touches_blocks_containing_queried_positions() {
        let dir = TempDir::new().unwrap();
        let base = dir.path().join("test");
        // Records spaced far enough to force multiple blocks (block sized
        // generously by default; here we generate enough variants that several
        // blocks are flushed).
        let records: Vec<AnnotationRecord> = (0..200_000)
            .map(|i| AnnotationRecord {
                chrom_idx: 0,
                position: 1000 + i,
                ref_allele: "A".into(),
                alt_allele: "G".into(),
                json: r#"{}"#.into(),
            })
            .collect();
        write_fixture(&base, records);

        let reader = SaReader::open(&base.with_extension("osa")).unwrap();
        // Preload a single position; only the containing block should warm up.
        reader.preload("chr1", &[1042]).unwrap();
        let ann = reader
            .annotate_position("chr1", 1042, "A", "G")
            .unwrap();
        assert!(ann.is_some());

        // Unknown chromosome must be a no-op rather than an error.
        reader.preload("chrUnknown", &[1, 2, 3]).unwrap();
    }

    #[test]
    fn missing_position_returns_none() {
        let dir = TempDir::new().unwrap();
        let base = dir.path().join("test");
        write_fixture(
            &base,
            vec![AnnotationRecord {
                chrom_idx: 0,
                position: 100,
                ref_allele: "A".into(),
                alt_allele: "G".into(),
                json: "{}".into(),
            }],
        );

        let reader = SaReader::open(&base.with_extension("osa")).unwrap();
        assert!(reader
            .annotate_position("chr1", 200, "A", "G")
            .unwrap()
            .is_none());
        assert!(reader
            .annotate_position("chr2", 100, "A", "G")
            .unwrap()
            .is_none());
    }
}
