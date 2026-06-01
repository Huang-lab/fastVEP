//! Reader for .osa position/allele-level annotation files.
//!
//! Uses memory-mapped I/O for the data file and binary search on the index
//! for O(log n) block lookups. Decompressed blocks are held in a thread-safe
//! byte-budgeted LRU cache shared across batches and across queries on the
//! same block.

use crate::block::{BlockEntry, SaBlock};
use crate::common::{chrom_aliases, OSA_MAGIC};
use crate::index::{BlockRef, SaIndex};
use anyhow::Result;
use fastvep_cache::annotation::{AnnotationProvider, AnnotationValue, SaMetadata};
use lru::LruCache;
use memmap2::Mmap;
use std::fs::File;
use std::num::NonZeroUsize;
use std::path::Path;
use std::sync::{Arc, Mutex};

/// Default per-reader byte budget for the decompressed-block cache (32 MiB).
///
/// Block sizes vary a lot by data source — clinvar/gnomad/spliceai are
/// ~10 MiB per block once decompressed, REVEL is ~25 MiB, PhyloP can be
/// ~40 MiB. A fixed *count* cap (e.g. 4 blocks) inflates total RSS by an
/// order of magnitude on PhyloP/REVEL-heavy stacks and OOMs on
/// full-genome inputs. A byte budget adapts: low-density readers cache
/// 2–3 blocks, high-density readers cache 1.
///
/// With ~100 readers in the full SA stack (one per chrom × DB), 32 MiB
/// caps cache memory at roughly 3.2 GiB worst case, leaving headroom on
/// a 12 GiB sandbox after the GFF3 cache + indexes (~3.5 GiB baseline).
/// Override via `FASTVEP_SA_CACHE_BYTES_PER_READER` (in bytes). The cache
/// is guaranteed to retain at least 1 block to avoid thrashing on a
/// single just-decompressed block under parallel queries.
const DEFAULT_CACHE_BYTES_PER_READER: usize = 32 * 1024 * 1024;

/// Soft upper bound on entries to prevent the underlying `LruCache`'s
/// capacity field from being a pathological size if blocks ever ended up
/// being tiny. The byte budget is the real gate.
const CACHE_MAX_ENTRIES: usize = 1024;

fn cache_bytes_per_reader() -> usize {
    std::env::var("FASTVEP_SA_CACHE_BYTES_PER_READER")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .filter(|&v| v > 0)
        .unwrap_or(DEFAULT_CACHE_BYTES_PER_READER)
}

/// Approximate in-memory footprint of a decompressed block: the BlockEntry
/// struct slots in the `Vec`, plus the heap storage backing each entry's
/// three `String`s. The `Vec` capacity is bounded by its length here
/// because the writer pre-sizes the allocation, so `len * size_of` is a
/// reasonable proxy for the slab.
fn block_bytes(entries: &[BlockEntry]) -> usize {
    let slot_bytes = std::mem::size_of::<BlockEntry>().saturating_mul(entries.len());
    let string_bytes: usize = entries
        .iter()
        .map(|e| e.ref_allele.len() + e.alt_allele.len() + e.json.len())
        .sum();
    slot_bytes.saturating_add(string_bytes)
}

/// LRU keyed by file_offset, with byte-based eviction. Evicts least-recently-
/// used entries until total cached bytes is within budget. Always keeps at
/// least one entry (the just-inserted one) so a single oversized block
/// doesn't keep falling out from under a parallel-worker batch.
///
/// Byte accounting goes through `pop_lru` (explicit) and `push` (which
/// returns the evicted entry on capacity overflow), so the inner `LruCache`
/// can never silently drop an entry without `total_bytes` reflecting it.
struct BlockCache {
    lru: LruCache<u64, (Arc<Vec<BlockEntry>>, usize)>,
    total_bytes: usize,
    budget_bytes: usize,
}

impl BlockCache {
    fn new(budget_bytes: usize) -> Self {
        let cap = NonZeroUsize::new(CACHE_MAX_ENTRIES).expect("non-zero");
        Self {
            lru: LruCache::new(cap),
            total_bytes: 0,
            budget_bytes,
        }
    }

    fn get(&mut self, offset: u64) -> Option<Arc<Vec<BlockEntry>>> {
        self.lru.get(&offset).map(|(arc, _)| Arc::clone(arc))
    }

    fn put(&mut self, offset: u64, value: Arc<Vec<BlockEntry>>, bytes: usize) {
        // Replace-in-place: drop the old bytes first so the budget loop
        // below sees the correct `total_bytes`.
        if let Some((_, old_bytes)) = self.lru.pop(&offset) {
            self.total_bytes = self.total_bytes.saturating_sub(old_bytes);
        }
        // Byte-budget eviction: free space for the new block, but always
        // leave room to insert it (cache may go above budget for a single
        // oversized block; parallel workers must not thrash that one).
        while self.total_bytes + bytes > self.budget_bytes && !self.lru.is_empty() {
            if let Some((_, (_, ev_bytes))) = self.lru.pop_lru() {
                self.total_bytes = self.total_bytes.saturating_sub(ev_bytes);
            } else {
                break;
            }
        }
        // `push` returns any entry the inner LruCache evicts to make room
        // (capacity overflow); subtract its bytes so `total_bytes` doesn't
        // drift if `CACHE_MAX_ENTRIES` ever bites before the byte budget.
        if let Some((_, (_, ev_bytes))) = self.lru.push(offset, (value, bytes)) {
            self.total_bytes = self.total_bytes.saturating_sub(ev_bytes);
        }
        self.total_bytes = self.total_bytes.saturating_add(bytes);
    }

    #[cfg(test)]
    fn len(&self) -> usize {
        self.lru.len()
    }
}

/// Reader for .osa annotation files.
///
/// Thread-safety: the reader is `Send + Sync`. Block decompression results are
/// cached in a `Mutex<LruCache>`; the lock is held only for the brief lookup
/// or insert and is released across the decompression and find_match steps.
pub struct SaReader {
    mmap: Mmap,
    index: SaIndex,
    metadata: SaMetadata,
    /// Byte-budgeted LRU cache of decompressed blocks, keyed by file offset.
    ///
    /// `Arc` lets a worker clone a reference to the block payload and drop the
    /// mutex before searching, so other workers can hit the cache concurrently.
    block_cache: Mutex<BlockCache>,
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

        Ok(Self {
            mmap,
            index,
            metadata,
            block_cache: Mutex::new(BlockCache::new(cache_bytes_per_reader())),
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

        // Cross-check the on-disk length prefix against the index. If they
        // disagree the `.osa` and `.osa.idx` are out of sync (corrupt or
        // mismatched files) and we'd otherwise silently decompress the wrong
        // byte range.
        // Bounds were just verified, but never `.expect()` on parsed bytes:
        // surface any unexpected slice shape as a typed error so debugging a
        // mismatched .osa/.osa.idx pair never produces a panic.
        let len_bytes: [u8; 4] = self.mmap[offset..offset + 4]
            .try_into()
            .map_err(|_| anyhow::anyhow!("expected 4-byte block length prefix at offset {}", offset))?;
        let on_disk_len = u32::from_le_bytes(len_bytes);
        if on_disk_len != compressed_len {
            anyhow::bail!(
                "Block length mismatch at offset {}: index says {} bytes, data file prefix says {}",
                file_offset,
                compressed_len,
                on_disk_len,
            );
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
            if let Some(arc) = cache.get(block_ref.file_offset) {
                return Ok(arc);
            }
        }

        // Slow path: decompress without holding the lock so other workers can
        // serve their own queries from the cache concurrently. If two threads
        // race on the same missing block they each decompress once; the second
        // `put` simply replaces an identical entry — acceptable for an LRU.
        let entries = self.decompress_block(block_ref.file_offset, block_ref.compressed_len)?;
        let bytes = block_bytes(&entries);
        let arc = Arc::new(entries);

        let mut cache = self
            .block_cache
            .lock()
            .map_err(|_| anyhow::anyhow!("SA block cache mutex poisoned"))?;
        cache.put(block_ref.file_offset, Arc::clone(&arc), bytes);
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

        // Honor the same chr*/bare/mitochondrial aliases as `find_blocks`
        // so a preload on `chr1` against an index built with `1` (or vice
        // versa) still primes the cache instead of silently no-op'ing.
        let blocks = chrom_aliases(chrom)
            .iter()
            .find_map(|alias| self.index.chromosomes.get(alias))
            .map(|v| v.as_slice());
        let blocks = match blocks {
            Some(b) => b,
            None => {
                // A chromosome the caller asked about that is not in the
                // index isn't necessarily an error (e.g., chrM absent from
                // ClinVar), but a typo would otherwise produce silently
                // empty annotations forever. Surface it at debug level so
                // operators can grep their logs without drowning in noise
                // on normal runs.
                log::debug!(
                    "SA preload: chromosome '{}' (and aliases) not present in {} index",
                    chrom,
                    self.metadata.name
                );
                return Ok(());
            }
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
        // Use a JSON payload big enough that the writer flushes multiple
        // 8 MiB blocks (entry_size accounting in SaBlock::add is
        // 4 + 2 + ref + 2 + alt + 4 + json bytes ≈ 13 + 200 ≈ 213 B; at
        // 100_000 entries that's ~21 MiB → at least 3 blocks).
        let big_json = "x".repeat(200);
        let records: Vec<AnnotationRecord> = (0..100_000)
            .map(|i| AnnotationRecord {
                chrom_idx: 0,
                position: 1000 + i,
                ref_allele: "A".into(),
                alt_allele: "G".into(),
                json: format!(r#"{{"i":{},"pad":"{}"}}"#, i, big_json),
            })
            .collect();
        write_fixture(&base, records);

        let reader = SaReader::open(&base.with_extension("osa")).unwrap();
        // Guard: the fixture must actually contain multiple blocks, otherwise
        // the assertion below is vacuous.
        let total_blocks: usize = reader
            .index
            .chromosomes
            .values()
            .map(|v| v.len())
            .sum();
        assert!(
            total_blocks >= 2,
            "test fixture should have ≥ 2 blocks, got {}",
            total_blocks
        );

        // Preload a single position; the cache should hold exactly the one
        // block that contains it, not the full chromosome.
        reader.preload("chr1", &[1042]).unwrap();
        let cached = reader.block_cache.lock().unwrap().len();
        assert_eq!(
            cached, 1,
            "preload of a single position should load exactly 1 block, got {}",
            cached
        );

        // The preloaded block must satisfy a real query.
        let ann = reader.annotate_position("chr1", 1042, "A", "G").unwrap();
        assert!(ann.is_some());

        // Unknown chromosome must be a no-op rather than an error.
        reader.preload("chrUnknown", &[1, 2, 3]).unwrap();
    }

    #[test]
    fn block_cache_evicts_lru_when_byte_budget_exceeded() {
        // Three "blocks" of 100 bytes each, budget of 250 bytes — the third
        // insert must evict the first to stay within budget.
        let mut cache = BlockCache::new(250);
        let mk = |i: u32| {
            Arc::new(vec![BlockEntry {
                position: i,
                ref_allele: "A".into(),
                alt_allele: "G".into(),
                json: "x".repeat(100),
            }])
        };
        cache.put(0, mk(0), 100);
        cache.put(1, mk(1), 100);
        cache.put(2, mk(2), 100); // evicts offset 0 (LRU)
        assert!(cache.get(0).is_none(), "offset 0 should have been evicted");
        assert!(cache.get(1).is_some());
        assert!(cache.get(2).is_some());
        assert!(cache.total_bytes <= cache.budget_bytes);
    }

    #[test]
    fn block_cache_retains_just_inserted_entry_even_if_oversized() {
        // A single block larger than the entire budget must still be cached;
        // otherwise concurrent workers querying the same oversized block
        // would each re-decompress it.
        let mut cache = BlockCache::new(50);
        let entry = Arc::new(vec![BlockEntry {
            position: 1,
            ref_allele: "A".into(),
            alt_allele: "G".into(),
            json: "x".repeat(1000),
        }]);
        cache.put(0, entry, 1000);
        assert!(cache.get(0).is_some(), "just-inserted block must be retained");
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
