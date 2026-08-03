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
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Arc, Mutex};

/// Soft upper bound on entries to prevent the underlying `LruCache`'s
/// capacity field from being a pathological size if blocks ever ended up
/// being tiny. The byte budget is the real gate.
const CACHE_MAX_ENTRIES: usize = 4096;

/// Monotonic id handed to each `SaReader` so its blocks are namespaced in the
/// shared cache (two shards can share a file offset without colliding). Never
/// reused, so a dropped reader's stale entries can't alias a new reader.
static NEXT_READER_ID: AtomicU64 = AtomicU64::new(0);

/// The one shared, process-wide block cache used by every `.osa` reader.
///
/// **Why one shared cache (issue #75):** the block cache used to be a fixed
/// 32 MiB *per reader*. On a dense, whole-genome source like SpliceAI — where
/// every 8 MiB block spans only ~13 kbp of genome — a batch of coordinate-
/// sorted variants touches far more than 4 blocks, and the `par_iter` annotate
/// phase runs `num_threads` workers that each walk a disjoint sub-range, so the
/// workers continually evict one another's in-flight block and re-decompress
/// the same blocks (the reported "adding spliceai.osa made annotate ~45×
/// slower"). A single shared, byte-budgeted, globally-LRU cache lets whichever
/// readers are hot for the current chromosome share the whole budget while
/// bounding total RAM regardless of how many per-chromosome shards are open —
/// see [`crate::common::sa_cache_budget_bytes`]. The v2 `.osa2` reader has its
/// own equivalent shared chunk cache.
static GLOBAL_BLOCK_CACHE: std::sync::LazyLock<Mutex<BlockCache>> =
    std::sync::LazyLock::new(|| {
        Mutex::new(BlockCache::new(crate::common::sa_cache_budget_bytes()))
    });

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

/// A cache key: which reader a block came from, plus its file offset. The
/// reader id namespaces entries so shards that happen to share an offset don't
/// collide in the single shared cache.
type CacheKey = (u64, u64);

/// LRU keyed by `(reader_id, file_offset)`, with byte-based eviction. Evicts
/// least-recently-used entries until total cached bytes is within budget.
/// Always keeps at least one entry (the just-inserted one) so a single
/// oversized block doesn't keep falling out from under a parallel-worker batch.
///
/// Byte accounting goes through `pop_lru` (explicit) and `push` (which
/// returns the evicted entry on capacity overflow), so the inner `LruCache`
/// can never silently drop an entry without `total_bytes` reflecting it.
struct BlockCache {
    lru: LruCache<CacheKey, (Arc<Vec<BlockEntry>>, usize)>,
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

    fn get(&mut self, key: CacheKey) -> Option<Arc<Vec<BlockEntry>>> {
        self.lru.get(&key).map(|(arc, _)| Arc::clone(arc))
    }

    fn put(&mut self, key: CacheKey, value: Arc<Vec<BlockEntry>>, bytes: usize) {
        // Replace-in-place: drop the old bytes first so the budget loop
        // below sees the correct `total_bytes`.
        if let Some((_, old_bytes)) = self.lru.pop(&key) {
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
        if let Some((_, (_, ev_bytes))) = self.lru.push(key, (value, bytes)) {
            self.total_bytes = self.total_bytes.saturating_sub(ev_bytes);
        }
        self.total_bytes = self.total_bytes.saturating_add(bytes);
    }

    /// Number of cached blocks belonging to a specific reader. Test-only:
    /// the shared cache holds blocks from every open reader, so counting must
    /// filter by id.
    #[cfg(test)]
    fn len_for_reader(&self, reader_id: u64) -> usize {
        self.lru
            .iter()
            .filter(|((rid, _), _)| *rid == reader_id)
            .count()
    }
}

/// Reader for .osa annotation files.
///
/// Thread-safety: the reader is `Send + Sync`. Decompressed blocks live in the
/// process-wide `GLOBAL_BLOCK_CACHE` (a `Mutex<LruCache>`); the lock is held
/// only for the brief lookup or insert and is released across the
/// decompression and find_match steps.
pub struct SaReader {
    mmap: Mmap,
    index: SaIndex,
    metadata: SaMetadata,
    /// Namespaces this reader's blocks in the shared cache.
    reader_id: u64,
    /// Optional private block cache with an explicit byte budget. `None` in
    /// production, where all readers share `GLOBAL_BLOCK_CACHE`; `Some` only
    /// via `open_with_cache_budget`, used by benchmarks and tests to pin a
    /// specific budget in isolation (the shared cache initializes its budget
    /// once, so it can't be varied in-process otherwise).
    local_cache: Option<Mutex<BlockCache>>,
    /// Count of block decompressions performed (cache misses). A diagnostic
    /// for cache thrashing: on a coordinate-sorted input each block should be
    /// decompressed ~once, so a count far above the number of distinct blocks
    /// touched indicates the block cache is too small for the parallel working
    /// set. Exposed via `decompress_count()`.
    decompress_count: AtomicU64,
}

impl SaReader {
    /// Open an .osa + .osa.idx file pair. Blocks are cached in the shared,
    /// byte-budgeted `GLOBAL_BLOCK_CACHE`.
    pub fn open(data_path: &Path) -> Result<Self> {
        Self::open_inner(data_path, None)
    }

    /// Open with a private block cache of exactly `budget_bytes`, bypassing the
    /// shared cache. Intended for benchmarks and tests that need to pin a
    /// specific budget in isolation; production code uses [`open`](Self::open).
    pub fn open_with_cache_budget(data_path: &Path, budget_bytes: usize) -> Result<Self> {
        Self::open_inner(data_path, Some(budget_bytes))
    }

    fn open_inner(data_path: &Path, local_budget: Option<usize>) -> Result<Self> {
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
            reader_id: NEXT_READER_ID.fetch_add(1, Ordering::Relaxed),
            local_cache: local_budget.map(|b| Mutex::new(BlockCache::new(b))),
            decompress_count: AtomicU64::new(0),
        })
    }

    /// Number of block decompressions performed since the reader was opened.
    /// Used by benchmarks/tests to detect cache thrashing.
    pub fn decompress_count(&self) -> u64 {
        self.decompress_count.load(Ordering::Relaxed)
    }

    /// The block cache this reader writes to: its private one if configured,
    /// otherwise the process-wide shared cache.
    fn cache(&self) -> &Mutex<BlockCache> {
        match &self.local_cache {
            Some(c) => c,
            None => &GLOBAL_BLOCK_CACHE,
        }
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
        let len_bytes: [u8; 4] = self.mmap[offset..offset + 4].try_into().map_err(|_| {
            anyhow::anyhow!("expected 4-byte block length prefix at offset {}", offset)
        })?;
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
    /// populating the shared LRU cache as needed.
    fn get_block(&self, block_ref: &BlockRef) -> Result<Arc<Vec<BlockEntry>>> {
        let key: CacheKey = (self.reader_id, block_ref.file_offset);
        let cache_mutex = self.cache();

        // Fast path: cache hit.
        {
            let mut cache = cache_mutex
                .lock()
                .map_err(|_| anyhow::anyhow!("SA block cache mutex poisoned"))?;
            if let Some(arc) = cache.get(key) {
                return Ok(arc);
            }
        }

        // Slow path: decompress without holding the lock so other workers can
        // serve their own queries from the cache concurrently. If two threads
        // race on the same missing block they each decompress once; the second
        // `put` simply replaces an identical entry — acceptable for an LRU.
        let entries = self.decompress_block(block_ref.file_offset, block_ref.compressed_len)?;
        self.decompress_count.fetch_add(1, Ordering::Relaxed);
        let bytes = block_bytes(&entries);
        let arc = Arc::new(entries);

        let mut cache = cache_mutex
            .lock()
            .map_err(|_| anyhow::anyhow!("SA block cache mutex poisoned"))?;
        cache.put(key, Arc::clone(&arc), bytes);
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
        let allele_ref = if self.metadata.match_by_allele {
            ref_allele
        } else {
            ""
        };
        let allele_alt = if self.metadata.match_by_allele {
            alt_allele
        } else {
            ""
        };

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
    use crate::common::{AnnotationRecord, DEFAULT_BLOCK_SIZE, SCHEMA_VERSION};
    use crate::index::IndexHeader;
    use crate::writer::SaWriter;
    use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
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
        let total_blocks: usize = reader.index.chromosomes.values().map(|v| v.len()).sum();
        assert!(
            total_blocks >= 2,
            "test fixture should have ≥ 2 blocks, got {}",
            total_blocks
        );

        // Preload a single position; exactly one block (the one containing it)
        // should be decompressed, not the full multi-block chromosome. The
        // per-reader decompress counter is used rather than the shared cache's
        // length so the assertion is immune to concurrent tests' entries and
        // to LRU eviction.
        reader.preload("chr1", &[1042]).unwrap();
        assert_eq!(
            reader.decompress_count(),
            1,
            "preload of a single position should decompress exactly 1 block"
        );

        // The preloaded block must satisfy a real query.
        let ann = reader.annotate_position("chr1", 1042, "A", "G").unwrap();
        assert!(ann.is_some());

        // Unknown chromosome must be a no-op rather than an error, and must
        // not decompress anything.
        reader.preload("chrUnknown", &[1, 2, 3]).unwrap();
    }

    /// Regression for issue #75: under parallel queries a too-small block
    /// cache makes workers re-decompress the same blocks over and over. With a
    /// cache large enough to hold the working set, each block is decompressed
    /// ~once. This is the whole fix, exercised end-to-end through the real
    /// query path. Uses `open_with_cache_budget` so the two budgets are pinned
    /// deterministically and in isolation from the shared cache.
    #[test]
    fn parallel_queries_do_not_thrash_with_adequate_cache() {
        let dir = TempDir::new().unwrap();
        let base = dir.path().join("dense");
        // ~200k dense records × ~213 B ≈ 42 MiB → ~5 blocks of 8 MiB.
        let pad = "x".repeat(200);
        let n: u32 = 200_000;
        let records: Vec<AnnotationRecord> = (0..n)
            .map(|i| AnnotationRecord {
                chrom_idx: 0,
                position: 1000 + i,
                ref_allele: "A".into(),
                alt_allele: "G".into(),
                json: format!(r#"{{"i":{},"pad":"{}"}}"#, i, pad),
            })
            .collect();
        write_fixture(&base, records);
        let osa = base.with_extension("osa");

        let total_blocks: usize = SaReader::open(&osa)
            .unwrap()
            .index
            .chromosomes
            .values()
            .map(|v| v.len())
            .sum();
        assert!(
            total_blocks > 2,
            "fixture must span several blocks for this test to be meaningful, got {}",
            total_blocks
        );

        // Sorted subset that still spans every block, mimicking a coordinate-
        // sorted VCF. Queried on a fixed 8-thread pool so the parallelism (and
        // thus the thrash) reproduces regardless of the host core count.
        let positions: Vec<u64> = (0..n).step_by(20).map(|i| (1000 + i) as u64).collect();
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(8)
            .build()
            .unwrap();
        // Faithful to `run_annotate`: warm the cache with a single sequential
        // preload, then run the queries in parallel. The preload both primes
        // the cache and pins the baseline — each distinct block is decompressed
        // exactly once here, with no first-touch races.
        let sweep = |reader: &SaReader| -> u64 {
            reader.preload("chr1", &positions).unwrap();
            pool.install(|| {
                positions.par_iter().for_each(|&p| {
                    // All queries must hit; a miss would mean the sweep isn't
                    // exercising the blocks and the counts would be meaningless.
                    assert!(reader
                        .annotate_position("chr1", p, "A", "G")
                        .unwrap()
                        .is_some());
                });
            });
            reader.decompress_count()
        };

        // Cache large enough for every block: the parallel phase adds nothing
        // on top of the preload — each block decompressed exactly once.
        let big = SaReader::open_with_cache_budget(&osa, (total_blocks + 4) * DEFAULT_BLOCK_SIZE)
            .unwrap();
        let big_decomps = sweep(&big);
        assert_eq!(
            big_decomps as usize, total_blocks,
            "adequate cache should decompress each block exactly once"
        );

        // Tiny 2-block cache: the preload's blocks are evicted before the
        // parallel phase, and workers then evict one another's in-flight block,
        // so the same blocks are decompressed many times over.
        let tiny = SaReader::open_with_cache_budget(&osa, 2 * DEFAULT_BLOCK_SIZE).unwrap();
        let tiny_decomps = sweep(&tiny);
        assert!(
            tiny_decomps >= big_decomps * 3,
            "small cache under parallelism must re-decompress heavily: tiny={} big={} blocks={}",
            tiny_decomps,
            big_decomps,
            total_blocks
        );
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
        cache.put((0, 0), mk(0), 100);
        cache.put((0, 1), mk(1), 100);
        cache.put((0, 2), mk(2), 100); // evicts offset 0 (LRU)
        assert!(
            cache.get((0, 0)).is_none(),
            "offset 0 should have been evicted"
        );
        assert!(cache.get((0, 1)).is_some());
        assert!(cache.get((0, 2)).is_some());
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
        cache.put((0, 0), entry, 1000);
        assert!(
            cache.get((0, 0)).is_some(),
            "just-inserted block must be retained"
        );
    }

    #[test]
    fn block_cache_namespaces_by_reader_id() {
        // Two readers with the same file offset must not alias in the shared
        // cache: reader 0's block and reader 1's block at offset 0 coexist.
        let mut cache = BlockCache::new(10_000);
        let mk = |tag: &str| {
            Arc::new(vec![BlockEntry {
                position: 1,
                ref_allele: "A".into(),
                alt_allele: "G".into(),
                json: tag.to_string(),
            }])
        };
        cache.put((0, 0), mk("reader0"), 100);
        cache.put((1, 0), mk("reader1"), 100);
        assert_eq!(cache.get((0, 0)).unwrap()[0].json, "reader0");
        assert_eq!(cache.get((1, 0)).unwrap()[0].json, "reader1");
        assert_eq!(cache.len_for_reader(0), 1);
        assert_eq!(cache.len_for_reader(1), 1);
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
