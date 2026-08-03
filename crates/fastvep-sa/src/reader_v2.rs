//! Reader for .osa2 format (ZIP-based chunked annotation files).
//!
//! Implements the `AnnotationProvider` trait with O(log n) lookups via
//! Var32 binary search on sorted genomic chunks.
//!
//! **I/O model (issue #75):** the ZIP central directory is parsed once at
//! `open` and every entry's compressed byte range is recorded up front, so a
//! chunk read is a lock-free `mmap` slice plus an in-thread inflate — no shared
//! `ZipArchive` mutex on the query path. Decompressed chunks live in a single
//! process-wide, byte-budgeted, globally-LRU cache shared by all `.osa2`
//! readers (mirror of the `.osa` block cache), so a dense whole-genome source
//! queried in parallel neither serializes on a lock nor thrashes a too-small
//! per-reader cache. See [`crate::common::sa_cache_budget_bytes`].

use crate::chunk::{delta_decode, Chunk};
use crate::common::chrom_aliases;
use crate::fields::{Field, FieldType};
use crate::kmer16::LongVariant;
use crate::var32;
use crate::writer_v2::{read_u32_array, Osa2Metadata};
use anyhow::{Context, Result};
use fastvep_cache::annotation::{AnnotationProvider, AnnotationValue, SaMetadata};
use flate2::read::DeflateDecoder;
use lru::LruCache;
use memmap2::Mmap;
use std::collections::HashMap;
use std::fs::File;
use std::io::Read;
use std::num::NonZeroUsize;
use std::path::Path;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Arc, Mutex};

/// Hard cap on a per-chunk zstd-decompressed JSON blob (256 MiB). Defends
/// against zstd bombs in maliciously crafted .osa2 files.
const MAX_JSON_BLOB_DECOMPRESSED: usize = 256 * 1024 * 1024;

/// Soft cap on cached chunk *entries*; the byte budget is the real gate.
const CHUNK_CACHE_MAX_ENTRIES: usize = 4096;

/// Monotonic id per reader so chunks are namespaced in the shared cache
/// (two shards can share a `chunk_id` without colliding). Never reused.
static NEXT_READER_ID: AtomicU64 = AtomicU64::new(0);

/// The one shared, process-wide chunk cache used by every `.osa2` reader. See
/// the module docs and [`crate::common::sa_cache_budget_bytes`] for why it is
/// shared and byte-budgeted rather than a fixed count per reader.
static GLOBAL_CHUNK_CACHE: std::sync::LazyLock<Mutex<ChunkCache>> =
    std::sync::LazyLock::new(|| {
        Mutex::new(ChunkCache::new(crate::common::sa_cache_budget_bytes()))
    });

/// A cache key: which reader a chunk came from, plus a chromosome-namespaced
/// chunk id (`(chrom_idx << 32) | chunk_id`).
type ChunkKey = (u64, u64);

/// Approximate resident footprint of a decompressed chunk: the parallel u32
/// arrays plus the heap behind any JSON-blob strings.
fn chunk_bytes(c: &Chunk) -> usize {
    let v32 = c.var32s.len() * std::mem::size_of::<u32>();
    let longs = c.longs.len() * std::mem::size_of::<LongVariant>();
    let vals: usize = c
        .values
        .iter()
        .map(|col| col.len() * std::mem::size_of::<u32>())
        .sum();
    let blobs: usize = c.json_blobs.as_ref().map_or(0, |b| {
        b.iter()
            .map(|s| s.len() + std::mem::size_of::<String>())
            .sum()
    });
    v32.saturating_add(longs)
        .saturating_add(vals)
        .saturating_add(blobs)
}

/// Byte-budgeted LRU of decompressed chunks. Same accounting discipline as the
/// `.osa` block cache: eviction goes through `pop_lru`/`push` so `total_bytes`
/// can never drift, and the just-inserted entry is always retained even if it
/// alone exceeds the budget.
struct ChunkCache {
    lru: LruCache<ChunkKey, (Arc<Chunk>, usize)>,
    total_bytes: usize,
    budget_bytes: usize,
}

impl ChunkCache {
    fn new(budget_bytes: usize) -> Self {
        let cap = NonZeroUsize::new(CHUNK_CACHE_MAX_ENTRIES).expect("non-zero");
        Self {
            lru: LruCache::new(cap),
            total_bytes: 0,
            budget_bytes,
        }
    }

    fn get(&mut self, key: ChunkKey) -> Option<Arc<Chunk>> {
        self.lru.get(&key).map(|(arc, _)| Arc::clone(arc))
    }

    fn put(&mut self, key: ChunkKey, value: Arc<Chunk>, bytes: usize) {
        if let Some((_, old_bytes)) = self.lru.pop(&key) {
            self.total_bytes = self.total_bytes.saturating_sub(old_bytes);
        }
        while self.total_bytes + bytes > self.budget_bytes && !self.lru.is_empty() {
            if let Some((_, (_, ev_bytes))) = self.lru.pop_lru() {
                self.total_bytes = self.total_bytes.saturating_sub(ev_bytes);
            } else {
                break;
            }
        }
        if let Some((_, (_, ev_bytes))) = self.lru.push(key, (value, bytes)) {
            self.total_bytes = self.total_bytes.saturating_sub(ev_bytes);
        }
        self.total_bytes = self.total_bytes.saturating_add(bytes);
    }
}

/// ZIP compression method for a single entry (only the two the writer emits).
#[derive(Clone, Copy, PartialEq, Eq)]
enum EntryMethod {
    Stored,
    Deflated,
}

/// Location of one ZIP entry's compressed data within the mmap, resolved once
/// at `open` so query-time reads never touch the `ZipArchive` (and thus never
/// take a shared lock).
#[derive(Clone, Copy)]
struct EntryLoc {
    data_start: u64,
    comp_size: u64,
    method: EntryMethod,
}

/// Reader for .osa2 annotation files.
pub struct Osa2Reader {
    /// The whole file, memory-mapped. Chunk reads slice directly into this.
    mmap: Mmap,
    /// Entry name → compressed byte range + method, from the central directory.
    entries: HashMap<String, EntryLoc>,
    /// Namespaces this reader's chunks in the shared cache.
    reader_id: u64,
    /// Optional private chunk cache with an explicit budget. `None` in
    /// production (all readers share `GLOBAL_CHUNK_CACHE`); `Some` only via
    /// `open_with_cache_budget`, for benchmarks/tests that pin a budget.
    local_cache: Option<Mutex<ChunkCache>>,
    /// Count of chunks built (cache misses) — a thrash diagnostic, exposed via
    /// `chunk_load_count()`.
    chunk_load_count: AtomicU64,
    metadata: Osa2Metadata,
    sa_metadata: SaMetadata,
    fields: Vec<Field>,
    /// Categorical string lookup tables per field.
    string_tables: Vec<Vec<String>>,
    /// On-disk chromosome name → small dense index. Serves two roles:
    /// `resolve_chrom` uses the keys to canonicalize a query name before any
    /// cache key is built (issue #37), and the index is folded into the chunk
    /// cache key so that `chr1`'s chunk 0 and `chr2`'s chunk 0 (same numeric
    /// `chunk_id`) never collide in the shared cache.
    chrom_index: HashMap<String, u32>,
}

impl Osa2Reader {
    /// Open an .osa2 file. Chunks are cached in the shared `GLOBAL_CHUNK_CACHE`.
    pub fn open(path: &Path) -> Result<Self> {
        Self::open_inner(path, None)
    }

    /// Open with a private chunk cache of exactly `budget_bytes`, bypassing the
    /// shared cache. For benchmarks/tests that need to pin a budget in
    /// isolation; production uses [`open`](Self::open).
    pub fn open_with_cache_budget(path: &Path, budget_bytes: usize) -> Result<Self> {
        Self::open_inner(path, Some(budget_bytes))
    }

    fn open_inner(path: &Path, local_budget: Option<usize>) -> Result<Self> {
        let file = File::open(path).with_context(|| format!("Opening {}", path.display()))?;
        // SAFETY: the file is opened read-only and the mmap is only read from.
        let mmap = unsafe { Mmap::map(&file)? };

        // Parse the central directory once to read the small metadata entries
        // and to record every entry's compressed byte range for the lock-free
        // read path. The archive is dropped at the end of this function.
        let mut archive = zip::ZipArchive::new(file)?;

        let metadata: Osa2Metadata = {
            let mut entry = archive
                .by_name("fastsa/metadata.json")
                .context("Missing fastsa/metadata.json")?;
            let mut buf = String::new();
            entry.read_to_string(&mut buf)?;
            serde_json::from_str(&buf)?
        };

        let fields: Vec<Field> = {
            let mut entry = archive
                .by_name("fastsa/config.json")
                .context("Missing fastsa/config.json")?;
            let mut buf = String::new();
            entry.read_to_string(&mut buf)?;
            serde_json::from_str(&buf)?
        };

        let mut string_tables: Vec<Vec<String>> = fields.iter().map(|_| Vec::new()).collect();
        for (i, field) in fields.iter().enumerate() {
            if field.ftype == FieldType::Categorical {
                let name = format!("fastsa/strings/{}.txt", field.alias);
                if let Ok(mut entry) = archive.by_name(&name) {
                    let mut buf = String::new();
                    entry.read_to_string(&mut buf)?;
                    string_tables[i] = buf.lines().map(|l| l.to_string()).collect();
                }
            }
        }

        // chunk_bits is used as a shift amount and as the within-chunk position
        // width in Var32 keys: 0 collapses every variant into chunk 0, and
        // values above var32::CHUNK_BITS would be truncated.
        if metadata.chunk_bits == 0 || metadata.chunk_bits > var32::CHUNK_BITS {
            anyhow::bail!(
                "Invalid chunk_bits {} (must be 1..={})",
                metadata.chunk_bits,
                var32::CHUNK_BITS
            );
        }

        // Single pass over the central directory: record each entry's data
        // range and assign each chromosome directory a dense index.
        let mut entries = HashMap::with_capacity(archive.len());
        let mut chrom_index: HashMap<String, u32> = HashMap::new();
        for i in 0..archive.len() {
            let zf = archive.by_index(i)?;
            let name = zf.name().to_string();
            let method = match zf.compression() {
                zip::CompressionMethod::Stored => EntryMethod::Stored,
                zip::CompressionMethod::Deflated => EntryMethod::Deflated,
                other => anyhow::bail!(
                    "Unsupported ZIP compression {:?} for entry '{}' in .osa2",
                    other,
                    name
                ),
            };
            let data_start = local_data_start(&mmap, zf.header_start())?;
            let comp_size = zf.compressed_size();
            entries.insert(
                name.clone(),
                EntryLoc {
                    data_start,
                    comp_size,
                    method,
                },
            );

            if let Some(rest) = name.strip_prefix("fastsa/") {
                if let Some((chrom, _)) = rest.split_once('/') {
                    if !matches!(chrom, "metadata.json" | "config.json" | "strings") {
                        let next = chrom_index.len() as u32;
                        chrom_index.entry(chrom.to_string()).or_insert(next);
                    }
                }
            }
        }
        drop(archive);

        let sa_metadata = SaMetadata {
            name: metadata.name.clone(),
            version: metadata.version.clone(),
            description: metadata.description.clone(),
            assembly: metadata.assembly.clone(),
            json_key: metadata.json_key.clone(),
            match_by_allele: metadata.match_by_allele,
            is_array: metadata.is_array,
            is_positional: metadata.is_positional,
        };

        Ok(Self {
            mmap,
            entries,
            reader_id: NEXT_READER_ID.fetch_add(1, Ordering::Relaxed),
            local_cache: local_budget.map(|b| Mutex::new(ChunkCache::new(b))),
            chunk_load_count: AtomicU64::new(0),
            metadata,
            sa_metadata,
            fields,
            string_tables,
            chrom_index,
        })
    }

    /// Number of chunks built (cache misses) since the reader was opened. Used
    /// by benchmarks/tests to detect chunk-cache thrashing.
    pub fn chunk_load_count(&self) -> u64 {
        self.chunk_load_count.load(Ordering::Relaxed)
    }

    /// The chunk cache this reader writes to: its private one if configured,
    /// otherwise the process-wide shared cache.
    fn cache(&self) -> &Mutex<ChunkCache> {
        match &self.local_cache {
            Some(c) => c,
            None => &GLOBAL_CHUNK_CACHE,
        }
    }

    /// Resolve a query chromosome name (`chr1`, `1`, `chrM`, `MT`, …) to the
    /// canonical on-disk name present in the archive. Returns the input
    /// unchanged when no alias matches — callers then naturally miss.
    fn resolve_chrom(&self, chrom: &str) -> String {
        if self.chrom_index.contains_key(chrom) {
            return chrom.to_string();
        }
        chrom_aliases(chrom)
            .into_iter()
            .find(|alias| self.chrom_index.contains_key(alias))
            .unwrap_or_else(|| chrom.to_string())
    }

    /// Read and decompress one ZIP entry straight from the mmap, with no lock.
    /// `Ok(None)` means the entry is absent (a legitimate "this chunk has no
    /// data for this field"); `Err` means the archive is corrupt/unreadable.
    fn read_entry(&self, name: &str) -> Result<Option<Vec<u8>>> {
        let Some(loc) = self.entries.get(name) else {
            return Ok(None);
        };
        let start = loc.data_start as usize;
        let end = start
            .checked_add(loc.comp_size as usize)
            .ok_or_else(|| anyhow::anyhow!("entry '{}' data range overflow", name))?;
        if end > self.mmap.len() {
            anyhow::bail!("entry '{}' extends beyond .osa2 file", name);
        }
        let raw = &self.mmap[start..end];
        match loc.method {
            EntryMethod::Stored => Ok(Some(raw.to_vec())),
            EntryMethod::Deflated => {
                let mut out = Vec::new();
                DeflateDecoder::new(raw)
                    .read_to_end(&mut out)
                    .with_context(|| format!("inflating ZIP entry '{}'", name))?;
                Ok(Some(out))
            }
        }
    }

    /// Build a chunk by reading its files from the mmap. Pure (no cache access)
    /// and lock-free, so many workers can build different chunks concurrently.
    ///
    /// Each sub-entry: absent → the empty/default case; any read/inflate error
    /// → propagated, since it means the .osa2 is corrupt (never silently turned
    /// into a false-negative lookup).
    fn build_chunk(&self, chrom: &str, chunk_id: u32) -> Result<Chunk> {
        // `chrom` is expected canonical (resolved at the public entry points).
        let prefix = format!("fastsa/{}/{}/", chrom, chunk_id);

        // Var32 keys. Absent ⇒ this chunk has no short variants, which implies
        // no long variants and no value arrays — short-circuit to empty.
        let var32s = match self.read_entry(&format!("{}var32.bin", prefix))? {
            Some(buf) => {
                let mut keys = read_u32_array(&buf)?;
                delta_decode(&mut keys);
                keys
            }
            None => return Ok(Chunk::empty()),
        };

        let longs: Vec<LongVariant> = match self.read_entry(&format!("{}too-long.enc", prefix))? {
            Some(buf) => bincode::deserialize(&buf).map_err(|e| {
                anyhow::anyhow!(
                    "failed to deserialize long-variant block for chunk {}/{}: {}",
                    chrom,
                    chunk_id,
                    e
                )
            })?,
            None => Vec::new(),
        };

        // Parallel value arrays (one per non-JsonBlob field, in field order).
        let mut values = Vec::new();
        for field in &self.fields {
            if field.ftype == FieldType::JsonBlob {
                continue;
            }
            match self.read_entry(&format!("{}{}.bin", prefix, field.alias))? {
                Some(buf) => values.push(read_u32_array(&buf)?),
                None => values.push(vec![field.missing_value; var32s.len()]),
            }
        }

        // JSON blobs, if any (content is zstd-compressed inside the ZIP entry).
        let json_blobs = match self.read_entry(&format!("{}json_blobs.zst", prefix))? {
            Some(buf) => {
                let mut decoder = zstd::stream::Decoder::new(buf.as_slice())?;
                let mut decompressed = Vec::new();
                (&mut decoder)
                    .take(MAX_JSON_BLOB_DECOMPRESSED as u64 + 1)
                    .read_to_end(&mut decompressed)?;
                if decompressed.len() > MAX_JSON_BLOB_DECOMPRESSED {
                    anyhow::bail!(
                        "JSON blob decompressed size exceeds limit ({} bytes)",
                        MAX_JSON_BLOB_DECOMPRESSED
                    );
                }
                let text = String::from_utf8(decompressed)?;
                Some(text.lines().map(|l| l.to_string()).collect())
            }
            None => None,
        };

        Ok(Chunk {
            var32s,
            longs,
            values,
            json_blobs,
        })
    }

    /// Return the chunk, hitting or populating the (shared) cache. Two workers
    /// racing on the same missing chunk each build once; the second `put`
    /// replaces an identical entry — acceptable for an LRU.
    fn get_chunk(&self, chrom: &str, chunk_id: u32) -> Result<Arc<Chunk>> {
        // Namespace the key by chromosome: `chr1`'s chunk 0 and `chr2`'s chunk
        // 0 share a numeric `chunk_id` but are different chunks. An unknown
        // chromosome (not on disk) maps to a shared sentinel — its chunks are
        // always empty, so collisions there are harmless.
        let chrom_idx = self.chrom_index.get(chrom).copied().unwrap_or(u32::MAX);
        let namespaced = ((chrom_idx as u64) << 32) | (chunk_id as u64);
        let key: ChunkKey = (self.reader_id, namespaced);
        let cache_mutex = self.cache();

        {
            let mut cache = cache_mutex
                .lock()
                .map_err(|_| anyhow::anyhow!("chunk cache mutex poisoned"))?;
            if let Some(arc) = cache.get(key) {
                return Ok(arc);
            }
        }

        // Build without holding the lock (lock-free mmap reads + inflate).
        let chunk = Arc::new(self.build_chunk(chrom, chunk_id)?);
        self.chunk_load_count.fetch_add(1, Ordering::Relaxed);
        let bytes = chunk_bytes(&chunk);

        let mut cache = cache_mutex
            .lock()
            .map_err(|_| anyhow::anyhow!("chunk cache mutex poisoned"))?;
        cache.put(key, Arc::clone(&chunk), bytes);
        Ok(chunk)
    }

    /// Query a variant in the chunk that would contain it.
    fn query(
        &self,
        chrom: &str,
        pos: u32,
        ref_allele: &[u8],
        alt_allele: &[u8],
    ) -> Result<Option<String>> {
        let chrom = self.resolve_chrom(chrom);
        let chunk_id = pos >> self.metadata.chunk_bits;
        let chunk = self.get_chunk(&chrom, chunk_id)?;

        if chunk.is_empty() {
            return Ok(None);
        }

        // chunk_bits validated in `open()` so the shift below is well-defined.
        let chunk_mask = (1u32 << self.metadata.chunk_bits) - 1;
        let within_pos = pos & chunk_mask;

        let idx = if self.metadata.is_positional {
            // Positional sources match by coordinate alone.
            chunk.find_short(var32::positional_key(within_pos))
        } else if var32::is_long(ref_allele.len(), alt_allele.len()) {
            chunk.find_long(pos, ref_allele, alt_allele)
        } else {
            var32::encode(within_pos, ref_allele, alt_allele).and_then(|key| chunk.find_short(key))
        };

        match idx {
            Some(i) => {
                // Defensive bounds check: `find_long` returns an index baked
                // into the on-disk record, so a corrupt chunk could yield an
                // index past the value columns and JSON-blob array. Without
                // this guard `reconstruct_json` would silently return `{}` and
                // the caller would treat it as a positive match. Take the max
                // across column lengths and the json_blobs length (all kept
                // parallel to the sorted record order by the writer).
                let data_len = chunk
                    .values
                    .iter()
                    .map(|c| c.len())
                    .max()
                    .unwrap_or(0)
                    .max(chunk.json_blobs.as_ref().map_or(0, |b| b.len()));
                if i >= data_len {
                    return Ok(None);
                }
                Ok(Some(chunk.reconstruct_json(
                    i,
                    &self.fields,
                    &self.string_tables,
                )))
            }
            None => Ok(None),
        }
    }
}

/// Compute the start offset of a ZIP entry's data from its local file header,
/// parsed directly out of the mmap. Robust against the crate's lazily-set
/// `data_start` (which panics if the entry hasn't been read) and avoids reading
/// any entry data at open. Local file header layout: 4-byte signature,
/// 22 bytes of fixed fields, then u16 name length and u16 extra length at
/// offsets 26/28, followed by name + extra, then the data.
fn local_data_start(mmap: &Mmap, header_start: u64) -> Result<u64> {
    let h = header_start as usize;
    let fixed_end = h
        .checked_add(30)
        .ok_or_else(|| anyhow::anyhow!("ZIP local header offset overflow"))?;
    if fixed_end > mmap.len() {
        anyhow::bail!("ZIP local header at {} extends beyond file", header_start);
    }
    if &mmap[h..h + 4] != b"PK\x03\x04" {
        anyhow::bail!("bad ZIP local header signature at {}", header_start);
    }
    let name_len = u16::from_le_bytes([mmap[h + 26], mmap[h + 27]]) as u64;
    let extra_len = u16::from_le_bytes([mmap[h + 28], mmap[h + 29]]) as u64;
    Ok(header_start + 30 + name_len + extra_len)
}

impl AnnotationProvider for Osa2Reader {
    fn name(&self) -> &str {
        &self.sa_metadata.name
    }

    fn json_key(&self) -> &str {
        &self.sa_metadata.json_key
    }

    fn metadata(&self) -> &SaMetadata {
        &self.sa_metadata
    }

    fn annotate_position(
        &self,
        chrom: &str,
        pos: u64,
        ref_allele: &str,
        alt_allele: &str,
    ) -> Result<Option<AnnotationValue>> {
        let pos_u32: u32 = pos
            .try_into()
            .map_err(|_| anyhow::anyhow!("Position {} exceeds u32::MAX", pos))?;
        match self.query(chrom, pos_u32, ref_allele.as_bytes(), alt_allele.as_bytes())? {
            Some(json) => {
                if self.sa_metadata.is_positional {
                    Ok(Some(AnnotationValue::Positional(json)))
                } else {
                    Ok(Some(AnnotationValue::Json(json)))
                }
            }
            None => Ok(None),
        }
    }

    fn preload(&self, chrom: &str, positions: &[u64]) -> Result<()> {
        if positions.is_empty() {
            return Ok(());
        }
        // Canonicalize once so preloaded chunks share cache keys with the
        // `annotate_position` calls that follow.
        let chrom = self.resolve_chrom(chrom);

        let mut chunk_ids: Vec<u32> = Vec::with_capacity(positions.len());
        for &p in positions {
            let p32: u32 = p
                .try_into()
                .map_err(|_| anyhow::anyhow!("Position {} exceeds u32::MAX", p))?;
            chunk_ids.push(p32 >> self.metadata.chunk_bits);
        }
        chunk_ids.sort_unstable();
        chunk_ids.dedup();

        for cid in chunk_ids {
            self.get_chunk(&chrom, cid)?;
        }
        Ok(())
    }
}
