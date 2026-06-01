//! Reader for .osa2 format (ZIP-based chunked annotation files).
//!
//! Implements the `AnnotationProvider` trait with O(log n) lookups via
//! Var32 binary search on sorted genomic chunks.

use crate::chunk::{delta_decode, Chunk};
use crate::common::chrom_aliases;
use crate::fields::{Field, FieldType};
use crate::kmer16::LongVariant;
use crate::var32;
use crate::writer_v2::{read_u32_array, Osa2Metadata};
use anyhow::{Context, Result};
use lru::LruCache;
use fastvep_cache::annotation::{AnnotationProvider, AnnotationValue, SaMetadata};
use std::collections::HashSet;
use std::fs::File;
use std::io::Read;
use std::num::NonZeroUsize;
use std::path::Path;
use std::sync::{Arc, Mutex};

/// Hard cap on a per-chunk zstd-decompressed JSON blob (256 MiB). Defends
/// against zstd bombs in maliciously crafted .osa2 files.
const MAX_JSON_BLOB_DECOMPRESSED: usize = 256 * 1024 * 1024;

/// Number of recently-used chunks held in the LRU cache.
const CHUNK_CACHE_SIZE: usize = 8;

/// Reader for .osa2 annotation files.
///
/// Loads genomic chunks on demand from a ZIP archive, caches recently used
/// chunks in an LRU cache, and performs binary search for variant lookups.
///
/// The chunk cache is guarded by a `Mutex` because `LruCache::get` is itself
/// a mutating operation (it reorders the LRU recency list), so concurrent
/// queries from multiple worker threads cannot share an `UnsafeCell`.
pub struct Osa2Reader {
    /// Path to the .osa2 ZIP file (re-opened for each chunk load).
    zip_path: std::path::PathBuf,
    metadata: Osa2Metadata,
    sa_metadata: SaMetadata,
    fields: Vec<Field>,
    /// Categorical string lookup tables per field.
    string_tables: Vec<Vec<String>>,
    /// Chromosome names present on-disk under `fastsa/<chrom>/...`.
    /// Used by `resolve_chrom` to canonicalize a query name (e.g. `chr1`
    /// → `1`) before any cache key is built, so a workload that mixes
    /// `chr*` and bare styles for the same physical contig never warms
    /// two cache slots for the same chunk. See issue #37.
    on_disk_chroms: HashSet<String>,
    /// LRU cache of loaded chunks, keyed by `"<canonical-chrom>/<chunk_id>"`.
    chunk_cache: Mutex<LruCache<String, Arc<Chunk>>>,
}

impl Osa2Reader {
    /// Open an .osa2 file.
    pub fn open(path: &Path) -> Result<Self> {
        let file = File::open(path).with_context(|| format!("Opening {}", path.display()))?;
        let mut archive = zip::ZipArchive::new(file)?;

        // Read metadata
        let metadata: Osa2Metadata = {
            let mut entry = archive.by_name("fastsa/metadata.json")
                .context("Missing fastsa/metadata.json")?;
            let mut buf = String::new();
            entry.read_to_string(&mut buf)?;
            serde_json::from_str(&buf)?
        };

        // Read field config
        let fields: Vec<Field> = {
            let mut entry = archive.by_name("fastsa/config.json")
                .context("Missing fastsa/config.json")?;
            let mut buf = String::new();
            entry.read_to_string(&mut buf)?;
            serde_json::from_str(&buf)?
        };

        // Read string tables
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

        // Validate metadata.chunk_bits before it's used as a shift amount
        // and as the within-chunk position width in Var32 keys.
        // A bit-count of 0 would put every variant in chunk 0, and values
        // above var32::CHUNK_BITS would be truncated by Var32 encoding.
        if metadata.chunk_bits == 0 || metadata.chunk_bits > var32::CHUNK_BITS {
            anyhow::bail!(
                "Invalid chunk_bits {} (must be 1..={})",
                metadata.chunk_bits,
                var32::CHUNK_BITS
            );
        }

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

        // Enumerate chromosome directories once so `resolve_chrom` can
        // pick the on-disk canonical name without re-scanning the ZIP per
        // query. Entries look like `fastsa/<chrom>/<chunk_id>/var32.bin`,
        // so we collect the second path segment under `fastsa/`.
        let mut on_disk_chroms = HashSet::new();
        for i in 0..archive.len() {
            let entry = archive.by_index(i)?;
            let name = entry.name();
            let Some(rest) = name.strip_prefix("fastsa/") else {
                continue;
            };
            let Some((chrom, _)) = rest.split_once('/') else {
                continue;
            };
            // Skip the sibling metadata/config/strings subtrees, which
            // share the `fastsa/` parent but aren't chromosome shards.
            if matches!(chrom, "metadata.json" | "config.json" | "strings") {
                continue;
            }
            on_disk_chroms.insert(chrom.to_string());
        }

        let cache_size = NonZeroUsize::new(CHUNK_CACHE_SIZE)
            .expect("CHUNK_CACHE_SIZE is a non-zero compile-time constant");

        Ok(Self {
            zip_path: path.to_path_buf(),
            metadata,
            sa_metadata,
            fields,
            string_tables,
            on_disk_chroms,
            chunk_cache: Mutex::new(LruCache::new(cache_size)),
        })
    }

    /// Resolve a query chromosome name (e.g. `chr1`, `1`, `chrM`, `MT`)
    /// to the canonical on-disk name actually present in the archive.
    /// Returns the input unchanged when no alias is present — callers
    /// then naturally produce empty results.
    fn resolve_chrom<'a>(&self, chrom: &'a str) -> String {
        if self.on_disk_chroms.contains(chrom) {
            return chrom.to_string();
        }
        chrom_aliases(chrom)
            .into_iter()
            .find(|alias| self.on_disk_chroms.contains(alias))
            .unwrap_or_else(|| chrom.to_string())
    }

    /// Build a chunk by reading its files from the ZIP archive. Pure (no cache
    /// access), so it can run with the cache mutex unheld and avoid blocking
    /// other readers during disk I/O.
    ///
    /// Each sub-entry is treated as follows:
    ///   * absent from the archive (`ZipError::FileNotFound`) — legitimate
    ///     "this chunk has no data for this field", continue with the empty
    ///     case;
    ///   * any other archive error — propagate, since this means the .osa2 is
    ///     corrupt or unreadable.
    /// Earlier revisions collapsed both cases together via `Err(_) => …`,
    /// which silently turned data corruption into false-negative lookups.
    fn build_chunk(&self, chrom: &str, chunk_id: u32) -> Result<Chunk> {
        let file = File::open(&self.zip_path)?;
        let mut archive = zip::ZipArchive::new(file)?;
        // `chrom` is expected to be canonical (resolved by `resolve_chrom`
        // at the public entry points), so no alias walk is needed here.
        let prefix = format!("fastsa/{}/{}/", chrom, chunk_id);

        // Read var32 keys. If absent, this chunk has no short variants at all,
        // which also implies no long variants and no parallel value arrays —
        // short-circuit to an empty chunk.
        let var32s = match archive.by_name(&format!("{}var32.bin", prefix)) {
            Ok(mut entry) => {
                let mut buf = Vec::new();
                entry.read_to_end(&mut buf)?;
                let mut keys = read_u32_array(&buf)?;
                delta_decode(&mut keys);
                keys
            }
            Err(zip::result::ZipError::FileNotFound) => return Ok(Chunk::empty()),
            Err(e) => {
                return Err(anyhow::anyhow!(
                    "failed to read var32 entry for chunk {}/{}: {}",
                    chrom,
                    chunk_id,
                    e
                ));
            }
        };

        // Read long variants.
        let longs: Vec<LongVariant> = match archive.by_name(&format!("{}too-long.enc", prefix)) {
            Ok(mut entry) => {
                let mut buf = Vec::new();
                entry.read_to_end(&mut buf)?;
                // Propagate bincode errors so a corrupt `too-long.enc` is
                // reported instead of silently masquerading as "no long
                // variants in this chunk".
                bincode::deserialize(&buf).map_err(|e| {
                    anyhow::anyhow!(
                        "failed to deserialize long-variant block for chunk {}/{}: {}",
                        chrom,
                        chunk_id,
                        e
                    )
                })?
            }
            Err(zip::result::ZipError::FileNotFound) => Vec::new(),
            Err(e) => {
                return Err(anyhow::anyhow!(
                    "failed to read long-variant entry for chunk {}/{}: {}",
                    chrom,
                    chunk_id,
                    e
                ));
            }
        };

        // Read parallel value arrays (one per non-JsonBlob field, in field order).
        let mut values = Vec::new();
        for field in &self.fields {
            if field.ftype == FieldType::JsonBlob {
                continue;
            }
            let name = format!("{}{}.bin", prefix, field.alias);
            match archive.by_name(&name) {
                Ok(mut entry) => {
                    let mut buf = Vec::new();
                    entry.read_to_end(&mut buf)?;
                    values.push(read_u32_array(&buf)?);
                }
                Err(zip::result::ZipError::FileNotFound) => {
                    values.push(vec![field.missing_value; var32s.len()]);
                }
                Err(e) => {
                    return Err(anyhow::anyhow!(
                        "failed to read value column '{}' for chunk {}/{}: {}",
                        field.alias,
                        chrom,
                        chunk_id,
                        e
                    ));
                }
            }
        }

        // Read JSON blobs if any
        let json_blobs = match archive.by_name(&format!("{}json_blobs.zst", prefix)) {
            Ok(mut entry) => {
                let mut buf = Vec::new();
                entry.read_to_end(&mut buf)?;
                // Bound the decompressed size to defend against zstd bombs.
                let mut decoder = zstd::stream::Decoder::new(buf.as_slice())?;
                let mut decompressed = Vec::new();
                let mut limited = (&mut decoder).take(MAX_JSON_BLOB_DECOMPRESSED as u64 + 1);
                limited.read_to_end(&mut decompressed)?;
                if decompressed.len() > MAX_JSON_BLOB_DECOMPRESSED {
                    anyhow::bail!(
                        "JSON blob decompressed size exceeds limit ({} bytes)",
                        MAX_JSON_BLOB_DECOMPRESSED
                    );
                }
                let text = String::from_utf8(decompressed)?;
                Some(text.lines().map(|l| l.to_string()).collect())
            }
            Err(zip::result::ZipError::FileNotFound) => None,
            Err(e) => {
                return Err(anyhow::anyhow!(
                    "failed to read json_blobs for chunk {}/{}: {}",
                    chrom,
                    chunk_id,
                    e
                ));
            }
        };

        Ok(Chunk { var32s, longs, values, json_blobs })
    }

    /// Ensure a chunk is in the LRU cache. Idempotent.
    fn load_chunk(&self, chrom: &str, chunk_id: u32) -> Result<()> {
        let cache_key = format!("{}/{}", chrom, chunk_id);
        {
            let cache = self
                .chunk_cache
                .lock()
                .map_err(|_| anyhow::anyhow!("chunk_cache mutex poisoned"))?;
            if cache.contains(&cache_key) {
                return Ok(());
            }
        }

        // Build the chunk without holding the lock to avoid blocking readers
        // during disk I/O. Two concurrent loads of the same chunk will each
        // build a chunk; the LRU keeps only the last `put`, which is acceptable.
        let chunk = Arc::new(self.build_chunk(chrom, chunk_id)?);

        let mut cache = self
            .chunk_cache
            .lock()
            .map_err(|_| anyhow::anyhow!("chunk_cache mutex poisoned"))?;
        cache.put(cache_key, chunk);
        Ok(())
    }

    /// Query a variant in the loaded chunks.
    fn query(&self, chrom: &str, pos: u32, ref_allele: &[u8], alt_allele: &[u8]) -> Result<Option<String>> {
        // Canonicalize before constructing the cache key so `chr1` and `1`
        // (same physical chunk) share a single LRU slot.
        let chrom = self.resolve_chrom(chrom);
        let chunk_id = pos >> self.metadata.chunk_bits;
        let cache_key = format!("{}/{}", chrom, chunk_id);

        // Ensure chunk is loaded
        self.load_chunk(&chrom, chunk_id)?;

        // `LruCache::get` mutates recency order, so lock only for lookup and
        // clone an `Arc` to release the mutex before search/reconstruction.
        let chunk = {
            let mut cache = self
                .chunk_cache
                .lock()
                .map_err(|_| anyhow::anyhow!("chunk_cache mutex poisoned"))?;
            match cache.get(&cache_key) {
                Some(c) => Arc::clone(c),
                None => return Ok(None),
            }
        };

        if chunk.is_empty() {
            return Ok(None);
        }

        // chunk_bits validated in `open()` so the shift below is well-defined.
        let chunk_mask = (1u32 << self.metadata.chunk_bits) - 1;
        let within_pos = pos & chunk_mask;

        let idx = if var32::is_long(ref_allele.len(), alt_allele.len()) {
            chunk.find_long(pos, ref_allele, alt_allele)
        } else {
            var32::encode(within_pos, ref_allele, alt_allele)
                .and_then(|key| chunk.find_short(key))
        };

        match idx {
            Some(i) => {
                // Defensive bounds check: `find_long` returns an index baked
                // into the on-disk `LongVariant` record, so an internally
                // inconsistent or corrupted chunk could yield an index past
                // the value columns and JSON-blob array. Without this guard
                // `reconstruct_json` would silently return `{}` and the caller
                // would treat it as a positive match.
                //
                // We take the *max* across the column lengths and the
                // json_blobs length (the writer keeps them parallel to the
                // sorted record order). A truly corrupt chunk with a value
                // column shorter than the others is caught by the per-column
                // `column.get(idx)` guard inside `reconstruct_json`.
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
                let json = chunk.reconstruct_json(i, &self.fields, &self.string_tables);
                Ok(Some(json))
            }
            None => Ok(None),
        }
    }
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
        let ref_bytes = ref_allele.as_bytes();
        let alt_bytes = alt_allele.as_bytes();

        let pos_u32: u32 = pos
            .try_into()
            .map_err(|_| anyhow::anyhow!("Position {} exceeds u32::MAX", pos))?;
        match self.query(chrom, pos_u32, ref_bytes, alt_bytes)? {
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

        // Canonicalize once so the chunks preloaded here share cache keys
        // with the subsequent `annotate_position` calls that follow,
        // regardless of which naming style each side uses.
        let chrom = self.resolve_chrom(chrom);

        // Determine which chunks need to be loaded. Reject positions that
        // overflow u32 rather than silently truncating into the wrong chunk.
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
            self.load_chunk(&chrom, cid)?;
        }

        Ok(())
    }
}
