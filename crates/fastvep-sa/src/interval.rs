//! Interval-based annotation reader/writer (.osi files).
//!
//! Used for structural variant databases (gnomAD SV, ClinGen dosage, DGV)
//! where annotations are regions rather than point positions.

use crate::common::{IntervalRecord, MAX_INDEX_PAYLOAD, OSI_MAGIC, SCHEMA_VERSION};
use anyhow::Result;
use fastvep_cache::annotation::{AnnotationProvider, AnnotationValue, SaMetadata};
use fastvep_core::chrom_alias_map;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::io::{Read, Write};
use std::path::Path;

/// Header for .osi files.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct IntervalHeader {
    pub schema_version: u16,
    pub json_key: String,
    pub name: String,
    pub version: String,
    pub assembly: String,
}

/// In-memory interval database loaded from an .osi file.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct IntervalIndex {
    pub header: IntervalHeader,
    /// Chromosome -> sorted list of intervals.
    pub intervals: HashMap<String, Vec<StoredInterval>>,
    /// Per-chromosome prefix maxima of `end`: `end_maxima[chrom][i] ==
    /// max(intervals[chrom][0..=i].end)`. Because it is monotonically
    /// non-decreasing it can be binary-searched to find the first interval
    /// that could possibly reach a query, which is what turns
    /// [`find_overlapping`](Self::find_overlapping) from a scan-from-index-0
    /// into two binary searches plus a bounded scan.
    ///
    /// Derived state, not part of the on-disk format: `#[serde(skip)]` keeps
    /// the `.osi` bytes byte-identical to files written before this field
    /// existed, and it is rebuilt by [`sort`](Self::sort) and
    /// [`read_from`](Self::read_from).
    #[serde(skip)]
    end_maxima: HashMap<String, Vec<u32>>,
    /// Every accepted spelling of every chromosome in the index → its key in
    /// `intervals` (issue #37). Same derived-state contract as `end_maxima`:
    /// `#[serde(skip)]`, rebuilt by [`sort`](Self::sort) and
    /// [`read_from`](Self::read_from).
    ///
    /// Replaces a per-query `chrom_aliases` call - which allocates a
    /// `Vec<String>` and probed the map once per alias - with a single
    /// allocation-free lookup.
    #[serde(skip)]
    chrom_lookup: HashMap<String, String>,
}

/// A stored interval with its annotation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StoredInterval {
    pub start: u32,
    pub end: u32,
    pub json: String,
}

impl IntervalIndex {
    /// Create a new empty interval index.
    pub fn new(header: IntervalHeader) -> Self {
        Self {
            header,
            intervals: HashMap::new(),
            end_maxima: HashMap::new(),
            chrom_lookup: HashMap::new(),
        }
    }

    /// Add an interval record.
    ///
    /// Invalidates the derived `end_maxima` for that chromosome; call
    /// [`sort`](Self::sort) once all records are added to rebuild it (and to
    /// establish the start-order invariant `find_overlapping` needs). The alias
    /// map is extended eagerly instead, so a query works even before `sort`.
    pub fn add(&mut self, record: IntervalRecord) {
        self.end_maxima.remove(&record.chrom);
        if !self.intervals.contains_key(&record.chrom) {
            for (alias, canonical) in chrom_alias_map([&record.chrom]) {
                self.chrom_lookup.entry(alias).or_insert(canonical);
            }
        }
        self.intervals
            .entry(record.chrom)
            .or_default()
            .push(StoredInterval {
                start: record.start,
                end: record.end,
                json: record.json,
            });
    }

    /// Sort all intervals by start position (call after adding all records).
    ///
    /// Also rebuilds the derived `end_maxima` and `chrom_lookup` structures
    /// [`find_overlapping`](Self::find_overlapping) relies on.
    pub fn sort(&mut self) {
        for intervals in self.intervals.values_mut() {
            intervals.sort_by_key(|i| i.start);
        }
        self.rebuild_derived();
    }

    /// Recompute every derived lookup structure from `intervals`.
    fn rebuild_derived(&mut self) {
        self.rebuild_end_maxima();
        self.chrom_lookup = chrom_alias_map(self.intervals.keys());
    }

    /// Resolve a query chromosome name to its key in `intervals`, tolerating
    /// `chr*` vs bare naming and mitochondrial aliases. `None` means the index
    /// has no contig by any accepted spelling of that name.
    pub fn resolve_chrom(&self, chrom: &str) -> Option<&str> {
        self.chrom_lookup.get(chrom).map(String::as_str)
    }

    /// Recompute the per-chromosome prefix maxima of `end` from `intervals`.
    /// Must run after any mutation of `intervals`, and only once they are
    /// sorted by `start`.
    fn rebuild_end_maxima(&mut self) {
        self.end_maxima.clear();
        self.end_maxima.reserve(self.intervals.len());
        for (chrom, intervals) in &self.intervals {
            let mut maxima = Vec::with_capacity(intervals.len());
            let mut running = 0u32;
            for iv in intervals {
                running = running.max(iv.end);
                maxima.push(running);
            }
            self.end_maxima.insert(chrom.clone(), maxima);
        }
    }

    /// The half-open `[lo, hi)` slice of `intervals[chrom]` that can contain an
    /// overlap of `[query_start, query_end]`. See
    /// [`find_overlapping`](Self::find_overlapping) for the reasoning; split out
    /// so the window width itself is testable rather than only its output.
    ///
    /// Resolves `chrom` itself rather than requiring a pre-resolved name, so it
    /// is correct in isolation; resolution is idempotent for a canonical name,
    /// so the caller re-resolving costs one hash lookup.
    fn candidate_window(&self, chrom: &str, query_start: u32, query_end: u32) -> (usize, usize) {
        let Some(chrom) = self.resolve_chrom(chrom) else {
            return (0, 0);
        };
        let intervals = match self.intervals.get(chrom) {
            Some(i) => i,
            None => return (0, 0),
        };

        // One past the last interval that begins at or before the query end.
        let hi = intervals.partition_point(|iv| iv.start <= query_end);

        // First interval that can reach back to `query_start`. Falls back to 0
        // when the derived maxima are absent or stale (an `IntervalIndex`
        // mutated via `add` without a following `sort`), which only costs the
        // old scan width and never changes the answer.
        let lo = match self.end_maxima.get(chrom) {
            Some(maxima) if maxima.len() == intervals.len() => {
                maxima.partition_point(|&max_end| max_end < query_start)
            }
            _ => 0,
        };

        (lo.min(hi), hi)
    }

    /// Find all intervals overlapping [query_start, query_end].
    ///
    /// Two binary searches bound the scan. A point query used to walk every
    /// interval on the chromosome starting at or before it, so cost grew with
    /// how far into the chromosome the query landed rather than with how many
    /// intervals actually overlap it.
    ///
    /// Measured on 620,754 real chr22 coordinates as 2 kb intervals, querying
    /// 15,437 real chr22 variants spread across the chromosome: 70.6 µs →
    /// 48.6 µs per query (1.45×). The gain is modest there because those
    /// intervals overlap heavily — ~25 contain any given point, so the window is
    /// legitimately wide and cloning those hits dominates. It grows with
    /// sparsity: `candidate_window`'s test builds a non-overlapping database and
    /// the window shrinks from the whole prefix to ≤4 intervals.
    ///
    /// * `hi` — one past the last interval with `start <= query_end`, found by
    ///   `partition_point` over the start-sorted list.
    /// * `lo` — the first interval that could reach `query_start`, found by
    ///   `partition_point` over the monotone `end_maxima`. Every interval
    ///   before `lo` ends strictly before `query_start`, so none can overlap.
    ///
    /// Only `[lo, hi)` is examined, and the per-interval `end >= query_start`
    /// filter inside that window is unchanged — so results are identical to a
    /// full scan, in the same (start-sorted) order.
    ///
    /// `chrom` may be given in any accepted spelling; it is resolved against the
    /// index's own contig set via [`resolve_chrom`](Self::resolve_chrom).
    pub fn find_overlapping(&self, chrom: &str, query_start: u32, query_end: u32) -> Vec<OverlapResult> {
        let Some(chrom) = self.resolve_chrom(chrom) else {
            return Vec::new();
        };
        let intervals = match self.intervals.get(chrom) {
            Some(i) => i,
            None => return Vec::new(),
        };
        let (lo, hi) = self.candidate_window(chrom, query_start, query_end);

        let mut results = Vec::new();
        if lo >= hi {
            return results;
        }
        for interval in &intervals[lo..hi] {
            if interval.end >= query_start {
                // Compute reciprocal overlap
                let overlap_start = interval.start.max(query_start);
                let overlap_end = interval.end.min(query_end);
                let overlap_len = (overlap_end as f64 - overlap_start as f64 + 1.0).max(0.0);
                let query_len = (query_end as f64 - query_start as f64 + 1.0).max(1.0);
                let interval_len = (interval.end as f64 - interval.start as f64 + 1.0).max(1.0);

                results.push(OverlapResult {
                    json: interval.json.clone(),
                    reciprocal_overlap: overlap_len / query_len.max(interval_len),
                    annotation_overlap: overlap_len / interval_len,
                });
            }
        }
        results
    }

    /// Serialize to a writer.
    pub fn write_to<W: Write>(&self, writer: &mut W) -> Result<()> {
        writer.write_all(OSI_MAGIC)?;
        writer.write_all(&SCHEMA_VERSION.to_le_bytes())?;
        let data = bincode::serialize(self)?;
        writer.write_all(&(data.len() as u64).to_le_bytes())?;
        writer.write_all(&data)?;
        Ok(())
    }

    /// Deserialize from a reader.
    pub fn read_from<R: Read>(reader: &mut R) -> Result<Self> {
        let mut magic = [0u8; 8];
        reader.read_exact(&mut magic)?;
        if &magic != OSI_MAGIC {
            anyhow::bail!("Invalid OSI magic");
        }
        let mut ver = [0u8; 2];
        reader.read_exact(&mut ver)?;
        if u16::from_le_bytes(ver) != SCHEMA_VERSION {
            anyhow::bail!("Unsupported OSI schema version");
        }
        let mut len_bytes = [0u8; 8];
        reader.read_exact(&mut len_bytes)?;
        let len_u64 = u64::from_le_bytes(len_bytes);
        if len_u64 > MAX_INDEX_PAYLOAD {
            anyhow::bail!(
                "OSI payload size {} exceeds limit {}",
                len_u64,
                MAX_INDEX_PAYLOAD
            );
        }
        let len: usize = len_u64
            .try_into()
            .map_err(|_| anyhow::anyhow!("OSI payload size {} exceeds usize", len_u64))?;
        let mut data = vec![0u8; len];
        reader.read_exact(&mut data)?;
        let mut index: IntervalIndex = bincode::deserialize(&data)?;
        // `find_overlapping` relies on per-chromosome intervals being sorted
        // by `start`. The writer enforces this via `sort()`, but a hand-
        // crafted or partially-corrupted .osi could violate the invariant
        // and produce silent missed overlaps. Sort on read as a safety net;
        // this is O(n log n) once per open and negligible at query time.
        let mut needed_sort = false;
        for intervals in index.intervals.values_mut() {
            if !intervals.windows(2).all(|w| w[0].start <= w[1].start) {
                intervals.sort_by_key(|i| i.start);
                needed_sort = true;
            }
        }
        if needed_sort {
            log::warn!(
                "OSI '{}': intervals were not stored in sorted order; \
                 sorted on load (writer should have called .sort())",
                index.header.name
            );
        }
        // The derived fields are `#[serde(skip)]`, so they arrive empty
        // regardless of what the file contained. Build them now, after the sort
        // above, so the very first query already gets the binary-searched path.
        index.rebuild_derived();
        Ok(index)
    }
}

/// Annotation-pipeline-facing wrapper around a loaded `.osi`. Holds the
/// SaMetadata view that the annotate pipeline expects from every SA
/// provider, plus the `IntervalIndex` itself for overlap queries.
pub struct OsiReader {
    pub index: IntervalIndex,
    metadata: SaMetadata,
}

impl OsiReader {
    /// Open a `.osi` file and wrap it as an `AnnotationProvider`.
    pub fn open(path: &Path) -> Result<Self> {
        let mut file = std::fs::File::open(path)?;
        let mut index = IntervalIndex::read_from(&mut file)?;
        // `find_overlapping` binary-searches both the start-sorted interval
        // list and the derived `end_maxima`, so the sort invariant is
        // load-bearing. `IntervalIndex::read_from` already sorts on read as a
        // safety net (see its docs); this second `sort()` is cheap on
        // already-sorted input and keeps the invariant explicit at the one
        // place the annotate pipeline enters from.
        index.sort();
        let metadata = SaMetadata {
            name: index.header.name.clone(),
            version: index.header.version.clone(),
            description: format!(
                "Interval annotation database from {}",
                path.display()
            ),
            assembly: index.header.assembly.clone(),
            json_key: index.header.json_key.clone(),
            // Intervals are inherently positional — overlap doesn't care
            // about REF/ALT. Setting these correctly lets the runtime
            // dispatch in the annotate pipeline skip allele matching.
            match_by_allele: false,
            is_array: true,
            is_positional: true,
        };
        Ok(Self { index, metadata })
    }
}

impl AnnotationProvider for OsiReader {
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
        _ref_allele: &str,
        _alt_allele: &str,
    ) -> Result<Option<AnnotationValue>> {
        // Point-query semantics: report every interval that contains `pos`.
        // `find_overlapping` resolves the chromosome through the index's alias
        // map, so a BED stored as `chrM` matches a VCF query of `MT` (and vice
        // versa), matching the behaviour of `SaReader`.
        let pos32: u32 = match u32::try_from(pos) {
            Ok(p) => p,
            Err(_) => return Ok(None),
        };
        let hits = self.index.find_overlapping(chrom, pos32, pos32);
        if hits.is_empty() {
            return Ok(None);
        }
        Ok(Some(AnnotationValue::Interval(
            hits.into_iter().map(|h| h.json).collect(),
        )))
    }
}

/// Result of an overlap query.
#[derive(Debug, Clone)]
pub struct OverlapResult {
    pub json: String,
    /// Overlap as fraction of the larger region.
    pub reciprocal_overlap: f64,
    /// Overlap as fraction of the annotation interval.
    pub annotation_overlap: f64,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_interval_round_trip() {
        let header = IntervalHeader {
            schema_version: SCHEMA_VERSION,
            json_key: "dgv".into(),
            name: "DGV".into(),
            version: "1.0".into(),
            assembly: "GRCh38".into(),
        };

        let mut index = IntervalIndex::new(header);
        index.add(IntervalRecord {
            chrom: "chr1".into(),
            start: 100,
            end: 500,
            json: r#"{"type":"DEL"}"#.into(),
        });
        index.add(IntervalRecord {
            chrom: "chr1".into(),
            start: 300,
            end: 800,
            json: r#"{"type":"DUP"}"#.into(),
        });
        index.sort();

        // Serialize and deserialize
        let mut buf = Vec::new();
        index.write_to(&mut buf).unwrap();
        let loaded = IntervalIndex::read_from(&mut std::io::Cursor::new(buf)).unwrap();

        assert_eq!(loaded.header.json_key, "dgv");
        assert_eq!(loaded.intervals["chr1"].len(), 2);
    }

    #[test]
    fn test_read_sorts_unsorted_intervals() {
        // Build an index with intervals deliberately out of order, serialize
        // it (skipping the writer's `sort()`), then verify read_from puts
        // them back in start-order so find_overlapping is correct.
        let header = IntervalHeader {
            schema_version: SCHEMA_VERSION,
            json_key: "test".into(),
            name: "Unsorted".into(),
            version: "1.0".into(),
            assembly: "GRCh38".into(),
        };
        let mut index = IntervalIndex::new(header);
        // Intentionally out of order:
        index.add(IntervalRecord { chrom: "chr1".into(), start: 500, end: 700, json: "{\"id\":\"B\"}".into() });
        index.add(IntervalRecord { chrom: "chr1".into(), start: 100, end: 200, json: "{\"id\":\"A\"}".into() });
        // Skip index.sort() so the on-disk layout is unsorted.

        let mut buf = Vec::new();
        index.write_to(&mut buf).unwrap();
        let loaded = IntervalIndex::read_from(&mut std::io::Cursor::new(buf)).unwrap();

        let intervals = &loaded.intervals["chr1"];
        assert_eq!(intervals[0].start, 100);
        assert_eq!(intervals[1].start, 500);

        // And find_overlapping returns both when the query spans them.
        let hits = loaded.find_overlapping("chr1", 50, 800);
        assert_eq!(hits.len(), 2);
    }

    #[test]
    fn test_find_overlapping() {
        let header = IntervalHeader {
            schema_version: SCHEMA_VERSION,
            json_key: "test".into(),
            name: "Test".into(),
            version: "1.0".into(),
            assembly: "GRCh38".into(),
        };

        let mut index = IntervalIndex::new(header);
        index.add(IntervalRecord {
            chrom: "chr1".into(),
            start: 100,
            end: 500,
            json: r#"{"id":"A"}"#.into(),
        });
        index.add(IntervalRecord {
            chrom: "chr1".into(),
            start: 400,
            end: 800,
            json: r#"{"id":"B"}"#.into(),
        });
        index.add(IntervalRecord {
            chrom: "chr1".into(),
            start: 1000,
            end: 1500,
            json: r#"{"id":"C"}"#.into(),
        });
        index.sort();

        // Query overlapping A and B
        let results = index.find_overlapping("chr1", 300, 600);
        assert_eq!(results.len(), 2);

        // Query overlapping only C
        let results = index.find_overlapping("chr1", 1200, 1300);
        assert_eq!(results.len(), 1);
        assert!(results[0].json.contains("\"C\""));

        // No overlap
        let results = index.find_overlapping("chr1", 900, 950);
        assert_eq!(results.len(), 0);
    }

    #[test]
    fn osi_reader_round_trip_and_provider_query() {
        // Write a `.osi`, reopen it as `OsiReader`, and exercise the
        // `AnnotationProvider` query path that the runtime annotate
        // pipeline actually uses.
        let header = IntervalHeader {
            schema_version: SCHEMA_VERSION,
            json_key: "myregions".into(),
            name: "MyRegions".into(),
            version: "1.0".into(),
            assembly: "GRCh38".into(),
        };
        let mut index = IntervalIndex::new(header);
        index.add(IntervalRecord {
            chrom: "chr1".into(),
            start: 100,
            end: 200,
            json: r#"{"name":"alpha"}"#.into(),
        });
        index.add(IntervalRecord {
            chrom: "chr1".into(),
            start: 150,
            end: 300,
            json: r#"{"name":"beta"}"#.into(),
        });
        index.sort();

        let tmp = tempfile::NamedTempFile::new().unwrap();
        let path = tmp.path().with_extension("osi");
        let mut file = std::fs::File::create(&path).unwrap();
        index.write_to(&mut file).unwrap();
        drop(file);

        let reader = OsiReader::open(&path).unwrap();
        assert_eq!(reader.json_key(), "myregions");
        assert!(reader.metadata().is_positional);

        // Position inside both intervals → returns both JSONs.
        let val = reader.annotate_position("chr1", 175, "", "").unwrap();
        match val {
            Some(AnnotationValue::Interval(v)) => {
                assert_eq!(v.len(), 2);
                assert!(v.iter().any(|s| s.contains("alpha")));
                assert!(v.iter().any(|s| s.contains("beta")));
            }
            other => panic!("expected Interval value, got {:?}", other),
        }

        // Bare-chromosome input matches the chr-prefixed BED-style storage.
        let val = reader.annotate_position("1", 175, "", "").unwrap();
        assert!(val.is_some(), "chr-prefix normalization should match");

        // Position outside any interval → None.
        let val = reader.annotate_position("chr1", 50, "", "").unwrap();
        assert!(val.is_none());
    }

    /// Reference implementation: the pre-optimization full scan. Any
    /// divergence between this and `find_overlapping` is a correctness bug.
    fn brute_force(
        intervals: &[StoredInterval],
        query_start: u32,
        query_end: u32,
    ) -> Vec<String> {
        intervals
            .iter()
            .filter(|iv| iv.start <= query_end && iv.end >= query_start)
            .map(|iv| iv.json.clone())
            .collect()
    }

    /// Deterministic xorshift so the differential test is reproducible without
    /// pulling in a `rand` dependency.
    fn next_rand(state: &mut u64) -> u64 {
        *state ^= *state << 13;
        *state ^= *state >> 7;
        *state ^= *state << 17;
        *state
    }

    fn test_header() -> IntervalHeader {
        IntervalHeader {
            schema_version: SCHEMA_VERSION,
            json_key: "t".into(),
            name: "T".into(),
            version: "1.0".into(),
            assembly: "GRCh38".into(),
        }
    }

    #[test]
    fn find_overlapping_matches_full_scan_on_randomized_intervals() {
        // Covers the shapes the binary-searched window has to survive:
        // duplicate starts, fully nested intervals, one very long interval
        // that reaches far past everything after it, and zero-length spans.
        let mut state = 0x243F_6A88_85A3_08D3u64;
        let mut index = IntervalIndex::new(test_header());
        for i in 0..600u32 {
            let start = (next_rand(&mut state) % 10_000) as u32;
            let len = match i % 7 {
                0 => 0,                                      // point interval
                1 => 5_000 + (next_rand(&mut state) % 4_000) as u32, // long reacher
                _ => (next_rand(&mut state) % 50) as u32,
            };
            index.add(IntervalRecord {
                chrom: "chr1".into(),
                start,
                end: start.saturating_add(len),
                json: format!("{{\"id\":{}}}", i),
            });
        }
        index.sort();
        let intervals = index.intervals["chr1"].clone();

        for q in 0..12_000u32 {
            // Point queries — the hot path from `annotate_position`.
            let got: Vec<String> = index
                .find_overlapping("chr1", q, q)
                .into_iter()
                .map(|h| h.json)
                .collect();
            let want = brute_force(&intervals, q, q);
            assert_eq!(got, want, "point query at {} diverged", q);
        }

        // Range queries, including ones that start past every interval.
        for _ in 0..2_000 {
            let a = (next_rand(&mut state) % 14_000) as u32;
            let b = a + (next_rand(&mut state) % 3_000) as u32;
            let got: Vec<String> = index
                .find_overlapping("chr1", a, b)
                .into_iter()
                .map(|h| h.json)
                .collect();
            assert_eq!(got, brute_force(&intervals, a, b), "range {}..{} diverged", a, b);
        }
    }

    #[test]
    fn candidate_window_does_not_scan_from_index_zero() {
        // The regression this guards: a point query near the end of a
        // chromosome used to walk every interval starting at or before it, so
        // cost grew with position rather than with overlap count. With the
        // `end_maxima` bound the window stays proportional to local density.
        let mut index = IntervalIndex::new(test_header());
        for i in 0..100_000u32 {
            let start = i * 10;
            index.add(IntervalRecord {
                chrom: "chr1".into(),
                start,
                end: start + 9,
                json: format!("{{\"id\":{}}}", i),
            });
        }
        index.sort();

        let (lo, hi) = index.candidate_window("chr1", 999_000, 999_000);
        assert!(
            hi - lo <= 4,
            "expected a tiny candidate window for a point query, got {}..{} ({} intervals)",
            lo,
            hi,
            hi - lo
        );
        assert!(lo > 99_000, "lower bound should skip the whole prefix, got {}", lo);

        // And it still finds the right interval.
        let hits = index.find_overlapping("chr1", 999_000, 999_000);
        assert_eq!(hits.len(), 1);
        assert!(hits[0].json.contains("99900"));
    }

    #[test]
    fn candidate_window_falls_back_to_full_scan_when_maxima_are_stale() {
        // `add` after `sort` deliberately drops that chromosome's maxima. The
        // window must widen to the whole prefix rather than silently trusting
        // a stale array and missing overlaps.
        let mut index = IntervalIndex::new(test_header());
        index.add(IntervalRecord { chrom: "chr1".into(), start: 10, end: 20, json: "{\"id\":\"a\"}".into() });
        index.add(IntervalRecord { chrom: "chr1".into(), start: 30, end: 40, json: "{\"id\":\"b\"}".into() });
        index.sort();
        assert_eq!(index.candidate_window("chr1", 35, 35).0, 1, "maxima should bound a fresh index");

        index.add(IntervalRecord { chrom: "chr1".into(), start: 50, end: 60, json: "{\"id\":\"c\"}".into() });
        assert_eq!(
            index.candidate_window("chr1", 35, 35).0,
            0,
            "a mutated index must fall back to scanning from 0"
        );
        // Correctness is preserved either way.
        let hits = index.find_overlapping("chr1", 35, 35);
        assert_eq!(hits.len(), 1);
        assert!(hits[0].json.contains("\"b\""));
    }

    #[test]
    fn osi_files_written_before_end_maxima_existed_still_load() {
        // `end_maxima` is `#[serde(skip)]`, so an `.osi` is still encoded as
        // exactly (header, intervals). Build that payload from a standalone
        // tuple — i.e. the old two-field layout — and require `read_from` to
        // accept it and to populate the maxima itself.
        let mut intervals: HashMap<String, Vec<StoredInterval>> = HashMap::new();
        intervals.insert(
            "chr1".into(),
            vec![
                StoredInterval { start: 100, end: 200, json: "{\"id\":\"a\"}".into() },
                StoredInterval { start: 150, end: 400, json: "{\"id\":\"b\"}".into() },
            ],
        );
        let legacy = (test_header(), intervals);

        let mut buf = Vec::new();
        buf.extend_from_slice(OSI_MAGIC);
        buf.extend_from_slice(&SCHEMA_VERSION.to_le_bytes());
        let payload = bincode::serialize(&legacy).unwrap();
        buf.extend_from_slice(&(payload.len() as u64).to_le_bytes());
        buf.extend_from_slice(&payload);

        let loaded = IntervalIndex::read_from(&mut std::io::Cursor::new(buf)).unwrap();
        assert_eq!(loaded.intervals["chr1"].len(), 2);
        assert_eq!(loaded.end_maxima["chr1"], vec![200, 400]);
        // Both derived fields are rebuilt, so the bounded path and the alias
        // resolution are live immediately after load.
        assert_eq!(loaded.resolve_chrom("1"), Some("chr1"));
        assert_eq!(loaded.candidate_window("chr1", 350, 350), (1, 2));
        assert_eq!(loaded.find_overlapping("chr1", 350, 350).len(), 1);
        assert_eq!(loaded.find_overlapping("1", 350, 350).len(), 1);
    }

    #[test]
    fn find_overlapping_accepts_any_chromosome_spelling() {
        // The alias map replaced a per-query `chrom_aliases` call; the accepted
        // spellings must be exactly the same set as before.
        let mut index = IntervalIndex::new(test_header());
        index.add(IntervalRecord { chrom: "chrM".into(), start: 100, end: 200, json: "{}".into() });
        index.add(IntervalRecord { chrom: "1".into(), start: 10, end: 20, json: "{}".into() });
        index.sort();

        for alias in ["chrM", "M", "MT", "chrMT"] {
            assert_eq!(
                index.find_overlapping(alias, 150, 150).len(),
                1,
                "mito spelling {alias} should resolve"
            );
        }
        // Bare-keyed contig reached from the chr-prefixed spelling and back.
        assert_eq!(index.find_overlapping("1", 15, 15).len(), 1);
        assert_eq!(index.find_overlapping("chr1", 15, 15).len(), 1);
        // A contig the index does not hold resolves to nothing.
        assert!(index.find_overlapping("chr9", 15, 15).is_empty());
        assert_eq!(index.resolve_chrom("chr9"), None);
    }

    #[test]
    fn add_registers_aliases_before_sort() {
        // `add` extends the alias map eagerly, so an index that was populated
        // but never sorted still resolves spellings (it falls back to the wide
        // scan window, which `candidate_window`'s own test covers).
        let mut index = IntervalIndex::new(test_header());
        index.add(IntervalRecord { chrom: "chr1".into(), start: 10, end: 20, json: "{}".into() });
        assert_eq!(index.resolve_chrom("1"), Some("chr1"));
        assert_eq!(index.find_overlapping("1", 15, 15).len(), 1);
    }

    #[test]
    fn write_to_output_is_unchanged_by_the_derived_field() {
        // Guards the on-disk format explicitly: the bytes `write_to` emits must
        // equal magic + version + len + bincode((header, intervals)).
        let mut index = IntervalIndex::new(test_header());
        index.add(IntervalRecord { chrom: "chr1".into(), start: 5, end: 9, json: "{}".into() });
        index.sort();

        let mut actual = Vec::new();
        index.write_to(&mut actual).unwrap();

        let legacy = (index.header.clone(), index.intervals.clone());
        let payload = bincode::serialize(&legacy).unwrap();
        let mut expected = Vec::new();
        expected.extend_from_slice(OSI_MAGIC);
        expected.extend_from_slice(&SCHEMA_VERSION.to_le_bytes());
        expected.extend_from_slice(&(payload.len() as u64).to_le_bytes());
        expected.extend_from_slice(&payload);

        assert_eq!(actual, expected);
    }

    #[test]
    fn osi_reader_matches_mito_chromosome_aliases() {
        // Index stores `chrM` (UCSC). Query with NCBI-style `MT` and bare
        // `M` must still match via the shared `chrom_aliases` helper.
        let header = IntervalHeader {
            schema_version: SCHEMA_VERSION,
            json_key: "mito".into(),
            name: "Mito".into(),
            version: "1.0".into(),
            assembly: "GRCh38".into(),
        };
        let mut index = IntervalIndex::new(header);
        index.add(IntervalRecord {
            chrom: "chrM".into(),
            start: 1000,
            end: 2000,
            json: r#"{"name":"mito_region"}"#.into(),
        });
        index.sort();

        let tmp = tempfile::NamedTempFile::new().unwrap();
        let path = tmp.path().with_extension("osi");
        let mut file = std::fs::File::create(&path).unwrap();
        index.write_to(&mut file).unwrap();
        drop(file);

        let reader = OsiReader::open(&path).unwrap();
        for alias in &["chrM", "M", "MT", "chrMT"] {
            let val = reader.annotate_position(alias, 1500, "", "").unwrap();
            assert!(val.is_some(), "mito query '{}' should match", alias);
        }
    }
}
