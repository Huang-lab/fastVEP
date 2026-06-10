//! Interval-based annotation reader/writer (.osi files).
//!
//! Used for structural variant databases (gnomAD SV, ClinGen dosage, DGV)
//! where annotations are regions rather than point positions.

use crate::common::{IntervalRecord, MAX_INDEX_PAYLOAD, OSI_MAGIC, SCHEMA_VERSION};
use anyhow::Result;
use fastvep_cache::annotation::{AnnotationProvider, AnnotationValue, SaMetadata};
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
        }
    }

    /// Add an interval record.
    pub fn add(&mut self, record: IntervalRecord) {
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
    pub fn sort(&mut self) {
        for intervals in self.intervals.values_mut() {
            intervals.sort_by_key(|i| i.start);
        }
    }

    /// Find all intervals overlapping [query_start, query_end].
    pub fn find_overlapping(&self, chrom: &str, query_start: u32, query_end: u32) -> Vec<OverlapResult> {
        let intervals = match self.intervals.get(chrom) {
            Some(i) => i,
            None => return Vec::new(),
        };

        // Binary search for first interval with start <= query_end
        // (since intervals are sorted by start, we need those where start <= query_end)
        let mut results = Vec::new();
        for interval in intervals {
            if interval.start > query_end {
                break; // Past the query range
            }
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
        // `find_overlapping` does a linear scan from the start of the
        // chromosome's interval vector, so the sort invariant is
        // load-bearing — `IntervalIndex::read_from` already sorts on read
        // as a safety net (see its docs), but if any caller mutates after
        // load they're responsible for re-sorting.
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
        // Point-query semantics: report every interval that contains
        // `pos`. We probe every standard alias for the chromosome so a
        // BED stored as `chrM` matches a VCF query of `MT` (and vice
        // versa), matching the behaviour of `SaReader`.
        let pos32: u32 = match u32::try_from(pos) {
            Ok(p) => p,
            Err(_) => return Ok(None),
        };
        let mut json_hits: Vec<String> = Vec::new();
        for alias in crate::common::chrom_aliases(chrom) {
            let hits = self.index.find_overlapping(&alias, pos32, pos32);
            if !hits.is_empty() {
                json_hits.extend(hits.into_iter().map(|h| h.json));
                // A `.osi` only stores each chromosome under one canonical
                // name, so the first alias that matches is the right one;
                // stop here to avoid double-counting if a (corrupted) file
                // somehow had two names for one contig.
                break;
            }
        }
        if json_hits.is_empty() {
            return Ok(None);
        }
        Ok(Some(AnnotationValue::Interval(json_hits)))
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
        assert_eq!(reader.metadata().is_positional, true);

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
