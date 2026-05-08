//! Interval-based annotation reader/writer (.osi files).
//!
//! Used for structural variant databases (gnomAD SV, ClinGen dosage, DGV)
//! where annotations are regions rather than point positions.

use crate::common::{IntervalRecord, MAX_INDEX_PAYLOAD, OSI_MAGIC, SCHEMA_VERSION};
use anyhow::Result;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::io::{Read, Write};

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
        if len_u64 > MAX_INDEX_PAYLOAD as u64 {
            anyhow::bail!(
                "OSI payload size {} exceeds limit {}",
                len_u64,
                MAX_INDEX_PAYLOAD
            );
        }
        let mut data = vec![0u8; len_u64 as usize];
        reader.read_exact(&mut data)?;
        let index: IntervalIndex = bincode::deserialize(&data)?;
        Ok(index)
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
}
