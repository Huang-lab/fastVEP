//! Genomic chunk: in-memory data structure for a ~1MB region.
//!
//! Each chunk contains sorted Var32 keys and parallel value arrays.
//! Chunks are loaded on demand from .osa2 ZIP archives.

use crate::fields::{Field, FieldType};
use crate::kmer16::LongVariant;

/// A loaded genomic chunk (~1MB region) with sorted variant keys and values.
pub struct Chunk {
    /// Sorted Var32-encoded variant keys for binary search.
    pub var32s: Vec<u32>,
    /// Long variants (ref+alt > 4 bases) sorted for binary search.
    pub longs: Vec<LongVariant>,
    /// Parallel value arrays: `values[field_idx][variant_idx]`.
    pub values: Vec<Vec<u32>>,
    /// Optional JSON blob strings for JsonBlob fields.
    pub json_blobs: Option<Vec<String>>,
}

impl Chunk {
    /// Create an empty chunk.
    pub fn empty() -> Self {
        Self { var32s: Vec::new(), longs: Vec::new(), values: Vec::new(), json_blobs: None }
    }

    /// Number of variants in this chunk.
    pub fn len(&self) -> usize {
        self.var32s.len()
    }

    pub fn is_empty(&self) -> bool {
        self.var32s.is_empty()
    }

    /// Look up a variant by Var32 key. Returns the index into value arrays.
    #[inline]
    pub fn find_short(&self, encoded: u32) -> Option<usize> {
        self.var32s.binary_search(&encoded).ok()
    }

    /// Look up a long variant. Returns the index into value arrays.
    pub fn find_long(&self, position: u32, ref_allele: &[u8], alt_allele: &[u8]) -> Option<usize> {
        let query = LongVariant {
            position,
            idx: 0,
            sequence: crate::kmer16::encode_var(ref_allele, alt_allele),
        };
        self.longs.binary_search(&query).ok().map(|i| self.longs[i].idx as usize)
    }

    /// Reconstruct a JSON string from parallel value arrays at the given index.
    ///
    /// `values` is parallel only to non-JsonBlob fields (in field-config order),
    /// so we maintain a separate `value_idx` that advances only for those.
    /// `strings` is parallel to all fields, so it uses the field index.
    pub fn reconstruct_json(
        &self,
        idx: usize,
        fields: &[Field],
        strings: &[Vec<String>],
    ) -> String {
        let mut parts = Vec::with_capacity(fields.len());
        let mut value_idx: usize = 0;

        for (fi, field) in fields.iter().enumerate() {
            if field.ftype == FieldType::JsonBlob {
                if let Some(ref blobs) = self.json_blobs {
                    if idx < blobs.len() && !blobs[idx].is_empty() {
                        parts.push(format!("\"{}\":{}", field.alias, blobs[idx]));
                    }
                }
                continue;
            }

            let column = match self.values.get(value_idx) {
                Some(c) => c,
                None => {
                    value_idx += 1;
                    continue;
                }
            };
            value_idx += 1;

            let stored = match column.get(idx) {
                Some(&s) => s,
                None => continue,
            };
            if stored == field.missing_value {
                continue; // Skip missing values in output
            }

            let val_str = crate::fields::format_value(field, stored, strings.get(fi).map(|v| v.as_slice()));
            if val_str != "null" {
                parts.push(format!("\"{}\":{}", field.alias, val_str));
            }
        }

        format!("{{{}}}", parts.join(","))
    }
}

/// Delta-encode a sorted u32 array in place. Returns the encoded array.
pub fn delta_encode(values: &[u32]) -> Vec<u32> {
    if values.is_empty() {
        return Vec::new();
    }
    let mut encoded = Vec::with_capacity(values.len());
    encoded.push(values[0]);
    for i in 1..values.len() {
        encoded.push(values[i].wrapping_sub(values[i - 1]));
    }
    encoded
}

/// Delta-decode a u32 array (cumulative sum). Modifies in place.
pub fn delta_decode(encoded: &mut [u32]) {
    for i in 1..encoded.len() {
        encoded[i] = encoded[i].wrapping_add(encoded[i - 1]);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::var32;

    #[test]
    fn test_delta_round_trip() {
        let original = vec![100, 105, 110, 200, 300];
        let encoded = delta_encode(&original);
        assert_eq!(encoded, vec![100, 5, 5, 90, 100]);

        let mut decoded = encoded;
        delta_decode(&mut decoded);
        assert_eq!(decoded, original);
    }

    #[test]
    fn test_delta_empty() {
        let encoded = delta_encode(&[]);
        assert!(encoded.is_empty());
    }

    #[test]
    fn test_delta_single() {
        let original = vec![42];
        let encoded = delta_encode(&original);
        assert_eq!(encoded, vec![42]);
        let mut decoded = encoded;
        delta_decode(&mut decoded);
        assert_eq!(decoded, original);
    }

    #[test]
    fn test_chunk_find_short() {
        let mut chunk = Chunk::empty();
        // Encode some variants
        let keys: Vec<u32> = (0..100)
            .filter_map(|i| var32::encode(i * 10, b"A", b"G"))
            .collect();
        chunk.var32s = keys;

        // Should find position 50 (i=5, pos=50)
        let query = var32::encode(50, b"A", b"G").unwrap();
        assert!(chunk.find_short(query).is_some());

        // Should NOT find position 55 (not in our set)
        let query = var32::encode(55, b"A", b"G").unwrap();
        assert!(chunk.find_short(query).is_none());
    }

    #[test]
    fn test_chunk_reconstruct_json_with_jsonblob_in_middle() {
        // Regression: previously `values[fi]` indexed by the field-config
        // position, which silently dropped any non-JsonBlob field that came
        // after a JsonBlob field. Verify the trailing Integer is emitted.
        let fields = vec![
            Field {
                field: "AF".into(), alias: "af".into(), ftype: FieldType::Float,
                multiplier: 1_000_000, zigzag: false, missing_value: u32::MAX,
                missing_string: ".".into(), description: String::new(),
            },
            Field {
                field: "blob".into(), alias: "blob".into(), ftype: FieldType::JsonBlob,
                multiplier: 1, zigzag: false, missing_value: u32::MAX,
                missing_string: ".".into(), description: String::new(),
            },
            Field {
                field: "AC".into(), alias: "ac".into(), ftype: FieldType::Integer,
                multiplier: 1, zigzag: false, missing_value: u32::MAX,
                missing_string: ".".into(), description: String::new(),
            },
        ];

        let mut chunk = Chunk::empty();
        chunk.var32s = vec![var32::encode(100, b"A", b"G").unwrap()];
        // Two non-JsonBlob columns, in field order: AF then AC.
        chunk.values = vec![vec![1234], vec![42]];
        chunk.json_blobs = Some(vec![r#"{"k":1}"#.to_string()]);

        let json = chunk.reconstruct_json(0, &fields, &[]);
        assert!(json.contains("\"af\":"), "missing af in: {}", json);
        assert!(json.contains("\"ac\":42"), "missing ac in: {}", json);
        assert!(json.contains("\"blob\":{\"k\":1}"), "missing blob in: {}", json);
    }

    #[test]
    fn test_chunk_reconstruct_json() {
        let fields = vec![
            Field {
                field: "AF".into(), alias: "allAf".into(), ftype: FieldType::Float,
                multiplier: 1_000_000, zigzag: false, missing_value: u32::MAX,
                missing_string: ".".into(), description: String::new(),
            },
            Field {
                field: "AC".into(), alias: "allAc".into(), ftype: FieldType::Integer,
                multiplier: 1, zigzag: false, missing_value: u32::MAX,
                missing_string: ".".into(), description: String::new(),
            },
        ];

        let mut chunk = Chunk::empty();
        chunk.var32s = vec![var32::encode(100, b"A", b"G").unwrap()];
        chunk.values = vec![
            vec![1234],   // AF * 1_000_000 = 0.001234
            vec![42],     // AC = 42
        ];

        let json = chunk.reconstruct_json(0, &fields, &[]);
        assert!(json.contains("\"allAf\":"));
        assert!(json.contains("\"allAc\":42"));
    }
}
