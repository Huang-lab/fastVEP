//! Genomic chunk: in-memory data structure for a ~1MB region.
//!
//! Each chunk contains sorted Var32 keys and parallel value arrays.
//! Chunks are loaded on demand from .osa2 ZIP archives.

use crate::fields::{Field, FieldType};
use crate::kmer16::{LongVariant, OtherVariant};

/// The JSON-blob payloads of a chunk: one buffer holding every row's text, plus
/// the offsets that delimit the rows.
///
/// These were a `Vec<String>` - one heap allocation and one copy per row. The
/// cost is paid per *stored* row, not per queried one, and a chunk holds every
/// row in its megabase of the genome. For a sparse source that is a handful; for
/// the densest ones it is the whole neighbourhood. SpliceAI scores three
/// alternates at essentially every position of a gene body, which puts on the
/// order of a million records in one chunk, so a single lookup allocated a
/// million strings, copied the entire decompressed blob into them, read one, and
/// dropped the rest (#101). Measured on a 3M-row-per-chunk source, 300 scattered
/// queries took 5.9s and 3.8 GB of resident memory.
///
/// Indexing the buffer instead costs eight bytes a row and no allocation, and
/// the text is kept exactly as it was decompressed, so what a query returns is
/// byte-for-byte what it returned before.
#[derive(Debug, Default, Clone)]
pub struct JsonBlobs {
    text: String,
    /// `(start, end)` byte offsets into `text`, one per row.
    ///
    /// `u32` is sound because `reader_v2::MAX_JSON_BLOB_DECOMPRESSED` rejects a
    /// chunk blob well before 4 GiB, and that check runs before this is built.
    rows: Vec<(u32, u32)>,
}

impl JsonBlobs {
    /// Index the newline-separated rows of a decompressed blob entry.
    ///
    /// Row boundaries come from [`str::lines`], so a `\r\n` archive is split the
    /// same way the `Vec<String>` build split it.
    pub fn from_text(text: String) -> Self {
        debug_assert!(
            text.len() <= u32::MAX as usize,
            "blob buffer exceeds the u32 offsets this type stores"
        );
        let base = text.as_ptr() as usize;
        let rows: Vec<(u32, u32)> = text
            .lines()
            .map(|line| {
                // `lines()` yields subslices of `text`, so the difference is an
                // offset into it. Integers only, so the borrow ends here.
                let start = line.as_ptr() as usize - base;
                (start as u32, (start + line.len()) as u32)
            })
            .collect();
        Self { text, rows }
    }

    /// Build from rows already in hand, for callers and tests that have them
    /// separated rather than as one buffer.
    pub fn from_rows<I, S>(rows: I) -> Self
    where
        I: IntoIterator<Item = S>,
        S: AsRef<str>,
    {
        let mut text = String::new();
        for (i, row) in rows.into_iter().enumerate() {
            if i > 0 {
                text.push('\n');
            }
            text.push_str(row.as_ref());
        }
        Self::from_text(text)
    }

    /// Number of rows.
    pub fn len(&self) -> usize {
        self.rows.len()
    }

    pub fn is_empty(&self) -> bool {
        self.rows.is_empty()
    }

    /// The blob for one row, borrowed out of the buffer.
    pub fn get(&self, idx: usize) -> Option<&str> {
        let &(start, end) = self.rows.get(idx)?;
        self.text.get(start as usize..end as usize)
    }

    /// Resident footprint, for the chunk cache's byte budget.
    pub fn heap_bytes(&self) -> usize {
        self.text.len() + self.rows.len() * std::mem::size_of::<(u32, u32)>()
    }
}

/// A loaded genomic chunk (~1MB region) with sorted variant keys and values.
pub struct Chunk {
    /// Sorted Var32-encoded variant keys for binary search.
    pub var32s: Vec<u32>,
    /// Long variants (ref+alt > 4 bases) sorted for binary search.
    pub longs: Vec<LongVariant>,
    /// Variants whose alleles are not 2-bit packable (contain `N`/IUPAC), kept
    /// verbatim and sorted for binary search. Almost always empty; real ClinVar
    /// puts ~0.015% of its records here. Absent from `.osa2` files written
    /// before this bucket existed, which simply leaves it empty.
    pub others: Vec<OtherVariant>,
    /// Parallel value arrays, one per non-JsonBlob field in field-config order:
    /// `values[non_jsonblob_position][variant_idx]`. JsonBlob fields do not
    /// contribute a column here (their payloads live in `json_blobs`), so this
    /// vector is shorter than `fields.len()` whenever any JsonBlob is present.
    /// `reconstruct_json` tracks a separate counter that advances only for
    /// non-JsonBlob fields when indexing into this array.
    pub values: Vec<Vec<u32>>,
    /// Optional JSON blob strings for JsonBlob fields.
    pub json_blobs: Option<JsonBlobs>,
}

impl Chunk {
    /// Create an empty chunk.
    pub fn empty() -> Self {
        Self {
            var32s: Vec::new(),
            longs: Vec::new(),
            others: Vec::new(),
            values: Vec::new(),
            json_blobs: None,
        }
    }

    /// Number of variants in this chunk (short + long + non-ACGT).
    pub fn len(&self) -> usize {
        self.var32s.len() + self.longs.len() + self.others.len()
    }

    /// A chunk is empty only when it holds no variants of any kind. Long-only
    /// chunks (all indels) must not be treated as empty, or every lookup
    /// against them would miss — and the same goes for a chunk that happens to
    /// hold only non-ACGT-allele variants.
    pub fn is_empty(&self) -> bool {
        self.var32s.is_empty() && self.longs.is_empty() && self.others.is_empty()
    }

    /// Look up a variant by Var32 key. Returns the index into value arrays.
    #[inline]
    pub fn find_short(&self, encoded: u32) -> Option<usize> {
        self.var32s.binary_search(&encoded).ok()
    }

    /// Look up a long variant. Returns the index into value arrays.
    ///
    /// Returns `None` if either allele contains a non-ACGT base — such a variant
    /// is filed in `others`, not here; see [`find_other`](Self::find_other).
    pub fn find_long(&self, position: u32, ref_allele: &[u8], alt_allele: &[u8]) -> Option<usize> {
        let sequence = crate::kmer16::encode_var(ref_allele, alt_allele)?;
        let query = LongVariant {
            position,
            idx: 0,
            sequence,
        };
        self.longs
            .binary_search(&query)
            .ok()
            .map(|i| self.longs[i].idx as usize)
    }

    /// Look up a variant whose alleles are not 2-bit packable. Returns the index
    /// into the value arrays.
    ///
    /// `others` is empty for every source without `N`/IUPAC alleles and for any
    /// `.osa2` written before the bucket existed, so this is a length check in
    /// the overwhelmingly common case.
    pub fn find_other(&self, position: u32, ref_allele: &[u8], alt_allele: &[u8]) -> Option<usize> {
        if self.others.is_empty() {
            return None;
        }
        let query = OtherVariant {
            position,
            idx: 0,
            ref_allele: ref_allele.to_vec(),
            alt_allele: alt_allele.to_vec(),
        };
        self.others
            .binary_search(&query)
            .ok()
            .map(|i| self.others[i].idx as usize)
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
                if let Some(blob) = self.json_blobs.as_ref().and_then(|b| b.get(idx)) {
                    if !blob.is_empty() {
                        // An empty alias marks a whole-record blob (see
                        // `writer_v2::raw_json_blob_fields`): the stored blob is
                        // the complete record object, so emit it verbatim rather
                        // than nesting it under a key. Such sources carry this
                        // as their sole field, so returning here is correct.
                        if field.alias.is_empty() {
                            return blob.to_string();
                        }
                        parts.push(format!("\"{}\":{}", field.alias, blob));
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

            let val_str =
                crate::fields::format_value(field, stored, strings.get(fi).map(|v| v.as_slice()));
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
                field: "AF".into(),
                alias: "af".into(),
                ftype: FieldType::Float,
                multiplier: 1_000_000,
                zigzag: false,
                missing_value: u32::MAX,
                missing_string: ".".into(),
                description: String::new(),
            },
            Field {
                field: "blob".into(),
                alias: "blob".into(),
                ftype: FieldType::JsonBlob,
                multiplier: 1,
                zigzag: false,
                missing_value: u32::MAX,
                missing_string: ".".into(),
                description: String::new(),
            },
            Field {
                field: "AC".into(),
                alias: "ac".into(),
                ftype: FieldType::Integer,
                multiplier: 1,
                zigzag: false,
                missing_value: u32::MAX,
                missing_string: ".".into(),
                description: String::new(),
            },
        ];

        let mut chunk = Chunk::empty();
        chunk.var32s = vec![var32::encode(100, b"A", b"G").unwrap()];
        // Two non-JsonBlob columns, in field order: AF then AC.
        chunk.values = vec![vec![1234], vec![42]];
        chunk.json_blobs = Some(JsonBlobs::from_rows([r#"{"k":1}"#]));

        let json = chunk.reconstruct_json(0, &fields, &[]);
        assert!(json.contains("\"af\":"), "missing af in: {}", json);
        assert!(json.contains("\"ac\":42"), "missing ac in: {}", json);
        assert!(
            json.contains("\"blob\":{\"k\":1}"),
            "missing blob in: {}",
            json
        );
    }

    #[test]
    fn test_chunk_reconstruct_whole_record_blob_verbatim() {
        // A single JsonBlob field with an empty alias means the stored blob is
        // the complete record object; reconstruct_json must return it verbatim
        // (not nested under a key), so v2 reproduces the v1 JSON exactly.
        let fields = vec![Field {
            field: String::new(),
            alias: String::new(),
            ftype: FieldType::JsonBlob,
            multiplier: 1,
            zigzag: false,
            missing_value: u32::MAX,
            missing_string: ".".into(),
            description: String::new(),
        }];
        let mut chunk = Chunk::empty();
        chunk.var32s = vec![var32::encode(100, b"A", b"G").unwrap()];
        let blob = r#"{"significance":["Pathogenic"],"reviewStatus":"criteria_provided"}"#;
        chunk.json_blobs = Some(JsonBlobs::from_rows([blob]));

        let json = chunk.reconstruct_json(0, &fields, &[]);
        assert_eq!(json, blob, "whole-record blob must round-trip verbatim");
    }

    #[test]
    fn test_chunk_reconstruct_json() {
        let fields = vec![
            Field {
                field: "AF".into(),
                alias: "allAf".into(),
                ftype: FieldType::Float,
                multiplier: 1_000_000,
                zigzag: false,
                missing_value: u32::MAX,
                missing_string: ".".into(),
                description: String::new(),
            },
            Field {
                field: "AC".into(),
                alias: "allAc".into(),
                ftype: FieldType::Integer,
                multiplier: 1,
                zigzag: false,
                missing_value: u32::MAX,
                missing_string: ".".into(),
                description: String::new(),
            },
        ];

        let mut chunk = Chunk::empty();
        chunk.var32s = vec![var32::encode(100, b"A", b"G").unwrap()];
        chunk.values = vec![
            vec![1234], // AF * 1_000_000 = 0.001234
            vec![42],   // AC = 42
        ];

        let json = chunk.reconstruct_json(0, &fields, &[]);
        assert!(json.contains("\"allAf\":"));
        assert!(json.contains("\"allAc\":42"));
    }
}
