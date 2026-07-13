//! Constants and common types for the fastSA binary annotation format.

/// Magic bytes for position/allele-level annotation files (.osa).
pub const OSA_MAGIC: &[u8; 8] = b"FSTSA_01";

/// Magic bytes for interval-level annotation files (.osi).
pub const OSI_MAGIC: &[u8; 8] = b"FSTSI_01";

/// Magic bytes for gene-level annotation files (.oga).
pub const OGA_MAGIC: &[u8; 8] = b"FSTGA_01";

/// Current schema version. Bump when the binary format changes.
pub const SCHEMA_VERSION: u16 = 1;

/// Default block size for compression (8 MiB).
pub const DEFAULT_BLOCK_SIZE: usize = 8 * 1024 * 1024;

/// Default zstd compression level (3 is a good speed/ratio tradeoff).
pub const ZSTD_LEVEL: i32 = 3;

/// Hard cap on a single bincode-serialized index payload (4 GiB). Used by
/// `.osa.idx`, `.osi`, and `.oga` readers to refuse malformed/malicious files
/// that claim absurd payload sizes, before allocating a buffer.
///
/// Stored as `u64` so the literal compiles on 32-bit targets (where
/// `usize` is 32 bits and cannot hold `2^32`). Readers must additionally
/// verify the value fits in `usize` before allocating.
pub const MAX_INDEX_PAYLOAD: u64 = 4 * 1024 * 1024 * 1024;

/// File extension for position/allele-level annotations.
pub const OSA_EXT: &str = "osa";

/// File extension for the index file.
pub const IDX_EXT: &str = "osa.idx";

/// File extension for interval-level annotations.
pub const OSI_EXT: &str = "osi";

/// File extension for gene-level annotations.
pub const OGA_EXT: &str = "oga";

/// A single annotation record ready for writing.
#[derive(Debug, Clone)]
pub struct AnnotationRecord {
    /// Chromosome index (numeric, mapped externally).
    pub chrom_idx: u16,
    /// 1-based genomic position.
    pub position: u32,
    /// Reference allele (empty string for positional annotations).
    pub ref_allele: String,
    /// Alternate allele (empty string for positional annotations).
    pub alt_allele: String,
    /// Pre-serialized JSON annotation string.
    pub json: String,
}

/// A single interval annotation record.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct IntervalRecord {
    /// Chromosome name.
    pub chrom: String,
    /// 1-based start position (inclusive).
    pub start: u32,
    /// 1-based end position (inclusive).
    pub end: u32,
    /// Pre-serialized JSON annotation string.
    pub json: String,
}

/// Chromosome name to index mapping for efficient lookups.
#[derive(Debug, Clone)]
pub struct ChromMap {
    name_to_idx: std::collections::HashMap<String, u16>,
}

impl ChromMap {
    /// Create a standard human chromosome mapping (chr1-22, chrX, chrY, chrM).
    pub fn standard_human() -> Self {
        let mut map = std::collections::HashMap::new();
        for i in 1..=22u16 {
            map.insert(format!("chr{}", i), i - 1);
            map.insert(format!("{}", i), i - 1);
        }
        map.insert("chrX".into(), 22);
        map.insert("X".into(), 22);
        map.insert("chrY".into(), 23);
        map.insert("Y".into(), 23);
        map.insert("chrM".into(), 24);
        map.insert("MT".into(), 24);
        map.insert("chrMT".into(), 24);
        Self { name_to_idx: map }
    }

    /// Look up a chromosome index by name.
    #[inline]
    pub fn get(&self, name: &str) -> Option<u16> {
        self.name_to_idx.get(name).copied()
    }
}

/// Equivalent on-disk names for a chromosome (chr↔bare, mitochondrial forms).
///
/// Re-exported from `fastvep-core` so SA readers and the cache builder share a
/// single implementation. See [`fastvep_core::chrom::chrom_aliases`].
pub use fastvep_core::chrom_aliases;

/// A single gene annotation record.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct GeneRecord {
    /// Gene symbol (e.g., "BRCA1").
    pub gene_symbol: String,
    /// Pre-serialized JSON annotation string.
    pub json: String,
}

/// Escape a string for embedding in a hand-built JSON string value.
///
/// Source parsers build their `.osa`/`.oga` JSON payloads by hand (not via
/// `serde_json`) for speed, so any field taken from an upstream file
/// (gene names, disease descriptions, COSMIC/ClinVar free-text fields,
/// etc.) must run through this before being interpolated into a `"..."`
/// value — otherwise a raw `"`, `\`, or control character in the source
/// data produces invalid JSON that silently corrupts that record for
/// every downstream `serde_json::from_str` consumer.
pub fn escape_json(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            // JSON requires every C0 control character (U+0000..=U+001F) to
            // be escaped; the common ones have short forms above, the rest
            // (backspace, form feed, vertical tab, NUL, etc.) must go out as
            // \u00XX or the record is invalid JSON.
            c if (c as u32) < 0x20 => {
                out.push_str(&format!("\\u{:04x}", c as u32));
            }
            c => out.push(c),
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn escape_json_produces_valid_json_for_quotes_backslashes_and_controls() {
        // A field carrying a quote, a backslash, and assorted control
        // characters (including ones without a short escape form: backspace,
        // form feed, vertical tab, NUL) must still round-trip as valid JSON.
        let raw = "a\"b\\c\n\r\t\x08\x0c\x0b\x00d";
        let json = format!(r#"{{"v":"{}"}}"#, escape_json(raw));
        let parsed: serde_json::Value =
            serde_json::from_str(&json).expect("escaped string must produce valid JSON");
        assert_eq!(parsed["v"], raw, "value must survive an escape/parse round-trip unchanged");
    }
}
