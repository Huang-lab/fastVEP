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

/// Equivalent on-disk names for a chromosome.
///
/// Different upstream sources (and different user VCFs) mix `chr*` and bare
/// styles, and historically the SA writer canonicalized to the bare form
/// (`"1"`, `"X"`, `"MT"`) while modern inputs use `chr*`. This returns the
/// query name first, then all known aliases — so the reader can satisfy a
/// `chr1` query against an index built with `1`, or vice versa, without
/// requiring users to rebuild databases. See issue #37.
///
/// The returned `Vec` is small (1–3 entries) so caller-side iteration is
/// cheap. The first element is always the input name unchanged.
pub fn chrom_aliases(chrom: &str) -> Vec<String> {
    let mut out = Vec::with_capacity(3);
    out.push(chrom.to_string());

    // chr1 <-> 1
    if let Some(stripped) = chrom.strip_prefix("chr") {
        if !stripped.is_empty() && stripped != chrom {
            out.push(stripped.to_string());
        }
    } else {
        out.push(format!("chr{}", chrom));
    }

    // Mitochondrial special case: chrM / M / MT / chrMT all refer to the
    // same contig but UCSC uses chrM, NCBI uses MT.
    let mito_set = ["chrM", "M", "MT", "chrMT"];
    if mito_set.contains(&chrom) {
        for alt in mito_set {
            if alt != chrom && !out.iter().any(|n| n == alt) {
                out.push(alt.to_string());
            }
        }
    }

    out
}

#[cfg(test)]
mod chrom_alias_tests {
    use super::chrom_aliases;

    #[test]
    fn chr_prefix_round_trips() {
        let aliases = chrom_aliases("chr1");
        assert!(aliases.iter().any(|n| n == "chr1"));
        assert!(aliases.iter().any(|n| n == "1"));
    }

    #[test]
    fn bare_form_round_trips() {
        let aliases = chrom_aliases("1");
        assert!(aliases.iter().any(|n| n == "1"));
        assert!(aliases.iter().any(|n| n == "chr1"));
    }

    #[test]
    fn mitochondrial_aliases_cover_all_four_forms() {
        for name in ["chrM", "M", "MT", "chrMT"] {
            let aliases = chrom_aliases(name);
            for form in ["chrM", "M", "MT", "chrMT"] {
                assert!(
                    aliases.iter().any(|n| n == form),
                    "{} should resolve {}",
                    name,
                    form
                );
            }
        }
    }

    #[test]
    fn unknown_contig_returns_just_self_and_chr_variant() {
        let aliases = chrom_aliases("HLA-A*01:01");
        assert!(aliases.iter().any(|n| n == "HLA-A*01:01"));
    }
}

/// A single gene annotation record.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct GeneRecord {
    /// Gene symbol (e.g., "BRCA1").
    pub gene_symbol: String,
    /// Pre-serialized JSON annotation string.
    pub json: String,
}
