//! Supplementary annotation provider traits and types.
//!
//! These traits define the interface for plugging external annotation sources
//! (ClinVar, gnomAD, dbSNP, conservation scores, etc.) into the OxiVEP pipeline.
//!
//! The simple data types ([`SupplementaryAnnotation`], [`GeneAnnotation`]) live in
//! `oxivep-core` to avoid circular dependencies. They are re-exported here for
//! convenience.

use anyhow::Result;

// Re-export data types from core for convenience.
pub use oxivep_core::{GeneAnnotation, SupplementaryAnnotation};

/// The value returned by an annotation provider for a single query.
#[derive(Debug, Clone)]
pub enum AnnotationValue {
    /// A JSON string representing one or more allele-specific annotations.
    Json(String),
    /// A positional annotation (same value regardless of allele, e.g., PhyloP score).
    Positional(String),
    /// Interval-based annotations (e.g., overlapping SV regions).
    Interval(Vec<String>),
}

/// Metadata about a supplementary annotation data source.
#[derive(Debug, Clone)]
pub struct SaMetadata {
    /// Human-readable name of the data source (e.g., "ClinVar").
    pub name: String,
    /// Version string (e.g., "2024-12-01").
    pub version: String,
    /// Release date or description.
    pub description: String,
    /// Genome assembly this data source is built for (e.g., "GRCh38").
    pub assembly: String,
    /// The JSON key used in output (e.g., "clinvar", "gnomad").
    pub json_key: String,
    /// Whether annotations are matched by allele (true) or positional (false).
    pub match_by_allele: bool,
    /// Whether the output is an array of annotations (true) or a single object (false).
    pub is_array: bool,
    /// Whether this is a positional annotation (same for all alleles at a position).
    pub is_positional: bool,
}

/// Trait for providing supplementary annotations at the variant level.
///
/// Implementations must be `Send + Sync` to support parallel annotation via rayon.
/// Each provider handles one data source (e.g., ClinVar, gnomAD, PhyloP).
pub trait AnnotationProvider: Send + Sync {
    /// Short name of this provider (e.g., "ClinVar").
    fn name(&self) -> &str;

    /// The JSON key used in structured output (e.g., "clinvar").
    fn json_key(&self) -> &str;

    /// Metadata about this annotation source.
    fn metadata(&self) -> &SaMetadata;

    /// Look up annotations for a specific variant position and alleles.
    ///
    /// Returns `None` if no annotation exists for this position/allele combination.
    /// For allele-specific sources (`match_by_allele = true`), `ref_allele` and
    /// `alt_allele` are used to match. For positional sources, only `chrom` and `pos`
    /// are used.
    fn annotate_position(
        &self,
        chrom: &str,
        pos: u64,
        ref_allele: &str,
        alt_allele: &str,
    ) -> Result<Option<AnnotationValue>>;

    /// Pre-load data for a batch of positions on a chromosome.
    ///
    /// Called before the parallel annotation phase to decompress and cache
    /// relevant blocks. The default implementation is a no-op.
    fn preload(&self, _chrom: &str, _positions: &[u64]) -> Result<()> {
        Ok(())
    }
}

/// Trait for providing gene-level annotations (OMIM, pLI scores, etc.).
///
/// Gene annotations are keyed by gene symbol rather than genomic position.
pub trait GeneAnnotationProvider: Send + Sync {
    /// Short name of this provider (e.g., "OMIM").
    fn name(&self) -> &str;

    /// The JSON key used in structured output (e.g., "omim").
    fn json_key(&self) -> &str;

    /// Look up annotations for a gene by its symbol (e.g., "BRCA1").
    ///
    /// Returns a JSON string if annotations exist, or `None`.
    fn annotate_gene(&self, gene_symbol: &str) -> Result<Option<String>>;
}
