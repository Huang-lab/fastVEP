//! Simple data types for supplementary and gene-level annotations.
//!
//! These live in oxivep-core so that both oxivep-io and oxivep-cache can use
//! them without creating circular dependencies.

/// A collected supplementary annotation ready for output.
#[derive(Debug, Clone)]
pub struct SupplementaryAnnotation {
    /// The JSON key (e.g., "clinvar", "gnomad").
    pub json_key: String,
    /// Whether this is an array annotation.
    pub is_array: bool,
    /// The pre-serialized JSON string.
    pub json_string: String,
}

/// A collected gene-level annotation ready for output.
#[derive(Debug, Clone)]
pub struct GeneAnnotation {
    /// The gene symbol this annotation applies to.
    pub gene_symbol: String,
    /// The JSON key (e.g., "omim").
    pub json_key: String,
    /// The pre-serialized JSON string.
    pub json_string: String,
}
