mod annotation_types;
pub mod chrom;
mod consequence;
mod position;

pub use annotation_types::{GeneAnnotation, SupplementaryAnnotation};
pub use chrom::{chrom_aliases, looks_like_refseq_accession, ChromSynonyms};
pub use consequence::{Consequence, Impact};
pub use position::{Allele, GenomicPosition, Strand, VariantType};
