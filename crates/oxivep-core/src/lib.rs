mod annotation_types;
mod consequence;
mod position;

pub use annotation_types::{GeneAnnotation, SupplementaryAnnotation};
pub use consequence::{Consequence, Impact};
pub use position::{Allele, GenomicPosition, Strand, VariantType};
