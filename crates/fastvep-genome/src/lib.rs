pub mod codon;
pub mod mitochondrial;
mod transcript;

pub use codon::CodonTable;
pub use mitochondrial::{
    is_mitochondrial, mitochondrial_codon_table, wrap_position, wrap_position_for, MT_LENGTH,
};
pub use transcript::{Exon, Gene, Transcript, Translation};
