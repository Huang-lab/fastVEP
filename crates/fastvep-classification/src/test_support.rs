//! Shared constructors for classification unit tests.
//!
//! [`ClassificationInput`] has grown to well over thirty fields, and every test
//! module used to build it with a full struct literal. Adding one field to the
//! struct then broke every one of those literals at compile time, which made
//! the cost of extending the input proportional to the number of tests instead
//! of to the change itself.
//!
//! Tests now start from [`minimal_input`] and override only the fields they are
//! actually exercising, so a new field is added in exactly one place. The
//! defaults here are deliberately inert: no annotations, no consequences, no
//! genotypes, so nothing fires unless a test asks for it.
//!
//! Every source of the reference stack is open and none has a record for the
//! variant: a complete database directory and a variant it does not know, which is
//! the situation the criteria were written for. A test of what a missing source does
//! overrides `sources`.

use crate::sa_extract::ClassificationInput;
use crate::sources::LoadedSources;
use fastvep_core::Impact;

/// A `ClassificationInput` with every annotation absent.
///
/// `gene_symbol` is `Some("TEST")` and `is_canonical` is true because criteria
/// that key off a named canonical transcript would otherwise take their
/// "nothing to evaluate" path in every test.
pub(crate) fn minimal_input() -> ClassificationInput {
    ClassificationInput {
        consequences: vec![],
        impact: Impact::Modifier,
        gene_symbol: Some("TEST".to_string()),
        is_canonical: true,
        amino_acids: None,
        protein_position: None,
        gnomad: None,
        clinvar: None,
        revel: None,
        splice_ai: None,
        dbnsfp: None,
        phylop: None,
        gerp: None,
        gene_constraints: None,
        omim: None,
        clinvar_protein: None,
        hgvs_c: None,
        predicted_nmd: None,
        nmd_escape_50nt: None,
        protein_truncation_pct: None,
        is_last_exon: None,
        in_critical_region: None,
        alt_start_codon_distance: None,
        same_splice_position_pathogenic: None,
        variant_coordinates: None,
        in_repeat_region: None,
        sources: LoadedSources::all(),
        is_pure_insertion: None,
        at_exon_edge: None,
        intronic_offset: None,
        hgvs_placement: None,
        proband_genotype: None,
        mother_genotype: None,
        father_genotype: None,
        companion_variants: vec![],
        functional_evidence: None,
    }
}
