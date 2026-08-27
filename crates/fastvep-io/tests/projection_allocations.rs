//! Pins the allocation behaviour of the VCF projection path.
//!
//! Every loaded supplementary source is offered every annotation, so
//! `format_supplementary_vcf_info` runs its per-source helpers once per
//! (variant x transcript x allele x source) whether or not that allele
//! carries a payload for that source. Those helpers used to render the
//! uploaded allele, split the VCF ALT column into a `Vec`, and build a
//! `HashSet` *before* checking whether any payload matched - so a run with no
//! `--sa-dir` at all still paid for all of it on every annotation, and threw
//! the results away. Each helper now establishes that a payload exists before
//! allocating anything.
//!
//! This file installs a counting global allocator, so it deliberately holds
//! exactly ONE test: `cargo test` runs the tests in a binary concurrently, and
//! a second test allocating on another thread would be counted here.

use fastvep_core::{Allele, Consequence, Impact, Strand, VariantType};
use fastvep_io::output::format_supplementary_vcf_info;
use fastvep_io::variant::{AlleleAnnotation, TranscriptVariation, VariationFeature};
use std::alloc::{GlobalAlloc, Layout, System};
use std::sync::atomic::{AtomicUsize, Ordering};

static ALLOCS: AtomicUsize = AtomicUsize::new(0);

struct Counting;

unsafe impl GlobalAlloc for Counting {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        ALLOCS.fetch_add(1, Ordering::Relaxed);
        unsafe { System.alloc(layout) }
    }
    unsafe fn dealloc(&self, ptr: *mut u8, layout: Layout) {
        unsafe { System.dealloc(ptr, layout) }
    }
    unsafe fn realloc(&self, ptr: *mut u8, layout: Layout, new_size: usize) -> *mut u8 {
        ALLOCS.fetch_add(1, Ordering::Relaxed);
        unsafe { System.realloc(ptr, layout, new_size) }
    }
}

#[global_allocator]
static ALLOCATOR: Counting = Counting;

fn allocations_during<F: FnOnce()>(f: F) -> usize {
    let before = ALLOCS.load(Ordering::Relaxed);
    f();
    ALLOCS.load(Ordering::Relaxed) - before
}

fn annotation(allele: &str) -> AlleleAnnotation {
    AlleleAnnotation {
        allele: Allele::Sequence(allele.as_bytes().to_vec()),
        consequences: vec![Consequence::MissenseVariant],
        impact: Impact::Moderate,
        cdna_position: Some((100, 100)),
        cds_position: Some((90, 90)),
        protein_position: Some((30, 30)),
        amino_acids: None,
        codons: None,
        exon: Some((2, 6)),
        intron: None,
        distance: None,
        protein_length: Some(300),
        escapes_nmd: None,
        hgvsc: None,
        hgvsp: None,
        hgvsg: None,
        hgvs_offset: None,
        existing_variation: Vec::new(),
        sift: None,
        polyphen: None,
        // The point of the test: no supplementary payload for any source.
        supplementary: Vec::new(),
        acmg_classification: None,
    }
}

fn transcript_variation(id: &str, alts: &[&str]) -> TranscriptVariation {
    TranscriptVariation {
        transcript_id: id.into(),
        gene_id: "ENSG00000000001".into(),
        gene_symbol: Some("GENE1".into()),
        biotype: "protein_coding".into(),
        allele_annotations: alts.iter().map(|a| annotation(a)).collect(),
        canonical: true,
        strand: Strand::Forward,
        source: None,
        protein_id: None,
        mane_select: None,
        mane_plus_clinical: None,
        tsl: None,
        appris: None,
        ccds: None,
        gencode_primary: false,
        symbol_source: None,
        hgnc_id: None,
        flags: Vec::new(),
    }
}

#[test]
fn projecting_a_variant_with_no_supplementary_payload_allocates_nothing() {
    // Twelve transcripts, two alleles each: the shape of a variant in a
    // gene-dense window, where the per-source helpers run most often.
    let alts = ["C", "G"];
    let vf = VariationFeature {
        position: fastvep_core::GenomicPosition::new("chr1", 1000, 1000, Strand::Forward),
        allele_string: "A/C/G".to_string(),
        ref_allele: Allele::Sequence(b"A".to_vec()),
        alt_alleles: alts
            .iter()
            .map(|a| Allele::Sequence(a.as_bytes().to_vec()))
            .collect(),
        variation_name: None,
        vcf_fields: None,
        transcript_variations: (0..12)
            .map(|i| transcript_variation(&format!("ENST{i:011}"), &alts))
            .collect(),
        existing_variants: Vec::new(),
        minimised: false,
        most_severe_consequence: Some(Consequence::MissenseVariant),
        variant_type: VariantType::Snv,
        sv_end: None,
        sv_len: None,
        supplementary_annotations: Vec::new(),
        gene_annotations: Vec::new(),
    };

    // Warm any lazily-initialised state so the measured call sees none of it.
    let _ = format_supplementary_vcf_info(&vf);

    let mut projected = Vec::new();
    let allocations = allocations_during(|| {
        projected = format_supplementary_vcf_info(&vf);
    });

    assert!(projected.is_empty(), "nothing to project: {projected:?}");
    assert_eq!(
        allocations, 0,
        "projecting a variant with no supplementary payload should not allocate; \
         each per-source helper must confirm a matching payload before it builds \
         the uploaded allele, the dedupe set, or the value vector"
    );
}
