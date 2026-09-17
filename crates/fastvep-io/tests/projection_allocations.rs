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
//! The same bar applies to the custom-source projection added for #116: it is
//! resolved per run from the loaded `json_key`s, not from a static descriptor,
//! but it runs in the same per-allele loop and so confirms a payload before it
//! allocates. The measured run has one loaded.
//!
//! The second property here is the *payload* path rather than the empty one.
//! `FV_CLINVAR_PROTEIN` renders a `&`-joined list of `pos:ref>alt:sig` records,
//! and ClinVar carries tens of them for a well-studied gene - 31 for FHL1 in
//! the report behind #123. Composing each record out of escaped `String`s and
//! a `format!` cost nine allocations per record, once per (variant x
//! transcript) for the tab writer; the records are written into one buffer
//! instead, so the count is bounded by that buffer's growth and does not scale
//! with the number of records.
//!
//! This file installs a counting global allocator, so it deliberately holds
//! exactly ONE test: `cargo test` runs the tests in a binary concurrently, and
//! a second test allocating on another thread would be counted here. Both
//! properties are therefore asserted in that one test, and its name says so -
//! a payload-budget regression failing under a name about the empty case
//! sends the next reader to the wrong function.

use fastvep_core::{Allele, Consequence, GeneAnnotation, Impact, Strand, VariantType};
use fastvep_io::output::{format_supplementary_vcf_info, LoadedSupplementarySpecs};
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
        intron_offset: None,
        shifted_intron_offset: None,
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
fn projection_allocations_stay_within_budget_with_and_without_a_payload() {
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

    // A run with a user's own `custom_vcf` database loaded. Its projection is
    // resolved from the source's `json_key` rather than a static descriptor
    // (#116), so it is the one helper that could not exist before the budget
    // was written - and it has to clear the same bar: no payload for this
    // allele, no allocation. Built outside the measured window because
    // resolving the loaded set is a once-per-run cost by design.
    let specs = LoadedSupplementarySpecs::new(&["my_panel".to_string()], &[]);

    // Warm any lazily-initialised state so the measured call sees none of it.
    let _ = format_supplementary_vcf_info(&vf, &specs);

    let mut projected = Vec::new();
    let allocations = allocations_during(|| {
        projected = format_supplementary_vcf_info(&vf, &specs);
    });

    assert!(projected.is_empty(), "nothing to project: {projected:?}");
    assert_eq!(
        allocations, 0,
        "projecting a variant with no supplementary payload should not allocate; \
         each per-source helper must confirm a matching payload before it builds \
         the uploaded allele, the dedupe set, or the value vector"
    );

    // Same variant, now carrying a gene-level ClinVar-protein payload with 40
    // records - the shape of a well-studied gene.
    let records: Vec<String> = (1..=40)
        .map(|i| format!(r#"{{"pos":{i},"refAa":"R","altAa":"H","sig":"Pathogenic"}}"#))
        .collect();
    let mut with_payload = vf;
    with_payload.gene_annotations = vec![GeneAnnotation {
        gene_symbol: "FHL1".to_string(),
        json_key: "clinvar_protein".to_string(),
        json_string: format!(r#"{{"proteinVariants":[{}]}}"#, records.join(",")),
    }];
    let gene_specs = LoadedSupplementarySpecs::new(&[], &["clinvar_protein".to_string()]);

    let mut rendered = Vec::new();
    let _ = format_supplementary_vcf_info(&with_payload, &gene_specs);
    let allocations = allocations_during(|| {
        rendered = format_supplementary_vcf_info(&with_payload, &gene_specs);
    });

    let value = &rendered
        .iter()
        .find(|(id, _)| id == "FV_CLINVAR_PROTEIN")
        .expect("the payload projects")
        .1;
    assert_eq!(
        value.matches('&').count(),
        39,
        "all 40 records render: {value}"
    );
    assert!(
        value.contains("1:R>H:Pathogenic"),
        "with literal delimiters: {value}"
    );

    // Measured on this payload: 663 when each record was composed with a
    // `format!` and escaped whole, 743 when the four leaves were escaped into
    // four `String`s to keep the delimiters literal (#123), and 345 once the
    // records were written into one buffer instead.
    //
    // What remains per record is `serde_json` parsing the payload into a
    // `Value` - a `Map` and a `String` per key - which this path has always
    // paid and which the budget does not try to hide. The bound is set below
    // every per-record-`String` spelling, so a return to one fails here rather
    // than in a profile six months later.
    eprintln!("ClinVar-protein payload projection: {allocations} allocations");
    assert!(
        allocations <= 400,
        "rendering 40 ClinVar-protein records took {allocations} allocations, \
         over the 400 budget; the records must be written into one buffer \
         rather than composed out of per-record `String`s"
    );
}
