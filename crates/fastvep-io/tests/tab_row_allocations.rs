//! Pins what the tab writer costs, in heap allocations, for each row it writes.
//!
//! The row loop runs once per (variant x transcript x allele): a variant in a gene-dense window has
//! dozens of rows and a callset has millions of variants, so a cell that is a `String` of its own,
//! collected into a `Vec` and joined, is built millions of times. The writer builds each row in
//! the one `String` that is then handed back, borrows every cell that already exists somewhere
//! (the transcript's names, the amino acids, the HGVS strings), and writes the ones it has to work
//! out (a position, a distance, a strand) straight into that `String`.
//!
//! So a row costs one allocation, its own `String`, and the rest of what is held here is the
//! per-variant constant: the line's `Vec`, `Uploaded_variation` and `Location`, which do not grow
//! with the number of rows.
//!
//! This file installs a counting global allocator that counts the allocations of the thread that is
//! measured and of no other: the test harness has a thread of its own that allocates while a test
//! runs, and a count of the whole process took those into the window, which failed a test of this
//! kind on a CI runner.

use fastvep_core::{
    Allele, Consequence, Impact, Strand, TranscriptLocation, TranscriptSpan, VariantType,
};
use fastvep_io::output::{format_tab_line_with, ColumnFlags, LoadedSupplementarySpecs, TabOptions};
use fastvep_io::variant::{AlleleAnnotation, TranscriptVariation, VariationFeature};
use std::alloc::{GlobalAlloc, Layout, System};
use std::cell::Cell;

thread_local! {
    // The allocations of the calling thread only. The test harness has a thread of its own that allocates
    // while a test runs, and a count of the whole process took those into the window that is measured.
    static ALLOCS: Cell<usize> = const { Cell::new(0) };
}

struct Counting;

unsafe impl GlobalAlloc for Counting {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        ALLOCS.with(|count| count.set(count.get() + 1));
        unsafe { System.alloc(layout) }
    }
    unsafe fn dealloc(&self, ptr: *mut u8, layout: Layout) {
        unsafe { System.dealloc(ptr, layout) }
    }
    unsafe fn realloc(&self, ptr: *mut u8, layout: Layout, new_size: usize) -> *mut u8 {
        ALLOCS.with(|count| count.set(count.get() + 1));
        unsafe { System.realloc(ptr, layout, new_size) }
    }
}

#[global_allocator]
static ALLOCATOR: Counting = Counting;

fn allocations_during<F: FnOnce()>(f: F) -> usize {
    let before = ALLOCS.with(Cell::get);
    f();
    ALLOCS.with(Cell::get) - before
}

/// A missense row as the annotation pass leaves it: positions, residues, codons and the three HGVS
/// strings all there to be written.
fn annotation(allele: &str) -> AlleleAnnotation {
    AlleleAnnotation {
        allele: Allele::Sequence(allele.as_bytes().to_vec()),
        consequences: vec![
            Consequence::MissenseVariant,
            Consequence::SpliceRegionVariant,
        ],
        impact: Impact::Moderate,
        location: TranscriptLocation {
            cdna: TranscriptSpan::new(Some(100), Some(100)),
            cds: TranscriptSpan::new(Some(90), Some(90)),
            protein: TranscriptSpan::new(Some(30), Some(30)),
            ..Default::default()
        },
        protein_position: Some((30, 30)),
        amino_acids: Some(("R".into(), "Q".into())),
        codons: Some(("cGg".into(), "cAg".into())),
        exon: Some((2, 6)),
        splice: Default::default(),
        distance: None,
        protein_length: Some(300),
        escapes_nmd: None,
        hgvsc: Some("ENST00000000001.1:c.90G>A".into()),
        hgvsp: Some("ENSP00000000001.1:p.Arg30Gln".into()),
        hgvsg: Some("chr1:g.1000A>C".into()),
        hgvs_offset: Some(2),
        existing_variation: Vec::new(),
        sift: None,
        polyphen: None,
        supplementary: Vec::new(),
        acmg_classification: None,
        pick: false,
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
        symbol_source: Some("HGNC".into()),
        hgnc_id: Some("HGNC:1".into()),
        flags: vec!["cds_start_NF".into()],
    }
}

fn variation(transcripts: usize, alts: &[&str]) -> VariationFeature {
    VariationFeature {
        position: fastvep_core::GenomicPosition::new("chr1", 1000, 1000, Strand::Forward),
        allele_string: format!("A/{}", alts.join("/")),
        ref_allele: Allele::Sequence(b"A".to_vec()),
        alt_alleles: alts
            .iter()
            .map(|a| Allele::Sequence(a.as_bytes().to_vec()))
            .collect(),
        variation_name: None,
        vcf_fields: None,
        transcript_variations: (0..transcripts)
            .map(|i| transcript_variation(&format!("ENST{i:011}"), alts))
            .collect(),
        known: Vec::new(),
        minimised: false,
        most_severe_consequence: Some(Consequence::MissenseVariant),
        variant_type: VariantType::Snv,
        sv_end: None,
        sv_len: None,
        supplementary_annotations: Vec::new(),
        gene_annotations: Vec::new(),
    }
}

#[test]
fn a_tab_row_costs_one_allocation_and_the_variant_a_few() {
    let flags = ColumnFlags {
        show_ref_allele: true,
        variant_class: true,
        symbol: true,
        biotype: true,
        canonical: true,
        numbers: true,
        hgvsg: true,
        ..Default::default()
    };
    // `--hgvs` adds HGVSc, HGVSp and HGVS_OFFSET to the optional columns
    let optional = flags.columns(true);
    let specs = LoadedSupplementarySpecs::new(&[], &[]);
    let opts = TabOptions {
        optional: &optional,
        ..Default::default()
    };

    let alts = ["C", "G"];
    let (small, large) = (variation(1, &alts), variation(30, &alts));
    // warm anything lazily initialised
    let _ = format_tab_line_with(&small, &specs, false, opts);

    let mut lines = Vec::new();
    let one = allocations_during(|| lines = format_tab_line_with(&small, &specs, false, opts));
    assert_eq!(lines.len(), 2);
    let many = allocations_during(|| lines = format_tab_line_with(&large, &specs, false, opts));
    assert_eq!(lines.len(), 60);
    // every cell is there: the 17 base columns, REF_ALLELE and the nine optional ones
    for line in &lines {
        assert_eq!(line.split('\t').count(), 17 + optional.len(), "{line}");
    }
    assert!(lines[0].contains("p.Arg30Gln\t2\t"), "{}", lines[0]);

    eprintln!(
        "tab writer: {one} allocations for 2 rows, {many} for 60 rows ({:.2} a row)",
        many as f64 / 60.0
    );
    // The rows: one `String` each, and the few of the variant (the `Vec` of lines, `Uploaded_variation`
    // and `Location`). A cell that is built into a `String` of its own is 60 more and fails here, and
    // so does a per-row `String` that is not the row's.
    assert!(
        many <= 60 + 8,
        "60 rows cost {many} allocations, more than a `String` each and a few for the variant: a \
         cell is being built into a `String` of its own"
    );
    assert!(
        many.saturating_sub(one) <= 58 + 2,
        "58 more rows cost {} more allocations: the per-row cost is more than its own `String`",
        many.saturating_sub(one)
    );
}
