//! Pins the cost of VEP's known-variant columns in the row writers.
//!
//! The cells (`AF`, `AFR_AF`, `gnomADg_NFE_AF`, `MAX_AF`, `CLIN_SIG`, ...) are
//! worked out once for each alt allele of a variant, and every transcript row of
//! that allele writes the same text. The row writers run once per (variant x
//! transcript x allele), and 33 columns that each built a `String` for the cell
//! would build 33 of them on every row: a variant in a gene-dense window has dozens
//! of rows, a callset has millions of variants.
//!
//! So the property held here is not a count but a difference: a variant that
//! carries all 33 cells allocates exactly as much in the tab writer as the same
//! variant carrying none (the cells are borrowed, and the `-` of an empty one is a
//! constant), and in the CSQ writer no more than the growth of the buffer that
//! the longer text needs.
//!
//! This file installs a counting global allocator that counts the allocations of the thread that is
//! measured and of no other: the test harness has a thread of its own that allocates while a test
//! runs, and a count of the whole process took those into the window, which failed a test of this
//! kind on a CI runner.
//! Both writers are asserted in that one test.

use fastvep_core::{
    Allele, Consequence, Impact, Strand, TranscriptLocation, TranscriptSpan, VariantType,
};
use fastvep_io::known::{AlleleKnown, KnownCells, KnownColumn, KnownFlags};
use fastvep_io::output::{
    csq_fields, format_csq_with, format_tab_line_with, ColumnFlags, CsqOptions,
    LoadedSupplementarySpecs, TabOptions,
};
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

fn annotation(allele: &str) -> AlleleAnnotation {
    AlleleAnnotation {
        allele: Allele::Sequence(allele.as_bytes().to_vec()),
        consequences: vec![Consequence::MissenseVariant],
        impact: Impact::Moderate,
        location: TranscriptLocation {
            cdna: TranscriptSpan::new(Some(100), Some(100)),
            cds: TranscriptSpan::new(Some(90), Some(90)),
            protein: TranscriptSpan::new(Some(30), Some(30)),
            ..Default::default()
        },
        protein_position: Some((30, 30)),
        amino_acids: None,
        codons: None,
        exon: Some((2, 6)),
        splice: Default::default(),
        distance: None,
        protein_length: Some(300),
        escapes_nmd: None,
        hgvsc: None,
        hgvsp: None,
        hgvsg: None,
        hgvs_offset: None,
        existing_variation: vec!["rs1".to_string()],
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
        symbol_source: None,
        hgnc_id: None,
        flags: Vec::new(),
    }
}

/// A cell for every one of the 33 columns.
fn full_cells() -> KnownCells {
    let mut cells = KnownCells::default();
    for (i, column) in KnownColumn::ALL.into_iter().enumerate() {
        cells.set(column, format!("0.{:04},{}", i + 1, i % 7));
    }
    cells
}

#[test]
fn known_columns_add_no_allocation_per_row_to_either_writer() {
    // Thirty transcripts, two alleles each: sixty rows of a variant in a
    // gene-dense window.
    let alts = ["C", "G"];
    let mut vf = VariationFeature {
        position: fastvep_core::GenomicPosition::new("chr1", 1000, 1000, Strand::Forward),
        allele_string: "A/C/G".to_string(),
        ref_allele: Allele::Sequence(b"A".to_vec()),
        alt_alleles: alts
            .iter()
            .map(|a| Allele::Sequence(a.as_bytes().to_vec()))
            .collect(),
        variation_name: None,
        vcf_fields: None,
        transcript_variations: (0..30)
            .map(|i| transcript_variation(&format!("ENST{i:011}"), &alts))
            .collect(),
        known: Vec::new(),
        minimised: false,
        most_severe_consequence: Some(Consequence::MissenseVariant),
        variant_type: VariantType::Snv,
        sv_end: None,
        sv_len: None,
        supplementary_annotations: Vec::new(),
        gene_annotations: Vec::new(),
    };

    let flags = ColumnFlags {
        known: KnownFlags::everything(),
        ..Default::default()
    };
    let optional = flags.columns(false);
    assert_eq!(optional.len(), 33, "every known column is asked for");
    let fields = csq_fields(false, flags);
    let specs = LoadedSupplementarySpecs::new(&[], &[]);
    let opts = TabOptions {
        optional: &optional,
        ..Default::default()
    };

    let measure = |vf: &VariationFeature| -> (usize, usize, Vec<String>, String) {
        // warm anything lazily initialised
        let _ = format_tab_line_with(vf, &specs, false, opts);
        let _ = format_csq_with(vf, &fields, CsqOptions::default());
        let mut lines = Vec::new();
        let tab = allocations_during(|| {
            lines = format_tab_line_with(vf, &specs, false, opts);
        });
        let mut csq = String::new();
        let csq_allocations = allocations_during(|| {
            csq = format_csq_with(vf, &fields, CsqOptions::default());
        });
        (tab, csq_allocations, lines, csq)
    };

    let (tab_without, csq_without, lines_without, csq_text_without) = measure(&vf);
    assert_eq!(lines_without.len(), 60);
    assert!(
        lines_without
            .iter()
            .all(|l| l.split('\t').filter(|c| *c == "-").count() >= 33),
        "without cells the 33 columns are `-`"
    );

    vf.known = alts
        .iter()
        .map(|a| AlleleKnown {
            allele: Allele::Sequence(a.as_bytes().to_vec()),
            cells: full_cells(),
        })
        .collect();
    let (tab_with, csq_with, lines_with, csq_text_with) = measure(&vf);
    assert_eq!(lines_with.len(), 60);
    assert!(
        lines_with.iter().all(|l| l.contains("0.0001,0\t")),
        "the cells reach every row"
    );

    eprintln!(
        "tab writer, 60 rows: {tab_without} allocations without known cells, {tab_with} with 33; \
         CSQ writer: {csq_without} and {csq_with}"
    );
    assert_eq!(
        tab_with,
        tab_without,
        "33 known-variant cells on each of 60 rows cost {} allocations in the tab writer; a cell \
         is borrowed from the variant and must not be built into a `String` for each row",
        tab_with as i64 - tab_without as i64
    );
    assert!(
        csq_with <= csq_without + 3,
        "the CSQ writer took {csq_with} allocations with 33 cells on each of 60 entries against \
         {csq_without} without: more than the growth of its buffer, so a cell is being built \
         into a `String` per entry"
    );
    assert!(csq_text_with.len() > csq_text_without.len() + 60 * 33);
}
