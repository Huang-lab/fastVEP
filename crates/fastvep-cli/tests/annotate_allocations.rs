//! Pins what `fastvep annotate` costs, in heap allocations, for each record and each row it writes.
//!
//! The per-variant loop runs once per (variant x transcript x allele), and about a quarter of the
//! instructions of a run over a gene-dense locus are in `malloc` and `free`: a `String` built in
//! that loop is built millions of times. A benchmark in a commit message does not keep it from
//! coming back, so the property is held here as a count.
//!
//! The count is a difference. A run pays for its start (the GFF3 or its cache, the FASTA, the
//! rayon pool) whatever the input is, and that is not what is held. The same world is annotated
//! twice, with 500 records and with 2,500 (the fixture's records over again, so that every one costs
//! what it costs), and what the 2,000 more records and the rows they wrote cost is divided out.
//! Both runs read the transcript cache the first one wrote, so neither has the GFF3 to parse.
//!
//! The budgets are what the loop allocates for each record and each row today, and a little over.
//! A change that adds a `String` to the loop adds a row's worth of them to the count and fails it;
//! one that removes some lowers the number it prints with `--nocapture`, and the budget comes down
//! with it.
//!
//! This file installs a counting global allocator, so it deliberately holds exactly ONE test:
//! `cargo test` runs the tests of a binary concurrently, and a second one allocating on another
//! thread would be counted here. The rayon pool is one thread, so that the count does not depend on
//! how the work was split.

use std::alloc::{GlobalAlloc, Layout, System};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};

use fastvep_cli::pipeline::{run_annotate, AnnotateConfig, ColumnFlags, PickFlags};
use tempfile::TempDir;

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

fn world(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("fixtures/vep_tab_world")
        .join(name)
}

/// The fixture's records, `copies` times over, under its own header.
fn vcf(copies: usize) -> String {
    let text = std::fs::read_to_string(world("in.vcf")).unwrap();
    let (header, records): (Vec<&str>, Vec<&str>) = text.lines().partition(|l| l.starts_with('#'));
    let mut out = header.join("\n");
    out.push('\n');
    for _ in 0..copies {
        out.push_str(&records.join("\n"));
        out.push('\n');
    }
    out
}

struct Run {
    allocations: usize,
    records: usize,
    rows: usize,
}

/// One run over `copies` copies of the fixture's records, with the transcript cache of `dir`.
fn run(dir: &Path, copies: usize, columns: ColumnFlags, hgvs: bool) -> Run {
    let input = dir.join(format!("in_{copies}.vcf"));
    std::fs::write(&input, vcf(copies)).unwrap();
    let out = dir.join(format!("out_{copies}.tab"));
    let config = AnnotateConfig {
        input: input.to_string_lossy().into(),
        output: out.to_string_lossy().into(),
        gff3: vec![dir.join("genes.gff3").to_string_lossy().into()],
        fasta: Some(dir.join("genome.fa").to_string_lossy().into()),
        output_format: "tab".into(),
        pick: PickFlags::default(),
        hgvs,
        columns,
        distance: 5000,
        cache_dir: None,
        transcript_cache: None,
        sa_dir: None,
        sa_only: false,
        acmg: false,
        acmg_config: None,
        pick_order: None,
        functional_evidence: None,
        proband: None,
        mother: None,
        father: None,
        gene_list: None,
        explicit_alleles: false,
        qc_rules: None,
        show_progress: false,
        compat: Default::default(),
        exclude: Default::default(),
    };
    let before = ALLOCS.load(Ordering::Relaxed);
    run_annotate(config).expect("annotation should succeed");
    let allocations = ALLOCS.load(Ordering::Relaxed) - before;
    let text = std::fs::read_to_string(&out).unwrap();
    Run {
        allocations,
        records: copies * (vcf(1).lines().filter(|l| !l.starts_with('#')).count()),
        rows: text.lines().filter(|l| !l.starts_with('#')).count(),
    }
}

/// What 2,000 more records cost, for each record and each row.
fn marginal(columns: ColumnFlags, hgvs: bool) -> (f64, f64) {
    let dir = TempDir::new().unwrap();
    std::fs::copy(world("genes.gff3"), dir.path().join("genes.gff3")).unwrap();
    std::fs::copy(world("genome.fa"), dir.path().join("genome.fa")).unwrap();
    // writes the transcript cache beside the GFF3 and the index beside the FASTA, which every run
    // after it reads
    let _ = run(dir.path(), 1, columns, hgvs);
    let small = run(dir.path(), 1, columns, hgvs);
    let large = run(dir.path(), 5, columns, hgvs);
    assert!(
        large.rows > small.rows * 4,
        "{} and {} rows",
        small.rows,
        large.rows
    );
    let more = (large.allocations - small.allocations) as f64;
    (
        more / (large.records - small.records) as f64,
        more / (large.rows - small.rows) as f64,
    )
}

#[test]
fn the_annotation_loop_allocates_what_its_budget_says_for_each_record_and_row() {
    // the pool is one thread, so that the count does not depend on how the work was split
    rayon::ThreadPoolBuilder::new()
        .num_threads(1)
        .build_global()
        .expect("the pool is built before any run uses it");

    // `--symbol --canonical` and no HGVS: the consequences and the row
    let plain = ColumnFlags {
        symbol: true,
        canonical: true,
        ..Default::default()
    };
    let (per_record_plain, per_row_plain) = marginal(plain, false);
    // `--hgvs --hgvsg` with every optional column: what the reporter's runs have
    let full = ColumnFlags {
        symbol: true,
        biotype: true,
        canonical: true,
        numbers: true,
        variant_class: true,
        show_ref_allele: true,
        hgvsg: true,
        ..Default::default()
    };
    let (per_record_full, per_row_full) = marginal(full, true);

    eprintln!(
        "plain: {per_record_plain:.1} allocations a record, {per_row_plain:.1} a row; \
         with HGVS and every column: {per_record_full:.1} a record, {per_row_full:.1} a row"
    );
    assert!(
        per_row_plain <= PLAIN_PER_ROW,
        "a row of a run with no HGVS costs {per_row_plain:.1} allocations; the budget is {PLAIN_PER_ROW}"
    );
    assert!(
        per_row_full <= FULL_PER_ROW,
        "a row of a run with HGVS and every column costs {per_row_full:.1} allocations; the budget is {FULL_PER_ROW}"
    );
}

/// What a row of a run with `--symbol --canonical` and no HGVS costs, over the whole path from the
/// VCF line to the written line.
const PLAIN_PER_ROW: f64 = 33.0;
/// What a row of a run with `--hgvs --hgvsg` and every optional column costs.
const FULL_PER_ROW: f64 = 58.0;
