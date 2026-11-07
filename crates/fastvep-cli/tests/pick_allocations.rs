//! Pins what a reducing `--pick*` costs, in heap allocations, when a variant has many transcripts
//! around it and the pick keeps one per gene.
//!
//! A variant of a run with `--distance 500000` has hundreds of transcripts around it, and every
//! allocation made for a transcript the pick then leaves out is made hundreds of times for a
//! variant and millions of times for a run. The pick is settled on what the predictor said, before
//! any row is built, and a transcript the variant is merely beside is not even predicted unless it
//! is chosen; so the work of a variant follows the transcripts that are kept, and not the ones that
//! are around. This holds that as a count: the same world is annotated twice, with a handful of
//! records and with a few hundred, and what the extra records cost is divided out.
//!
//! The world has 60 genes of 3 transcripts each within reach of every variant, so a variant has 180
//! transcripts around it and `--pick-allele-gene` keeps 60 rows of them. A variant is either inside
//! a gene, where its 3 transcripts are predicted in full and the other 177 are beside it, or in the
//! gap between two, where all 180 are beside it; each is a count of its own.
//!
//! This file installs a counting global allocator. It counts the allocations of the threads that do
//! the annotating, which are this one and the rayon pool's, and of no other: the test harness has a
//! thread of its own that allocates while a test runs, and a count of the whole process took those
//! into the window that is measured, which failed a test of this kind on a CI runner. The rayon pool
//! is one thread, so that the count does not depend on how the work was split.

use std::alloc::{GlobalAlloc, Layout, System};
use std::cell::Cell;
use std::fmt::Write as _;
use std::path::Path;
use std::sync::atomic::{AtomicUsize, Ordering};

use fastvep_cli::pipeline::{run_annotate, AnnotateConfig, ColumnFlags, PickFlags};
use tempfile::TempDir;

static ALLOCS: AtomicUsize = AtomicUsize::new(0);

thread_local! {
    // Whether this thread's allocations are counted. A `Cell<bool>` with a constant initial value
    // allocates nothing and has no destructor, so it is safe to read inside the allocator.
    static COUNTED: Cell<bool> = const { Cell::new(false) };
}

fn count() {
    if COUNTED.with(Cell::get) {
        ALLOCS.fetch_add(1, Ordering::Relaxed);
    }
}

struct Counting;

unsafe impl GlobalAlloc for Counting {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        count();
        unsafe { System.alloc(layout) }
    }
    unsafe fn dealloc(&self, ptr: *mut u8, layout: Layout) {
        unsafe { System.dealloc(ptr, layout) }
    }
    unsafe fn realloc(&self, ptr: *mut u8, layout: Layout, new_size: usize) -> *mut u8 {
        count();
        unsafe { System.realloc(ptr, layout, new_size) }
    }
}

#[global_allocator]
static ALLOCATOR: Counting = Counting;

const GENES: u64 = 60;
const TRANSCRIPTS: u64 = 3;

/// 60 genes along 120 kb, each of 3 coding transcripts that share their exon, and a genome to
/// read them from.
fn world(dir: &Path) {
    let mut gff3 = String::new();
    for g in 0..GENES {
        let start = 10_000 + g * 2_000;
        let end = start + 999;
        let strand = if g % 2 == 0 { '+' } else { '-' };
        writeln!(
            gff3,
            "1\ttest\tgene\t{start}\t{end}\t.\t{strand}\t.\tID=gene:ENSG{g:02};Name=GENE{g:02};biotype=protein_coding"
        )
        .unwrap();
        for t in 0..TRANSCRIPTS {
            let id = format!("ENST{g:02}_{t}");
            let canonical = if t == 0 { ";tag=Ensembl_canonical" } else { "" };
            writeln!(
                gff3,
                "1\ttest\tmRNA\t{start}\t{end}\t.\t{strand}\t.\tID=transcript:{id};Parent=gene:ENSG{g:02};biotype=protein_coding{canonical}"
            )
            .unwrap();
            writeln!(
                gff3,
                "1\ttest\texon\t{start}\t{end}\t.\t{strand}\t.\tID=exon:{id};Parent=transcript:{id};rank=1"
            )
            .unwrap();
            writeln!(
                gff3,
                "1\ttest\tCDS\t{start}\t{end}\t.\t{strand}\t0\tID=CDS:{id};Parent=transcript:{id}"
            )
            .unwrap();
        }
    }
    std::fs::write(dir.join("genes.gff3"), gff3).unwrap();

    let genome: String = (0..140_000u64)
        .map(|i| match (i * 7 + i / 3) % 4 {
            0 => 'A',
            1 => 'C',
            2 => 'G',
            _ => 'T',
        })
        .collect();
    let mut fasta = String::from(">1\n");
    for line in genome.as_bytes().chunks(60) {
        fasta.push_str(std::str::from_utf8(line).unwrap());
        fasta.push('\n');
    }
    std::fs::write(dir.join("genome.fa"), fasta).unwrap();
}

/// Where the substitutions of a run are.
#[derive(Clone, Copy, Debug)]
enum Shape {
    /// Inside a gene, in the exon its 3 transcripts share: those 3 are predicted in full, and the
    /// other 177 transcripts are beside the variant.
    InAGene,
    /// In the gap between two genes: all 180 transcripts are beside the variant.
    BetweenGenes,
}

/// `records` substitutions spread over the genes, each of the given shape.
fn vcf(dir: &Path, records: u64, shape: Shape) -> std::path::PathBuf {
    let genome = std::fs::read_to_string(dir.join("genome.fa")).unwrap();
    let genome: Vec<u8> = genome.lines().skip(1).flat_map(|l| l.bytes()).collect();
    let mut text =
        String::from("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
    for i in 0..records {
        // the gene of the record, and a place in it or in the 1,000 bases after it
        let gene = (i * 7) % GENES;
        let offset = (i * 997) % 1_000;
        let pos = 10_000
            + gene * 2_000
            + offset
            + match shape {
                Shape::InAGene => 0,
                Shape::BetweenGenes => 1_000,
            };
        let reference = genome[(pos - 1) as usize] as char;
        let alt = if reference == 'A' { 'C' } else { 'A' };
        writeln!(text, "1\t{pos}\t.\t{reference}\t{alt}\t.\t.\t.").unwrap();
    }
    let path = dir.join(format!("in_{records}.vcf"));
    std::fs::write(&path, text).unwrap();
    path
}

struct Run {
    allocations: usize,
    rows: usize,
}

fn run(dir: &Path, records: u64, shape: Shape, pick: PickFlags, hgvs: bool) -> Run {
    let input = vcf(dir, records, shape);
    let out = dir.join(format!("out_{records}.tab"));
    let config = AnnotateConfig {
        input: input.to_string_lossy().into(),
        output: out.to_string_lossy().into(),
        gff3: vec![dir.join("genes.gff3").to_string_lossy().into()],
        fasta: Some(dir.join("genome.fa").to_string_lossy().into()),
        output_format: "tab".into(),
        pick,
        hgvs,
        columns: ColumnFlags {
            symbol: true,
            canonical: true,
            ..Default::default()
        },
        distance: 500_000,
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
        rows: text.lines().filter(|l| !l.starts_with('#')).count(),
    }
}

/// What a record costs, in allocations: what the extra records of a larger run cost, divided out.
fn per_record(shape: Shape, pick: PickFlags, hgvs: bool) -> (f64, f64) {
    let dir = TempDir::new().unwrap();
    world(dir.path());
    // writes the transcript cache beside the GFF3, which every run after it reads
    let _ = run(dir.path(), 10, shape, pick, hgvs);
    let small = run(dir.path(), 20, shape, pick, hgvs);
    let large = run(dir.path(), 220, shape, pick, hgvs);
    assert!(
        large.rows > small.rows * 5,
        "{} and {} rows",
        small.rows,
        large.rows
    );
    let more = (large.allocations - small.allocations) as f64;
    (more / 200.0, (large.rows - small.rows) as f64 / 200.0)
}

#[test]
fn a_reducing_pick_allocates_for_the_transcripts_it_keeps_and_not_for_the_ones_around() {
    // the pool is one thread, so that the count does not depend on how the work was split, and
    // its thread is one of those that are counted
    rayon::ThreadPoolBuilder::new()
        .num_threads(1)
        .start_handler(|_| COUNTED.with(|counted| counted.set(true)))
        .build_global()
        .expect("the pool is built before any run uses it");
    COUNTED.with(|counted| counted.set(true));

    let pick = PickFlags {
        pick_allele_gene: true,
        ..Default::default()
    };
    for (shape, budget) in [
        (Shape::InAGene, IN_A_GENE),
        (Shape::BetweenGenes, BETWEEN_GENES),
    ] {
        let (per_record, rows_per_record) = per_record(shape, pick, true);
        // every transcript is around every variant: 180 of them, and one row for each of 60 genes
        assert!(
            (rows_per_record - GENES as f64).abs() < 1.0,
            "{shape:?}: {rows_per_record:.1} rows a record; the world is meant to keep {GENES}"
        );
        eprintln!(
            "--pick-allele-gene --hgvs, {shape:?}: {per_record:.1} allocations a record, {rows_per_record:.1} rows a record"
        );
        assert!(
            per_record <= budget,
            "{shape:?}: a record with {} transcripts around it and {GENES} rows kept costs {per_record:.1} allocations; the budget is {budget}",
            GENES * TRANSCRIPTS
        );
    }
}

/// What a record costs, over the whole path from the VCF line to the written lines: 897 today for a
/// record in a gene and 839 for one between genes, of which about 14 for each of the 60 rows kept.
/// Building a row for each of the 180 transcripts around the variant, as the pick used to be made
/// after the rows, was more than twice that. The budgets leave 6 per cent for what the standard
/// library's growth of a table or a vector does from one release to the next.
const IN_A_GENE: f64 = 950.0;
const BETWEEN_GENES: f64 = 890.0;
