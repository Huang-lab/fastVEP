//! The two per-variant loops, `fastvep annotate`'s and the library's (the one the web server runs),
//! write the same consequences and descriptions for the same input.
//!
//! They are two copies of one job, and they have drifted apart before: a rule added to one was
//! missing from the other until a probe found it. This runs the CLI over each of the committed
//! worlds and the library over the same file and gene models, and compares what each says of every
//! (variant, allele, transcript): the consequence terms, the impact, the residues and codons, the
//! positions, `HGVSc`, `HGVSp` and `HGVSg`.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use fastvep_annotate::AnnotationContext;
use fastvep_cli::pipeline::{run_annotate, AnnotateConfig, ColumnFlags, PickFlags};
use serde_json::Value;
use tempfile::TempDir;

/// What each loop writes for one (allele, transcript), the fields both write.
const FIELDS: &[&str] = &[
    "consequence_terms",
    "impact",
    "amino_acids",
    "codons",
    "cdna_start",
    "cdna_end",
    "cds_start",
    "cds_end",
    "protein_start",
    "protein_end",
    "exon",
    "intron",
    "hgvsc",
    "hgvsp",
    "hgvsg",
    "hgvs_offset",
    "distance",
    "gene_id",
    "gene_symbol",
    "biotype",
    "canonical",
    "strand",
    "protein_id",
];

fn fixture(world: &str, name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("fixtures")
        .join(world)
        .join(name)
}

type Key = (String, u64, String, String, String);

/// `(chromosome, start, allele_string, transcript, allele)` to the fields above, from the JSON
/// either loop writes.
fn by_row(variants: &[Value]) -> BTreeMap<Key, BTreeMap<&'static str, Value>> {
    let mut rows = BTreeMap::new();
    for v in variants {
        let (chrom, start, alleles) = (
            v["seq_region_name"].as_str().unwrap().to_string(),
            v["start"].as_u64().unwrap(),
            v["allele_string"].as_str().unwrap().to_string(),
        );
        for tc in v["transcript_consequences"]
            .as_array()
            .into_iter()
            .flatten()
        {
            let key = (
                chrom.clone(),
                start,
                alleles.clone(),
                tc["transcript_id"].as_str().unwrap_or("-").to_string(),
                tc["variant_allele"].as_str().unwrap_or("-").to_string(),
            );
            let fields = FIELDS
                .iter()
                .map(|f| (*f, tc.get(*f).cloned().unwrap_or(Value::Null)))
                .collect();
            rows.insert(key, fields);
        }
    }
    rows
}

fn cli_rows(
    world: &str,
    vcf: &Path,
    distance: u64,
) -> BTreeMap<Key, BTreeMap<&'static str, Value>> {
    let dir = TempDir::new().unwrap();
    let (gff3, fasta) = (dir.path().join("genes.gff3"), dir.path().join("genome.fa"));
    std::fs::copy(fixture(world, "genes.gff3"), &gff3).unwrap();
    std::fs::copy(fixture(world, "genome.fa"), &fasta).unwrap();
    let out = dir.path().join("out.json");
    run_annotate(AnnotateConfig {
        input: vcf.to_string_lossy().into(),
        output: out.to_string_lossy().into(),
        gff3: vec![gff3.to_string_lossy().into()],
        fasta: Some(fasta.to_string_lossy().into()),
        output_format: "json".into(),
        pick: PickFlags::default(),
        hgvs: true,
        columns: ColumnFlags {
            hgvsg: true,
            ..Default::default()
        },
        distance,
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
    })
    .expect("annotation should succeed");
    let json: Value = serde_json::from_str(&std::fs::read_to_string(&out).unwrap()).unwrap();
    by_row(json.as_array().expect("an array of variants"))
}

fn library_rows(
    world: &str,
    vcf_text: &str,
    distance: u64,
) -> BTreeMap<Key, BTreeMap<&'static str, Value>> {
    let dir = TempDir::new().unwrap();
    let (gff3, fasta) = (dir.path().join("genes.gff3"), dir.path().join("genome.fa"));
    std::fs::copy(fixture(world, "genes.gff3"), &gff3).unwrap();
    std::fs::copy(fixture(world, "genome.fa"), &fasta).unwrap();
    let mut ctx = AnnotationContext::new(
        Some(gff3.to_str().unwrap()),
        Some(fasta.to_str().unwrap()),
        None,
        distance,
    )
    .expect("the world should load");
    ctx.hgvs = true;
    let variants = ctx
        .annotate_vcf_text(vcf_text, None)
        .expect("annotation should succeed");
    by_row(&variants)
}

fn compare(world: &str, vcf: &Path, distance: u64, at_least: usize) {
    let text = std::fs::read_to_string(vcf).unwrap();
    let (cli, lib) = (
        cli_rows(world, vcf, distance),
        library_rows(world, &text, distance),
    );
    assert!(cli.len() >= at_least, "{world}: only {} rows", cli.len());
    let mut missing = Vec::new();
    let mut differ: BTreeMap<&str, Vec<String>> = BTreeMap::new();
    for (key, a) in &cli {
        let Some(b) = lib.get(key) else {
            missing.push(format!("{key:?}"));
            continue;
        };
        for f in FIELDS {
            if a[f] != b[f] {
                differ
                    .entry(f)
                    .or_default()
                    .push(format!("{key:?}: CLI {} library {}", a[f], b[f]));
            }
        }
    }
    let extra: Vec<_> = lib
        .keys()
        .filter(|k| !cli.contains_key(*k))
        .take(5)
        .collect();
    assert!(
        missing.is_empty(),
        "{world}: rows only the CLI writes, e.g. {:?}",
        &missing[..missing.len().min(5)]
    );
    assert!(
        extra.is_empty(),
        "{world}: rows only the library writes, e.g. {extra:?}"
    );
    for (field, rows) in &differ {
        eprintln!("{world} {field}: {} rows, e.g. {}", rows.len(), rows[0]);
    }
    assert!(
        differ.is_empty(),
        "{world}: the loops differ in {:?}",
        differ.keys()
    );
}

#[test]
fn the_two_loops_agree_on_alleles_that_are_not_plain_bases() {
    // `*`, `.`, `<NON_REF>`, ambiguity codes, symbolic alleles, several ALTs, on the generated world
    compare(
        "vep_tab_world",
        &fixture("vep_tab_alleles", "in.vcf"),
        5000,
        300,
    );
}

/// A VCF in lower case, which a soft-masked reference gives, with an ALT that repeats another
/// and a REF that is not the genome's, is read the same by both.
#[test]
fn the_two_loops_agree_on_a_vcf_that_is_not_what_a_caller_writes() {
    let dir = TempDir::new().unwrap();
    let source = std::fs::read_to_string(fixture("vep_tab_edges", "in.vcf")).unwrap();
    let mut out = String::new();
    for (n, line) in source.lines().enumerate() {
        if line.starts_with('#') {
            out.push_str(line);
            out.push('\n');
            continue;
        }
        let mut cells: Vec<String> = line.split('\t').map(String::from).collect();
        // lower case REF and ALT, as a soft-masked reference makes them
        cells[3] = cells[3].to_lowercase();
        cells[4] = cells[4].to_lowercase();
        if n % 7 == 0 {
            // an ALT that repeats the one before it
            cells[4] = format!("{},{}", cells[4], cells[4]);
        }
        if n % 11 == 0 && cells[3].len() == 1 {
            // a REF that is another base than the genome's
            cells[3] = if cells[3] == "a" {
                "c".into()
            } else {
                "a".into()
            };
        }
        out.push_str(&cells.join("\t"));
        out.push('\n');
    }
    let vcf = dir.path().join("odd.vcf");
    std::fs::write(&vcf, out).unwrap();
    compare("vep_tab_world", &vcf, 5000, 1000);
}

#[test]
fn the_two_loops_agree_on_the_world_of_frameshift_introns() {
    compare(
        "vep_tab_fsintron",
        &fixture("vep_tab_fsintron", "in.vcf"),
        5000,
        1000,
    );
}

#[test]
fn the_two_loops_agree_on_the_world_of_start_codons() {
    compare(
        "vep_tab_initiator",
        &fixture("vep_tab_initiator", "in.vcf"),
        0,
        1000,
    );
}

#[test]
fn the_two_loops_agree_on_the_world_of_repeats_over_introns() {
    compare(
        "vep_tab_intron_repeat",
        &fixture("vep_tab_intron_repeat", "in.vcf"),
        0,
        1000,
    );
}

#[test]
fn the_two_loops_agree_on_the_edges_of_the_generated_world() {
    compare(
        "vep_tab_world",
        &fixture("vep_tab_edges", "in.vcf"),
        5000,
        1000,
    );
}
