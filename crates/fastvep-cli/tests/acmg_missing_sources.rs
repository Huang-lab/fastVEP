//! What `--acmg` does when part of the reference stack is not loaded.
//!
//! A database directory that lacks a source does not make the run fail. Each criterion with no
//! data reports `evaluated: false`, and the run ends with a classification that rests on less
//! evidence and reads like any other. One criterion used to go the other way: PM2, "absent from
//! controls", read a variant with no gnomAD record as absent whether or not gnomAD was loaded, so
//! a run that never opened a gnomAD database met PM2 for every variant it classified. The
//! classifier is now told which sources are open (`LoadedSources`), and these tests drive the real
//! `run_annotate` and the real binary over a GFF3, a FASTA and a database built by `run_sa_build`,
//! so they fail if any link from the loaded database to the written criterion regresses.
//!
//! The reference is `ACGT` repeated and the gene is one coding exon, so the REF of each record
//! can be read off the pattern and a mismatch cannot muddy what the test is about.

use fastvep_annotate::AnnotationContext;
use fastvep_classification::{AcmgConfig, Source};
use fastvep_cli::pipeline::{run_annotate, run_sa_build, AnnotateConfig, PickFlags};
use serde_json::Value;
use std::fs;
use std::path::Path;
use std::process::Command;

const GFF3: &str = "##gff-version 3\n\
1\ttest\tgene\t1\t380\t.\t+\t.\tID=gene:ENSG_SRC;Name=SRCTEST;gene_name=SRCTEST;biotype=protein_coding\n\
1\ttest\tmRNA\t1\t380\t.\t+\t.\tID=transcript:ENST_SRC;Parent=gene:ENSG_SRC;biotype=protein_coding\n\
1\ttest\texon\t1\t380\t.\t+\t.\tID=exon:E1;Parent=transcript:ENST_SRC;rank=1\n\
1\ttest\tCDS\t1\t378\t.\t+\t0\tID=CDS:ENSP_SRC;Parent=transcript:ENST_SRC\n";

/// The base of the reference at a 1-based position.
fn reference_at(pos: usize) -> char {
    "ACGT".chars().nth((pos - 1) % 4).unwrap()
}

/// A missense change at `pos`: the base becomes an A (or a C where the reference has an A).
fn alt_at(pos: usize) -> char {
    if reference_at(pos) == 'A' {
        'C'
    } else {
        'A'
    }
}

struct Fixture {
    dir: tempfile::TempDir,
}

impl Fixture {
    /// A GFF3, a FASTA and a VCF with one variant at each of the given positions.
    fn new(positions: &[usize]) -> Self {
        let dir = tempfile::tempdir().unwrap();
        fs::write(dir.path().join("genes.gff3"), GFF3).unwrap();
        fs::write(
            dir.path().join("genome.fa"),
            format!(">1\n{}\n", "ACGT".repeat(95)),
        )
        .unwrap();
        let mut vcf =
            String::from("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
        for &p in positions {
            vcf.push_str(&format!(
                "1\t{p}\t.\t{}\t{}\t.\tPASS\t.\n",
                reference_at(p),
                alt_at(p)
            ));
        }
        fs::write(dir.path().join("in.vcf"), vcf).unwrap();
        Self { dir }
    }

    /// A gnomAD database holding one record, so that gnomAD is loaded and every other variant is
    /// one it does not know.
    fn with_gnomad_record(self, pos: usize, af: f64, ac: u64) -> Self {
        let sa = self.dir.path().join("sa");
        fs::create_dir_all(&sa).unwrap();
        let vcf = self.dir.path().join("gnomad.vcf");
        fs::write(
            &vcf,
            format!(
                "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
                 1\t{pos}\t.\t{}\t{}\t.\tPASS\tAF={af};AC={ac};AN=1000000;nhomalt=0\n",
                reference_at(pos),
                alt_at(pos)
            ),
        )
        .unwrap();
        run_sa_build(
            "gnomad",
            vcf.to_str().unwrap(),
            sa.join("gnomad").to_str().unwrap(),
            "GRCh38",
            None,
            &[],
            false,
        )
        .unwrap();
        self
    }

    fn path(&self, name: &str) -> String {
        self.dir.path().join(name).to_string_lossy().into()
    }

    fn annotate(&self, with_sa_dir: bool) -> Vec<Value> {
        let out = self.dir.path().join("out.json");
        run_annotate(AnnotateConfig {
            input: self.path("in.vcf"),
            output: out.to_string_lossy().into(),
            gff3: vec![self.path("genes.gff3")],
            fasta: Some(self.path("genome.fa")),
            output_format: "json".into(),
            pick: PickFlags::default(),
            hgvs: true,
            distance: 5000,
            cache_dir: None,
            transcript_cache: None,
            sa_dir: with_sa_dir.then(|| self.path("sa")),
            sa_only: false,
            acmg: true,
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
            columns: Default::default(),
            exclude: Default::default(),
        })
        .expect("annotation should succeed");
        serde_json::from_str(&fs::read_to_string(&out).unwrap()).unwrap()
    }
}

/// The criterion of a variant's first transcript row whose code starts with `prefix`.
fn criterion<'a>(variant: &'a Value, prefix: &str) -> &'a Value {
    let acmg = &variant["transcript_consequences"][0]["acmg"];
    acmg["criteria"]
        .as_array()
        .unwrap_or_else(|| panic!("an ACMG result with criteria: {variant}"))
        .iter()
        .find(|c| c["code"].as_str().is_some_and(|c| c.starts_with(prefix)))
        .unwrap_or_else(|| panic!("a criterion {prefix} in {acmg}"))
}

#[test]
fn with_no_gnomad_database_pm2_is_not_evaluated_and_does_not_read_a_missing_record_as_absence() {
    let fixture = Fixture::new(&[100]);
    let variants = fixture.annotate(false);
    let pm2 = criterion(&variants[0], "PM2");
    assert_eq!(pm2["evaluated"], false, "{pm2}");
    assert_eq!(pm2["met"], false, "{pm2}");
    let summary = pm2["summary"].as_str().unwrap();
    assert!(
        summary.contains("no gnomAD database is loaded"),
        "{summary}"
    );
    // BA1 cannot say anything either, and the class is what it was without PM2: VUS.
    assert_eq!(criterion(&variants[0], "BA1")["evaluated"], false);
    assert_eq!(
        variants[0]["transcript_consequences"][0]["acmg"]["classification"],
        "UncertainSignificance"
    );
}

#[test]
fn with_a_gnomad_database_loaded_a_variant_it_does_not_hold_is_absent_and_one_it_holds_is_common() {
    let fixture = Fixture::new(&[100, 200]).with_gnomad_record(200, 0.2, 200_000);
    let variants = fixture.annotate(true);

    let absent = criterion(&variants[0], "PM2");
    assert_eq!(absent["evaluated"], true, "{absent}");
    assert_eq!(absent["met"], true, "{absent}");
    assert!(
        absent["summary"]
            .as_str()
            .unwrap()
            .contains("Absent in gnomAD"),
        "{absent}"
    );

    let common = &variants[1];
    assert_eq!(criterion(common, "PM2")["met"], false);
    assert_eq!(criterion(common, "BA1")["met"], true);
    assert_eq!(
        common["transcript_consequences"][0]["acmg"]["classification"],
        "Benign"
    );
}

fn run_binary(fixture: &Fixture, extra: &[&str]) -> (std::process::Output, Vec<Value>) {
    let out = fixture.dir.path().join("cli_out.json");
    let output = Command::new(env!("CARGO_BIN_EXE_fastvep"))
        .arg("annotate")
        .arg("-i")
        .arg(fixture.path("in.vcf"))
        .arg("-o")
        .arg(&out)
        .arg("--gff3")
        .arg(fixture.path("genes.gff3"))
        .arg("--fasta")
        .arg(fixture.path("genome.fa"))
        .args(["--output-format", "json", "--no-progress"])
        .args(extra)
        .output()
        .expect("the fastvep binary runs");
    let text = fs::read_to_string(Path::new(&out)).unwrap_or_default();
    let variants = serde_json::from_str(&text).unwrap_or_default();
    (output, variants)
}

#[test]
fn the_command_line_says_which_sources_an_acmg_run_is_missing() {
    let fixture = Fixture::new(&[100]).with_gnomad_record(200, 0.2, 200_000);
    let sa = fixture.path("sa");
    let (output, _) = run_binary(&fixture, &["--acmg", "--sa-dir", &sa]);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(output.status.success(), "{stderr}");
    assert!(
        stderr.contains("--acmg is running on an incomplete reference stack"),
        "{stderr}"
    );
    // gnomAD is loaded, so it is not named; the eight that are not are.
    assert!(stderr.contains("8 of the 9 sources"), "{stderr}");
    assert!(!stderr.contains("gnomAD allele frequencies ("), "{stderr}");
    assert!(
        stderr.contains("REVEL missense scores (PP3, BP4)"),
        "{stderr}"
    );
}

#[test]
fn a_run_without_acmg_has_nothing_to_warn_about() {
    let fixture = Fixture::new(&[100]);
    let (output, variants) = run_binary(&fixture, &[]);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(output.status.success(), "{stderr}");
    assert!(!stderr.contains("incomplete reference stack"), "{stderr}");
    assert_eq!(variants.len(), 1);
}

/// The library behind both web servers has a per-variant loop of its own, and the two loops feed
/// one classifier, so the library gets the test the command line did.
#[test]
fn the_library_path_that_the_web_servers_use_reads_the_loaded_sources_as_the_command_line_does() {
    let fixture = Fixture::new(&[100]).with_gnomad_record(200, 0.2, 200_000);
    let annotate = |sa_dir: Option<String>| {
        let mut ctx = AnnotationContext::new(
            Some(&fixture.path("genes.gff3")),
            Some(&fixture.path("genome.fa")),
            sa_dir.as_deref(),
            5000,
        )
        .expect("the fixture should load");
        ctx.hgvs = true;
        let sources = ctx.loaded_sources();
        let acmg = AcmgConfig::default();
        let vcf = fs::read_to_string(fixture.path("in.vcf")).unwrap();
        let results = ctx
            .annotate_vcf_text_with_acmg(&vcf, None, Some(&acmg))
            .expect("annotation should succeed");
        (sources, results)
    };

    let (sources, results) = annotate(None);
    assert!(Source::ALL.iter().all(|s| !sources.has(*s)));
    let pm2 = criterion(&results[0], "PM2");
    assert_eq!(pm2["evaluated"], false, "{pm2}");
    assert!(pm2["summary"]
        .as_str()
        .unwrap()
        .contains("no gnomAD database is loaded"));

    let (sources, results) = annotate(Some(fixture.path("sa")));
    assert!(sources.has(Source::Gnomad));
    assert_eq!(sources.missing().count(), 8);
    let pm2 = criterion(&results[0], "PM2");
    assert_eq!(pm2["met"], true, "{pm2}");
}
