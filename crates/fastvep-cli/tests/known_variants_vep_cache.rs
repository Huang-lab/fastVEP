//! `--cache-dir` end to end: a VCF in, the tab file a user reads out, against a
//! variation cache laid out the way VEP converts it. Issue #139.
//!
//! Driven through `run_annotate` for the reason `flag_pick.rs` gives: the unit
//! tests build records by hand, so they cannot see a flag that never reaches the
//! lookup, a lookup whose answer is dropped before the writer, or a cache that
//! was refused in the library but ignored by the command.
//!
//! The cache is `fastvep-cache/fixtures/vep_cache/csi`: `all_vars.gz` bgzipped
//! and indexed with `tabix -C`, so it has a `.csi` and no `.tbi`, as a cache
//! converted by VEP does. Before the fix every `Existing_variation` cell here
//! was `-`, with no message.

use std::collections::BTreeMap;
use std::fs::File;
use std::io::Write;
use std::path::{Path, PathBuf};

use fastvep_cli::pipeline::{run_annotate, AnnotateConfig, PickFlags};
use tempfile::TempDir;

fn cache(kind: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../fastvep-cache/fixtures/vep_cache")
        .join(kind)
}

/// One lncRNA a little downstream of the test sites, so each chr10 variant gets
/// a transcript row (upstream_gene_variant) and `chrM` falls to an intergenic one.
fn write_gff3(dir: &Path) -> PathBuf {
    let path = dir.join("genes.gff3");
    std::fs::write(
        &path,
        "chr10\ttest\tgene\t11730\t14000\t.\t+\t.\tID=gene:ENSG_X;Name=GENEX;biotype=lncRNA\n\
         chr10\ttest\tlnc_RNA\t11730\t14000\t.\t+\t.\tID=transcript:ENST_X;Parent=gene:ENSG_X;biotype=lncRNA\n\
         chr10\ttest\texon\t11730\t12000\t.\t+\t.\tID=exon:E1;Parent=transcript:ENST_X;rank=1\n\
         chr10\ttest\texon\t13000\t14000\t.\t+\t.\tID=exon:E2;Parent=transcript:ENST_X;rank=2\n",
    )
    .unwrap();
    path
}

/// The shapes of issue #139, each as the VCF a caller writes (anchored).
fn write_vcf(dir: &Path) -> PathBuf {
    let path = dir.join("in.vcf");
    let mut f = File::create(&path).unwrap();
    writeln!(f, "##fileformat=VCFv4.2").unwrap();
    writeln!(f, "##contig=<ID=chr10,length=133797422>").unwrap();
    writeln!(f, "##contig=<ID=chrM,length=16569>").unwrap();
    writeln!(f, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO").unwrap();
    // AACCCT>A: ACCCT/- at 10401; the cache has it anchored
    writeln!(f, "chr10\t10400\t.\tAACCCT\tA\t.\t.\t.").unwrap();
    // TAACCC>T: AACCC/- at 10400; the cache has exactly that
    writeln!(f, "chr10\t10399\t.\tTAACCC\tT\t.\t.\t.").unwrap();
    // C>CT: -/T at 10409, end 10408; the cache stores it with end < start
    writeln!(f, "chr10\t10408\t.\tC\tCT\t.\t.\t.").unwrap();
    // G>GA: -/A at 10421; the cache has it anchored
    writeln!(f, "chr10\t10420\t.\tG\tGA\t.\t.\t.").unwrap();
    // the SNV site, as two records and as one two-allele record
    writeln!(f, "chr10\t10414\t.\tC\tA\t.\t.\t.").unwrap();
    writeln!(f, "chr10\t10414\t.\tC\tG\t.\t.\t.").unwrap();
    writeln!(f, "chr10\t10414\t.\tC\tA,G\t.\t.\t.").unwrap();
    // the mitochondrion, `chrM` in the input and `MT` in the cache
    writeln!(f, "chrM\t100\t.\tA\tG\t.\t.\t.").unwrap();
    path
}

fn config(input: &Path, out: &Path) -> AnnotateConfig {
    AnnotateConfig {
        input: input.to_string_lossy().into(),
        output: out.to_string_lossy().into(),
        gff3: vec![],
        fasta: None,
        output_format: "tab".into(),
        pick: PickFlags::default(),
        hgvs: false,
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
    }
}

fn run(dir: &Path, cache_dir: &Path) -> anyhow::Result<String> {
    let out = dir.join(format!(
        "out-{}.tab",
        cache_dir.file_name().unwrap().to_string_lossy()
    ));
    run_annotate(AnnotateConfig {
        gff3: vec![write_gff3(dir).to_string_lossy().into()],
        cache_dir: Some(cache_dir.to_string_lossy().into()),
        ..config(&write_vcf(dir), &out)
    })?;
    Ok(std::fs::read_to_string(&out).unwrap())
}

/// `(Location, Allele)` to the `Existing_variation` of every row with that key,
/// read by the declared column names.
fn existing_by_site(tab: &str) -> BTreeMap<(String, String), Vec<String>> {
    let header = tab
        .lines()
        .find(|l| l.starts_with("#Uploaded_variation"))
        .unwrap_or_else(|| panic!("no column header in:\n{tab}"));
    let col = |name: &str| {
        header
            .trim_start_matches('#')
            .split('\t')
            .position(|c| c == name)
            .unwrap_or_else(|| panic!("no {name} column in {header}"))
    };
    let (loc, allele, existing) = (col("Location"), col("Allele"), col("Existing_variation"));
    let mut by_site: BTreeMap<(String, String), Vec<String>> = BTreeMap::new();
    for row in tab.lines().filter(|l| !l.starts_with('#')) {
        let f: Vec<&str> = row.split('\t').collect();
        by_site
            .entry((f[loc].to_string(), f[allele].to_string()))
            .or_default()
            .push(f[existing].to_string());
    }
    by_site
}

fn expected() -> Vec<((&'static str, &'static str), &'static str)> {
    vec![
        // an anchored deletion in the cache, found from the minimised input
        (("chr10:10401-10405", "-"), "rsFIX_DEL_ANCH"),
        // a minimal deletion matches only its own record, not the anchored one
        (("chr10:10400-10404", "-"), "rsFIX_DEL_MIN"),
        // an insertion the cache stores with end below start
        (("chr10:10409-10408", "T"), "rsFIX_INS_MIN"),
        // an insertion the cache stores anchored
        (("chr10:10421-10420", "A"), "rsFIX_INS_ANCH"),
        // germline first, the somatic record without alleles last, the failed
        // record not at all
        (("chr10:10414", "A"), "rsFIX_SNV,rsFIX_SNV_MULTI,COSVFIX1"),
        (("chr10:10414", "G"), "rsFIX_SNV_MULTI,COSVFIX1"),
        // the input says chrM, the cache directory is MT
        (("chrM:100", "G"), "rsFIX_MT"),
    ]
}

#[test]
fn a_vep_converted_cache_gives_the_known_variants_of_every_shape_the_issue_reports() {
    let dir = TempDir::new().unwrap();
    let tab = run(dir.path(), &cache("csi")).expect("annotation should succeed");
    let by_site = existing_by_site(&tab);

    for ((location, allele), want) in expected() {
        let got = by_site
            .get(&(location.to_string(), allele.to_string()))
            .unwrap_or_else(|| panic!("no row for {location} {allele} in:\n{tab}"));
        // The two-allele record and the single-allele records give the same
        // answer for the same allele.
        for cell in got {
            assert_eq!(cell, want, "{location} {allele}\n{tab}");
        }
    }
}

#[test]
fn the_older_tbi_form_of_the_cache_gives_the_same_output() {
    let dir = TempDir::new().unwrap();
    let from_csi = run(dir.path(), &cache("csi")).unwrap();
    let from_tbi = run(dir.path(), &cache("tbi")).unwrap();
    // `Existing_variation` is what the index could change, but the whole file
    // is compared so that nothing else shifts either.
    assert_eq!(existing_by_site(&from_csi), existing_by_site(&from_tbi));
    assert_eq!(from_csi, from_tbi);
}

/// The cache in question cannot answer anything, so the run must say so and
/// stop. It used to run to the end and print `-` in every cell.
#[test]
fn a_cache_that_cannot_be_searched_stops_the_run() {
    let dir = TempDir::new().unwrap();
    let err = run(dir.path(), &cache("noindex")).expect_err("a cache with no index is refused");
    let message = format!("{err:#}");
    assert!(message.contains("--cache-dir"), "{message}");
    assert!(message.contains("no .csi or .tbi index"), "{message}");
    assert!(message.contains("tabix -C -s 1 -b 5 -e 5"), "{message}");
}

#[test]
fn an_unconverted_cache_directory_stops_the_run() {
    let dir = TempDir::new().unwrap();
    // The layout of an unconverted cache: info.txt and chromosome directories
    // holding the `{start}-{end}_var.gz` slices, and no all_vars.gz.
    let unconverted = dir.path().join("unconverted");
    std::fs::create_dir_all(unconverted.join("10")).unwrap();
    std::fs::write(unconverted.join("10/10001-1010000_var.gz"), b"").unwrap();
    std::fs::copy(cache("csi").join("info.txt"), unconverted.join("info.txt")).unwrap();

    let err = run(dir.path(), &unconverted).expect_err("an unconverted cache is refused");
    let message = format!("{err:#}");
    assert!(
        message.contains("holds no `{chr}/all_vars.gz`"),
        "{message}"
    );
    assert!(message.contains("unconverted"), "{message}");
}
