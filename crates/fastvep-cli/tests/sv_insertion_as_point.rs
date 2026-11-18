//! `--sv-insertion-as-point`: an insertion is the point between `POS` and `POS + 1`, and not the
//! interval of its `SVLEN` bases that Ensembl VEP reads.
//!
//! The interval can reach an exon that the point does not, so a long insertion in an intron has
//! terms for the coding sequence beside the intron's own. The option is off by default, so that the
//! rows are VEP's, and both per-variant loops (the CLI's and the library's, which the web server
//! runs) read the same way when it is on.

use std::path::{Path, PathBuf};

use fastvep_annotate::AnnotationContext;
use fastvep_cli::pipeline::{run_annotate, AnnotateConfig, ColumnFlags, CompatFlags, PickFlags};
use tempfile::TempDir;

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("fixtures/vep_tab_sv")
        .join(name)
}

/// An insertion of 100 bases in the intron of SVA between its first two exons (`1130-1299` and
/// `1300-1420`, the second of them coding), 50 bases from the second: as an interval it covers
/// `1251-1350`, and so 51 bases of that exon.
const VCF: &str = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
                   chr1\t1250\tintron_insertion\tN\t<INS>\t.\t.\tSVLEN=100\n\
                   chr1\t1250\tintron_mobile_element\tN\t<INS:ME:ALU>\t.\t.\tSVLEN=300;END=1250\n\
                   chr1\t1250\tduplication_beside\tN\t<DUP>\t.\t.\tEND=1350\n";

fn terms_of(consequence_terms: &[String]) -> String {
    consequence_terms.join(",")
}

/// The terms of ENST00000000401 for each record, by the CLI, with or without the option.
fn cli(as_point: bool) -> Vec<(String, String)> {
    let dir = TempDir::new().unwrap();
    let (gff3, fasta) = (dir.path().join("genes.gff3"), dir.path().join("genome.fa"));
    std::fs::copy(fixture("genes.gff3"), &gff3).unwrap();
    std::fs::copy(fixture("genome.fa"), &fasta).unwrap();
    let vcf = dir.path().join("in.vcf");
    std::fs::write(&vcf, VCF).unwrap();
    let out = dir.path().join("out.tab");
    run_annotate(AnnotateConfig {
        input: vcf.to_string_lossy().into(),
        output: out.to_string_lossy().into(),
        gff3: vec![gff3.to_string_lossy().into()],
        fasta: Some(fasta.to_string_lossy().into()),
        output_format: "tab".into(),
        pick: PickFlags::default(),
        hgvs: false,
        columns: ColumnFlags::default(),
        distance: 1000,
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
        compat: CompatFlags {
            sv_insertion_as_point: as_point,
            ..Default::default()
        },
        exclude: Default::default(),
    })
    .unwrap();
    std::fs::read_to_string(&out)
        .unwrap()
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| l.split('\t').collect::<Vec<_>>())
        .filter(|c| c[4] == "ENST00000000401")
        .map(|c| (c[0].to_string(), c[6].to_string()))
        .collect()
}

/// The same by the library loop.
fn library(as_point: bool) -> Vec<(String, String)> {
    let dir = TempDir::new().unwrap();
    let (gff3, fasta) = (dir.path().join("genes.gff3"), dir.path().join("genome.fa"));
    std::fs::copy(fixture("genes.gff3"), &gff3).unwrap();
    std::fs::copy(fixture("genome.fa"), &fasta).unwrap();
    let mut ctx = AnnotationContext::new(
        Some(gff3.to_str().unwrap()),
        Some(fasta.to_str().unwrap()),
        None,
        1000,
    )
    .unwrap();
    ctx.sv_insertion_as_point = as_point;
    let variants = ctx
        .annotate_vcf_text(VCF, PickFlags::default().resolve())
        .unwrap();
    variants
        .iter()
        .flat_map(|v| {
            v["transcript_consequences"]
                .as_array()
                .unwrap()
                .iter()
                .filter(|t| t["transcript_id"] == "ENST00000000401")
                .map(|t| {
                    let terms: Vec<String> = t["consequence_terms"]
                        .as_array()
                        .unwrap()
                        .iter()
                        .map(|c| c.as_str().unwrap().to_string())
                        .collect();
                    (v["id"].as_str().unwrap().to_string(), terms_of(&terms))
                })
                .collect::<Vec<_>>()
        })
        .collect()
}

#[test]
fn an_insertion_reaches_the_exon_by_its_interval_and_not_by_its_point() {
    let terms = |rows: &[(String, String)], id: &str| -> String {
        rows.iter()
            .find(|r| r.0 == id)
            .map(|r| r.1.clone())
            .unwrap()
    };
    let (interval, point) = (cli(false), cli(true));
    for id in ["intron_insertion", "intron_mobile_element"] {
        // VEP's reading: the 51 bases of the coding exon that the interval covers are the variant's
        assert!(
            terms(&interval, id).contains("coding_sequence_variant"),
            "{id}: {}",
            terms(&interval, id)
        );
        // the point is in the intron, 50 bases from the exon
        assert_eq!(terms(&point, id), "intron_variant", "{id}");
    }
    // a duplication is not an insertion: it keeps its interval
    assert_eq!(
        terms(&interval, "duplication_beside"),
        terms(&point, "duplication_beside")
    );
}

#[test]
fn the_library_loop_reads_an_insertion_as_the_cli_does_with_and_without_the_option() {
    for as_point in [false, true] {
        let (cli, library) = (cli(as_point), library(as_point));
        assert_eq!(cli.len(), 3, "{cli:?}");
        for (id, terms) in &cli {
            let other = library.iter().find(|r| &r.0 == id).map(|r| r.1.as_str());
            assert_eq!(other, Some(terms.as_str()), "{id} with as_point={as_point}");
        }
    }
}
