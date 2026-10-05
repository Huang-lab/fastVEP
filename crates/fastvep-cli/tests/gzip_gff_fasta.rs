//! Gzipped GFF3 and FASTA take the same path as a plain file.
//!
//! The VCF opener already decided by magic bytes, with the suffix as a
//! fallback. A GFF3 or FASTA that was only recognised by its `.gz` name
//! parsed as text when somebody renamed it, and a `.fa.gz` was never
//! decompressed at all.

use std::io::Write;
use std::path::{Path, PathBuf};

use fastvep_annotate::AnnotationContext;
use fastvep_cli::pipeline::{run_annotate, AnnotateConfig, PickFlags};
use flate2::write::GzEncoder;
use flate2::Compression;
use tempfile::TempDir;

const GFF3: &str = "\
1\ttest\tgene\t1\t12\t.\t+\t.\tID=gene:ENSG;Name=DEMO;biotype=protein_coding\n\
1\ttest\tmRNA\t1\t12\t.\t+\t.\tID=transcript:ENST;Parent=gene:ENSG;biotype=protein_coding\n\
1\ttest\texon\t1\t12\t.\t+\t.\tID=exon:ENSE;Parent=transcript:ENST;rank=1\n\
1\ttest\tCDS\t1\t12\t.\t+\t0\tID=CDS:ENSP;Parent=transcript:ENST\n";

const FASTA: &str = ">1\nATGTGGCGGTAA\n";

const VCF: &str = "\
##fileformat=VCFv4.2\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
1\t8\t.\tG\tA\t.\t.\t.\n";

fn gzip(text: &str) -> Vec<u8> {
    let mut enc = GzEncoder::new(Vec::new(), Compression::default());
    enc.write_all(text.as_bytes()).unwrap();
    enc.finish().unwrap()
}

fn write(dir: &Path, name: &str, bytes: &[u8]) -> PathBuf {
    let path = dir.join(name);
    std::fs::write(&path, bytes).unwrap();
    path
}

fn config(input: &Path, out: &Path) -> AnnotateConfig {
    AnnotateConfig {
        input: input.to_string_lossy().into(),
        output: out.to_string_lossy().into(),
        gff3: vec![],
        fasta: None,
        output_format: "vcf".into(),
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

fn variant_lines(annotated: &str) -> Vec<&str> {
    annotated
        .lines()
        .filter(|l| !l.starts_with('#') && !l.is_empty())
        .collect()
}

fn annotate(dir: &Path, gff3: &Path, fasta: &Path) -> String {
    let vcf = write(dir, "in.vcf", VCF.as_bytes());
    let out = dir.join("out.vcf");
    run_annotate(AnnotateConfig {
        // A fixed label, so the SOURCE column does not change with the filename.
        gff3: vec![format!("Demo={}", gff3.display())],
        fasta: Some(fasta.to_string_lossy().into()),
        ..config(&vcf, &out)
    })
    .expect("annotation should succeed");
    std::fs::read_to_string(&out).unwrap()
}

#[test]
fn gzipped_gff_and_fasta_match_the_plain_run() {
    let plain_dir = TempDir::new().unwrap();
    let plain = annotate(
        plain_dir.path(),
        &write(plain_dir.path(), "ann.gff3", GFF3.as_bytes()),
        &write(plain_dir.path(), "ref.fa", FASTA.as_bytes()),
    );

    let gz_dir = TempDir::new().unwrap();
    let named = annotate(
        gz_dir.path(),
        &write(gz_dir.path(), "ann.gff3.gz", &gzip(GFF3)),
        &write(gz_dir.path(), "ref.fa.gz", &gzip(FASTA)),
    );

    // The suffix is a fallback. These names claim to be plain text.
    let magic_dir = TempDir::new().unwrap();
    let magic = annotate(
        magic_dir.path(),
        &write(magic_dir.path(), "ann.gff3", &gzip(GFF3)),
        &write(magic_dir.path(), "ref.fa", &gzip(FASTA)),
    );

    assert_eq!(variant_lines(&named), variant_lines(&plain));
    assert_eq!(variant_lines(&magic), variant_lines(&plain));
    assert!(
        variant_lines(&plain).iter().any(|l| l.contains("missense")),
        "the plain run should have called the coding change, got:\n{plain}"
    );
}

#[test]
fn the_web_context_reads_a_gzipped_gff_and_fasta_whatever_they_are_named() {
    let dir = TempDir::new().unwrap();
    let gff3 = write(dir.path(), "ann.gff3", &gzip(GFF3));
    let fasta = write(dir.path(), "ref.fa", &gzip(FASTA));

    let ctx = AnnotationContext::new(
        Some(gff3.to_str().unwrap()),
        Some(fasta.to_str().unwrap()),
        None,
        5000,
    )
    .expect("gzipped inputs should load");

    assert_eq!(ctx.transcript_count(), 1);
    let seq = ctx
        .seq_provider
        .as_ref()
        .expect("fasta loaded")
        .fetch_sequence("1", 1, 4)
        .unwrap();
    assert_eq!(seq, b"ATGT");
}
