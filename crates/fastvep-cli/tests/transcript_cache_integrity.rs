//! A transcript cache that will not load must stop the run, not quietly change
//! the answer.
//!
//! `save_cache` used to create the destination and stream into it, so a second
//! annotation started while the first was still writing read a prefix: the
//! magic header and zstd frame header are written first, so it passed every
//! check and then failed inside deserialization. The CLI warned, fell back, and
//! with no GFF3 left to rebuild from it annotated every variant as intergenic at
//! exit code 0. On validation/human/vep_example_GRCh38.vcf that produced 173
//! well-formed intergenic calls in place of real consequences.

use std::fs::File;
use std::io::Write;
use std::path::{Path, PathBuf};

use fastvep_cli::pipeline::{run_annotate, AnnotateConfig};
use tempfile::TempDir;

const CDS: u64 = 1_001;

fn write_gff3(dir: &Path) -> PathBuf {
    let path = dir.join("ann.gff3");
    let end = CDS + 11;
    std::fs::write(
        &path,
        format!(
            "1\ttest\tgene\t{CDS}\t{end}\t.\t+\t.\tID=gene:ENSG_A;biotype=protein_coding\n\
             1\ttest\tmRNA\t{CDS}\t{end}\t.\t+\t.\tID=transcript:ENST_A;Parent=gene:ENSG_A;biotype=protein_coding;tag=Ensembl_canonical\n\
             1\ttest\texon\t{CDS}\t{end}\t.\t+\t.\tID=exon:ENSE_A;Parent=transcript:ENST_A;rank=1\n\
             1\ttest\tCDS\t{CDS}\t{end}\t.\t+\t0\tID=CDS:ENSP_A;Parent=transcript:ENST_A\n"
        ),
    )
    .unwrap();
    path
}

fn write_fasta(dir: &Path) -> PathBuf {
    let mut seq = vec![b'N'; 1_100];
    seq[(CDS - 1) as usize..(CDS - 1) as usize + 12].copy_from_slice(b"ATGTGGCGGTAA");
    let path = dir.join("ref.fa");
    let mut f = File::create(&path).unwrap();
    writeln!(f, ">1").unwrap();
    for chunk in seq.chunks(60) {
        f.write_all(chunk).unwrap();
        f.write_all(b"\n").unwrap();
    }
    path
}

fn write_vcf(dir: &Path) -> PathBuf {
    let path = dir.join("in.vcf");
    let mut f = File::create(&path).unwrap();
    writeln!(f, "##fileformat=VCFv4.2").unwrap();
    writeln!(f, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO").unwrap();
    writeln!(f, "1\t{}\t.\tG\tA\t.\t.\t.", CDS + 7).unwrap();
    path
}

fn config(input: &Path, out: &Path) -> AnnotateConfig {
    AnnotateConfig {
        input: input.to_string_lossy().into(),
        output: out.to_string_lossy().into(),
        gff3: vec![],
        fasta: None,
        output_format: "vcf".into(),
        pick: false,
        hgvs: true,
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

#[test]
fn a_truncated_explicit_transcript_cache_stops_the_run() {
    let dir = TempDir::new().unwrap();
    let gff3 = write_gff3(dir.path());
    let fasta = write_fasta(dir.path());
    let vcf = write_vcf(dir.path());

    // Build a real cache the ordinary way, and check the ordinary way works.
    let out_ok = dir.path().join("ok.vcf");
    run_annotate(AnnotateConfig {
        gff3: vec![gff3.to_string_lossy().into()],
        fasta: Some(fasta.to_string_lossy().into()),
        ..config(&vcf, &out_ok)
    })
    .expect("baseline annotation should succeed");
    let sidecar = dir.path().join("ann.gff3.fastvep.cache");
    assert!(sidecar.exists(), "expected a sidecar cache to be written");
    assert!(std::fs::read_to_string(&out_ok)
        .unwrap()
        .contains("missense_variant"));

    // Truncate it the way a concurrent writer or a full disk would.
    let full = std::fs::read(&sidecar).unwrap();
    let truncated = dir.path().join("truncated.cache");
    std::fs::write(&truncated, &full[..full.len() / 2]).unwrap();

    // Named explicitly and unusable: there is nothing to substitute, so the run
    // must fail rather than report the variant as intergenic.
    let out_bad = dir.path().join("bad.vcf");
    let err = run_annotate(AnnotateConfig {
        transcript_cache: Some(truncated.to_string_lossy().into()),
        fasta: Some(fasta.to_string_lossy().into()),
        ..config(&vcf, &out_bad)
    })
    .expect_err("a truncated --transcript-cache must not annotate anyway");
    let msg = err.to_string();
    assert!(
        msg.contains("could not be loaded") && msg.contains("truncated or corrupt"),
        "error should name the cache as the problem, got: {msg}"
    );
    assert!(
        !out_bad.exists()
            || !std::fs::read_to_string(&out_bad)
                .unwrap()
                .contains("intergenic_variant"),
        "a failed cache load must not leave an all-intergenic output behind"
    );
}

#[test]
fn a_run_with_no_usable_transcript_source_stops_instead_of_calling_everything_intergenic() {
    let dir = TempDir::new().unwrap();
    let fasta = write_fasta(dir.path());
    let vcf = write_vcf(dir.path());
    let out = dir.path().join("out.vcf");

    let err = run_annotate(AnnotateConfig {
        fasta: Some(fasta.to_string_lossy().into()),
        ..config(&vcf, &out)
    })
    .expect_err("no --gff3 and no cache must not silently produce intergenic calls");
    let msg = err.to_string();
    assert!(
        msg.contains("No transcript source is usable"),
        "unexpected error: {msg}"
    );
}

#[test]
fn a_sidecar_cache_that_will_not_load_is_rebuilt_from_the_gff3() {
    // The counterpart: when the GFF3 is right there, a corrupt sidecar is a
    // recoverable condition and the run should repair itself rather than fail.
    let dir = TempDir::new().unwrap();
    let gff3 = write_gff3(dir.path());
    let fasta = write_fasta(dir.path());
    let vcf = write_vcf(dir.path());

    let sidecar = dir.path().join("ann.gff3.fastvep.cache");
    std::fs::write(&sidecar, b"FSTVEP02 and then garbage that will not inflate").unwrap();
    // Make it look fresh against its source, so the run really does try to load
    // it rather than skipping it as stale.
    let future = std::time::SystemTime::now() + std::time::Duration::from_secs(3600);
    File::options()
        .write(true)
        .open(&sidecar)
        .unwrap()
        .set_times(std::fs::FileTimes::new().set_modified(future))
        .unwrap();

    let out = dir.path().join("out.vcf");
    run_annotate(AnnotateConfig {
        gff3: vec![gff3.to_string_lossy().into()],
        fasta: Some(fasta.to_string_lossy().into()),
        ..config(&vcf, &out)
    })
    .expect("a corrupt sidecar should be rebuilt from the GFF3, not fatal");
    assert!(std::fs::read_to_string(&out)
        .unwrap()
        .contains("missense_variant"));
}
