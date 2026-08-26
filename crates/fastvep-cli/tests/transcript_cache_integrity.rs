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
        msg.contains("cannot be used") && msg.contains("truncated or corrupt"),
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

/// Write a well-formed cache in the pre-#90 FSTVEP02 format, holding `n`
/// transcripts. Deliberately valid: the whole difficulty is that a poisoned
/// region-restricted cache and a whole-file one are the same bytes with a
/// different transcript count, and neither the count nor the mtime tells them
/// apart.
fn write_v2_cache(path: &Path, transcripts: &[fastvep_genome::Transcript]) {
    use std::io::{BufWriter, Write};
    let file = File::create(path).unwrap();
    let mut zst = zstd::Encoder::new(BufWriter::new(file), 1).unwrap();
    zst.write_all(b"FSTVEP02").unwrap();
    bincode::serialize_into(&mut zst, transcripts).unwrap();
    zst.finish().unwrap().flush().unwrap();
}

fn make_fresh(path: &Path) {
    // An hour into the future, so `cache_is_fresh` says yes against its source
    // and the run really attempts the load rather than skipping it as stale.
    let future = std::time::SystemTime::now() + std::time::Duration::from_secs(3600);
    File::options()
        .write(true)
        .open(path)
        .unwrap()
        .set_times(std::fs::FileTimes::new().set_modified(future))
        .unwrap();
}

#[test]
fn a_pre_90_sidecar_cache_is_rebuilt_rather_than_trusted() {
    // #90 stopped persisting region-restricted transcript sets, but every cache
    // already on disk from a build between the tabix read path (2026-06-10) and
    // that fix may be one. This is the regression that gate leaves open: a fresh
    // FSTVEP02 sidecar, holding a transcript set that is *valid* and *wrong for
    // this input*, next to the GFF3 it was keyed on.
    let dir = TempDir::new().unwrap();
    let gff3 = write_gff3(dir.path());
    let fasta = write_fasta(dir.path());
    let vcf = write_vcf(dir.path());

    // An empty transcript set stands in for "someone else's neighbourhood":
    // deserializes cleanly, and every variant comes back intergenic.
    let sidecar = dir.path().join("ann.gff3.fastvep.cache");
    write_v2_cache(&sidecar, &[]);
    make_fresh(&sidecar);

    let out = dir.path().join("out.vcf");
    run_annotate(AnnotateConfig {
        gff3: vec![gff3.to_string_lossy().into()],
        fasta: Some(fasta.to_string_lossy().into()),
        ..config(&vcf, &out)
    })
    .expect("a pre-#90 sidecar should be rebuilt from the GFF3, not fatal");

    let annotated = std::fs::read_to_string(&out).unwrap();
    assert!(
        annotated.contains("missense_variant"),
        "the run trusted a pre-#90 cache and lost the annotation; got:\n{annotated}"
    );
    assert!(
        !annotated.contains("intergenic_variant"),
        "expected the rebuilt transcript set, not an all-intergenic answer"
    );
    // And the rebuild republished in the current format, so the next run is fast.
    let republished = std::fs::read(&sidecar).unwrap();
    let mut magic = [0u8; 8];
    {
        use std::io::Read;
        zstd::Decoder::new(&republished[..])
            .unwrap()
            .read_exact(&mut magic)
            .unwrap();
    }
    assert_eq!(
        &magic, b"FSTVEP04",
        "rebuild should publish the current format"
    );
}

/// Write a well-formed cache in the pre-#98 FSTVEP03 format. Whole-file and
/// perfectly readable - the only thing wrong with it is that the parser which
/// filled it did not recognise `ncRNA_gene`, so its non-coding genes have no
/// symbol.
fn write_v3_cache(path: &Path, transcripts: &[fastvep_genome::Transcript]) {
    use std::io::{BufWriter, Write};
    let file = File::create(path).unwrap();
    let mut zst = zstd::Encoder::new(BufWriter::new(file), 1).unwrap();
    zst.write_all(b"FSTVEP03").unwrap();
    bincode::serialize_into(&mut zst, transcripts).unwrap();
    zst.finish().unwrap().flush().unwrap();
}

#[test]
fn a_pre_98_sidecar_cache_is_rebuilt_rather_than_trusted() {
    // The gap #98 would otherwise leave open. Its reporter runs one command; the
    // sidecar next to their GFF3 is newer than the GFF3, so upgrading to the
    // fixed parser changes nothing unless the stale cache is refused - the
    // parser simply never runs again. Same shape as the pre-#90 case above, one
    // format later.
    let dir = TempDir::new().unwrap();
    let gff3 = write_gff3(dir.path());
    let fasta = write_fasta(dir.path());
    let vcf = write_vcf(dir.path());

    let sidecar = dir.path().join("ann.gff3.fastvep.cache");
    write_v3_cache(&sidecar, &[]);
    make_fresh(&sidecar);

    let out = dir.path().join("out.vcf");
    run_annotate(AnnotateConfig {
        gff3: vec![gff3.to_string_lossy().into()],
        fasta: Some(fasta.to_string_lossy().into()),
        ..config(&vcf, &out)
    })
    .expect("a pre-#98 sidecar should be rebuilt from the GFF3, not fatal");

    let annotated = std::fs::read_to_string(&out).unwrap();
    assert!(
        annotated.contains("missense_variant"),
        "the run trusted a pre-#98 cache and lost the annotation; got:\n{annotated}"
    );
}

#[test]
fn a_pre_90_explicit_transcript_cache_stops_the_run() {
    // Named explicitly there is nothing to fall back to, so the only safe answer
    // is to refuse - and to say why, because "rebuild it" is not obvious from
    // "invalid cache".
    let dir = TempDir::new().unwrap();
    let fasta = write_fasta(dir.path());
    let vcf = write_vcf(dir.path());

    let stale = dir.path().join("old.cache");
    write_v2_cache(&stale, &[]);

    let out = dir.path().join("out.vcf");
    let err = run_annotate(AnnotateConfig {
        transcript_cache: Some(stale.to_string_lossy().into()),
        fasta: Some(fasta.to_string_lossy().into()),
        ..config(&vcf, &out)
    })
    .expect_err("a pre-#90 --transcript-cache must not annotate anyway");

    let msg = format!("{err:#}");
    assert!(
        msg.contains("FSTVEP02") || msg.contains("pre-0.3.1"),
        "the error should name the format, got: {msg}"
    );
    assert!(
        msg.contains("fastvep cache"),
        "the error should say how to recover, got: {msg}"
    );
    // One diagnosis, not two. The file is intact - sending the user to look for
    // a truncated write or a full disk would be a different problem than the
    // one they have, and the generic wording used to be appended to this one.
    assert!(
        !msg.contains("truncated or corrupt"),
        "a pre-#90 cache is intact; the error must not also blame corruption, got: {msg}"
    );
    assert!(
        !out.exists()
            || !std::fs::read_to_string(&out)
                .unwrap()
                .contains("intergenic_variant"),
        "a rejected cache must not leave an all-intergenic output behind"
    );
}
