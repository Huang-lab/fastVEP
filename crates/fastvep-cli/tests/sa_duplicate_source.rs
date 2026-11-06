//! Two databases that share a JSON key and both know a variant.
//!
//! An allele's annotation holds one value per key, so the pipeline keeps one and
//! drops the other. Which is decided by file name alone, because the databases
//! are read in path order, and for `gnomad` the loser can be the exomes frequency
//! behind every coding variant of a chromosome, which BA1, BS1, BS2 and PM2 all
//! read. Loading the exomes and the genomes of gnomAD side by side is the obvious
//! way to ask for both, and it is not a thing the classifier can honour: the
//! counts of two cohorts do not combine into a frequency without the per-group
//! counts the database does not hold. gnomAD's v4.1 joint release is the
//! combination, and the builder reads it.
//!
//! What the run owes the user is to say so. This lives in its own test binary
//! because the collector it reads is process-wide.

use fastvep_cli::pipeline::{run_annotate, run_sa_build, AnnotateConfig, PickFlags};
use serde_json::Value;
use std::fs;

const GFF3: &str = "##gff-version 3\n\
1\ttest\tgene\t1\t380\t.\t+\t.\tID=gene:ENSG_DUP;Name=DUPTEST;gene_name=DUPTEST;biotype=protein_coding\n\
1\ttest\tmRNA\t1\t380\t.\t+\t.\tID=transcript:ENST_DUP;Parent=gene:ENSG_DUP;biotype=protein_coding\n\
1\ttest\texon\t1\t380\t.\t+\t.\tID=exon:E1;Parent=transcript:ENST_DUP;rank=1\n\
1\ttest\tCDS\t1\t378\t.\t+\t0\tID=CDS:ENSP_DUP;Parent=transcript:ENST_DUP\n";

#[test]
fn two_gnomad_databases_for_one_chromosome_are_reported_and_the_last_by_path_is_kept() {
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();
    let sa = dir.join("sa");
    fs::create_dir(&sa).unwrap();

    fs::write(dir.join("genes.gff3"), GFF3).unwrap();
    fs::write(
        dir.join("genome.fa"),
        format!(">1\n{}\n", "ACGT".repeat(95)),
    )
    .unwrap();
    // The reference is `ACGT` repeated. Take the base at position 100 from it, so
    // that the record's REF agrees with the genome and a mismatch cannot muddy
    // what the test is about.
    let reference: String = "ACGT".repeat(95);
    let ref_at_100 = &reference[99..100];
    fs::write(
        dir.join("in.vcf"),
        format!(
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n1\t100\t.\t{ref_at_100}\tA\t.\tPASS\t.\n"
        ),
    )
    .unwrap();
    // `gnomad_a` sorts before `gnomad_b`, so `gnomad_b` is the one kept.
    for (name, af) in [("gnomad_a", 0.0001_f64), ("gnomad_b", 0.2_f64)] {
        let vcf = dir.join(format!("{name}.vcf"));
        fs::write(
            &vcf,
            format!(
                "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
                 1\t100\t.\t{ref_at_100}\tA\t.\tPASS\tAF={af};AC=10;AN=1000000;nhomalt=0\n"
            ),
        )
        .unwrap();
        run_sa_build(
            "gnomad",
            vcf.to_str().unwrap(),
            sa.join(name).to_str().unwrap(),
            "GRCh38",
            None,
            &[],
            false,
        )
        .unwrap();
    }

    let out = dir.join("out.json");
    run_annotate(AnnotateConfig {
        input: dir.join("in.vcf").to_string_lossy().into(),
        output: out.to_string_lossy().into(),
        gff3: vec![dir.join("genes.gff3").to_string_lossy().into()],
        fasta: Some(dir.join("genome.fa").to_string_lossy().into()),
        output_format: "json".into(),
        pick: PickFlags::default(),
        hgvs: false,
        distance: 5000,
        cache_dir: None,
        transcript_cache: None,
        sa_dir: Some(sa.to_string_lossy().into()),
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
        columns: Default::default(),
        exclude: Default::default(),
    })
    .expect("annotation should succeed");

    // The kept answer is the last by path: documented by the warning, pinned here.
    let variants: Vec<Value> = serde_json::from_str(&fs::read_to_string(&out).unwrap()).unwrap();
    let kept = variants[0]["transcript_consequences"][0]["gnomad"]["allAf"]
        .as_f64()
        .expect("a gnomad annotation");
    assert!((kept - 0.2).abs() < 1e-6, "allAf was {kept}");

    // And the run says that it chose.
    let lines = fastvep_annotate::sa_duplicate_answers().report_lines();
    assert!(
        lines
            .iter()
            .any(|l| l.contains("`gnomad`") && l.contains("1:100")),
        "{lines:?}"
    );
}
