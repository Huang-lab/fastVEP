//! End-to-end regression for issue #81: an in-frame insertion must not be
//! reported as a substitution in HGVSp.
//!
//! The CLI pipeline carries its own copy of the HGVSp routing block, separate
//! from the one in `fastvep-annotate`. Both had the same bug and both need
//! their own test, because nothing structurally keeps the two copies in step.
//! This one drives the real `run_annotate` entry point over a GFF3 + FASTA +
//! VCF on disk, so it fails if the routing regresses anywhere between parsing
//! the variant and writing the CSQ field.

use fastvep_cli::pipeline::{run_annotate, AnnotateConfig};
use std::fs::File;
use std::io::Write;
use std::path::Path;

// Single-exon coding transcript on contig `1`, CDS = the whole 12 bp exon.
// `ID=CDS:ENSP_INS1` supplies the protein ID that HGVSp is prefixed with.
const GFF3: &str = "##gff-version 3\n\
1\ttest\tgene\t1\t12\t.\t+\t.\tID=gene:ENSG_INS1;Name=INSTEST;gene_name=INSTEST;biotype=protein_coding\n\
1\ttest\tmRNA\t1\t12\t.\t+\t.\tID=transcript:ENST_INS1;Parent=gene:ENSG_INS1;biotype=protein_coding\n\
1\ttest\texon\t1\t12\t.\t+\t.\tID=exon:ENSE_INS1;Parent=transcript:ENST_INS1;rank=1\n\
1\ttest\tCDS\t1\t12\t.\t+\t0\tID=CDS:ENSP_INS1;Parent=transcript:ENST_INS1\n";

// ATG TGG CGG TAA = Met Trp Arg Ter.
const FASTA: &str = ">1\nATGTGGCGGTAA\n";

// Insert CGG after the Trp codon (bases 4-6), VCF-anchored at base 6:
// Met Trp Arg Ter -> Met Trp Arg Arg Ter, an in-frame duplication of Arg3.
const VCF: &str = "##fileformat=VCFv4.2\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
1\t6\t.\tG\tGCGG\t.\tPASS\t.\n";

fn write(path: &Path, contents: &str) {
    File::create(path)
        .unwrap()
        .write_all(contents.as_bytes())
        .unwrap();
}

#[test]
fn cli_inframe_insertion_hgvsp_is_a_delins_not_a_substitution() {
    let tmp = tempfile::tempdir().unwrap();
    let gff3 = tmp.path().join("test.gff3");
    let fasta = tmp.path().join("test.fa");
    let vcf_in = tmp.path().join("in.vcf");
    let vcf_out = tmp.path().join("out.vcf");
    write(&gff3, GFF3);
    write(&fasta, FASTA);
    write(&vcf_in, VCF);

    run_annotate(AnnotateConfig {
        input: vcf_in.to_string_lossy().into(),
        output: vcf_out.to_string_lossy().into(),
        gff3: vec![gff3.to_string_lossy().into()],
        fasta: Some(fasta.to_string_lossy().into()),
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
    })
    .expect("annotation should succeed");

    let out = std::fs::read_to_string(&vcf_out).unwrap();

    // Guard the premise: if this stops being classified as an in-frame
    // insertion, the HGVSp assertions below would pass for the wrong reason.
    assert!(
        out.contains("inframe_insertion"),
        "expected an inframe_insertion consequence, got:\n{}",
        out
    );

    // The bug: `hgvsp()` compared only the first residue of `R` vs `RR`,
    // found them equal, and emitted a synonymous call for a variant that
    // lengthens the protein.
    assert!(
        !out.contains("p.Arg3="),
        "in-frame insertion rendered as synonymous (issue #81):\n{}",
        out
    );
    assert!(
        out.contains("p.Arg3delinsArgArg"),
        "expected an un-normalised delins for the in-frame insertion, got:\n{}",
        out
    );
}
