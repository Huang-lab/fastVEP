//! End-to-end regression for issue #96: on the reverse strand the affected
//! residue span is anchored at its *last* residue, and when the reference run is
//! periodic the wrong end can look just as plausible as the right one.
//!
//! `protein_start` comes from the genomic left edge, which for a shrinking
//! change on a right-to-left transcript is the end of the affected range. #93
//! made `hgvsp_inframe_indel` try the other end when the peptide does not
//! corroborate the first, which fixed the cases where the first end is not a
//! match at all. It left the case where *both* ends match: `find()` takes the
//! first, and on the reverse strand the first is the wrong one.
//!
//! The gene below is built so that both ends match. The peptide is `MEGEGEA`
//! and the deletion removes `EGE` at residues 2-4 - a reference that also sits
//! at 4-6. Naming 4-6 is not another spelling of the same event: deleting 2-4
//! leaves `MGEA`, deleting 4-6 leaves `MEGA`, so the 3'-rule cannot merge them
//! and the description names residues this variant does not touch.
//!
//! Driven through `run_annotate` rather than the unit under test, because the
//! anchor convention is a property of the coordinates the pipeline hands over -
//! a unit test that passes `protein_start` by hand would be asserting my own
//! reading of the convention rather than the pipeline's.

use fastvep_cli::pipeline::{run_annotate, AnnotateConfig};
use std::fs::File;
use std::io::Write;
use std::path::Path;

/// Single-exon coding transcript on contig `1`, CDS = the whole 24 bp exon.
/// `ID=CDS:ENSP_REV1` supplies the HGVSp prefix.
fn gff3(strand: char) -> String {
    format!(
        "##gff-version 3\n\
         1\ttest\tgene\t1\t24\t.\t{strand}\t.\tID=gene:ENSG_REV1;Name=REVTEST;gene_name=REVTEST;biotype=protein_coding\n\
         1\ttest\tmRNA\t1\t24\t.\t{strand}\t.\tID=transcript:ENST_REV1;Parent=gene:ENSG_REV1;biotype=protein_coding;tag=Ensembl_canonical\n\
         1\ttest\texon\t1\t24\t.\t{strand}\t.\tID=exon:ENSE_REV1;Parent=transcript:ENST_REV1;rank=1\n\
         1\ttest\tCDS\t1\t24\t.\t{strand}\t0\tID=CDS:ENSP_REV1;Parent=transcript:ENST_REV1\n"
    )
}

// ATG GAA GGT GAA GGT GAA GCT TAA translates to Met Glu Gly Glu Gly Glu Ala Ter,
// so `EGE` sits at residues 2-4 and again at 4-6. The reverse-strand contig
// carries the reverse complement, so both strands describe the same protein and
// the same variant, and the only difference is which end of the affected span
// the pipeline hands to HGVSp.
const CDS_FORWARD: &str = ">1\nATGGAAGGTGAAGGTGAAGCTTAA\n";
const CDS_REVERSE: &str = ">1\nTTAAGCTTCACCTTCACCTTCCAT\n";

// Delete the codons for residues 2, 3 and 4 - CDS 4-12 either way, which is
// genomic 4-12 on the forward strand and genomic 13-21 on the reverse. Both
// VCF-anchored on the base before the deletion, as a caller's VCF would be.
const VCF_FORWARD: &str = "##fileformat=VCFv4.2\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
1\t3\t.\tGGAAGGTGAA\tG\t.\tPASS\t.\n";
const VCF_REVERSE: &str = "##fileformat=VCFv4.2\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
1\t12\t.\tCTTCACCTTC\tC\t.\tPASS\t.\n";

fn write(path: &Path, contents: &str) {
    File::create(path)
        .unwrap()
        .write_all(contents.as_bytes())
        .unwrap();
}

fn annotate(dir: &Path, strand: char) -> String {
    let gff3_path = dir.join("rev.gff3");
    let fasta = dir.join("rev.fa");
    let vcf_in = dir.join("in.vcf");
    let vcf_out = dir.join("out.vcf");
    let (contig, vcf) = if strand == '-' {
        (CDS_REVERSE, VCF_REVERSE)
    } else {
        (CDS_FORWARD, VCF_FORWARD)
    };
    write(&gff3_path, &gff3(strand));
    write(&fasta, contig);
    write(&vcf_in, vcf);

    run_annotate(AnnotateConfig {
        input: vcf_in.to_string_lossy().into(),
        output: vcf_out.to_string_lossy().into(),
        gff3: vec![gff3_path.to_string_lossy().into()],
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

    std::fs::read_to_string(&vcf_out).unwrap()
}

#[test]
fn cli_reverse_strand_deletion_names_the_residues_it_deletes() {
    let tmp = tempfile::tempdir().unwrap();
    let out = annotate(tmp.path(), '-');

    // Guard the premise twice over: a different consequence or a different
    // protein position would make the HGVSp assertion pass or fail for reasons
    // that have nothing to do with the anchor.
    assert!(
        out.contains("inframe_deletion"),
        "expected an inframe_deletion consequence, got:\n{out}"
    );
    assert!(
        out.contains("|2-4|"),
        "expected Protein_position 2-4 for the deleted residues, got:\n{out}"
    );

    assert!(
        out.contains("p.Glu2_Glu4del"),
        "expected the deletion of residues 2-4, got:\n{out}"
    );
    // The specific wrong answer, spelled out: `EGE` also sits at 4-6, and
    // taking the first corroborated anchor picks it.
    assert!(
        !out.contains("p.Glu4_Glu6del"),
        "HGVSp named residues 4-6, which this variant does not touch (issue #96):\n{out}"
    );
}

#[test]
fn cli_forward_strand_deletion_of_the_same_residues_is_unchanged() {
    // The mirror image, on a transcript that runs the other way. The protein,
    // the deleted residues and the expected description are identical; only the
    // end of the span the pipeline passes differs. Ordering the anchors by
    // strand is only correct if it leaves this case exactly as it was, so this
    // is the half of #96 that guards against overcorrecting.
    let tmp = tempfile::tempdir().unwrap();
    let out = annotate(tmp.path(), '+');

    assert!(
        out.contains("inframe_deletion"),
        "expected an inframe_deletion consequence, got:\n{out}"
    );
    assert!(
        out.contains("|2-4|"),
        "expected Protein_position 2-4 for the deleted residues, got:\n{out}"
    );
    assert!(
        out.contains("p.Glu2_Glu4del"),
        "expected the deletion of residues 2-4, got:\n{out}"
    );
}
