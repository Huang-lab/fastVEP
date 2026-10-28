//! The row of a transcript that no consequence term applies to is not `intergenic_variant`.
//!
//! A delins that reaches the last base of one exon, the one-base intron after it and the first base
//! of the next exon (`chr1:4556 AAG>CGA` over a reverse-strand transcript with a coding exon of six
//! bases) is neither inside an exon nor inside an intron, and none of Ensembl's predicates holds
//! for it. Ensembl VEP 116.3 falls back on its default overlap consequence
//! (`$cons = [$DEFAULT_OVERLAP_CONSEQUENCE] unless @$cons`) and writes `intergenic_variant` on the
//! row of the transcript, with `c.40_41delinsTCG` and `p.Leu14SerfsTer2` beside it (found by a
//! differential fuzz of the stack, run through the pinned 116.3 image). That says the change is
//! between genes where it is in one, and a filter that drops intergenic rows hides it, so this
//! writes a term of the transcript (`coding_transcript_variant`), and `docs/VEP_DIVERGENCE.md`
//! lists the difference. No transcript row of the 2,026,951 of VEP's own run over the Ensembl
//! cache's edge set, nor of the 212,240 over the RefSeq cache's, has the term.
//!
//! Run as a program, because the answer is the `Consequence` column.

use std::io::Write;
use std::path::Path;
use std::process::Command;

use tempfile::TempDir;

const GFF3: &str = "chr1\ttest\tgene\t4478\t4617\t.\t-\t.\tID=gene:ENSG00000000005;Name=G4;biotype=protein_coding;gene_id=ENSG00000000005;version=1\n\
chr1\ttest\tmRNA\t4478\t4617\t.\t-\t.\tID=transcript:ENST00000000005;Parent=gene:ENSG00000000005;Name=G4-001;biotype=protein_coding;tag=Ensembl_canonical;transcript_id=ENST00000000005;version=1\n\
chr1\ttest\texon\t4558\t4617\t.\t-\t.\tParent=transcript:ENST00000000005;Name=ENSE00000000015;exon_id=ENSE00000000015;rank=1;version=1\n\
chr1\ttest\texon\t4551\t4556\t.\t-\t.\tParent=transcript:ENST00000000005;Name=ENSE00000000016;exon_id=ENSE00000000016;rank=2;version=1\n\
chr1\ttest\texon\t4478\t4537\t.\t-\t.\tParent=transcript:ENST00000000005;Name=ENSE00000000017;exon_id=ENSE00000000017;rank=3;version=1\n\
chr1\ttest\tCDS\t4558\t4597\t.\t-\t0\tID=CDS:ENSP00000000005;Parent=transcript:ENST00000000005;protein_id=ENSP00000000005;version=1\n\
chr1\ttest\tCDS\t4551\t4556\t.\t-\t2\tID=CDS:ENSP00000000005;Parent=transcript:ENST00000000005;protein_id=ENSP00000000005;version=1\n\
chr1\ttest\tCDS\t4509\t4537\t.\t-\t2\tID=CDS:ENSP00000000005;Parent=transcript:ENST00000000005;protein_id=ENSP00000000005;version=1\n";

/// Bases 4,401 to 4,700 of the genome of the fuzz that found this, which is the transcript and
/// what is beside it: whether a splice term fires for the change depends on them. The 4,400 before
/// them are filler.
const BASES: &str = "TTAAAGTTGCAGCACCCATTAAGCCAGAGCTGATTGGCTCAGGCTTCGAAGAGAGCGGAATACCTAATACAACTTCTCTTCTAGGGGAACGGACAAGGGTGTGCCACATCAAAACCCGAAGAAGTTGGCTTTTTTCGGCACACTGCGCCCACTTAAAGGATGTCTCCAACTCCTGATGTCGTTGGGAAACGCCACATTCTCCTAAACTGCACAGGTGGAACCAGCCAGCGTCGGAGGAGTCGGGCCCTAATGGCACCTGTTGGGATGACGGTTACTTATGTCCGTGGTGTGAGCGGGCAG";

fn genome() -> String {
    format!("{}{BASES}", "A".repeat(4400))
}

/// `(variant id, Feature, Consequence)` of every row.
fn run(dir: &Path) -> Vec<(String, String, String)> {
    std::fs::write(dir.join("genes.gff3"), GFF3).unwrap();
    let fasta = dir.join("genome.fa");
    std::fs::write(&fasta, format!(">chr1\n{}\n", genome())).unwrap();
    let vcf = dir.join("in.vcf");
    let mut f = std::fs::File::create(&vcf).unwrap();
    writeln!(f, "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=4700>").unwrap();
    writeln!(f, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO").unwrap();
    let genome = genome();
    // the delins over the intron, with the genome's own bases as its reference
    let reference = &genome[4555..4558];
    writeln!(f, "chr1\t4556\tspans_the_intron\t{reference}\tCGA\t.\t.\t.").unwrap();
    // an ordinary change in the coding exon beside it, for the harness
    let b = &genome[4569..4570];
    let alt = if b == "A" { "G" } else { "A" };
    writeln!(f, "chr1\t4570\tinside_the_exon\t{b}\t{alt}\t.\t.\t.").unwrap();
    drop(f);
    let out = dir.join("out.tab");
    let output = Command::new(env!("CARGO_BIN_EXE_fastvep"))
        .args(["annotate", "-i"])
        .arg(&vcf)
        .arg("-o")
        .arg(&out)
        .arg("--gff3")
        .arg(dir.join("genes.gff3"))
        .arg("--fasta")
        .arg(&fasta)
        .args(["--output-format", "tab", "--hgvs", "--no-progress"])
        .output()
        .expect("the fastvep binary runs");
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let tab = std::fs::read_to_string(&out).unwrap();
    let header: Vec<&str> = tab
        .lines()
        .find(|l| l.starts_with("#Uploaded_variation"))
        .expect("a header")
        .split('\t')
        .collect();
    let col = |name: &str| header.iter().position(|c| *c == name).expect(name);
    let (id, feature, consequence) = (
        col("#Uploaded_variation"),
        col("Feature"),
        col("Consequence"),
    );
    tab.lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            (
                f[id].to_string(),
                f[feature].to_string(),
                f[consequence].to_string(),
            )
        })
        .collect()
}

#[test]
fn a_change_no_term_applies_to_is_a_transcript_variant_and_not_intergenic() {
    let dir = TempDir::new().unwrap();
    let rows = run(dir.path());
    let of = |id: &str| -> Vec<&str> {
        rows.iter()
            .filter(|(i, f, _)| i == id && f == "ENST00000000005")
            .map(|(_, _, c)| c.as_str())
            .collect()
    };
    assert_eq!(
        of("spans_the_intron"),
        ["coding_transcript_variant"],
        "VEP 116.3 writes intergenic_variant here, on the row of the transcript"
    );
    // a change inside the exon is an ordinary one, so the row above is not every row
    let inside = of("inside_the_exon");
    assert!(
        !inside.is_empty() && inside.iter().all(|c| !c.contains("intergenic")),
        "{inside:?}"
    );
}
