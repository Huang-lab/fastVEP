//! `HGVSg` of a change at the very edge of a contig does not name a base the contig does not have.
//!
//! An insertion before the first base was written `g.0dup` or `g.0_1insT` and one after the last
//! `g.100_101insC`, on every row of the record, transcript rows and upstream and downstream rows
//! alike, where Ensembl VEP 116.3 writes nothing: its `hgvs_genomic` returns nothing for a
//! variation feature that is not inside the slice it is on. A record of several alts passes that
//! check whatever its alleles are, so VEP writes the position that does not exist for one of
//! those too, and this does not.
//!
//! Run as a program, because the answer is the `HGVSg` column. The same variants without a FASTA:
//! the first base of the contig is known and the last is not.

use std::io::Write;
use std::path::Path;
use std::process::Command;

use tempfile::TempDir;

const CONTIG: &str = "chr10";
/// 100 bases: a `G` at the start and `TT` at the end, with a gene in the middle.
const GENOME: &str = "GCTAAGACAATTACATAACAATGGGGGAATTGTCGTCGAGACACAACAATTCGGCTCAATATGCGTCGCATTATCTCGAGGACCGACCGCAACAAGTATT";

fn variants() -> Vec<(&'static str, u64, &'static str, &'static str)> {
    vec![
        ("ins_before_1", 1, "G", "TG"),
        ("snv_1", 1, "G", "A"),
        ("ins_after_last", 100, "T", "TC"),
        ("dup_of_last", 100, "T", "TT"),
        ("del_last", 99, "TT", "T"),
        ("snv_last", 100, "T", "C"),
        ("ref_off_the_end", 99, "TTA", "T"),
        ("multi_ins_before_1", 1, "G", "TG,A"),
        ("multi_ins_after_last", 100, "T", "TC,A"),
        ("ins_inside_near_last", 99, "T", "TC"),
    ]
}

/// `(variant id, Allele, HGVSg)` of every row, one per allele and transcript.
fn run(dir: &Path, with_fasta: bool) -> Vec<(String, String, String)> {
    assert_eq!(GENOME.len(), 100);
    assert!(GENOME.starts_with('G') && GENOME.ends_with('T'));
    let gff3 = dir.join("genes.gff3");
    std::fs::write(
        &gff3,
        format!(
            "{CONTIG}\ttest\tgene\t30\t70\t.\t+\t.\tID=gene:G1;Name=G1;biotype=lncRNA\n\
             {CONTIG}\ttest\tlnc_RNA\t30\t70\t.\t+\t.\tID=transcript:T1;Parent=gene:G1;biotype=lncRNA\n\
             {CONTIG}\ttest\texon\t30\t70\t.\t+\t.\tID=exon:E1;Parent=transcript:T1;rank=1\n"
        ),
    )
    .unwrap();
    let fasta = dir.join("genome.fa");
    std::fs::write(&fasta, format!(">{CONTIG}\n{GENOME}\n")).unwrap();
    let vcf = dir.join("in.vcf");
    let mut f = std::fs::File::create(&vcf).unwrap();
    writeln!(f, "##fileformat=VCFv4.2").unwrap();
    writeln!(f, "##contig=<ID={CONTIG},length=100>").unwrap();
    writeln!(f, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO").unwrap();
    for (id, pos, r, a) in variants() {
        writeln!(f, "{CONTIG}\t{pos}\t{id}\t{r}\t{a}\t.\t.\t.").unwrap();
    }
    drop(f);
    let out = dir.join(if with_fasta {
        "with.tab"
    } else {
        "without.tab"
    });
    let mut command = Command::new(env!("CARGO_BIN_EXE_fastvep"));
    command
        .args(["annotate", "-i"])
        .arg(&vcf)
        .arg("-o")
        .arg(&out)
        .arg("--gff3")
        .arg(&gff3)
        .args(["--output-format", "tab", "--hgvsg", "--no-progress"]);
    if with_fasta {
        command.arg("--fasta").arg(&fasta);
    }
    let output = command.output().expect("the fastvep binary runs");
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
    let (id, allele, hgvsg) = (col("#Uploaded_variation"), col("Allele"), col("HGVSg"));
    tab.lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            (
                f[id].to_string(),
                f[allele].to_string(),
                f[hgvsg].to_string(),
            )
        })
        .collect()
}

/// What the rows of one variant write, once each: the description is the record's, so it is the
/// same on every row of it.
fn hgvsg_of(rows: &[(String, String, String)], id: &str, allele: &str) -> Vec<String> {
    let mut seen: Vec<String> = rows
        .iter()
        .filter(|(i, a, _)| i == id && a == allele)
        .map(|(_, _, h)| h.clone())
        .collect();
    seen.sort();
    seen.dedup();
    seen
}

#[test]
fn a_change_at_the_edge_of_a_contig_names_a_base_the_contig_has() {
    let dir = TempDir::new().unwrap();
    let rows = run(dir.path(), true);
    let one = |id: &str, allele: &str| {
        let h = hgvsg_of(&rows, id, allele);
        assert_eq!(h.len(), 1, "{id} {allele}: the rows disagree: {h:?}");
        h[0].clone()
    };
    // VEP 116.3 writes nothing for these: `crates/fastvep-annotate/fixtures/vep_hgvsg` holds the
    // same cases at both ends of two contigs, with what VEP itself wrote for each
    assert_eq!(one("ins_before_1", "T"), "-");
    assert_eq!(one("ins_after_last", "C"), "-");
    assert_eq!(
        one("dup_of_last", "T"),
        "-",
        "a dup of the last base: VEP writes none"
    );
    assert_eq!(one("ref_off_the_end", "-"), "-");
    // ... and the same for an allele of a record of several, where VEP writes a position that
    // does not exist (`chr10:g.0dup`, `chr10:g.100_101insC`)
    assert_eq!(one("multi_ins_before_1", "TG"), "-");
    assert_eq!(one("multi_ins_after_last", "TC"), "-");
    // what is inside the contig is described, at its very first and very last base too
    assert_eq!(one("multi_ins_before_1", "A"), "chr10:g.1G>A");
    assert_eq!(one("multi_ins_after_last", "A"), "chr10:g.100T>A");
    assert_eq!(one("snv_1", "A"), "chr10:g.1G>A");
    assert_eq!(one("snv_last", "C"), "chr10:g.100T>C");
    assert_eq!(one("del_last", "-"), "chr10:g.100del");
    assert_eq!(one("ins_inside_near_last", "C"), "chr10:g.99_100insC");
}

#[test]
fn without_a_reference_the_first_base_of_a_contig_is_still_known() {
    let dir = TempDir::new().unwrap();
    let rows = run(dir.path(), false);
    let one = |id: &str, allele: &str| {
        let h = hgvsg_of(&rows, id, allele);
        assert_eq!(h.len(), 1, "{id} {allele}: the rows disagree: {h:?}");
        h[0].clone()
    };
    assert_eq!(one("ins_before_1", "T"), "-");
    assert_eq!(one("multi_ins_before_1", "TG"), "-");
    assert_eq!(one("snv_1", "A"), "chr10:g.1G>A");
    // the length of the contig is not known, so the last base is not: the description stands
    assert_eq!(one("snv_last", "C"), "chr10:g.100T>C");
    assert_eq!(one("ins_after_last", "C"), "chr10:g.100_101insC");
}
