//! `VARIANT_CLASS` of a record whose ALT list repeats an ALT or the reference.
//!
//! Ensembl's class comes from its allele string, which holds such an ALT as it was written
//! (`T>C,C` is `T/C/C`, an `SNV`; `AGT>AGT,GA,ACT` is `AGT/AGT/GA/ACT`, an `indel`). The parser
//! reads the repeat as a placeholder, which gives it no row of its own (Ensembl writes none for
//! one that repeats the reference and one for the second of two that are alike) and then stood in
//! the allele string too, where no class has a `.`: 10 of the 14 records below were a
//! `sequence_alteration`. What each should be is what Ensembl VEP 116.3 wrote for the record,
//! run on all of them.
//!
//! Run as a program, because the class is a column of the output.

use std::io::Write;
use std::path::Path;
use std::process::Command;

use tempfile::TempDir;

/// `(id, REF, ALT, VARIANT_CLASS of Ensembl VEP 116.3)`
const RECORDS: [(&str, &str, &str, &str); 14] = [
    ("r01", "T", "C,C", "SNV"),
    ("r02", "T", "T,C", "SNV"),
    ("r03", "T", "C,T,C", "SNV"),
    ("r04", "TA", "T,T", "deletion"),
    ("r05", "T", "TA,TA", "insertion"),
    ("r06", "AGT", "AGT,GA,ACT", "indel"),
    ("r07", "GGC", "CGAGG,GGCGC,CGAGG", "indel"),
    ("r08", "TA", "TA,T", "sequence_alteration"),
    ("r09", "T", "C,C,G", "SNV"),
    ("r10", "AC", "AC,GT", "substitution"),
    ("r11", "A", "A,AT,AT", "sequence_alteration"),
    ("r12", "AC", "GT,GT", "substitution"),
    ("r13", "TA", "T,CA,T", "indel"),
    ("r15", "T", "C", "SNV"),
];

fn run(dir: &Path) -> Vec<(String, String, String)> {
    let gff3 = dir.join("genes.gff3");
    std::fs::write(
        &gff3,
        "1\ttest\tgene\t10000\t10100\t.\t+\t.\tID=gene:G1;Name=G1;biotype=lncRNA\n\
         1\ttest\tlnc_RNA\t10000\t10100\t.\t+\t.\tID=transcript:T1;Parent=gene:G1;biotype=lncRNA\n\
         1\ttest\texon\t10000\t10100\t.\t+\t.\tID=exon:E1;Parent=transcript:T1;rank=1\n",
    )
    .unwrap();
    let vcf = dir.join("in.vcf");
    let mut f = std::fs::File::create(&vcf).unwrap();
    writeln!(f, "##fileformat=VCFv4.2").unwrap();
    writeln!(f, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO").unwrap();
    for (n, (id, r, a, _)) in RECORDS.iter().enumerate() {
        writeln!(f, "1\t{}\t{id}\t{r}\t{a}\t.\t.\t.", 100 + 100 * n).unwrap();
    }
    drop(f);
    let out = dir.join("out.tab");
    let output = Command::new(env!("CARGO_BIN_EXE_fastvep"))
        .args(["annotate", "-i"])
        .arg(&vcf)
        .arg("-o")
        .arg(&out)
        .arg("--gff3")
        .arg(&gff3)
        .args(["--output-format", "tab", "--variant-class", "--no-progress"])
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
    let (id, allele, class) = (
        col("#Uploaded_variation"),
        col("Allele"),
        col("VARIANT_CLASS"),
    );
    tab.lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            (
                f[id].to_string(),
                f[allele].to_string(),
                f[class].to_string(),
            )
        })
        .collect()
}

#[test]
fn a_record_that_repeats_an_alt_or_the_reference_has_the_class_of_its_alleles() {
    let dir = TempDir::new().unwrap();
    let rows = run(dir.path());
    for (id, reference, alts, class) in RECORDS {
        let mine: Vec<&str> = rows
            .iter()
            .filter(|(i, _, _)| i == id)
            .map(|(_, _, c)| c.as_str())
            .collect();
        assert!(!mine.is_empty(), "{id} {reference}>{alts}: no row");
        assert!(
            mine.iter().all(|c| *c == class),
            "{id} {reference}>{alts}: VEP 116.3 writes {class}, here {mine:?}"
        );
    }
}

/// The repeat is still no allele of its own: one row for the first of two alike, and none for an
/// ALT that is the reference.
#[test]
fn the_rows_are_as_they_were() {
    let dir = TempDir::new().unwrap();
    let rows = run(dir.path());
    let alleles = |id: &str| {
        let mut a: Vec<&str> = rows
            .iter()
            .filter(|(i, _, _)| i == id)
            .map(|(_, a, _)| a.as_str())
            .collect();
        a.sort_unstable();
        a.dedup();
        a
    };
    assert_eq!(alleles("r01"), ["C"], "T>C,C");
    assert_eq!(alleles("r02"), ["C"], "T>T,C");
    assert_eq!(alleles("r06"), ["ACT", "GA"], "AGT>AGT,GA,ACT");
    let count = |id: &str| rows.iter().filter(|(i, _, _)| i == id).count();
    assert_eq!(count("r01"), count("r15"), "one row for T>C,C as for T>C");
}
