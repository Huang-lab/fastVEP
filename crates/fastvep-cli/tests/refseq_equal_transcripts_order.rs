//! Of the transcripts of a RefSeq gene that are equally long, the one marked `CANONICAL` is the
//! first in the file, and a file sorted for tabix with a `sort` that is not stable has lost that.
//!
//! Ensembl VEP's RefSeq cache marks the first of the equals in NCBI's order (48 of the 49 genes of
//! chromosomes 21 and 22 that have a tie for the longest). `sort -k1,1 -k4,4n -k5,5n`, which is the
//! command the tests of this repository and issue #87 sort a GFF3 with, puts records of one start
//! and end in the order of their text, and 25 of those 49 genes then have another transcript marked;
//! `sort -s -k1,1 -k4,4n` keeps NCBI's order and changes none. Nothing in the sorted file says which
//! was done, so a run says so when every group of transcripts of one start and end is in the order
//! of its text (505 of 505 groups after the plain sort, 122 of 505 in NCBI's order, measured).
//!
//! Run as a program, because the answer is on stderr and in the `CANONICAL` column.

use std::path::Path;
use std::process::Command;

use tempfile::TempDir;

const GENES: u64 = 25;

/// One RefSeq gene with two transcripts of one extent, written in `order`.
fn gene(n: u64, order: [u64; 2]) -> String {
    let (start, end) = (1000 * n, 1000 * n + 99);
    let mut text = format!(
        "1\tBestRefSeq\tgene\t{start}\t{end}\t.\t+\t.\tID=gene-G{n};Dbxref=GeneID:{n};Name=G{n};gbkey=Gene;gene=G{n};gene_biotype=protein_coding\n"
    );
    for k in order {
        let acc = format!("NM_{n:03}0{k}.1");
        let common =
            format!("Dbxref=GeneID:{n},GenBank:{acc};gbkey=mRNA;gene=G{n};transcript_id={acc}");
        text.push_str(&format!(
            "1\tBestRefSeq\tmRNA\t{start}\t{end}\t.\t+\t.\tID=rna-{acc};Parent=gene-G{n};{common};Name={acc}\n\
             1\tBestRefSeq\texon\t{start}\t{end}\t.\t+\t.\tID=exon-{acc}-1;Parent=rna-{acc};{common}\n"
        ));
    }
    text
}

/// `(stderr, [(transcript, CANONICAL)])` of a run over a file whose genes have their transcripts in
/// `order`, with a variant in each gene.
fn run(order: [u64; 2], dir: &Path) -> (String, Vec<(String, String)>) {
    let gff3 = dir.join("refseq.gff3");
    std::fs::write(
        &gff3,
        (1..=GENES).map(|n| gene(n, order)).collect::<String>(),
    )
    .unwrap();
    let vcf = dir.join("in.vcf");
    let mut text =
        String::from("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
    for n in 1..=GENES {
        text.push_str(&format!("1\t{}\t.\tA\tG\t.\t.\t.\n", 1000 * n + 50));
    }
    std::fs::write(&vcf, text).unwrap();
    let out = dir.join("out.tab");
    let output = Command::new(env!("CARGO_BIN_EXE_fastvep"))
        .args(["annotate", "-i"])
        .arg(&vcf)
        .arg("-o")
        .arg(&out)
        .arg("--gff3")
        .arg(&gff3)
        .args(["--output-format", "tab", "--canonical", "--no-progress"])
        .output()
        .expect("the fastvep binary runs");
    let stderr = String::from_utf8_lossy(&output.stderr).to_string();
    assert!(output.status.success(), "{stderr}");
    let tab = std::fs::read_to_string(&out).unwrap();
    let header: Vec<&str> = tab
        .lines()
        .find(|l| l.starts_with("#Uploaded_variation"))
        .expect("a header")
        .split('\t')
        .collect();
    let col = |name: &str| header.iter().position(|c| *c == name).expect(name);
    let (feature, canonical) = (col("Feature"), col("CANONICAL"));
    let rows = tab
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            (f[feature].to_string(), f[canonical].to_string())
        })
        .collect();
    (stderr, rows)
}

fn canonical(rows: &[(String, String)]) -> Vec<&str> {
    let mut c: Vec<&str> = rows
        .iter()
        .filter(|(_, c)| c == "YES")
        .map(|(t, _)| t.as_str())
        .collect();
    c.sort_unstable();
    // a transcript has a row for every variant within reach of it
    c.dedup();
    c
}

#[test]
fn a_file_whose_equal_transcripts_are_in_the_order_of_their_text_is_said_to_be() {
    let dir = TempDir::new().unwrap();
    // the order `sort` leaves them in without `-s`: NM_...01 before NM_...02
    let (stderr, rows) = run([1, 2], dir.path());
    assert!(
        stderr.contains(&format!(
            "are in the order of their text ({GENES} of {GENES} such groups)"
        )) && stderr.contains(&format!("For {GENES} genes")),
        "no warning in:\n{stderr}"
    );
    assert!(stderr.contains("sort -s -k1,1 -k4,4n"), "{stderr}");
    // the first in the file, which the warning is about
    let expected: Vec<String> = (1..=GENES).map(|n| format!("NM_{n:03}01.1")).collect();
    assert_eq!(canonical(&rows), expected);
}

#[test]
fn a_file_whose_equal_transcripts_are_not_in_that_order_is_not_warned_about() {
    let dir = TempDir::new().unwrap();
    // NM_...02 before NM_...01 in every gene: not the order `sort` would leave, and the first in
    // the file is marked all the same
    let (stderr, rows) = run([2, 1], dir.path());
    assert!(
        !stderr.contains("order of their text"),
        "a warning for a file that was not sorted by text:\n{stderr}"
    );
    let expected: Vec<String> = (1..=GENES).map(|n| format!("NM_{n:03}02.1")).collect();
    assert_eq!(canonical(&rows), expected);
}
