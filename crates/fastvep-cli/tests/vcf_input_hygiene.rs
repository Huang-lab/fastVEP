//! What `fastvep annotate` does with a VCF that is not what a pipeline hoped for, as a user meets
//! it: the program run over the small world of `fixtures/vep_tab_ablation` (one contig, `chr1`,
//! and eight transcripts), with its stderr and exit status read.
//!
//! The behaviour is Ensembl VEP's where VEP has one (issue #139): the alleles are annotated in
//! upper case, an ALT that is the reference or a repeat is no allele, a record that cannot be read
//! or is on a contig the gene models do not have is skipped with a warning and the run goes on.
//! The one thing added is that a run in which no record is on a known contig is an error, which
//! VEP only reports as a list of warnings and an empty file.

use std::path::{Path, PathBuf};
use std::process::{Command, Output};

use tempfile::TempDir;

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("fixtures/vep_tab_ablation")
        .join(name)
}

const HEADER: &str = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

/// A copy of the world in a temporary directory (the transcript cache a run writes goes beside the
/// GFF3, and neither it nor the FASTA's index belongs in the source tree) and the VCF to run.
struct Run {
    dir: TempDir,
}

impl Run {
    fn new() -> Self {
        let dir = TempDir::new().unwrap();
        for f in ["genes.gff3", "genome.fa"] {
            std::fs::copy(fixture(f), dir.path().join(f)).unwrap();
        }
        Self { dir }
    }

    /// `fastvep annotate` over `vcf` with `extra` flags: the process output and what it wrote.
    fn annotate(&self, vcf: &str, format: &str, extra: &[&str]) -> (Output, String) {
        let input = self.dir.path().join("in.vcf");
        std::fs::write(&input, vcf).unwrap();
        let out = self.dir.path().join(format!("out.{format}"));
        let _ = std::fs::remove_file(&out);
        let output = Command::new(env!("CARGO_BIN_EXE_fastvep"))
            .arg("annotate")
            .arg("-i")
            .arg(&input)
            .arg("-o")
            .arg(&out)
            .arg("--gff3")
            .arg(self.dir.path().join("genes.gff3"))
            .arg("--fasta")
            .arg(self.dir.path().join("genome.fa"))
            .args(["--output-format", format, "--no-progress", "--hgvs"])
            .args(extra)
            .output()
            .expect("the fastvep binary runs");
        let text = std::fs::read_to_string(&out).unwrap_or_default();
        (output, text)
    }
}

fn stderr(output: &Output) -> String {
    String::from_utf8_lossy(&output.stderr).into()
}

/// The data rows of a tab file, without the column that names the record, so that the rows of
/// two spellings of one record can be compared.
fn rows(tab: &str, drop_name: bool) -> Vec<String> {
    tab.lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| {
            if drop_name {
                l.split_once('\t').map_or("", |(_, rest)| rest).to_string()
            } else {
                l.to_string()
            }
        })
        .collect()
}

fn vcf(records: &[&str]) -> String {
    let mut text = HEADER.to_string();
    for r in records {
        text.push_str(r);
        text.push('\n');
    }
    text
}

/// A SNV, a deletion and an insertion inside the coding region of transcript A1 (`chr1:400-459`,
/// forward strand), which the rest of the tests spell in other ways.
fn reference_records(run: &Run) -> Vec<String> {
    let genome: String = std::fs::read_to_string(fixture("genome.fa"))
        .unwrap()
        .lines()
        .filter(|l| !l.starts_with('>'))
        .collect();
    let base = |p: usize| genome.as_bytes()[p - 1] as char;
    let _ = run;
    let other = |c: char| match c {
        'A' => 'C',
        'C' => 'G',
        'G' => 'T',
        _ => 'A',
    };
    vec![
        format!("chr1\t430\t.\t{}\t{}\t.\t.\t.", base(430), other(base(430))),
        format!(
            "chr1\t440\t.\t{}{}{}\t{}\t.\t.\t.",
            base(440),
            base(441),
            base(442),
            base(440)
        ),
        format!("chr1\t445\t.\t{}\t{}CG\t.\t.\t.", base(445), base(445)),
    ]
}

/// Alleles in lower case give VEP's answer for their upper-case twins: the same rows, and the
/// same `Allele`, `HGVSc` and `VARIANT_CLASS` (which was `sequence_alteration` for a lone SNV).
#[test]
fn a_lower_case_vcf_is_annotated_as_its_upper_case_twin() {
    let run = Run::new();
    let upper = reference_records(&run);
    let lower: Vec<String> = upper
        .iter()
        .map(|l| {
            let mut f: Vec<String> = l.split('\t').map(String::from).collect();
            f[3] = f[3].to_ascii_lowercase();
            f[4] = f[4].to_ascii_lowercase();
            f.join("\t")
        })
        .collect();
    let u: Vec<&str> = upper.iter().map(String::as_str).collect();
    let l: Vec<&str> = lower.iter().map(String::as_str).collect();
    let flags = ["--variant-class", "--hgvsg"];
    let (o1, tab_upper) = run.annotate(&vcf(&u), "tab", &flags);
    let (o2, tab_lower) = run.annotate(&vcf(&l), "tab", &flags);
    assert!(
        o1.status.success() && o2.status.success(),
        "{}",
        stderr(&o2)
    );
    assert!(rows(&tab_upper, true).len() >= 3);
    assert_eq!(rows(&tab_lower, true), rows(&tab_upper, true));
    // the record keeps the spelling it was written in, as VEP's `Uploaded_variation` does
    assert!(tab_lower.lines().any(|l| l.starts_with("chr1_430_")
        && l.contains("_430_")
        && l.split('\t')
            .next()
            .unwrap()
            .chars()
            .any(|c| c.is_ascii_lowercase())));
    assert!(!tab_lower.contains("sequence_alteration"));
}

/// An ALT that is the reference is no allele and one that repeats an earlier ALT is no other:
/// the rows are the single-ALT record's, and the VCF output has one line for each record.
#[test]
fn an_alt_that_repeats_the_reference_or_another_alt_writes_no_row_of_its_own() {
    let run = Run::new();
    let records = reference_records(&run);
    let snv = &records[0];
    let fields: Vec<&str> = snv.split('\t').collect();
    let (pos, reference, alt) = (fields[1], fields[3], fields[4]);
    let single = vcf(&[snv.as_str()]);
    let repeated = vcf(&[format!("chr1\t{pos}\t.\t{reference}\t{alt},{alt}\t.\t.\t.").as_str()]);
    let (o, tab_single) = run.annotate(&single, "tab", &[]);
    assert!(o.status.success(), "{}", stderr(&o));
    let (o, tab_repeated) = run.annotate(&repeated, "tab", &[]);
    assert!(o.status.success(), "{}", stderr(&o));
    assert!(!rows(&tab_single, false).is_empty());
    assert_eq!(rows(&tab_repeated, true), rows(&tab_single, true));

    // ALT equal to REF: no row at all, and the line is still in the VCF
    let same = vcf(&[format!("chr1\t{pos}\t.\t{reference}\t{reference}\t.\t.\t.").as_str()]);
    let (o, tab) = run.annotate(&same, "tab", &[]);
    assert!(o.status.success(), "{}", stderr(&o));
    assert!(rows(&tab, false).is_empty(), "{tab}");
    let (_, out_vcf) = run.annotate(&same, "vcf", &[]);
    assert_eq!(out_vcf.lines().filter(|l| !l.starts_with('#')).count(), 1);
    assert!(!out_vcf.contains("CSQ="), "{out_vcf}");
}

/// A record that cannot be read is skipped with its line number, the others are annotated as
/// they are without it, and the run says how many it skipped.
#[test]
fn a_record_that_cannot_be_read_is_skipped_and_counted() {
    let run = Run::new();
    let records = reference_records(&run);
    let clean = vcf(&records.iter().map(String::as_str).collect::<Vec<_>>());
    let (o, tab_clean) = run.annotate(&clean, "tab", &[]);
    assert!(o.status.success(), "{}", stderr(&o));

    let mut dirty = records.clone();
    dirty.insert(1, "chr1\t0\t.\tA\tC\t.\t.\t.".into());
    dirty.insert(2, "chr1\tx\t.\tA\tC\t.\t.\t.".into());
    dirty.push("too few fields".into());
    dirty.push("chr1\t450\t.\tA\t\t.\t.\t.".into());
    let (o, tab_dirty) = run.annotate(
        &vcf(&dirty.iter().map(String::as_str).collect::<Vec<_>>()),
        "tab",
        &[],
    );
    assert!(o.status.success(), "{}", stderr(&o));
    assert_eq!(rows(&tab_dirty, false), rows(&tab_clean, false));
    let err = stderr(&o);
    assert!(err.contains("skipped line 4 of the VCF"), "{err}");
    assert!(err.contains("skipped line 5 of the VCF"), "{err}");
    assert!(
        err.contains("4 record(s) of the VCF could not be read"),
        "{err}"
    );
}

/// A record on a contig that is not in the gene models or the reference is said to be, once for
/// the contig, with the contig of the run that is the same in another case when there is one, and is
/// still annotated as intergenic (a `--sa-dir` may have something for it); and a run in which none
/// is on a known contig is an error and not a file of intergenic rows.
#[test]
fn a_contig_the_gene_models_do_not_have_is_said_and_not_left_for_intergenic() {
    let run = Run::new();
    let records = reference_records(&run);
    let mut mixed = records.clone();
    mixed.push("chrZ\t100\t.\tA\tC\t.\t.\t.".into());
    mixed.push("chrZ\t200\t.\tA\tC\t.\t.\t.".into());
    mixed.push("CHR1\t430\t.\tA\tC\t.\t.\t.".into());
    let (o, tab) = run.annotate(
        &vcf(&mixed.iter().map(String::as_str).collect::<Vec<_>>()),
        "tab",
        &[],
    );
    assert!(o.status.success(), "{}", stderr(&o));
    let err = stderr(&o);
    assert!(
        err.matches("records on `chrZ`").count() == 1
            && err.matches("records on `CHR1`").count() == 1,
        "once for the contig: {err}"
    );
    assert!(
        err.contains("have `chr1`; contig names are case-sensitive"),
        "{err}"
    );
    assert!(
        err.contains(
            "3 record(s) are on contigs that the gene models and the reference do not have"
        ) && err.contains("CHR1 (1), chrZ (2)"),
        "{err}"
    );
    // the records of the contigs it does not know are intergenic, as they were
    let intergenic = rows(&tab, false)
        .into_iter()
        .filter(|r| r.contains("intergenic_variant"))
        .count();
    assert_eq!(intergenic, 3, "{tab}");

    // all of them on contigs it does not know: an error, whatever the output format
    let (o, _) = run.annotate(
        &vcf(&["Chr1\t430\t.\tA\tC\t.\t.\t.", "Chr1\t440\t.\tA\tC\t.\t.\t."]),
        "vcf",
        &[],
    );
    assert!(!o.status.success());
    let err = stderr(&o);
    assert!(
        err.contains("none of the 2 records of the VCF is on a contig"),
        "{err}"
    );
    assert!(err.contains("Chr1 (2)") && err.contains("chr1"), "{err}");
}

/// More contigs than are worth a line each are in the summary, not in a thousand warnings.
#[test]
fn a_vcf_with_many_unknown_contigs_does_not_write_a_warning_for_each() {
    let run = Run::new();
    let records = reference_records(&run);
    let mut many = records.clone();
    for i in 0..30 {
        many.push(format!("alt_{i}\t100\t.\tA\tC\t.\t.\t."));
    }
    let (o, _) = run.annotate(
        &vcf(&many.iter().map(String::as_str).collect::<Vec<_>>()),
        "tab",
        &[],
    );
    assert!(o.status.success(), "{}", stderr(&o));
    let err = stderr(&o);
    assert_eq!(
        err.matches("which is not a contig of the gene models")
            .count(),
        10,
        "{err}"
    );
    assert!(
        err.contains("30 record(s) are on contigs") && err.contains("and 22 more"),
        "{err}"
    );
}

/// A byte order mark in front of the first line, and a run of blank lines and comments among the
/// records, do not change what is annotated.
#[test]
fn a_byte_order_mark_and_stray_lines_change_nothing() {
    let run = Run::new();
    let records = reference_records(&run);
    let clean = vcf(&records.iter().map(String::as_str).collect::<Vec<_>>());
    let (o, tab_clean) = run.annotate(&clean, "tab", &[]);
    assert!(o.status.success(), "{}", stderr(&o));

    let mut odd = String::from("\u{feff}");
    odd.push_str(HEADER);
    odd.push_str(&records[0]);
    odd.push('\n');
    odd.push_str(&"\n#a comment\n".repeat(50_000));
    odd.push_str(&records[1]);
    odd.push_str("\r\n");
    odd.push_str(&records[2]);
    odd.push('\n');
    let (o, tab_odd) = run.annotate(&odd, "tab", &[]);
    assert!(o.status.success(), "{}", stderr(&o));
    assert_eq!(rows(&tab_odd, false), rows(&tab_clean, false));
}

/// A REF that is not the reference is described against the reference, as Ensembl VEP does: the
/// HGVSc of `A>C` at a base that is a `T` is that of `T>C`, which is what the `Codons` and the
/// HGVSg of the same row said already. The run says how many records disagree, with the first.
#[test]
fn a_ref_that_is_not_the_reference_is_described_against_the_reference_and_said() {
    let run = Run::new();
    let genome: String = std::fs::read_to_string(fixture("genome.fa"))
        .unwrap()
        .lines()
        .filter(|l| !l.starts_with('>'))
        .collect();
    let base = |p: usize| genome.as_bytes()[p - 1] as char;
    let wrong = |c: char| if c == 'A' { 'C' } else { 'A' };
    let alt = |c: char| match c {
        'A' => 'G',
        'C' => 'T',
        'G' => 'C',
        _ => 'A',
    };
    let right_snv = format!("chr1\t430\t.\t{}\t{}\t.\t.\t.", base(430), alt(base(430)));
    let wrong_snv = format!(
        "chr1\t430\t.\t{}\t{}\t.\t.\t.",
        wrong(base(430)),
        alt(base(430))
    );
    // a deletion of three bases whose last two REF bases are wrong
    let (b0, b1, b2) = (base(440), base(441), base(442));
    let right_del = format!(
        "chr1\t439\t.\t{}{b0}{b1}{b2}\t{}\t.\t.\t.",
        base(439),
        base(439)
    );
    let wrong_del = format!(
        "chr1\t439\t.\t{}{b0}{}{}\t{}\t.\t.\t.",
        base(439),
        wrong(b1),
        wrong(b2),
        base(439)
    );
    let flags = ["--hgvsg"];
    let (o, right) = run.annotate(
        &vcf(&[right_snv.as_str(), right_del.as_str()]),
        "tab",
        &flags,
    );
    assert!(o.status.success(), "{}", stderr(&o));
    assert!(!stderr(&o).contains("REF allele"), "{}", stderr(&o));
    let (o, tab) = run.annotate(
        &vcf(&[wrong_snv.as_str(), wrong_del.as_str()]),
        "tab",
        &flags,
    );
    assert!(o.status.success(), "{}", stderr(&o));

    let column = |text: &str, name: &str| -> Vec<String> {
        let header: Vec<&str> = text
            .lines()
            .find(|l| l.starts_with("#Uploaded_variation"))
            .unwrap()
            .trim_start_matches('#')
            .split('\t')
            .collect();
        let at = header.iter().position(|c| *c == name).unwrap();
        text.lines()
            .filter(|l| !l.starts_with('#'))
            .map(|l| l.split('\t').nth(at).unwrap().to_string())
            .collect()
    };
    for name in [
        "HGVSc",
        "HGVSp",
        "HGVSg",
        "Codons",
        "Amino_acids",
        "Consequence",
    ] {
        assert_eq!(column(&tab, name), column(&right, name), "{name}");
    }
    assert!(column(&tab, "HGVSc").iter().any(|c| c.contains("c.")));
    let err = stderr(&o);
    assert!(
        err.contains("2 record(s) have a REF allele that is not the reference's")
            && err.contains("chr1:430, where the VCF has REF")
            && err.contains("same assembly"),
        "{err}"
    );
}
