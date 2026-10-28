//! The CDS flags an Ensembl GFF3 export leaves out, read from the sequence, are not read on the
//! contig of an organelle genome.
//!
//! Ensembl's GFF3 names a transcript `transcript:ENST...` and has no `cds_start_NF` or
//! `cds_end_NF`, which Ensembl's database holds, so a run with a FASTA reads them off the codons
//! with the standard genetic code (`Transcript::read_ensembl_cds_state`). The rules were checked
//! against the nuclear genes of human chr21 and chr22. A gene of a yeast `Mito` contig that starts
//! with `ATA`, which that genome's code reads as a start and the standard one as Ile, came out
//! `cds_start_NF`, and was then annotated as an incomplete transcript is: a change in its first codon
//! was a `missense_variant` and not a `start_lost`.
//!
//! The same gene on two contigs, one nuclear and one named as the yeast, Drosophila and
//! Arabidopsis genomes of the test data name theirs, run as a program because the answer is in the
//! `FLAGS` field of the `CSQ`.

use std::io::Write;
use std::path::Path;
use std::process::Command;

use tempfile::TempDir;

/// The coding sequence of the gene: `ATA` with no UTR above it, then `AAA TGG TAA`.
const CDS: &str = "ATAAAATGGTAA";

/// One gene on `contig`, its CDS the whole of its one exon, from base 100.
fn gff3_line(contig: &str, n: usize) -> String {
    let (start, end) = (100, 100 + CDS.len() as u64 - 1);
    format!(
        "{contig}\tensembl\tgene\t{start}\t{end}\t.\t+\t.\tID=gene:ENSG_{n};Name=G{n};biotype=protein_coding;gene_id=ENSG_{n}\n\
         {contig}\tensembl\tmRNA\t{start}\t{end}\t.\t+\t.\tID=transcript:ENST_{n};Parent=gene:ENSG_{n};biotype=protein_coding;transcript_id=ENST_{n}\n\
         {contig}\tensembl\texon\t{start}\t{end}\t.\t+\t.\tParent=transcript:ENST_{n};rank=1\n\
         {contig}\tensembl\tCDS\t{start}\t{end}\t.\t+\t0\tID=CDS:ENSP_{n};Parent=transcript:ENST_{n};protein_id=ENSP_{n}\n"
    )
}

fn write_fasta(path: &Path, contigs: &[&str]) {
    let mut f = std::fs::File::create(path).unwrap();
    for contig in contigs {
        let mut seq = vec![b'C'; 140];
        seq[99..99 + CDS.len()].copy_from_slice(CDS.as_bytes());
        writeln!(f, ">{contig}").unwrap();
        f.write_all(&seq).unwrap();
        writeln!(f).unwrap();
    }
}

/// `(contig, FLAGS, Consequence)` of the transcript's row for a change of the first base of its CDS.
fn run(contigs: &[&str]) -> Vec<(String, String, String)> {
    let dir = TempDir::new().unwrap();
    let gff3 = dir.path().join("genes.gff3");
    std::fs::write(
        &gff3,
        contigs
            .iter()
            .enumerate()
            .map(|(n, c)| gff3_line(c, n))
            .collect::<String>(),
    )
    .unwrap();
    let fasta = dir.path().join("genome.fa");
    write_fasta(&fasta, contigs);
    let vcf = dir.path().join("in.vcf");
    let mut f = std::fs::File::create(&vcf).unwrap();
    writeln!(f, "##fileformat=VCFv4.2").unwrap();
    writeln!(f, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO").unwrap();
    for contig in contigs {
        writeln!(f, "{contig}\t100\t.\tA\tG\t.\t.\t.").unwrap();
    }
    let out = dir.path().join("out.vcf");
    let output = Command::new(env!("CARGO_BIN_EXE_fastvep"))
        .arg("annotate")
        .arg("-i")
        .arg(&vcf)
        .arg("-o")
        .arg(&out)
        .arg("--gff3")
        .arg(&gff3)
        .arg("--fasta")
        .arg(&fasta)
        .args(["--output-format", "vcf", "--no-progress"])
        .output()
        .expect("the fastvep binary runs");
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let vcf = std::fs::read_to_string(&out).unwrap();
    let format = vcf
        .lines()
        .find(|l| l.starts_with("##INFO=<ID=CSQ"))
        .and_then(|l| l.split("Format: ").nth(1))
        .expect("a CSQ header")
        .trim_end_matches("\">")
        .split('|')
        .map(String::from)
        .collect::<Vec<_>>();
    let col = |name: &str| format.iter().position(|c| c == name).expect(name);
    let (flags, consequence, feature) = (col("FLAGS"), col("Consequence"), col("Feature"));
    let mut rows = Vec::new();
    for line in vcf.lines().filter(|l| !l.starts_with('#')) {
        let fields: Vec<&str> = line.split('\t').collect();
        for entry in fields[7].trim_start_matches("CSQ=").split(',') {
            let cells: Vec<&str> = entry.split('|').collect();
            if cells[feature].starts_with("ENST_") {
                rows.push((
                    fields[0].to_string(),
                    cells[flags].to_string(),
                    cells[consequence].to_string(),
                ));
            }
        }
    }
    rows
}

#[test]
fn a_gene_of_an_organelle_contig_is_not_flagged_incomplete_by_the_standard_code() {
    for organelle in ["Mito", "mitochondrion_genome", "Mt", "Pt"] {
        let rows = run(&["1", organelle]);
        let on = |contig: &str| {
            rows.iter()
                .find(|(c, _, _)| c == contig)
                .unwrap_or_else(|| panic!("no row for {contig}: {rows:?}"))
        };
        // nuclear: `ATA` with no UTR above it is what the cache calls an incomplete start
        let nuclear = on("1");
        assert_eq!(nuclear.1, "cds_start_NF", "{rows:?}");
        assert!(
            nuclear.2.contains("missense_variant") && !nuclear.2.contains("start_lost"),
            "{rows:?}"
        );
        // the organelle genome is not read with the standard code, and its gene starts where it says
        let mine = on(organelle);
        assert_eq!(mine.1, "", "{organelle}: {rows:?}");
        assert!(mine.2.contains("start_lost"), "{organelle}: {rows:?}");
    }
}
