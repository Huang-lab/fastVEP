//! A gene Ensembl types `ncRNA_gene` must reach the CSQ column with its name.
//!
//! Reported as #98: `fastvep annotate --symbol` printed an empty SYMBOL where
//! VEP printed `WASH7P`, for every variant in a non-coding locus, while the
//! Gene column held the right stable ID. The ID survived because it comes from
//! the transcript's `Parent`; the name did not, because the GFF3 parser matched
//! only `gene` and `pseudogene` and dropped the whole `ncRNA_gene` record.
//!
//! The unit tests in `fastvep-cache::gff` pin the parser. This one pins what
//! the user actually sees: run the reporter's command shape end to end, and
//! read the SYMBOL out of the VCF that comes back.

use std::fs::File;
use std::io::Write;
use std::path::{Path, PathBuf};

use fastvep_cli::pipeline::{run_annotate, AnnotateConfig};
use tempfile::TempDir;

/// One non-coding locus and one protein-coding one, typed the way Ensembl types
/// them: `ncRNA_gene`/`lnc_RNA` for the first, `gene`/`mRNA` for the second. The
/// coding gene is the control — it was never affected, so if it also loses its
/// symbol the failure is somewhere other than the gene-type match.
fn write_gff3(dir: &Path) -> PathBuf {
    let path = dir.join("ann.gff3");
    std::fs::write(
        &path,
        "1\ttest\tncRNA_gene\t1001\t1100\t.\t+\t.\tID=gene:ENSG_NC;Name=WASH7P;biotype=lncRNA\n\
         1\ttest\tlnc_RNA\t1001\t1100\t.\t+\t.\tID=transcript:ENST_NC;Parent=gene:ENSG_NC;biotype=lncRNA;tag=Ensembl_canonical\n\
         1\ttest\texon\t1001\t1100\t.\t+\t.\tID=exon:ENSE_NC;Parent=transcript:ENST_NC;rank=1\n\
         1\ttest\tgene\t2001\t2012\t.\t+\t.\tID=gene:ENSG_PC;Name=CTRLGENE;biotype=protein_coding\n\
         1\ttest\tmRNA\t2001\t2012\t.\t+\t.\tID=transcript:ENST_PC;Parent=gene:ENSG_PC;biotype=protein_coding;tag=Ensembl_canonical\n\
         1\ttest\texon\t2001\t2012\t.\t+\t.\tID=exon:ENSE_PC;Parent=transcript:ENST_PC;rank=1\n\
         1\ttest\tCDS\t2001\t2012\t.\t+\t0\tID=CDS:ENSP_PC;Parent=transcript:ENST_PC\n",
    )
    .unwrap();
    path
}

fn write_fasta(dir: &Path) -> PathBuf {
    let mut seq = vec![b'N'; 2_100];
    seq[2000..2012].copy_from_slice(b"ATGTGGCGGTAA");
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
    // In the ncRNA exon, then in the coding exon.
    writeln!(f, "1\t1050\t.\tN\tA\t.\t.\t.").unwrap();
    writeln!(f, "1\t2008\t.\tG\tA\t.\t.\t.").unwrap();
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

/// The SYMBOL that the record at `pos` carries for transcript `feature`.
///
/// Keyed on the transcript because both loci fall inside the default 5 kb
/// window, so every record carries a CSQ entry for each of them.
fn symbol_for(annotated: &str, pos: &str, feature: &str) -> String {
    let fields = fastvep_io::output::DEFAULT_CSQ_FIELDS;
    let idx = |name: &str| {
        fields
            .iter()
            .position(|f| *f == name)
            .unwrap_or_else(|| panic!("{name} is a default CSQ field"))
    };
    let (symbol_idx, feature_idx) = (idx("SYMBOL"), idx("Feature"));

    let line = annotated
        .lines()
        .find(|l| !l.starts_with('#') && l.split('\t').nth(1) == Some(pos))
        .unwrap_or_else(|| panic!("no record at position {pos} in:\n{annotated}"));
    let info = line.split('\t').nth(7).expect("INFO column");
    let csq = info
        .split(';')
        .find_map(|kv| kv.strip_prefix("CSQ="))
        .unwrap_or_else(|| panic!("no CSQ on the record at {pos}: {line}"));

    csq.split(',')
        .map(|entry| entry.split('|').collect::<Vec<_>>())
        .find(|f| f.get(feature_idx) == Some(&feature))
        .map(|f| f.get(symbol_idx).copied().unwrap_or_default().to_string())
        .unwrap_or_else(|| panic!("no CSQ entry for {feature} at {pos}: {csq}"))
}

#[test]
fn an_ncrna_gene_reaches_the_csq_with_its_symbol() {
    let dir = TempDir::new().unwrap();
    let gff3 = write_gff3(dir.path());
    let fasta = write_fasta(dir.path());
    let vcf = write_vcf(dir.path());
    let out = dir.path().join("out.vcf");

    run_annotate(AnnotateConfig {
        gff3: vec![gff3.to_string_lossy().into()],
        fasta: Some(fasta.to_string_lossy().into()),
        ..config(&vcf, &out)
    })
    .expect("annotation should succeed");

    let annotated = std::fs::read_to_string(&out).unwrap();

    assert_eq!(
        symbol_for(&annotated, "1050", "ENST_NC"),
        "WASH7P",
        "an ncRNA_gene lost its symbol (#98), in:\n{annotated}"
    );

    // The control: a plain `gene` was never affected, so this failing too would
    // mean the symbol is being lost somewhere other than the gene-type match.
    assert_eq!(
        symbol_for(&annotated, "2008", "ENST_PC"),
        "CTRLGENE",
        "a protein_coding gene lost its symbol, in:\n{annotated}"
    );
}
