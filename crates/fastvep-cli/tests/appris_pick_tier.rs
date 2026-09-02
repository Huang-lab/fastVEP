//! `--pick`'s APPRIS tier has to be able to decide something.
//!
//! `pick.rs` scores eight tiers in VEP's `--pick_order`, and three unit tests
//! covered the APPRIS one. All three passed, and the tier had never decided a
//! pick in production: they built `TranscriptVariation` by hand, and the GFF3
//! parser wrote `appris: None` on every transcript it had ever produced. The
//! criterion returned the same constant for every candidate, so it could only
//! ever be a no-op, and the APPRIS column of every CSQ line was empty.
//!
//! The tier matters exactly where two genes overlap: each gene has its own
//! canonical transcript, so `canonical` cannot separate them, and if neither is
//! MANE the next tier that can is APPRIS. That is the case this file builds -
//! two overlapping protein-coding genes, both canonical, neither MANE, equal
//! TSL, equal biotype, equal consequence - so every tier above and below APPRIS
//! is tied and the pick is decided by APPRIS or by nothing.
//!
//! Driven through `run_annotate` rather than `pick_best_transcript_idx` on
//! purpose: the defect was in the parser, and a test that builds the struct
//! itself cannot see it. That is what the three unit tests were.

use std::fs::File;
use std::io::Write;
use std::path::{Path, PathBuf};

use fastvep_cli::pipeline::{run_annotate, AnnotateConfig};
use tempfile::TempDir;

/// Two overlapping protein-coding genes with GENCODE's APPRIS tags.
///
/// `ENST_ALT` sorts before `ENST_PRIN` alphabetically, which is `pick`'s final
/// tie-break, so with the APPRIS tier inert the wrong one wins. Reversing the
/// two IDs would make this test pass against the broken code.
fn write_gff3(dir: &Path) -> PathBuf {
    let path = dir.join("ann.gff3");
    std::fs::write(
        &path,
        "1\ttest\tgene\t1001\t1200\t.\t+\t.\tID=gene:ENSG_A;Name=GENEA;biotype=protein_coding\n\
         1\ttest\tmRNA\t1001\t1200\t.\t+\t.\tID=transcript:ENST_PRIN;Parent=gene:ENSG_A;biotype=protein_coding;tag=gencode_basic,Ensembl_canonical,appris_principal_1;transcript_support_level=1\n\
         1\ttest\texon\t1001\t1200\t.\t+\t.\tID=exon:E_PRIN;Parent=transcript:ENST_PRIN;rank=1\n\
         1\ttest\tCDS\t1001\t1200\t.\t+\t0\tID=CDS:P_PRIN;Parent=transcript:ENST_PRIN\n\
         1\ttest\tgene\t1001\t1200\t.\t+\t.\tID=gene:ENSG_B;Name=GENEB;biotype=protein_coding\n\
         1\ttest\tmRNA\t1001\t1200\t.\t+\t.\tID=transcript:ENST_ALT;Parent=gene:ENSG_B;biotype=protein_coding;tag=gencode_basic,Ensembl_canonical,appris_alternative_2;transcript_support_level=1\n\
         1\ttest\texon\t1001\t1200\t.\t+\t.\tID=exon:E_ALT;Parent=transcript:ENST_ALT;rank=1\n\
         1\ttest\tCDS\t1001\t1200\t.\t+\t0\tID=CDS:P_ALT;Parent=transcript:ENST_ALT\n",
    )
    .unwrap();
    path
}

fn write_fasta(dir: &Path) -> PathBuf {
    let mut seq = vec![b'N'; 1_300];
    // 200 coding bases starting at 1-based 1001, opening with ATG.
    seq[1000..1003].copy_from_slice(b"ATG");
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
    writeln!(f, "1\t1100\t.\tN\tA\t.\t.\t.").unwrap();
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
        hgvs: false,
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

fn csq_entries(annotated: &str) -> Vec<Vec<String>> {
    let line = annotated
        .lines()
        .find(|l| !l.starts_with('#'))
        .unwrap_or_else(|| panic!("no annotated record in:\n{annotated}"));
    let info = line.split('\t').nth(7).expect("INFO column");
    let csq = info
        .split(';')
        .find_map(|kv| kv.strip_prefix("CSQ="))
        .unwrap_or_else(|| panic!("no CSQ on:\n{line}"));
    csq.split(',')
        .map(|e| e.split('|').map(str::to_string).collect())
        .collect()
}

fn field_idx(name: &str) -> usize {
    fastvep_io::output::DEFAULT_CSQ_FIELDS
        .iter()
        .position(|f| *f == name)
        .unwrap_or_else(|| panic!("{name} is a default CSQ field"))
}

fn annotate(dir: &Path, pick: bool, pick_order: Option<&str>) -> String {
    let gff3 = write_gff3(dir);
    let fasta = write_fasta(dir);
    let vcf = write_vcf(dir);
    let out = dir.join(format!(
        "out-{}-{}.vcf",
        pick,
        pick_order.unwrap_or("default")
    ));
    run_annotate(AnnotateConfig {
        gff3: vec![gff3.to_string_lossy().into()],
        fasta: Some(fasta.to_string_lossy().into()),
        pick,
        pick_order: pick_order.map(str::to_string),
        ..config(&vcf, &out)
    })
    .expect("annotation should succeed");
    std::fs::read_to_string(&out).unwrap()
}

#[test]
fn a_gencode_appris_tag_reaches_the_csq_column() {
    let dir = TempDir::new().unwrap();
    let annotated = annotate(dir.path(), false, None);

    let (feature, appris) = (field_idx("Feature"), field_idx("APPRIS"));
    let entries = csq_entries(&annotated);

    let of = |tx: &str| {
        entries
            .iter()
            .find(|f| f.get(feature).map(String::as_str) == Some(tx))
            .unwrap_or_else(|| panic!("no CSQ entry for {tx} in:\n{annotated}"))
            .get(appris)
            .cloned()
            .unwrap_or_default()
    };

    // Empty for both before the parser read `tag=appris_*`.
    assert_eq!(of("ENST_PRIN"), "P1", "in:\n{annotated}");
    assert_eq!(of("ENST_ALT"), "ALT2", "in:\n{annotated}");
}

#[test]
fn pick_prefers_the_appris_principal_transcript() {
    let dir = TempDir::new().unwrap();
    let annotated = annotate(dir.path(), true, None);

    let feature = field_idx("Feature");
    let entries = csq_entries(&annotated);
    assert_eq!(entries.len(), 1, "--pick kept more than one: {entries:?}");
    assert_eq!(
        entries[0].get(feature).map(String::as_str),
        Some("ENST_PRIN"),
        "the APPRIS tier did not decide the pick; `ENST_ALT` is what the \
         alphabetical transcript-ID tie-break returns, in:\n{annotated}"
    );
}

#[test]
fn an_explicit_appris_first_pick_order_is_honoured() {
    let dir = TempDir::new().unwrap();
    // APPRIS ahead of everything. Before the parser populated the field this
    // was accepted and silently changed nothing, which is the behaviour
    // `parse_pick_order` refuses to allow for `length`.
    let annotated = annotate(dir.path(), true, Some("appris,mane_select,canonical,rank"));

    let feature = field_idx("Feature");
    let entries = csq_entries(&annotated);
    assert_eq!(entries.len(), 1, "--pick kept more than one: {entries:?}");
    assert_eq!(
        entries[0].get(feature).map(String::as_str),
        Some("ENST_PRIN"),
        "`--pick-order appris,...` did not change the pick, in:\n{annotated}"
    );
}
