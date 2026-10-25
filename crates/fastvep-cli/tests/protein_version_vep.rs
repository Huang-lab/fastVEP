//! The protein version in `HGVSp` is the GFF3's, as Ensembl VEP 116 reads it.
//!
//! Ensembl's GFF3 writes a protein's version on its CDS records,
//! `protein_id=ENSP00000363557;version=3`, and one protein in six of Ensembl 115's has one
//! that is not 1. VEP 116 puts it in the accession of every `HGVSp`
//! (`ENSP00000363557.3:p.Thr266Met`), with 1 where a record has none. VEP 115 gave every
//! translation version 1, and so did fastVEP, which parsed the version and then built every
//! transcript with `Some(1)`: on 20,241 records of GIAB HG002, 184 of the 844 `HGVSp` rows
//! differed from VEP 116.3 in nothing but that suffix.
//!
//! NCBI's GFF3 has no `version` attribute and carries the version in the id
//! (`protein_id=NP_000001.2`), which VEP writes as it is. fastVEP gave that protein version
//! 1 as well and appended it, `NP_000001.2.1`, an accession that does not exist.
//!
//! `fixtures/vep_protein_version` is four Ensembl-style coding genes whose proteins carry
//! versions 3, 1 (written out), none and 12, and one NCBI-style gene, with 29 variants in
//! them and the `HGVSp` Ensembl VEP 116.3 wrote for each (`regen.sh`, from the image pinned
//! by digest). This test runs fastVEP over the same files and holds every `HGVSp` to VEP's,
//! the whole string and not only the suffix, record by record: the NCBI-style transcript
//! is `rna-NM_000015.3` here and `NM_000015.3` in VEP, which is not what this holds.

use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};

use fastvep_cli::pipeline::{run_annotate, AnnotateConfig, PickFlags};
use tempfile::TempDir;

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("fixtures/vep_protein_version")
        .join(name)
}

/// fastVEP over the fixture's VCF and FASTA with the given GFF3 text, in a directory of its own
/// so that neither the transcript cache it writes beside the GFF3 nor the FASTA index lands in
/// the fixture.
fn annotate(gff3: &str, dir: &Path) -> String {
    let gff3_path = dir.join("genes.gff3");
    std::fs::write(&gff3_path, gff3).unwrap();
    let fasta = dir.join("genome.fa");
    std::fs::copy(fixture("genome.fa"), &fasta).unwrap();
    let out = dir.join("out.vcf");
    run_annotate(AnnotateConfig {
        input: fixture("in.vcf").to_string_lossy().into(),
        output: out.to_string_lossy().into(),
        gff3: vec![gff3_path.to_string_lossy().into()],
        fasta: Some(fasta.to_string_lossy().into()),
        output_format: "vcf".into(),
        pick: PickFlags::default(),
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
        columns: Default::default(),
    })
    .expect("annotation should succeed");
    std::fs::read_to_string(&out).unwrap()
}

/// record ID -> the HGVSp of every CSQ entry that has one.
fn hgvsp_by_record(annotated: &str) -> BTreeMap<String, BTreeSet<String>> {
    let hgvsp = fastvep_io::output::DEFAULT_CSQ_FIELDS
        .iter()
        .position(|f| *f == "HGVSp")
        .expect("HGVSp is a CSQ field");
    let mut rows: BTreeMap<String, BTreeSet<String>> = BTreeMap::new();
    for line in annotated.lines().filter(|l| !l.starts_with('#')) {
        let cols: Vec<&str> = line.split('\t').collect();
        let info = cols[7];
        let csq = info
            .split(';')
            .find_map(|kv| kv.strip_prefix("CSQ="))
            .unwrap_or_else(|| panic!("no CSQ on:\n{line}"));
        for entry in csq.split(',') {
            let fields: Vec<&str> = entry.split('|').collect();
            if !fields[hgvsp].is_empty() {
                rows.entry(cols[2].to_string())
                    .or_default()
                    .insert(fields[hgvsp].to_string());
            }
        }
    }
    rows
}

fn golden() -> BTreeMap<String, BTreeSet<String>> {
    let mut rows: BTreeMap<String, BTreeSet<String>> = BTreeMap::new();
    for line in std::fs::read_to_string(fixture("golden.tsv"))
        .unwrap()
        .lines()
        .filter(|l| !l.starts_with('#'))
    {
        let c: Vec<&str> = line.split('\t').collect();
        rows.entry(c[0].to_string())
            .or_default()
            .insert(c[2].to_string());
    }
    rows
}

#[test]
fn every_hgvsp_names_the_protein_version_vep_names() {
    let dir = TempDir::new().unwrap();
    let gff3 = std::fs::read_to_string(fixture("genes.gff3")).unwrap();
    let ours = hgvsp_by_record(&annotate(&gff3, dir.path()));
    let want = golden();
    // the fixture has what it says it has: versions 3, 1 and 12, and an id that carries its own
    let all: Vec<&String> = want.values().flatten().collect();
    for prefix in [
        "ENSP00000000011.3:",
        "ENSP00000000012.1:",
        "ENSP00000000014.12:",
        "NP_000015.2:",
    ] {
        assert!(
            all.iter().any(|h| h.starts_with(prefix)),
            "no {prefix} in VEP's answers"
        );
    }
    let wrong: Vec<String> = want
        .iter()
        .filter(|(record, hgvsp)| ours.get(*record) != Some(*hgvsp))
        .map(|(record, hgvsp)| {
            format!(
                "{record}: VEP {hgvsp:?}, here {:?}",
                ours.get(record).cloned().unwrap_or_default()
            )
        })
        .collect();
    assert!(
        wrong.is_empty(),
        "{} of {} records differ from VEP 116.3's HGVSp, first ten:\n{}",
        wrong.len(),
        want.len(),
        wrong
            .iter()
            .take(10)
            .cloned()
            .collect::<Vec<_>>()
            .join("\n")
    );
    assert_eq!(ours.len(), want.len(), "a record VEP writes no HGVSp for");
}

/// Both halves of the rule at once: with the `version` attributes taken off, every Ensembl
/// protein is `.1`, which is what VEP does with a record that has none.
#[test]
fn a_cds_without_a_version_is_version_one() {
    let dir = TempDir::new().unwrap();
    let gff = std::fs::read_to_string(fixture("genes.gff3")).unwrap();
    let gff: String = gff
        .lines()
        .map(|l| {
            if l.contains("\tCDS\t") {
                // `protein_id=ENSP...;version=N` becomes `protein_id=ENSP...`
                l.split(";version=").next().unwrap().to_string()
            } else {
                l.to_string()
            }
        })
        .collect::<Vec<_>>()
        .join("\n");
    let ours = hgvsp_by_record(&annotate(&(gff + "\n"), dir.path()));
    let ensembl: Vec<&String> = ours
        .values()
        .flatten()
        .filter(|h| h.starts_with("ENSP"))
        .collect();
    assert!(!ensembl.is_empty());
    assert!(
        ensembl.iter().all(|h| h.contains(".1:p.")),
        "a protein with no version is not `.1`: {ensembl:#?}"
    );
    // and an id that carries its version keeps that one, whatever the records say
    assert!(
        ours.values()
            .flatten()
            .any(|h| h.starts_with("NP_000015.2:p.")),
        "{ours:#?}"
    );
}
