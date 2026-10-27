//! A GFF3 with a damaged record, against what Ensembl VEP 116.3 does with it (issue #139).
//!
//! `fixtures/vep_gff_damaged` is a small RefSeq-style gene model (exons with no `rank`) with one
//! record damaged at a time, found by fuzzing the reader of fastVEP, and the tab file VEP writes
//! for each (`regen.sh` there runs `vep` from the 116.3 image in `--gff` mode). Nothing in a golden
//! was worked out here.
//!
//! - **A CDS phase of `.` or `-1`** is no phase, which VEP reads as a coding sequence that starts
//!   whole. It was read as `u64::MAX` bases to pad the first codon with: a run with a FASTA stopped
//!   in `capacity overflow`, one without in an overflow, and a release build wrapped.
//! - **A CDS phase that is not a phase** (`x`, `3`) stops VEP (`Bad value (x) for exon phase.
//!   Should only be -1,0,1,2`, exit 2, no rows), and stops this, naming the record.
//! - **Two exons of a transcript that overlap** leave the transcript out of VEP's run, with a
//!   warning (`Failed to add exon to transcript`). It stopped this one in a slice of the spliced
//!   sequence that started after it ended.
//! - **A reverse-strand CDS that ends in base 1** is read, as VEP reads it. Looking up its first
//!   codon in `transl_except` subtracted 2 from the base.

use std::collections::BTreeMap;
use std::io::Read;
use std::path::{Path, PathBuf};

use fastvep_cli::pipeline::{run_annotate, AnnotateConfig, ColumnFlags, PickFlags};
use flate2::read::GzDecoder;
use tempfile::TempDir;

fn here(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("fixtures/vep_gff_damaged")
        .join(name)
}

type Row = BTreeMap<String, String>;
/// Rows keyed by what names them: the record, the allele and the transcript.
type Table = BTreeMap<(String, String, String), Vec<Row>>;

fn table(text: &str) -> Table {
    let mut lines = text.lines().filter(|l| !l.starts_with("##"));
    let header: Vec<String> = lines
        .next()
        .expect("a header")
        .trim_start_matches('#')
        .split('\t')
        .map(String::from)
        .collect();
    let mut rows: Table = BTreeMap::new();
    for line in lines {
        let cells: Vec<&str> = line.split('\t').collect();
        assert_eq!(cells.len(), header.len(), "ragged row {line:?}");
        let row: Row = header
            .iter()
            .cloned()
            .zip(cells.into_iter().map(String::from))
            .collect();
        let key = (
            row["Uploaded_variation"].clone(),
            row["Allele"].clone(),
            row["Feature"].clone(),
        );
        rows.entry(key).or_default().push(row);
    }
    rows
}

fn golden(case: &str) -> Table {
    let mut text = String::new();
    GzDecoder::new(std::fs::File::open(here(&format!("golden_{case}.tab.gz"))).unwrap())
        .read_to_string(&mut text)
        .unwrap();
    table(&text)
}

/// What fastVEP does with the GFF3 `gff3` (a file of the fixture, or one this test wrote into
/// `dir`), the same flags as `regen.sh` gives VEP, with the FASTA or without one.
fn run(dir: &Path, gff3: &Path, vcf: &str, fasta: Option<&str>) -> anyhow::Result<(String, Table)> {
    // Copied, because the transcript cache a run writes goes beside its GFF3 and the index of the
    // FASTA beside it, and neither belongs in the source tree.
    let copy = |from: &Path, name: &str| {
        let to = dir.join(name);
        std::fs::copy(from, &to).unwrap();
        to
    };
    let gff3 = copy(gff3, "genes.gff3");
    let fasta = fasta.map(|f| copy(&here(f), "genome.fa"));
    let out = dir.join("out.tab");
    run_annotate(AnnotateConfig {
        input: here(vcf).to_string_lossy().into(),
        output: out.to_string_lossy().into(),
        gff3: vec![gff3.to_string_lossy().into()],
        fasta: fasta.map(|f| f.to_string_lossy().into()),
        output_format: "tab".into(),
        pick: PickFlags::default(),
        hgvs: true,
        columns: ColumnFlags {
            symbol: true,
            biotype: true,
            canonical: true,
            numbers: true,
            ..Default::default()
        },
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
        compat: Default::default(),
        exclude: Default::default(),
    })?;
    let text = std::fs::read_to_string(&out).unwrap();
    let rows = table(&text);
    Ok((text, rows))
}

/// The columns the reading of a gene model decides.
const COMPARED: &[&str] = &[
    "Location",
    "Consequence",
    "IMPACT",
    "Amino_acids",
    "Codons",
    "cDNA_position",
    "CDS_position",
    "Protein_position",
    "EXON",
    "INTRON",
    "DISTANCE",
    "STRAND",
    "HGVSc",
    "HGVSp",
];

/// Every row VEP wrote, fastVEP wrote, and no other, with every compared cell VEP's.
fn assert_is_vep(case: &str, gff3: &str, vcf: &str, fasta: &str) -> Table {
    let dir = TempDir::new().unwrap();
    let (_, ours) = run(dir.path(), &here(gff3), vcf, Some(fasta))
        .unwrap_or_else(|e| panic!("{case}: annotation should succeed: {e:#}"));
    let vep = golden(case);
    assert_eq!(
        vep.keys().collect::<Vec<_>>(),
        ours.keys().collect::<Vec<_>>(),
        "{case}: the rows differ"
    );
    let mut cells = 0;
    for (key, vep_rows) in &vep {
        let ours_rows = &ours[key];
        assert_eq!(vep_rows.len(), ours_rows.len(), "{case} {key:?}");
        for (v, o) in vep_rows.iter().zip(ours_rows) {
            for col in COMPARED {
                cells += 1;
                assert_eq!(v[*col], o[*col], "{case} {key:?} {col}");
            }
        }
    }
    assert!(cells > 100, "{case}: only {cells} cells compared");
    vep
}

#[test]
fn a_cds_with_no_phase_is_read_as_vep_reads_it() {
    for (case, gff3) in [
        ("phase_dot", "phase_dot.gff3"),
        ("phase_minus1", "phase_minus1.gff3"),
    ] {
        let vep = assert_is_vep(case, gff3, "phase.vcf", "genome.fa");
        // the golden is about what it says: the first codon of NM_090001.3 is where the CDS starts
        // whole, so a change in it is a start loss and one in the second codon a missense
        let consequences: Vec<&str> = vep
            .iter()
            .filter(|(k, _)| k.2 == "NM_090001.3")
            .flat_map(|(_, rows)| rows.iter().map(|r| r["Consequence"].as_str()))
            .collect();
        assert!(
            consequences.contains(&"start_lost"),
            "{case}: {consequences:?}"
        );
        assert!(
            consequences.contains(&"stop_gained"),
            "{case}: {consequences:?}"
        );
    }
}

/// With no FASTA the run has no sequence to build, and a phase of `u64::MAX` bases was added to a
/// CDS position instead.
#[test]
fn a_cds_with_no_phase_is_read_without_a_fasta_too() {
    let dir = TempDir::new().unwrap();
    let (_, rows) = run(dir.path(), &here("phase_dot.gff3"), "phase.vcf", None)
        .expect("annotation should succeed without a FASTA");
    let coding = |feature: &str| {
        rows.iter()
            .filter(|(k, _)| k.2 == feature)
            .flat_map(|(_, rows)| rows.iter())
            .filter(|r| r["CDS_position"] != "-")
            .count()
    };
    assert!(
        coding("NM_090001.3") > 5,
        "no CDS positions for the transcript"
    );
}

#[test]
fn a_transcript_whose_exons_overlap_is_left_out_as_vep_leaves_it_out() {
    let vep = assert_is_vep("overlap", "overlap.gff3", "overlap.vcf", "genome.fa");
    assert!(
        vep.keys().all(|k| k.2 != "NM_090002.2"),
        "VEP has a row for the transcript with overlapping exons"
    );
    assert!(
        vep.keys().any(|k| k.2 == "NM_090001.3"),
        "the other transcript is annotated"
    );
}

#[test]
fn a_cds_that_ends_in_base_one_is_read_as_vep_reads_it() {
    let vep = assert_is_vep("base1", "base1.gff3", "base1.vcf", "base1.fa");
    assert!(vep.keys().any(|k| k.2 == "NM_1.1"));
}

/// VEP stops on these (exit 2, no row), so there is no golden: fastVEP stops too, and says which
/// record and what a phase is, instead of panicking in a thread of the build.
#[test]
fn a_cds_phase_that_is_not_a_phase_stops_the_run_where_vep_stops() {
    let gff3 = std::fs::read_to_string(here("phase_dot.gff3")).unwrap();
    // the first CDS record of NM_090001.3, which has the phase `.`
    let record = "CDS\t1100\t1299\t.\t+\t.\t";
    assert!(gff3.contains(record));
    for bad in ["x", "3"] {
        let dir = TempDir::new().unwrap();
        let damaged = dir.path().join("damaged.gff3");
        std::fs::write(
            &damaged,
            gff3.replace(record, &format!("CDS\t1100\t1299\t.\t+\t{bad}\t")),
        )
        .unwrap();
        let err = run(dir.path(), &damaged, "phase.vcf", Some("genome.fa"))
            .expect_err("a CDS phase that is not a phase must stop the run");
        let msg = format!("{err:#}");
        assert!(
            msg.contains("chrR:1100-1299") && msg.contains("rna-NM_090001.3") && msg.contains(bad),
            "{bad}: {msg}"
        );
    }
}
