//! A lone indel written with bases the REF and the ALT share beyond the anchor, against what
//! Ensembl VEP 116.3 writes for it (issue #139).
//!
//! `fixtures/vep_tab_lone_indels` is 1,500 records, each with one ALT that is a deletion, an
//! insertion or a replacement of one length by another, written as a caller that does not trim
//! might write it: a prefix both alleles share (the anchor and up to two more bases), the bases
//! that change, and a suffix they share (none to three bases). `regen.sh` there runs `vep` from the
//! 116.3 image in `--gff` mode over them and keeps what it wrote; nothing in the golden was worked
//! out here.
//!
//! VEP minimises such a record from both ends before it annotates it (`minimise_alleles`), so its
//! `Allele`, `Location`, `REF_ALLELE` and `VARIANT_CLASS` are those of the minimal pair, and
//! so is everything that follows. Rows are paired on the record and the transcript, because the
//! allele is what differs when the two do not agree.
//!
//! The columns left out are the ones a later change makes exact, each of them where the change
//! reaches a start codon, a stop codon or the end of a transcript (docs/VEP_DIVERGENCE.md, Part 3).

use std::collections::BTreeMap;
use std::io::Read;
use std::path::{Path, PathBuf};

use fastvep_cli::pipeline::{run_annotate, AnnotateConfig, ColumnFlags, PickFlags};
use flate2::read::GzDecoder;
use tempfile::TempDir;

/// Columns that are not compared yet.
const NOT_YET: &[&str] = &["Consequence", "IMPACT", "HGVSc", "HGVSp", "HGVS_OFFSET"];

fn here(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("fixtures/vep_tab_lone_indels")
        .join(name)
}

fn world(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("fixtures/vep_tab_world")
        .join(name)
}

type Row = BTreeMap<String, String>;
/// Rows keyed by the record and the transcript.
type Table = BTreeMap<(String, String), Row>;

fn table(text: &str) -> Table {
    let mut lines = text.lines().filter(|l| !l.starts_with("##"));
    let header: Vec<String> = lines
        .next()
        .expect("a header")
        .trim_start_matches('#')
        .split('\t')
        .map(String::from)
        .collect();
    let mut rows = Table::new();
    for line in lines {
        let cells: Vec<&str> = line.split('\t').collect();
        assert_eq!(cells.len(), header.len(), "ragged row {line:?}");
        let row: Row = header
            .iter()
            .cloned()
            .zip(cells.into_iter().map(String::from))
            .collect();
        let key = (row["Uploaded_variation"].clone(), row["Feature"].clone());
        assert!(
            rows.insert(key.clone(), row).is_none(),
            "two rows for {key:?}"
        );
    }
    rows
}

fn golden() -> Table {
    let mut text = String::new();
    GzDecoder::new(std::fs::File::open(here("golden_flags.tab.gz")).unwrap())
        .read_to_string(&mut text)
        .unwrap();
    table(&text)
}

fn ours() -> Table {
    let dir = TempDir::new().unwrap();
    // Copied, because the transcript cache a run writes goes beside its GFF3 and the index of
    // the FASTA beside it, and neither belongs in the source tree.
    let (gff3, fasta) = (dir.path().join("genes.gff3"), dir.path().join("genome.fa"));
    std::fs::copy(world("genes.gff3"), &gff3).unwrap();
    std::fs::copy(world("genome.fa"), &fasta).unwrap();
    let out = dir.path().join("out.tab");
    run_annotate(AnnotateConfig {
        input: here("in.vcf").to_string_lossy().into(),
        output: out.to_string_lossy().into(),
        gff3: vec![gff3.to_string_lossy().into()],
        fasta: Some(fasta.to_string_lossy().into()),
        output_format: "tab".into(),
        pick: PickFlags::default(),
        hgvs: true,
        columns: ColumnFlags {
            show_ref_allele: true,
            variant_class: true,
            symbol: true,
            biotype: true,
            canonical: true,
            numbers: true,
            hgvsg: true,
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
    })
    .expect("annotation should succeed");
    table(&std::fs::read_to_string(&out).unwrap())
}

#[test]
fn every_compared_cell_is_vep_s_for_every_row() {
    let (vep, fv) = (golden(), ours());
    let missing: Vec<_> = vep
        .keys()
        .filter(|k| !fv.contains_key(*k))
        .take(5)
        .collect();
    let extra: Vec<_> = fv
        .keys()
        .filter(|k| !vep.contains_key(*k))
        .take(5)
        .collect();
    assert!(missing.is_empty(), "rows only VEP writes: {missing:?}");
    assert!(extra.is_empty(), "rows only fastVEP writes: {extra:?}");

    let (mut wrong, mut compared) = (Vec::new(), 0usize);
    for (key, want) in &vep {
        let have = &fv[key];
        for (column, v) in want {
            if NOT_YET.contains(&column.as_str()) {
                continue;
            }
            compared += 1;
            if &have[column] != v {
                wrong.push(format!(
                    "{} {} {column}: VEP {v:?}, here {:?}",
                    key.0, key.1, have[column]
                ));
            }
        }
    }
    assert!(compared > 70_000, "only {compared} cells compared");
    assert!(
        wrong.is_empty(),
        "{} of {compared} cells differ from VEP's, first ten:\n{}",
        wrong.len(),
        wrong[..wrong.len().min(10)].join("\n")
    );
}

/// The fixture is about what it says: most of its records are not minimal, and VEP rewrites them.
#[test]
fn the_records_are_not_minimal_and_vep_minimises_them() {
    let input = std::fs::read_to_string(here("in.vcf")).unwrap();
    let mut shared_beyond_the_anchor = 0;
    let mut records = 0;
    for line in input.lines().filter(|l| !l.starts_with('#')) {
        let c: Vec<&str> = line.split('\t').collect();
        let (r, a) = (c[3].as_bytes(), c[4].as_bytes());
        records += 1;
        let front = r.iter().zip(a).take_while(|(x, y)| x == y).count();
        let back = r
            .iter()
            .rev()
            .zip(a.iter().rev())
            .take_while(|(x, y)| x == y)
            .count();
        if front > 1 || back > 0 {
            shared_beyond_the_anchor += 1;
        }
    }
    assert_eq!(records, 1500);
    assert!(
        shared_beyond_the_anchor > 1300,
        "{shared_beyond_the_anchor} of {records} share more than the anchor"
    );
    // the Allele column of VEP's rows is the minimal one: never the whole ALT of the record
    let vep = golden();
    let rewritten = vep
        .values()
        .filter(|row| {
            let alt = row["Uploaded_variation"].rsplit('/').next().unwrap();
            row["Allele"] != alt
        })
        .count();
    assert!(rewritten > 2500, "only {rewritten} rows were rewritten");
}
