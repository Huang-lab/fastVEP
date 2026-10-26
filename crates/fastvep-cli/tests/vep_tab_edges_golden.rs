//! The tab output for a change at the edge of an exon, of the coding sequence or of a
//! transcript, against what Ensembl VEP 116.3 writes for the same change (issue #139).
//!
//! `fixtures/vep_tab_edges` is 2,271 single-ALT records, a sample of 17,672: every position
//! within ten bases of every exon and CDS boundary of the world in `fixtures/vep_tab_world`,
//! on both strands, in coding and non-coding transcripts: SNVs, MNVs, deletions of 1 to 30
//! bases, insertions (random and duplications), length-changing delins, and deletions of a
//! whole exon or of a whole intron. `regen.sh`
//! there runs `vep` from the 116.3 image in `--gff` mode over them and keeps what it wrote.
//! Nothing in the golden was worked out here.
//!
//! The columns compared are the ones that are exact; a column joins `COMPARED` when the
//! change that makes it exact lands. What is left is in `docs/VEP_DIVERGENCE.md` (Part 3),
//! measured on the full 17,672 records.

use std::collections::BTreeMap;
use std::io::Read;
use std::path::{Path, PathBuf};

use fastvep_cli::pipeline::{run_annotate, AnnotateConfig, ColumnFlags, PickFlags};
use flate2::read::GzDecoder;
use tempfile::TempDir;

/// The columns whose every cell, on every row, is VEP's.
const COMPARED: &[&str] = &[
    "Location",
    "cDNA_position",
    "CDS_position",
    "Protein_position",
    "EXON",
    "INTRON",
];

fn here(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("fixtures/vep_tab_edges")
        .join(name)
}

fn world(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("fixtures/vep_tab_world")
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

fn golden() -> Table {
    let mut text = String::new();
    GzDecoder::new(std::fs::File::open(here("golden_flags.tab.gz")).unwrap())
        .read_to_string(&mut text)
        .unwrap();
    table(&text)
}

fn ours() -> Table {
    let dir = TempDir::new().unwrap();
    // Copied, because the transcript cache a run writes goes beside its GFF3 and the index
    // of the FASTA beside it, and neither belongs in the source tree.
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

/// Every row VEP wrote, fastVEP wrote, and no other.
#[test]
fn both_tools_write_the_same_rows() {
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
}

/// The cells of the compared columns, row by row.
#[test]
fn the_compared_columns_are_vep_s_for_every_row() {
    let (vep, fv) = (golden(), ours());
    let mut wrong = Vec::new();
    let mut compared = 0usize;
    for (key, vep_rows) in &vep {
        let fv_rows = &fv[key];
        assert_eq!(vep_rows.len(), fv_rows.len(), "{key:?}");
        for (v, f) in vep_rows.iter().zip(fv_rows) {
            for col in COMPARED {
                compared += 1;
                if v[*col] != f[*col] {
                    wrong.push(format!(
                        "{} {} {col}: VEP {:?}, fastVEP {:?}",
                        key.0, key.2, v[*col], f[*col]
                    ));
                }
            }
        }
    }
    assert!(compared > 20_000, "only {compared} cells compared");
    assert!(
        wrong.is_empty(),
        "{} of {compared} cells differ from VEP, e.g.\n{}",
        wrong.len(),
        wrong[..wrong.len().min(12)].join("\n")
    );
}

/// The fixture exercises what it is for: both kinds of unknown end, an insertion beside an
/// exon, one beside the coding sequence, and a change over two exons.
#[test]
fn the_fixture_covers_the_shapes_the_columns_are_about() {
    let vep = golden();
    let cells = |col: &str| -> Vec<&str> {
        vep.values()
            .flatten()
            .map(|r| r[col].as_str())
            .collect::<Vec<_>>()
    };
    let cdna = cells("cDNA_position");
    assert!(cdna.iter().any(|c| c.starts_with("?-")), "no `?-N`");
    assert!(cdna.iter().any(|c| c.ends_with("-?")), "no `N-?`");
    assert!(
        cdna.iter().any(|c| {
            c.split_once('-').is_some_and(|(a, b)| {
                a.parse::<u64>()
                    .is_ok_and(|a| b.parse::<u64>() == Ok(a + 1))
            })
        }),
        "no insertion pair"
    );
    assert!(
        cells("CDS_position").contains(&"-") && cells("EXON").contains(&"-"),
        "no row without a CDS position or an exon"
    );
    assert!(
        cells("EXON")
            .iter()
            .any(|c| c.contains('-') && c.contains('/') && !c.starts_with('-')),
        "no change over several exons"
    );
    assert!(cells("INTRON").iter().any(|c| *c != "-"), "no INTRON");
}
