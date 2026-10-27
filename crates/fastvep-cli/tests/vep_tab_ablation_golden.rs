//! The tab output for a deletion that covers a whole transcript, against what Ensembl VEP 116.3
//! writes for the same change (issue #139).
//!
//! Ensembl's `transcript_ablation` is `feature_ablation`: a deletion (an alternate shorter than the
//! reference) whose span covers the whole transcript (`complete_overlap_feature`). It is the only
//! term such a change gets, because every other term's `include` excludes a complete overlap. A
//! predictor that asks the regions the span reaches one by one says `stop_retained_variant,
//! 5_prime_UTR_variant,3_prime_UTR_variant` for it, which is LOW, and a small non-coding
//! transcript was MODIFIER.
//!
//! `fixtures/vep_tab_ablation` is a world of eight transcripts of 21 to 130 bases, both strands,
//! coding and not, in one exon and in three, and 396 single-ALT records: a deletion at every
//! combination of start from 12 bases before the transcript to 3 bases in and end from 3 bases
//! before its last base to 12 after (so spans that cover it exactly, by one base, by many, and
//! that stop a base short are all there), and over each transcript a replacement by something
//! longer, one as long, one shorter, an insertion, and substitutions. `regen.sh` there runs `vep`
//! from the 116.3 image in `--gff` mode over them and keeps what it wrote. Nothing in the golden
//! was worked out here.
//!
//! Every column that is about the consequence and the position is compared, and no cell differs.

use std::collections::BTreeMap;
use std::io::Read;
use std::path::{Path, PathBuf};

use fastvep_cli::pipeline::{run_annotate, AnnotateConfig, ColumnFlags, PickFlags};
use flate2::read::GzDecoder;
use tempfile::TempDir;

/// The columns compared, cell by cell, on every row. `HGVSc`, `HGVS_OFFSET` and `HGVSp` are not:
/// a change that runs past an end of a transcript, and is then moved by the 3'-rule, is the subject
/// of `vep_tab_edges_golden.rs`, where the cases are listed, and in this world's deletions it is
/// most of what there is. (Ensembl writes `n.1_21del` for a deletion of all of a 21-base
/// transcript that the rule would carry a base past its end, and fastVEP nothing.)
const COMPARED: &[&str] = &[
    "Location",
    "cDNA_position",
    "CDS_position",
    "Protein_position",
    "EXON",
    "INTRON",
    "Consequence",
    "IMPACT",
    "Amino_acids",
    "Codons",
    "HGVSg",
];

fn here(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("fixtures/vep_tab_ablation")
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
    std::fs::copy(here("genes.gff3"), &gff3).unwrap();
    std::fs::copy(here("genome.fa"), &fasta).unwrap();
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
            ..Default::default()
        },
        distance: 1000,
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
    assert!(compared > 15_000, "only {compared} cells compared");
    assert!(
        wrong.is_empty(),
        "{} of {compared} cells differ from VEP, e.g.\n{}",
        wrong.len(),
        wrong[..wrong.len().min(12)].join("\n")
    );
}

/// The fixture shows what it is about, from VEP's golden alone: a deletion that covers a
/// transcript is `transcript_ablation` and nothing else, in both kinds of transcript and on both
/// strands, and a span that stops a base short is not.
#[test]
fn the_fixture_shows_what_vep_does() {
    let vep = golden();
    let rows = || vep.values().flatten();
    let ablated: Vec<&Row> = rows()
        .filter(|r| r["Consequence"].contains("transcript_ablation"))
        .collect();
    assert!(ablated.len() > 80, "{} ablations", ablated.len());
    assert!(ablated
        .iter()
        .all(|r| r["Consequence"] == "transcript_ablation" && r["IMPACT"] == "HIGH"));
    // a deletion, or a replacement by something shorter, which is `indel` and is a deletion to
    // Ensembl's `deletion` all the same
    assert!(ablated
        .iter()
        .all(|r| matches!(r["VARIANT_CLASS"].as_str(), "deletion" | "indel")));
    assert!(ablated.iter().any(|r| r["VARIANT_CLASS"] == "indel"));
    let has = |f: &dyn Fn(&Row) -> bool| ablated.iter().any(|r| f(r));
    assert!(has(&|r| r["STRAND"] == "1") && has(&|r| r["STRAND"] == "-1"));
    assert!(has(&|r| r["BIOTYPE"] == "lncRNA") && has(&|r| r["BIOTYPE"] == "protein_coding"));
    // a deletion that stops a base short of an end of the transcript is an exonic change
    let short = rows()
        .filter(|r| r["VARIANT_CLASS"] == "deletion" && r["EXON"] != "-")
        .filter(|r| !r["Consequence"].contains("transcript_ablation"))
        .count();
    assert!(short > 30, "{short} deletions that are not ablations");
}
