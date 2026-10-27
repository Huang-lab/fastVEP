//! The tab output for a change at the start codon, against what Ensembl VEP 116.3 writes for the
//! same change (issue #139).
//!
//! `fixtures/vep_tab_initiator` is a world of ten transcripts whose start codons are built to be
//! awkward: a second codon that repeats the last base of `ATG` (so that a deletion of codon 2 has a
//! leftmost form that touches the start codon), transcripts that begin at their initiator and so
//! have no 5' UTR, an initiator repeated (`ATG ATG`), a short UTR, and a start codon split by an
//! intron, on both strands; and 6,637 single-ALT records, every substitution, deletion of 1 to 6, 9
//! and 12 bases, insertion of 1 to 6 bases and replacement of two and three bases at every position
//! from 12 bases before the coding start to 14 after it. `regen.sh` there runs `vep` from the 116.3
//! image in `--gff` mode over them and keeps what it wrote. Nothing in the golden was worked out
//! here.
//!
//! What is compared is what the start codon decides: `Consequence`, `IMPACT`, `Amino_acids`,
//! `Codons`, the positions, `EXON`, `INTRON`, and `HGVSp` of the rows that are a start loss in
//! either tool. (`HGVSc` and `HGVSp` elsewhere, and the 3'-rule, are the subject of the other
//! `vep_tab_*_golden.rs` files; this world's repeats would only raise their listed differences
//! again.) Every difference is one of four documented classes, and each row that is one is checked
//! to have the shape the class says:
//!
//! - **A start loss that is not one** (`docs/VEP_DIVERGENCE.md`, divergence 14). VEP calls the
//!   initiator lost where the edited protein still begins with its Met: the loss of codon 2
//!   (`MA/M`, `start_lost,inframe_deletion,start_retained_variant`), a replacement over the start
//!   codon that leaves it as it was and changes codon 2, and any change on a transcript with no 5'
//!   UTR that leaves the ATG where it was, because the test that asks "is the ATG still at the
//!   coding start" is only made `if ($utr)`. fastVEP leaves out `start_lost` (and says
//!   `start_retained_variant` where VEP's test could not).
//! - **A start loss written as a range, or not written** (divergence 15). `p.Met1?` is what HGVS
//!   has for a change that loses the initiator, and VEP writes it for most of them; for a window of
//!   more than one residue it writes the residue names of the whole window sharing one position
//!   (`p.MetAla1_?2`), for an insertion in front of the initiator a range that runs backwards and
//!   names the last residue of the protein (`p.Pro1_?0`), and for a deletion that begins in the 5'
//!   UTR and runs into the start codon, which has no protein position to start from, nothing at
//!   all. fastVEP writes `p.Met1?` for every one.
//! - **A start loss that Ensembl does not see on a transcript that begins at its initiator**
//!   (divergence 14, the second way). The deletion of the whole ATG, and of the codons after it, is an
//!   `inframe_deletion` and `p.Met1del` there, because the test that asks whether the codon at the
//!   coding start is still an ATG is made `if ($utr)`. fastVEP asks it of a transcript with no UTR as
//!   it does of one that has one: `start_lost` beside the `inframe_deletion`, and `p.Met1?`.
//! - **A frameshift and a start loss for a protein that has neither** (divergence 14, the third
//!   way). An edit that touches the start codon and leaves the whole coding sequence behind it,
//!   with no ATG where the initiator was, leaves the protein as it was: VEP writes `start_lost`
//!   beside its own `start_retained_variant`, and `frameshift_variant` and `p.Met1?` as well, where
//!   fastVEP writes `start_retained_variant` and no notation.

use std::collections::BTreeMap;
use std::io::Read;
use std::path::{Path, PathBuf};

use fastvep_cli::pipeline::{run_annotate, AnnotateConfig, ColumnFlags, PickFlags};
use flate2::read::GzDecoder;
use tempfile::TempDir;

fn here(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("fixtures/vep_tab_initiator")
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
        distance: 0,
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

/// The transcripts of the world that begin at their initiator.
const NO_UTR: &[&str] = &["ENST00000000203", "ENST00000000204", "ENST00000000205"];

fn terms(row: &Row) -> Vec<&str> {
    row["Consequence"].split(',').collect()
}

fn has(row: &Row, term: &str) -> bool {
    terms(row).contains(&term)
}

/// The window of a row opens with its Met on both sides: `MA/M`, `M/MX`, `MAL/M`.
fn window_keeps_its_met(row: &Row) -> bool {
    row["Amino_acids"]
        .split_once('/')
        .is_some_and(|(r, a)| r.starts_with('M') && a.starts_with('M'))
}

#[derive(Debug, PartialEq, Eq, Clone, Copy, PartialOrd, Ord)]
enum Class {
    /// Divergence 14: `start_lost` where the window still opens with its Met.
    NotAStartLoss,
    /// Divergence 14, the second half: a transcript with no 5' UTR that keeps its ATG is
    /// `start_retained_variant`, as the same edit is on a transcript with one.
    NoUtr,
    /// Divergence 15: `p.Met1?`, not a range.
    StartLossNotation,
    /// Divergence 15, the second half: `p.Met1?` where VEP writes no HGVSp at all, for a deletion
    /// that begins in the 5' UTR and runs into the start codon.
    StartLossUnwritten,
    /// Divergence 14, the third way: the whole coding sequence is behind the edit and no ATG is at
    /// the old coding start, so the protein is the protein.
    CodingSequenceKept,
    /// Divergence 14, the second way once more: the whole ATG deleted in frame on a transcript with
    /// no 5' UTR, which VEP's test is not made for.
    NoUtrInitiatorDeleted,
}

/// Why a row that is not VEP's is not, or `None` when nothing documented explains it.
fn explain(v: &Row, f: &Row, differing: &[&str]) -> Option<Class> {
    let hgvsp = |r: &Row| r["HGVSp"].clone();
    // both call it a start loss and only the notation differs
    if has(v, "start_lost") && has(f, "start_lost") && differing == ["HGVSp"] {
        let protein = f["HGVSp"].split(':').next()?;
        let met1 = hgvsp(f) == format!("{protein}:p.Met1?");
        if v["HGVSp"] == "-" {
            return met1.then_some(Class::StartLossUnwritten);
        }
        let range = v["HGVSp"].contains("_?");
        return (range && met1 && v["HGVSp"].split(':').next() == Some(protein))
            .then_some(Class::StartLossNotation);
    }
    // the initiator deleted in frame, and the codons after it, on a transcript that begins at it:
    // VEP's terms and fastVEP's with `start_lost`, and `p.Met1?` for `p.Met1_Ala2del`
    if NO_UTR.contains(&v["Feature"].as_str())
        && !has(v, "start_lost")
        && has(f, "start_lost")
        && has(v, "inframe_deletion")
        && differing
            .iter()
            .all(|c| ["Consequence", "IMPACT", "HGVSp"].contains(c))
        && v["Amino_acids"]
            .split_once('/')
            .is_some_and(|(r, a)| r.starts_with('M') && a == "-")
        && f["HGVSp"].ends_with(":p.Met1?")
        && v["HGVSp"].contains(":p.Met1")
    {
        let mut want = terms(v);
        want.push("start_lost");
        want.sort_unstable();
        let mut got = terms(f);
        got.sort_unstable();
        if got == want {
            return Some(Class::NoUtrInitiatorDeleted);
        }
    }
    // VEP's start loss and frameshift for a protein that has neither: its own pair of
    // `start_lost` and `start_retained_variant` says its tests did not agree, and what is left is
    // `start_retained_variant`, no term about the protein changing, and no notation. No ATG
    // is where the initiator was, or the edit would be an insertion of one (`ATG ATG`).
    if has(v, "start_lost")
        && has(v, "start_retained_variant")
        && !has(f, "start_lost")
        && differing
            .iter()
            .all(|c| ["Consequence", "IMPACT", "HGVSp"].contains(c))
        && f["HGVSp"] == "-"
        && !f["Codons"]
            .split_once('/')
            .is_some_and(|(_, a)| a.to_ascii_uppercase().starts_with("ATG"))
    {
        let mut want: Vec<&str> = terms(v)
            .into_iter()
            .filter(|t| {
                !matches!(
                    *t,
                    "start_lost"
                        | "frameshift_variant"
                        | "inframe_insertion"
                        | "inframe_deletion"
                        | "protein_altering_variant"
                )
            })
            .collect();
        let mut got = terms(f);
        want.sort_unstable();
        got.sort_unstable();
        if got == want {
            return Some(Class::CodingSequenceKept);
        }
    }
    // VEP's start loss that is not one: the terms are VEP's less `start_lost`, and
    // `start_retained_variant` where VEP's test could not say it
    if has(v, "start_lost") && !has(f, "start_lost") {
        let allowed = ["Consequence", "IMPACT", "HGVSp"];
        let mut want: Vec<&str> = terms(v)
            .into_iter()
            .filter(|t| *t != "start_lost")
            .collect();
        if !want.contains(&"start_retained_variant") {
            want.push("start_retained_variant");
        }
        let (mut got, mut want_sorted) = (terms(f), want.clone());
        got.sort_unstable();
        want_sorted.sort_unstable();
        // the terms may also keep the in-frame term `start_lost` was hiding
        let only_hidden_terms_added = got.iter().all(|t| {
            want_sorted.contains(t)
                || matches!(
                    *t,
                    "inframe_insertion" | "inframe_deletion" | "missense_variant"
                )
        });
        // the window shows the Met, or it is the placeholder `X` of a frameshift on a transcript
        // that has no UTR, where the test that would say the ATG stayed was never made
        // or the coding sequence is all there behind inserted bases (`M/NM`): the edit is in front
        // of the initiator in all but spelling, and Ensembl's `_inv_start_altered` asks the place
        let cds_behind = has(v, "start_retained_variant")
            && v["Amino_acids"]
                .split_once('/')
                .is_some_and(|(r, a)| r.starts_with('M') && a.ends_with(r));
        let met_shown = window_keeps_its_met(v)
            || cds_behind
            || (NO_UTR.contains(&v["Feature"].as_str())
                && v["Amino_acids"]
                    .split_once('/')
                    .is_some_and(|(r, a)| r.starts_with('M') && a == "X"));
        return (differing.iter().all(|c| allowed.contains(c))
            && met_shown
            && !has(f, "start_lost")
            && only_hidden_terms_added)
            .then_some(Class::NotAStartLoss);
    }
    // an in-frame term on a transcript with no UTR, where VEP's test could not say the ATG stayed
    if differing == ["Consequence"]
        && NO_UTR.contains(&v["Feature"].as_str())
        && has(f, "start_retained_variant")
        && !has(v, "start_retained_variant")
    {
        let mut a = terms(v);
        a.push("start_retained_variant");
        a.sort_unstable();
        let mut b = terms(f);
        b.sort_unstable();
        return (a == b).then_some(Class::NoUtr);
    }
    None
}

/// The columns the start codon decides.
const CORE: &[&str] = &[
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
];

#[test]
fn the_start_codon_columns_are_vep_s_but_for_the_documented_classes() {
    let (vep, fv) = (golden(), ours());
    let mut classes: BTreeMap<Class, usize> = BTreeMap::new();
    let mut unexplained = Vec::new();
    let mut compared = 0usize;
    for (key, vep_rows) in &vep {
        let fv_rows = &fv[key];
        assert_eq!(vep_rows.len(), fv_rows.len(), "{key:?}");
        for (v, f) in vep_rows.iter().zip(fv_rows) {
            let start_loss = has(v, "start_lost") || has(f, "start_lost");
            let mut differing: Vec<&str> = Vec::new();
            for col in CORE.iter().chain(start_loss.then_some(&"HGVSp")) {
                compared += 1;
                if v[*col] != f[*col] {
                    differing.push(col);
                }
            }
            if differing.is_empty() {
                continue;
            }
            match explain(v, f, &differing) {
                Some(class) => *classes.entry(class).or_default() += 1,
                None => unexplained.push(format!(
                    "{} {} {differing:?}: VEP {:?} {:?} {:?}, fastVEP {:?} {:?} {:?}",
                    key.0,
                    key.2,
                    v["Consequence"],
                    v["Amino_acids"],
                    v["HGVSp"],
                    f["Consequence"],
                    f["Amino_acids"],
                    f["HGVSp"],
                )),
            }
        }
    }
    assert!(compared > 60_000, "only {compared} cells compared");
    assert!(
        unexplained.is_empty(),
        "{} rows differ from VEP in a way nothing documents, e.g.\n{}",
        unexplained.len(),
        unexplained[..unexplained.len().min(12)].join("\n")
    );
    // each class is there, so that a class that stopped happening is noticed too
    assert!(
        classes.get(&Class::StartLossNotation) > Some(&100),
        "{classes:?}"
    );
    assert!(
        classes.get(&Class::StartLossUnwritten) > Some(&100),
        "{classes:?}"
    );
    assert!(
        classes.get(&Class::NotAStartLoss) > Some(&20),
        "{classes:?}"
    );
    assert!(classes.get(&Class::NoUtr) > Some(&2), "{classes:?}");
    assert!(
        classes.get(&Class::CodingSequenceKept) > Some(&40),
        "{classes:?}"
    );
    assert!(
        classes.get(&Class::NoUtrInitiatorDeleted) > Some(&8),
        "{classes:?}"
    );
}

/// What the golden shows of VEP, alone, so that the differences above are about VEP and not about
/// the fixture: it calls the initiator lost and retained on one row (the two terms are opposites),
/// it writes a ranges for a start loss, and on a transcript with no UTR it calls a frameshift that
/// keeps the ATG a start loss where the same edit on a transcript with one is not.
#[test]
fn the_fixture_shows_what_vep_does() {
    let vep = golden();
    let rows = || vep.values().flatten();
    let both = rows()
        .filter(|r| {
            has(r, "start_lost") && has(r, "start_retained_variant") && window_keeps_its_met(r)
        })
        .count();
    assert!(
        both > 10,
        "{both} rows that are `start_lost` and `start_retained_variant` with a Met kept"
    );
    let ranges = rows().filter(|r| r["HGVSp"].contains("_?")).count();
    assert!(ranges > 100, "{ranges} start losses written as a range");
    let plain = rows().filter(|r| r["HGVSp"].ends_with(":p.Met1?")).count();
    assert!(plain > 50, "{plain} start losses written `p.Met1?`");
    // and for a protein that has no change it writes the frameshift and the loss as well
    let frameshift_and_loss_with_the_cds_behind = rows()
        .filter(|r| {
            has(r, "frameshift_variant")
                && has(r, "start_lost")
                && has(r, "start_retained_variant")
                && r["Amino_acids"]
                    .split_once('/')
                    .is_some_and(|(r, a)| r.starts_with('M') && a.ends_with('X'))
        })
        .count();
    assert!(
        frameshift_and_loss_with_the_cds_behind > 20,
        "{frameshift_and_loss_with_the_cds_behind}"
    );
    let frameshift_on_no_utr = rows()
        .filter(|r| NO_UTR.contains(&r["Feature"].as_str()))
        .filter(|r| has(r, "frameshift_variant") && has(r, "start_lost") && window_keeps_its_met(r))
        .count();
    assert!(frameshift_on_no_utr > 5, "{frameshift_on_no_utr}");
    let on = |strand: &str| rows().any(|r| r["STRAND"] == strand);
    assert!(on("1") && on("-1"));
}
