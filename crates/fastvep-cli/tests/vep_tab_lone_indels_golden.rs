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
//! Every column is compared. The cells where fastVEP is not VEP are listed in `KNOWN`, each
//! asserted to still differ (docs/VEP_DIVERGENCE.md, Part 3).

use std::collections::BTreeMap;
use std::io::Read;
use std::path::{Path, PathBuf};

use fastvep_cli::pipeline::{run_annotate, AnnotateConfig, ColumnFlags, PickFlags};
use flate2::read::GzDecoder;
use tempfile::TempDir;

/// Cells where stock VEP 116.3 is wrong, in a way `docs/VEP_DIVERGENCE.md` (Part 3) explains:
/// `(record, transcript, column, what fastVEP writes)`.
///
/// Eight are the `--gff` reader's UTR of a reverse-strand transcript, which puts a frameshift's
/// stop 12 residues short of the one the edited transcript has (the same rows, with that
/// accessor corrected, are fastVEP's); the ninth is a change that moves the frame at the stop
/// codon, where `frameshift` answers no and VEP's notation is `p.Ala200_Ter201delinsTer`.
const KNOWN: &[(&str, &str, &str, &str)] = &[
    (
        "chr10_10097_AGTTAAGCAA/AGTGGAAGCCA",
        "ENST00000000003",
        "HGVSp",
        "ENSP00000000003.1:p.Ter201SerextTer35",
    ),
    (
        "chr10_10098_GTTAAGCA/GGCA",
        "ENST00000000003",
        "HGVSp",
        "ENSP00000000003.1:p.Ter201%3D",
    ),
    (
        "chr10_10099_TTAAGCAAGA/TTAGATCTGCA",
        "ENST00000000003",
        "HGVSp",
        "ENSP00000000003.1:p.Leu199AlafsTer38",
    ),
    (
        "chr10_10100_TAAGCAAG/TATAAG",
        "ENST00000000003",
        "HGVSp",
        "ENSP00000000003.1:p.Ala200IlefsTer36",
    ),
    (
        "chr10_10102_AGCAAGACCCAAG/AGCAG",
        "ENST00000000003",
        "HGVSp",
        "ENSP00000000003.1:p.Leu197CysfsTer37",
    ),
    (
        "chr10_10108_ACCCAAGGTCCA/ACCA",
        "ENST00000000003",
        "HGVSp",
        "ENSP00000000003.1:p.Thr196SerfsTer38",
    ),
    (
        "chr10_10110_CCAAGG/CCACAGG",
        "ENST00000000003",
        "HGVSp",
        "ENSP00000000003.1:p.Leu197CysfsTer40",
    ),
    (
        "chr10_10111_CAAGGTCC/CAC",
        "ENST00000000003",
        "HGVSp",
        "ENSP00000000003.1:p.Trp195CysfsTer40",
    ),
    (
        "chr10_10113_AGGT/AGCAGGTATGT",
        "ENST00000000003",
        "HGVSp",
        "ENSP00000000003.1:p.Leu197TyrfsTer42",
    ),
    // a change that loses the initiator is `p.Met1?`, where VEP writes a range or two residue
    // names (`p.Met1_?2`, `p.MetCys1_?2`): not HGVS (divergence 14)
    (
        "chr10_1098_AGATGCAC/AGA",
        "ENST00000000001",
        "HGVSp",
        "ENSP00000000001.1:p.Met1?",
    ),
    (
        "chr10_1098_AGATGCAC/AGA",
        "ENST00000000002",
        "HGVSp",
        "ENSP00000000002.1:p.Met1?",
    ),
    (
        "chr10_1098_AGATGCAC/AGAC",
        "ENST00000000001",
        "HGVSp",
        "ENSP00000000001.1:p.Met1?",
    ),
    (
        "chr10_1098_AGATGCAC/AGAC",
        "ENST00000000002",
        "HGVSp",
        "ENSP00000000002.1:p.Met1?",
    ),
    (
        "chr10_11293_ATGCACATC/ATC",
        "ENST00000000003",
        "HGVSp",
        "ENSP00000000003.1:p.Met1?",
    ),
    // a start loss that VEP leaves without a notation, and an insertion in front of the initiator
    // that VEP calls a frameshift on a row that is only a 5' UTR variant (divergences 14 and 15)
    (
        "chr10_1093_AGCTGAGATGC/AGC",
        "ENST00000000001",
        "HGVSp",
        "ENSP00000000001.1:p.Met1?",
    ),
    (
        "chr10_1093_AGCTGAGATGC/AGC",
        "ENST00000000002",
        "HGVSp",
        "ENSP00000000002.1:p.Met1?",
    ),
    (
        "chr10_1095_CTGAGATG/CTG",
        "ENST00000000001",
        "HGVSp",
        "ENSP00000000001.1:p.Met1?",
    ),
    (
        "chr10_1095_CTGAGATG/CTG",
        "ENST00000000002",
        "HGVSp",
        "ENSP00000000002.1:p.Met1?",
    ),
    (
        "chr10_11295_GCACATCTGC/GCCGAAGC",
        "ENST00000000003",
        "HGVSp",
        "ENSP00000000003.1:p.Met1?",
    ),
    (
        "chr10_11299_ATCTG/ACCTGG",
        "ENST00000000003",
        "HGVSp",
        "ENSP00000000003.1:p.Met1?",
    ),
    (
        "chr10_11299_ATCTG/ATGTGACTTCTG",
        "ENST00000000003",
        "HGVSp",
        "-",
    ),
    // an insertion in the start codon that leaves the whole coding sequence behind it, with no ATG
    // at the old coding start (`c.2_3insCGAT` is `AT` `CGAT` `G`: the `ATG` and all that follows it
    // are four bases further along): VEP writes a frameshift and a start loss, and `p.Met1?`, for
    // a protein that has neither (divergence 14)
    (
        "chr10_11297_ACA/ACATCGA",
        "ENST00000000003",
        "Consequence",
        "start_retained_variant",
    ),
    ("chr10_11297_ACA/ACATCGA", "ENST00000000003", "HGVSp", "-"),
    (
        "chr10_11297_ACA/ACATCGA",
        "ENST00000000003",
        "IMPACT",
        "LOW",
    ),
    // an insertion at the edge of an exon that the 3'-rule's walk over the genome carries into the
    // intron: VEP writes no HGVSp, fastVEP walks the transcript instead (`vep_tab_intron_repeat`)
    (
        "chr10_10996_CAGCG/CAGCTCG",
        "ENST00000000003",
        "HGVSp",
        "ENSP00000000003.1:p.Pro101ArgfsTer4",
    ),
    (
        "chr10_10997_AGCG/AGCTCGGCCCG",
        "ENST00000000003",
        "HGVSp",
        "ENSP00000000003.1:p.Pro101ArgfsTer8",
    ),
];

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

    let (mut wrong, mut compared, mut known) = (Vec::new(), 0usize, 0usize);
    for (key, want) in &vep {
        let have = &fv[key];
        for (column, v) in want {
            compared += 1;
            if let Some((_, _, _, ours)) = KNOWN
                .iter()
                .find(|(rec, tr, c, _)| *rec == key.0 && *tr == key.1 && c == column)
            {
                assert_ne!(&have[column], v, "{} {column}: VEP now agrees", key.0);
                assert_eq!(&have[column], ours, "{} {column}", key.0);
                known += 1;
                continue;
            }
            if &have[column] != v {
                wrong.push(format!(
                    "{} {} {column}: VEP {v:?}, here {:?}",
                    key.0, key.1, have[column]
                ));
            }
        }
    }
    assert_eq!(
        known,
        KNOWN.len(),
        "a row `KNOWN` names is not in the golden"
    );
    assert!(compared > 90_000, "only {compared} cells compared");
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
