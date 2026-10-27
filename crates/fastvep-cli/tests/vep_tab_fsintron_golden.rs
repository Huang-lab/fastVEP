//! The tab output for a change in or near a transcript with a very short intron, against what
//! Ensembl VEP 116.3 writes for the same change (issue #139).
//!
//! Ensembl calls an intron of 13 bases or fewer a frameshift intron (`end - start <= 12`), and a
//! transcript that has one is read differently everywhere near its exons, not only in the intron:
//! a change in the intron earns no splice term and no `intron_variant` and, inside the coding
//! region, is a `coding_sequence_variant`; every exon of the transcript is 12 bases wider for the
//! "does it touch an exon" test, which takes `splice_polypyrimidine_tract_variant` from the 3 to
//! 12 bases before each of them and makes a change in the intron `coding` for `HGVSp`.
//!
//! `fixtures/vep_tab_fsintron` is a world of seven transcripts with introns of 1 to 60 bases, on
//! both strands, coding and not, two of them controls whose introns are all 14 bases or more,
//! and 2,879 single-ALT records around the exon edges and in the short introns: every substitution,
//! and every eighth of the other shapes. `regen.sh` there
//! runs `vep` from the 116.3 image in `--gff` mode over them and keeps what it wrote. Nothing in
//! the golden was worked out here. The same rules were checked against VEP's own cache on the
//! five transcripts of chr21 and chr22 that have such an intron (SHANK3 among them), and the
//! numbers are in `docs/VEP_DIVERGENCE.md`.
//!
//! Every column is compared, and the cells where fastVEP is not VEP are listed one by one, as
//! `vep_tab_edges_golden.rs` does. There are two, and they are one thing: a deletion of the last
//! base of an exon, in a repeat that the intron goes on with, which the 3'-rule's walk over the
//! genome carries out of the exon. VEP writes no HGVSp for it. fastVEP walks the sequence of the
//! transcript instead and writes the notation of the change the VCF makes in it (the next
//! fixture, `vep_tab_intron_repeat`, shows that this is what the same change is where the
//! transcript has no intron).

use std::collections::BTreeMap;
use std::io::Read;
use std::path::{Path, PathBuf};

use fastvep_cli::pipeline::{run_annotate, AnnotateConfig, ColumnFlags, PickFlags};
use flate2::read::GzDecoder;
use tempfile::TempDir;

/// The columns compared, cell by cell, on every row.
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
    "HGVSc",
    "HGVSp",
    "HGVS_OFFSET",
];

/// Cells where fastVEP and stock VEP 116.3 differ and VEP is the one that is wrong:
/// `(record, transcript, column, what fastVEP writes)`.
const KNOWN: &[(&str, &str, &str, &str)] = &[
    (
        "chr1_16094_ACTAGA/A",
        "ENST00000000106",
        "HGVSp",
        "ENSP00000000106.1:p.Thr19SerfsTer8",
    ),
    (
        "chr1_4276_GTGGGG/G",
        "ENST00000000102",
        "HGVSp",
        "ENSP00000000102.1:p.Pro56ArgfsTer22",
    ),
];

/// The transcripts whose introns include one of 13 bases or fewer, and the two that do not.
const FRAMESHIFT: &[&str] = &[
    "ENST00000000101",
    "ENST00000000102",
    "ENST00000000103",
    "ENST00000000104",
    "ENST00000000105",
];
const CONTROLS: &[&str] = &["ENST00000000106", "ENST00000000107"];

fn here(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("fixtures/vep_tab_fsintron")
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
    let mut seen = std::collections::BTreeSet::new();
    for (key, vep_rows) in &vep {
        let fv_rows = &fv[key];
        assert_eq!(vep_rows.len(), fv_rows.len(), "{key:?}");
        for (v, f) in vep_rows.iter().zip(fv_rows) {
            for col in COMPARED {
                compared += 1;
                let known = KNOWN
                    .iter()
                    .find(|(rec, tr, c, _)| *rec == key.0 && *tr == key.2 && c == col);
                if let Some((rec, tr, c, ours)) = known {
                    assert_ne!(v[*col], f[*col], "{rec} {tr} {c}: VEP now agrees");
                    assert_eq!(&f[*col], ours, "{rec} {tr} {c}");
                    seen.insert((*rec, *tr, *c));
                } else if v[*col] != f[*col] {
                    wrong.push(format!(
                        "{} {} {col}: VEP {:?}, fastVEP {:?}",
                        key.0, key.2, v[*col], f[*col]
                    ));
                }
            }
        }
    }
    assert_eq!(
        seen.len(),
        KNOWN.len(),
        "a row `KNOWN` names is not in the golden"
    );
    assert!(compared > 60_000, "only {compared} cells compared");
    assert!(
        wrong.is_empty(),
        "{} of {compared} cells differ from VEP, e.g.\n{}",
        wrong.len(),
        wrong[..wrong.len().min(12)].join("\n")
    );
}

/// How far into the intron, from the exon it precedes, an `HGVSc` of one position puts it
/// (`c.123-14A>G` is 14 before the exon, `c.123+4A>G` after the one before), `None` for a
/// description with no offset or a range.
fn intron_offset(hgvsc: &str) -> Option<i64> {
    let position = hgvsc.rsplit_once(':')?.1.get(2..)?;
    let position = position.trim_start_matches(['-', '*']);
    let digits = position.chars().take_while(char::is_ascii_digit).count();
    let rest = &position[digits..];
    let sign = rest.chars().next().filter(|c| matches!(c, '+' | '-'))?;
    let n: i64 = rest[1..]
        .chars()
        .take_while(char::is_ascii_digit)
        .collect::<String>()
        .parse()
        .ok()?;
    Some(if sign == '-' { -n } else { n })
}

/// The fixture exercises what it is for. This reads VEP's golden only, so it says what the
/// rule is, and the test above says fastVEP follows it: the polypyrimidine tract is 3 to 17
/// bases before an exon, and in a transcript with a frameshift intron an SNV in it is called
/// only 13 or more bases before one, where in the controls it is called all along it. A change
/// in a short intron inside the coding sequence is a `coding_sequence_variant` and nothing else.
#[test]
fn the_fixture_covers_the_rules_it_is_about() {
    let vep = golden();
    let tract = |transcripts: &[&str]| -> Vec<i64> {
        vep.iter()
            .filter(|((_, _, feature), _)| transcripts.contains(&feature.as_str()))
            .flat_map(|(_, rows)| rows)
            .filter(|r| r["VARIANT_CLASS"] == "SNV")
            .filter(|r| r["Consequence"].contains("splice_polypyrimidine_tract_variant"))
            .filter_map(|r| intron_offset(&r["HGVSc"]))
            .filter(|o| *o < 0)
            .collect()
    };
    let frameshift = tract(FRAMESHIFT);
    let controls = tract(CONTROLS);
    assert!(
        frameshift.len() > 8 && frameshift.iter().all(|o| (-17..=-13).contains(o)),
        "frameshift-intron transcripts: {frameshift:?}"
    );
    assert!(
        controls.len() > 30 && controls.iter().any(|o| (-12..=-3).contains(o)),
        "controls: {controls:?}"
    );

    let only_coding_sequence = vep
        .iter()
        .filter(|((_, _, feature), _)| FRAMESHIFT[..3].contains(&feature.as_str()))
        .flat_map(|(_, rows)| rows)
        .filter(|r| r["Consequence"] == "coding_sequence_variant")
        .count();
    assert!(
        only_coding_sequence > 50,
        "{only_coding_sequence} rows that are only `coding_sequence_variant`"
    );
    // both strands, and a transcript with no protein
    let on = |strand: &str| vep.values().flatten().any(|r| r["STRAND"] == strand);
    assert!(on("1") && on("-1"));
    assert!(vep.values().flatten().any(|r| r["Consequence"]
        .contains("non_coding_transcript_variant")
        && FRAMESHIFT.contains(&r["Feature"].as_str())));
    // a protein change for a change in a one-base intron that the 3'-rule carried into the exon
    assert!(
        vep.values()
            .flatten()
            .any(|r| r["Consequence"] == "coding_sequence_variant" && r["HGVSp"] != "-"),
        "no HGVSp for a change in a short intron"
    );
}
