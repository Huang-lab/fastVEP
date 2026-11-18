//! The tab output for structural variants, against what Ensembl VEP 116.3 writes for the same
//! records (issue #185).
//!
//! `fixtures/vep_tab_sv` is 6,692 records on a 14 kb contig of ten transcripts: every class VEP
//! has a rule for (`<DEL>`, `<DUP>`, `<INV>`, `<CNV>`, the tandem and mobile element classes, the copy
//! number alleles, `<INS>`, breakends in every notation), each laid over every place on every
//! transcript where one end of it can be on another kind of sequence than the other, and the records
//! that are about how a record is read: an `END` that `SVLEN` overrules, a type taken from `SVTYPE`,
//! several ALTs, a tandem repeat the INFO of which VEP turns into bases, and the records VEP skips.
//! `regen.sh` there runs `vep` from the 116.3 image in `--gff` mode over them and keeps what it wrote.
//! Nothing in the golden was worked out here.
//!
//! Twenty-three columns are compared cell by cell (`COMPARED`), and the rows are named by the record,
//! the allele and the feature, and come in VEP's order. `Existing_variation`, `FLAGS`,
//! `SYMBOL_SOURCE` and `HGNC_ID` are not compared: a `--gff` run of VEP has no variation source, no
//! `cds_start_NF` or `cds_end_NF`, and no symbol source or HGNC identifier of a GFF3 to write. The
//! cells where fastVEP is not VEP are listed in `KNOWN`, each with what fastVEP writes and why, and
//! each is asserted to still differ. `docs/STRUCTURAL_VARIANTS.md` has the rules, and what is left
//! out of the fixture because VEP's answer for it cannot be repeated.

use std::collections::BTreeMap;
use std::io::Read;
use std::path::{Path, PathBuf};

use fastvep_cli::pipeline::{run_annotate, AnnotateConfig, ColumnFlags, CompatFlags, PickFlags};
use flate2::read::GzDecoder;
use tempfile::TempDir;

/// The columns compared, cell by cell, on every row.
const COMPARED: &[&str] = &[
    "Location",
    "Consequence",
    "IMPACT",
    "cDNA_position",
    "CDS_position",
    "Protein_position",
    "Amino_acids",
    "Codons",
    "DISTANCE",
    "STRAND",
    "Gene",
    "Feature_type",
    "EXON",
    "INTRON",
    "VARIANT_CLASS",
    "REF_ALLELE",
    "HGVSc",
    "HGVSp",
    "HGVSg",
    "HGVS_OFFSET",
    "SYMBOL",
    "BIOTYPE",
    "CANONICAL",
];

/// Why a cell is not VEP's.
#[derive(Debug, Clone, Copy)]
#[allow(dead_code)] // the reason is documentation: it is read by whoever finds the row
enum Why {
    /// A tandem repeat with more than one allele is read as a small variant with several ALTs, and
    /// VEP 116 describes the alleles of such a record as the record, not as each would be alone:
    /// an insertion into the repeat beside a deletion allele is not moved in the 3' direction (it is
    /// `ins` where it was put, and fastVEP's is the `dup` of the 3' rule), and the `HGVSg` of a
    /// `delins` starts six bases after where the change does. `docs/VEP_DIVERGENCE.md`, Part 1,
    /// section 11, is the same difference in a plain VCF record.
    SeveralAlleles,
}

/// Cells where fastVEP and stock VEP 116.3 differ: `(record, allele, transcript, column, what
/// fastVEP writes, why)`.
const KNOWN: &[(&str, &str, &str, &str, &str, Why)] = &[
    (
        "tr_two_alleles",
        "CACACACACACACACACACACACACACACACACACACACACACA",
        "ENST00000000403",
        "HGVSc",
        "ENST00000000403.1:c.283+97_283+100dup",
        Why::SeveralAlleles,
    ),
    (
        "tr_two_alleles",
        "CACACACACACACACACACACACACACACACACACACACACACA",
        "ENST00000000403",
        "HGVS_OFFSET",
        "-40",
        Why::SeveralAlleles,
    ),
    (
        "tr_two_alleles",
        "CACACACACACACACACACACACACACACACACACA",
        "ENST00000000403",
        "HGVS_OFFSET",
        "-36",
        Why::SeveralAlleles,
    ),
    (
        "tr_two_alleles_no_rn",
        "CACACACACACACACACA",
        "ENST00000000403",
        "HGVS_OFFSET",
        "-18",
        Why::SeveralAlleles,
    ),
    (
        "tr_two_alleles_two_units",
        "CACACAGTGT",
        "ENST00000000403",
        "HGVSg",
        "chr1:g.3406_3439delinsGTGT",
        Why::SeveralAlleles,
    ),
];

fn here(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("fixtures/vep_tab_sv")
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

fn golden_text() -> String {
    let mut text = String::new();
    GzDecoder::new(std::fs::File::open(here("golden_flags.tab.gz")).unwrap())
        .read_to_string(&mut text)
        .unwrap();
    text
}

fn golden() -> Table {
    table(&golden_text())
}

/// The key of every row of a table, in the order the file has them.
fn order(text: &str) -> Vec<(String, String, String)> {
    let mut lines = text.lines().filter(|l| !l.starts_with("##"));
    let header: Vec<&str> = lines
        .next()
        .expect("a header")
        .trim_start_matches('#')
        .split('\t')
        .collect();
    let at = |name: &str| header.iter().position(|h| *h == name).expect(name);
    let (record, allele, feature) = (at("Uploaded_variation"), at("Allele"), at("Feature"));
    lines
        .map(|line| {
            let cells: Vec<&str> = line.split('\t').collect();
            (
                cells[record].to_string(),
                cells[allele].to_string(),
                cells[feature].to_string(),
            )
        })
        .collect()
}

fn ours() -> Table {
    table(&annotate(&here("in.vcf")))
}

/// What `fastvep annotate` writes for a VCF, over the fixture's genes and genome.
fn annotate(input: &Path) -> String {
    let dir = TempDir::new().unwrap();
    // Copied, because the transcript cache a run writes goes beside its GFF3 and the index of the
    // FASTA beside it, and neither belongs in the source tree.
    let (gff3, fasta) = (dir.path().join("genes.gff3"), dir.path().join("genome.fa"));
    std::fs::copy(here("genes.gff3"), &gff3).unwrap();
    std::fs::copy(here("genome.fa"), &fasta).unwrap();
    let out = dir.path().join("out.tab");
    run_annotate(AnnotateConfig {
        input: input.to_string_lossy().into(),
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
        compat: CompatFlags {
            total_length: true,
            ..Default::default()
        },
        exclude: Default::default(),
    })
    .expect("annotation should succeed");
    std::fs::read_to_string(&out).unwrap()
}

/// Every row VEP wrote, fastVEP wrote, and no other: a record VEP skips has none, and so does one
/// whose alleles are all the reference.
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

/// The rows come in the order VEP writes them: the records as the file has them, and for each one its
/// rows in VEP's order (the transcripts by their IDs, the alleles of a breakend as it lists them).
#[test]
fn the_rows_come_in_the_order_vep_writes_them() {
    let theirs = order(&golden_text());
    let ours = order(&annotate(&here("in.vcf")));
    assert_eq!(theirs.len(), ours.len(), "the number of rows");
    let first = theirs.iter().zip(&ours).position(|(vep, fv)| vep != fv);
    assert!(
        first.is_none(),
        "row {:?} is not where VEP puts it: VEP {:?}, fastVEP {:?}",
        first,
        first.map(|at| &theirs[at]),
        first.map(|at| &ours[at])
    );
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
                let known = KNOWN.iter().find(|(rec, allele, tr, c, _, _)| {
                    *rec == key.0 && *allele == key.1 && *tr == key.2 && c == col
                });
                if let Some((rec, allele, tr, c, ours, _)) = known {
                    assert_ne!(v[*col], f[*col], "{rec} {tr} {c}: VEP now agrees");
                    assert_eq!(&f[*col], ours, "{rec} {tr} {c}");
                    seen.insert((*rec, *allele, *tr, *c));
                } else if v[*col] != f[*col] {
                    wrong.push(format!(
                        "{} {} {} {col}: VEP {:?}, fastVEP {:?}",
                        key.0, key.1, key.2, v[*col], f[*col]
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
    assert!(compared > 300_000, "only {compared} cells compared");
    assert!(
        wrong.is_empty(),
        "{} of {compared} cells differ from VEP, e.g.\n{}",
        wrong.len(),
        wrong[..wrong.len().min(12)].join("\n")
    );
}

/// The fixture asks what it is for: every term a structural variant can have, every class of one,
/// the flank rows with their distances, the rows of a breakend's mate, and a tandem repeat that
/// became a small variant.
#[test]
fn the_fixture_covers_the_shapes_the_rules_are_about() {
    let vep = golden();
    let rows = || vep.values().flatten();
    let terms: std::collections::BTreeSet<&str> =
        rows().flat_map(|r| r["Consequence"].split(',')).collect();
    for term in [
        "transcript_ablation",
        "transcript_amplification",
        "feature_elongation",
        "feature_truncation",
        "frameshift_variant",
        "inframe_deletion",
        "stop_lost",
        "start_lost",
        "start_retained_variant",
        "coding_sequence_variant",
        "coding_transcript_variant",
        "5_prime_UTR_variant",
        "3_prime_UTR_variant",
        "intron_variant",
        "splice_polypyrimidine_tract_variant",
        "non_coding_transcript_exon_variant",
        "non_coding_transcript_variant",
        "NMD_transcript_variant",
        "upstream_gene_variant",
        "downstream_gene_variant",
        "intergenic_variant",
    ] {
        assert!(terms.contains(term), "no row has {term}");
    }
    let classes: std::collections::BTreeSet<&str> =
        rows().map(|r| r["VARIANT_CLASS"].as_str()).collect();
    for class in [
        "insertion",
        "deletion",
        "duplication",
        "tandem_duplication",
        "tandem_repeat",
        "copy_number_variation",
        "inversion",
        "chromosome_breakpoint",
        "Alu_insertion",
        "Alu_deletion",
        "LINE1_insertion",
        "LINE1_deletion",
        "SVA_insertion",
        "SVA_deletion",
        "HERV_insertion",
        "HERV_deletion",
        "mobile_element_insertion",
        "mobile_element_deletion",
    ] {
        assert!(classes.contains(class), "no row of the class {class}");
    }
    assert!(
        rows().any(|r| r["Allele"].contains('[') || r["Allele"].contains(']')),
        "no breakend allele"
    );
    assert!(
        rows().any(|r| r["Allele"].ends_with('.') && r["Allele"].len() > 1),
        "no breakend of the record's own"
    );
    assert!(
        rows().any(|r| r["DISTANCE"] != "-" && r["VARIANT_CLASS"] == "deletion"),
        "no flank row with a distance"
    );
    assert!(
        rows().any(|r| r["Feature"] == "-" && r["Consequence"] == "intergenic_variant"),
        "no intergenic row"
    );
    assert!(
        rows().any(|r| r["VARIANT_CLASS"] == "indel" && r["Uploaded_variation"].starts_with("tr_")),
        "no tandem repeat made a small variant"
    );
}

/// A record that VEP reads as a structural variant has no `REF_ALLELE`, no `Existing_variation`
/// and no HGVS in any row, and its `Allele` is the class or the breakend, never a base.
#[test]
fn a_structural_variant_row_has_no_bases_of_its_own() {
    let fv = ours();
    let mut checked = 0;
    for ((name, allele, _), rows) in &fv {
        // the records the world has that are small variants, or become one
        if ["tr_", "plain_", "bases_"]
            .iter()
            .any(|p| name.starts_with(p))
        {
            continue;
        }
        for r in rows {
            checked += 1;
            assert_eq!(r["REF_ALLELE"], "-", "{name} {allele}");
            for col in ["HGVSc", "HGVSp", "HGVSg"] {
                assert_eq!(r[col], "-", "{name} {allele} {col}");
            }
            assert!(
                !allele.chars().all(|c| "ACGTN".contains(c)),
                "{name}: the allele {allele} is bases"
            );
        }
    }
    assert!(checked > 10_000, "only {checked} rows checked");
}

/// A record whose only ALTs say nothing (`.`, `<NON_REF>`) has no row, over a transcript and away
/// from every one. VEP reads a `.` beside an `SVTYPE` as an allele of that class and writes rows for
/// it, which fastVEP does not, as it does not for the same ALT of a record without one (a `*` beside
/// an `SVTYPE` is an allele of the class in both).
#[test]
fn a_record_whose_only_alleles_are_placeholders_has_no_row() {
    let dir = TempDir::new().unwrap();
    let vcf = dir.path().join("placeholders.vcf");
    std::fs::write(
        &vcf,
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
         chr1\t1060\tdot_over_a_gene\tN\t.\t.\t.\tSVTYPE=DEL;END=1085\n\
         chr1\t13000\tdot_away\tN\t.\t.\t.\tSVTYPE=DEL;END=13050\n",
    )
    .unwrap();
    let written = annotate(&vcf);
    let rows: Vec<&str> = written.lines().filter(|l| !l.starts_with('#')).collect();
    assert!(rows.is_empty(), "rows for placeholder alleles: {rows:?}");
}
