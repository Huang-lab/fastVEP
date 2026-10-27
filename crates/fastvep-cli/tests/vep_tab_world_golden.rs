//! The tab and VCF output over a generated world, against what Ensembl VEP 116.3
//! writes for the same world (issue #139).
//!
//! `fixtures/vep_tab_world` is 500 records on both strands of a small genome with
//! coding and non-coding genes: SNVs, MNVs, minimal indels and multi-allelic records.
//! `regen.sh` there runs `vep` from the 116.3 image in `--gff` mode over it three times
//! and keeps what VEP wrote: the default tab file, the tab file with every optional
//! column fastVEP has, and the VCF with the same flags. Nothing in those goldens was
//! worked out here.
//!
//! Every column (field, for the VCF) both write is compared for every row, so the tests
//! hold the things the reporter's pipeline reads to VEP's own output and not to a reading
//! of its source: the columns and their order, `Uploaded_variation`, `Location`, the `-`
//! of a value a row does not have, `Amino_acids`, `FLAGS`, `CANONICAL`, the HGVS group and
//! `HGVS_OFFSET`, `HGVSg`, `VARIANT_CLASS` and `REF_ALLELE`.
//!
//! The one place the comparison is not exact is `HGVSp`, and only for a synonymous change
//! across two codons, which VEP writes `p.ProLeu145=` and HGVS writes `p.Pro145_Leu146=`
//! (`docs/VEP_DIVERGENCE.md`, divergence 3). `documented_hgvsp_difference` is the whole of that
//! exception, and the tests fail if it starts to cover more than a few rows.

use std::collections::{BTreeMap, HashMap};
use std::io::Read;
use std::path::{Path, PathBuf};

use fastvep_cli::pipeline::{run_annotate, AnnotateConfig, ColumnFlags, PickFlags};
use flate2::read::GzDecoder;
use tempfile::TempDir;

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("fixtures/vep_tab_world")
        .join(name)
}

fn read_gz(name: &str) -> String {
    let mut text = String::new();
    GzDecoder::new(std::fs::File::open(fixture(name)).unwrap())
        .read_to_string(&mut text)
        .unwrap();
    text
}

/// Every optional column fastVEP writes, which is what `golden_flags.*` was made with.
const ALL_FLAGS: ColumnFlags = ColumnFlags {
    show_ref_allele: true,
    variant_class: true,
    symbol: true,
    biotype: true,
    canonical: true,
    numbers: true,
    hgvsg: true,
};

fn run(dir: &Path, format: &str, columns: ColumnFlags, hgvs: bool) -> String {
    let out = dir.join(format!("out.{format}"));
    // Copied, because the transcript cache a run writes goes beside its GFF3 and the
    // index of the FASTA beside it, and neither belongs in the source tree.
    let (gff3, fasta) = (dir.join("genes.gff3"), dir.join("genome.fa"));
    std::fs::copy(fixture("genes.gff3"), &gff3).unwrap();
    std::fs::copy(fixture("genome.fa"), &fasta).unwrap();
    run_annotate(AnnotateConfig {
        input: fixture("in.vcf").to_string_lossy().into(),
        output: out.to_string_lossy().into(),
        gff3: vec![gff3.to_string_lossy().into()],
        fasta: Some(fasta.to_string_lossy().into()),
        output_format: format.into(),
        pick: PickFlags::default(),
        hgvs,
        columns,
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
    std::fs::read_to_string(&out).unwrap()
}

/// A tab file as a header and its rows, each row keyed by what names it: the record, the
/// allele and the transcript.
struct Table {
    header: Vec<String>,
    rows: BTreeMap<(String, String, String), Vec<HashMap<String, String>>>,
}

fn table(text: &str) -> Table {
    let mut lines = text.lines().filter(|l| !l.starts_with("##"));
    let header: Vec<String> = lines
        .next()
        .expect("a header")
        .trim_start_matches('#')
        .split('\t')
        .map(String::from)
        .collect();
    let mut rows: BTreeMap<_, Vec<HashMap<String, String>>> = BTreeMap::new();
    for line in lines {
        let cells: Vec<&str> = line.split('\t').collect();
        assert_eq!(cells.len(), header.len(), "ragged row {line:?}");
        let row: HashMap<String, String> = header
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
    Table { header, rows }
}

fn strip_digits(s: &str) -> String {
    s.chars().filter(|c| !c.is_ascii_digit()).collect()
}

/// Where fastVEP's `HGVSp` differs from VEP's on purpose, in this world: a synonymous change
/// across two codons, which VEP writes `p.SerGly76=` and HGVS writes `p.Ser76_Gly77=`
/// (`docs/VEP_DIVERGENCE.md`, divergence 3). The `=` is `%3D` in the tab file.
///
/// It used to be wider: the frameshifts and extensions whose numbers differ, and the delins that
/// are placed differently. With VEP's own notation transcribed (`fastvep_hgvs::protein`) the
/// first is what Ensembl's `--gff` reader does to the UTR of a reverse-strand transcript and
/// the second is gone from a world this small, and `vep_tab_edges_golden.rs` lists the rows
/// of both one by one.
fn documented_hgvsp_difference(vep: &str, ours: &str) -> bool {
    vep.ends_with("%3D")
        && ours.ends_with("%3D")
        && ours.contains('_')
        && strip_digits(&ours.replace('_', "")) == strip_digits(vep)
}

/// Compare two tab tables column by column over every row. Returns the HGVSp cells that
/// differ for a documented reason.
fn compare_tabs(ours: &Table, theirs: &Table) -> Vec<String> {
    assert_eq!(ours.header, theirs.header, "the column header");
    assert_eq!(
        ours.rows.keys().collect::<Vec<_>>(),
        theirs.rows.keys().collect::<Vec<_>>(),
        "the rows"
    );
    let mut documented = Vec::new();
    let mut wrong = Vec::new();
    for (key, want_rows) in &theirs.rows {
        let have_rows = &ours.rows[key];
        assert_eq!(have_rows.len(), want_rows.len(), "{key:?}");
        for (have, want) in have_rows.iter().zip(want_rows) {
            for column in &ours.header {
                let (a, b) = (&have[column], &want[column]);
                if a == b {
                    continue;
                }
                if column == "HGVSp" && documented_hgvsp_difference(b, a) {
                    documented.push(format!("{key:?} HGVSp VEP {b} here {a}"));
                    continue;
                }
                wrong.push(format!("{key:?} {column}: VEP {b:?}, here {a:?}"));
            }
        }
    }
    assert!(
        wrong.is_empty(),
        "{} cells differ from VEP's, first ten:\n{}",
        wrong.len(),
        wrong
            .iter()
            .take(10)
            .cloned()
            .collect::<Vec<_>>()
            .join("\n")
    );
    documented
}

#[test]
fn the_default_tab_file_agrees_with_vep_cell_for_cell() {
    let dir = TempDir::new().unwrap();
    let ours = table(&run(dir.path(), "tab", ColumnFlags::default(), false));
    let theirs = table(&read_gz("golden_default.tab.gz"));
    assert_eq!(ours.header.len(), 17);
    assert!(
        compare_tabs(&ours, &theirs).is_empty(),
        "no HGVSp column here"
    );
    // the world is not trivial: a thousand rows, coding and non-coding
    assert!(theirs.rows.len() > 1000);
}

#[test]
fn the_tab_file_with_every_optional_column_agrees_with_vep_cell_for_cell() {
    let dir = TempDir::new().unwrap();
    let ours = table(&run(dir.path(), "tab", ALL_FLAGS, true));
    let theirs = table(&read_gz("golden_flags.tab.gz"));
    assert_eq!(ours.header.len(), 30);
    let documented = compare_tabs(&ours, &theirs);
    // Only a few cells may be excused. This is what stops the exception from growing
    // into a way of not comparing HGVSp.
    assert!(
        documented.len() <= 3,
        "{} HGVSp cells differ for a documented reason:\n{}",
        documented.len(),
        documented.join("\n")
    );
}

/// The VCF's `CSQ`, field by field, for the fields VEP and fastVEP both write.
#[test]
fn the_vcf_csq_agrees_with_vep_field_for_field() {
    let dir = TempDir::new().unwrap();
    let ours = run(dir.path(), "vcf", ALL_FLAGS, true);
    let theirs = read_gz("golden_flags.vcf.gz");

    let csq_fields_of_ours: Vec<&str> = ours
        .lines()
        .find(|l| l.starts_with("##INFO=<ID=CSQ"))
        .and_then(|l| l.split("Format: ").nth(1))
        .and_then(|rest| rest.split("\">").next())
        .expect("a CSQ header")
        .split('|')
        .collect();
    let csq_fields_of_vep: Vec<&str> = theirs
        .lines()
        .find_map(|l| l.strip_prefix("##CSQ-fields: "))
        .expect("the golden declares its fields")
        .split('|')
        .collect();
    let shared: Vec<&str> = csq_fields_of_vep
        .iter()
        .copied()
        .filter(|f| csq_fields_of_ours.contains(f))
        .collect();
    // everything fastVEP is meant to share with VEP is shared
    for field in [
        "Allele",
        "Consequence",
        "IMPACT",
        "SYMBOL",
        "Gene",
        "Feature_type",
        "Feature",
        "BIOTYPE",
        "EXON",
        "INTRON",
        "HGVSc",
        "HGVSp",
        "cDNA_position",
        "CDS_position",
        "Protein_position",
        "Amino_acids",
        "Codons",
        "Existing_variation",
        "DISTANCE",
        "STRAND",
        "FLAGS",
        "SYMBOL_SOURCE",
        "HGNC_ID",
        "CANONICAL",
        "REF_ALLELE",
        "VARIANT_CLASS",
        "HGVS_OFFSET",
        "HGVSg",
    ] {
        assert!(shared.contains(&field), "{field} is not in both CSQs");
    }

    // (record, occurrence, allele, transcript) -> the entry as field name to value
    type Entries = BTreeMap<(String, usize, String, String), HashMap<String, String>>;
    let entries = |lines: Vec<(String, String)>, fields: &[&str]| -> Entries {
        let mut out = Entries::new();
        let mut seen: HashMap<String, usize> = HashMap::new();
        for (record, csq) in lines {
            let n = seen.entry(record.clone()).or_default();
            *n += 1;
            for entry in csq.split(',').filter(|e| !e.is_empty()) {
                let mut values: Vec<&str> = entry.split('|').collect();
                values.resize(fields.len(), "");
                let map: HashMap<String, String> = fields
                    .iter()
                    .map(|f| f.to_string())
                    .zip(values.iter().map(|v| v.to_string()))
                    .collect();
                out.insert(
                    (
                        record.clone(),
                        *n,
                        map["Allele"].clone(),
                        map["Feature"].clone(),
                    ),
                    map,
                );
            }
        }
        out
    };
    let ours_lines: Vec<(String, String)> = ours
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| {
            let c: Vec<&str> = l.split('\t').collect();
            let csq = c[7]
                .split(';')
                .find_map(|kv| kv.strip_prefix("CSQ="))
                .unwrap_or("")
                .to_string();
            (format!("{} {} {} {}", c[0], c[1], c[3], c[4]), csq)
        })
        .collect();
    let theirs_lines: Vec<(String, String)> = theirs
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| {
            let c: Vec<&str> = l.split('\t').collect();
            (
                format!("{} {} {} {}", c[0], c[1], c[3], c[4]),
                c[5].to_string(),
            )
        })
        .collect();
    let ours = entries(ours_lines, &csq_fields_of_ours);
    let theirs = entries(theirs_lines, &csq_fields_of_vep);

    assert_eq!(
        ours.keys().collect::<Vec<_>>(),
        theirs.keys().collect::<Vec<_>>(),
        "the CSQ entries"
    );
    let mut wrong = Vec::new();
    let mut documented = 0;
    for (key, want) in &theirs {
        let have = &ours[key];
        for field in &shared {
            let (a, b) = (&have[*field], &want[*field]);
            if a == b {
                continue;
            }
            if *field == "HGVSp"
                && documented_hgvsp_difference(if b.is_empty() { "-" } else { b }, a)
            {
                documented += 1;
                continue;
            }
            wrong.push(format!("{key:?} {field}: VEP {b:?}, here {a:?}"));
        }
    }
    assert!(
        wrong.is_empty(),
        "{} fields differ from VEP's, first ten:\n{}",
        wrong.len(),
        wrong
            .iter()
            .take(10)
            .cloned()
            .collect::<Vec<_>>()
            .join("\n")
    );
    assert!(documented <= 3, "{documented} HGVSp fields excused");
}
