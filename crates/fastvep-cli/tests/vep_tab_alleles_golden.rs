//! Alleles that are not plain bases, against what Ensembl VEP 116.3 writes for them (issue #139).
//!
//! `fixtures/vep_tab_alleles` is 171 records at nine kinds of place (a coding base on each
//! strand, a UTR base, an intronic one, a non-coding exon, an intergenic one) whose ALT or REF is
//! `*`, `.`, or holds an ambiguity code (`N`, `R`, `Y`), alone and beside a real ALT, in SNVs,
//! MNVs, insertions and replacements. `regen.sh` there runs `vep` from the 116.3 image in `--gff`
//! mode over them and keeps what it wrote, as a tab file and as the `CSQ` of a VCF.
//!
//! What VEP does, and what this holds fastVEP to:
//!
//! - `*` is a spanning deletion and carries no consequence of its own: no row for it, a record
//!   whose only ALT is `*` has no row at all (and, in a VCF, no `CSQ`), and the real ALTs beside
//!   it are annotated as they would be alone.
//! - `.` says no ALT was called: the record is not annotated.
//! - An allele that is not plain bases cannot be translated, so it has no `Amino_acids`, no
//!   `HGVSp` and no `HGVSc`, and a change that alters the reading frame is a `frameshift_variant`
//!   and anything else in the coding sequence is `coding_sequence_variant`, the same whether the
//!   unknown base is in the ALT or in the REF. fastVEP called an `N` a `missense_variant` with
//!   `p.Lys31Xaa`, and wrote a `frameshift_variant` for a `<*>` that sat beside a real ALT.
//! - Ambiguity codes are complemented as Ensembl does on the reverse strand (`R` is `Y`).
//!
//! `HGVSc`, `HGVS_OFFSET` and `HGVSg` of an allele in a record with several ALTs are not compared:
//! Ensembl 116 describes them differently from the allele alone (divergence 11 in
//! `docs/VEP_DIVERGENCE.md`), `*` makes a record multi-allelic, and they are the same
//! difference with or without it.
//!
//! `non_ref.vcf` holds the two placeholders VEP does not handle sensibly, `<NON_REF>` (it skips the
//! whole record) and `<*>` (it writes a `coding_sequence_variant` for it); what fastVEP does
//! instead is asserted on its own below.

use std::collections::{BTreeMap, HashMap};
use std::io::Read;
use std::path::{Path, PathBuf};

use fastvep_cli::pipeline::{run_annotate, AnnotateConfig, ColumnFlags, PickFlags};
use flate2::read::GzDecoder;
use tempfile::TempDir;

fn here(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("fixtures/vep_tab_alleles")
        .join(name)
}

fn world(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("fixtures/vep_tab_world")
        .join(name)
}

fn read_gz(name: &str) -> String {
    let mut text = String::new();
    GzDecoder::new(std::fs::File::open(here(name)).unwrap())
        .read_to_string(&mut text)
        .unwrap();
    text
}

fn run(input: &Path, format: &str) -> String {
    let dir = TempDir::new().unwrap();
    // Copied, because the transcript cache a run writes goes beside its GFF3 and the index of
    // the FASTA beside it, and neither belongs in the source tree.
    let (gff3, fasta) = (dir.path().join("genes.gff3"), dir.path().join("genome.fa"));
    std::fs::copy(world("genes.gff3"), &gff3).unwrap();
    std::fs::copy(world("genome.fa"), &fasta).unwrap();
    let out = dir.path().join(format!("out.{format}"));
    run_annotate(AnnotateConfig {
        input: input.to_string_lossy().into(),
        output: out.to_string_lossy().into(),
        gff3: vec![gff3.to_string_lossy().into()],
        fasta: Some(fasta.to_string_lossy().into()),
        output_format: format.into(),
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
        compat: Default::default(),
        exclude: Default::default(),
    })
    .expect("annotation should succeed");
    std::fs::read_to_string(&out).unwrap()
}

type Row = BTreeMap<String, String>;
type Table = BTreeMap<(String, String, String), Row>;

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
        let key = (
            row["Uploaded_variation"].clone(),
            row["Allele"].clone(),
            row["Feature"].clone(),
        );
        assert!(
            rows.insert(key.clone(), row).is_none(),
            "two rows for {key:?}"
        );
    }
    rows
}

/// Columns that describe an allele alone and so differ from VEP's for an allele in a record with
/// several ALTs, for the reason in the module comment.
const MULTI_ALT_DIVERGENT: &[&str] = &["HGVSc", "HGVS_OFFSET", "HGVSg"];

#[test]
fn the_tab_file_is_vep_s_for_every_row() {
    let ours = table(&run(&here("in.vcf"), "tab"));
    let theirs = table(&read_gz("golden_flags.tab.gz"));

    let extra: Vec<_> = ours
        .keys()
        .filter(|k| !theirs.contains_key(*k))
        .take(5)
        .collect();
    let missing: Vec<_> = theirs
        .keys()
        .filter(|k| !ours.contains_key(*k))
        .take(5)
        .collect();
    assert!(extra.is_empty(), "rows only fastVEP writes: {extra:?}");
    assert!(missing.is_empty(), "rows only VEP writes: {missing:?}");

    let mut wrong = Vec::new();
    for (key, want) in &theirs {
        let have = &ours[key];
        // `Uploaded_variation` of a record with several ALTs ends in the allele string
        let multi = key.0.matches('/').count() > 1;
        for (column, v) in want {
            if multi && MULTI_ALT_DIVERGENT.contains(&column.as_str()) {
                continue;
            }
            if &have[column] != v {
                wrong.push(format!(
                    "{} {} {column}: VEP {v:?}, here {:?}",
                    key.0, key.2, have[column]
                ));
            }
        }
    }
    assert!(
        wrong.is_empty(),
        "{} cells differ from VEP's, first ten:\n{}",
        wrong.len(),
        wrong[..wrong.len().min(10)].join("\n")
    );
}

/// The VCF: a record with a `*` and no other ALT has no CSQ, `.` is passed through unannotated,
/// and the entries of every other record are VEP's, allele by allele.
#[test]
fn the_vcf_csq_is_vep_s_entry_for_entry() {
    let ours = run(&here("in.vcf"), "vcf");
    let theirs = read_gz("golden_flags.vcf.gz");

    let fields_of = |text: &str, key: &str| -> Vec<String> {
        text.lines()
            .find_map(|l| l.strip_prefix(key))
            .map(|rest| rest.split('|').map(String::from).collect())
            .unwrap_or_default()
    };
    let ours_fields: Vec<String> = ours
        .lines()
        .find(|l| l.starts_with("##INFO=<ID=CSQ"))
        .and_then(|l| l.split("Format: ").nth(1))
        .and_then(|rest| rest.split("\">").next())
        .expect("a CSQ header")
        .split('|')
        .map(String::from)
        .collect();
    let vep_fields = fields_of(&theirs, "##CSQ-fields: ");
    let shared: Vec<&String> = vep_fields
        .iter()
        .filter(|f| ours_fields.contains(f))
        .collect();

    // record -> allele+transcript -> field -> value
    type Entries = BTreeMap<(String, String, String), HashMap<String, String>>;
    let entries = |lines: Vec<(String, String)>, fields: &[String]| -> Entries {
        let mut out = Entries::new();
        for (record, csq) in lines {
            for entry in csq.split(',').filter(|e| !e.is_empty()) {
                let mut values: Vec<&str> = entry.split('|').collect();
                values.resize(fields.len(), "");
                let map: HashMap<String, String> = fields
                    .iter()
                    .cloned()
                    .zip(values.iter().map(|v| v.to_string()))
                    .collect();
                out.insert(
                    (
                        record.clone(),
                        map["Allele"].clone(),
                        map["Feature"].clone(),
                    ),
                    map,
                );
            }
        }
        out
    };
    let lines = |text: &str, csq_column: usize, info: bool| -> Vec<(String, String)> {
        text.lines()
            .filter(|l| !l.starts_with('#'))
            .map(|l| {
                let c: Vec<&str> = l.split('\t').collect();
                let csq = if info {
                    c[7].split(';')
                        .find_map(|kv| kv.strip_prefix("CSQ="))
                        .unwrap_or("")
                        .to_string()
                } else {
                    c[csq_column].to_string()
                };
                (format!("{} {} {} {}", c[0], c[1], c[3], c[4]), csq)
            })
            .collect()
    };
    let ours_lines = lines(&ours, 7, true);
    let theirs_lines = lines(&theirs, 5, false);

    // `.` is not annotated by either, and VEP drops the record where fastVEP keeps the line
    let records_with_csq = |l: &[(String, String)]| -> Vec<String> {
        l.iter()
            .filter(|(_, c)| !c.is_empty())
            .map(|(r, _)| r.clone())
            .collect()
    };
    assert_eq!(
        records_with_csq(&ours_lines),
        records_with_csq(&theirs_lines),
        "records with a CSQ"
    );
    assert!(
        ours_lines
            .iter()
            .any(|(r, c)| r.ends_with(" *") && c.is_empty()),
        "a record whose only ALT is `*` keeps its line and has no CSQ"
    );

    let (ours, theirs) = (
        entries(ours_lines, &ours_fields),
        entries(theirs_lines, &vep_fields),
    );
    assert_eq!(
        ours.keys().collect::<Vec<_>>(),
        theirs.keys().collect::<Vec<_>>(),
        "the CSQ entries"
    );
    let mut wrong = Vec::new();
    for (key, want) in &theirs {
        let multi = key
            .0
            .rsplit(' ')
            .next()
            .is_some_and(|alt| alt.contains(','));
        for field in &shared {
            if multi && MULTI_ALT_DIVERGENT.contains(&field.as_str()) {
                continue;
            }
            if ours[key][*field] != want[*field] {
                wrong.push(format!(
                    "{key:?} {field}: VEP {:?}, here {:?}",
                    want[*field], ours[key][*field]
                ));
            }
        }
    }
    assert!(
        wrong.is_empty(),
        "{} fields differ from VEP's, first ten:\n{}",
        wrong.len(),
        wrong[..wrong.len().min(10)].join("\n")
    );
}

/// A record whose only ALT is `*` is in the JSON, with no transcript consequences, as in VEP's.
#[test]
fn json_keeps_a_star_only_record_and_gives_it_no_consequences() {
    let json = run(&here("in.vcf"), "json");
    let records: Vec<serde_json::Value> =
        serde_json::from_str(&json).expect("the JSON output is an array of records");
    let star_only: Vec<_> = records
        .iter()
        .filter(|r| {
            r["allele_string"]
                .as_str()
                .is_some_and(|a| a.ends_with("/*") && a.matches('/').count() == 1)
        })
        .collect();
    assert!(
        !star_only.is_empty(),
        "no `*`-only record in the JSON output"
    );
    for r in star_only {
        assert!(
            r["transcript_consequences"]
                .as_array()
                .is_none_or(|a| a.is_empty()),
            "{r}"
        );
    }
    // and a `*` beside a real ALT adds no consequence of its own
    for r in records.iter().filter(|r| {
        r["allele_string"]
            .as_str()
            .is_some_and(|a| a.ends_with("/*") && a.matches('/').count() > 1)
    }) {
        for tc in r["transcript_consequences"]
            .as_array()
            .into_iter()
            .flatten()
        {
            assert_ne!(tc["variant_allele"], "*", "{r}");
        }
    }
}

/// What fastVEP does with the placeholders VEP does not handle sensibly: no row for the
/// placeholder itself, and the real ALT beside it annotated as it would be alone.
#[test]
fn a_non_ref_or_star_placeholder_has_no_row_and_the_real_alt_beside_it_keeps_its_own() {
    let rows = table(&run(&here("non_ref.vcf"), "tab"));
    let alleles: Vec<&str> = rows.keys().map(|k| k.1.as_str()).collect();
    assert!(
        !alleles.iter().any(|a| a.starts_with('<')),
        "a row for a placeholder allele: {alleles:?}"
    );
    // the records with only a placeholder are not annotated at all
    assert!(
        !rows
            .keys()
            .any(|k| k.0 == "chr10_1220_G/<NON_REF>" || k.0 == "chr10_1240_G/<*>"),
        "{:?}",
        rows.keys().collect::<Vec<_>>()
    );
    // the real ALT beside a placeholder is annotated: G at 1230 and C at 1250
    for (record, allele) in [
        ("chr10_1230_C/G/<NON_REF>", "G"),
        ("chr10_1250_G/C/<*>", "C"),
    ] {
        let row = rows
            .iter()
            .find(|(k, _)| k.0 == record && k.1 == allele && k.2 == "ENST00000000001")
            .unwrap_or_else(|| panic!("no row for {record} {allele}"))
            .1;
        assert!(row["Consequence"].contains("missense_variant"), "{row:?}");
    }
}
