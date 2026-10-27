//! The known-variant columns, end to end, against Ensembl VEP 116.3 run as a
//! program: `AF`, `AFR_AF` .. `SAS_AF`, `gnomADe_*`, `gnomADg_*`, `MAX_AF`,
//! `MAX_AF_POPS`, `CLIN_SIG`, `SOMATIC`, `PHENO` and `PUBMED`. Issue #139: the
//! reporter runs `--af_1kg --check_existing`.
//!
//! `fastvep-cache/fixtures/vep_known_columns` is a cache with the 42 columns of
//! the real release-116 human cache, 446 input records (50 written by hand, each
//! a corner of what `add_colocated_variant_info` and
//! `add_colocated_frequency_data` do, and 396 generated, each with a random
//! payload) and the answers of VEP for eight sets of flags, one golden file each.
//! `regen.sh` there runs `vep --tab --cache --offline` from the pinned image and
//! keeps what it printed in these columns; nothing in a golden file was worked out
//! here.
//!
//! What this holds to account is the whole path of the program: that a flag
//! reaches the lookup, that the columns come in VEP's order, that a row for a
//! variant with no transcript (an intergenic one) carries them too, and that the
//! cells are VEP's, quirks and all (the `%.4f` of `AF`, the interpolated
//! frequency, `MAX_AF_POPS` repeating a population, the allele-specific
//! `CLIN_SIG`). The ports of the two Perl subs are unit-tested in
//! `fastvep-annotate/src/colocated.rs`; that VEP agrees with them is shown here.
//!
//! `CLIN_SIG` is compared as a set of its `;`-separated parts: when two known
//! variants give different significances for the allele VEP joins them in the
//! order of a Perl hash, which differs from one run to the next.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use fastvep_cli::pipeline::{run_annotate, AnnotateConfig, ColumnFlags, PickFlags};
use fastvep_io::known::KnownFlags;
use tempfile::TempDir;

fn world() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../fastvep-cache/fixtures/vep_known_columns")
}

fn config(input: &Path, out: &Path, format: &str, known: KnownFlags) -> AnnotateConfig {
    let world = world();
    AnnotateConfig {
        input: input.to_string_lossy().into(),
        output: out.to_string_lossy().into(),
        gff3: vec![world.join("genes.gff3").to_string_lossy().into()],
        fasta: None,
        output_format: format.into(),
        pick: PickFlags::default(),
        hgvs: false,
        columns: ColumnFlags {
            known,
            ..Default::default()
        },
        distance: 5000,
        cache_dir: Some(
            world
                .join("cache/homo_sapiens/116_GRCh38")
                .to_string_lossy()
                .into(),
        ),
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
    }
}

fn run(format: &str, known: KnownFlags) -> String {
    let dir = TempDir::new().unwrap();
    let out = dir.path().join(format!("out.{format}"));
    run_annotate(config(&world().join("in.vcf"), &out, format, known)).expect("annotation");
    std::fs::read_to_string(&out).unwrap()
}

/// The flags of each golden file, as `regen.sh` gave them to VEP.
fn sets() -> Vec<(&'static str, KnownFlags)> {
    let none = KnownFlags::default();
    vec![
        (
            "reporter",
            KnownFlags {
                af_1kg: true,
                check_existing: true,
                ..none
            },
        ),
        ("af", KnownFlags { af: true, ..none }),
        (
            "gnomad",
            KnownFlags {
                af_gnomade: true,
                af_gnomadg: true,
                ..none
            },
        ),
        (
            "max_af",
            KnownFlags {
                max_af: true,
                ..none
            },
        ),
        (
            "pubmed",
            KnownFlags {
                pubmed: true,
                ..none
            },
        ),
        (
            "existing",
            KnownFlags {
                check_existing: true,
                ..none
            },
        ),
        (
            "csa0",
            KnownFlags {
                check_existing: true,
                pubmed: true,
                no_clin_sig_allele: true,
                ..none
            },
        ),
        (
            "all",
            KnownFlags {
                af: true,
                af_1kg: true,
                af_gnomade: true,
                af_gnomadg: true,
                max_af: true,
                pubmed: true,
                check_existing: true,
                ..none
            },
        ),
    ]
}

/// The columns of `Existing_variation` and of the known variants: the ones a
/// golden file is about.
fn is_known_column(name: &str) -> bool {
    name == "Existing_variation" || fastvep_io::known::KnownColumn::from_name(name).is_some()
}

type Row = (String, String, String);

/// The golden file: its columns in VEP's order, and for every (record, allele,
/// feature) the cells that are not `-`.
fn golden(name: &str) -> (Vec<String>, BTreeMap<Row, BTreeMap<String, String>>) {
    let text = std::fs::read_to_string(world().join(format!("golden_{name}.tsv"))).unwrap();
    let mut columns = Vec::new();
    let mut rows = BTreeMap::new();
    for line in text.lines() {
        if let Some(rest) = line.strip_prefix("#columns\t") {
            columns = rest.split('\t').map(String::from).collect();
        } else if line.starts_with('#') {
            continue;
        } else {
            let f: Vec<&str> = line.split('\t').collect();
            let cells = f[3..]
                .iter()
                .map(|c| {
                    let (k, v) = c.split_once('=').expect("COLUMN=value");
                    (k.to_string(), v.to_string())
                })
                .collect();
            rows.insert((f[0].into(), f[1].into(), f[2].into()), cells);
        }
    }
    (columns, rows)
}

/// What fastVEP's tab file says in the known columns, in the golden file's shape.
fn parse_tab(tab: &str) -> (Vec<String>, BTreeMap<Row, BTreeMap<String, String>>) {
    let header: Vec<&str> = tab
        .lines()
        .find(|l| l.starts_with("#Uploaded_variation"))
        .unwrap_or_else(|| panic!("no column header in:\n{tab}"))
        .trim_start_matches('#')
        .split('\t')
        .collect();
    let at = |name: &str| header.iter().position(|c| *c == name).unwrap();
    let (id, allele, feature) = (at("Uploaded_variation"), at("Allele"), at("Feature"));
    let known: Vec<(usize, &str)> = header
        .iter()
        .enumerate()
        .filter(|(_, c)| is_known_column(c))
        .map(|(i, c)| (i, *c))
        .collect();
    let mut rows = BTreeMap::new();
    for line in tab.lines().filter(|l| !l.starts_with('#')) {
        let cells: Vec<&str> = line.split('\t').collect();
        let values = known
            .iter()
            .filter(|(i, _)| cells[*i] != "-")
            .map(|(i, c)| (c.to_string(), cells[*i].to_string()))
            .collect();
        // VEP upper-cases the alleles of the input. fastVEP writes a lower-case one
        // as it was given, which is not what these columns are about.
        rows.insert(
            (
                cells[id].to_string(),
                cells[allele].to_uppercase(),
                cells[feature].to_string(),
            ),
            values,
        );
    }
    (known.iter().map(|(_, c)| c.to_string()).collect(), rows)
}

fn normalised(column: &str, value: &str) -> String {
    if column == "CLIN_SIG" {
        let mut parts: Vec<&str> = value.split(';').collect();
        parts.sort_unstable();
        parts.join(";")
    } else {
        value.to_string()
    }
}

#[test]
fn every_set_of_flags_writes_the_cells_vep_wrote_in_the_order_it_wrote_them() {
    for (name, flags) in sets() {
        let (want_columns, want) = golden(name);
        let (got_columns, got) = parse_tab(&run("tab", flags));

        assert_eq!(
            got_columns, want_columns,
            "{name}: the known-variant columns, or their order, differ from VEP's"
        );
        assert_eq!(got.len(), want.len(), "{name}: rows written");

        let wrong: Vec<String> = want
            .iter()
            .filter_map(|(row, cells)| {
                let ours = got.get(row)?;
                let differing: Vec<String> = cells
                    .keys()
                    .chain(ours.keys())
                    .collect::<std::collections::BTreeSet<_>>()
                    .into_iter()
                    .filter(|c| {
                        normalised(c, cells.get(*c).map_or("-", String::as_str))
                            != normalised(c, ours.get(*c).map_or("-", String::as_str))
                    })
                    .map(|c| {
                        format!(
                            "{c}: VEP {:?}, here {:?}",
                            cells.get(c.as_str()),
                            ours.get(c.as_str())
                        )
                    })
                    .collect();
                (!differing.is_empty()).then(|| format!("{row:?} {}", differing.join("; ")))
            })
            .collect();
        let missing: Vec<&Row> = want.keys().filter(|k| !got.contains_key(*k)).collect();
        assert!(
            missing.is_empty() && wrong.is_empty(),
            "{name}: {} of {} rows differ from VEP's answer ({} not written), first ten:\n{}",
            wrong.len() + missing.len(),
            want.len(),
            missing.len(),
            wrong
                .iter()
                .take(10)
                .cloned()
                .collect::<Vec<_>>()
                .join("\n")
        );
    }
}

/// The same cells reach the VCF's `CSQ`, with `,` written `&` and `;` written `%3B`
/// as VEP's VCF writer does, and the extra fields in VEP's relative order.
#[test]
fn the_csq_carries_the_cells_of_the_tab_file() {
    let flags = sets().into_iter().find(|(n, _)| *n == "all").unwrap().1;
    let (_, tab_rows) = parse_tab(&run("tab", flags));
    let vcf = run("vcf", flags);

    let format = vcf
        .lines()
        .find(|l| l.starts_with("##INFO=<ID=CSQ"))
        .expect("a CSQ header");
    let fields: Vec<&str> = format
        .split("Format: ")
        .nth(1)
        .unwrap()
        .trim_end_matches("\">")
        .split('|')
        .collect();
    // the five fields the CSQ always had come first, where they were; the rest follow AF
    let position = |name: &str| fields.iter().position(|f| *f == name).unwrap();
    assert_eq!(position("AF") + 1, position("AFR_AF"));
    assert!(position("gnomADe_AF") < position("gnomADg_AF"));
    assert!(position("MAX_AF_POPS") + 1 == position("CLIN_SIG"));
    assert!(position("PHENO") < position("PUBMED"));

    let mut compared = 0;
    for line in vcf.lines().filter(|l| !l.starts_with('#')) {
        let cols: Vec<&str> = line.split('\t').collect();
        let id = cols[2];
        let csq = cols[7]
            .split(';')
            .find_map(|i| i.strip_prefix("CSQ="))
            .expect("a CSQ");
        for entry in csq.split(',') {
            let values: Vec<&str> = entry.split('|').collect();
            let field = |name: &str| values[position(name)];
            let row = (id.to_string(), field("Allele").to_uppercase(), {
                let feature = field("Feature");
                if feature.is_empty() { "-" } else { feature }.to_string()
            });
            let tab = &tab_rows[&row];
            for column in fastvep_io::known::KnownColumn::ALL {
                let name = column.name();
                let want = tab.get(name).map_or(String::new(), |v| {
                    v.replace(',', "&").replace(';', "%3B").replace(' ', "_")
                });
                assert_eq!(field(name), want, "{row:?} {name}");
                compared += 1;
            }
        }
    }
    assert!(compared > 10_000, "only {compared} cells compared");
}

/// A cache and no flag at all is VEP's `--check_existing` and nothing more: `CLIN_SIG`, `SOMATIC` and
/// `PHENO` are the cells of VEP's `--check_existing` run, and the frequency and PubMed columns, which
/// VEP only writes under their own flags, are empty. The CSQ keeps its 49 fields.
///
/// It used to carry a summary of the whole variant in these five fields (`AF` was the first
/// frequency of any known variant of the site, gnomAD's where it had one), which for a record with
/// several ALTs put one allele's number on every allele's row.
#[test]
fn a_cache_alone_is_check_existing_and_asks_for_no_frequency() {
    let none = KnownFlags::default();
    let (want_columns, want) = golden("existing");
    let (got_columns, got) = parse_tab(&run("tab", none));
    assert_eq!(got_columns, want_columns, "the columns of --check_existing");
    assert_eq!(got.len(), want.len());
    for (row, cells) in &want {
        let ours = &got[row];
        for column in ["CLIN_SIG", "SOMATIC", "PHENO", "Existing_variation"] {
            assert_eq!(
                normalised(column, cells.get(column).map_or("-", String::as_str)),
                normalised(column, ours.get(column).map_or("-", String::as_str)),
                "{row:?} {column}"
            );
        }
    }

    let vcf = run("vcf", none);
    let format = vcf
        .lines()
        .find(|l| l.starts_with("##INFO=<ID=CSQ"))
        .unwrap();
    let fields: Vec<&str> = format
        .split("Format: ")
        .nth(1)
        .unwrap()
        .trim_end_matches("\">")
        .split('|')
        .collect();
    assert_eq!(fields.len(), 49, "the CSQ gained or lost a field");
    assert!(!fields.contains(&"AFR_AF") && !fields.contains(&"MAX_AF"));
    let position = |name: &str| fields.iter().position(|f| *f == name).unwrap();

    let (mut with_clin_sig, mut entries) = (0, 0);
    for line in vcf.lines().filter(|l| !l.starts_with('#')) {
        let csq = line
            .split('\t')
            .nth(7)
            .unwrap()
            .split(';')
            .find_map(|i| i.strip_prefix("CSQ="))
            .expect("a CSQ");
        for entry in csq.split(',') {
            let values: Vec<&str> = entry.split('|').collect();
            entries += 1;
            assert_eq!(values[position("AF")], "", "no AF without --af: {entry}");
            assert_eq!(values[position("PUBMED")], "", "no PUBMED without --pubmed");
            with_clin_sig += usize::from(!values[position("CLIN_SIG")].is_empty());
        }
    }
    assert!(
        entries > 500 && with_clin_sig > 20,
        "{with_clin_sig} of {entries}"
    );

    // and with --af it is VEP's: t01 is 0.123456, which VEP writes `%.4f`
    let with_af = run(
        "vcf",
        KnownFlags {
            af: true,
            ..Default::default()
        },
    );
    let line = with_af
        .lines()
        .find(|l| l.split('\t').nth(2) == Some("t01"))
        .unwrap();
    assert!(line.contains("|0.1235|"), "{line}");
}

#[test]
fn a_flag_that_reads_the_cache_needs_one() {
    let dir = TempDir::new().unwrap();
    let out = dir.path().join("out.tab");
    for known in [
        KnownFlags {
            af_1kg: true,
            ..Default::default()
        },
        KnownFlags {
            check_existing: true,
            ..Default::default()
        },
        KnownFlags {
            pubmed: true,
            ..Default::default()
        },
    ] {
        let flags = known.requested().join(", ");
        let err = run_annotate(AnnotateConfig {
            cache_dir: None,
            ..config(&world().join("in.vcf"), &out, "tab", known)
        })
        .expect_err("no cache, no known variants");
        let message = format!("{err:#}");
        assert!(message.contains(&flags), "{message}");
        assert!(message.contains("--cache-dir"), "{message}");
    }
}
