//! VEP's `--total_length`, `--no_escape` and `--force_overwrite`, end to end (issue #139).
//!
//! The reporter's command line carries all three, and each was rejected by `fastvep annotate` with an
//! error, so a migrated command failed at once. They are accepted now, in VEP's spelling and in the
//! hyphenated one the rest of the CLI uses, and each does what VEP's does. What VEP's does was read
//! from its source and measured by running 116.3 (`vep_flags_compat_golden.rs` holds the output
//! against VEP's); this file is the behaviour, one flag at a time, on a RefSeq-style world.
//!
//! `fixtures/vep_flags_compat`: NB1 (`NM_090001.3`) and NB2 on the reverse strand (`NM_090002.2`),
//! coding with two exons each, NB3 (`NR_090003.1`, non-coding) and GENE A (`NM_090004.1`, one exon
//! and a symbol with a space).

use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};
use std::process::Command;

use fastvep_cli::pipeline::{run_annotate, AnnotateConfig, ColumnFlags, CompatFlags, PickFlags};
use tempfile::TempDir;

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("fixtures/vep_flags_compat")
        .join(name)
}

type Row = BTreeMap<String, String>;

fn table(text: &str) -> Vec<Row> {
    let mut lines = text.lines().filter(|l| !l.starts_with("##"));
    let header: Vec<String> = lines
        .next()
        .expect("a header")
        .trim_start_matches('#')
        .split('\t')
        .map(String::from)
        .collect();
    lines
        .map(|line| {
            let cells: Vec<&str> = line.split('\t').collect();
            assert_eq!(cells.len(), header.len(), "ragged row {line:?}");
            header
                .iter()
                .cloned()
                .zip(cells.into_iter().map(String::from))
                .collect()
        })
        .collect()
}

/// The world's tab output (every optional column) under `compat`.
fn tab(compat: CompatFlags) -> Vec<Row> {
    table(&run(compat, "tab"))
}

fn run(compat: CompatFlags, format: &str) -> String {
    let dir = TempDir::new().unwrap();
    // Copied, because the transcript cache a run writes goes beside its GFF3, and neither it nor
    // the FASTA's index belongs in the source tree.
    let (gff3, fasta) = (dir.path().join("genes.gff3"), dir.path().join("genome.fa"));
    std::fs::copy(fixture("genes.gff3"), &gff3).unwrap();
    std::fs::copy(fixture("genome.fa"), &fasta).unwrap();
    let out = dir.path().join(format!("out.{format}"));
    run_annotate(AnnotateConfig {
        input: fixture("in.vcf").to_string_lossy().into(),
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
            ..ColumnFlags::default()
        },
        compat,
        exclude: Default::default(),
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
    std::fs::read_to_string(&out).unwrap()
}

// --- --total_length ------------------------------------------------------------------------

#[test]
fn the_totals_are_the_lengths_of_the_cdna_the_cds_and_the_protein() {
    let rows = tab(CompatFlags {
        total_length: true,
        ..CompatFlags::default()
    });
    // (transcript, cDNA length, CDS length, protein length): the sums of the exons in
    // `genes.gff3`, the coding span and its codons without the stop
    let expected = [
        ("NM_090001.3", Some((600, 450, 149))),
        ("NM_090002.2", Some((600, 450, 149))),
        ("NM_090004.1", Some((300, 150, 49))),
        ("NR_090003.1", None),
    ];
    for (transcript, lengths) in expected {
        let mut exonic = 0;
        let mut coding = 0;
        for r in rows.iter().filter(|r| r["Feature"] == transcript) {
            let cdna = &r["cDNA_position"];
            if cdna == "-" {
                // an intron, or flank: no cell, and no total
                assert_eq!(r["CDS_position"], "-");
                assert_eq!(r["Protein_position"], "-");
                continue;
            }
            exonic += 1;
            let (cdna_len, cds_len, protein_len) = match lengths {
                Some(l) => (l.0, Some(l.1), Some(l.2)),
                None => (400, None, None),
            };
            assert!(cdna.ends_with(&format!("/{cdna_len}")), "{cdna}");
            match (r["CDS_position"].as_str(), cds_len) {
                ("-", _) => {} // a UTR: no CDS cell, and no CDS total
                (cell, Some(n)) => {
                    coding += 1;
                    assert!(cell.ends_with(&format!("/{n}")), "{cell}");
                    let protein = &r["Protein_position"];
                    assert!(
                        protein.ends_with(&format!("/{}", protein_len.unwrap())),
                        "{protein}"
                    );
                }
                (cell, None) => panic!("{transcript}: a CDS cell {cell}"),
            }
        }
        assert!(exonic >= 10, "{transcript}: {exonic} exonic rows");
        assert_eq!(coding >= 10, lengths.is_some(), "{transcript}: {coding}");
    }
}

#[test]
fn without_the_flag_the_cells_have_no_total() {
    for row in tab(CompatFlags::default()) {
        for col in ["cDNA_position", "CDS_position", "Protein_position"] {
            assert!(!row[col].contains('/'), "{col}: {}", row[col]);
        }
    }
}

#[test]
fn the_vcf_csq_carries_the_totals_too() {
    let vcf = run(
        CompatFlags {
            total_length: true,
            ..CompatFlags::default()
        },
        "vcf",
    );
    let format = vcf
        .lines()
        .find(|l| l.starts_with("##INFO=<ID=CSQ"))
        .and_then(|l| l.split("Format: ").nth(1))
        .map(|f| f.trim_end_matches("\">").split('|').collect::<Vec<_>>())
        .expect("a CSQ header");
    let at = |name: &str| format.iter().position(|f| *f == name).unwrap();
    let (cdna, cds, protein) = (
        at("cDNA_position"),
        at("CDS_position"),
        at("Protein_position"),
    );
    let (mut with_total, mut blank) = (0, 0);
    for line in vcf.lines().filter(|l| !l.starts_with('#')) {
        let info = line.split('\t').nth(7).unwrap();
        let Some(csq) = info.split(';').find_map(|kv| kv.strip_prefix("CSQ=")) else {
            continue;
        };
        for entry in csq.split(',') {
            let cells: Vec<&str> = entry.split('|').collect();
            if cells[cdna].is_empty() {
                // no cell and no total, and the CSQ leaves it empty, not `-`
                assert!(
                    cells[cds].is_empty() && cells[protein].is_empty(),
                    "{entry}"
                );
                blank += 1;
            } else {
                assert!(cells[cdna].contains('/'), "{}", cells[cdna]);
                with_total += 1;
            }
        }
    }
    assert!(with_total > 300 && blank > 500, "{with_total} {blank}");
}

// --- --no_escape ---------------------------------------------------------------------------

#[test]
fn a_synonymous_hgvsp_ends_in_an_equals_sign_only_with_no_escape() {
    let hgvsp = |compat: CompatFlags| -> BTreeSet<String> {
        tab(compat)
            .into_iter()
            .map(|r| r["HGVSp"].clone())
            .filter(|h| h.contains("%3D") || h.ends_with('='))
            .collect()
    };
    let escaped = hgvsp(CompatFlags::default());
    let plain = hgvsp(CompatFlags {
        no_escape: true,
        ..CompatFlags::default()
    });
    assert!(escaped.len() >= 4, "{escaped:?}");
    assert!(escaped.iter().all(|h| h.ends_with("%3D")), "{escaped:?}");
    assert!(plain.iter().all(|h| h.ends_with('=')), "{plain:?}");
    // the same notations, spelled the other way
    let unescaped: BTreeSet<String> = escaped.iter().map(|h| h.replace("%3D", "=")).collect();
    assert_eq!(unescaped, plain);
}

/// The CSQ entries of a VCF run as `(Feature, SYMBOL, HGVSp)`.
fn csq_of(compat: CompatFlags) -> Vec<(String, String, String)> {
    let vcf = run(compat, "vcf");
    let format: Vec<String> = vcf
        .lines()
        .find(|l| l.starts_with("##INFO=<ID=CSQ"))
        .and_then(|l| l.split("Format: ").nth(1))
        .map(|f| {
            f.trim_end_matches("\">")
                .split('|')
                .map(String::from)
                .collect()
        })
        .expect("a CSQ header");
    let at = |name: &str| format.iter().position(|f| f == name).unwrap();
    let (feature, symbol, hgvsp) = (at("Feature"), at("SYMBOL"), at("HGVSp"));
    let mut out = Vec::new();
    for line in vcf.lines().filter(|l| !l.starts_with('#')) {
        let info = line.split('\t').nth(7).unwrap();
        if let Some(csq) = info.split(';').find_map(|kv| kv.strip_prefix("CSQ=")) {
            for entry in csq.split(',') {
                let cells: Vec<&str> = entry.split('|').collect();
                out.push((
                    cells[feature].into(),
                    cells[symbol].into(),
                    cells[hgvsp].into(),
                ));
            }
        }
    }
    out
}

#[test]
fn in_a_vcf_the_space_of_a_symbol_is_percent_20_and_with_no_escape_an_underscore() {
    // `Name=GENE A`, VEP 116.3: `GENE%20A`, and `GENE_A` with `--no_escape`
    let symbols = |compat| -> BTreeSet<String> {
        csq_of(compat)
            .into_iter()
            .filter(|(f, _, _)| f == "NM_090004.1")
            .map(|(_, s, _)| s)
            .collect()
    };
    assert_eq!(
        symbols(CompatFlags::default()),
        BTreeSet::from(["GENE%20A".to_string()])
    );
    assert_eq!(
        symbols(CompatFlags {
            no_escape: true,
            ..CompatFlags::default()
        }),
        BTreeSet::from(["GENE_A".to_string()])
    );
    // the tab file never escaped it
    assert!(tab(CompatFlags::default())
        .iter()
        .any(|r| r["SYMBOL"] == "GENE A"));
}

#[test]
fn in_a_vcf_the_equals_sign_of_a_synonymous_hgvsp_follows_no_escape_too() {
    let hgvsp = |compat| -> BTreeSet<String> {
        csq_of(compat)
            .into_iter()
            .map(|(_, _, h)| h)
            .filter(|h| h.contains('='))
            .collect()
    };
    assert!(
        hgvsp(CompatFlags::default()).is_empty(),
        "an `=` that is not escaped"
    );
    let plain = hgvsp(CompatFlags {
        no_escape: true,
        ..CompatFlags::default()
    });
    assert!(
        plain.len() >= 4 && plain.iter().all(|h| h.ends_with('=')),
        "{plain:?}"
    );
}

// --- the command line ----------------------------------------------------------------------

fn fastvep(args: &[&str]) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_fastvep"))
        .args(args)
        .output()
        .expect("the fastvep binary runs")
}

#[test]
fn the_three_flags_are_accepted_in_vep_s_spelling_and_in_the_hyphenated_one() {
    let dir = TempDir::new().unwrap();
    let (gff3, fasta, out) = (
        dir.path().join("genes.gff3"),
        dir.path().join("genome.fa"),
        dir.path().join("out.tab"),
    );
    std::fs::copy(fixture("genes.gff3"), &gff3).unwrap();
    std::fs::copy(fixture("genome.fa"), &fasta).unwrap();
    let vcf = fixture("in.vcf");
    for flag in [
        "--total_length",
        "--total-length",
        "--no_escape",
        "--no-escape",
        "--force_overwrite",
        "--force-overwrite",
    ] {
        let result = fastvep(&[
            "annotate",
            "-i",
            vcf.to_str().unwrap(),
            "--gff3",
            gff3.to_str().unwrap(),
            "--fasta",
            fasta.to_str().unwrap(),
            "--output-format",
            "tab",
            "-o",
            out.to_str().unwrap(),
            "--no-progress",
            "--hgvs",
            flag,
        ]);
        assert!(
            result.status.success(),
            "{flag}: {}",
            String::from_utf8_lossy(&result.stderr)
        );
        let rows = table(&std::fs::read_to_string(&out).unwrap());
        assert!(rows.len() > 1000, "{flag}: {} rows", rows.len());
        // what each flag does shows in the file it wrote
        let cell = |col: &str| rows.iter().any(|r| r[col].contains('/'));
        let hgvsp: Vec<&str> = rows.iter().map(|r| r["HGVSp"].as_str()).collect();
        match flag {
            "--total_length" | "--total-length" => assert!(cell("cDNA_position"), "{flag}"),
            "--no_escape" | "--no-escape" => assert!(
                hgvsp.iter().any(|h| h.ends_with('=')) && !hgvsp.iter().any(|h| h.contains("%3D")),
                "{flag}"
            ),
            _ => {
                assert!(!cell("cDNA_position"), "{flag}");
                assert!(hgvsp.iter().any(|h| h.ends_with("%3D")), "{flag}");
            }
        }
    }
}

#[test]
fn force_overwrite_changes_nothing_because_fastvep_always_writes_over_its_output() {
    let dir = TempDir::new().unwrap();
    let (gff3, fasta, out) = (
        dir.path().join("genes.gff3"),
        dir.path().join("genome.fa"),
        dir.path().join("out.vcf"),
    );
    std::fs::copy(fixture("genes.gff3"), &gff3).unwrap();
    std::fs::copy(fixture("genome.fa"), &fasta).unwrap();
    let vcf = fixture("in.vcf");
    let annotate = |extra: &[&str]| -> String {
        let mut args = vec![
            "annotate",
            "-i",
            vcf.to_str().unwrap(),
            "--gff3",
            gff3.to_str().unwrap(),
            "--fasta",
            fasta.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
            "--no-progress",
        ];
        args.extend_from_slice(extra);
        // a file from an earlier run is there, which VEP refuses to write over without the flag
        std::fs::write(&out, "an earlier run\n").unwrap();
        let result = fastvep(&args);
        assert!(
            result.status.success(),
            "{extra:?}: {}",
            String::from_utf8_lossy(&result.stderr)
        );
        std::fs::read_to_string(&out).unwrap()
    };
    let without = annotate(&[]);
    assert!(without.starts_with("##fileformat"), "{without:.60}");
    assert_eq!(annotate(&["--force_overwrite"]), without);
    assert_eq!(annotate(&["--force-overwrite"]), without);
}

#[test]
fn the_help_says_what_each_flag_does() {
    let help = String::from_utf8(fastvep(&["annotate", "--help"]).stdout).unwrap();
    for flag in ["--total-length", "--no-escape", "--force-overwrite"] {
        assert!(help.contains(flag), "{flag} is not in the help");
    }
    assert!(help.contains("`vep` command line"), "{help}");
}
