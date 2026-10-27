//! The known-variant flags as a user types them: the `fastvep annotate` command,
//! run as a program over the generated world of
//! `fastvep-cache/fixtures/vep_known_columns`.
//!
//! `known_variant_columns_golden.rs` holds the columns to VEP's; this holds what
//! `main.rs` does with the command line, which no call of `run_annotate` can see: that
//! VEP's underscore spellings (`--af_1kg`, `--check_existing`, `--clin_sig_allele`) are
//! the same flags as the hyphen ones, that `--everything` brings the frequency flags
//! with it only where there is a cache to read them from, that a flag which reads the
//! cache without one stops the run, and that the VEP 116 flags this does not have are
//! refused and not taken for something else.

use std::path::{Path, PathBuf};
use std::process::{Command, Output};

use tempfile::TempDir;

fn world() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../fastvep-cache/fixtures/vep_known_columns")
}

fn cache() -> String {
    world()
        .join("cache/homo_sapiens/116_GRCh38")
        .to_string_lossy()
        .into()
}

/// `fastvep annotate` over the world, tab output, with `extra` after the flags every run has.
fn annotate(dir: &Path, extra: &[&str]) -> (Output, String) {
    annotate_as(dir, "tab", extra)
}

fn annotate_as(dir: &Path, format: &str, extra: &[&str]) -> (Output, String) {
    let out = dir.join(format!("out.{format}"));
    let output = Command::new(env!("CARGO_BIN_EXE_fastvep"))
        .arg("annotate")
        .arg("-i")
        .arg(world().join("in.vcf"))
        .arg("-o")
        .arg(&out)
        .arg("--gff3")
        .arg(world().join("genes.gff3"))
        .args(["--output-format", format, "--no-progress"])
        .args(extra)
        .output()
        .expect("the fastvep binary runs");
    let tab = std::fs::read_to_string(&out).unwrap_or_default();
    (output, tab)
}

fn header(tab: &str) -> Vec<String> {
    tab.lines()
        .find(|l| l.starts_with("#Uploaded_variation"))
        .map(|l| {
            l.trim_start_matches('#')
                .split('\t')
                .map(String::from)
                .collect()
        })
        .unwrap_or_default()
}

fn stderr(output: &Output) -> String {
    String::from_utf8_lossy(&output.stderr).into()
}

/// VEP's spelling of a flag is the same flag. Each pair is run and its two outputs are
/// the same file, and each differs from a run without the flags, so that an alias that
/// parsed and did nothing would not pass for one that worked.
#[test]
fn vep_s_underscore_spellings_are_the_hyphen_flags() {
    let cache = cache();
    let cases: Vec<(Vec<&str>, Vec<&str>)> = vec![
        (
            vec!["--af_1kg", "--check_existing"],
            vec!["--af-1kg", "--check-existing"],
        ),
        (
            vec!["--af_gnomade", "--af_gnomadg"],
            vec!["--af-gnomade", "--af-gnomadg"],
        ),
        (vec!["--af_gnomad"], vec!["--af-gnomade"]),
        (vec!["--af-gnomad"], vec!["--af-gnomade"]),
        (vec!["--max_af"], vec!["--max-af"]),
        (
            vec!["--check_existing", "--clin_sig_allele", "0"],
            vec!["--check-existing", "--clin-sig-allele", "0"],
        ),
        (
            vec!["--show_ref_allele", "--variant_class"],
            vec!["--show-ref-allele", "--variant-class"],
        ),
        (vec!["--flag_pick"], vec!["--flag-pick"]),
        (vec!["--flag_pick_allele"], vec!["--flag-pick-allele"]),
        (
            vec!["--flag_pick_allele_gene"],
            vec!["--flag-pick-allele-gene"],
        ),
        (vec!["--pick_allele"], vec!["--pick-allele"]),
        (vec!["--pick_allele_gene"], vec!["--pick-allele-gene"]),
        (
            vec!["--pick", "--pick_order", "rank,canonical"],
            vec!["--pick", "--pick-order", "rank,canonical"],
        ),
        (vec!["--buffer_size", "100"], vec!["--buffer-size", "100"]),
    ];
    let dir = TempDir::new().unwrap();
    let (output, plain) = annotate(dir.path(), &["--cache-dir", &cache]);
    assert!(output.status.success(), "{}", stderr(&output));

    for (underscore, hyphen) in cases {
        let with = |flags: &[&str]| {
            let mut args = vec!["--cache-dir", cache.as_str()];
            args.extend_from_slice(flags);
            args.iter().map(|a| a.to_string()).collect::<Vec<_>>()
        };
        let (u, h) = (with(&underscore), with(&hyphen));
        let (u_out, u_tab) = annotate(
            dir.path(),
            &u.iter().map(String::as_str).collect::<Vec<_>>(),
        );
        assert!(u_out.status.success(), "{underscore:?}: {}", stderr(&u_out));
        let (h_out, h_tab) = annotate(
            dir.path(),
            &h.iter().map(String::as_str).collect::<Vec<_>>(),
        );
        assert!(h_out.status.success(), "{hyphen:?}: {}", stderr(&h_out));
        assert_eq!(u_tab, h_tab, "{underscore:?} is not {hyphen:?}");
        // `--buffer_size` and `--pick_order` change no cell of this input, which is not what is
        // under test; the others change the columns or the rows
        if !matches!(underscore[0], "--buffer_size" | "--pick") {
            assert_ne!(u_tab, plain, "{underscore:?} did nothing");
        }
    }
}

/// `--everything` is VEP's set. With a cache it has the frequency flags in it, as in VEP; without
/// one it is the set fastVEP has always had, and not an error for a flag the user never typed.
#[test]
fn everything_brings_the_frequency_columns_where_there_is_a_cache_to_read_them_from() {
    let dir = TempDir::new().unwrap();

    let (output, tab) = annotate(dir.path(), &["--everything", "--cache-dir", &cache()]);
    assert!(output.status.success(), "{}", stderr(&output));
    let columns = header(&tab);
    for want in [
        "AF",
        "AFR_AF",
        "gnomADe_AF",
        "gnomADg_SAS_AF",
        "MAX_AF",
        "MAX_AF_POPS",
        "CLIN_SIG",
        "PUBMED",
    ] {
        assert!(
            columns.iter().any(|c| c == want),
            "--everything with a cache has no {want}: {columns:?}"
        );
    }

    let (output, tab) = annotate(dir.path(), &["--everything"]);
    assert!(output.status.success(), "{}", stderr(&output));
    let columns = header(&tab);
    assert!(
        !columns
            .iter()
            .any(|c| c == "AF" || c == "MAX_AF" || c == "CLIN_SIG"),
        "--everything without a cache wrote known-variant columns: {columns:?}"
    );
    assert!(columns.iter().any(|c| c == "VARIANT_CLASS"), "{columns:?}");
}

#[test]
fn a_flag_that_reads_the_cache_stops_a_run_that_has_none() {
    let dir = TempDir::new().unwrap();
    for flag in [
        "--af",
        "--af-1kg",
        "--af_gnomadg",
        "--max-af",
        "--pubmed",
        "--check_existing",
    ] {
        let (output, _) = annotate(dir.path(), &[flag]);
        assert!(!output.status.success(), "{flag} without a cache ran");
        let message = stderr(&output);
        assert!(message.contains("--cache-dir"), "{flag}: {message}");
    }
}

/// VEP 116 has no `--af_esp` and no `--af_exac`: its `Config.pm` lists `af`, `af_1kg`,
/// `af_gnomade`, `af_gnomad`, `af_gnomadg` and `max_af`. They are not taken for another flag.
#[test]
fn the_frequency_flags_vep_116_does_not_have_are_refused() {
    let dir = TempDir::new().unwrap();
    for flag in ["--af_esp", "--af-esp", "--af_exac", "--af-exac"] {
        let (output, _) = annotate(dir.path(), &["--cache-dir", &cache(), flag]);
        assert!(!output.status.success(), "{flag} was accepted");
        assert!(
            stderr(&output).contains("unexpected argument"),
            "{flag}: {}",
            stderr(&output)
        );
    }
}

/// A cache without the columns of a flag says so, once, where it is opened, and the columns are
/// empty as they are in VEP.
#[test]
fn a_cache_without_the_columns_of_a_flag_says_so() {
    let dir = TempDir::new().unwrap();
    let thin = dir.path().join("cache/homo_sapiens/116_GRCh38");
    std::fs::create_dir_all(thin.join("10")).unwrap();
    // the world's cache with its frequency columns removed from info.txt and the files: simplest is
    // a cache whose info.txt names only the first fifteen columns, which is a valid prefix
    // of every line (the reader looks columns up by position)
    let info =
        std::fs::read_to_string(world().join("cache/homo_sapiens/116_GRCh38/info.txt")).unwrap();
    let cut: String = info
        .lines()
        .map(|l| match l.strip_prefix("variation_cols\t") {
            Some(cols) => format!(
                "variation_cols\t{}",
                cols.split(',').take(15).collect::<Vec<_>>().join(",")
            ),
            None => l.to_string(),
        })
        .collect::<Vec<_>>()
        .join("\n");
    std::fs::write(thin.join("info.txt"), cut + "\n").unwrap();
    for name in ["all_vars.gz", "all_vars.gz.csi"] {
        std::fs::copy(
            world().join("cache/homo_sapiens/116_GRCh38/10").join(name),
            thin.join("10").join(name),
        )
        .unwrap();
    }
    let (output, tab) = annotate(
        dir.path(),
        &[
            "--cache-dir",
            &thin.to_string_lossy(),
            "--af-gnomadg",
            "--af-1kg",
        ],
    );
    assert!(output.status.success(), "{}", stderr(&output));
    let message = stderr(&output);
    assert!(
        message.contains("--af-gnomadg asks for the gnomAD genomes frequencies"),
        "no warning about the missing gnomAD genomes columns: {message}"
    );
    assert!(
        message.contains("--af-1kg asks for the 1000 Genomes frequencies"),
        "no warning about the missing 1000 Genomes columns: {message}"
    );
    // the columns are there, empty, as VEP writes them
    assert!(header(&tab).iter().any(|c| c == "gnomADg_AF"));
    let at = header(&tab).iter().position(|c| c == "AF").unwrap();
    assert!(tab
        .lines()
        .filter(|l| !l.starts_with('#'))
        .all(|l| l.split('\t').nth(at) == Some("-")));
}

/// JSON has no known-variant block, so a flag for one says it is ignored instead of being ignored.
#[test]
fn json_output_says_it_carries_no_known_variant_columns() {
    let dir = TempDir::new().unwrap();
    let (output, json) = annotate_as(dir.path(), "json", &["--cache-dir", &cache(), "--af-1kg"]);
    assert!(output.status.success(), "{}", stderr(&output));
    assert!(
        stderr(&output).contains(
            "--output-format json does not carry the known-variant columns; ignoring --af-1kg"
        ),
        "{}",
        stderr(&output)
    );
    assert!(!json.contains("\"AFR_AF\""), "JSON gained a field");
    // only the flag that was typed is named: `--cache-dir` is `--check-existing` for the formats
    // that have the columns, and JSON is not one of them
    assert!(
        !stderr(&output).contains("--check-existing"),
        "{}",
        stderr(&output)
    );

    let (output, _) = annotate_as(dir.path(), "json", &["--cache-dir", &cache()]);
    assert!(output.status.success(), "{}", stderr(&output));
    assert!(!stderr(&output).contains("ignoring"), "{}", stderr(&output));
}

/// Naming a cache is VEP's `--check_existing`: the columns come without the flag, in VEP's order,
/// and the frequencies, which VEP only writes under their flags, do not.
#[test]
fn a_cache_alone_adds_the_check_existing_columns_and_no_frequency() {
    let dir = TempDir::new().unwrap();
    let cache = cache();
    let (output, alone) = annotate(dir.path(), &["--cache-dir", &cache]);
    assert!(output.status.success(), "{}", stderr(&output));
    let (output, typed) = annotate(dir.path(), &["--cache-dir", &cache, "--check_existing"]);
    assert!(output.status.success(), "{}", stderr(&output));
    assert_eq!(alone, typed, "--cache-dir alone is --check_existing");
    let columns = header(&alone);
    for name in ["Existing_variation", "CLIN_SIG", "SOMATIC", "PHENO"] {
        assert!(columns.iter().any(|c| c == name), "{name} in {columns:?}");
    }
    for name in ["AF", "AFR_AF", "MAX_AF", "PUBMED"] {
        assert!(!columns.iter().any(|c| c == name), "{name} in {columns:?}");
    }
    assert!(
        !stderr(&output).contains("ignoring"),
        "no warning for a flag that was not typed: {}",
        stderr(&output)
    );
}
