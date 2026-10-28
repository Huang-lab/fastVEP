//! The flags of Ensembl VEP that `fastvep annotate` refuses, as a user who pastes a `vep` command
//! line meets them: the program run with each one.
//!
//! Each is `error: unexpected argument '--flag' found`, as `docs/VEP_MIGRATION.md` says, and clap
//! used to add a guess at what was meant that was wrong for most of them: `--assembly` was
//! "similar" to `--symbol`, `--json` to `--sa-only`, `--bam` to `--acmg`, `--no_stats` to
//! `--no_escape`, `--protein` to `--proband`, `--uniprot` to `--input`, `--use_transcript_ref` to
//! `--transcript-cache`, and the usage line under the error had the guess in it. Now the tip says
//! what there is instead, and the usage is the plain one.

use std::process::{Command, Output};

/// The flags the migration page lists as refused, each of which was run against the binary.
const REFUSED: &[&str] = &[
    "cache",
    "offline",
    "dir_cache",
    "refseq",
    "merged",
    "tab",
    "vcf",
    "json",
    "assembly",
    "species",
    "use_transcript_ref",
    "bam",
    "fields",
    "plugin",
    "custom",
    "fork",
    "no_stats",
    "quiet",
    "sift",
    "polyphen",
    "mane",
    "tsl",
    "appris",
    "ccds",
    "protein",
    "domains",
    "uniprot",
    "gene_phenotype",
    "regulatory",
    "nearest",
    "shift_3prime",
    "no_intergenic",
    "coding_only",
    "minimal",
    "allele_number",
    "individual",
    "vcf_info_field",
    "terms",
    "af_esp",
    "af_exac",
];

fn annotate(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_fastvep"))
        .args(["annotate", "-i", "in.vcf"])
        .args(args)
        .output()
        .expect("the fastvep binary runs")
}

fn stderr(output: &Output) -> String {
    String::from_utf8_lossy(&output.stderr).into()
}

#[test]
fn every_refused_flag_of_vep_says_what_there_is_instead_of_a_wrong_guess() {
    for name in REFUSED {
        for typed in [
            format!("--{name}"),
            format!("--{}", name.replace('_', "-")),
            format!("--{name}=1"),
        ] {
            let output = annotate(&[&typed]);
            let err = stderr(&output);
            let flag = typed.split('=').next().unwrap();
            assert_eq!(output.status.code(), Some(2), "{typed}: {err}");
            assert!(
                err.starts_with(&format!("error: unexpected argument '{flag}' found")),
                "{typed}: {err}"
            );
            assert!(err.contains("tip: "), "{typed}: {err}");
            assert!(
                err.contains("docs/VEP_MIGRATION.md"),
                "{typed}: the tip does not point at the page: {err}"
            );
            assert!(
                !err.contains("a similar argument exists"),
                "{typed}: clap's guess is back: {err}"
            );
            // and the usage line has no flag in it but the one every run needs
            let usage = err.lines().find(|l| l.starts_with("Usage:")).unwrap_or("");
            assert_eq!(
                usage.matches("--").count(),
                1,
                "{typed}: the usage line has a flag in it: {usage:?}"
            );
        }
    }
}

#[test]
fn the_tips_name_the_flag_to_use_in_the_cases_that_a_user_will_meet_first() {
    let tip = |flag: &str| stderr(&annotate(&[flag]));
    assert!(tip("--json").contains("--output-format"));
    assert!(tip("--tab").contains("--output-format"));
    assert!(tip("--refseq").contains("--gff3"));
    assert!(tip("--cache").contains("--cache-dir"));
    assert!(tip("--dir_cache").contains("--cache-dir"));
    assert!(tip("--fork").contains("RAYON_NUM_THREADS"));
    assert!(tip("--quiet").contains("--no-progress"));
    assert!(tip("--af_esp").contains("--af_gnomade"));
    assert!(tip("--plugin").contains("sa-build"));
}

/// A flag that is not one of VEP's is not touched, so a typo of one fastVEP has keeps the
/// suggestion that is right for it.
#[test]
fn a_flag_that_is_not_vep_s_keeps_clap_s_own_message() {
    let typo = stderr(&annotate(&["--hgv"]));
    assert!(
        typo.contains("a similar argument exists: '--hgvs'"),
        "{typo}"
    );
    let other = stderr(&annotate(&["--no-such-flag-anywhere"]));
    assert!(
        other.starts_with("error: unexpected argument '--no-such-flag-anywhere' found"),
        "{other}"
    );
    assert!(!other.contains("tip: "), "{other}");
}
