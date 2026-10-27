//! The scripts of `validation/` that put a command line together, run under the shell of the machine
//! they are for, with `docker`, `python3` and `fastvep` replaced by programs that write down what
//! they were given.
//!
//! `run_cache_oracle.sh` has arrays that are empty for one of its two kinds (`vep_extra` for
//! `ensembl`, `drop` for it too), and macOS's `/bin/bash` is 3.2, where expanding an empty array
//! under `set -u` is `vep_extra[@]: unbound variable`: the script stopped before it ran anything,
//! for `ensembl` and only there. Bash 4.4 and later expand it to nothing, so a machine with one
//! cannot show it, and the test holds the arguments each kind passes on whichever shell it has.

#![cfg(unix)]

use std::os::unix::fs::PermissionsExt;
use std::path::{Path, PathBuf};
use std::process::Command;

use tempfile::TempDir;

fn script() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../validation/run_cache_oracle.sh")
}

/// A program that writes its arguments, one to a line, to `<dir>/<name>.args`.
fn recorder(dir: &Path, name: &str) -> PathBuf {
    let path = dir.join(name);
    std::fs::write(
        &path,
        format!(
            "#!/bin/sh\nprintf '%s\\n' \"$@\" > '{}/{name}.args'\n",
            dir.display()
        ),
    )
    .unwrap();
    std::fs::set_permissions(&path, std::fs::Permissions::from_mode(0o755)).unwrap();
    path
}

/// What the script passed to each program, for a run of `kind`.
fn run(kind: &str) -> (Vec<String>, Vec<String>, Vec<String>) {
    let dir = TempDir::new().unwrap();
    let bin = dir.path().join("bin");
    std::fs::create_dir(&bin).unwrap();
    recorder(&bin, "docker");
    recorder(&bin, "python3");
    let fastvep = recorder(&bin, "fastvep");
    let (vcf, gff3, fasta) = (
        dir.path().join("in.vcf"),
        dir.path().join("genes.gff3"),
        dir.path().join("genome.fa"),
    );
    for f in [&vcf, &gff3, &fasta] {
        std::fs::write(f, "").unwrap();
    }
    let cache = dir.path().join("cache");
    std::fs::create_dir(&cache).unwrap();
    let path = format!("{}:{}", bin.display(), std::env::var("PATH").unwrap());
    let output = Command::new("/bin/bash")
        .arg(script())
        .args([
            kind,
            &vcf.to_string_lossy(),
            &dir.path().join("out").to_string_lossy(),
        ])
        .env("PATH", path)
        .env("VEP_CACHE", &cache)
        .env("GFF3", &gff3)
        .env("FASTA", &fasta)
        .env("FASTVEP", &fastvep)
        .output()
        .expect("bash runs");
    assert!(
        output.status.success(),
        "{kind}: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let args = |name: &str| -> Vec<String> {
        std::fs::read_to_string(bin.join(format!("{name}.args")))
            .unwrap_or_else(|e| panic!("{kind}: {name} was not run: {e}"))
            .lines()
            .map(String::from)
            .collect()
    };
    (args("docker"), args("fastvep"), args("python3"))
}

#[test]
fn the_ensembl_kind_runs_with_no_extra_vep_flag_and_no_dropped_transcripts() {
    let (docker, fastvep, python) = run("ensembl");
    assert!(docker.contains(&"--cache".to_string()), "{docker:?}");
    assert!(!docker.iter().any(|a| a == "--refseq"), "{docker:?}");
    for flag in ["--exclude-readthrough", "--exclude-artifact"] {
        assert!(fastvep.contains(&flag.to_string()), "{fastvep:?}");
    }
    // the script, the two tables and the prefix of the comparison, and nothing after
    assert_eq!(python.len(), 4, "{python:?}");
    assert!(python[0].ends_with("compare_cache_run.py"), "{python:?}");
}

#[test]
fn the_refseq_kind_leaves_out_predicted_transcripts_and_the_ones_vep_edits() {
    let (docker, fastvep, python) = run("refseq");
    for flag in ["--refseq", "--exclude_predicted"] {
        assert!(docker.contains(&flag.to_string()), "{docker:?}");
    }
    assert!(
        fastvep.contains(&"--exclude-predicted".to_string()),
        "{fastvep:?}"
    );
    assert_eq!(
        &python[4..],
        ["--drop-tx-where", "BAM_EDIT=OK"],
        "{python:?}"
    );
}
