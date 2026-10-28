//! What `fastvep annotate` writes does not depend on how many threads it ran on.
//!
//! The per-variant work runs in parallel and the rows are written in input order. A row out of
//! place, or a cell that one thread reads from another's half-built state, would show as a file
//! that differs between runs, and nothing else would say so. The program is run as a user runs
//! it, once on one thread and once on four, over a world of 4,753 changes, for each of the three
//! formats and for the flags that change the shape of a row.

use std::path::{Path, PathBuf};
use std::process::Command;

use tempfile::TempDir;

fn fixture(world: &str, name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("fixtures")
        .join(world)
        .join(name)
}

/// The file the program writes with `args` on `threads` threads.
fn run(dir: &Path, format: &str, args: &[&str], threads: u32) -> Vec<u8> {
    let out = dir.join(format!("out_{format}_{threads}"));
    let status = Command::new(env!("CARGO_BIN_EXE_fastvep"))
        .env("RAYON_NUM_THREADS", threads.to_string())
        .args(["annotate", "-i"])
        .arg(fixture("vep_tab_edges", "in.vcf"))
        .arg("--gff3")
        .arg(dir.join("genes.gff3"))
        .arg("--fasta")
        .arg(dir.join("genome.fa"))
        .args(["--output-format", format, "-o"])
        .arg(&out)
        .args(["--hgvs", "--no-progress"])
        .args(args)
        .output()
        .expect("the program runs");
    assert!(
        status.status.success(),
        "{}",
        String::from_utf8_lossy(&status.stderr)
    );
    std::fs::read(&out).unwrap()
}

fn same_on_one_thread_and_four(format: &str, args: &[&str]) {
    let dir = TempDir::new().unwrap();
    // copied, so that the transcript cache a run writes is not in the source tree and the runs do
    // not share one
    std::fs::copy(
        fixture("vep_tab_world", "genes.gff3"),
        dir.path().join("genes.gff3"),
    )
    .unwrap();
    std::fs::copy(
        fixture("vep_tab_world", "genome.fa"),
        dir.path().join("genome.fa"),
    )
    .unwrap();
    let one = run(dir.path(), format, args, 1);
    let four = run(dir.path(), format, args, 4);
    assert!(one.len() > 100_000, "{format}: only {} bytes", one.len());
    assert!(
        one == four,
        "{format} {args:?}: the file of one thread is not the file of four"
    );
}

#[test]
fn the_tab_file_is_the_same_on_one_thread_and_on_four() {
    same_on_one_thread_and_four(
        "tab",
        &[
            "--symbol",
            "--biotype",
            "--canonical",
            "--numbers",
            "--variant-class",
            "--show-ref-allele",
            "--hgvsg",
        ],
    );
}

#[test]
fn the_vcf_is_the_same_on_one_thread_and_on_four() {
    same_on_one_thread_and_four("vcf", &[]);
}

#[test]
fn the_json_is_the_same_on_one_thread_and_on_four() {
    same_on_one_thread_and_four("json", &[]);
}

#[test]
fn a_pick_is_the_same_on_one_thread_and_on_four() {
    same_on_one_thread_and_four("tab", &["--pick", "--symbol", "--canonical"]);
    same_on_one_thread_and_four("tab", &["--flag-pick-allele", "--symbol"]);
}
