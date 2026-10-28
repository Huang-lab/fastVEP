//! The known variants that `--cache-dir` attaches do not depend on how many threads looked them
//! up.
//!
//! The lookups of a batch are made in runs of consecutive records, each on whichever thread takes
//! it, and each run reads the cache's files from a place of its own. A run that started from the
//! wrong place, a window that one run left behind for another, or a run written out of its order
//! would show as a file that differs between one thread and several, and nothing else would say
//! so. The program is run as a user runs it, on one thread and on several (a number of threads
//! that does not divide the 446 records of the world, so that the runs end in different places),
//! with the flags that read the cache.

use std::path::{Path, PathBuf};
use std::process::Command;

use tempfile::TempDir;

fn world() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../fastvep-cache/fixtures/vep_known_columns")
}

/// The file the program writes with `args` on `threads` threads.
fn run(dir: &Path, args: &[&str], threads: u32) -> String {
    let out = dir.join(format!("out_{threads}.tab"));
    let status = Command::new(env!("CARGO_BIN_EXE_fastvep"))
        .env("RAYON_NUM_THREADS", threads.to_string())
        .args(["annotate", "-i"])
        .arg(world().join("in.vcf"))
        .arg("--gff3")
        .arg(dir.join("genes.gff3"))
        .arg("--cache-dir")
        .arg(world().join("cache/homo_sapiens/116_GRCh38"))
        .args(["--output-format", "tab", "--no-progress", "-o"])
        .arg(&out)
        .args(args)
        .output()
        .expect("the program runs");
    assert!(
        status.status.success(),
        "{}",
        String::from_utf8_lossy(&status.stderr)
    );
    std::fs::read_to_string(&out).unwrap()
}

fn same_on_every_number_of_threads(args: &[&str]) {
    let dir = TempDir::new().unwrap();
    // copied, so that the transcript cache a run writes is not in the source tree
    std::fs::copy(world().join("genes.gff3"), dir.path().join("genes.gff3")).unwrap();
    let one = run(dir.path(), args, 1);

    // that the lookups found something, or two empty files would be the same
    let header = one
        .lines()
        .find(|l| l.starts_with("#Uploaded_variation"))
        .expect("a header");
    let existing = header
        .split('\t')
        .position(|c| c == "Existing_variation")
        .expect("an Existing_variation column");
    let known = one
        .lines()
        .filter(|l| !l.starts_with('#'))
        .filter(|l| l.split('\t').nth(existing).is_some_and(|v| v != "-"))
        .count();
    assert!(known > 100, "only {known} rows have a known variant");

    for threads in [4, 7] {
        let many = run(dir.path(), args, threads);
        assert!(
            one == many,
            "{args:?}: the file of one thread is not the file of {threads}"
        );
    }
}

#[test]
fn the_reporters_flags_give_the_same_file_on_one_thread_and_on_several() {
    same_on_every_number_of_threads(&["--symbol", "--canonical", "--af-1kg", "--check-existing"]);
}

#[test]
fn every_flag_that_reads_the_cache_gives_the_same_file_on_one_thread_and_on_several() {
    same_on_every_number_of_threads(&[
        "--symbol",
        "--af",
        "--af-1kg",
        "--af-gnomade",
        "--af-gnomadg",
        "--max-af",
        "--pubmed",
        "--check-existing",
    ]);
}
