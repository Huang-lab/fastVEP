//! What `fastvep annotate` does when its output goes somewhere that is not a file, or that fails.
//!
//! The batches are written on a thread of their own, while the next is read and annotated, so a
//! write that fails is seen by one thread and the run is stopped by another. The unit tests of
//! that hand-over use an output of their own; this runs the program as a user runs it, with its
//! output on a pipe that is read to the end, on a pipe whose reader goes away, and on a device that
//! is always full.

use std::io::Read;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::time::{Duration, Instant};

use tempfile::TempDir;

fn fixture(world: &str, name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("fixtures")
        .join(world)
        .join(name)
}

/// A world in `dir`, copied so that the transcript cache a run writes is not in the source tree.
fn world() -> TempDir {
    let dir = TempDir::new().unwrap();
    for name in ["genes.gff3", "genome.fa"] {
        std::fs::copy(fixture("vep_tab_world", name), dir.path().join(name)).unwrap();
    }
    dir
}

/// `fastvep annotate` over the world of 2,271 changes (two batches and a part), writing `format` to
/// `output`.
fn annotate(dir: &Path, format: &str, output: &str) -> Command {
    let mut command = Command::new(env!("CARGO_BIN_EXE_fastvep"));
    command
        .args(["annotate", "-i"])
        .arg(fixture("vep_tab_edges", "in.vcf"))
        .arg("--gff3")
        .arg(dir.join("genes.gff3"))
        .arg("--fasta")
        .arg(dir.join("genome.fa"))
        .args(["--output-format", format, "-o", output])
        .args(["--hgvs", "--no-progress"]);
    command
}

#[test]
fn what_goes_to_stdout_is_what_goes_to_a_file() {
    let dir = world();
    for format in ["vcf", "tab", "json"] {
        let file = dir.path().join(format!("out.{format}"));
        let written = annotate(dir.path(), format, file.to_str().unwrap())
            .output()
            .unwrap();
        assert!(
            written.status.success(),
            "{}",
            String::from_utf8_lossy(&written.stderr)
        );
        let piped = annotate(dir.path(), format, "-").output().unwrap();
        assert!(
            piped.status.success(),
            "{}",
            String::from_utf8_lossy(&piped.stderr)
        );
        let in_the_file = std::fs::read(&file).unwrap();
        assert!(
            in_the_file.len() > 100_000,
            "{format}: {} bytes",
            in_the_file.len()
        );
        assert!(
            piped.stdout == in_the_file,
            "{format}: stdout and the file differ"
        );
    }
}

/// Waits for `child`, which has to end in `seconds`: a run that is stuck fails the test and is
/// killed, and does not stop the suite.
fn wait_for(child: &mut std::process::Child, seconds: u64) -> std::process::ExitStatus {
    let deadline = Instant::now() + Duration::from_secs(seconds);
    loop {
        if let Some(status) = child.try_wait().unwrap() {
            return status;
        }
        if Instant::now() > deadline {
            let _ = child.kill();
            panic!("the run was still going after {seconds} seconds");
        }
        std::thread::sleep(Duration::from_millis(10));
    }
}

/// What the run writes to stderr, read as it is written so that it never waits for a reader.
fn read_stderr(child: &mut std::process::Child) -> std::thread::JoinHandle<String> {
    let mut stderr = child.stderr.take().unwrap();
    std::thread::spawn(move || {
        let mut text = String::new();
        stderr.read_to_string(&mut text).unwrap();
        text
    })
}

#[test]
fn a_reader_that_goes_away_stops_the_run_with_an_error() {
    // `fastvep annotate ... -o - | head -c 64`. The output is megabytes, which is many times what
    // a pipe holds, so the run is still writing when the reader leaves; and it has to say that it
    // stopped, and not exit as though it had written what it was asked to.
    let dir = world();
    let mut child = annotate(dir.path(), "vcf", "-")
        .args(["--distance", "100000"])
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    let stderr = read_stderr(&mut child);
    let mut stdout = child.stdout.take().unwrap();
    let mut first = [0u8; 64];
    stdout.read_exact(&mut first).unwrap();
    assert!(
        first.starts_with(b"##fileformat=VCF"),
        "{:?}",
        String::from_utf8_lossy(&first)
    );
    drop(stdout);

    let status = wait_for(&mut child, 120);
    let stderr = stderr.join().unwrap();
    assert!(
        !status.success(),
        "the run succeeded, with nobody to read what it wrote"
    );
    assert!(stderr.contains("Broken pipe"), "{stderr}");
}

#[test]
#[cfg(target_os = "linux")]
fn a_full_disk_stops_the_run_with_an_error() {
    // `/dev/full` takes a write and answers that the disk is full
    let dir = world();
    let run = annotate(dir.path(), "vcf", "/dev/full").output().unwrap();
    let stderr = String::from_utf8_lossy(&run.stderr);
    assert!(!run.status.success(), "the run succeeded: {stderr}");
    assert!(stderr.contains("No space left on device"), "{stderr}");
}
