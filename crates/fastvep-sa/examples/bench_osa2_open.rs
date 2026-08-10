//! Reproduction/benchmark for issue #78: "after adding gnomad_chr*.osa2 to
//! --sa-dir the pre-startup time increased by 20+ minutes".
//!
//! Builds a set of per-chromosome gnomAD-schema `.osa2` shards (the exact
//! layout `fastvep sa-build --format osa2 --source gnomad` emits) into a
//! temporary "sa-dir", then times `Osa2Reader::open` for every shard the way
//! `load_sa_providers` does at annotate startup.
//!
//! The interesting number is not just wall time but `entries` and
//! `header_reads`: the old open path resolved every ZIP entry's data offset up
//! front, which is one random read per entry scattered across a multi-GB file.
//!
//! Usage:
//!   cargo run --release --example bench_osa2_open -p fastvep-sa
//!
//! Env knobs:
//!   OSA2_BENCH_CHROMS    number of chromosome shards to build (default 4)
//!   OSA2_BENCH_PER_CHROM records per shard                    (default 2_000_000)
//!   OSA2_BENCH_DIR       reuse/keep shards in this directory (default: temp dir)

use anyhow::Result;
use fastvep_sa::reader_v2::Osa2Reader;
use fastvep_sa::sources::gnomad::{gnomad_osa2_fields, gnomad_osa2_metadata};
use fastvep_sa::writer_v2::{Osa2Record, Osa2StreamWriter};
use std::env;
use std::fs::File;
use std::io::BufWriter;
use std::path::{Path, PathBuf};
use std::time::Instant;

const HUMAN_CHROMS: &[(&str, u32)] = &[
    ("chr1", 248_956_422),
    ("chr2", 242_193_529),
    ("chr3", 198_295_559),
    ("chr4", 190_214_555),
    ("chr5", 181_538_259),
    ("chr6", 170_805_979),
    ("chr7", 159_345_973),
    ("chr8", 145_138_636),
    ("chr9", 138_394_717),
    ("chr10", 133_797_422),
    ("chr11", 135_086_622),
    ("chr12", 133_275_309),
    ("chr13", 114_364_328),
    ("chr14", 107_043_718),
    ("chr15", 101_991_189),
    ("chr16", 90_338_345),
    ("chr17", 83_257_441),
    ("chr18", 80_373_285),
    ("chr19", 58_617_616),
    ("chr20", 64_444_167),
    ("chr21", 46_709_983),
    ("chr22", 50_818_468),
    ("chrX", 156_040_895),
];

fn env_usize(name: &str, default: usize) -> usize {
    env::var(name)
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(default)
}

/// Write one per-chromosome gnomAD `.osa2` shard with `n` evenly spread
/// records, using the streaming writer (the path `sa-build` uses).
fn build_shard(path: &Path, chrom: &str, length: u32, n: usize) -> Result<u64> {
    let fields = gnomad_osa2_fields();
    let metadata = gnomad_osa2_metadata("GRCh38");
    let n_fields = fields.len();
    let file = BufWriter::with_capacity(1 << 22, File::create(path)?);
    let mut writer = Osa2StreamWriter::new(file, &metadata, fields.clone())?;

    let step = (length as u64 / n.max(1) as u64).max(1);
    let alts = [b'C', b'G', b'T', b'A'];
    for j in 0..n {
        let pos = (1 + j as u64 * step).min(length as u64) as u32;
        let values: Vec<u32> = (0..n_fields)
            .map(|f| ((pos as u64 + f as u64 * 7) % 1_000_000) as u32)
            .collect();
        writer.push(Osa2Record {
            chrom: chrom.to_string(),
            position: pos,
            ref_allele: vec![b'A'],
            alt_allele: vec![alts[j % alts.len()]],
            values,
            json_blob: None,
        })?;
    }
    writer.finish()?;
    Ok(std::fs::metadata(path)?.len())
}

fn main() -> Result<()> {
    let n_chroms = env_usize("OSA2_BENCH_CHROMS", 4).min(HUMAN_CHROMS.len());
    let per_chrom = env_usize("OSA2_BENCH_PER_CHROM", 2_000_000);

    let (dir, _keep): (PathBuf, Option<tempfile::TempDir>) = match env::var("OSA2_BENCH_DIR") {
        Ok(d) => {
            std::fs::create_dir_all(&d)?;
            (PathBuf::from(d), None)
        }
        Err(_) => {
            let t = tempfile::TempDir::new()?;
            (t.path().to_path_buf(), Some(t))
        }
    };

    println!(
        "sa-dir: {}  ({} shards x {} records)",
        dir.display(),
        n_chroms,
        per_chrom
    );

    let mut paths = Vec::new();
    let mut total_bytes = 0u64;
    for (chrom, length) in HUMAN_CHROMS.iter().take(n_chroms) {
        let path = dir.join(format!("gnomad_{}.osa2", chrom));
        if path.exists() {
            total_bytes += std::fs::metadata(&path)?.len();
            paths.push(path);
            continue;
        }
        let t = Instant::now();
        let bytes = build_shard(&path, chrom, *length, per_chrom)?;
        total_bytes += bytes;
        println!(
            "  built {:<22} {:>8.1} MB in {:>6.1}s",
            path.file_name().unwrap().to_string_lossy(),
            bytes as f64 / 1e6,
            t.elapsed().as_secs_f64()
        );
        paths.push(path);
    }
    println!("total on-disk: {:.2} GB\n", total_bytes as f64 / 1e9);

    // Startup path: open every shard, exactly as load_sa_providers does.
    let t = Instant::now();
    let mut readers = Vec::new();
    for path in &paths {
        let t1 = Instant::now();
        let r = Osa2Reader::open(path)?;
        println!(
            "  open {:<22} entries={:>7}  header_reads={:>7}  {:>7.3}s",
            path.file_name().unwrap().to_string_lossy(),
            r.entry_count(),
            r.header_read_count(),
            t1.elapsed().as_secs_f64()
        );
        readers.push(r);
    }
    let open_secs = t.elapsed().as_secs_f64();
    let entries: usize = readers.iter().map(|r| r.entry_count()).sum();
    let header_reads: u64 = readers.iter().map(|r| r.header_read_count()).sum();
    println!(
        "\nOPEN TOTAL: {:.3}s  entries={}  header_reads={}",
        open_secs, entries, header_reads
    );

    // Sanity: a lookup must still work after open (guards the lazy-offset path).
    use fastvep_cache::annotation::AnnotationProvider;
    let mut hits = 0;
    for (i, r) in readers.iter().enumerate() {
        let (chrom, length) = HUMAN_CHROMS[i];
        let step = (length as u64 / per_chrom.max(1) as u64).max(1);
        let alts = [b'C', b'G', b'T', b'A'];
        for j in [0usize, per_chrom / 2, per_chrom.saturating_sub(1)] {
            let pos = (1 + j as u64 * step).min(length as u64);
            let alt = (alts[j % alts.len()] as char).to_string();
            if r.annotate_position(chrom, pos, "A", &alt)?.is_some() {
                hits += 1;
            }
        }
    }
    println!(
        "post-open lookups: {}/{} hit  header_reads_after={}",
        hits,
        readers.len() * 3,
        readers.iter().map(|r| r.header_read_count()).sum::<u64>()
    );

    Ok(())
}
