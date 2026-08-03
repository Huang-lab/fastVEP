//! Reproduces issue #75: adding a large, dense .osa (e.g. spliceai, 39 GB)
//! blows up annotate wall-time. Mimics the real pipeline access pattern:
//! coordinate-sorted variants processed batch-by-batch, each batch preloaded
//! sequentially then queried in *parallel* (rayon `par_iter`), exactly as
//! `run_annotate` does.
//!
//! The key diagnostic is `SaReader::decompress_count()`: on sorted input each
//! block should be decompressed ~once. A count far above the number of
//! distinct blocks means the per-reader block cache is too small for the
//! parallel working set (one in-flight block per worker thread), so threads
//! evict each other's blocks and re-decompress the same 8 MiB blocks over and
//! over — this is the spliceai slowdown.
//!
//! It reports OLD (fixed 32 MiB) vs NEW (thread-scaled) budgets side by side,
//! using `SaReader::open_with_cache_budget` so both run in one process.
//!
//! Usage (synthetic):
//!   cargo run --release --example bench_sa_parallel -p fastvep-sa
//! Usage (real data — every query hits a production .osa):
//!   SA_BENCH_OSA=spliceai_chr1.osa SA_BENCH_TSV=<dump_sa output> \
//!     cargo run --release --example bench_sa_parallel -p fastvep-sa
//!
//! Env knobs:
//!   SA_BENCH_OSA       real .osa to query (enables real-data mode)
//!   SA_BENCH_TSV       `dump_sa` TSV of chrom/pos/ref/alt to query
//!   SA_BENCH_SAMPLE    real-data: keep 1-in-N dump rows (default 4)
//!   SA_BENCH_POSITIONS synthetic: distinct positions ×3 alts (default 3_000_000)
//!   SA_BENCH_VARIANTS  synthetic: variants to query (default 100_000)
//!   SA_BENCH_STRIDE    synthetic: bp between consecutive variants (default 300)
//!   SA_BENCH_BATCH     preload batch size (default 1024)

use anyhow::Result;
use fastvep_cache::annotation::AnnotationProvider;
use fastvep_sa::common::{AnnotationRecord, SCHEMA_VERSION};
use fastvep_sa::index::IndexHeader;
use fastvep_sa::reader::SaReader;
use fastvep_sa::writer::SaWriter;
use rayon::prelude::*;
use std::env;
use std::fs::File;
use std::path::PathBuf;
use std::time::Instant;
use tempfile::TempDir;

fn env_usize(name: &str, default: usize) -> usize {
    env::var(name)
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(default)
}

/// Build a dense, spliceai-like .osa on a single contig: three alt records at
/// every consecutive position, with a fat JSON payload, so 8 MiB blocks each
/// span only ~13 kbp of genome — exactly like real precomputed SpliceAI, where
/// every genic base is annotated. `n_positions` distinct positions → `3 *
/// n_positions` records.
fn build_synthetic_osa(path: &std::path::Path, n_positions: usize) -> Result<(u64, u32)> {
    let mut records = Vec::with_capacity(n_positions * 3);
    let alts = ["G", "C", "T"];
    for p in 0..n_positions {
        let pos = (p + 1) as u32;
        for alt in alts {
            records.push(AnnotationRecord {
                chrom_idx: 0,
                position: pos,
                ref_allele: "A".into(),
                alt_allele: alt.into(),
                json: r#"{"gene":"GENEX","dsAg":0.01,"dsAl":0.0,"dsDg":0.85,"dsDl":0.0,"dpAg":5,"dpAl":-28,"dpDg":2,"dpDl":-13}"#.into(),
            });
        }
    }

    let header = IndexHeader {
        schema_version: SCHEMA_VERSION,
        json_key: "spliceai".into(),
        name: "Synthetic SpliceAI".into(),
        version: "bench".into(),
        description: "Synthetic dense fixture".into(),
        assembly: "GRCh38".into(),
        match_by_allele: true,
        is_array: false,
        is_positional: false,
    };

    let mut writer = SaWriter::new(header);
    writer.write_to_files(path, records.into_iter(), &["chr1".to_string()])?;
    let osa_path = path.with_extension("osa");
    Ok((std::fs::metadata(&osa_path)?.len(), n_positions as u32))
}

/// Sorted query positions on chr1, spaced `stride` bp apart — mimics a WGS VCF
/// whose consecutive variants are a few hundred bp apart, so a 1024-variant
/// batch spans many dense blocks.
fn build_query_workload(n_variants: usize, max_pos: u32, stride: u64) -> Vec<(String, u64)> {
    let mut variants = Vec::with_capacity(n_variants);
    let mut pos = 1u64;
    for _ in 0..n_variants {
        if pos > max_pos as u64 {
            break;
        }
        variants.push(("chr1".to_string(), pos));
        pos += stride;
    }
    variants
}

/// A query with an explicit ref/alt (for match_by_allele sources like SpliceAI).
type Query = (String, u64, String, String);

/// Real-data mode: load a `<file>.osa` and a TSV of `chrom\tpos\tref\talt[...]`
/// (e.g. produced by `dump_sa`), so every query is guaranteed to hit and the
/// benchmark exercises the real block layout of a production SpliceAI file.
fn load_real_workload(tsv: &str, sample: usize) -> Result<Vec<Query>> {
    use std::io::BufRead;
    let f = std::io::BufReader::new(File::open(tsv)?);
    let mut out = Vec::new();
    for (i, line) in f.lines().enumerate() {
        let line = line?;
        if i == 0 && line.starts_with("chrom") {
            continue; // header
        }
        if i % sample != 0 {
            continue;
        }
        let mut it = line.split('\t');
        let (Some(c), Some(p), Some(r), Some(a)) = (it.next(), it.next(), it.next(), it.next())
        else {
            continue;
        };
        if c == "?" {
            continue;
        }
        if let Ok(pos) = p.parse::<u64>() {
            out.push((c.to_string(), pos, r.to_string(), a.to_string()));
        }
    }
    out.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
    Ok(out)
}

fn main() -> Result<()> {
    let batch_size = env_usize("SA_BENCH_BATCH", 1024);
    let threads = rayon::current_num_threads();

    // Optional real-data mode: SA_BENCH_OSA=<file.osa> SA_BENCH_TSV=<dump.tsv>.
    let real_osa = env::var("SA_BENCH_OSA").ok();
    let (osa, workload, _dir_guard): (PathBuf, Vec<Query>, Option<TempDir>) =
        if let Some(osa_path) = real_osa {
            let tsv = env::var("SA_BENCH_TSV").expect("SA_BENCH_TSV required with SA_BENCH_OSA");
            let sample = env_usize("SA_BENCH_SAMPLE", 4);
            println!(
                "REAL data: osa={} tsv={} (1-in-{} sample)",
                osa_path, tsv, sample
            );
            let wl = load_real_workload(&tsv, sample)?;
            (PathBuf::from(osa_path), wl, None)
        } else {
            // ~3 M positions × 3 alts ≈ 9 M dense records → hundreds of 8 MiB blocks.
            let n_positions = env_usize("SA_BENCH_POSITIONS", 3_000_000);
            let n_variants = env_usize("SA_BENCH_VARIANTS", 100_000);
            let stride = env_usize("SA_BENCH_STRIDE", 300) as u64;
            println!(
                "SYNTHETIC: {} positions ×3 alts, {} variants @ {}bp stride",
                n_positions, n_variants, stride
            );
            let dir = TempDir::new()?;
            let base = dir.path().join("spliceai_bench");
            let t0 = Instant::now();
            let (osa_size, max_pos) = build_synthetic_osa(&base, n_positions)?;
            println!(
                "  built {} MB in {:.1}s",
                osa_size / (1024 * 1024),
                t0.elapsed().as_secs_f64()
            );
            let wl = build_query_workload(n_variants, max_pos, stride)
                .into_iter()
                .map(|(c, p)| (c, p, "A".to_string(), "G".to_string()))
                .collect();
            (base.with_extension("osa"), wl, Some(dir))
        };

    let n_variants = workload.len();
    println!(
        "config: {} queries, batch={}, rayon threads={}",
        n_variants, batch_size, threads
    );
    println!("  workload of {} sorted positions", workload.len());

    // Mimic run_annotate: sequential preload per batch, then PARALLEL query.
    let run = |reader: &SaReader| -> (f64, u64, u64) {
        let t0 = Instant::now();
        let mut hits = 0u64;
        for chunk in workload.chunks(batch_size) {
            let mut by_chrom: std::collections::HashMap<&str, Vec<u64>> =
                std::collections::HashMap::new();
            for (chrom, pos, _, _) in chunk {
                by_chrom.entry(chrom.as_str()).or_default().push(*pos);
            }
            for (chrom, positions) in &by_chrom {
                reader.preload(chrom, positions).unwrap();
            }
            hits += chunk
                .par_iter()
                .map(|(chrom, pos, r, a)| {
                    matches!(reader.annotate_position(chrom, *pos, r, a), Ok(Some(_))) as u64
                })
                .sum::<u64>();
        }
        (t0.elapsed().as_secs_f64(), hits, reader.decompress_count())
    };

    const MIB: usize = 1024 * 1024;
    // Old behavior (fixed 32 MiB) vs the new thread-scaled shared budget.
    let scaled = (threads * 3 * 8).clamp(32, 384); // blocks×8 MiB, matches reader clamp
    for (label, budget_mib) in [("OLD fixed 32 MiB", 32usize), ("NEW thread-scaled", scaled)] {
        let reader = SaReader::open_with_cache_budget(&osa, budget_mib * MIB)?;
        let (wall, hits, decomps) = run(&reader);
        println!("\n[{}]  budget={} MiB", label, budget_mib);
        println!(
            "  wall: {:.2}s  ({:.0} queries/s)",
            wall,
            n_variants as f64 / wall
        );
        println!("  hits: {}  block decompressions: {}", hits, decomps);
        println!(
            "  → {:.1} decompressions per 1k queries (lower is better)",
            decomps as f64 / (n_variants as f64 / 1000.0)
        );
    }
    Ok(())
}
