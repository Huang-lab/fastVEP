//! Synthetic benchmark for the .osa `SaReader` query path.
//!
//! Builds a ClinVar-scale .osa file, then drives the SA reader through a
//! realistic per-transcript × per-allele query pattern to measure the impact
//! of caching/preload changes. Standalone (no criterion dependency) so it
//! runs the same way before and after the optimisation.
//!
//! Usage:
//!   cargo run --release --example bench_sa_reader -p fastvep-sa
//!
//! Env knobs (optional):
//!   SA_BENCH_RECORDS    Number of SA records to write (default 1_000_000)
//!   SA_BENCH_VARIANTS   Number of user variants to query (default 50_000)
//!   SA_BENCH_TRANSCRIPTS Per-variant transcript fan-out (default 5)
//!   SA_BENCH_BATCH      Preload batch size (default 1024)

use anyhow::Result;
use fastvep_cache::annotation::AnnotationProvider;
use fastvep_sa::common::AnnotationRecord;
use fastvep_sa::index::IndexHeader;
use fastvep_sa::reader::SaReader;
use fastvep_sa::writer::SaWriter;
use std::env;
use std::path::PathBuf;
use std::time::Instant;
use tempfile::TempDir;

const HUMAN_CHROMS: &[(u16, &str, u32)] = &[
    (0, "chr1", 248_956_422),
    (1, "chr2", 242_193_529),
    (2, "chr3", 198_295_559),
    (3, "chr4", 190_214_555),
    (4, "chr5", 181_538_259),
    (5, "chr6", 170_805_979),
    (6, "chr7", 159_345_973),
    (7, "chr8", 145_138_636),
    (8, "chr9", 138_394_717),
    (9, "chr10", 133_797_422),
    (10, "chr11", 135_086_622),
    (11, "chr12", 133_275_309),
    (12, "chr13", 114_364_328),
    (13, "chr14", 107_043_718),
    (14, "chr15", 101_991_189),
    (15, "chr16", 90_338_345),
    (16, "chr17", 83_257_441),
    (17, "chr18", 80_373_285),
    (18, "chr19", 58_617_616),
    (19, "chr20", 64_444_167),
    (20, "chr21", 46_709_983),
    (21, "chr22", 50_818_468),
    (22, "chrX", 156_040_895),
];

fn env_usize(name: &str, default: usize) -> usize {
    env::var(name)
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(default)
}

fn chrom_map() -> Vec<String> {
    HUMAN_CHROMS.iter().map(|(_, n, _)| n.to_string()).collect()
}

fn build_synthetic_osa(path: &PathBuf, n_records: usize) -> Result<u64> {
    // Distribute records across chromosomes weighted by length.
    let total_len: u64 = HUMAN_CHROMS.iter().map(|(_, _, l)| *l as u64).sum();
    let mut records = Vec::with_capacity(n_records);

    let mut allocated = 0usize;
    for (i, (chrom_idx, _name, length)) in HUMAN_CHROMS.iter().enumerate() {
        let share = if i == HUMAN_CHROMS.len() - 1 {
            n_records - allocated
        } else {
            ((*length as u64 * n_records as u64) / total_len) as usize
        };
        allocated += share;

        if share == 0 {
            continue;
        }
        // Evenly-spaced positions on this chromosome, 1-based.
        let step = (*length as u64 / share as u64).max(1);
        for j in 0..share {
            let pos = (1 + j as u64 * step).min(*length as u64) as u32;
            records.push(AnnotationRecord {
                chrom_idx: *chrom_idx,
                position: pos,
                ref_allele: "A".into(),
                alt_allele: "G".into(),
                json: r#"{"sig":"Pathogenic","status":"reviewed_by_expert_panel"}"#.into(),
            });
        }
    }

    records.sort_by(|a, b| {
        a.chrom_idx
            .cmp(&b.chrom_idx)
            .then(a.position.cmp(&b.position))
    });

    let header = IndexHeader {
        schema_version: 1,
        json_key: "clinvar".into(),
        name: "Synthetic ClinVar".into(),
        version: "bench".into(),
        description: "Synthetic benchmark fixture".into(),
        assembly: "GRCh38".into(),
        match_by_allele: true,
        is_array: false,
        is_positional: false,
    };

    let mut writer = SaWriter::new(header);
    writer.write_to_files(path, records.into_iter(), &chrom_map())?;

    let osa_path = path.with_extension("osa");
    let size = std::fs::metadata(&osa_path)?.len();
    Ok(size)
}

fn build_query_workload(n_variants: usize) -> Vec<(String, u64)> {
    // Sorted positions across the genome, mimicking a normal VCF input.
    let total_len: u64 = HUMAN_CHROMS.iter().map(|(_, _, l)| *l as u64).sum();
    let mut variants = Vec::with_capacity(n_variants);
    for (_, name, length) in HUMAN_CHROMS {
        let share = ((*length as u64 * n_variants as u64) / total_len) as usize;
        if share == 0 {
            continue;
        }
        let step = (*length as u64 / share as u64).max(1);
        for j in 0..share {
            let pos = (1 + j as u64 * step).min(*length as u64);
            variants.push((name.to_string(), pos));
        }
    }
    variants
}

fn main() -> Result<()> {
    let n_records = env_usize("SA_BENCH_RECORDS", 1_000_000);
    let n_variants = env_usize("SA_BENCH_VARIANTS", 50_000);
    let transcripts = env_usize("SA_BENCH_TRANSCRIPTS", 5);
    let batch_size = env_usize("SA_BENCH_BATCH", 1024);

    println!(
        "Benchmark config: {} SA records, {} user variants, {}× transcript fan-out, batch={}",
        n_records, n_variants, transcripts, batch_size
    );

    let dir = TempDir::new()?;
    let base = dir.path().join("clinvar_bench");

    println!("Building synthetic .osa fixture...");
    let t0 = Instant::now();
    let osa_size = build_synthetic_osa(&base, n_records)?;
    println!(
        "  built {} MB .osa in {:.2}s",
        osa_size / (1024 * 1024),
        t0.elapsed().as_secs_f64()
    );

    let reader = SaReader::open(&base.with_extension("osa"))?;
    let workload = build_query_workload(n_variants);
    println!("  built workload of {} positions", workload.len());

    // --- Scenario A: preload + per-variant × per-transcript × per-allele queries ---
    println!("\n[scenario] Preload-per-batch + per-(transcript,allele) query pattern");
    let t0 = Instant::now();
    let mut hits = 0u64;
    let mut total_queries = 0u64;
    for chunk in workload.chunks(batch_size) {
        // Group positions by chromosome (sequential preload phase)
        let mut by_chrom: std::collections::HashMap<&str, Vec<u64>> =
            std::collections::HashMap::new();
        for (chrom, pos) in chunk {
            by_chrom.entry(chrom.as_str()).or_default().push(*pos);
        }
        for (chrom, positions) in &by_chrom {
            let _ = reader.preload(chrom, positions);
        }
        // Per-transcript / per-allele queries (parallel in real pipeline; here
        // we just measure total work).
        for (chrom, pos) in chunk {
            for _t in 0..transcripts {
                total_queries += 1;
                if let Ok(Some(_)) = reader.annotate_position(chrom, *pos, "A", "G") {
                    hits += 1;
                }
            }
        }
    }
    let elapsed = t0.elapsed().as_secs_f64();
    println!(
        "  scenario A: {:.2}s ({} queries, {} hits, {:.0} q/s)",
        elapsed,
        total_queries,
        hits,
        total_queries as f64 / elapsed
    );

    Ok(())
}
