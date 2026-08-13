//! Head-to-head benchmark: v1 `.osa` (zstd/JSON) vs v2 `.osa2` (chunked/u32).
//!
//! Builds the SAME synthetic gnomAD-scale dataset in both formats, then drives
//! an identical per-transcript × per-allele query workload through each reader.
//! Reports build time, on-disk size, and query throughput so the "is v2 better?"
//! question is answered with measurement rather than assertion.
//!
//! Usage:
//!   cargo run --release --example bench_v1_vs_v2 -p fastvep-sa
//!
//! Env knobs (optional):
//!   SA_BENCH_RECORDS     SA records to write         (default 2_000_000)
//!   SA_BENCH_VARIANTS    user variants to query      (default 100_000)
//!   SA_BENCH_TRANSCRIPTS per-variant transcript fan-out (default 5)
//!   SA_BENCH_BATCH       preload batch size          (default 1024)
//!   SA_BENCH_HIT_RATE    fraction of queries that hit an annotated site (default 0.5)

use anyhow::Result;
use fastvep_cache::annotation::AnnotationProvider;
use fastvep_sa::common::{AnnotationRecord, SCHEMA_VERSION};
use fastvep_sa::fields::{Field, FieldType};
use fastvep_sa::index::IndexHeader;
use fastvep_sa::reader::SaReader;
use fastvep_sa::reader_v2::Osa2Reader;
use fastvep_sa::writer::SaWriter;
use fastvep_sa::writer_v2::{Osa2Metadata, Osa2Record, Osa2Writer};
use std::env;
use std::path::{Path, PathBuf};
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

fn env_f64(name: &str, default: f64) -> f64 {
    env::var(name)
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(default)
}

fn chrom_map() -> Vec<String> {
    HUMAN_CHROMS.iter().map(|(_, n, _)| n.to_string()).collect()
}

/// One synthetic annotated site: position plus a deterministic gnomAD-style
/// payload (global AF, AN, AC). Shared by both format builders so the two
/// files hold byte-for-byte equivalent information.
struct SiteData {
    chrom_idx: u16,
    chrom: &'static str,
    position: u32,
    af: f64,
    an: i64,
    ac: i64,
}

fn build_sites(n_records: usize) -> Vec<SiteData> {
    let total_len: u64 = HUMAN_CHROMS.iter().map(|(_, _, l)| *l as u64).sum();
    let mut sites = Vec::with_capacity(n_records);
    let mut allocated = 0usize;
    for (i, (chrom_idx, name, length)) in HUMAN_CHROMS.iter().enumerate() {
        let share = if i == HUMAN_CHROMS.len() - 1 {
            n_records - allocated
        } else {
            ((*length as u64 * n_records as u64) / total_len) as usize
        };
        allocated += share;
        if share == 0 {
            continue;
        }
        let step = (*length as u64 / share as u64).max(1);
        for j in 0..share {
            let pos = (1 + j as u64 * step).min(*length as u64) as u32;
            // Deterministic pseudo-values so both formats encode identical data
            // without needing a RNG (Math.random-free, reproducible).
            let af = ((pos as u64 % 100_000) as f64) / 1_000_000.0; // 0.0 .. ~0.1
            let an = 150_000i64;
            let ac = (af * an as f64) as i64;
            sites.push(SiteData {
                chrom_idx: *chrom_idx,
                chrom: name,
                position: pos,
                af,
                an,
                ac,
            });
        }
    }
    sites
}

fn build_v1(path: &Path, sites: &[SiteData]) -> Result<u64> {
    let mut records: Vec<AnnotationRecord> = sites
        .iter()
        .map(|s| AnnotationRecord {
            chrom_idx: s.chrom_idx,
            position: s.position,
            ref_allele: "A".into(),
            alt_allele: "G".into(),
            json: format!(
                r#"{{"allAf":{:.6e},"allAn":{},"allAc":{}}}"#,
                s.af, s.an, s.ac
            ),
        })
        .collect();
    records.sort_by(|a, b| {
        a.chrom_idx
            .cmp(&b.chrom_idx)
            .then(a.position.cmp(&b.position))
    });

    let header = IndexHeader {
        schema_version: SCHEMA_VERSION,
        json_key: "gnomad".into(),
        name: "Synthetic gnomAD (v1)".into(),
        version: "bench".into(),
        description: "Synthetic benchmark fixture".into(),
        assembly: "GRCh38".into(),
        match_by_allele: true,
        is_array: false,
        is_positional: false,
    };
    let mut writer = SaWriter::new(header);
    writer.write_to_files(path, records.into_iter(), &chrom_map())?;
    Ok(std::fs::metadata(path.with_extension("osa"))?.len())
}

fn v2_fields() -> Vec<Field> {
    vec![
        Field {
            field: "AF".into(),
            alias: "allAf".into(),
            ftype: FieldType::Float,
            multiplier: 1_000_000,
            zigzag: false,
            missing_value: u32::MAX,
            missing_string: ".".into(),
            description: "Global allele frequency".into(),
        },
        Field {
            field: "AN".into(),
            alias: "allAn".into(),
            ftype: FieldType::Integer,
            multiplier: 1,
            zigzag: false,
            missing_value: u32::MAX,
            missing_string: ".".into(),
            description: "Total allele number".into(),
        },
        Field {
            field: "AC".into(),
            alias: "allAc".into(),
            ftype: FieldType::Integer,
            multiplier: 1,
            zigzag: false,
            missing_value: u32::MAX,
            missing_string: ".".into(),
            description: "Allele count".into(),
        },
    ]
}

fn build_v2(path: &Path, sites: &[SiteData]) -> Result<u64> {
    let fields = v2_fields();
    let mut records: Vec<Osa2Record> = sites
        .iter()
        .map(|s| Osa2Record {
            chrom: s.chrom.to_string(),
            position: s.position,
            ref_allele: b"A".to_vec(),
            alt_allele: b"G".to_vec(),
            values: vec![
                fields[0].encode_float(s.af),
                fields[1].encode_int(s.an),
                fields[2].encode_int(s.ac),
            ],
            json_blob: None,
        })
        .collect();
    records.sort_by(|a, b| a.chrom.cmp(&b.chrom).then(a.position.cmp(&b.position)));

    let metadata = Osa2Metadata {
        format_version: 2,
        name: "Synthetic gnomAD (v2)".into(),
        version: "bench".into(),
        assembly: "GRCh38".into(),
        json_key: "gnomad".into(),
        match_by_allele: true,
        is_array: false,
        is_positional: false,
        chunk_bits: 20,
        description: "Synthetic benchmark fixture".into(),
    };
    let writer = Osa2Writer::new(metadata, fields);
    let file = std::fs::File::create(path)?;
    writer.write_all(std::io::BufWriter::new(file), &records)?;
    Ok(std::fs::metadata(path)?.len())
}

/// Sorted query positions across the genome, mimicking a normal VCF input.
/// `hit_rate` controls how many land on an annotated site vs. a definite miss
/// (offset by +1 bp, which no synthetic record occupies given the spacing).
fn build_workload(
    sites: &[SiteData],
    n_variants: usize,
    hit_rate: f64,
) -> Vec<(&'static str, u64)> {
    if sites.is_empty() {
        return Vec::new();
    }
    let mut out = Vec::with_capacity(n_variants);
    // Stride through the annotated sites so queries stay genome-sorted.
    let stride = (sites.len() / n_variants.max(1)).max(1);
    let mut i = 0usize;
    let mut n = 0usize;
    while n < n_variants && i < sites.len() {
        let s = &sites[i];
        // Deterministic hit/miss split: first `hit_rate` share of each 100
        // queries hit, the rest are shifted to a guaranteed-empty position.
        let is_hit = ((n % 100) as f64) < hit_rate * 100.0;
        let pos = if is_hit {
            s.position as u64
        } else {
            s.position as u64 + 1
        };
        out.push((s.chrom, pos));
        i += stride;
        n += 1;
    }
    out
}

/// Run the preload-per-batch + per-(transcript,allele) query pattern.
fn run_workload(
    reader: &dyn AnnotationProvider,
    workload: &[(&'static str, u64)],
    transcripts: usize,
    batch_size: usize,
) -> Result<(f64, u64, u64)> {
    let t0 = Instant::now();
    let mut hits = 0u64;
    let mut total = 0u64;
    for chunk in workload.chunks(batch_size) {
        let mut by_chrom: std::collections::HashMap<&str, Vec<u64>> =
            std::collections::HashMap::new();
        for (chrom, pos) in chunk {
            by_chrom.entry(chrom).or_default().push(*pos);
        }
        for (chrom, positions) in &by_chrom {
            reader.preload(chrom, positions)?;
        }
        for (chrom, pos) in chunk {
            for _t in 0..transcripts {
                total += 1;
                if let Ok(Some(_)) = reader.annotate_position(chrom, *pos, "A", "G") {
                    hits += 1;
                }
            }
        }
    }
    Ok((t0.elapsed().as_secs_f64(), total, hits))
}

fn mb(bytes: u64) -> f64 {
    bytes as f64 / (1024.0 * 1024.0)
}

fn main() -> Result<()> {
    let n_records = env_usize("SA_BENCH_RECORDS", 2_000_000);
    let n_variants = env_usize("SA_BENCH_VARIANTS", 100_000);
    let transcripts = env_usize("SA_BENCH_TRANSCRIPTS", 5);
    let batch_size = env_usize("SA_BENCH_BATCH", 1024);
    let hit_rate = env_f64("SA_BENCH_HIT_RATE", 0.5);

    println!(
        "Config: {} records, {} variants, {}x fan-out, batch={}, hit_rate={:.2}\n",
        n_records, n_variants, transcripts, batch_size, hit_rate
    );

    let dir = TempDir::new()?;
    let base: PathBuf = dir.path().join("bench");

    println!("Generating {} synthetic sites...", n_records);
    let sites = build_sites(n_records);
    let workload = build_workload(&sites, n_variants, hit_rate);
    println!(
        "  {} sites, {} query positions\n",
        sites.len(),
        workload.len()
    );

    // --- Build ---
    let v1_base = base.with_extension("v1");
    let t0 = Instant::now();
    let v1_size = build_v1(&v1_base, &sites)?;
    let v1_build = t0.elapsed().as_secs_f64();

    let v2_path = base.with_extension("osa2");
    let t0 = Instant::now();
    let v2_size = build_v2(&v2_path, &sites)?;
    let v2_build = t0.elapsed().as_secs_f64();

    println!("BUILD");
    println!("  v1 (.osa) : {:.1} MB  in {:.2}s", mb(v1_size), v1_build);
    println!("  v2 (.osa2): {:.1} MB  in {:.2}s", mb(v2_size), v2_build);
    println!(
        "  -> v2 is {:.2}x the size, {:.2}x the build time\n",
        v2_size as f64 / v1_size as f64,
        v2_build / v1_build
    );

    // --- Query ---
    let v1_reader = SaReader::open(&v1_base.with_extension("osa"))?;
    let v2_reader = Osa2Reader::open(&v2_path)?;

    let (v1_t, v1_total, v1_hits) = run_workload(&v1_reader, &workload, transcripts, batch_size)?;
    let (v2_t, v2_total, v2_hits) = run_workload(&v2_reader, &workload, transcripts, batch_size)?;

    println!("QUERY (preload + per-(transcript,allele))");
    println!(
        "  v1 (.osa) : {:.2}s  {:.0} q/s  ({} hits / {} queries)",
        v1_t,
        v1_total as f64 / v1_t,
        v1_hits,
        v1_total
    );
    println!(
        "  v2 (.osa2): {:.2}s  {:.0} q/s  ({} hits / {} queries)",
        v2_t,
        v2_total as f64 / v2_t,
        v2_hits,
        v2_total
    );
    println!("  -> v2 query throughput is {:.2}x v1\n", v1_t / v2_t);

    if v1_hits != v2_hits {
        println!(
            "WARNING: hit counts differ (v1={}, v2={}) — formats are not returning \
             equivalent results!",
            v1_hits, v2_hits
        );
    } else {
        println!(
            "OK: both formats returned {} hits (equivalent lookups)",
            v1_hits
        );
    }

    Ok(())
}
