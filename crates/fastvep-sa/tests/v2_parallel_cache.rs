//! Issue #75 (v2 `.osa2`): the chunk reader must (a) return correct values
//! through its lock-free mmap read path and (b) not re-decompress chunks under
//! parallel queries when the shared cache is large enough. Mirrors the v1
//! `.osa` block-cache regression test.

use fastvep_cache::annotation::AnnotationProvider;
use fastvep_sa::fields::{Field, FieldType};
use fastvep_sa::reader_v2::Osa2Reader;
use fastvep_sa::writer_v2::{Osa2Metadata, Osa2Record, Osa2Writer};
use rayon::prelude::*;
use std::path::Path;
use tempfile::TempDir;

fn fields() -> Vec<Field> {
    vec![Field {
        field: "AF".into(),
        alias: "af".into(),
        ftype: FieldType::Float,
        multiplier: 1_000_000,
        zigzag: false,
        missing_value: u32::MAX,
        missing_string: ".".into(),
        description: "allele frequency".into(),
    }]
}

/// Dense fixture: one record per position across `n` positions, with a small
/// `chunk_bits` so the file spans many chunks. Returns (path base, n, chunk_bits).
fn build(dir: &Path, n: u32, chunk_bits: u32) -> std::path::PathBuf {
    let flds = fields();
    let records: Vec<Osa2Record> = (0..n)
        .map(|i| Osa2Record {
            chrom: "chr1".into(),
            position: 1000 + i,
            ref_allele: b"A".to_vec(),
            alt_allele: b"G".to_vec(),
            values: vec![flds[0].encode_float((i as f64 % 1000.0) / 1000.0)],
            json_blob: None,
        })
        .collect();
    let metadata = Osa2Metadata {
        format_version: 2,
        name: "dense v2".into(),
        version: "test".into(),
        assembly: "GRCh38".into(),
        json_key: "dense".into(),
        match_by_allele: true,
        is_array: false,
        is_positional: false,
        chunk_bits,
        description: String::new(),
    };
    let path = dir.join("dense.osa2");
    let writer = Osa2Writer::new(metadata, flds);
    let file = std::fs::File::create(&path).unwrap();
    writer
        .write_all(std::io::BufWriter::new(file), &records)
        .unwrap();
    path
}

#[test]
fn lock_free_reads_return_correct_values() {
    let dir = TempDir::new().unwrap();
    // chunk_bits=8 → 256 positions per chunk; 4000 positions → ~16 chunks.
    let path = build(dir.path(), 4000, 8);
    let reader = Osa2Reader::open(&path).unwrap();

    // Every stored position must resolve to its exact value; a +0.5 bp miss
    // (no record there) must return None.
    for i in [0u32, 1, 255, 256, 1234, 3999] {
        let pos = 1000 + i;
        let got = reader
            .annotate_position("chr1", pos as u64, "A", "G")
            .unwrap()
            .expect("stored position must hit");
        let json = match got {
            fastvep_cache::annotation::AnnotationValue::Json(j) => j,
            other => panic!("unexpected value kind: {:?}", other),
        };
        let expected_af = (i as f64 % 1000.0) / 1000.0;
        // Value survived the u32 round-trip (6-digit multiplier).
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let af = v["af"].as_f64().unwrap();
        assert!(
            (af - expected_af).abs() < 1e-6,
            "pos {} af {} != {}",
            pos,
            af,
            expected_af
        );
    }

    // A definite miss (alt not present) returns None, not a bogus hit.
    assert!(reader
        .annotate_position("chr1", 1000, "A", "T")
        .unwrap()
        .is_none());
    // chr/bare aliasing: querying bare "1" resolves to the on-disk "chr1".
    assert!(reader
        .annotate_position("1", 1000, "A", "G")
        .unwrap()
        .is_some());
}

#[test]
fn same_chunk_id_on_different_chromosomes_does_not_collide() {
    // `chr1` and `chr2` both have a record in numeric chunk 0. The chunk cache
    // must key on (chrom, chunk_id), not chunk_id alone — otherwise the first
    // chromosome's chunk 0 is served for the second's. Interleave the queries
    // so the second lookup is answered from a warm cache.
    let dir = TempDir::new().unwrap();
    let flds = fields();
    let records = vec![
        Osa2Record {
            chrom: "chr1".into(),
            position: 100,
            ref_allele: b"A".to_vec(),
            alt_allele: b"G".to_vec(),
            values: vec![flds[0].encode_float(0.111)],
            json_blob: None,
        },
        Osa2Record {
            chrom: "chr2".into(),
            position: 150,
            ref_allele: b"C".to_vec(),
            alt_allele: b"T".to_vec(),
            values: vec![flds[0].encode_float(0.222)],
            json_blob: None,
        },
    ];
    let metadata = Osa2Metadata {
        format_version: 2,
        name: "multi-chrom".into(),
        version: "test".into(),
        assembly: "GRCh38".into(),
        json_key: "mc".into(),
        match_by_allele: true,
        is_array: false,
        is_positional: false,
        chunk_bits: 20, // both positions land in chunk 0 of their chromosome
        description: String::new(),
    };
    let path = dir.path().join("mc.osa2");
    let writer = Osa2Writer::new(metadata, flds);
    writer
        .write_all(
            std::io::BufWriter::new(std::fs::File::create(&path).unwrap()),
            &records,
        )
        .unwrap();

    let reader = Osa2Reader::open(&path).unwrap();
    let af = |chrom, pos, r, a| -> f64 {
        let v = match reader.annotate_position(chrom, pos, r, a).unwrap().unwrap() {
            fastvep_cache::annotation::AnnotationValue::Json(j) => {
                serde_json::from_str::<serde_json::Value>(&j).unwrap()
            }
            _ => unreachable!(),
        };
        v["af"].as_f64().unwrap()
    };
    // Warm chr1/chunk0, then query chr2/chunk0 (same chunk_id) — must not alias.
    assert!((af("chr1", 100, "A", "G") - 0.111).abs() < 1e-6);
    assert!((af("chr2", 150, "C", "T") - 0.222).abs() < 1e-6);
    // And back again, from cache.
    assert!((af("chr1", 100, "A", "G") - 0.111).abs() < 1e-6);
}

#[test]
fn parallel_queries_do_not_thrash_with_adequate_cache() {
    let dir = TempDir::new().unwrap();
    // chunk_bits=8 → 256 positions/chunk (~2 KB each); 8000 positions → ~31
    // chunks so the 8-worker working set clearly exceeds a tiny cache.
    let n = 8000u32;
    let path = build(dir.path(), n, 8);

    let total_chunks = {
        let r = Osa2Reader::open(&path).unwrap();
        // Preload everything once, single-threaded, to count distinct chunks.
        let all: Vec<u64> = (0..n).map(|i| (1000 + i) as u64).collect();
        r.preload("chr1", &all).unwrap();
        r.chunk_load_count()
    };
    assert!(
        total_chunks > 4,
        "fixture must span several chunks, got {}",
        total_chunks
    );

    let positions: Vec<u64> = (0..n).step_by(3).map(|i| (1000 + i) as u64).collect();
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(8)
        .build()
        .unwrap();
    let sweep = |reader: &Osa2Reader| -> u64 {
        reader.preload("chr1", &positions).unwrap();
        pool.install(|| {
            positions.par_iter().for_each(|&p| {
                assert!(reader
                    .annotate_position("chr1", p, "A", "G")
                    .unwrap()
                    .is_some());
            });
        });
        reader.chunk_load_count()
    };

    // Chunk cache large enough for every chunk (512 MiB): the parallel phase
    // adds nothing on top of the preload — each chunk built exactly once.
    let big = Osa2Reader::open_with_cache_budget(&path, 512 * 1024 * 1024).unwrap();
    let big_loads = sweep(&big);
    assert_eq!(
        big_loads, total_chunks,
        "adequate cache should build each chunk exactly once"
    );

    // Tiny cache (holds ~1–2 of the ~2 KB chunks): the preloaded chunks are
    // evicted before the parallel phase and workers then evict one another's,
    // so chunks are rebuilt many times over.
    let tiny = Osa2Reader::open_with_cache_budget(&path, 4096).unwrap();
    let tiny_loads = sweep(&tiny);
    assert!(
        tiny_loads >= big_loads * 3,
        "small cache under parallelism must rebuild chunks: tiny={} big={} chunks={}",
        tiny_loads,
        big_loads,
        total_chunks
    );
}
