//! Issue #78: `--sa-dir` loading is parallel, so the contract that used to come
//! for free from a sequential `read_dir` loop now has to be asserted.
//!
//! Provider order decides the order of the supplementary output columns, so it
//! must be stable - and it must not depend on directory iteration order, which
//! differs between filesystems. These tests pin both: the order is the sorted
//! path order, and a source that fails to open is skipped rather than taking
//! the whole run down.

use fastvep_annotate::load_sa_providers;
use fastvep_sa::fields::{Field, FieldType};
use fastvep_sa::writer_v2::{Osa2Metadata, Osa2Record, Osa2Writer};
use std::path::Path;
use tempfile::TempDir;

fn write_osa2(dir: &Path, file_stem: &str, json_key: &str) {
    let fields = vec![Field {
        field: "AF".into(),
        alias: "af".into(),
        ftype: FieldType::Float,
        multiplier: 1_000_000,
        zigzag: false,
        missing_value: u32::MAX,
        missing_string: ".".into(),
        description: "allele frequency".into(),
    }];
    let metadata = Osa2Metadata {
        format_version: 2,
        name: json_key.into(),
        version: "test".into(),
        assembly: "GRCh38".into(),
        json_key: json_key.into(),
        match_by_allele: true,
        is_array: false,
        is_positional: false,
        chunk_bits: 20,
        description: String::new(),
    };
    let records = vec![Osa2Record {
        chrom: "chr1".into(),
        position: 1000,
        ref_allele: b"A".to_vec(),
        alt_allele: b"G".to_vec(),
        values: vec![fields[0].encode_float(0.25)],
        json_blob: None,
    }];
    let path = dir.join(format!("{file_stem}.osa2"));
    let file = std::fs::File::create(&path).unwrap();
    Osa2Writer::new(metadata, fields)
        .write_all(std::io::BufWriter::new(file), &records)
        .unwrap();
}

/// Written in an order that does not match the sorted order, so a loader that
/// simply preserved creation/iteration order would fail this.
#[test]
fn providers_are_ordered_by_path() {
    let dir = TempDir::new().unwrap();
    for (stem, key) in [
        ("gnomad_chr9", "k9"),
        ("gnomad_chr1", "k1"),
        ("gnomad_chr21", "k21"),
        ("gnomad_chr2", "k2"),
    ] {
        write_osa2(dir.path(), stem, key);
    }

    let providers = load_sa_providers(dir.path()).unwrap();
    let keys: Vec<String> = providers
        .iter()
        .map(|p| p.lock().unwrap().json_key().to_string())
        .collect();

    // Sorted by file name, which is lexicographic - chr1 < chr2 < chr21 < chr9.
    assert_eq!(keys, vec!["k1", "k2", "k21", "k9"]);
}

/// An unreadable source must be skipped with a warning, exactly as it was
/// before loading became parallel, and must not disturb the others' order.
#[test]
fn an_unopenable_source_is_skipped_not_fatal() {
    let dir = TempDir::new().unwrap();
    write_osa2(dir.path(), "a_good", "good_a");
    std::fs::write(dir.path().join("b_broken.osa2"), b"not a zip archive").unwrap();
    write_osa2(dir.path(), "c_good", "good_c");
    // Unrelated extensions in the same directory must be ignored.
    std::fs::write(dir.path().join("notes.txt"), b"ignore me").unwrap();

    let providers = load_sa_providers(dir.path()).unwrap();
    let keys: Vec<String> = providers
        .iter()
        .map(|p| p.lock().unwrap().json_key().to_string())
        .collect();
    assert_eq!(keys, vec!["good_a", "good_c"]);
}

#[test]
fn a_missing_directory_loads_nothing() {
    let dir = TempDir::new().unwrap();
    let providers = load_sa_providers(&dir.path().join("does_not_exist")).unwrap();
    assert!(providers.is_empty());
}
