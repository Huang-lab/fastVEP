//! End-to-end round-trip tests for SaWriter -> SaReader.

use fastvep_cache::annotation::AnnotationProvider;
use fastvep_sa::common::AnnotationRecord;
use fastvep_sa::index::IndexHeader;
use fastvep_sa::reader::SaReader;
use fastvep_sa::writer::SaWriter;

#[test]
fn test_writer_reader_round_trip() {
    let dir = tempfile::tempdir().unwrap();
    let base = dir.path().join("test_clinvar");

    let header = IndexHeader {
        schema_version: fastvep_sa::common::SCHEMA_VERSION,
        json_key: "clinvar".into(),
        name: "ClinVar".into(),
        version: "2024-12".into(),
        description: "Test ClinVar data".into(),
        assembly: "GRCh38".into(),
        match_by_allele: true,
        is_array: false,
        is_positional: false,
    };

    // Create records sorted by (chrom_idx, position)
    let chrom_map = vec!["chr1".to_string(), "chr2".to_string()];
    let records: Vec<AnnotationRecord> = (0..1000)
        .map(|i| AnnotationRecord {
            chrom_idx: if i < 700 { 0 } else { 1 },
            position: if i < 700 { 1000 + i } else { 1000 + (i - 700) },
            ref_allele: "A".into(),
            alt_allele: "G".into(),
            json: format!(r#"{{"sig":"pathogenic","id":{}}}"#, i),
        })
        .collect();

    let mut writer = SaWriter::new(header);
    writer
        .write_to_files(&base, records.into_iter(), &chrom_map)
        .unwrap();

    // Verify files exist
    assert!(base.with_extension("osa").exists());
    assert!(base.with_extension("osa.idx").exists());

    // Open reader
    let reader = SaReader::open(&base.with_extension("osa")).unwrap();
    assert_eq!(reader.name(), "ClinVar");
    assert_eq!(reader.json_key(), "clinvar");

    // Query a known record on chr1
    let result = reader.annotate_position("chr1", 1050, "A", "G").unwrap();
    assert!(result.is_some());
    let json = match result.unwrap() {
        fastvep_cache::annotation::AnnotationValue::Json(j) => j,
        other => panic!("Expected Json, got {:?}", other),
    };
    assert!(json.contains(r#""sig":"pathogenic""#));
    assert!(json.contains(r#""id":50"#));

    // Query a known record on chr2
    let result = reader.annotate_position("chr2", 1100, "A", "G").unwrap();
    assert!(result.is_some());

    // Query with wrong allele (should not match when match_by_allele=true)
    let result = reader.annotate_position("chr1", 1050, "A", "T").unwrap();
    assert!(result.is_none());

    // Query non-existent position
    let result = reader.annotate_position("chr1", 99999, "A", "G").unwrap();
    assert!(result.is_none());

    // Query non-existent chromosome
    let result = reader.annotate_position("chr99", 1000, "A", "G").unwrap();
    assert!(result.is_none());
}

#[test]
fn test_preload_then_query() {
    let dir = tempfile::tempdir().unwrap();
    let base = dir.path().join("test_gnomad");

    let header = IndexHeader {
        schema_version: fastvep_sa::common::SCHEMA_VERSION,
        json_key: "gnomad".into(),
        name: "gnomAD".into(),
        version: "4.0".into(),
        description: "Test gnomAD data".into(),
        assembly: "GRCh38".into(),
        match_by_allele: true,
        is_array: false,
        is_positional: false,
    };

    let chrom_map = vec!["chr1".to_string()];
    let records: Vec<AnnotationRecord> = (0..500)
        .map(|i| AnnotationRecord {
            chrom_idx: 0,
            position: 10000 + i * 10,
            ref_allele: "C".into(),
            alt_allele: "T".into(),
            json: format!(r#"{{"af":{:.6}}}"#, i as f64 / 10000.0),
        })
        .collect();

    let mut writer = SaWriter::new(header);
    writer
        .write_to_files(&base, records.into_iter(), &chrom_map)
        .unwrap();

    let reader = SaReader::open(&base.with_extension("osa")).unwrap();

    // Preload a range of positions
    let positions: Vec<u64> = (0..50).map(|i| 10000 + i * 10).collect();
    reader.preload("chr1", &positions).unwrap();

    // Query preloaded positions
    let result = reader.annotate_position("chr1", 10100, "C", "T").unwrap();
    assert!(result.is_some());
}

#[test]
fn test_positional_annotation() {
    let dir = tempfile::tempdir().unwrap();
    let base = dir.path().join("test_phylop");

    let header = IndexHeader {
        schema_version: fastvep_sa::common::SCHEMA_VERSION,
        json_key: "phylopScore".into(),
        name: "PhyloP".into(),
        version: "1.0".into(),
        description: "PhyloP conservation scores".into(),
        assembly: "GRCh38".into(),
        match_by_allele: false,
        is_array: false,
        is_positional: true,
    };

    let chrom_map = vec!["chr1".to_string()];
    let records: Vec<AnnotationRecord> = (0..100)
        .map(|i| AnnotationRecord {
            chrom_idx: 0,
            position: 1000 + i,
            ref_allele: String::new(),
            alt_allele: String::new(),
            json: format!("{:.3}", (i as f64 - 50.0) / 10.0),
        })
        .collect();

    let mut writer = SaWriter::new(header);
    writer
        .write_to_files(&base, records.into_iter(), &chrom_map)
        .unwrap();

    let reader = SaReader::open(&base.with_extension("osa")).unwrap();

    // Positional: alleles don't matter
    let result = reader.annotate_position("chr1", 1025, "A", "G").unwrap();
    assert!(result.is_some());
    match result.unwrap() {
        fastvep_cache::annotation::AnnotationValue::Positional(s) => {
            assert_eq!(s, "-2.500");
        }
        other => panic!("Expected Positional, got {:?}", other),
    }
}

// =============================================================================
// V2 Format (.osa2) Round-Trip Tests
// =============================================================================

#[test]
fn test_osa2_writer_reader_round_trip() {
    use fastvep_sa::fields::{Field, FieldType};
    use fastvep_sa::reader_v2::Osa2Reader;
    use fastvep_sa::writer_v2::{Osa2Metadata, Osa2Record, Osa2Writer};

    let dir = tempfile::tempdir().unwrap();
    let osa2_path = dir.path().join("test_gnomad.osa2");

    let metadata = Osa2Metadata {
        format_version: 2,
        name: "gnomAD".into(),
        version: "4.1".into(),
        assembly: "GRCh38".into(),
        json_key: "gnomad".into(),
        match_by_allele: true,
        is_array: false,
        is_positional: false,
        chunk_bits: 20,
        description: "Test gnomAD data".into(),
    };

    let fields = vec![
        Field {
            field: "AF".into(),
            alias: "allAf".into(),
            ftype: FieldType::Float,
            multiplier: 2_000_000,
            zigzag: false,
            missing_value: u32::MAX,
            missing_string: ".".into(),
            description: String::new(),
        },
        Field {
            field: "AC".into(),
            alias: "allAc".into(),
            ftype: FieldType::Integer,
            multiplier: 1,
            zigzag: false,
            missing_value: u32::MAX,
            missing_string: ".".into(),
            description: String::new(),
        },
    ];

    // Create records: 500 variants on chr1 across multiple chunks
    let records: Vec<Osa2Record> = (0..500)
        .map(|i| {
            let pos = 10_000 + i * 100; // Spread across positions
            let af = (i as f64 + 1.0) / 100_000.0;
            Osa2Record {
                chrom: "chr1".into(),
                position: pos,
                ref_allele: b"A".to_vec(),
                alt_allele: b"G".to_vec(),
                values: vec![
                    fields[0].encode_float(af), // AF
                    (i + 1) as u32,             // AC
                ],
                json_blob: None,
            }
        })
        .collect();

    let writer = Osa2Writer::new(metadata.clone(), fields.clone());
    let file = std::fs::File::create(&osa2_path).unwrap();
    writer.write_all(file, &records).unwrap();

    // Open and query
    let reader = Osa2Reader::open(&osa2_path).unwrap();
    assert_eq!(reader.name(), "gnomAD");
    assert_eq!(reader.json_key(), "gnomad");

    // Query a known position
    let result = reader
        .annotate_position("chr1", 10_000 + 50 * 100, "A", "G")
        .unwrap();
    assert!(
        result.is_some(),
        "Expected annotation at position {}",
        10_000 + 50 * 100
    );
    let json = match result.unwrap() {
        fastvep_cache::annotation::AnnotationValue::Json(j) => j,
        other => panic!("Expected Json, got {:?}", other),
    };
    // Should contain allAf and allAc fields
    assert!(json.contains("\"allAf\":"), "JSON missing allAf: {}", json);
    assert!(
        json.contains("\"allAc\":51"),
        "JSON missing allAc:51: {}",
        json
    );

    // Query with wrong allele should miss
    let result = reader
        .annotate_position("chr1", 10_000 + 50 * 100, "A", "T")
        .unwrap();
    assert!(result.is_none(), "Wrong allele should not match");

    // Query non-existent position
    let result = reader.annotate_position("chr1", 99999, "A", "G").unwrap();
    assert!(result.is_none());
}

#[test]
fn test_osa2_preload() {
    use fastvep_sa::fields::{Field, FieldType};
    use fastvep_sa::reader_v2::Osa2Reader;
    use fastvep_sa::writer_v2::{Osa2Metadata, Osa2Record, Osa2Writer};

    let dir = tempfile::tempdir().unwrap();
    let osa2_path = dir.path().join("test_preload.osa2");

    let metadata = Osa2Metadata {
        format_version: 2,
        name: "test".into(),
        version: "1.0".into(),
        assembly: "GRCh38".into(),
        json_key: "test".into(),
        match_by_allele: true,
        is_array: false,
        is_positional: false,
        chunk_bits: 20,
        description: String::new(),
    };

    let fields = vec![Field {
        field: "score".into(),
        alias: "score".into(),
        ftype: FieldType::Integer,
        multiplier: 1,
        zigzag: false,
        missing_value: u32::MAX,
        missing_string: ".".into(),
        description: String::new(),
    }];

    let records: Vec<Osa2Record> = (0..100)
        .map(|i| Osa2Record {
            chrom: "chr1".into(),
            position: 1000 + i * 10,
            ref_allele: b"C".to_vec(),
            alt_allele: b"T".to_vec(),
            values: vec![i as u32 + 1],
            json_blob: None,
        })
        .collect();

    let writer = Osa2Writer::new(metadata, fields);
    let file = std::fs::File::create(&osa2_path).unwrap();
    writer.write_all(file, &records).unwrap();

    let reader = Osa2Reader::open(&osa2_path).unwrap();

    // Preload positions
    let positions: Vec<u64> = (0..50).map(|i| 1000 + i * 10).collect();
    reader.preload("chr1", &positions).unwrap();

    // Query preloaded position
    let result = reader.annotate_position("chr1", 1050, "C", "T").unwrap();
    assert!(result.is_some());
}
