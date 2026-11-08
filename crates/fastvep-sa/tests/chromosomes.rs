//! A database can say which chromosomes it holds.
//!
//! A database of some chromosomes only opens without complaint and answers "no record" for every
//! variant of the others, which for gnomAD is read as "absent from the population". The classifier
//! can tell the two apart only if the database says what it holds, so both readers answer
//! `AnnotationProvider::chromosomes`, from what they hold and not from what they were meant to.

use fastvep_cache::annotation::AnnotationProvider;
use fastvep_sa::common::AnnotationRecord;
use fastvep_sa::fields::{Field, FieldType};
use fastvep_sa::index::IndexHeader;
use fastvep_sa::reader::SaReader;
use fastvep_sa::reader_v2::Osa2Reader;
use fastvep_sa::writer::SaWriter;
use fastvep_sa::writer_v2::{Osa2Metadata, Osa2Record, Osa2Writer};

#[test]
fn a_v1_database_holds_the_chromosomes_that_have_a_record() {
    let dir = tempfile::tempdir().unwrap();
    let base = dir.path().join("two_of_three");
    let header = IndexHeader {
        schema_version: fastvep_sa::common::SCHEMA_VERSION,
        json_key: "gnomad".into(),
        name: "gnomAD".into(),
        version: "test".into(),
        description: String::new(),
        assembly: "GRCh38".into(),
        match_by_allele: true,
        is_array: false,
        is_positional: false,
    };
    // The chromosome list has three names and the records use two of them: what the database
    // holds is what has a record.
    let chrom_map = vec!["chr1".to_string(), "chr2".to_string(), "chrX".to_string()];
    let records = vec![
        AnnotationRecord {
            chrom_idx: 0,
            position: 100,
            ref_allele: "A".into(),
            alt_allele: "G".into(),
            json: "{}".into(),
        },
        AnnotationRecord {
            chrom_idx: 2,
            position: 100,
            ref_allele: "A".into(),
            alt_allele: "G".into(),
            json: "{}".into(),
        },
    ];
    SaWriter::new(header)
        .write_to_files(&base, records.into_iter(), &chrom_map)
        .unwrap();
    let reader = SaReader::open(&base.with_extension("osa")).unwrap();
    assert_eq!(
        reader.chromosomes(),
        Some(vec!["chr1".to_string(), "chrX".to_string()])
    );
}

#[test]
fn a_v2_database_holds_the_chromosomes_that_have_a_chunk() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("shards.osa2");
    let field = Field {
        field: "AF".into(),
        alias: "af".into(),
        ftype: FieldType::Float,
        multiplier: 1_000_000,
        zigzag: false,
        missing_value: u32::MAX,
        missing_string: ".".into(),
        description: String::new(),
    };
    let metadata = Osa2Metadata {
        format_version: 2,
        name: "gnomAD".into(),
        version: "test".into(),
        assembly: "GRCh38".into(),
        json_key: "gnomad".into(),
        match_by_allele: true,
        is_array: false,
        is_positional: false,
        chunk_bits: 16,
        description: String::new(),
    };
    let records: Vec<Osa2Record> = ["chr17", "chr2", "chrX"]
        .iter()
        .map(|chrom| Osa2Record {
            chrom: (*chrom).into(),
            position: 1000,
            ref_allele: b"A".to_vec(),
            alt_allele: b"G".to_vec(),
            values: vec![field.encode_float(0.1)],
            json_blob: None,
        })
        .collect();
    // The writer takes its records in the order of the chromosomes of the file.
    let mut sorted = records;
    sorted.sort_by(|a, b| a.chrom.cmp(&b.chrom));
    let file = std::fs::File::create(&path).unwrap();
    Osa2Writer::new(metadata, vec![field])
        .write_all(std::io::BufWriter::new(file), &sorted)
        .unwrap();
    let reader = Osa2Reader::open(&path).unwrap();
    assert_eq!(
        reader.chromosomes(),
        Some(vec![
            "chr17".to_string(),
            "chr2".to_string(),
            "chrX".to_string()
        ])
    );
}
