//! A chunk that cannot be built is attempted once, not once per variant.
//!
//! Every lookup in a chunk's genomic width asked for the chunk, and a chunk
//! whose blobs fail to decompress paid the full decompression before failing
//! each time. One unreadable megabase of a dense database was enough to make a
//! run look hung rather than slow: 200,000 variants over three such chunks did
//! not finish in ten minutes, where remembering the verdict finishes in
//! seconds (#101).

use fastvep_cache::annotation::AnnotationProvider;
use fastvep_sa::fields::{Field, FieldType};
use fastvep_sa::reader_v2::Osa2Reader;
use fastvep_sa::writer_v2::{Osa2Metadata, Osa2Record, Osa2Writer};

fn blob_field() -> Vec<Field> {
    vec![Field {
        field: String::new(),
        alias: String::new(),
        ftype: FieldType::JsonBlob,
        multiplier: 1,
        zigzag: false,
        missing_value: u32::MAX,
        missing_string: ".".into(),
        description: String::new(),
    }]
}

/// Corrupt the payload of one ZIP entry in place, leaving every offset and
/// length in the archive untouched so the reader still finds the entry and only
/// fails once it tries to decompress it.
///
/// A local file header is: 4-byte signature, 22 bytes of fixed fields, then a
/// 2-byte name length and a 2-byte extra-field length, then the name, the extra
/// field, and the data.
fn corrupt_entry_payload(bytes: &mut [u8], name: &str) -> bool {
    let sig = b"PK\x03\x04";
    let name_bytes = name.as_bytes();
    let mut i = 0;
    while i + 30 + name_bytes.len() <= bytes.len() {
        if &bytes[i..i + 4] == sig {
            let name_len = u16::from_le_bytes([bytes[i + 26], bytes[i + 27]]) as usize;
            let extra_len = u16::from_le_bytes([bytes[i + 28], bytes[i + 29]]) as usize;
            let name_start = i + 30;
            if name_len == name_bytes.len()
                && &bytes[name_start..name_start + name_len] == name_bytes
            {
                let data = name_start + name_len + extra_len;
                // Flip enough leading bytes that neither a zstd frame header
                // nor a deflate stream can be read out of them.
                for b in bytes.iter_mut().skip(data).take(16) {
                    *b ^= 0xFF;
                }
                return true;
            }
        }
        i += 1;
    }
    false
}

#[test]
fn an_unreadable_chunk_is_attempted_once() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("blobs.osa2");

    // 64 records, all in chunk 0 of chr1.
    let records: Vec<Osa2Record> = (0..64u32)
        .map(|i| Osa2Record {
            chrom: "chr1".into(),
            position: 1_000 + i,
            ref_allele: b"A".to_vec(),
            alt_allele: b"G".to_vec(),
            values: Vec::new(),
            json_blob: Some(format!("{{\"score\":{i}}}")),
        })
        .collect();
    Osa2Writer::new(
        Osa2Metadata {
            format_version: 2,
            name: "Blobs".into(),
            version: "1".into(),
            assembly: "GRCh38".into(),
            json_key: "blobs".into(),
            match_by_allele: true,
            is_array: false,
            is_positional: false,
            chunk_bits: 20,
            description: "test".into(),
        },
        blob_field(),
    )
    .write_all(std::fs::File::create(&path).unwrap(), &records)
    .unwrap();

    // Sanity: intact, the record reads back.
    {
        let reader = Osa2Reader::open(&path).unwrap();
        assert!(reader
            .annotate_position("chr1", 1_010, "A", "G")
            .unwrap()
            .is_some());
        assert_eq!(reader.chunk_failure_count(), 0);
    }

    let mut bytes = std::fs::read(&path).unwrap();
    assert!(
        corrupt_entry_payload(&mut bytes, "fastsa/chr1/0/json_blobs.zst"),
        "the blob entry should be present to corrupt"
    );
    let broken = dir.path().join("broken.osa2");
    std::fs::write(&broken, &bytes).unwrap();

    let reader = Osa2Reader::open(&broken).unwrap();
    let mut errors = 0;
    for i in 0..500u32 {
        // Every one of these lands in the same unreadable chunk.
        if reader
            .annotate_position("chr1", (1_000 + i % 64) as u64, "A", "G")
            .is_err()
        {
            errors += 1;
        }
    }

    // Every lookup must still report the failure - it is not a miss.
    assert_eq!(errors, 500, "each lookup in a broken chunk must report it");
    // But the chunk itself is only ever attempted once.
    assert_eq!(
        reader.chunk_failure_count(),
        1,
        "the failed chunk should be remembered, not rebuilt per lookup"
    );
    // And nothing was cached as a successful load.
    assert_eq!(reader.chunk_load_count(), 0);
}
