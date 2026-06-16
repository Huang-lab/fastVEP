//! Integration test: a repeat-region indel that gnomAD stores in its
//! left-aligned representation is missed by a raw right-anchored query, but
//! matches once the query is run through `normalize_variant`. This is the
//! end-to-end proof for the PM2 frequency fix (gnomAD indels silently failing
//! to match → spurious "absent" → spurious PM2).

use anyhow::{anyhow, Result};
use fastvep_cache::annotation::{AnnotationProvider, AnnotationValue};
use fastvep_cache::normalize::normalize_variant;
use fastvep_cache::providers::SequenceProvider;
use fastvep_sa::common::{AnnotationRecord, SCHEMA_VERSION};
use fastvep_sa::index::IndexHeader;
use fastvep_sa::reader::SaReader;
use fastvep_sa::writer::SaWriter;

/// 1-based reference for a single contig: pos1=G then "TAAC" repeated.
struct StrRef;
impl SequenceProvider for StrRef {
    fn fetch_sequence(&self, chrom: &str, start: u64, end: u64) -> Result<Vec<u8>> {
        const SEQ: &[u8] = b"GTAACTAACTAACTAACG"; // pos1=G, then TAAC x4, then G
        if chrom != "chr2" {
            return Err(anyhow!("unknown contig"));
        }
        if start < 1 || end < start {
            return Err(anyhow!("bad range"));
        }
        let s0 = (start - 1) as usize;
        if s0 >= SEQ.len() {
            return Err(anyhow!("past contig end"));
        }
        let e = (end as usize).min(SEQ.len());
        Ok(SEQ[s0..e].to_vec())
    }
}

#[test]
fn gnomad_repeat_indel_matches_after_normalization() {
    let dir = tempfile::tempdir().unwrap();
    let base = dir.path().join("gnomad");

    let header = IndexHeader {
        schema_version: SCHEMA_VERSION,
        json_key: "gnomad".into(),
        name: "gnomAD".into(),
        version: "v4".into(),
        description: "Test gnomAD".into(),
        assembly: "GRCh38".into(),
        match_by_allele: true,
        is_array: false,
        is_positional: false,
    };

    // gnomAD stores the deletion in its LEFT-ALIGNED form: chr2:1 GTAAC>G.
    let records = vec![AnnotationRecord {
        chrom_idx: 0,
        position: 1,
        ref_allele: "GTAAC".into(),
        alt_allele: "G".into(),
        json: r#"{"allAf":0.012,"allAc":300}"#.into(),
    }];
    let chrom_map = vec!["chr2".to_string()];

    let mut writer = SaWriter::new(header);
    writer
        .write_to_files(&base, records.into_iter(), &chrom_map)
        .unwrap();
    let reader = SaReader::open(&base.with_extension("osa")).unwrap();

    // The SAME variant written right-anchored mid-repeat (chr2:5 CTAAC>C) — as
    // it might appear in an input VCF / ClinVar — does NOT match: today's bug.
    let raw = reader.annotate_position("chr2", 5, "CTAAC", "C").unwrap();
    assert!(raw.is_none(), "raw right-anchored query should miss gnomAD");

    // Normalize the raw query to gnomAD's minimal/left-aligned representation,
    // then query again — now it matches and we recover the real AF.
    let n = normalize_variant(&StrRef, "chr2", 5, "CTAAC", "C");
    assert_eq!((n.pos, n.ref_allele.as_str(), n.alt_allele.as_str()), (1, "GTAAC", "G"));

    let hit = reader
        .annotate_position("chr2", n.pos, &n.ref_allele, &n.alt_allele)
        .unwrap();
    assert!(hit.is_some(), "normalized query should match gnomAD");
    match hit.unwrap() {
        AnnotationValue::Json(j) => assert!(j.contains("\"allAf\":0.012")),
        other => panic!("expected Json, got {:?}", other),
    }
}
