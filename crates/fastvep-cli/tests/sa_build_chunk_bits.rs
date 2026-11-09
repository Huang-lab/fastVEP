//! `sa-build --chunk-bits`: the width of the chunks of an `.osa2`.
//!
//! A query inflates and indexes the whole chunk it falls in, so on a source with hundreds of
//! thousands of records in a megabase every lookup that misses the cache is dear: 1,000 scattered
//! GRCh38 variants took 24 s against the dbSNP of the hosted site, 23 ms a lookup, and 0.5 s without
//! it. Narrowing the chunks is the remedy, and these tests hold the two things it must not break: that
//! the file says what was asked, and that a query returns what it returned before.

use fastvep_annotate::load_sa_providers_reporting;
use fastvep_cli::pipeline::{run_sa_build_v2, with_chunk_bits};
use std::fs;
use std::path::Path;

const LONG_ALT: &str = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";

/// The records: positions that fall in different chunks at 16 bits (65,536 bases) and in the same
/// one at 20 (1,048,576), neighbours on both sides of a 16-bit boundary, a deletion that starts in
/// one 16-bit chunk and ends in the next, an insertion too long to be packed, and two alleles at one
/// position.
fn records() -> Vec<(u32, &'static str, &'static str)> {
    vec![
        (100, "A", "G"),
        (200, "A", "G"),
        (65_530, "ACGTACGTAC", "A"),
        (65_535, "A", "G"),
        (65_536, "A", "G"),
        (65_700, "A", "G"),
        (140_000, "A", "G"),
        (140_100, "A", LONG_ALT),
        (200_000, "A", "C"),
        (200_000, "A", "G"),
        (1_048_700, "A", "G"),
        (1_048_800, "A", "G"),
    ]
}

fn vcf() -> String {
    let mut out = String::from(
        "##fileformat=VCFv4.2\n##fileDate=2026-11-02\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
    );
    for (i, (pos, reference, alt)) in records().iter().enumerate() {
        out.push_str(&format!(
            "1\t{pos}\trs{}\t{reference}\t{alt}\t.\t.\t.\n",
            1000 + i
        ));
    }
    out
}

fn build(dir: &Path, bits: Option<u32>) -> Vec<u8> {
    let input = dir.join("dbsnp.vcf");
    fs::write(&input, vcf()).unwrap();
    let sa = dir.join("sa");
    fs::create_dir_all(&sa).unwrap();
    with_chunk_bits(bits, || {
        run_sa_build_v2(
            "dbsnp",
            input.to_str().unwrap(),
            sa.join("dbsnp").to_str().unwrap(),
            "GRCh38",
            None,
            &[],
            false,
        )
    })
    .unwrap();
    fs::read(sa.join("dbsnp.osa2")).unwrap()
}

fn has_entry(archive: &[u8], name: &str) -> bool {
    archive
        .windows(name.len())
        .any(|window| window == name.as_bytes())
}

fn answers(dir: &Path) -> Vec<String> {
    let (providers, skipped) = load_sa_providers_reporting(&dir.join("sa")).unwrap();
    assert!(skipped.is_empty(), "{skipped:?}");
    assert_eq!(providers.len(), 1);
    // Every record, and alleles and positions that hold nothing: the neighbours of the boundary,
    // the same position with another allele, and a position in a chunk nothing was written in.
    let mut probes: Vec<(u32, String, String)> = records()
        .iter()
        .map(|(pos, r, a)| (*pos, r.to_string(), a.to_string()))
        .collect();
    for (pos, r, a) in [
        (99, "A", "G"),
        (101, "A", "G"),
        (65_534, "A", "G"),
        (65_537, "A", "G"),
        (200_000, "A", "T"),
        (65_530, "ACGTACGTAC", "C"),
        (1_048_699, "A", "G"),
        (2_000_000, "A", "G"),
    ] {
        probes.push((pos, r.to_string(), a.to_string()));
    }
    probes
        .iter()
        .map(|(pos, r, a)| {
            let short = if a.len() > 12 { "<long>" } else { a.as_str() };
            format!(
                "{pos} {r}>{short}: {:?}",
                providers[0]
                    .annotate_position("1", u64::from(*pos), r, a)
                    .unwrap()
            )
        })
        .collect()
}

#[test]
fn a_narrower_chunk_is_what_the_file_says_and_every_query_answers_as_before() {
    let tmp = tempfile::tempdir().unwrap();
    let wide_dir = tmp.path().join("wide");
    let narrow_dir = tmp.path().join("narrow");
    fs::create_dir(&wide_dir).unwrap();
    fs::create_dir(&narrow_dir).unwrap();

    let wide = build(&wide_dir, None);
    let narrow = build(&narrow_dir, Some(16));

    // The source's own width is one megabase: the records of 1,048,700 and 1,048,800 are in chunk 1,
    // and there is no chunk 16.
    assert!(has_entry(&wide, "fastsa/chr1/0/var32.bin"));
    assert!(has_entry(&wide, "fastsa/chr1/1/var32.bin"));
    assert!(!has_entry(&wide, "fastsa/chr1/16/var32.bin"));
    // At 16 bits they are in chunk 16, 65,535 in chunk 0, 65,536 and 65,700 in chunk 1, 140,000 in 2
    // and 200,000 in 3.
    for chunk in [0, 1, 2, 3, 16] {
        assert!(
            has_entry(&narrow, &format!("fastsa/chr1/{chunk}/var32.bin")),
            "chunk {chunk}"
        );
    }
    assert!(!has_entry(&narrow, "fastsa/chr1/4/var32.bin"));

    let before = answers(&wide_dir);
    let after = answers(&narrow_dir);
    assert_eq!(before, after);
    // And they are answers: the records were found (the first twelve probes), and what was not
    // written was not.
    for (i, line) in before.iter().enumerate().take(records().len()) {
        assert!(line.contains(&format!("rs{}", 1000 + i)), "{line}");
    }
    for line in before.iter().skip(records().len()) {
        assert!(line.ends_with("None"), "{line}");
    }
}

#[test]
fn the_width_is_for_the_builds_inside_the_scope_alone() {
    let tmp = tempfile::tempdir().unwrap();
    let a = tmp.path().join("a");
    let b = tmp.path().join("b");
    fs::create_dir(&a).unwrap();
    fs::create_dir(&b).unwrap();
    let narrow = build(&a, Some(16));
    let after = build(&b, None);
    assert!(has_entry(&narrow, "fastsa/chr1/16/var32.bin"));
    assert!(
        !has_entry(&after, "fastsa/chr1/16/var32.bin"),
        "a build after the scope kept the narrow width"
    );
}
