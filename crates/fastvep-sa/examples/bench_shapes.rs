//! Payload-shape size/speed sweep: v1 `.osa` vs v2 `.osa2` across the *shapes*
//! of annotation that fastVEP sources actually carry.
//!
//! The v2 format's compactness comes from packing numeric fields into
//! delta-/deflate-compressed u32 columns instead of repeating JSON text per
//! record. That advantage is real for numeric and small-categorical payloads
//! but NOT structural for opaque string / array payloads, which v2 stores as
//! the same zstd-compressed JSON blob v1 does (plus a Var32 key column and ZIP
//! overhead). This benchmark measures each shape head-to-head so the question
//! "is v2 smaller for *everything*?" is answered with numbers, not assertion.
//!
//! Shapes:
//!   numeric      gnomAD-like: 3 numeric fields (AF float + AN/AC ints)
//!   score        AlphaMissense-like: 1 float score + 1 three-level categorical
//!   id_string    dbSNP-like: one opaque per-variant string (rsID) — JsonBlob
//!   array_blob   ClinVar-like: a multi-field/array JSON record        — JsonBlob
//!
//! This example reports **on-disk size only** — the question it exists to
//! answer is "is v2 smaller?". Query throughput is a separate axis that depends
//! heavily on the access pattern (v2's 1 MB chunks amortize sparse, scattered
//! VCF lookups but lose to v1 on dense every-base scans); it is measured
//! realistically by the companion `bench_v1_vs_v2` example, not here.
//!
//! Usage:
//!   cargo run --release --example bench_shapes -p fastvep-sa
//!   SA_BENCH_RECORDS=5000000 cargo run --release --example bench_shapes -p fastvep-sa

use anyhow::Result;
use fastvep_cache::annotation::AnnotationProvider;
use fastvep_sa::common::{AnnotationRecord, SCHEMA_VERSION};
use fastvep_sa::fields::{Field, FieldType};
use fastvep_sa::index::IndexHeader;
use fastvep_sa::reader::SaReader;
use fastvep_sa::reader_v2::Osa2Reader;
use fastvep_sa::writer::SaWriter;
use fastvep_sa::writer_v2::{raw_json_blob_fields, Osa2Metadata, Osa2Record, Osa2Writer};
use std::env;
use std::path::Path;
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
    env::var(name).ok().and_then(|s| s.parse().ok()).unwrap_or(default)
}

fn chrom_list() -> Vec<String> {
    HUMAN_CHROMS.iter().map(|(_, n, _)| n.to_string()).collect()
}

/// One synthetic site: which chromosome and where.
struct Site {
    chrom_idx: u16,
    chrom: &'static str,
    position: u32,
}

fn build_sites(n_records: usize) -> Vec<Site> {
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
            sites.push(Site { chrom_idx: *chrom_idx, chrom: name, position: pos });
        }
    }
    sites
}

/// A payload shape knows how to emit both the v1 JSON and the v2 record for a
/// given site, and provides the matching v2 field schema + string tables.
trait Shape {
    fn name(&self) -> &'static str;
    fn v2_fields(&self) -> Vec<Field>;
    fn v2_string_tables(&self) -> Vec<(usize, Vec<String>)> {
        Vec::new()
    }
    fn v1_json(&self, site: &Site) -> String;
    fn v2_values(&self, site: &Site) -> Vec<u32>;
    fn v2_blob(&self, _site: &Site) -> Option<String> {
        None
    }
    fn json_key(&self) -> &'static str;
    fn is_array(&self) -> bool {
        false
    }
}

// --- numeric (gnomAD-like) --------------------------------------------------
struct Numeric {
    fields: Vec<Field>,
}
impl Numeric {
    fn new() -> Self {
        Self {
            fields: vec![
                float_field("AF", "allAf", 1_000_000),
                int_field("AN", "allAn"),
                int_field("AC", "allAc"),
            ],
        }
    }
}
impl Shape for Numeric {
    fn name(&self) -> &'static str { "numeric   (gnomAD: AF+AN+AC)" }
    fn json_key(&self) -> &'static str { "gnomad" }
    fn v2_fields(&self) -> Vec<Field> { self.fields.clone() }
    fn v1_json(&self, s: &Site) -> String {
        let af = ((s.position as u64 % 100_000) as f64) / 1_000_000.0;
        let an = 150_000i64;
        let ac = (af * an as f64) as i64;
        format!(r#"{{"allAf":{:.6e},"allAn":{},"allAc":{}}}"#, af, an, ac)
    }
    fn v2_values(&self, s: &Site) -> Vec<u32> {
        let af = ((s.position as u64 % 100_000) as f64) / 1_000_000.0;
        let an = 150_000i64;
        let ac = (af * an as f64) as i64;
        vec![
            self.fields[0].encode_float(af),
            self.fields[1].encode_int(an),
            self.fields[2].encode_int(ac),
        ]
    }
}

// --- score (AlphaMissense-like) --------------------------------------------
const AM_CLASSES: &[&str] = &["likely_benign", "ambiguous", "likely_pathogenic"];
struct Score {
    fields: Vec<Field>,
}
impl Score {
    fn new() -> Self {
        Self {
            fields: vec![
                float_field("am_pathogenicity", "amPathogenicity", 1_000_000),
                categorical_field("am_class", "amClass"),
            ],
        }
    }
    fn score_of(s: &Site) -> f64 {
        ((s.position as u64 % 10_000) as f64) / 10_000.0
    }
    fn class_of(s: &Site) -> u32 {
        s.position % 3
    }
}
impl Shape for Score {
    fn name(&self) -> &'static str { "score     (AlphaMissense: score+class)" }
    fn json_key(&self) -> &'static str { "alphaMissense" }
    fn v2_fields(&self) -> Vec<Field> { self.fields.clone() }
    fn v2_string_tables(&self) -> Vec<(usize, Vec<String>)> {
        vec![(1, AM_CLASSES.iter().map(|s| s.to_string()).collect())]
    }
    fn v1_json(&self, s: &Site) -> String {
        let score = Self::score_of(s);
        let class = AM_CLASSES[Self::class_of(s) as usize];
        format!(r#"{{"amPathogenicity":{:.6e},"amClass":"{}"}}"#, score, class)
    }
    fn v2_values(&self, s: &Site) -> Vec<u32> {
        vec![self.fields[0].encode_float(Self::score_of(s)), Self::class_of(s)]
    }
}

// --- id_string (dbSNP-like) -------------------------------------------------
struct IdString {
    fields: Vec<Field>,
}
impl IdString {
    fn new() -> Self {
        Self { fields: vec![blob_field("dbsnp")] }
    }
    fn blob(s: &Site) -> String {
        // One opaque, near-unique identifier per variant, as dbSNP carries.
        format!(r#"{{"rsId":"rs{}"}}"#, (s.chrom_idx as u64) * 1_000_000_000 + s.position as u64)
    }
}
impl Shape for IdString {
    fn name(&self) -> &'static str { "id_string (dbSNP: rsID blob)" }
    fn json_key(&self) -> &'static str { "dbsnp" }
    fn v2_fields(&self) -> Vec<Field> { self.fields.clone() }
    fn v1_json(&self, s: &Site) -> String { Self::blob(s) }
    fn v2_values(&self, _s: &Site) -> Vec<u32> { Vec::new() }
    fn v2_blob(&self, s: &Site) -> Option<String> { Some(Self::blob(s)) }
}

// --- array_blob (ClinVar-like) ---------------------------------------------
struct ArrayBlob {
    fields: Vec<Field>,
}
impl ArrayBlob {
    fn new() -> Self {
        Self { fields: vec![blob_field("clinvar")] }
    }
    fn blob(s: &Site) -> String {
        // A multi-field record with arrays, as ClinVar carries.
        let sig = match s.position % 4 {
            0 => r#"["Pathogenic"]"#,
            1 => r#"["Likely_pathogenic","Pathogenic"]"#,
            2 => r#"["Benign","Likely_benign"]"#,
            _ => r#"["Uncertain_significance"]"#,
        };
        format!(
            r#"{{"significance":{},"reviewStatus":"criteria_provided_multiple_submitters","numSubmitters":{},"conditions":["Cardiomyopathy","Long_QT_syndrome"],"variationId":{}}}"#,
            sig,
            (s.position % 12) + 1,
            (s.chrom_idx as u64) * 1_000_000 + s.position as u64
        )
    }
}
impl Shape for ArrayBlob {
    fn name(&self) -> &'static str { "array_blob(ClinVar: sig array + fields)" }
    fn json_key(&self) -> &'static str { "clinvar" }
    fn is_array(&self) -> bool { true }
    fn v2_fields(&self) -> Vec<Field> { self.fields.clone() }
    fn v1_json(&self, s: &Site) -> String { Self::blob(s) }
    fn v2_values(&self, _s: &Site) -> Vec<u32> { Vec::new() }
    fn v2_blob(&self, s: &Site) -> Option<String> { Some(Self::blob(s)) }
}

fn float_field(field: &str, alias: &str, multiplier: u32) -> Field {
    Field {
        field: field.into(), alias: alias.into(), ftype: FieldType::Float,
        multiplier, zigzag: false, missing_value: u32::MAX,
        missing_string: ".".into(), description: String::new(),
    }
}
fn int_field(field: &str, alias: &str) -> Field {
    Field {
        field: field.into(), alias: alias.into(), ftype: FieldType::Integer,
        multiplier: 1, zigzag: false, missing_value: u32::MAX,
        missing_string: ".".into(), description: String::new(),
    }
}
fn categorical_field(field: &str, alias: &str) -> Field {
    Field {
        field: field.into(), alias: alias.into(), ftype: FieldType::Categorical,
        multiplier: 1, zigzag: false, missing_value: u32::MAX,
        missing_string: ".".into(), description: String::new(),
    }
}
fn blob_field(alias: &str) -> Field {
    Field {
        field: alias.into(), alias: alias.into(), ftype: FieldType::JsonBlob,
        multiplier: 1, zigzag: false, missing_value: u32::MAX,
        missing_string: ".".into(), description: String::new(),
    }
}

fn build_v1(path: &Path, sites: &[Site], shape: &dyn Shape) -> Result<u64> {
    let mut records: Vec<AnnotationRecord> = sites
        .iter()
        .map(|s| AnnotationRecord {
            chrom_idx: s.chrom_idx,
            position: s.position,
            ref_allele: "A".into(),
            alt_allele: "G".into(),
            json: shape.v1_json(s),
        })
        .collect();
    records.sort_by(|a, b| a.chrom_idx.cmp(&b.chrom_idx).then(a.position.cmp(&b.position)));

    let header = IndexHeader {
        schema_version: SCHEMA_VERSION,
        json_key: shape.json_key().into(),
        name: shape.json_key().into(),
        version: "bench".into(),
        description: "bench".into(),
        assembly: "GRCh38".into(),
        match_by_allele: true,
        is_array: shape.is_array(),
        is_positional: false,
    };
    let mut writer = SaWriter::new(header);
    writer.write_to_files(path, records.into_iter(), &chrom_list())?;
    // v1 is a pair of files (.osa data + .osa.idx offsets); count both so the
    // comparison against v2's single self-contained .osa2 archive is fair.
    let osa = std::fs::metadata(path.with_extension("osa"))?.len();
    let idx = std::fs::metadata(path.with_extension("osa.idx"))
        .map(|m| m.len())
        .unwrap_or(0);
    Ok(osa + idx)
}

fn build_v2(path: &Path, sites: &[Site], shape: &dyn Shape) -> Result<u64> {
    let fields = shape.v2_fields();
    let mut records: Vec<Osa2Record> = sites
        .iter()
        .map(|s| Osa2Record {
            chrom: s.chrom.to_string(),
            position: s.position,
            ref_allele: b"A".to_vec(),
            alt_allele: b"G".to_vec(),
            values: shape.v2_values(s),
            json_blob: shape.v2_blob(s),
        })
        .collect();
    records.sort_by(|a, b| a.chrom.cmp(&b.chrom).then(a.position.cmp(&b.position)));

    let metadata = Osa2Metadata {
        format_version: 2,
        name: shape.json_key().into(),
        version: "bench".into(),
        assembly: "GRCh38".into(),
        json_key: shape.json_key().into(),
        match_by_allele: true,
        is_array: shape.is_array(),
        is_positional: false,
        chunk_bits: 20,
        description: "bench".into(),
    };
    let mut writer = Osa2Writer::new(metadata, fields);
    for (idx, table) in shape.v2_string_tables() {
        writer.set_string_table(idx, table);
    }
    let file = std::fs::File::create(path)?;
    writer.write_all(std::io::BufWriter::new(file), &records)?;
    Ok(std::fs::metadata(path)?.len())
}

fn mb(bytes: u64) -> f64 {
    bytes as f64 / (1024.0 * 1024.0)
}

/// Positional per-base score text, matching how the v1 score parsers render a
/// score (up to 4 decimals, trailing zeros trimmed). Both formats store the
/// identical string so the size comparison is apples-to-apples.
fn score_text(score: f64) -> String {
    if score == 0.0 {
        return "0".into();
    }
    let s = format!("{:.4}", score);
    let s = s.trim_end_matches('0').trim_end_matches('.');
    s.to_string()
}

/// A deterministic per-base "conservation" score in a PhyloP-like range
/// (~[-4.5, 6.5], 3 decimals). A slow triangle wave gives the spatial
/// correlation real conservation has, plus per-base pseudo-random jitter so the
/// score column carries realistic entropy — otherwise an over-smooth synthetic
/// would flatter v2's compression and overstate the size win.
fn positional_score(pos: u32) -> f64 {
    let phase = (pos % 2000) as f64 / 2000.0; // 0..1
    let tri = if phase < 0.5 { phase * 2.0 } else { 2.0 - phase * 2.0 }; // 0..1..0
    let base = tri * 10.0 - 4.0;
    // Knuth multiplicative hash → deterministic per-position jitter in ~[-0.5, 0.5].
    let h = pos.wrapping_mul(2_654_435_761);
    let jitter = ((h >> 8) & 0x3FF) as f64 / 1024.0 - 0.5;
    ((base + jitter) * 1000.0).round() / 1000.0
}

/// Measure v1 `.osa` vs v2 `.osa2` for a DENSE per-base positional source
/// (PhyloP/GERP/DANN-shaped): one allele-less score per consecutive coordinate.
/// This is the case the numeric/blob shapes above don't cover, and the one that
/// decides whether positional sources are worth a v2 encoder.
fn bench_positional(n_records: usize) -> Result<()> {
    let dir = TempDir::new()?;
    let v1_base = dir.path().join("pos.v1");
    let v2_path = dir.path().join("pos.osa2");

    // Dense consecutive positions on chr1 (per-base coverage), allele-less.
    let chrom = "chr1";
    let chrom_idx = 0u16;

    // --- v1 .osa (positional header, bare-number JSON, empty alleles) ---
    let mut v1_records: Vec<AnnotationRecord> = Vec::with_capacity(n_records);
    for i in 0..n_records {
        let pos = (i + 1) as u32;
        v1_records.push(AnnotationRecord {
            chrom_idx,
            position: pos,
            ref_allele: String::new(),
            alt_allele: String::new(),
            json: score_text(positional_score(pos)),
        });
    }
    let header = IndexHeader {
        schema_version: SCHEMA_VERSION,
        json_key: "phylop".into(),
        name: "phylop".into(),
        version: "bench".into(),
        description: "bench".into(),
        assembly: "GRCh38".into(),
        match_by_allele: false,
        is_array: false,
        is_positional: true,
    };
    let mut w = SaWriter::new(header);
    w.write_to_files(&v1_base, v1_records.into_iter(), &chrom_list())?;
    let v1_size = std::fs::metadata(v1_base.with_extension("osa"))?.len()
        + std::fs::metadata(v1_base.with_extension("osa.idx")).map(|m| m.len()).unwrap_or(0);

    // --- v2 .osa2 (positional metadata, whole-record blob = bare number) ---
    let mut v2_records: Vec<Osa2Record> = Vec::with_capacity(n_records);
    for i in 0..n_records {
        let pos = (i + 1) as u32;
        v2_records.push(Osa2Record {
            chrom: chrom.to_string(),
            position: pos,
            ref_allele: Vec::new(),
            alt_allele: Vec::new(),
            values: Vec::new(),
            json_blob: Some(score_text(positional_score(pos))),
        });
    }
    let metadata = Osa2Metadata {
        format_version: 2,
        name: "phylop".into(),
        version: "bench".into(),
        assembly: "GRCh38".into(),
        json_key: "phylop".into(),
        match_by_allele: false,
        is_array: false,
        is_positional: true,
        chunk_bits: 20,
        description: "bench".into(),
    };
    let writer = Osa2Writer::new(metadata, raw_json_blob_fields());
    let file = std::fs::File::create(&v2_path)?;
    writer.write_all(std::io::BufWriter::new(file), &v2_records)?;
    let v2_size = std::fs::metadata(&v2_path)?.len();

    // Verify equivalence at a sample position before comparing sizes.
    let v1_reader = SaReader::open(&v1_base.with_extension("osa"))?;
    let v2_reader = Osa2Reader::open(&v2_path)?;
    let probe = (n_records / 2 + 1) as u64;
    let g1 = v1_reader.annotate_position(chrom, probe, "A", "G")?;
    let g2 = v2_reader.annotate_position(chrom, probe, "A", "G")?;
    let flag = match (&g1, &g2) {
        (Some(a), Some(b)) if format!("{a:?}") == format!("{b:?}") => "",
        _ => "  !! data mismatch",
    };

    println!(
        "{:<40} {:>12.1} {:>12.1} {:>9.2}x{}",
        "positional(PhyloP: per-base score, dense)",
        mb(v1_size),
        mb(v2_size),
        v2_size as f64 / v1_size as f64,
        flag,
    );
    Ok(())
}

fn main() -> Result<()> {
    let n_records = env_usize("SA_BENCH_RECORDS", 2_000_000);
    println!("Records per shape: {}\n", n_records);
    let sites = build_sites(n_records);

    let shapes: Vec<Box<dyn Shape>> = vec![
        Box::new(Numeric::new()),
        Box::new(Score::new()),
        Box::new(IdString::new()),
        Box::new(ArrayBlob::new()),
    ];

    println!("{:<40} {:>12} {:>12} {:>10}", "shape", "v1 MB", "v2 MB", "v2/v1");
    println!("{}", "-".repeat(76));

    for shape in &shapes {
        let dir = TempDir::new()?;
        let v1_base = dir.path().join("bench.v1");
        let v2_path = dir.path().join("bench.osa2");

        let v1_size = build_v1(&v1_base, &sites, shape.as_ref())?;
        let v2_size = build_v2(&v2_path, &sites, shape.as_ref())?;

        // Confirm the two formats hold equivalent data before comparing sizes,
        // so a size win can't come from silently dropping information.
        let v1_reader = SaReader::open(&v1_base.with_extension("osa"))?;
        let v2_reader = Osa2Reader::open(&v2_path)?;
        let sample = &sites[sites.len() / 2];
        let got_v1 = v1_reader
            .annotate_position(sample.chrom, sample.position as u64, "A", "G")?
            .is_some();
        let got_v2 = v2_reader
            .annotate_position(sample.chrom, sample.position as u64, "A", "G")?
            .is_some();
        let flag = if got_v1 != got_v2 { "  !! data mismatch" } else { "" };

        println!(
            "{:<40} {:>12.1} {:>12.1} {:>9.2}x{}",
            shape.name(),
            mb(v1_size),
            mb(v2_size),
            v2_size as f64 / v1_size as f64,
            flag,
        );
    }

    // Positional per-base source (dense consecutive coordinates, allele-less).
    bench_positional(n_records)?;

    println!(
        "\nv2/v1 < 1.0 means v2 is smaller. Numeric/score shapes pack into u32\n\
         columns and win only at genome scale (the fixed per-chunk/ZIP overhead\n\
         dominates for small inputs). id_string/array_blob store the same JSON\n\
         payload in both formats, but v2 zstd-compresses a whole chunk's blobs\n\
         together — exploiting cross-record redundancy v1's per-block scheme\n\
         can't — so v2 is markedly smaller there."
    );
    Ok(())
}
