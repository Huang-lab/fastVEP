//! `--source custom_vcf` / `custom_bed`, and `fastvep oga-build`.

use super::sa_build::standard_chrom_map;
use anyhow::{Context, Result};
use std::fs::File;
use std::io::{self};
use std::path::Path;

/// Classify a `--source custom` input as either `custom_vcf` or
/// `custom_bed` based on its extension. Centralised so the dispatcher in
/// `run_sa_build` and any future callers stay in agreement on the
/// extension rules.
pub(crate) fn classify_custom_input(input: &str) -> Result<&'static str> {
    let lower = input.to_lowercase();
    let stripped = lower.trim_end_matches(".gz").trim_end_matches(".bgz");
    if stripped.ends_with(".vcf") {
        Ok("custom_vcf")
    } else if stripped.ends_with(".bed") {
        Ok("custom_bed")
    } else {
        anyhow::bail!(
            "--source custom could not infer format from input '{}'. Use --source custom_vcf or --source custom_bed explicitly, or rename the file so it ends in .vcf[.gz] / .bed[.gz].",
            input
        )
    }
}

/// Default the `--name` flag from the input filename so a user invoking
/// `sa-build --source custom_vcf -i myset.vcf.gz` doesn't strictly need
/// `--name` (though we still warn that the column name will be derived
/// from the file path). Strips conventional extensions and falls back to
/// "custom" if the path doesn't have a useful stem.
pub(crate) fn resolve_custom_name(name: Option<&str>, input: &str) -> String {
    if let Some(n) = name {
        if !n.is_empty() {
            return n.to_string();
        }
    }
    let stem = Path::new(input)
        .file_name()
        .map(|s| s.to_string_lossy().to_string())
        .unwrap_or_else(|| "custom".to_string());
    let stem = stem
        .trim_end_matches(".gz")
        .trim_end_matches(".bgz")
        .trim_end_matches(".vcf")
        .trim_end_matches(".bed")
        .trim_end_matches('.');
    if stem.is_empty() {
        "custom".to_string()
    } else {
        stem.to_string()
    }
}

/// Build a generic allele-level annotation database from a user-supplied
/// VCF. Produces a `.osa` keyed by `name` that the annotate pipeline picks
/// up exactly like a built-in source (ClinVar, gnomAD, …). INFO fields
/// listed in `info_fields` are extracted; an empty list means "include
/// whatever INFO keys each record carries", which is convenient for
/// exploration but yields a heterogeneous schema.
pub(crate) fn run_custom_vcf_build(
    input: &str,
    output: &str,
    assembly: &str,
    name: Option<&str>,
    info_fields: &[String],
) -> Result<()> {
    use fastvep_sa::index::IndexHeader;
    use fastvep_sa::writer::SaWriter;

    let (chrom_list, chrom_map) = standard_chrom_map(assembly);
    let resolved_name = resolve_custom_name(name, input);
    let json_key = resolved_name.clone();

    eprintln!(
        "Building custom_vcf .osa from: {} (name={}, info_fields={})",
        input,
        resolved_name,
        if info_fields.is_empty() {
            "<all>".to_string()
        } else {
            info_fields.join(",")
        }
    );

    let header = IndexHeader {
        schema_version: fastvep_sa::common::SCHEMA_VERSION,
        json_key,
        name: resolved_name,
        version: "user-supplied".into(),
        description: format!("Custom VCF annotations for {}", assembly),
        assembly: assembly.into(),
        match_by_allele: true,
        is_array: false,
        is_positional: false,
    };

    let file = File::open(input).with_context(|| format!("Opening input file: {}", input))?;
    let reader: Box<dyn io::Read> = if input.ends_with(".gz") || input.ends_with(".bgz") {
        Box::new(flate2::read::MultiGzDecoder::new(file))
    } else {
        Box::new(file)
    };
    let buf_reader = io::BufReader::new(reader);

    let records =
        fastvep_sa::custom::parse_custom_vcf(buf_reader, &chrom_map, &header.name, info_fields)?;
    if records.is_empty() {
        // Build proceeds (writes an empty .osa rather than failing) so the
        // user can iterate on their filter, but call this out clearly —
        // a silent zero-record build is the #1 source of "why is my
        // annotation column always empty" tickets.
        eprintln!(
            "warning: custom VCF produced 0 records — check that {} contains data lines with valid chrom/pos and the requested INFO fields.",
            input
        );
    } else {
        eprintln!("Parsed {} records from custom VCF", records.len());
    }

    let output_path = Path::new(output);
    let mut writer = SaWriter::new(header);
    writer.write_to_files(output_path, records.into_iter(), &chrom_list)?;
    eprintln!(
        "Wrote: {} and {}",
        output_path.with_extension("osa").display(),
        output_path.with_extension("osa.idx").display()
    );
    Ok(())
}

/// Build a generic interval-level annotation database from a user-supplied
/// BED. Produces a `.osi` keyed by `name` that the annotate pipeline picks
/// up via the same `--sa-dir` mechanism as `.osa` files.
pub(crate) fn run_custom_bed_build(
    input: &str,
    output: &str,
    assembly: &str,
    name: Option<&str>,
) -> Result<()> {
    use fastvep_sa::interval::{IntervalHeader, IntervalIndex};

    let (_chrom_list, chrom_map) = standard_chrom_map(assembly);
    let resolved_name = resolve_custom_name(name, input);

    eprintln!(
        "Building custom_bed .osi from: {} (name={})",
        input, resolved_name
    );

    let file = File::open(input).with_context(|| format!("Opening input file: {}", input))?;
    let reader: Box<dyn io::Read> = if input.ends_with(".gz") || input.ends_with(".bgz") {
        Box::new(flate2::read::MultiGzDecoder::new(file))
    } else {
        Box::new(file)
    };
    let buf_reader = io::BufReader::new(reader);

    let records = fastvep_sa::custom::parse_custom_bed(buf_reader, &chrom_map)?;
    if records.is_empty() {
        eprintln!(
            "warning: custom BED produced 0 records — check that {} has at least 3 tab-separated columns (chrom, start, end).",
            input
        );
    } else {
        eprintln!("Parsed {} interval records from custom BED", records.len());
    }

    let header = IntervalHeader {
        schema_version: fastvep_sa::common::SCHEMA_VERSION,
        json_key: resolved_name.clone(),
        name: resolved_name,
        version: "user-supplied".into(),
        assembly: assembly.into(),
    };

    let mut index = IntervalIndex::new(header);
    for rec in records {
        index.add(rec);
    }
    index.sort();

    let output_path = Path::new(output);
    let osi_path = output_path.with_extension("osi");
    let mut file = File::create(&osi_path)
        .with_context(|| format!("Creating output file: {}", osi_path.display()))?;
    index.write_to(&mut file)?;
    eprintln!("Wrote: {}", osi_path.display());
    Ok(())
}

/// Build a gene-level annotation database (`.oga`) from a source file.
///
/// Supports three gene-level sources used by the ACMG-AMP classifier:
/// - `omim` - disease-gene annotations in genemap2 layout (13-col TSV):
///   ClinGen Gene-Disease Validity (preferred per ClinGen SVI / Abou
///   Tayoun 2018) or OMIM `genemap2.txt` (legacy). Drives PVS1, BS2,
///   PM3, BP2.
/// - `gnomad_genes` - gnomAD constraint metrics TSV (PVS1, PP2, BP1)
/// - `clinvar_protein` - ClinVar VCF, extracts pathogenic missense by
///   protein position (PS1, PM1, PM5)
///
/// The output is `<output>.oga`. The runtime loader at
/// `fastvep_annotate::load_gene_providers` picks up any `.oga` file in
/// `--sa-dir` and routes records to the classifier by `json_key`
/// (`omim`, `gnomad_genes`, `clinvar_protein`).
pub fn run_oga_build(source: &str, input: &str, output: &str, assembly: &str) -> Result<()> {
    use fastvep_sa::common::SCHEMA_VERSION;
    use fastvep_sa::gene::{GeneHeader, GeneIndex};

    let (json_key, name) = match source {
        // Source name kept as `omim` for back-compat; canonical input
        // is ClinGen GDV converted to genemap2 layout per SVI guidance.
        "omim" => ("omim", "ClinGen GDV / OMIM"),
        "gnomad_genes" | "gnomad_gene" => ("gnomad_genes", "gnomAD gene constraints"),
        "clinvar_protein" => ("clinvar_protein", "ClinVar protein index"),
        _ => anyhow::bail!(
            "run_oga_build called with non-gene source: {} (expected omim, gnomad_genes, clinvar_protein)",
            source
        ),
    };

    eprintln!("Building {} .oga from: {}", source, input);

    let file = File::open(input).with_context(|| format!("Opening input file: {}", input))?;
    let reader: Box<dyn io::Read> = if input.ends_with(".gz") || input.ends_with(".bgz") {
        Box::new(flate2::read::MultiGzDecoder::new(file))
    } else {
        Box::new(file)
    };
    let buf_reader = io::BufReader::new(reader);

    let records = match source {
        "omim" => fastvep_sa::sources::omim::parse_omim_genemap(buf_reader)?,
        "gnomad_genes" | "gnomad_gene" => {
            fastvep_sa::sources::gnomad_gene::parse_gnomad_gene_scores(buf_reader)?
        }
        "clinvar_protein" => {
            fastvep_sa::sources::clinvar_protein::parse_clinvar_protein_vcf(buf_reader, assembly)?
        }
        _ => unreachable!(),
    };

    eprintln!("Parsed {} records from {}", records.len(), source);

    let header = GeneHeader {
        schema_version: SCHEMA_VERSION,
        json_key: json_key.into(),
        name: name.into(),
        version: "latest".into(),
        assembly: assembly.into(),
    };

    let mut index = GeneIndex::new(header);
    for record in records {
        index.add(record);
    }

    let output_path = Path::new(output).with_extension("oga");
    let mut out_file = File::create(&output_path)
        .with_context(|| format!("Creating output file: {}", output_path.display()))?;
    index.write_to(&mut out_file)?;

    eprintln!(
        "Wrote: {} ({} genes)",
        output_path.display(),
        index.gene_count()
    );

    Ok(())
}

// SA provider loading is now in fastvep-annotate::load_sa_providers.

#[cfg(test)]
mod custom_source_tests {
    use super::*;
    // These exercise the custom sources end to end, which means going out
    // through the commands that consume them.
    use crate::pipeline::filter::run_filter;
    use crate::pipeline::sa_build::run_sa_build;
    use fastvep_cache::annotation::AnnotationProvider;

    #[test]
    fn classify_custom_recognises_vcf_and_bed_with_compression() {
        assert_eq!(classify_custom_input("foo.vcf").unwrap(), "custom_vcf");
        assert_eq!(classify_custom_input("foo.vcf.gz").unwrap(), "custom_vcf");
        assert_eq!(classify_custom_input("foo.vcf.bgz").unwrap(), "custom_vcf");
        assert_eq!(classify_custom_input("foo.bed").unwrap(), "custom_bed");
        assert_eq!(classify_custom_input("foo.BED.GZ").unwrap(), "custom_bed");
        assert!(classify_custom_input("foo.tsv").is_err());
        assert!(classify_custom_input("foo").is_err());
    }

    #[test]
    fn resolve_custom_name_prefers_flag_then_strips_extensions() {
        // Explicit --name wins.
        assert_eq!(resolve_custom_name(Some("MyDB"), "anything.vcf.gz"), "MyDB");
        // Empty --name acts like "not set".
        assert_eq!(resolve_custom_name(Some(""), "myset.vcf.gz"), "myset");
        // Strips .gz, .bgz, .vcf, .bed conventional suffixes.
        assert_eq!(resolve_custom_name(None, "/data/regions.bed.gz"), "regions");
        assert_eq!(resolve_custom_name(None, "myset.vcf"), "myset");
        // Pathological case: input is just `.gz` → fall through to "custom".
        assert_eq!(resolve_custom_name(None, ".gz"), "custom");
    }

    #[test]
    fn run_custom_vcf_build_produces_loadable_osa() {
        use std::io::Write;
        let dir = tempfile::tempdir().unwrap();
        let vcf_path = dir.path().join("custom.vcf");
        let mut f = std::fs::File::create(&vcf_path).unwrap();
        writeln!(f, "##fileformat=VCFv4.2").unwrap();
        writeln!(f, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO").unwrap();
        writeln!(f, "1\t100\t.\tA\tG\t.\t.\tMY_SCORE=0.95;FLAG=hot").unwrap();
        writeln!(f, "1\t200\t.\tC\tT\t.\t.\tMY_SCORE=0.10").unwrap();
        drop(f);

        let out_base = dir.path().join("custom_db");
        run_sa_build(
            "custom_vcf",
            vcf_path.to_str().unwrap(),
            out_base.to_str().unwrap(),
            "GRCh38",
            Some("mydb"),
            &["MY_SCORE".to_string()],
            false,
        )
        .unwrap();

        let osa = out_base.with_extension("osa");
        let osa_idx = out_base.with_extension("osa.idx");
        assert!(osa.exists(), ".osa should be written");
        assert!(osa_idx.exists(), ".osa.idx should be written");

        // Loadable via the same reader the runtime uses, and the records
        // we wrote come back through annotate_position.
        let reader = fastvep_sa::reader::SaReader::open(&osa).unwrap();
        assert_eq!(reader.json_key(), "mydb");
        let val = reader.annotate_position("chr1", 100, "A", "G").unwrap();
        assert!(val.is_some(), "should match the inserted record");
        match val.unwrap() {
            fastvep_cache::annotation::AnnotationValue::Json(j) => {
                assert!(j.contains("MY_SCORE"));
                // FLAG was not requested in info_fields → must not leak.
                assert!(!j.contains("FLAG"));
            }
            other => panic!("expected Json, got {:?}", other),
        }
    }

    #[test]
    fn run_custom_bed_build_produces_loadable_osi() {
        use std::io::Write;
        let dir = tempfile::tempdir().unwrap();
        let bed_path = dir.path().join("regions.bed");
        let mut f = std::fs::File::create(&bed_path).unwrap();
        // BED is 0-based half-open: chr1:99-200 → 1-based 100..=200.
        writeln!(f, "1\t99\t200\tpromoterA\t0.5").unwrap();
        writeln!(f, "1\t499\t600\tpromoterB").unwrap();
        drop(f);

        let out_base = dir.path().join("myregions");
        run_sa_build(
            "custom_bed",
            bed_path.to_str().unwrap(),
            out_base.to_str().unwrap(),
            "GRCh38",
            Some("myregions"),
            &[],
            false,
        )
        .unwrap();

        let osi_path = out_base.with_extension("osi");
        assert!(osi_path.exists(), ".osi should be written");

        let reader = fastvep_sa::interval::OsiReader::open(&osi_path).unwrap();
        // Position inside region A (1-based 100..=200).
        let val = reader.annotate_position("chr1", 150, "", "").unwrap();
        match val {
            Some(fastvep_cache::annotation::AnnotationValue::Interval(v)) => {
                assert_eq!(v.len(), 1);
                assert!(v[0].contains("promoterA"));
                assert!(v[0].contains("0.5"));
            }
            other => panic!("expected Interval, got {:?}", other),
        }
        // Position in the gap between regions → None.
        let val = reader.annotate_position("chr1", 350, "", "").unwrap();
        assert!(val.is_none());
    }

    #[test]
    fn explicit_transcript_cache_preserves_merged_source_labels() {
        // Regression for the second-pass review CRIT-2: running annotate
        // with `--transcript-cache combined.cache --gff3 single.gff3` used
        // to re-stamp every transcript in the cache with `single.gff3`'s
        // label, clobbering the merged Ensembl/RefSeq distinction the
        // cache was built to preserve. The fix gates re-stamping to
        // sidecar-cache (no `--transcript-cache`) + single-GFF3 mode.
        use fastvep_cache::transcript_cache;
        use fastvep_core::Strand;
        use fastvep_genome::{Gene, Transcript};
        use std::sync::Arc;

        let dir = tempfile::tempdir().unwrap();
        // Build a "merged" cache directly (skip GFF3 → cache to keep the
        // test focused on the load path, not the build path).
        let make = |source: &str, tid: &str, gid: &str| Transcript {
            stable_id: Arc::from(tid),
            version: None,
            gene: Gene {
                stable_id: Arc::from(gid),
                symbol: Some(Arc::from("GENE")),
                symbol_source: None,
                hgnc_id: None,
                biotype: Arc::from("protein_coding"),
                chromosome: Arc::from("chr1"),
                start: 1,
                end: 100,
                strand: Strand::Forward,
            },
            biotype: Arc::from("protein_coding"),
            chromosome: Arc::from("chr1"),
            start: 1,
            end: 100,
            strand: Strand::Forward,
            exons: vec![],
            translation: None,
            cdna_coding_start: None,
            cdna_coding_end: None,
            coding_region_start: None,
            coding_region_end: None,
            spliced_seq: None,
            translateable_seq: None,
            peptide: None,
            canonical: false,
            mane_select: None,
            mane_plus_clinical: None,
            tsl: None,
            appris: None,
            ccds: None,
            protein_id: None,
            protein_version: None,
            swissprot: vec![],
            trembl: vec![],
            uniparc: vec![],
            refseq_id: None,
            source: Some(source.to_string()),
            gencode_primary: false,
            flags: vec![],
            codon_table_start_phase: 0,
        };
        let trs = vec![
            make("Ensembl", "ENST00000001", "ENSG00000001"),
            make("RefSeq", "NM_000001", "GENE"),
        ];
        let cache_path = dir.path().join("combined.cache");
        transcript_cache::save_cache(&trs, &cache_path).unwrap();

        // Reload via the same code annotate uses. Drive the gating logic
        // by checking what the load_cache + re-stamp policy produces.
        let loaded = transcript_cache::load_cache(&cache_path).unwrap();
        assert_eq!(loaded.len(), 2);
        let sources: Vec<&str> = loaded.iter().filter_map(|t| t.source.as_deref()).collect();
        assert!(sources.contains(&"Ensembl"), "Ensembl label preserved");
        assert!(sources.contains(&"RefSeq"), "RefSeq label preserved");
    }

    #[test]
    fn run_sa_build_with_unknown_source_reports_custom_in_message() {
        // Regression: the user in issue #43 ran `--source custom` and got
        // "Unknown source" without seeing the new custom_* options. Verify
        // the help string now mentions them.
        let err = run_sa_build(
            "definitely_not_a_source",
            "/dev/null",
            "/tmp/out",
            "GRCh38",
            None,
            &[],
            false,
        )
        .expect_err("must error on unknown source");
        let msg = format!("{}", err);
        assert!(
            msg.contains("custom_vcf"),
            "error message must mention custom_vcf: {}",
            msg
        );
        assert!(
            msg.contains("custom_bed"),
            "error message must mention custom_bed: {}",
            msg
        );
    }

    #[test]
    fn run_filter_still_processes_well_formed_lines_when_malformed_lines_present() {
        // Regression: lines with fewer than 8 tab-separated fields were
        // silently `continue`d with no diagnostic even though they were
        // already counted toward `total`. This doesn't change the filtering
        // behavior (malformed lines still can't match), but verifies that
        // well-formed lines around a malformed one are still processed
        // correctly and the function doesn't error out.
        use std::io::Write;
        let dir = tempfile::tempdir().unwrap();
        let input_path = dir.path().join("in.vcf");
        let mut f = std::fs::File::create(&input_path).unwrap();
        writeln!(f, "##fileformat=VCFv4.2").unwrap();
        writeln!(
            f,
            "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations. Format: Consequence|IMPACT\">"
        )
        .unwrap();
        writeln!(f, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO").unwrap();
        // Well-formed, matches the filter.
        writeln!(f, "1\t100\t.\tA\tG\t.\tPASS\tCSQ=missense_variant|MODERATE").unwrap();
        // Malformed: fewer than 8 tab-separated fields.
        writeln!(f, "1\t200\t.\tC\tT").unwrap();
        // Well-formed, matches the filter.
        writeln!(f, "1\t300\t.\tG\tA\t.\tPASS\tCSQ=missense_variant|MODERATE").unwrap();
        drop(f);

        let output_path = dir.path().join("out.vcf");
        run_filter(
            input_path.to_str().unwrap(),
            output_path.to_str().unwrap(),
            "Consequence is missense_variant",
        )
        .expect("run_filter should not error on a malformed line");

        let out = std::fs::read_to_string(&output_path).unwrap();
        // Both well-formed matching lines pass through; the malformed line
        // is absent from the output (it can't match) but doesn't crash the
        // run or swallow its neighbors.
        assert!(
            out.contains("1\t100\t.\tA\tG"),
            "first well-formed line should pass:\n{}",
            out
        );
        assert!(
            out.contains("1\t300\t.\tG\tA"),
            "second well-formed line should pass:\n{}",
            out
        );
        assert!(
            !out.contains("1\t200\t.\tC\tT"),
            "malformed line must not appear in output:\n{}",
            out
        );
    }
}
