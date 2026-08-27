//! `fastvep sa-build` / `sa-convert`: build supplementary annotation files.

use super::custom::{
    classify_custom_input, resolve_custom_name, run_custom_bed_build, run_custom_vcf_build,
    run_oga_build,
};
use anyhow::{Context, Result};
use flate2::read::MultiGzDecoder;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufWriter, Read, Seek, SeekFrom};
use std::path::Path;

/// Standard chromosome ordering for SA builds.
///
/// Canonical names use the UCSC `chr*` style so the resulting `.osa.idx`
/// keys match the GRCh38 convention used by modern gnomAD / ClinVar
/// releases — and by most input VCFs we annotate against. The HashMap
/// resolves both `chr*` and bare forms (plus `MT`/`M`) so source parsers
/// can hand us either style. See issue #37.
pub(crate) fn standard_chrom_map(
    assembly: &str,
) -> (Vec<String>, std::collections::HashMap<String, u16>) {
    let chroms: Vec<String> = (1..=22)
        .map(|i| format!("chr{}", i))
        .chain(["chrX", "chrY", "chrM"].iter().map(|s| s.to_string()))
        .collect();
    let mut map: std::collections::HashMap<String, u16> = chroms
        .iter()
        .enumerate()
        .map(|(i, c)| (c.clone(), i as u16))
        .collect();
    // Accept bare-style aliases (e.g. NCBI / 1000G), mapping to the same
    // canonical index so the on-disk key still ends up `chr*`.
    for (i, c) in chroms.iter().enumerate() {
        if let Some(bare) = c.strip_prefix("chr") {
            map.insert(bare.to_string(), i as u16);
        }
    }
    // Mitochondrial aliases: canonical is `chrM`; accept `MT` and `chrMT`.
    if let Some(&mt_idx) = map.get("chrM") {
        map.insert("MT".to_string(), mt_idx);
        map.insert("chrMT".to_string(), mt_idx);
    }
    // RefSeq molecule accessions (`NC_000001.11`). NCBI's dbSNP VCF release
    // names contigs this way, so without these every record's chromosome misses
    // the map and zero records are parsed (issue #51). The accession→chromosome
    // relationship is assembly-specific and non-algorithmic, hence a table.
    if let Some(accessions) = fastvep_core::refseq_primary_accessions(assembly) {
        for (acc, chr_name) in accessions {
            if let Some(&idx) = map.get(*chr_name) {
                map.insert(acc.to_string(), idx);
            }
        }
    }
    (chroms, map)
}

/// Open an `sa-build` input file and transparently decompress gzip/bgzip.
///
/// Detection is by the gzip magic bytes (`0x1f 0x8b`) first, falling back to the
/// `.gz`/`.bgz` extension — mirroring `wrap_maybe_gzip_reader` on the annotate
/// path so a bgzipped source with a non-standard name (e.g. `gnomad.sites.vcf`)
/// still decodes instead of silently parsing to zero records. The file is
/// seekable, so we sniff and rewind rather than chaining the peeked bytes.
///
/// When `byte_counter` is supplied the raw (still-compressed) file is wrapped in
/// a `CountingReader` *before* decompression, so progress reflects compressed
/// bytes consumed against the on-disk file size.
fn open_sa_input(
    input: &str,
    byte_counter: Option<std::sync::Arc<std::sync::atomic::AtomicU64>>,
) -> Result<Box<dyn io::Read>> {
    let mut file = File::open(input).with_context(|| format!("Opening input file: {}", input))?;
    let mut magic = [0u8; 2];
    let n = file.read(&mut magic)?;
    file.seek(SeekFrom::Start(0))?;
    let is_gzip =
        (n == 2 && magic == [0x1f, 0x8b]) || input.ends_with(".gz") || input.ends_with(".bgz");

    let raw: Box<dyn io::Read> = match byte_counter {
        Some(bytes) => Box::new(crate::progress::CountingReader { inner: file, bytes }),
        None => Box::new(file),
    };
    if is_gzip {
        Ok(Box::new(MultiGzDecoder::new(raw)))
    } else {
        Ok(raw)
    }
}

/// Stream a coordinate-sorted source VCF straight into the .osa writer without
/// buffering every record in memory (issue #55: gnomAD/TOPMed/dbSNP releases
/// carry 100M+ records). The input is wrapped in a `CountingReader` so the
/// `ProgressMeter` can report byte-based % done + ETA for gz/bgz inputs whose
/// compressed size is known; plain or size-unknown inputs fall back to the
/// records-only meter.
///
/// `make_iter` adapts the (already gz-decoded) reader into a source-specific
/// streaming iterator; all callers must yield records in `(chrom_idx, position)`
/// order because the writer streams blocks and rejects out-of-order input.
fn run_streaming_sa_build<'a, I>(
    input: &str,
    output: &str,
    header: fastvep_sa::index::IndexHeader,
    chrom_map: &'a std::collections::HashMap<String, u16>,
    chrom_list: &[String],
    show_progress: bool,
    make_iter: impl FnOnce(
        io::BufReader<Box<dyn io::Read>>,
        &'a std::collections::HashMap<String, u16>,
    ) -> I,
) -> Result<()>
where
    I: Iterator<Item = Result<fastvep_sa::common::AnnotationRecord>>,
{
    use std::sync::atomic::AtomicU64;
    use std::sync::Arc;

    let output_path = Path::new(output);
    let mut writer = fastvep_sa::writer::SaWriter::new(header);

    let file_size = std::fs::metadata(input).map(|m| m.len()).unwrap_or(0);
    let byte_counter = Arc::new(AtomicU64::new(0));
    let reader = open_sa_input(input, Some(Arc::clone(&byte_counter)))?;
    let buf_reader = io::BufReader::new(reader);

    let mut meter = if show_progress && file_size > 0 {
        crate::progress::ProgressMeter::with_progress(true, file_size, byte_counter)
    } else {
        crate::progress::ProgressMeter::new(show_progress)
    };

    let records = make_iter(buf_reader, chrom_map).inspect(|record| {
        if record.is_ok() {
            meter.update();
        }
    });

    // The streaming writer creates and fills the .osa incrementally, so a
    // mid-stream failure (a parse error, or the writer's own out-of-order
    // `bail!`) leaves a truncated .osa and never writes the .idx. The old
    // buffer-then-write path parsed everything first and so left nothing behind
    // on failure; preserve that contract by removing the partial artifacts
    // rather than leaving a corrupt half-built database on disk.
    if let Err(e) = writer.write_results_to_files(output_path, records, chrom_list) {
        let _ = std::fs::remove_file(output_path.with_extension("osa"));
        let _ = std::fs::remove_file(output_path.with_extension("osa.idx"));
        return Err(e);
    }

    meter.finish();
    eprintln!(
        "Wrote: {} and {}",
        output_path.with_extension("osa").display(),
        output_path.with_extension("osa.idx").display()
    );

    Ok(())
}

/// PhyloP comes in two on-disk formats: UCSC fixed-step wig (`fixedStep
/// chrom=...`) and a simple `chrom\tpos\tscore` TSV (which is what we
/// emit when distilling PhyloP from gnomAD v4 INFO). Detect by peeking
/// the first non-blank line; everything else falls through to the TSV
/// parser since it tolerates BED-4 too. Streams either format straight
/// into the writer — PhyloP/GERP are per-base, genome-wide sources
/// (~3 billion positions for hg38) so buffering them as a Vec first
/// would exhaust memory long before a real build finishes.
fn iter_phylop_auto<'a, R: BufRead + 'a>(
    mut reader: R,
    chrom_to_idx: &'a HashMap<String, u16>,
) -> Box<dyn Iterator<Item = Result<fastvep_sa::common::AnnotationRecord>> + 'a> {
    let mut peek = [0u8; 16];
    let n = match reader.read(&mut peek) {
        Ok(n) => n,
        Err(e) => return Box::new(std::iter::once(Err(e).context("Peeking phyloP input"))),
    };
    let prefix = std::str::from_utf8(&peek[..n]).unwrap_or("");
    let chained = std::io::Cursor::new(peek[..n].to_vec()).chain(reader);
    let buf = io::BufReader::new(chained);
    if prefix.trim_start().starts_with("fixedStep") {
        Box::new(fastvep_sa::sources::scores::iter_wigfix(buf, chrom_to_idx))
    } else {
        Box::new(fastvep_sa::sources::scores::iter_score_tsv(
            buf,
            chrom_to_idx,
            false,
        ))
    }
}

/// Build a supplementary annotation .osa file from a source VCF.
pub fn run_sa_build(
    source: &str,
    input: &str,
    output: &str,
    assembly: &str,
    name: Option<&str>,
    info_fields: &[String],
    show_progress: bool,
) -> Result<()> {
    use fastvep_sa::index::IndexHeader;
    use fastvep_sa::writer::SaWriter;

    // Gene-level sources (.oga) — dispatched separately from variant-level (.osa).
    if matches!(
        source,
        "omim" | "gnomad_genes" | "gnomad_gene" | "clinvar_protein"
    ) {
        return run_oga_build(source, input, output, assembly);
    }

    // Custom VCF / BED — generic user-supplied annotation sources. `custom`
    // is an alias that auto-detects from the input extension. We dispatch
    // these before the built-in header-lookup so they don't have to thread
    // through every `match source` arm below.
    if matches!(source, "custom_vcf" | "custom_bed" | "custom") {
        let resolved = if source == "custom" {
            classify_custom_input(input)?
        } else {
            source
        };
        return match resolved {
            "custom_vcf" => run_custom_vcf_build(input, output, assembly, name, info_fields),
            "custom_bed" => run_custom_bed_build(input, output, assembly, name),
            _ => unreachable!("classify_custom_input only returns custom_vcf or custom_bed"),
        };
    }

    let (chrom_list, chrom_map) = standard_chrom_map(assembly);

    let header = match source {
        "clinvar" => IndexHeader {
            schema_version: fastvep_sa::common::SCHEMA_VERSION,
            json_key: "clinvar".into(),
            name: "ClinVar".into(),
            version: "latest".into(),
            description: format!("ClinVar annotations for {}", assembly),
            assembly: assembly.into(),
            match_by_allele: true,
            is_array: true,
            is_positional: false,
        },
        "gnomad" => IndexHeader {
            schema_version: fastvep_sa::common::SCHEMA_VERSION,
            json_key: "gnomad".into(),
            name: "gnomAD".into(),
            version: "latest".into(),
            description: format!("gnomAD population frequencies for {}", assembly),
            assembly: assembly.into(),
            match_by_allele: true,
            is_array: false,
            is_positional: false,
        },
        "dbsnp" => IndexHeader {
            schema_version: fastvep_sa::common::SCHEMA_VERSION,
            json_key: "dbsnp".into(),
            name: "dbSNP".into(),
            version: "latest".into(),
            description: format!("dbSNP RS IDs for {}", assembly),
            assembly: assembly.into(),
            match_by_allele: true,
            is_array: false,
            is_positional: false,
        },
        "phylop" | "gerp" | "dann" => {
            // Match the json_key the classifier looks for in
            // `fastvep_classification::sa_extract` (`phylop`, `gerp`).
            let json_key = match source {
                "phylop" => "phylop",
                "gerp" => "gerp",
                "dann" => "dann",
                _ => unreachable!(),
            };
            IndexHeader {
                schema_version: fastvep_sa::common::SCHEMA_VERSION,
                json_key: json_key.into(),
                name: source.to_uppercase(),
                version: "latest".into(),
                description: format!("{} conservation/prediction scores for {}", source, assembly),
                assembly: assembly.into(),
                match_by_allele: false,
                is_array: false,
                is_positional: true,
            }
        },
        "revel" => IndexHeader {
            schema_version: fastvep_sa::common::SCHEMA_VERSION,
            json_key: "revel".into(),
            name: "REVEL".into(),
            version: "latest".into(),
            description: format!("REVEL missense pathogenicity scores for {}", assembly),
            assembly: assembly.into(),
            match_by_allele: true,
            is_array: false,
            is_positional: false,
        },
        "spliceai" => IndexHeader {
            schema_version: fastvep_sa::common::SCHEMA_VERSION,
            json_key: "spliceAI".into(),
            name: "SpliceAI".into(),
            version: "latest".into(),
            description: format!("SpliceAI splice site predictions for {}", assembly),
            assembly: assembly.into(),
            match_by_allele: true,
            is_array: false,
            is_positional: false,
        },
        "primateai" => IndexHeader {
            schema_version: fastvep_sa::common::SCHEMA_VERSION,
            json_key: "primateAI".into(),
            name: "PrimateAI".into(),
            version: "latest".into(),
            description: format!("PrimateAI pathogenicity predictions for {}", assembly),
            assembly: assembly.into(),
            match_by_allele: true,
            is_array: false,
            is_positional: false,
        },
        "dbnsfp" => IndexHeader {
            schema_version: fastvep_sa::common::SCHEMA_VERSION,
            json_key: "dbnsfp".into(),
            name: "dbNSFP".into(),
            version: "latest".into(),
            description: format!("dbNSFP SIFT/PolyPhen predictions for {}", assembly),
            assembly: assembly.into(),
            match_by_allele: true,
            is_array: false,
            is_positional: false,
        },
        "cosmic" => IndexHeader {
            schema_version: fastvep_sa::common::SCHEMA_VERSION,
            json_key: "cosmic".into(),
            name: "COSMIC".into(),
            version: "latest".into(),
            description: format!("COSMIC somatic mutations for {}", assembly),
            assembly: assembly.into(),
            match_by_allele: true,
            is_array: false,
            is_positional: false,
        },
        "onekg" | "1000g" => IndexHeader {
            schema_version: fastvep_sa::common::SCHEMA_VERSION,
            json_key: "oneKg".into(),
            name: "1000 Genomes".into(),
            version: "latest".into(),
            description: format!("1000 Genomes population frequencies for {}", assembly),
            assembly: assembly.into(),
            match_by_allele: true,
            is_array: false,
            is_positional: false,
        },
        "topmed" => IndexHeader {
            schema_version: fastvep_sa::common::SCHEMA_VERSION,
            json_key: "topmed".into(),
            name: "TOPMed".into(),
            version: "latest".into(),
            description: format!("TOPMed population frequencies for {}", assembly),
            assembly: assembly.into(),
            match_by_allele: true,
            is_array: false,
            is_positional: false,
        },
        "alphamissense" => IndexHeader {
            schema_version: fastvep_sa::common::SCHEMA_VERSION,
            json_key: "alphaMissense".into(),
            name: "AlphaMissense".into(),
            version: "latest".into(),
            description: format!("AlphaMissense pathogenicity predictions for {}", assembly),
            assembly: assembly.into(),
            match_by_allele: true,
            is_array: false,
            is_positional: false,
        },
        "mitomap" => IndexHeader {
            schema_version: fastvep_sa::common::SCHEMA_VERSION,
            json_key: "mitomap".into(),
            name: "MitoMap".into(),
            version: "latest".into(),
            description: format!("MitoMap mitochondrial variants for {}", assembly),
            assembly: assembly.into(),
            match_by_allele: true,
            is_array: false,
            is_positional: false,
        },
        _ => anyhow::bail!(
            "Unknown source: {}. Supported: clinvar, gnomad, dbsnp, cosmic, onekg, topmed, mitomap, phylop, gerp, dann, revel, spliceai, primateai, dbnsfp, alphamissense, omim, gnomad_genes, clinvar_protein, custom_vcf, custom_bed, custom",
            source
        ),
    };

    eprintln!(
        "Building {} .osa: {} -> {}",
        source,
        input,
        Path::new(output).with_extension("osa").display()
    );

    // Large, coordinate-sorted population/score sources stream straight into
    // the writer instead of buffering every record (issue #55). They share one
    // helper that differs only in which per-source iterator adapts the reader.
    match source {
        "spliceai" => {
            return run_streaming_sa_build(
                input,
                output,
                header,
                &chrom_map,
                &chrom_list,
                show_progress,
                |r, m| fastvep_sa::sources::spliceai::iter_spliceai_vcf(r, m),
            );
        }
        "gnomad" => {
            return run_streaming_sa_build(
                input,
                output,
                header,
                &chrom_map,
                &chrom_list,
                show_progress,
                |r, m| fastvep_sa::sources::gnomad::iter_gnomad_vcf(r, m),
            );
        }
        "dbsnp" => {
            return run_streaming_sa_build(
                input,
                output,
                header,
                &chrom_map,
                &chrom_list,
                show_progress,
                |r, m| fastvep_sa::sources::dbsnp::iter_dbsnp_vcf(r, m),
            );
        }
        "topmed" => {
            return run_streaming_sa_build(
                input,
                output,
                header,
                &chrom_map,
                &chrom_list,
                show_progress,
                |r, m| fastvep_sa::sources::topmed::iter_topmed_vcf(r, m),
            );
        }
        "alphamissense" => {
            return run_streaming_sa_build(
                input,
                output,
                header,
                &chrom_map,
                &chrom_list,
                show_progress,
                |r, m| fastvep_sa::sources::alphamissense::iter_alphamissense_tsv(r, m),
            );
        }
        "phylop" => {
            return run_streaming_sa_build(
                input,
                output,
                header,
                &chrom_map,
                &chrom_list,
                show_progress,
                iter_phylop_auto,
            );
        }
        "gerp" | "dann" => {
            return run_streaming_sa_build(
                input,
                output,
                header,
                &chrom_map,
                &chrom_list,
                show_progress,
                |r, m| fastvep_sa::sources::scores::iter_score_tsv(r, m, false),
            );
        }
        _ => {}
    }

    let buf_reader = io::BufReader::new(open_sa_input(input, None)?);

    eprintln!(
        "INFO  ProgressMeter - Parsing '{}': {} -> {}",
        source,
        input,
        Path::new(output).with_extension("osa").display()
    );
    let t_parse = std::time::Instant::now();

    let records = match source {
        "clinvar" => fastvep_sa::sources::clinvar::parse_clinvar_vcf(buf_reader, &chrom_map)?,
        "cosmic" => fastvep_sa::sources::cosmic::parse_cosmic_vcf(buf_reader, &chrom_map)?,
        "onekg" | "1000g" => fastvep_sa::sources::onekg::parse_onekg_vcf(buf_reader, &chrom_map)?,
        "mitomap" => fastvep_sa::sources::mitomap::parse_mitomap(buf_reader, &chrom_map)?,
        "revel" => fastvep_sa::sources::revel::parse_revel(buf_reader, &chrom_map, 2)?,
        "primateai" => fastvep_sa::sources::primateai::parse_primateai(buf_reader, &chrom_map)?,
        "dbnsfp" => fastvep_sa::sources::dbnsfp::parse_dbnsfp(buf_reader, &chrom_map)?,
        _ => unreachable!(),
    };

    let n = records.len() as u64;
    if records.is_empty() {
        eprintln!(
            "Warning: 0 records parsed from {} — if the input is non-empty, this \
             usually means none of its chromosome names matched assembly '{}'. \
             NCBI releases (e.g. dbSNP) name contigs by RefSeq accession \
             (NC_000001.11); pass the --assembly (GRCh38 or GRCh37) matching the \
             build you downloaded.",
            source, assembly
        );
    }
    eprintln!(
        "INFO  ProgressMeter - Parsed {} records in {:.1} min. Writing ...",
        crate::progress::fmt_count(n),
        t_parse.elapsed().as_secs_f64() / 60.0
    );

    let output_path = Path::new(output);
    let t_write = std::time::Instant::now();
    let mut writer = SaWriter::new(header);
    writer.write_to_files(output_path, records.into_iter(), &chrom_list)?;

    eprintln!(
        "INFO  ProgressMeter - Done. Wrote {} records in {:.1} min.",
        crate::progress::fmt_count(n),
        t_write.elapsed().as_secs_f64() / 60.0
    );
    eprintln!(
        "Wrote: {} and {}",
        output_path.with_extension("osa").display(),
        output_path.with_extension("osa.idx").display()
    );

    Ok(())
}

/// The allele-level sources that build to the v2 `.osa2` format. Single source
/// of truth for the `--format auto` dispatch and the `--format osa2` error
/// message — adding a source to v2 means adding it here plus a match arm in
/// [`run_sa_build_v2`]. Includes the `1000g` alias of `onekg`.
///
/// This is every variant-level source, so `--format auto` builds v2 for all of
/// them and the v1 `.osa` writer is reachable only via an explicit
/// `--format osa`. Absent by nature, not by omission:
///
/// * gene-level sources (`omim`, `gnomad_genes`, `clinvar_protein`) are keyed by
///   gene symbol and build `.oga`;
/// * `custom_bed` is interval-keyed and builds `.osi`.
///
/// v2 is a variant-level container, so neither has a v2 form. Positional
/// per-base scores (PhyloP/GERP/DANN) ARE supported — keyed by coordinate via
/// `var32::positional_key`.
pub const OSA2_SUPPORTED_SOURCES: &[&str] = &[
    "gnomad",
    "onekg",
    "1000g",
    "topmed",
    "alphamissense",
    "dbsnp",
    "cosmic",
    "clinvar",
    "revel",
    "primateai",
    "dbnsfp",
    "phylop",
    "gerp",
    "dann",
    "spliceai",
    "mitomap",
    "custom_vcf",
];

/// The subset of [`OSA2_SUPPORTED_SOURCES`] whose v2 encoder decomposes the
/// payload into parallel u32 value columns (plus, for some, a categorical string
/// table). The rest store one whole-record JSON blob per variant, because their
/// payloads are high-cardinality ID strings or nested arrays that the numeric
/// layout cannot represent — see [`fastvep_sa::writer_v2::raw_json_blob_fields`].
///
/// The distinction matters to `sa-convert`: a `.osa` retains no field schema, so
/// a conversion can only ever produce the blob encoding. For a source in this
/// list, rebuilding from upstream yields a materially smaller file and the tool
/// says so; for the others, converting is equivalent to a rebuild.
pub const OSA2_DECOMPOSED_SOURCES: &[&str] = &[
    "gnomad",
    "onekg",
    "1000g",
    "topmed",
    "alphamissense",
    "spliceai",
];

/// Whether `--format osa2` (and `--format auto`) can build this source to v2.
pub fn source_supports_osa2(source: &str) -> bool {
    OSA2_SUPPORTED_SOURCES.contains(&source)
}

/// Whether this source's v2 encoder decomposes into numeric columns rather than
/// storing whole-record JSON blobs. See [`OSA2_DECOMPOSED_SOURCES`].
pub fn source_has_decomposed_osa2(source: &str) -> bool {
    OSA2_DECOMPOSED_SOURCES.contains(&source)
}

/// Recover the `--source` name from a database's `json_key`.
///
/// Most sources use their own name as the key, but several camel-case it
/// (`alphaMissense`, `oneKg`, `primateAI`, `spliceAI`), so the match is
/// case-insensitive. Returns `None` for a key that is not a built-in source —
/// notably a `custom_vcf` database, whose key is whatever the user passed to
/// `--name`.
pub fn source_from_json_key(json_key: &str) -> Option<&'static str> {
    let lowered = json_key.to_ascii_lowercase();
    OSA2_SUPPORTED_SOURCES
        .iter()
        .copied()
        .find(|s| *s == lowered.as_str())
}

/// Dispatch `sa-build` to the v1 or v2 builder according to `--format`:
///
/// * `auto` (the default) — build v2 for sources that support it (smaller and
///   faster to query at genome scale), v1 for everything else. This gives users
///   the higher-quality format per source without having to know which is which.
/// * `osa` / `v1` — force v1 `.osa`.
/// * `osa2` / `v2` — force v2 `.osa2` (errors if the source has no v2 encoder).
// Each argument is an independent coordinate, allele or flag with no
// natural grouping; bundling them into a struct would only move the
// argument list to the call site.
#[allow(clippy::too_many_arguments)]
pub fn run_sa_build_format(
    format: &str,
    source: &str,
    input: &str,
    output: &str,
    assembly: &str,
    name: Option<&str>,
    info_fields: &[String],
    show_progress: bool,
) -> Result<()> {
    match format {
        "osa" | "v1" => run_sa_build(
            source,
            input,
            output,
            assembly,
            name,
            info_fields,
            show_progress,
        ),
        "osa2" | "v2" => run_sa_build_v2(
            source,
            input,
            output,
            assembly,
            name,
            info_fields,
            show_progress,
        ),
        "auto" => {
            if source_supports_osa2(source) {
                run_sa_build_v2(
                    source,
                    input,
                    output,
                    assembly,
                    name,
                    info_fields,
                    show_progress,
                )
            } else {
                run_sa_build(
                    source,
                    input,
                    output,
                    assembly,
                    name,
                    info_fields,
                    show_progress,
                )
            }
        }
        other => anyhow::bail!(
            "Unknown --format '{}': expected 'auto' (default; best format per source), \
             'osa' (v1), or 'osa2' (v2)",
            other
        ),
    }
}

/// Build a v2 (`.osa2`) supplementary annotation database.
///
/// The v2 format (chunked ZIP, u32 value arrays, Var32 binary search) is
/// faster to query and smaller on disk at genome scale than v1 `.osa`. It is
/// wired for the allele-level sources in [`OSA2_SUPPORTED_SOURCES`]: numeric
/// payloads (gnomAD, 1000 Genomes, TOPMed, AlphaMissense) encode into parallel
/// u32 columns, while opaque string/array payloads (dbSNP, COSMIC, ClinVar)
/// are stored as whole-record JSON blobs. Other sources continue to build v1
/// `.osa` via [`run_sa_build`].
pub fn run_sa_build_v2(
    source: &str,
    input: &str,
    output: &str,
    assembly: &str,
    name: Option<&str>,
    info_fields: &[String],
    show_progress: bool,
) -> Result<()> {
    use fastvep_sa::sources::{
        alphamissense, clinvar, cosmic, dbnsfp, dbsnp, gnomad, mitomap, onekg, primateai, revel,
        scores, spliceai, topmed,
    };

    let (chrom_list, chrom_map) = standard_chrom_map(assembly);
    let out_path = Path::new(output).with_extension("osa2");

    match source {
        // gnomAD has a dedicated streaming encoder that avoids the v1 JSON
        // round-trip.
        "gnomad" => {
            eprintln!("Building gnomad .osa2: {} -> {}", input, out_path.display());
            let (buf_reader, meter) = open_sa_reader_with_meter(input, show_progress)?;
            let records = gnomad::iter_gnomad_osa2(buf_reader, &chrom_map);
            finish_osa2_build(
                &out_path,
                &gnomad::gnomad_osa2_metadata(assembly),
                gnomad::gnomad_osa2_fields(),
                &[],
                records,
                meter,
                input,
                assembly,
            )
        }
        // 1000 Genomes: the v1 parser buffers + sorts; bridge each record to
        // the v2 layout and stream it out.
        "onekg" | "1000g" => {
            eprintln!("Building onekg .osa2: {} -> {}", input, out_path.display());
            let (buf_reader, meter) = open_sa_reader_with_meter(input, show_progress)?;
            let fields = onekg::onekg_osa2_fields();
            // The bridge borrows the schema to re-encode each record; the
            // writer takes ownership of it. Keep separate copies so both can
            // live through the streaming loop.
            let bridge_fields = fields.clone();
            let mut v1 = onekg::parse_onekg_vcf(buf_reader, &chrom_map)?;
            v1.sort_by(|a, b| {
                a.chrom_idx
                    .cmp(&b.chrom_idx)
                    .then(a.position.cmp(&b.position))
            });
            let records = bridge_v1_records(v1.into_iter().map(Ok), &chrom_list, &bridge_fields);
            finish_osa2_build(
                &out_path,
                &onekg::onekg_osa2_metadata(assembly),
                fields,
                &[],
                records,
                meter,
                input,
                assembly,
            )
        }
        // TOPMed has a streaming v1 iterator (the full freeze is ~450M
        // records); bridge it record-by-record with bounded memory.
        "topmed" => {
            eprintln!("Building topmed .osa2: {} -> {}", input, out_path.display());
            let (buf_reader, meter) = open_sa_reader_with_meter(input, show_progress)?;
            let fields = topmed::topmed_osa2_fields();
            let bridge_fields = fields.clone();
            let records = bridge_v1_records(
                topmed::iter_topmed_vcf(buf_reader, &chrom_map),
                &chrom_list,
                &bridge_fields,
            );
            finish_osa2_build(
                &out_path,
                &topmed::topmed_osa2_metadata(assembly),
                fields,
                &[],
                records,
                meter,
                input,
                assembly,
            )
        }
        // AlphaMissense: numeric score + a 3-level categorical class. Streams
        // straight into v2 with a dedicated encoder (the generic v1 bridge
        // cannot encode categoricals) and a fixed 3-entry string table.
        "alphamissense" => {
            eprintln!(
                "Building alphamissense .osa2: {} -> {}",
                input,
                out_path.display()
            );
            let (buf_reader, meter) = open_sa_reader_with_meter(input, show_progress)?;
            let records =
                alphamissense::iter_alphamissense_osa2(buf_reader, &chrom_map, &chrom_list);
            finish_osa2_build(
                &out_path,
                &alphamissense::alphamissense_osa2_metadata(assembly),
                alphamissense::alphamissense_osa2_fields(),
                &alphamissense::alphamissense_string_tables(),
                records,
                meter,
                input,
                assembly,
            )
        }
        // SpliceAI: eight numeric columns (four delta scores, four delta
        // positions) plus a categorical gene symbol. The densest source fastVEP
        // ships against, so the column layout matters most here; the gene
        // vocabulary is only known once the input has streamed past, hence the
        // deferred string table.
        "spliceai" => {
            eprintln!(
                "Building spliceai .osa2: {} -> {}",
                input,
                out_path.display()
            );
            let (buf_reader, meter) = open_sa_reader_with_meter(input, show_progress)?;
            let genes = spliceai::GeneInterner::new();
            let records =
                spliceai::iter_spliceai_osa2(buf_reader, &chrom_map, &chrom_list, genes.clone());
            finish_osa2_build_deferred(
                &out_path,
                &spliceai::spliceai_osa2_metadata(assembly),
                spliceai::spliceai_osa2_fields(),
                records,
                meter,
                input,
                assembly,
                move || genes.string_tables(),
            )
        }
        // MitoMap: whole-record JSON blob (free-text disease + review status).
        // A few thousand records; wired to v2 for format uniformity rather than
        // speed. See `mitomap_osa2_metadata`.
        "mitomap" => {
            eprintln!(
                "Building mitomap .osa2: {} -> {}",
                input,
                out_path.display()
            );
            let (buf_reader, meter) = open_sa_reader_with_meter(input, show_progress)?;
            let v1 = mitomap::parse_mitomap(buf_reader, &chrom_map)?;
            let records = bridge_v1_raw_blobs(v1.into_iter().map(Ok), &chrom_list);
            finish_osa2_build(
                &out_path,
                &mitomap::mitomap_osa2_metadata(assembly),
                fastvep_sa::writer_v2::raw_json_blob_fields(),
                &[],
                records,
                meter,
                input,
                assembly,
            )
        }
        // Custom user VCF: whole-record JSON blob, since the INFO schema is
        // whatever the file carries. Keyed by the user's `--name`.
        "custom_vcf" => {
            let resolved_name = resolve_custom_name(name, input);
            eprintln!(
                "Building custom_vcf .osa2: {} -> {} (name={}, info_fields={})",
                input,
                out_path.display(),
                resolved_name,
                if info_fields.is_empty() {
                    "<all>".to_string()
                } else {
                    info_fields.join(",")
                }
            );
            let (buf_reader, meter) = open_sa_reader_with_meter(input, show_progress)?;
            let mut v1 = fastvep_sa::custom::parse_custom_vcf(
                buf_reader,
                &chrom_map,
                &resolved_name,
                info_fields,
            )?;
            // `parse_custom_vcf` already returns (chrom_idx, position) order, but
            // the v2 streaming writer *hard-fails* on a reopened chunk rather
            // than degrading, so re-establish the ordering at the boundary
            // instead of depending on a cross-crate invariant. A no-op on
            // already-sorted input; same belt-and-braces as the onekg arm above.
            v1.sort_by(|a, b| {
                a.chrom_idx
                    .cmp(&b.chrom_idx)
                    .then(a.position.cmp(&b.position))
            });
            let records = bridge_v1_raw_blobs(v1.into_iter().map(Ok), &chrom_list);
            finish_osa2_build(
                &out_path,
                &fastvep_sa::custom::custom_vcf_osa2_metadata(&resolved_name, assembly),
                fastvep_sa::writer_v2::raw_json_blob_fields(),
                &[],
                records,
                meter,
                input,
                assembly,
            )
        }
        // dbSNP: whole-record JSON blob per variant (RS ID + optional MAF).
        // Streams straight through (the full NCBI release is ~800M records).
        "dbsnp" => {
            eprintln!("Building dbsnp .osa2: {} -> {}", input, out_path.display());
            let (buf_reader, meter) = open_sa_reader_with_meter(input, show_progress)?;
            let records =
                bridge_v1_raw_blobs(dbsnp::iter_dbsnp_vcf(buf_reader, &chrom_map), &chrom_list);
            finish_osa2_build(
                &out_path,
                &dbsnp::dbsnp_osa2_metadata(assembly),
                fastvep_sa::writer_v2::raw_json_blob_fields(),
                &[],
                records,
                meter,
                input,
                assembly,
            )
        }
        // COSMIC: whole-record JSON blob per variant (COSV id + gene + count).
        // The coding-mutations file fits in memory; the v1 parser buffers and
        // sorts, then we bridge each record to a blob and stream it out.
        "cosmic" => {
            eprintln!("Building cosmic .osa2: {} -> {}", input, out_path.display());
            let (buf_reader, meter) = open_sa_reader_with_meter(input, show_progress)?;
            let v1 = cosmic::parse_cosmic_vcf(buf_reader, &chrom_map)?;
            let records = bridge_v1_raw_blobs(v1.into_iter().map(Ok), &chrom_list);
            finish_osa2_build(
                &out_path,
                &cosmic::cosmic_osa2_metadata(assembly),
                fastvep_sa::writer_v2::raw_json_blob_fields(),
                &[],
                records,
                meter,
                input,
                assembly,
            )
        }
        // ClinVar: whole-record JSON blob per variant (significance/phenotype
        // arrays, review status, population AFs). is_array=true is preserved in
        // the v2 metadata to match v1 exactly. Buffered+sorted by the v1 parser.
        "clinvar" => {
            eprintln!(
                "Building clinvar .osa2: {} -> {}",
                input,
                out_path.display()
            );
            let (buf_reader, meter) = open_sa_reader_with_meter(input, show_progress)?;
            let v1 = clinvar::parse_clinvar_vcf(buf_reader, &chrom_map)?;
            let records = bridge_v1_raw_blobs(v1.into_iter().map(Ok), &chrom_list);
            finish_osa2_build(
                &out_path,
                &clinvar::clinvar_osa2_metadata(assembly),
                fastvep_sa::writer_v2::raw_json_blob_fields(),
                &[],
                records,
                meter,
                input,
                assembly,
            )
        }
        // REVEL: single `{"score":..}` object per allele, stored as a
        // whole-record blob (its fixed-decimal score text rides through
        // untouched). The v1 parser buffers+sorts the CSV; column 2 is the
        // GRCh38 position, matching the v1 build path.
        "revel" => {
            eprintln!("Building revel .osa2: {} -> {}", input, out_path.display());
            let (buf_reader, meter) = open_sa_reader_with_meter(input, show_progress)?;
            let v1 = revel::parse_revel(buf_reader, &chrom_map, 2)?;
            let records = bridge_v1_raw_blobs(v1.into_iter().map(Ok), &chrom_list);
            finish_osa2_build(
                &out_path,
                &revel::revel_osa2_metadata(assembly),
                fastvep_sa::writer_v2::raw_json_blob_fields(),
                &[],
                records,
                meter,
                input,
                assembly,
            )
        }
        // Positional per-base scores (PhyloP/GERP/DANN): allele-less, keyed by
        // coordinate. Stream the score TSV/wig and store each bare-number score
        // as a whole-record blob. Genome-wide (~3B positions for PhyloP), so
        // these must stream — never buffer.
        "phylop" => {
            eprintln!("Building phylop .osa2: {} -> {}", input, out_path.display());
            let (buf_reader, meter) = open_sa_reader_with_meter(input, show_progress)?;
            let records =
                bridge_v1_raw_blobs(iter_phylop_auto(buf_reader, &chrom_map), &chrom_list);
            finish_osa2_build(
                &out_path,
                &scores::score_osa2_metadata("phylop", assembly),
                fastvep_sa::writer_v2::raw_json_blob_fields(),
                &[],
                records,
                meter,
                input,
                assembly,
            )
        }
        "gerp" | "dann" => {
            eprintln!(
                "Building {source} .osa2: {} -> {}",
                input,
                out_path.display()
            );
            let (buf_reader, meter) = open_sa_reader_with_meter(input, show_progress)?;
            let records = bridge_v1_raw_blobs(
                scores::iter_score_tsv(buf_reader, &chrom_map, false),
                &chrom_list,
            );
            finish_osa2_build(
                &out_path,
                &scores::score_osa2_metadata(source, assembly),
                fastvep_sa::writer_v2::raw_json_blob_fields(),
                &[],
                records,
                meter,
                input,
                assembly,
            )
        }
        // dbNSFP: composite SIFT/PolyPhen prediction strings per allele,
        // whole-record blob. The v1 parser auto-detects columns from the
        // header, buffers, and sorts; then we bridge to blobs.
        "dbnsfp" => {
            eprintln!("Building dbnsfp .osa2: {} -> {}", input, out_path.display());
            let (buf_reader, meter) = open_sa_reader_with_meter(input, show_progress)?;
            let v1 = dbnsfp::parse_dbnsfp(buf_reader, &chrom_map)?;
            let records = bridge_v1_raw_blobs(v1.into_iter().map(Ok), &chrom_list);
            finish_osa2_build(
                &out_path,
                &dbnsfp::dbnsfp_osa2_metadata(assembly),
                fastvep_sa::writer_v2::raw_json_blob_fields(),
                &[],
                records,
                meter,
                input,
                assembly,
            )
        }
        // PrimateAI: single `{"score":..}` object per allele, whole-record blob.
        "primateai" => {
            eprintln!(
                "Building primateai .osa2: {} -> {}",
                input,
                out_path.display()
            );
            let (buf_reader, meter) = open_sa_reader_with_meter(input, show_progress)?;
            let v1 = primateai::parse_primateai(buf_reader, &chrom_map)?;
            let records = bridge_v1_raw_blobs(v1.into_iter().map(Ok), &chrom_list);
            finish_osa2_build(
                &out_path,
                &primateai::primateai_osa2_metadata(assembly),
                fastvep_sa::writer_v2::raw_json_blob_fields(),
                &[],
                records,
                meter,
                input,
                assembly,
            )
        }
        _ => anyhow::bail!(
            "--format osa2 is currently supported for --source {} (got '{}'). Other \
             sources build v1 .osa; use --format auto (the default) to pick the best \
             format per source, or --format osa to force v1.",
            OSA2_SUPPORTED_SOURCES.join(", "),
            source
        ),
    }
}

/// Adapt an iterator of v1 `AnnotationRecord`s into whole-record-blob
/// `Osa2Record`s: each record's entire JSON is stored as the blob, so the v2
/// reader returns exactly the v1 bytes. Used for string/array payloads
/// (dbSNP, COSMIC, ClinVar) that don't decompose into numeric u32 columns.
fn bridge_v1_raw_blobs<'a, I>(
    records: I,
    chrom_list: &'a [String],
) -> impl Iterator<Item = Result<fastvep_sa::writer_v2::Osa2Record>> + 'a
where
    I: Iterator<Item = Result<fastvep_sa::common::AnnotationRecord>> + 'a,
{
    records.map(move |r| {
        let rec = r?;
        let chrom = chrom_list
            .get(rec.chrom_idx as usize)
            .cloned()
            .ok_or_else(|| anyhow::anyhow!("chrom_idx {} out of range", rec.chrom_idx))?;
        Ok(fastvep_sa::writer_v2::osa2_raw_blob_from_v1(&rec, chrom))
    })
}

/// Adapt an iterator of v1 `AnnotationRecord`s into `Osa2Record`s by
/// re-encoding each record's flat-object JSON against `fields`. Reuses the
/// existing, well-tested v1 source parsers as the front end.
fn bridge_v1_records<'a, I>(
    records: I,
    chrom_list: &'a [String],
    fields: &'a [fastvep_sa::fields::Field],
) -> impl Iterator<Item = Result<fastvep_sa::writer_v2::Osa2Record>> + 'a
where
    I: Iterator<Item = Result<fastvep_sa::common::AnnotationRecord>> + 'a,
{
    records.map(move |r| {
        let rec = r?;
        let chrom = chrom_list
            .get(rec.chrom_idx as usize)
            .cloned()
            .ok_or_else(|| anyhow::anyhow!("chrom_idx {} out of range", rec.chrom_idx))?;
        fastvep_sa::writer_v2::osa2_record_from_v1(&rec, chrom, fields)
    })
}

/// Open a (possibly gzipped) SA input and pair it with a byte-tracking
/// progress meter, matching the v1 streaming builder's setup.
fn open_sa_reader_with_meter(
    input: &str,
    show_progress: bool,
) -> Result<(
    io::BufReader<Box<dyn io::Read>>,
    crate::progress::ProgressMeter,
)> {
    use std::sync::atomic::AtomicU64;
    use std::sync::Arc;

    let file_size = std::fs::metadata(input).map(|m| m.len()).unwrap_or(0);
    let byte_counter = Arc::new(AtomicU64::new(0));
    let reader = open_sa_input(input, Some(Arc::clone(&byte_counter)))?;
    let buf_reader = io::BufReader::new(reader);
    let meter = if show_progress && file_size > 0 {
        crate::progress::ProgressMeter::with_progress(true, file_size, byte_counter)
    } else {
        crate::progress::ProgressMeter::new(show_progress)
    };
    Ok((buf_reader, meter))
}

/// Stream `records` into `out_path` as a `.osa2` file, then report. Mirrors
/// the v1 streaming builder's crash-safety contract: on any failure the
/// partial `.osa2` is removed so no corrupt database is left behind.
// Each argument is an independent coordinate, buffer or flag with no
// natural grouping; bundling them would only move the list to the
// call site.
#[allow(clippy::too_many_arguments)]
fn finish_osa2_build(
    out_path: &Path,
    metadata: &fastvep_sa::writer_v2::Osa2Metadata,
    fields: Vec<fastvep_sa::fields::Field>,
    string_tables: &[(usize, Vec<String>)],
    records: impl Iterator<Item = Result<fastvep_sa::writer_v2::Osa2Record>>,
    meter: crate::progress::ProgressMeter,
    input: &str,
    assembly: &str,
) -> Result<()> {
    let tables = string_tables.to_vec();
    finish_osa2_build_deferred(
        out_path,
        metadata,
        fields,
        records,
        meter,
        input,
        assembly,
        move || tables,
    )
}

/// As [`finish_osa2_build`], but the categorical string tables are produced by
/// `string_tables` *after* every record has streamed past.
///
/// Needed by sources whose categorical vocabulary is only discovered while
/// reading (SpliceAI's gene symbols), as opposed to a fixed table known up front
/// (AlphaMissense's three classes). The `.osa2` string table is global to the
/// archive and written during `finish`, so collecting it during the stream and
/// handing it over at the end is exactly the ordering the format wants.
#[allow(clippy::too_many_arguments)]
fn finish_osa2_build_deferred<F>(
    out_path: &Path,
    metadata: &fastvep_sa::writer_v2::Osa2Metadata,
    fields: Vec<fastvep_sa::fields::Field>,
    records: impl Iterator<Item = Result<fastvep_sa::writer_v2::Osa2Record>>,
    mut meter: crate::progress::ProgressMeter,
    input: &str,
    assembly: &str,
    string_tables: F,
) -> Result<()>
where
    F: FnOnce() -> Vec<(usize, Vec<String>)>,
{
    let n = match write_osa2_stream(
        out_path,
        metadata,
        fields,
        records,
        &mut meter,
        string_tables,
    ) {
        Ok(n) => n,
        Err(e) => {
            let _ = std::fs::remove_file(out_path);
            return Err(e);
        }
    };
    meter.finish();
    if n == 0 {
        eprintln!(
            "Warning: 0 records parsed from {} — if the input is non-empty, none of its \
             chromosome names matched assembly '{}'.",
            input, assembly
        );
    }
    eprintln!(
        "Wrote: {} ({} records)",
        out_path.display(),
        crate::progress::fmt_count(n)
    );
    Ok(())
}

/// Core streaming loop: create the `.osa2` writer, push every record, finish.
/// Returns the number of records written.
fn write_osa2_stream<F>(
    out_path: &Path,
    metadata: &fastvep_sa::writer_v2::Osa2Metadata,
    fields: Vec<fastvep_sa::fields::Field>,
    records: impl Iterator<Item = Result<fastvep_sa::writer_v2::Osa2Record>>,
    meter: &mut crate::progress::ProgressMeter,
    string_tables: F,
) -> Result<u64>
where
    F: FnOnce() -> Vec<(usize, Vec<String>)>,
{
    use fastvep_sa::writer_v2::Osa2StreamWriter;

    let out_file = std::fs::File::create(out_path)
        .with_context(|| format!("Creating {}", out_path.display()))?;
    let mut writer = Osa2StreamWriter::new(BufWriter::new(out_file), metadata, fields)?;

    let mut n = 0u64;
    for record in records {
        let record = record?;
        meter.update();
        writer.push(record)?;
        n += 1;
    }
    // Resolved after the stream so a source that interns its categorical
    // vocabulary while reading hands over the completed table. `finish` writes
    // the tables into the archive, so this still lands before finalization.
    for (field_idx, table) in string_tables() {
        writer.set_string_table(field_idx, table);
    }
    writer.finish()?;
    Ok(n)
}

/// Transcode an existing v1 `.osa` database into a v2 `.osa2`, in place of
/// re-downloading and re-parsing the upstream source.
///
/// Every variant-level source now builds v2 by default, but a `--sa-dir`
/// assembled before that change is full of `.osa` files whose upstream releases
/// may be large (dbSNP, SpliceAI), slow to fetch, or no longer available at the
/// exact version that was built. Rebuilding from source is not always an option,
/// so this reads the `.osa` and writes the equivalent `.osa2`.
///
/// **What is preserved:** the record set and every record's JSON payload, byte
/// for byte, plus the database's identity and matching semantics (`json_key`,
/// `name`, `version`, `description`, `assembly`, `match_by_allele`, `is_array`,
/// `is_positional`). Queries against the converted file return exactly what the
/// `.osa` returned.
///
/// **What is not:** the records are stored as whole-record JSON blobs, because
/// that is all a `.osa` retains — it has no field schema to recover a numeric
/// column layout from. A conversion therefore gets v2's chunked index, its
/// mmap-and-inflate read path, and its chunk-level zstd, but *not* the parallel
/// u32 columns that [`OSA2_DECOMPOSED_SOURCES`] build natively. Which encoding
/// ends up smaller is data-dependent (blob columns zstd well when adjacent
/// records are similar; value columns win when the payload is genuinely
/// numeric), so the tool points at the rebuild without promising a size win.
pub fn run_sa_convert(
    input: &str,
    output: &str,
    chunk_bits: u32,
    show_progress: bool,
) -> Result<()> {
    use fastvep_sa::reader::SaReader;
    use fastvep_sa::writer_v2::{Osa2Metadata, Osa2Record};

    let in_path = Path::new(input);
    match in_path.extension().and_then(|e| e.to_str()) {
        Some("osa") => {}
        Some("osa2") => anyhow::bail!(
            "{} is already a v2 .osa2 database; nothing to convert.",
            input
        ),
        Some(other @ ("osi" | "oga")) => anyhow::bail!(
            "sa-convert only converts variant-level .osa databases. '{}' is a .{} file \
             ({}), and v2 is a variant-level container with no equivalent form — keep \
             using it as-is.",
            input,
            other,
            if other == "osi" {
                "interval-level"
            } else {
                "gene-level"
            }
        ),
        _ => anyhow::bail!("sa-convert expects a v1 '.osa' input; got '{}'.", input),
    }

    let out_path = Path::new(output).with_extension("osa2");
    if out_path == in_path {
        anyhow::bail!("--output would overwrite the input ({})", in_path.display());
    }

    let reader = SaReader::open(in_path)
        .with_context(|| format!("Opening v1 database {}", in_path.display()))?;
    let header = reader.header().clone();

    eprintln!(
        "Converting {} -> {} (source={}, key={})",
        in_path.display(),
        out_path.display(),
        header.name,
        header.json_key
    );
    match source_from_json_key(&header.json_key) {
        Some(src) if source_has_decomposed_osa2(src) => eprintln!(
            "note: '{src}' has a column-oriented v2 encoder. This conversion preserves every \
             payload exactly, but can only store them as JSON blobs — a .osa retains no field \
             schema to rebuild the columns from. If you still have the upstream release, \
             `sa-build --source {src} --format osa2` gives you the column encoding instead; \
             which of the two is smaller depends on the data, so compare if size matters."
        ),
        Some(src) => eprintln!(
            "note: '{src}' stores whole-record JSON blobs in v2 as well, so this conversion is \
             equivalent to rebuilding from the upstream release."
        ),
        None => {}
    }

    let metadata = Osa2Metadata {
        format_version: 2,
        name: header.name.clone(),
        version: header.version.clone(),
        assembly: header.assembly.clone(),
        json_key: header.json_key.clone(),
        match_by_allele: header.match_by_allele,
        is_array: header.is_array,
        is_positional: header.is_positional,
        chunk_bits,
        description: header.description.clone(),
    };

    // The `.osa` is the progress denominator: `iter_records` walks it once, so
    // bytes-of-input is a faithful measure of how far along the conversion is.
    let total = std::fs::metadata(in_path).map(|m| m.len()).unwrap_or(0);
    let mut meter = if show_progress && total > 0 {
        crate::progress::ProgressMeter::new(true)
    } else {
        crate::progress::ProgressMeter::new(false)
    };

    let records = reader.iter_records().map(|r| {
        let (chrom, entry) = r?;
        Ok(Osa2Record {
            chrom: chrom.to_string(),
            position: entry.position,
            ref_allele: entry.ref_allele.into_bytes(),
            alt_allele: entry.alt_allele.into_bytes(),
            values: Vec::new(),
            json_blob: Some(entry.json),
        })
    });

    // Same crash-safety contract as the builders: a partial `.osa2` is removed
    // so a failed conversion never leaves a corrupt database behind.
    let n = match write_osa2_stream(
        &out_path,
        &metadata,
        fastvep_sa::writer_v2::raw_json_blob_fields(),
        records,
        &mut meter,
        Vec::new,
    ) {
        Ok(n) => n,
        Err(e) => {
            let _ = std::fs::remove_file(&out_path);
            return Err(e);
        }
    };
    meter.finish();

    if n == 0 {
        eprintln!(
            "warning: {} contained 0 records; wrote an empty .osa2.",
            in_path.display()
        );
    }
    let in_bytes = std::fs::metadata(in_path).map(|m| m.len()).unwrap_or(0)
        + std::fs::metadata(in_path.with_extension("osa.idx"))
            .map(|m| m.len())
            .unwrap_or(0);
    let out_bytes = std::fs::metadata(&out_path).map(|m| m.len()).unwrap_or(0);
    eprintln!(
        "Wrote: {} ({} records, {} -> {})",
        out_path.display(),
        crate::progress::fmt_count(n),
        crate::progress::fmt_bytes(in_bytes),
        crate::progress::fmt_bytes(out_bytes)
    );
    Ok(())
}
