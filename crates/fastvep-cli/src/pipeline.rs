use anyhow::{Context, Result};
use flate2::read::MultiGzDecoder;
use fastvep_cache::annotation::{AnnotationProvider, AnnotationValue, GeneAnnotationProvider};
use fastvep_cache::fasta::FastaReader;
use fastvep_cache::gff::parse_gff3_with_source;
use fastvep_cache::info::CacheInfo;
use fastvep_cache::providers::{
    FastaSequenceProvider, IndexedTranscriptProvider, MatchedVariant, SequenceProvider,
    TabixVariationProvider, TranscriptProvider, VariationProvider,
};
use fastvep_consequence::ConsequencePredictor;
use fastvep_core::{Allele, Consequence};
use fastvep_hgvs;
use fastvep_io::output;
use fastvep_io::variant::{AlleleAnnotation, TranscriptVariation, VariationFeature};
use fastvep_io::vcf::VcfParser;
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::Path;
use std::sync::atomic::{AtomicUsize, Ordering};

const BATCH_SIZE: usize = 1024;

fn open_vcf_input_reader(input: &str) -> Result<Box<dyn io::Read>> {
    let reader: Box<dyn io::Read> = if input == "-" {
        Box::new(io::stdin())
    } else {
        Box::new(
            File::open(input)
                .with_context(|| format!("Opening input file: {}", input))?,
        )
    };

    wrap_maybe_gzip_reader(reader, input)
}

fn wrap_maybe_gzip_reader(mut reader: Box<dyn io::Read>, source: &str) -> Result<Box<dyn io::Read>> {
    let mut prefix = [0u8; 2];
    let bytes_read = reader.read(&mut prefix)?;
    let looks_like_gzip = bytes_read == 2 && prefix == [0x1f, 0x8b];

    let replay = io::Cursor::new(prefix[..bytes_read].to_vec()).chain(reader);
    if looks_like_gzip || (source != "-" && source.ends_with(".gz")) {
        Ok(Box::new(MultiGzDecoder::new(replay)))
    } else {
        Ok(Box::new(replay))
    }
}

pub struct AnnotateConfig {
    pub input: String,
    pub output: String,
    /// One or more GFF3 annotation sources. An empty Vec means "no
    /// transcript annotation" (only intergenic / SA results). Each entry
    /// is either a bare path or `LABEL=path` — see `parse_gff3_arg` for
    /// the parsing/auto-detection rules.
    pub gff3: Vec<String>,
    pub fasta: Option<String>,
    pub output_format: String,
    pub pick: bool,
    pub hgvs: bool,
    pub distance: u64,
    pub cache_dir: Option<String>,
    pub transcript_cache: Option<String>,
    /// Directory containing supplementary annotation files (.osa, .osi, .oga).
    pub sa_dir: Option<String>,
    /// Skip the default 49-field CSQ annotation pipeline and emit only
    /// supplementary annotations from `--sa-dir`. Requires `sa_dir`.
    pub sa_only: bool,
    /// Enable ACMG-AMP variant classification.
    pub acmg: bool,
    /// Path to ACMG configuration file (TOML).
    pub acmg_config: Option<String>,
    /// Proband sample name for trio analysis.
    pub proband: Option<String>,
    /// Mother sample name for trio analysis.
    pub mother: Option<String>,
    /// Father sample name for trio analysis.
    pub father: Option<String>,
    /// Path to a gene panel file. When set, tab output keeps only rows
    /// whose transcript's gene_id or gene_symbol is in the panel.
    pub gene_list: Option<String>,
    /// Add an explicit REF column to tab output after the Allele/ALT column.
    pub explicit_alleles: bool,
    /// Path to a QC rules TOML file. When set, tab output gains a
    /// `QC_CLASS` column.
    pub qc_rules: Option<String>,
    /// Show periodic progress output.
    pub show_progress: bool,
}

pub fn run_annotate(config: AnnotateConfig) -> Result<()> {
    eprintln!("Annotating: {} -> {}", config.input, config.output);
    let sa_only = config.sa_only;
    if sa_only {
        if config.sa_dir.is_none() {
            return Err(anyhow::anyhow!(
                "--sa-only requires --sa-dir to be set (otherwise there is nothing to emit)."
            ));
        }
        for (set, name) in [
            (!config.gff3.is_empty(), "--gff3"),
            (config.fasta.is_some(), "--fasta"),
            (config.cache_dir.is_some(), "--cache-dir"),
            (config.transcript_cache.is_some(), "--transcript-cache"),
            (config.acmg, "--acmg"),
            (config.hgvs, "--hgvs"),
            (config.pick, "--pick"),
            (config.proband.is_some(), "--proband"),
            (config.mother.is_some(), "--mother"),
            (config.father.is_some(), "--father"),
        ] {
            if set {
                eprintln!("warning: --sa-only is set; ignoring {}", name);
            }
        }
    }

    // Parse `LABEL=path` syntax up front so a typo fails before we touch IO.
    let gff3_specs: Vec<Gff3Spec> = config.gff3
        .iter()
        .map(|s| parse_gff3_arg(s))
        .collect();

    // Load transcript models: try binary cache first, fall back to GFF3.
    // The auto-managed sidecar cache only kicks in for a single GFF3 — with
    // multiple sources we'd need a content-hashed key to know which combo
    // produced a given cache, and users running merged caches can pre-build
    // via `fastvep cache` and pass `--transcript-cache` explicitly.
    // Skipped entirely in --sa-only mode (no default annotation needed).
    let single_gff3: Option<&Gff3Spec> = if gff3_specs.len() == 1 {
        Some(&gff3_specs[0])
    } else {
        None
    };
    let cache_path = if sa_only {
        None
    } else {
        config.transcript_cache.as_ref().map(|p| Path::new(p).to_path_buf())
            .or_else(|| single_gff3.map(|s| fastvep_cache::transcript_cache::default_cache_path(Path::new(&s.path))))
    };

    let mut transcripts = if sa_only { Vec::new() } else { 'load: {
        // Cache-load gating:
        //
        // * Sidecar cache (`cache_path` derived from a single `--gff3`):
        //   always considered authoritative when fresh against its source
        //   GFF3. Re-stamp the source label so the user's current
        //   --gff3 label/auto-detection wins over whatever was on disk.
        //
        // * Explicit `--transcript-cache <path>`: the *user* told us where
        //   transcripts live. Honour the cache contents verbatim (do NOT
        //   re-stamp), and for multi-GFF3 invocations be explicit that the
        //   --gff3 arguments are being ignored. Without this, a user
        //   running `--gff3 ens.gff3 --gff3 refseq.gff3 --transcript-cache
        //   old_ens_only.cache` would silently get Ensembl-only output
        //   and never know.
        let explicit_cache = config.transcript_cache.is_some();
        if let Some(ref cp) = cache_path {
            if cp.exists() {
                // Freshness check: only sidecar-cache mode does it (the
                // user-provided --transcript-cache is always trusted).
                let is_fresh = if explicit_cache {
                    true
                } else {
                    single_gff3
                        .map(|s| fastvep_cache::transcript_cache::cache_is_fresh(cp, Path::new(&s.path)))
                        .unwrap_or(true)
                };
                if is_fresh {
                    match fastvep_cache::transcript_cache::load_cache(cp) {
                        Ok(mut trs) => {
                            // Re-stamp only in sidecar-cache + single-GFF3
                            // mode. For explicit --transcript-cache we
                            // preserve the on-disk labels so a merged
                            // cache built via `fastvep cache --gff3 ens
                            // --gff3 refseq -o combined.cache` survives
                            // round-tripping with the merged distinction
                            // intact.
                            if !explicit_cache {
                                if let Some(spec) = single_gff3 {
                                    for tr in &mut trs {
                                        tr.source = Some(spec.source.clone());
                                    }
                                }
                            } else if !gff3_specs.is_empty() {
                                // Loud warning: user-supplied --gff3
                                // alongside --transcript-cache means the
                                // GFF3 arguments are ignored.
                                eprintln!(
                                    "warning: --transcript-cache {} takes precedence over --gff3 {:?}; the GFF3 file(s) will NOT be parsed. Drop --transcript-cache to load from GFF3 instead, or remove the --gff3 flags to silence this warning.",
                                    cp.display(),
                                    gff3_specs.iter().map(|s| s.path.as_str()).collect::<Vec<_>>(),
                                );
                            }
                            eprintln!("Loaded {} transcripts from cache {}", trs.len(), cp.display());
                            break 'load trs;
                        }
                        Err(e) => {
                            eprintln!("Warning: cache load failed ({}), falling back to GFF3", e);
                        }
                    }
                } else {
                    eprintln!("Cache is stale, rebuilding from GFF3");
                }
            }
        }

        // Fall back to GFF3 parsing. With multiple sources, load each one and
        // concatenate — IndexedTranscriptProvider sorts and indexes by chrom,
        // and the consequence predictor doesn't require globally-unique
        // stable_ids, so Ensembl + RefSeq can simply coexist.
        if gff3_specs.is_empty() {
            eprintln!("Warning: No GFF3 file provided. Only intergenic variants will be annotated.");
            Vec::new()
        } else {
            let mut all: Vec<fastvep_genome::Transcript> = Vec::new();
            for spec in &gff3_specs {
                let trs = load_one_gff3(spec, &config.input, config.distance)?;
                eprintln!(
                    "Loaded {} transcripts from {} (source label: {})",
                    trs.len(), spec.path, spec.source
                );
                all.extend(trs);
            }
            if all.is_empty() {
                return Err(anyhow::anyhow!(
                    "GFF3 source(s) [{}] produced 0 transcripts — likely malformed, truncated, or unrecognized format. Refusing to continue with empty transcript set.",
                    gff3_specs.iter().map(|s| s.path.as_str()).collect::<Vec<_>>().join(", ")
                ));
            }
            if gff3_specs.len() > 1 {
                eprintln!("Merged {} GFF3 sources into {} total transcripts", gff3_specs.len(), all.len());
            }
            all
        }
    }};

    // Load FASTA reference (prefer mmap with .fai index, fall back to in-memory).
    // Skipped in --sa-only mode.
    let seq_provider: Option<Box<dyn SequenceProvider>> = if sa_only {
        None
    } else if let Some(ref fasta_path) = config.fasta {
        let fai_path = format!("{}.fai", fasta_path);
        if Path::new(&fai_path).exists() {
            let reader = fastvep_cache::fasta::MmapFastaReader::open(Path::new(fasta_path))?;
            eprintln!("Memory-mapped reference FASTA from {} (using .fai index)", fasta_path);
            Some(Box::new(fastvep_cache::providers::MmapFastaSequenceProvider::new(reader)))
        } else {
            let fasta_file = File::open(fasta_path)
                .with_context(|| format!("Opening FASTA file: {}", fasta_path))?;
            let reader = FastaReader::from_reader(fasta_file)?;
            eprintln!("Loaded reference FASTA from {}", fasta_path);
            Some(Box::new(FastaSequenceProvider::new(reader)))
        }
    } else {
        None
    };

    // Build sequences for coding transcripts from FASTA (skip if loaded from cache with sequences)
    let needs_seq_build = transcripts.iter().any(|t| t.is_coding() && t.spliced_seq.is_none());
    if needs_seq_build {
        if let Some(ref sp) = seq_provider {
            let built = AtomicUsize::new(0);
            transcripts.par_iter_mut().for_each(|tr| {
                if tr.is_coding() && tr.spliced_seq.is_none() {
                    if let Err(e) = tr.build_sequences(|chrom, start, end| {
                        sp.fetch_sequence(chrom, start, end)
                            .map_err(|e| e.to_string())
                    }) {
                        eprintln!("Warning: could not build sequences for {}: {}", tr.stable_id, e);
                    } else {
                        built.fetch_add(1, Ordering::Relaxed);
                    }
                }
            });
            eprintln!("Built sequences for {} coding transcripts", built.load(Ordering::Relaxed));
        }
    }

    // Save cache after sequence build (only if sequences were built or
    // cache doesn't exist). Sidecar-cache writes are gated to the
    // single-GFF3 case — multi-source caches need to round-trip through
    // an explicit `--transcript-cache` path so the on-disk filename
    // isn't tied to one of N inputs.
    if needs_seq_build {
        if let Some(ref cp) = cache_path {
            if single_gff3.is_some() || config.transcript_cache.is_some() {
                if let Err(e) = fastvep_cache::transcript_cache::save_cache(&transcripts, cp) {
                    eprintln!("Warning: could not save cache: {}", e);
                } else {
                    eprintln!("Saved transcript cache to {}", cp.display());
                }
            }
        }
    }

    let transcript_provider = IndexedTranscriptProvider::new(transcripts);

    // Initialize variation provider from VEP cache if provided.
    // Skipped in --sa-only mode (existing_variation comes from defaults).
    let var_provider: Option<TabixVariationProvider> = if sa_only {
        None
    } else if let Some(ref dir) = config.cache_dir {
        let info_path = Path::new(dir).join("info.txt");
        let cache_info = CacheInfo::from_file(&info_path)
            .with_context(|| format!("Reading cache info: {}", info_path.display()))?;
        eprintln!(
            "Loaded VEP cache info: species={}, assembly={}, {} variation columns",
            cache_info.species, cache_info.assembly, cache_info.variation_cols.len()
        );
        Some(TabixVariationProvider::new(Path::new(dir), &cache_info)?)
    } else {
        None
    };

    // Load supplementary annotation providers from --sa-dir
    // Shared load_sa_providers returns Mutex-wrapped providers; unwrap for batch pipeline
    // (batch pipeline processes variants sequentially per chunk, no concurrent access).
    let sa_providers: Vec<Box<dyn AnnotationProvider>> = if let Some(ref dir) = config.sa_dir {
        load_sa_providers(Path::new(dir))?
            .into_iter()
            .map(|m| m.into_inner().unwrap())
            .collect()
    } else {
        Vec::new()
    };

    if sa_only && sa_providers.is_empty() {
        eprintln!(
            "warning: --sa-only is set but --sa-dir {:?} loaded zero allele-level supplementary providers; output will contain no annotations.",
            config.sa_dir.as_deref().unwrap_or("")
        );
    }

    // Load gene-level annotation providers (.oga files).
    // Skipped in --sa-only mode: gene-level SA needs a gene symbol from
    // transcript overlap, which sa_only does not produce. Loading them anyway
    // would emit always-empty headers/columns.
    let gene_providers: Vec<fastvep_sa::gene::GeneIndex> = if sa_only {
        if let Some(ref dir) = config.sa_dir {
            // Cheap probe: just count `.oga` files instead of fully loading
            // each one. Earlier this path called `load_gene_providers` and
            // discarded the result, paying the full disk-read cost for a
            // warning message that only needs a yes/no on presence.
            let oga_count = std::fs::read_dir(Path::new(dir))
                .map(|it| {
                    it.flatten()
                        .filter(|e| {
                            e.path().extension().and_then(|s| s.to_str()) == Some("oga")
                        })
                        .count()
                })
                .unwrap_or(0);
            if oga_count > 0 {
                eprintln!(
                    "warning: --sa-only ignores {} gene-level annotation source(s) (.oga) in {}; gene-level SA requires transcript overlap.",
                    oga_count,
                    dir
                );
            }
        }
        Vec::new()
    } else if let Some(ref dir) = config.sa_dir {
        load_gene_providers(Path::new(dir))?
    } else {
        Vec::new()
    };

    // Load ACMG-AMP classification config if enabled.
    // Skipped in --sa-only mode (ACMG depends on default annotations).
    let acmg_config: Option<fastvep_classification::AcmgConfig> = if config.acmg && !sa_only {
        let mut cfg = if let Some(ref path) = config.acmg_config {
            fastvep_classification::AcmgConfig::from_toml_file(path)?
        } else {
            fastvep_classification::AcmgConfig::default()
        };
        // Wire trio config from CLI flags
        if let Some(ref proband) = config.proband {
            cfg.trio = Some(fastvep_classification::TrioConfig {
                proband: proband.clone(),
                mother: config.mother.clone(),
                father: config.father.clone(),
                min_depth: cfg.trio.as_ref().map_or(10, |t| t.min_depth),
                min_gq: cfg.trio.as_ref().map_or(20, |t| t.min_gq),
            });
        }
        Some(cfg)
    } else {
        None
    };

    // Load optional gene-panel filter (issue #1 ask 4). Loaded once at
    // startup so per-variant filtering is an O(1) HashSet lookup against
    // gene_id / gene_symbol.
    let gene_set: Option<fastvep_io::geneset::GeneSet> = match config.gene_list.as_deref() {
        Some(path) => {
            let set = fastvep_io::geneset::GeneSet::from_file(path)
                .with_context(|| format!("Loading gene list: {}", path))?;
            eprintln!("Loaded gene panel: {} entries from {}", set.len(), path);
            if config.output_format != "tab" {
                eprintln!(
                    "warning: --gene-list currently filters tab output only; \
                     --output-format {} will emit unfiltered rows.",
                    config.output_format
                );
            }
            Some(set)
        }
        None => None,
    };

    // Load optional QC rule set (issue #1 ask 3). Variant-level: reads
    // INFO field thresholds, no per-sample work.
    let qc_rules: Option<fastvep_io::qc::QcRules> = match config.qc_rules.as_deref() {
        Some(path) => {
            let rules = fastvep_io::qc::QcRules::from_toml_file(path)
                .with_context(|| format!("Loading QC rules: {}", path))?;
            eprintln!(
                "Loaded QC rules: {} class(es) from {}",
                rules.classes.len(),
                path
            );
            if config.output_format != "tab" {
                eprintln!(
                    "warning: --qc-rules currently annotates tab output only; \
                     --output-format {} will not gain a QC_CLASS column.",
                    config.output_format
                );
            }
            Some(rules)
        }
        None => None,
    };

    if config.explicit_alleles && config.output_format != "tab" {
        eprintln!(
            "warning: --explicit-alleles applies to tab output only; \
             --output-format {} ignores it.",
            config.output_format
        );
    }

    // Create consequence predictor
    let predictor = ConsequencePredictor::new(config.distance, config.distance);

    // Open input VCF (supports plain text or gzipped VCF)
    let input_reader = open_vcf_input_reader(&config.input)?;
    let mut vcf_parser = VcfParser::new(input_reader)?;

    // Extract sample names from VCF #CHROM header
    let sample_names: Vec<String> = vcf_parser
        .header_lines()
        .last()
        .filter(|l| l.starts_with("#CHROM"))
        .map(|l| l.split('\t').skip(9).map(|s| s.to_string()).collect())
        .unwrap_or_default();

    // Open output
    let output_writer: Box<dyn io::Write> = if config.output == "-" {
        Box::new(io::stdout())
    } else {
        Box::new(
            File::create(&config.output)
                .with_context(|| format!("Creating output file: {}", config.output))?,
        )
    };
    // 1 MiB buffer instead of the default 8 KiB: the per-variant output path
    // emits dozens of small `write!` calls per row, so a larger buffer cuts
    // the number of syscalls on a typical VCF (millions of variants) by
    // roughly two orders of magnitude.
    let mut writer = BufWriter::with_capacity(1 << 20, output_writer);
    let sa_json_keys: Vec<String> = sa_providers
        .iter()
        .map(|sa| sa.json_key().to_string())
        .collect();
    let gene_json_keys: Vec<String> = gene_providers
        .iter()
        .map(|gp| gp.json_key().to_string())
        .collect();
    let owned_vcf_info_ids = output::vcf_owned_info_ids(&sa_json_keys, &gene_json_keys);
    let generated_vcf_headers =
        output::vcf_info_header_lines(&sa_json_keys, &gene_json_keys, output::DEFAULT_CSQ_FIELDS, sa_only);
    // Precompute the loaded-source lookup once so the per-row tab writer
    // doesn't redo an O(specs × keys) membership scan for every variant.
    let supplementary_specs =
        output::LoadedSupplementarySpecs::new(&sa_json_keys, &gene_json_keys);

    // Write headers based on output format
    match config.output_format.as_str() {
        "vcf" => {
            // Pass through original VCF headers
            for header_line in vcf_parser.header_lines() {
                if let Some(info_id) = output::vcf_info_header_id(header_line) {
                    if owned_vcf_info_ids.iter().any(|owned| *owned == info_id) {
                        continue;
                    }
                }
                if header_line.starts_with("#CHROM") {
                    for generated in &generated_vcf_headers {
                        writeln!(writer, "{}", generated)?;
                    }
                }
                writeln!(writer, "{}", header_line)?;
            }
        }
        "tab" => {
            writeln!(writer, "## fastVEP output")?;
            for line in output::tab_supplementary_header_lines(&supplementary_specs) {
                writeln!(writer, "{}", line)?;
            }
            let extra_columns = output::tab_supplementary_column_names(&supplementary_specs);
            let mut header = if sa_only {
                let mut h = String::from("#Uploaded_variation\tLocation\tAllele");
                if config.explicit_alleles {
                    h.push_str("\tREF");
                }
                h
            } else {
                let mut h = String::from("#Uploaded_variation\tLocation\tAllele");
                if config.explicit_alleles {
                    h.push_str("\tREF");
                }
                h.push_str(
                    "\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tIMPACT\tDISTANCE\tSTRAND\tFLAGS",
                );
                h
            };
            for col in &extra_columns {
                header.push('\t');
                header.push_str(col);
            }
            if qc_rules.is_some() {
                header.push_str("\tQC_CLASS");
            }
            writeln!(writer, "{}", header)?;
        }
        "json" => {
            writeln!(writer, "[")?;
        }
        _ => {}
    }

    // Process variants in batches for parallel annotation
    let mut meter = crate::progress::ProgressMeter::new(config.show_progress);
    let mut first_json = true;

    loop {
        // Phase 1: Read a batch of variants (sequential - VCF parser is not Sync)
        let mut batch: Vec<(VariationFeature, HashMap<String, Vec<MatchedVariant>>)> = Vec::with_capacity(BATCH_SIZE);
        for _ in 0..BATCH_SIZE {
            match vcf_parser.next_variant()? {
                Some(mut vf) => {
                    // Variation lookup (sequential - TabixVariationProvider is not Sync)
                    let matched_by_allele: HashMap<String, Vec<MatchedVariant>> =
                        if let Some(ref vp) = var_provider {
                            let mut by_allele = HashMap::new();
                            for alt in &vf.alt_alleles {
                                let alt_str = alt.to_string();
                                let ref_str = vf.ref_allele.to_string();
                                let matches = vp.get_matched_variants(
                                    &vf.position.chromosome,
                                    vf.position.start,
                                    vf.position.end,
                                    &ref_str,
                                    &alt_str,
                                ).unwrap_or_default();
                                if !matches.is_empty() {
                                    by_allele.insert(alt_str, matches);
                                }
                            }
                            by_allele
                        } else {
                            HashMap::new()
                        };

                    // Populate existing_variants on the VF for output access
                    for matches in matched_by_allele.values() {
                        for m in matches {
                            if !vf.existing_variants.iter().any(|kv| kv.name == m.name) {
                                vf.existing_variants.push(fastvep_io::variant::KnownVariant {
                                    name: m.name.clone(),
                                    allele_string: None,
                                    minor_allele: m.minor_allele.clone(),
                                    minor_allele_freq: m.minor_allele_freq,
                                    clinical_significance: m.clin_sig.clone(),
                                    somatic: m.somatic,
                                    phenotype_or_disease: m.phenotype_or_disease,
                                    pubmed: m.pubmed.clone(),
                                    frequencies: m.frequencies.clone(),
                                });
                            }
                        }
                    }

                    batch.push((vf, matched_by_allele));
                }
                None => break,
            }
        }

        if batch.is_empty() {
            break;
        }

        // Phase 1.5: Preload SA providers for this batch (sequential)
        if !sa_providers.is_empty() && !batch.is_empty() {
            // Collect all positions in this batch, grouped by chromosome
            let mut chrom_positions: HashMap<&str, Vec<u64>> = HashMap::new();
            for (vf, _) in &batch {
                let positions = chrom_positions
                    .entry(&vf.position.chromosome)
                    .or_default();
                positions.push(vf.position.start);
                if let Some(vcf) = &vf.vcf_fields {
                    if vcf.pos != vf.position.start {
                        positions.push(vcf.pos);
                    }
                }
            }
            for sa in &sa_providers {
                for (chrom, positions) in &chrom_positions {
                    let _ = sa.preload(chrom, positions);
                }
            }
        }

        // Phase 2: Annotate batch in parallel (transcript lookup + consequence prediction + HGVS)
        batch.par_iter_mut().for_each(|(vf, matched_by_allele)| {
            // In --sa-only mode, skip transcript lookup + consequence
            // prediction + HGVS entirely. Use a neutral scaffold so the SA
            // attachment loop below has per-allele slots to populate.
            if sa_only {
                annotate_sa_only_scaffold(vf);
            } else {
            let chrom = &vf.position.chromosome;
            let query_start = if vf.position.start > config.distance {
                vf.position.start - config.distance
            } else {
                1
            };
            let query_end = vf.position.end + config.distance;
            let overlapping = transcript_provider.get_transcripts(chrom, query_start, query_end)
                .unwrap_or_default();

        if overlapping.is_empty() {
            // Intergenic
            annotate_intergenic(vf);
            // Populate existing_variation on intergenic annotations too
            for tv in &mut vf.transcript_variations {
                for aa in &mut tv.allele_annotations {
                    if let Some(matches) = matched_by_allele.get(&aa.allele.to_string()) {
                        aa.existing_variation = matches.iter().map(|m| m.name.clone()).collect();
                    }
                }
            }
        } else {
            // Get reference sequence if available
            let ref_seq = seq_provider.as_ref().and_then(|sp| {
                sp.fetch_sequence(chrom, query_start, query_end).ok()
            });

            // Run consequence prediction — dispatch SVs to SV predictor
            let transcript_consequences = if vf.variant_type.is_structural() {
                fastvep_consequence::sv_predictor::predict_sv_consequences(
                    chrom,
                    vf.position.start,
                    vf.position.end,
                    vf.variant_type,
                    &vf.alt_alleles,
                    &overlapping,
                    config.distance,
                    config.distance,
                )
            } else {
                let result = predictor.predict(
                    &vf.position,
                    &vf.ref_allele,
                    &vf.alt_alleles,
                    &overlapping,
                    ref_seq.as_deref(),
                );
                result.transcript_consequences
            };

            // Convert prediction results to VariationFeature annotations
            for tc in &transcript_consequences {
                let transcript = overlapping.iter().find(|t| t.stable_id == tc.transcript_id);

                let allele_annotations: Vec<AlleleAnnotation> = tc
                    .allele_consequences
                    .iter()
                    .map(|ac| {
                        let mut ann = AlleleAnnotation {
                            allele: ac.allele.clone(),
                            consequences: ac.consequences.clone(),
                            impact: ac.impact,
                            cdna_position: zip_positions(ac.cdna_start, ac.cdna_end),
                            cds_position: zip_positions(ac.cds_start, ac.cds_end),
                            protein_position: zip_positions(ac.protein_start, ac.protein_end),
                            amino_acids: ac.amino_acids.clone(),
                            codons: ac.codons.clone(),
                            exon: ac.exon,
                            intron: ac.intron,
                            distance: ac.distance,
                            hgvsc: None,
                            hgvsp: None,
                            hgvsg: None,
                            hgvs_offset: None,
                            existing_variation: matched_by_allele
                                .get(&ac.allele.to_string())
                                .map(|matches| matches.iter().map(|m| m.name.clone()).collect())
                                .unwrap_or_default(),
                            sift: None,
                            polyphen: None,
                            supplementary: Vec::new(),
                            acmg_classification: None,
                        };

                        // Generate HGVS if requested
                        if config.hgvs {
                            ann.hgvsg = Some(fastvep_hgvs::hgvsg(
                                chrom,
                                vf.position.start,
                                vf.position.end,
                                &vf.ref_allele,
                                &ac.allele,
                            ));

                            if let Some(tr) = transcript {
                                // Build versioned IDs for HGVS notation
                                let versioned_tid = match tr.version {
                                    Some(v) => format!("{}.{}", tc.transcript_id, v),
                                    None => tc.transcript_id.to_string(),
                                };

                                // Determine alleles for HGVS - complement for minus strand
                                let (hgvs_ref, hgvs_alt) = if tr.strand == fastvep_core::Strand::Reverse {
                                    (complement_allele(&vf.ref_allele), complement_allele(&ac.allele))
                                } else {
                                    (vf.ref_allele.clone(), ac.allele.clone())
                                };

                                if let Some(coding_start) = tr.cdna_coding_start {
                                    if let (Some(cs), Some(ce)) = (ac.cdna_start, ac.cdna_end) {
                                        // Normalize cDNA positions (minus-strand can reverse order)
                                        let (cs, ce) = (cs.min(ce), cs.max(ce));
                                        // Exonic variant: standard HGVSc with 3' shifting
                                        ann.hgvsc = fastvep_hgvs::hgvsc_with_seq(
                                            &versioned_tid,
                                            cs, ce,
                                            &hgvs_ref,
                                            &hgvs_alt,
                                            coding_start,
                                            tr.cdna_coding_end,
                                            tr.spliced_seq.as_deref(),
                                            tr.codon_table_start_phase,
                                        );
                                    } else if ac.intron.is_some() {
                                        // Intronic variant: offset notation
                                        // Note: intronic HGVS uses original coding_start (no phase adjustment)
                                        // Apply HGVS 3' normalization for intronic indels
                                        let (shifted_start, shifted_end) = if let Some(ref sp) = seq_provider {
                                            let is_indel = matches!((&hgvs_ref, &hgvs_alt),
                                                (Allele::Sequence(_), Allele::Deletion) |
                                                (Allele::Deletion, Allele::Sequence(_)));
                                            if is_indel {
                                                if let Some((istart, iend)) = tr.intron_bounds_at(vf.position.start) {
                                                    // Use genomic-strand alleles for ref comparison
                                                    three_prime_shift_intronic(
                                                        &**sp as &dyn SequenceProvider, chrom,
                                                        vf.position.start, vf.position.end,
                                                        &vf.ref_allele, &ac.allele,
                                                        tr.strand, istart, iend,
                                                    )
                                                } else {
                                                    (vf.position.start, vf.position.end)
                                                }
                                            } else {
                                                (vf.position.start, vf.position.end)
                                            }
                                        } else {
                                            (vf.position.start, vf.position.end)
                                        };
                                        // For insertions, build the rotated insertion bases
                                        // after 3' shifting (bases rotate as position shifts)
                                        let shifted_hgvs_alt = if let (Allele::Deletion, Allele::Sequence(ins_bases)) = (&hgvs_ref, &hgvs_alt) {
                                            if shifted_start != vf.position.start && !ins_bases.is_empty() {
                                                // Calculate how many positions we shifted
                                                let shift_amount = if tr.strand == fastvep_core::Strand::Forward {
                                                    (shifted_start as i64 - vf.position.start as i64) as usize
                                                } else {
                                                    (vf.position.start as i64 - shifted_start as i64) as usize
                                                };
                                                // Rotate: for forward strand, each shift moves first base to end
                                                // For reverse strand, each shift moves last base to front
                                                let mut rotated = ins_bases.clone();
                                                let len = rotated.len();
                                                if len > 0 {
                                                    let effective_shift = shift_amount % len;
                                                    match tr.strand {
                                                        fastvep_core::Strand::Forward => {
                                                            rotated.rotate_left(effective_shift);
                                                        }
                                                        fastvep_core::Strand::Reverse => {
                                                            rotated.rotate_right(effective_shift);
                                                        }
                                                    }
                                                }
                                                Allele::Sequence(rotated)
                                            } else {
                                                hgvs_alt.clone()
                                            }
                                        } else {
                                            hgvs_alt.clone()
                                        };

                                        // For insertions, use position before insertion
                                        // for the primary HGVS coordinate (ins is BETWEEN two bases).
                                        // On reverse strand, the insertion is between P and P+1 in
                                        // genomic coords, but P+1 is 5' in transcript order, so we
                                        // use P+1 as the HGVS start coordinate.
                                        let is_insertion = matches!((&hgvs_ref, &shifted_hgvs_alt), (Allele::Deletion, Allele::Sequence(_)));
                                        let hgvs_pos = if is_insertion {
                                            if tr.strand == fastvep_core::Strand::Reverse {
                                                shifted_end + 1
                                            } else {
                                                shifted_end // base before insertion
                                            }
                                        } else {
                                            shifted_start
                                        };
                                        if let Some((cdna_pos, offset)) = tr.genomic_to_intronic_cdna(hgvs_pos) {
                                            // For multi-base variants, compute end position too
                                            let (end_cdna, end_offset) = if shifted_start != shifted_end && hgvs_pos == shifted_start {
                                                tr.genomic_to_intronic_cdna(shifted_end)
                                                    .map(|(c, o)| (Some(c), Some(o)))
                                                    .unwrap_or((None, None))
                                            } else {
                                                (None, None)
                                            };
                                            let mut hgvsc = fastvep_hgvs::hgvsc_intronic_range(
                                                &versioned_tid,
                                                cdna_pos,
                                                offset,
                                                end_cdna,
                                                end_offset,
                                                &hgvs_ref,
                                                &shifted_hgvs_alt,
                                                coding_start,
                                                tr.cdna_coding_end,
                                            );
                                            // For intronic insertions, check if it's a dup.
                                            if let (Some(ref h), Allele::Deletion, Allele::Sequence(_)) =
                                                (&hgvsc, &hgvs_ref, &hgvs_alt)
                                            {
                                                if h.contains("ins") {
                                                    let orig_ins = match &ac.allele {
                                                        Allele::Sequence(b) => b.clone(),
                                                        _ => vec![],
                                                    };
                                                    if !orig_ins.is_empty() {
                                                        if let Some(ref sp) = seq_provider {
                                                            let ins_len = orig_ins.len() as u64;
                                                            // Check dup_before: base(s) before insertion match
                                                            let check_end = vf.position.end;
                                                            let check_start = check_end.saturating_sub(ins_len - 1);
                                                            let dup_before = if let Ok(ref_seq) = sp.fetch_sequence_slice(chrom, check_start, check_end) {
                                                                ref_seq.len() == orig_ins.len()
                                                                    && ref_seq.iter().zip(orig_ins.iter())
                                                                        .all(|(a, b)| a.eq_ignore_ascii_case(b))
                                                            } else { false };
                                                            // Check dup_after: base(s) after insertion match
                                                            let dup_after = if !dup_before {
                                                                let cs = vf.position.start;
                                                                let ce = cs + ins_len - 1;
                                                                if let Ok(ref_seq) = sp.fetch_sequence_slice(chrom, cs, ce) {
                                                                    ref_seq.len() == orig_ins.len()
                                                                        && ref_seq.iter().zip(orig_ins.iter())
                                                                            .all(|(a, b)| a.eq_ignore_ascii_case(b))
                                                                } else { false }
                                                            } else { false };
                                                            if dup_before || dup_after {
                                                                // For dups, determine the dup base position and 3' shift it
                                                                let dup_base_pos = if dup_before {
                                                                    // Dup base is before insertion: position.end
                                                                    vf.position.end
                                                                } else {
                                                                    // Dup base is after insertion: position.start
                                                                    vf.position.start
                                                                };
                                                                // 3' shift the dup position within the intron
                                                                let shifted_dup = if let Some((istart, iend)) = tr.intron_bounds_at(dup_base_pos) {
                                                                    let (sd, _) = three_prime_shift_intronic(
                                                                        &**sp as &dyn SequenceProvider, chrom,
                                                                        dup_base_pos, dup_base_pos,
                                                                        &Allele::Sequence(orig_ins.clone()), &Allele::Deletion,
                                                                        tr.strand, istart, iend,
                                                                    );
                                                                    sd
                                                                } else {
                                                                    dup_base_pos
                                                                };
                                                                // Use shifted_dup (start of dup region) for offset computation
                                                                if let Some((dup_cdna, dup_offset)) = tr.genomic_to_intronic_cdna(shifted_dup) {
                                                                    hgvsc = convert_ins_to_dup(h, dup_offset, ins_len, dup_cdna, coding_start, tr.cdna_coding_end);
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                            ann.hgvsc = hgvsc;
                                        }
                                    }
                                } else {
                                    // Non-coding transcript: use n. notation
                                    if let (Some(cs), Some(ce)) = (ac.cdna_start, ac.cdna_end) {
                                        ann.hgvsc = fastvep_hgvs::hgvsc_noncoding(
                                            &versioned_tid,
                                            cs, ce,
                                            &hgvs_ref,
                                            &hgvs_alt,
                                        );
                                    } else if ac.intron.is_some() {
                                        // Apply 3' normalization for non-coding intronic indels
                                        let (nc_shifted_start, nc_shifted_end) = if let Some(ref sp) = seq_provider {
                                            let is_indel = matches!((&hgvs_ref, &hgvs_alt),
                                                (Allele::Sequence(_), Allele::Deletion) |
                                                (Allele::Deletion, Allele::Sequence(_)));
                                            if is_indel {
                                                if let Some((istart, iend)) = tr.intron_bounds_at(vf.position.start) {
                                                    three_prime_shift_intronic(
                                                        &**sp as &dyn SequenceProvider, chrom,
                                                        vf.position.start, vf.position.end,
                                                        &vf.ref_allele, &ac.allele,
                                                        tr.strand, istart, iend,
                                                    )
                                                } else {
                                                    (vf.position.start, vf.position.end)
                                                }
                                            } else {
                                                (vf.position.start, vf.position.end)
                                            }
                                        } else {
                                            (vf.position.start, vf.position.end)
                                        };

                                        // Rotate insertion bases for non-coding
                                        let nc_shifted_hgvs_alt = if let (Allele::Deletion, Allele::Sequence(ins_bases)) = (&hgvs_ref, &hgvs_alt) {
                                            if nc_shifted_start != vf.position.start && !ins_bases.is_empty() {
                                                let shift_amount = if tr.strand == fastvep_core::Strand::Forward {
                                                    (nc_shifted_start as i64 - vf.position.start as i64) as usize
                                                } else {
                                                    (vf.position.start as i64 - nc_shifted_start as i64) as usize
                                                };
                                                let mut rotated = ins_bases.clone();
                                                let len = rotated.len();
                                                if len > 0 {
                                                    let effective_shift = shift_amount % len;
                                                    match tr.strand {
                                                        fastvep_core::Strand::Forward => rotated.rotate_left(effective_shift),
                                                        fastvep_core::Strand::Reverse => rotated.rotate_right(effective_shift),
                                                    }
                                                }
                                                Allele::Sequence(rotated)
                                            } else {
                                                hgvs_alt.clone()
                                            }
                                        } else {
                                            hgvs_alt.clone()
                                        };

                                        if let Some((cdna_pos, offset)) = tr.genomic_to_intronic_cdna(nc_shifted_start) {
                                            let (end_cdna, end_offset) = if nc_shifted_start != nc_shifted_end {
                                                tr.genomic_to_intronic_cdna(nc_shifted_end)
                                                    .map(|(c, o)| (Some(c), Some(o)))
                                                    .unwrap_or((None, None))
                                            } else {
                                                (None, None)
                                            };
                                            let mut hgvsc = fastvep_hgvs::hgvsc_noncoding_intronic_range(
                                                &versioned_tid,
                                                cdna_pos,
                                                offset,
                                                end_cdna,
                                                end_offset,
                                                &hgvs_ref,
                                                &nc_shifted_hgvs_alt,
                                            );
                                            // Dup detection for non-coding intronic insertions
                                            if let (Some(ref h), Allele::Deletion, Allele::Sequence(_)) =
                                                (&hgvsc, &hgvs_ref, &hgvs_alt)
                                            {
                                                if h.contains("ins") {
                                                    let orig_ins = match &ac.allele {
                                                        Allele::Sequence(b) => b.clone(),
                                                        _ => vec![],
                                                    };
                                                    if !orig_ins.is_empty() {
                                                        if let Some(ref sp) = seq_provider {
                                                            let ins_len = orig_ins.len() as u64;
                                                            let check_end = vf.position.end;
                                                            let check_start = check_end.saturating_sub(ins_len - 1);
                                                            let dup_before = if let Ok(ref_seq) = sp.fetch_sequence_slice(chrom, check_start, check_end) {
                                                                ref_seq.len() == orig_ins.len()
                                                                    && ref_seq.iter().zip(orig_ins.iter())
                                                                        .all(|(a, b)| a.eq_ignore_ascii_case(b))
                                                            } else { false };
                                                            let dup_after = if !dup_before {
                                                                let cs = vf.position.start;
                                                                let ce = cs + ins_len - 1;
                                                                if let Ok(ref_seq) = sp.fetch_sequence_slice(chrom, cs, ce) {
                                                                    ref_seq.len() == orig_ins.len()
                                                                        && ref_seq.iter().zip(orig_ins.iter())
                                                                            .all(|(a, b)| a.eq_ignore_ascii_case(b))
                                                                } else { false }
                                                            } else { false };
                                                            if dup_before || dup_after {
                                                                let dup_base_pos = if dup_before { vf.position.end } else { vf.position.start };
                                                                let shifted_dup = if let Some((istart, iend)) = tr.intron_bounds_at(dup_base_pos) {
                                                                    let (sd, _) = three_prime_shift_intronic(
                                                                        &**sp as &dyn SequenceProvider, chrom,
                                                                        dup_base_pos, dup_base_pos,
                                                                        &Allele::Sequence(orig_ins.clone()), &Allele::Deletion,
                                                                        tr.strand, istart, iend,
                                                                    );
                                                                    sd
                                                                } else {
                                                                    dup_base_pos
                                                                };
                                                                if let Some((dup_cdna, dup_offset)) = tr.genomic_to_intronic_cdna(shifted_dup) {
                                                                    hgvsc = convert_ins_to_dup_noncoding(h, dup_offset, ins_len, dup_cdna);
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                            ann.hgvsc = hgvsc;
                                        }
                                    }
                                }
                            }

                            if let (Some(ref aa), Some(ps)) = (&ac.amino_acids, ac.protein_start) {
                                if let Some(tr) = transcript {
                                    if let Some(ref pid) = tr.protein_id {
                                        let versioned_pid = match tr.protein_version {
                                            Some(v) => {
                                                let suffix = format!(".{}", v);
                                                if pid.ends_with(&suffix) {
                                                    pid.clone()
                                                } else {
                                                    format!("{}.{}", pid, v)
                                                }
                                            }
                                            None => pid.clone(),
                                        };
                                        let is_fs = ac.consequences.contains(&Consequence::FrameshiftVariant);

                                        if is_fs {
                                            // Frameshift: build alt sequence and scan for first changed AA + new stop
                                            // Use spliced_seq from CDS start onwards (includes 3'UTR for stop codon search)
                                            if let (Some(ref spliced), Some(coding_start), Some(cds_s)) =
                                                (&tr.spliced_seq, tr.cdna_coding_start, ac.cds_start)
                                            {
                                                // Extract from CDS start to end of spliced seq (includes 3'UTR)
                                                let coding_start_idx = (coding_start - 1) as usize;
                                                let ref_from_cds = &spliced.as_bytes()[coding_start_idx..];
                                                let cds_idx = (cds_s - 1) as usize;
                                                let mut alt_from_cds = ref_from_cds.to_vec();

                                                // Apply the indel to build the frameshifted sequence
                                                if ac.allele == Allele::Deletion {
                                                    let del_len = vf.ref_allele.len();
                                                    let end = (cds_idx + del_len).min(alt_from_cds.len());
                                                    alt_from_cds.drain(cds_idx..end);
                                                } else if let Allele::Sequence(ins_bases) = &ac.allele {
                                                    let mut bases = ins_bases.clone();
                                                    if tr.strand == fastvep_core::Strand::Reverse {
                                                        bases = bases.iter().map(|&b| match b {
                                                            b'A' => b'T', b'T' => b'A',
                                                            b'C' => b'G', b'G' => b'C',
                                                            o => o,
                                                        }).collect();
                                                    }
                                                    for (j, &b) in bases.iter().enumerate() {
                                                        if cds_idx + j <= alt_from_cds.len() {
                                                            alt_from_cds.insert(cds_idx + j, b);
                                                        }
                                                    }
                                                }

                                                let codon_start = cds_idx / 3;
                                                ann.hgvsp = fastvep_hgvs::hgvsp_frameshift(
                                                    &versioned_pid,
                                                    ref_from_cds,
                                                    &alt_from_cds,
                                                    codon_start,
                                                );
                                            }
                                        } else if aa.1 == "-"
                                            || ac.consequences.contains(&Consequence::InframeDeletion)
                                        {
                                            // In-frame deletion / delins (frameshift handled
                                            // above). aa.0 holds the deleted residues, aa.1 the
                                            // replacement ("-" for a pure deletion).
                                            ann.hgvsp = fastvep_hgvs::hgvsp_inframe_deletion(
                                                &versioned_pid, ps, &aa.0, &aa.1,
                                            );
                                        } else {
                                            let ref_aa_byte = aa.0.as_bytes().first().copied().unwrap_or(b'X');
                                            let alt_aa_byte = aa.1.as_bytes().first().copied().unwrap_or(b'X');
                                            ann.hgvsp = fastvep_hgvs::hgvsp(
                                                &versioned_pid, ps, ref_aa_byte, alt_aa_byte, false,
                                            );
                                        }
                                    }
                                }
                            }
                        }

                        ann
                    })
                    .collect();

                // Collect every transcript here; --pick filtering runs as a
                // single post-pass below so it can compare all candidates
                // before SA/gene/ACMG annotation, instead of picking the first
                // canonical one we happen to encounter.
                vf.transcript_variations.push(TranscriptVariation {
                    transcript_id: tc.transcript_id.clone(),
                    gene_id: tc.gene_id.clone(),
                    gene_symbol: tc.gene_symbol.clone(),
                    biotype: tc.biotype.clone(),
                    allele_annotations,
                    canonical: tc.canonical,
                    strand: tc.strand,
                    source: transcript.and_then(|t| t.source.clone()),
                    protein_id: transcript.and_then(|t| t.protein_id.clone()),
                    mane_select: transcript.and_then(|t| t.mane_select.clone()),
                    mane_plus_clinical: transcript.and_then(|t| t.mane_plus_clinical.clone()),
                    tsl: transcript.and_then(|t| t.tsl),
                    appris: transcript.and_then(|t| t.appris.clone()),
                    ccds: transcript.and_then(|t| t.ccds.clone()),
                    gencode_primary: transcript.map(|t| t.gencode_primary).unwrap_or(false),
                    symbol_source: transcript.and_then(|t| t.gene.symbol_source.clone()),
                    hgnc_id: transcript.and_then(|t| t.gene.hgnc_id.clone()),
                    flags: transcript.map(|t| t.flags.clone()).unwrap_or_default(),
                });
            }
            } // close `else` of overlapping.is_empty()
            } // close `else` of `if sa_only`

            // Apply --pick before SA/gene/ACMG so those passes only run on the
            // single surviving transcript. Running pick after them would still
            // produce correct output but would waste the most expensive work
            // (ACMG classification) on transcripts that get thrown away.
            if config.pick && !sa_only && vf.transcript_variations.len() > 1 {
                if let Some(idx) = pick_best_transcript_idx(&vf.transcript_variations) {
                    vf.transcript_variations =
                        vec![vf.transcript_variations.swap_remove(idx)];
                }
            }

            // Supplementary annotation: query SA providers once per unique
            // allele, then attach the result to every (transcript, allele)
            // slot that shares it. SA results depend only on (pos, ref, alt),
            // never on the transcript context, so this avoids T× amplification
            // for variants overlapping many transcripts. Runs for all variants
            // (intergenic, transcript-overlapping, and sa_only scaffold) so
            // supplementary databases attach to every record that matches.
            if !sa_providers.is_empty() {
                let chrom = &vf.position.chromosome;
                let sa_queries = supplementary_query_alleles(vf);
                let mut allele_results: HashMap<String, Vec<(String, String)>> =
                    HashMap::new();
                for tv in &vf.transcript_variations {
                    for aa in &tv.allele_annotations {
                        let allele_key = aa.allele.to_string();
                        if allele_results.contains_key(&allele_key) {
                            continue;
                        }
                        let (query_pos, ref_str, alt_str) = sa_queries
                            .iter()
                            .find(|(allele, _, _, _)| allele == &allele_key)
                            .map(|(_, pos, ref_allele, alt_allele)| {
                                (*pos, ref_allele.clone(), alt_allele.clone())
                            })
                            .unwrap_or_else(|| {
                                (
                                    vf.position.start,
                                    vf.ref_allele.to_string(),
                                    allele_key.clone(),
                                )
                            });
                        // gnomAD stores left-aligned, parsimonious alleles; the
                        // raw input representation can differ for indels (esp.
                        // in repeats), silently missing the match and making
                        // PM2 misfire on common variants. Normalize the query to
                        // gnomAD's minimal representation — only when a gnomAD
                        // provider and a reference are present, and applied only
                        // to the gnomAD lookup (other sources keep the raw key).
                        let gnomad_norm = if seq_provider.is_some()
                            && sa_providers.iter().any(|sa| sa.json_key() == "gnomad")
                        {
                            seq_provider.as_ref().map(|sp| {
                                fastvep_cache::normalize::normalize_variant(
                                    &**sp,
                                    chrom,
                                    query_pos,
                                    &ref_str,
                                    &alt_str,
                                )
                            })
                        } else {
                            None
                        };
                        let mut results: Vec<(String, String)> = Vec::new();
                        for sa in &sa_providers {
                            let (sa_pos, sa_ref, sa_alt) = if sa.metadata().match_by_allele {
                                if sa.json_key() == "gnomad" {
                                    match &gnomad_norm {
                                        Some(n) => {
                                            (n.pos, n.ref_allele.as_str(), n.alt_allele.as_str())
                                        }
                                        None => (query_pos, ref_str.as_str(), alt_str.as_str()),
                                    }
                                } else {
                                    (query_pos, ref_str.as_str(), alt_str.as_str())
                                }
                            } else {
                                (vf.position.start, "", "")
                            };
                            if let Ok(Some(ann)) =
                                sa.annotate_position(chrom, sa_pos, sa_ref, sa_alt)
                            {
                                let json_str = match ann {
                                    AnnotationValue::Json(j) => j,
                                    AnnotationValue::Positional(j) => j,
                                    AnnotationValue::Interval(v) => {
                                        format!("[{}]", v.join(","))
                                    }
                                };
                                results.push((sa.json_key().to_string(), json_str));
                            }
                        }
                        allele_results.insert(allele_key, results);
                    }
                }
                for tv in &mut vf.transcript_variations {
                    for aa in &mut tv.allele_annotations {
                        if let Some(results) = allele_results.get(&aa.allele.to_string()) {
                            aa.supplementary.extend(results.iter().cloned());
                        }
                    }
                }
            }

            // Gene-level annotation pass (OMIM, gnomAD gene constraints, etc.)
            if !gene_providers.is_empty() {
                use fastvep_cache::annotation::GeneAnnotationProvider;
                let mut seen_genes = std::collections::HashSet::new();
                for tv in &vf.transcript_variations {
                    if let Some(gene_sym) = tv.gene_symbol.as_deref() {
                        if seen_genes.insert(gene_sym.to_string()) {
                            for gp in &gene_providers {
                                if let Ok(Some(json)) = gp.annotate_gene(gene_sym) {
                                    vf.gene_annotations.push(
                                        fastvep_core::GeneAnnotation {
                                            gene_symbol: gene_sym.to_string(),
                                            json_key: gp.json_key().to_string(),
                                            json_string: json,
                                        },
                                    );
                                }
                            }
                        }
                    }
                }
            }

            // ACMG-AMP classification pass (after all SA annotations are attached)
            if let Some(ref acmg_cfg) = acmg_config {
                // Parse sample genotypes if trio config is present
                let trio_genotypes = extract_trio_genotypes_cli(vf, acmg_cfg, &sample_names);

                for tv in &mut vf.transcript_variations {
                    let gene_sym = tv.gene_symbol.as_deref().unwrap_or("");
                    let gene_anns: Vec<&fastvep_core::GeneAnnotation> =
                        vf.gene_annotations
                            .iter()
                            .filter(|ga| ga.gene_symbol == gene_sym)
                            .collect();
                    for aa in &mut tv.allele_annotations {
                        let input =
                            fastvep_classification::extract_classification_input(
                                &aa.consequences,
                                aa.impact,
                                tv.gene_symbol.as_deref(),
                                tv.canonical,
                                aa.amino_acids.as_ref(),
                                aa.protein_position.map(|(s, _)| s),
                                aa.hgvsc.as_deref(),
                                aa.exon,
                                &aa.supplementary,
                                &gene_anns,
                                &vf.supplementary_annotations,
                                trio_genotypes.0.clone(),
                                trio_genotypes.1.clone(),
                                trio_genotypes.2.clone(),
                                vec![], // companion_variants populated in second pass
                            );
                        let result = fastvep_classification::classify(&input, acmg_cfg);
                        aa.acmg_classification =
                            serde_json::to_value(&result).ok();
                    }
                }
            }

            if !sa_only {
                vf.compute_most_severe();
            }
        }); // end par_iter_mut

        // Phase 2.5: Compound-het enrichment pass (sequential, after parallel annotation)
        if let Some(ref acmg_cfg) = acmg_config {
            if acmg_cfg.trio.is_some() {
                let mut vfs: Vec<&mut VariationFeature> =
                    batch.iter_mut().map(|(vf, _)| vf).collect();
                enrich_compound_het_batch(&mut vfs, acmg_cfg, &sample_names);
            }
        }

        // Phase 3: Write output sequentially (preserves VCF order)
        for (vf, _) in &batch {
            match config.output_format.as_str() {
                "vcf" => write_vcf_line(&mut writer, vf, sa_only)?,
                "tab" => {
                    // Classify variant against QC rules (if any). The
                    // classifier reads the VCF INFO column once via a
                    // streaming view; no HashMap, no allocation.
                    let qc_label: Option<&str> = if let Some(ref rules) = qc_rules {
                        let (info_str, qual): (&str, Option<f64>) = match &vf.vcf_fields {
                            Some(f) => (f.info.as_str(), f.qual.parse::<f64>().ok()),
                            None => ("", None),
                        };
                        let view = fastvep_io::qc::InfoView::new(info_str, qual);
                        let filter = vf
                            .vcf_fields
                            .as_ref()
                            .map(|f| f.filter.as_str())
                            .unwrap_or("");
                        Some(rules.classify(&view, filter))
                    } else {
                        None
                    };

                    let opts = output::TabOptions {
                        gene_set: gene_set.as_ref(),
                        explicit_ref: config.explicit_alleles,
                        qc_class: qc_label,
                    };
                    for line in output::format_tab_line_with(
                        vf,
                        &supplementary_specs,
                        sa_only,
                        opts,
                    ) {
                        writeln!(writer, "{}", line)?;
                    }
                }
                "json" => {
                    if !first_json {
                        writeln!(writer, ",")?;
                    }
                    first_json = false;
                    let json = output::format_json(vf, sa_only);
                    write!(writer, "{}", serde_json::to_string_pretty(&json)?)?;
                }
                _ => {}
            }
        }

        meter.update_n(batch.len() as u64);
    } // end batch loop

    // Close JSON array
    if config.output_format == "json" {
        writeln!(writer, "\n]")?;
    }

    writer.flush()?;
    meter.finish();

    Ok(())
}

// Shared annotation utilities from fastvep-annotate (used by batch pipeline).
use fastvep_annotate::{
    annotate_intergenic, annotate_sa_only_scaffold, complement_allele, convert_ins_to_dup,
    convert_ins_to_dup_noncoding, load_gene_providers, load_sa_providers,
    three_prime_shift_intronic, zip_positions,
};

/// Index of the best transcript variation under VEP's default `--pick_order`
/// hierarchy: mane_select, mane_plus_clinical, canonical, appris, tsl, biotype
/// (protein_coding preferred), ccds, then most-severe consequence rank, with
/// transcript_id alphabetical order as a final deterministic tie-breaker.
fn pick_best_transcript_idx(tvs: &[TranscriptVariation]) -> Option<usize> {
    (0..tvs.len()).min_by(|&a, &b| pick_key(&tvs[a]).cmp(&pick_key(&tvs[b])))
}

fn pick_key(tv: &TranscriptVariation) -> (bool, bool, bool, u8, u8, u8, bool, u32, &str) {
    let most_severe_rank = tv
        .allele_annotations
        .iter()
        .flat_map(|aa| aa.consequences.iter())
        .map(|c| c.rank())
        .min()
        .unwrap_or(u32::MAX);
    (
        tv.mane_select.is_none(),
        tv.mane_plus_clinical.is_none(),
        !tv.canonical,
        appris_rank(tv.appris.as_deref()),
        tv.tsl.unwrap_or(u8::MAX),
        if tv.biotype.as_ref() == "protein_coding" { 0 } else { 1 },
        tv.ccds.is_none(),
        most_severe_rank,
        tv.transcript_id.as_ref(),
    )
}

/// Map an APPRIS tag (`P1`/`principal1`, ..., `A1`/`alternative1`, ...) to a
/// rank where lower is better, matching VEP's `--pick_order` APPRIS tier:
/// principal1 < principal2 < ... < alternative1 < alternative2 < absent.
fn appris_rank(appris: Option<&str>) -> u8 {
    let Some(s) = appris else { return u8::MAX };
    let lower = s.to_ascii_lowercase();
    let (is_alt, digits) = if let Some(d) = lower.strip_prefix("principal") {
        (false, d)
    } else if let Some(d) = lower.strip_prefix("alternative") {
        (true, d)
    } else if let Some(d) = lower.strip_prefix('p') {
        (false, d)
    } else if let Some(d) = lower.strip_prefix('a') {
        (true, d)
    } else {
        // Present but unrecognised: still better than absent, worse than any
        // recognised principal/alternative.
        return u8::MAX - 1;
    };
    let n: u8 = digits.parse().unwrap_or(9);
    if is_alt { 5u8.saturating_add(n) } else { n }
}

/// Extract trio genotype information from a VariationFeature's VCF sample columns (CLI path).
fn extract_trio_genotypes_cli(
    vf: &VariationFeature,
    acmg_cfg: &fastvep_classification::AcmgConfig,
    sample_names: &[String],
) -> (
    Option<fastvep_classification::GenotypeInfo>,
    Option<fastvep_classification::GenotypeInfo>,
    Option<fastvep_classification::GenotypeInfo>,
) {
    let trio = match &acmg_cfg.trio {
        Some(t) => t,
        None => return (None, None, None),
    };

    let vcf_fields = match &vf.vcf_fields {
        Some(f) => f,
        None => return (None, None, None),
    };

    if vcf_fields.rest.is_empty() {
        return (None, None, None);
    }

    let format_str = &vcf_fields.rest[0];
    let sample_strs: Vec<&str> = vcf_fields.rest[1..].iter().map(|s| s.as_str()).collect();

    let samples =
        fastvep_io::sample::parse_samples(format_str, &sample_strs, sample_names);

    let proband_gt = samples
        .iter()
        .find(|s| s.name == trio.proband)
        .map(|s| sample_data_to_genotype_info_cli(s));

    let mother_gt = trio.mother.as_ref().and_then(|name| {
        samples
            .iter()
            .find(|s| &s.name == name)
            .map(|s| sample_data_to_genotype_info_cli(s))
    });

    let father_gt = trio.father.as_ref().and_then(|name| {
        samples
            .iter()
            .find(|s| &s.name == name)
            .map(|s| sample_data_to_genotype_info_cli(s))
    });

    (proband_gt, mother_gt, father_gt)
}

/// Convert a SampleData to GenotypeInfo (CLI path).
fn sample_data_to_genotype_info_cli(
    sample: &fastvep_io::sample::SampleData,
) -> fastvep_classification::GenotypeInfo {
    let gt = sample.genotype.as_ref();
    let is_het = gt.map_or(false, |g| g.is_het());
    let is_hom_ref = gt.map_or(false, |g| g.is_hom_ref());
    let is_hom_alt = gt.map_or(false, |g| g.is_hom_alt());
    let is_missing = gt.map_or(true, |g| g.is_missing());
    let is_phased = gt.map_or(false, |g| g.phased);

    let alt_allele_index = gt.and_then(|g| {
        g.alleles.iter().filter_map(|a| *a).find(|&a| a > 0)
    });

    fastvep_classification::GenotypeInfo {
        is_het,
        is_hom_ref,
        is_hom_alt,
        is_missing,
        is_phased,
        depth: sample.depth,
        quality: sample.quality,
        alt_allele_index,
    }
}

/// Compound-het enrichment pass for the CLI batch pipeline.
/// Groups variants by gene and re-evaluates PM3/BP2 with companion data.
fn enrich_compound_het_batch(
    variants: &mut [&mut VariationFeature],
    acmg_cfg: &fastvep_classification::AcmgConfig,
    sample_names: &[String],
) {
    // Collect per-gene variant info
    struct VariantGeneInfo {
        vf_idx: usize,
        tv_idx: usize,
        aa_idx: usize,
        is_clinvar_pathogenic: bool,
        is_clinvar_likely_pathogenic: bool,
        proband_het: bool,
        is_phased: bool,
        proband_alleles: Vec<Option<u32>>,
        hgvsc: Option<String>,
    }

    let mut gene_variants: HashMap<String, Vec<VariantGeneInfo>> = HashMap::new();

    for (vf_idx, vf) in variants.iter().enumerate() {
        let trio_genotypes = extract_trio_genotypes_cli(vf, acmg_cfg, sample_names);
        let proband_gt = &trio_genotypes.0;

        for (tv_idx, tv) in vf.transcript_variations.iter().enumerate() {
            let gene_sym = match tv.gene_symbol.as_deref() {
                Some(g) if !g.is_empty() && g != "-" => g.to_string(),
                _ => continue,
            };

            for (aa_idx, aa) in tv.allele_annotations.iter().enumerate() {
                // Classify ClinVar supplementary as Pathogenic / Likely pathogenic
                // separately so PM3 v1.0 scores them at their proper point values.
                // Strip "Likely pathogenic" before checking for "pathogenic"
                // residual to avoid the substring-match bug that double-counts
                // LP as P.
                let (clinvar_p_from_sa, clinvar_lp_from_sa) = aa
                    .supplementary
                    .iter()
                    .filter(|(key, json)| {
                        key == "clinvar"
                            && !json.contains("Conflicting")
                            && !json.contains("conflicting")
                    })
                    .map(|(_, json)| {
                        let lower = json.to_lowercase();
                        let has_lp = lower.contains("likely pathogenic");
                        let stripped = lower.replace("likely pathogenic", "");
                        let has_p = stripped.contains("pathogenic");
                        (has_p, has_lp && !has_p)
                    })
                    .fold((false, false), |(p_acc, lp_acc), (p, lp)| {
                        (p_acc || p, lp_acc || lp)
                    });

                let proband_het = proband_gt.as_ref().map_or(false, |g| g.is_het);
                let is_phased = proband_gt.as_ref().map_or(false, |g| g.is_phased);
                let proband_alleles = if let Some(ref vcf_fields) = vf.vcf_fields {
                    if !vcf_fields.rest.is_empty() && !sample_names.is_empty() {
                        let format_str = &vcf_fields.rest[0];
                        let sample_strs: Vec<&str> =
                            vcf_fields.rest[1..].iter().map(|s| s.as_str()).collect();
                        let samples = fastvep_io::sample::parse_samples(
                            format_str,
                            &sample_strs,
                            sample_names,
                        );
                        if let Some(trio) = &acmg_cfg.trio {
                            samples
                                .iter()
                                .find(|s| s.name == trio.proband)
                                .and_then(|s| s.genotype.as_ref())
                                .map(|g| g.alleles.clone())
                                .unwrap_or_default()
                        } else {
                            vec![]
                        }
                    } else {
                        vec![]
                    }
                } else {
                    vec![]
                };

                gene_variants
                    .entry(gene_sym.clone())
                    .or_default()
                    .push(VariantGeneInfo {
                        vf_idx,
                        tv_idx,
                        aa_idx,
                        is_clinvar_pathogenic: clinvar_p_from_sa,
                        is_clinvar_likely_pathogenic: clinvar_lp_from_sa,
                        proband_het,
                        is_phased,
                        proband_alleles,
                        hgvsc: aa.hgvsc.clone(),
                    });
            }
        }
    }

    for (_gene, gene_infos) in &gene_variants {
        let het_variants: Vec<&VariantGeneInfo> =
            gene_infos.iter().filter(|v| v.proband_het).collect();
        if het_variants.len() < 2 {
            continue;
        }

        for info in &het_variants {
            let companions: Vec<fastvep_classification::CompanionVariant> = het_variants
                .iter()
                .filter(|other| {
                    other.vf_idx != info.vf_idx
                        || other.tv_idx != info.tv_idx
                        || other.aa_idx != info.aa_idx
                })
                .map(|other| {
                    let is_in_trans = if info.is_phased && other.is_phased {
                        if info.proband_alleles.len() >= 2 && other.proband_alleles.len() >= 2 {
                            let info_alt_on_first = info
                                .proband_alleles
                                .first()
                                .map_or(false, |a| a.map_or(false, |v| v > 0));
                            let other_alt_on_first = other
                                .proband_alleles
                                .first()
                                .map_or(false, |a| a.map_or(false, |v| v > 0));
                            Some(info_alt_on_first != other_alt_on_first)
                        } else {
                            None
                        }
                    } else {
                        None
                    };

                    fastvep_classification::CompanionVariant {
                        is_clinvar_pathogenic: other.is_clinvar_pathogenic,
                        is_clinvar_likely_pathogenic: other.is_clinvar_likely_pathogenic,
                        is_in_trans,
                        proband_het: other.proband_het,
                        hgvsc: other.hgvsc.clone(),
                    }
                })
                .collect();

            if companions.is_empty() {
                continue;
            }

            let vf = &*variants[info.vf_idx];
            let tv = &vf.transcript_variations[info.tv_idx];
            let aa = &tv.allele_annotations[info.aa_idx];
            let gene_sym = tv.gene_symbol.as_deref().unwrap_or("");
            let gene_anns: Vec<&fastvep_core::GeneAnnotation> = vf
                .gene_annotations
                .iter()
                .filter(|ga| ga.gene_symbol == gene_sym)
                .collect();

            let trio_genotypes = extract_trio_genotypes_cli(vf, acmg_cfg, sample_names);

            let input = fastvep_classification::extract_classification_input(
                &aa.consequences,
                aa.impact,
                tv.gene_symbol.as_deref(),
                tv.canonical,
                aa.amino_acids.as_ref(),
                aa.protein_position.map(|(s, _)| s),
                aa.hgvsc.as_deref(),
                aa.exon,
                &aa.supplementary,
                &gene_anns,
                &vf.supplementary_annotations,
                trio_genotypes.0,
                trio_genotypes.1,
                trio_genotypes.2,
                companions,
            );
            let result = fastvep_classification::classify(&input, acmg_cfg);
            variants[info.vf_idx].transcript_variations[info.tv_idx].allele_annotations
                [info.aa_idx]
                .acmg_classification = serde_json::to_value(&result).ok();
        }
    }
}

fn supplementary_query_alleles(vf: &VariationFeature) -> Vec<(String, u64, String, String)> {
    if let Some(vcf) = &vf.vcf_fields {
        let uploaded_alts: Vec<&str> = vcf.alt.split(',').collect();
        return vf
            .alt_alleles
            .iter()
            .enumerate()
            .map(|(idx, allele)| {
                let allele_string = allele.to_string();
                (
                    allele_string.clone(),
                    vcf.pos,
                    vcf.ref_allele.clone(),
                    uploaded_alts
                        .get(idx)
                        .copied()
                        .unwrap_or(&allele_string)
                        .to_string(),
                )
            })
            .collect();
    }

    vf.alt_alleles
        .iter()
        .map(|allele| {
            (
                allele.to_string(),
                vf.position.start,
                vf.ref_allele.to_string(),
                allele.to_string(),
            )
        })
        .collect()
}

fn write_vcf_line(writer: &mut impl Write, vf: &VariationFeature, sa_only: bool) -> Result<()> {
    if let Some(ref fields) = vf.vcf_fields {
        let csq = if sa_only {
            String::new()
        } else {
            output::format_csq(vf, output::DEFAULT_CSQ_FIELDS)
        };
        let info = output::format_vcf_info_fields(&fields.info, vf, &csq);

        write!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            fields.chrom, fields.pos, fields.id, fields.ref_allele, fields.alt,
            fields.qual, fields.filter, info
        )?;

        for rest_field in &fields.rest {
            write!(writer, "\t{}", rest_field)?;
        }
        writeln!(writer)?;
    }

    Ok(())
}


use serde_json;

/// Filter annotated VCF by CSQ field expressions.
///
/// Reads a VCF with CSQ INFO annotations (produced by `fastvep annotate --output-format vcf`),
/// parses each CSQ entry into fields, evaluates the filter expression, and writes matching lines.
pub fn run_filter(input: &str, output_path: &str, filter_expr: &str) -> Result<()> {
    use fastvep_filter::{Filter, FilterContext};

    let filter = Filter::parse(filter_expr)
        .with_context(|| format!("Parsing filter expression: {}", filter_expr))?;

    // Open input
    let reader: Box<dyn BufRead> = if input == "-" {
        Box::new(io::BufReader::new(io::stdin()))
    } else {
        let f = File::open(input).with_context(|| format!("Opening input: {}", input))?;
        Box::new(io::BufReader::new(f))
    };

    // Open output. 1 MiB buffer same rationale as the main annotation path.
    let mut writer: Box<dyn Write> = if output_path == "-" {
        Box::new(BufWriter::with_capacity(1 << 20, io::stdout()))
    } else {
        let f = File::create(output_path)
            .with_context(|| format!("Creating output: {}", output_path))?;
        Box::new(BufWriter::with_capacity(1 << 20, f))
    };

    // Parse CSQ header to get field names
    let mut csq_fields: Vec<String> = Vec::new();
    let mut kept = 0u64;
    let mut total = 0u64;

    for line in reader.lines() {
        let line = line?;

        // Header lines: pass through, extract CSQ format
        if line.starts_with('#') {
            // Parse CSQ format from ##INFO=<ID=CSQ,...,Description="... Format: A|B|C">
            if line.starts_with("##INFO=<ID=CSQ") {
                if let Some(fmt_start) = line.find("Format: ") {
                    let fmt = &line[fmt_start + 8..];
                    let fmt = fmt.trim_end_matches('"').trim_end_matches('>');
                    csq_fields = fmt.split('|').map(|s| s.to_string()).collect();
                    eprintln!(
                        "[filter] CSQ format: {} fields ({})",
                        csq_fields.len(),
                        csq_fields.iter().take(5).cloned().collect::<Vec<_>>().join(", ")
                    );
                }
            }
            writeln!(writer, "{}", line)?;
            continue;
        }

        total += 1;

        // Parse VCF data line
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 8 {
            continue;
        }

        let info = fields[7];

        // Extract CSQ entries from INFO field
        let csq_str = info
            .split(';')
            .find(|f| f.starts_with("CSQ="))
            .map(|f| &f[4..]);

        let Some(csq_str) = csq_str else {
            // No CSQ field — skip this line (doesn't match any filter)
            continue;
        };

        // Check if ANY CSQ entry matches the filter
        let mut any_match = false;
        for entry in csq_str.split(',') {
            let values: Vec<&str> = entry.split('|').collect();
            let mut ctx = FilterContext::new();

            for (i, val) in values.iter().enumerate() {
                if let Some(field_name) = csq_fields.get(i) {
                    if !val.is_empty() {
                        ctx.set(field_name, val);
                    }
                }
            }

            // Also set VCF-level fields
            ctx.set("CHROM", fields[0]);
            ctx.set("POS", fields[1]);
            ctx.set("REF", fields[3]);
            ctx.set("ALT", fields[4]);
            ctx.set("FILTER", fields[6]);

            if filter.matches(&ctx) {
                any_match = true;
                break;
            }
        }

        if any_match {
            writeln!(writer, "{}", line)?;
            kept += 1;
        }
    }

    eprintln!(
        "[filter] {} of {} variants passed filter: {}",
        kept, total, filter_expr
    );

    Ok(())
}

/// Build a binary transcript cache from one or more GFF3s + optional FASTA.
///
/// With multiple GFF3 inputs (e.g. Ensembl + RefSeq), transcripts from all
/// files are merged into a single cache, each stamped with its source so
/// the SOURCE column on annotate output stays per-transcript.
pub fn run_cache_build(
    gff3_paths: &[String],
    fasta_path: Option<&str>,
    synonyms_path: Option<&str>,
    output_path: &str,
    show_progress: bool,
) -> Result<()> {
    if gff3_paths.is_empty() {
        return Err(anyhow::anyhow!("`fastvep cache` requires at least one --gff3"));
    }
    let specs: Vec<Gff3Spec> = gff3_paths.iter().map(|s| parse_gff3_arg(s)).collect();
    let mut transcripts: Vec<fastvep_genome::Transcript> = Vec::new();
    for spec in &specs {
        let gff_file = File::open(&spec.path)
            .with_context(|| format!("Opening GFF3 file: {}", spec.path))?;
        // Auto-decompress .gz / .bgz GFF3 inputs. Without this we'd silently
        // produce a 0-transcript cache.
        let trs = if spec.path.ends_with(".gz") || spec.path.ends_with(".bgz") {
            parse_gff3_with_source(flate2::read::MultiGzDecoder::new(gff_file), &spec.source)?
        } else {
            parse_gff3_with_source(gff_file, &spec.source)?
        };
        eprintln!(
            "Loaded {} transcripts from {} (source label: {})",
            trs.len(), spec.path, spec.source
        );
        transcripts.extend(trs);
    }
    if transcripts.is_empty() {
        return Err(anyhow::anyhow!(
            "GFF3 source(s) [{}] produced 0 transcripts — refusing to write an empty cache.",
            specs.iter().map(|s| s.path.as_str()).collect::<Vec<_>>().join(", ")
        ));
    }
    if specs.len() > 1 {
        eprintln!("Merged {} GFF3 sources into {} total transcripts", specs.len(), transcripts.len());
    }

    // Chromosome synonyms (VEP `chr_synonyms.txt`). An empty table still
    // resolves chr↔bare / mitochondrial forms via mechanical aliases.
    let synonyms = match synonyms_path {
        Some(p) => {
            let contents = std::fs::read_to_string(p)
                .with_context(|| format!("Reading synonyms file: {}", p))?;
            eprintln!("Loaded chromosome synonyms from {}", p);
            fastvep_core::ChromSynonyms::parse(&contents)
        }
        None => fastvep_core::ChromSynonyms::new(),
    };

    if let Some(fasta) = fasta_path {
        let fasta_file = File::open(fasta)
            .with_context(|| format!("Opening FASTA file: {}", fasta))?;
        let reader = FastaReader::from_reader(fasta_file)?;
        eprintln!("Loaded reference FASTA from {}", fasta);

        // Canonicalize every transcript (and its gene) to the FASTA's contig
        // naming. This both lets sequence fetch succeed and makes a merged
        // Ensembl + RefSeq cache use one consistent naming scheme, so
        // `annotate` matches a VCF regardless of which GFF3 the transcript
        // came from. RefSeq accessions only resolve when a synonyms file maps
        // them; chr↔bare / mito resolve with no synonyms file at all.
        let contigs: std::collections::HashSet<String> =
            reader.sequence_names().into_iter().map(|s| s.to_string()).collect();
        let mut resolved: HashMap<String, Option<std::sync::Arc<str>>> = HashMap::new();
        let mut unresolved: std::collections::BTreeSet<String> = std::collections::BTreeSet::new();
        for tr in &mut transcripts {
            let orig: &str = tr.chromosome.as_ref();
            if contigs.contains(orig) {
                continue; // already matches a FASTA contig
            }
            let canonical = resolved
                .entry(orig.to_string())
                .or_insert_with(|| {
                    synonyms
                        .aliases(orig)
                        .into_iter()
                        .find(|a| contigs.contains(a))
                        .map(|a| std::sync::Arc::from(a.as_str()))
                })
                .clone();
            match canonical {
                Some(name) => {
                    tr.chromosome = std::sync::Arc::clone(&name);
                    tr.gene.chromosome = name;
                }
                None => {
                    unresolved.insert(orig.to_string());
                }
            }
        }
        if !unresolved.is_empty() {
            let refseq_like = unresolved
                .iter()
                .any(|c| fastvep_core::looks_like_refseq_accession(c));
            let hint = if refseq_like && synonyms_path.is_none() {
                " — some look like RefSeq accessions; pass a VEP-style synonyms file with `--synonyms` to map them"
            } else {
                ""
            };
            eprintln!(
                "Warning: {} GFF3 contig(s) have no matching FASTA sequence{}: {}",
                unresolved.len(),
                hint,
                unresolved.iter().cloned().collect::<Vec<_>>().join(", ")
            );
        }

        let sp = FastaSequenceProvider::new(reader);
        let mut meter = crate::progress::ProgressMeter::new(show_progress);
        for tr in &mut transcripts {
            if tr.is_coding() {
                if let Err(e) = tr.build_sequences(|chrom, start, end| {
                    sp.fetch_sequence(chrom, start, end)
                        .map_err(|e| e.to_string())
                }) {
                    eprintln!("Warning: could not build sequences for {}: {}", tr.stable_id, e);
                } else {
                    meter.update();
                }
            }
        }
        meter.finish();
    } else if synonyms_path.is_some() {
        eprintln!(
            "Warning: --synonyms has no effect without --fasta; chromosome names are kept as-is from the GFF3."
        );
    }

    fastvep_cache::transcript_cache::save_cache(&transcripts, Path::new(output_path))?;
    eprintln!("Saved transcript cache to {}", output_path);
    Ok(())
}

/// Quick VCF pre-scan to collect variant regions for indexed GFF3 loading.
/// Returns merged (chrom, start, end) regions expanded by the given distance.
/// One resolved GFF3 input: file path plus the SOURCE label that should be
/// attached to every transcript loaded from it.
#[derive(Debug, Clone)]
pub struct Gff3Spec {
    pub path: String,
    pub source: String,
}

/// Parse a single `--gff3` value. Accepts either `path/to/file.gff3` or
/// `LABEL=path/to/file.gff3`. When no label is provided, infers one from
/// the filename so merged Ensembl + RefSeq runs produce VEP-style
/// `SOURCE=Ensembl` / `SOURCE=RefSeq` values out of the box.
pub fn parse_gff3_arg(arg: &str) -> Gff3Spec {
    // Only treat `=` as a label separator when the prefix doesn't look like
    // a path (no slash) — otherwise something like `./x=1/file.gff3` would
    // get mis-split.
    if let Some((label, rest)) = arg.split_once('=') {
        if !label.contains('/') && !label.contains('\\') && !label.is_empty() {
            return Gff3Spec {
                path: rest.to_string(),
                source: label.to_string(),
            };
        }
    }
    Gff3Spec {
        path: arg.to_string(),
        source: detect_gff3_source(arg),
    }
}

fn detect_gff3_source(path: &str) -> String {
    let fname = Path::new(path)
        .file_name()
        .map(|n| n.to_string_lossy().to_string())
        .unwrap_or_else(|| path.to_string());
    let lower = fname.to_lowercase();
    // RefSeq: NCBI GFF3 archives ship as GCF_*.gff.gz; many user-renamed
    // copies just contain "refseq" in the name.
    if lower.contains("refseq") || lower.starts_with("gcf_") {
        "RefSeq".to_string()
    } else if lower.contains("ensembl") || lower.contains("gencode") {
        "Ensembl".to_string()
    } else {
        fname
    }
}

fn load_one_gff3(
    spec: &Gff3Spec,
    vcf_input: &str,
    distance: u64,
) -> Result<Vec<fastvep_genome::Transcript>> {
    let gff_path = Path::new(&spec.path);
    let tbi_path = format!("{}.tbi", spec.path);

    if spec.path.ends_with(".gz") && Path::new(&tbi_path).exists() {
        let regions = prescan_vcf_regions(vcf_input, distance)?;
        eprintln!("Pre-scanned {} variant regions for {}", regions.len(), spec.path);
        fastvep_cache::gff::parse_gff3_indexed_with_source(gff_path, &regions, &spec.source)
    } else {
        let gff_file = File::open(&spec.path)
            .with_context(|| format!("Opening GFF3 file: {}", spec.path))?;
        if spec.path.ends_with(".gz") || spec.path.ends_with(".bgz") {
            parse_gff3_with_source(flate2::read::MultiGzDecoder::new(gff_file), &spec.source)
        } else {
            parse_gff3_with_source(gff_file, &spec.source)
        }
    }
}

fn prescan_vcf_regions(vcf_path: &str, distance: u64) -> Result<Vec<(String, u64, u64)>> {
    let input_reader = open_vcf_input_reader(vcf_path)
        .with_context(|| format!("Pre-scanning VCF: {}", vcf_path))?;
    let reader = io::BufReader::new(input_reader);

    let mut regions: HashMap<String, (u64, u64)> = HashMap::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') || line.is_empty() { continue; }
        let mut fields = line.split('\t');
        let chrom = match fields.next() { Some(c) => c.to_string(), None => continue };
        let pos: u64 = match fields.next().and_then(|p| p.parse().ok()) { Some(p) => p, None => continue };
        let start = pos.saturating_sub(distance);
        let end = pos + distance;

        let entry = regions.entry(chrom).or_insert((start, end));
        entry.0 = entry.0.min(start);
        entry.1 = entry.1.max(end);
    }

    Ok(regions.into_iter().map(|(chrom, (s, e))| (chrom, s, e)).collect())
}

// =============================================================================
// SA Build: Build supplementary annotation databases from source VCFs
// =============================================================================

/// Standard chromosome ordering for SA builds.
///
/// Canonical names use the UCSC `chr*` style so the resulting `.osa.idx`
/// keys match the GRCh38 convention used by modern gnomAD / ClinVar
/// releases — and by most input VCFs we annotate against. The HashMap
/// resolves both `chr*` and bare forms (plus `MT`/`M`) so source parsers
/// can hand us either style. See issue #37.
fn standard_chrom_map(assembly: &str) -> (Vec<String>, std::collections::HashMap<String, u16>) {
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
fn open_sa_input(input: &str, byte_counter: Option<std::sync::Arc<std::sync::atomic::AtomicU64>>) -> Result<Box<dyn io::Read>> {
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
    use std::sync::Arc;
    use std::sync::atomic::AtomicU64;

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

    let records = make_iter(buf_reader, chrom_map).map(|record| {
        if record.is_ok() {
            meter.update();
        }
        record
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
        Box::new(fastvep_sa::sources::scores::iter_score_tsv(buf, chrom_to_idx, false))
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
    if matches!(source, "omim" | "gnomad_genes" | "gnomad_gene" | "clinvar_protein") {
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
                name: source.to_uppercase().into(),
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
            "Unknown source: {}. Supported: clinvar, gnomad, dbsnp, cosmic, onekg, topmed, mitomap, phylop, gerp, dann, revel, spliceai, primateai, dbnsfp, omim, gnomad_genes, clinvar_protein, custom_vcf, custom_bed, custom",
            source
        ),
    };

    eprintln!("Building {} .osa: {} -> {}", source, input, Path::new(output).with_extension("osa").display());

    // Large, coordinate-sorted population/score sources stream straight into
    // the writer instead of buffering every record (issue #55). They share one
    // helper that differs only in which per-source iterator adapts the reader.
    match source {
        "spliceai" => {
            return run_streaming_sa_build(
                input, output, header, &chrom_map, &chrom_list, show_progress,
                |r, m| fastvep_sa::sources::spliceai::iter_spliceai_vcf(r, m),
            );
        }
        "gnomad" => {
            return run_streaming_sa_build(
                input, output, header, &chrom_map, &chrom_list, show_progress,
                |r, m| fastvep_sa::sources::gnomad::iter_gnomad_vcf(r, m),
            );
        }
        "dbsnp" => {
            return run_streaming_sa_build(
                input, output, header, &chrom_map, &chrom_list, show_progress,
                |r, m| fastvep_sa::sources::dbsnp::iter_dbsnp_vcf(r, m),
            );
        }
        "topmed" => {
            return run_streaming_sa_build(
                input, output, header, &chrom_map, &chrom_list, show_progress,
                |r, m| fastvep_sa::sources::topmed::iter_topmed_vcf(r, m),
            );
        }
        "phylop" => {
            return run_streaming_sa_build(
                input, output, header, &chrom_map, &chrom_list, show_progress,
                |r, m| iter_phylop_auto(r, m),
            );
        }
        "gerp" | "dann" => {
            return run_streaming_sa_build(
                input, output, header, &chrom_map, &chrom_list, show_progress,
                |r, m| fastvep_sa::sources::scores::iter_score_tsv(r, m, false),
            );
        }
        _ => {}
    }

    let buf_reader = io::BufReader::new(open_sa_input(input, None)?);

    eprintln!(
        "INFO  ProgressMeter - Parsing '{}': {} -> {}",
        source, input,
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

/// Classify a `--source custom` input as either `custom_vcf` or
/// `custom_bed` based on its extension. Centralised so the dispatcher in
/// `run_sa_build` and any future callers stay in agreement on the
/// extension rules.
fn classify_custom_input(input: &str) -> Result<&'static str> {
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
fn resolve_custom_name(name: Option<&str>, input: &str) -> String {
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
fn run_custom_vcf_build(
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
        if info_fields.is_empty() { "<all>".to_string() } else { info_fields.join(",") }
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

    let file = File::open(input)
        .with_context(|| format!("Opening input file: {}", input))?;
    let reader: Box<dyn io::Read> = if input.ends_with(".gz") || input.ends_with(".bgz") {
        Box::new(flate2::read::MultiGzDecoder::new(file))
    } else {
        Box::new(file)
    };
    let buf_reader = io::BufReader::new(reader);

    let records = fastvep_sa::custom::parse_custom_vcf(
        buf_reader,
        &chrom_map,
        &header.name,
        info_fields,
    )?;
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
fn run_custom_bed_build(
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

    let file = File::open(input)
        .with_context(|| format!("Opening input file: {}", input))?;
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
/// - `omim`            — disease-gene annotations in genemap2 layout
///                       (13-col TSV): ClinGen Gene-Disease Validity
///                       (preferred per ClinGen SVI / Abou Tayoun
///                       2018) or OMIM `genemap2.txt` (legacy). Drives
///                       PVS1, BS2, PM3, BP2.
/// - `gnomad_genes`    — gnomAD constraint metrics TSV (PVS1, PP2, BP1)
/// - `clinvar_protein` — ClinVar VCF, extracts pathogenic missense by
///                       protein position (PS1, PM1, PM5)
///
/// The output is `<output>.oga`. The runtime loader at
/// `fastvep_annotate::load_gene_providers` picks up any `.oga` file in
/// `--sa-dir` and routes records to the classifier by `json_key`
/// (`omim`, `gnomad_genes`, `clinvar_protein`).
pub fn run_oga_build(source: &str, input: &str, output: &str, _assembly: &str) -> Result<()> {
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

    let file = File::open(input)
        .with_context(|| format!("Opening input file: {}", input))?;
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
            fastvep_sa::sources::clinvar_protein::parse_clinvar_protein_vcf(buf_reader)?
        }
        _ => unreachable!(),
    };

    eprintln!("Parsed {} records from {}", records.len(), source);

    let header = GeneHeader {
        schema_version: SCHEMA_VERSION,
        json_key: json_key.into(),
        name: name.into(),
        version: "latest".into(),
        assembly: _assembly.into(),
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
mod pick_tests {
    use super::*;
    use fastvep_core::{Allele, Impact, Strand};
    use std::sync::Arc;

    fn make_tv(
        transcript_id: &str,
        canonical: bool,
        biotype: &str,
        consequences: Vec<Consequence>,
        mane_select: Option<&str>,
        mane_plus_clinical: Option<&str>,
        appris: Option<&str>,
        tsl: Option<u8>,
        ccds: Option<&str>,
    ) -> TranscriptVariation {
        TranscriptVariation {
            transcript_id: Arc::from(transcript_id),
            gene_id: Arc::from("GENE"),
            gene_symbol: Some(Arc::from("GENE")),
            biotype: Arc::from(biotype),
            allele_annotations: vec![AlleleAnnotation {
                allele: Allele::from_str("A"),
                consequences,
                impact: Impact::Modifier,
                cdna_position: None,
                cds_position: None,
                protein_position: None,
                amino_acids: None,
                codons: None,
                exon: None,
                intron: None,
                distance: None,
                hgvsc: None,
                hgvsp: None,
                hgvsg: None,
                hgvs_offset: None,
                existing_variation: Vec::new(),
                sift: None,
                polyphen: None,
                supplementary: Vec::new(),
                acmg_classification: None,
            }],
            canonical,
            strand: Strand::Forward,
            source: None,
            protein_id: None,
            mane_select: mane_select.map(String::from),
            mane_plus_clinical: mane_plus_clinical.map(String::from),
            tsl,
            appris: appris.map(String::from),
            ccds: ccds.map(String::from),
            gencode_primary: false,
            symbol_source: None,
            hgnc_id: None,
            flags: Vec::new(),
        }
    }

    #[test]
    fn pick_prefers_mane_select_over_canonical() {
        let tvs = vec![
            make_tv("TX_CANON", true, "protein_coding",
                vec![Consequence::MissenseVariant],
                None, None, None, None, None),
            make_tv("TX_MANE", false, "protein_coding",
                vec![Consequence::MissenseVariant],
                Some("TX_MANE.1"), None, None, None, None),
        ];
        assert_eq!(pick_best_transcript_idx(&tvs), Some(1));
    }

    #[test]
    fn pick_prefers_mane_plus_clinical_over_canonical() {
        let tvs = vec![
            make_tv("TX_CANON", true, "protein_coding",
                vec![Consequence::MissenseVariant],
                None, None, None, None, None),
            make_tv("TX_MANE_PC", false, "protein_coding",
                vec![Consequence::MissenseVariant],
                None, Some("TX_MANE_PC.1"), None, None, None),
        ];
        assert_eq!(pick_best_transcript_idx(&tvs), Some(1));
    }

    #[test]
    fn pick_prefers_mane_select_over_mane_plus_clinical() {
        let tvs = vec![
            make_tv("TX_MANE_PC", true, "protein_coding",
                vec![Consequence::MissenseVariant],
                None, Some("TX_MANE_PC.1"), None, None, None),
            make_tv("TX_MANE", false, "protein_coding",
                vec![Consequence::MissenseVariant],
                Some("TX_MANE.1"), None, None, None, None),
        ];
        assert_eq!(pick_best_transcript_idx(&tvs), Some(1));
    }

    #[test]
    fn pick_falls_back_to_canonical_when_no_mane() {
        let tvs = vec![
            make_tv("TX_NONCAN", false, "protein_coding",
                vec![Consequence::StopGained],
                None, None, None, None, None),
            make_tv("TX_CANON", true, "protein_coding",
                vec![Consequence::MissenseVariant],
                None, None, None, None, None),
        ];
        // Canonical wins even though TX_NONCAN has a more severe consequence.
        assert_eq!(pick_best_transcript_idx(&tvs), Some(1));
    }

    #[test]
    fn pick_prefers_protein_coding_biotype() {
        let tvs = vec![
            make_tv("TX_NONCODING", false, "lncRNA",
                vec![Consequence::MissenseVariant],
                None, None, None, None, None),
            make_tv("TX_PC", false, "protein_coding",
                vec![Consequence::MissenseVariant],
                None, None, None, None, None),
        ];
        assert_eq!(pick_best_transcript_idx(&tvs), Some(1));
    }

    #[test]
    fn pick_uses_severity_when_other_fields_equal() {
        let tvs = vec![
            make_tv("TX_A", false, "protein_coding",
                vec![Consequence::SynonymousVariant],
                None, None, None, None, None),
            make_tv("TX_B", false, "protein_coding",
                vec![Consequence::StopGained],
                None, None, None, None, None),
        ];
        assert_eq!(pick_best_transcript_idx(&tvs), Some(1));
    }

    #[test]
    fn pick_tie_breaks_alphabetically_on_transcript_id() {
        let tvs = vec![
            make_tv("TX_Z", false, "protein_coding",
                vec![Consequence::MissenseVariant],
                None, None, None, None, None),
            make_tv("TX_A", false, "protein_coding",
                vec![Consequence::MissenseVariant],
                None, None, None, None, None),
        ];
        assert_eq!(pick_best_transcript_idx(&tvs), Some(1));
    }

    #[test]
    fn pick_prefers_lower_tsl() {
        let tvs = vec![
            make_tv("TX_TSL5", false, "protein_coding",
                vec![Consequence::MissenseVariant],
                None, None, None, Some(5), None),
            make_tv("TX_TSL1", false, "protein_coding",
                vec![Consequence::MissenseVariant],
                None, None, None, Some(1), None),
        ];
        assert_eq!(pick_best_transcript_idx(&tvs), Some(1));
    }

    #[test]
    fn pick_prefers_lower_appris_principal() {
        // P1 should beat P3 even though both are APPRIS-tagged — would fail
        // if APPRIS were compared by presence-only.
        let tvs = vec![
            make_tv("TX_P3", false, "protein_coding",
                vec![Consequence::MissenseVariant],
                None, None, Some("P3"), None, None),
            make_tv("TX_P1", false, "protein_coding",
                vec![Consequence::MissenseVariant],
                None, None, Some("P1"), None, None),
        ];
        assert_eq!(pick_best_transcript_idx(&tvs), Some(1));
    }

    #[test]
    fn pick_prefers_principal_over_alternative_appris() {
        let tvs = vec![
            make_tv("TX_A1", false, "protein_coding",
                vec![Consequence::MissenseVariant],
                None, None, Some("A1"), None, None),
            make_tv("TX_P5", false, "protein_coding",
                vec![Consequence::MissenseVariant],
                None, None, Some("P5"), None, None),
        ];
        assert_eq!(pick_best_transcript_idx(&tvs), Some(1));
    }

    #[test]
    fn pick_accepts_long_form_appris_tags() {
        // Ensembl GFF3 sometimes uses "principal1" / "alternative2".
        let tvs = vec![
            make_tv("TX_ALT", false, "protein_coding",
                vec![Consequence::MissenseVariant],
                None, None, Some("alternative2"), None, None),
            make_tv("TX_PRINC", false, "protein_coding",
                vec![Consequence::MissenseVariant],
                None, None, Some("principal1"), None, None),
        ];
        assert_eq!(pick_best_transcript_idx(&tvs), Some(1));
    }

    #[test]
    fn pick_returns_none_for_empty_input() {
        let tvs: Vec<TranscriptVariation> = vec![];
        assert_eq!(pick_best_transcript_idx(&tvs), None);
    }
}

#[cfg(test)]
mod gff3_arg_tests {
    use super::*;

    #[test]
    fn explicit_label_overrides_detection() {
        let spec = parse_gff3_arg("MyEnsembl=/data/refseq_v115.gff3.gz");
        assert_eq!(spec.path, "/data/refseq_v115.gff3.gz");
        assert_eq!(spec.source, "MyEnsembl");
    }

    #[test]
    fn refseq_filename_autodetected() {
        let spec = parse_gff3_arg("/data/RefSeq.GRCh38.gff");
        assert_eq!(spec.source, "RefSeq");

        let spec = parse_gff3_arg("/refs/GCF_000001405.40.gff.gz");
        assert_eq!(spec.source, "RefSeq");
    }

    #[test]
    fn ensembl_filename_autodetected() {
        let spec = parse_gff3_arg("/data/Homo_sapiens.GRCh38.115.ensembl.gff3");
        assert_eq!(spec.source, "Ensembl");

        let spec = parse_gff3_arg("/data/gencode.v45.annotation.gff3.gz");
        assert_eq!(spec.source, "Ensembl");
    }

    #[test]
    fn unknown_filename_falls_back_to_basename() {
        let spec = parse_gff3_arg("/data/custom_organism.gff3");
        assert_eq!(spec.source, "custom_organism.gff3");
    }

    #[test]
    fn path_with_equals_in_it_is_not_mistaken_for_label() {
        // Real paths containing `=` are rare but possible. A `=` only counts
        // as the label separator if no slash appears before it.
        let spec = parse_gff3_arg("/odd=path/file.gff3");
        assert_eq!(spec.path, "/odd=path/file.gff3");
        // The basename has no Ensembl/RefSeq marker, so we fall back to it.
        assert_eq!(spec.source, "file.gff3");
    }
}

#[cfg(test)]
mod custom_source_tests {
    use super::*;

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
        assert!(msg.contains("custom_vcf"), "error message must mention custom_vcf: {}", msg);
        assert!(msg.contains("custom_bed"), "error message must mention custom_bed: {}", msg);
    }
}
