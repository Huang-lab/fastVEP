//! `fastvep annotate`: the per-variant annotation pipeline.

use super::cache_build::{load_one_gff3, parse_gff3_arg, Gff3Spec};
use super::open_vcf_input_reader;
use super::pick::{parse_pick_order, pick_best_transcript_idx_with, DEFAULT_PICK_ORDER};
use anyhow::{Context, Result};
use fastvep_cache::annotation::{AnnotationProvider, AnnotationValue, GeneAnnotationProvider};
use fastvep_cache::fasta::FastaReader;
use fastvep_cache::info::CacheInfo;
use fastvep_cache::providers::{
    FastaSequenceProvider, IndexedTranscriptProvider, MatchedVariant, SequenceProvider,
    TabixVariationProvider, TranscriptProvider, VariationProvider,
};
use fastvep_cache::transcript_cache::StaleCacheFormat;
use fastvep_consequence::ConsequencePredictor;
use fastvep_core::{Allele, Consequence};
use fastvep_genome::Transcript;
use fastvep_hgvs;
use fastvep_io::output;
use fastvep_io::variant::{AlleleAnnotation, TranscriptVariation, VariationFeature};
use fastvep_io::vcf::VcfParser;
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;
use std::sync::atomic::{AtomicUsize, Ordering};

// Shared annotation utilities from fastvep-annotate (used by batch pipeline).
use fastvep_annotate::{
    annotate_intergenic, annotate_sa_only_scaffold, complement_allele, convert_ins_to_dup,
    convert_ins_to_dup_noncoding, load_gene_providers, load_sa_providers,
    three_prime_shift_intronic, zip_positions,
};

const BATCH_SIZE: usize = 1024;

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
    /// Order of criteria used by `--pick`, in VEP's `--pick_order` syntax.
    /// `None` uses VEP's default order.
    pub pick_order: Option<String>,
    /// Path to a curated functional-evidence TSV supplying PS3/BS3.
    pub functional_evidence: Option<String>,
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

pub fn run_annotate(mut config: AnnotateConfig) -> Result<()> {
    eprintln!("Annotating: {} -> {}", config.input, config.output);
    // Several ACMG criteria read HGVS c. notation rather than raw coordinates:
    // PVS1's canonical ±1/±2 splice gate and BP7's deep-intronic extension both
    // need the intronic offset, and BA1's Ghosh 2018 exception list is keyed by
    // `c.` string. Without `--hgvs` those signals are silently absent and the
    // criteria quietly fall back to weaker behaviour, so `--acmg` turns HGVS on
    // rather than degrading. The output field list is unchanged either way -
    // HGVSc/HGVSp are always in the CSQ header and were simply empty before.
    if config.acmg && !config.hgvs && !config.sa_only {
        config.hgvs = true;
    }
    let config = config;
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
    let gff3_specs: Vec<Gff3Spec> = config.gff3.iter().map(|s| parse_gff3_arg(s)).collect();

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
        config
            .transcript_cache
            .as_ref()
            .map(|p| Path::new(p).to_path_buf())
            .or_else(|| {
                single_gff3.map(|s| {
                    fastvep_cache::transcript_cache::default_cache_path(Path::new(&s.path))
                })
            })
    };

    // Set when any GFF3 source was read through the tabix path, which only
    // returns features overlapping this VCF's variant regions. Such a
    // transcript set is valid for *this* input and no other, so it must never
    // be persisted as a sidecar cache — see the save site below.
    let mut region_restricted = false;

    let mut transcripts = if sa_only {
        Vec::new()
    } else {
        'load: {
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
                            .map(|s| {
                                fastvep_cache::transcript_cache::cache_is_fresh(
                                    cp,
                                    Path::new(&s.path),
                                )
                            })
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
                                eprintln!(
                                    "Loaded {} transcripts from cache {}",
                                    trs.len(),
                                    cp.display()
                                );
                                break 'load trs;
                            }
                            Err(e) => {
                                // An explicit --transcript-cache is the user
                                // naming the transcript source. If it will not
                                // load we cannot quietly substitute something
                                // else: with no --gff3 to fall back to, the run
                                // continues with zero transcripts and calls
                                // every variant intergenic at exit code 0. A
                                // cache truncated by a concurrent writer or a
                                // full disk looked exactly like a small
                                // annotation set.
                                if explicit_cache {
                                    // A cache rejected for its format (pre-#90,
                                    // pre-#98) is intact, so "truncated or
                                    // corrupt" is the wrong diagnosis and would
                                    // send the user looking for a disk problem
                                    // they do not have. Say which of the two it
                                    // is, once, and let StaleCacheFormat supply
                                    // the reason for the format it found.
                                    let why = match e.downcast_ref::<StaleCacheFormat>() {
                                        Some(stale) => stale.to_string(),
                                        None => format!(
                                            "it could not be read ({e}), and is most likely \
                                             truncated or corrupt - delete it and rebuild with \
                                             `fastvep cache`."
                                        ),
                                    };
                                    return Err(anyhow::anyhow!(
                                        "Transcript cache {} cannot be used: {why} Refusing to \
                                         continue, because annotating without it would report \
                                         every variant as intergenic.",
                                        cp.display()
                                    ));
                                }
                                eprintln!(
                                    "Warning: sidecar cache {} could not be loaded ({e}); rebuilding from GFF3",
                                    cp.display()
                                );
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
                // Reachable only when a sidecar cache failed to load and there
                // is no GFF3 to rebuild from (an explicit cache failure has
                // already returned above). Annotating with zero transcripts
                // emits a complete, well-formed, all-intergenic VCF, which is
                // indistinguishable downstream from a real result.
                return Err(anyhow::anyhow!(
                    "No transcript source is usable: no --gff3 was given and no transcript cache \
                     could be loaded. Refusing to continue, because every variant would be \
                     reported as intergenic. Pass --gff3, or --sa-only for a \
                     transcript-free annotation."
                ));
            } else {
                let mut all: Vec<fastvep_genome::Transcript> = Vec::new();
                for spec in &gff3_specs {
                    let (trs, restricted) = load_one_gff3(spec, &config.input, config.distance)?;
                    region_restricted |= restricted;
                    // Say which kind of load this was. A tabix load returns a
                    // small fraction of the file's transcripts, and printing
                    // that count in the same shape as a whole-file count is
                    // what made a by-design optimisation look like data loss.
                    if restricted {
                        eprintln!(
                            "Loaded {} transcripts overlapping this input's variant regions from {} (source label: {}); the file as a whole holds more",
                            trs.len(),
                            spec.path,
                            spec.source
                        );
                    } else {
                        eprintln!(
                            "Loaded {} transcripts from {} (source label: {})",
                            trs.len(),
                            spec.path,
                            spec.source
                        );
                    }
                    all.extend(trs);
                }
                if all.is_empty() {
                    return Err(anyhow::anyhow!(
                    "GFF3 source(s) [{}] produced 0 transcripts — likely malformed, truncated, or unrecognized format. Refusing to continue with empty transcript set.",
                    gff3_specs.iter().map(|s| s.path.as_str()).collect::<Vec<_>>().join(", ")
                ));
                }
                if gff3_specs.len() > 1 {
                    eprintln!(
                        "Merged {} GFF3 sources into {} total transcripts",
                        gff3_specs.len(),
                        all.len()
                    );
                }
                all
            }
        }
    };

    // Load FASTA reference (prefer mmap with .fai index, fall back to in-memory).
    // Skipped in --sa-only mode.
    let seq_provider: Option<Box<dyn SequenceProvider>> = if sa_only {
        None
    } else if let Some(ref fasta_path) = config.fasta {
        let fai_path = format!("{}.fai", fasta_path);
        if Path::new(&fai_path).exists() {
            let reader = fastvep_cache::fasta::MmapFastaReader::open(Path::new(fasta_path))?;
            eprintln!(
                "Memory-mapped reference FASTA from {} (using .fai index)",
                fasta_path
            );
            Some(Box::new(
                fastvep_cache::providers::MmapFastaSequenceProvider::new(reader),
            ))
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
    let needs_seq_build = transcripts
        .iter()
        .any(|t| t.is_coding() && t.spliced_seq.is_none());
    if needs_seq_build {
        if let Some(ref sp) = seq_provider {
            let built = AtomicUsize::new(0);
            transcripts.par_iter_mut().for_each(|tr| {
                if tr.is_coding() && tr.spliced_seq.is_none() {
                    if let Err(e) = tr.build_sequences(|chrom, start, end| {
                        sp.fetch_sequence(chrom, start, end)
                            .map_err(|e| e.to_string())
                    }) {
                        eprintln!(
                            "Warning: could not build sequences for {}: {}",
                            tr.stable_id, e
                        );
                    } else {
                        built.fetch_add(1, Ordering::Relaxed);
                    }
                }
            });
            eprintln!(
                "Built sequences for {} coding transcripts",
                built.load(Ordering::Relaxed)
            );
        }
    }

    // Save cache after sequence build (only if sequences were built or
    // cache doesn't exist). Sidecar-cache writes are gated to the
    // single-GFF3 case - multi-source caches need to round-trip through
    // an explicit `--transcript-cache` path so the on-disk filename
    // isn't tied to one of N inputs.
    //
    // They are also gated on `!region_restricted`. A tabix load only returns
    // features overlapping this VCF's variants, so persisting it would leave a
    // sidecar cache that looks fresh to every later run against the same
    // --gff3 while holding only the first input's neighbourhood. The next run
    // on a different VCF then loads it, finds no transcripts near its own
    // variants, and reports them all as intergenic - exit code 0, no warning.
    // Measured before this gate: 171 of 173 variants in
    // validation/human/vep_example_GRCh38.vcf came back wrong that way,
    // including a missense call reduced to `intergenic_variant`.
    if needs_seq_build && !region_restricted {
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
            cache_info.species,
            cache_info.assembly,
            cache_info.variation_cols.len()
        );
        Some(TabixVariationProvider::new(Path::new(dir), &cache_info)?)
    } else {
        None
    };

    // Load supplementary annotation providers from --sa-dir.
    //
    // Reported on stderr like the other startup steps: this is the one that
    // scales with the number of files in --sa-dir, so a slow shared filesystem
    // shows up here rather than as a silent stall before the progress meter
    // (issue #78). `tracing` has no subscriber installed in the CLI.
    let sa_providers: Vec<Box<dyn AnnotationProvider>> = if let Some(ref dir) = config.sa_dir {
        let t0 = std::time::Instant::now();
        let loaded = load_sa_providers(Path::new(dir))?;
        eprintln!(
            "Loaded {} supplementary annotation source(s) from {} in {:.1}s",
            loaded.len(),
            dir,
            t0.elapsed().as_secs_f64()
        );
        loaded
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
                        .filter(|e| e.path().extension().and_then(|s| s.to_str()) == Some("oga"))
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
        let loaded = load_gene_providers(Path::new(dir))?;
        if !loaded.is_empty() {
            eprintln!(
                "Loaded {} gene-level annotation source(s) from {}",
                loaded.len(),
                dir
            );
        }
        loaded
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
    // Whether gnomAD queries get reference-normalized (see the SA attach loop
    // below) depends only on which providers loaded and whether a FASTA is
    // present — never on the variant. Resolved once here rather than re-scanning
    // every provider for every allele of every variant inside the parallel loop.
    let normalize_gnomad_queries =
        seq_provider.is_some() && sa_providers.iter().any(|sa| sa.json_key() == "gnomad");
    // BP3 must tell "this variant is not in a repeat" from "no repeat database
    // was loaded". An interval source only yields an annotation on a hit, so
    // those two are indistinguishable downstream unless the loaded-source list
    // is consulted here. Same matching the classifier uses to find the track.
    let repeat_db_loaded = sa_providers.iter().any(|sa| {
        let k = sa.json_key().to_lowercase();
        k.contains("repeat") || k.contains("repeatmasker") || k.contains("simple_repeat")
    });
    let gene_json_keys: Vec<String> = gene_providers
        .iter()
        .map(|gp| gp.json_key().to_string())
        .collect();

    // Curated functional-assay evidence for PS3/BS3. Parsed up front so a
    // malformed row fails before the run starts rather than after an hour of
    // annotation, and so the error names the offending line.
    // Resolved before annotation so a bad --pick-order fails immediately with
    // the offending criterion named, rather than after the run.
    let pick_order = match config.pick_order.as_deref() {
        Some(spec) => {
            let order = parse_pick_order(spec)?;
            eprintln!("Using --pick order: {}", spec);
            order
        }
        None => DEFAULT_PICK_ORDER.to_vec(),
    };

    let functional_index = match config.functional_evidence.as_deref() {
        Some(path) => {
            let idx = fastvep_classification::FunctionalEvidenceIndex::from_file(Path::new(path))
                .map_err(anyhow::Error::msg)?;
            eprintln!(
                "Loaded {} curated functional-evidence entries from {}",
                idx.len(),
                path
            );
            Some(idx)
        }
        None => None,
    };
    let functional_evidence = functional_index.as_ref();
    let owned_vcf_info_ids = output::vcf_owned_info_ids(&sa_json_keys, &gene_json_keys);
    let generated_vcf_headers = output::vcf_info_header_lines(
        &sa_json_keys,
        &gene_json_keys,
        output::DEFAULT_CSQ_FIELDS,
        sa_only,
    );
    // Precompute the loaded-source lookup once so the per-row tab writer
    // doesn't redo an O(specs × keys) membership scan for every variant.
    let supplementary_specs = output::LoadedSupplementarySpecs::new(&sa_json_keys, &gene_json_keys);

    // Write headers based on output format
    match config.output_format.as_str() {
        "vcf" => {
            // Pass through original VCF headers
            for header_line in vcf_parser.header_lines() {
                if let Some(info_id) = output::vcf_info_header_id(header_line) {
                    if owned_vcf_info_ids.contains(&info_id) {
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
        let mut batch: Vec<(VariationFeature, HashMap<String, Vec<MatchedVariant>>)> =
            Vec::with_capacity(BATCH_SIZE);
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
                                let matches = vp
                                    .get_matched_variants(
                                        &vf.position.chromosome,
                                        vf.position.start,
                                        vf.position.end,
                                        &ref_str,
                                        &alt_str,
                                    )
                                    .unwrap_or_default();
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
                                vf.existing_variants
                                    .push(fastvep_io::variant::KnownVariant {
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
                let positions = chrom_positions.entry(&vf.position.chromosome).or_default();
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
            for (i, tc) in transcript_consequences.iter().enumerate() {
                // `ConsequencePredictor::predict` emits one consequence per
                // transcript it was handed, in that order, so index `i` is the
                // match and the scan below never runs. The structural-variant
                // branch above makes no such promise, so a miss falls back to a
                // scan rather than losing the transcript. Before this, every
                // consequence rescanned the whole overlap set, and a variant in
                // a gene-dense window - a protocadherin or MHC cluster, where
                // `--distance` pulls in dozens of transcripts - paid for that
                // lookup in proportion to the square of the overlap.
                let transcript: Option<&Transcript> = overlapping
                    .get(i)
                    .copied()
                    .filter(|t| t.stable_id == tc.transcript_id)
                    .or_else(|| {
                        overlapping
                            .iter()
                            .copied()
                            .find(|t| t.stable_id == tc.transcript_id)
                    });

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
                            protein_position: ac.protein_range(),
                            amino_acids: ac.amino_acids.clone(),
                            codons: ac.codons.clone(),
                            exon: ac.exon,
                            intron: ac.intron,
                            distance: ac.distance,
                            protein_length: ac.protein_length,
                            escapes_nmd: ac.escapes_nmd,
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
                                                // Extract from CDS start to end of spliced seq (includes 3'UTR).
                                                // Guard against malformed/truncated GFF3-derived transcript data
                                                // where `coding_start` is inconsistent with the actual spliced
                                                // sequence length — skip HGVSp generation for this case rather
                                                // than panicking on an out-of-bounds slice.
                                                let coding_start_idx = (coding_start - 1) as usize;
                                                if coding_start >= 1 && coding_start_idx <= spliced.len() {
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
                                                let fs_codon_table =
                                                    if fastvep_genome::is_mitochondrial(&tr.chromosome) {
                                                        fastvep_genome::mitochondrial_codon_table()
                                                    } else {
                                                        fastvep_genome::CodonTable::standard()
                                                    };
                                                ann.hgvsp = fastvep_hgvs::hgvsp_frameshift(
                                                    &versioned_pid,
                                                    ref_from_cds,
                                                    &alt_from_cds,
                                                    codon_start,
                                                    &fs_codon_table,
                                                );
                                                }
                                            }
                                        } else if aa.1 == "-"
                                            || ac.consequences.contains(&Consequence::InframeDeletion)
                                            || ac.consequences.contains(&Consequence::InframeInsertion)
                                        {
                                            // In-frame indel / delins (frameshift handled
                                            // above). aa.0 holds the replaced residues, aa.1 the
                                            // replacement ("-" for a pure deletion).
                                            //
                                            // Insertions must route here too, not to `hgvsp()`:
                                            // that compares only the first residue of each side,
                                            // so `W/WR` reads as unchanged and renders
                                            // `p.Trp185=` for a variant that lengthens the
                                            // protein.
                                            //
                                            // The peptide lets `hgvsp_inframe_indel` apply the
                                            // HGVS 3'-rule and collapse a repeat to `dup`,
                                            // matching Ensembl VEP; it degrades to the
                                            // unshifted description without one.
                                            // `tr.strand` says which end of
                                            // `aa.0` the `ps` above names: on
                                            // the reverse strand a shrinking
                                            // change arrives anchored at the
                                            // end of its span (#89, #96).
                                            ann.hgvsp = fastvep_hgvs::hgvsp_inframe_indel(
                                                &versioned_pid,
                                                ps,
                                                &aa.0,
                                                &aa.1,
                                                tr.peptide.as_deref().map(str::as_bytes),
                                                tr.strand,
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
                if let Some(idx) =
                    pick_best_transcript_idx_with(&vf.transcript_variations, &pick_order)
                {
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
                let sa_queries = vf.query_alleles();
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
                        let gnomad_norm = if normalize_gnomad_queries {
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

                let functional_by_alt = resolve_functional_by_alt(functional_evidence, vf);
                // Resolved here rather than inside the classifier: the ClinVar
                // splice index is keyed by genomic coordinate, and
                // `extract_classification_input` is handed a transcript-level
                // view that carries none.
                let query_alleles = vf.query_alleles();

                for tv in &mut vf.transcript_variations {
                    let gene_sym = tv.gene_symbol.as_deref().unwrap_or("");
                    let gene_anns: Vec<&fastvep_core::GeneAnnotation> =
                        vf.gene_annotations
                            .iter()
                            .filter(|ga| ga.gene_symbol == gene_sym)
                            .collect();
                    for aa in &mut tv.allele_annotations {
                        let alt_idx = vf.alt_alleles.iter().position(|a| *a == aa.allele);
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
                                aa.protein_length,
                                aa.escapes_nmd,
                                repeat_db_loaded,
                                fastvep_annotate::splice_ps1_evidence(aa, &gene_anns, &query_alleles, alt_idx),
                                alt_idx.and_then(|i| query_alleles.get(i)).map(|(_, pos, r, a)| {
                                    (vf.position.chromosome.to_string(), *pos, r.clone(), a.clone())
                                }),
                                fastvep_classification::is_pure_insertion(&vf.ref_allele),
                                alt_idx.and_then(|i| functional_by_alt[i].clone()),
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
                enrich_compound_het_batch(
                    &mut vfs,
                    acmg_cfg,
                    &sample_names,
                    functional_evidence,
                    repeat_db_loaded,
                );
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
                    for line in
                        output::format_tab_line_with(vf, &supplementary_specs, sa_only, opts)
                    {
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

    let samples = fastvep_io::sample::parse_samples(format_str, &sample_strs, sample_names);

    let proband_gt = samples
        .iter()
        .find(|s| s.name == trio.proband)
        .map(sample_data_to_genotype_info_cli);

    let mother_gt = trio.mother.as_ref().and_then(|name| {
        samples
            .iter()
            .find(|s| &s.name == name)
            .map(sample_data_to_genotype_info_cli)
    });

    let father_gt = trio.father.as_ref().and_then(|name| {
        samples
            .iter()
            .find(|s| &s.name == name)
            .map(sample_data_to_genotype_info_cli)
    });

    (proband_gt, mother_gt, father_gt)
}

/// Convert a SampleData to GenotypeInfo (CLI path).
fn sample_data_to_genotype_info_cli(
    sample: &fastvep_io::sample::SampleData,
) -> fastvep_classification::GenotypeInfo {
    let gt = sample.genotype.as_ref();
    let is_het = gt.is_some_and(|g| g.is_het());
    let is_hom_ref = gt.is_some_and(|g| g.is_hom_ref());
    let is_hom_alt = gt.is_some_and(|g| g.is_hom_alt());
    let is_missing = gt.is_none_or(|g| g.is_missing());
    let is_phased = gt.is_some_and(|g| g.phased);

    let alt_allele_index = gt.and_then(|g| g.alleles.iter().filter_map(|a| *a).find(|&a| a > 0));

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
    functional_evidence: Option<&fastvep_classification::FunctionalEvidenceIndex>,
    repeat_db_loaded: bool,
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

                let proband_het = proband_gt.as_ref().is_some_and(|g| g.is_het);
                let is_phased = proband_gt.as_ref().is_some_and(|g| g.is_phased);
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

    for gene_infos in gene_variants.values() {
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
                                .is_some_and(|a| a.is_some_and(|v| v > 0));
                            let other_alt_on_first = other
                                .proband_alleles
                                .first()
                                .is_some_and(|a| a.is_some_and(|v| v > 0));
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
            let functional_by_alt = resolve_functional_by_alt(functional_evidence, vf);
            let query_alleles = vf.query_alleles();
            let alt_idx = vf.alt_alleles.iter().position(|a| *a == aa.allele);

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
                aa.protein_length,
                aa.escapes_nmd,
                repeat_db_loaded,
                fastvep_annotate::splice_ps1_evidence(aa, &gene_anns, &query_alleles, alt_idx),
                alt_idx
                    .and_then(|i| query_alleles.get(i))
                    .map(|(_, pos, r, a)| {
                        (
                            vf.position.chromosome.to_string(),
                            *pos,
                            r.clone(),
                            a.clone(),
                        )
                    }),
                fastvep_classification::is_pure_insertion(&vf.ref_allele),
                alt_idx.and_then(|i| functional_by_alt[i].clone()),
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
            fields.chrom,
            fields.pos,
            fields.id,
            fields.ref_allele,
            fields.alt,
            fields.qual,
            fields.filter,
            info
        )?;

        for rest_field in &fields.rest {
            write!(writer, "\t{}", rest_field)?;
        }
        writeln!(writer)?;
    }

    Ok(())
}

/// Curated functional evidence for each ALT allele of a record, if the run was
/// given a `--functional-evidence` file that names them.
///
/// Resolved for the whole record at once, before the mutable walk over
/// `transcript_variations` begins, because that walk borrows the record and a
/// per-allele lookup inside it would need the record again.
///
/// Keyed on the record's original VCF coordinates rather than fastVEP's
/// normalised alleles: the curated file is written by a human reading a VCF, so
/// that is the form the entry will be in. The result is positional, one slot
/// per ALT, so one allele's curated result is never applied to its neighbours.
fn resolve_functional_by_alt(
    index: Option<&fastvep_classification::FunctionalEvidenceIndex>,
    vf: &VariationFeature,
) -> Vec<Option<fastvep_classification::FunctionalEvidence>> {
    let empty = || vec![None; vf.alt_alleles.len()];
    let (Some(index), Some(vcf)) = (index, vf.vcf_fields.as_ref()) else {
        return empty();
    };
    (0..vf.alt_alleles.len())
        .map(|i| {
            index
                .for_vcf_allele(&vcf.chrom, vcf.pos, &vcf.ref_allele, &vcf.alt, i)
                .cloned()
        })
        .collect()
}
