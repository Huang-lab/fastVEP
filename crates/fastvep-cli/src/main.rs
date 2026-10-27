use anyhow::Result;
use clap::{Parser, Subcommand};

use fastvep_cli::{pipeline, webserver};

#[derive(Parser)]
#[command(name = "fastvep")]
#[command(about = "fastVEP - A high-performance variant effect predictor")]
#[command(version)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

// `Annotate` carries every CLI flag inline, so it is far larger than the other
// variants. That is the shape clap's derive wants, and exactly one of these is
// ever constructed, at startup, so the size difference costs nothing worth
// boxing a field to avoid.
#[allow(clippy::large_enum_variant)]
#[derive(Subcommand)]
enum Commands {
    /// Annotate variants with predicted consequences
    ///
    /// A flag that Ensembl VEP has is accepted in VEP's spelling as well as with
    /// hyphens: `--af_1kg` is `--af-1kg`, `--check_existing` is `--check-existing`.
    Annotate {
        /// Input file (VCF format). Use "-" for stdin.
        #[arg(short, long)]
        input: String,

        /// Output file. Use "-" for stdout.
        #[arg(short, long, default_value = "-")]
        output: String,

        /// GFF3 annotation file(s) for transcript models. Plain text or gzip
        /// (`.gz` / `.bgz`, or a gzip member under any name). A `.gz` with a
        /// sibling `.tbi` is read through tabix. May be repeated
        /// (`--gff3 a.gff3 --gff3 b.gff3`) or passed as a comma-separated
        /// list to replicate VEP's `--merged` cache (e.g. Ensembl + RefSeq
        /// in one annotation run). Each value may optionally be prefixed
        /// with `LABEL=` to control the SOURCE column for that file
        /// (e.g. `--gff3 Ensembl=ens.gff3 --gff3 RefSeq=refseq.gff3`); if
        /// no label is given, it is auto-detected from the filename
        /// (`refseq` / `gcf_` → RefSeq, `ensembl` / `gencode` → Ensembl,
        /// otherwise the basename is used).
        #[arg(long, num_args = 1.., value_delimiter = ',')]
        gff3: Vec<String>,

        /// Path to FASTA reference file. Plain text or gzip (`.gz` / `.bgz`,
        /// or a gzip member under any name). A `.fai` is used only for an
        /// uncompressed FASTA; a gzipped file is loaded into memory.
        #[arg(long)]
        fasta: Option<String>,

        /// Output format (vcf, tab, json)
        #[arg(long, alias = "output_format", default_value = "vcf")]
        output_format: String,

        /// Turn on VEP's `--everything` set, among the flags fastVEP has:
        /// `--hgvs --symbol --numbers --canonical --biotype --variant-class`, and
        /// with `--cache-dir` the frequency flags too (`--af --af-1kg --af-gnomade
        /// --af-gnomadg --max-af --pubmed`).
        /// VEP's set has no `--hgvsg` or `--show_ref_allele`, and neither does this.
        #[arg(long)]
        everything: bool,

        /// Number of variants to buffer
        #[arg(long, alias = "buffer_size", default_value_t = 5000)]
        buffer_size: usize,

        /// Pick one consequence per variant, by `--pick-order` (not by severity:
        /// severity is that order's last tie-break, as in VEP).
        ///
        /// Keeps the winning transcript's annotation for every alt allele.
        /// VEP's `--pick` emits one entry for the whole record, dropping the
        /// other alts; use `--pick-allele` for that, which agrees with VEP.
        #[arg(long, group = "pick_family")]
        pick: bool,

        /// Pick one consequence per (variant, allele). VEP's `--pick_allele`.
        #[arg(long, alias = "pick_allele", group = "pick_family")]
        pick_allele: bool,

        /// Pick one consequence per (variant, allele, gene).
        /// VEP's `--pick_allele_gene`.
        #[arg(long, alias = "pick_allele_gene", group = "pick_family")]
        pick_allele_gene: bool,

        /// Report every consequence, and mark the one `--pick` would have kept
        /// with `PICK=1`. VEP's `--flag_pick`.
        #[arg(long, alias = "flag_pick", group = "pick_family")]
        flag_pick: bool,

        /// Report every consequence, and mark the one chosen per allele.
        /// VEP's `--flag_pick_allele`.
        #[arg(long, alias = "flag_pick_allele", group = "pick_family")]
        flag_pick_allele: bool,

        /// Report every consequence, and mark the one chosen per allele and
        /// gene. VEP's `--flag_pick_allele_gene`.
        #[arg(long, alias = "flag_pick_allele_gene", group = "pick_family")]
        flag_pick_allele_gene: bool,

        /// Add `SYMBOL`, `SYMBOL_SOURCE` and `HGNC_ID` to tab output. VEP's
        /// `--symbol`. The VCF `CSQ` always has them.
        #[arg(long)]
        symbol: bool,

        /// Add HGVS notations: `HGVSc`, `HGVSp` and `HGVS_OFFSET` as tab columns,
        /// and their values in the VCF and JSON output. VEP's `--hgvs`.
        #[arg(long)]
        hgvs: bool,

        /// Add `HGVSg` (`chr10:g.10414C>A`) to tab and VCF output, and to
        /// JSON's `hgvsg`. VEP's `--hgvsg`.
        ///
        /// Written as VEP writes it: an insertion or deletion at its most 3'
        /// position along the genome, and an insertion that repeats the bases
        /// beside it as a `dup`. That reads the reference, so give `--fasta`;
        /// without one the description is the plain one, where the VCF put the
        /// variant. See docs/VEP_DIVERGENCE.md.
        #[arg(long)]
        hgvsg: bool,

        /// Add `CANONICAL` (`YES`) to tab output. VEP's `--canonical`. The VCF
        /// `CSQ` always has it.
        #[arg(long)]
        canonical: bool,

        /// Add `BIOTYPE` to tab output. VEP's `--biotype`. The VCF `CSQ`
        /// always has it.
        #[arg(long)]
        biotype: bool,

        /// Add `EXON` and `INTRON` (`3/10`) to tab output. VEP's `--numbers`.
        /// The VCF `CSQ` always has them.
        #[arg(long)]
        numbers: bool,

        /// Add `VARIANT_CLASS` (`SNV`, `insertion`, `deletion`, ...) to tab and
        /// VCF output. VEP's `--variant_class`.
        #[arg(long, alias = "variant_class")]
        variant_class: bool,

        /// Add `REF_ALLELE` to tab output. VEP's `--show_ref_allele`. The VCF
        /// `CSQ` always has it. (`--explicit-alleles` puts a `REF` column right
        /// after `Allele` instead.)
        #[arg(long, alias = "show_ref_allele")]
        show_ref_allele: bool,

        /// Write the length each position is counted in after `cDNA_position`,
        /// `CDS_position` and `Protein_position` (`14/203`), as VEP's
        /// `--total_length` does. The cDNA length goes on every row whose change
        /// touches an exon, the CDS and protein lengths on those that touch the
        /// coding sequence too, even when the cell itself is `-` (`-/1403`).
        /// Tab and VCF output; JSON writes its positions as numbers.
        #[arg(long, alias = "total_length")]
        total_length: bool,

        /// Do not percent-escape HGVS strings, as VEP's `--no_escape`: a synonymous
        /// `HGVSp` ends `=` and not `%3D`, and the space of a gene symbol in VCF
        /// output is `_` and not `%20`.
        #[arg(long, alias = "no_escape")]
        no_escape: bool,

        /// Accepted for the sake of a `vep` command line, and does nothing:
        /// `fastvep annotate` always writes over its output file, which is what
        /// VEP's `--force_overwrite` asks for.
        #[arg(long, alias = "force_overwrite")]
        force_overwrite: bool,

        /// Leave out RefSeq's predicted transcripts, `XM_` and `XR_`. VEP's
        /// `--exclude_predicted`, which is `not stable_id match ^X._`, and so is this.
        #[arg(long, alias = "exclude_predicted")]
        exclude_predicted: bool,

        /// Leave out transcripts the annotation tags `readthrough_transcript`, the ones that
        /// span two genes. VEP's Ensembl cache has none of them (it is built without them);
        /// a GFF3 has them, and VEP's `--gff` annotates them.
        #[arg(long, alias = "exclude_readthrough")]
        exclude_readthrough: bool,

        /// Leave out transcripts of biotype `artifact`. VEP's Ensembl cache has none of them
        /// (it is built without them); a GFF3 has them, and VEP's `--gff` annotates them.
        #[arg(long, alias = "exclude_artifact")]
        exclude_artifact: bool,

        /// Add `AF`, the 1000 Genomes combined allele frequency of the known
        /// variant, to the tab output and the CSQ. VEP's `--af`. Reads the
        /// variation cache, so it needs `--cache-dir`.
        ///
        /// Written as VEP writes it, per input allele and to four decimals.
        #[arg(long)]
        af: bool,

        /// Add `AF`, `AFR_AF`, `AMR_AF`, `EAS_AF`, `EUR_AF` and `SAS_AF`, the 1000
        /// Genomes frequencies. VEP's `--af_1kg`. Needs `--cache-dir`.
        #[arg(long, alias = "af_1kg")]
        af_1kg: bool,

        /// Add the gnomAD exomes frequencies, `gnomADe_AF` and one column for each
        /// of its populations. VEP's `--af_gnomade`. Needs `--cache-dir`.
        #[arg(long, alias = "af_gnomade", alias = "af-gnomad", alias = "af_gnomad")]
        af_gnomade: bool,

        /// Add the gnomAD genomes frequencies, `gnomADg_AF` and one column for
        /// each of its populations. VEP's `--af_gnomadg`. Needs `--cache-dir`.
        #[arg(long, alias = "af_gnomadg")]
        af_gnomadg: bool,

        /// Add `MAX_AF`, the highest frequency of any 1000 Genomes or gnomAD
        /// population (not the combined ones), and `MAX_AF_POPS`, the populations
        /// that have it. VEP's `--max_af`. Needs `--cache-dir`.
        #[arg(long, alias = "max_af")]
        max_af: bool,

        /// Add `CLIN_SIG`, `SOMATIC` and `PHENO` from the known variants of
        /// `--cache-dir`, to the tab output and the CSQ, as VEP's
        /// `--check_existing` does. Every `--af*` flag, `--max-af` and
        /// `--pubmed` imply it, and so does `--cache-dir`: a variation cache is
        /// all it is for here, where VEP's `--cache` is also the transcript
        /// cache, so the flag is only needed to be typed as VEP has it.
        #[arg(long, alias = "check_existing")]
        check_existing: bool,

        /// Add `PUBMED`, the publications that cite the known variants, and fill
        /// the CSQ's `PUBMED` with VEP's value. VEP's `--pubmed`. Needs
        /// `--cache-dir`.
        #[arg(long)]
        pubmed: bool,

        /// Whether `CLIN_SIG` is the significance a known variant gives for the
        /// input allele (1, VEP's default) or the list it gives for the variant
        /// (0). VEP's `--clin_sig_allele`.
        #[arg(
            long, alias = "clin_sig_allele",
            value_name = "0|1",
            default_value_t = 1,
            value_parser = clap::value_parser!(u8).range(0..=1)
        )]
        clin_sig_allele: u8,

        /// Upstream/downstream distance (bp)
        #[arg(long, default_value_t = 5000)]
        distance: u64,

        /// VEP's tabix variation cache, for `Existing_variation`, `CLIN_SIG`,
        /// `SOMATIC` and `PHENO` (VEP's `--check_existing`; the frequencies and
        /// the PubMed ids need `--af*`, `--max-af` and `--pubmed`): a directory
        /// of `{chr}/all_vars.gz` files, each with a `.csi` or `.tbi` index,
        /// which is what `convert_cache.pl` and the indexed-cache download
        /// provide.
        ///
        /// Only known variants come from here. Transcripts do not: gene models
        /// are read from `--gff3`, and VEP's `{start}-{end}.gz` transcript
        /// slices are never read. A directory with no variation file, or a file
        /// with no index, stops the run rather than annotating nothing.
        ///
        /// Contig names are matched through a `chr_synonyms.txt` in the
        /// directory, as VEP does, and through the `chr` and mitochondrial
        /// naming rules (`chrM`, `MT`).
        #[arg(long)]
        cache_dir: Option<String>,

        /// Path to binary transcript cache file (auto-generated if not specified)
        #[arg(long)]
        transcript_cache: Option<String>,

        /// Directory containing supplementary annotation files (.osa, .osi, .oga)
        #[arg(long)]
        sa_dir: Option<String>,

        /// Skip the default 49-field CSQ annotation pipeline (transcript
        /// consequence, HGVS, ACMG, VEP variation cache) and emit only
        /// supplementary annotations from --sa-dir. Requires --sa-dir.
        #[arg(long)]
        sa_only: bool,

        /// Enable ACMG-AMP variant classification
        #[arg(long)]
        acmg: bool,

        /// Path to ACMG configuration file (TOML) for custom thresholds
        #[arg(long)]
        acmg_config: Option<String>,

        /// Order of criteria used by --pick, VEP's `--pick_order` syntax.
        ///
        /// Default: mane_select,mane_plus_clinical,canonical,appris,tsl,biotype,ccds,rank
        /// (identical to VEP, where transcript status outranks consequence severity).
        /// For clinical reporting put `rank` first, so a variant that disrupts a
        /// transcript is never reported as merely neighbouring a MANE one.
        #[arg(long, alias = "pick_order")]
        pick_order: Option<String>,

        /// Path to a curated functional-evidence TSV supplying PS3/BS3.
        ///
        /// Columns: chrom, pos, ref, alt, criterion (PS3|BS3), strength, pmid, note.
        /// An entry also suppresses PP3/BP4/BP7 for that variant, since ClinGen SVI
        /// ranks experimental evidence above computational prediction.
        #[arg(long)]
        functional_evidence: Option<String>,

        /// Proband sample name for trio analysis (enables de novo / compound-het detection)
        #[arg(long)]
        proband: Option<String>,

        /// Mother sample name for trio analysis
        #[arg(long)]
        mother: Option<String>,

        /// Father sample name for trio analysis
        #[arg(long)]
        father: Option<String>,

        /// Path to a gene-panel file (one gene symbol or Ensembl gene ID per
        /// line; `#` comments and blank lines ignored). When set, tab output
        /// keeps only rows whose transcript belongs to a gene in the panel.
        #[arg(long)]
        gene_list: Option<String>,

        /// Add an explicit REF column to tab output (after the Allele/ALT
        /// column) so spreadsheets can see REF/ALT side-by-side without
        /// reparsing the Location string.
        #[arg(long)]
        explicit_alleles: bool,

        /// Path to a QC rules TOML file. When set, tab output gains a
        /// `QC_CLASS` column populated by the first class whose
        /// INFO-field thresholds the variant satisfies (variant-level,
        /// no per-sample parsing).
        #[arg(long)]
        qc_rules: Option<String>,

        /// Suppress periodic progress output
        #[arg(long, alias = "no_progress", default_value_t = false)]
        no_progress: bool,
    },

    /// Launch the web interface for interactive variant annotation
    Web {
        /// Port to listen on
        #[arg(long, default_value_t = 8080)]
        port: u16,

        /// GFF3 annotation file
        #[arg(long)]
        gff3: Option<String>,

        /// Path to FASTA reference file
        #[arg(long)]
        fasta: Option<String>,
    },

    /// Build a binary transcript cache for fast startup
    Cache {
        /// GFF3 annotation file(s). Plain text or gzip, same as `annotate`.
        /// May be repeated or comma-separated to
        /// build a merged cache (Ensembl + RefSeq); each value may be
        /// `LABEL=path` to control the SOURCE column.
        #[arg(long, num_args = 1.., value_delimiter = ',')]
        gff3: Vec<String>,

        /// Path to FASTA reference file (for pre-building sequences).
        /// Gzip is accepted, same as `annotate`.
        #[arg(long)]
        fasta: Option<String>,

        /// VEP-style chromosome synonyms file (e.g. Ensembl `chr_synonyms.txt`):
        /// one line per contig, whitespace-separated equivalent names. Used to
        /// reconcile GFF3 seqids (e.g. RefSeq `NC_000017.11`) with the FASTA's
        /// contig names. Only takes effect together with `--fasta`.
        #[arg(long)]
        synonyms: Option<String>,

        /// Output cache file path
        #[arg(short, long)]
        output: String,

        /// Suppress periodic progress output
        #[arg(long, default_value_t = false)]
        no_progress: bool,
    },

    /// Build a supplementary annotation database (.osa or .osi) from a source file
    SaBuild {
        /// Source type. Known sources (clinvar, gnomad, dbsnp, …) use their
        /// dedicated parsers; `custom_vcf` and `custom_bed` accept any
        /// well-formed VCF/BED file and produce a generic `.osa` or `.osi`
        /// keyed by `--name`. `custom` is an alias that auto-detects VCF vs
        /// BED from the input extension.
        #[arg(long)]
        source: String,

        /// Input file (VCF, typically gzipped)
        #[arg(short, long)]
        input: String,

        /// Output base path (will create .osa and .osa.idx, or .osi for BED)
        #[arg(short, long)]
        output: String,

        /// Genome assembly (e.g., GRCh38)
        #[arg(long, default_value = "GRCh38")]
        assembly: String,

        /// Display + JSON key name for custom_vcf / custom_bed sources.
        /// Optional — when omitted, the name is derived from the input
        /// filename (extensions stripped). Ignored for built-in sources.
        /// Becomes the `json_key` of the resulting database and the
        /// prefix of the column / INFO field on output.
        #[arg(long)]
        name: Option<String>,

        /// Comma-separated list of INFO fields to extract from a custom VCF.
        /// Empty (default) means "include every INFO key found on each record"
        /// — useful for quick exploration, but each record's JSON object will
        /// vary by which INFO keys it carries.
        #[arg(long, value_delimiter = ',')]
        info_fields: Vec<String>,

        /// On-disk format: `auto` (default), `osa` (v1), or `osa2` (v2).
        /// `auto` builds v2 for the sources that support it — smaller and
        /// faster to query at genome scale — and v1 for the rest, so you get
        /// the best format per source automatically. Pass `osa` to force v1
        /// (e.g. for a faster one-time build) or `osa2` to force v2.
        #[arg(long, default_value = "auto")]
        format: String,

        /// Suppress periodic progress output
        #[arg(long, default_value_t = false)]
        no_progress: bool,
    },

    /// Convert an existing v1 `.osa` supplementary annotation database to v2 `.osa2`
    SaConvert {
        /// Input v1 `.osa` file (its `.osa.idx` must sit alongside it)
        #[arg(short, long)]
        input: String,

        /// Output base path (`.osa2` is appended). Defaults to the input path
        /// with its extension replaced, i.e. `foo.osa` -> `foo.osa2`.
        #[arg(short, long)]
        output: Option<String>,

        /// Genomic width of one chunk, as a power of two (1..=20). A chunk is
        /// the unit a query loads, and a conversion stores whole-record JSON
        /// blobs, so the whole chunk's blobs are decompressed to read one
        /// record. The default of 20 (1 Mb) suits the sparse sources - dbSNP,
        /// ClinVar, COSMIC, REVEL - and keeps them compact. A source that
        /// scores nearly every base, SpliceAI above all, puts so much in one
        /// megabase that queries slow to a crawl and can exceed the blob
        /// decompression limit outright; pass 16 (65 kb) for those. Measured on
        /// a SpliceAI-density source, 300 scattered queries went from 9.7s
        /// (with most lookups failing) to 0.5s, for 12% more file size.
        #[arg(long, default_value_t = 20, value_parser = clap::value_parser!(u32).range(1..=20))]
        chunk_bits: u32,

        /// Suppress periodic progress output
        #[arg(long, default_value_t = false)]
        no_progress: bool,
    },

    /// Filter annotated VEP output
    Filter {
        /// Input file (VEP-annotated VCF)
        #[arg(short, long)]
        input: String,

        /// Output file
        #[arg(short, long, default_value = "-")]
        output: String,

        /// Filter expression
        #[arg(long)]
        filter: String,
    },
}

/// `Cli::parse`, except that a flag of Ensembl VEP that fastVEP does not have is told what there
/// is instead of what clap guesses it was mistaken for (`fastvep_cli::vep_flags`).
fn parse_cli() -> Cli {
    use clap::builder::StyledStr;
    use clap::error::{ContextKind, ContextValue, ErrorKind};

    match Cli::try_parse() {
        Ok(cli) => cli,
        Err(mut err) => {
            if err.kind() == ErrorKind::UnknownArgument {
                let tip = match err.get(ContextKind::InvalidArg) {
                    Some(ContextValue::String(arg)) => fastvep_cli::vep_flags::tip(arg),
                    _ => None,
                };
                if let Some(tip) = tip {
                    // The usage line clap prints under the error has the guess in it too
                    // (`Usage: fastvep annotate --input <INPUT> --symbol`), which says that
                    // `--symbol` is what to type instead. It is the usage of a flag that is not
                    // anything's near neighbour, which has no guess in it.
                    err.remove(ContextKind::SuggestedArg);
                    if let Some(usage) = plain_usage(&err) {
                        err.insert(ContextKind::Usage, ContextValue::StyledStr(usage));
                    }
                    err.insert(
                        ContextKind::Suggested,
                        ContextValue::StyledStrs(vec![StyledStr::from(tip)]),
                    );
                }
            }
            err.exit()
        }
    }
}

/// The usage line of an unknown-flag error that has no guess in it, for an error of the same
/// command as `err`.
fn plain_usage(err: &clap::Error) -> Option<clap::builder::StyledStr> {
    use clap::error::{ContextKind, ContextValue};

    let of_annotate = match err.get(ContextKind::Usage) {
        Some(ContextValue::StyledStr(usage)) => usage.to_string().contains("annotate"),
        _ => false,
    };
    if !of_annotate {
        return None;
    }
    let plain = Cli::try_parse_from(["fastvep", "annotate", "--no-such-flag-anywhere"]).err()?;
    match plain.get(ContextKind::Usage) {
        Some(ContextValue::StyledStr(usage)) => Some(usage.clone()),
        _ => None,
    }
}

fn main() -> Result<()> {
    env_logger::init();
    let cli = parse_cli();

    match cli.command {
        Commands::Annotate {
            input,
            output,
            gff3,
            fasta,
            output_format,
            everything,
            buffer_size: _,
            pick,
            pick_allele,
            pick_allele_gene,
            flag_pick,
            flag_pick_allele,
            flag_pick_allele_gene,
            symbol,
            hgvs,
            hgvsg,
            canonical,
            biotype,
            numbers,
            variant_class,
            show_ref_allele,
            total_length,
            no_escape,
            force_overwrite: _,
            exclude_predicted,
            exclude_readthrough,
            exclude_artifact,
            af,
            af_1kg,
            af_gnomade,
            af_gnomadg,
            max_af,
            check_existing,
            pubmed,
            clin_sig_allele,
            distance,
            cache_dir,
            transcript_cache,
            sa_dir,
            sa_only,
            acmg,
            acmg_config,
            pick_order,
            functional_evidence,
            proband,
            mother,
            father,
            gene_list,
            explicit_alleles,
            qc_rules,
            no_progress,
        } => {
            pipeline::run_annotate(pipeline::AnnotateConfig {
                input,
                output,
                gff3,
                fasta,
                output_format,
                pick: pipeline::PickFlags {
                    pick,
                    pick_allele,
                    pick_allele_gene,
                    flag_pick,
                    flag_pick_allele,
                    flag_pick_allele_gene,
                },
                // `--everything` is VEP's set of flags (Config.pm `OPTION_SETS`),
                // restricted to what fastVEP can write.
                hgvs: hgvs || everything,
                columns: pipeline::ColumnFlags {
                    show_ref_allele,
                    variant_class,
                    symbol,
                    biotype,
                    canonical,
                    numbers,
                    hgvsg,
                    known: fastvep_io::known::KnownFlags {
                        af,
                        af_1kg,
                        af_gnomade,
                        af_gnomadg,
                        max_af,
                        check_existing,
                        pubmed,
                        no_clin_sig_allele: clin_sig_allele == 0,
                    },
                }
                .union(if everything {
                    pipeline::ColumnFlags::everything()
                } else {
                    pipeline::ColumnFlags::default()
                })
                // `--everything` has the frequency flags too, and they read a
                // cache: with `--cache-dir` it is VEP's set, and without one it is
                // the set fastVEP has always had, instead of an error for a flag
                // the user did not type.
                .union(pipeline::ColumnFlags {
                    known: if everything && cache_dir.is_some() {
                        fastvep_io::known::KnownFlags::everything()
                    } else {
                        Default::default()
                    },
                    ..Default::default()
                }),
                compat: pipeline::CompatFlags {
                    total_length,
                    no_escape,
                },
                exclude: pipeline::TranscriptExclusions {
                    predicted: exclude_predicted,
                    readthrough: exclude_readthrough,
                    artifact: exclude_artifact,
                },
                distance,
                cache_dir,
                transcript_cache,
                sa_dir,
                sa_only,
                acmg,
                acmg_config,
                pick_order,
                functional_evidence,
                proband,
                mother,
                father,
                gene_list,
                explicit_alleles,
                qc_rules,
                show_progress: !no_progress,
            })?;
        }
        Commands::Cache {
            gff3,
            fasta,
            synonyms,
            output,
            no_progress,
        } => {
            pipeline::run_cache_build(
                &gff3,
                fasta.as_deref(),
                synonyms.as_deref(),
                &output,
                !no_progress,
            )?;
        }
        Commands::Web { port, gff3, fasta } => {
            webserver::run_server(port, gff3, fasta)?;
        }
        Commands::SaBuild {
            source,
            input,
            output,
            assembly,
            name,
            info_fields,
            format,
            no_progress,
        } => pipeline::run_sa_build_format(
            &format,
            &source,
            &input,
            &output,
            &assembly,
            name.as_deref(),
            &info_fields,
            !no_progress,
        )?,
        Commands::SaConvert {
            input,
            output,
            chunk_bits,
            no_progress,
        } => {
            // Default the output alongside the input: `foo.osa` -> `foo.osa2`.
            // `run_sa_convert` appends the extension, so hand it the stem.
            let output = output.unwrap_or_else(|| {
                std::path::Path::new(&input)
                    .with_extension("")
                    .to_string_lossy()
                    .into_owned()
            });
            pipeline::run_sa_convert(&input, &output, chunk_bits, !no_progress)?
        }
        Commands::Filter {
            input,
            output,
            filter,
        } => {
            pipeline::run_filter(&input, &output, &filter)?;
        }
    }

    Ok(())
}
