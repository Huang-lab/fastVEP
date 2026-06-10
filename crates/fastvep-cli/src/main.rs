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

#[derive(Subcommand)]
enum Commands {
    /// Annotate variants with predicted consequences
    Annotate {
        /// Input file (VCF format). Use "-" for stdin.
        #[arg(short, long)]
        input: String,

        /// Output file. Use "-" for stdout.
        #[arg(short, long, default_value = "-")]
        output: String,

        /// GFF3 annotation file(s) for transcript models. May be repeated
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

        /// Path to FASTA reference file
        #[arg(long)]
        fasta: Option<String>,

        /// Output format (vcf, tab, json)
        #[arg(long, default_value = "vcf")]
        output_format: String,

        /// Turn on all common annotation flags
        #[arg(long)]
        everything: bool,

        /// Number of variants to buffer
        #[arg(long, default_value_t = 5000)]
        buffer_size: usize,

        /// Pick one consequence per variant (most severe)
        #[arg(long)]
        pick: bool,

        /// Include gene symbol in output
        #[arg(long)]
        symbol: bool,

        /// Include HGVS notations
        #[arg(long)]
        hgvs: bool,

        /// Include canonical transcript flag
        #[arg(long)]
        canonical: bool,

        /// Upstream/downstream distance (bp)
        #[arg(long, default_value_t = 5000)]
        distance: u64,

        /// Path to VEP cache directory for known variant annotation
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
        /// GFF3 annotation file(s). May be repeated or comma-separated to
        /// build a merged cache (Ensembl + RefSeq); each value may be
        /// `LABEL=path` to control the SOURCE column.
        #[arg(long, num_args = 1.., value_delimiter = ',')]
        gff3: Vec<String>,

        /// Path to FASTA reference file (for pre-building sequences)
        #[arg(long)]
        fasta: Option<String>,

        /// Output cache file path
        #[arg(short, long)]
        output: String,
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

fn main() -> Result<()> {
    env_logger::init();
    let cli = Cli::parse();

    match cli.command {
        Commands::Annotate {
            input,
            output,
            gff3,
            fasta,
            output_format,
            everything: _,
            buffer_size: _,
            pick,
            symbol: _,
            hgvs,
            canonical: _,
            distance,
            cache_dir,
            transcript_cache,
            sa_dir,
            sa_only,
            acmg,
            acmg_config,
            proband,
            mother,
            father,
            gene_list,
            explicit_alleles,
            qc_rules,
        } => {
            pipeline::run_annotate(pipeline::AnnotateConfig {
                input,
                output,
                gff3,
                fasta,
                output_format,
                pick,
                hgvs,
                distance,
                cache_dir,
                transcript_cache,
                sa_dir,
                sa_only,
                acmg,
                acmg_config,
                proband,
                mother,
                father,
                gene_list,
                explicit_alleles,
                qc_rules,
            })?;
        }
        Commands::Cache { gff3, fasta, output } => {
            pipeline::run_cache_build(&gff3, fasta.as_deref(), &output)?;
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
        } => {
            pipeline::run_sa_build(
                &source,
                &input,
                &output,
                &assembly,
                name.as_deref(),
                &info_fields,
            )?;
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
