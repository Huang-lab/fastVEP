use clap::builder::BoolishValueParser;
use clap::Parser;
use std::path::PathBuf;
use std::sync::Arc;

use fastvep_web::build_router;
use fastvep_web::context::AnnotationContext;
use fastvep_web::genomes::default_genome_name;
use fastvep_web::handlers::{AppState, ServerConfig, SharedState};

#[derive(Parser)]
#[command(name = "fastvep-web")]
#[command(about = "fastVEP production web server for variant annotation")]
struct Cli {
    /// Port to listen on
    #[arg(long, default_value_t = 8080, env = "FASTVEP_PORT")]
    port: u16,

    /// Bind address
    #[arg(long, default_value = "0.0.0.0", env = "FASTVEP_BIND")]
    bind: String,

    /// GFF3 annotation file to load at startup
    #[arg(long, env = "FASTVEP_GFF3")]
    gff3: Option<String>,

    /// FASTA reference file (requires .fai index for mmap)
    #[arg(long, env = "FASTVEP_FASTA")]
    fasta: Option<String>,

    /// Supplementary annotation directory (.osa/.osa2 files)
    #[arg(long, env = "FASTVEP_SA_DIR")]
    sa_dir: Option<String>,

    /// Directory containing genome data (subdirs with GFF3 + FASTA per organism)
    ///
    /// Every genome in it can be named in the `genome` field of
    /// `POST /api/annotate`. Each is loaded the first time it is asked for and
    /// then kept, from its own GFF3, FASTA and `sa/` directory; nothing is
    /// shared with the gene model given by `--gff3`, and no request changes the
    /// answers of another.
    #[arg(long, env = "FASTVEP_DATA_DIR")]
    data_dir: Option<String>,

    /// Load every genome of `--data-dir` in the background at startup, instead
    /// of on its first request.
    ///
    /// Off by default: a data directory can hold more genomes than a small
    /// machine wants resident. With it on the server answers from the `--gff3`
    /// genome at once and the rest become available as the background load
    /// reaches them, so no visitor waits for a genome to load. A genome that
    /// fails to load is logged and fails again, with the reason, on its first
    /// request.
    #[arg(
        long,
        action = clap::ArgAction::SetTrue,
        value_parser = BoolishValueParser::new(),
        env = "FASTVEP_PRELOAD_GENOMES"
    )]
    preload_genomes: bool,

    /// Upstream/downstream distance in bp
    #[arg(long, default_value_t = 5000)]
    distance: u64,

    /// Maximum request body size in bytes
    #[arg(long, default_value_t = 10_485_760)]
    max_body_size: usize,

    /// Maximum concurrent annotation requests
    #[arg(long, default_value_t = 64)]
    max_concurrent: usize,

    /// Path to a file for persistent statistics tracking
    #[arg(long, default_value = "stats.json", env = "FASTVEP_STATS_FILE")]
    stats_file: PathBuf,

    /// Allow `/api/upload-gff3` and `/api/load-genome` to replace the active
    /// gene model.
    ///
    /// Off by default. A replacement is process-wide and permanent, so on a
    /// shared server one client trying out a gene model changes what every
    /// other client is answered against - with no error, and nothing in an
    /// `/api/annotate` response to say the model moved. Turn it on for a
    /// single-user desktop run, where the browser UI's own genome switching
    /// needs it; leave it off for anything several people share, and the gene
    /// model becomes a property of how the server was started.
    #[arg(
        long,
        action = clap::ArgAction::SetTrue,
        value_parser = BoolishValueParser::new(),
        env = "FASTVEP_ALLOW_MODEL_REPLACEMENT"
    )]
    allow_model_replacement: bool,

    /// Answer `"acmg": true` for a genome whose database directory lacks part of the ACMG
    /// reference stack.
    ///
    /// Off by default. The classifier reads nine supplementary sources
    /// (docs/ACMG_SETUP.md); without one, the criteria that read it report
    /// `evaluated: false` and the result is a classification that rests on less evidence than
    /// it appears to, almost always VUS. A server that holds ClinVar and dbSNP only would
    /// answer every request for ACMG that way, and nothing in a VUS says it was not given the
    /// evidence. So by default such a request is refused with HTTP 422 and the list of what is
    /// missing, and `GET /api/status` and `GET /api/genomes` say which genomes can be
    /// classified. Turn this on for a run where an incomplete stack is deliberate: every
    /// response then carries `acmg_stack` with the sources that were missing.
    #[arg(
        long,
        action = clap::ArgAction::SetTrue,
        value_parser = BoolishValueParser::new(),
        env = "FASTVEP_ACMG_ALLOW_PARTIAL"
    )]
    acmg_allow_partial: bool,
}

#[tokio::main]
async fn main() -> anyhow::Result<()> {
    tracing_subscriber::fmt()
        .with_env_filter(
            tracing_subscriber::EnvFilter::try_from_default_env()
                .unwrap_or_else(|_| "fastvep_web=info".into()),
        )
        .init();

    let cli = Cli::parse();

    tracing::info!("Loading annotation data...");
    let ctx = AnnotationContext::new(
        cli.gff3.as_deref(),
        cli.fasta.as_deref(),
        cli.sa_dir.as_deref(),
        cli.distance,
    )?;

    let data_dir = cli.data_dir.map(PathBuf::from);
    let sa_dir = cli.sa_dir.as_ref().map(PathBuf::from);
    if let Some(ref dir) = data_dir {
        tracing::info!("Data directory: {}", dir.display());
    }
    if let Some(ref dir) = sa_dir {
        tracing::info!("SA directory: {}", dir.display());
    }

    // Load persistent stats
    let mut initial_variants = 0;
    let mut initial_genomes = 0;
    if cli.stats_file.exists() {
        if let Ok(data) = std::fs::read_to_string(&cli.stats_file) {
            if let Ok(stats) = serde_json::from_str::<serde_json::Value>(&data) {
                initial_variants = stats["total_variants"].as_u64().unwrap_or(0);
                initial_genomes = stats["total_genomes"].as_u64().unwrap_or(0);
                tracing::info!(
                    "Loaded persistent stats: {} variants, {} genomes",
                    initial_variants,
                    initial_genomes
                );
            }
        }
    }

    if cli.acmg_allow_partial {
        tracing::warn!(
            "--acmg-allow-partial is set: ACMG classification is answered for genomes whose \
             database directory lacks part of the reference stack. Each response says which \
             sources were missing."
        );
    }
    if cli.allow_model_replacement {
        tracing::warn!(
            "--allow-model-replacement is set: any client may replace the gene model for every \
             other client. Intended for a single-user run."
        );
    }

    let default_genome = default_genome_name(cli.gff3.as_deref(), data_dir.as_deref());
    tracing::info!("Default genome: {default_genome}");

    let state: AppState = Arc::new(
        SharedState::new(
            ctx,
            default_genome,
            ServerConfig {
                data_dir,
                sa_dir,
                stats_file: Some(cli.stats_file),
                distance: cli.distance,
                allow_model_replacement: cli.allow_model_replacement,
                acmg_allow_partial: cli.acmg_allow_partial,
            },
        )
        .with_stats(initial_variants, initial_genomes),
    );

    if cli.preload_genomes {
        let preloading = Arc::clone(&state);
        tokio::spawn(async move {
            let loaded = preloading.preload_genomes().await;
            tracing::info!("Preloaded {loaded} more genomes");
        });
    }

    let app = build_router(state, cli.max_body_size, cli.max_concurrent);

    let addr = format!("{}:{}", cli.bind, cli.port);
    let listener = tokio::net::TcpListener::bind(&addr).await?;
    tracing::info!("fastVEP web server listening on http://{}", addr);

    axum::serve(listener, app)
        .with_graceful_shutdown(shutdown_signal())
        .await?;

    Ok(())
}

async fn shutdown_signal() {
    tokio::signal::ctrl_c()
        .await
        .expect("Failed to listen for ctrl-c");
    tracing::info!("Shutdown signal received, draining connections...");
}

#[cfg(test)]
mod tests {
    use super::*;

    // Both are plain switches. They read their environment variable through
    // `BoolishValueParser` so that `FASTVEP_PRELOAD_GENOMES=1`, which is what
    // anyone writes in a systemd unit, is accepted: clap's default for a bool
    // takes only `true` and `false`, and refuses `1` at startup with
    // "invalid value '1'".
    #[test]
    fn the_switches_are_flags_and_default_to_off() {
        let off = Cli::try_parse_from(["fastvep-web"]).unwrap();
        assert!(!off.preload_genomes && !off.allow_model_replacement && !off.acmg_allow_partial);

        let on = Cli::try_parse_from([
            "fastvep-web",
            "--preload-genomes",
            "--allow-model-replacement",
            "--acmg-allow-partial",
        ])
        .unwrap();
        assert!(on.preload_genomes && on.allow_model_replacement && on.acmg_allow_partial);
    }

    #[test]
    fn the_value_parser_accepts_the_usual_spellings_of_true_and_false() {
        use clap::builder::TypedValueParser;
        let parser = BoolishValueParser::new();
        let cmd = clap::Command::new("x");
        for yes in ["1", "true", "yes", "on", "y", "t"] {
            let v = parser.parse_ref(&cmd, None, std::ffi::OsStr::new(yes));
            assert!(v.unwrap(), "{yes}");
        }
        for no in ["0", "false", "no", "off", "n", "f"] {
            let v = parser.parse_ref(&cmd, None, std::ffi::OsStr::new(no));
            assert!(!v.unwrap(), "{no}");
        }
    }
}
