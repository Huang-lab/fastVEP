use axum::extract::State;
use axum::http::header;
use axum::response::IntoResponse;
use axum::Json;
use serde::{Deserialize, Serialize};
use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Arc, PoisonError, RwLock};
use std::time::Instant;

use fastvep_annotate::pick::{PickFlags, PickRequest};

use crate::context::AnnotationContext;
use crate::errors::AppError;
use crate::genomes::{resolve_genome_paths, scan_data_dir, GenomeRegistry};

#[derive(Serialize, Deserialize, Default)]
struct SavedStats {
    total_variants: u64,
    total_genomes: u64,
}

/// What a server is started with, besides the gene model itself.
pub struct ServerConfig {
    pub data_dir: Option<PathBuf>,
    pub sa_dir: Option<PathBuf>,
    pub stats_file: Option<PathBuf>,
    /// Upstream/downstream distance of the genomes loaded on request.
    pub distance: u64,
    pub allow_model_replacement: bool,
}

/// The name a gene model uploaded through `/api/upload-gff3` answers to. It is
/// also what `gff3_source` reads for such a model.
pub const UPLOADED_GENOME: &str = "user-upload";

/// Shared application state.
pub struct SharedState {
    /// The gene model the server started with: what a request that names no
    /// genome is answered from, and what `/api/load-genome` and
    /// `/api/upload-gff3` replace.
    pub ctx: RwLock<AnnotationContext>,
    /// The name `ctx` answers to, so that a request naming that genome is
    /// answered from `ctx` and a request naming any other is answered from the
    /// registry.
    ///
    /// It is written under the `ctx` write lock, by whatever replaces the model,
    /// and read under the `ctx` read lock, by whatever answers from it. That is
    /// what keeps a request for `hg38` from being answered from a mouse model
    /// that a replacement left behind under the old name. Always take `ctx`
    /// first; nothing takes this lock alone and then `ctx`.
    pub default_genome: RwLock<String>,
    /// The other genomes of the data directory, loaded when asked for.
    pub genomes: GenomeRegistry,
    pub data_dir: Option<PathBuf>,
    pub sa_dir: Option<PathBuf>,
    pub stats_file: Option<PathBuf>,
    pub total_variants: AtomicU64,
    pub total_genomes: AtomicU64,
    /// Whether `/api/upload-gff3` and `/api/load-genome` may replace the
    /// active gene model. Default off (`--allow-model-replacement`).
    ///
    /// The replacement is process-wide and permanent, so on the deployment
    /// `docs/API.md` describes - a lab server several colleagues share, behind
    /// a firewall - one client trying out a gene model silently changes what
    /// every other client is being answered against. Both clients are
    /// authorised and neither is doing anything wrong, which is why a firewall
    /// does not help and why this is off unless asked for. With it off, the
    /// gene model is a property of how the server was started (#113).
    ///
    /// Choosing a genome per request (`"genome"` in `/api/annotate`) replaces
    /// nothing and needs no flag.
    pub allow_model_replacement: bool,
    /// Bumped on every successful gene-model replacement.
    ///
    /// Reported alongside every annotation so a client can tell that two
    /// answers came from the same model. `gff3_source` alone cannot: every
    /// upload reports `user-upload`, and a transcript count can coincide. The
    /// counter is written under the `ctx` write lock and read under the read
    /// lock, so the lock - not the atomic's ordering - is what makes the value
    /// agree with the annotation it is reported with.
    pub gene_model_generation: AtomicU64,
}

/// Identity of the gene model an annotation was computed against.
///
/// Returned with every `/api/annotate` response. Before this, a client had no
/// way to tell: the answer was well-formed, HTTP 200, and computed against
/// whatever model happened to be loaded when it arrived. `/api/status` showed
/// the swap, but only to a client that polled it between every request.
#[derive(Serialize)]
struct GeneModelId {
    /// The genome that answered, by the name a request uses to ask for it. A
    /// client that chose a genome compares this with what it asked for.
    genome: String,
    gff3_source: Option<String>,
    transcripts: usize,
    generation: u64,
}

impl SharedState {
    /// A server answering from `ctx`, which is the genome called
    /// `default_genome`. Its counters start at zero; see [`Self::with_stats`].
    pub fn new(
        ctx: AnnotationContext,
        default_genome: impl Into<String>,
        config: ServerConfig,
    ) -> Self {
        Self {
            ctx: RwLock::new(ctx),
            default_genome: RwLock::new(default_genome.into()),
            genomes: GenomeRegistry::new(config.data_dir.clone(), config.distance),
            data_dir: config.data_dir,
            sa_dir: config.sa_dir,
            stats_file: config.stats_file,
            total_variants: AtomicU64::new(0),
            total_genomes: AtomicU64::new(0),
            allow_model_replacement: config.allow_model_replacement,
            gene_model_generation: AtomicU64::new(0),
        }
    }

    /// Start the usage counters from the values a previous run saved.
    pub fn with_stats(self, total_variants: u64, total_genomes: u64) -> Self {
        self.total_variants.store(total_variants, Ordering::Relaxed);
        self.total_genomes.store(total_genomes, Ordering::Relaxed);
        self
    }

    /// The name the gene model in `ctx` answers to. Read it while holding a
    /// `ctx` guard when it has to agree with what that guard is used for.
    fn current_default_genome(&self) -> String {
        self.default_genome
            .read()
            .unwrap_or_else(PoisonError::into_inner)
            .clone()
    }

    /// Read the default model's identity. Call while holding a `ctx` guard so
    /// the identity cannot be from a different model than the caller used.
    fn gene_model_id(&self, ctx: &AnnotationContext) -> GeneModelId {
        GeneModelId {
            genome: self.current_default_genome(),
            gff3_source: ctx.gff3_source.clone(),
            transcripts: ctx.transcript_count(),
            generation: self.gene_model_generation.load(Ordering::Relaxed),
        }
    }

    /// Refuse a gene-model replacement unless the server was started to allow
    /// one. The message names the flag rather than just saying no, because the
    /// browser UI's own workflow needs these endpoints and a single-user
    /// desktop run is exactly the case that should turn them on.
    fn check_model_replacement_allowed(&self) -> Result<(), AppError> {
        if self.allow_model_replacement {
            return Ok(());
        }
        Err(AppError::Forbidden(
            "Replacing the active gene model is disabled on this server. A replacement is \
             process-wide, so it changes what every other client is answered against. Start \
             fastvep-web with --allow-model-replacement to enable it."
                .into(),
        ))
    }

    /// Count one annotation run.
    ///
    /// `named` is whether the request chose its genome. Those runs are what the
    /// page's "Genome Analysis Sessions" counter counts: it used to count
    /// `/api/load-genome` calls, which the page made once per genome choice and
    /// which a server without `--allow-model-replacement` refuses.
    fn record_run(&self, variants: usize, named: bool) {
        self.total_variants
            .fetch_add(variants as u64, Ordering::Relaxed);
        if named {
            self.total_genomes.fetch_add(1, Ordering::Relaxed);
        }
        // Best-effort persistence, already on the blocking pool when this is
        // called - no need for a second spawn_blocking just for this fs::write.
        self.save_stats();
    }

    /// Load every genome of the data directory, one after another, so that no
    /// request has to wait for one. The genome the server started with is
    /// already loaded and is skipped. A genome that fails to load is reported
    /// and left to fail again, loudly, on its first request. Returns how many
    /// loaded.
    pub async fn preload_genomes(self: Arc<Self>) -> usize {
        let Some(dir) = self.data_dir.clone() else {
            return 0;
        };
        let Ok(listed) = tokio::task::spawn_blocking(move || scan_data_dir(&dir)).await else {
            return 0;
        };
        let default = self.current_default_genome();
        let mut loaded = 0;
        for listing in listed.into_iter().filter(|g| g.name != default) {
            match self.genomes.get_or_load(&listing.name).await {
                Ok(_) => loaded += 1,
                Err(e) => tracing::warn!("Could not preload genome {}: {e:?}", listing.name),
            }
        }
        loaded
    }
}

impl SharedState {
    pub fn save_stats(&self) {
        if let Some(ref path) = self.stats_file {
            let stats = SavedStats {
                total_variants: self.total_variants.load(Ordering::Relaxed),
                total_genomes: self.total_genomes.load(Ordering::Relaxed),
            };
            if let Ok(json) = serde_json::to_string(&stats) {
                let _ = std::fs::write(path, json);
            }
        }
    }
}

pub type AppState = Arc<SharedState>;

const INDEX_HTML: &str = include_str!("../../../web/index.html");
const LOGO_PNG: &[u8] = include_bytes!("../../../web/assets/logo.png");

pub async fn index_html() -> impl IntoResponse {
    (
        [(header::CONTENT_TYPE, "text/html; charset=utf-8")],
        INDEX_HTML,
    )
}

/// Serve the fastVEP logo PNG. Used as both the page logo and the browser
/// tab favicon (via <link rel="icon"> and <link rel="apple-touch-icon">).
pub async fn logo_png() -> impl IntoResponse {
    (
        [
            (header::CONTENT_TYPE, "image/png"),
            (header::CACHE_CONTROL, "public, max-age=86400"),
        ],
        LOGO_PNG,
    )
}

/// Enhanced status: reports transcript/genome/SA state so the SPA can
/// decide whether to use server-side annotation or upload example GFF3.
///
/// The model fields describe the default genome, as they always have.
/// `default_genome` is also how a client learns that this server answers per
/// request: an older server does not send it.
pub async fn status(State(state): State<AppState>) -> Json<serde_json::Value> {
    let (transcripts, gff3_source, has_fasta, sa_sources, generation, default_genome) = {
        let guard = state.ctx.read().unwrap();
        (
            guard.transcript_count(),
            guard.gff3_source.clone(),
            guard.seq_provider.is_some(),
            guard.sa_source_names(),
            state.gene_model_generation.load(Ordering::Relaxed),
            state.current_default_genome(),
        )
    };
    Json(serde_json::json!({
        "status": "ok",
        "backend": true,
        "version": env!("CARGO_PKG_VERSION"),
        "transcripts": transcripts,
        "gff3_source": gff3_source,
        "gene_model_generation": generation,
        "default_genome": default_genome,
        "has_fasta": has_fasta,
        "sa_sources": sa_sources,
        // So a UI can offer the model-replacing controls only when they will
        // work, instead of finding out with a 403 the user cannot interpret.
        "allow_model_replacement": state.allow_model_replacement,
        "total_variants": state.total_variants.load(Ordering::Relaxed),
        "total_genomes": state.total_genomes.load(Ordering::Relaxed),
    }))
}

/// The genomes a request can name: everything the data directory holds, plus
/// the one the server started with when that is not among them.
///
/// `loaded` says whether a request for the genome is answered at once or waits
/// for it to load, and `transcripts` is what it holds once it has.
pub async fn list_genomes(State(state): State<AppState>) -> Json<serde_json::Value> {
    let scanned = match state.data_dir.clone() {
        Some(dir) => tokio::task::spawn_blocking(move || scan_data_dir(&dir))
            .await
            .unwrap_or_default(),
        None => Vec::new(),
    };

    // Taken under one guard, so the entry and the name it carries agree.
    let (default_name, default_gff3, default_fasta, default_sa, default_transcripts) = {
        let guard = state.ctx.read().unwrap();
        (
            state.current_default_genome(),
            guard.gff3_source.clone().unwrap_or_default(),
            guard.seq_provider.is_some(),
            !guard.sa_source_names().is_empty(),
            guard.transcript_count(),
        )
    };

    let mut genomes: Vec<serde_json::Value> = Vec::with_capacity(scanned.len() + 1);
    if !scanned.iter().any(|g| g.name == default_name) {
        genomes.push(serde_json::json!({
            "name": default_name,
            "gff3": default_gff3,
            "has_fasta": default_fasta,
            "has_sa": default_sa,
            "loaded": true,
            "default": true,
            "transcripts": default_transcripts,
        }));
    }
    for listing in scanned {
        let is_default = listing.name == default_name;
        let transcripts = if is_default {
            Some(default_transcripts)
        } else {
            state
                .genomes
                .loaded(&listing.name)
                .map(|g| g.transcript_count())
        };
        genomes.push(serde_json::json!({
            "name": listing.name,
            "gff3": listing.gff3,
            "has_fasta": listing.has_fasta,
            "has_sa": listing.has_sa,
            "loaded": transcripts.is_some(),
            "default": is_default,
            "transcripts": transcripts,
        }));
    }

    Json(serde_json::json!({ "default": default_name, "genomes": genomes }))
}

/// Load a genome from the data directory by name.
#[derive(Deserialize)]
pub struct LoadGenomeRequest {
    name: String,
}

/// Replace the default gene model with another genome of the data directory.
///
/// Choosing a genome per request (`"genome"` in `/api/annotate`) does not
/// need this, and does not move anyone else's answers; this remains for the
/// single-user run that wants the default itself to change.
pub async fn load_genome(
    State(state): State<AppState>,
    Json(req): Json<LoadGenomeRequest>,
) -> Result<Json<serde_json::Value>, AppError> {
    state.check_model_replacement_allowed()?;
    let Some(ref data_dir) = state.data_dir else {
        return Err(AppError::BadRequest("No data directory configured".into()));
    };

    let paths = resolve_genome_paths(data_dir, &req.name)?;
    let gff3_path = paths.gff3;
    let fasta_path = paths.fasta;

    // Per-genome SA directory, else the global `--sa-dir`. The path comes back
    // from `resolve_genome_paths` rather than being re-joined here: rebuilding
    // it from `req.name` was safe only because the validation happened to run
    // first, and the traversal test did not cover the second join.
    //
    // Unlike a genome chosen per request, which never borrows another genome's
    // SA directory, this keeps the fallback: a single-user run with one
    // `--sa-dir` and several gene models is what it was written for.
    let sa_dir = paths.sa_dir.or_else(|| state.sa_dir.clone());

    let name = req.name.clone();
    let ctx = Arc::clone(&state);

    let start = Instant::now();
    let (transcripts, sa_sources, generation) = tokio::task::spawn_blocking(move || {
        let mut guard = ctx
            .ctx
            .write()
            .map_err(|e| anyhow::anyhow!("Lock poisoned: {}", e))?;
        let gff3_str = gff3_path
            .to_str()
            .ok_or_else(|| anyhow::anyhow!("GFF3 path is not valid UTF-8: {:?}", gff3_path))?;
        let fasta_str = fasta_path
            .as_ref()
            .map(|p| {
                p.to_str()
                    .ok_or_else(|| anyhow::anyhow!("FASTA path is not valid UTF-8: {:?}", p))
            })
            .transpose()?;
        let sa_str = sa_dir
            .as_ref()
            .map(|p| {
                p.to_str()
                    .ok_or_else(|| anyhow::anyhow!("SA dir path is not valid UTF-8: {:?}", p))
            })
            .transpose()?;
        guard.load_genome(gff3_str, fasta_str, sa_str)?;
        let tr_count = guard.transcript_count();
        let sa_names = guard.sa_source_names();
        // Under the write lock, as in `upload_gff3`: the model and the name it
        // answers to change together, so no reader sees one without the other.
        *ctx.default_genome
            .write()
            .map_err(|e| anyhow::anyhow!("Lock poisoned: {}", e))? = name;
        let generation = ctx.gene_model_generation.fetch_add(1, Ordering::Relaxed) + 1;
        ctx.total_genomes.fetch_add(1, Ordering::Relaxed);
        ctx.save_stats();
        Ok::<_, anyhow::Error>((tr_count, sa_names, generation))
    })
    .await??;

    let time_ms = start.elapsed().as_millis() as u64;
    Ok(Json(serde_json::json!({
        "name": req.name,
        "transcripts": transcripts,
        "sa_sources": sa_sources,
        "generation": generation,
        "time_ms": time_ms,
    })))
}

#[derive(Deserialize)]
pub struct AnnotateRequest {
    vcf: Option<String>,
    /// The genome to answer from, by the name `GET /api/genomes` lists.
    ///
    /// Left out, the request is answered from the genome the server started
    /// with. A name that is not a genome is refused; it is never answered from
    /// the default instead.
    genome: Option<String>,
    /// `pick`, `pick_allele`, `pick_allele_gene`, `flag_pick`,
    /// `flag_pick_allele`, `flag_pick_allele_gene` - flattened so the six
    /// names come from `PickFlags` itself and cannot drift from the CLI's or
    /// from the bundled server's.
    #[serde(flatten)]
    pick: PickFlags,
    /// Enable ACMG-AMP variant classification for this request.
    acmg: Option<bool>,
}

/// What the default model made of a request.
enum DefaultAnswer {
    Answered(Vec<serde_json::Value>, GeneModelId),
    /// The request named a genome and the default model is not it. The VCF text
    /// comes back so that a request of up to the body limit is not copied.
    NotTheDefault(String),
}

/// Annotate with one context, attaching ACMG when the request asked for it.
fn annotate_with(
    ctx: &AnnotationContext,
    vcf_text: &str,
    pick: Option<PickRequest>,
    acmg_requested: bool,
) -> anyhow::Result<Vec<serde_json::Value>> {
    // Borrow the existing config instead of cloning it; only build a
    // fresh default when none is loaded (no need to deep-clone
    // gene_overrides/ba1_exceptions on every ACMG-enabled request).
    let default_acmg;
    let acmg_config = if acmg_requested {
        Some(match ctx.acmg_config.as_ref() {
            Some(cfg) => cfg,
            None => {
                default_acmg = fastvep_classification::AcmgConfig::default();
                &default_acmg
            }
        })
    } else {
        None
    };
    ctx.annotate_vcf_text_with_acmg(vcf_text, pick, acmg_config)
}

/// Answer from a genome of the data directory other than the default one,
/// loading it first if this is the first request for it.
async fn answer_from_registry(
    state: &AppState,
    name: String,
    vcf_text: String,
    pick: Option<PickRequest>,
    acmg_requested: bool,
) -> Result<(Vec<serde_json::Value>, GeneModelId), AppError> {
    if state.genomes.data_dir().is_none() {
        return Err(AppError::BadRequest(format!(
            "Genome '{}' not found. This server was started without --data-dir and serves only '{}'.",
            name,
            state.current_default_genome()
        )));
    }
    let genome = state.genomes.get_or_load(&name).await?;
    let ctx = Arc::clone(state);
    let answered = tokio::task::spawn_blocking(move || {
        let results = annotate_with(&genome, &vcf_text, pick, acmg_requested)?;
        let gene_model = GeneModelId {
            genome: name,
            gff3_source: genome.gff3_source.clone(),
            transcripts: genome.transcript_count(),
            // A genome loaded on request is never replaced, so there is
            // nothing for a generation to count.
            generation: 0,
        };
        ctx.record_run(results.len(), true);
        Ok::<_, anyhow::Error>((results, gene_model))
    })
    .await??;
    Ok(answered)
}

pub async fn annotate(
    State(state): State<AppState>,
    Json(req): Json<AnnotateRequest>,
) -> Result<Json<serde_json::Value>, AppError> {
    let vcf_text = req.vcf.unwrap_or_default();
    if vcf_text.is_empty() {
        return Err(AppError::BadRequest("No VCF data provided".into()));
    }
    // An empty name is a client that meant to choose and did not, which is not
    // the same as one that left the field out.
    if req.genome.as_deref() == Some("") {
        return Err(AppError::BadRequest(
            "genome must not be empty: name one of GET /api/genomes, or leave it out for the \
             default"
                .into(),
        ));
    }

    let pick = req.pick.resolve();
    let acmg_requested = req.acmg.unwrap_or(false);
    let genome = req.genome;
    let wanted = genome.clone();
    let ctx = Arc::clone(&state);

    let start = Instant::now();
    let first = tokio::task::spawn_blocking(move || {
        // A read lock is sufficient: annotate_vcf_text_with_acmg only needs
        // &self, and the ACMG toggle is now passed as a per-call argument
        // instead of mutating the shared context. Previously this took a
        // write lock just to flip guard.acmg_config, which serialized every
        // concurrent request (including unrelated /api/status reads) behind
        // whichever annotation was running, and let one request's ACMG
        // preference clobber another's mid-flight.
        let guard = ctx
            .ctx
            .read()
            .map_err(|e| anyhow::anyhow!("Lock poisoned: {}", e))?;
        // Whether the default model is the genome asked for is decided under
        // the same guard that then answers from it. Decided earlier, a
        // replacement could land in between and the request would be answered
        // from a model that no longer has that name.
        if let Some(ref name) = wanted {
            if ctx.current_default_genome() != *name {
                return Ok::<_, anyhow::Error>(DefaultAnswer::NotTheDefault(vcf_text));
            }
        }
        let results = annotate_with(&guard, &vcf_text, pick, acmg_requested)?;
        // Read the model's identity before releasing the lock, so what the
        // response reports is the model that produced these results and not
        // whichever one a concurrent upload left behind (#113).
        let gene_model = ctx.gene_model_id(&guard);
        drop(guard);

        ctx.record_run(results.len(), wanted.is_some());

        Ok(DefaultAnswer::Answered(results, gene_model))
    })
    .await??;

    let (results, gene_model) = match first {
        DefaultAnswer::Answered(results, gene_model) => (results, gene_model),
        DefaultAnswer::NotTheDefault(vcf_text) => {
            let name = genome.ok_or_else(|| {
                AppError::Internal(anyhow::anyhow!(
                    "a request that named no genome was handed back by the default model"
                ))
            })?;
            answer_from_registry(&state, name, vcf_text, pick, acmg_requested).await?
        }
    };

    let time_ms = start.elapsed().as_millis() as u64;
    Ok(Json(serde_json::json!({
        "results": results,
        "count": results.len(),
        "time_ms": time_ms,
        "gene_model": gene_model,
    })))
}

pub async fn upload_gff3(
    State(state): State<AppState>,
    body: String,
) -> Result<Json<serde_json::Value>, AppError> {
    state.check_model_replacement_allowed()?;
    if body.is_empty() {
        return Err(AppError::BadRequest("No GFF3 data provided".into()));
    }

    let ctx = Arc::clone(&state);

    let start = Instant::now();
    let (genes, transcripts, generation) = tokio::task::spawn_blocking(move || {
        let mut guard = ctx
            .ctx
            .write()
            .map_err(|e| anyhow::anyhow!("Lock poisoned: {}", e))?;
        let (genes, transcripts) = guard.update_gff3_text(&body)?;
        // Bumped under the write lock, so a reader holding the read lock sees
        // a generation that matches the model it is about to use. The name
        // changes with it: the upload is no longer the genome it replaced.
        *ctx.default_genome
            .write()
            .map_err(|e| anyhow::anyhow!("Lock poisoned: {}", e))? = UPLOADED_GENOME.to_string();
        let generation = ctx.gene_model_generation.fetch_add(1, Ordering::Relaxed) + 1;
        Ok::<_, anyhow::Error>((genes, transcripts, generation))
    })
    .await??;

    let time_ms = start.elapsed().as_millis() as u64;
    Ok(Json(serde_json::json!({
        "genes": genes,
        "transcripts": transcripts,
        "generation": generation,
        "time_ms": time_ms,
    })))
}
