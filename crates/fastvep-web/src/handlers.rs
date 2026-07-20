use axum::extract::State;
use axum::http::header;
use axum::response::IntoResponse;
use axum::Json;
use serde::{Deserialize, Serialize};
use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Arc, RwLock};
use std::time::Instant;

use crate::context::AnnotationContext;
use crate::errors::AppError;

#[derive(Serialize, Deserialize, Default)]
struct SavedStats {
    total_variants: u64,
    total_genomes: u64,
}

/// Shared application state.
pub struct SharedState {
    pub ctx: RwLock<AnnotationContext>,
    pub data_dir: Option<PathBuf>,
    pub sa_dir: Option<PathBuf>,
    pub stats_file: Option<PathBuf>,
    pub total_variants: AtomicU64,
    pub total_genomes: AtomicU64,
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
pub async fn status(State(state): State<AppState>) -> Json<serde_json::Value> {
    let (transcripts, gff3_source, has_fasta, sa_sources) = {
        let guard = state.ctx.read().unwrap();
        (
            guard.transcript_count(),
            guard.gff3_source.clone(),
            guard.seq_provider.is_some(),
            guard.sa_source_names(),
        )
    };
    Json(serde_json::json!({
        "status": "ok",
        "backend": true,
        "transcripts": transcripts,
        "gff3_source": gff3_source,
        "has_fasta": has_fasta,
        "sa_sources": sa_sources,
        "total_variants": state.total_variants.load(Ordering::Relaxed),
        "total_genomes": state.total_genomes.load(Ordering::Relaxed),
    }))
}

/// List available genome GFF3 files from the data directory.
/// Scans for .gff3, .gff3.gz, and .fastvep.cache files.
pub async fn list_genomes(State(state): State<AppState>) -> Json<serde_json::Value> {
    let Some(ref data_dir) = state.data_dir else {
        return Json(serde_json::json!({ "genomes": [] }));
    };

    let mut genomes = Vec::new();

    // Scan subdirectories: each subdir is a genome with GFF3 + optional FASTA
    if let Ok(entries) = std::fs::read_dir(data_dir) {
        for entry in entries.flatten() {
            let path = entry.path();
            if path.is_dir() {
                let name = path
                    .file_name()
                    .map(|n| n.to_string_lossy().to_string())
                    .unwrap_or_default();
                let gff3 = find_file_with_ext(&path, &["gff3", "gff3.gz", "fastvep.cache"]);
                let fasta = find_file_with_ext(&path, &["fa", "fasta", "fa.gz", "fasta.gz"]);
                let has_sa = path.join("sa").is_dir()
                    && std::fs::read_dir(path.join("sa"))
                        .map(|rd| {
                            rd.flatten().any(|e| {
                                let n = e.file_name().to_string_lossy().to_string();
                                n.ends_with(".osa") || n.ends_with(".osa2")
                            })
                        })
                        .unwrap_or(false);
                if gff3.is_some() {
                    genomes.push(serde_json::json!({
                        "name": name,
                        "has_fasta": fasta.is_some(),
                        "has_sa": has_sa,
                    }));
                }
            }
        }
    }

    // Also scan top-level for loose GFF3 files
    if let Ok(entries) = std::fs::read_dir(data_dir) {
        for entry in entries.flatten() {
            let path = entry.path();
            if path.is_file() {
                let fname = path
                    .file_name()
                    .map(|n| n.to_string_lossy().to_string())
                    .unwrap_or_default();
                if fname.ends_with(".gff3")
                    || fname.ends_with(".gff3.gz")
                    || fname.ends_with(".fastvep.cache")
                {
                    let stem = fname
                        .trim_end_matches(".fastvep.cache")
                        .trim_end_matches(".gz")
                        .trim_end_matches(".gff3")
                        .to_string();
                    genomes.push(serde_json::json!({
                        "name": stem,
                        "has_fasta": false,
                    }));
                }
            }
        }
    }

    genomes.sort_by(|a, b| {
        a["name"]
            .as_str()
            .unwrap_or("")
            .cmp(b["name"].as_str().unwrap_or(""))
    });

    Json(serde_json::json!({ "genomes": genomes }))
}

/// Load a genome from the data directory by name.
#[derive(Deserialize)]
pub struct LoadGenomeRequest {
    name: String,
}

pub async fn load_genome(
    State(state): State<AppState>,
    Json(req): Json<LoadGenomeRequest>,
) -> Result<Json<serde_json::Value>, AppError> {
    let Some(ref data_dir) = state.data_dir else {
        return Err(AppError::BadRequest("No data directory configured".into()));
    };

    let (gff3_path, fasta_path) = resolve_genome_paths(data_dir, &req.name)?;

    // Check for per-genome SA directory (data_dir/<name>/sa/)
    let genome_sa_dir = data_dir.join(&req.name).join("sa");
    let sa_dir = if genome_sa_dir.is_dir() {
        Some(genome_sa_dir)
    } else {
        // Fall back to global --sa-dir
        state.sa_dir.clone()
    };

    let name = req.name.clone();
    let ctx = Arc::clone(&state);

    let start = Instant::now();
    let (transcripts, sa_sources) = tokio::task::spawn_blocking(move || {
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
        state.total_genomes.fetch_add(1, Ordering::Relaxed);
        state.save_stats();
        Ok::<_, anyhow::Error>((tr_count, sa_names))
    })
    .await??;

    let time_ms = start.elapsed().as_millis() as u64;
    Ok(Json(serde_json::json!({
        "name": name,
        "transcripts": transcripts,
        "sa_sources": sa_sources,
        "time_ms": time_ms,
    })))
}

#[derive(Deserialize)]
pub struct AnnotateRequest {
    vcf: Option<String>,
    pick: Option<bool>,
    /// Enable ACMG-AMP variant classification for this request.
    acmg: Option<bool>,
}

pub async fn annotate(
    State(state): State<AppState>,
    Json(req): Json<AnnotateRequest>,
) -> Result<Json<serde_json::Value>, AppError> {
    let vcf_text = req.vcf.unwrap_or_default();
    if vcf_text.is_empty() {
        return Err(AppError::BadRequest("No VCF data provided".into()));
    }

    let pick = req.pick.unwrap_or(false);
    let acmg_requested = req.acmg.unwrap_or(false);
    let ctx = Arc::clone(&state);

    let start = Instant::now();
    let results = tokio::task::spawn_blocking(move || {
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
        // Borrow the existing config instead of cloning it; only build a
        // fresh default when none is loaded (no need to deep-clone
        // gene_overrides/ba1_exceptions on every ACMG-enabled request).
        let default_acmg;
        let acmg_config = if acmg_requested {
            Some(match guard.acmg_config.as_ref() {
                Some(cfg) => cfg,
                None => {
                    default_acmg = fastvep_classification::AcmgConfig::default();
                    &default_acmg
                }
            })
        } else {
            None
        };
        let results = guard.annotate_vcf_text_with_acmg(&vcf_text, pick, acmg_config)?;
        drop(guard);

        ctx.total_variants
            .fetch_add(results.len() as u64, Ordering::Relaxed);
        // Best-effort persistence, already on the blocking pool here — no
        // need for a second spawn_blocking just for this fs::write.
        ctx.save_stats();

        Ok::<_, anyhow::Error>(results)
    })
    .await??;

    let time_ms = start.elapsed().as_millis() as u64;
    Ok(Json(serde_json::json!({
        "results": results,
        "count": results.len(),
        "time_ms": time_ms,
    })))
}

pub async fn upload_gff3(
    State(state): State<AppState>,
    body: String,
) -> Result<Json<serde_json::Value>, AppError> {
    if body.is_empty() {
        return Err(AppError::BadRequest("No GFF3 data provided".into()));
    }

    let ctx = Arc::clone(&state);

    let start = Instant::now();
    let (genes, transcripts) = tokio::task::spawn_blocking(move || {
        let mut guard = ctx
            .ctx
            .write()
            .map_err(|e| anyhow::anyhow!("Lock poisoned: {}", e))?;
        guard.update_gff3_text(&body)
    })
    .await??;

    let time_ms = start.elapsed().as_millis() as u64;
    Ok(Json(serde_json::json!({
        "genes": genes,
        "transcripts": transcripts,
        "time_ms": time_ms,
    })))
}

// --- helpers ---

fn find_file_with_ext(dir: &std::path::Path, extensions: &[&str]) -> Option<PathBuf> {
    // Search in order of extension priority (first ext = highest priority)
    for ext in extensions {
        let entries = std::fs::read_dir(dir).ok()?;
        for entry in entries.flatten() {
            let path = entry.path();
            if path.is_file() {
                let fname = path.file_name()?.to_string_lossy();
                if fname.ends_with(ext) && !fname.ends_with(&format!(".fastvep.cache.{}", ext)) {
                    return Some(path);
                }
            }
        }
    }
    None
}

fn resolve_genome_paths(
    data_dir: &std::path::Path,
    name: &str,
) -> Result<(PathBuf, Option<PathBuf>), AppError> {
    // Sanitize: reject any name containing path separators or traversal sequences
    if name.contains('/') || name.contains('\\') || name.contains("..") {
        return Err(AppError::BadRequest(format!(
            "Genome '{}' not found in data directory",
            name
        )));
    }

    // Check if it's a subdirectory
    let subdir = data_dir.join(name);
    if subdir.is_dir() {
        let gff3 = find_file_with_ext(&subdir, &["gff3", "gff3.gz", "fastvep.cache"])
            .ok_or_else(|| AppError::BadRequest(format!("No GFF3 found in genome '{}'", name)))?;
        let fasta = find_file_with_ext(&subdir, &["fa", "fasta", "fa.gz", "fasta.gz"]);
        return Ok((gff3, fasta));
    }

    // Check top-level files
    for ext in &["gff3", "gff3.gz", "fastvep.cache"] {
        let path = data_dir.join(format!("{}.{}", name, ext));
        if path.exists() {
            return Ok((path, None));
        }
    }

    Err(AppError::BadRequest(format!(
        "Genome '{}' not found in data directory",
        name
    )))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;

    fn tmp_dir(label: &str) -> PathBuf {
        let dir = std::env::temp_dir().join(format!(
            "fastvep-web-test-{}-{}-{:?}",
            label,
            std::process::id(),
            std::thread::current().id()
        ));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir).unwrap();
        dir
    }

    #[test]
    fn resolve_genome_paths_rejects_path_traversal_sequences() {
        let data_dir = tmp_dir("traversal");
        for name in ["../etc/passwd", "..\\windows", "a/../../b", "a/b", "a\\b"] {
            let err = resolve_genome_paths(&data_dir, name).unwrap_err();
            let resp = err.into_response();
            assert_eq!(resp.status(), axum::http::StatusCode::BAD_REQUEST);
        }
    }

    #[test]
    fn resolve_genome_paths_finds_gff3_in_subdirectory() {
        let data_dir = tmp_dir("subdir");
        let genome_dir = data_dir.join("hg38");
        fs::create_dir_all(&genome_dir).unwrap();
        fs::write(genome_dir.join("annotation.gff3"), "##gff-version 3\n").unwrap();
        fs::write(genome_dir.join("genome.fa"), ">chr1\nACGT\n").unwrap();

        let (gff3, fasta) = resolve_genome_paths(&data_dir, "hg38").unwrap();
        assert_eq!(gff3, genome_dir.join("annotation.gff3"));
        assert_eq!(fasta, Some(genome_dir.join("genome.fa")));
    }

    #[test]
    fn resolve_genome_paths_finds_top_level_loose_gff3() {
        let data_dir = tmp_dir("toplevel");
        fs::write(data_dir.join("mygenome.gff3"), "##gff-version 3\n").unwrap();

        let (gff3, fasta) = resolve_genome_paths(&data_dir, "mygenome").unwrap();
        assert_eq!(gff3, data_dir.join("mygenome.gff3"));
        assert_eq!(fasta, None);
    }

    #[test]
    fn resolve_genome_paths_errors_when_genome_not_found() {
        let data_dir = tmp_dir("missing");
        let err = resolve_genome_paths(&data_dir, "nonexistent").unwrap_err();
        let resp = err.into_response();
        assert_eq!(resp.status(), axum::http::StatusCode::BAD_REQUEST);
    }

    #[test]
    fn find_file_with_ext_prefers_higher_priority_extension() {
        let dir = tmp_dir("priority");
        fs::write(dir.join("genome.fa.gz"), b"").unwrap();
        fs::write(dir.join("genome.fa"), b"").unwrap();

        let found = find_file_with_ext(&dir, &["fa", "fa.gz"]).unwrap();
        assert_eq!(found, dir.join("genome.fa"));
    }

    #[test]
    fn find_file_with_ext_skips_fastvep_cache_shadow_files() {
        let dir = tmp_dir("cache-shadow");
        // A `.fastvep.cache.gff3` shadow file should not be mistaken for a
        // real `.gff3` when searching by the `gff3` extension.
        fs::write(dir.join("genome.fastvep.cache.gff3"), b"").unwrap();
        fs::write(dir.join("genome.gff3"), b"##gff-version 3\n").unwrap();

        let found = find_file_with_ext(&dir, &["gff3"]).unwrap();
        assert_eq!(found, dir.join("genome.gff3"));
    }

    #[test]
    fn find_file_with_ext_returns_none_when_absent() {
        let dir = tmp_dir("absent");
        assert!(find_file_with_ext(&dir, &["gff3", "gff3.gz"]).is_none());
    }
}
