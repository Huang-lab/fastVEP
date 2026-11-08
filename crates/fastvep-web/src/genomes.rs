//! The genomes a server can answer from, chosen per request.
//!
//! `fastvep-web` starts with one gene model, the one `--gff3` names. The only way
//! to answer from another genome used to be to replace that model for the whole
//! process (`/api/load-genome`), which changes what every other client is
//! answered against, so it is off unless `--allow-model-replacement` is given. A
//! server several people share therefore could not serve more than one genome.
//!
//! A request now names the genome it wants (`"genome"` in `POST /api/annotate`)
//! and is answered from that genome's own model, which the registry loads on
//! first use and keeps. Nothing is replaced, so no client can change what another
//! is answered against, and the response says which genome answered.
//!
//! A genome is a subdirectory of `--data-dir` holding a GFF3 (or its sidecar
//! transcript cache), optionally a FASTA, and optionally an `sa/` directory of
//! supplementary annotations.
//!
//! Two rules keep an answer from being plausible and wrong:
//!
//! - An unknown name is an error. It is never answered from the default genome,
//!   because a misspelled name that quietly returns human consequences for a
//!   mouse variant looks exactly like a correct answer.
//! - A genome's supplementary annotations come from its own `sa/` directory and
//!   from nowhere else. The `--sa-dir` the default genome was started with is
//!   never lent to another one: on the hosted site it holds human ClinVar and
//!   dbSNP, and a mouse variant that shares a coordinate with a human record
//!   would be annotated with the human record.

use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};
use std::time::Instant;

use tokio::sync::watch;

use crate::context::AnnotationContext;
use crate::errors::AppError;

const GFF3_EXTENSIONS: [&str; 3] = ["gff3", "gff3.gz", "fastvep.cache"];
const FASTA_EXTENSIONS: [&str; 4] = ["fa", "fasta", "fa.gz", "fasta.gz"];

/// One genome a server can be asked for, as found on disk.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GenomeListing {
    pub name: String,
    /// File name, without its directory, of the GFF3 the genome is read from, or
    /// of its sidecar cache when that is all there is. Ensembl names these
    /// `<Species>.<Assembly>.<release>.gff3`, which is how a client can say what
    /// the genome is without the server describing it.
    pub gff3: String,
    pub has_fasta: bool,
    pub has_sa: bool,
}

/// Every genome under the data directory: each subdirectory that holds a GFF3,
/// and each loose GFF3 (or cache) file at the top level. Sorted by name, one
/// entry per name, a directory winning over a loose file of the same name just
/// as [`resolve_genome_paths`] would resolve it.
pub fn scan_data_dir(data_dir: &Path) -> Vec<GenomeListing> {
    let mut genomes: Vec<GenomeListing> = Vec::new();
    let Ok(entries) = std::fs::read_dir(data_dir) else {
        return genomes;
    };
    let mut loose: Vec<GenomeListing> = Vec::new();
    for entry in entries.flatten() {
        let path = entry.path();
        let file_name = entry.file_name().to_string_lossy().into_owned();
        if path.is_dir() {
            if let Some(gff3) = find_file_with_ext(&path, &GFF3_EXTENSIONS) {
                genomes.push(GenomeListing {
                    name: file_name,
                    gff3: base_name(&gff3),
                    has_fasta: find_file_with_ext(&path, &FASTA_EXTENSIONS).is_some(),
                    has_sa: has_sa_files(&path.join("sa")),
                });
            }
        } else if path.is_file()
            && GFF3_EXTENSIONS
                .iter()
                .any(|e| file_name.ends_with(&format!(".{e}")))
        {
            let name = file_name
                .trim_end_matches(".fastvep.cache")
                .trim_end_matches(".gz")
                .trim_end_matches(".gff3")
                .to_string();
            loose.push(GenomeListing {
                name,
                gff3: file_name,
                has_fasta: false,
                has_sa: false,
            });
        }
    }
    for listing in loose {
        if !genomes.iter().any(|g| g.name == listing.name) {
            genomes.push(listing);
        }
    }
    genomes.sort_by(|a, b| a.name.cmp(&b.name));
    genomes
}

fn base_name(path: &Path) -> String {
    path.file_name()
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_default()
}

fn has_sa_files(sa_dir: &Path) -> bool {
    std::fs::read_dir(sa_dir)
        .map(|entries| {
            entries.flatten().any(|e| {
                let n = e.file_name().to_string_lossy().into_owned();
                n.ends_with(".osa") || n.ends_with(".osa2")
            })
        })
        .unwrap_or(false)
}

/// The name the gene model a server starts with answers to.
///
/// A model that lives in a genome directory of the data directory is that
/// genome, and answers to its directory name. That is what lets `GET
/// /api/genomes`, which scans the directory, and a request that names the
/// genome agree about which entry is the one already in memory, instead of
/// loading it a second time. Any other model is named for its file, and a
/// server started with no gene model at all answers to `default`.
pub fn default_genome_name(gff3: Option<&str>, data_dir: Option<&Path>) -> String {
    let Some(gff3) = gff3 else {
        return "default".to_string();
    };
    let path = Path::new(gff3);
    if let (Some(data_dir), Some(genome_dir)) = (data_dir, path.parent()) {
        let holds_it = genome_dir
            .parent()
            .is_some_and(|parent| same_directory(parent, data_dir));
        if let (true, Some(name)) = (holds_it, genome_dir.file_name()) {
            return name.to_string_lossy().into_owned();
        }
    }
    let file = base_name(path);
    let stem = file
        .trim_end_matches(".gz")
        .trim_end_matches(".gff3")
        .to_string();
    if stem.is_empty() {
        "default".to_string()
    } else {
        stem
    }
}

fn same_directory(a: &Path, b: &Path) -> bool {
    match (std::fs::canonicalize(a), std::fs::canonicalize(b)) {
        (Ok(a), Ok(b)) => a == b,
        _ => a == b,
    }
}

/// Every path a named genome contributes, all derived from the one validated
/// join.
///
/// The SA directory used to be re-joined at the call site from the raw request
/// name. That was safe, but by ordering rather than by construction - the
/// traversal check lives here, and nothing stopped a second join from being
/// written before it. Returning the path makes the sanitisation structural, so
/// the traversal test covers it.
#[derive(Debug)]
pub struct GenomePaths {
    pub gff3: PathBuf,
    pub fasta: Option<PathBuf>,
    /// `data_dir/<name>/sa`, when it exists. `None` for a top-level loose
    /// GFF3, which has no directory to hold one.
    pub sa_dir: Option<PathBuf>,
}

pub fn resolve_genome_paths(data_dir: &Path, name: &str) -> Result<GenomePaths, AppError> {
    // Sanitize: reject any name containing path separators or traversal sequences
    if name.contains('/') || name.contains('\\') || name.contains("..") {
        return Err(AppError::BadRequest(format!(
            "Genome '{}' not found in data directory",
            name
        )));
    }

    // A name has to be one that is on disk, letter for letter. On a
    // case-insensitive filesystem (macOS and Windows, where a self-hosted server
    // is often run) `data_dir.join("Alpha")` finds the directory `alpha`, and the
    // request would load a second copy of it under a name that `GET /api/genomes`
    // never listed.
    if !has_entry_named(data_dir, name) {
        return Err(AppError::BadRequest(format!(
            "Genome '{}' not found in data directory",
            name
        )));
    }

    // Check if it's a subdirectory
    let subdir = data_dir.join(name);
    if subdir.is_dir() {
        let gff3 = find_file_with_ext(&subdir, &GFF3_EXTENSIONS)
            .ok_or_else(|| AppError::BadRequest(format!("No GFF3 found in genome '{}'", name)))?;
        let fasta = find_file_with_ext(&subdir, &FASTA_EXTENSIONS);
        let sa = subdir.join("sa");
        let sa_dir = if sa.is_dir() { Some(sa) } else { None };
        return Ok(GenomePaths {
            gff3,
            fasta,
            sa_dir,
        });
    }

    // Check top-level files
    for ext in &GFF3_EXTENSIONS {
        let path = data_dir.join(format!("{}.{}", name, ext));
        if path.exists() {
            return Ok(GenomePaths {
                gff3: path,
                fasta: None,
                sa_dir: None,
            });
        }
    }

    Err(AppError::BadRequest(format!(
        "Genome '{}' not found in data directory",
        name
    )))
}

/// Whether the data directory has an entry called `name`, or a loose GFF3 file
/// for it, compared exactly rather than the way the filesystem would look it up.
fn has_entry_named(data_dir: &Path, name: &str) -> bool {
    let loose: Vec<String> = GFF3_EXTENSIONS
        .iter()
        .map(|ext| format!("{name}.{ext}"))
        .collect();
    std::fs::read_dir(data_dir).is_ok_and(|entries| {
        entries.flatten().any(|entry| {
            let entry_name = entry.file_name();
            let entry_name = entry_name.to_string_lossy();
            entry_name == name || loose.iter().any(|l| *l == entry_name)
        })
    })
}

pub fn find_file_with_ext(dir: &Path, extensions: &[&str]) -> Option<PathBuf> {
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

/// How a genome gets built from its files. A seam for tests, which need a load that takes a
/// known time; the server always uses [`load_context`].
type Loader =
    Arc<dyn Fn(&str, &GenomePaths, u64) -> anyhow::Result<AnnotationContext> + Send + Sync>;

/// What a load has to say so far: nothing while it runs, then the genome or why it could not
/// be built.
type Outcome = Option<Result<Arc<AnnotationContext>, Arc<str>>>;

/// The loads under way or finished, by genome name. Shared with the task that does a load so
/// that a failed one can take itself out, which lets the next request try again.
type Loads = Arc<Mutex<HashMap<String, watch::Receiver<Outcome>>>>;

/// The genomes loaded so far besides the one the server started with.
///
/// A genome is loaded the first time something asks for it and then kept, so
/// memory grows with the genomes actually used, bounded by the number of
/// directories under `--data-dir`. Concurrent first requests for one genome share
/// a single load.
///
/// The load is not owned by the request that started it. It runs as a task of its own, and a
/// request only waits for its outcome, so a client that gives up (the browser tab closes, the
/// proxy times out) does not stop the load or start a second one. A load that belonged to its
/// first request would be restarted by the next waiter each time that request went away, and
/// every restart builds the whole genome again: for human, 2.4 GiB of heap and 4.5 GiB while
/// its transcript cache is rebuilt, which a handful of abandoned requests multiply past what a
/// 16 GB host has.
pub struct GenomeRegistry {
    data_dir: Option<PathBuf>,
    distance: u64,
    loads: Loads,
    loads_started: AtomicUsize,
    loader: Loader,
}

impl GenomeRegistry {
    pub fn new(data_dir: Option<PathBuf>, distance: u64) -> Self {
        Self::with_loader(data_dir, distance, Arc::new(load_context))
    }

    fn with_loader(data_dir: Option<PathBuf>, distance: u64, loader: Loader) -> Self {
        Self {
            data_dir,
            distance,
            loads: Arc::new(Mutex::new(HashMap::new())),
            loads_started: AtomicUsize::new(0),
            loader,
        }
    }

    pub fn data_dir(&self) -> Option<&Path> {
        self.data_dir.as_deref()
    }

    /// How many loads have been started, which for a registry that works is
    /// the number of distinct genomes asked for. A test can assert that
    /// concurrent first requests shared one.
    pub fn loads_started(&self) -> usize {
        self.loads_started.load(Ordering::Relaxed)
    }

    /// The genome if it is already loaded; never waits and never loads.
    pub fn loaded(&self, name: &str) -> Option<Arc<AnnotationContext>> {
        let loads = self.loads.lock().unwrap_or_else(|p| p.into_inner());
        let outcome = loads.get(name)?.borrow().clone();
        outcome?.ok()
    }

    /// The genome, loading it first if this is the first request for it.
    ///
    /// The name is checked against the data directory before the registry
    /// remembers anything about it, so a request cannot grow the map by naming
    /// things that are not genomes.
    pub async fn get_or_load(&self, name: &str) -> Result<Arc<AnnotationContext>, AppError> {
        if let Some(ctx) = self.loaded(name) {
            return Ok(ctx);
        }
        let Some(data_dir) = self.data_dir.as_deref() else {
            return Err(AppError::BadRequest("No data directory configured".into()));
        };
        let paths = resolve_genome_paths(data_dir, name)?;
        let mut outcome = {
            let mut loads = self.loads.lock().unwrap_or_else(|p| p.into_inner());
            match loads.get(name) {
                Some(rx) => rx.clone(),
                None => {
                    let (tx, rx) = watch::channel(None);
                    loads.insert(name.to_string(), rx.clone());
                    self.loads_started.fetch_add(1, Ordering::Relaxed);
                    self.start_load(name.to_string(), paths, tx);
                    rx
                }
            }
        };
        loop {
            let answer = outcome.borrow().clone();
            match answer {
                Some(Ok(ctx)) => return Ok(ctx),
                Some(Err(why)) => {
                    return Err(AppError::Internal(anyhow::anyhow!(
                        "genome {name} could not be loaded: {why}"
                    )))
                }
                None => {}
            }
            outcome.changed().await.map_err(|_| {
                AppError::Internal(anyhow::anyhow!(
                    "the load of genome {name} ended without an answer"
                ))
            })?;
        }
    }

    /// Run the load of one genome on its own task and report the outcome on `tx`.
    fn start_load(&self, name: String, paths: GenomePaths, tx: watch::Sender<Outcome>) {
        let loader = Arc::clone(&self.loader);
        let loads = Arc::clone(&self.loads);
        let distance = self.distance;
        tokio::spawn(async move {
            let label = name.clone();
            let built = tokio::task::spawn_blocking(move || loader(&label, &paths, distance)).await;
            let outcome: Result<Arc<AnnotationContext>, Arc<str>> = match built {
                Ok(Ok(ctx)) => Ok(Arc::new(ctx)),
                Ok(Err(e)) => {
                    tracing::error!("Could not load genome {name}: {e:?}");
                    Err(format!("{e:#}").into())
                }
                Err(e) => {
                    tracing::error!("The load of genome {name} panicked: {e}");
                    Err(format!("{e}").into())
                }
            };
            if outcome.is_err() {
                // So that the next request tries again instead of finding the failure.
                loads
                    .lock()
                    .unwrap_or_else(|p| p.into_inner())
                    .remove(&name);
            }
            tx.send_replace(Some(outcome));
        });
    }
}

/// Build the context of one genome from its own files and nothing else.
///
/// `sa_dir` is the genome's `sa/` directory, never the global `--sa-dir`; see
/// the module documentation for why.
fn load_context(
    name: &str,
    paths: &GenomePaths,
    distance: u64,
) -> anyhow::Result<AnnotationContext> {
    let utf8 = |p: &Path| -> anyhow::Result<String> {
        p.to_str()
            .map(str::to_string)
            .ok_or_else(|| anyhow::anyhow!("path is not valid UTF-8: {:?}", p))
    };
    let gff3 = utf8(&paths.gff3)?;
    let fasta = paths.fasta.as_deref().map(utf8).transpose()?;
    let sa_dir = paths.sa_dir.as_deref().map(utf8).transpose()?;
    let started = Instant::now();
    let ctx = AnnotationContext::new(Some(&gff3), fasta.as_deref(), sa_dir.as_deref(), distance)?;
    tracing::info!(
        "Loaded genome {name}: {} transcripts in {:.1} s",
        ctx.transcript_count(),
        started.elapsed().as_secs_f64()
    );
    Ok(ctx)
}

#[cfg(test)]
mod tests {
    use super::*;
    use axum::response::IntoResponse;
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
        fs::create_dir_all(genome_dir.join("sa")).unwrap();

        let paths = resolve_genome_paths(&data_dir, "hg38").unwrap();
        assert_eq!(paths.gff3, genome_dir.join("annotation.gff3"));
        assert_eq!(paths.fasta, Some(genome_dir.join("genome.fa")));
        // The SA directory comes from the same validated join as the others,
        // rather than being rebuilt from the request name at the call site.
        assert_eq!(paths.sa_dir, Some(genome_dir.join("sa")));
    }

    #[test]
    fn resolve_genome_paths_reports_no_sa_dir_when_there_is_none() {
        let data_dir = tmp_dir("nosa");
        let genome_dir = data_dir.join("hg38");
        fs::create_dir_all(&genome_dir).unwrap();
        fs::write(genome_dir.join("annotation.gff3"), "##gff-version 3\n").unwrap();

        // `None` is what lets `/api/load-genome` fall back to the global
        // `--sa-dir`; an unconditional path would point the loader at a
        // directory that does not exist.
        assert_eq!(
            resolve_genome_paths(&data_dir, "hg38").unwrap().sa_dir,
            None
        );
    }

    #[test]
    fn resolve_genome_paths_finds_top_level_loose_gff3() {
        let data_dir = tmp_dir("toplevel");
        fs::write(data_dir.join("mygenome.gff3"), "##gff-version 3\n").unwrap();

        let paths = resolve_genome_paths(&data_dir, "mygenome").unwrap();
        assert_eq!(paths.gff3, data_dir.join("mygenome.gff3"));
        assert_eq!(paths.fasta, None);
        // A loose file has no directory to hold one.
        assert_eq!(paths.sa_dir, None);
    }

    #[test]
    fn resolve_genome_paths_errors_when_genome_not_found() {
        let data_dir = tmp_dir("missing");
        let err = resolve_genome_paths(&data_dir, "nonexistent").unwrap_err();
        let resp = err.into_response();
        assert_eq!(resp.status(), axum::http::StatusCode::BAD_REQUEST);
    }

    #[test]
    fn a_name_has_to_match_the_directory_letter_for_letter() {
        // On a case-insensitive filesystem the lookup alone would resolve
        // `Alpha` to `alpha`; on a case-sensitive one it is refused anyway. The
        // answer is the same on both, which is the point.
        let data_dir = tmp_dir("exact-name");
        fs::create_dir_all(data_dir.join("alpha")).unwrap();
        fs::write(data_dir.join("alpha").join("a.gff3"), "##gff-version 3\n").unwrap();
        fs::write(data_dir.join("loose.gff3"), "##gff-version 3\n").unwrap();

        assert!(resolve_genome_paths(&data_dir, "alpha").is_ok());
        assert!(resolve_genome_paths(&data_dir, "loose").is_ok());
        for name in ["Alpha", "ALPHA", "Loose"] {
            assert!(resolve_genome_paths(&data_dir, name).is_err(), "{name}");
        }
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

    #[test]
    fn scan_lists_each_genome_once_with_its_gff3_file_name() {
        let data_dir = tmp_dir("scan");
        let human = data_dir.join("hg38_ensembl_115");
        fs::create_dir_all(human.join("sa")).unwrap();
        fs::write(
            human.join("Homo_sapiens.GRCh38.115.gff3"),
            "##gff-version 3\n",
        )
        .unwrap();
        // The sidecar cache next to a GFF3 is the same genome, not a second one.
        fs::write(
            human.join("Homo_sapiens.GRCh38.115.gff3.fastvep.cache"),
            b"",
        )
        .unwrap();
        fs::write(human.join("genome.fa"), b">1\nA\n").unwrap();
        fs::write(human.join("sa").join("clinvar.osa2"), b"").unwrap();
        let mouse = data_dir.join("mouse");
        fs::create_dir_all(&mouse).unwrap();
        fs::write(
            mouse.join("Mus_musculus.GRCm39.115.gff3"),
            "##gff-version 3\n",
        )
        .unwrap();
        // A directory with no gene model is not a genome.
        fs::create_dir_all(data_dir.join("rat")).unwrap();
        // A loose file and a directory of the same name: the directory wins.
        fs::write(data_dir.join("mouse.gff3"), "##gff-version 3\n").unwrap();
        fs::write(data_dir.join("loose.gff3.gz"), b"").unwrap();

        let listed = scan_data_dir(&data_dir);
        let names: Vec<&str> = listed.iter().map(|g| g.name.as_str()).collect();
        assert_eq!(names, ["hg38_ensembl_115", "loose", "mouse"]);
        assert_eq!(listed[0].gff3, "Homo_sapiens.GRCh38.115.gff3");
        assert!(listed[0].has_fasta && listed[0].has_sa);
        assert_eq!(listed[1].gff3, "loose.gff3.gz");
        assert_eq!(listed[2].gff3, "Mus_musculus.GRCm39.115.gff3");
        assert!(!listed[2].has_fasta && !listed[2].has_sa);
    }

    #[test]
    fn a_gene_model_inside_the_data_directory_answers_to_its_directory_name() {
        let data_dir = tmp_dir("default-name");
        let genome_dir = data_dir.join("hg38_ensembl_115");
        fs::create_dir_all(&genome_dir).unwrap();
        let gff3 = genome_dir.join("Homo_sapiens.GRCh38.115.gff3");
        fs::write(&gff3, "##gff-version 3\n").unwrap();

        assert_eq!(
            default_genome_name(gff3.to_str(), Some(data_dir.as_path())),
            "hg38_ensembl_115"
        );
        // Spelled with a `.` segment, it is still the same directory.
        let roundabout = data_dir
            .join(".")
            .join("hg38_ensembl_115")
            .join("Homo_sapiens.GRCh38.115.gff3");
        assert_eq!(
            default_genome_name(roundabout.to_str(), Some(data_dir.as_path())),
            "hg38_ensembl_115"
        );
    }

    #[test]
    fn any_other_gene_model_is_named_for_its_file() {
        let data_dir = tmp_dir("default-name-outside");
        assert_eq!(
            default_genome_name(
                Some("/somewhere/else/mygenes.gff3"),
                Some(data_dir.as_path())
            ),
            "mygenes"
        );
        assert_eq!(
            default_genome_name(Some("x/Homo_sapiens.GRCh38.115.gff3.gz"), None),
            "Homo_sapiens.GRCh38.115"
        );
        assert_eq!(default_genome_name(None, None), "default");
    }

    const TINY_GFF3: &str = "##gff-version 3\n\
1\ttest\tgene\t1\t300\t.\t+\t.\tID=gene:G1;Name=G1;biotype=protein_coding\n\
1\ttest\tmRNA\t1\t300\t.\t+\t.\tID=transcript:T1;Parent=gene:G1;biotype=protein_coding\n\
1\ttest\texon\t1\t300\t.\t+\t.\tID=exon:E1;Parent=transcript:T1;rank=1\n\
1\ttest\tCDS\t1\t300\t.\t+\t0\tID=CDS:P1;Parent=transcript:T1\n";

    fn tiny_genome(label: &str) -> PathBuf {
        let data_dir = tmp_dir(label);
        let genome = data_dir.join("alpha");
        fs::create_dir_all(&genome).unwrap();
        fs::write(genome.join("alpha.gff3"), TINY_GFF3).unwrap();
        data_dir
    }

    /// A load that takes `delay`, and counts how many times it ran.
    fn slow_registry(
        data_dir: PathBuf,
        delay: std::time::Duration,
    ) -> (Arc<GenomeRegistry>, Arc<AtomicUsize>) {
        let ran = Arc::new(AtomicUsize::new(0));
        let counted = Arc::clone(&ran);
        let registry = GenomeRegistry::with_loader(
            Some(data_dir),
            5000,
            Arc::new(move |label, paths, distance| {
                counted.fetch_add(1, Ordering::Relaxed);
                std::thread::sleep(delay);
                load_context(label, paths, distance)
            }),
        );
        (Arc::new(registry), ran)
    }

    // The load belongs to the registry and not to the request that asked first. If it
    // belonged to that request, each one that went away would hand the load to the next
    // waiter, which would build the whole genome again while the abandoned build was still
    // running: a few clients that give up early multiply the memory of a human genome.
    #[tokio::test]
    async fn a_client_that_goes_away_neither_stops_a_load_nor_starts_another() {
        let (registry, ran) =
            slow_registry(tiny_genome("cancel"), std::time::Duration::from_millis(400));
        for _ in 0..6 {
            let r = Arc::clone(&registry);
            let request = tokio::spawn(async move { r.get_or_load("alpha").await.map(|_| ()) });
            tokio::time::sleep(std::time::Duration::from_millis(25)).await;
            request.abort();
            let _ = request.await;
        }
        let ctx = registry
            .get_or_load("alpha")
            .await
            .expect("the genome loads");
        assert_eq!(ctx.transcript_count(), 1);
        assert_eq!(ran.load(Ordering::Relaxed), 1, "the genome was built once");
        assert_eq!(registry.loads_started(), 1);
        assert!(registry.loaded("alpha").is_some());
    }

    #[tokio::test]
    async fn requests_that_arrive_while_a_genome_loads_share_the_load() {
        let (registry, ran) =
            slow_registry(tiny_genome("share"), std::time::Duration::from_millis(200));
        let waiting: Vec<_> = (0..8)
            .map(|_| {
                let r = Arc::clone(&registry);
                tokio::spawn(
                    async move { r.get_or_load("alpha").await.map(|c| c.transcript_count()) },
                )
            })
            .collect();
        for w in waiting {
            assert_eq!(w.await.unwrap().unwrap(), 1);
        }
        assert_eq!(ran.load(Ordering::Relaxed), 1);
    }

    // A genome that failed to load is not a genome that fails for ever: the failure is not
    // kept, so the next request builds it again, once whatever was wrong is fixed.
    #[tokio::test]
    async fn a_failed_load_is_not_remembered() {
        let data_dir = tiny_genome("retry");
        let attempts = Arc::new(AtomicUsize::new(0));
        let counted = Arc::clone(&attempts);
        let registry = GenomeRegistry::with_loader(
            Some(data_dir),
            5000,
            Arc::new(move |label, paths, distance| {
                if counted.fetch_add(1, Ordering::Relaxed) == 0 {
                    anyhow::bail!("the disk was not there yet");
                }
                load_context(label, paths, distance)
            }),
        );
        let first = registry.get_or_load("alpha").await;
        assert!(matches!(first, Err(AppError::Internal(_))));
        assert!(registry.loaded("alpha").is_none());
        let second = registry
            .get_or_load("alpha")
            .await
            .expect("the second try loads");
        assert_eq!(second.transcript_count(), 1);
        assert_eq!(registry.loads_started(), 2);
    }
}
