//! `fastvep cache`: build a binary transcript cache from GFF3 + FASTA.

use super::open_vcf_input_reader;
use anyhow::{Context, Result};
use fastvep_cache::fasta::FastaReader;
use fastvep_cache::gff::parse_gff3_with_source;
use fastvep_cache::providers::{FastaSequenceProvider, SequenceProvider};
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead};
use std::path::Path;

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
        return Err(anyhow::anyhow!(
            "`fastvep cache` requires at least one --gff3"
        ));
    }
    let specs: Vec<Gff3Spec> = gff3_paths.iter().map(|s| parse_gff3_arg(s)).collect();
    let mut transcripts: Vec<fastvep_genome::Transcript> = Vec::new();
    for spec in &specs {
        let gff_file =
            File::open(&spec.path).with_context(|| format!("Opening GFF3 file: {}", spec.path))?;
        // Auto-decompress .gz / .bgz GFF3 inputs. Without this we'd silently
        // produce a 0-transcript cache.
        let trs = if spec.path.ends_with(".gz") || spec.path.ends_with(".bgz") {
            parse_gff3_with_source(flate2::read::MultiGzDecoder::new(gff_file), &spec.source)?
        } else {
            parse_gff3_with_source(gff_file, &spec.source)?
        };
        eprintln!(
            "Loaded {} transcripts from {} (source label: {})",
            trs.len(),
            spec.path,
            spec.source
        );
        transcripts.extend(trs);
    }
    if transcripts.is_empty() {
        return Err(anyhow::anyhow!(
            "GFF3 source(s) [{}] produced 0 transcripts — refusing to write an empty cache.",
            specs
                .iter()
                .map(|s| s.path.as_str())
                .collect::<Vec<_>>()
                .join(", ")
        ));
    }
    if specs.len() > 1 {
        eprintln!(
            "Merged {} GFF3 sources into {} total transcripts",
            specs.len(),
            transcripts.len()
        );
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
        let fasta_file =
            File::open(fasta).with_context(|| format!("Opening FASTA file: {}", fasta))?;
        let reader = FastaReader::from_reader(fasta_file)?;
        eprintln!("Loaded reference FASTA from {}", fasta);

        // Canonicalize every transcript (and its gene) to the FASTA's contig
        // naming. This both lets sequence fetch succeed and makes a merged
        // Ensembl + RefSeq cache use one consistent naming scheme, so
        // `annotate` matches a VCF regardless of which GFF3 the transcript
        // came from. RefSeq accessions only resolve when a synonyms file maps
        // them; chr↔bare / mito resolve with no synonyms file at all.
        let contigs: std::collections::HashSet<String> = reader
            .sequence_names()
            .into_iter()
            .map(|s| s.to_string())
            .collect();
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
                    eprintln!(
                        "Warning: could not build sequences for {}: {}",
                        tr.stable_id, e
                    );
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

/// Load one GFF3 source, returning the transcripts and whether the load was
/// restricted to the input VCF's variant regions.
///
/// The `bool` is not cosmetic: a region-restricted set is only valid for the
/// VCF whose regions produced it, so callers must not persist it as a cache
/// keyed on the GFF3 path alone.
pub(crate) fn load_one_gff3(
    spec: &Gff3Spec,
    vcf_input: &str,
    distance: u64,
) -> Result<(Vec<fastvep_genome::Transcript>, bool)> {
    let gff_path = Path::new(&spec.path);
    let tbi_path = format!("{}.tbi", spec.path);

    if spec.path.ends_with(".gz") && Path::new(&tbi_path).exists() {
        let regions = prescan_vcf_regions(vcf_input, distance)?;
        eprintln!(
            "Pre-scanned {} variant regions for {}",
            regions.len(),
            spec.path
        );
        let trs =
            fastvep_cache::gff::parse_gff3_indexed_with_source(gff_path, &regions, &spec.source)?;
        Ok((trs, true))
    } else {
        let gff_file =
            File::open(&spec.path).with_context(|| format!("Opening GFF3 file: {}", spec.path))?;
        let trs = if spec.path.ends_with(".gz") || spec.path.ends_with(".bgz") {
            parse_gff3_with_source(flate2::read::MultiGzDecoder::new(gff_file), &spec.source)?
        } else {
            parse_gff3_with_source(gff_file, &spec.source)?
        };
        Ok((trs, false))
    }
}

fn prescan_vcf_regions(vcf_path: &str, distance: u64) -> Result<Vec<(String, u64, u64)>> {
    let input_reader = open_vcf_input_reader(vcf_path)
        .with_context(|| format!("Pre-scanning VCF: {}", vcf_path))?;
    let reader = io::BufReader::new(input_reader);

    let mut regions: HashMap<String, (u64, u64)> = HashMap::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') || line.is_empty() {
            continue;
        }
        let mut fields = line.split('\t');
        let chrom = match fields.next() {
            Some(c) => c.to_string(),
            None => continue,
        };
        let pos: u64 = match fields.next().and_then(|p| p.parse().ok()) {
            Some(p) => p,
            None => continue,
        };
        let start = pos.saturating_sub(distance);
        let end = pos + distance;

        let entry = regions.entry(chrom).or_insert((start, end));
        entry.0 = entry.0.min(start);
        entry.1 = entry.1.max(end);
    }

    Ok(regions
        .into_iter()
        .map(|(chrom, (s, e))| (chrom, s, e))
        .collect())
}

// =============================================================================
// SA Build: Build supplementary annotation databases from source VCFs
// =============================================================================

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
