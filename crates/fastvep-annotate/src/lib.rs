//! Shared annotation engine for fastVEP.
//!
//! Provides [`AnnotationContext`] which loads transcript models, reference sequences,
//! and supplementary annotation sources, then annotates VCF variants against them.
//!
//! Used by both `fastvep-web` (production HTTP server) and `fastvep-cli` (embedded web server).
//! The CLI batch pipeline (`run_annotate`) has its own streaming implementation but shares
//! the same underlying crates.

mod hgvs_normalize;
pub mod pick;

pub use hgvs_normalize::{
    convert_ins_to_dup_range, convert_ins_to_dup_range_noncoding, hgvsc_intronic_shifted,
    intronic_dup_span, intronic_ins_as_dup, three_prime_shift_intronic,
};

use anyhow::{Context, Result};
use fastvep_cache::annotation::{AnnotationProvider, AnnotationValue};
use fastvep_cache::fasta::FastaReader;
use fastvep_cache::gff::parse_gff3;
use fastvep_cache::providers::{
    FastaSequenceProvider, IndexedTranscriptProvider, SequenceProvider, TranscriptProvider,
};
use fastvep_consequence::ConsequencePredictor;
use fastvep_core::{Allele, Consequence};
use fastvep_genome::Transcript;
use fastvep_io::output;
use fastvep_io::variant::{AlleleAnnotation, TranscriptVariation, VariationFeature};
use fastvep_io::vcf::VcfParser;
use pick::{has_transcripts_to_pick, pick_best_transcript_idx_with, DEFAULT_PICK_ORDER};
use rayon::prelude::*;
use std::fs::File;
use std::path::Path;
use std::sync::atomic::{AtomicUsize, Ordering};

/// Pre-loaded annotation context shared by web and CLI.
///
/// Holds transcript models, a reference sequence provider, a consequence predictor,
/// and supplementary annotation providers (ClinVar, gnomAD, etc.).
pub struct AnnotationContext {
    pub transcript_provider: IndexedTranscriptProvider,
    pub seq_provider: Option<Box<dyn SequenceProvider + Send + Sync>>,
    pub predictor: ConsequencePredictor,
    pub gff3_source: Option<String>,
    pub distance: u64,
    pub hgvs: bool,
    /// Supplementary annotation providers (ClinVar, gnomAD, etc.).
    ///
    /// Held unwrapped: every `AnnotationProvider` method takes `&self` and the
    /// trait requires `Send + Sync`, so concurrent lookups need no external
    /// lock — the readers own the synchronization for their internal caches.
    /// These used to be `Mutex<Box<dyn AnnotationProvider>>`, which serialized
    /// every lookup of every provider across concurrent `/annotate` requests
    /// (and which the CLI pipeline worked around by immediately unwrapping).
    pub sa_providers: Vec<Box<dyn AnnotationProvider>>,
    /// Gene-level annotation providers (OMIM, gnomAD gene constraints, ClinVar protein index).
    pub gene_providers: Vec<fastvep_sa::gene::GeneIndex>,
    /// ACMG-AMP classification configuration (None = disabled).
    pub acmg_config: Option<fastvep_classification::AcmgConfig>,
    /// Curated functional-assay evidence supplying PS3/BS3 (None = not provided).
    ///
    /// PS3 and BS3 are the two criteria that cannot be computed from variant
    /// data. Without this the classifier reports them NotEvaluated, which is
    /// the honest answer rather than a silent absence of evidence.
    pub functional_evidence: Option<fastvep_classification::FunctionalEvidenceIndex>,
}

impl AnnotationContext {
    /// Build a context from GFF3, optional FASTA, and optional SA directory.
    pub fn new(
        gff3: Option<&str>,
        fasta: Option<&str>,
        sa_dir: Option<&str>,
        distance: u64,
    ) -> Result<Self> {
        let gff3_source: Option<String> = gff3.map(|p| {
            Path::new(p)
                .file_name()
                .map(|n| n.to_string_lossy().to_string())
                .unwrap_or_else(|| p.to_string())
        });

        let mut transcripts = if let Some(gff3_path) = gff3 {
            let cache_path =
                fastvep_cache::transcript_cache::default_cache_path(Path::new(gff3_path));
            let from_cache = if cache_path.exists() {
                let is_fresh = fastvep_cache::transcript_cache::cache_is_fresh(
                    &cache_path,
                    Path::new(gff3_path),
                );
                if is_fresh {
                    fastvep_cache::transcript_cache::load_cache(&cache_path).ok()
                } else {
                    None
                }
            } else {
                None
            };
            if let Some(trs) = from_cache {
                tracing::info!("Loaded {} transcripts from cache", trs.len());
                trs
            } else {
                let gff_file = File::open(gff3_path)
                    .with_context(|| format!("Opening GFF3 file: {}", gff3_path))?;
                // Auto-decompress gzipped GFF3. Without this, parse_gff3
                // reads binary gz bytes as text, yields zero transcripts,
                // and downstream silently produces empty annotations.
                let trs = if gff3_path.ends_with(".gz") || gff3_path.ends_with(".bgz") {
                    parse_gff3(flate2::read::MultiGzDecoder::new(gff_file))?
                } else {
                    parse_gff3(gff_file)?
                };
                if trs.is_empty() {
                    return Err(anyhow::anyhow!(
                        "GFF3 file {} produced 0 transcripts — likely malformed, truncated, or unrecognized format. Refusing to continue with empty transcript set.",
                        gff3_path
                    ));
                }
                tracing::info!("Loaded {} transcripts from {}", trs.len(), gff3_path);
                if let Err(e) = fastvep_cache::transcript_cache::save_cache(&trs, &cache_path) {
                    tracing::warn!("Could not save cache: {}", e);
                }
                trs
            }
        } else {
            Vec::new()
        };

        let seq_provider: Option<Box<dyn SequenceProvider + Send + Sync>> =
            if let Some(fasta_path) = fasta {
                let fai_path = format!("{}.fai", fasta_path);
                if Path::new(&fai_path).exists() {
                    let reader =
                        fastvep_cache::fasta::MmapFastaReader::open(Path::new(fasta_path))?;
                    tracing::info!("Memory-mapped FASTA from {}", fasta_path);
                    Some(Box::new(
                        fastvep_cache::providers::MmapFastaSequenceProvider::new(reader),
                    ))
                } else {
                    let fasta_file = File::open(fasta_path)
                        .with_context(|| format!("Opening FASTA: {}", fasta_path))?;
                    let reader = FastaReader::from_reader(fasta_file)?;
                    tracing::info!("Loaded FASTA from {}", fasta_path);
                    Some(Box::new(FastaSequenceProvider::new(reader)))
                }
            } else {
                None
            };

        // A non-coding transcript needs its spliced sequence only for HGVS - the
        // 3'-rule and `dup` collapsing are read off it - and building it for
        // every lncRNA and retained-intron transcript roughly doubles this step,
        // so it is built here where HGVS is always on and gated in the CLI.
        if let Some(ref sp) = seq_provider {
            let built = AtomicUsize::new(0);
            transcripts.par_iter_mut().for_each(|tr| {
                if tr.spliced_seq.is_none()
                    && tr
                        .build_sequences(|chrom, start, end| {
                            sp.fetch_sequence(chrom, start, end)
                                .map_err(|e| e.to_string())
                        })
                        .is_ok()
                {
                    built.fetch_add(1, Ordering::Relaxed);
                }
            });
            let built = built.load(Ordering::Relaxed);
            tracing::info!("Built sequences for {} transcripts", built);

            // Re-save cache with pre-built sequences so future startups skip this step
            if built > 0 {
                if let Some(gff3_path) = gff3 {
                    let cache_path =
                        fastvep_cache::transcript_cache::default_cache_path(Path::new(gff3_path));
                    match fastvep_cache::transcript_cache::save_cache(&transcripts, &cache_path) {
                        Ok(()) => tracing::info!("Updated cache with pre-built sequences"),
                        Err(e) => tracing::warn!("Could not update cache with sequences: {}", e),
                    }
                }
            }
        }

        let transcript_provider = IndexedTranscriptProvider::new(transcripts);
        let predictor = ConsequencePredictor::new(distance, distance);

        // Load supplementary annotation providers (.osa, .osa2 files)
        let sa_providers = if let Some(dir) = sa_dir {
            load_sa_providers(Path::new(dir))?
        } else {
            Vec::new()
        };

        let gene_providers = if let Some(dir) = sa_dir {
            load_gene_providers(Path::new(dir))?
        } else {
            Vec::new()
        };

        Ok(Self {
            functional_evidence: None,
            transcript_provider,
            seq_provider,
            predictor,
            gff3_source,
            distance,
            hgvs: true,
            sa_providers,
            gene_providers,
            acmg_config: None,
        })
    }

    pub fn transcript_count(&self) -> usize {
        self.transcript_provider.transcript_count()
    }

    /// Whether a repeat-interval source (RepeatMasker or equivalent) is loaded.
    ///
    /// BP3 needs this to tell "this variant is not in a repeat" from "no repeat
    /// database was loaded". An interval source only produces an annotation on
    /// a hit, so those two are indistinguishable from the annotation list alone,
    /// and reading a miss as the latter makes BP3 report a setup error on every
    /// in-frame indel outside a repeat.
    ///
    /// Matched on the same substrings the classifier uses to find the track, so
    /// a source named such that this returns `true` is exactly a source BP3 will
    /// actually read.
    pub fn repeat_db_loaded(&self) -> bool {
        self.sa_providers.iter().any(|p| {
            let k = p.json_key().to_lowercase();
            k.contains("repeat") || k.contains("repeatmasker") || k.contains("simple_repeat")
        })
    }

    /// Names of loaded supplementary annotation sources.
    pub fn sa_source_names(&self) -> Vec<String> {
        self.sa_providers
            .iter()
            .map(|p| p.name().to_string())
            .collect()
    }

    /// Load a genome from GFF3 path (+ optional FASTA + optional SA directory).
    /// Replaces transcripts, sequence provider, and SA providers.
    pub fn load_genome(
        &mut self,
        gff3_path: &str,
        fasta_path: Option<&str>,
        sa_dir: Option<&str>,
    ) -> Result<usize> {
        let new_ctx = Self::new(Some(gff3_path), fasta_path, sa_dir, self.distance)?;
        let tr_count = new_ctx.transcript_provider.transcript_count();
        self.transcript_provider = new_ctx.transcript_provider;
        self.seq_provider = new_ctx.seq_provider;
        self.predictor = new_ctx.predictor;
        self.gff3_source = new_ctx.gff3_source;
        self.sa_providers = new_ctx.sa_providers;
        self.gene_providers = new_ctx.gene_providers;
        Ok(tr_count)
    }

    /// Replace the transcript models by parsing GFF3 text uploaded from the browser.
    pub fn update_gff3_text(&mut self, gff3_text: &str) -> Result<(usize, usize)> {
        let mut transcripts = parse_gff3(gff3_text.as_bytes())?;
        let gene_count = {
            let mut genes = std::collections::HashSet::new();
            for t in &transcripts {
                genes.insert(t.gene.stable_id.clone());
            }
            genes.len()
        };

        if let Some(ref sp) = self.seq_provider {
            let mut built = 0usize;
            for tr in &mut transcripts {
                if tr.spliced_seq.is_none()
                    && tr
                        .build_sequences(|chrom, start, end| {
                            sp.fetch_sequence(chrom, start, end)
                                .map_err(|e| e.to_string())
                        })
                        .is_ok()
                {
                    built += 1;
                }
            }
            if built > 0 {
                tracing::info!("Built sequences for {} transcripts", built);
            }
        }

        let tr_count = transcripts.len();
        self.transcript_provider = IndexedTranscriptProvider::new(transcripts);
        self.gff3_source = Some("user-upload".to_string());
        tracing::info!(
            "Updated GFF3: {} genes, {} transcripts",
            gene_count,
            tr_count
        );
        Ok((gene_count, tr_count))
    }

    /// Annotate VCF text and return JSON results, using `self.acmg_config`.
    pub fn annotate_vcf_text(&self, vcf_text: &str, pick: bool) -> Result<Vec<serde_json::Value>> {
        self.annotate_vcf_text_with_acmg(vcf_text, pick, self.acmg_config.as_ref())
    }

    /// Annotate VCF text and return JSON results, using an explicit ACMG
    /// config for this call instead of `self.acmg_config`.
    ///
    /// This only needs `&self` (no mutation of shared context state), so
    /// callers serving concurrent requests over a shared `AnnotationContext`
    /// can take a read lock instead of a write lock — see fastvep-web's
    /// `annotate` handler, which used to mutate `self.acmg_config` per
    /// request under a write lock, serializing all concurrent traffic
    /// (including unrelated reads) behind whichever annotation was running.
    pub fn annotate_vcf_text_with_acmg(
        &self,
        vcf_text: &str,
        pick: bool,
        acmg_config: Option<&fastvep_classification::AcmgConfig>,
    ) -> Result<Vec<serde_json::Value>> {
        // Computed once: BP3 must distinguish "not in a repeat" from "no repeat
        // database loaded", and only the context knows which sources are open.
        let repeat_db_loaded = self.repeat_db_loaded();
        let mut vcf_parser = VcfParser::new(vcf_text.as_bytes())?;

        // Extract sample names from VCF #CHROM header
        let sample_names: Vec<String> = vcf_parser
            .header_lines()
            .last()
            .filter(|l| l.starts_with("#CHROM"))
            .map(|l| l.split('\t').skip(9).map(|s| s.to_string()).collect())
            .unwrap_or_default();

        let mut variants = vcf_parser.read_all()?;

        // gnomAD stores left-aligned, parsimonious alleles, so its queries are
        // normalized against the reference (see the SA block below). Whether that
        // applies depends only on the loaded context, not on the variant, so it
        // is resolved once here instead of re-scanning every provider for every
        // variant.
        let has_gnomad = self.seq_provider.is_some()
            && self.sa_providers.iter().any(|sa| sa.json_key() == "gnomad");

        let functional_evidence = self.functional_evidence.as_ref();

        for vf in &mut variants {
            let chrom = &vf.position.chromosome;
            let query_start = vf.position.start.saturating_sub(self.distance).max(1);
            let query_end = vf.position.end + self.distance;
            let overlapping = self
                .transcript_provider
                .get_transcripts(chrom, query_start, query_end)
                .unwrap_or_default();

            if overlapping.is_empty() {
                annotate_intergenic(vf);
            } else {
                let ref_seq = self
                    .seq_provider
                    .as_ref()
                    .and_then(|sp| sp.fetch_sequence(chrom, query_start, query_end).ok());

                let result = self.predictor.predict(
                    &vf.position,
                    &vf.ref_allele,
                    &vf.alt_alleles,
                    &overlapping,
                    ref_seq.as_deref(),
                );

                for (i, tc) in result.transcript_consequences.iter().enumerate() {
                    // `predict` returns one consequence per transcript it was
                    // handed, in that order, so index `i` is the match. The
                    // scan is the fallback for a predictor that ever stops
                    // promising that; without it the lookup cost grew with the
                    // square of the transcripts overlapping the variant.
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
                                existing_variation: vec![],
                                sift: None,
                                polyphen: None,
                                supplementary: Vec::new(),
                                acmg_classification: None,
                            };

                            if self.hgvs {
                                ann.hgvsg = Some(fastvep_hgvs::hgvsg(
                                    chrom,
                                    vf.position.start,
                                    vf.position.end,
                                    &vf.ref_allele,
                                    &ac.allele,
                                ));
                                if let Some(tr) = transcript {
                                    let versioned_tid = match tr.version {
                                        Some(v) => format!("{}.{}", tc.transcript_id, v),
                                        None => tc.transcript_id.to_string(),
                                    };
                                    let (hgvs_ref, hgvs_alt) =
                                        if tr.strand == fastvep_core::Strand::Reverse {
                                            (
                                                reverse_complement_allele(&vf.ref_allele),
                                                reverse_complement_allele(&ac.allele),
                                            )
                                        } else {
                                            (vf.ref_allele.clone(), ac.allele.clone())
                                        };
                                    if let Some(coding_start) = tr.cdna_coding_start {
                                        if let (Some(cs), Some(ce)) = (ac.cdna_start, ac.cdna_end) {
                                            let (cs, ce) = (cs.min(ce), cs.max(ce));
                                            ann.hgvsc = fastvep_hgvs::hgvsc_with_seq(
                                                &versioned_tid,
                                                cs,
                                                ce,
                                                &hgvs_ref,
                                                &hgvs_alt,
                                                coding_start,
                                                tr.cdna_coding_end,
                                                tr.spliced_seq.as_deref(),
                                                tr.codon_table_start_phase,
                                            );
                                        } else {
                                            // Not "is this intronic": a variant with
                                            // one end in an exon and the other in an
                                            // intron has no exonic cDNA pair, and
                                            // `intron_at` reads its first base only.
                                            // `hgvsc_intronic_shifted` returns
                                            // `None` for anything it cannot place,
                                            // which is the real guard.
                                            ann.hgvsc = hgvsc_intronic_shifted(
                                                self.seq_provider
                                                    .as_deref()
                                                    .map(|sp| sp as &dyn SequenceProvider),
                                                chrom,
                                                tr,
                                                &versioned_tid,
                                                vf.position.start,
                                                vf.position.end,
                                                &vf.ref_allele,
                                                &ac.allele,
                                                &hgvs_ref,
                                                &hgvs_alt,
                                                Some(coding_start),
                                                tr.cdna_coding_end,
                                            );
                                        }
                                    } else if let (Some(cs), Some(ce)) =
                                        (ac.cdna_start, ac.cdna_end)
                                    {
                                        ann.hgvsc = fastvep_hgvs::hgvsc_noncoding(
                                            &versioned_tid,
                                            cs,
                                            ce,
                                            &hgvs_ref,
                                            &hgvs_alt,
                                            tr.spliced_seq.as_deref(),
                                        );
                                    } else {
                                        // Not "is this intronic": a variant with
                                        // one end in an exon and the other in an
                                        // intron has no exonic cDNA pair, and
                                        // `intron_at` reads its first base only.
                                        // `hgvsc_intronic_shifted` returns
                                        // `None` for anything it cannot place,
                                        // which is the real guard.
                                        ann.hgvsc = hgvsc_intronic_shifted(
                                            self.seq_provider
                                                .as_deref()
                                                .map(|sp| sp as &dyn SequenceProvider),
                                            chrom,
                                            tr,
                                            &versioned_tid,
                                            vf.position.start,
                                            vf.position.end,
                                            &vf.ref_allele,
                                            &ac.allele,
                                            &hgvs_ref,
                                            &hgvs_alt,
                                            None,
                                            None,
                                        );
                                    }

                                    // HGVSp
                                    if let (Some(ref aa), Some(ps)) =
                                        (&ac.amino_acids, ac.protein_start)
                                    {
                                        // The residue pair is strand-ordered:
                                        // which of the two is the lower one
                                        // depends on the strand, and
                                        // `hgvsp_inframe_indel` needs both to
                                        // place an insertion between two
                                        // residues. A variant on an exon's edge
                                        // has one end outside the CDS and so may
                                        // have only one coordinate; using it for
                                        // both is what the single one meant
                                        // before.
                                        let pe = ac.protein_end.unwrap_or(ps);
                                        if let Some(ref pid) = tr.protein_id {
                                            let versioned_pid: String = match tr.protein_version {
                                                Some(v) => {
                                                    let suffix = format!(".{}", v);
                                                    if pid.ends_with(suffix.as_str()) {
                                                        pid.clone()
                                                    } else {
                                                        format!("{}.{}", pid, v)
                                                    }
                                                }
                                                None => pid.clone(),
                                            };
                                            let is_fs = ac
                                                .consequences
                                                .contains(&Consequence::FrameshiftVariant);
                                            if is_fs {
                                                if let (Some(spliced), Some(coding_start)) = (
                                                    tr.spliced_seq.as_deref(),
                                                    tr.cdna_coding_start,
                                                ) {
                                                    ann.hgvsp = cds_and_downstream(
                                                        tr,
                                                        spliced,
                                                        coding_start,
                                                    )
                                                    .and_then(|cds| {
                                                        fastvep_hgvs::hgvsp_frameshift_from_cds(
                                                            &versioned_pid,
                                                            &cds,
                                                            ac.cds_start,
                                                            ac.cds_end,
                                                            &vf.ref_allele,
                                                            &ac.allele,
                                                            tr.strand,
                                                            &frameshift_codon_table(tr),
                                                        )
                                                    });
                                                }
                                            } else if aa.1 == "-"
                                                || aa.0.len() != aa.1.len()
                                                || ac
                                                    .consequences
                                                    .contains(&Consequence::StartLost)
                                                || ac
                                                    .consequences
                                                    .contains(&Consequence::InframeDeletion)
                                                || ac
                                                    .consequences
                                                    .contains(&Consequence::InframeInsertion)
                                                // A window wider than one residue also routes
                                                // here even when the lengths match: `hgvsp()`
                                                // compares one residue per side, so `EP/ET`
                                                // reads as unchanged and renders `p.Glu153=`
                                                // for a change VEP calls `p.Pro154Thr` - and a
                                                // synonymous two-residue window needs its whole
                                                // span named.
                                            || aa.0.len() > 1
                                            {
                                                // In-frame indel / delins (frameshift
                                                // handled above). aa.0 holds the replaced
                                                // residues, aa.1 the replacement ("-" for a
                                                // pure deletion).
                                                //
                                                // Insertions must route here too, not to
                                                // `hgvsp()`: that compares only the first
                                                // residue of each side, so `W/WR` reads as
                                                // unchanged and renders `p.Trp185=` for a
                                                // variant that lengthens the protein.
                                                //
                                                // The residue counts decide, not the SO term.
                                                // A delins that replaces residues earns
                                                // `protein_altering_variant` or `stop_gained`
                                                // rather than either in-frame term, and
                                                // keying on the term alone dropped HGVSp
                                                // entirely for every one of them - 1,560
                                                // rows over the ClinVar 2-star in-frame
                                                // delins, where VEP names
                                                // `p.Lys666delinsAsnSer`.
                                                //
                                                // `start_lost` routes here whatever
                                                // its residue counts. Losing the
                                                // initiator makes the downstream
                                                // protein unknowable, so Ensembl
                                                // writes `p.Met1?` and not the
                                                // substitution that `hgvsp()` would
                                                // name - `p.Met1Asn` claims a protein
                                                // with an Asn at residue 1, and there
                                                // may be no protein at all.
                                                //
                                                // The peptide lets `hgvsp_inframe_indel`
                                                // apply the HGVS 3'-rule and collapse a
                                                // repeat to `dup`, matching Ensembl VEP; it
                                                // degrades to the unshifted description
                                                // without one.
                                                // `tr.strand` says which end
                                                // of `aa.0` the `ps` above
                                                // names: on the reverse strand a
                                                // shrinking change arrives
                                                // anchored at the end of its
                                                // span (#89, #96).
                                                ann.hgvsp = fastvep_hgvs::hgvsp_inframe_indel(
                                                    &versioned_pid,
                                                    ps,
                                                    pe,
                                                    &aa.0,
                                                    &aa.1,
                                                    tr.peptide.as_deref().map(str::as_bytes),
                                                    tr.strand,
                                                );
                                            } else {
                                                let ref_aa_byte =
                                                    aa.0.as_bytes()
                                                        .first()
                                                        .copied()
                                                        .unwrap_or(b'X');
                                                let alt_aa_byte =
                                                    aa.1.as_bytes()
                                                        .first()
                                                        .copied()
                                                        .unwrap_or(b'X');
                                                ann.hgvsp = fastvep_hgvs::hgvsp(
                                                    &versioned_pid,
                                                    ps,
                                                    ref_aa_byte,
                                                    alt_aa_byte,
                                                    false,
                                                );
                                            }
                                        }
                                    }
                                }
                            }
                            ann
                        })
                        .collect();

                    // Collect every transcript here; `pick` runs as a single
                    // post-pass below so it can compare all candidates, the
                    // same way `fastvep annotate --pick` does.
                    vf.transcript_variations.push(TranscriptVariation {
                        transcript_id: tc.transcript_id.clone(),
                        gene_id: tc.gene_id.clone(),
                        gene_symbol: tc.gene_symbol.clone(),
                        biotype: tc.biotype.clone(),
                        allele_annotations,
                        canonical: tc.canonical,
                        strand: tc.strand,
                        source: self.gff3_source.clone(),
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
            }

            // Apply `pick` before SA/ACMG so those passes only run on the
            // surviving transcript, and before `compute_most_severe` so the
            // summary term describes the transcript actually reported.
            // Shares `pick::` with `fastvep annotate --pick`: this path used to
            // keep "canonical, or the first transcript seen", which returned a
            // non-canonical first row next to the canonical one for a
            // single-gene variant.
            if pick && has_transcripts_to_pick(&vf.transcript_variations) {
                if let Some(idx) =
                    pick_best_transcript_idx_with(&vf.transcript_variations, DEFAULT_PICK_ORDER)
                {
                    vf.transcript_variations = vec![vf.transcript_variations.swap_remove(idx)];
                }
            }

            // Supplementary annotation: query SA providers once per unique
            // allele, then attach the result to each (transcript, allele)
            // slot. SA results don't depend on the transcript, so this avoids
            // the T× amplification across overlapping transcripts.
            if !self.sa_providers.is_empty() {
                let chrom = &vf.position.chromosome;
                let ref_str = vf.ref_allele.to_string();
                // Normalize the gnomAD query to its minimal representation so
                // indels (especially in repeats) match instead of silently
                // missing — which otherwise makes PM2 misfire on common
                // variants. Applied only to the gnomAD lookup; every other
                // source keeps the query unchanged.
                let mut allele_results: std::collections::HashMap<String, Vec<(String, String)>> =
                    std::collections::HashMap::new();
                for tv in &vf.transcript_variations {
                    for aa in &tv.allele_annotations {
                        let alt_str = aa.allele.to_string();
                        if allele_results.contains_key(&alt_str) {
                            continue;
                        }
                        let gnomad_norm = if has_gnomad {
                            self.seq_provider.as_ref().map(|sp| {
                                fastvep_cache::normalize::normalize_variant(
                                    &**sp,
                                    chrom,
                                    vf.position.start,
                                    &ref_str,
                                    &alt_str,
                                )
                            })
                        } else {
                            None
                        };
                        let mut results: Vec<(String, String)> = Vec::new();
                        for sa in &self.sa_providers {
                            let (q_pos, q_ref, q_alt) = if sa.json_key() == "gnomad" {
                                match &gnomad_norm {
                                    Some(n) => {
                                        (n.pos, n.ref_allele.as_str(), n.alt_allele.as_str())
                                    }
                                    None => (vf.position.start, ref_str.as_str(), alt_str.as_str()),
                                }
                            } else {
                                (vf.position.start, ref_str.as_str(), alt_str.as_str())
                            };
                            // A failed lookup is not a miss: see `SaLookupErrors`.
                            match sa.annotate_position(chrom, q_pos, q_ref, q_alt) {
                                Ok(Some(ann)) => {
                                    let json_str = match ann {
                                        AnnotationValue::Json(j) => j,
                                        AnnotationValue::Positional(j) => j,
                                        AnnotationValue::Interval(v) => {
                                            format!("[{}]", v.join(","))
                                        }
                                    };
                                    results.push((sa.json_key().to_string(), json_str));
                                }
                                Ok(None) => {}
                                Err(e) => sa_lookup_errors().record(sa.json_key(), &e),
                            }
                        }
                        allele_results.insert(alt_str, results);
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
            if !self.gene_providers.is_empty() {
                use fastvep_cache::annotation::GeneAnnotationProvider;
                let mut seen_genes = std::collections::HashSet::new();
                for tv in &vf.transcript_variations {
                    if let Some(gene_sym) = tv.gene_symbol.as_deref() {
                        if seen_genes.insert(gene_sym.to_string()) {
                            for gp in &self.gene_providers {
                                if let Ok(Some(json)) = gp.annotate_gene(gene_sym) {
                                    vf.gene_annotations.push(fastvep_core::GeneAnnotation {
                                        gene_symbol: gene_sym.to_string(),
                                        json_key: gp.json_key().to_string(),
                                        json_string: json,
                                    });
                                }
                            }
                        }
                    }
                }
            }

            // ACMG-AMP classification pass (after all SA annotations are attached)
            if let Some(acmg_cfg) = acmg_config {
                // Parse sample genotypes if trio config is present
                let trio_genotypes = extract_trio_genotypes(vf, acmg_cfg, &sample_names);

                let functional_by_alt = resolve_functional_by_alt(functional_evidence, vf);
                // Resolved here rather than inside the classifier: the ClinVar
                // splice index is keyed by genomic coordinate, and
                // `extract_classification_input` is handed a transcript-level
                // view that carries none.
                let query_alleles = vf.query_alleles();

                for tv in &mut vf.transcript_variations {
                    let gene_sym = tv.gene_symbol.as_deref().unwrap_or("");
                    let gene_anns: Vec<&fastvep_core::GeneAnnotation> = vf
                        .gene_annotations
                        .iter()
                        .filter(|ga| ga.gene_symbol == gene_sym)
                        .collect();
                    for aa in &mut tv.allele_annotations {
                        let alt_idx = vf.alt_alleles.iter().position(|a| *a == aa.allele);
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
                            splice_ps1_evidence(aa, &gene_anns, &query_alleles, alt_idx),
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
                            trio_genotypes.0.clone(),
                            trio_genotypes.1.clone(),
                            trio_genotypes.2.clone(),
                            vec![], // companion_variants populated in second pass
                        );
                        let result = fastvep_classification::classify(&input, acmg_cfg);
                        aa.acmg_classification = serde_json::to_value(&result).ok();
                    }
                }
            }

            vf.compute_most_severe();
        }

        // Compound-het enrichment pass: re-evaluate PM3/BP2 with companion variant data
        if let Some(acmg_cfg) = acmg_config {
            if acmg_cfg.trio.is_some() {
                enrich_compound_het(
                    &mut variants,
                    acmg_cfg,
                    &sample_names,
                    functional_evidence,
                    repeat_db_loaded,
                );
            }
        }

        Ok(variants
            .iter()
            .map(|vf| output::format_json(vf, false))
            .collect())
    }
}

/// Per-allele scaffold for `--sa-only` mode: creates a TranscriptVariation
/// per alt allele with empty consequences so the SA attachment loop has a
/// slot to populate while emitting no default-CSQ annotation.
pub fn annotate_sa_only_scaffold(vf: &mut VariationFeature) {
    for alt in &vf.alt_alleles {
        vf.transcript_variations.push(TranscriptVariation {
            transcript_id: "-".into(),
            gene_id: "-".into(),
            gene_symbol: None,
            biotype: "-".into(),
            allele_annotations: vec![AlleleAnnotation {
                allele: alt.clone(),
                consequences: vec![],
                impact: fastvep_core::Impact::Modifier,
                cdna_position: None,
                cds_position: None,
                protein_position: None,
                amino_acids: None,
                codons: None,
                exon: None,
                intron: None,
                distance: None,
                protein_length: None,
                escapes_nmd: None,
                hgvsc: None,
                hgvsp: None,
                hgvsg: None,
                hgvs_offset: None,
                existing_variation: vec![],
                sift: None,
                polyphen: None,
                supplementary: Vec::new(),
                acmg_classification: None,
            }],
            canonical: false,
            strand: fastvep_core::Strand::Forward,
            source: None,
            protein_id: None,
            mane_select: None,
            mane_plus_clinical: None,
            tsl: None,
            appris: None,
            ccds: None,
            gencode_primary: false,
            symbol_source: None,
            hgnc_id: None,
            flags: Vec::new(),
        });
    }
}

pub fn annotate_intergenic(vf: &mut VariationFeature) {
    for alt in &vf.alt_alleles {
        vf.transcript_variations.push(TranscriptVariation {
            transcript_id: "-".into(),
            gene_id: "-".into(),
            gene_symbol: None,
            biotype: "-".into(),
            allele_annotations: vec![AlleleAnnotation {
                allele: alt.clone(),
                consequences: vec![Consequence::IntergenicVariant],
                impact: fastvep_core::Impact::Modifier,
                cdna_position: None,
                cds_position: None,
                protein_position: None,
                amino_acids: None,
                codons: None,
                exon: None,
                intron: None,
                distance: None,
                protein_length: None,
                escapes_nmd: None,
                hgvsc: None,
                hgvsp: None,
                hgvsg: None,
                hgvs_offset: None,
                existing_variation: vec![],
                sift: None,
                polyphen: None,
                supplementary: Vec::new(),
                acmg_classification: None,
            }],
            canonical: false,
            strand: fastvep_core::Strand::Forward,
            source: None,
            protein_id: None,
            mane_select: None,
            mane_plus_clinical: None,
            tsl: None,
            appris: None,
            ccds: None,
            gencode_primary: false,
            symbol_source: None,
            hgnc_id: None,
            flags: Vec::new(),
        });
    }
    vf.most_severe_consequence = Some(Consequence::IntergenicVariant);
}

pub fn zip_positions(start: Option<u64>, end: Option<u64>) -> Option<(u64, u64)> {
    match (start, end) {
        (Some(s), Some(e)) => Some((s.min(e), s.max(e))),
        (Some(s), None) => Some((s, s)),
        (None, Some(e)) => Some((e, e)),
        _ => None,
    }
}

/// The `(cDNA anchor, offset)` pair one genomic position maps to on a
/// transcript, with an offset of 0 for an exonic base.
///
/// An indel that reaches from an exon into an intron has one end of each kind,
/// and a span written from its intronic end alone claims the change happens
/// where its far end does. That is not a cosmetic difference: PVS1 reads the
/// offset off HGVSc to decide whether a splice consequence reached the
/// canonical dinucleotide.
pub fn intronic_or_exonic_cdna(
    transcript: &fastvep_genome::Transcript,
    genomic: u64,
) -> Option<(u64, i64)> {
    transcript
        .genomic_to_intronic_cdna(genomic)
        .or_else(|| transcript.genomic_to_cdna(genomic).map(|c| (c, 0)))
}

/// The CDS and everything downstream of it, indexed by CDS position.
///
/// Byte `n - 1` is CDS position `n`, which is what every coordinate in an
/// `AlleleConsequenceResult` means, and the slice runs to the end of the
/// transcript rather than to the annotated terminator - a frameshift's new stop
/// is usually in what was the 3' UTR, and stopping at the old one would leave it
/// unfindable.
///
/// A CDS annotated as beginning part-way through a codon numbers its first
/// codon's missing bases anyway, so those are padded with `N` here exactly as
/// `Transcript::translateable_seq` pads them. Returns `None` when the
/// transcript's coding start does not fall inside its own spliced sequence,
/// which a truncated GFF3 can produce.
pub fn cds_and_downstream(
    transcript: &fastvep_genome::Transcript,
    spliced: &str,
    coding_start: u64,
) -> Option<Vec<u8>> {
    let start = usize::try_from(coding_start.checked_sub(1)?).ok()?;
    let bases = spliced.as_bytes();
    if start > bases.len() {
        return None;
    }
    let phase = transcript.codon_table_start_phase as usize;
    let mut cds = vec![b'N'; phase];
    cds.extend_from_slice(&bases[start..]);
    Some(cds)
}

/// The genetic code a frameshift's downstream peptide is read with.
///
/// Mitochondrial transcripts use NCBI table 2, where AGA and AGG terminate and
/// TGA is Trp; reading one with the nuclear table moves the new terminator.
pub fn frameshift_codon_table(
    transcript: &fastvep_genome::Transcript,
) -> fastvep_genome::CodonTable {
    if fastvep_genome::is_mitochondrial(&transcript.chromosome) {
        fastvep_genome::mitochondrial_codon_table()
    } else {
        fastvep_genome::CodonTable::standard()
    }
}

/// Supplementary-annotation lookups that failed, counted per source.
///
/// A provider error and a provider miss used to be the same thing at the call
/// site - both per-variant loops wrote `if let Ok(Some(ann))`, so a source that
/// *could not answer* produced a record indistinguishable from one it simply had
/// nothing for. That is the silent wrong answer this codebase is most exposed
/// to. A `.osa2` chunk whose JSON blob exceeds the decompression limit fails
/// every lookup in its megabase of the genome, and the only visible effect was a
/// missing column: on a dense converted SpliceAI database, 295 of 300 variants
/// came back with no score and the run reported success (#101).
///
/// Reporting per variant is not an option - the loop runs millions of times - so
/// the first message from each source is kept and the rest counted, for one
/// summary at the end of the run.
#[derive(Default)]
pub struct SaLookupErrors {
    /// `json_key -> (failures, first message)`.
    inner: std::sync::Mutex<std::collections::HashMap<String, (u64, String)>>,
}

impl SaLookupErrors {
    /// Record one failed lookup. Only ever called off the happy path, so the
    /// mutex is never contended by a successful annotation.
    pub fn record(&self, source: &str, err: &anyhow::Error) {
        let Ok(mut map) = self.inner.lock() else {
            return; // a poisoned counter must not take the run down with it
        };
        // Look the source up borrowed, and build the strings only when it is
        // new. A source that is unreadable rather than merely missing a record
        // fails on *every* variant, so this runs once per variant per broken
        // source across every worker; `entry()` would allocate the key and
        // format the error again each time.
        if let Some(seen) = map.get_mut(source) {
            seen.0 += 1;
            return;
        }
        map.insert(source.to_string(), (1, err.to_string()));
    }

    /// One line per source that failed, sorted so the output is reproducible.
    pub fn report_lines(&self) -> Vec<String> {
        let Ok(map) = self.inner.lock() else {
            return Vec::new();
        };
        let mut rows: Vec<_> = map.iter().collect();
        rows.sort_by(|a, b| a.0.cmp(b.0));
        rows.iter()
            .map(|(source, (count, first))| {
                format!(
                    "warning: {source} could not be read for {count} variant lookup(s); \
                     those variants carry no {source} annotation, which is not the same as \
                     having none. First error: {first}"
                )
            })
            .collect()
    }
}

/// The process-wide collector, shared by both annotation pipelines for the same
/// reason the SA caches are shared: the providers are reached from deep inside
/// two independent per-variant loops.
pub fn sa_lookup_errors() -> &'static SaLookupErrors {
    static ERRORS: std::sync::OnceLock<SaLookupErrors> = std::sync::OnceLock::new();
    ERRORS.get_or_init(SaLookupErrors::default)
}

/// Print the supplementary-annotation failure summary, if there is one.
pub fn report_sa_lookup_errors() {
    for line in sa_lookup_errors().report_lines() {
        eprintln!("{line}");
    }
}

/// Put an allele into transcript orientation for a reverse-strand transcript.
///
/// VCF alleles are given in genomic order, so reading them on the minus strand
/// means reversing them as well as complementing each base. Complementing in
/// place is indistinguishable for a single base, which is why it survived: it
/// showed up only once a multi-base alternate reached HGVS, where a
/// reverse-strand `ACG` was written `delinsTGC` instead of `delinsCGT`.
pub fn reverse_complement_allele(allele: &Allele) -> Allele {
    match allele {
        Allele::Sequence(bases) => {
            Allele::Sequence(fastvep_genome::codon::reverse_complement(bases))
        }
        other => other.clone(),
    }
}

/// Extract trio genotype information from a VariationFeature's VCF sample columns.
///
/// Returns (proband, mother, father) GenotypeInfo tuples.
fn extract_trio_genotypes(
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

    // rest[0] is FORMAT, rest[1..] are sample columns
    if vcf_fields.rest.is_empty() {
        return (None, None, None);
    }

    let format_str = &vcf_fields.rest[0];
    let sample_strs: Vec<&str> = vcf_fields.rest[1..].iter().map(|s| s.as_str()).collect();

    let samples = fastvep_io::sample::parse_samples(format_str, &sample_strs, sample_names);

    let proband_gt = samples
        .iter()
        .find(|s| s.name == trio.proband)
        .map(sample_data_to_genotype_info);

    let mother_gt = trio.mother.as_ref().and_then(|name| {
        samples
            .iter()
            .find(|s| &s.name == name)
            .map(sample_data_to_genotype_info)
    });

    let father_gt = trio.father.as_ref().and_then(|name| {
        samples
            .iter()
            .find(|s| &s.name == name)
            .map(sample_data_to_genotype_info)
    });

    (proband_gt, mother_gt, father_gt)
}

/// Convert a SampleData to GenotypeInfo.
fn sample_data_to_genotype_info(
    sample: &fastvep_io::sample::SampleData,
) -> fastvep_classification::GenotypeInfo {
    let gt = sample.genotype.as_ref();
    let is_het = gt.is_some_and(|g| g.is_het());
    let is_hom_ref = gt.is_some_and(|g| g.is_hom_ref());
    let is_hom_alt = gt.is_some_and(|g| g.is_hom_alt());
    let is_missing = gt.is_none_or(|g| g.is_missing());
    let is_phased = gt.is_some_and(|g| g.phased);

    // Determine which alt allele index is carried
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

/// Compound-het enrichment pass: after all variants are annotated,
/// group by gene and identify companion variant relationships,
/// then re-evaluate PM3/BP2 with companion data.
fn enrich_compound_het(
    variants: &mut [VariationFeature],
    acmg_cfg: &fastvep_classification::AcmgConfig,
    sample_names: &[String],
    functional_evidence: Option<&fastvep_classification::FunctionalEvidenceIndex>,
    repeat_db_loaded: bool,
) {
    use std::collections::HashMap;

    // Collect per-gene variant info: (variant_index, gene_symbol, ClinVar P/LP flags, proband_het, is_phased, hgvsc, allele_indices for phase)
    struct VariantGeneInfo {
        vf_idx: usize,
        tv_idx: usize,
        aa_idx: usize,
        is_clinvar_pathogenic: bool,
        is_clinvar_likely_pathogenic: bool,
        proband_het: bool,
        is_phased: bool,
        /// Proband's allele indices for phase comparison
        proband_alleles: Vec<Option<u32>>,
        hgvsc: Option<String>,
    }

    let mut gene_variants: HashMap<String, Vec<VariantGeneInfo>> = HashMap::new();

    for (vf_idx, vf) in variants.iter().enumerate() {
        let trio_genotypes = extract_trio_genotypes(vf, acmg_cfg, sample_names);
        let proband_gt = &trio_genotypes.0;

        for (tv_idx, tv) in vf.transcript_variations.iter().enumerate() {
            let gene_sym = match tv.gene_symbol.as_deref() {
                Some(g) if !g.is_empty() && g != "-" => g.to_string(),
                _ => continue,
            };

            for (aa_idx, aa) in tv.allele_annotations.iter().enumerate() {
                let is_clinvar_pathogenic = aa
                    .acmg_classification
                    .as_ref()
                    .and_then(|v| v.get("criteria"))
                    .and_then(|c| c.as_array())
                    .is_some_and(|criteria| {
                        // Check if this variant has ClinVar pathogenic data
                        criteria.iter().any(|c| {
                            c.get("code")
                                .and_then(|v| v.as_str())
                                .is_some_and(|code| code == "PP5" || code == "PS4")
                                && c.get("met").and_then(|v| v.as_bool()).unwrap_or(false)
                        })
                    });

                // Classify ClinVar supplementary as Pathogenic / Likely pathogenic
                // separately so PM3 v1.0 can score them at their proper point
                // values. A bare substring match on "pathogenic" matches both
                // "Pathogenic" and "Likely pathogenic" — this would over-score
                // LP companions as P. Strip "Likely pathogenic" first, then
                // see if any "pathogenic" remains: that residual signals true P.
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
                        is_clinvar_pathogenic: is_clinvar_pathogenic || clinvar_p_from_sa,
                        is_clinvar_likely_pathogenic: clinvar_lp_from_sa,
                        proband_het,
                        is_phased,
                        proband_alleles,
                        hgvsc: aa.hgvsc.clone(),
                    });
            }
        }
    }

    // For each gene with multiple het variants, build companion relationships and re-classify
    for gene_infos in gene_variants.values() {
        let het_variants: Vec<&VariantGeneInfo> =
            gene_infos.iter().filter(|v| v.proband_het).collect();
        if het_variants.len() < 2 {
            continue;
        }

        // For each het variant, build companion list from other het variants in the gene
        for info in &het_variants {
            let companions: Vec<fastvep_classification::CompanionVariant> = het_variants
                .iter()
                .filter(|other| {
                    other.vf_idx != info.vf_idx
                        || other.tv_idx != info.tv_idx
                        || other.aa_idx != info.aa_idx
                })
                .map(|other| {
                    // Determine trans/cis from phase information
                    let is_in_trans = if info.is_phased && other.is_phased {
                        // Both phased: check if they're on different haplotypes
                        // In a phased genotype like 0|1 vs 1|0, alleles at same index
                        // come from the same parent. So het 0|1 and 1|0 means they're
                        // on different haplotypes (trans).
                        if info.proband_alleles.len() >= 2 && other.proband_alleles.len() >= 2 {
                            let info_alt_on_first = info
                                .proband_alleles
                                .first()
                                .is_some_and(|a| a.is_some_and(|v| v > 0));
                            let other_alt_on_first = other
                                .proband_alleles
                                .first()
                                .is_some_and(|a| a.is_some_and(|v| v > 0));
                            // If alt alleles are on different haplotypes, they're in trans
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

            // Re-extract classification input with companion data and re-classify
            let vf = &variants[info.vf_idx];
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

            let trio_genotypes = extract_trio_genotypes(vf, acmg_cfg, sample_names);

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
                splice_ps1_evidence(aa, &gene_anns, &query_alleles, alt_idx),
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

/// Load supplementary annotation providers (.osa, .osa2, .osi files) from a
/// directory.
///
/// Opening is done in parallel and the results are ordered by path (issue
/// #78). A per-chromosome deployment puts 200+ shards in one `--sa-dir`, and
/// opening a shard is latency-bound rather than CPU-bound, so serial opens made
/// startup scale linearly with shard count for no reason. Sorting by path also
/// makes the provider order - and therefore the output column order - depend on
/// the file names rather than on directory iteration order.
pub fn load_sa_providers(sa_dir: &Path) -> Result<Vec<Box<dyn AnnotationProvider>>> {
    use fastvep_sa::interval::OsiReader;
    use fastvep_sa::reader::SaReader;
    use fastvep_sa::reader_v2::Osa2Reader;
    use rayon::prelude::*;

    if !sa_dir.is_dir() {
        tracing::warn!(
            "SA directory does not exist: {} (skipping)",
            sa_dir.display()
        );
        return Ok(Vec::new());
    }

    let paths = sorted_paths_with_extensions(sa_dir, &["osa2", "osa", "osi"])?;

    // A source that fails to open is skipped with a warning, as before; only
    // the directory walk itself is fatal.
    let providers: Vec<Box<dyn AnnotationProvider>> = paths
        .par_iter()
        .filter_map(|path| {
            let ext = path.extension().and_then(|e| e.to_str());
            let opened: Result<(&str, Box<dyn AnnotationProvider>)> = match ext {
                Some("osa2") => Osa2Reader::open(path).map(|r| ("SA v2", boxed(r))),
                Some("osa") => SaReader::open(path).map(|r| ("SA", boxed(r))),
                // Interval-level (.osi) - typically BED-derived custom sources.
                // Wired up alongside .osa so a directory with mixed file types
                // "just works" via --sa-dir; the OsiReader exposes the same
                // AnnotationProvider trait, returning AnnotationValue::Interval.
                Some("osi") => OsiReader::open(path).map(|r| ("SA interval", boxed(r))),
                _ => return None,
            };
            match opened {
                Ok((kind, provider)) => {
                    tracing::info!("Loaded {}: {} ({})", kind, provider.name(), path.display());
                    Some(provider)
                }
                Err(e) => {
                    tracing::warn!("Could not load {}: {}", path.display(), e);
                    None
                }
            }
        })
        .collect();

    Ok(providers)
}

fn boxed<P: AnnotationProvider + 'static>(provider: P) -> Box<dyn AnnotationProvider> {
    Box::new(provider)
}

/// Every file in `dir` whose extension is one of `exts`, sorted by path so the
/// caller's provider order is reproducible across machines and filesystems.
///
/// A failed directory entry is propagated rather than skipped: on a network
/// `--sa-dir` a transient error would otherwise drop a shard, and the run would
/// finish successfully with that source's annotations simply missing.
fn sorted_paths_with_extensions(dir: &Path, exts: &[&str]) -> Result<Vec<std::path::PathBuf>> {
    let mut paths: Vec<std::path::PathBuf> = Vec::new();
    for entry in std::fs::read_dir(dir)
        .with_context(|| format!("Listing annotation directory {}", dir.display()))?
    {
        let path = entry
            .with_context(|| format!("Reading an entry of {}", dir.display()))?
            .path();
        if path
            .extension()
            .and_then(|e| e.to_str())
            .is_some_and(|e| exts.contains(&e))
        {
            paths.push(path);
        }
    }
    paths.sort();
    Ok(paths)
}

/// Load gene-level annotation providers (.oga files) from a directory.
///
/// Same ordering and parallelism contract as [`load_sa_providers`].
pub fn load_gene_providers(sa_dir: &Path) -> Result<Vec<fastvep_sa::gene::GeneIndex>> {
    use rayon::prelude::*;

    if !sa_dir.is_dir() {
        return Ok(Vec::new());
    }

    let paths = sorted_paths_with_extensions(sa_dir, &["oga"])?;

    let providers: Vec<fastvep_sa::gene::GeneIndex> = paths
        .par_iter()
        .filter_map(|path| {
            match std::fs::File::open(path)
                .map_err(anyhow::Error::from)
                .and_then(|mut f| fastvep_sa::gene::GeneIndex::read_from(&mut f))
            {
                Ok(index) => {
                    tracing::info!(
                        "Loaded gene annotations: {} ({}, {} genes)",
                        index.header.name,
                        path.display(),
                        index.gene_count()
                    );
                    Some(index)
                }
                Err(e) => {
                    tracing::warn!("Could not load {}: {}", path.display(), e);
                    None
                }
            }
        })
        .collect();

    Ok(providers)
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
/// Whether PS1's splice path has a comparison variant for this allele.
///
/// A thin adapter: it pairs the allele with its VCF-form coordinates from
/// [`VariationFeature::query_alleles`] and hands them to the classifier, which
/// owns the rule. `None` whenever the allele has no coordinate to look up, which
/// is the same "cannot tell" the classifier returns for an unloaded index.
///
/// Public, and deliberately: `fastvep-cli` runs the same classification over
/// the same annotations by a different path, and a second copy of this adapter
/// there is how one of the two paths ends up quietly not applying a criterion.
pub fn splice_ps1_evidence(
    aa: &AlleleAnnotation,
    gene_anns: &[&fastvep_core::GeneAnnotation],
    query_alleles: &[(String, u64, String, String)],
    alt_idx: Option<usize>,
) -> Option<bool> {
    let (_, pos, ref_allele, alt) = query_alleles.get(alt_idx?)?;
    fastvep_classification::same_splice_position_pathogenic(
        &aa.consequences,
        gene_anns,
        aa.hgvsc.as_deref(),
        *pos,
        ref_allele,
        alt,
    )
}

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

#[cfg(test)]
mod tests {
    use super::*;

    /// fastvep-annotate is the shared engine both fastvep-cli and
    /// fastvep-web sit on top of, and had zero unit tests before this. This
    /// exercises the full annotate_vcf_text path without needing a real
    /// genome: an empty transcript set (gff3 = None) makes every variant
    /// intergenic, which is enough to verify the pipeline runs end-to-end
    /// and produces well-formed output.
    fn empty_context() -> AnnotationContext {
        AnnotationContext::new(None, None, None, 0).expect("empty context should build")
    }

    #[test]
    fn annotate_vcf_text_returns_one_result_per_variant() {
        let ctx = empty_context();
        let vcf = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
                   1\t100\t.\tA\tG\t.\tPASS\t.\n\
                   1\t200\t.\tC\tT\t.\tPASS\t.\n";
        let results = ctx.annotate_vcf_text(vcf, false).unwrap();
        assert_eq!(results.len(), 2);
        assert_eq!(
            results[0]["most_severe_consequence"],
            serde_json::json!("intergenic_variant")
        );
    }

    #[test]
    fn annotate_vcf_text_with_acmg_none_disables_classification_regardless_of_self_config() {
        let mut ctx = empty_context();
        // self.acmg_config is enabled...
        ctx.acmg_config = Some(fastvep_classification::AcmgConfig::default());
        let vcf = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
                   1\t100\t.\tA\tG\t.\tPASS\t.\n";

        // ...but an explicit None override (as fastvep-web now passes for
        // acmg_requested=false) must take precedence over self.acmg_config,
        // since concurrent requests share one AnnotationContext and must not
        // leak each other's ACMG preference.
        let results = ctx.annotate_vcf_text_with_acmg(vcf, false, None).unwrap();
        assert_eq!(results.len(), 1);
    }

    #[test]
    fn hgvsp_frameshift_skips_gracefully_on_inconsistent_spliced_seq_len() {
        // Regression: the HGVSp-frameshift branch slices
        // `spliced.as_bytes()[coding_start_idx..]` where `coding_start_idx`
        // comes from `tr.cdna_coding_start`. If `spliced_seq` is ever shorter
        // than `coding_start_idx` implies -- e.g. malformed/truncated GFF3-
        // or cache-derived transcript data -- this must not panic ("range
        // start index out of range"); it should just skip HGVSp generation.
        //
        // `cdna_coding_start`/`cdna_coding_end` (used by the predictor, via
        // `cdna_to_cds`, to classify the variant and via `translateable_seq`
        // to compute amino acids) are left internally consistent, so the
        // variant is still correctly classified as a frameshift with real
        // amino acids computed -- only `spliced_seq` itself (used solely by
        // fastvep-annotate's separate HGVSp-slicing code, not by the
        // predictor) is corrupted, isolating the guard under test.
        use fastvep_genome::{Exon, Gene, Transcript, Translation};

        // Two-exon transcript on chr1: exon1 is a 50bp 5' UTR (genomic
        // 1-50), exon2 is the fully-coding CDS (genomic 51-62): ATG AAA CCC
        // TAA (Met Lys Pro Stop).
        let mut transcript = Transcript {
            stable_id: "ENST_FS_TEST".into(),
            version: None,
            gene: Gene {
                stable_id: "ENSG_FS_TEST".into(),
                symbol: Some("FS-TEST".into()),
                symbol_source: None,
                hgnc_id: None,
                biotype: "protein_coding".into(),
                chromosome: "1".into(),
                start: 1,
                end: 62,
                strand: fastvep_core::Strand::Forward,
            },
            biotype: "protein_coding".into(),
            chromosome: "1".into(),
            start: 1,
            end: 62,
            strand: fastvep_core::Strand::Forward,
            exons: vec![
                Exon {
                    stable_id: "ENSE_FS_TEST_1".into(),
                    start: 1,
                    end: 50,
                    strand: fastvep_core::Strand::Forward,
                    phase: -1,
                    end_phase: 0,
                    rank: 1,
                },
                Exon {
                    stable_id: "ENSE_FS_TEST_2".into(),
                    start: 51,
                    end: 62,
                    strand: fastvep_core::Strand::Forward,
                    phase: 0,
                    end_phase: -1,
                    rank: 2,
                },
            ],
            translation: Some(Translation {
                stable_id: "ENSP_FS_TEST".into(),
                genomic_start: 51,
                genomic_end: 62,
                start_exon_rank: 2,
                start_exon_offset: 0,
                end_exon_rank: 2,
                end_exon_offset: 11,
            }),
            cdna_coding_start: Some(51),
            cdna_coding_end: Some(62),
            coding_region_start: Some(51),
            coding_region_end: Some(62),
            spliced_seq: None,
            translateable_seq: None,
            peptide: None,
            canonical: true,
            mane_select: None,
            mane_plus_clinical: None,
            tsl: None,
            appris: None,
            ccds: None,
            protein_id: Some("ENSP_FS_TEST".into()),
            protein_version: None,
            swissprot: vec![],
            trembl: vec![],
            uniparc: vec![],
            refseq_id: None,
            source: None,
            gencode_primary: false,
            flags: vec![],
            codon_table_start_phase: 0,
        };

        transcript
            .build_sequences(|_chrom, start, _end| {
                if start == 1 {
                    Ok(b"N".repeat(50))
                } else {
                    Ok(b"ATGAAACCCTAA".to_vec())
                }
            })
            .expect("build_sequences should succeed for a well-formed test transcript");
        assert_eq!(
            transcript.translateable_seq.as_deref(),
            Some("ATGAAACCCTAA")
        );
        assert_eq!(transcript.spliced_seq.as_deref().map(|s| s.len()), Some(62));

        // Simulate `spliced_seq` becoming shorter than `cdna_coding_start`
        // implies (e.g. a corrupted/truncated cache reload of just this
        // field) -- `cdna_coding_start` (51) is left untouched, so
        // coding_start_idx (50) now exceeds the truncated spliced_seq's
        // length (20).
        transcript.spliced_seq = transcript.spliced_seq.map(|s| s[..20].to_string());

        let mut ctx = empty_context();
        ctx.transcript_provider = IndexedTranscriptProvider::new(vec![transcript]);

        // 1bp deletion of the "G" at genomic position 53 (VCF-anchored at
        // 52): ATG -> AT_ + AAA... = frameshift.
        let vcf = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
                   1\t52\t.\tTG\tT\t.\tPASS\t.\n";

        // Must not panic despite the truncated spliced_seq.
        let results = ctx
            .annotate_vcf_text_with_acmg(vcf, false, None)
            .expect("annotation should succeed even with a truncated spliced_seq");
        assert_eq!(results.len(), 1);

        // Confirm this actually exercised the guarded path: the variant must
        // still be classified as a frameshift with real amino acids (proving
        // the truncated spliced_seq didn't derail classification, since that
        // uses translateable_seq/coding_region_start/coding_region_end
        // instead), but with no hgvsp emitted (proving the slicing guard
        // skipped it rather than panicking).
        let tc = &results[0]["transcript_consequences"][0];
        let terms = tc["consequence_terms"]
            .as_array()
            .expect("consequence_terms should be present");
        assert!(
            terms
                .iter()
                .any(|t| t.as_str() == Some("frameshift_variant")),
            "expected frameshift_variant, got: {:?}",
            terms
        );
        assert!(
            tc.get("amino_acids").is_some(),
            "amino acids should still be computed from the intact translateable_seq: {:?}",
            tc
        );
        assert!(
            tc.get("hgvsp").is_none(),
            "hgvsp should be omitted (skipped), not present, when spliced_seq is \
             shorter than coding_start_idx implies: {:?}",
            tc
        );
    }

    /// Single-CDS-exon transcript behind a 50 bp 5' UTR, built from a CDS given
    /// in transcript orientation. On the reverse strand the genomic sequence is
    /// its reverse complement and CDS base i sits at genomic
    /// `cds_end + 1 - i`, so the same peptide can be exercised from both
    /// directions with one description.
    fn cds_transcript(strand: fastvep_core::Strand, cds: &str) -> fastvep_genome::Transcript {
        use fastvep_genome::{Exon, Gene, Transcript, Translation};
        let cds_len = cds.len() as u64;
        let cds_start = 51;
        let cds_end = cds_start + cds_len - 1;
        let genomic: String = if strand == fastvep_core::Strand::Forward {
            cds.to_string()
        } else {
            cds.chars()
                .rev()
                .map(|c| match c {
                    'A' => 'T',
                    'T' => 'A',
                    'C' => 'G',
                    'G' => 'C',
                    other => other,
                })
                .collect()
        };

        let mut transcript = Transcript {
            stable_id: "ENST_ANCHOR".into(),
            version: None,
            gene: Gene {
                stable_id: "ENSG_ANCHOR".into(),
                symbol: Some("ANCHOR-TEST".into()),
                symbol_source: None,
                hgnc_id: None,
                biotype: "protein_coding".into(),
                chromosome: "1".into(),
                start: 1,
                end: cds_end,
                strand,
            },
            biotype: "protein_coding".into(),
            chromosome: "1".into(),
            start: 1,
            end: cds_end,
            strand,
            exons: vec![
                Exon {
                    stable_id: "ENSE_ANCHOR_1".into(),
                    start: 1,
                    end: 50,
                    strand,
                    phase: -1,
                    end_phase: 0,
                    rank: if strand == fastvep_core::Strand::Forward {
                        1
                    } else {
                        2
                    },
                },
                Exon {
                    stable_id: "ENSE_ANCHOR_2".into(),
                    start: cds_start,
                    end: cds_end,
                    strand,
                    phase: 0,
                    end_phase: -1,
                    rank: if strand == fastvep_core::Strand::Forward {
                        2
                    } else {
                        1
                    },
                },
            ],
            translation: Some(Translation {
                stable_id: "ENSP_ANCHOR".into(),
                genomic_start: cds_start,
                genomic_end: cds_end,
                start_exon_rank: if strand == fastvep_core::Strand::Forward {
                    2
                } else {
                    1
                },
                start_exon_offset: 0,
                end_exon_rank: if strand == fastvep_core::Strand::Forward {
                    2
                } else {
                    1
                },
                end_exon_offset: cds_len - 1,
            }),
            // cDNA runs 5'->3' along the transcript, so the CDS follows the
            // 50 bp UTR only on the forward strand. On the reverse strand the
            // genomically-rightmost exon comes first, putting the CDS at the
            // very start of the cDNA.
            cdna_coding_start: Some(if strand == fastvep_core::Strand::Forward {
                cds_start
            } else {
                1
            }),
            cdna_coding_end: Some(if strand == fastvep_core::Strand::Forward {
                cds_end
            } else {
                cds_len
            }),
            coding_region_start: Some(cds_start),
            coding_region_end: Some(cds_end),
            spliced_seq: None,
            translateable_seq: None,
            peptide: None,
            canonical: true,
            mane_select: None,
            mane_plus_clinical: None,
            tsl: None,
            appris: None,
            ccds: None,
            protein_id: Some("ENSP_ANCHOR".into()),
            protein_version: None,
            swissprot: vec![],
            trembl: vec![],
            uniparc: vec![],
            refseq_id: None,
            source: None,
            gencode_primary: false,
            flags: vec![],
            codon_table_start_phase: 0,
        };

        let genomic_bytes = genomic.into_bytes();
        transcript
            .build_sequences(move |_chrom, start, _end| {
                if start == 1 {
                    Ok(b"N".repeat(50))
                } else {
                    Ok(genomic_bytes.clone())
                }
            })
            .expect("build_sequences should succeed for a well-formed test transcript");
        transcript
    }

    /// Annotate one VCF record against one transcript and return its first
    /// transcript consequence.
    fn annotate_one(
        transcript: fastvep_genome::Transcript,
        pos: u64,
        r: &str,
        a: &str,
    ) -> serde_json::Value {
        let mut ctx = empty_context();
        ctx.transcript_provider = IndexedTranscriptProvider::new(vec![transcript]);
        let vcf = format!(
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
             1\t{pos}\t.\t{r}\t{a}\t.\tPASS\t.\n"
        );
        let results = ctx.annotate_vcf_text_with_acmg(&vcf, false, None).unwrap();
        assert_eq!(results.len(), 1, "expected one annotated record");
        results[0]["transcript_consequences"][0].clone()
    }

    /// Every residue an HGVSp names must actually be that residue in the
    /// peptide. This is the property that a wrong anchor breaks: the letters come
    /// from `Amino_acids` and the numbers from the anchor, so a mismatched pair
    /// yields `p.Ser1092_Phe1094delinsPhe` on a protein whose 1092-1094 is FQP.
    fn assert_hgvsp_agrees_with_peptide(hgvsp: &str, peptide: &str) {
        let three_to_one = |aa3: &str| -> Option<char> {
            const TABLE: &[(&str, char)] = &[
                ("Ala", 'A'),
                ("Arg", 'R'),
                ("Asn", 'N'),
                ("Asp", 'D'),
                ("Cys", 'C'),
                ("Gln", 'Q'),
                ("Glu", 'E'),
                ("Gly", 'G'),
                ("His", 'H'),
                ("Ile", 'I'),
                ("Leu", 'L'),
                ("Lys", 'K'),
                ("Met", 'M'),
                ("Phe", 'F'),
                ("Pro", 'P'),
                ("Ser", 'S'),
                ("Thr", 'T'),
                ("Trp", 'W'),
                ("Tyr", 'Y'),
                ("Val", 'V'),
                ("Sec", 'U'),
                ("Ter", '*'),
            ];
            TABLE.iter().find(|(k, _)| *k == aa3).map(|(_, v)| *v)
        };
        let body = hgvsp.split(":p.").nth(1).unwrap_or(hgvsp);
        // Only the positions left of `delins`/`ins` name reference residues; the
        // replacement residues after it are not in the reference peptide.
        let named = body.split("ins").next().unwrap_or(body);
        let mut checked = 0;
        let bytes = named.as_bytes();
        let mut i = 0;
        while i + 3 <= bytes.len() {
            let Some(one) = three_to_one(&named[i..i + 3]) else {
                i += 1;
                continue;
            };
            let digits: String = named[i + 3..]
                .chars()
                .take_while(|c| c.is_ascii_digit())
                .collect();
            if digits.is_empty() {
                i += 1;
                continue;
            }
            let pos: usize = digits.parse().unwrap();
            let actual = peptide.chars().nth(pos - 1);
            assert_eq!(
                actual,
                Some(one),
                "{hgvsp} names {} at residue {pos}, but the peptide has {:?} there ({peptide})",
                &named[i..i + 3],
                actual
            );
            checked += 1;
            i += 3 + digits.len();
        }
        assert!(checked > 0, "no residue positions found in {hgvsp}");
    }

    #[test]
    fn a_substitution_at_a_selenocysteine_codon_is_a_missense_not_a_stop_change() {
        // Ensembl annotates selenoprotein CDSs straight through their in-frame
        // UGA, which encodes selenocysteine. Translating that codon as a
        // terminator made the protein look like it ended there - SELENOW is 87
        // residues with its UGA at 13 - and made a substitution at the codon
        // look like a change to a stop rather than the missense it is.
        //
        // CDS ATG TGG TGA CGG TAA on the forward strand: M W U R *, with the
        // UGA at residue 3 and the real terminator at 5.
        use fastvep_core::Strand;
        let tr = cds_transcript(Strand::Forward, "ATGTGGTGACGGTAA");
        assert_eq!(
            tr.peptide.as_deref(),
            Some("MWUR*"),
            "the annotated CDS runs past the UGA, so residue 3 is selenocysteine"
        );

        // TGA -> TGC at the third base: Sec becomes Cys.
        let tc = annotate_one(tr, 59, "A", "C");
        let terms: Vec<&str> = tc["consequence_terms"]
            .as_array()
            .unwrap()
            .iter()
            .filter_map(|t| t.as_str())
            .collect();
        assert!(
            terms.contains(&"missense_variant"),
            "changing selenocysteine to cysteine is a missense, got {terms:?}"
        );
        assert!(
            !terms.iter().any(|t| t.contains("stop")),
            "must not be reported as any kind of stop change, got {terms:?}"
        );
        assert_eq!(tc["amino_acids"], serde_json::json!("U/C"));
        assert_eq!(
            tc["hgvsp"],
            serde_json::json!("ENSP_ANCHOR:p.Sec3Cys"),
            "selenocysteine renders as Sec, not as the \"???\" an unknown letter would give: {tc:?}"
        );
    }

    #[test]
    fn a_transcript_whose_cds_does_not_end_in_a_terminator_keeps_its_internal_stop() {
        // The guard that keeps readthrough from rewriting genuinely broken
        // annotations: without a terminator at the annotated end there is no
        // corroboration that the frame is right, so an internal UGA could as
        // easily be a CDS annotated in the wrong frame.
        use fastvep_core::Strand;
        let tr = cds_transcript(Strand::Forward, "ATGTGGTGACGGCGG");
        assert_eq!(
            tr.peptide.as_deref(),
            Some("MW*RR"),
            "no terminator at the end means the internal stop is left alone"
        );
    }

    #[test]
    fn hgvsp_inframe_insertion_is_normalised_on_the_reverse_strand_too() {
        // Asked for in the #86 review and not in the merged branch. The forward
        // case is covered below; a reverse-strand transcript exercises the
        // coordinate flip, where `protein_start` is derived from the genomic left
        // edge and so runs against the direction the transcript reads.
        //
        // CDS ATG TGG CGG TAA sits at genomic 51-62, so CDS base i is at genomic
        // 63 - i. Duplicating Arg3 means inserting after CDS 9 (genomic 54),
        // i.e. between genomic 53 and 54, with the bases complemented.
        use fastvep_core::Strand;
        let tr = cds_transcript(Strand::Reverse, "ATGTGGCGGTAA");
        assert_eq!(
            tr.peptide.as_deref(),
            Some("MWR*"),
            "test transcript premise"
        );

        let tc = annotate_one(tr, 53, "A", "ACCG");
        assert!(
            tc["consequence_terms"]
                .as_array()
                .unwrap()
                .iter()
                .any(|t| t.as_str() == Some("inframe_insertion")),
            "premise: should be an in-frame insertion, got {:?}",
            tc["consequence_terms"]
        );
        // See the note in the two-codon test: a codon-aligned insertion has no
        // reference codon, so the reference side is `-`.
        assert_eq!(tc["amino_acids"], serde_json::json!("-/R"));
        assert_eq!(
            tc["hgvsp"],
            serde_json::json!("ENSP_ANCHOR:p.Arg3dup"),
            "reverse-strand insertion repeating the residue it follows is a dup: {tc:?}"
        );
        assert_hgvsp_agrees_with_peptide(tc["hgvsp"].as_str().unwrap(), "MWR*");
    }

    #[test]
    fn hgvsp_inframe_insertion_of_two_codons_is_a_two_residue_dup() {
        // The second case asked for in the #86 review. A duplication spanning
        // more than one residue has to render as a range, `p.Arg3_Lys4dup`, not
        // as two separate events and not as an unshifted delins.
        //
        // CDS ATG TGG CGG AAG GAA TAA; inserting CGGAAG after CDS 12 repeats
        // Arg3-Lys4. The trailing Glu keeps the insertion off the terminator,
        // which is a separate unnormalised case (see the test below).
        use fastvep_core::Strand;
        let tr = cds_transcript(Strand::Forward, "ATGTGGCGGAAGGAATAA");
        assert_eq!(
            tr.peptide.as_deref(),
            Some("MWRKE*"),
            "test transcript premise"
        );

        let tc = annotate_one(tr, 62, "G", "GCGGAAG");
        // `-/RK`, not `E/RKE`: an insertion that falls exactly on a codon
        // boundary replaces no codon, so Ensembl's window is empty on the
        // reference side and holds only the inserted codons on the alternate.
        // This used to read `E/RKE`, which repeats the flanking residue on both
        // sides of a change that does not touch it.
        assert_eq!(tc["amino_acids"], serde_json::json!("-/RK"));
        assert_eq!(
            tc["hgvsp"],
            serde_json::json!("ENSP_ANCHOR:p.Arg3_Lys4dup"),
            "two-codon duplication should render as a residue range: {tc:?}"
        );
        assert_hgvsp_agrees_with_peptide(tc["hgvsp"].as_str().unwrap(), "MWRKE*");
    }

    #[test]
    fn hgvsp_names_real_residues_for_a_reverse_strand_shrinking_change() {
        // The bug behind #91. For a shrinking in-frame change the predictor's
        // `protein_start` is the genomic left edge, which on the reverse strand
        // is the *end* of the affected residue range. HGVSp took its residue
        // letters from `Amino_acids` and its numbers from that anchor, so it
        // emitted positions the protein does not have: on real ClinVar data,
        // `p.Ser1092_Phe1094delinsPhe` for a protein whose 1092-1094 is FQP,
        // with the SSF at 1090-1092. 17,654 of 80,679 in-frame indel
        // descriptions were affected.
        //
        // CDS ATG TGG TTT TTC AAG TAA = M W F F K * at genomic 51-68. Deleting
        // genomic 58-60 removes CDS 9-11, taking the FF pair to a single F, and
        // the anchor arrives as residue 4 - the end of the 3-4 range.
        use fastvep_core::Strand;
        let tr = cds_transcript(Strand::Reverse, "ATGTGGTTTTTCAAGTAA");
        assert_eq!(
            tr.peptide.as_deref(),
            Some("MWFFK*"),
            "test transcript premise"
        );

        let tc = annotate_one(tr, 57, "GAAA", "G");
        assert!(
            tc["consequence_terms"]
                .as_array()
                .unwrap()
                .iter()
                .any(|t| t.as_str() == Some("inframe_deletion")),
            "premise: should be an in-frame deletion, got {:?}",
            tc["consequence_terms"]
        );
        assert_eq!(tc["amino_acids"], serde_json::json!("FF/F"));

        // Before the fix: `p.Phe4_Phe5delinsPhe`, which calls residue 5 Phe when
        // the peptide has Lys there. The 3'-rule puts the deletion at the
        // C-terminal F of the pair, so residue 4.
        assert_eq!(
            tc["hgvsp"],
            serde_json::json!("ENSP_ANCHOR:p.Phe4del"),
            "shrinking reverse-strand change should name the residues it deletes: {tc:?}"
        );
        assert_hgvsp_agrees_with_peptide(tc["hgvsp"].as_str().unwrap(), "MWFFK*");
    }

    #[test]
    fn hgvsp_reads_a_periodic_reference_from_the_end_the_strand_determines() {
        // The residue #93 left behind, filed as #96. Trying the other end of the
        // span rescued the cases where the first end is not a match at all; where
        // the reference is periodic *both* ends match, and taking whichever came
        // first took the wrong one. The strand says which end applies, so it now
        // decides instead of the iteration order.
        //
        // CDS ATG GAA GGT GAA GGT GAA GCT TAA = M E G E G E A * at genomic 51-74,
        // so CDS base i sits at genomic 75 - i. `EGE` therefore sits at residues
        // 2-4 and again at 4-6. Deleting CDS 4-12 (genomic 63-71) removes
        // residues 2-4 and the anchor arrives as residue 4, the end of the range.
        use fastvep_core::Strand;
        let tr = cds_transcript(Strand::Reverse, "ATGGAAGGTGAAGGTGAAGCTTAA");
        assert_eq!(
            tr.peptide.as_deref(),
            Some("MEGEGEA*"),
            "test transcript premise"
        );

        let tc = annotate_one(tr, 62, "CTTCACCTTC", "C");
        assert!(
            tc["consequence_terms"]
                .as_array()
                .unwrap()
                .iter()
                .any(|t| t.as_str() == Some("inframe_deletion")),
            "premise: should be an in-frame deletion, got {:?}",
            tc["consequence_terms"]
        );
        assert_eq!(tc["amino_acids"], serde_json::json!("EGE/-"));
        assert_eq!(tc["protein_start"], serde_json::json!(2));
        assert_eq!(tc["protein_end"], serde_json::json!(4));

        // Before the fix: `p.Glu4_Glu6del`. Both descriptions are well-formed and
        // name real Glu residues, which is what makes this one dangerous - only
        // one of them describes this variant. Deleting 2-4 leaves MGEA; deleting
        // 4-6 leaves MEGA, so no later 3'-shift can reconcile them.
        assert_eq!(
            tc["hgvsp"],
            serde_json::json!("ENSP_ANCHOR:p.Glu2_Glu4del"),
            "periodic reference resolved to the wrong end of the span: {tc:?}"
        );
        assert_hgvsp_agrees_with_peptide(tc["hgvsp"].as_str().unwrap(), "MEGEGEA*");
    }

    #[test]
    fn hgvsp_reads_the_same_periodic_reference_from_the_start_on_the_forward_strand() {
        // The mirror: same peptide, same deleted residues, transcript running the
        // other way. Here `protein_start` already *is* the first affected residue,
        // so ordering the candidates by strand has to leave this untouched -
        // otherwise the fix for #96 would trade one wrong end for the other.
        //
        // Genomic 54-62 is CDS 4-12 on this strand, anchored at 53.
        use fastvep_core::Strand;
        let tr = cds_transcript(Strand::Forward, "ATGGAAGGTGAAGGTGAAGCTTAA");
        assert_eq!(
            tr.peptide.as_deref(),
            Some("MEGEGEA*"),
            "test transcript premise"
        );

        let tc = annotate_one(tr, 53, "GGAAGGTGAA", "G");
        assert_eq!(tc["amino_acids"], serde_json::json!("EGE/-"));
        assert_eq!(
            tc["hgvsp"],
            serde_json::json!("ENSP_ANCHOR:p.Glu2_Glu4del"),
            "forward-strand reading changed: {tc:?}"
        );
        assert_hgvsp_agrees_with_peptide(tc["hgvsp"].as_str().unwrap(), "MEGEGEA*");
    }

    #[test]
    fn hgvsp_inframe_insertion_is_normalised_not_a_substitution() {
        // Regression for issue #81: the HGVSp routing tested only
        // `aa.1 == "-" || InframeDeletion`, which no in-frame *insertion*
        // satisfies, so insertions fell through to `hgvsp()`. That compares
        // one byte per side, so `R/RR` looked unchanged and produced
        // `p.Arg3=` -- a synonymous call for a variant that lengthens the
        // protein by a residue. Unlike the malformed `p.Xaa123???` that #58
        // fixed, a well-formed `=` or missense is indistinguishable from a
        // genuine call downstream, which is what makes it dangerous.
        //
        // `CGG` inserted after the Trp codon repeats the Arg that follows it,
        // so the HGVS 3'-rule shifts the insertion one residue right and
        // writes the repeat as a duplication. This is the value Ensembl VEP
        // reports, and it exercises both halves of the guard: an insertion
        // never renders as a substitution, and the description it does render
        // is normalised rather than a literal `delinsArgArg`.
        use fastvep_genome::{Exon, Gene, Transcript, Translation};

        // Single-CDS-exon transcript on chr1 behind a 50 bp 5' UTR.
        // CDS (genomic 51-62) is ATG TGG CGG TAA = Met Trp Arg Ter.
        let mut transcript = Transcript {
            stable_id: "ENST_INS_TEST".into(),
            version: None,
            gene: Gene {
                stable_id: "ENSG_INS_TEST".into(),
                symbol: Some("INS-TEST".into()),
                symbol_source: None,
                hgnc_id: None,
                biotype: "protein_coding".into(),
                chromosome: "1".into(),
                start: 1,
                end: 62,
                strand: fastvep_core::Strand::Forward,
            },
            biotype: "protein_coding".into(),
            chromosome: "1".into(),
            start: 1,
            end: 62,
            strand: fastvep_core::Strand::Forward,
            exons: vec![
                Exon {
                    stable_id: "ENSE_INS_TEST_1".into(),
                    start: 1,
                    end: 50,
                    strand: fastvep_core::Strand::Forward,
                    phase: -1,
                    end_phase: 0,
                    rank: 1,
                },
                Exon {
                    stable_id: "ENSE_INS_TEST_2".into(),
                    start: 51,
                    end: 62,
                    strand: fastvep_core::Strand::Forward,
                    phase: 0,
                    end_phase: -1,
                    rank: 2,
                },
            ],
            translation: Some(Translation {
                stable_id: "ENSP_INS_TEST".into(),
                genomic_start: 51,
                genomic_end: 62,
                start_exon_rank: 2,
                start_exon_offset: 0,
                end_exon_rank: 2,
                end_exon_offset: 11,
            }),
            cdna_coding_start: Some(51),
            cdna_coding_end: Some(62),
            coding_region_start: Some(51),
            coding_region_end: Some(62),
            spliced_seq: None,
            translateable_seq: None,
            peptide: None,
            canonical: true,
            mane_select: None,
            mane_plus_clinical: None,
            tsl: None,
            appris: None,
            ccds: None,
            protein_id: Some("ENSP_INS_TEST".into()),
            protein_version: None,
            swissprot: vec![],
            trembl: vec![],
            uniparc: vec![],
            refseq_id: None,
            source: None,
            gencode_primary: false,
            flags: vec![],
            codon_table_start_phase: 0,
        };

        transcript
            .build_sequences(|_chrom, start, _end| {
                if start == 1 {
                    Ok(b"N".repeat(50))
                } else {
                    Ok(b"ATGTGGCGGTAA".to_vec())
                }
            })
            .expect("build_sequences should succeed for a well-formed test transcript");

        let mut ctx = empty_context();
        ctx.transcript_provider = IndexedTranscriptProvider::new(vec![transcript]);

        // Insert CGG after the Trp codon (genomic 54-56), anchored at 56:
        // Met Trp Arg Ter -> Met Trp Arg Arg Ter, an in-frame duplication
        // of Arg3.
        let vcf = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
                   1\t56\t.\tG\tGCGG\t.\tPASS\t.\n";

        let results = ctx.annotate_vcf_text_with_acmg(vcf, false, None).unwrap();
        assert_eq!(results.len(), 1);
        let tc = &results[0]["transcript_consequences"][0];

        // Guard the premise: if this stops being classified as an in-frame
        // insertion the assertion below would pass for the wrong reason.
        let terms = tc["consequence_terms"]
            .as_array()
            .expect("consequence_terms should be present");
        assert!(
            terms
                .iter()
                .any(|t| t.as_str() == Some("inframe_insertion")),
            "expected inframe_insertion, got: {:?}",
            terms
        );
        // See the note in `hgvsp_inframe_insertion_of_two_codons_is_a_two_residue_dup`.
        assert_eq!(tc["amino_acids"], serde_json::json!("-/R"));

        assert_eq!(
            tc["hgvsp"],
            serde_json::json!("ENSP_INS_TEST:p.Arg3dup"),
            "in-frame insertion must render as a normalised duplication, not a \
             substitution and not an unshifted delins: {:?}",
            tc
        );
    }

    /// A minus-strand transcript reads a VCF allele reverse-complemented, not
    /// complemented base-by-base. The two agree for a single base, which is why
    /// complementing in place survived until a multi-base alternate reached
    /// HGVS: `ACG` was written `TGC` instead of `CGT`.
    #[test]
    fn an_allele_in_transcript_orientation_is_reversed_as_well_as_complemented() {
        assert_eq!(
            reverse_complement_allele(&Allele::Sequence(b"ACG".to_vec())),
            Allele::Sequence(b"CGT".to_vec())
        );
        // A palindrome would not have caught the bug, and neither would an SNV.
        assert_eq!(
            reverse_complement_allele(&Allele::Sequence(b"A".to_vec())),
            Allele::Sequence(b"T".to_vec())
        );
        // Applying it twice is the identity, on any sequence.
        let original = Allele::Sequence(b"ACGGTTA".to_vec());
        assert_eq!(
            reverse_complement_allele(&reverse_complement_allele(&original)),
            original
        );
        // The placeholder alleles carry no sequence to orient.
        assert_eq!(
            reverse_complement_allele(&Allele::Deletion),
            Allele::Deletion
        );
        assert_eq!(reverse_complement_allele(&Allele::Missing), Allele::Missing);
    }
}
