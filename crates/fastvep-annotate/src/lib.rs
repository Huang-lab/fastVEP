//! Shared annotation engine for fastVEP.
//!
//! Provides [`AnnotationContext`] which loads transcript models, reference sequences,
//! and supplementary annotation sources, then annotates VCF variants against them.
//!
//! Used by both `fastvep-web` (production HTTP server) and `fastvep-cli` (embedded web server).
//! The CLI batch pipeline (`run_annotate`) has its own streaming implementation but shares
//! the same underlying crates.

mod hgvs_normalize;

pub use hgvs_normalize::{
    convert_ins_to_dup, convert_ins_to_dup_noncoding, convert_ins_to_dup_range,
    convert_ins_to_dup_range_noncoding, three_prime_shift_intronic,
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
use fastvep_io::output;
use fastvep_io::variant::{AlleleAnnotation, TranscriptVariation, VariationFeature};
use fastvep_io::vcf::VcfParser;
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

        // Build sequences for coding transcripts (parallel via rayon)
        if let Some(ref sp) = seq_provider {
            let built = AtomicUsize::new(0);
            transcripts.par_iter_mut().for_each(|tr| {
                if tr.is_coding() && tr.spliced_seq.is_none()
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
            tracing::info!("Built sequences for {} coding transcripts", built);

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
                if tr.is_coding() && tr.spliced_seq.is_none()
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
                tracing::info!("Built sequences for {} coding transcripts", built);
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
    pub fn annotate_vcf_text(
        &self,
        vcf_text: &str,
        pick: bool,
    ) -> Result<Vec<serde_json::Value>> {
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
            && self
                .sa_providers
                .iter()
                .any(|sa| sa.json_key() == "gnomad");

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
                let ref_seq = self.seq_provider.as_ref().and_then(|sp| {
                    sp.fetch_sequence(chrom, query_start, query_end).ok()
                });

                let result = self.predictor.predict(
                    &vf.position,
                    &vf.ref_allele,
                    &vf.alt_alleles,
                    &overlapping,
                    ref_seq.as_deref(),
                );

                for tc in &result.transcript_consequences {
                    let transcript =
                        overlapping.iter().find(|t| t.stable_id == tc.transcript_id);

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
                                protein_position: zip_positions(
                                    ac.protein_start,
                                    ac.protein_end,
                                ),
                                amino_acids: ac.amino_acids.clone(),
                                codons: ac.codons.clone(),
                                exon: ac.exon,
                                intron: ac.intron,
                                distance: ac.distance,
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
                                                complement_allele(&vf.ref_allele),
                                                complement_allele(&ac.allele),
                                            )
                                        } else {
                                            (vf.ref_allele.clone(), ac.allele.clone())
                                        };
                                    if let Some(coding_start) = tr.cdna_coding_start {
                                        if let (Some(cs), Some(ce)) =
                                            (ac.cdna_start, ac.cdna_end)
                                        {
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
                                        } else if ac.intron.is_some() {
                                            if let Some((cdna_pos, offset)) =
                                                tr.genomic_to_intronic_cdna(vf.position.start)
                                            {
                                                let (end_cdna, end_offset) =
                                                    if vf.position.start != vf.position.end {
                                                        tr.genomic_to_intronic_cdna(
                                                            vf.position.end,
                                                        )
                                                        .map(|(c, o)| (Some(c), Some(o)))
                                                        .unwrap_or((None, None))
                                                    } else {
                                                        (None, None)
                                                    };
                                                ann.hgvsc = fastvep_hgvs::hgvsc_intronic_range(
                                                    &versioned_tid,
                                                    cdna_pos,
                                                    offset,
                                                    end_cdna,
                                                    end_offset,
                                                    &hgvs_ref,
                                                    &hgvs_alt,
                                                    coding_start,
                                                    tr.cdna_coding_end,
                                                );
                                            }
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
                                        );
                                    } else if ac.intron.is_some() {
                                        if let Some((cdna_pos, offset)) =
                                            tr.genomic_to_intronic_cdna(vf.position.start)
                                        {
                                            let (end_cdna, end_offset) =
                                                if vf.position.start != vf.position.end {
                                                    tr.genomic_to_intronic_cdna(vf.position.end)
                                                        .map(|(c, o)| (Some(c), Some(o)))
                                                        .unwrap_or((None, None))
                                                } else {
                                                    (None, None)
                                                };
                                            ann.hgvsc =
                                                fastvep_hgvs::hgvsc_noncoding_intronic_range(
                                                    &versioned_tid,
                                                    cdna_pos,
                                                    offset,
                                                    end_cdna,
                                                    end_offset,
                                                    &hgvs_ref,
                                                    &hgvs_alt,
                                                );
                                        }
                                    }

                                    // HGVSp
                                    if let (Some(ref aa), Some(ps)) =
                                        (&ac.amino_acids, ac.protein_start)
                                    {
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
                                                if let (
                                                    Some(ref spliced),
                                                    Some(coding_start),
                                                    Some(cds_s),
                                                ) = (
                                                    &tr.spliced_seq,
                                                    tr.cdna_coding_start,
                                                    ac.cds_start,
                                                ) {
                                                    // Guard against malformed/truncated
                                                    // GFF3-derived transcript data where
                                                    // `coding_start` is inconsistent with the
                                                    // actual spliced sequence length — skip
                                                    // HGVSp generation for this case rather
                                                    // than panicking on an out-of-bounds slice.
                                                    let coding_start_idx =
                                                        (coding_start - 1) as usize;
                                                    let spliced_bytes: &[u8] =
                                                        spliced.as_bytes();
                                                    if coding_start >= 1
                                                        && coding_start_idx <= spliced_bytes.len()
                                                    {
                                                    let ref_from_cds =
                                                        &spliced_bytes[coding_start_idx..];
                                                    let cds_idx = (cds_s - 1) as usize;
                                                    let mut alt_from_cds =
                                                        ref_from_cds.to_vec();
                                                    if ac.allele == Allele::Deletion {
                                                        let del_len = vf.ref_allele.len();
                                                        let end = (cds_idx + del_len)
                                                            .min(alt_from_cds.len());
                                                        alt_from_cds.drain(cds_idx..end);
                                                    } else if let Allele::Sequence(ins_bases) =
                                                        &ac.allele
                                                    {
                                                        let mut bases = ins_bases.clone();
                                                        if tr.strand
                                                            == fastvep_core::Strand::Reverse
                                                        {
                                                            bases = bases
                                                                .iter()
                                                                .map(|&b| match b {
                                                                    b'A' => b'T',
                                                                    b'T' => b'A',
                                                                    b'C' => b'G',
                                                                    b'G' => b'C',
                                                                    o => o,
                                                                })
                                                                .collect();
                                                        }
                                                        for (j, &b) in
                                                            bases.iter().enumerate()
                                                        {
                                                            if cds_idx + j
                                                                <= alt_from_cds.len()
                                                            {
                                                                alt_from_cds
                                                                    .insert(cds_idx + j, b);
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
                                                    ann.hgvsp =
                                                        fastvep_hgvs::hgvsp_frameshift(
                                                            &versioned_pid,
                                                            ref_from_cds,
                                                            &alt_from_cds,
                                                            codon_start,
                                                            &fs_codon_table,
                                                        );
                                                    }
                                                }
                                            } else if aa.1 == "-"
                                                || ac
                                                    .consequences
                                                    .contains(&Consequence::InframeDeletion)
                                            {
                                                // In-frame deletion / delins (frameshift
                                                // handled above). aa.0 holds the deleted
                                                // residues, aa.1 the replacement ("-" for a
                                                // pure deletion).
                                                ann.hgvsp = fastvep_hgvs::hgvsp_inframe_deletion(
                                                    &versioned_pid,
                                                    ps,
                                                    &aa.0,
                                                    &aa.1,
                                                );
                                            } else {
                                                let ref_aa_byte = aa
                                                    .0
                                                    .as_bytes()
                                                    .first()
                                                    .copied()
                                                    .unwrap_or(b'X');
                                                let alt_aa_byte = aa
                                                    .1
                                                    .as_bytes()
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

                    let should_include =
                        !pick || tc.canonical || vf.transcript_variations.is_empty();
                    if should_include {
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
                            mane_plus_clinical: transcript
                                .and_then(|t| t.mane_plus_clinical.clone()),
                            tsl: transcript.and_then(|t| t.tsl),
                            appris: transcript.and_then(|t| t.appris.clone()),
                            ccds: transcript.and_then(|t| t.ccds.clone()),
                            gencode_primary: transcript
                                .map(|t| t.gencode_primary)
                                .unwrap_or(false),
                            symbol_source: transcript
                                .and_then(|t| t.gene.symbol_source.clone()),
                            hgnc_id: transcript.and_then(|t| t.gene.hgnc_id.clone()),
                            flags: transcript.map(|t| t.flags.clone()).unwrap_or_default(),
                        });
                    }
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
                let mut allele_results: std::collections::HashMap<
                    String,
                    Vec<(String, String)>,
                > = std::collections::HashMap::new();
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
                            if let Ok(Some(ann)) =
                                sa.annotate_position(chrom, q_pos, q_ref, q_alt)
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
            if let Some(acmg_cfg) = acmg_config {
                // Parse sample genotypes if trio config is present
                let trio_genotypes = extract_trio_genotypes(vf, acmg_cfg, &sample_names);

                let functional_by_alt = resolve_functional_by_alt(functional_evidence, vf);

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
                                fastvep_classification::is_pure_insertion(&vf.ref_allele),
                                vf.alt_alleles
                                    .iter()
                                    .position(|a| *a == aa.allele)
                                    .and_then(|i| functional_by_alt[i].clone()),
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
                );
            }
        }

        Ok(variants.iter().map(|vf| output::format_json(vf, false)).collect())
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

pub fn complement_allele(allele: &Allele) -> Allele {
    match allele {
        Allele::Sequence(bases) => {
            let comp: Vec<u8> = bases
                .iter()
                .map(|&b| match b {
                    b'A' | b'a' => b'T',
                    b'T' | b't' => b'A',
                    b'C' | b'c' => b'G',
                    b'G' | b'g' => b'C',
                    other => other,
                })
                .collect();
            Allele::Sequence(comp)
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

    let samples =
        fastvep_io::sample::parse_samples(format_str, &sample_strs, sample_names);

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
    let alt_allele_index = gt.and_then(|g| {
        g.alleles
            .iter()
            .filter_map(|a| *a)
            .find(|&a| a > 0)
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

/// Compound-het enrichment pass: after all variants are annotated,
/// group by gene and identify companion variant relationships,
/// then re-evaluate PM3/BP2 with companion data.
fn enrich_compound_het(
    variants: &mut [VariationFeature],
    acmg_cfg: &fastvep_classification::AcmgConfig,
    sample_names: &[String],
    functional_evidence: Option<&fastvep_classification::FunctionalEvidenceIndex>,
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
                                && c.get("met")
                                    .and_then(|v| v.as_bool())
                                    .unwrap_or(false)
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
                fastvep_classification::is_pure_insertion(&vf.ref_allele),
                vf.alt_alleles
                    .iter()
                    .position(|a| *a == aa.allele)
                    .and_then(|i| functional_by_alt[i].clone()),
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
pub fn load_sa_providers(
    sa_dir: &Path,
) -> Result<Vec<Box<dyn AnnotationProvider>>> {
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
pub fn load_gene_providers(
    sa_dir: &Path,
) -> Result<Vec<fastvep_sa::gene::GeneIndex>> {
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
        let results = ctx
            .annotate_vcf_text_with_acmg(vcf, false, None)
            .unwrap();
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
}
