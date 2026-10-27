//! `fastvep annotate`: the per-variant annotation pipeline.

use super::cache_build::{load_one_gff3, parse_gff3_arg, Gff3Spec};
use super::open_vcf_input_reader;
use anyhow::{Context, Result};
use fastvep_annotate::pick::{
    apply_pick, parse_pick_order, PickCriterion, PickFlags, PickRequest, DEFAULT_PICK_ORDER,
};
use fastvep_cache::annotation::{AnnotationProvider, AnnotationValue, GeneAnnotationProvider};
use fastvep_cache::fasta::ReferenceFasta;
use fastvep_cache::info::CacheInfo;
use fastvep_cache::providers::{
    FastaSequenceProvider, IndexedTranscriptProvider, KnownSite, MatchedVariant, SequenceProvider,
    TabixVariationProvider, TranscriptProvider, VariationProvider,
};
use fastvep_cache::transcript_cache::StaleCacheFormat;
use fastvep_consequence::ConsequencePredictor;
use fastvep_core::Allele;
use fastvep_genome::Transcript;
use fastvep_hgvs;
use fastvep_io::output::{self, ColumnFlags};
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
    annotate_intergenic, annotate_sa_only_scaffold, hgvs_allele_on_transcript, hgvsg_per_alt,
    load_gene_providers, load_sa_providers, splice_context, RowTranscripts,
};

const BATCH_SIZE: usize = 1024;

/// VEP flags that change what is written without being a column of their own.
///
/// `--force_overwrite` is not here: VEP refuses to write over an existing file without it
/// and fastVEP always writes over it, so the flag changes nothing and is only accepted.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct CompatFlags {
    /// `--total_length`: `cDNA_position`, `CDS_position` and `Protein_position` carry the length
    /// they are counted in, `14/203`.
    pub total_length: bool,
    /// `--no_escape`: the `=` of a synonymous `HGVSp` is not written `%3D` and the space of a VCF
    /// `SYMBOL` is `_` and not `%20`.
    pub no_escape: bool,
}

/// Transcripts to leave out of a run, which Ensembl VEP's cache does not hold.
///
/// A GFF3 has every transcript of its annotation. VEP's caches are built from a database, and the
/// build leaves some out (`AnnotationSource::Database::Transcript`): readthrough transcripts and
/// those of biotype `artifact`, from an Ensembl cache (109 of the 23,382 transcripts of chr21 and
/// chr22 in Ensembl's GFF3, and the only ones of them the cache lacks), and `--exclude_predicted`
/// asks VEP for RefSeq's models. None is on by default: dropping real transcripts from a user's
/// annotation is theirs to decide, where reading the file as it is is not a decision.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct TranscriptExclusions {
    /// `--exclude-predicted`: a stable id that matches `^X._` (`XM_`, `XR_`).
    pub predicted: bool,
    /// `--exclude-readthrough`: tagged `readthrough_transcript`.
    pub readthrough: bool,
    /// `--exclude-artifact`: biotype `artifact`.
    pub artifact: bool,
}

impl TranscriptExclusions {
    /// Whether any is on.
    pub fn any(self) -> bool {
        self != Self::default()
    }

    /// Whether `transcript` is one of the kinds this run leaves out.
    pub fn excludes(self, transcript: &Transcript) -> bool {
        // VEP's pattern: an `X`, any one character, an underscore.
        let predicted = matches!(transcript.stable_id.as_bytes(), [b'X', _, b'_', ..]);
        (self.predicted && predicted)
            || (self.readthrough && transcript.readthrough)
            || (self.artifact && &*transcript.biotype == "artifact")
    }
}

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
    /// The `--pick*` / `--flag-pick*` switch this run was given, if any.
    pub pick: PickFlags,
    pub hgvs: bool,
    /// The optional tab and CSQ columns VEP adds on request (`--hgvsg`,
    /// `--biotype`, `--numbers`, ...). `hgvs` is separate: it also decides
    /// whether the HGVS values are computed at all.
    pub columns: ColumnFlags,
    /// `--total_length` and `--no_escape`.
    pub compat: CompatFlags,
    /// The kinds of transcript to leave out ([`TranscriptExclusions`]).
    pub exclude: TranscriptExclusions,
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

/// Everything the per-variant work needs that is fixed for the whole run.
///
/// The annotation loop runs `annotate_variant` across rayon threads, so every
/// field here is a shared borrow of something built once during startup.
/// Grouping them is what lets that work be a named function: it used to be the
/// body of a 767-line closure nested inside `run_annotate`, indented far enough
/// that rustfmt left it alone.
struct AnnotationContext<'a> {
    config: &'a AnnotateConfig,
    /// `--sa-only`: emit supplementary annotations and nothing else.
    sa_only: bool,
    transcript_provider: &'a IndexedTranscriptProvider,
    seq_provider: Option<&'a dyn SequenceProvider>,
    predictor: &'a ConsequencePredictor,
    sa_providers: &'a [Box<dyn AnnotationProvider>],
    gene_providers: &'a [fastvep_sa::gene::GeneIndex],
    acmg_config: Option<&'a fastvep_classification::AcmgConfig>,
    functional_evidence: Option<&'a fastvep_classification::FunctionalEvidenceIndex>,
    pick_order: &'a [PickCriterion],
    /// `None` for a run that picks nothing.
    pick: Option<PickRequest>,
    sample_names: &'a [String],
    /// gnomAD queries need the VCF-style allele rather than fastVEP's
    /// normalised one; set when a gnomAD source is among `sa_providers`.
    normalize_gnomad_queries: bool,
    /// A tandem-repeat database is loaded, so PM4/BP3 can consult it.
    repeat_db_loaded: bool,
}

/// Annotate one variant: transcript overlap, consequence prediction, HGVS,
/// supplementary annotations, `--pick`, and ACMG classification.
///
/// Called once per variant from the parallel phase of `run_annotate`.
/// `known` carries the co-located known variants the sequential phase already
/// looked up for this record.
fn annotate_variant(ctx: &AnnotationContext<'_>, vf: &mut VariationFeature, known: &KnownLookup) {
    let matched_by_allele = &known.by_allele;
    // Destructured rather than used through `ctx.` so the body below reads the
    // same as it did when these were the enclosing function's locals.
    let &AnnotationContext {
        config,
        sa_only,
        transcript_provider,
        seq_provider,
        predictor,
        sa_providers,
        gene_providers,
        acmg_config,
        functional_evidence,
        pick_order,
        pick,
        sample_names,
        normalize_gnomad_queries,
        repeat_db_loaded,
    } = ctx;

    // VEP's known-variant columns, worked out once per alt allele: every
    // transcript row of an allele writes the same cells, so the rows borrow them.
    if config.columns.known.checks_existing() {
        vf.known = fastvep_annotate::known_cells(
            &known.site,
            &vf.ref_allele.to_string(),
            &known.alleles,
            config.columns.known,
        );
    }

    // Kept for the classification pass at the end, which runs after `--pick`
    // has deleted rows and needs the transcript each survivor came from.
    let mut overlapping = Vec::new();

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
        overlapping = transcript_provider
            .get_transcripts(chrom, query_start, query_end)
            .unwrap_or_default();

        if overlapping.is_empty() {
            // Intergenic
            annotate_intergenic(vf);
            // The variant's own HGVSg is there whether or not a transcript is:
            // VEP's intergenic output ends in the function that writes HGVSg
            // (`IntergenicVariationAllele_to_output_hash` returns
            // `VariationFeatureOverlapAllele_to_output_hash`). Only `--hgvsg`
            // asks for it here, as it does on the transcript rows; `fastvep-annotate`'s
            // loop has no flag of its own for it and writes it under `hgvs` on both.
            if config.columns.hgvsg {
                let hgvsg = hgvsg_per_alt(
                    seq_provider.map(|sp| sp as &dyn SequenceProvider),
                    &vf.position.chromosome,
                    vf.position.start,
                    vf.position.end,
                    &vf.ref_allele,
                    &vf.alt_alleles,
                );
                for tv in &mut vf.transcript_variations {
                    for aa in &mut tv.allele_annotations {
                        aa.hgvsg = vf
                            .alt_alleles
                            .iter()
                            .position(|alt| *alt == aa.allele)
                            .and_then(|i| hgvsg[i].clone());
                    }
                }
            }
            // Populate existing_variation on intergenic annotations too
            for tv in &mut vf.transcript_variations {
                for aa in &mut tv.allele_annotations {
                    if let Some(matches) = matched_by_allele.get(&aa.allele.to_string()) {
                        aa.existing_variation = matches.iter().map(|m| m.name.clone()).collect();
                    }
                }
            }
        } else {
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
                // No reference window is read for the predictor: it is given none to use (the
                // transcripts carry their own sequence), and reading the variant's window of
                // `--distance` bases either side for every variant was a few kilobytes copied out
                // of the FASTA per record for nothing.
                let result = predictor.predict(
                    &vf.position,
                    &vf.ref_allele,
                    &vf.alt_alleles,
                    &overlapping,
                    None,
                );
                result.transcript_consequences
            };

            // The reference an HGVS description is written against: the genome's, where the VCF
            // has `N` and where its REF is not the reference's, as Ensembl does for a
            // substitution (`reference_for_hgvs`). The `g.` description reads it too, so that the
            // two say the same thing about a record whose REF is wrong.
            let hgvs_ref = fastvep_annotate::reference_for_hgvs(
                seq_provider.map(|sp| sp as &dyn SequenceProvider),
                chrom,
                vf.position.start,
                vf.position.end,
                &vf.ref_allele,
            );

            // HGVSg is the variant's, the same on every transcript row, and reading
            // the reference to place it is not free: made once per record, and only
            // when an output carries it (the `HGVSg` column or field, or JSON's key,
            // which `--hgvs` has always filled).
            let hgvsg_by_alt =
                if config.columns.hgvsg || (config.hgvs && config.output_format == "json") {
                    hgvsg_per_alt(
                        seq_provider.map(|sp| sp as &dyn SequenceProvider),
                        chrom,
                        vf.position.start,
                        vf.position.end,
                        &hgvs_ref,
                        &vf.alt_alleles,
                    )
                } else {
                    Vec::new()
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

                // `*` and the placeholders (`.`, `<NON_REF>`, `<*>`) are not variants of
                // the transcript, and VEP writes no row for them: a record whose only ALT
                // is one has no row here, and the real ALTs beside one keep theirs.
                let allele_annotations: Vec<AlleleAnnotation> = tc
                    .allele_consequences
                    .iter()
                    .filter(|ac| ac.allele != Allele::Missing)
                    .map(|ac| {
                        let mut ann = AlleleAnnotation {
                            allele: ac.allele.clone(),
                            consequences: ac.consequences.clone(),
                            impact: ac.impact,
                            location: ac.location,
                            protein_position: ac.protein_range(),
                            amino_acids: ac.amino_acids.clone(),
                            codons: ac.codons.clone(),
                            exon: ac.exon,
                            splice: Default::default(),
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
                            pick: false,
                        };

                        // `--total_length` is a matter of the cells and of nothing that is
                        // worked out from them: `ac.location`, which the HGVS code reads, has
                        // no totals. A structural variant's cells are not written the way
                        // VEP writes them, with or without the flag, so it gets none.
                        if config.compat.total_length && !vf.variant_type.is_structural() {
                            if let Some(tr) = transcript {
                                tr.apply_total_length(
                                    &mut ann.location,
                                    vf.position.start,
                                    vf.position.end,
                                );
                            }
                        }

                        // HGVSg belongs to the variant and needs no transcript, so
                        // `--hgvsg` has it without `--hgvs`.
                        if !hgvsg_by_alt.is_empty() {
                            ann.hgvsg = vf
                                .alt_alleles
                                .iter()
                                .position(|alt| *alt == ac.allele)
                                .and_then(|i| hgvsg_by_alt[i].clone());
                        }

                        // How far the HGVS description moved 3' from the VCF's own position
                        let mut shift = 0u64;

                        // Generate HGVS if requested
                        if config.hgvs {
                            // An ALT that is not plain bases (`N`) has no HGVSc, as in
                            // Ensembl (`return undef if variation_feature_seq =~
                            // m/[^ACGT\-]/`).
                            if let Some(tr) = transcript.filter(|_| ac.allele.is_plain_bases()) {
                                // Build versioned IDs for HGVS notation
                                let versioned_tid = match tr.version {
                                    Some(v) => format!("{}.{}", tc.transcript_id, v),
                                    None => tc.transcript_id.to_string(),
                                };

                                // Read the way HGVS reads it: clipped of the
                                // bases the pair repeats at either end, and
                                // facing the transcript. A multi-allelic record
                                // is trimmed once for the whole site, so an
                                // allele that reaches here as a replacement may
                                // still be a plain deletion.
                                let hv = hgvs_allele_on_transcript(
                                    seq_provider.map(|sp| sp as &dyn SequenceProvider),
                                    chrom,
                                    tr,
                                    vf.position.start,
                                    vf.position.end,
                                    &hgvs_ref,
                                    &ac.allele,
                                    // `zip`: a pair with one end mapped is a variant with
                                    // one end in an exon and the other in an intron, which
                                    // has no cDNA span at all. Collapsing it to a point
                                    // names one base of a change that covers seventeen.
                                    // (The `cDNA_position` cell says `601-?` for it.)
                                    ac.cdna_start
                                        .zip(ac.cdna_end)
                                        .map(|(a, b)| (a.min(b), a.max(b))),
                                );
                                if let Some(hv) = hv {
                                    // Every indel is 3'-shifted over the genome, so
                                    // it routes through `hgvsc_shifted` whether or
                                    // not both its ends are exonic: the shift
                                    // crosses splice sites, and the description has
                                    // to follow it there. Without a reference
                                    // there is no walk to make, so those fall to the
                                    // cDNA renderer below, which shifts on the
                                    // spliced sequence instead - not the same answer
                                    // across a splice site, but the only one
                                    // available. `hgvsc_shifted` also returns `None`
                                    // for a variant it cannot place at all, and that
                                    // falls the same way.
                                    let shifted = (fastvep_annotate::is_shiftable_indel(
                                        &hv.hgvs_ref,
                                        &hv.hgvs_alt,
                                    ) && seq_provider.is_some())
                                    .then(|| {
                                        fastvep_annotate::hgvsc_shifted_outcome(
                                            seq_provider.map(|sp| sp as &dyn SequenceProvider),
                                            chrom,
                                            tr,
                                            &versioned_tid,
                                            hv.start,
                                            hv.end,
                                            hv.anchor_start,
                                            &hv.genomic_ref,
                                            &hv.genomic_alt,
                                            &hv.hgvs_ref,
                                            &hv.hgvs_alt,
                                            tr.cdna_coding_start,
                                            tr.cdna_coding_end,
                                            hv.shift,
                                        )
                                    });
                                    ann.hgvsc = match (shifted, tr.cdna_coding_start) {
                                        (
                                            Some(fastvep_annotate::ShiftedHgvsc::Written(h, moved)),
                                            _,
                                        ) => {
                                            shift = moved;
                                            Some(h)
                                        }
                                        // A description the 3'-rule carries off the end of the transcript, or one for a
                                        // change that is not a plain overlap of it: Ensembl writes nothing, and the cDNA
                                        // renderer below must not write what Ensembl would not.
                                        (
                                            Some(fastvep_annotate::ShiftedHgvsc::NoDescription),
                                            _,
                                        ) => None,
                                        (
                                            Some(fastvep_annotate::ShiftedHgvsc::CannotRender)
                                            | None,
                                            Some(coding_start),
                                        ) => match hv.cdna {
                                            Some((cs, ce)) => fastvep_hgvs::hgvsc_with_seq(
                                                &versioned_tid,
                                                cs,
                                                ce,
                                                &hv.hgvs_ref,
                                                &hv.hgvs_alt,
                                                coding_start,
                                                tr.cdna_coding_end,
                                                tr.spliced_seq.as_deref(),
                                                tr.codon_table_start_phase,
                                            ),
                                            // Not "is this intronic": a variant
                                            // with one end in an exon and the
                                            // other in an intron has no exonic
                                            // cDNA pair, and `intron_at` reads
                                            // its first base only.
                                            None => fastvep_annotate::hgvsc_shifted_with_shift(
                                                seq_provider.map(|sp| sp as &dyn SequenceProvider),
                                                chrom,
                                                tr,
                                                &versioned_tid,
                                                hv.start,
                                                hv.end,
                                                hv.anchor_start,
                                                &hv.genomic_ref,
                                                &hv.genomic_alt,
                                                &hv.hgvs_ref,
                                                &hv.hgvs_alt,
                                                Some(coding_start),
                                                tr.cdna_coding_end,
                                                hv.shift,
                                            )
                                            .map(|(h, moved)| {
                                                shift = moved;
                                                h
                                            }),
                                        },
                                        // Non-coding transcript: use n. notation
                                        (
                                            Some(fastvep_annotate::ShiftedHgvsc::CannotRender)
                                            | None,
                                            None,
                                        ) => match hv.cdna {
                                            Some((cs, ce)) => fastvep_hgvs::hgvsc_noncoding(
                                                &versioned_tid,
                                                cs,
                                                ce,
                                                &hv.hgvs_ref,
                                                &hv.hgvs_alt,
                                                tr.spliced_seq.as_deref(),
                                            ),
                                            None => fastvep_annotate::hgvsc_shifted_with_shift(
                                                seq_provider.map(|sp| sp as &dyn SequenceProvider),
                                                chrom,
                                                tr,
                                                &versioned_tid,
                                                hv.start,
                                                hv.end,
                                                hv.anchor_start,
                                                &hv.genomic_ref,
                                                &hv.genomic_alt,
                                                &hv.hgvs_ref,
                                                &hv.hgvs_alt,
                                                None,
                                                None,
                                                hv.shift,
                                            )
                                            .map(|(h, moved)| {
                                                shift = moved;
                                                h
                                            }),
                                        },
                                    };
                                }
                            }

                            // HGVSp
                            if let Some(tr) = transcript {
                                ann.hgvsp = fastvep_annotate::hgvsp_for_allele(
                                    &fastvep_annotate::HgvspRequest {
                                        predictor,
                                        transcript: tr,
                                        provider: seq_provider
                                            .map(|sp| sp as &dyn SequenceProvider),
                                        chrom,
                                        start: vf.position.start,
                                        end: vf.position.end,
                                        ref_allele: &vf.ref_allele,
                                        alt_allele: &ac.allele,
                                        location: &ac.location,
                                        consequences: &ac.consequences,
                                        amino_acids: ac
                                            .amino_acids
                                            .as_ref()
                                            .map(|(r, a)| (r.as_str(), a.as_str())),
                                        coding_sequence_kept: ac.coding_sequence_kept,
                                    },
                                );
                            }
                        }

                        // VEP's HGVS_OFFSET: the distance the description moved,
                        // signed by the strand (`offset * strand`), and only when
                        // there is a description to have moved.
                        if shift > 0 && (ann.hgvsc.is_some() || ann.hgvsp.is_some()) {
                            let strand = transcript.map_or(1, |t| t.strand.as_int());
                            ann.hgvs_offset = Some(shift as i64 * i64::from(strand));
                        }

                        ann
                    })
                    .collect();

                if allele_annotations.is_empty() {
                    continue;
                }

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

    // Apply the pick before SA/gene/ACMG so those passes only run on the
    // surviving transcripts. Running it after would still produce correct
    // output but would waste the most expensive work (ACMG classification)
    // on transcripts that get thrown away.
    //
    // That saving is exactly what a `--flag-pick*` run gives up: flagging
    // retains every transcript by definition, so classification runs on all
    // of them. The order still matters for the reducing modes, and putting
    // the two in the same place keeps the flag and the reduction deciding
    // from identical inputs.
    if let Some(request) = pick.filter(|_| !sa_only) {
        apply_pick(&mut vf.transcript_variations, &request.plan(pick_order));
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
        let mut allele_results: HashMap<String, Vec<(String, String)>> = HashMap::new();
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
                    seq_provider.map(|sp| {
                        fastvep_cache::normalize::normalize_variant(
                            sp, chrom, query_pos, &ref_str, &alt_str,
                        )
                    })
                } else {
                    None
                };
                let mut results: Vec<(String, String)> = Vec::new();
                for sa in sa_providers {
                    let (sa_pos, sa_ref, sa_alt) = if sa.metadata().match_by_allele {
                        if sa.json_key() == "gnomad" {
                            match &gnomad_norm {
                                Some(n) => (n.pos, n.ref_allele.as_str(), n.alt_allele.as_str()),
                                None => (query_pos, ref_str.as_str(), alt_str.as_str()),
                            }
                        } else {
                            (query_pos, ref_str.as_str(), alt_str.as_str())
                        }
                    } else {
                        (vf.position.start, "", "")
                    };
                    // A failed lookup is not a miss: see `SaLookupErrors`.
                    match sa.annotate_position(chrom, sa_pos, sa_ref, sa_alt) {
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
                        Err(e) => fastvep_annotate::sa_lookup_errors().record(sa.json_key(), &e),
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
                    for gp in gene_providers {
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
        let trio_genotypes = extract_trio_genotypes_cli(vf, acmg_cfg, sample_names);

        let functional_by_alt = resolve_functional_by_alt(functional_evidence, vf);
        // Resolved here rather than inside the classifier: the ClinVar
        // splice index is keyed by genomic coordinate, and
        // `extract_classification_input` is handed a transcript-level
        // view that carries none.
        let query_alleles = vf.query_alleles();
        let mut row_transcripts = RowTranscripts::new(&overlapping);

        for tv in &mut vf.transcript_variations {
            let gene_sym = tv.gene_symbol.as_deref().unwrap_or("");
            let gene_anns: Vec<&fastvep_core::GeneAnnotation> = vf
                .gene_annotations
                .iter()
                .filter(|ga| ga.gene_symbol == gene_sym)
                .collect();
            let transcript = row_transcripts.of(&tv.transcript_id);
            for aa in &mut tv.allele_annotations {
                let alt_idx = vf.alt_alleles.iter().position(|a| *a == aa.allele);
                if let Some(tr) = transcript {
                    aa.splice = splice_context(
                        seq_provider.map(|sp| sp as &dyn SequenceProvider),
                        &vf.position.chromosome,
                        tr,
                        vf.position.start,
                        vf.position.end,
                        &vf.ref_allele,
                        &aa.allele,
                        &aa.consequences,
                    );
                }
                let input = fastvep_classification::extract_classification_input(
                    &aa.consequences,
                    aa.impact,
                    tv.gene_symbol.as_deref(),
                    tv.canonical,
                    aa.amino_acids.as_ref(),
                    aa.protein_position.map(|(s, _)| s),
                    aa.hgvsc.as_deref(),
                    aa.splice,
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

    if !sa_only {
        vf.compute_most_severe();
    }
}

/// Everything between the `--gff3` / `--transcript-cache` arguments and a
/// transcript set in memory: the sidecar-cache decision, the fallback to
/// reading GFF3, and the merge when several sources were given.
///
/// The `bool` says whether any source came through the tabix path. A tabix read
/// returns only the features overlapping *this* VCF's variants, which makes the
/// result right for this input and wrong for every other one, so the caller
/// uses it to refuse to persist the set as a sidecar cache.
fn load_transcript_models(
    config: &AnnotateConfig,
    gff3_specs: &[Gff3Spec],
    single_gff3: Option<&Gff3Spec>,
    cache_path: Option<&Path>,
) -> Result<(Vec<Transcript>, bool)> {
    let mut region_restricted = false;
    let transcripts = 'load: {
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
        if let Some(cp) = cache_path {
            if cp.exists() {
                // Freshness check: only sidecar-cache mode does it (the
                // user-provided --transcript-cache is always trusted).
                let is_fresh = if explicit_cache {
                    true
                } else {
                    single_gff3
                        .map(|s| {
                            fastvep_cache::transcript_cache::cache_is_fresh(cp, Path::new(&s.path))
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
                            let stale = e.downcast_ref::<StaleCacheFormat>();
                            if explicit_cache {
                                // A cache rejected for its format (pre-#90,
                                // pre-#98) is intact, so "truncated or
                                // corrupt" is the wrong diagnosis and would
                                // send the user looking for a disk problem
                                // they do not have. Say which of the two it
                                // is, once, and let StaleCacheFormat supply
                                // the reason and the evidence for the format
                                // it found.
                                //
                                // The recovery clause is written here rather
                                // than inside StaleCacheFormat because it is
                                // specific to this path: the file was named on
                                // the command line, so deleting it produces no
                                // sidecar to take its place.
                                let why = match stale {
                                    Some(stale) => format!(
                                        "{stale} Rebuild it with `fastvep cache`, or drop \
                                         --transcript-cache and pass --gff3 so a fresh \
                                         sidecar is built alongside it."
                                    ),
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
                            // The sidecar rebuilds from the GFF3 on the next
                            // few lines, so this asks nothing of the user and
                            // says so in one line. A stale format is an
                            // expected consequence of upgrading, not a fault;
                            // the evidence behind it goes to the log for
                            // anyone who wants it.
                            match stale {
                                Some(stale) => {
                                    eprintln!(
                                        "Note: sidecar cache {} {}; rebuilding from GFF3",
                                        cp.display(),
                                        stale.summary()
                                    );
                                    log::debug!("{}", stale.detail);
                                }
                                None => eprintln!(
                                    "Warning: sidecar cache {} could not be loaded ({e}); rebuilding from GFF3",
                                    cp.display()
                                ),
                            }
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
            for spec in gff3_specs {
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
    };
    Ok((transcripts, region_restricted))
}

/// Open the reference FASTA named by `--fasta`. A `.fai` is used only when
/// the file is uncompressed; a gzipped FASTA is decompressed into memory.
fn open_sequence_provider(config: &AnnotateConfig) -> Result<Option<Box<dyn SequenceProvider>>> {
    let Some(fasta_path) = config.fasta.as_deref() else {
        return Ok(None);
    };
    Ok({
        match ReferenceFasta::open(Path::new(fasta_path))? {
            ReferenceFasta::Mapped(reader) => {
                eprintln!(
                    "Memory-mapped reference FASTA from {} (using .fai index)",
                    fasta_path
                );
                Some(Box::new(
                    fastvep_cache::providers::MmapFastaSequenceProvider::new(reader),
                ))
            }
            ReferenceFasta::Loaded {
                reader,
                ignored_fai,
            } => {
                if ignored_fai {
                    eprintln!(
                        "Ignoring {fasta_path}.fai: the FASTA is gzip-compressed, and a .fai addresses the uncompressed file. Loading it into memory."
                    );
                }
                eprintln!("Loaded reference FASTA from {}", fasta_path);
                Some(Box::new(FastaSequenceProvider::new(reader)))
            }
        }
    })
}

/// Fill in the spliced and protein sequences of transcripts that arrived
/// without them - from a cache written before sequences were built, or a GFF3
/// read on a run with no `--fasta`.
///
/// A non-coding transcript needs its spliced sequence only for HGVS: the
/// 3'-rule and `dup` collapsing are read off it, and without it an `n.`
/// coordinate is written where the VCF put it rather than where HGVS asks for
/// it - 2,800 rows over a 6,600-variant ClinVar sample. Nothing else in the
/// pipeline reads it, and building it for every lncRNA and retained-intron
/// transcript roughly doubles this step, so `for_hgvs` keeps that cost on the
/// runs that spend it.
///
/// Returns whether any transcript needed the work, which is what tells the
/// caller the sidecar cache on disk no longer matches what is in memory.
fn attach_transcript_sequences(
    transcripts: &mut [Transcript],
    seq_provider: Option<&dyn SequenceProvider>,
    for_hgvs: bool,
) -> bool {
    let wanted = |t: &Transcript| (for_hgvs || t.is_coding()) && t.spliced_seq.is_none();
    let needs_seq_build = transcripts.iter().any(&wanted);
    if needs_seq_build {
        if let Some(sp) = seq_provider {
            let built = AtomicUsize::new(0);
            transcripts.par_iter_mut().for_each(|tr| {
                if wanted(tr) {
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
        } else {
            // Without the CDS there is no way to tell a missense change from a
            // synonymous one, so `terms_for_window` reports the honest
            // `coding_sequence_variant` - which is MODIFIER, where a missense
            // call would be MODERATE. Anything filtering on IMPACT drops those
            // rows, so say so rather than letting the tier change go unseen.
            eprintln!(
                "Warning: no reference sequence (--fasta not given, and the transcript \
                 source carries none), so coding changes are reported as \
                 coding_sequence_variant with MODIFIER impact rather than \
                 missense/synonymous. Pass --fasta for coding consequences."
            );
        }
    }
    needs_seq_build
}

/// The contigs a run's gene models and reference are about, so that a record on any other is
/// said to be, and not taken for one in an intergenic region.
///
/// A VCF whose contigs are spelled another way than the GFF3's (`Chr10` and `chr10`: `chr1` and `1`
/// are matched, case is not) used to come out as one `intergenic_variant` row per record, at exit
/// code 0 and with nothing said, which is also the answer for a variant that is nowhere near a
/// gene, and reads the same. Ensembl VEP skips such a record and warns (`Chromosome chrZ not found
/// in annotation sources or synonyms`). This warns in the same way, once for each contig and with a
/// summary at the end, naming the contig of the run that is the same in another case when there is
/// one, and goes on as before with the record, because a `--sa-dir` or a variation cache can still
/// say something about a variant on a contig with no gene model, and VEP has neither. A run in which
/// not one record is on a known contig, and which has no other source to give it a use, is an error:
/// a run that annotates nothing is a naming mistake.
/// How many contigs of a VCF that the run does not know it says so for, one by one; the rest are in
/// the summary. (A GRCh38 analysis set has thousands of alternate contigs.)
const CONTIGS_SAID: usize = 10;

struct ContigCheck {
    known: std::collections::HashSet<String>,
    /// The contig of the last record and whether it is known: a VCF is mostly runs of one contig,
    /// and the aliases of a name are not worth building for each record.
    last: Option<(String, bool)>,
    unknown: std::collections::BTreeMap<String, u64>,
    on_known: u64,
}

impl ContigCheck {
    fn new(
        transcripts: &[Transcript],
        gff3_specs: &[Gff3Spec],
        seq_provider: Option<&dyn SequenceProvider>,
    ) -> Self {
        let mut known: std::collections::HashSet<String> = transcripts
            .iter()
            .map(|t| t.chromosome.to_string())
            .collect();
        // A tabix read keeps only the transcripts near the input's variants, so the ones it has are
        // not the contigs the file is about: its index is.
        for spec in gff3_specs {
            if Path::new(&format!("{}.tbi", spec.path)).exists() {
                if let Ok(names) = fastvep_cache::gff::tabix_sequence_names(Path::new(&spec.path)) {
                    known.extend(names);
                }
            }
        }
        if let Some(names) = seq_provider.and_then(|sp| sp.sequence_names()) {
            known.extend(names);
        }
        Self {
            known,
            last: None,
            unknown: Default::default(),
            on_known: 0,
        }
    }

    fn is_known(&self, chrom: &str) -> bool {
        fastvep_core::chrom_aliases(chrom)
            .iter()
            .any(|alias| self.known.contains(alias))
    }

    /// Note a record on `chrom`, and say so, once for the contig and for the first few, if the run
    /// does not know the contig.
    fn note(&mut self, chrom: &str) {
        let known = match &self.last {
            Some((last, known)) if last == chrom => *known,
            _ => {
                let known = self.is_known(chrom);
                self.last = Some((chrom.to_string(), known));
                known
            }
        };
        if known {
            self.on_known += 1;
            return;
        }
        let count = self.unknown.entry(chrom.to_string()).or_default();
        *count += 1;
        if *count == 1 && self.unknown.len() <= CONTIGS_SAID {
            let same_name = self
                .known
                .iter()
                .filter(|k| k.eq_ignore_ascii_case(chrom))
                .min()
                .map(|k| format!(" (the gene models and the reference have `{k}`; contig names are case-sensitive)"))
                .unwrap_or_default();
            eprintln!(
                "warning: the VCF has records on `{chrom}`, which is not a contig of the gene models or the reference{same_name}; they are annotated as intergenic"
            );
        }
    }

    /// The summary of a run, and the error of one that annotated nothing: `other_source` is whether
    /// the run has a `--sa-dir` or a variation cache, which can annotate a variant that no gene
    /// model reaches.
    fn finish(&self, other_source: bool) -> Result<()> {
        if self.unknown.is_empty() {
            return Ok(());
        }
        let list = |names: Vec<(&String, String)>| {
            let mut shown: Vec<String> = names
                .iter()
                .take(8)
                .map(|(n, extra)| format!("{n}{extra}"))
                .collect();
            if names.len() > 8 {
                shown.push(format!("and {} more", names.len() - 8));
            }
            shown.join(", ")
        };
        let skipped: u64 = self.unknown.values().sum();
        let unknown = list(
            self.unknown
                .iter()
                .map(|(n, c)| (n, format!(" ({c})")))
                .collect(),
        );
        if self.on_known == 0 && !other_source {
            let mut known: Vec<&String> = self.known.iter().collect();
            known.sort();
            let known = list(known.into_iter().map(|n| (n, String::new())).collect());
            return Err(anyhow::anyhow!(
                "none of the {skipped} records of the VCF is on a contig of the gene models or the reference: \
                 the VCF has {unknown}, and they have {known}. `chr1` and `1`, and `chrM` and `MT`, are \
                 matched; case is not. Rename the contigs of one of them."
            ));
        }
        eprintln!(
            "warning: {skipped} record(s) are on contigs that the gene models and the reference do not have, and were called intergenic: {unknown}"
        );
        Ok(())
    }
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
    // `--cache-dir` is the variation cache, and the known variants are all fastVEP reads from it,
    // so naming one is asking for them: it is VEP's `--check_existing`, which adds `CLIN_SIG`,
    // `SOMATIC` and `PHENO` beside `Existing_variation`. VEP needs the flag because its `--cache`
    // is also the transcript cache; a flag that cannot do anything without a second one would be
    // a trap here. The frequencies and PubMed ids still need theirs, as in VEP. JSON has no place
    // for these columns, and a flag the user did not type is not one to warn about.
    if config.cache_dir.is_some() && !config.sa_only && config.output_format != "json" {
        config.columns.known.check_existing = true;
    }
    let config = config;
    let sa_only = config.sa_only;
    if sa_only {
        if config.sa_dir.is_none() {
            return Err(anyhow::anyhow!(
                "--sa-only requires --sa-dir to be set (otherwise there is nothing to emit)."
            ));
        }
        let pick_switch = config.pick.requested_option();
        for (set, name) in [
            (!config.gff3.is_empty(), "--gff3"),
            (config.fasta.is_some(), "--fasta"),
            (config.cache_dir.is_some(), "--cache-dir"),
            (config.transcript_cache.is_some(), "--transcript-cache"),
            (config.acmg, "--acmg"),
            (config.hgvs, "--hgvs"),
            (config.columns.hgvsg, "--hgvsg"),
            (config.columns.symbol, "--symbol"),
            (config.columns.biotype, "--biotype"),
            (config.columns.canonical, "--canonical"),
            (config.columns.numbers, "--numbers"),
            (config.columns.variant_class, "--variant-class"),
            (config.columns.show_ref_allele, "--show-ref-allele"),
            (config.compat.total_length, "--total-length"),
            (config.compat.no_escape, "--no-escape"),
            // Names the switch the user actually typed, so the warning is not
            // about `--pick` when they passed `--flag-pick-allele-gene`.
            (pick_switch.is_some(), pick_switch.unwrap_or("--pick")),
            (config.proband.is_some(), "--proband"),
            (config.mother.is_some(), "--mother"),
            (config.father.is_some(), "--father"),
        ] {
            if set {
                eprintln!("warning: --sa-only is set; ignoring {}", name);
            }
        }
        for name in config.columns.known.requested() {
            eprintln!("warning: --sa-only is set; ignoring {}", name);
        }
    } else if config.columns.known.checks_existing() && config.cache_dir.is_none() {
        // The known variants are only in VEP's variation cache. Writing the
        // columns without it would print `-` for every variant, which reads as
        // "none of these is known" and is not an answer.
        return Err(anyhow::anyhow!(
            "{} read the known variants of VEP's variation cache, and no cache was given. \
             Pass --cache-dir DIR, the directory with a {{chr}}/all_vars.gz file for each \
             chromosome (the indexed cache VEP's own `--cache --offline` runs from).",
            config.columns.known.requested().join(", ")
        ));
    }

    // JSON has no place for them: VEP's own JSON carries the known variants in a
    // `colocated_variants` block that this does not write, and the tab file and the
    // VCF `CSQ` are where these columns are.
    if config.output_format == "json" {
        for name in config.columns.known.requested() {
            eprintln!(
                "warning: --output-format json does not carry the known-variant columns; ignoring {}",
                name
            );
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

    let (mut transcripts, region_restricted) = if sa_only {
        (Vec::new(), false)
    } else {
        load_transcript_models(&config, &gff3_specs, single_gff3, cache_path.as_deref())?
    };

    // Load FASTA reference (prefer mmap with .fai index, fall back to in-memory).
    // Skipped in --sa-only mode.
    let seq_provider: Option<Box<dyn SequenceProvider>> = if sa_only {
        None
    } else {
        open_sequence_provider(&config)?
    };

    let needs_seq_build =
        attach_transcript_sequences(&mut transcripts, seq_provider.as_deref(), config.hgvs);

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

    // After the sidecar cache is written, and never in the parser: the cache holds the file's
    // whole transcript set, so a run that leaves some out must not write what it kept as the
    // cache of the next one, which would not have them.
    if config.exclude.any() {
        let before = transcripts.len();
        transcripts.retain(|t| !config.exclude.excludes(t));
        eprintln!(
            "Left out {} of {} transcripts ({})",
            before - transcripts.len(),
            before,
            [
                (config.exclude.predicted, "--exclude-predicted"),
                (config.exclude.readthrough, "--exclude-readthrough"),
                (config.exclude.artifact, "--exclude-artifact"),
            ]
            .iter()
            .filter(|(on, _)| *on)
            .map(|(_, flag)| *flag)
            .collect::<Vec<_>>()
            .join(" "),
        );
        if transcripts.is_empty() {
            return Err(anyhow::anyhow!(
                "Every transcript was left out by the exclusion flags. Refusing to continue, because \
                 every variant would be reported as intergenic."
            ));
        }
    }

    let mut contig_check = if sa_only {
        None
    } else {
        Some(ContigCheck::new(
            &transcripts,
            &gff3_specs,
            seq_provider.as_deref(),
        ))
    };

    let transcript_provider = IndexedTranscriptProvider::new(transcripts);

    // Initialize variation provider from VEP cache if provided.
    // Skipped in --sa-only mode (existing_variation comes from defaults).
    //
    // Skipped for JSON as well, which has no place for what it finds: `--cache-dir` read every variant's
    // window of the cache and wrote nothing of it.
    let var_provider: Option<TabixVariationProvider> = if sa_only || config.output_format == "json"
    {
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
        warn_about_missing_columns(config.columns.known, &cache_info.variation_cols);
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

    if let Some(conflict) = fastvep_annotate::check_sa_assemblies(&sa_providers, &gene_providers) {
        eprintln!("warning: {conflict}");
    }

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

    // VEP's own JSON writer makes `cdna_start` and `cdna_end` by splitting the position cell at its
    // first `-`: with `--total_length` that cell is `5-7/1203`, and the `_end` it keeps is the
    // start, `5`, because `7/1203` is neither a number nor has a `?`, and a single position is
    // the string `"5/1203"` for both. Reproducing that would write a wrong end on purpose, so the
    // totals are left out of JSON, as the warning says. Measured on VEP 116.3 over the edge
    // fixture's 1,723 cDNA positions: without the flag 1,633 have an end that is not their
    // start, with it 1,525 of those have the start for an end, and the 90 single positions are
    // the string `"5/1403"`.
    if config.compat.total_length && config.output_format == "json" {
        eprintln!(
            "warning: --total-length applies to tab and VCF output only; \
             --output-format json writes its positions as numbers, without it."
        );
    }

    // Create consequence predictor
    let predictor = ConsequencePredictor::new(config.distance, config.distance);

    // Open input VCF (supports plain text or gzipped VCF)
    let input_reader = open_vcf_input_reader(&config.input)?;
    // A record that cannot be read is skipped with a message, as Ensembl VEP does, and counted at
    // the end: one malformed line of a whole-genome VCF is not a reason to lose the run.
    let mut vcf_parser = VcfParser::new(input_reader)?.lenient();

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

    // Resolved once, next to the order it will be paired with: the CSQ field
    // list, the tab header and the per-variant pick all have to agree about
    // whether this run flags, and deriving each from `config.pick`
    // independently is how they would come to disagree.
    let pick_request = config.pick.resolve();
    // `&& !sa_only` because the run has already told the user it is ignoring
    // the switch (see the warning above), and `annotate_variant` skips the
    // pick entirely in that mode. Declaring the column anyway would put a
    // `PICK` header over a column that is empty on every row, which reads as
    // "every entry lost the pick" rather than "no pick ran".
    let pick_flags_column = !sa_only && pick_request.is_some_and(|r| r.needs_pick_column());
    let csq_fields: Vec<&'static str> = output::csq_fields(pick_flags_column, config.columns);
    // The optional tab columns, in VEP's order. `--sa-only` writes no transcript
    // columns at all, and says so for each flag it ignores.
    let optional_columns: Vec<output::OptionalColumn> = if sa_only {
        Vec::new()
    } else {
        config.columns.columns(config.hgvs)
    };
    let csq_options = output::CsqOptions {
        no_escape: config.compat.no_escape,
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
    // Precompute the loaded-source lookup once so the per-row writers don't
    // redo an O(specs × keys) membership scan for every variant. The owned-ID
    // list and the header lines are derived from the same lookup, so the
    // headers a run declares and the fields it emits cannot disagree.
    let supplementary_specs = output::LoadedSupplementarySpecs::new(&sa_json_keys, &gene_json_keys);
    let owned_vcf_info_ids = output::vcf_owned_info_ids(&supplementary_specs);
    let generated_vcf_headers =
        output::vcf_info_header_lines(&supplementary_specs, &csq_fields, sa_only);

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
            let header = output::tab_header(
                sa_only,
                config.explicit_alleles,
                &optional_columns,
                pick_flags_column,
                &output::tab_supplementary_column_names(&supplementary_specs),
                qc_rules.is_some(),
            );
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
        let mut batch: Vec<(VariationFeature, KnownLookup)> = Vec::with_capacity(BATCH_SIZE);
        for _ in 0..BATCH_SIZE {
            match vcf_parser.next_variant()? {
                Some(vf) => {
                    if let Some(check) = contig_check.as_mut() {
                        check.note(&vf.position.chromosome);
                    }
                    // Variation lookup (sequential - it reads the cache's open files).
                    //
                    // A failed lookup stops the run. It used to be read as "no known
                    // variants", and a variant reported as novel because the cache
                    // could not be read cannot be told from one that is novel. The
                    // SA lookups are counted and reported at the end instead, but they
                    // run inside the parallel loop where an error cannot propagate;
                    // this one can, and the cache has already been validated at
                    // startup, so an error here is a real fault.
                    let known: KnownLookup = if let Some(ref vp) = var_provider {
                        known_variants(vp, &vf)?
                    } else {
                        KnownLookup::default()
                    };

                    batch.push((vf, known));
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

        let ctx = AnnotationContext {
            config: &config,
            sa_only,
            transcript_provider: &transcript_provider,
            seq_provider: seq_provider.as_deref(),
            predictor: &predictor,
            sa_providers: &sa_providers,
            gene_providers: &gene_providers,
            acmg_config: acmg_config.as_ref(),
            functional_evidence,
            pick_order: &pick_order,
            pick: pick_request,
            sample_names: &sample_names,
            normalize_gnomad_queries,
            repeat_db_loaded,
        };
        // Phase 2: Annotate batch in parallel (transcript lookup + consequence prediction + HGVS)
        batch
            .par_iter_mut()
            .for_each(|(vf, known)| annotate_variant(&ctx, vf, known));

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
                "vcf" => write_vcf_line(
                    &mut writer,
                    vf,
                    sa_only,
                    &supplementary_specs,
                    &csq_fields,
                    csq_options,
                )?,
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
                        pick_column: pick_flags_column,
                        optional: &optional_columns,
                        no_escape: config.compat.no_escape,
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
    // After the meter, so a run that lost annotations says so last and the
    // warning is not buried in progress output.
    fastvep_annotate::report_sa_lookup_errors();
    fastvep_annotate::report_ref_mismatches();
    if let Some(note) = var_provider.as_ref().and_then(|vp| vp.finish()) {
        eprintln!("{note}");
    }
    if vcf_parser.skipped() > 0 {
        eprintln!(
            "warning: {} record(s) of the VCF could not be read and were not annotated (the first are named above)",
            vcf_parser.skipped()
        );
    }
    if let Some(check) = &contig_check {
        check.finish(!sa_providers.is_empty() || var_provider.is_some())?;
    }

    Ok(())
}
/// Say so when the cache has none of the columns a flag asks for.
///
/// VEP answers with empty columns and says nothing, and so does this, because
/// that is what a script that works against VEP expects; but an empty `gnomADg_AF`
/// column over a cache that never had the data reads as "no variant is in gnomAD",
/// so the run says it once, here, where the cache is opened.
fn warn_about_missing_columns(known: fastvep_io::known::KnownFlags, variation_cols: &[String]) {
    let has = |name: &str| variation_cols.iter().any(|c| c == name);
    for group in &fastvep_io::known::FREQUENCY_GROUPS {
        if known.group(group.flag) && !group.keys.iter().any(|k| has(k.key)) {
            eprintln!(
                "warning: --{} asks for {} frequencies, and this cache's variation_cols has none \
                 of their columns ({}): they will be empty for every variant, as VEP writes them.",
                group.flag.replace('_', "-"),
                group.description,
                group
                    .keys
                    .iter()
                    .map(|k| k.key)
                    .collect::<Vec<_>>()
                    .join(", ")
            );
        }
    }
    if known.pubmed && !has("pubmed") {
        eprintln!(
            "warning: --pubmed asks for PUBMED, and this cache's variation_cols has no `pubmed` \
             column: it will be empty for every variant, as VEP writes it."
        );
    }
}

/// What the variation cache attached to one input record.
#[derive(Default)]
struct KnownLookup {
    /// The alt alleles the lookup was made with, symbolic ones and `*` left out:
    /// `AlleleMatch::alt_index` of the site counts into this list.
    alleles: Vec<Allele>,
    /// The known variants of each alt allele as the rows write it, in the order
    /// VEP prints `Existing_variation`.
    by_allele: HashMap<String, Vec<MatchedVariant>>,
    /// Every known variant attached to the record, in cache order, with the
    /// alleles each matched: what VEP's frequency columns are computed from.
    site: KnownSite,
}

/// The known variants VEP's `--check_existing` reports for one input record.
///
/// VEP leaves structural variants out of the check
/// (`VariationTabix::annotate_InputBuffer` keeps only plain
/// `VariationFeature`s), and a symbolic or missing allele has no sequence to
/// compare, so neither is looked up. Without the first rule a large `END=` would
/// read every cache record under the span.
fn known_variants(provider: &TabixVariationProvider, vf: &VariationFeature) -> Result<KnownLookup> {
    if vf.variant_type.is_structural() {
        return Ok(KnownLookup::default());
    }
    let alleles: Vec<Allele> = vf
        .alt_alleles
        .iter()
        .filter(|a| !matches!(a, Allele::Missing | Allele::Symbolic(_)))
        .cloned()
        .collect();
    if alleles.is_empty() {
        return Ok(KnownLookup::default());
    }
    let alts: Vec<String> = alleles.iter().map(ToString::to_string).collect();
    let site = provider.known_site(
        &vf.position.chromosome,
        vf.position.start,
        vf.position.end,
        &vf.ref_allele.to_string(),
        &alts,
    )?;
    let by_allele = alts
        .into_iter()
        .zip(site.per_alt(alleles.len()))
        .filter(|(_, matches)| !matches.is_empty())
        .collect();
    Ok(KnownLookup {
        alleles,
        by_allele,
        site,
    })
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
                aa.splice,
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

fn write_vcf_line(
    writer: &mut impl Write,
    vf: &VariationFeature,
    sa_only: bool,
    specs: &output::LoadedSupplementarySpecs,
    csq_fields: &[&str],
    csq_options: output::CsqOptions,
) -> Result<()> {
    if let Some(ref fields) = vf.vcf_fields {
        let csq = if sa_only {
            String::new()
        } else {
            output::format_csq_with(vf, csq_fields, csq_options)
        };
        let info = output::format_vcf_info_fields(&fields.info, vf, &csq, specs);

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
