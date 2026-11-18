//! Consequences of a structural variant, as Ensembl VEP assigns them.
//!
//! VEP has no model of what a structural variant does to a transcript. It has the class the record
//! was given ([`SvClass`]), the variant's interval, and a list of consequence definitions
//! (`%OVERLAP_CONSEQUENCES`, `Bio::EnsEMBL::Variation::Utils::Constants`) each with a predicate in
//! `Utils/VariationEffect.pm`. For each (transcript, allele) it asks the definitions in order of
//! their tier, skips those whose `include` conditions the pair does not meet, and keeps the ones
//! whose predicate holds (`BaseVariationFeatureOverlapAllele::get_all_OverlapConsequences`). This
//! module is that loop for the definitions a structural variant can reach, one function to each
//! predicate and named for it.
//!
//! What makes it a port and not a model is what it leaves as VEP has it:
//!
//! - **Tiers.** `transcript_ablation` and `transcript_amplification` are tier 1, and a tier-1 term
//!   ends the loop: a deletion that covers a transcript is that and nothing else.
//! - **The flags of `_bvfo_preds`.** A variant that covers the whole transcript is given only
//!   `complete_overlap`, `within_feature` and the biotype, no `exon`, `intron`, `utr` or `coding`,
//!   so none of the terms gated on them can be said of it, and it gets no positions. Every other
//!   overlap gets the full set, and `coding` is `exon` for a transcript the variant reaches the
//!   coding span of.
//! - **What a variant cannot be.** The terms about the protein that read the two alleles
//!   (`splice_*`, `stop_gained`, `missense_variant`, `synonymous_variant`, ...) are for small
//!   variants, `inframe_insertion` is "TO BE DONE" and always false, and `frameshift_variant`,
//!   `inframe_deletion`, `stop_lost`, `start_lost` and `start_retained_variant` have branches that
//!   read the interval alone.
//! - **The default.** A pair no definition holds for is `intergenic_variant`, on the transcript's
//!   row.
//!
//! `docs/STRUCTURAL_VARIANTS.md` lists the rules, and what VEP does in the cases that look like
//! defects.

use std::cell::OnceCell;
use std::sync::Arc;

use fastvep_core::{
    chrom::same_contig, Allele, AlleleNo, AlleleSite, Consequence, Impact, Strand,
    StructuralVariant, SvAllele, SvClass, TranscriptLocation,
};
use fastvep_genome::Transcript;

use crate::predictor::{AlleleConsequenceResult, TranscriptConsequence};
use crate::splice::{span_effects, SpanEffects, SpliceEffects};

/// How far from a transcript the position of an allele may be for VEP to make a row of it:
/// `MAX_DISTANCE_FROM_TRANSCRIPT`, a constant, which `--distance` does not change.
const CLOSE_TO_FEATURE: u64 = 5000;

/// The contigs a structural variant touches, by the names its record gives them: its own, then each
/// mate's, once.
///
/// A transcript is found by asking its provider for a contig by one of these names, so which of them
/// it is on is known from the question that found it, and no name of the transcript has to be
/// compared with one of the record's (which a chromosome synonyms file would make wrong).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Contigs(Vec<Arc<str>>);

impl Contigs {
    pub fn of(own: &str, sv: &StructuralVariant) -> Self {
        let mut names: Vec<Arc<str>> = vec![Arc::from(own)];
        for (name, _) in sv.mates() {
            if !names.iter().any(|n| same_contig(n, name)) {
                names.push(Arc::from(name));
            }
        }
        Self(names)
    }

    /// The number of the contig a mate is on.
    pub fn index_of(&self, name: &str) -> usize {
        self.0
            .iter()
            .position(|n| same_contig(n, name))
            .unwrap_or(0)
    }
}

/// A transcript a structural variant is paired with, and the contig of [`Contigs`] it was found on.
#[derive(Debug, Clone, Copy)]
pub struct SvTranscript<'a> {
    pub transcript: &'a Transcript,
    pub contig: usize,
}

/// What VEP says of one transcript, with the transcript it said it of.
#[derive(Debug, Clone)]
pub struct SvPrediction<'a> {
    pub consequence: TranscriptConsequence,
    pub transcript: &'a Transcript,
}

/// The consequences of a structural variant for each transcript it is paired with, in the order VEP
/// writes them (by the transcript's stable id).
///
/// `start` and `end` are the variant's interval, with `start = end + 1` for an insertion.
/// `alt_alleles` are the ALTs of the record, which the rows answer for by [`SvAllele::alt`], and
/// `distance` is `--distance`. A transcript none of whose alleles is close to it has no row, as in
/// VEP: it is paired and says nothing.
pub fn predict_sv<'a>(
    sv: &StructuralVariant,
    alt_alleles: &[Allele],
    start: u64,
    end: u64,
    contigs: &Contigs,
    transcripts: &[SvTranscript<'a>],
    distance: u64,
) -> Vec<SvPrediction<'a>> {
    let mut out: Vec<SvPrediction<'a>> = Vec::with_capacity(transcripts.len());
    for found in transcripts {
        let tr = found.transcript;
        let alleles: Vec<(usize, &SvAllele)> = sv
            .alleles
            .iter()
            .enumerate()
            .filter(|(_, a)| is_close_to_feature(a, found, contigs, start, end))
            .collect();
        if alleles.is_empty() {
            continue;
        }
        let ctx = Ctx::new(tr, sv.class, start, end, distance);
        let location = ctx.cells();
        let allele_consequences = alleles
            .into_iter()
            .map(|(number, allele)| {
                let breakend_in_feature = ctx.breakend_in_feature(allele, found, contigs);
                let consequences = ctx.consequences(breakend_in_feature);
                let impact = consequences
                    .iter()
                    .map(|c| c.impact())
                    .min()
                    .unwrap_or(Impact::Modifier);
                let beside = consequences.iter().any(|c| {
                    matches!(
                        c,
                        Consequence::UpstreamGeneVariant | Consequence::DownstreamGeneVariant
                    )
                });
                AlleleConsequenceResult {
                    allele: alt_alleles
                        .get(allele.alt)
                        .cloned()
                        .unwrap_or(Allele::Missing),
                    sv_allele: AlleleNo::of(number),
                    consequences,
                    impact,
                    cdna_start: None,
                    cdna_end: None,
                    cds_start: None,
                    cds_end: None,
                    protein_start: None,
                    protein_end: None,
                    amino_acids: None,
                    codons: None,
                    exon: None,
                    intron: None,
                    location,
                    distance: beside.then(|| ctx.distance_to_transcript()),
                    protein_length: tr.peptide_length(),
                    escapes_nmd: None,
                    coding_sequence_kept: false,
                }
            })
            .collect();
        out.push(SvPrediction {
            consequence: TranscriptConsequence {
                transcript_id: Arc::clone(&tr.stable_id),
                gene_id: tr.gene.stable_id.clone(),
                gene_symbol: tr.gene.symbol.clone(),
                biotype: Arc::clone(&tr.biotype),
                allele_consequences,
                canonical: tr.canonical,
                strand: tr.strand,
            },
            transcript: tr,
        });
    }
    // `get_all_TranscriptStructuralVariations` returns `sort keys` of the transcripts' stable ids
    out.sort_by(|a, b| {
        a.consequence
            .transcript_id
            .cmp(&b.consequence.transcript_id)
    });
    out
}

/// `_close_to_feature` (`StructuralVariationOverlap.pm`), which decides which alleles a breakend gets
/// a row for: those on the transcript's contig and within [`CLOSE_TO_FEATURE`] bases of it. An allele
/// with no position of its own is always in.
fn is_close_to_feature(
    allele: &SvAllele,
    found: &SvTranscript<'_>,
    contigs: &Contigs,
    start: u64,
    end: u64,
) -> bool {
    let tr = found.transcript;
    let near = |s: u64, e: u64| {
        // `overlap($start, $end, $slice->start, $slice->end)` of the transcript expanded
        e.saturating_add(CLOSE_TO_FEATURE) >= tr.start
            && s <= tr.end.saturating_add(CLOSE_TO_FEATURE)
    };
    match &allele.site {
        AlleleSite::None => true,
        AlleleSite::Own => found.contig == 0 && near(start, end),
        AlleleSite::Mate {
            chromosome,
            position,
        } => found.contig == contigs.index_of(chromosome) && near(*position, *position),
    }
}

/// The flags `_bvfo_preds` leaves for a variant that reaches a transcript without covering it.
struct Pre {
    /// `_overlapped_exons` is not empty ([`Transcript::overlaps_exon_flag`])
    exon: bool,
    /// `_overlapped_introns` is not empty
    intron: bool,
    /// the variant reaches the UTR: it is wholly outside the coding region or not wholly inside it
    utr: bool,
    /// a transcript with a coding region, which the variant reaches, and an exon with it
    coding: bool,
}

/// One structural variant against one transcript: the facts every predicate reads, worked out once.
struct Ctx<'a> {
    tr: &'a Transcript,
    class: SvClass,
    s: u64,
    e: u64,
    distance: u64,
    /// `complete_overlap_feature`
    complete: bool,
    /// `within_feature`: a plain overlap, so an insertion needs both its flanks inside
    wf: bool,
    span: OnceCell<SpanEffects>,
    location: OnceCell<TranscriptLocation>,
}

impl<'a> Ctx<'a> {
    fn new(tr: &'a Transcript, class: SvClass, s: u64, e: u64, distance: u64) -> Self {
        Self {
            tr,
            class,
            s,
            e,
            distance,
            complete: s <= tr.start && e >= tr.end,
            wf: e >= tr.start && s <= tr.end,
            span: OnceCell::new(),
            location: OnceCell::new(),
        }
    }

    fn lo(&self) -> u64 {
        self.s.min(self.e)
    }

    fn hi(&self) -> u64 {
        self.s.max(self.e)
    }

    /// Ensembl's insertion: `start = end + 1`
    fn is_insertion(&self) -> bool {
        self.s == self.e + 1
    }

    fn length(&self) -> u64 {
        self.e + 1 - self.s
    }

    /// `_intron_effects` of the whole span
    fn span(&self) -> &SpanEffects {
        self.span
            .get_or_init(|| span_effects(self.tr, self.s, self.e))
    }

    fn splice(&self) -> &SpliceEffects {
        &self.span().splice
    }

    /// Where the span sits in the transcript, as the exon mapper reports it
    fn mapped(&self) -> &TranscriptLocation {
        self.location
            .get_or_init(|| self.tr.location(self.s, self.e))
    }

    fn coding_region(&self) -> Option<(u64, u64)> {
        Some((self.tr.coding_region_start?, self.tr.coding_region_end?))
    }

    fn pre(&self) -> Pre {
        let exon = self.tr.overlaps_exon_flag(self.s, self.e);
        let mut pre = Pre {
            exon,
            intron: self.span().near_intron,
            utr: false,
            coding: false,
        };
        if let Some((cs, ce)) = self.coding_region() {
            if self.hi() < cs || self.lo() > ce {
                pre.utr = true;
            } else {
                pre.utr = self.lo() < cs || self.hi() > ce;
                pre.coding = exon;
            }
        }
        pre
    }

    // ---- the predicates of VariationEffect.pm ------------------------------------------------

    /// `within_cdna`: the mapper has a coordinate for the variant that lies in the transcript's
    /// cDNA (`end > 0` and `start <= length`), or it is inside a frameshift intron.
    ///
    /// An insertion is the pair `(c + 1, c)` that `Mapper::map_insert` leaves for the one flank that
    /// is in an exon, so it is not in the cDNA when that flank is the transcript's own first or last
    /// base and the insertion is on the outside of it: `(1, 0)` has no `end`, `(L + 1, L)` has a
    /// `start` past the length. Ensembl VEP 116.3 says only `upstream_gene_variant` of such an
    /// insertion, where the exon test alone would add `feature_elongation`.
    fn within_cdna(&self) -> bool {
        let mapped = if self.is_insertion() {
            let cdna = self.mapped().cdna;
            matches!((cdna.start, cdna.end),
                (Some(start), Some(end)) if end > 0 && start <= self.tr.cdna_length())
        } else {
            // some exon has a base of the variant, and the mapper has a coordinate for that piece
            self.tr
                .exons
                .iter()
                .any(|x| self.s <= x.end && self.e >= x.start)
        };
        mapped || (self.wf && self.splice().within_frameshift_intron)
    }

    /// `within_cds`
    fn within_cds(&self) -> bool {
        let Some((cs, ce)) = self.coding_region().filter(|_| self.tr.is_coding()) else {
            return false;
        };
        let in_cds = if self.is_insertion() {
            // the pair of bases it sits between maps to a CDS coordinate or does not
            !self.mapped().cds.is_none()
        } else {
            let (lo, hi) = (self.lo(), self.hi());
            self.tr
                .exons
                .iter()
                .any(|x| lo.max(x.start).max(cs) <= hi.min(x.end).min(ce))
        };
        in_cds || (self.splice().within_frameshift_intron && self.e >= cs && self.s <= ce)
    }

    /// some exon contains the whole variant (`complete_within_feature` of an exon)
    fn in_one_exon(&self) -> bool {
        self.tr
            .exons
            .iter()
            .any(|x| self.s >= x.start && self.e <= x.end)
    }

    /// `frameshift`, the branch for a structural variant: a deletion inside one exon, whose length
    /// is not a multiple of three. It does not ask whether the exon is coding.
    fn frameshift(&self) -> bool {
        self.class.is_deletion() && self.in_one_exon() && !self.length().is_multiple_of(3)
    }

    /// `inframe_deletion`: a deletion inside the coding part of one exon (the mapper gives it one
    /// CDS coordinate and no gap), a multiple of three long
    fn inframe_deletion(&self) -> bool {
        let Some((cs, ce)) = self.coding_region() else {
            return false;
        };
        self.class.is_deletion()
            && self
                .tr
                .exons
                .iter()
                .any(|x| self.s >= x.start.max(cs) && self.e <= x.end.min(ce))
            && self.length().is_multiple_of(3)
    }

    /// `partial_codon`: the variant starts in the last codon of a coding sequence that ends in a
    /// part of one
    fn partial_codon(&self) -> bool {
        let (Some(p), Some(len)) = (self.mapped().protein.start, self.tr.cds_length()) else {
            return false;
        };
        let last_codon = len as i64 - (p as i64 * 3 - 3);
        (1..3).contains(&last_codon)
    }

    /// `_overlaps_start_codon`: the cDNA positions of the variant's two ends meet the first three
    /// bases of the coding sequence. It is false when either end is outside an exon (the mapper has
    /// no coordinate for it) and for a transcript whose start is not known (`cds_start_NF`).
    fn overlaps_start_codon(&self) -> bool {
        if self.tr.flags.iter().any(|f| f == "cds_start_NF") {
            return false;
        }
        let cdna = self.mapped().cdna;
        let (Some(first), Some(last), Some(coding_start)) =
            (cdna.start, cdna.end, self.tr.cdna_coding_start)
        else {
            return false;
        };
        first != 0 && last != 0 && last >= coding_start && first <= coding_start + 2
    }

    /// `start_lost`: it overlaps the start codon by the mapper, and the genomic bases of that codon
    /// by the interval
    fn start_lost(&self) -> bool {
        let Some((cs, ce)) = self.coding_region() else {
            return false;
        };
        self.overlaps_start_codon()
            && match self.tr.strand {
                Strand::Forward => cs + 2 >= self.s && cs <= self.e,
                Strand::Reverse => ce >= self.s && ce.saturating_sub(2) <= self.e,
            }
    }

    /// `stop_lost`: a deletion over the genomic bases of the stop codon, unless it starts in a last
    /// codon that is partial
    fn stop_lost(&self) -> bool {
        let Some((cs, ce)) = self.coding_region() else {
            return false;
        };
        self.class.is_deletion()
            && !self.partial_codon()
            && match self.tr.strand {
                Strand::Forward => ce >= self.s && ce.saturating_sub(2) <= self.e,
                Strand::Reverse => cs + 2 >= self.s && cs <= self.e,
            }
    }

    /// `_before_coding`, which is the 5' side on the forward strand and the 3' side on the reverse
    fn before_coding(&self) -> bool {
        let Some((cs, _)) = self.coding_region().filter(|_| self.tr.is_coding()) else {
            return false;
        };
        if self.is_insertion() && self.s == cs {
            return true;
        }
        self.e >= self.tr.start && cs > 0 && self.s < cs
    }

    /// `_after_coding`
    fn after_coding(&self) -> bool {
        let Some((_, ce)) = self.coding_region().filter(|_| self.tr.is_coding()) else {
            return false;
        };
        if self.is_insertion() && self.e == ce {
            return true;
        }
        self.e > ce && self.s <= self.tr.end
    }

    /// `non_coding_exon_variant`: not a complete overlap, and a plain overlap with an exon of a
    /// transcript that has no protein
    fn non_coding_exon_variant(&self) -> bool {
        !self.complete
            && !self.tr.is_coding()
            && self
                .tr
                .exons
                .iter()
                .any(|x| self.e >= x.start && self.s <= x.end)
    }

    /// `_before_start`
    fn before_start(&self, distance: u64) -> bool {
        self.e.saturating_add(distance) >= self.tr.start && self.e < self.tr.start
    }

    /// `_after_end`
    fn after_end(&self, distance: u64) -> bool {
        self.s <= self.tr.end.saturating_add(distance) && self.s > self.tr.end
    }

    /// `feature_truncation` of a breakpoint: the allele's own position is inside the transcript.
    ///
    /// A breakend that has no mate (`A.`, `.A`, a bare `<BND>`) has no position of its own to ask
    /// about, and `within_feature` then asks about the variant's, which is the same place. Ensembl VEP
    /// writes the term for it in a normal run and not in a run with `--no_stats`, where the allele has
    /// had `{breakend}` created on it by `OutputFactory.pm` before the term is worked out (the hash is
    /// empty and has no chromosome, so nothing is inside anything); the first is the term the code is
    /// written to give, so it is the one here.
    fn breakend_in_feature(
        &self,
        allele: &SvAllele,
        found: &SvTranscript<'_>,
        contigs: &Contigs,
    ) -> bool {
        match &allele.site {
            AlleleSite::None | AlleleSite::Own => {
                found.contig == 0 && self.e >= self.tr.start && self.s <= self.tr.end
            }
            AlleleSite::Mate {
                chromosome,
                position,
            } => {
                found.contig == contigs.index_of(chromosome)
                    && *position >= self.tr.start
                    && *position <= self.tr.end
            }
        }
    }

    // ---- the loop of get_all_OverlapConsequences ----------------------------------------------

    /// The terms of this variant for this transcript, in rank order. `breakend_in_feature` is the
    /// one thing that depends on the allele, for a breakpoint.
    fn consequences(&self, breakend_in_feature: bool) -> Vec<Consequence> {
        let class = self.class;

        // tier 1: a complete overlap by a deletion or a gain is that and nothing else
        if self.complete {
            if class.is_deletion() {
                return vec![Consequence::TranscriptAblation];
            }
            if class.is_copy_number_gain() {
                return vec![Consequence::TranscriptAmplification];
            }
        }

        let tr = self.tr;
        let mut out: Vec<Consequence> = Vec::with_capacity(6);
        let protein_coding = &*tr.biotype == "protein_coding";

        if self.wf {
            if !self.complete {
                let pre = self.pre();
                if pre.coding {
                    let frameshift = self.frameshift();
                    let inframe_deletion = self.inframe_deletion();
                    if frameshift {
                        out.push(Consequence::FrameshiftVariant);
                    }
                    if self.stop_lost() {
                        out.push(Consequence::StopLost);
                    }
                    if self.start_lost() {
                        out.push(Consequence::StartLost);
                    }
                    if inframe_deletion {
                        out.push(Consequence::InframeDeletion);
                    }
                    if self.overlaps_start_codon() {
                        out.push(Consequence::StartRetainedVariant);
                    }
                    // `coding_unknown`
                    if self.within_cds() && !(inframe_deletion || frameshift) {
                        out.push(Consequence::CodingSequenceVariant);
                    }
                }
                if pre.exon && pre.utr && self.within_cdna() {
                    if tr.strand == Strand::Forward && self.before_coding()
                        || tr.strand == Strand::Reverse && self.after_coding()
                    {
                        out.push(Consequence::FivePrimeUtrVariant);
                    }
                    if tr.strand == Strand::Forward && self.after_coding()
                        || tr.strand == Strand::Reverse && self.before_coding()
                    {
                        out.push(Consequence::ThreePrimeUtrVariant);
                    }
                }
                if pre.exon && !protein_coding && self.non_coding_exon_variant() {
                    out.push(Consequence::NonCodingTranscriptExonVariant);
                }
                if pre.intron {
                    let splice = self.splice();
                    if splice.intronic {
                        out.push(Consequence::IntronVariant);
                    }
                    if splice.polypyrimidine_tract {
                        out.push(Consequence::SplicePolypyrimidineTractVariant);
                    }
                }
            }
            if &*tr.biotype == "nonsense_mediated_decay" {
                out.push(Consequence::NmdTranscriptVariant);
            }
            if !protein_coding && !tr.is_coding() && !self.non_coding_exon_variant() {
                out.push(Consequence::NonCodingTranscriptVariant);
            }
            if protein_coding && self.complete && tr.is_coding() {
                out.push(Consequence::CodingTranscriptVariant);
            }
        }

        // `feature_elongation`: a gain or an insertion inside the transcript that touches an exon
        if class.is_insertion() && self.within_cdna() && self.s >= tr.start && self.e <= tr.end {
            out.push(Consequence::FeatureElongation);
        }

        // `feature_truncation`
        if class.is_breakpoint() {
            if breakend_in_feature {
                out.push(Consequence::FeatureTruncation);
            }
        } else if class.is_deletion()
            && self.within_cdna()
            && ((self.wf && !self.complete && (self.e > tr.end || self.s < tr.start))
                || (self.s >= tr.start && self.e <= tr.end))
        {
            out.push(Consequence::FeatureTruncation);
        }

        if !self.wf {
            let (up, down) = (self.distance, self.distance);
            let (upstream, downstream) = match tr.strand {
                Strand::Forward => (self.before_start(up), self.after_end(down)),
                Strand::Reverse => (self.after_end(up), self.before_start(down)),
            };
            if upstream {
                out.push(Consequence::UpstreamGeneVariant);
            }
            if downstream {
                out.push(Consequence::DownstreamGeneVariant);
            }
        }

        // the definition no predicate holds for is the default one
        if out.is_empty() {
            out.push(Consequence::IntergenicVariant);
        }
        out.sort_by_key(|c| c.rank());
        out
    }

    /// `distance_to_transcript`
    fn distance_to_transcript(&self) -> i64 {
        let d = |a: u64, b: u64| a.abs_diff(b) as i64;
        d(self.s, self.tr.start)
            .min(d(self.s, self.tr.end))
            .min(d(self.e, self.tr.start))
            .min(d(self.e, self.tr.end))
    }

    /// The `cDNA_position`, `CDS_position`, `Protein_position`, `EXON` and `INTRON` cells, as
    /// `TranscriptStructuralVariationAllele_to_output_hash` writes them: positions when the variant
    /// is `exon` (and CDS and protein when it is `coding`), the exon numbers then too, the intron
    /// numbers when it is `intron`, and nothing for a variant that covers the transcript.
    fn cells(&self) -> TranscriptLocation {
        if self.complete || !self.wf {
            return TranscriptLocation::default();
        }
        let pre = self.pre();
        let mapped = self.mapped();
        let mut cells = TranscriptLocation::default();
        if pre.exon {
            cells.cdna = mapped.cdna;
            cells.exon = mapped.exon;
            if pre.coding {
                cells.cds = mapped.cds;
                cells.protein = mapped.protein;
            }
        }
        if pre.intron {
            cells.intron = mapped.intron;
        }
        cells
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fastvep_core::SvClass::*;
    use fastvep_genome::{Exon, Gene, Translation};

    /// A transcript on `chr1`: exons `(start, end)` in genomic order and the coding span.
    fn transcript(
        strand: Strand,
        biotype: &str,
        exons: &[(u64, u64)],
        cds: Option<(u64, u64)>,
    ) -> Transcript {
        let n = exons.len();
        let in_order: Vec<(u64, u64)> = if strand == Strand::Forward {
            exons.to_vec()
        } else {
            exons.iter().rev().copied().collect()
        };
        let mut t = Transcript {
            stable_id: Arc::from("ENST_T"),
            version: Some(1),
            gene: Gene {
                stable_id: Arc::from("ENSG_T"),
                symbol: None,
                symbol_source: None,
                hgnc_id: None,
                biotype: Arc::from(biotype),
                chromosome: Arc::from("chr1"),
                start: exons[0].0,
                end: exons[n - 1].1,
                strand,
            },
            biotype: Arc::from(biotype),
            chromosome: Arc::from("chr1"),
            start: exons[0].0,
            end: exons[n - 1].1,
            strand,
            exons: in_order
                .iter()
                .enumerate()
                .map(|(i, &(start, end))| Exon {
                    stable_id: format!("E{i}"),
                    start,
                    end,
                    strand,
                    phase: 0,
                    end_phase: 0,
                    rank: i as u32 + 1,
                })
                .collect(),
            translation: None,
            cdna_coding_start: None,
            cdna_coding_end: None,
            coding_region_start: None,
            coding_region_end: None,
            spliced_seq: None,
            translateable_seq: None,
            peptide: None,
            canonical: true,
            mane_select: None,
            mane_plus_clinical: None,
            tsl: None,
            appris: None,
            ccds: None,
            protein_id: None,
            protein_version: None,
            swissprot: vec![],
            trembl: vec![],
            uniparc: vec![],
            refseq_id: None,
            source: None,
            gencode_primary: false,
            flags: vec![],
            codon_table_start_phase: 0,
            initiator_met: false,
            cds_flags_from_sequence: false,
            readthrough: false,
        };
        if let Some((s, e)) = cds {
            let (first, last) = if strand == Strand::Forward {
                (s, e)
            } else {
                (e, s)
            };
            t.cdna_coding_start = t.genomic_to_cdna(first);
            t.cdna_coding_end = t.genomic_to_cdna(last);
            t.coding_region_start = Some(s);
            t.coding_region_end = Some(e);
            t.translation = Some(Translation {
                stable_id: "ENSP_T".into(),
                genomic_start: s,
                genomic_end: e,
                start_exon_rank: 1,
                start_exon_offset: 0,
                end_exon_rank: n as u32,
                end_exon_offset: 0,
            });
        }
        t
    }

    /// The transcript `A5` of the world the real VEP 116.3 was run over to write the expected
    /// terms below: forward strand, three exons, coding from 2009 to 2135.
    fn a5() -> Transcript {
        transcript(
            Strand::Forward,
            "protein_coding",
            &[(2000, 2039), (2060, 2089), (2110, 2139)],
            Some((2009, 2135)),
        )
    }

    /// `A6`: reverse strand, three exons, coding from 2404 to 2530.
    fn a6() -> Transcript {
        transcript(
            Strand::Reverse,
            "protein_coding",
            &[(2400, 2439), (2460, 2489), (2510, 2539)],
            Some((2404, 2530)),
        )
    }

    /// `A3`, a non-coding transcript of one exon.
    fn a3() -> Transcript {
        transcript(Strand::Forward, "lncRNA", &[(1200, 1239)], None)
    }

    fn terms(tr: &Transcript, class: SvClass, s: u64, e: u64) -> Vec<&'static str> {
        let ctx = Ctx::new(tr, class, s, e, 1000);
        ctx.consequences(false)
            .into_iter()
            .map(Consequence::so_term)
            .collect()
    }

    #[test]
    fn a_deletion_over_a_transcript_is_its_ablation_and_nothing_else() {
        assert_eq!(terms(&a5(), Deletion, 1991, 2150), ["transcript_ablation"]);
        assert_eq!(terms(&a3(), Deletion, 1191, 1250), ["transcript_ablation"]);
        // exactly the transcript is a complete overlap too
        assert_eq!(terms(&a5(), Deletion, 2000, 2139), ["transcript_ablation"]);
        // a base short of it is not
        assert_ne!(terms(&a5(), Deletion, 2001, 2139), ["transcript_ablation"]);
    }

    #[test]
    fn a_duplication_over_a_transcript_is_its_amplification() {
        assert_eq!(
            terms(&a5(), Duplication, 1991, 2150),
            ["transcript_amplification"]
        );
        assert_eq!(
            terms(&a3(), TandemDuplication, 1191, 1250),
            ["transcript_amplification"]
        );
        // an insertion is not a copy number gain, and neither is a copy number variation
        assert_eq!(
            terms(&a5(), Insertion, 1991, 2150),
            ["coding_transcript_variant"]
        );
        assert_eq!(
            terms(&a5(), CopyNumberVariation, 1991, 2150),
            ["coding_transcript_variant"]
        );
    }

    #[test]
    fn what_covers_a_transcript_without_ablating_or_amplifying_it_has_only_the_terms_of_the_whole()
    {
        // no exon, intron, UTR or coding flag is set for a complete overlap
        assert_eq!(
            terms(&a5(), Inversion, 1991, 2150),
            ["coding_transcript_variant"]
        );
        assert_eq!(
            terms(&a3(), Inversion, 1191, 1250),
            ["non_coding_transcript_variant"]
        );
        assert_eq!(
            terms(&a3(), CopyNumberVariation, 1191, 1250),
            ["non_coding_transcript_variant"]
        );
    }

    #[test]
    fn a_deletion_that_reaches_a_transcript_truncates_it_and_names_what_it_reaches() {
        // exon 1 to exon 2 over the intron between: coding, and the intron
        assert_eq!(
            terms(&a5(), Deletion, 2031, 2070),
            [
                "feature_truncation",
                "coding_sequence_variant",
                "intron_variant"
            ]
        );
        // wholly in an intron, in the polypyrimidine tract: not a truncation (no exon is touched)
        assert_eq!(
            terms(&a5(), Deletion, 2046, 2055),
            ["splice_polypyrimidine_tract_variant", "intron_variant"]
        );
        // 5' UTR into the start codon
        assert_eq!(
            terms(&a5(), Deletion, 2004, 2015),
            [
                "start_lost",
                "feature_truncation",
                "start_retained_variant",
                "coding_sequence_variant",
                "5_prime_UTR_variant"
            ]
        );
        // over the stop codon and past the end of the transcript
        assert_eq!(
            terms(&a5(), Deletion, 2131, 2140),
            [
                "stop_lost",
                "feature_truncation",
                "coding_sequence_variant",
                "3_prime_UTR_variant"
            ]
        );
        // a deletion that starts upstream of the transcript and ends in its UTR
        assert_eq!(
            terms(&a5(), Deletion, 1996, 2005),
            ["feature_truncation", "5_prime_UTR_variant"]
        );
    }

    #[test]
    fn a_deletion_inside_one_exon_is_in_frame_or_a_frameshift_by_its_length() {
        // 2013..2030 is 18 bases of the coding sequence
        assert_eq!(
            terms(&a5(), Deletion, 2013, 2030),
            ["feature_truncation", "inframe_deletion"]
        );
        // 17 bases
        assert_eq!(
            terms(&a5(), Deletion, 2013, 2029),
            ["frameshift_variant", "feature_truncation"]
        );
        // the other classes of deletion read the same
        assert_eq!(
            terms(&a5(), AluDeletion, 2013, 2030),
            ["feature_truncation", "inframe_deletion"]
        );
        // `frameshift_variant` is gated on the `coding` flag, which a deletion wholly in the UTR does
        // not have, so its length is not asked
        assert_eq!(
            terms(&a5(), Deletion, 2002, 2005),
            ["feature_truncation", "5_prime_UTR_variant"]
        );
    }

    #[test]
    fn a_duplication_or_insertion_inside_a_transcript_elongates_it() {
        assert_eq!(
            terms(&a5(), Duplication, 2013, 2030),
            ["feature_elongation", "coding_sequence_variant"]
        );
        assert_eq!(
            terms(&a5(), Duplication, 2061, 2100),
            [
                "feature_elongation",
                "coding_sequence_variant",
                "intron_variant"
            ]
        );
        // a tandem duplication is one, a copy number variation is not
        assert_eq!(
            terms(&a5(), TandemDuplication, 2013, 2030),
            ["feature_elongation", "coding_sequence_variant"]
        );
        assert_eq!(
            terms(&a5(), CopyNumberVariation, 2013, 2030),
            ["coding_sequence_variant"]
        );
        // an insertion is the interval of the bases it sits between (start = end + 1)
        assert_eq!(
            terms(&a5(), Insertion, 2013, 2012),
            ["feature_elongation", "coding_sequence_variant"]
        );
        // one that runs out of the transcript is not an elongation
        assert_eq!(
            terms(&a5(), AluInsertion, 2013, 2312),
            [
                "coding_sequence_variant",
                "3_prime_UTR_variant",
                "intron_variant"
            ]
        );
        // a duplication of non-coding exon
        assert_eq!(
            terms(&a3(), Duplication, 1211, 1220),
            ["feature_elongation", "non_coding_transcript_exon_variant"]
        );
    }

    #[test]
    fn an_inversion_or_a_copy_number_variation_names_only_where_it_is() {
        assert_eq!(
            terms(&a5(), Inversion, 2013, 2030),
            ["coding_sequence_variant"]
        );
        assert_eq!(
            terms(&a5(), Inversion, 2061, 2100),
            ["coding_sequence_variant", "intron_variant"]
        );
        assert_eq!(
            terms(&a3(), Inversion, 1211, 1220),
            ["non_coding_transcript_exon_variant"]
        );
    }

    #[test]
    fn the_reverse_strand_is_read_in_transcript_order() {
        // from the 3' UTR over the stop codon of A6 (stop codon at 2404..2406)
        assert_eq!(
            terms(&a6(), Deletion, 2401, 2410),
            [
                "frameshift_variant",
                "stop_lost",
                "feature_truncation",
                "3_prime_UTR_variant"
            ]
        );
        // 5' end, into the coding sequence: exon 1 ends in the 5' UTR at 2539
        assert_eq!(
            terms(&a6(), Deletion, 2501, 2515),
            [
                "feature_truncation",
                "coding_sequence_variant",
                "intron_variant"
            ]
        );
    }

    #[test]
    fn a_variant_that_does_not_reach_a_transcript_is_upstream_or_downstream_of_it() {
        let t = a5();
        // the forward strand: before the start is upstream, after the end is downstream
        assert_eq!(terms(&t, Deletion, 1701, 1750), ["upstream_gene_variant"]);
        assert_eq!(terms(&t, Deletion, 2151, 2200), ["downstream_gene_variant"]);
        // beyond the distance it is the default term
        assert_eq!(terms(&t, Deletion, 100, 150), ["intergenic_variant"]);
        // the reverse strand swaps them
        let r = a6();
        assert_eq!(terms(&r, Deletion, 2201, 2250), ["downstream_gene_variant"]);
        assert_eq!(terms(&r, Deletion, 2551, 2600), ["upstream_gene_variant"]);
    }

    #[test]
    fn a_breakend_truncates_a_transcript_only_where_its_own_allele_is_in_it() {
        let tr = a5();
        let ctx = Ctx::new(&tr, ChromosomeBreakpoint, 2021, 2021, 1000);
        let names = |inside: bool| -> Vec<_> {
            ctx.consequences(inside)
                .into_iter()
                .map(Consequence::so_term)
                .collect()
        };
        assert_eq!(
            names(true),
            ["feature_truncation", "coding_sequence_variant"]
        );
        assert_eq!(names(false), ["coding_sequence_variant"]);
    }

    #[test]
    fn nothing_that_reads_two_alleles_is_said_of_a_structural_variant() {
        for class in [
            Deletion,
            Duplication,
            Insertion,
            Inversion,
            CopyNumberVariation,
        ] {
            for (s, e) in [
                (2031, 2070),
                (2013, 2030),
                (2036, 2041),
                (2058, 2061),
                (2004, 2015),
            ] {
                for t in terms(&a5(), class, s, e) {
                    assert!(
                        !t.starts_with("splice_acceptor")
                            && !t.starts_with("splice_donor")
                            && t != "splice_region_variant"
                            && t != "stop_gained"
                            && t != "missense_variant"
                            && t != "synonymous_variant"
                            && t != "inframe_insertion",
                        "{class:?} {s}-{e}: {t}"
                    );
                }
            }
        }
    }

    #[test]
    fn the_cells_are_written_for_what_reaches_an_exon_and_not_for_what_covers_the_transcript() {
        let t = a5();
        let cells = |s, e| Ctx::new(&t, Deletion, s, e, 1000).cells();
        // a complete overlap has none
        assert!(cells(1991, 2150).cdna.is_none());
        // one that touches an exon has positions, with `?` where an end is not in one
        let c = cells(2031, 2070);
        assert_eq!((c.cdna.start, c.cdna.end), (Some(32), Some(51)));
        assert_eq!((c.cds.start, c.cds.end), (Some(23), Some(42)));
        assert_eq!((c.protein.start, c.protein.end), (Some(8), Some(14)));
        // one wholly in an intron has an intron number and nothing else
        let c = cells(2046, 2055);
        assert!(c.cdna.is_none() && c.exon.is_none());
        assert!(c.intron.is_some());
    }

    #[test]
    fn the_distance_is_the_nearest_of_the_four_ends() {
        let t = a5();
        let ctx = Ctx::new(&t, Deletion, 1701, 1750, 1000);
        assert_eq!(ctx.distance_to_transcript(), 250);
    }

    #[test]
    fn rows_are_in_the_order_of_the_transcripts_stable_ids_and_an_allele_far_away_has_none() {
        let mut a = a5();
        a.stable_id = Arc::from("ENST_B");
        let mut b = a3();
        b.stable_id = Arc::from("ENST_A");
        let sv = StructuralVariant {
            class: Deletion,
            alleles: vec![SvAllele {
                label: Arc::from("deletion"),
                site: AlleleSite::None,
                alt: 0,
            }],
            tandem: None,
        };
        let contigs = Contigs::of("chr1", &sv);
        let found = [
            SvTranscript {
                transcript: &a,
                contig: 0,
            },
            SvTranscript {
                transcript: &b,
                contig: 0,
            },
        ];
        let alts = [Allele::Symbolic("<DEL>".into())];
        let rows = predict_sv(&sv, &alts, 2031, 2070, &contigs, &found, 1000);
        let ids: Vec<&str> = rows.iter().map(|r| &*r.consequence.transcript_id).collect();
        assert_eq!(ids, ["ENST_A", "ENST_B"]);
        let row = &rows[1].consequence.allele_consequences[0];
        assert_eq!(row.sv_allele.and_then(|n| sv.label(n)), Some("deletion"));
        assert_eq!(row.allele, alts[0]);
        // each row names the transcript it is of
        assert_eq!(&*rows[0].transcript.stable_id, "ENST_A");
    }

    #[test]
    fn a_mate_is_close_to_a_transcript_within_five_thousand_bases_and_on_its_contig() {
        let t = a5();
        let sv = StructuralVariant {
            class: ChromosomeBreakpoint,
            alleles: vec![
                SvAllele {
                    label: Arc::from("A."),
                    site: AlleleSite::Own,
                    alt: 0,
                },
                SvAllele {
                    label: Arc::from("A[chr1:2500["),
                    site: AlleleSite::Mate {
                        chromosome: Arc::from("chr1"),
                        position: 2500,
                    },
                    alt: 0,
                },
                SvAllele {
                    label: Arc::from("A[chr9:2100["),
                    site: AlleleSite::Mate {
                        chromosome: Arc::from("chr9"),
                        position: 2100,
                    },
                    alt: 0,
                },
            ],
            tandem: None,
        };
        let contigs = Contigs::of("chr1", &sv);
        assert_eq!(contigs.index_of("chr1"), 0);
        assert_eq!(contigs.index_of("1"), 0);
        let on_own = SvTranscript {
            transcript: &t,
            contig: 0,
        };
        let on_mate = SvTranscript {
            transcript: &t,
            contig: 1,
        };
        let close = |found: &SvTranscript<'_>| -> Vec<String> {
            sv.alleles
                .iter()
                .filter(|a| is_close_to_feature(a, found, &contigs, 2021, 2021))
                .map(|a| a.label.to_string())
                .collect()
        };
        assert_eq!(close(&on_own), ["A.", "A[chr1:2500["]);
        assert_eq!(close(&on_mate), ["A[chr9:2100["]);
    }
}
