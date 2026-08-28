use fastvep_core::{Allele, Consequence, GenomicPosition, Impact, Strand};
use fastvep_genome::{is_mitochondrial, mitochondrial_codon_table, CodonTable, Transcript};
use std::sync::Arc;

use crate::splice;

/// Result of consequence prediction for a variant against a transcript.
#[derive(Debug, Clone)]
pub struct TranscriptConsequence {
    pub transcript_id: Arc<str>,
    pub gene_id: Arc<str>,
    pub gene_symbol: Option<Arc<str>>,
    pub biotype: Arc<str>,
    pub allele_consequences: Vec<AlleleConsequenceResult>,
    pub canonical: bool,
    pub strand: Strand,
}

/// Consequence result for a single allele against a transcript.
#[derive(Debug, Clone)]
pub struct AlleleConsequenceResult {
    pub allele: Allele,
    pub consequences: Vec<Consequence>,
    pub impact: Impact,
    /// Transcript-relative coordinate pairs, in *strand order*: `start` is the
    /// position of the genomically leftmost affected base and `end` the
    /// rightmost. On the reverse strand a transcript runs right-to-left, so
    /// `start > end` there and the pair is not `(lo, hi)`.
    ///
    /// Use [`AlleleConsequenceResult::protein_range`] (and sort the others the
    /// same way) whenever a span is wanted. A bare `protein_start` is only the
    /// first affected residue on the forward strand; the exported
    /// `Protein_position` matches Ensembl VEP because the output layer sorts the
    /// pair before printing it.
    pub cdna_start: Option<u64>,
    pub cdna_end: Option<u64>,
    pub cds_start: Option<u64>,
    pub cds_end: Option<u64>,
    pub protein_start: Option<u64>,
    pub protein_end: Option<u64>,
    /// `(replaced residues, replacement residues)`.
    ///
    /// The replaced residues are **not** guaranteed to begin at
    /// `protein_start`. They are built from the lower of the two CDS
    /// coordinates for deletions and from `cds_start` otherwise, so for a
    /// shrinking in-frame change on the reverse strand `protein_start` is the
    /// *end* of the affected range and the residues begin one codon earlier.
    /// Over 37,122 in-frame ClinVar rows the replaced residues sat at
    /// `protein_start` in 71.1% of cases, at `protein_start - 1` in 11.9%, and
    /// at neither in 17.0%.
    ///
    /// Anything that needs the true anchor - HGVSp, in particular - has to
    /// corroborate it against the transcript peptide rather than trust either
    /// scalar. See `fastvep_hgvs::protein` and issue #89.
    pub amino_acids: Option<(String, String)>,
    pub codons: Option<(String, String)>,
    pub exon: Option<(u32, u32)>,
    pub intron: Option<(u32, u32)>,
    pub distance: Option<i64>,
    /// Full-length peptide of this transcript, in residues. `None` for
    /// non-coding transcripts.
    pub protein_length: Option<u64>,
    /// Whether a premature termination codon at this position is predicted to
    /// escape nonsense-mediated decay (the 50-nt rule; see
    /// [`Transcript::escapes_nmd`]). Populated for exonic positions in coding
    /// transcripts; `None` when the variant is intronic or the transcript has
    /// no usable exon structure.
    ///
    /// For a frameshift the true stop codon lies somewhat downstream of the
    /// variant; the variant's own position is used as the stand-in, which is
    /// the approximation the published PVS1 implementations make.
    pub escapes_nmd: Option<bool>,
}

/// Resolve a codon-table residue against the transcript's own peptide.
///
/// The codon table cannot see readthrough. A selenoprotein's in-frame UGA
/// translates to `*` here, while `build_sequences` has already resolved it to
/// selenocysteine from the annotated CDS extent, so a substitution at that codon
/// reads as `stop_lost` and a synonymous change reads as `stop_retained`.
///
/// The peptide is consulted only when the table says terminator and the peptide
/// disagrees, which is exactly the readthrough case; every other position still
/// decides by the codon table alone.
///
/// This describes the *reference* protein, so it must never be applied to an alt
/// residue that the variant changed - doing so would rewrite a genuine
/// `stop_gained` into the reference residue. Callers pass alt residues through
/// only for codons the variant leaves untouched.
/// Whether the transcript's annotation claims a complete initiator codon.
///
/// A `cds_start_NF` transcript does not, and neither does one whose CDS is
/// annotated as beginning part-way through a codon - the three bases at
/// `cdna_coding_start` are then the tail of a codon that starts before them.
/// Ensembl checks the flag before every start-codon predicate; the phase is the
/// same fact in a different field, and without it every length-changing variant
/// reaching that end would be `start_lost`, because those three bases never
/// read as ATG.
fn start_codon_known(transcript: &Transcript) -> bool {
    transcript.codon_table_start_phase == 0 && !transcript.flags.iter().any(|f| f == "cds_start_NF")
}

fn resolve_readthrough_residue(transcript: &Transcript, codon_index: usize, translated: u8) -> u8 {
    if translated != b'*' {
        return translated;
    }
    transcript
        .peptide
        .as_deref()
        .and_then(|p| p.as_bytes().get(codon_index).copied())
        .filter(|&b| b != b'*')
        .unwrap_or(translated)
}

/// Drop the common prefix and the common suffix of two sequences, Ensembl's
/// `Bio::EnsEMBL::Variation::Utils::Sequence::trim_sequences`.
///
/// Used to ask whether one sequence is an interior slice of the other: if
/// nothing is left of the shorter side, everything it contained matched in
/// place.
fn trim_common_ends<'a>(a: &'a [u8], b: &'a [u8]) -> (&'a [u8], &'a [u8]) {
    let mut lo = 0;
    while lo < a.len() && lo < b.len() && a[lo] == b[lo] {
        lo += 1;
    }
    let mut hi = 0;
    while hi < a.len() - lo && hi < b.len() - lo && a[a.len() - 1 - hi] == b[b.len() - 1 - hi] {
        hi += 1;
    }
    (&a[lo..a.len() - hi], &b[lo..b.len() - hi])
}

impl AlleleConsequenceResult {
    /// The affected residue span as `(lo, hi)`, whichever way the transcript
    /// runs.
    ///
    /// `protein_start`/`protein_end` are strand-ordered, so on the reverse
    /// strand the raw pair arrives reversed - a 3 bp deletion spanning residues
    /// 3 and 4 comes back as `protein_start = 4`, `protein_end = 3`. Every
    /// consumer that wants a span therefore has to sort the pair, and each one
    /// that open-coded it was one `min`/`max` away from an off-by-one that
    /// nothing downstream could detect.
    ///
    /// A single present coordinate describes a single residue.
    pub fn protein_range(&self) -> Option<(u64, u64)> {
        match (self.protein_start, self.protein_end) {
            (Some(s), Some(e)) => Some((s.min(e), s.max(e))),
            (Some(s), None) => Some((s, s)),
            (None, Some(e)) => Some((e, e)),
            (None, None) => None,
        }
    }
}

/// Full prediction result for a variant.
#[derive(Debug, Clone)]
pub struct PredictionResult {
    pub transcript_consequences: Vec<TranscriptConsequence>,
    pub most_severe: Option<Consequence>,
}

/// The consequence prediction engine.
/// A ref/alt pair rendered for display: amino acids (`("T", "M")`) or codons
/// (`("aCg", "aTg")`). `None` when the change does not produce one.
pub type DisplayPair = Option<(String, String)>;

/// What a change inside the CDS resolves to.
///
/// `additional` carries a second SO term for the changes that earn two. Ensembl
/// evaluates every predicate in `VariationEffect.pm` independently and keeps all
/// that hold, so one delins can be both `stop_gained` and
/// `protein_altering_variant` - 164 of the 1,803 coding rows the ClinVar 2-star
/// in-frame delins produce are exactly that pair. One slot rather than a `Vec`
/// because this value is built once per (variant x transcript x allele) and a
/// `Vec` there is an allocation on the hot path; no shape in that set yields
/// more than two coding terms.
pub struct CodingChange {
    pub consequence: Consequence,
    pub additional: Option<Consequence>,
    pub amino_acids: DisplayPair,
    pub codons: DisplayPair,
}

/// The codon window a CDS-internal change edits, translated on both sides.
///
/// This is the whole of Ensembl's coding model. Every term in
/// `VariationEffect.pm` from `start_lost` to `synonymous_variant` is decided
/// by comparing these two peptides and these two codon strings, and every
/// `Amino_acids` / `Codons` value VEP prints is one of them. Building it once
/// and reading the terms off it is not a refactor of convenience: the four
/// shapes that used to have their own arithmetic (SNV/MNV, pure insertion,
/// pure deletion, delins) each got a different subset of the model right,
/// and the disagreements were exactly where the subsets differed.
struct CodonWindow {
    /// The window's residues before and after the change. `X` stands for a
    /// trailing codon the window does not complete, which is how a frameshift's
    /// alternate peptide ends. Either side is empty when the window is - an
    /// insertion between two codons replaces no residues - and stays empty here
    /// rather than becoming `-`, because every predicate below tests it as a
    /// string and `-` would be a residue to them. Ensembl draws the same line:
    /// `_get_peptide_alleles` maps its `-` to `''` before any predicate sees it.
    ref_aas: String,
    alt_aas: String,
    ref_window: Vec<u8>,
    alt_window: Vec<u8>,
    /// Rendered with the replaced span uppercase; `-` for an empty side.
    ref_codons: String,
    alt_codons: String,
    /// CDS bases the reference covers, and bases the alternate supplies.
    ref_len: usize,
    alt_len: usize,
    /// 1-based protein positions of the window's first and last codons.
    /// Ensembl's `translation_start` / `translation_end`, so an empty window -
    /// an insertion between two codons - has `tl_end == tl_start - 1`.
    tl_start: usize,
    tl_end: usize,
    /// Whether the codon the change starts in runs off the end of the CDS -
    /// Ensembl's `partial_codon`. It is a property of the window's *first*
    /// codon, not its last: a window that merely ends on an incomplete codon is
    /// still a codon edit, and Ensembl reports `inframe_deletion` for one.
    partial_codon: bool,
}

impl CodingChange {
    /// The common case: one term, with whatever ref/alt pairs go with it.
    fn single(consequence: Consequence, amino_acids: DisplayPair, codons: DisplayPair) -> Self {
        Self {
            consequence,
            additional: None,
            amino_acids,
            codons,
        }
    }
}

pub struct ConsequencePredictor {
    pub upstream_distance: u64,
    pub downstream_distance: u64,
    codon_table: CodonTable,
    /// Vertebrate mitochondrial codon table (NCBI translation table 2), used
    /// instead of `codon_table` whenever the transcript being predicted is on
    /// the mitochondrial chromosome (see [`is_mitochondrial`]). Built once
    /// up front so per-allele translation doesn't reconstruct the table.
    mt_codon_table: CodonTable,
}

impl ConsequencePredictor {
    pub fn new(upstream_distance: u64, downstream_distance: u64) -> Self {
        Self {
            upstream_distance,
            downstream_distance,
            codon_table: CodonTable::standard(),
            mt_codon_table: mitochondrial_codon_table(),
        }
    }

    /// Select the codon table to use for a given transcript: the vertebrate
    /// mitochondrial table (NCBI table 2) for MT transcripts, the standard
    /// nuclear table otherwise. AGA/AGG (Arg->Stop), ATA (Ile->Met), and TGA
    /// (Stop->Trp) all differ between the two, so using the wrong table on an
    /// MT variant silently mis-predicts stop_gained/missense/synonymous.
    fn codon_table_for(&self, transcript: &Transcript) -> &CodonTable {
        if is_mitochondrial(&transcript.chromosome) {
            &self.mt_codon_table
        } else {
            &self.codon_table
        }
    }

    /// Predict consequences of a variant against a set of transcripts.
    pub fn predict(
        &self,
        position: &GenomicPosition,
        ref_allele: &Allele,
        alt_alleles: &[Allele],
        transcripts: &[&Transcript],
        ref_seq: Option<&[u8]>,
    ) -> PredictionResult {
        let mut transcript_consequences = Vec::new();

        for transcript in transcripts {
            let tc =
                self.predict_transcript(position, ref_allele, alt_alleles, transcript, ref_seq);
            transcript_consequences.push(tc);
        }

        let all_consequences: Vec<Consequence> = transcript_consequences
            .iter()
            .flat_map(|tc| {
                tc.allele_consequences
                    .iter()
                    .flat_map(|ac| ac.consequences.iter().copied())
            })
            .collect();

        let most_severe = Consequence::most_severe(&all_consequences);

        PredictionResult {
            transcript_consequences,
            most_severe,
        }
    }

    fn predict_transcript(
        &self,
        position: &GenomicPosition,
        ref_allele: &Allele,
        alt_alleles: &[Allele],
        transcript: &Transcript,
        ref_seq: Option<&[u8]>,
    ) -> TranscriptConsequence {
        let allele_consequences: Vec<AlleleConsequenceResult> = alt_alleles
            .iter()
            .map(|alt| self.predict_allele(position, ref_allele, alt, transcript, ref_seq))
            .collect();

        TranscriptConsequence {
            transcript_id: transcript.stable_id.clone(),
            gene_id: transcript.gene.stable_id.clone(),
            gene_symbol: transcript.gene.symbol.clone(),
            biotype: transcript.biotype.clone(),
            allele_consequences,
            canonical: transcript.canonical,
            strand: transcript.strand,
        }
    }

    fn predict_allele(
        &self,
        position: &GenomicPosition,
        ref_allele: &Allele,
        alt_allele: &Allele,
        transcript: &Transcript,
        _ref_seq: Option<&[u8]>,
    ) -> AlleleConsequenceResult {
        // `Allele::Missing` covers the VCF placeholder alleles that assert no
        // alternate sequence at this site: `*` (spanning upstream deletion)
        // and the non-variant forms `.` / `<NON_REF>` / `<*>` that the VCF
        // reader maps here. There is no alternate sequence to translate, so
        // emitting any consequence for them is meaningless. VEP likewise
        // gives `*` no consequence.
        if *alt_allele == Allele::Missing {
            return AlleleConsequenceResult {
                allele: alt_allele.clone(),
                consequences: Vec::new(),
                impact: Impact::Modifier,
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
                distance: None,
                protein_length: None,
                escapes_nmd: None,
            };
        }

        let var_start = position.start;
        let var_end = position.end;
        let tr_start = transcript.start;
        let tr_end = transcript.end;

        let mut consequences = Vec::new();
        let mut cds_start = None;
        let mut cds_end = None;
        let mut protein_start = None;
        let mut protein_end = None;
        let mut amino_acids = None;
        let mut codons = None;
        let mut distance = None;

        // 1. Check if variant overlaps the transcript at all
        let overlaps = var_start <= tr_end && var_end >= tr_start;

        if !overlaps {
            // Check upstream/downstream
            let dist = self.distance_to_transcript(var_start, var_end, transcript);
            if let Some(d) = dist {
                distance = Some(d);
                let abs_dist = d.unsigned_abs();
                if abs_dist <= self.upstream_distance {
                    match transcript.strand {
                        Strand::Forward => {
                            if var_end < tr_start {
                                consequences.push(Consequence::UpstreamGeneVariant);
                            } else {
                                consequences.push(Consequence::DownstreamGeneVariant);
                            }
                        }
                        Strand::Reverse => {
                            if var_start > tr_end {
                                consequences.push(Consequence::UpstreamGeneVariant);
                            } else {
                                consequences.push(Consequence::DownstreamGeneVariant);
                            }
                        }
                    }
                }
            }

            if consequences.is_empty() {
                consequences.push(Consequence::IntergenicVariant);
            }

            let impact = Consequence::worst_impact(&consequences).unwrap_or(Impact::Modifier);
            return AlleleConsequenceResult {
                allele: alt_allele.clone(),
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
                distance,
                protein_length: None,
                escapes_nmd: None,
            };
        }

        // 2. Map to cDNA coordinates
        let cdna_start = transcript.genomic_to_cdna(var_start);
        let cdna_end = transcript.genomic_to_cdna(var_end);

        // 3. Determine exon/intron location
        // Use range-based overlap for exon detection to handle large indels
        let exon_info = transcript
            .exon_at(var_start)
            .or_else(|| transcript.exon_overlapping(var_start, var_end))
            .map(|(i, t)| (i as u32 + 1, t as u32));
        let intron_info = transcript
            .intron_at(var_start)
            .map(|(i, t)| (i as u32 + 1, t as u32));

        let in_exon = exon_info.is_some();
        let in_intron = intron_info.is_some();

        // 4. Check splice sites (always check regardless of coding status)
        //
        // Every term comes from one pass over the introns, because Ensembl's
        // suppression rules are not a fallback chain: a donor variant that also
        // hits the 5th base earns both terms, and so does an acceptor variant in
        // the polypyrimidine tract. Only `splice_region_variant` is displaced by
        // the more specific terms, and `splice_effects` applies that itself.
        //
        // Measured over a 6,600-variant ClinVar sample against real VEP 115.1,
        // this took the splice-term error from 1,503 rows of 150,725 to none.
        let splice = splice::splice_effects(transcript, var_start, var_end, ref_allele, alt_allele);
        if splice.donor {
            consequences.push(Consequence::SpliceDonorVariant);
        }
        if splice.acceptor {
            consequences.push(Consequence::SpliceAcceptorVariant);
        }
        if splice.donor_5th_base {
            consequences.push(Consequence::SpliceDonorFifthBaseVariant);
        }
        if splice.donor_region {
            consequences.push(Consequence::SpliceDonorRegionVariant);
        }
        if splice.polypyrimidine_tract {
            consequences.push(Consequence::SplicePolypyrimidineTractVariant);
        }
        if splice.region {
            consequences.push(Consequence::SpliceRegionVariant);
        }

        // 5. Coding vs non-coding transcript
        if transcript.is_coding() {
            let coding_start = transcript.coding_region_start.unwrap_or(0);
            let coding_end = transcript.coding_region_end.unwrap_or(0);

            // Map to CDS coordinates
            if let Some(cs) = cdna_start {
                cds_start = transcript.cdna_to_cds(cs);
            }
            if let Some(ce) = cdna_end {
                cds_end = transcript.cdna_to_cds(ce);
            }
            if let Some(cds_s) = cds_start {
                protein_start = Some(Transcript::cds_to_protein(cds_s));
            }
            if let Some(cds_e) = cds_end {
                protein_end = Some(Transcript::cds_to_protein(cds_e));
            }

            // Which regions the variant *reaches*, not which region it starts
            // in. A variant may span the CDS/UTR boundary, and then it belongs
            // to both: VEP reports the pair, `stop_lost&3_prime_UTR_variant`.
            //
            // Testing only `var_start` - and choosing one branch from the
            // result - made the answer depend on the strand, because the
            // genomic start of a reverse-strand variant is its *last* base in
            // transcript order. A delins over FOXL2's stop codon came back as
            // `3_prime_UTR_variant` alone, MODIFIER where VEP says HIGH, and
            // the mirror case on the forward strand lost `start_lost` the same
            // way (#100). Each region is now its own additive test.
            let hits_coding =
                in_exon && self.is_in_coding_region(var_start, var_end, coding_start, coding_end);
            let hits_5_utr = in_exon && self.overlaps_5_utr(var_start, var_end, transcript);
            let hits_3_utr = in_exon && self.overlaps_3_utr(var_start, var_end, transcript);

            if hits_5_utr {
                consequences.push(Consequence::FivePrimeUtrVariant);
            }
            if hits_3_utr {
                consequences.push(Consequence::ThreePrimeUtrVariant);
            }
            if hits_coding {
                // Coding exonic variant - determine coding consequence
                let coding_conseq = self.predict_coding_consequence(
                    ref_allele, alt_allele, transcript, cds_start, cds_end, cdna_start, cdna_end,
                );
                if let Some(change) = coding_conseq {
                    consequences.push(change.consequence);
                    if let Some(extra) = change.additional {
                        consequences.push(extra);
                    }
                    amino_acids = change.amino_acids;
                    codons = change.codons;
                } else {
                    consequences.push(Consequence::CodingSequenceVariant);
                }
            }
            // `intron_variant` is not a fallback. Ensembl reports it whenever
            // the variant reaches the intron's interior, beside whatever exonic
            // or splice term it also earned - `splice_donor_variant,intron_variant`
            // is its usual answer for a deletion over a donor site. Reserving it
            // for a variant that earned nothing else dropped it from 1,294 rows
            // of a 6,600-variant ClinVar sample.
            if splice.intronic {
                consequences.push(Consequence::IntronVariant);
            }
        } else {
            // Non-coding transcript
            if in_exon {
                consequences.push(Consequence::NonCodingTranscriptExonVariant);
            } else if in_intron {
                consequences.push(Consequence::NonCodingTranscriptVariant);
            }
            if splice.intronic {
                consequences.push(Consequence::IntronVariant);
            }
        }

        // If still no consequences, add catch-all
        if consequences.is_empty() {
            if transcript.is_coding() {
                consequences.push(Consequence::CodingTranscriptVariant);
            } else {
                consequences.push(Consequence::NonCodingTranscriptVariant);
            }
        }

        // Add NMD_transcript_variant modifier for nonsense_mediated_decay transcripts
        if &*transcript.biotype == "nonsense_mediated_decay" {
            consequences.push(Consequence::NmdTranscriptVariant);
        }

        // Deduplicate
        consequences.sort_by_key(|c| c.rank());
        consequences.dedup();

        let impact = Consequence::worst_impact(&consequences).unwrap_or(Impact::Modifier);

        // PVS1 decision-tree inputs (Abou Tayoun 2018). Both are properties of
        // the transcript and the position, so they are derived here where the
        // exon structure is in hand rather than re-derived downstream.
        let protein_length = transcript
            .is_coding()
            .then(|| transcript.peptide_length())
            .flatten();
        let escapes_nmd = if transcript.is_coding() {
            cdna_start.and_then(|p| transcript.escapes_nmd(p))
        } else {
            None
        };

        AlleleConsequenceResult {
            allele: alt_allele.clone(),
            consequences,
            impact,
            cdna_start,
            cdna_end,
            cds_start,
            cds_end,
            protein_start,
            protein_end,
            amino_acids,
            codons,
            exon: exon_info,
            intron: intron_info,
            distance,
            protein_length,
            escapes_nmd,
        }
    }

    /// Predict the coding consequence (missense, synonymous, frameshift, etc.)
    // Each argument is an independent coordinate or allele; grouping them into
    // a struct would only move the argument list to the call site.
    #[allow(clippy::too_many_arguments)]
    fn predict_coding_consequence(
        &self,
        ref_allele: &Allele,
        alt_allele: &Allele,
        transcript: &Transcript,
        cds_start: Option<u64>,
        cds_end: Option<u64>,
        cdna_start: Option<u64>,
        cdna_end: Option<u64>,
    ) -> Option<CodingChange> {
        // A variant reaching past either end of the CDS is not a codon edit:
        // the bases falling in the UTR belong to no codon, and the frame
        // arithmetic below would count them anyway - which is how a delins
        // over FOXL2's stop codon came back `frameshift_variant` on the
        // forward strand and never reached this function at all on the reverse
        // one. VEP reaches the same conclusion structurally: `cds_start` /
        // `cds_end` are undefined when that end of the variant maps outside
        // the coding region, and `frameshift` returns 0 unless both are
        // defined, so the term comes from its start/stop codon predicates.
        //
        // Test the cDNA spans rather than `cds_start.is_none()`: a CDS
        // coordinate is also None for a flank that falls in an intron, and an
        // indel at an exon edge is an ordinary indel.
        let reaches_past_cds = match (
            cdna_start,
            cdna_end,
            transcript.cdna_coding_start,
            transcript.cdna_coding_end,
        ) {
            (Some(s), Some(e), Some(coding_s), Some(coding_e)) => {
                s.min(e) < coding_s || s.max(e) > coding_e
            }
            _ => false,
        };
        if reaches_past_cds {
            return self.predict_cds_boundary_consequence(
                ref_allele, alt_allele, transcript, cdna_start, cdna_end,
            );
        }
        // One codon window answers every shape. A CDS-internal change is a
        // codon edit whatever its allele string looks like: an SNV, a
        // multi-nucleotide substitution, a pure insertion, a pure deletion and a
        // delins all differ only in how many bases each side of the window
        // contributes. The four separate code paths this replaced disagreed with
        // real VEP 115.1 on 16,000 `Amino_acids` values and 18,000 `Codons`
        // values over a 6,600-variant ClinVar sample, each in its own way.
        let Some(window) =
            self.build_codon_window(transcript, ref_allele, alt_allele, cds_start, cds_end)
        else {
            // The window ran off the end of an incomplete CDS: Ensembl's
            // `partial_codon`, which suppresses every codon term and reports
            // `incomplete_terminal_codon_variant` beside
            // `coding_sequence_variant`. Anything else that has no window - a
            // change not contiguous in CDS space - is `coding_sequence_variant`
            // alone, which is what `None` gives the caller.
            let partial_codon = transcript
                .translateable_seq
                .as_deref()
                .zip(cds_start)
                .is_some_and(|(seq, cds_s)| {
                    let first = cds_s.saturating_sub(1) as usize;
                    (first / 3) * 3 + 3 > seq.len()
                });
            return partial_codon.then_some(CodingChange {
                consequence: Consequence::CodingSequenceVariant,
                additional: Some(Consequence::IncompleteTerminalCodonVariant),
                amino_acids: None,
                codons: None,
            });
        };

        // Ensembl's `_overlaps_start_codon`: a cDNA-space overlap with the three
        // bases at `cdna_coding_start`. A `cds_start_NF` transcript, or one whose
        // CDS begins part-way through a codon, does not claim to carry an
        // initiator, so no term about one is available for it.
        //
        // Ensembl's CDS coordinates run in transcript order on both strands;
        // ours are strand-ordered, so a reverse-strand deletion arrives with its
        // ends swapped and the overlap has to be asked of the sorted pair. An
        // insertion is the exception both ways round: it is the zero-length
        // interval `end = start - 1`, and keeping it inverted is what makes an
        // insertion count only where it sits *inside* the initiator.
        let overlaps_initiator = start_codon_known(transcript)
            && match (cds_start, cds_end) {
                (Some(s), Some(e)) => {
                    let (lo, hi) = (s.min(e), s.max(e));
                    let (s, e) = if *ref_allele == Allele::Deletion {
                        (hi, lo)
                    } else {
                        (lo, hi)
                    };
                    s <= 3 && e >= 1
                }
                _ => false,
            };

        // Ensembl's `_peptide` is the protein without its terminator; ours
        // carries it, because that is what the annotation's own translation
        // ends with.
        let peptide_len = transcript
            .peptide
            .as_deref()
            .map(|p| p.strip_suffix('*').unwrap_or(p).len());
        Some(self.terms_for_window(&window, overlaps_initiator, peptide_len))
    }

    /// The SO terms Ensembl derives from a codon window.
    ///
    /// A port of the predicate graph in `VariationEffect.pm` release/115. The
    /// order matters and is not a severity order: `stop_retained` gates both
    /// `frameshift` (l. 1435) and `stop_gained` (l. 1208), `start_lost` gates
    /// `inframe_insertion` (l. 1100) and `protein_altering_variant` (l. 375),
    /// and `protein_altering_variant` defers to almost everything. Reading them
    /// as a severity chain - pick the worst and stop - loses the second term on
    /// the ~1,100 rows per 6,600 ClinVar variants that earn two.
    ///
    /// `overlaps_initiator` is Ensembl's `_overlaps_start_codon`, which the
    /// caller computes because it is a coordinate question rather than a
    /// peptide one.
    fn terms_for_window(
        &self,
        w: &CodonWindow,
        overlaps_initiator: bool,
        peptide_len: Option<usize>,
    ) -> CodingChange {
        let (ref_pep, alt_pep) = (w.ref_aas.as_str(), w.alt_aas.as_str());
        let extends = |pep: &str| pep.starts_with(ref_pep) || pep.ends_with(ref_pep);

        // `stop_lost` and `stop_retained` are asked first because everything
        // below defers to one or the other.
        let stop_lost = ref_pep.contains('*') && !alt_pep.contains('*');
        // `ref_eq_alt_sequence` (l. 1321). It matters because `frameshift` and
        // `stop_gained` both defer to `stop_retained`: an insertion that leaves
        // the terminator where it was is an `inframe_insertion` to Ensembl
        // however many bases it adds.
        //
        // Two of Ensembl's three clauses are reproduced. The first is not, and
        // that is a deliberate divergence rather than an omission. It reads
        // `$ref_pep eq substr($alt_pep, 0, 1) && $alt_pep =~ /\*/` - the
        // replacement keeps the residue it starts on and introduces a terminator
        // *anywhere* - which does not test whether the annotated terminator
        // survived, and fires on any frameshift whose new stop lands in the
        // first codon after the insertion point. Real VEP 115.1 therefore calls
        // BRCA1 `c.5030_5033dup` `inframe_insertion,stop_retained_variant`,
        // MODERATE, where ClinVar 3-star says Pathogenic and PVS1 applies; the
        // same for `c.1499_1508dup`, BRCA2 `c.3205_3206ins…`, TP53
        // `c.895_919dup` and 30 more in the ClinVar 2-star set. Reproducing it
        // would report a known truncating variant as a moderate in-frame
        // insertion, which is the failure this codebase exists to avoid, so
        // these keep `stop_gained` and `frameshift_variant`.
        //
        // The remaining two do test the terminator. One asks whether it sits at
        // the same residue on both sides; the other whether the edited *protein*
        // still matches the reference over the reference's own length and grows
        // by fewer than three residues past it - only possible when nothing
        // follows the window, so it needs the peptide's length. Deriving that
        // from the CDS length instead assumes the CDS is exactly the peptide
        // plus a terminator, and a `cds_end_NF` transcript is not.
        let grows_past_the_last_residue = |pep_len: usize| {
            // Nothing after the window means the edit cannot displace anything.
            if w.tl_end < pep_len || w.tl_start == 0 {
                return false;
            }
            let before = w.tl_start - 1;
            let grown = before + alt_pep.len();
            // The residues of the window that lie inside the peptide have to
            // survive unchanged at the front of the replacement.
            let kept = pep_len.saturating_sub(before).min(ref_pep.len());
            grown > pep_len && grown - pep_len < 3 && alt_pep.starts_with(&ref_pep[..kept])
        };
        // `stop_retained` also declines on an incomplete terminal codon.
        let stop_retained = !stop_lost
            && !w.partial_codon
            && !alt_pep.is_empty()
            && ((ref_pep.contains('*') && ref_pep.find('*') == alt_pep.find('*'))
                || peptide_len.is_some_and(grows_past_the_last_residue));

        // `frameshift` and `inframe_deletion` both decline when the codon the
        // change starts in is incomplete; `protein_altering_variant` does not.
        let frameshift =
            !stop_retained && !w.partial_codon && !w.alt_len.abs_diff(w.ref_len).is_multiple_of(3);
        let stop_gained = !stop_retained && alt_pep.contains('*') && !ref_pep.contains('*');

        // Ensembl withholds the codon strings from a frameshift
        // (`_get_codon_alleles` returns nothing for one), so neither in-frame
        // term can hold there.
        let (mut inframe_insertion, inframe_deletion) = if frameshift {
            (false, false)
        } else if w.alt_window.len() > w.ref_window.len() {
            // `inframe_insertion` cuts the alternate peptide back to its first
            // terminator before the prefix/suffix test, so an insertion that
            // also introduces a stop is still an insertion.
            let alt_to_stop = match alt_pep.find('*') {
                Some(i) => &alt_pep[..=i],
                None => alt_pep,
            };
            (extends(alt_to_stop), false)
        } else if w.alt_window.len() < w.ref_window.len() {
            // `inframe_deletion` tests the codon strings rather than the
            // peptides: the replacement counts as a deletion only when it
            // reproduces a prefix, a suffix, or a whole-codon interior slice of
            // what it replaced.
            let interior = {
                let (r, a) = trim_common_ends(&w.ref_window, &w.alt_window);
                a.is_empty() && r.len().is_multiple_of(3)
            };
            let matched = !w.partial_codon
                && (w.ref_window.starts_with(&w.alt_window)
                    || w.ref_window.ends_with(&w.alt_window)
                    || interior);
            (false, matched)
        } else {
            (false, false)
        };

        // `_ins_del_start_altered` (l. 1028) edits the transcript and asks
        // whether the initiator survived. Inside a window that starts at the CDS
        // start, that is the same question as whether the edited window still
        // opens on a start codon - unless the edit made the sequence shorter, in
        // which case Ensembl answers "altered" without looking.
        let start_altered = w.alt_len < w.ref_len
            || match w.alt_window.get(..3) {
                Some([a, b, c]) => !CodonTable::is_start(&[*a, *b, *c]),
                _ => true,
            };
        let length_changed = w.ref_len != w.alt_len;
        let start_lost = overlaps_initiator
            && if length_changed {
                start_altered && !(inframe_insertion || inframe_deletion)
            } else {
                // The peptide test, which Ensembl skips when either side is
                // missing or the alternate is the frameshift placeholder.
                !ref_pep.is_empty() && !alt_pep.is_empty() && alt_pep != "X" && !extends(alt_pep)
            };
        // `start_retained_variant` is the complement of the same question, and
        // only for a length change: for a substitution Ensembl reaches it
        // through `_snp_start_altered`, which this window does not model.
        let start_retained = overlaps_initiator && length_changed && !start_altered;
        inframe_insertion &= !start_lost;

        let protein_altering = ref_pep.len() != alt_pep.len()
            && !ref_pep.starts_with('*')
            && !alt_pep.starts_with('*')
            && !extends(alt_pep)
            && !inframe_deletion
            && !start_lost
            && !frameshift;

        let mut terms: Vec<Consequence> = [
            (stop_gained, Consequence::StopGained),
            (frameshift, Consequence::FrameshiftVariant),
            (stop_lost, Consequence::StopLost),
            (start_lost, Consequence::StartLost),
            (inframe_insertion, Consequence::InframeInsertion),
            (inframe_deletion, Consequence::InframeDeletion),
            (protein_altering, Consequence::ProteinAlteringVariant),
            (start_retained, Consequence::StartRetainedVariant),
            (stop_retained, Consequence::StopRetainedVariant),
        ]
        .into_iter()
        .filter_map(|(fired, term)| fired.then_some(term))
        .collect();

        // A same-length replacement that earned nothing above changed the
        // protein or it did not - unless a residue is unknown, in which case
        // nobody can say. Ensembl's `coding_unknown` (l. 1507) reports
        // `coding_sequence_variant` whenever a peptide carries an `X` and no
        // other term holds, and the `X` here comes from an ambiguous reference
        // base: a CDS padded with `N` for an incomplete first codon translates
        // to `X`, and calling that synonymous claims the protein is unchanged
        // when the reference residue was never known.
        let unresolved = ref_pep.contains('X') || alt_pep.contains('X');
        if terms.is_empty() && !length_changed {
            terms.push(if unresolved {
                Consequence::CodingSequenceVariant
            } else if ref_pep == alt_pep {
                Consequence::SynonymousVariant
            } else {
                Consequence::MissenseVariant
            });
        }
        // `incomplete_terminal_codon_variant` sits beside whatever else held,
        // and beside `coding_sequence_variant` when nothing else did.
        if w.partial_codon {
            if terms.is_empty() {
                terms.push(Consequence::CodingSequenceVariant);
            }
            terms.push(Consequence::IncompleteTerminalCodonVariant);
        }
        terms.sort_by_key(|c| c.rank());

        // `-` is how an absent side is written, for residues as for codons.
        let shown = |pep: &str| {
            if pep.is_empty() {
                "-".to_string()
            } else {
                pep.to_string()
            }
        };
        let amino_acids = Some((shown(&w.ref_aas), shown(&w.alt_aas)));
        let codons = Some((w.ref_codons.clone(), w.alt_codons.clone()));
        match terms.first() {
            Some(&first) => CodingChange {
                consequence: first,
                additional: terms.get(1).copied(),
                amino_acids,
                codons,
            },
            // Nothing held: a length change whose peptides Ensembl cannot
            // classify is `coding_sequence_variant`, which the caller supplies.
            None => CodingChange::single(Consequence::CodingSequenceVariant, amino_acids, codons),
        }
    }

    /// Build the codon window for a change whose reference bases are contiguous
    /// in CDS space.
    ///
    /// The window is every codon the reference allele touches. For an insertion
    /// the reference covers no bases at all, so the window is the single codon
    /// it falls inside - and *nothing*, when it falls exactly on a codon
    /// boundary, which is why VEP writes `-/HENKTKGD` for a codon-aligned
    /// insertion rather than repeating the flanking residue on both sides.
    ///
    /// Returns `None` when no window describes the change: it is not contiguous
    /// in the CDS (a delins across a splice junction covers fewer CDS bases than
    /// it has reference bases), or the CDS has no sequence loaded.
    fn build_codon_window(
        &self,
        transcript: &Transcript,
        ref_allele: &Allele,
        alt_allele: &Allele,
        cds_start: Option<u64>,
        cds_end: Option<u64>,
    ) -> Option<CodonWindow> {
        let seq = transcript.translateable_seq.as_deref()?.as_bytes();
        // One end may be missing for an insertion at an exon edge; the branch
        // below reads the pair itself rather than through these.
        let (cds_lo, cds_hi) = match (cds_start, cds_end) {
            (Some(s), Some(e)) => (s.min(e), s.max(e)),
            _ => (0, 0),
        };

        // `first` is the 0-based CDS offset the change starts at; for an
        // insertion it is the offset it is inserted *before*, so a change of no
        // reference bases still has a position.
        let (first, ref_len) = if *ref_allele == Allele::Deletion {
            // Ensembl's zero-length interval: `cds_start == cds_end + 1`, so the
            // insertion goes after the lower of the two.
            //
            // An insertion on an exon's edge has one end in the intron, and so
            // only one CDS coordinate. The insertion point is still determined:
            // it abuts the exonic base, on whichever side the strand puts the
            // intron. `cds_start` comes from the genomic left edge and `cds_end`
            // from the right, so the surviving coordinate says which. Without
            // this an insertion at an exon's first base was
            // `coding_sequence_variant` - 39 rows of a 6,600-variant ClinVar
            // sample, LOW where VEP says HIGH, every one a frameshift.
            let point = match (cds_start, cds_end) {
                (Some(_), Some(_)) if cds_hi == cds_lo + 1 => cds_lo,
                (Some(s), None) if transcript.strand == Strand::Reverse => s,
                (Some(s), None) => s.checked_sub(1)?,
                (None, Some(e)) if transcript.strand == Strand::Forward => e,
                (None, Some(e)) => e.checked_sub(1)?,
                _ => return None,
            };
            (point as usize, 0usize)
        } else {
            let ref_len = ref_allele.len();
            if cds_lo < 1 || cds_hi - cds_lo + 1 != ref_len as u64 {
                return None;
            }
            ((cds_lo - 1) as usize, ref_len)
        };

        // Alt bases in transcript orientation. On the reverse strand the VCF
        // gives them in genomic order, so they are reversed *as well as*
        // complemented: complementing in place reported HPS4
        // c.1060_1061delTCinsAG as `stop_gained` when it is synonymous.
        let alt_cds: Vec<u8> = match alt_allele {
            Allele::Sequence(bases) => match transcript.strand {
                Strand::Forward => bases.clone(),
                Strand::Reverse => bases.iter().rev().map(|&b| complement(b)).collect(),
            },
            _ => Vec::new(),
        };

        let win_start = (first / 3) * 3;
        let win_end = if ref_len == 0 {
            // An insertion on a codon boundary belongs to no codon.
            if first.is_multiple_of(3) {
                win_start
            } else {
                win_start + 3
            }
        } else {
            ((first + ref_len - 1) / 3 + 1) * 3
        };
        // The last codon of an incomplete CDS is short, and the window ends
        // wherever the sequence does: Ensembl reports `aAGA/a` for a deletion
        // there, translating the leftover as `X` rather than declining to
        // describe the change. Refusing the window instead left 3 rows per 230
        // IMPACT-changing ClinVar variants as `coding_sequence_variant` where
        // VEP says `inframe_deletion` or `frameshift_variant`.
        let win_end = win_end.min(seq.len());
        if win_start > win_end {
            return None;
        }

        let ref_window = seq[win_start..win_end].to_vec();
        let lead = first - win_start; // unchanged bases before the change
        let trail = win_end.saturating_sub(first + ref_len); // unchanged bases after it
        let mut alt_window = Vec::with_capacity(lead + alt_cds.len() + trail);
        if lead > ref_window.len() {
            return None;
        }
        alt_window.extend_from_slice(&ref_window[..lead]);
        alt_window.extend_from_slice(&alt_cds);
        alt_window.extend_from_slice(&ref_window[ref_window.len() - trail..]);

        let table = self.codon_table_for(transcript);
        // A codon the edit happens to leave byte-identical to the reference
        // codon at the same offset keeps whatever the reference resolved to, so
        // a selenoprotein's readthrough UGA is not re-read as a new terminator.
        // Every other codon is translated plainly, which is what makes a
        // genuinely introduced stop a stop.
        let translate = |window: &[u8], anchored: bool| -> String {
            let mut pep: String = window
                .chunks_exact(3)
                .enumerate()
                .map(|(i, codon)| {
                    let [a, b, c] = [codon[0], codon[1], codon[2]];
                    let translated = table.translate(&[a, b, c]);
                    let index = win_start / 3 + i;
                    if anchored || ref_window.get(i * 3..i * 3 + 3) == Some(codon) {
                        resolve_readthrough_residue(transcript, index, translated) as char
                    } else {
                        translated as char
                    }
                })
                .collect();
            // Ensembl marks a codon the window does not complete with `X`, but
            // only when the codons before it did not already end translation:
            // `if($partial_codon && $pep ne '*') { $pep .= 'X' }`
            // (`TranscriptVariationAllele.pm` release/115). So `tAac` is `*` and
            // `gttg` is `VX`, and a peptide that merely *contains* a terminator
            // still gets its `X` - `SIFNYIITLFQ*YSFIPYX`.
            if !window.len().is_multiple_of(3) && pep != "*" {
                pep.push('X');
            }
            pep
        };

        // VEP's codon rendering for an edit of any length: the unchanged flanks
        // of the window stay lowercase and exactly the replaced or inserted
        // bases are uppercase, on both sides - `gGg/gTCCCg` for one base
        // replaced by four, and `TAT/CAC` for a whole codon replaced even though
        // its middle base did not change. That is not the same as "uppercase
        // wherever the two differ", which is what produced `TaT/CaC`.
        let render = |window: &[u8], upper_len: usize| -> String {
            if window.is_empty() {
                return "-".to_string();
            }
            window
                .iter()
                .enumerate()
                .map(|(i, &b)| {
                    if i >= lead && i < lead + upper_len {
                        (b as char).to_ascii_uppercase()
                    } else {
                        (b as char).to_ascii_lowercase()
                    }
                })
                .collect()
        };

        Some(CodonWindow {
            ref_aas: translate(&ref_window, true),
            alt_aas: translate(&alt_window, false),
            ref_codons: render(&ref_window, ref_len),
            alt_codons: render(&alt_window, alt_cds.len()),
            tl_start: win_start / 3 + 1,
            tl_end: win_end.div_ceil(3),
            partial_codon: win_start + 3 > seq.len(),
            ref_window,
            alt_window,
            ref_len,
            alt_len: alt_cds.len(),
        })
    }

    /// Consequence for a coding variant that reaches past one end of the CDS.
    ///
    /// Such a variant almost always touches the initiator or the terminator:
    /// reaching past `cdna_coding_start` while still overlapping the CDS means
    /// covering part of `[cdna_coding_start, cdna_coding_start + 2]`, and the
    /// same holds at the other end. So the question is usually only whether
    /// that codon survives the change.
    ///
    /// It is not a guarantee, because the overlap that put us here is measured
    /// in genomic space and can be satisfied by intronic bases alone. Returning
    /// `None` then is correct rather than merely safe: the caller falls back to
    /// `coding_sequence_variant`, which is what a coding change nobody can
    /// resolve any further is.
    ///
    /// The amino-acid and codon pair is deliberately left unset. A variant of
    /// this shape has no single replaced codon to name - part of it is not in
    /// the CDS at all - and the alternative is inventing one: the forward-strand
    /// case used to reach the frameshift formatter and report
    /// `p.Tyr175TerfsTer1` for a change that does not shift any frame.
    fn predict_cds_boundary_consequence(
        &self,
        ref_allele: &Allele,
        alt_allele: &Allele,
        transcript: &Transcript,
        cdna_start: Option<u64>,
        cdna_end: Option<u64>,
    ) -> Option<CodingChange> {
        let coding_start = transcript.cdna_coding_start?;
        let coding_end = transcript.cdna_coding_end?;
        let (s, e) = (cdna_start?, cdna_end?);
        let (cdna_lo, cdna_hi) = (s.min(e), s.max(e));
        if coding_end < coding_start + 2 {
            return None; // no room for a codon at either end
        }

        // VEP's `_overlaps_stop_codon` / `_overlaps_start_codon`: a cDNA-space
        // overlap with the three bases of the terminator or the initiator. A
        // `cds_end_NF` / `cds_start_NF` transcript is one whose annotation does
        // not claim to carry that codon, so neither term is available for it -
        // VEP checks the same two flags before either predicate.
        let flagged = |flag: &str| transcript.flags.iter().any(|f| f == flag);
        let overlaps_stop =
            cdna_lo <= coding_end && cdna_hi >= coding_end - 2 && !flagged("cds_end_NF");
        let overlaps_start =
            cdna_lo <= coding_start + 2 && cdna_hi >= coding_start && start_codon_known(transcript);

        // A variant that swallows the whole CDS destroys both codons; report
        // the more severe of the two terms, which is `stop_lost`.
        if overlaps_stop {
            let still_a_stop = self
                .edited_codon(transcript, cdna_lo, cdna_hi, alt_allele, coding_end - 2)
                .map(|codon| self.codon_table_for(transcript).translate(&codon) == b'*');
            return Some(match still_a_stop {
                Some(true) => CodingChange::single(Consequence::StopRetainedVariant, None, None),
                // Either the terminator now reads as something else, or the
                // edited transcript no longer reaches the position it sat at.
                // Both mean the annotated stop is gone from where it was.
                _ => CodingChange::single(Consequence::StopLost, None, None),
            });
        }

        if overlaps_start {
            // The initiator is not decided by re-reading its codon, and the
            // asymmetry with the terminator above is VEP's, not an oversight.
            // Past the stop codon the 3' UTR simply supplies the next bases and
            // translation runs on, so whatever now sits at that offset *is* the
            // terminator. An initiator instead has to be the first ATG in its
            // context, so `_ins_del_start_altered` spares a change only when it
            // leaves the 5' UTR byte-identical *and* still finds ATG at the
            // coding start; a length change reaching across the junction fails
            // that by construction, and VEP calls it `start_lost`.
            //
            // Reading the codon alone is too permissive here. PLA2G6's coding
            // sequence begins `ATGGATGTCA`, and a `c.-2_3delins` brings its
            // second in-frame ATG onto the coding start: the codon still reads
            // ATG while the CDS has lost four bases.
            //
            // Read the length change off the alleles, not off the cDNA span:
            // an insertion's span is its two flanking bases while it deletes
            // nothing, and a variant spanning an intron has more reference
            // bases than cDNA positions. VEP's `increase_length` /
            // `decrease_length` are likewise properties of the allele pair.
            if ref_allele.len() != alt_allele.len() {
                return Some(CodingChange::single(Consequence::StartLost, None, None));
            }
            // A same-length change cannot shift the frame, so the codon it
            // leaves at the coding start is the whole question.
            let still_a_start = self
                .edited_codon(transcript, cdna_lo, cdna_hi, alt_allele, coding_start)
                .map(|codon| CodonTable::is_start(&codon));
            if still_a_start != Some(true) {
                return Some(CodingChange::single(Consequence::StartLost, None, None));
            }
        }

        None
    }

    /// The three bases at cDNA position `codon_start` *after* `alt_allele`
    /// replaces `[cdna_lo, cdna_hi]`, in transcript orientation.
    ///
    /// VEP decides whether an indel took the terminator or the initiator away
    /// by editing the transcript's own sequence and re-reading one codon at a
    /// fixed offset (`_ins_del_stop_altered`, `_ins_del_start_altered`) rather
    /// than by arithmetic on the variant's length, and it has to: length
    /// arithmetic cannot tell a delins that removes the stop codon from one
    /// that happens to rebuild it, and those are `stop_lost` and
    /// `stop_retained_variant` respectively.
    ///
    /// Read the codon in place instead of materialising the edited sequence.
    /// This runs inside the per-variant loop and the sequence being edited is
    /// the whole transcript, so a copy here would be a copy per variant.
    ///
    /// `None` when the run loaded no sequences, or when the edited cDNA is too
    /// short to reach the codon.
    fn edited_codon(
        &self,
        transcript: &Transcript,
        cdna_lo: u64,
        cdna_hi: u64,
        alt_allele: &Allele,
        codon_start: u64,
    ) -> Option<[u8; 3]> {
        let seq = transcript.spliced_seq.as_deref()?.as_bytes();
        let alt: &[u8] = match alt_allele {
            Allele::Sequence(bases) => bases,
            Allele::Deletion => &[],
            // `*` and the non-variant placeholders assert no alternate
            // sequence, and a symbolic `<DEL>`/`<DUP>` names one without
            // spelling it, so neither has an edited codon to read. Structural
            // alleles are the SV predictor's job.
            Allele::Missing | Allele::Symbolic(_) => return None,
        };
        let alt_len = alt.len() as u64;

        let mut codon = [0u8; 3];
        for (i, slot) in codon.iter_mut().enumerate() {
            let pos = codon_start + i as u64; // 1-based, in the *edited* cDNA
            let base = if pos < cdna_lo {
                *seq.get((pos - 1) as usize)?
            } else if pos < cdna_lo + alt_len {
                // Inside the replacement. The alternate bases arrive in
                // genomic order, so a reverse-strand transcript reads them
                // reverse-complemented, not merely complemented in place.
                let offset = (pos - cdna_lo) as usize;
                match transcript.strand {
                    Strand::Forward => alt[offset],
                    Strand::Reverse => complement(alt[alt.len() - 1 - offset]),
                }
            } else {
                // Past the replacement, so the next original base is the one
                // after the replaced span rather than the one at `pos`.
                let original = cdna_hi + (pos - (cdna_lo + alt_len)) + 1;
                *seq.get((original - 1) as usize)?
            };
            *slot = base.to_ascii_uppercase();
        }
        Some(codon)
    }

    /// Compute amino acids and codons affected by an indel variant.
    /// Returns (amino_acids, codons) tuples.
    /// For frameshifts: ref codon with VEP-style case formatting, truncated alt codon.
    fn distance_to_transcript(
        &self,
        var_start: u64,
        var_end: u64,
        transcript: &Transcript,
    ) -> Option<i64> {
        // For insertions (end < start), use start for distance calculation
        // since start represents the actual insertion position
        let effective_start = var_start.min(var_end);
        let effective_end = var_start.max(var_end);
        if effective_end < transcript.start {
            Some((transcript.start - effective_end) as i64)
        } else if effective_start > transcript.end {
            Some((effective_start - transcript.end) as i64)
        } else {
            Some(0)
        }
    }

    fn is_in_coding_region(
        &self,
        var_start: u64,
        var_end: u64,
        coding_start: u64,
        coding_end: u64,
    ) -> bool {
        var_start <= coding_end && var_end >= coding_start
    }

    /// Does the variant's span reach the untranslated region on the low
    /// genomic side of the coding region - `[transcript.start, cds_start - 1]`?
    ///
    /// Insertions widen the span by the base on either side of the insertion
    /// point, which is what makes an insertion sitting exactly on the coding
    /// boundary a UTR variant. VEP carries the same rule as the two special
    /// cases at the top of `_before_coding` / `_after_coding`.
    ///
    /// The span tested is genomic, so it includes any intron lying between the
    /// transcript's edge and the coding region, and a variant whose only bases
    /// there are intronic still counts. That is VEP's rule rather than an
    /// oversight - `_before_coding` and `_after_coding` are plain `overlap`
    /// calls against the same genomic interval, gated only on the variant
    /// touching an exon somewhere - and matching it keeps the UTR terms
    /// comparable with VEP's. It can only add a MODIFIER term next to a
    /// correctly-derived one, never change the reported impact.
    fn reaches_low_side_utr(&self, var_start: u64, var_end: u64, transcript: &Transcript) -> bool {
        let Some(coding_start) = transcript.coding_region_start else {
            return false;
        };
        if coding_start <= transcript.start {
            return false; // no UTR on this side
        }
        let (lo, hi) = (var_start.min(var_end), var_start.max(var_end));
        lo < coding_start && hi >= transcript.start
    }

    /// The same test on the high genomic side - `[cds_end + 1, transcript.end]`.
    fn reaches_high_side_utr(&self, var_start: u64, var_end: u64, transcript: &Transcript) -> bool {
        let Some(coding_end) = transcript.coding_region_end else {
            return false;
        };
        if coding_end >= transcript.end {
            return false; // no UTR on this side
        }
        let (lo, hi) = (var_start.min(var_end), var_start.max(var_end));
        hi > coding_end && lo <= transcript.end
    }

    /// Does the variant's span reach the 5' UTR? Which genomic side that is
    /// depends on the strand; the overlap test itself does not.
    fn overlaps_5_utr(&self, var_start: u64, var_end: u64, transcript: &Transcript) -> bool {
        if transcript.coding_region_start.is_none() || transcript.coding_region_end.is_none() {
            return false;
        }
        match transcript.strand {
            Strand::Forward => self.reaches_low_side_utr(var_start, var_end, transcript),
            Strand::Reverse => self.reaches_high_side_utr(var_start, var_end, transcript),
        }
    }

    /// Does the variant's span reach the 3' UTR?
    fn overlaps_3_utr(&self, var_start: u64, var_end: u64, transcript: &Transcript) -> bool {
        if transcript.coding_region_start.is_none() || transcript.coding_region_end.is_none() {
            return false;
        }
        match transcript.strand {
            Strand::Forward => self.reaches_high_side_utr(var_start, var_end, transcript),
            Strand::Reverse => self.reaches_low_side_utr(var_start, var_end, transcript),
        }
    }
}

fn complement(base: u8) -> u8 {
    match base {
        b'A' | b'a' => b'T',
        b'T' | b't' => b'A',
        b'C' | b'c' => b'G',
        b'G' | b'g' => b'C',
        other => other,
    }
}

impl Default for ConsequencePredictor {
    fn default() -> Self {
        Self::new(5000, 5000)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn protein_range_sorts_a_reverse_strand_pair() {
        // The trap this exists to remove: on the reverse strand the raw pair
        // arrives with start past end, so a consumer reading `protein_start` as
        // "the first affected residue" is one codon off, and `Protein_position`
        // would print backwards.
        let mut ac = bare_result();

        ac.protein_start = Some(4);
        ac.protein_end = Some(3);
        assert_eq!(ac.protein_range(), Some((3, 4)));

        // Forward strand: already ordered, unchanged.
        ac.protein_start = Some(3);
        ac.protein_end = Some(4);
        assert_eq!(ac.protein_range(), Some((3, 4)));

        // A single coordinate describes a single residue, from either side.
        ac.protein_start = Some(7);
        ac.protein_end = None;
        assert_eq!(ac.protein_range(), Some((7, 7)));
        ac.protein_start = None;
        ac.protein_end = Some(7);
        assert_eq!(ac.protein_range(), Some((7, 7)));

        ac.protein_end = None;
        assert_eq!(ac.protein_range(), None);
    }

    fn bare_result() -> AlleleConsequenceResult {
        AlleleConsequenceResult {
            allele: Allele::from_str("A"),
            consequences: vec![],
            impact: Impact::Modifier,
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
            distance: None,
            protein_length: None,
            escapes_nmd: None,
        }
    }

    use fastvep_genome::{Exon, Gene, Translation};

    fn make_coding_transcript() -> Transcript {
        // A simple protein-coding transcript on forward strand:
        // Exon1: 1000-1200 (UTR: 1000-1049, CDS: 1050-1200)
        // Intron: 1201-1999
        // Exon2: 2000-2300 (all CDS)
        // Intron: 2301-3999
        // Exon3: 4000-5000 (CDS: 4000-4500, UTR: 4501-5000)
        //
        // CDS length: 151 + 301 + 501 = 953 bases
        // translateable_seq: from cDNA pos 51 to 953+50=1003
        let translateable = "ATGGCTTCAAAGCCC".to_string() + &"A".repeat(938); // starts with ATG

        Transcript {
            stable_id: "ENST00000001".into(),
            version: None,
            gene: Gene {
                stable_id: "ENSG00000001".into(),
                symbol: Some("TESTGENE".into()),
                symbol_source: Some("HGNC".into()),
                hgnc_id: None,
                biotype: "protein_coding".into(),
                chromosome: "chr1".into(),
                start: 1000,
                end: 5000,
                strand: Strand::Forward,
            },
            biotype: "protein_coding".into(),
            chromosome: "chr1".into(),
            start: 1000,
            end: 5000,
            strand: Strand::Forward,
            exons: vec![
                Exon {
                    stable_id: "E1".into(),
                    start: 1000,
                    end: 1200,
                    strand: Strand::Forward,
                    phase: -1,
                    end_phase: 0,
                    rank: 1,
                },
                Exon {
                    stable_id: "E2".into(),
                    start: 2000,
                    end: 2300,
                    strand: Strand::Forward,
                    phase: 0,
                    end_phase: 1,
                    rank: 2,
                },
                Exon {
                    stable_id: "E3".into(),
                    start: 4000,
                    end: 5000,
                    strand: Strand::Forward,
                    phase: 1,
                    end_phase: -1,
                    rank: 3,
                },
            ],
            translation: Some(Translation {
                stable_id: "ENSP00000001".into(),
                genomic_start: 1050,
                genomic_end: 4500,
                start_exon_rank: 1,
                start_exon_offset: 50,
                end_exon_rank: 3,
                end_exon_offset: 500,
            }),
            cdna_coding_start: Some(51),
            cdna_coding_end: Some(1003),
            coding_region_start: Some(1050),
            coding_region_end: Some(4500),
            spliced_seq: None,
            translateable_seq: Some(translateable),
            peptide: None,
            canonical: true,
            mane_select: None,
            mane_plus_clinical: None,
            tsl: Some(1),
            appris: None,
            ccds: None,
            protein_id: Some("ENSP00000001".into()),
            protein_version: None,
            swissprot: vec![],
            trembl: vec![],
            uniparc: vec![],
            refseq_id: None,
            source: None,
            gencode_primary: false,
            flags: vec![],
            codon_table_start_phase: 0,
        }
    }

    fn make_noncoding_transcript() -> Transcript {
        Transcript {
            stable_id: "ENST_NC".into(),
            version: None,
            gene: Gene {
                stable_id: "ENSG_NC".into(),
                symbol: Some("NCRNA1".into()),
                symbol_source: None,
                hgnc_id: None,
                biotype: "lncRNA".into(),
                chromosome: "chr1".into(),
                start: 10000,
                end: 12000,
                strand: Strand::Forward,
            },
            biotype: "lncRNA".into(),
            chromosome: "chr1".into(),
            start: 10000,
            end: 12000,
            strand: Strand::Forward,
            exons: vec![
                Exon {
                    stable_id: "E1".into(),
                    start: 10000,
                    end: 10500,
                    strand: Strand::Forward,
                    phase: -1,
                    end_phase: -1,
                    rank: 1,
                },
                Exon {
                    stable_id: "E2".into(),
                    start: 11500,
                    end: 12000,
                    strand: Strand::Forward,
                    phase: -1,
                    end_phase: -1,
                    rank: 2,
                },
            ],
            translation: None,
            cdna_coding_start: None,
            cdna_coding_end: None,
            coding_region_start: None,
            coding_region_end: None,
            spliced_seq: None,
            translateable_seq: None,
            peptide: None,
            canonical: false,
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
        }
    }

    #[test]
    fn test_upstream_variant() {
        let predictor = ConsequencePredictor::default();
        let tr = make_coding_transcript();
        let pos = GenomicPosition::new("chr1", 500, 500, Strand::Forward);
        let result = predictor.predict(
            &pos,
            &Allele::from_str("A"),
            &[Allele::from_str("G")],
            &[&tr],
            None,
        );

        assert_eq!(result.transcript_consequences.len(), 1);
        let tc = &result.transcript_consequences[0];
        assert_eq!(tc.allele_consequences.len(), 1);
        assert!(tc.allele_consequences[0]
            .consequences
            .contains(&Consequence::UpstreamGeneVariant));
        assert_eq!(tc.allele_consequences[0].distance, Some(500));
    }

    #[test]
    fn test_downstream_variant() {
        let predictor = ConsequencePredictor::default();
        let tr = make_coding_transcript();
        let pos = GenomicPosition::new("chr1", 5500, 5500, Strand::Forward);
        let result = predictor.predict(
            &pos,
            &Allele::from_str("A"),
            &[Allele::from_str("G")],
            &[&tr],
            None,
        );

        let ac = &result.transcript_consequences[0].allele_consequences[0];
        assert!(ac
            .consequences
            .contains(&Consequence::DownstreamGeneVariant));
        assert_eq!(ac.distance, Some(500));
    }

    #[test]
    fn test_intergenic_variant() {
        let predictor = ConsequencePredictor::default();
        let tr = make_coding_transcript();
        // Very far away
        let pos = GenomicPosition::new("chr1", 100000, 100000, Strand::Forward);
        let result = predictor.predict(
            &pos,
            &Allele::from_str("A"),
            &[Allele::from_str("G")],
            &[&tr],
            None,
        );

        let ac = &result.transcript_consequences[0].allele_consequences[0];
        assert!(ac.consequences.contains(&Consequence::IntergenicVariant));
    }

    #[test]
    fn test_5_prime_utr_variant() {
        let predictor = ConsequencePredictor::default();
        let tr = make_coding_transcript();
        // Position 1020 is in exon1 (1000-1200), before CDS start (1050)
        let pos = GenomicPosition::new("chr1", 1020, 1020, Strand::Forward);
        let result = predictor.predict(
            &pos,
            &Allele::from_str("A"),
            &[Allele::from_str("G")],
            &[&tr],
            None,
        );

        let ac = &result.transcript_consequences[0].allele_consequences[0];
        assert!(ac.consequences.contains(&Consequence::FivePrimeUtrVariant));
    }

    #[test]
    fn test_3_prime_utr_variant() {
        let predictor = ConsequencePredictor::default();
        let tr = make_coding_transcript();
        // Position 4600 is in exon3 (4000-5000), after CDS end (4500)
        let pos = GenomicPosition::new("chr1", 4600, 4600, Strand::Forward);
        let result = predictor.predict(
            &pos,
            &Allele::from_str("A"),
            &[Allele::from_str("G")],
            &[&tr],
            None,
        );

        let ac = &result.transcript_consequences[0].allele_consequences[0];
        assert!(ac.consequences.contains(&Consequence::ThreePrimeUtrVariant));
    }

    #[test]
    fn test_intron_variant() {
        let predictor = ConsequencePredictor::default();
        let tr = make_coding_transcript();
        // Position 1500 is in intron1 (1201-1999), away from splice sites
        let pos = GenomicPosition::new("chr1", 1500, 1500, Strand::Forward);
        let result = predictor.predict(
            &pos,
            &Allele::from_str("A"),
            &[Allele::from_str("G")],
            &[&tr],
            None,
        );

        let ac = &result.transcript_consequences[0].allele_consequences[0];
        assert!(ac.consequences.contains(&Consequence::IntronVariant));
    }

    #[test]
    fn test_splice_donor_variant() {
        let predictor = ConsequencePredictor::default();
        let tr = make_coding_transcript();
        // Position 1201 is first base of intron1 → splice donor
        let pos = GenomicPosition::new("chr1", 1201, 1201, Strand::Forward);
        let result = predictor.predict(
            &pos,
            &Allele::from_str("G"),
            &[Allele::from_str("A")],
            &[&tr],
            None,
        );

        let ac = &result.transcript_consequences[0].allele_consequences[0];
        assert!(ac.consequences.contains(&Consequence::SpliceDonorVariant));
        assert_eq!(ac.impact, Impact::High);
    }

    #[test]
    fn test_splice_acceptor_variant() {
        let predictor = ConsequencePredictor::default();
        let tr = make_coding_transcript();
        // Position 1999 is last base of intron1 → splice acceptor
        let pos = GenomicPosition::new("chr1", 1999, 1999, Strand::Forward);
        let result = predictor.predict(
            &pos,
            &Allele::from_str("G"),
            &[Allele::from_str("A")],
            &[&tr],
            None,
        );

        let ac = &result.transcript_consequences[0].allele_consequences[0];
        assert!(ac
            .consequences
            .contains(&Consequence::SpliceAcceptorVariant));
        assert_eq!(ac.impact, Impact::High);
    }

    #[test]
    fn test_synonymous_variant() {
        let predictor = ConsequencePredictor::default();
        let tr = make_coding_transcript();
        // First CDS position is at genomic 1050, which is cDNA pos 51, CDS pos 1
        // translateable_seq starts with "ATG" (Met)
        // CDS pos 3 (third base of first codon) - change G to A: ATA still codes for... wait
        // ATG -> Met. Let's change position 3 of ATG from G to something that's still Met: not possible
        // Let's use a different codon. CDS pos 4-6 is "GCT" (Ala). GCC also codes for Ala.
        // Genomic pos for CDS pos 4 = 1050 + 3 = 1053
        // Change T at CDS pos 6 to C: GCT -> GCC both = Ala
        let pos = GenomicPosition::new("chr1", 1055, 1055, Strand::Forward);
        let result = predictor.predict(
            &pos,
            &Allele::from_str("T"),
            &[Allele::from_str("C")],
            &[&tr],
            None,
        );

        let ac = &result.transcript_consequences[0].allele_consequences[0];
        assert!(
            ac.consequences.contains(&Consequence::SynonymousVariant),
            "Expected synonymous, got: {:?}",
            ac.consequences
        );
        assert_eq!(ac.impact, Impact::Low);
    }

    #[test]
    fn test_missense_variant() {
        let predictor = ConsequencePredictor::default();
        let tr = make_coding_transcript();
        // CDS pos 4-6 is "GCT" (Ala). Change first base G to T: TCT = Ser (different!)
        // Genomic pos for CDS pos 4 = 1050 + 3 = 1053
        let pos = GenomicPosition::new("chr1", 1053, 1053, Strand::Forward);
        let result = predictor.predict(
            &pos,
            &Allele::from_str("G"),
            &[Allele::from_str("T")],
            &[&tr],
            None,
        );

        let ac = &result.transcript_consequences[0].allele_consequences[0];
        assert!(
            ac.consequences.contains(&Consequence::MissenseVariant),
            "Expected missense, got: {:?}",
            ac.consequences
        );
        assert_eq!(ac.impact, Impact::Moderate);
    }

    #[test]
    fn test_stop_gained() {
        let predictor = ConsequencePredictor::default();
        let tr = make_coding_transcript();
        // CDS pos 4-6 is "GCT". Change to "TAA" (stop) → need to change pos 4,5,6
        // For simplicity, change CDS pos 4: G->T, pos 5: C->A, pos 6: T->A
        // But our predictor works one SNV at a time. Let's pick a codon that's one base
        // away from a stop. "TCA" (Ser) → change C→A: TAA (stop). But we'd need that codon.
        // Actually, translateable_seq[6..9] = "TCA" (positions 7-9 in 1-based)
        // CDS pos 7 is at genomic 1050+6 = 1056
        // Change T to T (no), we need C at pos 8 to become something.
        // Let's just use translateable[3..6] = "GCT" and change pos 4 (G) to T: "TCT" = Ser
        // That's missense, not stop. Let's try another approach.
        // translateable[9..12] = "AAG" (Lys). Change A at pos 10 to T: TAG = stop!
        // CDS pos 10 is at genomic 1050+9 = 1059
        let pos = GenomicPosition::new("chr1", 1059, 1059, Strand::Forward);
        let result = predictor.predict(
            &pos,
            &Allele::from_str("A"),
            &[Allele::from_str("T")],
            &[&tr],
            None,
        );

        let ac = &result.transcript_consequences[0].allele_consequences[0];
        assert!(
            ac.consequences.contains(&Consequence::StopGained),
            "Expected stop_gained, got: {:?}. translateable[9..12]={:?}",
            ac.consequences,
            &tr.translateable_seq.as_ref().unwrap()[9..12]
        );
        assert_eq!(ac.impact, Impact::High);
    }

    #[test]
    fn test_frameshift_variant() {
        let predictor = ConsequencePredictor::default();
        let tr = make_coding_transcript();
        // Deletion of 1 base in CDS → frameshift
        // CDS pos 4 at genomic 1053
        let pos = GenomicPosition::new("chr1", 1053, 1053, Strand::Forward);
        let result = predictor.predict(
            &pos,
            &Allele::from_str("G"),
            &[Allele::Deletion],
            &[&tr],
            None,
        );

        let ac = &result.transcript_consequences[0].allele_consequences[0];
        assert!(
            ac.consequences.contains(&Consequence::FrameshiftVariant),
            "Expected frameshift, got: {:?}",
            ac.consequences
        );
        assert_eq!(ac.impact, Impact::High);
    }

    #[test]
    fn test_inframe_deletion() {
        let predictor = ConsequencePredictor::default();
        let tr = make_coding_transcript();
        // Deletion of 3 bases in CDS → inframe deletion
        let pos = GenomicPosition::new("chr1", 1053, 1055, Strand::Forward);
        let result = predictor.predict(
            &pos,
            &Allele::from_str("GCT"),
            &[Allele::Deletion],
            &[&tr],
            None,
        );

        let ac = &result.transcript_consequences[0].allele_consequences[0];
        assert!(
            ac.consequences.contains(&Consequence::InframeDeletion),
            "Expected inframe_deletion, got: {:?}",
            ac.consequences
        );
        assert_eq!(ac.impact, Impact::Moderate);
    }

    #[test]
    fn test_noncoding_exon_variant() {
        let predictor = ConsequencePredictor::default();
        let tr = make_noncoding_transcript();
        let pos = GenomicPosition::new("chr1", 10100, 10100, Strand::Forward);
        let result = predictor.predict(
            &pos,
            &Allele::from_str("A"),
            &[Allele::from_str("G")],
            &[&tr],
            None,
        );

        let ac = &result.transcript_consequences[0].allele_consequences[0];
        assert!(ac
            .consequences
            .contains(&Consequence::NonCodingTranscriptExonVariant));
    }

    #[test]
    fn test_start_lost() {
        let predictor = ConsequencePredictor::default();
        let tr = make_coding_transcript();
        // CDS pos 1 (first base of ATG) is at genomic 1050
        // Change A to G: GTG is not a standard start codon
        let pos = GenomicPosition::new("chr1", 1050, 1050, Strand::Forward);
        let result = predictor.predict(
            &pos,
            &Allele::from_str("A"),
            &[Allele::from_str("G")],
            &[&tr],
            None,
        );

        let ac = &result.transcript_consequences[0].allele_consequences[0];
        assert!(
            ac.consequences.contains(&Consequence::StartLost),
            "Expected start_lost, got: {:?}",
            ac.consequences
        );
    }

    #[test]
    fn test_multiple_transcripts() {
        let predictor = ConsequencePredictor::default();
        let tr1 = make_coding_transcript();
        let tr2 = make_noncoding_transcript();

        // Position in tr1's intron, not overlapping tr2
        let pos = GenomicPosition::new("chr1", 1500, 1500, Strand::Forward);
        let result = predictor.predict(
            &pos,
            &Allele::from_str("A"),
            &[Allele::from_str("G")],
            &[&tr1, &tr2],
            None,
        );

        assert_eq!(result.transcript_consequences.len(), 2);
        // tr1: intron variant
        assert!(result.transcript_consequences[0].allele_consequences[0]
            .consequences
            .contains(&Consequence::IntronVariant));
        // tr2: 8500bp away (>5000), so intergenic
        assert!(result.transcript_consequences[1].allele_consequences[0]
            .consequences
            .contains(&Consequence::IntergenicVariant));
    }

    #[test]
    fn test_most_severe_across_transcripts() {
        let predictor = ConsequencePredictor::default();
        let tr = make_coding_transcript();

        // Splice donor is more severe than intron variant
        let pos = GenomicPosition::new("chr1", 1201, 1201, Strand::Forward);
        let result = predictor.predict(
            &pos,
            &Allele::from_str("G"),
            &[Allele::from_str("A")],
            &[&tr],
            None,
        );

        assert_eq!(result.most_severe, Some(Consequence::SpliceDonorVariant));
    }

    #[test]
    fn test_multi_allelic() {
        let predictor = ConsequencePredictor::default();
        let tr = make_coding_transcript();

        // Two alt alleles at a coding position
        let pos = GenomicPosition::new("chr1", 1053, 1053, Strand::Forward);
        let result = predictor.predict(
            &pos,
            &Allele::from_str("G"),
            &[Allele::from_str("T"), Allele::from_str("C")],
            &[&tr],
            None,
        );

        let tc = &result.transcript_consequences[0];
        assert_eq!(tc.allele_consequences.len(), 2);
    }

    // ---- variants spanning either end of the CDS (#100) ----

    /// A single-exon transcript carrying its own sequence, so the terminator
    /// and initiator tests can be decided rather than guessed.
    ///
    /// cDNA position `p` sits at genomic `999 + p` on the forward strand and
    /// `1201 - p` on the reverse one. The layout either way:
    ///
    ///   cDNA   1..10    5' UTR
    ///   cDNA  11..40    CDS - `ATG`, eight `GCT`, `TAA`
    ///   cDNA  41..201   3' UTR - `CC`, then `TGA`, then filler
    ///
    /// The `TGA` three bases into the 3' UTR is deliberate: it lands on the old
    /// terminator's offset after exactly five cDNA bases are removed, which is
    /// what separates `stop_lost` from `stop_retained_variant` below.
    fn make_boundary_transcript(strand: Strand) -> Transcript {
        let five_utr = "C".repeat(10);
        let cds = format!("ATG{}TAA", "GCT".repeat(8));
        let three_utr = format!("CCTGA{}", "G".repeat(156));
        let spliced = format!("{five_utr}{cds}{three_utr}");
        assert_eq!(spliced.len(), 201);

        let (coding_region_start, coding_region_end) = match strand {
            // CDS is cDNA 11..40
            Strand::Forward => (1010, 1039),
            Strand::Reverse => (1161, 1190),
        };

        Transcript {
            stable_id: "ENST_BOUND".into(),
            version: None,
            gene: Gene {
                stable_id: "ENSG_BOUND".into(),
                symbol: Some("BOUNDGENE".into()),
                symbol_source: Some("HGNC".into()),
                hgnc_id: None,
                biotype: "protein_coding".into(),
                chromosome: "chr1".into(),
                start: 1000,
                end: 1200,
                strand,
            },
            biotype: "protein_coding".into(),
            chromosome: "chr1".into(),
            start: 1000,
            end: 1200,
            strand,
            exons: vec![Exon {
                stable_id: "E1".into(),
                start: 1000,
                end: 1200,
                strand,
                phase: -1,
                end_phase: -1,
                rank: 1,
            }],
            translation: Some(Translation {
                stable_id: "ENSP_BOUND".into(),
                genomic_start: coding_region_start,
                genomic_end: coding_region_end,
                start_exon_rank: 1,
                start_exon_offset: 10,
                end_exon_rank: 1,
                end_exon_offset: 40,
            }),
            cdna_coding_start: Some(11),
            cdna_coding_end: Some(40),
            coding_region_start: Some(coding_region_start),
            coding_region_end: Some(coding_region_end),
            spliced_seq: Some(spliced),
            translateable_seq: Some(cds),
            peptide: None,
            canonical: true,
            mane_select: None,
            mane_plus_clinical: None,
            tsl: Some(1),
            appris: None,
            ccds: None,
            protein_id: Some("ENSP_BOUND".into()),
            protein_version: None,
            swissprot: vec![],
            trembl: vec![],
            uniparc: vec![],
            refseq_id: None,
            source: None,
            gencode_primary: false,
            flags: vec![],
            codon_table_start_phase: 0,
        }
    }

    /// Genomic span of a cDNA range on the given strand, as `(start, end)`.
    fn genomic_span(strand: Strand, cdna_lo: u64, cdna_hi: u64) -> (u64, u64) {
        match strand {
            Strand::Forward => (999 + cdna_lo, 999 + cdna_hi),
            Strand::Reverse => (1201 - cdna_hi, 1201 - cdna_lo),
        }
    }

    fn consequences_for(
        strand: Strand,
        cdna_lo: u64,
        cdna_hi: u64,
        ref_allele: &Allele,
        alt_allele: &Allele,
    ) -> Vec<Consequence> {
        let predictor = ConsequencePredictor::default();
        let tr = make_boundary_transcript(strand);
        let (g_start, g_end) = genomic_span(strand, cdna_lo, cdna_hi);
        let pos = GenomicPosition::new("chr1", g_start, g_end, Strand::Forward);
        let result = predictor.predict(
            &pos,
            ref_allele,
            std::slice::from_ref(alt_allele),
            &[&tr],
            None,
        );
        result.transcript_consequences[0].allele_consequences[0]
            .consequences
            .clone()
    }

    /// The reported case: a delins over the terminator and on into the 3' UTR
    /// is `stop_lost&3_prime_UTR_variant`, on either strand.
    ///
    /// It used to depend on the strand. The region was chosen from the
    /// variant's genomic *start*, which on the reverse strand is its last base
    /// in transcript order, so FOXL2's `c.1127_*4delinsCG` matched the 3' UTR
    /// first and never reached the coding branch at all - MODIFIER where VEP
    /// says HIGH. On the forward strand the same shape matched the coding
    /// branch first and lost the UTR term instead.
    #[test]
    fn delins_over_the_stop_codon_is_stop_lost_on_either_strand() {
        for strand in [Strand::Forward, Strand::Reverse] {
            // cDNA 38..44: the whole terminator plus four 3' UTR bases.
            let got = consequences_for(
                strand,
                38,
                44,
                &Allele::from_str("TAACCTG"),
                &Allele::from_str("AC"),
            );
            assert!(
                got.contains(&Consequence::StopLost),
                "{strand:?}: expected stop_lost, got {got:?}"
            );
            assert!(
                got.contains(&Consequence::ThreePrimeUtrVariant),
                "{strand:?}: expected the 3'UTR term alongside it, got {got:?}"
            );
            assert_eq!(Consequence::worst_impact(&got), Some(Impact::High));
        }
    }

    /// A deletion reaching past the terminator is `stop_retained_variant` when
    /// the bases that move up into its place still read as a stop, and
    /// `stop_lost` when they do not. One base of difference decides it, which
    /// is why the codon is re-read from the sequence rather than inferred from
    /// the variant's length.
    #[test]
    fn a_deletion_past_the_stop_codon_distinguishes_lost_from_retained() {
        for strand in [Strand::Forward, Strand::Reverse] {
            // Removing cDNA 38..42 slides the 3' UTR's `TGA` onto the
            // terminator's offset: still a stop.
            let retained = consequences_for(
                strand,
                38,
                42,
                &Allele::from_str("TAACC"),
                &Allele::Deletion,
            );
            assert!(
                retained.contains(&Consequence::StopRetainedVariant),
                "{strand:?}: expected stop_retained_variant, got {retained:?}"
            );
            assert!(!retained.contains(&Consequence::StopLost));

            // One base less, and `CTG` lands there instead: the stop is gone.
            let lost =
                consequences_for(strand, 38, 41, &Allele::from_str("TAAC"), &Allele::Deletion);
            assert!(
                lost.contains(&Consequence::StopLost),
                "{strand:?}: expected stop_lost, got {lost:?}"
            );
            assert!(!lost.contains(&Consequence::StopRetainedVariant));
        }
    }

    /// The mirror case at the other end of the CDS, which had the same bug with
    /// the strands swapped: a length change reaching across the initiator is
    /// `start_lost&5_prime_UTR_variant`.
    #[test]
    fn delins_over_the_start_codon_is_start_lost_on_either_strand() {
        for strand in [Strand::Forward, Strand::Reverse] {
            // cDNA 9..13: two 5' UTR bases and the first three of the CDS.
            let got = consequences_for(
                strand,
                9,
                13,
                &Allele::from_str("CCATG"),
                &Allele::from_str("TT"),
            );
            assert!(
                got.contains(&Consequence::StartLost),
                "{strand:?}: expected start_lost, got {got:?}"
            );
            assert!(
                got.contains(&Consequence::FivePrimeUtrVariant),
                "{strand:?}: expected the 5'UTR term alongside it, got {got:?}"
            );
            assert_eq!(Consequence::worst_impact(&got), Some(Impact::High));
        }
    }

    /// A length change across the initiator is `start_lost` even when an ATG
    /// still happens to land on the coding start, because the reading frame
    /// behind it has moved. PLA2G6's coding sequence begins `ATGGATGTCA`, so a
    /// `c.-2_3delins` there leaves its second in-frame ATG at the coding start
    /// while the CDS is four bases shorter.
    #[test]
    fn start_lost_does_not_depend_on_the_codon_that_lands_at_the_coding_start() {
        // Deleting cDNA 9..13 leaves cDNA 14 (`GCT`'s `T`) onward at the coding
        // start; insert `ATG` so an initiator sits there again.
        for strand in [Strand::Forward, Strand::Reverse] {
            let got = consequences_for(
                strand,
                9,
                13,
                &Allele::from_str("CCATG"),
                &Allele::from_str("ATG"),
            );
            assert!(
                got.contains(&Consequence::StartLost),
                "{strand:?}: expected start_lost despite an ATG at the coding start, got {got:?}"
            );
        }
    }

    /// A variant that stays inside one region keeps exactly the term it had, on
    /// both strands - the additive region tests must not add a second one.
    #[test]
    fn variants_inside_a_single_region_gain_no_extra_term() {
        for strand in [Strand::Forward, Strand::Reverse] {
            // Wholly in the 3' UTR.
            let utr3 =
                consequences_for(strand, 60, 62, &Allele::from_str("GGG"), &Allele::Deletion);
            assert_eq!(
                utr3,
                vec![Consequence::ThreePrimeUtrVariant],
                "{strand:?}: 3'UTR deletion"
            );

            // Wholly in the 5' UTR.
            let utr5 = consequences_for(strand, 3, 5, &Allele::from_str("CCC"), &Allele::Deletion);
            assert_eq!(
                utr5,
                vec![Consequence::FivePrimeUtrVariant],
                "{strand:?}: 5'UTR deletion"
            );

            // Wholly inside the CDS, clear of both terminator and initiator.
            let cds = consequences_for(strand, 20, 22, &Allele::from_str("TGC"), &Allele::Deletion);
            assert_eq!(
                cds,
                vec![Consequence::InframeDeletion],
                "{strand:?}: in-CDS deletion"
            );
        }
    }

    /// A transcript whose annotation does not claim to carry the terminator
    /// cannot lose it. VEP checks `cds_end_NF` before its stop predicates and
    /// `cds_start_NF` before its start ones.
    #[test]
    fn a_cds_end_nf_transcript_does_not_report_stop_lost() {
        let predictor = ConsequencePredictor::default();
        let mut tr = make_boundary_transcript(Strand::Forward);
        tr.flags = vec!["cds_end_NF".to_string()];
        let (g_start, g_end) = genomic_span(Strand::Forward, 38, 44);
        let pos = GenomicPosition::new("chr1", g_start, g_end, Strand::Forward);
        let result = predictor.predict(
            &pos,
            &Allele::from_str("TAACCTG"),
            &[Allele::from_str("AC")],
            &[&tr],
            None,
        );
        let got = &result.transcript_consequences[0].allele_consequences[0].consequences;
        assert!(
            !got.contains(&Consequence::StopLost),
            "cds_end_NF transcript should not report stop_lost, got {got:?}"
        );
        assert!(got.contains(&Consequence::ThreePrimeUtrVariant));
    }

    /// The codon at `codon_start` read the obvious way: materialise the whole
    /// edited cDNA, then slice. `edited_codon` must agree with this in every
    /// case; the only reason it does not work this way is that the sequence
    /// being edited is the whole transcript and the loop runs per variant.
    fn naive_edited_codon(
        seq: &str,
        strand: Strand,
        cdna_lo: u64,
        cdna_hi: u64,
        alt: &Allele,
        codon_start: u64,
    ) -> Option<[u8; 3]> {
        let bytes = seq.as_bytes();
        let alt_tx: Vec<u8> = match alt {
            Allele::Sequence(b) => match strand {
                Strand::Forward => b.clone(),
                Strand::Reverse => b.iter().rev().map(|&x| complement(x)).collect(),
            },
            Allele::Deletion => Vec::new(),
            _ => return None,
        };
        let lo = (cdna_lo - 1) as usize;
        let hi = cdna_hi as usize; // exclusive
        if hi > bytes.len() || lo > hi {
            return None;
        }
        let mut edited: Vec<u8> = Vec::with_capacity(bytes.len());
        edited.extend_from_slice(&bytes[..lo]);
        edited.extend_from_slice(&alt_tx);
        edited.extend_from_slice(&bytes[hi..]);

        let start = (codon_start - 1) as usize;
        if start + 3 > edited.len() {
            return None;
        }
        let mut out = [0u8; 3];
        for (i, slot) in out.iter_mut().enumerate() {
            *slot = edited[start + i].to_ascii_uppercase();
        }
        Some(out)
    }

    /// `edited_codon` reads a codon out of the edited cDNA without building it.
    /// Check it against the build-it-and-slice version over a wide sweep of
    /// shapes, on both strands, including the ones the boundary path never
    /// generates - the arithmetic should not depend on that.
    #[test]
    fn edited_codon_agrees_with_building_the_edited_sequence() {
        let predictor = ConsequencePredictor::default();
        // A deterministic xorshift, so a failure is reproducible.
        let mut state: u64 = 0x2545_F491_4F6C_DD1D;
        let mut next = move || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };

        let bases = *b"ACGT";
        let mut checked = 0usize;
        for strand in [Strand::Forward, Strand::Reverse] {
            let mut tr = make_boundary_transcript(strand);
            let seq = tr.spliced_seq.clone().unwrap();
            let seq_len = seq.len() as u64;

            for _ in 0..4_000 {
                let lo = 1 + next() % (seq_len - 1);
                let span = next() % 12; // 0..11 extra bases
                let hi = (lo + span).min(seq_len);
                let alt_len = (next() % 6) as usize; // 0 => a pure deletion
                let alt = if alt_len == 0 {
                    Allele::Deletion
                } else {
                    Allele::Sequence((0..alt_len).map(|_| bases[(next() % 4) as usize]).collect())
                };
                let codon_start = 1 + next() % (seq_len - 2);

                let got = predictor.edited_codon(&tr, lo, hi, &alt, codon_start);
                let want = naive_edited_codon(&seq, strand, lo, hi, &alt, codon_start);
                assert_eq!(
                    got, want,
                    "{strand:?}: lo={lo} hi={hi} alt={alt:?} codon_start={codon_start}"
                );
                checked += 1;
            }

            // And with no sequence loaded there is nothing to read.
            tr.spliced_seq = None;
            assert_eq!(
                predictor.edited_codon(&tr, 10, 12, &Allele::Deletion, 11),
                None,
                "a transcript with no sequence has no edited codon"
            );
        }
        assert_eq!(checked, 8_000);
    }

    /// A CDS annotated as beginning part-way through a codon has no complete
    /// initiator, so nothing at its start can be `start_lost`. Without the
    /// phase guard every length-changing variant reaching that end would be,
    /// because the three bases there never read as ATG.
    #[test]
    fn a_phase_offset_transcript_does_not_report_start_lost() {
        let predictor = ConsequencePredictor::default();
        let mut tr = make_boundary_transcript(Strand::Forward);
        tr.codon_table_start_phase = 2;
        let (g_start, g_end) = genomic_span(Strand::Forward, 9, 13);
        let pos = GenomicPosition::new("chr1", g_start, g_end, Strand::Forward);
        let result = predictor.predict(
            &pos,
            &Allele::from_str("CCATG"),
            &[Allele::from_str("TT")],
            &[&tr],
            None,
        );
        let got = &result.transcript_consequences[0].allele_consequences[0].consequences;
        assert!(
            !got.contains(&Consequence::StartLost),
            "a phase-offset CDS has no initiator to lose, got {got:?}"
        );
        // The same transcript with a complete initiator does report it, so the
        // guard is what makes the difference and not the variant.
        let tr0 = make_boundary_transcript(Strand::Forward);
        let result = predictor.predict(
            &pos,
            &Allele::from_str("CCATG"),
            &[Allele::from_str("TT")],
            &[&tr0],
            None,
        );
        assert!(result.transcript_consequences[0].allele_consequences[0]
            .consequences
            .contains(&Consequence::StartLost));
    }

    /// An insertion sitting exactly on either coding boundary is a UTR variant
    /// and nothing more. Its span is the two bases it sits between, so it never
    /// deletes a codon, and the boundary path must not claim one was lost.
    #[test]
    fn an_insertion_on_a_coding_boundary_reports_only_the_utr_term() {
        let predictor = ConsequencePredictor::default();
        for strand in [Strand::Forward, Strand::Reverse] {
            let tr = make_boundary_transcript(strand);
            // The initiator is cDNA 11 and the terminator ends at cDNA 40.
            for (cdna_flank, forbidden) in
                [(11u64, Consequence::StartLost), (41, Consequence::StopLost)]
            {
                let genomic = match strand {
                    Strand::Forward => 999 + cdna_flank,
                    Strand::Reverse => 1201 - cdna_flank,
                };
                // Ensembl's zero-length interval: end = start - 1.
                let pos = GenomicPosition::new("chr1", genomic, genomic - 1, Strand::Forward);
                let result = predictor.predict(
                    &pos,
                    &Allele::Deletion,
                    &[Allele::from_str("TTT")],
                    &[&tr],
                    None,
                );
                let got = &result.transcript_consequences[0].allele_consequences[0].consequences;
                assert!(
                    !got.contains(&forbidden),
                    "{strand:?} insertion at cDNA {cdna_flank}: unexpected {forbidden:?} in {got:?}"
                );
            }
        }
    }
    /// A delins in CDS terms, run through the full predictor.
    ///
    /// `cds_lo` is the CDS coordinate of the first replaced base in transcript
    /// order; `reference` and `replacement` are read in transcript order too, so
    /// the reverse-strand case is written the same way and the helper supplies
    /// the reverse complement the VCF would carry.
    fn delins_at(
        strand: Strand,
        cds_lo: u64,
        reference: &str,
        replacement: &str,
    ) -> AlleleConsequenceResult {
        let tr = make_boundary_transcript(strand);
        // CDS n is cDNA n + 10 on this transcript.
        let (lo, hi) = genomic_span(
            strand,
            cds_lo + 10,
            cds_lo + 10 + reference.len() as u64 - 1,
        );
        let orient = |s: &str| -> Allele {
            Allele::Sequence(match strand {
                Strand::Forward => s.as_bytes().to_vec(),
                Strand::Reverse => fastvep_genome::codon::reverse_complement(s.as_bytes()),
            })
        };
        let pos = GenomicPosition::new("chr1", lo, hi, Strand::Forward);
        let alt = orient(replacement);
        let result = ConsequencePredictor::default().predict(
            &pos,
            &orient(reference),
            std::slice::from_ref(&alt),
            &[&tr],
            None,
        );
        result.transcript_consequences[0].allele_consequences[0].clone()
    }

    /// A delins that replaces residues without preserving the reference ones at
    /// either end of the replacement is `protein_altering_variant`, not an
    /// in-frame indel.
    ///
    /// Choosing the term from the direction of the length change - which is what
    /// this code did - called all of these `inframe_deletion` or
    /// `inframe_insertion`. Real VEP 115.1 over the 156 in-frame delins in the
    /// ClinVar 2-star set gives `protein_altering_variant` on 1,231 transcript
    /// rows, `inframe_insertion` on 81 and `inframe_deletion` on none.
    #[test]
    fn a_delins_that_replaces_residues_is_protein_altering_on_either_strand() {
        for strand in [Strand::Forward, Strand::Reverse] {
            // CDS 4-9 is Ala2 Ala3 (`GCTGCT`), replaced by `TGG` (Trp).
            let ac = delins_at(strand, 4, "GCTGCT", "TGG");
            assert!(
                ac.consequences
                    .contains(&Consequence::ProteinAlteringVariant),
                "{strand:?}: expected protein_altering_variant, got {:?}",
                ac.consequences
            );
            assert_eq!(
                ac.amino_acids,
                Some(("AA".to_string(), "W".to_string())),
                "{strand:?}"
            );
            // Every base of the window is replaced, so none of it stays lower.
            assert_eq!(
                ac.codons,
                Some(("GCTGCT".to_string(), "TGG".to_string())),
                "{strand:?}"
            );
        }
    }

    /// The replacement preserving the reference residues at one end of itself is
    /// what makes a delins an in-frame indel. The codon rendering keeps the
    /// unchanged flanks of the window lowercase and uppercases exactly the
    /// replaced and inserted bases - `gGg/gTCCCg`, not "uppercase from the first
    /// difference", which would lose the trailing lowercase base.
    #[test]
    fn a_delins_that_extends_the_reference_residues_is_an_inframe_insertion() {
        for strand in [Strand::Forward, Strand::Reverse] {
            // CDS 5 is the middle base of Ala2's `gCt`. Replacing `C` with
            // `CTGC` rebuilds the window as `GCTGCT`, so the alt peptide still
            // starts with the reference residue and gains one.
            let ac = delins_at(strand, 5, "C", "CTGC");
            assert!(
                ac.consequences.contains(&Consequence::InframeInsertion),
                "{strand:?}: expected inframe_insertion, got {:?}",
                ac.consequences
            );
            assert_eq!(
                ac.amino_acids,
                Some(("A".to_string(), "AA".to_string())),
                "{strand:?}"
            );
            assert_eq!(
                ac.codons,
                Some(("gCt".to_string(), "gCTGCt".to_string())),
                "{strand:?}"
            );
        }
    }

    /// A delins that introduces a terminator earns both terms. Ensembl evaluates
    /// each predicate independently and keeps all that hold, so the new stop
    /// gives `stop_gained` and the changed residue count gives
    /// `protein_altering_variant` - 164 of the 1,803 coding rows in the ClinVar
    /// in-frame delins set are that pair in real VEP.
    ///
    /// `Amino_acids` keeps the whole translated window including what follows
    /// the new terminator, which is what VEP reports (`SL/MEP*S`).
    #[test]
    fn a_delins_introducing_a_terminator_reports_both_terms() {
        for strand in [Strand::Forward, Strand::Reverse] {
            let ac = delins_at(strand, 4, "GCTGCT", "GCCTAGGCC");
            for expected in [Consequence::StopGained, Consequence::ProteinAlteringVariant] {
                assert!(
                    ac.consequences.contains(&expected),
                    "{strand:?}: expected {expected:?}, got {:?}",
                    ac.consequences
                );
            }
            assert_eq!(ac.impact, Impact::High, "{strand:?}");
            assert_eq!(
                ac.amino_acids,
                Some(("AA".to_string(), "A*A".to_string())),
                "{strand:?}"
            );
        }
    }

    /// A replacement whose peptide *begins* with the terminator is `stop_gained`
    /// alone: `protein_altering_variant` declines when either peptide starts
    /// with `*`. That is what separates VEP's `HQ/*` rows, which carry one term,
    /// from its `SL/MEP*S` rows, which carry two.
    #[test]
    fn a_delins_whose_peptide_begins_with_a_terminator_is_stop_gained_alone() {
        for strand in [Strand::Forward, Strand::Reverse] {
            let ac = delins_at(strand, 4, "GCTGCT", "TAG");
            assert_eq!(
                ac.amino_acids,
                Some(("AA".to_string(), "*".to_string())),
                "{strand:?}"
            );
            assert!(
                ac.consequences.contains(&Consequence::StopGained),
                "{strand:?}: got {:?}",
                ac.consequences
            );
            assert!(
                !ac.consequences
                    .contains(&Consequence::ProteinAlteringVariant),
                "{strand:?}: protein_altering_variant must decline, got {:?}",
                ac.consequences
            );
        }
    }

    /// A delins over the initiator is `start_lost`, and that outranks the length
    /// change. Ensembl asks whether the reference residues survived at either
    /// end of the replacement, not whether some ATG still sits at the coding
    /// start.
    #[test]
    fn a_delins_over_the_initiator_is_start_lost() {
        for strand in [Strand::Forward, Strand::Reverse] {
            let ac = delins_at(strand, 1, "ATGGCT", "CCC");
            assert!(
                ac.consequences.contains(&Consequence::StartLost),
                "{strand:?}: got {:?}",
                ac.consequences
            );
            assert_eq!(ac.impact, Impact::High, "{strand:?}");
        }
    }

    /// A change whose reference bases are not contiguous in CDS space has no
    /// codon window, and building one from the low CDS coordinate and the allele
    /// length translates codons the variant never touched. VEP reports
    /// `coding_sequence_variant` with no residues for these, which is what
    /// resolving to no coding term produces.
    ///
    /// This is what a change straddling a splice junction looks like from here:
    /// its two ends are both in the CDS, but further apart than it has bases.
    /// Over a 6,600-variant ClinVar sample it was 156 rows of invented residue
    /// change, ten of them a false `stop_gained`.
    #[test]
    fn a_change_that_is_not_contiguous_in_the_cds_names_no_residues() {
        let predictor = ConsequencePredictor::default();
        let tr = make_boundary_transcript(Strand::Forward);
        for (reference, replacement) in [
            ("GCTGCT", "TGG"),    // delins
            ("GCTGCT", "TGGCCC"), // equal-length MNV
        ] {
            let change = predictor.predict_coding_consequence(
                &Allele::Sequence(reference.as_bytes().to_vec()),
                &Allele::Sequence(replacement.as_bytes().to_vec()),
                &tr,
                Some(4),
                Some(20), // 17 CDS positions for six reference bases
                Some(14),
                Some(30),
            );
            assert!(
                change.is_none(),
                "{reference}/{replacement}: a non-contiguous change must resolve to \
                 no coding term, got {:?}",
                change.map(|c| c.consequence)
            );
        }
    }

    /// An insertion that falls exactly between two codons replaces no codon, so
    /// Ensembl's window is empty on the reference side.
    ///
    /// It reports `-/HENKTKGD` and `-/CATGAG...`, not the flanking residue
    /// repeated on both sides. Anchoring to the preceding codon instead was 791
    /// of 1,808 in-frame insertion rows over a 6,600-variant ClinVar sample, and
    /// it also cost the `dup` collapsing downstream: with the flanking residue
    /// in the way, `p.His553_Asp560dup` came out as an eight-residue `ins`.
    #[test]
    fn a_codon_aligned_insertion_names_no_reference_residue() {
        for strand in [Strand::Forward, Strand::Reverse] {
            // Insert three bases after CDS 3, on a codon boundary.
            let ac = insertion_at(strand, 3, "GGG");
            assert_eq!(
                ac.amino_acids,
                Some(("-".to_string(), "G".to_string())),
                "{strand:?}: {:?}",
                ac.amino_acids
            );
            assert_eq!(
                ac.codons,
                Some(("-".to_string(), "GGG".to_string())),
                "{strand:?}"
            );
            assert!(
                ac.consequences.contains(&Consequence::InframeInsertion),
                "{strand:?}: got {:?}",
                ac.consequences
            );

            // One base further in, the insertion sits inside a codon and the
            // window is that codon.
            let ac = insertion_at(strand, 4, "GGG");
            let (r, a) = ac.amino_acids.clone().unwrap();
            assert_eq!(r.len(), 1, "{strand:?}: {r}/{a}");
            assert_eq!(a.len(), 2, "{strand:?}: {r}/{a}");
        }
    }

    /// Ensembl marks a codon the window does not complete with `X`, but only
    /// when the codons before it did not already end translation.
    ///
    /// `if($partial_codon && $pep ne '*') { $pep .= 'X' }`
    /// (`TranscriptVariationAllele.pm` release/115). Appending it
    /// unconditionally reported `Y/*X` where VEP writes `Y/*` - 288 rows of a
    /// 6,600-variant ClinVar sample.
    #[test]
    fn a_partial_codon_after_a_terminator_adds_no_placeholder() {
        // The fixture's CDS begins ATG GCT GCT ...; inserting one base after
        // CDS 3 shifts the frame from codon 2 on.
        let ac = insertion_at(Strand::Forward, 3, "T");
        let (_, alt) = ac.amino_acids.clone().unwrap();
        assert!(
            alt.ends_with('X') || alt == "*",
            "a frameshift window ends in X unless it ended in a terminator: {alt}"
        );
        assert!(
            !alt.ends_with("*X"),
            "nothing is translated past a terminator: {alt}"
        );
    }

    /// An insertion on an exon's edge has one end in the intron, so only one of
    /// its two CDS coordinates exists - and it is still a frameshift.
    ///
    /// Requiring both coordinates left 39 rows of a 6,600-variant ClinVar sample
    /// as `coding_sequence_variant`, LOW where VEP says HIGH, every one of them
    /// a frameshift at the first or last base of an exon.
    #[test]
    fn an_insertion_on_an_exon_edge_is_still_a_codon_edit() {
        for strand in [Strand::Forward, Strand::Reverse] {
            let tr = make_boundary_transcript(strand);
            // The fixture's first exon ends at cDNA 20, so an insertion between
            // cDNA 20 and 21 has one end in the intron.
            let (edge, _) = genomic_span(strand, 20, 20);
            let (next, _) = genomic_span(strand, 21, 21);
            let (lo, hi) = (edge.min(next), edge.max(next));
            let pos = GenomicPosition::new("chr1", hi, lo, Strand::Forward);
            let alt = Allele::Sequence(b"T".to_vec());
            let result = ConsequencePredictor::default().predict(
                &pos,
                &Allele::Deletion,
                std::slice::from_ref(&alt),
                &[&tr],
                None,
            );
            let got = &result.transcript_consequences[0].allele_consequences[0].consequences;
            assert!(
                got.contains(&Consequence::FrameshiftVariant),
                "{strand:?}: got {got:?}"
            );
        }
    }

    /// An insertion of `bases` after CDS position `cds_lo`, on either strand.
    fn insertion_at(strand: Strand, cds_lo: u64, bases: &str) -> AlleleConsequenceResult {
        let tr = make_boundary_transcript(strand);
        // CDS n is cDNA n + 10 on this transcript, and an insertion is the
        // zero-length interval between two adjacent bases.
        let (a, _) = genomic_span(strand, cds_lo + 10, cds_lo + 10);
        let (b, _) = genomic_span(strand, cds_lo + 11, cds_lo + 11);
        let pos = GenomicPosition::new("chr1", a.max(b), a.min(b), Strand::Forward);
        let alt = Allele::Sequence(match strand {
            Strand::Forward => bases.as_bytes().to_vec(),
            Strand::Reverse => fastvep_genome::codon::reverse_complement(bases.as_bytes()),
        });
        let result = ConsequencePredictor::default().predict(
            &pos,
            &Allele::Deletion,
            std::slice::from_ref(&alt),
            &[&tr],
            None,
        );
        result.transcript_consequences[0].allele_consequences[0].clone()
    }
}
