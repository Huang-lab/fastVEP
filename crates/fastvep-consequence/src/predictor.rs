use fastvep_core::{Allele, Consequence, GenomicPosition, Impact, Strand};
use fastvep_genome::codon::format_codon_window;
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
        // A site is matched against the variant's whole span, because that is
        // what Ensembl does: every test in `_intron_effects` is an `overlap`
        // against `$r_start, $r_end`. Reading `var_start` alone missed any
        // variant whose *second* base landed on the site.
        //
        // `is_splice_donor_5th_base`, `is_splice_donor_region` and
        // `is_splice_polypyrimidine_tract` below still read `var_start` alone.
        // That is a known remaining gap, not a judgement that they differ: they
        // are LOW terms, so missing one cannot drop a HIGH the way missing a
        // donor or acceptor does, and each needs its own check against VEP.
        //
        // Measured over a 6,600-variant ClinVar sample against real VEP 115.1,
        // the span change took donor/acceptor errors from 1,341 to 170 and the
        // extended-term errors from 2,935 to 1,906.
        if splice::is_splice_donor(transcript, var_start, var_end) {
            consequences.push(Consequence::SpliceDonorVariant);
        }
        if splice::is_splice_acceptor(transcript, var_start, var_end) {
            consequences.push(Consequence::SpliceAcceptorVariant);
        }

        // Only add extended splice consequences if not already a donor/acceptor
        let is_essential_splice = consequences.iter().any(|c| {
            matches!(
                c,
                Consequence::SpliceDonorVariant | Consequence::SpliceAcceptorVariant
            )
        });

        if !is_essential_splice {
            let is_donor_5th = splice::is_splice_donor_5th_base(transcript, var_start);
            let is_donor_region = splice::is_splice_donor_region(transcript, var_start);
            if is_donor_5th {
                consequences.push(Consequence::SpliceDonorFifthBaseVariant);
            } else if is_donor_region {
                consequences.push(Consequence::SpliceDonorRegionVariant);
            }
            if splice::is_splice_polypyrimidine_tract(transcript, var_start) {
                consequences.push(Consequence::SplicePolypyrimidineTractVariant);
            }
            // VEP excludes splice_region_variant when a more specific splice term is present:
            // splice_donor_region_variant or splice_donor_5th_base_variant
            if !is_donor_5th
                && !is_donor_region
                && splice::is_splice_region(transcript, var_start, var_end)
            {
                consequences.push(Consequence::SpliceRegionVariant);
            }
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
            // An exonic region term already describes the variant, so the
            // intron term stays reserved for a variant that reached no exonic
            // region at all - the condition the `else if` chain used to encode.
            if !(hits_5_utr || hits_3_utr || hits_coding) && in_intron && !is_essential_splice {
                // VEP excludes intron_variant for positions at splice donor/acceptor sites
                consequences.push(Consequence::IntronVariant);
            }
        } else {
            // Non-coding transcript
            if in_exon {
                consequences.push(Consequence::NonCodingTranscriptExonVariant);
            } else if in_intron {
                if !is_essential_splice {
                    consequences.push(Consequence::IntronVariant);
                }
                consequences.push(Consequence::NonCodingTranscriptVariant);
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
        let cds_pos_start = cds_start?;

        let ref_len = ref_allele.len();
        let alt_len = alt_allele.len();

        // Check if this is a frameshift or in-frame indel
        let is_deletion = *ref_allele != Allele::Deletion && *alt_allele == Allele::Deletion;
        let is_insertion = *ref_allele == Allele::Deletion && *alt_allele != Allele::Missing;
        let is_indel = is_deletion || is_insertion || ref_len != alt_len;

        // A delins - both sides real sequence, lengths differing by a multiple
        // of three - is resolved entirely by `predict_delins_consequence`,
        // including the `None` that means "no codon window describes this".
        // Falling through to the length-only branch below would put back the
        // `inframe_deletion` / `(residue, "X")` answer it exists to replace.
        let is_inframe_delins = matches!(
            (ref_allele, alt_allele),
            (Allele::Sequence(_), Allele::Sequence(_))
        ) && ref_len != alt_len
            && ((ref_len as i64 - alt_len as i64).unsigned_abs() as usize).is_multiple_of(3);
        if is_inframe_delins {
            return self.predict_delins_consequence(
                ref_allele, alt_allele, transcript, cds_start, cds_end,
            );
        }

        if is_indel {
            let (consequence, is_frameshift) = if is_deletion || is_insertion {
                let indel_len = if is_deletion { ref_len } else { alt_len };
                if indel_len % 3 != 0 {
                    (Consequence::FrameshiftVariant, true)
                } else if is_insertion {
                    (Consequence::InframeInsertion, false)
                } else {
                    (Consequence::InframeDeletion, false)
                }
            } else {
                let len_diff = (ref_len as i64 - alt_len as i64).unsigned_abs() as usize;
                if !len_diff.is_multiple_of(3) {
                    (Consequence::FrameshiftVariant, true)
                } else if ref_len > alt_len {
                    (Consequence::InframeDeletion, false)
                } else {
                    (Consequence::InframeInsertion, false)
                }
            };

            // For deletions on reverse strand, cds_start maps to the end of the
            // deletion in CDS space. Use the lower CDS position as the start.
            let cds_pos = if is_deletion {
                match cds_end {
                    Some(ce) => cds_pos_start.min(ce),
                    None => cds_pos_start,
                }
            } else {
                cds_pos_start
            };

            // Try to compute amino acids and codons from translateable_seq
            let (aa_pair, codon_pair) = self.compute_indel_amino_acids(
                transcript,
                cds_pos,
                ref_allele,
                alt_allele,
                is_frameshift,
            );

            return Some(CodingChange::single(consequence, aa_pair, codon_pair));
        }

        // Same length substitution (SNV or MNV).
        //
        // A multi-nucleotide substitution is NOT a sequence of independent
        // single-base changes: it must be translated as one block, over every
        // codon it touches. Two things go wrong if it is treated as a
        // single-codon edit anchored at `cds_pos_start`:
        //
        //  1. A change straddling a codon boundary loses the bases that fall
        //     past the first codon (MUTYH c.1164_1165delinsAT was reported as
        //     `CTg/TAg` -> stop_gained when the real change is synonymous).
        //  2. On the reverse strand `cds_pos_start` is the CDS coordinate of
        //     the *last* changed base, and the alt bases arrive in genomic
        //     order, so they must be reverse-complemented, not merely
        //     complemented in place (HPS4 c.1060_1061delTCinsAG was reported
        //     as `tCC/tGA` -> stop_gained when the real change is p.Ser354=).
        if let Some(ref translateable_seq) = transcript.translateable_seq {
            let seq_bytes = translateable_seq.as_bytes();

            // Alt bases in transcript orientation.
            let alt_bases_cds: Vec<u8> = match alt_allele {
                Allele::Sequence(bases) => match transcript.strand {
                    Strand::Forward => bases.clone(),
                    Strand::Reverse => bases.iter().rev().map(|&b| complement(b)).collect(),
                },
                _ => Vec::new(),
            };

            // CDS coordinate of the first changed base in transcript
            // orientation. On the reverse strand the genomic start maps to the
            // higher CDS coordinate, so take the lower of the two ends.
            let cds_lo = match cds_end {
                Some(ce) => cds_pos_start.min(ce),
                None => cds_pos_start,
            };

            // The changed bases have to be contiguous in CDS space before a
            // codon window can describe them. An MNV straddling a splice
            // junction covers fewer CDS positions than it has bases, and a
            // window built from `cds_lo` and the allele length alone then
            // translates codons the variant never touched: 156 rows over a
            // 6,600-variant ClinVar sample, ten of them a false `stop_gained`.
            // VEP reports `coding_sequence_variant` with no amino acids for
            // these, which is what falling through to `None` produces.
            let contiguous_in_cds = match (cds_start, cds_end) {
                (Some(s), Some(e)) => s.max(e) - s.min(e) + 1 == alt_bases_cds.len() as u64,
                _ => false,
            };

            let n = alt_bases_cds.len();
            if n > 0 && cds_lo >= 1 && contiguous_in_cds {
                let first_changed = (cds_lo - 1) as usize;
                let last_changed = first_changed + n - 1;
                let win_start = (first_changed / 3) * 3;
                let win_end = (last_changed / 3 + 1) * 3;

                if win_end <= seq_bytes.len() {
                    let ref_window = seq_bytes[win_start..win_end].to_vec();
                    let mut alt_window = ref_window.clone();
                    for (i, &base) in alt_bases_cds.iter().enumerate() {
                        alt_window[first_changed - win_start + i] = base;
                    }

                    let table = self.codon_table_for(transcript);
                    let mut ref_aas = String::new();
                    let mut alt_aas = String::new();
                    for c in (0..ref_window.len()).step_by(3) {
                        let r: [u8; 3] = [ref_window[c], ref_window[c + 1], ref_window[c + 2]];
                        let a: [u8; 3] = [alt_window[c], alt_window[c + 1], alt_window[c + 2]];
                        let codon_index = (win_start + c) / 3;
                        let ref_res = resolve_readthrough_residue(
                            transcript,
                            codon_index,
                            table.translate(&r),
                        );
                        // An untouched codon keeps whatever the reference
                        // resolved to; a changed one is translated plainly, so a
                        // newly-introduced stop is still a stop.
                        let alt_res = if a == r { ref_res } else { table.translate(&a) };
                        ref_aas.push(ref_res as char);
                        alt_aas.push(alt_res as char);
                    }

                    let (ref_codon_str, alt_codon_str) =
                        format_codon_window(&ref_window, &alt_window);
                    let codon_pair = Some((ref_codon_str, alt_codon_str));
                    let aa_pair = Some((ref_aas.clone(), alt_aas.clone()));

                    // Start codon is the first codon of the window only when
                    // the window begins at CDS position 1.
                    if win_start == 0 {
                        let first_alt: [u8; 3] = [alt_window[0], alt_window[1], alt_window[2]];
                        if !CodonTable::is_start(&first_alt) {
                            return Some(CodingChange::single(
                                Consequence::StartLost,
                                aa_pair,
                                codon_pair,
                            ));
                        }
                    }

                    if ref_aas == alt_aas {
                        if ref_aas.contains('*') {
                            return Some(CodingChange::single(
                                Consequence::StopRetainedVariant,
                                aa_pair,
                                codon_pair,
                            ));
                        }
                        return Some(CodingChange::single(
                            Consequence::SynonymousVariant,
                            aa_pair,
                            codon_pair,
                        ));
                    }

                    // A stop introduced anywhere in the window is stop_gained;
                    // a reference stop that the change removes is stop_lost.
                    let gained_stop = ref_aas
                        .chars()
                        .zip(alt_aas.chars())
                        .any(|(r, a)| a == '*' && r != '*');
                    if gained_stop {
                        return Some(CodingChange::single(
                            Consequence::StopGained,
                            aa_pair,
                            codon_pair,
                        ));
                    }
                    let lost_stop = ref_aas
                        .chars()
                        .zip(alt_aas.chars())
                        .any(|(r, a)| r == '*' && a != '*');
                    if lost_stop {
                        return Some(CodingChange::single(
                            Consequence::StopLost,
                            aa_pair,
                            codon_pair,
                        ));
                    }

                    return Some(CodingChange::single(
                        Consequence::MissenseVariant,
                        aa_pair,
                        codon_pair,
                    ));
                }
            }
        }

        // Nothing above could read a codon window: the change is not contiguous
        // in CDS space, or the window runs off the end of an incomplete CDS.
        //
        // `missense_variant` used to be returned here for a same-length change,
        // on the reasoning that a substitution in the CDS changes a residue. It
        // does not follow - the change could equally be synonymous, and we have
        // just established that we cannot tell. VEP reports
        // `coding_sequence_variant`, which is what a coding change nobody can
        // resolve any further is, and `None` is how the caller reaches it. Over a
        // 6,600-variant ClinVar sample this was 156 rows of MODERATE where VEP
        // says MODIFIER, each naming an amino-acid change built from codons the
        // variant never touched.
        None
    }

    /// Consequence for a delins inside the CDS: a change where both sides are
    /// real sequence and the lengths differ.
    ///
    /// Such a change is neither a pure insertion nor a pure deletion, and the
    /// terms Ensembl gives it do not follow from the direction of the length
    /// change - they follow from the peptides. `protein_altering_variant`
    /// (`VariationEffect.pm` release/115, l. 375) is what a replacement earns
    /// when it is not a prefix or suffix extension of what it replaces;
    /// `inframe_insertion` (l. 1100) and `inframe_deletion` (l. 1160) are
    /// reserved for the cases where it is. Choosing by length alone - which is
    /// what this code did until now - got the term wrong far more often than
    /// right: over the 156 in-frame delins in the ClinVar 2-star set, real VEP
    /// 115.1 reports `protein_altering_variant` on 1,231 transcript rows,
    /// `inframe_insertion` on 81 and `inframe_deletion` on none.
    ///
    /// It also had nothing to say about the residues. The old path fell through
    /// to `(ref_residue, "X")` with no codons at all, so 1,697 rows carried an
    /// `Amino_acids` of `E/X` and an empty `Codons` - and HGVSp built
    /// `p.Glu721delinsXaa` on top of it, which reads as "one residue replaced
    /// by one unknown residue" for a variant that inserts twelve.
    ///
    /// Returns `None` when the change cannot be read as a codon window: it runs
    /// across a splice junction (so it covers fewer CDS bases than reference
    /// bases), or the window would run off the end of the CDS. VEP reports
    /// `coding_sequence_variant` for those instead of naming a codon change,
    /// and returning `None` is how the caller reaches that term.
    fn predict_delins_consequence(
        &self,
        ref_allele: &Allele,
        alt_allele: &Allele,
        transcript: &Transcript,
        cds_start: Option<u64>,
        cds_end: Option<u64>,
    ) -> Option<CodingChange> {
        let (Allele::Sequence(ref_bases), Allele::Sequence(alt_bases)) = (ref_allele, alt_allele)
        else {
            return None;
        };
        let seq_bytes = transcript.translateable_seq.as_deref()?.as_bytes();

        let (cds_s, cds_e) = (cds_start?, cds_end?);
        let (cds_lo, cds_hi) = (cds_s.min(cds_e), cds_s.max(cds_e));
        if cds_lo < 1 {
            return None;
        }
        // The replaced bases have to be contiguous in CDS space for a codon
        // window to describe them. A delins reaching across a splice junction
        // covers fewer CDS positions than it has reference bases, and one end
        // of it is not in the CDS at all.
        if cds_hi - cds_lo + 1 != ref_bases.len() as u64 {
            return None;
        }

        // Alt bases in transcript orientation. On the reverse strand the VCF
        // gives them in genomic order, so they are reversed *as well as*
        // complemented: complementing in place reported HPS4
        // c.1060_1061delTCinsAG as `stop_gained` when it is synonymous.
        let alt_cds: Vec<u8> = match transcript.strand {
            Strand::Forward => alt_bases.clone(),
            Strand::Reverse => alt_bases.iter().rev().map(|&b| complement(b)).collect(),
        };

        // The affected window is the whole codons the replaced bases touch, so
        // every codon in it carries at least one changed base.
        let first = (cds_lo - 1) as usize;
        let last = first + ref_bases.len() - 1;
        let win_start = (first / 3) * 3;
        let win_end = (last / 3 + 1) * 3;
        if win_end > seq_bytes.len() {
            return None; // the last codon is partial - VEP's `partial_codon`
        }
        let ref_window = &seq_bytes[win_start..win_end];
        let lead = first - win_start; // unchanged bases before the change
        let trail = win_end - 1 - last; // unchanged bases after it

        let mut alt_window = Vec::with_capacity(lead + alt_cds.len() + trail);
        alt_window.extend_from_slice(&ref_window[..lead]);
        alt_window.extend_from_slice(&alt_cds);
        alt_window.extend_from_slice(&ref_window[ref_window.len() - trail..]);
        // In frame by construction: the window is whole codons and the length
        // difference is a multiple of three, which is what routed us here.
        debug_assert!(alt_window.len().is_multiple_of(3));

        let table = self.codon_table_for(transcript);
        let mut ref_aas = String::with_capacity(ref_window.len() / 3);
        for c in (0..ref_window.len()).step_by(3) {
            let codon = [ref_window[c], ref_window[c + 1], ref_window[c + 2]];
            ref_aas.push(resolve_readthrough_residue(
                transcript,
                (win_start + c) / 3,
                table.translate(&codon),
            ) as char);
        }
        let mut alt_aas = String::with_capacity(alt_window.len() / 3);
        for c in (0..alt_window.len()).step_by(3) {
            let codon = [alt_window[c], alt_window[c + 1], alt_window[c + 2]];
            let translated = table.translate(&codon);
            // A codon the edit happens to leave byte-identical to the reference
            // codon at the same offset keeps whatever the reference resolved to,
            // so a selenoprotein's readthrough UGA is not re-read as a new
            // terminator. Every other codon is translated plainly, which is what
            // makes a genuinely introduced stop a stop.
            let residue = if ref_window.get(c..c + 3) == Some(&codon[..]) {
                resolve_readthrough_residue(transcript, (win_start + c) / 3, translated)
            } else {
                translated
            };
            alt_aas.push(residue as char);
        }

        // VEP's codon rendering for a length-changing edit: the unchanged flanks
        // of the window stay lowercase and exactly the replaced or inserted
        // bases are uppercase, on both sides - `gGg/gTCCCg` for one base
        // replaced by four. That is not the same as "uppercase from the first
        // difference", which would lose the lowercase tail.
        let render = |window: &[u8], upper_len: usize| -> String {
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
        let codons = Some((
            render(ref_window, ref_bases.len()),
            render(&alt_window, alt_cds.len()),
        ));
        let amino_acids = Some((ref_aas.clone(), alt_aas.clone()));

        // `inframe_insertion` cuts the alt peptide back to its first terminator
        // before the prefix/suffix test (l. 1124), so an insertion that also
        // introduces a stop is still an insertion. The other predicates read the
        // untruncated peptide.
        let alt_to_stop = match alt_aas.find('*') {
            Some(i) => &alt_aas[..=i],
            None => alt_aas.as_str(),
        };
        let extends = |pep: &str| pep.starts_with(&ref_aas) || pep.ends_with(&ref_aas);

        // `start_lost` (l. 851) asks whether the initiator's residue survived at
        // either end of the replacement, not whether some ATG still sits at the
        // coding start; it takes precedence over every term below, and
        // `protein_altering_variant` defers to it explicitly.
        if win_start == 0 && !extends(&alt_aas) {
            return Some(CodingChange::single(
                Consequence::StartLost,
                amino_acids,
                codons,
            ));
        }

        let ref_has_stop = ref_aas.contains('*');
        let alt_has_stop = alt_aas.contains('*');
        if ref_has_stop {
            // The window carried the annotated terminator. Whether the edited
            // window still ends translation there is the whole question, and it
            // outranks the length change.
            let consequence = if alt_has_stop {
                Consequence::StopRetainedVariant
            } else {
                Consequence::StopLost
            };
            return Some(CodingChange::single(consequence, amino_acids, codons));
        }

        let inframe = if alt_window.len() > ref_window.len() {
            extends(alt_to_stop).then_some(Consequence::InframeInsertion)
        } else {
            // `inframe_deletion` tests the codon strings rather than the
            // peptides: the replacement counts as a deletion only when it
            // reproduces a prefix, a suffix, or a whole-codon interior slice of
            // what it replaced.
            let interior = {
                let (r, a) = trim_common_ends(ref_window, &alt_window);
                a.is_empty() && r.len().is_multiple_of(3)
            };
            (ref_window.starts_with(&alt_window) || ref_window.ends_with(&alt_window) || interior)
                .then_some(Consequence::InframeDeletion)
        };

        // `protein_altering_variant` (l. 375) is the term for a replacement that
        // changes the residue count without preserving the reference residues at
        // either end. It defers to `inframe_deletion` but not to
        // `inframe_insertion`, and its own peptide test is the untruncated one -
        // so an insertion carrying a new stop earns both terms. It also defers
        // to a peptide *starting* with a terminator, which is why `HQ/*` is
        // `stop_gained` alone while `SL/MEP*S` is the pair.
        let protein_altering = ref_aas.len() != alt_aas.len()
            && !ref_aas.starts_with('*')
            && !alt_aas.starts_with('*')
            && !extends(&alt_aas)
            && inframe != Some(Consequence::InframeDeletion);

        let stop_gained = alt_has_stop;
        Some(match (stop_gained, protein_altering, inframe) {
            (true, pa, _) => CodingChange {
                consequence: Consequence::StopGained,
                additional: pa.then_some(Consequence::ProteinAlteringVariant),
                amino_acids,
                codons,
            },
            (false, _, Some(term)) => CodingChange::single(term, amino_acids, codons),
            (false, true, None) => {
                CodingChange::single(Consequence::ProteinAlteringVariant, amino_acids, codons)
            }
            // Same residue count, reference residues preserved, no stop
            // involved: the replacement produced the same protein.
            (false, false, None) => CodingChange::single(
                if ref_aas == alt_aas {
                    Consequence::SynonymousVariant
                } else {
                    Consequence::MissenseVariant
                },
                amino_acids,
                codons,
            ),
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
        // A non-zero start phase says the annotated CDS begins part-way through
        // a codon, so this transcript has no complete initiator to lose and the
        // three bases at `cdna_coding_start` are not one. Without this, every
        // length-changing variant reaching that end of such a transcript would
        // be `start_lost`, because `is_start` can never hold there.
        let start_codon_known = !flagged("cds_start_NF") && transcript.codon_table_start_phase == 0;
        let overlaps_start =
            cdna_lo <= coding_start + 2 && cdna_hi >= coding_start && start_codon_known;

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
    fn compute_indel_amino_acids(
        &self,
        transcript: &Transcript,
        cds_pos: u64,
        ref_allele: &Allele,
        alt_allele: &Allele,
        is_frameshift: bool,
    ) -> (DisplayPair, DisplayPair) {
        let translateable_seq = match transcript.translateable_seq.as_ref() {
            Some(s) => s,
            None => return (None, None),
        };
        let seq_bytes = translateable_seq.as_bytes();
        let cds_idx = (cds_pos - 1) as usize;

        if cds_idx >= seq_bytes.len() {
            return (None, None);
        }

        // Get the codon at the affected position
        let codon_number = cds_idx / 3;
        let codon_offset = cds_idx % 3;
        let codon_start = codon_number * 3;

        if codon_start + 3 > seq_bytes.len() {
            return (None, None);
        }

        let ref_codon = [
            seq_bytes[codon_start],
            seq_bytes[codon_start + 1],
            seq_bytes[codon_start + 2],
        ];
        let ref_aa = resolve_readthrough_residue(
            transcript,
            codon_number,
            self.codon_table_for(transcript).translate(&ref_codon),
        );
        let ref_aa_str = String::from(ref_aa as char);

        if is_frameshift {
            // A frameshift's alt residue is `X` and its alt codon is only the
            // remains of the codon the indel starts in, so the whole edited CDS
            // is never read here. It used to be built anyway - a clone of
            // `translateable_seq` plus a splice, per frameshift variant, thrown
            // away unused. Only the allele shapes still need checking.
            if !matches!(
                (ref_allele, alt_allele),
                (Allele::Sequence(_), Allele::Deletion)
                    | (Allele::Deletion, Allele::Sequence(_))
                    | (Allele::Sequence(_), Allele::Sequence(_))
            ) {
                return (None, None);
            }

            // Build codon display: VEP style with deleted base uppercase
            // ref codon: lowercase bases, uppercase at the deleted position(s)
            let mut ref_codon_display = String::with_capacity(3);
            for (i, &base) in ref_codon.iter().enumerate() {
                if i == codon_offset {
                    ref_codon_display.push((base as char).to_ascii_uppercase());
                } else {
                    ref_codon_display.push((base as char).to_ascii_lowercase());
                }
            }

            // alt codon: show only the remaining bases of the original codon after the indel
            // For a deletion at offset 2 in a 3-base codon: show only the 2 remaining bases
            let alt_codon_display: String = {
                let mut original_codon: Vec<u8> = ref_codon.to_vec();
                match (ref_allele, alt_allele) {
                    (Allele::Sequence(_), Allele::Deletion) => {
                        // Remove the deleted base(s) from the codon
                        let del_len = ref_allele.len().min(3 - codon_offset);
                        let end = (codon_offset + del_len).min(original_codon.len());
                        original_codon.drain(codon_offset..end);
                    }
                    (Allele::Deletion, Allele::Sequence(ins_bases)) => {
                        // Insert bases into the codon at the offset
                        let mut bases = ins_bases.clone();
                        if transcript.strand == Strand::Reverse {
                            bases = bases.iter().map(|&b| complement(b)).collect();
                        }
                        for (j, &b) in bases.iter().enumerate() {
                            original_codon.insert(codon_offset + j, b);
                        }
                    }
                    _ => {}
                }
                original_codon
                    .iter()
                    .map(|&b| (b as char).to_ascii_lowercase())
                    .collect()
            };

            // For frameshifts, alt amino acid is always X (unknown/frameshift)
            // For pure insertions, VEP uses "-" for ref amino acid/codon
            // and only the inserted bases for alt codon
            let (fs_ref_aa, fs_ref_codon, fs_alt_codon) = if *ref_allele == Allele::Deletion {
                let ins_codon = if let Allele::Sequence(ins_bases) = alt_allele {
                    let mut bases = ins_bases.clone();
                    if transcript.strand == Strand::Reverse {
                        bases = bases.iter().map(|&b| complement(b)).collect();
                    }
                    bases
                        .iter()
                        .map(|&b| (b as char).to_ascii_uppercase())
                        .collect::<String>()
                } else {
                    alt_codon_display
                };
                ("-".to_string(), "-".to_string(), ins_codon)
            } else {
                (ref_aa_str, ref_codon_display, alt_codon_display)
            };
            let aa_pair = Some((fs_ref_aa, "X".to_string()));
            let codon_pair = Some((fs_ref_codon, fs_alt_codon));
            (aa_pair, codon_pair)
        } else {
            // In-frame indel: build alt sequence and translate affected codons
            let mut alt_seq: Vec<u8> = seq_bytes.to_vec();
            match (ref_allele, alt_allele) {
                (Allele::Sequence(_), Allele::Deletion) => {
                    // In-frame deletion: remove bases and compare ref/alt amino acids
                    let del_len = ref_allele.len();
                    let end = (cds_idx + del_len).min(alt_seq.len());
                    alt_seq.drain(cds_idx..end);

                    // Number of complete codons deleted
                    let del_codons = del_len / 3;

                    if codon_offset == 0 {
                        // Deletion starts at codon boundary: VEP shows deleted AAs vs "-"
                        let ref_end = (codon_start + del_codons * 3).min(seq_bytes.len());
                        let ref_region = &seq_bytes[codon_start..ref_end];
                        let ref_aas: String = ref_region
                            .chunks(3)
                            .filter(|c| c.len() == 3)
                            .map(|c| {
                                self.codon_table_for(transcript)
                                    .translate(&[c[0], c[1], c[2]])
                                    as char
                            })
                            .collect();
                        let ref_codons: String = ref_region
                            .iter()
                            .map(|&b| (b as char).to_uppercase().next().unwrap())
                            .collect();
                        let aa_pair = Some((ref_aas, "-".to_string()));
                        let codon_pair = Some((ref_codons, "-".to_string()));
                        return (aa_pair, codon_pair);
                    } else {
                        // Deletion within a codon: show affected codons ref and alt
                        let n_ref_codons = del_codons + 1;
                        let ref_end = (codon_start + n_ref_codons * 3).min(seq_bytes.len());
                        let ref_region = &seq_bytes[codon_start..ref_end];
                        let ref_aas: String = ref_region
                            .chunks(3)
                            .filter(|c| c.len() == 3)
                            .map(|c| {
                                self.codon_table_for(transcript)
                                    .translate(&[c[0], c[1], c[2]])
                                    as char
                            })
                            .collect();
                        let alt_codon_end = (codon_start + 3).min(alt_seq.len());
                        let alt_region = &alt_seq[codon_start..alt_codon_end];
                        let alt_aas: String = if alt_region.len() == 3 {
                            String::from(self.codon_table_for(transcript).translate(&[
                                alt_region[0],
                                alt_region[1],
                                alt_region[2],
                            ]) as char)
                        } else {
                            "-".to_string()
                        };
                        let ref_codons: String = ref_region
                            .iter()
                            .map(|&b| (b as char).to_uppercase().next().unwrap())
                            .collect();
                        let alt_codons: String = if alt_aas == "-" {
                            "-".to_string()
                        } else {
                            alt_region
                                .iter()
                                .map(|&b| (b as char).to_uppercase().next().unwrap())
                                .collect()
                        };
                        let aa_pair = Some((ref_aas, alt_aas));
                        let codon_pair = Some((ref_codons, alt_codons));
                        return (aa_pair, codon_pair);
                    }
                }
                (Allele::Deletion, Allele::Sequence(ins_bases)) => {
                    // In-frame insertion: reverse-complement for reverse strand
                    let mut bases: Vec<u8> = ins_bases.clone();
                    if transcript.strand == Strand::Reverse {
                        bases = bases.iter().rev().map(|&b| complement(b)).collect();
                    }
                    // For reverse strand, the VCF insertion point maps to one base
                    // earlier in CDS space, so shift the insertion index by 1
                    let ins_idx = if transcript.strand == Strand::Reverse {
                        cds_idx + 1
                    } else {
                        cds_idx
                    };
                    for (i, &b) in bases.iter().enumerate() {
                        if ins_idx + i <= alt_seq.len() {
                            alt_seq.insert(ins_idx + i, b);
                        }
                    }

                    // Ref: the single codon at the insertion point
                    let ref_codon_str: String = ref_codon
                        .iter()
                        .map(|&b| (b as char).to_lowercase().next().unwrap())
                        .collect();

                    // Alt: translate codons spanning the insertion
                    let ins_codons = (bases.len() / 3) + 1;
                    let alt_end = (codon_start + ins_codons * 3).min(alt_seq.len());
                    let alt_region = &alt_seq[codon_start..alt_end];
                    let alt_aas: String = alt_region
                        .chunks(3)
                        .filter(|c| c.len() == 3)
                        .map(|c| {
                            self.codon_table_for(transcript)
                                .translate(&[c[0], c[1], c[2]]) as char
                        })
                        .collect();

                    // Build alt codon string: original bases lowercase, inserted uppercase
                    let ins_offset_in_codon = if transcript.strand == Strand::Reverse {
                        codon_offset + 1
                    } else {
                        codon_offset
                    };
                    let mut alt_codon_display = String::new();
                    for (i, &b) in alt_region.iter().enumerate() {
                        let is_original = if i < ins_offset_in_codon {
                            true
                        } else {
                            i >= ins_offset_in_codon + bases.len()
                        };
                        if is_original {
                            alt_codon_display.push((b as char).to_lowercase().next().unwrap());
                        } else {
                            alt_codon_display.push((b as char).to_uppercase().next().unwrap());
                        }
                    }

                    let aa_pair = Some((ref_aa_str, alt_aas));
                    let codon_pair = Some((ref_codon_str, alt_codon_display));
                    return (aa_pair, codon_pair);
                }
                _ => {}
            }
            (Some((ref_aa_str, "X".to_string())), None)
        }
    }

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
}
