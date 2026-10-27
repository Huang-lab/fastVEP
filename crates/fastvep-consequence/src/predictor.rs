use fastvep_core::{
    Allele, Consequence, GenomicPosition, Impact, Strand, TranscriptLocation, TranscriptSpan,
};
use fastvep_genome::{is_mitochondrial, mitochondrial_codon_table, CodonTable, Transcript};
use std::sync::Arc;

use crate::edits::{oriented, Edit, EditTests};
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
    /// The `cDNA_position`, `CDS_position`, `Protein_position`, `EXON` and `INTRON` cells as
    /// Ensembl VEP writes them, which differ from the fields above wherever the change has
    /// an end that is not in an exon (`?-7`) or is an insertion (`601-602`); see
    /// [`Transcript::location`]. The fields above stay what the criteria read: the
    /// strand-ordered coordinates the consequence and HGVS code work in.
    pub location: TranscriptLocation,
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
/// The codon table cannot see readthrough. A selenoprotein's in-frame UGA translates to `*` here,
/// while `build_sequences` has already resolved it to selenocysteine from the annotated CDS
/// extent, so a substitution at that codon reads as `stop_lost` and a synonymous change reads as
/// `stop_retained`.
///
/// The peptide is consulted only when the table says terminator and the peptide disagrees, which
/// is exactly the readthrough case; every other position still decides by the codon table alone.
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

/// The reference residue of the first codon, where the annotation reads the initiator as Met
/// though its codon is not `ATG` (`Transcript::initiator_met`: Ensembl's `initial_met` edit,
/// RefSeq's `transl_except ... aa:Met`).
///
/// Ensembl VEP applies the edit to the *reference* peptide of the codon
/// (`TranscriptVariationAllele::peptide`, under `is_reference`) and to no other, so `Amino_acids`
/// reads `M/L` for a CTG start and a change to the codon is a `start_lost`; an alternate codon
/// is translated plainly, even one that is the same codon. Only codon 1 is recoded.
fn resolve_initiator_residue(transcript: &Transcript, codon_index: usize, translated: u8) -> u8 {
    if codon_index == 0 && transcript.initiator_met && translated != b'M' {
        b'M'
    } else {
        translated
    }
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

/// Ensembl's `within_cds`: some piece of the change lies in the coding sequence.
///
/// An end with a CDS coordinate is such a piece. A change that starts and ends outside the CDS
/// can still cover some of it, though: a deletion of a whole exon, from the intron on one side to
/// beyond the transcript on the other, has coding bases in the middle, and for Ensembl that is
/// `within_cds` and so `coding_sequence_variant`. An insertion has no middle.
fn within_cds(
    transcript: &Transcript,
    location: &TranscriptLocation,
    var_start: u64,
    var_end: u64,
) -> bool {
    if location.cds.start.is_some() || location.cds.end.is_some() {
        return true;
    }
    let (Some(coding_start), Some(coding_end)) =
        (transcript.coding_region_start, transcript.coding_region_end)
    else {
        return false;
    };
    var_start <= var_end
        && transcript.exons.iter().any(|exon| {
            let (lo, hi) = (exon.start.max(coding_start), exon.end.min(coding_end));
            lo <= hi && var_start <= hi && var_end >= lo
        })
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
/// `consequences` is every SO term that holds, in rank order, and it can be empty: Ensembl
/// evaluates each predicate in `VariationEffect.pm` independently and keeps all that hold, so one
/// delins can be both `stop_gained` and `protein_altering_variant` - 164 of the 1,803 coding rows
/// the ClinVar 2-star in-frame delins produce are exactly that pair - and one that replaces the
/// initiator, introduces a stop and shifts the frame fires `stop_gained`, `frameshift_variant`
/// *and* `start_lost`. A slot for "the" term, and then for a second, dropped the lowest-ranked of
/// those without a diagnostic.
pub struct CodingChange {
    pub consequences: Vec<Consequence>,
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

/// The two peptides Ensembl's `hgvs_protein` starts from, for one allele at one location.
///
/// Both are the window of codons the change touches, the reference's and the allele's, with `-`
/// for a window that holds no residue (an insertion between two codons has no reference residue).
pub struct ProteinWindow {
    pub ref_peptide: String,
    pub alt_peptide: String,
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
                location: TranscriptLocation::default(),
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
                location: TranscriptLocation::default(),
                distance,
                protein_length: None,
                escapes_nmd: None,
            };
        }

        // Where the change falls on this transcript, as Ensembl's mapper reports it. The coding
        // terms are questions put to those coordinates, so they are worked out once, here.
        let location = transcript.location(var_start, var_end);

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
            //
            // `within_cdna`, which gates the UTR terms, has a second way to hold that
            // `within_cds` has for the CDS: a change inside a frameshift intron is within the
            // transcript, so one in a short intron of the UTR is a UTR variant. It was
            // `coding_transcript_variant` here, a term Ensembl gives a change that covers a whole
            // transcript.
            let in_cdna = in_exon || splice.within_frameshift_intron;
            let hits_coding =
                in_exon && self.is_in_coding_region(var_start, var_end, coding_start, coding_end);
            let hits_5_utr = in_cdna && self.overlaps_5_utr(var_start, var_end, transcript);
            let hits_3_utr = in_cdna && self.overlaps_3_utr(var_start, var_end, transcript);

            if hits_5_utr {
                consequences.push(Consequence::FivePrimeUtrVariant);
            }
            if hits_3_utr {
                consequences.push(Consequence::ThreePrimeUtrVariant);
            }
            if hits_coding {
                // Coding exonic variant - determine coding consequence
                let change = self.predict_coding_consequence(
                    ref_allele, alt_allele, transcript, &location, var_start, var_end,
                );
                consequences.extend(change.consequences);
                amino_acids = change.amino_acids;
                codons = change.codons;
            } else if splice.within_frameshift_intron
                && var_start <= coding_end
                && var_end >= coding_start
            {
                // `within_cds` has a second way to hold: a change inside a frameshift intron
                // that lies within the coding region (`VariationEffect.pm`, "we also need to
                // check if the vf is in a frameshift intron within the CDS"). Such an intron
                // was skipped by the splice terms, and the change has no CDS coordinate and so
                // no peptide, which is all `coding_unknown` asks for. Without this the row fell
                // to the catch-all below and was `coding_transcript_variant`, a term Ensembl
                // gives only a change that covers the whole transcript.
                consequences.push(Consequence::CodingSequenceVariant);
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
            // Non-coding transcript. `non_coding_exon_variant` is a plain `overlap` of the
            // variant with an exon, so an insertion has to have both of the bases it sits between
            // inside it; one on an exon's first base has only one, and Ensembl calls it a
            // `non_coding_transcript_variant` (`within_non_coding_gene` is whatever remains).
            if location.exon.is_some() {
                consequences.push(Consequence::NonCodingTranscriptExonVariant);
            } else if in_intron || in_exon {
                consequences.push(Consequence::NonCodingTranscriptVariant);
            }
            if splice.intronic {
                consequences.push(Consequence::IntronVariant);
            }
        }

        // A change whose span covers the whole transcript (`complete_overlap_feature`) is read
        // differently from one that reaches part of it, in three ways (`VariationEffect.pm`):
        //
        // - a deletion (`deletion`: the alternate is shorter than the reference, `_get_alleles`)
        //   is `transcript_ablation` (`feature_ablation`), and nothing else: every other term's
        //   `include` excludes a complete overlap, so the 5' UTR, 3' UTR, `stop_retained_variant`
        //   and the splice terms the span reaches on the way are not said, and the positions,
        //   `EXON` and `HGVSc` are written as for any other change;
        // - any other change is not `coding_unknown` (it returns 0 on a complete overlap) and so
        //   not a `coding_sequence_variant`, but is a `coding_transcript_variant` beside whatever
        //   else it earned, if the transcript has a protein;
        // - and it is not a `non_coding_exon_variant` either, so on a transcript without a protein
        //   it is `non_coding_transcript_variant`.
        //
        // This predictor used to ask the regions the span reaches one by one: a deletion of all of
        // a 21-base transcript was `stop_retained_variant,5_prime_UTR_variant,3_prime_UTR_variant`,
        // LOW, a non-coding transcript the same change covered was MODIFIER, and only the
        // structural-variant predictor said `transcript_ablation` at all.
        if var_start <= tr_start && var_end >= tr_end {
            let deletion = matches!(
                (ref_allele.known_length(), alt_allele.known_length()),
                (Some(r), Some(a)) if a < r
            );
            if deletion {
                consequences.clear();
                consequences.push(Consequence::TranscriptAblation);
                amino_acids = None;
                codons = None;
            } else if transcript.is_coding() {
                consequences.retain(|c| *c != Consequence::CodingSequenceVariant);
                consequences.push(Consequence::CodingTranscriptVariant);
            } else if let Some(at) = consequences
                .iter()
                .position(|c| *c == Consequence::NonCodingTranscriptExonVariant)
            {
                consequences[at] = Consequence::NonCodingTranscriptVariant;
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
            location,
            distance,
            protein_length,
            escapes_nmd,
        }
    }

    /// Predict the coding consequence (missense, synonymous, frameshift, etc.)
    fn predict_coding_consequence(
        &self,
        ref_allele: &Allele,
        alt_allele: &Allele,
        transcript: &Transcript,
        location: &TranscriptLocation,
        var_start: u64,
        var_end: u64,
    ) -> CodingChange {
        // Ensembl translates an allele only when it and the reference are plain bases (see
        // `Allele::is_plain_bases`), and every test of the start and the stop codon reads
        // the edited sequence or the peptide, so for one that is not there is none to make.
        let plain = ref_allele.is_plain_bases() && alt_allele.is_plain_bases();
        // `ref_length` of `_pre_consequence_predicates`: the span of the change in the genome,
        // `0` for an insertion (`end = start - 1`). The length classes (`snp`, `insertion`,
        // `deletion`) are decided by it and not by the reference allele or the CDS span.
        let ref_length = (var_end as i64 - var_start as i64 + 1).max(0) as usize;
        let edits = if plain {
            EditTests::compute(&Edit {
                transcript,
                table: self.codon_table_for(transcript),
                location,
                start: var_start,
                end: var_end,
                alt: alt_allele.as_bytes(),
                ref_len: ref_length,
                alt_len: alt_allele.len(),
            })
        } else {
            EditTests::default()
        };
        let changes_length = ref_length != alt_allele.len();

        // A change with an end that has no CDS coordinate - it reaches into the UTR, or has a
        // flank in the intron beside an exon - is not a codon edit: the bases outside the CDS
        // belong to no codon, and the frame arithmetic would count them anyway, which is how a
        // delins over FOXL2's stop codon came back `frameshift_variant`. Ensembl reaches the same
        // conclusion structurally: `codon` returns nothing unless `cds_start` and `cds_end` are
        // both defined, so there is no peptide and every term that compares one is out; what is
        // left is what the edited sequence says about the start and the stop codon, and
        // `coding_sequence_variant` for a change that has a base in the CDS.
        if location.cds.start.is_none() || location.cds.end.is_none() {
            return self.terms_without_peptides(
                transcript,
                &edits,
                changes_length,
                within_cds(transcript, location, var_start, var_end),
                location,
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
            self.build_codon_window(transcript, alt_allele, location.cds, Some(ref_length))
        else {
            // The CDS has no sequence loaded, so there is no window to read
            return self.terms_without_peptides(transcript, &edits, changes_length, true, location);
        };

        // An allele that is not plain bases has no peptide in Ensembl (`peptide` is undef
        // unless `seq_is_unambiguous_dna`), so no `Amino_acids` and no term about the protein.
        // What is left is `frameshift`, which only needs lengths (`seq_length` counts the
        // ambiguity codes) and declines on a partial codon, and `coding_unknown`'s
        // `coding_sequence_variant` for everything else. An `N` was a `missense_variant`
        // with `p.Lys31Xaa`, and a `<*>` beside a real ALT a `frameshift_variant`.
        if !plain {
            let frameshift = !window.partial_codon
                && alt_allele
                    .known_length()
                    .is_some_and(|len| len.abs_diff(window.ref_len) % 3 != 0);
            let mut consequences = vec![if frameshift {
                Consequence::FrameshiftVariant
            } else {
                Consequence::CodingSequenceVariant
            }];
            if window.partial_codon {
                consequences.push(Consequence::IncompleteTerminalCodonVariant);
            }
            return CodingChange {
                consequences,
                amino_acids: None,
                codons: Some((window.ref_codons.clone(), window.alt_codons.clone())),
            };
        }

        // Ensembl's `_peptide` is the protein without its terminator; ours
        // carries it, because that is what the annotation's own translation
        // ends with.
        let peptide = transcript
            .peptide
            .as_deref()
            .map(|p| p.strip_suffix('*').unwrap_or(p));
        self.terms_for_window(&window, &edits, peptide, location.cds, ref_length)
    }

    /// The SO terms Ensembl derives from a codon window.
    ///
    /// A port of `VariationEffect.pm` and of the `include` pre-filters of `Config.pm` in release
    /// 116 (read in full, not recalled): that is not the graph 115 had, which is why this does not
    /// look like what stood here. `stop_lost` and `stop_retained` consult the edited sequence
    /// when the alternate peptide holds an `X`; `ref_eq_alt_sequence` asks that what follows the
    /// edited protein *opens with a stop* and no longer that it is shorter than three residues;
    /// `frameshift` declines when the change starts on the terminator; an insertion or a deletion
    /// of a terminator is neither `inframe_insertion` nor `inframe_deletion`.
    ///
    /// Two things make it a graph and not a ladder. The predicates call one another, and the
    /// order is not a severity order: `stop_retained` gates both `frameshift` and `stop_gained`,
    /// `start_lost` gates `inframe_insertion` and `protein_altering_variant`, and
    /// `protein_altering_variant` defers to almost everything. And every term has an `include`
    /// that decides whether its predicate is asked at all: `frameshift_variant` only of a change
    /// that moves the length (`snp => 0`), `inframe_insertion` and `inframe_deletion` only of one
    /// that adds or removes bases, `missense_variant` only of one that does neither. The predicate
    /// calls inside other predicates do not go through the `include`, so each is computed in its
    /// raw form here and gated where it becomes a term; without the gate a frameshift whose
    /// alternate peptide is as long as the reference (`K/X`) is a `missense_variant`.
    ///
    /// `e` carries the answers to the questions Ensembl puts to the edited transcript, which are
    /// coordinate and sequence questions and not peptide ones; `cds` is Ensembl's CDS position
    /// of the change, `?` where an end is not coding.
    fn terms_for_window(
        &self,
        w: &CodonWindow,
        e: &EditTests,
        peptide: Option<&str>,
        cds: TranscriptSpan,
        ref_length: usize,
    ) -> CodingChange {
        let (r, a) = (w.ref_aas.as_str(), w.alt_aas.as_str());
        let partial = w.partial_codon;
        // `snp`, `insertion` and `deletion` of `_pre_consequence_predicates`: the lengths of the
        // two alleles in the genome, which the codon strings carry no sign of when the change is a
        // frameshift and which the CDS span does not give for a change that crosses an intron
        let changes_length = ref_length != w.alt_len;
        let (insertion, deletion) = (w.alt_len > ref_length, ref_length > w.alt_len);
        let starts_or_ends_with_ref = |pep: &str| pep.starts_with(r) || pep.ends_with(r);

        // `stop_lost`: the peptides say so, unless the alternate holds an `X` (a frameshift's
        // unfinished codon), when only the edited sequence can
        let stop_lost = !partial
            && if a.contains('X') {
                e.ins_del_stop_altered
            } else {
                !a.contains('*') && r.contains('*')
            };

        // `ref_eq_alt_sequence`. Its first clause asks whether the reference protein with the
        // alternate residues put in its place, cut back to the reference's own length, *is* the
        // reference protein and what spilled past that length opens with a terminator: an
        // insertion of a stop just before the stop. 115 asked only that the spill be shorter than
        // three residues, which called a 4-base frameshift near the end of a protein an
        // in-frame insertion that retains the stop. The edited protein is compared in place and
        // not built: this runs once per (variant x transcript x allele), and a `String` here
        // is a String per row.
        let ref_eq_alt_sequence = |peptide: &str| -> bool {
            if r == "X" && a == "X" {
                return false; // an incomplete terminal codon on both sides
            }
            let pep = peptide.as_bytes();
            let (len, alt) = (pep.len(), a.as_bytes());
            // on the terminator, one past the protein's last residue, with a replacement that
            // opens with a stop
            if w.tl_start > len && alt.first() == Some(&b'*') {
                return true;
            }
            // both peptides hold a terminator at the same place
            if r.contains('*') && r.find('*') == a.find('*') {
                return true;
            }
            if w.tl_start == 0 {
                return false;
            }
            // `substr($mut_seq, $tl_start - 1, $tl_end - $tl_start + 1) = $alt_pep`, which clips
            // at the end: a window sitting *on* the terminator appends
            let before = (w.tl_start - 1).min(len);
            let after = (before + (w.tl_end + 1).saturating_sub(w.tl_start)).min(len);
            let edited_len = before + alt.len() + (len - after);
            if edited_len <= len {
                return false; // nothing spilled past the reference's length
            }
            let at = |i: usize| -> Option<u8> {
                if i < before {
                    pep.get(i).copied()
                } else if i < before + alt.len() {
                    Some(alt[i - before])
                } else {
                    pep.get(after + (i - before - alt.len())).copied()
                }
            };
            at(len) == Some(b'*') && (0..len).all(|i| at(i) == Some(pep[i]))
        };

        // `stop_retained`: a peptide that is whole is judged by the comparison above, one with an
        // `X` or none at all by the edited sequence, with the length of an insertion counted
        let stop_retained = !partial
            && !stop_lost
            && if !a.is_empty() && !a.contains('X') {
                peptide.is_some_and(ref_eq_alt_sequence)
            } else {
                changes_length && e.overlaps_stop_cil && !e.ins_del_stop_altered_cil
            };

        // `frameshift` needs both ends of the change in the CDS, declines on a partial codon, on
        // `stop_retained`, and on a change that starts on the terminator (it does not move the
        // frame there); otherwise it is the length arithmetic on the CDS span
        let cds_span = match (cds.start, cds.end) {
            (Some(s), Some(en)) => Some(en as i64 - s as i64 + 1),
            _ => None,
        };
        let frameshift = !partial
            && !stop_retained
            && !r.starts_with('*')
            && cds_span.is_some_and(|var_len| (w.alt_len as i64 - var_len).rem_euclid(3) != 0);
        let stop_gained = !stop_retained && !stop_lost && a.contains('*') && !r.contains('*');

        // The start side. `start_retained_variant` asks the edited sequence whether the start
        // survived: for a change that keeps the length, whether the codon at the start is still
        // an ATG, and otherwise whether the whole coding sequence is still there behind whatever
        // was added (atg -> aCGAtg)
        let overlaps_initiator = e.overlaps_initiator;
        let start_retained = overlaps_initiator
            && if changes_length {
                !e.ins_del_start_altered
            } else {
                !e.snp_start_altered
            };

        // Ensembl evaluates `inframe_insertion` with `start_lost` still at its default of 0 while
        // `start_lost` is being decided, and `inframe_insertion` itself declines once
        // `start_lost` holds, so the pair is decided in two steps. `_get_codon_alleles` returns
        // nothing for a frameshift, and neither in-frame term can hold without codons.
        let inframe_insertion_before_start_lost = !frameshift
            && w.alt_window.len() > w.ref_window.len()
            && !(start_retained && a.ends_with(r))
            && !(r == "*" && a == "*")
            && {
                // `inframe_insertion` cuts the alternate peptide back to its first terminator
                // before the prefix/suffix test, so an insertion that also introduces a stop
                // is still an insertion
                let to_stop = match a.find('*') {
                    Some(i) => &a[..=i],
                    None => a,
                };
                starts_or_ends_with_ref(to_stop)
            };
        let inframe_deletion =
            !partial && !frameshift && w.alt_window.len() < w.ref_window.len() && r != "*" && {
                // the codon strings, not the peptides: the replacement counts as a deletion
                // only when it reproduces a prefix, a suffix, or a whole-codon interior slice
                // of what it replaced
                let interior = {
                    let (rr, aa) = trim_common_ends(&w.ref_window, &w.alt_window);
                    aa.is_empty() && rr.len().is_multiple_of(3)
                };
                w.ref_window.starts_with(&w.alt_window)
                    || w.ref_window.ends_with(&w.alt_window)
                    || interior
            };

        // `start_lost`: three ways. An edit that alters the initiator and is neither in-frame
        // term (an in-frame deletion that removes it falls through to the next two); the first
        // codon after the 5' UTR no longer reading ATG, which is judged for every kind of change
        // and not only inversions; or the peptide test, which Ensembl skips when either side is
        // missing or the alternate is the frameshift placeholder.
        //
        // None of the three holds when the window still opens with its Met after the edit: a
        // protein that begins with Met has not lost its initiator. Ensembl's peptide test (the
        // third way) cannot see that for a deletion, whose alternate is shorter than the reference
        // it is compared with, so the loss of residue 2 (`MA/M`) is
        // `start_lost,inframe_deletion,start_retained_variant`, a pair that contradicts itself,
        // and `p.Ala2?`, which says residue 2 is unknown when it is gone. Only a Met counts: a
        // transcript that starts on `CTG` has a Leu for its first residue, and any change to that
        // codon, `CTG` to `CTA` as much as the rest, is a start loss in the cache of Ensembl
        // itself.
        let initiator_kept = w.tl_start == 1 && r.starts_with('M') && a.starts_with('M');
        let start_lost = overlaps_initiator
            && !initiator_kept
            && ((e.ins_del_start_altered
                && !(inframe_insertion_before_start_lost || inframe_deletion))
                || e.inv_start_altered
                || (!r.is_empty()
                    && !a.is_empty()
                    && a != "X"
                    && w.tl_start == 1
                    && !a.ends_with(r)
                    && !a.starts_with(r)));
        let inframe_insertion = inframe_insertion_before_start_lost && !start_lost;

        let protein_altering = r.len() != a.len()
            && !r.starts_with('*')
            && !a.starts_with('*')
            && !starts_or_ends_with_ref(a)
            && !inframe_deletion
            && !start_lost
            && !frameshift;

        // `missense_variant`, `synonymous_variant` and `coding_sequence_variant`, each as Ensembl
        // asks it. An unknown residue silences `synonymous_variant` (it tests both peptides for an
        // `X`) but not `missense_variant`, and `coding_unknown` does not defer to either: so a
        // window that resolves one residue and not the next earns *both* `missense_variant` and
        // `coding_sequence_variant` (`QX/HX` on ENST00000515665).
        let unresolved = r.contains('X') || a.contains('X');
        let missense = r != a
            && r.len() == a.len()
            && !start_lost
            && !stop_lost
            && !stop_gained
            && !partial
            && !stop_retained;
        let synonymous =
            !r.is_empty() && a == r && !stop_retained && !start_retained && !unresolved;
        let coding_unknown = unresolved
            && !(frameshift
                || inframe_deletion
                || protein_altering
                || start_retained
                || start_lost
                || stop_retained
                || stop_lost);

        // The `include` of each term, then every term that holds, in rank order
        let consequences = [
            (stop_gained, Consequence::StopGained),
            (changes_length && frameshift, Consequence::FrameshiftVariant),
            (stop_lost, Consequence::StopLost),
            (start_lost, Consequence::StartLost),
            (
                insertion && inframe_insertion,
                Consequence::InframeInsertion,
            ),
            (deletion && inframe_deletion, Consequence::InframeDeletion),
            (!changes_length && missense, Consequence::MissenseVariant),
            (protein_altering, Consequence::ProteinAlteringVariant),
            (partial, Consequence::IncompleteTerminalCodonVariant),
            (start_retained, Consequence::StartRetainedVariant),
            (stop_retained, Consequence::StopRetainedVariant),
            (synonymous, Consequence::SynonymousVariant),
            (coding_unknown, Consequence::CodingSequenceVariant),
        ]
        .into_iter()
        .filter_map(|(fired, term)| fired.then_some(term))
        .collect();

        // `-` is how an absent side is written, for residues as for codons.
        let shown = |pep: &str| {
            if pep.is_empty() {
                "-".to_string()
            } else {
                pep.to_string()
            }
        };
        CodingChange {
            consequences,
            amino_acids: Some((shown(&w.ref_aas), shown(&w.alt_aas))),
            codons: Some((w.ref_codons.clone(), w.alt_codons.clone())),
        }
    }

    /// The terms of a change in the coding sequence that has no peptide, because one end of it has
    /// no CDS coordinate (it reaches into the UTR or the intron beside an exon): only the
    /// questions put to the edited sequence can be answered, and everything else Ensembl asks
    /// needs a peptide. What is left is `start_lost`, `start_retained_variant`, `stop_lost`,
    /// `stop_retained_variant`, `incomplete_terminal_codon_variant` and `coding_sequence_variant`,
    /// and `_ins_del_stop_altered_cil` giving up when the change has no CDS start is what makes a
    /// deletion that starts in the 3' UTR and covers the terminator `stop_retained_variant` and
    /// not `stop_lost`.
    ///
    /// `coding_sequence_variant` is `coding_unknown`: it needs a base inside the CDS (`within_cds`,
    /// which an insertion in the UTR next to the start codon does not have) and defers to every
    /// start and stop term, but not to `incomplete_terminal_codon_variant`.
    fn terms_without_peptides(
        &self,
        transcript: &Transcript,
        e: &EditTests,
        changes_length: bool,
        within_cds: bool,
        location: &TranscriptLocation,
    ) -> CodingChange {
        // `partial_codon` is asked of the first residue the change touches, when it has one, and
        // both stop predicates decline when it holds
        let partial = location.protein.start.is_some_and(|p| {
            transcript
                .translateable_seq
                .as_deref()
                .is_some_and(|seq| p >= 1 && seq.len().saturating_sub((p as usize - 1) * 3) < 3)
        });
        let stop_lost = !partial && e.ins_del_stop_altered;
        let stop_retained = !partial
            && !stop_lost
            && changes_length
            && e.overlaps_stop_cil
            && !e.ins_del_stop_altered_cil;
        let start_lost = e.overlaps_initiator && (e.ins_del_start_altered || e.inv_start_altered);
        let start_retained = e.overlaps_initiator
            && if changes_length {
                !e.ins_del_start_altered
            } else {
                !e.snp_start_altered
            };
        let coding_unknown =
            within_cds && !(start_retained || start_lost || stop_retained || stop_lost);
        let consequences = [
            (stop_lost, Consequence::StopLost),
            (start_lost, Consequence::StartLost),
            (partial, Consequence::IncompleteTerminalCodonVariant),
            (start_retained, Consequence::StartRetainedVariant),
            (stop_retained, Consequence::StopRetainedVariant),
            (coding_unknown, Consequence::CodingSequenceVariant),
        ]
        .into_iter()
        .filter_map(|(fired, term)| fired.then_some(term))
        .collect();
        CodingChange {
            consequences,
            amino_acids: None,
            codons: None,
        }
    }

    /// Build the codon window of a change from the CDS span Ensembl's mapper gives it.
    ///
    /// The window is every codon the span touches. For an insertion the span covers no bases at
    /// all (`cds` is the pair `(c + 1, c)`), so the window is the single codon it falls inside -
    /// and *nothing*, when it falls exactly on a codon boundary, which is why VEP writes
    /// `-/HENKTKGD` for a codon-aligned insertion rather than repeating the flanking residue on
    /// both sides.
    ///
    /// The span is counted in the CDS (`cds_end - cds_start + 1`), not in the genome, and the whole
    /// alternate allele replaces it. That is what Ensembl does (`vf_nt_len` in `codon`), and it
    /// matters for a change that crosses an intron: a deletion of 17 CDS bases and the 600 intron
    /// bases between them is a 17-base deletion of the coding sequence, a frameshift, and the
    /// window is the six codons around it (`LLLPSI/X`). A change that is not in one piece was
    /// given no residues at all here before, and so `coding_sequence_variant`.
    ///
    /// `ref_bases` is how many bases of the reference codons are written in uppercase, which is
    /// the length of the reference allele *in the genome* (`display_codon` uppercases
    /// `length($self->feature_seq)` from the codon position on), and not the CDS span the codons
    /// are cut with. They differ for a change with an intron between its ends: a deletion of
    /// `CACAC` over a one-base intron replaces four bases of the CDS and Ensembl writes
    /// `CACCGc`, five. `None` is the CDS span, for a caller that does not show the codons.
    ///
    /// Returns `None` when the CDS has no sequence loaded or the span has no CDS coordinate.
    fn build_codon_window(
        &self,
        transcript: &Transcript,
        alt_allele: &Allele,
        cds: TranscriptSpan,
        ref_bases: Option<usize>,
    ) -> Option<CodonWindow> {
        let seq = transcript.translateable_seq.as_deref()?.as_bytes();
        let (cds_start, cds_end) = (cds.start?, cds.end?);
        if cds_start == 0 || cds_end == 0 {
            return None;
        }

        // `first` is the 0-based CDS offset the change starts at; for an insertion it is the
        // offset it is inserted *before*, so a change of no reference bases still has a position.
        let ref_len = (cds_end as i64 - cds_start as i64 + 1).max(0) as usize;
        let first = if ref_len == 0 {
            cds_end as usize
        } else {
            (cds_start - 1) as usize
        };

        // Alt bases in transcript orientation (reversed as well as complemented on the reverse
        // strand; see `oriented`)
        let alt_cds: Vec<u8> = oriented(alt_allele.as_bytes(), transcript.strand).into_owned();

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
        let full_win_end = win_end;
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
        // Ensembl reads the alternate codon out of the edited CDS with the 3' UTR appended
        // (`_get_alternate_cds`) and the reference codon out of the CDS alone, so a change in a
        // last codon that the CDS leaves short completes it from the UTR on the alternate side
        // only: `gA/gCt` and `X/A` for ENST00000512526 `17:50167975 A>C`, where the reference
        // codon is the two bases the CDS has and the alternate is three.
        let missing = full_win_end - win_end;
        if missing > 0 {
            if let (Some(spliced), Some(coding_end)) = (
                transcript.spliced_seq.as_deref(),
                transcript.cdna_coding_end,
            ) {
                let utr = spliced.as_bytes().get(coding_end as usize..).unwrap_or(&[]);
                alt_window.extend_from_slice(&utr[..missing.min(utr.len())]);
            }
        }

        let table = self.codon_table_for(transcript);
        // A codon the edit happens to leave byte-identical to the reference
        // codon at the same offset keeps whatever the reference resolved to, so
        // a selenoprotein's readthrough UGA is not re-read as a new terminator.
        // Every other codon is translated plainly, which is what makes a
        // genuinely introduced stop a stop.
        let translate = |window: &[u8], anchored: bool| -> String {
            let mut pep: String = window
                .as_chunks::<3>()
                .0
                .iter()
                .enumerate()
                .map(|(i, codon)| {
                    let translated = table.translate(codon);
                    let index = win_start / 3 + i;
                    if anchored {
                        let residue = resolve_readthrough_residue(transcript, index, translated);
                        resolve_initiator_residue(transcript, index, residue) as char
                    } else if ref_window.get(i * 3..i * 3 + 3) == Some(&codon[..]) {
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
            ref_codons: render(&ref_window, ref_bases.unwrap_or(ref_len)),
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

    /// The two peptides `hgvs_protein` starts from for `alt_allele` at `location`.
    ///
    /// This is the window the consequence code reads, and for the same reason: Ensembl's HGVSp
    /// is built from the peptides of the *same* codon window `Amino_acids` is, only at the
    /// location of the 3'-shifted change, so the caller passes the location it shifted to.
    ///
    /// `None` when the CDS has no sequence loaded or the span has no CDS coordinate.
    pub fn protein_window(
        &self,
        transcript: &Transcript,
        alt_allele: &Allele,
        location: &TranscriptLocation,
    ) -> Option<ProteinWindow> {
        let w = self.build_codon_window(transcript, alt_allele, location.cds, None)?;
        let or_dash = |aas: String| if aas.is_empty() { "-".to_string() } else { aas };
        Some(ProteinWindow {
            ref_peptide: or_dash(w.ref_aas),
            alt_peptide: or_dash(w.alt_aas),
        })
    }

    /// The translation of the CDS with `alt_allele` in place of the CDS span `cds` and the 3' UTR
    /// after it, which is `_get_alternate_cds()->translate()` of Ensembl: what a frameshift or
    /// a lost terminator is counted against. Only whole codons are translated, as BioPerl does.
    ///
    /// The span is the one the mapper gives the change, counted in the CDS, and the whole allele
    /// replaces it (see [`Self::build_codon_window`]). `None` when the CDS has no sequence or the
    /// span has no CDS coordinate.
    pub fn edited_translation(
        &self,
        transcript: &Transcript,
        alt_allele: &Allele,
        cds: TranscriptSpan,
    ) -> Option<String> {
        let seq = transcript.translateable_seq.as_deref()?.as_bytes();
        let (cds_start, cds_end) = (cds.start?, cds.end?);
        if cds_start == 0 || cds_end == 0 {
            return None;
        }
        let ref_len = (cds_end as i64 - cds_start as i64 + 1).max(0) as usize;
        let first = if ref_len == 0 {
            cds_end as usize
        } else {
            (cds_start - 1) as usize
        };
        let alt_cds = oriented(alt_allele.as_bytes(), transcript.strand);
        // A transcript whose spliced sequence ends before its own coding region does is not one
        // to describe a frameshift of: what follows the change would be read from nothing.
        let utr: &[u8] = match (
            transcript.spliced_seq.as_deref(),
            transcript.cdna_coding_end,
        ) {
            (Some(spliced), Some(coding_end)) => spliced.as_bytes().get(coding_end as usize..)?,
            _ => &[],
        };
        let mut edited = Vec::with_capacity(seq.len() + utr.len() + alt_cds.len());
        edited.extend_from_slice(&seq[..first.min(seq.len())]);
        edited.extend_from_slice(&alt_cds);
        edited.extend_from_slice(seq.get(first + ref_len..).unwrap_or(&[]));
        edited.extend_from_slice(utr);

        let table = self.codon_table_for(transcript);
        Some(
            edited
                .as_chunks::<3>()
                .0
                .iter()
                .map(|codon| table.translate(codon) as char)
                .collect(),
        )
    }

    /// Distance from the change to the nearest end of the transcript, `0` when it overlaps.
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

    /// Whether the change reaches the genomic extent of the coding region, introns between the
    /// coding exons included. An insertion counts as the two bases it sits between (`min_vf` and
    /// `max_vf` of `_bvfo_preds`), which is what opens the `coding` pre-filter for one on the
    /// first coding base: the predicates behind it then decline for their own reasons, and the
    /// UTR term it gets there is not a substitute for them.
    fn is_in_coding_region(
        &self,
        var_start: u64,
        var_end: u64,
        coding_start: u64,
        coding_end: u64,
    ) -> bool {
        let (lo, hi) = (var_start.min(var_end), var_start.max(var_end));
        lo <= coding_end && hi >= coding_start
    }

    /// Ensembl's `_before_coding`: does the change overlap the genomic interval from the
    /// transcript's first base to the one before the coding region - the UTR on the low side?
    ///
    /// Three details of the Perl are kept, because rows depend on each of them.
    ///
    /// An insertion on the first coding base is special-cased (`bvf_s == bvf_e + 1 && bvf_s ==
    /// cds_s`): the interval test would ask that both bases it sits between are in the UTR, and
    /// one is not, yet Ensembl calls it a UTR variant.
    ///
    /// The interval is genomic, so it includes any intron lying between the transcript's edge and
    /// the coding region, and a variant whose only bases there are intronic still counts. That is
    /// Ensembl's rule rather than an oversight - the term is gated only on the variant touching an
    /// exon somewhere - and matching it keeps the UTR terms comparable with VEP's. It can only add
    /// a MODIFIER term next to a correctly-derived one, never change the reported impact.
    ///
    /// And `overlap` does not ask that the interval hold a base. For a transcript whose coding
    /// region reaches its own edge the interval is `[start, start - 1]`, and a change that spans
    /// that edge still overlaps it, so a deletion of a stop codon that runs on past the last base
    /// of the transcript is a `3_prime_UTR_variant`: nine of the 145 ClinVar rows whose
    /// consequence moved when the coding terms were ported are that. Guarding on "there is a UTR
    /// on this side", as this did, is what left the term off.
    fn before_coding(&self, var_start: u64, var_end: u64, transcript: &Transcript) -> bool {
        let Some(coding_start) = transcript.coding_region_start else {
            return false;
        };
        // `overlap($bvf_s, $bvf_e, $t_s, $cds_s - 1)` is `e >= t_start && s <= cds_start - 1`, which
        // is `s < cds_start`; i64 because an insertion is `start = end + 1`
        let (s, e) = (var_start as i64, var_end as i64);
        let (t_start, cds_start) = (transcript.start as i64, coding_start as i64);
        (s == e + 1 && s == cds_start) || (e >= t_start && s < cds_start)
    }

    /// Ensembl's `_after_coding`, the same on the high side: `[cds_end + 1, transcript.end]`.
    fn after_coding(&self, var_start: u64, var_end: u64, transcript: &Transcript) -> bool {
        let Some(coding_end) = transcript.coding_region_end else {
            return false;
        };
        let (s, e) = (var_start as i64, var_end as i64);
        let (t_end, cds_end) = (transcript.end as i64, coding_end as i64);
        (s == e + 1 && e == cds_end) || (e > cds_end && s <= t_end)
    }

    /// Does the variant's span reach the 5' UTR? Which genomic side that is
    /// depends on the strand; the overlap test itself does not.
    fn overlaps_5_utr(&self, var_start: u64, var_end: u64, transcript: &Transcript) -> bool {
        if transcript.coding_region_start.is_none() || transcript.coding_region_end.is_none() {
            return false;
        }
        match transcript.strand {
            Strand::Forward => self.before_coding(var_start, var_end, transcript),
            Strand::Reverse => self.after_coding(var_start, var_end, transcript),
        }
    }

    /// Does the variant's span reach the 3' UTR?
    fn overlaps_3_utr(&self, var_start: u64, var_end: u64, transcript: &Transcript) -> bool {
        if transcript.coding_region_start.is_none() || transcript.coding_region_end.is_none() {
            return false;
        }
        match transcript.strand {
            Strand::Forward => self.after_coding(var_start, var_end, transcript),
            Strand::Reverse => self.before_coding(var_start, var_end, transcript),
        }
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
    use crate::edits::complement;

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

    /// An in-frame deletion that removes the initiator is a start loss too.
    ///
    /// The terms of a hand-built window: `edits` is what the edited sequence says, `cds` the CDS
    /// span of the change, `ref_length` its length in the genome.
    fn window_terms(
        window: &CodonWindow,
        edits: EditTests,
        peptide: &str,
        cds: (u64, u64),
        ref_length: usize,
    ) -> Vec<Consequence> {
        ConsequencePredictor::default()
            .terms_for_window(
                window,
                &edits,
                Some(peptide),
                TranscriptSpan::new(Some(cds.0), Some(cds.1)),
                ref_length,
            )
            .consequences
    }

    /// Ensembl's `start_lost` (VariationEffect.pm release/116, l. 851) answers on the edited
    /// sequence first and falls through to the peptide test when that declines because the edit is
    /// an in-frame indel. Stopping at the first test reported KCNA2 `c.3_11del` on
    /// ENST00000639048 - `MNII/I`, the initiator among the residues it removes - as a MODERATE
    /// in-frame deletion where VEP 115.1 calls it `start_lost,inframe_deletion`, HIGH. Four rows
    /// over the 400-variant ClinVar in-frame deletion set, and all four are IMPACT.
    #[test]
    fn an_inframe_deletion_that_removes_the_initiator_is_a_start_loss() {
        let window = CodonWindow {
            ref_aas: "MNII".into(),
            alt_aas: "I".into(),
            ref_window: b"ATGAACATCATT".to_vec(),
            alt_window: b"ATT".to_vec(),
            ref_codons: "atGAACATCATt".into(),
            alt_codons: "att".into(),
            ref_len: 9,
            alt_len: 0,
            tl_start: 1,
            tl_end: 4,
            partial_codon: false,
        };
        // the edited sequence says the start was altered, and `start_lost` defers to the in-frame
        // deletion for that: it is the peptide test, which sees `I` where `MNII` was, that holds
        let edits = EditTests {
            overlaps_initiator: true,
            ins_del_start_altered: true,
            ..EditTests::default()
        };
        let terms = window_terms(&window, edits, "MNIIKEL", (3, 11), 9);
        assert_eq!(
            terms,
            vec![Consequence::StartLost, Consequence::InframeDeletion]
        );
    }

    /// Release 116 changed `ref_eq_alt_sequence`, the test behind `stop_retained_variant` for a
    /// peptide that is whole (l. 1326 of VariationEffect.pm; the 115 file asked only that what
    /// spilled past the protein's length was shorter than three residues). It now asks whether the
    /// protein with the alternate residues put in its place *is* the reference protein and what
    /// spilled opens with a terminator, or whether both peptides hold a terminator at the same
    /// place, or whether the window sits on the terminator and the replacement opens with one.
    ///
    /// So an insertion that merely repeats what follows it is an `inframe_insertion` and no longer
    /// also `stop_retained_variant` (MSH6 `c.4106_4108dup` on ENST00000936511 was both in 115.1),
    /// and one residue added in front of the stop is the same. The edge probe of 17,672 changes at
    /// every exon, CDS and transcript boundary agrees with VEP 116.3 on every row of that shape.
    #[test]
    fn stop_retained_follows_release_116_s_ref_eq_alt_sequence() {
        let peptide = "MSRQKEL";
        let on =
            |tl_start: usize, ref_aas: &str, alt_aas: &str, ref_codons: &str, alt_codons: &str| {
                CodonWindow {
                    ref_aas: ref_aas.into(),
                    alt_aas: alt_aas.into(),
                    ref_window: ref_codons.to_ascii_uppercase().into_bytes(),
                    alt_window: alt_codons.to_ascii_uppercase().into_bytes(),
                    ref_codons: ref_codons.into(),
                    alt_codons: alt_codons.into(),
                    ref_len: 0,
                    alt_len: alt_codons.len() - ref_codons.len(),
                    tl_start,
                    tl_end: tl_start + ref_codons.len() / 3 - 1,
                    partial_codon: false,
                }
            };
        let terms = |w: &CodonWindow, cds: (u64, u64)| {
            window_terms(w, EditTests::default(), peptide, cds, 0)
        };

        // A Leu in front of the Leu already there repeats what follows it: the protein is the
        // same and one residue longer, and what spilled is `L`, not a stop
        let repeat = CodonWindow {
            ref_aas: "E".into(),
            alt_aas: "EL".into(),
            ref_window: b"GAA".to_vec(),
            alt_window: b"GAATTA".to_vec(),
            ref_codons: "gaa".into(),
            alt_codons: "gaATTa".into(),
            ref_len: 3,
            alt_len: 6,
            tl_start: 6,
            tl_end: 6,
            partial_codon: false,
        };
        let got = window_terms(&repeat, EditTests::default(), peptide, (16, 18), 3);
        assert!(
            !got.contains(&Consequence::StopRetainedVariant),
            "release 116 does not retain the stop here: {got:?}"
        );

        // A whole stop codon inserted in front of the terminator: the window is empty and sits on
        // the terminator (one past the protein's last residue), and the replacement opens with a
        // stop. That is the first clause, and it holds.
        let stop_before_stop = on(8, "", "*", "", "TAA");
        let got = terms(&stop_before_stop, (21, 20));
        assert!(
            got.contains(&Consequence::StopRetainedVariant),
            "an inserted stop in front of the stop retains it: {got:?}"
        );

        // One residue in front of the stop is not a stop, so it does not
        let one_residue = on(8, "", "F", "", "TTT");
        let got = terms(&one_residue, (21, 20));
        assert!(
            !got.contains(&Consequence::StopRetainedVariant),
            "a residue in front of the stop does not retain it: {got:?}"
        );
        assert!(got.contains(&Consequence::InframeInsertion), "{got:?}");

        // Both peptides hold a terminator at the same place: the second clause
        let same_place = CodonWindow {
            ref_aas: "*".into(),
            alt_aas: "*".into(),
            ref_window: b"TAA".to_vec(),
            alt_window: b"TAG".to_vec(),
            ref_codons: "taa".into(),
            alt_codons: "taG".into(),
            ref_len: 1,
            alt_len: 1,
            tl_start: 8,
            tl_end: 8,
            partial_codon: false,
        };
        let got = window_terms(&same_place, EditTests::default(), peptide, (22, 22), 1);
        assert_eq!(got, vec![Consequence::StopRetainedVariant]);
    }

    /// A window that resolves one residue and not the next earns both terms.
    ///
    /// Ensembl's `missense_variant` (l. 1085) asks only that the two peptides differ and that the
    /// codon is complete; the `X` silences `synonymous_variant` (l. 1077) and not it, and
    /// `coding_unknown` (l. 1616) defers to neither. Reporting `coding_sequence_variant` alone made a
    /// real missense MODIFIER, which anything filtering on IMPACT drops: IL7R `5:35860925 GC>TT`
    /// on ENST00000515665, `QX/HX`, is `missense_variant,coding_sequence_variant` to real VEP 115.1.
    #[test]
    fn an_unknown_residue_beside_a_changed_one_keeps_the_missense() {
        let window = |ref_aas: &str, alt_aas: &str, partial_codon: bool| CodonWindow {
            ref_aas: ref_aas.into(),
            alt_aas: alt_aas.into(),
            ref_window: b"CAGC".to_vec(),
            alt_window: b"CATT".to_vec(),
            ref_codons: "caGC".into(),
            alt_codons: "caTT".into(),
            ref_len: 2,
            alt_len: 2,
            tl_start: 52,
            tl_end: 53,
            partial_codon,
        };
        let peptide = "M".repeat(100);
        let terms =
            |w: &CodonWindow| window_terms(w, EditTests::default(), &peptide, (154, 155), 2);

        assert_eq!(
            terms(&window("QX", "HX", false)),
            vec![
                Consequence::MissenseVariant,
                Consequence::CodingSequenceVariant
            ]
        );

        // Unchanged where nothing resolves: an `X` on both sides is not a
        // synonymous change, because the reference residue was never known.
        assert_eq!(
            terms(&window("X", "X", false)),
            vec![Consequence::CodingSequenceVariant]
        );

        // And `missense_variant` still declines on a partial codon, as Ensembl's does (l. 1093),
        // leaving the terms that describe one.
        assert_eq!(
            terms(&window("QX", "HX", true)),
            vec![
                Consequence::IncompleteTerminalCodonVariant,
                Consequence::CodingSequenceVariant
            ]
        );
    }

    /// The `include` of `missense_variant` (Config.pm) keeps it off any change that moves the
    /// length, and `frameshift_variant`'s (`snp => 0`) off any that does not. Without the first, a
    /// frameshift whose alternate peptide is as long as the reference (`K/X`) is a missense
    /// variant: the predicate itself only compares two peptides.
    #[test]
    fn the_length_gates_of_missense_and_frameshift_hold() {
        // `AAG` (K) with a base removed: the window shrinks to `AG`, which translates as `X`
        let frameshift = CodonWindow {
            ref_aas: "K".into(),
            alt_aas: "X".into(),
            ref_window: b"AAG".to_vec(),
            alt_window: b"AG".to_vec(),
            ref_codons: "aAg".into(),
            alt_codons: "ag".into(),
            ref_len: 1,
            alt_len: 0,
            tl_start: 5,
            tl_end: 5,
            partial_codon: false,
        };
        let peptide = "M".repeat(100);
        let got = window_terms(&frameshift, EditTests::default(), &peptide, (14, 14), 1);
        assert_eq!(got, vec![Consequence::FrameshiftVariant], "got {got:?}");
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
            location: TranscriptLocation::default(),
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
            initiator_met: false,
            cds_flags_from_sequence: false,
            readthrough: false,
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
            initiator_met: false,
            cds_flags_from_sequence: false,
            readthrough: false,
        }
    }

    /// Nine predicates run over the codon window and every one that holds is a
    /// term Ensembl would report. Keeping two dropped the third silently.
    ///
    /// A change that replaces the initiator, introduces a stop and shifts the
    /// frame holds `stop_gained`, `frameshift_variant` and `start_lost` at once.
    /// `start_lost` ranks below the other two, so it was the one lost.
    ///
    /// PVS1 happens to be unaffected here - `NonsenseOrFrameshift` outranks
    /// `StartLost` when both are present, so the criterion takes the same branch
    /// either way - but the reported consequence set was still short a term
    /// Ensembl's model holds.
    #[test]
    fn every_coding_term_that_holds_is_reported() {
        let predictor = ConsequencePredictor::default();
        let tr = make_coding_transcript();
        // CDS 1-3 is the initiator.
        let pos = GenomicPosition::new("chr1", 1050, 1052, Strand::Forward);
        let result = predictor.predict(
            &pos,
            &Allele::from_str("ATG"),
            &[Allele::from_str("TAAA")],
            &[&tr],
            None,
        );
        let got = &result.transcript_consequences[0].allele_consequences[0].consequences;
        for term in [
            Consequence::StopGained,
            Consequence::FrameshiftVariant,
            Consequence::StartLost,
        ] {
            assert!(got.contains(&term), "{term:?} missing from {got:?}");
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

    /// What Ensembl 116.3 writes for an allele that is not plain bases, checked in
    /// `fixtures/vep_tab_alleles`: no peptide, so no amino acids and no term about the protein.
    /// `N` in an SNV was `missense_variant,coding_sequence_variant` with `K/X`.
    #[test]
    fn an_allele_with_an_unknown_base_has_no_protein_change() {
        let predictor = ConsequencePredictor::default();
        let tr = make_coding_transcript();
        // CDS 1050.. is ATG GCT TCA AAG CCC: 1053 is the first base of GCT
        let predict = |start: u64, end: u64, reference: &str, alt: &str| {
            let pos = GenomicPosition::new("chr1", start, end, Strand::Forward);
            let result = predictor.predict(
                &pos,
                &Allele::from_str(reference),
                &[Allele::from_str(alt)],
                &[&tr],
                None,
            );
            result.transcript_consequences[0].allele_consequences[0].clone()
        };

        let snv = predict(1053, 1053, "G", "N");
        assert_eq!(snv.consequences, vec![Consequence::CodingSequenceVariant]);
        assert_eq!(snv.impact, Impact::Modifier);
        assert_eq!(snv.amino_acids, None);
        // the codons are still shown, as VEP shows `Aaa/Naa`
        assert_eq!(snv.codons, Some(("Gct".to_string(), "Nct".to_string())));

        // an unknown REF is the same: the reference is not a base to compare against
        let unknown_ref = predict(1053, 1053, "N", "A");
        assert_eq!(
            unknown_ref.consequences,
            vec![Consequence::CodingSequenceVariant]
        );
        assert_eq!(unknown_ref.amino_acids, None);

        // lengths are still known, so the frame is: one base is a frameshift, three are not
        let one = predict(1054, 1053, "-", "N");
        assert_eq!(one.consequences, vec![Consequence::FrameshiftVariant]);
        assert_eq!(one.amino_acids, None);
        let three = predict(1056, 1055, "-", "NNN");
        assert_eq!(three.consequences, vec![Consequence::CodingSequenceVariant]);
        let replaced = predict(1053, 1055, "GCT", "NN");
        assert_eq!(replaced.consequences, vec![Consequence::FrameshiftVariant]);

        // a symbolic allele has no length either, so it cannot move the frame
        let symbolic = predict(1053, 1053, "G", "<*>");
        assert_eq!(
            symbolic.consequences,
            vec![Consequence::CodingSequenceVariant]
        );
        assert_eq!(symbolic.amino_acids, None);

        // and a plain allele is untouched
        let plain = predict(1053, 1053, "G", "C");
        assert!(plain.amino_acids.is_some());
        assert_ne!(plain.consequences, vec![Consequence::CodingSequenceVariant]);
    }

    #[test]
    fn the_complement_of_an_ambiguity_code_is_ensembls() {
        // `tr/acgtrymkswhbvdnACGTRYMKSWHBVDN/tgcayrkmswdvbhnTGCAYRKMSWDVBHN/`
        for (a, b) in [
            (b'R', b'Y'),
            (b'Y', b'R'),
            (b'M', b'K'),
            (b'K', b'M'),
            (b'S', b'S'),
            (b'W', b'W'),
            (b'H', b'D'),
            (b'D', b'H'),
            (b'B', b'V'),
            (b'V', b'B'),
            (b'N', b'N'),
            (b'a', b't'),
            (b'r', b'y'),
        ] {
            assert_eq!(complement(a), b, "{}", a as char);
        }
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
            initiator_met: false,
            cds_flags_from_sequence: false,
            readthrough: false,
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

    /// A first codon that is not `ATG` and that the annotation reads as Met (`initiator_met`) has
    /// `M` for its reference residue: `Amino_acids` is `M/L` for a change that keeps the codon Leu,
    /// and the change is a `start_lost` and not also a `synonymous_variant`. That is what Ensembl
    /// VEP 116.3's cache writes for a complete CTG start (`M/L`, `Ctg/Ttg`, `start_lost`) and for
    /// RefSeq's GTG one (`M/V`). Without the edit (RefSeq's CTG starts have none) the residue is
    /// the code table's, `L`, and the cache writes `L`, `ctG/ctT`, `start_lost,synonymous_variant`:
    /// a change to a start codon is `start_lost` whenever the codon it leaves is not `ATG`, and
    /// the two peptides are equal. The alternate codon is translated plainly either way.
    #[test]
    fn an_initiator_met_makes_the_reference_residue_of_codon_one_m() {
        let predictor = ConsequencePredictor::default();
        for strand in [Strand::Forward, Strand::Reverse] {
            let amino_acids_and_terms = |initiator_met: bool| {
                let mut tr = make_boundary_transcript(strand);
                let cds = format!("CTG{}TAA", "GCT".repeat(8));
                tr.translateable_seq = Some(cds.clone());
                tr.peptide = Some(format!(
                    "{}{}*",
                    if initiator_met { "M" } else { "L" },
                    "A".repeat(8)
                ));
                tr.initiator_met = initiator_met;
                // the third base of codon 1 (cDNA 13): G to A keeps it Leu
                let (lo, hi) = genomic_span(strand, 13, 13);
                let (r, a) = match strand {
                    Strand::Forward => ("G", "A"),
                    Strand::Reverse => ("C", "T"),
                };
                let pos = GenomicPosition::new("chr1", lo, hi, Strand::Forward);
                let alt = Allele::from_str(a);
                let result = predictor.predict(
                    &pos,
                    &Allele::from_str(r),
                    std::slice::from_ref(&alt),
                    &[&tr],
                    None,
                );
                let ac = &result.transcript_consequences[0].allele_consequences[0];
                (ac.amino_acids.clone(), ac.consequences.clone())
            };
            let (aa, terms) = amino_acids_and_terms(true);
            assert_eq!(aa, Some(("M".to_string(), "L".to_string())), "{strand:?}");
            assert_eq!(terms, [Consequence::StartLost], "{strand:?}");
            let (aa, terms) = amino_acids_and_terms(false);
            assert_eq!(aa, Some(("L".to_string(), "L".to_string())), "{strand:?}");
            assert_eq!(
                terms,
                [Consequence::StartLost, Consequence::SynonymousVariant],
                "{strand:?}"
            );
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

    /// [`make_boundary_transcript`] with a CDS of its own and `utr5` bases of 5' UTR in front of
    /// it (the 3' UTR is the filler that keeps the transcript 201 bases long), so that a change
    /// at the start codon can be decided against the edited sequence.
    fn make_start_codon_transcript(strand: Strand, utr5: usize, cds: &str) -> Transcript {
        assert_eq!(cds.len(), 30);
        let mut tr = make_boundary_transcript(strand);
        let spliced = format!(
            "{}{}{}",
            "C".repeat(utr5),
            cds,
            "G".repeat(201 - utr5 - cds.len())
        );
        let (crs, cre) = match strand {
            Strand::Forward => (1000 + utr5 as u64, 999 + utr5 as u64 + 30),
            Strand::Reverse => (1171 - utr5 as u64, 1200 - utr5 as u64),
        };
        tr.spliced_seq = Some(spliced);
        tr.translateable_seq = Some(cds.to_string());
        tr.cdna_coding_start = Some(utr5 as u64 + 1);
        tr.cdna_coding_end = Some(utr5 as u64 + 30);
        tr.coding_region_start = Some(crs);
        tr.coding_region_end = Some(cre);
        if let Some(t) = tr.translation.as_mut() {
            t.genomic_start = crs;
            t.genomic_end = cre;
            t.start_exon_offset = utr5 as u64;
            t.end_exon_offset = utr5 as u64 + 30;
        }
        tr
    }

    fn terms_at_the_start_codon(
        strand: Strand,
        utr5: usize,
        cds: &str,
        cdna_lo: u64,
        cdna_hi: u64,
        ref_allele: &str,
        alt_allele: &Allele,
    ) -> Vec<Consequence> {
        let tr = make_start_codon_transcript(strand, utr5, cds);
        let (g_start, g_end) = genomic_span(strand, cdna_lo, cdna_hi);
        let pos = GenomicPosition::new("chr1", g_start, g_end, Strand::Forward);
        let result = ConsequencePredictor::default().predict(
            &pos,
            &Allele::from_str(ref_allele),
            std::slice::from_ref(alt_allele),
            &[&tr],
            None,
        );
        result.transcript_consequences[0].allele_consequences[0]
            .consequences
            .clone()
    }

    /// A protein that still begins with its Met has not lost its initiator. `GGC` at `c.3_5` of
    /// `ATG GCG GCT ...` is the deletion of codon 2 written at its leftmost position, which touches
    /// the start codon and leaves the edited transcript starting `ATGGCT`: Ensembl's peptide test,
    /// which compares the window `MA` with `M`, calls it `start_lost` beside `start_retained_variant`.
    #[test]
    fn the_deletion_of_codon_two_is_not_a_start_loss() {
        let cds = format!("ATGGCG{}TAA", "GCT".repeat(7));
        for strand in [Strand::Forward, Strand::Reverse] {
            let got = terms_at_the_start_codon(strand, 10, &cds, 13, 15, "GGC", &Allele::Deletion);
            assert!(
                !got.contains(&Consequence::StartLost),
                "{strand:?}: the Met is still there, got {got:?}"
            );
            assert!(
                got.contains(&Consequence::InframeDeletion)
                    && got.contains(&Consequence::StartRetainedVariant),
                "{strand:?}: {got:?}"
            );
            assert_eq!(Consequence::worst_impact(&got), Some(Impact::Moderate));
        }
    }

    /// What is still a start loss: the initiator deleted, a base of it changed, and a change that
    /// leaves a first residue that is not Met.
    #[test]
    fn a_change_that_takes_the_met_away_is_still_a_start_loss() {
        let cds = format!("ATGGCG{}TAA", "GCT".repeat(7));
        for strand in [Strand::Forward, Strand::Reverse] {
            // the initiator codon, in frame: codon 2 is the first one now
            let got = terms_at_the_start_codon(strand, 10, &cds, 11, 13, "ATG", &Allele::Deletion);
            assert!(got.contains(&Consequence::StartLost), "{strand:?}: {got:?}");
            // the two codons: the window `MA` becomes `A`
            let got =
                terms_at_the_start_codon(strand, 10, &cds, 11, 16, "ATGGCG", &Allele::Deletion);
            assert!(got.contains(&Consequence::StartLost), "{strand:?}: {got:?}");
            // two of its bases, a frameshift whose first codon is not ATG
            let got = terms_at_the_start_codon(strand, 10, &cds, 12, 13, "TG", &Allele::Deletion);
            assert!(got.contains(&Consequence::StartLost), "{strand:?}: {got:?}");
            // a substitution
            let (r, a) = match strand {
                Strand::Forward => ("G", "C"),
                Strand::Reverse => ("C", "G"),
            };
            let got = terms_at_the_start_codon(strand, 10, &cds, 13, 13, r, &Allele::from_str(a));
            assert!(got.contains(&Consequence::StartLost), "{strand:?}: {got:?}");
        }
    }

    /// An insertion inside the start codon that keeps its ATG (`atg` to `atGg`) is the same edit
    /// whether the transcript has a 5' UTR or not. Ensembl asks "is the ATG still at the coding
    /// start" only for a transcript that has a UTR, so on one that begins at the initiator the
    /// frameshift is also `start_lost` (DYRK1A `c.3dup`: `start_retained_variant` on the
    /// transcripts with a UTR, `start_lost` on ENST00000398956).
    #[test]
    fn a_transcript_with_no_utr_reads_an_insertion_that_keeps_the_atg_as_one_with_a_utr_does() {
        let cds = format!("ATGGCT{}TAA", "GCT".repeat(7));
        for strand in [Strand::Forward, Strand::Reverse] {
            // a G between cDNA 2 and 3 of the start codon, written as T>TG in transcript terms
            let at = |utr5: usize| {
                let lo = utr5 as u64 + 2;
                let (r, a) = match strand {
                    Strand::Forward => ("T", "TG"),
                    Strand::Reverse => ("A", "CA"),
                };
                terms_at_the_start_codon(strand, utr5, &cds, lo, lo, r, &Allele::from_str(a))
            };
            let with_utr = at(10);
            let without_utr = at(0);
            assert!(
                with_utr.contains(&Consequence::StartRetainedVariant)
                    && with_utr.contains(&Consequence::FrameshiftVariant),
                "{strand:?}: {with_utr:?}"
            );
            assert!(!with_utr.contains(&Consequence::StartLost), "{strand:?}");
            assert_eq!(
                without_utr, with_utr,
                "{strand:?}: the UTR changed the answer"
            );
        }
    }

    /// A deletion that covers the whole transcript is `transcript_ablation` and nothing else, on
    /// either strand, for a coding transcript and a non-coding one. A span one base short of
    /// either end of it is not.
    #[test]
    fn a_deletion_that_covers_the_transcript_is_a_transcript_ablation() {
        let ablated = |tr: &Transcript, lo: u64, hi: u64, reference: &str, alt: &Allele| {
            let pos = GenomicPosition::new("chr1", lo, hi, Strand::Forward);
            let result = ConsequencePredictor::default().predict(
                &pos,
                &Allele::from_str(reference),
                std::slice::from_ref(alt),
                &[tr],
                None,
            );
            let ac = &result.transcript_consequences[0].allele_consequences[0];
            (ac.consequences.clone(), ac.impact)
        };
        for strand in [Strand::Forward, Strand::Reverse] {
            let coding = make_boundary_transcript(strand);
            let non_coding = {
                let mut t = make_boundary_transcript(strand);
                t.biotype = "lncRNA".into();
                t.translation = None;
                t.cdna_coding_start = None;
                t.cdna_coding_end = None;
                t.coding_region_start = None;
                t.coding_region_end = None;
                t.translateable_seq = None;
                t
            };
            let long = "A".repeat(230);
            for tr in [&coding, &non_coding] {
                for (lo, hi) in [(990, 1210), (1000, 1200), (1000, 1210), (990, 1200)] {
                    let reference = &long[..(hi - lo + 1) as usize];
                    let (got, impact) = ablated(tr, lo, hi, reference, &Allele::Deletion);
                    assert_eq!(
                        got,
                        [Consequence::TranscriptAblation],
                        "{strand:?} {lo}..{hi}"
                    );
                    assert_eq!(impact, Impact::High);
                }
                // a shorter replacement is a deletion too, as Ensembl's `deletion` has it
                let (got, _) = ablated(tr, 990, 1210, &long[..221], &Allele::from_str("AA"));
                assert_eq!(got, [Consequence::TranscriptAblation], "{strand:?} delins");
                // one base short of an end: an exonic deletion like any other
                for (lo, hi) in [(1001, 1210), (990, 1199)] {
                    let reference = &long[..(hi - lo + 1) as usize];
                    let (got, _) = ablated(tr, lo, hi, reference, &Allele::Deletion);
                    assert!(
                        !got.contains(&Consequence::TranscriptAblation),
                        "{strand:?} {lo}..{hi}: {got:?}"
                    );
                }
                // a change that is not shorter than what it replaces is not a deletion: a
                // transcript with a protein is a `coding_transcript_variant`, and one without is
                // `non_coding_transcript_variant`, not an exon variant
                let (got, _) =
                    ablated(tr, 990, 1210, &long[..221], &Allele::from_str(&long[..221]));
                assert!(!got.contains(&Consequence::TranscriptAblation), "{got:?}");
                assert!(
                    !got.contains(&Consequence::CodingSequenceVariant),
                    "{got:?}"
                );
                assert!(
                    !got.contains(&Consequence::NonCodingTranscriptExonVariant),
                    "{got:?}"
                );
                if tr.is_coding() {
                    assert!(
                        got.contains(&Consequence::CodingTranscriptVariant),
                        "{got:?}"
                    );
                } else {
                    assert!(
                        got.contains(&Consequence::NonCodingTranscriptVariant),
                        "{got:?}"
                    );
                }
            }
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
    ///
    /// "Exactly on" is between cDNA 10 and 11 for the initiator and between 40 and 41 for the
    /// terminator, which is a different genomic base on the two strands: the insertion is the
    /// zero-length interval that ends one before the base that follows it *in the genome*. One
    /// base further in is inside the codon, and that is another matter.
    #[test]
    fn an_insertion_on_a_coding_boundary_reports_only_the_utr_term() {
        let predictor = ConsequencePredictor::default();
        for strand in [Strand::Forward, Strand::Reverse] {
            let tr = make_boundary_transcript(strand);
            // The initiator is cDNA 11 and the terminator ends at cDNA 40; `cdna_flank` is the
            // base the insertion sits in front of, in transcript order.
            for (cdna_flank, forbidden) in
                [(11u64, Consequence::StartLost), (41, Consequence::StopLost)]
            {
                let genomic = match strand {
                    Strand::Forward => 999 + cdna_flank,
                    // the base after the insertion point in genomic order is the one before it in
                    // transcript order
                    Strand::Reverse => 1202 - cdna_flank,
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
    /// The same insertion one base further in, between the first and second base of the
    /// initiator, destroys it: `ATG` becomes `ATTTTG`, and the first codon is no longer an ATG.
    /// Ensembl's `_inv_start_altered` reads the codon that now sits where the initiator did, so
    /// this is `start_lost` on either strand, and it must not depend on which genomic side the
    /// transcript's 5' end is on.
    #[test]
    fn an_insertion_inside_the_initiator_is_a_start_loss_on_either_strand() {
        let predictor = ConsequencePredictor::default();
        for strand in [Strand::Forward, Strand::Reverse] {
            let tr = make_boundary_transcript(strand);
            // between cDNA 11 and 12
            let genomic = match strand {
                Strand::Forward => 999 + 12,
                Strand::Reverse => 1202 - 12,
            };
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
                got.contains(&Consequence::StartLost),
                "{strand:?} insertion inside the initiator: expected start_lost, got {got:?}"
            );
        }
    }

    /// The boundary transcript cut off at its terminator: the coding region reaches the 3' end of
    /// the transcript, so there is no 3' UTR to speak of.
    fn make_transcript_ending_at_its_stop(strand: Strand) -> Transcript {
        let mut tr = make_boundary_transcript(strand);
        // cDNA 1..=40 are the 5' UTR (10 bases) and the CDS (30)
        let (start, end) = match strand {
            Strand::Forward => (1000, 1039),
            Strand::Reverse => (1161, 1200),
        };
        tr.start = start;
        tr.end = end;
        tr.gene.start = start;
        tr.gene.end = end;
        tr.exons[0].start = start;
        tr.exons[0].end = end;
        tr.spliced_seq = tr.spliced_seq.map(|s| s[..40].to_string());
        tr
    }

    /// `_before_coding` and `_after_coding` are plain `overlap` calls against the UTR interval,
    /// and `overlap` does not ask that the interval hold a base. For a transcript whose coding
    /// region reaches its own end the interval is empty (`[end + 1, end]`), and a change that
    /// spans that end - the stop codon and on past the last base - still overlaps it, so it is a
    /// `3_prime_UTR_variant`. VEP 116.3 writes `stop_retained_variant,3_prime_UTR_variant` for
    /// ClinVar's 265214 (KRIT1-215, ENST00000452773, `7:92241100-92241103 ACTT>-`, whose stop
    /// codon is the last three bases of the transcript) and for eight more like it; a rule that
    /// asked first whether there is a UTR on that side left the term off all nine.
    #[test]
    fn a_change_spanning_the_end_of_a_transcript_with_no_utr_is_still_a_utr_variant() {
        let predictor = ConsequencePredictor::default();
        for strand in [Strand::Forward, Strand::Reverse] {
            let tr = make_transcript_ending_at_its_stop(strand);
            let run = |cdna_lo: u64, cdna_hi: u64| {
                let (lo, hi) = genomic_span(strand, cdna_lo, cdna_hi);
                let pos = GenomicPosition::new("chr1", lo, hi, Strand::Forward);
                let reference = Allele::Sequence(vec![b'A'; (hi - lo + 1) as usize]);
                predictor
                    .predict(&pos, &reference, &[Allele::Deletion], &[&tr], None)
                    .transcript_consequences[0]
                    .allele_consequences[0]
                    .consequences
                    .clone()
            };
            // the last two bases of the stop codon (cDNA 39, 40) and two bases beyond the transcript
            let spanning = run(39, 42);
            assert!(
                spanning.contains(&Consequence::ThreePrimeUtrVariant),
                "{strand:?}: a deletion over the end of the transcript, got {spanning:?}"
            );
            // the same deletion without the two bases beyond it has no UTR term
            let inside = run(39, 40);
            assert!(
                !inside.contains(&Consequence::ThreePrimeUtrVariant),
                "{strand:?}: a deletion inside the transcript, got {inside:?}"
            );
        }
    }

    /// Ensembl reads the alternate codon out of the edited CDS with the 3' UTR appended
    /// (`_get_alternate_cds`) and the reference codon out of the CDS alone. A last codon that the
    /// CDS leaves short is therefore two bases on the reference side and three on the alternate
    /// side, the third from the UTR: ClinVar's 1011962, ENST00000512526 `17:50167975 A>C`, is
    /// `X/A` and `gA/gCt` in VEP 116.3, and was `X` and `gA/gC` here.
    #[test]
    fn a_last_codon_the_cds_leaves_short_is_completed_from_the_utr_on_the_alternate_side() {
        let predictor = ConsequencePredictor::default();
        let mut tr = make_boundary_transcript(Strand::Forward);
        // the CDS now ends one base early, at cDNA 39, so its last codon `TA` is short
        tr.cdna_coding_end = Some(39);
        tr.coding_region_end = Some(1038);
        tr.translateable_seq = tr.translateable_seq.map(|s| s[..29].to_string());
        let (g, _) = genomic_span(Strand::Forward, 39, 39);
        let pos = GenomicPosition::new("chr1", g, g, Strand::Forward);
        let result = predictor.predict(
            &pos,
            &Allele::from_str("A"),
            &[Allele::from_str("G")],
            &[&tr],
            None,
        );
        let got = &result.transcript_consequences[0].allele_consequences[0];
        // `TGA` is a stop
        assert_eq!(got.amino_acids, Some(("X".into(), "*".into())));
        assert_eq!(got.codons, Some(("tA".into(), "tGa".into())));
        assert!(got
            .consequences
            .contains(&Consequence::IncompleteTerminalCodonVariant));
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

    /// A change that crosses an intron is measured in the CDS, as Ensembl measures it.
    ///
    /// `codon` takes the CDS span of the change (`cds_end - cds_start + 1`, the introns between
    /// the coding bases not counted) as the length of what the alternate allele replaces, and
    /// `frameshift` does the same arithmetic. So a deletion of three bases at the end of one exon,
    /// the intron and three bases at the start of the next is a six-base deletion of the coding
    /// sequence, an in-frame one, with the two codons it removed named; one base longer it is a
    /// frameshift. Real VEP 116.3 on the edge probe: `LLLPSI/X`, `p.Leu98HisfsTer3` for a deletion
    /// of 17 coding bases across an intron of 600. This used to resolve to no residues and
    /// `coding_sequence_variant`, on the grounds that the change was "not contiguous", which was
    /// true of the genome and not of the coding sequence.
    #[test]
    fn a_change_across_an_intron_counts_only_its_coding_bases() {
        let predictor = ConsequencePredictor::default();
        let tr = make_coding_transcript();
        // exon 1 ends at 1200 and exon 2 starts at 2000; CDS position of 1198 is 149, of 2000 is 152
        for (end, want, coding_bases) in [
            (2002u64, Consequence::InframeDeletion, 6),
            (2003, Consequence::FrameshiftVariant, 7),
        ] {
            let pos = GenomicPosition::new("chr1", 1198, end, Strand::Forward);
            let reference = Allele::Sequence(vec![b'A'; (end - 1198 + 1) as usize]);
            let result = predictor.predict(&pos, &reference, &[Allele::Deletion], &[&tr], None);
            let got = &result.transcript_consequences[0].allele_consequences[0];
            assert!(
                got.consequences.contains(&want),
                "{coding_bases} coding bases: expected {want:?}, got {:?}",
                got.consequences
            );
            assert!(
                got.amino_acids.is_some(),
                "{coding_bases} coding bases: the residues it changes are named"
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

    /// The coding fixture with its first intron cut to 13 bases (`1201-1213`), which makes it a
    /// frameshift intron, and the exon after it moved to follow.
    fn make_transcript_with_a_frameshift_intron() -> Transcript {
        let mut tr = make_coding_transcript();
        tr.exons[1].start = 1214;
        tr.exons[1].end = 1514;
        tr
    }

    #[test]
    fn a_change_in_a_frameshift_intron_inside_the_cds_is_a_coding_sequence_variant() {
        // `within_cds` holds for a change in a frameshift intron that lies within the coding
        // region, so `coding_unknown` does: no splice term, no `intron_variant`, and not the
        // `coding_transcript_variant` the catch-all would have given it.
        let tr = make_transcript_with_a_frameshift_intron();
        let pos = GenomicPosition::new("chr1", 1205, 1205, Strand::Forward);
        let result = ConsequencePredictor::default().predict(
            &pos,
            &Allele::from_str("A"),
            &[Allele::from_str("G")],
            &[&tr],
            None,
        );
        let got = &result.transcript_consequences[0].allele_consequences[0].consequences;
        assert_eq!(got, &vec![Consequence::CodingSequenceVariant]);
    }

    #[test]
    fn the_same_change_in_an_ordinary_intron_is_an_intron_variant() {
        let tr = make_coding_transcript();
        let pos = GenomicPosition::new("chr1", 1500, 1500, Strand::Forward);
        let result = ConsequencePredictor::default().predict(
            &pos,
            &Allele::from_str("A"),
            &[Allele::from_str("G")],
            &[&tr],
            None,
        );
        let got = &result.transcript_consequences[0].allele_consequences[0].consequences;
        assert_eq!(got, &vec![Consequence::IntronVariant]);
    }

    /// `within_cdna`, which gates the UTR terms, holds for a change in a frameshift intron as
    /// `within_cds` does for the CDS: one in a short intron of the UTR is a UTR variant, on the side of
    /// the coding region it is on, and not the `coding_transcript_variant` the catch-all gave it. An
    /// intron that is not short is `intron_variant`, as before.
    #[test]
    fn a_frameshift_intron_in_a_utr_is_a_utr_variant() {
        for strand in [Strand::Forward, Strand::Reverse] {
            // a 9-base intron inside the 5' UTR (on the low side) and one inside the 3' UTR (high)
            let mut tr = make_coding_transcript();
            tr.strand = strand;
            tr.exons = vec![
                (1000, 1020),
                (1030, 1200),
                (2000, 2300),
                (4000, 4600),
                (4610, 5000),
            ]
            .into_iter()
            .enumerate()
            .map(|(i, (start, end))| Exon {
                stable_id: format!("E{i}"),
                start,
                end,
                strand,
                phase: 0,
                end_phase: 0,
                rank: i as u32 + 1,
            })
            .collect();
            let terms = |tr: &Transcript, at: u64| {
                let result = ConsequencePredictor::default().predict(
                    &GenomicPosition::new("chr1", at, at, Strand::Forward),
                    &Allele::from_str("A"),
                    &[Allele::from_str("G")],
                    &[tr],
                    None,
                );
                result.transcript_consequences[0].allele_consequences[0]
                    .consequences
                    .clone()
            };
            let (low, high) = match strand {
                Strand::Forward => (
                    Consequence::FivePrimeUtrVariant,
                    Consequence::ThreePrimeUtrVariant,
                ),
                Strand::Reverse => (
                    Consequence::ThreePrimeUtrVariant,
                    Consequence::FivePrimeUtrVariant,
                ),
            };
            assert_eq!(terms(&tr, 1025), [low], "{strand:?} in the low UTR");
            assert_eq!(terms(&tr, 4605), [high], "{strand:?} in the high UTR");

            // an intron of the same place that is long is an intron, not a UTR term
            tr.exons[1].start = 1100;
            tr.exons[0].end = 1020;
            assert_eq!(
                terms(&tr, 1050),
                [Consequence::IntronVariant],
                "{strand:?}: an ordinary intron"
            );
        }
    }

    #[test]
    fn a_frameshift_intron_before_the_coding_region_is_not_a_coding_sequence_variant() {
        // `within_cds` asks that the change overlap the coding region too.
        let mut tr = make_transcript_with_a_frameshift_intron();
        tr.coding_region_start = Some(1300);
        let pos = GenomicPosition::new("chr1", 1205, 1205, Strand::Forward);
        let result = ConsequencePredictor::default().predict(
            &pos,
            &Allele::from_str("A"),
            &[Allele::from_str("G")],
            &[&tr],
            None,
        );
        let got = &result.transcript_consequences[0].allele_consequences[0].consequences;
        assert!(
            !got.contains(&Consequence::CodingSequenceVariant),
            "{got:?}"
        );
    }

    #[test]
    fn a_deletion_over_a_one_base_intron_writes_the_reference_codons_in_uppercase_by_genome_length()
    {
        // `display_codon` uppercases `length($self->feature_seq)` of the reference allele from
        // the codon position on, and that is the five bases of the genome the deletion removes,
        // though only four of them are in the CDS. The codons are cut with the CDS span (four,
        // CDS 49 to 52, a window of two codons); the case is not.
        let mut tr = make_coding_transcript();
        // exon 1 ends at 1099 (CDS 50 is its last base) and exon 2 starts at 1101: the base
        // between them is a one-base intron
        tr.exons[0].end = 1099;
        tr.exons[1].start = 1101;
        tr.exons[1].end = 1400;
        let result = ConsequencePredictor::default().predict(
            &GenomicPosition::new("chr1", 1098, 1102, Strand::Forward),
            &Allele::from_str("ACGTA"),
            &[Allele::Deletion],
            &[&tr],
            None,
        );
        let ac = &result.transcript_consequences[0].allele_consequences[0];
        let (ref_codons, _) = ac.codons.clone().expect("a codon window");
        assert_eq!(ref_codons.len(), 6, "{ref_codons}");
        let upper: Vec<bool> = ref_codons.chars().map(|c| c.is_ascii_uppercase()).collect();
        assert_eq!(upper, [true, true, true, true, true, false], "{ref_codons}");
    }
}
