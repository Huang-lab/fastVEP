use crate::codon::{reverse_complement, CodonTable};
use fastvep_core::Strand;
use serde::{Deserialize, Serialize};
use std::fmt::Write;
use std::sync::Arc;

/// The exons of a transcript in transcript order ([`Transcript::sorted_exons`]).
pub(crate) enum ExonsInOrder<'a> {
    /// `exons` as they are stored
    Stored(&'a [Exon]),
    /// `exons` from the last to the first
    Backwards(&'a [Exon]),
    /// neither, so sorted
    Sorted(Vec<&'a Exon>),
}

impl<'a> ExonsInOrder<'a> {
    pub(crate) fn len(&self) -> usize {
        match self {
            Self::Stored(e) | Self::Backwards(e) => e.len(),
            Self::Sorted(e) => e.len(),
        }
    }

    pub(crate) fn get(&self, i: usize) -> Option<&'a Exon> {
        match self {
            Self::Stored(e) => e.get(i),
            Self::Backwards(e) => e.len().checked_sub(i + 1).and_then(|j| e.get(j)),
            Self::Sorted(e) => e.get(i).copied(),
        }
    }

    pub(crate) fn iter(&self) -> impl Iterator<Item = &'a Exon> + '_ {
        (0..self.len()).filter_map(move |i| self.get(i))
    }

    pub(crate) fn last(&self) -> Option<&'a Exon> {
        self.len().checked_sub(1).and_then(|i| self.get(i))
    }
}

impl std::ops::Index<usize> for ExonsInOrder<'_> {
    type Output = Exon;

    fn index(&self, i: usize) -> &Exon {
        self.get(i).expect("an exon of the transcript")
    }
}

/// A gene model.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Gene {
    pub stable_id: Arc<str>,
    pub symbol: Option<Arc<str>>,
    pub symbol_source: Option<String>,
    pub hgnc_id: Option<String>,
    pub biotype: Arc<str>,
    pub chromosome: Arc<str>,
    pub start: u64,
    pub end: u64,
    pub strand: Strand,
}

/// A transcript model with all data needed for consequence prediction.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Transcript {
    pub stable_id: Arc<str>,
    /// Version number from annotation (e.g., 7 for ENST00000348295.7)
    pub version: Option<u32>,
    pub gene: Gene,
    pub biotype: Arc<str>,
    pub chromosome: Arc<str>,
    pub start: u64,
    pub end: u64,
    pub strand: Strand,
    pub exons: Vec<Exon>,
    pub translation: Option<Translation>,

    // Pre-computed fields for consequence prediction
    /// Start of the coding region in cDNA coordinates (1-based).
    pub cdna_coding_start: Option<u64>,
    /// End of the coding region in cDNA coordinates (1-based).
    pub cdna_coding_end: Option<u64>,
    /// Start of the coding region in genomic coordinates.
    pub coding_region_start: Option<u64>,
    /// End of the coding region in genomic coordinates.
    pub coding_region_end: Option<u64>,

    // Spliced and translated sequences
    pub spliced_seq: Option<String>,
    pub translateable_seq: Option<String>,
    pub peptide: Option<String>,

    // Annotation metadata
    pub canonical: bool,
    pub mane_select: Option<String>,
    pub mane_plus_clinical: Option<String>,
    pub tsl: Option<u8>,
    pub appris: Option<String>,
    pub ccds: Option<String>,
    pub protein_id: Option<String>,
    pub protein_version: Option<u32>,
    pub swissprot: Vec<String>,
    pub trembl: Vec<String>,
    pub uniparc: Vec<String>,
    pub refseq_id: Option<String>,
    pub source: Option<String>,
    pub gencode_primary: bool,
    /// Flags like "cds_end_NF", "cds_start_NF" from annotation tags.
    pub flags: Vec<String>,
    /// Ensembl phase offset for CDS numbering (converted from GFF3 first CDS phase).
    /// Added to the raw CDS position to account for incomplete CDS starts.
    #[serde(default)]
    pub codon_table_start_phase: u64,
    /// The annotation reads the first codon as Met although it is not `ATG`, so residue 1 is `M`
    /// where the codon table says `L` or `V`.
    ///
    /// Ensembl's database records it as an `initial_met` edit of the translation and RefSeq's
    /// GFF3 as `transl_except=(pos:...,aa:Met)`; Ensembl VEP applies either to the reference
    /// residue and to nothing else, so `Amino_acids` reads `M/L` and `HGVSp` `p.Met1?` for a CTG
    /// start. Set from the `transl_except` by the GFF3 parser, and for an Ensembl GFF3, which has
    /// no such record, by [`Transcript::build_sequences`] (see `cds_flags_from_sequence`).
    #[serde(default)]
    pub initiator_met: bool,
    /// The GFF3 this transcript came from is Ensembl's export, which leaves out two things
    /// Ensembl's database (and so VEP's cache) holds: the `cds_start_NF` / `cds_end_NF` flags and
    /// the `initial_met` edit. [`Transcript::build_sequences`] reads both off the sequence, once
    /// the CDS is in hand, and sets `flags` and `initiator_met`.
    ///
    /// False for a GFF3 that carries its own flags (GENCODE's), which are taken as they are, and
    /// for RefSeq's, whose cache has none.
    #[serde(default)]
    pub cds_flags_from_sequence: bool,
    /// The annotation tags it a readthrough transcript, one that spans two genes. Ensembl VEP's
    /// cache leaves these out (`AnnotationSource::Database::Transcript`), a GFF3 does not.
    #[serde(default)]
    pub readthrough: bool,
}

/// Reinterpret an in-frame UGA that the annotation places *inside* the coding
/// sequence as selenocysteine.
///
/// A GFF3 CDS is an assertion about where translation starts and ends. When it
/// runs past an in-frame stop, the annotation is saying that codon is read
/// through, and the only codon the genetic code reads through this way is UGA,
/// which encodes selenocysteine. Translating it as a terminator instead makes
/// the protein look like it ends partway: SELENOW is 87 residues and its UGA is
/// at residue 13, so consumers that stop at the first `*` see a 13-residue
/// protein. That is what blocked HGVSp normalisation past the Sec position, and
/// it would make a substitution at the Sec codon look like `stop_lost` rather
/// than the missense it is.
///
/// Three conditions, all read from the annotation rather than from a list of
/// gene names:
///
///   1. The coding sequence is a whole number of codons.
///   2. Its final codon is a terminator, so the annotated end is where
///      translation really stops and the frame is corroborated.
///   3. The internal codon is literally `TGA`, and the active codon table calls
///      `TGA` a stop. Selenocysteine has no other codon, so `TAA`/`TAG` inside a
///      CDS are left as terminators - they signal an annotation problem rather
///      than readthrough, and silently translating them would invent residues.
///
/// A CDS failing 1 or 2 is untouched, which leaves incomplete and
/// internally-inconsistent annotations exactly as they were. The vertebrate
/// mitochondrial table reads `TGA` as Trp, so condition 3 never fires there.
fn resolve_readthrough_selenocysteine(peptide: &mut [u8], cds: &[u8], table: &CodonTable) {
    if !table.is_stop(b"TGA") {
        return;
    }
    // A partial final codon means the annotated CDS is not a whole protein, so
    // its interior is not trustworthy either.
    if !cds.len().is_multiple_of(3) {
        return;
    }
    // The last residue must be the terminator the annotation promises. Without
    // it we cannot tell a readthrough codon from a CDS annotated in the wrong
    // frame, and reading through would extend a protein that does end early.
    match peptide.last() {
        Some(b'*') => {}
        _ => return,
    }

    let last = peptide.len() - 1;
    for i in 0..last {
        if peptide[i] != b'*' {
            continue;
        }
        let codon = &cds[i * 3..i * 3 + 3];
        if codon.eq_ignore_ascii_case(b"TGA") {
            peptide[i] = b'U';
        }
    }
}

/// Where a genomic span sits against a transcript's exons and introns. See
/// [`Transcript::span_site`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SpanSite {
    /// Wholly inside one exon.
    Exon,
    /// Wholly inside one intron. `offset` is the signed distance from the exon
    /// boundary of the nearest base the span covers: `+N` counted from the
    /// donor, `-N` from the acceptor.
    Intron { offset: i64 },
    /// Anything else: it reaches an exon and an intron, crosses an exon between
    /// two introns, covers a whole intron, or does not touch the transcript.
    Spans,
}

/// One intron and the exon boundaries on either side of it.
struct IntronBounds {
    start: u64,
    end: u64,
    /// Genomic position of the last base of the upstream exon.
    donor: u64,
    /// Genomic position of the first base of the downstream exon.
    acceptor: u64,
    strand: Strand,
}

impl IntronBounds {
    /// The signed offset of a position inside this intron: `+N` from the donor
    /// when that is the nearer boundary, `-N` from the acceptor otherwise.
    fn offset_of(&self, pos: u64) -> i64 {
        let (from_donor, from_acceptor) = match self.strand {
            Strand::Forward => (pos - self.donor, self.acceptor - pos),
            Strand::Reverse => (self.donor - pos, pos - self.acceptor),
        };
        if from_donor <= from_acceptor {
            from_donor as i64
        } else {
            -(from_acceptor as i64)
        }
    }
}

/// How far apart an intron's two ends may be for Ensembl to call it a frameshift intron:
/// `abs($intron_end - $intron_start) <= 12` in `BaseTranscriptVariation.pm`, which is an intron
/// of 13 bases or fewer. The comment above that test says 12, and the test is the one that runs.
pub const FRAMESHIFT_INTRON_SPAN: u64 = 12;

/// How many bases Ensembl widens every exon by, for the question "does this change touch an
/// exon", in a transcript that has a frameshift intron (`$stretch` in `_overlapped_exons`).
pub const FRAMESHIFT_EXON_STRETCH: u64 = 12;

impl Transcript {
    /// Whether this transcript is protein-coding.
    pub fn is_coding(&self) -> bool {
        self.translation.is_some()
    }

    /// Whether some intron is short enough for Ensembl to call it a frameshift intron: it is
    /// there to make the transcript's sequence agree with an assembly that has an extra or a
    /// missing base, and Ensembl treats a change in one as exonic ([`FRAMESHIFT_INTRON_SPAN`]).
    pub fn has_frameshift_intron(&self) -> bool {
        self.introns()
            .any(|intron| intron.end - intron.start <= FRAMESHIFT_INTRON_SPAN)
    }

    /// Ensembl's "overlaps an exon" (`_overlapped_exons` is not empty), the `exon` flag that
    /// `_bvfo_preds` caches and that the polypyrimidine term and the `coding` flag are gated on.
    ///
    /// It is not the same question as "is inside an exon", twice over. The pair it asks about
    /// is the change's two flanking bases sorted (`min_vf`, `max_vf`), so an insertion touches
    /// the exon if either of the bases it sits between is in it, where the exon an insertion is
    /// *in* (`non_coding_exon_variant`, which re-checks with a plain `overlap`) needs both. And
    /// in a transcript with a frameshift intron every exon is [`FRAMESHIFT_EXON_STRETCH`] bases
    /// wider on each side, so that a change in or next to the tiny intron "is exonic": which also
    /// makes every change within 12 bases of any exon of that transcript touch one.
    pub fn overlaps_exon_flag(&self, start: u64, end: u64) -> bool {
        let (lo, hi) = (start.min(end), start.max(end));
        let reaches = |stretch: u64| {
            self.exons.iter().any(|exon| {
                lo <= exon.end.saturating_add(stretch) && hi >= exon.start.saturating_sub(stretch)
            })
        };
        reaches(0) || (reaches(FRAMESHIFT_EXON_STRETCH) && self.has_frameshift_intron())
    }

    /// `within_feature` of `_bvfo_preds`: a plain `overlap` of the change with the transcript, so
    /// an insertion (`end = start - 1`) needs both of the bases it sits between inside it.
    pub fn within_feature(&self, start: u64, end: u64) -> bool {
        end >= self.start && start <= self.end
    }

    /// `$pre->{exon}` of `_bvfo_preds`: the change is within the transcript and touches an exon
    /// ([`Self::overlaps_exon_flag`]). It decides whether Ensembl writes a `cDNA_position` at all.
    pub fn pre_exon(&self, start: u64, end: u64) -> bool {
        self.within_feature(start, end) && self.overlaps_exon_flag(start, end)
    }

    /// `$pre->{coding}` of `_bvfo_preds`, which `hgvs_protein` and the `CDS_position` and
    /// `Protein_position` cells ask of the change *as it stands*, before the 3'-rule moves it.
    ///
    /// It is a test of extent and not of coordinates. The change has to be `within_feature`, a
    /// plain `overlap` with the transcript, so an insertion needs both of the bases it sits
    /// between inside it; it has to reach the coding region, and an exon, where an insertion
    /// counts as the two bases it sits between (`min_vf` and `max_vf`). Whether the mapper then
    /// finds a CDS coordinate for it is not asked: an insertion in front of the initiator has
    /// none and is `coding` all the same, and the 3'-rule can carry it onto the initiator, which
    /// is how Ensembl writes `p.Met1AsnfsTer?` for one, and `-/654` is what `--total_length`
    /// writes for the `CDS_position` of a deletion that spans the whole of an exon from intron
    /// to intron.
    ///
    /// "An exon" is [`Self::overlaps_exon_flag`], which in a transcript with a frameshift
    /// intron is every base within 12 of an exon: a deletion inside the intron is `coding`, and
    /// the 3'-rule can carry it into the exon, where Ensembl then writes the protein change of it
    /// (`p.Glu56SerfsTer136` for a deletion in SHANK3's two-base intron, `c.166del`).
    pub fn pre_coding(&self, start: u64, end: u64) -> bool {
        if !self.within_feature(start, end) {
            return false;
        }
        let (Some(coding_start), Some(coding_end)) =
            (self.coding_region_start, self.coding_region_end)
        else {
            return false;
        };
        let (lo, hi) = (start.min(end), start.max(end));
        lo <= coding_end && hi >= coding_start && self.overlaps_exon_flag(lo, hi)
    }

    /// The length of the coding sequence as Ensembl counts it for `--total_length`:
    /// `length(translateable_seq)`, which includes the stop codon and, for a transcript whose first
    /// codon is partial, the `N`s in front of it. `None` for a transcript with no coding sequence.
    ///
    /// Without a FASTA there is no sequence to measure, and the length is worked out from the
    /// coding span and the phase instead, which is the same number.
    pub fn cds_length(&self) -> Option<u64> {
        if let Some(seq) = self.translateable_seq.as_deref().filter(|s| !s.is_empty()) {
            return Some(seq.len() as u64);
        }
        if !self.is_coding() {
            return None;
        }
        let span = self.cdna_coding_end?.checked_sub(self.cdna_coding_start?)? + 1;
        Some(span + self.codon_table_start_phase)
    }

    /// The protein's accession as `HGVSp` writes it: `ENSP00000363557.3`.
    ///
    /// Ensembl's GFF3 keeps the version beside the id (`protein_id=ENSP00000363557;version=3`),
    /// so it is appended. NCBI's carries it in the id (`protein_id=NP_000001.2`), and Ensembl
    /// VEP writes that as it is (`NP_000001.2:p.Ala6Thr`): appending the version to it would
    /// make `NP_000001.2.1`, an accession that does not exist, so a version the id already
    /// has stands. (The old check compared the two versions and only caught an id that ended
    /// in `.1`, which is the one case where the answer was right anyway.)
    pub fn versioned_protein_id(&self) -> Option<String> {
        let id = self.protein_id.as_ref()?;
        let has_version = id
            .rsplit_once('.')
            .is_some_and(|(_, v)| !v.is_empty() && v.bytes().all(|b| b.is_ascii_digit()));
        Some(match self.protein_version {
            // sized for what it holds: `format!` grows from nothing, and this runs once per
            // (variant x transcript x allele) with `--hgvs`
            Some(v) if !has_version => {
                let mut out = String::with_capacity(id.len() + 11);
                out.push_str(id);
                out.push('.');
                let _ = write!(out, "{v}");
                out
            }
            _ => id.clone(),
        })
    }

    /// Total number of exons.
    pub fn exon_count(&self) -> usize {
        self.exons.len()
    }

    /// Total number of introns (exons - 1, minimum 0).
    pub fn intron_count(&self) -> usize {
        self.exons.len().saturating_sub(1)
    }

    /// Get the cDNA length (sum of exon lengths).
    pub fn cdna_length(&self) -> u64 {
        self.exons.iter().map(|e| e.end - e.start + 1).sum()
    }

    /// Map a genomic position to a cDNA position.
    /// Returns None if the position is not in an exon.
    pub fn genomic_to_cdna(&self, genomic_pos: u64) -> Option<u64> {
        let mut cdna_pos = 0u64;
        let sorted_exons = self.sorted_exons();

        for exon in sorted_exons.iter() {
            if genomic_pos >= exon.start && genomic_pos <= exon.end {
                return match self.strand {
                    Strand::Forward => Some(cdna_pos + (genomic_pos - exon.start) + 1),
                    Strand::Reverse => Some(cdna_pos + (exon.end - genomic_pos) + 1),
                };
            }
            cdna_pos += exon.end - exon.start + 1;
        }
        None
    }

    /// Map a cDNA position to a CDS position.
    /// Returns None if the position is not in the coding region.
    pub fn cdna_to_cds(&self, cdna_pos: u64) -> Option<u64> {
        let coding_start = self.cdna_coding_start?;
        let coding_end = self.cdna_coding_end?;
        if cdna_pos >= coding_start && cdna_pos <= coding_end {
            Some(cdna_pos - coding_start + 1 + self.codon_table_start_phase)
        } else {
            None
        }
    }

    /// Map a CDS position (1-based) to a protein position (1-based).
    pub fn cds_to_protein(cds_pos: u64) -> u64 {
        (cds_pos - 1) / 3 + 1
    }

    /// Length of the full-length peptide in residues, excluding the stop codon.
    ///
    /// Prefers the translated peptide when sequences have been built, and
    /// otherwise derives it from the coding cDNA span, so the value is
    /// available whether or not the run was given a FASTA. `None` for
    /// non-coding transcripts.
    pub fn peptide_length(&self) -> Option<u64> {
        if let Some(pep) = self.peptide.as_deref() {
            let residues = pep.trim_end_matches('*').chars().count() as u64;
            if residues > 0 {
                return Some(residues);
            }
        }
        let cds_len = self.cdna_coding_end?.checked_sub(self.cdna_coding_start?)?
            + 1
            + self.codon_table_start_phase;
        // A complete CDS ends in a stop codon, which is not a residue.
        (cds_len / 3).checked_sub(1).filter(|n| *n > 0)
    }

    /// Whether a premature termination codon at `cdna_pos` (1-based cDNA
    /// coordinate) is predicted to escape nonsense-mediated decay.
    ///
    /// This is the 50-nucleotide rule the ClinGen SVI's PVS1 decision tree
    /// runs on (Abou Tayoun 2018, Hum Mutat 39:1517): the surveillance complex
    /// is deposited at exon-exon junctions, so decay is triggered only when the
    /// PTC sits more than 50 nt upstream of the *final* junction. A PTC in the
    /// last exon, in the last 50 nt of the penultimate exon, or in a
    /// single-exon transcript therefore escapes.
    ///
    /// Returns `None` when the transcript has no exon structure to measure
    /// against, which the caller must treat as "unknown" rather than "escapes".
    pub fn escapes_nmd(&self, cdna_pos: u64) -> Option<bool> {
        let sorted = self.sorted_exons();
        let last = sorted.last()?;
        if sorted.len() < 2 {
            // No junction anywhere in the transcript, so nothing marks the
            // message for decay however early the PTC falls.
            return Some(true);
        }
        let last_exon_len = last.end.checked_sub(last.start)? + 1;
        // cDNA coordinate of the final base before the last exon-exon junction.
        let last_junction = self.cdna_length().checked_sub(last_exon_len)?;
        // Escape window is the 50 bases ending at the junction, plus everything
        // downstream of it.
        Some(cdna_pos + 50 > last_junction)
    }

    /// Build spliced_seq, translateable_seq, and peptide from a FASTA reference.
    ///
    /// The `fetch_seq` closure takes (chrom, start, end) using 1-based inclusive
    /// coordinates and returns the genomic sequence as uppercase bytes.
    pub fn build_sequences<F>(&mut self, fetch_seq: F) -> Result<(), String>
    where
        F: Fn(&str, u64, u64) -> Result<Vec<u8>, String>,
    {
        let sorted = self.sorted_exons_owned();

        // Build spliced sequence by concatenating exon sequences in transcript order
        let mut spliced = Vec::new();
        for exon in &sorted {
            let seq = fetch_seq(&self.chromosome, exon.start, exon.end)?;
            match self.strand {
                Strand::Forward => spliced.extend_from_slice(&seq),
                Strand::Reverse => spliced.extend_from_slice(&reverse_complement(&seq)),
            }
        }

        let spliced_str = String::from_utf8(spliced)
            .map_err(|e| format!("Invalid UTF-8 in spliced sequence: {}", e))?;
        self.spliced_seq = Some(spliced_str.clone());

        // Extract translateable sequence (coding portion of spliced cDNA)
        // When codon_table_start_phase > 0, CDS numbering includes bases before
        // cdna_coding_start. Start extraction earlier by the phase offset so that
        // CDS position N maps to translateable_seq index N-1.
        if let (Some(cs), Some(ce)) = (self.cdna_coding_start, self.cdna_coding_end) {
            let phase = self.codon_table_start_phase as usize;
            let cs_idx = (cs - 1) as usize;
            let ce_idx = ce as usize;
            // A coding region that starts after it ends has nothing to be read out of the spliced
            // sequence: exons that are out of order or overlap put it there, and slicing it
            // panics, which in the parallel build of a whole GFF3 takes the run down with it.
            if cs_idx > ce_idx {
                return Err(format!(
                    "the coding region of {} starts at cDNA {} and ends at {}, before it starts",
                    self.stable_id, cs, ce
                ));
            }
            if ce_idx <= spliced_str.len() {
                // For incomplete CDS starts (phase > 0), pad with N bases so
                // CDS position 1 maps to translateable index phase (aligning codons).
                let padding = "N".repeat(phase);
                let raw_translateable = &spliced_str[cs_idx..ce_idx];
                let translateable = format!("{}{}", padding, raw_translateable);
                self.translateable_seq = Some(translateable.to_string());

                // Translate to peptide. Mitochondrial transcripts use the
                // vertebrate mitochondrial genetic code (NCBI table 2) —
                // AGA/AGG are stop codons, ATA is Met, and TGA is Trp there,
                // unlike the standard nuclear table.
                let codon_table = if crate::mitochondrial::is_mitochondrial(&self.chromosome) {
                    crate::mitochondrial::mitochondrial_codon_table()
                } else {
                    CodonTable::standard()
                };
                let mut peptide_bytes = codon_table.translate_seq(translateable.as_bytes());
                resolve_readthrough_selenocysteine(
                    &mut peptide_bytes,
                    translateable.as_bytes(),
                    &codon_table,
                );
                // Mitochondrial genes are left alone: their initiators are `ATA` and `ATT` as
                // often as `ATG` and they have no UTR, which is what the reading below takes for
                // an unconfirmed start.
                if self.cds_flags_from_sequence
                    && !crate::mitochondrial::is_mitochondrial(&self.chromosome)
                {
                    self.read_ensembl_cds_state(translateable.as_bytes(), &codon_table);
                }
                // `Transcript::translate` writes a first codon the table lists as a start (`CTG` is
                // one, for the standard table) as `M`, and the translation's `initial_met` edit
                // does the same for the codons it is given to (`initiator_met`). They are two
                // things in Ensembl and they are two things here: the first is the protein
                // `hgvs_protein` reads a frameshift against (`p.Met1?` for a CTG start), the second
                // is also what `Amino_acids` reads for codon 1 (`resolve_readthrough_residue`).
                let first_codon_is_start = translateable
                    .as_bytes()
                    .first_chunk::<3>()
                    .is_some_and(|codon| codon_table.lists_start(codon));
                if first_codon_is_start || self.initiator_met {
                    if let Some(first) = peptide_bytes.first_mut() {
                        *first = b'M';
                    }
                }
                self.peptide = Some(String::from_utf8_lossy(&peptide_bytes).to_string());
            }
        }

        Ok(())
    }

    /// Set `flags` and `initiator_met` the way Ensembl's database has them, for a transcript whose
    /// GFF3 is Ensembl's export and so says neither (`cds_flags_from_sequence`).
    ///
    /// Read against VEP 116.3's own cache for chr21 and chr22 (12,867 coding transcripts, the
    /// Ensembl 116 GFF3 and the genome), the rules below give the cache's `cds_start_NF` on 12,859
    /// and its `cds_end_NF` on 12,858:
    ///
    /// - `cds_end_NF` when the coding sequence does not end in a stop codon in frame. All 621 the
    ///   cache flags are found, and 9 it does not flag are too (a CDS of 698 bases, say, that the
    ///   database calls complete).
    /// - `cds_start_NF` when the first CDS exon starts part-way through a codon (all 229 of
    ///   those carry the flag), or when its first codon is not `ATG` and there is no 5' UTR above
    ///   it. A transcript that begins with another codon *and* has a UTR is complete in the
    ///   cache (18 of 18); one without a UTR is not (128 of 130). The six it flags that start with
    ///   `ATG` cannot be told from the rest.
    /// - `initial_met`, the edit that makes residue 1 `M`, for a complete transcript that starts
    ///   with `CTG`, `GTG` or `ACG`. The cache gives it to all 14 complete transcripts that start
    ///   so, and not to the six complete ones that start with something else (`GCT`, `GAA`);
    ///   twelve of the 14 have a UTR and are found, the other two are read as unconfirmed above.
    ///   `TTG` is the fourth codon of the usual list, and none of the cache's complete transcripts
    ///   starts with it, so it is not assumed.
    fn read_ensembl_cds_state(&mut self, translateable: &[u8], table: &CodonTable) {
        let first = translateable.first_chunk::<3>();
        let has_five_prime_utr = self.cdna_coding_start.is_some_and(|start| start > 1);
        let start_unconfirmed =
            self.codon_table_start_phase > 0 || (first != Some(b"ATG") && !has_five_prime_utr);
        let end_in_frame_stop = translateable.len().is_multiple_of(3)
            && translateable
                .last_chunk::<3>()
                .is_some_and(|codon| table.is_stop(codon));
        if start_unconfirmed {
            self.flags.push("cds_start_NF".to_string());
        }
        if !end_in_frame_stop {
            self.flags.push("cds_end_NF".to_string());
        }
        if !start_unconfirmed && matches!(first, Some(b"CTG" | b"GTG" | b"ACG")) {
            self.initiator_met = true;
        }
    }

    /// Get exons sorted by transcript order, returning owned Exons.
    fn sorted_exons_owned(&self) -> Vec<Exon> {
        let mut exons: Vec<Exon> = self.exons.clone();
        match self.strand {
            Strand::Forward => exons.sort_by_key(|e| e.start),
            Strand::Reverse => exons.sort_by_key(|e| std::cmp::Reverse(e.start)),
        }
        exons
    }

    /// Get exons sorted by transcript order (forward for +, reverse for -).
    ///
    /// The GFF3 reader stores them that way, so this is a view of `exons` and not a copy that is
    /// sorted: it runs for every (variant x transcript x allele), a few times each, and the copy
    /// was an allocation and a sort of the exons of the transcript every time.
    pub(crate) fn sorted_exons(&self) -> ExonsInOrder<'_> {
        let exons = &self.exons[..];
        let ascending = exons.windows(2).all(|w| w[0].start <= w[1].start);
        let descending = exons.windows(2).all(|w| w[0].start >= w[1].start);
        let strictly_ascending = || exons.windows(2).all(|w| w[0].start < w[1].start);
        let strictly_descending = || exons.windows(2).all(|w| w[0].start > w[1].start);
        match self.strand {
            // a stable sort of what is already in order is what is in order, and of what is
            // strictly the other way round is the list read from its end
            Strand::Forward if ascending => ExonsInOrder::Stored(exons),
            Strand::Forward if strictly_descending() => ExonsInOrder::Backwards(exons),
            Strand::Reverse if descending => ExonsInOrder::Stored(exons),
            Strand::Reverse if strictly_ascending() => ExonsInOrder::Backwards(exons),
            _ => {
                let mut sorted: Vec<&Exon> = exons.iter().collect();
                match self.strand {
                    Strand::Forward => sorted.sort_by_key(|e| e.start),
                    Strand::Reverse => sorted.sort_by_key(|e| std::cmp::Reverse(e.start)),
                }
                ExonsInOrder::Sorted(sorted)
            }
        }
    }

    /// Determine which exon (0-indexed rank) a genomic position falls in.
    /// Returns (exon_index, total_exons) or None if intronic/outside.
    pub fn exon_at(&self, genomic_pos: u64) -> Option<(usize, usize)> {
        let sorted = self.sorted_exons();
        for (i, exon) in sorted.iter().enumerate() {
            if genomic_pos >= exon.start && genomic_pos <= exon.end {
                return Some((i, sorted.len()));
            }
        }
        None
    }

    /// Check if a genomic range [start, end] overlaps any exon.
    /// Returns the first overlapping exon info, or None.
    pub fn exon_overlapping(&self, range_start: u64, range_end: u64) -> Option<(usize, usize)> {
        let (start, end) = (range_start.min(range_end), range_start.max(range_end));
        let sorted = self.sorted_exons();
        for (i, exon) in sorted.iter().enumerate() {
            if start <= exon.end && end >= exon.start {
                return Some((i, sorted.len()));
            }
        }
        None
    }

    /// Map a genomic position in an intron to HGVSc offset notation.
    ///
    /// Returns `(nearest_exon_boundary_cdna_pos, signed_offset)`:
    /// - Positive offset: relative to donor end of upstream exon (e.g., c.151+5)
    /// - Negative offset: relative to acceptor start of downstream exon (e.g., c.152-3)
    pub fn genomic_to_intronic_cdna(&self, genomic_pos: u64) -> Option<(u64, i64)> {
        let sorted = self.sorted_exons();
        let n_introns = sorted.len().saturating_sub(1);

        for i in 0..n_introns {
            // Intron genomic range depends on strand-based sort order
            let (intron_start, intron_end) = match self.strand {
                Strand::Forward => (sorted[i].end + 1, sorted[i + 1].start - 1),
                Strand::Reverse => (sorted[i + 1].end + 1, sorted[i].start - 1),
            };

            if genomic_pos >= intron_start && genomic_pos <= intron_end {
                // "Upstream" in transcript order = sorted[i], "downstream" = sorted[i+1]
                // For forward: donor is sorted[i].end, acceptor is sorted[i+1].start
                // For reverse: donor is sorted[i].start (lower genomic), acceptor is sorted[i+1].end (higher genomic)
                let dist_from_donor = match self.strand {
                    Strand::Forward => genomic_pos - sorted[i].end,
                    Strand::Reverse => sorted[i].start - genomic_pos,
                };
                let dist_from_acceptor = match self.strand {
                    Strand::Forward => sorted[i + 1].start - genomic_pos,
                    Strand::Reverse => genomic_pos - sorted[i + 1].end,
                };

                if dist_from_donor <= dist_from_acceptor {
                    // Closer to donor (upstream exon end): c.X+offset
                    let donor_genomic = match self.strand {
                        Strand::Forward => sorted[i].end,
                        Strand::Reverse => sorted[i].start,
                    };
                    let donor_cdna = self.genomic_to_cdna(donor_genomic)?;
                    return Some((donor_cdna, dist_from_donor as i64));
                } else {
                    // Closer to acceptor (downstream exon start): c.X-offset
                    let acceptor_genomic = match self.strand {
                        Strand::Forward => sorted[i + 1].start,
                        Strand::Reverse => sorted[i + 1].end,
                    };
                    let acceptor_cdna = self.genomic_to_cdna(acceptor_genomic)?;
                    return Some((acceptor_cdna, -(dist_from_acceptor as i64)));
                }
            }
        }
        None
    }

    /// The signed offset from the nearest exon boundary that a genomic span
    /// reaches: `+N` on the donor side, `-N` on the acceptor side, `None` for a
    /// span that touches no intronic base.
    ///
    /// What is returned for a range is the smallest |offset| it *covers*, not
    /// the smaller of its two endpoints, matching
    /// `fastvep_core::parse_intronic_offset` - the two have to agree, because
    /// the criteria that read this also read offsets parsed out of ClinVar's
    /// own HGVS strings. A span with one end in an exon and the other in an
    /// intron runs continuously into it, so the nearest intronic base it
    /// reaches is the first one, `+1` or `-1`, whatever its far end reads.
    ///
    /// This is the same number the criteria used to recover by parsing a `+N` /
    /// `-N` token out of the HGVSc, and it is not the same answer: `c.` is a
    /// display form, 3'-shifted, so the offset can be shifted away entirely.
    pub fn intronic_offset_covered(&self, start: u64, end: u64) -> Option<i64> {
        let (lo, hi) = (start.min(end), start.max(end));
        let low = self.intronic_offset_at(lo);
        // One base reads the same at both ends, and almost every variant is one.
        let high = if hi == lo {
            low
        } else {
            self.intronic_offset_at(hi)
        };
        match (low, high) {
            (None, None) => None,
            (None, Some(o)) | (Some(o), None) => Some(o.signum()),
            (Some(x), Some(y)) => Some(if x.abs() <= y.abs() { x } else { y }),
        }
    }

    /// Where a genomic span sits against the exons and introns of this
    /// transcript, for a caller that needs to know whether it is *wholly* inside
    /// one of them.
    ///
    /// [`Self::intronic_offset_covered`] cannot say: it reads `+1` for a span
    /// that runs from an exon into an intron and the nearer end's offset for one
    /// that crosses an exon between two introns, and it reads `None` for both a
    /// span inside an exon and one that covers a whole intron. This separates
    /// the three outcomes a criterion has to tell apart. An insertion is the
    /// span of the two bases it sits between, so one that sits on an exon
    /// boundary is [`SpanSite::Spans`].
    ///
    /// Reads the exons in the order they are stored, as the other per-allele
    /// helpers do, and allocates nothing.
    pub fn span_site(&self, start: u64, end: u64) -> SpanSite {
        let (lo, hi) = (start.min(end), start.max(end));
        if self.exons.iter().any(|e| e.start <= lo && hi <= e.end) {
            return SpanSite::Exon;
        }
        match self.introns().find(|i| i.start <= lo && hi <= i.end) {
            Some(intron) => {
                let (a, b) = (intron.offset_of(lo), intron.offset_of(hi));
                SpanSite::Intron {
                    offset: if a.abs() <= b.abs() { a } else { b },
                }
            }
            None => SpanSite::Spans,
        }
    }

    /// Whether a change touches the first base of an exon or the last three, for
    /// an exon that has a splice site on that side.
    ///
    /// These are the positions Walker 2023 (ClinGen SVI Splicing Subgroup)
    /// excludes from BP7, for the reason BP7's own doc comment gives, that
    /// SpliceAI can miss splice impact there.
    /// The first base of a transcript's first exon and the last bases of its last
    /// are not next to a splice site, so they are not edges. Counted in
    /// transcript order, so on the reverse strand the first base of an exon is
    /// its higher coordinate.
    pub fn touches_splice_exon_edge(&self, start: u64, end: u64) -> bool {
        let (lo, hi) = (start.min(end), start.max(end));
        let lowest = self.exons.iter().map(|e| e.start).min().unwrap_or(0);
        let highest = self.exons.iter().map(|e| e.end).max().unwrap_or(0);
        self.exons.iter().any(|exon| {
            // The transcript's own first exon has no acceptor and its last no
            // donor; which of the outermost exons that is depends on the strand.
            let (has_acceptor, has_donor) = match self.strand {
                Strand::Forward => (exon.start != lowest, exon.end != highest),
                Strand::Reverse => (exon.end != highest, exon.start != lowest),
            };
            let (first_base, final_three) = match self.strand {
                Strand::Forward => (
                    (exon.start, exon.start),
                    (exon.end.saturating_sub(2).max(exon.start), exon.end),
                ),
                Strand::Reverse => (
                    (exon.end, exon.end),
                    (exon.start, (exon.start + 2).min(exon.end)),
                ),
            };
            let overlaps = |(a, b): (u64, u64)| lo <= b && a <= hi;
            (has_acceptor && overlaps(first_base)) || (has_donor && overlaps(final_three))
        })
    }

    /// The signed distance from `pos` to the nearest exon boundary of the intron
    /// it sits in, `None` for a position in no intron.
    ///
    /// The same number [`Self::genomic_to_intronic_cdna`] returns beside the cDNA
    /// anchor, without the anchor: that call maps the boundary through the
    /// spliced coordinates and, like every accessor that sorts, allocates doing
    /// it - per end, per allele, per overlapping transcript, to produce a value
    /// the only caller here then throws away.
    fn intronic_offset_at(&self, pos: u64) -> Option<i64> {
        self.introns()
            .find(|i| i.start <= pos && pos <= i.end)
            .map(|i| i.offset_of(pos))
    }

    /// The introns of this transcript.
    ///
    /// Adjacent exons bound one intron, whichever way the list runs, so pairing
    /// them directly - rather than through `sorted_exons` - is what keeps the
    /// callers allocation-free. The cost is the invariant that the list is
    /// ascending or descending, asserted here and shared with
    /// `fastvep-consequence`'s splice code, which makes the same bargain for the
    /// same reason (`gff.rs` stores exons in transcript order, which is one of
    /// the two). Which exon is the donor side is the strand's to say, not the
    /// list's: the lower exon is upstream on the forward strand and downstream on
    /// the reverse.
    fn introns(&self) -> impl Iterator<Item = IntronBounds> + '_ {
        debug_assert!(
            self.exons.windows(2).all(|p| p[0].start <= p[1].start)
                || self.exons.windows(2).all(|p| p[0].start >= p[1].start),
            "exons of {} are neither ascending nor descending",
            self.stable_id
        );
        self.exons.windows(2).filter_map(move |pair| {
            let (lower, upper) = if pair[0].start <= pair[1].start {
                (&pair[0], &pair[1])
            } else {
                (&pair[1], &pair[0])
            };
            let (start, end) = (lower.end + 1, upper.start.checked_sub(1)?);
            let (donor, acceptor) = match self.strand {
                Strand::Forward => (lower.end, upper.start),
                Strand::Reverse => (upper.start, lower.end),
            };
            // A transcript with overlapping exons yields start > end.
            (start <= end).then_some(IntronBounds {
                start,
                end,
                donor,
                acceptor,
                strand: self.strand,
            })
        })
    }

    /// Get the genomic boundaries (start, end) of the intron containing the given position.
    /// Returns None if the position is not in an intron.
    pub fn intron_bounds_at(&self, genomic_pos: u64) -> Option<(u64, u64)> {
        let sorted = self.sorted_exons();
        let n_introns = sorted.len().saturating_sub(1);
        for i in 0..n_introns {
            let (intron_start, intron_end) = match self.strand {
                Strand::Forward => (sorted[i].end + 1, sorted[i + 1].start - 1),
                Strand::Reverse => (sorted[i + 1].end + 1, sorted[i].start - 1),
            };
            if genomic_pos >= intron_start && genomic_pos <= intron_end {
                return Some((intron_start, intron_end));
            }
        }
        None
    }

    /// Determine which intron (0-indexed rank) a genomic position falls in.
    /// Returns (intron_index, total_introns) or None if exonic/outside.
    pub fn intron_at(&self, genomic_pos: u64) -> Option<(usize, usize)> {
        let sorted = self.sorted_exons();
        let n_introns = sorted.len().saturating_sub(1);
        for i in 0..n_introns {
            // For forward strand: sorted by ascending start, intron is between sorted[i].end and sorted[i+1].start
            // For reverse strand: sorted by descending start, intron is between sorted[i+1].end and sorted[i].start
            let (intron_start, intron_end) = match self.strand {
                Strand::Forward => (sorted[i].end + 1, sorted[i + 1].start - 1),
                Strand::Reverse => (sorted[i + 1].end + 1, sorted[i].start - 1),
            };
            if genomic_pos >= intron_start && genomic_pos <= intron_end {
                return Some((i, n_introns));
            }
        }
        None
    }
}

/// An exon within a transcript.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Exon {
    pub stable_id: String,
    pub start: u64,
    pub end: u64,
    pub strand: Strand,
    /// Reading frame phase at start of exon (-1, 0, 1, 2).
    pub phase: i8,
    /// Reading frame phase at end of exon.
    pub end_phase: i8,
    /// 1-based rank of the exon in the transcript.
    pub rank: u32,
}

impl Exon {
    pub fn length(&self) -> u64 {
        self.end - self.start + 1
    }
}

/// Translation metadata for a protein-coding transcript.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Translation {
    pub stable_id: String,
    /// Genomic start of the translation (start codon first base).
    pub genomic_start: u64,
    /// Genomic end of the translation (stop codon last base).
    pub genomic_end: u64,
    /// The exon rank where translation starts.
    pub start_exon_rank: u32,
    /// Offset within the start exon (0-based).
    pub start_exon_offset: u64,
    /// The exon rank where translation ends.
    pub end_exon_rank: u32,
    /// Offset within the end exon (0-based).
    pub end_exon_offset: u64,
}

#[cfg(test)]
mod tests {
    /// Translate then resolve, the way `build_sequences` does, so the tests
    /// exercise the same pairing of peptide and CDS.
    fn translated(cds: &str) -> String {
        let table = CodonTable::standard();
        let mut pep = table.translate_seq(cds.as_bytes());
        super::resolve_readthrough_selenocysteine(&mut pep, cds.as_bytes(), &table);
        String::from_utf8(pep).unwrap()
    }

    #[test]
    fn an_internal_tga_inside_a_complete_cds_is_selenocysteine() {
        // The shape of every real selenoprotein in Ensembl 115: a UGA partway
        // through a CDS the annotation carries on past, and a real stop at the
        // end. SELENOW is exactly this, with the UGA at residue 13 of 87.
        // ATG AAA TGA CGG TAA -> M K * R *  ->  M K U R *
        assert_eq!(translated("ATGAAATGACGGTAA"), "MKUR*");

        // Several of them, as in TXNRD2 and the multi-Sec transcripts.
        // ATG TGA AAA TGA CGG TAA
        assert_eq!(translated("ATGTGAAAATGACGGTAA"), "MUKUR*");

        // Case-insensitively, since soft-masked reference sequence is lowercase.
        assert_eq!(translated("atgaaatgacggtaa"), "MKUR*");
    }

    #[test]
    fn only_tga_reads_through_and_only_inside_a_complete_cds() {
        // TAA and TAG are not selenocysteine codons. Inside a CDS they signal an
        // annotation problem, and inventing a residue there would be a guess.
        // ATG AAA TAA CGG TAA / ATG AAA TAG CGG TAA
        assert_eq!(translated("ATGAAATAACGGTAA"), "MK*R*");
        assert_eq!(translated("ATGAAATAGCGGTAA"), "MK*R*");

        // The terminator the annotation promises is never itself read through,
        // even though it is a TGA.
        // ATG AAA CGG TGA
        assert_eq!(translated("ATGAAACGGTGA"), "MKR*");

        // No terminator at the end: the frame is uncorroborated, so an internal
        // TGA could equally be a CDS annotated in the wrong frame. Left alone.
        // ATG AAA TGA CGG CGG
        assert_eq!(translated("ATGAAATGACGGCGG"), "MK*RR");

        // A partial final codon means the annotated CDS is not a whole protein.
        // ATG AAA TGA CGG TA
        assert_eq!(translated("ATGAAATGACGGTA"), "MK*R");
    }

    #[test]
    fn the_mitochondrial_table_has_no_readthrough_to_resolve() {
        // NCBI table 2 reads TGA as Trp, so it never produces the internal
        // terminator this resolves, and AGA/AGG stops must not be touched.
        let table = crate::mitochondrial::mitochondrial_codon_table();
        let cds = "ATGAAATGACGGAGA";
        let mut pep = table.translate_seq(cds.as_bytes());
        let before = pep.clone();
        super::resolve_readthrough_selenocysteine(&mut pep, cds.as_bytes(), &table);
        assert_eq!(pep, before, "mitochondrial translation must be untouched");
        assert!(
            !pep.contains(&b'U'),
            "no selenocysteine in the mitochondrial table: {:?}",
            String::from_utf8_lossy(&pep)
        );
    }

    use super::*;

    /// The exons in transcript order are what a stable sort of them by start gives (ascending for
    /// the forward strand, descending for the reverse), however they were stored: in order, the
    /// other way round, shuffled, with starts that repeat.
    #[test]
    fn the_exons_in_transcript_order_are_the_sorted_exons_however_they_were_stored() {
        let mut seed = 0x2545_f491_4f6c_dd1du64;
        let mut next = move |bound: u64| {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            seed % bound
        };
        let mut stored_as_is = 0;
        let mut stored_backwards = 0;
        let mut sorted_copy = 0;
        for round in 0..4000 {
            let mut t = make_test_transcript();
            t.strand = if round % 2 == 0 {
                Strand::Forward
            } else {
                Strand::Reverse
            };
            let n = next(8) as usize;
            // starts from a small range, so that some repeat; `rank` says where each began
            let mut exons: Vec<Exon> = (0..n)
                .map(|rank| Exon {
                    stable_id: format!("E{rank}"),
                    start: 100 + next(12) * 10,
                    end: 0,
                    strand: t.strand,
                    phase: -1,
                    end_phase: -1,
                    rank: rank as u32,
                })
                .collect();
            for e in &mut exons {
                e.end = e.start + 5;
            }
            // the shapes a reader can leave them in
            match round % 4 {
                0 => exons.sort_by_key(|e| e.start),
                1 => exons.sort_by_key(|e| std::cmp::Reverse(e.start)),
                2 => exons.reverse(),
                _ => {}
            }
            t.exons = exons;
            let mut expected: Vec<&Exon> = t.exons.iter().collect();
            match t.strand {
                Strand::Forward => expected.sort_by_key(|e| e.start),
                Strand::Reverse => expected.sort_by_key(|e| std::cmp::Reverse(e.start)),
            }
            let got = t.sorted_exons();
            match got {
                ExonsInOrder::Stored(_) => stored_as_is += 1,
                ExonsInOrder::Backwards(_) => stored_backwards += 1,
                ExonsInOrder::Sorted(_) => sorted_copy += 1,
            }
            assert_eq!(got.len(), expected.len());
            let ranks: Vec<u32> = got.iter().map(|e| e.rank).collect();
            let wanted: Vec<u32> = expected.iter().map(|e| e.rank).collect();
            assert_eq!(ranks, wanted, "{:?} {:?}", t.strand, t.exons);
            assert_eq!(got.last().map(|e| e.rank), wanted.last().copied());
            for (i, e) in expected.iter().enumerate() {
                assert_eq!(got[i].rank, e.rank);
            }
        }
        // all three ways of answering were taken
        assert!(stored_as_is > 100 && stored_backwards > 100 && sorted_copy > 100);
    }

    fn make_test_transcript() -> Transcript {
        Transcript {
            stable_id: "ENST00000001".into(),
            version: None,
            gene: Gene {
                stable_id: "ENSG00000001".into(),
                symbol: Some("TESTGENE".into()),
                symbol_source: Some("HGNC".into()),
                hgnc_id: Some("HGNC:1".into()),
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
                    stable_id: "ENSE00000001".into(),
                    start: 1000,
                    end: 1200,
                    strand: Strand::Forward,
                    phase: -1,
                    end_phase: 0,
                    rank: 1,
                },
                Exon {
                    stable_id: "ENSE00000002".into(),
                    start: 2000,
                    end: 2300,
                    strand: Strand::Forward,
                    phase: 0,
                    end_phase: 1,
                    rank: 2,
                },
                Exon {
                    stable_id: "ENSE00000003".into(),
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
            cdna_coding_end: Some(952),
            coding_region_start: Some(1050),
            coding_region_end: Some(4500),
            spliced_seq: None,
            translateable_seq: None,
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

    #[test]
    fn test_is_coding() {
        let tr = make_test_transcript();
        assert!(tr.is_coding());
    }

    #[test]
    fn test_peptide_length_from_coding_span() {
        // cDNA coding span 51..=952 is 902 bases; the trailing stop codon is
        // not a residue, so 902 / 3 - 1 = 299.
        assert_eq!(make_test_transcript().peptide_length(), Some(299));
    }

    #[test]
    fn test_peptide_length_prefers_the_translated_peptide() {
        let mut tr = make_test_transcript();
        tr.peptide = Some("MKVLA*".to_string());
        assert_eq!(
            tr.peptide_length(),
            Some(5),
            "the stop codon is not a residue"
        );
    }

    #[test]
    fn test_peptide_length_is_none_for_noncoding() {
        let mut tr = make_test_transcript();
        tr.translation = None;
        tr.cdna_coding_start = None;
        tr.cdna_coding_end = None;
        assert_eq!(tr.peptide_length(), None);
    }

    #[test]
    fn test_escapes_nmd_50nt_rule() {
        // Exons are 201 + 301 + 1001 = 1503 cDNA bases, so the final
        // exon-exon junction sits at cDNA position 502.
        let tr = make_test_transcript();
        assert_eq!(
            tr.escapes_nmd(452),
            Some(false),
            "51 nt upstream of the junction decays"
        );
        assert_eq!(
            tr.escapes_nmd(453),
            Some(true),
            "50 nt upstream is the escape window"
        );
        assert_eq!(
            tr.escapes_nmd(502),
            Some(true),
            "the junction itself escapes"
        );
        assert_eq!(tr.escapes_nmd(900), Some(true), "the last exon escapes");
        assert_eq!(tr.escapes_nmd(10), Some(false), "far upstream decays");
    }

    #[test]
    fn test_escapes_nmd_single_exon_transcript() {
        // With no exon-exon junction the surveillance complex is never
        // deposited, so no PTC triggers decay however early it falls.
        let mut tr = make_test_transcript();
        tr.exons.truncate(1);
        assert_eq!(tr.escapes_nmd(10), Some(true));
    }

    #[test]
    fn test_escapes_nmd_on_the_reverse_strand() {
        // Transcript order is reversed, so the last exon is the 1000..=1200
        // one (201 bases) and the junction lands at 1503 - 201 = 1302.
        let mut tr = make_test_transcript();
        tr.strand = Strand::Reverse;
        for e in &mut tr.exons {
            e.strand = Strand::Reverse;
        }
        assert_eq!(tr.escapes_nmd(1252), Some(false));
        assert_eq!(tr.escapes_nmd(1253), Some(true));
    }

    #[test]
    fn test_exon_count() {
        let tr = make_test_transcript();
        assert_eq!(tr.exon_count(), 3);
        assert_eq!(tr.intron_count(), 2);
    }

    #[test]
    fn test_genomic_to_cdna_forward() {
        let tr = make_test_transcript();
        // Position in first exon
        assert_eq!(tr.genomic_to_cdna(1000), Some(1));
        assert_eq!(tr.genomic_to_cdna(1200), Some(201));
        // Position in second exon
        assert_eq!(tr.genomic_to_cdna(2000), Some(202));
        // Position in intron
        assert_eq!(tr.genomic_to_cdna(1500), None);
    }

    #[test]
    fn test_exon_at() {
        let tr = make_test_transcript();
        assert_eq!(tr.exon_at(1100), Some((0, 3)));
        assert_eq!(tr.exon_at(2100), Some((1, 3)));
        assert_eq!(tr.exon_at(4500), Some((2, 3)));
        assert_eq!(tr.exon_at(1500), None); // intron
    }

    #[test]
    fn test_intron_at() {
        let tr = make_test_transcript();
        assert_eq!(tr.intron_at(1500), Some((0, 2)));
        assert_eq!(tr.intron_at(3000), Some((1, 2)));
        assert_eq!(tr.intron_at(1100), None); // exon
    }

    /// The three-exon fixture read the other way round: the same coordinates, on
    /// the reverse strand, with the exons stored in transcript order, which is
    /// descending. Exon 1 is 4000-5000, exon 2 2000-2300, exon 3 1000-1200.
    fn reverse_test_transcript() -> Transcript {
        let mut tr = make_test_transcript();
        tr.strand = Strand::Reverse;
        tr.gene.strand = Strand::Reverse;
        tr.exons.reverse();
        for (i, e) in tr.exons.iter_mut().enumerate() {
            e.strand = Strand::Reverse;
            e.rank = i as u32 + 1;
        }
        tr
    }

    /// A reverse-strand transcript whose exons are stored ascending, as several
    /// fixtures across the workspace build them, reads the same as one stored in
    /// transcript order: which exon is the donor side is the strand's to say.
    #[test]
    fn the_exon_order_in_the_list_does_not_change_what_the_helpers_report() {
        let descending = reverse_test_transcript();
        let mut ascending = descending.clone();
        ascending.exons.reverse();
        for pos in (900..5100).step_by(7) {
            assert_eq!(
                ascending.intronic_offset_at(pos),
                descending.intronic_offset_at(pos),
                "{pos}"
            );
            assert_eq!(
                ascending.touches_splice_exon_edge(pos, pos),
                descending.touches_splice_exon_edge(pos, pos),
                "{pos}"
            );
        }
        assert_eq!(
            ascending.span_site(3999, 3999),
            descending.span_site(3999, 3999)
        );
    }

    /// The offset the allocation-free lookup reports is the one
    /// `genomic_to_intronic_cdna` always did, base for base, on both strands.
    #[test]
    fn the_allocation_free_intronic_offset_agrees_with_the_mapping_it_replaces() {
        for tr in [make_test_transcript(), reverse_test_transcript()] {
            for pos in 900..5100 {
                assert_eq!(
                    tr.intronic_offset_at(pos),
                    tr.genomic_to_intronic_cdna(pos).map(|(_, off)| off),
                    "strand {:?}, position {pos}",
                    tr.strand
                );
            }
        }
    }

    /// What a span is, against the exons and introns, on both strands. Introns
    /// are 1201-1999 and 2301-3999, whose first bases are `+1` from the exon
    /// above them in the forward transcript and `-1` from the exon below.
    #[test]
    fn a_span_is_inside_an_exon_inside_an_intron_or_neither() {
        let fwd = make_test_transcript();
        assert_eq!(fwd.span_site(1100, 1150), SpanSite::Exon);
        assert_eq!(fwd.span_site(1300, 1310), SpanSite::Intron { offset: 100 });
        // The nearer end decides inside one intron: 1990 is -10 from exon 2.
        assert_eq!(fwd.span_site(1500, 1990), SpanSite::Intron { offset: -10 });
        // A single base against a boundary.
        assert_eq!(fwd.span_site(1201, 1201), SpanSite::Intron { offset: 1 });
        assert_eq!(fwd.span_site(1999, 1999), SpanSite::Intron { offset: -1 });
        // An exon and an intron, an exon between two introns, a whole intron,
        // an insertion sitting on a boundary, and a change off the transcript.
        assert_eq!(fwd.span_site(1195, 1205), SpanSite::Spans);
        assert_eq!(fwd.span_site(1990, 2310), SpanSite::Spans);
        assert_eq!(fwd.span_site(1100, 2100), SpanSite::Spans);
        assert_eq!(fwd.span_site(1200, 1201), SpanSite::Spans);
        assert_eq!(fwd.span_site(10, 20), SpanSite::Spans);
        // An insertion is its two flanking bases, whichever order they arrive in.
        assert_eq!(fwd.span_site(1305, 1304), SpanSite::Intron { offset: 104 });

        let rev = reverse_test_transcript();
        assert_eq!(rev.span_site(4100, 4200), SpanSite::Exon);
        // Intron 1 is 2301-3999, entered from exon 1 (4000-5000) at 3999.
        assert_eq!(rev.span_site(3999, 3999), SpanSite::Intron { offset: 1 });
        assert_eq!(rev.span_site(2301, 2301), SpanSite::Intron { offset: -1 });
        assert_eq!(rev.span_site(3900, 3910), SpanSite::Intron { offset: 90 });
        assert_eq!(rev.span_site(3990, 4005), SpanSite::Spans);
        // 2310 is -10 from exon 2 and 3995 is +5 from exon 1: the nearer end.
        assert_eq!(rev.span_site(2310, 3995), SpanSite::Intron { offset: 5 });
    }

    /// Walker 2023: the first base of an exon and its last three, where there is
    /// a splice site on that side. Exon 2 here is 2000-2300 and has one on each.
    #[test]
    fn the_first_base_and_last_three_of_an_exon_are_its_splice_edges() {
        let fwd = make_test_transcript();
        let edge = |p: u64| fwd.touches_splice_exon_edge(p, p);
        assert!(edge(2000), "first base");
        assert!(!edge(2001), "second base is not an edge");
        assert!(!edge(2297), "fourth from the end is not an edge");
        assert!(edge(2298) && edge(2299) && edge(2300), "last three");
        // The first exon has no acceptor, the last no donor: their outer ends
        // are not next to a splice site.
        assert!(!edge(1000) && !edge(1001), "first base of the first exon");
        assert!(
            edge(1198) && edge(1200),
            "last bases of the first exon have a donor"
        );
        assert!(!edge(5000) && !edge(4998), "last bases of the last exon");
        assert!(edge(4000), "first base of the last exon has an acceptor");
        // A span counts if any base of it does.
        assert!(fwd.touches_splice_exon_edge(2290, 2298));
        assert!(!fwd.touches_splice_exon_edge(2100, 2200));

        // Reverse strand: transcript order puts an exon's first base at its
        // higher coordinate. Exon 2 is 2000-2300 again.
        let rev = reverse_test_transcript();
        let edge = |p: u64| rev.touches_splice_exon_edge(p, p);
        assert!(edge(2300), "first base, transcript order");
        assert!(!edge(2299));
        assert!(
            edge(2000) && edge(2001) && edge(2002),
            "last three, transcript order"
        );
        assert!(!edge(2003));
        assert!(!edge(5000), "first base of the first exon");
        assert!(!edge(1000), "last base of the last exon");
    }

    #[test]
    fn test_cds_to_protein() {
        assert_eq!(Transcript::cds_to_protein(1), 1);
        assert_eq!(Transcript::cds_to_protein(3), 1);
        assert_eq!(Transcript::cds_to_protein(4), 2);
        assert_eq!(Transcript::cds_to_protein(6), 2);
        assert_eq!(Transcript::cds_to_protein(7), 3);
    }

    #[test]
    fn test_build_sequences_forward() {
        let mut tr = make_test_transcript();
        // Mock FASTA: return known sequences for each exon
        let result = tr.build_sequences(|_chrom, start, end| {
            let len = (end - start + 1) as usize;
            match start {
                1000 => {
                    // Exon 1: 1000-1200 (201 bases)
                    // UTR: first 50 bases, CDS starts at offset 50
                    let mut seq = vec![b'N'; len];
                    // Put ATG at offset 50 (genomic 1050)
                    seq[50] = b'A';
                    seq[51] = b'T';
                    seq[52] = b'G';
                    // Fill rest of CDS with GCT (Ala)
                    for (i, base) in seq.iter_mut().enumerate().skip(53) {
                        *base = b"GCT"[(i - 53) % 3];
                    }
                    Ok(seq)
                }
                2000 => {
                    // Exon 2: 2000-2300 (301 bases), all CDS
                    let mut seq = vec![b'A'; len];
                    for (i, base) in seq.iter_mut().enumerate() {
                        *base = b"AAG"[i % 3]; // Lys
                    }
                    Ok(seq)
                }
                4000 => {
                    // Exon 3: 4000-5000 (1001 bases)
                    let mut seq = vec![b'N'; len];
                    // CDS: first 501 bases (4000-4500)
                    for (i, base) in seq.iter_mut().enumerate().take(501) {
                        *base = b"GGT"[i % 3]; // Gly
                    }
                    Ok(seq)
                }
                _ => Ok(vec![b'N'; len]),
            }
        });

        assert!(result.is_ok());
        assert!(tr.spliced_seq.is_some());
        assert!(tr.translateable_seq.is_some());
        assert!(tr.peptide.is_some());

        let ts = tr.translateable_seq.as_ref().unwrap();
        // translateable_seq should start with ATG
        assert!(
            ts.starts_with("ATG"),
            "translateable_seq starts with: {}",
            &ts[..6]
        );
        // Length should be cdna_coding_end - cdna_coding_start + 1 = 952 - 51 + 1 = 902
        // Actually: cdna_coding_end = 952 in the test fixture
        let expected_len =
            (tr.cdna_coding_end.unwrap() - tr.cdna_coding_start.unwrap() + 1) as usize;
        // Wait, the slice is [cs-1..ce] which is ce - (cs-1) = ce - cs + 1 items
        assert_eq!(ts.len(), expected_len);

        let peptide = tr.peptide.as_ref().unwrap();
        // Peptide should start with M (Met)
        assert!(
            peptide.starts_with('M'),
            "peptide starts with: {}",
            &peptide[..1]
        );
    }
    /// `HGVSp` names the protein as `id.version`, and not twice: Ensembl's GFF3 keeps the version
    /// beside the id, NCBI's in it.
    #[test]
    fn the_protein_accession_carries_one_version() {
        let accession = |id: Option<&str>, version: Option<u32>| {
            let mut tr = make_test_transcript();
            tr.protein_id = id.map(String::from);
            tr.protein_version = version;
            tr.versioned_protein_id()
        };
        // Ensembl: the version is beside the id
        assert_eq!(
            accession(Some("ENSP00000363557"), Some(3)).as_deref(),
            Some("ENSP00000363557.3")
        );
        assert_eq!(
            accession(Some("ENSP00000363557"), Some(12)).as_deref(),
            Some("ENSP00000363557.12")
        );
        // NCBI: it is in the id, and the version GFF3 parsing supplies (1 where a record has none)
        // must not be put after it. `NP_000001.2.1` is not an accession; `.1.1` was not one either.
        assert_eq!(
            accession(Some("NP_000001.2"), Some(1)).as_deref(),
            Some("NP_000001.2")
        );
        assert_eq!(
            accession(Some("NP_000001.1"), Some(1)).as_deref(),
            Some("NP_000001.1")
        );
        assert_eq!(
            accession(Some("NP_000001.12"), Some(1)).as_deref(),
            Some("NP_000001.12")
        );
        // no version known, none added; no protein, no accession
        assert_eq!(
            accession(Some("ENSP00000363557"), None).as_deref(),
            Some("ENSP00000363557")
        );
        assert_eq!(accession(None, Some(3)), None);
        // a dot that is not a version is not one
        assert_eq!(
            accession(Some("PROT.abc"), Some(2)).as_deref(),
            Some("PROT.abc.2")
        );
    }

    /// The three-exon fixture (`1000-1200`, `2000-2300`, `4000-5000`) with its third exon moved
    /// to start `gap` bases after the second ends, which makes the second intron `gap` long.
    fn with_second_intron(gap: u64) -> Transcript {
        let mut tr = make_test_transcript();
        tr.exons[2].start = 2301 + gap;
        tr.exons[2].end = tr.exons[2].start + 999;
        tr
    }

    #[test]
    fn an_intron_is_a_frameshift_intron_up_to_thirteen_bases() {
        // `abs($intron_end - $intron_start) <= 12` is an intron whose last base is 12 past its
        // first, which is 13 bases long: the comment beside that test says 12.
        assert!(with_second_intron(1).has_frameshift_intron());
        assert!(with_second_intron(13).has_frameshift_intron());
        assert!(!with_second_intron(14).has_frameshift_intron());
        assert!(!make_test_transcript().has_frameshift_intron());
        // the exons of a reverse-strand transcript are stored descending
        let mut reverse = with_second_intron(13);
        reverse.strand = Strand::Reverse;
        reverse.exons.reverse();
        assert!(reverse.has_frameshift_intron());
    }

    #[test]
    fn overlapping_exons_make_no_intron_to_be_short() {
        // A RefSeq transcript with exons that overlap has an "intron" with start > end, which
        // Ensembl skips before it asks how long it is.
        let mut tr = make_test_transcript();
        tr.exons[2].start = 2250;
        tr.exons[2].end = 3000;
        assert!(!tr.has_frameshift_intron());
    }

    #[test]
    fn the_exon_flag_widens_every_exon_only_in_a_transcript_with_a_frameshift_intron() {
        let plain = make_test_transcript();
        assert!(plain.overlaps_exon_flag(1200, 1200));
        assert!(!plain.overlaps_exon_flag(1201, 1201));
        // third exon at 2306 to 3305: the intron 2301-2305 is five bases
        let short = with_second_intron(5);
        assert!(short.overlaps_exon_flag(1212, 1212), "12 past exon 1");
        assert!(!short.overlaps_exon_flag(1213, 1213), "13 past exon 1");
        assert!(short.overlaps_exon_flag(1988, 1988), "12 before exon 2");
        assert!(!short.overlaps_exon_flag(1987, 1987), "13 before exon 2");
        // the far side of the same exon is widened too, and by the same amount
        assert!(short.overlaps_exon_flag(3317, 3317));
        assert!(!short.overlaps_exon_flag(3318, 3318));
    }

    #[test]
    fn an_insertion_touches_an_exon_when_either_base_it_sits_between_does() {
        // An insertion arrives as `end = start - 1`. `_overlapped_exons` sorts the pair, so the
        // base on either side counts, where "inside the exon" needs both.
        let tr = make_test_transcript();
        assert!(tr.overlaps_exon_flag(2000, 1999), "between 1999 and 2000");
        assert!(tr.overlaps_exon_flag(1201, 1200), "between 1200 and 1201");
        assert!(!tr.overlaps_exon_flag(1501, 1500), "between 1500 and 1501");
    }

    /// A one-exon transcript made of `utr5` bases, `cds` and ten bases of 3' UTR, built the way a
    /// run with a FASTA builds it.
    fn built(cds: &str, utr5: usize, phase: u64, from_sequence: bool, chrom: &str) -> Transcript {
        let seq = format!("{}{}{}", "C".repeat(utr5), cds, "C".repeat(10));
        let mut t = make_test_transcript();
        t.chromosome = chrom.into();
        t.exons = vec![Exon {
            stable_id: "E1".into(),
            start: 1,
            end: seq.len() as u64,
            strand: Strand::Forward,
            phase: -1,
            end_phase: -1,
            rank: 1,
        }];
        t.start = 1;
        t.end = seq.len() as u64;
        t.cdna_coding_start = Some(utr5 as u64 + 1);
        t.cdna_coding_end = Some((utr5 + cds.len()) as u64);
        t.codon_table_start_phase = phase;
        t.flags = vec![];
        t.initiator_met = false;
        t.cds_flags_from_sequence = from_sequence;
        t.build_sequences(|_chrom, start, end| {
            Ok(seq.as_bytes()[start as usize - 1..end as usize].to_vec())
        })
        .unwrap();
        t
    }

    /// Ensembl's GFF3 export has no `cds_start_NF` or `cds_end_NF`, and the cache VEP reads them
    /// from does. Read off the sequence they agree with it on 12,858 of the 12,867 coding
    /// transcripts of chr21 and chr22 (the end) and 12,859 (the start): see
    /// [`Transcript::read_ensembl_cds_state`] for the rules and what they miss.
    #[test]
    fn the_cds_flags_of_an_ensembl_export_are_read_from_the_sequence() {
        // complete: ATG, a 5' UTR, a stop codon in frame at the end
        let t = built("ATGAAATGGTAA", 5, 0, true, "chr1");
        assert!(t.flags.is_empty(), "{:?}", t.flags);

        // no stop codon at the end
        let t = built("ATGAAATGGTTT", 5, 0, true, "chr1");
        assert_eq!(t.flags, ["cds_end_NF"]);

        // ends in a stop codon that is not in frame: the CDS is 13 bases
        let t = built("ATGAAATGGCTAA", 5, 0, true, "chr1");
        assert_eq!(
            t.flags,
            ["cds_end_NF"],
            "a stop out of frame is not the end"
        );

        // the first exon of the CDS starts part-way through a codon (two bases of it are missing,
        // so the 7 that are there make nine with them)
        let t = built("GTGGTAA", 5, 2, true, "chr1");
        assert_eq!(t.flags, ["cds_start_NF"]);

        // a first codon that is not ATG with nothing above it...
        let t = built("GCTAAATGGTAA", 0, 0, true, "chr1");
        assert_eq!(t.flags, ["cds_start_NF"]);
        // ...and with a UTR above it, which is how the cache has the 18 that are complete
        let t = built("GCTAAATGGTAA", 5, 0, true, "chr1");
        assert!(t.flags.is_empty(), "{:?}", t.flags);

        // both, in VEP's order
        let t = built("GCTAAATGGTTT", 0, 0, true, "chr1");
        assert_eq!(t.flags, ["cds_start_NF", "cds_end_NF"]);

        // a transcript whose GFF3 carries its own flags is left to them
        let t = built("GCTAAATGGTTT", 0, 0, false, "chr1");
        assert!(t.flags.is_empty(), "{:?}", t.flags);

        // mitochondrial genes have no UTR and start with ATA and ATT as often as ATG
        let t = built("ATTAAATGGTAA", 0, 0, true, "chrM");
        assert!(t.flags.is_empty(), "{:?}", t.flags);
    }

    /// Residue 1 is `M` in two ways and Ensembl VEP keeps them apart. `Transcript::translate`
    /// writes a first codon the code table lists as a start (`TTG`, `CTG` and `ATG` for table 1) as
    /// `M`, which is the protein `hgvs_protein` reads a frameshift against; the translation's
    /// `initial_met` edit makes `M` of the codons it is given to (`CTG`, `GTG` and `ACG` of the
    /// complete transcripts of chr21 and chr22), and that is what `Amino_acids` reads.
    #[test]
    fn residue_one_is_met_by_the_start_codon_rule_and_by_the_initial_met_edit() {
        // complete, with a UTR: the edit applies, and the table's rule too for a CTG
        let t = built("CTGAAATGGTAA", 5, 0, true, "chr1");
        assert!(t.initiator_met);
        assert_eq!(t.peptide.as_deref(), Some("MKW*"));
        let t = built("GTGAAATGGTAA", 5, 0, true, "chr1");
        assert!(
            t.initiator_met,
            "GTG is not in the table's list; the edit is"
        );
        assert_eq!(t.peptide.as_deref(), Some("MKW*"));
        let t = built("ACGAAATGGTAA", 5, 0, true, "chr1");
        assert!(t.initiator_met);
        assert_eq!(t.peptide.as_deref(), Some("MKW*"));

        // not a codon the edit is given to
        let t = built("GCTAAATGGTAA", 5, 0, true, "chr1");
        assert!(!t.initiator_met);
        assert_eq!(t.peptide.as_deref(), Some("AKW*"));

        // incomplete at the start: no edit, but `translate` still writes a CTG as M (the cache
        // has `pep1=M` and `Amino_acids` `L` for these)
        let t = built("CTGAAATGGTAA", 0, 0, true, "chr1");
        assert!(!t.initiator_met);
        assert_eq!(t.peptide.as_deref(), Some("MKW*"));
        let t = built("GTGAAATGGTAA", 0, 0, true, "chr1");
        assert!(!t.initiator_met);
        assert_eq!(
            t.peptide.as_deref(),
            Some("VKW*"),
            "GTG is not a start of table 1"
        );

        // RefSeq says so itself with `transl_except ... aa:Met`
        let mut t = make_test_transcript();
        t.initiator_met = true;
        t.cds_flags_from_sequence = false;
        let seq = "GTGAAATGGTAA";
        t.cdna_coding_start = Some(1);
        t.cdna_coding_end = Some(12);
        t.codon_table_start_phase = 0;
        t.exons = vec![Exon {
            stable_id: "E1".into(),
            start: 1,
            end: 12,
            strand: Strand::Forward,
            phase: 0,
            end_phase: -1,
            rank: 1,
        }];
        t.build_sequences(|_c, s, e| Ok(seq.as_bytes()[s as usize - 1..e as usize].to_vec()))
            .unwrap();
        assert_eq!(t.peptide.as_deref(), Some("MKW*"));
    }
}
