//! The structural variants of a run: which transcripts Ensembl VEP pairs a variant with, and what
//! it says of each pair. Both per-variant loops (`fastvep annotate` and the library the web server
//! runs) call this, so that a structural variant gets one answer whichever way it is asked.
//!
//! The rules are in `fastvep_consequence::sv_predictor`; what is here is the lookup they need that
//! only a transcript provider can do, and the rows of a variant that no transcript is near.

use fastvep_cache::providers::{SequenceProvider, TranscriptProvider};
use fastvep_consequence::sv_predictor::{predict_sv, Contigs, SvPrediction, SvTranscript};
use fastvep_core::{
    Allele, AlleleNo, Consequence, GenomicPosition, Impact, Strand, StructuralVariant,
    TranscriptLocation,
};
use fastvep_genome::Transcript;
use fastvep_io::variant::{AlleleAnnotation, TranscriptVariation, VariationFeature};
use fastvep_io::vcf::small_variant;

/// The transcripts a structural variant is paired with (`annotate_InputBuffer` of
/// `AnnotationType/Transcript.pm`): those whose `--distance` window holds the variant's interval
/// or, for a breakend, the position of one of its mates, on the contig that mate names.
///
/// Each is found by asking the provider for the contig by the name the record gives it, and is kept
/// with the number of that contig ([`Contigs`]), so that no name of a transcript is compared with
/// the record's.
pub fn paired_transcripts<'a>(
    provider: &'a dyn TranscriptProvider,
    vf: &VariationFeature,
    sv: &StructuralVariant,
    contigs: &Contigs,
    distance: u64,
) -> Vec<SvTranscript<'a>> {
    let mut found: Vec<SvTranscript<'a>> = Vec::new();
    let mut ask = |chrom: &str, from: u64, to: u64, contig: usize| {
        let window_start = from.saturating_sub(distance).max(1);
        let window_end = to.saturating_add(distance);
        for transcript in provider
            .get_transcripts(chrom, window_start, window_end)
            .unwrap_or_default()
        {
            if !found.iter().any(|f| std::ptr::eq(f.transcript, transcript)) {
                found.push(SvTranscript { transcript, contig });
            }
        }
    };
    ask(
        &vf.position.chromosome,
        vf.position.start,
        vf.position.end,
        0,
    );
    for (chrom, position) in sv.mates() {
        ask(chrom, position, position, contigs.index_of(chrom));
    }
    found
}

/// The transcripts a variant is annotated against, and for a structural variant what Ensembl VEP says
/// of each. Both per-variant loops ask this, so that a variant is paired with the same transcripts
/// whichever way it is asked.
pub struct Overlap<'a> {
    /// The transcripts, in the order of the rows: for a small variant those within `--distance` of
    /// it, and for a structural variant those it is paired with ([`paired_transcripts`]).
    pub transcripts: Vec<&'a Transcript>,
    /// Set for a structural variant that has a transcript: what VEP says of each of `transcripts`.
    /// A structural variant with none is intergenic ([`annotate_sv_intergenic`]).
    pub predictions: Option<Vec<SvPrediction<'a>>>,
}

pub fn overlap<'a>(
    provider: &'a dyn TranscriptProvider,
    vf: &VariationFeature,
    distance: u64,
) -> Overlap<'a> {
    if vf.sv.is_some() {
        return match annotate_sv(provider, vf, distance) {
            Some(predictions) => Overlap {
                transcripts: predictions.iter().map(|p| p.transcript).collect(),
                predictions: Some(predictions),
            },
            None => Overlap {
                transcripts: Vec::new(),
                predictions: None,
            },
        };
    }
    let start = vf.position.start.saturating_sub(distance).max(1);
    let end = vf.position.end.saturating_add(distance);
    Overlap {
        transcripts: provider
            .get_transcripts(&vf.position.chromosome, start, end)
            .unwrap_or_default(),
        predictions: None,
    }
}

/// What Ensembl VEP says of a structural variant against the transcripts around it: `None` when no
/// transcript is paired with it (the variant is intergenic).
pub fn annotate_sv<'a>(
    provider: &'a dyn TranscriptProvider,
    vf: &VariationFeature,
    distance: u64,
) -> Option<Vec<SvPrediction<'a>>> {
    let sv = vf.sv.as_deref()?;
    let contigs = Contigs::of(&vf.position.chromosome, sv);
    let found = paired_transcripts(provider, vf, sv, &contigs, distance);
    if found.is_empty() {
        return None;
    }
    Some(predict_sv(
        sv,
        &vf.alt_alleles,
        vf.position.start,
        vf.position.end,
        &contigs,
        &found,
        distance,
    ))
}

/// The rows of a structural variant that no transcript is paired with
/// (`IntergenicStructuralVariation`): one `intergenic_variant` for each of its alleles, a
/// breakend's own and each mate's included, with no feature.
pub fn annotate_sv_intergenic(vf: &mut VariationFeature) {
    let Some(sv) = vf.sv.as_deref() else {
        return;
    };
    let rows: Vec<TranscriptVariation> = sv
        .alleles
        .iter()
        .enumerate()
        // an ALT that asserts nothing (`.`, `<NON_REF>`) has no row here, as it has none over a
        // transcript
        .filter(|(_, allele)| {
            !matches!(vf.alt_alleles.get(allele.alt), None | Some(Allele::Missing))
        })
        .map(|(number, allele)| TranscriptVariation {
            transcript_id: TranscriptVariation::NO_FEATURE.into(),
            gene_id: "-".into(),
            gene_symbol: None,
            biotype: "-".into(),
            allele_annotations: vec![AlleleAnnotation {
                allele: vf
                    .alt_alleles
                    .get(allele.alt)
                    .cloned()
                    .unwrap_or(Allele::Missing),
                sv_allele: AlleleNo::of(number),
                consequences: vec![Consequence::IntergenicVariant],
                impact: Impact::Modifier,
                splice: Default::default(),
                location: TranscriptLocation::default(),
                protein_position: None,
                amino_acids: None,
                codons: None,
                exon: None,
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
                pick: false,
            }],
            canonical: false,
            strand: Strand::Forward,
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
        })
        .collect();
    if !rows.is_empty() {
        vf.most_severe_consequence = Some(Consequence::IntergenicVariant);
    }
    vf.transcript_variations.extend(rows);
}

/// A tandem repeat that has the INFO to be read as bases (`RUS` and `RUC` or `RB`) is annotated as the
/// small variant Ensembl VEP makes of it (`create_StructuralVariationFeatures`): the reference is the
/// bases of its interval, its alleles are the repeat unit, repeated, and one that is not as long as
/// the reference is trimmed of what it shares with it. Without a reference to read, or when a repeat
/// is longer than `max_sv_size`, it stays the structural variant it was.
///
/// A record that names no repeat (a `<DUP:TANDEM>` with the INFO but without `RN`, which counts only
/// the `CNV:TR` alleles) or whose repeats are the reference has no allele to write a row for, and
/// has none here either.
pub fn expand_tandem_repeat(
    vf: &mut VariationFeature,
    sequence: &dyn SequenceProvider,
    max_sv_size: Option<u64>,
) {
    let Some(tandem) = vf.sv.as_deref().and_then(|sv| sv.tandem.as_ref()) else {
        return;
    };
    let (start, end) = (vf.position.start, vf.position.end);
    let Some(repeats) = tandem.alleles(max_sv_size) else {
        return;
    };
    let Ok(bases) = sequence.fetch_sequence(&vf.position.chromosome, start, end) else {
        return;
    };
    if end < start || bases.len() as u64 != end - start + 1 {
        return;
    }
    let reference = String::from_utf8_lossy(&bases).to_ascii_uppercase();
    let small = small_variant(start, &reference, &repeats);
    vf.position = GenomicPosition::new(
        vf.position.chromosome.clone(),
        small.start,
        small.end,
        Strand::Forward,
    );
    vf.ref_allele = Allele::from_str(&small.reference);
    vf.alt_alleles = small.alt_alleles();
    vf.allele_string = small.allele_string;
    vf.variant_type = small.variant_type;
    vf.sv = None;
}
