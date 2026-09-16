//! The HGVSc normalisation path: 3'-shifting an indel over the genome, and
//! naming a shifted insertion as the duplication it is.
//!
//! [`hgvsc_shifted`] is the entry point, and both annotation loops call
//! it. They used to carry a copy each - the CLI's shifted, the library's did
//! not - so the same intronic duplication came out normalised from
//! `fastvep annotate` and unnormalised from the server.

use fastvep_cache::providers::SequenceProvider;

/// Convert intronic insertion to dup notation with explicit start/end positions
/// (coding).
///
/// Each end carries its own `(exon anchor, offset)` pair. A span deep inside one
/// intron shares an anchor, but one running past the intron's midpoint is
/// written from the exon on either side - `c.5044+27_5045-47dup` - so the two
/// ends are not interchangeable.
pub fn convert_ins_to_dup_range(
    hgvsc: &str,
    start: (u64, i64),
    end: (u64, i64),
    coding_start: u64,
    coding_end: Option<u64>,
) -> Option<String> {
    let prefix_end = hgvsc
        .find(":c.")
        .map(|i| i + 3)
        .or_else(|| hgvsc.find(":n.").map(|i| i + 3))?;
    let prefix = &hgvsc[..prefix_end];

    let build_pos = |cdna: u64, off: i64| -> String {
        let raw = cdna as i64 - coding_start as i64 + 1;
        let cp = if raw <= 0 { raw - 1 } else { raw };
        // An offset of 0 is an exonic end, written as the anchor alone. Folding
        // it into the negative arm renders `c.21` as `c.210`.
        let anchor = match coding_end.filter(|&ce| cp >= 0 && cdna > ce) {
            Some(ce) => format!("*{}", cdna - ce),
            None => format!("{}", cp),
        };
        match off.cmp(&0) {
            std::cmp::Ordering::Greater => format!("{}+{}", anchor, off),
            std::cmp::Ordering::Less => format!("{}{}", anchor, off),
            std::cmp::Ordering::Equal => anchor,
        }
    };

    if start == end {
        Some(format!("{}{}dup", prefix, build_pos(start.0, start.1)))
    } else {
        Some(format!(
            "{}{}_{}dup",
            prefix,
            build_pos(start.0, start.1),
            build_pos(end.0, end.1)
        ))
    }
}

/// Convert intronic insertion to dup notation with explicit start/end positions
/// (non-coding). See [`convert_ins_to_dup_range`] for why each end carries its
/// own anchor.
pub fn convert_ins_to_dup_range_noncoding(
    hgvsc: &str,
    start: (u64, i64),
    end: (u64, i64),
) -> Option<String> {
    let prefix_end = hgvsc
        .find(":n.")
        .map(|i| i + 3)
        .or_else(|| hgvsc.find(":c.").map(|i| i + 3))?;
    let prefix = &hgvsc[..prefix_end];

    let build_pos = |cdna: u64, off: i64| -> String {
        match off.cmp(&0) {
            std::cmp::Ordering::Greater => format!("{}+{}", cdna, off),
            std::cmp::Ordering::Less => format!("{}{}", cdna, off),
            std::cmp::Ordering::Equal => format!("{}", cdna),
        }
    };

    if start == end {
        Some(format!("{}{}dup", prefix, build_pos(start.0, start.1)))
    } else {
        Some(format!(
            "{}{}_{}dup",
            prefix,
            build_pos(start.0, start.1),
            build_pos(end.0, end.1)
        ))
    }
}

/// A window of reference the shift walks through, refilled in blocks rather than
/// fetched a base at a time.
///
/// The walk used to call `fetch_sequence` once per base, and twice per base for
/// a deletion, each returning an owned `Vec` - so a variant sliding 4,000 bases
/// down a repeat cost 4,000 heap allocations, or 8,000 as a deletion, once per
/// (variant x transcript x allele).
///
/// The block **grows**, and that is the point: almost every shift stops within a
/// base or two, so a window that opened at its full size would read and copy
/// hundreds of bases to answer two questions - work done before knowing it is
/// needed, which measured 11 % slower over a real indel-only callset than the
/// per-base reads it replaced. Starting small and doubling keeps the common case
/// at one short read and still collapses a long walk to a handful.
struct RefWindow<'a> {
    provider: &'a dyn SequenceProvider,
    chrom: &'a str,
    /// 1-based genomic position of `bases[0]`. Zero until the first fill.
    origin: u64,
    bases: Vec<u8>,
    next_block: u64,
}

impl<'a> RefWindow<'a> {
    /// Enough for a shift that stops immediately, which is nearly all of them.
    const FIRST_BLOCK: u64 = 16;
    /// Past this a walk is in a long repeat and the reads are already amortised.
    const MAX_BLOCK: u64 = 1024;

    fn new(provider: &'a dyn SequenceProvider, chrom: &'a str) -> Self {
        Self {
            provider,
            chrom,
            origin: 0,
            bases: Vec::new(),
            next_block: Self::FIRST_BLOCK,
        }
    }

    /// The base at `pos`, uppercased, or `None` past the contig.
    fn base(&mut self, pos: u64) -> Option<u8> {
        if pos == 0 {
            return None;
        }
        if self.origin == 0 || pos < self.origin || pos >= self.origin + self.bases.len() as u64 {
            let block = self.next_block;
            self.next_block = (block * 2).min(Self::MAX_BLOCK);
            // Centre the block on the request so a walk in either direction has
            // room; the caller's direction is not known here.
            let origin = pos.saturating_sub(block / 2).max(1);
            self.bases = self
                .provider
                .fetch_sequence_slice(self.chrom, origin, origin + block - 1)
                .ok()?;
            self.origin = origin;
            if self.bases.is_empty() {
                return None;
            }
        }
        self.bases
            .get((pos - self.origin) as usize)
            .map(|b| b.to_ascii_uppercase())
    }
}

/// 3' shift an indel along the transcript direction, over the reference.
///
/// HGVS describes a change at the most 3' position its reference sequence
/// allows, and the reference a `c.` description is numbered against is the
/// gene's genomic sequence - introns included. So the walk is over the genome in
/// transcript orientation, and `bound_start`/`bound_end` are the caller's
/// limits on it: in production the transcript's own extent, so a change can
/// travel out of the intron it was written in and into the exon beyond, which is
/// where HGVS puts it and where VEP writes it.
///
/// Bounding the walk by the *intron* instead was worth 529 of the 602 HGVSc rows
/// disagreeing with real VEP 115.1 over a 6,600-variant ClinVar sample: a
/// deletion of the last intronic base came out `c.1792-1del` where the exonic
/// base it repeats makes it `c.1793del`.
///
/// Returns the shifted genomic start and end positions.
// Each argument is an independent coordinate, allele or flag with no
// natural grouping; bundling them into a struct would only move the
// argument list to the call site.
#[allow(clippy::too_many_arguments)]
pub fn three_prime_shift_genomic(
    seq_provider: &dyn SequenceProvider,
    chrom: &str,
    start: u64,
    end: u64,
    ref_allele: &fastvep_core::Allele,
    alt_allele: &fastvep_core::Allele,
    strand: fastvep_core::Strand,
    bound_start: u64,
    bound_end: u64,
) -> (u64, u64) {
    use fastvep_core::Allele;

    match (ref_allele, alt_allele) {
        // A deletion slides 3' while the base past its end repeats the base at
        // its start, which leaves the deleted sequence unchanged.
        //
        // The two cursors sit a whole deletion apart, so they get a window each:
        // sharing one would refetch on every step of any deletion longer than a
        // block, which is worse than the per-base reads this replaced.
        (Allele::Sequence(ref_bases), Allele::Deletion) if !ref_bases.is_empty() => {
            let (mut s, mut e) = (start, end);
            let mut ahead = RefWindow::new(seq_provider, chrom);
            let mut behind = RefWindow::new(seq_provider, chrom);
            match strand {
                fastvep_core::Strand::Forward => loop {
                    let next = e + 1;
                    if next > bound_end {
                        break;
                    }
                    match (ahead.base(next), behind.base(s)) {
                        (Some(a), Some(b)) if a == b => {
                            s += 1;
                            e += 1;
                        }
                        _ => break,
                    }
                },
                fastvep_core::Strand::Reverse => loop {
                    if s == 0 || s - 1 < bound_start {
                        break;
                    }
                    match (ahead.base(s - 1), behind.base(e)) {
                        (Some(a), Some(b)) if a == b => {
                            s -= 1;
                            e -= 1;
                        }
                        _ => break,
                    }
                },
            }
            (s, e)
        }
        // An insertion slides over a base that repeats the next base of the
        // inserted sequence, which rotates the sequence by one each step.
        (Allele::Deletion, Allele::Sequence(ins_bases)) if !ins_bases.is_empty() => {
            let ins_len = ins_bases.len();
            let mut pos = start;
            let mut shift = 0usize;
            let mut window = RefWindow::new(seq_provider, chrom);
            match strand {
                fastvep_core::Strand::Forward => loop {
                    if pos > bound_end {
                        break;
                    }
                    let expected = ins_bases[shift % ins_len].to_ascii_uppercase();
                    match window.base(pos) {
                        Some(b) if b == expected => {
                            pos += 1;
                            shift += 1;
                        }
                        _ => break,
                    }
                },
                fastvep_core::Strand::Reverse => loop {
                    if pos == 0 || pos - 1 < bound_start {
                        break;
                    }
                    let expected = ins_bases[ins_len - 1 - (shift % ins_len)].to_ascii_uppercase();
                    match window.base(pos - 1) {
                        Some(b) if b == expected => {
                            pos -= 1;
                            shift += 1;
                        }
                        _ => break,
                    }
                },
            }
            (pos, pos.saturating_sub(1))
        }
        _ => (start, end),
    }
}

/// The block a 3'-shifted intronic insertion duplicates, in genomic coordinates.
///
/// After a maximal 3'-shift the duplicated copy can only sit immediately 5' of
/// the insertion point *in transcript orientation*; anything further 3' would
/// have been shifted over. The inserted string rotates one base per position
/// shifted, so the block is compared against that rotated form and not against
/// the string the VCF carried.
///
/// `shifted_start` follows the same convention [`three_prime_shift_genomic`]
/// returns: the insertion sits between genomic `shifted_start - 1` and
/// `shifted_start`.
///
/// The dup anchor used to be derived by re-shifting the *unshifted* insertion
/// point through `three_prime_shift_genomic` over a single position, which
/// walks while the next base repeats the current one - a homopolymer test. A
/// `TG` insertion in a `TGTGTG…` repeat therefore never moved at all and named
/// the copy 16 bases 5' of the one HGVS asks for. That accounted for 1,390 of
/// the 1,538 HGVSc rows disagreeing with real Ensembl VEP 115.1 over a
/// genome-wide HG002 sample, and it is why the two tools split by strand: this
/// walk stays where the VCF put the variant, so fastVEP's anchor was always the
/// genomically-left one whichever way the transcript ran.
pub fn intronic_dup_span(
    seq_provider: &dyn SequenceProvider,
    chrom: &str,
    shifted_start: u64,
    ins_bases: &[u8],
    shift: u64,
    strand: fastvep_core::Strand,
) -> Option<(u64, u64)> {
    let len = ins_bases.len();
    if len == 0 {
        return None;
    }
    // One base of rotation per position travelled, in the direction of travel:
    // shifting 3' on the forward strand walks the first base to the end, and on
    // the reverse strand - where 3' runs towards lower coordinates - the last
    // base walks to the front.
    let rot = (shift % len as u64) as usize;
    let mut rotated = ins_bases.to_vec();
    match strand {
        fastvep_core::Strand::Forward => rotated.rotate_left(rot),
        fastvep_core::Strand::Reverse => rotated.rotate_right(rot),
    }

    // 5' of the insertion point is genomically before it on the forward strand
    // and after it on the reverse.
    let (lo, hi) = match strand {
        fastvep_core::Strand::Forward => {
            let hi = shifted_start.checked_sub(1)?;
            (hi.checked_sub(len as u64 - 1)?, hi)
        }
        fastvep_core::Strand::Reverse => (shifted_start, shifted_start + len as u64 - 1),
    };
    if lo == 0 {
        return None;
    }

    let block = seq_provider.fetch_sequence_slice(chrom, lo, hi).ok()?;
    if block.len() != len {
        return None;
    }
    block
        .iter()
        .zip(rotated.iter())
        .all(|(a, b)| a.eq_ignore_ascii_case(b))
        .then_some((lo, hi))
}

/// Is this change one the HGVS 3'-rule can move - a pure deletion or a pure
/// insertion?
///
/// A substitution has nowhere to go, and neither tool shifts a delins, so the
/// callers use this to decide what routes through [`hgvsc_shifted`].
pub fn is_shiftable_indel(
    hgvs_ref: &fastvep_core::Allele,
    hgvs_alt: &fastvep_core::Allele,
) -> bool {
    use fastvep_core::Allele;
    matches!(
        (hgvs_ref, hgvs_alt),
        (Allele::Sequence(_), Allele::Deletion) | (Allele::Deletion, Allele::Sequence(_))
    )
}

/// One allele of a site, clipped of the bases its reference and alternate
/// repeat at either end and turned to face the transcript.
///
/// A VCF record carries one position for the whole site, so a multi-allelic
/// indel can only have the one base *every* allele shares stripped from it -
/// that is what the reader does, and it is what Ensembl's parser does too
/// (`Parser/VCF.pm`, release/115: `if(scalar keys %first_bases == 1)`). An
/// allele that is a clean deletion therefore arrives as a replacement. TTC28
/// `22:28225368 AAAGAAG>AAAG,A` reaches the annotator as `AAGAAG` -> `AAG` and
/// came out `c.553-61775_553-61770delinsCTT`, where VEP writes
/// `c.553-61772_553-61770del`: the same six bases, read as the three-base
/// deletion they are.
///
/// Ensembl closes the gap twice over. `Parser.pm`'s `post_process_vfs` sends
/// any record whose alleles differ in length through `minimise_alleles`, and
/// `InputBuffer::split_variants` then turns it into one VariationFeature per
/// ALT, each trimmed against the reference on its own and annotated separately
/// before being rejoined for output; `_clip_alleles` in
/// `TranscriptVariationAllele.pm` clips again while building the notation.
/// Neither is gated on `--minimal`.
///
/// This clips where the description is built, which is the second of the two.
/// A description carries its own coordinates, so getting it right needs nothing
/// else to move: the consequence caller reads the differing region out of the
/// untrimmed pair and agrees with VEP on every row of all three samples, and
/// splitting the site into one variant per allele to match the first would put
/// that agreement at risk for no field that is wrong today. What stays
/// untrimmed is the *reported allele* - `AAG` here where VEP prints `-` - and
/// the positional fields beside it, which describe the site's span.
///
/// The description was 801 of the 824 HGVSc rows disagreeing with real VEP
/// 115.1 over a 1-in-200 sample of the GIAB HG002 callset.
pub struct HgvsAllele {
    /// Genomic span of the clipped change. An insertion leaves `start` one past
    /// `end`, which is Ensembl's zero-length interval.
    pub start: u64,
    pub end: u64,
    /// The clipped pair as the genome reads it.
    pub genomic_ref: fastvep_core::Allele,
    pub genomic_alt: fastvep_core::Allele,
    /// The same pair as the transcript reads it.
    pub hgvs_ref: fastvep_core::Allele,
    pub hgvs_alt: fastvep_core::Allele,
    /// cDNA span of the clipped change, for the caller that had one before the
    /// clip. Both ends of a cDNA pair are exonic, so the bases clipped off
    /// either end are exonic too and the span moves with them.
    pub cdna: Option<(u64, u64)>,
}

/// How many bases a pair repeats at its front and at its back, clipping the
/// front first or the back first.
///
/// The order decides the answer whenever one side runs out - `AAGAAG` against
/// `AAG` is the front three bases or the back three, never both - and Ensembl
/// takes the front of the sequence *as the transcript reads it*, which is the
/// back of the genomic one on the reverse strand. Where neither side runs out
/// the two orders agree: the leading run and the trailing run cannot overlap
/// while a differing base separates them.
fn shared_ends(r: &[u8], a: &[u8], back_first: bool) -> (usize, usize) {
    let (mut front, mut back) = (0usize, 0usize);
    let clip_front = |front: &mut usize, back: usize| {
        while *front + back < r.len()
            && *front + back < a.len()
            && r[*front].eq_ignore_ascii_case(&a[*front])
        {
            *front += 1;
        }
    };
    let clip_back = |back: &mut usize, front: usize| {
        while front + *back < r.len()
            && front + *back < a.len()
            && r[r.len() - 1 - *back].eq_ignore_ascii_case(&a[a.len() - 1 - *back])
        {
            *back += 1;
        }
    };
    if back_first {
        clip_back(&mut back, front);
        clip_front(&mut front, back);
    } else {
        clip_front(&mut front, back);
        clip_back(&mut back, front);
    }
    (front, back)
}

/// `g.` for one allele, clipped of the bases its reference and alternate repeat
/// at either end.
///
/// Ensembl clips here too, and in genomic orientation: `hgvs_genomic`
/// (`VariationFeature.pm`, release/115) calls `trim_sequences` without a strand,
/// so the front of the *sequence* gives way first whichever way the gene runs.
/// The same allele of TTC28 `22:28225368 AAAGAAG>AAAG,A` that reaches
/// [`hgvs_allele`] as a replacement was written `22:g.28225369_28225374delinsAAG`
/// here, beside an HGVSc that had already been read as the deletion it is.
///
/// This does not close the field: Ensembl also applies the 3'-rule over the
/// genome, which nothing here does, so `9:905488 C>CTGTGTGTG` stays
/// `g.905488_905489insTGTGTGTG` where VEP writes `g.905507_905514dup`. Over a
/// 1-in-200 sample of the GIAB HG002 callset the two disagree on 11 of 13,535
/// records; clipping accounts for part of that and the genomic shift for the
/// rest.
pub fn hgvsg_clipped(
    chrom: &str,
    start: u64,
    end: u64,
    genomic_ref: &fastvep_core::Allele,
    genomic_alt: &fastvep_core::Allele,
) -> String {
    use fastvep_core::Allele;
    let (front, back) = match (genomic_ref, genomic_alt) {
        (Allele::Sequence(r), Allele::Sequence(a)) => shared_ends(r, a, false),
        _ => (0, 0),
    };
    // Nothing repeats, which is every ordinary variant: answer from the pair
    // the caller already holds rather than building a copy of it to answer from.
    if front + back == 0 {
        return fastvep_hgvs::hgvsg(chrom, start, end, genomic_ref, genomic_alt);
    }
    let clipped = |bases: &[u8]| match &bases[front..bases.len() - back] {
        [] => Allele::Deletion,
        kept => Allele::Sequence(kept.to_vec()),
    };
    let (Allele::Sequence(r), Allele::Sequence(a)) = (genomic_ref, genomic_alt) else {
        unreachable!("a non-zero clip needs two sequences")
    };
    fastvep_hgvs::hgvsg(
        chrom,
        start + front as u64,
        end - back as u64,
        &clipped(r),
        &clipped(a),
    )
}

/// The genomic span one allele actually changes, clipped of the bases its
/// reference and alternate repeat at either end.
///
/// [`hgvs_allele`] is the same clip carrying the alleles along; this is the span
/// alone, for a caller that needs to know *where* the change is and not what to
/// call it. It allocates nothing.
pub fn clipped_span(
    strand: fastvep_core::Strand,
    start: u64,
    end: u64,
    genomic_ref: &fastvep_core::Allele,
    genomic_alt: &fastvep_core::Allele,
) -> (u64, u64) {
    use fastvep_core::{Allele, Strand};
    let (front, back) = match (genomic_ref, genomic_alt) {
        (Allele::Sequence(r), Allele::Sequence(a)) => shared_ends(r, a, strand == Strand::Reverse),
        _ => (0, 0),
    };
    (start + front as u64, end - back as u64)
}

/// Read one allele of a site as HGVS describes it. See [`HgvsAllele`].
///
/// `start`/`end` and the pair are genomic; `cdna` is the span the predictor
/// mapped for the unclipped reference, low end first.
pub fn hgvs_allele(
    strand: fastvep_core::Strand,
    start: u64,
    end: u64,
    genomic_ref: &fastvep_core::Allele,
    genomic_alt: &fastvep_core::Allele,
    cdna: Option<(u64, u64)>,
) -> HgvsAllele {
    use fastvep_core::{Allele, Strand};

    let reverse = strand == Strand::Reverse;
    let (front, back) = match (genomic_ref, genomic_alt) {
        // Only a pair of sequences can repeat anything. A `-` on either side is
        // already minimal, and `*` or a symbolic allele names no bases to
        // compare.
        (Allele::Sequence(r), Allele::Sequence(a)) => shared_ends(r, a, reverse),
        _ => (0, 0),
    };

    let clipped = |bases: &[u8]| -> Allele {
        match &bases[front..bases.len() - back] {
            [] => Allele::Deletion,
            kept => Allele::Sequence(kept.to_vec()),
        }
    };
    let (genomic_ref, genomic_alt) = match (genomic_ref, genomic_alt) {
        (Allele::Sequence(r), Allele::Sequence(a)) if front + back > 0 => (clipped(r), clipped(a)),
        (r, a) => (r.clone(), a.clone()),
    };
    let (hgvs_ref, hgvs_alt) = if reverse {
        (
            crate::reverse_complement_allele(&genomic_ref),
            crate::reverse_complement_allele(&genomic_alt),
        )
    } else {
        (genomic_ref.clone(), genomic_alt.clone())
    };

    // cDNA runs in transcript order, so the front of the genomic pair is its
    // low end only on the forward strand.
    let (cdna_front, cdna_back) = if reverse {
        (back, front)
    } else {
        (front, back)
    };
    HgvsAllele {
        start: start + front as u64,
        end: end - back as u64,
        genomic_ref,
        genomic_alt,
        hgvs_ref,
        hgvs_alt,
        // Saturating because a clip is bounded by the shorter allele, not by
        // the transcript: a change at cDNA position 1 whose back clips away can
        // land at 0, which is not a position, and the renderer refuses it.
        cdna: cdna.map(|(lo, hi)| (lo + cdna_front as u64, hi.saturating_sub(cdna_back as u64))),
    }
}

/// Do two genomic positions of one transcript sit in the same exon or intron, or
/// in an exon and an intron that touch?
///
/// That is the widest span a single `c.` range can name without an offset having
/// to run through zero - see [`intronic_ins_as_dup`], which is the only caller.
fn one_region_or_adjacent(transcript: &fastvep_genome::Transcript, first: u64, last: u64) -> bool {
    let (a, b) = (
        transcript.intron_bounds_at(first),
        transcript.intron_bounds_at(last),
    );
    match (a, b) {
        // Same intron, or both exonic.
        (Some(x), Some(y)) => x == y,
        (None, None) => {
            // Both exonic: the same exon, since a block that left one would have
            // an intronic base between.
            let (lo, hi) = (first.min(last), first.max(last));
            transcript
                .exons
                .iter()
                .any(|e| e.start <= lo && hi <= e.end)
        }
        // One of each: the intron has to touch the exon the other end is in.
        (Some((s, e)), None) | (None, Some((s, e))) => {
            let exonic = if a.is_some() { last } else { first };
            transcript.exons.iter().any(|x| {
                x.start <= exonic && exonic <= x.end && (x.end + 1 == s || e + 1 == x.start)
            })
        }
    }
}

/// Rewrite a 3'-shifted intronic insertion as a duplication, when it is one.
///
/// `hgvsc` is the insertion notation already built for the shifted position, and
/// is returned rewritten. `None` means the insertion does not duplicate an
/// adjacent block, or the block leaves the intron it started in - a range
/// written across that boundary names bases in the wrong intron, so the
/// insertion notation it already carries is the correct description.
///
/// `coding_start` is `None` for a transcript numbered from its first base, which
/// selects `n.` numbering.
// The arguments are independent coordinates, alleles and transcript state with
// no natural grouping; a struct would only move the list to the call site.
#[allow(clippy::too_many_arguments)]
pub fn intronic_ins_as_dup(
    seq_provider: &dyn SequenceProvider,
    chrom: &str,
    transcript: &fastvep_genome::Transcript,
    hgvsc: &str,
    shifted_start: u64,
    ins_bases: &[u8],
    shift: u64,
    coding_start: Option<u64>,
    coding_end: Option<u64>,
) -> Option<String> {
    let (lo, hi) = intronic_dup_span(
        seq_provider,
        chrom,
        shifted_start,
        ins_bases,
        shift,
        transcript.strand,
    )?;
    // Transcript order, not genomic: on the reverse strand the block's 5' end is
    // its higher coordinate.
    let (first, last) = match transcript.strand {
        fastvep_core::Strand::Forward => (lo, hi),
        fastvep_core::Strand::Reverse => (hi, lo),
    };
    let start = crate::intronic_or_exonic_cdna(transcript, first)?;
    let end = crate::intronic_or_exonic_cdna(transcript, last)?;
    // The two ends have to be writable as one range, and an offset counts from
    // its own exon and does not run through zero: stepping back from `+1` does
    // not arrive at `-2`, it names bases in the intron before. So a block may
    // span one splice site at most - the exon and the intron on one side of it,
    // where Ensembl writes `c.1176_1183+7dup` - and never two. Writing a
    // crossing range anyway put PVS1's offset gate at `-2` for a duplication
    // sitting on the donor and called two ClinVar-benign MSH6 and DSP variants
    // likely pathogenic.
    //
    // A span running past the intron's own midpoint is fine, and Ensembl writes
    // it from the exon on either side: `c.5044+27_5045-47dup`. Requiring one
    // shared anchor instead of one shared region left 160 such rows as `ins`.
    if !one_region_or_adjacent(transcript, first, last) {
        return None;
    }
    match coding_start {
        Some(cs) => convert_ins_to_dup_range(hgvsc, start, end, cs, coding_end),
        None => convert_ins_to_dup_range_noncoding(hgvsc, start, end),
    }
}

/// Build the HGVSc for an indel: 3'-shifted over the genome, and written as a
/// duplication where the shifted insertion sits against the block it copies.
///
/// Every indel routes here, exonic or not, because the shift does not stop at a
/// splice site and so neither can the description. Each end is mapped on its
/// own and an exonic one is written as its anchor alone, which is how a change
/// that travels out of an exon and into the intron beyond comes out
/// `c.220+1del` and one travelling the other way `c.1793del`.
///
/// Shifting on the *spliced* sequence instead - which is what the exonic path
/// did - is wrong in both directions at once. It cannot follow a repeat into the
/// intron, and it will happily walk a deletion over a splice junction, naming a
/// block that is contiguous in the mRNA and not in the DNA the description is
/// numbered against: `c.71_81del` for BRCA1 `c.70_80del`, where the eleventh
/// base named is the first base of the next exon.
///
/// Both annotation loops call this. They had drifted: the CLI's copy shifted and
/// converted to `dup`, the library's did neither, so the same intronic
/// duplication came out normalised from `fastvep annotate` and unnormalised from
/// the server and the web UI.
///
/// `genomic_ref` and `genomic_alt` are the alleles as the VCF carried them;
/// `hgvs_ref` and `hgvs_alt` are the same pair in transcript orientation. The
/// shift reads the reference, so it needs the genomic pair; the notation is
/// written from the transcript pair.
///
/// `coding_start` is `None` for a transcript numbered from its first base, which
/// selects `n.` numbering.
// The arguments are independent coordinates, alleles and transcript state with
// no natural grouping; a struct would only move the list to the call site.
#[allow(clippy::too_many_arguments)]
pub fn hgvsc_shifted(
    seq_provider: Option<&dyn SequenceProvider>,
    chrom: &str,
    transcript: &fastvep_genome::Transcript,
    versioned_tid: &str,
    var_start: u64,
    var_end: u64,
    genomic_ref: &fastvep_core::Allele,
    genomic_alt: &fastvep_core::Allele,
    hgvs_ref: &fastvep_core::Allele,
    hgvs_alt: &fastvep_core::Allele,
    coding_start: Option<u64>,
    coding_end: Option<u64>,
) -> Option<String> {
    use fastvep_core::{Allele, Strand};

    let is_insertion = matches!(
        (hgvs_ref, hgvs_alt),
        (Allele::Deletion, Allele::Sequence(_))
    );
    let is_indel = is_insertion
        || matches!(
            (hgvs_ref, hgvs_alt),
            (Allele::Sequence(_), Allele::Deletion)
        );

    // A change that does not touch the transcript has no position on it, and the
    // 3'-rule must not give it one: an upstream deletion sitting in a repeat
    // that runs into the first exon would otherwise shift in and be described as
    // `n.1_3del`, naming bases of a transcript it never reaches. VEP writes no
    // HGVSc for those, and neither did this before every indel started routing
    // through here.
    let (var_lo, var_hi) = (var_start.min(var_end), var_start.max(var_end));
    if var_hi < transcript.start || var_lo > transcript.end {
        return None;
    }
    // The walk is bounded by the *transcript*, not by the exon or intron the
    // variant starts in. HGVS shifts over the reference sequence a description
    // is numbered against, and for `c.` that sequence is genomic - so a change
    // whose repeat runs on past the splice site travels with it.
    let (shifted_start, shifted_end) = match seq_provider.filter(|_| is_indel) {
        Some(sp) => three_prime_shift_genomic(
            sp,
            chrom,
            var_start,
            var_end,
            genomic_ref,
            genomic_alt,
            transcript.strand,
            transcript.start,
            transcript.end,
        ),
        None => (var_start, var_end),
    };
    // Distance travelled, which is also how far the inserted string rotated.
    let shift = match transcript.strand {
        Strand::Forward => shifted_start.saturating_sub(var_start),
        Strand::Reverse => var_start.saturating_sub(shifted_start),
    };
    // `hgvs_alt` is already in transcript orientation, where the shift always
    // travels 3' whichever way the transcript runs, so the string always rotates
    // left - the strand is spent before this point, on reverse-complementing it.
    //
    // Rotating right on the reverse strand instead, which is what this did,
    // named the right position with the wrong bases: `c.21-21_21-20insGTC` for a
    // variant the same transcript calls `insCGT` when the VCF spells it one base
    // over. Genomically that rotation is correct and [`intronic_dup_span`] keeps
    // it, because that function reads the reference rather than the transcript.
    let shifted_alt = match hgvs_alt {
        Allele::Sequence(ins) if is_insertion && shift > 0 && !ins.is_empty() => {
            let mut rotated = ins.clone();
            let k = (shift % rotated.len() as u64) as usize;
            rotated.rotate_left(k);
            Allele::Sequence(rotated)
        }
        other => other.clone(),
    };

    // An insertion is written over the two bases it sits between, and both are
    // mapped: `shifted_end` and `shifted_start` are the pair, in transcript
    // order on the forward strand and reversed on the reverse.
    //
    // Letting the renderer infer the second coordinate as `offset + 1` instead
    // breaks across the middle of an intron, where `+n` counts from one exon and
    // `-m` from the next: an insertion between `c.20+30` and `c.21-30` came out
    // `c.20+30_20+31ins…`, naming a base past the half the donor-side offsets
    // reach.
    let (anchor_pos, second_pos) = if is_insertion {
        match transcript.strand {
            Strand::Forward => (shifted_end, Some(shifted_end + 1)),
            Strand::Reverse => (shifted_end + 1, Some(shifted_end)),
        }
    } else {
        (
            shifted_start,
            (shifted_start != shifted_end).then_some(shifted_end),
        )
    };
    let (cdna_pos, offset) = crate::intronic_or_exonic_cdna(transcript, anchor_pos)?;
    let (end_cdna, end_offset) = second_pos
        .and_then(|p| crate::intronic_or_exonic_cdna(transcript, p))
        .map(|(c, o)| (Some(c), Some(o)))
        .unwrap_or((None, None));

    let hgvsc = match coding_start {
        Some(cs) => fastvep_hgvs::hgvsc_intronic_range(
            versioned_tid,
            cdna_pos,
            offset,
            end_cdna,
            end_offset,
            hgvs_ref,
            &shifted_alt,
            cs,
            coding_end,
        ),
        None => fastvep_hgvs::hgvsc_noncoding_intronic_range(
            versioned_tid,
            cdna_pos,
            offset,
            end_cdna,
            end_offset,
            hgvs_ref,
            &shifted_alt,
        ),
    }?;

    if let (true, Allele::Sequence(ins), Some(sp)) = (is_insertion, genomic_alt, seq_provider) {
        if hgvsc.contains("ins") && !ins.is_empty() {
            if let Some(dup) = intronic_ins_as_dup(
                sp,
                chrom,
                transcript,
                &hgvsc,
                shifted_start,
                ins,
                shift,
                coding_start,
                coding_end,
            ) {
                return Some(dup);
            }
        }
    }
    Some(hgvsc)
}

#[cfg(test)]
mod tests {
    use super::*;
    use anyhow::{anyhow, Result};
    use fastvep_core::{Allele, Strand};
    use fastvep_genome::{Exon, Gene, Transcript};

    /// Minimal `SequenceProvider` over a 1-based reference string for one contig,
    /// mirroring the real readers' contract: 1-based inclusive, `Err` past the end.
    struct StrRef(&'static str);
    impl SequenceProvider for StrRef {
        fn fetch_sequence(&self, _chrom: &str, start: u64, end: u64) -> Result<Vec<u8>> {
            if start < 1 || end < start {
                return Err(anyhow!("bad range"));
            }
            let b = self.0.as_bytes();
            let s0 = (start - 1) as usize;
            if s0 >= b.len() {
                return Err(anyhow!("past contig end"));
            }
            Ok(b[s0..(end as usize).min(b.len())].to_vec())
        }
    }

    /// Two exons on `strand` with one intron between them, so an intronic
    /// position has an anchor on either side. Exon 1 is 1..=20, exon 2 is
    /// 81..=100, and the intron is 21..=80.
    fn transcript(strand: Strand) -> Transcript {
        let exon = |start: u64, end: u64, rank: u32| Exon {
            stable_id: format!("ENSE{}", rank),
            start,
            end,
            strand,
            phase: 0,
            end_phase: 0,
            rank,
        };
        Transcript {
            stable_id: "ENST00000000001".into(),
            version: Some(1),
            gene: Gene {
                stable_id: "ENSG00000000001".into(),
                symbol: Some("TEST".into()),
                symbol_source: None,
                hgnc_id: None,
                biotype: "protein_coding".into(),
                chromosome: "1".into(),
                start: 1,
                end: 100,
                strand,
            },
            biotype: "protein_coding".into(),
            chromosome: "1".into(),
            start: 1,
            end: 100,
            strand,
            exons: vec![exon(1, 20, 1), exon(81, 100, 2)],
            translation: None,
            cdna_coding_start: Some(1),
            cdna_coding_end: Some(40),
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
        }
    }

    /// 20 exonic bases, then a `TG` repeat filling the intron, then exon 2. An
    /// insertion of `TG` anywhere in that repeat is the same variant.
    const TG_REPEAT: &str = "AAAAAAAAAAAAAAAAAAAA\
                             TGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTG\
                             CCCCCCCCCCCCCCCCCCCC";

    /// The duplicated block sits immediately 5' of the *shifted* insertion
    /// point, so a `TG` insertion that travelled the whole repeat names the last
    /// two bases of it, not the two it was written against.
    #[test]
    fn dup_span_follows_the_shifted_insertion_not_the_vcf_position() {
        let r = StrRef(TG_REPEAT);
        // Insertion point right after exon 1; the repeat runs 21..=80, so the
        // maximal 3' shift on the forward strand travels its full 60 bases.
        let span = intronic_dup_span(&r, "1", 81, b"TG", 60, Strand::Forward);
        assert_eq!(span, Some((79, 80)));
    }

    /// The inserted string rotates one base per position travelled. After an odd
    /// shift of a two-base insert the block is `GT`, and reading the unrotated
    /// `TG` against the reference would reject a duplication that is real.
    #[test]
    fn dup_span_compares_against_the_rotated_insert() {
        let r = StrRef(TG_REPEAT);
        // An odd shift lands the insertion point between a G and a T.
        let span = intronic_dup_span(&r, "1", 80, b"TG", 59, Strand::Forward);
        assert_eq!(span, Some((78, 79)));
        let block = r.fetch_sequence("1", 78, 79).unwrap();
        assert_eq!(&block, b"GT", "the block is the rotated form, not `TG`");
    }

    /// 3' on the reverse strand runs towards lower coordinates, so the block a
    /// reverse-strand insertion duplicates lies *after* the insertion point.
    #[test]
    fn dup_span_reads_the_other_side_on_the_reverse_strand() {
        let r = StrRef(TG_REPEAT);
        let span = intronic_dup_span(&r, "1", 21, b"TG", 0, Strand::Reverse);
        assert_eq!(span, Some((21, 22)));
    }

    /// A one-base walk through a homopolymer is not the repeat test: this is the
    /// shape the old dup anchor used, and it stopped after one base of `TGTG…`.
    #[test]
    fn dup_span_is_none_when_the_block_does_not_repeat() {
        let r = StrRef(TG_REPEAT);
        // Inside the exon's poly-A run, a `TG` insert duplicates nothing.
        assert_eq!(
            intronic_dup_span(&r, "1", 15, b"TG", 0, Strand::Forward),
            None
        );
    }

    #[test]
    fn dup_span_declines_rather_than_reading_past_the_contig() {
        let r = StrRef(TG_REPEAT);
        assert_eq!(
            intronic_dup_span(&r, "1", 1, b"TG", 0, Strand::Forward),
            None
        );
        assert_eq!(
            intronic_dup_span(&r, "1", 0, b"TG", 0, Strand::Forward),
            None
        );
        assert_eq!(
            intronic_dup_span(&r, "1", 81, b"", 0, Strand::Forward),
            None
        );
    }

    /// The whole rewrite: a shifted intronic insertion becomes a `dup` naming the
    /// block at its 3' end, written from the anchor its own side of the intron.
    #[test]
    fn intronic_insertion_becomes_a_dup_at_its_shifted_position() {
        let r = StrRef(TG_REPEAT);
        let tr = transcript(Strand::Forward);
        let out = intronic_ins_as_dup(
            &r,
            "1",
            &tr,
            "ENST00000000001.1:c.20+59_20+60insTG",
            81,
            b"TG",
            60,
            Some(1),
            Some(40),
        );
        // 79 and 80 are the last two intronic bases, one and two before exon 2,
        // so they are written from the *downstream* anchor: c.21-2_21-1.
        assert_eq!(out.as_deref(), Some("ENST00000000001.1:c.21-2_21-1dup"));
    }

    /// A duplicated block that sits inside the exon is named there, in exonic
    /// coordinates. It reads as an intronic case only from the insertion
    /// notation the caller brings; the block itself is two exonic bases, and
    /// refusing it - which is what mapping both ends through
    /// `genomic_to_intronic_cdna` did, since that returns `None` for an exonic
    /// position - left a duplication written as a plain insertion.
    #[test]
    fn a_block_inside_the_exon_is_named_in_exonic_coordinates() {
        let r = StrRef(TG_REPEAT);
        let tr = transcript(Strand::Reverse);
        // Genomic 81-82 are the first two bases of exon 2, which on the reverse
        // strand are c.20 and c.19.
        let out = intronic_ins_as_dup(
            &r,
            "1",
            &tr,
            "ENST00000000001.1:c.21-1_21insCC",
            81,
            b"CC",
            0,
            Some(1),
            Some(40),
        );
        assert_eq!(out.as_deref(), Some("ENST00000000001.1:c.19_20dup"));
    }

    /// One `c.` range can name a span that crosses one splice site and no more.
    /// An offset counts from its own exon and does not run through zero, so a
    /// block reaching from one exon across an intron into the next cannot be
    /// written as a range at all: stepping back from `+1` does not arrive at
    /// `-2`, it names bases in the intron before. Writing the crossing range
    /// anyway put PVS1's offset gate at `-2` for a duplication sitting on the
    /// donor and called two ClinVar-benign MSH6 and DSP variants likely
    /// pathogenic.
    #[test]
    fn a_range_may_cross_one_splice_site_and_no_more() {
        let tr = transcript(Strand::Forward);
        // Exon 1 is 1..=20, the intron 21..=80, exon 2 is 81..=100.
        assert!(one_region_or_adjacent(&tr, 5, 10), "one exon");
        assert!(one_region_or_adjacent(&tr, 30, 40), "one intron");
        assert!(
            one_region_or_adjacent(&tr, 18, 25),
            "exon into the intron after it"
        );
        assert!(
            one_region_or_adjacent(&tr, 75, 85),
            "intron into the exon after it"
        );
        assert!(
            !one_region_or_adjacent(&tr, 10, 85),
            "across a whole intron"
        );
        assert!(one_region_or_adjacent(&tr, 5, 19), "the whole of one exon");
    }

    /// The property the whole shift exists to provide, and the one that broke:
    /// every spelling of the same insertion inside a repeat must name the same
    /// block. Walking every position of the `TG` repeat, each with the rotation
    /// that spelling carries, must land on one answer.
    ///
    /// Against real Ensembl VEP 115.1 this was 0 of 168 transcript rows agreeing
    /// across the two spellings of six genome-wide variants; VEP agreed on all
    /// 168.
    #[test]
    fn every_spelling_of_one_insertion_names_the_same_block() {
        let r = StrRef(TG_REPEAT);
        // The repeat runs 21..=80. An insertion written at offset `k` into it
        // carries the insert rotated by `k`, and has `60 - k` left to travel.
        let answers: std::collections::HashSet<_> = (0..60)
            .map(|k| {
                let mut ins = b"TG".to_vec();
                ins.rotate_left(k % 2);
                intronic_dup_span(&r, "1", 81, &ins, (60 - k) as u64, Strand::Forward)
            })
            .collect();
        assert_eq!(
            answers,
            [Some((79, 80))].into_iter().collect(),
            "every spelling must name the last two bases of the repeat"
        );
    }

    /// A span running past the intron's own midpoint is still one intron, and is
    /// written from the exon on either side. Requiring one shared anchor instead
    /// left 160 such rows as `ins` where Ensembl writes `c.5044+27_5045-47dup`.
    #[test]
    fn a_span_crossing_the_intron_midpoint_is_still_one_dup() {
        // A 60-base insert filling the whole intron duplicates all of it.
        let r = StrRef(TG_REPEAT);
        let tr = transcript(Strand::Forward);
        let ins: Vec<u8> = TG_REPEAT.as_bytes()[20..80].to_vec();
        let out = intronic_ins_as_dup(
            &r,
            "1",
            &tr,
            "ENST00000000001.1:c.20+60_21-0ins…",
            81,
            &ins,
            0,
            Some(1),
            Some(40),
        );
        assert_eq!(out.as_deref(), Some("ENST00000000001.1:c.20+1_21-1dup"));
    }

    /// The offset a criterion reads is measured on the transcript, over the span
    /// the allele actually changes. The fixture's intron runs 21..80, so the
    /// donor's own base is 21.
    #[test]
    fn the_offset_is_measured_over_the_change_and_not_over_the_record() {
        let tr = transcript(Strand::Forward);
        let offset = |r: &str, a: &str, start: u64, end: u64| {
            let (lo, hi) = clipped_span(
                Strand::Forward,
                start,
                end,
                &Allele::Sequence(r.as_bytes().to_vec()),
                &Allele::Sequence(a.as_bytes().to_vec()),
            );
            tr.intronic_offset_covered(lo, hi)
        };
        // Six bases from the donor's own base, of which the first three are
        // unchanged: the record reaches `+1`, the deletion sits at `+4`.
        assert_eq!(offset("AAGAAG", "AAG", 21, 26), Some(4));
        // The same record with nothing to clip is the `+1` it looks like.
        assert_eq!(offset("AAGAAG", "C", 21, 26), Some(1));
        // Wholly exonic reaches no intronic base at all.
        assert_eq!(offset("AA", "C", 10, 11), None);
        // Running out of the exon into the intron reaches the first base of it,
        // whatever the far end reads - the same rule the HGVS string parser
        // follows, because the two answers are compared against each other.
        assert_eq!(offset("AAAAA", "C", 19, 23), Some(1));
        // Deep in the intron, the nearer boundary wins.
        assert_eq!(offset("AA", "C", 76, 77), Some(-4));
    }

    /// Which boundary an intronic base counts from is a property of the
    /// transcript's direction, so the two strands mirror each other.
    #[test]
    fn the_strand_decides_which_boundary_the_offset_counts_from() {
        let (fwd, rev) = (transcript(Strand::Forward), transcript(Strand::Reverse));
        // 21 is the base after exon 1 and 80 the base before exon 2.
        assert_eq!(fwd.intronic_offset_covered(21, 21), Some(1));
        assert_eq!(fwd.intronic_offset_covered(80, 80), Some(-1));
        assert_eq!(rev.intronic_offset_covered(80, 80), Some(1));
        assert_eq!(rev.intronic_offset_covered(21, 21), Some(-1));
        // A span reaching both exons covers the whole intron and names no
        // boundary, which is what the HGVS string for it says too.
        assert_eq!(fwd.intronic_offset_covered(10, 90), None);
    }

    /// The rule is "clip what the two repeat at either end", not "strip the one
    /// base the VCF anchored on": what is left of `AAGAAG` -> `AAG` is a
    /// deletion, and naming it a replacement invents three bases of change.
    #[test]
    fn a_pair_that_repeats_at_one_end_clips_down_to_the_change_it_is() {
        let hv = hgvs_allele(
            Strand::Forward,
            100,
            105,
            &Allele::Sequence(b"AAGAAG".to_vec()),
            &Allele::Sequence(b"AAG".to_vec()),
            Some((10, 15)),
        );
        assert_eq!((hv.start, hv.end), (103, 105));
        assert_eq!(hv.genomic_ref, Allele::Sequence(b"AAG".to_vec()));
        assert_eq!(hv.genomic_alt, Allele::Deletion);
        assert_eq!(hv.cdna, Some((13, 15)));
    }

    /// Which end is clipped when only one can be decides where the change
    /// lands, and Ensembl clips the front of the sequence *as the transcript
    /// reads it* - the far end of the genomic one on the reverse strand.
    #[test]
    fn the_strand_decides_which_repeated_end_gives_way() {
        let call = |strand| {
            let hv = hgvs_allele(
                strand,
                100,
                102,
                &Allele::Sequence(b"ACA".to_vec()),
                &Allele::Sequence(b"A".to_vec()),
                Some((10, 12)),
            );
            (hv.start, hv.end, hv.cdna)
        };
        // Forward: the front two bases match, so the back two go.
        assert_eq!(call(Strand::Forward), (101, 102, Some((11, 12))));
        // Reverse: `TGT` against `T` matches at its front, which is the genomic
        // back, so the front two go instead - and cDNA, which runs the other
        // way, still loses its first two.
        assert_eq!(call(Strand::Reverse), (100, 101, Some((11, 12))));
    }

    /// A clip that empties the reference leaves Ensembl's zero-length interval,
    /// `start` one past `end`, which is what every insertion downstream reads.
    #[test]
    fn a_pair_whose_reference_clips_away_becomes_an_insertion() {
        let hv = hgvs_allele(
            Strand::Forward,
            100,
            100,
            &Allele::Sequence(b"T".to_vec()),
            &Allele::Sequence(b"TT".to_vec()),
            Some((10, 10)),
        );
        assert_eq!((hv.start, hv.end), (101, 100));
        assert_eq!(hv.genomic_ref, Allele::Deletion);
        assert_eq!(hv.genomic_alt, Allele::Sequence(b"T".to_vec()));
        assert!(is_shiftable_indel(&hv.hgvs_ref, &hv.hgvs_alt));
        assert_eq!(hv.cdna, Some((11, 10)));
    }

    /// Both ends may repeat at once. An equal-length pair clips to the same core
    /// whichever end goes first, so the strand cannot move it.
    #[test]
    fn a_pair_repeating_at_both_ends_clips_to_the_same_core_either_way() {
        for strand in [Strand::Forward, Strand::Reverse] {
            let hv = hgvs_allele(
                strand,
                100,
                103,
                &Allele::Sequence(b"ACGT".to_vec()),
                &Allele::Sequence(b"ATGT".to_vec()),
                Some((10, 13)),
            );
            assert_eq!((hv.start, hv.end), (101, 101));
            assert_eq!(hv.genomic_ref, Allele::Sequence(b"C".to_vec()));
            assert_eq!(hv.genomic_alt, Allele::Sequence(b"T".to_vec()));
        }
    }

    /// Nothing repeats, nothing moves. This is every ordinary variant.
    #[test]
    fn a_pair_that_repeats_nothing_is_left_where_it_was() {
        for (r, a) in [
            (
                Allele::Sequence(b"A".to_vec()),
                Allele::Sequence(b"G".to_vec()),
            ),
            (
                Allele::Sequence(b"AC".to_vec()),
                Allele::Sequence(b"GT".to_vec()),
            ),
            (Allele::Sequence(b"AC".to_vec()), Allele::Deletion),
            (Allele::Deletion, Allele::Sequence(b"AC".to_vec())),
            (Allele::Sequence(b"A".to_vec()), Allele::Missing),
        ] {
            let hv = hgvs_allele(Strand::Forward, 100, 101, &r, &a, Some((10, 11)));
            assert_eq!((hv.start, hv.end, hv.cdna), (100, 101, Some((10, 11))));
            assert_eq!((hv.genomic_ref, hv.genomic_alt), (r, a));
        }
    }

    /// The reverse-strand pair is complemented as well as clipped, and the two
    /// orientations have to agree about which bases are left.
    #[test]
    fn the_transcript_sees_the_clipped_pair_complemented() {
        let hv = hgvs_allele(
            Strand::Reverse,
            100,
            105,
            &Allele::Sequence(b"AAGAAG".to_vec()),
            &Allele::Sequence(b"AAG".to_vec()),
            None,
        );
        // `CTTCTT` against `CTT` repeats at its front, and the transcript's
        // front is the genome's back, so the last three bases go.
        assert_eq!((hv.start, hv.end), (100, 102));
        assert_eq!(hv.genomic_ref, Allele::Sequence(b"AAG".to_vec()));
        assert_eq!(hv.hgvs_ref, Allele::Sequence(b"CTT".to_vec()));
        assert_eq!(hv.hgvs_alt, Allele::Deletion);
    }
}
