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
    /// Where `start` would be had the front of the pair been clipped first, which
    /// is where the VCF anchors the change (the base every allele shares at the
    /// front is the anchor) and what `HGVS_OFFSET` is measured from.
    ///
    /// It is `start` except on the reverse strand, where the pair is clipped from
    /// the transcript's end of it, and then only where the two clips overlap: `CTGCAC`
    /// to `C` deletes `TGCAC` from the front and `CTGCA` from the back, one base
    /// apart. A site that carries a single ALT is clipped of its anchor by the VCF
    /// parser, so for it the two are the same; a site with several ALTs is not, and
    /// without this the same deletion would be reported one base short of the shift
    /// it was given alone.
    pub anchor_start: u64,
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
        // Saturating, as in `hgvs_allele`: `end` need not belong to these alleles.
        end.saturating_sub(back as u64),
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
    // Saturating, as in `hgvs_allele`: `end` need not belong to these alleles.
    (start + front as u64, end.saturating_sub(back as u64))
}

/// Where HGVS places one allele against a transcript: the most 3' position its
/// sequence context allows, read for *whether it touches an exon*.
///
/// [`Transcript::intronic_offset_covered`] over the record's own span answers
/// "how far into the intron is this change written". PVS1, PM4 and BP4 need the
/// other question, and it has a different answer for a change in a repeat. All
/// the positions such a change can be written at produce the same edited
/// sequence, so one that *can* be written wholly inside an intron, clear of
/// `+1` and `+2`, alters no exon base and no canonical dinucleotide when it is
/// written anywhere else. PTEN `10:87931087 GAGGT>G` is `splice_donor_variant`
/// and `c.253+4_253+7del` in both tools, and the intron still opens `GTATGA`
/// afterwards. MSH6 `2:47806358 G>GGTAT` is `frameshift_variant` because the
/// record anchors on the exon's last base, and `c.3801+3_3801+6dup` once HGVS
/// has placed it.
///
/// On a donor the 3' end of the range is the one furthest into the intron, so
/// the placement answers it. **On an acceptor the 3' end is the one nearest the
/// exon**, because the shift runs toward it: a change written deeper in the
/// intron can arrive wholly inside the exon, as MLH1 `TAGG>T` does at an exon
/// that opens `GAGC` (`c.1732_1734del`) and FOXP2 `CCAG>C` does down a `CAG`
/// repeat that runs across the acceptor (`c.618_620del`), and that is
/// [`HgvsPlacement::Exonic`]. That asymmetry is HGVS's, not this function's.
///
/// The allele is clipped of the bases its reference and alternate repeat before
/// anything else, exactly as the HGVSc is, so a multi-allelic record whose
/// alleles arrive as replacements is placed the way it is described. The walk
/// is bounded by the transcript, as the description's is, and a change that does
/// not touch the transcript is not given a position on it.
///
/// The reference is read only to take the shift. Without one - or for a change
/// that is not a pure insertion or deletion - the change is placed where it sits.
#[allow(clippy::too_many_arguments)]
pub fn hgvs_placement(
    seq_provider: Option<&dyn SequenceProvider>,
    chrom: &str,
    transcript: &fastvep_genome::Transcript,
    var_start: u64,
    var_end: u64,
    genomic_ref: &fastvep_core::Allele,
    genomic_alt: &fastvep_core::Allele,
) -> fastvep_core::HgvsPlacement {
    use fastvep_core::HgvsPlacement;
    use fastvep_genome::SpanSite;

    let hv = hgvs_allele(
        transcript.strand,
        var_start,
        var_end,
        genomic_ref,
        genomic_alt,
        None,
    );
    // A span that ends before it starts by more than the one base a pure
    // insertion's does belongs to no change on this transcript: the record's
    // `end` is another allele's. There is no placement to work out, so the
    // record's terms decide, as they do where the placement is not worked out.
    if hv.start > hv.end.saturating_add(1) {
        return HgvsPlacement::Spans;
    }
    let (lo, hi) = (hv.start.min(hv.end), hv.start.max(hv.end));
    if hi < transcript.start || lo > transcript.end {
        return HgvsPlacement::Spans;
    }
    let (start, end) = match seq_provider.filter(|_| is_shiftable_indel(&hv.hgvs_ref, &hv.hgvs_alt))
    {
        Some(sp) => three_prime_shift_genomic(
            sp,
            chrom,
            hv.start,
            hv.end,
            &hv.genomic_ref,
            &hv.genomic_alt,
            transcript.strand,
            transcript.start,
            transcript.end,
        ),
        None => (hv.start, hv.end),
    };
    match transcript.span_site(start, end) {
        SpanSite::Intron { offset } => HgvsPlacement::Intronic { offset },
        SpanSite::Exon => HgvsPlacement::Exonic {
            in_frame: length_change(&hv.genomic_ref, &hv.genomic_alt).is_some_and(|n| n % 3 == 0),
        },
        SpanSite::Spans => HgvsPlacement::Spans,
    }
}

/// The net change in length of a pair of alleles, or `None` where a symbolic or
/// missing allele names no bases to count.
fn length_change(
    genomic_ref: &fastvep_core::Allele,
    genomic_alt: &fastvep_core::Allele,
) -> Option<i64> {
    use fastvep_core::Allele;
    let len = |a: &Allele| match a {
        Allele::Sequence(b) => Some(b.len() as i64),
        Allele::Deletion => Some(0),
        Allele::Missing | Allele::Symbolic(_) => None,
    };
    Some(len(genomic_alt)? - len(genomic_ref)?)
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

    let anchor_front = match (genomic_ref, genomic_alt) {
        (Allele::Sequence(r), Allele::Sequence(a)) => shared_ends(r, a, false).0,
        _ => 0,
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
        anchor_start: start + anchor_front as u64,
        // Saturating because `end` need not belong to these alleles: a sequence
        // ALT beside a symbolic one in the same record takes its `end` from INFO
        // END, which describes the symbolic allele and can be smaller than the
        // clip. A subtraction that wraps aborts a debug or test build.
        end: end.saturating_sub(back as u64),
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
    hgvsc_shifted_with_shift(
        seq_provider,
        chrom,
        transcript,
        versioned_tid,
        var_start,
        var_end,
        var_start,
        genomic_ref,
        genomic_alt,
        hgvs_ref,
        hgvs_alt,
        coding_start,
        coding_end,
    )
    .map(|(hgvsc, _)| hgvsc)
}

/// [`hgvsc_shifted`], and how many bases the change travelled 3' to be written
/// where HGVS asks for it. That distance is VEP's `HGVS_OFFSET` (times the
/// strand): it is what a reader needs to put the description back where the
/// caller's VCF had the variant.
///
/// `anchor_start` is where that is: [`HgvsAllele::anchor_start`], which is
/// `var_start` unless a clip from the transcript's end of the pair left the change
/// a base nearer to where it travels to than the VCF had it.
// The arguments of `hgvsc_shifted` and the anchor, for the same reason.
#[allow(clippy::too_many_arguments)]
pub fn hgvsc_shifted_with_shift(
    seq_provider: Option<&dyn SequenceProvider>,
    chrom: &str,
    transcript: &fastvep_genome::Transcript,
    versioned_tid: &str,
    var_start: u64,
    var_end: u64,
    anchor_start: u64,
    genomic_ref: &fastvep_core::Allele,
    genomic_alt: &fastvep_core::Allele,
    hgvs_ref: &fastvep_core::Allele,
    hgvs_alt: &fastvep_core::Allele,
    coding_start: Option<u64>,
    coding_end: Option<u64>,
) -> Option<(String, u64)> {
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
    // What is reported counts from where the VCF anchors the change, which on the
    // reverse strand can be a base short of where the clip left it; the rotation
    // above must not, because `hgvs_alt` is the string as the clip left it.
    let reported = match transcript.strand {
        Strand::Forward => shifted_start.saturating_sub(anchor_start),
        Strand::Reverse => anchor_start.saturating_sub(shifted_start),
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
                return Some((dup, reported));
            }
        }
    }
    Some((hgvsc, reported))
}

// ---------------------------------------------------------------------------
// HGVS g., as Ensembl writes it
// ---------------------------------------------------------------------------
//
// `VariationFeature::hgvs_genomic` (ensembl-variation, release/116), transcribed
// and then held to a run of VEP 116.3 on every allele of a generated world (see
// `crates/fastvep-cli/fixtures/vep_tab_world`). What it adds to the plain `g.`
// description is the 3'-rule over the genome: an insertion or deletion is written at
// its most 3' position, and an insertion that repeats the bases beside it is a `dup`.

/// `trim_sequences($ref, $alt, $start, $end)` of Ensembl's `Utils::Sequence`, with
/// its defaults: the left first, then the right, and an allele left empty stays empty.
fn trim_pair(reference: &str, alternate: &str, start: i64, end: i64) -> (String, String, i64, i64) {
    let (mut r, mut a) = (reference.as_bytes(), alternate.as_bytes());
    let (mut start, mut end) = (start, end);
    while let (Some((x, r_rest)), Some((y, a_rest))) = (r.split_first(), a.split_first()) {
        if x != y {
            break;
        }
        (r, a) = (r_rest, a_rest);
        start += 1;
    }
    while let (Some((x, r_rest)), Some((y, a_rest))) = (r.split_last(), a.split_last()) {
        if x != y {
            break;
        }
        (r, a) = (r_rest, a_rest);
        end -= 1;
    }
    (
        String::from_utf8_lossy(r).into_owned(),
        String::from_utf8_lossy(a).into_owned(),
        start,
        end,
    )
}

/// Perl's `substr($s, $offset, $length)` for the cases the HGVS code reaches: a
/// negative offset counts from the end, and whatever lies outside the string is
/// simply not there.
fn perl_substr(s: &[u8], offset: i64, length: Option<i64>) -> &[u8] {
    let total = s.len() as i64;
    let first = if offset < 0 { total + offset } else { offset };
    let last = match length {
        Some(l) => first + l,
        None => total,
    };
    let (first, last) = (first.clamp(0, total), last.clamp(0, total));
    if first >= last {
        &[]
    } else {
        &s[first as usize..last as usize]
    }
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum NotationKind {
    Sub,
    Del,
    Ins,
    Dup,
    Inv,
    Delins,
    /// `[n]`, an allele that is a whole multiple of the reference
    Repeat(usize),
    /// ref and alt equal after clipping
    Same,
}

struct Notation {
    start: i64,
    end: i64,
    reference: String,
    alt: String,
    kind: NotationKind,
}

fn reverse_complement_bases(seq: &str) -> String {
    seq.bytes()
        .rev()
        .map(|b| match b {
            b'A' => 'T',
            b'C' => 'G',
            b'G' => 'C',
            b'T' => 'A',
            other => other as char,
        })
        .collect()
}

/// `hgvs_variant_notation` of Ensembl's `Utils::Sequence`: what kind of change the
/// alt is against the reference bases of `ref_start..=ref_end` (1-based, in
/// `ref_seq`), displayed at `display_start..=display_end`. `None` when they are the
/// same bases.
fn hgvs_variant_notation(
    alt_allele: &str,
    ref_seq: &[u8],
    ref_start: i64,
    ref_end: i64,
    display_start: i64,
    display_end: i64,
) -> Option<Notation> {
    let ref_length = (ref_end - ref_start + 1).max(0);
    let alt: String = alt_allele.chars().filter(|&c| c != '-').collect();
    let alt_length = alt.len() as i64;
    let ref_allele =
        String::from_utf8_lossy(perl_substr(ref_seq, ref_start - 1, Some(ref_length))).into_owned();
    if ref_allele == alt {
        return None;
    }
    let mut n = Notation {
        start: display_start,
        end: display_end,
        reference: ref_allele.clone(),
        alt: alt.clone(),
        kind: NotationKind::Delins,
    };
    if alt_length == 0 {
        n.kind = NotationKind::Del;
        return Some(n);
    }
    if ref_length == alt_length {
        n.kind = if ref_length == 1 {
            NotationKind::Sub
        } else if alt == reverse_complement_bases(&ref_allele) {
            NotationKind::Inv
        } else {
            NotationKind::Delins
        };
        return Some(n);
    }
    if ref_length == 0 {
        // An insertion that repeats the bases beside it is a duplication: those after
        // it are looked at first, then those before it. (VEP reverses that for a feature
        // on the other strand, which a VCF record never is, and for an insertion that
        // has been shifted 3' only the bases before it can match anyway.)
        for direction in [1i64, -1] {
            let neighbour = if direction == 1 {
                perl_substr(ref_seq, ref_end, Some(alt_length))
            } else {
                perl_substr(ref_seq, ref_end - alt_length, Some(alt_length))
            };
            if neighbour == alt.as_bytes() {
                if direction == 1 {
                    n.end = display_start + alt_length - 1;
                } else {
                    n.start = display_end - alt_length + 1;
                }
                n.kind = NotationKind::Dup;
                return Some(n);
            }
        }
        n.start = display_end;
        n.end = display_start;
        n.kind = NotationKind::Ins;
        return Some(n);
    }
    if alt_length % ref_length == 0 {
        let multiple = (alt_length / ref_length) as usize;
        if alt == ref_allele.repeat(multiple) {
            n.kind = if multiple == 2 {
                NotationKind::Dup
            } else {
                NotationKind::Repeat(multiple)
            };
            return Some(n);
        }
    }
    Some(n)
}

/// `_clip_alleles` of `VariationFeature.pm`: a delins whose alleles share bases at
/// either end is clipped of them, and may turn out to be something else.
fn clip_alleles(mut n: Notation) -> Notation {
    let mut preseq = String::new();
    while let (Some(r), Some(a)) = (n.reference.chars().next(), n.alt.chars().next()) {
        if r != a {
            break;
        }
        preseq.push(r);
        n.reference.remove(0);
        n.alt.remove(0);
        n.start += 1;
    }
    while let (Some(r), Some(a)) = (n.reference.chars().last(), n.alt.chars().last()) {
        if r != a {
            break;
        }
        n.reference.pop();
        n.alt.pop();
        n.end -= 1;
    }
    if n.reference.is_empty() {
        n.reference = "-".into();
    }
    if n.alt.is_empty() {
        n.alt = "-".into();
    }
    if n.reference == n.alt {
        n.kind = NotationKind::Same;
    } else if n.reference != "-" && n.reference.len() == 1 && n.alt.len() == 1 {
        n.kind = NotationKind::Sub;
    } else if n.reference == "-" && !n.alt.is_empty() {
        let keep = preseq.len().saturating_sub(n.alt.len());
        if n.alt == preseq[keep..] {
            n.kind = NotationKind::Dup;
            n.start -= n.alt.len() as i64;
        } else {
            n.kind = NotationKind::Ins;
            n.start -= 1;
            n.end = n.start + 1;
        }
    } else if !n.reference.is_empty() && n.alt == "-" {
        n.kind = NotationKind::Del;
    }
    n
}

/// `format_hgvs_string` for a genomic description.
fn format_hgvs_g(chrom: &str, n: &Notation) -> Option<String> {
    let numbering = if chrom == "MT" { 'm' } else { 'g' };
    let coordinates = if n.start == n.end {
        n.start.to_string()
    } else {
        format!("{}_{}", n.start, n.end)
    };
    let body = match n.kind {
        NotationKind::Sub => format!("{}{}>{}", n.start, n.reference, n.alt),
        NotationKind::Inv if n.reference.len() == 1 => {
            format!("{}{}>{}", n.start, n.reference, n.alt)
        }
        NotationKind::Del => format!("{coordinates}del"),
        NotationKind::Inv => format!("{coordinates}inv"),
        NotationKind::Dup => format!("{coordinates}dup"),
        NotationKind::Delins => format!("{coordinates}delins{}", n.alt),
        NotationKind::Ins => format!("{coordinates}ins{}", n.alt),
        NotationKind::Repeat(multiple) => format!("{coordinates}[{multiple}]"),
        // VEP warns and writes nothing for a type it cannot name
        NotationKind::Same => return None,
    };
    Some(format!("{chrom}:{numbering}.{body}"))
}

/// The HGVS `g.` description of one allele of a record, as Ensembl VEP writes it.
///
/// `start`, `end`, `reference` and `alts` are the record as the VCF parser leaves
/// it. A lone indel alt is minimised first, as VEP does before it annotates one
/// (`Parser::minimise_alleles`, the same step the known-variant lookup takes).
///
/// The description is read off the reference sequence, not off the VCF's REF: where
/// they differ Ensembl describes the genome, and when the alt is the genome's own
/// base it writes nothing: `None`. Where the reference cannot be read (no FASTA, a
/// contig it lacks) the plain description of [`hgvsg_clipped`] stands in, which is
/// what every run gave before.
pub fn hgvsg_vep(
    seq_provider: &dyn SequenceProvider,
    chrom: &str,
    start: u64,
    end: u64,
    reference: &fastvep_core::Allele,
    alts: &[fastvep_core::Allele],
    alt_index: usize,
) -> Option<String> {
    use fastvep_core::Allele;
    let text = |a: &Allele| -> Option<String> {
        match a {
            Allele::Deletion => Some("-".to_string()),
            Allele::Sequence(bases) => Some(String::from_utf8_lossy(bases).to_ascii_uppercase()),
            _ => None,
        }
    };
    let mut ref_text = text(reference)?;
    let (record_start, record_end) = (start, end);
    let mut alt_texts: Vec<String> = alts.iter().map(|a| text(a).unwrap_or_default()).collect();
    let (mut start, mut end) = (start as i64, end as i64);

    // a lone indel alt is the minimal pair VEP annotates
    let length = |s: &str| if s == "-" { 0 } else { s.len() };
    if let [only] = alt_texts.as_slice() {
        if length(only) != length(&ref_text) {
            let bare = |s: &str| {
                if s == "-" {
                    String::new()
                } else {
                    s.to_string()
                }
            };
            let (r, a, s, e) = trim_pair(&bare(&ref_text), &bare(only), start, end);
            ref_text = if r.is_empty() { "-".into() } else { r };
            alt_texts = vec![if a.is_empty() { "-".into() } else { a }];
            (start, end) = (s, e);
        }
    }
    let original_allele = alt_texts.get(alt_index)?.clone();
    // Only the alt is looked at, as Ensembl does: a VCF REF of `N` still has a genome
    // base to describe, and it is the genome's base that the description names.
    let alt_ok = |s: &str| s.chars().all(|c| matches!(c, 'A' | 'C' | 'G' | 'T' | '-'));
    if original_allele.is_empty() || !alt_ok(&original_allele) {
        return None;
    }

    // `_get_flank_seq`: 100 bases either side, more when an allele is long, to
    // find where the 3'-rule stops and to look for the bases a dup would copy.
    let mut add_length = 100i64;
    for al in std::iter::once(&ref_text).chain(alt_texts.iter()) {
        if al.len() as i64 > 50 {
            add_length = al.len() as i64 + 50;
        }
    }
    let mut seq_start = start - add_length;
    let seq_end = end + add_length;
    let mut ref_start = add_length;
    let mut ref_end = add_length + end - start;
    if seq_start < 0 {
        ref_start += seq_start;
        ref_end += seq_start;
        seq_start = 0;
    }
    let Ok(ref_seq) =
        seq_provider.fetch_sequence(chrom, (seq_start + 1) as u64, seq_end.max(1) as u64)
    else {
        return Some(hgvsg_clipped(
            chrom,
            record_start,
            record_end,
            reference,
            &alts[alt_index],
        ));
    };

    let multi = alt_texts.len() > 1;
    let (mut chr_start, mut chr_end) = (start, end);
    let mut allele = original_allele;
    let mut ref_allele = ref_text.clone();
    let (mut r_start, mut r_end) = (ref_start, ref_end);
    if multi {
        let bare = |s: &str| {
            if s == "-" {
                String::new()
            } else {
                s.to_string()
            }
        };
        let (r, a, s, e) = trim_pair(&bare(&ref_text), &bare(&allele), chr_start, chr_end);
        r_start = ref_start + (s - chr_start);
        r_end = ref_end + (e - chr_end);
        (chr_start, chr_end) = (s, e);
        allele = if a.is_empty() { "-".into() } else { a };
        ref_allele = if r.is_empty() { "-".into() } else { r };
    }

    // only insertions and deletions are 3'-shifted
    let is_deletion = allele == "-" && ref_allele != "-";
    let is_insertion = ref_allele == "-" && allele != "-";
    let (offset, check_allele) = if is_deletion || is_insertion {
        let mut seq_to_check = if is_deletion {
            ref_allele.clone().into_bytes()
        } else {
            allele.clone().into_bytes()
        };
        let downstream = perl_substr(&ref_seq, r_end, None);
        // `get_3prime_seq_offset`
        let mut offset = 0i64;
        for n in 0..=downstream.len() {
            let next_ref = downstream.get(n).copied();
            if seq_to_check.first().copied() != next_ref {
                break;
            }
            offset += 1;
            seq_to_check.remove(0);
            seq_to_check.push(next_ref.expect("equal to a base"));
        }
        let shifted = if is_insertion {
            String::from_utf8_lossy(&seq_to_check).into_owned()
        } else {
            "-".to_string()
        };
        (offset, shifted)
    } else {
        (0, allele.clone())
    };

    let mut notation = hgvs_variant_notation(
        &check_allele,
        &ref_seq,
        r_start + offset,
        r_end + offset,
        chr_start + offset,
        chr_end + offset,
    )?;
    if notation.kind == NotationKind::Delins {
        notation = clip_alleles(notation);
    }
    format_hgvs_g(chrom, &notation)
}

/// The `g.` description of each alt allele of a record, which belongs to the
/// variant and is the same on every transcript row, so it is made once per record
/// and not once per row.
///
/// Ensembl's when there is a reference to read ([`hgvsg_vep`]), the plain
/// description of [`hgvsg_clipped`] when there is not.
///
/// An allele that is not plain bases has none, with or without a reference: `*` (a
/// deletion that began upstream), a symbolic allele, a base that is not one of
/// `ACGT`. Ensembl skips those (`next if $allele =~ m/[^ACGT\-]/ig` in
/// `hgvs_genomic`), and the plain description had nothing to say for them but
/// `chr10:g.1400?`, which is not HGVS.
pub fn hgvsg_per_alt(
    seq_provider: Option<&dyn SequenceProvider>,
    chrom: &str,
    start: u64,
    end: u64,
    reference: &fastvep_core::Allele,
    alts: &[fastvep_core::Allele],
) -> Vec<Option<String>> {
    use fastvep_core::Allele;
    let describable = |alt: &Allele| match alt {
        Allele::Deletion => true,
        Allele::Sequence(bases) => bases
            .iter()
            .all(|b| matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T')),
        Allele::Missing | Allele::Symbolic(_) => false,
    };
    alts.iter()
        .enumerate()
        .map(|(i, alt)| {
            if !describable(alt) {
                return None;
            }
            match seq_provider {
                Some(sp) => hgvsg_vep(sp, chrom, start, end, reference, alts, i),
                None => Some(hgvsg_clipped(chrom, start, end, reference, alt)),
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use anyhow::{anyhow, Result};
    use fastvep_core::HgvsPlacement;
    use fastvep_core::{Allele, Strand};
    use fastvep_genome::SpanSite;
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

    /// 20 exonic bases, an intron opening `AGGT AGGT` and then not repeating,
    /// then exon 2. A four-base deletion at the donor can be written at `+1` or
    /// at `+5`, and both edit the sequence the same way. Exactly as long as the
    /// transcript it is read against, which `the_fixtures_match_the_transcript`
    /// holds it to: a reference that drifts from the coordinates it stands for
    /// makes a placement test pass or fail on filler.
    const DONOR_REPEAT: &str = "AAAAAAAAAAAAAAAAAAAA\
                                AGGTAGGTCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\
                                GGGGGGGGGGGGGGGGGGGG";

    /// The exon ends `AGGT` and the intron opens `AGGT`, so the repeat straddles
    /// the boundary and a duplication of its unit can be anchored inside the
    /// exon or on its last base.
    const EXON_END_REPEAT: &str = "AAAAAAAAAAAAAAAAAGGT\
                                   AGGTCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\
                                   GGGGGGGGGGGGGGGGGGGG";

    /// An intron whose last bases are `CTAG`, then an exon that opens `GAGGCC`:
    /// the MLH1 / RYR1 / ATM acceptor shape. The intron is 21..=80 and `TAG`
    /// is 78..=80; exon 2 is 81..=100.
    const ACCEPTOR_REPEAT: &str = "AAAAAAAAAAAAAAAAAAAA\
                                   CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCTAG\
                                   GAGGCCCCCCCCCCCCCCCC";

    #[test]
    fn the_fixtures_match_the_transcript() {
        for fixture in [DONOR_REPEAT, EXON_END_REPEAT, ACCEPTOR_REPEAT] {
            assert_eq!(fixture.len(), 100, "{fixture}");
        }
    }

    /// Where the change *is* and where HGVS *places* it are two different
    /// questions, and PVS1 needs the second: a deletion that can be written at
    /// `+5` can be written without touching `+1` or `+2`, so the canonical `GT`
    /// survives it whichever alignment is chosen.
    #[test]
    fn the_placement_is_the_furthest_the_change_reaches_on_a_donor() {
        let r = StrRef(DONOR_REPEAT);
        let tr = transcript(Strand::Forward);
        let (del_ref, del_alt) = (Allele::Sequence(b"AGGT".to_vec()), Allele::Deletion);
        // Written at the donor's own first base.
        assert_eq!(tr.intronic_offset_covered(21, 24), Some(1));
        // Placed four bases further in, clear of the dinucleotide.
        assert_eq!(
            hgvs_placement(Some(&r), "1", &tr, 21, 24, &del_ref, &del_alt),
            HgvsPlacement::Intronic { offset: 5 }
        );
    }

    /// Without a reference there is no walk to make, so the change is placed
    /// where it sits - the same fallback the description itself takes.
    #[test]
    fn the_placement_falls_back_to_the_position_when_there_is_no_reference() {
        let tr = transcript(Strand::Forward);
        let (del_ref, del_alt) = (Allele::Sequence(b"AGGT".to_vec()), Allele::Deletion);
        assert_eq!(
            hgvs_placement(None, "1", &tr, 21, 24, &del_ref, &del_alt),
            HgvsPlacement::Intronic { offset: 1 }
        );
    }

    /// A change inside an exon is exonic, one off the transcript is given no
    /// position on it, and a substitution has nowhere to go.
    /// A sequence ALT beside a symbolic one takes its `end` from INFO END, which
    /// belongs to the symbolic allele and can be smaller than the clip. The
    /// subtraction used to underflow, which aborts a debug or test build on a
    /// line of a VCF (`1 30 . CTCCTCCT GATGCT,<DEL> . . END=1;SVTYPE=DEL`).
    #[test]
    fn a_span_that_ends_before_the_clip_does_not_underflow() {
        let r = Allele::Sequence(b"CTCCTCCT".to_vec());
        let a = Allele::Sequence(b"GATGCT".to_vec());
        // The back clip is two bases, and `end` is 1.
        assert_eq!(clipped_span(Strand::Forward, 30, 1, &r, &a), (30, 0));
        let hv = hgvs_allele(Strand::Forward, 30, 1, &r, &a, None);
        assert_eq!((hv.start, hv.end), (30, 0));
        let _ = hgvsg_clipped("1", 30, 1, &r, &a);
        // There is no change on the transcript to place, so the terms decide.
        for strand in [Strand::Forward, Strand::Reverse] {
            assert_eq!(
                hgvs_placement(None, "1", &transcript(strand), 30, 1, &r, &a),
                HgvsPlacement::Spans
            );
        }
        // A pure insertion's span ends one base before it starts, and is placed.
        let ins = (
            Allele::Sequence(b"G".to_vec()),
            Allele::Sequence(b"GTT".to_vec()),
        );
        assert_eq!(
            hgvs_placement(
                None,
                "1",
                &transcript(Strand::Forward),
                10,
                10,
                &ins.0,
                &ins.1
            ),
            HgvsPlacement::Exonic { in_frame: false }
        );
    }

    #[test]
    fn the_placement_names_an_exonic_change_and_declines_one_off_the_transcript() {
        let r = StrRef(DONOR_REPEAT);
        let tr = transcript(Strand::Forward);
        let snv = (
            Allele::Sequence(b"A".to_vec()),
            Allele::Sequence(b"G".to_vec()),
        );
        assert_eq!(
            hgvs_placement(Some(&r), "1", &tr, 10, 10, &snv.0, &snv.1),
            HgvsPlacement::Exonic { in_frame: true }
        );
        assert_eq!(
            hgvs_placement(Some(&r), "1", &tr, 500, 500, &snv.0, &snv.1),
            HgvsPlacement::Spans
        );
    }

    /// The MSH6 `c.3801+3_3801+6dup` shape. The record inserts `AGGT` after the
    /// exon's last base, where the intron already opens `AGGT AGGT`: it is
    /// written at the exon boundary and placed past the second copy, `+8`, with
    /// the exon never touched.
    #[test]
    fn an_insertion_anchored_on_the_last_exon_base_is_placed_inside_the_intron() {
        let r = StrRef(DONOR_REPEAT);
        let tr = transcript(Strand::Forward);
        let (ins_ref, ins_alt) = (Allele::Deletion, Allele::Sequence(b"AGGT".to_vec()));
        // Between exon base 20 and intron base 21: on the boundary, where it sits.
        assert_eq!(tr.span_site(20, 21), SpanSite::Spans);
        assert_eq!(
            hgvs_placement(Some(&r), "1", &tr, 21, 20, &ins_ref, &ins_alt),
            HgvsPlacement::Intronic { offset: 8 }
        );
    }

    /// The same duplication written at either end of a repeat that straddles the
    /// exon boundary. A caller that left-aligns anchors inside the exon, as
    /// NBEAL2 `c.6560+4_6560+20dup` arrives, and one that does not anchors on
    /// the exon's last base, as MSH6 `c.3801+3_3801+6dup` does. The record
    /// anchored inside the exon reaches no intronic base the way it is written,
    /// so nothing about it says where HGVS will put it; HGVS names one place for
    /// both.
    #[test]
    fn the_placement_does_not_depend_on_which_end_of_the_repeat_the_record_anchors_at() {
        let r = StrRef(EXON_END_REPEAT);
        let tr = transcript(Strand::Forward);
        let (ins_ref, ins_alt) = (Allele::Deletion, Allele::Sequence(b"AGGT".to_vec()));
        assert_eq!(tr.intronic_offset_covered(17, 16), None);
        for (start, end) in [(21, 20), (17, 16)] {
            assert_eq!(
                hgvs_placement(Some(&r), "1", &tr, start, end, &ins_ref, &ins_alt),
                HgvsPlacement::Intronic { offset: 4 },
                "anchored between {end} and {start}"
            );
        }
    }

    /// A multi-allelic record is trimmed once for the whole site, so one allele
    /// can arrive as a replacement that clips to a plain deletion. It is placed
    /// the way its HGVSc is described, not where the unclipped pair would put it.
    #[test]
    fn a_replacement_that_clips_to_a_deletion_is_placed_like_the_deletion() {
        let r = StrRef(DONOR_REPEAT);
        let tr = transcript(Strand::Forward);
        // `AAGGT>AA`... clips the shared leading `A`s, leaving `AGGT>-` at 21..24.
        let (rep_ref, rep_alt) = (
            Allele::Sequence(b"AAGGT".to_vec()),
            Allele::Sequence(b"A".to_vec()),
        );
        assert_eq!(
            hgvs_placement(Some(&r), "1", &tr, 20, 24, &rep_ref, &rep_alt),
            HgvsPlacement::Intronic { offset: 5 }
        );
    }

    /// `TAGG>T` deletes the intron's `AG` and the exon's first base, so the
    /// record is a splice acceptor variant. Exon 2 opens `GAGG`, so the same
    /// sequence is the exon's `AGG` deleted, in frame, with the acceptor `AG`
    /// untouched - and the shift runs toward the exon, so HGVS writes it there.
    #[test]
    fn an_acceptor_deletion_that_the_shift_carries_into_the_exon_is_exonic_and_in_frame() {
        let r = StrRef(ACCEPTOR_REPEAT);
        let tr = transcript(Strand::Forward);
        let (del_ref, del_alt) = (Allele::Sequence(b"AGG".to_vec()), Allele::Deletion);
        // Where it sits: the acceptor's last two bases and the exon's first.
        assert_eq!(tr.intronic_offset_covered(79, 81), Some(-1));
        assert_eq!(
            hgvs_placement(Some(&r), "1", &tr, 79, 81, &del_ref, &del_alt),
            HgvsPlacement::Exonic { in_frame: true }
        );
        // One base fewer cannot slide anywhere: the exon's `G` is not an `A`.
        let (one_ref, one_alt) = (Allele::Sequence(b"AG".to_vec()), Allele::Deletion);
        assert_eq!(
            hgvs_placement(Some(&r), "1", &tr, 79, 80, &one_ref, &one_alt),
            HgvsPlacement::Intronic { offset: -1 }
        );
    }

    /// The reverse strand reads the same repeat from the other end: the donor is
    /// at the intron's high coordinate and the 3' shift runs down it. The
    /// reference is the forward fixture reverse-complemented and every position
    /// mirrored, so this is DONOR_REPEAT as the reverse transcript sees it.
    #[test]
    fn the_placement_is_the_same_on_the_reverse_strand() {
        fn mirrored(s: &str) -> &'static str {
            let rc: String = s
                .bytes()
                .rev()
                .map(|b| match b {
                    b'A' => 'T',
                    b'T' => 'A',
                    b'C' => 'G',
                    b'G' => 'C',
                    other => other as char,
                })
                .collect();
            Box::leak(rc.into_boxed_str())
        }
        let r = StrRef(mirrored(DONOR_REPEAT));
        let tr = transcript(Strand::Reverse);
        // AGGT at 21..=24 forward is ACCT at 77..=80 here, written at the donor.
        let (del_ref, del_alt) = (Allele::Sequence(b"ACCT".to_vec()), Allele::Deletion);
        assert_eq!(tr.intronic_offset_covered(77, 80), Some(1));
        assert_eq!(
            hgvs_placement(Some(&r), "1", &tr, 77, 80, &del_ref, &del_alt),
            HgvsPlacement::Intronic { offset: 5 }
        );
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

    /// A multi-allelic record's allele that trims to a substitution is written as the
    /// substitution it is.
    ///
    /// VEP 116.3 writes `g.12delinsAA` for the `TAA` allele of `TGA>TG,TAA`: its
    /// `hgvs_genomic` trims the reference span to the one base that differs and then
    /// writes the *untrimmed* alt over it, which is a different change (it inserts a base
    /// the variant does not have). `g.12G>A` is the change, and what is written here.
    /// `docs/VEP_DIVERGENCE.md` lists it, and the generated worlds leave such records out.
    #[test]
    fn an_allele_that_trims_to_a_substitution_is_a_substitution() {
        let reference = StrRef("CCCCCCCCCCTGACCCCCCCCCC");
        // as the VCF parser leaves it: the shared first base taken off, so GA > G, AA
        let alts = [
            fastvep_core::Allele::Sequence(b"G".to_vec()),
            fastvep_core::Allele::Sequence(b"AA".to_vec()),
        ];
        let got = hgvsg_per_alt(
            Some(&reference),
            "chr1",
            12,
            13,
            &fastvep_core::Allele::Sequence(b"GA".to_vec()),
            &alts,
        );
        assert_eq!(
            got,
            vec![
                Some("chr1:g.13del".to_string()),
                Some("chr1:g.12G>A".to_string())
            ]
        );
    }

    /// The description is read off the genome, not off the VCF's REF: where they differ
    /// Ensembl describes the genome, and when the alt is the genome's own base it
    /// writes nothing.
    #[test]
    fn hgvsg_describes_the_genome_and_writes_nothing_for_its_own_base() {
        let reference = StrRef("ACGTACGTACGTACGT");
        let alt_t = [fastvep_core::Allele::Sequence(b"T".to_vec())];
        // position 3 is G in the genome whatever the VCF says its REF is
        let got = hgvsg_per_alt(
            Some(&reference),
            "chr1",
            3,
            3,
            &fastvep_core::Allele::Sequence(b"A".to_vec()),
            &alt_t,
        );
        assert_eq!(got, vec![Some("chr1:g.3G>T".to_string())]);
        // position 4 is T already
        let got = hgvsg_per_alt(
            Some(&reference),
            "chr1",
            4,
            4,
            &fastvep_core::Allele::Sequence(b"A".to_vec()),
            &alt_t,
        );
        assert_eq!(got, vec![None]);
        // and a contig named MT is numbered `m.`, as VEP numbers it
        let got = hgvsg_per_alt(
            Some(&reference),
            "MT",
            3,
            3,
            &fastvep_core::Allele::Sequence(b"G".to_vec()),
            &alt_t,
        );
        assert_eq!(got, vec![Some("MT:m.3G>T".to_string())]);
    }

    /// Without a reference the description is the plain one, where the VCF put the
    /// variant: what every run gave before the genome was read.
    #[test]
    fn hgvsg_without_a_reference_is_the_plain_description() {
        let alts = [fastvep_core::Allele::Sequence(b"T".to_vec())];
        let got = hgvsg_per_alt(
            None,
            "chr1",
            3,
            3,
            &fastvep_core::Allele::Sequence(b"A".to_vec()),
            &alts,
        );
        assert_eq!(got, vec![Some("chr1:g.3A>T".to_string())]);
    }

    /// `HGVS_OFFSET` of a deletion counts from where the VCF anchors it, whichever
    /// way the pair is clipped.
    ///
    /// On the reverse strand the pair is clipped from the transcript's end, so
    /// `CTGCAC` to `C` is the deletion of `CTGCA`, not `TGCAC`, and a record that
    /// carries a second ALT (which keeps the anchor `C`, where a lone ALT has had it
    /// taken off) used to report the shift one base short of the one its lone twin
    /// got: 2 against 3, for the same description. Found by running both through
    /// Ensembl VEP 116.3, which writes `-` for the multi-allelic one and 3 for the
    /// lone one.
    #[test]
    fn hgvs_offset_counts_from_where_the_vcf_anchors_the_change() {
        use fastvep_core::Allele;
        // 36..=43 is `CACTGCAC`, the base before it a `T` that ends the repeat.
        let mut genome = vec![b'A'; 100];
        genome[34] = b'T';
        genome[35..43].copy_from_slice(b"CACTGCAC");
        let genome = Box::leak(String::from_utf8(genome).unwrap().into_boxed_str());
        let reference = StrRef(genome);
        let tr = transcript(Strand::Reverse);
        let seq = |s: &[u8]| Allele::Sequence(s.to_vec());

        let shift_of = |start: u64, end: u64, r: Allele, a: Allele| {
            let hv = hgvs_allele(Strand::Reverse, start, end, &r, &a, None);
            let got = hgvsc_shifted_with_shift(
                Some(&reference as &dyn SequenceProvider),
                "1",
                &tr,
                "T.1",
                hv.start,
                hv.end,
                hv.anchor_start,
                &hv.genomic_ref,
                &hv.genomic_alt,
                &hv.hgvs_ref,
                &hv.hgvs_alt,
                Some(1),
                Some(40),
            )
            .expect("a description");
            (hv.start, hv.anchor_start, got)
        };

        // the lone ALT: the parser has taken the anchor off, `TGCAC` at 39..=43
        let (start, anchor, (lone_hgvs, lone_shift)) =
            shift_of(39, 43, seq(b"TGCAC"), Allele::Deletion);
        assert_eq!((start, anchor, lone_shift), (39, 39, 3));
        // a record with several ALTs: `CTGCAC` to `C` at 38..=43, clipped from the back
        let (start, anchor, (multi_hgvs, multi_shift)) =
            shift_of(38, 43, seq(b"CTGCAC"), seq(b"C"));
        assert_eq!((start, anchor), (38, 39), "clipped at 38, anchored at 39");
        assert_eq!(multi_hgvs, lone_hgvs, "the same description");
        assert_eq!(multi_shift, lone_shift, "and the same distance travelled");
    }

    /// An allele that is not plain bases has no HGVSg, as in Ensembl, whether or not
    /// there is a reference to read: the plain description would have been `g.3?`.
    /// The alleles beside it keep theirs.
    ///
    /// Ensembl 116.3 on the HG002 sample writes no row at all for a `*` allele, and
    /// `HGVSg` for the allele beside it as usual (`6:31318928 A>T,*`).
    #[test]
    fn an_allele_that_is_not_plain_bases_has_no_hgvsg() {
        use fastvep_core::Allele;
        let reference = Allele::Sequence(b"A".to_vec());
        let alts = [
            Allele::Sequence(b"T".to_vec()),
            Allele::Missing,
            Allele::Symbolic("<DEL>".into()),
            Allele::Sequence(b"N".to_vec()),
        ];
        let genome = StrRef("CCAGG");
        // a VCF REF of `N` is not a reason to say nothing: Ensembl names the genome's base
        // (`g.1600C>A` for the REF `N` of a record in which the genome has a `C`), and a
        // run with no reference can only say what the VCF said
        let n_ref = Allele::Sequence(b"N".to_vec());
        let only_t = [Allele::Sequence(b"T".to_vec())];
        assert_eq!(
            hgvsg_per_alt(
                Some(&genome as &dyn SequenceProvider),
                "chr1",
                3,
                3,
                &n_ref,
                &only_t
            ),
            vec![Some("chr1:g.3A>T".to_string())]
        );
        assert_eq!(
            hgvsg_per_alt(None, "chr1", 3, 3, &n_ref, &only_t),
            vec![Some("chr1:g.3N>T".to_string())]
        );
        for provider in [None, Some(&genome as &dyn SequenceProvider)] {
            let got = hgvsg_per_alt(provider, "chr1", 3, 3, &reference, &alts);
            assert_eq!(
                got,
                vec![Some("chr1:g.3A>T".to_string()), None, None, None],
                "{}",
                if provider.is_some() {
                    "with a reference"
                } else {
                    "without one"
                }
            );
        }
    }
}
