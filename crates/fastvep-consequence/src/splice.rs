use fastvep_core::{Allele, Strand};
use fastvep_genome::Transcript;

// Splice site boundaries relative to the intron, in genomic coordinates.
// `intron_start` is the intron's first base and `intron_end` its last, so
// every range below is written the way Ensembl writes it in `_intron_effects`
// (`BaseTranscriptVariationAllele.pm` release/115) and the strand enters only
// when the flags are turned into SO terms:
//
//   start_splice_site       intron_start   .. intron_start+1   (donor on +)
//   end_splice_site         intron_end-1   .. intron_end       (acceptor on +)
//   fifth_base              intron_start+4                     (donor 5th on +)
//   donor_region            intron_start+2 .. intron_start+5
//   polypyrimidine          intron_end-16  .. intron_end-2
//   fifth_base_reverse      intron_end-4                       (donor 5th on -)
//   donor_region_reverse    intron_end-5   .. intron_end-2
//   polypyrimidine_reverse  intron_start+2 .. intron_start+16

/// The splice terms one variant's span earns against one transcript.
///
/// Ensembl computes these together rather than as independent predicates, and
/// the difference matters: the suppression rules are not a hierarchy. A
/// `splice_donor_variant` and a `splice_donor_5th_base_variant` coexist on the
/// same variant, and so do a `splice_polypyrimidine_tract_variant` and a
/// `splice_region_variant`; only `splice_region_variant` is displaced by the
/// more specific terms, and only `splice_donor_region_variant` is displaced by
/// the 5th base. Suppressing all four extended terms behind an essential site -
/// which is what reading them as a fallback chain amounts to - dropped 812 of
/// them over a 6,600-variant ClinVar sample. With every rule below in place,
/// the splice terms match real VEP 115.1 on all 150,725 rows of that sample.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct SpliceEffects {
    /// Ensembl's `within_intron`: the variant reaches the intron's interior,
    /// past the two bases of each splice site. It is not a splice term, but it
    /// comes from the same pass and answers the same kind of question, and it is
    /// subject to none of the suppression below - `intron_variant` sits beside a
    /// `splice_donor_variant` and beside an exonic term alike.
    pub intronic: bool,
    pub donor: bool,
    pub acceptor: bool,
    pub donor_5th_base: bool,
    pub donor_region: bool,
    pub polypyrimidine_tract: bool,
    pub region: bool,
}

/// Whether `[var_start, var_end]` overlaps `[site_start, site_end]`.
///
/// Ensembl's `overlap`, and the insertion convention falls out of it rather than
/// needing a case: an insertion arrives as the zero-length interval
/// `end = start - 1`, so the test narrows to `site_start < start <= site_end` -
/// an insertion counts only where it sits *inside* the site, not where it abuts
/// it. Do not sort the pair before calling; the order carries that meaning.
fn overlaps(var_start: u64, var_end: u64, site_start: u64, site_end: u64) -> bool {
    var_start <= site_end && var_end >= site_start
}

/// Every splice term Ensembl gives `[var_start, var_end]` against `transcript`.
///
/// A port of `_intron_effects` plus the SO-term predicates in
/// `Utils/VariationEffect.pm` that read its flags. Two properties of that code
/// are load-bearing and easy to lose by rewriting it as separate predicates.
///
/// Every test is an overlap rather than an equality, so a variant whose second
/// base lands on a site is on that site. Reading `var_start` alone made BRCA2
/// `c.9256_9256+1delinsTA` - the exon's last base and the intron's first,
/// changed together - a `splice_region_variant`, LOW where VEP says HIGH.
///
/// And what is matched against a site is not the variant's span but the parts
/// of it that differ; see [`for_each_differing_region`].
pub fn splice_effects(
    transcript: &Transcript,
    var_start: u64,
    var_end: u64,
    ref_allele: &Allele,
    alt_allele: &Allele,
) -> SpliceEffects {
    // Which introns are looked at is decided by the variant's *whole* span, not
    // by the differing regions: Ensembl selects the intron lists once from
    // `$vf->{start}..{end}` and then applies the per-region tests only to what
    // that selection returned. The two are not interchangeable, because the
    // padding in `differing_regions` can push a region past the span. Reading
    // the regions for selection too reported `splice_donor_variant` for a
    // frameshift five bases inside an exon, on the strength of padding that
    // reached the donor.
    //
    // `splice_polypyrimidine_tract_variant` is gated once more, before any
    // predicate runs: `include => { exon => 0, intron => 1 }` in
    // `Utils/Config.pm`, so a delins that starts in an exon and runs into the
    // tract earns an acceptor term and no tract term.
    let mut touches_exon = false;
    for exon in &transcript.exons {
        touches_exon |= overlaps(var_start, var_end, exon.start, exon.end);
    }

    let mut start_site = false;
    let mut end_site = false;
    let mut fifth = false;
    let mut fifth_rev = false;
    let mut donor_region = false;
    let mut donor_region_rev = false;
    let mut polypy = false;
    let mut polypy_rev = false;
    let mut region = false;
    let mut intronic = false;

    // Regions outside, introns inside - the nesting Ensembl uses, and it is
    // observable. `splice_region` is *assigned* rather than or-assigned on each
    // pass, so the last (region, intron) pair decides it: a change whose first
    // differing base is in the splice region and whose last one is not comes out
    // with no splice term at all. That reads like a bug, but it is the answer
    // VEP gives, on 63 rows over a 6,600-variant ClinVar sample, and guessing
    // the tidier rule instead is how a `missense_variant` acquires a
    // `splice_region_variant` VEP never reports.
    for_each_differing_region(
        var_start,
        var_end,
        ref_allele,
        alt_allele,
        |r_start, r_end| {
            // Ensembl's zero-length interval for an insertion; see `overlaps`.
            let insertion = r_end + 1 == r_start;
            let (lo, hi) = (r_start.min(r_end), r_start.max(r_end));
            for_each_intron(transcript, |intron_start, intron_end| {
                // An intron short enough to be a frameshift intron - 12 bases or
                // fewer - is skipped whole: Ensembl treats a variant inside one as
                // exonic rather than spliced.
                if intron_end - intron_start < 12
                    && overlaps(r_start, r_end, intron_start, intron_end)
                {
                    return;
                }
                let back = |n: u64| intron_end.saturating_sub(n);
                // `_overlapped_introns`: the intron plus a 3-base flank. Only the
                // polypyrimidine tract and `intronic` are read from this set.
                if overlaps(
                    var_start,
                    var_end,
                    intron_start.saturating_sub(3),
                    intron_end + 3,
                ) {
                    intronic |= overlaps(r_start, r_end, intron_start + 2, back(2))
                        || (insertion && (r_start == intron_start + 2 || r_end == back(2)));
                    // The two polypyrimidine tests are the only ones Ensembl runs
                    // against the *sorted* pair - `_intron_effects` swaps `$start`
                    // and `$end` before them and nowhere else. For an insertion,
                    // whose interval is `end = start - 1`, sorting turns the test
                    // from "strictly inside" into "abuts or inside", so an insertion
                    // sitting on the tract's first base is in the tract. That is a
                    // real distinction rather than an accident: the tract is a
                    // stretch of sequence, and inserting at either edge lengthens
                    // it. 302 rows over a 6,600-variant ClinVar sample.
                    polypy |= overlaps(lo, hi, back(16), back(2));
                    polypy_rev |= overlaps(lo, hi, intron_start + 2, intron_start + 16);
                }
                // `_overlapped_introns_boundary`: the two ends of the intron. Every
                // other flag is read from this set. Which introns are in either set
                // is decided by the variant's *whole* span, not by the differing
                // region, because Ensembl selects the lists once from
                // `$vf->{start}..{end}` and then applies the per-region tests to
                // what that selection returned. Reading the region for selection too
                // reported `splice_donor_variant` for a frameshift five bases inside
                // an exon, on the strength of padding that reached the donor.
                if !(overlaps(
                    var_start,
                    var_end,
                    intron_start.saturating_sub(3),
                    intron_start + 7,
                ) || overlaps(
                    var_start,
                    var_end,
                    intron_end.saturating_sub(7),
                    intron_end + 3,
                )) {
                    return;
                }
                start_site |= overlaps(r_start, r_end, intron_start, intron_start + 1);
                end_site |= overlaps(r_start, r_end, back(1), intron_end);
                fifth |= overlaps(r_start, r_end, intron_start + 4, intron_start + 4);
                fifth_rev |= overlaps(r_start, r_end, back(4), back(4));
                donor_region |= overlaps(r_start, r_end, intron_start + 2, intron_start + 5);
                donor_region_rev |= overlaps(r_start, r_end, back(5), back(2));
                if !start_site && !end_site {
                    region = intron_overlap(r_start, r_end, intron_start, intron_end, insertion);
                }
            });
        },
    );

    let forward = transcript.strand == Strand::Forward;
    let donor = if forward { start_site } else { end_site };
    let acceptor = if forward { end_site } else { start_site };
    let donor_5th_base = if forward { fifth } else { fifth_rev };
    // `splice_donor_region_variant` defers to the 5th base and to nothing else.
    let donor_region = !donor_5th_base
        && if forward {
            donor_region
        } else {
            donor_region_rev
        };
    let polypyrimidine_tract = !touches_exon && if forward { polypy } else { polypy_rev };
    // `splice_region_variant` is the only term the others displace, and it is
    // displaced by all four of them.
    let region = region && !donor && !acceptor && !donor_5th_base && !donor_region;

    SpliceEffects {
        intronic,
        donor,
        acceptor,
        donor_5th_base,
        donor_region,
        polypyrimidine_tract,
        region,
    }
}

/// The genomic stretches of a variant that Ensembl actually tests against a
/// splice site: the runs of bases where the reference and alternate strings
/// differ, not the reference allele's whole span.
///
/// A port of `_get_differing_regions` (`VariationFeatureOverlapAllele.pm`
/// release/115). The two strings are compared position by position with the
/// shorter padded, so `CCT/TCTG` yields two one-base regions - offsets 0 and 3 -
/// rather than one three-base span. That is not a refinement of the span: it is
/// both narrower (the matching interior is not tested) and *wider* (the padding
/// makes every base past the reference allele's end differ, so `TAC/GCAAAAA`
/// is tested over seven bases, four of them past the reference). Reading the
/// plain span instead reported `splice_acceptor_variant` for a delins whose
/// changed bases straddle the site without touching it, and missed the donor
/// on a delins whose insertion runs into it.
///
/// A variant with one side of length 1 or less - every SNV, every pure
/// insertion, every pure deletion - is a single region covering the reference
/// span, which for an insertion is Ensembl's zero-length `end = start - 1`.
fn for_each_differing_region<F>(
    var_start: u64,
    var_end: u64,
    ref_allele: &Allele,
    alt_allele: &Allele,
    mut visit: F,
) where
    F: FnMut(u64, u64),
{
    let (ref_bases, alt_bases) = (ref_allele.as_bytes(), alt_allele.as_bytes());
    // Ensembl measures the reference side from the coordinates, not the string,
    // so a zero-length insertion interval reports 0 here.
    let ref_len = (var_end + 1).saturating_sub(var_start) as usize;
    if ref_len <= 1 || alt_bases.len() <= 1 {
        visit(var_start, var_end);
        return;
    }

    // A run at a time, so nothing is collected: this runs once per
    // (variant x transcript x allele), where a `Vec` would be one allocation
    // per SNV for the single region every SNV has.
    let mut run: Option<(u64, u64)> = None;
    for i in 0..ref_len.max(alt_bases.len()) {
        if ref_bases.get(i) == alt_bases.get(i) {
            continue;
        }
        let pos = var_start + i as u64;
        match run {
            Some((s, e)) if e + 1 == pos => run = Some((s, pos)),
            Some((s, e)) => {
                visit(s, e);
                run = Some((pos, pos));
            }
            None => run = Some((pos, pos)),
        }
    }
    if let Some((s, e)) = run {
        visit(s, e);
    }
}

/// SO:0001630, "within 1-3 bases of the exon or 3-8 bases of the intron".
///
/// A port of Ensembl's `_intron_overlap` (`Utils/VariationEffect.pm`). It is
/// purely genomic - the four ranges are symmetric about each intron, so the
/// strand never enters. And it has an insertion clause: an insertion sitting
/// exactly on an intron edge, or on the inner edge of a donor or acceptor
/// dinucleotide, is in the splice region even though its zero-length interval
/// overlaps none of the four ranges. Without that clause a `c.830+2dup` lands
/// in no splice term at all - MODIFIER where VEP says LOW.
fn intron_overlap(
    var_start: u64,
    var_end: u64,
    intron_start: u64,
    intron_end: u64,
    insertion: bool,
) -> bool {
    overlaps(var_start, var_end, intron_start + 2, intron_start + 7)
        || overlaps(
            var_start,
            var_end,
            intron_end.saturating_sub(7),
            intron_end.saturating_sub(2),
        )
        || overlaps(
            var_start,
            var_end,
            intron_start.saturating_sub(3),
            intron_start.saturating_sub(1),
        )
        || overlaps(var_start, var_end, intron_end + 1, intron_end + 3)
        || (insertion
            && (var_start == intron_start
                || var_end == intron_end
                || var_start == intron_start + 2
                || var_end == intron_end.saturating_sub(2)))
}

/// Call `visit` with the genomic bounds of every intron, in transcript order.
///
/// An intron shorter than the 17-base windows above makes those windows overlap
/// each other; that is Ensembl's arithmetic too, and the flags stay meaningful
/// because every range is clamped by the overlap test itself. The offsets that
/// reach backwards from `intron_end` saturate rather than wrap, which is what a
/// negative lower bound does in Perl.
fn for_each_intron<F>(transcript: &Transcript, mut visit: F)
where
    F: FnMut(u64, u64),
{
    // Exons are stored in transcript order (`gff.rs` sorts them on load), so
    // adjacent entries bound one intron whichever way the transcript runs.
    // Reading the pair through `min`/`max` rather than sorting a copy keeps this
    // free of the per-(variant x transcript) allocation a `Vec` would cost.
    for pair in transcript.exons.windows(2) {
        let (intron_start, intron_end) = (
            pair[0].end.min(pair[1].end) + 1,
            pair[0].start.max(pair[1].start) - 1,
        );
        // A RefSeq transcript with overlapping exons yields start > end.
        if intron_start > intron_end {
            continue;
        }
        visit(intron_start, intron_end);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fastvep_genome::{Exon, Gene, Transcript, Translation};

    // One-site views of `splice_effects`, so each test reads as the question it
    // is asking. A single-base position is the degenerate span `(p, p)`.
    /// A single-base substitution over `(s, e)`, so each test reads as the
    /// question it is asking rather than as an allele pair. An SNV and an
    /// insertion both give one differing region covering the span, which is what
    /// these fixtures mean to exercise.
    fn effects(tr: &Transcript, s: u64, e: u64) -> SpliceEffects {
        let n = (e + 1).saturating_sub(s).max(1) as usize;
        splice_effects(
            tr,
            s,
            e,
            &Allele::Sequence(vec![b'A'; n]),
            &Allele::Sequence(vec![b'T'; n]),
        )
    }
    fn donor(tr: &Transcript, s: u64, e: u64) -> bool {
        effects(tr, s, e).donor
    }
    fn acceptor(tr: &Transcript, s: u64, e: u64) -> bool {
        effects(tr, s, e).acceptor
    }
    fn region(tr: &Transcript, s: u64, e: u64) -> bool {
        effects(tr, s, e).region
    }
    fn polypy(tr: &Transcript, p: u64) -> bool {
        effects(tr, p, p).polypyrimidine_tract
    }

    fn make_forward_transcript() -> Transcript {
        // Exon1: 1000-1200, Intron1: 1201-1999, Exon2: 2000-2300
        Transcript {
            stable_id: "ENST_TEST".into(),
            version: None,
            gene: Gene {
                stable_id: "ENSG_TEST".into(),
                symbol: None,
                symbol_source: None,
                hgnc_id: None,
                biotype: "protein_coding".into(),
                chromosome: "chr1".into(),
                start: 1000,
                end: 2300,
                strand: Strand::Forward,
            },
            biotype: "protein_coding".into(),
            chromosome: "chr1".into(),
            start: 1000,
            end: 2300,
            strand: Strand::Forward,
            exons: vec![
                Exon {
                    stable_id: "E1".into(),
                    start: 1000,
                    end: 1200,
                    strand: Strand::Forward,
                    phase: 0,
                    end_phase: 0,
                    rank: 1,
                },
                Exon {
                    stable_id: "E2".into(),
                    start: 2000,
                    end: 2300,
                    strand: Strand::Forward,
                    phase: 0,
                    end_phase: 0,
                    rank: 2,
                },
            ],
            translation: Some(Translation {
                stable_id: "P1".into(),
                genomic_start: 1000,
                genomic_end: 2300,
                start_exon_rank: 1,
                start_exon_offset: 0,
                end_exon_rank: 2,
                end_exon_offset: 300,
            }),
            cdna_coding_start: Some(1),
            cdna_coding_end: Some(502),
            coding_region_start: Some(1000),
            coding_region_end: Some(2300),
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
    fn test_splice_donor() {
        let tr = make_forward_transcript();
        // Intron: 1201-1999. Donor = first 2 bases: 1201, 1202
        assert!(donor(&tr, 1201, 1201));
        assert!(donor(&tr, 1202, 1202));
        assert!(!donor(&tr, 1203, 1203));
        assert!(!donor(&tr, 1200, 1200)); // exonic
    }

    #[test]
    fn test_splice_acceptor() {
        let tr = make_forward_transcript();
        // Intron: 1201-1999. Acceptor = last 2 bases: 1998, 1999
        assert!(acceptor(&tr, 1998, 1998));
        assert!(acceptor(&tr, 1999, 1999));
        assert!(!acceptor(&tr, 1997, 1997));
        assert!(!acceptor(&tr, 2000, 2000)); // exonic
    }

    /// `splice_region_variant` is what is left after the more specific terms
    /// have taken their bases, not a band that overlaps them.
    ///
    /// Ensembl's `splice_region` returns 0 for a variant that already earned
    /// `splice_donor_region_variant` or `splice_donor_5th_base_variant`
    /// (`VariationEffect.pm` release/115, l. 616), so on the intronic donor side
    /// only bases 6 and 7 are left for it - bases 2 to 5 are the donor region and
    /// base 4 is the 5th base. This used to assert `splice_region` at
    /// `intron_start + 2`, which is the donor region.
    #[test]
    fn test_splice_region() {
        let tr = make_forward_transcript();
        // Exon 1 ends at 1200, intron 1 is 1201..1999, exon 2 starts at 2000.
        for pos in [1198u64, 1199, 1200] {
            assert!(region(&tr, pos, pos), "last 3 bases of the exon: {pos}");
        }
        for pos in [2000u64, 2001, 2002] {
            assert!(
                region(&tr, pos, pos),
                "first 3 bases of the next exon: {pos}"
            );
        }
        for pos in [1207u64, 1208] {
            assert!(region(&tr, pos, pos), "intron bases 7 and 8: {pos}");
        }
        for pos in [1203u64, 1204, 1206] {
            let e = effects(&tr, pos, pos);
            assert!(e.donor_region && !e.region, "donor region takes {pos}");
        }
        let fifth = effects(&tr, 1205, 1205);
        assert!(
            fifth.donor_5th_base && !fifth.donor_region && !fifth.region,
            "the 5th base takes 1205 from both"
        );
        assert!(!region(&tr, 1209, 1209), "past the splice region");
        assert!(!region(&tr, 1500, 1500), "mid-intron");
    }

    #[test]
    fn test_polypyrimidine_forward() {
        let tr = make_forward_transcript();
        // Intron: 1201-1999. Acceptor at end (1999).
        // Polypyrimidine tract: 3-17 bases from acceptor = positions 1982-1996
        // (intron_end - 17 = 1982, intron_end - 3 = 1996)
        // VEP definition: intron_end-16 to intron_end-2 (1983-1997)

        // Check boundaries
        for pos in 1980..=2000 {
            let in_ppt = polypy(&tr, pos);
            eprintln!(
                "  pos {} (dist from end = {}): ppt={}",
                pos,
                1999u64.saturating_sub(pos),
                in_ppt
            );
        }

        // Distance 17 from intron_end(1999) = 1982 = intron_end - 17
        // Distance 3 from intron_end(1999) = 1996 = intron_end - 3
        // But VEP measures from exon boundary (2000):
        //   dist 17 from exon = 2000-17 = 1983 = intron_end - 16
        //   dist 3 from exon = 2000-3 = 1997 = intron_end - 2
        assert!(
            polypy(&tr, 1983),
            "pos 1983 (dist 17 from exon) should be PPT"
        );
        assert!(
            polypy(&tr, 1997),
            "pos 1997 (dist 3 from exon) should be PPT"
        );
        assert!(
            !polypy(&tr, 1982),
            "pos 1982 (dist 18 from exon) should NOT be PPT"
        );
        assert!(
            !polypy(&tr, 1998),
            "pos 1998 (dist 2 from exon) should NOT be PPT - it's the acceptor site"
        );
    }

    fn make_reverse_transcript() -> Transcript {
        // Reverse strand: exons sorted in descending genomic order
        // Exon1 (rank 1, 5'): 2000-2300 (higher coords)
        // Exon2 (rank 2, 3'): 1000-1200 (lower coords)
        // Intron: 1201-1999
        // For reverse: donor at intron_end (1999), acceptor at intron_start (1201)
        Transcript {
            stable_id: "ENST_REV".into(),
            version: None,
            gene: Gene {
                stable_id: "ENSG_REV".into(),
                symbol: None,
                symbol_source: None,
                hgnc_id: None,
                biotype: "protein_coding".into(),
                chromosome: "chr1".into(),
                start: 1000,
                end: 2300,
                strand: Strand::Reverse,
            },
            biotype: "protein_coding".into(),
            chromosome: "chr1".into(),
            start: 1000,
            end: 2300,
            strand: Strand::Reverse,
            exons: vec![
                Exon {
                    stable_id: "E1".into(),
                    start: 2000,
                    end: 2300,
                    strand: Strand::Reverse,
                    phase: 0,
                    end_phase: 0,
                    rank: 1,
                },
                Exon {
                    stable_id: "E2".into(),
                    start: 1000,
                    end: 1200,
                    strand: Strand::Reverse,
                    phase: 0,
                    end_phase: 0,
                    rank: 2,
                },
            ],
            translation: Some(Translation {
                stable_id: "P1".into(),
                genomic_start: 1000,
                genomic_end: 2300,
                start_exon_rank: 1,
                start_exon_offset: 0,
                end_exon_rank: 2,
                end_exon_offset: 200,
            }),
            cdna_coding_start: Some(1),
            cdna_coding_end: Some(502),
            coding_region_start: Some(1000),
            coding_region_end: Some(2300),
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
    fn test_polypyrimidine_reverse() {
        let tr = make_reverse_transcript();
        // Intron: 1201-1999. On reverse strand, acceptor at intron_start (1201).
        // Polypyrimidine tract: 3-17 bases from acceptor
        // VEP: distance measured from exon boundary (1200)
        //   dist 3: 1200+3 = 1203 = intron_start + 2
        //   dist 17: 1200+17 = 1217 = intron_start + 16

        for pos in 1199..=1220 {
            let in_ppt = polypy(&tr, pos);
            eprintln!(
                "  REV pos {} (dist from intron_start=1201: {}): ppt={}",
                pos,
                pos as i64 - 1201,
                in_ppt
            );
        }

        // c.X-17 = 17 bases from exon boundary = 1200+17 = 1217 = intron_start + 16
        assert!(
            polypy(&tr, 1217),
            "pos 1217 (c.X-17, dist 16 from intron_start) should be PPT"
        );
        // c.X-3 = 3 bases from exon boundary = 1200+3 = 1203 = intron_start + 2
        assert!(
            polypy(&tr, 1203),
            "pos 1203 (c.X-3, dist 2 from intron_start) should be PPT"
        );
        // c.X-18 = 18 bases = 1218 = intron_start + 17
        assert!(!polypy(&tr, 1218), "pos 1218 (c.X-18) should NOT be PPT");
        // c.X-2 = 2 bases = 1202 = intron_start + 1 (acceptor site)
        assert!(
            !polypy(&tr, 1202),
            "pos 1202 (c.X-2, acceptor site) should NOT be PPT"
        );
    }
    /// A variant whose *second* base lands on the donor dinucleotide is a
    /// `splice_donor_variant`, not merely a `splice_region_variant`.
    ///
    /// Ensembl tests every splice site by overlap with the variant's whole span
    /// (`_intron_effects`, `BaseTranscriptVariationAllele.pm`). Reading
    /// `var_start` alone made BRCA2 `c.9256_9256+1delinsTA` - the last base of an
    /// exon and the first of the intron, changed together - LOW where VEP says
    /// HIGH. It also reported `stop_gained` for it, from a codon window built as
    /// if the two bases were adjacent in the CDS.
    #[test]
    fn a_variant_reaching_onto_the_donor_dinucleotide_is_a_donor_variant() {
        let tr = make_forward_transcript();
        // The first intron starts at 1201, so 1201-1202 is the donor.
        assert!(
            donor(&tr, 1200, 1201),
            "a change over the exon/intron junction must reach the donor"
        );
        assert!(
            !donor(&tr, 1198, 1199),
            "a change wholly inside the exon must not"
        );
        assert!(
            donor(&tr, 1202, 1210),
            "a change starting on the donor's second base still reaches it"
        );
        assert!(
            !donor(&tr, 1203, 1210),
            "a change starting past the donor must not"
        );

        // An insertion arrives as the zero-length interval `end = start - 1`, and
        // counts only where it sits inside the site rather than abutting it -
        // which is what Ensembl's `overlap` yields for that convention.
        assert!(
            donor(&tr, 1202, 1201),
            "an insertion between the two donor bases is inside the site"
        );
        assert!(
            !donor(&tr, 1201, 1200),
            "an insertion at the exon/intron boundary abuts the site, not inside it"
        );
        assert!(
            !donor(&tr, 1203, 1202),
            "an insertion just past the donor abuts it, not inside it"
        );
    }
    /// An insertion on an intron edge or on the inner edge of a donor or
    /// acceptor dinucleotide is in the splice region, even though its
    /// zero-length interval overlaps none of the four ranges that define one.
    ///
    /// This is Ensembl's own special case (`_intron_overlap`), and without it a
    /// duplication like `c.830+2dup` collects no splice term at all - MODIFIER
    /// where real VEP 115.1 reports `splice_region_variant`, LOW. Ten ClinVar
    /// 2-star pathogenic splice-site duplications land here.
    #[test]
    fn an_insertion_on_a_splice_boundary_is_in_the_splice_region() {
        let tr = make_forward_transcript();
        // First intron: 1201..1999 inclusive on this fixture.
        // An insertion arrives as `end = start - 1`.
        for (start, end, expected, what) in [
            (1201u64, 1200u64, true, "on the intron's first base"),
            (1203, 1202, true, "on the donor dinucleotide's inner edge"),
            (2000, 1999, true, "on the intron's last base"),
            (
                1998,
                1997,
                true,
                "on the acceptor dinucleotide's inner edge",
            ),
            (1600, 1599, false, "deep inside the intron"),
        ] {
            assert_eq!(
                region(&tr, start, end),
                expected,
                "insertion {what} ({start},{end})"
            );
        }
    }

    /// The bases Ensembl tests are the ones that differ, not the reference
    /// allele's span - and that is both a narrowing and a widening.
    #[test]
    fn only_the_differing_bases_are_matched_against_a_site() {
        let tr = make_forward_transcript();
        // Exon 1 ends at 1200; intron 1 is 1201..1999, so the donor is 1201-1202.
        let seq = |b: &str| Allele::Sequence(b.as_bytes().to_vec());

        // `CCT/TCTG` differs at offsets 0 and 3 only, so a change starting at
        // 1200 leaves the donor untouched even though its span covers it - and
        // is instead tested at 1203, in the donor region.
        let split = splice_effects(&tr, 1200, 1202, &seq("CCT"), &seq("TCTG"));
        assert!(!split.donor, "the matching interior is not tested");
        assert!(split.donor_region, "the padded tail past it is");

        // The same span with every base changed does reach the donor.
        let solid = splice_effects(&tr, 1200, 1202, &seq("CCT"), &seq("TTA"));
        assert!(solid.donor, "a change over the donor reaches it");

        // Padding makes every base past the reference allele differ, so a
        // two-base replacement by seven is tested over seven bases.
        let padded = splice_effects(&tr, 1198, 1199, &seq("CC"), &seq("TTAAAAA"));
        assert!(
            padded.donor,
            "the inserted tail runs into the donor and is tested there"
        );
    }

    /// An insertion on the polypyrimidine tract's first base is in the tract,
    /// even though the same insertion on a donor's first base is not in the
    /// donor.
    ///
    /// Ensembl runs the two polypyrimidine tests against the sorted pair and
    /// every other test against the unsorted one, so the zero-length insertion
    /// interval means "abuts or inside" for the tract and "strictly inside"
    /// everywhere else.
    #[test]
    fn an_insertion_abutting_the_polypyrimidine_tract_is_in_it() {
        let tr = make_forward_transcript();
        // Intron 1 is 1201..1999, so the tract is 1983..1997.
        let ins = |at: u64| {
            splice_effects(
                &tr,
                at + 1,
                at,
                &Allele::Deletion,
                &Allele::Sequence(b"T".to_vec()),
            )
        };
        assert!(
            ins(1982).polypyrimidine_tract,
            "an insertion on the tract's first base lengthens the tract"
        );
        assert!(ins(1990).polypyrimidine_tract, "and one inside it");
        assert!(
            !ins(1981).polypyrimidine_tract,
            "one a base further out does not"
        );
    }

    /// The extended splice terms are not a fallback for the essential ones.
    ///
    /// Ensembl reports `splice_donor_variant,splice_donor_5th_base_variant` for
    /// a deletion covering both, and `splice_acceptor_variant,splice_polypyrimidine_tract_variant`
    /// for one covering those. Only `splice_region_variant` is displaced.
    #[test]
    fn an_essential_site_does_not_swallow_the_extended_terms() {
        let tr = make_forward_transcript();
        let del = |s: u64, e: u64| {
            let n = (e - s + 1) as usize;
            splice_effects(
                &tr,
                s,
                e,
                &Allele::Sequence(vec![b'A'; n]),
                &Allele::Deletion,
            )
        };
        // Intron 1 is 1201..1999: donor 1201-1202, 5th base 1205, tract 1983..1997,
        // acceptor 1998-1999.
        let donor_side = del(1201, 1205);
        assert!(donor_side.donor && donor_side.donor_5th_base);
        assert!(
            !donor_side.region,
            "splice_region is the only term the others displace"
        );
        let acceptor_side = del(1995, 1999);
        assert!(acceptor_side.acceptor && acceptor_side.polypyrimidine_tract);
    }

    /// A change whose *first* differing base is in the splice region and whose
    /// last one is not reports no splice region at all.
    ///
    /// Ensembl assigns rather than or-assigns `splice_region` on each pass of
    /// its region loop, so the last differing region decides it. This reads like
    /// a bug and is not one to work around: it is the answer real VEP 115.1
    /// gives, and or-assigning instead put a `splice_region_variant` on 63 rows
    /// of a 6,600-variant ClinVar sample that VEP calls plain missense.
    #[test]
    fn the_last_differing_region_decides_the_splice_region() {
        let tr = make_forward_transcript();
        let seq = |b: &str| Allele::Sequence(b.as_bytes().to_vec());
        // Exon 2 starts at 2000, so 2000..2002 is the exonic splice region.
        // `ATAT/TAAG` differs at offsets 0, 1 and 3: 2000, 2001 and 2003.
        let trailing = splice_effects(&tr, 2000, 2003, &seq("ATAT"), &seq("TAAG"));
        assert!(
            !trailing.region,
            "the last differing region is past the splice region"
        );
        // `ATAT/TAAT` differs at offsets 0 and 1 only, both inside it.
        let leading = splice_effects(&tr, 2000, 2003, &seq("ATAT"), &seq("TAAT"));
        assert!(leading.region, "every differing region is inside it");
    }
}
