//! The location cells of a transcript row (`cDNA_position`, `CDS_position`,
//! `Protein_position`, `EXON`, `INTRON`), computed the way Ensembl VEP computes them.

use crate::transcript::Transcript;
use fastvep_core::{Ordinals, Strand, TranscriptLocation, TranscriptSpan};

impl Transcript {
    /// Where a change at genomic `start..=end` sits in this transcript, as the cells of
    /// VEP's tab output.
    ///
    /// `start` and `end` are Ensembl's own for a variant feature: the first and last base
    /// a deletion or substitution touches, and for an insertion `start = end + 1`, the two
    /// bases it sits between.
    ///
    /// **Positions.** Ensembl maps the whole span through the transcript's exons
    /// (`Mapper::map_coordinates`, in transcript order), which gives a run of coordinates
    /// and gaps, and takes the `start` of the first piece and the `end` of the last one
    /// (`BaseTranscriptVariation::cdna_start`). A piece that is a gap, which is any
    /// stretch of the span that is not exon, leaves that end undefined, and
    /// `format_coords` writes it `?`. So an end has a coordinate exactly when the base at
    /// that end, the 5'-most for the start and the 3'-most for the end, is in an exon;
    /// the CDS and the protein are the same test against the coding sequence
    /// (`TranscriptMapper::genomic2cds` makes the UTR part of a piece a gap).
    ///
    /// An insertion is mapped by `Mapper::map_insert`, which maps the two bases beside it
    /// and drops the gaps. A flank in an exon contributes its coordinate and the one
    /// across the boundary, whether or not that one exists: an insertion in the intron
    /// right after an exon's last base `601` is `601-602`. For the CDS the piece
    /// `(c + 1, c)` is a gap when `c` lies before the first coding base or `c + 1` after
    /// the last, so an insertion immediately before the initiator or after the terminator
    /// has no CDS or protein position.
    ///
    /// **Numbers.** `EXON` and `INTRON` list every exon, or intron, that the span overlaps
    /// by `overlap($vf_start, $vf_end, $feature_start, $feature_end)`, which for an
    /// insertion (`start > end`) is true only when the insertion sits between two bases of
    /// the feature. Introns are the gaps between consecutive exons in transcript order.
    ///
    /// Checked against Ensembl VEP 116.3 over 17,657 changes at every exon, CDS and
    /// transcript edge of a generated gene set (`fixtures/vep_tab_edges`).
    pub fn location(&self, start: u64, end: u64) -> TranscriptLocation {
        let n = self.exons.len();
        if n == 0 {
            return TranscriptLocation::default();
        }
        let forward = self.strand == Strand::Forward;
        let insertion = start == end + 1;

        // The exons in transcript order, without sorting a copy: the list runs one way
        // or the other, and which way is the same bargain `introns` makes (`gff.rs`
        // stores exons in transcript order). This runs once per allele per overlapping
        // transcript.
        debug_assert!(
            self.exons.windows(2).all(|p| p[0].start <= p[1].start)
                || self.exons.windows(2).all(|p| p[0].start >= p[1].start),
            "exons of {} are neither ascending nor descending",
            self.stable_id
        );
        let flip = (self.exons[0].start <= self.exons[n - 1].start) != forward;
        let at = |i: usize| &self.exons[if flip { n - 1 - i } else { i }];

        let cdna_of = |g: u64| -> Option<u64> {
            let mut before = 0u64;
            for i in 0..n {
                let e = at(i);
                if g >= e.start && g <= e.end {
                    let within = if forward { g - e.start } else { e.end - g };
                    return Some(before + within + 1);
                }
                before += e.end - e.start + 1;
            }
            None
        };

        // The ends in transcript order. For an insertion they are its two flanks, `end`
        // being the base to its left and `start` the one to its right.
        let (five, three) = match (insertion, forward) {
            (false, true) => (start, end),
            (false, false) => (end, start),
            (true, true) => (end, start),
            (true, false) => (start, end),
        };
        let (c5, c3) = (cdna_of(five), cdna_of(three));
        let cdna = if insertion {
            match (c5, c3) {
                (Some(c), _) => TranscriptSpan::new(Some(c + 1), Some(c)),
                (None, Some(c)) => TranscriptSpan::new(Some(c), Some(c - 1)),
                (None, None) => TranscriptSpan::NONE,
            }
        } else {
            TranscriptSpan::new(c5, c3)
        };

        let (cds, protein) = match (self.cdna_coding_start, self.cdna_coding_end) {
            (Some(cs), Some(ce)) if self.is_coding() => {
                let phase = self.codon_table_start_phase;
                let to_cds = |c: u64| c - cs + 1 + phase;
                let coding = |c: u64| c >= cs && c <= ce;
                let cds = match (insertion, cdna.start, cdna.end) {
                    (true, Some(s), Some(e)) if !(e < cs || s > ce) => {
                        TranscriptSpan::new(Some(to_cds(s)), Some(to_cds(e)))
                    }
                    (true, _, _) => TranscriptSpan::NONE,
                    (false, s, e) => TranscriptSpan::new(
                        s.filter(|&c| coding(c)).map(to_cds),
                        e.filter(|&c| coding(c)).map(to_cds),
                    ),
                };
                let residue = |p: u64| p.div_ceil(3);
                (
                    cds,
                    TranscriptSpan::new(cds.start.map(residue), cds.end.map(residue)),
                )
            }
            _ => (TranscriptSpan::NONE, TranscriptSpan::NONE),
        };

        // Ensembl's `overlap`: `(f1_end >= f2_start) and (f1_start <= f2_end)`
        let overlaps = |s: u64, e: u64| end >= s && start <= e;
        let exon = ordinals(
            (0..n)
                .filter(|&i| overlaps(at(i).start, at(i).end))
                .map(|i| i as u32 + 1),
            n,
        );
        let intron = ordinals(
            (0..n - 1)
                .filter(|&i| {
                    // the intron between the i-th exon and the next runs from the end of
                    // the lower one to the start of the upper one, whichever way the
                    // transcript does; coordinates are 1-based, so `start - 1` is safe
                    let (a, b) = (at(i), at(i + 1));
                    let (lower, upper) = if a.start <= b.start { (a, b) } else { (b, a) };
                    let (s, e) = (lower.end + 1, upper.start - 1);
                    // overlapping exons (seen in RefSeq) leave an empty intron
                    s <= e && overlaps(s, e)
                })
                .map(|i| i as u32 + 1),
            n - 1,
        );

        TranscriptLocation {
            cdna,
            cds,
            protein,
            exon,
            intron,
        }
    }
}

/// The first and last of the numbers `hits` yields, with the count they are out of.
fn ordinals(hits: impl Iterator<Item = u32>, total: usize) -> Option<Ordinals> {
    let (mut first, mut last) = (u32::MAX, 0u32);
    for n in hits {
        first = first.min(n);
        last = last.max(n);
    }
    (last != 0).then_some(Ordinals {
        first,
        last,
        total: total as u32,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::transcript::{Exon, Gene, Translation};

    /// A transcript of the generated gene set the Ensembl VEP 116.3 comparison runs on
    /// (`fixtures/vep_tab_world/gen_world.py`): exons as `(start, end)` ascending, and
    /// the coding span as genomic `(start, end)` when there is one.
    fn transcript(strand: Strand, exons: &[(u64, u64)], cds: Option<(u64, u64)>) -> Transcript {
        let mut t = Transcript {
            stable_id: "ENST".into(),
            version: Some(1),
            gene: Gene {
                stable_id: "ENSG".into(),
                symbol: None,
                symbol_source: None,
                hgnc_id: None,
                biotype: "protein_coding".into(),
                chromosome: "chr10".into(),
                start: exons[0].0,
                end: exons[exons.len() - 1].1,
                strand,
            },
            biotype: "protein_coding".into(),
            chromosome: "chr10".into(),
            start: exons[0].0,
            end: exons[exons.len() - 1].1,
            strand,
            exons: exons
                .iter()
                .enumerate()
                .map(|(i, &(start, end))| Exon {
                    stable_id: format!("ENSE{i}"),
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
        };
        if let Some((s, e)) = cds {
            // the first coding base is the low end on the forward strand, the high end on
            // the reverse one
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
                stable_id: "ENSP".into(),
                genomic_start: s,
                genomic_end: e,
                start_exon_rank: 1,
                start_exon_offset: 0,
                end_exon_rank: exons.len() as u32,
                end_exon_offset: 0,
            });
        }
        t
    }

    fn cell(span: TranscriptSpan) -> String {
        let mut s = String::new();
        span.write_to(&mut s);
        if s.is_empty() {
            "-".into()
        } else {
            s
        }
    }

    fn ordinal(o: Option<Ordinals>) -> String {
        let mut s = String::new();
        match o {
            Some(o) => o.write_to(&mut s),
            None => s.push('-'),
        }
        s
    }

    /// `(cDNA_position, CDS_position, Protein_position, EXON, INTRON)` as the tab file
    /// writes them.
    fn cells(t: &Transcript, start: u64, end: u64) -> [String; 5] {
        let l = t.location(start, end);
        [
            cell(l.cdna),
            cell(l.cds),
            cell(l.protein),
            ordinal(l.exon),
            ordinal(l.intron),
        ]
    }

    fn forward_coding() -> Transcript {
        // ENST00000000001 of the generated gene set: coding 1100..3251, cDNA 401..1054
        transcript(
            Strand::Forward,
            &[(700, 1300), (2000, 2200), (3000, 3600)],
            Some((1100, 3251)),
        )
    }

    fn reverse_coding() -> Transcript {
        // ENST00000000003: coding 10099..11300 (cDNA 401..1003), running right to left
        transcript(
            Strand::Reverse,
            &[(9700, 10400), (11000, 11700)],
            Some((10099, 11300)),
        )
    }

    fn reverse_noncoding() -> Transcript {
        // ENST00000000005
        transcript(Strand::Reverse, &[(2500, 2600), (3700, 3800)], None)
    }

    // Every expected value below is what Ensembl VEP 116.3 wrote for that change
    // (`chr10_POS_REF/ALT`, in `fixtures/vep_tab_edges/golden_flags.tab.gz`).

    #[test]
    fn a_deletion_that_ends_in_an_intron_has_no_coordinate_for_that_end() {
        // chr10_1291_CGCACAACCATTG/C deletes 1292..=1303, exon 1 ends at 1300
        assert_eq!(
            cells(&forward_coding(), 1292, 1303),
            ["593-?", "193-?", "65-?", "1/3", "1/2"]
        );
    }

    #[test]
    fn a_deletion_that_starts_in_an_intron_has_none_for_its_start() {
        // chr10_1990_GACATGTCTTGCG/G deletes 1991..=2002, exon 2 starts at 2000
        assert_eq!(
            cells(&forward_coding(), 1991, 2002),
            ["?-604", "?-204", "?-68", "2/3", "1/2"]
        );
    }

    #[test]
    fn the_unknown_end_is_the_transcripts_not_the_genomes() {
        // reverse strand: the genomically-low end of the span is the 3' end of the change
        // chr10_2490_TTGACCGCAGCTG/T deletes 2491..=2502; exon 2 (transcript order) is 2500..=2600
        assert_eq!(
            cells(&reverse_noncoding(), 2491, 2502),
            ["200-?", "-", "-", "2/2", "-"]
        );
        // chr10_3791_AGCTGGTACAGGT/A deletes 3792..=3803; exon 1 is 3700..=3800
        assert_eq!(
            cells(&reverse_noncoding(), 3792, 3803),
            ["?-9", "-", "-", "1/2", "-"]
        );
    }

    #[test]
    fn a_change_inside_an_exon_that_leaves_the_coding_sequence_has_an_unknown_end() {
        // chr10_1090_AGTAGCTGAGATG/A deletes 1091..=1102, the coding sequence starts at 1100
        assert_eq!(
            cells(&forward_coding(), 1091, 1102),
            ["392-403", "?-3", "?-1", "1/3", "-"]
        );
        // reverse strand: chr10_10089_CAGAAATCAGTTA/C deletes 10090..=10101, coding ends at 10099
        assert_eq!(
            cells(&reverse_coding(), 10090, 10101),
            ["1001-1012", "601-?", "201-?", "2/2", "-"]
        );
    }

    #[test]
    fn an_insertion_is_the_pair_of_positions_it_sits_between() {
        let t = forward_coding();
        // chr10_1100_A/AA, inside the initiator: between 1100 and 1101
        assert_eq!(cells(&t, 1101, 1100), ["401-402", "1-2", "1", "1/3", "-"]);
        // chr10_1300_A/AA, in the intron right after exon 1's last base
        assert_eq!(
            cells(&t, 1301, 1300),
            ["601-602", "201-202", "67-68", "-", "-"]
        );
        // chr10_1999_T/TACATGTCTT, in the intron right before exon 2's first base
        assert_eq!(
            cells(&t, 2000, 1999),
            ["601-602", "201-202", "67-68", "-", "-"]
        );
    }

    #[test]
    fn an_insertion_beside_the_coding_sequence_has_no_cds_position() {
        let t = forward_coding();
        // chr10_1099_G/GA: between the last UTR base and the initiator
        assert_eq!(cells(&t, 1100, 1099), ["400-401", "-", "-", "1/3", "-"]);
        // chr10_3251_A/AG: between the terminator's last base and the first UTR base
        assert_eq!(cells(&t, 3252, 3251), ["1054-1055", "-", "-", "3/3", "-"]);
    }

    #[test]
    fn an_insertion_on_the_reverse_strand_is_read_from_the_other_flank() {
        // chr10_10400_G/GA: after exon 2's first base in transcript order, in the intron
        assert_eq!(
            cells(&reverse_coding(), 10401, 10400),
            ["701-702", "301-302", "101", "-", "-"]
        );
        // chr10_3699_C/CAAGGC on a non-coding transcript, before exon 1's first base
        assert_eq!(
            cells(&reverse_noncoding(), 3700, 3699),
            ["101-102", "-", "-", "-", "-"]
        );
    }

    #[test]
    fn a_change_wholly_outside_the_transcript_has_no_location() {
        let t = forward_coding();
        assert_eq!(cells(&t, 100, 100), ["-", "-", "-", "-", "-"]);
    }

    #[test]
    fn the_direction_the_exon_list_runs_does_not_change_the_cells() {
        // `location` reads the list in place instead of sorting a copy, so the list may
        // run either way on either strand
        let spans = [
            (1292, 1303),
            (1991, 2002),
            (1101, 1100),
            (1301, 1300),
            (3252, 3251),
            (1291, 2001),
            (100, 100),
        ];
        let mut flipped = forward_coding();
        flipped.exons.reverse();
        for (s, e) in spans {
            assert_eq!(flipped.location(s, e), forward_coding().location(s, e));
        }
        let mut flipped = reverse_coding();
        flipped.exons.reverse();
        for (s, e) in [(10090, 10101), (10401, 10400), (10399, 11001), (9000, 9000)] {
            assert_eq!(flipped.location(s, e), reverse_coding().location(s, e));
        }
    }

    #[test]
    fn a_deletion_across_two_exons_names_both_and_the_intron_between() {
        // 1291..=2001 runs from exon 1 through the whole of intron 1 into exon 2
        let l = forward_coding().location(1291, 2001);
        assert_eq!(ordinal(l.exon), "1-2/3");
        assert_eq!(ordinal(l.intron), "1/2");
    }
}
