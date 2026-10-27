//! The protein notation of a change that loses the initiator, and of one that touches only the UTR
//! (issue #139).
//!
//! A transcript is built from its exons and intron in the order of the transcript and placed on
//! either strand, and a change is given where it is in the *transcript*, so that the same test
//! covers both strands. `vep_tab_intron_repeat_golden.rs` does the same at the scale of 3,001
//! changes against Ensembl VEP; the cases here are the ones worth reading.

use fastvep_annotate::{hgvsp_for_allele, HgvspRequest};
use fastvep_cache::fasta::FastaReader;
use fastvep_cache::providers::{FastaSequenceProvider, SequenceProvider};
use fastvep_consequence::ConsequencePredictor;
use fastvep_core::{Allele, GenomicPosition, Strand};
use fastvep_genome::{Exon, Gene, Transcript, Translation};

const COMP: [(u8, u8); 4] = [(b'A', b'T'), (b'C', b'G'), (b'G', b'C'), (b'T', b'A')];

fn revcomp(s: &str) -> String {
    s.bytes()
        .rev()
        .map(|b| COMP.iter().find(|(x, _)| *x == b).unwrap().1 as char)
        .collect()
}

/// A transcript on a genome of its own, and the genome to read the 3'-rule's repeats from.
struct World {
    transcript: Transcript,
    provider: FastaSequenceProvider,
    strand: Strand,
    /// the length of the genome, which is the place of the last base of the transcript's 5' end
    /// on the reverse strand
    genome_len: u64,
    /// genomic coordinate of each base of the spliced transcript, in transcript order
    place: Vec<u64>,
}

impl World {
    /// `exons` are the sequences of the exons in transcript order and `intron` is between the
    /// first two (the rest follow each other); the CDS is the stretch of the spliced sequence
    /// from `cds.0` to `cds.1`, 0-based and half open.
    fn new(exons: &[&str], intron: &str, cds: (usize, usize), strand: Strand) -> World {
        // the sequence of the transcript's own strand, introns included
        let mut stranded = String::new();
        let mut lens = Vec::new();
        for (i, e) in exons.iter().enumerate() {
            stranded.push_str(e);
            lens.push(e.len());
            if i == 0 {
                stranded.push_str(intron);
            }
        }
        let genome = match strand {
            Strand::Forward => stranded.clone(),
            Strand::Reverse => revcomp(&stranded),
        };
        let genome_len = genome.len() as u64;
        // where each base of the transcript is in the genome, 1-based
        let mut place = Vec::new();
        let mut at = 0usize; // offset into `stranded`
        for (i, n) in lens.iter().enumerate() {
            for k in 0..*n {
                let o = at + k;
                place.push(match strand {
                    Strand::Forward => o as u64 + 1,
                    Strand::Reverse => genome_len - o as u64,
                });
            }
            at += n + if i == 0 { intron.len() } else { 0 };
        }
        let mut ascending: Vec<(u64, u64)> = Vec::new();
        let mut cursor = 0u64;
        for (i, n) in lens.iter().enumerate() {
            let (lo, hi) = (cursor + 1, cursor + *n as u64);
            ascending.push(match strand {
                Strand::Forward => (lo, hi),
                Strand::Reverse => (genome_len - hi + 1, genome_len - lo + 1),
            });
            cursor = hi + if i == 0 { intron.len() as u64 } else { 0 };
        }
        ascending.sort_unstable();
        let exon_list: Vec<Exon> = ascending
            .iter()
            .enumerate()
            .map(|(i, (s, e))| Exon {
                stable_id: format!("E{i}"),
                start: *s,
                end: *e,
                strand,
                phase: -1,
                end_phase: -1,
                rank: i as u32 + 1,
            })
            .collect();
        let (start, end) = (ascending[0].0, ascending.last().unwrap().1);
        let (coding_lo, coding_hi) = (place[cds.0], place[cds.1 - 1]);
        let (coding_region_start, coding_region_end) =
            (coding_lo.min(coding_hi), coding_lo.max(coding_hi));
        let mut transcript = Transcript {
            stable_id: "ENST_WALK".into(),
            version: Some(1),
            gene: Gene {
                stable_id: "ENSG_WALK".into(),
                symbol: None,
                symbol_source: None,
                hgnc_id: None,
                biotype: "protein_coding".into(),
                chromosome: "1".into(),
                start,
                end,
                strand,
            },
            biotype: "protein_coding".into(),
            chromosome: "1".into(),
            start,
            end,
            strand,
            exons: exon_list,
            translation: Some(Translation {
                stable_id: "ENSP_WALK".into(),
                genomic_start: coding_region_start,
                genomic_end: coding_region_end,
                start_exon_rank: 1,
                start_exon_offset: cds.0 as u64,
                end_exon_rank: 1,
                end_exon_offset: cds.1 as u64,
            }),
            cdna_coding_start: Some(cds.0 as u64 + 1),
            cdna_coding_end: Some(cds.1 as u64),
            coding_region_start: Some(coding_region_start),
            coding_region_end: Some(coding_region_end),
            spliced_seq: None,
            translateable_seq: None,
            peptide: None,
            canonical: true,
            mane_select: None,
            mane_plus_clinical: None,
            tsl: None,
            appris: None,
            ccds: None,
            protein_id: Some("ENSP_WALK".into()),
            protein_version: Some(1),
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
        let fasta = format!(">1\n{genome}\n");
        let provider =
            FastaSequenceProvider::new(FastaReader::from_reader(fasta.as_bytes()).unwrap());
        transcript
            .build_sequences(|chrom, s, e| {
                provider
                    .fetch_sequence(chrom, s, e)
                    .map_err(|err| err.to_string())
            })
            .unwrap();
        World {
            transcript,
            provider,
            strand,
            genome_len,
            place,
        }
    }

    /// The change that inserts `bases` (read along the transcript) between the bases `after` and
    /// `after + 1` of the spliced transcript, 1-based, as a VCF has it: where `start = end + 1`.
    fn insertion(&self, after: usize, bases: &str) -> Change {
        let (left, right) = (self.place[after - 1], self.place[after]);
        match self.strand {
            Strand::Forward => Change {
                start: left + 1,
                end: left,
                ref_allele: Allele::Deletion,
                alt_allele: Allele::from_str(bases),
            },
            Strand::Reverse => Change {
                // the bases are read along the other strand
                start: left,
                end: right,
                ref_allele: Allele::Deletion,
                alt_allele: Allele::from_str(&revcomp(bases)),
            },
        }
    }

    /// The change that deletes the bases `first` to `last` of the spliced transcript, 1-based.
    fn deletion(&self, first: usize, last: usize) -> Change {
        let (a, b) = (self.place[first - 1], self.place[last - 1]);
        let (lo, hi) = (a.min(b), a.max(b));
        let genome = self.genome();
        Change {
            start: lo,
            end: hi,
            ref_allele: Allele::from_str(&genome[lo as usize - 1..hi as usize]),
            alt_allele: Allele::Deletion,
        }
    }

    fn genome(&self) -> String {
        let bytes = self
            .provider
            .fetch_sequence("1", 1, self.genome_len)
            .unwrap();
        String::from_utf8(bytes).unwrap()
    }

    /// The notation of the change, with the genome to walk over when `walk` is true.
    fn hgvsp(&self, change: &Change, walk: bool) -> Option<String> {
        let predictor = ConsequencePredictor::new(5000, 5000);
        let position = GenomicPosition::new("1", change.start, change.end, Strand::Forward);
        let result = predictor.predict(
            &position,
            &change.ref_allele,
            std::slice::from_ref(&change.alt_allele),
            &[&self.transcript],
            None,
        );
        let ac = &result.transcript_consequences[0].allele_consequences[0];
        hgvsp_for_allele(&HgvspRequest {
            predictor: &predictor,
            transcript: &self.transcript,
            provider: walk.then_some(&self.provider as &dyn SequenceProvider),
            chrom: "1",
            start: change.start,
            end: change.end,
            ref_allele: &change.ref_allele,
            alt_allele: &change.alt_allele,
            location: &ac.location,
            consequences: &ac.consequences,
            amino_acids: ac
                .amino_acids
                .as_ref()
                .map(|(r, a)| (r.as_str(), a.as_str())),
        })
    }
}

struct Change {
    start: u64,
    end: u64,
    ref_allele: Allele,
    alt_allele: Allele,
}

const UTR5: &str = "ACGTACGTA";
/// The same transcript with no intron in it: `ATG GCT AAG`, six bases of a repeat, `CACAC G TGG TAA`
/// and a 3' UTR.
const JOINED: &str = "ACGTACGTAATGGCTAAGACACACCACACGTGGTAAACGTAC";

fn one_exon(strand: Strand) -> World {
    World::new(&[JOINED], "", (UTR5.len(), 36), strand)
}

/// A deletion that begins in the 5' UTR and runs into the start codon has no protein position for
/// Ensembl to start from, so it writes nothing; the initiator is lost all the same and HGVS names
/// that `p.Met1?`.
#[test]
fn a_change_that_loses_the_initiator_is_met1_where_ensembl_writes_nothing() {
    for strand in [Strand::Forward, Strand::Reverse] {
        let world = one_exon(strand);
        // bases 8 to 11 of the transcript: the last two of the UTR and `AT` of the initiator
        let change = world.deletion(8, 11);
        assert_eq!(
            world.hgvsp(&change, true).as_deref(),
            Some("ENSP_WALK.1:p.Met1?"),
            "{strand:?}"
        );
        // a change in the UTR alone is not one
        let utr = world.deletion(3, 4);
        assert_eq!(world.hgvsp(&utr, true), None, "{strand:?} UTR");
    }
}

/// An insertion between the last base of the 5' UTR and the initiator that repeats the codon
/// beside it: its consequence is `5_prime_UTR_variant` and nothing else, and Ensembl writes the
/// change of the protein it would be a base along (`p.Met1_Ala2insMet`) beside it. A row with no
/// term for the coding sequence has no protein notation.
#[test]
fn a_change_that_is_only_in_the_utr_has_no_protein_notation() {
    for strand in [Strand::Forward, Strand::Reverse] {
        let world = one_exon(strand);
        let change = world.insertion(9, "ATG");
        assert_eq!(world.hgvsp(&change, true), None, "{strand:?}");
    }
}
