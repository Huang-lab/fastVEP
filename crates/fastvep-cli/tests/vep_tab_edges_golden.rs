//! The tab output for a change at the edge of an exon, of the coding sequence or of a
//! transcript, against what Ensembl VEP 116.3 writes for the same change (issue #139).
//!
//! `fixtures/vep_tab_edges` is 2,271 single-ALT records, a sample of 17,672: every position
//! within ten bases of every exon and CDS boundary of the world in `fixtures/vep_tab_world`,
//! on both strands, in coding and non-coding transcripts: SNVs, MNVs, deletions of 1 to 30
//! bases, insertions (random and duplications), length-changing delins, and deletions of a
//! whole exon or of a whole intron. `regen.sh`
//! there runs `vep` from the 116.3 image in `--gff` mode over them and keeps what it wrote.
//! Nothing in the golden was worked out here.
//!
//! Every column is compared. The cells where fastVEP is not VEP are listed in `KNOWN`, each with
//! what fastVEP writes and why, and each is asserted to still differ: a list that only excused would
//! stop meaning anything the day VEP changes. `docs/VEP_DIVERGENCE.md` (Part 3) has the reasons
//! and the measurements on the full 17,672 records.

use std::collections::BTreeMap;
use std::io::Read;
use std::path::{Path, PathBuf};

use fastvep_cli::pipeline::{run_annotate, AnnotateConfig, ColumnFlags, PickFlags};
use flate2::read::GzDecoder;
use tempfile::TempDir;

/// The columns compared, cell by cell, on every row.
const COMPARED: &[&str] = &[
    "Location",
    "cDNA_position",
    "CDS_position",
    "Protein_position",
    "EXON",
    "INTRON",
    "Consequence",
    "IMPACT",
    "Amino_acids",
    "Codons",
    "HGVSc",
    "HGVSp",
    "HGVS_OFFSET",
];

/// Why a cell is not VEP's. Each is a section of `docs/VEP_DIVERGENCE.md`, Part 3.
#[derive(Debug, Clone, Copy)]
#[allow(dead_code)] // the reason is documentation: it is read by whoever finds the row
enum Why {
    /// VEP's `--gff` reader builds the UTR of a reverse-strand transcript from the wrong place
    /// (`BaseTranscriptVariation::_generic_utr` transfers the transcript to its feature slice),
    /// so a stop counted into the 3' UTR, and `start_retained_variant` and `stop_lost`, are
    /// answers about bases that are not there. With that accessor reading `spliced_seq` as the
    /// core API's own do, VEP and fastVEP agree on every one of these cells.
    Utr,
    /// An insertion on a reverse-strand transcript that the 3'-rule moves further than it is long:
    /// `shift_feature_seqs` turns the string `seq_length - shift_length` times, which is a
    /// negative number of times and so none, and the edited CDS holds an insertion that is not
    /// the variant's. Rows that agree once that is reproduced.
    Rotation,
    /// A change that moves the frame at a stop codon, where `frameshift` answers no and the
    /// notation is built from a window that ends mid-codon: `p.Ter218=` beside `stop_lost`, and
    /// `p.Ter218delinsTyrTer` where the edited translation has a residue that is not a stop.
    StopCodon,
    /// An insertion past the 3' end of a reverse-strand transcript, written with its two
    /// positions the wrong way round.
    ReversedPair,
    /// A change that loses the initiator is `p.Met1?`. For a window of more than one residue VEP
    /// writes the residue names of the whole window sharing one position (`p.MetCys1_?2`) or the
    /// first one and a range (`p.Met1_?2`), for an insertion in front of the initiator a range
    /// that runs backwards and names the last residue of the protein (`p.Ala1_?0`), and for a
    /// deletion that begins in the 5' UTR and runs into the start codon, which has no protein start
    /// (`?-1`), nothing: none is HGVS or says that the initiator is gone.
    StartLoss,
    /// An insertion between the last base of the 5' UTR and the initiator that the 3'-rule carries
    /// onto the first codon is written `p.Met1IlefsTer?` on a row whose only term is
    /// `5_prime_UTR_variant`: the consequence is of the change where the VCF puts it and the
    /// notation of the change where the 3'-rule moved it, and the row says both. fastVEP writes
    /// no protein change for a row with no term for the coding sequence.
    UtrOnly,
    /// A change in a repeat that goes on across an intron, which the 3'-rule's walk over the genome
    /// carries out of the exon: VEP then writes no HGVSp (the walk ended in the intron, which is
    /// where an insertion at the edge of an exon goes when the intron repeats its bases), or one
    /// for the residue after the one the change is in (it turned the inserted string once for each
    /// base of the genome it passed, the intron's included). fastVEP walks the sequence of the
    /// transcript instead, and writes the notation `vep_tab_intron_repeat` shows is the same as for
    /// the same change in a transcript that has no intron there.
    IntronWalk,
    /// An edit that leaves the whole coding sequence behind it, one UTR shorter or longer
    /// (`c.-8_1del` inside a repeat, or bases inserted inside the start codon that repeat the one
    /// before): Ensembl's `_inv_start_altered` asks the codon at the old coding start and writes
    /// `start_lost` (HIGH) beside `start_retained_variant`. The edited sequence decides here, and
    /// the protein is the protein.
    CdsKept,
}

/// Cells where fastVEP and stock VEP 116.3 differ and VEP is the one that is wrong:
/// `(record, transcript, column, what fastVEP writes, why)`.
const KNOWN: &[(&str, &str, &str, &str, Why)] = &[
    (
        "chr10_10098_GTTAAGCAAGACCCAA/G",
        "ENST00000000003",
        "HGVSp",
        "ENSP00000000003.1:p.Leu197_Ter201delextTer2",
        Why::Utr,
    ),
    (
        "chr10_10099_TT/GAA",
        "ENST00000000003",
        "HGVSp",
        "ENSP00000000003.1:p.Ter201PheextTer35",
        Why::StopCodon,
    ),
    (
        "chr10_10099_TTAAG/T",
        "ENST00000000003",
        "HGVSp",
        "ENSP00000000003.1:p.Ala200AspfsTer2",
        Why::StopCodon,
    ),
    (
        "chr10_10100_T/G",
        "ENST00000000003",
        "HGVSp",
        "ENSP00000000003.1:p.Ter201SerextTer7",
        Why::Utr,
    ),
    (
        "chr10_10100_TAA/CG",
        "ENST00000000003",
        "HGVSp",
        "ENSP00000000003.1:p.Ter201AspextTer1",
        Why::Utr,
    ),
    (
        "chr10_10102_AGCA/CC",
        "ENST00000000003",
        "HGVSp",
        "ENSP00000000003.1:p.Ala200ValfsTer36",
        Why::Utr,
    ),
    (
        "chr10_10103_G/GCAAG",
        "ENST00000000003",
        "HGVSp",
        "ENSP00000000003.1:p.Ter201CysextTer36",
        Why::Utr,
    ),
    (
        "chr10_10103_GC/TAA",
        "ENST00000000003",
        "HGVSp",
        "ENSP00000000003.1:p.Ala200LeufsTer37",
        Why::Utr,
    ),
    (
        "chr10_10103_GCAAGA/G",
        "ENST00000000003",
        "HGVSp",
        "ENSP00000000003.1:p.Ala200ThrfsTer35",
        Why::Utr,
    ),
    (
        "chr10_10104_C/CA",
        "ENST00000000003",
        "HGVSp",
        "ENSP00000000003.1:p.Ala200CysfsTer37",
        Why::Utr,
    ),
    (
        "chr10_10106_AGAC/CA",
        "ENST00000000003",
        "HGVSp",
        "ENSP00000000003.1:p.Gly198ValfsTer38",
        Why::Utr,
    ),
    (
        "chr10_10108_A/AA",
        "ENST00000000003",
        "HGVSp",
        "ENSP00000000003.1:p.Leu199SerfsTer38",
        Why::Utr,
    ),
    (
        "chr10_10108_ACC/A",
        "ENST00000000003",
        "HGVSp",
        "ENSP00000000003.1:p.Gly198SerfsTer38",
        Why::Utr,
    ),
    (
        "chr10_10109_C/CCCAA",
        "ENST00000000003",
        "HGVSp",
        "ENSP00000000003.1:p.Gly198ValfsTer40",
        Why::Utr,
    ),
    (
        "chr10_10395_T/TATCGGTTAT",
        "ENST00000000003",
        "HGVSp",
        "ENSP00000000003.1:p.Ile103_Asp105dup",
        Why::Rotation,
    ),
    (
        "chr10_10395_T/TTAT",
        "ENST00000000003",
        "HGVSp",
        "ENSP00000000003.1:p.Ile103dup",
        Why::Rotation,
    ),
    (
        "chr10_10399_A/ATGGAA",
        "ENST00000000003",
        "HGVSp",
        "ENSP00000000003.1:p.Thr104ProfsTer2",
        Why::Rotation,
    ),
    (
        "chr10_10400_G/GGAAG",
        "ENST00000000003",
        "HGVSp",
        "ENSP00000000003.1:p.Ile103PhefsTer5",
        Why::Rotation,
    ),
    (
        "chr10_11294_T/TAT",
        "ENST00000000003",
        "HGVSp",
        "ENSP00000000003.1:p.Ile4Ter",
        Why::Rotation,
    ),
    (
        "chr10_11295_G/GATTATG",
        "ENST00000000003",
        "HGVSp",
        "ENSP00000000003.1:p.Ile3_Ile4dup",
        Why::Rotation,
    ),
    (
        "chr10_11299_A/ATGCACA",
        "ENST00000000003",
        "HGVSp",
        "ENSP00000000003.1:p.Met1_Cys2dup",
        Why::Rotation,
    ),
    (
        "chr10_3248_ATG/A",
        "ENST00000000001",
        "HGVSp",
        "ENSP00000000001.1:p.Ter218ArgextTer32",
        Why::StopCodon,
    ),
    (
        "chr10_3248_ATG/A",
        "ENST00000000002",
        "HGVSp",
        "ENSP00000000002.1:p.Ter151ArgextTer32",
        Why::StopCodon,
    ),
    (
        "chr10_3249_T/TAT",
        "ENST00000000001",
        "HGVSp",
        "ENSP00000000001.1:p.Ter218TyrextTer24",
        Why::StopCodon,
    ),
    (
        "chr10_3249_T/TAT",
        "ENST00000000002",
        "HGVSp",
        "ENSP00000000002.1:p.Ter151TyrextTer24",
        Why::StopCodon,
    ),
    (
        "chr10_3249_T/TT",
        "ENST00000000001",
        "HGVSp",
        "ENSP00000000001.1:p.Ter218LeuextTer33",
        Why::StopCodon,
    ),
    (
        "chr10_3249_T/TT",
        "ENST00000000002",
        "HGVSp",
        "ENSP00000000002.1:p.Ter151LeuextTer33",
        Why::StopCodon,
    ),
    (
        "chr10_3249_TGA/CG",
        "ENST00000000001",
        "HGVSp",
        "ENSP00000000001.1:p.Ter218ArgextTer23",
        Why::StopCodon,
    ),
    (
        "chr10_3249_TGA/CG",
        "ENST00000000002",
        "HGVSp",
        "ENSP00000000002.1:p.Ter151ArgextTer23",
        Why::StopCodon,
    ),
    (
        "chr10_3250_G/GCATG",
        "ENST00000000001",
        "HGVSp",
        "ENSP00000000001.1:p.Ter218CysextTer34",
        Why::StopCodon,
    ),
    (
        "chr10_3250_G/GCATG",
        "ENST00000000002",
        "HGVSp",
        "ENSP00000000002.1:p.Ter151CysextTer34",
        Why::StopCodon,
    ),
    (
        "chr10_3250_G/GTTTA",
        "ENST00000000001",
        "HGVSp",
        "ENSP00000000001.1:p.Ter218CysextTer34",
        Why::StopCodon,
    ),
    (
        "chr10_3250_G/GTTTA",
        "ENST00000000002",
        "HGVSp",
        "ENSP00000000002.1:p.Ter151CysextTer34",
        Why::StopCodon,
    ),
    (
        "chr10_3600_ACCT/TA",
        "ENST00000000001",
        "HGVSc",
        "ENST00000000001.1:c.*348_*349insT",
        Why::ReversedPair,
    ),
    (
        "chr10_3600_ACCT/TA",
        "ENST00000000002",
        "HGVSc",
        "ENST00000000002.1:c.*348_*349insT",
        Why::ReversedPair,
    ),
    (
        "chr10_11299_ATCTGCTCTG/A",
        "ENST00000000003",
        "Consequence",
        "start_retained_variant,5_prime_UTR_variant",
        Why::CdsKept,
    ),
    (
        "chr10_11299_ATCTGCTCTG/A",
        "ENST00000000003",
        "IMPACT",
        "LOW",
        Why::CdsKept,
    ),
    (
        "chr10_1099_GATGCACGTACGTGTCGCTACGAGGGGATGT/G",
        "ENST00000000001",
        "HGVSp",
        "ENSP00000000001.1:p.Met1?",
        Why::StartLoss,
    ),
    (
        "chr10_1099_GATGCACGTACGTGTCGCTACGAGGGGATGT/G",
        "ENST00000000002",
        "HGVSp",
        "ENSP00000000002.1:p.Met1?",
        Why::StartLoss,
    ),
    (
        "chr10_1100_ATGCA/A",
        "ENST00000000001",
        "HGVSp",
        "ENSP00000000001.1:p.Met1?",
        Why::StartLoss,
    ),
    (
        "chr10_1100_ATGCA/A",
        "ENST00000000002",
        "HGVSp",
        "ENSP00000000002.1:p.Met1?",
        Why::StartLoss,
    ),
    (
        "chr10_1101_TGCACG/T",
        "ENST00000000001",
        "HGVSp",
        "ENSP00000000001.1:p.Met1?",
        Why::StartLoss,
    ),
    (
        "chr10_1101_TGCACG/T",
        "ENST00000000002",
        "HGVSp",
        "ENSP00000000002.1:p.Met1?",
        Why::StartLoss,
    ),
    (
        "chr10_11293_ATGCAC/A",
        "ENST00000000003",
        "HGVSp",
        "ENSP00000000003.1:p.Met1?",
        Why::StartLoss,
    ),
    (
        "chr10_11294_TGCACAT/T",
        "ENST00000000003",
        "HGVSp",
        "ENSP00000000003.1:p.Met1?",
        Why::StartLoss,
    ),
    (
        "chr10_11296_CACA/AG",
        "ENST00000000003",
        "HGVSp",
        "ENSP00000000003.1:p.Met1?",
        Why::StartLoss,
    ),
    (
        "chr10_11297_AC/TCG",
        "ENST00000000003",
        "HGVSp",
        "ENSP00000000003.1:p.Met1?",
        Why::StartLoss,
    ),
    (
        "chr10_11298_C/CATCTGCTCT",
        "ENST00000000003",
        "HGVSp",
        "-",
        Why::StartLoss,
    ),
    (
        "chr10_1091_GTAGCTGAGA/G",
        "ENST00000000001",
        "HGVSp",
        "ENSP00000000001.1:p.Met1?",
        Why::StartLoss,
    ),
    (
        "chr10_1091_GTAGCTGAGA/G",
        "ENST00000000002",
        "HGVSp",
        "ENSP00000000002.1:p.Met1?",
        Why::StartLoss,
    ),
    (
        "chr10_1092_TAGCTGAGATGCA/T",
        "ENST00000000001",
        "HGVSp",
        "ENSP00000000001.1:p.Met1?",
        Why::StartLoss,
    ),
    (
        "chr10_1092_TAGCTGAGATGCA/T",
        "ENST00000000002",
        "HGVSp",
        "ENSP00000000002.1:p.Met1?",
        Why::StartLoss,
    ),
    (
        "chr10_1093_AGCTGAGATGCACGTACGTGT/A",
        "ENST00000000001",
        "HGVSp",
        "ENSP00000000001.1:p.Met1?",
        Why::StartLoss,
    ),
    (
        "chr10_1093_AGCTGAGATGCACGTACGTGT/A",
        "ENST00000000002",
        "HGVSp",
        "ENSP00000000002.1:p.Met1?",
        Why::StartLoss,
    ),
    (
        "chr10_1094_GCTGAGATGCACGTACGTGTCGCTACGAGGG/G",
        "ENST00000000001",
        "HGVSp",
        "ENSP00000000001.1:p.Met1?",
        Why::StartLoss,
    ),
    (
        "chr10_1094_GCTGAGATGCACGTACGTGTCGCTACGAGGG/G",
        "ENST00000000002",
        "HGVSp",
        "ENSP00000000002.1:p.Met1?",
        Why::StartLoss,
    ),
    (
        "chr10_1096_TGAGATG/T",
        "ENST00000000001",
        "HGVSp",
        "ENSP00000000001.1:p.Met1?",
        Why::StartLoss,
    ),
    (
        "chr10_1096_TGAGATG/T",
        "ENST00000000002",
        "HGVSp",
        "ENSP00000000002.1:p.Met1?",
        Why::StartLoss,
    ),
    (
        "chr10_1097_GAGATGCACG/G",
        "ENST00000000001",
        "HGVSp",
        "ENSP00000000001.1:p.Met1?",
        Why::StartLoss,
    ),
    (
        "chr10_1097_GAGATGCACG/G",
        "ENST00000000002",
        "HGVSp",
        "ENSP00000000002.1:p.Met1?",
        Why::StartLoss,
    ),
    (
        "chr10_1098_AGA/TC",
        "ENST00000000001",
        "HGVSp",
        "ENSP00000000001.1:p.Met1?",
        Why::StartLoss,
    ),
    (
        "chr10_1098_AGA/TC",
        "ENST00000000002",
        "HGVSp",
        "ENSP00000000002.1:p.Met1?",
        Why::StartLoss,
    ),
    (
        "chr10_1098_AGATGCACGTACGTGT/A",
        "ENST00000000001",
        "HGVSp",
        "ENSP00000000001.1:p.Met1?",
        Why::StartLoss,
    ),
    (
        "chr10_1098_AGATGCACGTACGTGT/A",
        "ENST00000000002",
        "HGVSp",
        "ENSP00000000002.1:p.Met1?",
        Why::StartLoss,
    ),
    (
        "chr10_1099_GA/CT",
        "ENST00000000001",
        "HGVSp",
        "ENSP00000000001.1:p.Met1?",
        Why::StartLoss,
    ),
    (
        "chr10_1099_GA/CT",
        "ENST00000000002",
        "HGVSp",
        "ENSP00000000002.1:p.Met1?",
        Why::StartLoss,
    ),
    (
        "chr10_11291_TTATGCACATCTG/T",
        "ENST00000000003",
        "HGVSp",
        "ENSP00000000003.1:p.Met1?",
        Why::StartLoss,
    ),
    (
        "chr10_11292_TATGCACATCTGCTCTGTATAGTATCCCCTA/T",
        "ENST00000000003",
        "HGVSp",
        "ENSP00000000003.1:p.Met1?",
        Why::StartLoss,
    ),
    (
        "chr10_11295_GCACATCTGC/G",
        "ENST00000000003",
        "HGVSp",
        "ENSP00000000003.1:p.Met1?",
        Why::StartLoss,
    ),
    (
        "chr10_11296_CACATCTGCTCTGTATAGTAT/C",
        "ENST00000000003",
        "HGVSp",
        "ENSP00000000003.1:p.Met1?",
        Why::StartLoss,
    ),
    (
        "chr10_11297_ACATCT/A",
        "ENST00000000003",
        "HGVSp",
        "ENSP00000000003.1:p.Met1?",
        Why::StartLoss,
    ),
    (
        "chr10_11298_CATCTGC/C",
        "ENST00000000003",
        "HGVSp",
        "ENSP00000000003.1:p.Met1?",
        Why::StartLoss,
    ),
    (
        "chr10_11299_ATC/CA",
        "ENST00000000003",
        "HGVSp",
        "ENSP00000000003.1:p.Met1?",
        Why::StartLoss,
    ),
    (
        "chr10_1099_G/GAT",
        "ENST00000000001",
        "HGVSp",
        "-",
        Why::UtrOnly,
    ),
    (
        "chr10_1099_G/GAT",
        "ENST00000000002",
        "HGVSp",
        "-",
        Why::UtrOnly,
    ),
    (
        "chr10_11300_T/TAT",
        "ENST00000000003",
        "HGVSp",
        "-",
        Why::UtrOnly,
    ),
    (
        "chr10_1299_C/CATTGC",
        "ENST00000000001",
        "HGVSp",
        "ENSP00000000001.1:p.Ala68LeufsTer7",
        Why::IntronWalk,
    ),
    (
        "chr10_1299_C/CATTGC",
        "ENST00000000002",
        "HGVSp",
        "ENSP00000000002.1:p.Asp68LeufsTer6",
        Why::IntronWalk,
    ),
    (
        "chr10_1300_A/ATTGCGA",
        "ENST00000000001",
        "HGVSp",
        "ENSP00000000001.1:p.Pro67_Ala68insLeuArg",
        Why::IntronWalk,
    ),
    (
        "chr10_1300_A/ATTGCGA",
        "ENST00000000002",
        "HGVSp",
        "ENSP00000000002.1:p.Pro67_Asp68insLeuArg",
        Why::IntronWalk,
    ),
    (
        "chr10_2195_A/ATACAC",
        "ENST00000000001",
        "HGVSp",
        "ENSP00000000001.1:p.Asp135TyrfsTer6",
        Why::IntronWalk,
    ),
    (
        "chr10_2196_T/TACACTC",
        "ENST00000000001",
        "HGVSp",
        "ENSP00000000001.1:p.His134_Asp135insSerHis",
        Why::IntronWalk,
    ),
    (
        "chr10_2198_C/CACTCT",
        "ENST00000000001",
        "HGVSp",
        "ENSP00000000001.1:p.Asp135SerfsTer6",
        Why::IntronWalk,
    ),
    (
        "chr10_2200_C/CTACAC",
        "ENST00000000001",
        "HGVSp",
        "ENSP00000000001.1:p.Asp135TyrfsTer6",
        Why::IntronWalk,
    ),
    (
        "chr10_2200_C/CTG",
        "ENST00000000001",
        "HGVSp",
        "ENSP00000000001.1:p.Asp135TrpfsTer5",
        Why::IntronWalk,
    ),
    (
        "chr10_10999_C/CC",
        "ENST00000000003",
        "HGVSp",
        "ENSP00000000003.1:p.Pro101ArgfsTer6",
        Why::IntronWalk,
    ),
    (
        "chr10_11000_G/GCCATCG",
        "ENST00000000003",
        "HGVSp",
        "ENSP00000000003.1:p.Leu100_Pro101insArgTrp",
        Why::IntronWalk,
    ),
    (
        "chr10_11002_A/AAGCGTA",
        "ENST00000000003",
        "HGVSp",
        "ENSP00000000003.1:p.Leu100_Pro101insArgLeu",
        Why::IntronWalk,
    ),
    (
        "chr10_11298_C/CATCTGCTCT",
        "ENST00000000003",
        "Consequence",
        "start_retained_variant",
        Why::CdsKept,
    ),
    (
        "chr10_11298_C/CATCTGCTCT",
        "ENST00000000003",
        "IMPACT",
        "LOW",
        Why::CdsKept,
    ),
];

fn here(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("fixtures/vep_tab_edges")
        .join(name)
}

fn world(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("fixtures/vep_tab_world")
        .join(name)
}

type Row = BTreeMap<String, String>;
/// Rows keyed by what names them: the record, the allele and the transcript.
type Table = BTreeMap<(String, String, String), Vec<Row>>;

fn table(text: &str) -> Table {
    let mut lines = text.lines().filter(|l| !l.starts_with("##"));
    let header: Vec<String> = lines
        .next()
        .expect("a header")
        .trim_start_matches('#')
        .split('\t')
        .map(String::from)
        .collect();
    let mut rows: Table = BTreeMap::new();
    for line in lines {
        let cells: Vec<&str> = line.split('\t').collect();
        assert_eq!(cells.len(), header.len(), "ragged row {line:?}");
        let row: Row = header
            .iter()
            .cloned()
            .zip(cells.into_iter().map(String::from))
            .collect();
        let key = (
            row["Uploaded_variation"].clone(),
            row["Allele"].clone(),
            row["Feature"].clone(),
        );
        rows.entry(key).or_default().push(row);
    }
    rows
}

fn golden() -> Table {
    let mut text = String::new();
    GzDecoder::new(std::fs::File::open(here("golden_flags.tab.gz")).unwrap())
        .read_to_string(&mut text)
        .unwrap();
    table(&text)
}

fn ours() -> Table {
    let dir = TempDir::new().unwrap();
    // Copied, because the transcript cache a run writes goes beside its GFF3 and the index
    // of the FASTA beside it, and neither belongs in the source tree.
    let (gff3, fasta) = (dir.path().join("genes.gff3"), dir.path().join("genome.fa"));
    std::fs::copy(world("genes.gff3"), &gff3).unwrap();
    std::fs::copy(world("genome.fa"), &fasta).unwrap();
    let out = dir.path().join("out.tab");
    run_annotate(AnnotateConfig {
        input: here("in.vcf").to_string_lossy().into(),
        output: out.to_string_lossy().into(),
        gff3: vec![gff3.to_string_lossy().into()],
        fasta: Some(fasta.to_string_lossy().into()),
        output_format: "tab".into(),
        pick: PickFlags::default(),
        hgvs: true,
        columns: ColumnFlags {
            show_ref_allele: true,
            variant_class: true,
            symbol: true,
            biotype: true,
            canonical: true,
            numbers: true,
            hgvsg: true,
            ..Default::default()
        },
        distance: 5000,
        cache_dir: None,
        transcript_cache: None,
        sa_dir: None,
        sa_only: false,
        acmg: false,
        acmg_config: None,
        pick_order: None,
        functional_evidence: None,
        proband: None,
        mother: None,
        father: None,
        gene_list: None,
        explicit_alleles: false,
        qc_rules: None,
        show_progress: false,
        compat: Default::default(),
        exclude: Default::default(),
    })
    .expect("annotation should succeed");
    table(&std::fs::read_to_string(&out).unwrap())
}

/// Every row VEP wrote, fastVEP wrote, and no other.
#[test]
fn both_tools_write_the_same_rows() {
    let (vep, fv) = (golden(), ours());
    let missing: Vec<_> = vep
        .keys()
        .filter(|k| !fv.contains_key(*k))
        .take(5)
        .collect();
    let extra: Vec<_> = fv
        .keys()
        .filter(|k| !vep.contains_key(*k))
        .take(5)
        .collect();
    assert!(missing.is_empty(), "rows only VEP writes: {missing:?}");
    assert!(extra.is_empty(), "rows only fastVEP writes: {extra:?}");
}

/// The cells of the compared columns, row by row.
#[test]
fn the_compared_columns_are_vep_s_for_every_row() {
    let (vep, fv) = (golden(), ours());
    let mut wrong = Vec::new();
    let mut compared = 0usize;
    let mut seen = std::collections::BTreeSet::new();
    for (key, vep_rows) in &vep {
        let fv_rows = &fv[key];
        assert_eq!(vep_rows.len(), fv_rows.len(), "{key:?}");
        for (v, f) in vep_rows.iter().zip(fv_rows) {
            for col in COMPARED {
                compared += 1;
                let known = KNOWN
                    .iter()
                    .find(|(rec, tr, c, _, _)| *rec == key.0 && *tr == key.2 && c == col);
                if let Some((rec, tr, c, ours, _)) = known {
                    assert_ne!(v[*col], f[*col], "{rec} {tr} {c}: VEP now agrees");
                    assert_eq!(&f[*col], ours, "{rec} {tr} {c}");
                    seen.insert((*rec, *tr, *c));
                } else if v[*col] != f[*col] {
                    wrong.push(format!(
                        "{} {} {col}: VEP {:?}, fastVEP {:?}",
                        key.0, key.2, v[*col], f[*col]
                    ));
                }
            }
        }
    }
    assert_eq!(
        seen.len(),
        KNOWN.len(),
        "a row `KNOWN` names is not in the golden"
    );
    assert!(compared > 50_000, "only {compared} cells compared");
    assert!(
        wrong.is_empty(),
        "{} of {compared} cells differ from VEP, e.g.\n{}",
        wrong.len(),
        wrong[..wrong.len().min(12)].join("\n")
    );
}

/// The fixture exercises what it is for: both kinds of unknown end, an insertion beside an
/// exon, one beside the coding sequence, and a change over two exons.
#[test]
fn the_fixture_covers_the_shapes_the_columns_are_about() {
    let vep = golden();
    let cells = |col: &str| -> Vec<&str> {
        vep.values()
            .flatten()
            .map(|r| r[col].as_str())
            .collect::<Vec<_>>()
    };
    let cdna = cells("cDNA_position");
    assert!(cdna.iter().any(|c| c.starts_with("?-")), "no `?-N`");
    assert!(cdna.iter().any(|c| c.ends_with("-?")), "no `N-?`");
    assert!(
        cdna.iter().any(|c| {
            c.split_once('-').is_some_and(|(a, b)| {
                a.parse::<u64>()
                    .is_ok_and(|a| b.parse::<u64>() == Ok(a + 1))
            })
        }),
        "no insertion pair"
    );
    assert!(
        cells("CDS_position").contains(&"-") && cells("EXON").contains(&"-"),
        "no row without a CDS position or an exon"
    );
    assert!(
        cells("EXON")
            .iter()
            .any(|c| c.contains('-') && c.contains('/') && !c.starts_with('-')),
        "no change over several exons"
    );
    assert!(cells("INTRON").iter().any(|c| *c != "-"), "no INTRON");
}
