//! End-to-end regression tests for GitHub issue #68: the vertebrate
//! mitochondrial codon table (NCBI translation table 2) and circular-genome
//! coordinate wrapping existed in `fastvep_genome::mitochondrial` but were
//! never wired into the production consequence-prediction / HGVS-protein
//! pipeline. These tests build a minimal MT transcript + FASTA and exercise
//! the exact production call chain (`Transcript::build_sequences` ->
//! `fastvep_cache` sequence provider -> `ConsequencePredictor::predict` ->
//! `fastvep_hgvs::hgvsp`) rather than calling the codon table directly, so a
//! regression in any one of those wiring points fails the test.

use fastvep_cache::fasta::FastaReader;
use fastvep_cache::providers::{FastaSequenceProvider, SequenceProvider};
use fastvep_consequence::ConsequencePredictor;
use fastvep_core::{Allele, Consequence, GenomicPosition, Strand};
use fastvep_genome::{Exon, Gene, Transcript, Translation, MT_LENGTH};

/// Build a single-exon, fully-coding MT transcript whose CDS is exactly the
/// concatenation of `exon_start..=exon_end` (using the "extended", possibly
/// origin-wrapping, genomic numbering described on [`fetch_circular`] in
/// fastvep-cache: coordinates may run past `MT_LENGTH` to represent a region
/// that physically continues from position 1).
fn mt_transcript(exon_start: u64, exon_end: u64, cds_len: u64) -> Transcript {
    Transcript {
        stable_id: "ENST_MT_TEST".into(),
        version: None,
        gene: Gene {
            stable_id: "ENSG_MT_TEST".into(),
            symbol: Some("MT-TEST".into()),
            symbol_source: None,
            hgnc_id: None,
            biotype: "protein_coding".into(),
            chromosome: "MT".into(),
            start: exon_start,
            end: exon_end,
            strand: Strand::Forward,
        },
        biotype: "protein_coding".into(),
        chromosome: "MT".into(),
        start: exon_start,
        end: exon_end,
        strand: Strand::Forward,
        exons: vec![Exon {
            stable_id: "ENSE_MT_TEST".into(),
            start: exon_start,
            end: exon_end,
            strand: Strand::Forward,
            phase: -1,
            end_phase: -1,
            rank: 1,
        }],
        translation: Some(Translation {
            stable_id: "ENSP_MT_TEST".into(),
            genomic_start: exon_start,
            genomic_end: exon_end,
            start_exon_rank: 1,
            start_exon_offset: 0,
            end_exon_rank: 1,
            end_exon_offset: cds_len - 1,
        }),
        cdna_coding_start: Some(1),
        cdna_coding_end: Some(cds_len),
        coding_region_start: Some(exon_start),
        coding_region_end: Some(exon_end),
        spliced_seq: None,
        translateable_seq: None,
        peptide: None,
        canonical: true,
        mane_select: None,
        mane_plus_clinical: None,
        tsl: None,
        appris: None,
        ccds: None,
        protein_id: Some("ENSP_MT_TEST".into()),
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

/// Build sequences on `transcript` using a `FastaSequenceProvider` over
/// `fasta_text`, exactly the way `fastvep-annotate`/`fastvep-cli` do in
/// production (so an MT-origin-spanning fetch goes through the real
/// `fetch_circular` wraparound logic in `fastvep-cache`, not a test double).
fn build_sequences_from_fasta(transcript: &mut Transcript, fasta_text: &str) {
    let reader = FastaReader::from_reader(fasta_text.as_bytes()).expect("valid test FASTA");
    let provider = FastaSequenceProvider::new(reader);
    transcript
        .build_sequences(|chrom, start, end| {
            provider
                .fetch_sequence(chrom, start, end)
                .map_err(|e| e.to_string())
        })
        .expect("build_sequences should succeed for a well-formed MT test transcript");
}

/// TGA is a stop codon in the standard nuclear code but Trp in the
/// vertebrate mitochondrial code. A Trp(TGG)->TGA substitution on an MT
/// transcript must therefore be `synonymous_variant`, not `stop_gained`
/// (which is what fastVEP predicted before the vertebrate mitochondrial
/// table was wired into `ConsequencePredictor`/`Transcript::build_sequences`).
#[test]
fn mt_tga_reads_as_tryptophan_not_stop() {
    // CDS (positions 1-12 on a 12bp toy "MT" contig): ATG TGG AAA TAA
    // (Met, Trp, Lys, Stop) -- codon 2 (protein position 2) is the target.
    let fasta = ">MT\nATGTGGAAATAA\n";
    let mut transcript = mt_transcript(1, 12, 12);
    build_sequences_from_fasta(&mut transcript, fasta);
    assert_eq!(transcript.translateable_seq.as_deref(), Some("ATGTGGAAATAA"));
    // Peptide translation (Transcript::build_sequences) must also use the
    // mitochondrial table: MTAK* would be wrong; the correct read is M W K *.
    assert_eq!(transcript.peptide.as_deref(), Some("MWK*"));

    // Variant: position 6 (last base of codon 2), G -> A: TGG -> TGA.
    let position = GenomicPosition::new("MT", 6, 6, Strand::Forward);
    let ref_allele = Allele::Sequence(b"G".to_vec());
    let alt_allele = Allele::Sequence(b"A".to_vec());

    let predictor = ConsequencePredictor::new(5000, 5000);
    let result = predictor.predict(&position, &ref_allele, &[alt_allele], &[&transcript], None);

    let tc = &result.transcript_consequences[0];
    let ac = &tc.allele_consequences[0];

    assert_eq!(
        ac.consequences,
        vec![Consequence::SynonymousVariant],
        "TGA must read as Trp (synonymous with Trp->Trp), not stop_gained, on an MT transcript"
    );
    assert_eq!(ac.protein_start, Some(2));
    let (ref_aa, alt_aa) = ac.amino_acids.clone().expect("amino acids computed");
    assert_eq!((ref_aa.as_str(), alt_aa.as_str()), ("W", "W"));

    let hgvsp = fastvep_hgvs::hgvsp(
        "ENSP_MT_TEST",
        ac.protein_start.unwrap(),
        ref_aa.as_bytes()[0],
        alt_aa.as_bytes()[0],
        false,
    );
    assert_eq!(hgvsp, Some("ENSP_MT_TEST:p.Trp2=".to_string()));
}

/// An MT allele spanning the circular origin (human rCRS position
/// 16569 -> 1) must be stitched into a single correct codon instead of being
/// silently truncated or erroring. This variant also exercises the AGA
/// codon: Arg in the standard code, a stop in the vertebrate mitochondrial
/// code, so getting the wrong table *or* the wrong (unwrapped) sequence both
/// independently corrupt the call -- only the fully-wired path gets both
/// pieces right and reports stop_gained.
#[test]
fn mt_origin_spanning_allele_produces_stop_gained() {
    assert_eq!(MT_LENGTH, 16569);

    // Build a 16,569 base circular "MT" contig. Only the bases this test
    // actually reads are meaningful; everything else is filler ('N').
    let mut seq = vec![b'N'; MT_LENGTH as usize];
    // Codon 1 (protein pos 1, start codon), genomic 16566-16568: "ATG".
    seq[16565] = b'A'; // pos 16566
    seq[16566] = b'T'; // pos 16567
    seq[16567] = b'G'; // pos 16568
    // Codon 2 (protein pos 2), genomic 16569, then wraps to 1, 2: "CGA".
    seq[16568] = b'C'; // pos 16569 (last base of the contig)
    seq[0] = b'G'; // pos 1 (wrapped)
    seq[1] = b'A'; // pos 2 (wrapped)
    // Codon 3 (protein pos 3), wrapped genomic 3, 4, 5: "AAA" (Lys).
    seq[2] = b'A';
    seq[3] = b'A';
    seq[4] = b'A';
    // Codon 4 (protein pos 4), wrapped genomic 6, 7, 8: "TAA" (Stop, both tables).
    seq[5] = b'T';
    seq[6] = b'A';
    seq[7] = b'A';
    let fasta = format!(">MT\n{}\n", String::from_utf8(seq).unwrap());

    // Exon spans genomic 16566 .. 16577 in "extended" numbering (12 bases:
    // 16566-16569 physically, then wrapping to 1-8) -- matching how
    // fastvep-io computes an unwrapped `end` for an MT allele/feature whose
    // footprint runs past the last physical base.
    let mut transcript = mt_transcript(16566, 16577, 12);
    build_sequences_from_fasta(&mut transcript, &fasta);
    assert_eq!(
        transcript.translateable_seq.as_deref(),
        Some("ATGCGAAAATAA"),
        "codon spanning the origin (16569->1->2 = CGA) must be assembled correctly, not truncated"
    );
    assert_eq!(transcript.peptide.as_deref(), Some("MRK*"));

    // Variant: genomic 16569-16570 (extended; 16570 wraps to physical
    // position 1), REF "CG" -> ALT "AG": codon 2 CGA -> AGA.
    let position = GenomicPosition::new("MT", 16569, 16570, Strand::Forward);
    let ref_allele = Allele::Sequence(b"CG".to_vec());
    let alt_allele = Allele::Sequence(b"AG".to_vec());

    let predictor = ConsequencePredictor::new(5000, 5000);
    let result = predictor.predict(&position, &ref_allele, &[alt_allele], &[&transcript], None);

    let tc = &result.transcript_consequences[0];
    let ac = &tc.allele_consequences[0];

    assert_eq!(
        ac.consequences,
        vec![Consequence::StopGained],
        "AGA is a stop under the vertebrate mitochondrial code -- an allele \
         spanning the MT origin that creates it must be classified as \
         stop_gained (a standard-table read would call this synonymous \
         Arg->Arg and silently miss a real stop gain)"
    );
    assert_eq!(ac.protein_start, Some(2));
    let (ref_aa, alt_aa) = ac.amino_acids.clone().expect("amino acids computed");
    assert_eq!((ref_aa.as_str(), alt_aa.as_str()), ("R", "*"));

    let hgvsp = fastvep_hgvs::hgvsp(
        "ENSP_MT_TEST",
        ac.protein_start.unwrap(),
        ref_aa.as_bytes()[0],
        alt_aa.as_bytes()[0],
        false,
    );
    assert_eq!(hgvsp, Some("ENSP_MT_TEST:p.Arg2Ter".to_string()));
}
