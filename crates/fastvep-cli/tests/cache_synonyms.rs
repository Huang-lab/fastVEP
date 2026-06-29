//! End-to-end tests for merged cache building across chromosome naming
//! schemes (issue #47): Ensembl bare names (`17`) + RefSeq accessions
//! (`NC_000017.11`) reconciled against a single FASTA via `--synonyms`.

use fastvep_cli::pipeline::run_cache_build;
use std::fs::File;
use std::io::Write;
use std::path::Path;

// Ensembl-style GFF3: bare seqid `17`.
const ENSEMBL_GFF3: &str = "##gff-version 3\n\
17\tensembl\tgene\t1\t30\t.\t+\t.\tID=gene:ENSG1;Name=ENSG1;gene_name=ENSG1\n\
17\tensembl\tmRNA\t1\t30\t.\t+\t.\tID=transcript:ENST1;Parent=gene:ENSG1;biotype=protein_coding\n\
17\tensembl\texon\t1\t30\t.\t+\t.\tID=exon:EE1;Parent=transcript:ENST1;rank=1\n\
17\tensembl\tCDS\t1\t30\t.\t+\t0\tID=CDS:EP1;Parent=transcript:ENST1\n";

// RefSeq-style GFF3: accession seqid `NC_000017.11`.
const REFSEQ_GFF3: &str = "##gff-version 3\n\
NC_000017.11\tBestRefSeq\tgene\t1\t30\t.\t+\t.\tID=gene:RSG1;Name=RSG1;gene_name=RSG1\n\
NC_000017.11\tBestRefSeq\tmRNA\t1\t30\t.\t+\t.\tID=transcript:RST1;Parent=gene:RSG1;biotype=protein_coding\n\
NC_000017.11\tBestRefSeq\texon\t1\t30\t.\t+\t.\tID=exon:RE1;Parent=transcript:RST1;rank=1\n\
NC_000017.11\tBestRefSeq\tCDS\t1\t30\t.\t+\t0\tID=CDS:RP1;Parent=transcript:RST1\n";

// FASTA uses the Ensembl bare contig name only.
const FASTA: &str = ">17\nACGTACGTACGTACGTACGTACGTACGTAC\n";

// VEP chr_synonyms.txt mapping the three equivalent names.
const SYNONYMS: &str = "17\tchr17\tNC_000017.11\n";

fn write(path: &Path, contents: &str) {
    File::create(path).unwrap().write_all(contents.as_bytes()).unwrap();
}

#[test]
fn merged_cache_with_synonyms_canonicalizes_to_fasta_names() {
    let tmp = tempfile::tempdir().unwrap();
    let ens = tmp.path().join("ensembl.gff3");
    let rs = tmp.path().join("refseq.gff3");
    let fa = tmp.path().join("ref.fa");
    let syn = tmp.path().join("chr_synonyms.txt");
    let out = tmp.path().join("combined.cache");
    write(&ens, ENSEMBL_GFF3);
    write(&rs, REFSEQ_GFF3);
    write(&fa, FASTA);
    write(&syn, SYNONYMS);

    run_cache_build(
        &[ens.to_string_lossy().into(), rs.to_string_lossy().into()],
        Some(fa.to_str().unwrap()),
        Some(syn.to_str().unwrap()),
        out.to_str().unwrap(),
        false,
    )
    .expect("merged cache build should succeed with a synonyms file");

    let transcripts = fastvep_cache::transcript_cache::load_cache(&out).unwrap();
    assert_eq!(transcripts.len(), 2, "both sources should contribute a transcript");
    for tr in &transcripts {
        // (b) every transcript canonicalized to the FASTA contig name.
        assert_eq!(&*tr.chromosome, "17", "transcript {} not canonicalized", tr.stable_id);
        assert_eq!(&*tr.gene.chromosome, "17", "gene of {} not canonicalized", tr.stable_id);
        // (c) coding sequences were built for both sources.
        assert!(
            tr.spliced_seq.is_some(),
            "sequence not built for {} — FASTA fetch must have matched",
            tr.stable_id
        );
    }
}

#[test]
fn refseq_unresolved_without_synonyms_but_build_still_succeeds() {
    let tmp = tempfile::tempdir().unwrap();
    let ens = tmp.path().join("ensembl.gff3");
    let rs = tmp.path().join("refseq.gff3");
    let fa = tmp.path().join("ref.fa");
    let out = tmp.path().join("combined.cache");
    write(&ens, ENSEMBL_GFF3);
    write(&rs, REFSEQ_GFF3);
    write(&fa, FASTA);

    // No synonyms file: the RefSeq accession has no mechanical relationship to
    // `17`, so it can't be canonicalized. The build must still succeed.
    run_cache_build(
        &[ens.to_string_lossy().into(), rs.to_string_lossy().into()],
        Some(fa.to_str().unwrap()),
        None,
        out.to_str().unwrap(),
        false,
    )
    .expect("build should not fail just because one contig is unresolved");

    let transcripts = fastvep_cache::transcript_cache::load_cache(&out).unwrap();
    let ens_tr = transcripts.iter().find(|t| &*t.stable_id == "ENST1").unwrap();
    let rs_tr = transcripts.iter().find(|t| &*t.stable_id == "RST1").unwrap();
    // Ensembl transcript resolves and gets a sequence.
    assert_eq!(&*ens_tr.chromosome, "17");
    assert!(ens_tr.spliced_seq.is_some());
    // RefSeq transcript keeps its accession and has no sequence.
    assert_eq!(&*rs_tr.chromosome, "NC_000017.11");
    assert!(rs_tr.spliced_seq.is_none());
}
