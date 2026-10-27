//! A transcript cache that will not load must stop the run, not quietly change
//! the answer.
//!
//! `save_cache` used to create the destination and stream into it, so a second
//! annotation started while the first was still writing read a prefix: the
//! magic header and zstd frame header are written first, so it passed every
//! check and then failed inside deserialization. The CLI warned, fell back, and
//! with no GFF3 left to rebuild from it annotated every variant as intergenic at
//! exit code 0. On validation/human/vep_example_GRCh38.vcf that produced 173
//! well-formed intergenic calls in place of real consequences.

use std::fs::File;
use std::io::Write;
use std::path::{Path, PathBuf};

use fastvep_cli::pipeline::{run_annotate, AnnotateConfig, PickFlags};
use tempfile::TempDir;

const CDS: u64 = 1_001;

fn write_gff3(dir: &Path) -> PathBuf {
    let path = dir.join("ann.gff3");
    let end = CDS + 11;
    std::fs::write(
        &path,
        format!(
            "1\ttest\tgene\t{CDS}\t{end}\t.\t+\t.\tID=gene:ENSG_A;biotype=protein_coding\n\
             1\ttest\tmRNA\t{CDS}\t{end}\t.\t+\t.\tID=transcript:ENST_A;Parent=gene:ENSG_A;biotype=protein_coding;tag=Ensembl_canonical\n\
             1\ttest\texon\t{CDS}\t{end}\t.\t+\t.\tID=exon:ENSE_A;Parent=transcript:ENST_A;rank=1\n\
             1\ttest\tCDS\t{CDS}\t{end}\t.\t+\t0\tID=CDS:ENSP_A;Parent=transcript:ENST_A\n"
        ),
    )
    .unwrap();
    path
}

fn write_fasta(dir: &Path) -> PathBuf {
    let mut seq = vec![b'N'; 1_100];
    seq[(CDS - 1) as usize..(CDS - 1) as usize + 12].copy_from_slice(b"ATGTGGCGGTAA");
    let path = dir.join("ref.fa");
    let mut f = File::create(&path).unwrap();
    writeln!(f, ">1").unwrap();
    for chunk in seq.chunks(60) {
        f.write_all(chunk).unwrap();
        f.write_all(b"\n").unwrap();
    }
    path
}

fn write_vcf(dir: &Path) -> PathBuf {
    let path = dir.join("in.vcf");
    let mut f = File::create(&path).unwrap();
    writeln!(f, "##fileformat=VCFv4.2").unwrap();
    writeln!(f, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO").unwrap();
    writeln!(f, "1\t{}\t.\tG\tA\t.\t.\t.", CDS + 7).unwrap();
    path
}

fn config(input: &Path, out: &Path) -> AnnotateConfig {
    AnnotateConfig {
        input: input.to_string_lossy().into(),
        output: out.to_string_lossy().into(),
        gff3: vec![],
        fasta: None,
        output_format: "vcf".into(),
        pick: PickFlags::default(),
        hgvs: true,
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
        columns: Default::default(),
    }
}

#[test]
fn a_truncated_explicit_transcript_cache_stops_the_run() {
    let dir = TempDir::new().unwrap();
    let gff3 = write_gff3(dir.path());
    let fasta = write_fasta(dir.path());
    let vcf = write_vcf(dir.path());

    // Build a real cache the ordinary way, and check the ordinary way works.
    let out_ok = dir.path().join("ok.vcf");
    run_annotate(AnnotateConfig {
        gff3: vec![gff3.to_string_lossy().into()],
        fasta: Some(fasta.to_string_lossy().into()),
        ..config(&vcf, &out_ok)
    })
    .expect("baseline annotation should succeed");
    let sidecar = dir.path().join("ann.gff3.fastvep.cache");
    assert!(sidecar.exists(), "expected a sidecar cache to be written");
    assert!(std::fs::read_to_string(&out_ok)
        .unwrap()
        .contains("missense_variant"));

    // Truncate it the way a concurrent writer or a full disk would.
    let full = std::fs::read(&sidecar).unwrap();
    let truncated = dir.path().join("truncated.cache");
    std::fs::write(&truncated, &full[..full.len() / 2]).unwrap();

    // Named explicitly and unusable: there is nothing to substitute, so the run
    // must fail rather than report the variant as intergenic.
    let out_bad = dir.path().join("bad.vcf");
    let err = run_annotate(AnnotateConfig {
        transcript_cache: Some(truncated.to_string_lossy().into()),
        fasta: Some(fasta.to_string_lossy().into()),
        ..config(&vcf, &out_bad)
    })
    .expect_err("a truncated --transcript-cache must not annotate anyway");
    let msg = err.to_string();
    assert!(
        msg.contains("cannot be used") && msg.contains("truncated or corrupt"),
        "error should name the cache as the problem, got: {msg}"
    );
    assert!(
        !out_bad.exists()
            || !std::fs::read_to_string(&out_bad)
                .unwrap()
                .contains("intergenic_variant"),
        "a failed cache load must not leave an all-intergenic output behind"
    );
}

#[test]
fn a_run_with_no_usable_transcript_source_stops_instead_of_calling_everything_intergenic() {
    let dir = TempDir::new().unwrap();
    let fasta = write_fasta(dir.path());
    let vcf = write_vcf(dir.path());
    let out = dir.path().join("out.vcf");

    let err = run_annotate(AnnotateConfig {
        fasta: Some(fasta.to_string_lossy().into()),
        ..config(&vcf, &out)
    })
    .expect_err("no --gff3 and no cache must not silently produce intergenic calls");
    let msg = err.to_string();
    assert!(
        msg.contains("No transcript source is usable"),
        "unexpected error: {msg}"
    );
}

#[test]
fn a_sidecar_cache_that_will_not_load_is_rebuilt_from_the_gff3() {
    // The counterpart: when the GFF3 is right there, a corrupt sidecar is a
    // recoverable condition and the run should repair itself rather than fail.
    let dir = TempDir::new().unwrap();
    let gff3 = write_gff3(dir.path());
    let fasta = write_fasta(dir.path());
    let vcf = write_vcf(dir.path());

    let sidecar = dir.path().join("ann.gff3.fastvep.cache");
    std::fs::write(&sidecar, b"FSTVEP02 and then garbage that will not inflate").unwrap();
    // Make it look fresh against its source, so the run really does try to load
    // it rather than skipping it as stale.
    let future = std::time::SystemTime::now() + std::time::Duration::from_secs(3600);
    File::options()
        .write(true)
        .open(&sidecar)
        .unwrap()
        .set_times(std::fs::FileTimes::new().set_modified(future))
        .unwrap();

    let out = dir.path().join("out.vcf");
    run_annotate(AnnotateConfig {
        gff3: vec![gff3.to_string_lossy().into()],
        fasta: Some(fasta.to_string_lossy().into()),
        ..config(&vcf, &out)
    })
    .expect("a corrupt sidecar should be rebuilt from the GFF3, not fatal");
    assert!(std::fs::read_to_string(&out)
        .unwrap()
        .contains("missense_variant"));
}

/// Write a well-formed cache in the pre-#90 FSTVEP02 format, holding `n`
/// transcripts. Deliberately valid: the whole difficulty is that a poisoned
/// region-restricted cache and a whole-file one are the same bytes with a
/// different transcript count, and neither the count nor the mtime tells them
/// apart.
fn write_v2_cache(path: &Path, transcripts: &[fastvep_genome::Transcript]) {
    use std::io::{BufWriter, Write};
    let file = File::create(path).unwrap();
    let mut zst = zstd::Encoder::new(BufWriter::new(file), 1).unwrap();
    zst.write_all(b"FSTVEP02").unwrap();
    bincode::serialize_into(&mut zst, transcripts).unwrap();
    zst.finish().unwrap().flush().unwrap();
}

fn make_fresh(path: &Path) {
    // An hour into the future, so `cache_is_fresh` says yes against its source
    // and the run really attempts the load rather than skipping it as stale.
    let future = std::time::SystemTime::now() + std::time::Duration::from_secs(3600);
    File::options()
        .write(true)
        .open(path)
        .unwrap()
        .set_times(std::fs::FileTimes::new().set_modified(future))
        .unwrap();
}

#[test]
fn a_pre_90_sidecar_cache_is_rebuilt_rather_than_trusted() {
    // #90 stopped persisting region-restricted transcript sets, but every cache
    // already on disk from a build between the tabix read path (2026-06-10) and
    // that fix may be one. This is the regression that gate leaves open: a fresh
    // FSTVEP02 sidecar, holding a transcript set that is *valid* and *wrong for
    // this input*, next to the GFF3 it was keyed on.
    let dir = TempDir::new().unwrap();
    let gff3 = write_gff3(dir.path());
    let fasta = write_fasta(dir.path());
    let vcf = write_vcf(dir.path());

    // An empty transcript set stands in for "someone else's neighbourhood":
    // deserializes cleanly, and every variant comes back intergenic.
    let sidecar = dir.path().join("ann.gff3.fastvep.cache");
    write_v2_cache(&sidecar, &[]);
    make_fresh(&sidecar);

    let out = dir.path().join("out.vcf");
    run_annotate(AnnotateConfig {
        gff3: vec![gff3.to_string_lossy().into()],
        fasta: Some(fasta.to_string_lossy().into()),
        ..config(&vcf, &out)
    })
    .expect("a pre-#90 sidecar should be rebuilt from the GFF3, not fatal");

    let annotated = std::fs::read_to_string(&out).unwrap();
    assert!(
        annotated.contains("missense_variant"),
        "the run trusted a pre-#90 cache and lost the annotation; got:\n{annotated}"
    );
    assert!(
        !annotated.contains("intergenic_variant"),
        "expected the rebuilt transcript set, not an all-intergenic answer"
    );
    // And the rebuild republished in the current format, so the next run is fast.
    let republished = std::fs::read(&sidecar).unwrap();
    let mut magic = [0u8; 8];
    {
        use std::io::Read;
        zstd::Decoder::new(&republished[..])
            .unwrap()
            .read_exact(&mut magic)
            .unwrap();
    }
    assert_eq!(
        &magic, b"FSTVEP08",
        "rebuild should publish the current format"
    );
}

/// Write a well-formed cache in the pre-#98 FSTVEP03 format. Whole-file and
/// perfectly readable - the only thing wrong with it is that the parser which
/// filled it did not recognise `ncRNA_gene`, so its non-coding genes have no
/// symbol.
fn write_v3_cache(path: &Path, transcripts: &[fastvep_genome::Transcript]) {
    use std::io::{BufWriter, Write};
    let file = File::create(path).unwrap();
    let mut zst = zstd::Encoder::new(BufWriter::new(file), 1).unwrap();
    zst.write_all(b"FSTVEP03").unwrap();
    bincode::serialize_into(&mut zst, transcripts).unwrap();
    zst.finish().unwrap().flush().unwrap();
}

#[test]
fn a_pre_98_sidecar_cache_is_rebuilt_rather_than_trusted() {
    // The gap #98 would otherwise leave open. Its reporter runs one command; the
    // sidecar next to their GFF3 is newer than the GFF3, so upgrading to the
    // fixed parser changes nothing unless the stale cache is refused - the
    // parser simply never runs again. Same shape as the pre-#90 case above, one
    // format later.
    let dir = TempDir::new().unwrap();
    let gff3 = write_gff3(dir.path());
    let fasta = write_fasta(dir.path());
    let vcf = write_vcf(dir.path());

    let sidecar = dir.path().join("ann.gff3.fastvep.cache");
    write_v3_cache(&sidecar, &[]);
    make_fresh(&sidecar);

    let out = dir.path().join("out.vcf");
    run_annotate(AnnotateConfig {
        gff3: vec![gff3.to_string_lossy().into()],
        fasta: Some(fasta.to_string_lossy().into()),
        ..config(&vcf, &out)
    })
    .expect("a pre-#98 sidecar should be rebuilt from the GFF3, not fatal");

    let annotated = std::fs::read_to_string(&out).unwrap();
    assert!(
        annotated.contains("missense_variant"),
        "the run trusted a pre-#98 cache and lost the annotation; got:\n{annotated}"
    );
}

/// Write a well-formed cache in the pre-APPRIS FSTVEP04 format. Whole-file,
/// every gene named - the only thing wrong with it is that the parser which
/// filled it did not read GENCODE's `appris_*` tags, so no transcript in it
/// carries an APPRIS call and `--pick`'s APPRIS tier has nothing to compare.
fn write_v4_cache(path: &Path, transcripts: &[fastvep_genome::Transcript]) {
    use std::io::{BufWriter, Write};
    let file = File::create(path).unwrap();
    let mut zst = zstd::Encoder::new(BufWriter::new(file), 1).unwrap();
    zst.write_all(b"FSTVEP04").unwrap();
    bincode::serialize_into(&mut zst, transcripts).unwrap();
    zst.finish().unwrap().flush().unwrap();
}

#[test]
fn a_pre_appris_sidecar_cache_is_rebuilt_rather_than_trusted() {
    // Third instance of the same gap: the sidecar is newer than its GFF3, so a
    // user who upgrades and re-runs the identical command never re-parses the
    // file, and the APPRIS column stays empty and the pick stays wrong. The
    // cache loads cleanly; nothing in it distinguishes "GENCODE model whose
    // APPRIS tags were dropped" from "Ensembl model that never had any".
    let dir = TempDir::new().unwrap();
    let gff3 = write_gff3(dir.path());
    let fasta = write_fasta(dir.path());
    let vcf = write_vcf(dir.path());

    let sidecar = dir.path().join("ann.gff3.fastvep.cache");
    write_v4_cache(&sidecar, &[]);
    make_fresh(&sidecar);

    let out = dir.path().join("out.vcf");
    run_annotate(AnnotateConfig {
        gff3: vec![gff3.to_string_lossy().into()],
        fasta: Some(fasta.to_string_lossy().into()),
        ..config(&vcf, &out)
    })
    .expect("a pre-APPRIS sidecar should be rebuilt from the GFF3, not fatal");

    let annotated = std::fs::read_to_string(&out).unwrap();
    assert!(
        annotated.contains("missense_variant"),
        "the run trusted a pre-APPRIS cache and lost the annotation; got:\n{annotated}"
    );
}

/// Write a well-formed cache in the pre-protein-version FSTVEP05 format. Whole-file, every
/// gene named, APPRIS read - the only thing wrong with it is that the parser which filled it
/// gave every protein version 1, whatever the GFF3's CDS records said.
fn write_v5_cache(path: &Path, transcripts: &[fastvep_genome::Transcript]) {
    use std::io::{BufWriter, Write};
    let file = File::create(path).unwrap();
    let mut zst = zstd::Encoder::new(BufWriter::new(file), 1).unwrap();
    zst.write_all(b"FSTVEP05").unwrap();
    bincode::serialize_into(&mut zst, transcripts).unwrap();
    zst.finish().unwrap().flush().unwrap();
}

#[test]
fn a_pre_protein_version_sidecar_cache_is_rebuilt_rather_than_trusted() {
    // Fourth instance of the same gap, and the one that shows it on the answer: the GFF3's
    // CDS records say the protein is version 3, the cache the old parser wrote says 1, and
    // the sidecar is newer than its GFF3, so a user who upgrades and re-runs the identical
    // command would keep `ENSP_A.1:p.` for ever if the stale cache were read.
    let dir = TempDir::new().unwrap();
    let gff3 = write_gff3(dir.path());
    // `ID=CDS:ENSP_A;Parent=...` is the last line: give the CDS a version
    let text = std::fs::read_to_string(&gff3).unwrap();
    std::fs::write(&gff3, format!("{};version=3\n", text.trim_end())).unwrap();
    let fasta = write_fasta(dir.path());
    let vcf = write_vcf(dir.path());
    let run = |out: &Path| {
        run_annotate(AnnotateConfig {
            gff3: vec![gff3.to_string_lossy().into()],
            fasta: Some(fasta.to_string_lossy().into()),
            ..config(&vcf, out)
        })
        .expect("annotation should succeed");
        std::fs::read_to_string(out).unwrap()
    };

    // the ordinary way: the version is the GFF3's
    let first = run(&dir.path().join("first.vcf"));
    assert!(
        first.contains("ENSP_A.3:p."),
        "the protein version is the CDS record's; got:\n{first}"
    );

    // what the old parser would have left beside the GFF3: the same transcripts, version 1
    let sidecar = dir.path().join("ann.gff3.fastvep.cache");
    let mut transcripts = fastvep_cache::transcript_cache::load_cache(&sidecar).unwrap();
    assert!(transcripts.iter().all(|t| t.protein_version == Some(3)));
    for t in &mut transcripts {
        t.protein_version = Some(1);
    }
    write_v5_cache(&sidecar, &transcripts);
    make_fresh(&sidecar);

    let second = run(&dir.path().join("second.vcf"));
    assert!(
        second.contains("ENSP_A.3:p.") && !second.contains("ENSP_A.1:p."),
        "the run trusted a pre-protein-version cache; got:\n{second}"
    );
}

/// Write a well-formed cache in the pre-RefSeq-naming FSTVEP06 format. Whole-file, APPRIS
/// read, protein versions right - the only thing wrong with it is that the parser which
/// filled it named a RefSeq transcript by its ID, its gene by its own ID and its biotype by
/// its SO term.
fn write_v6_cache(path: &Path, transcripts: &[fastvep_genome::Transcript]) {
    use std::io::{BufWriter, Write};
    let file = File::create(path).unwrap();
    let mut zst = zstd::Encoder::new(BufWriter::new(file), 1).unwrap();
    zst.write_all(b"FSTVEP06").unwrap();
    bincode::serialize_into(&mut zst, transcripts).unwrap();
    zst.finish().unwrap().flush().unwrap();
}

#[test]
fn a_pre_refseq_names_sidecar_cache_is_rebuilt_rather_than_trusted() {
    // Fifth instance of the same gap, and again it shows on the answer: a RefSeq GFF3 names a
    // transcript `rna-NM_000001.3` by its ID and `NM_000001.3` by its `transcript_id`, its
    // gene `gene-TP` and `7157` by its Dbxref, and Ensembl VEP writes the second of each. The
    // sidecar is newer than its GFF3, so a user who upgrades and re-runs the identical
    // command would keep the old names for ever if the stale cache were read.
    let dir = TempDir::new().unwrap();
    let gff3 = dir.path().join("refseq.gff3");
    let end = CDS + 11;
    std::fs::write(
        &gff3,
        format!(
            "1\tBestRefSeq\tgene\t{CDS}\t{end}\t.\t+\t.\tID=gene-TP;Dbxref=GeneID:7157;Name=TP;gbkey=Gene;gene=TP;gene_biotype=protein_coding\n\
             1\tBestRefSeq\tmRNA\t{CDS}\t{end}\t.\t+\t.\tID=rna-NM_000001.3;Parent=gene-TP;Dbxref=GeneID:7157,GenBank:NM_000001.3;gbkey=mRNA;gene=TP;transcript_id=NM_000001.3\n\
             1\tBestRefSeq\texon\t{CDS}\t{end}\t.\t+\t.\tID=exon-NM_000001.3-1;Parent=rna-NM_000001.3;transcript_id=NM_000001.3\n\
             1\tBestRefSeq\tCDS\t{CDS}\t{end}\t.\t+\t0\tID=cds-NP_000001.2;Parent=rna-NM_000001.3;gbkey=CDS;gene=TP;protein_id=NP_000001.2\n"
        ),
    )
    .unwrap();
    let fasta = write_fasta(dir.path());
    let vcf = write_vcf(dir.path());
    let run = |out: &Path| {
        run_annotate(AnnotateConfig {
            gff3: vec![gff3.to_string_lossy().into()],
            fasta: Some(fasta.to_string_lossy().into()),
            ..config(&vcf, out)
        })
        .expect("annotation should succeed");
        std::fs::read_to_string(out).unwrap()
    };

    // the ordinary way: Ensembl VEP's names
    let first = run(&dir.path().join("first.vcf"));
    assert!(
        first.contains("|7157|Transcript|NM_000001.3|protein_coding|"),
        "a RefSeq transcript is named as VEP names it; got:\n{first}"
    );
    assert!(first.contains("NM_000001.3:c.") && first.contains("NP_000001.2:p."));

    // what the old parser would have left beside the GFF3
    let sidecar = dir.path().join("refseq.gff3.fastvep.cache");
    let mut transcripts = fastvep_cache::transcript_cache::load_cache(&sidecar).unwrap();
    for t in &mut transcripts {
        t.stable_id = format!("rna-{}", t.stable_id).into();
        t.gene.stable_id = "gene-TP".into();
        t.biotype = "mRNA".into();
    }
    write_v6_cache(&sidecar, &transcripts);
    make_fresh(&sidecar);

    let second = run(&dir.path().join("second.vcf"));
    assert!(
        second.contains("|7157|Transcript|NM_000001.3|protein_coding|")
            && !second.contains("rna-NM_"),
        "the run trusted a pre-RefSeq-names cache; got:\n{second}"
    );
}

/// Write a well-formed cache in the pre-parity FSTVEP07 format. Whole-file, names right, protein
/// versions right - the only thing wrong with it is that the parser which filled it left a RefSeq
/// gene with no symbol source and no canonical transcript.
fn write_v7_cache(path: &Path, transcripts: &[fastvep_genome::Transcript]) {
    use std::io::{BufWriter, Write};
    let file = File::create(path).unwrap();
    let mut zst = zstd::Encoder::new(BufWriter::new(file), 1).unwrap();
    zst.write_all(b"FSTVEP07").unwrap();
    bincode::serialize_into(&mut zst, transcripts).unwrap();
    zst.finish().unwrap().flush().unwrap();
}

#[test]
fn a_pre_cache_parity_sidecar_cache_is_rebuilt_rather_than_trusted() {
    // Sixth instance of the same gap, and again it shows on the answer: Ensembl VEP's RefSeq
    // cache writes `EntrezGene` for the symbol source of a RefSeq gene and marks its MANE Select
    // transcript canonical, and the parser that wrote a V7 cache did neither. The sidecar is newer
    // than its GFF3, so a user who upgrades and re-runs the identical command would keep an empty
    // symbol source for ever if the stale cache were read.
    let dir = TempDir::new().unwrap();
    let gff3 = dir.path().join("refseq.gff3");
    let end = CDS + 11;
    std::fs::write(
        &gff3,
        format!(
            "1\tBestRefSeq\tgene\t{CDS}\t{end}\t.\t+\t.\tID=gene-TP;Dbxref=GeneID:7157,HGNC:HGNC:11998;Name=TP;gbkey=Gene;gene=TP;gene_biotype=protein_coding\n\
             1\tBestRefSeq\tmRNA\t{CDS}\t{end}\t.\t+\t.\tID=rna-NM_000001.3;Parent=gene-TP;Dbxref=Ensembl:ENST00000269305.9,GeneID:7157,GenBank:NM_000001.3;gbkey=mRNA;gene=TP;tag=MANE Select;transcript_id=NM_000001.3\n\
             1\tBestRefSeq\texon\t{CDS}\t{end}\t.\t+\t.\tID=exon-NM_000001.3-1;Parent=rna-NM_000001.3;transcript_id=NM_000001.3\n\
             1\tBestRefSeq\tCDS\t{CDS}\t{end}\t.\t+\t0\tID=cds-NP_000001.2;Parent=rna-NM_000001.3;gbkey=CDS;gene=TP;protein_id=NP_000001.2\n"
        ),
    )
    .unwrap();
    let fasta = write_fasta(dir.path());
    let vcf = write_vcf(dir.path());
    let run = |out: &Path| {
        run_annotate(AnnotateConfig {
            gff3: vec![gff3.to_string_lossy().into()],
            fasta: Some(fasta.to_string_lossy().into()),
            ..config(&vcf, out)
        })
        .expect("annotation should succeed");
        std::fs::read_to_string(out).unwrap()
    };
    // `SYMBOL_SOURCE`, `CANONICAL`, `MANE`, then `MANE_SELECT` follow `HGNC_ID` in the CSQ
    let wanted = "|YES|EntrezGene||MANE_Select|ENST00000269305.9|";

    let first = run(&dir.path().join("first.vcf"));
    assert!(
        first.contains(wanted),
        "a RefSeq gene is EntrezGene, its MANE Select transcript canonical; got:\n{first}"
    );

    // what the old parser would have left beside the GFF3
    let sidecar = dir.path().join("refseq.gff3.fastvep.cache");
    let mut transcripts = fastvep_cache::transcript_cache::load_cache(&sidecar).unwrap();
    for t in &mut transcripts {
        t.gene.symbol_source = None;
        t.canonical = false;
        t.mane_select = None;
    }
    write_v7_cache(&sidecar, &transcripts);
    make_fresh(&sidecar);

    let second = run(&dir.path().join("second.vcf"));
    assert!(
        second.contains(wanted),
        "the run trusted a pre-parity cache; got:\n{second}"
    );
}

#[test]
fn a_pre_90_explicit_transcript_cache_stops_the_run() {
    // Named explicitly there is nothing to fall back to, so the only safe answer
    // is to refuse - and to say why, because "rebuild it" is not obvious from
    // "invalid cache".
    let dir = TempDir::new().unwrap();
    let fasta = write_fasta(dir.path());
    let vcf = write_vcf(dir.path());

    let stale = dir.path().join("old.cache");
    write_v2_cache(&stale, &[]);

    let out = dir.path().join("out.vcf");
    let err = run_annotate(AnnotateConfig {
        transcript_cache: Some(stale.to_string_lossy().into()),
        fasta: Some(fasta.to_string_lossy().into()),
        ..config(&vcf, &out)
    })
    .expect_err("a pre-#90 --transcript-cache must not annotate anyway");

    let msg = format!("{err:#}");
    assert!(
        msg.contains("FSTVEP02") || msg.contains("pre-0.3.1"),
        "the error should name the format, got: {msg}"
    );
    assert!(
        msg.contains("fastvep cache"),
        "the error should say how to recover, got: {msg}"
    );
    // The recovery has to be one this path can actually perform. The cache was
    // named on the command line, so there is no sidecar to fall back to and
    // deleting the file conjures none - advice the message used to give because
    // a single fixed clause served both call sites.
    assert!(
        !msg.contains("let the sidecar rebuild itself"),
        "an explicitly named cache has no sidecar to fall back to, got: {msg}"
    );
    // One diagnosis, not two. The file is intact - sending the user to look for
    // a truncated write or a full disk would be a different problem than the
    // one they have, and the generic wording used to be appended to this one.
    assert!(
        !msg.contains("truncated or corrupt"),
        "a pre-#90 cache is intact; the error must not also blame corruption, got: {msg}"
    );
    assert!(
        !out.exists()
            || !std::fs::read_to_string(&out)
                .unwrap()
                .contains("intergenic_variant"),
        "a rejected cache must not leave an all-intergenic output behind"
    );
}
