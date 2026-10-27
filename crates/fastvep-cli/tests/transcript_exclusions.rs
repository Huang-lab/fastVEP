//! `--exclude-predicted`, `--exclude-readthrough` and `--exclude-artifact`: the transcripts Ensembl
//! VEP's caches do not hold, left out of a run on request.
//!
//! A GFF3 has every transcript of its annotation. VEP's Ensembl cache is built without the
//! readthrough ones and the ones of biotype `artifact` (`AnnotationSource::Database::Transcript`:
//! 109 of the 23,382 transcripts of chr21 and chr22 in Ensembl 116's GFF3, and the only ones of
//! them the cache lacks), and `--exclude_predicted` leaves RefSeq's `XM_` and `XR_` out of the
//! RefSeq one. None of it is on unless asked for.

use std::fs::File;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;

use fastvep_cli::pipeline::{run_annotate, AnnotateConfig, PickFlags, TranscriptExclusions};
use tempfile::TempDir;

/// One gene of each kind, each with one exon that a variant of the VCF falls in.
fn write_gff3(dir: &Path) -> PathBuf {
    let path = dir.join("mixed.gff3");
    let mut gff = String::new();
    // RefSeq: a reference transcript and a predicted one
    gff.push_str("1\tBestRefSeq\tgene\t1000\t1239\t.\t+\t.\tID=gene-A;Dbxref=GeneID:1;Name=A;gbkey=Gene;gene=A;gene_biotype=lncRNA\n");
    gff.push_str("1\tBestRefSeq\tlnc_RNA\t1000\t1239\t.\t+\t.\tID=rna-NR_000001.1;Parent=gene-A;Dbxref=GeneID:1;gbkey=ncRNA;gene=A;transcript_id=NR_000001.1\n");
    gff.push_str("1\tBestRefSeq\texon\t1000\t1239\t.\t+\t.\tID=exon-NR_000001.1-1;Parent=rna-NR_000001.1;transcript_id=NR_000001.1\n");
    gff.push_str("1\tGnomon\tgene\t2000\t2239\t.\t+\t.\tID=gene-B;Dbxref=GeneID:2;Name=B;gbkey=Gene;gene=B;gene_biotype=lncRNA\n");
    gff.push_str("1\tGnomon\tlnc_RNA\t2000\t2239\t.\t+\t.\tID=rna-XR_000002.1;Parent=gene-B;Dbxref=GeneID:2;gbkey=ncRNA;gene=B;transcript_id=XR_000002.1\n");
    gff.push_str("1\tGnomon\texon\t2000\t2239\t.\t+\t.\tID=exon-XR_000002.1-1;Parent=rna-XR_000002.1;transcript_id=XR_000002.1\n");
    // Ensembl: a readthrough transcript, an artifact, and an ordinary one
    for (n, start, extra) in [
        (
            "RT",
            3000,
            "biotype=lncRNA;tag=gencode_basic,readthrough_transcript",
        ),
        ("ART", 4000, "biotype=artifact"),
        ("OK", 5000, "biotype=lncRNA;tag=gencode_basic"),
    ] {
        let end = start + 239;
        gff.push_str(&format!(
            "1\tensembl\tncRNA_gene\t{start}\t{end}\t.\t+\t.\tID=gene:ENSG_{n};Name=G{n};biotype=lncRNA;gene_id=ENSG_{n}\n\
             1\tensembl\tlnc_RNA\t{start}\t{end}\t.\t+\t.\tID=transcript:ENST_{n};Parent=gene:ENSG_{n};{extra};transcript_id=ENST_{n}\n\
             1\tensembl\texon\t{start}\t{end}\t.\t+\t.\tParent=transcript:ENST_{n};rank=1\n"
        ));
    }
    // a coding transcript, so that a run has a sequence to build and writes the sidecar cache
    gff.push_str("1\tBestRefSeq\tgene\t6000\t6011\t.\t+\t.\tID=gene-C;Dbxref=GeneID:3;Name=C;gbkey=Gene;gene=C;gene_biotype=protein_coding\n");
    gff.push_str("1\tBestRefSeq\tmRNA\t6000\t6011\t.\t+\t.\tID=rna-NM_000003.1;Parent=gene-C;Dbxref=GeneID:3;gbkey=mRNA;gene=C;transcript_id=NM_000003.1\n");
    gff.push_str("1\tBestRefSeq\texon\t6000\t6011\t.\t+\t.\tID=exon-NM_000003.1-1;Parent=rna-NM_000003.1;transcript_id=NM_000003.1\n");
    gff.push_str("1\tBestRefSeq\tCDS\t6000\t6011\t.\t+\t0\tID=cds-NP_000003.1;Parent=rna-NM_000003.1;gbkey=CDS;gene=C;protein_id=NP_000003.1\n");
    std::fs::write(&path, gff).unwrap();
    path
}

fn write_fasta(dir: &Path) -> PathBuf {
    let mut seq = vec![b'N'; 6_100];
    seq[5_999..6_011].copy_from_slice(b"ATGTGGCGGTAA");
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
    for pos in [1100, 2100, 3100, 4100, 5100, 6005] {
        writeln!(f, "1\t{pos}\tv{pos}\tG\tA\t.\t.\t.").unwrap();
    }
    path
}

fn config(input: &Path, out: &Path, gff3: &Path, exclude: TranscriptExclusions) -> AnnotateConfig {
    let fasta = input.with_file_name("ref.fa");
    AnnotateConfig {
        input: input.to_string_lossy().into(),
        output: out.to_string_lossy().into(),
        gff3: vec![gff3.to_string_lossy().into()],
        fasta: fasta.exists().then(|| fasta.to_string_lossy().into()),
        output_format: "tab".into(),
        pick: PickFlags::default(),
        hgvs: false,
        distance: 0,
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
        exclude,
    }
}

/// The transcripts the run wrote rows for, `-` for the row of a variant that has none (the
/// intergenic one that a variant in a transcript left out becomes).
fn features(
    dir: &Path,
    name: &str,
    gff3: &Path,
    vcf: &Path,
    exclude: TranscriptExclusions,
) -> Vec<String> {
    let out = dir.join(name);
    run_annotate(config(vcf, &out, gff3, exclude)).expect("annotation should succeed");
    let tab = std::fs::read_to_string(&out).unwrap();
    let mut header = None;
    let mut features = Vec::new();
    for line in tab.lines().filter(|l| !l.starts_with("##")) {
        let cells: Vec<&str> = line.split('\t').collect();
        match header {
            None => {
                header = cells
                    .iter()
                    .position(|c| *c == "#Uploaded_variation")
                    .map(|_| cells.iter().position(|c| *c == "Feature").unwrap())
            }
            Some(i) => features.push(cells[i].to_string()),
        }
    }
    features.sort();
    features
}

#[test]
fn the_exclusion_flags_leave_out_what_vep_s_caches_do_not_hold() {
    let dir = TempDir::new().unwrap();
    let gff3 = write_gff3(dir.path());
    write_fasta(dir.path());
    let vcf = write_vcf(dir.path());
    let run = |name: &str, exclude| features(dir.path(), name, &gff3, &vcf, exclude);

    // The first run is the one that writes the sidecar cache, and it leaves everything out: what
    // it kept must not be what the next run reads.
    let everything = TranscriptExclusions {
        predicted: true,
        readthrough: true,
        artifact: true,
    };
    assert_eq!(
        run("all.tab", everything),
        ["-", "-", "-", "ENST_OK", "NM_000003.1", "NR_000001.1"]
    );
    let sidecar = dir.path().join("mixed.gff3.fastvep.cache");
    let cached = fastvep_cache::transcript_cache::load_cache(&sidecar).unwrap();
    assert_eq!(cached.len(), 6, "the cache is the whole file, flags or not");

    let all = [
        "ENST_ART",
        "ENST_OK",
        "ENST_RT",
        "NM_000003.1",
        "NR_000001.1",
        "XR_000002.1",
    ];
    assert_eq!(run("none.tab", TranscriptExclusions::default()), all);

    // `^X._`, which is VEP's `--exclude_predicted`
    let predicted = TranscriptExclusions {
        predicted: true,
        ..Default::default()
    };
    assert_eq!(
        run("predicted.tab", predicted),
        [
            "-",
            "ENST_ART",
            "ENST_OK",
            "ENST_RT",
            "NM_000003.1",
            "NR_000001.1"
        ]
    );
    let readthrough = TranscriptExclusions {
        readthrough: true,
        ..Default::default()
    };
    assert_eq!(
        run("readthrough.tab", readthrough),
        [
            "-",
            "ENST_ART",
            "ENST_OK",
            "NM_000003.1",
            "NR_000001.1",
            "XR_000002.1"
        ]
    );
    let artifact = TranscriptExclusions {
        artifact: true,
        ..Default::default()
    };
    assert_eq!(
        run("artifact.tab", artifact),
        [
            "-",
            "ENST_OK",
            "ENST_RT",
            "NM_000003.1",
            "NR_000001.1",
            "XR_000002.1"
        ]
    );
}

#[test]
fn leaving_out_every_transcript_stops_the_run() {
    // A run that has no transcript left would call every variant intergenic at exit code 0, the
    // failure an empty transcript source is refused for.
    let dir = TempDir::new().unwrap();
    let gff3 = dir.path().join("only_predicted.gff3");
    std::fs::write(
        &gff3,
        "1\tGnomon\tgene\t2000\t2239\t.\t+\t.\tID=gene-B;Dbxref=GeneID:2;Name=B;gbkey=Gene;gene=B;gene_biotype=lncRNA\n\
         1\tGnomon\tlnc_RNA\t2000\t2239\t.\t+\t.\tID=rna-XR_000002.1;Parent=gene-B;Dbxref=GeneID:2;gbkey=ncRNA;gene=B;transcript_id=XR_000002.1\n\
         1\tGnomon\texon\t2000\t2239\t.\t+\t.\tID=exon-XR_000002.1-1;Parent=rna-XR_000002.1;transcript_id=XR_000002.1\n",
    )
    .unwrap();
    let vcf = write_vcf(dir.path());
    let out = dir.path().join("out.tab");
    let err = run_annotate(config(
        &vcf,
        &out,
        &gff3,
        TranscriptExclusions {
            predicted: true,
            ..Default::default()
        },
    ))
    .expect_err("nothing left to annotate with");
    assert!(format!("{err:#}").contains("left out"), "{err:#}");
}

/// `--sa-only` loads no transcripts, so the exclusion flags have nothing to act on, and it says so
/// as it does for every other flag about transcripts. It used to count `0 of 0` and refuse the run
/// because "every transcript was left out", which is not what happened and sent the user looking
/// for a GFF3 they had not given.
#[test]
fn sa_only_ignores_the_exclusion_flags_and_says_so() {
    let dir = TempDir::new().unwrap();
    let vcf = write_vcf(dir.path());
    let sa = dir.path().join("sa");
    std::fs::create_dir(&sa).unwrap();
    for flag in [
        "--exclude-predicted",
        "--exclude-readthrough",
        "--exclude-artifact",
    ] {
        let out = dir.path().join("out.vcf");
        let output = Command::new(env!("CARGO_BIN_EXE_fastvep"))
            .arg("annotate")
            .arg("-i")
            .arg(&vcf)
            .arg("-o")
            .arg(&out)
            .arg("--sa-only")
            .arg("--sa-dir")
            .arg(&sa)
            .args(["--output-format", "vcf", "--no-progress", flag])
            .output()
            .expect("the fastvep binary runs");
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert!(output.status.success(), "{flag}: {stderr}");
        assert!(
            stderr.contains(&format!("--sa-only is set; ignoring {flag}")),
            "{flag}: {stderr}"
        );
        assert!(!stderr.contains("left out"), "{flag}: {stderr}");
    }
}
