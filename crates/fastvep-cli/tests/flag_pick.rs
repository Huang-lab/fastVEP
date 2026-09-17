//! `--flag-pick*` and `--pick-allele*`: VEP's pick family beyond `--pick`.
//!
//! Driven through `run_annotate` rather than `apply_pick`, for the reason
//! `appris_pick_tier.rs` gives at length: the unit tests build
//! `TranscriptVariation` by hand, so they cannot see a flag that never reaches
//! the writer, a `PICK` column missing from the header, or a field list that
//! disagrees with the values written under it. All three are how this feature
//! would arrive looking finished and be useless.
//!
//! The fixture is two overlapping protein-coding genes, each with its own
//! canonical transcript, and a two-alt site inside both - the smallest shape
//! where the three scopes give three different answers.

use std::fs::File;
use std::io::Write;
use std::path::{Path, PathBuf};

use fastvep_cli::pipeline::{run_annotate, AnnotateConfig, PickFlags};
use tempfile::TempDir;

/// Two overlapping genes, two transcripts each.
///
/// Within each gene one transcript is `Ensembl_canonical` and the other is
/// not, so the `canonical` tier decides the within-gene pick. Across genes
/// every tier ties, so the cross-gene pick falls to the transcript-ID
/// tie-break: `ENST_A1` < `ENST_B1`.
fn write_gff3(dir: &Path) -> PathBuf {
    let path = dir.join("ann.gff3");
    std::fs::write(
        &path,
        "1\ttest\tgene\t1001\t1200\t.\t+\t.\tID=gene:ENSG_A;Name=GENEA;biotype=protein_coding\n\
         1\ttest\tmRNA\t1001\t1200\t.\t+\t.\tID=transcript:ENST_A1;Parent=gene:ENSG_A;biotype=protein_coding;tag=Ensembl_canonical;transcript_support_level=1\n\
         1\ttest\texon\t1001\t1200\t.\t+\t.\tID=exon:E_A1;Parent=transcript:ENST_A1;rank=1\n\
         1\ttest\tCDS\t1001\t1200\t.\t+\t0\tID=CDS:P_A1;Parent=transcript:ENST_A1\n\
         1\ttest\tmRNA\t1001\t1200\t.\t+\t.\tID=transcript:ENST_A2;Parent=gene:ENSG_A;biotype=protein_coding;transcript_support_level=1\n\
         1\ttest\texon\t1001\t1200\t.\t+\t.\tID=exon:E_A2;Parent=transcript:ENST_A2;rank=1\n\
         1\ttest\tCDS\t1001\t1200\t.\t+\t0\tID=CDS:P_A2;Parent=transcript:ENST_A2\n\
         1\ttest\tgene\t1001\t1200\t.\t+\t.\tID=gene:ENSG_B;Name=GENEB;biotype=protein_coding\n\
         1\ttest\tmRNA\t1001\t1200\t.\t+\t.\tID=transcript:ENST_B1;Parent=gene:ENSG_B;biotype=protein_coding;tag=Ensembl_canonical;transcript_support_level=1\n\
         1\ttest\texon\t1001\t1200\t.\t+\t.\tID=exon:E_B1;Parent=transcript:ENST_B1;rank=1\n\
         1\ttest\tCDS\t1001\t1200\t.\t+\t0\tID=CDS:P_B1;Parent=transcript:ENST_B1\n\
         1\ttest\tmRNA\t1001\t1200\t.\t+\t.\tID=transcript:ENST_B2;Parent=gene:ENSG_B;biotype=protein_coding;transcript_support_level=1\n\
         1\ttest\texon\t1001\t1200\t.\t+\t.\tID=exon:E_B2;Parent=transcript:ENST_B2;rank=1\n\
         1\ttest\tCDS\t1001\t1200\t.\t+\t0\tID=CDS:P_B2;Parent=transcript:ENST_B2\n",
    )
    .unwrap();
    path
}

fn write_fasta(dir: &Path) -> PathBuf {
    let mut seq = vec![b'A'; 1_300];
    seq[1000..1003].copy_from_slice(b"ATG");
    let path = dir.join("ref.fa");
    let mut f = File::create(&path).unwrap();
    writeln!(f, ">1").unwrap();
    for chunk in seq.chunks(60) {
        f.write_all(chunk).unwrap();
        f.write_all(b"\n").unwrap();
    }
    path
}

/// One record, two alt alleles, inside both genes' single exon.
fn write_vcf(dir: &Path) -> PathBuf {
    let path = dir.join("in.vcf");
    let mut f = File::create(&path).unwrap();
    writeln!(f, "##fileformat=VCFv4.2").unwrap();
    writeln!(f, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO").unwrap();
    writeln!(f, "1\t1100\t.\tA\tC,G\t.\t.\t.").unwrap();
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
        hgvs: false,
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
    }
}

fn annotate(dir: &Path, flags: PickFlags, format: &str) -> String {
    let gff3 = write_gff3(dir);
    let fasta = write_fasta(dir);
    let vcf = write_vcf(dir);
    let out = dir.join(format!(
        "out-{}-{}",
        flags.requested_option().unwrap_or("none").trim_matches('-'),
        format
    ));
    run_annotate(AnnotateConfig {
        gff3: vec![gff3.to_string_lossy().into()],
        fasta: Some(fasta.to_string_lossy().into()),
        pick: flags,
        output_format: format.into(),
        ..config(&vcf, &out)
    })
    .expect("annotation should succeed");
    std::fs::read_to_string(&out).unwrap()
}

/// The CSQ field names this run declared, read out of the header it wrote
/// rather than assumed - the point of several of these assertions is that the
/// header and the values agree.
fn declared_fields(annotated: &str) -> Vec<String> {
    let line = annotated
        .lines()
        .find(|l| l.starts_with("##INFO=<ID=CSQ"))
        .unwrap_or_else(|| panic!("no CSQ header in:\n{annotated}"));
    let fmt = line
        .split("Format: ")
        .nth(1)
        .and_then(|rest| rest.split("\">").next())
        .expect("CSQ header declares a format");
    fmt.split('|').map(str::to_string).collect()
}

/// `(allele, symbol, transcript, pick)` per CSQ entry, with `pick` read by the
/// column's declared position.
fn entries(annotated: &str) -> Vec<(String, String, String, String)> {
    let fields = declared_fields(annotated);
    let idx = |name: &str| fields.iter().position(|f| f == name);
    let line = annotated
        .lines()
        .find(|l| !l.starts_with('#'))
        .unwrap_or_else(|| panic!("no annotated record in:\n{annotated}"));
    let csq = line
        .split('\t')
        .nth(7)
        .expect("INFO column")
        .split(';')
        .find_map(|kv| kv.strip_prefix("CSQ="))
        .unwrap_or_else(|| panic!("no CSQ on:\n{line}"));
    csq.split(',')
        .map(|entry| {
            let v: Vec<&str> = entry.split('|').collect();
            let get = |name: &str| {
                idx(name)
                    .and_then(|i| v.get(i))
                    .map(|s| s.to_string())
                    .unwrap_or_default()
            };
            (
                get("Allele"),
                get("SYMBOL"),
                get("Feature"),
                idx("PICK")
                    .and_then(|i| v.get(i))
                    .map(|s| s.to_string())
                    .unwrap_or_default(),
            )
        })
        .collect()
}

fn flags_for(option: &str) -> PickFlags {
    let mut f = PickFlags::default();
    match option {
        "--pick" => f.pick = true,
        "--pick-allele" => f.pick_allele = true,
        "--pick-allele-gene" => f.pick_allele_gene = true,
        "--flag-pick" => f.flag_pick = true,
        "--flag-pick-allele" => f.flag_pick_allele = true,
        "--flag-pick-allele-gene" => f.flag_pick_allele_gene = true,
        other => panic!("unknown option {other}"),
    }
    f
}

#[test]
fn a_plain_run_reports_every_transcript_and_declares_no_pick_column() {
    let dir = TempDir::new().unwrap();
    let out = annotate(dir.path(), PickFlags::default(), "vcf");
    // Four transcripts x two alleles.
    assert_eq!(entries(&out).len(), 8, "in:\n{out}");
    assert!(
        !declared_fields(&out).contains(&"PICK".to_string()),
        "a run that picks nothing should not declare PICK: {:?}",
        declared_fields(&out)
    );
}

#[test]
fn each_reducing_scope_leaves_the_entries_its_name_promises() {
    for (option, expected) in [
        // One transcript, both of its alleles - fastVEP's documented
        // divergence from VEP's `--pick`, which would leave one entry.
        (
            "--pick",
            vec![("C", "GENEA", "ENST_A1"), ("G", "GENEA", "ENST_A1")],
        ),
        // One per allele; every tier ties across the genes, so the
        // transcript-ID tie-break picks GENEA for both.
        (
            "--pick-allele",
            vec![("C", "GENEA", "ENST_A1"), ("G", "GENEA", "ENST_A1")],
        ),
        // One per allele per gene: the canonical transcript of each.
        (
            "--pick-allele-gene",
            vec![
                ("C", "GENEA", "ENST_A1"),
                ("G", "GENEA", "ENST_A1"),
                ("C", "GENEB", "ENST_B1"),
                ("G", "GENEB", "ENST_B1"),
            ],
        ),
    ] {
        let dir = TempDir::new().unwrap();
        let out = annotate(dir.path(), flags_for(option), "vcf");
        let got: Vec<(String, String, String)> = entries(&out)
            .into_iter()
            .map(|(a, s, t, _)| (a, s, t))
            .collect();
        let want: Vec<(String, String, String)> = expected
            .into_iter()
            .map(|(a, s, t)| (a.into(), s.into(), t.into()))
            .collect();
        // Order within a record is not part of the contract; membership is.
        let mut got_sorted = got.clone();
        got_sorted.sort();
        let mut want_sorted = want.clone();
        want_sorted.sort();
        assert_eq!(got_sorted, want_sorted, "{option} in:\n{out}");
        assert!(
            !declared_fields(&out).contains(&"PICK".to_string()),
            "{option} reduces, so it should not declare PICK"
        );
    }
}

#[test]
fn each_flagging_scope_keeps_everything_and_marks_what_it_would_have_kept() {
    for option in [
        "--flag-pick",
        "--flag-pick-allele",
        "--flag-pick-allele-gene",
    ] {
        let reducing = option
            .strip_prefix("--flag")
            .map(|r| format!("-{r}"))
            .unwrap();

        let dir = TempDir::new().unwrap();
        let flagged = annotate(dir.path(), flags_for(option), "vcf");
        let dir2 = TempDir::new().unwrap();
        let reduced = annotate(dir2.path(), flags_for(&reducing), "vcf");

        assert!(
            declared_fields(&flagged).contains(&"PICK".to_string()),
            "{option} must declare the PICK column it writes"
        );
        let all = entries(&flagged);
        assert_eq!(all.len(), 8, "{option} dropped entries in:\n{flagged}");

        let mut marked: Vec<(String, String, String)> = all
            .iter()
            .filter(|(_, _, _, pick)| pick == "1")
            .map(|(a, s, t, _)| (a.clone(), s.clone(), t.clone()))
            .collect();
        let mut kept: Vec<(String, String, String)> = entries(&reduced)
            .into_iter()
            .map(|(a, s, t, _)| (a, s, t))
            .collect();
        marked.sort();
        kept.sort();
        assert_eq!(
            marked, kept,
            "{option} flagged a different set than {reducing} keeps"
        );

        // Every other entry says so explicitly rather than being blank-by-accident.
        assert!(
            all.iter().any(|(_, _, _, pick)| pick.is_empty()),
            "{option} should leave the column empty on the entries it did not pick"
        );
    }
}

#[test]
fn the_tab_writer_carries_the_same_flag_under_a_declared_column() {
    let dir = TempDir::new().unwrap();
    let out = annotate(dir.path(), flags_for("--flag-pick-allele-gene"), "tab");

    let header = out
        .lines()
        .find(|l| l.starts_with("#Uploaded_variation"))
        .unwrap_or_else(|| panic!("no tab header in:\n{out}"));
    let cols: Vec<&str> = header.split('\t').collect();
    let pick_col = cols
        .iter()
        .position(|c| *c == "PICK")
        .unwrap_or_else(|| panic!("tab header should declare PICK: {header}"));

    let rows: Vec<Vec<&str>> = out
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| l.split('\t').collect())
        .collect();
    assert_eq!(rows.len(), 8, "every row is retained in:\n{out}");
    for row in &rows {
        assert_eq!(
            row.len(),
            cols.len(),
            "row has {} cells for {} declared columns: {row:?}",
            row.len(),
            cols.len()
        );
    }
    let flagged = rows.iter().filter(|r| r[pick_col] == "1").count();
    assert_eq!(flagged, 4, "one per allele per gene in:\n{out}");
    // Unflagged rows carry `-`, like every other empty cell in this format.
    assert!(
        rows.iter()
            .all(|r| r[pick_col] == "1" || r[pick_col] == "-"),
        "PICK cells should be 1 or -, in:\n{out}"
    );
}

#[test]
fn the_json_writer_marks_only_the_picked_consequences() {
    let dir = TempDir::new().unwrap();
    let out = annotate(dir.path(), flags_for("--flag-pick-allele-gene"), "json");
    let parsed: serde_json::Value = serde_json::from_str(&out).expect("valid JSON array");
    let consequences = parsed[0]["transcript_consequences"]
        .as_array()
        .expect("transcript_consequences array");
    assert_eq!(consequences.len(), 8, "every consequence is retained");
    let picked = consequences
        .iter()
        .filter(|tc| tc.get("pick").and_then(|v| v.as_u64()) == Some(1))
        .count();
    assert_eq!(picked, 4, "one per allele per gene in:\n{out}");

    // A plain run leaves the key off entirely rather than writing `0`.
    let dir2 = TempDir::new().unwrap();
    let plain = annotate(dir2.path(), PickFlags::default(), "json");
    let parsed: serde_json::Value = serde_json::from_str(&plain).unwrap();
    assert!(
        parsed[0]["transcript_consequences"]
            .as_array()
            .unwrap()
            .iter()
            .all(|tc| tc.get("pick").is_none()),
        "a run that picks nothing should not mention pick:\n{plain}"
    );
}
