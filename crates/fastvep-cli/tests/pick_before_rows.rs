//! A reducing `--pick*` chooses among a variant's transcripts before it builds their rows, and
//! builds only the rows it chooses. The flagging `--flag-pick*` builds every row and marks the
//! ones that `apply_pick` chooses. They are two paths to one answer, which this holds in worlds made
//! at random: for every scope, the rows a reducing run writes are the rows its flagging
//! counterpart marks, column for column and in order.
//!
//! The fixed fixtures of `flag_pick.rs` have two genes and one site, the smallest shape in which the
//! scopes differ. This has up to twenty genes of up to four transcripts, every attribute a pick
//! reads varied across them (canonical, MANE Select, MANE Plus Clinical, APPRIS, TSL, CCDS and
//! the biotype), and indels, multi-allelic records, an ALT that repeats another, a `*` allele and a
//! symbolic deletion, at distances from a few kilobases to enough that every transcript of the
//! world is beside the variant.
//!
//! The library's loop, the one the web server runs, builds every row and deletes the ones the pick
//! did not choose, so another test holds the CLI to the rows it keeps.
//!
//! What a run attaches to its rows from a supplementary annotation database (ClinVar here, built
//! from the variants of the world) is attached after the pick, to the rows it kept, so the same
//! comparisons are made with one.

use std::collections::BTreeMap;
use std::fmt::Write as _;
use std::path::{Path, PathBuf};

use fastvep_annotate::AnnotationContext;
use fastvep_cli::pipeline::{run_annotate, run_sa_build, AnnotateConfig, ColumnFlags, PickFlags};
use serde_json::Value;
use tempfile::TempDir;

/// A small deterministic generator, so that a failure is the same failure on every machine.
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }
    fn below(&mut self, n: u64) -> u64 {
        self.next() % n
    }
    fn chance(&mut self, one_in: u64) -> bool {
        self.below(one_in) == 0
    }
}

const CHROM_LEN: u64 = 220_000;

struct World {
    gff3: String,
    fasta: String,
    vcf: String,
}

fn world(seed: u64) -> World {
    let mut rng = Rng(seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1);
    for _ in 0..8 {
        rng.next();
    }

    // the genome
    let genome: Vec<u8> = (0..CHROM_LEN + 10_000)
        .map(|_| b"ACGT"[rng.below(4) as usize])
        .collect();
    let mut fasta = String::from(">1\n");
    for line in genome.chunks(60) {
        fasta.push_str(std::str::from_utf8(line).unwrap());
        fasta.push('\n');
    }

    // the genes: up to twenty, each of up to four transcripts that differ in what a pick reads
    let mut gff3 = String::new();
    let mut spans = Vec::new();
    for g in 0..5 + rng.below(16) {
        let start = 1_000 + rng.below(CHROM_LEN - 20_000);
        let end = start + 2_000 + rng.below(12_000);
        let strand = if rng.chance(2) { '+' } else { '-' };
        spans.push((start, end));
        writeln!(
            gff3,
            "1\ttest\tgene\t{start}\t{end}\t.\t{strand}\t.\tID=gene:ENSG{g:02};Name=GENE{g:02};biotype=protein_coding"
        )
        .unwrap();
        for t in 0..=rng.below(4) {
            let id = format!("ENST{g:02}_{t}");
            let biotype = match rng.below(10) {
                0..=6 => "protein_coding",
                7 | 8 => "lncRNA",
                _ => "processed_pseudogene",
            };
            let mut tags = Vec::new();
            if rng.chance(3) {
                tags.push("Ensembl_canonical");
            }
            if rng.chance(6) {
                tags.push("MANE_Select");
            }
            if rng.chance(5) {
                tags.push("MANE_Plus_Clinical");
            }
            if rng.chance(3) {
                tags.push(
                    [
                        "appris_principal_1",
                        "appris_principal_2",
                        "appris_alternative_1",
                    ][rng.below(3) as usize],
                );
            }
            let mut attributes =
                format!("ID=transcript:{id};Parent=gene:ENSG{g:02};biotype={biotype}");
            if !tags.is_empty() {
                write!(attributes, ";tag={}", tags.join(",")).unwrap();
            }
            if rng.chance(3) {
                write!(attributes, ";ccdsid=CCDS{g}{t}").unwrap();
            }
            if rng.chance(2) {
                write!(attributes, ";transcript_support_level={}", 1 + rng.below(5)).unwrap();
            }
            let tstart = start + rng.below(1_000);
            let mut tend = tstart;
            let mut exons = Vec::new();
            for _ in 0..=rng.below(3) {
                let exon_start = tend + 100 + rng.below(600);
                let exon_end = exon_start + 3 * (30 + rng.below(90)) - 1;
                exons.push((exon_start, exon_end));
                tend = exon_end;
            }
            let tstart = exons[0].0;
            writeln!(
                gff3,
                "1\ttest\tmRNA\t{tstart}\t{tend}\t.\t{strand}\t.\t{attributes}"
            )
            .unwrap();
            for (i, (exon_start, exon_end)) in exons.iter().enumerate() {
                writeln!(
                    gff3,
                    "1\ttest\texon\t{exon_start}\t{exon_end}\t.\t{strand}\t.\tID=exon:{id}_{i};Parent=transcript:{id};rank={}",
                    i + 1
                )
                .unwrap();
                if biotype == "protein_coding" {
                    writeln!(
                        gff3,
                        "1\ttest\tCDS\t{exon_start}\t{exon_end}\t.\t{strand}\t0\tID=CDS:{id}_{i};Parent=transcript:{id}"
                    )
                    .unwrap();
                }
            }
        }
    }

    // the variants: near the genes, mostly, so that they overlap transcripts and sit beside
    // others; indels and multi-allelic records among the substitutions
    let base = |pos: u64| genome[(pos - 1) as usize] as char;
    let mut records: Vec<(u64, String)> = Vec::new();
    for _ in 0..40 + rng.below(40) {
        let pos = if rng.chance(5) {
            2 + rng.below(CHROM_LEN)
        } else {
            let (s, e) = spans[rng.below(spans.len() as u64) as usize];
            s.saturating_sub(300) + rng.below(e - s + 600)
        }
        .max(2);
        let other = |rng: &mut Rng, avoid: char| loop {
            let b = "ACGT".chars().nth(rng.below(4) as usize).unwrap();
            if b != avoid {
                return b;
            }
        };
        let reference = base(pos);
        let (ref_allele, alts): (String, Vec<String>) = match rng.below(26) {
            0..=10 => (
                reference.to_string(),
                vec![other(&mut rng, reference).to_string()],
            ),
            11..=12 => {
                // a deletion of one to six bases, anchored
                let n = 1 + rng.below(6);
                let r: String = (0..=n).map(|i| base(pos + i)).collect();
                (r, vec![reference.to_string()])
            }
            13..=14 => {
                // an insertion of one to six bases
                let ins: String = (0..=rng.below(6)).map(|_| other(&mut rng, 'X')).collect();
                (reference.to_string(), vec![format!("{reference}{ins}")])
            }
            15..=17 => {
                let first = other(&mut rng, reference);
                let second = loop {
                    let b = other(&mut rng, reference);
                    if b != first {
                        break b;
                    }
                };
                let mut alts = vec![first.to_string(), second.to_string()];
                if rng.chance(3) {
                    alts.push("*".into());
                }
                (reference.to_string(), alts)
            }
            18 => (
                reference.to_string(),
                vec![
                    format!("{reference}{}", other(&mut rng, 'X')),
                    other(&mut rng, reference).to_string(),
                ],
            ),
            19 => (reference.to_string(), vec!["<DEL>".to_string()]),
            // the placeholders that assert no alternate sequence, alone and beside a real allele:
            // a record whose only ALT is one has no row for any transcript
            20 | 21 => (reference.to_string(), vec!["*".to_string()]),
            22 => (reference.to_string(), vec!["<NON_REF>".to_string()]),
            // an ALT that repeats another, which is one allele of the site that two columns of the
            // record name
            24 | 25 => {
                let alt = other(&mut rng, reference).to_string();
                (reference.to_string(), vec![alt.clone(), alt])
            }
            _ => (
                reference.to_string(),
                vec![
                    other(&mut rng, reference).to_string(),
                    "<NON_REF>".to_string(),
                ],
            ),
        };
        let info = if alts[0] == "<DEL>" {
            format!("SVTYPE=DEL;END={}", pos + 500 + rng.below(8_000))
        } else {
            ".".to_string()
        };
        records.push((
            pos,
            format!(
                "1\t{pos}\t.\t{ref_allele}\t{}\t.\t.\t{info}",
                alts.join(",")
            ),
        ));
    }
    records.sort_by_key(|(pos, _)| *pos);
    let mut vcf =
        String::from("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
    for (_, line) in records {
        vcf.push_str(&line);
        vcf.push('\n');
    }
    World { gff3, fasta, vcf }
}

fn config(input: &Path, out: &Path) -> AnnotateConfig {
    AnnotateConfig {
        input: input.to_string_lossy().into(),
        output: out.to_string_lossy().into(),
        gff3: vec![],
        fasta: None,
        output_format: "tab".into(),
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
        columns: ColumnFlags {
            symbol: true,
            canonical: true,
            biotype: true,
            numbers: true,
            hgvsg: true,
            ..Default::default()
        },
        exclude: Default::default(),
    }
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

struct Files {
    _dir: TempDir,
    gff3: PathBuf,
    fasta: PathBuf,
    vcf: PathBuf,
    out: PathBuf,
    /// A directory of supplementary annotation databases, for a run that reads them.
    sa: Option<PathBuf>,
}

fn files(world: &World) -> Files {
    let dir = TempDir::new().unwrap();
    let (gff3, fasta, vcf) = (
        dir.path().join("genes.gff3"),
        dir.path().join("genome.fa"),
        dir.path().join("in.vcf"),
    );
    std::fs::write(&gff3, &world.gff3).unwrap();
    std::fs::write(&fasta, &world.fasta).unwrap();
    std::fs::write(&vcf, &world.vcf).unwrap();
    let out = dir.path().join("out.tab");
    Files {
        _dir: dir,
        gff3,
        fasta,
        vcf,
        out,
        sa: None,
    }
}

/// A ClinVar database for the variants of a world: its plain substitutions, every other one, with
/// a significance and a review status that differ from record to record, so that what is attached
/// to a row is not one value for them all.
fn clinvar(dir: &Path, vcf: &str) -> PathBuf {
    let significance = [
        "Pathogenic",
        "Likely_pathogenic",
        "Uncertain_significance",
        "Likely_benign",
        "Benign",
    ];
    let review = [
        "criteria_provided,_multiple_submitters,_no_conflicts",
        "criteria_provided,_single_submitter",
        "no_assertion_criteria_provided",
    ];
    let mut text =
        String::from("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
    let mut kept = 0;
    for (i, line) in vcf.lines().filter(|l| !l.starts_with('#')).enumerate() {
        let f: Vec<&str> = line.split('\t').collect();
        let plain = f[3].len() == 1 && f[4].len() == 1 && "ACGT".contains(f[4]);
        if !plain || i % 2 == 1 {
            continue;
        }
        writeln!(
            text,
            "{}\t{}\t{i}\t{}\t{}\t.\t.\tCLNSIG={};CLNREVSTAT={};CLNDN=disease_{i};CLNVC=single_nucleotide_variant;CLNVCSO=SO:0001483",
            f[0],
            f[1],
            f[3],
            f[4],
            significance[i % significance.len()],
            review[i % review.len()]
        )
        .unwrap();
        kept += 1;
    }
    assert!(kept >= 3, "{kept} records of ClinVar");
    let source = dir.join("clinvar.vcf");
    std::fs::write(&source, text).unwrap();
    let sa = dir.join("sa");
    std::fs::create_dir_all(&sa).unwrap();
    run_sa_build(
        "clinvar",
        source.to_str().unwrap(),
        sa.join("clinvar").to_str().unwrap(),
        "GRCh38",
        None,
        &[],
        false,
    )
    .unwrap();
    sa
}

/// A world's files with the database of [`clinvar`] to read.
fn files_with_clinvar(world: &World) -> Files {
    let mut files = files(world);
    files.sa = Some(clinvar(files._dir.path(), &world.vcf));
    files
}

/// The header columns and the rows of a tab output.
fn read(files: &Files) -> (Vec<String>, Vec<Vec<String>>) {
    let text = std::fs::read_to_string(&files.out).unwrap();
    let mut header = Vec::new();
    let mut rows = Vec::new();
    for line in text.lines() {
        if line.starts_with("##") {
            continue;
        }
        let cells: Vec<String> = line.split('\t').map(str::to_string).collect();
        if line.starts_with('#') {
            header = cells;
        } else {
            rows.push(cells);
        }
    }
    (header, rows)
}

fn annotate_as(
    files: &Files,
    option: &str,
    distance: u64,
    order: Option<&str>,
    format: &str,
    acmg: bool,
) {
    run_annotate(AnnotateConfig {
        gff3: vec![files.gff3.to_string_lossy().into()],
        fasta: Some(files.fasta.to_string_lossy().into()),
        pick: flags_for(option),
        distance,
        pick_order: order.map(str::to_string),
        output_format: format.into(),
        acmg,
        sa_dir: files.sa.as_ref().map(|dir| dir.to_string_lossy().into()),
        ..config(&files.vcf, &files.out)
    })
    .unwrap_or_else(|e| {
        panic!("{option} at distance {distance}, order {order:?}, {format}, acmg {acmg}: {e:#}")
    });
}

const SCOPES: [(&str, &str); 3] = [
    ("--pick", "--flag-pick"),
    ("--pick-allele", "--flag-pick-allele"),
    ("--pick-allele-gene", "--flag-pick-allele-gene"),
];

/// Annotates with the flagging option and with its reducing counterpart, and holds the rows of the
/// second to the ones the first marked, column for column and in order. Returns the number of rows
/// the flagging run wrote and the number the reducing run did.
fn reducing_matches_flagging(
    files: &Files,
    (reduce, flag): (&str, &str),
    distance: u64,
    order: Option<&str>,
    acmg: bool,
    context: &str,
) -> (usize, usize) {
    annotate_as(files, flag, distance, order, "tab", acmg);
    let (flag_header, flag_rows) = read(files);
    let pick = flag_header
        .iter()
        .position(|c| c == "PICK")
        .expect("the flagging run declares a PICK column");
    let feature = flag_header
        .iter()
        .position(|c| c == "Feature")
        .expect("the output has a Feature column");
    let without_pick = |cells: &[String]| -> Vec<String> {
        cells
            .iter()
            .enumerate()
            .filter(|(i, _)| *i != pick)
            .map(|(_, c)| c.clone())
            .collect()
    };
    let expected_header = without_pick(&flag_header);
    // A variant that overlaps no transcript has a row with no feature, which no pick chooses or
    // deletes: the flagging run does not mark it and the reducing run keeps it.
    let expected_rows: Vec<Vec<String>> = flag_rows
        .iter()
        .filter(|row| row[pick] == "1" || row[feature] == "-")
        .map(|row| without_pick(row))
        .collect();

    annotate_as(files, reduce, distance, order, "tab", acmg);
    let (header, rows) = read(files);

    let context = format!("{context}, {reduce}, distance {distance}, order {order:?}, acmg {acmg}");
    assert_eq!(header, expected_header, "{context}: the header");

    assert_eq!(
        rows.len(),
        expected_rows.len(),
        "{context}: the number of rows"
    );
    for (i, (row, expected)) in rows.iter().zip(&expected_rows).enumerate() {
        assert_eq!(row, expected, "{context}: row {i}");
    }
    (flag_rows.len(), rows.len())
}

#[test]
fn what_a_reducing_pick_writes_is_what_its_flagging_counterpart_marks() {
    const ORDERS: [Option<&str>; 5] = [
        None,
        Some("canonical,biotype"),
        Some("rank,mane_select,canonical"),
        Some("tsl,appris,ccds,mane_select"),
        Some("mane_plus_clinical,canonical"),
    ];
    const DISTANCES: [u64; 4] = [0, 2_000, 30_000, 250_000];

    let (mut rows_before, mut rows_after) = (0usize, 0usize);
    for seed in 0..24 {
        let w = world(seed);
        let files = files(&w);
        let mut rng = Rng(seed + 1000);
        for round in 0..4 {
            let distance = DISTANCES[(seed as usize + round) % DISTANCES.len()];
            let order = ORDERS[rng.below(ORDERS.len() as u64) as usize];
            for scope in SCOPES {
                let (before, after) = reducing_matches_flagging(
                    &files,
                    scope,
                    distance,
                    order,
                    false,
                    &format!("seed {seed}"),
                );
                rows_before += before;
                rows_after += after;
            }
        }
    }
    // not vacuous: the picks left out most of what a plain run writes at these distances
    assert!(
        rows_after * 3 < rows_before,
        "{rows_after} rows kept of {rows_before}: the worlds were meant to give a pick much to choose from"
    );
}

/// The entries of a reducing run in JSON, ACMG classification included, are the entries that its
/// flagging counterpart marks, in the worlds of `seeds` with the files `files_for` makes of them.
/// Returns the number of entries compared, those that have a classification, and those with a
/// ClinVar record.
fn json_entries_match(files_for: fn(&World) -> Files, seeds: u64) -> (usize, usize, usize) {
    let entries = |files: &Files| -> Vec<Value> {
        let text = std::fs::read_to_string(&files.out).unwrap();
        serde_json::from_str(&text).unwrap()
    };
    let (mut compared, mut classified, mut with_clinvar) = (0usize, 0usize, 0usize);
    for seed in 0..seeds {
        let w = world(seed);
        let files = files_for(&w);
        for (reduce, flag) in SCOPES {
            for distance in [2_000, 250_000] {
                let context = format!("seed {seed}, {reduce}, distance {distance}");
                annotate_as(&files, flag, distance, None, "json", true);
                let flagged = entries(&files);
                annotate_as(&files, reduce, distance, None, "json", true);
                let reduced = entries(&files);

                assert_eq!(flagged.len(), reduced.len(), "{context}: variants");
                for (flagged, reduced) in flagged.iter().zip(&reduced) {
                    let kept: Vec<Value> = flagged["transcript_consequences"]
                        .as_array()
                        .into_iter()
                        .flatten()
                        // a variant that overlaps no transcript has a row with none, which no
                        // pick chooses or deletes
                        .filter(|tc| tc.get("pick").is_some() || tc["transcript_id"] == "-")
                        .map(|tc| {
                            let mut tc = tc.clone();
                            tc.as_object_mut().unwrap().remove("pick");
                            tc
                        })
                        .collect();
                    let got = reduced["transcript_consequences"]
                        .as_array()
                        .cloned()
                        .unwrap_or_default();
                    let at = format!(
                        "{context}: {}:{}",
                        flagged["seq_region_name"], flagged["start"]
                    );
                    assert_eq!(kept, got, "{at}: the transcript consequences");
                    assert_eq!(
                        flagged["intergenic_consequences"], reduced["intergenic_consequences"],
                        "{at}: the intergenic consequences"
                    );
                    compared += got.len();
                    classified += got.iter().filter(|tc| tc.get("acmg").is_some()).count();
                    with_clinvar += got.iter().filter(|tc| tc.get("clinvar").is_some()).count();
                }
            }
        }
    }
    (compared, classified, with_clinvar)
}

/// The ACMG classification reads what is left of a variant after the pick: its transcripts, its
/// genes, the transcript each row came from. A reducing run that picks before it builds the rows
/// has to leave all of it as a run that picks after it does, so each entry the reducing run
/// writes, ACMG classification included, is the entry its flagging counterpart marks.
///
/// In JSON, which is where a run writes the classification of each row.
#[test]
fn the_classification_of_the_rows_a_reducing_pick_keeps_is_the_one_its_flagging_counterpart_gives()
{
    let (compared, classified, _) = json_entries_match(files, 4);
    eprintln!("{classified} of {compared} entries have a classification");
    assert!(compared > 200, "{compared} entries compared");
    assert!(
        classified > 20,
        "{classified} of {compared} entries have a classification: the comparison is meant to include it"
    );
}

/// The same with ClinVar read: its record of a variant is attached to the rows the pick kept, and
/// the classification reads that record (PS1, PM5, PP5 and BP6 look at what ClinVar says), so it
/// is in what is compared.
#[test]
fn the_clinvar_records_and_classification_of_the_rows_a_reducing_pick_keeps_are_its_flagging_counterparts(
) {
    let (compared, classified, with_clinvar) = json_entries_match(files_with_clinvar, 4);
    eprintln!("{with_clinvar} of {compared} entries have ClinVar, {classified} a classification");
    assert!(compared > 200, "{compared} entries compared");
    assert!(
        with_clinvar > 20,
        "{with_clinvar} of {compared} entries have a ClinVar record: the comparison is meant to include them"
    );
    assert!(classified > 20, "{classified} classified");
}

/// What a tab output carries of ClinVar, the cells of its column of them that hold a value.
fn clinvar_cells(files: &Files) -> usize {
    let (header, rows) = read(files);
    let column = header
        .iter()
        .position(|c| c.eq_ignore_ascii_case("FV_CLINVAR"))
        .expect("a run that reads ClinVar has a column of it");
    rows.iter()
        .filter(|row| {
            row.get(column)
                .is_some_and(|cell| !cell.is_empty() && cell != "-")
        })
        .count()
}

/// In tab, the rows of a reducing pick with ClinVar read are the rows its flagging counterpart
/// marks, the column of ClinVar included, at every scope.
#[test]
fn the_clinvar_column_of_the_rows_a_reducing_pick_keeps_is_its_flagging_counterparts() {
    let (mut rows, mut cells) = (0usize, 0usize);
    for seed in 0..8 {
        let w = world(seed);
        let files = files_with_clinvar(&w);
        for scope in SCOPES {
            for distance in [2_000, 250_000] {
                let (_, after) = reducing_matches_flagging(
                    &files,
                    scope,
                    distance,
                    None,
                    false,
                    &format!("seed {seed}, ClinVar"),
                );
                rows += after;
                // the reducing run is the last that wrote the file
                cells += clinvar_cells(&files);
            }
        }
    }
    assert!(rows > 300, "{rows} rows compared");
    assert!(
        cells > 50,
        "{cells} of {rows} cells hold a ClinVar record: the comparison is meant to include them"
    );
}

/// The JSON and the VCF of a reducing pick have the rows of its tab output.
#[test]
fn a_reducing_pick_writes_the_same_rows_in_json_and_in_vcf_as_in_tab() {
    let mut compared = 0usize;
    for seed in 0..8 {
        let w = world(seed);
        let files = files(&w);
        for (reduce, _) in SCOPES {
            annotate_as(&files, reduce, 30_000, None, "tab", false);
            let tab = read(&files).1.len();

            annotate_as(&files, reduce, 30_000, None, "json", false);
            let text = std::fs::read_to_string(&files.out).unwrap();
            let json: usize = serde_json::from_str::<Vec<Value>>(&text)
                .unwrap()
                .iter()
                .map(|v| {
                    ["transcript_consequences", "intergenic_consequences"]
                        .iter()
                        .map(|key| v[*key].as_array().map_or(0, Vec::len))
                        .sum::<usize>()
                })
                .sum();

            annotate_as(&files, reduce, 30_000, None, "vcf", false);
            let text = std::fs::read_to_string(&files.out).unwrap();
            let vcf: usize = text
                .lines()
                .filter(|line| !line.starts_with('#'))
                .filter_map(|line| line.split('\t').nth(7))
                .filter_map(|info| info.split(';').find_map(|f| f.strip_prefix("CSQ=")))
                .map(|csq| csq.split(',').count())
                .sum();

            assert_eq!(
                json, tab,
                "seed {seed}, {reduce}: JSON entries and tab rows"
            );
            assert_eq!(vcf, tab, "seed {seed}, {reduce}: CSQ entries and tab rows");
            compared += tab;
        }
    }
    assert!(compared > 100, "{compared} rows compared");
}

/// What each loop writes for one (allele, transcript) of a variant that a pick can disturb.
const FIELDS: [&str; 7] = [
    "consequence_terms",
    "impact",
    "distance",
    "gene_id",
    "biotype",
    "hgvsc",
    "clinvar",
];

/// `(chromosome, start, allele_string, transcript, allele)`: the rows of a variant that repeats an
/// ALT have the one key twice.
type Key = (String, u64, String, String, String);

/// The transcript consequences of the JSON either loop writes, by row.
fn rows_of(variants: &[Value]) -> BTreeMap<Key, Vec<Vec<Value>>> {
    let mut rows: BTreeMap<Key, Vec<Vec<Value>>> = BTreeMap::new();
    for v in variants {
        for tc in v["transcript_consequences"]
            .as_array()
            .into_iter()
            .flatten()
        {
            let key = (
                v["seq_region_name"].as_str().unwrap().to_string(),
                v["start"].as_u64().unwrap(),
                v["allele_string"].as_str().unwrap().to_string(),
                tc["transcript_id"].as_str().unwrap_or("-").to_string(),
                tc["variant_allele"].as_str().unwrap_or("-").to_string(),
            );
            rows.entry(key)
                .or_default()
                .push(FIELDS.iter().map(|f| tc[*f].clone()).collect());
        }
    }
    rows
}

/// The library's loop, which the web server runs, builds every row of every transcript and deletes
/// the ones a reducing pick did not choose. The CLI chooses before it builds any, on what the
/// predictor said. They share the scoring and nothing else of the pick, which is the shape drift
/// takes, so the rows the CLI keeps are held to the rows the library keeps in the same worlds as
/// above, where the attributes a pick reads differ from transcript to transcript.
///
/// Without the symbolic deletions: the library's loop has no structural-variant predictor, so it
/// answers a `<DEL>` as it would a small variant, which is not what the CLI says of it and is
/// nothing to do with the pick.
#[test]
fn the_library_loop_keeps_the_rows_the_cli_keeps() {
    let (mut kept, mut with_clinvar) = (0usize, 0usize);
    for seed in 0..16 {
        let mut w = world(seed);
        w.vcf = w
            .vcf
            .lines()
            .filter(|line| !line.contains("<DEL>"))
            .map(|line| format!("{line}\n"))
            .collect();
        // with ClinVar, which both loops attach to the rows they kept
        let files = files_with_clinvar(&w);
        for option in ["--pick", "--pick-allele", "--pick-allele-gene"] {
            for distance in [2_000, 250_000] {
                annotate_as(&files, option, distance, None, "json", false);
                let text = std::fs::read_to_string(&files.out).unwrap();
                let cli = rows_of(&serde_json::from_str::<Vec<Value>>(&text).unwrap());

                let mut context = AnnotationContext::new(
                    Some(files.gff3.to_str().unwrap()),
                    Some(files.fasta.to_str().unwrap()),
                    files.sa.as_ref().map(|dir| dir.to_str().unwrap()),
                    distance,
                )
                .expect("the world should load");
                context.hgvs = true;
                let library = rows_of(
                    &context
                        .annotate_vcf_text(&w.vcf, flags_for(option).resolve())
                        .expect("annotation should succeed"),
                );

                let context = format!("seed {seed}, {option}, distance {distance}");
                let only = |a: &BTreeMap<Key, Vec<Vec<Value>>>,
                            b: &BTreeMap<Key, Vec<Vec<Value>>>| {
                    a.keys()
                        .filter(|k| !b.contains_key(*k))
                        .take(3)
                        .cloned()
                        .collect::<Vec<_>>()
                };
                assert_eq!(
                    only(&cli, &library),
                    Vec::<Key>::new(),
                    "{context}: rows only the CLI keeps"
                );
                assert_eq!(
                    only(&library, &cli),
                    Vec::<Key>::new(),
                    "{context}: rows only the library keeps"
                );
                for (key, rows) in &cli {
                    assert_eq!(rows, &library[key], "{context}: {key:?}");
                }
                kept += cli.len();
                with_clinvar += cli
                    .values()
                    .flatten()
                    .filter(|row| !row[FIELDS.len() - 1].is_null())
                    .count();
            }
        }
    }
    assert!(kept > 1_000, "{kept} rows compared");
    assert!(
        with_clinvar > 100,
        "{with_clinvar} of {kept} rows have a ClinVar record: the comparison is meant to include them"
    );
}
