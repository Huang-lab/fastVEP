//! The optional columns of VEP's tab output, asked for by VEP's own flag names
//! (`--hgvs`, `--hgvsg`, `--symbol`, `--biotype`, `--canonical`, `--numbers`,
//! `--variant-class`, `--show-ref-allele`), end to end (issue #139).
//!
//! VEP adds a tab column only when its flag is given, and in the order of the
//! `FLAG_FIELDS` table of `Constants.pm` (ensembl-vep release/116). fastVEP's
//! tab output had 17 fixed columns and accepted `--hgvs`, `--symbol` and
//! `--canonical` without writing anything for them. The reporter's pipeline
//! reads `BIOTYPE`, `EXON`, `INTRON`, `HGVSc`, `HGVSp`, `HGVS_OFFSET`, `HGVSg`,
//! `VARIANT_CLASS` and `REF_ALLELE` from the tab file.
//!
//! The tab and the VCF writers share the value of each of these columns, so a
//! test below holds one against the other: that is how they would come to
//! disagree.

use std::collections::HashMap;
use std::fs::File;
use std::io::Write;
use std::path::{Path, PathBuf};

use fastvep_cli::pipeline::{run_annotate, AnnotateConfig, ColumnFlags, PickFlags};
use tempfile::TempDir;

/// One canonical lncRNA, `chr10:11730-14000`, exons 11730-12000 and 13000-14000.
fn write_gff3(dir: &Path) -> PathBuf {
    let path = dir.join("genes.gff3");
    std::fs::write(
        &path,
        "chr10\ttest\tgene\t11730\t14000\t.\t+\t.\tID=gene:ENSG_X;Name=GENEX;biotype=lncRNA\n\
         chr10\ttest\tlnc_RNA\t11730\t14000\t.\t+\t.\tID=transcript:ENST_X;Parent=gene:ENSG_X;biotype=lncRNA;tag=Ensembl_canonical\n\
         chr10\ttest\texon\t11730\t12000\t.\t+\t.\tID=exon:E1;Parent=transcript:ENST_X;rank=1\n\
         chr10\ttest\texon\t13000\t14000\t.\t+\t.\tID=exon:E2;Parent=transcript:ENST_X;rank=2\n",
    )
    .unwrap();
    path
}

fn write_vcf(dir: &Path) -> PathBuf {
    let path = dir.join("in.vcf");
    let mut f = File::create(&path).unwrap();
    writeln!(f, "##fileformat=VCFv4.2").unwrap();
    writeln!(f, "##contig=<ID=chr10,length=133797422>").unwrap();
    writeln!(f, "##contig=<ID=chrM,length=16569>").unwrap();
    writeln!(f, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO").unwrap();
    // an SNV in the first exon
    writeln!(f, "chr10\t11800\t.\tG\tT\t.\t.\t.").unwrap();
    // a deletion in the first exon (GA>G is A/- at 11901)
    writeln!(f, "chr10\t11900\t.\tGA\tG\t.\t.\t.").unwrap();
    // an insertion in the first exon
    writeln!(f, "chr10\t11950\t.\tC\tCTT\t.\t.\t.").unwrap();
    // an SNV in the intron between the two exons
    writeln!(f, "chr10\t12500\t.\tA\tC\t.\t.\t.").unwrap();
    // no transcript near
    writeln!(f, "chrM\t100\t.\tA\tG\t.\t.\t.").unwrap();
    path
}

fn run(dir: &Path, format: &str, columns: ColumnFlags, hgvs: bool) -> String {
    run_vcf(dir, &write_vcf(dir), format, columns, hgvs)
}

fn run_vcf(dir: &Path, vcf: &Path, format: &str, columns: ColumnFlags, hgvs: bool) -> String {
    let out = dir.join(format!("out-{}.{format}", dir.read_dir().unwrap().count()));
    run_annotate(AnnotateConfig {
        input: vcf.to_string_lossy().into(),
        output: out.to_string_lossy().into(),
        gff3: vec![write_gff3(dir).to_string_lossy().into()],
        fasta: None,
        output_format: format.into(),
        pick: PickFlags::default(),
        hgvs,
        columns,
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
    })
    .expect("annotation should succeed");
    std::fs::read_to_string(&out).unwrap()
}

fn header(tab: &str) -> Vec<String> {
    tab.lines()
        .find(|l| l.starts_with("#Uploaded_variation"))
        .unwrap_or_else(|| panic!("no column header in:\n{tab}"))
        .trim_start_matches('#')
        .split('\t')
        .map(String::from)
        .collect()
}

fn rows(tab: &str) -> Vec<HashMap<String, String>> {
    let names = header(tab);
    tab.lines()
        .filter(|l| !l.starts_with('#'))
        .map(|row| {
            let cells: Vec<&str> = row.split('\t').collect();
            assert_eq!(cells.len(), names.len(), "ragged row: {row:?}");
            names
                .iter()
                .cloned()
                .zip(cells.into_iter().map(String::from))
                .collect()
        })
        .collect()
}

fn find<'a>(rows: &'a [HashMap<String, String>], uploaded: &str) -> &'a HashMap<String, String> {
    rows.iter()
        .find(|r| r["Uploaded_variation"] == uploaded)
        .unwrap_or_else(|| panic!("no row {uploaded} in {rows:#?}"))
}

const BASE_17: [&str; 17] = [
    "Uploaded_variation",
    "Location",
    "Allele",
    "Gene",
    "Feature",
    "Feature_type",
    "Consequence",
    "cDNA_position",
    "CDS_position",
    "Protein_position",
    "Amino_acids",
    "Codons",
    "Existing_variation",
    "IMPACT",
    "DISTANCE",
    "STRAND",
    "FLAGS",
];

/// With no flag, the tab file is the 17 columns it always was.
#[test]
fn without_flags_the_tab_file_keeps_its_17_columns() {
    let dir = TempDir::new().unwrap();
    let tab = run(dir.path(), "tab", ColumnFlags::default(), false);
    assert_eq!(header(&tab), BASE_17);
}

/// Each flag adds the columns VEP's `FLAG_FIELDS` gives it, in that table's order:
/// REF_ALLELE before IMPACT, then (after FLAGS) VARIANT_CLASS, the SYMBOL group,
/// BIOTYPE, CANONICAL, EXON and INTRON, the HGVS group, HGVSg. That is the header
/// of a VEP 116.3 run with all of them.
#[test]
fn the_flags_add_vep_s_columns_in_vep_s_order() {
    let dir = TempDir::new().unwrap();
    let all = ColumnFlags {
        show_ref_allele: true,
        variant_class: true,
        symbol: true,
        biotype: true,
        canonical: true,
        numbers: true,
        hgvsg: true,
    };
    let tab = run(dir.path(), "tab", all, true);
    let mut want: Vec<&str> = BASE_17[..13].to_vec();
    want.push("REF_ALLELE");
    want.extend(&BASE_17[13..]);
    want.extend([
        "VARIANT_CLASS",
        "SYMBOL",
        "SYMBOL_SOURCE",
        "HGNC_ID",
        "BIOTYPE",
        "CANONICAL",
        "EXON",
        "INTRON",
        "HGVSc",
        "HGVSp",
        "HGVS_OFFSET",
        "HGVSg",
    ]);
    assert_eq!(header(&tab), want);

    // and one flag adds only its own
    let dir = TempDir::new().unwrap();
    let only_biotype = ColumnFlags {
        biotype: true,
        ..ColumnFlags::default()
    };
    let tab = run(dir.path(), "tab", only_biotype, false);
    let mut want: Vec<&str> = BASE_17.to_vec();
    want.push("BIOTYPE");
    assert_eq!(header(&tab), want);
}

#[test]
fn hgvs_adds_the_hgvs_group_and_not_hgvsg() {
    let dir = TempDir::new().unwrap();
    let tab = run(dir.path(), "tab", ColumnFlags::default(), true);
    let mut want: Vec<&str> = BASE_17.to_vec();
    want.extend(["HGVSc", "HGVSp", "HGVS_OFFSET"]);
    assert_eq!(header(&tab), want);

    let rows = rows(&tab);
    let snv = find(&rows, "chr10_11800_G/T");
    // base 71 of a non-coding transcript
    assert_eq!(snv["HGVSc"], "ENST_X:n.71G>T");
    assert_eq!(snv["HGVSp"], "-");
    assert_eq!(snv["HGVS_OFFSET"], "-");
}

#[test]
fn hgvsg_is_written_on_every_row_including_the_intergenic_one() {
    let dir = TempDir::new().unwrap();
    let flags = ColumnFlags {
        hgvsg: true,
        ..ColumnFlags::default()
    };
    // `--hgvsg` alone: no HGVSc, and still HGVSg
    let tab = run(dir.path(), "tab", flags, false);
    assert!(header(&tab).contains(&"HGVSg".to_string()));
    assert!(!header(&tab).contains(&"HGVSc".to_string()));
    let rows = rows(&tab);
    assert_eq!(find(&rows, "chr10_11800_G/T")["HGVSg"], "chr10:g.11800G>T");
    assert_eq!(find(&rows, "chrM_100_A/G")["HGVSg"], "chrM:g.100A>G");
}

/// `--hgvs` has never put an HGVSg on the intergenic row, and the web
/// annotator (`fastvep-annotate`) still does not, so `--hgvs` alone must not
/// start here: that is one loop changing without the other. `--hgvsg` is what
/// asks for it.
#[test]
fn only_hgvsg_puts_a_genomic_hgvs_on_the_intergenic_row() {
    let intergenic_hgvsg = |columns: ColumnFlags, hgvs: bool| -> Option<String> {
        let dir = TempDir::new().unwrap();
        let json = run(dir.path(), "json", columns, hgvs);
        let parsed: serde_json::Value = serde_json::from_str(&json).expect("a JSON array");
        let variant = parsed
            .as_array()
            .unwrap()
            .iter()
            .find(|v| v["seq_region_name"] == "chrM")
            .unwrap_or_else(|| panic!("no chrM variant in:\n{json}"));
        variant["transcript_consequences"][0]
            .get("hgvsg")
            .map(|h| h.as_str().unwrap().to_string())
    };

    assert_eq!(intergenic_hgvsg(ColumnFlags::default(), true), None);
    let asked = ColumnFlags {
        hgvsg: true,
        ..ColumnFlags::default()
    };
    assert_eq!(
        intergenic_hgvsg(asked, false).as_deref(),
        Some("chrM:g.100A>G")
    );
}

/// A `*` allele (a deletion that began upstream) gets no row at all, and the allele beside it
/// has its own `HGVSg`: Ensembl VEP 116.3 writes no row for the `*` and `chr10:g.11800G>T` for
/// the `T`. An earlier fastVEP wrote a row for the `*` with an empty `Consequence`, and before
/// that one with `chr10:g.11800?`, which is not HGVS.
#[test]
fn a_star_allele_has_no_hgvsg_and_the_allele_beside_it_keeps_its_own() {
    let dir = TempDir::new().unwrap();
    let vcf = dir.path().join("star.vcf");
    std::fs::write(
        &vcf,
        "##fileformat=VCFv4.2\n##contig=<ID=chr10,length=133797422>\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
         chr10\t11800\tstar\tG\tT,*\t.\t.\t.\n\
         chr10\t11900\tn_alt\tG\tN\t.\t.\t.\n",
    )
    .unwrap();
    write_gff3(dir.path());
    let columns = ColumnFlags {
        hgvsg: true,
        ..ColumnFlags::default()
    };
    let tab = run_vcf(dir.path(), &vcf, "tab", columns, true);
    let rows = rows(&tab);
    let hgvsg = |uploaded: &str, allele: &str| -> String {
        rows.iter()
            .find(|r| r["Uploaded_variation"] == uploaded && r["Allele"] == allele)
            .unwrap_or_else(|| panic!("no row {uploaded} {allele} in {rows:#?}"))["HGVSg"]
            .clone()
    };
    assert_eq!(hgvsg("star", "T"), "chr10:g.11800G>T");
    assert!(
        !rows
            .iter()
            .any(|r| r["Uploaded_variation"] == "star" && r["Allele"] == "*"),
        "no row is written for a `*` allele: {rows:#?}"
    );
    assert_eq!(hgvsg("n_alt", "N"), "-");
}

#[test]
fn variant_class_is_veps_term_for_each_shape() {
    let dir = TempDir::new().unwrap();
    let flags = ColumnFlags {
        variant_class: true,
        ..ColumnFlags::default()
    };
    let rows = rows(&run(dir.path(), "tab", flags, false));
    assert_eq!(find(&rows, "chr10_11800_G/T")["VARIANT_CLASS"], "SNV");
    assert_eq!(find(&rows, "chr10_11900_GA/G")["VARIANT_CLASS"], "deletion");
    assert_eq!(
        find(&rows, "chr10_11950_C/CTT")["VARIANT_CLASS"],
        "insertion"
    );
    // a variant class belongs to the variant, so the intergenic row has it too
    assert_eq!(find(&rows, "chrM_100_A/G")["VARIANT_CLASS"], "SNV");
}

#[test]
fn the_transcript_columns_are_dashes_on_an_intergenic_row() {
    let dir = TempDir::new().unwrap();
    let flags = ColumnFlags {
        show_ref_allele: true,
        symbol: true,
        biotype: true,
        canonical: true,
        numbers: true,
        ..ColumnFlags::default()
    };
    let rows = rows(&run(dir.path(), "tab", flags, false));

    let exonic = find(&rows, "chr10_11800_G/T");
    assert_eq!(exonic["REF_ALLELE"], "G");
    assert_eq!(exonic["SYMBOL"], "GENEX");
    assert_eq!(exonic["BIOTYPE"], "lncRNA");
    assert_eq!(exonic["CANONICAL"], "YES");
    assert_eq!(exonic["EXON"], "1/2");
    assert_eq!(exonic["INTRON"], "-");

    let intronic = find(&rows, "chr10_12500_A/C");
    assert_eq!(intronic["EXON"], "-");
    assert_eq!(intronic["INTRON"], "1/1");

    // an insertion has `-` as its reference allele
    assert_eq!(find(&rows, "chr10_11950_C/CTT")["REF_ALLELE"], "-");

    let intergenic = find(&rows, "chrM_100_A/G");
    assert_eq!(intergenic["REF_ALLELE"], "A");
    for col in ["SYMBOL", "BIOTYPE", "CANONICAL", "EXON", "INTRON"] {
        assert_eq!(intergenic[col], "-", "{col} of an intergenic row");
    }
}

/// VEP's FLAGS column holds the transcript's flags (`cds_start_NF`,
/// `cds_end_NF`); canonical status is its own column. The tab writer put the word
/// `canonical` in FLAGS.
#[test]
fn flags_does_not_carry_the_canonical_marker() {
    let dir = TempDir::new().unwrap();
    let rows = rows(&run(dir.path(), "tab", ColumnFlags::default(), false));
    for row in &rows {
        assert_ne!(row["FLAGS"], "canonical", "{row:?}");
    }
}

/// `--everything` is VEP's set (Config.pm `OPTION_SETS`), restricted to the
/// columns fastVEP can write: it turns on HGVS, symbol, numbers, canonical,
/// biotype and variant class, and not `--hgvsg` or `--show_ref_allele`.
#[test]
fn everything_is_vep_s_set_of_the_columns_fastvep_writes() {
    let dir = TempDir::new().unwrap();
    let flags = ColumnFlags::everything();
    let tab = run(dir.path(), "tab", flags, true);
    let names = header(&tab);
    for col in [
        "VARIANT_CLASS",
        "SYMBOL",
        "BIOTYPE",
        "CANONICAL",
        "EXON",
        "INTRON",
        "HGVSc",
    ] {
        assert!(names.contains(&col.to_string()), "{col} in {names:?}");
    }
    for col in ["HGVSg", "REF_ALLELE"] {
        assert!(!names.contains(&col.to_string()), "{col} in {names:?}");
    }
}

/// The VCF writer carries `HGVSg` and `VARIANT_CLASS` only when asked, so the
/// 49-field CSQ every existing reader parses does not change.
#[test]
fn csq_gains_hgvsg_and_variant_class_only_when_asked() {
    let csq_fields = |vcf: &str| -> Vec<String> {
        vcf.lines()
            .find(|l| l.starts_with("##INFO=<ID=CSQ"))
            .and_then(|l| l.split("Format: ").nth(1))
            .and_then(|rest| rest.split("\">").next())
            .map(|f| f.split('|').map(String::from).collect())
            .expect("a CSQ header")
    };

    let dir = TempDir::new().unwrap();
    let plain = csq_fields(&run(dir.path(), "vcf", ColumnFlags::default(), false));
    assert_eq!(plain.len(), 49);
    assert!(!plain.contains(&"HGVSg".to_string()));
    assert!(!plain.contains(&"VARIANT_CLASS".to_string()));

    let dir = TempDir::new().unwrap();
    let flags = ColumnFlags {
        hgvsg: true,
        variant_class: true,
        ..ColumnFlags::default()
    };
    let asked = csq_fields(&run(dir.path(), "vcf", flags, false));
    assert_eq!(asked.len(), 51);
    // every field of the plain CSQ is still there, in order
    let without_new: Vec<&String> = asked
        .iter()
        .filter(|f| *f != "HGVSg" && *f != "VARIANT_CLASS")
        .collect();
    assert_eq!(without_new, plain.iter().collect::<Vec<_>>());
    // HGVSg follows HGVS_OFFSET, as in VEP
    let at = |name: &str| asked.iter().position(|f| f == name).unwrap();
    assert_eq!(at("HGVSg"), at("HGVS_OFFSET") + 1);
}

/// What the tab and the VCF writers say for the same column must be the same
/// value. Compared by name over every row both write: a tab `-` is a VCF empty.
#[test]
fn tab_and_csq_agree_on_every_shared_column() {
    let flags = ColumnFlags {
        show_ref_allele: true,
        variant_class: true,
        symbol: true,
        biotype: true,
        canonical: true,
        numbers: true,
        hgvsg: true,
    };
    let dir = TempDir::new().unwrap();
    let tab = rows(&run(dir.path(), "tab", flags, true));
    let dir = TempDir::new().unwrap();
    let vcf = run(dir.path(), "vcf", flags, true);

    let fields: Vec<String> = vcf
        .lines()
        .find(|l| l.starts_with("##INFO=<ID=CSQ"))
        .and_then(|l| l.split("Format: ").nth(1))
        .and_then(|rest| rest.split("\">").next())
        .map(|f| f.split('|').map(String::from).collect())
        .expect("a CSQ header");

    // (Location, Allele, Feature) -> CSQ field -> value, over every entry
    let mut csq: HashMap<(String, String, String), HashMap<String, String>> = HashMap::new();
    for record in vcf.lines().filter(|l| !l.starts_with('#')) {
        let cols: Vec<&str> = record.split('\t').collect();
        let info = cols[7];
        let Some(value) = info.split(';').find_map(|kv| kv.strip_prefix("CSQ=")) else {
            continue;
        };
        for entry in value.split(',') {
            let by_name: HashMap<String, String> = fields
                .iter()
                .cloned()
                .zip(entry.split('|').map(String::from))
                .collect();
            let location = format!("{}:{}", cols[0], cols[1]);
            let key = (
                location,
                by_name["Allele"].clone(),
                by_name["Feature"].clone(),
            );
            csq.insert(key, by_name);
        }
    }

    let shared = [
        "REF_ALLELE",
        "VARIANT_CLASS",
        "SYMBOL",
        "SYMBOL_SOURCE",
        "HGNC_ID",
        "BIOTYPE",
        "CANONICAL",
        "EXON",
        "INTRON",
        "HGVSc",
        "HGVSp",
        "HGVS_OFFSET",
        "HGVSg",
    ];
    let mut compared = 0;
    for row in &tab {
        // the tab Location of an SNV is `chr:pos`; the VCF record sits at its own POS
        let (chrom, pos) = row["Location"].split_once(':').unwrap();
        let pos = pos.split('-').next().unwrap();
        // an indel's VCF POS is the anchor base, one before the minimised start
        let candidates = [
            pos.to_string(),
            (pos.parse::<u64>().unwrap() - 1).to_string(),
        ];
        let feature = if row["Feature"] == "-" {
            String::new()
        } else {
            row["Feature"].clone()
        };
        let entry = candidates.iter().find_map(|p| {
            csq.get(&(
                format!("{chrom}:{p}"),
                row["Allele"].clone(),
                feature.clone(),
            ))
        });
        let Some(entry) = entry else { continue };
        for col in shared {
            let in_csq = entry[col].replace('&', ",");
            // A tab `-` is a value the row does not have, which the VCF leaves
            // empty. It is also a value: the reference allele of an insertion is
            // `-` in both.
            let agree = row[col] == in_csq || (row[col] == "-" && in_csq.is_empty());
            assert!(
                agree,
                "{col} of {} differs between the writers: tab {:?}, VCF {:?}",
                row["Uploaded_variation"], row[col], in_csq
            );
        }
        compared += 1;
    }
    assert!(compared >= 4, "only {compared} rows could be compared");
}
