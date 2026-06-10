//! End-to-end tests for supplementary annotation CLI paths.
//!
//! Each `sa-build` test writes a small fixture for the source, calls
//! `run_sa_build` (the same entrypoint the CLI uses), and reads the resulting
//! database back to confirm the round-trip.

use fastvep_cli::pipeline::{run_annotate, run_sa_build, AnnotateConfig};
use fastvep_sa::gene::GeneIndex;
use std::fs::{self, File};
use std::io::Write;

const SPLICEAI_SOURCE_VCF: &str = include_str!("../fixtures/spliceai/spliceai-mini.vcf");
const SPLICEAI_INDEL_SOURCE_VCF: &str =
    include_str!("../fixtures/spliceai/spliceai-indel-mini.vcf");
const GNOMAD_SOURCE_VCF: &str = include_str!("../fixtures/spliceai/gnomad-mini.vcf");
const INPUT_NO_SPLICEAI_INFO_VCF: &str =
    include_str!("../fixtures/spliceai/input-no-spliceai-info.vcf");
const MINI_GFF3: &str = include_str!("../fixtures/spliceai/mini.gff3");

fn read_oga(path: &std::path::Path) -> GeneIndex {
    let mut f = File::open(path).expect("open .oga");
    GeneIndex::read_from(&mut f).expect("parse .oga")
}

#[test]
fn sa_build_omim_writes_oga_with_records() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("genemap2.txt");
    let output = tmp.path().join("omim");

    // Minimal genemap2.txt fixture — column layout matches the real OMIM
    // file format (gene symbol at index 5, MIM at index 8, phenotypes at 12).
    let fixture = "# Generated\n\
                   # Copyright OMIM\n\
                   1\tp36.33\t1:10001-20000\tGene1\t\tBRCA1\tprotein\t\t113705\t\t\t\tBreast cancer, 114480 (3), Autosomal dominant; Ovarian cancer, 167000 (3)\n\
                   1\tp36.33\t1:30001-40000\tGene2\t\tTP53\tprotein\t\t191170\t\t\t\tLi-Fraumeni syndrome 1, 151623 (3), Autosomal dominant\n";
    File::create(&input)
        .unwrap()
        .write_all(fixture.as_bytes())
        .unwrap();

    run_sa_build(
        "omim",
        input.to_str().unwrap(),
        output.to_str().unwrap(),
        "GRCh38",
        None,
        &[],
    )
    .unwrap();

    let oga_path = output.with_extension("oga");
    assert!(oga_path.exists(), ".oga file should be written");

    let index = read_oga(&oga_path);
    assert_eq!(index.header.json_key, "omim");
    assert_eq!(index.header.name, "ClinGen GDV / OMIM");
    assert!(index.gene_count() >= 2);

    let brca1 = index.get("BRCA1").expect("BRCA1 should be present");
    let json = brca1.first().unwrap();
    assert!(
        json.contains("113705"),
        "BRCA1 should have MIM 113705: {}",
        json
    );
    assert!(json.contains("Breast cancer"));
}

#[test]
fn sa_build_gnomad_genes_writes_oga_with_records() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("constraint.tsv");
    let output = tmp.path().join("gnomad_genes");

    let fixture = "\
gene\ttranscript\tobs_lof\texp_lof\toe_lof\toe_lof_upper\tpLI\tmis_z\tsyn_z
BRCA1\tENST00000357654\t0\t50.2\t0.00\t0.03\t1.0000\t3.45\t0.12
TP53\tENST00000269305\t0\t25.1\t0.00\t0.05\t0.9999\t5.67\t-0.34
";
    File::create(&input)
        .unwrap()
        .write_all(fixture.as_bytes())
        .unwrap();

    run_sa_build(
        "gnomad_genes",
        input.to_str().unwrap(),
        output.to_str().unwrap(),
        "GRCh38",
        None,
        &[],
    )
    .unwrap();

    let oga_path = output.with_extension("oga");
    let index = read_oga(&oga_path);
    assert_eq!(index.header.json_key, "gnomad_genes");
    assert_eq!(index.gene_count(), 2);

    let brca1 = index.get("BRCA1").unwrap();
    let json = brca1.first().unwrap();
    assert!(json.contains("\"pLI\":1.0000"));
    assert!(json.contains("\"loeuf\":0.0300"));
    assert!(json.contains("\"misZ\":3.45"));
}

#[test]
fn sa_build_gnomad_gene_alias_routes_to_same_key() {
    // Both `gnomad_genes` (plural) and `gnomad_gene` (singular) are accepted
    // for the CLI; both must produce a database with json_key=gnomad_genes
    // so the classifier picks them up consistently.
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("constraint.tsv");
    let output = tmp.path().join("gnomad_gene_alias");
    let fixture = "\
gene\ttranscript\tobs_lof\texp_lof\toe_lof\toe_lof_upper\tpLI\tmis_z\tsyn_z
BRCA1\tENST00000357654\t0\t50.2\t0.00\t0.03\t1.0000\t3.45\t0.12
";
    File::create(&input)
        .unwrap()
        .write_all(fixture.as_bytes())
        .unwrap();

    run_sa_build(
        "gnomad_gene",
        input.to_str().unwrap(),
        output.to_str().unwrap(),
        "GRCh38",
        None,
        &[],
    )
    .unwrap();

    let index = read_oga(&output.with_extension("oga"));
    assert_eq!(index.header.json_key, "gnomad_genes");
}

#[test]
fn sa_build_clinvar_protein_writes_oga_with_records() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("clinvar.vcf");
    let output = tmp.path().join("clinvar_protein");

    // Minimal ClinVar VCF: two pathogenic missense entries with protein
    // change in CLNHGVS, one rejected (Benign).
    let fixture = "\
##fileformat=VCFv4.1
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
17\t7676154\t12345\tG\tA\t.\t.\tCLNSIG=Pathogenic;MC=SO:0001583|missense_variant;GENEINFO=TP53:7157;CLNHGVS=NP_000537.3:p.Arg175His
17\t7676156\t12346\tT\tC\t.\t.\tCLNSIG=Likely_pathogenic;MC=SO:0001583|missense_variant;GENEINFO=TP53:7157;CLNHGVS=NP_000537.3:p.Arg248Trp
17\t7676160\t12347\tG\tA\t.\t.\tCLNSIG=Benign;MC=SO:0001583|missense_variant;GENEINFO=TP53:7157;CLNHGVS=NP_000537.3:p.Pro72Leu
";
    File::create(&input)
        .unwrap()
        .write_all(fixture.as_bytes())
        .unwrap();

    run_sa_build(
        "clinvar_protein",
        input.to_str().unwrap(),
        output.to_str().unwrap(),
        "GRCh38",
        None,
        &[],
    )
    .unwrap();

    let index = read_oga(&output.with_extension("oga"));
    assert_eq!(index.header.json_key, "clinvar_protein");
    assert!(index.gene_count() >= 1);

    let tp53 = index.get("TP53").expect("TP53 should be present");
    let json = tp53.first().unwrap();
    // The two pathogenic entries should make it through; the Benign one shouldn't.
    assert!(
        json.contains("\"pos\":175"),
        "should include p.Arg175His: {}",
        json
    );
    assert!(
        json.contains("\"pos\":248"),
        "should include p.Arg248Trp: {}",
        json
    );
    assert!(
        !json.contains("\"pos\":72"),
        "Benign p.Pro72Leu must NOT be in index: {}",
        json
    );
}

#[test]
fn sa_build_unknown_source_errors_with_full_supported_list() {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("nope.txt");
    File::create(&input).unwrap().write_all(b"").unwrap();

    let err = run_sa_build(
        "definitely_not_a_source",
        input.to_str().unwrap(),
        "out",
        "GRCh38",
        None,
        &[],
    )
    .expect_err("must error on unknown source");
    let msg = format!("{}", err);
    // Error should list the new gene-level sources alongside variant-level ones.
    assert!(msg.contains("omim"), "error should mention omim: {}", msg);
    assert!(msg.contains("gnomad_genes"));
    assert!(msg.contains("clinvar_protein"));
}

#[test]
fn annotate_vcf_emits_spliceai_from_fastsa() {
    let tmp = tempfile::tempdir().unwrap();
    let spliceai_source = tmp.path().join("spliceai-mini.vcf");
    let input_vcf = tmp.path().join("input-no-spliceai-info.vcf");
    let gff3 = tmp.path().join("mini.gff3");
    let output_base = tmp.path().join("spliceai-mini");
    let output_vcf = tmp.path().join("annotated.vcf");
    let transcript_cache = tmp.path().join("mini.fastvep.cache");

    fs::write(&spliceai_source, SPLICEAI_SOURCE_VCF).unwrap();
    fs::write(&input_vcf, INPUT_NO_SPLICEAI_INFO_VCF).unwrap();
    fs::write(&gff3, MINI_GFF3).unwrap();

    run_sa_build(
        "spliceai",
        spliceai_source.to_str().unwrap(),
        output_base.to_str().unwrap(),
        "GRCh38",
        None,
        &[],
    )
    .unwrap();

    run_annotate(AnnotateConfig {
        input: input_vcf.to_string_lossy().into_owned(),
        output: output_vcf.to_string_lossy().into_owned(),
        gff3: vec![gff3.to_string_lossy().into_owned()],
        fasta: None,
        output_format: "vcf".into(),
        pick: false,
        hgvs: false,
        distance: 0,
        cache_dir: None,
        transcript_cache: Some(transcript_cache.to_string_lossy().into_owned()),
        sa_dir: Some(tmp.path().to_string_lossy().into_owned()),
        sa_only: false,
        acmg: false,
        acmg_config: None,
        proband: None,
        mother: None,
        father: None,
        gene_list: None,
        explicit_alleles: false,
        qc_rules: None,
    })
    .unwrap();

    let annotated = fs::read_to_string(output_vcf).unwrap();

    assert!(
        annotated.contains("##INFO=<ID=SpliceAI,Number=.,Type=String,Description=\"SpliceAI annotations. Format: ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL\">"),
        "VCF output should declare the SpliceAI INFO header:\n{}",
        annotated
    );
    assert!(
        annotated.contains("SpliceAI=G|GENE1|0.01|0.00|0.85|0.00|5|-28|2|-13"),
        "VCF output should emit SpliceAI from fastSA for single-alt records:\n{}",
        annotated
    );
    assert!(
        annotated.contains(
            "SpliceAI=T|GENE2|0.00|0.10|0.00|0.92|3|-5|10|-2,A|GENE2|0.50|0.00|0.00|0.00|7|0|0|0"
        ),
        "VCF output should emit one SpliceAI value per matching alternate allele:\n{}",
        annotated
    );
}

#[test]
fn annotate_vcf_emits_spliceai_for_uploaded_indel_alleles() {
    let tmp = tempfile::tempdir().unwrap();
    let spliceai_source = tmp.path().join("spliceai-indel-mini.vcf");
    let input_vcf = tmp.path().join("input-no-spliceai-info.vcf");
    let gff3 = tmp.path().join("mini.gff3");
    let output_base = tmp.path().join("spliceai-indel-mini");
    let phylop_source = tmp.path().join("phylop-indel.tsv");
    let phylop_output_base = tmp.path().join("phylop-indel");
    let output_vcf = tmp.path().join("annotated.vcf");
    let transcript_cache = tmp.path().join("mini.fastvep.cache");

    fs::write(&spliceai_source, SPLICEAI_INDEL_SOURCE_VCF).unwrap();
    fs::write(&gff3, MINI_GFF3).unwrap();
    fs::write(
        &phylop_source,
        "\
chr1\t26001\t3.14
chr1\t26011\t2.71
",
    )
    .unwrap();
    fs::write(
        &input_vcf,
        "\
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
1\t26000\t.\tGA\tG\t.\t.\t.
1\t26010\t.\tG\tGA\t.\t.\t.
",
    )
    .unwrap();

    run_sa_build(
        "spliceai",
        spliceai_source.to_str().unwrap(),
        output_base.to_str().unwrap(),
        "GRCh38",
        None,
        &[],
    )
    .unwrap();
    run_sa_build(
        "phylop",
        phylop_source.to_str().unwrap(),
        phylop_output_base.to_str().unwrap(),
        "GRCh38",
        None,
        &[],
    )
    .unwrap();

    run_annotate(AnnotateConfig {
        input: input_vcf.to_string_lossy().into_owned(),
        output: output_vcf.to_string_lossy().into_owned(),
        gff3: vec![gff3.to_string_lossy().into_owned()],
        fasta: None,
        output_format: "vcf".into(),
        pick: false,
        hgvs: false,
        distance: 0,
        cache_dir: None,
        transcript_cache: Some(transcript_cache.to_string_lossy().into_owned()),
        sa_dir: Some(tmp.path().to_string_lossy().into_owned()),
        sa_only: false,
        acmg: false,
        acmg_config: None,
        proband: None,
        mother: None,
        father: None,
        gene_list: None,
        explicit_alleles: false,
        qc_rules: None,
    })
    .unwrap();

    let annotated = fs::read_to_string(output_vcf).unwrap();
    assert!(
        annotated.contains("SpliceAI=G|GENE1|0.10|0.00|0.00|0.00|4|7|27|17"),
        "deletion should use uploaded ALT allele from SpliceAI source:\n{}",
        annotated
    );
    assert!(
        annotated.contains("SpliceAI=GA|GENE1|0.00|0.20|0.00|0.00|1|7|27|17"),
        "insertion should use uploaded ALT allele from SpliceAI source:\n{}",
        annotated
    );
    assert!(
        annotated.contains("FV_PHYLOP=G|3.14"),
        "positional scores should still query fastVEP's normalized indel position:\n{}",
        annotated
    );
    assert!(!annotated.contains("SpliceAI=-|"), "{annotated}");
}

#[test]
fn annotate_vcf_replaces_existing_fastvep_info() {
    let tmp = tempfile::tempdir().unwrap();
    let spliceai_source = tmp.path().join("spliceai-mini.vcf");
    let input_vcf = tmp.path().join("input-existing-info.vcf");
    let gff3 = tmp.path().join("mini.gff3");
    let output_base = tmp.path().join("spliceai-mini");
    let output_vcf = tmp.path().join("annotated.vcf");
    let transcript_cache = tmp.path().join("mini.fastvep.cache");

    fs::write(&spliceai_source, SPLICEAI_SOURCE_VCF).unwrap();
    fs::write(&gff3, MINI_GFF3).unwrap();
    fs::write(
        &input_vcf,
        "\
##fileformat=VCFv4.2
##INFO=<ID=CSQ,Number=.,Type=String,Description=\"stale CSQ\">
##INFO=<ID=SpliceAI,Number=.,Type=String,Description=\"stale SpliceAI\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
1\t25000\t.\tA\tG\t.\t.\tDP=12;CSQ=old;SpliceAI=old
",
    )
    .unwrap();

    run_sa_build(
        "spliceai",
        spliceai_source.to_str().unwrap(),
        output_base.to_str().unwrap(),
        "GRCh38",
        None,
        &[],
    )
    .unwrap();

    run_annotate(AnnotateConfig {
        input: input_vcf.to_string_lossy().into_owned(),
        output: output_vcf.to_string_lossy().into_owned(),
        gff3: vec![gff3.to_string_lossy().into_owned()],
        fasta: None,
        output_format: "vcf".into(),
        pick: false,
        hgvs: false,
        distance: 0,
        cache_dir: None,
        transcript_cache: Some(transcript_cache.to_string_lossy().into_owned()),
        sa_dir: Some(tmp.path().to_string_lossy().into_owned()),
        sa_only: false,
        acmg: false,
        acmg_config: None,
        proband: None,
        mother: None,
        father: None,
        gene_list: None,
        explicit_alleles: false,
        qc_rules: None,
    })
    .unwrap();

    let annotated = fs::read_to_string(output_vcf).unwrap();
    assert_eq!(annotated.matches("##INFO=<ID=CSQ,").count(), 1, "{annotated}");
    assert_eq!(
        annotated.matches("##INFO=<ID=SpliceAI,").count(),
        1,
        "{annotated}"
    );
    assert!(annotated.contains("DP=12;SpliceAI=G|GENE1|0.01|0.00|0.85|0.00|5|-28|2|-13;CSQ=G|"));
    assert!(!annotated.contains("CSQ=old"), "{annotated}");
    assert!(!annotated.contains("SpliceAI=old"), "{annotated}");
    assert!(!annotated.contains("stale CSQ"), "{annotated}");
    assert!(!annotated.contains("stale SpliceAI"), "{annotated}");
}

#[test]
fn annotate_vcf_emits_fastsa_projection_for_gnomad() {
    let tmp = tempfile::tempdir().unwrap();
    let gnomad_source = tmp.path().join("gnomad-mini.vcf");
    let input_vcf = tmp.path().join("input.vcf");
    let gff3 = tmp.path().join("mini.gff3");
    let output_base = tmp.path().join("gnomad-mini");
    let output_vcf = tmp.path().join("annotated.vcf");
    let transcript_cache = tmp.path().join("mini.fastvep.cache");

    fs::write(&gnomad_source, GNOMAD_SOURCE_VCF).unwrap();
    fs::write(&input_vcf, INPUT_NO_SPLICEAI_INFO_VCF).unwrap();
    fs::write(&gff3, MINI_GFF3).unwrap();

    run_sa_build(
        "gnomad",
        gnomad_source.to_str().unwrap(),
        output_base.to_str().unwrap(),
        "GRCh38",
        None,
        &[],
    )
    .unwrap();

    run_annotate(AnnotateConfig {
        input: input_vcf.to_string_lossy().into_owned(),
        output: output_vcf.to_string_lossy().into_owned(),
        gff3: vec![gff3.to_string_lossy().into_owned()],
        fasta: None,
        output_format: "vcf".into(),
        pick: false,
        hgvs: false,
        distance: 0,
        cache_dir: None,
        transcript_cache: Some(transcript_cache.to_string_lossy().into_owned()),
        sa_dir: Some(tmp.path().to_string_lossy().into_owned()),
        sa_only: false,
        acmg: false,
        acmg_config: None,
        proband: None,
        mother: None,
        father: None,
        gene_list: None,
        explicit_alleles: false,
        qc_rules: None,
    })
    .unwrap();

    let annotated = fs::read_to_string(output_vcf).unwrap();
    assert!(
        annotated.contains("##INFO=<ID=FV_GNOMAD,Number=.,Type=String"),
        "{annotated}"
    );
    assert!(
        annotated.contains("FV_GNOMAD=G|0.00012|12|100000|0|0.00021"),
        "{annotated}"
    );
    assert!(!annotated.contains("FV_GNOMAD={"), "{annotated}");
}

#[test]
fn annotate_tab_emits_fastsa_columns_for_clinvar_and_gnomad() {
    // Regression test for issue #30: tab output silently dropped every
    // supplementary annotation source. After fix, each loaded source must
    // produce one extra tab column with the same pipe-delimited value used
    // for the VCF `FV_*` INFO field.
    let tmp = tempfile::tempdir().unwrap();
    let gnomad_source = tmp.path().join("gnomad-mini.vcf");
    let input_vcf = tmp.path().join("input.vcf");
    let gff3 = tmp.path().join("mini.gff3");
    let gnomad_base = tmp.path().join("gnomad-mini");
    let output_tab = tmp.path().join("annotated.tab");
    let transcript_cache = tmp.path().join("mini.fastvep.cache");

    write_clinvar_fixture(tmp.path());
    fs::write(&gnomad_source, GNOMAD_SOURCE_VCF).unwrap();
    fs::write(&input_vcf, INPUT_NO_SPLICEAI_INFO_VCF).unwrap();
    fs::write(&gff3, MINI_GFF3).unwrap();

    run_sa_build(
        "gnomad",
        gnomad_source.to_str().unwrap(),
        gnomad_base.to_str().unwrap(),
        "GRCh38",
        None,
        &[],
    )
    .unwrap();

    run_annotate(AnnotateConfig {
        input: input_vcf.to_string_lossy().into_owned(),
        output: output_tab.to_string_lossy().into_owned(),
        gff3: vec![gff3.to_string_lossy().into_owned()],
        fasta: None,
        output_format: "tab".into(),
        pick: false,
        hgvs: false,
        distance: 0,
        cache_dir: None,
        transcript_cache: Some(transcript_cache.to_string_lossy().into_owned()),
        sa_dir: Some(tmp.path().to_string_lossy().into_owned()),
        sa_only: false,
        acmg: false,
        acmg_config: None,
        proband: None,
        mother: None,
        father: None,
        gene_list: None,
        explicit_alleles: false,
        qc_rules: None,
    })
    .unwrap();

    let annotated = fs::read_to_string(output_tab).unwrap();

    // Tab schema header must document the pipe format for every loaded source.
    assert!(
        annotated.contains("## COLUMN=<ID=FV_CLINVAR,Description=\"fastVEP ClinVar annotations. Format: ALLELE|SIGNIFICANCE|REVIEW_STATUS|PHENOTYPES|VARIANT_CLASS|SO_ACCESSION\">"),
        "tab output must declare FV_CLINVAR schema:\n{}",
        annotated
    );
    assert!(
        annotated.contains("## COLUMN=<ID=FV_GNOMAD,Description="),
        "tab output must declare FV_GNOMAD schema:\n{}",
        annotated
    );

    // The column header line must end with the FV_* columns.
    let column_header = annotated
        .lines()
        .find(|l| l.starts_with("#Uploaded_variation"))
        .expect("missing tab column header");
    assert!(
        column_header.ends_with("\tFV_CLINVAR\tFV_GNOMAD"),
        "extra columns must append after FLAGS in spec order: {}",
        column_header
    );

    // At least one data row carries a populated FV_CLINVAR value at position 25000.
    let data_rows: Vec<&str> = annotated
        .lines()
        .filter(|l| !l.starts_with('#'))
        .collect();
    assert!(!data_rows.is_empty(), "tab output must contain data rows");
    let pos25k = data_rows
        .iter()
        .find(|r| r.contains("1:25000\t"))
        .expect("expected an annotated row for chr1:25000");
    let cols: Vec<&str> = pos25k.split('\t').collect();
    assert_eq!(cols.len(), 17 + 2);
    assert_eq!(
        cols[17], "G|Pathogenic|criteria_provided%2C_multiple_submitters%2C_no_conflicts|Breast_cancer|SNV|SO%3A0001483",
        "FV_CLINVAR tab column must match the VCF pipe schema: {}",
        pos25k
    );
    // The fixture only populates a subset of gnomAD fields; absent positions
    // render as empty between pipes (matches VCF behavior).
    assert_eq!(
        cols[18], "G|0.00012|12|100000|0|0.00021||||||0.00009|||",
        "FV_GNOMAD tab column must match the VCF pipe schema: {}",
        pos25k
    );

    // Sanity: no raw JSON leaked into the tab file.
    assert!(!annotated.contains('{'), "tab output must not contain raw JSON:\n{}", annotated);
    assert!(!annotated.contains('}'), "tab output must not contain raw JSON:\n{}", annotated);
}

/// Build a minimal ClinVar SA database in a temp dir and return its path.
/// Used by the --sa-only tests below.
fn write_clinvar_fixture(tmp: &std::path::Path) {
    let clinvar_source = tmp.join("clinvar-mini.vcf");
    let clinvar_base = tmp.join("clinvar-mini");
    let clinvar_fixture = "\
##fileformat=VCFv4.1
##INFO=<ID=CLNSIG,Number=.,Type=String>
##INFO=<ID=CLNREVSTAT,Number=.,Type=String>
##INFO=<ID=CLNDN,Number=.,Type=String>
##INFO=<ID=CLNVC,Number=.,Type=String>
##INFO=<ID=CLNVCSO,Number=.,Type=String>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
1\t25000\trs1\tA\tG\t.\t.\tCLNSIG=Pathogenic;CLNREVSTAT=criteria_provided,_multiple_submitters,_no_conflicts;CLNDN=Breast_cancer;CLNVC=SNV;CLNVCSO=SO:0001483
";
    fs::write(&clinvar_source, clinvar_fixture).unwrap();
    run_sa_build(
        "clinvar",
        clinvar_source.to_str().unwrap(),
        clinvar_base.to_str().unwrap(),
        "GRCh38",
        None,
        &[],
    )
    .unwrap();
}

#[test]
fn sa_only_vcf_omits_csq_and_default_pipeline() {
    // --sa-only must skip the default 49-field CSQ annotation entirely:
    // no ##INFO=<ID=CSQ> header, no CSQ= INFO field on any data row.
    // FV_CLINVAR (and any other --sa-dir-loaded source) still emits.
    let tmp = tempfile::tempdir().unwrap();
    let input_vcf = tmp.path().join("input.vcf");
    let output_vcf = tmp.path().join("annotated.vcf");
    fs::write(&input_vcf, INPUT_NO_SPLICEAI_INFO_VCF).unwrap();
    write_clinvar_fixture(tmp.path());

    run_annotate(AnnotateConfig {
        input: input_vcf.to_string_lossy().into_owned(),
        output: output_vcf.to_string_lossy().into_owned(),
        gff3: vec![],
        fasta: None,
        output_format: "vcf".into(),
        pick: false,
        hgvs: false,
        distance: 0,
        cache_dir: None,
        transcript_cache: None,
        sa_dir: Some(tmp.path().to_string_lossy().into_owned()),
        sa_only: true,
        acmg: false,
        acmg_config: None,
        proband: None,
        mother: None,
        father: None,
        gene_list: None,
        explicit_alleles: false,
        qc_rules: None,
    })
    .unwrap();

    let annotated = fs::read_to_string(&output_vcf).unwrap();

    assert!(
        !annotated.contains("##INFO=<ID=CSQ"),
        "--sa-only must not emit the CSQ INFO header:\n{}",
        annotated
    );
    assert!(
        annotated.contains("##INFO=<ID=FV_CLINVAR,"),
        "--sa-only must still emit FV_CLINVAR header:\n{}",
        annotated
    );
    for row in annotated.lines().filter(|l| !l.starts_with('#')) {
        assert!(
            !row.contains("CSQ="),
            "--sa-only data row must not contain CSQ=: {}",
            row
        );
    }
    assert!(
        annotated.contains("FV_CLINVAR=G|Pathogenic"),
        "--sa-only must emit FV_CLINVAR on chr1:25000 data row:\n{}",
        annotated
    );
}

#[test]
fn sa_only_tab_emits_minimal_columns() {
    // --sa-only tab layout: header is Uploaded_variation, Location, Allele
    // followed by one column per loaded SA source (FV_CLINVAR here).
    let tmp = tempfile::tempdir().unwrap();
    let input_vcf = tmp.path().join("input.vcf");
    let output_tab = tmp.path().join("annotated.tab");
    fs::write(&input_vcf, INPUT_NO_SPLICEAI_INFO_VCF).unwrap();
    write_clinvar_fixture(tmp.path());

    run_annotate(AnnotateConfig {
        input: input_vcf.to_string_lossy().into_owned(),
        output: output_tab.to_string_lossy().into_owned(),
        gff3: vec![],
        fasta: None,
        output_format: "tab".into(),
        pick: false,
        hgvs: false,
        distance: 0,
        cache_dir: None,
        transcript_cache: None,
        sa_dir: Some(tmp.path().to_string_lossy().into_owned()),
        sa_only: true,
        acmg: false,
        acmg_config: None,
        proband: None,
        mother: None,
        father: None,
        gene_list: None,
        explicit_alleles: false,
        qc_rules: None,
    })
    .unwrap();

    let annotated = fs::read_to_string(&output_tab).unwrap();
    let column_header = annotated
        .lines()
        .find(|l| l.starts_with("#Uploaded_variation"))
        .expect("missing tab column header");
    assert_eq!(
        column_header, "#Uploaded_variation\tLocation\tAllele\tFV_CLINVAR",
        "sa-only tab header must be minimal: {}",
        column_header
    );

    let data_rows: Vec<&str> = annotated.lines().filter(|l| !l.starts_with('#')).collect();
    assert!(!data_rows.is_empty(), "expected sa-only tab data rows");
    for row in &data_rows {
        let cols: Vec<&str> = row.split('\t').collect();
        assert_eq!(cols.len(), 4, "sa-only tab row must have 4 columns: {}", row);
    }
    let pos25k = data_rows
        .iter()
        .find(|r| r.contains("1:25000\t"))
        .expect("expected a row at 1:25000");
    let cols: Vec<&str> = pos25k.split('\t').collect();
    assert_eq!(cols[2], "G");
    assert_eq!(
        cols[3],
        "G|Pathogenic|criteria_provided%2C_multiple_submitters%2C_no_conflicts|Breast_cancer|SNV|SO%3A0001483"
    );
}

#[test]
fn sa_only_json_omits_transcript_consequences() {
    // --sa-only JSON: top-level "alleles" array carries per-allele SA payloads;
    // no transcript_consequences / most_severe_consequence keys.
    let tmp = tempfile::tempdir().unwrap();
    let input_vcf = tmp.path().join("input.vcf");
    let output_json = tmp.path().join("annotated.json");
    fs::write(&input_vcf, INPUT_NO_SPLICEAI_INFO_VCF).unwrap();
    write_clinvar_fixture(tmp.path());

    run_annotate(AnnotateConfig {
        input: input_vcf.to_string_lossy().into_owned(),
        output: output_json.to_string_lossy().into_owned(),
        gff3: vec![],
        fasta: None,
        output_format: "json".into(),
        pick: false,
        hgvs: false,
        distance: 0,
        cache_dir: None,
        transcript_cache: None,
        sa_dir: Some(tmp.path().to_string_lossy().into_owned()),
        sa_only: true,
        acmg: false,
        acmg_config: None,
        proband: None,
        mother: None,
        father: None,
        gene_list: None,
        explicit_alleles: false,
        qc_rules: None,
    })
    .unwrap();

    let annotated = fs::read_to_string(&output_json).unwrap();
    let parsed: serde_json::Value = serde_json::from_str(&annotated).expect("valid JSON array");
    let arr = parsed.as_array().expect("top-level JSON must be an array");
    assert!(!arr.is_empty(), "expected at least one variant record");
    for record in arr {
        let obj = record.as_object().expect("record must be an object");
        assert!(
            !obj.contains_key("transcript_consequences"),
            "sa-only JSON must omit transcript_consequences: {}",
            record
        );
        assert!(
            !obj.contains_key("most_severe_consequence"),
            "sa-only JSON must omit most_severe_consequence: {}",
            record
        );
        assert!(
            obj.contains_key("alleles"),
            "sa-only JSON must include alleles array: {}",
            record
        );
    }

    let first = arr[0].as_object().unwrap();
    let alleles = first["alleles"].as_array().unwrap();
    let g = alleles
        .iter()
        .find(|a| a["allele"].as_str() == Some("G"))
        .expect("expected an allele:G entry on chr1:25000");
    let clinvar = g.get("clinvar").expect("clinvar key on allele:G");
    // The fixture's CLNSIG=Pathogenic round-trips as a single-element array
    // of strings. Assert the exact shape so a regression that flips array vs
    // bare-string representation is caught.
    let significance = clinvar["significance"]
        .as_array()
        .expect("clinvar.significance should be an array");
    assert_eq!(
        significance.iter().filter_map(|v| v.as_str()).collect::<Vec<_>>(),
        vec!["Pathogenic"],
        "clinvar.significance should be exactly [\"Pathogenic\"]: {}",
        clinvar
    );
    assert_eq!(
        clinvar["reviewStatus"].as_str(),
        Some("criteria_provided,_multiple_submitters,_no_conflicts"),
        "clinvar.reviewStatus mismatch: {}",
        clinvar
    );
    let phenotypes = clinvar["phenotypes"]
        .as_array()
        .expect("clinvar.phenotypes should be an array");
    assert_eq!(
        phenotypes.iter().filter_map(|v| v.as_str()).collect::<Vec<_>>(),
        vec!["Breast_cancer"],
        "clinvar.phenotypes mismatch: {}",
        clinvar
    );
}

#[test]
fn sa_only_requires_sa_dir() {
    let tmp = tempfile::tempdir().unwrap();
    let input_vcf = tmp.path().join("input.vcf");
    let output_vcf = tmp.path().join("annotated.vcf");
    fs::write(&input_vcf, INPUT_NO_SPLICEAI_INFO_VCF).unwrap();

    let err = run_annotate(AnnotateConfig {
        input: input_vcf.to_string_lossy().into_owned(),
        output: output_vcf.to_string_lossy().into_owned(),
        gff3: vec![],
        fasta: None,
        output_format: "vcf".into(),
        pick: false,
        hgvs: false,
        distance: 0,
        cache_dir: None,
        transcript_cache: None,
        sa_dir: None,
        sa_only: true,
        acmg: false,
        acmg_config: None,
        proband: None,
        mother: None,
        father: None,
        gene_list: None,
        explicit_alleles: false,
        qc_rules: None,
    })
    .expect_err("--sa-only without --sa-dir must error");
    assert!(
        err.to_string().contains("--sa-only requires --sa-dir"),
        "error message should mention --sa-dir requirement: {}",
        err
    );
}

#[test]
fn sa_only_multi_allelic_emits_per_alt_rows_with_independent_sa_columns() {
    // Regression: --sa-only mode must lookup supplementary annotations
    // independently for each ALT of a multi-allelic site. The input has
    // 1:30000 C>T,A. A ClinVar fixture is seeded with a match for C>T only.
    // We expect exactly two rows for 1:30000 — one per alt — with
    // FV_CLINVAR populated only on the T row.
    let tmp = tempfile::tempdir().unwrap();
    let input_vcf = tmp.path().join("input.vcf");
    let output_tab = tmp.path().join("annotated.tab");
    fs::write(&input_vcf, INPUT_NO_SPLICEAI_INFO_VCF).unwrap();

    // Custom ClinVar fixture: pathogenic for C>T at chr1:30000, nothing at A.
    let clinvar_source = tmp.path().join("clinvar-mini.vcf");
    let clinvar_base = tmp.path().join("clinvar-mini");
    let clinvar_fixture = "\
##fileformat=VCFv4.1
##INFO=<ID=CLNSIG,Number=.,Type=String>
##INFO=<ID=CLNREVSTAT,Number=.,Type=String>
##INFO=<ID=CLNDN,Number=.,Type=String>
##INFO=<ID=CLNVC,Number=.,Type=String>
##INFO=<ID=CLNVCSO,Number=.,Type=String>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
1\t30000\trs2\tC\tT\t.\t.\tCLNSIG=Likely_pathogenic;CLNREVSTAT=criteria_provided;CLNDN=Test;CLNVC=SNV;CLNVCSO=SO:0001483
";
    fs::write(&clinvar_source, clinvar_fixture).unwrap();
    run_sa_build(
        "clinvar",
        clinvar_source.to_str().unwrap(),
        clinvar_base.to_str().unwrap(),
        "GRCh38",
        None,
        &[],
    )
    .unwrap();

    run_annotate(AnnotateConfig {
        input: input_vcf.to_string_lossy().into_owned(),
        output: output_tab.to_string_lossy().into_owned(),
        gff3: vec![],
        fasta: None,
        output_format: "tab".into(),
        pick: false,
        hgvs: false,
        distance: 0,
        cache_dir: None,
        transcript_cache: None,
        sa_dir: Some(tmp.path().to_string_lossy().into_owned()),
        sa_only: true,
        acmg: false,
        acmg_config: None,
        proband: None,
        mother: None,
        father: None,
        gene_list: None,
        explicit_alleles: false,
        qc_rules: None,
    })
    .unwrap();

    let annotated = fs::read_to_string(&output_tab).unwrap();
    let rows_30k: Vec<&str> = annotated
        .lines()
        .filter(|l| !l.starts_with('#') && l.contains("1:30000\t"))
        .collect();
    assert_eq!(rows_30k.len(), 2, "expected one row per ALT, got: {:?}", rows_30k);

    let mut t_row = None;
    let mut a_row = None;
    for row in rows_30k {
        let cols: Vec<&str> = row.split('\t').collect();
        assert_eq!(cols.len(), 4, "sa-only tab row must be 4 cols: {}", row);
        match cols[2] {
            "T" => t_row = Some(cols[3].to_string()),
            "A" => a_row = Some(cols[3].to_string()),
            other => panic!("unexpected allele: {}", other),
        }
    }
    let t_fv = t_row.expect("missing T row");
    let a_fv = a_row.expect("missing A row");
    assert!(
        t_fv.starts_with("T|Likely_pathogenic"),
        "T row must carry ClinVar match: {}",
        t_fv
    );
    assert_eq!(a_fv, "-", "A row must have no ClinVar match: {}", a_fv);
}

#[test]
fn sa_only_strips_preexisting_csq_from_input_info() {
    // Regression: --sa-only must strip any stale CSQ=... already present in
    // the input VCF's INFO column. Otherwise the output has CSQ= data rows
    // with no matching ##INFO=<ID=CSQ> header (we drop the header in sa_only).
    let tmp = tempfile::tempdir().unwrap();
    let input_vcf = tmp.path().join("input_with_csq.vcf");
    let output_vcf = tmp.path().join("annotated.vcf");
    fs::write(
        &input_vcf,
        "##fileformat=VCFv4.2\n\
         ##INFO=<ID=CSQ,Number=.,Type=String,Description=\"stale\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
         1\t25000\t.\tA\tG\t.\t.\tCSQ=stale_value;DP=10\n",
    )
    .unwrap();
    write_clinvar_fixture(tmp.path());

    run_annotate(AnnotateConfig {
        input: input_vcf.to_string_lossy().into_owned(),
        output: output_vcf.to_string_lossy().into_owned(),
        gff3: vec![],
        fasta: None,
        output_format: "vcf".into(),
        pick: false,
        hgvs: false,
        distance: 0,
        cache_dir: None,
        transcript_cache: None,
        sa_dir: Some(tmp.path().to_string_lossy().into_owned()),
        sa_only: true,
        acmg: false,
        acmg_config: None,
        proband: None,
        mother: None,
        father: None,
        gene_list: None,
        explicit_alleles: false,
        qc_rules: None,
    })
    .unwrap();

    let annotated = fs::read_to_string(&output_vcf).unwrap();
    assert!(
        !annotated.contains("##INFO=<ID=CSQ"),
        "--sa-only must drop the CSQ header from input:\n{}",
        annotated
    );
    let data_row = annotated
        .lines()
        .find(|l| !l.starts_with('#'))
        .expect("expected a data row");
    assert!(
        !data_row.contains("CSQ="),
        "--sa-only must strip stale CSQ= from input INFO: {}",
        data_row
    );
    assert!(
        data_row.contains("DP=10"),
        "non-CSQ INFO fields must pass through: {}",
        data_row
    );
    assert!(
        data_row.contains("FV_CLINVAR="),
        "FV_CLINVAR must still be added: {}",
        data_row
    );
}

#[test]
fn sa_only_strips_csq_when_in_middle_of_info_field() {
    // CSQ stripping must work even when CSQ is sandwiched between other INFO
    // keys. Earlier this was only tested with CSQ at the leading position;
    // make sure the middle case also leaves neighbouring fields intact.
    let tmp = tempfile::tempdir().unwrap();
    let input_vcf = tmp.path().join("input_csq_middle.vcf");
    let output_vcf = tmp.path().join("annotated.vcf");
    fs::write(
        &input_vcf,
        "##fileformat=VCFv4.2\n\
         ##INFO=<ID=CSQ,Number=.,Type=String,Description=\"stale\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
         1\t25000\t.\tA\tG\t.\t.\tAC=1;CSQ=stale_middle;AF=0.5\n",
    )
    .unwrap();
    write_clinvar_fixture(tmp.path());

    run_annotate(AnnotateConfig {
        input: input_vcf.to_string_lossy().into_owned(),
        output: output_vcf.to_string_lossy().into_owned(),
        gff3: vec![],
        fasta: None,
        output_format: "vcf".into(),
        pick: false,
        hgvs: false,
        distance: 0,
        cache_dir: None,
        transcript_cache: None,
        sa_dir: Some(tmp.path().to_string_lossy().into_owned()),
        sa_only: true,
        acmg: false,
        acmg_config: None,
        proband: None,
        mother: None,
        father: None,
        gene_list: None,
        explicit_alleles: false,
        qc_rules: None,
    })
    .unwrap();

    let annotated = fs::read_to_string(&output_vcf).unwrap();
    let row = annotated
        .lines()
        .find(|l| !l.starts_with('#'))
        .expect("expected a data row");
    assert!(!row.contains("CSQ=stale"), "stale CSQ in middle not stripped: {}", row);
    assert!(row.contains("AC=1"), "AC must be preserved: {}", row);
    assert!(row.contains("AF=0.5"), "AF must be preserved: {}", row);
    assert!(row.contains("FV_CLINVAR="), "FV_CLINVAR must still be added: {}", row);
}

#[test]
fn intergenic_variant_with_sa_dir_in_default_mode_emits_fv_clinvar() {
    // Documents a behavior change in this PR: when running in default mode
    // (no --sa-only) with --sa-dir loaded, intergenic variants now also
    // receive supplementary annotations. Previously the SA attachment block
    // was nested inside the `else` of `if overlapping.is_empty()` and so
    // ran only for variants overlapping at least one transcript.
    let tmp = tempfile::tempdir().unwrap();
    let input_vcf = tmp.path().join("input.vcf");
    let output_vcf = tmp.path().join("annotated.vcf");
    let gff3 = tmp.path().join("mini.gff3");
    let transcript_cache = tmp.path().join("mini.fastvep.cache");

    // The fixture transcript covers chr1:25000-30000 region. We pick chr2
    // (no transcripts) to force the intergenic path.
    fs::write(
        &input_vcf,
        "##fileformat=VCFv4.2\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
         2\t25000\t.\tA\tG\t.\t.\t.\n",
    )
    .unwrap();
    fs::write(&gff3, MINI_GFF3).unwrap();

    // ClinVar fixture covers chr2:25000.
    let clinvar_source = tmp.path().join("clinvar-mini.vcf");
    let clinvar_base = tmp.path().join("clinvar-mini");
    fs::write(
        &clinvar_source,
        "##fileformat=VCFv4.1\n\
         ##INFO=<ID=CLNSIG,Number=.,Type=String>\n\
         ##INFO=<ID=CLNREVSTAT,Number=.,Type=String>\n\
         ##INFO=<ID=CLNDN,Number=.,Type=String>\n\
         ##INFO=<ID=CLNVC,Number=.,Type=String>\n\
         ##INFO=<ID=CLNVCSO,Number=.,Type=String>\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
         2\t25000\trs1\tA\tG\t.\t.\tCLNSIG=Pathogenic;CLNREVSTAT=criteria_provided,_single_submitter;CLNDN=Disease;CLNVC=SNV;CLNVCSO=SO:0001483\n",
    )
    .unwrap();
    run_sa_build(
        "clinvar",
        clinvar_source.to_str().unwrap(),
        clinvar_base.to_str().unwrap(),
        "GRCh38",
        None,
        &[],
    )
    .unwrap();

    run_annotate(AnnotateConfig {
        input: input_vcf.to_string_lossy().into_owned(),
        output: output_vcf.to_string_lossy().into_owned(),
        gff3: vec![gff3.to_string_lossy().into_owned()],
        fasta: None,
        output_format: "vcf".into(),
        pick: false,
        hgvs: false,
        distance: 0,
        cache_dir: None,
        transcript_cache: Some(transcript_cache.to_string_lossy().into_owned()),
        sa_dir: Some(tmp.path().to_string_lossy().into_owned()),
        sa_only: false,
        acmg: false,
        acmg_config: None,
        proband: None,
        mother: None,
        father: None,
        gene_list: None,
        explicit_alleles: false,
        qc_rules: None,
    })
    .unwrap();

    let annotated = fs::read_to_string(&output_vcf).unwrap();
    let data_row = annotated
        .lines()
        .find(|l| l.starts_with("2\t25000\t"))
        .expect("expected an intergenic data row");
    assert!(
        data_row.contains("FV_CLINVAR=G|Pathogenic"),
        "intergenic variant should now receive FV_CLINVAR: {}",
        data_row
    );
    // The CSQ field is still emitted with the intergenic_variant consequence.
    assert!(
        data_row.contains("CSQ=G|intergenic_variant|"),
        "intergenic variant should still emit CSQ: {}",
        data_row
    );
}

#[test]
fn annotate_tab_gene_list_filters_to_panel_genes() {
    // Issue #1 ask #4: --gene-list keeps only rows whose transcript matches
    // a gene in the panel. Variants in non-panel genes drop out entirely.
    let tmp = tempfile::tempdir().unwrap();
    let input_vcf = tmp.path().join("input.vcf");
    let gff3 = tmp.path().join("mini.gff3");
    let panel = tmp.path().join("panel.txt");
    let output_tab = tmp.path().join("annotated.tab");
    let transcript_cache = tmp.path().join("mini.fastvep.cache");

    fs::write(&input_vcf, INPUT_NO_SPLICEAI_INFO_VCF).unwrap();
    fs::write(&gff3, MINI_GFF3).unwrap();

    // Panel includes GENE1 (in the GFF3) but not OTHER_GENE.
    fs::write(&panel, "# DDR panel\nGENE1\nOTHER_GENE\n").unwrap();

    run_annotate(AnnotateConfig {
        input: input_vcf.to_string_lossy().into_owned(),
        output: output_tab.to_string_lossy().into_owned(),
        gff3: vec![gff3.to_string_lossy().into_owned()],
        fasta: None,
        output_format: "tab".into(),
        pick: false,
        hgvs: false,
        distance: 0,
        cache_dir: None,
        transcript_cache: Some(transcript_cache.to_string_lossy().into_owned()),
        sa_dir: None,
        sa_only: false,
        acmg: false,
        acmg_config: None,
        proband: None,
        mother: None,
        father: None,
        gene_list: Some(panel.to_string_lossy().into_owned()),
        explicit_alleles: false,
        qc_rules: None,
    })
    .unwrap();

    let annotated = fs::read_to_string(&output_tab).unwrap();
    let data_rows: Vec<&str> = annotated.lines().filter(|l| !l.starts_with('#')).collect();
    assert!(!data_rows.is_empty(), "panel-matching rows should remain");
    for row in &data_rows {
        let cols: Vec<&str> = row.split('\t').collect();
        // Gene id is column 3 (0-based); the fixture uses gene id "GENE1".
        assert_eq!(
            cols[3], "GENE1",
            "every emitted row should belong to a panel gene: {}",
            row
        );
    }
}

#[test]
fn annotate_tab_explicit_alleles_inserts_ref_column() {
    // Issue #1 ask #1: --explicit-alleles adds a REF column right after the
    // Allele column so spreadsheets can read REF/ALT side-by-side.
    let tmp = tempfile::tempdir().unwrap();
    let input_vcf = tmp.path().join("input.vcf");
    let gff3 = tmp.path().join("mini.gff3");
    let output_tab = tmp.path().join("annotated.tab");
    let transcript_cache = tmp.path().join("mini.fastvep.cache");

    fs::write(&input_vcf, INPUT_NO_SPLICEAI_INFO_VCF).unwrap();
    fs::write(&gff3, MINI_GFF3).unwrap();

    run_annotate(AnnotateConfig {
        input: input_vcf.to_string_lossy().into_owned(),
        output: output_tab.to_string_lossy().into_owned(),
        gff3: vec![gff3.to_string_lossy().into_owned()],
        fasta: None,
        output_format: "tab".into(),
        pick: false,
        hgvs: false,
        distance: 0,
        cache_dir: None,
        transcript_cache: Some(transcript_cache.to_string_lossy().into_owned()),
        sa_dir: None,
        sa_only: false,
        acmg: false,
        acmg_config: None,
        proband: None,
        mother: None,
        father: None,
        gene_list: None,
        explicit_alleles: true,
        qc_rules: None,
    })
    .unwrap();

    let annotated = fs::read_to_string(&output_tab).unwrap();
    let column_header = annotated
        .lines()
        .find(|l| l.starts_with("#Uploaded_variation"))
        .expect("missing tab column header");
    assert!(
        column_header.starts_with("#Uploaded_variation\tLocation\tAllele\tREF\tGene"),
        "REF column must come right after Allele: {}",
        column_header
    );

    let pos25k = annotated
        .lines()
        .filter(|l| !l.starts_with('#'))
        .find(|r| r.contains("1:25000\t"))
        .expect("expected a row at 1:25000");
    let cols: Vec<&str> = pos25k.split('\t').collect();
    // Base 17 columns + 1 REF column = 18.
    assert_eq!(cols.len(), 18, "row should carry 18 cols, got: {:?}", cols);
    assert_eq!(cols[2], "G", "Allele column holds ALT");
    assert_eq!(cols[3], "A", "REF column holds REF allele");
}

#[test]
fn annotate_tab_qc_rules_classifies_from_info_fields() {
    // Issue #1 ask #3: --qc-rules adds a QC_CLASS column populated from
    // INFO fields. No per-sample work; rules are evaluated on each row.
    let tmp = tempfile::tempdir().unwrap();
    let input_vcf = tmp.path().join("input.vcf");
    let gff3 = tmp.path().join("mini.gff3");
    let qc_rules_path = tmp.path().join("qc.toml");
    let output_tab = tmp.path().join("annotated.tab");
    let transcript_cache = tmp.path().join("mini.fastvep.cache");

    // Two variants with different INFO/DP values.
    fs::write(
        &input_vcf,
        "##fileformat=VCFv4.2\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
         1\t25000\t.\tA\tG\t.\tPASS\tDP=50;QD=30;MQ=60\n\
         1\t30000\t.\tC\tT\t.\tPASS\tDP=5;QD=2\n",
    )
    .unwrap();
    fs::write(&gff3, MINI_GFF3).unwrap();
    fs::write(
        &qc_rules_path,
        r#"
[[class]]
name = "HIGH_QC"
min_dp = 15
min_qd = 20

[[class]]
name = "LOW_QC"
min_dp = 8
"#,
    )
    .unwrap();

    run_annotate(AnnotateConfig {
        input: input_vcf.to_string_lossy().into_owned(),
        output: output_tab.to_string_lossy().into_owned(),
        gff3: vec![gff3.to_string_lossy().into_owned()],
        fasta: None,
        output_format: "tab".into(),
        pick: false,
        hgvs: false,
        distance: 0,
        cache_dir: None,
        transcript_cache: Some(transcript_cache.to_string_lossy().into_owned()),
        sa_dir: None,
        sa_only: false,
        acmg: false,
        acmg_config: None,
        proband: None,
        mother: None,
        father: None,
        gene_list: None,
        explicit_alleles: false,
        qc_rules: Some(qc_rules_path.to_string_lossy().into_owned()),
    })
    .unwrap();

    let annotated = fs::read_to_string(&output_tab).unwrap();
    let column_header = annotated
        .lines()
        .find(|l| l.starts_with("#Uploaded_variation"))
        .expect("missing tab column header");
    assert!(
        column_header.ends_with("\tQC_CLASS"),
        "QC_CLASS column should be appended last: {}",
        column_header
    );

    let row_25k = annotated
        .lines()
        .filter(|l| !l.starts_with('#'))
        .find(|r| r.contains("1:25000\t"))
        .expect("expected a row at 1:25000");
    assert!(
        row_25k.ends_with("\tHIGH_QC"),
        "DP=50, QD=30 should classify as HIGH_QC: {}",
        row_25k
    );

    let row_30k = annotated
        .lines()
        .filter(|l| !l.starts_with('#'))
        .find(|r| r.contains("1:30000\t"))
        .expect("expected a row at 1:30000");
    assert!(
        row_30k.ends_with("\tFAIL_QC"),
        "DP=5 falls below LOW_QC threshold (min_dp=8) → fallback: {}",
        row_30k
    );
}
