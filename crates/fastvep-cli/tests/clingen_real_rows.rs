//! Real ClinGen curations, through the converter and the `.oga` builder, read by the classifier.
//!
//! `fixtures/clingen/genemap.tsv` is what `clingen_gdv_to_oga.py` made of 263 real curations in 185
//! genes of ClinGen's export of 11 October 2026 (`make_clingen_parse_fixture.py` says how they were
//! picked and why), and `fixtures/clingen/expected.tsv` is what each gene must read as, worked out from
//! the curations' own columns and not from the converter's output. The fixture holds every mode of
//! inheritance the file uses, genes with an autosomal dominant and an autosomal recessive disease, genes
//! whose curations are all weak, the three labels that contain a `;` (which the builder splits into
//! pieces, one of them without the tag) and the three X-linked disorders whose label says "dominant".
//!
//! The unit tests of the classifier hand it phenotype strings written for the purpose. These are the
//! strings the pipeline really produces, so a change to the converter's format, to the builder's
//! splitting, or to the classifier's reading of the tag shows up here as a gene that stopped reading as
//! its curation says.

use fastvep_classification::sa_extract::OmimData;
use fastvep_cli::pipeline::run_sa_build;
use fastvep_sa::gene::GeneIndex;
use std::collections::BTreeMap;
use std::fs;

const GENEMAP: &str = include_str!("../fixtures/clingen/genemap.tsv");
const EXPECTED: &str = include_str!("../fixtures/clingen/expected.tsv");

struct Expect {
    established: bool,
    no_valid: bool,
    dominant: bool,
    recessive: bool,
    contexts: String,
}

fn expectations() -> BTreeMap<String, Expect> {
    EXPECTED
        .lines()
        .skip(1)
        .map(|line| {
            let f: Vec<&str> = line.split('\t').collect();
            (
                f[0].to_string(),
                Expect {
                    established: f[1] == "true",
                    no_valid: f[2] == "true",
                    dominant: f[3] == "true",
                    recessive: f[4] == "true",
                    contexts: f.get(5).copied().unwrap_or("").to_string(),
                },
            )
        })
        .collect()
}

/// The `.oga` the CLI builds from the fixture, read back gene by gene.
fn build() -> (tempfile::TempDir, GeneIndex) {
    let tmp = tempfile::tempdir().unwrap();
    let input = tmp.path().join("genemap.txt");
    let output = tmp.path().join("omim");
    fs::write(&input, GENEMAP).unwrap();
    run_sa_build(
        "omim",
        input.to_str().unwrap(),
        output.to_str().unwrap(),
        "GRCh38",
        None,
        &[],
        false,
    )
    .unwrap();
    let index = GeneIndex::read_from(&mut fs::File::open(output.with_extension("oga")).unwrap())
        .expect("the .oga reads back");
    (tmp, index)
}

fn annotation(index: &GeneIndex, gene: &str) -> OmimData {
    let records = index
        .get(gene)
        .unwrap_or_else(|| panic!("{gene} is not in the .oga"));
    serde_json::from_str(records.first().expect("a record")).unwrap()
}

fn contexts_of(omim: &OmimData) -> String {
    omim.disorder_contexts()
        .iter()
        .map(|c| format!("{}:{}", c.inheritance_code(), c.entries.len()))
        .collect::<Vec<_>>()
        .join("|")
}

#[test]
fn every_gene_reads_as_its_curations_say() {
    let (_tmp, index) = build();
    let expected = expectations();
    assert_eq!(
        expected.len(),
        185,
        "the fixture changed: regenerate its expectations"
    );
    assert_eq!(index.gene_count(), expected.len());

    let mut wrong = Vec::new();
    for (gene, want) in &expected {
        let omim = annotation(&index, gene);
        let got = (
            omim.has_established_relationship(),
            omim.has_no_valid_relationship(),
            omim.has_dominant_inheritance(),
            omim.has_recessive_inheritance(),
            contexts_of(&omim),
        );
        let want = (
            want.established,
            want.no_valid,
            want.dominant,
            want.recessive,
            want.contexts.clone(),
        );
        if got != want {
            wrong.push(format!("{gene}: read {got:?}, curations say {want:?}"));
        }
    }
    assert!(
        wrong.is_empty(),
        "{} genes:\n{}",
        wrong.len(),
        wrong.join("\n")
    );
}

/// The three labels of the file that contain a `;`. The builder splits phenotypes on it, so the gene
/// arrives as a piece of the label with no tag beside the piece that has one.
#[test]
fn a_label_with_a_semicolon_is_one_disease_and_not_two() {
    let (_tmp, index) = build();
    for gene in ["MAPK8IP3", "SCAPER", "NSMCE3"] {
        let omim = annotation(&index, gene);
        let phenotypes = omim.phenotypes.as_ref().unwrap();
        assert!(
            phenotypes.len() >= 2,
            "{gene}: the builder splits the label, which is what this guards: {phenotypes:?}"
        );
        assert_eq!(
            omim.entries().count(),
            1,
            "{gene}: the untagged piece is not a disease: {phenotypes:?}"
        );
    }
    // NSMCE3's only curation is Limited, and it used to escape the validity gate by the untagged piece.
    assert!(annotation(&index, "NSMCE3").has_no_valid_relationship());
    // The other two are Definitive and stay established, with the inheritance of the curation.
    assert!(annotation(&index, "MAPK8IP3").has_dominant_inheritance());
    assert!(annotation(&index, "SCAPER").has_recessive_inheritance());
}

/// Three X-linked disorders have "dominant" in the name. The label is not the curated field.
#[test]
fn an_x_linked_disorder_called_dominant_is_x_linked() {
    let (_tmp, index) = build();
    let mut checked = 0;
    for (gene, want) in expectations() {
        let omim = annotation(&index, &gene);
        let named_dominant = omim.entries().any(|e| {
            e.label.to_lowercase().contains("dominant")
                && e.clingen.is_some_and(|t| {
                    t.moi == Some(fastvep_classification::ModeOfInheritance::XLinked)
                })
        });
        if named_dominant && want.contexts == "XL:1" {
            assert!(!omim.has_dominant_inheritance(), "{gene}");
            checked += 1;
        }
    }
    assert!(checked >= 1, "the fixture holds no such gene");
}

/// The review's examples: a recessive gene whose disease name does not say so, one that is dominant and
/// recessive, and the semidominant ABCC8 with both besides.
#[test]
fn the_genes_the_review_named_read_as_curated() {
    let (_tmp, index) = build();
    let pah = annotation(&index, "PAH");
    assert!(pah.has_recessive_inheritance() && !pah.has_dominant_inheritance());
    assert!(!pah.has_several_disorder_contexts());

    let brca1 = annotation(&index, "BRCA1");
    assert!(brca1.has_dominant_inheritance() && brca1.has_recessive_inheritance());
    assert!(brca1.has_several_disorder_contexts());
    assert_eq!(contexts_of(&brca1), "AD:1|AR:1");

    let abcc8 = annotation(&index, "ABCC8");
    assert_eq!(contexts_of(&abcc8), "AD:2|AR:1|SD:1");
}
