//! The supplementary sources the classifier reads, and what a run that lacks one cannot say.
//!
//! A criterion with no data reports `evaluated: false` and says nothing else, so a database
//! directory that lacks a source produces no error. It produces a classification that rests on
//! less evidence, and nothing in a VUS tells a reader whether the evidence was weighed and found
//! balanced or never looked at. This module is where the classifier learns which is which.
//!
//! * [`REFERENCE_STACK`] names the nine sources the criteria read, once, with the criteria each
//!   one drives and what a run loses without it. The web server answers "can this genome be
//!   classified?" from it, the command line warns from it, and the documents quote it, so none
//!   of them keeps a list of its own that can drift from the code.
//! * [`LoadedSources`] says which of them a run has open. It exists because a source that is
//!   not loaded and a source that holds no record for a variant look the same in the
//!   annotation: neither adds a field. The one criterion that reads a missing record as
//!   evidence, PM2 ("absent from controls"), must be able to tell them apart, or it fires for
//!   every variant of a run that never opened gnomAD.

use serde::Serialize;

/// One of the nine sources of the reference stack.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum Source {
    /// gnomAD allele counts and frequencies.
    Gnomad,
    /// ClinVar, keyed by allele.
    Clinvar,
    /// REVEL missense scores.
    Revel,
    /// SpliceAI delta scores.
    SpliceAi,
    /// PhyloP conservation.
    PhyloP,
    /// gnomAD gene constraint (pLI, LOEUF, missense and synonymous Z).
    GnomadGenes,
    /// The disease annotation of a gene: ClinGen gene-disease validity, or OMIM.
    Omim,
    /// ClinVar's classified missense variants by protein position, and its canonical splice
    /// variants by genomic position.
    ClinvarProtein,
    /// A repeat-interval track such as RepeatMasker.
    Repeats,
}

/// How a source is keyed.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum SourceScope {
    Allele,
    Gene,
    Interval,
}

/// What the classifier knows about a source.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
pub struct SourceInfo {
    /// Not serialised: `id` names the source for anything outside this crate.
    #[serde(skip)]
    pub source: Source,
    /// The `json_key` its database registers under, which is what `--sa-dir` is read by: the key
    /// comes from the header of the file, and the name of the file does not matter.
    pub id: &'static str,
    pub label: &'static str,
    pub scope: SourceScope,
    /// The criteria that read it.
    pub drives: &'static [&'static str],
    /// What a classification loses without it.
    pub without: &'static str,
}

/// The reference stack, in the order `docs/ACMG_SETUP.md` builds it.
pub const REFERENCE_STACK: [SourceInfo; 9] = [
    SourceInfo {
        source: Source::Gnomad,
        id: "gnomad",
        label: "gnomAD allele frequencies",
        scope: SourceScope::Allele,
        drives: &["BA1", "BS1", "BS2", "PM2"],
        without: "no population frequency: BA1, BS1, BS2 and PM2 cannot be evaluated, so nothing is called too common to be pathogenic",
    },
    SourceInfo {
        source: Source::Clinvar,
        id: "clinvar",
        label: "ClinVar",
        scope: SourceScope::Allele,
        drives: &["PM2", "BA1", "BS1", "BS2", "PM3", "BP2", "PP5", "BP6"],
        without: "PM2 loses the check against ClinVar's own population frequencies, BA1, BS1 and BS2 lose the gate on ClinVar's low-penetrance and risk-allele terms, and PM3 and BP2 cannot read the ClinVar class of a companion variant",
    },
    SourceInfo {
        source: Source::Revel,
        id: "revel",
        label: "REVEL missense scores",
        scope: SourceScope::Allele,
        drives: &["PP3", "BP4"],
        without: "no computational evidence for missense variants: PP3 and BP4 cannot be evaluated for them",
    },
    SourceInfo {
        source: Source::SpliceAi,
        id: "spliceai",
        label: "SpliceAI splice predictions",
        scope: SourceScope::Allele,
        drives: &["PVS1", "PP3", "BP4", "BP7"],
        without: "no splice prediction: PVS1 cannot check a canonical splice site against it, and PP3, BP4 and BP7 cannot be evaluated for splicing",
    },
    SourceInfo {
        source: Source::PhyloP,
        id: "phylop",
        label: "PhyloP conservation",
        scope: SourceScope::Allele,
        drives: &["BP7"],
        without: "BP7 cannot tell a conserved position from one that is not",
    },
    SourceInfo {
        source: Source::GnomadGenes,
        id: "gnomad_genes",
        label: "gnomAD gene constraint",
        scope: SourceScope::Gene,
        drives: &["PVS1", "PP2", "BP1"],
        without: "no gene constraint: PP2 and BP1 cannot be evaluated, and PVS1 loses its constraint evidence",
    },
    SourceInfo {
        source: Source::Omim,
        id: "omim",
        label: "gene-disease validity and inheritance (ClinGen, or OMIM)",
        scope: SourceScope::Gene,
        drives: &["PVS1", "PP2", "PM1", "PM2", "BS2", "PM3", "BP2"],
        without: "no gene-disease validity or inheritance: the validity gate on PVS1, PP2 and PM1 never blocks, and every gene is read as of unknown inheritance by PM2, BS2, PM3 and BP2",
    },
    SourceInfo {
        source: Source::ClinvarProtein,
        id: "clinvar_protein",
        label: "ClinVar protein and splice index",
        scope: SourceScope::Gene,
        drives: &["PS1", "PM1", "PM5", "BP1"],
        without: "no protein or splice index: PS1, PM1, PM5 and BP1 cannot be evaluated",
    },
    SourceInfo {
        source: Source::Repeats,
        id: "repeatmasker",
        label: "RepeatMasker repeat intervals",
        scope: SourceScope::Interval,
        drives: &["BP3"],
        without: "BP3 cannot tell an in-frame change inside a repeat from one outside",
    },
];

impl Source {
    pub const ALL: [Source; 9] = [
        Source::Gnomad,
        Source::Clinvar,
        Source::Revel,
        Source::SpliceAi,
        Source::PhyloP,
        Source::GnomadGenes,
        Source::Omim,
        Source::ClinvarProtein,
        Source::Repeats,
    ];

    /// Which source a database registers as, from its `json_key`, or `None` for a key the
    /// classifier does not read (`dbsnp`, `dbnsfp`, `gerp`, ...).
    ///
    /// The spellings are the ones `extract_classification_input` accepts, because a database
    /// that is loaded under a key the extraction ignores is a database the classifier never
    /// reads; `the_keys_that_name_a_source_are_the_keys_the_extraction_reads` holds the two
    /// together. SpliceAI ships under several capitalisations across builders, and the repeat
    /// track is found by any key that contains `repeat`.
    pub fn from_json_key(key: &str) -> Option<Source> {
        Some(match key {
            "gnomad" => Source::Gnomad,
            "clinvar" => Source::Clinvar,
            "revel" => Source::Revel,
            "spliceai" | "spliceAi" | "spliceAI" | "splice_ai" => Source::SpliceAi,
            "phylop" | "phyloP" => Source::PhyloP,
            "gnomad_genes" | "gnomad_gene" => Source::GnomadGenes,
            "omim" => Source::Omim,
            "clinvar_protein" => Source::ClinvarProtein,
            _ if contains_ignore_ascii_case(key, "repeat") => Source::Repeats,
            _ => return None,
        })
    }

    pub fn info(self) -> &'static SourceInfo {
        &REFERENCE_STACK[self as usize]
    }
}

/// `haystack.to_lowercase().contains(needle)` for a lowercase ASCII needle, without building the
/// lowercase string: this runs for every supplementary entry of every allele.
fn contains_ignore_ascii_case(haystack: &str, needle: &str) -> bool {
    haystack
        .as_bytes()
        .windows(needle.len())
        .any(|w| w.eq_ignore_ascii_case(needle.as_bytes()))
}

/// Which sources of the reference stack a run has open.
///
/// `Default` is "none": the safe reading of a value nobody filled in is that nothing was
/// loaded, so that a criterion that needs a source declines instead of firing on its absence.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct LoadedSources {
    bits: u16,
}

impl LoadedSources {
    /// Every source of the stack, as the benchmark of the paper ran with.
    pub fn all() -> Self {
        Self {
            bits: (1 << Source::ALL.len()) - 1,
        }
    }

    /// The sources a set of loaded databases amounts to, from their `json_key`s. Keys that name
    /// no source of the stack are ignored.
    pub fn from_json_keys<'a>(keys: impl IntoIterator<Item = &'a str>) -> Self {
        let mut loaded = Self::default();
        for key in keys {
            if let Some(source) = Source::from_json_key(key) {
                loaded.bits |= 1 << source as u16;
            }
        }
        loaded
    }

    pub fn has(self, source: Source) -> bool {
        self.bits & (1 << source as u16) != 0
    }

    /// The sources of the stack that are not loaded, in the order of [`REFERENCE_STACK`].
    pub fn missing(self) -> impl Iterator<Item = &'static SourceInfo> {
        REFERENCE_STACK.iter().filter(move |s| !self.has(s.source))
    }

    /// The sources of the stack that are loaded, in the order of [`REFERENCE_STACK`].
    pub fn present(self) -> impl Iterator<Item = &'static SourceInfo> {
        REFERENCE_STACK.iter().filter(move |s| self.has(s.source))
    }

    /// Whether every source of the stack is loaded, which is what the benchmark of the paper
    /// had and what a classification is only comparable with when it has too.
    pub fn is_complete(self) -> bool {
        self.missing().next().is_none()
    }

    /// One sentence naming what is missing and what that costs, for a log line or a message, or
    /// `None` when the stack is complete.
    pub fn shortfall(self) -> Option<String> {
        let missing: Vec<String> = self
            .missing()
            .map(|s| format!("{} ({})", s.label, s.drives.join(", ")))
            .collect();
        if missing.is_empty() {
            return None;
        }
        Some(format!(
            "{} of the {} sources the classifier reads are not loaded: {}",
            missing.len(),
            REFERENCE_STACK.len(),
            missing.join("; ")
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_stack_is_indexed_by_its_enum() {
        // `Source::info` indexes the table by discriminant, so the two must stay in step.
        for (i, info) in REFERENCE_STACK.iter().enumerate() {
            assert_eq!(info.source as usize, i, "{}", info.id);
            assert_eq!(Source::ALL[i], info.source);
            assert_eq!(info.source.info().id, info.id);
        }
    }

    #[test]
    fn nothing_is_loaded_by_default_and_everything_is_in_all() {
        let none = LoadedSources::default();
        assert!(Source::ALL.iter().all(|s| !none.has(*s)));
        assert_eq!(none.missing().count(), 9);
        let all = LoadedSources::all();
        assert!(Source::ALL.iter().all(|s| all.has(*s)));
        assert!(all.is_complete());
        assert_eq!(all.shortfall(), None);
    }

    #[test]
    fn databases_are_recognised_by_the_key_in_their_header() {
        let loaded = LoadedSources::from_json_keys([
            "gnomad",
            "clinvar",
            "spliceAI",
            "phyloP",
            "gnomad_gene",
            "simple_repeat",
            "dbsnp",
            "dbnsfp",
        ]);
        for s in [
            Source::Gnomad,
            Source::Clinvar,
            Source::SpliceAi,
            Source::PhyloP,
            Source::GnomadGenes,
            Source::Repeats,
        ] {
            assert!(loaded.has(s), "{s:?}");
        }
        for s in [Source::Revel, Source::Omim, Source::ClinvarProtein] {
            assert!(!loaded.has(s), "{s:?}");
        }
        // The databases the classifier does not read are not part of the stack.
        assert_eq!(Source::from_json_key("dbsnp"), None);
        assert_eq!(Source::from_json_key("gerp"), None);
        assert_eq!(Source::from_json_key(""), None);
    }

    #[test]
    fn a_repeat_track_is_found_by_any_key_that_contains_repeat_in_any_case() {
        for key in [
            "repeatmasker",
            "RepeatMasker",
            "simple_repeat",
            "my_REPEATS_v2",
        ] {
            assert_eq!(Source::from_json_key(key), Some(Source::Repeats), "{key}");
        }
        assert_eq!(Source::from_json_key("repea"), None);
    }

    #[test]
    fn a_run_without_gnomad_evaluates_none_of_the_four_criteria_that_read_it() {
        use crate::config::AcmgConfig;
        use crate::criteria::evaluate_all_criteria;
        use crate::test_support::minimal_input;
        let mut input = minimal_input();
        input.sources = LoadedSources::from_json_keys(["clinvar"]);
        let criteria = evaluate_all_criteria(&input, &AcmgConfig::default());
        for code in ["BA1", "BS1", "BS2", "PM2"] {
            let c = criteria
                .iter()
                .find(|c| c.code.starts_with(code))
                .unwrap_or_else(|| panic!("{code} is in the list"));
            assert!(!c.evaluated && !c.met, "{code}: {c:?}");
        }
        // With the database open and no record of the variant, PM2 is back: that is absence.
        input.sources = LoadedSources::all();
        let criteria = evaluate_all_criteria(&input, &AcmgConfig::default());
        let pm2 = criteria.iter().find(|c| c.code.starts_with("PM2")).unwrap();
        assert!(pm2.evaluated && pm2.met, "{pm2:?}");
    }

    #[test]
    fn the_shortfall_names_each_missing_source_and_what_it_drives() {
        let loaded = LoadedSources::from_json_keys(["clinvar"]);
        let text = loaded.shortfall().unwrap();
        assert!(text.starts_with("8 of the 9 sources"), "{text}");
        assert!(
            text.contains("gnomAD allele frequencies (BA1, BS1, BS2, PM2)"),
            "{text}"
        );
        assert!(!text.contains("ClinVar ("), "{text}");
        assert_eq!(
            loaded.present().map(|s| s.id).collect::<Vec<_>>(),
            ["clinvar"]
        );
    }
}

/// What each source of [`REFERENCE_STACK`] is claimed to drive is checked against the criteria
/// themselves, not against the documents: a handful of inputs rich enough that every criterion of
/// the stack has something to read are classified with the source and without it, and the
/// criteria whose result differs are compared with the source's `drives`.
///
/// The check that matters is the one-sided one. A criterion that starts reading a source and is
/// not listed under it would make the message "this is what the missing source costs" wrong, so no
/// criterion outside `drives` (and outside the few that the reconciliation pass sets aside beside
/// another, listed under `KNOCK_ON`) may change when a source is blanked. The other side is that
/// blanking a source must change *something* it drives, or the entry is describing a source the
/// inputs never exercise. The same lists were measured on the benchmark, with the nine
/// leave-one-source-out runs of the paper (`analysis/acmg_benchmark/results/v37_manuscript/runs.tsv`),
/// which is where the entries of ClinVar, REVEL, the protein index and the constraint file come from.
#[cfg(test)]
mod what_each_source_drives {
    use super::*;
    use crate::config::AcmgConfig;
    use crate::criteria::evaluate_all_criteria;
    use crate::sa_extract::{
        ClassificationInput, ClinvarData, ClinvarProteinData, ClinvarProteinVariant, GnomadData,
        GnomadGeneData, OmimData, RevelData, SpliceAiData,
    };
    use crate::test_support::minimal_input;
    use fastvep_core::{Consequence, Impact};
    use std::collections::BTreeSet;
    use std::sync::Arc;

    /// Criteria that change only because another criterion that reads the source was set aside or
    /// let through, not because they read it: PP3 is set aside beside PS1 or PM5, which the protein
    /// index feeds. The reconciliation pass has more of these (PM1 beside a Strong PP3, which moves
    /// when REVEL does), and the benchmark shows them at scale (`runs.tsv`); an input here that
    /// brings another one out is for whoever adds it to decide, which is why the list holds only
    /// what these inputs exercise.
    const KNOCK_ON: &[(Source, &[&str])] = &[(Source::ClinvarProtein, &["PP3"])];

    /// `PM2_Supporting` and `PM2` are one criterion at two strengths.
    fn base(code: &str) -> &str {
        for suffix in [
            "_Supporting",
            "_Moderate",
            "_Strong",
            "_Very_Strong",
            "_Very",
        ] {
            if let Some(stem) = code.strip_suffix(suffix) {
                return stem;
            }
        }
        code
    }

    fn state(input: &ClassificationInput) -> BTreeSet<(String, bool, bool)> {
        evaluate_all_criteria(input, &AcmgConfig::default())
            .into_iter()
            .map(|c| (base(&c.code).to_string(), c.evaluated, c.met))
            .collect()
    }

    /// The input as a run without `source` would have built it: the record gone, the fact that
    /// the database is open gone, and whatever the extraction derives from that database gone.
    fn without(input: &ClassificationInput, source: Source) -> ClassificationInput {
        let mut blank = input.clone();
        let mut keys = Vec::new();
        for s in Source::ALL {
            if s != source && input.sources.has(s) {
                keys.push(s.info().id);
            }
        }
        blank.sources = LoadedSources::from_json_keys(keys);
        match source {
            Source::Gnomad => blank.gnomad = None,
            Source::Clinvar => blank.clinvar = None,
            Source::Revel => blank.revel = None,
            Source::SpliceAi => blank.splice_ai = None,
            Source::PhyloP => blank.phylop = None,
            Source::GnomadGenes => blank.gene_constraints = None,
            Source::Omim => blank.omim = None,
            Source::ClinvarProtein => {
                blank.clinvar_protein = None;
                blank.in_critical_region = None;
                blank.same_splice_position_pathogenic = None;
            }
            Source::Repeats => blank.in_repeat_region = None,
        }
        blank
    }

    fn constrained() -> GnomadGeneData {
        GnomadGeneData {
            pli: Some(1.0),
            loeuf: Some(0.1),
            mis_z: Some(4.0),
            syn_z: Some(0.0),
        }
    }

    fn omim(inheritance: &str) -> OmimData {
        OmimData {
            mim_number: Some(1),
            phenotypes: Some(vec![format!(
                "a test disorder (ClinGen Definitive/{inheritance}, MONDO:0000001)"
            )]),
        }
    }

    fn protein_index() -> Arc<ClinvarProteinData> {
        let entry = |alt: &str, n| ClinvarProteinVariant {
            pos: 100,
            ref_aa: "A".into(),
            alt_aa: alt.into(),
            sig: "Pathogenic".into(),
            n,
        };
        let neighbour = |pos| ClinvarProteinVariant {
            pos,
            ref_aa: "G".into(),
            alt_aa: "R".into(),
            sig: "Pathogenic".into(),
            n: 2,
        };
        Arc::new(ClinvarProteinData {
            // PS1 (same change, two independent records), PM5 (another change at the residue),
            // PM1 (a cluster around it) and BP1 (a gene with few pathogenic missense records
            // elsewhere) all read this one index.
            protein_variants: vec![
                entry("D", 3),
                entry("V", 2),
                neighbour(98),
                neighbour(99),
                neighbour(101),
                neighbour(102),
                neighbour(103),
            ],
            benign_indexed: true,
            splice_indexed: true,
            splice_assembly: Some("GRCh38".into()),
            splice_positions: vec![],
        })
    }

    fn splice(score: f64) -> SpliceAiData {
        SpliceAiData {
            ds_ag: Some(score),
            ds_al: Some(score),
            ds_dg: Some(score),
            ds_dl: Some(score),
            ..Default::default()
        }
    }

    /// A missense change in a constrained gene whose disease is dominant, with every record that
    /// points at pathogenic and none at benign.
    fn pathogenic_looking_missense() -> ClassificationInput {
        ClassificationInput {
            consequences: vec![Consequence::MissenseVariant],
            impact: Impact::Moderate,
            amino_acids: Some(("A".into(), "D".into())),
            protein_position: Some(100),
            gnomad: Some(GnomadData {
                all_ac: Some(0),
                all_af: Some(0.0),
                all_an: Some(1_000_000),
                ..Default::default()
            }),
            clinvar: Some(ClinvarData::default()),
            revel: Some(RevelData { score: Some(0.95) }),
            splice_ai: Some(splice(0.01)),
            phylop: Some(4.0),
            gene_constraints: Some(constrained()),
            omim: Some(omim("AD")),
            clinvar_protein: Some(protein_index()),
            ..minimal_input()
        }
    }

    /// A missense change that is common, scores as tolerated, and sits in a gene where the disease
    /// is truncating, with a ClinVar record that calls it low penetrance.
    fn benign_looking_missense() -> ClassificationInput {
        ClassificationInput {
            consequences: vec![Consequence::MissenseVariant],
            impact: Impact::Moderate,
            amino_acids: Some(("A".into(), "T".into())),
            protein_position: Some(500),
            gnomad: Some(GnomadData {
                all_ac: Some(10_000),
                all_af: Some(0.01),
                all_an: Some(1_000_000),
                all_hc: Some(3),
                ..Default::default()
            }),
            clinvar: Some(ClinvarData {
                significance: Some(vec!["Pathogenic, low penetrance".into()]),
                // The gate reads a record only from two stars up.
                review_status: Some("criteria_provided,_multiple_submitters,_no_conflicts".into()),
                af_tgp: Some(0.01),
                ..Default::default()
            }),
            revel: Some(RevelData { score: Some(0.01) }),
            splice_ai: Some(splice(0.01)),
            phylop: Some(0.1),
            gene_constraints: Some(constrained()),
            omim: Some(omim("AR")),
            clinvar_protein: Some(protein_index()),
            ..minimal_input()
        }
    }

    fn synonymous() -> ClassificationInput {
        ClassificationInput {
            consequences: vec![Consequence::SynonymousVariant],
            impact: Impact::Low,
            gnomad: Some(GnomadData {
                all_ac: Some(100),
                all_af: Some(0.0001),
                all_an: Some(1_000_000),
                ..Default::default()
            }),
            clinvar: Some(ClinvarData::default()),
            splice_ai: Some(splice(0.01)),
            phylop: Some(0.1),
            at_exon_edge: Some(false),
            ..minimal_input()
        }
    }

    fn frameshift() -> ClassificationInput {
        ClassificationInput {
            consequences: vec![Consequence::FrameshiftVariant],
            impact: Impact::High,
            gnomad: Some(GnomadData {
                all_ac: Some(0),
                all_af: Some(0.0),
                all_an: Some(1_000_000),
                ..Default::default()
            }),
            clinvar: Some(ClinvarData::default()),
            gene_constraints: Some(constrained()),
            omim: Some(omim("AD")),
            predicted_nmd: Some(true),
            nmd_escape_50nt: Some(false),
            protein_truncation_pct: Some(0.5),
            is_last_exon: Some(false),
            ..minimal_input()
        }
    }

    fn splice_donor() -> ClassificationInput {
        ClassificationInput {
            consequences: vec![Consequence::SpliceDonorVariant],
            impact: Impact::High,
            gnomad: Some(GnomadData {
                all_ac: Some(0),
                all_af: Some(0.0),
                all_an: Some(1_000_000),
                ..Default::default()
            }),
            clinvar: Some(ClinvarData::default()),
            splice_ai: Some(splice(0.9)),
            gene_constraints: Some(constrained()),
            omim: Some(omim("AD")),
            predicted_nmd: Some(true),
            nmd_escape_50nt: Some(false),
            protein_truncation_pct: Some(0.5),
            is_last_exon: Some(false),
            same_splice_position_pathogenic: Some(true),
            clinvar_protein: Some(protein_index()),
            ..minimal_input()
        }
    }

    fn inframe_deletion_in_a_repeat() -> ClassificationInput {
        ClassificationInput {
            consequences: vec![Consequence::InframeDeletion],
            impact: Impact::Moderate,
            gnomad: Some(GnomadData {
                all_ac: Some(100),
                all_af: Some(0.0001),
                all_an: Some(1_000_000),
                ..Default::default()
            }),
            in_repeat_region: Some(true),
            ..minimal_input()
        }
    }

    /// A missense change with a high REVEL score and nothing else to read, so that PP3 is not set
    /// aside beside PS1 or PM5.
    fn missense_with_only_a_revel_score() -> ClassificationInput {
        ClassificationInput {
            consequences: vec![Consequence::MissenseVariant],
            impact: Impact::Moderate,
            amino_acids: Some(("A".into(), "D".into())),
            protein_position: Some(100),
            revel: Some(RevelData { score: Some(0.95) }),
            ..minimal_input()
        }
    }

    /// A frameshift whose only route to loss of function is the constraint file.
    fn frameshift_known_by_constraint_alone() -> ClassificationInput {
        ClassificationInput {
            omim: None,
            ..frameshift()
        }
    }

    /// A frameshift whose only route to loss of function is ClinGen's curation.
    fn frameshift_known_by_clingen_alone() -> ClassificationInput {
        ClassificationInput {
            gene_constraints: None,
            ..frameshift()
        }
    }

    /// A canonical splice site that SpliceAI predicts is still intact.
    fn splice_site_predicted_intact() -> ClassificationInput {
        ClassificationInput {
            splice_ai: Some(splice(0.01)),
            same_splice_position_pathogenic: None,
            clinvar_protein: None,
            ..splice_donor()
        }
    }

    /// A recessive gene, at a frequency between the dominant PM2 bar (4e-5) and the recessive one
    /// (7e-5).
    fn recessive_gene_between_the_pm2_bars() -> ClassificationInput {
        ClassificationInput {
            consequences: vec![Consequence::MissenseVariant],
            impact: Impact::Moderate,
            gnomad: Some(GnomadData {
                all_ac: Some(5),
                all_af: Some(0.00005),
                all_an: Some(100_000),
                ..Default::default()
            }),
            omim: Some(omim("AR")),
            ..minimal_input()
        }
    }

    /// No gnomAD record, and ClinVar's own frequencies say the variant is not rare.
    fn no_gnomad_record_and_a_common_clinvar_frequency() -> ClassificationInput {
        ClassificationInput {
            consequences: vec![Consequence::MissenseVariant],
            impact: Impact::Moderate,
            gnomad: None,
            clinvar: Some(ClinvarData {
                af_tgp: Some(0.01),
                ..Default::default()
            }),
            ..minimal_input()
        }
    }

    /// A missense change in a gene that truncating variants cause disease in, with and without a
    /// protein index that holds enough pathogenic missense records to veto BP1.
    fn missense_in_a_truncating_gene(
        index: Option<Arc<ClinvarProteinData>>,
    ) -> ClassificationInput {
        ClassificationInput {
            consequences: vec![Consequence::MissenseVariant],
            impact: Impact::Moderate,
            amino_acids: Some(("A".into(), "T".into())),
            protein_position: Some(500),
            gene_constraints: Some(GnomadGeneData {
                pli: Some(1.0),
                loeuf: Some(0.1),
                mis_z: Some(0.5),
                syn_z: Some(0.0),
            }),
            clinvar_protein: index,
            ..minimal_input()
        }
    }

    fn scenarios() -> Vec<(&'static str, ClassificationInput)> {
        vec![
            (
                "missense with only a REVEL score",
                missense_with_only_a_revel_score(),
            ),
            (
                "frameshift known by constraint alone",
                frameshift_known_by_constraint_alone(),
            ),
            (
                "frameshift known by ClinGen alone",
                frameshift_known_by_clingen_alone(),
            ),
            (
                "splice site predicted intact",
                splice_site_predicted_intact(),
            ),
            (
                "recessive gene between the PM2 bars",
                recessive_gene_between_the_pm2_bars(),
            ),
            (
                "no gnomAD record and a common ClinVar frequency",
                no_gnomad_record_and_a_common_clinvar_frequency(),
            ),
            (
                "missense in a truncating gene",
                missense_in_a_truncating_gene(None),
            ),
            (
                "missense in a truncating gene with a missense mechanism in the index",
                missense_in_a_truncating_gene(Some(protein_index())),
            ),
            ("pathogenic-looking missense", pathogenic_looking_missense()),
            ("benign-looking missense", benign_looking_missense()),
            ("synonymous", synonymous()),
            ("frameshift", frameshift()),
            ("splice donor", splice_donor()),
            (
                "inframe deletion in a repeat",
                inframe_deletion_in_a_repeat(),
            ),
        ]
    }

    /// The criteria that differ when `source` is blanked, over every scenario that has it.
    fn changed_by(source: Source) -> BTreeSet<String> {
        let mut changed = BTreeSet::new();
        for (name, input) in scenarios() {
            let mut full = input;
            full.sources = LoadedSources::all();
            let before = state(&full);
            let after = state(&without(&full, source));
            for (code, ..) in before.symmetric_difference(&after) {
                eprintln!("{source:?}: {code} changes in the {name}");
                changed.insert(code.clone());
            }
        }
        changed
    }

    #[test]
    fn blanking_a_source_changes_only_criteria_it_is_said_to_drive() {
        for info in REFERENCE_STACK {
            let knock_on: &[&str] = KNOCK_ON
                .iter()
                .find(|(s, _)| *s == info.source)
                .map(|(_, k)| *k)
                .unwrap_or(&[]);
            let changed = changed_by(info.source);
            let undeclared: Vec<&String> = changed
                .iter()
                .filter(|c| !info.drives.contains(&c.as_str()) && !knock_on.contains(&c.as_str()))
                .collect();
            assert!(
                undeclared.is_empty(),
                "{}: criteria {:?} changed when it was removed and are not in its `drives` \
                 ({:?}): add them there, or to KNOCK_ON if another criterion's firing moved them",
                info.id,
                undeclared,
                info.drives
            );
        }
    }

    #[test]
    fn blanking_a_source_changes_something_it_is_said_to_drive() {
        for info in REFERENCE_STACK {
            let changed = changed_by(info.source);
            assert!(
                info.drives.iter().any(|d| changed.contains(*d)),
                "{}: none of {:?} changed when the source was removed (changed: {:?}); either its \
                 entry describes a criterion that never reads it, or the inputs here do not \
                 exercise it",
                info.id,
                info.drives,
                changed
            );
        }
    }
}
