use fastvep_core::{GeneAnnotation, SupplementaryAnnotation};
use fastvep_core::{Consequence, Impact};
use serde::Deserialize;

/// gnomAD population frequency data.
#[derive(Debug, Clone, Default, Deserialize)]
pub struct GnomadData {
    #[serde(rename = "allAf")]
    pub all_af: Option<f64>,
    #[serde(rename = "allAn")]
    pub all_an: Option<u64>,
    #[serde(rename = "allAc")]
    pub all_ac: Option<u64>,
    #[serde(rename = "allHc")]
    pub all_hc: Option<u64>,
    #[serde(rename = "afrAf")]
    pub afr_af: Option<f64>,
    #[serde(rename = "nfeAf")]
    pub nfe_af: Option<f64>,
    #[serde(rename = "easAf")]
    pub eas_af: Option<f64>,
    #[serde(rename = "amrAf")]
    pub amr_af: Option<f64>,
    #[serde(rename = "asjAf")]
    pub asj_af: Option<f64>,
    #[serde(rename = "finAf")]
    pub fin_af: Option<f64>,
    #[serde(rename = "othAf")]
    pub oth_af: Option<f64>,
    #[serde(rename = "sasAf")]
    pub sas_af: Option<f64>,
}

impl GnomadData {
    /// Maximum allele frequency across all populations.
    pub fn max_pop_af(&self) -> Option<f64> {
        [
            self.all_af, self.afr_af, self.nfe_af, self.eas_af, self.amr_af, self.asj_af,
            self.fin_af, self.oth_af, self.sas_af,
        ]
        .into_iter()
        .flatten()
        .reduce(f64::max)
    }
}

/// ClinVar clinical significance data.
#[derive(Debug, Clone, Default, Deserialize)]
pub struct ClinvarData {
    pub significance: Option<Vec<String>>,
    #[serde(rename = "reviewStatus")]
    pub review_status: Option<String>,
    pub phenotypes: Option<Vec<String>>,
    #[serde(rename = "variantClass")]
    pub variant_class: Option<String>,
}

impl ClinvarData {
    /// Check if any significance term contains a pathogenic classification.
    pub fn has_pathogenic(&self) -> bool {
        self.significance.as_ref().map_or(false, |sigs| {
            sigs.iter().any(|s| {
                let lower = s.to_lowercase();
                lower.contains("pathogenic") && !lower.contains("conflicting")
            })
        })
    }

    /// Check if any significance term contains a benign classification.
    pub fn has_benign(&self) -> bool {
        self.significance.as_ref().map_or(false, |sigs| {
            sigs.iter().any(|s| {
                let lower = s.to_lowercase();
                lower.contains("benign") && !lower.contains("conflicting")
            })
        })
    }

    /// Returns the review star level (0-4).
    pub fn review_stars(&self) -> u8 {
        match self.review_status.as_deref() {
            Some(s) if s.contains("practice_guideline") || s.contains("practice guideline") => 4,
            Some(s)
                if s.contains("reviewed_by_expert_panel")
                    || s.contains("reviewed by expert panel") =>
            {
                3
            }
            Some(s) if s.contains("multiple_submitters") || s.contains("multiple submitters") => 2,
            Some(s)
                if (s.contains("criteria_provided") || s.contains("criteria provided"))
                    && !s.contains("no_assertion") && !s.contains("no assertion") =>
            {
                1
            }
            _ => 0,
        }
    }
}

/// REVEL missense pathogenicity score.
#[derive(Debug, Clone, Deserialize)]
pub struct RevelData {
    pub score: Option<f64>,
}

/// SpliceAI delta scores.
#[derive(Debug, Clone, Default, Deserialize)]
pub struct SpliceAiData {
    pub gene: Option<String>,
    #[serde(rename = "dsAg")]
    pub ds_ag: Option<f64>,
    #[serde(rename = "dsAl")]
    pub ds_al: Option<f64>,
    #[serde(rename = "dsDg")]
    pub ds_dg: Option<f64>,
    #[serde(rename = "dsDl")]
    pub ds_dl: Option<f64>,
    #[serde(rename = "dpAg")]
    pub dp_ag: Option<i32>,
    #[serde(rename = "dpAl")]
    pub dp_al: Option<i32>,
    #[serde(rename = "dpDg")]
    pub dp_dg: Option<i32>,
    #[serde(rename = "dpDl")]
    pub dp_dl: Option<i32>,
}

impl SpliceAiData {
    /// Maximum delta score across all four splice effects.
    pub fn max_delta_score(&self) -> Option<f64> {
        [self.ds_ag, self.ds_al, self.ds_dg, self.ds_dl]
            .into_iter()
            .flatten()
            .reduce(f64::max)
    }
}

/// dbNSFP SIFT/PolyPhen predictions.
#[derive(Debug, Clone, Default, Deserialize)]
pub struct DbNsfpData {
    pub sift: Option<String>,
    pub polyphen: Option<String>,
}

/// Parsed prediction result from dbNSFP format strings.
#[derive(Debug, Clone)]
pub struct PredictionResult {
    pub prediction: String,
    pub score: Option<f64>,
}

impl DbNsfpData {
    /// Parse SIFT prediction from format "deleterious(0.001)" or "tolerated(0.123)".
    pub fn parse_sift(&self) -> Option<PredictionResult> {
        self.sift.as_ref().and_then(|s| parse_prediction_string(s))
    }

    /// Parse PolyPhen prediction from format "probably_damaging(0.998)".
    pub fn parse_polyphen(&self) -> Option<PredictionResult> {
        self.polyphen.as_ref().and_then(|s| parse_prediction_string(s))
    }
}

fn parse_prediction_string(s: &str) -> Option<PredictionResult> {
    if let Some(paren_pos) = s.find('(') {
        let prediction = s[..paren_pos].to_string();
        let score = s[paren_pos + 1..]
            .trim_end_matches(')')
            .parse::<f64>()
            .ok();
        Some(PredictionResult { prediction, score })
    } else {
        Some(PredictionResult {
            prediction: s.to_string(),
            score: None,
        })
    }
}

/// gnomAD gene-level constraint data.
#[derive(Debug, Clone, Default, Deserialize)]
pub struct GnomadGeneData {
    #[serde(rename = "pLI")]
    pub pli: Option<f64>,
    pub loeuf: Option<f64>,
    #[serde(rename = "misZ")]
    pub mis_z: Option<f64>,
    #[serde(rename = "synZ")]
    pub syn_z: Option<f64>,
}

/// OMIM gene-disease data.
#[derive(Debug, Clone, Default, Deserialize)]
pub struct OmimData {
    #[serde(rename = "mimNumber")]
    pub mim_number: Option<u32>,
    pub phenotypes: Option<Vec<String>>,
}

impl OmimData {
    /// Check if any phenotype suggests autosomal dominant inheritance.
    pub fn has_dominant_inheritance(&self) -> bool {
        self.phenotypes.as_ref().map_or(false, |ps| {
            ps.iter().any(|p| {
                let lower = p.to_lowercase();
                lower.contains("autosomal dominant")
                    || lower.contains("{ad}")
                    || (lower.contains("dominant") && !lower.contains("recessive"))
            })
        })
    }

    /// Check if any phenotype suggests autosomal recessive inheritance.
    pub fn has_recessive_inheritance(&self) -> bool {
        self.phenotypes.as_ref().map_or(false, |ps| {
            ps.iter().any(|p| {
                let lower = p.to_lowercase();
                lower.contains("autosomal recessive")
                    || lower.contains("{ar}")
                    || (lower.contains("recessive") && !lower.contains("dominant"))
            })
        })
    }
}

/// ClinVar pathogenic variants indexed by protein position (from .oga).
#[derive(Debug, Clone, Default, Deserialize)]
pub struct ClinvarProteinData {
    #[serde(rename = "proteinVariants")]
    pub protein_variants: Vec<ClinvarProteinVariant>,
}

/// A single ClinVar pathogenic variant at a protein position.
#[derive(Debug, Clone, Deserialize)]
pub struct ClinvarProteinVariant {
    pub pos: u64,
    #[serde(rename = "refAa")]
    pub ref_aa: String,
    #[serde(rename = "altAa")]
    pub alt_aa: String,
    pub sig: String,
}

/// Genotype information for a sample at a specific variant.
#[derive(Debug, Clone)]
pub struct GenotypeInfo {
    pub is_het: bool,
    pub is_hom_ref: bool,
    pub is_hom_alt: bool,
    pub is_missing: bool,
    pub is_phased: bool,
    pub depth: Option<u32>,
    pub quality: Option<u32>,
    /// Which alt allele index this genotype carries (1-based). None if hom_ref or missing.
    pub alt_allele_index: Option<u32>,
}

impl GenotypeInfo {
    /// Returns true if genotype passes depth and quality thresholds.
    pub fn passes_quality(&self, min_depth: u32, min_gq: u32) -> bool {
        let depth_ok = self.depth.map_or(false, |d| d >= min_depth);
        let gq_ok = self.quality.map_or(false, |q| q >= min_gq);
        depth_ok && gq_ok
    }

    /// Returns true if the sample carries the variant allele (het or hom_alt).
    pub fn carries_variant(&self) -> bool {
        self.is_het || self.is_hom_alt
    }
}

/// Information about another variant in the same gene for compound-het analysis.
#[derive(Debug, Clone)]
pub struct CompanionVariant {
    /// Whether the companion variant is ClinVar pathogenic
    pub is_clinvar_pathogenic: bool,
    /// Whether the companion variant is ClinVar likely pathogenic.
    /// Used by PM3 v1.0 points scoring (PR7) — LP companions earn fewer
    /// points than P companions. Defaults to false to preserve back-compat.
    pub is_clinvar_likely_pathogenic: bool,
    /// Whether variants are in trans (different haplotypes). None = unphased.
    pub is_in_trans: Option<bool>,
    /// Whether proband is heterozygous for the companion variant
    pub proband_het: bool,
    /// HGVSc of the companion variant for reporting
    pub hgvsc: Option<String>,
}

/// All data needed for ACMG-AMP classification, extracted from pipeline annotations.
#[derive(Debug, Clone)]
pub struct ClassificationInput {
    pub consequences: Vec<Consequence>,
    pub impact: Impact,
    pub gene_symbol: Option<String>,
    pub is_canonical: bool,
    pub amino_acids: Option<(String, String)>,
    pub protein_position: Option<u64>,
    pub gnomad: Option<GnomadData>,
    pub clinvar: Option<ClinvarData>,
    pub revel: Option<RevelData>,
    pub splice_ai: Option<SpliceAiData>,
    pub dbnsfp: Option<DbNsfpData>,
    pub phylop: Option<f64>,
    pub gerp: Option<f64>,
    pub gene_constraints: Option<GnomadGeneData>,
    pub omim: Option<OmimData>,
    /// ClinVar pathogenic variants at protein positions for this gene (from .oga).
    pub clinvar_protein: Option<ClinvarProteinData>,
    /// HGVS c. notation for the variant (e.g. "c.845G>A"). Used by the BA1
    /// exception list lookup (Ghosh 2018) to identify well-known high-AF
    /// pathogenic variants that must not call BA1. `None` if the pipeline
    /// did not produce HGVS — BA1 then proceeds with its default behavior.
    pub hgvs_c: Option<String>,
    // ── PVS1 decision-tree signals (Abou Tayoun 2018) ────────────────────
    // All Optional. The PVS1 evaluator uses them to grade the strength of
    // null-variant evidence (PVS1 / PVS1_Strong / PVS1_Moderate /
    // PVS1_Supporting / N/A). When unpopulated, PVS1 falls back to its
    // legacy binary behavior (full Very Strong if a null variant is in a
    // LOF-intolerant gene). The pipeline (fastvep-cli) computes these from
    // cached transcript exon coordinates + ClinVar protein index.
    /// True if the predicted premature termination is expected to undergo
    /// nonsense-mediated decay (NOT in 3'-most exon AND NOT in last 50 nt
    /// of penultimate exon).
    pub predicted_nmd: Option<bool>,
    /// Fraction of the protein removed by the variant (0.0–1.0).
    pub protein_truncation_pct: Option<f64>,
    /// True if the variant lies in the 3'-most (last) exon.
    pub is_last_exon: Option<bool>,
    /// True if downstream pathogenic variants exist past the truncation point
    /// (used as a proxy for "critical functional region").
    pub in_critical_region: Option<bool>,
    /// Distance (in codons) to the next downstream Met for start-loss
    /// variants. None if no alternative start codon exists.
    pub alt_start_codon_distance: Option<i64>,
    /// PS1 splice path (Walker 2023): set to true when this canonical ±1/2
    /// splice variant matches a known pathogenic splice variant predicted to
    /// produce the same RNA outcome (e.g. same intron, same direction of
    /// splice loss). The pipeline computes this from a position-indexed
    /// ClinVar splice catalog. None when the data isn't available.
    pub same_splice_position_pathogenic: Option<bool>,
    /// Whether variant overlaps a repeat region (from RepeatMasker .osi).
    pub in_repeat_region: Option<bool>,
    /// Whether the variant sits at the first base or last 3 bases of an exon
    /// (the canonical splice region). Per Walker 2023 (ClinGen SVI Splicing
    /// Subgroup), BP7 must NOT fire for synonymous variants at these positions
    /// because SpliceAI can miss splice impact in this region. `None` means
    /// the pipeline didn't populate this signal (BP7 falls back to legacy
    /// behavior — fire if SpliceAI low + PhyloP low).
    pub at_exon_edge: Option<bool>,
    /// Signed offset in bp from the nearest exon boundary, for intronic
    /// variants. Convention: positive after the donor (e.g. c.123+5 → +5),
    /// negative before the acceptor (e.g. c.123-15 → -15). Per Walker 2023,
    /// BP7 may extend to intronic variants outside the standard splice region:
    /// donor-side `offset ≥ 7` or acceptor-side `offset ≤ -21`. `None` for
    /// non-intronic variants or when the pipeline can't compute it.
    pub intronic_offset: Option<i64>,
    /// Proband genotype information (from trio VCF)
    pub proband_genotype: Option<GenotypeInfo>,
    /// Mother genotype information (from trio VCF)
    pub mother_genotype: Option<GenotypeInfo>,
    /// Father genotype information (from trio VCF)
    pub father_genotype: Option<GenotypeInfo>,
    /// Other variants in the same gene that the proband carries (for compound-het)
    pub companion_variants: Vec<CompanionVariant>,
}

/// Extract classification input from pipeline annotation data.
///
/// Parses the pre-serialized JSON strings from supplementary annotations into typed structs.
pub fn extract_classification_input(
    consequences: &[Consequence],
    impact: Impact,
    gene_symbol: Option<&str>,
    is_canonical: bool,
    amino_acids: Option<&(String, String)>,
    protein_position: Option<u64>,
    allele_supplementary: &[(String, String)],
    gene_annotations: &[&GeneAnnotation],
    variant_supplementary: &[SupplementaryAnnotation],
    proband_genotype: Option<GenotypeInfo>,
    mother_genotype: Option<GenotypeInfo>,
    father_genotype: Option<GenotypeInfo>,
    companion_variants: Vec<CompanionVariant>,
) -> ClassificationInput {
    let mut gnomad = None;
    let mut clinvar = None;
    let mut revel = None;
    let mut splice_ai = None;
    let mut dbnsfp = None;

    // Parse per-allele supplementary annotations
    for (key, json_str) in allele_supplementary {
        match key.as_str() {
            "gnomad" => {
                gnomad = serde_json::from_str(json_str).ok();
            }
            "clinvar" => {
                clinvar = serde_json::from_str(json_str).ok();
            }
            "revel" => {
                revel = serde_json::from_str(json_str).ok();
            }
            "spliceai" | "spliceAi" | "splice_ai" => {
                splice_ai = serde_json::from_str(json_str).ok();
            }
            "dbnsfp" => {
                dbnsfp = serde_json::from_str(json_str).ok();
            }
            _ => {}
        }
    }

    // Parse positional supplementary annotations (PhyloP, GERP)
    let mut phylop = None;
    let mut gerp = None;
    for sa in variant_supplementary {
        match sa.json_key.as_str() {
            "phylop" | "phyloP" => {
                phylop = sa.json_string.trim_matches('"').parse::<f64>().ok();
            }
            "gerp" | "GERP" => {
                gerp = sa.json_string.trim_matches('"').parse::<f64>().ok();
            }
            _ => {}
        }
    }

    // Parse gene-level annotations
    let mut gene_constraints = None;
    let mut omim = None;
    let mut clinvar_protein = None;
    for ga in gene_annotations {
        match ga.json_key.as_str() {
            "gnomad_genes" | "gnomad_gene" => {
                gene_constraints = serde_json::from_str(&ga.json_string).ok();
            }
            "omim" => {
                omim = serde_json::from_str(&ga.json_string).ok();
            }
            "clinvar_protein" => {
                clinvar_protein = serde_json::from_str(&ga.json_string).ok();
            }
            _ => {}
        }
    }

    // Check if variant overlaps a repeat region (from interval SA)
    let in_repeat_region = {
        let has_repeat = allele_supplementary.iter().any(|(key, _)| {
            let k = key.to_lowercase();
            k.contains("repeat") || k.contains("repeatmasker") || k.contains("simple_repeat")
        });
        if has_repeat { Some(true) } else { None }
    };

    ClassificationInput {
        consequences: consequences.to_vec(),
        impact,
        gene_symbol: gene_symbol.map(|s| s.to_string()),
        is_canonical,
        amino_acids: amino_acids.cloned(),
        protein_position,
        gnomad,
        clinvar,
        revel,
        splice_ai,
        dbnsfp,
        phylop,
        gerp,
        gene_constraints,
        omim,
        clinvar_protein,
        // hgvs_c is populated by the pipeline once HGVS is computed; remains
        // None for now and BA1 falls back to default behavior.
        hgvs_c: None,
        // PVS1 decision-tree signals — populated once the pipeline plumbing
        // (transcript exon coords + ClinVar protein index) lands. Until
        // then, PVS1 falls back to its legacy binary rule.
        predicted_nmd: None,
        protein_truncation_pct: None,
        is_last_exon: None,
        in_critical_region: None,
        alt_start_codon_distance: None,
        same_splice_position_pathogenic: None,
        in_repeat_region,
        // BP7 exon-edge / deep-intronic signals (Walker 2023). The pipeline
        // populates these once per-transcript exon coordinates are wired in;
        // until then they remain None and BP7 falls back to legacy behavior.
        at_exon_edge: None,
        intronic_offset: None,
        proband_genotype,
        mother_genotype,
        father_genotype,
        companion_variants,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gnomad_max_pop_af() {
        let g = GnomadData {
            all_af: Some(0.001),
            afr_af: Some(0.05),
            nfe_af: Some(0.0001),
            ..Default::default()
        };
        assert!((g.max_pop_af().unwrap() - 0.05).abs() < 1e-10);
    }

    #[test]
    fn test_gnomad_max_pop_af_none() {
        let g = GnomadData::default();
        assert!(g.max_pop_af().is_none());
    }

    #[test]
    fn test_clinvar_pathogenic() {
        let c = ClinvarData {
            significance: Some(vec!["Pathogenic".to_string()]),
            ..Default::default()
        };
        assert!(c.has_pathogenic());
        assert!(!c.has_benign());
    }

    #[test]
    fn test_clinvar_likely_pathogenic() {
        let c = ClinvarData {
            significance: Some(vec!["Likely_pathogenic".to_string()]),
            ..Default::default()
        };
        assert!(c.has_pathogenic());
    }

    #[test]
    fn test_clinvar_conflicting_not_pathogenic() {
        let c = ClinvarData {
            significance: Some(vec![
                "Conflicting_interpretations_of_pathogenicity".to_string(),
            ]),
            ..Default::default()
        };
        assert!(!c.has_pathogenic());
    }

    #[test]
    fn test_clinvar_review_stars() {
        assert_eq!(
            ClinvarData {
                review_status: Some(
                    "criteria_provided,_multiple_submitters,_no_conflicts".to_string()
                ),
                ..Default::default()
            }
            .review_stars(),
            2
        );
        assert_eq!(
            ClinvarData {
                review_status: Some("reviewed_by_expert_panel".to_string()),
                ..Default::default()
            }
            .review_stars(),
            3
        );
        assert_eq!(
            ClinvarData {
                review_status: Some("criteria_provided,_single_submitter".to_string()),
                ..Default::default()
            }
            .review_stars(),
            1
        );
        assert_eq!(
            ClinvarData {
                review_status: Some("no_assertion_criteria_provided".to_string()),
                ..Default::default()
            }
            .review_stars(),
            0
        );
    }

    #[test]
    fn test_spliceai_max_delta() {
        let s = SpliceAiData {
            ds_ag: Some(0.01),
            ds_al: Some(0.85),
            ds_dg: Some(0.02),
            ds_dl: Some(0.10),
            ..Default::default()
        };
        assert!((s.max_delta_score().unwrap() - 0.85).abs() < 1e-10);
    }

    #[test]
    fn test_parse_prediction_string() {
        let r = parse_prediction_string("deleterious(0.001)").unwrap();
        assert_eq!(r.prediction, "deleterious");
        assert!((r.score.unwrap() - 0.001).abs() < 1e-10);

        let r = parse_prediction_string("probably_damaging(0.998)").unwrap();
        assert_eq!(r.prediction, "probably_damaging");
        assert!((r.score.unwrap() - 0.998).abs() < 1e-10);

        let r = parse_prediction_string("tolerated").unwrap();
        assert_eq!(r.prediction, "tolerated");
        assert!(r.score.is_none());
    }

    #[test]
    fn test_gnomad_deserialization() {
        let json = r#"{"allAf":1.234e-03,"allAn":150000,"allAc":150,"allHc":2,"afrAf":2.0e-03,"nfeAf":5.0e-04}"#;
        let g: GnomadData = serde_json::from_str(json).unwrap();
        assert!((g.all_af.unwrap() - 0.001234).abs() < 1e-8);
        assert_eq!(g.all_an.unwrap(), 150000);
        assert!((g.afr_af.unwrap() - 0.002).abs() < 1e-8);
    }

    #[test]
    fn test_clinvar_deserialization() {
        let json = r#"{"significance":["Pathogenic"],"reviewStatus":"criteria_provided,_multiple_submitters,_no_conflicts","phenotypes":["Breast_cancer"],"variantClass":"SNV"}"#;
        let c: ClinvarData = serde_json::from_str(json).unwrap();
        assert!(c.has_pathogenic());
        assert_eq!(c.review_stars(), 2);
    }
}
