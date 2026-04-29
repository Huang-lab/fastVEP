use serde::{Deserialize, Serialize};
use std::collections::HashMap;

use crate::types::EvidenceStrength;

/// Trio configuration for de novo and compound heterozygote analysis.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TrioConfig {
    /// Sample name of the proband (required)
    pub proband: String,
    /// Sample name of the mother (optional)
    pub mother: Option<String>,
    /// Sample name of the father (optional)
    pub father: Option<String>,
    /// Minimum read depth for reliable genotype call (default: 10)
    #[serde(default = "default_min_depth")]
    pub min_depth: u32,
    /// Minimum genotype quality for reliable genotype call (default: 20)
    #[serde(default = "default_min_gq")]
    pub min_gq: u32,
}

fn default_min_depth() -> u32 {
    10
}
fn default_min_gq() -> u32 {
    20
}

/// Configuration for ACMG-AMP classification thresholds and behavior.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AcmgConfig {
    // ── Frequency thresholds ──
    /// BA1: allele frequency threshold for standalone benign (default: 0.05)
    #[serde(default = "default_ba1")]
    pub ba1_af_threshold: f64,
    /// BS1: allele frequency threshold for strong benign (default: 0.01)
    #[serde(default = "default_bs1")]
    pub bs1_af_threshold: f64,
    /// PM2: legacy single allele frequency threshold (default: 0.0001).
    ///
    /// Retained for backward compatibility with configs predating PR4. The
    /// classifier prefers the inheritance-aware fields below
    /// (`pm2_ad_af_threshold` / `pm2_ar_af_threshold`); this field is no
    /// longer consulted by the default code path but remains so that existing
    /// TOML configs continue to deserialize.
    #[serde(default = "default_pm2")]
    pub pm2_af_threshold: f64,
    /// PM2: AF threshold for autosomal-dominant or unknown-inheritance genes
    /// per ClinGen SVI v1.0 (Sept 2020). Default 0.0 → strict absence (AC = 0).
    #[serde(default = "default_pm2_ad")]
    pub pm2_ad_af_threshold: f64,
    /// PM2: AF threshold for autosomal-recessive genes per ClinGen SVI v1.0
    /// (Sept 2020). Default 0.00007 (0.007%).
    #[serde(default = "default_pm2_ar")]
    pub pm2_ar_af_threshold: f64,

    // ── REVEL thresholds (ClinGen SVI calibrated) ──
    /// PP3 supporting threshold (default: 0.644)
    #[serde(default = "default_pp3_supporting")]
    pub pp3_revel_supporting: f64,
    /// PP3 moderate threshold (default: 0.773)
    #[serde(default = "default_pp3_moderate")]
    pub pp3_revel_moderate: f64,
    /// PP3 strong threshold (default: 0.932)
    #[serde(default = "default_pp3_strong")]
    pub pp3_revel_strong: f64,
    /// BP4 supporting threshold (default: 0.290)
    #[serde(default = "default_bp4_supporting")]
    pub bp4_revel_supporting: f64,
    /// BP4 moderate threshold (default: 0.183)
    #[serde(default = "default_bp4_moderate")]
    pub bp4_revel_moderate: f64,
    /// BP4 strong threshold (default: 0.016)
    #[serde(default = "default_bp4_strong")]
    pub bp4_revel_strong: f64,
    /// BP4 very strong threshold (default: 0.003) — Pejaver 2022 endorses
    /// REVEL ≤ 0.003 as Very Strong benign evidence; only REVEL reaches this band.
    #[serde(default = "default_bp4_very_strong")]
    pub bp4_revel_very_strong: f64,

    // ── SpliceAI thresholds ──
    /// SpliceAI delta score threshold for PP3 pathogenic splice impact (default: 0.2).
    /// Per Walker 2023 SVI Splicing Subgroup, ≥ 0.2 yields PP3 at *Supporting* strength only.
    /// SpliceAI alone does not reach Strong — that requires experimental RNA assay (PVS1_RNA).
    #[serde(default = "default_spliceai_pathogenic")]
    pub spliceai_pathogenic: f64,
    /// SpliceAI delta score upper bound for BP4 benign splice impact (default: 0.1).
    /// Per Walker 2023 SVI Splicing Subgroup, ≤ 0.1 yields BP4 at Supporting strength.
    /// Scores between 0.1 and 0.2 are uninformative.
    #[serde(default = "default_spliceai_benign")]
    pub spliceai_benign: f64,

    // ── Conservation thresholds ──
    /// PhyloP threshold for conserved position (default: 2.0)
    #[serde(default = "default_phylop")]
    pub phylop_conserved: f64,

    // ── Gene constraint thresholds ──
    /// pLI threshold for LOF-intolerant gene (default: 0.9)
    #[serde(default = "default_pli")]
    pub pli_lof_intolerant: f64,
    /// LOEUF threshold for LOF-intolerant gene (default: 0.35)
    #[serde(default = "default_loeuf")]
    pub loeuf_lof_intolerant: f64,
    /// Missense Z-score threshold for PP2 (default: 3.09)
    #[serde(default = "default_misz")]
    pub pp2_misz_threshold: f64,

    // ── PM1 hotspot detection thresholds ──
    /// Window size (in amino acid positions) for hotspot detection (default: 5)
    #[serde(default = "default_pm1_window")]
    pub pm1_hotspot_window: u64,
    /// Minimum pathogenic variants in window to call hotspot (default: 3)
    #[serde(default = "default_pm1_threshold")]
    pub pm1_hotspot_min_pathogenic: u32,

    // ── ClinGen SVI behavior modifications ──
    /// Downgrade PM2 from Moderate to Supporting (ClinGen SVI recommendation)
    #[serde(default = "default_true")]
    pub pm2_downgrade_to_supporting: bool,
    /// Enable PP5/BP6 criteria (disabled by default per ClinGen SVI)
    #[serde(default)]
    pub use_pp5_bp6: bool,
    /// Opt back into the legacy PS4 proxy: ClinVar pathogenic with ≥3 review
    /// stars → PS4. ClinGen SVI considers this proxy invalid (true PS4
    /// requires case-control statistics), so it is disabled by default. Set
    /// `true` only for backward-comparable benchmarks.
    #[serde(default)]
    pub use_clinvar_stars_as_ps4_proxy: bool,

    /// Variants exempted from BA1 per the ClinGen SVI updated recommendation
    /// (Ghosh et al. 2018, Hum Mutat). These are well-known high-AF variants
    /// whose pathogenicity is established despite exceeding the 5% AF threshold.
    /// Defaults to the original 9-variant list; users may add VCEP-specific
    /// entries via TOML.
    #[serde(default = "default_ba1_exceptions")]
    pub ba1_exceptions: Vec<Ba1Exception>,

    /// Minimum allele number (AN) required for BA1 / BS1 to fire (gnomAD v4
    /// guidance, ClinGen SVI March 2024). With v4's massive expansion (807k
    /// exomes, 76k genomes), the overall AN should be ≥ 2000 before a
    /// frequency-based call is reliable. When the AN drops below this
    /// threshold, BA1/BS1 are marked NotEvaluated rather than firing on
    /// noisy frequency estimates. Default 2000.
    #[serde(default = "default_min_an")]
    pub min_an_for_frequency_criteria: u64,

    // ── Gene-specific overrides ──
    #[serde(default)]
    pub gene_overrides: HashMap<String, GeneOverride>,

    // ── Trio configuration ──
    /// Trio configuration for de novo and compound heterozygote analysis.
    #[serde(default)]
    pub trio: Option<TrioConfig>,
}

/// One entry on the BA1 exception list. A variant is exempt from BA1 when its
/// `(gene_symbol, hgvs_c)` matches an entry, regardless of allele frequency.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Ba1Exception {
    pub gene: String,
    /// HGVS c. notation, e.g. "c.845G>A". Compared case-insensitively.
    pub hgvs_c: String,
    /// Optional human-readable reason — surfaced in the criterion `summary`.
    #[serde(default)]
    pub reason: Option<String>,
}

/// Gene-specific overrides for ACMG-AMP criteria.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GeneOverride {
    /// Disease mechanism: "LOF", "GOF", "LOF_and_GOF"
    pub mechanism: Option<String>,
    /// Override BS1 allele frequency threshold
    pub bs1_af_threshold: Option<f64>,
    /// Override PM2 allele frequency threshold
    pub pm2_af_threshold: Option<f64>,
    /// Criteria codes to disable for this gene
    #[serde(default)]
    pub disabled_criteria: Vec<String>,
    /// Criteria strength overrides (code -> new strength)
    #[serde(default)]
    pub strength_overrides: HashMap<String, EvidenceStrength>,
    /// Per-disorder thresholds for genes associated with multiple disorders
    /// (ClinGen SVI guidance July 2025). The classifier consumes whichever
    /// disorder context is active for the call; in the absence of explicit
    /// disorder selection, this scaffold is currently informational only —
    /// the active disorder selection mechanism is part of a follow-up PR.
    #[serde(default)]
    pub disorders: HashMap<String, DisorderOverride>,
}

/// Per-disorder override values within a multi-disorder gene.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct DisorderOverride {
    /// Inheritance for this disorder ("AD", "AR", or "AD_AR").
    pub inheritance: Option<String>,
    /// Override BS1 AF threshold for this disorder.
    pub bs1_af_threshold: Option<f64>,
    /// Override PM2 AF threshold for this disorder.
    pub pm2_af_threshold: Option<f64>,
}

impl Default for AcmgConfig {
    fn default() -> Self {
        Self {
            ba1_af_threshold: 0.05,
            bs1_af_threshold: 0.01,
            pm2_af_threshold: 0.0001,
            pm2_ad_af_threshold: 0.0,
            pm2_ar_af_threshold: 0.00007,
            pp3_revel_supporting: 0.644,
            pp3_revel_moderate: 0.773,
            pp3_revel_strong: 0.932,
            bp4_revel_supporting: 0.290,
            bp4_revel_moderate: 0.183,
            bp4_revel_strong: 0.016,
            bp4_revel_very_strong: 0.003,
            spliceai_pathogenic: 0.2,
            spliceai_benign: 0.1,
            phylop_conserved: 2.0,
            pli_lof_intolerant: 0.9,
            loeuf_lof_intolerant: 0.35,
            pp2_misz_threshold: 3.09,
            pm1_hotspot_window: 5,
            pm1_hotspot_min_pathogenic: 3,
            pm2_downgrade_to_supporting: true,
            use_pp5_bp6: false,
            ba1_exceptions: default_ba1_exceptions(),
            use_clinvar_stars_as_ps4_proxy: false,
            min_an_for_frequency_criteria: 2000,
            gene_overrides: HashMap::new(),
            trio: None,
        }
    }
}

impl AcmgConfig {
    /// Load configuration from a TOML file.
    pub fn from_toml_file(path: &str) -> anyhow::Result<Self> {
        let content = std::fs::read_to_string(path)?;
        let config: Self = toml::from_str(&content)?;
        Ok(config)
    }

    /// Get gene-specific override, if any.
    pub fn gene_override(&self, gene: &str) -> Option<&GeneOverride> {
        self.gene_overrides.get(gene)
    }

    /// Check if a criterion is disabled for a given gene.
    pub fn is_criterion_disabled(&self, gene: &str, criterion_code: &str) -> bool {
        self.gene_overrides
            .get(gene)
            .map_or(false, |o| o.disabled_criteria.iter().any(|c| c == criterion_code))
    }

    /// Get effective BS1 threshold for a gene (gene-specific or default).
    pub fn effective_bs1_threshold(&self, gene: Option<&str>) -> f64 {
        gene.and_then(|g| {
            self.gene_overrides
                .get(g)
                .and_then(|o| o.bs1_af_threshold)
        })
        .unwrap_or(self.bs1_af_threshold)
    }

}

// Default value functions for serde
fn default_ba1() -> f64 { 0.05 }
fn default_bs1() -> f64 { 0.01 }
fn default_pm2() -> f64 { 0.0001 }
fn default_pm2_ad() -> f64 { 0.0 }
fn default_pm2_ar() -> f64 { 0.00007 }
fn default_pp3_supporting() -> f64 { 0.644 }
fn default_pp3_moderate() -> f64 { 0.773 }
fn default_pp3_strong() -> f64 { 0.932 }
fn default_bp4_supporting() -> f64 { 0.290 }
fn default_bp4_moderate() -> f64 { 0.183 }
fn default_bp4_strong() -> f64 { 0.016 }
fn default_bp4_very_strong() -> f64 { 0.003 }
fn default_spliceai_pathogenic() -> f64 { 0.2 }
fn default_spliceai_benign() -> f64 { 0.1 }
fn default_phylop() -> f64 { 2.0 }
fn default_pli() -> f64 { 0.9 }
fn default_loeuf() -> f64 { 0.35 }
fn default_misz() -> f64 { 3.09 }
fn default_pm1_window() -> u64 { 5 }
fn default_pm1_threshold() -> u32 { 3 }
fn default_true() -> bool { true }
fn default_min_an() -> u64 { 2000 }

/// Default BA1 exception list (Ghosh et al. 2018, Hum Mutat — 9 variants).
fn default_ba1_exceptions() -> Vec<Ba1Exception> {
    let mk = |gene: &str, hgvs: &str, reason: &str| Ba1Exception {
        gene: gene.to_string(),
        hgvs_c: hgvs.to_string(),
        reason: Some(reason.to_string()),
    };
    vec![
        mk("ACAD9", "c.-44_-41dupTAAG", "Ghosh 2018 BA1 exception (VUS)"),
        mk("GJB2", "c.109G>A", "Ghosh 2018 BA1 exception (Pathogenic) — DFNB1 hearing loss"),
        mk("HFE", "c.187C>G", "Ghosh 2018 BA1 exception (Pathogenic) — hereditary hemochromatosis"),
        mk("HFE", "c.845G>A", "Ghosh 2018 BA1 exception (Pathogenic) — hereditary hemochromatosis (p.Cys282Tyr)"),
        mk("MEFV", "c.1105C>T", "Ghosh 2018 BA1 exception (VUS)"),
        mk("MEFV", "c.1223G>A", "Ghosh 2018 BA1 exception (VUS)"),
        mk("PIBF1", "c.1214G>A", "Ghosh 2018 BA1 exception (VUS)"),
        mk("ACADS", "c.511C>T", "Ghosh 2018 BA1 exception (VUS)"),
        mk("BTD", "c.1330G>C", "Ghosh 2018 BA1 exception (Pathogenic) — biotinidase deficiency"),
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default_config() {
        let cfg = AcmgConfig::default();
        assert!((cfg.ba1_af_threshold - 0.05).abs() < 1e-10);
        assert!(cfg.pm2_downgrade_to_supporting);
        assert!(!cfg.use_pp5_bp6);
    }

    #[test]
    fn test_gene_override() {
        let mut cfg = AcmgConfig::default();
        cfg.gene_overrides.insert(
            "BRCA1".to_string(),
            GeneOverride {
                mechanism: Some("LOF".to_string()),
                bs1_af_threshold: Some(0.001),
                pm2_af_threshold: None,
                disabled_criteria: vec![],
                strength_overrides: HashMap::new(),
                disorders: HashMap::new(),
            },
        );
        assert_eq!(cfg.effective_bs1_threshold(Some("BRCA1")), 0.001);
        assert_eq!(cfg.effective_bs1_threshold(Some("TP53")), 0.01);
        assert_eq!(cfg.effective_bs1_threshold(None), 0.01);
    }

    #[test]
    fn test_toml_deserialization() {
        let toml_str = r#"
ba1_af_threshold = 0.05
bs1_af_threshold = 0.005
pm2_af_threshold = 0.00005

[gene_overrides.BRCA1]
mechanism = "LOF"
bs1_af_threshold = 0.001
disabled_criteria = ["BP1"]
"#;
        let cfg: AcmgConfig = toml::from_str(toml_str).unwrap();
        assert!((cfg.bs1_af_threshold - 0.005).abs() < 1e-10);
        assert!(cfg.is_criterion_disabled("BRCA1", "BP1"));
        assert!(!cfg.is_criterion_disabled("BRCA1", "PM2"));
    }
}
