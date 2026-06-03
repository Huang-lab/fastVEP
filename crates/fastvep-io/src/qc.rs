//! Rules-based QC classification from variant-level INFO fields.
//!
//! Reads a TOML file describing ordered classes. Each class names a set of
//! numeric thresholds against VCF INFO fields (e.g. `DP`, `QD`, `MQ`, `AF`,
//! `FS`, `SOR`, …). A variant is assigned the *first* class whose every
//! threshold is satisfied; otherwise it falls through to `fallback` (default
//! `"FAIL_QC"`).
//!
//! The classifier is variant-level by design: it does *not* parse per-sample
//! FORMAT columns. This keeps it O(1) extra work per variant regardless of
//! how many samples the VCF carries.
//!
//! ## Example
//!
//! ```toml
//! fallback = "FAIL_QC"
//!
//! [[class]]
//! name = "HIGH_QC"
//! min_dp = 15
//! min_qd = 20
//! min_mq = 40
//! max_fs = 60.0
//!
//! [[class]]
//! name = "LOW_QC"
//! min_dp = 8
//! min_qd = 10
//! ```
//!
//! Any field can also be addressed by raw INFO ID via `[class.min]` /
//! `[class.max]` tables, e.g. `min = { DP = 15, "MQRankSum" = -12.5 }`.

use anyhow::{Context, Result};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::path::Path;

/// Default class name returned when no rule matches.
const DEFAULT_FALLBACK: &str = "FAIL_QC";

/// One ordered class in the rule set.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QcClass {
    /// Class label emitted in the `QC_CLASS` column.
    pub name: String,

    // ── Common GATK-style shortcuts. All optional. ──
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub min_dp: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub max_dp: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub min_qd: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub min_mq: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub min_af: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub max_af: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub max_fs: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub max_sor: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub min_mqranksum: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub min_readposranksum: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub min_qual: Option<f64>,

    /// Free-form `field -> threshold` map. Values are *minimum* required
    /// values for the named INFO field. Use to extend beyond the shortcuts.
    #[serde(default, skip_serializing_if = "HashMap::is_empty")]
    pub min: HashMap<String, f64>,

    /// Free-form `field -> threshold` map. Values are *maximum* allowed
    /// values for the named INFO field.
    #[serde(default, skip_serializing_if = "HashMap::is_empty")]
    pub max: HashMap<String, f64>,

    /// Require the VCF FILTER column to equal one of these values
    /// (case-sensitive). When unset, FILTER is ignored.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub require_filter: Vec<String>,
}

/// Loaded QC rule set.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QcRules {
    /// Ordered classes. The first one whose thresholds all pass wins.
    #[serde(rename = "class", default)]
    pub classes: Vec<QcClass>,
    /// Label assigned when no class matches. Defaults to `"FAIL_QC"`.
    #[serde(default = "default_fallback")]
    pub fallback: String,
}

fn default_fallback() -> String {
    DEFAULT_FALLBACK.to_string()
}

impl QcRules {
    /// Load a QC rule set from a TOML file.
    pub fn from_toml_file(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref();
        let content = std::fs::read_to_string(path)
            .with_context(|| format!("Reading QC rules: {}", path.display()))?;
        let rules: Self = toml::from_str(&content)
            .with_context(|| format!("Parsing QC rules TOML: {}", path.display()))?;
        rules.validate()?;
        Ok(rules)
    }

    /// Sanity check: class names must be non-empty.
    fn validate(&self) -> Result<()> {
        for c in &self.classes {
            if c.name.trim().is_empty() {
                anyhow::bail!("QC class is missing a non-empty `name`");
            }
        }
        Ok(())
    }

    /// Classify a variant. `info` is the parsed VCF INFO column.
    /// `filter` is the VCF FILTER column (e.g. `"PASS"`).
    /// Returns the matched class name, or the fallback when none match.
    pub fn classify<'a>(&'a self, info: &InfoView<'_>, filter: &str) -> &'a str {
        for class in &self.classes {
            if class_matches(class, info, filter) {
                return &class.name;
            }
        }
        &self.fallback
    }
}

fn class_matches(class: &QcClass, info: &InfoView<'_>, filter: &str) -> bool {
    if !class.require_filter.is_empty() && !class.require_filter.iter().any(|f| f == filter) {
        return false;
    }

    let checks: &[(&str, Option<f64>, Cmp)] = &[
        ("DP", class.min_dp, Cmp::Min),
        ("DP", class.max_dp, Cmp::Max),
        ("QD", class.min_qd, Cmp::Min),
        ("MQ", class.min_mq, Cmp::Min),
        ("AF", class.min_af, Cmp::Min),
        ("AF", class.max_af, Cmp::Max),
        ("FS", class.max_fs, Cmp::Max),
        ("SOR", class.max_sor, Cmp::Max),
        ("MQRankSum", class.min_mqranksum, Cmp::Min),
        ("ReadPosRankSum", class.min_readposranksum, Cmp::Min),
        ("QUAL", class.min_qual, Cmp::Min),
    ];

    for (key, threshold, cmp) in checks {
        let Some(threshold) = *threshold else {
            continue;
        };
        let observed = if *key == "QUAL" {
            info.qual
        } else {
            info.numeric(key)
        };
        if !passes(observed, threshold, *cmp) {
            return false;
        }
    }

    for (key, threshold) in &class.min {
        if !passes(info.numeric(key), *threshold, Cmp::Min) {
            return false;
        }
    }
    for (key, threshold) in &class.max {
        if !passes(info.numeric(key), *threshold, Cmp::Max) {
            return false;
        }
    }

    true
}

#[derive(Clone, Copy)]
enum Cmp {
    Min,
    Max,
}

fn passes(observed: Option<f64>, threshold: f64, cmp: Cmp) -> bool {
    let Some(observed) = observed else {
        // A missing field is a non-match. Without the value, we can't claim
        // the threshold is satisfied, so the rule fails.
        return false;
    };
    match cmp {
        Cmp::Min => observed >= threshold,
        Cmp::Max => observed <= threshold,
    }
}

/// Cheap, allocation-free view over a VCF INFO column. Built once per row;
/// parsing is lazy: numeric() scans the INFO string in O(len).
///
/// Carries the parsed QUAL as a side channel so QC rules can reference it
/// alongside INFO fields (QUAL is its own VCF column, not INFO).
pub struct InfoView<'a> {
    info: &'a str,
    qual: Option<f64>,
}

impl<'a> InfoView<'a> {
    pub fn new(info: &'a str, qual: Option<f64>) -> Self {
        Self { info, qual }
    }

    /// Look up a numeric value by INFO key. Returns `None` when the key is
    /// missing, valueless (flag), or unparseable. For multi-valued fields
    /// (`A,B,C`), the *first* value is taken.
    pub fn numeric(&self, key: &str) -> Option<f64> {
        if self.info.is_empty() || self.info == "." {
            return None;
        }
        // INFO is `KEY=VAL;KEY=VAL;FLAG;...`. Hand-roll the scan to avoid
        // building a HashMap per row.
        let mut rest = self.info;
        while !rest.is_empty() {
            let (chunk, tail) = match rest.find(';') {
                Some(i) => (&rest[..i], &rest[i + 1..]),
                None => (rest, ""),
            };
            rest = tail;
            let (k, v) = match chunk.split_once('=') {
                Some(kv) => kv,
                None => continue,
            };
            if k != key {
                continue;
            }
            let first = v.split(',').next()?;
            return first.parse::<f64>().ok();
        }
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn parse(toml: &str) -> QcRules {
        let r: QcRules = toml::from_str(toml).unwrap();
        r.validate().unwrap();
        r
    }

    #[test]
    fn high_qc_when_all_thresholds_pass() {
        let rules = parse(
            r#"
            [[class]]
            name = "HIGH_QC"
            min_dp = 15
            min_qd = 20
            min_mq = 40
            "#,
        );
        let info = InfoView::new("DP=30;QD=25;MQ=60", None);
        assert_eq!(rules.classify(&info, "PASS"), "HIGH_QC");
    }

    #[test]
    fn falls_through_to_low_then_fail() {
        let rules = parse(
            r#"
            [[class]]
            name = "HIGH_QC"
            min_dp = 15
            min_qd = 20

            [[class]]
            name = "LOW_QC"
            min_dp = 8
            "#,
        );
        // QD too low → not HIGH_QC; DP ≥ 8 → LOW_QC.
        let info = InfoView::new("DP=10;QD=5", None);
        assert_eq!(rules.classify(&info, "PASS"), "LOW_QC");

        // DP too low → no class matches.
        let info = InfoView::new("DP=5;QD=5", None);
        assert_eq!(rules.classify(&info, "PASS"), "FAIL_QC");
    }

    #[test]
    fn custom_fallback() {
        let rules = parse(
            r#"
            fallback = "DROP"
            [[class]]
            name = "PASS_QC"
            min_dp = 15
            "#,
        );
        let info = InfoView::new("DP=1", None);
        assert_eq!(rules.classify(&info, "PASS"), "DROP");
    }

    #[test]
    fn max_threshold_enforced() {
        let rules = parse(
            r#"
            [[class]]
            name = "GOOD"
            max_fs = 60.0
            "#,
        );
        assert_eq!(rules.classify(&InfoView::new("FS=30.0", None), "PASS"), "GOOD");
        assert_eq!(rules.classify(&InfoView::new("FS=65.0", None), "PASS"), "FAIL_QC");
    }

    #[test]
    fn missing_field_means_no_match() {
        let rules = parse(
            r#"
            [[class]]
            name = "GOOD"
            min_dp = 10
            "#,
        );
        // No DP in INFO → cannot prove threshold passed.
        assert_eq!(rules.classify(&InfoView::new("AF=0.5", None), "PASS"), "FAIL_QC");
    }

    #[test]
    fn require_filter_pass() {
        let rules = parse(
            r#"
            [[class]]
            name = "GOOD"
            require_filter = ["PASS"]
            min_dp = 10
            "#,
        );
        let info = InfoView::new("DP=20", None);
        assert_eq!(rules.classify(&info, "PASS"), "GOOD");
        assert_eq!(rules.classify(&info, "LowQual"), "FAIL_QC");
    }

    #[test]
    fn freeform_min_max_tables() {
        let rules = parse(
            r#"
            [[class]]
            name = "GOOD"
            [class.min]
            "MQRankSum" = -12.5
            [class.max]
            "FS" = 60.0
            "#,
        );
        let info = InfoView::new("MQRankSum=-1.0;FS=10.0", None);
        assert_eq!(rules.classify(&info, "PASS"), "GOOD");
        let info = InfoView::new("MQRankSum=-20.0;FS=10.0", None);
        assert_eq!(rules.classify(&info, "PASS"), "FAIL_QC");
    }

    #[test]
    fn multi_value_takes_first() {
        let rules = parse(
            r#"
            [[class]]
            name = "GOOD"
            min_af = 0.2
            "#,
        );
        let info = InfoView::new("AF=0.5,0.05", None);
        assert_eq!(rules.classify(&info, "PASS"), "GOOD");
    }

    #[test]
    fn qual_threshold() {
        let rules = parse(
            r#"
            [[class]]
            name = "GOOD"
            min_qual = 30.0
            "#,
        );
        assert_eq!(
            rules.classify(&InfoView::new(".", Some(99.0)), "PASS"),
            "GOOD"
        );
        assert_eq!(
            rules.classify(&InfoView::new(".", Some(10.0)), "PASS"),
            "FAIL_QC"
        );
    }
}
