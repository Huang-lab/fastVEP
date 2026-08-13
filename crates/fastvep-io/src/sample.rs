//! Multi-sample VCF FORMAT/genotype parsing.
//!
//! Parses FORMAT and sample columns from VCF to extract per-sample
//! genotype, depth, quality, and allele depth information.

use std::collections::HashMap;

/// Parsed genotype for a sample.
#[derive(Debug, Clone)]
pub struct Genotype {
    /// Allele indices (0 = ref, 1 = first alt, etc.). None = missing ".".
    pub alleles: Vec<Option<u32>>,
    /// Whether the genotype is phased (true = '|', false = '/').
    pub phased: bool,
}

impl Genotype {
    /// Returns true if heterozygous (two different non-missing alleles).
    pub fn is_het(&self) -> bool {
        let present: Vec<u32> = self.alleles.iter().filter_map(|a| *a).collect();
        present.len() >= 2 && present.iter().min() != present.iter().max()
    }

    /// Returns true if homozygous reference (all alleles = 0).
    pub fn is_hom_ref(&self) -> bool {
        self.alleles.iter().all(|a| *a == Some(0))
    }

    /// Returns true if homozygous alternate (all non-missing alleles same, > 0).
    pub fn is_hom_alt(&self) -> bool {
        let present: Vec<u32> = self.alleles.iter().filter_map(|a| *a).collect();
        !present.is_empty()
            && present.iter().all(|a| *a > 0)
            && present.iter().min() == present.iter().max()
    }

    /// Returns true if all alleles are missing.
    pub fn is_missing(&self) -> bool {
        self.alleles.iter().all(|a| a.is_none())
    }

    /// Short string representation: "het", "hom_ref", "hom_alt", "missing".
    pub fn class(&self) -> &'static str {
        if self.is_missing() {
            "missing"
        } else if self.is_hom_ref() {
            "hom_ref"
        } else if self.is_hom_alt() {
            "hom_alt"
        } else if self.is_het() {
            "het"
        } else {
            "other"
        }
    }
}

/// Per-sample data extracted from VCF FORMAT columns.
#[derive(Debug, Clone)]
pub struct SampleData {
    /// Sample name (from #CHROM header).
    pub name: String,
    /// Parsed genotype (from GT field).
    pub genotype: Option<Genotype>,
    /// Read depth (from DP field).
    pub depth: Option<u32>,
    /// Genotype quality (from GQ field).
    pub quality: Option<u32>,
    /// Per-allele read depths (from AD field).
    pub allele_depths: Vec<u32>,
    /// All FORMAT fields as raw strings.
    pub raw_fields: HashMap<String, String>,
}

/// Parse the FORMAT string and sample columns into SampleData.
///
/// - `format_str`: the FORMAT column (e.g., "GT:DP:GQ:AD")
/// - `sample_strs`: each sample column (e.g., "0/1:30:99:15,15")
/// - `sample_names`: sample names from the #CHROM header
pub fn parse_samples(
    format_str: &str,
    sample_strs: &[&str],
    sample_names: &[String],
) -> Vec<SampleData> {
    let format_keys: Vec<&str> = format_str.split(':').collect();

    sample_strs
        .iter()
        .enumerate()
        .map(|(i, sample_str)| {
            let name = sample_names
                .get(i)
                .cloned()
                .unwrap_or_else(|| format!("SAMPLE_{}", i));

            let values: Vec<&str> = sample_str.split(':').collect();
            let mut raw_fields = HashMap::new();
            for (j, key) in format_keys.iter().enumerate() {
                if let Some(val) = values.get(j) {
                    raw_fields.insert(key.to_string(), val.to_string());
                }
            }

            let genotype = raw_fields.get("GT").map(|gt| parse_genotype(gt));
            let depth = raw_fields.get("DP").and_then(|v| v.parse().ok());
            let quality = raw_fields.get("GQ").and_then(|v| v.parse().ok());
            let allele_depths = raw_fields
                .get("AD")
                .map(|v| v.split(',').filter_map(|s| s.parse().ok()).collect())
                .unwrap_or_default();

            SampleData {
                name,
                genotype,
                depth,
                quality,
                allele_depths,
                raw_fields,
            }
        })
        .collect()
}

/// Parse a GT string like "0/1", "1|1", "./.", "0/1/2".
fn parse_genotype(gt: &str) -> Genotype {
    let phased = gt.contains('|');
    let sep = if phased { '|' } else { '/' };
    let alleles = gt
        .split(sep)
        .map(|s| if s == "." { None } else { s.parse().ok() })
        .collect();
    Genotype { alleles, phased }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_genotype() {
        let gt = parse_genotype("0/1");
        assert!(gt.is_het());
        assert!(!gt.phased);

        let gt = parse_genotype("0|0");
        assert!(gt.is_hom_ref());
        assert!(gt.phased);

        let gt = parse_genotype("1/1");
        assert!(gt.is_hom_alt());

        let gt = parse_genotype("./.");
        assert!(gt.is_missing());
    }

    #[test]
    fn test_parse_samples() {
        let names = vec!["proband".into(), "mother".into(), "father".into()];
        let samples = parse_samples(
            "GT:DP:GQ:AD",
            &["0/1:30:99:15,15", "0/0:25:80:25,0", "0/0:28:90:28,0"],
            &names,
        );

        assert_eq!(samples.len(), 3);
        assert_eq!(samples[0].name, "proband");
        assert!(samples[0].genotype.as_ref().unwrap().is_het());
        assert_eq!(samples[0].depth, Some(30));
        assert_eq!(samples[0].quality, Some(99));
        assert_eq!(samples[0].allele_depths, vec![15, 15]);

        assert_eq!(samples[1].name, "mother");
        assert!(samples[1].genotype.as_ref().unwrap().is_hom_ref());

        assert_eq!(samples[2].name, "father");
        assert!(samples[2].genotype.as_ref().unwrap().is_hom_ref());
    }

    #[test]
    fn test_genotype_class() {
        assert_eq!(parse_genotype("0/1").class(), "het");
        assert_eq!(parse_genotype("0/0").class(), "hom_ref");
        assert_eq!(parse_genotype("1/1").class(), "hom_alt");
        assert_eq!(parse_genotype("./.").class(), "missing");
    }
}
