use anyhow::{Context, Result};
use std::collections::HashMap;
use std::fs;
use std::path::Path;

/// Parsed VEP cache info.txt metadata.
#[derive(Debug, Clone)]
pub struct CacheInfo {
    pub species: String,
    pub assembly: String,
    pub version: Option<String>,
    pub sift: Option<String>,
    pub polyphen: Option<String>,
    pub regulatory: bool,
    pub variation_cols: Vec<String>,
    pub var_type: Option<String>,
    pub serialiser_type: Option<String>,
    pub valid_chromosomes: Vec<String>,
    pub version_data: HashMap<String, String>,
    pub raw: HashMap<String, String>,
}

impl CacheInfo {
    /// Parse an info.txt file from a VEP cache directory.
    pub fn from_file(path: &Path) -> Result<Self> {
        let contents =
            fs::read_to_string(path).with_context(|| format!("Reading {}", path.display()))?;
        Self::parse(&contents)
    }

    /// Parse info.txt content.
    pub fn parse(contents: &str) -> Result<Self> {
        let mut raw = HashMap::new();

        for line in contents.lines() {
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            if let Some((key, value)) = line.split_once('\t') {
                raw.insert(key.to_string(), value.to_string());
            }
        }

        let species = raw.get("species").cloned().unwrap_or_default();
        let assembly = raw.get("assembly").cloned().unwrap_or_default();
        let version = raw.get("version").cloned();
        let sift = raw.get("sift").cloned();
        let polyphen = raw.get("polyphen").cloned();
        let regulatory = raw.get("regulatory").map(|v| v == "1").unwrap_or(false);
        let var_type = raw.get("var_type").cloned();
        let serialiser_type = raw.get("serialiser_type").cloned();

        let variation_cols = raw
            .get("variation_cols")
            .map(|v| v.split(',').map(|s| s.to_string()).collect())
            .unwrap_or_default();

        let valid_chromosomes = raw
            .get("valid_chromosomes")
            .map(|v| v.split(',').map(|s| s.trim().to_string()).collect())
            .unwrap_or_default();

        let mut version_data = HashMap::new();
        for (key, value) in &raw {
            if key.starts_with("source_") {
                version_data.insert(key.clone(), value.clone());
            }
        }

        Ok(Self {
            species,
            assembly,
            version,
            sift,
            polyphen,
            regulatory,
            variation_cols,
            var_type,
            serialiser_type,
            valid_chromosomes,
            version_data,
            raw,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_info() {
        let content = "\
species\thomo_sapiens
assembly\tGRCh38
version\t115
sift\tb
polyphen\tb
regulatory\t1
var_type\ttabix
serialiser_type\tstorable
variation_cols\tvariation_name,failed,somatic,start,end,allele_string,strand,minor_allele,minor_allele_freq
valid_chromosomes\t1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y,MT
source_gencode\tGENCODE 46
source_assembly\tGRCh38.p14";

        let info = CacheInfo::parse(content).unwrap();
        assert_eq!(info.species, "homo_sapiens");
        assert_eq!(info.assembly, "GRCh38");
        assert_eq!(info.version.as_deref(), Some("115"));
        assert_eq!(info.sift.as_deref(), Some("b"));
        assert_eq!(info.polyphen.as_deref(), Some("b"));
        assert!(info.regulatory);
        assert_eq!(info.var_type.as_deref(), Some("tabix"));
        assert_eq!(info.variation_cols.len(), 9);
        assert_eq!(info.variation_cols[0], "variation_name");
        assert_eq!(info.valid_chromosomes.len(), 25);
        assert_eq!(info.valid_chromosomes[0], "1");
        assert_eq!(
            info.version_data.get("source_gencode").unwrap(),
            "GENCODE 46"
        );
    }
}
