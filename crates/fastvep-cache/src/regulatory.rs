//! Regulatory region annotations.
//!
//! Parses Ensembl regulatory GFF3 features and provides interval-based
//! lookup for regulatory_region_variant consequences.

use anyhow::Result;
use std::collections::HashMap;
use std::io::BufRead;

/// A regulatory feature from Ensembl's regulatory build.
#[derive(Debug, Clone)]
pub struct RegulatoryFeature {
    pub stable_id: String,
    pub feature_type: RegulatoryType,
    pub chromosome: String,
    pub start: u64,
    pub end: u64,
}

/// Types of regulatory features.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RegulatoryType {
    Promoter,
    Enhancer,
    CtcfBindingSite,
    TfBindingSite,
    OpenChromatinRegion,
    PromoterFlankingRegion,
    Other,
}

impl RegulatoryType {
    pub fn from_str(s: &str) -> Self {
        match s.to_lowercase().as_str() {
            "promoter" => Self::Promoter,
            "enhancer" => Self::Enhancer,
            "ctcf_binding_site" | "ctcf" => Self::CtcfBindingSite,
            "tf_binding_site" | "tfbs" => Self::TfBindingSite,
            "open_chromatin_region" | "open_chromatin" => Self::OpenChromatinRegion,
            "promoter_flanking_region" => Self::PromoterFlankingRegion,
            _ => Self::Other,
        }
    }

    pub fn as_str(self) -> &'static str {
        match self {
            Self::Promoter => "promoter",
            Self::Enhancer => "enhancer",
            Self::CtcfBindingSite => "CTCF_binding_site",
            Self::TfBindingSite => "TF_binding_site",
            Self::OpenChromatinRegion => "open_chromatin_region",
            Self::PromoterFlankingRegion => "promoter_flanking_region",
            Self::Other => "regulatory_region",
        }
    }
}

/// In-memory regulatory feature index for interval-based lookup.
pub struct RegulatoryProvider {
    /// Chromosome -> sorted list of regulatory features.
    by_chrom: HashMap<String, Vec<RegulatoryFeature>>,
}

impl RegulatoryProvider {
    /// Build from a list of features.
    pub fn new(mut features: Vec<RegulatoryFeature>) -> Self {
        let mut by_chrom: HashMap<String, Vec<RegulatoryFeature>> = HashMap::new();
        for f in features.drain(..) {
            by_chrom.entry(f.chromosome.clone()).or_default().push(f);
        }
        for feats in by_chrom.values_mut() {
            feats.sort_by_key(|f| f.start);
        }
        Self { by_chrom }
    }

    /// Find all regulatory features overlapping [start, end] on a chromosome.
    pub fn find_overlapping(&self, chrom: &str, start: u64, end: u64) -> Vec<&RegulatoryFeature> {
        let feats = match self.by_chrom.get(chrom) {
            Some(f) => f,
            None => return Vec::new(),
        };

        let mut results = Vec::new();
        for f in feats {
            if f.start > end {
                break;
            }
            if f.end >= start {
                results.push(f);
            }
        }
        results
    }

    pub fn feature_count(&self) -> usize {
        self.by_chrom.values().map(|v| v.len()).sum()
    }
}

/// Parse Ensembl regulatory build GFF3 into RegulatoryFeatures.
pub fn parse_regulatory_gff3<R: BufRead>(reader: R) -> Result<Vec<RegulatoryFeature>> {
    let mut features = Vec::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') || line.is_empty() {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 9 {
            continue;
        }

        let feature_type = fields[2];
        if feature_type != "regulatory_region"
            && feature_type != "enhancer"
            && feature_type != "promoter"
            && !feature_type.contains("binding_site")
        {
            continue;
        }

        let chrom = fields[0].to_string();
        let start: u64 = match fields[3].parse() {
            Ok(s) => s,
            Err(_) => continue,
        };
        let end: u64 = match fields[4].parse() {
            Ok(e) => e,
            Err(_) => continue,
        };
        // GFF3 coordinates are 1-based and inclusive: start=0 or an inverted
        // range parses fine as u64 but is invalid, and (as with the
        // transcript GFF3 parser in gff.rs) risks colliding with everything
        // at position 0 or underflowing downstream `- 1` offset math.
        if start == 0 || start > end {
            log::warn!(
                "Regulatory GFF3: skipping line with invalid coordinates (start={}, end={})",
                start,
                end
            );
            continue;
        }

        // Parse attributes for ID and feature_type
        let attrs = fields[8];
        let mut stable_id = String::new();
        let mut reg_type = RegulatoryType::Other;

        for attr in attrs.split(';') {
            let attr = attr.trim();
            if let Some(id) = attr.strip_prefix("ID=") {
                stable_id = id.to_string();
            } else if let Some(ft) = attr.strip_prefix("feature_type=") {
                reg_type = RegulatoryType::from_str(ft);
            } else if let Some(bt) = attr.strip_prefix("biotype=") {
                if reg_type == RegulatoryType::Other {
                    reg_type = RegulatoryType::from_str(bt);
                }
            }
        }

        if reg_type == RegulatoryType::Other && feature_type != "regulatory_region" {
            reg_type = RegulatoryType::from_str(feature_type);
        }

        features.push(RegulatoryFeature {
            stable_id,
            feature_type: reg_type,
            chromosome: chrom,
            start,
            end,
        });
    }

    Ok(features)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_regulatory_provider() {
        let features = vec![
            RegulatoryFeature {
                stable_id: "ENSR00001".into(),
                feature_type: RegulatoryType::Promoter,
                chromosome: "chr1".into(),
                start: 1000,
                end: 2000,
            },
            RegulatoryFeature {
                stable_id: "ENSR00002".into(),
                feature_type: RegulatoryType::Enhancer,
                chromosome: "chr1".into(),
                start: 3000,
                end: 4000,
            },
        ];

        let provider = RegulatoryProvider::new(features);
        assert_eq!(provider.feature_count(), 2);

        let hits = provider.find_overlapping("chr1", 1500, 1600);
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].stable_id, "ENSR00001");

        let hits = provider.find_overlapping("chr1", 2500, 3500);
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].stable_id, "ENSR00002");

        let hits = provider.find_overlapping("chr1", 5000, 6000);
        assert_eq!(hits.len(), 0);
    }

    #[test]
    fn test_parse_regulatory_gff3() {
        let gff = "\
##gff-version 3
chr1\tensembl_regulation\tregulatory_region\t1000\t2000\t.\t.\t.\tID=ENSR00001;feature_type=promoter
chr1\tensembl_regulation\tregulatory_region\t3000\t4000\t.\t.\t.\tID=ENSR00002;feature_type=enhancer
";
        let features = parse_regulatory_gff3(gff.as_bytes()).unwrap();
        assert_eq!(features.len(), 2);
        assert_eq!(features[0].feature_type, RegulatoryType::Promoter);
        assert_eq!(features[1].feature_type, RegulatoryType::Enhancer);
    }

    #[test]
    fn test_parse_regulatory_gff3_skips_invalid_coordinates() {
        // start=0 and start>end parse fine as u64 but are invalid 1-based
        // GFF3 coordinates (same bug class as fastvep-cache::gff.rs).
        let gff = "\
##gff-version 3
chr1\tensembl_regulation\tregulatory_region\t0\t2000\t.\t.\t.\tID=ENSR00001;feature_type=promoter
chr1\tensembl_regulation\tregulatory_region\t4000\t3000\t.\t.\t.\tID=ENSR00002;feature_type=enhancer
chr1\tensembl_regulation\tregulatory_region\t5000\t6000\t.\t.\t.\tID=ENSR00003;feature_type=promoter
";
        let features = parse_regulatory_gff3(gff.as_bytes()).unwrap();
        assert_eq!(features.len(), 1);
        assert_eq!(features[0].stable_id, "ENSR00003");
    }
}
