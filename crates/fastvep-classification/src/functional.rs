//! Curated functional-assay evidence for PS3 and BS3.
//!
//! PS3 and BS3 are the two criteria fastVEP cannot compute. They ask whether a
//! well-established in vitro or in vivo assay showed the variant to damage, or
//! not damage, gene function - a question answered in a paper, not in a
//! database. fastVEP will not mine literature for them: a wrong PMID in a
//! clinical report is worse than no evidence at all. What it will do is
//! consume a curated file, which is how every VCEP pipeline works.
//!
//! The file is a TSV keyed by genomic coordinate:
//!
//! ```text
//! #chrom  pos       ref  alt  criterion  strength    pmid      note
//! chr15   88855485  A    G    PS3        Strong      29625052  minigene shows exon skipping
//! chr2    47478343  G    A    BS3        Supporting  31391288  normal MMR activity in vitro
//! ```
//!
//! `strength` is optional and defaults to Strong, the criteria's own default.
//! It is settable because Brnich 2020 (ClinGen SVI functional-evidence
//! recommendation) makes assay strength a judgement about the assay - its
//! validity, controls and dynamic range - rather than a constant. A curator
//! who has read the paper is the one who knows whether it supports Strong or
//! only Supporting.
//!
//! ## Why this outranks the predictors
//!
//! The SVI's ordering is explicit: functional evidence is stronger than
//! computational prediction, and PP3/BP4/BP7 must not be used to argue against
//! sound experimental data. So an entry here does more than add a criterion -
//! it suppresses the computational criteria for that variant, which is handled
//! in the reconciliation pass. The round-2 review's OCA2 c.1503A>G is the case
//! that motivated it: a synonymous variant with published splice-defect data,
//! where fastVEP was firing BP7 and BP4 and calling it benign.

use std::collections::HashMap;
use std::io::BufRead;

use crate::types::EvidenceStrength;

/// One curated functional-assay result.
#[derive(Debug, Clone, PartialEq)]
pub struct FunctionalEvidence {
    /// `PS3` (damaging) or `BS3` (not damaging).
    pub criterion: FunctionalCriterion,
    /// Strength the curator assigned to the assay (Brnich 2020).
    pub strength: EvidenceStrength,
    /// Supporting publication, surfaced in the criterion summary so a reader
    /// can go straight to the source.
    pub pmid: Option<String>,
    /// Free-text note from the curator, also surfaced in the summary.
    pub note: Option<String>,
}

/// Which of the two functional criteria an entry asserts.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FunctionalCriterion {
    /// Assay shows a damaging effect.
    Ps3,
    /// Assay shows no damaging effect.
    Bs3,
}

impl FunctionalCriterion {
    pub fn code(&self) -> &'static str {
        match self {
            Self::Ps3 => "PS3",
            Self::Bs3 => "BS3",
        }
    }
}

/// Genomic key for a functional-evidence entry.
///
/// Coordinates rather than HGVS deliberately: HGVS is transcript-relative, so
/// the same curated result would need re-keying whenever the pipeline picks a
/// different transcript, and normalising it back to a genomic position is the
/// work this file exists to avoid.
type VariantKey = (String, u64, String, String);

/// A loaded functional-evidence file, ready for per-variant lookup.
#[derive(Debug, Clone, Default)]
pub struct FunctionalEvidenceIndex {
    entries: HashMap<VariantKey, FunctionalEvidence>,
}

/// Strip a leading `chr` so `chr7` and `7` are the same contig. VCFs disagree
/// about the prefix and a curated file should not have to guess which
/// convention the caller's VCF uses.
fn normalize_chrom(chrom: &str) -> String {
    chrom
        .strip_prefix("chr")
        .unwrap_or(chrom)
        .to_ascii_uppercase()
}

fn parse_strength(s: &str) -> Result<EvidenceStrength, String> {
    match s.to_ascii_lowercase().replace(['_', '-', ' '], "").as_str() {
        "supporting" => Ok(EvidenceStrength::Supporting),
        "moderate" => Ok(EvidenceStrength::Moderate),
        "strong" => Ok(EvidenceStrength::Strong),
        "verystrong" => Ok(EvidenceStrength::VeryStrong),
        "standalone" => Ok(EvidenceStrength::Standalone),
        other => Err(format!(
            "unknown strength {:?} (expected Supporting, Moderate, Strong, VeryStrong or Standalone)",
            other
        )),
    }
}

impl FunctionalEvidenceIndex {
    /// Parse a functional-evidence TSV.
    ///
    /// Blank lines and `#` comments are skipped, so the documented header line
    /// is valid input. Every other malformed line is an error rather than a
    /// silent skip: this file is hand-curated clinical evidence, and a typo
    /// that quietly dropped a row would show up as a variant mysteriously
    /// missing its PS3 rather than as a message.
    pub fn from_reader<R: BufRead>(reader: R) -> Result<Self, String> {
        let mut entries: HashMap<VariantKey, FunctionalEvidence> = HashMap::new();

        for (lineno, line) in reader.lines().enumerate() {
            let line = line.map_err(|e| format!("line {}: {}", lineno + 1, e))?;
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let f: Vec<&str> = line.split('\t').map(str::trim).collect();
            let at = |i: usize| -> &str { f.get(i).copied().unwrap_or("") };
            if f.len() < 5 {
                return Err(format!(
                    "line {}: expected at least 5 tab-separated columns (chrom, pos, ref, alt, criterion), found {}",
                    lineno + 1,
                    f.len()
                ));
            }

            let pos: u64 = at(1)
                .parse()
                .map_err(|_| format!("line {}: position {:?} is not a number", lineno + 1, at(1)))?;

            let criterion = match at(4).to_ascii_uppercase().as_str() {
                "PS3" => FunctionalCriterion::Ps3,
                "BS3" => FunctionalCriterion::Bs3,
                other => {
                    return Err(format!(
                        "line {}: criterion {:?} is not PS3 or BS3. fastVEP consumes functional-assay evidence only; other criteria are computed or come from the trio/VCF inputs",
                        lineno + 1,
                        other
                    ))
                }
            };

            let strength = match at(5) {
                "" | "." => EvidenceStrength::Strong,
                s => parse_strength(s).map_err(|e| format!("line {}: {}", lineno + 1, e))?,
            };

            let opt = |s: &str| match s {
                "" | "." => None,
                v => Some(v.to_string()),
            };

            let key = (
                normalize_chrom(at(0)),
                pos,
                at(2).to_ascii_uppercase(),
                at(3).to_ascii_uppercase(),
            );
            let entry = FunctionalEvidence {
                criterion,
                strength,
                pmid: opt(at(6)),
                note: opt(at(7)),
            };

            // A file asserting both PS3 and BS3 for one variant is a curation
            // error that would otherwise resolve by whichever line came last.
            if let Some(prev) = entries.get(&key) {
                if prev.criterion != entry.criterion {
                    return Err(format!(
                        "line {}: {}:{} {}>{} is given both {} and {}; a variant cannot have functional evidence in both directions",
                        lineno + 1,
                        key.0,
                        key.1,
                        key.2,
                        key.3,
                        prev.criterion.code(),
                        entry.criterion.code(),
                    ));
                }
            }
            entries.insert(key, entry);
        }

        Ok(Self { entries })
    }

    /// Load from a path.
    pub fn from_file(path: &std::path::Path) -> Result<Self, String> {
        let file = std::fs::File::open(path)
            .map_err(|e| format!("opening functional evidence file {}: {}", path.display(), e))?;
        Self::from_reader(std::io::BufReader::new(file))
            .map_err(|e| format!("{}: {}", path.display(), e))
    }

    /// Resolve the entry for one ALT allele of a VCF record.
    ///
    /// Keyed on the record's *original* VCF coordinates, not fastVEP's
    /// normalised ones: the curated file is written by a human reading a VCF,
    /// so `chr1 55039878 CACTGCTGCT C` is what they will have typed, and
    /// silently requiring the trimmed form would make entries fail to match
    /// with no diagnostic.
    ///
    /// `alt_index` selects from a multi-allelic ALT column, so the curated
    /// result for one allele of a site is not applied to its neighbours.
    pub fn for_vcf_allele(
        &self,
        chrom: &str,
        pos: u64,
        reference: &str,
        alt_column: &str,
        alt_index: usize,
    ) -> Option<&FunctionalEvidence> {
        let alt = alt_column.split(',').nth(alt_index)?;
        self.get(chrom, pos, reference, alt)
    }

    /// Look up one variant. `chrom` may carry a `chr` prefix or not.
    pub fn get(&self, chrom: &str, pos: u64, reference: &str, alt: &str) -> Option<&FunctionalEvidence> {
        self.entries.get(&(
            normalize_chrom(chrom),
            pos,
            reference.to_ascii_uppercase(),
            alt.to_ascii_uppercase(),
        ))
    }

    pub fn len(&self) -> usize {
        self.entries.len()
    }

    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn parse(s: &str) -> Result<FunctionalEvidenceIndex, String> {
        FunctionalEvidenceIndex::from_reader(s.as_bytes())
    }

    #[test]
    fn test_parses_the_documented_format() {
        let idx = parse(
            "#chrom\tpos\tref\talt\tcriterion\tstrength\tpmid\tnote\n\
             chr15\t88855485\tA\tG\tPS3\tStrong\t29625052\tminigene shows exon skipping\n\
             2\t47478343\tG\tA\tBS3\tSupporting\t31391288\tnormal MMR activity\n",
        )
        .expect("must parse");
        assert_eq!(idx.len(), 2);

        let ps3 = idx.get("chr15", 88855485, "A", "G").expect("PS3 entry");
        assert_eq!(ps3.criterion, FunctionalCriterion::Ps3);
        assert_eq!(ps3.strength, EvidenceStrength::Strong);
        assert_eq!(ps3.pmid.as_deref(), Some("29625052"));

        let bs3 = idx.get("chr2", 47478343, "G", "A").expect("BS3 entry");
        assert_eq!(bs3.criterion, FunctionalCriterion::Bs3);
        assert_eq!(bs3.strength, EvidenceStrength::Supporting);
    }

    #[test]
    fn test_chrom_prefix_and_case_do_not_matter() {
        // The curated file and the user's VCF are written by different people.
        let idx = parse("7\t100\ta\tg\tPS3\n").unwrap();
        assert!(idx.get("chr7", 100, "A", "G").is_some());
        assert!(idx.get("7", 100, "A", "G").is_some());
        assert!(idx.get("chr7", 100, "a", "g").is_some());
    }

    #[test]
    fn test_strength_defaults_to_strong_and_is_optional() {
        let idx = parse("1\t10\tA\tT\tPS3\n1\t20\tA\tT\tBS3\t.\n").unwrap();
        assert_eq!(idx.get("1", 10, "A", "T").unwrap().strength, EvidenceStrength::Strong);
        assert_eq!(idx.get("1", 20, "A", "T").unwrap().strength, EvidenceStrength::Strong);
    }

    #[test]
    fn test_curator_can_downgrade_a_weak_assay() {
        // Brnich 2020: assay strength is a judgement about the assay, so a
        // curator must be able to say Supporting.
        let idx = parse("1\t10\tA\tT\tPS3\tModerate\n").unwrap();
        assert_eq!(idx.get("1", 10, "A", "T").unwrap().strength, EvidenceStrength::Moderate);
    }

    #[test]
    fn test_malformed_lines_are_errors_not_silent_skips() {
        for (bad, expect) in [
            ("1\t10\tA\tT\n", "at least 5"),
            ("1\tNOTANUMBER\tA\tT\tPS3\n", "not a number"),
            ("1\t10\tA\tT\tPM3\n", "not PS3 or BS3"),
            ("1\t10\tA\tT\tPS3\tVeryWeak\n", "unknown strength"),
        ] {
            let err = parse(bad).expect_err("must reject");
            assert!(err.contains(expect), "{bad:?} gave {err:?}, wanted {expect:?}");
        }
    }

    #[test]
    fn test_contradictory_entries_for_one_variant_are_rejected() {
        let err = parse("1\t10\tA\tT\tPS3\n1\t10\tA\tT\tBS3\n")
            .expect_err("PS3 and BS3 for one variant is a curation error");
        assert!(err.contains("both directions"), "got: {err}");
    }

    #[test]
    fn test_comments_and_blank_lines_are_skipped() {
        let idx = parse("# a header\n\n1\t10\tA\tT\tPS3\n\n# trailing\n").unwrap();
        assert_eq!(idx.len(), 1);
    }

    #[test]
    fn test_absent_variant_looks_up_to_none() {
        let idx = parse("1\t10\tA\tT\tPS3\n").unwrap();
        assert!(idx.get("1", 11, "A", "T").is_none());
        assert!(idx.get("2", 10, "A", "T").is_none());
        assert!(idx.get("1", 10, "A", "C").is_none());
    }
}
