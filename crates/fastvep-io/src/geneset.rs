//! Gene-set membership filter.
//!
//! Loads a panel of gene symbols and/or Ensembl gene IDs from a text file
//! (one identifier per line; blank lines and `#` comments ignored) into a
//! `HashSet`. The annotation pipeline uses it to drop transcript rows whose
//! gene is not in the set — the lookup is O(1) per transcript and the set
//! is loaded once at startup, so the per-variant hot path is unaffected
//! beyond a single membership check.

use anyhow::{Context, Result};
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// Set of gene symbols and Ensembl gene IDs to include in output. A row is
/// kept when its transcript's symbol OR gene ID matches.
#[derive(Debug, Clone, Default)]
pub struct GeneSet {
    members: HashSet<String>,
}

impl GeneSet {
    /// Build a `GeneSet` from an iterator of identifier strings.
    pub fn from_iter<I, S>(iter: I) -> Self
    where
        I: IntoIterator<Item = S>,
        S: Into<String>,
    {
        Self {
            members: iter
                .into_iter()
                .map(|s| s.into())
                .filter(|s| !s.is_empty())
                .collect(),
        }
    }

    /// Load a gene set from a text file. Each non-empty, non-`#`-prefixed
    /// line contributes one identifier. Whitespace is trimmed. Tab-separated
    /// files are supported: only the first column is read so panels exported
    /// from spreadsheets work directly.
    pub fn from_file(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref();
        let file = File::open(path)
            .with_context(|| format!("Opening gene list: {}", path.display()))?;
        let reader = BufReader::new(file);
        let mut members = HashSet::new();
        for line in reader.lines() {
            let line = line?;
            let trimmed = line.trim();
            if trimmed.is_empty() || trimmed.starts_with('#') {
                continue;
            }
            let first = trimmed.split('\t').next().unwrap_or("").trim();
            if first.is_empty() {
                continue;
            }
            members.insert(first.to_string());
        }
        Ok(Self { members })
    }

    /// Number of identifiers in the set.
    pub fn len(&self) -> usize {
        self.members.len()
    }

    /// True when the set has no members. Callers can short-circuit filtering.
    pub fn is_empty(&self) -> bool {
        self.members.is_empty()
    }

    /// Returns true when either the gene symbol or the Ensembl gene ID
    /// (with or without a trailing `.<version>` suffix) is in the set.
    pub fn contains_gene(&self, gene_id: &str, gene_symbol: Option<&str>) -> bool {
        if let Some(sym) = gene_symbol {
            if self.members.contains(sym) {
                return true;
            }
        }
        if self.members.contains(gene_id) {
            return true;
        }
        // Strip Ensembl version suffix and try again.
        if let Some(base) = gene_id.split('.').next() {
            if base != gene_id && self.members.contains(base) {
                return true;
            }
        }
        false
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn contains_by_symbol_or_gene_id() {
        let set = GeneSet::from_iter(["BRCA1", "ENSG00000139618"]);
        assert!(set.contains_gene("ENSG00000012048", Some("BRCA1")));
        assert!(set.contains_gene("ENSG00000139618", Some("BRCA2")));
        assert!(!set.contains_gene("ENSG00000141510", Some("TP53")));
    }

    #[test]
    fn version_suffix_stripped() {
        let set = GeneSet::from_iter(["ENSG00000012048"]);
        assert!(set.contains_gene("ENSG00000012048.23", None));
        assert!(set.contains_gene("ENSG00000012048", None));
    }

    #[test]
    fn from_file_ignores_blanks_and_comments() {
        let mut tmp = tempfile_for_test();
        writeln!(tmp, "# panel: DDR genes").unwrap();
        writeln!(tmp).unwrap();
        writeln!(tmp, "BRCA1").unwrap();
        writeln!(tmp, "BRCA2\tsome_note").unwrap();
        writeln!(tmp, "  TP53  ").unwrap();
        tmp.flush().unwrap();

        let set = GeneSet::from_file(tmp.path()).unwrap();
        assert_eq!(set.len(), 3);
        assert!(set.contains_gene("X", Some("BRCA1")));
        assert!(set.contains_gene("X", Some("BRCA2")));
        assert!(set.contains_gene("X", Some("TP53")));
        assert!(!set.contains_gene("X", Some("EGFR")));
    }

    fn tempfile_for_test() -> tempfile::NamedTempFile {
        tempfile::NamedTempFile::new().expect("tempfile")
    }
}
