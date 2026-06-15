//! Chromosome name aliasing and synonym resolution.
//!
//! Genome resources disagree about how to name contigs: Ensembl uses bare
//! names (`1`, `X`, `MT`), UCSC prefixes `chr` (`chr1`, `chrX`, `chrM`), and
//! NCBI RefSeq uses accessions (`NC_000001.11`). The first two differ by a
//! mechanical rule and are handled by [`chrom_aliases`]; the RefSeq accessions
//! have no algorithmic relationship to the others and can only be reconciled
//! with a mapping table — a VEP-style `chr_synonyms.txt`, parsed by
//! [`ChromSynonyms`].

use std::collections::HashMap;

/// Equivalent on-disk names for a chromosome, derived by mechanical rules.
///
/// Different upstream sources (and different user VCFs) mix `chr*` and bare
/// styles, and historically the SA writer canonicalized to the bare form
/// (`"1"`, `"X"`, `"MT"`) while modern inputs use `chr*`. This returns the
/// query name first, then all known aliases — so the reader can satisfy a
/// `chr1` query against an index built with `1`, or vice versa, without
/// requiring users to rebuild databases. See issue #37.
///
/// The returned `Vec` is small (1–3 entries) so caller-side iteration is
/// cheap. The first element is always the input name unchanged.
pub fn chrom_aliases(chrom: &str) -> Vec<String> {
    let mut out = Vec::with_capacity(3);
    out.push(chrom.to_string());

    // An empty input is most likely a programming error upstream — return
    // it unchanged so the caller still sees a single (empty) "alias" and
    // can surface a meaningful miss, rather than synthesizing a bogus
    // `"chr"` lookup.
    if chrom.is_empty() {
        return out;
    }

    // chr1 <-> 1
    if let Some(stripped) = chrom.strip_prefix("chr") {
        if !stripped.is_empty() && stripped != chrom {
            out.push(stripped.to_string());
        }
    } else {
        out.push(format!("chr{}", chrom));
    }

    // Mitochondrial special case: chrM / M / MT / chrMT all refer to the
    // same contig but UCSC uses chrM, NCBI uses MT.
    let mito_set = ["chrM", "M", "MT", "chrMT"];
    if mito_set.contains(&chrom) {
        for alt in mito_set {
            if alt != chrom && !out.iter().any(|n| n == alt) {
                out.push(alt.to_string());
            }
        }
    }

    out
}

/// Heuristic: does this name look like an NCBI RefSeq molecule accession
/// (`NC_000017.11`, `NW_…`, `NT_…`, `NG_…`)?
///
/// Used only to enrich error messages — when a contig like `NC_000017.11`
/// can't be found in a FASTA built with Ensembl bare names, the failure is
/// almost certainly a naming mismatch that `--synonyms` would fix, so we say
/// so. The pattern is `^N[CWTG]_\d+\.\d+$`; implemented by hand to avoid
/// pulling a regex dependency into the foundational crate.
pub fn looks_like_refseq_accession(name: &str) -> bool {
    let bytes = name.as_bytes();
    if bytes.len() < 6 {
        return false;
    }
    if bytes[0] != b'N' || !matches!(bytes[1], b'C' | b'W' | b'T' | b'G') || bytes[2] != b'_' {
        return false;
    }
    let rest = &name[3..];
    let Some((digits, version)) = rest.split_once('.') else {
        return false;
    };
    !digits.is_empty()
        && digits.bytes().all(|b| b.is_ascii_digit())
        && !version.is_empty()
        && version.bytes().all(|b| b.is_ascii_digit())
}

/// Bidirectional chromosome synonym table, loaded from a VEP-style
/// `chr_synonyms.txt`.
///
/// Each non-empty line of the file lists names that all refer to the same
/// contig, separated by whitespace/tabs (e.g. `1\tchr1\tNC_000001.11`). Every
/// name on a line maps to the full set of its line-mates, so a lookup from any
/// synonym yields all the others regardless of column order.
#[derive(Debug, Clone, Default)]
pub struct ChromSynonyms {
    /// name -> the other names sharing its line (the name itself excluded).
    map: HashMap<String, Vec<String>>,
}

impl ChromSynonyms {
    /// An empty table. Lookups still return the mechanical [`chrom_aliases`]
    /// fallbacks, so chr↔bare and mitochondrial forms resolve with zero
    /// configuration.
    pub fn new() -> Self {
        Self::default()
    }

    /// Parse a VEP `chr_synonyms.txt`. Lines are split on any whitespace;
    /// blank lines and single-token lines are ignored (nothing to map).
    pub fn parse(contents: &str) -> Self {
        let mut map: HashMap<String, Vec<String>> = HashMap::new();
        for line in contents.lines() {
            let names: Vec<&str> = line.split_whitespace().collect();
            if names.len() < 2 {
                continue;
            }
            for &name in &names {
                let entry = map.entry(name.to_string()).or_default();
                for &other in &names {
                    if other != name && !entry.iter().any(|n| n == other) {
                        entry.push(other.to_string());
                    }
                }
            }
        }
        Self { map }
    }

    /// True when no synonym lines were loaded.
    pub fn is_empty(&self) -> bool {
        self.map.is_empty()
    }

    /// All names equivalent to `query`, most-specific first: the query itself,
    /// then its file synonyms, then the mechanical [`chrom_aliases`]
    /// (chr↔bare, mito) of all of those. De-duplicated, query first.
    pub fn aliases(&self, query: &str) -> Vec<String> {
        let mut out: Vec<String> = Vec::new();
        let push = |name: String, out: &mut Vec<String>| {
            if !out.iter().any(|n| n == &name) {
                out.push(name);
            }
        };

        push(query.to_string(), &mut out);
        if let Some(syns) = self.map.get(query) {
            for s in syns {
                push(s.clone(), &mut out);
            }
        }
        // Fold in mechanical aliases for the query and each synonym so a
        // synonyms file that only lists `NC_000001.11 1` still satisfies a
        // `chr1` FASTA.
        for base in out.clone() {
            for alias in chrom_aliases(&base) {
                push(alias, &mut out);
            }
        }
        out
    }
}

#[cfg(test)]
mod chrom_alias_tests {
    use super::chrom_aliases;

    #[test]
    fn chr_prefix_round_trips() {
        let aliases = chrom_aliases("chr1");
        assert!(aliases.iter().any(|n| n == "chr1"));
        assert!(aliases.iter().any(|n| n == "1"));
    }

    #[test]
    fn bare_form_round_trips() {
        let aliases = chrom_aliases("1");
        assert!(aliases.iter().any(|n| n == "1"));
        assert!(aliases.iter().any(|n| n == "chr1"));
    }

    #[test]
    fn mitochondrial_aliases_cover_all_four_forms() {
        for name in ["chrM", "M", "MT", "chrMT"] {
            let aliases = chrom_aliases(name);
            for form in ["chrM", "M", "MT", "chrMT"] {
                assert!(
                    aliases.iter().any(|n| n == form),
                    "{} should resolve {}",
                    name,
                    form
                );
            }
        }
    }

    #[test]
    fn unknown_contig_returns_just_self_and_chr_variant() {
        let aliases = chrom_aliases("HLA-A*01:01");
        // Exactly two: the input plus its `chr`-prefixed form. Pinning the
        // length here keeps a future over-eager alias expansion from
        // silently broadening the lookup set.
        assert_eq!(aliases.len(), 2, "unexpected aliases: {:?}", aliases);
        assert_eq!(aliases[0], "HLA-A*01:01");
        assert_eq!(aliases[1], "chrHLA-A*01:01");
    }

    #[test]
    fn empty_input_does_not_synthesize_bogus_chr_alias() {
        // Regression: an earlier version pushed `format!("chr{}", "")` =
        // `"chr"` for empty input, which then collided with the synthetic
        // chr-strip case and could match unrelated index keys.
        let aliases = chrom_aliases("");
        assert_eq!(aliases, vec![String::new()]);
    }
}

#[cfg(test)]
mod synonym_tests {
    use super::{looks_like_refseq_accession, ChromSynonyms};

    #[test]
    fn refseq_accession_detected() {
        assert!(looks_like_refseq_accession("NC_000017.11"));
        assert!(looks_like_refseq_accession("NW_009646201.1"));
        assert!(looks_like_refseq_accession("NT_187361.1"));
        assert!(looks_like_refseq_accession("NG_012232.1"));
    }

    #[test]
    fn non_refseq_names_rejected() {
        for n in ["17", "chr17", "NC_00017", "NC_000017", "X", "", "NX_000017.11"] {
            assert!(!looks_like_refseq_accession(n), "{} misclassified", n);
        }
    }

    #[test]
    fn synonyms_map_bidirectionally() {
        let syn = ChromSynonyms::parse("17\tchr17\tNC_000017.11\n1 chr1\n");
        // From the RefSeq accession we reach the Ensembl bare name.
        assert!(syn.aliases("NC_000017.11").iter().any(|n| n == "17"));
        // From the bare name we reach the accession.
        assert!(syn.aliases("17").iter().any(|n| n == "NC_000017.11"));
        // Query itself comes first.
        assert_eq!(syn.aliases("NC_000017.11")[0], "NC_000017.11");
    }

    #[test]
    fn empty_table_still_yields_mechanical_aliases() {
        let syn = ChromSynonyms::new();
        assert!(syn.is_empty());
        assert!(syn.aliases("chr1").iter().any(|n| n == "1"));
    }

    #[test]
    fn single_token_lines_ignored() {
        let syn = ChromSynonyms::parse("17\n\nchr17\n");
        assert!(syn.is_empty());
    }

    #[test]
    fn partial_synonyms_fold_in_chr_prefix() {
        // File maps accession <-> bare only; chr-prefixed FASTA should still
        // resolve via the mechanical fallback.
        let syn = ChromSynonyms::parse("NC_000001.11\t1\n");
        assert!(syn.aliases("NC_000001.11").iter().any(|n| n == "chr1"));
    }
}
