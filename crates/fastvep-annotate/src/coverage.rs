//! Which chromosomes the loaded databases hold, so that "no record" is not read as "absent".
//!
//! A database of some chromosomes only (one shard of a per-chromosome build, a region-restricted
//! build, a chromosome whose download failed and was not noticed) opens without complaint and answers
//! "no record" for every variant of the others. For gnomAD that is read as "absent from the
//! population", and PM2 fires: the `deploy-full.sh` of this project once built a gnomAD of chromosome
//! X alone, and a variant of chromosome 1 was called absent from a database that had never held
//! chromosome 1. The stack counted as complete, because a gnomAD database was loaded.
//!
//! The classifier is told, for each variant, whether the gnomAD databases that are loaded hold its
//! chromosome (`LoadedSources::covers_variant`), and the server decides whether a genome can be
//! classified from whether they hold the chromosomes that the genome's gene model has transcripts on.

use fastvep_cache::annotation::AnnotationProvider;
use fastvep_classification::{LoadedSources, Source};
use std::collections::HashSet;

/// A chromosome name as it is compared: without a `chr` prefix, upper-case, and the mitochondrial
/// genome as `MT`, so that `chr1`, `Chr1` and `1` are one chromosome and `M` and `MT` are one.
pub fn normalize_chromosome(name: &str) -> String {
    let bare = name
        .strip_prefix("chr")
        .or_else(|| name.strip_prefix("Chr"))
        .or_else(|| name.strip_prefix("CHR"))
        .unwrap_or(name);
    let upper = bare.to_ascii_uppercase();
    if upper == "M" {
        "MT".to_string()
    } else {
        upper
    }
}

/// Whether a chromosome is one of those a population database is expected to hold: an autosome
/// (a number) or X, by its normalized name. Y and the mitochondrial genome are left out of what a
/// genome must be covered for (the gnomAD of GRCh37 has no Y), and the variants on them are
/// declined one by one, as any variant on a chromosome that the database lacks is.
pub fn is_expected_chromosome(normalized: &str) -> bool {
    normalized == "X"
        || (!normalized.is_empty()
            && normalized.len() <= 3
            && normalized.bytes().all(|b| b.is_ascii_digit()))
}

/// The chromosomes the loaded gnomAD databases hold.
#[derive(Debug, Clone, Default)]
pub struct SourceCoverage {
    /// The normalized names of the chromosomes of every gnomAD database that is loaded, together.
    /// `None` when no gnomAD database is loaded (which is not a coverage question) or when one of
    /// them cannot say what it holds, which is read as holding anything.
    gnomad: Option<HashSet<String>>,
}

impl SourceCoverage {
    /// What the loaded databases cover. The shards of a per-chromosome build are many databases of
    /// one source, and what they cover is the union.
    pub fn of(providers: &[Box<dyn AnnotationProvider>]) -> Self {
        let mut covered = HashSet::new();
        let mut any = false;
        let mut unknown = false;
        for provider in providers {
            let key = [provider.json_key()];
            if !LoadedSources::from_json_keys(key).has(Source::Gnomad) {
                continue;
            }
            any = true;
            match provider.chromosomes() {
                Some(names) => covered.extend(names.iter().map(|n| normalize_chromosome(n))),
                None => unknown = true,
            }
        }
        Self {
            gnomad: (any && !unknown).then_some(covered),
        }
    }

    /// `loaded` for a variant on `chromosome`: the same sources, with gnomAD marked as holding
    /// nothing of it when the databases that are loaded do not hold it.
    pub fn sources_for(&self, loaded: LoadedSources, chromosome: &str) -> LoadedSources {
        match &self.gnomad {
            Some(covered) if !covered.contains(&normalize_chromosome(chromosome)) => {
                loaded.without_coverage(Source::Gnomad)
            }
            _ => loaded,
        }
    }

    /// The chromosomes of `wanted` (normalized) that the loaded gnomAD databases do not hold, in
    /// the order of `wanted`. Empty when they hold them all, when none is loaded, and when they
    /// cannot say.
    pub fn gnomad_lacks(&self, wanted: &[String]) -> Vec<String> {
        match &self.gnomad {
            Some(covered) => wanted
                .iter()
                .filter(|c| !covered.contains(*c))
                .cloned()
                .collect(),
            None => Vec::new(),
        }
    }
}

/// The chromosomes of a gene model that a population database has to cover for the genome to be
/// classified: the autosomes and X it has transcripts on, normalized, in natural order (1, 2, ...
/// 22, X).
pub fn expected_chromosomes<'a>(names: impl IntoIterator<Item = &'a str>) -> Vec<String> {
    let mut expected: Vec<String> = names
        .into_iter()
        .map(normalize_chromosome)
        .filter(|n| is_expected_chromosome(n))
        .collect();
    expected.sort_by_key(|n| match n.parse::<u32>() {
        Ok(number) => (0, number, String::new()),
        Err(_) => (1, 0, n.clone()),
    });
    expected.dedup();
    expected
}

/// A list of chromosomes for a sentence: `1, 2, 3, ... 22 and X`, with a run of numbers written as
/// a range once it is longer than three.
pub fn describe_chromosomes(chromosomes: &[String]) -> String {
    let numbers: Vec<u32> = chromosomes.iter().filter_map(|c| c.parse().ok()).collect();
    let others: Vec<&String> = chromosomes
        .iter()
        .filter(|c| c.parse::<u32>().is_err())
        .collect();
    let mut parts: Vec<String> = Vec::new();
    let mut i = 0;
    while i < numbers.len() {
        let mut j = i;
        while j + 1 < numbers.len() && numbers[j + 1] == numbers[j] + 1 {
            j += 1;
        }
        if j - i >= 3 {
            parts.push(format!("{}-{}", numbers[i], numbers[j]));
        } else {
            parts.extend(numbers[i..=j].iter().map(|n| n.to_string()));
        }
        i = j + 1;
    }
    parts.extend(others.iter().map(|c| c.to_string()));
    match parts.len() {
        0 => String::new(),
        1 => parts.remove(0),
        _ => {
            let last = parts.pop().unwrap_or_default();
            format!("{} and {}", parts.join(", "), last)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn names_are_compared_without_the_prefix_the_database_happens_to_use() {
        for (a, b) in [
            ("chr1", "1"),
            ("Chr1", "1"),
            ("chrX", "x"),
            ("M", "MT"),
            ("chrM", "MT"),
        ] {
            assert_eq!(normalize_chromosome(a), normalize_chromosome(b), "{a} {b}");
        }
        assert_ne!(normalize_chromosome("chr1"), normalize_chromosome("chr11"));
    }

    #[test]
    fn the_autosomes_and_x_are_what_a_genome_is_covered_for() {
        let expected = expected_chromosomes([
            "chr1",
            "2",
            "chr10",
            "X",
            "Y",
            "MT",
            "KI270728.1",
            "GL000194.1",
            "2L",
            "chr1",
        ]);
        assert_eq!(expected, ["1", "2", "10", "X"]);
    }

    #[test]
    fn chromosomes_are_listed_as_a_sentence() {
        let all: Vec<String> = (1..=22)
            .map(|n| n.to_string())
            .chain(["X".to_string()])
            .collect();
        assert_eq!(describe_chromosomes(&all), "1-22 and X");
        assert_eq!(
            describe_chromosomes(&["1".into(), "2".into(), "X".into()]),
            "1, 2 and X"
        );
        assert_eq!(describe_chromosomes(&["7".into()]), "7");
        assert_eq!(describe_chromosomes(&[]), "");
    }
}
