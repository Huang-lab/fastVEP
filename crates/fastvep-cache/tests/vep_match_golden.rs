//! The known-variant matcher against Ensembl VEP's own answers.
//!
//! `fixtures/vep_match/golden.tsv` is what VEP's `get_matched_variant_alleles`
//! (ensembl-variation release/115, run as-is) returns for 2,000-odd pairs of an
//! input site and a cache record: the shapes issue #139 reports, the corners of
//! the Perl, and a seeded sample of anchored, shifted, multi-allelic,
//! other-strand and random alleles. The matcher has to reproduce each answer,
//! including *which* input allele matched *which* record allele and in what
//! order.
//!
//! The committed sample is small enough to run everywhere. To run the same
//! check on a much larger set, make one and point the second test at it. That
//! test is `#[ignore]`d, so that CI stays hermetic and says it did not run it,
//! and with `--ignored` and no file it stops and says what to set, where it used
//! to pass without having checked anything:
//!
//! ```text
//! crates/fastvep-cache/fixtures/vep_match/regen.sh big 200000 /tmp/pairs.tsv
//! FASTVEP_MATCH_CASES=/tmp/pairs.tsv cargo test -p fastvep-cache --test vep_match_golden -- --ignored
//! ```

use fastvep_cache::variation::match_alleles;
use std::path::Path;

/// `(a_allele, b_allele)` pairs, in the order VEP reports them.
type Pairs = Vec<(String, String)>;

/// `a_allele,b_allele` pairs, in order; empty for "-".
fn parse_expected(field: &str) -> Pairs {
    if field == "-" {
        return Vec::new();
    }
    field
        .split(';')
        .map(|pair| {
            let (a, b) = pair
                .split_once(',')
                .unwrap_or_else(|| panic!("malformed expected pair {pair:?}"));
            (a.to_string(), b.to_string())
        })
        .collect()
}

/// Our answer to one line, in the same shape as the expected one.
fn answer(line: &str) -> (Pairs, Pairs) {
    let f: Vec<&str> = line.split('\t').collect();
    assert_eq!(f.len(), 7, "a pair has seven fields: {line:?}");

    let mut a = f[0].split('/');
    let a_ref = a.next().expect("input reference");
    let a_alts: Vec<&str> = a.collect();
    let a_pos: u64 = f[1].parse().expect("input position");
    assert_eq!(f[2], "1", "the input is always forward-strand");
    let b_pos: u64 = f[4].parse().expect("record position");
    let b_strand: i8 = f[5].parse().expect("record strand");

    let ours = match_alleles(a_ref, &a_alts, a_pos, f[3], b_pos, b_strand)
        .into_iter()
        .map(|m| (a_alts[m.alt_index].to_string(), m.record_allele))
        .collect();
    (ours, parse_expected(f[6]))
}

fn check_file(path: &Path) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    let mut checked = 0;
    let mut wrong = Vec::new();
    for line in text
        .lines()
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
    {
        let (ours, theirs) = answer(line);
        if ours != theirs && wrong.len() < 10 {
            wrong.push(format!(
                "{line}\n    VEP says {theirs:?}\n    we say  {ours:?}"
            ));
        }
        checked += 1;
    }
    assert!(
        wrong.is_empty(),
        "{} of the first mismatches out of {checked} pairs:\n{}",
        wrong.len(),
        wrong.join("\n")
    );
    checked
}

#[test]
fn the_matcher_reproduces_vep_on_the_committed_pairs() {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("fixtures/vep_match/golden.tsv");
    let checked = check_file(&path);
    // The file must keep both kinds of answer, or a matcher that always says
    // "no" (or always "yes") would pass.
    let text = std::fs::read_to_string(&path).unwrap();
    let rows: Vec<&str> = text.lines().filter(|l| !l.starts_with('#')).collect();
    let matches = rows.iter().filter(|l| !l.ends_with("\t-")).count();
    assert!(checked > 1_500, "golden.tsv shrank to {checked} pairs");
    assert!(
        matches > 500 && rows.len() - matches > 500,
        "{matches} matches of {}",
        rows.len()
    );
}

#[test]
#[ignore = "needs FASTVEP_MATCH_CASES, a file made by fixtures/vep_match/regen.sh (see the module docs)"]
fn the_matcher_reproduces_vep_on_a_large_generated_set() {
    let path = std::env::var_os("FASTVEP_MATCH_CASES").expect(
        "set FASTVEP_MATCH_CASES to a file of pairs made by \
         crates/fastvep-cache/fixtures/vep_match/regen.sh big N FILE",
    );
    let checked = check_file(Path::new(&path));
    eprintln!("{checked} pairs agree with VEP");
}
