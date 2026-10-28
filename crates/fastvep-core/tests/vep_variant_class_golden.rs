//! `VARIANT_CLASS` against Ensembl VEP's own answers.
//!
//! `fixtures/vep_class/golden.tsv` is the class VEP's `SO_variation_class`
//! (ensembl-variation release/115, run as it is) gives each allele string: every
//! handwritten corner of the Perl, then a seeded sample with at least 150 of
//! each class VEP answers with. Our port has to give the same class for every
//! one.
//!
//! To run the same check on a much larger set, make one and point the second test
//! at it; without the variable that test does nothing, so CI stays hermetic:
//!
//! ```text
//! crates/fastvep-core/fixtures/vep_class/regen.sh big 100000 /tmp/classes.tsv
//! FASTVEP_CLASS_CASES=/tmp/classes.tsv cargo test -p fastvep-core --test vep_variant_class_golden
//! ```

use fastvep_core::so_variation_class;
use std::collections::BTreeSet;
use std::path::Path;

fn check_file(path: &Path) -> (usize, BTreeSet<String>) {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    let mut checked = 0;
    let mut classes = BTreeSet::new();
    let mut wrong = Vec::new();
    for line in text.lines().filter(|l| !l.starts_with('#')) {
        // An allele string can be empty, so only the tab splits a line.
        let (alleles, expected) = line
            .split_once('\t')
            .unwrap_or_else(|| panic!("no tab in {line:?}"));
        let ours = so_variation_class(alleles);
        if ours != expected && wrong.len() < 10 {
            wrong.push(format!("{alleles:?}: VEP says {expected}, we say {ours}"));
        }
        classes.insert(expected.to_string());
        checked += 1;
    }
    assert!(
        wrong.is_empty(),
        "mismatches among {checked} strings:\n{}",
        wrong.join("\n")
    );
    (checked, classes)
}

#[test]
fn the_port_reproduces_vep_on_the_committed_strings() {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("fixtures/vep_class/golden.tsv");
    let (checked, classes) = check_file(&path);
    assert!(checked > 1_000, "golden.tsv shrank to {checked} strings");
    // Every class VEP can answer with is in the file, or a port that never says
    // one of them would pass.
    let expected: BTreeSet<String> = [
        "SNV",
        "substitution",
        "insertion",
        "deletion",
        "indel",
        "tandem_repeat",
        "sequence_alteration",
    ]
    .iter()
    .map(ToString::to_string)
    .collect();
    assert_eq!(classes, expected);
}

#[test]
fn the_port_reproduces_vep_on_a_large_generated_set() {
    let Some(path) = std::env::var_os("FASTVEP_CLASS_CASES") else {
        eprintln!("FASTVEP_CLASS_CASES not set; skipping (see the module docs)");
        return;
    };
    let (checked, _) = check_file(Path::new(&path));
    eprintln!("{checked} strings agree with VEP");
}
