//! `VARIANT_CLASS`: the Sequence Ontology term VEP reports for a variant's
//! alleles under `--variant_class`.
//!
//! VEP gets it from `VariationFeature::class_SO_term`, which is
//! `SO_variation_class($allele_string)` of
//! `Bio::EnsEMBL::Variation::Utils::Sequence` (ensembl-variation release/115):
//! a function of the allele string alone, so the class of a multi-allelic site
//! depends on all of its alleles and not on any one of them. This is that
//! function, branch for branch, checked against VEP's own Perl in
//! `tests/vep_variant_class_golden.rs`.
//!
//! It is the class of the allele string *as VEP holds it*: after the VCF
//! parser has stripped the shared anchor base, so `-/T` and not `C/CT`.
//! Structural variants take their class from the parser in VEP and are not
//! covered here.

const SNV: &str = "SNV";
const SUBSTITUTION: &str = "substitution";
const INSERTION: &str = "insertion";
const DELETION: &str = "deletion";
const INDEL: &str = "indel";
const TANDEM_REPEAT: &str = "tandem_repeat";
const SEQUENCE_ALTERATION: &str = "sequence_alteration";

/// `/^[A-Z]+$/`: a non-empty run of capital letters. The Perl's `[A-Z]` is
/// every capital letter, so `N`, `R` and `Y` count as bases.
fn is_run(s: &str) -> bool {
    !s.is_empty() && s.bytes().all(|b| b.is_ascii_uppercase())
}

/// `/^[A-Z](\/[A-Z])+$/`: single letters, two or more, joined by `/`.
fn is_list_of_single_letters(s: &str) -> bool {
    let b = s.as_bytes();
    b.len() >= 3
        && b.len() % 2 == 1
        && b.iter().enumerate().all(|(i, &c)| {
            if i % 2 == 0 {
                c.is_ascii_uppercase()
            } else {
                c == b'/'
            }
        })
}

/// `/^[A-Z]+(\/[A-Z]+)+$/`: runs of capital letters, two or more, joined by `/`.
fn is_list_of_runs(s: &str) -> bool {
    let mut parts = s.split('/');
    let first_ok = parts.next().is_some_and(is_run);
    let mut rest = 0;
    for p in parts {
        if !is_run(p) {
            return false;
        }
        rest += 1;
    }
    first_ok && rest >= 1
}

/// `/\)\d+/`: a closing parenthesis followed by a digit, as in `(CAG)8`.
fn has_repeat_count(s: &str) -> bool {
    s.as_bytes()
        .windows(2)
        .any(|w| w[0] == b')' && w[1].is_ascii_digit())
}

/// The SO term for an Ensembl allele string (`REF/ALT1/ALT2`, `-` for an empty
/// allele), as `SO_variation_class` returns it with its default of the first
/// allele being the reference.
pub fn so_variation_class(alleles: &str) -> &'static str {
    // Anything not recognised stays this.
    let mut class = SEQUENCE_ALTERATION;

    if is_list_of_single_letters(alleles) {
        // A/T, A/T/G
        return SNV;
    }
    if has_repeat_count(alleles) {
        // (CAG)8/(CAG)9
        return TANDEM_REPEAT;
    }

    // Perl's `split` drops trailing empty fields and keeps leading ones.
    let mut parts: Vec<&str> = alleles.split('/').collect();
    while parts.last() == Some(&"") {
        parts.pop();
    }

    match parts.as_slice() {
        [] => {}
        [only] => {
            if only.contains("INS") {
                // (LARGEINSERTION)
                class = INSERTION;
            } else if only.contains("DEL") {
                // (308 BP DELETION)
                class = DELETION;
            }
        }
        [reference, alts @ ..] => {
            if *reference == "-" {
                if alts.len() == 1 && alts[0].contains("DEL") {
                    // -/(LARGEDELETION) (rather oddly!)
                    class = DELETION;
                }
                // A separate test in the Perl, not an else: it can override the
                // one above, though for these inputs it never disagrees.
                if alts.iter().all(|a| is_run(a) || a.contains("INS")) {
                    // -/ATT, -/(LARGEINSERTION)
                    class = INSERTION;
                }
                // else a mix of insertion and deletion: sequence_alteration
            } else if is_run(reference) {
                if alts.iter().all(|a| a.contains('-') || a.contains("DEL")) {
                    // A/-, A/(LARGEDELETION)
                    class = DELETION;
                } else if is_list_of_runs(alleles) {
                    // AA/TT is a substitution, AA/TTTT an indel
                    let same_size = alts.iter().all(|a| a.len() == reference.len());
                    class = if same_size { SUBSTITUTION } else { INDEL };
                }
            } else if reference.contains("DEL") && alts.iter().all(|a| a.contains('-')) {
                // (LARGEDELETION)/-, (2345 BP DELETION)/-
                class = DELETION;
            }
        }
    }
    class
}

#[cfg(test)]
mod tests {
    use super::so_variation_class as class;

    #[test]
    fn the_common_shapes() {
        assert_eq!(class("C/A"), "SNV");
        assert_eq!(class("C/A/G"), "SNV");
        assert_eq!(class("AA/TT"), "substitution");
        assert_eq!(class("ACCCT/-"), "deletion");
        assert_eq!(class("-/T"), "insertion");
        assert_eq!(class("-/ACGT"), "insertion");
        assert_eq!(class("AA/TTTT"), "indel");
        assert_eq!(class("(CAG)8/(CAG)9"), "tandem_repeat");
    }

    /// A multi-allelic site has one class, from all its alleles together.
    #[test]
    fn a_multiallelic_site_has_one_class() {
        // a base, a longer allele: an indel
        assert_eq!(class("A/G/ACG"), "indel");
        // a deletion and a base: neither pure
        assert_eq!(class("A/-/T"), "sequence_alteration");
        assert_eq!(class("-/A/T"), "insertion");
    }

    /// The Perl only calls capital letters bases.
    #[test]
    fn lowercase_is_not_recognised() {
        assert_eq!(class("a/g"), "sequence_alteration");
        assert_eq!(class("A/g"), "sequence_alteration");
    }

    #[test]
    fn malformed_strings_fall_through() {
        for s in ["A", "-", "-/-", "A/", "/A", "A//G", ""] {
            assert_eq!(class(s), "sequence_alteration", "{s:?}");
        }
    }
}
