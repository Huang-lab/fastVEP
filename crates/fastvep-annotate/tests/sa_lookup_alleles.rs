//! What the supplementary-annotation lookup asks each database, and in what form.
//!
//! Every database of alleles (`.osa`, `.osa2`, the ClinVar indices) is keyed on the REF and ALT the
//! VCF wrote, anchor base included, at the position the VCF gave. The parser trims an indel the way
//! VEP does for the consequence step: an insertion `A>AT` at 100 becomes `->T` at 101. The
//! annotation loop of the library asked with the trimmed form and the one of the command line with
//! the VCF form, so the command line found an indel in gnomAD and ClinVar and both web servers
//! missed it, which reads as "absent from gnomAD" and lets PM2 fire on a common variant. The loops
//! now share `supplementary_for_alleles`, and these tests hold what it asks.

use anyhow::Result;
use fastvep_annotate::{annotate_intergenic, sa_duplicate_answers, supplementary_for_alleles};
use fastvep_cache::annotation::{AnnotationProvider, AnnotationValue, SaMetadata};
use fastvep_io::variant::VariationFeature;
use fastvep_io::vcf::VcfParser;
use std::sync::Mutex;

/// A database that holds the records it was given, keyed as the real ones are, and remembers every
/// question it was asked.
struct Keyed {
    metadata: SaMetadata,
    holds: Vec<(&'static str, u64, &'static str, &'static str)>,
    asked: Mutex<Vec<(u64, String, String)>>,
}

impl Keyed {
    fn new(
        json_key: &str,
        match_by_allele: bool,
        holds: Vec<(&'static str, u64, &'static str, &'static str)>,
    ) -> Self {
        Self {
            metadata: SaMetadata {
                name: json_key.into(),
                version: "1".into(),
                description: "test".into(),
                assembly: "GRCh38".into(),
                json_key: json_key.into(),
                match_by_allele,
                is_array: false,
                is_positional: !match_by_allele,
            },
            holds,
            asked: Mutex::new(Vec::new()),
        }
    }
}

impl AnnotationProvider for Keyed {
    fn name(&self) -> &str {
        &self.metadata.name
    }
    fn json_key(&self) -> &str {
        &self.metadata.json_key
    }
    fn metadata(&self) -> &SaMetadata {
        &self.metadata
    }
    fn annotate_position(
        &self,
        chrom: &str,
        pos: u64,
        ref_allele: &str,
        alt_allele: &str,
    ) -> Result<Option<AnnotationValue>> {
        self.asked
            .lock()
            .unwrap()
            .push((pos, ref_allele.to_string(), alt_allele.to_string()));
        let found = if self.metadata.match_by_allele {
            self.holds.iter().any(|(c, p, r, a)| {
                *c == chrom && *p == pos && *r == ref_allele && *a == alt_allele
            })
        } else {
            self.holds
                .iter()
                .any(|(c, p, _, _)| *c == chrom && *p == pos)
        };
        Ok(found
            .then(|| AnnotationValue::Json(format!(r#"{{"from":"{}"}}"#, self.metadata.json_key))))
    }
}

/// One record of a VCF, parsed as the annotation loops parse it and given the placeholder rows
/// that stand for a variant with no overlapping transcript, which is all the lookup needs.
fn variant(chrom: &str, pos: u64, reference: &str, alt: &str) -> VariationFeature {
    let vcf = format!(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n{chrom}\t{pos}\t.\t{reference}\t{alt}\t.\t.\t.\n"
    );
    let mut vfs = VcfParser::new(vcf.as_bytes()).unwrap().read_all().unwrap();
    let mut vf = vfs.remove(0);
    annotate_intergenic(&mut vf);
    vf
}

fn ask(
    vf: &VariationFeature,
    providers: Vec<Box<dyn AnnotationProvider>>,
) -> std::collections::HashMap<String, Vec<(String, String)>> {
    supplementary_for_alleles(vf, &providers, None, false)
}

#[test]
fn an_insertion_is_looked_up_as_the_vcf_wrote_it() {
    let vf = variant("1", 100, "A", "AT");
    // What makes this a test of the form: the parser has trimmed the allele the databases hold.
    assert_eq!(vf.position.start, 101);
    assert_eq!(vf.ref_allele.to_string(), "-");

    let gnomad = Keyed::new("gnomad", true, vec![("1", 100, "A", "AT")]);
    let results = ask(&vf, vec![Box::new(gnomad)]);
    let found = &results["T"];
    assert_eq!(
        found.len(),
        1,
        "gnomAD holds `1:100 A>AT` and was not asked for it: {results:?}"
    );
    assert_eq!(found[0].0, "gnomad");
}

#[test]
fn a_deletion_is_looked_up_as_the_vcf_wrote_it() {
    let vf = variant("1", 200, "CAG", "C");
    assert_eq!(vf.position.start, 201);
    let clinvar = Keyed::new("clinvar", true, vec![("1", 200, "CAG", "C")]);
    let results = ask(&vf, vec![Box::new(clinvar)]);
    assert_eq!(
        results.values().map(Vec::len).sum::<usize>(),
        1,
        "ClinVar holds `1:200 CAG>C` and was not asked for it: {results:?}"
    );
}

#[test]
fn a_substitution_is_asked_as_it_always_was() {
    let vf = variant("1", 300, "A", "G");
    let db = Keyed::new("gnomad", true, vec![("1", 300, "A", "G")]);
    let results = ask(&vf, vec![Box::new(db)]);
    assert_eq!(results["G"].len(), 1);
}

/// A site with a substitution and an insertion is not minimised by the parser as a single-ALT record
/// is, and each of its alleles is asked as the VCF wrote it.
#[test]
fn each_allele_of_a_multiallelic_site_is_looked_up_as_the_vcf_wrote_it() {
    let vf = variant("1", 400, "A", "G,AT");
    let db = Keyed::new(
        "gnomad",
        true,
        vec![("1", 400, "A", "G"), ("1", 400, "A", "AT")],
    );
    let results = ask(&vf, vec![Box::new(db)]);
    assert_eq!(results.len(), 2, "{results:?}");
    for (allele, found) in &results {
        assert_eq!(found.len(), 1, "allele {allele} found nothing: {results:?}");
    }
}

/// A source that is not keyed by allele is asked by position alone, as the command line has always asked it.
#[test]
fn a_source_not_keyed_by_allele_is_asked_by_position_alone() {
    let vf = variant("1", 100, "A", "AT");
    let phylop = std::sync::Arc::new(Keyed::new("phylop", false, vec![("1", 101, "", "")]));
    struct Shared(std::sync::Arc<Keyed>);
    impl AnnotationProvider for Shared {
        fn name(&self) -> &str {
            self.0.name()
        }
        fn json_key(&self) -> &str {
            self.0.json_key()
        }
        fn metadata(&self) -> &SaMetadata {
            self.0.metadata()
        }
        fn annotate_position(
            &self,
            c: &str,
            p: u64,
            r: &str,
            a: &str,
        ) -> Result<Option<AnnotationValue>> {
            self.0.annotate_position(c, p, r, a)
        }
    }
    let results = ask(&vf, vec![Box::new(Shared(phylop.clone()))]);
    assert_eq!(results["T"].len(), 1, "{results:?}");
    let asked = phylop.asked.lock().unwrap();
    assert_eq!(asked.as_slice(), [(101, String::new(), String::new())]);
}

/// Two databases of one source that both hold an allele are recorded at the position the VCF gave,
/// which is the position the databases were asked at, and not at the trimmed one.
#[test]
fn a_source_answered_twice_is_recorded_at_the_position_it_was_asked_at() {
    let vf = variant("1", 100, "A", "AT");
    let first = Keyed::new("twice_for_the_test", true, vec![("1", 100, "A", "AT")]);
    let second = Keyed::new("twice_for_the_test", true, vec![("1", 100, "A", "AT")]);
    ask(&vf, vec![Box::new(first), Box::new(second)]);
    let lines = sa_duplicate_answers().report_lines();
    let line = lines
        .iter()
        .find(|l| l.contains("`twice_for_the_test`"))
        .unwrap_or_else(|| panic!("the duplicate was not reported: {lines:?}"));
    assert!(line.contains("at 1:100)"), "{line}");
}
