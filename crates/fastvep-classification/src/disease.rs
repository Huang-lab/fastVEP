//! The classifier's typed reading of a gene's disease annotation.
//!
//! The `omim` gene database carries a list of phenotype strings per gene. Two
//! sources fill it, and they say different things:
//!
//! * **ClinGen Gene-Disease Validity**, converted by `clingen_gdv_to_oga.py`.
//!   Every string ends in a tag the converter writes from the curation's own
//!   columns: `<disease label> (ClinGen <classification>/<MOI>, <MONDO id>)`,
//!   for instance `Bohring-Opitz syndrome (ClinGen Definitive/AD, MONDO:0011510)`.
//!   The classification and the mode of inheritance are *curated fields*, one
//!   value per curation, from a closed list.
//! * **OMIM `genemap2.txt`**, the legacy source. Free text, no tag: the
//!   inheritance, where it is stated at all, is a word in the phenotype.
//!
//! Everything the criteria need from those strings is read here, once, and
//! handed over typed. Before this module the criteria each searched the text
//! for the words "dominant" and "recessive" and never looked at the tag, so the
//! curated inheritance of a ClinGen curation reached the classifier only when
//! the disease's *name* happened to repeat it: 39 of the 740 genes whose
//! established ClinGen curations are all AD, and 60 of the 1,285 that are all
//! AR, among the genes with a variant in the ClinVar 2-star+ benchmark
//! (`docs/ACMG.md`, limitation 20).
//!
//! Parsing borrows from the stored strings and allocates nothing, because the
//! classifier asks these questions several times per variant, transcript and
//! allele.

use crate::sa_extract::OmimData;

/// ClinGen's mode-of-inheritance codes, as the Gene-Disease Validity export
/// spells them. The six codes below are the only ones in the export of
/// 11 October 2026 (3,705 curations; the export is not in the repository, and
/// the 263 curations of `crates/fastvep-cli/fixtures/clingen/` are).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum ModeOfInheritance {
    /// `AD`: autosomal dominant.
    AutosomalDominant,
    /// `AR`: autosomal recessive.
    AutosomalRecessive,
    /// `XL`: X-linked. ClinGen's code does not say dominant or recessive.
    XLinked,
    /// `SD`: semidominant. Heterozygotes are affected, and homozygotes more
    /// severely, so it is dominant and recessive at once: a criterion that needs
    /// a gene to be dominant only or recessive only does not apply to it (PM2
    /// takes the strict threshold, BS2 the arm for homozygotes and hemizygotes, and
    /// BP2's trans arm does not apply), and PM3, which needs a recessive
    /// disease, does.
    Semidominant,
    /// `MT`: mitochondrial.
    Mitochondrial,
    /// `UD`: undetermined. ClinGen says it does not know.
    Undetermined,
}

impl ModeOfInheritance {
    /// The mode a ClinGen code names, `None` for a code this module does not
    /// know. An unknown code is not guessed at: it reads as no inheritance
    /// information, which is what every criterion already does for a gene it
    /// has nothing about.
    pub fn from_code(code: &str) -> Option<Self> {
        match code.trim() {
            "AD" => Some(Self::AutosomalDominant),
            "AR" => Some(Self::AutosomalRecessive),
            "XL" => Some(Self::XLinked),
            "SD" => Some(Self::Semidominant),
            "MT" => Some(Self::Mitochondrial),
            "UD" => Some(Self::Undetermined),
            _ => None,
        }
    }

    /// ClinGen's code for this mode.
    pub fn code(self) -> &'static str {
        match self {
            Self::AutosomalDominant => "AD",
            Self::AutosomalRecessive => "AR",
            Self::XLinked => "XL",
            Self::Semidominant => "SD",
            Self::Mitochondrial => "MT",
            Self::Undetermined => "UD",
        }
    }

    /// This mode's bit in the set `has_several_disorder_contexts` keeps. An
    /// exhaustive match, so that a mode added to the enum has to be given a bit
    /// here and cannot inherit one that is already taken.
    const fn bit(self) -> u8 {
        match self {
            Self::AutosomalDominant => 1 << 0,
            Self::AutosomalRecessive => 1 << 1,
            Self::XLinked => 1 << 2,
            Self::Semidominant => 1 << 3,
            Self::Mitochondrial => 1 << 4,
            Self::Undetermined => 1 << 5,
        }
    }
}

/// The bit of the entries that carry no curated mode: none of theirs is a context of a mode.
const NO_MODE_BIT: u8 = 1 << 6;

/// How a gene, or one of its disorders, is inherited, as far as the criteria
/// that branch on it can tell.
///
/// Two flags and not an enum, because a gene is routinely both and the criteria
/// ask two separate questions: PM2 and BS2 branch on "is it recessive only" and
/// "is it dominant only", and a gene that is neither (X-linked, mitochondrial,
/// undetermined, or not annotated) takes the arm for unknown inheritance, which
/// is the cautious one for each of them.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct Inheritance {
    /// Heterozygotes are affected.
    pub dominant: bool,
    /// Homozygotes and compound heterozygotes are affected.
    pub recessive: bool,
}

impl Inheritance {
    /// Both flags of `other` set here too.
    fn merge(&mut self, other: Inheritance) {
        self.dominant |= other.dominant;
        self.recessive |= other.recessive;
    }

    /// What a curated mode of inheritance says. X-linked, mitochondrial and
    /// undetermined say neither: an X-linked disorder is neither the
    /// heterozygote's nor the homozygote's, and ClinGen's `XL` does not split
    /// dominant from recessive, so a guess here would be the classifier's.
    fn from_mode(mode: ModeOfInheritance) -> Self {
        match mode {
            ModeOfInheritance::AutosomalDominant => Self {
                dominant: true,
                recessive: false,
            },
            ModeOfInheritance::AutosomalRecessive => Self {
                dominant: false,
                recessive: true,
            },
            ModeOfInheritance::Semidominant => Self {
                dominant: true,
                recessive: true,
            },
            ModeOfInheritance::XLinked
            | ModeOfInheritance::Mitochondrial
            | ModeOfInheritance::Undetermined => Self::default(),
        }
    }

    /// What the words of an untagged phenotype say, for the legacy OMIM source,
    /// which has no curated field: an OMIM-style `Autosomal dominant`, a
    /// `{AD}` marker, or "dominant" (or "recessive") in the name of the
    /// disease. A phenotype that names both is read through the explicit forms
    /// only (`Autosomal dominant`, `{AD}`), so "dominant and recessive" in a
    /// name alone ends up as neither.
    fn from_text(text: &str) -> Self {
        let dominant = contains_ignore_case(text, "dominant");
        let recessive = contains_ignore_case(text, "recessive");
        Self {
            dominant: contains_ignore_case(text, "autosomal dominant")
                || contains_ignore_case(text, "{ad}")
                || (dominant && !recessive),
            recessive: contains_ignore_case(text, "autosomal recessive")
                || contains_ignore_case(text, "{ar}")
                || (recessive && !dominant),
        }
    }
}

/// ASCII-case-insensitive substring test. The needles are ASCII, so folding
/// ASCII is enough, and it does not allocate a lower-cased copy of every
/// phenotype on every call the way `to_lowercase().contains(..)` did.
fn contains_ignore_case(haystack: &str, needle: &str) -> bool {
    let (h, n) = (haystack.as_bytes(), needle.as_bytes());
    n.is_empty() || h.windows(n.len()).any(|w| w.eq_ignore_ascii_case(n))
}

/// The ClinGen tag of a phenotype: `(ClinGen <classification>/<MOI>, <MONDO>)`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ClingenTag<'a> {
    /// The classification as ClinGen words it: `Definitive`, `Strong`,
    /// `Moderate`, `Limited`, `Disputed`, `Refuted` or
    /// `No Known Disease Relationship`.
    pub classification: &'a str,
    /// The curated mode of inheritance, `None` where the code is blank or not
    /// one this module knows.
    pub moi: Option<ModeOfInheritance>,
    /// The MONDO identifier of the disease, `None` for the converter's `.`
    /// placeholder.
    pub mondo: Option<&'a str>,
}

impl ClingenTag<'_> {
    /// Whether the classification falls short of the SVI evidence threshold -
    /// `Limited`, `Disputed`, `Refuted` or `No Known Disease Relationship`.
    ///
    /// Only these four are weak. A classification this module has never seen
    /// counts as established, as it always has, because the alternative is to
    /// withhold PVS1, PP2 and PM1 from a gene on the strength of a word nobody
    /// has checked.
    pub fn is_weak(&self) -> bool {
        matches!(
            self.classification,
            "Limited" | "Disputed" | "Refuted" | "No Known Disease Relationship"
        )
    }
}

/// One phenotype string of a gene, read into its parts.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DiseaseEntry<'a> {
    /// The string as the gene database carries it.
    pub raw: &'a str,
    /// The disease label, without ClinGen's tag.
    pub label: &'a str,
    /// ClinGen's classification, inheritance and disease id; `None` for an OMIM
    /// `genemap2` phenotype, which carries none, and for a stray fragment of a
    /// label that was split on `;` (see [`OmimData::entries`]).
    pub clingen: Option<ClingenTag<'a>>,
}

impl<'a> DiseaseEntry<'a> {
    /// Read one phenotype string.
    ///
    /// The tag is the *last* `(ClinGen ` in the string and has to close the
    /// string, so a label that has parentheses of its own, or a comma or a
    /// slash, cannot be mistaken for it. Anything that does not end in a
    /// well-formed tag is an untagged entry and is read as text.
    pub fn parse(raw: &'a str) -> Self {
        let trimmed = raw.trim();
        const OPEN: &str = "(ClinGen ";
        let tagged = trimmed.rfind(OPEN).and_then(|start| {
            let inner = trimmed[start + OPEN.len()..].strip_suffix(')')?;
            let (head, mondo) = match inner.split_once(',') {
                Some((head, mondo)) => (head, Some(mondo.trim())),
                None => (inner, None),
            };
            let (classification, moi) = match head.split_once('/') {
                Some((class, moi)) => (class.trim(), ModeOfInheritance::from_code(moi)),
                None => (head.trim(), None),
            };
            let tag = ClingenTag {
                classification,
                moi,
                mondo: mondo.filter(|m| !m.is_empty() && *m != "."),
            };
            Some((trimmed[..start].trim_end(), tag))
        });
        match tagged {
            Some((label, tag)) => Self {
                raw,
                label,
                clingen: Some(tag),
            },
            None => Self {
                raw,
                label: trimmed,
                clingen: None,
            },
        }
    }

    /// Whether the disease meets the SVI evidence threshold. An untagged entry
    /// does: OMIM listing a phenotype is the assertion, and the legacy source
    /// has no weaker tier to fall back on.
    pub fn is_established(&self) -> bool {
        !self.clingen.is_some_and(|t| t.is_weak())
    }

    /// What the entry says about inheritance: the curated mode for a ClinGen
    /// curation, whatever the label says, and the words of the label for an
    /// OMIM phenotype.
    ///
    /// The curated field wins outright because it is the more reliable of the
    /// two and the two do not disagree where both speak: in the export of
    /// 11 October 2026 no label says dominant on a curation tagged AR or
    /// recessive on one tagged AD. The only mismatch is three X-linked
    /// disorders whose label says "dominant" (GJB1, PDK3 and PHEX, in the
    /// exports of 30 April and of 11 October 2026), which the tag reads,
    /// correctly, as X-linked.
    pub fn inheritance(&self) -> Inheritance {
        match self.clingen {
            Some(tag) => tag.moi.map(Inheritance::from_mode).unwrap_or_default(),
            None => Inheritance::from_text(self.label),
        }
    }
}

/// The established diseases of a gene that share one mode of inheritance.
///
/// A gene with an autosomal dominant disease and an autosomal recessive one -
/// 124 genes in the ClinGen export of 11 October 2026, ATM, ABCC8, ACTA1, ALPL
/// and ANO5 among them - is two questions to the criteria that branch on
/// inheritance, not one: PM2's frequency threshold, which of BS2's arms applies,
/// whether PM3 or BP2's trans arm can apply. ClinGen's guidance on genes with
/// several disorders
/// (January 2024) is to answer each against its own disease, and this is the
/// unit it is answered for.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DisorderContext<'a> {
    /// The mode of inheritance the diseases share; `None` for the untagged
    /// entries of a source that curates none.
    pub moi: Option<ModeOfInheritance>,
    /// The established diseases, in the order the gene lists them.
    pub entries: Vec<DiseaseEntry<'a>>,
}

impl DisorderContext<'_> {
    /// ClinGen's code for the shared mode, or `unspecified` for the entries that
    /// carry no curated mode: a curation whose code ClinGen left blank or this
    /// module does not know, and the untagged entries of the legacy OMIM source.
    /// Those all share one context, so a gene of that source is one context and
    /// has no list of its own.
    pub fn inheritance_code(&self) -> &'static str {
        self.moi.map_or("unspecified", ModeOfInheritance::code)
    }

    /// The gene's disease annotation reduced to this context alone, which is
    /// what every criterion that reads `ClassificationInput::omim` then
    /// evaluates against.
    pub fn restricted_annotation(&self) -> OmimData {
        OmimData {
            mim_number: None,
            phenotypes: Some(self.entries.iter().map(|e| e.raw.to_string()).collect()),
        }
    }
}

impl OmimData {
    /// The gene's phenotype strings, read.
    ///
    /// **A string without a tag is dropped from a gene that has any tagged
    /// string.** The two sources never mix inside one gene: a ClinGen build
    /// tags every phenotype it writes and an OMIM build tags none. An untagged
    /// string beside tagged ones is therefore a piece of a label that held a
    /// `;`. The converter joins phenotypes with `;` and the reader splits them
    /// on it, so `x; y (ClinGen Limited/AD, MONDO:1)` arrives as `x` and
    /// `y (ClinGen Limited/AD, MONDO:1)`, and the piece without the tag counted
    /// as an established disease. NSMCE3, whose only curation is Limited, was
    /// the one weak-only gene that escaped the validity gate that way
    /// (`docs/ACMG.md`, limitation 25). MAPK8IP3 and SCAPER carry the same
    /// split in a Definitive curation, where it does no harm.
    pub fn entries(&self) -> impl Iterator<Item = DiseaseEntry<'_>> + '_ {
        let phenotypes: &[String] = self.phenotypes.as_deref().unwrap_or(&[]);
        let tagged_gene = phenotypes
            .iter()
            .any(|p| DiseaseEntry::parse(p).clingen.is_some());
        phenotypes
            .iter()
            .map(|p| DiseaseEntry::parse(p))
            .filter(move |e| !e.raw.trim().is_empty() && (!tagged_gene || e.clingen.is_some()))
    }

    /// True when at least one listed disease meets the SVI evidence threshold -
    /// a ClinGen classification of Definitive, Strong or Moderate.
    ///
    /// The `.oga` carries every classification, including the weak and negative
    /// ones, because the distinction this method draws is only possible if they
    /// are present. Reading the classification out of the phenotype string is
    /// what lets `has_no_valid_relationship` mean "ClinGen looked and found
    /// little or nothing" rather than "not in the file".
    ///
    /// A gene from an OMIM `genemap2` source carries no ClinGen classification
    /// at all. Those entries count as established: OMIM listing a phenotype is
    /// the assertion, and the legacy source has no weaker tier to fall back on.
    pub fn has_established_relationship(&self) -> bool {
        self.entries().any(|e| e.is_established())
    }

    /// True when ClinGen curated this gene and concluded there is little or no
    /// evidence for any of its proposed disease relationships.
    ///
    /// This is positive evidence against a gene-disease relationship, and it is
    /// what the validity gate blocks on. Absence from the file is *not* this:
    /// ClinGen has reached roughly 2,400 genes, so a gene it has not curated is
    /// overwhelmingly one nobody got to yet. Treating absence as this cost 1,497
    /// pathogenic calls in run v10 - SPAST, ABCB11, FLG and LAMB3 among them.
    pub fn has_no_valid_relationship(&self) -> bool {
        let mut any = false;
        for entry in self.entries() {
            if entry.is_established() {
                return false;
            }
            any = true;
        }
        any
    }

    /// How the gene's established diseases are inherited.
    ///
    /// Considers only diseases that met the evidence threshold: the inheritance
    /// of a relationship ClinGen refuted is not a fact about this gene, and PM2,
    /// BS2, PM3 and BP2 all branch on it. A gene with an autosomal dominant
    /// disease and an autosomal recessive one is both. PM2's threshold, BS2's arm
    /// and BP2's trans arm treat that as neither "dominant only" nor
    /// "recessive only", which is the cautious reading for each; PM3 asks only
    /// that one established disease be recessive, because the observation it
    /// counts is evidence about the variant whatever else the gene causes.
    /// [`Self::disorder_contexts`] is how a caller asks the question one
    /// disease at a time.
    pub fn inheritance(&self) -> Inheritance {
        let mut all = Inheritance::default();
        for entry in self.entries().filter(|e| e.is_established()) {
            all.merge(entry.inheritance());
        }
        all
    }

    /// Check if any established disease is dominantly inherited: ClinGen's
    /// curated `AD` or `SD`, or, for an OMIM phenotype, the words that say so.
    pub fn has_dominant_inheritance(&self) -> bool {
        self.inheritance().dominant
    }

    /// Check if any established disease is recessively inherited: ClinGen's
    /// curated `AR` or `SD`, or, for an OMIM phenotype, the words that say so.
    pub fn has_recessive_inheritance(&self) -> bool {
        self.inheritance().recessive
    }

    /// Whether the established diseases fall under more than one mode of
    /// inheritance, so that [`Self::disorder_contexts`] is worth building.
    ///
    /// Asked of every variant of every gene that has a disease annotation, so it
    /// is answered with a bit set and no allocation.
    pub fn has_several_disorder_contexts(&self) -> bool {
        // One bit per code, and one for every entry without a curated mode, so
        // that two contexts are two bits - the same grouping `disorder_contexts`
        // makes, which this must agree with.
        let mut seen = 0u8;
        for entry in self.entries().filter(|e| e.is_established()) {
            let moi = entry.clingen.and_then(|t| t.moi);
            seen |= moi.map_or(NO_MODE_BIT, ModeOfInheritance::bit);
            if seen.count_ones() > 1 {
                return true;
            }
        }
        false
    }

    /// The established diseases grouped by mode of inheritance, in ClinGen's
    /// code order (AD, AR, XL, SD, MT, UD) with the entries that carry no
    /// curated mode last.
    ///
    /// A curation whose code ClinGen left blank or this module does not know is
    /// grouped with the entries that carry none: it says nothing about
    /// inheritance, and a context of its own would be a context nobody asked
    /// for.
    pub fn disorder_contexts(&self) -> Vec<DisorderContext<'_>> {
        let mut contexts: Vec<DisorderContext<'_>> = Vec::new();
        for entry in self.entries().filter(|e| e.is_established()) {
            let moi = entry.clingen.and_then(|t| t.moi);
            match contexts.iter_mut().find(|c| c.moi == moi) {
                Some(context) => context.entries.push(entry),
                None => contexts.push(DisorderContext {
                    moi,
                    entries: vec![entry],
                }),
            }
        }
        contexts.sort_by_key(|c| c.moi.map_or(u8::MAX, |m| m as u8));
        contexts
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn gene(phenotypes: &[&str]) -> OmimData {
        OmimData {
            mim_number: None,
            phenotypes: Some(phenotypes.iter().map(|s| s.to_string()).collect()),
        }
    }

    // ── The tag ──

    #[test]
    fn test_parse_reads_classification_inheritance_and_disease_id() {
        let e =
            DiseaseEntry::parse("Bohring-Opitz syndrome (ClinGen Definitive/AD, MONDO:0011510)");
        assert_eq!(e.label, "Bohring-Opitz syndrome");
        let tag = e.clingen.expect("a tagged entry");
        assert_eq!(tag.classification, "Definitive");
        assert_eq!(tag.moi, Some(ModeOfInheritance::AutosomalDominant));
        assert_eq!(tag.mondo, Some("MONDO:0011510"));
        assert!(e.is_established());
    }

    #[test]
    fn test_parse_reads_every_code_the_export_uses() {
        for (code, mode) in [
            ("AD", ModeOfInheritance::AutosomalDominant),
            ("AR", ModeOfInheritance::AutosomalRecessive),
            ("XL", ModeOfInheritance::XLinked),
            ("SD", ModeOfInheritance::Semidominant),
            ("MT", ModeOfInheritance::Mitochondrial),
            ("UD", ModeOfInheritance::Undetermined),
        ] {
            let raw = format!("a disease (ClinGen Strong/{code}, MONDO:0000001)");
            let e = DiseaseEntry::parse(&raw);
            assert_eq!(e.clingen.unwrap().moi, Some(mode), "{code}");
            assert_eq!(mode.code(), code);
        }
    }

    #[test]
    fn test_parse_label_with_its_own_punctuation_is_not_mistaken_for_the_tag() {
        // Real labels carry commas, slashes and parentheses.
        let e = DiseaseEntry::parse(
            "Charcot-Marie-Tooth disease, type 2N (axonal, AD/AR variants) (ClinGen Moderate/AR, MONDO:0013212)",
        );
        assert_eq!(
            e.label,
            "Charcot-Marie-Tooth disease, type 2N (axonal, AD/AR variants)"
        );
        let tag = e.clingen.unwrap();
        assert_eq!(tag.classification, "Moderate");
        assert_eq!(tag.moi, Some(ModeOfInheritance::AutosomalRecessive));
    }

    #[test]
    fn test_parse_multiword_classification() {
        let e = DiseaseEntry::parse(
            "congenital heart disease (ClinGen No Known Disease Relationship/UD, MONDO:0005453)",
        );
        let tag = e.clingen.unwrap();
        assert_eq!(tag.classification, "No Known Disease Relationship");
        assert!(tag.is_weak());
        assert!(!e.is_established());
    }

    #[test]
    fn test_parse_placeholder_disease_id_is_none() {
        let e = DiseaseEntry::parse("a disease (ClinGen Definitive/AD, .)");
        assert_eq!(e.clingen.unwrap().mondo, None);
    }

    #[test]
    fn test_parse_an_unknown_code_reads_as_no_inheritance_and_stays_tagged() {
        let e = DiseaseEntry::parse("a disease (ClinGen Definitive/ZZ, MONDO:0000001)");
        let tag = e.clingen.unwrap();
        assert_eq!(tag.moi, None);
        assert_eq!(e.inheritance(), Inheritance::default());
    }

    #[test]
    fn test_parse_text_after_the_tag_leaves_the_entry_untagged() {
        // A tag has to close the string. Text after it means this is not the
        // converter's output, so it is read as the free text it is.
        let e = DiseaseEntry::parse("x (ClinGen Definitive/AD, MONDO:1) and more");
        assert!(e.clingen.is_none());
    }

    #[test]
    fn test_unrecognised_classification_counts_as_established() {
        // The behaviour before the typed reading, kept: only the four named
        // classifications are weak.
        let e = DiseaseEntry::parse("a disease (ClinGen Supportive/AD, MONDO:1)");
        assert!(e.is_established());
    }

    // ── What the tag says about inheritance ──

    #[test]
    fn test_curated_inheritance_does_not_need_the_disease_name_to_repeat_it() {
        // The defect this module exists to remove. Neither label says a word
        // about inheritance, and the old reading found nothing in either.
        let ad = gene(&["Bohring-Opitz syndrome (ClinGen Definitive/AD, MONDO:0011510)"]);
        assert!(ad.has_dominant_inheritance());
        assert!(!ad.has_recessive_inheritance());

        let ar = gene(&["ataxia-telangiectasia (ClinGen Definitive/AR, MONDO:0008840)"]);
        assert!(ar.has_recessive_inheritance());
        assert!(!ar.has_dominant_inheritance());
    }

    #[test]
    fn test_semidominant_is_both() {
        let sd = gene(&["a disease (ClinGen Definitive/SD, MONDO:0000001)"]);
        assert_eq!(
            sd.inheritance(),
            Inheritance {
                dominant: true,
                recessive: true
            }
        );
    }

    #[test]
    fn test_x_linked_mitochondrial_and_undetermined_are_neither() {
        for code in ["XL", "MT", "UD"] {
            let g = gene(&[&format!(
                "a disease (ClinGen Definitive/{code}, MONDO:0000001)"
            )]);
            assert_eq!(g.inheritance(), Inheritance::default(), "{code}");
        }
    }

    #[test]
    fn test_the_tag_beats_the_label() {
        // Three real X-linked curations have "dominant" in the label. The tag
        // is the curated field, so the gene is X-linked and not dominant.
        let g =
            gene(&["X-linked dominant chondrodysplasia (ClinGen Definitive/XL, MONDO:0000001)"]);
        assert!(!g.has_dominant_inheritance());
        assert!(!g.has_recessive_inheritance());
    }

    #[test]
    fn test_weak_curations_do_not_state_the_genes_inheritance() {
        // A relationship ClinGen refuted is not a fact about this gene.
        let g = gene(&[
            "a refuted disease (ClinGen Refuted/AD, MONDO:0000001)",
            "a real disease (ClinGen Definitive/AR, MONDO:0000002)",
        ]);
        assert!(!g.has_dominant_inheritance());
        assert!(g.has_recessive_inheritance());
    }

    #[test]
    fn test_a_gene_with_a_dominant_and_a_recessive_disease_is_both() {
        let g = gene(&[
            "breast cancer (ClinGen Definitive/AD, MONDO:0000001)",
            "ataxia-telangiectasia (ClinGen Definitive/AR, MONDO:0000002)",
        ]);
        assert!(g.has_dominant_inheritance());
        assert!(g.has_recessive_inheritance());
    }

    // ── The legacy OMIM source keeps reading words ──

    #[test]
    fn test_omim_phenotypes_are_read_as_text_as_before() {
        for (text, dominant, recessive) in [
            ("Breast cancer, 114480 (3), Autosomal dominant", true, false),
            (
                "Cystic fibrosis, 219700 (3), Autosomal recessive",
                false,
                true,
            ),
            ("a disease {AD}", true, false),
            ("a disease {ar}", false, true),
            ("dominant disorder", true, false),
            ("recessive disorder", false, true),
            ("X-linked recessive disorder", false, true),
            // "dominant" and "recessive" both present, neither in the explicit form
            ("dominant and recessive forms", false, false),
            // the explicit form wins over the other word
            ("Autosomal dominant, with recessive forms", true, false),
            ("a disease with no inheritance stated", false, false),
        ] {
            let g = gene(&[text]);
            assert_eq!(
                (g.has_dominant_inheritance(), g.has_recessive_inheritance()),
                (dominant, recessive),
                "{text}"
            );
        }
    }

    #[test]
    fn test_text_matching_ignores_ascii_case() {
        let g = gene(&["AUTOSOMAL DOMINANT"]);
        assert!(g.has_dominant_inheritance());
    }

    // ── Validity ──

    #[test]
    fn test_established_and_no_valid_relationship() {
        let weak = gene(&[
            "x (ClinGen Disputed/AD, MONDO:1)",
            "y (ClinGen Limited/AR, MONDO:2)",
        ]);
        assert!(!weak.has_established_relationship());
        assert!(weak.has_no_valid_relationship());

        let mixed = gene(&[
            "x (ClinGen Disputed/AD, MONDO:1)",
            "y (ClinGen Moderate/AR, MONDO:2)",
        ]);
        assert!(mixed.has_established_relationship());
        assert!(!mixed.has_no_valid_relationship());

        let omim = gene(&["Breast cancer, 114480 (3), Autosomal dominant"]);
        assert!(omim.has_established_relationship());
        assert!(!omim.has_no_valid_relationship());

        for none in [
            OmimData {
                mim_number: Some(0),
                phenotypes: None,
            },
            gene(&[]),
            gene(&["", "  "]),
        ] {
            assert!(!none.has_established_relationship());
            assert!(!none.has_no_valid_relationship(), "{none:?}");
        }
    }

    #[test]
    fn test_a_fragment_of_a_split_label_is_not_an_established_disease() {
        // NSMCE3's only curation is Limited, and its label ends in ';'. The
        // converter joins phenotypes with ';' and the reader splits on it, so
        // the string reaches the classifier as an untagged piece and a tagged
        // one. The untagged piece used to count as an established disease and
        // the gene escaped the validity gate.
        let nsmce3 = gene(&[
            "immunodeficiency 49 and neutropenia;",
            "NSMCE3-related (ClinGen Limited/AR, MONDO:0044118)",
        ]);
        assert!(!nsmce3.has_established_relationship());
        assert!(nsmce3.has_no_valid_relationship());

        // The same split in a Definitive curation changes nothing it should not.
        let split = gene(&[
            "neurodevelopmental disorder with",
            "a variable phenotype (ClinGen Definitive/AD, MONDO:0000001)",
        ]);
        assert!(split.has_established_relationship());
        assert!(split.has_dominant_inheritance());
        assert_eq!(split.entries().count(), 1);
    }

    // ── Contexts ──

    #[test]
    fn test_one_inheritance_is_not_several_contexts() {
        let g = gene(&[
            "a (ClinGen Definitive/AD, MONDO:1)",
            "b (ClinGen Strong/AD, MONDO:2)",
            "c (ClinGen Disputed/AR, MONDO:3)",
        ]);
        assert!(!g.has_several_disorder_contexts());
        let contexts = g.disorder_contexts();
        assert_eq!(contexts.len(), 1);
        assert_eq!(contexts[0].inheritance_code(), "AD");
        assert_eq!(
            contexts[0].entries.len(),
            2,
            "the disputed AR one is not established"
        );
    }

    #[test]
    fn test_a_dominant_and_a_recessive_disease_are_two_contexts() {
        let g = gene(&[
            "ataxia-telangiectasia (ClinGen Definitive/AR, MONDO:1)",
            "breast cancer (ClinGen Definitive/AD, MONDO:2)",
            "another recessive (ClinGen Moderate/AR, MONDO:3)",
        ]);
        assert!(g.has_several_disorder_contexts());
        let contexts = g.disorder_contexts();
        let codes: Vec<_> = contexts.iter().map(|c| c.inheritance_code()).collect();
        assert_eq!(
            codes,
            ["AD", "AR"],
            "in ClinGen's code order, not the gene's"
        );
        assert_eq!(contexts[1].entries.len(), 2);

        // Each context, restricted, reads as that inheritance alone.
        let ad = contexts[0].restricted_annotation();
        assert!(ad.has_dominant_inheritance() && !ad.has_recessive_inheritance());
        let ar = contexts[1].restricted_annotation();
        assert!(ar.has_recessive_inheritance() && !ar.has_dominant_inheritance());
    }

    #[test]
    fn test_contexts_ignore_weak_curations_and_stray_fragments() {
        let g = gene(&[
            "a stray piece",
            "a (ClinGen Definitive/AD, MONDO:1)",
            "b (ClinGen Refuted/AR, MONDO:2)",
        ]);
        assert!(!g.has_several_disorder_contexts());
        assert_eq!(g.disorder_contexts().len(), 1);
    }

    #[test]
    fn test_untagged_entries_form_their_own_context() {
        let g = gene(&[
            "Breast cancer, 114480 (3), Autosomal dominant",
            "Ovarian cancer, 167000 (3)",
        ]);
        assert!(!g.has_several_disorder_contexts());
        let contexts = g.disorder_contexts();
        assert_eq!(contexts.len(), 1);
        assert_eq!(contexts[0].inheritance_code(), "unspecified");
        assert_eq!(contexts[0].entries.len(), 2);
    }

    #[test]
    fn test_a_gene_with_no_annotation_has_no_contexts() {
        assert!(gene(&[]).disorder_contexts().is_empty());
        assert!(!gene(&[]).has_several_disorder_contexts());
    }
    /// `has_several_disorder_contexts` keeps one bit for each mode and one for the entries that carry none: two
    /// modes that shared a bit would be one context to it and two to `disorder_contexts`.
    #[test]
    fn test_every_mode_has_a_bit_of_its_own() {
        let mut seen = NO_MODE_BIT;
        for code in ["AD", "AR", "XL", "SD", "MT", "UD"] {
            let mode = ModeOfInheritance::from_code(code).unwrap();
            assert_eq!(seen & mode.bit(), 0, "{code} shares a bit");
            seen |= mode.bit();
        }
        assert_eq!(seen.count_ones(), 7);
    }
}
