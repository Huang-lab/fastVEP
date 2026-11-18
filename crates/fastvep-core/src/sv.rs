//! Structural variants as Ensembl VEP holds them.
//!
//! VEP does not predict the effect of a structural variant from its alleles. It reads one fact, the
//! class the parser gave the record (`class_SO_term` of `StructuralVariationFeature`), and every
//! consequence rule of `Utils/VariationEffect.pm` is a test of that class and of where the variant's
//! interval lies against a transcript. [`SvClass`] is that class, and [`StructuralVariant`] is the
//! rest of what the rules and the output read: the alleles VEP writes a row for, each with the name
//! the `Allele` column gives it.
//!
//! The rules themselves are in `fastvep-consequence` (`sv_predictor`) and the reading of a VCF
//! record into this form is in `fastvep-io` (`sv_record`).

use std::num::NonZeroU32;
use std::sync::Arc;

/// The class VEP gives a structural variant: the Sequence Ontology term of `%SO_TERMS` in
/// `Bio::EnsEMBL::Variation::Utils::Config`.
///
/// The set is closed. A record whose type is not one of these is not annotated at all (VEP writes
/// "X is not a supported structural variant type" and goes on to the next line), so no rule has a
/// class it does not know. That is why the rules below are a `match` and not a pattern on the term's
/// text, which is what the Perl does (`class_SO_term =~ /deletion|loss/i`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum SvClass {
    Insertion,
    MobileElementInsertion,
    AluInsertion,
    HervInsertion,
    Line1Insertion,
    SvaInsertion,
    Deletion,
    MobileElementDeletion,
    AluDeletion,
    HervDeletion,
    Line1Deletion,
    SvaDeletion,
    TandemRepeat,
    TandemDuplication,
    Duplication,
    CopyNumberVariation,
    Inversion,
    ChromosomeBreakpoint,
}

impl SvClass {
    /// The term as VEP writes it: `VARIANT_CLASS`, and the `Allele` column for every class but a
    /// breakend with mates.
    pub const fn so_term(self) -> &'static str {
        match self {
            Self::Insertion => "insertion",
            Self::MobileElementInsertion => "mobile_element_insertion",
            Self::AluInsertion => "Alu_insertion",
            Self::HervInsertion => "HERV_insertion",
            Self::Line1Insertion => "LINE1_insertion",
            Self::SvaInsertion => "SVA_insertion",
            Self::Deletion => "deletion",
            Self::MobileElementDeletion => "mobile_element_deletion",
            Self::AluDeletion => "Alu_deletion",
            Self::HervDeletion => "HERV_deletion",
            Self::Line1Deletion => "LINE1_deletion",
            Self::SvaDeletion => "SVA_deletion",
            Self::TandemRepeat => "tandem_repeat",
            Self::TandemDuplication => "tandem_duplication",
            Self::Duplication => "duplication",
            Self::CopyNumberVariation => "copy_number_variation",
            Self::Inversion => "inversion",
            Self::ChromosomeBreakpoint => "chromosome_breakpoint",
        }
    }

    /// `deletion` of `VariationEffect.pm`: a copy number loss, or a class named deletion. It is also
    /// the `deletion` flag of `_bvf_preds` (`/deletion|loss/i`), which is the same set: no class
    /// has `loss` in its name.
    pub const fn is_deletion(self) -> bool {
        matches!(
            self,
            Self::Deletion
                | Self::MobileElementDeletion
                | Self::AluDeletion
                | Self::HervDeletion
                | Self::Line1Deletion
                | Self::SvaDeletion
        )
    }

    /// `copy_number_gain`: a duplication, tandem or not, or a tandem repeat (`/gain/i` names no
    /// class).
    pub const fn is_copy_number_gain(self) -> bool {
        matches!(
            self,
            Self::Duplication | Self::TandemDuplication | Self::TandemRepeat
        )
    }

    /// `insertion`: a copy number gain, or a class named insertion. It is also the
    /// `increase_length` flag of `_bvf_preds` (`/insertion|duplication|gain|tandem/i`), the same
    /// set.
    pub const fn is_insertion(self) -> bool {
        self.is_copy_number_gain()
            || matches!(
                self,
                Self::Insertion
                    | Self::MobileElementInsertion
                    | Self::AluInsertion
                    | Self::HervInsertion
                    | Self::Line1Insertion
                    | Self::SvaInsertion
            )
    }

    /// `chromosome_breakpoint`: the class of a breakend.
    pub const fn is_breakpoint(self) -> bool {
        matches!(self, Self::ChromosomeBreakpoint)
    }

    /// The `decrease_length` flag of `_bvf_preds`, which `feature_truncation` is gated on: a deletion,
    /// or a breakpoint.
    pub const fn decreases_length(self) -> bool {
        self.is_deletion() || self.is_breakpoint()
    }

    /// Ensembl's `get_SO_term` (`modules/Bio/EnsEMBL/VEP/Parser.pm`): the class of the type a record
    /// gives, which is its `ALT` (the alleles joined with `/`) or, where the `ALT` does not start
    /// with a known type, its `SVTYPE`. `Err` carries the abbreviation VEP names in its warning
    /// (`X is not a supported structural variant type`).
    ///
    /// Its steps are tried in this order and the first to match decides:
    ///
    /// 1. a type repeated (`<INS>,<INS>`, `INS/INS`) is the type once;
    /// 2. `INS:ME` and `DEL:ME`, with `:ALU`, `:HERV`, `:LINE1` (or `:L1`) or `:SVA` after them;
    /// 3. `DEL` together with `DUP`: a copy number variation;
    /// 4. `DUP:TANDEM`: a tandem duplication;
    /// 5. `CNV:TR`: a tandem repeat;
    /// 6. `CN` and a digit (`<CN0>`): a deletion for `CN0`, a duplication for `CN2`, a copy number
    ///    variation for any other number;
    /// 7. `CNV`: a copy number variation;
    /// 8. a bracket anywhere, or a `.` at either end: a breakpoint;
    /// 9. `<` at the start or `>` at the end: the name between them, up to its first `:`;
    /// 10. the type as it is.
    ///
    /// The abbreviation is then looked up in the table of [`Self::from_abbreviation`], which is
    /// case-sensitive, so `<del>` is not a deletion.
    pub fn from_type(type_text: &str) -> Result<Self, String> {
        let mut text = type_text;
        // 1. `$type = $parts[0]` when every part, stripped of spaces and one pair of angle
        // brackets and upper-cased, is the first one. Perl's `split` drops trailing empty fields.
        if text.contains(['/', ',']) {
            let mut parts: Vec<&str> = text.split(['/', ',']).collect();
            while parts.last().is_some_and(|p| p.is_empty()) {
                parts.pop();
            }
            if let Some(first) = parts.first() {
                let key = |part: &str| -> String {
                    let part = part.trim();
                    let part = part.strip_prefix('<').unwrap_or(part);
                    let part = part.strip_suffix('>').unwrap_or(part);
                    part.to_ascii_uppercase()
                };
                let first_key = key(first);
                if parts.iter().all(|p| key(p) == first_key) {
                    text = first;
                }
            }
        }

        let abbreviation: String = if let Some(me) = mobile_element(text) {
            me
        } else if contains_ci(text, "DEL") && contains_ci(text, "DUP") {
            "CNV".into()
        } else if contains_ci(text, "DUP:TANDEM") {
            "TDUP".into()
        } else if contains_ci(text, "CNV:TR") {
            "TREP".into()
        } else if has_copy_number_allele(text) {
            // `/^<?CN=?0>?$/` and `/^<?CN=?2>?$/` carry no `i`: they are upper case only.
            if is_exact_copy_number(text, b'0') {
                "DEL".into()
            } else if is_exact_copy_number(text, b'2') {
                "DUP".into()
            } else {
                "CNV".into()
            }
        } else if contains_ci(text, "CNV") {
            "CNV".into()
        } else if text.contains(['[', ']']) || text.starts_with('.') || text.ends_with('.') {
            "BND".into()
        } else if text.starts_with('<') || text.ends_with('>') {
            let bare: String = text.chars().filter(|c| *c != '<' && *c != '>').collect();
            // `s/\:.+//g`: from the first colon that has something after it
            match bare
                .char_indices()
                .find(|&(i, c)| c == ':' && i + 1 < bare.len())
            {
                Some((at, _)) => bare[..at].to_string(),
                None => bare,
            }
        } else {
            text.to_string()
        };

        Self::from_abbreviation(&abbreviation).ok_or(abbreviation)
    }

    /// The table of `%SO_TERMS`, by the abbreviation `get_SO_term` reduces a type to.
    fn from_abbreviation(abbreviation: &str) -> Option<Self> {
        Some(match abbreviation {
            "INS" => Self::Insertion,
            "INS_ME" => Self::MobileElementInsertion,
            "INS_ALU" => Self::AluInsertion,
            "INS_HERV" => Self::HervInsertion,
            "INS_LINE1" => Self::Line1Insertion,
            "INS_SVA" => Self::SvaInsertion,
            "DEL" => Self::Deletion,
            "DEL_ME" => Self::MobileElementDeletion,
            "DEL_ALU" => Self::AluDeletion,
            "DEL_HERV" => Self::HervDeletion,
            "DEL_LINE1" => Self::Line1Deletion,
            "DEL_SVA" => Self::SvaDeletion,
            "TREP" => Self::TandemRepeat,
            "TDUP" => Self::TandemDuplication,
            "DUP" => Self::Duplication,
            "CNV" => Self::CopyNumberVariation,
            "INV" => Self::Inversion,
            "BND" => Self::ChromosomeBreakpoint,
            _ => return None,
        })
    }
}

/// `/needle/i` for an ASCII needle.
fn contains_ci(haystack: &str, needle: &str) -> bool {
    let (h, n) = (haystack.as_bytes(), needle.as_bytes());
    h.len() >= n.len() && h.windows(n.len()).any(|w| w.eq_ignore_ascii_case(n))
}

/// `/(INS|DEL):(ME):?([A-Z0-9]+)?/i` of `get_SO_term`, read the way that code reads its captures:
/// `INS` or `DEL`, then the element (`ALU`, `HERV`, `LINE1` with `L1` as its other name, `SVA`) if
/// there is one, or plain `ME`.
fn mobile_element(text: &str) -> Option<String> {
    let bytes = text.as_bytes();
    let found = (0..bytes.len()).find_map(|i| {
        let rest = bytes.get(i..)?;
        let kind = if rest.len() >= 3 && rest[..3].eq_ignore_ascii_case(b"INS") {
            "INS"
        } else if rest.len() >= 3 && rest[..3].eq_ignore_ascii_case(b"DEL") {
            "DEL"
        } else {
            return None;
        };
        (rest.get(3) == Some(&b':') && rest.len() >= 6 && rest[4..6].eq_ignore_ascii_case(b"ME"))
            .then(|| (kind, &rest[6..]))
    })?;
    let (kind, after) = found;
    let after = after.strip_prefix(b":").unwrap_or(after);
    let element: String = after
        .iter()
        .take_while(|b| b.is_ascii_alphanumeric())
        .map(|b| b.to_ascii_uppercase() as char)
        .collect();
    let element = if element == "L1" {
        "LINE1"
    } else {
        element.as_str()
    };
    let subtype = match element {
        "ALU" | "HERV" | "LINE1" | "SVA" => element,
        _ => "ME",
    };
    Some(format!("{kind}_{subtype}"))
}

/// `/CN=?[0-9]/i`: a copy number allele such as `<CN0>` or `<CN=3>`.
fn has_copy_number_allele(text: &str) -> bool {
    let b = text.as_bytes();
    (0..b.len().saturating_sub(2)).any(|i| {
        b[i..i + 2].eq_ignore_ascii_case(b"CN")
            && (b[i + 2].is_ascii_digit()
                || (b[i + 2] == b'=' && b.get(i + 3).is_some_and(u8::is_ascii_digit)))
    })
}

/// `/^<?CN=?D>?$/` for one digit `D`, upper case only.
fn is_exact_copy_number(text: &str, digit: u8) -> bool {
    let t = text.strip_prefix('<').unwrap_or(text);
    let t = t.strip_suffix('>').unwrap_or(t);
    let t = t.strip_prefix("CN").unwrap_or(t);
    let t = t.strip_prefix('=').unwrap_or(t);
    // the prefix must really have been there
    text.contains("CN") && t.as_bytes() == [digit]
}

/// Where an allele of a structural variant is, for the one rule that is about a position of its
/// own: `feature_truncation` of a breakpoint holds when the allele's breakend lies in the
/// transcript.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AlleleSite {
    /// An allele with no position of its own: every allele of a variant that is not a breakend, and
    /// a single breakend. For a single breakend the question of `feature_truncation` is asked of the
    /// variant's own position, which is the same place as [`Self::Own`].
    None,
    /// The variant itself, the first allele of a breakend with mates: its own interval.
    Own,
    /// The far end of a mate (`A[chr2:321682[` is `chr2` position `321682`).
    Mate { chromosome: Arc<str>, position: u64 },
}

/// One allele of a structural variant: a row for each transcript it is close to.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SvAllele {
    /// What the `Allele` column says: the class (`deletion`), or for a breakend with mates its own
    /// string (`A.`, `A[chr2:321682[`).
    pub label: Arc<str>,
    pub site: AlleleSite,
    /// The `ALT` of the record this allele answers for, which is what the supplementary
    /// annotations of a variant are keyed by.
    pub alt: usize,
}

/// Which allele of a structural variant a row answers for: a position in
/// [`StructuralVariant::alleles`], where the allele's label is kept once for the variant.
///
/// A row carries this, four bytes, and not the label, which is an `Arc<str>` of sixteen. Every
/// row of every transcript of every variant has the field, and only a structural variant's rows
/// have a label, so the label's width is what every small variant would pay for: with the label in
/// the row, `AlleleAnnotation` was 648 bytes instead of 632 and a one-thread run over 202,417
/// variants with VCF output took 2.3% more user time (3.48 s to 3.56 s, medians of 8 runs).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct AlleleNo(NonZeroU32);

impl AlleleNo {
    /// The number of the allele at `index` of [`StructuralVariant::alleles`]; `None` for an index
    /// that does not fit, which a variant cannot have: the alleles of one are bounded long before.
    pub fn of(index: usize) -> Option<Self> {
        let number = u32::try_from(index).ok()?.checked_add(1)?;
        NonZeroU32::new(number).map(Self)
    }

    /// The index in [`StructuralVariant::alleles`].
    pub fn index(self) -> usize {
        self.0.get() as usize - 1
    }
}

/// A tandem repeat that VEP turns into alleles it can read as bases: what its INFO says of the
/// repeat unit and how often it occurs (`RUS`, `RUC` or `RB`, `RN`), kept as written because the
/// reference it is placed on is not known until the record is annotated.
///
/// `Bio::EnsEMBL::VEP::Parser::VCF::_expand_tandem_repeat_allele_string` builds, for each record of
/// `RN`, the unit repeated `RUC` times (or `RB` divided by the unit's length), and annotates the
/// repeat as a small variant: its reference is the bases of the interval, its alleles are these.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TandemRepeat {
    /// `RUS`: the repeat units, one for each unit of each record, `,` between them
    pub units: String,
    /// `RUC`: how many times each unit occurs, or empty
    pub counts: String,
    /// `RB`: how many bases each unit makes, or empty
    pub bases: String,
    /// `RN`: how many units each allele has, `.` for an allele that has none; absent for one
    pub records: Option<String>,
    /// How many times `CNV:TR` is in the type, which is how many alleles there are with no `RN`
    pub types: usize,
}

impl TandemRepeat {
    /// The most alleles a repeat is read as bases for. VEP has no limit, and a record whose `RN`
    /// lists a million counts is a million alleles for each transcript it meets: 2 MB of INFO was
    /// 1.1 GB of rows in a fuzz of the program, and with distinct units the search for an allele that
    /// repeats an earlier one grows with the square of the count (30,000 distinct units took a second
    /// and 250,000 took 81). A request to the web server holds its plain ALT list to 256, and a record
    /// that lists more repeats than this stays the structural variant it is.
    pub const MAX_ALLELES: usize = 1_024;

    /// How many alleles the repeat has: one for each record of `RN` that is not `.`, and with no `RN`
    /// one for each `CNV:TR` of the type.
    pub fn allele_count(&self) -> usize {
        match &self.records {
            Some(rn) => rn.split(',').filter(|record| *record != ".").count(),
            None => self.types,
        }
    }

    /// The repeat alleles (`_expand_tandem_repeat_allele_string`), or `None` when one is longer than
    /// `max_size` bases (`None` for none: VEP's `--max_sv_size -1`, for which every repeat is too
    /// long), when the repeat has more than [`Self::MAX_ALLELES`], or when the INFO gives `RB` and no
    /// unit to divide it by. The alleles are in upper case, as VEP writes them whatever case the unit
    /// was given in.
    ///
    /// VEP bounds each run of a unit and nothing else, so a record of sixty alleles of ten megabases
    /// each is six hundred megabytes of repeat. The alleles together are bounded here by the same
    /// limit, which no record that VEP reads as bases comes near.
    pub fn alleles(&self, max_size: Option<u64>) -> Option<Vec<String>> {
        if self.allele_count() > Self::MAX_ALLELES {
            return None;
        }
        let mut units = self.units.split(',');
        let mut counts = self.counts.split(',').filter(|c| !c.is_empty());
        let mut bases = self.bases.split(',').filter(|c| !c.is_empty());
        let has_counts = !self.counts.is_empty();
        let records: Vec<&str> = match &self.records {
            Some(rn) => rn.split(',').collect(),
            None => vec!["1"; self.types],
        };
        let too_long = |size: u128| max_size.is_none_or(|limit| size > u128::from(limit));
        let mut alleles = Vec::with_capacity(records.len());
        let mut total = 0u128;
        for record in records {
            if record == "." {
                continue;
            }
            let repeats = perl_count(record);
            // with no limit every run is too long, the first one that is made included
            if repeats > 0 && max_size.is_none() {
                return None;
            }
            let mut allele = String::new();
            for _ in 0..repeats {
                // A list that has run out gives Perl's `undef`: no unit, or a count of none, so
                // the repeat adds nothing (`RN=2` over one `RUC` is the first repeat alone), and
                // nothing after it can add anything either, since the units do not come back.
                // Only `RB`, which is divided by the unit's length, cannot do without a unit:
                // there Perl divides by zero and VEP stops, and here the record stays what it was.
                let Some(unit) = units.next() else {
                    break;
                };
                let unit = if unit == "." { "N" } else { unit };
                let number = if has_counts {
                    counts.next().map_or(0, perl_count)
                } else if unit.is_empty() {
                    return None;
                } else {
                    bases.next().map_or(0, perl_count) / unit.len()
                };
                let size = unit.len() as u128 * number as u128;
                if too_long(size) || too_long(total + allele.len() as u128 + size) {
                    return None;
                }
                allele.push_str(&unit.to_ascii_uppercase().repeat(number));
            }
            total += allele.len() as u128;
            alleles.push(allele);
        }
        Some(alleles)
    }
}

/// A whole number out of a string, as Perl reads it where a count is wanted: its leading digits, or
/// nothing (zero), and past the largest number that fits, that one.
fn perl_count(text: &str) -> usize {
    text.trim()
        .bytes()
        .take_while(u8::is_ascii_digit)
        .fold(0usize, |n, digit| {
            n.saturating_mul(10)
                .saturating_add(usize::from(digit - b'0'))
        })
}

/// A structural variant: its class, and the alleles VEP writes rows for.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StructuralVariant {
    pub class: SvClass,
    pub alleles: Vec<SvAllele>,
    /// Set for a tandem repeat that has the INFO to be read as bases ([`TandemRepeat`])
    pub tandem: Option<TandemRepeat>,
}

impl StructuralVariant {
    /// What the `Allele` column says for the allele a row answers for.
    pub fn label(&self, allele: AlleleNo) -> Option<&str> {
        self.alleles.get(allele.index()).map(|a| &*a.label)
    }

    /// The positions the variant has besides its own, for the transcripts that must be looked for
    /// there: a breakend's mates.
    pub fn mates(&self) -> impl Iterator<Item = (&str, u64)> + '_ {
        self.alleles.iter().filter_map(|a| match &a.site {
            AlleleSite::Mate {
                chromosome,
                position,
            } => Some((&**chromosome, *position)),
            _ => None,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn class(text: &str) -> Result<&'static str, String> {
        SvClass::from_type(text).map(SvClass::so_term)
    }

    #[test]
    fn an_allele_number_is_the_index_it_was_made_from_and_costs_a_row_four_bytes() {
        assert_eq!(AlleleNo::of(0).map(AlleleNo::index), Some(0));
        assert_eq!(AlleleNo::of(41).map(AlleleNo::index), Some(41));
        // the last index that has a number after it does not exist: no variant has that many alleles
        assert_eq!(AlleleNo::of(u32::MAX as usize), None);
        assert_eq!(std::mem::size_of::<Option<AlleleNo>>(), 4);
    }

    #[test]
    fn the_label_of_a_number_is_that_of_the_allele_at_its_index() {
        let allele = |label: &str| SvAllele {
            label: Arc::from(label),
            site: AlleleSite::None,
            alt: 0,
        };
        let sv = StructuralVariant {
            class: SvClass::Deletion,
            alleles: vec![allele("A."), allele("deletion")],
            tandem: None,
        };
        assert_eq!(sv.label(AlleleNo::of(0).unwrap()), Some("A."));
        assert_eq!(sv.label(AlleleNo::of(1).unwrap()), Some("deletion"));
        assert_eq!(sv.label(AlleleNo::of(2).unwrap()), None);
    }

    #[test]
    fn the_plain_types_are_their_terms() {
        assert_eq!(class("<DEL>"), Ok("deletion"));
        assert_eq!(class("<DUP>"), Ok("duplication"));
        assert_eq!(class("<INS>"), Ok("insertion"));
        assert_eq!(class("<INV>"), Ok("inversion"));
        assert_eq!(class("<CNV>"), Ok("copy_number_variation"));
        assert_eq!(class("<BND>"), Ok("chromosome_breakpoint"));
        assert_eq!(class("DEL"), Ok("deletion"));
    }

    #[test]
    fn a_breakend_in_any_notation_is_a_breakpoint() {
        for alt in ["A[chr1:2500[", "]chr1:2500]A", "A]2:5]", ".A", "A."] {
            assert_eq!(class(alt), Ok("chromosome_breakpoint"), "{alt}");
        }
        // a name with a dot in it is a name, not a single breakend
        assert_eq!(class("<DEL>"), Ok("deletion"));
    }

    #[test]
    fn a_mobile_element_names_its_family_and_unknown_ones_are_plain_me() {
        assert_eq!(class("<INS:ME:ALU>"), Ok("Alu_insertion"));
        assert_eq!(class("<INS:ME:Alu>"), Ok("Alu_insertion"));
        assert_eq!(class("<DEL:ME:L1>"), Ok("LINE1_deletion"));
        assert_eq!(class("<INS:ME:HERV>"), Ok("HERV_insertion"));
        assert_eq!(class("<INS:ME:SVA>"), Ok("SVA_insertion"));
        assert_eq!(class("<INS:ME>"), Ok("mobile_element_insertion"));
        assert_eq!(class("<DEL:ME:XYZ>"), Ok("mobile_element_deletion"));
    }

    #[test]
    fn copy_number_alleles_follow_the_number() {
        assert_eq!(class("<CN0>"), Ok("deletion"));
        assert_eq!(class("<CN=0>"), Ok("deletion"));
        assert_eq!(class("<CN2>"), Ok("duplication"));
        assert_eq!(class("<CN1>"), Ok("copy_number_variation"));
        assert_eq!(class("<CN3>"), Ok("copy_number_variation"));
        assert_eq!(class("<CN9>"), Ok("copy_number_variation"));
        // several alleles together are not CN0 or CN2 alone
        assert_eq!(class("<CN0>,<CN2>"), Ok("copy_number_variation"));
    }

    #[test]
    fn tandem_types_come_before_their_neighbours() {
        assert_eq!(class("<DUP:TANDEM>"), Ok("tandem_duplication"));
        assert_eq!(class("<CNV:TR>"), Ok("tandem_repeat"));
    }

    #[test]
    fn a_deletion_beside_a_duplication_is_a_copy_number_variation() {
        assert_eq!(class("<DEL>/<DUP>"), Ok("copy_number_variation"));
        assert_eq!(class("<DEL>,<DUP>"), Ok("copy_number_variation"));
    }

    #[test]
    fn a_type_that_is_repeated_is_the_type_once() {
        assert_eq!(class("<DEL>/<DEL>"), Ok("deletion"));
        assert_eq!(class("<INS>,<INS>"), Ok("insertion"));
        // but not when what is beside it is a base: `G/<DEL>` is judged on its own text, which is
        // why such a record takes its SVTYPE instead
        assert_eq!(class("G/<DEL>"), Err("G/DEL".to_string()));
    }

    #[test]
    fn a_type_that_vep_does_not_know_is_refused_with_its_name() {
        assert_eq!(class("<CPX>"), Err("CPX".into()));
        assert_eq!(class("<NON_REF>"), Err("NON_REF".into()));
        assert_eq!(class("<CTX>"), Err("CTX".into()));
        assert_eq!(class("<STR>"), Err("STR".into()));
        // the lookup is case-sensitive
        assert_eq!(class("<del>"), Err("del".into()));
        // what follows the first colon is dropped from the name
        assert_eq!(class("<CPX:FOO>"), Err("CPX".into()));
    }

    #[test]
    fn the_rules_read_the_class_the_way_the_perl_regexes_do() {
        use SvClass::*;
        for c in [
            Deletion,
            MobileElementDeletion,
            AluDeletion,
            HervDeletion,
            Line1Deletion,
            SvaDeletion,
        ] {
            assert!(c.is_deletion() && c.decreases_length() && !c.is_insertion());
        }
        for c in [
            Insertion,
            MobileElementInsertion,
            AluInsertion,
            HervInsertion,
            Line1Insertion,
            SvaInsertion,
            Duplication,
            TandemDuplication,
            TandemRepeat,
        ] {
            assert!(c.is_insertion() && !c.is_deletion() && !c.decreases_length());
        }
        // a copy number gain is a duplication or a tandem repeat, and an insertion is not one
        assert!(Duplication.is_copy_number_gain() && TandemRepeat.is_copy_number_gain());
        assert!(!Insertion.is_copy_number_gain() && !AluInsertion.is_copy_number_gain());
        // a breakpoint shortens, and a bare copy number variation or inversion does neither
        assert!(ChromosomeBreakpoint.decreases_length() && !ChromosomeBreakpoint.is_deletion());
        for c in [CopyNumberVariation, Inversion] {
            assert!(!c.decreases_length() && !c.is_insertion() && !c.is_deletion());
        }
    }

    fn repeat(units: &str, counts: &str, records: Option<&str>) -> super::TandemRepeat {
        super::TandemRepeat {
            units: units.into(),
            counts: counts.into(),
            bases: String::new(),
            records: records.map(String::from),
            types: 0,
        }
    }

    #[test]
    fn a_repeat_is_its_unit_as_many_times_as_it_is_counted() {
        let tr = repeat("CA,GT", "3,2", Some("2"));
        assert_eq!(tr.alleles(Some(1000)), Some(vec!["CACACAGTGT".to_string()]));
        // a count the INFO does not have is no repeat, and a unit it does not have adds nothing
        assert_eq!(
            repeat("CA,GT", "3", Some("2")).alleles(Some(1000)),
            Some(vec!["CACACA".to_string()])
        );
        assert_eq!(
            repeat("CA", "3", Some("2")).alleles(Some(1000)),
            Some(vec!["CACACA".to_string()])
        );
    }

    #[test]
    fn a_count_that_does_not_fit_in_a_number_is_too_long_and_not_a_panic() {
        for count in [
            "9223372036854775807",
            "18446744073709551615",
            "18446744073709551616",
            "99999999999999999999999999",
        ] {
            assert_eq!(
                repeat("CA", count, Some("1")).alleles(Some(10_000_000)),
                None,
                "{count}"
            );
        }
    }

    #[test]
    fn more_units_than_there_are_ends_with_the_ones_there_are() {
        // Perl would go round this loop 18 quintillion times; the units are gone after one
        let tr = repeat("CA", "3", Some("18446744073709551615"));
        assert_eq!(tr.alleles(Some(1000)), Some(vec!["CACACA".to_string()]));
    }

    #[test]
    fn the_alleles_together_are_bounded_by_the_limit() {
        let units = vec!["ACGT"; 60].join(",");
        let counts = vec!["25"; 60].join(",");
        let records = vec!["1"; 60].join(",");
        let tr = repeat(&units, &counts, Some(&records));
        // each allele is 100 bases and the sixty are 6,000
        assert_eq!(tr.alleles(Some(6000)).map(|a| a.len()), Some(60));
        assert_eq!(tr.alleles(Some(5999)), None);
    }

    #[test]
    fn with_no_limit_every_repeat_is_too_long() {
        assert_eq!(repeat("CA", "3", Some("1")).alleles(None), None);
        // and a record that has no repeat in it has none to be too long
        assert_eq!(repeat("CA", "3", Some(".")).alleles(None), Some(vec![]));
    }

    /// VEP upper-cases the allele string it builds (`Parser.pm`), so a unit written in lower case
    /// is the same repeat as its upper-case twin.
    #[test]
    fn a_unit_is_written_in_upper_case() {
        assert_eq!(
            repeat("ca", "3", Some("1")).alleles(Some(1000)),
            Some(vec!["CACACA".to_string()])
        );
    }

    #[test]
    fn a_repeat_of_more_than_the_most_alleles_stays_what_it_is() {
        let most = super::TandemRepeat::MAX_ALLELES;
        let many = |n: usize| {
            repeat(
                &vec!["A"; n].join(","),
                &vec!["1"; n].join(","),
                Some(&vec!["1"; n].join(",")),
            )
        };
        assert_eq!(many(most + 1).allele_count(), most + 1);
        assert_eq!(many(most + 1).alleles(Some(10_000_000)), None);
        assert_eq!(
            many(most).alleles(Some(10_000_000)).map(|a| a.len()),
            Some(most)
        );
        // an allele that has no records (`.`) is not one
        assert_eq!(repeat("CA", "3", Some(".,.,1")).allele_count(), 1);
    }
}
