//! A VCF record that Ensembl VEP reads as a structural variant.
//!
//! This is the part of VEP that decides *what a record is*, before any consequence is asked of it:
//! whether it is a structural variant at all (`create_VariationFeatures` of `Parser/VCF.pm`), its
//! class (`get_SO_term`, in [`fastvep_core::SvClass::from_type`]), its interval (`get_start` and
//! `get_end` of `Bio::EnsEMBL::IO::Parser::BaseVCF4`), the alleles it writes a row for
//! (`StructuralVariationOverlap.pm`) and the records it refuses (`create_StructuralVariationFeatures`
//! and `validate_vf`). Every step is a transcription, and where a Perl regex is what decides, the
//! function that stands for it names it.
//!
//! The interval is the thing to know. For every symbolic class, `<DUP>` and `<INV>` and `<INS>`
//! included, the start is `POS + 1`; the end is `START + |SVLEN| - 1` when there is an `SVLEN`
//! (which wins over `END`), else `END`, else the start. A deletion needs one of them or it is
//! skipped.
//!
//! One deliberate difference: the placeholders `.`, `<NON_REF>` and `<*>` that a gVCF puts beside
//! (or in place of) an ALT are not read as alleles of the record, so `C,<NON_REF>` is the SNV it
//! always was here. VEP reads that record as a structural variant of unknown type and drops it.

use std::sync::Arc;

use fastvep_core::{AlleleSite, StructuralVariant, SvAllele, SvClass, TandemRepeat};

/// A record Ensembl VEP skips with a warning, as opposed to one it cannot read. The reader reports
/// it and goes on, whether it is lenient or not.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Skipped(pub String);

impl std::fmt::Display for Skipped {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for Skipped {}

/// A record that is a structural variant, as VEP holds it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SvRecord {
    /// `start` and `end` of the `StructuralVariationFeature`. An insertion has `start = end + 1`.
    pub start: u64,
    pub end: u64,
    /// `allele_string`: the ALTs joined with `/` (`<DEL>/<DUP>`), or for a breakend the reference
    /// and the breakends (`A/A[chr1:2500[`), which is what VEP's JSON calls the `allele_string`.
    pub allele_string: String,
    pub sv: StructuralVariant,
}

/// The ALTs that are alleles: those that assert something. A placeholder is not one.
pub fn is_placeholder(alt: &str) -> bool {
    matches!(alt, "." | "<NON_REF>" | "<*>")
}

/// `/[<\[\]][^\*]+[>\]\[]|^\.\w+|\w+\.$/` of `create_VariationFeatures`, the test of the joined ALTs
/// that makes a record an SV when it has no `SVTYPE`.
///
/// The first branch wants a `<`, `[` or `]`, then one or more characters that are not `*`, then a
/// `>`, `]` or `[`: `<DEL>` and `A[chr1:2500[` have it, `<*>` does not (nothing but `*` between its
/// brackets). It is read in one pass, keeping the earliest opener since the last `*`.
fn looks_symbolic(alts: &str) -> bool {
    let b = alts.as_bytes();
    let mut opener: Option<usize> = None;
    for (j, &c) in b.iter().enumerate() {
        if c == b'*' {
            opener = None;
            continue;
        }
        if matches!(c, b'>' | b']' | b'[') && opener.is_some_and(|i| i + 2 <= j) {
            return true;
        }
        if matches!(c, b'<' | b'[' | b']') && opener.is_none() {
            opener = Some(j);
        }
    }
    let is_word = |c: u8| c.is_ascii_alphanumeric() || c == b'_';
    // `^\.\w+`
    if b.first() == Some(&b'.') && b.get(1).is_some_and(|&c| is_word(c)) {
        return true;
    }
    // `\w+\.$`
    b.len() >= 2 && b[b.len() - 1] == b'.' && is_word(b[b.len() - 2])
}

/// Perl's `split('=', $piece)` read as `my ($key, $value) = ...`: the key, and the value unless the
/// piece has no `=` at all (a flag). Assigning to two scalars splits into three fields at most, and
/// a split with a limit keeps its trailing empty fields, so `KEY=` has a value, the empty string,
/// which is defined where a flag's is not (checked with perl 5.34: `SVLEN=` is defined and empty,
/// so a deletion with it has a length of nothing and is skipped as incomplete, where one without it
/// keeps its `END`).
fn key_value(piece: &str) -> (&str, Option<&str>) {
    let mut parts = piece.split('=');
    (parts.next().unwrap_or(""), parts.next())
}

/// The value of an INFO key as VEP's `get_info` holds it: the last occurrence wins, and a key with
/// no `=` (a flag) is undefined; `KEY=` is defined and empty.
fn info_value<'a>(info: &'a str, key: &str) -> Option<&'a str> {
    let mut found = None;
    for piece in info.split(';') {
        let (k, v) = key_value(piece);
        if k == key {
            found = v;
        }
    }
    found
}

/// Perl truth of a defined value: neither empty nor `0`.
fn truthy(value: Option<&str>) -> bool {
    value.is_some_and(|v| !v.is_empty() && v != "0")
}

/// A Perl number out of a string: its leading integer, or 0 when there is none (`abs('.')`), and
/// past the largest that fits, that one: Perl reads `99999999999999999999` as a number that is
/// larger than any size limit, and a length of nothing is not what it says. A value with a fraction
/// or an exponent is read to its integer part.
fn perl_integer(text: &str) -> i64 {
    let t = text.trim_start();
    let (sign, digits) = match t.as_bytes().first() {
        Some(b'-') => (-1, &t[1..]),
        Some(b'+') => (1, &t[1..]),
        _ => (1, t),
    };
    let magnitude = digits
        .bytes()
        .take_while(u8::is_ascii_digit)
        .fold(0i64, |n, digit| {
            n.saturating_mul(10).saturating_add(i64::from(digit - b'0'))
        });
    sign * magnitude
}

/// `$alt !~ /^<?(DEL|INS|DUP|INV|CNV|CN=?[0-9]+)/`, which decides whether `SVTYPE` replaces the ALT
/// as the type: it does not when the ALT already starts with a type, in upper case.
fn starts_with_a_type(alt: &str) -> bool {
    let t = alt.strip_prefix('<').unwrap_or(alt);
    ["DEL", "INS", "DUP", "INV", "CNV"]
        .iter()
        .any(|name| t.starts_with(name))
        || t.strip_prefix("CN").is_some_and(|rest| {
            let rest = rest.strip_prefix('=').unwrap_or(rest);
            rest.as_bytes().first().is_some_and(u8::is_ascii_digit)
        })
}

/// Whether Ensembl VEP reads this record as a structural variant (`create_VariationFeatures`).
///
/// `join("", REF, ALTs)` must not be bases alone, and then either `SVTYPE` is set (true in Perl's
/// sense: not empty, not `0`) or the ALTs joined with `,` look symbolic ([`looks_symbolic`]).
/// `alts` are the alleles that assert something; see [`is_placeholder`].
pub fn is_structural<'a>(
    ref_given: &'a str,
    alts: impl Iterator<Item = &'a str> + Clone,
    info: &str,
) -> bool {
    if is_acgt_joined(ref_given, alts.clone()) {
        return false;
    }
    // Most records that are not bases alone are `*`, `N`, lower case or IUPAC alleles, which have
    // neither an `SVTYPE` nor anything symbolic: they are told from the structural ones without
    // reading the INFO piece by piece or building the joined ALTs.
    if info.contains("SVTYPE") && truthy(info_value(info, "SVTYPE")) {
        return true;
    }
    let (first, last) = (alts.clone().next(), alts.clone().last());
    let may_look_symbolic = alts
        .clone()
        .any(|alt| alt.bytes().any(|c| matches!(c, b'<' | b'[' | b']')))
        || first.is_some_and(|alt| alt.starts_with('.'))
        || last.is_some_and(|alt| alt.ends_with('.'));
    may_look_symbolic && looks_symbolic(&alts.collect::<Vec<_>>().join(","))
}

/// `join("", $ref, @$alts) =~ /^[ACGT]+$/` over the ALT column as it is written, the alleles with a
/// comma between them, in one pass over the bytes of the record and nothing else. Nearly every
/// record is bases alone, and is told from the rest here, before any iterator over its alleles is
/// made: this runs once for every record of a run, on the one thread that reads the file.
pub fn is_bases_only(ref_given: &str, alt_field: &str) -> bool {
    let is_base = |b: u8| matches!(b, b'A' | b'C' | b'G' | b'T');
    ref_given.bytes().all(is_base) && alt_field.bytes().all(|b| is_base(b) || b == b',')
}

/// `join("", $ref, @$alts) =~ /^[ACGT]+$/` without building the string.
fn is_acgt_joined<'a>(ref_given: &'a str, alts: impl Iterator<Item = &'a str>) -> bool {
    let mut any = false;
    for part in std::iter::once(ref_given).chain(alts) {
        if !part.bytes().all(|b| matches!(b, b'A' | b'C' | b'G' | b'T')) {
            return false;
        }
        any |= !part.is_empty();
    }
    any
}

/// One mate of a breakend, as `_parse_breakends` reads a piece of the allele string with
/// `([A-Za-z]*)\s*[\[\]]\s*(\S+)\s*:\s*([0-9]+)\s*[\[\]]\s*([(A-Za-z)]*)`: the contig and the position
/// between its brackets, and whether the piece begins with a bracket (a mate placed to the left,
/// which is `$alt_string =~ '^(\[|\])'`).
///
/// Read without spaces, which no VCF allele holds. The contig is the longest run that can precede a
/// `:digits` and a closing bracket, so a contig with colons in it (`HLA-A*01:01:01:01`) keeps them.
fn parse_mate(piece: &str) -> Option<(&str, u64, bool)> {
    let b = piece.as_bytes();
    if b.iter().any(u8::is_ascii_whitespace) {
        return None;
    }
    let is_bracket = |c: u8| c == b'[' || c == b']';
    // The contig runs to the last colon that has digits and then a bracket after it (`\S+` takes as
    // much as it can), and it is the same colon for every bracket that could open the piece, so it
    // is found once, from the right: a piece of `[[[[...` is not read again from each bracket.
    let colon = (0..b.len()).rev().find(|&c| {
        b[c] == b':' && {
            let digits = b[c + 1..].iter().take_while(|d| d.is_ascii_digit()).count();
            digits > 0 && b.get(c + 1 + digits).copied().is_some_and(is_bracket)
        }
    })?;
    // The first bracket opens the piece, the letters before it being the reference base(s), and
    // needs at least one character of contig between it and the colon.
    let open = b.iter().position(|&c| is_bracket(c))?;
    if colon < open + 2 {
        return None;
    }
    let digits = b[colon + 1..]
        .iter()
        .take_while(|d| d.is_ascii_digit())
        .count();
    let position = piece[colon + 1..colon + 1 + digits].parse::<u64>().ok()?;
    Some((&piece[open + 1..colon], position, is_bracket(b[0])))
}

/// `^[\.]?([ACGTN]+)[\.]?$`: the reference base(s) of a breakend's own piece.
fn reference_of_piece(piece: &str) -> Option<&str> {
    let t = piece.strip_prefix('.').unwrap_or(piece);
    let t = t.strip_suffix('.').unwrap_or(t);
    (!t.is_empty()
        && t.bytes()
            .all(|b| matches!(b, b'A' | b'C' | b'G' | b'N' | b'T')))
    .then_some(t)
}

/// Read a record that [`is_structural`] says is one: `alts` are `(index in the ALT column, text)`
/// of the alleles that assert something, in order.
///
/// `Err` is a record VEP skips: a type it does not know, a deletion with no extent, an interval that
/// ends before it starts, a coordinate that is not a number.
///
/// `insertion_as_point` reads an insertion as the point between `POS` and `POS + 1` (the zero-length
/// interval `start = end + 1`) whatever its `SVLEN` and `END` say, which is not what VEP does: see
/// [`fastvep_io::vcf::ParseOptions::insertion_as_point`](crate::vcf::ParseOptions).
pub fn read(
    pos: u64,
    ref_given: &str,
    alts: &[(usize, &str)],
    info: &str,
    insertion_as_point: bool,
) -> Result<SvRecord, Skipped> {
    let alt_slash: String = alts.iter().map(|a| a.1).collect::<Vec<_>>().join("/");
    let (class, type_text) = class_of(info, &alt_slash)?;

    // `get_start`: the base before the variant is the padding of the record, for anything that is not
    // bases alone
    let sv_like =
        info.contains("SVTYPE") || alts.iter().any(|a| a.1.contains(['<', '[', ']', '>']));
    let start = if sv_like || alts.iter().any(|a| a.1.len() != ref_given.len()) {
        pos.saturating_add(1)
    } else {
        pos
    };
    let Interval {
        mut end,
        has_extent,
    } = end_of(class, pos, ref_given, start, sv_like, info)?;
    if insertion_as_point && class.is_insertion() && !class.is_copy_number_gain() {
        end = i128::from(start) - 1;
    }
    let allele_string = allele_string_of(class, &alt_slash, ref_given, info);

    // The sums are made in `i128`, which holds any position and any length a record can name, and
    // the interval is judged there
    let start_wide = i128::from(start);
    if class.is_deletion() && (start_wide > end || !has_extent) {
        return Err(Skipped("deletion looks incomplete".into()));
    }
    if start_wide > end.saturating_add(1) {
        return Err(Skipped(format!(
            "start > end+1 : (START={start}, END={end})"
        )));
    }
    // and then brought into the range of a coordinate: at most what `i64` holds, so that no sum of
    // one with a distance or a length can wrap. The record that names more is the web server's to
    // refuse, and has no transcript here.
    let start = start.min(i64::MAX as u64);
    let end = end.clamp(0, i128::from(i64::MAX)) as u64;

    let alleles = alleles_of(class, &allele_string, alts);
    Ok(SvRecord {
        start,
        end,
        allele_string,
        sv: StructuralVariant {
            class,
            alleles,
            tandem: tandem_of(class, info, type_text),
        },
    })
}

/// The class of a record, and the type it was read from. The type is the ALT, unless `SVTYPE` is set
/// and the ALT does not already start with one; a type VEP has no class for is a skipped record.
fn class_of<'a>(info: &'a str, alt_slash: &'a str) -> Result<(SvClass, &'a str), Skipped> {
    let type_text = match info_value(info, "SVTYPE") {
        Some(svtype) if truthy(Some(svtype)) && !starts_with_a_type(alt_slash) => svtype,
        _ => alt_slash,
    };
    let class = SvClass::from_type(type_text).map_err(|abbreviation| {
        Skipped(format!(
            "{abbreviation} is not a supported structural variant type"
        ))
    })?;
    Ok((class, type_text))
}

/// The end of a record's interval, and whether the record said how far it goes at all.
struct Interval {
    end: i128,
    /// an `SVLEN` or an `END`: a deletion that has neither is "incomplete"
    has_extent: bool,
}

/// `get_end`. `SVLEN` wins over `END`, then `END`, then the one base after the padding of a record
/// that is not bases alone, then the bases of the reference.
///
/// A breakpoint's INFO/END is Manta's mate position and is taken out first, if it is true: Perl tests
/// the value, so an `END=0` stays where it is, and is the end.
fn end_of(
    class: SvClass,
    pos: u64,
    ref_given: &str,
    start: u64,
    sv_like: bool,
    info: &str,
) -> Result<Interval, Skipped> {
    let start_wide = i128::from(start);
    let svlen = info_value(info, "SVLEN");
    let end_info = info_value(info, "END");
    let end_text = if class.is_breakpoint() && truthy(end_info) {
        None
    } else {
        end_info
    };
    let end = if let Some(len) = svlen {
        // A tandem repeat takes the longest of its `SVLEN` values for its length, as it is, and
        // anything else the size of the first one.
        let length = if is_tandem(class) {
            len.split(',').map(perl_integer).max().unwrap_or(0)
        } else {
            perl_integer(len.split(',').next().unwrap_or("")).abs()
        };
        start_wide + i128::from(length) - 1
    } else if let Some(text) = end_text {
        // a number VEP reads as one, however large: it is the size limit's to refuse
        text.trim().parse::<i128>().map_err(|_| {
            Skipped(format!(
                "Invalid start '{start}' or end '{text}' coordinate"
            ))
        })?
    } else if sv_like {
        start_wide
    } else {
        i128::from(pos) + ref_given.len() as i128 - 1
    };
    Ok(Interval {
        end,
        has_extent: svlen.is_some() || end_text.is_some(),
    })
}

/// The classes that VEP reads as bases when the INFO has what it takes.
fn is_tandem(class: SvClass) -> bool {
    matches!(class, SvClass::TandemRepeat | SvClass::TandemDuplication)
}

/// The allele string of the structural variant feature: the ALTs joined with `/`, and for a
/// breakpoint the reference and the breakends, with the mate that Manta's `CHR2` names written in.
fn allele_string_of(class: SvClass, alt_slash: &str, ref_given: &str, info: &str) -> String {
    let mut allele_string = alt_slash.to_string();
    if !class.is_breakpoint() {
        return allele_string;
    }
    if let Some(chr2) = info_value(info, "CHR2") {
        let end2 = info_value(info, "END2");
        // Manta puts the mate's position in END, which a breakpoint's own interval does not use
        let mate_position = if truthy(end2) {
            end2
        } else {
            info_value(info, "END")
        };
        if let Some(position) = mate_position {
            let is_bnd = {
                let t = allele_string.strip_prefix('<').unwrap_or(&allele_string);
                let t = t.strip_suffix('>').unwrap_or(t);
                t.eq_ignore_ascii_case("BND")
            };
            let base = if is_bnd {
                "N".to_string()
            } else {
                format!("{allele_string}/N")
            };
            allele_string = format!("{base}[{chr2}:{position}[");
        }
    }
    if !(allele_string.starts_with('.') || allele_string.ends_with('.')) {
        allele_string = format!("{ref_given}/{allele_string}");
    }
    allele_string
}

/// The repeat that Ensembl VEP reads as bases: `if ($info->{RUS} && ($info->{RUC} || $info->{RB}))`.
fn tandem_of(class: SvClass, info: &str, type_text: &str) -> Option<TandemRepeat> {
    // `($info->{RUC} || '')` and `($info->{RB} || '')`: a count of 0 is false in Perl, so it is no
    // count, and `RB` is what is read where `RUC=0` stands.
    let true_or_empty = |key: &str| {
        info_value(info, key)
            .filter(|value| truthy(Some(value)))
            .unwrap_or("")
            .to_string()
    };
    (is_tandem(class)
        && truthy(info_value(info, "RUS"))
        && (truthy(info_value(info, "RUC")) || truthy(info_value(info, "RB"))))
    .then(|| TandemRepeat {
        units: info_value(info, "RUS").unwrap_or("").to_string(),
        counts: true_or_empty("RUC"),
        bases: true_or_empty("RB"),
        records: info_value(info, "RN")
            .filter(|rn| !rn.is_empty() && *rn != "0")
            .map(String::from),
        types: type_text.matches("CNV:TR").count(),
    })
}

/// The alleles VEP writes rows for (`StructuralVariationOverlap::new`).
///
/// A class that is not a breakpoint has one per piece of the allele string, each named by the class
/// and numbered from 1, so a record with two ALTs has two identical rows. A breakpoint with mates
/// has the variant itself (numbered 0, named by its reference base and a dot on the side the mates
/// are not) and one per mate, named by its own text. A breakpoint with no mates (`.A`, `A.`) is
/// alone, with no position of its own.
fn alleles_of(class: SvClass, allele_string: &str, alts: &[(usize, &str)]) -> Vec<SvAllele> {
    let term: Arc<str> = Arc::from(class.so_term());
    let alt_of = |i: usize| alts.get(i).or_else(|| alts.last()).map_or(0, |a| a.0);
    let pieces: Vec<&str> = allele_string.split('/').collect();

    let plain = |pieces: &[&str]| -> Vec<SvAllele> {
        (0..pieces.len())
            .map(|i| SvAllele {
                label: Arc::clone(&term),
                site: AlleleSite::None,
                alt: alt_of(i),
            })
            .collect()
    };
    if !class.is_breakpoint() {
        return plain(&pieces);
    }

    // `_parse_breakends`
    let mut reference: Option<&str> = None;
    let mut left = false;
    let mut mates: Vec<(&str, &str, u64)> = Vec::new();
    for piece in &pieces {
        match parse_mate(piece) {
            Some((contig, position, placed_left)) => {
                left |= placed_left;
                mates.push((piece, contig, position));
            }
            None => reference = reference_of_piece(piece),
        }
    }
    if mates.is_empty() {
        return plain(&pieces);
    }
    let base = reference.unwrap_or("N");
    let own = if left {
        format!(".{base}")
    } else {
        format!("{base}.")
    };
    let mut out = vec![SvAllele {
        label: Arc::from(own),
        site: AlleleSite::Own,
        alt: alt_of(0),
    }];
    for (n, (piece, contig, position)) in mates.into_iter().enumerate() {
        out.push(SvAllele {
            label: Arc::from(piece),
            site: AlleleSite::Mate {
                chromosome: Arc::from(contig),
                position,
            },
            alt: alt_of(n),
        });
    }
    out
}

/// What a record with these ALTs is as far as the reader needs to know before it has read the
/// alleles: used by the small-variant path to leave a record to the structural one.
pub fn real_alts(alt_field: &str) -> impl Iterator<Item = (usize, &str)> + Clone {
    // `split(',')` drops the empty fields at the end, so `<DEL>,` is the one allele and `,,` none
    let alts = alt_field.trim_end_matches(',');
    alts.split(',')
        .enumerate()
        .filter(move |(_, alt)| !alts.is_empty() && !is_placeholder(alt))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn read_one(pos: u64, reference: &str, alt: &str, info: &str) -> Result<SvRecord, Skipped> {
        let alts: Vec<(usize, &str)> = real_alts(alt).collect();
        read(pos, reference, &alts, info, false)
    }

    fn read_one_as_point(pos: u64, reference: &str, alt: &str, info: &str) -> SvRecord {
        let alts: Vec<(usize, &str)> = real_alts(alt).collect();
        read(pos, reference, &alts, info, true).unwrap()
    }

    fn routed(reference: &str, alt: &str, info: &str) -> bool {
        is_structural(reference, real_alts(alt).map(|a| a.1), info)
    }

    #[test]
    fn symbolic_and_breakend_records_are_structural_and_bases_are_not() {
        assert!(routed("A", "<DEL>", "END=10"));
        assert!(routed("A", "<INS>", "."));
        assert!(routed("A", "A[chr1:2500[", "."));
        assert!(routed("A", ".A", "."));
        assert!(routed("A", "A.", "."));
        assert!(routed("A", "G,<DEL>", "."));
        for (reference, alt) in [("A", "G"), ("AT", "A"), ("A", "ACGT,T"), ("C", "G,T")] {
            assert!(!routed(reference, alt, "SVTYPE=DEL"), "{reference}>{alt}");
        }
    }

    #[test]
    fn svtype_makes_a_record_structural_when_an_allele_is_not_bases() {
        assert!(routed("A", "N", "SVTYPE=DEL"));
        assert!(!routed("A", "N", "."));
        // a flag, or a zero, is not true
        assert!(!routed("A", "N", "SVTYPE"));
        assert!(!routed("A", "N", "SVTYPE=0"));
    }

    #[test]
    fn a_star_between_brackets_is_not_a_symbolic_allele_and_the_placeholders_are_not_alleles() {
        assert!(!looks_symbolic("<*>"));
        assert!(looks_symbolic("<*>,<DEL>"));
        assert!(!routed("A", "<*>", "END=10"));
        assert!(!routed("A", "<NON_REF>", "END=10"));
        assert!(!routed("A", ".", "."));
        assert!(!routed("A", "*", "."));
        // beside a real ALT the placeholder is not an allele, so the record is its real ALT's
        assert!(!routed("C", "G,<NON_REF>", "."));
        assert!(routed("A", "<DEL>,<NON_REF>", "END=10"));
    }

    #[test]
    fn the_start_is_the_base_after_pos_for_every_symbolic_class() {
        for alt in ["<DEL>", "<DUP>", "<INV>", "<CNV>", "<INS>"] {
            let r = read_one(100, "A", alt, "END=200").unwrap();
            assert_eq!((r.start, r.end), (101, 200), "{alt}");
        }
    }

    #[test]
    fn svlen_wins_over_end_and_an_insertion_spans_its_length() {
        let r = read_one(100, "A", "<DEL>", "END=900;SVLEN=-19").unwrap();
        assert_eq!((r.start, r.end), (101, 119));
        let r = read_one(100, "A", "<INS>", "SVLEN=50").unwrap();
        assert_eq!((r.start, r.end), (101, 150));
        // `END=POS` is a zero-length insertion: start is end + 1
        let r = read_one(100, "A", "<INS>", "END=100").unwrap();
        assert_eq!((r.start, r.end), (101, 100));
        // neither: the one base after POS
        let r = read_one(100, "A", "<INS>", ".").unwrap();
        assert_eq!((r.start, r.end), (101, 101));
        // an SVLEN that is not a number is zero long
        let r = read_one(100, "A", "<INS>", "SVLEN=.").unwrap();
        assert_eq!((r.start, r.end), (101, 100));
    }

    #[test]
    fn a_deletion_without_an_extent_or_that_ends_before_it_starts_is_skipped() {
        let skipped = |info: &str| read_one(100, "A", "<DEL>", info).unwrap_err().0;
        assert_eq!(skipped("."), "deletion looks incomplete");
        assert_eq!(skipped("END=100"), "deletion looks incomplete");
        assert!(skipped("END=50").starts_with("deletion looks incomplete"));
        // a duplication with no extent is the one base
        assert!(read_one(100, "A", "<DUP>", ".").is_ok());
        // an end more than one base before the start is refused for any class
        assert!(read_one(100, "A", "<INV>", "END=50")
            .unwrap_err()
            .0
            .starts_with("start > end+1"));
    }

    #[test]
    fn a_type_vep_does_not_know_is_skipped_with_its_name() {
        assert_eq!(
            read_one(100, "A", "<CPX>", "END=200").unwrap_err().0,
            "CPX is not a supported structural variant type"
        );
        // SVTYPE stands for an ALT that does not start with a type
        assert_eq!(
            read_one(100, "A", "<X>", "SVTYPE=DEL;END=200")
                .unwrap()
                .sv
                .class,
            SvClass::Deletion
        );
        // and does not for one that does
        assert_eq!(
            read_one(100, "A", "<DEL>", "SVTYPE=DUP;END=200")
                .unwrap()
                .sv
                .class,
            SvClass::Deletion
        );
    }

    #[test]
    fn several_alts_are_one_allele_each_with_the_same_name() {
        let r = read_one(100, "A", "<DEL>,<DUP>", "END=200").unwrap();
        assert_eq!(r.sv.class, SvClass::CopyNumberVariation);
        assert_eq!(r.allele_string, "<DEL>/<DUP>");
        let labels: Vec<_> = r.sv.alleles.iter().map(|a| (&*a.label, a.alt)).collect();
        assert_eq!(
            labels,
            [("copy_number_variation", 0), ("copy_number_variation", 1)]
        );
        // a base ALT beside an SV one takes the SVTYPE and still gets a row of its own
        let r = read_one(100, "G", "C,<DEL>", "END=200;SVTYPE=DEL").unwrap();
        assert_eq!(r.sv.class, SvClass::Deletion);
        assert_eq!(r.sv.alleles.len(), 2);
    }

    #[test]
    fn a_breakend_with_a_mate_has_its_own_allele_and_one_for_the_mate() {
        let r = read_one(2020, "A", "A[chr1:2500[", "SVTYPE=BND").unwrap();
        assert_eq!((r.start, r.end), (2021, 2021));
        assert_eq!(r.allele_string, "A/A[chr1:2500[");
        let a = &r.sv.alleles;
        assert_eq!(a.len(), 2);
        assert_eq!((&*a[0].label, &a[0].site), ("A.", &AlleleSite::Own));
        assert_eq!(&*a[1].label, "A[chr1:2500[");
        assert_eq!(
            a[1].site,
            AlleleSite::Mate {
                chromosome: Arc::from("chr1"),
                position: 2500
            }
        );
        // a mate to the left puts the dot on the left
        let r = read_one(2020, "A", "]chr1:2500]A", "SVTYPE=BND").unwrap();
        assert_eq!(&*r.sv.alleles[0].label, ".A");
    }

    #[test]
    fn a_single_breakend_is_one_allele_named_by_its_class_with_no_position() {
        for alt in [".A", "A."] {
            let r = read_one(2020, "A", alt, "SVTYPE=BND").unwrap();
            assert_eq!((r.start, r.end), (2021, 2021), "{alt}");
            assert_eq!(r.allele_string, alt);
            assert_eq!(r.sv.alleles.len(), 1);
            assert_eq!(&*r.sv.alleles[0].label, "chromosome_breakpoint");
            assert_eq!(r.sv.alleles[0].site, AlleleSite::None);
        }
    }

    #[test]
    fn a_symbolic_breakend_with_chr2_is_rewritten_to_the_mate_it_names() {
        for info in [
            "SVTYPE=BND;CHR2=chr1;END2=2500",
            "SVTYPE=BND;CHR2=chr1;END=2500",
        ] {
            let r = read_one(2020, "A", "<BND>", info).unwrap();
            assert_eq!(r.allele_string, "A/N[chr1:2500[", "{info}");
            // the END that was the mate is not the end of the variant
            assert_eq!((r.start, r.end), (2021, 2021), "{info}");
            assert_eq!(&*r.sv.alleles[1].label, "N[chr1:2500[");
        }
    }

    #[test]
    fn a_mate_on_a_contig_with_colons_keeps_them() {
        assert_eq!(
            parse_mate("A[HLA-A*01:01:01:01:55["),
            Some(("HLA-A*01:01:01:01", 55, false))
        );
        assert_eq!(parse_mate("]2:321682]A"), Some(("2", 321682, true)));
        assert_eq!(parse_mate("A"), None);
        assert_eq!(parse_mate(".A"), None);
    }

    #[test]
    fn the_alt_index_of_each_allele_is_that_of_the_alt_in_the_column() {
        // a placeholder before the SV does not shift the index the allele answers for
        let alts: Vec<(usize, &str)> = real_alts("<NON_REF>,<DEL>").collect();
        let r = read(100, "A", &alts, "END=200", false).unwrap();
        assert_eq!(r.sv.alleles[0].alt, 1);
    }

    #[test]
    fn an_info_piece_is_read_as_perl_splits_it() {
        for (piece, expected) in [
            ("A=b", ("A", Some("b"))),
            ("A=b=c", ("A", Some("b"))),
            ("A=", ("A", Some(""))),
            ("A==", ("A", Some(""))),
            ("A=b=", ("A", Some("b"))),
            ("A==c", ("A", Some(""))),
            ("A", ("A", None)),
            ("=v", ("", Some("v"))),
            ("", ("", None)),
        ] {
            assert_eq!(key_value(piece), expected, "{piece:?}");
        }
    }

    #[test]
    fn a_piece_of_a_million_brackets_is_read_once() {
        // each bracket was a place to start from, and each start read the rest of the piece
        let piece = format!("G{}", "[".repeat(1_000_000));
        let started = std::time::Instant::now();
        assert_eq!(parse_mate(&piece), None);
        let piece = format!("G{}chr1:100[", "[".repeat(1_000_000));
        assert_eq!(
            parse_mate(&piece).map(|(contig, position, left)| (contig.len(), position, left)),
            Some((999_999 + 4, 100, false))
        );
        assert!(started.elapsed().as_secs() < 5);
    }

    #[test]
    fn a_trailing_comma_is_not_another_allele() {
        let alts: Vec<_> = real_alts("<DEL>,").collect();
        assert_eq!(alts, vec![(0, "<DEL>")]);
        assert_eq!(real_alts(",,").count(), 0);
        assert_eq!(real_alts("<DEL>,<DUP>").count(), 2);
        // a gap in the middle is an allele, as it is for Perl
        assert_eq!(real_alts("<DEL>,,<DUP>").count(), 3);
    }

    #[test]
    fn an_end_of_zero_on_a_breakend_is_the_end() {
        // Perl takes the END out of a breakend's INFO only when it is true, and `0` is not
        let r = read_one(100, "A", "A[chr1:2500[", "SVTYPE=BND;END=0");
        assert_eq!(
            r.unwrap_err().0,
            "start > end+1 : (START=101, END=0)".to_string()
        );
        let r = read_one(100, "A", "A[chr1:2500[", "SVTYPE=BND;END=2500").unwrap();
        assert_eq!((r.start, r.end), (101, 101));
    }

    #[test]
    fn lengths_and_positions_beyond_the_range_of_a_coordinate_do_not_wrap() {
        // the largest a length can be, which `abs` and the sum over a start made wrap
        for svlen in ["SVLEN=-9223372036854775807", "SVLEN=9223372036854775807"] {
            let r = read_one(100, "A", "<DUP>", svlen).unwrap();
            assert_eq!(r.start, 101);
            assert!(r.end >= r.start && r.end <= i64::MAX as u64, "{svlen}");
        }
        // a position beyond what `i64` holds is a coordinate no transcript has, and not a negative one
        let r = read_one(i64::MAX as u64, "A", "<DUP>", "SVLEN=5").unwrap();
        assert!(r.start <= i64::MAX as u64 && r.end <= i64::MAX as u64);
        let r = read_one(u64::MAX - 1, "A", "<INS>", ".").unwrap();
        assert_eq!((r.start, r.end), (i64::MAX as u64, i64::MAX as u64));
    }

    /// `KEY=` is a value, the empty string (`key_value`), so an empty `SVLEN` is a length of nothing
    /// and an empty `END` is no coordinate, as `get_end` reads them (`defined($info->{SVLEN})`).
    #[test]
    fn an_empty_svlen_or_end_is_a_value_and_not_an_absent_one() {
        let skipped =
            |pos, reference, alt, info| read_one(pos, reference, alt, info).unwrap_err().0;
        assert_eq!(
            skipped(1060, "C", "<DEL>", "SVLEN=;END=1085"),
            "deletion looks incomplete"
        );
        assert!(skipped(1061, "T", "<DUP>", "END=").starts_with("Invalid start"));
        // an inversion with the empty length is the interval of no bases after its padding
        let r = read_one(1062, "T", "<INV>", "SVLEN=;END=1090").unwrap();
        assert_eq!((r.start, r.end), (1063, 1062));
    }

    #[test]
    fn a_length_past_the_largest_number_is_that_number() {
        assert_eq!(perl_integer("99999999999999999999"), i64::MAX);
        assert_eq!(perl_integer("-99999999999999999999"), -i64::MAX);
        assert_eq!(perl_integer(" +12abc"), 12);
        assert_eq!(perl_integer("."), 0);
        // so the record is as large as any, which the size limit refuses, and not an insertion of
        // no bases that an annotation then writes rows for
        let r = read_one(1060, "C", "<INS>", "SVLEN=99999999999999999999").unwrap();
        assert!(r.end > 10_000_000_000, "{}", r.end);
    }

    /// `RUC=0` is false in Perl (`$info->{RUC} || ''`): it is no count, and `RB` is read in its place.
    #[test]
    fn a_count_of_zero_is_no_count_and_the_bases_are_read() {
        let r = read_one(
            3399,
            "C",
            "<CNV:TR>",
            "SVTYPE=CNV;END=3439;SVLEN=40;RN=1;RUS=CA;RUC=0;RB=10",
        )
        .unwrap();
        let tandem = r.sv.tandem.expect("a tandem repeat");
        assert_eq!((tandem.counts.as_str(), tandem.bases.as_str()), ("", "10"));
        assert_eq!(
            tandem.alleles(Some(1000)),
            Some(vec!["CACACACACA".to_string()])
        );
        // and with neither a count nor the bases it is not read as bases at all
        let r = read_one(
            3399,
            "C",
            "<CNV:TR>",
            "SVTYPE=CNV;END=3439;SVLEN=40;RN=1;RUS=CA;RUC=0",
        )
        .unwrap();
        assert!(r.sv.tandem.is_none());
    }

    #[test]
    fn an_insertion_is_the_point_between_pos_and_the_next_base_when_it_is_asked_to_be() {
        // as VEP reads it: the bases after `POS`, as many as `SVLEN` says
        let r = read_one(100, "A", "<INS:ME:ALU>", "SVLEN=300").unwrap();
        assert_eq!((r.start, r.end), (101, 400));
        // as a point: `start = end + 1`, whatever the length and the `END` are
        for info in ["SVLEN=300", "END=500;SVLEN=300", "END=500", "."] {
            let r = read_one_as_point(100, "A", "<INS:ME:ALU>", info);
            assert_eq!((r.start, r.end), (101, 100), "{info}");
        }
        for alt in ["<INS>", "<INS:ME>", "<INS:ME:LINE1>", "<INS:ME:SVA>"] {
            let r = read_one_as_point(100, "A", alt, "SVLEN=40");
            assert_eq!((r.start, r.end), (101, 100), "{alt}");
        }
        // and nothing else is an insertion: a duplication, a deletion and a tandem repeat keep their
        // intervals
        for (alt, info) in [
            ("<DUP>", "END=200"),
            ("<DEL>", "END=200"),
            ("<DUP:TANDEM>", "END=200"),
            ("<INV>", "END=200"),
        ] {
            let r = read_one_as_point(100, "A", alt, info);
            assert_eq!((r.start, r.end), (101, 200), "{alt}");
        }
    }
}
