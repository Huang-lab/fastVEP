//! HGVS protein nomenclature, as Ensembl VEP writes it.
//!
//! This is a transcription of `hgvs_protein` and the routines it calls in
//! `Bio::EnsEMBL::Variation::TranscriptVariationAllele` (release 116), one Rust function to each
//! Perl sub and in the same order, so that a difference between the two can be found by reading
//! them side by side. Ensembl's notation is not derived from the protein the change produces;
//! it is the result of a sequence of string operations on two short peptides (the window the
//! change edits, before and after), and what it writes at the ends of a protein - a start codon,
//! a stop codon, the last residue - is the sum of those operations and nothing else. Rules
//! inferred from its outputs reproduced the common cases and missed those, which is what this
//! replaces.
//!
//! What the Perl reads from the transcript and the variant is in [`ProteinChange`] and
//! [`EditedTranslation`]. Both are *inputs*: this module knows no coordinates, no strand and no
//! genome, which is what makes it testable against Ensembl's own rows one at a time.
//!
//! Four places differ from Ensembl on purpose. Each is marked where it happens, and each is a case
//! where what Ensembl writes is not what the sequence says or is not HGVS:
//!
//! 1. **`*` and `X` stay different.** Ensembl rewrites `*` to `X` before it tests the
//!    peptides and then rewrites every `Xaa` to `Ter`, so an unresolved residue is written as
//!    a terminator: a change whose alternate window ends on an incomplete codon is
//!    `p.Ter201delinsTyrTer` there and a deletion in a stop codon is `p.Ter583=` next to a
//!    `stop_lost` consequence. Here `*` is `Ter` and `X` is `Xaa`, and a change that moves
//!    the frame at a stop codon, where the window ends mid-codon, is handed to this module as a
//!    frameshift so that it is described from the edited translation (the driver,
//!    `fastvep_annotate::hgvsp_for_allele`, decides).
//! 2. **A shift reaches the end of the protein.** Ensembl stops short of it in two ways. The
//!    residues after a change come from `_get_surrounding_peptides`, which returns nothing when
//!    `length <= position`, so a change that ends one residue before the end is not shifted
//!    at all, and the walk in `_shift_3prime` is bounded by `length(post_seq) - n`, which
//!    halts it `n - 1` residues short for a change of `n` residues. That puts the NT5C2
//!    deletion at `p.Glu559_Glu560del` on a protein that ends in a run of seven Glu. The
//!    3'-rule has no such bound (#94).
//! 3. **A synonymous window names its span.** `p.LeuIle2672=` is two residues sharing one
//!    position, which HGVS does not define; `p.Leu2672_Ile2673=` is what it means.
//! 4. **A change that loses the initiator is `p.Met1?`.** Ensembl writes that for most of them,
//!    and for a window of more than one residue two residue names sharing one position
//!    (`p.MetAsn1_?2`), and for an insertion in front of the initiator a range that runs
//!    backwards and names the last residue of the protein as if it were the first
//!    (`p.Asp1_?0`), because `substr` reads the residue in front of residue 1 from the end of
//!    the string. Neither is HGVS and the second is wrong.

use fastvep_genome::codon::aa_one_to_three;
use std::fmt::Write;

/// What Ensembl's `hgvs_protein` reads from a `TranscriptVariationAllele` once its 3' shift has
/// been applied.
///
/// `translation_start` and `translation_end` are the variant's own protein coordinates
/// (`translation_start`, `translation_end` of the `TranscriptVariation`), of the 3'-shifted
/// variant for an insertion or a deletion; they are not the span of the window, although the
/// window is every codon between them. For an insertion between two codons `translation_end`
/// is one less than `translation_start`.
#[derive(Debug, Clone, Copy)]
pub struct ProteinChange<'a> {
    /// `ENSP00000000003.1`: the translation's name with its version.
    pub ref_name: &'a str,
    pub translation_start: u64,
    pub translation_end: u64,
    /// `peptide` of the reference allele and of this one, over the window. `-` for a window with
    /// no residues, which is what Ensembl's `peptide` returns for one. `*` is a terminator and
    /// `X` a residue that could not be resolved.
    pub ref_peptide: &'a str,
    pub alt_peptide: &'a str,
    /// The reference CDS translated: the protein, and a final `*` when the CDS ends in a stop
    /// codon (`_translateable_seq` through `translate`). `_peptide` is this without that `*`.
    pub reference_translation: &'a str,
    /// `frameshift`, `start_lost` and `stop_lost` of `VariationEffect.pm`, as predicates and not as
    /// terms: `hgvs_protein` calls them directly, so the `include` filters that decide which
    /// terms Ensembl asks about do not apply.
    pub frameshift: bool,
    pub start_lost: bool,
    pub stop_lost: bool,
}

impl ProteinChange<'_> {
    /// `_peptide`: the protein without its terminator.
    fn protein(&self) -> &str {
        self.reference_translation
            .strip_suffix('*')
            .unwrap_or(self.reference_translation)
    }
}

/// The part of Ensembl's `hgvs_protein` that needs the edited transcript, which is the
/// translation of `_get_alternate_cds`: the CDS with the change applied and the 3' UTR appended.
///
/// A trait and not a slice because most notations never ask for it - a substitution is
/// decided by the two residues in the window - and building it costs an allocation and a
/// translation of everything downstream of the change. `None` is a CDS Ensembl cannot edit
/// (`_get_alternate_cds` returns undef), and the notation is then not written.
pub trait EditedTranslation {
    fn alt_translation(&self) -> Option<&str>;
}

/// Perl's `substr($s, $offset, $length)` over ASCII.
///
/// A negative offset counts from the end, a negative length leaves that many off it, and a
/// range partly outside the string is clipped to the part inside. A read that starts past the
/// end is undef in Perl; every use here treats that as the empty string, which is how Perl's
/// string operators treat it.
fn substr(s: &str, offset: i64, length: Option<i64>) -> &str {
    let n = s.len() as i64;
    let first = if offset < 0 { n + offset } else { offset };
    if first > n {
        return "";
    }
    let last = match length {
        None => n,
        Some(l) if l >= 0 => first + l,
        Some(l) => n + l,
    };
    let (lo, hi) = (first.clamp(0, n), last.clamp(0, n));
    if hi <= lo {
        return "";
    }
    &s[lo as usize..hi as usize]
}

/// `Bio::SeqUtils::seq3`.
fn seq3(peptide: &str) -> String {
    peptide.bytes().map(aa_one_to_three).collect()
}

/// Ensembl's `$hgvs_notation` hash.
#[derive(Debug, Clone, Default)]
struct Notation {
    start: i64,
    end: i64,
    reference: String,
    alt: String,
    /// `_clip_alleles` keeps the reference as it was before clipping, for the routines that
    /// look past the end of the protein for a terminator.
    original_ref: Option<String>,
    /// What clipping took off the front.
    preseq: Option<String>,
    /// `$hgvs_notation->{type}`; `None` is Perl's undef.
    kind: Option<Kind>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Kind {
    Fs,
    Ins,
    Del,
    Sub,
    Delins,
    Dup,
    Same,
    /// The empty string, which `_get_hgvs_peptides` assigns to a change that loses the initiator.
    Unknown,
}

/// `hgvs_protein`: the HGVSp of one allele, `None` where Ensembl writes nothing.
pub fn hgvs_protein(change: &ProteinChange<'_>, edited: &dyn EditedTranslation) -> Option<String> {
    // The shape nearly every row has: one residue for another, or one that is not changed, and
    // nothing about either that the rest of the routine would look at. What it does with such a
    // pair is `Arg41Lys` or `Leu41=`, without clipping or shifting or a second look at the
    // protein, and writing that directly is a `String` where the general path is a dozen.
    if let ([r], [a]) = (change.ref_peptide.as_bytes(), change.alt_peptide.as_bytes()) {
        let plain = |b: &u8| b.is_ascii_uppercase() && *b != b'X';
        if !(change.frameshift || change.start_lost || change.stop_lost) && plain(r) && plain(a) {
            // one allocation of the size it needs, which `format!` does not promise
            let mut out = String::with_capacity(change.ref_name.len() + 20);
            out.push_str(change.ref_name);
            out.push_str(":p.");
            out.push_str(aa_one_to_three(*r));
            let _ = write!(out, "{}", change.translation_start);
            out.push_str(if r == a { "=" } else { aa_one_to_three(*a) });
            return Some(out);
        }
    }
    let mut n = Notation {
        start: change.translation_start as i64,
        end: change.translation_end as i64,
        reference: change.ref_peptide.to_string(),
        alt: change.alt_peptide.to_string(),
        ..Notation::default()
    };
    // `return undef unless $hgvs_notation->{ref}`
    if n.reference.is_empty() {
        return None;
    }
    if n.alt != n.reference {
        clip_alleles(&mut n);
    }
    get_hgvs_protein_type(change, &mut n);
    get_hgvs_peptides(change, edited, &mut n)?;
    Some(get_hgvs_protein_format(change, edited, n))
}

/// `_clip_alleles` for protein numbering: take what the two peptides share off both ends.
fn clip_alleles(n: &mut Notation) {
    let (mut check_alt, mut check_ref) = (n.alt.clone(), n.reference.clone());
    let (mut check_start, mut check_end) = (n.start, n.end);
    n.original_ref = Some(n.reference.clone());

    let mut preseq = String::new();
    // The bound is the length of the *original* reference, as in the Perl; when the two run out
    // together every remaining pass compares two empty strings and moves the start on.
    for _ in 0..n.reference.len() {
        let next_ref = substr(&check_ref, 0, Some(1)).to_string();
        let next_alt = substr(&check_alt, 0, Some(1)).to_string();

        // a stop re-created by the change: no change to the protein
        if next_ref == "*" && next_alt == "*" {
            n.kind = Some(Kind::Same);
            return;
        }
        if next_ref == next_alt {
            check_start += 1;
            check_ref = substr(&check_ref, 1, None).to_string();
            check_alt = substr(&check_alt, 1, None).to_string();
            preseq.push_str(&next_ref);
        } else {
            break;
        }
    }

    let len = check_ref.len();
    for _ in 0..len {
        let next_ref = substr(&check_ref, -1, Some(1)).to_string();
        let next_alt = substr(&check_alt, -1, Some(1)).to_string();
        if next_ref == next_alt {
            check_ref.pop();
            check_alt.pop();
            check_end -= 1;
        } else {
            break;
        }
    }

    n.alt = check_alt;
    n.reference = check_ref;
    n.start = check_start;
    n.end = check_end;
    n.preseq = Some(preseq);

    // The type set here is not read again - `_get_hgvs_protein_type` recomputes it from the
    // clipped peptides - but it is kept because the Perl sets it.
    if n.reference == n.alt {
        n.kind = Some(Kind::Same);
    } else if n.reference != "-" && n.reference.len() == 1 && n.alt.len() == 1 {
        n.kind = Some(Kind::Sub);
    } else if n.reference.is_empty() && !n.alt.is_empty() {
        n.kind = Some(Kind::Ins);
    } else if !n.reference.is_empty() && n.alt.is_empty() {
        n.kind = Some(Kind::Del);
    }
}

/// `_get_hgvs_protein_type`, for a window whose two peptides are both known (the only one this
/// is called with: Ensembl's other branch is for an allele with no peptide, and none is
/// written for one).
fn get_hgvs_protein_type(change: &ProteinChange<'_>, n: &mut Notation) {
    if change.frameshift {
        n.kind = Some(Kind::Fs);
        return;
    }
    // Ensembl rewrites the first `*` of each peptide to `X` here, which is how a terminator
    // becomes a residue, and so the same string as an unresolved one later on. Only the
    // lengths are read below, so leaving `*` alone changes nothing but that (#1 of the module
    // doc).
    let (r, a) = (n.reference.len(), n.alt.len());
    n.kind = Some(if n.reference == "-" || n.reference.is_empty() {
        Kind::Ins
    } else if n.alt.is_empty() || n.alt == "-" {
        Kind::Del
    } else if r == 1 && a == 1 {
        Kind::Sub
    } else if (a > 0 && r > 0 && a != r) || (a > 1 && r > 1) {
        Kind::Delins
    } else {
        Kind::Sub
    });
}

/// `_get_hgvs_peptides`: the reference and alternate peptides as the notation names them, in
/// three-letter code. `None` where Ensembl returns undef.
fn get_hgvs_peptides(
    change: &ProteinChange<'_>,
    edited: &dyn EditedTranslation,
    n: &mut Notation,
) -> Option<()> {
    match n.kind {
        Some(Kind::Fs) => {
            // the Ensembl alt and ref peptides are not the HGVS ones for a frameshift
            get_fs_peptides(change, edited, n)?;
        }
        Some(Kind::Ins) => {
            // is what follows the insertion a repeat of it?
            check_peptides_post_var(change, n);
            // and does it repeat what precedes it?
            if !n.alt.contains('*') {
                check_for_peptide_duplication(change, n);
            }
            if n.kind == Some(Kind::Dup) {
                return Some(());
            }
            // the reference of an insertion is the pair of residues around it
            let min = n.start.min(n.end);
            let flank = get_surrounding_peptides(change, min, n.original_ref.as_deref(), Some(2))?;
            if flank.is_empty() {
                return None;
            }
            n.reference = flank;
        }
        Some(Kind::Del) => check_peptides_post_var(change, n),
        _ => {}
    }

    if n.reference != "-" {
        n.reference = seq3(&n.reference);
    }
    if n.alt != "-" {
        n.alt = seq3(&n.alt);
    }
    if n.alt == "-" {
        n.alt = "del".to_string();
    }

    if change.start_lost {
        // the initiator is gone, so there is probably no translation at all: the alternate
        // allele is unknown, and what was lost is residue 1 whatever the window is (#4 of the
        // module doc)
        let first = substr(change.protein(), 0, Some(1));
        n.reference = if first.is_empty() {
            "Met".to_string()
        } else {
            seq3(first)
        };
        n.start = 1;
        n.end = 1;
        n.alt = "?".to_string();
        n.kind = Some(Kind::Unknown);
    } else if n.kind == Some(Kind::Del) {
        // `if ($hgvs_notation->{ref} =~ /\w+/) { alt = "del" } else { _get_del_peptides }`. A
        // deletion has a reference that is neither empty nor `-` (that is what makes it a
        // deletion), so the second arm cannot be reached from here.
        n.alt = "del".to_string();
    } else if n.kind == Some(Kind::Fs) {
        // only the first reference residue is quoted for a frameshift
        n.reference = substr(&n.reference, 0, Some(3)).to_string();
    }
    Some(())
}

/// `_get_fs_peptides`: the first residue the frameshift changes, found by comparing the
/// reference protein with the translation of the edited CDS from the variant's own start.
fn get_fs_peptides(
    change: &ProteinChange<'_>,
    edited: &dyn EditedTranslation,
    n: &mut Notation,
) -> Option<()> {
    let alt_trans = edited.alt_translation()?;
    if alt_trans.is_empty() {
        return None;
    }
    // the reference protein with its terminator, for checking purposes
    let ref_trans = format!("{}*", change.protein());
    n.start = change.translation_start as i64;
    let alt_len = alt_trans.len() as i64;

    // the terminator is deleted and nothing follows it in the alternate sequence
    if n.start > alt_len {
        n.alt = "del".to_string();
        n.kind = Some(Kind::Del);
        return Some(());
    }
    while n.start <= alt_len {
        n.reference = substr(&ref_trans, n.start - 1, Some(1)).to_string();
        n.alt = substr(alt_trans, n.start - 1, Some(1)).to_string();
        if n.reference == "*" && n.alt == "*" {
            // the change is at the stop codon but keeps it: synonymous
            n.kind = Some(Kind::Same);
            return Some(());
        }
        // a frameshift may leave the residue unchanged for a while
        if n.reference != n.alt {
            break;
        }
        n.start += 1;
    }
    Some(())
}

/// `_get_surrounding_peptides`: the residues of the reference protein from `ref_pos`, or `None`
/// where there are none there.
fn get_surrounding_peptides(
    change: &ProteinChange<'_>,
    ref_pos: i64,
    original_ref: Option<&str>,
    length: Option<i64>,
) -> Option<String> {
    let mut ref_trans = change.protein().to_string();
    // a change that starts on the terminator has one to read
    if let Some(original) = original_ref.filter(|o| o.starts_with('*')) {
        ref_trans.push_str(original);
    }
    // can't find a residue after the end
    if (ref_trans.len() as i64) <= ref_pos {
        return None;
    }
    Some(substr(&ref_trans, ref_pos - 1, length).to_string())
}

/// `_check_peptides_post_var`: HGVS counts from the first residue that differs, so look at what
/// follows the change and move it along while it repeats what the change is made of.
fn check_peptides_post_var(change: &ProteinChange<'_>, n: &mut Notation) {
    let post_pos = n.end + 1;
    // The residues from `post_pos` on. Ensembl asks `_get_surrounding_peptides`, which answers
    // nothing when `length <= post_pos`, so a change that ends one residue short of the end of
    // the protein has the last residue hidden from it and is not shifted at all (#2 of the module
    // doc): that function is written for the residues around an insertion, and compares
    // like one. The last residue is a residue.
    let mut protein = change.protein().to_string();
    if let Some(original) = n.original_ref.as_deref().filter(|o| o.starts_with('*')) {
        protein.push_str(original);
    }
    if post_pos < 1 || post_pos > protein.len() as i64 {
        // a stop is deleted and there is no sequence beyond it to check
        return;
    }
    shift_3prime(n, substr(&protein, post_pos - 1, None));
}

/// `_shift_3prime`: HGVS aligns a change as far 3' as the sequence allows,
/// `AAC[TG]TAT => AACT[GT]AT`.
fn shift_3prime(n: &mut Notation, post_seq: &str) {
    let mut seq_to_check = match n.kind {
        Some(Kind::Ins) => n.alt.clone(),
        Some(Kind::Del) => n.reference.clone(),
        _ => return,
    };
    if seq_to_check.is_empty() {
        return;
    }
    // Ensembl's bound is `length($post_seq) - length($seq_to_check)`, which halts the walk n - 1
    // residues before the end of the protein for a change of n residues while comparing a single
    // residue per step (#2 of the module doc; #94).
    for next in post_seq.bytes() {
        let first = seq_to_check.as_bytes()[0];
        if first != next {
            break;
        }
        // move the change along, and rotate what it is made of
        n.start += 1;
        n.end += 1;
        seq_to_check.remove(0);
        seq_to_check.push(first as char);
    }
    match n.kind {
        Some(Kind::Ins) => n.alt = seq_to_check,
        _ => n.reference = seq_to_check,
    }
}

/// `_check_for_peptide_duplication`: an insertion that repeats the residues just before it is a
/// duplication.
fn check_for_peptide_duplication(change: &ProteinChange<'_>, n: &mut Notation) {
    let reference_trans = change.reference_translation;
    let mut upstream = substr(reference_trans, 0, Some(n.start - 1)).to_string();
    // add back on anything previously chopped off the reference
    if let Some(pre) = &n.preseq {
        upstream.push_str(pre);
    }

    // does the alternate match the reference just before the change?
    let alt_len = n.alt.len() as i64;
    let test_new_start = n.start - alt_len - 1;
    if upstream.len() as i64 >= test_new_start + alt_len
        && test_new_start >= 0
        && substr(&upstream, test_new_start, Some(alt_len)) == n.alt
    {
        n.kind = Some(Kind::Dup);
        n.end = n.start - 1;
        n.start -= alt_len;
        n.alt = seq3(&n.alt);
    }
}

/// `_stop_loss_extra_AA`: how many residues the protein runs on for, counted to the first stop of
/// the edited translation. `None` where there is none to count to.
///
/// `test` is `"fs"` for a frameshift, which counts from the first changed residue; anything
/// else counts from the end of the reference protein.
fn stop_loss_extra_aa(
    change: &ProteinChange<'_>,
    edited: &dyn EditedTranslation,
    ref_var_pos: i64,
    test: Option<&str>,
) -> Option<i64> {
    // `return undef unless $ref_var_pos`
    if ref_var_pos == 0 {
        return None;
    }
    let alt_trans = edited.alt_translation()?;
    // `$+[0]`: one past the index of the first stop, which is its 1-based position
    let first_stop = alt_trans.find('*')? as i64 + 1;
    let extra = if test == Some("fs") {
        first_stop - ref_var_pos
    } else {
        first_stop - 1 - change.protein().len() as i64
    };
    // a stop in the first place counted is not an extension
    (extra > 0).then_some(extra)
}

/// `_get_hgvs_protein_format`.
fn get_hgvs_protein_format(
    change: &ProteinChange<'_>,
    edited: &dyn EditedTranslation,
    mut n: Notation,
) -> String {
    let mut hgvs = format!("{}:p.", change.ref_name);
    let kind = n.kind;

    // The synonymous way to write it since release 15.11. A window of more than one residue is
    // written as the span it is (#3 of the module doc).
    if n.reference == n.alt && kind != Some(Kind::Fs) && kind != Some(Kind::Ins) {
        if n.reference.len() > 3 && n.reference.len().is_multiple_of(3) {
            let first = substr(&n.reference, 0, Some(3));
            let last = substr(&n.reference, -3, Some(3));
            let lo = n.start.min(n.end);
            let hi = n.start.max(n.end);
            hgvs.push_str(&format!("{first}{lo}_{last}{hi}="));
        } else {
            hgvs.push_str(&format!("{}{}=", n.reference, n.start));
        }
        return hgvs;
    }

    // a stop lost to a deletion or a substitution, whatever the cause: the residue it became and
    // how far the protein now runs
    if change.stop_lost && matches!(kind, Some(Kind::Del) | Some(Kind::Sub)) {
        n.alt = substr(&n.alt, 0, Some(3)).to_string();
        // `?` for a new stop that cannot be found
        let aa_til_stop = stop_loss_extra_aa(change, edited, n.start - 1, None)
            .map_or_else(|| "?".to_string(), |v| v.to_string());
        n.alt.push_str(&format!("extTer{aa_til_stop}"));

        if n.reference.len() > 3 && kind == Some(Kind::Del) {
            let first = substr(&n.reference, 0, Some(3));
            let last = substr(&n.reference, -3, Some(3));
            hgvs.push_str(&format!("{first}{}_{last}{}{}", n.start, n.end, n.alt));
        } else {
            hgvs.push_str(&format!("{}{}{}", n.reference, n.start, n.alt));
        }
    } else if kind == Some(Kind::Dup) {
        if n.start < n.end {
            // only the first and last of a long duplicated run
            let first = substr(&n.alt, 0, Some(3));
            let last = substr(&n.alt, -3, Some(3));
            hgvs.push_str(&format!("{first}{}_{last}{}dup", n.start, n.end));
        } else {
            hgvs.push_str(&format!("{}{}dup", n.alt, n.start));
        }
    } else if kind == Some(Kind::Sub) {
        hgvs.push_str(&format!("{}{}{}", n.reference, n.start, n.alt));
    } else if matches!(kind, Some(Kind::Delins) | Some(Kind::Ins)) {
        // don't report other peptides after a stop is found
        if let Some(at) = n.alt.find("Ter").filter(|&at| at + 3 < n.alt.len()) {
            n.alt.truncate(at + 3);
        }
        // the first and last residue of the reference only
        let first = substr(&n.reference, 0, Some(3)).to_string();
        let last = substr(&n.reference, -3, Some(3)).to_string();
        let type_name = if kind == Some(Kind::Ins) {
            "ins"
        } else {
            "delins"
        };
        if n.start == n.end && kind == Some(Kind::Delins) {
            hgvs.push_str(&format!("{first}{}{type_name}{}", n.start, n.alt));
        } else {
            // put the ends in order if they are not
            if n.start > n.end {
                std::mem::swap(&mut n.start, &mut n.end);
            }
            hgvs.push_str(&format!(
                "{first}{}_{last}{}{type_name}{}",
                n.start, n.end, n.alt
            ));
        }
    } else if kind == Some(Kind::Fs) {
        if n.alt == "Ter" {
            // a stop gained: a substitution, since the stop comes at once
            hgvs.push_str(&format!("{}{}{}", n.reference, n.start, n.alt));
        } else {
            // count the residues to the next stop; `?` shows that none is predicted
            let aa_til_stop = stop_loss_extra_aa(change, edited, n.start - 1, Some("fs"));
            if n.reference == "Ter" {
                // the stop is the first residue lost, which is an extension, and an extension
                // does not count the stop itself
                let count = aa_til_stop.map_or_else(|| "?".to_string(), |v| (v - 1).to_string());
                hgvs.push_str(&format!("{}{}{}extTer{count}", n.reference, n.start, n.alt));
            } else {
                let count = aa_til_stop.map_or_else(|| "?".to_string(), |v| v.to_string());
                hgvs.push_str(&format!("{}{}{}fsTer{count}", n.reference, n.start, n.alt));
            }
        }
    } else if kind == Some(Kind::Del) {
        n.alt = "del".to_string();
        if n.reference.len() > 3 {
            let first = substr(&n.reference, 0, Some(3));
            let last = substr(&n.reference, -3, Some(3));
            hgvs.push_str(&format!("{first}{}_{last}{}{}", n.start, n.end, n.alt));
        } else {
            hgvs.push_str(&format!("{}{}{}", n.reference, n.start, n.alt));
        }
    } else {
        // what is left is a change that loses the initiator: `p.Met1?`
        hgvs.push_str(&format!("{}{}{}", n.reference, n.start, n.alt));
    }
    hgvs
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The edited translation of a test, when it needs one.
    struct Edited(Option<&'static str>);

    impl EditedTranslation for Edited {
        fn alt_translation(&self) -> Option<&str> {
            self.0
        }
    }

    const NONE: Edited = Edited(None);

    /// A change to residues `start..=end` of a protein, as Ensembl sees it: the two peptides of
    /// the window and the reference translation, `*` and all.
    fn change<'a>(
        ref_peptide: &'a str,
        alt_peptide: &'a str,
        start: u64,
        end: u64,
        reference_translation: &'a str,
    ) -> ProteinChange<'a> {
        ProteinChange {
            ref_name: "P.1",
            translation_start: start,
            translation_end: end,
            ref_peptide,
            alt_peptide,
            reference_translation,
            frameshift: false,
            start_lost: false,
            stop_lost: false,
        }
    }

    fn write(change: &ProteinChange<'_>, edited: &Edited) -> String {
        hgvs_protein(change, edited).unwrap_or_else(|| "-".to_string())
    }

    #[test]
    fn perl_substr_reads_from_the_end_and_clips() {
        assert_eq!(substr("MAKL", 1, Some(2)), "AK");
        assert_eq!(substr("MAKL", -1, Some(2)), "L");
        assert_eq!(substr("MAKL", -3, Some(3)), "AKL");
        assert_eq!(substr("MAK", -3, Some(3)), "MAK");
        assert_eq!(substr("MK", -3, Some(3)), "MK");
        assert_eq!(substr("MAKL", 2, None), "KL");
        assert_eq!(substr("MAKL", 4, None), "");
        assert_eq!(substr("MAKL", 9, Some(1)), "");
        assert_eq!(substr("MAKL", 0, Some(-1)), "MAK");
    }

    #[test]
    fn a_substitution_names_both_residues() {
        let c = change("R", "Q", 41, 41, "MAKRPL*");
        assert_eq!(write(&c, &NONE), "P.1:p.Arg41Gln");
    }

    #[test]
    fn a_synonymous_change_is_an_equals_sign() {
        let c = change("L", "L", 41, 41, "MAKRPL*");
        assert_eq!(write(&c, &NONE), "P.1:p.Leu41=");
    }

    /// Ensembl writes `p.LeuIle2672=`, two residues sharing one position. The span is what it
    /// means (#3 of the module doc).
    #[test]
    fn a_synonymous_window_of_two_residues_names_its_span() {
        let c = change("LI", "LI", 2672, 2673, "MAKRPL*");
        assert_eq!(write(&c, &NONE), "P.1:p.Leu2672_Ile2673=");
    }

    #[test]
    fn a_stop_gained_is_ter() {
        let c = change("R", "*", 41, 41, "MAKRPL*");
        assert_eq!(write(&c, &NONE), "P.1:p.Arg41Ter");
    }

    /// An unresolved residue is `Xaa` and not `Ter` (#1 of the module doc).
    #[test]
    fn an_unresolved_residue_is_not_a_terminator() {
        let c = change("R", "X", 41, 41, "MAKRPL*");
        assert_eq!(write(&c, &NONE), "P.1:p.Arg41Xaa");
    }

    /// The terminator becomes a residue and the protein runs on to the next stop: the residues
    /// after the old end are counted, the one that replaces it included.
    #[test]
    fn a_lost_terminator_names_the_extension_it_causes() {
        let mut c = change("*", "Q", 4, 4, "MKP*");
        c.stop_lost = true;
        // MKP, then Q, A, A, and a stop in 7th place: 3 residues past the old end of MKP
        let edited = Edited(Some("MKPQAA*"));
        assert_eq!(write(&c, &edited), "P.1:p.Ter4GlnextTer3");
        // no stop in what follows
        assert_eq!(write(&c, &Edited(Some("MKPQAAL"))), "P.1:p.Ter4GlnextTer?");
    }

    #[test]
    fn a_deleted_terminator_is_a_deletion_that_extends() {
        let mut c = change("*", "-", 4, 4, "MAK*");
        c.stop_lost = true;
        assert_eq!(
            write(&c, &Edited(Some("MAKLLLLL*"))),
            "P.1:p.Ter4delextTer5"
        );
    }

    #[test]
    fn a_frameshift_names_the_first_residue_that_changes_and_counts_to_the_stop() {
        let mut c = change("K", "X", 3, 3, "MAKPLQRE*");
        c.frameshift = true;
        // the frame now reads M A Q * : Lys3 is the first residue that differs, and the stop is
        // the second residue of the new frame
        assert_eq!(write(&c, &Edited(Some("MAQ*"))), "P.1:p.Lys3GlnfsTer2");
        // the first residue that differs is a stop: that is a nonsense change and not a frameshift
        assert_eq!(write(&c, &Edited(Some("MA*"))), "P.1:p.Lys3Ter");
        // no stop to count to
        assert_eq!(write(&c, &Edited(Some("MAQLLL"))), "P.1:p.Lys3GlnfsTer?");
        // a residue the frameshift leaves alone is skipped over
        assert_eq!(write(&c, &Edited(Some("MAKLLQ*"))), "P.1:p.Pro4LeufsTer4");
        // an alternate CDS Ensembl cannot build is no notation
        assert_eq!(hgvs_protein(&c, &NONE), None);
    }

    /// A frameshift that starts on the terminator is an extension, and an extension does not count
    /// the stop it ends in.
    #[test]
    fn a_frameshift_that_loses_the_terminator_is_an_extension() {
        let mut c = change("-", "X", 4, 3, "MAK*");
        c.frameshift = true;
        assert_eq!(write(&c, &Edited(Some("MAKLLLL*"))), "P.1:p.Ter4LeuextTer4");
        // and one that keeps it at the same place changes nothing
        assert_eq!(write(&c, &Edited(Some("MAK*"))), "P.1:p.Ter4=");
    }

    #[test]
    fn a_deletion_is_shifted_to_its_3_prime_end() {
        // K at 3, 4 and 5: the HGVS position of one deleted K is the last
        let c = change("K", "-", 4, 4, "MAKKKPL*");
        assert_eq!(write(&c, &NONE), "P.1:p.Lys5del");
        let c = change("KK", "-", 3, 4, "MAKKKPL*");
        assert_eq!(write(&c, &NONE), "P.1:p.Lys4_Lys5del");
    }

    /// Ensembl stops short of the end of the protein twice over: it cannot see the last residue
    /// from a change that ends one residue before it, and its walk is bounded by `n - 1` residues
    /// for a change of `n`. That puts NT5C2's deletion of two Glu from a run of seven, in a
    /// protein that ends in the run, at `Glu559_Glu560` and not at the end of it (#2 of the
    /// module doc).
    #[test]
    fn a_deletion_in_a_run_at_the_end_of_the_protein_reaches_the_end() {
        // the walk is bounded: 5 residues after the change, 2 of them in the way
        let c = change("EE", "-", 3, 4, "MAEEEEEEE*");
        assert_eq!(write(&c, &NONE), "P.1:p.Glu8_Glu9del");
        // and the last residue is hidden: the change ends at 8 of 9
        let c = change("EE", "-", 7, 8, "MAEEEEEEE*");
        assert_eq!(write(&c, &NONE), "P.1:p.Glu8_Glu9del");
        // a change that ends at the last residue has nothing after it to shift over
        let c = change("EE", "-", 8, 9, "MAEEEEEEE*");
        assert_eq!(write(&c, &NONE), "P.1:p.Glu8_Glu9del");
    }

    #[test]
    fn a_replacement_is_a_delins_and_is_not_shifted() {
        let c = change("KK", "R", 3, 4, "MAKKKPL*");
        assert_eq!(write(&c, &NONE), "P.1:p.Lys3_Lys4delinsArg");
        let c = change("K", "RW", 3, 3, "MAKKKPL*");
        assert_eq!(write(&c, &NONE), "P.1:p.Lys3delinsArgTrp");
    }

    #[test]
    fn an_insertion_that_repeats_what_precedes_it_is_a_duplication() {
        // between residues 4 and 5 of MAKKKPL: after the 3' shift it follows the last K
        let c = change("-", "K", 5, 4, "MAKKKPL*");
        assert_eq!(write(&c, &NONE), "P.1:p.Lys5dup");
        let c = change("-", "KK", 5, 4, "MAKKKPL*");
        assert_eq!(write(&c, &NONE), "P.1:p.Lys4_Lys5dup");
    }

    #[test]
    fn an_insertion_is_written_between_the_residues_it_sits_between() {
        let c = change("-", "W", 4, 3, "MAKKKPL*");
        assert_eq!(write(&c, &NONE), "P.1:p.Lys3_Lys4insTrp");
    }

    /// Ensembl looks for the residue after an insertion in the protein without its terminator, so
    /// one after the last residue has none to look at and gets no notation.
    #[test]
    fn an_insertion_after_the_last_residue_is_not_described() {
        let c = change("-", "W", 7, 6, "MAKKKP*");
        assert_eq!(hgvs_protein(&c, &NONE), None);
    }

    /// A change that loses the initiator is `p.Met1?` whatever the window is. Ensembl writes the
    /// residue names of the whole window for one of more than one residue (`p.MetAla1_?2`), and for
    /// a frameshift whose first changed residue is a later one a range that runs backwards
    /// (`p.Ala2_?1`).
    #[test]
    fn a_change_that_loses_the_initiator_is_met1_unknown() {
        let mut c = change("M", "I", 1, 1, "MAKKKPL*");
        c.start_lost = true;
        assert_eq!(write(&c, &NONE), "P.1:p.Met1?");
        let mut c = change("MA", "I", 1, 2, "MAKKKPL*");
        c.start_lost = true;
        assert_eq!(write(&c, &NONE), "P.1:p.Met1?");
        // a frameshift that does it reads the same
        let mut c = change("M", "X", 1, 1, "MAKKKPL*");
        c.start_lost = true;
        c.frameshift = true;
        assert_eq!(write(&c, &Edited(Some("NKKL*"))), "P.1:p.Met1?");
        assert_eq!(write(&c, &Edited(Some("MIKKL*"))), "P.1:p.Met1?");
        // and so does a deletion of the whole window
        let mut c = change("MAK", "-", 1, 3, "MAKKKPL*");
        c.start_lost = true;
        assert_eq!(write(&c, &NONE), "P.1:p.Met1?");
    }

    /// An insertion in front of residue 1 has no residue before it, and Ensembl reads the one that
    /// is not there from the end of the protein, so its `p.Ala1_?0` names the last residue of the
    /// protein as if it were the first. The residue that was lost is the first one.
    #[test]
    fn an_insertion_in_front_of_the_initiator_names_residue_1_and_not_the_last() {
        let mut c = change("M", "IAEM", 1, 1, "MKPA*");
        c.start_lost = true;
        assert_eq!(write(&c, &NONE), "P.1:p.Met1?");
    }

    /// What residue 1 is called is read off the protein, so a protein whose first residue is not a
    /// Met (an initiator the annotation gives as something else) is named as it is.
    #[test]
    fn the_residue_that_was_lost_is_the_first_of_the_protein() {
        let mut c = change("L", "I", 1, 1, "LAKKKPL*");
        c.start_lost = true;
        assert_eq!(write(&c, &NONE), "P.1:p.Leu1?");
    }

    #[test]
    fn a_change_with_no_reference_residue_is_not_described() {
        let c = change("", "K", 3, 3, "MAKKKPL*");
        assert_eq!(hgvs_protein(&c, &NONE), None);
    }

    /// A stop that comes back where it was is nothing, whatever the rest of the window did.
    #[test]
    fn a_window_that_keeps_its_stop_ends_in_one() {
        // the stop is first in both, and what follows it is not protein
        let c = change("*", "*LX", 201, 201, "MAK*");
        assert_eq!(write(&c, &NONE), "P.1:p.Ter201delinsTer");
        let c = change("A*", "V*X", 200, 201, "MAK*");
        assert_eq!(write(&c, &NONE), "P.1:p.Ala200_Ter201delinsValTer");
    }

    /// Deterministic pseudo-random source: a seeded LCG rather than a `rand` dependency, because
    /// the sweep below has to fail the same way twice or a failure cannot be investigated.
    struct Lcg(u64);

    impl Lcg {
        fn below(&mut self, n: usize) -> usize {
            self.0 = self
                .0
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            ((self.0 >> 33) as usize) % n
        }
    }

    /// Apply a description back to the reference protein, checking on the way that every residue
    /// it names is the residue the protein has there.
    ///
    /// This is the inverse of the description, which is what makes it a check worth having: a
    /// span that is well formed, names real residues and is wrong still reconstructs the wrong
    /// protein.
    fn apply_to_protein(description: &str, protein: &str) -> String {
        let one = |aa3: &str| -> char {
            match aa3 {
                "Ala" => 'A',
                "Glu" => 'E',
                "Gly" => 'G',
                "Lys" => 'K',
                other => panic!("unexpected residue {other} in {description}"),
            }
        };
        let body = description
            .split(":p.")
            .nth(1)
            .unwrap_or_else(|| panic!("no :p. in {description}"));
        let residues = |text: &str| -> String {
            text.as_bytes()
                .chunks(3)
                .map(|c| one(std::str::from_utf8(c).unwrap()))
                .collect()
        };
        // the span, and what replaces it
        let (span, replacement, inserted_between) =
            if let Some((span, inserted)) = body.split_once("delins") {
                (span, residues(inserted), false)
            } else if let Some(span) = body.strip_suffix("del") {
                (span, String::new(), false)
            } else if let Some((span, inserted)) = body.split_once("ins") {
                (span, residues(inserted), true)
            } else if let Some(span) = body.strip_suffix("dup") {
                (span, String::new(), false)
            } else {
                // a substitution: `Gly18Lys`
                let (from, rest) = body.split_at(3);
                let (digits, to) = rest.split_at(rest.len() - 3);
                let pos: usize = digits
                    .parse()
                    .unwrap_or_else(|_| panic!("unparsable substitution {description}"));
                assert_eq!(
                    protein.chars().nth(pos - 1),
                    Some(one(from)),
                    "{description} names {from} at {pos}"
                );
                return format!("{}{}{}", &protein[..pos - 1], one(to), &protein[pos..]);
            };

        let mut bounds = Vec::new();
        for part in span.split('_') {
            let (aa3, digits) = part.split_at(3);
            let pos: usize = digits
                .parse()
                .unwrap_or_else(|_| panic!("unparsable position in {description}"));
            assert_eq!(
                protein.chars().nth(pos - 1),
                Some(one(aa3)),
                "{description} names {aa3} at {pos}, the protein has {:?} ({protein})",
                protein.chars().nth(pos - 1)
            );
            bounds.push(pos);
        }
        let (lo, hi) = (bounds[0], *bounds.last().unwrap());
        assert!(lo <= hi, "inverted span in {description}");
        if body.ends_with("dup") {
            let copy = &protein[lo - 1..hi];
            format!("{}{}{}", &protein[..hi], copy, &protein[hi..])
        } else if inserted_between {
            assert_eq!(
                hi,
                lo + 1,
                "an insertion sits between neighbours: {description}"
            );
            format!("{}{}{}", &protein[..lo], replacement, &protein[lo..])
        } else {
            format!("{}{}{}", &protein[..lo - 1], replacement, &protein[hi..])
        }
    }

    /// The property that matters, over 6,000 in-frame changes on a four-residue alphabet chosen
    /// to make repeats common: apply the description back to the reference protein and you must
    /// get the protein the variant produces. A description can be well formed, name real
    /// residues, and still fail this; that is what #96 was, and what #91 was before it.
    ///
    /// Changes start at residue 2 or later, because one that reaches the initiator is deliberately
    /// described as unknown (`?`), and end before the last residue, because one at the very end
    /// of a protein is where Ensembl's 3'-rule and the 3'-rule part company (#2 of the module
    /// doc) and is pinned on its own.
    #[test]
    fn every_description_reconstructs_the_protein_the_variant_produces() {
        const ALPHABET: &[u8] = b"AEGK";
        let mut rng = Lcg(0x5EED_1234_9ABC_DEF0);
        let mut checked = 0;
        for _ in 0..6_000 {
            let len = 8 + rng.below(14);
            let protein: String = (0..len)
                .map(|_| ALPHABET[rng.below(ALPHABET.len())] as char)
                .collect();
            let translation = format!("{protein}*");
            // residues lo..=hi replaced by `replacement` (nothing for a deletion)
            let n = 1 + rng.below(4);
            let lo = 2 + rng.below(len - n - 2);
            let hi = lo + n - 1;
            let m = rng.below(n + 3);
            let replacement: String = (0..m)
                .map(|_| ALPHABET[rng.below(ALPHABET.len())] as char)
                .collect();
            let expected = format!("{}{}{}", &protein[..lo - 1], replacement, &protein[hi..]);
            if expected == protein {
                continue;
            }

            let reference = &protein[lo - 1..hi];
            let alt = if replacement.is_empty() {
                "-"
            } else {
                &replacement
            };
            let got = write(
                &change(reference, alt, lo as u64, hi as u64, &translation),
                &NONE,
            );
            // an insertion is a window with no reference residue, between lo - 1 and lo
            assert_eq!(
                apply_to_protein(&got, &protein),
                expected,
                "{got} does not reconstruct {expected} from {protein} ({reference}/{alt} at {lo}-{hi})"
            );
            checked += 1;

            // and a pure insertion at the same place
            if !replacement.is_empty() {
                let got = write(
                    &change("-", &replacement, lo as u64, lo as u64 - 1, &translation),
                    &NONE,
                );
                let expected = format!(
                    "{}{}{}",
                    &protein[..lo - 1],
                    replacement,
                    &protein[lo - 1..]
                );
                if expected != protein {
                    assert_eq!(
                        apply_to_protein(&got, &protein),
                        expected,
                        "{got} does not reconstruct {expected} from {protein} (insertion of {replacement} before {lo})"
                    );
                    checked += 1;
                }
            }
        }
        assert!(checked > 6_000, "sweep covered only {checked} cases");
    }
}
