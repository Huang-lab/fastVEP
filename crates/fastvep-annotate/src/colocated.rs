//! VEP's known-variant columns, computed as VEP computes them.
//!
//! This is a transcription of two subs of `Bio::EnsEMBL::VEP::OutputFactory`
//! (ensembl-vep release/116): `add_colocated_variant_info` (`CLIN_SIG`,
//! `SOMATIC`, `PHENO`, `PUBMED`) and `add_colocated_frequency_data` (`AF` and the
//! population frequencies, `MAX_AF`, `MAX_AF_POPS`), one Rust function to each and
//! in the same order, so that a difference between the two can be found by reading
//! them side by side. What they do is not what one would write from a description
//! of the columns, and a plausible reading is wrong in ways that print
//! plausible numbers:
//!
//! - `AF` is formatted `%.4f` and no other population is;
//! - the frequency of an allele the record does not list is *interpolated* as
//!   `1 - total`, for `AF` of a two-allele record only;
//! - the frequencies are read from every known variant of the site for every
//!   allele, in the order the cache holds them, including a record that matched
//!   only another allele of the site (the interpolated frequency of its other
//!   allele can then be this allele's);
//! - `MAX_AF` leaves out the combined populations (`AF`, `gnomADe`, `gnomADg`) and
//!   lists every population that ties for it, once per known variant that had it;
//! - an allele-specific `CLIN_SIG` replaces the list of variant-level ones, and
//!   when the record's reference allele is not the input's the alleles of
//!   `clin_sig_allele` are reverse-complemented before the input allele is looked
//!   up in them.
//!
//! The columns are per input allele, and the same on every transcript row of it,
//! so [`known_cells`] works them out once per (variant x alt allele) and the
//! writers borrow the text. VEP works them out per row.
//!
//! Every behaviour here was checked against VEP 116.3 run on the same cache
//! (`crates/fastvep-cache/fixtures/vep_known_columns`), and `golden` tests hold it
//! to those runs.

use fastvep_cache::providers::{KnownRecord, KnownSite};
use fastvep_cache::variation::{existing_variation_rank, reverse_complement, Matched};
use fastvep_core::Allele;
use fastvep_io::known::{AlleleKnown, KnownCells, KnownColumn, KnownFlags, FREQUENCY_GROUPS};
use std::borrow::Cow;

/// Perl's idea of a string and a number, which is what the Perl being ported
/// computes with.
mod perl {
    /// Whether Perl takes a string as true: anything but `""` and `"0"`. So
    /// `"0.0000"` is true and `"0"` is not, which is why a population whose
    /// whole column reads `0` is skipped and one that reads `A:0` is not.
    pub fn is_true(s: &str) -> bool {
        !s.is_empty() && s != "0"
    }

    /// A string as a number: the longest numeric prefix, 0 if there is none
    /// (`"1e-05"` is 1e-05, `"0.5abc"` is 0.5, `"abc"` is 0).
    pub fn number(s: &str) -> f64 {
        let s = s.trim_start();
        let bytes = s.as_bytes();
        let mut i = 0;
        if matches!(bytes.first(), Some(b'+' | b'-')) {
            i += 1;
        }
        let digits = |from: usize| {
            bytes[from..]
                .iter()
                .take_while(|b| b.is_ascii_digit())
                .count()
        };
        let int_digits = digits(i);
        i += int_digits;
        let mut frac_digits = 0;
        if bytes.get(i) == Some(&b'.') {
            frac_digits = digits(i + 1);
            if int_digits > 0 || frac_digits > 0 {
                i += 1 + frac_digits;
            }
        }
        if int_digits == 0 && frac_digits == 0 {
            return 0.0;
        }
        if matches!(bytes.get(i), Some(b'e' | b'E')) {
            let mut j = i + 1;
            if matches!(bytes.get(j), Some(b'+' | b'-')) {
                j += 1;
            }
            let exp_digits = digits(j);
            if exp_digits > 0 {
                i = j + exp_digits;
            }
        }
        s[..i].parse().unwrap_or(0.0)
    }

    /// `sprintf("%.4f", $x)`.
    pub fn four_decimals(x: f64) -> String {
        format!("{x:.4}")
    }

    /// A number as Perl prints it: `%.15g`, the shortest of fixed and scientific
    /// notation with trailing zeros dropped (`1 - 0.3` is `0.7`, not
    /// `0.7000000000000001`).
    pub fn string(x: f64) -> String {
        if x == 0.0 {
            return "0".to_string();
        }
        if !x.is_finite() {
            return if x.is_nan() {
                "NaN".to_string()
            } else if x > 0.0 {
                "Inf".to_string()
            } else {
                "-Inf".to_string()
            };
        }
        // 15 significant digits, correctly rounded: d.dddddddddddddde[-]x
        let sci = format!("{x:.14e}");
        let (mantissa, exponent) = sci.split_once('e').expect("exponent form");
        let exponent: i32 = exponent.parse().expect("integer exponent");
        let (sign, mantissa) = match mantissa.strip_prefix('-') {
            Some(m) => ("-", m),
            None => ("", mantissa),
        };
        let digits: String = mantissa.chars().filter(|c| *c != '.').collect();
        let digits = digits.trim_end_matches('0');
        let digits = if digits.is_empty() { "0" } else { digits };
        // `%g` uses scientific notation when the exponent is below -4 or at least
        // the precision.
        if !(-4..15).contains(&exponent) {
            let mut out = String::from(sign);
            out.push_str(&digits[..1]);
            if digits.len() > 1 {
                out.push('.');
                out.push_str(&digits[1..]);
            }
            out.push('e');
            out.push(if exponent < 0 { '-' } else { '+' });
            out.push_str(&format!("{:02}", exponent.abs()));
            return out;
        }
        let mut out = String::from(sign);
        if exponent < 0 {
            out.push_str("0.");
            out.push_str(&"0".repeat((-exponent - 1) as usize));
            out.push_str(digits);
        } else {
            let int_len = exponent as usize + 1;
            if digits.len() <= int_len {
                out.push_str(digits);
                out.push_str(&"0".repeat(int_len - digits.len()));
            } else {
                out.push_str(&digits[..int_len]);
                out.push('.');
                out.push_str(&digits[int_len..]);
            }
        }
        out
    }

    /// `split(',', $s)`: Perl drops the empty fields at the end, so a string of
    /// nothing but separators (or nothing) has no fields at all.
    pub fn split(s: &str, sep: char) -> impl Iterator<Item = &str> {
        let trimmed = s.trim_end_matches(sep);
        (!trimmed.is_empty())
            .then(|| trimmed.split(sep))
            .into_iter()
            .flatten()
    }
}

/// What `add_colocated_variant_info` and `add_colocated_frequency_data` write into
/// VEP's output hash for one row, restricted to the keys they own.
#[derive(Default)]
struct Hash {
    /// `$tmp->{CLIN_SIG}`: the variant-level values, pushed in order.
    clin_sig: Vec<String>,
    /// `$tmp->{CLIN_SIG}` as `join(';', keys %clin_sigs)`, which replaces the list
    /// when any allele-specific value was found.
    clin_sig_allele: Option<String>,
    pubmed: Vec<String>,
    somatic: Vec<bool>,
    pheno: Vec<bool>,
    /// `$hash->{"${key}_AF"}`: one list per output column, values unique.
    frequencies: Vec<(KnownColumn, Vec<String>)>,
    /// `$hash->{MAX_AF}`: the text of the cache value, or `0`.
    max_af: Option<String>,
    /// `$hash->{MAX_AF_POPS}`: cache column names, repeats and all.
    max_af_pops: Vec<&'static str>,
}

/// The known variants of a site in the order `Existing_variation` is sorted into:
/// germline before somatic, then by the prefix of the name, ties in cache order.
/// `add_colocated_variant_info` walks them this way and
/// `add_colocated_frequency_data` does not.
fn sorted_by_name(site: &KnownSite) -> Vec<&KnownRecord> {
    let mut sorted: Vec<&KnownRecord> = site.records.iter().collect();
    sorted.sort_by_key(|k| existing_variation_rank(&k.record.variation_name, k.record.somatic));
    sorted
}

/// The input allele at `alt_index` has a known variant that matched it, or the
/// known variant has no alleles at all (`COSMIC_MUTATION`) and so no matched
/// alleles to check (`next unless (grep {$_->{a_allele} eq $this_allele} @$matched)`).
fn matched_for_allele(known: &KnownRecord, alt_index: usize) -> bool {
    known.matched.covers(alt_index)
}

/// `add_colocated_variant_info`.
///
/// `known` is the site's known variants sorted as above, `this_allele` the
/// allele of the row and `ref_allele` the reference allele VEP prints for it
/// (`-` for an insertion).
fn colocated_variant_info(
    hash: &mut Hash,
    known: &[&KnownRecord],
    alt_index: usize,
    this_allele: &str,
    ref_allele: &str,
    flags: KnownFlags,
) {
    let mut clin_sig_allele_exists = false;
    // `%clin_sigs`: Perl's hash order is not defined, and VEP's own output for two
    // distinct values differs from run to run. First seen is the one order that
    // can be reproduced.
    let mut clin_sigs: Vec<String> = Vec::new();

    for ex in known {
        // check allele match
        if !matched_for_allele(ex, alt_index) {
            continue;
        }
        let ex = &*ex.record;

        if let (Some(clin_sig_allele), false) = (&ex.clin_sig_allele, flags.no_clin_sig_allele) {
            // The significances are written against the record's reference
            // allele. When the input's reference is another allele, VEP takes the
            // record to be on the other strand.
            let allele_match = !ex
                .clin_sig_ref_allele
                .as_deref()
                .is_some_and(|r| perl::is_true(r) && r != ref_allele);

            // `%cs_hash`: the values given for each allele, joined by `,`
            let mut cs_hash: Vec<(String, String)> = Vec::new();
            for cs in perl::split(clin_sig_allele, ';') {
                // `split(':', $cs)`: the allele, then its value (which is empty
                // when the part has no `:`, as `.= undef` appends nothing)
                let mut parts = perl::split(cs, ':');
                let Some(allele) = parts.next() else {
                    continue;
                };
                let value = parts.next().unwrap_or("");
                let allele = if allele_match {
                    allele.to_string()
                } else {
                    reverse_complement(allele)
                };
                match cs_hash.iter_mut().find(|(a, _)| *a == allele) {
                    Some((_, existing)) => {
                        if !existing.is_empty() {
                            existing.push(',');
                        }
                        existing.push_str(value);
                    }
                    None => cs_hash.push((allele, value.to_string())),
                }
            }
            if let Some((_, value)) = cs_hash.into_iter().find(|(a, _)| a == this_allele) {
                if !clin_sigs.contains(&value) {
                    clin_sigs.push(value);
                }
            }
            clin_sig_allele_exists = true;
        }

        // clin sig and pubmed?
        if let Some(clin_sig) = ex.clin_sig.as_deref() {
            if perl::is_true(clin_sig) && !clin_sig_allele_exists {
                hash.clin_sig
                    .extend(perl::split(clin_sig, ',').map(String::from));
            }
        }
        // The cache holds `pubmed` split already; a field that was `0` is no
        // field to Perl, and splits to one `0` here.
        if flags.pubmed && !(ex.pubmed.is_empty() || ex.pubmed == ["0"]) {
            hash.pubmed.extend(ex.pubmed.iter().cloned());
        }

        // somatic? phenotype or disease?
        hash.somatic.push(ex.somatic);
        hash.pheno.push(ex.phenotype_or_disease);
    }

    if !clin_sigs.is_empty() && !flags.no_clin_sig_allele {
        hash.clin_sig_allele = Some(clin_sigs.join(";"));
    }
}

/// Add `f` to the list of an output column unless it is there (`merge_arrays`).
fn merge_into(hash: &mut Hash, column: KnownColumn, f: &str) {
    let list = match hash.frequencies.iter().position(|(c, _)| *c == column) {
        Some(i) => &mut hash.frequencies[i].1,
        None => {
            hash.frequencies.push((column, Vec::new()));
            &mut hash.frequencies.last_mut().expect("just pushed").1
        }
    };
    if !list.iter().any(|v| v == f) {
        list.push(f.to_string());
    }
}

/// `add_colocated_frequency_data`, for one known variant of the site.
fn colocated_frequency_data(
    hash: &mut Hash,
    known: &KnownRecord,
    alt_index: usize,
    this_allele: &str,
    flags: KnownFlags,
) {
    let any_group = FREQUENCY_GROUPS.iter().any(|g| flags.group(g.flag));
    if !(any_group || flags.max_af) {
        return;
    }
    let ex = &*known.record;
    let ex_alleles: Vec<&str> = ex.allele_string.split('/').collect();

    // `my ($matched_allele) = grep {...} @{$ex->{matched_alleles} || []}`: the
    // first of the record's alleles that matched this input allele.
    let matched_allele: Option<&str> = match &known.matched {
        Matched::Alleles(matches) => matches
            .iter()
            .find(|m| m.alt_index == alt_index)
            .map(|m| m.record_allele.as_str()),
        Matched::AllAlleles => None,
    };

    // `return $hash unless $matched_allele || (grep {$_ eq 'af'} @keys)`, where
    // `@keys` is every group under `--max_af` and the groups asked for otherwise
    if matched_allele.is_none() && !(flags.max_af || flags.af) {
        return;
    }

    let mut max_af = 0.0_f64;
    let mut max_af_text: Option<String> = None;
    let mut max_af_pops: Vec<&'static str> = Vec::new();
    // Scratch for each population read, kept across them: this runs for every
    // known variant of every variant that has one.
    let mut freq_data: Vec<(&str, Cow<'_, str>)> = Vec::new();
    let mut remaining: Vec<&str> = Vec::new();

    for group in &FREQUENCY_GROUPS {
        let recorded = flags.group(group.flag);
        // `@keys = grep {$self->{$_}} @keys unless $self->{max_af}`
        if !recorded && !flags.max_af {
            continue;
        }
        for frequency_key in group.keys {
            let Some(data) = ex
                .frequencies
                .get(frequency_key.key)
                .filter(|v| perl::is_true(v))
            else {
                continue;
            };

            // get the frequencies for each allele
            freq_data.clear();
            let mut total = 0.0;
            // `%remaining = map {$_ => 1} @ex_alleles`: which alleles have no
            // frequency of their own
            remaining.clear();
            for a in &ex_alleles {
                if !remaining.contains(a) {
                    remaining.push(a);
                }
            }
            for pair in perl::split(data, ',') {
                let mut parts = perl::split(pair, ':');
                let Some(a) = parts.next() else {
                    continue;
                };
                let mut f = Cow::Borrowed(parts.next().unwrap_or(""));
                if frequency_key.key == "AF" {
                    // this format is just to keep old compatibility with dbSNP import
                    f = Cow::Owned(perl::four_decimals(perl::number(&f)));
                }
                total += perl::number(&f);
                match freq_data.iter_mut().find(|(k, _)| *k == a) {
                    Some(slot) => slot.1 = f,
                    None => freq_data.push((a, f)),
                }
                remaining.retain(|r| *r != a);
            }

            // interpolate the frequency for the remaining allele if there is only 1;
            // that can be done reliably for AF only, as only the minor AF is stored
            let interpolated =
                ex_alleles.len() == 2 && remaining.len() == 1 && frequency_key.key == "AF";
            if interpolated {
                freq_data.push((remaining[0], Cow::Owned(perl::string(1.0 - total))));
            }

            let lookup = |allele: &str| {
                freq_data
                    .iter()
                    .rev()
                    .find(|(k, _)| *k == allele)
                    .map(|(_, f)| &**f)
            };
            let from_matched = matched_allele.and_then(lookup);
            let from_interpolation = if interpolated {
                lookup(this_allele).filter(|f| perl::is_true(f))
            } else {
                None
            };
            let Some(f) = from_matched.or(from_interpolation) else {
                continue;
            };

            // record the frequency if requested
            if recorded {
                if let Some(column) = frequency_key.column {
                    merge_into(hash, column, f);
                }
            }

            // update max_af data if required; no combined-level population counts
            if flags.max_af
                && !matches!(
                    frequency_key.key,
                    "AF" | "ExAC" | "ExAC_Adj" | "gnomADe" | "gnomADg"
                )
            {
                let value = perl::number(f);
                if value > max_af {
                    max_af = value;
                    max_af_text = Some(f.to_string());
                    max_af_pops = vec![frequency_key.key];
                } else if value == max_af && !max_af_pops.contains(&frequency_key.key) {
                    max_af_pops.push(frequency_key.key);
                }
            }
        }
    }

    // add/update max_af info
    if flags.max_af && !max_af_pops.is_empty() {
        // `my $current_max = $hash->{MAX_AF} ||= 0`
        if !hash.max_af.as_deref().is_some_and(perl::is_true) {
            hash.max_af = Some("0".to_string());
        }
        let current_max = hash.max_af.as_deref().map_or(0.0, perl::number);
        if max_af > current_max {
            hash.max_af = Some(max_af_text.unwrap_or_else(|| "0".to_string()));
            hash.max_af_pops.clear();
        }
        if max_af >= current_max {
            hash.max_af_pops.extend(max_af_pops);
        }
    }
}

impl Hash {
    /// The text each column has in VEP's tab file (`convert_arrayref`: a list is
    /// joined with `,`), for the columns that have a value.
    fn into_cells(self, flags: KnownFlags) -> KnownCells {
        let mut cells = KnownCells::with_capacity(self.frequencies.len() + 8);
        for (column, mut values) in self.frequencies {
            // a list of one is the value itself: no need to build the join of it
            match values.len() {
                0 => {}
                1 => cells.set(column, values.pop().expect("one value")),
                _ => cells.set(column, values.join(",")),
            }
        }
        if let Some(max_af) = self.max_af {
            cells.set(KnownColumn::MaxAf, max_af);
        }
        if !self.max_af_pops.is_empty() {
            cells.set(KnownColumn::MaxAfPops, self.max_af_pops.join(","));
        }
        // post-process to remove all-0, e.g. SOMATIC
        if let Some(allele_specific) = self.clin_sig_allele {
            if !flags.no_clin_sig_allele {
                cells.set(KnownColumn::ClinSig, allele_specific);
            }
        } else if self.clin_sig.iter().any(|v| perl::is_true(v)) {
            cells.set(KnownColumn::ClinSig, self.clin_sig.join(","));
        }
        let flags_cell = |list: &[bool]| {
            list.iter().any(|b| *b).then(|| {
                list.iter()
                    .map(|b| if *b { "1" } else { "0" })
                    .collect::<Vec<_>>()
                    .join(",")
            })
        };
        if let Some(cell) = flags_cell(&self.somatic) {
            cells.set(KnownColumn::Somatic, cell);
        }
        if let Some(cell) = flags_cell(&self.pheno) {
            cells.set(KnownColumn::Pheno, cell);
        }
        if self.pubmed.iter().any(|v| perl::is_true(v)) {
            cells.set(KnownColumn::Pubmed, self.pubmed.join(","));
        }
        cells
    }
}

/// VEP's known-variant columns for each alt allele of one input site.
///
/// `site` is what the variation cache attached to the site, `alleles` the alt
/// alleles in the order the lookup was made with (`AlleleMatch::alt_index` counts
/// into them) and `ref_allele` the site's reference allele as VEP writes it. Only
/// alleles with something to show are returned.
pub fn known_cells(
    site: &KnownSite,
    ref_allele: &str,
    alleles: &[Allele],
    flags: KnownFlags,
) -> Vec<AlleleKnown> {
    if site.is_empty() || !flags.checks_existing() {
        return Vec::new();
    }
    let by_name = sorted_by_name(site);
    let mut out = Vec::new();
    for (alt_index, allele) in alleles.iter().enumerate() {
        let this_allele = allele.to_string();
        let mut hash = Hash::default();
        colocated_variant_info(
            &mut hash,
            &by_name,
            alt_index,
            &this_allele,
            ref_allele,
            flags,
        );
        // in the order of the cache, as `foreach my $ex(@{$vf->{existing}})`
        for known in &site.records {
            colocated_frequency_data(&mut hash, known, alt_index, &this_allele, flags);
        }
        let cells = hash.into_cells(flags);
        if !cells.is_empty() {
            out.push(AlleleKnown {
                allele: allele.clone(),
                cells,
            });
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use fastvep_cache::variation::{AlleleMatch, VariationCacheRecord};
    use std::collections::HashMap;
    use std::sync::Arc;

    fn record(name: &str, allele_string: &str) -> VariationCacheRecord {
        VariationCacheRecord {
            variation_name: name.into(),
            failed: false,
            somatic: false,
            start: 100,
            end: 100,
            allele_string: allele_string.into(),
            strand: 1,
            minor_allele: None,
            minor_allele_freq: None,
            clin_sig: None,
            clin_sig_allele: None,
            clin_sig_ref_allele: None,
            phenotype_or_disease: false,
            pubmed: vec![],
            frequencies: HashMap::new(),
        }
    }

    fn with_freq(mut r: VariationCacheRecord, key: &str, value: &str) -> VariationCacheRecord {
        r.frequencies.insert(key.into(), value.into());
        r
    }

    /// A known variant that matched input allele `alt_index` through its allele `b`.
    fn matched(r: VariationCacheRecord, alt_index: usize, b: &str) -> KnownRecord {
        KnownRecord {
            record: Arc::new(r),
            matched: Matched::Alleles(vec![AlleleMatch {
                alt_index,
                record_allele: b.into(),
            }]),
        }
    }

    fn alleles(alts: &[&str]) -> Vec<Allele> {
        alts.iter()
            .map(|a| {
                if *a == "-" {
                    Allele::Deletion
                } else {
                    Allele::Sequence(a.as_bytes().to_vec())
                }
            })
            .collect()
    }

    fn cells(
        records: Vec<KnownRecord>,
        reference: &str,
        alts: &[&str],
        flags: KnownFlags,
    ) -> Vec<Vec<(&'static str, String)>> {
        let site = KnownSite { records };
        let found = known_cells(&site, reference, &alleles(alts), flags);
        (0..alts.len())
            .map(|i| {
                found
                    .iter()
                    .find(|k| k.allele == alleles(alts)[i])
                    .map(|k| {
                        k.cells
                            .iter()
                            .map(|(c, v)| (c.name(), v.to_string()))
                            .collect()
                    })
                    .unwrap_or_default()
            })
            .collect()
    }

    fn get<'a>(row: &'a [(&'static str, String)], column: &str) -> Option<&'a str> {
        row.iter()
            .find(|(c, _)| *c == column)
            .map(|(_, v)| v.as_str())
    }

    fn af() -> KnownFlags {
        KnownFlags {
            af: true,
            ..Default::default()
        }
    }

    fn af_1kg() -> KnownFlags {
        KnownFlags {
            af_1kg: true,
            ..Default::default()
        }
    }

    fn max_af() -> KnownFlags {
        KnownFlags {
            max_af: true,
            ..Default::default()
        }
    }

    // ---- Perl's numbers and strings -------------------------------------------

    #[test]
    fn a_number_is_printed_as_perl_prints_it() {
        // `%.15g`, which drops the noise of a subtraction
        assert_eq!(perl::string(1.0 - 0.3), "0.7");
        assert_eq!(perl::string(1.0 - 0.8765), "0.1235");
        assert_eq!(perl::string(1.0), "1");
        assert_eq!(perl::string(0.0), "0");
        assert_eq!(perl::string(1.0 / 3.0), "0.333333333333333");
        assert_eq!(perl::string(2.0 / 3.0), "0.666666666666667");
        assert_eq!(perl::string(1e-5), "1e-05");
        assert_eq!(perl::string(0.0001), "0.0001");
        assert_eq!(perl::string(1.5e-7), "1.5e-07");
        assert_eq!(perl::string(123456789012345.0), "123456789012345");
        assert_eq!(perl::string(1e15), "1e+15");
        assert_eq!(perl::string(-0.25), "-0.25");
        assert_eq!(perl::string(0.1 + 0.2), "0.3");
    }

    #[test]
    fn a_string_is_a_number_as_far_as_it_reads_as_one() {
        assert_eq!(perl::number("0.5"), 0.5);
        assert_eq!(perl::number("1e-05"), 1e-5);
        assert_eq!(perl::number("8.291e-06"), 8.291e-6);
        assert_eq!(perl::number("  .25"), 0.25);
        assert_eq!(perl::number("0.5abc"), 0.5);
        assert_eq!(perl::number("2."), 2.0);
        assert_eq!(perl::number("abc"), 0.0);
        assert_eq!(perl::number(""), 0.0);
        assert_eq!(perl::number("-1.5e2x"), -150.0);
        assert_eq!(
            perl::number("1e"),
            1.0,
            "an exponent with no digits is not read"
        );
    }

    #[test]
    fn four_decimals_rounds_the_way_printf_does() {
        assert_eq!(perl::four_decimals(0.123456), "0.1235");
        assert_eq!(
            perl::four_decimals(0.00005),
            "0.0001",
            "5.0000000000000002e-05 is above the tie"
        );
        assert_eq!(perl::four_decimals(1e-5), "0.0000");
        assert_eq!(perl::four_decimals(0.9999), "0.9999");
        // an exact tie in binary (1/32) goes to the even digit, as glibc's does
        assert_eq!(perl::four_decimals(0.03125), "0.0312");
        assert_eq!(perl::four_decimals(0.09375), "0.0938");
    }

    #[test]
    fn truth_is_perls() {
        assert!(!perl::is_true(""));
        assert!(!perl::is_true("0"));
        assert!(perl::is_true("0.0"));
        assert!(perl::is_true("0.0000"));
        assert!(perl::is_true("A:0"));
        assert_eq!(perl::split("a,b,,", ',').collect::<Vec<_>>(), ["a", "b"]);
        assert_eq!(perl::split("a,,b", ',').collect::<Vec<_>>(), ["a", "", "b"]);
        assert_eq!(perl::split("", ',').count(), 0);
        assert_eq!(perl::split(",,", ',').count(), 0);
    }

    // ---- add_colocated_frequency_data -----------------------------------------

    /// `AF` is `%.4f`, and `AFR` (like every other key) is printed as the cache has it.
    #[test]
    fn af_is_formatted_to_four_decimals_and_other_populations_are_not() {
        let r = with_freq(
            with_freq(record("rs1", "C/A"), "AF", "A:0.123456"),
            "AFR",
            "A:0.00123456",
        );
        let flags = KnownFlags {
            af_1kg: true,
            ..Default::default()
        };
        let out = cells(vec![matched(r, 0, "A")], "C", &["A"], flags);
        assert_eq!(get(&out[0], "AF"), Some("0.1235"));
        assert_eq!(get(&out[0], "AFR_AF"), Some("0.00123456"));
    }

    /// Only the minor allele's AF is stored: the other allele of a two-allele record
    /// gets `1 - total`, printed as Perl prints a number.
    #[test]
    fn af_of_the_other_allele_of_a_two_allele_record_is_interpolated() {
        let r = with_freq(record("rs1", "C/A"), "AF", "C:0.8765");
        let out = cells(vec![matched(r, 0, "A")], "C", &["A"], af());
        assert_eq!(get(&out[0], "AF"), Some("0.1235"));
    }

    /// Interpolation is for two alleles: a record with three has no AF for the one it lists nothing for.
    #[test]
    fn af_of_a_three_allele_record_is_not_interpolated() {
        let r = with_freq(record("rs1", "C/A/G"), "AF", "A:0.2");
        let out = cells(
            vec![KnownRecord {
                record: Arc::new(r),
                matched: Matched::Alleles(vec![
                    AlleleMatch {
                        alt_index: 0,
                        record_allele: "A".into(),
                    },
                    AlleleMatch {
                        alt_index: 1,
                        record_allele: "G".into(),
                    },
                ]),
            }],
            "C",
            &["A", "G"],
            af(),
        );
        assert_eq!(get(&out[0], "AF"), Some("0.2000"));
        assert_eq!(get(&out[1], "AF"), None);
    }

    /// Only AF is interpolated; every other population must list the allele.
    #[test]
    fn another_population_is_not_interpolated() {
        let r = with_freq(record("rs1", "C/A"), "AFR", "C:0.7");
        let out = cells(vec![matched(r, 0, "A")], "C", &["A"], af_1kg());
        assert_eq!(get(&out[0], "AFR_AF"), None);
    }

    /// The key is looked up by the record's own allele (`b_allele`), which for a
    /// record on the other strand is not the input's.
    #[test]
    fn the_frequency_is_looked_up_by_the_alleles_of_the_record() {
        let r = with_freq(record("rs1", "G/T"), "AFR", "T:0.25");
        // the input is C>A, on the other strand of the record
        let out = cells(vec![matched(r, 0, "T")], "C", &["A"], af_1kg());
        assert_eq!(get(&out[0], "AFR_AF"), Some("0.25"));
    }

    /// A column that is just `0` is false to Perl and skipped; `A:0` is a frequency of 0.
    #[test]
    fn a_population_of_0_is_skipped_and_one_of_allele_colon_0_is_not() {
        let r = with_freq(with_freq(record("rs1", "C/A"), "AFR", "0"), "AMR", "A:0");
        let out = cells(vec![matched(r, 0, "A")], "C", &["A"], af_1kg());
        assert_eq!(get(&out[0], "AFR_AF"), None);
        assert_eq!(get(&out[0], "AMR_AF"), Some("0"));
    }

    /// Two known variants with the same frequency print it once; with different ones, in cache order.
    #[test]
    fn frequencies_of_several_known_variants_are_merged_in_cache_order() {
        let a = with_freq(record("rsB", "C/A"), "AF", "A:0.2");
        let b = with_freq(record("COSV1", "C/A"), "AF", "A:0.9");
        let mut b = b;
        b.somatic = true;
        let c = with_freq(record("rsA", "C/A"), "AF", "A:0.1");
        let d = with_freq(record("rsD", "C/A"), "AF", "A:0.2");
        let out = cells(
            vec![
                matched(a, 0, "A"),
                matched(b, 0, "A"),
                matched(c, 0, "A"),
                matched(d, 0, "A"),
            ],
            "C",
            &["A"],
            af(),
        );
        // cache order, not Existing_variation order (germline first): 0.2, 0.9, 0.1, and 0.2 again
        assert_eq!(get(&out[0], "AF"), Some("0.2000,0.9000,0.1000"));
    }

    /// Only the groups asked for are written; `--max_af` reads the others without writing them.
    #[test]
    fn only_the_groups_asked_for_are_written() {
        let r = with_freq(
            with_freq(record("rs1", "C/A"), "AFR", "A:0.3"),
            "gnomADg_NFE",
            "A:0.4",
        );
        let flags = KnownFlags {
            af_gnomadg: true,
            ..Default::default()
        };
        let out = cells(vec![matched(r, 0, "A")], "C", &["A"], flags);
        assert_eq!(get(&out[0], "AFR_AF"), None);
        assert_eq!(get(&out[0], "gnomADg_NFE_AF"), Some("0.4"));
    }

    /// A known variant attached through one allele of the site is read for the
    /// others: the interpolated frequency of its other allele is theirs when they
    /// are that allele. (Which records the cache attaches is not under test here;
    /// VEP does not attach `C/A` to `G>A,C`, so this is the port's own reading of
    /// `add_colocated_frequency_data`, which no known-variant record can reach for
    /// a single-base site.)
    #[test]
    fn a_record_that_matched_another_allele_still_gives_the_interpolated_frequency() {
        // the record C/A attached through A; for the allele C, which is the record's
        // own reference allele, AF is 1 - 0.3
        let r = with_freq(record("rs1", "C/A"), "AF", "A:0.3000");
        let out = cells(vec![matched(r, 0, "A")], "G", &["A", "C"], af());
        assert_eq!(get(&out[0], "AF"), Some("0.3000"));
        assert_eq!(get(&out[1], "AF"), Some("0.7"));
    }

    /// ... but only for `AF`, so `--af_1kg` alone does not read an unmatched record at all.
    #[test]
    fn an_unmatched_record_is_read_only_when_af_is_among_the_groups() {
        let r = with_freq(record("rs1", "C/A"), "AF", "A:0.3000");
        let out = cells(vec![matched(r.clone(), 0, "A")], "G", &["A", "C"], af_1kg());
        assert_eq!(get(&out[1], "AF"), None, "af_1kg is not `af`");
        let flags = KnownFlags {
            af_1kg: true,
            max_af: true,
            ..Default::default()
        };
        let out = cells(vec![matched(r, 0, "A")], "G", &["A", "C"], flags);
        assert_eq!(
            get(&out[1], "AF"),
            Some("0.7"),
            "--max_af puts `af` among the keys"
        );
    }

    /// A record with no alleles has no allele to look a frequency up by.
    #[test]
    fn a_record_with_no_alleles_has_no_frequency() {
        let r = with_freq(record("COSV1", "COSMIC_MUTATION"), "AF", "A:0.3");
        let out = cells(
            vec![KnownRecord {
                record: Arc::new(r),
                matched: Matched::AllAlleles,
            }],
            "C",
            &["A"],
            af(),
        );
        assert_eq!(get(&out[0], "AF"), None);
    }

    // ---- MAX_AF ----------------------------------------------------------------

    #[test]
    fn max_af_is_the_largest_population_and_leaves_out_the_combined_ones() {
        let r = with_freq(
            with_freq(
                with_freq(
                    with_freq(record("rs1", "C/A"), "AF", "A:0.5"),
                    "gnomADe",
                    "A:0.9",
                ),
                "AFR",
                "A:0.01",
            ),
            "gnomADe_NFE",
            "A:0.015",
        );
        let out = cells(vec![matched(r, 0, "A")], "C", &["A"], max_af());
        assert_eq!(get(&out[0], "MAX_AF"), Some("0.015"));
        assert_eq!(get(&out[0], "MAX_AF_POPS"), Some("gnomADe_NFE"));
    }

    /// Populations that tie are all listed, in the order VEP reads the groups.
    #[test]
    fn max_af_lists_every_population_that_ties() {
        let r = with_freq(
            with_freq(
                with_freq(record("rs1", "C/A"), "EAS", "A:0.02"),
                "AMR",
                "A:0.02",
            ),
            "gnomADg_AMR",
            "A:0.02",
        );
        let out = cells(vec![matched(r, 0, "A")], "C", &["A"], max_af());
        assert_eq!(get(&out[0], "MAX_AF"), Some("0.02"));
        assert_eq!(get(&out[0], "MAX_AF_POPS"), Some("AMR,EAS,gnomADg_AMR"));
    }

    /// A population of 0 is a population: `MAX_AF` is 0 and lists the zeros.
    #[test]
    fn max_af_of_only_zeros_is_zero_and_lists_them() {
        let r = with_freq(
            with_freq(record("rs1", "C/A"), "AFR", "A:0"),
            "gnomADg_NFE",
            "A:0",
        );
        let out = cells(vec![matched(r, 0, "A")], "C", &["A"], max_af());
        assert_eq!(get(&out[0], "MAX_AF"), Some("0"));
        assert_eq!(get(&out[0], "MAX_AF_POPS"), Some("AFR,gnomADg_NFE"));
    }

    /// The same population at the same maximum in two known variants is listed twice.
    #[test]
    fn max_af_pops_repeats_a_population_across_known_variants() {
        let a = with_freq(record("rs1", "C/A"), "AFR", "A:0.3");
        let b = with_freq(record("rs2", "C/A"), "AFR", "A:0.3");
        let out = cells(
            vec![matched(a, 0, "A"), matched(b, 0, "A")],
            "C",
            &["A"],
            max_af(),
        );
        assert_eq!(get(&out[0], "MAX_AF_POPS"), Some("AFR,AFR"));
    }

    /// A later known variant with a higher maximum replaces the list; a lower one does not touch it.
    #[test]
    fn a_higher_maximum_in_a_later_known_variant_replaces_the_list() {
        let a = with_freq(record("rs1", "C/A"), "AFR", "A:0.3");
        let b = with_freq(record("rs2", "C/A"), "EUR", "A:0.5");
        let c = with_freq(record("rs3", "C/A"), "SAS", "A:0.1");
        let out = cells(
            vec![matched(a, 0, "A"), matched(b, 0, "A"), matched(c, 0, "A")],
            "C",
            &["A"],
            max_af(),
        );
        assert_eq!(get(&out[0], "MAX_AF"), Some("0.5"));
        assert_eq!(get(&out[0], "MAX_AF_POPS"), Some("EUR"));
    }

    /// The value is the cache's own text.
    #[test]
    fn max_af_is_written_as_the_cache_has_it() {
        let r = with_freq(record("rs1", "C/A"), "gnomADg_NFE", "A:1e-05");
        let out = cells(vec![matched(r, 0, "A")], "C", &["A"], max_af());
        assert_eq!(get(&out[0], "MAX_AF"), Some("1e-05"));
    }

    // ---- add_colocated_variant_info -----------------------------------------------

    fn checks() -> KnownFlags {
        KnownFlags {
            check_existing: true,
            ..Default::default()
        }
    }

    #[test]
    fn clin_sig_of_a_record_is_split_on_commas() {
        let mut r = record("rs1", "C/A");
        r.clin_sig = Some("likely_benign,benign".into());
        let out = cells(vec![matched(r, 0, "A")], "C", &["A"], checks());
        assert_eq!(get(&out[0], "CLIN_SIG"), Some("likely_benign,benign"));
    }

    #[test]
    fn clin_sig_of_several_records_is_one_list() {
        let mut a = record("rs1", "C/A");
        a.clin_sig = Some("pathogenic".into());
        let mut b = record("rs2", "C/A");
        b.clin_sig = Some("benign,likely_benign".into());
        let out = cells(
            vec![matched(a, 0, "A"), matched(b, 0, "A")],
            "C",
            &["A"],
            checks(),
        );
        assert_eq!(
            get(&out[0], "CLIN_SIG"),
            Some("pathogenic,benign,likely_benign")
        );
    }

    /// The significance given for the allele is the one that is written, instead of the list.
    #[test]
    fn clin_sig_allele_picks_the_value_of_this_allele() {
        let mut r = record("rs1", "C/A/G");
        r.clin_sig = Some("pathogenic,likely_pathogenic,benign".into());
        r.clin_sig_allele = Some("A:pathogenic;A:likely_pathogenic;G:benign".into());
        r.clin_sig_ref_allele = Some("C".into());
        let known = KnownRecord {
            record: Arc::new(r),
            matched: Matched::Alleles(vec![
                AlleleMatch {
                    alt_index: 0,
                    record_allele: "A".into(),
                },
                AlleleMatch {
                    alt_index: 1,
                    record_allele: "G".into(),
                },
            ]),
        };
        let out = cells(vec![known], "C", &["A", "G"], checks());
        assert_eq!(
            get(&out[0], "CLIN_SIG"),
            Some("pathogenic,likely_pathogenic")
        );
        assert_eq!(get(&out[1], "CLIN_SIG"), Some("benign"));
    }

    /// `--clin_sig_allele 0` goes back to the list of the variant.
    #[test]
    fn clin_sig_allele_zero_reads_the_variants_list() {
        let mut r = record("rs1", "C/A");
        r.clin_sig = Some("pathogenic".into());
        r.clin_sig_allele = Some("A:benign".into());
        r.clin_sig_ref_allele = Some("C".into());
        let flags = KnownFlags {
            check_existing: true,
            no_clin_sig_allele: true,
            ..Default::default()
        };
        let out = cells(vec![matched(r, 0, "A")], "C", &["A"], flags);
        assert_eq!(get(&out[0], "CLIN_SIG"), Some("pathogenic"));
    }

    /// An allele-specific entry that is not for this allele leaves nothing, and still
    /// keeps the variant-level list of the records after it out (`$clin_sig_allele_exists`).
    #[test]
    fn clin_sig_allele_for_another_allele_gives_no_clin_sig_and_blocks_later_lists() {
        let mut a = record("rsX1", "C/A");
        a.clin_sig = Some("benign".into());
        let mut b = record("rsX2", "C/A");
        b.clin_sig = Some("pathogenic".into());
        b.clin_sig_allele = Some("G:pathogenic".into());
        b.clin_sig_ref_allele = Some("C".into());
        let mut c = record("rsX3", "C/A");
        c.clin_sig = Some("likely_benign".into());
        let out = cells(
            vec![matched(a, 0, "A"), matched(b, 0, "A"), matched(c, 0, "A")],
            "C",
            &["A"],
            checks(),
        );
        // rsX1 is read before the allele-specific record; rsX3 after it is not
        assert_eq!(get(&out[0], "CLIN_SIG"), Some("benign"));
    }

    /// The allele-specific entry comes first and is for this allele: it replaces the list.
    #[test]
    fn an_allele_specific_clin_sig_replaces_the_list_whatever_the_order() {
        let mut a = record("rsW1", "C/A");
        a.clin_sig = Some("pathogenic".into());
        a.clin_sig_allele = Some("A:pathogenic".into());
        a.clin_sig_ref_allele = Some("C".into());
        let mut b = record("rsW2", "C/A");
        b.clin_sig = Some("likely_benign".into());
        let out = cells(
            vec![matched(a, 0, "A"), matched(b, 0, "A")],
            "C",
            &["A"],
            checks(),
        );
        assert_eq!(get(&out[0], "CLIN_SIG"), Some("pathogenic"));
    }

    /// When the record's reference allele is not the input's, the record's alleles
    /// are reverse-complemented before the input allele is looked up in them.
    #[test]
    fn clin_sig_allele_is_reverse_complemented_when_the_reference_alleles_differ() {
        let mut r = record("rs1", "G/T");
        r.strand = -1;
        r.clin_sig_allele = Some("T:pathogenic".into());
        r.clin_sig_ref_allele = Some("G".into());
        // input C>A: the record's G/T is its reverse complement
        let out = cells(vec![matched(r, 0, "T")], "C", &["A"], checks());
        assert_eq!(get(&out[0], "CLIN_SIG"), Some("pathogenic"));
    }

    /// A deletion's `clin_sig_allele` names the deleted bases, which is not `-`,
    /// the allele of the row.
    #[test]
    fn clin_sig_allele_naming_the_deleted_bases_finds_no_value_for_a_deletion() {
        let mut r = record("rs1", "AAGT/-");
        r.clin_sig = Some("pathogenic".into());
        r.clin_sig_allele = Some("AAGT:pathogenic".into());
        r.clin_sig_ref_allele = Some("AAGT".into());
        let out = cells(vec![matched(r, 0, "-")], "AAGT", &["-"], checks());
        assert_eq!(get(&out[0], "CLIN_SIG"), None);
    }

    #[test]
    fn two_distinct_allele_specific_values_are_joined_with_a_semicolon_in_the_order_met() {
        let mut a = record("rsU1", "C/A");
        a.clin_sig_allele = Some("A:pathogenic".into());
        a.clin_sig_ref_allele = Some("C".into());
        let mut b = record("rsU2", "C/A");
        b.clin_sig_allele = Some("A:benign".into());
        b.clin_sig_ref_allele = Some("C".into());
        let out = cells(
            vec![matched(a, 0, "A"), matched(b, 0, "A")],
            "C",
            &["A"],
            checks(),
        );
        assert_eq!(get(&out[0], "CLIN_SIG"), Some("pathogenic;benign"));
    }

    /// SOMATIC and PHENO list one value per known variant, and only when one is set.
    #[test]
    fn somatic_and_pheno_are_lists_that_are_kept_when_any_is_set() {
        let mut a = record("COSV100", "C/A");
        a.somatic = true;
        a.phenotype_or_disease = true;
        let b = record("rsT07a", "C/A");
        let mut c = record("rsT07b", "C/A");
        c.phenotype_or_disease = true;
        let out = cells(
            vec![matched(a, 0, "A"), matched(b, 0, "A"), matched(c, 0, "A")],
            "C",
            &["A"],
            checks(),
        );
        // germline first (rsT07a, rsT07b), then the somatic one
        assert_eq!(get(&out[0], "SOMATIC"), Some("0,0,1"));
        assert_eq!(get(&out[0], "PHENO"), Some("0,1,1"));
    }

    #[test]
    fn somatic_and_pheno_of_all_zeros_are_no_value() {
        let out = cells(
            vec![matched(record("rs1", "C/A"), 0, "A")],
            "C",
            &["A"],
            checks(),
        );
        assert!(out[0].is_empty());
    }

    /// PUBMED is not deduplicated, and is only read when asked for.
    #[test]
    fn pubmed_is_every_id_of_every_known_variant_and_only_when_asked_for() {
        let mut a = record("rsT13a", "C/A");
        a.pubmed = vec!["111".into(), "222".into()];
        let mut b = record("rsT13b", "C/A");
        b.pubmed = vec!["222".into(), "333".into()];
        let records = || vec![matched(a.clone(), 0, "A"), matched(b.clone(), 0, "A")];
        let asked = KnownFlags {
            pubmed: true,
            ..Default::default()
        };
        let out = cells(records(), "C", &["A"], asked);
        assert_eq!(get(&out[0], "PUBMED"), Some("111,222,222,333"));
        let out = cells(records(), "C", &["A"], checks());
        assert_eq!(get(&out[0], "PUBMED"), None);
    }

    /// A known variant that matched only another allele is not reported for this one.
    #[test]
    fn a_known_variant_of_another_allele_adds_nothing() {
        let mut r = record("rs1", "C/G");
        r.clin_sig = Some("benign".into());
        r.phenotype_or_disease = true;
        let out = cells(vec![matched(r, 1, "G")], "C", &["A", "G"], checks());
        assert!(out[0].is_empty());
        assert_eq!(get(&out[1], "CLIN_SIG"), Some("benign"));
        assert_eq!(get(&out[1], "PHENO"), Some("1"));
    }

    /// A record with no alleles is reported for every allele of the site.
    #[test]
    fn a_record_with_no_alleles_is_reported_for_every_allele() {
        let mut r = record("COSV14", "COSMIC_MUTATION");
        r.somatic = true;
        r.phenotype_or_disease = true;
        let out = cells(
            vec![KnownRecord {
                record: Arc::new(r),
                matched: Matched::AllAlleles,
            }],
            "C",
            &["A", "G"],
            checks(),
        );
        for row in &out {
            assert_eq!(get(row, "SOMATIC"), Some("1"));
            assert_eq!(get(row, "PHENO"), Some("1"));
        }
    }

    #[test]
    fn nothing_is_computed_for_a_run_that_asks_for_none_of_it() {
        let site = KnownSite {
            records: vec![matched(
                with_freq(record("rs1", "C/A"), "AF", "A:0.5"),
                0,
                "A",
            )],
        };
        assert!(known_cells(&site, "C", &alleles(&["A"]), KnownFlags::default()).is_empty());
        assert!(known_cells(&KnownSite::default(), "C", &alleles(&["A"]), af()).is_empty());
    }
}
