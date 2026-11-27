//! Extra INFO columns for a built-in VCF source (`sa-build --extra-info`).
//!
//! A built-in source stores a fixed set of columns. `--extra-info` adds INFO keys of
//! the input to that set without a code change. Each entry is `KEY[:ENCODING][=ALIAS]`:
//!
//! - `KEY` is the INFO key exactly as the input spells it. An `##INFO` header line
//!   must declare it, because that line gives the column its type and its number of values.
//! - `ENCODING` says how a `Type=Float` value is stored in an `.osa2` column, and a
//!   Float key must have one: a preset of [`FLOAT_PRESETS`] or a multiplier, with a
//!   trailing `z` for signed values. A `Type=Integer` key takes `z` alone, for signed
//!   values. Other types take none.
//! - `ALIAS` is the JSON key of the value; the INFO key when it is not given.
//!
//! The encoding of a Float is not guessed: a header says nothing about the range or the
//! precision of a Float. A value outside the range of the encoding stops a v2 build, because
//! stored it would read back as 0 or as the top of the range. `docs/SUPPLEMENTARY_ANNOTATIONS.md`
//! lists the encoding to use for common keys.

use crate::fields::{Field, FieldType};
use anyhow::{bail, Context, Result};
use std::cell::RefCell;
use std::collections::{HashMap, HashSet};
use std::io::BufRead;

/// Named encodings for a `Type=Float` key, as (name, multiplier, zigzag, use).
///
/// A value is stored as `value * multiplier`, truncated toward zero, in a u32, so the multiplier sets
/// both the smallest step (`1 / multiplier`) and the largest value stored as itself.
pub const FLOAT_PRESETS: &[(&str, u32, bool, &str)] = &[
    (
        "af",
        2_000_000,
        false,
        "allele frequencies: 0 to 2147, in steps of 5e-7 (the built-in gnomAD AF columns)",
    ),
    (
        "score",
        1_000_000,
        false,
        "scores and probabilities: 0 to 4294, in steps of 1e-6",
    ),
    (
        "signed",
        10_000,
        true,
        "signed values: -214748 to 214748, in steps of 1e-4",
    ),
];

fn preset_names() -> String {
    FLOAT_PRESETS
        .iter()
        .map(|(name, ..)| *name)
        .collect::<Vec<_>>()
        .join(", ")
}

/// One `--extra-info` entry before the header is consulted.
#[derive(Debug, Clone, PartialEq, Eq)]
struct Entry {
    key: String,
    encoding: Option<String>,
    alias: String,
}

fn parse_entry(raw: &str) -> Result<Entry> {
    let raw = raw.trim();
    let (lhs, alias) = match raw.split_once('=') {
        Some((lhs, alias)) => (lhs, Some(alias)),
        None => (raw, None),
    };
    let (key, encoding) = match lhs.split_once(':') {
        Some((key, encoding)) => (key, Some(encoding)),
        None => (lhs, None),
    };
    if key.is_empty() {
        bail!("--extra-info entry '{raw}' names no INFO key");
    }
    if encoding == Some("") {
        bail!("--extra-info {key}: the encoding after ':' is empty");
    }
    let alias = alias.unwrap_or(key);
    // The alias is written into the JSON of every record unescaped, as the
    // built-in aliases are, so it is held to characters that need no escaping.
    if alias.is_empty()
        || !alias
            .chars()
            .all(|c| c.is_ascii_alphanumeric() || matches!(c, '_' | '.' | '-'))
    {
        bail!(
            "--extra-info {key}: the JSON key '{alias}' may hold only letters, digits, '_', '.' and '-'"
        );
    }
    Ok(Entry {
        key: key.to_string(),
        encoding: encoding.map(str::to_string),
        alias: alias.to_string(),
    })
}

/// What an `##INFO` header line declares about a key.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct InfoDecl {
    pub id: String,
    pub number: String,
    pub ty: String,
}

/// Read ID, Number and Type from an `##INFO=<...>` line, in any attribute order.
///
/// A `Description` may hold commas and escaped quotes, so the attributes are split
/// on the commas outside quotes only.
pub fn parse_info_decl(line: &str) -> Option<InfoDecl> {
    let body = line.strip_prefix("##INFO=<")?.strip_suffix('>')?;
    let (mut id, mut number, mut ty) = (None, None, None);
    for (key, value) in header_attributes(body) {
        match key {
            "ID" => id = Some(value),
            "Number" => number = Some(value),
            "Type" => ty = Some(value),
            _ => {}
        }
    }
    Some(InfoDecl {
        id: id?.to_string(),
        number: number?.to_string(),
        ty: ty?.to_string(),
    })
}

fn header_attributes(body: &str) -> Vec<(&str, &str)> {
    let bytes = body.as_bytes();
    let mut out = Vec::new();
    let (mut start, mut i, mut in_quotes) = (0, 0, false);
    while i <= bytes.len() {
        let split = i == bytes.len() || (bytes[i] == b',' && !in_quotes);
        if !split {
            match bytes[i] {
                b'\\' if in_quotes => i += 1,
                b'"' => in_quotes = !in_quotes,
                _ => {}
            }
            i += 1;
            continue;
        }
        if let Some((key, value)) = body[start..i].split_once('=') {
            let value = value
                .strip_prefix('"')
                .and_then(|v| v.strip_suffix('"'))
                .unwrap_or(value);
            out.push((key.trim(), value));
        }
        i += 1;
        start = i;
    }
    out
}

/// The header lines of a VCF, up to the first record.
pub fn read_header_lines<R: BufRead>(reader: R) -> Result<Vec<String>> {
    let mut lines = Vec::new();
    for line in reader.lines() {
        let line = line.context("Reading the VCF header")?;
        if !line.starts_with('#') {
            break;
        }
        lines.push(line);
    }
    Ok(lines)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Kind {
    Integer,
    Float,
    Flag,
    Text,
}

/// Which of the comma-separated values of a key belongs to an allele.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Arity {
    /// `Number=A`: one per alternate allele.
    Allele,
    /// `Number=R`: the reference first, then one per alternate allele.
    AlleleAfterRef,
    /// `Number=1`: one for the site.
    Site,
    /// `Type=Flag`: set or not set.
    Flag,
}

#[derive(Debug)]
struct Column {
    key: String,
    kind: Kind,
    arity: Arity,
    field: Field,
}

impl Column {
    fn new(entry: &Entry, decl: &InfoDecl) -> Result<Self> {
        let key = &entry.key;
        let kind = match decl.ty.as_str() {
            "Integer" => Kind::Integer,
            "Float" => Kind::Float,
            "Flag" => Kind::Flag,
            "String" | "Character" => Kind::Text,
            other => bail!(
                "--extra-info {key}: the header declares Type={other}, which is not a VCF type"
            ),
        };
        let arity = match (kind, decl.number.as_str()) {
            (Kind::Flag, _) => Arity::Flag,
            (_, "A") => Arity::Allele,
            (_, "R") => Arity::AlleleAfterRef,
            (_, "1") => Arity::Site,
            (_, other) => bail!(
                "--extra-info {key}: the header declares Number={other}; a column holds one value \
                 for each allele, so only Number=A, R or 1 can be stored"
            ),
        };
        let (multiplier, zigzag) = match (kind, entry.encoding.as_deref()) {
            (Kind::Float, None) => bail!(
                "--extra-info {key}: a Type=Float key needs an encoding, as in {key}:af. The header \
                 does not say the range or the precision of its values, and a scale that does not \
                 fit them stores wrong values. Presets: {}; or a multiplier, with a trailing z for \
                 signed values (such as 1000z).",
                preset_names()
            ),
            (Kind::Float, Some(encoding)) => parse_float_encoding(key, encoding)?,
            (Kind::Integer, None) => (1, false),
            (Kind::Integer, Some("z")) => (1, true),
            (Kind::Integer, Some(other)) => bail!(
                "--extra-info {key}: a Type=Integer key takes only the encoding z (signed), not '{other}'"
            ),
            (_, None) => (1, false),
            (_, Some(other)) => bail!(
                "--extra-info {key}: a Type={} key takes no encoding, not '{other}'",
                decl.ty
            ),
        };
        let ftype = match kind {
            Kind::Integer => FieldType::Integer,
            Kind::Float => FieldType::Float,
            Kind::Flag => FieldType::Flag,
            Kind::Text => FieldType::Categorical,
        };
        Ok(Self {
            key: key.clone(),
            kind,
            arity,
            field: Field {
                field: key.clone(),
                alias: entry.alias.clone(),
                ftype,
                multiplier,
                zigzag,
                missing_value: u32::MAX,
                missing_string: ".".into(),
                description: format!("INFO/{key}, added with --extra-info"),
            },
        })
    }

    /// The raw value of this key for one alternate allele, or `None` when the
    /// record does not carry it or carries VCF's missing value.
    fn raw<'m>(&self, info_map: &'m HashMap<String, String>, allele_idx: usize) -> Option<&'m str> {
        let value = info_map.get(&self.key)?;
        let one = match self.arity {
            Arity::Allele => value.split(',').nth(allele_idx),
            Arity::AlleleAfterRef => value.split(',').nth(allele_idx + 1),
            Arity::Site | Arity::Flag => Some(value.as_str()),
        };
        one.filter(|v| !v.is_empty() && *v != ".")
    }

    /// Whether `value` lies in the range this column's encoding can store.
    ///
    /// `Field::encode_float` clamps instead of failing, which is right for the built-in
    /// columns whose ranges are known and wrong for one whose range the user guessed: a
    /// negative value without `z` would be stored as 0, and a value above the range as its
    /// top, and a query could not tell either from a measurement.
    fn float_fits(&self, value: f64) -> Result<()> {
        let field = &self.field;
        let scaled = value * f64::from(field.multiplier);
        // Truncation toward zero is the encoding's own step, so (-1, 0) still stores as 0.
        let (low, high) = if field.zigzag {
            (f64::from(i32::MIN) + 1.0, f64::from(i32::MAX) + 1.0)
        } else {
            (0.0, f64::from(u32::MAX))
        };
        if scaled > low - 1.0 && scaled < high {
            return Ok(());
        }
        bail!(
            "--extra-info {}: the value {value} does not fit the encoding (multiplier {}{}), which \
             stores {} to {}. Choose an encoding whose range holds every value of the key.",
            self.key,
            field.multiplier,
            if field.zigzag { ", signed" } else { "" },
            low / f64::from(field.multiplier),
            (high - 1.0) / f64::from(field.multiplier),
        )
    }

    fn int(&self, raw: &str) -> Option<i64> {
        // A value the v2 column cannot hold would be dropped or saturated there; v1
        // leaves it out too, so the formats agree. u32::MAX and its zigzag twin
        // i32::MIN are the missing value.
        let (min, max) = if self.field.zigzag {
            (i64::from(i32::MIN) + 1, i64::from(i32::MAX))
        } else {
            (0, i64::from(u32::MAX) - 1)
        };
        raw.parse::<i64>().ok().filter(|n| (min..=max).contains(n))
    }
}

fn parse_float_encoding(key: &str, encoding: &str) -> Result<(u32, bool)> {
    if let Some((_, multiplier, zigzag, _)) =
        FLOAT_PRESETS.iter().find(|(name, ..)| *name == encoding)
    {
        return Ok((*multiplier, *zigzag));
    }
    let (digits, zigzag) = match encoding.strip_suffix('z') {
        Some(digits) => (digits, true),
        None => (encoding, false),
    };
    match digits.parse::<u32>() {
        Ok(multiplier) if multiplier > 0 => Ok((multiplier, zigzag)),
        _ => bail!(
            "--extra-info {key}: '{encoding}' is not an encoding. Use a preset ({}) or a positive \
             multiplier, with a trailing z for signed values (such as 1000z).",
            preset_names()
        ),
    }
}

/// The strings a `Type=String` column has met, in the order it met them. The
/// `.osa2` stores the index and the archive carries the table.
#[derive(Debug, Default)]
struct Vocabulary {
    index: HashMap<String, u32>,
    strings: Vec<String>,
}

impl Vocabulary {
    fn intern(&mut self, value: &str) -> u32 {
        if let Some(&i) = self.index.get(value) {
            return i;
        }
        let i = self.strings.len() as u32;
        self.index.insert(value.to_string(), i);
        self.strings.push(value.to_string());
        i
    }
}

/// The extra columns of one build, resolved against the input's header.
#[derive(Debug)]
pub struct ExtraInfo {
    columns: Vec<Column>,
    /// One per column; only the `Type=String` ones are ever filled. Interior
    /// mutability, because the iterators that encode records share this with the closure that
    /// hands the tables to the writer after the last record.
    vocabularies: RefCell<Vec<Vocabulary>>,
}

impl ExtraInfo {
    /// Resolve the `--extra-info` entries against the header lines of the input.
    ///
    /// `reserved` holds the JSON keys the source writes itself, which an alias may not
    /// take: a record with one key twice reads as whichever value a JSON parser keeps.
    pub fn resolve(
        entries: &[String],
        header_lines: &[String],
        reserved: &[String],
    ) -> Result<Self> {
        let decls: HashMap<String, InfoDecl> = header_lines
            .iter()
            .filter_map(|line| parse_info_decl(line))
            .map(|decl| (decl.id.clone(), decl))
            .collect();
        let reserved: HashSet<&str> = reserved.iter().map(String::as_str).collect();
        let (mut keys, mut aliases) = (HashSet::new(), HashSet::new());
        let mut columns = Vec::with_capacity(entries.len());
        for raw in entries {
            let entry = parse_entry(raw)?;
            if !keys.insert(entry.key.clone()) {
                bail!("--extra-info {}: the key is listed twice", entry.key);
            }
            if reserved.contains(entry.alias.as_str()) {
                bail!(
                    "--extra-info {}: the source already writes the JSON key '{}'; give another \
                     one, as in {}={}Extra",
                    entry.key,
                    entry.alias,
                    entry.key,
                    entry.alias
                );
            }
            if !aliases.insert(entry.alias.clone()) {
                bail!(
                    "--extra-info: two entries write the JSON key '{}'",
                    entry.alias
                );
            }
            let decl = decls.get(&entry.key).with_context(|| {
                format!(
                    "--extra-info {key}: the input declares no ##INFO=<ID={key}>, so its type and \
                     number are unknown. Keys are matched exactly as the input spells them.",
                    key = entry.key
                )
            })?;
            columns.push(Column::new(&entry, decl)?);
        }
        let vocabularies = (0..columns.len()).map(|_| Vocabulary::default()).collect();
        Ok(Self {
            columns,
            vocabularies: RefCell::new(vocabularies),
        })
    }

    pub fn is_empty(&self) -> bool {
        self.columns.is_empty()
    }

    /// The `.osa2` fields of the extra columns, in the order [`Self::encode`] fills them.
    pub fn fields(&self) -> Vec<Field> {
        self.columns.iter().map(|c| c.field.clone()).collect()
    }

    /// The u32 values of the extra columns for one alternate allele.
    ///
    /// A Float the column's encoding cannot hold is an error, not a clamped value: the
    /// `.osa2` would answer with 0 or the top of the range, which reads as data (see
    /// [`Column::float_fits`]). An Integer the column cannot hold is stored as missing,
    /// as `--extra-info` documents, because that one is visibly absent in every query.
    pub fn encode(
        &self,
        info_map: &HashMap<String, String>,
        allele_idx: usize,
    ) -> Result<Vec<u32>> {
        let mut vocabularies = self.vocabularies.borrow_mut();
        self.columns
            .iter()
            .zip(vocabularies.iter_mut())
            .map(|(column, vocabulary)| {
                let field = &column.field;
                if column.kind == Kind::Flag {
                    return Ok(u32::from(info_map.contains_key(&column.key)));
                }
                let Some(raw) = column.raw(info_map, allele_idx) else {
                    return Ok(field.missing_value);
                };
                Ok(match column.kind {
                    Kind::Integer => column
                        .int(raw)
                        .map_or(field.missing_value, |n| field.encode_int(n)),
                    Kind::Float => match raw.parse::<f64>() {
                        // NaN and infinities are the missing value, as for the built-in columns.
                        Ok(f) if !f.is_finite() => field.missing_value,
                        Ok(f) => {
                            column.float_fits(f)?;
                            field.encode_float(f)
                        }
                        Err(_) => field.missing_value,
                    },
                    Kind::Text => vocabulary.intern(raw),
                    Kind::Flag => unreachable!("a flag returned above"),
                })
            })
            .collect()
    }

    /// The JSON members of the extra columns for one alternate allele, for a v1 `.osa`.
    /// An absent value and an unset flag write nothing, as the built-in columns do.
    pub fn json_members(
        &self,
        info_map: &HashMap<String, String>,
        allele_idx: usize,
        out: &mut Vec<String>,
    ) {
        for column in &self.columns {
            let alias = &column.field.alias;
            if column.kind == Kind::Flag {
                if info_map.contains_key(&column.key) {
                    out.push(format!("\"{alias}\":true"));
                }
                continue;
            }
            let Some(raw) = column.raw(info_map, allele_idx) else {
                continue;
            };
            match column.kind {
                Kind::Integer => {
                    if let Some(n) = column.int(raw) {
                        out.push(format!("\"{alias}\":{n}"));
                    }
                }
                Kind::Float => {
                    if let Some(f) = raw.parse::<f64>().ok().filter(|f| f.is_finite()) {
                        out.push(format!("\"{alias}\":{f:.6e}"));
                    }
                }
                Kind::Text => {
                    let quoted = serde_json::to_string(raw).expect("a string serializes");
                    out.push(format!("\"{alias}\":{quoted}"));
                }
                Kind::Flag => unreachable!("a flag continued above"),
            }
        }
    }

    /// The string tables of the `Type=String` columns, for a writer whose extra columns
    /// start at field index `offset`. Call it after the last record.
    pub fn string_tables(&self, offset: usize) -> Vec<(usize, Vec<String>)> {
        let vocabularies = self.vocabularies.borrow();
        self.columns
            .iter()
            .zip(vocabularies.iter())
            .enumerate()
            .filter(|(_, (column, _))| column.kind == Kind::Text)
            .map(|(i, (_, vocabulary))| (offset + i, vocabulary.strings.clone()))
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const HEADER: &[&str] = &[
        "##fileformat=VCFv4.2",
        r#"##INFO=<ID=AC_nfe,Number=A,Type=Integer,Description="Alternate allele count, nfe">"#,
        r#"##INFO=<ID=AN_nfe,Number=1,Type=Integer,Description="Total allele number, nfe">"#,
        r#"##INFO=<ID=fafmax_faf99_max,Number=A,Type=Float,Description="Max FAF, 99% CI">"#,
        r#"##INFO=<Description="Group of the max FAF, \"grpmax\", by ancestry",Type=String,Number=A,ID=fafmax_faf95_max_gen_anc>"#,
        r#"##INFO=<ID=InbreedingCoeff,Number=A,Type=Float,Description="Inbreeding coefficient">"#,
        r#"##INFO=<ID=AD_ish,Number=R,Type=Integer,Description="Per allele, reference first">"#,
        r#"##INFO=<ID=lcr,Number=0,Type=Flag,Description="Low-complexity region">"#,
        r#"##INFO=<ID=vep,Number=.,Type=String,Description="VEP annotations">"#,
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    ];

    fn header() -> Vec<String> {
        HEADER.iter().map(|s| s.to_string()).collect()
    }

    fn entries(list: &[&str]) -> Vec<String> {
        list.iter().map(|s| s.to_string()).collect()
    }

    fn resolve(list: &[&str]) -> Result<ExtraInfo> {
        ExtraInfo::resolve(&entries(list), &header(), &["allAf".to_string()])
    }

    fn info(pairs: &[(&str, &str)]) -> HashMap<String, String> {
        pairs
            .iter()
            .map(|(k, v)| (k.to_string(), v.to_string()))
            .collect()
    }

    fn json(
        extra: &ExtraInfo,
        info_map: &HashMap<String, String>,
        allele_idx: usize,
    ) -> serde_json::Value {
        let mut members = Vec::new();
        extra.json_members(info_map, allele_idx, &mut members);
        serde_json::from_str(&format!("{{{}}}", members.join(","))).unwrap()
    }

    #[test]
    fn test_header_declaration_is_read_in_any_order_past_quoted_commas() {
        let decl = parse_info_decl(HEADER[4]).unwrap();
        assert_eq!(decl.id, "fafmax_faf95_max_gen_anc");
        assert_eq!(decl.number, "A");
        assert_eq!(decl.ty, "String");
    }

    #[test]
    fn test_header_lines_stop_at_the_first_record() {
        let vcf = "##fileformat=VCFv4.2\n#CHROM\tPOS\n1\t100\n##not_header\n";
        let lines = read_header_lines(vcf.as_bytes()).unwrap();
        assert_eq!(lines, vec!["##fileformat=VCFv4.2", "#CHROM\tPOS"]);
    }

    #[test]
    fn test_type_and_number_come_from_the_header() {
        let extra = resolve(&[
            "AC_nfe=nfeAc",
            "AN_nfe=nfeAn",
            "fafmax_faf99_max:af=faf99Max",
            "fafmax_faf95_max_gen_anc=faf95MaxGenAnc",
            "lcr",
        ])
        .unwrap();
        let fields = extra.fields();
        let types: Vec<_> = fields.iter().map(|f| (f.alias.as_str(), f.ftype)).collect();
        assert_eq!(
            types,
            vec![
                ("nfeAc", FieldType::Integer),
                ("nfeAn", FieldType::Integer),
                ("faf99Max", FieldType::Float),
                ("faf95MaxGenAnc", FieldType::Categorical),
                ("lcr", FieldType::Flag),
            ]
        );
        assert_eq!(fields[2].multiplier, 2_000_000);
        assert!(!fields[2].zigzag);
    }

    #[test]
    fn test_a_float_key_without_an_encoding_is_refused() {
        let err = resolve(&["fafmax_faf99_max"]).unwrap_err().to_string();
        assert!(err.contains("needs an encoding"), "{err}");
        assert!(err.contains("af, score, signed"), "{err}");
    }

    #[test]
    fn test_float_encodings_are_presets_or_multipliers() {
        assert_eq!(parse_float_encoding("k", "af").unwrap(), (2_000_000, false));
        assert_eq!(
            parse_float_encoding("k", "score").unwrap(),
            (1_000_000, false)
        );
        assert_eq!(parse_float_encoding("k", "signed").unwrap(), (10_000, true));
        assert_eq!(parse_float_encoding("k", "1000").unwrap(), (1000, false));
        assert_eq!(parse_float_encoding("k", "1000z").unwrap(), (1000, true));
        for bad in ["0", "z", "-5", "1e6", "afz", "x"] {
            assert!(
                parse_float_encoding("k", bad).is_err(),
                "{bad} was accepted"
            );
        }
    }

    #[test]
    fn test_an_integer_key_takes_only_z() {
        let extra = resolve(&["AC_nfe:z"]).unwrap();
        assert!(extra.fields()[0].zigzag);
        assert!(resolve(&["AC_nfe:af"]).is_err());
        assert!(resolve(&["lcr:af"]).is_err());
        assert!(resolve(&["fafmax_faf95_max_gen_anc:z"]).is_err());
    }

    #[test]
    fn test_keys_the_input_does_not_declare_are_refused() {
        let err = resolve(&["AC_joint_nfe"]).unwrap_err().to_string();
        assert!(
            err.contains("declares no ##INFO=<ID=AC_joint_nfe>"),
            "{err}"
        );
    }

    #[test]
    fn test_a_number_a_column_cannot_hold_is_refused() {
        let err = resolve(&["vep"]).unwrap_err().to_string();
        assert!(err.contains("Number=."), "{err}");
    }

    #[test]
    fn test_duplicate_keys_and_aliases_are_refused() {
        assert!(resolve(&["AC_nfe", "AC_nfe=other"]).is_err());
        assert!(resolve(&["AC_nfe=n", "AN_nfe=n"]).is_err());
        let err = resolve(&["AC_nfe=allAf"]).unwrap_err().to_string();
        assert!(err.contains("already writes the JSON key 'allAf'"), "{err}");
    }

    #[test]
    fn test_an_alias_needing_json_escapes_is_refused() {
        assert!(resolve(&["AC_nfe=a\"b"]).is_err());
        assert!(resolve(&["AC_nfe=a b"]).is_err());
        assert!(resolve(&["AC_nfe="]).is_err());
        assert!(resolve(&[":af"]).is_err());
    }

    #[test]
    fn test_values_are_sliced_per_allele_by_number() {
        let extra = resolve(&["AC_nfe", "AN_nfe", "AD_ish"]).unwrap();
        let record = info(&[
            ("AC_nfe", "30,7"),
            ("AN_nfe", "60000"),
            ("AD_ish", "100,30,7"),
        ]);
        assert_eq!(extra.encode(&record, 0).unwrap(), vec![30, 60000, 30]);
        assert_eq!(extra.encode(&record, 1).unwrap(), vec![7, 60000, 7]);
    }

    #[test]
    fn test_absent_and_dot_values_are_missing() {
        let extra = resolve(&["AC_nfe", "fafmax_faf99_max:af", "lcr"]).unwrap();
        let record = info(&[("AC_nfe", ".")]);
        assert_eq!(
            extra.encode(&record, 0).unwrap(),
            vec![u32::MAX, u32::MAX, 0]
        );
        assert_eq!(json(&extra, &record, 0), serde_json::json!({}));
    }

    #[test]
    fn test_v1_and_v2_carry_the_same_values() {
        let extra = resolve(&[
            "AC_nfe=nfeAc",
            "fafmax_faf99_max:af=faf99Max",
            "InbreedingCoeff:signed=inbreeding",
            "fafmax_faf95_max_gen_anc=faf95MaxGenAnc",
            "lcr",
        ])
        .unwrap();
        let record = info(&[
            ("AC_nfe", "30"),
            ("fafmax_faf99_max", "0.0003"),
            ("InbreedingCoeff", "-0.25"),
            ("fafmax_faf95_max_gen_anc", "nfe"),
            ("lcr", ""),
        ]);
        let v = json(&extra, &record, 0);
        assert_eq!(v["nfeAc"], 30);
        assert!((v["faf99Max"].as_f64().unwrap() - 0.0003).abs() < 1e-12);
        assert!((v["inbreeding"].as_f64().unwrap() + 0.25).abs() < 1e-12);
        assert_eq!(v["faf95MaxGenAnc"], "nfe");
        assert_eq!(v["lcr"], true);

        let values = extra.encode(&record, 0).unwrap();
        let fields = extra.fields();
        assert_eq!(values[0], 30);
        assert!((fields[1].decode_float(values[1]) - 0.0003).abs() < 1e-6);
        assert!((fields[2].decode_float(values[2]) + 0.25).abs() < 1e-4);
        assert_eq!(values[4], 1);
        let tables = extra.string_tables(50);
        assert_eq!(tables, vec![(53, vec!["nfe".to_string()])]);
        assert_eq!(values[3], 0);
    }

    #[test]
    fn test_string_values_share_one_table() {
        let extra = resolve(&["fafmax_faf95_max_gen_anc"]).unwrap();
        let codes: Vec<u32> = ["nfe", "afr", "nfe"]
            .iter()
            .map(|g| {
                extra
                    .encode(&info(&[("fafmax_faf95_max_gen_anc", g)]), 0)
                    .unwrap()[0]
            })
            .collect();
        assert_eq!(codes, vec![0, 1, 0]);
        assert_eq!(
            extra.string_tables(0),
            vec![(0, vec!["nfe".to_string(), "afr".to_string()])]
        );
    }

    #[test]
    fn test_a_negative_integer_is_missing_unless_signed() {
        let unsigned = resolve(&["AC_nfe"]).unwrap();
        let signed = resolve(&["AC_nfe:z"]).unwrap();
        let record = info(&[("AC_nfe", "-3")]);
        assert_eq!(unsigned.encode(&record, 0).unwrap(), vec![u32::MAX]);
        assert_eq!(json(&unsigned, &record, 0), serde_json::json!({}));
        let stored = signed.encode(&record, 0).unwrap()[0];
        assert_eq!(signed.fields()[0].decode_int(stored), -3);
        assert_eq!(json(&signed, &record, 0)["AC_nfe"], -3);
    }

    #[test]
    fn test_an_integer_the_column_cannot_hold_is_missing_in_both_formats() {
        let unsigned = resolve(&["AC_nfe"]).unwrap();
        let signed = resolve(&["AC_nfe:z"]).unwrap();
        for (extra, value) in [
            (&unsigned, "4294967295"),
            (&unsigned, "5000000000"),
            (&signed, "3000000000"),
            (&signed, "-2147483648"),
        ] {
            let record = info(&[("AC_nfe", value)]);
            assert_eq!(extra.encode(&record, 0).unwrap(), vec![u32::MAX], "{value}");
            assert_eq!(json(extra, &record, 0), serde_json::json!({}), "{value}");
        }
        let top = info(&[("AC_nfe", "4294967294")]);
        assert_eq!(json(&unsigned, &top, 0)["AC_nfe"], 4294967294_u64);
        assert_eq!(unsigned.encode(&top, 0).unwrap(), vec![u32::MAX - 1]);
    }

    #[test]
    fn test_a_float_the_encoding_cannot_hold_stops_the_build() {
        let unsigned = resolve(&["InbreedingCoeff:af"]).unwrap();
        let signed = resolve(&["InbreedingCoeff:signed"]).unwrap();
        let encode = |extra: &ExtraInfo, v: &str| extra.encode(&info(&[("InbreedingCoeff", v)]), 0);

        // A negative value would be stored as 0, and a huge one as the top of the range.
        let err = encode(&unsigned, "-0.25").unwrap_err().to_string();
        assert!(
            err.contains("InbreedingCoeff") && err.contains("-0.25"),
            "{err}"
        );
        assert!(encode(&unsigned, "5000").is_err());
        assert!(encode(&signed, "300000").is_err());
        assert!(encode(&signed, "-300000").is_err());

        // What fits is stored, including values below one step, which truncate to 0.
        assert!(encode(&unsigned, "0.0000001").is_ok());
        assert!(encode(&unsigned, "-0.0000001").is_ok());
        assert!(encode(&signed, "-0.25").is_ok());
        // Missing values are not out of range.
        assert_eq!(encode(&unsigned, ".").unwrap(), vec![u32::MAX]);
        assert_eq!(encode(&unsigned, "nan").unwrap(), vec![u32::MAX]);
    }

    #[test]
    fn test_a_number_1_value_is_read_whole() {
        // VCF does not permit a comma in a Number=1 value; one that has it is kept whole
        // rather than cut at the comma with no sign.
        let extra = ExtraInfo::resolve(
            &entries(&["site_text"]),
            &[r#"##INFO=<ID=site_text,Number=1,Type=String,Description="x">"#.to_string()],
            &[],
        )
        .unwrap();
        let record = info(&[("site_text", "a,b")]);
        assert_eq!(json(&extra, &record, 1)["site_text"], "a,b");
    }

    #[test]
    fn test_a_backslash_at_the_end_of_a_quoted_value_keeps_the_attributes() {
        let decl = parse_info_decl(r#"##INFO=<ID=k,Number=1,Type=String,Description="ends in ">"#);
        let decl = decl.expect("declaration");
        assert_eq!(
            (decl.id.as_str(), decl.number.as_str(), decl.ty.as_str()),
            ("k", "1", "String")
        );
        let tail = header_attributes(r#"ID=k,Description="x"#);
        assert_eq!(tail.len(), 2, "{tail:?}");
    }

    #[test]
    fn test_a_text_value_is_escaped_in_json() {
        let extra = resolve(&["fafmax_faf95_max_gen_anc"]).unwrap();
        let record = info(&[("fafmax_faf95_max_gen_anc", "a\"b\\c")]);
        assert_eq!(
            json(&extra, &record, 0)["fafmax_faf95_max_gen_anc"],
            "a\"b\\c"
        );
    }
}
