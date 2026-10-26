use crate::variant::{AlleleAnnotation, TranscriptVariation, VariationFeature};
use fastvep_core::{so_variation_class, Allele, Consequence, Impact, Ordinals, TranscriptSpan};
use serde_json::Value;
use std::fmt::Write as FmtWrite;

/// Format a VCF CSQ INFO field value from a VariationFeature.
///
/// Fields match the standard VEP CSQ format:
/// Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|
/// EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|
/// Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS
pub fn format_csq(vf: &VariationFeature, fields: &[&str]) -> String {
    let mut result = String::with_capacity(1024);

    // Resolve the column names once for the whole variant. A variant in a
    // gene-dense window produces dozens of entries, and every one of them used
    // to re-match all 49 header strings.
    let fields: Vec<CsqField> = fields.iter().copied().map(CsqField::from_name).collect();

    let mut first = true;
    for tv in &vf.transcript_variations {
        for aa in &tv.allele_annotations {
            if !first {
                result.push(',');
            }
            first = false;
            format_csq_entry_into(vf, tv, aa, &fields, &mut result);
        }
    }

    result
}

/// One CSQ column, resolved from its header name.
///
/// The per-entry writer used to `match` the column's header string, and it ran
/// once per column per entry - 49 string matches for every
/// (variant x transcript x allele). A profile of the annotate path put 7% of
/// the program's instructions in `slice::cmp` underneath it. The names are
/// fixed for a whole run, so they are resolved once and the writer dispatches
/// on a discriminant.
///
/// `Unknown` covers a header fastVEP does not populate - `MOTIF_NAME` and the
/// rest of the regulatory columns are in [`DEFAULT_CSQ_FIELDS`] but have no
/// writer - and renders as an empty column, which is what the old catch-all
/// arm did. Naming it means the writer's `match` is exhaustive, so a column
/// added here without a writer is a compile error rather than a silently blank
/// field.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum CsqField {
    Allele,
    Consequence,
    Impact,
    Symbol,
    Gene,
    FeatureType,
    Feature,
    Biotype,
    Exon,
    Intron,
    Hgvsg,
    Hgvsc,
    Hgvsp,
    CdnaPosition,
    CdsPosition,
    ProteinPosition,
    AminoAcids,
    Codons,
    ExistingVariation,
    RefAllele,
    UploadedAllele,
    Distance,
    Strand,
    Flags,
    Pick,
    Canonical,
    SymbolSource,
    HgncId,
    Mane,
    ManeSelect,
    ManePlusClinical,
    Tsl,
    Appris,
    Ccds,
    GencodePrimary,
    Ensp,
    Sift,
    PolyPhen,
    Af,
    ClinSig,
    Somatic,
    Pheno,
    Pubmed,
    Source,
    HgvsOffset,
    VariantClass,
    Acmg,
    AcmgCriteria,
    Unknown,
}

impl CsqField {
    fn from_name(name: &str) -> Self {
        match name {
            "Allele" => Self::Allele,
            "Consequence" => Self::Consequence,
            "IMPACT" => Self::Impact,
            "SYMBOL" => Self::Symbol,
            "Gene" => Self::Gene,
            "Feature_type" => Self::FeatureType,
            "Feature" => Self::Feature,
            "BIOTYPE" => Self::Biotype,
            "EXON" => Self::Exon,
            "INTRON" => Self::Intron,
            "HGVSg" => Self::Hgvsg,
            "HGVSc" => Self::Hgvsc,
            "HGVSp" => Self::Hgvsp,
            "cDNA_position" => Self::CdnaPosition,
            "CDS_position" => Self::CdsPosition,
            "Protein_position" => Self::ProteinPosition,
            "Amino_acids" => Self::AminoAcids,
            "Codons" => Self::Codons,
            "Existing_variation" => Self::ExistingVariation,
            "REF_ALLELE" => Self::RefAllele,
            "UPLOADED_ALLELE" => Self::UploadedAllele,
            "DISTANCE" => Self::Distance,
            "STRAND" => Self::Strand,
            "FLAGS" => Self::Flags,
            "PICK" => Self::Pick,
            "CANONICAL" => Self::Canonical,
            "SYMBOL_SOURCE" => Self::SymbolSource,
            "HGNC_ID" => Self::HgncId,
            "MANE" => Self::Mane,
            "MANE_SELECT" => Self::ManeSelect,
            "MANE_PLUS_CLINICAL" => Self::ManePlusClinical,
            "TSL" => Self::Tsl,
            "APPRIS" => Self::Appris,
            "CCDS" => Self::Ccds,
            "GENCODE_PRIMARY" => Self::GencodePrimary,
            "ENSP" => Self::Ensp,
            "SIFT" => Self::Sift,
            "PolyPhen" => Self::PolyPhen,
            "AF" => Self::Af,
            "CLIN_SIG" => Self::ClinSig,
            "SOMATIC" => Self::Somatic,
            "PHENO" => Self::Pheno,
            "PUBMED" => Self::Pubmed,
            "SOURCE" => Self::Source,
            "HGVS_OFFSET" => Self::HgvsOffset,
            "VARIANT_CLASS" => Self::VariantClass,
            "ACMG" => Self::Acmg,
            "ACMG_CRITERIA" => Self::AcmgCriteria,
            _ => Self::Unknown,
        }
    }
}

/// Write a single CSQ entry directly into the output buffer, avoiding intermediate String allocations.
fn format_csq_entry_into(
    vf: &VariationFeature,
    tv: &TranscriptVariation,
    aa: &AlleleAnnotation,
    fields: &[CsqField],
    buf: &mut String,
) {
    for (i, field) in fields.iter().enumerate() {
        if i > 0 {
            buf.push('|');
        }
        // Write each field value directly into buf via escape_csq_str, avoiding
        // temporary String allocations for most fields.
        match *field {
            CsqField::Allele => escape_csq_str(&aa.allele.to_string(), buf),
            CsqField::Consequence => {
                for (j, c) in aa.consequences.iter().enumerate() {
                    if j > 0 {
                        buf.push('&');
                    }
                    buf.push_str(c.so_term());
                }
            }
            CsqField::Impact => buf.push_str(aa.impact.as_str()),
            CsqField::Symbol => escape_csq_str(tv.gene_symbol.as_deref().unwrap_or_default(), buf),
            // A row with no feature (an intergenic variant) leaves these empty, as
            // VEP does; the placeholder's own `-` and `Transcript` are not values.
            CsqField::Gene if !tv.is_featureless() => escape_csq_str(&tv.gene_id, buf),
            CsqField::FeatureType if !tv.is_featureless() => buf.push_str("Transcript"),
            CsqField::Feature if !tv.is_featureless() => escape_csq_str(&tv.transcript_id, buf),
            CsqField::Biotype if !tv.is_featureless() => escape_csq_str(&tv.biotype, buf),
            CsqField::Gene | CsqField::FeatureType | CsqField::Feature | CsqField::Biotype => {}
            CsqField::Exon => {
                if let Some(o) = aa.location.exon {
                    o.write_to(buf);
                }
            }
            CsqField::Intron => {
                if let Some(o) = aa.location.intron {
                    o.write_to(buf);
                }
            }
            CsqField::Hgvsg => escape_csq_str(aa.hgvsg.as_deref().unwrap_or_default(), buf),
            CsqField::Hgvsc => escape_csq_str(aa.hgvsc.as_deref().unwrap_or_default(), buf),
            CsqField::Hgvsp => escape_csq_str(aa.hgvsp.as_deref().unwrap_or_default(), buf),
            CsqField::CdnaPosition => aa.location.cdna.write_to(buf),
            CsqField::CdsPosition => aa.location.cds.write_to(buf),
            CsqField::ProteinPosition => aa.location.protein.write_to(buf),
            CsqField::AminoAcids => {
                if let Some((ref r, ref a)) = aa.amino_acids {
                    escape_csq_str(r, buf);
                    if r != a {
                        buf.push('/');
                        escape_csq_str(a, buf);
                    }
                }
            }
            CsqField::Codons => {
                if let Some((ref r, ref a)) = aa.codons {
                    escape_csq_str(r, buf);
                    buf.push('/');
                    escape_csq_str(a, buf);
                }
            }
            CsqField::ExistingVariation => {
                for (j, ev) in aa.existing_variation.iter().enumerate() {
                    if j > 0 {
                        buf.push('&');
                    }
                    escape_csq_str(ev, buf);
                }
            }
            // VEP's VCF writer blanks a `-` in every field but Allele
            // (`output_hash_to_vcf_info_chunk`), so an insertion has no reference
            // allele here, where the tab file writes `-`.
            CsqField::RefAllele => {
                let reference = vf.ref_allele.to_string();
                if reference != "-" {
                    escape_csq_str(&reference, buf);
                }
            }
            CsqField::UploadedAllele => {
                if let Some(ref vcf) = vf.vcf_fields {
                    escape_csq_str(&vcf.ref_allele, buf);
                    buf.push('/');
                    escape_csq_str(&vcf.alt, buf);
                } else {
                    escape_csq_str(&vf.ref_allele.to_string(), buf);
                    buf.push('/');
                    escape_csq_str(&aa.allele.to_string(), buf);
                }
            }
            CsqField::Distance => {
                if let Some(d) = aa.distance {
                    let _ = write!(buf, "{}", d);
                }
            }
            CsqField::Strand => {
                if !tv.is_featureless() {
                    let _ = write!(buf, "{}", tv.strand.as_int());
                }
            }
            CsqField::Flags => {
                for (j, f) in tv.flags.iter().enumerate() {
                    if j > 0 {
                        buf.push('&');
                    }
                    buf.push_str(f);
                }
            }
            // `1` or empty, as VEP writes it - measured on 115.1 with
            // `--flag_pick_allele_gene`, which flags one entry per (allele,
            // gene) and leaves the column empty on the rest.
            CsqField::Pick => {
                if aa.pick {
                    buf.push('1');
                }
            }
            CsqField::Canonical => {
                if tv.canonical {
                    buf.push_str("YES");
                }
            }
            CsqField::SymbolSource => {
                escape_csq_str(tv.symbol_source.as_deref().unwrap_or_default(), buf)
            }
            CsqField::HgncId => escape_csq_str(tv.hgnc_id.as_deref().unwrap_or_default(), buf),
            CsqField::Mane => {
                if tv.mane_select.is_some() {
                    buf.push_str("MANE_Select");
                } else if tv.mane_plus_clinical.is_some() {
                    buf.push_str("MANE_Plus_Clinical");
                }
            }
            CsqField::ManeSelect => {
                escape_csq_str(tv.mane_select.as_deref().unwrap_or_default(), buf)
            }
            CsqField::ManePlusClinical => {
                escape_csq_str(tv.mane_plus_clinical.as_deref().unwrap_or_default(), buf)
            }
            CsqField::Tsl => {
                if let Some(t) = tv.tsl {
                    let _ = write!(buf, "{}", t);
                }
            }
            CsqField::Appris => escape_csq_str(tv.appris.as_deref().unwrap_or_default(), buf),
            CsqField::Ccds => escape_csq_str(tv.ccds.as_deref().unwrap_or_default(), buf),
            CsqField::GencodePrimary => {
                if tv.gencode_primary {
                    buf.push_str("YES");
                }
            }
            CsqField::Ensp => escape_csq_str(tv.protein_id.as_deref().unwrap_or_default(), buf),
            CsqField::Sift => escape_csq_str(aa.sift.as_deref().unwrap_or_default(), buf),
            CsqField::PolyPhen => escape_csq_str(aa.polyphen.as_deref().unwrap_or_default(), buf),
            CsqField::Af => {
                if let Some(f) = vf.existing_variants.iter().find_map(|kv| {
                    kv.frequencies
                        .get("gnomAD")
                        .or_else(|| kv.frequencies.get("gnomADe"))
                        .or_else(|| kv.frequencies.get("minor_allele_freq"))
                }) {
                    let _ = write!(buf, "{}", f);
                }
            }
            CsqField::ClinSig => {
                if let Some(cs) = vf
                    .existing_variants
                    .iter()
                    .find_map(|kv| kv.clinical_significance.as_deref())
                {
                    escape_csq_str(cs, buf);
                }
            }
            CsqField::Somatic => {
                if vf.existing_variants.iter().any(|kv| kv.somatic) {
                    buf.push('1');
                }
            }
            CsqField::Pheno => {
                if vf
                    .existing_variants
                    .iter()
                    .any(|kv| kv.phenotype_or_disease)
                {
                    buf.push('1');
                }
            }
            CsqField::Pubmed => {
                let mut first_pub = true;
                for kv in &vf.existing_variants {
                    for p in &kv.pubmed {
                        if !first_pub {
                            buf.push('&');
                        }
                        first_pub = false;
                        buf.push_str(p);
                    }
                }
            }
            CsqField::Source => escape_csq_str(tv.source.as_deref().unwrap_or_default(), buf),
            CsqField::HgvsOffset => {
                if let Some(o) = aa.hgvs_offset {
                    let _ = write!(buf, "{}", o);
                }
            }
            CsqField::VariantClass => {
                // a structural variant's class comes from the parser in VEP
                if !vf.variant_type.is_structural() {
                    buf.push_str(so_variation_class(&vf.allele_string));
                }
            }
            CsqField::Acmg => {
                if let Some(ref acmg) = aa.acmg_classification {
                    if let Some(sh) = acmg.get("shorthand").and_then(|v| v.as_str()) {
                        buf.push_str(sh);
                    }
                }
            }
            CsqField::AcmgCriteria => {
                if let Some(ref acmg) = aa.acmg_classification {
                    if let Some(criteria) = acmg.get("criteria").and_then(|v| v.as_array()) {
                        let mut first = true;
                        for c in criteria {
                            if c.get("met").and_then(|v| v.as_bool()).unwrap_or(false) {
                                if let Some(code) = c.get("code").and_then(|v| v.as_str()) {
                                    if !first {
                                        buf.push('&');
                                    }
                                    first = false;
                                    buf.push_str(code);
                                }
                            }
                        }
                    }
                }
            }
            CsqField::Unknown => {}
        }
    }
}

/// Escape special characters in CSQ field values, appending to an existing
/// buffer.
///
/// Five of these substitutions are the ones Ensembl VEP makes, measured
/// against 115.1 by `validation/run_escaping_probe.sh` rather than assumed:
/// `,` and `|` become `&`, `;` becomes `%3B`, `=` becomes `%3D` (VEP writes a
/// synonymous HGVSp as `ENSP00000269305.1:p.Thr125%3D`), and a space becomes
/// `_`. Note that VEP is not self-consistent on `=` - the value it appends to
/// CSQ for a `--custom` annotation keeps a literal one - so this follows what
/// it does to CSQ's own fields, which is also the only reading VCF 4.3
/// permits.
///
/// TAB / CR / LF are encoded without a VEP measurement behind them; see the
/// comment at that arm.
///
/// `:` is *not* escaped, by VEP or here: every HGVSc carries one.
fn escape_csq_str(value: &str, buf: &mut String) {
    // Nearly every value that reaches here - a transcript ID, a gene symbol,
    // an SO term, an HGVS string, a run of bases - contains none of the five
    // characters that need substituting, and copying such a value whole skips
    // a UTF-8 encode and a capacity check per character. All five are ASCII,
    // so scanning bytes decides this exactly; only a value that really needs
    // escaping falls through to the character loop.
    if !value
        .bytes()
        .any(|b| matches!(b, b',' | b'|' | b';' | b'=' | b' ' | b'\t' | b'\r' | b'\n'))
    {
        buf.push_str(value);
        return;
    }
    for c in value.chars() {
        match c {
            ',' | '|' => buf.push('&'),
            ';' => buf.push_str("%3B"),
            '=' => buf.push_str("%3D"),
            // VEP substitutes rather than encodes, so a reader cannot recover
            // an original `_` from it either. Matching VEP is worth more than
            // reversibility in a column whose whole point is to be diffable
            // against VEP's, and no CSQ field carries free text.
            ' ' => buf.push('_'),
            // Not measured against VEP, because no source could deliver one
            // to measure with: a tab or a newline here does not mangle a
            // value, it ends the column or the record, and there is no
            // reading of "agrees with VEP" that justifies writing a VCF a
            // reader will silently truncate. `push_escaped_subfield_char`
            // encodes the same three for the same reason.
            '\t' => buf.push_str("%09"),
            '\r' => buf.push_str("%0D"),
            '\n' => buf.push_str("%0A"),
            _ => buf.push(c),
        }
    }
}

/// `escape_csq_str` as an expression, for the tests that assert one value at a
/// time.
///
/// This used to be a second implementation with its own `replace` chain, and
/// it had drifted: it substituted `_` for a space and the writer above did
/// not, so `test_escape_csq_value` asserted VEP's behaviour against a function
/// no output path called. A wrapper cannot drift.
#[cfg(test)]
fn escape_csq_value(value: &str) -> String {
    let mut buf = String::with_capacity(value.len());
    escape_csq_str(value, &mut buf);
    buf
}

/// A position cell of the tab format. VEP writes `-` for every value it does
/// not have (`OutputFactory::Tab::output_hash_to_line`), so a span with neither
/// end is `-` here and not the empty string the CSQ field uses.
fn span_cell(span: TranscriptSpan) -> String {
    let mut cell = String::new();
    span.write_to(&mut cell);
    if cell.is_empty() {
        cell.push('-');
    }
    cell
}

/// VEP's `Uploaded_variation`: the first ID the record carries, otherwise a name
/// built from the input line itself, `{chr}_{pos}_{ref}/{alt}`.
///
/// For VCF input that is the line as it was written: the POS column, REF, and
/// ALT with the commas of several alts turned into slashes. `chr10_150_CA/C` for
/// `chr10 150 . CA C`, and `chr10_109_CA/C/TA` for `CA C,TA`, not the minimal
/// `chr10_151_A/-` of the annotation (`OutputFactory::VariationFeature_to_output_hash`,
/// VEP 116). A reader that joins on this name, as the reporter of #139 does, gets
/// the key the caller's own file gives it.
///
/// VEP 115 mixed the two, the trimmed start with the original alleles
/// (`chr10_151_CA/C`), and 116 is the release this follows.
///
/// Input that did not come from a VCF line has only the annotation to name it:
/// `{chr}_{start}_{allele_string}`.
pub fn uploaded_variation(vf: &VariationFeature) -> String {
    let id = vf
        .variation_name
        .as_deref()
        .and_then(|names| names.split(';').next())
        .filter(|id| !id.is_empty() && *id != ".");
    if let Some(id) = id {
        return id.to_string();
    }
    match &vf.vcf_fields {
        Some(line) => {
            let mut name = String::with_capacity(
                line.chrom.len() + line.ref_allele.len() + line.alt.len() + 16,
            );
            name.push_str(&line.chrom);
            name.push('_');
            name.push_str(&line.pos.to_string());
            name.push('_');
            name.push_str(&line.ref_allele);
            for alt in line.alt.split(',') {
                name.push('/');
                name.push_str(alt);
            }
            name
        }
        None => format!(
            "{}_{}_{}",
            vf.position.chromosome, vf.position.start, vf.allele_string
        ),
    }
}

/// VEP's `pep_allele_string`: `L/P`, or just `L` when the protein does not
/// change. The VCF `CSQ` writer already did this; the tab and JSON writers wrote
/// `L/L`.
fn amino_acids_text(reference: &str, alternate: &str) -> String {
    if reference == alternate {
        reference.to_string()
    } else {
        format!("{reference}/{alternate}")
    }
}

/// Default CSQ fields matching Ensembl VEP's extended output format.
///
/// Includes all standard VEP fields (CANONICAL, CCDS, ENSP, SOURCE,
/// HGVS_OFFSET) plus extended annotations (MANE, SIFT, PolyPhen, etc.).
pub const DEFAULT_CSQ_FIELDS: &[&str] = &[
    "Allele",
    "Consequence",
    "IMPACT",
    "SYMBOL",
    "Gene",
    "Feature_type",
    "Feature",
    "BIOTYPE",
    "EXON",
    "INTRON",
    "HGVSc",
    "HGVSp",
    "cDNA_position",
    "CDS_position",
    "Protein_position",
    "Amino_acids",
    "Codons",
    "Existing_variation",
    "REF_ALLELE",
    "UPLOADED_ALLELE",
    "DISTANCE",
    "STRAND",
    "FLAGS",
    "CANONICAL",
    "SYMBOL_SOURCE",
    "HGNC_ID",
    "MANE",
    "MANE_SELECT",
    "MANE_PLUS_CLINICAL",
    "TSL",
    "APPRIS",
    "CCDS",
    "ENSP",
    "SOURCE",
    "HGVS_OFFSET",
    "SIFT",
    "PolyPhen",
    "AF",
    "CLIN_SIG",
    "SOMATIC",
    "PHENO",
    "PUBMED",
    "MOTIF_NAME",
    "MOTIF_POS",
    "HIGH_INF_POS",
    "MOTIF_SCORE_CHANGE",
    "TRANSCRIPTION_FACTORS",
    "ACMG",
    "ACMG_CRITERIA",
];

/// [`DEFAULT_CSQ_FIELDS`] with `PICK` spliced in after `FLAGS`, for a run that
/// flags a pick instead of reducing to it.
///
/// A column that is always present but always empty is worse than no column:
/// it widens every record of every run for a flag almost nobody passes, and it
/// tells a reader nothing about whether a pick was even requested. VEP adds
/// the column only under `--flag_pick*` too, and puts it in this position.
pub fn csq_fields_with_pick() -> Vec<&'static str> {
    csq_fields(true, ColumnFlags::default())
}

/// The CSQ field list of a run: [`DEFAULT_CSQ_FIELDS`], plus what VEP writes only
/// on request. `PICK` goes after `FLAGS` for a flagging run (see above), then
/// `VARIANT_CLASS` for `--variant-class`; `HGVSg` goes after `HGVS_OFFSET` for
/// `--hgvsg`. That is VEP's relative order, and with neither flag the list is
/// the 49 fields every reader of the CSQ already parses.
pub fn csq_fields(pick: bool, flags: ColumnFlags) -> Vec<&'static str> {
    let mut fields = Vec::with_capacity(DEFAULT_CSQ_FIELDS.len() + 3);
    for &name in DEFAULT_CSQ_FIELDS {
        fields.push(name);
        match name {
            "FLAGS" => {
                if pick {
                    fields.push("PICK");
                }
                if flags.variant_class {
                    fields.push("VARIANT_CLASS");
                }
            }
            "HGVS_OFFSET" if flags.hgvsg => fields.push("HGVSg"),
            _ => {}
        }
    }
    fields
}

/// Generate the VCF INFO header line for CSQ.
pub fn csq_header_line(fields: &[&str]) -> String {
    format!(
        "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from fastVEP. Format: {}\">",
        fields.join("|")
    )
}

/// Description string shared by SpliceAI's VCF `##INFO=` header and the
/// tab `## COLUMN=<...>` prologue so both formats document the same schema.
pub const SPLICEAI_DESCRIPTION: &str =
    "SpliceAI annotations. Format: ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL";

/// Generate the VCF INFO header line for SpliceAI annotations emitted from fastSA.
pub fn spliceai_header_line() -> String {
    format!(
        "##INFO=<ID=SpliceAI,Number=.,Type=String,Description=\"{}\">",
        SPLICEAI_DESCRIPTION
    )
}

#[derive(Clone, Copy)]
struct VcfProjectionSpec {
    json_key: &'static str,
    info_id: &'static str,
    description: &'static str,
    fields: &'static [(&'static str, &'static str)],
    kind: VcfProjectionKind,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum VcfProjectionKind {
    AlleleObject,
    AlleleScalar,
    GeneObject,
    ClinvarProtein,
}

const CLINVAR_FIELDS: &[(&str, &str)] = &[
    ("SIGNIFICANCE", "significance"),
    ("REVIEW_STATUS", "reviewStatus"),
    ("PHENOTYPES", "phenotypes"),
    ("VARIANT_CLASS", "variantClass"),
    ("SO_ACCESSION", "soAccession"),
];
const GNOMAD_FIELDS: &[(&str, &str)] = &[
    ("ALL_AF", "allAf"),
    ("ALL_AC", "allAc"),
    ("ALL_AN", "allAn"),
    ("ALL_HC", "allHc"),
    ("AFR_AF", "afrAf"),
    ("AMR_AF", "amrAf"),
    ("ASJ_AF", "asjAf"),
    ("EAS_AF", "easAf"),
    ("FIN_AF", "finAf"),
    ("MID_AF", "midAf"),
    ("NFE_AF", "nfeAf"),
    ("OTH_AF", "othAf"),
    ("REMAINING_AF", "remainingAf"),
    ("SAS_AF", "sasAf"),
];
const DBSNP_FIELDS: &[(&str, &str)] = &[("ID", "id"), ("GLOBAL_MAF", "globalMaf")];
const COSMIC_FIELDS: &[(&str, &str)] = &[("ID", "id"), ("GENE", "gene"), ("COUNT", "count")];
const ONEKG_FIELDS: &[(&str, &str)] = &[
    ("ALL_AF", "allAf"),
    ("AFR_AF", "afrAf"),
    ("AMR_AF", "amrAf"),
    ("EAS_AF", "easAf"),
    ("EUR_AF", "eurAf"),
    ("SAS_AF", "sasAf"),
];
const TOPMED_FIELDS: &[(&str, &str)] = &[
    ("ALL_AF", "allAf"),
    ("ALL_AC", "allAc"),
    ("ALL_AN", "allAn"),
];
const MITOMAP_FIELDS: &[(&str, &str)] = &[("DISEASE", "disease"), ("STATUS", "status")];
const SCORE_FIELDS: &[(&str, &str)] = &[("SCORE", "")];
const SCORE_OBJECT_FIELDS: &[(&str, &str)] = &[("SCORE", "score")];
const DBNSFP_FIELDS: &[(&str, &str)] = &[("SIFT", "sift"), ("POLYPHEN", "polyphen")];
const ALPHAMISSENSE_FIELDS: &[(&str, &str)] =
    &[("PATHOGENICITY", "amPathogenicity"), ("CLASS", "amClass")];
const OMIM_FIELDS: &[(&str, &str)] = &[("MIM_NUMBER", "mimNumber"), ("PHENOTYPES", "phenotypes")];
const GNOMAD_GENE_FIELDS: &[(&str, &str)] = &[
    ("PLI", "pLI"),
    ("LOEUF", "loeuf"),
    ("MIS_Z", "misZ"),
    ("SYN_Z", "synZ"),
];
const CLINVAR_PROTEIN_FIELDS: &[(&str, &str)] = &[("PROTEIN_VARIANTS", "proteinVariants")];

const VCF_PROJECTION_SPECS: &[VcfProjectionSpec] = &[
    VcfProjectionSpec {
        json_key: "clinvar",
        info_id: "FV_CLINVAR",
        description: "fastVEP ClinVar annotations. Format: ALLELE|SIGNIFICANCE|REVIEW_STATUS|PHENOTYPES|VARIANT_CLASS|SO_ACCESSION",
        fields: CLINVAR_FIELDS,
        kind: VcfProjectionKind::AlleleObject,
    },
    VcfProjectionSpec {
        json_key: "gnomad",
        info_id: "FV_GNOMAD",
        description: "fastVEP gnomAD annotations. Format: ALLELE|ALL_AF|ALL_AC|ALL_AN|ALL_HC|AFR_AF|AMR_AF|ASJ_AF|EAS_AF|FIN_AF|MID_AF|NFE_AF|OTH_AF|REMAINING_AF|SAS_AF",
        fields: GNOMAD_FIELDS,
        kind: VcfProjectionKind::AlleleObject,
    },
    VcfProjectionSpec {
        json_key: "dbsnp",
        info_id: "FV_DBSNP",
        description: "fastVEP dbSNP annotations. Format: ALLELE|ID|GLOBAL_MAF",
        fields: DBSNP_FIELDS,
        kind: VcfProjectionKind::AlleleObject,
    },
    VcfProjectionSpec {
        json_key: "cosmic",
        info_id: "FV_COSMIC",
        description: "fastVEP COSMIC annotations. Format: ALLELE|ID|GENE|COUNT",
        fields: COSMIC_FIELDS,
        kind: VcfProjectionKind::AlleleObject,
    },
    VcfProjectionSpec {
        json_key: "oneKg",
        info_id: "FV_1KG",
        description: "fastVEP 1000 Genomes annotations. Format: ALLELE|ALL_AF|AFR_AF|AMR_AF|EAS_AF|EUR_AF|SAS_AF",
        fields: ONEKG_FIELDS,
        kind: VcfProjectionKind::AlleleObject,
    },
    VcfProjectionSpec {
        json_key: "topmed",
        info_id: "FV_TOPMED",
        description: "fastVEP TOPMed annotations. Format: ALLELE|ALL_AF|ALL_AC|ALL_AN",
        fields: TOPMED_FIELDS,
        kind: VcfProjectionKind::AlleleObject,
    },
    VcfProjectionSpec {
        json_key: "mitomap",
        info_id: "FV_MITOMAP",
        description: "fastVEP MitoMap annotations. Format: ALLELE|DISEASE|STATUS",
        fields: MITOMAP_FIELDS,
        kind: VcfProjectionKind::AlleleObject,
    },
    VcfProjectionSpec {
        json_key: "phylop",
        info_id: "FV_PHYLOP",
        description: "fastVEP PhyloP annotations. Format: ALLELE|SCORE",
        fields: SCORE_FIELDS,
        kind: VcfProjectionKind::AlleleScalar,
    },
    VcfProjectionSpec {
        json_key: "gerp",
        info_id: "FV_GERP",
        description: "fastVEP GERP annotations. Format: ALLELE|SCORE",
        fields: SCORE_FIELDS,
        kind: VcfProjectionKind::AlleleScalar,
    },
    VcfProjectionSpec {
        json_key: "dann",
        info_id: "FV_DANN",
        description: "fastVEP DANN annotations. Format: ALLELE|SCORE",
        fields: SCORE_FIELDS,
        kind: VcfProjectionKind::AlleleScalar,
    },
    VcfProjectionSpec {
        json_key: "revel",
        info_id: "FV_REVEL",
        description: "fastVEP REVEL annotations. Format: ALLELE|SCORE",
        fields: SCORE_OBJECT_FIELDS,
        kind: VcfProjectionKind::AlleleObject,
    },
    VcfProjectionSpec {
        json_key: "primateAI",
        info_id: "FV_PRIMATEAI",
        description: "fastVEP PrimateAI annotations. Format: ALLELE|SCORE",
        fields: SCORE_OBJECT_FIELDS,
        kind: VcfProjectionKind::AlleleObject,
    },
    VcfProjectionSpec {
        json_key: "dbnsfp",
        info_id: "FV_DBNSFP",
        description: "fastVEP dbNSFP annotations. Format: ALLELE|SIFT|POLYPHEN",
        fields: DBNSFP_FIELDS,
        kind: VcfProjectionKind::AlleleObject,
    },
    VcfProjectionSpec {
        json_key: "alphaMissense",
        info_id: "FV_ALPHAMISSENSE",
        description: "fastVEP AlphaMissense annotations. Format: ALLELE|PATHOGENICITY|CLASS",
        fields: ALPHAMISSENSE_FIELDS,
        kind: VcfProjectionKind::AlleleObject,
    },
    VcfProjectionSpec {
        json_key: "omim",
        info_id: "FV_OMIM",
        description: "fastVEP OMIM annotations. Format: SYMBOL|MIM_NUMBER|PHENOTYPES",
        fields: OMIM_FIELDS,
        kind: VcfProjectionKind::GeneObject,
    },
    VcfProjectionSpec {
        json_key: "gnomad_genes",
        info_id: "FV_GNOMAD_GENE",
        description: "fastVEP gnomAD gene constraint annotations. Format: SYMBOL|PLI|LOEUF|MIS_Z|SYN_Z",
        fields: GNOMAD_GENE_FIELDS,
        kind: VcfProjectionKind::GeneObject,
    },
    VcfProjectionSpec {
        json_key: "clinvar_protein",
        info_id: "FV_CLINVAR_PROTEIN",
        description: "fastVEP ClinVar protein annotations. Format: SYMBOL|PROTEIN_VARIANTS",
        fields: CLINVAR_PROTEIN_FIELDS,
        kind: VcfProjectionKind::ClinvarProtein,
    },
];

/// A loaded supplementary source with no entry in [`VCF_PROJECTION_SPECS`].
///
/// `VCF_PROJECTION_SPECS` is a closed table of `&'static` descriptors, so it
/// can only name sources that existed when it was written. A `custom_vcf` /
/// `custom_bed` database is named by `--name` at build time and carries
/// whatever fields its source file had, so it matches no entry - and the VCF
/// and tab writers walk that table. The annotation was therefore dropped from
/// both formats with no header, no column and no warning, reaching only
/// `--output-format json`, which iterates the record's own keys (#116).
///
/// The projection is self-describing (`KEY=VALUE&KEY=VALUE`) rather than the
/// fixed pipe list the built-ins use, because the key set genuinely varies:
/// with `--info-fields` unset the builder stores every INFO key each record
/// happens to carry, so there is no layout to declare in the header. Naming
/// each value makes a heterogeneous source readable and keeps the one failure
/// this repository cares about off the table - a reader can never mistake
/// which field a position holds, because position does not carry meaning.
struct CustomProjection {
    /// The source's `json_key`: what `aa.supplementary` / `gene_annotations`
    /// file the record under, and the key the JSON output already used.
    json_key: String,
    /// VCF INFO ID and tab column name, derived from `json_key`.
    info_id: String,
    /// Shared by the `##INFO=` header line and the tab `## COLUMN=` prologue,
    /// so the two formats document the same schema - the same reason the
    /// built-in specs keep one `description` for both.
    description: String,
    /// `true` for an `.oga` gene-level source, whose annotations hang off a
    /// gene symbol rather than an allele.
    gene_keyed: bool,
}

/// Derive a VCF INFO ID for a source with no static spec.
///
/// VCF requires an INFO ID to match `[A-Za-z_][0-9A-Za-z_.]*`, and `--name` is
/// free text, so every character outside that set becomes `_`. The `FV_`
/// prefix marks the field fastVEP-owned - which is what makes re-annotating a
/// VCF replace it rather than append a second copy, via
/// [`vcf_owned_info_ids`] - and guarantees the leading character whatever the
/// name started with.
fn custom_info_id(json_key: &str) -> String {
    let mut id = String::with_capacity(3 + json_key.len());
    id.push_str("FV_");
    for c in json_key.chars() {
        if c.is_ascii_alphanumeric() || c == '_' || c == '.' {
            id.push(c.to_ascii_uppercase());
        } else {
            id.push('_');
        }
    }
    id
}

/// Render a source name for a `Description="..."` attribute.
///
/// The name reaches here from `--name`, so it can hold the quote that would
/// end the attribute or a newline that would end the header line. Neither can
/// be escaped inside a VCF header, so both are replaced.
fn describable_name(json_key: &str) -> String {
    json_key
        .chars()
        .map(|c| match c {
            '"' => '\'',
            '\\' => '/',
            c if c.is_control() => ' ',
            c => c,
        })
        .collect()
}

/// Format fastSA SpliceAI annotations as a VCF-compatible INFO field value.
pub fn format_spliceai_info(vf: &VariationFeature) -> Option<String> {
    format_spliceai_projection(vf)
}

/// Return VCF INFO IDs that fastVEP owns for this run's loaded sources.
///
/// Reads the same [`LoadedSupplementarySpecs`] the header lines and the
/// per-row values come from. Deriving the three lists separately let them
/// disagree, and the disagreement that matters is silent: an emitted `FV_*`
/// field whose `##INFO=` header was never written is an invalid VCF that most
/// readers accept anyway.
pub fn vcf_owned_info_ids(specs: &LoadedSupplementarySpecs) -> Vec<&str> {
    let mut ids = vec!["CSQ"];
    if specs.spliceai {
        ids.push("SpliceAI");
    }
    for spec in specs.loaded_specs() {
        ids.push(spec.info_id);
    }
    for custom in &specs.custom {
        ids.push(&custom.info_id);
    }
    ids
}

/// Generate fastVEP-owned VCF INFO header lines for the loaded sources.
pub fn vcf_info_header_lines(
    specs: &LoadedSupplementarySpecs,
    csq_fields: &[&str],
    sa_only: bool,
) -> Vec<String> {
    let mut headers = if sa_only {
        Vec::new()
    } else {
        vec![csq_header_line(csq_fields)]
    };
    if specs.spliceai {
        headers.push(spliceai_header_line());
    }
    for spec in specs.loaded_specs() {
        headers.push(format!(
            "##INFO=<ID={},Number=.,Type=String,Description=\"{}\">",
            spec.info_id, spec.description
        ));
    }
    for custom in &specs.custom {
        headers.push(format!(
            "##INFO=<ID={},Number=.,Type=String,Description=\"{}\">",
            custom.info_id, custom.description
        ));
    }
    headers
}

/// Parse a VCF INFO header ID from a structured INFO header line.
pub fn vcf_info_header_id(line: &str) -> Option<&str> {
    let rest = line.strip_prefix("##INFO=<ID=")?;
    let end = rest.find([',', '>'])?;
    Some(&rest[..end])
}

/// Precomputed "which supplementary projections are loaded" view, derived
/// once from the loaded SA / gene-annotation keys. Reused by the tab output
/// helpers so each row avoids a fresh O(specs × keys) scan.
pub struct LoadedSupplementarySpecs {
    /// Whether the SpliceAI fastSA source is loaded; emitted as its own
    /// un-namespaced INFO ID / tab column.
    spliceai: bool,
    /// Parallel to `VCF_PROJECTION_SPECS`: `true` when that spec's source is
    /// loaded for this run.
    per_spec: Vec<bool>,
    /// Loaded sources that no `VCF_PROJECTION_SPECS` entry describes - the
    /// user's own `custom_vcf` / `custom_bed` databases. Resolved once here
    /// rather than rediscovered per variant: the set is a property of which
    /// providers opened, and the per-variant loop runs once per
    /// (variant x transcript x allele x source).
    custom: Vec<CustomProjection>,
    /// Total number of extra tab columns (`spliceai` + number of `true`
    /// entries in `per_spec` + `custom`), cached at construction so the
    /// per-row hot path doesn't re-scan for the column count.
    column_count: usize,
}

impl LoadedSupplementarySpecs {
    /// Build the lookup from the SA-provider and gene-provider JSON keys
    /// computed by the annotation pipeline.
    pub fn new(sa_keys: &[String], gene_keys: &[String]) -> Self {
        let spliceai = sa_keys.iter().any(|key| key == "spliceAI");
        let per_spec: Vec<bool> = VCF_PROJECTION_SPECS
            .iter()
            .map(|spec| match spec.kind {
                VcfProjectionKind::GeneObject | VcfProjectionKind::ClinvarProtein => {
                    gene_keys.iter().any(|key| key == spec.json_key)
                }
                VcfProjectionKind::AlleleObject | VcfProjectionKind::AlleleScalar => {
                    sa_keys.iter().any(|key| key == spec.json_key)
                }
            })
            .collect();

        // Anything loaded that the static table does not describe. Seeded with
        // the IDs already in play so a derived ID cannot collide with one of
        // them - `--name gnomad_gene` sanitises to `FV_GNOMAD_GENE`, which the
        // built-in `gnomad_genes` spec owns - and two custom names that
        // sanitise alike cannot collide with each other either. Two INFO
        // fields sharing an ID is a malformed VCF; the suffix keeps the run
        // going and the `Description` names the source it came from.
        let mut taken: std::collections::HashSet<&str> = ["CSQ", "SpliceAI"].into_iter().collect();
        for (spec, loaded) in VCF_PROJECTION_SPECS.iter().zip(per_spec.iter()) {
            if *loaded {
                taken.insert(spec.info_id);
            }
        }
        let mut taken: std::collections::HashSet<String> =
            taken.into_iter().map(str::to_string).collect();
        let mut custom = Vec::new();
        let described = |keys: &[String], gene_keyed: bool| -> Vec<(String, bool)> {
            keys.iter()
                .filter(|key| {
                    *key != "spliceAI"
                        && !VCF_PROJECTION_SPECS
                            .iter()
                            .any(|spec| spec.json_key == key.as_str())
                })
                .map(|key| (key.clone(), gene_keyed))
                .collect()
        };
        for (json_key, gene_keyed) in described(sa_keys, false)
            .into_iter()
            .chain(described(gene_keys, true))
        {
            let base = custom_info_id(&json_key);
            let mut info_id = base.clone();
            let mut n = 2;
            while !taken.insert(info_id.clone()) {
                info_id = format!("{base}_{n}");
                n += 1;
            }
            let subject = if gene_keyed { "SYMBOL" } else { "ALLELE" };
            let description = format!(
                "fastVEP custom annotations from {}. Format: {subject}|KEY=VALUE&KEY=VALUE, \
                 naming the fields the source record carries (key-sorted)",
                describable_name(&json_key)
            );
            custom.push(CustomProjection {
                json_key,
                info_id,
                description,
                gene_keyed,
            });
        }

        let column_count =
            spliceai as usize + per_spec.iter().filter(|b| **b).count() + custom.len();
        Self {
            spliceai,
            per_spec,
            custom,
            column_count,
        }
    }

    /// Number of extra tab columns this lookup will produce (including
    /// SpliceAI). True O(1) — value is cached at construction.
    pub fn column_count(&self) -> usize {
        self.column_count
    }

    /// True when no supplementary source is loaded; callers can skip the
    /// per-row work entirely in that case.
    pub fn is_empty(&self) -> bool {
        self.column_count == 0
    }

    fn loaded_specs(&self) -> impl Iterator<Item = &'static VcfProjectionSpec> + '_ {
        VCF_PROJECTION_SPECS
            .iter()
            .zip(self.per_spec.iter())
            .filter_map(|(spec, loaded)| if *loaded { Some(spec) } else { None })
    }
}

/// Return the supplementary projection column names (info IDs) for tab
/// output, in the same order as `vcf_info_header_lines` emits them. SpliceAI
/// is included first when it's loaded, followed by every loaded
/// `VCF_PROJECTION_SPECS` entry.
pub fn tab_supplementary_column_names(specs: &LoadedSupplementarySpecs) -> Vec<&str> {
    let mut names = Vec::with_capacity(specs.column_count());
    if specs.spliceai {
        names.push("SpliceAI");
    }
    for spec in specs.loaded_specs() {
        names.push(spec.info_id);
    }
    for custom in &specs.custom {
        names.push(&custom.info_id);
    }
    names
}

/// Format `##` header comment lines describing the tab pipe-format for each
/// loaded supplementary source. Mirrors `vcf_info_header_lines` content so
/// the tab and VCF outputs document the same schema.
pub fn tab_supplementary_header_lines(specs: &LoadedSupplementarySpecs) -> Vec<String> {
    let mut lines = Vec::with_capacity(specs.column_count());
    if specs.spliceai {
        lines.push(format!(
            "## COLUMN=<ID=SpliceAI,Description=\"{}\">",
            SPLICEAI_DESCRIPTION
        ));
    }
    for spec in specs.loaded_specs() {
        lines.push(format!(
            "## COLUMN=<ID={},Description=\"{}\">",
            spec.info_id, spec.description
        ));
    }
    for custom in &specs.custom {
        lines.push(format!(
            "## COLUMN=<ID={},Description=\"{}\">",
            custom.info_id, custom.description
        ));
    }
    lines
}

/// Render the per-allele supplementary projection columns for a single tab
/// row, in the same order as `tab_supplementary_column_names`. Missing
/// values render as `-` so the column count is fixed.
pub fn format_supplementary_tab_columns_for_allele(
    vf: &VariationFeature,
    tv: &TranscriptVariation,
    aa: &AlleleAnnotation,
    specs: &LoadedSupplementarySpecs,
) -> Vec<String> {
    let mut cols = Vec::with_capacity(specs.column_count());
    if specs.spliceai {
        cols.push(format_spliceai_for_allele(vf, aa).unwrap_or_else(|| "-".to_string()));
    }
    for spec in specs.loaded_specs() {
        let value = match spec.kind {
            VcfProjectionKind::AlleleObject | VcfProjectionKind::AlleleScalar => {
                format_allele_projection_for_aa(vf, aa, spec)
            }
            VcfProjectionKind::GeneObject => format_gene_projection_for_tv(vf, tv, spec),
            VcfProjectionKind::ClinvarProtein => {
                format_clinvar_protein_projection_for_tv(vf, tv, spec)
            }
        };
        cols.push(value.unwrap_or_else(|| "-".to_string()));
    }
    for custom in &specs.custom {
        let value = if custom.gene_keyed {
            format_custom_gene_projection_for_tv(vf, tv, custom)
        } else {
            format_custom_allele_projection_for_aa(vf, aa, custom)
        };
        cols.push(value.unwrap_or_else(|| "-".to_string()));
    }
    cols
}

/// Format all supplementary VCF INFO projections for an annotated variant.
///
/// `specs` supplies the run's custom sources; the built-in specs are walked in
/// full and yield nothing for a source this variant carries no payload for, as
/// before.
pub fn format_supplementary_vcf_info(
    vf: &VariationFeature,
    specs: &LoadedSupplementarySpecs,
) -> Vec<(String, String)> {
    let mut projected = Vec::new();

    if let Some(value) = format_spliceai_projection(vf) {
        projected.push(("SpliceAI".to_string(), value));
    }

    for spec in VCF_PROJECTION_SPECS {
        let value = match spec.kind {
            VcfProjectionKind::AlleleObject | VcfProjectionKind::AlleleScalar => {
                format_allele_projection(vf, spec)
            }
            VcfProjectionKind::GeneObject => format_gene_projection(vf, spec),
            VcfProjectionKind::ClinvarProtein => format_clinvar_protein_projection(vf, spec),
        };
        if let Some(value) = value {
            projected.push((spec.info_id.to_string(), value));
        }
    }

    for custom in &specs.custom {
        let value = if custom.gene_keyed {
            format_custom_gene_projection(vf, custom)
        } else {
            format_custom_allele_projection(vf, custom)
        };
        if let Some(value) = value {
            projected.push((custom.info_id.clone(), value));
        }
    }

    projected
}

/// Format final VCF INFO by replacing fastVEP-owned fields and appending current projections.
///
/// CSQ is always stripped from `original_info` (it is fastVEP-owned), then
/// re-added only when `csq` is non-empty. This prevents stale `CSQ=` from a
/// re-annotated VCF from leaking through in `--sa-only` mode where the CSQ
/// header has already been removed.
pub fn format_vcf_info_fields(
    original_info: &str,
    vf: &VariationFeature,
    csq: &str,
    specs: &LoadedSupplementarySpecs,
) -> String {
    let mut projections = format_supplementary_vcf_info(vf, specs);
    if !csq.is_empty() {
        projections.push(("CSQ".to_string(), csq.to_string()));
    }

    let mut fields: Vec<String> = if original_info == "." || original_info.is_empty() {
        Vec::new()
    } else {
        original_info
            .split(';')
            .filter(|field| {
                let key = field.split_once('=').map_or(*field, |(key, _)| key);
                if key == "CSQ" {
                    return false;
                }
                !projections.iter().any(|(id, _)| id == key)
            })
            .map(ToOwned::to_owned)
            .collect()
    };

    for (id, value) in projections {
        fields.push(format!("{id}={value}"));
    }

    if fields.is_empty() {
        ".".into()
    } else {
        fields.join(";")
    }
}

fn format_spliceai_projection(vf: &VariationFeature) -> Option<String> {
    // O(n) order-preserving dedup (was O(n^2) via `Vec::contains`).
    let mut values: Vec<String> = Vec::new();
    let mut seen: std::collections::HashSet<String> = std::collections::HashSet::new();
    for tv in &vf.transcript_variations {
        for aa in &tv.allele_annotations {
            let Some(joined) = format_spliceai_for_allele(vf, aa) else {
                continue;
            };
            // The per-aa helper already comma-joins multiple entries for one
            // allele; split them back out so variant-level dedupe stays
            // entry-by-entry rather than across whole allele payloads.
            for value in joined.split(',') {
                // Allocate the owned String once; clone it into the HashSet
                // and only retain it in `values` on first-seen. Earlier code
                // called `value.to_string()` twice per unique entry, doubling
                // allocations in this hot path.
                let owned = value.to_string();
                if seen.insert(owned.clone()) {
                    values.push(owned);
                }
            }
        }
    }
    if values.is_empty() {
        None
    } else {
        Some(values.join(","))
    }
}

fn format_spliceai_entry(allele: &str, json_str: &str) -> Option<String> {
    let value: Value = serde_json::from_str(json_str).ok()?;
    let obj = value.as_object()?;
    let gene = obj.get("gene")?.as_str()?;

    Some(format!(
        "{}|{}|{}|{}|{}|{}|{}|{}|{}|{}",
        allele,
        escape_spliceai_field(gene),
        format_spliceai_float(obj.get("dsAg")?)?,
        format_spliceai_float(obj.get("dsAl")?)?,
        format_spliceai_float(obj.get("dsDg")?)?,
        format_spliceai_float(obj.get("dsDl")?)?,
        format_spliceai_int(obj.get("dpAg")?)?,
        format_spliceai_int(obj.get("dpAl")?)?,
        format_spliceai_int(obj.get("dpDg")?)?,
        format_spliceai_int(obj.get("dpDl")?)?,
    ))
}

fn format_spliceai_float(value: &Value) -> Option<String> {
    value.as_f64().map(|v| format!("{:.2}", v))
}

fn format_spliceai_int(value: &Value) -> Option<String> {
    if let Some(v) = value.as_i64() {
        Some(v.to_string())
    } else {
        value.as_f64().map(|v| format!("{:.0}", v))
    }
}

fn escape_spliceai_field(value: &str) -> String {
    escape_vcf_subfield(value)
}

fn format_allele_projection(vf: &VariationFeature, spec: &VcfProjectionSpec) -> Option<String> {
    // Preserve first-seen order while deduplicating in O(n) instead of the
    // earlier O(n^2) `Vec::contains` scan. A variant overlapping many
    // transcripts with the same SA payload can otherwise quadruple-walk
    // every projection string in a hot loop.
    let mut values: Vec<String> = Vec::new();
    let mut seen: std::collections::HashSet<String> = std::collections::HashSet::new();
    for tv in &vf.transcript_variations {
        for aa in &tv.allele_annotations {
            let Some(joined) = format_allele_projection_for_aa(vf, aa, spec) else {
                continue;
            };
            for value in joined.split(',') {
                // Single allocation per unique value (see `format_spliceai_projection`).
                let owned = value.to_string();
                if seen.insert(owned.clone()) {
                    values.push(owned);
                }
            }
        }
    }
    if values.is_empty() {
        None
    } else {
        Some(values.join(","))
    }
}

/// `GeneIndex::annotate_gene` returns either a single JSON object for a gene
/// with one record or a JSON array of objects when several records hit the
/// same gene symbol. This helper flattens both shapes to a slice of object
/// `Value`s so the projection helpers can produce one pipe entry per record
/// without losing data on multi-record genes.
fn iter_gene_objects(parsed: &Value) -> Vec<&Value> {
    match parsed {
        Value::Object(_) => vec![parsed],
        Value::Array(items) => items.iter().filter(|v| v.is_object()).collect(),
        _ => Vec::new(),
    }
}

fn format_gene_projection(vf: &VariationFeature, spec: &VcfProjectionSpec) -> Option<String> {
    let mut values: Vec<String> = Vec::new();
    let mut seen: std::collections::HashSet<String> = std::collections::HashSet::new();
    for ga in &vf.gene_annotations {
        if ga.json_key != spec.json_key {
            continue;
        }
        let Ok(parsed) = serde_json::from_str::<Value>(&ga.json_string) else {
            continue;
        };
        for object in iter_gene_objects(&parsed) {
            let mut parts = Vec::with_capacity(spec.fields.len() + 1);
            parts.push(escape_vcf_subfield(&ga.gene_symbol));
            if let Some(obj) = object.as_object() {
                for (_, json_key) in spec.fields {
                    parts.push(
                        obj.get(*json_key)
                            .map(json_value_to_vcf)
                            .unwrap_or_default(),
                    );
                }
            }
            let value = parts.join("|");
            if seen.insert(value.clone()) {
                values.push(value);
            }
        }
    }
    if values.is_empty() {
        None
    } else {
        Some(values.join(","))
    }
}

fn format_spliceai_for_allele(vf: &VariationFeature, aa: &AlleleAnnotation) -> Option<String> {
    // Establish that this allele carries a payload for the source before
    // building anything. Every annotation is offered to every loaded source,
    // so most calls match nothing and the whole body below was previously
    // paid for to reach an empty `values` and return `None`.
    if !aa.supplementary.iter().any(|(key, _)| key == "spliceAI") {
        return None;
    }
    let allele = uploaded_allele_for_annotation(vf, &aa.allele);
    let mut values: Vec<String> = Vec::new();
    let mut seen: std::collections::HashSet<String> = std::collections::HashSet::new();
    for (key, json_str) in &aa.supplementary {
        if key != "spliceAI" {
            continue;
        }
        if let Some(value) = format_spliceai_entry(&allele, json_str) {
            if seen.insert(value.clone()) {
                values.push(value);
            }
        }
    }
    if values.is_empty() {
        None
    } else {
        Some(values.join(","))
    }
}

fn format_allele_projection_for_aa(
    vf: &VariationFeature,
    aa: &AlleleAnnotation,
    spec: &VcfProjectionSpec,
) -> Option<String> {
    // As in `format_spliceai_for_allele`: no payload for this source means no
    // column, and finding that out costs a key comparison rather than an
    // allele render, a `HashSet`, and two `Vec`s.
    if !aa.supplementary.iter().any(|(key, _)| key == spec.json_key) {
        return None;
    }
    let allele = uploaded_allele_for_annotation(vf, &aa.allele);
    let mut values: Vec<String> = Vec::new();
    let mut seen: std::collections::HashSet<String> = std::collections::HashSet::new();
    for (key, json_str) in &aa.supplementary {
        if key != spec.json_key {
            continue;
        }
        // Truncated or otherwise non-JSON supplementary payloads used to fall
        // back to `Value::String(json_str.clone())` here, which silently
        // pushed the raw blob into the projection. Now we log at debug level
        // and skip this entry — callers see a missing column instead of a
        // mis-shaped one.
        let parsed = match serde_json::from_str::<Value>(json_str) {
            Ok(v) => v,
            Err(e) => {
                log::debug!(
                    "Skipping non-JSON supplementary payload for key='{}': {} (payload snippet: {})",
                    spec.json_key,
                    e,
                    json_str.chars().take(80).collect::<String>(),
                );
                continue;
            }
        };
        let entries = match spec.kind {
            VcfProjectionKind::AlleleScalar => {
                vec![format!(
                    "{}|{}",
                    escape_vcf_subfield(&allele),
                    json_value_to_vcf(&parsed)
                )]
            }
            _ => format_object_projection_entries(&allele, &parsed, spec.fields),
        };
        for value in entries {
            if seen.insert(value.clone()) {
                values.push(value);
            }
        }
    }
    if values.is_empty() {
        None
    } else {
        Some(values.join(","))
    }
}

/// Append one `KEY=VALUE` pair, `&`-separated from what is already there.
///
/// Skips a value with no scalar rendering. `escape_vcf_subfield` percent-encodes
/// `=`, `&`, `|` and `,`, so neither a key nor a value can contain a separator
/// and the pair list splits back apart exactly.
fn push_custom_pair(out: &mut String, key: &str, value: &Value) {
    let rendered = match value {
        Value::Null | Value::Array(_) | Value::Object(_) => return,
        leaf => escape_vcf_subfield(&json_leaf_to_string(leaf)),
    };
    if !out.is_empty() {
        out.push('&');
    }
    out.push_str(key);
    out.push('=');
    out.push_str(&rendered);
}

/// Render a custom source's record object as `KEY=VALUE&KEY=VALUE`.
///
/// `serde_json::Map` is a `BTreeMap` here (the `preserve_order` feature is
/// off), so the pairs come out key-sorted and every record of a source lists
/// its fields in the same order - which is what makes the column greppable
/// even when the records are heterogeneous.
///
/// An array value repeats its key (`AF=0.1&AF=0.2`) rather than joining the
/// elements the way [`json_value_to_vcf`] does, because `&` already separates
/// pairs: joining would produce `AF=0.1&0.2`, whose second element has no key
/// and reads as a malformed pair. Nested objects are dropped - neither custom
/// builder emits one, and flattening would have to invent a separator for the
/// compound key.
fn format_custom_pairs(object: &serde_json::Map<String, Value>) -> String {
    let mut out = String::new();
    for (key, value) in object {
        let key = escape_vcf_subfield(key);
        match value {
            Value::Array(items) => {
                for item in items {
                    push_custom_pair(&mut out, &key, item);
                }
            }
            leaf => push_custom_pair(&mut out, &key, leaf),
        }
    }
    out
}

/// One projection entry per record object: `SUBJECT|KEY=VALUE&KEY=VALUE`.
///
/// `subject` is the allele or the gene symbol, matching the built-in specs'
/// leading field. The pair list is emitted even when empty, so the entry is
/// always two `|`-separated fields: a three-column BED carries no attributes
/// at all, and `region|` still says the allele hit it - the same way the
/// built-in specs render a value they have no data for.
fn format_custom_entries(subject: &str, parsed: &Value) -> Vec<String> {
    let records: Vec<&Value> = match parsed {
        Value::Array(items) => items.iter().collect(),
        single => vec![single],
    };
    records
        .into_iter()
        .filter_map(|record| {
            let subject = escape_vcf_subfield(subject);
            match record {
                Value::Object(object) => Some(format!("{subject}|{}", format_custom_pairs(object))),
                // A source whose record is a bare scalar rather than an object
                // has nothing to name, so it projects like the built-in
                // `AlleleScalar` specs do.
                Value::Null | Value::Array(_) => None,
                leaf => Some(format!(
                    "{subject}|{}",
                    escape_vcf_subfield(&json_leaf_to_string(leaf))
                )),
            }
        })
        .collect()
}

/// Project a custom allele-level source for one allele annotation (tab rows).
fn format_custom_allele_projection_for_aa(
    vf: &VariationFeature,
    aa: &AlleleAnnotation,
    custom: &CustomProjection,
) -> Option<String> {
    // Establish a payload before building anything, as the built-in helpers
    // do: every loaded source is offered every annotation, so most calls match
    // nothing and the allele render, the dedupe set and the value vector would
    // all be paid for to reach an empty result.
    if !aa
        .supplementary
        .iter()
        .any(|(key, _)| *key == custom.json_key)
    {
        return None;
    }
    let allele = uploaded_allele_for_annotation(vf, &aa.allele);
    let mut values: Vec<String> = Vec::new();
    let mut seen: std::collections::HashSet<String> = std::collections::HashSet::new();
    for (key, json_str) in &aa.supplementary {
        if *key != custom.json_key {
            continue;
        }
        // Same policy as `format_allele_projection_for_aa`: a truncated or
        // otherwise non-JSON payload is skipped with a debug line rather than
        // pushed through raw, so a caller sees a missing column instead of a
        // mis-shaped one.
        let parsed = match serde_json::from_str::<Value>(json_str) {
            Ok(v) => v,
            Err(e) => {
                log::debug!(
                    "Skipping non-JSON supplementary payload for custom key='{}': {} (payload snippet: {})",
                    custom.json_key,
                    e,
                    json_str.chars().take(80).collect::<String>(),
                );
                continue;
            }
        };
        for value in format_custom_entries(&allele, &parsed) {
            if seen.insert(value.clone()) {
                values.push(value);
            }
        }
    }
    if values.is_empty() {
        None
    } else {
        Some(values.join(","))
    }
}

/// Project a custom allele-level source across every allele (VCF INFO).
fn format_custom_allele_projection(
    vf: &VariationFeature,
    custom: &CustomProjection,
) -> Option<String> {
    let mut values: Vec<String> = Vec::new();
    let mut seen: std::collections::HashSet<String> = std::collections::HashSet::new();
    for tv in &vf.transcript_variations {
        for aa in &tv.allele_annotations {
            let Some(joined) = format_custom_allele_projection_for_aa(vf, aa, custom) else {
                continue;
            };
            for value in joined.split(',') {
                let owned = value.to_string();
                if seen.insert(owned.clone()) {
                    values.push(owned);
                }
            }
        }
    }
    if values.is_empty() {
        None
    } else {
        Some(values.join(","))
    }
}

/// Project a custom gene-level source for one transcript (tab rows).
fn format_custom_gene_projection_for_tv(
    vf: &VariationFeature,
    tv: &TranscriptVariation,
    custom: &CustomProjection,
) -> Option<String> {
    format_custom_gene_entries(vf, custom, tv.gene_symbol.as_deref())
}

/// Project a custom gene-level source across every gene (VCF INFO).
fn format_custom_gene_projection(
    vf: &VariationFeature,
    custom: &CustomProjection,
) -> Option<String> {
    format_custom_gene_entries(vf, custom, None)
}

/// Shared body for the two gene-level custom projections.
///
/// `symbol_filter` restricts the annotations to one gene for a tab row that
/// has a symbol; `None` emits every match, which is both the VCF behaviour and
/// what a tab row for a transcript with no symbol needs so it does not lose
/// data. Same split as `format_gene_projection` / `format_gene_projection_for_tv`.
fn format_custom_gene_entries(
    vf: &VariationFeature,
    custom: &CustomProjection,
    symbol_filter: Option<&str>,
) -> Option<String> {
    if !vf
        .gene_annotations
        .iter()
        .any(|ga| ga.json_key == custom.json_key)
    {
        return None;
    }
    let mut values: Vec<String> = Vec::new();
    let mut seen: std::collections::HashSet<String> = std::collections::HashSet::new();
    for ga in &vf.gene_annotations {
        if ga.json_key != custom.json_key {
            continue;
        }
        if let Some(want) = symbol_filter {
            if ga.gene_symbol != want {
                continue;
            }
        }
        let Ok(parsed) = serde_json::from_str::<Value>(&ga.json_string) else {
            continue;
        };
        for value in format_custom_entries(&ga.gene_symbol, &parsed) {
            if seen.insert(value.clone()) {
                values.push(value);
            }
        }
    }
    if values.is_empty() {
        None
    } else {
        Some(values.join(","))
    }
}

fn format_gene_projection_for_tv(
    vf: &VariationFeature,
    tv: &TranscriptVariation,
    spec: &VcfProjectionSpec,
) -> Option<String> {
    // Gene annotations are keyed by HGNC symbol. When the transcript has a
    // known symbol, filter to that gene; otherwise emit every matching
    // annotation for this variant — same dedupe semantics the VCF emitter
    // uses — so tab rows for transcripts with no symbol don't lose data.
    let symbol_filter = tv.gene_symbol.as_deref();
    let mut values: Vec<String> = Vec::new();
    let mut seen: std::collections::HashSet<String> = std::collections::HashSet::new();
    for ga in &vf.gene_annotations {
        if ga.json_key != spec.json_key {
            continue;
        }
        if let Some(want) = symbol_filter {
            if ga.gene_symbol != want {
                continue;
            }
        }
        let Ok(parsed) = serde_json::from_str::<Value>(&ga.json_string) else {
            continue;
        };
        for object in iter_gene_objects(&parsed) {
            let mut parts = Vec::with_capacity(spec.fields.len() + 1);
            parts.push(escape_vcf_subfield(&ga.gene_symbol));
            if let Some(obj) = object.as_object() {
                for (_, json_key) in spec.fields {
                    parts.push(
                        obj.get(*json_key)
                            .map(json_value_to_vcf)
                            .unwrap_or_default(),
                    );
                }
            }
            let value = parts.join("|");
            if seen.insert(value.clone()) {
                values.push(value);
            }
        }
    }
    if values.is_empty() {
        None
    } else {
        Some(values.join(","))
    }
}

fn format_clinvar_protein_projection_for_tv(
    vf: &VariationFeature,
    tv: &TranscriptVariation,
    spec: &VcfProjectionSpec,
) -> Option<String> {
    // See `format_gene_projection_for_tv`: fall back to variant-level
    // dedupe when the transcript has no associated gene symbol.
    let symbol_filter = tv.gene_symbol.as_deref();
    let mut out = String::new();
    // `(start, end)` of each accepted value inside `out`, for the dedupe.
    let mut seen: Vec<(usize, usize)> = Vec::new();
    for ga in &vf.gene_annotations {
        if ga.json_key != spec.json_key {
            continue;
        }
        if let Some(want) = symbol_filter {
            if ga.gene_symbol != want {
                continue;
            }
        }
        let Ok(parsed) = serde_json::from_str::<Value>(&ga.json_string) else {
            continue;
        };
        // Written straight into the output, with the `,` that joins two
        // genes' entries, so neither the value nor the dedupe key is ever a
        // second copy of the rendered list.
        let start = out.len();
        if !out.is_empty() {
            out.push(',');
        }
        let value_at = out.len();
        if !push_clinvar_protein_value(&mut out, &ga.gene_symbol, &parsed) {
            out.truncate(start);
            continue;
        }
        // One entry per overlapping gene, so this list is one to a handful
        // long: a scan of it beats hashing a value that can be kilobytes, and
        // it compares in place rather than cloning to build a key.
        if seen
            .iter()
            .any(|&(from, to)| out[from..to] == out[value_at..])
        {
            out.truncate(start);
            continue;
        }
        seen.push((value_at, out.len()));
    }
    if out.is_empty() {
        None
    } else {
        Some(out)
    }
}

/// Pull every `proteinVariants[*]` entry out of a ClinVar-protein JSON
/// payload and render the `&`-joined `pos:ref>alt:sig` list. Accepts either
/// a single object payload or the JSON array `GeneIndex::annotate_gene`
/// emits when a gene has multiple records.
/// Append every `proteinVariants` record in `parsed` to `out`, `&`-separated.
///
/// Returns whether anything was written, which is also what tells the caller
/// whether the `SYMBOL|` prefix it wrote first has a value to belong to.
fn push_clinvar_protein_records(out: &mut String, parsed: &Value) -> bool {
    let mut wrote = false;
    for obj in iter_gene_objects(parsed) {
        let Some(vars) = obj.get("proteinVariants").and_then(|v| v.as_array()) else {
            continue;
        };
        for v in vars {
            // The separator belongs to the record that follows it, so a
            // skipped record cannot leave a trailing `&`.
            wrote |= push_clinvar_protein_record(out, v, wrote);
        }
    }
    wrote
}

/// Append `SYMBOL|pos:ref>alt:sig&...` for one gene annotation, or leave `out`
/// exactly as it was if the payload holds no renderable record.
///
/// The value is built in place rather than returned, because the tab writer
/// reaches this once per (variant x transcript x allele): composing it with a
/// `format!` and then cloning it into a dedupe set copied the whole rendered
/// list twice on each of those, and the list grows with how many records
/// ClinVar has for the gene.
fn push_clinvar_protein_value(out: &mut String, symbol: &str, parsed: &Value) -> bool {
    let start = out.len();
    escape_vcf_subfield_into(symbol, out);
    out.push('|');
    if push_clinvar_protein_records(out, parsed) {
        return true;
    }
    out.truncate(start);
    false
}

/// Append one `&`-separated `pos:ref>alt:sig` record to `out`, or nothing if
/// the entry does not name a substitution.
///
/// Writes into the caller's buffer rather than returning a `String`, the way
/// `format_csq_entry_into` does and for the same reason: with a gene-level
/// source loaded this runs once per (variant x transcript) for the tab writer,
/// and ClinVar carries tens of protein variants for a well-studied gene - FHL1
/// has 31 in the report that prompted #123. Composing each record through four
/// escaped `String`s and a `format!` built nine allocations per record; this
/// builds none beyond the growth of `out` itself.
fn push_clinvar_protein_record(out: &mut String, v: &Value, needs_separator: bool) -> bool {
    // An entry missing any of the three coordinates is not a substitution this
    // field can describe. A *present* but null one renders empty, which is what
    // the previous `json_leaf_to_string` spelling did, so an odd record still
    // keeps its place in the list rather than shifting the ones after it.
    let (Some(pos), Some(ref_aa), Some(alt_aa)) = (v.get("pos"), v.get("refAa"), v.get("altAa"))
    else {
        return false;
    };

    if needs_separator {
        out.push('&');
    }
    push_clinvar_protein_leaf(out, pos);
    out.push(':');
    push_clinvar_protein_leaf(out, ref_aa);
    out.push('>');
    push_clinvar_protein_leaf(out, alt_aa);
    out.push(':');
    if let Some(sig) = v.get("sig") {
        push_clinvar_protein_leaf(out, sig);
    }
    true
}

/// Append one leaf of a `pos:ref>alt:sig` record, escaped.
///
/// `:` and `>` are this field's own delimiters and no other `FV_*` field uses
/// them, so [`escape_vcf_subfield`] leaves both alone and this encodes them
/// for the values that sit *between* them. Order matters: the shared escaper
/// runs on the same character, so a literal `%` becomes `%25` and the result
/// stays exactly one decode pass deep.
///
/// No payload can reach here containing either character today - `pos` is a
/// number, `refAa` / `altAa` are single residues, and `sig` is one of the four
/// `&'static str` constants `normalise_significance` returns. This is what
/// keeps that true if a fifth term or a three-letter residue code ever
/// arrives, rather than silently splitting a record in two.
fn push_clinvar_protein_leaf(out: &mut String, value: &Value) {
    match value {
        // A number renders as digits, `-`, `.` and `e`, none of which either
        // escaper touches, so it is written straight through.
        Value::Number(n) => {
            let _ = write!(out, "{}", n);
        }
        Value::Bool(b) => out.push_str(if *b { "true" } else { "false" }),
        Value::String(text) => {
            for c in text.chars() {
                match c {
                    ':' => out.push_str("%3A"),
                    '>' => out.push_str("%3E"),
                    _ => push_escaped_subfield_char(c, out),
                }
            }
        }
        // `json_leaf_to_string` renders these empty too.
        Value::Null | Value::Array(_) | Value::Object(_) => {}
    }
}

fn format_clinvar_protein_projection(
    vf: &VariationFeature,
    spec: &VcfProjectionSpec,
) -> Option<String> {
    let mut out = String::new();
    // `(start, end)` of each accepted value inside `out`, for the dedupe.
    let mut seen: Vec<(usize, usize)> = Vec::new();
    for ga in &vf.gene_annotations {
        if ga.json_key != spec.json_key {
            continue;
        }
        let Ok(parsed) = serde_json::from_str::<Value>(&ga.json_string) else {
            continue;
        };
        // Written straight into the output, with the `,` that joins two
        // genes' entries, so neither the value nor the dedupe key is ever a
        // second copy of the rendered list.
        let start = out.len();
        if !out.is_empty() {
            out.push(',');
        }
        let value_at = out.len();
        if !push_clinvar_protein_value(&mut out, &ga.gene_symbol, &parsed) {
            out.truncate(start);
            continue;
        }
        // One entry per overlapping gene, so this list is one to a handful
        // long: a scan of it beats hashing a value that can be kilobytes, and
        // it compares in place rather than cloning to build a key.
        if seen
            .iter()
            .any(|&(from, to)| out[from..to] == out[value_at..])
        {
            out.truncate(start);
            continue;
        }
        seen.push((value_at, out.len()));
    }
    if out.is_empty() {
        None
    } else {
        Some(out)
    }
}

fn format_object_projection_entries(
    allele: &str,
    value: &Value,
    fields: &[(&str, &str)],
) -> Vec<String> {
    let values: Vec<&Value> = match value {
        Value::Array(items) => items.iter().collect(),
        _ => vec![value],
    };
    values
        .into_iter()
        .filter_map(|value| {
            let obj = value.as_object()?;
            let mut parts = Vec::with_capacity(fields.len() + 1);
            parts.push(escape_vcf_subfield(allele));
            for (_, json_key) in fields {
                if json_key.is_empty() {
                    parts.push(json_value_to_vcf(value));
                } else {
                    parts.push(
                        obj.get(*json_key)
                            .map(json_value_to_vcf)
                            .unwrap_or_default(),
                    );
                }
            }
            Some(parts.join("|"))
        })
        .collect()
}

fn json_value_to_vcf(value: &Value) -> String {
    match value {
        Value::Null => String::new(),
        Value::Array(items) => items
            .iter()
            .map(json_value_to_vcf)
            .filter(|s| !s.is_empty())
            .collect::<Vec<_>>()
            .join("&"),
        Value::Object(_) => String::new(),
        _ => escape_vcf_subfield(&json_leaf_to_string(value)),
    }
}

fn json_leaf_to_string(value: &Value) -> String {
    match value {
        Value::Null => String::new(),
        Value::Bool(v) => v.to_string(),
        Value::Number(v) => v.to_string(),
        Value::String(v) => v.clone(),
        Value::Array(_) | Value::Object(_) => String::new(),
    }
}

/// Percent-encode the characters that would otherwise be read as structure
/// inside an `FV_*` pipe field.
///
/// The set is "illegal, one of our own delimiters, or a character that breaks
/// a reader of the column". VCF 4.3 restricts only three inside an INFO value
/// (`;` and `=` are not permitted, `,` only as a list delimiter); `|` and `&`
/// are the delimiters these fields are built from; `%` is the escape
/// introducer, so it has to be encoded for a single decode pass to be exact;
/// CR / LF / TAB would end the value or the record.
///
/// **Exactly one character was dropped from this set**, the colon, because it
/// is reserved by nothing and encoding it only cost legibility - `#123`
/// reported an `FV_CLINVAR_PROTEIN` value with sixty `%3A` in it. Ensembl VEP
/// 115.1 writes `:` through literally (`validation/run_escaping_probe.sh`
/// measures this), and so does fastVEP's own CSQ writer - every HGVSc carries
/// one - so encoding it here contradicted the same record's CSQ column.
///
/// Two characters the spec does *not* require are still encoded, and both for
/// the same reason: matching VEP is worth nothing on a field VEP does not
/// emit, while surviving the tools that read these columns is worth a lot.
///
/// - **Space** (`%20`): a whitespace-free INFO column survives an `awk` or
///   `cut -d' '` pipeline, and `%20` is reversible where VEP's `_`
///   substitution is not - `Breast_cancer` cannot be told back from
///   `Breast cancer`. The single-decode-pass promise in
///   docs/SUPPLEMENTARY_ANNOTATIONS.md depends on that.
/// - **Double quote** (`%22`): VEP writes it literally, and dropping it here
///   to match cost a real invariant - `R`'s `read.delim` and Python's `csv`
///   both default to `"` as a quote character, so one inch mark in a disease
///   name or a `--custom` BED field silently swallows the rest of the record.
///   The design spec in docs/superpowers/ states the property, and the
///   `"VCF INFO must not contain JSON quotes"` assertion below is what holds
///   it; neither survives a reader that only has VEP's behaviour to go on.
fn escape_vcf_subfield(value: &str) -> String {
    let mut escaped = String::with_capacity(value.len());
    escape_vcf_subfield_into(value, &mut escaped);
    escaped
}

/// [`escape_vcf_subfield`] appending to a caller-owned buffer.
fn escape_vcf_subfield_into(value: &str, out: &mut String) {
    for c in value.chars() {
        push_escaped_subfield_char(c, out);
    }
}

/// One character of a pipe field, percent-encoded if it would otherwise be
/// read as structure. The single definition of the set described on
/// [`escape_vcf_subfield`]; `push_clinvar_protein_leaf` shares it so a field
/// with delimiters of its own cannot drift from the shared rules.
#[inline]
fn push_escaped_subfield_char(c: char, out: &mut String) {
    match c {
        ';' => out.push_str("%3B"),
        '=' => out.push_str("%3D"),
        '%' => out.push_str("%25"),
        ',' => out.push_str("%2C"),
        '\r' => out.push_str("%0D"),
        '\n' => out.push_str("%0A"),
        '\t' => out.push_str("%09"),
        ' ' => out.push_str("%20"),
        '"' => out.push_str("%22"),
        '|' => out.push_str("%7C"),
        '&' => out.push_str("%26"),
        _ => out.push(c),
    }
}

fn uploaded_allele_for_annotation(vf: &VariationFeature, allele: &Allele) -> String {
    // Rendering the allele is the fallback, not the starting point: it used to
    // be built on entry and then thrown away on the common path, once per
    // (variant x transcript x allele x source). Splitting the ALT column
    // lazily rather than collecting it drops the other allocation.
    let Some(vcf) = &vf.vcf_fields else {
        return allele.to_string();
    };

    vf.alt_alleles
        .iter()
        .position(|alt| alt == allele)
        .and_then(|idx| vcf.alt.split(',').nth(idx))
        .map_or_else(|| allele.to_string(), str::to_string)
}

/// The optional columns of VEP's tab output, by the VEP flag that asks for each.
///
/// VEP adds a tab column only when its flag is given, and orders the added
/// columns by the `FLAG_FIELDS` table of `Constants.pm` (ensembl-vep
/// release/115), not by the order of the flags on the command line.
/// [`ColumnFlags::columns`] does the same. `--hgvs` is not here: the run
/// already carries it for the HGVS values themselves.
///
/// A run without any of these has the 17 columns the tab output always had.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct ColumnFlags {
    /// `--show-ref-allele`: `REF_ALLELE`.
    pub show_ref_allele: bool,
    /// `--variant-class`: `VARIANT_CLASS`.
    pub variant_class: bool,
    /// `--symbol`: `SYMBOL`, `SYMBOL_SOURCE`, `HGNC_ID`.
    pub symbol: bool,
    /// `--biotype`: `BIOTYPE`.
    pub biotype: bool,
    /// `--canonical`: `CANONICAL`.
    pub canonical: bool,
    /// `--numbers`: `EXON`, `INTRON`.
    pub numbers: bool,
    /// `--hgvsg`: `HGVSg`.
    pub hgvsg: bool,
}

impl ColumnFlags {
    /// What `--everything` turns on, among the columns fastVEP writes: VEP's set
    /// (`OPTION_SETS` in `Config.pm`) has `--symbol`, `--numbers`, `--canonical`,
    /// `--biotype` and `--variant_class` (and `--hgvs`, which the caller adds)
    /// but neither `--hgvsg` nor `--show_ref_allele`.
    pub fn everything() -> Self {
        Self {
            variant_class: true,
            symbol: true,
            biotype: true,
            canonical: true,
            numbers: true,
            ..Self::default()
        }
    }

    /// Every flag set in either.
    pub fn union(self, other: Self) -> Self {
        Self {
            show_ref_allele: self.show_ref_allele || other.show_ref_allele,
            variant_class: self.variant_class || other.variant_class,
            symbol: self.symbol || other.symbol,
            biotype: self.biotype || other.biotype,
            canonical: self.canonical || other.canonical,
            numbers: self.numbers || other.numbers,
            hgvsg: self.hgvsg || other.hgvsg,
        }
    }

    /// The columns to write, in VEP's order. `hgvs` says whether the run has
    /// `--hgvs`, which adds `HGVSc`, `HGVSp` and `HGVS_OFFSET`.
    pub fn columns(self, hgvs: bool) -> Vec<OptionalColumn> {
        let mut columns = Vec::new();
        if self.show_ref_allele {
            columns.push(OptionalColumn::RefAllele);
        }
        if self.variant_class {
            columns.push(OptionalColumn::VariantClass);
        }
        if self.symbol {
            columns.extend([
                OptionalColumn::Symbol,
                OptionalColumn::SymbolSource,
                OptionalColumn::HgncId,
            ]);
        }
        if self.biotype {
            columns.push(OptionalColumn::Biotype);
        }
        if self.canonical {
            columns.push(OptionalColumn::Canonical);
        }
        if self.numbers {
            columns.extend([OptionalColumn::Exon, OptionalColumn::Intron]);
        }
        if hgvs {
            columns.extend([
                OptionalColumn::Hgvsc,
                OptionalColumn::Hgvsp,
                OptionalColumn::HgvsOffset,
            ]);
        }
        if self.hgvsg {
            columns.push(OptionalColumn::Hgvsg);
        }
        columns
    }

    /// Whether any column is asked for.
    pub fn any(self) -> bool {
        self != Self::default()
    }
}

/// One of VEP's optional columns. The tab writer and the CSQ writer give each the
/// same value; `tab_columns_vep_flags.rs` holds one against the other.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OptionalColumn {
    RefAllele,
    VariantClass,
    Symbol,
    SymbolSource,
    HgncId,
    Biotype,
    Canonical,
    Exon,
    Intron,
    Hgvsc,
    Hgvsp,
    HgvsOffset,
    Hgvsg,
}

impl OptionalColumn {
    /// The header, which is VEP's name for the column.
    pub fn name(self) -> &'static str {
        match self {
            Self::RefAllele => "REF_ALLELE",
            Self::VariantClass => "VARIANT_CLASS",
            Self::Symbol => "SYMBOL",
            Self::SymbolSource => "SYMBOL_SOURCE",
            Self::HgncId => "HGNC_ID",
            Self::Biotype => "BIOTYPE",
            Self::Canonical => "CANONICAL",
            Self::Exon => "EXON",
            Self::Intron => "INTRON",
            Self::Hgvsc => "HGVSc",
            Self::Hgvsp => "HGVSp",
            Self::HgvsOffset => "HGVS_OFFSET",
            Self::Hgvsg => "HGVSg",
        }
    }

    /// The cell of a transcript row. `-` for a value the row does not have, and
    /// the same text the CSQ writer gives the field.
    fn cell(
        self,
        vf: &VariationFeature,
        tv: &TranscriptVariation,
        aa: &AlleleAnnotation,
    ) -> String {
        fn or_dash(v: Option<&str>) -> String {
            v.filter(|s| !s.is_empty()).unwrap_or("-").to_string()
        }
        fn ordinal(v: Option<Ordinals>) -> String {
            v.map_or_else(
                || "-".to_string(),
                |o| {
                    let mut cell = String::new();
                    o.write_to(&mut cell);
                    cell
                },
            )
        }
        match self {
            Self::RefAllele | Self::VariantClass | Self::Hgvsg => self.variant_cell(vf, Some(aa)),
            Self::Symbol => or_dash(tv.gene_symbol.as_deref()),
            Self::SymbolSource => or_dash(tv.symbol_source.as_deref()),
            Self::HgncId => or_dash(tv.hgnc_id.as_deref()),
            Self::Biotype => or_dash(Some(&tv.biotype)),
            Self::Canonical => if tv.canonical { "YES" } else { "-" }.to_string(),
            Self::Exon => ordinal(aa.location.exon),
            Self::Intron => ordinal(aa.location.intron),
            Self::Hgvsc => or_dash(aa.hgvsc.as_deref()),
            // VEP encodes `=` as `%3D` in the tab file as well as in the VCF
            // (`OutputFactory`, "URI encode =", unless --no_escape)
            Self::Hgvsp => or_dash(
                aa.hgvsp
                    .as_deref()
                    .map(|h| h.replace('=', "%3D"))
                    .as_deref(),
            ),
            Self::HgvsOffset => aa
                .hgvs_offset
                .map_or_else(|| "-".to_string(), |o| o.to_string()),
        }
    }

    /// The cell of a row that has no allele annotation at all (the row written for
    /// a variant with no transcript row): the variant's own values, `-` for the
    /// rest. A row that does have one passes it, for `HGVSg`.
    fn variant_cell(self, vf: &VariationFeature, aa: Option<&AlleleAnnotation>) -> String {
        match self {
            Self::RefAllele => vf.ref_allele.to_string(),
            // a structural variant's class comes from the parser in VEP
            Self::VariantClass if vf.variant_type.is_structural() => "-".to_string(),
            Self::VariantClass => so_variation_class(&vf.allele_string).to_string(),
            Self::Hgvsg => aa
                .and_then(|aa| aa.hgvsg.as_deref())
                .filter(|s| !s.is_empty())
                .unwrap_or("-")
                .to_string(),
            _ => "-".to_string(),
        }
    }
}

/// The first 17 columns of the tab output, whatever else a run adds: VEP's 13
/// base columns (`DEFAULT_OUTPUT_COLS` in `Constants.pm`) and the four every
/// fastVEP tab file has always had after them. The header is written from this
/// list, and the row writer's 17 cells are held to it by a test.
pub const TAB_BASE_COLUMNS: [&str; 17] = [
    "Uploaded_variation",
    "Location",
    "Allele",
    "Gene",
    "Feature",
    "Feature_type",
    "Consequence",
    "cDNA_position",
    "CDS_position",
    "Protein_position",
    "Amino_acids",
    "Codons",
    "Existing_variation",
    "IMPACT",
    "DISTANCE",
    "STRAND",
    "FLAGS",
];

/// How many of [`TAB_BASE_COLUMNS`] come before `IMPACT`: VEP's
/// `DEFAULT_OUTPUT_COLS`. `FLAG_FIELDS` orders `REF_ALLELE` between them and
/// `IMPACT`, and `PICK` straight after `FLAGS`.
const COLUMNS_BEFORE_IMPACT: usize = 13;

/// The header line of the tab file, `#` included, for a run with these columns.
///
/// The row writer ([`format_tab_line_with`]) writes its cells in this order, and
/// a test holds the two to each other for every combination of options: the order
/// is VEP's `FLAG_FIELDS` (`Constants.pm`, release/116), with the columns fastVEP
/// adds on its own (`REF` after `Allele`, the supplementary sources, `QC_CLASS`)
/// where it has always put them.
///
/// `--sa-only` writes no transcript columns, so it has none of VEP's.
pub fn tab_header(
    sa_only: bool,
    explicit_ref: bool,
    optional: &[OptionalColumn],
    pick_column: bool,
    supplementary: &[&str],
    qc_class: bool,
) -> String {
    let mut names: Vec<&str> = Vec::new();
    if sa_only {
        names.extend(&TAB_BASE_COLUMNS[..3]);
        if explicit_ref {
            names.push("REF");
        }
    } else {
        names.extend(&TAB_BASE_COLUMNS[..3]);
        if explicit_ref {
            names.push("REF");
        }
        names.extend(&TAB_BASE_COLUMNS[3..COLUMNS_BEFORE_IMPACT]);
        let mut optional = optional;
        if let [OptionalColumn::RefAllele, rest @ ..] = optional {
            names.push(OptionalColumn::RefAllele.name());
            optional = rest;
        }
        names.extend(&TAB_BASE_COLUMNS[COLUMNS_BEFORE_IMPACT..]);
        if pick_column {
            names.push("PICK");
        }
        names.extend(optional.iter().map(|c| c.name()));
    }
    names.extend(supplementary);
    if sa_only && pick_column {
        names.push("PICK");
    }
    if qc_class {
        names.push("QC_CLASS");
    }
    format!("#{}", names.join("\t"))
}

/// Optional knobs that extend the base tab schema. All fields are
/// no-op by default so the hot path is unchanged unless a caller opts in.
#[derive(Default, Clone, Copy)]
pub struct TabOptions<'a> {
    /// Drop rows whose transcript's gene_id / gene_symbol is not in this
    /// set. Variant-level rows (intergenic) are dropped entirely; in
    /// `sa_only` mode the filter is ignored because there is no transcript
    /// context.
    pub gene_set: Option<&'a crate::geneset::GeneSet>,
    /// Insert an explicit `REF` column right after `Allele`.
    pub explicit_ref: bool,
    /// Append a `QC_CLASS` column carrying the supplied class name. When
    /// `None`, the column is omitted.
    pub qc_class: Option<&'a str>,
    /// Write a `PICK` column carrying `1` on the flagged (transcript, allele)
    /// pairs and `-` on the rest, straight after `FLAGS` as VEP does. Set by
    /// `--flag-pick*`; omitted otherwise, for the reason in
    /// [`csq_fields_with_pick`].
    pub pick_column: bool,
    /// VEP's optional columns in `FLAG_FIELDS` order ([`ColumnFlags::columns`]):
    /// `REF_ALLELE` before `IMPACT`, the rest after `FLAGS` (and `PICK`), and
    /// before the columns of loaded supplementary sources.
    pub optional: &'a [OptionalColumn],
}

/// `1` when flagged, `-` when not, so the column is never blank - every other
/// empty cell in this format is `-` too.
fn pick_cell(flagged: bool) -> String {
    if flagged { "1" } else { "-" }.to_string()
}

/// Format a VariationFeature as a tab-delimited VEP output line using the
/// default schema. Equivalent to `format_tab_line_with(vf, specs, sa_only,
/// TabOptions::default())`.
///
/// `specs` enumerates which supplementary annotation sources are loaded for
/// this run (build once via `LoadedSupplementarySpecs::new`). For each
/// loaded source, one extra tab column is appended carrying the same
/// pipe-delimited value emitted into the corresponding VCF `FV_*` INFO
/// field (or the `SpliceAI` INFO field for the SpliceAI source). The
/// VCF-side `<ID>=` prefix is *not* included — the column header line
/// already names each column. The column order matches
/// `tab_supplementary_column_names`. Missing values render as `-`.
pub fn format_tab_line(
    vf: &VariationFeature,
    specs: &LoadedSupplementarySpecs,
    sa_only: bool,
) -> Vec<String> {
    format_tab_line_with(vf, specs, sa_only, TabOptions::default())
}

/// Tab formatter with optional schema extensions. See [`TabOptions`].
pub fn format_tab_line_with(
    vf: &VariationFeature,
    specs: &LoadedSupplementarySpecs,
    sa_only: bool,
    opts: TabOptions<'_>,
) -> Vec<String> {
    let mut lines = Vec::new();

    // `format_coords` of VEP's `Utils.pm`: the lower coordinate first. An
    // insertion has its end below its start (it sits between two bases), and VEP
    // writes it as `chr10:11950-11951`, not `chr10:11951-11950`.
    let (low, high) = (
        vf.position.start.min(vf.position.end),
        vf.position.start.max(vf.position.end),
    );
    let location = if low == high {
        format!("{}:{}", vf.position.chromosome, low)
    } else {
        format!("{}:{}-{}", vf.position.chromosome, low, high)
    };

    let uploaded_variation = uploaded_variation(vf);

    let extra_count = specs.column_count();
    let ref_str: Option<String> = if opts.explicit_ref {
        Some(vf.ref_allele.to_string())
    } else {
        None
    };
    let qc_extra = if opts.qc_class.is_some() { 1 } else { 0 };
    let ref_extra = if opts.explicit_ref { 1 } else { 0 };
    let pick_extra = usize::from(opts.pick_column);

    if sa_only {
        for tv in &vf.transcript_variations {
            for aa in &tv.allele_annotations {
                let mut parts: Vec<String> =
                    Vec::with_capacity(3 + ref_extra + extra_count + qc_extra + pick_extra);
                parts.push(uploaded_variation.clone());
                parts.push(location.clone());
                parts.push(aa.allele.to_string());
                if let Some(ref r) = ref_str {
                    parts.push(r.clone());
                }
                if extra_count > 0 {
                    parts.extend(format_supplementary_tab_columns_for_allele(
                        vf, tv, aa, specs,
                    ));
                }
                if opts.pick_column {
                    parts.push(pick_cell(aa.pick));
                }
                if let Some(cls) = opts.qc_class {
                    parts.push(cls.to_string());
                }
                lines.push(parts.join("\t"));
            }
        }
        // Defensive fallback: if transcript_variations is empty (the scaffold
        // should prevent this, but guard against direct callers), still emit
        // one row per alt allele padded with `-` so column count matches the
        // header.
        if lines.is_empty() {
            for alt in &vf.alt_alleles {
                let mut parts: Vec<String> =
                    Vec::with_capacity(3 + ref_extra + extra_count + qc_extra + pick_extra);
                parts.push(uploaded_variation.clone());
                parts.push(location.clone());
                parts.push(alt.to_string());
                if let Some(ref r) = ref_str {
                    parts.push(r.clone());
                }
                if extra_count > 0 {
                    parts.extend(vec!["-".to_string(); extra_count]);
                }
                if opts.pick_column {
                    parts.push(pick_cell(false));
                }
                if let Some(cls) = opts.qc_class {
                    parts.push(cls.to_string());
                }
                lines.push(parts.join("\t"));
            }
        }
        return lines;
    }

    let row_capacity = 17 + ref_extra + opts.optional.len() + extra_count + qc_extra + pick_extra;

    for tv in &vf.transcript_variations {
        if let Some(set) = opts.gene_set {
            if !set.contains_gene(&tv.gene_id, tv.gene_symbol.as_deref()) {
                continue;
            }
        }
        for aa in &tv.allele_annotations {
            let consequence_str = aa
                .consequences
                .iter()
                .map(|c| c.so_term())
                .collect::<Vec<_>>()
                .join(",");

            let impact_str = aa.impact.as_str().to_string();
            let distance_str = aa
                .distance
                .map(|d| d.to_string())
                .unwrap_or_else(|| "-".to_string());
            // No feature, no strand: VEP prints `-`, not the placeholder's `1`.
            let featureless = tv.is_featureless();
            let strand_str = if featureless {
                "-".to_string()
            } else {
                format!("{}", tv.strand.as_int())
            };
            // VEP's FLAGS are the transcript's own (`cds_start_NF`, `cds_end_NF`);
            // whether it is canonical is the CANONICAL column. This wrote the word
            // `canonical` here and never the flags.
            let flags_str = if tv.flags.is_empty() {
                "-".to_string()
            } else {
                tv.flags.join(",")
            };

            let mut parts: Vec<String> = Vec::with_capacity(row_capacity);
            parts.push(uploaded_variation.clone());
            parts.push(location.clone());
            parts.push(aa.allele.to_string());
            if let Some(ref r) = ref_str {
                parts.push(r.clone());
            }
            parts.push(tv.gene_id.to_string());
            parts.push(tv.transcript_id.to_string());
            parts.push(if featureless { "-" } else { "Transcript" }.to_string());
            parts.push(consequence_str);
            parts.push(span_cell(aa.location.cdna));
            parts.push(span_cell(aa.location.cds));
            parts.push(span_cell(aa.location.protein));
            parts.push(
                aa.amino_acids
                    .as_ref()
                    .map(|(r, a)| amino_acids_text(r, a))
                    .unwrap_or_else(|| "-".to_string()),
            );
            parts.push(
                aa.codons
                    .as_ref()
                    .map(|(r, a)| format!("{}/{}", r, a))
                    .unwrap_or_else(|| "-".to_string()),
            );
            parts.push(if aa.existing_variation.is_empty() {
                "-".to_string()
            } else {
                aa.existing_variation.join(",")
            });
            // VEP's order (`FLAG_FIELDS`): REF_ALLELE, then IMPACT to FLAGS, then
            // PICK, then the other optional columns.
            let mut optional = opts.optional;
            if let [OptionalColumn::RefAllele, rest @ ..] = optional {
                parts.push(OptionalColumn::RefAllele.cell(vf, tv, aa));
                optional = rest;
            }
            parts.push(impact_str);
            parts.push(distance_str);
            parts.push(strand_str);
            parts.push(flags_str);
            if opts.pick_column {
                parts.push(pick_cell(aa.pick));
            }
            for col in optional {
                parts.push(col.cell(vf, tv, aa));
            }

            if extra_count > 0 {
                parts.extend(format_supplementary_tab_columns_for_allele(
                    vf, tv, aa, specs,
                ));
            }

            if let Some(cls) = opts.qc_class {
                parts.push(cls.to_string());
            }

            lines.push(parts.join("\t"));
        }
    }

    // If no transcript annotations, still output the variant with intergenic
    // — unless a gene-set filter is active, in which case an unannotated
    // variant has nothing to match and is dropped.
    if vf.transcript_variations.is_empty() && opts.gene_set.is_none() {
        for alt in &vf.alt_alleles {
            let mut parts: Vec<String> = Vec::with_capacity(row_capacity);
            parts.push(uploaded_variation.clone());
            parts.push(location.clone());
            parts.push(alt.to_string());
            if let Some(ref r) = ref_str {
                parts.push(r.clone());
            }
            parts.extend(vec!["-".to_string(); 3]);
            parts.push(Consequence::IntergenicVariant.so_term().to_string());
            // cDNA, CDS and protein position, amino acids, codons, existing variation
            parts.extend(vec!["-".to_string(); 6]);
            let mut optional = opts.optional;
            if let [OptionalColumn::RefAllele, rest @ ..] = optional {
                parts.push(OptionalColumn::RefAllele.variant_cell(vf, None));
                optional = rest;
            }
            // An intergenic variant is MODIFIER, as the placeholder row says; the
            // cells after it (DISTANCE, STRAND, FLAGS) belong to a feature.
            parts.push(Impact::Modifier.as_str().to_string());
            parts.extend(vec!["-".to_string(); 3]);
            if opts.pick_column {
                parts.push(pick_cell(false));
            }
            for col in optional {
                parts.push(col.variant_cell(vf, None));
            }
            if extra_count > 0 {
                parts.extend(vec!["-".to_string(); extra_count]);
            }
            if let Some(cls) = opts.qc_class {
                parts.push(cls.to_string());
            }
            lines.push(parts.join("\t"));
        }
    }

    lines
}

/// Format a VariationFeature as JSON.
///
/// When `sa_only` is true (used by `--sa-only` mode), the
/// `transcript_consequences`, `most_severe_consequence`, and per-allele
/// ACMG metadata blocks are omitted, and per-allele supplementary
/// annotations are surfaced under a top-level `alleles` array
/// (`[{"allele": .., "<sa_key>": <value>, ...}]`). Variant-level
/// `supplementary_annotations` and `gene_annotations` (if any) still emit.
pub fn format_json(vf: &VariationFeature, sa_only: bool) -> serde_json::Value {
    let mut obj = serde_json::Map::new();

    obj.insert("id".into(), json_str(&vf.variation_name));
    obj.insert(
        "seq_region_name".into(),
        serde_json::Value::String(vf.position.chromosome.clone()),
    );
    obj.insert(
        "start".into(),
        serde_json::Value::Number(vf.position.start.into()),
    );
    obj.insert(
        "end".into(),
        serde_json::Value::Number(vf.position.end.into()),
    );
    obj.insert(
        "allele_string".into(),
        serde_json::Value::String(vf.allele_string.clone()),
    );
    obj.insert(
        "strand".into(),
        serde_json::Value::Number(vf.position.strand.as_int().into()),
    );

    if sa_only {
        let alleles: Vec<serde_json::Value> = vf
            .transcript_variations
            .iter()
            .flat_map(|tv| {
                tv.allele_annotations.iter().map(|aa| {
                    let mut a = serde_json::Map::new();
                    a.insert(
                        "allele".into(),
                        serde_json::Value::String(aa.allele.to_string()),
                    );
                    for (key, json_str) in &aa.supplementary {
                        if let Ok(val) = serde_json::from_str::<serde_json::Value>(json_str) {
                            a.insert(key.clone(), val);
                        }
                    }
                    serde_json::Value::Object(a)
                })
            })
            .collect();
        if !alleles.is_empty() {
            obj.insert("alleles".into(), serde_json::Value::Array(alleles));
        }

        for sa in &vf.supplementary_annotations {
            if let Ok(val) = serde_json::from_str::<serde_json::Value>(&sa.json_string) {
                obj.insert(sa.json_key.clone(), val);
            }
        }
        if !vf.gene_annotations.is_empty() {
            let mut genes_map = serde_json::Map::new();
            for ga in &vf.gene_annotations {
                if let Ok(val) = serde_json::from_str::<serde_json::Value>(&ga.json_string) {
                    genes_map
                        .entry(ga.gene_symbol.clone())
                        .or_insert_with(|| serde_json::Value::Object(serde_json::Map::new()))
                        .as_object_mut()
                        .map(|obj| obj.insert(ga.json_key.clone(), val));
                }
            }
            if !genes_map.is_empty() {
                obj.insert("genes".into(), serde_json::Value::Object(genes_map));
            }
        }

        return serde_json::Value::Object(obj);
    }

    if let Some(ref msq) = vf.most_severe_consequence {
        obj.insert(
            "most_severe_consequence".into(),
            serde_json::Value::String(msq.so_term().to_string()),
        );
    }

    let transcript_consequences: Vec<serde_json::Value> = vf
        .transcript_variations
        .iter()
        .flat_map(|tv| {
            tv.allele_annotations.iter().map(move |aa| {
                let mut tc = serde_json::Map::new();
                tc.insert(
                    "gene_id".into(),
                    serde_json::Value::String(tv.gene_id.to_string()),
                );
                tc.insert(
                    "transcript_id".into(),
                    serde_json::Value::String(tv.transcript_id.to_string()),
                );
                tc.insert(
                    "biotype".into(),
                    serde_json::Value::String(tv.biotype.to_string()),
                );
                if let Some(ref sym) = tv.gene_symbol {
                    tc.insert(
                        "gene_symbol".into(),
                        serde_json::Value::String(sym.to_string()),
                    );
                }
                tc.insert(
                    "consequence_terms".into(),
                    serde_json::Value::Array(
                        aa.consequences
                            .iter()
                            .map(|c| serde_json::Value::String(c.so_term().to_string()))
                            .collect(),
                    ),
                );
                tc.insert(
                    "impact".into(),
                    serde_json::Value::String(aa.impact.as_str().to_string()),
                );
                tc.insert(
                    "variant_allele".into(),
                    serde_json::Value::String(aa.allele.to_string()),
                );
                tc.insert(
                    "strand".into(),
                    serde_json::Value::Number(tv.strand.as_int().into()),
                );
                if tv.canonical {
                    tc.insert("canonical".into(), serde_json::Value::Number(1.into()));
                }
                // Present only when set, like `canonical` beside it and like
                // VEP's own JSON: a `"pick": 0` on every consequence of every
                // run would say a pick happened and this one lost, which is
                // false for a run that never picked.
                if aa.pick {
                    tc.insert("pick".into(), serde_json::Value::Number(1.into()));
                }
                if let Some(ref ms) = tv.mane_select {
                    tc.insert("mane_select".into(), serde_json::Value::String(ms.clone()));
                }
                if let Some(ref mpc) = tv.mane_plus_clinical {
                    tc.insert(
                        "mane_plus_clinical".into(),
                        serde_json::Value::String(mpc.clone()),
                    );
                }
                if let Some(t) = tv.tsl {
                    tc.insert("tsl".into(), serde_json::Value::Number(t.into()));
                }
                if let Some(ref a) = tv.appris {
                    tc.insert("appris".into(), serde_json::Value::String(a.clone()));
                }
                if let Some(ref c) = tv.ccds {
                    tc.insert("ccds".into(), serde_json::Value::String(c.clone()));
                }
                if tv.gencode_primary {
                    tc.insert(
                        "gencode_primary".into(),
                        serde_json::Value::Number(1.into()),
                    );
                }
                if let Some(ref pid) = tv.protein_id {
                    tc.insert("protein_id".into(), serde_json::Value::String(pid.clone()));
                }
                // VEP's JSON splits the position cell and leaves out an end that is `?`
                for (prefix, span) in [
                    ("cdna", aa.location.cdna),
                    ("cds", aa.location.cds),
                    ("protein", aa.location.protein),
                ] {
                    let (s, e) = span.json_bounds();
                    if let Some(s) = s {
                        tc.insert(
                            format!("{prefix}_start"),
                            serde_json::Value::Number(s.into()),
                        );
                    }
                    if let Some(e) = e {
                        tc.insert(format!("{prefix}_end"), serde_json::Value::Number(e.into()));
                    }
                }
                if let Some(ref aas) = aa.amino_acids {
                    tc.insert(
                        "amino_acids".into(),
                        serde_json::Value::String(amino_acids_text(&aas.0, &aas.1)),
                    );
                }
                if let Some(ref cdns) = aa.codons {
                    tc.insert(
                        "codons".into(),
                        serde_json::Value::String(format!("{}/{}", cdns.0, cdns.1)),
                    );
                }
                for (key, ordinals) in [("exon", aa.location.exon), ("intron", aa.location.intron)]
                {
                    if let Some(o) = ordinals {
                        let mut cell = String::new();
                        o.write_to(&mut cell);
                        tc.insert(key.into(), serde_json::Value::String(cell));
                    }
                }
                if let Some(ref h) = aa.hgvsg {
                    tc.insert("hgvsg".into(), serde_json::Value::String(h.clone()));
                }
                if let Some(ref h) = aa.hgvsc {
                    tc.insert("hgvsc".into(), serde_json::Value::String(h.clone()));
                }
                if let Some(ref h) = aa.hgvsp {
                    tc.insert("hgvsp".into(), serde_json::Value::String(h.clone()));
                }
                if let Some(d) = aa.distance {
                    tc.insert("distance".into(), serde_json::Value::Number(d.into()));
                }
                // Per-allele supplementary annotations (ClinVar, gnomAD, etc.)
                for (key, json_str) in &aa.supplementary {
                    if let Ok(val) = serde_json::from_str::<serde_json::Value>(json_str) {
                        tc.insert(key.clone(), val);
                    }
                }
                // ACMG-AMP classification
                if let Some(ref acmg) = aa.acmg_classification {
                    tc.insert("acmg".into(), acmg.clone());
                }
                serde_json::Value::Object(tc)
            })
        })
        .collect();

    obj.insert(
        "transcript_consequences".into(),
        serde_json::Value::Array(transcript_consequences),
    );

    // Variant type (for SVs)
    if vf.variant_type != fastvep_core::VariantType::Unknown {
        obj.insert(
            "variant_type".into(),
            serde_json::Value::String(format!("{:?}", vf.variant_type)),
        );
    }
    if let Some(sv_end) = vf.sv_end {
        obj.insert("sv_end".into(), serde_json::Value::Number(sv_end.into()));
    }
    if let Some(sv_len) = vf.sv_len {
        obj.insert("sv_len".into(), serde_json::Value::Number(sv_len.into()));
    }

    // Supplementary annotations (pre-serialized JSON from SA providers)
    for sa in &vf.supplementary_annotations {
        if let Ok(val) = serde_json::from_str::<serde_json::Value>(&sa.json_string) {
            obj.insert(sa.json_key.clone(), val);
        }
    }

    // Gene-level annotations
    if !vf.gene_annotations.is_empty() {
        let mut genes_map = serde_json::Map::new();
        for ga in &vf.gene_annotations {
            if let Ok(val) = serde_json::from_str::<serde_json::Value>(&ga.json_string) {
                genes_map
                    .entry(ga.gene_symbol.clone())
                    .or_insert_with(|| serde_json::Value::Object(serde_json::Map::new()))
                    .as_object_mut()
                    .map(|obj| obj.insert(ga.json_key.clone(), val));
            }
        }
        if !genes_map.is_empty() {
            obj.insert("genes".into(), serde_json::Value::Object(genes_map));
        }
    }

    serde_json::Value::Object(obj)
}

/// Format a variant as Nirvana-style structured JSON.
///
/// This is a richer format with nested sections:
/// - `position`: chromosome, position coordinates
/// - `variants`: array of per-allele annotations with supplementary data
/// - `genes`: gene-level annotations keyed by symbol
pub fn format_nirvana_json(vf: &VariationFeature) -> serde_json::Value {
    let mut obj = serde_json::Map::new();

    // Position section
    obj.insert(
        "chromosome".into(),
        serde_json::Value::String(vf.position.chromosome.clone()),
    );
    obj.insert(
        "position".into(),
        serde_json::Value::Number(vf.position.start.into()),
    );
    obj.insert(
        "end".into(),
        serde_json::Value::Number(vf.position.end.into()),
    );

    let ref_str = vf.ref_allele.to_string();
    obj.insert("refAllele".into(), serde_json::Value::String(ref_str));

    let alt_strs: Vec<serde_json::Value> = vf
        .alt_alleles
        .iter()
        .map(|a| serde_json::Value::String(a.to_string()))
        .collect();
    obj.insert("altAlleles".into(), serde_json::Value::Array(alt_strs));

    if vf.variant_type != fastvep_core::VariantType::Unknown {
        obj.insert(
            "variantType".into(),
            serde_json::Value::String(format!("{:?}", vf.variant_type)),
        );
    }

    // Variants section: one per alt allele with all transcript consequences
    let mut variants = Vec::new();
    for alt in &vf.alt_alleles {
        let alt_str = alt.to_string();
        let mut var_obj = serde_json::Map::new();
        var_obj.insert(
            "altAllele".into(),
            serde_json::Value::String(alt_str.clone()),
        );

        // Collect transcript consequences for this allele
        let mut transcripts = Vec::new();
        for tv in &vf.transcript_variations {
            for aa in &tv.allele_annotations {
                if aa.allele.to_string() == alt_str {
                    let mut tc = serde_json::Map::new();
                    tc.insert(
                        "transcriptId".into(),
                        serde_json::Value::String(tv.transcript_id.to_string()),
                    );
                    tc.insert(
                        "geneId".into(),
                        serde_json::Value::String(tv.gene_id.to_string()),
                    );
                    if let Some(ref sym) = tv.gene_symbol {
                        tc.insert(
                            "geneSymbol".into(),
                            serde_json::Value::String(sym.to_string()),
                        );
                    }
                    tc.insert(
                        "biotype".into(),
                        serde_json::Value::String(tv.biotype.to_string()),
                    );
                    tc.insert(
                        "consequences".into(),
                        serde_json::Value::Array(
                            aa.consequences
                                .iter()
                                .map(|c| serde_json::Value::String(c.so_term().to_string()))
                                .collect(),
                        ),
                    );
                    tc.insert(
                        "impact".into(),
                        serde_json::Value::String(aa.impact.as_str().to_string()),
                    );
                    if tv.canonical {
                        tc.insert("isCanonical".into(), serde_json::Value::Bool(true));
                    }
                    if let Some(ref ms) = tv.mane_select {
                        tc.insert("maneSelect".into(), serde_json::Value::String(ms.clone()));
                    }
                    if let Some(ref mpc) = tv.mane_plus_clinical {
                        tc.insert(
                            "manePlusClinical".into(),
                            serde_json::Value::String(mpc.clone()),
                        );
                    }
                    if let Some(t) = tv.tsl {
                        tc.insert("tsl".into(), serde_json::Value::Number(t.into()));
                    }
                    if let Some(ref a) = tv.appris {
                        tc.insert("appris".into(), serde_json::Value::String(a.clone()));
                    }
                    if let Some(ref c) = tv.ccds {
                        tc.insert("ccds".into(), serde_json::Value::String(c.clone()));
                    }
                    if tv.gencode_primary {
                        tc.insert("isGencodePrimary".into(), serde_json::Value::Bool(true));
                    }
                    if let Some(ref pid) = tv.protein_id {
                        tc.insert("proteinId".into(), serde_json::Value::String(pid.clone()));
                    }
                    if let Some(ref h) = aa.hgvsg {
                        tc.insert("hgvsg".into(), serde_json::Value::String(h.clone()));
                    }
                    if let Some(ref h) = aa.hgvsc {
                        tc.insert("hgvsc".into(), serde_json::Value::String(h.clone()));
                    }
                    if let Some(ref h) = aa.hgvsp {
                        tc.insert("hgvsp".into(), serde_json::Value::String(h.clone()));
                    }

                    // Per-allele supplementary annotations
                    for (key, json_str) in &aa.supplementary {
                        if let Ok(val) = serde_json::from_str::<serde_json::Value>(json_str) {
                            tc.insert(key.clone(), val);
                        }
                    }

                    transcripts.push(serde_json::Value::Object(tc));
                }
            }
        }
        var_obj.insert("transcripts".into(), serde_json::Value::Array(transcripts));
        variants.push(serde_json::Value::Object(var_obj));
    }
    obj.insert("variants".into(), serde_json::Value::Array(variants));

    // Gene annotations
    if !vf.gene_annotations.is_empty() {
        let mut genes_map = serde_json::Map::new();
        for ga in &vf.gene_annotations {
            if let Ok(val) = serde_json::from_str::<serde_json::Value>(&ga.json_string) {
                genes_map
                    .entry(ga.gene_symbol.clone())
                    .or_insert_with(|| serde_json::Value::Object(serde_json::Map::new()))
                    .as_object_mut()
                    .map(|obj| obj.insert(ga.json_key.clone(), val));
            }
        }
        if !genes_map.is_empty() {
            obj.insert("genes".into(), serde_json::Value::Object(genes_map));
        }
    }

    serde_json::Value::Object(obj)
}

fn json_str(opt: &Option<String>) -> serde_json::Value {
    match opt {
        Some(s) => serde_json::Value::String(s.clone()),
        None => serde_json::Value::Null,
    }
}

#[cfg(test)]
mod tests {
    /// A `LoadedSupplementarySpecs` with no custom sources, for the tests that
    /// only exercise the built-in projections. Those walk
    /// `VCF_PROJECTION_SPECS` in full regardless of what loaded, so an empty
    /// lookup is the right neutral value.
    fn no_custom_sources() -> super::LoadedSupplementarySpecs {
        super::LoadedSupplementarySpecs::new(&[], &[])
    }

    use super::*;
    use crate::variant::{AlleleAnnotation, TranscriptVariation, VariationFeature, VcfFields};
    use fastvep_core::{
        Allele, Consequence, GeneAnnotation, GenomicPosition, Impact, Strand,
        SupplementaryAnnotation, TranscriptLocation, VariantType,
    };
    use std::sync::Arc;

    /// Every default column resolves to its own writer, and the five that do
    /// not are named.
    ///
    /// `CsqField::from_name` is a second copy of the column vocabulary, so a
    /// typo in it - `"HGVSc" => Self::Hgvsg` - would silently write the wrong
    /// value into a column, or blank it. Nothing else in the suite would fail
    /// for a column no test happens to assert.
    #[test]
    fn every_default_csq_column_resolves_to_a_distinct_writer() {
        // fastVEP emits the regulatory columns' headers for VEP header
        // compatibility but has no writer for them, so they render empty.
        const NO_WRITER: &[&str] = &[
            "MOTIF_NAME",
            "MOTIF_POS",
            "HIGH_INF_POS",
            "MOTIF_SCORE_CHANGE",
            "TRANSCRIPTION_FACTORS",
        ];

        let unresolved: Vec<&str> = DEFAULT_CSQ_FIELDS
            .iter()
            .copied()
            .filter(|name| CsqField::from_name(name) == CsqField::Unknown)
            .collect();
        assert_eq!(
            unresolved, NO_WRITER,
            "a default column stopped resolving to a writer, or gained one"
        );

        let mut seen: Vec<(CsqField, &str)> = Vec::new();
        for name in DEFAULT_CSQ_FIELDS.iter().copied() {
            let field = CsqField::from_name(name);
            if field == CsqField::Unknown {
                continue;
            }
            if let Some((_, other)) = seen.iter().find(|(f, _)| *f == field) {
                panic!("{name} and {other} both resolve to {field:?}");
            }
            seen.push((field, name));
        }
        assert_eq!(seen.len(), DEFAULT_CSQ_FIELDS.len() - NO_WRITER.len());
    }

    #[test]
    fn test_escape_csq_value() {
        assert_eq!(escape_csq_value("hello,world"), "hello&world");
        assert_eq!(escape_csq_value("a;b"), "a%3Bb");
        assert_eq!(escape_csq_value("a|b"), "a&b");
        assert_eq!(escape_csq_value("a b"), "a_b");
        assert_eq!(escape_csq_value("p.Leu153="), "p.Leu153%3D");
        // A space is the only substitution the byte pre-scan had not been
        // told about, so a value whose *only* special character is a space
        // took the copy-whole path and reached the INFO column unescaped.
        assert_eq!(escape_csq_value("Breast cancer"), "Breast_cancer");
        // Left alone: VEP writes it, and every HGVSc has one.
        assert_eq!(
            escape_csq_value("ENST00000370131.3:c.452_454del"),
            "ENST00000370131.3:c.452_454del"
        );
        // These do not mangle a value, they end the column or the record.
        assert_eq!(escape_csq_value("a\tb"), "a%09b");
        assert_eq!(escape_csq_value("a\r\nb"), "a%0D%0Ab");
    }

    #[test]
    fn test_escape_vcf_subfield_encodes_only_structure() {
        // Illegal in an INFO value, or one of the pipe field's own delimiters.
        assert_eq!(escape_vcf_subfield("a;b"), "a%3Bb");
        assert_eq!(escape_vcf_subfield("a=b"), "a%3Db");
        assert_eq!(escape_vcf_subfield("a,b"), "a%2Cb");
        assert_eq!(escape_vcf_subfield("a|b"), "a%7Cb");
        assert_eq!(escape_vcf_subfield("a&b"), "a%26b");
        assert_eq!(escape_vcf_subfield("a%b"), "a%25b");
        assert_eq!(escape_vcf_subfield("a\tb"), "a%09b");

        // Reserved by nothing, and VEP 115.1 writes it through literally, as
        // fastVEP's own CSQ column already did - every HGVSc carries one. The
        // single character #123 is about.
        assert_eq!(
            escape_vcf_subfield("NM_001159702.3:c.-101+5113dup"),
            "NM_001159702.3:c.-101+5113dup"
        );

        // Legal unencoded, kept encoded on purpose. VEP writes a literal space
        // as `_` and a literal `"` as itself; neither helps here, because no
        // `FV_*` field exists in VEP to be diffed against, and both cost a
        // reader of the column - `awk` splits on the space, and R's
        // `read.delim` and Python's `csv` treat the quote as an opening quote.
        assert_eq!(escape_vcf_subfield("Breast cancer"), "Breast%20cancer");
        assert_eq!(escape_vcf_subfield("6\" of tissue"), "6%22%20of%20tissue");
    }

    /// `push_clinvar_protein_records` as an expression. A wrapper, not a
    /// second implementation: the production path writes into a buffer it
    /// already owns, and a test that re-implemented that would stop testing
    /// it.
    fn rendered_records(parsed: &serde_json::Value) -> String {
        let mut out = String::new();
        push_clinvar_protein_records(&mut out, parsed);
        out
    }

    #[test]
    fn test_clinvar_protein_record_keeps_its_delimiters_literal() {
        let parsed: serde_json::Value = serde_json::from_str(
            r#"{"proteinVariants":[{"pos":175,"refAa":"R","altAa":"H","sig":"Pathogenic"},
                                   {"pos":248,"refAa":"R","altAa":"W","sig":"Likely_pathogenic"}]}"#,
        )
        .unwrap();
        assert_eq!(
            rendered_records(&parsed),
            "175:R>H:Pathogenic&248:R>W:Likely_pathogenic"
        );
    }

    #[test]
    fn test_clinvar_protein_leaf_cannot_forge_a_delimiter() {
        // No source reaches this today - `sig` is one of four constants and the
        // residues are single letters. If one ever does, the record still
        // splits into the three leaves a reader expects rather than five.
        let parsed: serde_json::Value = serde_json::from_str(
            r#"{"proteinVariants":[{"pos":175,"refAa":"R","altAa":"H","sig":"odd:term>here"}]}"#,
        )
        .unwrap();
        let rendered = rendered_records(&parsed);
        assert_eq!(rendered, "175:R>H:odd%3Aterm%3Ehere");

        let leaves: Vec<&str> = rendered.split(':').collect();
        assert_eq!(leaves, vec!["175", "R>H", "odd%3Aterm%3Ehere"]);
        assert_eq!(leaves[1].split('>').count(), 2);
    }

    #[test]
    fn test_csq_header() {
        let header = csq_header_line(&["Allele", "Consequence"]);
        assert!(header.contains("Format: Allele|Consequence"));
    }

    #[test]
    fn test_span_cell() {
        let cell = |s, e| span_cell(TranscriptSpan::new(s, e));
        assert_eq!(cell(Some(100), Some(100)), "100");
        assert_eq!(cell(Some(100), Some(200)), "100-200");
        assert_eq!(cell(None, None), "-");
        // a span with an end outside the transcript, and an insertion
        assert_eq!(cell(None, Some(7)), "?-7");
        assert_eq!(cell(Some(601), None), "601-?");
        assert_eq!(cell(Some(602), Some(601)), "601-602");
    }

    fn projection_test_variant() -> VariationFeature {
        let supplementary = vec![
            ("clinvar".into(), r#"{"significance":["Pathogenic","Likely_pathogenic"],"reviewStatus":"criteria_provided,_multiple_submitters,_no_conflicts","phenotypes":["Breast,cancer","Ovarian|cancer"],"variantClass":"SNV","soAccession":"SO:0001483"}"#.into()),
            ("gnomad".into(), r#"{"allAf":1.2e-4,"allAc":12,"allAn":100000,"allHc":0,"afrAf":2.1e-4,"nfeAf":9.0e-5}"#.into()),
            ("dbsnp".into(), r#"{"id":"rs123","globalMaf":0.042}"#.into()),
            ("cosmic".into(), r#"{"id":"COSV123","gene":"GENE1","count":7}"#.into()),
            ("oneKg".into(), r#"{"allAf":0.01,"afrAf":0.02,"amrAf":0.03,"easAf":0.04,"eurAf":0.05,"sasAf":0.06}"#.into()),
            ("topmed".into(), r#"{"allAf":0.001,"allAc":4,"allAn":20000}"#.into()),
            ("mitomap".into(), r#"{"disease":"MELAS;like","status":"Reported"}"#.into()),
            ("phylop".into(), "3.14".into()),
            ("gerp".into(), "-1.5".into()),
            ("dann".into(), "0.99".into()),
            ("revel".into(), r#"{"score":0.8123}"#.into()),
            ("spliceAI".into(), r#"{"gene":"GENE|1","dsAg":0.01,"dsAl":0.0,"dsDg":0.85,"dsDl":0.0,"dpAg":5,"dpAl":-28,"dpDg":2,"dpDl":-13}"#.into()),
            ("primateAI".into(), r#"{"score":0.4567}"#.into()),
            ("dbnsfp".into(), r#"{"sift":"deleterious(0.010)","polyphen":"probably_damaging(0.980)"}"#.into()),
        ];

        VariationFeature {
            position: GenomicPosition::new("1", 25000, 25000, Strand::Forward),
            allele_string: "A/G".into(),
            ref_allele: Allele::from_str("A"),
            alt_alleles: vec![Allele::from_str("G")],
            variation_name: None,
            vcf_fields: None,
            transcript_variations: vec![TranscriptVariation {
                transcript_id: Arc::from("TX1"),
                gene_id: Arc::from("GENE1"),
                gene_symbol: Some(Arc::from("GENE1")),
                biotype: Arc::from("protein_coding"),
                allele_annotations: vec![AlleleAnnotation {
                    allele: Allele::from_str("G"),
                    consequences: vec![Consequence::MissenseVariant],
                    impact: Impact::Moderate,
                    location: TranscriptLocation::default(),
                    protein_position: None,
                    amino_acids: None,
                    codons: None,
                    exon: None,
                    splice: Default::default(),
                    distance: None,
                    protein_length: None,
                    escapes_nmd: None,
                    hgvsc: None,
                    hgvsp: None,
                    hgvsg: None,
                    hgvs_offset: None,
                    existing_variation: Vec::new(),
                    sift: None,
                    polyphen: None,
                    supplementary,
                    acmg_classification: None,
                    pick: false,
                }],
                canonical: false,
                strand: Strand::Forward,
                source: None,
                protein_id: None,
                mane_select: None,
                mane_plus_clinical: None,
                tsl: None,
                appris: None,
                ccds: None,
                gencode_primary: false,
                symbol_source: None,
                hgnc_id: None,
                flags: Vec::new(),
            }],
            existing_variants: Vec::new(),
            minimised: false,
            most_severe_consequence: None,
            variant_type: VariantType::Snv,
            sv_end: None,
            sv_len: None,
            supplementary_annotations: vec![SupplementaryAnnotation {
                json_key: "customVariant".into(),
                is_array: false,
                json_string: r#"{"note":"top-level"}"#.into(),
            }],
            gene_annotations: vec![
                GeneAnnotation {
                    gene_symbol: "GENE1".into(),
                    json_key: "omim".into(),
                    json_string: r#"{"mimNumber":113705,"phenotypes":["Breast cancer","Ovarian,cancer"]}"#.into(),
                },
                GeneAnnotation {
                    gene_symbol: "GENE1".into(),
                    json_key: "gnomad_genes".into(),
                    json_string: r#"{"pLI":1.0,"loeuf":0.03,"misZ":3.45,"synZ":0.12}"#.into(),
                },
                GeneAnnotation {
                    gene_symbol: "GENE1".into(),
                    json_key: "clinvar_protein".into(),
                    json_string: r#"{"proteinVariants":[{"pos":175,"refAa":"R","altAa":"H","sig":"Pathogenic"}]}"#.into(),
                },
            ],
        }
    }

    /// A variant carrying only user-supplied custom payloads, so nothing below
    /// can be satisfied by a built-in projection.
    fn custom_source_variant(
        allele_payloads: &[(&str, &str)],
        gene_payloads: &[(&str, &str)],
    ) -> VariationFeature {
        let mut vf = projection_test_variant();
        vf.supplementary_annotations.clear();
        vf.gene_annotations = gene_payloads
            .iter()
            .map(|(key, json)| GeneAnnotation {
                gene_symbol: "GENE1".into(),
                json_key: (*key).into(),
                json_string: (*json).into(),
            })
            .collect();
        for tv in &mut vf.transcript_variations {
            for aa in &mut tv.allele_annotations {
                aa.supplementary = allele_payloads
                    .iter()
                    .map(|(key, json)| ((*key).to_string(), (*json).to_string()))
                    .collect();
            }
        }
        vf
    }

    #[test]
    fn format_spliceai_info_returns_the_spliceai_projection_alone() {
        // Public API with no in-tree caller, so nothing else would notice if
        // it broke. It used to build every projection and then filter for the
        // one entry it wanted; it now asks for that entry directly, which is
        // the same value `format_supplementary_vcf_info` files under
        // `SpliceAI`.
        let vf = projection_test_variant();
        let direct = format_spliceai_info(&vf).expect("SpliceAI payload is in the fixture");
        assert_eq!(direct, "G|GENE%7C1|0.01|0.00|0.85|0.00|5|-28|2|-13");

        let via_all = format_supplementary_vcf_info(&vf, &no_custom_sources())
            .into_iter()
            .find_map(|(id, value)| (id == "SpliceAI").then_some(value));
        assert_eq!(via_all, Some(direct));

        // No payload, no value - not an empty string.
        let bare = custom_source_variant(&[], &[]);
        assert_eq!(format_spliceai_info(&bare), None);
    }

    #[test]
    fn a_custom_source_reaches_the_vcf_and_the_tab_column() {
        // The bug in #116: a `custom_vcf` database annotated `--output-format
        // json` and nothing else, because the VCF and tab writers projected
        // through a static table that could not name it. The annotation was
        // dropped with no header, no column and no warning - a complete-looking
        // VCF missing the field the run was for.
        let vf = custom_source_variant(
            &[("mypanel", r#"{"AC":"7","AF":"0.005814","AN":"1204"}"#)],
            &[],
        );
        let specs = LoadedSupplementarySpecs::new(&["mypanel".to_string()], &[]);

        let projections = format_supplementary_vcf_info(&vf, &specs);
        assert_eq!(
            projections,
            vec![(
                "FV_MYPANEL".to_string(),
                "G|AC=7&AF=0.005814&AN=1204".to_string()
            )]
        );

        // ... and the same value in the tab column, without the `ID=` prefix.
        let columns = tab_supplementary_column_names(&specs);
        assert_eq!(columns, vec!["FV_MYPANEL"]);
        let rows = format_tab_line(&vf, &specs, false);
        assert!(
            rows.iter()
                .all(|row| row.ends_with("\tG|AC=7&AF=0.005814&AN=1204")),
            "tab rows lost the custom column: {rows:?}"
        );
    }

    #[test]
    fn every_emitted_custom_field_is_declared_in_the_header() {
        // A `FV_*` field with no `##INFO=` line is an invalid VCF that most
        // readers accept anyway, so the mismatch would be silent. The header
        // list and the value list are now derived from one
        // `LoadedSupplementarySpecs`; this is what that buys.
        let vf = custom_source_variant(
            &[
                ("panel one", r#"{"AF":"0.01"}"#),
                ("panel/one", r#"{"AF":"0.02"}"#),
            ],
            &[("my_genes", r#"{"tier":"1"}"#)],
        );
        let specs = LoadedSupplementarySpecs::new(
            &["panel one".to_string(), "panel/one".to_string()],
            &["my_genes".to_string()],
        );

        let declared: Vec<String> = vcf_info_header_lines(&specs, DEFAULT_CSQ_FIELDS, false)
            .iter()
            .filter_map(|line| vcf_info_header_id(line).map(str::to_string))
            .collect();
        for (id, _) in format_supplementary_vcf_info(&vf, &specs) {
            assert!(
                declared.contains(&id),
                "{id} was emitted but never declared; declared = {declared:?}"
            );
        }
        // And fastVEP claims ownership of each, so re-annotating replaces
        // rather than duplicates them.
        let owned = vcf_owned_info_ids(&specs);
        for id in &declared {
            assert!(
                owned.contains(&id.as_str()),
                "{id} is declared but not owned"
            );
        }
    }

    #[test]
    fn a_separator_inside_a_custom_value_cannot_break_the_pair_list() {
        // The pair list is only unambiguous because `escape_vcf_subfield`
        // percent-encodes every separator it uses. A raw `=` or `&` in a value
        // would silently invent a key or truncate one.
        let vf = custom_source_variant(&[("p", r#"{"NOTE":"a=b&c|d,e;f","OK":"1"}"#)], &[]);
        let specs = LoadedSupplementarySpecs::new(&["p".to_string()], &[]);
        let (_, value) = format_supplementary_vcf_info(&vf, &specs)
            .into_iter()
            .next()
            .expect("custom projection missing");
        assert_eq!(value, "G|NOTE=a%3Db%26c%7Cd%2Ce%3Bf&OK=1");

        // Splitting the way the header documents recovers exactly two fields
        // and exactly two pairs.
        let (_allele, pairs) = value.split_once('|').unwrap();
        let keys: Vec<&str> = pairs
            .split('&')
            .map(|pair| pair.split_once('=').unwrap().0)
            .collect();
        assert_eq!(keys, vec!["NOTE", "OK"]);
    }

    #[test]
    fn a_custom_record_with_no_fields_still_reports_the_hit() {
        // A three-column BED carries no attributes at all. "In a listed
        // region" is the whole annotation, and dropping the entry would be
        // indistinguishable from a miss.
        let vf = custom_source_variant(&[("regions", "{}")], &[]);
        let specs = LoadedSupplementarySpecs::new(&["regions".to_string()], &[]);
        assert_eq!(
            format_supplementary_vcf_info(&vf, &specs),
            vec![("FV_REGIONS".to_string(), "G|".to_string())]
        );
    }

    #[test]
    fn a_gene_level_custom_source_projects_on_the_symbol() {
        let vf = custom_source_variant(&[], &[("my_genes", r#"{"panel":"cardio","tier":2}"#)]);
        let specs = LoadedSupplementarySpecs::new(&[], &["my_genes".to_string()]);
        assert_eq!(
            format_supplementary_vcf_info(&vf, &specs),
            vec![(
                "FV_MY_GENES".to_string(),
                "GENE1|panel=cardio&tier=2".to_string()
            )]
        );
        // The header says which subject the leading field holds.
        let header = vcf_info_header_lines(&specs, DEFAULT_CSQ_FIELDS, false)
            .into_iter()
            .find(|line| line.contains("FV_MY_GENES"))
            .unwrap();
        assert!(header.contains("Format: SYMBOL|KEY=VALUE"), "{header}");
    }

    #[test]
    fn a_custom_array_value_repeats_its_key() {
        // `&` already separates pairs, so joining the elements the way
        // `json_value_to_vcf` does would leave the second element keyless.
        let vf = custom_source_variant(&[("p", r#"{"AF":[0.1,0.2],"ID":"x"}"#)], &[]);
        let specs = LoadedSupplementarySpecs::new(&["p".to_string()], &[]);
        let (_, value) = format_supplementary_vcf_info(&vf, &specs)
            .into_iter()
            .next()
            .unwrap();
        assert_eq!(value, "G|AF=0.1&AF=0.2&ID=x");
    }

    #[test]
    fn custom_info_ids_are_valid_vcf_and_never_collide() {
        // `--name` is free text; a VCF INFO ID is not. Anything outside
        // `[0-9A-Za-z_.]` becomes `_`, which means distinct names can sanitise
        // alike - and two INFO fields sharing an ID is a malformed VCF.
        assert_eq!(custom_info_id("my-panel"), "FV_MY_PANEL");
        assert_eq!(custom_info_id("my panel!"), "FV_MY_PANEL_");
        assert_eq!(custom_info_id("v1.2"), "FV_V1.2");

        let specs = LoadedSupplementarySpecs::new(
            &[
                "my-panel".to_string(),
                "my_panel".to_string(),
                "my/panel".to_string(),
                // Sanitises onto the built-in `gnomad_genes` spec's ID.
                "gnomad_gene".to_string(),
            ],
            &["gnomad_genes".to_string()],
        );
        let ids = vcf_owned_info_ids(&specs);
        let unique: std::collections::HashSet<&&str> = ids.iter().collect();
        assert_eq!(ids.len(), unique.len(), "duplicate INFO ID in {ids:?}");
        assert!(ids.contains(&"FV_MY_PANEL"));
        assert!(ids.contains(&"FV_MY_PANEL_2"));
        assert!(ids.contains(&"FV_MY_PANEL_3"));
        // The built-in keeps its ID; the custom source that wanted it yields.
        assert!(ids.contains(&"FV_GNOMAD_GENE"));
        assert!(ids.contains(&"FV_GNOMAD_GENE_2"));
    }

    #[test]
    fn a_quote_in_a_source_name_cannot_end_the_header_attribute() {
        // The name reaches the `Description="..."` attribute from `--name`,
        // and a VCF header has no escape for a quote inside one.
        let specs = LoadedSupplementarySpecs::new(&["pa\"nel\nx".to_string()], &[]);
        let header = vcf_info_header_lines(&specs, DEFAULT_CSQ_FIELDS, false)
            .into_iter()
            .find(|line| line.contains("FV_PA_NEL_X"))
            .expect("custom header missing");
        assert_eq!(
            header.matches('"').count(),
            2,
            "unbalanced quotes: {header}"
        );
        assert!(!header.contains('\n'));
    }

    #[test]
    fn vcf_projection_emits_supported_fastsa_sources_without_json_payloads() {
        let vf = projection_test_variant();
        let projections = format_supplementary_vcf_info(&vf, &no_custom_sources());
        let ids: Vec<&str> = projections.iter().map(|(id, _)| id.as_str()).collect();

        for expected in [
            "FV_CLINVAR",
            "FV_GNOMAD",
            "FV_DBSNP",
            "FV_COSMIC",
            "FV_1KG",
            "FV_TOPMED",
            "FV_MITOMAP",
            "FV_PHYLOP",
            "FV_GERP",
            "FV_DANN",
            "FV_REVEL",
            "SpliceAI",
            "FV_PRIMATEAI",
            "FV_DBNSFP",
            "FV_OMIM",
            "FV_GNOMAD_GENE",
            "FV_CLINVAR_PROTEIN",
        ] {
            assert!(
                ids.contains(&expected),
                "missing projection {expected}: {projections:?}"
            );
        }

        let info = projections
            .iter()
            .map(|(id, value)| format!("{id}={value}"))
            .collect::<Vec<_>>()
            .join(";");
        assert!(
            !info.contains('{'),
            "VCF INFO must not contain JSON objects: {info}"
        );
        assert!(
            !info.contains('}'),
            "VCF INFO must not contain JSON objects: {info}"
        );
        assert!(
            !info.contains('"'),
            "VCF INFO must not contain JSON quotes: {info}"
        );
        assert!(info.contains("SpliceAI=G|GENE%7C1|0.01|0.00|0.85|0.00|5|-28|2|-13"));
        assert!(info.contains("FV_CLINVAR=G|Pathogenic&Likely_pathogenic|criteria_provided%2C_multiple_submitters%2C_no_conflicts|Breast%2Ccancer&Ovarian%7Ccancer|SNV|SO:0001483"));
        assert!(info.contains("FV_PHYLOP=G|3.14"));
        assert!(info.contains("FV_REVEL=G|0.8123"));
        assert!(info.contains("FV_PRIMATEAI=G|0.4567"));
        assert!(info.contains("FV_OMIM=GENE1|113705|Breast%20cancer&Ovarian%2Ccancer"));
        assert!(info.contains("FV_CLINVAR_PROTEIN=GENE1|175:R>H:Pathogenic"));
    }

    #[test]
    fn vcf_info_replaces_existing_fastvep_owned_fields() {
        let vf = projection_test_variant();
        let csq = "G|missense_variant|MODERATE";
        let info = format_vcf_info_fields(
            "DP=12;CSQ=old;SpliceAI=old;FV_CLINVAR=old;KEEP=1",
            &vf,
            csq,
            &no_custom_sources(),
        );

        assert!(info.contains("DP=12"));
        assert!(info.contains("KEEP=1"));
        assert!(info.contains("CSQ=G|missense_variant|MODERATE"));
        assert!(info.contains("SpliceAI=G|GENE%7C1|0.01|0.00|0.85|0.00|5|-28|2|-13"));
        assert!(info.contains("FV_CLINVAR=G|Pathogenic&Likely_pathogenic"));
        assert!(!info.contains("CSQ=old"));
        assert!(!info.contains("SpliceAI=old"));
        assert!(!info.contains("FV_CLINVAR=old"));
    }

    fn all_supplementary_sa_keys() -> Vec<String> {
        vec![
            "clinvar".into(),
            "gnomad".into(),
            "dbsnp".into(),
            "cosmic".into(),
            "oneKg".into(),
            "topmed".into(),
            "mitomap".into(),
            "phylop".into(),
            "gerp".into(),
            "dann".into(),
            "revel".into(),
            "primateAI".into(),
            "dbnsfp".into(),
            "spliceAI".into(),
        ]
    }

    fn all_supplementary_gene_keys() -> Vec<String> {
        vec![
            "omim".into(),
            "gnomad_genes".into(),
            "clinvar_protein".into(),
        ]
    }

    #[test]
    fn tab_supplementary_column_names_match_vcf_header_order() {
        let sa_keys = all_supplementary_sa_keys();
        let gene_keys = all_supplementary_gene_keys();
        let specs = LoadedSupplementarySpecs::new(&sa_keys, &gene_keys);
        let cols = tab_supplementary_column_names(&specs);
        // SpliceAI is emitted first to mirror the VCF header ordering, then
        // every loaded VCF_PROJECTION_SPECS entry in declaration order.
        assert_eq!(cols.first().copied(), Some("SpliceAI"));
        for expected in [
            "FV_CLINVAR",
            "FV_GNOMAD",
            "FV_DBSNP",
            "FV_COSMIC",
            "FV_1KG",
            "FV_TOPMED",
            "FV_MITOMAP",
            "FV_PHYLOP",
            "FV_GERP",
            "FV_DANN",
            "FV_REVEL",
            "FV_PRIMATEAI",
            "FV_DBNSFP",
            "FV_OMIM",
            "FV_GNOMAD_GENE",
            "FV_CLINVAR_PROTEIN",
        ] {
            assert!(
                cols.contains(&expected),
                "missing column {expected}: {cols:?}"
            );
        }
    }

    #[test]
    fn tab_supplementary_header_lines_document_pipe_format() {
        let sa_keys = vec!["clinvar".into()];
        let gene_keys: Vec<String> = vec![];
        let specs = LoadedSupplementarySpecs::new(&sa_keys, &gene_keys);
        let lines = tab_supplementary_header_lines(&specs);
        let joined = lines.join("\n");
        assert!(
            joined.contains("## COLUMN=<ID=FV_CLINVAR,Description=\"fastVEP ClinVar annotations. Format: ALLELE|SIGNIFICANCE|REVIEW_STATUS|PHENOTYPES|VARIANT_CLASS|SO_ACCESSION\">"),
            "header lines should document FV_CLINVAR pipe format: {joined}"
        );
    }

    #[test]
    fn tab_line_emits_supplementary_columns_for_loaded_sources() {
        let vf = projection_test_variant();
        let sa_keys = all_supplementary_sa_keys();
        let gene_keys = all_supplementary_gene_keys();
        let specs = LoadedSupplementarySpecs::new(&sa_keys, &gene_keys);
        let lines = format_tab_line(&vf, &specs, false);
        assert_eq!(lines.len(), 1, "test variant has exactly one TV×AA row");
        let row = &lines[0];
        let cols: Vec<&str> = row.split('\t').collect();

        let extra_cols = tab_supplementary_column_names(&specs);
        assert_eq!(
            cols.len(),
            17 + extra_cols.len(),
            "tab row should carry 17 base columns plus one per loaded source"
        );

        let vcf_projections: std::collections::HashMap<String, String> =
            format_supplementary_vcf_info(&vf, &no_custom_sources())
                .into_iter()
                .collect();

        for (i, info_id) in extra_cols.iter().enumerate() {
            let tab_value = cols[17 + i];
            let expected = vcf_projections
                .get(*info_id)
                .cloned()
                .unwrap_or_else(|| "-".into());
            assert_eq!(
                tab_value, expected,
                "tab column {info_id} should equal VCF INFO value"
            );
        }

        assert!(!row.contains('{'), "tab row must not contain JSON: {row}");
        assert!(!row.contains('}'), "tab row must not contain JSON: {row}");
        assert!(
            !row.contains('"'),
            "tab row must not contain JSON quotes: {row}"
        );
    }

    #[test]
    fn tab_line_emits_dashes_when_no_supplementary_loaded_but_columns_requested() {
        let mut vf = projection_test_variant();
        for tv in &mut vf.transcript_variations {
            for aa in &mut tv.allele_annotations {
                aa.supplementary.clear();
            }
        }
        vf.gene_annotations.clear();

        let specs =
            LoadedSupplementarySpecs::new(&["clinvar".to_string(), "dbnsfp".to_string()], &[]);
        let lines = format_tab_line(&vf, &specs, false);
        let cols: Vec<&str> = lines[0].split('\t').collect();
        assert_eq!(cols.len(), 17 + 2);
        assert_eq!(cols[17], "-", "FV_CLINVAR column should be '-' when absent");
        assert_eq!(cols[18], "-", "FV_DBNSFP column should be '-' when absent");
    }

    #[test]
    fn tab_line_has_no_extra_columns_when_no_sa_loaded() {
        let vf = projection_test_variant();
        let specs = LoadedSupplementarySpecs::new(&[], &[]);
        assert!(specs.is_empty());
        let lines = format_tab_line(&vf, &specs, false);
        let cols: Vec<&str> = lines[0].split('\t').collect();
        assert_eq!(
            cols.len(),
            17,
            "without any loaded SA sources, the tab row keeps the original 17-column shape"
        );
    }

    #[test]
    fn tab_line_intergenic_row_carries_full_column_shape_when_sources_loaded() {
        // Regression test: when --sa-dir is loaded but a variant has no
        // transcript hits, the intergenic-allele row must still match the
        // header's column count so downstream parsers don't see ragged rows.
        let mut vf = projection_test_variant();
        vf.transcript_variations.clear();
        let sa_keys = vec!["clinvar".to_string(), "dbnsfp".to_string()];
        let gene_keys = vec!["omim".to_string()];
        let specs = LoadedSupplementarySpecs::new(&sa_keys, &gene_keys);

        let lines = format_tab_line(&vf, &specs, false);
        assert_eq!(
            lines.len(),
            1,
            "single alt allele yields one intergenic row"
        );
        let cols: Vec<&str> = lines[0].split('\t').collect();
        assert_eq!(
            cols.len(),
            17 + 3,
            "intergenic row must include FV_* columns"
        );
        assert_eq!(cols[6], "intergenic_variant");
        // All FV_* columns should be `-` for an intergenic row.
        for tail in &cols[17..] {
            assert_eq!(*tail, "-", "intergenic row should leave FV_* columns empty");
        }
    }

    #[test]
    fn tab_line_intergenic_row_keeps_17_columns_without_sa() {
        // Backwards-compatibility check for users not passing --sa-dir.
        let mut vf = projection_test_variant();
        vf.transcript_variations.clear();
        let specs = LoadedSupplementarySpecs::new(&[], &[]);
        let lines = format_tab_line(&vf, &specs, false);
        let cols: Vec<&str> = lines[0].split('\t').collect();
        assert_eq!(cols.len(), 17);
    }

    #[test]
    fn tab_line_emits_per_alt_values_for_multi_alt_variant() {
        // Per-tab-row dedupe: a variant with two alts and a distinct
        // supplementary payload per alt must produce two rows whose
        // FV_CLINVAR cells differ.
        use fastvep_core::Allele;

        let mut vf = projection_test_variant();
        // Replace single-alt transcript-variation with one TV per alt.
        let base_tv = vf.transcript_variations[0].clone();
        vf.alt_alleles = vec![Allele::from_str("G"), Allele::from_str("T")];
        vf.allele_string = "A/G/T".into();
        let mut tv_g = base_tv.clone();
        tv_g.allele_annotations[0].allele = Allele::from_str("G");
        tv_g.allele_annotations[0].supplementary = vec![(
            "clinvar".into(),
            r#"{"significance":["Pathogenic"],"reviewStatus":"","phenotypes":[],"variantClass":"SNV","soAccession":""}"#.into(),
        )];
        let mut tv_t = base_tv.clone();
        tv_t.allele_annotations[0].allele = Allele::from_str("T");
        tv_t.allele_annotations[0].supplementary = vec![(
            "clinvar".into(),
            r#"{"significance":["Benign"],"reviewStatus":"","phenotypes":[],"variantClass":"SNV","soAccession":""}"#.into(),
        )];
        vf.transcript_variations = vec![tv_g, tv_t];

        let specs = LoadedSupplementarySpecs::new(&["clinvar".to_string()], &[]);
        let lines = format_tab_line(&vf, &specs, false);
        assert_eq!(lines.len(), 2);
        assert!(
            lines[0]
                .split('\t')
                .next_back()
                .unwrap()
                .starts_with("G|Pathogenic|"),
            "first row should carry the ALT-G clinvar value: {}",
            lines[0]
        );
        assert!(
            lines[1]
                .split('\t')
                .next_back()
                .unwrap()
                .starts_with("T|Benign|"),
            "second row should carry the ALT-T clinvar value: {}",
            lines[1]
        );
    }

    #[test]
    fn tab_line_gene_projection_works_when_transcript_has_no_symbol() {
        // Regression for the reviewer's Major #3: when a transcript has no
        // gene_symbol, the per-tv gene filter previously fell back to
        // gene_id which never matched any GeneAnnotation. Now we drop the
        // symbol filter in that case so tab output carries the same value
        // as VCF.
        let mut vf = projection_test_variant();
        vf.transcript_variations[0].gene_symbol = None;
        let specs = LoadedSupplementarySpecs::new(&[], &["omim".to_string()]);
        let lines = format_tab_line(&vf, &specs, false);
        let cols: Vec<&str> = lines[0].split('\t').collect();
        let fv_omim = cols.last().unwrap();
        assert!(
            fv_omim.starts_with("GENE1|113705|"),
            "FV_OMIM cell should still carry the gene annotation: {}",
            fv_omim
        );
    }

    #[test]
    fn gene_projections_flatten_json_array_gene_payloads() {
        // GeneIndex::annotate_gene wraps multi-record genes in a JSON
        // array (`[{...},{...}]`). The variant- and per-tv-level helpers
        // must flatten that into one pipe entry per record instead of
        // emitting just `GENE1` with no fields (or dropping the value
        // entirely for clinvar_protein). Regression test for Copilot
        // review of PR #31.
        let mut vf = projection_test_variant();
        vf.gene_annotations = vec![
            GeneAnnotation {
                gene_symbol: "GENE1".into(),
                json_key: "omim".into(),
                json_string: r#"[{"mimNumber":113705,"phenotypes":["Breast cancer"]},{"mimNumber":167000,"phenotypes":["Ovarian cancer"]}]"#.into(),
            },
            GeneAnnotation {
                gene_symbol: "GENE1".into(),
                json_key: "clinvar_protein".into(),
                json_string: r#"[{"proteinVariants":[{"pos":175,"refAa":"R","altAa":"H","sig":"Pathogenic"}]},{"proteinVariants":[{"pos":248,"refAa":"R","altAa":"W","sig":"Pathogenic"}]}]"#.into(),
            },
        ];

        // VCF projection — should carry both OMIM entries comma-joined.
        let info: std::collections::HashMap<String, String> =
            format_supplementary_vcf_info(&vf, &no_custom_sources())
                .into_iter()
                .collect();
        let omim = info.get("FV_OMIM").expect("FV_OMIM should be present");
        assert!(
            omim.contains("GENE1|113705|Breast%20cancer"),
            "FV_OMIM should carry first record: {omim}"
        );
        assert!(
            omim.contains("GENE1|167000|Ovarian%20cancer"),
            "FV_OMIM should carry second record: {omim}"
        );
        let cvp = info
            .get("FV_CLINVAR_PROTEIN")
            .expect("FV_CLINVAR_PROTEIN should be present even for array payload");
        assert!(
            cvp.contains("175:R>H:Pathogenic"),
            "FV_CLINVAR_PROTEIN should carry first proteinVariants entry: {cvp}"
        );
        assert!(
            cvp.contains("248:R>W:Pathogenic"),
            "FV_CLINVAR_PROTEIN should carry second proteinVariants entry: {cvp}"
        );

        // Tab projection — same pipe values land in the per-row cells.
        let specs = LoadedSupplementarySpecs::new(
            &[],
            &["omim".to_string(), "clinvar_protein".to_string()],
        );
        let lines = format_tab_line(&vf, &specs, false);
        let row = &lines[0];
        let cols: Vec<&str> = row.split('\t').collect();
        let fv_omim = cols[17];
        let fv_cvp = cols[18];
        assert!(
            fv_omim.contains("113705") && fv_omim.contains("167000"),
            "tab FV_OMIM cell should carry both records: {fv_omim}"
        );
        assert!(
            fv_cvp.contains("175:R>H") && fv_cvp.contains("248:R>W"),
            "tab FV_CLINVAR_PROTEIN cell should carry both proteinVariants: {fv_cvp}"
        );
    }

    #[test]
    fn tab_line_gracefully_handles_malformed_supplementary_json() {
        // Garbage inside `aa.supplementary` must not panic and must not
        // bleed into the tab output as raw bytes.
        let mut vf = projection_test_variant();
        for tv in &mut vf.transcript_variations {
            for aa in &mut tv.allele_annotations {
                aa.supplementary = vec![
                    ("clinvar".into(), "{not valid json".into()),
                    ("dbnsfp".into(), "".into()),
                ];
            }
        }
        vf.gene_annotations.clear();
        let specs =
            LoadedSupplementarySpecs::new(&["clinvar".to_string(), "dbnsfp".to_string()], &[]);
        let lines = format_tab_line(&vf, &specs, false);
        let row = &lines[0];
        assert!(!row.contains('{'), "malformed JSON must not leak: {row}");
        // Whatever the renderer does with the broken payload, the cell is
        // a single tab-delimited field — not a multi-tab eruption.
        let cols: Vec<&str> = row.split('\t').collect();
        assert_eq!(cols.len(), 17 + 2);
    }

    #[test]
    fn supplementary_annotations_doc_lists_every_projection_spec() {
        // The user-facing `docs/SUPPLEMENTARY_ANNOTATIONS.md` is the single
        // place we point users at for the FV_* schema. If a new
        // VCF_PROJECTION_SPECS entry is added without a matching doc entry,
        // downstream parsers won't know its pipe format. This test fails CI
        // in that case so the doc cannot drift from the code.
        let doc = include_str!("../../../docs/SUPPLEMENTARY_ANNOTATIONS.md");
        for spec in VCF_PROJECTION_SPECS {
            assert!(
                doc.contains(spec.info_id),
                "docs/SUPPLEMENTARY_ANNOTATIONS.md is missing INFO ID {}",
                spec.info_id
            );
            // Each spec's Description ends with "Format: <pipe layout>".
            // Verify the pipe layout itself appears verbatim in the doc so
            // grep-style downstream consumers can pull the schema from
            // either source. (The "fastVEP X annotations. Format: " prose
            // intentionally does not need to appear verbatim.)
            let format = spec
                .description
                .split_once("Format: ")
                .map(|(_, layout)| layout)
                .unwrap_or(spec.description);
            assert!(
                doc.contains(format),
                "docs/SUPPLEMENTARY_ANNOTATIONS.md is missing pipe format for {}: expected `{}`",
                spec.info_id,
                format
            );
        }
        let spliceai_format = "ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL";
        assert!(
            doc.contains(spliceai_format),
            "docs/SUPPLEMENTARY_ANNOTATIONS.md must document the SpliceAI pipe format"
        );
    }

    #[test]
    fn vcf_projection_uses_uploaded_indel_alt_allele() {
        let mut vf = projection_test_variant();
        vf.position.start = 26001;
        vf.position.end = 26001;
        vf.ref_allele = Allele::from_str("A");
        vf.alt_alleles = vec![Allele::Deletion];
        vf.allele_string = "A/-".into();
        vf.vcf_fields = Some(crate::variant::VcfFields {
            chrom: "1".into(),
            pos: 26000,
            id: ".".into(),
            ref_allele: "GA".into(),
            alt: "G".into(),
            qual: ".".into(),
            filter: ".".into(),
            info: ".".into(),
            rest: Vec::new(),
        });
        let aa = &mut vf.transcript_variations[0].allele_annotations[0];
        aa.allele = Allele::Deletion;
        aa.supplementary = vec![(
            "spliceAI".into(),
            r#"{"gene":"GENE1","dsAg":0.1,"dsAl":0.0,"dsDg":0.0,"dsDl":0.0,"dpAg":4,"dpAl":7,"dpDg":27,"dpDl":17}"#.into(),
        )];

        let info = format_supplementary_vcf_info(&vf, &no_custom_sources())
            .into_iter()
            .map(|(id, value)| format!("{id}={value}"))
            .collect::<Vec<_>>()
            .join(";");
        assert!(
            info.contains("SpliceAI=G|GENE1|0.10|0.00|0.00|0.00|4|7|27|17"),
            "{info}"
        );
        assert!(!info.contains("SpliceAI=-|"), "{info}");
    }

    #[test]
    fn tab_options_explicit_ref_adds_ref_column() {
        let vf = projection_test_variant();
        let specs = LoadedSupplementarySpecs::new(&[], &[]);
        let opts = TabOptions {
            explicit_ref: true,
            ..TabOptions::default()
        };
        let lines = format_tab_line_with(&vf, &specs, false, opts);
        assert_eq!(lines.len(), 1);
        let cols: Vec<&str> = lines[0].split('\t').collect();
        // Base 17 columns + 1 REF column inserted after Allele (index 2).
        assert_eq!(cols.len(), 18, "expected 18 cols, got: {:?}", cols);
        assert_eq!(cols[2], "G", "Allele column should still hold ALT");
        assert_eq!(cols[3], "A", "REF column should follow Allele");
    }

    #[test]
    fn tab_options_qc_class_adds_trailing_column() {
        let vf = projection_test_variant();
        let specs = LoadedSupplementarySpecs::new(&[], &[]);
        let opts = TabOptions {
            qc_class: Some("HIGH_QC"),
            ..TabOptions::default()
        };
        let lines = format_tab_line_with(&vf, &specs, false, opts);
        let cols: Vec<&str> = lines[0].split('\t').collect();
        assert_eq!(cols.len(), 18);
        assert_eq!(cols.last().copied(), Some("HIGH_QC"));
    }

    #[test]
    fn tab_options_gene_set_drops_non_matching_transcripts() {
        let vf = projection_test_variant();
        let specs = LoadedSupplementarySpecs::new(&[], &[]);

        let included = crate::geneset::GeneSet::from_iter(["GENE1"]);
        let opts_match = TabOptions {
            gene_set: Some(&included),
            ..TabOptions::default()
        };
        let lines = format_tab_line_with(&vf, &specs, false, opts_match);
        assert_eq!(lines.len(), 1, "GENE1 in panel → row kept");

        let excluded = crate::geneset::GeneSet::from_iter(["OTHER"]);
        let opts_no_match = TabOptions {
            gene_set: Some(&excluded),
            ..TabOptions::default()
        };
        let lines = format_tab_line_with(&vf, &specs, false, opts_no_match);
        assert!(lines.is_empty(), "GENE1 not in panel → row dropped");
    }

    #[test]
    fn tab_options_combined_ref_and_qc_class_preserve_order() {
        let vf = projection_test_variant();
        let specs = LoadedSupplementarySpecs::new(&[], &[]);
        let opts = TabOptions {
            explicit_ref: true,
            qc_class: Some("HIGH_QC"),
            ..TabOptions::default()
        };
        let lines = format_tab_line_with(&vf, &specs, false, opts);
        let cols: Vec<&str> = lines[0].split('\t').collect();
        // Base 17 + REF + QC = 19. REF at index 3, QC last.
        assert_eq!(cols.len(), 19);
        assert_eq!(cols[3], "A");
        assert_eq!(cols.last().copied(), Some("HIGH_QC"));
    }

    #[test]
    fn tab_line_default_unchanged_when_options_empty() {
        // Regression guard: the no-option call path must match the original
        // 17-column shape exactly.
        let vf = projection_test_variant();
        let specs = LoadedSupplementarySpecs::new(&[], &[]);
        let baseline = format_tab_line(&vf, &specs, false);
        let with_default = format_tab_line_with(&vf, &specs, false, TabOptions::default());
        assert_eq!(baseline, with_default);
        assert_eq!(baseline[0].split('\t').count(), 17);
    }

    // -- the cell rules VEP's tab and CSQ writers follow (issue #139) ----------

    /// What `annotate_intergenic` pushes when no transcript overlaps.
    fn intergenic_variant() -> VariationFeature {
        let mut vf = projection_test_variant();
        let tv = &mut vf.transcript_variations[0];
        tv.transcript_id = Arc::from("-");
        tv.gene_id = Arc::from("-");
        tv.gene_symbol = None;
        tv.biotype = Arc::from("-");
        let aa = &mut tv.allele_annotations[0];
        aa.consequences = vec![Consequence::IntergenicVariant];
        aa.impact = Impact::Modifier;
        aa.supplementary.clear();
        vf
    }

    fn tab_cells(vf: &VariationFeature) -> Vec<String> {
        let specs = LoadedSupplementarySpecs::new(&[], &[]);
        let rows = format_tab_line(vf, &specs, false);
        assert_eq!(rows.len(), 1, "one allele on one row: {rows:?}");
        rows[0].split('\t').map(String::from).collect()
    }

    fn csq_by_name(vf: &VariationFeature) -> std::collections::HashMap<&'static str, String> {
        let csq = format_csq(vf, DEFAULT_CSQ_FIELDS);
        let entry = csq.split(',').next().expect("one entry");
        DEFAULT_CSQ_FIELDS
            .iter()
            .copied()
            .zip(entry.split('|').map(String::from))
            .collect()
    }

    /// `Uploaded_variation` for a VCF record without an ID is the line as it was
    /// written, `{CHROM}_{POS}_{REF}/{ALT}`, whatever the annotation trimmed it
    /// to. The expected names are what Ensembl VEP 116.3 printed for these lines,
    /// run in Docker (`crates/fastvep-cli/fixtures/vep_tab_world`).
    #[test]
    fn uploaded_variation_of_a_vcf_record_is_the_line_as_written() {
        let line = |chrom: &str, pos: u64, reference: &str, alt: &str| VcfFields {
            chrom: chrom.into(),
            pos,
            id: ".".into(),
            ref_allele: reference.into(),
            alt: alt.into(),
            qual: ".".into(),
            filter: ".".into(),
            info: ".".into(),
            rest: vec![],
        };
        // (the line's POS, REF, ALT) against the minimal coordinates and
        // alleles the annotation works with
        for (fields, start, end, minimal, want) in [
            (
                line("chr10", 10400, "AACCCT", "A"),
                10401,
                10405,
                "ACCCT/-",
                "chr10_10400_AACCCT/A",
            ),
            (
                line("chr10", 10408, "C", "CT"),
                10409,
                10408,
                "-/T",
                "chr10_10408_C/CT",
            ),
            (
                line("chr10", 10414, "C", "A"),
                10414,
                10414,
                "C/A",
                "chr10_10414_C/A",
            ),
            (
                line("chr10", 10414, "C", "A,G"),
                10414,
                10414,
                "C/A/G",
                "chr10_10414_C/A/G",
            ),
            (
                line("chr10", 109, "CA", "C,TA"),
                109,
                110,
                "CA/C/TA",
                "chr10_109_CA/C/TA",
            ),
            (line("chrM", 100, "A", "G"), 100, 100, "A/G", "chrM_100_A/G"),
        ] {
            let mut vf = projection_test_variant();
            vf.position = GenomicPosition::new(fields.chrom.as_str(), start, end, Strand::Forward);
            vf.allele_string = minimal.into();
            vf.vcf_fields = Some(fields);
            assert_eq!(uploaded_variation(&vf), want);
            assert_eq!(
                tab_cells(&vf)[0],
                want,
                "and the tab row's first cell is that name"
            );
        }
    }

    /// An ID wins, the first of several; `.` is no ID.
    #[test]
    fn uploaded_variation_prefers_the_first_id() {
        let mut vf = projection_test_variant();
        vf.variation_name = Some("rs710322;rs187914".into());
        assert_eq!(uploaded_variation(&vf), "rs710322");
        vf.variation_name = Some(".".into());
        assert_eq!(uploaded_variation(&vf), "1_25000_A/G");
    }

    /// Input that is not a VCF line has only the annotation to name it: the
    /// start alone, not a range, and the alleles it has.
    #[test]
    fn uploaded_variation_without_a_vcf_line_is_named_from_the_annotation() {
        let mut vf = projection_test_variant();
        assert!(vf.vcf_fields.is_none());
        assert_eq!(uploaded_variation(&vf), "1_25000_A/G");

        vf.position = GenomicPosition::new("chr10", 10401, 10405, Strand::Forward);
        vf.allele_string = "ACCCT/-".into();
        assert_eq!(uploaded_variation(&vf), "chr10_10401_ACCCT/-");
    }

    /// `pep_allele_string`: one letter when the protein does not change. The CSQ
    /// writer had it; the tab and JSON writers wrote `L/L`. VEP 116.3 prints `L`
    /// in all three.
    #[test]
    fn identical_amino_acids_are_written_once() {
        assert_eq!(amino_acids_text("L", "L"), "L");
        assert_eq!(amino_acids_text("L", "P"), "L/P");
        assert_eq!(amino_acids_text("*", "*"), "*");

        let mut vf = projection_test_variant();
        vf.transcript_variations[0].allele_annotations[0].amino_acids =
            Some(("T".into(), "T".into()));
        assert_eq!(tab_cells(&vf)[10], "T");
        vf.transcript_variations[0].allele_annotations[0].amino_acids =
            Some(("T".into(), "M".into()));
        assert_eq!(tab_cells(&vf)[10], "T/M");
    }

    /// VEP's `format_coords` puts the lower coordinate first, so an insertion
    /// (end below start) is `chr10:11950-11951`. It was written `11951-11950`.
    #[test]
    fn location_is_written_low_to_high_as_vep_does() {
        let mut vf = projection_test_variant();
        for (start, end, want) in [
            (25000, 25000, "1:25000"),
            (25000, 25004, "1:25000-25004"),
            (25001, 25000, "1:25000-25001"), // an insertion
        ] {
            vf.position = GenomicPosition::new("1", start, end, Strand::Forward);
            assert_eq!(tab_cells(&vf)[1], want, "{start}..{end}");
        }
    }

    #[test]
    fn uploaded_variation_is_the_first_id_when_the_record_has_one() {
        let mut vf = projection_test_variant();
        for (id, want) in [
            (Some("rs123"), "rs123"),
            // VEP keeps `$ids->[0]` of the semicolon-separated column
            (Some("rs1;rs2"), "rs1"),
            // no usable ID: the generated name
            (Some("."), "1_25000_A/G"),
            (Some(""), "1_25000_A/G"),
            (None, "1_25000_A/G"),
        ] {
            vf.variation_name = id.map(String::from);
            assert_eq!(uploaded_variation(&vf), want, "ID {id:?}");
        }
    }

    /// VEP's tab writer prints `-` for every value it does not have
    /// (`Tab::output_hash_to_line`). The position cells printed an empty string.
    #[test]
    fn an_absent_value_is_a_dash_in_every_tab_cell() {
        let cells = tab_cells(&projection_test_variant());
        assert_eq!(cells.len(), 17);
        assert!(
            cells.iter().all(|c| !c.is_empty()),
            "an empty cell where VEP writes `-`: {cells:?}"
        );
        // cDNA, CDS and protein position
        assert_eq!(cells[7..=9].to_vec(), vec!["-", "-", "-"]);

        let mut vf = projection_test_variant();
        let location = &mut vf.transcript_variations[0].allele_annotations[0].location;
        location.cdna = TranscriptSpan::new(Some(100), Some(100));
        location.cds = TranscriptSpan::new(Some(50), Some(52));
        location.protein = TranscriptSpan::new(Some(17), Some(17));
        assert_eq!(tab_cells(&vf)[7..=9].to_vec(), vec!["100", "50-52", "17"]);

        // an end the transcript has no coordinate for is a `?`, as in VEP
        let location = &mut vf.transcript_variations[0].allele_annotations[0].location;
        location.cdna = TranscriptSpan::new(None, Some(7));
        location.cds = TranscriptSpan::new(Some(601), None);
        assert_eq!(tab_cells(&vf)[7..=8].to_vec(), vec!["?-7", "601-?"]);
    }

    /// VEP's JSON writer splits the position cell on `-` and leaves out an end that is
    /// `?` (`OutputFactory::JSON`), so a change that begins outside an exon has a
    /// `cdna_end` and no `cdna_start`, and an insertion's pair is written low to high.
    #[test]
    fn json_leaves_out_the_end_of_a_span_that_has_no_coordinate() {
        let mut vf = projection_test_variant();
        let location = &mut vf.transcript_variations[0].allele_annotations[0].location;
        location.cdna = TranscriptSpan::new(None, Some(7));
        location.cds = TranscriptSpan::new(Some(601), None);
        location.protein = TranscriptSpan::new(Some(68), Some(67));
        location.exon = Some(Ordinals {
            first: 2,
            last: 3,
            total: 10,
        });
        let json = format_json(&vf, false);
        let tc = &json["transcript_consequences"][0];
        assert!(tc.get("cdna_start").is_none(), "{tc}");
        assert_eq!(tc["cdna_end"], 7);
        assert_eq!(tc["cds_start"], 601);
        assert!(tc.get("cds_end").is_none(), "{tc}");
        assert_eq!(tc["protein_start"], 67);
        assert_eq!(tc["protein_end"], 68);
        assert_eq!(tc["exon"], "2-3/10");
        assert!(tc.get("intron").is_none(), "{tc}");
    }

    #[test]
    fn a_transcript_row_keeps_its_feature_cells() {
        let cells = tab_cells(&projection_test_variant());
        assert_eq!(cells[3], "GENE1");
        assert_eq!(cells[4], "TX1");
        assert_eq!(cells[5], "Transcript");
        assert_eq!(cells[15], "1");
    }

    /// An intergenic row has no feature, so VEP never sets `Feature_type` or
    /// `STRAND` for it. fastVEP printed `Transcript` and `1`, which also made
    /// `Feature_type is Transcript` select intergenic rows.
    #[test]
    fn an_intergenic_row_has_no_feature_cells_in_tab() {
        let cells = tab_cells(&intergenic_variant());
        // Gene, Feature, Feature_type
        assert_eq!(cells[3..=5].to_vec(), vec!["-", "-", "-"]);
        assert_eq!(cells[6], "intergenic_variant");
        assert_eq!(cells[13], "MODIFIER");
        // DISTANCE, STRAND, FLAGS
        assert_eq!(cells[14..=16].to_vec(), vec!["-", "-", "-"]);
        assert!(cells.iter().all(|c| !c.is_empty()), "{cells:?}");
    }

    /// The row written when a variant has no transcript rows at all printed `-`
    /// for IMPACT, where an intergenic variant is MODIFIER.
    #[test]
    fn the_row_for_a_variant_with_no_transcript_rows_is_modifier() {
        let mut vf = projection_test_variant();
        vf.transcript_variations.clear();
        let cells = tab_cells(&vf);
        assert_eq!(cells.len(), 17);
        assert_eq!(cells[6], "intergenic_variant");
        assert_eq!(cells[13], "MODIFIER");
        assert_eq!(cells[3..=5].to_vec(), vec!["-", "-", "-"]);
        assert_eq!(cells[14..=16].to_vec(), vec!["-", "-", "-"]);
        assert!(cells.iter().all(|c| !c.is_empty()), "{cells:?}");
    }

    #[test]
    fn an_intergenic_csq_entry_leaves_the_feature_fields_empty() {
        let f = csq_by_name(&intergenic_variant());
        for name in [
            "SYMBOL",
            "Gene",
            "Feature_type",
            "Feature",
            "BIOTYPE",
            "DISTANCE",
            "STRAND",
        ] {
            assert_eq!(f[name], "", "{name} of an intergenic entry");
        }
        assert_eq!(f["Consequence"], "intergenic_variant");
        assert_eq!(f["IMPACT"], "MODIFIER");

        let g = csq_by_name(&projection_test_variant());
        assert_eq!(g["Gene"], "GENE1");
        assert_eq!(g["Feature_type"], "Transcript");
        assert_eq!(g["Feature"], "TX1");
        assert_eq!(g["BIOTYPE"], "protein_coding");
        assert_eq!(g["STRAND"], "1");
    }

    #[test]
    fn only_the_placeholder_is_featureless() {
        assert!(intergenic_variant().transcript_variations[0].is_featureless());
        assert!(!projection_test_variant().transcript_variations[0].is_featureless());
    }

    // -- VEP's optional columns (issue #139) ----------------------------------

    fn names(columns: &[OptionalColumn]) -> Vec<&'static str> {
        columns.iter().map(|c| c.name()).collect()
    }

    /// The order is VEP's `FLAG_FIELDS` table (Constants.pm), not the order the
    /// flags were typed in, and a flag adds exactly its own columns.
    #[test]
    fn the_flags_give_vep_s_columns_in_vep_s_order() {
        assert!(ColumnFlags::default().columns(false).is_empty());
        assert!(!ColumnFlags::default().any());

        let all = ColumnFlags {
            show_ref_allele: true,
            variant_class: true,
            symbol: true,
            biotype: true,
            canonical: true,
            numbers: true,
            hgvsg: true,
        };
        assert_eq!(
            names(&all.columns(true)),
            [
                "REF_ALLELE",
                "VARIANT_CLASS",
                "SYMBOL",
                "SYMBOL_SOURCE",
                "HGNC_ID",
                "BIOTYPE",
                "CANONICAL",
                "EXON",
                "INTRON",
                "HGVSc",
                "HGVSp",
                "HGVS_OFFSET",
                "HGVSg",
            ]
        );
        // --hgvs without --hgvsg has no HGVSg, and --hgvsg without --hgvs no HGVSc
        let only_g = ColumnFlags {
            hgvsg: true,
            ..ColumnFlags::default()
        };
        assert_eq!(names(&only_g.columns(false)), ["HGVSg"]);
        assert_eq!(
            names(&ColumnFlags::default().columns(true)),
            ["HGVSc", "HGVSp", "HGVS_OFFSET"]
        );
    }

    /// VEP's `--everything` (Config.pm OPTION_SETS) has these, among the columns
    /// fastVEP writes, and neither `--hgvsg` nor `--show_ref_allele`.
    #[test]
    fn everything_is_vep_s_set_among_the_columns_fastvep_writes() {
        let e = ColumnFlags::everything();
        assert_eq!(
            names(&e.columns(true)),
            [
                "VARIANT_CLASS",
                "SYMBOL",
                "SYMBOL_SOURCE",
                "HGNC_ID",
                "BIOTYPE",
                "CANONICAL",
                "EXON",
                "INTRON",
                "HGVSc",
                "HGVSp",
                "HGVS_OFFSET",
            ]
        );
        // union keeps what either has
        let asked = ColumnFlags {
            hgvsg: true,
            show_ref_allele: true,
            ..ColumnFlags::default()
        };
        let both = e.union(asked);
        assert!(both.hgvsg && both.show_ref_allele && both.symbol && both.variant_class);
    }

    #[test]
    fn the_csq_field_list_only_grows_on_request() {
        // nothing asked: the 49 fields every reader already parses
        assert_eq!(
            csq_fields(false, ColumnFlags::default()),
            DEFAULT_CSQ_FIELDS.to_vec()
        );
        let both = ColumnFlags {
            variant_class: true,
            hgvsg: true,
            ..ColumnFlags::default()
        };
        let fields = csq_fields(true, both);
        let at = |name: &str| fields.iter().position(|f| *f == name).unwrap();
        assert_eq!(fields.len(), DEFAULT_CSQ_FIELDS.len() + 3);
        // PICK then VARIANT_CLASS after FLAGS, HGVSg after HGVS_OFFSET: VEP's order
        assert_eq!(at("PICK"), at("FLAGS") + 1);
        assert_eq!(at("VARIANT_CLASS"), at("PICK") + 1);
        assert_eq!(at("HGVSg"), at("HGVS_OFFSET") + 1);
        // and the pick-only list is what it was
        assert_eq!(
            csq_fields_with_pick(),
            csq_fields(true, ColumnFlags::default())
        );
    }

    fn cells_with(vf: &VariationFeature, flags: ColumnFlags, hgvs: bool) -> Vec<String> {
        let specs = LoadedSupplementarySpecs::new(&[], &[]);
        let optional = flags.columns(hgvs);
        let opts = TabOptions {
            optional: &optional,
            ..TabOptions::default()
        };
        let rows = format_tab_line_with(vf, &specs, false, opts);
        assert_eq!(rows.len(), 1, "{rows:?}");
        rows[0].split('\t').map(String::from).collect()
    }

    /// The optional cells follow FLAGS, except `REF_ALLELE`, which VEP puts before
    /// `IMPACT`, and come before everything fastVEP-specific.
    #[test]
    fn optional_cells_sit_after_flags_and_before_the_extras() {
        let mut vf = projection_test_variant();
        {
            let tv = &mut vf.transcript_variations[0];
            tv.canonical = true;
            tv.symbol_source = Some("HGNC".into());
            tv.hgnc_id = Some("HGNC:5".into());
            let aa = &mut tv.allele_annotations[0];
            aa.location.exon = Some(Ordinals {
                first: 3,
                last: 3,
                total: 10,
            });
            aa.hgvsc = Some("TX1:c.100A>G".into());
            aa.hgvsp = Some("ENSP1:p.Ile34Val".into());
            aa.hgvs_offset = Some(2);
            aa.hgvsg = Some("1:g.25000A>G".into());
        }
        let flags = ColumnFlags {
            show_ref_allele: true,
            variant_class: true,
            symbol: true,
            biotype: true,
            canonical: true,
            numbers: true,
            hgvsg: true,
        };
        let cells = cells_with(&vf, flags, true);
        assert_eq!(cells.len(), 17 + 13);
        assert_eq!(
            cells[13], "A",
            "REF_ALLELE is the first column after the 13"
        );
        assert_eq!(
            cells[18..].to_vec(),
            vec![
                "SNV",
                "GENE1",
                "HGNC",
                "HGNC:5",
                "protein_coding",
                "YES",
                "3/10",
                "-",
                "TX1:c.100A>G",
                "ENSP1:p.Ile34Val",
                "2",
                "1:g.25000A>G",
            ]
        );
    }

    #[test]
    fn a_row_with_none_of_the_optional_values_has_dashes() {
        let vf = projection_test_variant();
        let flags = ColumnFlags {
            symbol: true,
            canonical: true,
            numbers: true,
            ..ColumnFlags::default()
        };
        let cells = cells_with(&vf, flags, true);
        // SYMBOL is GENE1 (the test gene has a symbol); the rest has nothing
        assert_eq!(
            cells[17..].to_vec(),
            vec!["GENE1", "-", "-", "-", "-", "-", "-", "-", "-"]
        );
    }

    /// An intergenic row has the variant's own values and no feature's.
    #[test]
    fn an_intergenic_row_keeps_the_variant_level_columns_only() {
        let mut vf = intergenic_variant();
        vf.transcript_variations[0].allele_annotations[0].hgvsg = Some("1:g.25000A>G".into());
        let flags = ColumnFlags {
            show_ref_allele: true,
            variant_class: true,
            symbol: true,
            biotype: true,
            canonical: true,
            numbers: true,
            hgvsg: true,
        };
        let cells = cells_with(&vf, flags, false);
        assert_eq!(cells[13], "A", "REF_ALLELE");
        assert_eq!(
            cells[18..].to_vec(),
            // VARIANT_CLASS, SYMBOL x3, BIOTYPE, CANONICAL, EXON, INTRON, HGVSg
            vec!["SNV", "-", "-", "-", "-", "-", "-", "-", "1:g.25000A>G"]
        );
    }

    /// The row written for a variant with no transcript row at all keeps them too.
    #[test]
    fn the_fallback_row_keeps_the_variant_level_columns_only() {
        let mut vf = projection_test_variant();
        vf.transcript_variations.clear();
        let flags = ColumnFlags {
            show_ref_allele: true,
            variant_class: true,
            biotype: true,
            hgvsg: true,
            ..ColumnFlags::default()
        };
        let cells = cells_with(&vf, flags, false);
        assert_eq!(cells.len(), 17 + 4);
        assert_eq!(cells[13], "A", "REF_ALLELE");
        assert_eq!(cells[18..].to_vec(), vec!["SNV", "-", "-"]);
    }

    /// A structural variant's class comes from the parser in VEP, not from its
    /// allele string, so there is nothing to write for it here.
    #[test]
    fn a_structural_variant_has_no_variant_class() {
        let mut vf = projection_test_variant();
        vf.variant_type = VariantType::CopyNumberLoss;
        let flags = ColumnFlags {
            variant_class: true,
            ..ColumnFlags::default()
        };
        assert_eq!(cells_with(&vf, flags, false)[17], "-");
        let f = csq_by_name_with(&vf, flags);
        assert_eq!(f["VARIANT_CLASS"], "");
    }

    fn csq_by_name_with(
        vf: &VariationFeature,
        flags: ColumnFlags,
    ) -> std::collections::HashMap<&'static str, String> {
        let fields = csq_fields(false, flags);
        let csq = format_csq(vf, &fields);
        let entry = csq.split(',').next().expect("one entry");
        fields
            .iter()
            .copied()
            .zip(entry.split('|').map(String::from))
            .collect()
    }

    #[test]
    fn csq_writes_variant_class_and_hgvsg_when_they_are_in_the_field_list() {
        let mut vf = projection_test_variant();
        vf.transcript_variations[0].allele_annotations[0].hgvsg = Some("1:g.25000A>G".into());
        let flags = ColumnFlags {
            variant_class: true,
            hgvsg: true,
            ..ColumnFlags::default()
        };
        let f = csq_by_name_with(&vf, flags);
        assert_eq!(f["VARIANT_CLASS"], "SNV");
        assert_eq!(f["HGVSg"], "1:g.25000A>G");
    }

    /// FLAGS is the transcript's own flags, as in VEP; the word `canonical` it
    /// used to carry is CANONICAL's job.
    #[test]
    fn the_flags_cell_holds_transcript_flags_and_not_the_canonical_marker() {
        let mut vf = projection_test_variant();
        vf.transcript_variations[0].canonical = true;
        assert_eq!(tab_cells(&vf)[16], "-");

        vf.transcript_variations[0].flags = vec!["cds_start_NF".into(), "cds_end_NF".into()];
        assert_eq!(tab_cells(&vf)[16], "cds_start_NF,cds_end_NF");
    }

    /// The header and the row writer are two pieces of code that have to agree on
    /// the base columns; the header is written from this list.
    #[test]
    fn the_base_columns_are_the_ones_the_row_writer_writes() {
        assert_eq!(TAB_BASE_COLUMNS.len(), 17);
        assert_eq!(
            tab_cells(&projection_test_variant()).len(),
            TAB_BASE_COLUMNS.len()
        );
        // the first 13 are VEP's DEFAULT_OUTPUT_COLS (Constants.pm)
        assert_eq!(
            TAB_BASE_COLUMNS[..13].to_vec(),
            [
                "Uploaded_variation",
                "Location",
                "Allele",
                "Gene",
                "Feature",
                "Feature_type",
                "Consequence",
                "cDNA_position",
                "CDS_position",
                "Protein_position",
                "Amino_acids",
                "Codons",
                "Existing_variation",
            ]
        );
    }

    /// What a real VEP 116.3 run writes after the 13 default columns when every
    /// flag fastVEP has is given: `FLAG_FIELDS` order, `REF_ALLELE` ahead of
    /// `IMPACT`, `PICK` straight after `FLAGS`.
    #[test]
    fn the_header_follows_veps_flag_fields_order() {
        let optional = ColumnFlags {
            show_ref_allele: true,
            variant_class: true,
            symbol: true,
            biotype: true,
            canonical: true,
            numbers: true,
            hgvsg: true,
        }
        .columns(true);
        let header = tab_header(false, false, &optional, true, &[], false);
        let names: Vec<&str> = header.trim_start_matches('#').split('\t').collect();
        assert_eq!(names[..13], TAB_BASE_COLUMNS[..13]);
        assert_eq!(
            names[13..].to_vec(),
            [
                "REF_ALLELE",
                "IMPACT",
                "DISTANCE",
                "STRAND",
                "FLAGS",
                "PICK",
                "VARIANT_CLASS",
                "SYMBOL",
                "SYMBOL_SOURCE",
                "HGNC_ID",
                "BIOTYPE",
                "CANONICAL",
                "EXON",
                "INTRON",
                "HGVSc",
                "HGVSp",
                "HGVS_OFFSET",
                "HGVSg",
            ]
        );
    }

    /// The header and the rows are written by two pieces of code. For every
    /// combination of options, each column's header names the cell under it.
    #[test]
    fn every_column_of_the_header_has_its_own_cell_under_it() {
        let mut vf = projection_test_variant();
        {
            let tv = &mut vf.transcript_variations[0];
            tv.canonical = true;
            tv.symbol_source = Some("HGNC".into());
            tv.hgnc_id = Some("HGNC:5".into());
            let aa = &mut tv.allele_annotations[0];
            aa.pick = true;
            aa.location.exon = Some(Ordinals {
                first: 3,
                last: 3,
                total: 10,
            });
            aa.hgvsc = Some("TX1:c.100A>G".into());
            aa.hgvsp = Some("ENSP1:p.Ile34Val".into());
            aa.hgvs_offset = Some(2);
            aa.hgvsg = Some("1:g.25000A>G".into());
        }
        let expected = |name: &str| -> Option<&'static str> {
            Some(match name {
                "Allele" => "G",
                "REF" | "REF_ALLELE" => "A",
                "IMPACT" => "MODERATE",
                "FLAGS" => "-",
                "PICK" => "1",
                "QC_CLASS" => "pass",
                "VARIANT_CLASS" => "SNV",
                "SYMBOL" => "GENE1",
                "SYMBOL_SOURCE" => "HGNC",
                "HGNC_ID" => "HGNC:5",
                "BIOTYPE" => "protein_coding",
                "CANONICAL" => "YES",
                "EXON" => "3/10",
                "HGVSc" => "TX1:c.100A>G",
                "HGVSp" => "ENSP1:p.Ile34Val",
                "HGVS_OFFSET" => "2",
                "HGVSg" => "1:g.25000A>G",
                _ => return None,
            })
        };
        let specs = LoadedSupplementarySpecs::new(&[], &[]);
        let mut combinations = 0;
        for bits in 0u32..(1 << 10) {
            let on = |i: u32| bits & (1 << i) != 0;
            let flags = ColumnFlags {
                show_ref_allele: on(0),
                variant_class: on(1),
                symbol: on(2),
                biotype: on(3),
                canonical: on(4),
                numbers: on(5),
                hgvsg: on(6),
            };
            let (hgvs, pick, explicit_ref) = (on(7), on(8), on(9));
            let optional = flags.columns(hgvs);
            for qc in [None, Some("pass")] {
                let opts = TabOptions {
                    explicit_ref,
                    qc_class: qc,
                    pick_column: pick,
                    optional: &optional,
                    ..TabOptions::default()
                };
                let header = tab_header(false, explicit_ref, &optional, pick, &[], qc.is_some());
                let names: Vec<&str> = header.trim_start_matches('#').split('\t').collect();
                let rows = format_tab_line_with(&vf, &specs, false, opts);
                assert_eq!(rows.len(), 1);
                let cells: Vec<&str> = rows[0].split('\t').collect();
                assert_eq!(names.len(), cells.len(), "{names:?} against {cells:?}");
                for (name, cell) in names.iter().zip(&cells) {
                    if let Some(want) = expected(name) {
                        assert_eq!(*cell, want, "column {name} of {names:?}");
                    }
                }
                combinations += 1;
            }
        }
        assert_eq!(combinations, 2048);
    }

    /// VEP encodes `=` as `%3D` in the HGVSp of the tab file, as it does in the
    /// VCF: a synonymous change is `p.Thr11%3D` in both.
    #[test]
    fn the_tab_hgvsp_encodes_an_equals_sign_as_vep_does() {
        let mut vf = projection_test_variant();
        vf.transcript_variations[0].allele_annotations[0].hgvsp = Some("ENSP1:p.Thr11=".into());
        let cells = cells_with(&vf, ColumnFlags::default(), true);
        assert_eq!(cells[18], "ENSP1:p.Thr11%3D");
        // and the CSQ, which already did
        let flags = ColumnFlags::default();
        assert_eq!(csq_by_name_with(&vf, flags)["HGVSp"], "ENSP1:p.Thr11%3D");
    }

    /// An insertion's reference allele is `-` in the tab file and empty in the
    /// CSQ, where VEP's VCF writer blanks every `-` but the one in `Allele`.
    #[test]
    fn an_insertions_ref_allele_is_a_dash_in_tab_and_empty_in_csq() {
        let mut vf = projection_test_variant();
        vf.position = GenomicPosition::new("1", 25001, 25000, Strand::Forward);
        vf.allele_string = "-/T".into();
        vf.ref_allele = Allele::Deletion;
        vf.alt_alleles = vec![Allele::from_str("T")];
        vf.transcript_variations[0].allele_annotations[0].allele = Allele::from_str("T");
        let flags = ColumnFlags {
            show_ref_allele: true,
            ..ColumnFlags::default()
        };
        assert_eq!(cells_with(&vf, flags, false)[13], "-");
        assert_eq!(csq_by_name_with(&vf, flags)["REF_ALLELE"], "");
        // a base is a base in both
        let snv = projection_test_variant();
        assert_eq!(cells_with(&snv, flags, false)[13], "A");
        assert_eq!(csq_by_name_with(&snv, flags)["REF_ALLELE"], "A");
    }
}
