use crate::variant::{AlleleAnnotation, TranscriptVariation, VariationFeature};
use fastvep_core::{Allele, Consequence};
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

    let mut first = true;
    for tv in &vf.transcript_variations {
        for aa in &tv.allele_annotations {
            if !first {
                result.push(',');
            }
            first = false;
            format_csq_entry_into(vf, tv, aa, fields, &mut result);
        }
    }

    result
}

/// Write a single CSQ entry directly into the output buffer, avoiding intermediate String allocations.
fn format_csq_entry_into(
    vf: &VariationFeature,
    tv: &TranscriptVariation,
    aa: &AlleleAnnotation,
    fields: &[&str],
    buf: &mut String,
) {
    for (i, field) in fields.iter().enumerate() {
        if i > 0 {
            buf.push('|');
        }
        // Write each field value directly into buf via escape_csq_str, avoiding
        // temporary String allocations for most fields.
        match *field {
            "Allele" => escape_csq_str(&aa.allele.to_string(), buf),
            "Consequence" => {
                for (j, c) in aa.consequences.iter().enumerate() {
                    if j > 0 { buf.push('&'); }
                    buf.push_str(c.so_term());
                }
            }
            "IMPACT" => buf.push_str(aa.impact.as_str()),
            "SYMBOL" => escape_csq_str(tv.gene_symbol.as_deref().unwrap_or_default(), buf),
            "Gene" => escape_csq_str(&tv.gene_id, buf),
            "Feature_type" => buf.push_str("Transcript"),
            "Feature" => escape_csq_str(&tv.transcript_id, buf),
            "BIOTYPE" => escape_csq_str(&tv.biotype, buf),
            "EXON" => {
                if let Some((n, t)) = aa.exon {
                    let _ = write!(buf, "{}/{}", n, t);
                }
            }
            "INTRON" => {
                if let Some((n, t)) = aa.intron {
                    let _ = write!(buf, "{}/{}", n, t);
                }
            }
            "HGVSg" => escape_csq_str(aa.hgvsg.as_deref().unwrap_or_default(), buf),
            "HGVSc" => escape_csq_str(aa.hgvsc.as_deref().unwrap_or_default(), buf),
            "HGVSp" => escape_csq_str(aa.hgvsp.as_deref().unwrap_or_default(), buf),
            "cDNA_position" => write_position_range(aa.cdna_position, buf),
            "CDS_position" => write_position_range(aa.cds_position, buf),
            "Protein_position" => write_position_range(aa.protein_position, buf),
            "Amino_acids" => {
                if let Some((ref r, ref a)) = aa.amino_acids {
                    escape_csq_str(r, buf);
                    if r != a {
                        buf.push('/');
                        escape_csq_str(a, buf);
                    }
                }
            }
            "Codons" => {
                if let Some((ref r, ref a)) = aa.codons {
                    escape_csq_str(r, buf);
                    buf.push('/');
                    escape_csq_str(a, buf);
                }
            }
            "Existing_variation" => {
                for (j, ev) in aa.existing_variation.iter().enumerate() {
                    if j > 0 { buf.push('&'); }
                    escape_csq_str(ev, buf);
                }
            }
            "REF_ALLELE" => escape_csq_str(&vf.ref_allele.to_string(), buf),
            "UPLOADED_ALLELE" => {
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
            "DISTANCE" => {
                if let Some(d) = aa.distance {
                    let _ = write!(buf, "{}", d);
                }
            }
            "STRAND" => { let _ = write!(buf, "{}", tv.strand.as_int()); }
            "FLAGS" => {
                for (j, f) in tv.flags.iter().enumerate() {
                    if j > 0 { buf.push('&'); }
                    buf.push_str(f);
                }
            }
            "CANONICAL" => { if tv.canonical { buf.push_str("YES"); } }
            "SYMBOL_SOURCE" => escape_csq_str(tv.symbol_source.as_deref().unwrap_or_default(), buf),
            "HGNC_ID" => escape_csq_str(tv.hgnc_id.as_deref().unwrap_or_default(), buf),
            "MANE" => {
                if tv.mane_select.is_some() {
                    buf.push_str("MANE_Select");
                } else if tv.mane_plus_clinical.is_some() {
                    buf.push_str("MANE_Plus_Clinical");
                }
            }
            "MANE_SELECT" => escape_csq_str(tv.mane_select.as_deref().unwrap_or_default(), buf),
            "MANE_PLUS_CLINICAL" => escape_csq_str(tv.mane_plus_clinical.as_deref().unwrap_or_default(), buf),
            "TSL" => { if let Some(t) = tv.tsl { let _ = write!(buf, "{}", t); } }
            "APPRIS" => escape_csq_str(tv.appris.as_deref().unwrap_or_default(), buf),
            "CCDS" => escape_csq_str(tv.ccds.as_deref().unwrap_or_default(), buf),
            "GENCODE_PRIMARY" => { if tv.gencode_primary { buf.push_str("YES"); } }
            "ENSP" => escape_csq_str(tv.protein_id.as_deref().unwrap_or_default(), buf),
            "SIFT" => escape_csq_str(aa.sift.as_deref().unwrap_or_default(), buf),
            "PolyPhen" => escape_csq_str(aa.polyphen.as_deref().unwrap_or_default(), buf),
            "AF" => {
                if let Some(f) = vf.existing_variants.iter().find_map(|kv| {
                    kv.frequencies.get("gnomAD")
                        .or_else(|| kv.frequencies.get("gnomADe"))
                        .or_else(|| kv.frequencies.get("minor_allele_freq"))
                }) {
                    let _ = write!(buf, "{}", f);
                }
            }
            "CLIN_SIG" => {
                if let Some(cs) = vf.existing_variants.iter()
                    .find_map(|kv| kv.clinical_significance.as_deref())
                {
                    escape_csq_str(cs, buf);
                }
            }
            "SOMATIC" => { if vf.existing_variants.iter().any(|kv| kv.somatic) { buf.push('1'); } }
            "PHENO" => { if vf.existing_variants.iter().any(|kv| kv.phenotype_or_disease) { buf.push('1'); } }
            "PUBMED" => {
                let mut first_pub = true;
                for kv in &vf.existing_variants {
                    for p in &kv.pubmed {
                        if !first_pub { buf.push('&'); }
                        first_pub = false;
                        buf.push_str(p);
                    }
                }
            }
            "SOURCE" => escape_csq_str(tv.source.as_deref().unwrap_or_default(), buf),
            "HGVS_OFFSET" => {
                if let Some(o) = aa.hgvs_offset {
                    let _ = write!(buf, "{}", o);
                }
            }
            "ACMG" => {
                if let Some(ref acmg) = aa.acmg_classification {
                    if let Some(sh) = acmg.get("shorthand").and_then(|v| v.as_str()) {
                        buf.push_str(sh);
                    }
                }
            }
            "ACMG_CRITERIA" => {
                if let Some(ref acmg) = aa.acmg_classification {
                    if let Some(criteria) = acmg.get("criteria").and_then(|v| v.as_array()) {
                        let mut first = true;
                        for c in criteria {
                            if c.get("met").and_then(|v| v.as_bool()).unwrap_or(false) {
                                if let Some(code) = c.get("code").and_then(|v| v.as_str()) {
                                    if !first { buf.push('&'); }
                                    first = false;
                                    buf.push_str(code);
                                }
                            }
                        }
                    }
                }
            }
            _ => {}
        }
    }
}

/// Escape special characters in CSQ field values, appending to an existing buffer.
fn escape_csq_str(value: &str, buf: &mut String) {
    for c in value.chars() {
        match c {
            ',' | '|' => buf.push('&'),
            ';' => buf.push_str("%3B"),
            '=' => buf.push_str("%3D"),
            _ => buf.push(c),
        }
    }
}

/// Escape special characters in CSQ field values.
#[cfg(test)]
fn escape_csq_value(value: &str) -> String {
    value
        .replace(',', "&")
        .replace(';', "%3B")
        .replace('=', "%3D")
        .replace('|', "&")
        .replace(' ', "_")
}

fn write_position_range(pos: Option<(u64, u64)>, buf: &mut String) {
    match pos {
        Some((start, end)) if start == end => { let _ = write!(buf, "{}", start); }
        Some((start, end)) => { let _ = write!(buf, "{}-{}", start, end); }
        None => {}
    }
}

fn format_position_range(pos: Option<(u64, u64)>) -> String {
    match pos {
        Some((start, end)) if start == end => start.to_string(),
        Some((start, end)) => format!("{}-{}", start, end),
        None => String::new(),
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
const TOPMED_FIELDS: &[(&str, &str)] = &[("ALL_AF", "allAf"), ("ALL_AC", "allAc"), ("ALL_AN", "allAn")];
const MITOMAP_FIELDS: &[(&str, &str)] = &[("DISEASE", "disease"), ("STATUS", "status")];
const SCORE_FIELDS: &[(&str, &str)] = &[("SCORE", "")];
const SCORE_OBJECT_FIELDS: &[(&str, &str)] = &[("SCORE", "score")];
const DBNSFP_FIELDS: &[(&str, &str)] = &[("SIFT", "sift"), ("POLYPHEN", "polyphen")];
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

/// Format fastSA SpliceAI annotations as a VCF-compatible INFO field value.
pub fn format_spliceai_info(vf: &VariationFeature) -> Option<String> {
    format_supplementary_vcf_info(vf)
        .into_iter()
        .find_map(|(id, value)| if id == "SpliceAI" { Some(value) } else { None })
}

/// Return VCF INFO IDs that fastVEP owns for the given loaded sources.
pub fn vcf_owned_info_ids(sa_keys: &[String], gene_keys: &[String]) -> Vec<&'static str> {
    let mut ids = vec!["CSQ"];
    if sa_keys.iter().any(|key| key == "spliceAI") {
        ids.push("SpliceAI");
    }
    for spec in VCF_PROJECTION_SPECS {
        let loaded = match spec.kind {
            VcfProjectionKind::GeneObject | VcfProjectionKind::ClinvarProtein => {
                gene_keys.iter().any(|key| key == spec.json_key)
            }
            VcfProjectionKind::AlleleObject | VcfProjectionKind::AlleleScalar => {
                sa_keys.iter().any(|key| key == spec.json_key)
            }
        };
        if loaded {
            ids.push(spec.info_id);
        }
    }
    ids
}

/// Generate fastVEP-owned VCF INFO header lines for the loaded sources.
pub fn vcf_info_header_lines(
    sa_keys: &[String],
    gene_keys: &[String],
    csq_fields: &[&str],
) -> Vec<String> {
    let mut headers = vec![csq_header_line(csq_fields)];
    if sa_keys.iter().any(|key| key == "spliceAI") {
        headers.push(spliceai_header_line());
    }
    for spec in VCF_PROJECTION_SPECS {
        let loaded = match spec.kind {
            VcfProjectionKind::GeneObject | VcfProjectionKind::ClinvarProtein => {
                gene_keys.iter().any(|key| key == spec.json_key)
            }
            VcfProjectionKind::AlleleObject | VcfProjectionKind::AlleleScalar => {
                sa_keys.iter().any(|key| key == spec.json_key)
            }
        };
        if loaded {
            headers.push(format!(
                "##INFO=<ID={},Number=.,Type=String,Description=\"{}\">",
                spec.info_id, spec.description
            ));
        }
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
    /// Total number of extra tab columns (`spliceai` + number of `true`
    /// entries in `per_spec`), cached at construction so the per-row hot
    /// path doesn't re-scan `per_spec` for the column count.
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
        let column_count = spliceai as usize + per_spec.iter().filter(|b| **b).count();
        Self { spliceai, per_spec, column_count }
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
pub fn tab_supplementary_column_names(specs: &LoadedSupplementarySpecs) -> Vec<&'static str> {
    let mut names = Vec::with_capacity(specs.column_count());
    if specs.spliceai {
        names.push("SpliceAI");
    }
    for spec in specs.loaded_specs() {
        names.push(spec.info_id);
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
    cols
}

/// Format all supplementary VCF INFO projections for an annotated variant.
pub fn format_supplementary_vcf_info(vf: &VariationFeature) -> Vec<(String, String)> {
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

    projected
}

/// Format final VCF INFO by replacing fastVEP-owned fields and appending current projections.
pub fn format_vcf_info_fields(original_info: &str, vf: &VariationFeature, csq: &str) -> String {
    let mut projections = format_supplementary_vcf_info(vf);
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
    let mut values: Vec<String> = Vec::new();
    for tv in &vf.transcript_variations {
        for aa in &tv.allele_annotations {
            let Some(joined) = format_spliceai_for_allele(vf, aa) else {
                continue;
            };
            // The per-aa helper already comma-joins multiple entries for one
            // allele; split them back out so variant-level dedupe stays
            // entry-by-entry rather than across whole allele payloads.
            for value in joined.split(',') {
                let owned = value.to_string();
                if !values.contains(&owned) {
                    values.push(owned);
                }
            }
        }
    }
    if values.is_empty() { None } else { Some(values.join(",")) }
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
    let mut values: Vec<String> = Vec::new();
    for tv in &vf.transcript_variations {
        for aa in &tv.allele_annotations {
            let Some(joined) = format_allele_projection_for_aa(vf, aa, spec) else {
                continue;
            };
            for value in joined.split(',') {
                let owned = value.to_string();
                if !values.contains(&owned) {
                    values.push(owned);
                }
            }
        }
    }
    if values.is_empty() { None } else { Some(values.join(",")) }
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
                    parts.push(obj.get(*json_key).map(json_value_to_vcf).unwrap_or_default());
                }
            }
            let value = parts.join("|");
            if !values.contains(&value) {
                values.push(value);
            }
        }
    }
    if values.is_empty() { None } else { Some(values.join(",")) }
}

fn format_spliceai_for_allele(vf: &VariationFeature, aa: &AlleleAnnotation) -> Option<String> {
    let allele = uploaded_allele_for_annotation(vf, &aa.allele);
    let mut values: Vec<String> = Vec::new();
    for (key, json_str) in &aa.supplementary {
        if key != "spliceAI" {
            continue;
        }
        if let Some(value) = format_spliceai_entry(&allele, json_str) {
            if !values.contains(&value) {
                values.push(value);
            }
        }
    }
    if values.is_empty() { None } else { Some(values.join(",")) }
}

fn format_allele_projection_for_aa(
    vf: &VariationFeature,
    aa: &AlleleAnnotation,
    spec: &VcfProjectionSpec,
) -> Option<String> {
    let allele = uploaded_allele_for_annotation(vf, &aa.allele);
    let mut values: Vec<String> = Vec::new();
    for (key, json_str) in &aa.supplementary {
        if key != spec.json_key {
            continue;
        }
        let parsed = serde_json::from_str::<Value>(json_str)
            .unwrap_or_else(|_| Value::String(json_str.clone()));
        let entries = match spec.kind {
            VcfProjectionKind::AlleleScalar => {
                vec![format!("{}|{}", escape_vcf_subfield(&allele), json_value_to_vcf(&parsed))]
            }
            _ => format_object_projection_entries(&allele, &parsed, spec.fields),
        };
        for value in entries {
            if !values.contains(&value) {
                values.push(value);
            }
        }
    }
    if values.is_empty() { None } else { Some(values.join(",")) }
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
                    parts.push(obj.get(*json_key).map(json_value_to_vcf).unwrap_or_default());
                }
            }
            let value = parts.join("|");
            if !values.contains(&value) {
                values.push(value);
            }
        }
    }
    if values.is_empty() { None } else { Some(values.join(",")) }
}

fn format_clinvar_protein_projection_for_tv(
    vf: &VariationFeature,
    tv: &TranscriptVariation,
    spec: &VcfProjectionSpec,
) -> Option<String> {
    // See `format_gene_projection_for_tv`: fall back to variant-level
    // dedupe when the transcript has no associated gene symbol.
    let symbol_filter = tv.gene_symbol.as_deref();
    let mut values: Vec<String> = Vec::new();
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
        let variants = collect_clinvar_protein_variants(&parsed);
        if variants.is_empty() {
            continue;
        }
        let value = format!("{}|{}", escape_vcf_subfield(&ga.gene_symbol), variants);
        if !values.contains(&value) {
            values.push(value);
        }
    }
    if values.is_empty() { None } else { Some(values.join(",")) }
}

/// Pull every `proteinVariants[*]` entry out of a ClinVar-protein JSON
/// payload and render the `&`-joined `pos:ref>alt:sig` list. Accepts either
/// a single object payload or the JSON array `GeneIndex::annotate_gene`
/// emits when a gene has multiple records.
fn collect_clinvar_protein_variants(parsed: &Value) -> String {
    iter_gene_objects(parsed)
        .into_iter()
        .filter_map(|obj| obj.get("proteinVariants").and_then(|v| v.as_array()))
        .flat_map(|vars| vars.iter())
        .filter_map(|v| {
            let pos = v.get("pos").map(json_leaf_to_string)?;
            let ref_aa = v.get("refAa").map(json_leaf_to_string)?;
            let alt_aa = v.get("altAa").map(json_leaf_to_string)?;
            let sig = v.get("sig").map(json_leaf_to_string).unwrap_or_default();
            Some(escape_vcf_subfield(&format!("{pos}:{ref_aa}>{alt_aa}:{sig}")))
        })
        .collect::<Vec<_>>()
        .join("&")
}

fn format_clinvar_protein_projection(
    vf: &VariationFeature,
    spec: &VcfProjectionSpec,
) -> Option<String> {
    let mut values: Vec<String> = Vec::new();
    for ga in &vf.gene_annotations {
        if ga.json_key != spec.json_key {
            continue;
        }
        let Ok(parsed) = serde_json::from_str::<Value>(&ga.json_string) else {
            continue;
        };
        let variants = collect_clinvar_protein_variants(&parsed);
        if variants.is_empty() {
            continue;
        }
        let value = format!("{}|{}", escape_vcf_subfield(&ga.gene_symbol), variants);
        if !values.contains(&value) {
            values.push(value);
        }
    }
    if values.is_empty() { None } else { Some(values.join(",")) }
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
                    parts.push(obj.get(*json_key).map(json_value_to_vcf).unwrap_or_default());
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

fn escape_vcf_subfield(value: &str) -> String {
    let mut escaped = String::with_capacity(value.len());
    for c in value.chars() {
        match c {
            ':' => escaped.push_str("%3A"),
            ';' => escaped.push_str("%3B"),
            '=' => escaped.push_str("%3D"),
            '%' => escaped.push_str("%25"),
            ',' => escaped.push_str("%2C"),
            '\r' => escaped.push_str("%0D"),
            '\n' => escaped.push_str("%0A"),
            '\t' => escaped.push_str("%09"),
            ' ' => escaped.push_str("%20"),
            '"' => escaped.push_str("%22"),
            '|' => escaped.push_str("%7C"),
            '&' => escaped.push_str("%26"),
            _ => escaped.push(c),
        }
    }
    escaped
}

fn uploaded_allele_for_annotation(vf: &VariationFeature, allele: &Allele) -> String {
    let allele_string = allele.to_string();
    let Some(vcf) = &vf.vcf_fields else {
        return allele_string;
    };

    let uploaded_alts: Vec<&str> = vcf.alt.split(',').collect();
    vf.alt_alleles
        .iter()
        .position(|alt| alt == allele)
        .and_then(|idx| uploaded_alts.get(idx).copied())
        .unwrap_or(&allele_string)
        .to_string()
}

/// Format a VariationFeature as a tab-delimited VEP output line.
///
/// `specs` enumerates which supplementary annotation sources are loaded for
/// this run (build once via `LoadedSupplementarySpecs::new`). For each
/// loaded source, one extra tab column is appended carrying the same
/// pipe-delimited value emitted into the corresponding VCF `FV_*` INFO
/// field (or the `SpliceAI` INFO field for the SpliceAI source). The
/// VCF-side `<ID>=` prefix is *not* included — the column header line
/// already names each column. The column order matches
/// `tab_supplementary_column_names`. Missing values render as `-`.
pub fn format_tab_line(vf: &VariationFeature, specs: &LoadedSupplementarySpecs) -> Vec<String> {
    let mut lines = Vec::new();

    let location = if vf.position.start == vf.position.end {
        format!("{}:{}", vf.position.chromosome, vf.position.start)
    } else {
        format!(
            "{}:{}-{}",
            vf.position.chromosome, vf.position.start, vf.position.end
        )
    };

    let uploaded_variation = vf
        .variation_name
        .clone()
        .unwrap_or_else(|| format!("{}_{}", location, vf.allele_string));

    let extra_count = specs.column_count();
    let row_capacity = 17 + extra_count;

    for tv in &vf.transcript_variations {
        for aa in &tv.allele_annotations {
            let consequence_str = aa
                .consequences
                .iter()
                .map(|c| c.so_term())
                .collect::<Vec<_>>()
                .join(",");

            let impact_str = aa.impact.as_str().to_string();
            let distance_str = aa.distance.map(|d| d.to_string()).unwrap_or("-".to_string());
            let strand_str = format!("{}", tv.strand.as_int());
            let flags_str = if tv.canonical { "canonical" } else { "-" };

            let mut parts: Vec<String> = Vec::with_capacity(row_capacity);
            parts.push(uploaded_variation.clone());
            parts.push(location.clone());
            parts.push(aa.allele.to_string());
            parts.push(tv.gene_id.to_string());
            parts.push(tv.transcript_id.to_string());
            parts.push("Transcript".to_string());
            parts.push(consequence_str);
            parts.push(format_position_range(aa.cdna_position));
            parts.push(format_position_range(aa.cds_position));
            parts.push(format_position_range(aa.protein_position));
            parts.push(
                aa.amino_acids
                    .as_ref()
                    .map(|(r, a)| format!("{}/{}", r, a))
                    .unwrap_or("-".to_string()),
            );
            parts.push(
                aa.codons
                    .as_ref()
                    .map(|(r, a)| format!("{}/{}", r, a))
                    .unwrap_or("-".to_string()),
            );
            parts.push(if aa.existing_variation.is_empty() {
                "-".to_string()
            } else {
                aa.existing_variation.join(",")
            });
            parts.push(impact_str);
            parts.push(distance_str);
            parts.push(strand_str);
            parts.push(flags_str.to_string());

            if extra_count > 0 {
                parts.extend(format_supplementary_tab_columns_for_allele(vf, tv, aa, specs));
            }

            lines.push(parts.join("\t"));
        }
    }

    // If no transcript annotations, still output the variant with intergenic.
    // Pad the row to the full 17 + extra column shape so all tab rows share
    // a header — even when only intergenic alleles are present.
    if vf.transcript_variations.is_empty() {
        for alt in &vf.alt_alleles {
            let mut parts: Vec<String> = Vec::with_capacity(row_capacity);
            parts.push(uploaded_variation.clone());
            parts.push(location.clone());
            parts.push(alt.to_string());
            parts.extend(vec!["-".to_string(); 3]);
            parts.push(Consequence::IntergenicVariant.so_term().to_string());
            parts.extend(vec!["-".to_string(); 10]);
            if extra_count > 0 {
                parts.extend(vec!["-".to_string(); extra_count]);
            }
            lines.push(parts.join("\t"));
        }
    }

    lines
}

/// Format a VariationFeature as JSON.
pub fn format_json(vf: &VariationFeature) -> serde_json::Value {
    let mut obj = serde_json::Map::new();

    obj.insert("id".into(), json_str(&vf.variation_name));
    obj.insert(
        "seq_region_name".into(),
        serde_json::Value::String(vf.position.chromosome.clone()),
    );
    obj.insert("start".into(), serde_json::Value::Number(vf.position.start.into()));
    obj.insert("end".into(), serde_json::Value::Number(vf.position.end.into()));
    obj.insert(
        "allele_string".into(),
        serde_json::Value::String(vf.allele_string.clone()),
    );
    obj.insert("strand".into(), serde_json::Value::Number(vf.position.strand.as_int().into()));

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
                if let Some(ref ms) = tv.mane_select {
                    tc.insert("mane_select".into(), serde_json::Value::String(ms.clone()));
                }
                if let Some(ref mpc) = tv.mane_plus_clinical {
                    tc.insert("mane_plus_clinical".into(), serde_json::Value::String(mpc.clone()));
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
                    tc.insert("gencode_primary".into(), serde_json::Value::Number(1.into()));
                }
                if let Some(ref pid) = tv.protein_id {
                    tc.insert("protein_id".into(), serde_json::Value::String(pid.clone()));
                }
                if let Some((s, e)) = aa.cdna_position {
                    tc.insert("cdna_start".into(), serde_json::Value::Number(s.into()));
                    tc.insert("cdna_end".into(), serde_json::Value::Number(e.into()));
                }
                if let Some((s, e)) = aa.cds_position {
                    tc.insert("cds_start".into(), serde_json::Value::Number(s.into()));
                    tc.insert("cds_end".into(), serde_json::Value::Number(e.into()));
                }
                if let Some((s, e)) = aa.protein_position {
                    tc.insert("protein_start".into(), serde_json::Value::Number(s.into()));
                    tc.insert("protein_end".into(), serde_json::Value::Number(e.into()));
                }
                if let Some(ref aas) = aa.amino_acids {
                    tc.insert("amino_acids".into(),
                        serde_json::Value::String(format!("{}/{}", aas.0, aas.1)));
                }
                if let Some(ref cdns) = aa.codons {
                    tc.insert("codons".into(),
                        serde_json::Value::String(format!("{}/{}", cdns.0, cdns.1)));
                }
                if let Some((n, t)) = aa.exon {
                    tc.insert("exon".into(),
                        serde_json::Value::String(format!("{}/{}", n, t)));
                }
                if let Some((n, t)) = aa.intron {
                    tc.insert("intron".into(),
                        serde_json::Value::String(format!("{}/{}", n, t)));
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
    obj.insert("chromosome".into(), serde_json::Value::String(vf.position.chromosome.clone()));
    obj.insert("position".into(), serde_json::Value::Number(vf.position.start.into()));
    obj.insert("end".into(), serde_json::Value::Number(vf.position.end.into()));

    let ref_str = vf.ref_allele.to_string();
    obj.insert("refAllele".into(), serde_json::Value::String(ref_str));

    let alt_strs: Vec<serde_json::Value> = vf
        .alt_alleles
        .iter()
        .map(|a| serde_json::Value::String(a.to_string()))
        .collect();
    obj.insert("altAlleles".into(), serde_json::Value::Array(alt_strs));

    if vf.variant_type != fastvep_core::VariantType::Unknown {
        obj.insert("variantType".into(), serde_json::Value::String(format!("{:?}", vf.variant_type)));
    }

    // Variants section: one per alt allele with all transcript consequences
    let mut variants = Vec::new();
    for alt in &vf.alt_alleles {
        let alt_str = alt.to_string();
        let mut var_obj = serde_json::Map::new();
        var_obj.insert("altAllele".into(), serde_json::Value::String(alt_str.clone()));

        // Collect transcript consequences for this allele
        let mut transcripts = Vec::new();
        for tv in &vf.transcript_variations {
            for aa in &tv.allele_annotations {
                if aa.allele.to_string() == alt_str {
                    let mut tc = serde_json::Map::new();
                    tc.insert("transcriptId".into(), serde_json::Value::String(tv.transcript_id.to_string()));
                    tc.insert("geneId".into(), serde_json::Value::String(tv.gene_id.to_string()));
                    if let Some(ref sym) = tv.gene_symbol {
                        tc.insert("geneSymbol".into(), serde_json::Value::String(sym.to_string()));
                    }
                    tc.insert("biotype".into(), serde_json::Value::String(tv.biotype.to_string()));
                    tc.insert("consequences".into(), serde_json::Value::Array(
                        aa.consequences.iter().map(|c| serde_json::Value::String(c.so_term().to_string())).collect()
                    ));
                    tc.insert("impact".into(), serde_json::Value::String(aa.impact.as_str().to_string()));
                    if tv.canonical { tc.insert("isCanonical".into(), serde_json::Value::Bool(true)); }
                    if let Some(ref ms) = tv.mane_select { tc.insert("maneSelect".into(), serde_json::Value::String(ms.clone())); }
                    if let Some(ref mpc) = tv.mane_plus_clinical { tc.insert("manePlusClinical".into(), serde_json::Value::String(mpc.clone())); }
                    if let Some(t) = tv.tsl { tc.insert("tsl".into(), serde_json::Value::Number(t.into())); }
                    if let Some(ref a) = tv.appris { tc.insert("appris".into(), serde_json::Value::String(a.clone())); }
                    if let Some(ref c) = tv.ccds { tc.insert("ccds".into(), serde_json::Value::String(c.clone())); }
                    if tv.gencode_primary { tc.insert("isGencodePrimary".into(), serde_json::Value::Bool(true)); }
                    if let Some(ref pid) = tv.protein_id { tc.insert("proteinId".into(), serde_json::Value::String(pid.clone())); }
                    if let Some(ref h) = aa.hgvsg { tc.insert("hgvsg".into(), serde_json::Value::String(h.clone())); }
                    if let Some(ref h) = aa.hgvsc { tc.insert("hgvsc".into(), serde_json::Value::String(h.clone())); }
                    if let Some(ref h) = aa.hgvsp { tc.insert("hgvsp".into(), serde_json::Value::String(h.clone())); }

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
    use super::*;
    use crate::variant::{AlleleAnnotation, TranscriptVariation, VariationFeature};
    use fastvep_core::{
        Allele, Consequence, GeneAnnotation, GenomicPosition, Impact, Strand,
        SupplementaryAnnotation, VariantType,
    };
    use std::sync::Arc;

    #[test]
    fn test_escape_csq_value() {
        assert_eq!(escape_csq_value("hello,world"), "hello&world");
        assert_eq!(escape_csq_value("a;b"), "a%3Bb");
        assert_eq!(escape_csq_value("a|b"), "a&b");
        assert_eq!(escape_csq_value("a b"), "a_b");
        assert_eq!(escape_csq_value("p.Leu153="), "p.Leu153%3D");
    }

    #[test]
    fn test_csq_header() {
        let header = csq_header_line(&["Allele", "Consequence"]);
        assert!(header.contains("Format: Allele|Consequence"));
    }

    #[test]
    fn test_format_position_range() {
        assert_eq!(format_position_range(Some((100, 100))), "100");
        assert_eq!(format_position_range(Some((100, 200))), "100-200");
        assert_eq!(format_position_range(None), "");
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
                    cdna_position: None,
                    cds_position: None,
                    protein_position: None,
                    amino_acids: None,
                    codons: None,
                    exon: None,
                    intron: None,
                    distance: None,
                    hgvsc: None,
                    hgvsp: None,
                    hgvsg: None,
                    hgvs_offset: None,
                    existing_variation: Vec::new(),
                    sift: None,
                    polyphen: None,
                    supplementary,
                    acmg_classification: None,
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

    #[test]
    fn vcf_projection_emits_supported_fastsa_sources_without_json_payloads() {
        let vf = projection_test_variant();
        let projections = format_supplementary_vcf_info(&vf);
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
            assert!(ids.contains(&expected), "missing projection {expected}: {projections:?}");
        }

        let info = projections
            .iter()
            .map(|(id, value)| format!("{id}={value}"))
            .collect::<Vec<_>>()
            .join(";");
        assert!(!info.contains('{'), "VCF INFO must not contain JSON objects: {info}");
        assert!(!info.contains('}'), "VCF INFO must not contain JSON objects: {info}");
        assert!(!info.contains('"'), "VCF INFO must not contain JSON quotes: {info}");
        assert!(info.contains("SpliceAI=G|GENE%7C1|0.01|0.00|0.85|0.00|5|-28|2|-13"));
        assert!(info.contains("FV_CLINVAR=G|Pathogenic&Likely_pathogenic|criteria_provided%2C_multiple_submitters%2C_no_conflicts|Breast%2Ccancer&Ovarian%7Ccancer|SNV|SO%3A0001483"));
        assert!(info.contains("FV_PHYLOP=G|3.14"));
        assert!(info.contains("FV_REVEL=G|0.8123"));
        assert!(info.contains("FV_PRIMATEAI=G|0.4567"));
        assert!(info.contains("FV_OMIM=GENE1|113705|Breast%20cancer&Ovarian%2Ccancer"));
        assert!(info.contains("FV_CLINVAR_PROTEIN=GENE1|175%3AR>H%3APathogenic"));
    }

    #[test]
    fn vcf_info_replaces_existing_fastvep_owned_fields() {
        let vf = projection_test_variant();
        let csq = "G|missense_variant|MODERATE";
        let info = format_vcf_info_fields(
            "DP=12;CSQ=old;SpliceAI=old;FV_CLINVAR=old;KEEP=1",
            &vf,
            csq,
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
        vec!["omim".into(), "gnomad_genes".into(), "clinvar_protein".into()]
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
            assert!(cols.contains(&expected), "missing column {expected}: {cols:?}");
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
        let lines = format_tab_line(&vf, &specs);
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
            format_supplementary_vcf_info(&vf)
                .into_iter()
                .collect();

        for (i, info_id) in extra_cols.iter().enumerate() {
            let tab_value = cols[17 + i];
            let expected = vcf_projections.get(*info_id).cloned().unwrap_or_else(|| "-".into());
            assert_eq!(
                tab_value, expected,
                "tab column {info_id} should equal VCF INFO value"
            );
        }

        assert!(!row.contains('{'), "tab row must not contain JSON: {row}");
        assert!(!row.contains('}'), "tab row must not contain JSON: {row}");
        assert!(!row.contains('"'), "tab row must not contain JSON quotes: {row}");
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

        let specs = LoadedSupplementarySpecs::new(
            &["clinvar".to_string(), "dbnsfp".to_string()],
            &[],
        );
        let lines = format_tab_line(&vf, &specs);
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
        let lines = format_tab_line(&vf, &specs);
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

        let lines = format_tab_line(&vf, &specs);
        assert_eq!(lines.len(), 1, "single alt allele yields one intergenic row");
        let cols: Vec<&str> = lines[0].split('\t').collect();
        assert_eq!(cols.len(), 17 + 3, "intergenic row must include FV_* columns");
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
        let lines = format_tab_line(&vf, &specs);
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
        let lines = format_tab_line(&vf, &specs);
        assert_eq!(lines.len(), 2);
        assert!(
            lines[0].split('\t').last().unwrap().starts_with("G|Pathogenic|"),
            "first row should carry the ALT-G clinvar value: {}",
            lines[0]
        );
        assert!(
            lines[1].split('\t').last().unwrap().starts_with("T|Benign|"),
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
        let lines = format_tab_line(&vf, &specs);
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
            format_supplementary_vcf_info(&vf).into_iter().collect();
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
            cvp.contains("175%3AR>H%3APathogenic"),
            "FV_CLINVAR_PROTEIN should carry first proteinVariants entry: {cvp}"
        );
        assert!(
            cvp.contains("248%3AR>W%3APathogenic"),
            "FV_CLINVAR_PROTEIN should carry second proteinVariants entry: {cvp}"
        );

        // Tab projection — same pipe values land in the per-row cells.
        let specs = LoadedSupplementarySpecs::new(
            &[],
            &["omim".to_string(), "clinvar_protein".to_string()],
        );
        let lines = format_tab_line(&vf, &specs);
        let row = &lines[0];
        let cols: Vec<&str> = row.split('\t').collect();
        let fv_omim = cols[17];
        let fv_cvp = cols[18];
        assert!(
            fv_omim.contains("113705") && fv_omim.contains("167000"),
            "tab FV_OMIM cell should carry both records: {fv_omim}"
        );
        assert!(
            fv_cvp.contains("175%3AR>H") && fv_cvp.contains("248%3AR>W"),
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
        let specs = LoadedSupplementarySpecs::new(
            &["clinvar".to_string(), "dbnsfp".to_string()],
            &[],
        );
        let lines = format_tab_line(&vf, &specs);
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

        let info = format_supplementary_vcf_info(&vf)
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
}
