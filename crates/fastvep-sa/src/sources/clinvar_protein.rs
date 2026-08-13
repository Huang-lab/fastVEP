//! ClinVar protein-position index builder for .oga files.
//!
//! Parses ClinVar VCF to extract protein-level data for classified missense
//! variants, enabling PS1, PM5, and PM1 (hotspot) ACMG criteria evaluation.
//!
//! Both directions of assertion are indexed. The pathogenic entries are what
//! PS1, PM5 and PM1's hotspot count read. The benign ones exist for the second
//! half of PM1's own definition - Richards 2015 asks for a hotspot or critical
//! domain "**without benign variation**" - which cannot be tested from an index
//! that only ever recorded one direction. `benignIndexed` marks a file built this
//! way, so a classifier reading an older index can tell "no benign variation
//! here" apart from "this file never carried any".

use crate::common::{escape_json, GeneRecord};
use anyhow::{Context, Result};
use std::collections::HashMap;
use std::io::BufRead;

/// How ClinVar classifies a missense variant, reduced to the four terms the
/// index carries. Conflicting and uncertain records are not indexed at all:
/// they are evidence of disagreement, not of either direction.
fn normalise_significance(clnsig_lower: &str) -> Option<&'static str> {
    if clnsig_lower.contains("conflicting") {
        return None;
    }
    let pathogenic = clnsig_lower.contains("pathogenic");
    let benign = clnsig_lower.contains("benign");
    // "Likely_pathogenic" contains neither "benign" nor a conflict marker, but
    // a term asserting both directions at once is not usable evidence.
    match (pathogenic, benign) {
        (true, false) if clnsig_lower.contains("likely") => Some("Likely_pathogenic"),
        (true, false) => Some("Pathogenic"),
        (false, true) if clnsig_lower.contains("likely") => Some("Likely_benign"),
        (false, true) => Some("Benign"),
        _ => None,
    }
}

/// Sentinel used while deduplicating, never written to disk: marks a residue
/// substitution that ClinVar asserts in both directions.
const CONFLICTED: &str = "";

/// Whether an indexed significance term is a benign assertion. The four terms
/// the index carries are closed, so this is a total function over them.
fn is_benign(sig: &str) -> bool {
    sig.eq_ignore_ascii_case("Benign") || sig.eq_ignore_ascii_case("Likely_benign")
}

/// A classified missense variant at a specific protein position.
#[derive(Debug, Clone)]
struct ProteinVariant {
    pos: u64,
    ref_aa: String,
    alt_aa: String,
    sig: String,
    /// The nucleotide-level change that produced this amino-acid change, e.g.
    /// `c.5074G>C`. Several distinct nucleotide changes can produce the same
    /// residue substitution, and PS1 is defined precisely on that distinction
    /// ("same amino acid change ... regardless of the nucleotide change").
    /// Empty when the source did not expose a c. token.
    cdna: String,
}

/// Serialize one gene's protein variants into a `GeneRecord`.
///
/// Dedups by `(pos, ref_aa, alt_aa)`, then sorts before serializing: HashMap
/// iteration order is unspecified, so without the sort the `proteinVariants`
/// array order (and thus the on-disk `.oga` bytes) would vary run-to-run for
/// identical input, breaking build reproducibility/diffing. Shared by both the
/// VCF and variant_summary parse paths so the two stay byte-for-byte identical.
///
/// Each entry also carries `n`, the number of *distinct nucleotide changes*
/// that produce this amino-acid change. PS1 needs it: an entry with `n = 1`
/// backed only by the variant being classified is that variant's own record,
/// not independent evidence, and firing PS1 off it is circular. `n` is omitted
/// when it is 1 so older `.oga` files (which default it to 1) stay compatible.
fn serialize_gene_record(gene: String, variants: &[ProteinVariant]) -> GeneRecord {
    use std::collections::BTreeSet;
    let mut unique: HashMap<(u64, String, String), (String, BTreeSet<String>)> = HashMap::new();
    for v in variants {
        let e = unique
            .entry((v.pos, v.ref_aa.clone(), v.alt_aa.clone()))
            .or_insert_with(|| (v.sig.clone(), BTreeSet::new()));
        // Distinct alleles producing the same residue substitution can carry
        // assertions in opposite directions. That is a conflict at the level
        // this index is keyed on, so it is dropped for the same reason a
        // `Conflicting_interpretations` record is: it supports neither side.
        // Marked rather than removed here so a later record cannot resurrect it.
        if is_benign(&e.0) != is_benign(&v.sig) {
            e.0 = CONFLICTED.to_string();
        }
        if !v.cdna.is_empty() {
            e.1.insert(v.cdna.clone());
        }
    }
    unique.retain(|_, (sig, _)| sig != CONFLICTED);

    let mut unique: Vec<_> = unique.into_iter().collect();
    unique.sort_by(|((pos_a, ref_a, alt_a), _), ((pos_b, ref_b, alt_b), _)| {
        pos_a.cmp(pos_b).then(ref_a.cmp(ref_b)).then(alt_a.cmp(alt_b))
    });

    let variant_jsons: Vec<String> = unique
        .iter()
        .map(|((pos, ref_aa, alt_aa), (sig, cdnas))| {
            let n = cdnas.len().max(1);
            let n_field = if n > 1 { format!(r#","n":{}"#, n) } else { String::new() };
            format!(
                r#"{{"pos":{},"refAa":"{}","altAa":"{}","sig":"{}"{}}}"#,
                pos, escape_json(ref_aa), escape_json(alt_aa), escape_json(sig), n_field
            )
        })
        .collect();

    GeneRecord {
        gene_symbol: gene,
        json: format!(
            r#"{{"benignIndexed":true,"proteinVariants":[{}]}}"#,
            variant_jsons.join(",")
        ),
    }
}

/// Parse a ClinVar source (VCF or variant_summary.txt.gz) and produce
/// gene-level records of classified missense variants indexed by protein
/// position. Auto-detects format from the header line:
///
/// - **VCF (`clinvar.vcf.gz`)**: header begins with `#`. Per-record protein
///   change is rarely available in the VCF (ClinVar's MC field is just an SO
///   term, and CLNHGVS is genomic), so the VCF path will yield very few
///   records — the variant_summary path is preferred.
/// - **variant_summary.txt (`variant_summary.txt.gz`)**: header begins with
///   `#AlleleID` and contains a `Name` column with full HGVS like
///   `NM_007294.4(BRCA1):c.5074G>C (p.Asp1692His)`. Protein changes are
///   extracted from the parenthesised `p.` block.
///
/// Output JSON per gene:
/// `{"benignIndexed":true,"proteinVariants":[{"pos":175,"refAa":"R","altAa":"H","sig":"Pathogenic"}, ...]}`
pub fn parse_clinvar_protein_vcf<R: BufRead>(reader: R) -> Result<Vec<GeneRecord>> {
    // Buffer the first line to detect format. Both formats start with `#`,
    // but the variant_summary header begins with `#AlleleID\t...`.
    let mut iter = reader.lines();
    let first = match iter.next() {
        Some(l) => l.context("Reading first line of ClinVar source")?,
        None => return Ok(Vec::new()),
    };
    if first.starts_with("#AlleleID") || first.contains("\tName\t") {
        return parse_variant_summary_inner(first, iter);
    }
    // Otherwise treat as VCF — re-prepend the buffered first line.
    let chained = std::iter::once(Ok(first)).chain(iter);
    parse_clinvar_vcf_inner(chained)
}

fn parse_clinvar_vcf_inner<I>(lines: I) -> Result<Vec<GeneRecord>>
where
    I: Iterator<Item = std::io::Result<String>>,
{
    // Collect classified missense variants per gene, both directions.
    let mut gene_variants: HashMap<String, Vec<ProteinVariant>> = HashMap::new();

    for line in lines {
        let line = line.context("Reading ClinVar VCF line")?;
        if line.starts_with('#') {
            continue;
        }

        let fields: Vec<&str> = line.splitn(9, '\t').collect();
        if fields.len() < 8 {
            continue;
        }

        let info = fields[7];
        let info_map = parse_info(info);

        let clnsig = match info_map.get("CLNSIG") {
            Some(sig) => sig.to_lowercase(),
            None => continue,
        };
        let Some(sig_clean) = normalise_significance(&clnsig) else {
            continue;
        };

        // Check molecular consequence for missense
        let mc = info_map.get("MC").map(|s| s.as_str()).unwrap_or("");
        if !mc.contains("missense_variant") {
            continue;
        }

        // Extract gene symbol from GENEINFO=GENE:ID
        let gene_symbol = match info_map.get("GENEINFO") {
            Some(gi) => {
                if let Some(colon_pos) = gi.find(':') {
                    gi[..colon_pos].to_string()
                } else {
                    gi.to_string()
                }
            }
            None => continue,
        };

        // Try to extract protein change from MC field or CLNHGVS
        // MC format: SO:0001583|missense_variant  (no protein info)
        // CLNHGVS format: NC_000017.11:g.7676154G>A (genomic, not protein)
        // We need to parse protein change from alternative fields

        // Try CLNHGVS for protein notation (some entries have multiple, pipe-separated)
        let mut protein_variant = None;

        // Parse from MC field which sometimes includes protein change info
        // Format can be: "SO:0001583|missense_variant|NP_000537.3:p.Arg175His"
        for part in mc.split(',') {
            if let Some(prot) = extract_protein_from_mc(part) {
                protein_variant = Some(prot);
                break;
            }
        }

        // Also try parsing from CLNHGVS (may contain p. notation)
        if protein_variant.is_none() {
            if let Some(hgvs) = info_map.get("CLNHGVS") {
                for part in hgvs.split(',') {
                    if let Some(prot) = parse_protein_hgvs(part) {
                        protein_variant = Some(prot);
                        break;
                    }
                }
            }
        }

        if let Some(pv) = protein_variant {
            let sig_clean = sig_clean.to_string();
            gene_variants
                .entry(gene_symbol)
                .or_default()
                .push(ProteinVariant {
                    pos: pv.0,
                    ref_aa: pv.1,
                    alt_aa: pv.2,
                    sig: sig_clean,
                    // The VCF path has no c. token (CLNHGVS is genomic), so
                    // distinct nucleotide changes cannot be counted here.
                    // variant_summary is the preferred source and does carry it.
                    cdna: String::new(),
                });
        }
    }

    // Convert to GeneRecords
    let mut records: Vec<GeneRecord> = gene_variants
        .into_iter()
        .map(|(gene, variants)| serialize_gene_record(gene, &variants))
        .collect();

    records.sort_by(|a, b| a.gene_symbol.cmp(&b.gene_symbol));
    Ok(records)
}

/// Parse ClinVar `variant_summary.txt` (or .gz) into protein-position records.
///
/// Header columns of interest (1-based per ClinVar docs, 0-based here):
///  - col 2: `Type` (e.g. "single nucleotide variant")
///  - col 3: `Name` — full HGVS, e.g. `NM_007294.4(BRCA1):c.5074G>C (p.Asp1692His)`
///  - col 5: `GeneSymbol`
///  - col 7: `ClinicalSignificance`
///  - col 25: `ReviewStatus`
///  - col 26: `Assembly` (we keep all rows; protein change is independent of build)
fn parse_variant_summary_inner<I>(header: String, lines: I) -> Result<Vec<GeneRecord>>
where
    I: Iterator<Item = std::io::Result<String>>,
{
    let header_fields: Vec<&str> = header.trim_start_matches('#').split('\t').collect();
    let find = |needle: &str| header_fields.iter().position(|f| f.eq_ignore_ascii_case(needle));
    let i_name = find("Name").context("variant_summary missing Name column")?;
    let i_gene = find("GeneSymbol").context("variant_summary missing GeneSymbol column")?;
    let i_sig = find("ClinicalSignificance").context("variant_summary missing ClinicalSignificance column")?;

    let mut gene_variants: HashMap<String, Vec<ProteinVariant>> = HashMap::new();
    for line in lines {
        let line = line.context("Reading variant_summary line")?;
        if line.starts_with('#') || line.is_empty() {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() <= i_name.max(i_gene).max(i_sig) {
            continue;
        }
        let sig = cols[i_sig].to_lowercase();
        let Some(sig_clean) = normalise_significance(&sig) else {
            continue;
        };
        let gene = cols[i_gene].trim();
        if gene.is_empty() || gene == "-" {
            continue;
        }
        let name = cols[i_name];
        // Extract the parenthesised "p." block from the Name column.
        let pv = match parse_protein_from_summary_name(name) {
            Some(pv) => pv,
            None => continue,
        };

        let sig_clean = sig_clean.to_string();
        // The nucleotide change, e.g. `c.5074G>C` from
        // `NM_007294.4(BRCA1):c.5074G>C (p.Asp1692His)`. variant_summary lists
        // one row per assembly, so deduping on this also collapses the
        // GRCh37/GRCh38 duplicate rows for the same allele.
        let cdna = extract_cdna_token(name);

        // GeneSymbol may be a semicolon-delimited list; split and emit per gene.
        for g in gene.split(';').map(str::trim).filter(|g| !g.is_empty()) {
            gene_variants
                .entry(g.to_string())
                .or_default()
                .push(ProteinVariant {
                    pos: pv.0,
                    ref_aa: pv.1.clone(),
                    alt_aa: pv.2.clone(),
                    sig: sig_clean.clone(),
                    cdna: cdna.clone(),
                });
        }
    }

    let mut records: Vec<GeneRecord> = gene_variants
        .into_iter()
        .map(|(gene, variants)| serialize_gene_record(gene, &variants))
        .collect();
    records.sort_by(|a, b| a.gene_symbol.cmp(&b.gene_symbol));
    Ok(records)
}

/// Pull the `c.` token out of a variant_summary `Name` value like
/// `NM_007294.4(BRCA1):c.5074G>C (p.Asp1692His)`. Returns an empty string when
/// there is no `c.` token to key on.
fn extract_cdna_token(name: &str) -> String {
    let Some(idx) = name.find(":c.") else {
        return String::new();
    };
    let rest = &name[idx + 1..];
    let end = rest.find(' ').unwrap_or(rest.len());
    rest[..end].to_string()
}

/// Pull a `(pos, ref, alt)` protein change from a Name string like
/// `"NM_007294.4(BRCA1):c.5074G>C (p.Asp1692His)"`. Returns None for
/// anything that isn't a missense (silent / synonymous / stop / frameshift /
/// non-protein names).
fn parse_protein_from_summary_name(name: &str) -> Option<(u64, String, String)> {
    // Find the trailing `(p.` block.
    let p_open = name.rfind("(p.")?;
    let after = &name[p_open + 1..];
    let p_close = after.find(')')?;
    parse_protein_hgvs(&after[..p_close])
}

/// Extract protein position and amino acid change from MC field component.
/// MC can contain entries like: "SO:0001583|missense_variant" (no protein info usually)
/// but some ClinVar entries have extended format.
fn extract_protein_from_mc(mc_part: &str) -> Option<(u64, String, String)> {
    // Look for p. notation in the MC field
    if let Some(p_idx) = mc_part.find(":p.") {
        return parse_protein_hgvs(&mc_part[p_idx + 1..]);
    }
    None
}

/// Parse a protein HGVS expression like "p.Arg175His" or "p.R175H"
/// Returns (position, ref_aa, alt_aa) using single-letter codes.
fn parse_protein_hgvs(hgvs: &str) -> Option<(u64, String, String)> {
    let p_str = if let Some(idx) = hgvs.find("p.") {
        &hgvs[idx + 2..]
    } else {
        return None;
    };

    // Try three-letter codes first: "Arg175His"
    if let Some(result) = parse_three_letter_protein(p_str) {
        return Some(result);
    }

    // Try single-letter codes: "R175H"
    if p_str.len() >= 3 {
        let first = p_str.chars().next()?;
        if first.is_ascii_uppercase() {
            // Extract digits
            let digits: String = p_str[1..].chars().take_while(|c| c.is_ascii_digit()).collect();
            if let Ok(pos) = digits.parse::<u64>() {
                let rest = &p_str[1 + digits.len()..];
                if let Some(alt_aa) = rest.chars().next() {
                    if alt_aa.is_ascii_uppercase() && alt_aa != '*' {
                        return Some((pos, first.to_string(), alt_aa.to_string()));
                    }
                }
            }
        }
    }

    None
}

/// Parse three-letter amino acid protein change like "Arg175His"
fn parse_three_letter_protein(s: &str) -> Option<(u64, String, String)> {
    let aa_map: HashMap<&str, &str> = [
        ("Ala", "A"), ("Arg", "R"), ("Asn", "N"), ("Asp", "D"), ("Cys", "C"),
        ("Gln", "Q"), ("Glu", "E"), ("Gly", "G"), ("His", "H"), ("Ile", "I"),
        ("Leu", "L"), ("Lys", "K"), ("Met", "M"), ("Phe", "F"), ("Pro", "P"),
        ("Ser", "S"), ("Thr", "T"), ("Trp", "W"), ("Tyr", "Y"), ("Val", "V"),
        ("Sec", "U"), ("Pyl", "O"), ("Ter", "*"),
    ]
    .iter()
    .copied()
    .collect();

    // Find ref AA (first 3 chars)
    if s.len() < 4 {
        return None;
    }
    let ref_three = s.get(..3)?;
    let ref_aa = aa_map.get(ref_three)?;

    // Extract position digits
    let rest = &s[3..];
    let digits: String = rest.chars().take_while(|c| c.is_ascii_digit()).collect();
    let pos = digits.parse::<u64>().ok()?;

    // Find alt AA
    let after_digits = &rest[digits.len()..];
    if after_digits.len() < 3 {
        return None;
    }
    let alt_three = after_digits.get(..3)?;
    let alt_aa = aa_map.get(alt_three)?;

    // Skip stop codon/terminator variants (not missense)
    if *alt_aa == "*" {
        return None;
    }

    Some((pos, ref_aa.to_string(), alt_aa.to_string()))
}

fn parse_info(info: &str) -> HashMap<String, String> {
    let mut map = HashMap::new();
    for part in info.split(';') {
        if let Some(eq_pos) = part.find('=') {
            let key = &part[..eq_pos];
            let val = &part[eq_pos + 1..];
            map.insert(key.to_string(), val.to_string());
        }
    }
    map
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_protein_hgvs_three_letter() {
        let result = parse_protein_hgvs("p.Arg175His").unwrap();
        assert_eq!(result, (175, "R".to_string(), "H".to_string()));
    }

    #[test]
    fn test_parse_protein_hgvs_single_letter() {
        let result = parse_protein_hgvs("p.R175H").unwrap();
        assert_eq!(result, (175, "R".to_string(), "H".to_string()));
    }

    #[test]
    fn test_parse_protein_hgvs_with_prefix() {
        let result = parse_protein_hgvs("NP_000537.3:p.Arg175His").unwrap();
        assert_eq!(result, (175, "R".to_string(), "H".to_string()));
    }

    #[test]
    fn test_parse_protein_hgvs_stop_codon_rejected() {
        assert!(parse_protein_hgvs("p.Arg175Ter").is_none());
    }

    #[test]
    fn test_parse_three_letter_protein() {
        let result = parse_three_letter_protein("Cys315Met").unwrap();
        assert_eq!(result, (315, "C".to_string(), "M".to_string()));
    }

    #[test]
    fn test_parse_three_letter_protein_multibyte_ref_no_panic() {
        // A multi-byte UTF-8 character straddling byte offset 3 in untrusted
        // ClinVar free text must not panic on the ref-AA slice ("byte index
        // is not a char boundary") — it should just fail to parse.
        // "C" (1 byte) + "€" (3 bytes, occupying byte offsets 1..4) + "s315Met":
        // byte offset 3 lands inside the € character, not on a boundary.
        let s = "C\u{20AC}s315Met";
        assert!(s.len() >= 4);
        assert!(parse_three_letter_protein(s).is_none());
    }

    #[test]
    fn test_parse_three_letter_protein_multibyte_alt_no_panic() {
        // Same hazard on the alt-AA slice: "Cys315" + "M€t" puts a multi-byte
        // character straddling byte offset 3 of `after_digits`.
        let s = "Cys315M\u{20AC}t";
        assert!(parse_three_letter_protein(s).is_none());
    }

    #[test]
    fn test_parse_protein_from_summary_name_typical() {
        let n = "NM_007294.4(BRCA1):c.5074G>C (p.Asp1692His)";
        assert_eq!(
            parse_protein_from_summary_name(n).unwrap(),
            (1692, "D".to_string(), "H".to_string())
        );
    }

    #[test]
    fn test_parse_protein_from_summary_name_silent_skipped() {
        // p.= and p.? are not missense; should yield None
        assert!(parse_protein_from_summary_name("NM_x:c.1A>G (p.=)").is_none());
        assert!(parse_protein_from_summary_name("NM_x:c.1+5G>T (p.?)").is_none());
    }

    #[test]
    fn test_parse_clinvar_protein_variant_summary_format() {
        // Synthetic variant_summary with two pathogenic missense, one silent
        let data = "\
#AlleleID\tType\tName\tGeneID\tGeneSymbol\tHGNC_ID\tClinicalSignificance\tClinSigSimple\tLastEvaluated\tRS# (dbSNP)\tnsv/esv (dbVar)\tRCVaccession\tPhenotypeIDS\tPhenotypeList\tOrigin\tOriginSimple\tAssembly\tChromosomeAccession\tChromosome\tStart\tStop\tReferenceAllele\tAlternateAllele\tCytogenetic\tReviewStatus\tNumberSubmitters\tGuidelines\tTestedInGTR\tOtherIDs\tSubmitterCategories\tVariationID\tPositionVCF\tReferenceAlleleVCF\tAlternateAlleleVCF\tSomaticClinicalImpact\tSomaticClinicalImpactLastEvaluated\tReviewStatusClinicalImpact\tOncogenicity\tOncogenicityLastEvaluated\tReviewStatusOncogenicity
1\tsingle nucleotide variant\tNM_007294.4(BRCA1):c.5074G>C (p.Asp1692His)\t672\tBRCA1\tHGNC:1100\tPathogenic\t1\t-\t-\t-\tRCV000031208\t-\t-\tgermline\tgermline\tGRCh38\tNC_000017.11\t17\t43057105\t43057105\tG\tC\t17q21.31\tcriteria provided, multiple submitters, no conflicts\t5\t-\t-\t-\t1\t1\t43057105\tG\tC\t-\t-\t-\t-\t-\t-
2\tsingle nucleotide variant\tNM_000546.6(TP53):c.524G>A (p.Arg175His)\t7157\tTP53\tHGNC:11998\tPathogenic\t1\t-\t-\t-\tRCV000\t-\t-\tgermline\tgermline\tGRCh38\tNC_000017.11\t17\t7674220\t7674220\tG\tA\t17p13.1\tcriteria provided, multiple submitters, no conflicts\t5\t-\t-\t-\t1\t2\t7674220\tG\tA\t-\t-\t-\t-\t-\t-
3\tsingle nucleotide variant\tNM_000218.3(KCNQ1):c.123C>T (p.=)\t3784\tKCNQ1\tHGNC:6294\tBenign\t0\t-\t-\t-\tRCV000\t-\t-\tgermline\tgermline\tGRCh38\tNC_000011.10\t11\t1\t1\tC\tT\t11p15.5-p15.4\tcriteria provided, single submitter\t1\t-\t-\t-\t1\t3\t1\tC\tT\t-\t-\t-\t-\t-\t-
";
        let records = parse_clinvar_protein_vcf(data.as_bytes()).unwrap();
        // Two genes (BRCA1, TP53) - the KCNQ1 entry is silent (`p.=`), so it
        // has no residue substitution to index in either direction.
        assert_eq!(records.len(), 2);
        let brca1 = records.iter().find(|r| r.gene_symbol == "BRCA1").unwrap();
        assert!(brca1.json.contains("\"pos\":1692"));
        assert!(brca1.json.contains("\"refAa\":\"D\""));
        assert!(brca1.json.contains("\"altAa\":\"H\""));
        assert!(brca1.json.contains("\"sig\":\"Pathogenic\""));
    }

    #[test]
    fn test_parse_clinvar_protein_variant_order_is_deterministic() {
        // A gene with several distinct protein variants must always
        // serialize proteinVariants in the same (sorted-by-position) order
        // regardless of HashMap iteration order — otherwise the .oga bytes
        // would vary run-to-run for identical input.
        let header = "#AlleleID\tType\tName\tGeneID\tGeneSymbol\tHGNC_ID\tClinicalSignificance\tClinSigSimple\tLastEvaluated\tRS# (dbSNP)\tnsv/esv (dbVar)\tRCVaccession\tPhenotypeIDS\tPhenotypeList\tOrigin\tOriginSimple\tAssembly\tChromosomeAccession\tChromosome\tStart\tStop\tReferenceAllele\tAlternateAllele\tCytogenetic\tReviewStatus\tNumberSubmitters\tGuidelines\tTestedInGTR\tOtherIDs\tSubmitterCategories\tVariationID\tPositionVCF\tReferenceAlleleVCF\tAlternateAlleleVCF\tSomaticClinicalImpact\tSomaticClinicalImpactLastEvaluated\tReviewStatusClinicalImpact\tOncogenicity\tOncogenicityLastEvaluated\tReviewStatusOncogenicity";
        let mut rows = vec![header.to_string()];
        // Deliberately out of position order in the input.
        for (variation_id, pos_aa) in [(1, "p.Arg175His"), (2, "p.Gly245Ser"), (3, "p.Cys135Tyr"), (4, "p.Arg273Cys")] {
            rows.push(format!(
                "{vid}\tsingle nucleotide variant\tNM_000546.6(TP53):c.1A>T ({pos_aa})\t7157\tTP53\tHGNC:11998\tPathogenic\t1\t-\t-\t-\tRCV000\t-\t-\tgermline\tgermline\tGRCh38\tNC_000017.11\t17\t1\t1\tG\tA\t17p13.1\tcriteria provided, multiple submitters, no conflicts\t5\t-\t-\t-\t1\t{vid}\t1\tG\tA\t-\t-\t-\t-\t-\t-",
                vid = variation_id,
                pos_aa = pos_aa,
            ));
        }
        let data = rows.join("\n") + "\n";

        let first = parse_clinvar_protein_vcf(data.as_bytes()).unwrap();
        let second = parse_clinvar_protein_vcf(data.as_bytes()).unwrap();
        assert_eq!(
            first[0].json, second[0].json,
            "identical input must produce byte-identical proteinVariants ordering across runs"
        );

        let tp53 = first.iter().find(|r| r.gene_symbol == "TP53").unwrap();
        let pos_135 = tp53.json.find("\"pos\":135").unwrap();
        let pos_175 = tp53.json.find("\"pos\":175").unwrap();
        let pos_245 = tp53.json.find("\"pos\":245").unwrap();
        let pos_273 = tp53.json.find("\"pos\":273").unwrap();
        assert!(
            pos_135 < pos_175 && pos_175 < pos_245 && pos_245 < pos_273,
            "proteinVariants must be sorted by position: {}",
            tp53.json
        );
    }

    /// One variant_summary row, parameterised on the fields these tests vary.
    fn summary_rows(rows: &[(&str, &str, &str)]) -> String {
        let header = "#AlleleID\tType\tName\tGeneID\tGeneSymbol\tHGNC_ID\tClinicalSignificance\tClinSigSimple\tLastEvaluated\tRS# (dbSNP)\tnsv/esv (dbVar)\tRCVaccession\tPhenotypeIDS\tPhenotypeList\tOrigin\tOriginSimple\tAssembly\tChromosomeAccession\tChromosome\tStart\tStop\tReferenceAllele\tAlternateAllele\tCytogenetic\tReviewStatus\tNumberSubmitters\tGuidelines\tTestedInGTR\tOtherIDs\tSubmitterCategories\tVariationID\tPositionVCF\tReferenceAlleleVCF\tAlternateAlleleVCF\tSomaticClinicalImpact\tSomaticClinicalImpactLastEvaluated\tReviewStatusClinicalImpact\tOncogenicity\tOncogenicityLastEvaluated\tReviewStatusOncogenicity";
        let mut out = vec![header.to_string()];
        for (i, (cdna, protein, sig)) in rows.iter().enumerate() {
            out.push(format!(
                "{i}\tsingle nucleotide variant\tNM_000546.6(TP53):{cdna} ({protein})\t7157\tTP53\tHGNC:11998\t{sig}\t1\t-\t-\t-\tRCV000\t-\t-\tgermline\tgermline\tGRCh38\tNC_000017.11\t17\t1\t1\tG\tA\t17p13.1\tcriteria provided, multiple submitters, no conflicts\t5\t-\t-\t-\t1\t{i}\t1\tG\tA\t-\t-\t-\t-\t-\t-"
            ));
        }
        out.join("\n") + "\n"
    }

    #[test]
    fn test_index_carries_both_directions() {
        // PM1's "without benign variation" half cannot be tested from an index
        // that only ever recorded pathogenic assertions.
        let data = summary_rows(&[
            ("c.524G>A", "p.Arg175His", "Pathogenic"),
            ("c.733G>A", "p.Gly245Ser", "Likely pathogenic"),
            ("c.404G>A", "p.Cys135Tyr", "Benign"),
            ("c.817C>T", "p.Arg273Cys", "Likely benign"),
        ]);
        let records = parse_clinvar_protein_vcf(data.as_bytes()).unwrap();
        let tp53 = &records[0];
        assert!(tp53.json.contains(r#""benignIndexed":true"#));
        for sig in ["Pathogenic", "Likely_pathogenic", "Benign", "Likely_benign"] {
            assert!(
                tp53.json.contains(&format!(r#""sig":"{sig}""#)),
                "missing {sig} in {}",
                tp53.json
            );
        }
    }

    #[test]
    fn test_index_drops_uncertain_and_conflicting_records() {
        let data = summary_rows(&[
            ("c.524G>A", "p.Arg175His", "Uncertain significance"),
            ("c.733G>A", "p.Gly245Ser", "Conflicting classifications of pathogenicity"),
        ]);
        let records = parse_clinvar_protein_vcf(data.as_bytes()).unwrap();
        assert!(records.is_empty(), "got {records:?}");
    }

    #[test]
    fn test_a_residue_asserted_both_ways_is_dropped() {
        // Two distinct alleles produce the same substitution and ClinVar
        // classifies them oppositely. That is a conflict at the level this
        // index is keyed on, and it supports neither PM1's hotspot count nor
        // its benign-variation test.
        let data = summary_rows(&[
            ("c.524G>A", "p.Arg175His", "Pathogenic"),
            ("c.524G>C", "p.Arg175His", "Benign"),
            ("c.733G>A", "p.Gly245Ser", "Pathogenic"),
        ]);
        let records = parse_clinvar_protein_vcf(data.as_bytes()).unwrap();
        let tp53 = &records[0];
        assert!(!tp53.json.contains(r#""pos":175"#), "got {}", tp53.json);
        assert!(tp53.json.contains(r#""pos":245"#), "got {}", tp53.json);
    }
}
