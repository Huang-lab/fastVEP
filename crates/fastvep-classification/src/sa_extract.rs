use fastvep_core::{GeneAnnotation, SupplementaryAnnotation};
use fastvep_core::{is_canonical_dinucleotide_offset, parse_intronic_offset};
use fastvep_core::{Consequence, Impact};
use serde::Deserialize;

/// gnomAD population frequency data.
#[derive(Debug, Clone, Default, Deserialize)]
pub struct GnomadData {
    #[serde(rename = "allAf")]
    pub all_af: Option<f64>,
    #[serde(rename = "allAn")]
    pub all_an: Option<u64>,
    #[serde(rename = "allAc")]
    pub all_ac: Option<u64>,
    #[serde(rename = "allHc")]
    pub all_hc: Option<u64>,
    #[serde(rename = "afrAf")]
    pub afr_af: Option<f64>,
    #[serde(rename = "nfeAf")]
    pub nfe_af: Option<f64>,
    #[serde(rename = "easAf")]
    pub eas_af: Option<f64>,
    #[serde(rename = "amrAf")]
    pub amr_af: Option<f64>,
    #[serde(rename = "asjAf")]
    pub asj_af: Option<f64>,
    #[serde(rename = "finAf")]
    pub fin_af: Option<f64>,
    #[serde(rename = "midAf")]
    pub mid_af: Option<f64>,
    #[serde(rename = "othAf")]
    pub oth_af: Option<f64>,
    #[serde(rename = "remainingAf")]
    pub remaining_af: Option<f64>,
    #[serde(rename = "sasAf")]
    pub sas_af: Option<f64>,

    // ── Extended QC / stratified fields ──
    //
    // Every field below is absent from databases built before the gnomAD
    // builder emitted it. `None` / `false` therefore has to mean "this
    // database cannot answer the question", and every consumer must fall back
    // to the behaviour it had before the field existed. Do not read an absent
    // flag as an assertion that the site is clean.
    /// FILTER=AC0: allele count is zero once low-confidence genotypes are removed.
    #[serde(rename = "filterAc0", default)]
    pub filter_ac0: bool,
    /// FILTER=AS_VQSR: the site failed the allele-specific variant-quality model.
    #[serde(rename = "filterVqsr", default)]
    pub filter_vqsr: bool,
    /// FILTER=InbreedingCoeff: genotypes are distributed unlike real population data.
    #[serde(rename = "filterInbreeding", default)]
    pub filter_inbreeding: bool,
    /// The site falls in a low-complexity region.
    #[serde(rename = "lcr", default)]
    pub lcr: bool,
    /// The site falls in a segmental duplication.
    #[serde(rename = "segdup", default)]
    pub segdup: bool,
    /// A sex-chromosome site outside a pseudoautosomal region, where XY samples
    /// are hemizygous and `all_ac_xy` is a count of individuals.
    #[serde(rename = "nonPar", default)]
    pub non_par: bool,
    /// Alternate allele count in XY samples.
    #[serde(rename = "allAcXY")]
    pub all_ac_xy: Option<u64>,
    /// Total allele number in XY samples.
    #[serde(rename = "allAnXY")]
    pub all_an_xy: Option<u64>,
    /// Filtering allele frequency across all samples (Whiffin 2017).
    #[serde(rename = "faf95")]
    pub faf95: Option<f64>,
    /// Maximum filtering allele frequency across genetic-ancestry groups.
    #[serde(rename = "faf95Max")]
    pub faf95_max: Option<f64>,
}

impl GnomadData {
    /// Maximum allele frequency across all populations. Includes both
    /// gnomAD v2.1 codes (`oth`) and v4.1 codes (`mid`, `remaining`).
    pub fn max_pop_af(&self) -> Option<f64> {
        [
            self.all_af, self.afr_af, self.nfe_af, self.eas_af, self.amr_af, self.asj_af,
            self.fin_af, self.mid_af, self.oth_af, self.remaining_af, self.sas_af,
        ]
        .into_iter()
        .flatten()
        .reduce(f64::max)
    }

    /// The name of any FILTER entry that fired on this record, or `None` when
    /// the record passed (or the database predates these fields).
    ///
    /// A non-PASS gnomAD record is not evidence about a real population
    /// frequency, so benign frequency criteria must not be evaluated on one.
    pub fn failed_filter(&self) -> Option<&'static str> {
        if self.filter_ac0 {
            Some("AC0")
        } else if self.filter_vqsr {
            Some("AS_VQSR")
        } else if self.filter_inbreeding {
            Some("InbreedingCoeff")
        } else {
            None
        }
    }

    /// The name of any region flag marking this site as one where short-read
    /// allele frequencies are systematically unreliable.
    ///
    /// This is the per-site form of the curated homologous-gene list
    /// (Mandelker 2016, PMID 27228465): a segmental duplication is exactly the
    /// context in which reads from a paralogue mismap onto the gene of
    /// interest and inflate its apparent frequency.
    pub fn unreliable_region(&self) -> Option<&'static str> {
        if self.segdup {
            Some("segmental duplication")
        } else if self.lcr {
            Some("low-complexity region")
        } else {
            None
        }
    }

    /// Number of hemizygous individuals observed, for a non-PAR sex-chromosome
    /// site.
    ///
    /// Outside the pseudoautosomal regions an XY sample carries a single X (or
    /// Y) allele, so `AC_XY` counts individuals, not chromosomes. `None` for
    /// autosomes, for PAR sites, and for databases built before the XY columns
    /// existed. Returning `None` rather than 0 keeps "we cannot tell" distinct
    /// from "we looked and there are none".
    pub fn hemizygote_count(&self) -> Option<u64> {
        if self.non_par {
            self.all_ac_xy
        } else {
            None
        }
    }

    /// Number of XY individuals surveyed at a non-PAR sex-chromosome site.
    pub fn hemizygote_individuals(&self) -> Option<u64> {
        if self.non_par {
            self.all_an_xy
        } else {
            None
        }
    }

    /// The filtering allele frequency to test benign frequency criteria
    /// against: the maximum across genetic-ancestry groups where available,
    /// falling back to the global FAF.
    ///
    /// The FAF is the lower bound of the 95 % confidence interval on the
    /// frequency (Whiffin 2017, PMID 28518168), so it is robust both to a
    /// founder variant that is common in one population and rare elsewhere,
    /// and to a frequency estimated from very few alleles. Both were raised in
    /// the round-2 medical-genetics review. `None` when the database predates
    /// these columns, in which case callers fall back to the point estimate.
    pub fn filtering_af(&self) -> Option<f64> {
        match (self.faf95_max, self.faf95) {
            (Some(m), Some(g)) => Some(m.max(g)),
            (Some(m), None) => Some(m),
            (None, g) => g,
        }
    }
}

/// ClinVar clinical significance data.
#[derive(Debug, Clone, Default, Deserialize)]
pub struct ClinvarData {
    pub significance: Option<Vec<String>>,
    #[serde(rename = "reviewStatus")]
    pub review_status: Option<String>,
    pub phenotypes: Option<Vec<String>>,
    #[serde(rename = "variantClass")]
    pub variant_class: Option<String>,
    /// ClinVar-distributed population allele frequencies (ExAC / 1000G / ESP).
    /// Absent in caches built before these were emitted → `None`.
    #[serde(rename = "afExac")]
    pub af_exac: Option<f64>,
    #[serde(rename = "afTgp")]
    pub af_tgp: Option<f64>,
    #[serde(rename = "afEsp")]
    pub af_esp: Option<f64>,
}

impl ClinvarData {
    /// Check if any significance term contains a pathogenic classification.
    pub fn has_pathogenic(&self) -> bool {
        self.significance.as_ref().is_some_and(|sigs| {
            sigs.iter().any(|s| {
                let lower = s.to_lowercase();
                lower.contains("pathogenic") && !lower.contains("conflicting")
            })
        })
    }

    /// Check if any significance term contains a benign classification.
    pub fn has_benign(&self) -> bool {
        self.significance.as_ref().is_some_and(|sigs| {
            sigs.iter().any(|s| {
                let lower = s.to_lowercase();
                lower.contains("benign") && !lower.contains("conflicting")
            })
        })
    }

    /// Maximum of the ClinVar-distributed population allele frequencies
    /// (ExAC / 1000G / ESP). `None` when none are present (e.g. caches built
    /// before these fields were emitted). Used as a PM2 frequency backstop
    /// when gnomAD has no record at the variant.
    pub fn max_pop_af(&self) -> Option<f64> {
        [self.af_exac, self.af_tgp, self.af_esp]
            .into_iter()
            .flatten()
            .reduce(f64::max)
    }

    /// The ClinVar significance term marking this variant as low-penetrance or
    /// as a risk allele, if any. ClinVar's controlled vocabulary carries
    /// `Pathogenic, low penetrance`, `Likely pathogenic, low penetrance`,
    /// `Established risk allele`, `Likely risk allele` and `risk factor`.
    ///
    /// Such a variant is expected to be frequent in the population and is
    /// outside BS2's "full penetrance expected at an early age" precondition,
    /// so the frequency-based benign criteria must not fire on it. SERPINA1
    /// PI*Z (1,236 gnomAD homozygotes) and the F2 3'UTR risk allele were both
    /// called Benign this way in the round-2 review.
    pub fn low_penetrance_term(&self) -> Option<String> {
        let sigs = self.significance.as_ref()?;
        sigs.iter()
            .find(|s| {
                let lower = s.to_lowercase().replace('_', " ");
                lower.contains("low penetrance")
                    || lower.contains("risk allele")
                    || lower.contains("risk factor")
            })
            .cloned()
    }

    /// Returns the review star level (0-4).
    pub fn review_stars(&self) -> u8 {
        match self.review_status.as_deref() {
            Some(s) if s.contains("practice_guideline") || s.contains("practice guideline") => 4,
            Some(s)
                if s.contains("reviewed_by_expert_panel")
                    || s.contains("reviewed by expert panel") =>
            {
                3
            }
            Some(s) if s.contains("multiple_submitters") || s.contains("multiple submitters") => 2,
            Some(s)
                if (s.contains("criteria_provided") || s.contains("criteria provided"))
                    && !s.contains("no_assertion") && !s.contains("no assertion") =>
            {
                1
            }
            _ => 0,
        }
    }
}

/// REVEL missense pathogenicity score.
#[derive(Debug, Clone, Deserialize)]
pub struct RevelData {
    pub score: Option<f64>,
}

/// SpliceAI delta scores.
#[derive(Debug, Clone, Default, Deserialize)]
pub struct SpliceAiData {
    pub gene: Option<String>,
    #[serde(rename = "dsAg")]
    pub ds_ag: Option<f64>,
    #[serde(rename = "dsAl")]
    pub ds_al: Option<f64>,
    #[serde(rename = "dsDg")]
    pub ds_dg: Option<f64>,
    #[serde(rename = "dsDl")]
    pub ds_dl: Option<f64>,
    #[serde(rename = "dpAg")]
    pub dp_ag: Option<i32>,
    #[serde(rename = "dpAl")]
    pub dp_al: Option<i32>,
    #[serde(rename = "dpDg")]
    pub dp_dg: Option<i32>,
    #[serde(rename = "dpDl")]
    pub dp_dl: Option<i32>,
}

impl SpliceAiData {
    /// Maximum delta score across all four splice effects.
    pub fn max_delta_score(&self) -> Option<f64> {
        [self.ds_ag, self.ds_al, self.ds_dg, self.ds_dl]
            .into_iter()
            .flatten()
            .reduce(f64::max)
    }
}

/// dbNSFP SIFT/PolyPhen predictions.
#[derive(Debug, Clone, Default, Deserialize)]
pub struct DbNsfpData {
    pub sift: Option<String>,
    pub polyphen: Option<String>,
}

/// Parsed prediction result from dbNSFP format strings.
#[derive(Debug, Clone)]
pub struct PredictionResult {
    pub prediction: String,
    pub score: Option<f64>,
}

impl DbNsfpData {
    /// Parse SIFT prediction from format "deleterious(0.001)" or "tolerated(0.123)".
    pub fn parse_sift(&self) -> Option<PredictionResult> {
        self.sift.as_ref().and_then(|s| parse_prediction_string(s))
    }

    /// Parse PolyPhen prediction from format "probably_damaging(0.998)".
    pub fn parse_polyphen(&self) -> Option<PredictionResult> {
        self.polyphen.as_ref().and_then(|s| parse_prediction_string(s))
    }
}

fn parse_prediction_string(s: &str) -> Option<PredictionResult> {
    if let Some(paren_pos) = s.find('(') {
        let prediction = s[..paren_pos].to_string();
        let score = s[paren_pos + 1..]
            .trim_end_matches(')')
            .parse::<f64>()
            .ok();
        Some(PredictionResult { prediction, score })
    } else {
        Some(PredictionResult {
            prediction: s.to_string(),
            score: None,
        })
    }
}

/// gnomAD gene-level constraint data.
#[derive(Debug, Clone, Default, Deserialize)]
pub struct GnomadGeneData {
    #[serde(rename = "pLI")]
    pub pli: Option<f64>,
    pub loeuf: Option<f64>,
    #[serde(rename = "misZ")]
    pub mis_z: Option<f64>,
    #[serde(rename = "synZ")]
    pub syn_z: Option<f64>,
}

/// OMIM gene-disease data.
#[derive(Debug, Clone, Default, Deserialize)]
pub struct OmimData {
    #[serde(rename = "mimNumber")]
    pub mim_number: Option<u32>,
    pub phenotypes: Option<Vec<String>>,
}

impl OmimData {
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
        self.phenotypes.as_ref().is_some_and(|ps| {
            ps.iter().any(|p| {
                let t = p.trim();
                !t.is_empty() && !is_weak_clingen_classification(t)
            })
        })
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
        self.phenotypes
            .as_ref()
            .is_some_and(|ps| !ps.is_empty() && !self.has_established_relationship())
    }

    /// Check if any phenotype suggests autosomal dominant inheritance.
    ///
    /// Considers only diseases that met the evidence threshold: the inheritance
    /// of a relationship ClinGen refuted is not a fact about this gene, and PM2
    /// and BS2 both branch on it.
    pub fn has_dominant_inheritance(&self) -> bool {
        self.phenotypes.as_ref().is_some_and(|ps| {
            ps.iter().filter(|p| !is_weak_clingen_classification(p)).any(|p| {
                let lower = p.to_lowercase();
                lower.contains("autosomal dominant")
                    || lower.contains("{ad}")
                    || (lower.contains("dominant") && !lower.contains("recessive"))
            })
        })
    }

    /// Check if any phenotype suggests autosomal recessive inheritance.
    ///
    /// Established diseases only, for the same reason as
    /// [`Self::has_dominant_inheritance`].
    pub fn has_recessive_inheritance(&self) -> bool {
        self.phenotypes.as_ref().is_some_and(|ps| {
            ps.iter().filter(|p| !is_weak_clingen_classification(p)).any(|p| {
                let lower = p.to_lowercase();
                lower.contains("autosomal recessive")
                    || lower.contains("{ar}")
                    || (lower.contains("recessive") && !lower.contains("dominant"))
            })
        })
    }
}

/// True for a ClinGen classification that does not meet the SVI evidence
/// threshold: `Limited`, `Disputed`, `Refuted`, `No Known Disease
/// Relationship`. Matched inside the `(ClinGen <Class>/<MOI>, MONDO:...)`
/// suffix the builder writes.
fn is_weak_clingen_classification(phenotype: &str) -> bool {
    let Some(rest) = phenotype.split("(ClinGen ").nth(1) else {
        // No ClinGen tag: an OMIM genemap2 entry. Counts as established.
        return false;
    };
    let class = rest.split('/').next().unwrap_or("").trim();
    matches!(
        class,
        "Limited" | "Disputed" | "Refuted" | "No Known Disease Relationship"
    )
}

/// ClinVar variants indexed by protein position (from .oga).
#[derive(Debug, Clone, Default, Deserialize)]
pub struct ClinvarProteinData {
    #[serde(rename = "proteinVariants")]
    pub protein_variants: Vec<ClinvarProteinVariant>,
    /// Whether the index that produced this record carries benign as well as
    /// pathogenic assertions. PM1's "without benign variation" test needs the
    /// distinction: in an index built before this field existed, the absence of
    /// benign entries is a property of the builder, not of the gene, and
    /// reading it as "no benign variation here" would fire PM1 on exactly the
    /// evidence the criterion asks us to rule out.
    #[serde(rename = "benignIndexed", default)]
    pub benign_indexed: bool,
    /// Whether the index that produced this record carries the canonical
    /// splice-site variants PS1's splice path reads. Same three-state reason as
    /// [`Self::benign_indexed`]: in a file built before the splice pass existed,
    /// an empty [`Self::splice_positions`] says nothing about the gene, so PS1
    /// must report "no data" rather than "no match".
    #[serde(rename = "spliceIndexed", default)]
    pub splice_indexed: bool,
    /// Reference assembly the [`Self::splice_positions`] coordinates are on.
    /// Protein positions are assembly-independent; genomic ones are not, and
    /// silently comparing GRCh38 index positions against GRCh37 calls would
    /// match the wrong sites. Surfaced in PS1's `details` so the build a call
    /// rests on is visible in the output.
    #[serde(rename = "spliceAssembly", default)]
    pub splice_assembly: Option<String>,
    /// Classified ClinVar variants on a canonical splice dinucleotide of this
    /// gene, keyed by genomic position.
    #[serde(rename = "splicePositions", default)]
    pub splice_positions: Vec<ClinvarSpliceVariant>,
}

/// A single classified ClinVar variant sitting on a canonical splice
/// dinucleotide (`+1`, `+2`, `-1` or `-2`).
#[derive(Debug, Clone, Deserialize)]
pub struct ClinvarSpliceVariant {
    /// 1-based genomic position, on the assembly named by
    /// [`ClinvarProteinData::splice_assembly`].
    pub pos: u64,
    #[serde(rename = "ref")]
    pub ref_allele: String,
    pub alt: String,
    /// Signed intronic offset, one of `+1`, `+2`, `-1`, `-2`.
    pub off: i64,
    /// `Pathogenic` or `Likely_pathogenic`.
    pub sig: String,
}

/// Whether a ClinVar `Pathogenic` variant *other than this one* sits on the
/// same canonical splice dinucleotide - PS1's splice path (Walker 2023, ClinGen
/// SVI Splicing Subgroup, PMID 37352859, Table 3).
///
/// Returns `None` when the question cannot be asked: the variant is not a
/// canonical ±1/2 splice variant, no `clinvar_protein` gene record is loaded,
/// or the record predates the splice pass. `Some(false)` is a real answer -
/// the index was consulted and holds no qualifying comparison variant.
///
/// Three deliberate restrictions, each straight out of Table 3:
///
/// - **Comparison variant must be `Pathogenic`, not `Likely_pathogenic`.**
///   Table 3's rows for a canonical-dinucleotide variant under assessment read
///   `N/A` in the LP column. The subgroup's published response to feedback
///   (22 March 2024, item 7c) gives the reason: "since it is so easy for a
///   ±1,2 dinucleotide variant to reach likely pathogenic, we placed
///   constraints on using these variants as reference to make sure there
///   actually was clinical evidence informing that pathogenic classification".
///   Whether an LP call rests on real clinical evidence is a curator's
///   judgement, not something an index can answer.
/// - **Comparison variant must be within the *same* canonical dinucleotide.**
///   Table 3's prerequisite is that the predicted RNA event of the two variants
///   "precisely match". Two variants on the same dinucleotide abolish the same
///   donor or acceptor, so that holds by construction. It does not hold for a
///   comparison variant elsewhere in the splice motif, which is why the rows
///   covering those are left to curators rather than implemented here.
/// - **The variant's own ClinVar record never counts.** A comparison variant is
///   by definition another variant. This is exact here - the index carries
///   ref/alt, so the self entry is identified rather than guessed - which is why
///   it does not consult
///   [`AcmgConfig::exclude_self_from_clinvar_evidence`](crate::config::AcmgConfig::exclude_self_from_clinvar_evidence);
///   that knob exists for the protein index, which cannot tell which entry is
///   the variant being classified.
///
/// The resulting strength is *not* Strong. See
/// [`crate::criteria`] for the PVS1-dependent grading Table 3 specifies.
pub fn same_splice_position_pathogenic(
    consequences: &[Consequence],
    gene_annotations: &[&GeneAnnotation],
    hgvs_c: Option<&str>,
    pos: u64,
    ref_allele: &str,
    alt_allele: &str,
) -> Option<bool> {
    let is_canonical_splice = consequences.iter().any(|c| {
        matches!(
            c,
            Consequence::SpliceAcceptorVariant | Consequence::SpliceDonorVariant
        )
    });
    if !is_canonical_splice {
        return None;
    }

    let cpd: ClinvarProteinData = gene_annotations
        .iter()
        .find(|ga| ga.json_key == "clinvar_protein")
        .and_then(|ga| serde_json::from_str(&ga.json_string).ok())?;
    if !cpd.splice_indexed {
        return None;
    }

    let own_offset = hgvs_c.and_then(parse_intronic_offset);
    let met = cpd.splice_positions.iter().any(|v| {
        v.sig.eq_ignore_ascii_case("Pathogenic")
            && !is_same_allele(v, pos, ref_allele, alt_allele)
            && same_canonical_dinucleotide(pos, own_offset, v.pos, v.off)
    });
    Some(met)
}

/// Whether an index entry is the variant being classified.
fn is_same_allele(v: &ClinvarSpliceVariant, pos: u64, ref_allele: &str, alt_allele: &str) -> bool {
    v.pos == pos
        && v.ref_allele.eq_ignore_ascii_case(ref_allele)
        && v.alt.eq_ignore_ascii_case(alt_allele)
}

/// Whether two canonical splice-site variants sit on the same donor or acceptor
/// dinucleotide.
///
/// Same genomic position settles it outright. Otherwise the two bases of a
/// dinucleotide are genomically adjacent - on either strand, `+1` and `+2` are
/// neighbours, as are `-1` and `-2` - and the offsets must agree in sign, which
/// keeps a donor from pairing with an acceptor. Introns are far longer than the
/// 4 bp that would be needed for a donor and an acceptor dinucleotide to be
/// adjacent, so the sign test is a guard rather than the load-bearing part.
///
/// Without the variant's own offset (no HGVS on the call) only the
/// same-position case can be decided, and the function conservatively declines
/// the adjacent one.
fn same_canonical_dinucleotide(
    vua_pos: u64,
    vua_offset: Option<i64>,
    cmp_pos: u64,
    cmp_offset: i64,
) -> bool {
    if !is_canonical_dinucleotide_offset(cmp_offset) {
        return false;
    }
    if vua_pos == cmp_pos {
        return true;
    }
    let Some(off) = vua_offset.filter(|o| is_canonical_dinucleotide_offset(*o)) else {
        return false;
    };
    vua_pos.abs_diff(cmp_pos) == 1 && off.signum() == cmp_offset.signum()
}

/// A single classified ClinVar variant at a protein position.
#[derive(Debug, Clone, Deserialize)]
pub struct ClinvarProteinVariant {
    pub pos: u64,
    #[serde(rename = "refAa")]
    pub ref_aa: String,
    #[serde(rename = "altAa")]
    pub alt_aa: String,
    pub sig: String,
    /// Number of distinct nucleotide changes in ClinVar that produce this
    /// amino-acid change. PS1 is defined on exactly this distinction ("same
    /// amino acid change ... regardless of the nucleotide change"), so an
    /// entry with `n == 1` backed only by the variant under classification is
    /// that variant's own record rather than independent evidence. Older
    /// `.oga` files omit the field and default to 1.
    #[serde(default = "default_variant_count")]
    pub n: u32,
}

fn default_variant_count() -> u32 {
    1
}

/// Genotype information for a sample at a specific variant.
#[derive(Debug, Clone)]
pub struct GenotypeInfo {
    pub is_het: bool,
    pub is_hom_ref: bool,
    pub is_hom_alt: bool,
    pub is_missing: bool,
    pub is_phased: bool,
    pub depth: Option<u32>,
    pub quality: Option<u32>,
    /// Which alt allele index this genotype carries (1-based). None if hom_ref or missing.
    pub alt_allele_index: Option<u32>,
}

impl GenotypeInfo {
    /// Returns true if genotype passes depth and quality thresholds.
    pub fn passes_quality(&self, min_depth: u32, min_gq: u32) -> bool {
        let depth_ok = self.depth.is_some_and(|d| d >= min_depth);
        let gq_ok = self.quality.is_some_and(|q| q >= min_gq);
        depth_ok && gq_ok
    }

    /// Returns true if the sample carries the variant allele (het or hom_alt).
    pub fn carries_variant(&self) -> bool {
        self.is_het || self.is_hom_alt
    }
}

/// Information about another variant in the same gene for compound-het analysis.
#[derive(Debug, Clone)]
pub struct CompanionVariant {
    /// Whether the companion variant is ClinVar pathogenic
    pub is_clinvar_pathogenic: bool,
    /// Whether the companion variant is ClinVar likely pathogenic.
    /// Used by PM3 v1.0 points scoring (PR7) — LP companions earn fewer
    /// points than P companions. Defaults to false to preserve back-compat.
    pub is_clinvar_likely_pathogenic: bool,
    /// Whether variants are in trans (different haplotypes). None = unphased.
    pub is_in_trans: Option<bool>,
    /// Whether proband is heterozygous for the companion variant
    pub proband_het: bool,
    /// HGVSc of the companion variant for reporting
    pub hgvsc: Option<String>,
}

/// All data needed for ACMG-AMP classification, extracted from pipeline annotations.
#[derive(Debug, Clone)]
pub struct ClassificationInput {
    pub consequences: Vec<Consequence>,
    pub impact: Impact,
    pub gene_symbol: Option<String>,
    pub is_canonical: bool,
    pub amino_acids: Option<(String, String)>,
    pub protein_position: Option<u64>,
    pub gnomad: Option<GnomadData>,
    pub clinvar: Option<ClinvarData>,
    pub revel: Option<RevelData>,
    pub splice_ai: Option<SpliceAiData>,
    pub dbnsfp: Option<DbNsfpData>,
    pub phylop: Option<f64>,
    pub gerp: Option<f64>,
    pub gene_constraints: Option<GnomadGeneData>,
    pub omim: Option<OmimData>,
    /// ClinVar pathogenic variants at protein positions for this gene (from .oga).
    pub clinvar_protein: Option<ClinvarProteinData>,
    /// HGVS c. notation for the variant (e.g. "c.845G>A"). Used by the BA1
    /// exception list lookup (Ghosh 2018) to identify well-known high-AF
    /// pathogenic variants that must not call BA1. `None` if the pipeline
    /// did not produce HGVS — BA1 then proceeds with its default behavior.
    pub hgvs_c: Option<String>,
    // ── PVS1 decision-tree signals (Abou Tayoun 2018) ────────────────────
    // All Optional. The PVS1 evaluator uses them to grade the strength of
    // null-variant evidence (PVS1 / PVS1_Strong / PVS1_Moderate /
    // PVS1_Supporting / N/A). When unpopulated, PVS1 falls back to its
    // legacy binary behavior (full Very Strong if a null variant is in a
    // LOF-intolerant gene). The pipeline (fastvep-cli) computes these from
    // cached transcript exon coordinates + ClinVar protein index.
    /// True if the predicted premature termination is expected to undergo
    /// nonsense-mediated decay, judged by the last-exon proxy: a PTC escapes
    /// only in the 3'-most exon.
    ///
    /// This is the coarser of the two NMD signals. It misses the escape window
    /// at the 3' end of the penultimate exon, and so reports "will decay" for
    /// PTCs that will not. [`Self::nmd_escape_50nt`] measures that properly;
    /// which of the two PVS1 consults is
    /// [`AcmgConfig::pvs1_nmd_50nt_rule`](crate::config::AcmgConfig::pvs1_nmd_50nt_rule).
    pub predicted_nmd: Option<bool>,
    /// Whether the PTC escapes nonsense-mediated decay under the exact 50-nt
    /// rule (Abou Tayoun 2018): escape in the last exon, in the last 50 nt of
    /// the penultimate exon, or in a single-exon transcript.
    ///
    /// `None` when the pipeline could not measure it - an intronic variant with
    /// no cDNA coordinate, a non-coding transcript, or a caller that does not
    /// populate the field. PVS1 then falls back to [`Self::predicted_nmd`].
    pub nmd_escape_50nt: Option<bool>,
    /// Fraction of the protein removed by the variant (0.0–1.0).
    pub protein_truncation_pct: Option<f64>,
    /// True if the variant lies in the 3'-most (last) exon.
    pub is_last_exon: Option<bool>,
    /// True if downstream pathogenic variants exist past the truncation point
    /// (used as a proxy for "critical functional region").
    pub in_critical_region: Option<bool>,
    /// Distance (in codons) to the next downstream Met for start-loss
    /// variants. None if no alternative start codon exists.
    pub alt_start_codon_distance: Option<i64>,
    /// PS1 splice path (Walker 2023): true when a ClinVar `Pathogenic` variant
    /// other than this one sits on the same canonical splice dinucleotide, and
    /// so is predicted to abolish the same donor or acceptor. `None` when the
    /// question could not be asked. Computed by the caller from the ClinVar
    /// splice index - see [`same_splice_position_pathogenic`], whose doc
    /// comment carries the Table 3 reasoning.
    pub same_splice_position_pathogenic: Option<bool>,
    /// Whether variant overlaps a repeat region (from RepeatMasker .osi).
    pub in_repeat_region: Option<bool>,
    /// The variant's VCF-form position as `(chrom, pos, ref, alt)`, on the
    /// build the run is annotating against.
    ///
    /// Most criteria reason about a transcript-level view and need no
    /// coordinate. BA1's exception list does: it is a list of specific
    /// variants, and an HGVS `c.` token only identifies one relative to the
    /// transcript it was written for. `None` when the caller did not supply
    /// one, in which case the exception list falls back to matching on
    /// gene plus HGVS.
    pub variant_coordinates: Option<(String, u64, String, String)>,
    /// True when the variant is a pure insertion or duplication (VCF REF is a
    /// single anchor base that the ALT extends). Used by PVS1's canonical
    /// splice track: an insertion adjacent to, or inside, the ±1/±2
    /// dinucleotide does not necessarily destroy it. `c.802-2dupA` inserts the
    /// same base the acceptor already has, so the intron still ends in AG.
    /// `None` when the pipeline did not supply allele shape.
    pub is_pure_insertion: Option<bool>,
    /// Whether the variant sits at the first base or last 3 bases of an exon
    /// (the canonical splice region). Per Walker 2023 (ClinGen SVI Splicing
    /// Subgroup), BP7 must NOT fire for synonymous variants at these positions
    /// because SpliceAI can miss splice impact in this region. `None` means
    /// the pipeline didn't populate this signal (BP7 falls back to legacy
    /// behavior — fire if SpliceAI low + PhyloP low).
    pub at_exon_edge: Option<bool>,
    /// Signed offset in bp from the nearest exon boundary, for intronic
    /// variants. Convention: positive after the donor (e.g. c.123+5 → +5),
    /// negative before the acceptor (e.g. c.123-15 → -15). Per Walker 2023,
    /// BP7 may extend to intronic variants outside the standard splice region:
    /// donor-side `offset ≥ 7` or acceptor-side `offset ≤ -21`. `None` for
    /// non-intronic variants or when the pipeline can't compute it.
    pub intronic_offset: Option<i64>,
    /// Proband genotype information (from trio VCF)
    pub proband_genotype: Option<GenotypeInfo>,
    /// Mother genotype information (from trio VCF)
    pub mother_genotype: Option<GenotypeInfo>,
    /// Father genotype information (from trio VCF)
    pub father_genotype: Option<GenotypeInfo>,
    /// Other variants in the same gene that the proband carries (for compound-het)
    pub companion_variants: Vec<CompanionVariant>,
    /// Curated functional-assay result for this variant, if the run was given a
    /// `--functional-evidence` file that lists it.
    ///
    /// Resolved by the caller rather than here, because the lookup is by
    /// genomic coordinate and this function is handed a transcript-level view
    /// that has no coordinates in it. Keeping the file I/O and the coordinate
    /// matching in the pipeline also leaves the classifier consuming typed
    /// evidence, the same as every other criterion.
    pub functional_evidence: Option<crate::functional::FunctionalEvidence>,
}

/// True when the reference allele is the empty/"-" allele, i.e. the variant is
/// a pure insertion or duplication. PVS1's canonical splice track consults this
/// because inserting bases beside the ±1/±2 dinucleotide does not necessarily
/// destroy it (a duplication of the acceptor's own base leaves the intron still
/// ending in AG).
pub fn is_pure_insertion(ref_allele: &fastvep_core::Allele) -> Option<bool> {
    Some(*ref_allele == fastvep_core::Allele::Deletion)
}

/// Extract classification input from pipeline annotation data.
///
/// Parses the pre-serialized JSON strings from supplementary annotations into typed structs.
#[allow(clippy::too_many_arguments)]
pub fn extract_classification_input(
    consequences: &[Consequence],
    impact: Impact,
    gene_symbol: Option<&str>,
    is_canonical: bool,
    amino_acids: Option<&(String, String)>,
    protein_position: Option<u64>,
    hgvs_c: Option<&str>,
    exon: Option<(u32, u32)>,
    protein_length: Option<u64>,
    escapes_nmd: Option<bool>,
    repeat_db_loaded: bool,
    same_splice_position_pathogenic: Option<bool>,
    variant_coordinates: Option<(String, u64, String, String)>,
    is_pure_insertion: Option<bool>,
    functional_evidence: Option<crate::functional::FunctionalEvidence>,
    allele_supplementary: &[(String, String)],
    gene_annotations: &[&GeneAnnotation],
    variant_supplementary: &[SupplementaryAnnotation],
    proband_genotype: Option<GenotypeInfo>,
    mother_genotype: Option<GenotypeInfo>,
    father_genotype: Option<GenotypeInfo>,
    companion_variants: Vec<CompanionVariant>,
) -> ClassificationInput {
    let mut gnomad = None;
    let mut clinvar = None;
    let mut revel = None;
    let mut splice_ai = None;
    let mut dbnsfp = None;

    // Parse per-allele supplementary annotations
    for (key, json_str) in allele_supplementary {
        match key.as_str() {
            "gnomad" => {
                gnomad = serde_json::from_str(json_str).ok();
            }
            "clinvar" => {
                clinvar = serde_json::from_str(json_str).ok();
            }
            "revel" => {
                revel = serde_json::from_str(json_str).ok();
            }
            // SpliceAI ships under several capitalisations across builders
            // (`spliceAI` from the SaWriter pipeline, `spliceai` lowercased,
            // `splice_ai` snake_case). All three resolve to the same struct.
            "spliceai" | "spliceAi" | "spliceAI" | "splice_ai" => {
                splice_ai = serde_json::from_str(json_str).ok();
            }
            "dbnsfp" => {
                dbnsfp = serde_json::from_str(json_str).ok();
            }
            _ => {}
        }
    }

    // Parse positional supplementary annotations (PhyloP, GERP). The CLI
    // pipeline attaches positional SAs into `aa.supplementary` (allele-level)
    // alongside the per-allele SAs above, so we look in both places —
    // whichever fires first wins and the other becomes a no-op.
    let mut phylop = None;
    let mut gerp = None;
    for (key, json_str) in allele_supplementary {
        match key.as_str() {
            "phylop" | "phyloP" => {
                if phylop.is_none() {
                    phylop = json_str.trim_matches('"').parse::<f64>().ok();
                }
            }
            "gerp" | "GERP" => {
                if gerp.is_none() {
                    gerp = json_str.trim_matches('"').parse::<f64>().ok();
                }
            }
            _ => {}
        }
    }
    for sa in variant_supplementary {
        match sa.json_key.as_str() {
            "phylop" | "phyloP" => {
                if phylop.is_none() {
                    phylop = sa.json_string.trim_matches('"').parse::<f64>().ok();
                }
            }
            "gerp" | "GERP" => {
                if gerp.is_none() {
                    gerp = sa.json_string.trim_matches('"').parse::<f64>().ok();
                }
            }
            _ => {}
        }
    }

    // Parse gene-level annotations
    let mut gene_constraints = None;
    let mut omim = None;
    let mut clinvar_protein: Option<ClinvarProteinData> = None;
    for ga in gene_annotations {
        match ga.json_key.as_str() {
            "gnomad_genes" | "gnomad_gene" => {
                gene_constraints = serde_json::from_str(&ga.json_string).ok();
            }
            "omim" => {
                omim = serde_json::from_str(&ga.json_string).ok();
            }
            "clinvar_protein" => {
                clinvar_protein = serde_json::from_str(&ga.json_string).ok();
            }
            _ => {}
        }
    }

    // Whether the variant overlaps a repeat region, from an interval `.osi`.
    //
    // Three-state on purpose. An interval source only yields an annotation when
    // the position falls inside an interval, so a miss and an unloaded database
    // look identical from `allele_supplementary` alone - which is why the
    // caller has to say whether a repeat source was loaded at all. Without that
    // distinction BP3 reported "repeat region data not available (load
    // RepeatMasker .osi)" for every in-frame indel that simply was not in a
    // repeat, which reads as a setup error rather than as the answer.
    //
    // `None` still means "cannot tell", and BP3 declines to evaluate on it.
    let in_repeat_region = {
        let hit = allele_supplementary.iter().any(|(key, _)| {
            let k = key.to_lowercase();
            k.contains("repeat") || k.contains("repeatmasker") || k.contains("simple_repeat")
        });
        match (hit, repeat_db_loaded) {
            (true, _) => Some(true),
            (false, true) => Some(false),
            (false, false) => None,
        }
    };

    // ── PVS1 / BP7 decision-tree signals derived from data the pipeline
    // already computes (exon rank/total + HGVS c. notation). ────────────────
    // `is_last_exon`: the variant sits in the 3'-most exon (rank == total).
    let is_last_exon = exon.and_then(|(rank, total)| (total > 0).then_some(rank == total));
    // `intronic_offset`: distance from the nearest exon boundary, parsed from
    // the HGVS `+N`/`-N` token. None for purely exonic variants.
    let intronic_offset = hgvs_c.and_then(parse_intronic_offset);
    // Both NMD signals are carried, and PVS1 picks between them per config.
    // The last-exon proxy is the historical one; the 50-nt measurement is the
    // rule Abou Tayoun 2018 actually states. They disagree only for a PTC in
    // the 3' end of the penultimate exon, where the proxy claims decay and the
    // measurement finds escape.
    let predicted_nmd = is_last_exon.map(|last| !last);

    // `protein_truncation_pct`: the fraction of the peptide the truncation
    // removes, counting the first lost residue. The tree branches at 10%.
    let protein_truncation_pct = match (protein_position, protein_length) {
        (Some(pos), Some(len)) if len > 0 && pos <= len => {
            Some((len - pos + 1) as f64 / len as f64)
        }
        _ => None,
    };

    // `in_critical_region`: the tree asks whether the removed portion of the
    // protein is known to matter. The .oga ClinVar protein index answers it
    // directly - a pathogenic variant at a residue downstream of the truncation
    // point is positive evidence that the lost region is functional. Left
    // `None` when the gene has no index entry, so absence of the file never
    // reads as "the region is unimportant".
    let in_critical_region = match (protein_position, clinvar_protein.as_ref()) {
        (Some(pos), Some(idx)) => Some(idx.protein_variants.iter().any(|v| {
            let sig = v.sig.to_ascii_lowercase();
            v.pos > pos && sig.contains("pathogenic") && !sig.contains("conflicting")
        })),
        _ => None,
    };

    ClassificationInput {
        consequences: consequences.to_vec(),
        impact,
        gene_symbol: gene_symbol.map(|s| s.to_string()),
        is_canonical,
        amino_acids: amino_acids.cloned(),
        protein_position,
        gnomad,
        clinvar,
        revel,
        splice_ai,
        dbnsfp,
        phylop,
        gerp,
        gene_constraints,
        omim,
        clinvar_protein,
        // Threaded from the caller (typically `aa.hgvsc` from the annotation
        // context). When the pipeline doesn't compute HGVS — i.e. the user
        // didn't pass `--hgvs` — this stays `None` and BA1 falls back to its
        // default behavior (no exception-list lookup).
        hgvs_c: hgvs_c.map(|s| s.to_string()),
        // PVS1 decision-tree signals, all derived above.
        // `alt_start_codon_distance` (start-loss track) still awaits dedicated
        // plumbing and stays `None`, so that branch keeps its conservative
        // default. `same_splice_position_pathogenic` is resolved by the caller,
        // because the ClinVar splice index is keyed by genomic coordinate and
        // this function is handed a transcript-level view with no coordinates
        // in it - the same division of labour as `functional_evidence`.
        predicted_nmd,
        nmd_escape_50nt: escapes_nmd,
        protein_truncation_pct,
        is_last_exon,
        in_critical_region,
        alt_start_codon_distance: None,
        same_splice_position_pathogenic,
        variant_coordinates,
        in_repeat_region,
        is_pure_insertion,
        // BP7 exon-edge / deep-intronic signals (Walker 2023). `intronic_offset`
        // is derived above; `at_exon_edge` awaits exon-coordinate plumbing and
        // stays None (BP7 exon-edge exclusion falls back to legacy behavior).
        at_exon_edge: None,
        intronic_offset,
        proband_genotype,
        mother_genotype,
        father_genotype,
        companion_variants,
        functional_evidence,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gnomad_max_pop_af() {
        let g = GnomadData {
            all_af: Some(0.001),
            afr_af: Some(0.05),
            nfe_af: Some(0.0001),
            ..Default::default()
        };
        assert!((g.max_pop_af().unwrap() - 0.05).abs() < 1e-10);
    }

    #[test]
    fn test_gnomad_max_pop_af_none() {
        let g = GnomadData::default();
        assert!(g.max_pop_af().is_none());
    }

    #[test]
    fn test_clinvar_pathogenic() {
        let c = ClinvarData {
            significance: Some(vec!["Pathogenic".to_string()]),
            ..Default::default()
        };
        assert!(c.has_pathogenic());
        assert!(!c.has_benign());
    }

    #[test]
    fn test_clinvar_likely_pathogenic() {
        let c = ClinvarData {
            significance: Some(vec!["Likely_pathogenic".to_string()]),
            ..Default::default()
        };
        assert!(c.has_pathogenic());
    }

    #[test]
    fn test_clinvar_conflicting_not_pathogenic() {
        let c = ClinvarData {
            significance: Some(vec![
                "Conflicting_interpretations_of_pathogenicity".to_string(),
            ]),
            ..Default::default()
        };
        assert!(!c.has_pathogenic());
    }

    #[test]
    fn test_clinvar_review_stars() {
        assert_eq!(
            ClinvarData {
                review_status: Some(
                    "criteria_provided,_multiple_submitters,_no_conflicts".to_string()
                ),
                ..Default::default()
            }
            .review_stars(),
            2
        );
        assert_eq!(
            ClinvarData {
                review_status: Some("reviewed_by_expert_panel".to_string()),
                ..Default::default()
            }
            .review_stars(),
            3
        );
        assert_eq!(
            ClinvarData {
                review_status: Some("criteria_provided,_single_submitter".to_string()),
                ..Default::default()
            }
            .review_stars(),
            1
        );
        assert_eq!(
            ClinvarData {
                review_status: Some("no_assertion_criteria_provided".to_string()),
                ..Default::default()
            }
            .review_stars(),
            0
        );
    }

    #[test]
    fn test_spliceai_max_delta() {
        let s = SpliceAiData {
            ds_ag: Some(0.01),
            ds_al: Some(0.85),
            ds_dg: Some(0.02),
            ds_dl: Some(0.10),
            ..Default::default()
        };
        assert!((s.max_delta_score().unwrap() - 0.85).abs() < 1e-10);
    }

    #[test]
    fn test_parse_prediction_string() {
        let r = parse_prediction_string("deleterious(0.001)").unwrap();
        assert_eq!(r.prediction, "deleterious");
        assert!((r.score.unwrap() - 0.001).abs() < 1e-10);

        let r = parse_prediction_string("probably_damaging(0.998)").unwrap();
        assert_eq!(r.prediction, "probably_damaging");
        assert!((r.score.unwrap() - 0.998).abs() < 1e-10);

        let r = parse_prediction_string("tolerated").unwrap();
        assert_eq!(r.prediction, "tolerated");
        assert!(r.score.is_none());
    }

    #[test]
    fn test_gnomad_deserialization() {
        let json = r#"{"allAf":1.234e-03,"allAn":150000,"allAc":150,"allHc":2,"afrAf":2.0e-03,"nfeAf":5.0e-04}"#;
        let g: GnomadData = serde_json::from_str(json).unwrap();
        assert!((g.all_af.unwrap() - 0.001234).abs() < 1e-8);
        assert_eq!(g.all_an.unwrap(), 150000);
        assert!((g.afr_af.unwrap() - 0.002).abs() < 1e-8);
    }

    #[test]
    fn test_clinvar_deserialization() {
        let json = r#"{"significance":["Pathogenic"],"reviewStatus":"criteria_provided,_multiple_submitters,_no_conflicts","phenotypes":["Breast_cancer"],"variantClass":"SNV"}"#;
        let c: ClinvarData = serde_json::from_str(json).unwrap();
        assert!(c.has_pathogenic());
        assert_eq!(c.review_stars(), 2);
        // No AF keys → max_pop_af is None (back-compat with old caches).
        assert_eq!(c.max_pop_af(), None);
    }

    #[test]
    fn test_clinvar_population_af_parsing() {
        let json = r#"{"significance":["Pathogenic"],"afExac":0.0113,"afEsp":0.002}"#;
        let c: ClinvarData = serde_json::from_str(json).unwrap();
        assert_eq!(c.af_exac, Some(0.0113));
        assert_eq!(c.af_esp, Some(0.002));
        assert_eq!(c.af_tgp, None);
        assert_eq!(c.max_pop_af(), Some(0.0113));
    }

    /// Build a `ClassificationInput` through the real entry point, varying only
    /// the PVS1 tree inputs. Everything else is inert.
    fn extract(
        protein_position: Option<u64>,
        protein_length: Option<u64>,
        escapes_nmd: Option<bool>,
        exon: Option<(u32, u32)>,
        clinvar_protein_json: Option<&str>,
    ) -> ClassificationInput {
        let gene_anns: Vec<GeneAnnotation> = clinvar_protein_json
            .map(|j| {
                vec![GeneAnnotation {
                    gene_symbol: "TEST".to_string(),
                    json_key: "clinvar_protein".to_string(),
                    json_string: j.to_string(),
                }]
            })
            .unwrap_or_default();
        let refs: Vec<&GeneAnnotation> = gene_anns.iter().collect();
        extract_classification_input(
            &[Consequence::StopGained],
            Impact::High,
            Some("TEST"),
            true,
            None,
            protein_position,
            None,
            exon,
            protein_length,
            escapes_nmd,
            false, // repeat_db_loaded: no interval source in these unit tests
            None,  // same_splice_position_pathogenic: resolved by the caller
            None,  // variant_coordinates: these unit tests are transcript-level
            None,
            None,
            &[],
            &refs,
            &[],
            None,
            None,
            None,
            vec![],
        )
    }

    #[test]
    fn test_both_nmd_signals_are_carried_separately() {
        // Exon 3 of 10 reads as NMD-competent under the last-exon proxy, while
        // the measured 50-nt rule found the position inside the escape window
        // at the end of the penultimate exon. Both are reported; PVS1 picks
        // between them per `pvs1_nmd_50nt_rule`, so neither may overwrite the
        // other here.
        let input = extract(None, None, Some(true), Some((3, 10)), None);
        assert_eq!(input.predicted_nmd, Some(true), "last-exon proxy");
        assert_eq!(input.nmd_escape_50nt, Some(true), "measured 50-nt rule");
        assert_eq!(input.is_last_exon, Some(false));
    }

    #[test]
    fn test_predicted_nmd_is_the_last_exon_proxy() {
        assert_eq!(extract(None, None, None, Some((10, 10)), None).predicted_nmd, Some(false));
        assert_eq!(extract(None, None, None, Some((3, 10)), None).predicted_nmd, Some(true));
        assert_eq!(extract(None, None, None, None, None).predicted_nmd, None);
    }

    #[test]
    fn test_protein_truncation_pct_counts_the_first_lost_residue() {
        // Truncating at residue 91 of 100 removes residues 91..=100, i.e. 10 %.
        let input = extract(Some(91), Some(100), None, None, None);
        assert_eq!(input.protein_truncation_pct, Some(0.10));
        // Out-of-range or absent inputs produce no claim rather than a wrong one.
        assert_eq!(extract(Some(101), Some(100), None, None, None).protein_truncation_pct, None);
        assert_eq!(extract(Some(50), None, None, None, None).protein_truncation_pct, None);
    }

    #[test]
    fn test_in_critical_region_reads_downstream_clinvar_pathogenic_variants() {
        let idx = r#"{"proteinVariants":[
            {"pos":900,"refAa":"R","altAa":"Q","sig":"Pathogenic","n":2},
            {"pos":100,"refAa":"A","altAa":"T","sig":"Pathogenic","n":2}
        ]}"#;
        // Truncating at 500 loses residue 900, which ClinVar calls pathogenic.
        assert_eq!(extract(Some(500), None, None, None, Some(idx)).in_critical_region, Some(true));
        // Truncating at 950 loses nothing ClinVar has an entry for.
        assert_eq!(extract(Some(950), None, None, None, Some(idx)).in_critical_region, Some(false));
    }

    #[test]
    fn test_in_critical_region_is_unknown_without_an_index() {
        // Absence of the .oga must never read as "the lost region is unimportant".
        assert_eq!(extract(Some(500), None, None, None, None).in_critical_region, None);
    }

    #[test]
    fn test_in_critical_region_ignores_conflicting_clinvar_entries() {
        let idx = r#"{"proteinVariants":[
            {"pos":900,"refAa":"R","altAa":"Q","sig":"Conflicting_interpretations_of_pathogenicity","n":1}
        ]}"#;
        assert_eq!(extract(Some(500), None, None, None, Some(idx)).in_critical_region, Some(false));
    }

    // ── PS1 splice path (Walker 2023 Table 3) ────────────────────────────

    /// A gene record carrying the splice index, with `rows` as the raw JSON of
    /// `splicePositions`.
    fn splice_index(rows: &str) -> GeneAnnotation {
        GeneAnnotation {
            gene_symbol: "TEST".to_string(),
            json_key: "clinvar_protein".to_string(),
            json_string: format!(
                r#"{{"benignIndexed":true,"proteinVariants":[],"spliceIndexed":true,"spliceAssembly":"GRCh38","splicePositions":[{rows}]}}"#
            ),
        }
    }

    fn splice_ps1(index: Option<&GeneAnnotation>, hgvs: Option<&str>, pos: u64, r: &str, a: &str) -> Option<bool> {
        let anns: Vec<&GeneAnnotation> = index.into_iter().collect();
        same_splice_position_pathogenic(
            &[Consequence::SpliceAcceptorVariant],
            &anns,
            hgvs,
            pos,
            r,
            a,
        )
    }

    #[test]
    fn test_ps1_splice_fires_on_a_different_allele_at_the_same_position() {
        let idx = splice_index(r#"{"pos":100,"ref":"A","alt":"G","off":-2,"sig":"Pathogenic"}"#);
        assert_eq!(splice_ps1(Some(&idx), Some("c.376-2A>T"), 100, "A", "T"), Some(true));
    }

    #[test]
    fn test_ps1_splice_never_fires_off_the_variants_own_record() {
        // The index is built from ClinVar, so a variant that is itself ClinVar
        // Pathogenic finds its own entry. A comparison variant is by
        // definition another variant, and firing off yourself is circular.
        let idx = splice_index(r#"{"pos":100,"ref":"A","alt":"G","off":-2,"sig":"Pathogenic"}"#);
        assert_eq!(splice_ps1(Some(&idx), Some("c.376-2A>G"), 100, "A", "G"), Some(false));
    }

    #[test]
    fn test_ps1_splice_reaches_the_other_base_of_the_same_dinucleotide() {
        // c.376-1 and c.376-2 are the two bases of one acceptor, genomically
        // adjacent on either strand.
        let idx = splice_index(r#"{"pos":101,"ref":"G","alt":"A","off":-1,"sig":"Pathogenic"}"#);
        assert_eq!(splice_ps1(Some(&idx), Some("c.376-2A>T"), 100, "A", "T"), Some(true));
    }

    #[test]
    fn test_ps1_splice_does_not_cross_from_a_donor_to_an_acceptor() {
        // Opposite-signed offsets are different splice sites even if the
        // coordinates were somehow adjacent.
        let idx = splice_index(r#"{"pos":101,"ref":"G","alt":"A","off":1,"sig":"Pathogenic"}"#);
        assert_eq!(splice_ps1(Some(&idx), Some("c.376-2A>T"), 100, "A", "T"), Some(false));
    }

    #[test]
    fn test_ps1_splice_rejects_a_likely_pathogenic_comparison() {
        // Table 3 reads N/A in the LP column for a canonical-dinucleotide
        // variant under assessment: an LP call at a ±1,2 position is too easy
        // to reach for its clinical evidence to be borrowable unchecked.
        let idx = splice_index(r#"{"pos":100,"ref":"A","alt":"G","off":-2,"sig":"Likely_pathogenic"}"#);
        assert_eq!(splice_ps1(Some(&idx), Some("c.376-2A>T"), 100, "A", "T"), Some(false));
    }

    #[test]
    fn test_ps1_splice_without_an_index_is_unknown_not_absent() {
        assert_eq!(splice_ps1(None, Some("c.376-2A>T"), 100, "A", "T"), None);

        // A record built before the splice pass carries no `spliceIndexed`,
        // and its silence says nothing about the gene.
        let old = GeneAnnotation {
            gene_symbol: "TEST".to_string(),
            json_key: "clinvar_protein".to_string(),
            json_string: r#"{"benignIndexed":true,"proteinVariants":[]}"#.to_string(),
        };
        assert_eq!(splice_ps1(Some(&old), Some("c.376-2A>T"), 100, "A", "T"), None);
    }

    #[test]
    fn test_ps1_splice_declines_a_non_canonical_variant() {
        let idx = splice_index(r#"{"pos":100,"ref":"A","alt":"G","off":-2,"sig":"Pathogenic"}"#);
        let anns = vec![&idx];
        assert_eq!(
            same_splice_position_pathogenic(
                &[Consequence::MissenseVariant],
                &anns,
                Some("c.524G>A"),
                100,
                "A",
                "T"
            ),
            None
        );
    }

    #[test]
    fn test_ps1_splice_without_hgvs_still_matches_the_same_position() {
        // No HGVS means no offset for the variant being classified, so the
        // adjacent base cannot be reasoned about - but the same position is
        // the same dinucleotide whatever the offset is.
        let idx = splice_index(
            r#"{"pos":100,"ref":"A","alt":"G","off":-2,"sig":"Pathogenic"},
               {"pos":101,"ref":"G","alt":"A","off":-1,"sig":"Pathogenic"}"#,
        );
        assert_eq!(splice_ps1(Some(&idx), None, 100, "A", "T"), Some(true));
        assert_eq!(splice_ps1(Some(&idx), None, 102, "C", "T"), Some(false));
    }
}
