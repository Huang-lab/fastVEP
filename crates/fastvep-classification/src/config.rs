use serde::{Deserialize, Serialize};
use std::collections::HashMap;

use crate::types::EvidenceStrength;

/// Trio configuration for de novo and compound heterozygote analysis.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TrioConfig {
    /// Sample name of the proband (required)
    pub proband: String,
    /// Sample name of the mother (optional)
    pub mother: Option<String>,
    /// Sample name of the father (optional)
    pub father: Option<String>,
    /// Minimum read depth for reliable genotype call (default: 10)
    #[serde(default = "default_min_depth")]
    pub min_depth: u32,
    /// Minimum genotype quality for reliable genotype call (default: 20)
    #[serde(default = "default_min_gq")]
    pub min_gq: u32,
}

fn default_min_depth() -> u32 {
    10
}
fn default_min_gq() -> u32 {
    20
}

/// Configuration for ACMG-AMP classification thresholds and behavior.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AcmgConfig {
    // ── Frequency thresholds ──
    /// BA1: allele frequency threshold for standalone benign (default: 0.05)
    #[serde(default = "default_ba1")]
    pub ba1_af_threshold: f64,
    /// BS1: allele frequency above which a variant is more common than the
    /// disorder allows. Default 0.005.
    ///
    /// Richards 2015 words BS1 as "allele frequency is greater than expected
    /// for disorder", which is a per-disease quantity; 1 % was a placeholder
    /// for it, not a derivation. Whiffin 2017's maximum credible population
    /// allele frequency lands well below 1 % for essentially every Mendelian
    /// disorder, and published VCEP specifications cluster between 0.1 % and
    /// 0.5 %.
    ///
    /// Swept alone over a 1-in-10 sample of the benchmark, which is what picked
    /// the value:
    ///
    /// | bar | benign recall | false-pathogenic | false-benign |
    /// |---|---:|---:|---:|
    /// | 0.01 | 56.3 % | 4 | 1 |
    /// | **0.005** | **58.5 %** | **4** | **1** |
    /// | 0.001 | 61.9 % | 3 | 2 |
    /// | 0.0005 | 63.1 % | 3 | 7 |
    /// | 0.0001 | 65.9 % | 3 | 18 |
    ///
    /// On the sample 0.005 looked free. **It is not.** Confirming it on the
    /// full 673,660-variant set (v11 to v12) put benign recall at 65.3 % to
    /// 68.4 % and false-benign at 10 to 14: roughly 4,850 additional correct
    /// benign calls for 4 additional missed diagnoses, about 1,200 to 1. That
    /// is a good trade by the same exchange-rate reasoning used for the BS2
    /// prevalence bar, but it is a trade, and a sample of this size cannot
    /// resolve a four-count difference.
    ///
    /// Below 0.005 the cost climbs steeply - past 0.0005 false-benign roughly
    /// triples - so the remainder is left to a per-gene-disease threshold,
    /// which is what Richards' "expected for disorder" actually asks for.
    #[serde(default = "default_bs1")]
    pub bs1_af_threshold: f64,
    /// PM2: legacy single allele frequency threshold (default: 0.0001).
    ///
    /// Retained for backward compatibility with configs predating PR4. The
    /// classifier prefers the inheritance-aware fields below
    /// (`pm2_ad_af_threshold` / `pm2_ar_af_threshold`); this field is no
    /// longer consulted by the default code path but remains so that existing
    /// TOML configs continue to deserialize.
    #[serde(default = "default_pm2")]
    pub pm2_af_threshold: f64,
    /// PM2: AF threshold for autosomal-dominant or unknown-inheritance genes.
    /// Default 4e-5, i.e. "extremely rare" rather than strictly absent.
    ///
    /// Richards 2015 words PM2 as "absent from controls (or at extremely low
    /// frequency if recessive)", and fastVEP read that literally: `0.0`, so a
    /// dominant variant seen even once anywhere in gnomAD failed PM2. That
    /// reading does not survive the change in denominator. The 2015 text was
    /// written against ExAC's 60,706 exomes; gnomAD v4 carries 730,947 exomes
    /// plus 76,215 genomes. "Absent" from a cohort twelve times larger is a
    /// far stricter test than the authors were specifying, and it is the wrong
    /// test: a singleton among 800,000 people *is* the "not seen in the general
    /// population" that PM2 was asking about.
    ///
    /// ClinGen SVI moved the same way, downgrading PM2 to Supporting in 2020
    /// because it had been over-weighted, and the VCEP specifications written
    /// since then give dominant genes explicit non-zero bars, clustered between
    /// 1e-5 and 1e-4 (the Cardiomyopathy VCEP's MYH7 specification uses 4e-5;
    /// the PTEN VCEP uses 1e-5). 4e-5 sits inside that published range rather
    /// than being a number of our own.
    ///
    /// Measured over the ClinVar 2-star+ benchmark, sweeping this key alone:
    ///
    /// | bar | pathogenic recall | false-pathogenic | benign recall |
    /// |---|---:|---:|---:|
    /// | 0 (strict absence) | 37.8 % | 1 | 56.3 % |
    /// | 1e-6 | 46.2 % | 1 | 56.3 % |
    /// | 1e-5 | 54.0 % | 2 | 56.3 % |
    /// | **4e-5** | **56.8 %** | **2** | **56.3 %** |
    /// | 1e-4 | 57.5 % | 2 | 56.3 % |
    /// | 2e-4 | 57.6 % | 3 | 56.3 % |
    ///
    /// The curve flattens after 4e-5 - the last 0.8 pp costs a 2.5× loosening -
    /// and benign recall never moves, PM2 being pathogenic-direction only.
    /// Going further than 4e-5 would also take us past what any VCEP has
    /// published for a dominant gene, to buy a fraction of a point.
    #[serde(default = "default_pm2_ad")]
    pub pm2_ad_af_threshold: f64,
    /// PM2: AF threshold for autosomal-recessive genes per ClinGen SVI v1.0
    /// (Sept 2020). Default 0.00007 (0.007%).
    #[serde(default = "default_pm2_ar")]
    pub pm2_ar_af_threshold: f64,

    // ── REVEL thresholds (ClinGen SVI calibrated) ──
    /// PP3 supporting threshold (default: 0.644)
    #[serde(default = "default_pp3_supporting")]
    pub pp3_revel_supporting: f64,
    /// PP3 moderate threshold (default: 0.773)
    #[serde(default = "default_pp3_moderate")]
    pub pp3_revel_moderate: f64,
    /// PP3 strong threshold (default: 0.932)
    #[serde(default = "default_pp3_strong")]
    pub pp3_revel_strong: f64,
    /// BP4 supporting threshold (default: 0.290)
    #[serde(default = "default_bp4_supporting")]
    pub bp4_revel_supporting: f64,
    /// BP4 moderate threshold (default: 0.183)
    #[serde(default = "default_bp4_moderate")]
    pub bp4_revel_moderate: f64,
    /// BP4 strong threshold (default: 0.016)
    #[serde(default = "default_bp4_strong")]
    pub bp4_revel_strong: f64,
    /// BP4 very strong threshold (default: 0.003) — Pejaver 2022 endorses
    /// REVEL ≤ 0.003 as Very Strong benign evidence; only REVEL reaches this band.
    #[serde(default = "default_bp4_very_strong")]
    pub bp4_revel_very_strong: f64,

    // ── SpliceAI thresholds ──
    /// SpliceAI delta score threshold for PP3 pathogenic splice impact (default: 0.2).
    /// Per Walker 2023 SVI Splicing Subgroup, ≥ 0.2 yields PP3 at *Supporting* strength only.
    /// SpliceAI alone does not reach Strong — that requires experimental RNA assay (PVS1_RNA).
    #[serde(default = "default_spliceai_pathogenic")]
    pub spliceai_pathogenic: f64,
    /// SpliceAI delta score upper bound for BP4 benign splice impact (default: 0.1).
    /// Per Walker 2023 SVI Splicing Subgroup, ≤ 0.1 yields BP4 at Supporting strength.
    /// Scores between 0.1 and 0.2 are uninformative.
    #[serde(default = "default_spliceai_benign")]
    pub spliceai_benign: f64,

    // ── Conservation thresholds ──
    /// PhyloP threshold for conserved position (default: 2.0)
    #[serde(default = "default_phylop")]
    pub phylop_conserved: f64,

    // ── Gene constraint thresholds ──
    /// pLI threshold for LOF-intolerant gene (default: 0.9)
    #[serde(default = "default_pli")]
    pub pli_lof_intolerant: f64,
    /// LOEUF threshold for LOF-intolerant gene (default: 0.35)
    #[serde(default = "default_loeuf")]
    pub loeuf_lof_intolerant: f64,
    /// Missense Z-score threshold for PP2 (default: 3.09)
    #[serde(default = "default_misz")]
    pub pp2_misz_threshold: f64,

    // ── PM1 hotspot detection thresholds ──
    /// Window size (in amino acid positions) for hotspot detection (default: 5)
    #[serde(default = "default_pm1_window")]
    pub pm1_hotspot_window: u64,
    /// Minimum pathogenic variants in window to call hotspot (default: 3)
    #[serde(default = "default_pm1_threshold")]
    pub pm1_hotspot_min_pathogenic: u32,

    // ── ClinGen SVI behavior modifications ──
    /// Downgrade PM2 from Moderate to Supporting (ClinGen SVI recommendation)
    #[serde(default = "default_true")]
    pub pm2_downgrade_to_supporting: bool,
    /// When `input.gnomad` is `None` (no gnomAD record at all for the
    /// variant), treat the variant as absent from gnomAD and fire PM2.
    /// Per ClinGen SVI v1.0, "absent or extremely rare in population
    /// databases" is the PM2 trigger; if a variant is not in the loaded
    /// gnomAD `.osa`, the natural interpretation is that gnomAD never
    /// observed it (i.e. it IS absent). Default `true`.
    ///
    /// Set `false` to keep the strict-coverage stance (PM2 NotEvaluated
    /// when no record present) — useful when gnomAD data was loaded for
    /// only a subset of input regions and you want PM2 silenced outside
    /// that coverage.
    #[serde(default = "default_true")]
    pub pm2_absent_when_no_record: bool,
    /// Minimum allele count (AC) required to fire BS2 on an autosomal-
    /// dominant or X-linked-dominant gene from heterozygous gnomAD
    /// observations alone. Singleton / doubleton observations of a
    /// novel allele in a 100K-cohort are not sufficient evidence that
    /// the variant is tolerated in healthy adults - Richards 2015 BS2
    /// requires "observed in unaffected adult". Recessive and X-linked
    /// genes take the separate null-individual test instead, governed by
    /// [`Self::bs2_ar_min_hom`] and
    /// [`Self::bs2_hom_prevalence_threshold`]. Default `5` mirrors common
    /// ClinGen VCEP practice (e.g. Hereditary Cancer / Lynch).
    #[serde(default = "default_bs2_ad_min_ac")]
    pub bs2_ad_min_ac: u64,
    /// Absolute floor on homozygous observations before BS2 can fire on a
    /// recessive (or unknown-inheritance) gene. A single homozygote is not
    /// evidence of tolerance: 20 of the 52 BS2 misfires in the round-2
    /// medical-genetics review had exactly one, and 6 had two. Default 2.
    #[serde(default = "default_bs2_ar_min_hom")]
    pub bs2_ar_min_hom: u64,
    /// Maximum disease prevalence, as a fraction of individuals, that BS2's
    /// homozygote test is asked to rule out. BS2 fires only when the 95 %
    /// lower confidence bound on the homozygote frequency in gnomAD exceeds
    /// this value, i.e. when there are demonstrably more homozygotes than the
    /// disorder itself could account for.
    ///
    /// This is what makes the criterion scale with cohort size rather than
    /// with a bare count: one homozygote among 730 K individuals and one among
    /// 5 K are very different observations, and Richards 2015 BS2 asks for an
    /// observation "in a healthy adult", not for any observation at all.
    ///
    /// **Default 1e-3, chosen from measurement rather than convention.** A
    /// sweep of the full 673,660-variant ClinVar 2-star+ benchmark
    /// (`analysis/acmg_benchmark/scripts/sweep_acmg_thresholds.py`) gives:
    ///
    /// | bar | false-benign calls | correct benign calls |
    /// |---|---|---|
    /// | 1e-6 | 54 | 139,270 |
    /// | 1e-5 | 45 | 136,014 |
    /// | 1e-4 | 40 | 133,407 |
    /// | 1e-3 | 38 | 132,815 |
    ///
    /// The curve has no knee, so the data alone does not pick a value; the
    /// choice rests on what the parameter means and on which error is worse.
    /// It is a *maximum credible disease prevalence*, so it has to cover the
    /// most prevalent Mendelian conditions BS2 is applied to, not the typical
    /// one. Hearing loss, alpha-1 antitrypsin deficiency and familial
    /// Mediterranean fever in high-prevalence populations all sit near
    /// 1 in 1,000, and a 1e-5 bar is two orders of magnitude too tight for
    /// them - which is exactly the failure a medical geneticist raised in the
    /// round-2 review. A false-benign call is a missed diagnosis, whereas a
    /// lost benign call becomes a VUS and costs triage effort, so the
    /// asymmetry favours the safer bar. The step from 1e-4 to 1e-3 is also the
    /// cheapest on the curve, at 296 correct benign calls per false-benign
    /// avoided against 542 for the step before it.
    ///
    /// Raise or lower it freely: a lab willing to trade the other way sets a
    /// smaller value, and a gene-specific prevalence from a VCEP
    /// specification should override it per gene once that table lands.
    #[serde(default = "default_bs2_prevalence")]
    pub bs2_hom_prevalence_threshold: f64,
    /// Genes where population allele frequencies cannot be trusted because
    /// paralogues, pseudogenes or segmental duplications cause systematic
    /// mismapping (Mandelker et al. 2016, PMID 27228465). BA1, BS1, BS2 and
    /// PM2 are all reported NotEvaluated for these genes rather than fired on
    /// unreliable counts.
    #[serde(default = "default_homology_unreliable_genes")]
    pub homology_unreliable_genes: Vec<String>,
    /// Suppress BS1/BS2 when ClinVar (2 stars or better) labels the variant
    /// with a low-penetrance or risk-allele term. Such a variant is outside
    /// BS2's "full penetrance expected at an early age" precondition by
    /// definition, and its population frequency is expected to be high.
    ///
    /// This gate reads ClinVar, so it must be disabled when measuring
    /// concordance against ClinVar itself. Default `true`.
    #[serde(default = "default_true")]
    pub clinvar_low_penetrance_blocks_benign_frequency: bool,
    /// Treat gnomAD's `segdup` and `lcr` region flags as making a site's
    /// allele frequency unusable, in either direction.
    ///
    /// This is the per-site form of [`Self::homology_unreliable_genes`]: a
    /// segmental duplication is exactly the context in which reads from a
    /// paralogue mismap onto the gene of interest and inflate its apparent
    /// frequency. It is the more aggressive of the two, since it fires on
    /// individual sites inside otherwise well-behaved genes, so it is
    /// switchable. gnomAD's own FILTER verdict is not covered by this flag and
    /// is never ignored. Default `true`; a no-op against annotation databases
    /// built before the flags were extracted.
    #[serde(default = "default_true")]
    pub gnomad_region_flags_block_frequency: bool,
    /// Test BA1 and BS1 against gnomAD's filtering allele frequency (the 95 %
    /// CI lower bound, maximised over genetic-ancestry groups) rather than the
    /// population-maximum point estimate.
    ///
    /// This is the ClinGen/Whiffin 2017 recommendation: a point estimate makes
    /// a frequency measured from a few hundred alleles look as solid as one
    /// measured from hundreds of thousands. Default `true`; a no-op against
    /// annotation databases built before the FAF columns were extracted.
    #[serde(default = "default_true")]
    pub use_filtering_af: bool,
    /// Combine criteria with the ClinGen SVI point system (Tavtigian 2020)
    /// instead of the Richards 2015 combining table.
    ///
    /// The two are encodings of the same Bayesian model, and where they
    /// disagree the point system is the SVI's current recommendation and what
    /// VCEP specifications are written against. The disagreement that matters
    /// most is a lone PVS1: 8 points, inside the Likely Pathogenic band, but no
    /// row of Table 5 matches it, so the table returns Uncertain Significance.
    ///
    /// Measured on the ClinVar 2-star+ benchmark, run v10 put **2,319
    /// truth-pathogenic variants in VUS on `PVS1` and nothing else** - nonsense
    /// and frameshift variants in haploinsufficient disease genes.
    #[serde(default = "default_true")]
    pub use_point_system: bool,
    /// Require an established gene-disease relationship before the criteria
    /// that presuppose one - PVS1, PP2 and PM1 - may fire.
    ///
    /// Every pathogenic criterion is an argument that a variant disrupts a gene
    /// *in a way that causes a known disease*. PVS1 makes that explicit: Abou
    /// Tayoun 2018 opens by requiring that loss of function be an established
    /// mechanism for the gene, which presupposes the gene causes disease at
    /// all. PP2 and PM1 are the same argument at lower strength - "missense is
    /// how this gene causes disease", "this residue matters for the disease".
    /// Applied to a gene with no established relationship, all three assert a
    /// disease association that nobody has demonstrated.
    ///
    /// The gate reads the loaded gene-disease source (ClinGen Gene-Disease
    /// Validity, filtered to Definitive/Strong/Moderate; or OMIM `genemap2` for
    /// installations still on the legacy source). A gene absent from it takes
    /// the three criteria to NotEvaluated with
    /// `no_established_gene_disease_relationship`, which lands the variant in
    /// VUS - the honest answer when the gene itself is not established.
    ///
    /// **Degrades to a no-op when no gene-disease source is loaded**, since
    /// "absent from a database that was never opened" is not evidence. Default
    /// `true`.
    #[serde(default = "default_true")]
    pub require_gene_disease_validity: bool,
    /// Let a curated disease mechanism decide whether PVS1 applies, instead of
    /// only letting it enable PVS1.
    ///
    /// Abou Tayoun 2018 opens the PVS1 decision tree by asking whether loss of
    /// function is a known mechanism for the gene. fastVEP used to read
    /// [`GeneOverride::mechanism`] in one direction only - a `"LOF"` statement
    /// could switch PVS1 on for a gene that gnomAD constraint had not flagged,
    /// but a gain-of-function statement could not switch it off. With this
    /// enabled a mechanism that excludes loss of function (`"GOF"`,
    /// `"DOMINANT_NEGATIVE"`) takes PVS1 to NotApplicable.
    ///
    /// Constraint scores cannot substitute for this. A gain-of-function gene
    /// under strong purifying selection has a high pLI for the same reason a
    /// haploinsufficient one does, so pLI alone will happily carry PVS1 into a
    /// gene where a null allele is the harmless outcome - PCSK9 being the
    /// clearest case, where loss of function lowers LDL and is the mechanism
    /// of an approved drug class.
    ///
    /// Only genes carrying an explicit mechanism statement are affected;
    /// unknown mechanism is never read as gain of function. Default `true`.
    #[serde(default = "default_true")]
    pub mechanism_gates_pvs1: bool,
    /// Use the measured 50-nucleotide rule for PVS1's NMD prediction instead of
    /// the last-exon proxy. Default `false`, and the reason is a measurement
    /// rather than a preference.
    ///
    /// Abou Tayoun 2018 predicts nonsense-mediated decay by asking whether the
    /// premature termination codon sits more than 50 nt upstream of the final
    /// exon-exon junction. fastVEP can now answer that exactly
    /// ([`Transcript::escapes_nmd`]); before, it approximated it with "is the
    /// variant in the last exon?", which misses the escape window at the 3' end
    /// of the penultimate exon and so keeps those PTCs at full Very Strong.
    ///
    /// Turning the exact rule on is the guideline-faithful setting, and it is
    /// off by default because of what it costs against ClinVar. On a
    /// 4,000-variant sample of PVS1-carrying 2-star records it moved 58 calls
    /// from Likely Pathogenic to Uncertain: **54 of them ClinVar calls
    /// Pathogenic**, 4 Uncertain, and none Benign. Extrapolated to the full
    /// PVS1 population that is roughly 825 correct pathogenic calls given up to
    /// correct a single false-pathogenic one.
    ///
    /// The trade is not a defect in either the rule or the tree - both give the
    /// right answer. A PTC 20 nt before the last junction really does escape
    /// decay, and Abou Tayoun really does grade an NMD-escaping truncation in a
    /// critical region at PVS1_Strong, which is 4 points and Uncertain on its
    /// own. ClinVar's curators call those variants Pathogenic because they also
    /// had segregation, case and functional evidence that no annotator can
    /// compute. Enabling this makes fastVEP more faithful to the guideline and
    /// less concordant with ClinVar at the same time, so the choice is the
    /// user's rather than a default.
    #[serde(default)]
    pub pvs1_nmd_50nt_rule: bool,
    /// Optional ceiling on the strength PP3 may reach from computational
    /// evidence alone. `None` (the default) means uncapped.
    ///
    /// Pejaver 2022 is a ClinGen SVI product and explicitly calibrates
    /// REVEL >= 0.932 to Strong, so capping by default would put fastVEP
    /// outside the guideline. The round-2 medical-genetics review holds the
    /// stricter view that a predictor should not reach Strong on its own, and
    /// several VCEP specifications agree. This knob lets a lab following that
    /// convention configure it (`pp3_max_strength = "Moderate"`) instead of
    /// patching the classifier.
    #[serde(default)]
    pub pp3_max_strength: Option<EvidenceStrength>,
    /// Exclude the variant being classified from the ClinVar-derived evidence
    /// that PS1 and PM1 read. PS1 means "same amino acid change as a
    /// *previously established* pathogenic variant, regardless of the
    /// nucleotide change"; matching a variant against its own ClinVar record
    /// is circular and inflates any ClinVar-based benchmark. When the variant
    /// itself is ClinVar pathogenic/likely-pathogenic, its own entry is
    /// discounted from the protein-index match count. Default `true`.
    #[serde(default = "default_true")]
    pub exclude_self_from_clinvar_evidence: bool,
    /// Maximum number of pathogenic/likely-pathogenic missense variants a gene
    /// may have in the ClinVar protein index before BP1 is ruled out. At or
    /// above this count, missense is an established disease mechanism for the
    /// gene and BP1 cannot apply. Default 3.
    #[serde(default = "default_bp1_max_pathogenic_missense")]
    pub bp1_max_pathogenic_missense: u32,
    /// Far boundary, in nucleotides from the nearest exon edge, of BP7's
    /// intronic extension. Default 300.
    ///
    /// Both ends of the range are meaningful, so no sentinel is needed: a value
    /// below Walker's near boundary of 7 turns the intronic extension off
    /// entirely, and one larger than any human intron restores it unbounded.
    ///
    /// Walker 2023 gives BP7's intronic extension a *near* boundary (donor
    /// `+7`, acceptor `-21`) and states no far one. That is consistent for the
    /// SVI, whose recommendations are written for a curator applying BP7 to a
    /// single variant in a gene they know. An automated pipeline applies it to
    /// every intronic position in the genome, and out in the deep intron the
    /// only evidence BP7 has left is a SpliceAI score in the regime where
    /// SpliceAI is weakest: pseudoexon activation, which is also the dominant
    /// mechanism by which a deep-intronic variant is pathogenic at all.
    ///
    /// So the far boundary is set by measurement. Over the full ClinVar
    /// 2-star-or-better set, the extension's exchange rate - correct benign
    /// firings per missed diagnosis - collapses with intron depth:
    ///
    /// | bound | correct benign calls | missed diagnoses | marginal cost |
    /// |---|---:|---:|---|
    /// | 50 | 176,132 | 25 | 248 correct benign per diagnosis recovered |
    /// | 100 | 177,122 | 29 | 204 |
    /// | 200 | 177,735 | 32 | 155 |
    /// | **300** | **177,890** | **33** | **40** |
    /// | 500 | 177,930 | 34 | 26 |
    /// | 1000 | 177,981 | 36 | 31 |
    /// | unbounded | 178,420 | 50 | - |
    ///
    /// 300 is the knee. Bounding an unbounded extension at 300 recovers 17 of
    /// its 50 missed diagnoses for 530 correct benign calls, about 31 each;
    /// tightening further costs 155 and then 204 and then 248 per diagnosis.
    /// For scale, the BS2 prevalence bar this classifier already ships trades
    /// ~296 correct benign calls per missed diagnosis avoided, so everything
    /// down to 300 is cheap by the project's own standard and everything below
    /// it is not. The sweep is reproducible with `sweep_acmg_thresholds.py
    /// --sweep "bp7_max_intron_offset=..."` and is recorded in `docs/ACMG.md`.
    ///
    /// This bounds only the Walker 2023 intronic extension. A synonymous
    /// variant is exonic and is unaffected.
    #[serde(default = "default_bp7_max_intron_offset")]
    pub bp7_max_intron_offset: u64,
    /// How many benign or likely-benign ClinVar missense variants may sit in
    /// PM1's window before the window stops counting as a hotspot. Default 0,
    /// which is Richards 2015 read literally: PM1 wants a mutational hot spot
    /// or critical domain "without benign variation".
    ///
    /// The knob exists because "without benign variation" is a judgement in
    /// practice - a single benign submission at the edge of a well-established
    /// active site does not dissolve the domain - and a lab following a VCEP
    /// specification that tolerates one or two can say so here. The default
    /// stays at the literal reading.
    ///
    /// Has no effect against a `clinvar_protein.oga` built before benign
    /// assertions were indexed; see [`ClinvarProteinData::benign_indexed`].
    #[serde(default)]
    pub pm1_max_benign_in_window: u32,
    /// Enable PP5/BP6 criteria (disabled by default per ClinGen SVI)
    #[serde(default)]
    pub use_pp5_bp6: bool,
    /// Opt back into the legacy PS4 proxy: ClinVar pathogenic with ≥3 review
    /// stars → PS4. ClinGen SVI considers this proxy invalid (true PS4
    /// requires case-control statistics), so it is disabled by default. Set
    /// `true` only for backward-comparable benchmarks.
    #[serde(default)]
    pub use_clinvar_stars_as_ps4_proxy: bool,

    /// Variants exempted from BA1 per the ClinGen SVI updated recommendation
    /// (Ghosh et al. 2018, Hum Mutat). These are well-known high-AF variants
    /// whose pathogenicity is established despite exceeding the 5% AF threshold.
    /// Defaults to the original 9-variant list; users may add VCEP-specific
    /// entries via TOML.
    #[serde(default = "default_ba1_exceptions")]
    pub ba1_exceptions: Vec<Ba1Exception>,

    /// Minimum allele number (AN) required for BA1 / BS1 to fire (gnomAD v4
    /// guidance, ClinGen SVI March 2024). With v4's massive expansion (807k
    /// exomes, 76k genomes), the overall AN should be ≥ 2000 before a
    /// frequency-based call is reliable. When the AN drops below this
    /// threshold, BA1/BS1 are marked NotEvaluated rather than firing on
    /// noisy frequency estimates. Default 2000.
    #[serde(default = "default_min_an")]
    pub min_an_for_frequency_criteria: u64,

    /// Curated disease mechanism per gene, consulted by PVS1 when the gene has
    /// no entry in [`Self::gene_overrides`].
    ///
    /// Values are `"LOF"`, `"GOF"`, `"DOMINANT_NEGATIVE"` or `"LOF_and_GOF"`,
    /// matched case-insensitively. This is shipped data rather than user
    /// policy, which is why it is a separate map: a TOML that sets one
    /// `[gene_overrides.X]` block would otherwise silently replace the whole
    /// curated table. `gene_overrides` still wins per gene where both name the
    /// same one.
    #[serde(default = "default_gene_mechanisms")]
    pub gene_mechanisms: HashMap<String, String>,

    // ── Gene-specific overrides ──
    #[serde(default)]
    pub gene_overrides: HashMap<String, GeneOverride>,

    // ── Trio configuration ──
    /// Trio configuration for de novo and compound heterozygote analysis.
    #[serde(default)]
    pub trio: Option<TrioConfig>,
}

/// One entry on the BA1 exception list. A variant is exempt from BA1 when it
/// matches an entry, regardless of allele frequency.
///
/// Matching is by genomic coordinate first and HGVS second, and the order
/// matters. An HGVS `c.` token is only meaningful against the transcript it was
/// written for: Ghosh 2018 lists BTD `c.1330G>C` on NM_000060.4, and fastVEP
/// reports the same variant as `c.1270G>C` on ENST00000643237.3. The same
/// variant can also be spelled two valid ways - ACAD9's `c.-44_-41dupTAAG` is
/// fastVEP's `c.-45_-44insTAAG`. A coordinate has neither problem.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Ba1Exception {
    pub gene: String,
    /// HGVS c. notation, e.g. "c.845G>A".
    ///
    /// Compared case-insensitively, and against the part of the call's HGVS
    /// after the last `:` - the pipeline emits `ENST00000357618.10:c.845G>A`,
    /// so a whole-string comparison against a bare `c.` token never matches.
    pub hgvs_c: String,
    /// Optional human-readable reason - surfaced in the criterion `summary`.
    #[serde(default)]
    pub reason: Option<String>,
    /// Chromosome, without a `chr` prefix. Set together with the other three
    /// coordinate fields or not at all.
    #[serde(default)]
    pub chrom: Option<String>,
    /// 1-based VCF position on [`Self::assembly`].
    #[serde(default)]
    pub pos: Option<u64>,
    /// VCF-form reference allele.
    #[serde(default, rename = "ref")]
    pub ref_allele: Option<String>,
    /// VCF-form alternate allele.
    #[serde(default)]
    pub alt: Option<String>,
    /// Assembly the coordinate is on. Defaults to GRCh38; an entry on another
    /// build is matched by HGVS only, because comparing its position against a
    /// GRCh38 call would silently exempt the wrong variant.
    #[serde(default = "default_exception_assembly")]
    pub assembly: String,
    /// ClinVar variation ID, for provenance. Not used for matching - fastVEP's
    /// `clinvar.osa` does not carry the ID - but it is what makes an entry
    /// checkable against the source table.
    #[serde(default)]
    pub clinvar_id: Option<String>,
}

fn default_exception_assembly() -> String {
    "GRCh38".to_string()
}

/// The part of an HGVS expression after the last `:`, which is the `c.` token
/// itself. Returns the whole string when there is no prefix.
fn hgvs_suffix(hgvs: &str) -> &str {
    hgvs.rsplit(':').next().unwrap_or(hgvs)
}

impl Ba1Exception {
    /// Whether this entry exempts the variant described by the arguments.
    ///
    /// `coordinates` is `(chrom, pos, ref, alt)` in VCF form on GRCh38, or
    /// `None` when the caller could not supply one.
    pub fn matches(
        &self,
        gene: Option<&str>,
        hgvs_c: Option<&str>,
        coordinates: Option<(&str, u64, &str, &str)>,
    ) -> bool {
        if let (Some((chrom, pos, ref_allele, alt)), Some(e_chrom), Some(e_pos), Some(e_ref), Some(e_alt)) =
            (coordinates, self.chrom.as_deref(), self.pos, self.ref_allele.as_deref(), self.alt.as_deref())
        {
            if self.assembly.eq_ignore_ascii_case("GRCh38")
                && e_chrom.trim_start_matches("chr").eq_ignore_ascii_case(chrom.trim_start_matches("chr"))
                && e_pos == pos
                && e_ref.eq_ignore_ascii_case(ref_allele)
                && e_alt.eq_ignore_ascii_case(alt)
            {
                return true;
            }
        }
        match (gene, hgvs_c) {
            (Some(g), Some(h)) => {
                self.gene.eq_ignore_ascii_case(g)
                    && hgvs_suffix(&self.hgvs_c).eq_ignore_ascii_case(hgvs_suffix(h))
            }
            _ => false,
        }
    }
}

/// Gene-specific overrides for ACMG-AMP criteria.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GeneOverride {
    /// Disease mechanism: `"LOF"`, `"GOF"`, `"DOMINANT_NEGATIVE"` or
    /// `"LOF_and_GOF"`, matched case-insensitively. Overrides the curated
    /// [`AcmgConfig::gene_mechanisms`] entry for this gene. A mechanism that
    /// excludes loss of function takes PVS1 to NotApplicable; see
    /// [`AcmgConfig::mechanism_gates_pvs1`].
    pub mechanism: Option<String>,
    /// Override the BA1 standalone-benign allele frequency threshold.
    ///
    /// Where a ClinGen VCEP has published a bar for this gene it outranks the
    /// global default, which is measured across all genes and cannot know the
    /// disorder's prevalence, penetrance or allelic heterogeneity. Published
    /// bars run in both directions: ABCA4's is 0.163, three times looser than
    /// the 0.05 default, and CDKL5's is 8.3e-5, six hundred times tighter.
    pub ba1_af_threshold: Option<f64>,
    /// Override BS1 allele frequency threshold
    pub bs1_af_threshold: Option<f64>,
    /// Override PM2 allele frequency threshold
    pub pm2_af_threshold: Option<f64>,
    /// Criteria codes to disable for this gene
    #[serde(default)]
    pub disabled_criteria: Vec<String>,
    /// Criteria strength overrides (code -> new strength)
    #[serde(default)]
    pub strength_overrides: HashMap<String, EvidenceStrength>,
    /// Per-disorder thresholds for genes associated with multiple disorders
    /// (ClinGen SVI guidance July 2025). The classifier consumes whichever
    /// disorder context is active for the call; in the absence of explicit
    /// disorder selection, this scaffold is currently informational only —
    /// the active disorder selection mechanism is part of a follow-up PR.
    #[serde(default)]
    pub disorders: HashMap<String, DisorderOverride>,
}

/// Per-disorder override values within a multi-disorder gene.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct DisorderOverride {
    /// Inheritance for this disorder ("AD", "AR", or "AD_AR").
    pub inheritance: Option<String>,
    /// Override BA1 AF threshold for this disorder.
    #[serde(default)]
    pub ba1_af_threshold: Option<f64>,
    /// Override BS1 AF threshold for this disorder.
    pub bs1_af_threshold: Option<f64>,
    /// Override PM2 AF threshold for this disorder.
    pub pm2_af_threshold: Option<f64>,
}

impl Default for AcmgConfig {
    fn default() -> Self {
        Self {
            ba1_af_threshold: 0.05,
            bs1_af_threshold: default_bs1(),
            pm2_af_threshold: 0.0001,
            pm2_ad_af_threshold: default_pm2_ad(),
            pm2_ar_af_threshold: 0.00007,
            pp3_revel_supporting: 0.644,
            pp3_revel_moderate: 0.773,
            pp3_revel_strong: 0.932,
            bp4_revel_supporting: 0.290,
            bp4_revel_moderate: 0.183,
            bp4_revel_strong: 0.016,
            bp4_revel_very_strong: 0.003,
            spliceai_pathogenic: 0.2,
            spliceai_benign: 0.1,
            phylop_conserved: 2.0,
            pli_lof_intolerant: 0.9,
            loeuf_lof_intolerant: 0.35,
            pp2_misz_threshold: 3.09,
            pm1_hotspot_window: 5,
            pm1_hotspot_min_pathogenic: 3,
            pm2_downgrade_to_supporting: true,
            pm2_absent_when_no_record: true,
            bs2_ad_min_ac: 5,
            bs2_ar_min_hom: 2,
            bs2_hom_prevalence_threshold: 1e-3,
            homology_unreliable_genes: default_homology_unreliable_genes(),
            clinvar_low_penetrance_blocks_benign_frequency: true,
            gnomad_region_flags_block_frequency: true,
            use_filtering_af: true,
            use_point_system: true,
            require_gene_disease_validity: true,
            mechanism_gates_pvs1: true,
            pvs1_nmd_50nt_rule: false,
            gene_mechanisms: default_gene_mechanisms(),
            pp3_max_strength: None,
            exclude_self_from_clinvar_evidence: true,
            bp1_max_pathogenic_missense: 3,
            bp7_max_intron_offset: default_bp7_max_intron_offset(),
            pm1_max_benign_in_window: 0,
            use_pp5_bp6: false,
            ba1_exceptions: default_ba1_exceptions(),
            use_clinvar_stars_as_ps4_proxy: false,
            min_an_for_frequency_criteria: 2000,
            gene_overrides: HashMap::new(),
            trio: None,
        }
    }
}

impl AcmgConfig {
    /// Load configuration from a TOML file.
    pub fn from_toml_file(path: &str) -> anyhow::Result<Self> {
        let content = std::fs::read_to_string(path)?;
        let config: Self = toml::from_str(&content)?;
        Ok(config)
    }

    /// Get gene-specific override, if any.
    pub fn gene_override(&self, gene: &str) -> Option<&GeneOverride> {
        self.gene_overrides.get(gene)
    }

    /// Check if a criterion is disabled for a given gene.
    pub fn is_criterion_disabled(&self, gene: &str, criterion_code: &str) -> bool {
        self.gene_overrides
            .get(gene)
            .is_some_and(|o| o.disabled_criteria.iter().any(|c| c == criterion_code))
    }

    /// The disease mechanism to use for a gene: an explicit
    /// [`GeneOverride::mechanism`] if the user set one, otherwise the curated
    /// [`Self::gene_mechanisms`] entry, otherwise `None` for "not curated".
    pub fn effective_mechanism(&self, gene: Option<&str>) -> Option<&str> {
        let gene = gene?;
        self.gene_overrides
            .get(gene)
            .and_then(|o| o.mechanism.as_deref())
            .or_else(|| self.gene_mechanisms.get(gene).map(String::as_str))
    }

    /// True when the gene's population frequencies are unreliable because of
    /// paralogue / pseudogene / segmental-duplication mismapping.
    pub fn is_homology_unreliable(&self, gene: Option<&str>) -> bool {
        match gene {
            Some(g) => self
                .homology_unreliable_genes
                .iter()
                .any(|h| h.eq_ignore_ascii_case(g)),
            None => false,
        }
    }

    /// Get effective BA1 threshold for a gene (gene-specific or default).
    pub fn effective_ba1_threshold(&self, gene: Option<&str>) -> f64 {
        gene.and_then(|g| {
            self.gene_overrides
                .get(g)
                .and_then(|o| o.ba1_af_threshold)
        })
        .unwrap_or(self.ba1_af_threshold)
    }

    /// Get effective BS1 threshold for a gene (gene-specific or default).
    pub fn effective_bs1_threshold(&self, gene: Option<&str>) -> f64 {
        gene.and_then(|g| {
            self.gene_overrides
                .get(g)
                .and_then(|o| o.bs1_af_threshold)
        })
        .unwrap_or(self.bs1_af_threshold)
    }

}

// Default value functions for serde
fn default_ba1() -> f64 { 0.05 }
fn default_bs1() -> f64 { 0.005 }
fn default_pm2() -> f64 { 0.0001 }
fn default_pm2_ad() -> f64 { 0.00004 }
fn default_pm2_ar() -> f64 { 0.00007 }
fn default_bs2_ad_min_ac() -> u64 { 5 }
fn default_bs2_ar_min_hom() -> u64 { 2 }
fn default_bs2_prevalence() -> f64 { 1e-3 }
fn default_bp1_max_pathogenic_missense() -> u32 { 3 }
fn default_bp7_max_intron_offset() -> u64 { 300 }

/// Genes whose gnomAD frequencies are unreliable because of paralogue,
/// pseudogene or segmental-duplication mismapping (Mandelker et al. 2016,
/// PMID 27228465, "Navigating highly homologous genes in a molecular
/// diagnostic setting"). The reviewer flagged CYP21A2 and STRC directly in
/// the round-2 review and cited this reference; HBA1/HBA2 rows in the same
/// file have the same problem. Users may replace the list via TOML.
fn default_homology_unreliable_genes() -> Vec<String> {
    [
        // Reviewer-flagged in the round-2 discordance review.
        "CYP21A2", "STRC", "HBA1", "HBA2",
        // Canonical members of the Mandelker 2016 set that also appear in
        // clinical panels and carry a pseudogene or near-identical paralogue.
        "SMN1", "SMN2", "PMS2", "NEB", "OTOA", "GBA", "IKBKG", "CFC1", "NCF1",
        "TTN", "CYP11B1", "CYP11B2", "HBB", "HBD", "SBDS", "FCGR3A", "FCGR3B",
        "CR1", "C4A", "C4B", "TUBB8", "MOCS1", "OPN1LW", "OPN1MW", "GTF2I",
    ]
    .iter()
    .map(|s| s.to_string())
    .collect()
}
/// Curated disease mechanisms, for the PVS1 question "is loss of function a
/// known mechanism for this gene?" (Abou Tayoun 2018).
///
/// Deliberately small. Every entry is a gene where the mechanism is settled
/// enough that a ClinGen VCEP specification, or an approved therapy built on
/// it, depends on the distinction. Genes whose mechanism is merely unstudied
/// are absent, and absence never blocks PVS1.
///
/// The `LOF_and_GOF` entries are inert - they leave PVS1 exactly as it was -
/// and are listed anyway, because each is a gene someone will reasonably
/// propose adding as GOF. Recording why they are not is cheaper than
/// re-litigating them.
fn default_gene_mechanisms() -> HashMap<String, String> {
    [
        // ── Gain of function: a null allele is not the disease mechanism ──
        // LoF lowers LDL and is protective; it is the mechanism of the
        // PCSK9-inhibitor drug class. Flagged in the round-2 review.
        ("PCSK9", "GOF"),
        // Constitutive MDA5 signalling → Aicardi-Goutières 7 / Singleton-Merten.
        ("IFIH1", "GOF"),
        // Recurrent de novo missense (p.Arg87) → DEE65 by altered function.
        ("CYFIP2", "GOF"),
        // Cardiac conduction disease with dilated cardiomyopathy.
        ("TNNI3K", "GOF"),
        // Inflammasome activation → cryopyrin-associated periodic syndromes.
        ("NLRP3", "GOF"),
        // ALS by toxic gain of function; heterozygous null carriers are healthy.
        ("SOD1", "GOF"),
        // Prion disease by misfolding; heterozygous nulls are healthy.
        ("PRNP", "GOF"),
        // RASopathies: activating missense in the RAS-MAPK cascade. The ClinGen
        // RASopathy VCEP does not apply PVS1 to these genes.
        ("HRAS", "GOF"),
        ("KRAS", "GOF"),
        ("NRAS", "GOF"),
        ("BRAF", "GOF"),
        ("RAF1", "GOF"),
        ("MAP2K1", "GOF"),
        ("MAP2K2", "GOF"),
        ("SOS1", "GOF"),
        ("RIT1", "GOF"),
        ("SHOC2", "GOF"),
        ("PTPN11", "GOF"),
        // Primary open-angle glaucoma is caused by mutant myocilin misfolding
        // and accumulating; whole-gene deletions and null alleles do not cause
        // it. Surfaced by run v11 as a truncating variant ClinVar calls benign
        // that had collected full PVS1.
        ("MYOC", "GOF"),
        // SCA17 is a CAG/polyglutamine expansion - a gain of function in the
        // repeat, not loss of the protein. Same run, same shape.
        ("TBP", "GOF"),
        // ── Dominant negative: also not loss of function ──
        // ClinGen MYH7 specification (Kelly 2018): PVS1 is not applicable,
        // truncating variants are not established as causing HCM.
        ("MYH7", "DOMINANT_NEGATIVE"),
        // ── Both mechanisms: PVS1 still applies, listed for the record ──
        // Malignant hyperthermia is GoF; the congenital myopathies include
        // recessive LoF.
        ("RYR1", "LOF_and_GOF"),
        // Brugada is LoF, LQT3 is GoF.
        ("SCN5A", "LOF_and_GOF"),
        // LQT1 is LoF, short-QT/atrial fibrillation is GoF.
        ("KCNQ1", "LOF_and_GOF"),
        // McCune-Albright is GoF, pseudohypoparathyroidism is LoF.
        ("GNAS", "LOF_and_GOF"),
        // The achondroplasia group is GoF, CATSHL syndrome is LoF.
        ("FGFR3", "LOF_and_GOF"),
    ]
    .iter()
    .map(|(g, m)| (g.to_string(), m.to_string()))
    .collect()
}

fn default_pp3_supporting() -> f64 { 0.644 }
fn default_pp3_moderate() -> f64 { 0.773 }
fn default_pp3_strong() -> f64 { 0.932 }
fn default_bp4_supporting() -> f64 { 0.290 }
fn default_bp4_moderate() -> f64 { 0.183 }
fn default_bp4_strong() -> f64 { 0.016 }
fn default_bp4_very_strong() -> f64 { 0.003 }
fn default_spliceai_pathogenic() -> f64 { 0.2 }
fn default_spliceai_benign() -> f64 { 0.1 }
fn default_phylop() -> f64 { 2.0 }
fn default_pli() -> f64 { 0.9 }
fn default_loeuf() -> f64 { 0.35 }
fn default_misz() -> f64 { 3.09 }
fn default_pm1_window() -> u64 { 5 }
fn default_pm1_threshold() -> u32 { 3 }
fn default_true() -> bool { true }
fn default_min_an() -> u64 { 2000 }

/// Default BA1 exception list (Ghosh et al. 2018, Hum Mutat — 9 variants).
/// The ClinGen SVI BA1 exception list (Ghosh 2018, PMID 30311378, Table 1).
///
/// Nine variants that sit above 5 % in at least one population and must not be
/// called standalone-benign on that basis. Some are outright pathogenic (HFE
/// C282Y, GJB2 p.Val37Ile, BTD p.Asp444His); the rest are common in one
/// ancestry group for reasons that have nothing to do with being harmless.
///
/// Coordinates are GRCh38, lifted from the GRCh37 positions in the paper's
/// Table 1 via the ClinVar variation IDs it also lists, which are carried here
/// so any entry can be checked against the source. Matching prefers the
/// coordinate precisely because two of these nine cannot be matched on their
/// published HGVS: BTD's `c.1330G>C` is NM_000060.4 numbering and fastVEP
/// reports `c.1270G>C` on the Ensembl canonical transcript, and ACAD9's
/// `c.-44_-41dupTAAG` is the same variant fastVEP spells `c.-45_-44insTAAG`.
fn default_ba1_exceptions() -> Vec<Ba1Exception> {
    let mk = |gene: &str, hgvs: &str, clinvar_id: &str, chrom: &str, pos: u64, r: &str, a: &str, reason: &str| {
        Ba1Exception {
            gene: gene.to_string(),
            hgvs_c: hgvs.to_string(),
            reason: Some(reason.to_string()),
            chrom: Some(chrom.to_string()),
            pos: Some(pos),
            ref_allele: Some(r.to_string()),
            alt: Some(a.to_string()),
            assembly: default_exception_assembly(),
            clinvar_id: Some(clinvar_id.to_string()),
        }
    };
    vec![
        mk("ACAD9", "c.-44_-41dupTAAG", "1018", "3", 128_879_647, "C", "CTAAG",
           "Ghosh 2018 BA1 exception (VUS); 12.6 % in African/African American"),
        mk("GJB2", "c.109G>A", "17023", "13", 20_189_473, "C", "T",
           "Ghosh 2018 BA1 exception (Pathogenic, p.Val37Ile); 7.2 % in East Asian, autosomal recessive deafness"),
        mk("HFE", "c.187C>G", "10", "6", 26_090_951, "C", "G",
           "Ghosh 2018 BA1 exception (Pathogenic, p.His63Asp); 13.7 % in non-Finnish European, hereditary hemochromatosis"),
        mk("HFE", "c.845G>A", "9", "6", 26_092_913, "G", "A",
           "Ghosh 2018 BA1 exception (Pathogenic, p.Cys282Tyr); 5.1 % in non-Finnish European, hereditary hemochromatosis"),
        mk("MEFV", "c.1105C>T", "2551", "16", 3_249_586, "G", "A",
           "Ghosh 2018 BA1 exception (VUS, p.Pro369Ser); 7.2 % in East Asian, familial Mediterranean fever"),
        mk("MEFV", "c.1223G>A", "2552", "16", 3_249_468, "C", "T",
           "Ghosh 2018 BA1 exception (VUS, p.Arg408Gln); 5.4 % in East Asian, familial Mediterranean fever"),
        mk("PIBF1", "c.1214G>A", "217689", "13", 72_835_359, "G", "A",
           "Ghosh 2018 BA1 exception (VUS, p.Arg405Gln); 9.9 % in Latino, Joubert syndrome"),
        mk("ACADS", "c.511C>T", "3830", "12", 120_737_875, "C", "T",
           "Ghosh 2018 BA1 exception (VUS, p.Arg171Trp); 6.6 % in Finnish only, SCAD deficiency"),
        mk("BTD", "c.1330G>C", "1900", "3", 15_645_186, "G", "C",
           "Ghosh 2018 BA1 exception (Pathogenic, p.Asp444His); 5.4 % in Finnish only, biotinidase deficiency"),
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default_config() {
        let cfg = AcmgConfig::default();
        assert!((cfg.ba1_af_threshold - 0.05).abs() < 1e-10);
        assert!(cfg.pm2_downgrade_to_supporting);
        assert!(!cfg.use_pp5_bp6);
    }

    #[test]
    fn test_gene_override() {
        let mut cfg = AcmgConfig::default();
        cfg.gene_overrides.insert(
            "BRCA1".to_string(),
            GeneOverride {
                mechanism: Some("LOF".to_string()),
                ba1_af_threshold: None,
                bs1_af_threshold: Some(0.001),
                pm2_af_threshold: None,
                disabled_criteria: vec![],
                strength_overrides: HashMap::new(),
                disorders: HashMap::new(),
            },
        );
        assert_eq!(cfg.effective_bs1_threshold(Some("BRCA1")), 0.001);
        // Genes without an override fall back to the measured default.
        assert_eq!(cfg.effective_bs1_threshold(Some("TP53")), 0.005);
        assert_eq!(cfg.effective_bs1_threshold(None), 0.005);
    }

    #[test]
    fn test_vcep_frequency_bars_load_from_toml_in_both_directions() {
        // Both real, both from the ClinGen CSpec Registry: ABCA4's published
        // BA1 is three times *looser* than fastVEP's global default and CDKL5's
        // is six hundred times tighter. A per-gene BA1 has to express both, and
        // an unknown key here would be silently dropped by serde - which is why
        // this asserts the values took effect rather than that the file parsed.
        let toml = r#"
[gene_overrides.ABCA4]
ba1_af_threshold = 0.163
bs1_af_threshold = 0.0163
pm2_af_threshold = 0.0001

[gene_overrides.CDKL5]
ba1_af_threshold = 8.3e-5
"#;
        let cfg: AcmgConfig = toml::from_str(toml).expect("VCEP threshold TOML must parse");
        assert_eq!(cfg.effective_ba1_threshold(Some("ABCA4")), 0.163);
        assert_eq!(cfg.effective_bs1_threshold(Some("ABCA4")), 0.0163);
        assert_eq!(cfg.effective_ba1_threshold(Some("CDKL5")), 8.3e-5);
        // CDKL5 sets no BS1, so it keeps the global bar rather than inheriting
        // its own BA1.
        assert_eq!(cfg.effective_bs1_threshold(Some("CDKL5")), 0.005);
        assert_eq!(cfg.effective_ba1_threshold(Some("TP53")), 0.05);
    }

    #[test]
    fn test_the_shipped_vcep_threshold_table_parses() {
        // The generated table is the artefact this whole path exists for, and
        // a TOML that no longer deserializes would otherwise only be noticed at
        // a user's command line.
        let path = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../analysis/acmg_benchmark/data/vcep_thresholds.toml"
        );
        let Ok(text) = std::fs::read_to_string(path) else {
            // Not present in a packaged crate; nothing to check.
            return;
        };
        let cfg: AcmgConfig =
            toml::from_str(&text).expect("generated VCEP threshold table must parse");
        assert!(
            cfg.gene_overrides.len() > 50,
            "expected the full table, got {} genes",
            cfg.gene_overrides.len()
        );
        for (gene, over) in &cfg.gene_overrides {
            for (label, value) in [
                ("BA1", over.ba1_af_threshold),
                ("BS1", over.bs1_af_threshold),
                ("PM2", over.pm2_af_threshold),
            ] {
                if let Some(v) = value {
                    assert!(
                        (0.0..=0.5).contains(&v),
                        "{gene} {label} = {v} is not an allele frequency"
                    );
                }
            }
            if let (Some(ba1), Some(bs1)) = (over.ba1_af_threshold, over.bs1_af_threshold) {
                assert!(ba1 >= bs1, "{gene}: BA1 {ba1} below BS1 {bs1}");
            }
            if let (Some(bs1), Some(pm2)) = (over.bs1_af_threshold, over.pm2_af_threshold) {
                assert!(bs1 >= pm2, "{gene}: BS1 {bs1} below PM2 {pm2}");
            }
        }
    }

    #[test]
    fn test_curated_mechanisms_ship_by_default() {
        let cfg = AcmgConfig::default();
        assert_eq!(cfg.effective_mechanism(Some("PCSK9")), Some("GOF"));
        assert_eq!(cfg.effective_mechanism(Some("MYH7")), Some("DOMINANT_NEGATIVE"));
        assert_eq!(cfg.effective_mechanism(Some("RYR1")), Some("LOF_and_GOF"));
        assert_eq!(cfg.effective_mechanism(Some("BRCA1")), None);
        assert_eq!(cfg.effective_mechanism(None), None);
        assert!(cfg.require_gene_disease_validity);
        assert!(cfg.mechanism_gates_pvs1);
    }

    #[test]
    fn test_gene_overrides_mechanism_wins_over_the_curated_table() {
        let cfg: AcmgConfig = toml::from_str(
            "[gene_overrides.PCSK9]\nmechanism = \"LOF\"\n",
        )
        .unwrap();
        assert_eq!(cfg.effective_mechanism(Some("PCSK9")), Some("LOF"));
        // ... and the rest of the shipped table survives that override.
        assert_eq!(cfg.effective_mechanism(Some("MYH7")), Some("DOMINANT_NEGATIVE"));
    }

    #[test]
    fn test_gene_mechanisms_table_can_be_replaced_wholesale() {
        // Documented behaviour in ACMG.md: naming the table in TOML replaces
        // it, exactly as `homology_unreliable_genes` and `ba1_exceptions` do.
        let cfg: AcmgConfig = toml::from_str(
            "[gene_mechanisms]\nMYGENE = \"GOF\"\n",
        )
        .unwrap();
        assert_eq!(cfg.effective_mechanism(Some("MYGENE")), Some("GOF"));
        assert_eq!(cfg.effective_mechanism(Some("PCSK9")), None);
    }

    #[test]
    fn test_toml_deserialization() {
        let toml_str = r#"
ba1_af_threshold = 0.05
bs1_af_threshold = 0.005
pm2_af_threshold = 0.00005

[gene_overrides.BRCA1]
mechanism = "LOF"
bs1_af_threshold = 0.001
disabled_criteria = ["BP1"]
"#;
        let cfg: AcmgConfig = toml::from_str(toml_str).unwrap();
        assert!((cfg.bs1_af_threshold - 0.005).abs() < 1e-10);
        assert!(cfg.is_criterion_disabled("BRCA1", "BP1"));
        assert!(!cfg.is_criterion_disabled("BRCA1", "PM2"));
    }
}
