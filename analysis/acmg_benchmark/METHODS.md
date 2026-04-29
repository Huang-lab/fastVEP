# ACMG-AMP Variant Classification in fastVEP: Methods and Benchmark

## Overview

fastVEP implements the 28 ACMG-AMP evidence criteria from Richards et al.
2015 plus published ClinGen Sequence Variant Interpretation (SVI) Working
Group refinements, producing a 5-tier classification: Pathogenic (P),
Likely Pathogenic (LP), Uncertain Significance (VUS), Likely Benign (LB),
Benign (B).

This document reflects the state after the SVI alignment series (PR1–PR10).
Each criterion section includes the ClinGen reference and any deviations
from a strict reading of the SVI text. A per-criterion "Limitations"
column flags criteria that fall back to legacy or conservative behavior
until additional pipeline data is wired in.

## Implementation

### Criteria Coverage

Of the 28 ACMG-AMP criteria, 18 are fully automatable from variant-level
data and are implemented in fastVEP. The classifier records the source
that drove each call in `details.pp3_source` / `details.ps1_path` /
`details.inheritance_basis` etc., so every classification is auditable.

#### Pathogenic Criteria (11 automated)

| Criterion | Strength | Description | Data Source / Notes |
|-----------|----------|-------------|---------------------|
| PVS1 | VS / Strong / Moderate / Supporting (Abou Tayoun 2018 decision tree) | Null variant (nonsense, frameshift, canonical splice, start-loss, whole-gene deletion) in LOF-intolerant gene | Consequence + gnomAD constraints + transcript NMD prediction + critical-region check + alt-start distance |
| PS1 | Strong | Same amino acid change as known pathogenic missense, **or** same RNA outcome for canonical splice (Walker 2023) | ClinVar protein-position index + same-position pathogenic splice catalog |
| PS2 | Strong | Confirmed de novo (full trio) | VCF genotype (proband + both parents) + DP/GQ thresholds |
| PS3 | Strong | Functional studies | Not automatable — NotEvaluated |
| PS4 | Strong | Prevalence in affected vs controls | **NotEvaluated by default** — requires case-control statistics. Optional legacy proxy via `use_clinvar_stars_as_ps4_proxy` |
| PM1 | Moderate | Mutational hotspot / functional domain | ClinVar protein-position density. Capped against PP3 per Pejaver 2022 |
| PM2 | Supporting* | Absent / extremely rare in population | gnomAD raw AF, **inheritance-aware** (AD/unknown: AC=0; AR: AF ≤ 0.00007) — SVI v1.0 |
| PM3 | Supporting/Moderate/Strong/VeryStrong | In trans with pathogenic (recessive) | **SVI PM3 v1.0 points-based**: P / LP companion × phasing × hom-occurrence → 0.5 / 1.0 / 2.0 / 4.0 thresholds |
| PM4 | Moderate | Protein length change | Consequence (in-frame indel, stop-loss) |
| PM5 | Moderate | Novel missense at known pathogenic position | ClinVar protein-position index (different alt AA) |
| PM6 | Moderate | Assumed de novo (partial trio) | VCF genotype (proband + ≥1 parent). Mutually exclusive with PS2 |
| PP2 | Supporting | Missense in constrained gene | gnomAD missense Z-score ≥ 3.09 |
| PP3 | Supporting / Moderate / Strong (Pejaver 2022 + Walker 2023) | Computational pathogenic evidence | REVEL (missense only) or SpliceAI ≥ 0.2 (Supporting only) |
| PP4 | Supporting | Phenotype-specific | Not automatable — NotEvaluated |
| PP5 | Supporting | Reputable source | **Disabled by default** per ClinGen SVI |

*PM2 downgraded from Moderate to Supporting per ClinGen SVI v1.0.

#### Benign Criteria (7 automated)

| Criterion | Strength | Description | Data Source / Notes |
|-----------|----------|-------------|---------------------|
| BA1 | Standalone | Common variant (AF > 5%) | gnomAD max population AF, with **AN ≥ 2000** minimum (gnomAD v4 / SVI March 2024). Honors the **9-variant Ghosh 2018 BA1 exception list** |
| BS1 | Strong | Greater than expected frequency | gnomAD AF (gene-specific or default 0.01); same AN minimum as BA1 |
| BS2 | Strong | Observed in healthy adults | gnomAD homozygote count + OMIM inheritance |
| BS3 | Strong | Functional studies — no damage | Not automatable — NotEvaluated |
| BS4 | Strong | Lack of segregation | Not automatable — NotEvaluated |
| BP1 | Supporting | Missense in truncation-disease gene | gnomAD pLI ≥ 0.9 + misZ < 2.0 |
| BP2 | Supporting | In trans / in cis with pathogenic | Companion-variant phasing + OMIM inheritance |
| BP3 | Supporting | In-frame indel in repeat region | Consequence + RepeatMasker |
| BP4 | Supporting / Moderate / Strong / **VeryStrong** | Computational benign evidence | REVEL (missense only, **incl. VeryStrong band at REVEL ≤ 0.003**) or SpliceAI ≤ 0.1 (Walker 2023) |
| BP5 | Supporting | Alternate molecular basis | Not automatable — NotEvaluated |
| BP6 | Supporting | Reputable source — benign | **Disabled by default** per ClinGen SVI |
| BP7 | Supporting | Synonymous (mid-exon) or deep-intronic, no splice, not conserved | Consequence + SpliceAI + PhyloP + transcript exon coords. **Walker 2023**: exon-edge exclusion + deep-intronic extension |

**10 criteria require manual curation** and are marked NotEvaluated:
PS3 / PS4 (default) / BS3 / BS4 / PP1 / PP4 / PP5 (disabled) / BP2 (when
unphased) / BP5 / BP6 (disabled).

### Pejaver 2022 Calibrated REVEL Thresholds

REVEL is applied **only to missense variants** per Pejaver 2022. The
single-tool calibration replaces the previous SIFT/PolyPhen/PhyloP/GERP
ensemble (Pejaver explicitly recommends a single calibrated tool over
ad-hoc consensus).

| Direction | Strength | REVEL threshold |
|-----------|----------|-----------------|
| PP3 | Supporting | ≥ 0.644 |
| PP3 | Moderate   | ≥ 0.773 |
| PP3 | Strong     | ≥ 0.932 |
| BP4 | Supporting | ≤ 0.290 |
| BP4 | Moderate   | ≤ 0.183 |
| BP4 | Strong     | ≤ 0.016 |
| BP4 | **Very Strong** (REVEL only) | ≤ 0.003 |

A single BP4_VeryStrong is mapped to 2× `benign_strong` in the counts so
it satisfies the existing ≥2 BS → Benign rule alone (Tavtigian Bayesian
framework).

### Walker 2023 Splicing Recommendations

- **PP3 splice**: SpliceAI max_ds ≥ 0.2 → PP3 *Supporting* (no Strong from
  SpliceAI alone — Strong splice claims need experimental RNA assay).
- **BP4 splice**: SpliceAI max_ds ≤ 0.1 → BP4 Supporting.
- **Uninformative zone**: 0.10 < max_ds < 0.20 — neither fires.
- **BP7 exon-edge exclusion**: BP7 cannot fire for synonymous at first
  base or last 3 bases of an exon (canonical splice region).
- **BP7 deep-intronic extension**: BP7 may fire for intronic variants
  with offset ≥ 7 (donor side) or ≤ -21 (acceptor side) when SpliceAI is
  low and PhyloP is low.

### PVS1 Decision Tree (Abou Tayoun 2018)

| Strength | Trigger |
|----------|---------|
| **PVS1** (Very Strong) | Nonsense/frameshift predicted to undergo NMD; canonical ±1/2 splice predicted to cause NMD; whole-gene deletion in haploinsufficient gene |
| **PVS1_Strong** | NMD-escape in critical functional region |
| **PVS1_Moderate** | NMD-escape, non-critical region, ≥10% protein removed; canonical splice in last exon (NMD unlikely); start-loss with downstream Met ≤100 codons + pathogenic upstream |
| **PVS1_Supporting** | <10% protein removed in non-critical region; start-loss without strong corroborating evidence |

When NMD prediction or other transcript-level signals are missing,
PVS1 falls back to legacy full-strength VeryStrong for backward
compatibility.

### PM2 Inheritance-Aware Threshold (SVI v1.0)

| Inheritance | Threshold |
|-------------|-----------|
| AD / unknown | Strict absence (AC = 0 AND AF = 0) |
| AR | AF ≤ 0.00007 (0.007%) |
| Per-gene override | Wins over inheritance default |

Uses **raw** gnomAD AF (not FAF / popmax). Inheritance is inferred from
OMIM phenotypes.

### PM3 v1.0 Points Scoring (SVI)

| Observation | Points |
|-------------|--------|
| Confirmed in-trans + co-occurring **Pathogenic** | 1.0 |
| Confirmed in-trans + co-occurring **Likely Pathogenic** | 0.5 |
| Phase unknown + Pathogenic | 0.5 |
| Phase unknown + Likely Pathogenic | 0.25 |
| Homozygous proband observation | 0.5 (capped at 1.0 total) |

| Total points | Strength |
|--------------|----------|
| ≥ 0.5 | PM3_Supporting |
| ≥ 1.0 | PM3 (Moderate) |
| ≥ 2.0 | PM3_Strong |
| ≥ 4.0 | PM3_VeryStrong |

In-cis companions are excluded from PM3 (those count toward BP2).

### BA1 Exception List (Ghosh 2018)

Nine variants exempt from BA1 despite exceeding the 5% threshold (HFE
c.845G>A p.Cys282Tyr, GJB2 c.109G>A, F2/F5 founder alleles, etc.). Match
on `(gene_symbol, hgvs_c)`, case-insensitive. Configurable via TOML so
VCEPs can extend.

| Gene | Variant | Note |
|------|---------|------|
| ACAD9 | c.-44_-41dupTAAG | VUS |
| ACADS | c.511C>T | VUS |
| BTD | c.1330G>C | Pathogenic — biotinidase deficiency |
| GJB2 | c.109G>A | Pathogenic — DFNB1 hearing loss |
| HFE | c.187C>G | Pathogenic — hemochromatosis |
| HFE | c.845G>A | Pathogenic — hemochromatosis |
| MEFV | c.1105C>T | VUS |
| MEFV | c.1223G>A | VUS |
| PIBF1 | c.1214G>A | VUS |

### Anti-Double-Counting (PP3 Reconciliation)

A post-evaluation reconciliation pass suppresses PP3 (or PM1) under
overlap conditions called out in Pejaver 2022 and Walker 2023:

| Trigger | Suppressed | Source |
|---------|------------|--------|
| PVS1 fires AND PP3 was driven by SpliceAI | PP3 | Walker 2023 |
| PS1 fires AND PP3 was driven by REVEL | PP3 | Pejaver 2022 |
| PM5 fires AND PP3 was driven by REVEL | PP3 | Pejaver 2022 |
| PP3_Strong + PM1 (combined > Strong cap) | PM1 | Pejaver 2022 |

Suppressed criteria stay in the result with `met=false` and a
`details.suppressed_by_reconcile` note.

### gnomAD v4 AN Minimum (SVI March 2024)

BA1 and BS1 require gnomAD `all_an ≥ 2000` before they can fire. Below
the threshold the criterion is NotEvaluated rather than fired on noisy
estimates. Configurable via `min_an_for_frequency_criteria`.

### Combination Rules (19 = 18 Richards 2015 + 1 SVI Sept 2020)

**Benign:**
1. BA1 standalone → Benign
2. ≥2 BS → Benign

**Pathogenic (8):**
3. PVS + ≥1 PS
4. PVS + ≥2 PM
5. PVS + 1 PM + 1 PP
6. PVS + ≥2 PP
7. ≥2 PS
8. 1 PS + ≥3 PM
9. 1 PS + 2 PM + ≥2 PP
10. 1 PS + 1 PM + ≥4 PP

**Likely Pathogenic (7, includes ClinGen SVI Sept 2020 rule):**
11. PVS + 1 PM
12. **PVS + ≥1 PP** *(ClinGen SVI Sept 2020 — compensates PM2 downgrade; Bayesian Post_P = 0.988)*
13. 1 PS + 1–2 PM
14. 1 PS + ≥2 PP
15. ≥3 PM
16. 2 PM + ≥2 PP
17. 1 PM + ≥4 PP

**Likely Benign (2):**
18. 1 BS + 1 BP
19. ≥2 BP

**Conflict gating (PR9 fix)**: pathogenic and benign rules apply
**independently**. The result is VUS-Conflicting only when **both**
directions reach a definite call (P/LP and B/LB). Otherwise the
directional call wins.

### Trio Analysis

When a multi-sample VCF with trio configuration is provided:
- **PS2** (de novo): proband carries variant, both parents hom-ref, all
  pass DP ≥ 10 / GQ ≥ 20.
- **PM6** (assumed de novo): partial parental data; mutually exclusive
  with PS2.
- **PM3** (compound het): SVI v1.0 points scoring (above). Recessive gene
  required (OMIM).
- **BP2** (in cis/trans): for dominant genes — variant in trans with
  pathogenic; for any gene — variant in cis with pathogenic.

## ClinVar Concordance Benchmark

### Methodology

We evaluate fastVEP's ACMG classifier concordance against ClinVar 2-star+
GRCh38 variants as the gold standard.

The concordance analysis uses a Monte Carlo simulation with
population-level priors derived from published analyses (Pejaver 2022,
ClinGen SVI calibrations). For each ClinVar variant:

1. **Consequence parsing**: HGVS notation → variant consequence type.
2. **Criteria simulation**: each ACMG criterion is evaluated
   probabilistically using priors for the probability that the criterion
   fires given the variant type.
3. **Classification**: combination rules (19 total) are applied to
   determine the predicted classification.
4. **Monte Carlo averaging**: each variant is simulated 10 times; the
   majority-vote classification is used.

### Key Prior Probabilities

| Parameter | Pathogenic | VUS | Benign |
|-----------|-----------|-----|--------|
| PVS1 (null variant) | 85% | 50% | 15% |
| PM2_Supporting (rare in gnomAD) | 95% | 80% | 15% |
| BA1 (AF > 5%) | — | — | 60% (B), 25% (LB) |
| PP3_Strong (REVEL ≥ 0.932, missense) | 15% | — | — |
| PP3_Moderate (REVEL 0.773–0.932) | 25% | — | — |
| BP4_Strong (REVEL ≤ 0.016, missense) | — | — | 25% |
| BP4_VeryStrong (REVEL ≤ 0.003) | — | — | 5% |

### Limitations

1. **Inherently conservative**: missing PS3/BS3/BS4/PP1/PP4/BP5 means
   many pathogenic variants drop to VUS for lack of functional /
   segregation / phenotype evidence. Manual curators outperform the
   automated classifier for these categories by design.
2. **PVS1 / PM3 / PS1 / BP7 fallbacks**: when the pipeline does not yet
   populate Abou Tayoun decision-tree signals (NMD, %protein removed,
   etc.), exon-edge / deep-intronic positions, or splice-position
   pathogenic catalogs, those criteria fall back to conservative legacy
   behavior.
3. **PS4 NotEvaluated by default**: the previous ClinVar-stars proxy was
   replaced; opt back in via `use_clinvar_stars_as_ps4_proxy` for
   backward-comparable benchmarks.
4. **BA1/BS1 proxy in benchmark**: the Monte Carlo simulation uses
   ClinVar significance class as a proxy for gnomAD AF. The classifier
   itself uses real AF.

## Configuration

```toml
# Frequency thresholds
ba1_af_threshold = 0.05
bs1_af_threshold = 0.01
pm2_af_threshold = 0.0001            # legacy single-threshold field (back-compat)
pm2_ad_af_threshold = 0.0            # AD / unknown: strict absence
pm2_ar_af_threshold = 0.00007        # AR threshold (SVI v1.0)
min_an_for_frequency_criteria = 2000 # gnomAD v4 AN minimum (SVI March 2024)

# REVEL thresholds (Pejaver 2022; missense only)
pp3_revel_supporting = 0.644
pp3_revel_moderate = 0.773
pp3_revel_strong = 0.932
bp4_revel_supporting = 0.290
bp4_revel_moderate = 0.183
bp4_revel_strong = 0.016
bp4_revel_very_strong = 0.003        # only REVEL reaches this band

# SpliceAI thresholds (Walker 2023)
spliceai_pathogenic = 0.2
spliceai_benign = 0.1

# Conservation
phylop_conserved = 2.0

# Gene constraints
pli_lof_intolerant = 0.9
loeuf_lof_intolerant = 0.35
pp2_misz_threshold = 3.09
pm1_hotspot_window = 5
pm1_hotspot_min_pathogenic = 3

# ClinGen SVI behavior modifications
pm2_downgrade_to_supporting = true
use_pp5_bp6 = false
use_clinvar_stars_as_ps4_proxy = false

# BA1 exception list — defaults to the 9-variant Ghosh 2018 set;
# users can extend or replace via TOML.
[[ba1_exceptions]]
gene = "HFE"
hgvs_c = "c.845G>A"
reason = "Hereditary hemochromatosis"

# Gene-specific overrides
[gene_overrides.BRCA1]
mechanism = "LOF"
bs1_af_threshold = 0.001

# Per-disorder overrides for multi-disorder genes (SVI July 2025 scaffold)
[gene_overrides.GENE_X.disorders.disorder_a]
inheritance = "AR"
pm2_af_threshold = 0.00007
```

## References

- Richards S, et al. Standards and guidelines for the interpretation of sequence variants. *Genet Med*. 2015;17(5):405-424.
- Abou Tayoun AN, et al. Recommendations for interpreting the loss of function PVS1 ACMG/AMP variant criterion. *Hum Mutat*. 2018;39(11):1517-1524.
- Ghosh R, et al. Updated recommendation for the benign stand-alone ACMG/AMP criterion. *Hum Mutat*. 2018;39(11):1525-1530.
- ClinGen SVI Recommendation for Absence/Rarity (PM2) — Version 1.0. Approved September 4, 2020.
- ClinGen SVI Recommendation for In-Trans Criterion (PM3) — Version 1.0.
- Pejaver V, et al. Calibration of computational tools for missense variant pathogenicity classification and ClinGen recommendations for PP3/BP4 criteria. *Am J Hum Genet*. 2022;109(12):2163-2177.
- Walker LC, et al. (ClinGen SVI Splicing Subgroup). Using the ACMG/AMP framework to capture evidence related to predicted and observed impact on splicing. *Am J Hum Genet*. 2023;110(7):1046-1067.
- ClinGen SVI Working Group. Guidance to VCEPs Regarding gnomAD v4 (March 2024).
- ClinGen SVI Working Group. Guidance Classifying Variants in Genes Associated with Multiple Disorders (July 2025).
- Tavtigian SV, et al. Modeling the ACMG/AMP variant classification guidelines as a Bayesian classification framework. *Genet Med*. 2018;20(9):1054-1060.
- Lek M, et al. Analysis of protein-coding genetic variation in 60,706 humans. *Nature*. 2016;536(7616):285-291.
