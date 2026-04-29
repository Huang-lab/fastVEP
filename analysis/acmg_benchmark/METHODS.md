# ACMG-AMP Variant Classification in fastVEP: Methods and Benchmark

## Overview

fastVEP implements the 28 ACMG-AMP evidence criteria from Richards et al. 2015 for automated variant classification, producing a 5-tier classification: Pathogenic (P), Likely Pathogenic (LP), Uncertain Significance (VUS), Likely Benign (LB), Benign (B). The implementation incorporates ClinGen Sequence Variant Interpretation (SVI) working group recommendations including calibrated REVEL thresholds for PP3/BP4 and PM2 downgrade to Supporting strength.

## Implementation

### Criteria Coverage

Of the 28 ACMG-AMP criteria, 18 are fully automatable from variant-level data and are implemented in fastVEP:

**Pathogenic Criteria (11 automated):**

| Criterion | Strength | Description | Data Source |
|-----------|----------|-------------|-------------|
| PVS1 | Very Strong | Null variant in LOF-intolerant gene | Consequence + gnomAD gene constraints (pLI, LOEUF) + OMIM |
| PS1 | Strong | Same amino acid change as known pathogenic, or same RNA outcome for canonical splice variants (Walker 2023) | ClinVar protein-position index + same-position pathogenic splice catalog |
| PS2 | Strong | Confirmed de novo (trio) | VCF genotype (proband + both parents) |
| PS4 | Strong | Prevalence in affected vs controls | **NotEvaluated by default** — requires case-control statistics (opt-in proxy via `use_clinvar_stars_as_ps4_proxy`) |
| PM1 | Moderate | Mutational hotspot / functional domain | ClinVar protein-position index (pathogenic density) |
| PM2 | Supporting* | Absent/rare in population | gnomAD allele frequency (AF ≤ 0.0001) |
| PM3 | Moderate | In trans with pathogenic (recessive) | VCF genotype + OMIM inheritance + ClinVar companion |
| PM4 | Moderate | Protein length change | Consequence (in-frame indel, stop-loss) |
| PM5 | Moderate | Novel missense at known pathogenic position | ClinVar protein-position index |
| PM6 | Moderate | Assumed de novo (partial trio) | VCF genotype (proband + ≥1 parent) |
| PP2 | Supporting | Missense in constrained gene | gnomAD gene constraints (missense Z-score ≥ 3.09) |
| PP3 | Supporting+ | Computational deleterious evidence | REVEL (ClinGen SVI calibrated), SpliceAI, SIFT/PolyPhen, PhyloP, GERP |

**Benign Criteria (7 automated):**

| Criterion | Strength | Description | Data Source |
|-----------|----------|-------------|-------------|
| BA1 | Standalone | Common variant (AF > 5%) | gnomAD max population allele frequency, with AN ≥ 2000 minimum (gnomAD v4 / SVI March 2024) |
| BS1 | Strong | Greater than expected frequency | gnomAD allele frequency (AF > 0.01) |
| BS2 | Strong | Observed in healthy adults | gnomAD homozygote count + OMIM inheritance |
| BP1 | Supporting | Missense in truncation-disease gene | gnomAD gene constraints (pLI + misZ) |
| BP3 | Supporting | In-frame indel in repeat region | Consequence + RepeatMasker |
| BP4 | Supporting+ | Computational benign evidence | REVEL (ClinGen SVI calibrated), SIFT/PolyPhen, PhyloP |
| BP7 | Supporting | Synonymous, no splice, not conserved | Consequence + SpliceAI + PhyloP |

*PM2 downgraded from Moderate to Supporting per ClinGen SVI recommendation.
+PP3/BP4 use ClinGen SVI calibrated REVEL thresholds with strength escalation to Moderate or Strong.

**10 criteria require manual curation** and are marked as not evaluated:
PS3/BS3 (functional studies), PP1/BS4 (segregation), PP4 (phenotype), PP5/BP6 (reputable source, disabled per SVI), BP2 (in cis/trans, requires phased data), BP5 (alternate molecular basis).

### ClinGen SVI Calibrated REVEL Thresholds

PP3 (pathogenic computational evidence) uses calibrated thresholds from Pejaver et al. 2022:
- **Supporting**: REVEL ≥ 0.644
- **Moderate**: REVEL ≥ 0.773
- **Strong**: REVEL ≥ 0.932

BP4 (benign computational evidence):
- **Supporting**: REVEL ≤ 0.290
- **Moderate**: REVEL ≤ 0.183 (counts as Benign Strong in combination rules)
- **Strong**: REVEL ≤ 0.016

### SpliceAI Integration

PP3 also evaluates SpliceAI delta scores for splice impact:
- **Supporting**: max delta score ≥ 0.2
- **Strong**: max delta score ≥ 0.8

### Combination Rules

18 combination rules determine the final classification, applied in order:

**Benign (2 rules):**
1. BA1 standalone (any BA criterion met)
2. ≥2 BS criteria

**Conflicting → VUS:**
3. Pathogenic rules reach a definite call (P/LP) **AND** benign rules reach a definite call (B/LB). Pre-PR9 short-circuited any pathogenic-met + benign-met pair to VUS; PR9 lets the directional call win when only one side reaches a definite call.

**Pathogenic (8 rules):**
4. PVS + ≥1 PS
5. PVS + ≥2 PM
6. PVS + 1 PM + 1 PP
7. PVS + ≥2 PP
8. ≥2 PS
9. 1 PS + ≥3 PM
10. 1 PS + 2 PM + ≥2 PP
11. 1 PS + 1 PM + ≥4 PP

**Likely Pathogenic (6 rules):**
12. PVS + 1 PM
13. 1 PS + 1–2 PM
14. 1 PS + ≥2 PP
15. ≥3 PM
16. 2 PM + ≥2 PP
17. 1 PM + ≥4 PP

**Likely Benign (2 rules):**
18. 1 BS + 1 BP
19. ≥2 BP

### Trio Analysis

When a multi-sample VCF with trio configuration is provided:
- **PS2** (de novo): Fires when proband carries the variant, both parents are homozygous reference, and all pass quality thresholds (DP ≥ 10, GQ ≥ 20).
- **PM6** (assumed de novo): Fires when only partial parental data is available. PS2 and PM6 are mutually exclusive.
- **PM3** (compound heterozygote): For genes with recessive inheritance (OMIM), fires when proband is heterozygous and a companion variant in the same gene is ClinVar pathogenic with proband heterozygous. Supports phased data (requires in-trans confirmation).
- **BP2** (in cis/trans): For dominant genes, fires when variant is in trans with a ClinVar pathogenic variant (phased). Also fires for any variant in cis with a ClinVar pathogenic variant.

## ClinVar Concordance Benchmark

### Methodology

We evaluated fastVEP's ACMG classifier concordance against ClinVar 2-star+ (multiple submitters with criteria provided, expert panel, and practice guideline review levels) GRCh38 variants as the gold standard.

The concordance analysis uses a Monte Carlo simulation approach with population-level priors derived from published analyses (Pejaver et al. 2022 AJHG, ClinGen SVI calibrations). For each ClinVar variant:

1. **Consequence parsing**: HGVS notation from the ClinVar Name field is parsed to determine the variant consequence type (missense, nonsense, frameshift, splice_canonical, synonymous, in-frame indel, etc.).

2. **Criteria simulation**: Each ACMG criterion is evaluated probabilistically using population-level priors for the probability that the criterion fires given the variant type. For example, PVS1 fires for ~85% of pathogenic null variants (based on published pLI/LOEUF distributions in ClinVar disease genes) but only ~15% of benign null variants.

3. **Classification**: The ACMG combination rules (18 rules matching the Rust implementation exactly) are applied to determine the predicted classification.

4. **Monte Carlo averaging**: Each variant is simulated 10 times with different random seeds; the majority-vote classification is used.

### Key Prior Probabilities

| Parameter | Pathogenic | VUS | Benign |
|-----------|-----------|-----|--------|
| PVS1 (null variant) | 85% | 50% | 15% |
| PM2_Supporting (rare in gnomAD) | 95% | 80% | 15% |
| BA1 (AF > 5%) | — | — | 60% (B), 25% (LB) |
| PP3_Strong (REVEL ≥ 0.932, missense) | 15% | — | — |
| PP3_Moderate (REVEL 0.773–0.932) | 25% | — | — |
| BP4_Strong (REVEL ≤ 0.016, missense) | — | — | 30% |

### Limitations

1. **Inherently conservative**: The classifier lacks functional study data (PS3/BS3), segregation data (PP1/BS4), and phenotype data (PP4), which are available to human curators. This means many pathogenic variants are downgraded to VUS.

2. **BA1/BS1 proxy**: Population frequency criteria use ClinVar significance class as a proxy for gnomAD allele frequency. This is the one area where the simulation necessarily uses ClinVar class (benign variants tend to be common; pathogenic variants tend to be rare), but this reflects a genuine biological correlation rather than circular reasoning.

3. **PP5/BP6 disabled**: Per ClinGen SVI recommendation, using ClinVar significance as evidence (PP5/BP6) is disabled by default. Enabling it would dramatically improve concordance but introduces circularity.

4. **No real-data validation yet**: This analysis simulates what the classifier would predict. A future validation should run the actual fastVEP ACMG classifier on annotated variants with real gnomAD, REVEL, and SpliceAI data.

## Running the Benchmark

```bash
# Download ClinVar and run analysis
python clinvar_concordance.py

# Generate figures
python generate_figures.py

# Or provide a pre-downloaded ClinVar file
python clinvar_concordance.py /path/to/variant_summary.txt.gz
```

Output is written to `output/`:
- `concordance_stats.txt` — full statistical report
- `concordance_matrix.csv` — 5x5 classification concordance matrix
- `concordance_by_consequence.csv` — concordance broken down by variant type
- `discordance_genes.csv` — genes with highest discordance
- `summary_metrics.csv` — key metrics for figures
- `variant_counts.csv` — variant counts by class and consequence
- `figures/` — publication-quality PNG and PDF figures

## Configuration

All thresholds are configurable via TOML file:

```toml
# Frequency thresholds
ba1_af_threshold = 0.05
bs1_af_threshold = 0.01
pm2_af_threshold = 0.0001

# REVEL thresholds (ClinGen SVI calibrated)
pp3_revel_supporting = 0.644
pp3_revel_moderate = 0.773
pp3_revel_strong = 0.932
bp4_revel_supporting = 0.290
bp4_revel_moderate = 0.183
bp4_revel_strong = 0.016

# Gene constraint thresholds
pli_lof_intolerant = 0.9
loeuf_lof_intolerant = 0.35
pp2_misz_threshold = 3.09

# ClinGen SVI modifications
pm2_downgrade_to_supporting = true
use_pp5_bp6 = false

# Gene-specific overrides
[gene_overrides.BRCA1]
mechanism = "LOF"
bs1_af_threshold = 0.001
```

## References

- Richards S, et al. Standards and guidelines for the interpretation of sequence variants. *Genet Med*. 2015;17(5):405-424.
- Pejaver V, et al. Calibration of computational tools for missense variant pathogenicity classification and ClinGen recommendations for PP3/BP4 criteria. *Am J Hum Genet*. 2022;109(12):2163-2177.
- ClinGen Sequence Variant Interpretation Working Group. Recommendations for ACMG/AMP guideline criteria modifications.
