# v7 — Current run (full SA stack + 5 fixes; BS1/BS2 tightened)

Builds on v6 with two ClinGen-SVI correctness fixes for the BS-tier
benign criteria.

## Code fixes added in v7 (vs v6)

1. **BS1 uses max-population AF, not cohort `all_af`**. ClinGen SVI
   applies BS1 against the max-pop AF (mirroring BA1) so a 5%-AF
   variant in a single subpopulation can't slip under a 1% threshold
   when the global cohort dilutes it. Effect: BS1 fires went **6.4×**
   higher (4,104 → 26,291), and Benign exact-match calls went up by
   **+14,485**.
2. **BS2 requires AC ≥ 5 for AD genes** (configurable via
   `bs2_ad_min_ac`). Richards 2015 BS2 says "observed in healthy
   adult", which a single heterozygote in a 100K cohort doesn't
   satisfy — singletons are sequencing-noise plausibility, not
   evidence of tolerance. ClinGen VCEPs commonly use AC ≥ 5 for AD
   carrier counts. Effect: BS2 false-positives on Pathogenic ClinVar
   variants cut by **52%** (809 → 389).

## SA stack (unchanged from v6)

| SA source | Loaded? | Notes |
|-----------|:-------:|-------|
| ClinVar (.osa)            | ✅ | 4,402,501 records |
| ClinVar protein (.oga)    | ✅ | 4,554 genes |
| gnomAD v4.1 exomes (.osa, per-chrom)         | ✅ | 25 chromosomes |
| gnomAD v4.1 gene constraints (.oga)         | ✅ | 18,173 genes |
| REVEL v1.3 (.osa, per-chrom)        | ✅ | 24 chromosomes |
| **PhyloP** (.osa, per-chrom)        | ✅ | distilled from gnomAD v4 INFO |
| **SpliceAI** (.osa, per-chrom)      | ✅ | distilled from gnomAD v4 INFO |
| **ClinGen GDV (.oga)**     | ✅ | 2,419 Definitive/Strong/Moderate genes |

## Code fixes vs v1 (full set in v7)

1. SpliceAI camelCase mismatch in `sa_extract.rs` *(v4)*
2. PhyloP routing fix in `sa_extract.rs` *(v4)*
3. Indel allele matching in concordance script *(v5)*
4. PM2 fires when variant absent from gnomAD *(v6)*
5. BP4-splice gated to non-PVS1 consequences *(v6)*
6. **BS1 uses max-pop AF** *(v7 — new)*
7. **BS2 AD requires AC ≥ 5** *(v7 — new)*

## Headline metrics

| Metric | Value | Δ vs v1 |
|--------|------:|--------:|
| Same-direction concordance | **70.8%** | **+16.1 pp** |
| Exact match | **58.7%** | +6.0 pp |
| Opposite direction | 0.0% | (-0.005 pp) |
| Pathogenic recall | **64.0%** | +48 pp |
| Likely_pathogenic recall | **52.0%** | +31 pp |
| Likely_benign recall | **42.7%** | +39 pp |
| Benign recall | **59.0%** | +26 pp |

## Files

- `clinvar_2star.fastvep.vcf.gz` — bgzipped VCF with ACMG in CSQ INFO
- `concordance_matrix.csv` — 5-class truth × predicted matrix
- `concordance_summary.txt` — text rollup
- `concordance_by_chrom.csv`, `concordance_by_consequence.csv`
- `criterion_firing_rates.csv` — per-criterion fire counts by truth class
- `rule_distribution.csv` — top criteria-set signatures
- `discrepancies.tsv` — opposite-direction calls
- `figures/` — 6 PNG + PDF figures (incl. v1 vs v7 comparison panels)

For per-version SA stack and code-fix diff, see `../RUN_VERSIONS.md`.
