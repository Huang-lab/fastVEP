# v1 — Baseline run

This is the **baseline run** preserved for the v1 → v7 comparison.

## What was loaded

| SA source | Loaded? | Notes |
|-----------|:-------:|-------|
| ClinVar (.osa)            | ✅ | from `clinvar.vcf.gz` |
| ClinVar protein (.oga)    | ✅ | from `variant_summary.txt.gz` |
| gnomAD v4.1 exomes (.osa, per-chrom)         | ✅ | 25 chromosomes |
| gnomAD v4.1 gene constraints (.oga)         | ✅ | 18,173 genes |
| REVEL v1.3 (.osa, per-chrom)        | ✅ | 24 chromosomes |
| **PhyloP** (.osa)                   | ❌ | NOT loaded |
| **SpliceAI** (.osa)                 | ❌ | NOT loaded |
| **ClinGen Gene-Disease Validity (.oga)**     | ❌ | NOT loaded |

## Headline metrics

| Metric | Value |
|--------|------:|
| Same-direction concordance | 54.7 % |
| Exact match | 52.7 % |
| Opposite direction | 0.005 % |
| Likely_benign recall | **3.2 %** ← collapsed because BP7 had no PhyloP/SpliceAI to gate on |
| Benign recall | 33.2 % |
| Pathogenic recall | 15.7 % |
| Likely_pathogenic recall | 20.9 % |

## Files

- `concordance_matrix.csv` — 5-class truth × predicted matrix
- `concordance_summary.txt` — text rollup
- (No raw VCF.gz; the annotation phase was JSON output that has since been
  superseded. The matrix above was reconstructed from documentation.)

For the v7 run with full SA stack and bug fixes, see `../output_v7/`.
For the per-version SA stack and code-fix delta, see `../RUN_VERSIONS.md`.
