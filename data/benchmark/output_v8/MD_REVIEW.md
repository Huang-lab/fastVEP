# Medical-geneticist review queue (v8 benchmark)

This document curates the variants where fastVEP and ClinVar 2-star+
disagree most strongly — i.e. where the classifier commits to a
directional call that contradicts the expert-panel review. It excludes
cases where fastVEP calls VUS due to insufficient automatable evidence
(those need data, not curation).

v8 adds the three ACMG criterion-execution fixes from the previous
medical-genetics review (PVS1 grading, BP4 deep-exonic-missense gate,
PM2 gnomAD normalization + ClinVar-AF backstop; merged in #50).

## Top-line numbers (v8, same April truth set as v7)

| Outcome | Count | What it is |
|---------|------:|------------|
| **Opposite-direction discrepancies** | **122** | Variants where fastVEP calls P/LP and ClinVar says LB/B (or vice versa). The candidate set for MD review. Down from **313** in v7 (−61%). |
| Same-class differences (P/LP/B/LB ↔ VUS) | ~10,000 | *Not* opposite-direction; reflect missing manual-curation evidence (PS3 functional, PP1 segregation, etc.) by design. |
| Same-direction calls | 496,076 | Classifier and ClinVar agree on direction (73.6%, up from 70.8%). |
| Truth-not-annotated (NoCall) | 0 | v7 had 1 (a non-coding regulatory variant); v8 has none. |

## Important: this list is now complete (v7's was truncated)

The v7 review export the geneticist worked from listed only **101**
opposite-direction variants, even though v7's true opposite count was
**313**. The concordance script capped its discrepancy log at 10,000
rows and hit the cap at chr5, silently dropping every opposite call on
chr6+. v8 fixes that (opposite-direction calls are never truncated), so
this `discrepancies_for_md_review.tsv` contains the **complete 122**.

Relative to the 101 the geneticist already reviewed:

- **72 are now resolved** by the v8 fixes (no longer opposite);
- **29 remain** opposite (the "held" cases — gene mechanism, founder /
  penetrance exceptions, gene-specific ACMG rules, malformed ClinVar
  records, or local-data disagreement; see
  `../discordance_review/what_changed_and_held.tsv`);
- **93 are new to review** — they were always opposite but the old
  truncated export (chr1–5 only) never surfaced them.

## Curated review file

`discrepancies_for_md_review.tsv` ranks all 122 opposite-direction
cases by:

1. **3-star (expert panel) > 2-star** — the truth is more reliable, so
   classifier disagreement is more interesting.
2. **Extreme reversals (Pathogenic ↔ Benign) > one-step (LP ↔ LB)**.
3. **More criteria fired** — higher classifier confidence in the
   disagreement.

Each row carries the full annotation panel (ClinVar HGVS / phenotype /
review status / population AFs; fastVEP HGVSc/p, consequence, MANE,
REVEL, SpliceAI components, PhyloP, gnomAD AFs, gene constraints,
ClinGen/OMIM phenotypes, ClinVar-protein neighbours, and the fired ACMG
criteria) so each call can be adjudicated without leaving the table.

## ClinVar release-refresh check

The versioned chain uses the frozen April truth set for comparability.
As a robustness check, the same v8 engine re-run against a fresh ClinVar
**2026-06-27** truth set (684,160 variants) gave 61.2% exact / 73.6%
same-direction / 130 opposite — essentially identical. See `README.md`.
