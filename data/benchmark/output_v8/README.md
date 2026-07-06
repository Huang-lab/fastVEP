# v8 — ClinVar discordance-review fixes (PVS1 grading, BP4 gate, PM2 normalization)

Builds on v7 with the three ACMG criterion-execution fixes surfaced by the
medical-genetics discordance review (see `../../../crates/fastvep-classification`,
merged in #50). Same 673,660-variant April truth set as v1–v7, so the numbers
are directly comparable down the chain.

## Code fixes added in v8 (vs v7)

1. **PVS1 is graded, not always Very-Strong.** The decision-tree inputs
   (`is_last_exon`, `intronic_offset`, `predicted_nmd`) were hard-coded empty,
   so every null variant fired at full Very-Strong. Now wired from exon
   rank/total and HGVS intronic offset: last-exon (NMD-escape) PTCs downgrade to
   Moderate; deep-intronic variants mislabeled as canonical splice no longer
   fire PVS1.
2. **BP4 no longer fires on deep-exonic missense.** The SpliceAI benign path ran
   for any non-PVS1 consequence, so a missense with SpliceAI = 0 collected BP4
   even with a clearly pathogenic REVEL. Now gated to splice-relevant
   consequences.
3. **PM2 stops calling common indels "absent."** Repeat-region indels are stored
   left-aligned in gnomAD; a differently-anchored query silently missed them →
   spurious PM2. Queries are now reference-normalized (vt-style left-alignment)
   before the gnomAD lookup, with ClinVar's ExAC/1000G/ESP AFs as a backstop.

## Headline metrics (vs v7, same April truth set)

| Metric | v7 | **v8** | Δ |
|--------|---:|-------:|--:|
| Exact match | 58.7% | **61.3%** | +2.6 pp (+17,465) |
| Same-direction | 70.8% | **73.6%** | +2.8 pp (+18,972) |
| Opposite-direction | 313 | **122** | −191 (−61%) |
| No-call | 1 | **0** | — |

Per-class exact gains: Likely_benign +7,009, VUS +6,467, Benign +2,340,
Likely_pathogenic +1,394, Pathogenic +255. One conservative trade-off:
Pathogenic *same-direction* −582 (genuine pathogenics the tighter PVS1/PM2 moved
to VUS), dwarfed by the gains and by Pathogenic opposite-errors falling 168→67.

## SA stack (unchanged from v7)

Same as v7: ClinVar + ClinVar-protein + gnomAD exomes/gene-constraints + REVEL +
PhyloP + SpliceAI + ClinGen GDV. The `clinvar.osa` here carries ClinVar's
AF_EXAC/AF_TGP/AF_ESP frequencies (the PM2 backstop, new in v8).

## Files

- `clinvar_2star.fastvep.vcf.gz` — bgzipped VCF with ACMG in the CSQ INFO field
- `concordance_summary.txt` — text rollup
- `concordance_matrix.csv` — 5-class truth × predicted matrix
- `concordance_by_chrom.csv`, `concordance_by_consequence.csv`
- `criterion_firing_rates.csv`, `rule_distribution.csv`
- `discrepancies.tsv` — all opposite-direction calls (uncapped) + up to 10k soft ones
- **`discrepancies_for_md_review.tsv`** — the 122 opposite-direction calls enriched
  with the full annotation panel, in the same format the geneticist reviewed
  (the medical-geneticist review table)
- `opposite_direction.vcf` / `.fastvep.vcf` / `.fastvep.json.gz` — the 122-variant
  subset and its HGVS/JSON re-annotation that feed the review table

## Note on the review table (important)

The `discrepancies_for_md_review.tsv` in this run has **122** opposite-direction
variants — the *complete* set. The earlier v7 export the geneticist reviewed had
only ~101, because `03_evaluate_concordance.py` capped the discrepancy log at
10,000 rows and hit the cap at chr5, silently dropping every opposite call on
chr6+. v8 fixes that (opposite-direction calls are never truncated), so:

- of the 101 she reviewed, **72 are now resolved** and **29 remain** opposite;
- **93** of the current 122 were never in her export (mostly chr6+) — new to review.

A per-variant "what changed / what held" view of her original 101 is in
`../discordance_review/what_changed_and_held.tsv`.

## ClinVar release-refresh check (not a chain entry)

The versioned chain uses the frozen April truth set. As a robustness check, the
same v8 engine was re-run against a freshly regenerated truth set from ClinVar
**2026-06-27** (684,160 variants): 61.2% exact / 73.6% same-direction / 130
opposite — essentially identical, so the gain isn't a snapshot artifact. The
source `clinvar.vcf.gz` is the 2026-06-27 release; regenerate the truth set to it
with `01_extract_clinvar_2star.py` when you want a release-matched benchmark.

For the per-version SA stack and code-fix diff, see `../RUN_VERSIONS.md`.
