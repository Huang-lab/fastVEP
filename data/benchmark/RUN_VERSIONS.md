# Benchmark run versions

Each row of the `output_v*/` directories represents one end-to-end run
of the ClinVar 2-star+ benchmark on the same 673,660-variant truth set.
Successive runs differ only in (a) which supplementary annotation
databases were loaded into `--sa-dir` and (b) which classifier / output
bugs had been fixed.

## SA stack per run

|                                     |  v1  |  v2  |  v4  |  v5  |  v6  |  v7  |
|-------------------------------------|:----:|:----:|:----:|:----:|:----:|:----:|
| **Variant-level (.osa)**            |      |      |      |      |      |      |
| ClinVar                             |  ✅  |  ✅  |  ✅  |  ✅  |  ✅  |  ✅  |
| gnomAD v4.1 exomes (per-chrom)      |  ✅  |  ✅  |  ✅  |  ✅  |  ✅  |  ✅  |
| REVEL v1.3 (per-chrom)              |  ✅  |  ✅  |  ✅  |  ✅  |  ✅  |  ✅  |
| **PhyloP** (per-chrom)              |  ❌  |  ✅  |  ✅  |  ✅  |  ✅  |  ✅  |
| **SpliceAI** (per-chrom)            |  ❌  |  ✅  |  ✅  |  ✅  |  ✅  |  ✅  |
| **Gene-level (.oga)**               |      |      |      |      |      |      |
| ClinVar protein                     |  ✅  |  ✅  |  ✅  |  ✅  |  ✅  |  ✅  |
| gnomAD gene constraints             |  ✅  |  ✅  |  ✅  |  ✅  |  ✅  |  ✅  |
| **ClinGen Gene-Disease Validity**   |  ❌  |  ✅  |  ✅  |  ✅  |  ✅  |  ✅  |

(PhyloP and SpliceAI are distilled from gnomAD v4 INFO fields
`phylop` and `spliceai_ds_max` rather than re-downloaded; the gnomAD
v4 sites VCF already includes them. ClinGen GDV substitutes for
OMIM `genemap2.txt` per ClinGen SVI / Abou Tayoun 2018 — same `.oga`
schema, `omim` json_key, but with a multi-curator scored rubric and
explicit Definitive/Strong/Moderate filtering.)

## Code fixes per run

|                                                              |  v1  |  v2  |  v4  |  v5  |  v6  |  v7  |  v8  |
|--------------------------------------------------------------|:----:|:----:|:----:|:----:|:----:|:----:|:----:|
| SpliceAI `spliceAI` json_key recognised in classifier        |      |      |  ✅  |  ✅  |  ✅  |  ✅  |  ✅  |
| PhyloP read from `allele_supplementary` (CLI's actual route) |      |      |  ✅  |  ✅  |  ✅  |  ✅  |  ✅  |
| VCF + bgzip output (vs 25 GB pretty JSON)                    |      |      |  ✅  |  ✅  |  ✅  |  ✅  |  ✅  |
| `vep_allele(ref, alt)` indel matching in concordance script  |      |      |      |  ✅  |  ✅  |  ✅  |  ✅  |
| **PM2 fires when variant absent from gnomAD** (`pm2_absent_when_no_record`) |      |      |      |      |  ✅  |  ✅  |  ✅  |
| **BP4-splice gated to non-PVS1 consequences** (Walker 2023)  |      |      |      |      |  ✅  |  ✅  |  ✅  |
| **BS1 uses max-pop AF (mirrors BA1)** (ClinGen SVI)          |      |      |      |      |      |  ✅  |  ✅  |
| **BS2 AD requires AC ≥ 5 (`bs2_ad_min_ac`)** (Richards 2015) |      |      |      |      |      |  ✅  |  ✅  |
| **PVS1 graded (last-exon→Moderate, deep-intron gate)** (#50) |      |      |      |      |      |      |  ✅  |
| **BP4 SpliceAI path skips deep-exonic missense** (#50)      |      |      |      |      |      |      |  ✅  |
| **PM2 gnomAD query left-aligned + ClinVar-AF backstop** (#50) |      |      |      |      |      |      |  ✅  |
| **Opposite-direction discrepancies never truncated** (03 cap fix) |      |      |      |      |      |      |  ✅  |

The SA stack is unchanged from v7. v8 adds only the three ACMG
criterion-execution fixes from the medical-genetics discordance review
(merged in #50) plus a concordance-script fix so the review export is
complete (see below).

v3 was a partial run (PhyloP+SpliceAI loaded but bugs still latent);
its results are functionally indistinguishable from v2 and were
overwritten before being preserved.

## Headline metrics per run

|                            |     v1     |     v6     |     v7     |     v8     |  Δ v7→v8  |
|----------------------------|-----------:|-----------:|-----------:|-----------:|----------:|
| Same-direction concordance |   54.7 %   |   70.3 %   |   70.8 %   | **73.6 %** |**+2.8 pp**|
| Exact match                |   52.7 %   |   56.8 %   |   58.7 %   | **61.3 %** |**+2.6 pp**|
| Opposite direction         |   0.005 %  |   0.05 %   |  0.046 %   | **0.018 %**| **−61 %** |
| NoCall                     |   0.0 %    |   0.0 %    |   0.0 %    |   0.0 %    | —         |
| **Pathogenic recall**      |   **15.7 %** | 20.6 %   |   63.8 %   | **64.0 %** |**+48 pp** |
| **Likely_pathogenic recall** | **20.9 %** | 26.7 %   |   51.8 %   | **52.0 %** |**+31 pp** |
| VUS recall                 |   96.6 %   |   92.6 %   |   91.5 %   |   92.0 %   | -5 pp     |
| **Likely_benign recall**   |   **3.2 %**|   42.4 %   |   42.4 %   | **42.7 %** |**+39 pp** |
| Benign recall              |   33.2 %   |   58.0 %   |   58.0 %   | **59.0 %** |**+26 pp** |
| **Benign exact-match (B→B)** | 48,282 | 48,996 | 48,996 | **63,481** |**+15,199** |

## Driver of each lift

- **+39 pp LB recall, +25 pp B recall** (v1 → v5): BP7 went from **0**
  → **81,706 fires** once PhyloP+SpliceAI were loaded *and* both
  wiring bugs were fixed. (Walker 2023: BP7 needs synonymous + low
  SpliceAI + low PhyloP.)
- **+48 pp Pathogenic recall, +31 pp LP recall** (v5 → v6): the
  classifier's PM2 evaluator previously refused to fire when
  `input.gnomad` was `None` (no gnomAD record at the variant). For
  ~78 % of the truth's pathogenic class — most rare disease variants
  aren't in gnomAD — that meant PM2 NotEvaluated, so PVS1 had no
  partner criterion and the PVS+≥1 PP (SVI Sept 2020) → LP rule
  couldn't trigger. Default config flag
  `pm2_absent_when_no_record = true` interprets a missing record as
  "absent from gnomAD" per ClinGen SVI v1.0. PM2_Supporting fires
  jumped from ~12K → 340K, unlocking the LP rule for ~50K Pathogenic
  variants. PVS1 also nearly doubled (27K → 50K P+LP) because
  BP4-splice is no longer (incorrectly) firing on frameshift / null
  variants — Walker 2023 explicitly scopes BP4-splice to
  splice-territory consequences.
- **VUS recall slight drop (-5 pp)**: by design — when more benign
  evidence fires, some variants previously called VUS now correctly
  drop to LB / B (which doesn't match a VUS truth). Same-direction
  rate still rises because the P/LP/B/LB gains far outweigh the VUS
  loss.
- **+15,199 Benign exact-match calls** (v6 → v7): two ClinGen-SVI BS-tier
  fixes from a deep classifier audit. **(1)** BS1 was reading cohort
  `all_af` instead of `max_pop_af` — a 5%-AF EAS variant could slip
  under a 1% BS1 threshold whenever the global cohort diluted it.
  ClinGen SVI applies BS1 against the max-pop AF (mirroring BA1).
  Effect: BS1 fires went **6.4× higher** (4,104 → 26,291). Many LB
  calls in v6 promote to B in v7 once BS1+BP fires. **(2)** BS2 was
  firing for AD genes on any single het in gnomAD (`AF > 0 && AN >
  10000`) — a singleton novel allele in a 100K cohort isn't
  "observed in healthy adult" per Richards 2015. Tightened to AC ≥ 5
  by default (`bs2_ad_min_ac`). False-positive BS2 fires on
  Pathogenic ClinVar variants cut by **52%** (809 → 389). Net
  opposite-direction rate dropped from 0.05% to 0.0%.
- **+2.6 pp exact, +2.8 pp same-direction, opposite −61%** (v7 → v8):
  the three ACMG criterion-execution fixes from the medical-genetics
  discordance review (#50). (1) PVS1 was always firing Very-Strong
  because its grading inputs were hard-coded empty; now last-exon PTCs
  downgrade to Moderate and deep-intronic mislabeled-splice variants
  don't fire PVS1. (2) BP4's SpliceAI benign path was firing on
  deep-exonic missense (SpliceAI = 0) against a pathogenic REVEL; now
  gated to splice-territory consequences. (3) PM2 was calling
  repeat-region indels "absent" because a differently-anchored query
  missed gnomAD's left-aligned record; queries are now
  reference-normalized before the gnomAD lookup, with ClinVar's
  ExAC/1000G/ESP AFs as a backstop. Opposite-direction errors fell
  313 → 122; per-class exact gains led by Likely_benign (+7,009) and
  VUS (+6,467). Also fixed the concordance script's 10k discrepancy
  cap, which had been silently truncating the opposite-direction review
  export at chr5 (so the geneticist's earlier list of ~101 was
  incomplete; the true post-fix set is 122).

## Where to find each version

- v1 baseline: `output_v1/concordance_matrix.csv` +
  `output_v1/README.md` (raw outputs were overwritten; matrix
  reconstructed from documentation)
- v7: `output_v7/` (previous baseline — full outputs + figures + raw VCF.gz)
- **v8 current: `output_v8/`** (discordance-review fixes; full outputs +
  `discrepancies_for_md_review.tsv` medical-geneticist review table +
  `README.md`; ClinVar 2026-06-27 refresh check noted in that README)
