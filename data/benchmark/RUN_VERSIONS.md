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

## v9: round-2 medical-genetics review

v9 responds to the second review round (`discordant 122.xlsx` plus the
reviewer's covering email). Every change traces to a specific point she raised.
Full write-up in [`../../docs/ACMG_EXPERT_REVIEW_ROUND2.md`](../../docs/ACMG_EXPERT_REVIEW_ROUND2.md);
run-specific notes in [`output_v9/README.md`](output_v9/README.md).

The SA stack is unchanged except for a rebuild of `clinvar_protein.oga`, which
now carries `n` per entry: the number of distinct nucleotide changes producing
that amino-acid change. PS1 needs it to tell an independent precedent from the
variant's own ClinVar record.

### Annotation bugs (found while reproducing her cases)

| Bug | Evidence |
|---|---|
| `ALT=.` parsed as a real allele | ClinVar reference-agreement records produced `T/.` amino acids, `cgA/cg.` codons and full Likely Pathogenic calls. 382 records. |
| MNV/delins translated in the wrong codon frame | A change straddling a codon boundary dropped the bases past the first codon, and reverse-strand alt bases were complemented without being reversed. MUTYH `GC>AT` and HPS4 `GA>CT` were reported `stop_gained`; both are synonymous. |
| PVS1 fired on splice-site indels that leave the dinucleotide intact | `PTEN c.802-2dupA` and `BRIP1 c.2258-2dup` insert the base the acceptor already carries. ATM `GTAATC>G` had SpliceAI 0.00. |

### Criteria changes

| Change | Criterion firings (v8 → v9) |
|---|---|
| BP1 blocked when the gene has an established pathogenic-missense spectrum | 22,567 → 6,471 |
| BS2 homozygote test scaled to cohort size (95 % lower bound vs a prevalence bar, floor of 2) | 101,853 → 67,615 |
| PM1 restricted to missense/in-frame and suppressed under PVS1 | 41,845 → 23,323 |
| BP4 not applied to in-frame indels / stop-lost / protein-altering | included above |
| BA1/BS1/BS2/PM2 NotEvaluated for homology-confounded genes (PMID 27228465) | included above |
| BS1/BS2 NotEvaluated for ClinVar low-penetrance / risk alleles | included above |
| Strong-or-above evidence on one side blocks a definite call on the other | combiner rule |
| PS1 requires an independent nucleotide change | 23,022 → 1,274 |

### The v9 honesty correction

`clinvar_protein.oga` is built from ClinVar and contained the variant being
classified, so PS1 was matching each variant against its own ClinVar record.
**94 % of v8's PS1 firings were self-matches** (23,022 → 1,274 once excluded).
PM1's hotspot count had the same problem, with the variant contributing one of
the three pathogenic neighbours it needed.

Reported both ways so the effect is visible:

- **Criteria fixes alone** (v8 → v9 CI): same-direction 73.6 % → 73.0 %, opposite-direction **122 → 41 (−66 %)**. The reviewer's feedback costs 0.6 pp of same-direction agreement and removes two-thirds of the contradictions.
- **Removing the self-reference** (v9 CI → v9 LOO): same-direction 73.0 % → 71.0 %, Pathogenic recall 60.8 % → 48.1 %. That 2 pp is the part of v8's number that was self-fulfilling.

Any ClinVar-derived evidence path makes ClinVar concordance partly circular.
Publish the leave-one-out number as the headline and the ClinVar-informed number
as the "what a lab pipeline actually does" comparison.

The SA stack is unchanged from v7. v8 adds only the three ACMG
criterion-execution fixes from the medical-genetics discordance review
(merged in #50) plus a concordance-script fix so the review export is
complete (see below).

v3 was a partial run (PhyloP+SpliceAI loaded but bugs still latent);
its results are functionally indistinguishable from v2 and were
overwritten before being preserved.

## v10: gnomAD v4 QC, hemizygote and filtering-AF columns (Tier B)

v10 is a **data-layer** run, not a criteria-logic run. The gnomAD builder had
never captured nine fields the ACMG frequency criteria need, so B2, B3 and B5
of the round-2 review plan were blocked on a re-extraction rather than on any
classification work. Write-up in
[`docs/ACMG_EXPERT_REVIEW_ROUND2.md`](../../docs/ACMG_EXPERT_REVIEW_ROUND2.md).

### New gnomAD columns

| Column | Source | Used by |
|---|---|---|
| `filterAc0`, `filterVqsr`, `filterInbreeding` | VCF FILTER | BA1/BS1/BS2/PM2 report NotEvaluated on a non-PASS record |
| `lcr`, `segdup` | INFO flags | same, per site rather than per gene |
| `nonPar` | INFO flag | required to read `AC_XY` as hemizygotes |
| `allAcXY`, `allAnXY` | `AC_XY`, `AN_XY` | BS2 counts hemizygotes |
| `faf95`, `faf95Max` | `faf95`, `fafmax_faf95_max` | BA1/BS1 test the filtering AF |

The database is rebuilt as `.osa2`, which is **smaller than the v1 `.osa` it
replaces even with the nine extra columns**: chr22 is 9.3 MB against 16.2 MB for
852,255 records.

`nhomalt_XY` was in the plan and is deliberately **not** extracted. gnomAD calls
XY samples haploid outside the PAR and so never records one as homozygous:
across all 6,955 non-PAR chrX records in the IDS region it is zero, including at
a site with 109,916 XY carriers. `AC_XY` is the hemizygote count.

### Criteria changes

- BA1/BS1/BS2/PM2 share one precondition gate (`criteria::frequency_gate`)
  rather than repeating it. A record gnomAD itself rejected, or a site it flags
  `segdup`/`lcr`, blocks all four: when a frequency cannot be believed, neither
  its presence nor its absence is evidence.
- BS2 counts individuals with no functional copy - homozygotes plus, on a
  non-PAR sex-chromosome site, hemizygotes. Cohort size becomes
  `(AN - AN_XY)/2 + AN_XY`, which is exactly `AN/2` with no XY columns, so
  autosomal results are unchanged.
- BA1/BS1 test the filtering allele frequency (Whiffin 2017) rather than a point
  estimate, from a shared helper so the two cannot disagree.

### Backward compatibility

Every new field is optional and every gate degrades to its previous behaviour
when absent, so upgrading fastVEP without rebuilding the annotation database
cannot change a call. Two switches cover the judgement calls:
`gnomad_region_flags_block_frequency` and `use_filtering_af`.

## Headline metrics per run

|                            |     v1     |     v6     |     v7     |     v8     |  v9 (CI)   |  v9 (LOO)  |
|----------------------------|-----------:|-----------:|-----------:|-----------:|-----------:|-----------:|
| Same-direction concordance |   54.7 %   |   70.3 %   |   70.8 %   |   73.6 %   |   73.0 %   | **71.0 %** |
| Exact match                |   52.7 %   |   56.8 %   |   58.7 %   |   61.3 %   |   60.7 %   | **59.9 %** |
| Opposite direction (count) |     36     |    338     |    311     |    122     |   **41**   |   **46**   |
| NoCall                     |   0.0 %    |   0.0 %    |   0.0 %    |   0.0 %    |   0.1 %    |   0.1 %    |
| Pathogenic recall          |   15.7 %   |   20.6 %   |   63.8 %   |   63.3 %   |   60.8 %   |   48.1 %   |
| Likely_pathogenic recall   |   20.9 %   |   26.7 %   |   51.8 %   |   62.5 %   |   55.5 %   |   31.1 %   |
| VUS recall                 |   96.6 %   |   92.6 %   |   91.5 %   |   94.2 %   | **97.3 %** | **97.3 %** |
| Likely_benign recall       |    3.2 %   |   42.4 %   |   42.4 %   |   48.4 %   |   45.5 %   |   45.5 %   |
| Benign recall              |   33.2 %   |   58.0 %   |   58.0 %   |   61.8 %   |   57.2 %   |   57.2 %   |

**v9 is reported in two modes.** `v9 (CI)` is ClinVar-informed and directly
comparable to v8; `v9 (LOO)` is leave-one-out, where PS1 and PM1 discount the
variant's own ClinVar record. See "The v9 honesty correction" below. The v9
NoCall column is the 382 ClinVar reference-agreement records (`ALT=.`) that v8
and earlier were silently annotating as real variants.

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
