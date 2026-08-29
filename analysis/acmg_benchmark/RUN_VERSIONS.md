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

## v10: gnomAD v4 QC columns, the gene-level gates, and two measured thresholds

v10 closes every item of the round-2 medical-genetics review except B1, which
needs a curated table rather than code. Write-up in
[`docs/ACMG_EXPERT_REVIEW_ROUND2.md`](../../docs/ACMG_EXPERT_REVIEW_ROUND2.md).

### Headline

Against v9 leave-one-out, the like-for-like mode (both exclude the variant's own
ClinVar record):

| Metric | v9 (LOO) | **v10** | change |
|---|---:|---:|---:|
| Exact match | 59.9 % | **61.2 %** | +1.3 pp |
| Same-direction | 71.0 % | **74.1 %** | +3.1 pp |
| **Opposite-direction** | **46** | **25** | **−46 %** |
| Pathogenic recall | 48.1 % | **58.5 %** | +10.4 pp |
| Benign recall | 57.2 % | **65.3 %** | +8.1 pp |
| VUS recall | 97.3 % | 97.0 % | −0.3 pp |
| NoCall | 382 | 382 | — |

Every headline metric improved except VUS recall, by 0.3 pp.

### What is in it

The data layer (B2, B3, B5) plus the gene-level gates (B6, B7), the curated
functional-evidence input (B8), the two Tier C leftovers (C5, C6), and two
thresholds moved from convention to measurement (BS2 prevalence, PM2 for
dominant genes).

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

## Threshold sweeps

Thresholds in the classifier are guideline values wherever a guideline states one
(REVEL bands from Pejaver 2022, SpliceAI from Walker 2023, BA1 at 5%). Where no
guideline states one, the value should come from measurement rather than from
convention. `sweep_bs2_thresholds.py` runs the full benchmark once per setting
(78 s each) and reports what each choice costs.

### BS2 prevalence bar (swept 2026-08-12)

| bar | BS2 fires | false-benign | correct benign | marginal cost per false-benign avoided |
|---|---:|---:|---:|---|
| 1e-6 | 63,285 | 54 | 139,270 | |
| 1e-5 (old default) | 51,875 | 45 | 136,014 | 362 |
| 5e-5 | 42,219 | 41 | 133,949 | 516 |
| 1e-4 | 37,808 | 40 | 133,407 | 542 |
| **1e-3 (new default)** | **27,290** | **38** | **132,815** | **296** |

No knee, so the data fixes the exchange rate but not the choice. The bar is a
maximum credible disease prevalence, so it has to cover the most prevalent
Mendelian conditions rather than the typical one; hearing loss, alpha-1
antitrypsin deficiency and familial Mediterranean fever sit near 1 in 1,000. A
false-benign call is a missed diagnosis and a lost benign call is a VUS, so the
asymmetry favours the safer bar, and the last step is also the cheapest.
Reasoning in full in [`docs/ACMG.md`](../../docs/ACMG.md#choosing-the-bs2-prevalence-bar).

Two knobs measured at the same time and found inert:

- `bs2_ad_min_ac`: 5 to 100 removes 45 BS2 firings on pathogenic truth and
  changes false-benign and opposite-direction counts by **zero**. Firing is not
  deciding, and an earlier claim here that this was the weakest threshold was
  based on firing counts and was wrong.
- `bs2_ar_min_hom`: 1 versus 10 is identical on every metric. Subsumed by the
  prevalence test; retained only as a floor for configs that lower the bar.

Swept against `sa_db_noclinvar`, which was stable while the gnomAD rebuild ran
and removes the ClinVar-derived criteria, so the fit is less circular. Absolute
numbers therefore differ from the v9 table below; the curve is the point. To be
re-confirmed on v10.

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
- v8: `output_v8/` (discordance-review fixes; full outputs +
  `discrepancies_for_md_review.tsv` medical-geneticist review table +
  `README.md`; ClinVar 2026-06-27 refresh check noted in that README)
- **v14 current: `output_v14/`**, produced by `run_full_benchmark_v14.sh`.
  v14 is the first run whose `clinvar_protein.oga` indexes benign as well as
  pathogenic missense; a run against an older index is still valid but PM1's
  benign-variation half is inert in it.

## v11: SVI point system on the pathogenic side, and B7 corrected

Two changes, both found by mining v10's own discrepancy table rather than from the review list.

| Metric | v10 | **v11** | change |
|---|---:|---:|---:|
| Exact match | 61.2 % | **61.4 %** | +0.2 pp |
| Same-direction | 74.1 % | **75.0 %** | +0.9 pp |
| Pathogenic recall | 58.5 % | **65.0 %** | +6.5 pp |
| Likely-pathogenic recall | 37.9 % | **47.9 %** | +10.0 pp |
| Benign recall | 65.3 % | 65.3 % | — |
| VUS recall | 97.0 % | 96.9 % | −0.1 pp |
| **False-benign (missed diagnoses)** | **10** | **10** | **—** |
| False-pathogenic | 15 | 46 | +31 |
| Opposite-direction (total) | 25 | 56 | +31 |

**The combining rules.** v10 left 2,319 truth-pathogenic variants in VUS on `PVS1` and nothing else -
a lone Very Strong criterion, which is 8 points and inside Likely Pathogenic under the ClinGen SVI
point system, but matches no row of the Richards 2015 table. fastVEP now scores the pathogenic side
by points and keeps the 2015 table for the benign side; the measurement behind that split is in
[`docs/ACMG.md`](../../docs/ACMG.md#which-combining-rules-apply).

**B7 corrected.** The gene-disease validity gate blocked on absence from ClinGen GDV, which cost
1,497 truth-pathogenic PVS1 firings in genes ClinGen has not curated (SPAST, ABCB11, FLG, LAMB3).
It now blocks only where every ClinGen curation of the gene is Limited/Disputed/Refuted/No Known
Disease Relationship - positive evidence rather than a gap in coverage.

**Read the error columns separately.** The two directions did not move together. False-benign - a
missed diagnosis - is unchanged at 10. The whole increase is false-pathogenic, 15 to 46, and it is
the direct consequence of a lone PVS1 now reaching Likely Pathogenic. The trade is roughly 5,200
additional correct pathogenic calls for 31 additional false-pathogenic ones, in the direction that
gets scrutinised rather than filed away.

## v12: BS1 threshold from measurement

One change: `bs1_af_threshold` from 1 % to 0.5 %.

Richards 2015 words BS1 as "allele frequency is greater than expected for disorder", a
per-disease quantity; 1 % was a placeholder for it rather than a derivation, and far too
permissive for a dominant early-onset condition. The reviewer's MTR (7,000 hets), KMT2C (1,721),
TP63 (57), CHD7 (30) and PTCH1 (15) notes all turn on this.

| Metric | v11 | **v12** |
|---|---:|---:|
| Exact match | 61.4 % | **61.8 %** |
| Same-direction | 75.0 % | **76.0 %** |
| Benign recall | 65.3 % | **68.4 %** |
| Likely-benign recall | 45.7 % | **46.8 %** |
| Pathogenic recall | 65.0 % | 65.0 % |
| False-benign | 10 | 14 |
| False-pathogenic | 46 | 45 |

The sweep that picked 0.5 % ran on a 1-in-10 sample, where it looked free - benign recall up,
neither error count moving. The full set says otherwise: +3.1 pp of benign recall for 4 additional
missed diagnoses, roughly 1,200 to 1. A good trade by the same exchange-rate reasoning used for the
BS2 prevalence bar, but a trade, and a sample of that size cannot resolve a four-count difference.
Recorded here because the sample result was quoted before it was confirmed.

## v13: `--acmg` did not imply `--hgvs`, and the PVS1 tree got its missing inputs

The largest single item is a plumbing bug rather than a criterion.

**`--acmg` now turns HGVS on.**
Three criteria read HGVS `c.` notation rather than raw coordinates: PVS1's canonical ±1/±2 splice gate, BP7's Walker 2023 intronic extension, and BA1's Ghosh 2018 exception list, which is keyed by `c.` string.
Without `--hgvs` all three saw an empty string and quietly fell back to weaker behaviour, and no benchmark run from v9 to v12 passed the flag.
The library path was unaffected - `AnnotationContext` already defaults `hgvs: true` - so this was a CLI-only gap, which is why it survived four review rounds.
The output field list is unchanged either way; HGVSc and HGVSp were always in the CSQ header and were simply empty.

**The PVS1 decision tree now has real inputs.**
`protein_truncation_pct` is derived from the CDS position and the transcript's peptide length, and `in_critical_region` from ClinVar pathogenic entries downstream of the truncation point.
Before this, every graded branch that needed a magnitude fell through to a single fallback.
PVS1 grading on truth-pathogenic variants went from `PVS1_Moderate` 3,528 / `PVS1_Strong` 0 to `PVS1_Moderate` 359 / `PVS1_Strong` 2,826.

**BP7 no longer treats a missing SpliceAI score as a prediction of no impact.**
Richards 2015 asks that "splicing prediction algorithms predict no impact"; the absence of a prediction is not one.
This is correct in principle and it changed nothing measurable, which is worth recording: the deep-intronic false-benign calls it was meant to explain all *have* SpliceAI scores, reading 0.00 to 0.14.
That diagnosis was wrong and the real cause is in v14 below.

| Metric | v12 | **v13** |
|---|---:|---:|
| Exact match | 61.8 % | **62.7 %** |
| Same-direction | 76.0 % | **77.6 %** |
| Benign recall | 68.4 % | **71.6 %** |
| Likely-benign recall | 46.8 % | **51.8 %** |
| Pathogenic recall | 65.0 % | 65.0 % |
| False-benign | 14 | **50** |
| False-pathogenic | 45 | 44 |

**Read the false-benign column.** The whole benign-side gain is BP7's intronic extension switching on for the first time, and so is the whole cost: 44,803 additional correct BP7 firings on benign-truth variants against 36 additional missed diagnoses, every one of them deep-intronic.
That ratio is fine in aggregate and terrible in the tail, which is what v14 addresses.

Two changes shipped in v13 but off by default, both because measurement said so:

- **`pvs1_nmd_50nt_rule`.** The exact 50-nucleotide NMD rule Abou Tayoun 2018 states, rather than the last-exon proxy. On a 4,000-variant sample it moved 58 calls from Likely Pathogenic to VUS - 54 of them ClinVar Pathogenic - so it gives up roughly 825 correct pathogenic calls to correct one false-pathogenic. Shipped, documented, default `false`.
- **A tempering rule for lone PVS1 where gnomAD frequency was withheld as unreliable.** Implemented, measured, and **deleted**: 307 correct pathogenic calls lost against 7 corrected. The premise was wrong - `AC0` and `AS_VQSR` are the normal state of a genuinely rare pathogenic variant, not a mismapping marker.

## v14: a far boundary for BP7, and the other half of PM1

Two changes, both taking a criterion from half-implemented to complete. Strictly dominant: 18 opposite-direction calls fixed, none introduced.

**BP7's intronic extension gets a far boundary at 300 nt** (`bp7_max_intron_offset`).
Walker 2023 gives the extension a near boundary (donor +7, acceptor -21) and no far one, which is coherent for a curator working one variant at a time and not for a pipeline applying it to every intronic position in the genome.
Out in the deep intron the only evidence BP7 has is a SpliceAI score in the pseudoexon-activation regime where SpliceAI is weakest - and pseudoexon activation is essentially the only way a deep-intronic variant is pathogenic.
On all 36 deep-intronic ClinVar-pathogenic variants v13 called benign, SpliceAI is present and reads 0.00 to 0.14: the criterion is not missing its evidence out there, it is being misled by it.
The boundary is set by sweep, at the knee where the marginal cost jumps from 31 correct benign calls per diagnosis recovered to 155; the table is in [`docs/ACMG.md`](../../docs/ACMG.md#choosing-bp7s-far-boundary).
The same sweep answers the question the bound raises - whether to keep Walker's extension at all - and the answer is yes: turning it off drops missed diagnoses to 14 but gives up 3.8 pp of benign recall, about 570 correct benign calls per diagnosis, which is stricter than either BS1 or BS2 is held to here.

**PM1 now tests both halves of its own definition.**
Richards 2015 asks for a mutational hot spot or critical domain "without benign variation", and fastVEP counted only the pathogenic half - because `clinvar_protein.oga` only ever indexed pathogenic missense.
The index now carries both directions, marked with `benignIndexed` so an older file is skipped rather than read as "no benign variation here", and a residue ClinVar asserts in both directions is dropped as conflicting.
The same leave-one-out that PS1 uses applies: a variant ClinVar calls benign is not allowed to be its own benign neighbour.
Controls hold - TP53 p.Arg248 keeps PM1 on 23 pathogenic and 0 benign neighbours - while MSH2 p.Gly315Val loses it and falls to VUS, which is the reviewer's cspec GN137 point ("PM1 doesn't apply for MSH2").

| Metric | v13 | **v14** |
|---|---:|---:|
| Exact match | 62.7 % | 62.7 % |
| Same-direction | 77.6 % | 77.5 % |
| **False-benign** | **50** | **33** |
| False-pathogenic | 44 | **43** |
| **Opposite-direction (total)** | **94** | **76** |
| Pathogenic recall | 65.0 % | 64.6 % |
| Benign recall | 71.6 % | 71.5 % |
| VUS recall | 96.7 % | **96.8 %** |

**What each change cost.**
The BP7 boundary is close to free: 17 missed diagnoses recovered for 530 correct benign calls, 0.1 pp of benign recall.
PM1's benign half is not: it removes 1,523 PM1 firings on truth-pathogenic variants and costs 0.4 pp of pathogenic recall, buying one false-pathogenic call and a criterion that matches its own definition.
That is a deliberate choice of correctness over the headline number, and it is the only change in four rounds where those two point in opposite directions.

**One thing the reviewer's rows did *not* get fixed by this.**
The PM1 change resolves MSH2, but not CHD7 or PTCH1, where the objection was gene-level ("variants are spread out in the gene", "no clear hotspots for PTCH1") rather than about benign variation in a ±5-residue window.
Both still have 3 pathogenic and 0 benign neighbours there. A gene-level hotspot-density test would be needed, and no guideline text asks for one.

## v15: BP3 switched on

BP3 has been implemented since v1 and inert since v1, for want of a data file. It reads a
positional interval database and none was ever built.

The source is the UCSC RepeatMasker track - public, no registration, 155 MB - converted by
[`repeatmasker_to_bed.py`](scripts/sa_sources/repeatmasker_to_bed.py) and built as
`custom_bed --name repeatmasker`. 5.3 M primary-assembly intervals, ~240 MB as `.osi`.
`download_sa_sources.sh` now fetches it alongside everything else.

**What it bought: essentially nothing on the headline, and that is the honest headline.**

| Metric | v14 | **v15** |
|---|---:|---:|
| Exact match | 62.67 % | 62.67 % |
| Same-direction | 77.50 % | 77.50 % |
| Benign recall | 71.47 % | **71.49 %** |
| Likely-benign recall | 51.57 % | **51.58 %** |
| False-benign | 33 | 33 |
| False-pathogenic | 43 | 43 |
| Opposite-direction | 76 | 76 |

BP3 fires on 1,307 variants, split 688 on benign-truth against 56 on pathogenic-truth -
directionally correct at roughly 12:1. It moves benign recall by 0.02 pp and changes no
error count at all, which is what a lone Supporting criterion should do: at 1 point it
cannot reach a classification by itself, and the combining rules require it to be joined by
other benign evidence. The value is that a curator reading the evidence now sees the repeat
context, not that the aggregate moved.

**A defect the data exposed.** An interval source only produces an annotation when the
position falls inside an interval, so "this variant is not in a repeat" and "no repeat
database is loaded" reached the classifier identically. BP3 therefore answered *every*
in-frame indel outside a repeat with "repeat region data not available (load RepeatMasker
.osi)" - a setup error, where the truth was that it had looked and found nothing. The
pipeline now tells the classifier whether a repeat source is loaded, so the three states are
distinct. Same bug class as BP7 reading a missing SpliceAI score as a prediction of no
impact, fixed in v13.

Verified that the fix changes messaging only: re-running v15 with and without it produces
byte-identical concordance matrices.

## v16: PS1's splice path, and the pile-up at nine points

PS1's splice track (Walker 2023, PMID 37352859) had never fired.
`same_splice_position_pathogenic` was a field nothing populated, so every canonical ±1/2 splice variant got `evaluated: false` on PS1.
Unlike BP3 in v15, the missing piece was not a download: the evidence is ClinVar, which the benchmark has always loaded.
What was missing was an index keyed the way the criterion asks.

`clinvar_protein.oga` now carries a second index alongside the protein one: every ClinVar (Likely) Pathogenic variant whose HGVS `c.` token puts it on a canonical splice dinucleotide, keyed by genomic position and carrying its REF/ALT and signed offset.
The file grew 13.4 MB to 17.9 MB and still builds in 11 s from the same `variant_summary.txt.gz`.

**This is the largest single move any criterion has produced in this series.**

| Metric | v15 | **v16** |
|---|---:|---:|
| Exact match | 62.67 % | **63.11 %** |
| Pathogenic exact | 0.17 % | **4.38 %** |
| Likely-pathogenic exact | 47.17 % | 44.23 % |
| Same-direction | 77.50 % | 77.50 % |
| Pathogenic recall | 64.61 % | 64.61 % |
| Benign recall | 71.49 % | 71.49 % |
| False-benign | 33 | 33 |
| False-pathogenic | 43 | 43 |
| Opposite-direction | 76 | 76 |

Nothing on the benign side moved by a single variant, and no new opposite-direction call appeared.
The entire effect is 3,793 variants moving `Likely_pathogenic` to `Pathogenic`, and the signature table shows the mechanism exactly:

```
-3793  PM2_Supporting+PVS1
+3793  PM2_Supporting+PS1_Supporting+PVS1
```

PVS1 is 8 Tavtigian points and PM2_Supporting is 1.
Nine points is Likely pathogenic; Pathogenic needs ten.
Half the ClinVar 2-star+ pathogenic set was sitting on exactly that signature, one point short, which is why adding a Supporting criterion to a subset of it moved 3,793 calls.

Of those 3,793, **3,363 are true Pathogenic** (a gain), 411 are true Likely_pathogenic and 19 are true VUS (both losses, both same-direction).
That is 7.8 right moves for every wrong one.
The 411 are worth stating plainly rather than hiding in the net: on those, ClinVar's submitters landed on LP and fastVEP now says P, and the disagreement is exactly whether Walker 2023's PS1 code was applied.

Across the whole set PS1's splice path fires 4,776 times, of which 4,719 (98.8 %) are on pathogenic-side truth and 11 are on benign-side truth.

### What was implemented, and what was deliberately not

Table 3 of Walker 2023 grades PS1 for splicing on two axes: where the comparison variant sits relative to the variant under assessment, and what PVS1 code the variant already carries.
Only the rows where the table's own prerequisite holds by construction are implemented.

That prerequisite is that "the predicted event of the VUA must precisely match the predicted event of the comparison (Likely) Pathogenic variant".
Two variants on the *same canonical dinucleotide* abolish the same donor or acceptor, so it holds automatically.
For a comparison variant elsewhere in the splice motif it does not, and no index can settle it - those rows stay with curators.

Three restrictions follow from the table and from the subgroup's [published response to feedback](https://clinicalgenome.org/docs/clingen-svi-splicing-subgroup-response-to-feedback/) (22 March 2024):

- **The comparison variant must be `Pathogenic`, not `Likely_pathogenic`.** Table 3 reads `N/A` in the LP column for a canonical-dinucleotide variant under assessment. The subgroup's reason (item 7c): "since it is so easy for a ±1,2 dinucleotide variant to reach likely pathogenic, we placed constraints on using these variants as reference to make sure there actually was clinical evidence informing that pathogenic classification". Whether an LP call rests on real clinical evidence is a curator's judgement.
- **The strength is Supporting, not Strong.** Table 3 row 3. The reduction exists "to prevent overweighting of the VUA compared to the original (Likely) Pathogenic comparison variant". Row 5 raises it to full PS1 when the variant's own PVS1 was downgraded, which is implemented in `reconcile_evidence` and **fired zero times on this dataset** - a downgraded PVS1 and a catalogued same-dinucleotide neighbour did not co-occur once in 673,660 variants.
- **The variant's own ClinVar record never counts.** The index carries REF/ALT, so the self entry is identified exactly rather than guessed at, which is why this does not consult `exclude_self_from_clinvar_evidence` - that knob exists for the protein index, which cannot tell which entry is the variant being classified. Verified end to end on two variants that are the only catalogued Pathogenic allele on their dinucleotide, CYB5A c.130-2A>G and DARS2 c.492+2T>C: both report "the ClinVar splice index holds no other Pathogenic variant on this dinucleotide" rather than firing off themselves.

447 of the 4,776 firings are on canonical splice variants carrying no PVS1 at all - a gene with no established LOF mechanism, typically.
Table 3 has no row for those, so they keep the conservative Supporting.

Combining PS1 with PVS1 is not double counting here, and the subgroup says so directly (item 7b): "if there is a relevant pathogenic variant with the same predicted impact as the variant under assessment, then you can use PS1 as well as PVS1(RNA)".

### Two smaller things the work changed

The splice index is keyed by *genomic* position, where the protein index is assembly-independent.
The builder therefore filters `variant_summary` to one assembly (`--assembly`, default GRCh38, which the flag already carried and the builder previously ignored) and stamps the name into every gene record, so a call's provenance is visible rather than assumed.

A gene whose only classified variants are canonical splice ones has no residue substitution to index, and used to fall out of `clinvar_protein.oga` entirely.
It now gets a record with an empty `proteinVariants`.
That is a real behaviour change beyond PS1: `in_critical_region` for such a gene goes from "cannot tell" to "no downstream pathogenic variant indexed", which moved exactly one variant from `PVS1_Moderate` to `PVS1_Supporting`.

## v17: published VCEP frequency bars, and why they do not ship on by default

fastVEP's BA1, BS1 and PM2 bars are single global numbers, measured across all genes.
A ClinGen Variant Curation Expert Panel has done what a global measurement cannot: looked up one disorder's prevalence, penetrance and allelic heterogeneity, and derived a bar for one gene.
Where a panel has published one it outranks anything measured across everything.

[`build_vcep_thresholds.py`](scripts/build_vcep_thresholds.py) reads them out of the [CSpec Registry](https://cspec.genome.network/cspec/ui/svi/)'s JSON-LD API and emits [`data/vcep_thresholds.toml`](data/vcep_thresholds.toml), loadable with `--acmg-config`.
117 genes, 304 bars.
It also emits [`data/vcep_thresholds_audit.tsv`](data/vcep_thresholds_audit.tsv), which carries the verbatim sentence every number was read out of, plus every rejection and its reason, so the parse can be checked rather than trusted.

`GeneOverride` gained `ba1_af_threshold`, which it had been missing: it could express a per-gene BS1 and PM2 but not a per-gene BA1, and BA1 is where most of the published specification lives.

**The result is a good trade that is not yet a safe default.**

| Metric | v16 | **v17** |
|---|---:|---:|
| Exact match | 63.11 % | **63.37 %** |
| Benign recall | 71.49 % | **72.19 %** |
| Pathogenic recall | 64.61 % | 64.60 % |
| False-benign | 33 | 37 |
| False-pathogenic | 43 | **42** |
| Opposite-direction | 76 | 79 |

2,028 more true-Benign variants are called Benign, and 1,780 more calls land exactly right.
The cost is four new false-benign calls, all four of them a *tightened* VCEP bar firing BA1 on a variant ClinVar calls pathogenic:

| Variant | Gene | Truth | VCEP BA1 | fastVEP default |
|---|---|---|---:|---:|
| 1:173914872 A>T | SERPINC1 | Pathogenic | 0.002 | 0.05 |
| 11:22262138 G>A | ANO5 | Pathogenic | 0.003 | 0.05 |
| 15:42403721 C>G | CAPN3 | Pathogenic | 0.005 | 0.05 |
| 20:44413714 C>T | HNF4A | Likely pathogenic | 0.0001 | 0.05 |

### The half of the specification that is not machine-readable

**Two of those four genes carry a founder-variant caveat in the very sentence the bar was read from.** ANO5 and CAPN3 are both limb-girdle muscular dystrophy genes with well-known common founder alleles, and the panels wrote the exclusion into the rule text.
32 of the 286 usable bars, across 18 genes, say something like "and variant is excluded as founder pathogenic variant" - MLH1, MSH2, MSH6, PMS2, BRCA1, BRCA2, TP53, LDLR, GALT, five sarcoglycans, DYSF, FOXN1, RMRP, ANO5, CAPN3.

fastVEP has the mechanism for that already: `ba1_exceptions`, the Ghosh 2018 list.
What it does not have is the exception *contents*, because the panels publish them as prose and, in several cases, as an appendix behind a link.
Applying a tightened bar without its exceptions is applying half a specification, and the half left out is the half that protects pathogenic founder variants.

The other two are different and worth separating.
SERPINC1 and HNF4A carry no caveat: their bar and ClinVar's classification simply disagree, and HNF4A/MODY1 is a textbook reduced-penetrance disorder where a genuinely pathogenic allele is genuinely common.
`clinvar_low_penetrance_blocks_benign_frequency` exists for exactly this and did not catch it, because ClinVar does not carry a low-penetrance term on that record.

**So the table ships as an opt-in file rather than as a default.** Turning it on is one flag and buys 2,028 correct benign calls for four missed pathogenic ones - a ratio of about 500:1, far better than BP7's far-boundary trade. But three of those four are avoidable rather than intrinsic, and the right order of work is exceptions first, default second.

### The gene-disease attribute table

[`build_gene_disease_attributes.py`](scripts/build_gene_disease_attributes.py) assembles the rest of B1 from Orphanet, joined to the VCEP bars: [`data/gene_disease_attributes.tsv`](data/gene_disease_attributes.tsv), 6,857 gene-disorder pairs across 4,199 genes.

| Field | Coverage | Source |
|---|---:|---|
| Inheritance | 96.2 % | Orphanet `en_product9_ages` |
| Age of onset (earliest recorded) | 91.4 % | Orphanet `en_product9_ages` |
| Prevalence class | 84.8 % | Orphanet `en_product9_prev` |
| Prevalence as a numeric upper bound | 61.8 % | derived from the class |
| Published BA1 / BS1 / PM2 | 4.5 % / 4.2 % / 2.7 % | ClinGen CSpec |
| **Penetrance** | **0 %** | **not sourced anywhere** |

Penetrance is the field BS2 most wants - its precondition is "full penetrance expected at an early age" - and no public resource publishes it per gene-disorder pair in machine-readable form.
It is emitted blank rather than guessed at.
The same goes for allelic contribution and genetic heterogeneity, the other two inputs to a Whiffin 2017 maximum-credible-AF calculation: VCEPs compute them per gene and publish only the resulting bar.

### What the extractor refuses to do

304 bars parsed; 107 criterion entries were rejected rather than guessed at, each with its reason in the audit table: 37 PM2 rules that specify absence rather than a bar, 27 that state separate bars for dominant and recessive disease (which one per-gene number cannot hold), 21 with no inequality attached to a frequency, 13 whose exponent was a superscript that did not survive markup stripping, 4 genuinely contradictory, 2 outside any believable range, 3 marked Not Applicable by the panel.
Three genes are dropped because two panels state different bars for them, and `gene_overrides` is keyed by gene alone.

Two guards earned their place during development, and both stay:

- **`0.1%` parsed as `0.1`, not `0.001`.** Python's regex alternation is leftmost-first, not longest-match, so a decimal alternative placed before the percentage alternative swallowed the digits and dropped the `%`. A thousand-fold error, in the benign direction, silent. The ordering in `NUMBER` is now load-bearing and commented as such.
- **BA1 ≥ BS1 ≥ PM2 is checked per gene, and a violation drops the gene.** This is what caught PM2 bars of 2 % read out of "the most frequent pathogenic variant accounts for no more than 2% of..." across eight cardiomyopathy genes, and BA1 bars of 4 and 6 read out of superscript exponents. Both are blocked upstream now; the check stays for the next mis-parse, which will be one nobody predicted.

## v18: the BA1 exception list had never fired, and the benchmark cannot see it

The reviewer supplied Table 1 of Ghosh 2018 (PMID 30311378), the ClinGen SVI list of nine variants that must not be called standalone-benign on frequency.
All nine were already shipped as defaults, cited to the right paper.

**Not one of them had ever matched.**

The matcher compared the full HGVS string the pipeline emits, `ENST00000357618.10:c.187C>G`, against the bare `c.187C>G` the list holds.
The unit tests passed because they fed the matcher bare `c.` tokens the pipeline does not produce - the same defect shape as `--acmg` not implying `--hgvs` in v13, and BP3's never-built database in v15.
Four of the nine were being called Benign on frequency alone:

| Variant | Was | Now |
|---|---|---|
| HFE c.187C>G (p.His63Asp) | Benign | VUS |
| MEFV c.1105C>T (p.Pro369Ser) | Benign | VUS |
| PIBF1 c.1214G>A (p.Arg405Gln) | Benign | Likely benign |
| ACAD9 c.-44_-41dupTAAG | Benign | Likely benign |

Entries are now keyed on **genomic coordinate**, with HGVS as a fallback and the ClinVar variation IDs carried for provenance.
Coordinate matching is not belt-and-braces here, it is required: two of the nine cannot be matched on their published HGVS at all.
BTD is `c.1330G>C` on NM_000060.4 and `c.1270G>C` on the Ensembl canonical, and ACAD9's `c.-44_-41dupTAAG` is the same variant fastVEP spells `c.-45_-44insTAAG`.

### The benchmark is structurally blind to this

| Metric | v16 | **v18** |
|---|---:|---:|
| Every cell of the concordance matrix | | **identical** |

Not one call changed across 673,660 variants, and that is the explanation for how a dead list survived this long rather than a disappointing result.

Seven of the nine are `Conflicting classifications of pathogenicity` in ClinVar, so they carry no consensus label and never enter the truth set.
Of the two that are scored - GJB2 p.Val37Ile and ACADS p.Arg171Trp - neither call moves.
The exception list exists precisely for variants where expert opinion is divided, or where a variant is common *and* pathogenic, which is exactly the population a consensus-truth benchmark discards.

**No aggregate metric computed against ClinVar consensus can validate this criterion.** Only case review can, which is why the nine are now covered by an end-to-end test that annotates all of them through the real binary rather than by a unit test over synthetic inputs.

### What the reviewer's own rows show

Her review flagged GJB2, ACADS, MEFV and BTD - and in every case a *different allele* from the one Ghosh lists.
The mechanism now works and is empty of the entries her cases need.
There is a second reason it is thin for founder alleles: BA1 and BS1 test gnomAD's **grpmax** filtering allele frequency, and gnomAD deliberately excludes the Finnish, Ashkenazi and remaining groups from grpmax because they are bottlenecked.
A Finnish or Ashkenazi founder allele is therefore invisible to that statistic by construction, which is why Ghosh had to list ACADS and BTD by hand - both annotated in the source table as ">5 % MAF only in Finnish".

## v19: a data-quality veto that closed only the benign half of the ledger

The round-5 medical-genetics review is about frequency evidence: which disorders it may be applied to, and which alleles it may be applied to.
Working through the 76 opposite-direction calls to answer that turned up a defect nobody had named, in 21 of the 76.

When gnomAD rejects a site - its own FILTER, a low-complexity tract, a segmental duplication - the classifier withholds BA1, BS1, BS2 **and** PM2.
That is even-handed in intention and badly one-sided in effect.
The benign side loses up to Standalone plus two Strong criteria; the pathogenic side loses one Supporting; and everything left standing, PVS1 above all, is read from the very alignments just declared untrustworthy and is kept at full strength.

The result was a class of confident calls made with the benign half of the ledger closed by the pipeline itself:

| Variant | gnomAD | Veto | ClinVar | fastVEP v18 |
|---|---:|---|---|---|
| RAI1 `c.840del` | AF 0.326, **48,739 hom** | low-complexity region | Benign | Likely pathogenic |
| RAI1 `c.837_838del` | AF 0.387, 10,183 hom | low-complexity region | Benign | Likely pathogenic |
| GIGYF2 `c.3630_3631insG` | AF 0.120, 1,775 hom | low-complexity region | Benign | Likely pathogenic |
| ANAPC1 `c.1393C>T` | AF 0.282 | gnomAD FILTER=AS_VQSR | Benign | Likely pathogenic |
| DIAPH1 `c.3149-1G>T` | AF 0.071 | gnomAD FILTER=AS_VQSR | Likely benign | Likely pathogenic |

A lone Likely pathogenic call now becomes Uncertain when the frequency criteria were withheld for a **site-level** reason, with the veto named in the summary.
Two distinctions in that sentence are load-bearing and both were established by measurement.

**Site-level, not gene-level.** A gnomAD FILTER or a region flag says the alignments here are wrong, which undermines the frameshift call as much as the frequency.
Curated homology (Mandelker 2016) and a ClinVar low-penetrance label say something narrower - that the *frequency* is confounded - and leave the consequence call intact.

**Only where the suppressed frequency mattered.** Vetoing on the site flag alone cost 2,630 correct pathogenic calls to remove 18 wrong ones, because most flagged sites carry a frequency far too low to have supported any benign criterion, so withholding it cost the benign side nothing.
Requiring that the suppressed frequency would itself have met BA1 or BS1 keeps the rule on the variants it was written for: **14 wrong calls removed for 6 correct ones**.

## v20: the exception list reaches BS1 and BS2, for hypomorphic alleles

The reviewer's sharpest round-5 point is a class the frequency criteria are structurally wrong about rather than a threshold that needs tuning.
Where a variant causes disease only in trans with a rarer null allele, its being common is what the disease model *predicts*, and its homozygotes being unaffected is also what the model predicts.
Reading either as evidence against pathogenicity inverts the mechanism.

`Ba1Exception` now carries a `blocks` field naming which frequency criteria an entry suppresses, defaulting to `["BA1"]`, which is the scope Ghosh 2018 wrote for.
Three hypomorphic alleles ship with `blocks = ["BA1", "BS1", "BS2"]`, all three taken from the discordance table under review:

| Variant | gnomAD AF | Hom | v19 | **v20** |
|---|---:|---:|---|---|
| GAA `c.-32-13T>G`, late-onset Pompe | 5.4e-3 | 23 | Likely benign | **VUS** |
| CFTR `c.1210-11T>G`, the 5T allele | 9.8e-3 | 27 | Likely benign | **VUS** |
| SPTA1 `c.4339-99C>T`, alpha-LELY | 6.6e-3 | 40 | Likely benign | **VUS** |

None reaches 5 %, which is why a BA1-only list never touched them.
None reaches Likely pathogenic either, and should not: the deciding evidence is PM3, the second allele in trans, which no annotator computes from a single variant.

A smaller hole closed on the way. BS1 stands down where BA1 already covers the frequency, but it compared against the global 5 % default rather than `effective_ba1_threshold(gene)`, so in a gene with a looser published bar - ABCA4's BA1 is 0.163 - a 6 % allele collected no frequency evidence in either direction.

### v18 to v20

| Metric | v18 | v19 | **v20** |
|---|---:|---:|---:|
| Exact match | 63.1 % | 63.1 % | **63.1 %** |
| Same-direction | 77.5 % | 77.5 % | **77.5 %** |
| Opposite-direction | 76 | 62 | **59** |
| Called pathogenic, truth benign | 43 | 29 | **29** |
| Called benign, truth pathogenic | 33 | 33 | **30** |

Seventeen opposite-direction errors removed with the headline rates unmoved, which is what a fix to a specific defect should look like.

## Progression, v9 to v14

| Metric | v9 | v10 | v11 | v12 | v13 | v14 |
|---|---:|---:|---:|---:|---:|---:|
| Exact match | 59.9 % | 61.2 % | 61.4 % | 61.8 % | 62.7 % | **62.7 %** |
| Same-direction | 71.0 % | 74.1 % | 75.0 % | 76.0 % | 77.6 % | **77.5 %** |
| Pathogenic recall | 48.1 % | 58.5 % | 65.0 % | 65.0 % | 65.0 % | **64.6 %** |
| Benign recall | 57.2 % | 65.3 % | 65.3 % | 68.4 % | 71.6 % | **71.5 %** |
| False-benign | 18 | 10 | 10 | 14 | 50 | **33** |
| False-pathogenic | 28 | 15 | 46 | 45 | 44 | **43** |

The v9-to-v12 columns are quoted from their own sections above; v13 and v14 are computed
from `concordance_matrix.csv` on the same definitions.

### Is v12 better than v14?

On opposite-direction calls alone, yes: 59 against 76. That column is worth reading before
concluding anything from it, because the difference is not really between versions.

Everything separating them is one criterion, BP7's Walker 2023 intronic extension, which v12
was not running at all - it was inert because `--acmg` did not imply `--hgvs`. Setting
`bp7_max_intron_offset = 0` on v14's code turns it back off, and that configuration beats
v12 on every error count while keeping the rest of v13 and v14:

| Metric | v12 | v14 | v14, BP7 extension off |
|---|---:|---:|---:|
| Exact match | 61.8 % | **62.7 %** | 61.8 % |
| Same-direction | 76.0 % | **77.5 %** | 75.9 % |
| Benign recall | 68.4 % | **71.5 %** | 68.4 % |
| Pathogenic recall | **65.0 %** | 64.6 % | 64.6 % |
| False-benign | **14** | 33 | **14** |
| False-pathogenic | 45 | 43 | **43** |
| Opposite-direction | 59 | 76 | **57** |

So the question is not which run to prefer but whether to apply an SVI recommendation the
tool had never actually been running. The trade and the reasoning are in
[`docs/ACMG.md`](../../docs/ACMG.md#choosing-bp7s-far-boundary); the short version is that
it costs 19 missed diagnoses for ~10,800 correct benign calls, and that switching a
published recommendation off entirely is a deviation from guidance in a way that bounding it
at a measured knee is not.

### What the medical-geneticist review can and cannot settle

Her 122 adjudicated rows cannot arbitrate between v12, v13 and v14. Scored against her own
call on the 74 rows where she stated one, all three versions are **identical**: 7 agree, 56
where fastVEP now says VUS, 11 opposite. Not one row moved between them. Every difference
between these versions is in variants she never saw - 53 % of v12's error set, 70 % of
v13's and 64 % of v14's have been reviewed by nobody.

The directionality is worth recording because it points at what would help next:

| v14 against her call, 74 adjudicated rows | n |
|---|---:|
| fastVEP VUS, she called it pathogenic | 46 |
| fastVEP VUS, she called it benign | 9 |
| fastVEP over-calls pathogenic | 8 |
| Agree | 7 |
| fastVEP over-calls benign | 3 |
| She called VUS, fastVEP definite | 1 |

Among hard errors we over-call pathogenic about 3:1, though 3 of those 8 are the CBS rows
where her colour coding (ClinVar's Benign is correct) contradicts her own note ("not sure
why it is called as Benign in ClinVar ... this insertion is a well studied one"), so read
that as 5 to 8. Among soft errors the ratio runs 5:1 the other way: 46 variants she would
call pathogenic where fastVEP defers to VUS. That gap is not a threshold problem - her notes
cite segregation, in-trans phase and functional data that no annotator computes - which is
why B1 and further adjudication, not further tuning, are what move it.

## v23: the annotation layer, not the criteria

v23 is the first run in this series whose differences come entirely from consequence
calling. No criterion changed, no threshold moved, and the supplementary-annotation stack is
byte-identical to v22. The criteria fire differently because they are reading a different
consequence.

The baseline is v22, which is numerically identical to the documented v20 (v21 and v22 were
re-runs that moved nothing, which is why neither has a section here).

Three annotation fixes, all of them cases where fastVEP disagreed with Ensembl VEP:

- **A splice site is matched against the whole variant, not its first base.** Every splice
  predicate tested only the leftmost base, so a deletion or delins whose *second* base landed
  on the donor or acceptor dinucleotide was missed. BRCA2 `c.9256_9256+1delinsTA` was
  `splice_region_variant`, LOW, where VEP says HIGH. This accounts for almost all the
  movement below.
- **A `delins` is read from its peptides.** A change that both deletes and inserts bases was
  typed from the direction of its length change alone, with `Amino_acids` of the form `E/X`
  and no `Codons` at all. VEP calls most of them `protein_altering_variant`.
- **A codon window is refused when the variant is not contiguous in the CDS.** A change
  straddling a splice junction was translated as though its bases were adjacent in the coding
  sequence.

### v22 to v23

| Metric | v22 | **v23** | change |
|---|---:|---:|---:|
| Exact match | 425,145 | **425,235** | +90 |
| Same-direction | 522,101 | **522,297** | +196 |
| Opposite-direction | 59 | **61** | +2 |
| Exact-match rate | 63.1 % | 63.1 % | - |
| Same-direction rate | 77.5 % | 77.5 % | - |

243 variants changed their call. The rates are unmoved because 243 out of 673,660 does not
reach a decimal place; the interesting number is not the rate but which way the 243 went.

### Scored against real VEP

Every one of the 243 was re-run through Ensembl VEP 115.1 (Docker
`ensemblorg/ensembl-vep:release_115.1`, the same GFF3 and FASTA, `--pick`). VEP picked the
same transcript in all 243 cases, so the comparison is like for like:

| | v22 | **v23** |
|---|---:|---:|
| Changed calls whose IMPACT tier matches VEP | 5 / 243 | **243 / 243** |

That is the whole argument for this run. The calls that moved were calls where fastVEP's own
consequence disagreed with VEP, and they now agree.

| Movement | n | ClinVar 2-star+ says |
|---|---:|---|
| Gained P/LP | 193 | P/LP - recovered true positives |
| Gained P/LP | 9 | VUS |
| Gained P/LP | 3 | B/LB - **new opposite-direction** |
| Lost P/LP | 21 | VUS |
| Lost P/LP | 10 | P/LP - **new false-benign** |
| Lost P/LP | 1 | B/LB |
| LP to P | 3 | P/LP |
| Other | 3 | - |

Both flagged buckets are the same class of variant seen from two sides, and both are
written up for the next medical-genetics round in
[`results/v23/MD_REVIEW.md`](results/v23/MD_REVIEW.md):

- The **10 new false-benign** are all `c.N+2dup` single-base duplications two bases into the
  intron, in MECR, COL11A1, PCDH15, SLX4, BRIP1, SLC7A9, ERCC2, NEB, GLB1 and NPHP3. VEP
  agrees with the new `splice_region_variant` call on all ten: duplicating the base after
  `+2` leaves the `GT` intact. This is the class the round-2 review already reached in the
  *classifier* (PVS1 stands down for an insertion at the site without SpliceAI support);
  what changed is that the consequence itself is now right, rather than PVS1's gate carrying
  the whole burden.
- The **3 new opposite-direction** calls are deletions that genuinely remove part of a donor
  or acceptor dinucleotide - NFKB2, FOXP2, LMX1B - and VEP calls all three HIGH. The
  disagreement is with ClinVar's classification rather than with the annotation. Two sit in
  repeat-rich introns where the deletion 3'-shifts far from the site in HGVS notation, which
  is a plausible reason a submitter read them as benign.

## v26: the same annotation layer, finished

v23 fixed three consequence-calling bugs. v26 finishes the job: it replaces the last of the
hand-rolled annotation logic with a port of Ensembl's own model, and corrects two places
where the *classifier* was misreading the (now more complete) annotation. No criterion logic
changed, no threshold moved, and the supplementary-annotation stack is identical to v23.

Measured against real Ensembl VEP 115.1 over a 6,600-variant stratified ClinVar sample -
150,725 transcript rows, the same GFF3 and FASTA both tools read:

| Field | v23 | **v26** |
|---|---:|---:|
| Consequence terms (coding rows) | 98.0 % | **100 %** |
| `Amino_acids` | 77.7 % | **100 %** |
| `Codons` | 72.5 % | **100 %** |
| HGVSc (coding rows) | 87.7 % | **99.5 %** |
| HGVSp (coding rows) | 83.4 % | **98.8 %** |
| Splice terms (all rows) | 99.00 % | **100 %** |
| Whole consequence set (all rows) | 97.96 % | **99.92 %** |
| IMPACT (all rows) | 99.76 % | **99.93 %** |

### What moved

- **One codon window for every shape.** Four code paths - substitution, pure insertion, pure
  deletion, delins - each had a partial version of `VariationEffect.pm`'s rules. They are now
  one window, translated on both sides, with every term read off the peptide and codon pair.
- **Splice terms from one pass over the introns.** `_intron_effects` tests the parts of a
  variant that *differ*, selects its intron lists from the whole span, sorts the pair for the
  polypyrimidine tests and nowhere else, and assigns rather than or-assigns `splice_region`.
  The extended terms coexist with donor and acceptor rather than hiding behind them, and so
  does `intron_variant`.
- **HGVS applies the 3'-rule before naming anything.** An insertion is shifted before it is
  tested for duplication; an equal-length replacement by the reverse complement is `inv`; a
  span with one end in an exon and the other in an intron is written over both ends.
- **PVS1 reads the offset a span covers, not its endpoint.** `c.764_771+9delins…` reaches the
  canonical dinucleotide on its way; reading `+9` stood PVS1 down on 14 ClinVar-pathogenic
  donor and acceptor deletions. Those spans had no HGVSc at all before, which is why the gate
  only started firing once HGVSc improved.
- **A duplication is not written across an exon boundary.** Offsets do not run through zero.

### v23 to v26

| Metric | v23 | **v26** | change |
|---|---:|---:|---:|
| Exact match | 425,235 | **425,322** | +87 |
| Same-direction | 522,297 | **522,473** | +176 |
| Opposite-direction | 61 | 64 | +3 |
| P/LP recall | 62.22 % | **62.42 %** | +187 |
| False-benign (ClinVar P/LP, fastVEP B/LB) | 30 | 30 | 0 |
| False-pathogenic (ClinVar B/LB, fastVEP P/LP) | 31 | 34 | +3 |

290 calls changed. Every one was re-run through Ensembl VEP 115.1 with `--pick`; VEP chose the
same transcript in all 290.

| | v23 | **v26** |
|---|---:|---:|
| Changed calls whose IMPACT tier matches VEP | 232 / 290 | **276 / 290** |

### The one deliberate disagreement

The 14 that do not match are a single class, and fastVEP disagrees with VEP on purpose.
`ref_eq_alt_sequence`'s first clause asks whether a replacement keeps the residue it starts on
and introduces a terminator *anywhere* - not whether the annotated terminator survived - and
`frameshift_variant` and `stop_gained` both defer to it. Real VEP 115.1 therefore calls BRCA1
`c.5030_5033dup` `inframe_insertion,stop_retained_variant`, MODERATE, where ClinVar 3-star
says Pathogenic; the same for BRCA1 `c.1499_1508dup`, BRCA2 `c.3205_3206ins…`, TP53
`c.895_919dup` and 30 more in this set. fastVEP keeps `stop_gained` and `frameshift_variant`.

That clause accounts for 101 of the 117 remaining consequence-set differences and 369 of the
370 IMPACT differences. Every one is fastVEP reporting HIGH for a variant that truncates the
protein.

### Review queue

Two documents. [`results/v26/MD_REVIEW.md`](results/v26/MD_REVIEW.md) is the v23-to-v26 delta with
the full table in `call_changes_v23_to_v26.tsv`.
[`results/v26/MD_REVIEW_ROUND6.md`](results/v26/MD_REVIEW_ROUND6.md) is the standalone note for the
reviewer, covering **v9 to v26** - the run her last comments were written against through this one -
because round 5's table and note went out and have not come back. It folds in the answers to her six
points, the gnomAD-rejection defect those comments uncovered, both annotation rounds and the two
classifier defects, and carries round 5's 59 rows forward into
`discordant_64_round6.xlsx`, built by `scripts/07_build_round6_table.py`. Three buckets need her:

- **Three insertions at an exon's last base** (MSH6 ×2, DSP) that are frameshifts where the
  VCF puts them and intronic duplications after the 3'-shift. ClinVar calls all three benign.
  PVS1's "a pure insertion may not do what its position suggests" gate exists on the splice
  track and not on the frameshift track; whether it should is a criteria question.
- **One inversion over an acceptor** (NBN, ClinVar Pathogenic) that loses PVS1 because
  Ensembl matches only the bases that differ, and an inversion leaves its palindromic
  positions unchanged. VEP agrees it is MODIFIER.
- The deliberate divergence above, if the reviewer would rather match VEP.

v24 and v25 were intermediate runs of this work. v24 exposed both classifier bugs above -
34 ClinVar-pathogenic calls lost PVS1 in it - and v25 fixed one of the two. Neither is a
shipped state; they have no section here and no results directory.

## v27: the 3'-rule reaches the splice junction, and PVS1 stops firing on an intact donor

No criterion logic changed, no threshold moved, and the supplementary-annotation stack is
identical to v23 and v26.
The binary is `master` through #104 plus the then-open #106.
Every call that moved, moved because an HGVSc string moved.

#104 let the HGVS 3'-shift apply to a span with one end in an exon and the other in an intron.
v26 stopped such a span at the junction; v27 shifts it the way Ensembl does.
541 variants got a different HGVSc, and consequence terms and IMPACT are **byte-identical on all
673,660 variants** - which is also the measurement that says #106 contributed nothing here, since
adding coding terms is the only thing #106 does.

### v26 to v27

| Metric | v26 | **v27** | change |
|---|---:|---:|---:|
| Exact match | 425,322 | **425,334** | +12 |
| Same-direction | 522,473 | **522,474** | +1 |
| Opposite-direction | 64 | **59** | -5 |
| False-benign (ClinVar P/LP, fastVEP B/LB) | 30 | 30 | 0 |
| False-pathogenic (ClinVar B/LB, fastVEP P/LP) | 34 | **29** | -5 |
| P/LP recall | 62.42 % | 62.41 % | -15 variants |

56 calls changed: 34 `LP -> VUS`, 11 `VUS -> LB`, 11 `LB -> VUS`.

### Every one of the 56 was checked against real Ensembl VEP 115.1

Docker `ensemblorg/ensembl-vep:release_115.1`, the same GFF3 and FASTA fastVEP reads, `--pick`.

| HGVSc VEP agrees with | count |
|---|---:|
| v27 | **55** |
| v26 | 0 |
| neither | 1 |

The one VEP matches neither on is VHL `c.341-25_370dup`, a duplication whose span runs from the
intron into the exon.
fastVEP cannot write that span in either run, and it is the same boundary-crossing family already
recorded against #104.

### The 34 that lost PVS1

Each is a deletion that starts in an exon and runs through the donor dinucleotide *as the VCF spells
it*, and is a purely intronic deletion at +3 or later once shifted.
CDH1 is the shape: `c.2438_2439+2del` in v26, `c.2439+5_2439+8del` in v27 and in VEP.

PVS1's splice track reads the offset off the HGVSc, so it stood down on all 34.
**That is the right answer, and v26's PVS1 was a false positive.**
Verified rather than argued: for each of the 34 the maximal 3'-shift was recomputed directly from
the FASTA, independently of both tools.
It lands on the offset VEP and fastVEP report in all 34, every donor is a canonical `GT`, and in
none of them does the shifted span touch +1 or +2.
The shifted spelling is the one that preserves the longest prefix of the reference, so the sequence
through +2 is untouched and the donor survives in the variant transcript.

The consequence stays `splice_donor_variant`, in v27 and in VEP both, because Ensembl does not
3'-shift before calling consequences.
A HIGH consequence next to a VUS call looks like a contradiction and is not one: the two fields
answer questions about different positions, by design.

### The 22 that moved on BP7

BP7's intron-offset gate reads the same shifted offset.
The shift moves a donor-side variant away from the site and an acceptor-side variant toward it, so
11 variants cleared the gate and 11 no longer do.
Both directions follow from the same rule, because the shifted spelling always puts the disturbed
sequence downstream of the shifted position.
VEP agrees with v27 on 21 of the 22.

### Review queue: nothing new

| | v26 | **v27** |
|---|---:|---:|
| Discordant rows for the reviewer | 64 | **59** |
| Resolved since v26 | | 5 |
| **New since v26** | | **0** |

The five that resolved are ClinVar-benign variants fastVEP had been calling likely pathogenic on the
false-positive PVS1 above: AGRN, NFKB2, PTEN, LZTR1, LMX1B.
No variant entered the queue, so round 6 stands as sent and needs no reissue.

One question is worth folding into round 7 when her comments arrive, and it is a criteria question
rather than a defect.
15 of the 34 are ClinVar P/LP that are now VUS.
Their canonical `GT` is intact, so ClinGen SVI's PVS1 splice track does not reach them, but the
deletion still removes `+5`/`+6` of the donor consensus in several.
Whether a near-donor deletion should carry PVS1 at reduced strength, or be left to SpliceAI through
PP3, is hers to rule on.

`scripts/08_diff_calls.py` produces the call-changes table; the v23-to-v26 one was built by hand,
which is why no script existed for it.
It reproduces that table exactly - 290 changed calls - when run over v23 and v26.
