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
- **v8 current: `output_v8/`** (discordance-review fixes; full outputs +
  `discrepancies_for_md_review.tsv` medical-geneticist review table +
  `README.md`; ClinVar 2026-06-27 refresh check noted in that README)

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

## Progression, v9 to v12

| Metric | v9 | v10 | v11 | v12 |
|---|---:|---:|---:|---:|
| Exact match | 59.9 % | 61.2 % | 61.4 % | **61.8 %** |
| Same-direction | 71.0 % | 74.1 % | 75.0 % | **76.0 %** |
| Pathogenic recall | 48.1 % | 58.5 % | 65.0 % | **65.0 %** |
| Benign recall | 57.2 % | 65.3 % | 65.3 % | **68.4 %** |
| False-benign | 18 | 10 | 10 | 14 |
| False-pathogenic | 28 | 15 | 46 | 45 |
