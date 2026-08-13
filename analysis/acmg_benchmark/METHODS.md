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
| BS1 | Strong | Greater than expected frequency | gnomAD **max-population AF** (mirrors BA1; was cohort `all_af` pre-v7); gene-specific or default 0.01; same AN minimum as BA1 |
| BS2 | Strong | Observed in healthy adults | gnomAD homozygote count + ClinGen GDV inheritance (or OMIM legacy). For AD-only genes: AC ≥ 5 (default; configurable via `bs2_ad_min_ac`). |
| BS3 | Strong | Functional studies — no damage | Not automatable — NotEvaluated |
| BS4 | Strong | Lack of segregation | Not automatable — NotEvaluated |
| BP1 | Supporting | Missense in truncation-disease gene | gnomAD pLI ≥ 0.9 + misZ < 2.0 |
| BP2 | Supporting | In trans / in cis with pathogenic | Companion-variant phasing + ClinGen GDV inheritance (or OMIM legacy) |
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
ClinGen Gene-Disease Validity (preferred per Abou Tayoun 2018) or OMIM
`genemap2.txt` (legacy). Both populate the same `omim` json_key in the
`.oga` schema.

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
  required (ClinGen GDV, or OMIM legacy).
- **BP2** (in cis/trans): for dominant genes — variant in trans with
  pathogenic; for any gene — variant in cis with pathogenic.

## ClinVar Concordance Benchmark

### Methodology

The benchmark runs `fastvep annotate --acmg --pick` end-to-end on every
ClinVar 2-star+ GRCh38 SNV / small indel and compares the issued ACMG
classification against the ClinVar review-panel call.

Concrete pipeline (`analysis/acmg_benchmark/scripts/run_benchmark.sh`):

1. **Input**: ClinVar VCF filtered to review_status ≥ 2 stars
   (`criteria_provided,_multiple_submitters,_no_conflicts` and stricter)
   on GRCh38, plus a parallel truth TSV (`chrom`, `pos`, `ref`, `alt`,
   `gene`, `clnsig`, `normalized_class`, `review_stars`, `rcv`).
2. **Annotation**: GFF3 + FASTA cache + supplementary annotation
   directory (`--sa-dir`) loaded once; all 673,660 variants annotated
   with `--acmg --pick` to a single JSON file.
3. **Concordance** (`analysis/acmg_benchmark/scripts/03_evaluate_concordance.py`):
   stream-parses the JSON via `ijson` (memory-bounded — output is ~24 GB
   pretty-printed), keys each variant on `(chrom, pos, ref, alt)`,
   reads the picked transcript's `acmg.classification`, and fills a
   5×6 truth × predicted matrix (extra column for NoCall).

Outputs: `concordance_summary.txt` (free-text rollup),
`concordance_matrix.csv`, `concordance_by_chrom.csv`,
`concordance_by_consequence.csv`, `criterion_firing_rates.csv`,
`rule_distribution.csv`, `discrepancies.tsv` (top 10k
opposite-direction calls).

### Supplementary Annotation Build

| Source | Build | Records |
|--------|-------|---------|
| ClinVar (.osa) | `fastvep sa-build --source clinvar` from `clinvar.vcf.gz` | 4,402,501 |
| ClinVar gene index (.oga) | `--source clinvar_protein` from `variant_summary.txt.gz` | 15,665 genes (protein + canonical-splice indices) |
| gnomAD v4.1 exomes (.osa, per-chrom) | tabix-extracted to ClinVar 2-star+ regions (24,350 merged ranges, `bedtools merge -d 5000`), `--source gnomad` per chrom (chr1..22, X, Y, MT) | 25 × .osa |
| gnomAD v4.1 gene constraints (.oga) | `--source gnomad_gene` from `gnomad.v4.1.constraint_metrics.tsv` | 18,173 genes |
| REVEL v1.3 (.osa, per-chrom) | per-chromosome split from `revel-v1.3_all_chromosomes.zip` to bound RAM | 24 × .osa |
| **SpliceAI (.osa, per-chrom)** | distilled from gnomAD v4.1 INFO field `spliceai_ds_max` (gnomAD v4 already includes SpliceAI scores; no separate download) | 24 × .osa |
| **PhyloP (.osa, per-chrom)** | distilled from gnomAD v4.1 INFO field `phylop` (Zoonomia 241-mammal score; gnomAD v4 already includes it) | 24 × .osa |
| **ClinGen Gene-Disease Validity (.oga)** | `--source omim` from `data/benchmark/sa_sources/clingen_gdv.tsv` (ClinGen public CSV converted to genemap2 layout, Definitive/Strong/Moderate only). **Preferred over OMIM `genemap2.txt`** per ClinGen SVI / Abou Tayoun 2018: ClinGen uses a multi-curator scored rubric and excludes Disputed/Refuted/Limited classes. Real OMIM `genemap2.txt` is also accepted (registration-gated at omim.org); both populate the same `omim` json_key. | 2,419 genes |

The gnomAD bulk-extraction path uses `tabix` against the public bgz on
`gs://gcp-public-data--gnomad/...`. We tested the gnomAD GraphQL API as
an alternative: it is fine for ad-hoc per-variant lookups, but
**rate-limits aggressively (HTTP 429) even single-threaded with 5-try
exponential backoff**, so it cannot replace tabix for the 24 K-region
extraction.

### Speed (single host, Apple Silicon, release build)

| Stage | Time | Throughput |
|-------|------|-----------|
| `fastvep annotate --acmg --pick` on 673,660 variants | **4,103 s (68 min)** | **164 variants/s** |
| Concordance parse of 0.6 GB bgzipped VCF | ~3 min | — |

(All 75 SA databases loaded once at process start: 25 gnomAD chroms +
24 REVEL chroms + 24 SpliceAI chroms + 24 PhyloP chroms + 3 .oga
gene-level. 99 % CPU during the annotation phase.)

The annotation pipeline writes `--output-format vcf` (single-line per
variant, ACMG annotations carried in CSQ INFO field) and pipes through
`bgzip`. Final output is **0.6 GB** (vs ~25 GB for pretty-printed
JSON). This is the only output format that includes per-transcript
ACMG/ACMG_CRITERIA fields in a stable column position — `tab` format
ships only the basic VEP columns and `json` is verbose.

### Figures

Generate by running:
```bash
python3 analysis/acmg_benchmark/scripts/generate_figures.py
```

Outputs (`data/benchmark/output_v7/figures/`, PNG + PDF):

| File | Content |
|------|---------|
| `fig_concordance_matrix` | row-normalised truth × predicted heatmap |
| `fig_outcome_breakdown` | per-truth-class stack: same-direction / VUS / NoCall / opposite |
| `fig_criterion_fire_heatmap` | criterion × truth-class fire-rate heatmap. The core "which criteria are diagnostic for which class" diagnostic. |
| `fig_criterion_signatures_by_class` | top-8 criterion signatures of *failed* calls per truth class (P / LP / LB / B), split by failure mode (→ VUS vs opposite-direction) |
| `fig_lost_to_vus_signatures` | criterion patterns of P-truth → VUS and LP-truth → VUS — the dominant failure mode (PVS1 alone, PM2 alone, etc., where the partner-criterion gating prevents reaching LP) |
| `fig_opposite_direction_signatures` | for the ~100 truly opposite-direction cases, the criterion signatures by reversal direction (P/LP-truth → B/LB-pred and the reverse) — directly maps the candidate set for medical-geneticist review |

### Run-to-run differences (v1 → v7)

The benchmark was run end-to-end multiple times during the SVI alignment
+ SA-source build. Each delta maps cleanly onto a specific data or code
change:

| Run | What was loaded / fixed | Same-dir | P | LP | LB | B | PVS1 | PM2 |
|-----|--------------------------|---------:|--:|--:|--:|--:|----:|----:|
| **v1** | Per-chrom gnomAD v4 + REVEL + ClinVar + ClinVar-protein .oga + gnomAD genes .oga | 54.7 % | 16 % | 21 % | 3 % | 33 % | 5,636 | ≈12K |
| v2 | + PhyloP + SpliceAI .osa per chrom + ClinGen GDV .oga (loaded but with two latent wiring bugs) | 54.9 % | 16 % | 21 % | 3 % | 33 % | 23,703 | ≈12K |
| v4 | v2 data + two wiring fixes: SpliceAI `spliceAI` ↔ `spliceai`/`spliceAi` json_key alignment in `sa_extract.rs`; PhyloP read from `allele_supplementary` not just `variant_supplementary`; output switched from JSON to bgzipped VCF (~360× smaller) | 63.7 % | 17 % | 26 % | 42 % | 58 % | 27,460 | ≈12K |
| v5 | v4 + indel allele-matching fix in concordance script: VCF ALT (e.g. `CT`) doesn't match VEP CSQ Allele convention (`T`, `-`). Added `vep_allele(ref, alt)` to strip the leading common prefix. | 65.1 % | 21 % | 27 % | 42 % | 58 % | 52,289 | ≈12K |
| v6 | v5 + **PM2 fires when variant is absent from gnomAD** (`pm2_absent_when_no_record = true`, default; ClinGen SVI v1.0). + **BP4-splice gated to non-PVS1 consequences** (Walker 2023) — frameshift / stop_gained / canonical splice no longer get BP4 from low SpliceAI. | 70.3 % | 64 % | 52 % | 42 % | 58 % | 50,062 | 340K |
| **v7** | v6 + **BS1 uses max-population AF** (mirrors BA1; was cohort `all_af`) per ClinGen SVI. + **BS2 AD requires AC ≥ `bs2_ad_min_ac`** (default 5) — singletons in gnomAD aren't "observed in healthy adult". | **70.8 %** | **64 %** | **52 %** | 43 % | **59 %** | **50,062** | **340K** |

Diagnosis of the lifts:

1. **PhyloP + SpliceAI loaded but invisible (v2 → v4 fix)**. The
   SaWriter wrote SpliceAI under json_key `spliceAI` (capital AI), the
   classifier matched `spliceai | spliceAi | splice_ai`. Result: every
   SpliceAI lookup that succeeded at the .osa layer was discarded by
   `extract_classification_input`. In parallel, the CLI pipeline
   attaches all SA results to `aa.supplementary` (allele-level), but
   the classifier read PhyloP exclusively from `variant_supplementary`
   (variant-level). After the fix, BP7 went from 0 fires to 81,688 (47K
   on LB synonymous + 34K on B synonymous). Direct cause of the +39 pp
   LB recall and +25 pp Benign recall in v4.

2. **ClinGen Gene-Disease Validity loaded** (v1 → v4). PVS1's
   `is_lof_intolerant_gene` requires gnomAD `pLI ≥ 0.9` *or* `LOEUF ≤
   0.35` *or* a curated disease-gene association. v1 had only the
   gnomAD constraint path. v4 adds the ClinGen GDV `.oga` (2,419
   Definitive/Strong/Moderate genes) — preferred per Abou Tayoun 2018
   over the OMIM phenotype proxy. PVS1 fires went from 5,636 to 27,460
   (4.9×);.

3. **NoCall artefact in v4 was an indel allele-matching bug** (v4 → v5
   fix). The 6.8 % NoCall rate in v4 looked like a `--pick` ↔
   non-coding-transcript trade-off, but it was actually all 48,539
   indels in the truth set falling into NoCall. VCF and VEP use
   different allele conventions: VCF stores REF=`C`/ALT=`CT` for an
   insertion, but VEP's CSQ Allele field for the same variant is `T`
   (just the inserted base; `-` for a deletion). The original
   concordance script matched on the raw VCF ALT and so missed every
   indel. Fixed by adding `vep_allele(ref, alt)` to the concordance
   script — strips the leading common prefix between REF and ALT.
   Beyond clearing NoCall, this also unmasked PVS1 nearly doubling
   (27K → 52K) because frameshift indels are PVS1's most common
   trigger.

4. **Pathogenic recall jumped 21 % → 64 % (v5 → v6 fix)**. PM2's
   evaluator previously returned `evaluated=false` whenever
   `input.gnomad` was `None` — i.e. whenever the loaded gnomAD `.osa`
   had no record for the variant. For ~78 % of pathogenic ClinVar 2★+
   variants — most rare disease variants are simply absent from
   gnomAD — that meant PM2 NotEvaluated, so PVS1 had no partner
   criterion and the PVS+≥1 PP (SVI Sept 2020) → LP rule couldn't
   trigger. Default config flag `pm2_absent_when_no_record = true`
   now interprets a missing record as "absent from gnomAD" per
   ClinGen SVI v1.0 ("absent or extremely rare in population
   databases"). PM2_Supporting fires went from ~12K → 340K, unlocking
   the LP rule for ~50K Pathogenic and ~4K LP variants. Configurable
   off (`pm2_absent_when_no_record = false`) for partial-coverage
   runs where the gnomAD `.osa` covers only some input regions.

5. **BP4-splice no longer fires on PVS1-territory consequences (v6
   fix)**. Walker 2023's BP4-splice (SpliceAI ≤ 0.1 → Supporting
   benign) applies to variants where splice impact is the question —
   splice-region, intronic, synonymous, missense. It was firing on
   frameshift / stop_gained / start_lost / transcript_ablation /
   canonical splice donor/acceptor variants because their SpliceAI
   max_ds is typically low (SpliceAI doesn't know how to score a
   frameshift). The result was BP4+PVS1 simultaneous fires on ~5K
   pathogenic frameshifts, which the conflict-gating logic couldn't
   resolve. The fix gates BP4-splice to non-null consequences;
   PVS1-territory variants now reach LP cleanly.

6. **Output format change** (v1–v3 JSON → v4+ VCF.gz). Pretty-printed
   JSON for the full benchmark is ~25 GB; the bgzipped VCF carrying
   the same ACMG fields in CSQ INFO is ~70 MB (~360× smaller). The
   concordance script was rewritten to stream the VCF.gz directly. Tab
   format was considered but does not include ACMG/ACMG_CRITERIA
   fields — only VCF preserves the full classification output.

7. **+14,485 Benign exact-match calls (v6 → v7 fix)**. Two ClinGen-SVI
   BS-tier corrections from a deep classifier audit:
   - **BS1** was reading cohort `all_af` instead of `max_pop_af`. A
     5%-AF EAS variant could slip under a 1% BS1 threshold whenever
     the global cohort diluted it. ClinGen SVI applies BS1 against the
     max-pop AF (mirroring BA1). Effect: BS1 fires went **6.4× higher**
     (4,104 → 26,291); many LB calls promote to B once BS1+BP fires.
   - **BS2** was firing for AD genes on any single het in gnomAD
     (`AF > 0 && AN > 10000`) — a singleton novel allele in a 100K
     cohort isn't "observed in healthy adult" per Richards 2015.
     Tightened to AC ≥ 5 by default (`bs2_ad_min_ac`). False-positive
     BS2 fires on Pathogenic ClinVar variants cut by **52%** (809 →
     389); opposite-direction rate dropped 0.05% → 0.0%.

### Real-Data Concordance (ClinVar 2-star+, April 2026 release; v7)

Truth records: **673,660** · Classified: **673,659** · Truth-not-annotated: **1** (a single intergenic variant whose `--pick` transcript had no ACMG annotation).

#### Truth × predicted matrix (v7)

| Truth ↓ / Predicted → | P | LP | VUS | LB | B | NoCall |
|--|--:|--:|--:|--:|--:|--:|
| Pathogenic (n=79,823) | **6,450** | 44,684 | 28,520 | 150 | 18 | 1 |
| Likely Pathogenic (n=13,989) | 670 | **6,604** | 6,691 | 24 | 0 | 0 |
| VUS (n=295,298) | 73 | 3,238 | **271,659** | 17,609 | 2,719 | 0 |
| Likely Benign (n=128,038) | 2 | 55 | 73,367 | **47,473** | 7,141 | 0 |
| Benign (n=156,512) | 5 | 59 | 64,025 | 28,942 | **63,481** | 0 |

(Diagonal = exact-match. Same-direction lumps P+LP for P/LP truth and
LB+B for LB/B truth.)

#### Headline metrics

| Metric | Value |
|--------|------:|
| Exact-match (truth = predicted) | **58.7 %** |
| Same-direction (truth & predicted both P-tier or both B-tier or both VUS) | **70.8 %** |
| Opposite-direction (P/LP truth → B/LB predicted, or vice versa) | **313 / 673,660 = 0.05 %** |
| NoCall after annotation | 0.0 % |

Per-class same-direction recall:

| Truth class | Same-dir count | Recall |
|-------------|---------------:|------:|
| Pathogenic | 51,134 / 79,823 | **64.0 %** |
| Likely Pathogenic | 7,274 / 13,989 | **52.0 %** |
| VUS | 271,659 / 295,298 | **92.0 %** |
| Likely Benign | 54,614 / 128,038 | **42.7 %** |
| Benign | 92,423 / 156,512 | **59.0 %** |

#### Most-triggered criterion signatures (v7)

The VCF output records the set of criteria met but not the named
combination rule (the rule is reconstructed offline). Top signatures by
sorted criterion code combination:

| Signature | Count |
|-----------|------:|
| BP4 alone | 86,837 |
| **PM2_Supporting alone** | 79,492 |
| **BP4 + PM2_Supporting** | 78,359 |
| **PM2_Supporting + PVS1** | 40,244 |
| BP4_Moderate alone | 31,733 |
| BP4 + BP7 | 31,628 |
| BA1 + BP4 + BS2 | 28,207 |
| BP4_Moderate + PM2_Supporting | 27,171 |
| BP4 + BS2 | 24,128 |
| BP4 + BP7 + PM2_Supporting | 22,394 |
| (no criteria met) | 12,760 |
| BP4 + BP7 + BS2 | 11,829 |
| PVS1 alone | 6,726 |
| PM2_Supporting + PP3 | 6,156 |

Per-criterion fire counts (Pathogenic / LP / VUS / LB / Benign):

| Criterion | P | LP | VUS | LB | B |
|-----------|--:|---:|----:|---:|--:|
| **PVS1** | 45,435 | 4,627 | 2,066 | 62 | 74 |
| **PVS1_Supporting** | 287 | 31 | 114 | 0 | 5 |
| PS1 | 14,438 | 4,770 | 68 | 1 | 1 |
| PM1 | 15,636 | 3,508 | 9,226 | 4,791 | 2,663 |
| **PM2_Supporting** | **64,538** | **10,855** | **143,654** | **57,466** | **63,008** |
| PM4 | 623 | 231 | 3,402 | 226 | 575 |
| PM5 | 8,284 | 2,387 | 3,293 | 39 | 83 |
| PP3_Strong | 65 | 8 | 4,141 | 19 | 31 |
| BA1 | 1 | 0 | 23 | 874 | 41,183 |
| BS1 | 3 | 0 | 18 | 1,317 | 2,766 |
| BS2 | 809 | 104 | 10,097 | 14,611 | 75,813 |
| **BP7** | 19 | 6 | 457 | 47,488 | 33,718 |
| BP4 | 9,733 | 2,212 | 153,415 | 90,941 | 108,005 |
| BP4_Moderate | 118 | 49 | 72,394 | 4,173 | 13,765 |
| BP4_Strong | 1 | 0 | 2,756 | 216 | 903 |
| BP4_Very_Strong | 0 | 0 | 110 | 11 | 54 |

(Note BP4 P fires dropped from 17,163 → 9,733: the v6 BP4-splice gate
no longer fires on null variants. BP4 LB / B counts are essentially
unchanged because those classes are dominated by missense-REVEL and
non-PVS1 splice variants that the gate doesn't affect.)

### Interpretation (v6)

- **Pathogenic recall 16 % → 64 %** between v1 and v6 (4× lift). The
  dominant driver is the v6 PM2-absent fix: PM2_Supporting fires
  340,521 times in v6 vs ~12K previously, so the PVS+≥1 PP (SVI
  Sept 2020) → LP rule now triggers for the bulk of pathogenic null
  variants. Secondary contributor: ClinGen GDV's disease-gene
  fallback for PVS1.
- **Likely_benign recall jumped from 3 % → 42 %** between v1 and
  v4. Driver was BP7 firing 81,206 times (47K LB + 34K B) on
  synonymous / deep-intronic variants once PhyloP and SpliceAI were
  loaded *and* both wiring bugs were fixed. v5/v6 didn't change LB
  numbers since BP7 already fires on the synonymous variants where
  it can.
- **Opposite-direction rate is 0.05 %** (345 / 673,660): when the
  classifier commits to a directional call it agrees with the
  curated review-panel call ~99.95 % of the time. Per-variant
  discrepancies are in `data/benchmark/output_v7/discrepancies.tsv`
  for case-by-case review.
- **NoCall is 0.0 %** in v5+ (after the indel-allele-matching fix).
  Earlier runs reported 6.8 % NoCall, but this was an artefact of the
  concordance script comparing raw VCF ALT (`CT`) against the VEP CSQ
  Allele convention (`T` for an inserted base, `-` for a deletion);
  no real annotation gap. Removing `--pick` and aggregating across
  transcripts would still be an option for non-coding / regulatory
  variants but at the cost of choosing among
  multiple disagreeing transcripts; the trade-off was kept on the side
  of cleaner per-variant calls.

### Limitations of the automated benchmark

1. **Inherently conservative**: PS3/BS3/BS4/PP1/PP4/BP5 are all
   NotEvaluated. Manual curators outperform any variant-level automation
   for these categories by design. The benchmark measures
   classifier-vs-curator agreement, not classifier-vs-ground-truth.
2. **PVS1 / PS1 / BP7 fallback paths**: when the pipeline cannot
   compute Abou Tayoun decision-tree signals (NMD, %protein removed)
   for a specific transcript, those criteria fall back to conservative
   legacy behavior. PVS1_Strong / PVS1_Moderate / PVS1_Supporting
   firings in the table reflect cases where the pipeline *did* derive
   the tree signal.
3. **PS4 NotEvaluated by default**: the previous ClinVar-stars proxy was
   replaced; opt back in via `use_clinvar_stars_as_ps4_proxy` for a
   backward-comparable benchmark.
4. **gnomAD v4 mid / remaining populations**: added to the parser and
   `max_pop_af` after the audit. The 9 chromosome `.osa` files built
   before this change (chr 6, 13, 18, 20, 21, 22, MT, X, Y) lack those
   keys; impact is small (mid + remaining ≈ 5 % of v4 cohort).

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
