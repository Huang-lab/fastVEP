# Medical-geneticist review queue (v7 benchmark)

This document curates the variants where fastVEP and ClinVar 2-star+
disagree most strongly — i.e. where the classifier commits to a
directional call that contradicts the expert-panel review. It excludes
cases where fastVEP calls VUS due to insufficient automatable evidence
(those need data, not curation).

## Top-line numbers

| Outcome | Count | What it is |
|---------|------:|------------|
| **Truth-not-annotated (NoCall)** | **1** | A single regulatory variant chr6:99593111 G>C (PRDM13, North Carolina macular dystrophy, OMIM 616842.0002). Falls upstream of any transcript so VEP picks `intergenic_variant` and no ACMG block is emitted. Pathogenicity here is via a non-coding regulatory mechanism that the standard ACMG framework doesn't address — flag as "out-of-scope-of-ACMG" and review manually. |
| **Opposite-direction discrepancies** | **101** | Variants where fastVEP calls P/LP and ClinVar says LB/B (or vice versa). The candidate set for MD review. Down from 112 in v6 — the BS2 AC ≥ 5 fix dropped 11 false-positive Pathogenic→LB calls. |
| Same-class differences (P/LP/B/LB ↔ VUS) | 9,899 | These are *not* opposite-direction; they reflect missing manual-curation evidence (PS3 functional, PP1 segregation, etc.) by design. |
| Same-direction calls | 477,104 | Classifier and ClinVar agree on direction. |

The single NoCall is **not** a classifier bug — it's a recognized gap
between coding-variant ACMG and non-coding regulatory pathogenicity.
Recommendation: fastVEP should ship a flag `--report-noncoding` that
emits a stub ACMG block for upstream/3'UTR/intergenic calls so they
appear in downstream review queues rather than vanishing.

## Curated review file

`discrepancies_for_md_review.tsv` ranks all 101 opposite-direction
cases by:

1. **3-star (expert panel) > 2-star** — the truth is more reliable, so
   classifier disagreement is more interesting.
2. **Pathogenic↔Benign (full reversal) > LP↔LB** — extremity weighted
   higher.
3. **Number of criteria fired** — more evidence on the disagreeing
   side = more confident classifier disagreement = stronger candidate
   for re-curation.

Built by `analysis/acmg_benchmark/real_data/04_build_md_review_table.py`
which joins three sources:

- `data/benchmark/output_v7/discrepancies.tsv` — the raw concordance
  output
- `data/benchmark/clinvar_2star.vcf` — for the full ClinVar INFO
- `data/benchmark/output_v7/opposite_direction.fastvep.vcf` — a
  re-annotation of the 112 discrepancy variants with `--hgvs` so the
  HGVSc / HGVSp fields are populated (the full-benchmark VCF.gz omits
  HGVS to stay small)

### TSV columns (43 total)

**Identity / scoring (7 columns)**

| Column | Description |
|--------|-------------|
| `priority_score` | sorted descending; higher = more reviewable |
| `chrom`, `pos`, `ref`, `alt` | GRCh38 coordinates (un-prefixed chrom) |
| `gene` | ClinVar's primary gene symbol |
| `stars` | ClinVar review-status star level (2 or 3) |

**Classification (5 columns)**

| Column | Description |
|--------|-------------|
| `truth_class` | ClinVar normalised class (Pathogenic / Likely_pathogenic / Likely_benign / Benign) |
| `fastvep_class` | fastVEP ACMG-AMP call (P / LP / LB / B) |
| `n_criteria_met` | count of fastVEP criteria that fired |
| `fastvep_met_criteria` | semicolon-separated criterion codes (e.g. `PS1;BS2;BP4`) |
| `consequence_top` | top SO term from VEP picked transcript |

**ClinVar INFO fields (16 columns, prefixed `clinvar_`)**

| Column | Description |
|--------|-------------|
| `clinvar_ALLELEID` | ClinVar Allele identifier |
| `clinvar_CLNHGVS` | top-level HGVS expression (e.g. `NC_000002.12:g.26463969C>G`) |
| `clinvar_CLNDN` | preferred disease name(s), pipe-separated per assertion |
| `clinvar_CLNDISDB` | disease-database tag-value pairs (MedGen, OMIM, Orphanet, MONDO) |
| `clinvar_CLNREVSTAT` | review status (`reviewed_by_expert_panel`, `criteria_provided,_multiple_submitters,_no_conflicts`, etc.) |
| `clinvar_CLNSIG` | significance (`Pathogenic`, `Likely_benign`, etc.) |
| `clinvar_CLNSIGCONF` | conflicting-significance details when present |
| `clinvar_CLNSIGSCV` | SCV accession of the panel call |
| `clinvar_CLNVC` | variant class (single_nucleotide_variant, etc.) |
| `clinvar_CLNVI` | external identifier list (ClinGen CA-id, OMIM, UniProt, etc.) |
| `clinvar_GENEINFO` | gene symbol(s) + Entrez ID at this position |
| `clinvar_MC` | molecular-consequence SO term(s) |
| `clinvar_ORIGIN` | allele origin code (germline / somatic / etc.) |
| `clinvar_AF_EXAC`, `clinvar_AF_TGP`, `clinvar_AF_ESP` | older population AFs from ClinVar's frozen sources |

**fastVEP CSQ fields (14 columns, prefixed `fastvep_`)**

| Column | Description |
|--------|-------------|
| `fastvep_HGVSc` | coding-sequence HGVS on the picked transcript |
| `fastvep_HGVSp` | protein HGVS on the picked transcript |
| `fastvep_BIOTYPE` | transcript biotype (protein_coding, lncRNA, etc.) |
| `fastvep_EXON` / `fastvep_INTRON` | feature/total |
| `fastvep_MANE_SELECT` | MANE Select transcript ID when available |
| `fastvep_CANONICAL` | YES if Ensembl canonical |
| `fastvep_ENSP` | Ensembl protein ID |
| `fastvep_CCDS` | CCDS identifier |
| `fastvep_HGNC_ID` | HGNC numeric ID (empty unless an HGNC xref source is loaded) |
| `fastvep_SIFT` / `fastvep_PolyPhen` | empty unless dbNSFP `.osa` loaded |
| `fastvep_ACMG` | predicted classification short code (P / LP / VUS / LB / B) |
| `fastvep_ACMG_CRITERIA` | `&`-joined met criteria codes |

**fastVEP per-variant score panel (24 columns from JSON, prefixed `fastvep_`)**

These are the actual values that drive the ACMG criteria decisions —
distilled from the same JSON that fastVEP's classifier consumes. Empty
values mean the data source isn't loaded or the variant is absent from it.

| Column | Drives | Source |
|--------|--------|--------|
| `fastvep_revel_score` | PP3 / BP4 (missense) per Pejaver 2022 | REVEL v1.3 `.osa` |
| `fastvep_phylop` | BP7 conservation gate per Walker 2023 | gnomAD v4 INFO `phylop` (Zoonomia 241-mammal score) distilled to PhyloP `.osa` |
| `fastvep_gerp` | unused (placeholder); empty in v7 | — |
| `fastvep_spliceai_dsAg` / `_dsAl` / `_dsDg` / `_dsDl` | acceptor-gain / acceptor-loss / donor-gain / donor-loss delta scores | SpliceAI `.osa` (distilled from gnomAD v4 INFO `spliceai_ds_max`) |
| `fastvep_spliceai_max_ds` | PP3 splice (≥0.2) / BP4 splice (≤0.1) per Walker 2023 | computed = max of the four ds heads |
| `fastvep_spliceai_gene` | gene context for the SpliceAI score | (always `"gnomad"` since we distilled from gnomAD INFO) |
| `fastvep_gnomad_allAf` / `_allAc` / `_allAn` / `_allHc` | PM2 / BA1 / BS1 / BS2 (BS2 is `_allHc > 0` for AD genes) | gnomAD v4.1 exomes `.osa` |
| `fastvep_gnomad_<pop>Af` (afr, amr, asj, eas, fin, mid, nfe, remaining, sas) | per-population AF; drives `max_pop_af` for BA1/BS1 | gnomAD v4.1 |
| `fastvep_gnomad_max_pop_af` | BA1 (>5 %) / BS1 (>1 % default) | computed = max across all populations |
| `fastvep_acmg_classification` | full ACMG call (Pathogenic / Likely_pathogenic / UncertainSignificance / Likely_benign / Benign) | classifier output |
| `fastvep_acmg_triggered_rule` | which Richards 2015 + SVI combination rule fired (e.g. `PVS + ≥1 PP (SVI)` → LP) | classifier output |

**fastVEP gene-level fields (6 columns, prefixed `fastvep_gene_`)**

These come from the `.oga` gene-level annotation databases and inform
the gene-context criteria.

| Column | Drives | Source |
|--------|--------|--------|
| `fastvep_gene_pLI` | PVS1 (≥0.9 → LOF-intolerant), BP1 (missense in pLI≥0.9 truncation gene = sub-threshold) | gnomAD v4.1 constraint metrics `.oga` |
| `fastvep_gene_LOEUF` | PVS1 (≤0.35 → LOF-intolerant) | same |
| `fastvep_gene_misZ` | PP2 (≥3.09 → missense-constrained) | same |
| `fastvep_gene_synZ` | sanity check on misZ calibration | same |
| `fastvep_gene_omim_phenotypes` | PVS1 disease-gene fallback; PM3 / BP2 inheritance (AR / AD); BS2 hom requires AR | ClinGen Gene-Disease Validity `.oga` (preferred per Abou Tayoun 2018) or OMIM `genemap2.txt` (legacy) |
| `fastvep_gene_clinvar_protein_n` | count of pathogenic protein variants registered for this gene; PS1 / PM5 / PM1 driver source | ClinVar `variant_summary.txt.gz` `.oga` |

**Helper (1 column)**

| Column | Description |
|--------|-------------|
| `review_question` | pre-formatted prompt: "Why does fastVEP call X when ClinVar says Y? Inspect: <HGVS>; criteria fired = <list>" |

### Coverage on the 112 discrepancies

Empty cells reflect either (a) the data source doesn't apply to the
variant (REVEL is missense-only; SpliceAI doesn't fire on far-from-
splice variants) or (b) the variant is truly absent from the source
(gnomAD `_allAf` empty = variant not in gnomAD = PM2_Supporting fires
in the classifier per the v7 fix).

| Score | Populated | Note |
|-------|----------:|------|
| REVEL | 60 / 112 | Only missense variants; matches missense subset |
| PhyloP | 74 / 112 | Variants where the position has a gnomAD record (PhyloP is distilled from gnomAD INFO) |
| SpliceAI max_ds | 71 / 112 | Same gnomAD-distillation coverage |
| gnomAD allAf | 59 / 112 | The other 53 are absent from gnomAD — PM2_Supporting fires |

## Top 7 (3-star expert-panel opposite calls)

These are the **highest-priority** cases — the ClinVar review-panel
call is from a Variant Curation Expert Panel.

| # | Gene | Variant | HGVSp | ClinVar (3*) | fastVEP | Criteria fired | Review hypothesis |
|---|------|---------|-------|------------------|---------|----------------|-------------------|
| 1 | **OTOF**  | chr2:26463969 C>G | p.Glu1700Gln | Pathogenic | LB | PS1, BS2, BP4 | Same-AA-as-known-pathogenic + healthy-homozygote conflict. Likely a hypomorphic allele where heterozygotes are healthy but compound-het with another pathogenic allele causes auditory neuropathy. ACMG combiner correctly hits VUS-conflicting; ClinVar's panel has functional/segregation evidence we don't have. |
| 2 | **CYP1B1** | chr2:38075207 C>T | p.Gly61Glu | Pathogenic | LB | PS1, BS2, BP4 | Same pattern as OTOF — known pathogenic-by-AA but observed in healthy homozygotes. Recessive primary congenital glaucoma; this is OMIM 601771.0003. Clinical variability is well-documented. Re-evaluate whether BS2 should be silenced for genes where homozygous healthy adults *can* exist (incomplete penetrance). |
| 3 | **SCN2A**  | chr2:165344715 A>G | p.Lys908Arg | **Benign** | **LP** | PM1, PM2_Supporting, PM5, PP2, BP4 | Five criteria all firing in the pathogenic direction (with BP4 from REVEL ≤ 0.290). ClinVar says benign. **This is the strongest classifier-disagrees-with-panel case.** Worth checking the actual ClinVar submission rationale — possible that the variant is in a benign domain despite the hotspot context. |
| 4 | **MSH2**   | chr2:47416297 G>A | p.Gly315Val | Likely Benign | LP | PS1, PM1, BP4 | Lynch syndrome locus — same-AA-as-pathogenic call. ClinVar has hypothesized hypomorphic / VUS reclassification; check current InSiGHT classification. |
| 5 | **MYOC**   | chr1:171652476 G>A | p.Arg46Ter | Likely Benign | LP | PVS1, PM2_Supporting | Stop_gained in MYOC (juvenile open-angle glaucoma, OMIM 601652). PVS1+PM2 → LP per ClinGen SVI rule. Glaucoma panel may have downgraded based on incomplete penetrance / variable expressivity in healthy elderly. Worth checking gene-specific PVS1 calibration (is MYOC truly haploinsufficient?). |
| 6 | **MSH6**   | chr2:47806652 GTAAC>G | (splice) | Likely Benign | LP | PVS1, PM2_Supporting | Splice donor variant in Lynch-syndrome locus. PVS1 fires on canonical splice. Likely a panel re-classification based on RNA studies showing tolerated transcript or in-frame skip. |
| 7 | **MSH6**   | chr2:47806842 T>TTTGA | p.Lys1358AspfsTer2 | Likely Benign | LP | PVS1, PM2_Supporting | Frameshift variant in MSH6. Same pattern — panel may have downgraded based on functional / segregation evidence. |

## Top 20 (all 2-star+ opposite-direction)

See `discrepancies_for_md_review.tsv` for the full ranked list.
Notable patterns:

- **PS1 + BS2** conflict (~30 % of opposites): a same-AA-as-pathogenic
  match co-occurring with healthy-adult homozygote observations. This
  is the *literal* signature of incomplete-penetrance / hypomorphic
  alleles. ACMG-AMP doesn't have a clean rule for these — manual
  curation needs to choose which evidence dominates.
- **PVS1 + BS2** (~10 %): null-variant pathogenic by mechanism but
  ClinVar has BS2 from healthy adult homozygotes. Likely needs RNA /
  protein functional data the classifier can't see.
- **MSH6, MSH2, MUTYH** appear repeatedly — Lynch-syndrome and
  MAP-related genes with active expert-panel re-curation. These are
  worth a sweep through current InSiGHT classifications.
- **PCSK9** (7 cases) — gain-of-function disease gene where the
  ClinGen Familial Hypercholesterolemia VCEP recently re-curated many
  variants; the v7 calls predate the latest VCEP guidance.

## Recommended workflow

1. Open `discrepancies_for_md_review.tsv` in a spreadsheet.
2. Filter `stars == 3` (top 7 cases) — review first.
3. For each case, look up `rcv` in ClinVar and read the assertion
   evidence (especially PS3 functional and PP1 segregation rationale).
4. If the panel call is well-supported, consider whether fastVEP needs
   a per-gene config override (e.g. silence BS2 for genes with
   well-documented incomplete penetrance, or downgrade PVS1 for genes
   where some null variants are tolerated).
5. If the panel call seems out-of-date, flag for ClinVar re-submission.
6. For the remaining 105 2-star cases, batch-review by gene
   (the file is sorted by gene within priority).

## What's NOT in this list (and why)

- **All 9,899 P/LP → VUS cases**: not a wrong call, just missing
  evidence. Fix is data (PS3 functional, PP1 segregation, PP4
  phenotype-specific) not curation.
- **The 1 NoCall (PRDM13 regulatory variant)**: out of ACMG scope.
- **Any same-direction call**: classifier and panel agree on direction
  (P↔LP or B↔LB) — these are concordant even if not exact-match.
