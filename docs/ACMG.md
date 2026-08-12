# ACMG-AMP Variant Classification in fastVEP

fastVEP implements automated ACMG-AMP variant classification based on the standards published by Richards et al. 2015 (*Genet Med* 17:405-424), with ClinGen Sequence Variant Interpretation (SVI) working group recommendations incorporated.

## Overview

The classifier evaluates all 28 ACMG-AMP evidence criteria for each variant-allele-transcript combination and produces a 5-tier classification:

| Classification | Shorthand | Color (Web UI) |
|---|---|---|
| Pathogenic | P | Red |
| Likely Pathogenic | LP | Orange |
| Uncertain Significance | VUS | Gray |
| Likely Benign | LB | Blue |
| Benign | B | Green |

## Quick Start

### CLI

```bash
# Basic ACMG classification
fastvep annotate \
  --input variants.vcf \
  --gff3 genes.gff3 \
  --fasta ref.fa \
  --sa-dir ./sa/ \
  --acmg \
  --output-format json

# With custom thresholds
fastvep annotate \
  --input variants.vcf \
  --gff3 genes.gff3 \
  --fasta ref.fa \
  --sa-dir ./sa/ \
  --acmg \
  --acmg-config acmg_config.toml \
  --output-format json

# With trio analysis (de novo detection)
fastvep annotate \
  --input trio.vcf \
  --gff3 genes.gff3 \
  --fasta ref.fa \
  --sa-dir ./sa/ \
  --acmg \
  --proband CHILD01 \
  --mother MOTHER01 \
  --father FATHER01 \
  --output-format json
```

### Web UI

1. Check the **ACMG-AMP Classification** checkbox in the options panel
2. Click **Annotate**
3. The results table shows an **ACMG** column with color-coded classification badges
4. Click any badge to view the full evidence detail modal showing all 28 criteria
5. The **Summary** tab includes an ACMG classification distribution chart

## CLI Parameters

| Flag | Description | Default |
|---|---|---|
| `--acmg` | Enable ACMG-AMP classification | disabled |
| `--acmg-config <FILE>` | Path to TOML configuration file for custom thresholds | built-in defaults |
| `--functional-evidence <FILE>` | Curated functional-assay TSV supplying PS3/BS3 | none |
| `--proband <SAMPLE>` | Proband sample name in multi-sample VCF (enables PS2/PM6 de novo detection) | none |
| `--mother <SAMPLE>` | Mother sample name for trio analysis | none |
| `--father <SAMPLE>` | Father sample name for trio analysis | none |

## Configuration File

All thresholds are configurable via a TOML file passed with `--acmg-config`. Any omitted field uses its default value.

```toml
# ── Allele frequency thresholds ──
ba1_af_threshold = 0.05        # BA1: benign standalone (>5% in any population)
bs1_af_threshold = 0.01        # BS1: benign strong (greater than expected for disorder)

# ── PM2 inheritance-aware thresholds (ClinGen SVI v1.0, Sept 2020) ──
pm2_ad_af_threshold = 0.00004  # PM2: AD/unknown — "extremely rare", not strictly absent.
                               # Chosen by measurement and sits inside the range published
                               # VCEP specs use for dominant genes; see "Choosing the PM2
                               # bar for dominant genes" below. Set 0.0 for the literal
                               # Richards 2015 strict-absence reading.
pm2_ar_af_threshold = 0.00007  # PM2: AR — AF ≤ 0.00007 (0.007%)
pm2_af_threshold = 0.0001      # Legacy single-threshold field, retained for back-compat
                               # with pre-PR4 configs; not consulted in the default path.

# ── REVEL thresholds (ClinGen SVI calibrated, Pejaver 2022) ──
pp3_revel_supporting = 0.644   # PP3 Supporting
pp3_revel_moderate = 0.773     # PP3 Moderate
pp3_revel_strong = 0.932       # PP3 Strong
bp4_revel_supporting = 0.290   # BP4 Supporting
bp4_revel_moderate = 0.183     # BP4 Moderate
bp4_revel_strong = 0.016       # BP4 Strong
bp4_revel_very_strong = 0.003  # BP4 Very Strong (only REVEL reaches this band)

# ── SpliceAI thresholds (Walker 2023, ClinGen SVI Splicing Subgroup) ──
spliceai_pathogenic = 0.2      # PP3 Supporting cap (SpliceAI alone never reaches Strong)
spliceai_benign = 0.1          # BP4 Supporting threshold; 0.1–0.2 is uninformative

# ── Conservation thresholds ──
phylop_conserved = 2.0         # PhyloP score above which position is "conserved"

# ── Gene constraint thresholds ──
pli_lof_intolerant = 0.9       # pLI score for LOF-intolerant gene
loeuf_lof_intolerant = 0.35    # LOEUF upper bound for LOF-intolerant gene
pp2_misz_threshold = 3.09      # Missense Z-score threshold for PP2

# ── PM1 hotspot detection ──
pm1_hotspot_window = 5         # Window size (amino acid positions) for hotspot scan
pm1_hotspot_min_pathogenic = 3 # Minimum pathogenic variants in window to call hotspot

# ── gnomAD v4 AN minimum (ClinGen SVI March 2024) ──
min_an_for_frequency_criteria = 2000  # BA1/BS1 require AN ≥ this; below → NotEvaluated

# ── BS2 null-individual test (Richards 2015 "full penetrance ... at an early age") ──
bs2_ar_min_hom = 2                      # Absolute floor on observed individuals with no
                                        # functional copy (homozygotes, plus hemizygotes on
                                        # a non-PAR sex-chromosome site)
bs2_hom_prevalence_threshold = 1e-3     # BS2 fires only when the 95% lower bound on their
                                        # frequency exceeds this prevalence bar, so the
                                        # criterion scales with cohort size rather than with
                                        # a bare count. Chosen by measurement, not convention
                                        # - see "Choosing the BS2 prevalence bar" below.

# NOTE: `bs2_ar_min_hom` is a belt-and-braces floor, not an active knob. A sweep of the
# full benchmark found it inert at the default prevalence bar: values from 1 to 10 give
# identical false-benign and opposite-direction counts, because the prevalence test
# subsumes it. It is kept so that a config lowering the prevalence bar still cannot fire
# BS2 on a single observation.

# ── Genes where population frequencies cannot be trusted (Mandelker 2016, PMID 27228465) ──
# BA1/BS1/BS2/PM2 all report NotEvaluated for these. Specifying the list in TOML REPLACES
# the built-in one (CYP21A2, STRC, HBA1/2, SMN1/2, PMS2, NEB, OTOA, GBA, ...).
# homology_unreliable_genes = ["CYP21A2", "STRC"]

# ── gnomAD record quality (requires a database built with the v4 QC columns) ──
gnomad_region_flags_block_frequency = true  # A site gnomAD flags `segdup` or `lcr` reports
                                            # NotEvaluated for BA1/BS1/BS2/PM2: the per-site
                                            # form of the homologous-gene list above. gnomAD's
                                            # own FILTER verdict (AC0 / AS_VQSR /
                                            # InbreedingCoeff) is never ignored and is not
                                            # covered by this switch.
use_filtering_af = true                     # Test BA1/BS1 against the filtering allele
                                            # frequency (95% CI lower bound, max across
                                            # ancestry groups; Whiffin 2017) rather than the
                                            # population-maximum point estimate
# Both are no-ops against a database built before those columns were extracted.

# ── Gene-level preconditions (require omim.oga / ClinGen GDV) ──
require_gene_disease_validity = true   # PVS1, PP2 and PM1 report NotEvaluated with
                                       # `no_established_gene_disease_relationship` for a gene
                                       # absent from the loaded gene-disease validity source.
                                       # A no-op when no such source is loaded, since "absent
                                       # from a file nobody opened" is not evidence.
mechanism_gates_pvs1 = true            # A curated mechanism that excludes loss of function
                                       # ("GOF", "DOMINANT_NEGATIVE") takes PVS1 to
                                       # NotApplicable, instead of mechanism only ever being
                                       # able to switch PVS1 on (Abou Tayoun 2018)

# ── ClinGen SVI behavior ──
pm2_downgrade_to_supporting = true     # Downgrade PM2 from Moderate to Supporting (SVI)
bp1_max_pathogenic_missense = 3        # BP1 blocked once the gene has this many pathogenic
                                       # missense variants in ClinVar (missense is then an
                                       # established mechanism, whatever the constraint says)
exclude_self_from_clinvar_evidence = true   # PS1/PM1 discount the variant's own ClinVar
                                            # record. Set false for a ClinVar-informed run.
clinvar_low_penetrance_blocks_benign_frequency = true  # BS1/BS2 NotEvaluated when ClinVar
                                            # labels the variant low-penetrance / risk allele
# pp3_max_strength = "Moderate"        # Optional ceiling on PP3 from computational
                                       # evidence alone. Unset (uncapped) by default:
                                       # Pejaver 2022 calibrates REVEL >= 0.932 to
                                       # Strong. Set it if your lab follows the
                                       # stricter "a predictor alone never reaches
                                       # Strong" convention.
use_pp5_bp6 = false                    # Enable PP5/BP6 (disabled by default per SVI)
use_clinvar_stars_as_ps4_proxy = false # Opt back into the legacy ClinVar-stars PS4 proxy
                                       # (true PS4 needs case-control statistics, so off by default)

# ── BA1 exception list (Ghosh 2018, 9 known-pathogenic high-AF variants) ──
# Specifying ba1_exceptions in TOML REPLACES the default 9-variant list.
# Include the defaults plus your additions if you want to retain them.
[[ba1_exceptions]]
gene = "HFE"
hgvs_c = "c.845G>A"
reason = "p.Cys282Tyr — hereditary hemochromatosis"

# ── Trio analysis ──
[trio]
proband = "CHILD01"
mother = "MOTHER01"
father = "FATHER01"
min_depth = 10                 # Minimum read depth for reliable genotype
min_gq = 20                    # Minimum genotype quality

# ── Curated disease mechanisms (shipped table; PCSK9, IFIH1, MYH7, RASopathies, ...) ──
# Specifying gene_mechanisms in TOML REPLACES the shipped table.
# Kept separate from gene_overrides precisely so that setting one [gene_overrides.X]
# block below does not silently discard it.
[gene_mechanisms]
PCSK9 = "GOF"          # LoF lowers LDL — the mechanism of an approved drug class
MYH7 = "DOMINANT_NEGATIVE"   # ClinGen MYH7 spec: PVS1 not applicable
RYR1 = "LOF_and_GOF"   # MH is GoF, the congenital myopathies are LoF → PVS1 still applies

# ── Gene-specific overrides ──
# A `mechanism` here wins over the gene_mechanisms table for that gene.
[gene_overrides.BRCA1]
mechanism = "LOF"
bs1_af_threshold = 0.001
# pm2_af_threshold = 0.00005

[gene_overrides.TP53]
mechanism = "LOF_and_GOF"
disabled_criteria = ["BP1"]

# [gene_overrides.GENE.strength_overrides]
# PM2 = "Moderate"   # Override strength for specific criteria
```

## Choosing the BS2 prevalence bar

`bs2_hom_prevalence_threshold` is the one threshold in the classifier that was set by
measurement rather than by convention, so it is worth showing the working.

BS2 fires on a recessive or X-linked gene only when the 95% lower confidence bound on the
frequency of individuals with no functional copy exceeds this bar.
The bar is therefore a **maximum credible disease prevalence**: the point above which the
observed null individuals cannot all be explained by the disorder itself.

Sweeping it across the full 673,660-variant ClinVar 2-star+ benchmark
([`sweep_acmg_thresholds.py`](../data/benchmark/scripts/sweep_acmg_thresholds.py)) gives:

| bar | BS2 fires | false-benign calls | correct benign calls | marginal cost |
|---|---:|---:|---:|---|
| 1e-6 | 63,285 | 54 | 139,270 | |
| 1e-5 | 51,875 | 45 | 136,014 | 362 correct benign lost per false-benign avoided |
| 5e-5 | 42,219 | 41 | 133,949 | 516 |
| 1e-4 | 37,808 | 40 | 133,407 | 542 |
| **1e-3** | **27,290** | **38** | **132,815** | **296** |

Two things follow.

**The data does not pick a value by itself.** The curve is smooth, with no knee, so every
setting is a defensible point on a trade. What the sweep does establish is the exchange
rate, and that no setting of this or any other BS2 knob drives false-benign calls below 38;
BS2 appears in only 9 of the 45 at the old default, so the residual is not a threshold
problem at all.

**The choice rests on the parameter's meaning and on which error is worse.** As a maximum
credible prevalence, the bar has to cover the most prevalent Mendelian conditions BS2 is
applied to, not the typical one. Hearing loss, alpha-1 antitrypsin deficiency and familial
Mediterranean fever in high-prevalence populations all sit near 1 in 1,000, so a 1e-5 bar is
two orders of magnitude too tight for exactly the disorders a medical geneticist flagged in
review. A false-benign call is a missed diagnosis; a lost benign call becomes a VUS and
costs triage effort. That asymmetry favours the safer bar, and the step to 1e-3 is also the
cheapest on the curve.

Hence the 1e-3 default. It is a plain config value: a lab that weighs VUS burden more
heavily sets a smaller one, and a gene-specific prevalence from a VCEP specification should
override it per gene once that table exists.

Two related knobs were measured at the same time and found inert, which is worth knowing
before spending effort on them:

- `bs2_ad_min_ac` (dominant branch) barely affects the final call. Raising it from 5 to 100
  removes 45 BS2 firings on pathogenic truth but changes false-benign and
  opposite-direction counts by zero. A criterion firing is not a criterion deciding.
- `bs2_ar_min_hom` is fully subsumed by the prevalence test: 1 versus 10 gives identical
  results. It is retained only as a floor for configs that lower the prevalence bar.

## Choosing the PM2 bar for dominant genes

Richards 2015 words PM2 as "absent from controls (or at extremely low frequency if
recessive)", and fastVEP read that literally: `pm2_ad_af_threshold = 0.0`, so a variant in
a dominant or uncharacterised-inheritance gene failed PM2 if gnomAD had seen it even once.

That reading does not survive the change in denominator.
The 2015 text was written against ExAC's 60,706 exomes.
gnomAD v4 carries 730,947 exomes plus 76,215 genomes, so "absent" now means absent from a
cohort more than twelve times larger, which is a far stricter test than the authors were
specifying.
It is also the wrong test: a singleton among 800,000 people **is** the "not seen in the
general population" that PM2 was asking about.

ClinGen SVI moved the same way, downgrading PM2 to Supporting in 2020 because it had been
over-weighted, and the VCEP specifications written since give dominant genes explicit
non-zero bars, clustered between 1e-5 and 1e-4 (the Cardiomyopathy VCEP's MYH7
specification uses 4e-5; the PTEN VCEP uses 1e-5).

Sweeping the key alone over the ClinVar 2-star+ benchmark:

| bar | pathogenic recall | false-pathogenic calls | benign recall |
|---|---:|---:|---:|
| 0 (strict absence) | 37.8 % | 1 | 56.3 % |
| 1e-6 | 46.2 % | 1 | 56.3 % |
| 5e-6 | 51.8 % | 2 | 56.3 % |
| 1e-5 | 54.0 % | 2 | 56.3 % |
| 2e-5 | 55.7 % | 2 | 56.3 % |
| **4e-5** | **56.8 %** | **2** | **56.3 %** |
| 1e-4 | 57.5 % | 2 | 56.3 % |
| 2e-4 | 57.6 % | 3 | 56.3 % |

Unlike the BS2 curve, this one has a knee.
Recall climbs 19 points between strict absence and 4e-5 and then flattens; the last 0.8
points cost a further 2.5-fold loosening.
Benign recall never moves, PM2 being pathogenic-direction only, which is the sanity check
that the knob is doing what it claims.

**4e-5** is the choice: it is at the knee, it costs one additional false-pathogenic call
on the sweep sample, and it sits inside the range VCEPs have actually published for
dominant genes rather than being a number of our own.
Going further would buy a fraction of a point while taking us past any published VCEP
specification.

The literal 2015 reading remains one line away for a lab that wants it
(`pm2_ad_af_threshold = 0.0`).

**Caveat, the same one that applies to BS2.** ClinVar significance is itself assigned
partly from gnomAD frequency, so a frequency threshold tuned against ClinVar is partly
fitted to its own input. The curve's shape is the useful output, not its absolute level.

## Functional evidence (PS3 / BS3)

PS3 and BS3 are the two criteria fastVEP cannot compute.
They ask whether a well-established in vitro or in vivo assay showed the variant to damage, or not damage, gene function - a question answered in a paper, not a database.
fastVEP will not mine literature for them: a wrong PMID in a clinical report is worse than no evidence at all.
What it will do is consume a curated file, which is how every VCEP pipeline works.

```bash
fastvep annotate ... --acmg --functional-evidence functional.tsv
```

The file is a TSV keyed by genomic coordinate:

```text
#chrom  pos       ref  alt  criterion  strength    pmid      note
chr15   88855485  A    G    PS3        Strong      29625052  minigene shows exon skipping
2       47478343  G    A    BS3        Supporting  31391288  normal MMR activity in vitro
```

`chrom` may carry a `chr` prefix or not, and coordinates are the ones in your VCF rather than
fastVEP's normalised alleles, so an entry can be copied straight off the record it describes.
`strength`, `pmid` and `note` are optional; `strength` defaults to Strong.
Blank lines and `#` comments are skipped, so the header above is valid input.
Every other malformed row is a hard error naming the line, because a typo that silently dropped a
row would surface as a variant mysteriously missing its PS3 rather than as a message.

**Strength is a curator's judgement, not a constant.** Brnich et al. 2020 (the ClinGen SVI
functional-evidence recommendation) grades assay strength on validity, controls and dynamic range,
so the column exists for a curator who has read the paper to say Supporting where the assay only
supports Supporting.

**An entry outranks the predictors.** The SVI ordering is explicit: functional evidence is stronger
than computational prediction, and PP3/BP4/BP7 must not be used to argue against sound experimental
data. So a PS3 or BS3 entry also suppresses PP3, BP4 and BP7 for that variant, each with
`Superseded by functional evidence` in its summary. This runs in both directions - an assay showing
no damaging effect outranks a high REVEL exactly as a damaging assay outranks a low one.

The round-2 review's OCA2 case is what this is for. Before, a synonymous OCA2 variant with published
splice-defect data collected BP7 and BP4 from a SpliceAI score of 0.00 and came out Benign:

```text
before   call=B     BA1 + BS2 + BP4 + BP7
after    call=VUS   PS3 (PMID 26637981) + BA1 + BS2, BP4 and BP7 superseded
```

Without an entry, PS3 and BS3 stay `evaluated: false` - no assertion in either direction, which is
the honest answer when nobody has run the assay.

## Required Data Sources

ACMG classification draws on multiple supplementary annotation (SA) sources. Place `.osa`/`.osa2` (allele-level) and `.oga` (gene-level) files in the SA directory:

### Allele-Level Sources (`.osa` / `.osa2`)

| Source | SA Key | Used By | Description |
|---|---|---|---|
| **gnomAD** | `gnomad` | BA1, BS1, BS2, PM2 | Per-population allele frequencies + AN + homozygote counts (BA1/BS1 max-pop AF; BA1/BS1 require AN ≥ 2,000) |
| **ClinVar** | `clinvar` | PP5, BP6 (off by default per SVI); PS4 only when `use_clinvar_stars_as_ps4_proxy = true` | Clinical significance, review status, phenotypes; companion-variant lookups for PM3 / BP2 |
| **REVEL** | `revel` | PP3 (missense), BP4 (missense, including Very Strong band ≤ 0.003) | Missense pathogenicity score (0-1); not applied to non-missense per Pejaver 2022 |
| **SpliceAI** | `spliceai` | PP3 (caps at Supporting per Walker 2023), BP4 (≤ 0.1 → Supporting), BP7 | Splice site delta scores |
| **dbNSFP** | `dbnsfp` | Transparency-only (SIFT / PolyPhen surfaced in `details`) | The pre-PR1 ≥3-of-4 consensus path was removed per Pejaver 2022. |
| **1000 Genomes** | `onekg` | PM2 (supplement) | Population frequencies |
| **TOPMed** | `topmed` | PM2 (supplement) | Population frequencies |

### Positional Sources (`.osa`)

| Source | SA Key | Used By | Description |
|---|---|---|---|
| **PhyloP** | `phylop` | BP7 (conservation tier) | Conservation scores |

### Gene-Level Sources (`.oga`)

| Source | SA Key | Used By | Description |
|---|---|---|---|
| **gnomAD Gene Constraints** | `gnomad_genes` | PVS1, PP2, BP1 | pLI, LOEUF, misZ, synZ |
| **ClinGen GDV** (or OMIM) | `omim` | PVS1, PP2, PM1, BS2, PM3, BP2 | Gene-disease validity, disease associations, inheritance patterns |
| **ClinVar Protein Index** | `clinvar_protein` | PS1, PM1, PM5 | Pathogenic missense by protein position |

### Optional Sources

| Source | SA Key | Used By | Description |
|---|---|---|---|
| **RepeatMasker** | `repeatmasker` | BP3 | Repeat region intervals (`.osi` format) |

### Gene-level (.oga) sources

`fastvep sa-build` supports three gene-level sources, each producing a `.oga` file that the runtime picks up automatically from `--sa-dir`:

```bash
fastvep sa-build --source omim -i genemap2.txt -o sa/omim --assembly GRCh38
fastvep sa-build --source gnomad_genes -i gnomad.v4.1.constraint_metrics.tsv -o sa/gnomad_genes --assembly GRCh38
fastvep sa-build --source clinvar_protein -i clinvar.vcf.gz -o sa/clinvar_protein --assembly GRCh38
```

When a `.oga` is missing, dependent criteria (PVS1, PS1, PM1, PM5, PM3, BP1, BP2, PP2, BS2) degrade gracefully to `evaluated: false` rather than misfiring.
The gene-disease validity gate degrades the other way round, and deliberately: with no `omim.oga` loaded it does not fire at all, because a gene missing from a file nobody opened is not a gene without a disease. Loading the source is therefore what switches the gate on, and PVS1/PP2/PM1 become *stricter* when it is present, not more permissive. See [ACMG_SETUP.md](ACMG_SETUP.md) for download URLs, expected file sizes, and end-to-end verification.

## Evidence Criteria Reference

### Pathogenic Criteria

| Code | Strength | Description | Data Source | Automatable |
|---|---|---|---|---|
| **PVS1** | Very Strong → Supporting* | Null variant in LOF-intolerant gene; graded per Abou Tayoun 2018 decision tree (NMD prediction, %protein removed, critical region, alt start codon, last exon). Two gene-level preconditions run first: the gene must have an established gene-disease relationship, and its curated mechanism must not exclude loss of function | Consequence + pLI/LOEUF/ClinGen GDV + curated mechanism + transcript context | Yes |
| **PS1** | Strong | Same AA change as established pathogenic missense **or** canonical ±1/2 splice predicted to produce same RNA outcome (Walker 2023) | ClinVar protein index + splice catalog | Yes (with .oga) |
| **PS2** | Strong | De novo with confirmed parents | Trio VCF genotypes | Yes (with trio) |
| **PS3** | Strong | Functional studies show damaging. Read from `--functional-evidence`; suppresses PP3/BP4/BP7 for that variant (SVI ranks experimental evidence above prediction) | Curated TSV | Yes (with file) |
| **PS4** | Strong | Prevalence in affected >> controls | Case-control statistics (NotEvaluated by default; ClinVar review-stars proxy is invalid per SVI, opt in via `use_clinvar_stars_as_ps4_proxy`) | No (Partial in proxy mode) |
| **PM1** | Moderate | Mutational hotspot / critical domain. NotEvaluated where the gene has no established gene-disease relationship: a cluster of assertions in an uncurated gene is not evidence of a critical region | ClinVar protein density + ClinGen GDV | Yes (with .oga) |
| **PM2** | Supporting* | Absent/rare in population databases - inheritance-aware: AD/unknown fires at AF ≤ 0.00004, AR at AF ≤ 0.00007 (SVI v1.0 plus published VCEP practice; see "Choosing the PM2 bar for dominant genes"). NotEvaluated wherever BA1/BS1/BS2 are: absence is only evidence when the database could have seen the variant | gnomAD AF + QC flags + OMIM inheritance | Yes |
| **PM3** | Supporting → Very Strong* | In trans with pathogenic (recessive); points-based per SVI v1.0 (in-trans/P=1.0pt, in-trans/LP=0.5pt, unphased/P=0.5pt, unphased/LP=0.25pt, hom=0.5pt cap 1.0) | Phased VCF + ClinVar | Yes (with trio) |
| **PM4** | Moderate | Protein length change (in-frame/stop-loss) | Consequence type | Yes |
| **PM5** | Moderate | Different pathogenic missense at same residue | ClinVar protein index | Yes (with .oga) |
| **PM6** | Moderate | Assumed de novo (partial confirmation) | Partial trio VCF | Yes (with trio) |
| **PP1** | Supporting | Co-segregation in family | Pedigree data | No |
| **PP2** | Supporting | Missense in constrained gene. NotEvaluated where the gene has no established gene-disease relationship: missense constraint measures tolerance to variation, not that missense is how the gene causes disease | Gene misZ score + ClinGen GDV | Yes |
| **PP3** | Supporting-Strong | Computational evidence (deleterious) — REVEL **missense-only** (Pejaver 2022) + SpliceAI ≥ 0.2 caps at Supporting (Walker 2023). Ensemble SIFT/PolyPhen/PhyloP/GERP consensus path removed (Pejaver 2022 endorses single calibrated tool only). | REVEL + SpliceAI | Yes |
| **PP4** | Supporting | Phenotype-specific for single-gene disease | HPO phenotype data | No |
| **PP5** | Supporting | Reputable source reports pathogenic | ClinVar (disabled by default per SVI) | Partial |

*\*PM2 is downgraded from Moderate to Supporting per ClinGen SVI recommendation. PVS1 and PM3 are escalated/de-escalated by graded subcodes (e.g. `PVS1_Strong`, `PM3_Supporting`).*

### Benign Criteria

| Code | Strength | Description | Data Source | Automatable |
|---|---|---|---|---|
| **BA1** | Standalone | AF > 5%, AN ≥ 2,000 (gnomAD v4 SVI March 2024). Tested against the **filtering allele frequency** (95% CI lower bound, max across ancestry groups; Whiffin 2017) where the database provides it, else the population maximum. 9-variant Ghosh 2018 exception list (HFE c.845G>A, MEFV common, BTD c.1330G>C, etc.) blocks BA1 regardless of AF | gnomAD population AFs + FAF + AN + HGVSc | Yes |
| **BS1** | Strong | Cross-population AF > expected, read from the same statistic as BA1 (filtering AF where available, else population maximum). NotEvaluated for homology-confounded genes, for gnomAD-flagged `segdup`/`lcr` sites, for non-PASS gnomAD records, and for ClinVar low-penetrance / risk alleles | gnomAD per-population AFs + FAF + QC flags + ClinVar terms | Yes |
| **BS2** | Strong | Observed in healthy adult - AD/X-linked-D requires AC ≥ 5 (`bs2_ad_min_ac`). AR / X-linked counts individuals with no functional copy (homozygotes **plus hemizygotes** on a non-PAR sex-chromosome site) and requires ≥ `bs2_ar_min_hom` of them **and** a 95% lower bound on their frequency above `bs2_hom_prevalence_threshold` (default 1e-3, set by measurement - see [Choosing the BS2 prevalence bar](#choosing-the-bs2-prevalence-bar)), so the test scales with cohort size. Richards 2015's "full penetrance expected at an early age" qualifier is what this implements: one or two such individuals in a 730 K cohort is what a late-onset or reduced-penetrance disorder looks like, not tolerance | gnomAD hom + hemizygote counts + AN + inheritance | Yes |
| **BS3** | Strong | Functional studies show no damage. Read from `--functional-evidence`; suppresses PP3/BP4/BP7 for that variant | Curated TSV | Yes (with file) |
| **BS4** | Strong | Lack of segregation | Pedigree data | No |
| **BP1** | Supporting | Missense in a gene where primarily truncating variants cause disease. Requires both the constraint signature (high pLI, low misZ) **and** a mutation spectrum without an established missense mechanism (< `bp1_max_pathogenic_missense` pathogenic missense in the ClinVar protein index) | pLI + misZ + ClinVar protein index | Yes |
| **BP2** | Supporting | In trans/cis with pathogenic | Phased VCF + ClinVar | Yes (with trio) |
| **BP3** | Supporting | In-frame indel in repeat region | Consequence + RepeatMasker | Yes (with .osi) |
| **BP4** | Supporting-Very Strong | Computational evidence (benign) - REVEL **missense-only** with Very Strong band at ≤ 0.003 (Pejaver 2022) + SpliceAI ≤ 0.1 → Supporting (Walker 2023). BP4-splice is withheld from PVS1-territory consequences, from deep-exonic missense, and from in-frame indels / stop-lost / protein-altering variants, for which no calibrated benign predictor exists. | REVEL + SpliceAI | Yes |
| **BP5** | Supporting | Alternate molecular basis in case | Case-level analysis | No |
| **BP6** | Supporting | Reputable source reports benign | ClinVar (disabled by default per SVI) | Partial |
| **BP7** | Supporting | Synonymous + no splice impact + not conserved. Per Walker 2023: must NOT fire for synonymous at first base / last 3 bases of an exon (`at_exon_edge`); extends to deep-intronic offsets ≥ 7 (donor) or ≤ -21 (acceptor). | Consequence + SpliceAI + PhyloP + exon position | Yes |

### PP3/BP4 Strength Elevation (ClinGen SVI Calibrated)

PP3 and BP4 can be elevated beyond Supporting based on REVEL score (Pejaver 2022). BP4 reaches Very Strong only via REVEL — none of the other 12 calibrated tools reach that band.

| REVEL Score | PP3 Strength | BP4 Strength |
|---|---|---|
| ≥ 0.932 | Strong | — |
| ≥ 0.773 | Moderate | — |
| ≥ 0.644 | Supporting | — |
| ≤ 0.290 | — | Supporting |
| ≤ 0.183 | — | Moderate |
| ≤ 0.016 | — | Strong |
| ≤ 0.003 | — | **Very Strong** |

### PVS1 Strength Grading (Abou Tayoun 2018)

PVS1 is graded by a decision tree over null-variant context. The output code carries the strength suffix (e.g. `PVS1_Moderate`).

| Variant context | PVS1 strength |
|---|---|
| Nonsense / frameshift, NMD predicted | **Very Strong** (`PVS1`) |
| Canonical ±1/2 splice, NMD predicted | **Very Strong** (`PVS1`) |
| Whole-gene deletion in haploinsufficient gene | **Very Strong** (`PVS1`) |
| NMD-escape **in critical functional region** | Strong (`PVS1_Strong`) |
| NMD-escape, non-critical, **≥ 10 % protein removed** | Moderate (`PVS1_Moderate`) |
| Canonical splice in last exon | Moderate (`PVS1_Moderate`) |
| Start-loss with downstream Met ≤ 100 codons + pathogenic upstream | Moderate (`PVS1_Moderate`) |
| NMD-escape, non-critical, < 10 % protein removed | Supporting (`PVS1_Supporting`) |
| Start-loss with no corroborating evidence | Supporting (`PVS1_Supporting`) |

When the pipeline does not populate the tree signals (`predicted_nmd`, `protein_truncation_pct`, `is_last_exon`, `in_critical_region`, `alt_start_codon_distance`), PVS1 falls back to legacy Very Strong for any null variant in a LOF-intolerant gene. The graded result is exposed in `details` for transparency.

### PM3 Strength Grading (SVI v1.0 Points)

PM3 sums points across compound-het / homozygous observations and maps the total to a strength tier:

| Scenario | Points |
|---|---|
| Confirmed in-trans + Pathogenic companion | 1.0 |
| Confirmed in-trans + Likely Pathogenic companion | 0.5 |
| Phase unknown + Pathogenic | 0.5 |
| Phase unknown + Likely Pathogenic | 0.25 |
| Homozygous occurrence (proband hom-alt) | 0.5 each, capped at 1.0 total |

| Total points | PM3 strength |
|---|---|
| ≥ 0.5 | `PM3_Supporting` |
| ≥ 1.0 | `PM3` (Moderate) |
| ≥ 2.0 | `PM3_Strong` |
| ≥ 4.0 | `PM3_Very_Strong` |

In-cis companions are excluded (they count toward BP2 instead). Requires AR inheritance from OMIM.

### Anti-Double-Counting (Reconciliation Pass)

After per-criterion evaluation, a reconciliation pass suppresses computational evidence (PP3) that would double-count molecular signals already captured by other criteria. Suppressed criteria appear with `met: false` and `details.suppressed_by_reconcile` explaining why.

| Trigger | Action | Reference |
|---|---|---|
| PVS1 fires + PM1 fires | Suppress PM1 | Abou Tayoun 2018 (PVS1 already counts loss of the region; PM1 is residue-level evidence for missense) |
| PVS1 fires + PP3 driven by SpliceAI | Suppress PP3 | Walker 2023 (PVS1 already counts splice signal) |
| PS1 fires + PP3 driven by REVEL (missense) | Suppress PP3 | Pejaver 2022 (PS1 covers residue) |
| PM5 fires + PP3 driven by REVEL (missense) | Suppress PP3 | Pejaver 2022 (PM5 covers residue) |
| PP3 elevated to Strong + PM1 fires | Suppress PM1 (cap combined at Strong = 4 Tavtigian points) | Pejaver 2022 |

## Classification Combination Rules

### Pathogenic (8 rules)
1. PVS >= 1 AND PS >= 1
2. PVS >= 1 AND PM >= 2
3. PVS >= 1 AND PM >= 1 AND PP >= 1
4. PVS >= 1 AND PP >= 2
5. PS >= 2
6. PS >= 1 AND PM >= 3
7. PS >= 1 AND PM >= 2 AND PP >= 2
8. PS >= 1 AND PM >= 1 AND PP >= 4

### Likely Pathogenic (7 rules)
1. PVS >= 1 AND PM >= 1
2. **PVS >= 1 AND PP >= 1** (ClinGen SVI Sept 2020 — compensates for PM2 downgrade; Bayesian Post_P = 0.988)
3. PS >= 1 AND PM = 1-2
4. PS >= 1 AND PP >= 2
5. PM >= 3
6. PM >= 2 AND PP >= 2
7. PM >= 1 AND PP >= 4

### Benign (2 rules)
1. BA >= 1 (standalone)
2. BS >= 2

### Likely Benign (2 rules)
1. BS >= 1 AND BP >= 1
2. BP >= 2

### Conflicting Evidence

The combiner evaluates the pathogenic and benign rule sets **independently**. Conflicting → VUS only when **both directions reach a definite call** (P/LP and B/LB simultaneously); otherwise the directional call wins. So:

- PVS1 + PM2_Supporting + BP4_Supporting → LP (PVS+PP rule fires; the lone BP4_Supporting doesn't reach LB)
- PVS1 + 2 BS → VUS (Conflicting — both sides definite)
- PVS1 + BS1 alone → plain VUS (no benign rule fires; no "Conflicting" label)

This replaces the pre-PR9 behavior, which short-circuited any pathogenic-met + benign-met combination to VUS — over-zealous, since sub-threshold benign evidence shouldn't override a definite pathogenic call.

## Trio Analysis

When a multi-sample VCF is provided with `--proband`, `--mother`, and `--father` flags:

### De Novo Detection (PS2 / PM6)
- **PS2** (Strong): Proband carries variant, both parents homozygous reference, all three pass quality thresholds (DP >= 10, GQ >= 20 by default)
- **PM6** (Moderate): Partial trio — only one parent available or one parent fails quality, but available parent(s) are homozygous reference

### Compound Heterozygote Detection (PM3 / BP2)
After individual variant classification, a second pass groups variants by gene:
- **PM3** (Supporting → Very Strong): In a recessive-inheritance gene, the proband is het or hom for the variant. Companions in trans / phase-unknown contribute points (1.0 / 0.5 / 0.25 depending on phasing × ClinVar significance), homozygous occurrences add 0.5 (capped at 1.0). Total points ≥ 0.5 / 1.0 / 2.0 / 4.0 map to `PM3_Supporting` / `PM3` / `PM3_Strong` / `PM3_Very_Strong`. Phase-aware: uses phased genotypes (VCF `|` separator) when available.
- **BP2** (Supporting): Variant is in cis with a ClinVar pathogenic variant, or in trans with pathogenic in a dominant gene. Requires phased data.

## Output Format

### JSON (Web API and CLI `--output-format json`)

Each transcript consequence includes an `acmg` field:

```json
{
  "transcript_consequences": [{
    "gene_symbol": "BRCA1",
    "consequence_terms": ["frameshift_variant"],
    "impact": "HIGH",
    "acmg": {
      "classification": "Likely_pathogenic",
      "shorthand": "LP",
      "triggered_rule": "PVS + PM",
      "criteria": [
        {
          "code": "PVS1",
          "direction": "Pathogenic",
          "strength": "VeryStrong",
          "met": true,
          "evaluated": true,
          "summary": "Null variant in LOF-intolerant gene BRCA1 (pLI=1.00, LOEUF=0.03)"
        },
        {
          "code": "PM2_Supporting",
          "direction": "Pathogenic",
          "strength": "Supporting",
          "met": true,
          "evaluated": true,
          "summary": "Absent from gnomAD"
        }
      ],
      "counts": {
        "pathogenic_very_strong": 1,
        "pathogenic_strong": 0,
        "pathogenic_moderate": 0,
        "pathogenic_supporting": 1,
        "benign_standalone": 0,
        "benign_strong": 0,
        "benign_supporting": 0
      }
    }
  }]
}
```

### VCF CSQ Field

Two fields appended to the CSQ INFO annotation:
- `ACMG`: Classification shorthand (P, LP, VUS, LB, B)
- `ACMG_CRITERIA`: Met criteria codes joined by `&` (e.g., `PVS1&PM2_Supporting`)

### TSV Output

Two columns added after IMPACT:
- `ACMG`: Classification shorthand
- `ACMG_CRITERIA`: Met criteria codes (comma-separated)

## Architecture

### Crate: `fastvep-classification`

```
crates/fastvep-classification/src/
  lib.rs              # Public API: classify(), extract_classification_input()
  types.rs            # EvidenceStrength, EvidenceCriterion, AcmgClassification, AcmgResult
  sa_extract.rs       # ClassificationInput, typed SA deserialization, GenotypeInfo, CompanionVariant
  config.rs           # AcmgConfig, TrioConfig, GeneOverride, TOML loading
  combiner.rs         # 19 classification combination rules
  criteria/
    mod.rs            # evaluate_all_criteria() orchestrator
    pvs1.rs           # PVS1: null variant in LOF gene
    pathogenic_strong.rs    # PS1, PS2, PS3, PS4
    pathogenic_moderate.rs  # PM1, PM2, PM3, PM4, PM5, PM6
    pathogenic_supporting.rs # PP1, PP2, PP3, PP4, PP5
    benign_standalone.rs    # BA1
    benign_strong.rs        # BS1, BS2, BS3, BS4
    benign_supporting.rs    # BP1, BP2, BP3, BP4, BP5, BP6, BP7
```

### Data Flow

```
VCF Input
  |
  v
Consequence Prediction (fastvep-consequence)
  |
  v
Supplementary Annotation (fastvep-sa)
  |  Per-allele: ClinVar, gnomAD, REVEL, SpliceAI, dbNSFP
  |  Positional: PhyloP, GERP
  v
Gene-Level Annotation (fastvep-sa .oga)
  |  OMIM, gnomAD gene constraints, ClinVar protein index
  v
Sample Genotype Extraction (if trio configured)
  |  Parse FORMAT/GT/DP/GQ from VCF sample columns
  v
ACMG Classification Pass (fastvep-classification)
  |  1. extract_classification_input() -> ClassificationInput
  |  2. evaluate_all_criteria() -> Vec<EvidenceCriterion>
  |  3. combine() -> (AcmgClassification, triggered_rule)
  |  4. Store AcmgResult as serde_json::Value on AlleleAnnotation
  v
Compound-Het Enrichment Pass (if trio configured)
  |  Group variants by gene, detect companion relationships
  |  Re-evaluate PM3/BP2 with companion data
  v
Output (JSON / VCF CSQ / TSV)
```

### Integration Points

| File | Role |
|---|---|
| `crates/fastvep-annotate/src/lib.rs` | Web engine: loads .oga, runs gene annotation pass, ACMG classification, compound-het enrichment |
| `crates/fastvep-cli/src/pipeline.rs` | CLI batch: same pipeline with parallel variant processing |
| `crates/fastvep-io/src/variant.rs` | `AlleleAnnotation.acmg_classification: Option<serde_json::Value>` |
| `crates/fastvep-io/src/output.rs` | ACMG in JSON, VCF CSQ (`ACMG`, `ACMG_CRITERIA`), TSV |
| `crates/fastvep-web/src/handlers.rs` | Web API `acmg` request field |
| `web/index.html` | ACMG column, badges, evidence detail modal, summary chart |

## Limitations

1. **PS3/BS3** (functional studies): Cannot be automated — requires curated functional assay databases
2. **PP1/BS4** (segregation): Requires multi-generation pedigree with affection status beyond trio
3. **PP4** (phenotype specificity): Requires patient HPO phenotype terms
4. **BP5** (alternate molecular basis): Requires case-level multi-gene analysis
5. **PS4** is `NotEvaluated` by default — true PS4 needs case-control statistics; the legacy ClinVar-stars proxy is invalid per SVI. Opt back in via `use_clinvar_stars_as_ps4_proxy = true` for backward-comparable benchmarks.
6. **PS1/PM5/PM1** require the ClinVar protein index `.oga` file to be built and loaded. PS1 splice-RNA path requires the pipeline to populate `same_splice_position_pathogenic`; without it, splice PS1 is `evaluated: false`.
7. **PVS1 grading** uses Abou Tayoun 2018 signals (`predicted_nmd`, `protein_truncation_pct`, `is_last_exon`, `in_critical_region`, `alt_start_codon_distance`). When the pipeline cannot derive these for a transcript, PVS1 falls back to legacy Very Strong on any null variant in an LOF-intolerant gene.
8. **BP7 exon-edge / deep-intronic extension** uses optional `at_exon_edge` / `intronic_offset` fields. When unset, BP7 falls back to the legacy synonymous-only rule.
9. **BP3** requires RepeatMasker interval `.osi` file to be built and loaded
10. **PS2/PM6/PM3/BP2** require a multi-sample VCF with trio sample names configured
11. Compound heterozygote detection (PM3/BP2) works per-batch in the CLI; variants in different batches within the same gene may not be cross-referenced
12. **Multi-disorder genes** (SVI July 2025): the per-disorder override schema (`gene_overrides[GENE].disorders[DISORDER]`) is in place but the active-disorder selection mechanism is informational scaffolding pending a follow-up PR.

## References

- Richards S, et al. Standards and guidelines for the interpretation of sequence variants. *Genet Med*. 2015;17(5):405-424.
- Abou Tayoun AN, et al. Recommendations for interpreting the loss of function PVS1 ACMG/AMP variant criterion. *Hum Mutat*. 2018;39(11):1517-1524.
- Ghosh R, et al. Updated recommendation for the benign stand-alone ACMG/AMP criterion. *Hum Mutat*. 2018;39(11):1525-1530.
- Tavtigian SV, et al. Modeling the ACMG/AMP variant classification guidelines as a Bayesian classification framework. *Genet Med*. 2018;20(9):1054-1060.
- ClinGen SVI Recommendation for Absence/Rarity (PM2) — Version 1.0. Approved September 4, 2020.
- ClinGen SVI Recommendation for In-Trans Criterion (PM3) — Version 1.0.
- Pejaver V, et al. Calibration of computational tools for missense variant pathogenicity classification and ClinGen recommendations for PP3/BP4 criteria. *Am J Hum Genet*. 2022;109(12):2163-2177.
- Walker LC, et al. (ClinGen SVI Splicing Subgroup). Using the ACMG/AMP framework to capture evidence related to predicted and observed impact on splicing. *Am J Hum Genet*. 2023;110(7):1046-1067.
- ClinGen SVI Working Group. Guidance to VCEPs Regarding gnomAD v4 (March 2024).
- ClinGen SVI Working Group. Guidance Classifying Variants in Genes Associated with Multiple Disorders (July 2025).
- ClinGen Sequence Variant Interpretation (SVI) recommendations: https://clinicalgenome.org/tools/clingen-variant-classification-guidance/
