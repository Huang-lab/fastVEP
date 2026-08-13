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
| `--pick-order <LIST>` | Criteria order used by `--pick`, VEP `--pick_order` syntax | VEP's default |
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

# ── BP7's intronic extension (Walker 2023) ──
# Near boundary is Walker's: donor +7, acceptor -21. The far boundary is not in the
# paper and is set by measurement - see "Choosing BP7's far boundary" below.
# Both ends of the range are meaningful: a value below 7 turns the intronic extension
# off entirely, and one larger than any intron (1_000_000) restores it unbounded.
bp7_max_intron_offset = 300    # nucleotides from the nearest exon edge

# ── Conservation thresholds ──
phylop_conserved = 2.0         # PhyloP score above which position is "conserved"

# ── Gene constraint thresholds ──
pli_lof_intolerant = 0.9       # pLI score for LOF-intolerant gene
loeuf_lof_intolerant = 0.35    # LOEUF upper bound for LOF-intolerant gene
pp2_misz_threshold = 3.09      # Missense Z-score threshold for PP2

# ── PM1 hotspot detection ──
# Richards 2015 asks for a hot spot or critical domain "without benign variation", so
# both halves are tested. The benign half needs a clinvar_protein.oga built after benign
# assertions were indexed; against an older file it is skipped rather than assumed.
pm1_hotspot_window = 5         # Window size (amino acid positions) for hotspot scan
pm1_hotspot_min_pathogenic = 3 # Minimum pathogenic variants in window to call hotspot
pm1_max_benign_in_window = 0   # Benign/likely-benign ClinVar missense tolerated in the
                               # same window before it stops counting as a hotspot

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

# ── Combining rules ──
use_point_system = true                # Score the PATHOGENIC side with the ClinGen SVI point
                                       # system (Tavtigian 2020) instead of the Richards 2015
                                       # table. The benign side keeps the 2015 rules either
                                       # way - see "Which combining rules apply" below for the
                                       # measurement behind that split.

# ── Gene-level preconditions (require omim.oga / ClinGen GDV) ──
require_gene_disease_validity = true   # PVS1, PP2 and PM1 report NotEvaluated with
                                       # `no_valid_gene_disease_relationship` when EVERY ClinGen
                                       # curation of the gene is Limited/Disputed/Refuted/No Known
                                       # Disease Relationship. A gene ClinGen has not curated is
                                       # unaffected: absence is a fact about curation coverage,
                                       # not about the gene.
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

# ── Curated frequency-exception list ──
# Ships with Ghosh 2018's nine BA1 exceptions plus three hypomorphic alleles.
# Specifying ba1_exceptions in TOML REPLACES the whole shipped list; include
# the defaults alongside your additions if you want to keep them.
#
# `blocks` names which frequency criteria the entry suppresses. It defaults to
# ["BA1"], which is the scope Ghosh 2018 wrote for. Widen it for an allele whose
# frequency is what the disease model predicts rather than evidence against it:
# a hypomorph that causes disease only in trans with a rarer null allele is
# expected to be common AND expected to be seen homozygous in unaffected people,
# so BS1 and BS2 have to be blocked too.
#
# Matching is by coordinate first and HGVS second. A `c.` token identifies a
# variant only relative to the transcript it was written for.
[[ba1_exceptions]]
gene = "HFE"
hgvs_c = "c.845G>A"
reason = "p.Cys282Tyr - hereditary hemochromatosis"

[[ba1_exceptions]]
gene = "GAA"
hgvs_c = "c.-32-13T>G"
chrom = "17"
pos = 80104542
ref = "T"
alt = "G"
blocks = ["BA1", "BS1", "BS2"]
reason = "leaky splice allele of late-onset Pompe disease, pathogenic in trans with a null allele"

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
([`sweep_acmg_thresholds.py`](../analysis/acmg_benchmark/scripts/sweep_acmg_thresholds.py)) gives:

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

## Choosing BP7's far boundary

Walker 2023 extends BP7 from synonymous variants to intronic ones and gives that extension
a **near** boundary: donor-side `+7`, acceptor-side `-21`.
It gives no far boundary, and that is coherent for the SVI, whose recommendations are
written for a curator applying BP7 to one variant in a gene they know.
An automated pipeline applies it to every intronic position in the genome, and the two
situations are not the same.

Out in the deep intron the only evidence BP7 has left is a SpliceAI score, in precisely the
regime where SpliceAI is weakest.
A deep-intronic variant is pathogenic almost exclusively by activating a pseudoexon, and
pseudoexon activation is SpliceAI's documented blind spot.
This is not a hypothesis about our data: on all 36 deep-intronic ClinVar-pathogenic
variants that v13 called benign, SpliceAI is present and reads between 0.00 and 0.14.
The criterion is not missing its evidence out there - it is being actively misled by it.

So the far boundary is set by measurement
([`sweep_acmg_thresholds.py`](../analysis/acmg_benchmark/scripts/sweep_acmg_thresholds.py), full
673,660-variant benchmark, this key swept alone):

| `bp7_max_intron_offset` | correct benign calls | missed diagnoses | opposite-direction | marginal cost |
|---|---:|---:|---:|---|
| 50 | 176,132 | 25 | 69 | 248 correct benign per diagnosis recovered |
| 100 | 177,122 | 29 | 73 | 204 |
| 200 | 177,735 | 32 | 76 | 155 |
| **300** | **177,890** | **33** | **77** | **40** |
| 500 | 177,930 | 34 | 78 | 26 |
| 1000 | 177,981 | 36 | 80 | 31 |
| unbounded | 178,420 | 50 | 94 | - |

Pathogenic recall is unchanged at every setting; this knob moves only the benign side.

**Is the extension worth having at all?** That is the question the table above cannot answer,
because none of its rows turn the extension off. Scored separately, on the same concordance
definitions the benchmark's other tables use:

| | extension off (`0`) | **default (`300`)** |
|---|---:|---:|
| Exact match | 61.8 % | **62.7 %** |
| Same-direction | 75.9 % | **77.5 %** |
| Benign recall | 68.4 % | **71.5 %** |
| Likely-benign recall | 46.8 % | **51.6 %** |
| Pathogenic recall | 64.6 % | 64.6 % |
| **False-benign (missed diagnoses)** | **14** | 33 |
| False-pathogenic | 43 | 43 |
| **Opposite-direction (total)** | **57** | 76 |

Neither column dominates, and which one is "better" depends on the metric.
Off minimises hard errors; on maximises agreement.
Going from off to 300 buys 3.1 pp of benign recall - about 10,800 correct benign calls -
for 19 missed diagnoses, roughly **570 correct benign calls per diagnosis**.
That sits between the two frequency trades this classifier already makes (BS2's prevalence
bar at ~296, BS1's threshold at ~1,200), so declining the extension would mean holding it to
a stricter standard than the criteria beside it.

The deciding argument is not the exchange rate, though, it is provenance: Walker 2023 is a
published ClinGen SVI recommendation, and switching it off entirely is a deviation from
guidance while bounding it at a measured knee is not.
So 300 ships, and a lab that wants the conservative setting has it in one line:
`bp7_max_intron_offset = 0`.

Read the last column, not the first two.
Bounding the extension at 300 recovers 17 of its 50 missed diagnoses for 530 correct benign
calls, about **31 correct benign calls per diagnosis**.
Tightening past 300 costs 155, then 204, then 248 per diagnosis - a five- to eight-fold
step, and the knee is unambiguous.

For scale against a trade this classifier already makes: the BS2 prevalence bar above gives
up roughly **296** correct benign calls per false-benign call avoided, and the BS1 bar about
**1,200**.
By the project's own established exchange rate, every diagnosis recovered down to 300 is
cheap, and every one below it is not.

Both ends of the range stay meaningful, so there is no sentinel value:
`bp7_max_intron_offset = 0` turns the intronic extension off entirely (nothing can be below
Walker's near boundary of 7), and `1_000_000` restores it unbounded.
The bound applies only to the intronic extension; a synonymous variant is exonic and is
unaffected.

## PS1 for splicing

Richards 2015 defines PS1 on the protein: the same amino acid change as a previously established pathogenic variant.
Walker 2023 (ClinGen SVI Splicing Subgroup, PMID 37352859) extends it to splicing, where "the same change" means the same predicted RNA event rather than the same residue.
fastVEP implements the part of that extension an index can settle on its own.

### The rule as implemented

A canonical ±1/2 splice variant gets PS1 when **another** ClinVar `Pathogenic` variant sits on the **same canonical dinucleotide**.
Both variants abolish the same donor or acceptor, so Table 3's prerequisite - that the two predicted RNA events "precisely match" - holds by construction rather than by assumption.

The strength follows Table 3, and is graded on the PVS1 code the variant already carries:

| Baseline PVS1 on the variant | PS1 code |
|---|---|
| `PVS1` (Very Strong) | `PS1_Supporting` |
| `PVS1_Strong` / `_Moderate` / `_Supporting` | `PS1` (Strong) |
| No PVS1 at all | `PS1_Supporting` (no row of Table 3 covers this; the conservative one stands) |

The reduction under a full-strength PVS1 exists "to prevent overweighting of the VUA compared to the original (Likely) Pathogenic comparison variant".
It is applied in the reconciliation pass, because criteria are evaluated independently and PVS1's graded outcome is not visible to PS1 itself.

Combining PS1 with PVS1 is sanctioned rather than double counting.
The subgroup's [response to feedback](https://clinicalgenome.org/docs/clingen-svi-splicing-subgroup-response-to-feedback/) (22 March 2024, item 7b) is explicit: "if there is a relevant pathogenic variant with the same predicted impact as the variant under assessment, then you can use PS1 as well as PVS1(RNA)".

### Three restrictions, and why each is there

**The comparison variant must be `Pathogenic`, not `Likely_pathogenic`.**
Table 3 reads `N/A` in the LP column for a canonical-dinucleotide variant under assessment.
Item 7c of the same response document gives the reason: "since it is so easy for a ±1,2 dinucleotide variant to reach likely pathogenic, we placed constraints on using these variants as reference to make sure there actually was clinical evidence informing that pathogenic classification".
Whether an LP call rests on real clinical evidence is a curator's judgement, not something an index can answer, so LP entries are indexed but never fire.

**The comparison variant must be on the same dinucleotide, not merely in the same splice motif.**
Table 3 does cover a comparison variant elsewhere in the donor or acceptor motif, and does cover variants under assessment outside the canonical dinucleotide.
Neither is implemented, because for those the prerequisite has to be *established* - a variant at -3 may produce a different aberration from one at -1 - and establishing it needs a per-variant SpliceAI comparison the index does not carry.
Those rows are left to curators.

**The variant's own ClinVar record never counts.**
A comparison variant is by definition another variant, and the index is built from ClinVar, so a variant that is itself ClinVar Pathogenic will find its own entry.
The splice index carries REF and ALT, so the self entry is matched exactly rather than guessed at; this is why the splice path does not consult `exclude_self_from_clinvar_evidence`, which exists for the protein index precisely because that one *cannot* tell which entry is the variant being classified.

### The data it reads

`clinvar_protein.oga` carries the index, alongside the protein one it has always had:

```json
{"benignIndexed": true,
 "proteinVariants": [...],
 "spliceIndexed": true,
 "spliceAssembly": "GRCh38",
 "splicePositions": [{"pos": 7676000, "ref": "A", "alt": "G", "off": -2, "sig": "Pathogenic"}]}
```

It is built from `variant_summary.txt.gz`, whose `Name` column carries the HGVS `c.` token that identifies a `±1`/`±2` position, and whose `PositionVCF` / `ReferenceAlleleVCF` / `AlternateAlleleVCF` columns give the coordinate to match a call against.
Rebuild with:

```bash
fastvep sa-build --source clinvar_protein \
  -i variant_summary.txt.gz -o <sa-dir>/clinvar_protein
```

Two consequences of keying on genomic position rather than protein position are worth knowing.
The index is **assembly-specific**, so the builder filters to one assembly (`--assembly`, default GRCh38) and stamps the name into every gene record.
And `spliceIndexed` is what separates "this gene has no catalogued comparison variant" from "this file was built before the splice pass existed" - on an older file PS1 reports `evaluated: false` rather than a negative call.

### What it moved

On the 673,660-variant ClinVar 2-star+ benchmark, PS1's splice path fires 4,776 times, 98.8 % of them on pathogenic-side truth, and lifts the exact-match rate from 62.67 % to 63.11 % without moving a single benign-side call or adding an opposite-direction error.
The whole effect is 3,793 variants going from Likely pathogenic to Pathogenic, of which 3,363 are true Pathogenic.
The mechanism is arithmetic: `PVS1` (8 points) + `PM2_Supporting` (1) is 9, one short of the 10 that Pathogenic requires, and a large fraction of the pathogenic set was sitting on exactly that signature.
Full numbers in [`analysis/acmg_benchmark/RUN_VERSIONS.md`](../analysis/acmg_benchmark/RUN_VERSIONS.md).

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

## Which transcript `--pick` reports

`--pick` reduces a variant to one consequence block, and which one it keeps is decided by
`--pick-order`.
The default is VEP's, exactly:

```text
mane_select,mane_plus_clinical,canonical,appris,tsl,biotype,ccds,rank
```

Note where `rank` - consequence severity - sits: **last**.
Transcript status outranks how badly the variant damages the transcript.
At a locus where genes overlap, that produces a result which is correct VEP behaviour and wrong
clinical reporting, and the round-2 review caught it: CYP21A2 variants were reported on **C4B** as
`downstream_gene_variant`, and STRC-region variants on **TIMM9** as `upstream_gene_variant`, because
the neighbouring gene's MANE transcript outranked the disrupted non-MANE one.
Stock VEP makes the same choice; this is not a fastVEP deviation.

For clinical reporting, put `rank` first:

```bash
fastvep annotate ... --pick \
  --pick-order rank,mane_select,mane_plus_clinical,canonical,appris,tsl,biotype,ccds
```

Measured over 300 ClinVar variants in the CYP21A2 / STRC / KIAA0586 regions, the gene `--pick`
reports:

| gene reported | VEP default order | `rank` first |
|---|---:|---:|
| KIAA0586 | 187 | 217 |
| CYP21A2 | 34 | 77 |
| C4B (neighbour) | 43 | 0 |
| TIMM9 (neighbour) | 32 | 2 |
| TNXB | 4 | 4 |

Severity first, then the same tie-breaks as before, so equal-consequence choices still resolve to
MANE - the ordering changes which question is asked first, not which transcripts are eligible.

Criteria you leave out are not consulted at all rather than appended, so `--pick-order rank` really
does mean "severity, then transcript ID". Unknown or repeated criteria are rejected at startup with
the offending name. `length` (VEP's final tie-break) is rejected rather than silently ignored:
fastVEP's annotation record does not carry transcript length.

**The default stays VEP-identical** so a default run of each tool picks the same transcript, which
is what the head-to-head concordance numbers rest on. This is a flag, not a behaviour change.

## Gene-disease validity: curated-negative, never absent

`require_gene_disease_validity` blocks PVS1, PP2 and PM1 on **positive evidence against** a
gene-disease relationship, never on absence from the file.

ClinGen Gene-Disease Validity classifies each proposed gene-disease pair as Definitive, Strong,
Moderate, Limited, Disputed, Refuted or No Known Disease Relationship. The `.oga` carries all seven,
and the gate fires only when *every* curation of a gene falls in the bottom four. A gene ClinGen has
simply not reached is unaffected.

That distinction is not academic. ClinGen has curated roughly 2,400 genes to
Definitive/Strong/Moderate, so absence is overwhelmingly "nobody got to it yet". The first version
of this gate blocked on absence, and run v10 measured the cost: **1,497 truth-pathogenic null
variants lost PVS1**, in genes including SPAST, ABCB11, FLG and LAMB3 - all of which cause disease
and none of which are in the file. Worse, the genes the gate was written to catch (RYK, GIGYF2) are
absent from ClinGen too, so blocking on absence was not selecting for gene-disease validity at all.
It was selecting for curation coverage.

The consequence, stated plainly: fastVEP cannot automatically decline PVS1 for a gene a curator
believes is not disease-associated unless ClinGen has said so. RYK and GIGYF2 come back as
false-pathogenic in v11 for exactly this reason. `gene_overrides.<GENE>.disabled_criteria` is the
manual route.

## Which combining rules apply

Two schemes exist for turning met criteria into a classification, and they are encodings of the
same Bayesian model: the Richards 2015 combining table, and the ClinGen SVI point system
(Tavtigian et al. 2020, *Genet Med* 22:1735).
fastVEP uses **the point system for the pathogenic direction and the 2015 table for the benign
direction**, which is not a compromise but a measured choice.

Points are 1 / 2 / 4 / 8 for Supporting / Moderate / Strong / Very Strong, summed over met
pathogenic criteria: `>= 10` Pathogenic, `6..=9` Likely Pathogenic.

**Why points on the pathogenic side.** The table has a documented gap at a lone PVS1. It scores 8
points, squarely inside Likely Pathogenic, but matches no row of Table 5, so the table returns
Uncertain Significance. Run v10 put **2,319 truth-pathogenic variants in VUS on `PVS1` and nothing
else** - nonsense and frameshift variants in haploinsufficient disease genes, which is close to the
least uncertain thing a classifier sees.

**Why the table on the benign side.** Tavtigian's Likely Benign band opens at `-1`, so a single BP4
is enough for a benign call, where Richards requires two benign supporting criteria or a strong plus
a supporting. Measured on the benchmark:

| scheme | pathogenic recall | benign recall | false-benign |
|---|---:|---:|---:|
| 2015 table, both sides | 56.8 % | 56.3 % | 1 |
| full point system | 63.6 % | 71.7 % | **36** |
| full points, benign floor tightened to −2 | 63.6 % | **45.6 %** | 5 |
| **points pathogenic + table benign** | **63.6 %** | **56.3 %** | **1** |

The full point system buys 15 pp of benign recall for 36 false-benign calls, and **22 of those 36
are a lone BP4**. Tightening the band instead removes them but costs more benign recall than the
table had, because a lone BP4 carries a great many correct benign calls too. Using each scheme where
it is the safer one takes the pathogenic gain and pays nothing: a missed diagnosis and a variant
left uncertain are not comparable errors.

Set `use_point_system = false` for the 2015 table on both sides.

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
| **ClinVar gene index** | `clinvar_protein` | PS1, PM1, PM5, BP1 | Two indices in one file. Classified missense by protein position, both directions: pathogenic entries drive PS1's missense path, PM5 and PM1's hotspot count, benign entries drive PM1's "without benign variation" test. Plus canonical splice-dinucleotide variants by genomic position, driving [PS1's splice path](#ps1-for-splicing) |

### Optional Sources

| Source | SA Key | Used By | Description |
|---|---|---|---|
| **RepeatMasker** | `repeatmasker` | BP3 | Repeat region intervals (`.osi`). Built from the UCSC track - see [Building the RepeatMasker track](#building-the-repeatmasker-track-bp3) |

### Building the RepeatMasker track (BP3)

BP3 is "in-frame deletions/insertions in a repetitive region without a known
function" (Richards 2015), and it reads a positional interval database. The
source is the UCSC RepeatMasker track, public and unregistered:

```bash
curl -O https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz
python3 analysis/acmg_benchmark/scripts/sa_sources/repeatmasker_to_bed.py rmsk.txt.gz > repeatmasker.bed
fastvep sa-build --source custom_bed --name repeatmasker \
    -i repeatmasker.bed -o sa_databases/repeatmasker
```

`--name repeatmasker` is not cosmetic. The classifier finds the track by looking
for a supplementary key containing `repeat`, so a database built under another
name loads successfully and is then silently ignored.

5.3 M intervals on the primary assembly, about 240 MB as `.osi`. All repeat
classes are emitted, which is the literal reading of "repetitive region"; the
class is carried in the BED name column, because an in-frame indel inside a
simple tandem repeat and one inside an exonized Alu are not equally good
arguments for benignity and the distinction is invisible if it is dropped.

**What it is worth.** Measured over the full ClinVar 2-star-or-better set, BP3
fires on 1,307 variants, split 688 on benign-truth against 56 on
pathogenic-truth - directionally right, roughly 12:1. Its effect on the headline
is nil: benign recall moves 71.47 % to 71.49 %, and every error count is
unchanged. That is the expected behaviour of a lone Supporting criterion, which
cannot reach a classification by itself. BP3 earns its place by being visible to
a curator reading the evidence, not by moving an aggregate.

**Three states, not two.** An interval source only produces an annotation when
the position falls inside an interval, so "not in a repeat" and "no database
loaded" arrive at the classifier identically. The pipeline therefore tells the
classifier whether a repeat source was loaded at all, and BP3 reports
`evaluated: true, met: false` for an in-frame indel that was checked and found
outside a repeat, against `evaluated: false` when there was nothing to check
against. Without that distinction BP3 answered every in-frame indel outside a
repeat with "load RepeatMasker .osi", which reads as a setup error rather than
as the answer.

### Gene-level (.oga) sources

`fastvep sa-build` supports three gene-level sources, each producing a `.oga` file that the runtime picks up automatically from `--sa-dir`:

```bash
fastvep sa-build --source omim -i genemap2.txt -o sa/omim --assembly GRCh38
fastvep sa-build --source gnomad_genes -i gnomad.v4.1.constraint_metrics.tsv -o sa/gnomad_genes --assembly GRCh38
fastvep sa-build --source clinvar_protein -i variant_summary.txt.gz -o sa/clinvar_protein --assembly GRCh38
```

When a `.oga` is missing, dependent criteria (PVS1, PS1, PM1, PM5, PM3, BP1, BP2, PP2, BS2) degrade gracefully to `evaluated: false` rather than misfiring.
`clinvar_protein.oga` files built before benign assertions were indexed stay readable, and PM1 skips its benign-variation half against them rather than reading a structurally empty count as "no benign variation here": the record carries `benignIndexed` to tell the two apart.
`spliceIndexed` does the same job for PS1's splice path.
Rebuilding the source is what switches either half on.
Prefer `variant_summary.txt.gz` over `clinvar.vcf.gz` as the input: the VCF exposes neither a protein change for most records nor an HGVS `c.` token at all, so it yields very few protein entries and no splice index whatsoever.
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
| **PM1** | Moderate | Mutational hotspot / critical domain **without benign variation** - both halves of Richards 2015's wording are tested, so a window that also holds benign ClinVar missense does not count. NotEvaluated where the gene has no established gene-disease relationship: a cluster of assertions in an uncurated gene is not evidence of a critical region | ClinVar protein density (both directions) + ClinGen GDV | Yes (with .oga) |
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
| **BA1** | Standalone | AF > 5%, AN ≥ 2,000 (gnomAD v4 SVI March 2024). Tested against the **filtering allele frequency** (95% CI lower bound, max across ancestry groups; Whiffin 2017) where the database provides it, else the population maximum. The curated exception list blocks BA1 regardless of AF - Ghosh 2018's nine plus three hypomorphic alleles; see [Curated frequency-exception list](#curated-frequency-exception-list) | gnomAD population AFs + FAF + AN + HGVSc | Yes |
| **BS1** | Strong | Cross-population AF > expected, read from the same statistic as BA1 (filtering AF where available, else population maximum). NotEvaluated for homology-confounded genes, for gnomAD-flagged `segdup`/`lcr` sites, for non-PASS gnomAD records, and for ClinVar low-penetrance / risk alleles. Stands down where the gene's own BA1 bar already covers the frequency. Blocked outright by a curated exception whose `blocks` names BS1 | gnomAD per-population AFs + FAF + QC flags + ClinVar terms | Yes |
| **BS2** | Strong | Observed in healthy adult - AD/X-linked-D requires AC ≥ 5 (`bs2_ad_min_ac`). AR / X-linked counts individuals with no functional copy (homozygotes **plus hemizygotes** on a non-PAR sex-chromosome site) and requires ≥ `bs2_ar_min_hom` of them **and** a 95% lower bound on their frequency above `bs2_hom_prevalence_threshold` (default 1e-3, set by measurement - see [Choosing the BS2 prevalence bar](#choosing-the-bs2-prevalence-bar)), so the test scales with cohort size. Richards 2015's "full penetrance expected at an early age" qualifier is what this implements: one or two such individuals in a 730 K cohort is what a late-onset or reduced-penetrance disorder looks like, not tolerance Blocked outright by a curated exception whose `blocks` names BS2, which is how a hypomorphic allele's unaffected homozygotes are kept from reading as tolerance | gnomAD hom + hemizygote counts + AN + inheritance | Yes |
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

One branch remains unreachable because its signal has no plumbing yet: `alt_start_codon_distance` (start-loss), which therefore always lands at `PVS1_Supporting`.

#### Which NMD prediction the tree runs on

Abou Tayoun 2018 predicts nonsense-mediated decay with the **50-nucleotide rule**: a premature termination codon triggers decay only when it sits more than 50 nt upstream of the final exon-exon junction.
fastVEP can compute that exactly, and also carries the coarser **last-exon proxy** it has always used.
The two disagree in one place - a PTC in the 3' end of the penultimate exon, where the proxy claims decay and the measurement finds escape.

`pvs1_nmd_50nt_rule` chooses between them, and it defaults to `false` (the proxy) on a measurement rather than a preference.

On a 4,000-variant sample of PVS1-carrying ClinVar 2-star records, switching the exact rule on moved 58 calls from Likely Pathogenic to Uncertain: **54 of them ClinVar calls Pathogenic**, 4 Uncertain, none Benign.
Extrapolated to the full PVS1 population that is roughly 825 correct pathogenic calls given up to correct a single false-pathogenic one.

Nothing in that trade is a defect in either the rule or the tree.
A PTC 20 nt before the last junction really does escape decay, and Abou Tayoun really does grade an NMD-escaping truncation in a critical region at `PVS1_Strong` - 4 points, and Uncertain on its own.
ClinVar's curators call those variants Pathogenic because they also had segregation, case and functional evidence that no annotator can compute.
So enabling the flag makes fastVEP more faithful to the guideline and less concordant with ClinVar at the same time, and the choice belongs to the user rather than to a default.

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
6. **PS1's splice path covers only the same canonical dinucleotide.** Walker 2023 Table 3 also allows a comparison variant elsewhere in the same splice motif, and allows PS1 on variants *outside* the canonical dinucleotide. Neither is implemented, because both need the table's prerequisite - that the two variants' predicted RNA events "precisely match" - to be established rather than assumed. See [PS1 for splicing](#ps1-for-splicing). Requires a `clinvar_protein.oga` built by a version that carries `spliceIndexed`; older files make PS1 report `evaluated: false` on canonical splice variants rather than assuming no comparison exists.
7. **PVS1 grading** uses Abou Tayoun 2018 signals (`predicted_nmd`, `protein_truncation_pct`, `is_last_exon`, `in_critical_region`, `alt_start_codon_distance`). When the pipeline cannot derive these for a transcript, PVS1 falls back to legacy Very Strong on any null variant in an LOF-intolerant gene.
8. **BP7 exon-edge / deep-intronic extension** uses optional `at_exon_edge` / `intronic_offset` fields. When unset, BP7 falls back to the legacy synonymous-only rule.
9. **BP3** requires a RepeatMasker interval `.osi`. Build it from the UCSC track - see [Building the RepeatMasker track](#building-the-repeatmasker-track-bp3). Without one BP3 reports `evaluated: false` rather than assuming a variant is not in a repeat.
10. **PS2/PM6/PM3/BP2** require a multi-sample VCF with trio sample names configured
11. Compound heterozygote detection (PM3/BP2) works per-batch in the CLI; variants in different batches within the same gene may not be cross-referenced
12. **Multi-disorder genes** (SVI July 2025): the per-disorder override schema (`gene_overrides[GENE].disorders[DISORDER]`) is in place and the generated VCEP table populates it, but the active-disorder selection mechanism is informational scaffolding pending a follow-up PR - the gene-level values are what the classifier reads.
13. **The frequency dead zone between PM2 and BS1.** See below - this is the largest known gap on the frequency side and it is open by decision, not oversight.
14. **Penetrance is unsourced.** BS2's precondition is "full penetrance expected at an early age" and no public resource publishes penetrance per gene-disorder pair. Prevalence, onset and inheritance are covered; see [Per-gene frequency bars](#per-gene-frequency-bars-clingen-vcep-specifications).
15. **Published VCEP bars are applied without their founder-variant exceptions.** The exclusions are prose in the specification text, and `ba1_exceptions` carries no panel-specific entries beyond Ghosh 2018's nine, which is why the generated threshold table is opt-in rather than a default.
16. **The ClinVar-concordance benchmark cannot validate the BA1 exception list.** Seven of Ghosh 2018's nine variants are `Conflicting classifications` in ClinVar and carry no consensus truth label, so they are not in the truth set; the criterion exists for exactly the population a consensus benchmark discards. It is covered by an end-to-end test over all nine instead. Any future exception entries need the same treatment.
17. **The hypomorphic-allele list is three entries long.** GAA `c.-32-13T>G`, CFTR `c.1210-11T>G` and SPTA1 `c.4339-99C>T` ship with `blocks = ["BA1", "BS1", "BS2"]`, because an allele that causes disease only in trans with a rarer null is expected to be common and expected to be tolerated homozygous. The mechanism is general and the list is not; it grows by curation only. See [The curated frequency-exception list](#curated-frequency-exception-list).
18. **A variant at a site where the population data is vetoed cannot reach Likely pathogenic alone.** Where gnomAD's FILTER or a low-complexity / segmental-duplication flag makes the frequency unusable, BA1, BS1, BS2 and PM2 are all withheld while PVS1 - read from the same alignments - is not, so a lone PVS1 call would rest on a one-sided ledger. Such calls are reported Uncertain with the veto named. The rule applies only when the suppressed frequency would itself have met BA1 or BS1; applying it to every flagged site cost 2,630 correct pathogenic calls to remove 18 wrong ones.

## The frequency dead zone between PM2 and BS1

Both frequency bars are set by measurement, and between them is a range where fastVEP
returns no frequency evidence in either direction.

For a dominant or uncharacterised-inheritance gene, PM2 requires a filtering allele
frequency below `pm2_ad_af_threshold` (4e-5) and BS1 fires above `bs1_af_threshold` (5e-3).
A variant sitting between those - roughly 1 in 25,000 to 1 in 200 - is too common to be
"extremely rare" and not common enough to be "greater than expected for disorder", so
neither criterion applies and the variant carries no frequency evidence at all. That band is
about two orders of magnitude wide and it is not a narrow edge case: it is where a great
many uncertain missense variants in dominant genes land.

The gap is real in both directions but it hurts most on the benign side. A variant at 1e-3
in a dominant, early-onset, fully penetrant condition is almost certainly too common to
cause it, and fastVEP says nothing. The reviewer's CHD7 (30 heterozygotes) and PTCH1 (15)
rows are exactly this shape, and neither is resolved at any *global* BS1 threshold that does
not also cost missed diagnoses elsewhere - the sweep in
[Choosing the BS2 prevalence bar](#choosing-the-bs2-prevalence-bar) shows the same behaviour
for its own bar.

**Why it is not closed by moving a number.** Narrowing the band means lowering BS1, and the
BS1 sweep recorded in [`bs1_af_threshold`](../crates/fastvep-classification/src/config.rs)
shows the cost climbing steeply below 0.5 %: past 0.05 % the missed-diagnosis count roughly
triples. A single global value cannot be right for both a dominant early-onset condition and
a common recessive one, which is what Richards 2015's actual wording - "greater than
expected **for disorder**" - has been asking for all along.

**What would close it** is a per-gene-disease attribute table: prevalence, penetrance class,
onset class and any published VCEP BA1/BS1 threshold, keyed by gene and disorder. Whiffin
2017 gives the formula. Most of that table now exists - see
[Per-gene frequency bars](#per-gene-frequency-bars-clingen-vcep-specifications) - and what
remains missing is penetrance, which nothing publishes in machine-readable form.

## Per-gene frequency bars (ClinGen VCEP specifications)

Where a ClinGen Variant Curation Expert Panel has published a BA1, BS1 or PM2 bar for a gene,
it outranks fastVEP's global default. The panel has looked up that disorder's prevalence,
penetrance and allelic heterogeneity; a number measured across all genes cannot.

`gene_overrides` carries them, one per gene:

```toml
[gene_overrides.ABCA4]
ba1_af_threshold = 0.163      # the ABCA4 VCEP's published bar, 3x looser than the default
bs1_af_threshold = 0.0163
pm2_af_threshold = 0.0001
```

A generated table of 304 bars across 117 genes ships at
[`analysis/acmg_benchmark/data/vcep_thresholds.toml`](../analysis/acmg_benchmark/data/vcep_thresholds.toml),
read out of the [ClinGen CSpec Registry](https://cspec.genome.network/cspec/ui/svi/) by
[`build_vcep_thresholds.py`](../analysis/acmg_benchmark/scripts/build_vcep_thresholds.py).
Load it with `--acmg-config`. Alongside it is an audit table carrying the verbatim sentence
every number was read out of, and every rejection with its reason.

**It is opt-in, and that is deliberate.** On the ClinVar 2-star+ benchmark it buys 2,028 more
correct Benign calls and lifts exact match by 0.26 pp, at a cost of four new false-benign
calls. Two of those four are founder variants - ANO5 and CAPN3, both limb-girdle muscular
dystrophy - that the panels explicitly excluded, in prose the generator cannot read. 32 of
the 286 usable bars carry a caveat of that shape. `ba1_exceptions` is the right home for
them and is currently empty of panel-specific entries, so a tightened bar today is applied
without the exclusions that were published with it. Filling that list is what has to come
before these bars become a default.

The wider attribute table -
[`gene_disease_attributes.tsv`](../analysis/acmg_benchmark/data/gene_disease_attributes.tsv),
6,857 gene-disorder pairs - carries inheritance (96 % coverage), earliest age of onset
(91 %) and prevalence class (85 %) from Orphanet. **Penetrance is 0 % and is emitted blank**,
because no public resource publishes it per gene-disorder pair in machine-readable form. It
is the field BS2's "full penetrance expected at an early age" precondition most wants, and
closing it needs a clinician rather than a download.

## Curated frequency-exception list

A per-gene bar answers "how common is too common in this gene".
It cannot answer "this particular allele's frequency is not evidence about it at all", and
two classes of allele need exactly that answer.

**Common and pathogenic.** Ghosh 2018 (PMID 30311378) lists nine variants that sit above 5 %
in at least one population and must not be called standalone-benign on that basis. All nine
ship as defaults with `blocks = ["BA1"]`, the scope the paper wrote for.

**Common because the disease model says so.** A hypomorphic allele that causes disease only
in trans with a rarer null allele is *expected* to be frequent, and is *expected* to be seen
homozygous in unaffected people. Reading either as evidence against pathogenicity inverts
the mechanism, so an entry of this kind sets `blocks = ["BA1", "BS1", "BS2"]`. Three ship as
defaults, all raised by the round-5 medical-genetics review:

| Variant | gnomAD AF | Homozygotes | Mechanism |
|---|---:|---:|---|
| GAA `c.-32-13T>G` | 5.4e-3 | 23 | leaky splicing; late-onset Pompe in trans with a null allele |
| CFTR `c.1210-11T>G` | 9.8e-3 | 27 | the 5T poly-pyrimidine tract; CBAVD or CF in trans with a CF-causing variant |
| SPTA1 `c.4339-99C>T` | 6.6e-3 | 40 | alpha-LELY; elliptocytosis in trans with an SPTA1 mutation |

None of the three reaches 5 %, which is why a BA1-only list never touched them, and none
reaches Likely pathogenic even with the exception applied - the deciding evidence is PM3,
the second allele in trans, which no annotator computes from a single variant.

Matching is by genomic coordinate first and HGVS second, and the order is load-bearing
rather than a preference: a `c.` token identifies a variant only relative to the transcript
it was written for, and two of Ghosh's nine cannot be matched on their published HGVS at all.

The list grows by curation. There is no automated route to it, and an entry without a stated
mechanism is not one worth adding.

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
