# ACMG-AMP Variant Classification in fastVEP

fastVEP implements automated ACMG-AMP variant classification based on the standards published by Richards et al. 2015 (*Genet Med* 17:405-424), with ClinGen Sequence Variant Interpretation (SVI) working group recommendations incorporated.

## Overview

The classifier implements all 28 ACMG-AMP evidence criteria and evaluates them for each variant-allele-transcript combination (a default run reports 26, because PP5 and BP6 are off unless `use_pp5_bp6` is set) and produces a 5-tier classification:

| Classification | Shorthand | Color (Web UI) |
|---|---|---|
| Pathogenic | P | Red |
| Likely Pathogenic | LP | Orange |
| Uncertain Significance | VUS | Gray |
| Likely Benign | LB | Blue |
| Benign | B | Green |

## Quick Start

> **First, build the databases.** Every `--sa-dir ./sa/` below stands for the nine-source reference stack.
> A source that is missing does not stop a run: each criterion that reads it reports `evaluated: false` instead of firing, so a half-built `--sa-dir` produces VUS rather than an error.
> The command line says which sources are missing, once, on stderr (see [Missing sources](#missing-sources)).
> [ACMG_SETUP.md](ACMG_SETUP.md) has the verified download URL, build command and measured size for each source, plus `scripts/check_acmg_stack.py` to confirm each one answers a query before you trust a run.

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

The checkbox is available only for a genome whose database directory holds all nine sources of the reference stack.
For any other genome it is disabled and the page names the sources that are missing, because a classification made without them would rest on less evidence than it appears to (see [Missing sources](#missing-sources)).
`POST /api/annotate` with `"acmg": true` for such a genome is refused with HTTP 422 and the same list, and `--acmg-allow-partial` answers it anyway and says so in the response ([docs/API.md](API.md#acmg-needs-the-reference-stack)).

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
bs1_af_threshold = 0.005       # BS1: benign strong (greater than expected for disorder)

# ── PM2 inheritance-aware thresholds (not from ClinGen SVI PM2 v1.0, which sets only the strength) ──
pm2_ad_af_threshold = 0.00004  # PM2: AD/unknown — "extremely rare", not strictly absent.
                               # Chosen by measurement and sits inside the range published
                               # VCEP specs use for dominant genes; see "Choosing the PM2
                               # threshold for dominant genes" below. Set 0.0 for the literal
                               # Richards 2015 strict-absence reading.
pm2_ar_af_threshold = 0.00007  # PM2: recessive-only - AF ≤ 0.00007 (0.007%); equals the Hearing Loss VCEP's threshold, swept for the manuscript: flat
pm2_af_threshold = 0.0001      # Legacy single-threshold field, retained for back-compat
                               # with pre-PR4 configs; not consulted in the default path.
pm2_absent_when_no_record = true  # PM2: a variant with no gnomAD record is read as absent (unless ClinVar's own
                               # population frequencies say it is not extremely rare). Coverage is not read, so
                               # this includes sites that gnomAD could not see; false leaves PM2 not evaluated
                               # wherever there is no record. A run with no gnomAD database loaded at all never
                               # reads a missing record as absence, whatever this says. See limitation 30.

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
# The near end of the range is Walker's: donor +7, acceptor -21. The other end, the intronic
# distance, is not in the paper and is set by measurement - see "Choosing BP7's intronic distance" below.
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

# ── Allele-number floor for BA1/BS1 (from the dataset-size condition of the BA1 recommendation,
#    Ghosh 2018; the March 2024 gnomAD v4 communication states no floor) ──
min_an_for_frequency_criteria = 2000  # BA1/BS1 require AN ≥ this; below → NotEvaluated
                                      # (BS2 only needs an AN to exist)

# ── BS2 (Richards 2015 "full penetrance ... at an early age") ──
# Both branches hold a 95% lower bound on the frequency of the individuals the disorder
# would affect against a maximum credible prevalence, so the criterion scales with cohort
# size rather than with a bare count. Both thresholds were chosen by measurement, not convention
# - see "Choosing the BS2 prevalence threshold" below.
bs2_ad_min_ac = 5                       # Dominant: absolute floor on carrier alleles
bs2_ad_prevalence_threshold = 3e-4      # Dominant: the lower bound on the carrier frequency
                                        # (AC - nhomalt individuals) must exceed this. 0
                                        # switches the threshold off and leaves the floor as the
                                        # whole test
bs2_ar_min_hom = 2                      # Recessive / X-linked: absolute floor on observed
                                        # individuals with no functional copy (homozygotes,
                                        # plus hemizygotes on a non-PAR sex-chromosome site)
bs2_hom_prevalence_threshold = 1e-3     # Recessive / X-linked: the lower bound on their
                                        # frequency must exceed this

# NOTE: `bs2_ar_min_hom` is a belt-and-braces floor, not an active knob. A sweep of the
# full benchmark found it inert at the default prevalence threshold: values from 1 to 10 give
# identical reference-discordant benign and opposite-direction counts, because the prevalence test
# subsumes it. It is kept so that a config lowering the prevalence threshold still cannot fire
# BS2 on a single observation. `bs2_ad_min_ac` is not inert: see "The dominant threshold".

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
# Ships with Ghosh 2018's nine BA1 exceptions, three hypomorphic alleles and three
# alleles the round-6 review named.
# Specifying ba1_exceptions in TOML REPLACES the whole default list; include
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

# ── Curated disease mechanisms (default table; PCSK9, IFIH1, MYH7, RASopathies, ...) ──
# Specifying gene_mechanisms in TOML REPLACES the default table.
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

# Cap or raise a criterion for one gene, which is what a VCEP specification
# routinely does. Keyed by the base code (`PP3` catches PP3_Strong and
# PP3_Moderate too) or by an exact graded code, exact winning. A restrengthed
# criterion is renamed to match, so `PP3_Strong` capped to Supporting is
# reported as `PP3`, never as a code and a strength that disagree.
[gene_overrides.MYH7.strength_overrides]
PP3 = "Moderate"
```

Every key above is checked on load. An unrecognised one is an error, not a
shrug: a file like this exists precisely to depart from the defaults, and a
misspelling that silently returns them tells the reader nothing. This
classifier has included three separate criteria that never fired, and each was
invisible for want of exactly that objection.

## Choosing the BS2 prevalence threshold

`bs2_hom_prevalence_threshold` is the one threshold in the classifier that was set by
measurement rather than by convention, so it is worth showing the working.

BS2 fires on a recessive or X-linked gene only when the 95% lower confidence bound on the
frequency of individuals with no functional copy exceeds this threshold.
The threshold is therefore a **maximum credible disease prevalence**: the point above which the
observed null individuals cannot all be explained by the disorder itself.

Sweeping it across the full 673,660-variant ClinVar 2-star+ benchmark
([`sweep_acmg_thresholds.py`](../analysis/acmg_benchmark/scripts/sweep_acmg_thresholds.py)) gives:

| threshold | BS2 fires | reference-discordant benign calls | correct benign calls | marginal cost |
|---|---:|---:|---:|---|
| 1e-6 | 63,285 | 54 | 139,270 | |
| 1e-5 | 51,875 | 45 | 136,014 | 362 correct benign lost per reference-discordant benign call avoided |
| 5e-5 | 42,219 | 41 | 133,949 | 516 |
| 1e-4 | 37,808 | 40 | 133,407 | 542 |
| **1e-3** | **27,290** | **38** | **132,815** | **296** |

Two things follow.

**The data does not pick a value by itself.** The curve is smooth, with no knee, so every
setting is a defensible point on a trade. What the sweep does establish is the exchange
rate, and that no setting of this or any other BS2 knob drives reference-discordant benign calls below 38;
BS2 appears in only 9 of the 45 at the old default, so the residual is not a threshold
problem at all.

**The choice rests on the parameter's meaning and on which error is worse.** As a maximum
credible prevalence, the threshold has to cover the most prevalent Mendelian conditions BS2 is
applied to, not the typical one. Hearing loss, alpha-1 antitrypsin deficiency and familial
Mediterranean fever in high-prevalence populations all sit near 1 in 1,000, so a 1e-5 threshold is
two orders of magnitude too tight for exactly the disorders a medical geneticist flagged in
review. If ClinVar is right, a benign call on a pathogenic variant is the costlier error; a lost benign call becomes a VUS and
costs triage effort. That asymmetry favours the safer threshold, and the step to 1e-3 is also the
cheapest on the curve.

Hence the 1e-3 default. It is a plain config value: a lab that weighs VUS burden more
heavily sets a smaller one, and a gene-specific prevalence from a VCEP specification should
override it per gene once that table exists.

### The dominant threshold

The dominant branch used to fire at `AC >= bs2_ad_min_ac` (5), whatever the size of the cohort.
Nine carriers among gnomAD v4's 1.4 million alleles is not "observed in a healthy adult" for a disorder with full penetrance at an early age, because a count that small is about what a rare dominant disease would produce by itself.
SNRNP200 p.Ser1087Leu, which ClinVar classifies Pathogenic, is the round-6 reviewer's case: nine carriers, and BS2 called it tolerated.

The branch now asks the question the recessive branch asks.
It counts carriers, `AC - nhomalt`, so that a homozygote counts once.
It requires the 95% lower bound on their frequency to exceed `bs2_ad_prevalence_threshold`, and the floor as well.
A dominant-only gene is decided in that branch alone.
Master had no threshold: a dominant-only variant that cleared the floor fired at once, and one below the floor that had a homozygote fell into the recessive branch and could be fired there.
Neither can reach the recessive branch now, so a variant that clears the floor but not the threshold cannot be fired by its homozygotes.

Each row below is a whole run of the final classifier over the 673,660-variant benchmark ([`run_benchmark.sh`](../analysis/acmg_benchmark/scripts/run_benchmark.sh), one `ACMG_CONFIG` per row, the configurations in [`results/v37/configs`](../analysis/acmg_benchmark/results/v37/configs)), so the rows from "threshold off" down differ in the BS2 setting alone.
The master row is v36 and differs from the threshold-off row by 534 calls that are the other changes of this round.
The columns count calls in the same direction as ClinVar's 2-star+ classification, which is a reference and not ground truth.

| Run | Exact | Same-direction | Opposite | P same | LP same | LB same | B same |
|---|---:|---:|---:|---:|---:|---:|---:|
| v36 (master) | 425,328 | 522,542 | 59 | 51,796 | 6,750 | 66,046 | 111,965 |
| threshold off (floor 5) | 425,289 | 522,388 | 49 | 51,793 | 6,747 | 65,825 | 111,855 |
| 3e-5 | 425,773 | 522,873 | 49 | 51,839 | 6,753 | 65,610 | 111,837 |
| 1e-4 | 425,877 | 522,938 | 49 | 51,852 | 6,755 | 65,500 | 111,803 |
| **3e-4 (default)** | 425,869 | 522,886 | 49 | 51,855 | 6,755 | 65,461 | 111,745 |
| 1e-3 | 425,753 | 522,819 | 49 | 51,855 | 6,755 | 65,437 | 111,693 |
| 3e-3 | 425,653 | 522,775 | 49 | 51,855 | 6,755 | 65,427 | 111,652 |

Against the threshold switched off, 3e-4 recovers 70 Pathogenic and Likely pathogenic calls (58,540 to 58,610 same-direction) and gives up 474 benign-direction ones (177,680 to 177,206), about 7 to 1.
Nothing above 3e-4 recovers another pathogenic call, and every step above it only gives up more benign ones.
Below it the recovery is incomplete: 58,592 at 3e-5 and 58,607 at 1e-4.
That is the rule the default follows, the smallest threshold on the swept grid at which every Pathogenic and Likely pathogenic call the threshold can recover is recovered.
It is a statement about this benchmark and about which error is worse, not a claim about any disease's prevalence.
The data alone do not pick it: exact and same-direction agreement both peak at 1e-4, three Pathogenic calls short of it.
The asymmetry the recessive threshold was chosen on applies here as well, since, if ClinVar is right, a benign call on a pathogenic variant costs more than a benign call lost to VUS.

Against the threshold switched off, 1,518 calls change.
1,007 of them move toward ClinVar's classification and 438 away, and the other 73 stay as far from it as they were.
Of the 70 recovered, 69 go from VUS to Likely pathogenic and one to Pathogenic.
68 of them carry a PVS1 that BS2 had been cancelling (BEST1, ALG8, NEK8, APOB and COL11A1 among them) and the other two are CASR missense variants.
The same change moves 32 variants ClinVar calls VUS to Likely pathogenic, which is the opposite error type and part of what the trade costs.
22 of them carry a PVS1 that BS2 had been cancelling, and the other 10 are missense variants on PP3 or PM1.
474 benign-direction calls are lost to VUS, 364 Likely benign and 110 Benign, and 7 Benign calls become Likely benign, which stays benign-direction.

A flat floor does the same work with the same count in a cohort of any size:

| Run | Exact | Same-direction | Opposite | P same | LP same | LB same | B same |
|---|---:|---:|---:|---:|---:|---:|---:|
| floor 5 (the old rule) | 425,289 | 522,388 | 49 | 51,793 | 6,747 | 65,825 | 111,855 |
| floor 20 | 425,721 | 522,827 | 49 | 51,833 | 6,753 | 65,656 | 111,844 |
| floor 100 | 425,872 | 522,918 | 49 | 51,853 | 6,755 | 65,488 | 111,786 |
| floor 1,000 | 425,726 | 522,798 | 49 | 51,855 | 6,755 | 65,434 | 111,670 |
| threshold 3e-4, floor 5 | 425,869 | 522,886 | 49 | 51,855 | 6,755 | 65,461 | 111,745 |

A floor of 1,000 reaches the same pathogenic recovery as the threshold, 58,610, and gives up 576 benign-direction calls against the threshold's 474.
A floor of 100 has a slightly better exchange rate and recovers two fewer.
What the threshold adds is that it scales: 1,000 carriers is a large share of a cohort of a few thousand and a small one of 730,000 individuals.

**This corrects an earlier note.**
`bs2_ad_min_ac` was described here as barely affecting the final call.
That was true of the reference-discordant benign and opposite-direction counts, which is all the sweep behind it looked at, and it is not true of agreement: with the threshold off, raising the floor from 5 to 20, 100 and 1,000 recovers 46, 68 and 70 Pathogenic and Likely pathogenic calls.

**Applying BA1's and BS1's minimum AN to BS2 was tried and rejected.**
We read the SVI guidance behind that minimum as being about the noise of a frequency estimate, and a count test over a small cohort is already as strict as the cohort is small.
With the threshold off it moved 693 calls, 690 of them away from ClinVar's classification.
BS2 is not evaluated for a record with no AN or an AN of zero, since the cohort size is then unknown, and takes no other minimum.

**The sweeps above are of the dominant arm as it was reached before v38.**
Until then a gene was dominant here only if its disease name said so, which was 50 of the 2,924 genes of the ClinGen file the benchmark read.
v38 reads ClinGen's curated mode of inheritance instead, and 746 genes are dominant-only (limitation 20), so the tables above describe the arm on a small and not typical set of genes, and the threshold and the floor were swept again on the genes it now reaches (`analysis/acmg_benchmark/RUN_VERSIONS.md`, "v38").
The carrier threshold is no longer flat there: each step from 3e-3 to the default 3e-4 adds benign-side calls and no reference-discordant benign call, and below it a step buys 715, 678 and 315 correct benign-side calls per added reference-discordant benign call (3, 7 and 23 more at 1e-4, 3e-5 and no threshold), so the default threshold is the last step that adds no error.
With the threshold off, a floor of 20 alleles gives 33 reference-discordant benign calls, 100 gives 26 and 1,000 gives 24.
It also stops counting carriers beside BS1 and BA1, which read the same frequency: see the last row of the [reconciliation table](#anti-double-counting-reconciliation-pass).

One other knob was measured and found inert on the reference-discordant benign count:

- `bs2_ar_min_hom` is fully subsumed by the prevalence test: 1 versus 10 gives identical results.
  It is retained only as a floor for configs that lower the prevalence threshold.

## Choosing the PM2 threshold for dominant genes

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
non-zero thresholds, clustered between 1e-5 and 1e-4 (the Cardiomyopathy VCEP's MYH7
specification uses 4e-5; the PTEN VCEP uses 1e-5).

Sweeping the key alone over the ClinVar 2-star+ benchmark:

| threshold | pathogenic recall | reference-discordant pathogenic calls | benign recall |
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

**4e-5** is the choice: it is at the knee, it costs one additional reference-discordant pathogenic call
on the sweep sample, and it sits inside the range VCEPs have actually published for
dominant genes rather than being a number of our own.
Going further would buy a fraction of a point while taking us past any published VCEP
specification.

The literal 2015 reading remains one line away for a lab that wants it
(`pm2_ad_af_threshold = 0.0`).

**Caveat, the same one that applies to BS2.** ClinVar significance is itself assigned
partly from gnomAD frequency, so a frequency threshold tuned against ClinVar is partly
fitted to its own input. The curve's shape is the useful output, not its absolute level.

**The knee above belongs to the 2015 combining table, which the sweep was made under.**
With the point system that the classifier now uses, a lone PVS1 is already Likely
pathogenic without PM2, and the same sweep is nearly flat: recall of ClinVar's
pathogenic-side variants moves from 61.97% at strict absence to 62.49% at 5e-4 on v37, and
from 61.98% to 62.48% on v38, which reads inheritance from ClinGen's curated field.
Under the 2015 table, one switch away, the knee is back (38.8% at strict absence and 56.7%
at 4e-5).
A threshold can move only the calls that need the criterion it gates, and the combining rule
decides how many do (`manuscript/acmg`, Section 3.3 and Section 3.8).

## Choosing BP7's intronic distance

Walker 2023 extends BP7 from synonymous variants to intronic ones and gives that extension
a **near** boundary: donor-side `+7`, acceptor-side `-21`.
It gives no intronic distance, and that is coherent for the SVI, whose recommendations are
written for a curator applying BP7 to one variant in a gene they know.
An automated pipeline applies it to every intronic position in the genome, and the two
situations are not the same.

Two other departures from the paper, found when it was read against the code in October 2026.
Walker 2023 applies BP7 after BP4 has been assigned, which for SpliceAI is a score of 0.1 or
below, and cautions against a conservation test for BP7 unless its tool, depth and cut point
have an empirical justification.
This classifier takes a score below 0.2 as "no impact" and keeps Richards 2015's "not highly
conserved" as a PhyloP below `phylop_conserved`.
Limitation 21 below gives what aligning it would change.

Out in the deep intron the only evidence BP7 has left is a SpliceAI score, in precisely the
regime where SpliceAI is weakest.
A deep-intronic variant is pathogenic almost exclusively by activating a pseudoexon, and
pseudoexon activation is SpliceAI's documented blind spot.
This is not a hypothesis about our data: on all 36 deep-intronic ClinVar-pathogenic
variants that v13 called benign, SpliceAI is present and reads between 0.00 and 0.14.
The criterion is not missing its evidence out there - it is being actively misled by it.

So the intronic distance is set by measurement
([`sweep_acmg_thresholds.py`](../analysis/acmg_benchmark/scripts/sweep_acmg_thresholds.py), full
673,660-variant benchmark, this key swept alone):

| `bp7_max_intron_offset` | correct benign calls | reference-discordant benign calls | opposite-direction | marginal cost |
|---|---:|---:|---:|---|
| 50 | 176,132 | 25 | 69 | 248 correct benign per reference-discordant benign call avoided |
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
| **Reference-discordant benign** | **14** | 33 |
| Reference-discordant pathogenic | 43 | 43 |
| **Opposite-direction (total)** | **57** | 76 |

Neither column dominates, and which one is "better" depends on the metric.
Off minimises hard errors; on maximises agreement.
Going from off to 300 buys 3.1 pp of benign recall - about 10,800 correct benign calls -
for 19 reference-discordant benign calls, roughly **570 correct benign calls per reference-discordant benign call avoided**.
That sits between the two frequency trades this classifier already makes (BS2's prevalence
threshold at ~296, BS1's threshold at ~1,200), so declining the extension would mean holding it to
a stricter standard than the criteria beside it.

The deciding argument is not the exchange rate, though, it is provenance: Walker 2023 is a
published ClinGen SVI recommendation, and switching it off entirely is a deviation from
guidance while bounding it at a measured knee is not.
So 300 ships, and a lab that wants the conservative setting has it in one line:
`bp7_max_intron_offset = 0`.

Read the last column, not the first two.
Bounding the extension at 300 avoids 17 of its 50 reference-discordant benign calls for 530 correct benign
calls, about **31 correct benign calls per reference-discordant benign call avoided**.
Tightening past 300 costs 155, then 204, then 248 per reference-discordant benign call avoided - a five- to eight-fold
step, and the knee is unambiguous.

For scale against a trade this classifier already makes: the BS2 prevalence threshold above gives
up roughly **296** correct benign calls per reference-discordant benign call avoided, and the BS1 threshold about
**1,200**.
By the project's own established exchange rate, every reference-discordant benign call avoided down to 300 is
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

### Four restrictions, and why each is there

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

**The variant must abolish the dinucleotide where HGVS places it.**
A record anchored on the dinucleotide can describe a change that HGVS places wholly inside the intron past ±2, or wholly inside an exon in frame, and that change abolishes no dinucleotide.
PS1's splice path stands down for it, as PVS1's splice track does, and records the placement; see [Where HGVS places an indel](#where-hgvs-places-an-indel).

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

On the 673,660-variant ClinVar 2-star+ benchmark, PS1's splice path fires 4,776 times, 98.8 % of them on variants ClinVar classifies as Pathogenic or Likely pathogenic, and lifts the exact-match rate from 62.67 % to 63.11 % without moving a single benign-side call or adding an opposite-direction error.
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

## Where HGVS places an indel

A consequence term is decided from the position a VCF record anchors a change at.
An indel in a repeat can be written at several positions that all edit the sequence the same way, and HGVS's 3' rule names the most 3' one; [VEP_DIVERGENCE.md](VEP_DIVERGENCE.md) has how the shift is done and checked.
The two can disagree about whether the change touches an exon at all, and PVS1, PS1's splice path, PM1, PM4, BP3 and BP4's splice path used to follow the record.

| Record | Term the record earns | HGVS | What the edited sequence is |
|---|---|---|---|
| MSH6 `2:47806358 G>GGTAT` | `frameshift_variant` | `c.3801+3_3801+6dup` | four intronic bases inserted; exon 8 untouched, donor `GT` intact |
| MLH1 `3:37047516 TAGG>T` at an acceptor | `splice_acceptor_variant` | `c.1732_1734del` | the exon's own first three bases (`GAG`) deleted, in frame; the acceptor `AG` intact |
| FOXP2 `7:114631524 CCAG>C` | `splice_acceptor_variant` | `c.618_620del` | one unit of a `CAG` repeat that runs across the acceptor, in frame, inside the exon; the acceptor `AG` intact |
| NBEAL2 `3:47005320 C>CGGTGCGGCCCAGGGGCT` | `frameshift_variant`, anchored inside the exon | `c.6560+4_6560+20dup` | seventeen intronic bases inserted |

The consequence terms are unchanged: they are still what fastVEP's own predictor reports, and this does not touch them.
What changed is what the criteria ask.
PVS1 asks whether the *edited sequence* loses function, and a term that describes how the record was anchored cannot answer that.

Each reported (transcript, allele) row now carries a placement, worked out after `--pick` for the rows that survive it.
The alleles are clipped of the bases the pair repeats, as the HGVSc is, shifted 3' over the reference and bounded by the transcript, and the shifted span is classed against the transcript's exons and introns.
An insertion is the two bases it sits between.
A substitution has one way of being written and is never walked, and without a reference sequence every change is placed where it sits.

| HGVS places the change | PVS1 | PS1, splice path | PM1, PM4, BP3 | BP4 splice path |
|---|---|---|---|---|
| wholly inside an intron, beyond ±2 | does not apply | does not apply | no residue for PM1, no length change for PM4 and BP3 | open, as for any intronic variant |
| wholly inside an exon, in frame, and the record's only null term is a splice one | the splice track does not apply | does not apply | unchanged | unchanged |
| wholly inside an exon, out of frame | the record's terms decide | the record's terms decide | unchanged | unchanged |
| across an exon and an intron, or across an exon between two introns | the record's terms decide | the record's terms decide | unchanged | unchanged |
| within ±2 of an exon, or not worked out | the record's terms decide | the record's terms decide | unchanged | unchanged |

A coding-null term beside the splice term keeps PVS1 in the second row, because the placement says nothing about whether that term is real.
PS1's splice path asks whether the variant abolishes a canonical splice dinucleotide, and it stands down for the placements PVS1's splice track stands down for, so that a row with only a splice term is never a splice variant to one criterion and not to the other.
It stands down whatever else the row carries: a coding-null term beside the splice one keeps PVS1 for that term, and does not keep PS1, which compares the dinucleotide the placement shows intact.
A change wholly inside an exon and out of frame is left to the record's terms by both.
It leaves the dinucleotide intact as well, but it is a frameshift whose record anchors on the splice site.
The 46 acceptor deletions and insertions of the benchmark that HGVS places in the exon out of frame, NF1 `c.4174del` among them, are PVS1 variants that also carry `PS1_Supporting`.
Taking PS1 alone away from them moved 32 calls from Pathogenic to Likely pathogenic, 29 of them on variants ClinVar calls Pathogenic, for a distinction PVS1 does not make, and it was not kept.
Each stand-down is recorded on the criterion: `details.hgvs_placement` carries the placement, and the summary says why the record's term was set aside.
This replaces an earlier rule that re-graded a splice term to a co-occurring frameshift one when the offset was deep.
That term describes the same anchor, so there was nothing to re-grade to.

The placement walks the reference, so it is taken only for an allele whose consequence terms a criterion reads it for (`fastvep_classification::reads_hgvs_placement`).
The intronic offset, which several criteria read, and the exon edge, which BP7 reads for a synonymous change, are worked out for every classified row from the transcript's exon list in place.
A counting-allocator test pins that neither allocates, and that the walk costs a handful of allocations whatever distance the change slides.
Both pipelines call one `splice_context` in `fastvep-annotate`, and look each surviving row's transcript up by the identity of its id, so two models sharing an id are told apart.

### What it measured

On the 673,660-variant benchmark, 40 variants lose PVS1 to this rule and 36 of them change call.
25 are variants ClinVar calls VUS, 5 are Benign or Likely benign, and 6 are Pathogenic or Likely pathogenic.
30 of the 36 move toward ClinVar's classification and 5 away.
The reviewer's MSH6 and DSP rows are among them, and so is FOXP2.

The six Pathogenic and Likely pathogenic variants are the price of the rule, and it is not zero:

| Variant | HGVS | ClinVar | Call before, after |
|---|---|---|---|
| ALPL | `c.69_74del`, in-frame, at an acceptor | Pathogenic/Likely pathogenic | LP, VUS |
| KCNQ1 | `c.478_480del`, in-frame, at an acceptor | Likely pathogenic | P, VUS |
| RYR1 | `c.958_960del`, in-frame, at an acceptor | Likely pathogenic | LP, VUS |
| LDLR | `c.190+2_190+3dup` | Pathogenic/Likely pathogenic | LP, VUS |
| DYSF | `c.1267_1276+4dup` | Pathogenic | LP, VUS |
| DYSF | `c.1272_1276+3dup` | Likely pathogenic | LP, VUS |

In each, the edit leaves the exon and the canonical dinucleotide in place, so PVS1's premise does not hold.
ClinVar's classification rests on evidence that an annotator cannot compute from the variant, such as functional, case and segregation data.
They now report VUS with the stand-down recorded in PVS1's `details.hgvs_placement`, which is where a curator should look.

PS1's splice path, PM1 and BP3 follow the same placement, and BP3 now records why it set an in-frame term aside.
Adding them changed no call in the benchmark.
It changed the criteria of 43 variants, each of which lost a `PS1_Supporting` that borrowed a Pathogenic variant's evidence for a dinucleotide HGVS places the change away from (`c.1136+3del` and `c.3556+3_3556+7dup` are two), and they keep their calls: 42 are VUS and one is Benign.

### BP7 and the exon edge

Walker 2023 excludes BP7 for a synonymous change at an exon's first base or its last three.
The criterion read `at_exon_edge`, this document described the exclusion, and the pipeline never set the field, so the exclusion did not exist.
It is now read off the transcript's exon coordinates, and an exon's outer end, which has no splice site beside it, is not an edge.

706 variants lose BP7 and 475 of them change call: 330 are variants ClinVar calls Benign or Likely benign, 143 are VUS, and 2 are Pathogenic or Likely pathogenic.
145 move toward ClinVar's classification and 330 away.
That is the measured price of following the guideline, and it is paid in benign calls that become VUS.
All 475 of those sit at an exon edge in the Ensembl GFF3, which was checked independently of fastVEP's own exon code.
Of 89,213 synonymous variants that still carry BP7, a random 39,999 were checked the same way and none is at an edge.

BP7's deep-intronic extension reads where the record sits, and not where HGVS places the change.
That keeps it from calling a variant deep-intronic while the row beside it says `splice_donor_region_variant`, and it is the one place the two offsets disagree in the less cautious direction.
141 of the 1,267 `intron_variant` indels that carry BP7 are numbered by HGVS inside the splice region while the record's own offset is outside it.
ClinVar calls 116 of them Benign, 17 Likely benign, 8 VUS and none Pathogenic or Likely pathogenic, so requiring both offsets to clear the region would give up 133 benign calls ClinVar agrees with and recover no pathogenic call.
That is why it was left.

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

The `appris` tier is only as good as the gene model.
APPRIS reaches fastVEP through GENCODE's `tag=appris_principal_1` / `tag=appris_alternative_2`
attributes; **Ensembl's own GFF3 carries no APPRIS tag at all**, so on an Ensembl model every
transcript scores the same at that tier and the pick falls through to `tsl`.
That is the same answer VEP gives from the same file, but it is worth knowing which tier is
actually deciding: at a locus where two genes' canonical transcripts overlap and neither is MANE,
an Ensembl model decides the pick on TSL, and a GENCODE model decides it on APPRIS.

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

**That advice is for reporting consequences, and it is not free for an ACMG call.**
Run over the whole ClinVar 2-star+ benchmark, the classifier gives:

| Run | Exact | Same-direction | Opposite | P same | LP same | LB same | B same |
|---|---:|---:|---:|---:|---:|---:|---:|
| default order | 425,869 | 522,886 | 49 | 51,855 | 6,755 | 65,461 | 111,745 |
| MANE Plus Clinical first | 425,645 | 522,484 | 47 | 51,689 | 6,655 | 65,355 | 111,722 |
| severity (`rank`) first | 423,481 | 519,120 | 452 | 52,821 | 6,895 | 63,171 | 109,754 |

Severity first gains 1,106 Pathogenic and Likely pathogenic same-direction calls and loses 4,281 benign-direction ones, and opposite-direction calls rise from 49 to 452.
It picks the transcript on which the variant looks worst, and an ACMG call made there is a call about the worst reading and not about the variant.
MANE Plus Clinical first resolves the reviewer's DST row (`c.1779-222C>T` is intronic on MANE Select and `c.64C>T p.Arg22Ter` on the Plus Clinical transcript) and costs 224 exact matches and 402 same-direction ones overall.
A Plus Clinical transcript is added for variants the Select one cannot report (as we read MANE's definition), and is not a better reading of every variant in its gene.
The supportable form is to report its call beside the MANE Select one when they differ, which is not implemented.

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
of this gate blocked on absence, and run v10 measured the cost: **1,497 ClinVar-Pathogenic null
variants lost PVS1**, in genes including SPAST, ABCB11, FLG and LAMB3 - all of which cause disease
and none of which are in the file. Worse, the genes the gate was written to catch (RYK, GIGYF2) are
absent from ClinGen too, so blocking on absence was not selecting for gene-disease validity at all.
It was selecting for curation coverage.

The consequence, stated plainly: fastVEP cannot automatically decline PVS1 for a gene a curator
believes is not disease-associated unless ClinGen has said so. RYK and GIGYF2 come back as
reference-discordant pathogenic in v11 for exactly this reason. `gene_overrides.<GENE>.disabled_criteria` is the
manual route.

## Which combining rules apply

Two schemes exist for turning met criteria into a classification, and they are encodings of the
same Bayesian model: the Richards 2015 combining table, and the ClinGen SVI point system
(Tavtigian et al. 2020, *Hum Mutat* 41:1734-1737).
fastVEP uses **the point system for the pathogenic direction and the 2015 table for the benign
direction**, which is not a compromise but a measured choice.

Points are 1 / 2 / 4 / 8 for Supporting / Moderate / Strong / Very Strong, summed over met
pathogenic criteria: `>= 10` Pathogenic, `6..=9` Likely Pathogenic.

**Why points on the pathogenic side.** The table has a documented gap at a lone PVS1. It scores 8
points, squarely inside Likely Pathogenic, but matches no row of Table 5, so the table returns
Uncertain Significance. Run v10 put **2,319 ClinVar-Pathogenic variants in VUS on `PVS1` and nothing
else** - nonsense and frameshift variants in haploinsufficient disease genes, which is close to the
least uncertain thing a classifier sees.

**Why the table on the benign side.** Tavtigian's Likely Benign band opens at `-1`, so a single BP4
is enough for a benign call, where Richards requires two benign supporting criteria or a strong plus
a supporting. Measured on the benchmark:

| scheme | pathogenic recall | benign recall | reference-discordant benign |
|---|---:|---:|---:|
| 2015 table, both sides | 56.8 % | 56.3 % | 1 |
| full point system | 63.6 % | 71.7 % | **36** |
| full points, benign floor tightened to −2 | 63.6 % | **45.6 %** | 5 |
| **points pathogenic + table benign** | **63.6 %** | **56.3 %** | **1** |

The full point system buys 15 pp of benign recall for 36 reference-discordant benign calls, and **22 of those 36
are a lone BP4**. Tightening the band instead removes them but costs more benign recall than the
table had, because a lone BP4 carries a great many correct benign calls too. Using each scheme where
it is the safer one takes the pathogenic gain and pays nothing: if ClinVar is right, a benign call on a pathogenic variant and a variant
left uncertain are not comparable errors.

Set `use_point_system = false` for the 2015 table on both sides.

## Required Data Sources

ACMG classification draws on multiple supplementary annotation (SA) sources. Place `.osa`/`.osa2` (allele-level) and `.oga` (gene-level) files in the SA directory:

### Allele-Level Sources (`.osa` / `.osa2`)

| Source | SA Key | Used By | Description |
|---|---|---|---|
| **gnomAD** | `gnomad` | BA1, BS1, BS2, PM2 | Per-population allele frequencies + AN + homozygote counts (BA1/BS1 max-pop AF; BA1/BS1 require AN ≥ 2,000) |
| **ClinVar** | `clinvar` | PM2 (frequency backstop), BA1, BS1 and BS2 (low-penetrance gate), PM3 and BP2 (companion lookup); PP5, BP6 (off by default per SVI); PS4 only when `use_clinvar_stars_as_ps4_proxy = true` | Clinical significance, review status, phenotypes, population frequencies (ExAC, 1000 Genomes, ESP); companion-variant lookups for PM3 / BP2 |
| **REVEL** | `revel` | PP3 (missense), BP4 (missense, including Very Strong band ≤ 0.003) | Missense pathogenicity score (0-1); not applied to non-missense per Pejaver 2022 |
| **SpliceAI** | `spliceai` | PP3 (caps at Supporting per Walker 2023), BP4 (≤ 0.1 → Supporting), BP7 | Splice site delta scores |
| **dbNSFP** | `dbnsfp` | Transparency-only (SIFT / PolyPhen surfaced in `details`) | The pre-PR1 ≥3-of-4 consensus path was removed per Pejaver 2022. |

The 1000 Genomes (`onekg`), TOPMed (`topmed`), AlphaMissense, PrimateAI, dbSNP and COSMIC databases can be built and loaded, and their values appear in the annotation, but no criterion reads them ([ACMG_SETUP.md](ACMG_SETUP.md#what-is-deliberately-not-in-the-stack)).

### Positional Sources (`.osa` / `.osa2`)

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

### Missing sources

A run that lacks a source of the reference stack is not an error, and it is not silent either.
The command line prints, once and on stderr, which of the nine sources are not loaded and which criteria read each of them (`warning: --acmg is running on an incomplete reference stack: ...`).
The output files are not changed by it, so a script that discards stderr does not see it.
A criterion that reads a source that is not loaded reports `evaluated: false`.
PM2 used to be the exception: it read a variant with no gnomAD record as absent from controls whether or not a gnomAD database was loaded, so a run without one met PM2 for every variant it classified.
It now declines when no gnomAD database is loaded, whatever `pm2_absent_when_no_record` says, and so do BA1, BS1 and BS2.
That setting is about a database that is loaded and holds no record of the variant (limitation 30), which is a different statement.

A source can also be loaded and hold only some of the genome, and for gnomAD that is the same mistake.
A database built for some chromosomes (one shard of a per-chromosome build, a region-restricted build, a chromosome that failed to build and was not noticed) opens without complaint and says "no record" for every variant of the others, which PM2 reads as absence from the population.
A deployment script of this project once built a gnomAD of chromosome X alone.
The classifier is told, for each variant, whether the gnomAD databases that are loaded hold its chromosome, and declines PM2, BA1, BS1 and BS2 for a variant on one that they do not (`holds no record of chromosome 1`), whatever the genome; chromosomes are compared without their `chr` prefix, and the shards of one source count as the union of what they hold.
The command line says which chromosomes of the gene model the loaded gnomAD lacks (`warning: --acmg: the gnomAD databases that are loaded hold nothing of chromosomes ...`), and the web server counts a gnomAD that lacks a chromosome of the gene model's autosomes or X as missing, so the genome is not classified until it is built (Y and the mitochondrial genome are left out of that, because the gnomAD of GRCh37 has no Y: a variant there is declined one by one).
A database that cannot say which chromosomes it holds is read as holding all of them.

`fastvep-web` goes further than a warning, because the people who read its answers are not the ones who built its databases: it refuses to classify a genome that lacks any of the nine sources, and `GET /api/status` and `GET /api/genomes` say which genomes can be classified.
The nine sources are named once, in `crates/fastvep-classification/src/sources.rs` (`REFERENCE_STACK`), with the criteria that read each of them and what a classification loses without it.
The command line, the web server and this document all read from that list.
A test classifies a set of inputs with each source and without it and fails if a criterion that is not listed under a source changes when the source is removed, so the list cannot drift from the criteria.
It found one that the documents had missed: ClinVar's low-penetrance gate also withholds BA1.

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

**What it is worth.** Measured over the full ClinVar 2-star-or-better set when it was
switched on (v15 in `RUN_VERSIONS.md`), BP3 fired on 1,307 variants, split 688 on
ClinVar-benign variants against 56 on ClinVar-pathogenic ones - directionally right,
roughly 12:1; on v37 it fires on 1,296, 684 against 55. Its effect on the headline
was nil: benign recall moved 71.47 % to 71.49 %, and every error count was
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
| **PVS1** | Very Strong → Supporting* | Null variant in LOF-intolerant gene; graded per Abou Tayoun 2018 decision tree (NMD prediction, %protein removed, critical region, alt start codon, last exon). Two gene-level preconditions run first: the gene must have an established gene-disease relationship, and its curated mechanism must not exclude loss of function. A change HGVS places wholly inside an intron beyond ±2 is not a null variant whatever term its record's anchor produced, and neither is an in-frame exonic change under a splice-only term - see [Where HGVS places an indel](#where-hgvs-places-an-indel) | Consequence + HGVS placement + pLI/LOEUF/ClinGen GDV + curated mechanism + transcript context | Yes |
| **PS1** | Strong | Same AA change as established pathogenic missense **or** canonical ±1/2 splice predicted to produce same RNA outcome (Walker 2023) | ClinVar protein index + splice catalog | Yes (with .oga) |
| **PS2** | Strong | De novo with confirmed parents | Trio VCF genotypes | Yes (with trio) |
| **PS3** | Strong | Functional studies show damaging. Read from `--functional-evidence`; suppresses PP3/BP4/BP7 for that variant (SVI ranks experimental evidence above prediction) | Curated TSV | Yes (with file) |
| **PS4** | Strong | Prevalence in affected >> controls | Case-control statistics (NotEvaluated by default; ClinVar review-stars proxy is invalid per SVI, opt in via `use_clinvar_stars_as_ps4_proxy`) | No (Partial in proxy mode) |
| **PM1** | Moderate | Mutational hotspot / critical domain **without benign variation** - both halves of Richards 2015's wording are tested, so a window that also holds benign ClinVar missense does not count. NotEvaluated where the gene has no established gene-disease relationship: a cluster of assertions in an uncurated gene is not evidence of a critical region | ClinVar protein density (both directions) + ClinGen GDV | Yes (with .oga) |
| **PM2** | Supporting* | Absent/rare in population databases - inheritance-aware: AD/unknown fires at AF ≤ 0.00004 (published VCEP practice and a sweep; see "Choosing the PM2 threshold for dominant genes"), recessive-only at AF ≤ 0.00007 (equal to the Hearing Loss VCEP's recessive threshold, swept for the ACMG manuscript and flat; SVI v1.0 lowers PM2 to Supporting and gives no frequency). NotEvaluated wherever BA1/BS1/BS2 are: absence is only evidence when the database could have seen the variant | gnomAD AF + QC flags + inheritance (ClinGen's curated mode, or the words of an OMIM phenotype) | Yes |
| **PM3** | Supporting → Very Strong* | In trans with pathogenic (recessive); points-based per SVI v1.0 (in-trans/P=1.0pt, in-trans/LP=0.5pt, unphased/P=0.5pt, unphased/LP=0.25pt, hom=0.5pt cap 1.0) | Phased VCF + ClinVar | Yes (with trio) |
| **PM4** | Moderate | Protein length change (in-frame/stop-loss). Not applied where HGVS places the change wholly inside an intron beyond ±2: no exon base is edited, so the protein is unchanged | Consequence type + HGVS placement | Yes |
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
| **BA1** | Standalone | AF > 5%, and AN ≥ 2,000 at the site (the BA1 recommendation of Ghosh 2018 asks for a dataset of at least 2,000 observed alleles; applying it to the site's AN is ours, and the March 2024 gnomAD v4 communication states no floor). Tested against gnomAD's **grpmax filtering allele frequency** (the lower bound of the 95% CI of the allele frequency in the genetic ancestry group where it is highest, as the ClinGen VCEP Review Committee specifies for BA1 and BS1; the filtering allele frequency is the concept of Whiffin 2017) where the database provides it, else the population maximum. The curated exception list blocks BA1 regardless of AF - Ghosh 2018's nine, three hypomorphic alleles and three the round-6 review named; see [Curated frequency-exception list](#curated-frequency-exception-list) | gnomAD population AFs + FAF + AN + HGVSc | Yes |
| **BS1** | Strong | Cross-population AF > expected, read from the same statistic as BA1 (filtering AF where available, else population maximum). NotEvaluated for homology-confounded genes, for gnomAD-flagged `segdup`/`lcr` sites, for non-PASS gnomAD records, and for ClinVar low-penetrance / risk alleles. Stands down where the gene's own BA1 threshold already covers the frequency. Blocked outright by a curated exception whose `blocks` names BS1 | gnomAD per-population AFs + FAF + QC flags + ClinVar terms | Yes |
| **BS2** | Strong | Observed in healthy adult, for a disorder with full penetrance expected at an early age (Richards 2015). Both branches ask whether more healthy individuals have the allele than the disorder could account for, by holding a 95% lower bound on their frequency against a prevalence threshold, so the test scales with cohort size. **Dominant** genes count carriers (`AC - nhomalt`, a homozygote once) and require AC ≥ `bs2_ad_min_ac` (5) **and** a lower bound above `bs2_ad_prevalence_threshold` (default 3e-4, set by measurement). **Recessive / X-linked** genes count individuals with no functional copy (homozygotes **plus hemizygotes** on a non-PAR sex-chromosome site) and require ≥ `bs2_ar_min_hom` of them **and** a lower bound above `bs2_hom_prevalence_threshold` (1e-3). Both thresholds are set by measurement - see [Choosing the BS2 prevalence threshold](#choosing-the-bs2-prevalence-threshold). A gene is dominant when an established ClinGen curation of it is tagged `AD` or `SD`, the curated mode of inheritance, or, for an OMIM phenotype, when its words say so (limitation 20 has the history); a gene that is both, or neither, takes the recessive / X-linked arm. The dominant arm's carrier count is not added beside BS1 or BA1, which read the same gnomAD frequency (the last rule of the reconciliation pass). Not evaluated where gnomAD AN is missing, since there is then no cohort size; it does not take BA1's and BS1's minimum AN, which we read as SVI guidance about the noise of a frequency estimate, and the lower bound already accounts for a small cohort. Blocked outright by a curated exception whose `blocks` names BS2, which is how a hypomorphic allele's unaffected homozygotes are kept from reading as tolerance | gnomAD AC + hom + hemizygote counts + AN + inheritance | Yes |
| **BS3** | Strong | Functional studies show no damage. Read from `--functional-evidence`; suppresses PP3/BP4/BP7 for that variant | Curated TSV | Yes (with file) |
| **BS4** | Strong | Lack of segregation | Pedigree data | No |
| **BP1** | Supporting | Missense in a gene where primarily truncating variants cause disease. Requires both the constraint signature (high pLI, low misZ) **and** a mutation spectrum without an established missense mechanism (< `bp1_max_pathogenic_missense` pathogenic missense in the ClinVar protein index) | pLI + misZ + ClinVar protein index | Yes |
| **BP2** | Supporting | In trans/cis with pathogenic | Phased VCF + ClinVar | Yes (with trio) |
| **BP3** | Supporting | In-frame indel in repeat region. Not applied where HGVS places the change wholly inside an intron beyond ±2 | Consequence + HGVS placement + RepeatMasker | Yes (with .osi) |
| **BP4** | Supporting-Very Strong | Computational evidence (benign) - REVEL **missense-only** with Very Strong band at ≤ 0.003 (Pejaver 2022) + SpliceAI ≤ 0.1 → Supporting (Walker 2023). BP4-splice is withheld from PVS1-territory consequences, from deep-exonic missense, and from in-frame indels / stop-lost / protein-altering variants, for which no calibrated benign predictor exists - except where HGVS places the change wholly inside an intron beyond ±2, which makes it as splice-assessable as any other intronic variant. A splice BP4 is also set aside beside a PP3 met from REVEL (reconciliation rule below). | REVEL + SpliceAI + HGVS placement | Yes |
| **BP5** | Supporting | Alternate molecular basis in case | Case-level analysis | No |
| **BP6** | Supporting | Reputable source reports benign | ClinVar (disabled by default per SVI) | Partial |
| **BP7** | Supporting | Synonymous + no splice impact + not conserved. Per Walker 2023: must NOT fire for synonymous at first base / last 3 bases of an exon (`at_exon_edge`, read off the transcript's exon coordinates; an exon's outer end, which has no splice site beside it, is not an edge); extends to deep-intronic offsets ≥ 7 (donor) or ≤ -21 (acceptor). | Consequence + SpliceAI + PhyloP + exon position | Yes |

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

### PVS1 Strength Grading (after Abou Tayoun 2018)

PVS1 is graded by a decision tree over null-variant context, after the published one. The output code carries the strength suffix (e.g. `PVS1_Moderate`).

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

**Where this table departs from the published tree** (Figure 1 of Abou Tayoun 2018, read against the paper in October 2026):
the tree grades an NMD-escaping truncation outside a critical region `PVS1_Strong` when more than 10 % of the protein is removed and `PVS1_Moderate` when less is, so the two rows above for that case are one step lower than the paper's;
and it has no last-exon branch for a canonical splice variant, grading by critical region and by the fraction of protein removed (Strong, Strong, Moderate), so the `PVS1_Moderate` row for it is ours.
The tree has no rule for a truncation in the first coding exon, where translation can re-initiate downstream, and this classifier has none either.
Regrading the first two would change 7 of 673,660 calls on the ClinVar 2-star+ benchmark (a replay; `manuscript/acmg`, Section 3.2), so it has not been done.

#### Which NMD prediction the tree runs on

Abou Tayoun 2018 predicts nonsense-mediated decay with the **50-nucleotide rule**: a premature termination codon triggers decay only when it sits more than 50 nt upstream of the final exon-exon junction.
fastVEP can compute that exactly, and also carries the coarser **last-exon proxy** it has always used.
The two disagree in one place - a PTC in the 3' end of the penultimate exon, where the proxy claims decay and the measurement finds escape.

`pvs1_nmd_50nt_rule` chooses between them, and it defaults to `false` (the proxy) on a measurement rather than a preference.

On a 4,000-variant sample of PVS1-carrying ClinVar 2-star records, switching the exact rule on moved 58 calls from Likely Pathogenic to Uncertain: **54 of them ClinVar calls Pathogenic**, 4 Uncertain, none Benign.
Extrapolated to the full PVS1 population that is roughly 825 correct pathogenic calls given up to remove a single reference-discordant pathogenic call.
On the full benchmark, in a run of the final classifier with only this setting changed (`manuscript/acmg`, Section 3.4), the exact rule loses 603 calls in direction and 0.73 points of pathogenic-side recall to avoid one reference-discordant pathogenic call.

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

In-cis companions are excluded (they count toward BP2 instead).
Requires recessive inheritance in the gene's disease annotation: ClinGen's curated `AR` or `SD` on an established curation, or the words of an OMIM phenotype.

### Anti-Double-Counting (Reconciliation Pass)

After per-criterion evaluation, a reconciliation pass suppresses a criterion that counts a signal another criterion has already counted.
Suppressed criteria appear with `met: false` and `details.suppressed_by_reconcile` explaining why.
Each rule is a function of the criteria that fired and nothing else; the BS2 rule and the PS1 upgrade run first, and touch no criterion that the other rules read.

| Trigger | Action | Reference |
|---|---|---|
| PVS1 fires + PM1 fires | Suppress PM1 | Abou Tayoun 2018 (PVS1 already counts loss of the region; PM1 is residue-level evidence for missense) |
| PS3 or BS3 fires (curated functional evidence) | Suppress PP3, BP4 and BP7 | ClinGen SVI: experimental evidence outranks prediction, in either direction |
| PVS1 fires + PP3 driven by SpliceAI | Suppress PP3 | Walker 2023 (PVS1 already counts splice signal) |
| PS1 fires + PP3 driven by REVEL (missense) | Suppress PP3 | Pejaver 2022 (PS1 covers residue) |
| PM5 fires + PP3 driven by REVEL (missense) | Suppress PP3 | Pejaver 2022 (PM5 covers residue) |
| PP3 elevated to Strong + PM1 fires | Suppress PM1 (cap combined at Strong = 4 Tavtigian points) | Pejaver 2022 |
| PP3 driven by REVEL (missense) + BP4 driven by SpliceAI | Suppress BP4 | Walker 2023, as quoted from the medRxiv version because the published text could not be checked here (the protein and splicing effects are assessed separately and the more deleterious stands; SpliceAI saying splicing is unaffected is not a vote against a damaging missense prediction) |
| BS1 or BA1 fires + BS2 decided on carriers (a dominant gene) | Suppress BS2 | One gnomAD observation counted twice: a carrier count is the allele count that BS1 and BA1 turn into a frequency, and two Strong benign criteria read as Benign |

One effect is an upgrade and not a suppression: PS1's splice path is graded against the PVS1 the variant already carries (Walker 2023, Table 3).

**BS2's carriers beside BS1 and BA1.**
BS2's dominant arm counts the individuals who carry the allele, `AC - nhomalt`, and tests a lower bound on their frequency against `bs2_ad_prevalence_threshold`, 3e-4.
That is the same gnomAD observation that BS1 and BA1 test against thresholds of their own, and the thresholds are far apart: BS1's is 5e-3 on a filtering allele frequency, and the carrier threshold is 3e-4 on the frequency of carriers, which is about 1.5e-4 as an allele frequency.
In the genes that ClinGen curates as autosomal dominant only, the arm fires on 7,234 of the 7,328 variants that meet BS1 and on all 8,869 that meet BA1, so it adds a second Strong criterion for the same observation nearly everywhere it can.
The first time the classifier read ClinGen's curated inheritance and so reached the carrier arm for the genes it had missed, that took two ClinVar Pathogenic variants to Benign: TTR p.Val142Ile (722 gnomAD alleles, 3 of them homozygous) and VWF p.Arg854Gln (7,432 alleles), each at VUS with BS1 and each Benign with BS1 and BS2.
With the arm withheld beside BS1 and BA1 both stay at VUS, and the arm still acts below BS1's threshold, which is the only place it adds anything.
Compared directly, the unguarded and the adopted reading differ on 6,494 calls, every one of them Benign in the unguarded reading: 5,192 that the adopted one calls Likely benign and 1,302 that it leaves at VUS, because under the Richards rules two Strong benign criteria are Benign and one Strong benign criterion alone is not a call.
That is the cost in exact agreement of not counting one frequency twice, and no ClinVar pathogenic-side call is among them.
The recessive arm, which counts individuals with no working copy, is not reconciled with BS1: a homozygote is a genotype that an allele frequency does not give.
The hemizygotes that the arm counts on the non-pseudoautosomal X are part of the allele count, so there the same doubling remains (limitation 29).

**What the pass does not do.**
PP2, PM1 and PP3 on one missense variant are counted as three independent criteria, as the point system counts them, except where the PP3_Strong rule above removes PM1, although they can be correlated.
That is limitation 27, which has what a combined cap at Strong would do.

### Genes with more than one mode of inheritance

ClinGen's guidance on genes associated with several disorders (version 1, January 2024) is to classify a variant against each disease on its own, because a variant can be benign for one and pathogenic for another.
The criteria that depend on the disease are the ones that branch on its mode of inheritance (PM2's frequency threshold, which arm of BS2 applies, whether PM3 can apply, BP2) and on its validity (whether PVS1, PP2 and PM1 are withheld).
A gene with an autosomal dominant disease and an autosomal recessive one is two questions to them, and a single call has to read the gene's diseases together.
PM2's threshold, BS2's arm and BP2's trans arm take the reading for a gene that is neither dominant only nor recessive only, which is the cautious one for each and answers neither disease.
PM3 asks only that one established disease be recessive, because the observation it counts, a variant in trans with a pathogenic one, is evidence about the variant whatever else the gene causes.
BP2's trans arm is the argument of a fully penetrant dominant disease, which a recessive disease of the same gene explains away, so it needs dominant and none recessive.

When the established diseases of a gene fall under more than one mode of inheritance, the variant is classified once more for each mode, with the gene's disease annotation cut down to the established diseases that have it.
Nothing else about the input changes: the variant, its frequencies and the gene's constraint are not disease-specific.
The result keeps the call for the gene as a whole exactly as it was, and adds the list.

In JSON:

```json
"acmg": {
  "classification": "Likely_pathogenic", "shorthand": "LP", "criteria": [ ... ],
  "per_disorder": [
    { "inheritance": "AD",
      "disorders": [ { "label": "BRCA1-related cancer predisposition", "mondo": "MONDO:0700268", "validity": "Definitive" } ],
      "classification": "Likely_benign", "shorthand": "LB", "criteria": [ ... ], "counts": { ... } },
    { "inheritance": "AR", "disorders": [ ... ], "classification": "VUS", "shorthand": "VUS", "criteria": [ ... ], "counts": { ... } }
  ]
}
```

In a VCF run with `--acmg` the last CSQ column is `ACMG_DISORDERS`, one `INHERITANCE:CLASS:CRITERIA` entry for each mode joined by `&`, with the criteria that were met joined by `+`: `AD:LB:BS2+BP4&AR:VUS:BP4`.
It is empty for the gene that has one mode, which is nearly all of them, and it is not part of the 49 default fields, so a run without `--acmg` has the columns it always had.
The web evidence modal lists the same entries under "By mode of inheritance".

`inheritance` is ClinGen's code (`AD`, `AR`, `XL`, `SD`, `MT`, `UD`), or `unspecified` for the entries that carry no curated mode: a curation whose code ClinGen left blank or this version does not know, and the untagged phenotypes of an OMIM `genemap2.txt` source.
The untagged phenotypes of an OMIM gene all share one context, so such a gene has no list of its own, and what the classifier made of it is the call for the gene as a whole, which reads each phenotype by its words.
`SD` (semidominant) reads as both dominant and recessive, and `XL`, `MT` and `UD` as neither, so they take the arm for unknown inheritance.
Only established diseases make a mode: a gene whose only curation is Limited, Disputed, Refuted or No Known Disease Relationship has none, and is withheld from PVS1, PP2 and PM1 as before.

On the ClinVar 2-star+ benchmark, 128 genes have more than one mode, and 87,599 of the 673,278 scored variants (13.0%) are in them.
For 482 of those variants the modes give different classes, and for 435 of them different directions; of the 20,754 ClinVar pathogenic-side variants among them, 15 differ in class and 14 in direction.
The effect on a benchmark is small because only two criteria ever differ between the modes there: BS2, for 2,976 variants, and PM2, for 1,461 (PM3 and BP2 need a proband, which the benchmark lacks), and the criteria that the modes meet differ at all for 4,437 of the 87,599 variants (5.1%).
For a variant in ABCC8, ATM, BRCA1, BRCA2 or PKD1 that is the call to read before reporting.

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
- **BP2** (Supporting): Variant is in cis with a ClinVar pathogenic variant, or in trans with pathogenic in a gene whose established diseases are dominant and none recessive.
  Requires phased data.

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

A run with `--acmg` has a third and last field, `ACMG_DISORDERS`, which is empty except for a gene whose established diseases have more than one mode of inheritance: `AD:LB:BS2+BP4&AR:VUS:BP4` (see [Genes with more than one mode of inheritance](#genes-with-more-than-one-mode-of-inheritance)).
It is not one of the 49 default fields, so a run without `--acmg` has the columns it always had.

### TSV Output

The tab output is VEP's tab layout and has no ACMG column: a run with `--acmg --output-format tab` writes the columns that a run without `--acmg` writes.
Use the VCF, where `ACMG`, `ACMG_CRITERIA` and `ACMG_DISORDERS` are in the CSQ field, or JSON, for the classification.

## Architecture

### Crate: `fastvep-classification`

```
crates/fastvep-classification/src/
  lib.rs              # Public API: classify(), extract_classification_input()
  types.rs            # EvidenceStrength, EvidenceCriterion, AcmgClassification, AcmgResult
  sa_extract.rs       # ClassificationInput, typed SA deserialization, GenotypeInfo, CompanionVariant
  disease.rs          # a gene's disease annotation, read: ClinGen's classification, curated inheritance and MONDO id; modes of inheritance
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
| `crates/fastvep-cli/src/pipeline/annotate.rs` | CLI batch: same pipeline with parallel variant processing |
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
8. **BP7 exon-edge / deep-intronic extension** uses optional `at_exon_edge` / `intronic_offset` fields. Both pipelines populate them from the transcript's exon coordinates; a caller of the classification crate that does not gets the legacy synonymous-only rule.
9. **BP3** requires a RepeatMasker interval `.osi`. Build it from the UCSC track - see [Building the RepeatMasker track](#building-the-repeatmasker-track-bp3). Without one BP3 reports `evaluated: false` rather than assuming a variant is not in a repeat.
10. **PS2/PM6/PM3/BP2** require a multi-sample VCF with trio sample names configured
11. Compound heterozygote detection (PM3/BP2) works per-batch in the CLI; variants in different batches within the same gene may not be cross-referenced
12. **A gene with several disorders is classified once for each mode of inheritance, not once for each disorder, and no disorder is chosen for you.**
    A gene whose established ClinGen curations fall under more than one mode of inheritance (128 genes of the ClinVar benchmark, 87,599 of its 673,278 scored variants) is also classified against each, and the result carries a `per_disorder` list: in JSON an array of entries, each naming the diseases it was made against with their MONDO identifiers, and in a VCF run with `--acmg` a last CSQ column, `ACMG_DISORDERS` (see [Genes with more than one mode of inheritance](#genes-with-more-than-one-mode-of-inheritance)).
    The call for the gene as a whole reads the modes together, which is the cautious reading for PM2's threshold, BS2's arm and BP2's trans arm (PM3 asks only that one established disease be recessive), and is unchanged by the list.
    The modes are grouped by inheritance and not by disease because every criterion that branches on a disease branches on its mode of inheritance (PM2's threshold, BS2's arm, PM3, BP2) or on its validity, and two diseases of one gene that share a mode are classified identically.
    Three things are not done.
    The `disorders` blocks of the generated VCEP table are still not read, because the generator writes the same thresholds under every disorder of a gene and nothing in it says which disorder a threshold was meant for.
    The `--disorder MONDO:...` selector that the round-2 plan named is not implemented, so the user reads the list and chooses.
    And a gene whose curations are all weak has no mode to classify against, which is limitation 24.
13. **The frequency dead zone between PM2 and BS1.** See below - this is the largest known gap on the frequency side and it is open by decision, not oversight.
14. **Penetrance is unsourced.** BS2's precondition is "full penetrance expected at an early age" and no public resource publishes penetrance per gene-disorder pair. Prevalence, onset and inheritance are covered; see [Per-gene frequency thresholds](#per-gene-frequency-thresholds-clingen-vcep-specifications).
15. **Published VCEP thresholds are applied without their founder-variant exceptions.** The exclusions are prose in the specification text, and `ba1_exceptions` carries no panel-specific entries beyond Ghosh 2018's nine, which is why the generated threshold table is opt-in rather than a default.
16. **The ClinVar-concordance benchmark cannot validate the BA1 exception list.** Seven of Ghosh 2018's nine variants are not in the April 2026 ClinVar 2-star+ set: six are `Conflicting classifications` in ClinVar and carry no consensus classification, and the seventh (HFE p.Cys282Tyr) was not in the April set and is in the October one; the criterion exists for exactly the population a consensus benchmark discards. It is covered by an end-to-end test over all nine instead. Any future exception entries need the same treatment.
17. **The hypomorphic-allele list is three entries long, and the reviewer-named list three more.** GAA `c.-32-13T>G`, CFTR `c.1210-11T>G` and SPTA1 `c.4339-99C>T` ship with `blocks = ["BA1", "BS1", "BS2"]`, because an allele that causes disease only in trans with a rarer null is expected to be common and expected to be tolerated homozygous. The mechanism is general and the list is not; it grows by curation only. G6PD `c.563C>T`, G6PD `c.202G>A` and SERPINC1 `c.439A>G` ship the same way on the round-6 reviewer's word, with no PMID, which each entry says. See [The curated frequency-exception list](#curated-frequency-exception-list).
18. **A variant at a site where the population data is vetoed cannot reach Likely pathogenic alone.** Where gnomAD's FILTER or a low-complexity / segmental-duplication flag makes the frequency unusable, BA1, BS1, BS2 and PM2 are all withheld while PVS1 - read from the same alignments - is not, so a lone PVS1 call would rest on a one-sided ledger. Such calls are reported Uncertain with the veto named. The rule applies only when the suppressed frequency would itself have met BA1 or BS1; applying it to every flagged site cost 2,630 correct pathogenic calls to remove 18 wrong ones.
19. **HGVS placement reaches PVS1, PS1's splice path, PM1, PM4, BP3 and BP4's splice path, and no further.**
    A placement withdraws the criteria that read a protein or splice consequence it shows to be only the record's anchor, and for a null term placed wholly in the intron it opens BP4's splice path.
    That added BP4 to 28 variants in the benchmark and changed 11 calls, 10 toward ClinVar's classification and one away, and one of the 28 is a ClinVar Likely pathogenic variant.
    An in-frame exonic change under a splice-only term loses PVS1's splice track but earns no PM4 or BP3, because the record carries no in-frame term: FOXP2's Q209del in a `CAG` repeat is an example.
    BP7 reads the record's own intronic offset and not the placed one, on purpose, and that is the one place the two disagree in the less cautious direction.
    On an acceptor the shift runs toward the exon, so an intronic indel in a pyrimidine tract that the record anchors deeper than `-21` can be numbered by HGVS inside the splice region.
    Of the 1,267 `intron_variant` indels that carry BP7, 141 are numbered that way.
    ClinVar calls 116 of them Benign, 17 Likely benign, 8 VUS and none Pathogenic, so asking both offsets to clear the region would trade 133 benign calls ClinVar agrees with for no pathogenic call recovered, and it was not done.
    PP3's splice path is not gated on placement.
    A change that spans an exon and an intron, or crosses an exon between two introns, is left to the record's terms.
    So is an exonic frameshift under a splice term, which PVS1 keeps as a null variant and PS1 keeps with it.
    On an acceptor, a deletion of the final `AG` of `...AGAG|exon` is placed on the dinucleotide by the 3' rule although an equivalent placement leaves the intron ending in `AG`, so it keeps PVS1; this was reasoned and not measured.
    Without a reference sequence every change is placed where it sits.
    At a transcript's end the same edit can be classed as spanning the transcript or as exonic, depending on how the record is anchored.
    Spanning leaves the record's terms in charge, so the cost is a stand-down that is not made; this was found in review by comparing the placement with a brute-force oracle, and not measured on the benchmark.
    Six ClinVar Pathogenic or Likely pathogenic variants fall to VUS under the rule; see [Where HGVS places an indel](#where-hgvs-places-an-indel).
20. **"Dominant" and "recessive" were read from disease names, not from ClinGen's inheritance tag, until v38.**
    `has_dominant_inheritance` matched OMIM-style text and "dominant" in the disease name, and ignored the curated mode of inheritance, one of `AD`, `AR`, `XL`, `SD`, `MT` and `UD`, in `(ClinGen Definitive/AD, MONDO:...)`.
    Of the 740 genes with a variant in the benchmark whose established ClinGen curations are all AD, 39 registered as dominant, and of the 1,285 that are all AR, 60 registered as recessive; of the 2,924 genes of the file, 50 were read as dominant-only, 67 as recessive-only, 12 as both and 2,795 as neither.
    BS2's dominant arm, PM2's recessive threshold, PM3 and BP2 share the predicate, so each reached roughly one gene in twenty.
    v38 reads the curated field: the 2,924 genes read as 746 dominant-only, 1,310 recessive-only, 169 both and 699 neither, and the genes with a benchmark variant as 740, 1,285 and 167.
    An OMIM `genemap2.txt` phenotype has no tag and is still read by its words.
    Against the commit before, on the 673,660 records, it adds 149 exact matches and 1,490 calls in the same direction as ClinVar's classification and moves 25 ClinVar pathogenic-side variants.
    Twenty-one go from Likely pathogenic to Pathogenic, because PM2's recessive threshold now reaches recessive genes that v37 did not read as recessive (PAH, DHCR7, CYP27A1 and others).
    Four go from Likely pathogenic to VUS: NEK1 p.Ser1036Ter, DDX41 p.Asp140GlyfsTer2, SRCAP p.Pro1879ThrfsTer21 and PROKR2 p.His20MetfsTer24, each a lone PVS1 before and PVS1 with BS2 now, on 295 to 324 gnomAD carriers.
    Reference-discordant benign calls (ClinVar pathogenic-side called benign-side) stay at 24, and reference-discordant pathogenic calls go from 25 to 24.
    Reading the tag for BS2 had been measured on v37 and not done, because it moved six ClinVar pathogenic-side variants toward benign: those four, and TTR p.Val142Ile and VWF p.Arg854Gln, which went from VUS to Benign.
    The last two are what the reconciliation pass now prevents: BS2's carrier count and BS1 read one frequency, and counting both was what made them Benign.
    v38 reads the field although v37 declined to.
    The unread field made the classifier answer a question about the name of a disease where the source curates the answer, which is a defect in what is read and not a threshold, and the rule that prevents the two Benign calls is stated by what it counts and not by a list of variants (it changes 6,494 calls), although it was written after the two were seen.
    The other four move to VUS and not to a benign-side call, and are accepted as the price of a carrier arm that applies to every gene read as dominant, because there is no penetrance to read (limitation 14).
    `manuscript/acmg`, Section 3.8, gives the reasons and the decision about the label guard.
    The tag cannot supply the penetrance that BS2 asks for, and ClinGen does not publish it (limitation 14).
    ClinGen's labels for NEK1 and DDX41 are a susceptibility and a predisposition syndrome, and withholding BS2's carrier arm from such labels recovers those two for 65 fewer correct benign-side calls (53 fewer calls in direction on balance), which was measured and not adopted because it reads the name of a disease again; why SRCAP and PROKR2 are affected was not established.
    `RUN_VERSIONS.md`, "v38", gives the runs, and `manuscript/acmg`, Section 3.8, the same on the variants of ClinVar's release of 4 October 2026.
21. **BP7 keeps the 2015 conservation test and takes SpliceAI below 0.2 as no impact.**
    Walker 2023 applies BP7 after BP4, which for SpliceAI is a score of 0.1 or below, and cautions against a conservation test unless its tool, depth and cut point have an empirical justification.
    `phylop_conserved` (2.0) is from the first commit with no source recorded, and nothing but BP7 reads it.
    On the ClinVar 2-star+ benchmark, switching the test off adds 16,482 correct benign-side calls and 2 reference-discordant benign calls (variants of NF1 and MAPT), and asking for SpliceAI of at most 0.1 as well removes 402 of them and no error (`manuscript/acmg`, Section 3.3).
    Neither has been changed.
22. **BP1's veto counts the variant's own ClinVar record.**
    `exclude_self_from_clinvar_evidence` discounts it for PS1 and PM1, and the veto of `bp1_max_pathogenic_missense` reads the same index without discounting it.
    Leaving the record out adds 3 reference-discordant benign calls on the benchmark (FBXO28, PCGF2 and SBF1, each a ClinVar-pathogenic missense variant whose own entry is the third of the three that the veto needs; `manuscript/acmg`, Section 3.4).
    It has not been changed.
23. **The ClinGen file is a download, and the benchmark's was incomplete.**
    The export of 30 April 2026 lacked 104 curations that ClinGen dates before that day, among them CFTR (Definitive, autosomal recessive, 2022); seven further rows of the current export only re-label curations that the April file lists, and one is dated after that file was downloaded.
    A gene with no curation is not blocked from PVS1, but with no other route to loss of function it is not licensed either, so all 558 null variants of CFTR in the benchmark were VUS.
    With the 112 rows put back, 542 of them are called pathogenic-side (`scripts/15_clingen_complete_export.py`, `manuscript/acmg`, Section 3.5).
    The check that finds this is mechanical: list the genes in which ClinVar's pathogenic stop-gained and frameshift variants carry no PVS1, and compare them with an independent listing.
24. **PVS1's gene-level gates cannot see which disease a variant is pathogenic for.**
    The validity gate withholds PVS1 where every ClinGen curation of the gene is weak, and the mechanism gate where the curated mechanism is gain of function, whatever disease the variant belongs to.
    RASA1's only curation in the file is a Disputed one for Noonan syndrome, while ClinVar names capillary malformation-arteriovenous malformation for 39 of its 40 truncating variants that carry no PVS1.
    Item 12 narrows this only for a gene that has several established curations: RASA1 and LRP2 have none that is established, so there is no mode to classify against, and what they need is the disease named by the variant's own ClinVar record, which a VCF does not carry.
25. **The validity gate used to miss a weak-only gene whose ClinGen label contains a semicolon, and no longer does.**
    The gate reads the tag of each phenotype, and the builder splits phenotypes on `;`, so NSMCE3 (Limited), whose label ends in `;`, left an untagged fragment that counted as established and the gene escaped the gate.
    The typed reading (`OmimData::entries`) drops an untagged piece from a gene that has any tagged one, since a ClinGen build tags every phenotype it writes and an OMIM build tags none.
    NSMCE3 was the only one of the 505 weak-only genes in the file of 30 April 2026 that escaped, and it has no variant in the benchmark, so no call changed.
    Two Definitive curations (MAPK8IP3 and SCAPER) contain a `;` as well and are unaffected.
    `crates/fastvep-cli/tests/clingen_real_rows.rs` covers all three, from real curations through the converter and the builder.
26. **BP1's veto count of 3 was never swept, and a stricter count calls fewer ClinVar-VUS variants benign-side.**
    The veto withholds BP1 in a gene that has at least `bp1_max_pathogenic_missense` pathogenic missense records in the ClinVar protein index, and the count was set with the veto and not varied.
    On the ClinVar 2-star+ benchmark a count of 1 gains 1,716 calls in direction (1,865 ClinVar-VUS variants no longer called benign-side, less 149 correct benign-side calls given up) and removes 2 reference-discordant benign calls (IQSEC1 and PRR12, the two BP1 calls of the residual).
    The removal of PRR12 rests on the variant's own record: the index holds one pathogenic missense record for PRR12, the variant itself.
    The removal of IQSEC1 does not at a count of 1, because the index holds a second record for it (p.Thr343Met) besides the variant, and does at a count of 2.
    With the own record left out as well (item 22) a count of 1 removes IQSEC1 only and a count of 2 neither; a count of 2 alone gains 849 calls and removes one.
    Above 3, a count of 5 adds 3 reference-discordant benign calls, 10 adds 5 and 25 adds 12 (`manuscript/acmg`, Section 3.3).
    It has not been changed.

27. **PP2, PM1 and PP3 on one missense variant are counted as three independent criteria.**
    The point system counts criteria as independent, although these three can be correlated, and the reconciliation pass removes only PM1 beside a Strong PP3.
    On the 673,278 scored variants, 1,591 have all three, 913 of them at Likely pathogenic.
    A combined cap of 4 points (Strong) on the three, by a real run, changes 910 calls, every one from Likely pathogenic to VUS: 320 variants ClinVar calls pathogenic-side, 587 it calls VUS and 3 it calls benign-side.
    It raises agreement in direction by 267 calls and exact agreement by 512, lowers recall of ClinVar's pathogenic-side variants from 62.5% to 62.1%, and avoids 3 reference-discordant pathogenic calls, giving up about 107 correct pathogenic-side calls for each.
    It was not adopted, because it changes the combining rule and does not correct a defect, because the benchmark does not separate the two readings (the 587 ClinVar-VUS variants it keeps at VUS are not shown to be errors), and because 107 calls per contradiction is at the level that the manuscript's knee rule (Section 2.5) leaves undecided.
    The run is `v38_alt_triple_cap`, built from `manuscript/acmg/replay_validation/inheritance_alternatives_env.patch`; on the October release it changes 1,018 calls and gives up 354 correct pathogenic-side calls for 4 contradictions avoided, and on the new variants 37 for 1 (`manuscript/acmg`, Section 3.8).

28. **The frequency criteria read gnomAD's exomes in the benchmark, and the joint release adds the genomes at a measured price and without the region flags.**
    The benchmark and the setup guide build the `gnomad` database from the exomes release.
    gnomAD's joint release of exomes and genomes can be built instead (`analysis/acmg_benchmark/scripts/sa_sources/build_gnomad_release_per_chrom.sh joint`), and the builder detects its `_joint` field names, but in the headers read on 4 November 2026 the release declares no `lcr`, `segdup` or `non_par` flag.
    A database built from it therefore cannot withhold BA1, BS1, BS2 and PM2 in a low-complexity region or a segmental duplication (`gnomad_region_flags_block_frequency`), and cannot say which X-chromosome variants lie outside the pseudoautosomal regions.
    The builder derives `non_par` for GRCh38 from the pseudoautosomal coordinates (chrX 10,001 to 2,781,479 and 155,701,383 to 156,030,895; chrY 10,001 to 2,781,479 and 56,887,903 to 57,217,415), which agree with the flag in gnomAD's genomes file at the edges that could be read (exactly at the end of the first region of the X chromosome, and within the gaps between the records at the others), and `sa-build` prints a warning to standard error that the other two flags are missing.
    A joint database built before v38 has no `nonPar` flag in it, and an absent flag reads as "cannot say", so BS2 counts no hemizygotes from it: rebuild it with the current builder.
    On every chromosome (673,278 benchmark variants) the joint release holds a record for every variant that the exomes hold, with an allele count, an allele number and a count of homozygotes that are never lower, and a record for 17,382 more.
    A database from it changes 13,302 calls (12,983 on ClinVar benign-side variants, 209 on VUS and 110 on pathogenic-side variants) and gives 10,190 more exact calls and 9,321 more in direction than the exomes, for one more reference-discordant benign call (24 to 25) and two more reference-discordant pathogenic calls (24 to 26).
    4,075 and 3,473 of those are what the exomes give when the region flags withhold nothing, with one more reference-discordant pathogenic call (24 to 25), so the genomes themselves add 6,115 and 5,848, with one more reference-discordant benign call and one more reference-discordant pathogenic call.
    A release without the flags gains agreement by giving up the safeguard, and the three added errors and the one removed are variants at a threshold: KMT2C p.Cys391Ter and RAI1 p.Gln280HisfsTer83 become Likely pathogenic because their filtering allele frequency falls from 0.522% and 0.537% to 0.495% and 0.498%, below BS1's threshold of 0.5%, so the site-level frequency veto no longer holds them at VUS; TXNL4A's deletion upstream of the gene (ClinVar Pathogenic) becomes Likely benign because BS1 fires at 0.560%, where the exomes, which cover it with 364,572 alleles against the joint release's 516,266, give 0.384%; NF1 p.Met895Ile becomes Likely pathogenic because the segmental-duplication flag no longer withholds PM2; and PIEZO1's splice-acceptor deletion (ClinVar Benign) stops being Likely pathogenic, at 0.947% against 0.198%.
    Thirteen ClinVar pathogenic-side variants lose a pathogenic-side call and 39 gain one.
    Ten of the 13 lose PM2 because their record carries a filter of gnomAD's quality control (AC0 for eight, AS_VQSR for one and both for one): five have a clean record in the exomes, with one to four alleles, and the filter is the genomes', which the database reports because it reports a filter that fired in either callset, and five are not in the exomes at all; three gain BS2 because the low-complexity flag of the exomes had withheld it (MYO6 p.Gln918ThrfsTer24 and PRRT2 p.Arg217ProfsTer8 and p.Arg217GlufsTer12).
    The either-callset rule is deliberate, as the conservative direction for a guard against trusting a frequency, and where the filter is AC0 and the combined count is positive, 741 of the benchmark's variants lose PM2 this way against the exomes without the flags, and 11 calls move, 9 of them on ClinVar pathogenic-side variants.
    On this benchmark the joint release is the better-measured source for classification, with about 3,100 calls in direction for each added opposite-direction call where the manuscript's knee rule (Section 2.5) is at about a hundred, and what it gives up is the flags, whose protection a benchmark of ClinVar's classified variants cannot measure.
    The benchmark of the manuscript still reads the exomes, because every other number in it was measured on them; `ACMG_SETUP.md` says which to build.
    A second gnomAD database in the annotation directory used to answer in place of the first without a message, the later path winning; both annotation loops now record a source that two databases answer, and the command line says so at the end of the run.
    The extraction of the benchmark's 24,350 regions from gnomAD's public files took about 6 hours on one machine, limited by the network (the last 19 chromosomes took 4 hours 46 minutes at about 10 MB/s), and the 24 databases total 610 MB, against 369 MB for the exomes (`manuscript/acmg`, Section 3.8).

29. **BS2's arm for homozygotes and hemizygotes is not set aside beside BS1 and BA1, and on the X chromosome it counts the same frequency twice.**
    The rule that sets BS2's carrier count aside beside BS1 and BA1 does not reach the arm for homozygotes and hemizygotes.
    A homozygote is a genotype that an allele frequency does not give, but the hemizygotes of the non-pseudoautosomal X are part of the allele count that BS1 and BA1 read.
    On the benchmark BS2 is met beside BS1 or BA1 through that arm on 36,818 variants, 1,598 of them on the X chromosome.
    BA1 alone is Benign, so the doubling changes a call only beside BS1 without BA1: 1,301 variants (443 on the X chromosome), all called Benign, of which 1,298 are ClinVar benign-side variants and 3 are ClinVar VUS variants, and none is a ClinVar pathogenic-side variant.
    It was left as it was because the benchmark shows no harm from it, and setting the arm aside would only lower the call of variants that ClinVar also calls benign-side; an X-linked gene with a pathogenic variant at such a frequency would be the case that decides it (`manuscript/acmg`, Section 3.8).

30. **PM2 reads a variant with no gnomAD record as absent, whether or not gnomAD could have seen the site.**
    The classifier reads no coverage: a database holds the records of the variants that gnomAD called, and a variant that is not among them has no record and no allele number to say whether its site was covered.
A run with no gnomAD database loaded at all is not this case: PM2 is not evaluated for it, because no variant has a record there to be missing.
    With `pm2_absent_when_no_record = true`, the default, PM2 is met for such a variant unless ClinVar's own population frequencies show that it is not extremely rare, and that includes sites that the exomes do not cover.
    The reviewer of the benchmark gave a site covered in under half of gnomAD's exomes (PRODH) as a reason for one of the PVS1 calls on a ClinVar benign-side variant.
    Setting it to `false` leaves PM2 not evaluated wherever there is no record, which on the benchmark costs 0.22 points of agreement and 2.9 points of pathogenic-side recall (`manuscript/acmg`, Section 3.4 and Figure 3A).
    The joint release does not change the rule: a variant that is in neither callset is still read as absent, and nothing says whether its site was covered.
    A PM2 that asked for coverage would need gnomAD's coverage tracks as a source, which is not built (`manuscript/acmg`, Section 4.4).

## The frequency dead zone between PM2 and BS1

Both frequency thresholds are set by measurement, and between them is a range where fastVEP
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
not also cost reference-discordant benign calls elsewhere - the sweep in
[Choosing the BS2 prevalence threshold](#choosing-the-bs2-prevalence-threshold) shows the same behaviour
for its own threshold.

**Why it is not closed by moving a number.** Narrowing the band means lowering BS1, and the
BS1 sweep recorded in [`bs1_af_threshold`](../crates/fastvep-classification/src/config.rs)
shows the cost climbing steeply below 0.5 %: past 0.05 % the reference-discordant benign count roughly
triples. A single global value cannot be right for both a dominant early-onset condition and
a common recessive one, which is what Richards 2015's actual wording - "greater than
expected **for disorder**" - has been asking for all along.

**What would close it** is a per-gene-disease attribute table: prevalence, penetrance class,
onset class and any published VCEP BA1/BS1 threshold, keyed by gene and disorder. Whiffin
2017 gives the formula. Most of that table now exists - see
[Per-gene frequency thresholds](#per-gene-frequency-thresholds-clingen-vcep-specifications) - and what
remains missing is penetrance, which nothing publishes in machine-readable form.

## Per-gene frequency thresholds (ClinGen VCEP specifications)

Where a ClinGen Variant Curation Expert Panel has published a BA1, BS1 or PM2 threshold for a gene,
it outranks fastVEP's global default. The panel has looked up that disorder's prevalence,
penetrance and allelic heterogeneity; a number measured across all genes cannot.

`gene_overrides` carries them, one per gene:

```toml
[gene_overrides.ABCA4]
ba1_af_threshold = 0.163      # the ABCA4 VCEP's published threshold, 3x looser than the default
bs1_af_threshold = 0.0163
pm2_af_threshold = 0.0001
```

A generated table of 304 thresholds across 117 genes ships at
[`analysis/acmg_benchmark/data/vcep_thresholds.toml`](../analysis/acmg_benchmark/data/vcep_thresholds.toml),
read out of the [ClinGen CSpec Registry](https://cspec.genome.network/cspec/ui/svi/) by
[`build_vcep_thresholds.py`](../analysis/acmg_benchmark/scripts/build_vcep_thresholds.py).
Load it with `--acmg-config`. Alongside it is an audit table carrying the verbatim sentence
every number was read out of, and every rejection with its reason.

**It is opt-in, and that is deliberate.** On the ClinVar 2-star+ benchmark it buys 2,028 more
correct Benign calls and lifts exact match by 0.26 pp, at a cost of four new reference-discordant benign
calls. Two of those four are founder variants - ANO5 and CAPN3, both limb-girdle muscular
dystrophy - that the panels explicitly excluded, in prose the generator cannot read. 32 of
the 286 usable thresholds carry a caveat of that shape. `ba1_exceptions` is the right home for
them and is currently empty of panel-specific entries, so a tightened threshold today is applied
without the exclusions that were published with it. Filling that list is what has to come
before these thresholds become a default.

The wider attribute table -
[`gene_disease_attributes.tsv`](../analysis/acmg_benchmark/data/gene_disease_attributes.tsv),
6,857 gene-disorder pairs - carries inheritance (96 % coverage), earliest age of onset
(91 %) and prevalence class (85 %) from Orphanet. **Penetrance is 0 % and is emitted blank**,
because no public resource publishes it per gene-disorder pair in machine-readable form. It
is the field BS2's "full penetrance expected at an early age" precondition most wants, and
closing it needs a clinician rather than a download.

## Curated frequency-exception list

A per-gene threshold answers "how common is too common in this gene".
It cannot answer "this particular allele's frequency is not evidence about it at all", and
two classes of allele need exactly that answer.

**Common and pathogenic.** Ghosh 2018 (PMID 30311383) lists nine variants that sit above 5 %
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

**Common, and a medical geneticist says it is not evidence.**
Three more ship with `blocks = ["BA1", "BS1", "BS2"]` because the round-6 reviewer named them.
Each entry quotes the reviewer and gives the frequencies read from the gnomAD v4 database the benchmark runs against.
No PMID came with any of them, so none is cited, and no mechanism is asserted that the review did not state.

| Variant | gnomAD AF | Population maximum | Homozygotes | Reviewer |
|---|---:|---:|---:|---|
| G6PD `c.563C>T` (p.Ser188Phe) | 1.5e-3 | 4.4e-2 | 21 | "Mediterranean deficiency variant; G6PD variants are common;exceptions to Benign population criteria" |
| G6PD `c.202G>A` (p.Val68Met) | 3.7e-3 | 0.131 | 173 | "G6PD variants are common;exceptions to Benign population criteria" |
| SERPINC1 `c.439A>G` (p.Thr147Ala) | 1.8e-4 | 7.0e-3 | 3 | "founder variant in the African population; exceptions to Benign population criteria" |

The two G6PD alleles are not a ruling on the gene.
"G6PD variants are common" suggests that every G6PD variant is exempt, and only these two were put forward.
Either route to an entry, the coordinate or the gene and HGVS token, is enough for a match, so with both supplied a typo in one is hidden by the other, and an entry could ship half inert.
The benchmark that would show it is not part of CI, so each route is exercised on its own, in the form the pipeline produces, with a BA1 case for G6PD `c.202G>A`, whose population maximum crosses 5 %.

None of the first three reaches 5 %, which is why a BA1-only list never touched them, and none
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
- ClinGen VCEP Review Committee. Communication to ClinGen VCEPs regarding gnomAD v4, version 1.0 (8 March 2024); versions of August 2024 and June 2025 keep the recommendation to use the grpmax filtering allele frequency.
- ClinGen. Guidance: Classifying Variants in Genes Associated with Multiple Disorders, version 1 (January 2024).
- ClinGen Sequence Variant Interpretation (SVI) recommendations: https://clinicalgenome.org/tools/clingen-variant-classification-guidance/
