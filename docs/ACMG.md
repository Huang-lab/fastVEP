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
pm2_ad_af_threshold = 0.0      # PM2: AD/unknown — strict absence (AC=0 AND AF=0)
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

# ── ClinGen SVI behavior ──
pm2_downgrade_to_supporting = true     # Downgrade PM2 from Moderate to Supporting (SVI)
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

# ── Gene-specific overrides ──
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

## Required Data Sources

ACMG classification draws on multiple supplementary annotation (SA) sources. Place `.osa`/`.osa2` (allele-level) and `.oga` (gene-level) files in the SA directory:

### Allele-Level Sources (`.osa` / `.osa2`)

| Source | SA Key | Used By | Description |
|---|---|---|---|
| **gnomAD** | `gnomad` | BA1, BS1, BS2, PM2 | Population allele frequencies + homozygote counts |
| **ClinVar** | `clinvar` | PS4, PP5, BP6 | Clinical significance, review status, phenotypes |
| **REVEL** | `revel` | PP3, BP4 | Missense pathogenicity score (0-1) |
| **SpliceAI** | `spliceai` | PP3, BP7 | Splice site delta scores |
| **dbNSFP** | `dbnsfp` | PP3, BP4 | SIFT/PolyPhen predictions |
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
| **OMIM** | `omim` | PVS1, BS2, PM3, BP2 | Disease associations, inheritance patterns |
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

When a `.oga` is missing, dependent criteria (PVS1, PS1, PM1, PM5, PM3, BP1, BP2, PP2, BS2) degrade gracefully to `evaluated: false` rather than misfiring. See [ACMG_SETUP.md](ACMG_SETUP.md) for download URLs, expected file sizes, and end-to-end verification.

## Evidence Criteria Reference

### Pathogenic Criteria

| Code | Strength | Description | Data Source | Automatable |
|---|---|---|---|---|
| **PVS1** | Very Strong | Null variant in LOF-intolerant gene | Consequence + pLI/LOEUF/OMIM | Yes |
| **PS1** | Strong | Same AA change as established pathogenic | ClinVar protein index | Yes (with .oga) |
| **PS2** | Strong | De novo with confirmed parents | Trio VCF genotypes | Yes (with trio) |
| **PS3** | Strong | Functional studies show damaging | External data | No |
| **PS4** | Strong | Prevalence in affected >> controls | ClinVar expert panel | Partial |
| **PM1** | Moderate | Mutational hotspot / critical domain | ClinVar protein density | Yes (with .oga) |
| **PM2** | Supporting* | Absent/rare in population databases | gnomAD AF | Yes |
| **PM3** | Moderate | In trans with pathogenic (recessive) | Phased VCF + ClinVar | Yes (with trio) |
| **PM4** | Moderate | Protein length change (in-frame/stop-loss) | Consequence type | Yes |
| **PM5** | Moderate | Different pathogenic missense at same residue | ClinVar protein index | Yes (with .oga) |
| **PM6** | Moderate | Assumed de novo (partial confirmation) | Partial trio VCF | Yes (with trio) |
| **PP1** | Supporting | Co-segregation in family | Pedigree data | No |
| **PP2** | Supporting | Missense in constrained gene | Gene misZ score | Yes |
| **PP3** | Supporting-Strong | Computational evidence (deleterious) | REVEL (missense, Pejaver 2022 calibrated) + SpliceAI (caps at Supporting per Walker 2023) | Yes |
| **PP4** | Supporting | Phenotype-specific for single-gene disease | HPO phenotype data | No |
| **PP5** | Supporting | Reputable source reports pathogenic | ClinVar (disabled by default per SVI) | Partial |

*\*PM2 is downgraded from Moderate to Supporting per ClinGen SVI recommendation.*

### Benign Criteria

| Code | Strength | Description | Data Source | Automatable |
|---|---|---|---|---|
| **BA1** | Standalone | AF > 5% in any population | gnomAD population AFs | Yes |
| **BS1** | Strong | AF > expected for disorder | gnomAD global AF | Yes |
| **BS2** | Strong | Observed in healthy adult | gnomAD hom count + OMIM inheritance | Yes |
| **BS3** | Strong | Functional studies show no damage | External data | No |
| **BS4** | Strong | Lack of segregation | Pedigree data | No |
| **BP1** | Supporting | Missense in LOF-only gene | pLI + misZ | Yes |
| **BP2** | Supporting | In trans/cis with pathogenic | Phased VCF + ClinVar | Yes (with trio) |
| **BP3** | Supporting | In-frame indel in repeat region | Consequence + RepeatMasker | Yes (with .osi) |
| **BP4** | Supporting-Very Strong | Computational evidence (benign) | REVEL (missense, Pejaver 2022 calibrated; ≤ 0.003 reaches Very Strong) + SpliceAI ≤ 0.1 (Walker 2023) | Yes |
| **BP5** | Supporting | Alternate molecular basis in case | Case-level analysis | No |
| **BP6** | Supporting | Reputable source reports benign | ClinVar (disabled by default per SVI) | Partial |
| **BP7** | Supporting | Synonymous + no splice impact + not conserved | Consequence + SpliceAI + PhyloP | Yes |

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
- **PM3** (Moderate): In a recessive-inheritance gene, proband is heterozygous for two variants, and the companion variant is ClinVar pathogenic. Phase-aware: uses phased genotypes (VCF `|` separator) when available.
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
  combiner.rs         # 18 classification combination rules
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
5. **PS1/PM5/PM1** require the ClinVar protein index `.oga` file to be built and loaded
6. **BP3** requires RepeatMasker interval `.osi` file to be built and loaded
7. **PS2/PM6/PM3/BP2** require a multi-sample VCF with trio sample names configured
8. Compound heterozygote detection (PM3/BP2) works per-batch in the CLI; variants in different batches within the same gene may not be cross-referenced

## References

- Richards S, et al. Standards and guidelines for the interpretation of sequence variants. *Genet Med*. 2015;17(5):405-424.
- Tavtigian SV, et al. Modeling the ACMG/AMP variant classification guidelines as a Bayesian classification framework. *Genet Med*. 2018;20(9):1054-1060.
- ClinGen Sequence Variant Interpretation (SVI) recommendations: https://clinicalgenome.org/tools/clingen-variant-classification-guidance/
- Pejaver V, et al. Calibration of computational tools for missense variant pathogenicity classification and ClinGen recommendations for PP3/BP4 criteria. *Am J Hum Genet*. 2022;109(12):2163-2177.
