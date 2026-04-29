# ACMG Setup Guide — Building Supplementary Annotation Databases

This guide walks through downloading and building every supplementary annotation source the ACMG-AMP classifier needs. Run-of-the-mill annotation needs only a subset; full ACMG classification needs more.

If you just want to *use* the classifier, see [ACMG.md](ACMG.md) for the criteria reference and configuration. This document is about getting the data in place.

> **Critical concept:** `fastvep sa-build` is a *converter*, not a downloader. It reads a source VCF/TSV/wig file you already downloaded and produces a `.osa` (data) + `.osa.idx` (index) pair. **Building the database does not download the source.** If you skip the download, `sa-build` may succeed on an empty/placeholder file and produce a `.osa` that contains zero records, which then yields blank annotations at runtime. (See [Issue #4](https://github.com/Huang-lab/fastVEP/issues/4).)

## Quick verification

After building any source, run this one-liner to confirm the database actually has data:

```bash
ls -la sa_databases/<source>.osa sa_databases/<source>.osa.idx
# Both files should be present. .osa should be many MB to many GB depending on source.
# A clinvar.osa under ~5 MB usually means an empty or partial build.
```

You can also test annotation with one variant from the source VCF and grep for the expected SA key:

```bash
echo -e '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t69134\trs1\tA\tG\t.\t.\t.' > /tmp/probe.vcf
fastvep annotate -i /tmp/probe.vcf --gff3 data/Homo_sapiens.GRCh38.115.gff3 \
  --sa-dir sa_databases/ --output-format json | grep -i 'clinvar\|gnomad'
```

## Minimum viable ACMG setup

Roughly 90% of ACMG criteria fire with just three sources:

| Source | Used by | Download size | Build time |
|---|---|---|---|
| **gnomAD v4** | BA1, BS1, BS2, PM2 | ~30–60 GB (genomes), ~5 GB (exomes) | ~20–60 min |
| **ClinVar** | PS4, PP5, BP6 | ~50 MB | ~30 sec |
| **REVEL** | PP3, BP4 (missense) | ~3 GB | ~5 min |

That leaves the splicing path (PP3/BP4 splice, BP7 — needs SpliceAI), conservation-only paths (PP3/BP4/BP7 fallback — needs PhyloP), and gene-level criteria (PVS1, PS1, PM1, PM5 — see [Gene-level sources](#gene-level-sources-oga) below).

## Source-by-source recipes

Every recipe assumes you start from a fresh `sa_databases/` directory:

```bash
mkdir -p sa_databases && cd sa_databases
```

### ClinVar — clinical significance

Used by PS4, PP5, BP6. Updated weekly.

```bash
# 1. Download
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz

# 2. Verify the file is real (~30–60 MB, contains actual VCF records)
zcat clinvar.vcf.gz | head -30 | grep -c '^#CHROM'  # → 1
zcat clinvar.vcf.gz | grep -v '^#' | wc -l           # → ~3M records

# 3. Build
fastvep sa-build --source clinvar -i clinvar.vcf.gz -o clinvar --assembly GRCh38

# Expected output:
#   sa_databases/clinvar.osa      (~80–120 MB)
#   sa_databases/clinvar.osa.idx  (~5–10 MB)
```

### gnomAD — population allele frequencies

Used by BA1, BS1, BS2, PM2. The single largest source by far.

gnomAD v4 is split per chromosome. You can build per-chromosome `.osa` files and put them all in `--sa-dir`, or merge first and build once. Per-chromosome is recommended for incremental builds.

```bash
# 1. Download (per-chromosome, ~30 GB total for genomes v4.0)
# Browse https://gnomad.broadinstitute.org/downloads and pick:
#   gnomad.genomes.v4.0.sites.chr{1..22,X,Y}.vcf.bgz
# Or use gsutil:
gsutil -m cp gs://gcp-public-data--gnomad/release/4.0/vcf/genomes/gnomad.genomes.v4.0.sites.chr*.vcf.bgz .

# 2. Build per chromosome
for chr in {1..22} X Y; do
  fastvep sa-build --source gnomad \
    -i gnomad.genomes.v4.0.sites.chr${chr}.vcf.bgz \
    -o gnomad_chr${chr} \
    --assembly GRCh38
done

# Expected output:
#   sa_databases/gnomad_chr*.osa      (totals ~5–8 GB)
#   sa_databases/gnomad_chr*.osa.idx
```

> **Note on naming:** all `.osa` files in `--sa-dir` are loaded; their JSON keys come from the source type (set when building), not the filename. So `gnomad_chr1.osa` and `gnomad_chr2.osa` both register under the `gnomad` key.

### REVEL — missense pathogenicity

Used by PP3 and BP4 (missense path, ClinGen SVI calibrated).

```bash
# 1. Download (one-time, ~3 GB; REVEL has no version updates after 2023)
wget https://rothsj06.dmz.hpc.mssm.edu/revel-v1.3_all_chromosomes.zip
unzip revel-v1.3_all_chromosomes.zip
# Produces: revel_with_transcript_ids (one TSV per chromosome, no header reuse)

# 2. Build
fastvep sa-build --source revel -i revel_with_transcript_ids -o revel --assembly GRCh38

# Expected output:
#   sa_databases/revel.osa      (~2 GB)
#   sa_databases/revel.osa.idx
```

### SpliceAI — splice site predictions

Used by PP3 and BP4 (splice path) and BP7. Requires Illumina account for download.

```bash
# 1. Download from https://basespace.illumina.com/s/otSPW8hnhaZR
#    (look for spliceai_scores.masked.snv.hg38.vcf.gz)

# 2. Build
fastvep sa-build --source spliceai -i spliceai_scores.masked.snv.hg38.vcf.gz -o spliceai --assembly GRCh38

# Expected output:
#   sa_databases/spliceai.osa      (~10 GB)
#   sa_databases/spliceai.osa.idx
```

### dbNSFP — SIFT / PolyPhen / metaSVM

Used by PP3 and BP4 (consensus signals; in transparency-only mode after PR1, no longer drives PP3/BP4 firing).

```bash
# 1. Download from https://sites.google.com/site/jpopgen/dbNSFP
#    (latest is dbNSFP4.x; ~30 GB unpacked)

# 2. Build
fastvep sa-build --source dbnsfp -i dbNSFP4.5a.zip -o dbnsfp --assembly GRCh38
```

### PhyloP — conservation scores

Used by PP3, BP4, BP7. Position-level, not allele-level.

```bash
# 1. Download (UCSC)
wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/phyloP100way/hg38.phyloP100way.wigFix.gz

# 2. Build
fastvep sa-build --source phylop -i hg38.phyloP100way.wigFix.gz -o phylop --assembly GRCh38
```

### GERP — evolutionary rate

Used by PP3 (consensus path).

```bash
# 1. Download from https://hgdownload.soe.ucsc.edu/gbdb/hg38/bbi/
#    Convert to TSV: chrom, pos, score (one row per position)

# 2. Build
fastvep sa-build --source gerp -i gerp_scores.tsv -o gerp --assembly GRCh38
```

## Gene-level sources (`.oga`)

`fastvep sa-build` supports three gene-level sources. The output is a `.oga` file (gene-keyed annotation index); place it in the same `--sa-dir` as your `.osa` files and the runtime will pick it up automatically.

| Source | Builder key | json_key (runtime) | Used by |
|---|---|---|---|
| OMIM `genemap2.txt` | `omim` | `omim` | PVS1, BS2, PM3, BP2 |
| gnomAD gene constraints | `gnomad_genes` | `gnomad_genes` | PVS1, PP2, BP1 |
| ClinVar protein index | `clinvar_protein` | `clinvar_protein` | PS1, PM1, PM5 |

The classifier degrades gracefully when these are absent — every criterion that depends on a missing `.oga` is marked `evaluated: false` rather than firing — so it's fine to start with allele-level sources only and layer these in later.

### OMIM — gene-phenotype map

Used by PVS1 (LOF-intolerance proxy from disease association), BS2 (inheritance pattern), PM3 (recessive disorder check), BP2 (recessive cis check).

```bash
# 1. Download genemap2.txt from OMIM (account required, free for academic use)
#    https://www.omim.org/downloads
#    File: genemap2.txt (~5 MB, tab-separated)

# 2. Build
fastvep sa-build --source omim -i genemap2.txt -o sa_databases/omim --assembly GRCh38

# Expected output:
#   sa_databases/omim.oga (~2–4 MB)
```

### gnomAD gene constraints

Used by PVS1 (pLI / LOEUF), PP2 (mis_z), BP1 (LOF-only gene check).

```bash
# 1. Download the constraint TSV
#    https://gnomad.broadinstitute.org/downloads → "Gene constraint metrics"
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv

# 2. Build
fastvep sa-build --source gnomad_genes -i gnomad.v4.1.constraint_metrics.tsv -o sa_databases/gnomad_genes --assembly GRCh38

# Expected output:
#   sa_databases/gnomad_genes.oga (~4–6 MB)
```

### ClinVar protein index

Used by PS1 (same-AA pathogenic match), PM5 (different-AA at same position), PM1 (hotspot density).

This source is built from the **same** ClinVar VCF you used for the allele-level `clinvar.osa` — but `clinvar_protein` extracts a different view (pathogenic missense indexed by protein position) and writes a separate `.oga`.

```bash
# Reuse the clinvar.vcf.gz you already downloaded for the .osa build
fastvep sa-build --source clinvar_protein -i clinvar.vcf.gz -o sa_databases/clinvar_protein --assembly GRCh38

# Expected output:
#   sa_databases/clinvar_protein.oga (~5–10 MB)
```

### Verifying gene-level annotations

`.oga` files don't show up in the standard `ls` size sanity check the same way `.osa` files do — they're more compact (single MB range). To verify:

```bash
# Each .oga is loaded at startup; --acmg run will log the gene count.
fastvep annotate -i tests/test.vcf --gff3 tests/test.gff3 \
  --sa-dir sa_databases/ --acmg --output-format json 2>&1 \
  | grep -E 'Loaded gene annotations'
# Expected:
#   Loaded gene annotations: OMIM (... N genes)
#   Loaded gene annotations: gnomAD gene constraints (... N genes)
#   Loaded gene annotations: ClinVar protein index (... N genes)
```

## Putting it together

```bash
fastvep annotate \
  -i your_variants.vcf \
  -o annotated.vcf \
  --gff3 data/Homo_sapiens.GRCh38.115.gff3 \
  --fasta data/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  --sa-dir sa_databases/ \
  --acmg \
  --output-format json \
  --hgvs
```

`--output-format json` is recommended while you're verifying the setup — it shows the full ACMG result block with every criterion's evaluation, so you can confirm the SA data is being picked up:

```json
{
  "acmg": {
    "classification": "LikelyPathogenic",
    "shorthand": "LP",
    "triggered_rule": "PVS + >=1 PP (SVI)",
    "criteria": [
      { "code": "PM2_Supporting", "met": true, "summary": "Absent in gnomAD ..." },
      ...
    ]
  }
}
```

If `criteria` shows everything as `met: false` or `evaluated: false`, double-check:

1. **The `.osa` files are non-trivially sized** (gnomAD genomes should be several GB, not a few KB)
2. **Contig naming matches** between your input VCF and the source VCF used to build the database (`chr1` vs `1` is the most common mismatch)
3. **Assembly matches** (GRCh38 throughout — mixing GRCh37 sources will silently produce no matches)
4. **You passed `--acmg`** — without it, the classifier doesn't run

## Configuration

Per-criterion thresholds, gene-specific overrides, and trio settings all live in a single TOML file. See [ACMG.md § Configuration File](ACMG.md#configuration-file) for the complete schema.

```bash
fastvep annotate ... --acmg --acmg-config my_thresholds.toml
```
