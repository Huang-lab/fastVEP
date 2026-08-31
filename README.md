# fastVEP

A high-performance Variant Effect Predictor written in Rust. fastVEP predicts the functional consequences of genomic variants (SNPs, insertions, deletions, structural variants) on genes, transcripts, and protein sequences, with direct integration of clinical and population databases.

fastVEP is inspired by and aims to be compatible with [Ensembl VEP](https://www.ensembl.org/info/docs/tools/vep/index.html) and [Illumina Nirvana](https://github.com/Illumina/Nirvana), while delivering significantly better performance through Rust's zero-cost abstractions and native parallelism.

**Try it now:** A hosted web server is available at [fastVEP.org](https://fastVEP.org) — paste VCF data and get annotated results instantly, no installation required. As of July 2026 it has served **2,372 genome analysis sessions** and annotated **19,013 variants**.

The hosted instance is for interactive use in the browser. For programmatic access, run your own server: see the **[REST API guide](docs/API.md)**.

## Features

- **Variant Consequence Prediction** — Classifies variants using 49 [Sequence Ontology](http://www.sequenceontology.org/) terms (missense, frameshift, splice donor, copy_number_change, transcript_ablation, etc.)
- **Structural Variant Support** — Full SV pipeline: `<DEL>`, `<DUP>`, `<INV>`, `<CNV>`, `<BND>`, `<INS>`, `<STR>` with SV-specific consequence prediction
- **Supplementary Annotations** — Direct integration with ClinVar, gnomAD, dbSNP, COSMIC, 1000 Genomes, TOPMed, MitoMap via the native fastSA format (v1: zstd block compression; v2: echtvar-inspired chunked ZIP with Var32 encoding, parallel u32 value arrays, delta encoding). Both use a shared, thread-scaled, byte-budgeted LRU cache of decompressed blocks/chunks and lock-free reads, so a dense whole-genome source (e.g. SpliceAI, dbSNP) queried in parallel stays fast without re-decompressing the same data
- **Prediction Scores** — PhyloP, GERP, REVEL, SpliceAI, PrimateAI, DANN, AlphaMissense conservation and pathogenicity scores; SIFT/PolyPhen via dbNSFP
- **Gene-Level Annotations** — OMIM phenotypes, gnomAD gene constraint (pLI, LOEUF), ClinGen gene-disease validity
- **Filter Engine** — Expression-based filtering compatible with VEP's filter_vep syntax
- **HGVS Nomenclature** — Generates HGVSg, HGVSc, and HGVSp notations with 3' normalization
- **Multiple Output Formats** — VCF (with 49-field CSQ), tab-delimited, JSON (including Nirvana-style structured output)
- **Multi-Sample Support** — Parse FORMAT/GT/DP/GQ/AD fields per sample with genotype classification
- **Regulatory Region Detection** — Promoters, enhancers, CTCF binding sites, TF binding sites from Ensembl regulatory build
- **Mitochondrial Support** — Circular coordinate handling, vertebrate mitochondrial codon table (NCBI table 2)
- **Custom Annotations** — User-provided VCF (`--source custom_vcf`) and BED (`--source custom_bed`) files; `.osi` interval databases load alongside `.osa` via `--sa-dir`
- **ACMG-AMP Classification** — `--acmg` runs the full Richards 2015 + ClinGen SVI rule set (28 criteria, configurable thresholds, trio / compound-het support via `--proband`/`--mother`/`--father`)
- **VEP `--merged` Cache** — `--gff3` is repeatable on `annotate` and `cache`; combine Ensembl + RefSeq in a single run with per-transcript SOURCE labels
- **`--sa-only` Mode** — Skip the default CSQ pipeline and emit only supplementary annotations, useful for re-annotating already-annotated VCFs
- **Gzipped VCF Input** — `annotate` auto-detects `.vcf.gz` / `.vcf.bgz` (no upstream decompression needed)
- **Web Interface** — Built-in web GUI for interactive variant annotation
- **GFF3 Annotation Support** — Load gene models from standard GFF3 files (any organism)

## Quick Start

### 1. Install Rust (if you don't have it)

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
source "$HOME/.cargo/env"
```

### 2. Build and install fastVEP

```bash
git clone https://github.com/Huang-lab/fastVEP.git
cd fastVEP

# Build and install both binaries to ~/.cargo/bin/
cargo install --path crates/fastvep-cli   # fastvep (CLI annotator)
cargo install --path crates/fastvep-web   # fastvep-web (production web server)

# Verify it works
fastvep --version
```

> **Note:** `cargo install` places the binary in `~/.cargo/bin/`. If `fastvep` is not found after install, run `source "$HOME/.cargo/env"` or add this line to your `~/.zshrc` (or `~/.bashrc`):
> ```bash
> source "$HOME/.cargo/env"
> ```

#### Alternative: build a conda package

> fastVEP is **not yet published on bioconda** — `conda install -c bioconda fastvep` will not work. The recipe below is bioconda-shaped and builds locally in the meantime.

The repo ships a recipe under `conda/recipe/` that builds both `fastvep` and `fastvep-web` into a local conda package (Linux and macOS):

```bash
# One-time: tools for building conda packages
conda install -n base -c conda-forge conda-build cargo-bundle-licenses

# Build the package from the repo root (builds the tagged release tarball)
conda build conda/recipe

# Install into a fresh environment
conda create -n fastvep -c local fastvep
conda activate fastvep
fastvep --version
```

To build your unreleased working tree instead, swap the `url:`/`sha256:` pair in
[`conda/recipe/meta.yaml`](conda/recipe/meta.yaml) for the commented-out `path: ../..` line.

### 3. Try it — annotate the included test data

fastVEP ships with a small test VCF and GFF3 so you can try it immediately:

```bash
# Annotate 12 test variants covering SNVs, indels, splice sites, UTRs, and intergenic regions
fastvep annotate -i tests/test.vcf --gff3 tests/test.gff3 --hgvs --output-format tab
```

### 4. Build supplementary annotation databases

```bash
# Every allele-level source builds the smaller, faster v2 `.osa2` format
# automatically (`--format auto` is the default); pass `--format osa` to force
# v1 `.osa`. Already have `.osa` files? `fastvep sa-convert` upgrades them in
# place of a re-download.

# Build ClinVar annotation database (writes clinvar.osa2)
fastvep sa-build --source clinvar --input clinvar.vcf.gz --output clinvar

# Build gnomAD population frequency database (writes gnomad.osa2)
fastvep sa-build --source gnomad --input gnomad.genomes.v4.vcf.bgz --output gnomad

# Build AlphaMissense pathogenicity predictions (writes alphamissense.osa2).
# Download AlphaMissense_hg38.tsv.gz from Zenodo record 8208688.
fastvep sa-build --source alphamissense --input AlphaMissense_hg38.tsv.gz --output alphamissense

# Build PhyloP conservation scores (see docs/ACMG_SETUP.md for how to
# obtain hg38.phyloP100way.wigFix.gz — UCSC ships it as one file per
# chromosome, not a single combined download)
fastvep sa-build --source phylop --input hg38.phyloP100way.wigFix.gz --output phylop

# Build SpliceAI predictions (writes spliceai.osa2). The densest source
# fastVEP ships against, and the one v2 helps most: its eight delta
# scores/positions become numeric columns and the gene symbol a string-table
# index, for 0.14x-0.56x the v1 size depending on how much the scores vary.
fastvep sa-build --source spliceai --input spliceai_scores.vcf.gz --output spliceai

# Upgrade a v1 database built before .osa2 became the default, without
# re-downloading the upstream release (writes dbsnp.osa2 alongside it)
fastvep sa-convert --input sa_databases/dbsnp.osa
```

### 5. Annotate with supplementary databases

```bash
# Annotate with all databases in a directory
fastvep annotate \
  -i your_variants.vcf \
  -o annotated.vcf \
  --gff3 Homo_sapiens.GRCh38.112.gff3 \
  --fasta Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  --sa-dir /path/to/annotation_databases/ \
  --hgvs
```

### 6. Filter annotated variants

```bash
# Filter for high-impact or rare missense variants
fastvep filter \
  -i annotated.vcf \
  --filter "IMPACT is HIGH or (Consequence in missense_variant and AF < 0.001)"
```

### 7. Launch the web interface

```bash
# Quick start — uses a built-in example gene model (OR4F5, chr1)
fastvep-web

# With your own data
fastvep-web --gff3 Homo_sapiens.GRCh38.115.gff3 --fasta Homo_sapiens.GRCh38.dna.primary_assembly.fa

# With supplementary annotations (ClinVar, gnomAD, etc.)
fastvep-web --gff3 genes.gff3 --fasta ref.fa --sa-dir /path/to/sa_databases/
```

Open http://localhost:8080 in your browser. The web interface lets you paste VCF data, switch gene models, and view results in an interactive table.

> **Note:** `fastvep-web` is a separate production-quality binary (axum/tokio, async, multi-connection). The legacy `fastvep web` command still works but is single-threaded.

### REST API (self-hosted)

The same binary serves a JSON API over the same annotation engine, so anything the web interface can do a script can do:

```bash
curl -s -X POST http://localhost:8080/api/annotate \
  -H 'Content-Type: application/json' \
  -d '{"vcf": "17\t43124027\t.\tG\tA\t50\tPASS\t.", "acmg": true}'
```

It is designed for self-hosting on a workstation or local network, which is also the only supported way to use it: please do not point scripts at fastVEP.org, which is a single small machine sized for interactive browser use. Bulk VCFs belong in `fastvep annotate`, not in an HTTP request.

Full setup, endpoint reference, and response shapes: **[docs/API.md](docs/API.md)**.

## Local Setup Guide

This section walks through setting up fastVEP with full annotation capabilities (gene models, reference sequence, and supplementary databases like ClinVar and gnomAD).

### Step 1: Download reference data

```bash
mkdir -p data && cd data

# Gene models (GFF3) — pick your organism
# Human GRCh38
wget https://ftp.ensembl.org/pub/release-115/gff3/homo_sapiens/Homo_sapiens.GRCh38.115.gff3.gz
gunzip Homo_sapiens.GRCh38.115.gff3.gz

# Reference FASTA (needed for HGVS and sequence context)
wget https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# Create FASTA index (enables memory-mapped access — important for large genomes)
samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa
```

### Step 2: Build supplementary annotation databases

Each supplementary database (ClinVar, gnomAD, etc.) is built in **two steps** — *download the source file*, then run `fastvep sa-build` to convert it into a fastSA database (a `.osa2` file by default; a `.osa` + `.osa.idx` pair under `--format osa`). **`sa-build` is a converter, not a downloader; if you skip the download, the resulting database will be empty and your annotations will silently come back blank.** After each build, check that the file size matches the expected magnitude (column below); a few-KB database is the tell that the source file wasn't real.

```bash
mkdir -p sa_databases

# ── ClinVar — clinical variant significance ──
# Download (~50 MB)
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
# Build (expect ~80–120 MB .osa)
fastvep sa-build --source clinvar -i clinvar.vcf.gz -o sa_databases/clinvar --assembly GRCh38

# ── gnomAD v4 — population allele frequencies ──
# Download per-chromosome from https://gnomad.broadinstitute.org/downloads
# (~30–60 GB total for genomes v4.0)
fastvep sa-build --source gnomad -i gnomad.genomes.v4.0.sites.vcf.bgz -o sa_databases/gnomad --assembly GRCh38

# ── dbSNP — variant identifiers ──
wget https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz
fastvep sa-build --source dbsnp -i GCF_000001405.40.gz -o sa_databases/dbsnp --assembly GRCh38

# ── COSMIC — somatic mutations (requires license) ──
# https://cancer.sanger.ac.uk/cosmic/download
fastvep sa-build --source cosmic -i CosmicCodingMuts.vcf.gz -o sa_databases/cosmic --assembly GRCh38
```

**Verify before moving on:**

```bash
ls -la sa_databases/*.osa
# Expected: clinvar ~100 MB; gnomad several GB; dbsnp ~5 GB.
# Anything < 1 MB usually means an empty build — re-check the source file.
```

> For ACMG-AMP classification specifically (REVEL, SpliceAI, PhyloP, dbNSFP, OMIM, ClinVar protein index, etc.), see the dedicated **[ACMG Setup Guide](docs/ACMG_SETUP.md)** — it walks through every source the classifier needs with download URLs, build commands, expected disk sizes, and a verification recipe.

### Step 3: Run the CLI annotator

```bash
fastvep annotate \
  -i your_variants.vcf \
  -o annotated.vcf \
  --gff3 data/Homo_sapiens.GRCh38.115.gff3 \
  --fasta data/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  --sa-dir sa_databases/ \
  --hgvs
```

### Step 4: Run the web server

```bash
# Install the web server binary
cargo install --path crates/fastvep-web

# Run with all annotation sources
fastvep-web \
  --gff3 data/Homo_sapiens.GRCh38.115.gff3 \
  --fasta data/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  --sa-dir sa_databases/ \
  --port 8080
```

All flags also accept environment variables (`FASTVEP_GFF3`, `FASTVEP_FASTA`, `FASTVEP_SA_DIR`, `FASTVEP_PORT`) for container deployments.

### Multi-organism setup

To serve multiple genomes from the web interface, organize data into subdirectories and use `--data-dir`. Each subdirectory is one genome — the server auto-detects GFF3, FASTA, and SA files inside.

**Directory layout:**

```
genomes/
  human_grch38/
    Homo_sapiens.GRCh38.115.gff3       # gene models (required)
    Homo_sapiens.GRCh38.dna.primary_assembly.fa   # reference (optional, for HGVS)
    Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai
    sa/                                 # supplementary annotations (optional)
      clinvar.osa2
      gnomad.osa2
      dbsnp.osa2
  mouse_grcm39/
    Mus_musculus.GRCm39.115.gff3
    mouse.fa
    mouse.fa.fai
  zebrafish/
    Danio_rerio.GRCz11.115.gff3
```

**Setup:**

```bash
mkdir -p genomes/human_grch38/sa genomes/mouse_grcm39 genomes/zebrafish

# Human: GFF3 + FASTA + SA databases
cp data/Homo_sapiens.GRCh38.115.gff3 genomes/human_grch38/
cp data/Homo_sapiens.GRCh38.dna.primary_assembly.fa* genomes/human_grch38/
# Copy every supplementary database, whichever format it built to
# (.osa2/.osa allele-level, .osi interval, .oga gene-level).
cp sa_databases/*.osa2 sa_databases/*.osa sa_databases/*.osi sa_databases/*.oga \
   genomes/human_grch38/sa/ 2>/dev/null   # ClinVar, gnomAD, PhyloP, OMIM, etc.

# Mouse
wget -O- https://ftp.ensembl.org/pub/release-115/gff3/mus_musculus/Mus_musculus.GRCm39.115.gff3.gz | gunzip > genomes/mouse_grcm39/Mus_musculus.GRCm39.115.gff3

# Zebrafish
wget -O- https://ftp.ensembl.org/pub/release-115/gff3/danio_rerio/Danio_rerio.GRCz11.115.gff3.gz | gunzip > genomes/zebrafish/Danio_rerio.GRCz11.115.gff3
```

**Run:**

```bash
fastvep-web --data-dir genomes/
```

Users can switch between organisms from the dropdown in the web UI. When a genome has a `sa/` subdirectory, its SA databases are automatically loaded. The dropdown shows "(FASTA + SA)" labels for genomes that have these resources.

`--sa-dir` is optional — if provided, it serves as a fallback for genomes that don't have their own `sa/` folder. If the directory doesn't exist, the server starts without SA (no error).

fastVEP works with any organism — just provide the matching GFF3 (and optionally FASTA for HGVS).

### Merged annotation (Ensembl + RefSeq, à la VEP `--merged`)

`--gff3` accepts multiple values, so a single run can annotate against
Ensembl and RefSeq side-by-side — fastVEP's analog of VEP's `--merged`
cache. The SOURCE column of each CSQ entry records which file produced
that transcript.

```bash
# Auto-detected SOURCE labels (filename contains "ensembl"/"gencode" → Ensembl;
# "refseq" or GCF_ prefix → RefSeq; otherwise the basename).
fastvep annotate -i variants.vcf \
  --gff3 Homo_sapiens.GRCh38.115.gff3 \
  --gff3 GCF_000001405.40.gff.gz \
  --fasta GRCh38.fa --hgvs

# Or set the labels explicitly with `LABEL=path` (also accepts comma-separated):
fastvep annotate -i variants.vcf \
  --gff3 Ensembl=ensembl.gff3,RefSeq=refseq.gff3 \
  --fasta GRCh38.fa
```

Each transcript carries its source through to all output formats:
`SOURCE` field in the VCF CSQ string, the `source` key in JSON, and the
SOURCE column in tab output (with `--canonical` / extended fields).
Transcripts from both sources are queried independently — overlap
detection, HGVS, and supplementary annotation all run per-transcript,
so RefSeq NM_… and Ensembl ENST… records appear as separate CSQ entries
for the same variant.

The auto-managed sidecar cache (`<gff3>.fastvep.cache`) only kicks in
for single-GFF3 runs. For merged workflows, pre-build a combined binary
cache with `fastvep cache --gff3 ensembl.gff3 -o combined.cache` once
per source and pass `--transcript-cache combined.cache` on subsequent
runs — though GFF3 parsing is fast enough that this is rarely needed.

#### Chromosome naming across sources (`--synonyms`)

Ensembl GFF3/FASTA use bare contig names (`17`), UCSC prefixes `chr`
(`chr17`), and NCBI RefSeq GFF3 uses accessions (`NC_000017.11`). When
you build a merged cache against one FASTA, those names must be
reconciled or sequence-building fails with
`Chromosome 'NC_000017.11' not found in FASTA index`.

- **`chr` ↔ bare and mitochondrial (`M`/`MT`/`chrM`/`chrMT`) forms resolve
  automatically** — no configuration needed.
- **RefSeq accessions need a mapping table.** Pass a VEP-style synonyms
  file with `--synonyms` (only takes effect together with `--fasta`):

  ```bash
  fastvep cache \
    --gff3 Homo_sapiens.GRCh38.115.gff3 \
    --gff3 GCF_000001405.40_GRCh38.p14_genomic.gff \
    --fasta GRCh38.fa \
    --synonyms chr_synonyms.txt \
    -o combined.cache
  ```

  The file format matches Ensembl VEP's `chr_synonyms.txt`: one line per
  contig, whitespace/tab-separated equivalent names
  (e.g. `17  chr17  NC_000017.11`). Ensembl ships one with each release;
  you can also generate a two-column map from the NCBI
  `*_assembly_report.txt` (`Sequence-Name` ↔ `RefSeq-Accn` columns).

Every transcript is **canonicalized to the FASTA's contig names at build
time**, so the merged cache uses one consistent naming scheme and
`annotate` matches your VCF regardless of which GFF3 each transcript came
from. Contigs with no FASTA match are reported in a warning and skipped
(no sequences built for them).

## Supplementary Annotation Sources

fastVEP supports direct integration with clinical and population databases through its native fastSA binary format. Build once with `fastvep sa-build`, then use `--sa-dir` to annotate:

| Source | Type | Description | Build Command |
|--------|------|-------------|---------------|
| **ClinVar** | Allele-specific | Clinical significance, review status, phenotypes | `--source clinvar` |
| **gnomAD** | Allele-specific | Population frequencies (8 populations), allele counts | `--source gnomad` |
| **dbSNP** | Allele-specific | RS IDs, global minor allele frequency | `--source dbsnp` |
| **COSMIC** | Allele-specific | Somatic mutations, gene, sample counts | `--source cosmic` |
| **1000 Genomes** | Allele-specific | Population frequencies (AFR, AMR, EAS, EUR, SAS) | `--source onekg` |
| **TOPMed** | Allele-specific | Population frequencies, allele counts | `--source topmed` |
| **MitoMap** | Allele-specific | Mitochondrial disease associations | `--source mitomap` |
| **PhyloP** | Positional | Phylogenetic conservation scores | `--source phylop` |
| **GERP** | Positional | Evolutionary rate profiling | `--source gerp` |
| **DANN** | Positional | Deleterious annotation scores | `--source dann` |
| **REVEL** | Allele-specific | Missense pathogenicity predictions | `--source revel` |
| **SpliceAI** | Allele-specific | Splice site effect predictions (delta scores) | `--source spliceai` |
| **PrimateAI** | Allele-specific | Primate-based pathogenicity | `--source primateai` |
| **dbNSFP** | Allele-specific | SIFT/PolyPhen predictions | `--source dbnsfp` |
| **AlphaMissense** | Allele-specific | Missense pathogenicity score + class | `--source alphamissense` |
| **OMIM / ClinGen GDV** | Gene-level (`.oga`) | Disease-gene annotations driving PVS1, BS2, PM3, BP2 in ACMG | `--source omim` |
| **gnomAD constraint** | Gene-level (`.oga`) | gnomAD constraint metrics (pLI, LOEUF) for PVS1, PP2, BP1 | `--source gnomad_genes` |
| **ClinVar protein index** | Gene-level (`.oga`) | Pathogenic missense by protein position (PS1, PM1, PM5) | `--source clinvar_protein` |
| **Custom VCF** | Allele-specific (`.osa`) | Any user-supplied VCF, INFO fields become the JSON object | `--source custom_vcf` |
| **Custom BED** | Interval (`.osi`) | Any user-supplied BED, name/score columns become the JSON object | `--source custom_bed` |
| **Custom (auto-detect)** | VCF or BED | Auto-detects format from `.vcf[.gz]` / `.bed[.gz]` extension | `--source custom` |

For the per-source VCF `FV_*` / tab column / JSON-key schema (pipe formats,
escaping rules, identifiers), see
[`docs/SUPPLEMENTARY_ANNOTATIONS.md`](docs/SUPPLEMENTARY_ANNOTATIONS.md).

### Custom annotation sources

You don't have to wait for a built-in parser to plug in your own data —
`sa-build` accepts arbitrary VCFs and BEDs via `--source custom_vcf`,
`--source custom_bed`, or `--source custom` (auto-detects from the
input extension). The `--name` flag becomes the JSON / column key for
the resulting database, so it shows up in output exactly like a
first-class source.

```bash
# Custom allele-level VCF — select which INFO fields to keep
fastvep sa-build --source custom_vcf \
  --name clinical --info-fields CLIN_LABEL,CLIN_SCORE \
  -i my_clinical.vcf.gz -o sa_databases/clinical

# Custom interval-level BED — score/name columns flow through automatically
fastvep sa-build --source custom_bed \
  --name myregions \
  -i my_regions.bed -o sa_databases/myregions

# Annotate as usual — both .osa and .osi in --sa-dir are picked up
fastvep annotate -i variants.vcf --gff3 genes.gff3 \
  --sa-dir sa_databases/ --output-format json
```

Allele-level custom VCFs produce a `.osa` and attach to records whose
`(pos, ref, alt)` matches. Interval-level custom BEDs produce a `.osi`
and attach via positional overlap (returning every interval that
contains the variant). Omit `--info-fields` to capture every INFO key
on every record — convenient for exploration, but the resulting JSON
objects will be heterogeneous.

## Command Reference

### `fastvep annotate`

| Flag | Description | Default |
|------|-------------|---------|
| `-i, --input` | Input VCF file (`-` for stdin; `.vcf.gz` auto-detected) | *required* |
| `-o, --output` | Output file (`-` for stdout) | `-` |
| `--gff3` | GFF3 gene annotation file. May be repeated to replicate VEP's `--merged` cache (Ensembl + RefSeq in a single run); each value may be `LABEL=path` to control the SOURCE column. | -- |
| `--fasta` | Reference FASTA file | -- |
| `--output-format` | `vcf`, `tab`, or `json` | `vcf` |
| `--hgvs` | Include HGVS notations | off |
| `--pick` | Report only the most severe consequence per variant | off |
| `--symbol` | Include gene symbol in output | off |
| `--canonical` | Include canonical-transcript flag in output | off |
| `--everything` | Turn on all common annotation flags | off |
| `--distance` | Upstream/downstream distance in bp | `5000` |
| `--buffer-size` | Variants buffered per parallel batch | `5000` |
| `--sa-dir` | Directory containing `.osa` / `.osa2` / `.osi` / `.oga` supplementary annotations | -- |
| `--sa-only` | Skip the default CSQ annotation and emit only supplementary annotations from `--sa-dir` (requires `--sa-dir`) | off |
| `--cache-dir` | Path to VEP cache directory for known-variant lookup | -- |
| `--transcript-cache` | Path to binary transcript cache file (overrides the auto-managed `<gff3>.fastvep.cache` sidecar) | -- |
| `--acmg` | Run ACMG-AMP classification (Richards 2015 + ClinGen SVI); adds `ACMG` + `ACMG_CRITERIA` to CSQ | off |
| `--acmg-config` | TOML file with custom ACMG thresholds | built-in defaults |
| `--proband` / `--mother` / `--father` | Sample names for trio analysis — enables PS2 (de novo), PM6, BP2 | -- |
| `--gene-list` | Path to a gene-panel file (one symbol or Ensembl gene ID per line). Tab output drops rows whose transcript isn't on the panel. | -- |
| `--explicit-alleles` | Add an explicit `REF` column to tab output after `Allele` | off |
| `--qc-rules` | TOML file of QC class rules; populates a `QC_CLASS` column in tab output | -- |

### `fastvep sa-build`

| Flag | Description | Default |
|------|-------------|---------|
| `--source` | Source type (clinvar, gnomad, dbsnp, cosmic, onekg, topmed, mitomap, phylop, gerp, dann, revel, spliceai, primateai, dbnsfp, omim, gnomad_genes, clinvar_protein, custom_vcf, custom_bed, custom) | *required* |
| `-i, --input` | Input file (VCF/TSV/wigFix/BED, supports .gz) | *required* |
| `-o, --output` | Output base path (creates .osa2 by default; .osa + .osa.idx under `--format osa`, .osi for BED, .oga for gene-level sources) | *required* |
| `--assembly` | Genome assembly | `GRCh38` |
| `--name` | Display + JSON-key name for `custom_*` sources | derived from input filename |
| `--info-fields` | Comma-separated INFO keys to extract for `custom_vcf` | all INFO keys |
| `--format` | On-disk format: `auto` (default), `osa` (v1), or `osa2` (v2). Every allele-level source has a v2 encoder, so `auto` builds the smaller, faster `.osa2` for all of them; gene- and interval-level sources have no v2 form and build `.oga`/`.osi` regardless. See [Choosing v1 vs v2](docs/SUPPLEMENTARY_ANNOTATIONS.md#choosing-v1-vs-v2---format) for when to force `osa`. | `auto` |

### `fastvep sa-convert`

Upgrades an existing v1 `.osa` database to v2 `.osa2` without re-downloading and
re-parsing the upstream source. Preserves the record set, every record's JSON
payload byte for byte, and the database's identity and matching semantics — see
[Converting an existing v1 database](docs/SUPPLEMENTARY_ANNOTATIONS.md#converting-an-existing-v1-database-sa-convert)
for the one caveat about column-oriented sources.

| Flag | Description | Default |
|------|-------------|---------|
| `-i, --input` | Input v1 `.osa` file (its `.osa.idx` must sit alongside it) | *required* |
| `-o, --output` | Output base path (`.osa2` appended) | input path with the extension replaced |
| `--no-progress` | Suppress periodic progress output | off |

### `fastvep filter`

| Flag | Description | Default |
|------|-------------|---------|
| `-i, --input` | Input VEP-annotated VCF | *required* |
| `-o, --output` | Output file | `-` |
| `--filter` | Filter expression (filter_vep-compatible syntax) | *required* |

Filter syntax examples:
```
IMPACT is HIGH
Consequence in missense_variant,stop_gained,frameshift_variant
AF < 0.001
IMPACT is HIGH and AF < 0.01
(IMPACT is HIGH or IMPACT is MODERATE) and not Consequence is synonymous_variant
```

### `fastvep-web` (production web server)

| Flag | Description | Default |
|------|-------------|---------|
| `--gff3` | GFF3 gene annotation file | -- |
| `--fasta` | Reference FASTA file | -- |
| `--sa-dir` | Directory containing `.osa` / `.osa2` / `.osi` / `.oga` supplementary annotations | -- |
| `--data-dir` | Directory of genome subdirectories (for multi-organism switching) | -- |
| `--port` | HTTP port (also `FASTVEP_PORT` env) | `8080` |
| `--bind` | Bind address (also `FASTVEP_BIND` env) | `0.0.0.0` |
| `--distance` | Upstream/downstream distance in bp | `5000` |
| `--max-body-size` | Max request body in bytes | `10485760` |
| `--max-concurrent` | Max concurrent annotation requests | `64` |
| `--stats-file` | JSON file to write per-request stats to (also `FASTVEP_STATS_FILE` env) | -- |

### `fastvep cache`

Pre-builds a binary transcript cache for fast startup. Accepts the same
multi-`--gff3` / `LABEL=path` syntax as `annotate`, so a merged Ensembl
+ RefSeq cache can be built once and reused via `--transcript-cache`.

| Flag | Description | Default |
|------|-------------|---------|
| `--gff3` | GFF3 annotation file(s). Repeatable; each value may be `LABEL=path`. | *required* |
| `--fasta` | Reference FASTA (for pre-building spliced sequences) | -- |
| `--synonyms` | VEP-style `chr_synonyms.txt` mapping equivalent contig names (e.g. RefSeq `NC_000017.11` ↔ `17`). Reconciles mixed Ensembl/RefSeq naming against the FASTA. Only effective with `--fasta`. | -- |
| `-o, --output` | Output cache file path | *required* |

## Output Formats

### VCF Output

Consequence annotations are added as a `CSQ` field in the INFO column with 49 pipe-delimited fields matching Ensembl VEP's extended format, plus fastVEP-specific `ACMG` and `ACMG_CRITERIA` slots when `--acmg` is set. When supplementary annotation databases are loaded with `--sa-dir`, fastVEP also emits VCF-compatible INFO projections for supported fastSA sources: standard `SpliceAI` for SpliceAI databases, and fastVEP-specific `FV_*` fields such as `FV_CLINVAR`, `FV_GNOMAD`, `FV_DBSNP`, `FV_REVEL`, and gene-level `FV_OMIM`.

The VCF output never embeds raw JSON in INFO values. Use `--output-format json` for the richest structured representation of all supplementary annotation objects.

### Tab Output

One line per variant-transcript-allele combination with 17 columns.

### JSON Output

Structured JSON with `transcript_consequences` array per variant, including supplementary annotations from SA providers (ClinVar, gnomAD, etc.) and gene-level annotations.

## Consequence Types

fastVEP predicts 49 consequence types organized by impact:

| Impact | Consequences |
|--------|-------------|
| **HIGH** | transcript_ablation, splice_acceptor_variant, splice_donor_variant, stop_gained, frameshift_variant, stop_lost, start_lost, transcript_amplification, TFBS_ablation, regulatory_region_ablation |
| **MODERATE** | inframe_insertion, inframe_deletion, missense_variant, protein_altering_variant, regulatory_region_amplification, TFBS_amplification |
| **LOW** | splice_region_variant, splice_donor_5th_base_variant, splice_donor_region_variant, splice_polypyrimidine_tract_variant, synonymous_variant, start_retained_variant, stop_retained_variant, incomplete_terminal_codon_variant |
| **MODIFIER** | coding_sequence_variant, 5_prime_UTR_variant, 3_prime_UTR_variant, non_coding_transcript_exon_variant, intron_variant, upstream_gene_variant, downstream_gene_variant, intergenic_variant, copy_number_change, copy_number_increase, copy_number_decrease, short_tandem_repeat_change, transcript_variant, and others |

## Architecture

```
crates/
  fastvep-core/         # Core types: Consequence (49 SO terms), VariantType, Allele, Impact
  fastvep-genome/       # Transcript, Exon, Gene, CodonTable, mitochondrial codon table
  fastvep-cache/        # GFF3 parser, FASTA reader, annotation providers, regulatory regions
  fastvep-consequence/  # Consequence prediction: small variants + SV predictor
  fastvep-hgvs/         # HGVS nomenclature generation (c., p., g.)
  fastvep-io/           # VCF parser (incl. SVs), output formatters, multi-sample parsing
  fastvep-filter/       # Filter engine: lexer, parser, evaluator (filter_vep-compatible)
  fastvep-sa/           # Supplementary annotation format (fastSA):
                       #   v1 (.osa): zstd block compression, binary search,
                       #     shared byte-budgeted LRU cache of decompressed blocks
                       #     (thread-scaled; override via FASTVEP_SA_CACHE_BYTES)
                       #   v2 (.osa2): echtvar-inspired chunked ZIP with Var32 encoding,
                       #     parallel u32 value arrays, delta encoding, lock-free mmap
                       #     reads + shared byte-budgeted LRU chunk cache
                       #   .osi: interval-level annotations (BED-derived), positional overlap
                       #   .oga: gene-level annotations (OMIM, gnomAD constraint, ClinVar
                       #     protein index)
                       # Source parsers: ClinVar, gnomAD (incl. v4.1 joint), dbSNP, COSMIC,
                       # 1000G, TOPMed, MitoMap, PhyloP, GERP, DANN, REVEL, SpliceAI,
                       # PrimateAI, dbNSFP, plus user-supplied custom_vcf / custom_bed.
  fastvep-annotate/     # Shared annotation pipeline (used by CLI batch and web server):
                       #   variant overlap, consequence prediction, HGVS, SA/gene SA
                       #   provider loading
  fastvep-classification/ # ACMG-AMP variant classification engine (Richards 2015 +
                       #   ClinGen SVI). 28 criteria, trio/compound-het support,
                       #   configurable thresholds via TOML
  fastvep-cli/          # CLI binary. `src/pipeline/` holds one module per
                       #   subcommand: annotate, cache_build, sa_build, filter,
                       #   custom, plus pick (the --pick criteria)
  fastvep-web/          # Production web server (axum/tokio): async, multi-connection,
                       #   genome switching, SA integration, rate limiting
web/                   # Web GUI (HTML/CSS/JS, embedded in both server binaries)
tests/                 # Test data: chr1 (OR4F5) and chr17 (BRCA1) VCF + GFF3
```

## Running Tests

```bash
cargo test --workspace          # ~1,000 tests
cargo test -p fastvep-consequence  # Consequence prediction tests (incl. SV)
cargo test -p fastvep-filter       # Filter engine tests
cargo test -p fastvep-sa           # Supplementary annotation format tests
```

## Performance Benchmarks

Benchmarked on Apple M-series (ARM64), release build with LTO. Median of 3 runs, full Ensembl annotations with FASTA and HGVS.

### Multi-Organism Throughput (Gold-Standard Datasets)

| Organism | Transcripts | Variants | Source | Time | Throughput |
|----------|-------------|----------|--------|------|------------|
| Yeast (R64, full genome) | 7,036 | 260,526 | Ensembl/SGD | 3.0s | **85,934 v/s** |
| Drosophila (BDGP6, full) | 35,442 | 4,438,427 | DGRP2 | 57.3s | **77,486 v/s** |
| Arabidopsis (TAIR10, full) | 54,013 | 12,883,854 | 1001 Genomes | 168.7s | **76,378 v/s** |
| Mouse (GRCm39, full genome) | 142,626 | 26,062,054 | MGP CAST/EiJ | 338.0s | **77,113 v/s** |
| Human full WGS (GRCh38) | 508,530 | 4,048,342 | GIAB HG002 | 86.3s | **46,917 v/s** |

### vs. Ensembl VEP v115.1 (head-to-head, single-thread, GIAB HG002)

Both tools single-threaded on identical GRCh38 Ensembl 115 GFF3 + FASTA with `--hgvs --symbol --canonical`. fastVEP uses its binary transcript cache; Ensembl VEP runs via Docker (`ensemblorg/ensembl-vep:release_115.1`) in `--gff` mode against a bgzip+tabix-indexed GFF3. VEP's full-genome time is measured per chromosome and summed — a single-process whole-genome run OOMs above ~2.7M variants at the 7.65 GiB Docker default, and the per-chromosome sum keeps VEP's memory bounded (its best case).

**Full human WGS (4,048,342 variants), single-thread:**

| Annotation | fastVEP | Ensembl VEP | Speedup |
|-----------|---------|-------------|---------|
| Consequence only | 197.8s | 4,621s (~77 min) | **23.4×** |
| + ClinVar | 197.4s | 4,803s | **24.3×** |
| + gnomAD + ClinVar | 218.5s | 4,905s | **22.4×** |

With fastVEP's default multi-threading (10 cores), the same full WGS completes in **93s** (consequence-only) to **104.5s** (+gnomAD+ClinVar) on this system.

**chr22 scaling, single-thread** — fastVEP carries a fixed ~2.7s binary-cache load, so VEP (which tabix-fetches only overlapping GFF regions) is faster below ~1–2K variants; beyond that fastVEP's per-variant throughput dominates:

| Variants | fastVEP | VEP | Speedup |
|----------|---------|-----|---------|
| 1,000 | 2.68s | 0.45s | 0.17× (VEP faster) |
| 5,000 | 2.86s | 6.29s | 2.2× |
| 10,000 | 3.25s | 13.41s | 4.1× |
| 50,284 | 5.09s | 93.95s | **18.5×** |

| Resource | Ensembl VEP | fastVEP |
|----------|-------------|---------|
| Peak memory (100K variants) | ~500 MB | **2.8 MB** |
| Binary / install size | ~200 MB installed | **3.3 MB** |
| Dependencies | Perl 5.22+, DBI, 10+ CPAN modules | **None** |

### Supplementary annotation format (fastSA v2)

The `.osa2` format (echtvar-inspired chunked binary; `--format auto` default) stores real annotation databases compactly, with annotation output **byte-identical** to the v1 `.osa` format (md5-verified):

| Source | v1 (.osa) | v2 (.osa2) | Ratio |
|--------|-----------|-----------|-------|
| ClinVar (4.44M records) | 48.1 MB | 39.0 MB | 0.81× |
| REVEL (chr1, 8.25M records) | 40.5 MB | 19.0 MB | 0.47× |
| REVEL (chr22, 1.78M records) | 8.7 MB | 3.8 MB | 0.44× |
| SpliceAI (chr22, 828K records) | 4.8 MB | 2.7 MB | 0.56× |
| gnomAD (chr22, 852K records) | 15.4 MB | 12.9 MB | 0.84× |
| PhyloP (chr22, 852K records) | 3.1 MB | 1.4 MB | 0.45× |

Annotating 308,740 real chr22 variants against all five chr22/ClinVar databases at once (8.3M records total): **12.7 s from the `.osa` set, 8.5 s from the `.osa2` set — 1.50×**, md5-identical output. Every record of all five databases was queried through both readers and compared byte for byte (8,318,903 keys, zero differences).

## VEP Concordance

fastVEP's consequence and HGVS output is a port of Ensembl VEP's own model, validated against
**real Ensembl VEP v115.1** (Docker `ensemblorg/ensembl-vep:release_115.1`, `--gff` mode, the same
GRCh38 Ensembl 115 GFF3 and FASTA) on three datasets, because they exercise different code paths.

**Field-level, all annotation types.**
All 23 compared CSQ fields agree on 100% of 2,340 shared (allele, transcript) pairs from the
173-variant VEP example set (`validation/run_validation.sh`), across all 12 consequence types the
set contains.

**Consequences and HGVS under indels, MNVs and splice sites.**
The example set is SNV-only, so the harder shapes are validated on a 6,600-variant stratified
sample of the ClinVar 2-star+ set - **150,725 matched (variant, allele, transcript) rows**, 72,632
of them coding. That sample is built to be hard: 54.5% of its variants are not SNVs, against 7.2%
of the ClinVar 2-star+ set it is drawn from, so these are the rates on the shapes that disagree.

| Field | Scope | Agreement |
|---|---|---:|
| `Amino_acids` | coding rows | **100 %** |
| `Codons` | coding rows | **100 %** |
| Splice terms | all rows | **100 %** |
| Consequence terms | coding rows | 99.86 % |
| Whole consequence set | all rows | 99.92 % |
| `IMPACT` | all rows | 99.93 % |
| `HGVSc` | all rows | 99.60 % |
| `HGVSp` | all rows | 99.43 % |

**Genome-wide, on an ordinary callset.**
The same comparison over a systematic 1-in-200 sample of the GIAB HG002 WGS callset - 20,241
variants, **122,317 matched rows** - is the other end of that range: the consequence set and
`IMPACT` agree on **100.000 %** of rows, `HGVSp` on **99.985 %**, `HGVSc` on **99.879 %**.
Most of the 148 remaining `HGVSc` rows come from multi-allelic VCF records (1.18% of this
callset), whose per-allele trimming against the reference is a known gap - see
[docs/VEP_DIVERGENCE.md](docs/VEP_DIVERGENCE.md).

**HGVSp under in-frame indels.**
Protein-level normalisation is additionally checked on 400 ClinVar in-frame deletions run through
both tools: 99.17% exact string agreement on 5,192 (variant, transcript) pairs.

### Divergences

Agreement is the evidence that the port is faithful, not the objective. The output is a prediction
a clinician may act on, so where Ensembl is demonstrably wrong in a way that changes a call,
fastVEP is right instead.

There are five such places, and one of them matters clinically: Ensembl reports a frameshift that
introduces a premature stop as `inframe_insertion,stop_retained_variant`, MODERATE. BRCA1
`c.5030_5033dup` - ClinVar 3-star Pathogenic - is one of 34 such variants in the ClinVar 2-star+
set. fastVEP reports `stop_gained,frameshift_variant`, HIGH, so those 34 keep PVS1.

Every divergence, its cause in the Ensembl source, its row count, and the list of fastVEP's own
remaining gaps in the other direction are in **[docs/VEP_DIVERGENCE.md](docs/VEP_DIVERGENCE.md)**.

## Citation

If you use fastVEP in your research, please cite:

**fastVEP: A Fast, Comprehensive Variant Effect Predictor Written in Rust**  
Kuan-lin Huang  
*bioRxiv* (2026)  
doi: [https://doi.org/10.64898/2026.04.14.718452](https://doi.org/10.64898/2026.04.14.718452)  
URL: [https://www.biorxiv.org/content/10.64898/2026.04.14.718452v1](https://www.biorxiv.org/content/10.64898/2026.04.14.718452v1)

## License

Apache License 2.0

## Acknowledgements

fastVEP is inspired by [Ensembl VEP](https://www.ensembl.org/info/docs/tools/vep/index.html) by EMBL-EBI and [Illumina Nirvana](https://github.com/Illumina/Nirvana). The consequence prediction logic follows the Sequence Ontology term definitions and the Ensembl variant annotation framework. The supplementary annotation system (fastSA v2) incorporates algorithms and encoding strategies from [echtvar](https://github.com/brentp/echtvar).
