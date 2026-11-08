# fastVEP

A high-performance Variant Effect Predictor written in Rust. fastVEP predicts the functional consequences of genomic variants (SNPs, insertions, deletions, structural variants) on genes, transcripts, and protein sequences, with direct integration of clinical and population databases.

fastVEP is inspired by and aims to be compatible with [Ensembl VEP](https://www.ensembl.org/info/docs/tools/vep/index.html) and [Illumina Nirvana](https://github.com/Illumina/Nirvana), while delivering significantly better performance through Rust's zero-cost abstractions and native parallelism.

**Try it now:** A hosted web server is available at [fastVEP.org](https://fastVEP.org) - paste VCF data and get annotated results instantly, no installation required.
On 7 November 2026 it had served **3,427 genome analysis sessions** and annotated **160,983 variants**.
Since 8 November 2026 it runs a build of the 0.5.0 release candidate and answers from nine genomes: human (GRCh38 and GRCh37), mouse, rat, zebrafish, fruit fly, *C. elegans*, yeast and Arabidopsis.

The hosted instance is for interactive use in the browser. For programmatic access, run your own server: see the **[REST API guide](docs/API.md)**.

**Questions and discussion:** join the [fastVEP Slack](https://join.slack.com/t/fastvep/shared_invite/zt-4g20tqjvo-wGyIEnN~Ncrc4zXiu_1XJA).
Bug reports and feature requests belong in [GitHub issues](https://github.com/Huang-lab/fastVEP/issues).

## Features

- **Variant Consequence Prediction** — Classifies variants using [Sequence Ontology](http://www.sequenceontology.org/) terms: 49 are defined and the predictors emit 38 of them (missense, frameshift, splice donor, copy_number_change, transcript_ablation, etc.; `manuscript/data/so_terms.csv` lists each)
- **Structural Variant Support** — Full SV pipeline: `<DEL>`, `<DUP>`, `<INV>`, `<CNV>`, `<BND>`, `<INS>`, `<STR>` with SV-specific consequence prediction
- **Supplementary Annotations** — Direct integration with ClinVar, gnomAD, dbSNP, COSMIC, 1000 Genomes, TOPMed, MitoMap via the native fastSA format (v1: zstd block compression; v2: echtvar-inspired chunked ZIP with Var32 encoding, parallel u32 value arrays, delta encoding). Both use a shared, thread-scaled, byte-budgeted LRU cache of decompressed blocks/chunks and lock-free reads, so a dense whole-genome source (e.g. SpliceAI, dbSNP) queried in parallel stays fast without re-decompressing the same data
- **Prediction Scores** — PhyloP, GERP, REVEL, SpliceAI, PrimateAI, DANN, AlphaMissense conservation and pathogenicity scores; SIFT/PolyPhen via dbNSFP
- **Gene-Level Annotations** — OMIM phenotypes, gnomAD gene constraint (pLI, LOEUF), ClinGen gene-disease validity
- **Filter Engine** — Expression-based filtering compatible with VEP's filter_vep syntax
- **HGVS Nomenclature** — Generates HGVSg, HGVSc, and HGVSp notations with 3' normalization
- **Multiple Output Formats** — VCF (with 49-field CSQ), tab-delimited, JSON
- **Multi-Sample Support** — Parse FORMAT/GT/DP/GQ/AD fields per sample with genotype classification
- **Regulatory regions** — not annotated yet: the library has a parser and an interval index for Ensembl regulatory-build GFF3 (`fastvep-cache::regulatory`), but no command uses them, `--regulatory` is refused, and the five regulatory columns of the CSQ are empty
- **Mitochondrial Support** — Circular coordinate handling, vertebrate mitochondrial codon table (NCBI table 2)
- **Custom Annotations** — User-provided VCF (`--source custom_vcf`) and BED (`--source custom_bed`) files; `.osi` interval databases load alongside `.osa` via `--sa-dir`
- **ACMG-AMP Classification** — `--acmg` runs the full Richards 2015 + ClinGen SVI rule set (28 criteria, configurable thresholds, trio / compound-het support via `--proband`/`--mother`/`--father`)
- **VEP `--merged` Cache** — `--gff3` is repeatable on `annotate` and `cache`; combine Ensembl + RefSeq in a single run with per-transcript SOURCE labels
- **`--sa-only` Mode** — Skip the default CSQ pipeline and emit only supplementary annotations, useful for re-annotating already-annotated VCFs
- **Gzipped inputs** — `annotate` auto-detects gzip-compressed VCF, GFF3, and FASTA by magic bytes, with a `.gz` / `.bgz` filename fallback. A `.gff3.gz` that has a `.tbi` is still read through tabix. A `.fai` applies only to an uncompressed FASTA, because its offsets address that file; a gzipped FASTA is loaded into memory
- **Web Interface** — Built-in web GUI for interactive variant annotation
- **GFF3 Annotation Support** — Load gene models from standard GFF3 files (any organism)

## Quick Start

### 1. Install Rust (if you don't have it)

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
source "$HOME/.cargo/env"
```

Building also needs a C compiler, since the `zstd` and `mimalloc` crates compile C: the command line tools of Xcode on macOS (`xcode-select --install`), `build-essential` on Debian and Ubuntu.

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

#### Alternative: download a prebuilt binary

Each [release](https://github.com/Huang-lab/fastVEP/releases) carries `fastvep` and `fastvep-web` in one archive per platform, with a `.sha256` file beside it:

| Platform | Archive |
| --- | --- |
| Linux, x86_64 | `fastVEP-<version>-x86_64-unknown-linux-gnu.tar.gz` |
| Linux, aarch64 (releases after v0.4.0) | `fastVEP-<version>-aarch64-unknown-linux-gnu.tar.gz` |
| macOS, Apple silicon | `fastVEP-<version>-aarch64-apple-darwin.tar.gz` |

```bash
VERSION=v0.4.0                    # a release tag
TARGET=x86_64-unknown-linux-gnu   # or aarch64-unknown-linux-gnu, aarch64-apple-darwin
BASE=https://github.com/Huang-lab/fastVEP/releases/download/$VERSION
curl -LO "$BASE/fastVEP-$VERSION-$TARGET.tar.gz"
curl -LO "$BASE/fastVEP-$VERSION-$TARGET.tar.gz.sha256"
sha256sum -c "fastVEP-$VERSION-$TARGET.tar.gz.sha256"   # shasum -a 256 -c on macOS
tar -xzf "fastVEP-$VERSION-$TARGET.tar.gz"
"fastVEP-$VERSION-$TARGET/fastvep" --version
```

The Linux binaries need glibc 2.34 or newer: Ubuntu 22.04, Debian 12, RHEL 9 and its rebuilds (Rocky Linux, AlmaLinux), and anything newer.
The release workflow runs each Linux archive in Rocky Linux 9 before it publishes it, so that a build which asks for a newer glibc stops the release.
On an older system (Ubuntu 20.04, Debian 11, RHEL 8), build from source as above.

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
# Annotate the 8 test variants (coding, 5' UTR, 3' UTR, intronic and downstream variants, and a deletion, in BRCA1 and TP53)
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
fastvep sa-build --source gnomad --input gnomad.genomes.v4.1.sites.chr1.vcf.bgz --output gnomad_chr1

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
# Starts with no gene model: every variant is an intergenic_variant until one is loaded
fastvep-web

# With your own data
fastvep-web --gff3 Homo_sapiens.GRCh38.115.gff3 --fasta Homo_sapiens.GRCh38.dna.primary_assembly.fa

# With supplementary annotations (ClinVar, gnomAD, etc.)
fastvep-web --gff3 genes.gff3 --fasta ref.fa --sa-dir /path/to/sa_databases/
```

Open http://localhost:8080 in your browser. The web interface lets you paste VCF data and view results in an interactive table. With `--data-dir` the page offers every genome in that directory, and each analysis names the genome it wants, so choosing one needs no flag and changes nobody else's answers.
ACMG classification is offered for a genome only when its database directory holds all nine sources the classifier reads; for any other the checkbox is disabled and names what is missing, and the API refuses with `422` (`--acmg-allow-partial` answers anyway and says so in every response; [docs/API.md](docs/API.md#acmg-needs-the-reference-stack)).
Uploading your own GFF3 in the page replaces the server's gene model for every client, so that needs `--allow-model-replacement`, which is off by default.

> **Note:** `fastvep-web` is a separate production-quality binary (axum/tokio, async, multi-connection). The legacy `fastvep web` command still works but is single-threaded.

### REST API (self-hosted)

The same binary serves a JSON API over the same annotation engine, so anything the web interface can do a script can do:

```bash
curl -s -X POST http://localhost:8080/api/annotate \
  -H 'Content-Type: application/json' \
  -d '{"vcf": "17\t43124025\t.\tA\tG\t50\tPASS\t.", "acmg": true}'
```

It is designed for self-hosting on a workstation or local network, which is also the only supported way to use it: please do not point scripts at fastVEP.org, which is a single small machine sized for interactive browser use. Bulk VCFs belong in `fastvep annotate`, not in an HTTP request.

Full setup, endpoint reference, and response shapes: **[docs/API.md](docs/API.md)**.

## Local Setup Guide

This section walks through setting up fastVEP with full annotation capabilities (gene models, reference sequence, and supplementary databases like ClinVar and gnomAD).

### Step 1: Download reference data

```bash
mkdir -p data && cd data

# Gene models (GFF3) — pick your organism
# Human GRCh38. `annotate` reads the `.gz` directly; gunzip is optional.
wget https://ftp.ensembl.org/pub/release-115/gff3/homo_sapiens/Homo_sapiens.GRCh38.115.gff3.gz

# Reference FASTA (needed for HGVS and sequence context).
# `annotate` reads `.fa.gz` directly, into memory. For a whole genome,
# gunzip and index it so the run memory-maps the file instead.
wget https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# Create FASTA index (enables memory-mapped access — important for large genomes)
samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa
```

### Step 2: Build supplementary annotation databases

For every source, for both GRCh38 and GRCh37, with download URLs, sizes, which ones need an account, and a script that does the downloading and building, see [docs/SA_DATABASES.md](docs/SA_DATABASES.md).
The quickest way to try supplementary annotation is `scripts/build-sa-databases.sh --assembly GRCh38`, which builds the small open sources in a few minutes.

Each supplementary database (ClinVar, gnomAD, etc.) is built in **two steps** — *download the source file*, then run `fastvep sa-build` to convert it into a fastSA database (a `.osa2` file by default; a `.osa` + `.osa.idx` pair under `--format osa`). **`sa-build` is a converter, not a downloader.** Handed a missing, truncated or wrong-build file it finds no record to write, and it stops with an error that names the usual causes instead of writing a database of nothing (`--allow-empty` writes one anyway, for the shard of a chromosome that has no record in your regions). After each build, check that the file size matches the expected magnitude (column below); a few-KB database is the tell that the source file wasn't real. A release that dates itself (a VCF with a `##fileDate=` line, as ClinVar's does) is stamped with that date, which `GET /api/status` reports, so that a database says how old it is.

```bash
mkdir -p sa_databases

# ── ClinVar - clinical variant significance ──
# Download (193 MB as of 2026-08-31; re-released weekly and growing)
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
# Build (expect ~40-60 MB)
fastvep sa-build --source clinvar -i clinvar.vcf.gz -o sa_databases/clinvar --assembly GRCh38

# ── gnomAD v4.1 — population allele frequencies ──
# Per-chromosome sites VCFs. Measured 2026-08-31, chr1-22 X Y: the exomes
# release is 198 GB and the genomes release is 563 GB, so most setups want
# tabix-extracted target regions rather than the whole thing. See
# docs/ACMG_SETUP.md for the region-restricted recipe, and for the joint release
# of the exomes and the genomes, which gives more calls in direction on the
# ACMG benchmark than the exomes alone.
GNOMAD=https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/exomes
wget "$GNOMAD/gnomad.exomes.v4.1.sites.chr1.vcf.bgz"
fastvep sa-build --source gnomad -i gnomad.exomes.v4.1.sites.chr1.vcf.bgz -o sa_databases/gnomad_chr1 --assembly GRCh38

# ── dbSNP — variant identifiers ──
wget https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz
fastvep sa-build --source dbsnp -i GCF_000001405.40.gz -o sa_databases/dbsnp --assembly GRCh38

# ── COSMIC — somatic mutations (requires license) ──
# https://cancer.sanger.ac.uk/cosmic/download
fastvep sa-build --source cosmic -i CosmicCodingMuts.vcf.gz -o sa_databases/cosmic --assembly GRCh38
```

**Verify before moving on:**

```bash
ls -la sa_databases/*.osa2 sa_databases/*.osa
# Expected: clinvar ~50 MB; gnomad hundreds of MB to several GB depending on
# how much of the release you built; dbsnp ~5 GB.
# Anything < 1 MB usually means an empty build - re-check the source file.
```

> **Setting up for ACMG-AMP classification?** Read the **[ACMG Setup Guide](docs/ACMG_SETUP.md)** instead of this section.
> The classifier reads a specific nine-source stack: gnomAD, ClinVar, REVEL, SpliceAI, PhyloP, ClinGen Gene-Disease Validity, gnomAD gene constraints, the ClinVar protein/splice index, and RepeatMasker.
> All of it is public and none of it needs an account, because **SpliceAI and PhyloP are distilled from gnomAD's own `spliceai_ds_max` and `phylop` INFO fields rather than downloaded from Illumina BaseSpace and UCSC**, and the disease-gene source is ClinGen Gene-Disease Validity rather than OMIM.
> The guide gives the verified download URL, build command and measured size for each, explains what silently stops working when one is missing, and ends with `scripts/check_acmg_stack.py`, which annotates eight probe variants and tells you which sources actually answered.
> An incomplete `--sa-dir` does not error; it turns pathogenic variants into VUS.
> `fastvep annotate --acmg` prints which sources are missing, once, on stderr, and `fastvep-web` refuses to classify a genome that lacks any of them.

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

Fourteen flags also take an environment variable, for container deployments: `FASTVEP_GFF3`, `FASTVEP_FASTA`, `FASTVEP_SA_DIR`, `FASTVEP_DATA_DIR`, `FASTVEP_PORT`, `FASTVEP_BIND`, `FASTVEP_STATS_FILE`, `FASTVEP_PRELOAD_GENOMES`, `FASTVEP_ACMG_ALLOW_PARTIAL`, `FASTVEP_ALLOW_MODEL_REPLACEMENT`, `FASTVEP_MAX_PARALLEL`, `FASTVEP_QUEUE_TIMEOUT_SECONDS`, `FASTVEP_MAX_RESPONSE_BYTES` and `FASTVEP_MAX_REQUEST_SECONDS`.
The three switches accept `true`, `1`, `yes` or `on`.

### Multi-organism setup

To serve multiple genomes from the web interface, organize data into subdirectories and use `--data-dir`.
Each subdirectory is one genome, and the server detects the GFF3, FASTA and SA files inside it.

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

Each analysis names its genome, from the dropdown in the web interface or in the `genome` field of `POST /api/annotate`, and is answered from that genome's own model.
Nothing is replaced, so choosing a genome needs no `--allow-model-replacement` and cannot change what any other client is answered against.
A genome is loaded the first time it is asked for and then kept; `--preload-genomes` loads them all in the background at startup instead, so no visitor waits.
An unknown name is refused rather than answered from the default genome.
Details and response shapes: [docs/API.md](docs/API.md#choosing-the-genome).

Each genome is read from its own GFF3, FASTA and `sa/` subdirectory.
`--sa-dir` belongs to the genome given by `--gff3`: it is never lent to another genome of `--data-dir`, so a mouse variant is never annotated with a human ClinVar record.
If that directory doesn't exist, the server starts without SA (no error).

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

Each transcript carries its source through to the VCF and JSON output:
the `SOURCE` field in the VCF CSQ string and the `source` key in JSON.
(Tab output has no `SOURCE` column.)
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
| **gnomAD** | Allele-specific | Population frequencies (ten groups: AFR, AMR, ASJ, EAS, FIN, MID, NFE, OTH, SAS and the remainder), allele counts | `--source gnomad` |
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
| **OMIM / ClinGen GDV** | Gene-level (`.oga`) | Disease-gene validity and mode of inheritance, driving PVS1, PP2, PM1 (validity) and PM2, BS2, PM3, BP2 (inheritance) in ACMG | `--source omim` |
| **gnomAD constraint** | Gene-level (`.oga`) | gnomAD constraint metrics (pLI, LOEUF) for PVS1, PP2, BP1 | `--source gnomad_genes` |
| **ClinVar protein index** | Gene-level (`.oga`) | Pathogenic missense by protein position (PS1, PM1, PM5, BP1) | `--source clinvar_protein` |
| **Custom VCF** | Allele-specific (`.osa2`) | Any user-supplied VCF, INFO fields become the JSON object | `--source custom_vcf` |
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
  --sa-dir sa_databases/ --output-format vcf
```

Allele-level custom VCFs produce a `.osa2` and attach to records whose
`(pos, ref, alt)` matches. Interval-level custom BEDs produce a `.osi`
and attach via positional overlap (returning every interval that
contains the variant). Omit `--info-fields` to capture every INFO key
on every record — convenient for exploration, but the resulting JSON
objects will be heterogeneous.

Every output format carries the annotation. JSON files it under the
`--name` key; VCF and tab project it under `FV_<NAME>`, as a
self-describing `ALLELE|KEY=VALUE&KEY=VALUE` value rather than the fixed
pipe layout a built-in source uses, because a custom source's field set
is whatever its records carry:

```console
$ fastvep sa-build --source custom_vcf --name mypanel -i panel.vcf -o sa/mypanel
$ fastvep annotate -i variants.vcf --gff3 genes.gff3 --sa-dir sa/ --output-format vcf | grep FV_MYPANEL
##INFO=<ID=FV_MYPANEL,Number=.,Type=String,Description="fastVEP custom annotations from mypanel. ...">
17  43124090  .  A  G  30  PASS  FV_MYPANEL=G|AC=7&AF=0.005814&AN=1204;CSQ=...
```

A local frequency panel needs nothing beyond this: `AC` / `AN` / `AF`
come through by name, and `bcftools query -f '%INFO/FV_MYPANEL\n'` reads
them back. See
[docs/SUPPLEMENTARY_ANNOTATIONS.md](docs/SUPPLEMENTARY_ANNOTATIONS.md#custom-sources)
for the parsing rules and the escaping table.

## Command Reference

### `fastvep annotate`

| Flag | Description | Default |
|------|-------------|---------|
| `-i, --input` | Input VCF file (`-` for stdin; `.vcf.gz` auto-detected) | *required* |
| `-o, --output` | Output file (`-` for stdout) | `-` |
| `--gff3` | GFF3 gene annotation file. May be repeated to replicate VEP's `--merged` cache (Ensembl + RefSeq in a single run); each value may be `LABEL=path` to control the SOURCE column. | -- |
| `--fasta` | Reference FASTA file | -- |
| `--output-format` | `vcf`, `tab`, or `json` | `vcf` |
| `--hgvs` | Add `HGVSc`, `HGVSp` and `HGVS_OFFSET`: as columns in tab output, as values in the VCF and JSON output. VEP's `--hgvs` | off |
| `--hgvsg` | Add `HGVSg` (`chr10:g.10414C>A`) to tab and VCF output, and the `hgvsg` key of JSON. VEP's `--hgvsg`: an insertion or deletion is written at its most 3' position along the genome and an insertion that repeats its neighbours is a `dup`, which reads the reference, so give `--fasta`. See [docs/VEP_DIVERGENCE.md](docs/VEP_DIVERGENCE.md#11-an-allele-of-a-record-with-several-alts-is-described-as-it-would-be-alone) | off |
| `--pick` | Report one consequence per variant, chosen by `--pick-order`, keeping that transcript's annotation for every alt allele | off |
| `--pick-allele` | One consequence per (variant, allele). VEP's `--pick_allele` | off |
| `--pick-allele-gene` | One consequence per (variant, allele, gene). VEP's `--pick_allele_gene` | off |
| `--flag-pick` | Report *every* consequence and mark the one `--pick` would have kept with `PICK=1`. VEP's `--flag_pick` | off |
| `--flag-pick-allele` | As `--flag-pick`, marking the pick per allele | off |
| `--flag-pick-allele-gene` | As `--flag-pick`, marking the pick per allele and gene | off |
| `--pick-order` | Criteria order used by `--pick`, VEP `--pick_order` syntax. Severity (`rank`) is last in the default, so this is not "most severe" -- see [docs/ACMG.md](docs/ACMG.md#which-transcript---pick-reports) | VEP's default |
| `--symbol` | Add `SYMBOL`, `SYMBOL_SOURCE` and `HGNC_ID` to tab output (the VCF `CSQ` always has them). VEP's `--symbol` | off |
| `--canonical` | Add `CANONICAL` (`YES`) to tab output (the VCF `CSQ` always has it). VEP's `--canonical` | off |
| `--biotype` | Add `BIOTYPE` to tab output (the VCF `CSQ` always has it). VEP's `--biotype` | off |
| `--numbers` | Add `EXON` and `INTRON` (`3/10`) to tab output (the VCF `CSQ` always has them). VEP's `--numbers` | off |
| `--variant-class` | Add `VARIANT_CLASS` (`SNV`, `insertion`, `deletion`, ...) to tab and VCF output. VEP's `--variant_class` | off |
| `--show-ref-allele` | Add `REF_ALLELE` to tab output (the VCF `CSQ` always has it). VEP's `--show_ref_allele`. (`--explicit-alleles` puts `REF` right after `Allele` instead) | off |
| `--total-length` | Write the length each position is counted in after `cDNA_position`, `CDS_position` and `Protein_position`, `14/203`, in tab and VCF output (JSON writes numbers). VEP's `--total_length` | off |
| `--no-escape` | Do not percent-escape HGVS strings: a synonymous `HGVSp` ends in `=` and not `%3D`, and the space of a gene symbol in the VCF `CSQ` is `_` and not `%20`. VEP's `--no_escape` | off |
| `--force-overwrite` | Accepted so that a `vep` command line is not refused, and does nothing: `fastvep annotate` always writes over its output file, which is what VEP's `--force_overwrite` asks for | -- |
| `--exclude-predicted` | Leave out RefSeq's predicted transcripts, `XM_` and `XR_`. VEP's `--exclude_predicted`, which is `not stable_id match ^X._`, and so is this. A variant left with no transcript is an `intergenic_variant` | off |
| `--exclude-readthrough` | Leave out the transcripts the annotation tags `readthrough_transcript`, which span two genes. VEP's Ensembl cache is built without them and a GFF3 has them, so give it to match `vep --cache` | off |
| `--exclude-artifact` | Leave out transcripts of biotype `artifact`. VEP's Ensembl cache is built without them and a GFF3 has them, so give it to match `vep --cache` | off |
| `--af` | Add `AF`, the 1000 Genomes combined allele frequency of the known variant, to the tab output and the CSQ. VEP's `--af`. Written as VEP writes it: per input allele, to four decimals, with the other allele of a two-allele record taken as `1 - AF`. Reads the variation cache, so it needs `--cache-dir`, as do the six flags below | off |
| `--af-1kg` | Add `AF`, `AFR_AF`, `AMR_AF`, `EAS_AF`, `EUR_AF` and `SAS_AF`. VEP's `--af_1kg` | off |
| `--af-gnomade` | Add `gnomADe_AF` and the gnomAD exomes populations (`gnomADe_AFR_AF` ... `gnomADe_SAS_AF`). VEP's `--af_gnomade`; its old name `--af-gnomad` is accepted | off |
| `--af-gnomadg` | Add `gnomADg_AF` and the gnomAD genomes populations (`gnomADg_AFR_AF` ... `gnomADg_SAS_AF`). VEP's `--af_gnomadg` | off |
| `--max-af` | Add `MAX_AF`, the highest frequency of any 1000 Genomes or gnomAD population (the combined ones are left out), and `MAX_AF_POPS`, the populations that have it. VEP's `--max_af` | off |
| `--check-existing` | Add `CLIN_SIG`, `SOMATIC` and `PHENO` from the known variants, to the tab output and the CSQ, as VEP's `--check_existing` does. Every frequency flag and `--pubmed` imply it, and so does `--cache-dir`, which has no other use here, so the flag is only needed to be typed as VEP has it. `Existing_variation` is filled whenever `--cache-dir` is given | on with `--cache-dir` |
| `--pubmed` | Add `PUBMED`, the publications that cite the known variants, to the tab output and the CSQ. VEP's `--pubmed` | off |
| `--clin-sig-allele` | `1`: `CLIN_SIG` is what a known variant says about the input allele, where it says it; `0`: what it says about the variant. VEP's `--clin_sig_allele` | `1` |
| `--everything` | Turn on `--hgvs --symbol --numbers --canonical --biotype --variant-class`: VEP's `--everything`, restricted to the columns fastVEP writes (not `--hgvsg` or `--show-ref-allele`, which VEP's does not include either). With `--cache-dir` it also turns on `--af --af-1kg --af-gnomade --af-gnomadg --max-af --pubmed`, as VEP's does | off |
| `--distance` | Upstream/downstream distance in bp | `5000` |
| `--buffer-size` | Accepted for VEP compatibility and ignored: records are annotated in batches of 1,024 | ignored |
| `--sa-dir` | Directory containing `.osa` / `.osa2` / `.osi` / `.oga` supplementary annotations | -- |
| `--sa-only` | Skip the default CSQ annotation and emit only supplementary annotations from `--sa-dir` (requires `--sa-dir`) | off |
| `--cache-dir` | VEP's tabix variation cache, for `Existing_variation`: a directory of `{chr}/all_vars.gz` files, each with a `.csi` or `.tbi` index (what `convert_cache.pl` and the indexed-cache download provide). Only known variants come from it: gene models come from `--gff3`, and VEP's `{start}-{end}.gz` transcript slices are never read. Naming a cache is VEP's `--check_existing`: `Existing_variation`, `CLIN_SIG`, `SOMATIC` and `PHENO`; the frequencies and the PubMed ids need `--af*`, `--max-af` and `--pubmed`, as in VEP. Matching follows VEP's (anchored and minimal indels match, as in `--check_existing`). A directory with no variation file, or a file with no index, stops the run. Contigs are matched through a `chr_synonyms.txt` in the directory and the `chr`/`MT` naming rules | -- |
| `--transcript-cache` | Path to binary transcript cache file (overrides the auto-managed `<gff3>.fastvep.cache` sidecar) | -- |
| `--acmg` | Run ACMG-AMP classification (Richards 2015 + ClinGen SVI); adds `ACMG` + `ACMG_CRITERIA` to CSQ, and a last `ACMG_DISORDERS` column with the call against each mode of inheritance of a gene that has more than one | off |
| `--acmg-config` | TOML file with custom ACMG thresholds | built-in defaults |
| `--functional-evidence` | Curated functional-evidence TSV supplying PS3/BS3 (columns: chrom, pos, ref, alt, criterion, strength, pmid, note); an entry also suppresses PP3/BP4/BP7 for that variant | -- |
| `--proband` / `--mother` / `--father` | Sample names for trio analysis — enables PS2 (de novo), PM6, BP2 | -- |
| `--gene-list` | Path to a gene-panel file (one symbol or Ensembl gene ID per line). Tab output drops rows whose transcript isn't on the panel. | -- |
| `--explicit-alleles` | Add an explicit `REF` column to tab output after `Allele` | off |
| `--qc-rules` | TOML file of QC class rules; populates a `QC_CLASS` column in tab output | -- |
| `--no-progress` | Suppress periodic progress output | off |

### `fastvep sa-build`

| Flag | Description | Default |
|------|-------------|---------|
| `--source` | Source type (clinvar, gnomad, dbsnp, cosmic, onekg, topmed, mitomap, phylop, gerp, dann, revel, spliceai, primateai, dbnsfp, alphamissense, omim, gnomad_genes, clinvar_protein, custom_vcf, custom_bed, custom) | *required* |
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
| `--chunk-bits` | Genomic width of one chunk as a power of two (1 to 20). The default of 20 (1 Mb) suits sparse sources (dbSNP, ClinVar, COSMIC, REVEL); a source that scores nearly every base, SpliceAI above all, needs 16 (65 kb) | `20` |
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
| `--data-dir` | Directory of genome subdirectories, one per genome, that a request can name | -- |
| `--preload-genomes` | Load every genome of `--data-dir` in the background at startup instead of on its first request (also `FASTVEP_PRELOAD_GENOMES` env) | off |
| `--port` | HTTP port (also `FASTVEP_PORT` env) | `8080` |
| `--bind` | Bind address (also `FASTVEP_BIND` env) | `0.0.0.0` |
| `--distance` | Upstream/downstream distance in bp | `5000` |
| `--max-body-size` | Max request body in bytes | `10485760` |
| `--max-concurrent` | Max concurrent requests of any kind | `64` |
| `--max-parallel` | Annotations that run at once; the next request waits `--queue-timeout-seconds` for a place, then gets `503` with `Retry-After` (also `FASTVEP_MAX_PARALLEL` env) | number of cores |
| `--queue-timeout-seconds` | How long a request waits for a place to annotate (also `FASTVEP_QUEUE_TIMEOUT_SECONDS` env) | `30` |
| `--max-response-bytes` | Largest JSON answer one request may produce, built in memory a batch of variants at a time; beyond it `413` (also `FASTVEP_MAX_RESPONSE_BYTES` env) | `134217728` |
| `--max-request-seconds` | Longest one annotation may run; beyond it `413` (also `FASTVEP_MAX_REQUEST_SECONDS` env) | `90` |
| `--stats-file` | JSON file in which the running totals of variants annotated and genomes loaded are kept (also `FASTVEP_STATS_FILE` env) | `stats.json` |
| `--allow-model-replacement` | Let `/api/upload-gff3` and `/api/load-genome` replace the gene model of the running server, for every client; choosing a genome per request does not need it (also `FASTVEP_ALLOW_MODEL_REPLACEMENT` env) | off |
| `--acmg-allow-partial` | Answer a request for ACMG classification for a genome whose database directory lacks part of the nine-source reference stack, saying so in every response, instead of refusing it with `422` (also `FASTVEP_ACMG_ALLOW_PARTIAL` env) | off |

### `fastvep cache`

Pre-builds a binary transcript cache for fast startup. Accepts the same
multi-`--gff3` / `LABEL=path` syntax as `annotate`, so a merged Ensembl
+ RefSeq cache can be built once and reused via `--transcript-cache`.

The cache is a sequence of zstd frames of about 8 MB, each with a checksum, and a table of them at
the end. A run with more than two threads decodes the frames side by side; a cache of one frame,
which earlier builds wrote and which has no table, is read as a stream, as is every cache by a run
of one or two threads (which does not look at the table). A cache that fails a checksum, that does
not end where its list does, or whose table does not describe it, is rebuilt from the GFF3 when it
is the auto-managed sidecar, and stops the run when it is named with `--transcript-cache`.

| Flag | Description | Default |
|------|-------------|---------|
| `--gff3` | GFF3 annotation file(s). Repeatable; each value may be `LABEL=path`. | *required* |
| `--fasta` | Reference FASTA (for pre-building spliced sequences) | -- |
| `--synonyms` | VEP-style `chr_synonyms.txt` mapping equivalent contig names (e.g. RefSeq `NC_000017.11` ↔ `17`). Reconciles mixed Ensembl/RefSeq naming against the FASTA. Only effective with `--fasta`. | -- |
| `-o, --output` | Output cache file path | *required* |

## Output Formats

### VCF Output

Consequence annotations are added as a `CSQ` field in the INFO column with 49 pipe-delimited fields matching Ensembl VEP's extended format, of which the last two, fastVEP's `ACMG` and `ACMG_CRITERIA`, are filled when `--acmg` is set (a run with `--acmg` also ends with an `ACMG_DISORDERS` field, empty unless the gene has more than one mode of inheritance; see [docs/ACMG.md](docs/ACMG.md#genes-with-more-than-one-mode-of-inheritance)). When supplementary annotation databases are loaded with `--sa-dir`, fastVEP also emits VCF-compatible INFO projections for supported fastSA sources: standard `SpliceAI` for SpliceAI databases, and fastVEP-specific `FV_*` fields such as `FV_CLINVAR`, `FV_GNOMAD`, `FV_DBSNP`, `FV_REVEL`, and gene-level `FV_OMIM`.

The VCF output never embeds raw JSON in INFO values. Use `--output-format json` for the richest structured representation of all supplementary annotation objects.

### Tab Output

One line per variant, transcript and allele.
The 17 columns are VEP's 13 base columns (`Uploaded_variation`, `Location`, `Allele`, `Gene`, `Feature`, `Feature_type`, `Consequence`, `cDNA_position`, `CDS_position`, `Protein_position`, `Amino_acids`, `Codons`, `Existing_variation`) and `IMPACT`, `DISTANCE`, `STRAND`, `FLAGS`.
Around them go `REF` (with `--explicit-alleles`, right after `Allele`), `PICK` (with `--flag-pick*`, right after `FLAGS`, as in VEP), one column per loaded supplementary source and `QC_CLASS` (with `--qc-rules`), the last two at the end.

The cells follow VEP's tab writer:

- A value a row does not have is `-`, never an empty cell.
- `Location` runs from the lower coordinate to the higher, so an insertion, which sits between two bases, is `chr10:11950-11951`.
- `Uploaded_variation` is the VCF ID (the first, if there are several), or the line itself, `{CHROM}_{POS}_{REF}/{ALT}`, when the record has none: `chr10_10414_C/A`, `chr10_10400_AACCCT/A`, and `chr10_109_CA/C/TA` for a record with two alts.
  It is not built from the minimal alleles of the annotation, as it was.
- `Amino_acids` is `L/P`, or just `L` when the protein does not change, as in the VCF `CSQ` and in JSON.
- A record with one ALT that is an indel is minimised from both ends before it is annotated, as VEP does: `CAGT>CT` is the deletion of `AG` and its `Allele` is `-`, whatever the caller left in.
- An ALT of `*` (a spanning deletion) or one of the placeholders a gVCF carries (`.`, `<NON_REF>`, `<*>`) gets no row, and a record with none but those has none; the real ALTs beside one are annotated as they would be alone.
  An ALT or REF that holds an ambiguity code (`N`) cannot be translated: in the coding sequence it is `coding_sequence_variant` (or `frameshift_variant` when its length moves the frame), with no `Amino_acids`, `HGVSp` or `HGVSc`, as in VEP.
- `cDNA_position`, `CDS_position` and `Protein_position` are the positions of the first and the last base of the change in transcript order, and an end that is not in an exon (for the cDNA) or not in the coding sequence (for the CDS and the protein) has no coordinate: `?-7` for a deletion that starts in an intron or upstream of the transcript, `601-?` for one that runs out of the exon.
  An insertion has the two positions it sits between (`601-602`), also when one of them is across an exon boundary, and no CDS or protein position when it sits just before the initiator or just after the terminator.
  `EXON` and `INTRON` list every exon and intron the change overlaps, `2-3/10` for a change over two, and an insertion at an exon's edge overlaps neither.
  The same cells are written to the VCF `CSQ`, and JSON has `cdna_start` and `cdna_end` (and the CDS and protein pairs) for the ends that have one.
  With `--total-length` each cell is followed by the length it is counted in, as in VEP: `451-?/1403`, `51/654`, `17/217` (the cDNA, the coding sequence with its stop codon, the protein without it).
  The cDNA length goes on every row whose change touches an exon, and the other two on the rows that touch the coding sequence too, whatever the cell holds, so a deletion that covers an exon from intron to intron is `-/1403`, `-/654` and `-/217`.
- A variant with no overlapping transcript has no feature, so `Gene`, `Feature`, `Feature_type`, `DISTANCE`, `STRAND` and `FLAGS` are `-`; `IMPACT` is still `MODIFIER`.
  A VCF `CSQ` entry for it leaves those fields empty.

VEP adds further tab columns only when asked, and so does fastVEP, with VEP's flag names (spelled with hyphens here): `--show-ref-allele`, `--variant-class`, `--symbol`, `--biotype`, `--canonical`, `--numbers`, `--hgvs`, `--hgvsg` and the known-variant flags below.
They are written in the order VEP's own table gives them: `REF_ALLELE` between `Existing_variation` and `IMPACT`, then, after `FLAGS` (and `PICK`), `VARIANT_CLASS`, `SYMBOL`, `SYMBOL_SOURCE`, `HGNC_ID`, `BIOTYPE`, `CANONICAL`, `EXON`, `INTRON`, `HGVSc`, `HGVSp`, `HGVS_OFFSET` and `HGVSg`, and before any supplementary-source column.
A flag you do not give adds nothing, so a run without any of them has the 17 columns above.
`FLAGS` holds the transcript's flags (`cds_start_NF`), as in VEP, and `CANONICAL` says whether the transcript is canonical.
`HGVSp` writes `=` as `%3D`, as VEP does in its tab file and in its VCF (`--no-escape` leaves it an `=`), and `HGVS_OFFSET` is how far the HGVS description moved 3' from where the VCF had the variant, signed by the strand.
The underscore spellings VEP uses for those three flags (`--total_length`, `--no_escape`, `--force_overwrite`) are accepted beside the hyphenated ones.

The columns of the known variants (`--af`, `--af-1kg`, `--af-gnomade`, `--af-gnomadg`, `--max-af`, `--check-existing`, `--pubmed`) read the variation cache of `--cache-dir` and are worked out as VEP works them out, for each input allele and from every known variant of the site.
They come after `HGVSg`, in VEP's order: `AF` and the 1000 Genomes populations, the gnomAD exomes and genomes populations, `MAX_AF`, `MAX_AF_POPS`, `CLIN_SIG`, `SOMATIC`, `PHENO`, `PUBMED`.
A flag that reads the cache without one is an error rather than a column of `-`.
The VCF `CSQ` keeps `AF`, `CLIN_SIG`, `SOMATIC`, `PHENO` and `PUBMED` among its 49 default fields, so that a reader that finds a field by its position keeps finding it, and the columns of the other flags follow `AF`.
They hold VEP's per-allele value or nothing: `CLIN_SIG`, `SOMATIC` and `PHENO` whenever there is a `--cache-dir`, `AF` under `--af` or `--af-1kg`, `PUBMED` under `--pubmed`.
They used to hold a summary of the whole variant (the first frequency of any known variant of the site, gnomAD's where it had one, which is not VEP's `AF`), and for a record with several ALTs that put one allele's number on every allele's row.
A reader of that `AF` now gives `--af-1kg` (VEP's `AF`) or `--af-gnomade`.
How the two tools differ, and what was checked against VEP 116.3, is in [docs/VEP_DIVERGENCE.md](docs/VEP_DIVERGENCE.md).

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
  fastvep-cache/        # GFF3 parser, FASTA reader, annotation providers, regulatory regions,
                       #   the binary transcript cache (zstd frames, decoded side by side)
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
  fastvep-annotate/     # Shared annotation engine (AnnotationContext, used by both web servers; the
                       #   `fastvep annotate` loop is in fastvep-cli and calls the same helpers):
                       #   variant overlap, consequence prediction, HGVS, SA/gene SA
                       #   provider loading, and the pick (the --pick criteria)
  fastvep-classification/ # ACMG-AMP variant classification engine (Richards 2015 +
                       #   ClinGen SVI). 28 criteria, trio/compound-het support,
                       #   configurable thresholds via TOML
  fastvep-cli/          # CLI binary (mimalloc is its allocator). `src/pipeline/`
                       #   holds one module per subcommand: annotate, cache_build,
                       #   sa_build, filter, custom
  fastvep-web/          # Production web server (axum/tokio): async, multi-connection,
                       #   genome switching, SA integration, a cap on concurrent requests
web/                   # Web GUI (HTML/CSS/JS, embedded in both server binaries)
tests/                 # Test data: chr1 (OR4F5) and chr17 (BRCA1) VCF + GFF3
```

## Running Tests

```bash
cargo test --workspace          # 1,731 tests on macOS, 1,732 on Linux, all must pass
cargo test -p fastvep-consequence  # Consequence prediction tests (incl. SV)
cargo test -p fastvep-filter       # Filter engine tests
cargo test -p fastvep-sa           # Supplementary annotation format tests
```

## Performance Benchmarks

fastVEP's numbers below were measured on 2026-11-07 on an Apple M4 (10 cores, 4 performance and 6 efficiency, ARM64, macOS 26.5.1), with the release build (LTO) of 0.5.0.
`benchmarks/results/v0.5.0/` holds every run, what it ran on and the hash of the binary; `benchmarks/run_head_to_head.py` and `benchmarks/run_benchmark.sh` are the commands that made them.
Each is the median of 5 runs, with the full Ensembl annotations, the FASTA and HGVS.
The comparison with Ensembl VEP 116.3 was made on the same machine, with VEP in Docker (`benchmarks/run_vep_side.py`), and the fastVEP times that it is set against were measured in the same session (`benchmarks/results/v0.5.0/with_vep/`): a laptop runs one binary at speeds that differ by several per cent from one session to another, so each speedup is the time of VEP over the time of fastVEP of one session.
A laptop is a noisy instrument: the 5 consequence-only whole-genome runs on one thread took 41 to 46 seconds in that session (43 to 55 in the earlier one), and those on every core 11.0 to 11.3 (12 to 17).

### Multi-Organism Throughput

Default threads, `--hgvs`, the output written to a file, the binary transcript cache already built (the first session of 2026-11-07).

| Organism | Transcripts | Variants | Source | Time | Throughput |
|----------|-------------|----------|--------|------|------------|
| Yeast (R64, full genome) | 7,024 | 260,526 | Ensembl/SGD | 0.5s | **553,839 v/s** |
| Drosophila (BDGP6, full) | 41,600 | 4,438,427 | DGRP2 (dm3 coordinates) | 9.6s | **461,567 v/s** |
| Arabidopsis (TAIR10, full) | 54,013 | 12,883,854 | Ensembl Plants 113 | 27.0s | **477,136 v/s** |
| Mouse (GRCm39, full genome) | 278,396 | 26,062,054 | MGP CAST/EiJ (GRCm38 coordinates) | 65.4s | **398,455 v/s** |
| Human full WGS (GRCh38) | 509,650 | 4,048,342 | GIAB HG002 | 12.4s | **326,930 v/s** |

Transcripts are those that fastVEP reads from each GFF3.
The Drosophila and mouse variant files are in the coordinates of dm3 and GRCm38, against gene models of BDGP6.54 and GRCm39: those two rows measure throughput, and the consequences they predict are not those of the variants.

### vs. Ensembl VEP 116.3 (head-to-head, single-thread, GIAB HG002)

Both tools single-threaded on identical GRCh38 Ensembl 115 GFF3 + FASTA with `--hgvs --symbol --canonical`.
fastVEP uses its binary transcript cache; Ensembl VEP runs via Docker (`ensemblorg/ensembl-vep:release_116.3`, pinned by its digest) in `--gff` mode against a bgzip+tabix-indexed GFF3.
VEP's full-genome time is measured per chromosome and summed, one run at a time: a single-process whole-genome run of VEP 115.1 ran out of the 7.65 GiB Docker default above ~2.7M variants (July 2026, not repeated), and the per-chromosome sum keeps VEP's memory bounded (its best case).
Release 116.3 takes 1.3 to 1.9 times as long as the release 115.1 of the earlier version of this table: on the chr22 inputs below, on the same day, 115.1 took 0.46, 5.98, 12.88, 83.65 s where 116.3 took 0.64, 7.91, 24.07, 112.51 s (`benchmarks/results/v0.5.0/with_vep/vep_115_1/`).

**Full human WGS (4,048,342 variants), single-thread:**

| Annotation | fastVEP | Ensembl VEP | Speedup |
|-----------|---------|-------------|---------|
| Consequence only | 42.2s | 6,516s (~109 min) | **154×** |
| + ClinVar | 44.5s | 6,662s | **150×** |
| + gnomAD + ClinVar | 54.6s | 7,334s | **134×** |

With fastVEP's default multi-threading (10 cores), the same full WGS completes in **11.1s** (consequence-only) to **17.0s** (+gnomAD+ClinVar) on this system.
A whole-genome run on all cores writes 11.1 GB (the multi-organism table above has 12.4s for the same genome with `--hgvs` alone, measured in the first session).
Both tools read the same ClinVar and gnomAD files in the two annotation rows (`data/benchmark/sa_sources`: the ClinVar VCF of 27 June 2026, and the extracts of the joint release of gnomAD v4.1 restricted to the regions of its ClinVar 2-star+ set), fastVEP as the fastSA databases built from them (`data/benchmark/sa_db` and `sa_db_joint`) and VEP as `--custom` sources; on chromosome 22 they found the same 1,057 ClinVar records and 2,623 gnomAD records (`benchmarks/check_custom_matches.py`).

**chr22 scaling, single-thread** — fastVEP reads its binary transcript cache of 509,650 transcripts on every run: a run of one variant takes about 0.9s on one thread and 0.4s on eight, nearly all of it reading the cache.
VEP (which tabix-fetches only overlapping GFF regions) is therefore faster at 1,000 variants, and fastVEP from a few thousand on:

| Variants | fastVEP | VEP | Speedup |
|----------|---------|-----|---------|
| 1,000 | 0.88s | 0.64s | 0.7× (VEP faster) |
| 5,000 | 0.94s | 7.91s | 8.4× |
| 10,000 | 1.05s | 24.07s | 23× |
| 50,284 | 1.54s | 112.51s | **73×** |

| Resource | Ensembl VEP | fastVEP |
|----------|-------------|---------|
| Peak memory (100K variants) | 0.94 GiB (the peak of the VEP process) | 2.6 GiB (the whole gene model, with its sequences, is in memory) |
| Binary / install size | 0.97 GB (the Docker image, unpacked; 0.31 GB to pull) | **5.8 MB** |
| Dependencies | Perl (5.22 or later recommended), the Perl libraries Archive::Zip and DBI, a C compiler, and what its `INSTALL.pl` fetches (the README of ensembl-vep release/115.2) | **None** at run time (a C compiler to build) |

### Threads and memory

fastVEP uses every core unless `RAYON_NUM_THREADS=N` says otherwise.
The peak memory of a run is the gene model and not the VCF or the threads: 2,700 MiB on one thread and 2,702 MiB on eight for 100,000 variants against the human Ensembl GFF3 with sequences (macOS), and 1,170 MiB and 1,188 MiB with a cache built without a FASTA, which has no sequences (consequence only, no `--hgvs`).
On Linux (arm64, glibc 2.39, in a container) a run of 200 variants peaks at 2,648 MB on 4 threads and 2,677 MB on 64, where the build before #166 peaked at 2,520 MB on either, and a run of 100,000 variants peaks at 2,673 MB on 4 threads, 2,720 MB on 32 and 2,834 MB on 64, against 2,519 MB: the peak does not grow with the number of variants, and the `fastvep` binary's allocator (mimalloc) holds 4 to 9% more of it than glibc's does on one to 32 threads and 12% more on 64, and nothing more on macOS.
Under an address-space limit (`ulimit -v`, the `h_vmem` of some schedulers) a run of 100,000 variants on 8 threads needs about 1.2 GiB more than before: it ran down to a limit between 7,000 and 7,200 MiB, against 5,800 to 6,000 MiB, and below that it stops with `Cannot allocate memory`.

### Supplementary annotation format (fastSA v2)

The `.osa2` format (echtvar-inspired chunked binary; `--format auto` default) stores real annotation databases compactly, with annotation output **byte-identical** to the v1 `.osa` format (md5-verified).
Measured on 2026-11-07 with 0.5.0 by `benchmarks/run_sa_formats.py` (sizes in megabytes of 10^6 bytes; the `.osa` with its `.osa.idx`):

| Source | Records | v1 (.osa) | v2 (.osa2) | Ratio | Build v1 | Build v2 |
|--------|---------|-----------|-----------|-------|----------|----------|
| ClinVar | 4.44M | 50.5 MB | 40.9 MB | 0.81× | 10.9s | 12.2s |
| REVEL (chr1) | 8.25M | 42.5 MB | 18.8 MB | 0.44× | 3.5s | 5.4s |
| REVEL (chr22) | 1.78M | 9.1 MB | 4.0 MB | 0.44× | 0.8s | 1.2s |

Annotating the 673,660 variants of the ClinVar 2-star+ set against ClinVar alone takes 7.1 s from the `.osa` and 7.0 s from the `.osa2`, and against ClinVar with both REVEL chromosomes 7.6 s and 7.1 s (median of 3 runs); the four outputs of each pair are md5-identical, and 38,537 of their records carry a REVEL score.
The sizes of SpliceAI, gnomAD and PhyloP databases, and the annotation against five databases at once, are in [`docs/SUPPLEMENTARY_ANNOTATIONS.md`](docs/SUPPLEMENTARY_ANNOTATIONS.md); their sources are not among the files that the repository's scripts fetch, and they were not measured again for this release.

## VEP Concordance

fastVEP's consequence and HGVS output is a port of Ensembl VEP's own model, and it is compared with **real Ensembl VEP 116.3**, the release that most users run now.
VEP ran in Docker (`ensemblorg/ensembl-vep:release_116.3`) in `--gff` mode over the same GRCh38 Ensembl 115 GFF3 and FASTA that fastVEP reads, with `--tab --symbol --biotype --canonical --numbers --variant_class --show_ref_allele --hgvs --hgvsg`, and `validation/compare_tab.py` pairs the rows of the two outputs and compares them cell by cell (27 columns with the HGVS columns and 23 without them: a GFF3 has no symbol source, HGNC identifier or CDS flags, so VEP's `--gff` reader writes `-` where fastVEP, like VEP's own cache, fills them).
The counts below were measured on 7 November 2026 with the binary of 0.5.0 (`validation/measure_vep116.py` makes every row of the table in one command), and are the counts that commit `5048fc5` gave on 27 October 2026, to the cell.
The two builds write the same bytes on every input of the comparison that can be run without VEP, as a VCF and as the tab output of this comparison (11 runs, the whole HG002 call set among them; `validation/compare_builds.py`, [`benchmarks/results/v0.5.0/identity_vs_5048fc5/`](benchmarks/results/v0.5.0/identity_vs_5048fc5/)).

| Input | Records | Paired rows | Cells compared | Cells that differ |
|---|---:|---:|---:|---:|
| GIAB HG002 WGS, a systematic 1-in-200 sample | 20,241 | 186,055 | 5,023,485 | 3,851 (0.077 %) |
| ClinVar 2-star+, a systematic 1-in-5 sample, without the HGVS columns | 134,732 | 3,136,295 | 72,134,785 | 112 (0.00016 %) |
| The same sample with `--hgvs --hgvsg` | 134,732 | 3,136,295 | 84,679,965 | 2,999 (0.0035 %) |
| A generated world built to be hard (both strands, multi-ALT records) | 3,000 | 7,213 | 216,390 | 954 (0.44 %) |
| Its edges: every shape of change within ten bases of every exon, CDS and transcript boundary | 17,672 | 36,958 | 1,108,740 | 845 (0.076 %) |
| Its lone indels, written with shared bases | 1,500 | 3,176 | 95,280 | 26 (0.027 %) |
| The 173 variants of VEP's own example set | 173 | 2,977 | 80,379 | 0 |
| The 400 ClinVar in-frame deletions | 400 | 9,465 | 255,555 | 53 (0.021 %) |

**An ordinary callset.**
On the HG002 sample `Consequence`, `IMPACT`, `Amino_acids`, `Codons`, `cDNA_position`, `CDS_position`, `Protein_position`, `Location`, `REF_ALLELE`, `VARIANT_CLASS`, `DISTANCE`, `STRAND`, `EXON`, `INTRON`, `SYMBOL` and `BIOTYPE` agree on every paired row.
The 3,851 cells that differ are all in four HGVS columns: `HGVS_OFFSET` 2,070, `HGVSc` 1,729, `HGVSg` 35 and `HGVSp` 17.
Each has a documented cause: 3,822 are alleles of records with several ALTs, which fastVEP describes as each would be alone and Ensembl 116 describes by what else the record holds (divergence 11), 12 are the `*` coordinate that Ensembl writes for an intronic base after the terminator (divergence 4), and 17 are a change that the 3'-rule carries over an intron, for which Ensembl writes no protein change (divergence 17).
On a RefSeq GFF3 (chr21 and chr22 of GRCh38.p14, the 527 variants of the sample on them) 818 of 73,359 cells differ: 135 are the multi-allelic difference of divergence 11, and 683 are the canonical flag and the biotype, which a RefSeq GFF3 does not carry (VEP's `--gff` reader guesses the biotype and has no canonical flag, and fastVEP reads them as VEP's RefSeq cache does); against the output of VEP 116.3's own Ensembl cache (`vep --cache --offline`, chromosomes 21 and 22; measured on 27 October 2026 with `5048fc5`, and not repeated, because the extracted caches are no longer on disk) fastVEP over the matching GFF3 differs in one `Consequence` cell of the 168,813 paired rows of the ClinVar sample.

**The shapes that disagree.**
A real callset holds few of them, so the generated worlds are built to hold many, and their rates are the rates on those shapes.
Smaller worlds of the same kinds, with the output that VEP 116.3 wrote for them, are fixtures of the test suite (20,363 records in nine fixtures), and the 26 tests that compare fastVEP's output with them pass.

**HGVSp under in-frame indels.**
For a two-residue deletion in the poly-glutamate C-terminus of NT5C2 (`10:103089678 TTCCTCC>T`) fastVEP writes `p.Glu560_Glu561del`, the 3'-most placement that HGVS asks for, and Ensembl `p.Glu559_Glu560del`: Ensembl's 3'-shift stops short of the end of the protein in two ways, both of which are in release 116 as they were in 115 (divergence 9).

An earlier comparison, against release 115.1 and the build of 15 September 2026 (v0.3.0), is in [docs/VEP_DIVERGENCE.md](docs/VEP_DIVERGENCE.md#what-is-measured); its output differs from this release's on every input, so its numbers are not quoted here.
The comparisons need Docker: `validation/measure_vep116.py` makes the 116.3 one in one command (about an hour on an Apple M4, most of it Ensembl VEP over the ClinVar sample), `validation/run_divergence.sh` makes the 115.1 one again (a few minutes), and Part 3 of that document gives the commands of the 116.3 one one by one.

### Divergences

Agreement is the evidence that the port is faithful, not the objective. The output is a prediction
a clinician may act on, so where Ensembl is demonstrably wrong in a way that changes a call,
fastVEP is right instead.

Against release 115.1 there were eleven such places, and one of them mattered clinically.
Ensembl reported a change that introduces a premature stop as `inframe_insertion,stop_retained_variant`, MODERATE.
BRCA1 `c.5030_5033dup` - ClinVar Pathogenic - is one of 52 such variants in the ClinVar 2-star+ set, 38 of them Pathogenic or Likely pathogenic.
fastVEP reported `stop_gained,frameshift_variant`, HIGH, so those 52 kept PVS1.
Release 116 rewrote that clause, and Ensembl VEP 116.3 writes `stop_gained,frameshift_variant` for them too.
Against 116.3 two of the eleven are gone, that one and the second (a start codon a deletion removes), which fastVEP now reads as VEP does.
Four more (the count or the first residue of a frameshift or of an extension, 5 to 8) turned out not to be differences with VEP but with its `--gff` reader, which takes the 3' UTR of a reverse-strand transcript from the wrong place and which a VEP cache does not do.
Two more were found in the protein notation, one in VEP's 3'-shifting of an insertion on the reverse strand and one in how it reads an unresolved residue (12 and 13).
Release 116.3 itself turned up six more, found by reading the start codon, the UTR, the introns and the VCF the way a clinician's pipeline meets them: a start codon is not lost where the protein still begins with its Met or where the whole coding sequence is behind the edit (14), a change that loses the initiator is `p.Met1?` whatever Ensembl writes (15), a row that is only a UTR variant has no protein change (16), a change that the 3'-rule carries over an intron is described as the change it is in the transcript (17), a record the run cannot place is said so and a REF that is not the genome's is described against the genome (18), and an insertion in front of residue 1 has no `HGVSp` (19).

Every divergence, its cause in the Ensembl source, its row count, and the list of fastVEP's own
remaining gaps in the other direction are in **[docs/VEP_DIVERGENCE.md](docs/VEP_DIVERGENCE.md)**.

Moving a pipeline from `vep`? What fastVEP reads from a VEP installation, a flag-by-flag table, and what is
not covered (plugins, `--custom`, RefSeq BAM edits) are in **[docs/VEP_MIGRATION.md](docs/VEP_MIGRATION.md)**.

## Community

Join the [fastVEP Slack](https://join.slack.com/t/fastvep/shared_invite/zt-4g20tqjvo-wGyIEnN~Ncrc4zXiu_1XJA) for questions and discussion.
Bug reports and feature requests belong in [GitHub issues](https://github.com/Huang-lab/fastVEP/issues).

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
