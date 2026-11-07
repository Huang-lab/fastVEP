# fastVEP benchmark suite

This directory is the entry point for **all** benchmarking in the repo.
There are three full-data benchmarks (plus a synthetic microbench),
each with its own runner, data, and output location. Everything documented
below has been **URL-verified** against the upstream hosts — re-verify
before a fresh download with:

```bash
bash benchmarks/check_urls.sh           # all
bash benchmarks/check_urls.sh speed     # multi-organism throughput only
bash benchmarks/check_urls.sh acmg      # ACMG concordance only
```

| Benchmark | Question it answers | Runner | Inputs | Outputs |
|-----------|---------------------|--------|--------|---------|
| **1. Speed** (multi-organism, and the fastVEP side of the head-to-head with Ensembl VEP) | How fast is annotation across genome sizes, and how much does it cost in memory? | [`benchmarks/run_benchmark.sh`](run_benchmark.sh) (five organisms), [`benchmarks/run_head_to_head.py`](run_head_to_head.py) (chr22 scaling, whole genome, start-up, memory, formats) | `test_data/organisms/<org>/` (~30 GB, not in the repository) | [`benchmarks/results/`](results/), one directory per build measured |
| **2. ACMG concordance** (clinical)   | How well does `--acmg` agree with curated ClinVar 2-star+ calls? | [`analysis/acmg_benchmark/scripts/run_benchmark.sh`](../analysis/acmg_benchmark/scripts/run_benchmark.sh) | `data/benchmark/` (~30 GB SA dbs + sources, not in the repository) | [`analysis/acmg_benchmark/results/`](../analysis/acmg_benchmark/results/), [`manuscript/acmg/data/`](../manuscript/acmg/data/) |
| **3. VEP validation** (correctness)  | Does fastVEP agree with Ensembl VEP on the same input? | [`validation/run_validation.sh`](../validation/run_validation.sh), [`validation/run_divergence.sh`](../validation/run_divergence.sh) | [`validation/human/`](../validation/human/), [`validation/mouse/`](../validation/mouse/), the human inputs of benchmark 1 | the counts of [`docs/VEP_DIVERGENCE.md`](../docs/VEP_DIVERGENCE.md); the runs write `validation/results/` (not in the repository) |

A fourth, narrower benchmark — the synthetic `.osa` reader microbench
([`crates/fastvep-sa/examples/bench_sa_reader.rs`](../crates/fastvep-sa/examples/bench_sa_reader.rs))
— is invoked via `cargo run --release --example bench_sa_reader -p fastvep-sa`
and needs no external data.

---

## Small test data (lives in-repo, 28 KB, gitignored exceptions)

These tiny fixtures are tracked in git so the README quickstart and the docs
walkthroughs work without downloading anything:

| Path | Contents | Used by |
|------|----------|---------|
| [`tests/test.vcf`](../tests/test.vcf)               | 8 variants on chr17, 7 in BRCA1 and 1 in TP53 (coding, 5' UTR, 3' UTR, intronic, downstream, one deletion) | README quickstart, [`docs/ACMG_SETUP.md`](../docs/ACMG_SETUP.md), manual smoke runs |
| [`tests/test.gff3`](../tests/test.gff3)             | BRCA1 and TP53, 49 transcripts                      | README quickstart, [`docs/ACMG_SETUP.md`](../docs/ACMG_SETUP.md), manual smoke runs |
| [`tests/test_chr1.vcf`](../tests/test_chr1.vcf)     | 11 chr1 variants                                    | manual smoke runs for chr1-specific code paths |
| [`tests/test_chr1.gff3`](../tests/test_chr1.gff3)   | chr1 transcript fragment                            | manual smoke runs for chr1-specific code paths |

The `tests/` directory itself is gitignored (large generated artefacts can
land there during a run); these specific files are force-added with
`git add -f`. **Do not delete them** without checking what depends on them.

Note: the Rust test suite does *not* read from `tests/` — every `cargo test`
fixture is either inline in the test source or under
`crates/*/fixtures/`, so `cargo test` passes on a clone that lacks this
directory entirely.

---

## 1. Speed benchmark — multi-organism throughput

**What it does.** Runs `fastvep annotate --hgvs` end-to-end on a real
genome-scale VCF for each of yeast, drosophila, arabidopsis, mouse, and
human (3 runs, median reported, binary GFF3 cache pre-warmed), with every
core and the output written to a file. A run in which `fastvep` exits with
an error stops the script: it is not timed. `FASTVEP_BIN=/path/to/fastvep`
benchmarks another binary (a build of a branch, or a release archive) and
skips the build, and `RUNS=5` takes the median of five.
The cache of each gene model is warmed with a run of one variant before its timed runs, also when a cache is already there:
a cache that an older version wrote is rejected and rebuilt by the first run that meets it, and that run must not be a timed one.
The script writes `benchmark_environment.tsv` beside its CSV (date, system, compiler, the binary's version and hash, runs).

**Run it:**

```bash
bash benchmarks/download_data.sh --all     # ~12 GB raw download, ~30 GB decompressed
bash benchmarks/run_benchmark.sh           # writes benchmarks/output/benchmark_results.csv
```

`benchmarks/output/` is not in the repository.
The measurements that the README and the manuscript quote are copied from it into `benchmarks/results/<build>/`, which is.

Or one organism at a time: `bash benchmarks/download_data.sh yeast` etc.

### Real benchmark inputs (sources, URLs, sizes)

All URLs verified working. `download_data.sh` writes each file to the
location `run_benchmark.sh` expects.

| Organism | Assembly | Reference (GFF3 + FASTA) — Ensembl 115 | Benchmark VCF | Variants | Source URL |
|----------|----------|----------------------------------------|---------------|---------:|------------|
| Yeast       | R64-1-1   | [Ensembl 115 yeast](https://ftp.ensembl.org/pub/release-115/gff3/saccharomyces_cerevisiae/)        | `yeast_ensembl_full.vcf`         | ~260 K  | [Ensembl 115 saccharomyces_cerevisiae.vcf.gz](https://ftp.ensembl.org/pub/release-115/variation/vcf/saccharomyces_cerevisiae/saccharomyces_cerevisiae.vcf.gz) |
| Drosophila  | BDGP6.54 (dm6) | [Ensembl 115 drosophila](https://ftp.ensembl.org/pub/release-115/gff3/drosophila_melanogaster/) | `drosophila_dgrp2_full.vcf`      | ~4.4 M  | [DGRP2 dm6-lifted (Zenodo)](https://zenodo.org/record/155396/files/dgrp2_dm6_dbSNP.vcf.gz) — Aerts Lab deposit (Zenodo DOI [10.5281/zenodo.155396](https://doi.org/10.5281/zenodo.155396)). The original NCSU host (`dgrp2.gnets.ncsu.edu`) is offline. |
| Arabidopsis | TAIR10    | [Ensembl Plants 60 arabidopsis](http://ftp.ensemblgenomes.org/pub/plants/release-60/gff3/arabidopsis_thaliana/) | `arabidopsis_1001g_full.vcf`     | ~12.9 M | [Ensembl Plants 60 arabidopsis_thaliana.vcf.gz](http://ftp.ensemblgenomes.org/pub/plants/release-60/variation/vcf/arabidopsis_thaliana/arabidopsis_thaliana.vcf.gz) (TAIR10, working equivalent to the original 1001 Genomes VCF — see note below) |
| Mouse       | GRCm39    | [Ensembl 115 mouse](https://ftp.ensembl.org/pub/release-115/gff3/mus_musculus/)                    | `mouse_ensembl_1m.vcf` (head -n 1M of MGP REL-2112) | 1 M | [MGP REL-2112 mgp_REL2021_snps.vcf.gz](https://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-2112-v8-SNPs_Indels/mgp_REL2021_snps.vcf.gz) |
| Human       | GRCh38    | [Ensembl 115 human](https://ftp.ensembl.org/pub/release-115/gff3/homo_sapiens/)                    | `human_giab_hg002_full.vcf` (+ high-confidence BED) | ~4.05 M | [GIAB HG002 v4.2.1 benchmark VCF](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz) + [`_noinconsistent.bed`](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed) |

**Optional sixth organism (opt in with `download_data.sh elegans`)**: C. elegans
WBcel235 + [CaeNDR WI.20250625 hard-filter isotype VCF](https://caendr-open-access-data-bucket.s3.us-east-2.amazonaws.com/dataset_release/c_elegans/20250625/variation/WI.20250625.hard-filter.isotype.vcf.gz)
(~11 GB; set `CAENDR_SUBSET=1` to subset to 1 M variants instead).

### Substitutions made vs the manuscript's original sources

| Dataset                       | Original URL (broken)                                              | What we use now |
|-------------------------------|--------------------------------------------------------------------|-----------------|
| **DGRP2** freeze-2 drosophila (~4.4M) | `http://dgrp2.gnets.ncsu.edu/data/website/dgrp2.vcf.gz` (timeout) | **Resolved.** [Aerts Lab Zenodo deposit](https://zenodo.org/records/155396) hosts the dm6-lifted DGRP2 + dbSNP IDs (`dgrp2_dm6_dbSNP.vcf.gz`, 173 MB). dm6 == BDGP6, so it matches the Ensembl reference. |
| **MGP CAST/EiJ** mouse (~26M) | `https://ftp.ebi.ac.uk/.../CAST_EiJ.mgp.v8.snps.dbSNP142.vcf.gz` (404), and `ftp.mousegenomes.org/pub/release/current_snps/` unreachable | Aggregated [MGP REL-2112 SNPs on EBI](https://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-2112-v8-SNPs_Indels/mgp_REL2021_snps.vcf.gz) — fewer variants than the strain-specific file but on the official EBI mirror. `head -n 1,000,000` for the `mouse_ensembl_1m.vcf` slot. |
| **1001 Genomes** Arabidopsis (~12.9M) | `1001genomes.org/data/GMI-MPI/releases/v3.1/SNPs_all_methods/*.vcf.gz` (404; site is now a SPA) | The `intersection_snp_short_indel_vcf/` directory contains only per-strain VCFs, not a combined VCF. We use Ensembl Plants 60's `arabidopsis_thaliana.vcf.gz` instead — same TAIR10 coordinates. |

If you need to exactly reproduce the manuscript variant counts in
`benchmarks/output/benchmark_results.csv`,
the preserved on-disk VCFs in `test_data/organisms/{drosophila,mouse}/`
are the source of truth: `run_benchmark.sh` uses `mouse_cast_eij_full.vcf`
(the 26M-variant file the README's table reports) when it is there, and the
1M-variant `mouse_ensembl_1m.vcf` otherwise.

Two of the preserved variant files are not in the assembly of the gene model they are annotated against.
The header of `drosophila_dgrp2_full.vcf` says `##reference=dm3`, against the BDGP6.54 (dm6) gene model, and that of `mouse_cast_eij_full.vcf` names `GRCm38_68.fa`, against the GRCm39 gene model.
The throughput of those two runs is measured all the same, and the consequences they predict are not those of the variants.
The Zenodo file in the table above is lifted to dm6 and would replace the first; a GRCm39 call set of the Mouse Genomes Project would replace the second (neither has been run).
The yeast file is of Ensembl 115, like its gene model, and the human file is GRCh38; the header of the Arabidopsis file says Ensembl release 113 (`##source=ensembl;version=113`) and its gene model is Ensembl Plants 60, and both are of the TAIR10 assembly.

### Head-to-head with Ensembl VEP (the fastVEP side), start-up, memory and formats

The README's comparison with Ensembl VEP v115.1 runs each tool on one thread
over the whole of GIAB HG002 (4,048,342 variants, GRCh38 Ensembl 115 GFF3 and
FASTA) with `--hgvs --symbol --canonical`, and on the first 1,000, 5,000 and
10,000 variants of its chromosome 22 and all 50,284 of them.
[`run_head_to_head.py`](run_head_to_head.py) is the fastVEP side of that, and of the resource table beside it, in one command:

```bash
python3 benchmarks/run_head_to_head.py --bin target/release/fastvep \
    --built-from "v0.5.0-rc.1 (c12d9f9)" --out benchmarks/results/v0.5.0-rc.1
python3 benchmarks/run_head_to_head.py --bin ... --only scaling,memory --runs 3    # a part of it
```

It measures, with the FASTA and the cache already built: the first run on a gene model (parse the GFF3 and write the cache), start-up (a run of one variant),
the chr22 series, the peak memory of 100,000 variants, the three output formats, the whole genome on one thread and on every core
in three annotation conditions (consequences only, with ClinVar, with gnomAD and ClinVar), and the size of the binary.
Every cell is run once per round and the rounds are repeated, so that a slow minute falls on all cells and not on one, and the median is reported.
It writes every run (`raw.tsv`), the medians (`summary.tsv`), what it ran on (`environment.tsv`: the machine, the binary's hash, the SHA-256 of every input, the load)
and whether the multi-threaded outputs are the single-threaded ones byte for byte (`identity.tsv`), and it stops when they are not.
The two annotation conditions read the ClinVar and gnomAD databases of the ACMG benchmark (`data/benchmark/sa_db`, see below); `--sa-db DIR` names others.
Nothing else should run on the machine while it does: a browser, a backup or another build moves the medians by more than the differences measured.

`consequence_distribution.py` counts the consequence calls of an annotated VCF (Table 8 of the manuscript), from a file or a pipe, so that the 10 GB output of a whole genome need not be kept.
`check_lean_inputs.py` measures what the cut of the chromosome 22 inputs to five columns is worth to fastVEP (about 3% of the time of 50,284 records on one thread), because the inputs that Ensembl VEP read in July are not recorded.
`run_sa_formats.py` measures the fastSA v1 against v2 table (sizes, build times, and the annotation of the ClinVar 2-star+ set against each format, which must come out as the same bytes).
`manuscript/build_data.py` turns a results directory into the CSV files of `manuscript/data/`, and checks that they still agree with it; `manuscript/verify_numbers.py` checks that the numbers of the text are those of the CSV files.
`benchmarks/run_vep_side.py` is the Ensembl VEP side (the chromosome 22 series, start-up and peak memory in Docker); it was written for the release candidate and has not been run, because Docker could not be started on the day (`benchmarks/results/v0.5.0-rc.1/README.md`).

Ensembl VEP's side (the Docker image, `--gff` mode against a bgzip and tabix copy of the same
GFF3, the whole genome as the sum of 22 per-chromosome runs because a single
process runs out of the Docker default memory above about 2.7M variants) is
described in the manuscript's methods; its times were measured on
2026-07-29 and are kept in `manuscript/data/vep_baseline_2026-07-29.csv`.

---

## 2. ACMG concordance benchmark — clinical correctness

**What it does.** Annotates every ClinVar 2-star+ variant (~674 K SNV/small indels) with
`fastvep annotate --acmg --pick`, then scores per-class concordance
against the ClinVar review-panel call.
ClinVar's classifications are a reference and not ground truth.
The classifier of this release is v38 ([`analysis/acmg_benchmark/RUN_VERSIONS.md`](../analysis/acmg_benchmark/RUN_VERSIONS.md), Section 3.8 of the manuscript in [`manuscript/acmg/`](../manuscript/acmg/)).
Of the 673,278 variants of the April 2026 ClinVar release that it scores it agrees exactly with 426,022 (63.3 %), in direction with 524,380 (77.9 %),
and calls 48 in the opposite direction (0.007 %): 24 ClinVar pathogenic or likely pathogenic variants called likely benign or benign, and 24 the other way.
The scored matrix, the firing counts and the provenance of every run of the manuscript are in
[`analysis/acmg_benchmark/results/v37_manuscript/`](../analysis/acmg_benchmark/results/v37_manuscript/).

**Run it:**

```bash
# 1. Download SA source files (~1.4 GB)
bash analysis/acmg_benchmark/scripts/sa_sources/download_sa_sources.sh

# 2. Build the supplementary annotation .osa / .oga databases
#    (see scripts in analysis/acmg_benchmark/scripts/sa_sources/ - `build_gnomad_per_chrom.sh`,
#    `build_spliceai_phylop.sh`, and `fastvep sa-build` for the rest).

# 3. Extract the ClinVar 2-star+ classifications
python3 analysis/acmg_benchmark/scripts/01_extract_clinvar_2star.py \
    data/benchmark/sa_sources/clinvar.vcf.gz \
    data/benchmark/

# 4. Annotate + score concordance
bash analysis/acmg_benchmark/scripts/run_benchmark.sh v1

# 5. (optional) Regenerate the manuscript figures
python3 analysis/acmg_benchmark/scripts/generate_figures.py
```

### Real benchmark inputs (sources, URLs)

All URLs verified working. Downloaded by
[`analysis/acmg_benchmark/scripts/sa_sources/download_sa_sources.sh`](../analysis/acmg_benchmark/scripts/sa_sources/download_sa_sources.sh).

| Source | Drives | Size | URL |
|--------|--------|-----:|-----|
| ClinVar GRCh38 VCF                        | `clinvar.osa` (per-allele ClinSig)                                | 190 MB | [ftp.ncbi.nlm.nih.gov/.../clinvar.vcf.gz](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz) |
| ClinVar `variant_summary.txt.gz`          | `clinvar_protein.oga` (PS1/PM1/PM5 protein-position catalog)      | 419 MB | [ftp.ncbi.nlm.nih.gov/.../variant_summary.txt.gz](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz) |
| gnomAD v4.1 exomes per-chrom VCFs         | `gnomad_chrN.osa2` (PM2/BA1/BS1/BS2) — tabix-extracted to ClinVar 2-star+ regions | ~12 GB raw / chrom (only the ClinVar regions are kept) | [storage.googleapis.com/.../gnomad.exomes.v4.1.sites.chr22.vcf.bgz](https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/exomes/gnomad.exomes.v4.1.sites.chr22.vcf.bgz) (one per chrom) |
| gnomAD v4.1 constraint TSV                | `gnomad_genes.oga` (pLI/LOEUF/misZ for PVS1/PP2/BP1)              | 95 MB  | [storage.googleapis.com/.../gnomad.v4.1.constraint_metrics.tsv](https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv) |
| REVEL v1.3 all chromosomes                | per-chrom `revel_chrN.osa` (PP3/BP4 missense)                     | 637 MB | [rothsj06.dmz.hpc.mssm.edu/.../revel-v1.3_all_chromosomes.zip](https://rothsj06.dmz.hpc.mssm.edu/revel-v1.3_all_chromosomes.zip) |
| SpliceAI scores                           | per-chrom `spliceai_chrN.osa` (PP3/BP4/BP7 splice)                | 0 (distilled from gnomAD v4 `INFO/spliceai_ds_max`) | - see [`build_spliceai_phylop.sh`](../analysis/acmg_benchmark/scripts/sa_sources/build_spliceai_phylop.sh) |
| PhyloP (Zoonomia 241-mammal)              | per-chrom `phylop_chrN.osa` (BP7 conservation gate)               | 0 (distilled from gnomAD v4 `INFO/phylop`) | - see [`build_spliceai_phylop.sh`](../analysis/acmg_benchmark/scripts/sa_sources/build_spliceai_phylop.sh) |
| ClinGen Gene-Disease Validity (CSV)       | `omim.oga` (PVS1 disease-gene fallback)                           | 1 MB   | [search.clinicalgenome.org/kb/gene-validity/download](https://search.clinicalgenome.org/kb/gene-validity/download) |

OMIM `genemap2.txt` is **also** accepted (registration-gated at
omim.org); ClinGen GDV is preferred per Abou Tayoun 2018 and is what the
runs of the manuscript use. Both populate the same `omim` json_key in the `.oga` schema.

A fresh `sa-build` writes the v2 format (`.osa2`) for the sources that have it, and the reader opens either.
The databases of the benchmark itself are a mix, made on different days: ClinVar, REVEL, SpliceAI and PhyloP are v1 (`.osa`) and gnomAD is v2 (`.osa2`).
The head-to-head of benchmark 1 reads the same ClinVar and gnomAD files.

### Reading the output

- [`analysis/acmg_benchmark/STATUS.md`](../analysis/acmg_benchmark/STATUS.md) - on-disk inventory (what's downloaded, what's built)
- [`analysis/acmg_benchmark/RUN_VERSIONS.md`](../analysis/acmg_benchmark/RUN_VERSIONS.md) - v1 → v38: what changed between runs and its measured effect, and the record of the reruns
- [`analysis/acmg_benchmark/METHODS.md`](../analysis/acmg_benchmark/METHODS.md) — full methodology & per-criterion implementation
- [`analysis/acmg_benchmark/results/`](../analysis/acmg_benchmark/results/) - the concordance matrix and criterion firing rates of the runs that are kept (v1, v7, v8, v9, v23, v26, v27, v37 and the runs of the manuscript)
- [`manuscript/acmg/`](../manuscript/acmg/) - the manuscript, and its tables and figures with the checks that keep them in step

---

## 3. VEP validation

**What it does.** Runs both fastVEP and Ensembl VEP (via Docker) on the
same VCF + GFF3 + FASTA and diffs the outputs. Lives in
[`validation/`](../validation/) — see
[`validation/run_validation.sh`](../validation/run_validation.sh) (the VEP example variants and a mouse chromosome)
and [`validation/run_divergence.sh`](../validation/run_divergence.sh) (the three inputs whose counts are in
[`docs/VEP_DIVERGENCE.md`](../docs/VEP_DIVERGENCE.md): a ClinVar sample built to be hard, one in 200 of the HG002 variants, and 400 in-frame deletions).
Both compare against VEP 115.1; the comparisons with VEP 116.3 and with VEP's own cache (Part 3 of that document) have their own scripts, `run_tab_world.sh`, `run_cache_oracle.sh` and `run_known_columns.sh`.
Ensembl VEP's output does not depend on the build of fastVEP, so a comparison that was measured with one build describes any build that writes the same bytes:
[`validation/compare_builds.py`](../validation/compare_builds.py) runs two builds over the inputs of these comparisons, as a VCF and as the tab output of Part 3, and compares what they write byte for byte.
That is how the counts of Part 3, measured on 27 October 2026 with commit `5048fc5`, are the release candidate's: `benchmarks/results/v0.5.0-rc.1/identity_vs_5048fc5/` has the hashes of both.

External data needed:
- Mouse chr19 FASTA: [Ensembl 115 chromosome.19.fa.gz](https://ftp.ensembl.org/pub/release-115/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.chromosome.19.fa.gz) — auto-downloaded by the script
- Mouse chr19 GFF3: [Ensembl 115 chromosome.19.gff3.gz](https://ftp.ensembl.org/pub/release-115/gff3/mus_musculus/Mus_musculus.GRCm39.115.chromosome.19.gff3.gz) — auto-downloaded by the script
- Human chr22 1KGP VCF: tracked in-repo at [`validation/human/chr22_1kgp.vcf`](../validation/human/chr22_1kgp.vcf)
- Human VEP example VCF: tracked in-repo at [`validation/human/vep_example_GRCh38.vcf`](../validation/human/vep_example_GRCh38.vcf)
- Docker image: `ensemblorg/ensembl-vep:release_115.1`

---

## 4. SA reader microbench (developer-only)

Synthetic; no external data. Useful for catching regressions in the
`.osa` query path:

```bash
cargo run --release --example bench_sa_reader -p fastvep-sa
# env knobs: SA_BENCH_RECORDS, SA_BENCH_VARIANTS, SA_BENCH_TRANSCRIPTS, SA_BENCH_BATCH
```

Source: [`crates/fastvep-sa/examples/bench_sa_reader.rs`](../crates/fastvep-sa/examples/bench_sa_reader.rs).
