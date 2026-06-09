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
| **1. Speed** (multi-organism)        | How fast is annotation across genome sizes? | [`benchmarks/run_benchmark.sh`](run_benchmark.sh)                                               | [`test_data/organisms/<org>/`](../test_data/organisms/) (~30 GB)                  | [`benchmarks/output/benchmark_results.csv`](output/benchmark_results.csv) |
| **2. ACMG concordance** (clinical)   | How well does `--acmg` agree with curated ClinVar 2-star+ calls? | [`data/benchmark/run_full_benchmark.sh`](../data/benchmark/run_full_benchmark.sh)               | [`data/benchmark/`](../data/benchmark/) (~30 GB SA dbs + sources)                 | [`data/benchmark/output_v7/`](../data/benchmark/output_v7/) |
| **3. VEP validation** (correctness)  | Does fastVEP agree with Ensembl VEP on the same input? | [`validation/run_validation.sh`](../validation/run_validation.sh)                                | [`validation/human/`](../validation/human/), [`validation/mouse/`](../validation/mouse/) | [`validation/results/`](../validation/results/) |

A fourth, narrower benchmark — the synthetic `.osa` reader microbench
([`crates/fastvep-sa/examples/bench_sa_reader.rs`](../crates/fastvep-sa/examples/bench_sa_reader.rs))
— is invoked via `cargo run --release --example bench_sa_reader -p fastvep-sa`
and needs no external data.

---

## Small test data (lives in-repo, ~40 KB, gitignored exceptions)

These tiny fixtures are tracked in git so `cargo test`, the README quickstart,
and CI all work without downloading anything:

| Path | Contents | Used by |
|------|----------|---------|
| [`tests/test.vcf`](../tests/test.vcf)               | 12 BRCA1 variants (SNV/indel/splice/UTR/intergenic) | README quickstart, integration tests |
| [`tests/test.gff3`](../tests/test.gff3)             | BRCA1 region, ~10 transcripts                       | README quickstart, integration tests |
| [`tests/test.fa`](../tests/test.fa)                 | BRCA1 region FASTA                                  | HGVSp + amino-acid prediction tests |
| [`tests/test_chr1.vcf`](../tests/test_chr1.vcf)     | A few chr1 variants                                 | smoke tests for chr1-specific code paths |
| [`tests/test_chr1.gff3`](../tests/test_chr1.gff3)   | chr1 transcript fragment                            | smoke tests for chr1-specific code paths |

The `tests/` directory itself is gitignored (large generated artefacts can
land there during a run); these specific files are force-added with
`git add -f`. **Do not delete them** without checking what depends on them.

---

## 1. Speed benchmark — multi-organism throughput

**What it does.** Runs `fastvep annotate --hgvs` end-to-end on a real
genome-scale VCF for each of yeast, drosophila, arabidopsis, mouse, and
human (3 runs, median reported, binary GFF3 cache pre-warmed).

**Run it:**

```bash
bash benchmarks/download_data.sh --all     # ~12 GB raw download, ~30 GB decompressed
bash benchmarks/run_benchmark.sh           # writes benchmarks/output/benchmark_results.csv
```

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
[`benchmarks/output/benchmark_results.csv`](output/benchmark_results.csv),
the preserved on-disk VCFs in `test_data/organisms/{drosophila,mouse}/`
are the source of truth.

---

## 2. ACMG concordance benchmark — clinical correctness

**What it does.** Annotates every ClinVar 2-star+ variant (~674 K SNV/small indels) with
`fastvep annotate --acmg --pick`, then scores per-class concordance
against the ClinVar review-panel call. The current state of the art is in
[`data/benchmark/output_v7/`](../data/benchmark/output_v7/) — **70.8 % same-direction
concordance, 0.0 % opposite-direction**.

**Run it:**

```bash
# 1. Download SA source files (~1.4 GB)
bash data/benchmark/sa_sources/download_sa_sources.sh

# 2. Build the supplementary annotation .osa / .oga databases
#    (see scripts in data/benchmark/sa_sources/ — `build_gnomad_per_chrom.sh`,
#    `build_spliceai_phylop.sh`, and `fastvep sa-build` for the rest).

# 3. Generate the ClinVar 2-star+ truth set
python3 analysis/acmg_benchmark/real_data/01_extract_clinvar_2star.py \
    data/benchmark/sa_sources/clinvar.vcf.gz \
    data/benchmark/

# 4. Annotate + score concordance
bash data/benchmark/run_full_benchmark.sh

# 5. (optional) Regenerate the manuscript figures
python3 analysis/acmg_benchmark/real_data/generate_figures.py
```

### Real benchmark inputs (sources, URLs)

All URLs verified working. Downloaded by
[`data/benchmark/sa_sources/download_sa_sources.sh`](../data/benchmark/sa_sources/download_sa_sources.sh).

| Source | Drives | Size | URL |
|--------|--------|-----:|-----|
| ClinVar GRCh38 VCF                        | `clinvar.osa` (per-allele ClinSig)                                | 190 MB | [ftp.ncbi.nlm.nih.gov/.../clinvar.vcf.gz](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz) |
| ClinVar `variant_summary.txt.gz`          | `clinvar_protein.oga` (PS1/PM1/PM5 protein-position catalog)      | 419 MB | [ftp.ncbi.nlm.nih.gov/.../variant_summary.txt.gz](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz) |
| gnomAD v4.1 exomes per-chrom VCFs         | `gnomad_chrN.osa` (PM2/BA1/BS1/BS2) — tabix-extracted to ClinVar 2-star+ regions | ~12 GB raw / chrom (only the ClinVar regions are kept) | [storage.googleapis.com/.../gnomad.exomes.v4.1.sites.chr22.vcf.bgz](https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/exomes/gnomad.exomes.v4.1.sites.chr22.vcf.bgz) (one per chrom) |
| gnomAD v4.1 constraint TSV                | `gnomad_genes.oga` (pLI/LOEUF/misZ for PVS1/PP2/BP1)              | 95 MB  | [storage.googleapis.com/.../gnomad.v4.1.constraint_metrics.tsv](https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv) |
| REVEL v1.3 all chromosomes                | per-chrom `revel_chrN.osa` (PP3/BP4 missense)                     | 637 MB | [rothsj06.dmz.hpc.mssm.edu/.../revel-v1.3_all_chromosomes.zip](https://rothsj06.dmz.hpc.mssm.edu/revel-v1.3_all_chromosomes.zip) |
| SpliceAI scores                           | per-chrom `spliceai_chrN.osa` (PP3/BP4/BP7 splice)                | 0 (distilled from gnomAD v4 `INFO/spliceai_ds_max`) | — see [`build_spliceai_phylop.sh`](../data/benchmark/sa_sources/build_spliceai_phylop.sh) |
| PhyloP (Zoonomia 241-mammal)              | per-chrom `phylop_chrN.osa` (BP7 conservation gate)               | 0 (distilled from gnomAD v4 `INFO/phylop`) | — see [`build_spliceai_phylop.sh`](../data/benchmark/sa_sources/build_spliceai_phylop.sh) |
| ClinGen Gene-Disease Validity (CSV)       | `omim.oga` (PVS1 disease-gene fallback)                           | 1 MB   | [search.clinicalgenome.org/kb/gene-validity/download](https://search.clinicalgenome.org/kb/gene-validity/download) |

OMIM `genemap2.txt` is **also** accepted (registration-gated at
omim.org); ClinGen GDV is preferred per Abou Tayoun 2018 and is what the
v7 run used. Both populate the same `omim` json_key in the `.oga` schema.

### Reading the output

- [`data/benchmark/STATUS.md`](../data/benchmark/STATUS.md) — on-disk inventory (what's downloaded, what's built)
- [`data/benchmark/RUN_VERSIONS.md`](../data/benchmark/RUN_VERSIONS.md) — v1 → v7 delta (per-run SA stack + code-fix diff)
- [`analysis/acmg_benchmark/METHODS.md`](../analysis/acmg_benchmark/METHODS.md) — full methodology & per-criterion implementation
- [`data/benchmark/output_v7/`](../data/benchmark/output_v7/) — current concordance matrix, criterion fire rates, figures

---

## 3. VEP validation

**What it does.** Runs both fastVEP and Ensembl VEP (via Docker) on the
same VCF + GFF3 + FASTA and diffs the outputs. Lives in
[`validation/`](../validation/) — see
[`validation/run_validation.sh`](../validation/run_validation.sh).

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
