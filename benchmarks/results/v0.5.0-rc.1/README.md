# Measurements of fastVEP v0.5.0-rc.1

Everything in this directory was measured on 7 November 2026 (UTC) with the binary built from the tag `v0.5.0-rc.1` (commit `c12d9f9`, whose crates are those of `master` at `aa5e5ee`; SHA-256 `0063e01bd143...`, 5,834,640 bytes), on an Apple M4 (10 cores: 4 performance and 6 efficiency; 24 GiB; macOS 26.5.1; rustc 1.98.0).
What the README and `manuscript/manuscript.md` quote about speed, memory, the five organisms, the consequence distribution, the fastSA formats and the tests is in these files, and `manuscript/build_data.py` turns them into `manuscript/data/`.
`environment.tsv` holds the machine, the load, the hash of the binary and of every input.

## What is here

| File | Made by | What it holds |
|------|---------|---------------|
| `raw.tsv`, `summary.tsv`, `environment.tsv`, `identity.tsv`, `gene_models.tsv` | `benchmarks/run_head_to_head.py` | every run of Tables 2 to 5 and of the whole genome (round 0 is the untimed warm-up), their medians, what it ran on, whether the multi-threaded outputs are the single-threaded ones byte for byte, and how many transcripts each gene model holds |
| `multi_organism.csv`, `multi_organism_environment.tsv` | `benchmarks/run_benchmark.sh` (`RUNS=5`, under `caffeinate -i`) | Table 7: the five organisms, every core, `--hgvs`, median of five runs |
| `consequence_most_severe.csv`, `consequence_every_term.csv`, `consequence_run.txt` | `benchmarks/consequence_distribution.py`, over the output of the whole HG002 genome | Table 8: 37,214,646 calls of 4,048,342 variants, each under its most severe term and under every term |
| `code_metrics.csv`, `so_terms.csv`, `crate_dependencies.csv` | `manuscript/code_metrics.py` | the lines, tests, consequence terms and crate dependencies that the manuscript counts |
| `tests_executed.csv` | `cargo test --workspace --no-fail-fast --locked` | what ran: 1,731 passed, 0 failed, 1 ignored, in 121 test targets |
| `sa_formats.csv`, `sa_formats_annotation.csv`, `sa_formats_environment.tsv` | `benchmarks/run_sa_formats.py` | Table 6: fastSA v1 against v2 |
| `lean_inputs.tsv` | `benchmarks/check_lean_inputs.py` | the chromosome 22 inputs cut to five columns, against the records as they are |
| `identity_vs_5048fc5/`, `identity_vs_ac2e2b6/` | `validation/compare_builds.py` | whether this build writes the same bytes as the build of 27 October 2026 (`5048fc5`) and of 15 September 2026 (`ac2e2b6`) on the inputs of the comparisons with Ensembl VEP |

The ACMG runs of the release candidate are in `analysis/acmg_benchmark/results/v38_rc1/` (`RUN_VERSIONS.md` there, "The release candidate reproduces v38").

## Three measurements of one output

The VCF that `annotate --hgvs --symbol --canonical` writes for the whole HG002 genome (4,048,342 records, 11,133,276,060 bytes) has the MD5 `418eb82f968435d459db01609fddf3a5` in all three places that hash it:
the single-threaded and the all-core runs of the head-to-head (`identity.tsv`), the stream that `consequence_distribution.py` counted (`consequence_run.txt`), and the run of the build of 27 October 2026 (`identity_vs_5048fc5/identity.tsv`).

## Runs that were made and are not kept

A laptop is a noisy instrument, and four runs were set aside because something other than the program was running.
Each was repeated, and the repetition is what is kept.

- The first campaign of `run_head_to_head.py` ran while audit agents and commands of mine were running: its single-threaded whole-genome runs took 43.6 to 61.6 seconds, against 43.3 to 55.1 in the one that is kept, made with nothing else running (the 55.1 is round 1: its CPU time, 44.6 seconds, is that of the other four, 43.6 to 44.3, so it spent about ten seconds waiting).
- The first run of the five organisms overlapped the sleeps of the machine: macOS idle-sleeps between commands and suspends a background job while it does (`pmset -g log` shows idle sleeps from 00:23 to 01:32 EST, with only brief wakes), so human took 18.7 seconds against the 12.4 kept (yeast 0.50 against 0.47, Drosophila 9.28 against 9.62, Arabidopsis 26.2 against 27.0, mouse 65.9 against 65.4).
  The kept run is under `caffeinate -i`.
- The first run of `run_sa_formats.py` shared the machine with a `docker desktop start` that kept polling for an engine that never came up (about half a core): its builds took the same time to within 0.2 seconds and its annotation 7.4 seconds against 7.1.
  A second run, before REVEL was added to the annotation and the two formats were made to alternate, took 6.9 and 6.8 seconds.
- The first timing of the ACMG annotation (`19_time_annotation.py`) overlapped a search of the file system of mine during its one-thread run: 84.9 seconds against the 74.4 kept.

## What was not measured

Ensembl VEP.
Docker Desktop could not be started on the day (its virtual machine was already dead when the work began, a `docker desktop restart` did not bring it back, and the engine did not come up when the application was started by hand; a staged update of Docker Desktop, build 244584, has been `in_progress` since 26 October).
So the Ensembl VEP times of the README and the manuscript are those of 29 July 2026 (`manuscript/data/vep_baseline_2026-07-29.csv`); its start-up and peak memory are not measured, and the accuracy counts are those of 27 October 2026, which `identity_vs_5048fc5/` shows this build to reproduce.
`benchmarks/run_vep_side.py` is the Ensembl VEP side of the speed tables, written for the day Docker is back; it has not been run.
To measure the accuracy again, run `validation/run_tab_world.sh` and the commands of Part 3 of `docs/VEP_DIVERGENCE.md`, and replace the counts of `manuscript/data/vep116_agreement.csv` and `vep116_hg002_columns.csv`.

## To repeat it

```bash
cargo build --release -p fastvep-cli                       # or: FASTVEP_BIN=... for a build made elsewhere
python3 benchmarks/run_head_to_head.py --bin target/release/fastvep --built-from "v0.5.0-rc.1 (c12d9f9)" --out benchmarks/results/v0.5.0-rc.1
RUNS=5 FASTVEP_BIN=target/release/fastvep caffeinate -i bash benchmarks/run_benchmark.sh
target/release/fastvep annotate -i test_data/organisms/human/human_giab_hg002_full.vcf --gff3 ... --fasta ... --hgvs --symbol --canonical -o - \
    | python3 benchmarks/consequence_distribution.py - --out-prefix benchmarks/results/v0.5.0-rc.1/consequence
python3 manuscript/code_metrics.py --out benchmarks/results/v0.5.0-rc.1
python3 benchmarks/run_sa_formats.py --bin target/release/fastvep --out benchmarks/results/v0.5.0-rc.1
python3 benchmarks/check_lean_inputs.py --bin target/release/fastvep --out benchmarks/results/v0.5.0-rc.1
python3 validation/compare_builds.py --base OLD --base-label 5048fc5 --new target/release/fastvep --new-label v0.5.0-rc.1 --out benchmarks/results/v0.5.0-rc.1/identity_vs_5048fc5
python3 manuscript/build_data.py --results benchmarks/results/v0.5.0-rc.1 && python3 manuscript/verify_numbers.py
```

Run nothing else while a timing script runs, and hold a laptop awake with `caffeinate -i`.
