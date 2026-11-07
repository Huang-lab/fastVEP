# Measurements of fastVEP 0.5.0

Everything in this directory was measured on 7 November 2026 (UTC) with the binary built from the tag `v0.5.0-rc.1` (commit `c12d9f9`, whose crates are those of `master` at `aa5e5ee`; SHA-256 `0063e01bd143...`, 5,834,640 bytes), on an Apple M4 (10 cores: 4 performance and 6 efficiency; 24 GiB; macOS 26.5.1; rustc 1.98.0).
What the README and `manuscript/manuscript.md` quote about speed, memory, the five organisms, the consequence distribution, the fastSA formats and the tests is in these files, and `manuscript/build_data.py` turns them into `manuscript/data/`.
`environment.tsv` holds the machine, the load, the hash of the binary and of every input.
The binary prints `fastvep 0.5.0-rc.1`, because it was built from the tag of the release candidate.
The numbers are those of the release if it is cut from the same crates, which `git diff --stat v0.5.0-rc.1 v0.5.0 -- crates Cargo.lock` shows once the tag `v0.5.0` exists (it is empty if they are); if it is not empty, the commands under "To repeat it" are to be run again.
`web/index.html` is embedded in the binary, and the description in its head, which quotes the speed-up, is the one difference to expect there.

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
| `with_vep/raw.tsv`, `summary.tsv`, `environment.tsv`, `identity.tsv` | `benchmarks/run_head_to_head.py --only startup,scaling,memory,formats,wgs,wgs_joint` | the parts that are compared with Ensembl VEP, run again in the session of the VEP runs: Tables 2 to 5, and the whole genome with ClinVar and the joint release of gnomAD (Table 3's third row) |
| `with_vep/vep_raw.tsv`, `vep_summary.tsv`, `vep_genome.tsv`, `vep_genome_sums.csv`, `vep_baseline.csv`, `vep_resources.csv`, `vep_environment.tsv` | `benchmarks/run_vep_side.py` | Ensembl VEP 116.3 on the same machine: every run with the CPU time of VEP beside the wall time, the medians, one run per chromosome and condition and their sums, what the tables quote, and the image, its size and the `--custom` arguments |
| `with_vep/custom_matches_chr22.tsv` | `benchmarks/check_custom_matches.py` | whether the two tools find the same records in the same ClinVar and gnomAD files, on chromosome 22 |
| `vep116_agreement.csv`, `vep116_hg002_columns.csv`, `vep116_refseq_columns.csv`, `vep116_run.tsv`, `vep116_diff/` | `validation/measure_vep116.py` | the counts of the comparison with Ensembl VEP 116.3 on every input, the HG002 sample by column and by cause, what ran, and every row that differs |
| `identity_vs_5048fc5/`, `identity_vs_ac2e2b6/` | `validation/compare_builds.py` | whether this build writes the same bytes as the build of 27 October 2026 (`5048fc5`) and of 15 September 2026 (`ac2e2b6`) on the inputs of the comparisons with Ensembl VEP |

The ACMG runs of 0.5.0 are in `analysis/acmg_benchmark/results/v38_v0.5.0/` (`RUN_VERSIONS.md` there, "Release 0.5.0 reproduces v38").

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

## Ensembl VEP, and the second session

The first session of these measurements (the night of 6 to 7 November) ran without Docker: its virtual machine was already dead when the work began, a `docker desktop restart` from the shell did not bring it back, and the engine did not come up when the application was started by hand.
After Docker Desktop was opened from its window on 7 November, Ensembl VEP 116.3 (the image `ensemblorg/ensembl-vep@sha256:2a46b0ba4ac5...`) was measured on the same machine, and the fastVEP measurements that are compared with it were made again beside it (times are UTC).

- **Accuracy** (`validation/measure_vep116.py`, 12:59 to 13:59, with the binary of the first session): every count of the comparison with 116.3 is in `vep116_agreement.csv` and `vep116_hg002_columns.csv`, the rows that differ are in `vep116_diff/`, and `vep116_run.tsv` says what ran.
  The counts that `docs/VEP_DIVERGENCE.md` already had are the ones that commit `5048fc5` gave on 27 October 2026, to the cell.
  The comparison of RefSeq's GFF3 is the exception: it is new here, and its count differs from the 135 of 81,510 cells that the document had, which was made before fastVEP read a transcript as VEP's cache holds it (`vep116_refseq_columns.csv` has the 818 of 73,359 cells by column).
- **Speed, second session** (`with_vep/`): `run_head_to_head.py` ran start-up, scaling, memory, formats, `wgs` and `wgs_joint` again (14:10 to 14:39), and `run_vep_side.py` ran Ensembl VEP: the chromosome 22 series, start-up, peak memory and the size of the image, then the whole genome, one run at a time, from 15:06 to 20:48.
  The second session was needed because the first ran slower than the second: the medians of the consequence-only genome on one thread were 44.7 s and 42.2 s, and those of the condition with gnomAD's exomes 61.6 s and 51.1 s.
  A ratio of a time of the one to a time of the other would carry that difference, so every ratio of the tables is of one session.
  A first run of the joint condition, started at 14:00 and stopped after four of its five rounds (`with_vep/stopped_run_joint.log`), had shown the difference: its consequence-only control took 40.9 to 44.6 s.
- **Ensembl VEP 116.3 against 115.1**: release 116.3 takes 1.3 to 1.9 times as long as release 115.1 on the four inputs of the chromosome 22 series, measured on the same day with the same flags (`with_vep/vep_115_1/`); the 13.41 s that 115.1 took for the 10,000 variants in July 2026 is reproduced (12.9 s on 7 November).
- **One run at a time**: two VEP runs at a time cost 10% each and three at a time 16% (chromosome 22, 50,284 variants: 112.5 s alone, 124 s two at a time, 131 s three at a time, CPU time inside the container up by the same share), so the whole-genome runs were made one after the other, once for each chromosome and condition.
  Their wall times are `vep_genome.tsv`, their sums `vep_genome_sums.csv` (6,516 s, 6,662 s and 7,334 s for the three conditions), and the CPU time of VEP inside the container is 99.6% of the wall time.
  `custom_matches_chr22.tsv` says that both tools found the same 1,057 ClinVar records and 2,623 gnomAD records in the same files on chromosome 22.

What was not measured again: the comparisons with VEP's own Ensembl and RefSeq caches (chromosomes 21 and 22), because the extracted copies of the caches, taken from the 26 to 28 GB archives of Ensembl's FTP, are no longer on disk.
Their numbers in `docs/VEP_DIVERGENCE.md` are those of 27 October 2026 with commit `5048fc5`.

## To repeat it

```bash
cargo build --release -p fastvep-cli                       # or: FASTVEP_BIN=... for a build made elsewhere
python3 benchmarks/run_head_to_head.py --bin target/release/fastvep --built-from "v0.5.0-rc.1 (c12d9f9)" --out benchmarks/results/v0.5.0
python3 benchmarks/run_head_to_head.py --bin target/release/fastvep --built-from "v0.5.0-rc.1 (c12d9f9)" --out benchmarks/results/v0.5.0/with_vep \
    --only startup,scaling,memory,formats,wgs,wgs_joint                                      # the second session, beside Ensembl VEP
python3 benchmarks/run_vep_side.py --out benchmarks/results/v0.5.0/with_vep                    # needs Docker; about seven hours, one run at a time
python3 benchmarks/run_vep_side.py --out benchmarks/results/v0.5.0/with_vep/vep_115_1 --image ensemblorg/ensembl-vep:release_115.1 --only scaling,startup
python3 benchmarks/check_custom_matches.py --bin target/release/fastvep --chrom 22 --out benchmarks/results/v0.5.0/with_vep
python3 validation/measure_vep116.py --bin target/release/fastvep --out benchmarks/results/v0.5.0 --refseq-gff3 GCF_000001405.40_GRCh38.p14_genomic.gff.gz
RUNS=5 FASTVEP_BIN=target/release/fastvep caffeinate -i bash benchmarks/run_benchmark.sh
target/release/fastvep annotate -i test_data/organisms/human/human_giab_hg002_full.vcf --gff3 ... --fasta ... --hgvs --symbol --canonical -o - \
    | python3 benchmarks/consequence_distribution.py - --out-prefix benchmarks/results/v0.5.0/consequence
python3 manuscript/code_metrics.py --out benchmarks/results/v0.5.0
python3 benchmarks/run_sa_formats.py --bin target/release/fastvep --out benchmarks/results/v0.5.0
python3 benchmarks/check_lean_inputs.py --bin target/release/fastvep --out benchmarks/results/v0.5.0
python3 validation/compare_builds.py --base OLD --base-label 5048fc5 --new target/release/fastvep --new-label v0.5.0-rc.1 --out benchmarks/results/v0.5.0/identity_vs_5048fc5
python3 manuscript/build_data.py --results benchmarks/results/v0.5.0 && python3 manuscript/verify_numbers.py
```

Run nothing else while a timing script runs, and hold a laptop awake with `caffeinate -i`.
