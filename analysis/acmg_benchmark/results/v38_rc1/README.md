# The runs of v38 made again by the release candidate

These are the 38 runs of Section 3.8 of the ACMG manuscript on the April 2026 release of ClinVar and the 4 on the October release, made again on 7 November 2026 by binaries built from the tag `v0.5.0-rc.1` (`c12d9f9`), to check that the code that is merged gives the calls the manuscript reports.
`analysis/acmg_benchmark/RUN_VERSIONS.md`, "The release candidate reproduces v38", is the account of it.

| File | What it holds |
|------|---------------|
| `runs_v38.tsv`, `holdout_runs_v38.tsv`, `run_matrix_v38.tsv` | the concordance matrix, the firing counts and the description of every run; they are the files of `../v37_manuscript/` line for line |
| `provenance/` | what made each run: the flags, the hash of the binary and of every configuration, the inputs read, the signature of the database directory, the commit of the checkout and the time |
| `configs/` | the `--acmg-config` of each run that needs one |
| `rerun_comparison.csv` | `scripts/16_compare_runs.py` against `data/benchmark/manuscript_runs/`: 42 of 42 runs identical in every scored file, and the 15 annotated VCFs that both keep identical record by record |
| `annotation_time.tsv` | `scripts/19_time_annotation.py`: the annotation of the benchmark on every core and on one thread, and its scoring |

`runs.tsv` and `holdout_runs.tsv` are here as the harness writes them and hold only their header, because the runs of the v37 family were not repeated; `run_matrix.tsv` describes every run of the harness and is the file of `../v37_manuscript/` line for line.
The annotated VCFs (about 250 MB each; the per-variant tables were not made for these runs) are not in the repository; they are in `data/benchmark/manuscript_runs_rc1/` on the machine that made them.
