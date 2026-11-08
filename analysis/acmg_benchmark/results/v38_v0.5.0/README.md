# The 42 runs of Section 3.8 made a second time

These are the 38 runs of Section 3.8 of the ACMG manuscript on the April 2026 release of ClinVar and the 4 on the October release, made a second time on 8 November 2026 by the same binaries as the first time (the tag `v0.5.0-rc.1`, `c12d9f9`, for the adopted change and its alternatives, and the build of the commit before the change for `v38_base` and `v38_oct_base`), with the complete databases, to check that the runs reproduce.
`analysis/acmg_benchmark/RUN_VERSIONS.md`, "The benchmark is re-based on the complete databases", is the account of it.
The repetition of 7 November 2026 by release 0.5.0 of the first drafts, which was made with the first databases, is "Release 0.5.0 reproduces v38" in the same file, and its outputs are in `data/benchmark/manuscript_runs_v0.5.0_april_databases/`.

| File | What it holds |
|------|---------------|
| `runs_v38.tsv`, `holdout_runs_v38.tsv`, `run_matrix_v38.tsv` | the concordance matrix, the firing counts and the description of every run; they are the files of `../v37_manuscript/` line for line |
| `provenance/` | what made each run: the flags, the hash of the binary and of every configuration, the inputs read, the signature of the database directory, the commit of the checkout and the time |
| `configs/` | the `--acmg-config` of each run that needs one |
| `rerun_comparison.csv` | `scripts/16_compare_runs.py` against `data/benchmark/manuscript_runs/`: 42 of 42 runs identical in every scored file, and the 15 annotated VCFs that both keep identical record by record |
| `annotation_time.tsv` | `scripts/19_time_annotation.py`: the annotation of the benchmark on every core and on one thread, and its scoring, measured alone on 7 November 2026 with the first databases and not repeated |

`runs.tsv` and `holdout_runs.tsv` are here as the harness writes them and hold only their header, because the runs of the v37 family were not repeated; `run_matrix.tsv` describes every run of the harness and is the file of `../v37_manuscript/` line for line.
The annotated VCFs of the repetition (about 250 MB each; the per-variant tables were not made for these runs) were deleted after `rerun_comparison.csv` had shown them identical to the first runs; the scored files and the provenance of every run are in `data/benchmark/manuscript_runs_v0.5.0/` on the machine that made them and in this directory.
