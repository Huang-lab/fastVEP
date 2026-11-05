# ACMG concordance benchmark

How well does `fastvep annotate --acmg` agree with expert-curated ClinVar calls?

This directory is the single home for that benchmark's **code, documentation and
results**.
Everything it reads and writes at scale - input VCFs, the ~30 GB of supplementary
annotation databases, and the full per-run output directories - lives under
`data/benchmark/`, which is entirely gitignored.

That split is the point of the layout, and it is worth stating why, because the
repository spent several months with the opposite arrangement.
Git cannot re-include a file whose parent directory is excluded.
While `.gitignore` excluded `/data/` and tracked files lived beneath it, every
one of them needed `git add -f`, and three were silently missed - a benchmark
runner and two analysis write-ups that simply never entered the repository.
Nothing under `data/` is tracked now, and nothing needs to be.

**If it is code, documentation or a small result, it belongs here. If it is an
input, a database or a bulk output, it belongs in `data/benchmark/`.**

## Layout

| Path | What |
|---|---|
| [`METHODS.md`](METHODS.md) | What the benchmark measures, how the ClinVar 2-star+ set is built, what the numbers mean (ClinVar's classifications are a reference, not ground truth) |
| [`RUN_VERSIONS.md`](RUN_VERSIONS.md) | Per-run history: what changed between v1 and v38, with measured effect, and the record of the full rerun of 4 November 2026 |
| [`STATUS.md`](STATUS.md) | On-disk inventory - what is downloaded and built locally |
| [`scripts/`](scripts/) | Everything executable, plus [its own README](scripts/README.md) |
| [`scripts/sa_sources/`](scripts/sa_sources/) | Building the supplementary annotation databases from public sources |
| [`results/`](results/) | Small tracked artifacts per run: concordance matrices, summaries, figures |
| [`data/`](data/) | Small tracked *derived* data: curated tables generated from public sources, not bulk inputs |
| [`notes/`](notes/) | Analysis write-ups that are not tied to a single run |

For the *speed* and VEP-equivalence benchmarks, see [`benchmarks/`](../../benchmarks/)
instead. This directory is only about clinical concordance.

## Running it

```bash
# One-time: fetch sources (~30 GB, into data/benchmark/sa_sources/) and build databases
bash analysis/acmg_benchmark/scripts/sa_sources/download_sa_sources.sh

# One-time: extract ClinVar's 2-star+ classifications
python3 analysis/acmg_benchmark/scripts/01_extract_clinvar_2star.py

# Each run - the version string names the output directory
bash analysis/acmg_benchmark/scripts/run_benchmark.sh v15

# Or with a config, to measure one change against the same ClinVar 2-star+ set.
# A variant of a run is a run of its own, so it gets its own version name.
ACMG_CONFIG=analysis/acmg_benchmark/data/vcep_thresholds.toml \
  bash analysis/acmg_benchmark/scripts/run_benchmark.sh v15_vcep

# Any other `fastvep annotate` flag, the same way
EXTRA_ARGS="--pick-order mane_plus_clinical,mane_select,canonical,appris,tsl,biotype,ccds,rank" \
  bash analysis/acmg_benchmark/scripts/run_benchmark.sh v15_mpc

# Any number of runs against a baseline, call by call
python3 analysis/acmg_benchmark/scripts/09_compare_runs.py \
  --baseline v15 --runs v15_vcep v15_mpc --out /tmp/v15_comparison
```

The run script takes about 2 minutes to annotate 673,660 variants plus a couple
of minutes to score, and writes to `data/benchmark/output_v15/`.
It records what produced the annotation in `run_provenance.txt` beside it: the
flags, the hash of the binary and of every `--acmg-config` and
`--functional-evidence` file, which inputs were read (by name and size, not
content) and the commit.
It reuses an existing annotated VCF only when the same binary, flags,
configuration and inputs made it, and stops with exit code 3 otherwise, so a
second run under a name already used cannot be scored against the first run's
annotation.
Record what changed, and its measured effect, in [`RUN_VERSIONS.md`](RUN_VERSIONS.md).

To vary one threshold instead of one version, use the sweep harness rather than a
hand-edited config:

```bash
python3 analysis/acmg_benchmark/scripts/sweep_acmg_thresholds.py \
  --sa-dir data/benchmark/sa_db --out-dir data/benchmark/sweeps \
  --sweep "bp7_max_intron_offset=50,100,300,1000"
```

## The ACMG manuscript

[`manuscript/acmg/`](../../manuscript/acmg/) holds a manuscript on how the classifier was built, reviewed and measured, and on its discordance with ClinVar 2-star+ classifications in its final version.
Its tables are produced from the scored runs by the scripts below, in this order, and its text is checked against them.
Every run is made by a binary that an environment variable names, and not by the checkout's own `target/release/fastvep`.
The checkout's classifier is v38, and a run called `v37_...` made by it would be a run of v38 under v37's name, which `11_manuscript_runs.py` and `14_temporal_holdout.py` refuse to do.
The provenance file of a run records the commit and the number of uncommitted files of the checkout that ran the harness, which is not the commit of the binary; `results/v37_manuscript/binaries.tsv` says what each binary was built from.

```bash
# 0. The six binaries: v37; the two measurement builds of v37 (replay_validation/combiner_cutoff_env.patch and
#    bp_hooks_env.patch); the commit before the change of Section 3.8; the change; and the change with the alternatives
#    (replay_validation/inheritance_alternatives_env.patch). Each is built from a clean export of its commit in a
#    target directory of its own, and the script stops if a patched build does not contain its patch.
#    V38_TREE=. builds the checkout as it is, uncommitted files included; V38_REF=<commit> builds a commit.
#    It prints the six `export` lines that the steps below read their binaries from.
V38_TREE=. analysis/acmg_benchmark/scripts/build_manuscript_binaries.sh /path/to/binaries
export FASTVEP_V37_BIN=/path/to/binaries/fastvep_v37 FASTVEP_CUTOFF_BIN=/path/to/binaries/fastvep_cut \
       FASTVEP_HOOKS_BIN=/path/to/binaries/fastvep_hook FASTVEP_V38_BASE_BIN=/path/to/binaries/fastvep_base \
       FASTVEP_V38_BIN=/path/to/binaries/fastvep_v38 FASTVEP_V38_ALT_BIN=/path/to/binaries/fastvep_alt

# 0b. v37 itself, which every run of the April set differs from by one setting, and the ten runs of round 6 that are
#     v37 with a configuration or a flag changed (the configurations are in results/v37/configs)
FASTVEP="$FASTVEP_V37_BIN" bash analysis/acmg_benchmark/scripts/run_benchmark.sh v37
for n in bs2off bs2prev_3e-5 bs2prev_1e-4 bs2prev_1e-3 bs2prev_3e-3 bs2ac_20 bs2ac_100 bs2ac_1000; do
  ACMG_CONFIG=analysis/acmg_benchmark/results/v37/configs/$n.toml FASTVEP="$FASTVEP_V37_BIN" \
    bash analysis/acmg_benchmark/scripts/run_benchmark.sh v37_$n
done
FASTVEP="$FASTVEP_V37_BIN" EXTRA_ARGS="--pick-order mane_plus_clinical,mane_select,canonical,appris,tsl,biotype,ccds,rank" \
  bash analysis/acmg_benchmark/scripts/run_benchmark.sh v37_mpc
FASTVEP="$FASTVEP_V37_BIN" EXTRA_ARGS="--pick-order rank,mane_select,mane_plus_clinical,canonical,appris,tsl,biotype,ccds" \
  bash analysis/acmg_benchmark/scripts/run_benchmark.sh v37_rank

# 1. The three runs that change an input need the ClinGen file completed first
#    (the commands are in the docstring of 15_clingen_complete_export.py)

# 2. The sweeps, ablations, interaction and gap runs, and the runs of the two measurement builds, which are skipped
#    where their variable is not set (about 3 minutes each; --workers 3 is comfortable)
python3 analysis/acmg_benchmark/scripts/11_manuscript_runs.py --workers 3

# 3. The later ClinVar release (this also writes its per-variant table), and the seven runs of the changes that
#    Section 4.5 of the manuscript proposes, made on the variants of that release
python3 analysis/acmg_benchmark/scripts/14_temporal_holdout.py all
python3 analysis/acmg_benchmark/scripts/11_manuscript_runs.py --holdout --workers 3

# 3b. The per-variant table of v37 and of each run whose variants the statistics read one by one
python3 analysis/acmg_benchmark/scripts/12_per_variant_table.py \
  --run data/benchmark/output_v37 --out data/benchmark/manuscript_runs/v37_per_variant.tsv.gz
for run in v37_gap_clingen_complete v37_gap_orphanet_lof v37_gap_both v37_hk_bp1_loo v37_sw_phylop_1000 \
           v37_hk_bp7_walker v37_sw_bp1veto_1 v37_sw_bp1veto_2 v37_hk_bp1_loo_veto1 v37_hk_bp1_loo_veto2; do
  python3 analysis/acmg_benchmark/scripts/12_per_variant_table.py \
    --run data/benchmark/manuscript_runs/$run --out data/benchmark/manuscript_runs/${run}_per_variant.tsv.gz
done
for run in v37_oct_clingen_complete v37_oct_orphanet_lof v37_oct_both v37_oct_bp7_walker v37_oct_bp1_loo v37_oct_bp1_veto1 v37_oct_no_an_floor; do
  python3 analysis/acmg_benchmark/scripts/12_per_variant_table.py \
    --run data/benchmark/manuscript_runs/$run --labels data/benchmark/clinvar_20261004_covered/clinvar_2star_labels.tsv \
    --out data/benchmark/manuscript_runs/${run}_per_variant.tsv.gz
done

# 3c. Why BA1 and BS1 did not fire on the Benign variants of 2% or more that are not called benign-side (needs the tables above)
python3 analysis/acmg_benchmark/scripts/14_temporal_holdout.py blockers

# 4. The runs of Section 3.8 (what changed after v37: ClinGen's curated mode of inheritance is read), which need the
#    three binaries of step 0 that are not v37's. The databases of gnomAD's joint exomes and genomes release come first,
#    for all 24 chromosomes: about 6 hours from gnomAD's public files, limited by the network and not by the script
#    (JOBS=8 runs eight chromosomes at a time, the default is four), and resumable, so that a run that was stopped
#    continues after its last finished batch. 11_manuscript_runs.py makes `v38_joint` once all 24 exist; the pilot on
#    chr21 and chr22 that was made first reads the same two databases, and its input, the benchmark's variants on those
#    two chromosomes, is written by 11_manuscript_runs.py from the benchmark's own files
FASTVEP="$FASTVEP_V38_BIN" JOBS=8 analysis/acmg_benchmark/scripts/sa_sources/build_gnomad_release_per_chrom.sh joint
python3 analysis/acmg_benchmark/scripts/11_manuscript_runs.py --v38 --workers 3
python3 analysis/acmg_benchmark/scripts/11_manuscript_runs.py --v38 --holdout --workers 3

# 4c. The per-variant table of each run whose variants the statistics read one by one
for run in v38 v38_base v38_alt_unguarded v38_alt_label_guard v38_alt_bs2_as_v37 v38_alt_triple_cap \
           v38_noflags v38_joint v38_chr2122_exomes v38_chr2122_exomes_noflags v38_chr2122_joint; do
  python3 analysis/acmg_benchmark/scripts/12_per_variant_table.py \
    --run data/benchmark/manuscript_runs/$run --out data/benchmark/manuscript_runs/${run}_per_variant.tsv.gz
done
for run in v38_oct_base v38_oct_shipped v38_oct_alt_unguarded v38_oct_alt_triple_cap; do
  python3 analysis/acmg_benchmark/scripts/12_per_variant_table.py \
    --run data/benchmark/manuscript_runs/$run --labels data/benchmark/clinvar_20261004_covered/clinvar_2star_labels.tsv \
    --out data/benchmark/manuscript_runs/${run}_per_variant.tsv.gz
done

# 4d. The variants that the joint release turns into errors, or out of a pathogenic-side call, annotated one at a time
#     through both databases (the filtering allele frequency and the flags are in the json and not in the VCF)
FASTVEP_V38_BIN="$FASTVEP_V38_BIN" python3 analysis/acmg_benchmark/scripts/17_joint_boundary_cases.py

# 4e. The status of the 28 criteria against the engine's json output, on every 40th record of the benchmark (the VCF lists
#     only the criteria that were met; the json says whether each could be evaluated). It stops, writing nothing, if the engine
#     disagrees with the groups of criteria that the script states
FASTVEP_V38_BIN="$FASTVEP_V38_BIN" python3 analysis/acmg_benchmark/scripts/18_criteria_status.py

# 5. Every table, then the figures, then the checks
python3 analysis/acmg_benchmark/scripts/13_manuscript_stats.py
python3 manuscript/acmg/generate_figures.py
python3 manuscript/acmg/build_tables.py --write
python3 manuscript/acmg/build_tables.py --check
python3 manuscript/acmg/build_supplement.py --write
python3 manuscript/acmg/verify_numbers.py
```

`run_benchmark.sh` and the scripts that start it read the environment, so a run that is made by hand needs `FASTVEP` set the same way.
The tables and the checks are deterministic: after the statistics script had run under each of two hash seeds (`PYTHONHASHSEED=0` and `1`) the 89 files of `manuscript/acmg/data/` were the same, and the figures are byte for byte the same on a rerun.
Some files are not written by the statistics script: `change_ledger.csv` and `round6_rulings.tsv` are typed from the project's record, `clingen_missing_curations.csv` comes from `15_clingen_complete_export.py`, `temporal_high_af_blockers_v37.csv` from `14_temporal_holdout.py blockers`, `v38_joint_errors.csv` and `v38_joint_pm2_lost.csv` from `17_joint_boundary_cases.py`, and `criteria_status.csv` from `18_criteria_status.py`.

The inputs are dated and none of them is in the repository, and the download scripts fetch the current files and not the dated ones: the ClinVar weekly releases of 4 April 2026 (the benchmark and its labels) and of 4 October 2026, the ClinGen Gene-Disease Validity exports of 30 April 2026 (the file the benchmark reads) and of 11 October 2026, and the ClinVar `variant_summary` and gnomAD files of Section 2.2 of the manuscript.
A rerun on other releases of them is a rerun of the method and not of the paper's numbers.
The run records are in [`results/v37_manuscript/`](results/v37_manuscript/).
`scripts/16_compare_runs.py` compares a rerun with the runs it repeats, file by file and record by record, and writes `results/v37_manuscript/rerun_comparison.csv`; use it when any step above is made again.

## The classifier itself

Thresholds, criteria and the reasoning behind every value set by measurement are
documented in [`docs/ACMG.md`](../../docs/ACMG.md), not here.
This directory measures the classifier; it does not define it.

The medical-genetics review that drove most of the criteria changes is in
[`docs/ACMG_EXPERT_REVIEW_ROUND2.md`](../../docs/ACMG_EXPERT_REVIEW_ROUND2.md).
What each subsequent run did about it is in
[`RUN_VERSIONS.md`](RUN_VERSIONS.md); the runs between reviews are ours, not hers.

## Curated tables in `data/`

Two tables are generated from public sources and tracked, because they are small,
because they are inputs to classification rather than outputs of a run, and
because their provenance needs to be reviewable line by line.

| File | What | Regenerate with |
|---|---|---|
| [`data/vcep_thresholds.toml`](data/vcep_thresholds.toml) | Published ClinGen VCEP BA1/BS1/PM2 bars for 117 genes, as `--acmg-config` input | [`scripts/build_vcep_thresholds.py`](scripts/build_vcep_thresholds.py) |
| [`data/vcep_thresholds_audit.tsv`](data/vcep_thresholds_audit.tsv) | Every bar and every rejection, with the verbatim sentence it came from | same script |
| [`data/gene_disease_attributes.tsv`](data/gene_disease_attributes.tsv) | Prevalence, onset and inheritance per gene-disorder pair, joined to the VCEP bars | [`scripts/build_gene_disease_attributes.py`](scripts/build_gene_disease_attributes.py) |

The audit table is the point. A threshold read out of prose by a regex is a
claim about a clinical guideline, and the only way to make that claim checkable
is to keep the sentence next to the number. Every rejection is in there too,
with its reason, so "this gene has no bar" and "this generator would not guess"
stay distinguishable.

`vcep_thresholds.toml` is **not** on by default. It is a measured improvement
that also introduces four new false-benign calls, two of which are founder
variants the panels excluded in prose the generator cannot read. See the v17
section of [`RUN_VERSIONS.md`](RUN_VERSIONS.md) for the full trade.
