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
| [`METHODS.md`](METHODS.md) | What the benchmark measures, how the truth set is built, what the numbers mean |
| [`RUN_VERSIONS.md`](RUN_VERSIONS.md) | Per-run history: what changed between v1 and v17, with measured effect |
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

# One-time: build the truth set
python3 analysis/acmg_benchmark/scripts/01_extract_clinvar_2star.py

# Each run - the version string names the output directory
bash analysis/acmg_benchmark/scripts/run_benchmark.sh v15

# Or with a config, to measure one change against the same truth set
ACMG_CONFIG=analysis/acmg_benchmark/data/vcep_thresholds.toml \
  bash analysis/acmg_benchmark/scripts/run_benchmark.sh v17
```

The run script takes about 90 seconds to annotate 673,660 variants plus a couple
of minutes to score, reuses an existing annotated VCF if one is present, and
writes to `data/benchmark/output_v15/`.
Record what changed, and its measured effect, in [`RUN_VERSIONS.md`](RUN_VERSIONS.md).

To vary one threshold instead of one version, use the sweep harness rather than a
hand-edited config:

```bash
python3 analysis/acmg_benchmark/scripts/sweep_acmg_thresholds.py \
  --sa-dir data/benchmark/sa_db --out-dir data/benchmark/sweeps \
  --sweep "bp7_max_intron_offset=50,100,300,1000"
```

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
