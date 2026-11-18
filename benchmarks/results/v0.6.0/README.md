# Measurements of fastVEP 0.6.0

What depends on the build that the release is cut from is here, measured with the binary that `environment.tsv` describes.
Every speed, memory and agreement number that the README and `manuscript/manuscript.md` quote and that this directory does not hold is in [`../v0.5.0/`](../v0.5.0/), measured with the release candidate of 0.5.0.
Those numbers describe 0.6.0 as well, for two reasons that are measured:

- 0.6.0 writes the same bytes as master and as a rebuild of the release candidate on the inputs of the comparison with Ensembl VEP 116.3 (`identity_vs_master/`, `identity_vs_rc1/`) and on every input and with every flag of the benchmarks (25 cases, [`../identity_sweep/`](../identity_sweep/)), apart from what follows Ensembl VEP 116.3 and is listed there: the `PICK` cell of the `intergenic_variant` rows under `--flag-pick*`, and the `most_severe_consequence` of 3 JSON records of HG002.
- It takes the same time: the fastest runs of 16 workloads differ by -1.0% to +1.4% from master's, and the instructions of one thread by -0.39% to +0.23% in 6 of the 8 workloads ([`../small_variant_regression/`](../small_variant_regression/)).

The speed of structural variants, which 0.5.0 did not annotate as Ensembl VEP does, is in [`../sv_head_to_head/`](../sv_head_to_head/), and what Ensembl VEP 116.3 itself wrote for a generated world, for the three records above and for the flagging forms of `--pick` is in [`../vep_oracle/`](../vep_oracle/).

## What is here

| File | Made by | What it holds |
|------|---------|---------------|
| `environment.tsv` | `shasum`, `rustc --version`, `sw_vers` | the binary (`fastvep 0.6.0`, built from the release commit with `cargo build --release --locked -p fastvep-cli`), its SHA-256 and its size, which the manuscript quotes as the size of the executable |
| `code_metrics.csv`, `so_terms.csv`, `crate_dependencies.csv` | `manuscript/code_metrics.py --out benchmarks/results/v0.6.0` | the lines, tests, consequence terms (41, of which 33 are emitted) and crate dependencies that the manuscript counts |
| `tests_executed.csv` | `cargo test --workspace --no-fail-fast --locked` | what ran: 1,889 passed, 0 failed, 1 ignored, in 130 test targets; the suite was run three times on the release commit, once alone and twice while the machine did other work, with the same result |
| `identity_vs_master/`, `identity_vs_rc1/` | `validation/compare_builds.py` | whether this build writes the same bytes as master and as a rebuild of `v0.5.0-rc.1` on the 11 inputs of the comparison with Ensembl VEP 116.3 that can be run without it |

`manuscript/build_data.py --results benchmarks/results/v0.5.0 --build benchmarks/results/v0.6.0` turns these and the runs of `v0.5.0` into `manuscript/data/`, and `manuscript/verify_numbers.py` checks the text against them.

## To repeat it

```sh
cargo build --release --locked -p fastvep-cli
shasum -a 256 target/release/fastvep; ls -l target/release/fastvep
manuscript/code_metrics.py --out benchmarks/results/v0.6.0
cargo test --workspace --no-fail-fast --locked
validation/compare_builds.py ...      # see its header: master, the release candidate and this build
```
