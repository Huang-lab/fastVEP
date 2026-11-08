# The complete databases of the benchmark

The benchmark reads nine evidence sources, and three of them, gnomAD, SpliceAI and PhyloP, are built from the gnomAD v4.1 exomes sites file.
Since 8 November 2026 all three are built from the whole of that file, and the benchmark reads those databases.
This directory records the 72 files, how they were built, and how they were checked.

## What the files are

- `gnomad_chr*.osa2`, `spliceai_chr*.osa2` and `phylop_chr*.osa2` for chromosomes 1 to 22, X and Y: 72 files, 2,718,429,803 bytes.
- All are `.osa2` files of format version 2 for GRCh38, with 16-bit chunks (`chunk_bits` is 16 in each file's `fastsa/metadata.json`).
  The first databases of the benchmark have 20-bit chunks for gnomAD and are `.osa` files of the first format for SpliceAI and PhyloP.
  The reader of the tag `v0.5.0-rc.1`, and the readers of the commits 9d84be0 and 5d8004a, read `chunk_bits` from each file.
- `manifest.tsv` lists each file with its size, its SHA-256 and the fields of its metadata, and `SHA256SUMS` is the same checksums in the form `sha256sum -c` reads.

## How they were built

They were built on the host that serves fastvep.org, on 8 November 2026, from the whole of the gnomAD v4.1 exomes sites files and not from regions of them.
SpliceAI and PhyloP were distilled from the `spliceai_ds_max` and `phylop` fields of those files by `scripts/extract_gnomad_scores.py` without its `--regions` option.
All three were built with `fastvep sa-build --source gnomad`, `--source spliceai` and `--source phylop` with `--chunk-bits 16`, by the binary `0.5.0-rc.1+web.1ea790c`.
This is the record of whoever built them: the build was not repeated on the machine of the benchmark, and the checksums are what ties the copies here to the files of the host.
The copies were compared with the host's files by their SHA-256 on 8 November 2026, and every one of the 72 was identical.

## How the benchmark reads them

`scripts/sa_sources/make_complete_stack.sh` makes `data/benchmark/sa_db_complete/`, a directory of links to every file of the first stack `data/benchmark/sa_db/` except its SpliceAI, PhyloP and gnomAD files, and to these 72 files in their place.
It stops if any of the 72 files lacks the checksum of `SHA256SUMS`.
REVEL, ClinVar, the ClinVar protein index, ClinGen, gnomAD's gene constraints and RepeatMasker are the same files as in the first stack, because their dates are part of the design of the benchmark.
`run_benchmark.sh` and `11_manuscript_runs.py` read `sa_db_complete` by default, and `SA_DIR` and `BENCH_SA_DB` name another stack.

## Why the first databases were replaced

The first SpliceAI and PhyloP databases were distilled in April 2026 from regions of the same gnomAD file, and they lack records that the file carries.
On the 673,660 variants of the benchmark, 486,580 have a SpliceAI value under the first databases and 529,992 under the complete ones, and 532,531 have a PhyloP value under the first and 580,460 under the complete.
No variant has a value only under the first databases, and where both have a value the values are the same.
gnomAD and REVEL have the same variants with a value under both stacks, and the same values.
`complete_vs_april_databases.md` in the directory above sets what the runs say under each stack side by side.

## The check of gnomAD

The complete gnomAD was adopted for every run only after the April set had been run with and without it.
With the complete SpliceAI and PhyloP databases and the tag's binary, the benchmark's own regional gnomAD databases and the complete gnomAD gave the same ACMG call, the same criteria and the same SpliceAI, PhyloP, gnomAD, REVEL and ClinVar fields for each of the 673,660 variants (8 November 2026, `fastvep 0.5.0-rc.1` of `data/benchmark/manuscript_bins_v38/fastvep_v38`).
