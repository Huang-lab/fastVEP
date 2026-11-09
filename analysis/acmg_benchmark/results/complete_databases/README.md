# The complete databases of the benchmark

The benchmark reads nine evidence sources, and three of them, gnomAD, SpliceAI and PhyloP, are built from the gnomAD v4.1 exomes sites file.
Since 8 November 2026 all three are built from the whole of that file, and the benchmark reads those databases.
This directory records the 72 files, how they were built, and how they were checked.

## What the files are

- `gnomad_chr*.osa2`, `spliceai_chr*.osa2` and `phylop_chr*.osa2` for chromosomes 1 to 22, X and Y: 72 files, 2,718,429,803 bytes.
- All are `.osa2` files of format version 2 for GRCh38, with 16-bit chunks (`chunk_bits` is 16 in each file's `fastsa/metadata.json`).
  The databases of April 2026 have 20-bit chunks for gnomAD and are `.osa` files of format version 1 for SpliceAI and PhyloP.
  The reader of the tag `v0.5.0-rc.1`, and the readers of the commits 9d84be0 and 5d8004a, read `chunk_bits` from each file.
- `manifest.tsv` lists each file with its size, its SHA-256 and the fields of its metadata, and `SHA256SUMS` is the same checksums in the form `sha256sum -c` reads.

## How they were built

They were built on the host that serves fastvep.org, on 8 November 2026, from the whole of the gnomAD v4.1 exomes sites files and not from regions of them.
SpliceAI and PhyloP were distilled from the `spliceai_ds_max` and `phylop` fields of those files by `scripts/extract_gnomad_scores.py` without its `--regions` option.
All three were built with `fastvep sa-build --source gnomad`, `--source spliceai` and `--source phylop` with `--chunk-bits 16`, by the binary `0.5.0-rc.1+web.1ea790c`.
This is the record of whoever built them: the build was not repeated on the machine of the benchmark, and the checksums are what ties the copies here to the files of the host.
The copies were compared with the host's files by their SHA-256 on 8 November 2026, and every one of the 72 was identical.

## How the benchmark reads them

`scripts/sa_sources/make_complete_stack.sh` makes `data/benchmark/sa_db_complete/`, a directory of links to every file of the stack of April 2026, `data/benchmark/sa_db/`, except its SpliceAI, PhyloP and gnomAD files, and to these 72 files in their place.
It stops if any of the 72 files lacks the checksum of `SHA256SUMS`.
REVEL, ClinVar, the ClinVar protein index, ClinGen, gnomAD's gene constraints and RepeatMasker are the same files as in the stack of April 2026, because their dates are part of the design of the benchmark.
`run_benchmark.sh` and `11_manuscript_runs.py` read `sa_db_complete` by default, and `SA_DIR` and `BENCH_SA_DB` name another stack.

## What the databases of April 2026 lack

The SpliceAI and PhyloP databases of April 2026 were distilled from regions of the same gnomAD file, and they lack records that the file carries.
Of the 673,278 variants of the benchmark that have an alternative allele (the other 382 are reference-agreement records, which get no call under either stack), 486,580 have a SpliceAI value under the databases of April 2026 and 529,992 under the complete ones, and 532,531 have a PhyloP value under the databases of April 2026 and 580,460 under the complete ones.
That is 8.2% of the SpliceAI values and 8.3% of the PhyloP values that exist.
No variant has a value only under the databases of April 2026, and where both have a value the values are the same.
gnomAD and REVEL have the same variants with a value under both stacks, and the same values.
`../complete_vs_april_databases.md` sets what the runs say under each stack side by side, and `../complete_vs_april_databases/databases_*.csv` hold the numbers.

## Where the databases of April 2026 lack their records

`scripts/sa_sources/extract_spliceai_phylop_from_gnomad.py` read the gnomAD file in batches of 200 consecutive regions of `data/benchmark/sa_sources/clinvar_2star_regions.bed` (the regions of one chromosome at a time), wrote SpliceAI and PhyloP from the same `tabix` stream, ignored the exit status of `tabix` and sent its messages to /dev/null.
The variants that only the complete databases have lie in 14 batches, the same 14 for SpliceAI and for PhyloP.
In each of them the databases of April 2026 have values up to a point and none after it, while the complete databases have values after it, and the first lost position is the same in both sources in all 14.
That is 43,412 SpliceAI records and 47,472 PhyloP records.
The other 457 PhyloP records are those of indels whose second base (position + 1) lies outside every region of the extraction.
No other record is missing from either source.
`../complete_vs_april_databases/databases_lost_blocks.csv` lists the 14 batches, and `databases_causes.csv` the counts.
The pattern is what a stream that ends early gives, in a script that does not look at the exit status and does not keep the messages.
It is not a record of such a stream: the extraction was not repeated, and the messages of April were not kept.
The regional gnomAD databases of the benchmark were built by `build_gnomad_per_chrom.sh`, which also ignores the exit status of `tabix` and also reads the regions 200 at a time, and they hold every record that the complete ones do for the benchmark's variants.

## The check of gnomAD

The complete gnomAD was adopted for every run only after the April set had been run with and without it.
With the complete SpliceAI and PhyloP databases and the tag's binary, the benchmark's own regional gnomAD databases and the complete gnomAD gave the same ACMG call, the same criteria and the same SpliceAI, PhyloP, gnomAD, REVEL and ClinVar fields for each of the 673,660 variants (8 November 2026, `fastvep 0.5.0-rc.1` of `data/benchmark/manuscript_bins_v38/fastvep_v38`).
