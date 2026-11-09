# What the complete databases changed

The benchmark was first built with SpliceAI and PhyloP databases that lack records which the gnomAD file they were distilled from carries (`results/complete_databases/README.md`).
Every run of the paper was then made again with databases built from the whole of that file, by the same classifiers on the same inputs, and the runs on the first databases were kept.
This file sets the two side by side.
It is written by `scripts/21_compare_databases.py` from the per-variant tables of the runs, and the numbers are in `results/complete_vs_april_databases/databases_*.csv`.
REVEL, ClinVar, ClinGen, the ClinVar protein index and RepeatMasker are the same files under both stacks.

## What the first databases lacked

| Source | Variants with a value under the first databases | Under the complete databases | Only under the complete | Only under the first | Differing where both have one |
|:---|---:|---:|---:|---:|---:|
| SpliceAI | 486,580 | 529,992 | 43,412 | 0 | 0 |
| PhyloP | 532,531 | 580,460 | 47,929 | 0 | 0 |
| gnomAD | 537,038 | 537,038 | 0 | 0 | 0 |
| REVEL | 347,292 | 347,292 | 0 | 0 | 0 |

The table is over the 673,278 variants of the April set that have an alternative allele (the other 382 are reference-agreement records, which get no call).
The complete SpliceAI database has a value for 529,992 of them and the first for 486,580, so the first lacks 8.2% of the records that exist; for PhyloP the shares are 580,460 and 532,531, and 8.3%.
No variant has a value only under the first databases, and where both have one the values are the same, so the complete databases add records and change none.
gnomAD and REVEL have the same variants with a value, and the same values, under both stacks.

## Where the missing records are

The extraction of April read gnomAD in batches of 200 regions of `clinvar_2star_regions.bed`, a chromosome at a time, ignored the exit status of `tabix` and sent its messages to /dev/null.
The 43,412 SpliceAI records that only the complete database holds lie in 14 stretches, each inside one batch and reaching its end: in each of those batches the first database has values up to a point and none after it, and the complete database has values after it.
For PhyloP the same is true of 47,472 of the 47,929 records, in 14 batches.
The other 457 are all indels, and the second base (position + 1) of each lies outside every region of the extraction.
The batches are the same for the two sources (14 batches lose records in both of 14), and the first lost position is the same in 14 of them, as it would be if one stream per batch fed both files, which is how the extraction script works (one `tabix` stream is read once and written to both).
The pattern is what a stream that ends early gives, and it is not a record of a stream that did: the extraction was not repeated, and the script threw its messages away.
`databases_lost_blocks.csv` lists each batch, with the variants of each source that have a value before the break and those lost after it.

## The April set

| Measure | v37, first databases | v37, complete | v38, first databases | v38, complete |
|:---|---:|---:|---:|---:|
| Exact agreement | 425,869 | 431,555 | 426,022 | 431,727 |
| Agreement in direction | 522,886 | 532,688 | 524,380 | 534,335 |
| Opposite-direction calls | 49 | 52 | 48 | 51 |
| Reference-discordant benign calls | 24 | 27 | 24 | 27 |
| Reference-discordant pathogenic calls | 25 | 25 | 24 | 24 |
| Recall of ClinVar pathogenic-side variants | 62.5% | 62.5% | 62.5% | 62.5% |
| Recall of ClinVar benign-side variants | 62.3% | 65.7% | 62.9% | 66.4% |
| Recall of ClinVar VUS | 97.2% | 97.2% | 97.2% | 97.1% |
| BP4 fires on | 359,427 | 378,930 | 359,436 | 378,939 |
| BP7 fires on | 132,615 | 144,948 | 132,615 | 144,948 |
| PP3 fires on | 68,012 | 69,133 | 68,012 | 69,133 |
| PVS1 fires on | 61,341 | 61,342 | 61,334 | 61,335 |
| PVS1 fires on canonical splice variants | 11,541 | 11,542 | 11,541 | 11,542 |

There are 13,449 canonical splice donor and acceptor variants in the set.
Under v38 the complete databases give +5,705 exact calls (+0.85 points) and +9,955 calls in direction (+1.48 points), and the opposite-direction calls go from 48 to 51.
Recall of ClinVar benign-side variants goes from 62.9% to 66.4% and recall of pathogenic-side variants from 62.5% to 62.5%.
Under v37 the complete databases give +5,686 exact calls (+0.84 points) and +9,802 calls in direction (+1.46 points), and the opposite-direction calls go from 49 to 52.
Recall of ClinVar benign-side variants goes from 62.3% to 65.7% and recall of pathogenic-side variants from 62.5% to 62.5%.

### The opposite-direction calls that differ

Under v38, 3 variants enter the opposite-direction calls, 0 leave them and 0 stay in them with another call.

- TMEM216 11:61392563 G>T (ClinVar Pathogenic, upstream_gene_variant): VUS (BS1) becomes Likely_benign (BS1&BP4); SpliceAI none to 0.07, PhyloP none to -0.668.
- VWF 12:6022888 G>A (ClinVar Pathogenic, synonymous_variant): VUS (no criterion) becomes Likely_benign (BP4&BP7); SpliceAI none to 0.00, PhyloP none to -0.105.
- GTPBP3 19:17339046 GGGGA>G (ClinVar Pathogenic, intron_variant): VUS (PM2_Supporting) becomes Likely_benign (PM2_Supporting&BP4&BP7); SpliceAI none to 0.00, PhyloP none to 0.312.

Under v37, 3 variants enter the opposite-direction calls, 0 leave them and 0 stay in them with another call.

- TMEM216 11:61392563 G>T (ClinVar Pathogenic, upstream_gene_variant): VUS (BS1) becomes Likely_benign (BS1&BP4); SpliceAI none to 0.07, PhyloP none to -0.668.
- VWF 12:6022888 G>A (ClinVar Pathogenic, synonymous_variant): VUS (no criterion) becomes Likely_benign (BP4&BP7); SpliceAI none to 0.00, PhyloP none to -0.105.
- GTPBP3 19:17339046 GGGGA>G (ClinVar Pathogenic, intron_variant): VUS (PM2_Supporting) becomes Likely_benign (PM2_Supporting&BP4&BP7); SpliceAI none to 0.00, PhyloP none to 0.312.

## The October release

The first design kept the 772,585 variants of the release that a gnomAD extract covered, and the complete databases need no such rule, so the release is scored whole (781,135 variants at two stars or better).
The rows below set the two stacks side by side over the variants that both hold, and the last rows give the complete stack on the whole release.

| Set | Measure | v37, first | v37, complete | v38, first | v38, complete |
|:---|:---|---:|---:|---:|---:|
| October release, the variants of the first design | Variants | 772,585 | 772,585 | 772,585 | 772,585 |
| October release, the variants of the first design | Agreement in direction | 591,576 | 604,246 | 593,343 | 606,193 |
| October release, the variants of the first design | Opposite-direction calls | 64 | 68 | 63 | 67 |
| October release, new variants of the first design | Variants | 107,528 | 107,528 | 107,528 | 107,528 |
| October release, new variants of the first design | Agreement in direction | 75,348 | 78,228 | 75,607 | 78,512 |
| October release, new variants of the first design | Opposite-direction calls | 18 | 19 | 18 | 19 |
| October release, all variants | Variants | - | 781,135 | - | 781,135 |
| October release, all variants | Agreement in direction | - | 609,923 | - | 611,945 |
| October release, all variants | Opposite-direction calls | - | 69 | - | 68 |

## The sweeps and the effects of the family

The 12 sweeps of the paper keep their shape but for 4, and 0 of them differs in its knee (BS1's is at 0.001 on both stacks, the far boundary of BP7's is at 300 on both stacks).

- BS1 bar: settings 0.0005 and 0.002 change order
- BS2 carrier bar, dominant-only genes: settings 0 and 0.001; 0 and 0.003 change order
- BS2 carrier floor, bar off: settings 5 and 1000 change order
- BP7 PhyloP bar: the cost of the step from 2.0 to 3 goes from zero to positive

Of the 245 effects of the switches, sources, hooks and gap runs on 7 measures, 7 changed sign, appeared or vanished, or changed by more than 20%, and 6 of those are of fewer than 10 calls.

- v38_sw_no_lowpen_gate false_benign: +1 on the first databases, +2 on the complete ones (fewer than 10 calls)
- v38_sw_no_lowpen_gate opposite: +1 on the first databases, +2 on the complete ones (fewer than 10 calls)
- v38_sw_no_homology_list same_pathogenic: +0 on the first databases, +1 on the complete ones (fewer than 10 calls)
- v38_src_nospliceai opposite: -13 on the first databases, -16 on the complete ones
- v38_src_nognomad false_benign: -2 on the first databases, -3 on the complete ones (fewer than 10 calls)
- v38_hk_bp7_walker false_benign: +2 on the first databases, +3 on the complete ones (fewer than 10 calls)
- v38_hk_bp7_walker opposite: +2 on the first databases, +3 on the complete ones (fewer than 10 calls)

## What this file does not say

The first databases were not rebuilt, and the run of the extraction that left them short was not repeated, so the cause above is read from the pattern of the records and not observed.
The earlier numbered versions v1 to v36 were not repeated on the complete databases, and their results are on the first ones.
A comparison of two stacks says what the missing records changed in these runs, and it does not say which stack is closer to the truth, which only ClinVar's classifications can be compared with.
