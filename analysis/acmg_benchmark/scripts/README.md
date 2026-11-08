# Real-data ACMG benchmark — scripts

End-to-end concordance harness against ClinVar 2-star+ variants. Runs the
**actual classifier** through `fastvep annotate --acmg` against
real allele- and gene-level annotations — no simulation, no priors.

The runner (`run_benchmark.sh`) and the analysis scripts are in this directory.
Every input and every bulk output is under `data/benchmark/`, which git ignores.

## Pipeline

```bash
# 1. Extract ClinVar's 2-star+ classifications from a ClinVar VCF (one-time)
python3 01_extract_clinvar_2star.py \
    /path/to/clinvar.vcf.gz \
    ../../../data/benchmark/

# 2. Build the SA databases (one-time per source).
#    See sa_sources/ for the build scripts:
#      - build_gnomad_per_chrom.sh                (gnomAD v4 exomes per chrom)
#      - build_spliceai_phylop.sh                 (distilled from gnomAD INFO)
#      - clingen_gdv_to_oga.py                    (ClinGen Gene-Disease Validity)
#    Plus the generic `fastvep sa-build` for ClinVar / ClinVar-protein /
#    REVEL / gnomAD gene constraints.

# 3. Annotate + score concordance
bash run_benchmark.sh v1

# 4. Render figures
python3 generate_figures.py
```

## Files in this directory

| File | Purpose |
|------|---------|
| `01_extract_clinvar_2star.py` | Filter ClinVar `clinvar.vcf.gz` → 2-star+ VCF + parallel TSV of ClinVar's classifications (`clinvar_2star_labels.tsv`, formerly `clinvar_2star_truth.tsv`) |
| `03_evaluate_concordance.py`  | Stream the bgzipped VCF output of `fastvep annotate --acmg` and emit concordance matrix / per-chrom / per-consequence / per-criterion CSVs |
| `04_build_md_review_table.py` | Build the review table of the medical geneticist's first round: the opposite-direction calls with the ClinVar record and the annotation fields they need |
| `05_append_current_calls.py`  | Append fastVEP's current call to a reviewer's own annotated workbook, leaving her notes, colours and row order as they are |
| `06_build_round5_table.py`, `07_build_round6_table.py` | Build the round-5 and round-6 review workbooks: one sheet, every call that still contradicts ClinVar 2-star+, with the reviewer's earlier verdicts carried onto the rows they belong to |
| `08_diff_calls.py`            | Diff two runs' ACMG calls variant by variant into a `call_changes_vN_to_vM.tsv`, reusing `03`'s own `pick_csq` so both describe the same transcript |
| `09_compare_runs.py`          | Compare runs with a baseline call by call: the headline counts, the two error types, and every changed call grouped by the criteria it lost and gained |
| `10_replay_bs2_guards.py`     | Replay a BS2 widening from the annotated VCFs of two runs and report what a guard on the new BS2 would have kept (limitation 20 in `docs/ACMG.md`) |
| `11_manuscript_runs.py`       | Run the threshold sweeps, the single-switch and leave-one-source-out ablations, the PM2 interaction under the 2015 table, the runs of two departures from the sources (BP1's veto, BP7), the runs that change an input (ClinGen, Orphanet) and (with a patched build) the point-cutoff runs behind the ACMG manuscript, each as a variant of `v37` through `run_benchmark.sh`, and with `--holdout` the six proposed changes on the October 2026 release; writes `results/v37_manuscript/` |
| `12_per_variant_table.py`     | Join a run's annotated VCF to ClinVar's labels, one row per variant (call, criteria, consequence, REVEL, SpliceAI, PhyloP, gnomAD counts, constraint, ClinGen tags), reusing `03`'s `pick_csq` |
| `13_manuscript_stats.py`      | Every table the ACMG manuscript reports about v37 and v38, written to `manuscript/acmg/data/`, including the offline replay of the combining rules and its validation against real runs |
| `14_temporal_holdout.py`      | Score v37 on a later ClinVar release and on the variants that were not in the benchmark set, and (`blockers`) say why BA1 and BS1 did not fire on its high-frequency Benign variants that were not called benign-side |
| `15_clingen_complete_export.py` | Put back the ClinGen curations that the export of 30 April 2026 did not list, from a current export, so that the `gap` runs of `11_manuscript_runs.py` can read the complete file |
| `16_compare_runs.py`          | Compare a rerun with the runs it repeats, file by file and record by record (the annotated VCFs by record, the per-variant tables by column name), and write `results/v37_manuscript/rerun_comparison.csv` |
| `17_joint_boundary_cases.py`  | Annotate, one at a time and through both databases with `--output-format json`, the variants that gnomAD's joint release turns into an error or out of a pathogenic-side call (the json carries the filtering allele frequency and the flags, which the VCF does not), and count the variants that lose PM2 by the reason, into `manuscript/acmg/data/v38_joint_errors.csv` and `v38_joint_pm2_lost.csv` |
| `18_criteria_status.py`       | Annotate every 40th record of the benchmark with `--output-format json`, which reports for each criterion whether it could be evaluated and whether it was met, and write the status of all 28 criteria (what each needs, whether a sites-only input can supply it, and how often the engine evaluated and met it) into `manuscript/acmg/data/criteria_status.csv`; it stops without writing if the engine disagrees with the groups it states. Needs `FASTVEP_V38_BIN` |
| `19_time_annotation.py`       | Time the annotation of the benchmark (`annotate --acmg --pick` over the 673,660 variants, piped to `bgzip`) on every core in rounds and on one thread, with its peak memory, and the scoring after it; writes `annotation_time.tsv` |
| `run_benchmark.sh`            | Annotate the ClinVar 2-star+ set with `--acmg --pick` and score it (`run_benchmark.sh vNN`); the flags, the hashes of the binary and the configurations, and the inputs it read are written to `run_provenance.txt`, and an annotated VCF is reused only when the same binary, flags, configuration and inputs made it |
| `sweep_acmg_thresholds.py`    | Run the benchmark once per value of one numeric configuration key and report what each value costs and buys |
| `build_bs2_review_table.py`, `make_bs2_review_xlsx.py` | Build the BS2 gene-disease review table and its workbook |
| `sa_sources/`                 | Download the sources and build the supplementary databases: `download_sa_sources.sh`, `build_gnomad_per_chrom.sh` (exomes) and `build_gnomad_release_per_chrom.sh` (exomes, genomes or joint), `build_spliceai_phylop.sh`, `clingen_gdv_to_oga.py`, `repeatmasker_to_bed.py`, `make_clingen_parse_fixture.py` (the fixture of real ClinGen rows that the classifier's parsing test reads) |
| `build_manuscript_binaries.sh` | Build the six binaries the manuscript's runs are made by (v37, its two patched builds, the commit before v38, v38 and its alternatives build), each from a clean export in a target directory of its own, and print the environment variables the run scripts read them from |
| `generate_figures.py`         | Read the v7 outputs at `data/benchmark/output_v7/` and emit 6 PDF + PNG panels (incl. v1 vs v7 comparisons) |
| `build_vcep_thresholds.py`    | Read published ClinGen VCEP BA1/BS1/PM2 thresholds out of the CSpec Registry API into `../data/vcep_thresholds.toml`, plus an audit table carrying the source sentence behind every number and the reason behind every rejection |
| `build_gene_disease_attributes.py` | Join Orphanet prevalence, age of onset and inheritance to those thresholds into `../data/gene_disease_attributes.tsv`, one row per gene-disorder pair |

## Data sources

| Source | Used by | Approx size |
|--------|---------|-------------|
| Ensembl GFF3 v115 | transcript models | 1.2 GB |
| Ensembl primary-assembly FASTA | amino acid prediction (PS1/PM5/PM1/PP3/BP4) | 3.2 GB |
| ClinVar `clinvar.vcf.gz` | PS4, PP5, BP6 (allele-level) | 190 MB |
| ClinVar `variant_summary.txt.gz` | **PS1 / PM1 / PM5** (protein-position index) | 438 MB |
| gnomAD v4.1 constraint metrics TSV | PVS1, PP2, BP1 (gene constraints) | 95 MB |
| gnomAD v4.1 exomes per-chrom VCF | **PM2 / BA1 / BS1 / BS2** + (distilled) PhyloP & SpliceAI | ~12 GB raw / chrom |
| REVEL `revel-v1.3_all_chromosomes.zip` | **PP3 / BP4 missense** | 637 MB |
| ClinGen Gene-Disease Validity CSV | PVS1 disease-gene fallback (preferred per Abou Tayoun 2018; OMIM `genemap2.txt` accepted as legacy) | 1 MB |
| ClinGen CSpec Registry (JSON-LD API) | published per-gene **BA1 / BS1 / PM2** thresholds | 10 MB, fetched |
| Orphanet `en_product6` / `en_product9_prev` / `en_product9_ages` | disorder prevalence, onset and inheritance per gene | 46 MB |

The classifier degrades gracefully — criteria with missing inputs are
marked `evaluated: false` rather than firing on noisy data.

## Bugs surfaced by this benchmark

The end-to-end real-data run surfaced several correctness issues that
unit tests had not caught. Each is fixed with a regression test.

### Classifier-side (`crates/fastvep-*`)

1. **`gnomad_genes` parser produced 0 records on gnomAD v4.x**. Old
   parser knew only v2.1 column names (`pLI`, `oe_lof_upper`, `mis_z`,
   `syn_z`); v4 uses dotted namespaces (`lof.pLI`, `lof.oe_ci.upper`,
   `mis.z_score`, `syn.z_score`). Fix: schema auto-detect; canonical /
   MANE-select preferred when v4 emits one row per transcript. Test:
   `test_parse_gnomad_gene_scores_v41_format`.
2. **`clinvar_protein` parser produced 0 records on ClinVar VCF**.
   ClinVar's `MC` field is just an SO term — no `p.` notation. Fix:
   parser auto-detects format and now accepts `variant_summary.txt.gz`,
   the ClinVar tab dump that exposes full HGVS Names with the
   `(p.Asp1692His)` block. Test:
   `test_parse_clinvar_protein_variant_summary_format`.
3. **GFF3 `.gz` input silently produced 0 transcripts**. `parse_gff3`
   was not gzip-aware. Fix: detect `.gz` / `.bgz` and wrap in
   `MultiGzDecoder` in three call sites; empty transcript output now a
   hard error.
4. **PM2 fired on every variant when gnomAD wasn't loaded**.
   `evaluate_pm2` returned `met=true, evaluated=true` when
   `input.gnomad` was `None`, conflating "absent in gnomAD" with "not
   loaded". Fix: returns `met=false, evaluated=false` when no record;
   PM2 only fires on a true `AC=0, AF=0` record. Tests:
   `test_pm2_no_gnomad_data_not_evaluated`,
   `test_pm2_truly_absent_with_gnomad_record_fires`.
5. **gnomAD v4 `mid` and `remaining` populations not parsed**.
   `parse_gnomad_vcf` and `GnomadData::max_pop_af` covered v2 codes
   only. Fix: added both.
6. **SpliceAI camelCase mismatch**. SaWriter wrote json_key `spliceAI`
   but classifier matched `spliceai | spliceAi | splice_ai`. Fix: added
   `spliceAI` to the match arms.
7. **PhyloP routing**. PhyloP attached at allele-level
   (`aa.supplementary`) but classifier read it from variant-level
   (`variant_supplementary`). Fix: read from both.

### Concordance-script-side (`03_evaluate_concordance.py`)

8. **Indel allele matching**. VCF stores REF=`C`/ALT=`CT` for an
   insertion; VEP CSQ Allele uses the inserted-base form (`T`) and `-`
   for deletions. The original script matched raw VCF ALT against CSQ
   Allele and so missed all 48,539 indels in the ClinVar 2-star+ set (counted as
   NoCall). Fix: `vep_allele(ref, alt)` strips the leading common
   prefix.

See [`../RUN_VERSIONS.md`](../RUN_VERSIONS.md) for the per-run impact of each fix.
