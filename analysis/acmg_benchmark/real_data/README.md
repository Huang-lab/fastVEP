# Real-data ACMG benchmark

End-to-end concordance harness against ClinVar 2-star+ variants. This complements the Monte-Carlo simulation in `../clinvar_concordance.py` by running the **actual classifier** through `fastvep annotate --acmg` against real allele- and gene-level annotations.

## What you need

| Source | Used by | Approx size |
|--------|---------|-------------|
| Ensembl GFF3 v115 | transcript models | 50 MB compressed |
| Ensembl primary-assembly FASTA | amino acid prediction (PS1/PM5/PM1/PP3/BP4) | 880 MB compressed |
| ClinVar `clinvar.vcf.gz` | PS4, PP5, BP6 (allele-level) | 190 MB |
| ClinVar `variant_summary.txt.gz` | **PS1 / PM1 / PM5** (protein-position index) | 110 MB |
| gnomAD v4.1 constraint metrics TSV | PVS1, PP2, BP1 (gene constraints) | 95 MB |
| gnomAD v4.x exomes/genomes per-chrom VCF | **PM2 / BA1 / BS1 / BS2** (allele-level AF) | 12 GB / chr (exomes) |
| REVEL `revel-v1.3_all_chromosomes.zip` | **PP3 / BP4 missense** | 3 GB |
| SpliceAI masked.snv VCF (optional) | PP3/BP4 splice + BP7 | 10 GB |

The classifier degrades gracefully — criteria with missing inputs are marked `evaluated: false` rather than firing on noisy data.

## Pipeline

```bash
# 1. Extract 2-star+ ClinVar truth set (~680k variants of P/LP/VUS/LB/B)
python3 01_extract_clinvar_2star.py /path/to/clinvar.vcf.gz .

# 2. Build SA databases (one-time per source — see docs/ACMG_SETUP.md)
fastvep sa-build --source clinvar          -i clinvar.vcf.gz             -o sa_db/clinvar
fastvep sa-build --source clinvar_protein  -i variant_summary.txt.gz     -o sa_db/clinvar_protein
fastvep sa-build --source gnomad_genes     -i gnomad.v4.1.constraint_metrics.tsv -o sa_db/gnomad_genes
fastvep sa-build --source gnomad           -i gnomad.exomes.v4.1.sites.chr<N>.vcf.bgz -o sa_db/gnomad_chr<N>
fastvep sa-build --source revel            -i revel_with_transcript_ids  -o sa_db/revel

# 3. Run the classifier on the truth VCF
bash 02_run_classifier.sh

# 4. Compute concordance vs ClinVar truth
python3 03_evaluate_concordance.py \
    --truth clinvar_2star_truth.tsv \
    --predictions clinvar_2star.fastvep.json \
    --out output/
```

## Bugs surfaced by this benchmark (and fixed in PR #fix/acmg-bugs-from-real-clinvar-benchmark)

The first end-to-end run on real data — even a chr17-only subset of 47k variants — surfaced four real correctness issues that none of the unit tests had caught. Each is now fixed with a focused regression test.

### 1. `gnomad_genes` parser produced 0 records on gnomAD v4.x

The parser accepted only the v2.1 column naming (`pLI`, `oe_lof_upper`, `mis_z`, `syn_z`). gnomAD v4.x uses dotted namespaces — `lof.pLI`, `lof.oe_ci.upper`, `mis.z_score`, `syn.z_score` — so feeding the latest constraint-metrics TSV silently produced an empty `.oga`. PVS1 then had no constraint data, PP2 always missed, BP1 always missed.

**Fix**: parser now resolves both schemas and prefers the canonical / MANE-select transcript when the TSV emits one row per transcript. New regression test: `test_parse_gnomad_gene_scores_v41_format`.

### 2. `clinvar_protein` parser produced 0 records on ClinVar VCF

ClinVar's `MC` field is just an SO term (`SO:0001583|missense_variant`) — it does not encode `p.` notation. `CLNHGVS` is genomic. So the parser, which scanned `MC` and `CLNHGVS` for `p.` blocks, never matched. PS1/PM1/PM5 were perpetually `evaluated: false`.

**Fix**: parser auto-detects format from the header line and now also accepts `variant_summary.txt.gz` — the canonical ClinVar tab dump that exposes full HGVS Names with the `(p.Asp1692His)` block. New regression test: `test_parse_clinvar_protein_variant_summary_format`.

### 3. GFF3 `.gz` input silently produced 0 transcripts

`File::open(gff3_path)` was passed to `parse_gff3` without gzip detection. With a gzipped GFF3 the parser saw raw deflate bytes, found no GFF lines, and returned an empty Vec. The annotator then proceeded with zero transcripts — every variant was annotated as intergenic, and `--acmg` produced criteria scored against an empty transcript context.

**Fix**: detect `.gz` / `.bgz` extensions in three call sites (annotate-context loader, CLI batch path, cache-build) and wrap in `MultiGzDecoder`. Empty transcript output is now treated as a hard error rather than a silent success.

### 4. PM2 fired on every variant when gnomAD wasn't loaded

`evaluate_pm2` returned `met: true, evaluated: true` whenever `input.gnomad` was `None`, with summary "Absent from gnomAD (no record)". This conflated two genuinely different states:

- The variant exists in a gnomAD record showing AC=0 (truly absent — PM2 should fire)
- gnomAD wasn't loaded for this run, or this variant's chromosome isn't covered by the loaded subset (PM2 cannot be evaluated)

The bug pushed PM2_Supporting onto every variant in the chr17 benchmark, which then triggered the SVI PVS+PP rule for every null variant and inflated the pathogenic call rate. The pre-fix concordance run showed PM2_Supporting firing for **100%** of B/LB/VUS records — clearly impossible.

**Fix**: when no gnomAD record is present at all, PM2 returns `met: false, evaluated: false` with a summary explaining the data is missing. The user has to load a gnomAD `.osa` for PM2 to fire. Updated unit tests (`test_pm2_no_gnomad_data_not_evaluated`, `test_pm2_truly_absent_with_gnomad_record_fires`) to lock in the corrected semantics. The two `lib.rs` integration tests that depended on the old behavior were re-anchored to use a truly-absent gnomAD record (AC=0, AF=0).

## Initial chr17 benchmark snapshot (post-fixes)

47,591 ClinVar 2-star+ variants on chr17. **No** allele-level gnomAD database loaded (so PM2 / BA1 / BS1 / BS2 cannot fire), **no** REVEL (so PP3/BP4 missense cannot fire). What's left: PVS1, PS1, PM1, PM5, PP2, BP1, and PM4.

| truth | n | exact | same-direction | opposite | no_call |
|-------|---|-------|----------------|----------|---------|
| Pathogenic | 4,775 | 0 | 859 (18%) | 0 | 0 |
| Likely_pathogenic | 888 | 219 | 219 (25%) | 0 | 0 |
| VUS | 18,015 | 18,011 | 18,011 (>99%) | 0 | 0 |
| Likely_benign | 8,022 | 0 | 0 | 0 | 0 |
| Benign | 10,539 | 0 | 0 | 0 | 0 |

What this tells us:

- The classifier is now working **correctly** for the criteria it has data for. Per-class criterion firing rates show good specificity (e.g. PVS1 fires 23.5% of true Pathogenic, 0% of true Benign).
- The 0% benign call rate is a direct consequence of missing gnomAD allele-level data (BA1, BS1, BS2) and missing REVEL (BP4). When those sources are loaded, BA1 alone should immediately reclassify the bulk of Benign variants.
- The pathogenic same-direction rate of 18% is constrained by missing REVEL (~70% of pathogenic missense need PP3 to reach LP).

The benchmark provides a clear yardstick: this run + REVEL + gnomAD allele-level should land in the 60–80% same-direction range for pathogenic and ≥80% for benign, matching the published ACMG concordance literature.
