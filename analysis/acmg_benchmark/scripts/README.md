# Real-data ACMG benchmark — scripts

End-to-end concordance harness against ClinVar 2-star+ variants. Runs the
**actual classifier** through `fastvep annotate --acmg` against
real allele- and gene-level annotations — no simulation, no priors.

The orchestrating runner and all intermediate / final outputs live
under `data/benchmark/`. This directory holds only the analysis
scripts.

## Pipeline

```bash
# 1. Build the truth set from a ClinVar VCF (one-time)
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
| `01_extract_clinvar_2star.py` | Filter ClinVar `clinvar.vcf.gz` → 2-star+ truth VCF + parallel TSV |
| `03_evaluate_concordance.py`  | Stream the bgzipped VCF output of `fastvep annotate --acmg` and emit concordance matrix / per-chrom / per-consequence / per-criterion CSVs |
| `generate_figures.py`         | Read the v7 outputs at `data/benchmark/output_v7/` and emit 6 PDF + PNG panels (incl. v1 vs v7 comparisons) |
| `build_vcep_thresholds.py`    | Read published ClinGen VCEP BA1/BS1/PM2 bars out of the CSpec Registry API into `../data/vcep_thresholds.toml`, plus an audit table carrying the source sentence behind every number and the reason behind every rejection |
| `build_gene_disease_attributes.py` | Join Orphanet prevalence, age of onset and inheritance to those bars into `../data/gene_disease_attributes.tsv`, one row per gene-disorder pair |

## Data sources

| Source | Used by | Approx size |
|--------|---------|-------------|
| Ensembl GFF3 v115 | transcript models | 1.2 GB |
| Ensembl primary-assembly FASTA | amino acid prediction (PS1/PM5/PM1/PP3/BP4) | 3.2 GB |
| ClinVar `clinvar.vcf.gz` | PS4, PP5, BP6 (allele-level) | 190 MB |
| ClinVar `variant_summary.txt.gz` | **PS1 / PM1 / PM5** (protein-position index) | 110 MB |
| gnomAD v4.1 constraint metrics TSV | PVS1, PP2, BP1 (gene constraints) | 95 MB |
| gnomAD v4.1 exomes per-chrom VCF | **PM2 / BA1 / BS1 / BS2** + (distilled) PhyloP & SpliceAI | ~12 GB raw / chrom |
| REVEL `revel-v1.3_all_chromosomes.zip` | **PP3 / BP4 missense** | 637 MB |
| ClinGen Gene-Disease Validity CSV | PVS1 disease-gene fallback (preferred per Abou Tayoun 2018; OMIM `genemap2.txt` accepted as legacy) | 1 MB |
| ClinGen CSpec Registry (JSON-LD API) | published per-gene **BA1 / BS1 / PM2** bars | 10 MB, fetched |
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
   Allele and so missed all 48,539 indels in the truth (counted as
   NoCall). Fix: `vep_allele(ref, alt)` strips the leading common
   prefix.

See [`../RUN_VERSIONS.md`](../RUN_VERSIONS.md) for the per-run impact of each fix.
