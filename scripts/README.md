# scripts/

Operational helpers that are not part of the build.
Benchmark-specific scripts live in [`analysis/acmg_benchmark/scripts/`](../analysis/acmg_benchmark/scripts/) and [`benchmarks/`](../benchmarks/) instead.

| Script | What it does |
|---|---|
| `check_acmg_stack.py` | Annotates eight real ClinVar variants and reports which supplementary annotation sources actually answered. Exits non-zero when one answers nothing. |
| `extract_gnomad_scores.py` | Distils SpliceAI and PhyloP databases out of the gnomAD v4.1 sites VCF, so neither needs its own download. |
| `deploy-minimal.sh` | Human GRCh38 gene models and reference sequence only. |
| `deploy-clinical.sh` | The clinical deployment: human plus ClinVar and dbSNP. |
| `deploy-full.sh` | Every organism and every supplementary source. |
| `release.sh` | Cuts a tagged release. |

## Setting up for ACMG classification

Read [`docs/ACMG_SETUP.md`](../docs/ACMG_SETUP.md).
It builds the nine-source stack the classifier reads and ends with `check_acmg_stack.py`.

The reason that verification step exists, rather than a size check on the built files, is that an incomplete `--sa-dir` is not an error condition.
`fastvep annotate --acmg` runs to completion against a half-built stack and classifies every variant; the criteria with no data behind them report `evaluated: false` and drop out of the evidence.
A missing REVEL database does not fail the run, it turns pathogenic missense variants into VUS.
So the thing worth checking is not whether `sa-build` succeeded but whether each source answers a query:

```bash
python3 scripts/check_acmg_stack.py \
    --gff3 data/Homo_sapiens.GRCh38.115.gff3 \
    --fasta data/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    --sa-dir sa_databases/
```

RepeatMasker is the one source checked indirectly, through BP3's `evaluated` flag.
An interval database only yields an annotation when the position falls inside an interval, so a variant that simply is not in a repeat looks exactly like an unloaded file from the annotation alone.
BP3 carries that distinction and the script reads it there.

## SpliceAI and PhyloP come from gnomAD

`extract_gnomad_scores.py` exists because gnomAD v4.1 already publishes both scores in its INFO column, as `spliceai_ds_max` and `phylop`.
Taking them from there removes the Illumina BaseSpace account that SpliceAI otherwise requires and the 5.2 GB UCSC phyloP download, and the same tabix fetch that supplies allele frequencies supplies all three databases.

```bash
python3 scripts/extract_gnomad_scores.py --chrom 17 --regions targets.bed --out-dir extracts/
```

The script's docstring records exactly what that costs.
In short: gnomAD stores only SpliceAI's maximum delta score, which is the only field any ACMG criterion reads, so the classification is unaffected; the four DP position fields have no gnomAD equivalent and are written as placeholders.
gnomAD's phyloP is the Zoonomia 241-mammal score rather than UCSC's 100-way vertebrate score, and fastVEP's `phylop_conserved` default was measured against the Zoonomia one.

## Checking upstream before a download

```bash
bash benchmarks/check_urls.sh acmg     # every URL in the ACMG guide
bash benchmarks/check_urls.sh speed    # every URL in the speed benchmark
bash benchmarks/check_urls.sh sizes    # re-derive the guide's download sizes
```
