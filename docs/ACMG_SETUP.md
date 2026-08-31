# ACMG setup: building the reference annotation stack

This guide builds the exact set of supplementary annotation databases the ACMG-AMP classifier reads.
It is the stack the ClinVar concordance benchmark has run against since v15, when the last of the nine sources was added, so following it end to end reproduces [`analysis/acmg_benchmark/`](../analysis/acmg_benchmark/)'s published behaviour rather than approximating it.

For what the criteria mean and how to retune their thresholds, see [ACMG.md](ACMG.md).
For the on-disk format of the databases you are about to build, see [SUPPLEMENTARY_ANNOTATIONS.md](SUPPLEMENTARY_ANNOTATIONS.md).

**Every source in this stack is public, and no registration, licence or account is required for any of it.**

Two upstream files that normally do gate an ACMG setup are not used here, and it is worth being explicit about where their data comes from instead:

- **SpliceAI comes out of the gnomAD sites VCF, not from Illumina BaseSpace.** gnomAD v4.1 publishes Illumina's score in its own INFO column as `spliceai_ds_max`, so the same gnomAD fetch that gives you allele frequencies also gives you the splice predictions. No Illumina account is involved anywhere in this guide. The BaseSpace file is listed as an alternative in [its section](#spliceai-and-phylop-both-distilled-from-the-gnomad-vcf-allele-level), together with exactly what the gnomAD version does and does not carry.
- **The disease-gene source is ClinGen Gene-Disease Validity, not OMIM.** ClinGen SVI prefers it on the merits, and it downloads with a plain `curl`. Real OMIM `genemap2.txt` still parses if you have a licence.

PhyloP comes out of the same gnomAD INFO column (`phylop`), which also removes the 5.2 GB UCSC download.
That one carries a caveat worth reading before you rely on it, because gnomAD's phyloP is not UCSC's phyloP.

> **The failure mode to design against.** An incomplete `--sa-dir` does not produce an error.
> `fastvep annotate --acmg` runs, classifies every variant, and silently omits the evidence it has no data for.
> A missing REVEL database does not fail the run; it turns pathogenic missense variants into VUS.
> That is why this guide ends with [a verification step](#verify-the-stack-answers) rather than with a build command, and why `sa-build` succeeding is not the thing to check.

## What the classifier actually reads

Nine sources, in three scopes.
The "silent effect if absent" column is the one that matters: fastVEP marks a criterion `evaluated: false` rather than firing it on missing data, so a gap here narrows the evidence rather than announcing itself.

| Source | `--source` | Scope | Drives | Silent effect if absent |
|---|---|---|---|---|
| gnomAD v4.1 sites | `gnomad` | allele | BA1, BS1, BS2, PM2 | No frequency evidence at all: nothing is ever called too common to be pathogenic, and PM2 falls back to ClinVar's frequencies |
| ClinVar | `clinvar` | allele | PM2 frequency backstop, PM3/BP2 companion lookup, PP5/BP6 | PM2 loses its backstop, so a common variant absent from gnomAD reads as absent from the population |
| REVEL v1.3 | `revel` | allele | PP3, BP4 (missense) | No computational evidence on any missense variant |
| SpliceAI, **from gnomAD's `spliceai_ds_max`** | `spliceai` | allele | PVS1 splice branch, PP3, BP4 (splice), BP7 | No splice prediction: PVS1 cannot check a predicted-intact site, BP7 cannot clear a synonymous variant |
| PhyloP, **from gnomAD's `phylop`** (Zoonomia) | `phylop` | allele | BP7 conservation tier | BP7 stops distinguishing conserved from unconserved positions |
| ClinGen Gene-Disease Validity | `omim` | gene | PVS1 gate, PP2/PM1 gate | The gene-validity gate never blocks, so PVS1 can fire in genes ClinGen curated and refuted |
| gnomAD gene constraints | `gnomad_genes` | gene | PVS1 (pLI/LOEUF), PP2, BP1 | PVS1 loses its constraint evidence and PP2/BP1 stop evaluating |
| ClinVar protein + splice index | `clinvar_protein` | gene | PS1, PM1, PM5, BP1 | Four criteria stop evaluating entirely |
| RepeatMasker | `custom_bed` | interval | BP3 | BP3 reports "cannot tell" for every in-frame indel |

Three of ClinVar's consumers are off unless you switch them on.
PP5 and BP6 need `use_pp5_bp6 = true`, and PS4 is `NotEvaluated` by default because ClinVar review stars are not the case-control statistics PS4 asks for; `use_clinvar_stars_as_ps4_proxy = true` opts into the older proxy behaviour.
Both defaults follow ClinGen SVI.
ClinVar still earns its place in the stack without them, through the PM2 backstop, the PM3/BP2 companion lookup, and PS1's comparison variant.

Two of those rows are less obvious than they look.

`omim` is the json_key, not the source.
The canonical content is **ClinGen Gene-Disease Validity**, which ClinGen SVI and Abou Tayoun 2018 prefer over OMIM: it is a multi-curator scored rubric that explicitly classifies Limited, Disputed and Refuted relationships instead of listing every reported association.
The key kept its historical name.
Real OMIM `genemap2.txt` still parses into the same `.oga` and remains supported.

`clinvar_protein` is a second, differently shaped read of ClinVar, not a duplicate of the allele database.
It holds two indices: classified missense keyed by protein position, and canonical splice-dinucleotide variants keyed by genomic position.
It must be built from `variant_summary.txt.gz`, not from the allele VCF, for reasons given in [its section](#clinvar-protein-and-splice-index-gene-level).

### What is deliberately not in the stack

`dbnsfp` and `gerp` build and load, and their values appear in the JSON `details` block, but **no criterion fires on them**.
An earlier PP3/BP4 path took a 3-of-4 vote across SIFT, PolyPhen, PhyloP and GERP; it was removed per Pejaver 2022, whose recommendation is a single calibrated tool, and that tool is REVEL.
SIFT and PolyPhen are surfaced for a reviewer to read and nothing more.
Build them if you want the predictions in your output, skip them if you do not.
`alphamissense`, `primateai`, `dbsnp`, `cosmic`, `topmed` and `onekg` are likewise supported by `sa-build` and are not read by the classifier.

## Before the databases: gene models and reference sequence

The classifier is downstream of consequence prediction, so a gap here is a gap in everything.

```bash
mkdir -p data && cd data

# Gene models. Ensembl 115, GRCh38.
curl -LO https://ftp.ensembl.org/pub/release-115/gff3/homo_sapiens/Homo_sapiens.GRCh38.115.gff3.gz
gunzip Homo_sapiens.GRCh38.115.gff3.gz

# Reference sequence, and its index. Without the FASTA there are no coding
# sequences, so there are no missense calls, so PS1, PM5, PP2, BP1 and the
# whole missense half of PP3/BP4 have nothing to act on.
curl -LO https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa
```

`--acmg` turns `--hgvs` on for you.
Several criteria read the HGVS `c.` string rather than raw coordinates: PVS1's canonical +/-1,2 splice gate and BP7's deep-intronic extension both need the intronic offset, and BA1's Ghosh 2018 exception list is keyed by `c.` notation.
Passing `--hgvs` explicitly changes nothing.

## Build the stack

Every command below is a two-step pattern: download the upstream release, then convert it with `fastvep sa-build`.
**`sa-build` is a converter, not a downloader.**
Handed a missing or truncated file it can succeed and write a database containing nothing, which then annotates nothing, which then classifies everything as VUS.
Check each download before building it.

```bash
REPO=/path/to/fastVEP          # this checkout, for the helper scripts below
mkdir -p sa_databases sa_sources && cd sa_sources
```

### ClinVar (allele level)

Download the NCBI release specifically.
It carries `AF_EXAC`, `AF_TGP` and `AF_ESP`, which become the PM2 frequency backstop consulted when gnomAD has no matching record.
A stripped or third-party ClinVar VCF may not have them, and PM2 then treats a common variant as absent from the population.

```bash
curl -LO https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
curl -LO https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.md5

# Verify before building. A truncated download builds a partial database
# without erroring; this file once arrived ending mid-chr4.
gzip -t clinvar.vcf.gz && echo "gzip OK"
md5 -q clinvar.vcf.gz            # macOS; use md5sum elsewhere. Compare to the .md5
zcat clinvar.vcf.gz | grep -m1 '^##INFO=<ID=AF_EXAC'   # the backstop field must be present

fastvep sa-build --source clinvar -i clinvar.vcf.gz \
    -o ../sa_databases/clinvar --assembly GRCh38
```

ClinVar is re-released weekly, so this is the one source worth refreshing on a schedule.
For a clinical setup, pin a dated release instead of tracking the moving file, so a re-run reproduces its calls:

```bash
# https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/weekly/ lists them by date
curl -LO https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/weekly/clinvar_20260822.vcf.gz
```

Rebuild `clinvar_protein.oga` from the matching `variant_summary.txt.gz` whenever you refresh, so the two ClinVar views agree with each other.

### ClinVar protein and splice index (gene level)

Build this from `variant_summary.txt.gz`.
The allele VCF also parses, but it is much the weaker input: its `MC` field is a bare Sequence Ontology term rather than a protein change, and `CLNHGVS` is genomic, so it carries no HGVS `c.` token at all.
From the VCF you get few protein entries and **no splice index whatsoever**, which silently removes PS1's splice path.

```bash
curl -LO https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz
gzip -t variant_summary.txt.gz

fastvep sa-build --source clinvar_protein -i variant_summary.txt.gz \
    -o ../sa_databases/clinvar_protein --assembly GRCh38
```

`--assembly` selects which rows feed the splice index and is stamped into every gene record.
Protein positions are assembly-independent and come from all rows; genomic splice positions are not.

### gnomAD sites (allele level)

This is the one source whose full download is a real cost.
Measured on 2026-08-31 across chr1-22, X and Y:

| Release | Total download |
|---|---|
| `gnomad.exomes.v4.1.sites` | **198.4 GB** |
| `gnomad.genomes.v4.1.sites` | **563.0 GB** |

You almost certainly do not need all of it.
gnomAD's VCFs are tabix-indexed on a public bucket, so you can pull only the regions your variants occupy and never download the rest.
For a panel or exome workflow that is a few hundred megabytes instead of a few hundred gigabytes.

```bash
# Region-restricted (recommended). targets.bed is your capture target, panel,
# or the merged regions of the VCFs you plan to annotate.
BASE=https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/exomes/gnomad.exomes.v4.1.sites
for chr in {1..22} X Y; do
  url="${BASE}.chr${chr}.vcf.bgz"
  tabix -h "$url" "chr${chr}:1-1000" 2>/dev/null | grep '^#' > gnomad_chr${chr}.vcf
  # gnomAD is chr-prefixed; accept a BED written either way and normalise.
  awk -v c="chr${chr}" -v n="${chr}" '$1==c || $1==n {print c":"$2+1"-"$3}' targets.bed \
    | xargs -n 200 tabix "$url" >> gnomad_chr${chr}.vcf
  bgzip -f gnomad_chr${chr}.vcf
  fastvep sa-build --source gnomad -i gnomad_chr${chr}.vcf.gz \
      -o ../sa_databases/gnomad_chr${chr} --assembly GRCh38
  rm -f gnomad_chr${chr}.vcf.gz
done
```

To build whole-genome coverage instead, download each `.vcf.bgz` in full and run the same `sa-build` on it.
Use `.../vcf/genomes/gnomad.genomes.v4.1.sites.chr${chr}.vcf.bgz` for the genomes release.

Every `.osa`/`.osa2` in `--sa-dir` is loaded, and the JSON key comes from the source type rather than the filename, so `gnomad_chr1` and `gnomad_chr2` both register under `gnomad`.
Per-chromosome files are the easier shape: they make the build resumable and let you add a chromosome later.

### gnomAD gene constraints (gene level)

```bash
curl -LO https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv

fastvep sa-build --source gnomad_genes -i gnomad.v4.1.constraint_metrics.tsv \
    -o ../sa_databases/gnomad_genes --assembly GRCh38
```

### ClinGen Gene-Disease Validity (gene level)

The gate this feeds blocks only on positive evidence against a gene, never on absence.
ClinGen has curated roughly 3,000 genes out of a genome, so "not in the file" means "nobody has reached it yet" and must not suppress PVS1.
Run v10 of the benchmark measured what the opposite policy costs: blocking on absence stripped PVS1 from 1,497 truth-pathogenic null variants in genes like SPAST, ABCB11, FLG and LAMB3, all of which cause disease and none of which ClinGen has curated.

The converter therefore emits every classification, including Limited, Disputed, Refuted and No Known Disease Relationship, so the classifier can tell "curated and found wanting" apart from "not curated".

```bash
curl -L -o clingen_gene_validity.csv https://search.clinicalgenome.org/kb/gene-validity/download

# Two positional arguments: input CSV, output TSV. Not stdin/stdout.
python3 "$REPO/analysis/acmg_benchmark/scripts/sa_sources/clingen_gdv_to_oga.py" \
    clingen_gene_validity.csv clingen_gdv.tsv

fastvep sa-build --source omim -i clingen_gdv.tsv \
    -o ../sa_databases/omim --assembly GRCh38
```

Expect roughly 3,000 genes and an `omim.oga` under 1 MB.
The run on 2026-08-31 produced 3,030.

**Using real OMIM instead.** `genemap2.txt` from <https://www.omim.org/downloads> populates the same `.oga` and the same json_key; pass it straight to `--source omim`.
Access requires a registered account, and an unauthenticated request to that page is refused, so it cannot be scripted the way the ClinGen download can.

### REVEL (allele level)

```bash
curl -LO https://rothsj06.dmz.hpc.mssm.edu/revel-v1.3_all_chromosomes.zip
unzip revel-v1.3_all_chromosomes.zip
```

The archive is 667 MB and holds exactly one member, `revel_with_transcript_ids`, which unpacks to 6.50 GB.
It is a **CSV**, not a TSV, with this header:

```
chr,hg19_pos,grch38_pos,ref,alt,aaref,aaalt,REVEL,Ensembl_transcriptid
```

fastVEP reads column 2, `grch38_pos`.
Rows whose GRCh38 position is empty (liftover failures) are skipped rather than misplaced.

```bash
fastvep sa-build --source revel -i revel_with_transcript_ids \
    -o ../sa_databases/revel --assembly GRCh38
```

Splitting the CSV per chromosome first and building 24 databases bounds peak memory, which is what the benchmark does.
The result is identical either way.

### SpliceAI and PhyloP, both distilled from the gnomAD VCF (allele level)

**There is no separate SpliceAI download in this stack.**
gnomAD v4.1 publishes Illumina's SpliceAI score in its own INFO column, under `spliceai_ds_max`, and describes it in the VCF header as "Illumina's SpliceAI max delta score".
It publishes phyloP alongside it under `phylop`.
So the tabix fetch you already ran for allele frequencies is the same fetch that produces the splice and conservation databases, and the Illumina BaseSpace file never enters the picture.

Three databases, one source.
This is also what removes the 5.2 GB UCSC phyloP download from the setup.

```bash
for chr in {1..22} X Y; do
  python3 "$REPO/scripts/extract_gnomad_scores.py" --chrom $chr --regions targets.bed --out-dir .
  fastvep sa-build --source spliceai -i spliceai_chr${chr}.vcf \
      -o ../sa_databases/spliceai_chr${chr} --assembly GRCh38
  fastvep sa-build --source phylop -i phylop_chr${chr}.tsv \
      -o ../sa_databases/phylop_chr${chr} --assembly GRCh38
  rm -f spliceai_chr${chr}.vcf phylop_chr${chr}.tsv
done
```

Three things about this path are worth knowing before you rely on it.

**The SpliceAI loss does not reach the classification.**
gnomAD stores only the maximum delta score, so the extractor writes that max into all four DS fields.
Every ACMG consumer of SpliceAI in fastVEP reads `max_delta_score()` and nothing else, in PVS1's splice branch, PP3, BP4 and BP7, so the maximum of four copies of the maximum is the same number the real file would give.
What is lost is the four DP position fields, which have no gnomAD equivalent and are written as `0`.
The `SpliceAI` INFO field fastVEP emits therefore carries placeholder positions; do not read them downstream.

**gnomAD's phyloP is a different score from UCSC's.**
gnomAD ships the Zoonomia 241-placental-mammal score.
UCSC's `phyloP100way` is a 100-way vertebrate alignment.
They are not interchangeable, and fastVEP's `phylop_conserved` default of 2.0 was measured against the Zoonomia score, so this is the one the shipped threshold expects.

**Coverage is gnomAD's coverage.**
The exomes release covers capture targets and their flanks, so a deep-intronic variant has no record there and BP7's deep-intronic extension has nothing to read.
Use the genomes release if that matters to you, or the upstream files below.

**Upstream alternatives.**
Illumina's precomputed SpliceAI (all four components, real DP positions) is at <https://basespace.illumina.com/s/otSPW8hnhaZR> and needs an Illumina account:

```bash
fastvep sa-build --source spliceai -i spliceai_scores.masked.snv.hg38.vcf.gz \
    -o ../sa_databases/spliceai --assembly GRCh38
```

UCSC's phyloP is public, one file per chromosome, 5.2 GB in total:

```bash
UCSC=https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP100way/hg38.100way.phyloP100way
for c in {1..22} X Y M; do
  curl -LO "$UCSC/chr${c}.phyloP100way.wigFix.gz"
done

# Concatenate IN CHROMOSOME ORDER into one multi-member gzip stream. fastVEP
# reads gzip with MultiGzDecoder, so concatenating compressed members is safe.
# `cat chr*.wigFix.gz` sorts lexicographically (chr1, chr10, chr11, ..., chr2),
# and sa-build streams a positional source straight into the index, so a
# lexicographic glob builds happily and then aborts with a "not sorted" error.
: > hg38.phyloP100way.wigFix.gz
for c in {1..22} X Y M; do cat "chr${c}.phyloP100way.wigFix.gz" >> hg38.phyloP100way.wigFix.gz; done

fastvep sa-build --source phylop -i hg38.phyloP100way.wigFix.gz \
    -o ../sa_databases/phylop --assembly GRCh38
```

If you take this route, re-check `phylop_conserved` against your own data before trusting BP7: the threshold shipped with fastVEP was measured on the Zoonomia score.

### RepeatMasker (interval level)

```bash
curl -LO https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz

python3 "$REPO/analysis/acmg_benchmark/scripts/sa_sources/repeatmasker_to_bed.py" \
    rmsk.txt.gz > repeatmasker.bed

fastvep sa-build --source custom_bed --name repeatmasker \
    -i repeatmasker.bed -o ../sa_databases/repeatmasker --assembly GRCh38
```

`--name repeatmasker` is load-bearing.
The classifier locates this track by looking for a supplementary key containing `repeat`, so a database built under any other name loads and is then ignored.

The converter reports what it kept, and on the 2026-08-31 release that is:

```
kept 5,317,291 intervals; skipped 366,399 on non-primary contigs, 0 by class filter, 0 malformed
```

About 240 MB as a `.osi`.

Without it BP3 reports `evaluated: false` for every in-frame indel rather than assuming the variant is not in a repeat.

## Verify the stack answers

`ls` on the output directory tells you the files exist, which is not the question.
The question is whether each source answers a query, and the only way to know is to annotate something and look.

```bash
python3 scripts/check_acmg_stack.py \
    --gff3 data/Homo_sapiens.GRCh38.115.gff3 \
    --fasta data/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    --sa-dir sa_databases/
```

It annotates eight real ClinVar variants spanning the consequence classes the criteria branch on, then reports which of the nine sources produced an annotation and which criteria never evaluated.
It exits non-zero when a source answers nothing, so it can gate a deployment.

```
source                 key                    status    evidence
gnomAD                 gnomad                 OK        7/8 annotated                 BA1 BS1 BS2 PM2
ClinVar                clinvar                OK        8/8 annotated                 PM2 backstop, PP5/BP6, PM3/BP2
REVEL                  revel                  OK        3/8 annotated                 PP3 BP4 (missense)
SpliceAI               spliceAI               OK        6/8 annotated                 PVS1 splice, PP3 BP4 (splice), BP7
PhyloP                 phylop                 OK        6/8 annotated                 BP7 conservation tier
ClinGen GDV / OMIM     omim                   OK        7/8 variants' genes           PVS1 gate, PP2/PM1 gate
gnomAD constraint      gnomad_genes           OK        8/8 variants' genes           PVS1 pLI/LOEUF, PP2, BP1
ClinVar protein index  clinvar_protein        OK        8/8 variants' genes           PS1 PM1 PM5 BP1
RepeatMasker           (BP3 evaluated)        OK        BP3 evaluated 8/8             BP3

Complete: every source in the reference stack answered at least one query.
```

Partial counts are expected and are not failures.
REVEL scores only missense positions, not every variant is in gnomAD, and a region-restricted build answers only inside its regions.
A **zero** is the signal.

Add `--fastvep target/release/fastvep` if the binary is not on your `PATH`, and `--acmg-config` to check the stack against your own thresholds.

RepeatMasker is checked indirectly, through BP3's `evaluated` flag.
An interval database yields an annotation only when the position falls inside an interval, so no key appears for a variant that simply is not in a repeat, and that is indistinguishable from an unloaded file by looking at the annotation alone.
BP3 carries the distinction instead.

When something reads zero, the three usual causes are:

1. **Contig naming.** `chr1` in one file and `1` in the other. This is the most common one by a distance.
2. **Assembly.** A GRCh37 source against GRCh38 input produces no matches and no error.
3. **An empty build.** Check the source file, not the database.

Criteria reported as never evaluated because no supplementary annotation supplies them are working as designed: PS2, PM6, PP1 and BS4 need trio or family data (`--proband`, `--mother`, `--father`), PS3 and BS3 need `--functional-evidence`, and PS4 is `NotEvaluated` by default per ClinGen SVI because review stars are not case-control statistics.

## Run it

```bash
fastvep annotate \
  -i your_variants.vcf \
  -o annotated.vcf \
  --gff3 data/Homo_sapiens.GRCh38.115.gff3 \
  --fasta data/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  --sa-dir sa_databases/ \
  --acmg --pick \
  --output-format vcf
```

Use `--output-format json` while you are still verifying the setup.
It carries the full ACMG block for every criterion, including the ones that did not fire and why:

```json
{
  "acmg": {
    "classification": "LikelyPathogenic",
    "shorthand": "LP",
    "triggered_rule": "PVS + >=1 PP (SVI)",
    "criteria": [
      { "code": "PM2_Supporting", "met": true, "evaluated": true, "summary": "Absent in gnomAD ..." }
    ]
  }
}
```

`vcf` is the format to use once you are past setup.
It is the only format that carries per-transcript ACMG and ACMG_CRITERIA fields in a stable column position, and it is roughly 40x smaller than pretty-printed JSON on a large cohort.

For trio analysis, add `--proband`, `--mother` and `--father` with the sample names, which turns on PS2 and PM6.
For curated experimental evidence, `--functional-evidence` takes a TSV of chrom, pos, ref, alt, criterion, strength, pmid, note, and supplies PS3/BS3.
An entry there also suppresses PP3/BP4/BP7 for that variant, because ClinGen SVI ranks experimental evidence above computational prediction.

## Sizes, measured

Re-derive this table at any time rather than trusting it, since upstream releases grow:

```bash
bash benchmarks/check_urls.sh sizes
```

It reads each `Content-Length` with a HEAD request and downloads nothing.
The figures below are what it printed on 2026-08-31:

| File | Size |
|---|---|
| `clinvar.vcf.gz` | 193.4 MB |
| `variant_summary.txt.gz` | 442.5 MB |
| `gnomad.v4.1.constraint_metrics.tsv` | 95.5 MB |
| ClinGen gene-validity CSV | 1.1 MB |
| `revel-v1.3_all_chromosomes.zip` | 667.1 MB (6.50 GB unpacked) |
| `rmsk.txt.gz` | 155.6 MB |
| gnomAD v4.1 exomes, chr1-22 X Y | 198.4 GB (full release; a region-restricted extract is a fraction of it) |
| gnomAD v4.1 genomes, chr1-22 X Y | 563.0 GB (same) |
| SpliceAI | 0, it comes out of the gnomAD fetch above |
| PhyloP | 0, same |
| UCSC phyloP100way, 25 files | 5.2 GB, only if you choose it over gnomAD's |
| Illumina SpliceAI (BaseSpace) | not used; account required |

Built, from the benchmark's own `sa_db`, totalling 1.3 GB:

| Database | Built size | Coverage |
|---|---|---|
| `clinvar.osa` | 50 MB | whole release |
| `clinvar_protein.oga` | 18 MB | 15,665 genes |
| `gnomad_chr*.osa2` | 369 MB | region-restricted |
| `gnomad_genes.oga` | 1.5 MB | 18,173 genes |
| `omim.oga` | 0.4 MB | 2,924 genes |
| `revel_chr*.osa` | 421 MB | genome-wide |
| `spliceai_chr*.osa` | 151 MB | region-restricted |
| `phylop_chr*.osa` | 115 MB | region-restricted |
| `repeatmasker.osi` | 239 MB | whole assembly |

The three region-restricted rows are scoped to the benchmark's ClinVar 2-star+ regions.
Yours will scale with your own target file.

## Keeping the URLs honest

Every download URL in this guide is checked by:

```bash
bash benchmarks/check_urls.sh acmg     # liveness, 13 URLs
bash benchmarks/check_urls.sh sizes    # download sizes
```

Run the first before a fresh build.
Upstream hosts do move, and this guide has already outlived two moves: `hgdownload.cse.ucsc.edu` stopped resolving entirely, taking both the RepeatMasker table and phyloP with it to `hgdownload.soe.ucsc.edu`, and Ensembl rotated the fly assembly from BDGP6.46 to BDGP6.54 and the rat from mRatBN7.2 to GRCr8 between releases.
A URL that 404s is the good case, because you notice.
The dangerous one is the guide that still resolves and quietly hands you a stale assembly.

Two URLs are deliberately outside that check.
`https://www.omim.org/downloads` refuses an unauthenticated request, so it returns 403 whether the page is healthy or not, and there is nothing useful to assert about it.
The dbNSFP and GERP sources are not read by the classifier, so they are not part of this stack.

## Configuration

Per-criterion thresholds, gene-specific mechanism overrides and trio settings live in one TOML file.
See [ACMG.md, Configuration File](ACMG.md#configuration-file) for the full schema.

```bash
fastvep annotate ... --acmg --acmg-config my_thresholds.toml
```
