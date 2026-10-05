# Supplementary annotation databases: the complete manifest

This is the list of every source `fastvep sa-build` understands, for both GRCh38 and GRCh37.
For each one it gives where the file comes from, what it costs to download, the build command, and whether it can be fetched without an account.
It answers the question "is there a pre-built package, or a recommended set with URLs and build commands".

**There is no pre-built package.**
The databases are built from the upstream releases by `fastvep sa-build`, and the biggest of those releases (gnomAD, dbSNP) are hundreds of gigabytes.
What exists instead is a script that does the download, the integrity check and the build for every source that is publicly downloadable: [`scripts/build-sa-databases.sh`](../scripts/build-sa-databases.sh).

For what each database emits, see [SUPPLEMENTARY_ANNOTATIONS.md](SUPPLEMENTARY_ANNOTATIONS.md).
For the ACMG subset with the reasoning behind each source, see [ACMG_SETUP.md](ACMG_SETUP.md).
This page does not repeat that reasoning.

## Start here

You do not need every source to test that supplementary annotation works.
The `small` tier is about 1 GB of downloads, takes a few minutes, and gives you ClinVar, ClinVar protein/splice, gnomAD gene constraints, ClinGen gene-disease validity and RepeatMasker:

```bash
scripts/build-sa-databases.sh --assembly GRCh38 --out sa_databases           # sources: small
scripts/build-sa-databases.sh --assembly GRCh37 --out sa_databases_grch37
```

Add the allele-level frequency and score sources one chromosome at a time with `--chroms`, which keeps the first run short:

```bash
scripts/build-sa-databases.sh --assembly GRCh38 --out sa_databases \
    --sources gnomad,revel,alphamissense --chroms 22
```

Then annotate with `--sa-dir sa_databases`.
Every `.osa2`, `.osa`, `.osi` and `.oga` in the directory is loaded, and a database's key in the output comes from its source rather than its filename, so `gnomad_chr1` and `gnomad_chr2` both feed `FV_GNOMAD`.

**Use one `--sa-dir` per assembly.**
fastVEP warns when one `--sa-dir` holds databases built for different assemblies, and names which databases are on which.
It cannot warn when the whole directory is for the wrong assembly, because a GFF3 does not declare one: a GRCh37 directory used on a GRCh38 VCF loads and answers with the data of a different base or nothing.
Keep `sa_databases/` and `sa_databases_grch37/` apart, as above.

## The manifest

"Script" means `build-sa-databases.sh` handles it.
Sizes were read from the servers' `Content-Length` on 2026-10-05, except the GRCh38 gnomAD totals, which are the 2026-08-31 measurement in [ACMG_SETUP.md](ACMG_SETUP.md).

| Source (`--source`) | Access | GRCh38 | GRCh37 | Script |
|---|---|---|---|---|
| ClinVar (`clinvar`) | open | NCBI `vcf_GRCh38/clinvar.vcf.gz`, 193 MB | NCBI `vcf_GRCh37/clinvar.vcf.gz`, 200 MB | `clinvar` |
| ClinVar protein and splice (`clinvar_protein`) | open | NCBI `variant_summary.txt.gz`, 438 MB. One file serves both assemblies; `--assembly` picks the splice rows | same file | `clinvar_protein` |
| gnomAD gene constraints (`gnomad_genes`) | open | v4.1 `constraint_metrics.tsv`, 95 MB | v2.1.1 `lof_metrics.by_gene.txt.bgz`, 4 MB | `gnomad_genes` |
| ClinGen gene-disease validity (`omim`) | open | `search.clinicalgenome.org/kb/gene-validity/download`, 1 MB. Gene-level, so the same file serves both | same file | `clingen` |
| RepeatMasker (`custom_bed --name repeatmasker`) | open | UCSC `hg38/database/rmsk.txt.gz`, 155 MB | UCSC `hg19/database/rmsk.txt.gz`, 148 MB | `repeatmasker` |
| gnomAD sites (`gnomad`) | open | v4.1. Exomes 198 GB, genomes 563 GB (chr1-22, X, Y) | v2.1.1. Exomes 63 GB, genomes 494 GB (chr1-22, X) | `gnomad` |
| 1000 Genomes (`onekg`) | open | NYGC 30x, 445 MB for chr22 | Phase 3, 205 MB for chr22 | `onekg` |
| AlphaMissense (`alphamissense`) | open | Zenodo 8208688 `AlphaMissense_hg38.tsv.gz`, 642 MB | `AlphaMissense_hg19.tsv.gz`, 622 MB | `alphamissense` |
| REVEL (`revel`) | open | `revel-v1.3_all_chromosomes.zip`, 667 MB. One file carries both assemblies | same file | `revel` |
| PhyloP (`phylop`) | open | UCSC `hg38.100way.phyloP100way`, 5.2 GB | UCSC `hg19.100way.phyloP100way`, 5.4 GB | `phylop` |
| dbSNP (`dbsnp`) | open | `GCF_000001405.40.gz`, 29.6 GB | `GCF_000001405.25.gz`, 28.2 GB | `dbsnp` |
| SpliceAI (`spliceai`) | GRCh38: open, derived from gnomAD. GRCh37: account | gnomAD v4.1 `spliceai_ds_max`, via `scripts/extract_gnomad_scores.py` | Illumina BaseSpace only (gnomAD v2.1.1 has no SpliceAI field) | no |
| COSMIC (`cosmic`) | account | `CosmicCodingMuts.vcf.gz` from cancer.sanger.ac.uk | the file for that build from the same page | no |
| PrimateAI (`primateai`) | account | the provider's `chr pos ref alt score` file for that build | same | no |
| TOPMed (`topmed`) | account | the TOPMed freeze VCF | not located | no |
| dbNSFP (`dbnsfp`) | not fetched | dbNSFP 4.x | the same file, read through its `hg19_*` columns | no |
| GERP, DANN (`gerp`, `dann`) | see note | a `chrom pos score` TSV you prepare | the same, GRCh37 coordinates | no |
| MitoMap (`mitomap`) | see note | a `pos ref alt disease status` TSV you prepare | same | no |

"Account" means a file I could not fetch without logging in or accepting terms, or could not locate at a stable public URL.
That is a statement about what could be fetched from a script, not about the provider's licence, so check the provider's terms before redistributing a database built from one.

**Notes on specific rows.**

- **gnomAD GRCh37 and GRCh38 are different releases, not the same release twice.** GRCh38 is v4.1 and GRCh37 is v2.1.1. Both parse, with the field naming detected from the header, but their per-population columns and the number of samples differ, so frequencies are not comparable across the two.
- **1000 Genomes GRCh38 is the NYGC 30x call set, GRCh37 is Phase 3.** The two spell their per-population frequency fields differently (`AF_EAS` against `EAS_AF`). Both are read. Before this page was written, the GRCh38 file built a database whose population columns were all empty, which looked exactly like a population with no data.
- **REVEL and dbNSFP carry both assemblies in one file.** `--assembly` chooses the coordinate column (`hg19_pos` or `grch38_pos` for REVEL; `hg19_chr`/`hg19_pos(1-based)` or `#chr`/`pos(1-based)` for dbNSFP). Before this page was written both always read the GRCh38 column, so a GRCh37 build indexed GRCh38 positions under a GRCh37 label and returned the score of a different base without any error.
- **chrX and chrY for 1000 Genomes.** The script builds chr1-22 only, because the chrX files are named differently on both servers. Fetch them by hand and build with `--source onekg`.
- **GERP and DANN** are read as `chrom pos score` (or BED-like `chrom start end score`) TSVs sorted by position. I did not verify an upstream distribution against the parser, so prepare the TSV from whichever release you use.
- **MitoMap** does not publish the TSV the parser reads, and mitomap.org returned HTTP 403 to scripted requests, so the file has to be assembled by hand from its tables.

## What was tested

Everything marked "open" was built from the real upstream file on 2026-10-05 on macOS, with the release binary and `scripts/build-sa-databases.sh`.
Where a whole download was not worth repeating, the cell says what was built instead.

| Source | GRCh38 | GRCh37 |
|---|---|---|
| ClinVar | whole file, 4,587,125 records | whole file, 4,587,491 records |
| ClinVar protein and splice | whole file, 15,904 genes | whole file, 15,904 genes |
| gnomAD gene constraints | v4.1, 18,173 genes | v2.1.1, 19,658 genes |
| ClinGen gene-disease validity | 3,057 genes | 3,057 genes |
| RepeatMasker | 5,317,291 intervals | 5,232,241 intervals |
| gnomAD sites | exomes v4.1, chr22 1 Mb region, 238,961 records | exomes v2.1.1, chr22 1 Mb region, 16,093 records |
| 1000 Genomes | chr22 whole, 1,066,557 records | chr22 whole, 1,110,240 records |
| AlphaMissense | whole file, 71,697,556 records | whole file, 69,716,655 records |
| REVEL | chr22, 1,776,286 records | chr22, 1,785,378 records |
| PhyloP | UCSC chr22, 36,218,018 records | UCSC chr22, 34,480,521 records |
| dbSNP | first 300,000 lines of the file, 368,117 records | first 300,000 lines of the file, 368,979 records |
| SpliceAI from gnomAD | chr22 slice, 55,428 records | not available |
| dbNSFP | unit tests only: no copy of the file was available | unit tests only |
| COSMIC, PrimateAI, TOPMed, GERP, DANN, MitoMap | not tested against an upstream file | not tested against an upstream file |

dbSNP is the one open source not built from the whole download.
It is 28-30 GB, and the script's `dbsnp` step downloads and builds it unchanged; the first 300,000 lines exercise the same parser and the same RefSeq contig mapping.

Building is not the same as answering, so each assembly was also annotated end to end with the databases the script produced:

- **GRCh37.** A chr22 VCF in GRCh37 coordinates (1,883 ClinVar records from a 1 Mb window) was annotated against the Ensembl GRCh37 GFF3, read gzipped, with the script's databases. Every source answered: ClinVar 1,883, PhyloP 1,883, gnomAD gene constraints 1,844, ClinVar protein 1,816, REVEL 1,315, gnomAD sites 1,137, AlphaMissense 1,078, 1000 Genomes 269, RepeatMasker 78. ClinGen did not answer in that window because none of its genes are curated; a probe of 20 ClinVar records in SLC25A1, which is, returned `FV_OMIM=SLC25A1|0|mitochondrial disease (ClinGen Definitive/AR ...)`.
- **GRCh38.** The 1,000 chr22 variants in `validation/human/chr22_1kgp.vcf` against the Ensembl 115 GFF3. 1000 Genomes answered 687, PhyloP 1,000, RepeatMasker 605, ClinVar protein and gnomAD gene constraints 402, AlphaMissense 3, REVEL 3, ClinVar 1. gnomAD sites answered none, because that database was deliberately built from a 1 Mb region the VCF does not occupy.
- **REVEL, both assemblies.** The same chr22 CSV was built twice. The GRCh38 database answers `22:16590880 G>C` with 0.136 and not `22:17071770 G>C`. The GRCh37 database answers `22:17071770 G>C` with 0.136 and not `22:16590880 G>C`. Those are the same REVEL row in its two coordinate systems.

## Check it answers

`sa-build` printing a record count proves the file parsed, not that the database matches your VCF.
Two failures leave a database that builds and then annotates nothing: an assembly mismatch, and a truncated download.
The script refuses a build of zero records, but a plausible-looking count from the wrong assembly passes that test.

So run one variant you know is in each source, with `--sa-only`, which skips consequence prediction and prints only the supplementary fields:

```bash
printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n22\t16590880\t.\tG\tC\t.\t.\t.\n' > probe.vcf
fastvep annotate -i probe.vcf -o probe.out.vcf --sa-only --sa-dir sa_databases --output-format vcf
grep -v '^##' probe.out.vcf
```

An empty INFO column means that database did not answer.
For the nine ACMG sources there is a ready-made version: `scripts/check_acmg_stack.py`, described in [scripts/README.md](../scripts/README.md).

## Per-source commands

The script is these commands in a loop.
Run them by hand when you want one source or a different release.
`sa-build` reads gzip by its magic bytes, so downloaded `.gz`/`.bgz` files go in as they are.

```bash
A=GRCh38   # or GRCh37

# ClinVar. Use the NCBI release: it carries AF_EXAC/AF_TGP/AF_ESP, the PM2 backstop.
curl -LO https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_$A/clinvar.vcf.gz
fastvep sa-build --source clinvar -i clinvar.vcf.gz -o sa_databases/clinvar --assembly $A

# ClinVar protein and splice index. Build it from variant_summary, not from the VCF.
curl -LO https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz
fastvep sa-build --source clinvar_protein -i variant_summary.txt.gz -o sa_databases/clinvar_protein --assembly $A

# gnomAD gene constraints.
curl -LO https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv       # GRCh38
curl -LO https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz # GRCh37
fastvep sa-build --source gnomad_genes -i <that file> -o sa_databases/gnomad_genes --assembly $A

# gnomAD sites, one chromosome at a time. GRCh37 contigs are bare (22), GRCh38 are chr22.
curl -LO https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/exomes/gnomad.exomes.v4.1.sites.chr22.vcf.bgz         # GRCh38
curl -LO https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.22.vcf.bgz      # GRCh37
fastvep sa-build --source gnomad -i <that file> -o sa_databases/gnomad_chr22 --assembly $A

# AlphaMissense. hg38 for GRCh38, hg19 for GRCh37.
curl -LO https://zenodo.org/records/8208688/files/AlphaMissense_hg38.tsv.gz
fastvep sa-build --source alphamissense -i AlphaMissense_hg38.tsv.gz -o sa_databases/alphamissense --assembly GRCh38

# REVEL. The archive holds one CSV; --assembly picks the coordinate column.
curl -LO https://rothsj06.dmz.hpc.mssm.edu/revel-v1.3_all_chromosomes.zip
unzip -p revel-v1.3_all_chromosomes.zip | awk -F, 'NR==1||$1==22' > revel_chr22.csv   # per chromosome bounds memory
fastvep sa-build --source revel -i revel_chr22.csv -o sa_databases/revel_chr22 --assembly $A

# PhyloP, per chromosome. Use hg19 in both places for GRCh37.
curl -LO https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP100way/hg38.100way.phyloP100way/chr22.phyloP100way.wigFix.gz
fastvep sa-build --source phylop -i chr22.phyloP100way.wigFix.gz -o sa_databases/phylop_chr22 --assembly GRCh38

# dbSNP. Contigs are RefSeq accessions (NC_000001.11); --assembly maps them.
curl -LO https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz            # GRCh38
curl -LO https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.25.gz            # GRCh37
fastvep sa-build --source dbsnp -i GCF_000001405.40.gz -o sa_databases/dbsnp --assembly GRCh38

# 1000 Genomes. Only columns 1-8 are used, so drop the genotypes on the way down.
curl -L https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr22.filtered.SNV_INDEL_SV_phased_panel.vcf.gz \
    | gzip -dc | cut -f1-8 > onekg_chr22.vcf                                                # GRCh38
curl -L https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz \
    | gzip -dc | cut -f1-8 > onekg_chr22.vcf                                                # GRCh37
fastvep sa-build --source onekg -i onekg_chr22.vcf -o sa_databases/onekg_chr22 --assembly $A

# RepeatMasker. --name repeatmasker is load-bearing: the classifier finds the track by it.
curl -LO https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz            # hg19 for GRCh37
python3 analysis/acmg_benchmark/scripts/sa_sources/repeatmasker_to_bed.py rmsk.txt.gz > repeatmasker.bed
fastvep sa-build --source custom_bed --name repeatmasker -i repeatmasker.bed -o sa_databases/repeatmasker --assembly $A

# ClinGen gene-disease validity (json_key omim).
curl -L -o clingen_gene_validity.csv https://search.clinicalgenome.org/kb/gene-validity/download
python3 analysis/acmg_benchmark/scripts/sa_sources/clingen_gdv_to_oga.py clingen_gene_validity.csv clingen_gdv.tsv
fastvep sa-build --source omim -i clingen_gdv.tsv -o sa_databases/omim --assembly $A
```

For the account-gated sources the command is the same shape once you have the file: `--source cosmic`, `primateai`, `topmed` or `dbnsfp`, with `-i` the file you downloaded and `--assembly` the build it is in.
SpliceAI for GRCh38 without an account is covered in [ACMG_SETUP.md](ACMG_SETUP.md#spliceai-and-phylop-both-distilled-from-the-gnomad-vcf-allele-level), including what the gnomAD-derived version does not carry.
