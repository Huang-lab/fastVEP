# Moving from `vep` to `fastvep`

fastVEP is a standalone annotator.
Its consequence and HGVS engine is a port of Ensembl VEP's, and [VEP_DIVERGENCE.md](VEP_DIVERGENCE.md) measures where the two agree.
It is not a drop-in for the `vep` command line.
This page says what carries over from a VEP installation, what does not, and how to check your own sites before you switch.
It describes fastVEP with the whole of the #139 series merged: the known-variant reader, the tab columns, the protein version and the RefSeq names.

Each statement below was checked against fastVEP's code, the output of its `--help`, the file it names, the source of Ensembl VEP (release/116), or a run of Ensembl VEP 116.3, the program, on the same input (Docker `ensemblorg/ensembl-vep:release_116.3`).
Where something was not checked, it says so.

## What fastVEP reads from a VEP installation

| In the VEP cache directory | Read? | How |
|---|---|---|
| `{chr}/all_vars.gz` with its `.csi` or `.tbi`, `info.txt`, `chr_synonyms.txt` | yes | `--cache-dir`, for `Existing_variation` only, and only in `fastvep annotate` |
| `{chr}/{start}-{end}.gz` transcript slices (Perl `Storable`) | **no** | gene models come from GFF3 |
| `{chr}/{start}-{end}_reg.gz` regulatory slices | **no** | |
| `{chr}/{start}-{end}_var.gz`, the unconverted variation slices | **no**, and the run stops with an error | convert with VEP's `convert_cache.pl`, or use the indexed cache download |
| plugin data files | no | |

So the transcript models are always rebuilt from a GFF3, Ensembl's or NCBI's:

```bash
fastvep annotate -i in.vcf --gff3 genes.gff3 --fasta genome.fa --cache-dir ~/.vep/homo_sapiens_refseq/116_GRCh38 ...
# or compile the GFF3 once and reuse it
fastvep cache --gff3 genes.gff3 --fasta genome.fa -o genes.cache
fastvep annotate -i in.vcf --transcript-cache genes.cache ...
```

`--cache-dir` takes the VEP variation cache as it is.
VEP's converter indexes it with `tabix -C`, so the index is a `.csi`; fastVEP reads that and the older `.tbi`.
The contig names of the input are matched through the `chr_synonyms.txt` in the directory, if it has one, and the `chr`/`MT` rules.
A contig the cache does not hold is counted, and the run ends with a warning that names it.

## Flags

fastVEP spells a multi-word flag with a hyphen (`--pick-allele`), where VEP uses an underscore (`--pick_allele`).
The underscore spelling is not accepted.

| VEP | fastVEP | Note |
|---|---|---|
| `--cache --offline --dir_cache DIR` | `--cache-dir DIR` | known variants only, see above; there is no database mode |
| `--fasta FILE` | `--fasta FILE` | |
| `--tab`, `--vcf`, `--json` | `--output-format tab`, `vcf`, `json` | `vcf` is the default |
| `--check_existing` | implied by `--cache-dir` | |
| `--refseq`, `--merged` | `--gff3` | pass the NCBI GFF3, or both, with `LABEL=path`; the `SOURCE` field of the CSQ says which |
| `--hgvs` | `--hgvs` | adds `HGVSc`, `HGVSp` and `HGVS_OFFSET` |
| `--hgvsg` | `--hgvsg` | written as VEP writes it, an indel 3'-shifted over the genome and a repeat a `dup`, which reads the reference, so give `--fasta` |
| `--symbol`, `--biotype`, `--canonical`, `--numbers`, `--variant_class`, `--show_ref_allele` | `--symbol`, `--biotype`, `--canonical`, `--numbers`, `--variant-class`, `--show-ref-allele` | tab columns, in VEP's order; the VCF `CSQ` has most of them always |
| `--everything` | `--everything` | `--hgvs --symbol --numbers --canonical --biotype --variant-class`; VEP's set has more |
| `--pick`, `--pick_allele`, `--pick_allele_gene`, `--flag_pick`, `--flag_pick_allele`, `--flag_pick_allele_gene`, `--pick_order` | the same, with hyphens | |
| `--distance` | `--distance` | |

Not accepted, with nothing equivalent: `--af`, `--af_1kg`, `--af_gnomad` and the other `--af_*` and `--max_af`; `--exclude_predicted`; `--use_transcript_ref` and `--bam`; `--total_length`; `--fields`; `--no_escape`; `--plugin`; `--custom`; `--assembly` and `--species`; `--force_overwrite` (the `-o` file is overwritten as it is).
fastVEP rejects each of them with an error, so a command line that uses one fails at once.

Population frequencies come from supplementary databases you build with `sa-build` (gnomAD, 1000 Genomes, TOPMed), and appear as `FV_*` fields and columns, not as VEP's `AF` columns.

### A worked example

A VEP command line of the kind a production pipeline has:

```bash
vep --cache --offline --force_overwrite --tab \
    --show_ref_allele --variant_class --biotype --af_1kg --check_existing \
    --hgvs --hgvsg --refseq --exclude_predicted --no_escape \
    --use_transcript_ref --numbers --total_length \
    --fasta $ref --dir_cache $dir_cache --assembly GRCh38 --bam $bam \
    --plugin REVEL,... --custom file=phylop.bw,format=bigwig
```

and what carries over to fastVEP:

```bash
fastvep annotate -i in.vcf -o out.tab --output-format tab \
    --gff3 refseq.gff3 --fasta $ref --cache-dir $dir_cache \
    --show-ref-allele --variant-class --biotype --hgvs --hgvsg --numbers
```

What is dropped, and why: `--af_1kg` (no VEP-style frequency columns), `--exclude_predicted` (filter the GFF3), `--no_escape` and `--force_overwrite` (not needed), `--use_transcript_ref` and `--bam` (no BAM-edited transcripts), `--total_length` (not implemented), `--assembly` (the assembly is whatever the GFF3 and FASTA are), `--plugin` and `--custom` (see below).
`--refseq` becomes the choice of GFF3, and `--check_existing` is implied by `--cache-dir`.

## Plugins and custom tracks

There is no plugin API and no `--custom`.
`fastvep sa-build` is the extension point.
Its dedicated sources, from the code, are `gnomad`, `onekg` (`1000g`), `topmed`, `alphamissense`, `dbsnp`, `cosmic`, `clinvar`, `revel`, `primateai`, `dbnsfp`, `phylop`, `gerp`, `dann`, `spliceai` and `mitomap`, plus the gene-level `omim`, `gnomad_genes` and `clinvar_protein`, and `custom_vcf` and `custom_bed` for any VCF or BED.
See [SA_DATABASES.md](SA_DATABASES.md) for how each is built.

For the data sets people usually run VEP plugins for:

| Data set | In fastVEP |
|---|---|
| AlphaMissense | dedicated source `alphamissense` |
| dbNSFP | dedicated source `dbnsfp` (which of its columns it keeps was not checked here) |
| REVEL | dedicated source `revel` |
| SpliceAI | dedicated source `spliceai` |
| PhyloP | dedicated source `phylop`, built from wigFix files such as UCSC's; bigWig is not read directly |
| dbscSNV, Condel | no dedicated source; whether they can be converted to a `custom_vcf` was not checked |

A plugin you wrote yourself (`RefSeqBamRepair`, `RefSeqIntronAnchorFix`, `Selective3PrimeConsequence`, `SpliceMotifVariant`) has no equivalent hook.

## RefSeq

Use the NCBI GFF3 for the assembly (the whole-genome `genomic.gff` of `GCF_000001405.40_GRCh38.p14`, for instance).
Its sequence names are accessions (`NC_000010.11`) and the FASTA's are `chr10` or `10`, so build a transcript cache with a synonyms file and pass it to `annotate`:

```bash
fastvep cache --gff3 GCF_000001405.40_GRCh38.p14_genomic.gff --fasta GRCh38.fa \
  --synonyms chr_synonyms.txt -o refseq.cache
fastvep annotate -i in.vcf --transcript-cache refseq.cache --fasta GRCh38.fa ...
```

The synonyms file is in the format VEP reads from a cache directory's `chr_synonyms.txt`, one line of equivalent names per contig.
[README.md](../README.md#chromosome-naming-across-sources---synonyms) has the details.

Two things are not covered, and were not measured:

- VEP's RefSeq cache carries transcripts whose sequence differs from the genome, repaired from a BAM (`--bam`, `--use_transcript_ref`).
  fastVEP builds transcript sequence from the GFF3 and the FASTA and has no equivalent, so its agreement with VEP on those transcripts is unknown.
- There is no `--exclude_predicted`.
  Filter the GFF3 to the records you want before you build the cache.

fastVEP names what a RefSeq GFF3 holds as VEP's GFF reader does: a transcript is its `transcript_id` (`NM_000001.3`, not the ID `rna-NM_000001.3`), a gene is the `GeneID` of its Dbxref (`7157`, not `gene-TP53`), and a biotype the record does not state is guessed from its type and its `gbkey` (`protein_coding` for an mRNA, `ncRNA` for an lnc_RNA, `misc_RNA`, `rRNA`, `tRNA`, ...).
The protein in an `HGVSp` is the id NCBI gives, `NP_000001.2`.

[VEP_DIVERGENCE.md, Part 3](VEP_DIVERGENCE.md#a-refseq-gff3) compares the two tools on a RefSeq GFF3: chr21 and chr22 of GRCh38.p14 and the 527 HG002 variants on them, run through VEP 116.3 `--gff` and fastVEP.
135 of 81,510 cells differ, all of them the multi-allelic difference described under [What will differ](#what-will-differ).
What has not been compared is a run against VEP's RefSeq *cache*, which is what `--refseq` reads: VEP builds it from RefSeq by its own pipeline.

## The output

The tab file and the VCF `CSQ` follow the cells of VEP 116.3, which was run to check them: `Uploaded_variation` (the VCF line itself for a record with no ID), `Location`, `-` for an absent value, `Amino_acids` as a single letter where the residue is unchanged, `HGVSp` with `=` as `%3D`, and no feature on an intergenic row.
The columns are the ones you ask for, in VEP's order.
[README.md](../README.md#tab-output) lists the rules.
The columns VEP computes per allele and per known variant (`AF` and its populations, `CLIN_SIG`, `SOMATIC`, `PHENO`) are not written as VEP writes them.
The `CSQ` has fields with those names that summarise the whole variant, and its `AF` takes gnomAD first where VEP's is 1000 Genomes.

### What will differ

On a whole-genome callset, 20,241 records of GIAB HG002 annotated with the same GRCh38 Ensembl 115 GFF3 and FASTA by VEP 116.3 and by fastVEP, 3,876 of 5,581,650 cells of the tab files differ (0.07 %), and every column but five (`HGVS_OFFSET`, `HGVSc`, `HGVSp`, `HGVSg`, `cDNA_position`) is the same on every row.
[VEP_DIVERGENCE.md, Part 3](VEP_DIVERGENCE.md#part-3---the-tab-output-against-vep-1163) lists each class with its size and says which tool is right.
The ones you are likely to meet:

- **A record with several ALT alleles.**
  VEP 116 no longer splits it into one variant per allele, and then writes `HGVSc` of an insertion or deletion allele without the 3' shift and without `HGVS_OFFSET`, and `HGVSg` of an allele that trims to a substitution as a `delins` that does not reproduce the variant.
  fastVEP writes every allele as it would be alone.
  On HG002, 36 % of the rows of such records differ in `HGVSc`, and nothing else does because of it.
  A `*` ALT counts as another allele: `TC>T,*` is shifted by fastVEP and not by VEP.
- **A `*` allele.**
  VEP writes no row for it, and fastVEP writes one, with an empty `Consequence`.
- **A change that runs off the end of a transcript or an exon.**
  VEP writes `?-7` or `199-?` for its position, `INTRON` beside `EXON`, and an `HGVSc` of the part inside the transcript; fastVEP writes the one coordinate it has.
- **An indel written with bases the alleles share at the back.**
  VEP minimises it from both ends and fastVEP from the front only, so `Allele`, `Location` and `REF_ALLELE` differ, and a consequence term can.
  None of the 2,818 indels of HG002 is like that.
- **Rows only one tool writes.**
  fastVEP reads every transcript of the GFF3 and VEP's `--gff` reader keeps only some feature types (`unconfirmed_transcript` and `scRNA` are not among them), so a variant near one is `intergenic_variant` in VEP and an upstream or downstream variant in fastVEP: 127 rows of 186,182 on HG002, and 7 the other way.

- **A RefSeq GFF3.**
  `SOURCE` is `RefSeq` where VEP writes the name of the GFF3, and fastVEP writes rows for the mature microRNA records (`rna-MIR185`), which hang off a transcript and which VEP writes nothing for.
  Where VEP can name no biotype it drops the transcript, and fastVEP keeps it, named by its SO term.

A RefSeq run against VEP's RefSeq cache has not been compared, and so none of this is known for it.

## Before you switch

Run both tools on a sample of your own sites, with the flags you use, and compare row by row.
`validation/run_divergence.sh` and `validation/compare_rows.py` do that for the Ensembl GFF3 and write one file per disagreeing field, which is the quickest way to see whether a disagreement touches the fields you read.
