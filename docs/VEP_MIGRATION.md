# Moving from `vep` to `fastvep`

fastVEP is a standalone annotator.
Its consequence and HGVS engine is a port of Ensembl VEP's, and [VEP_DIVERGENCE.md](VEP_DIVERGENCE.md) measures where the two agree.
It is not a drop-in for the `vep` command line: there is no plugin API, no `--custom`, no database mode and no BAM-edited transcript, and gene models come from a GFF3 and not from a VEP cache.
What it does write is VEP's output: the tab cells and columns, the columns VEP computes from its variation cache (`AF` and the other frequencies, `MAX_AF`, `CLIN_SIG`, `SOMATIC`, `PHENO`, `PUBMED`), and most of the flags that ask for them.
This page says what carries over from a VEP installation, what does not, and how to check your own sites before you switch.
It describes fastVEP with the whole of the #139 series merged: the known-variant reader and columns, the tab columns, the protein version, the RefSeq names, the positions and terms at the edge of a transcript, the flags below, the changes that make a GFF3 run read a transcript as VEP's cache holds it, and what a review of the series found at the start codon, in the protein notation, in the VCF a run reads and in the speed of a run.

Each statement below was checked against fastVEP's code, the output of its `--help`, the file it names, the source of Ensembl VEP (release/116), or a run of Ensembl VEP 116.3, the program, on the same input (Docker `ensemblorg/ensembl-vep:release_116.3`).
Where something was not checked, it says so.

## What fastVEP reads from a VEP installation

| In the VEP cache directory | Read? | How |
|---|---|---|
| `{chr}/all_vars.gz` with its `.csi` or `.tbi`, `info.txt`, `chr_synonyms.txt` | yes | `--cache-dir`, for `Existing_variation` and for the columns VEP computes from the same records, and only in `fastvep annotate` |
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

fastVEP spells a multi-word flag with a hyphen (`--pick-allele`).
Every flag it shares with VEP is also accepted in VEP's underscore spelling (`--pick_allele`), so those can be pasted as they are.

| VEP | fastVEP | Note |
|---|---|---|
| `--cache --offline --dir_cache DIR` | `--cache-dir DIR` | known variants only, see above; there is no database mode |
| `--fasta FILE` | `--fasta FILE` | |
| `--tab`, `--vcf`, `--json` | `--output-format tab`, `vcf`, `json` | `vcf` is the default |
| `--refseq`, `--merged` | `--gff3` | pass the NCBI GFF3, or both, with `LABEL=path`; the `SOURCE` field of the CSQ says which |
| `--check_existing` | `--check_existing`, or `--cache-dir` alone | adds `CLIN_SIG`, `SOMATIC` and `PHENO`, per allele as VEP computes them, and fills `Existing_variation`; needs `--cache-dir`, which is the variation cache and has no other use here, so naming one is `--check_existing`. The frequencies and the PubMed ids need their own flags, as in VEP |
| `--af`, `--af_1kg`, `--af_gnomade`, `--af_gnomadg`, `--max_af`, `--pubmed`, `--clin_sig_allele` | the same | the columns of VEP's `FLAG_FIELDS`, in VEP's order, from the records of the variation cache; every frequency flag, `--max_af` and `--pubmed` imply `--check_existing`, as in VEP; they need `--cache-dir` and stop the run without it. `--af_gnomad`, which VEP keeps as another name for `--af_gnomade`, is accepted |
| `--hgvs` | `--hgvs` | adds `HGVSc`, `HGVSp` and `HGVS_OFFSET` |
| `--hgvsg` | `--hgvsg` | written as VEP writes it, an indel 3'-shifted over the genome and a repeat a `dup`, which reads the reference, so give `--fasta` |
| `--symbol`, `--biotype`, `--canonical`, `--numbers`, `--variant_class`, `--show_ref_allele` | the same | tab columns, in VEP's order; the VCF `CSQ` has most of them always |
| `--total_length` | `--total_length` | `cDNA_position`, `CDS_position` and `Protein_position` as `14/203`, in tab and VCF output; JSON writes numbers and gives no totals |
| `--no_escape` | `--no_escape` | a synonymous `HGVSp` ends in `=` and not `%3D`; the space of a gene symbol in the VCF `CSQ` is `_` and not `%20` |
| `--force_overwrite` | `--force_overwrite` | accepted and does nothing: `fastvep annotate` always writes over its output file |
| `--exclude_predicted` | `--exclude_predicted` | leaves out `XM_` and `XR_` transcripts (VEP's filter is `not stable_id match ^X._`); a variant left with no transcript is an `intergenic_variant` |
| (VEP's Ensembl cache is built without them) | `--exclude-readthrough`, `--exclude-artifact` | leave out readthrough transcripts and those of biotype `artifact`, which an Ensembl GFF3 has and VEP's cache does not; give them to get the cache's transcript set. None of the three exclusion flags is on by default |
| `--everything` | `--everything` | `--hgvs --symbol --numbers --canonical --biotype --variant_class`, and with `--cache-dir` the frequency flags too; VEP's set has more |
| `--pick`, `--pick_allele`, `--pick_allele_gene`, `--flag_pick`, `--flag_pick_allele`, `--flag_pick_allele_gene`, `--pick_order` | the same | |
| `--distance` | `--distance` | |

Any other VEP flag is refused, with `error: unexpected argument`, so a command line that uses one fails at once, and for each flag of the list below the message says what there is instead (`--cache` and `--offline`: the gene models are a GFF3 and the known variants `--cache-dir`; `--refseq`: pass the GFF3; `--tab`: `--output-format tab`; `--plugin` and `--custom`: `sa-build`; `--fork`: `RAYON_NUM_THREADS`).
Of the ones a production command line is likely to carry, that is `--cache`, `--offline` and `--dir_cache`, `--refseq` and `--merged`, `--tab`, `--vcf` and `--json`, `--assembly` and `--species`, `--use_transcript_ref` and `--bam`, `--fields`, `--plugin`, `--custom`, `--fork`, `--no_stats` and `--quiet`, and the annotation flags that have no counterpart here (`--sift`, `--polyphen`, `--mane`, `--tsl`, `--appris`, `--ccds`, `--protein`, `--domains`, `--uniprot`, `--gene_phenotype`, `--regulatory`, `--nearest`, `--shift_3prime`, `--no_intergenic`, `--coding_only`, `--minimal`, `--allele_number`, `--individual`, `--vcf_info_field` and `--terms`).
`--af_esp` and `--af_exac` are refused as well, and VEP 116.3 has neither.
That list was run, each flag on its own, against the binary.

Population frequencies beyond what a variation cache carries come from supplementary databases you build with `sa-build` (gnomAD, 1000 Genomes, TOPMed), and appear as `FV_*` fields and columns.
The columns above are VEP's own, computed from the cache's records.

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

and what carries over to fastVEP, the flags typed as VEP spells them:

```bash
fastvep annotate -i in.vcf -o out.tab --output-format tab \
    --gff3 refseq.gff3 --fasta $ref --cache-dir $dir_cache \
    --show_ref_allele --variant_class --biotype --af_1kg --check_existing \
    --hgvs --hgvsg --exclude_predicted --no_escape --numbers --total_length \
    --force_overwrite
```

What is dropped, and why: `--cache` and `--offline` (fastVEP never connects to a database), `--refseq` (the choice of GFF3), `--assembly` (the assembly is whatever the GFF3 and FASTA are), `--use_transcript_ref` and `--bam` (no BAM-edited transcripts), `--plugin` and `--custom` (see below).
`--dir_cache` becomes `--cache-dir`.

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
fastvep annotate -i in.vcf --transcript-cache refseq.cache --fasta GRCh38.fa --exclude_predicted ...
```

The synonyms file is in the format VEP reads from a cache directory's `chr_synonyms.txt`, one line of equivalent names per contig.
[README.md](../README.md#chromosome-naming-across-sources---synonyms) has the details.
Give `--exclude_predicted` as you would to `vep --refseq`: the GFF3 holds the predicted `XM_` and `XR_` models, and without the flag they are annotated.
If you bgzip and tabix the GFF3 so that a run reads only the regions of its input, sort it with a stable sort (`sort -s -k1,1 -k4,4n`, then `bgzip` and `tabix -p gff`): of equally long transcripts the first in the file is the one marked `CANONICAL`, as in VEP's cache for 48 of the 49 genes of chr21 and chr22 that have a tie, and a plain `sort -k1,1 -k4,4n -k5,5n` puts records of one start and end in the order of their text and changes the pick for 25 of those 49.
A run says so when it sees a file in that order.
A read through tabix marks the same transcripts as a read of the whole file, which it did not before (4,266 of 212,240 rows differed).

fastVEP reads a RefSeq GFF3 as VEP's RefSeq cache holds it.
A transcript is its `transcript_id` (`NM_000001.3`), a gene is the `GeneID` of its Dbxref (`7157`), and the protein in an `HGVSp` is the id NCBI gives (`NP_000001.2`).
`SYMBOL_SOURCE` is `EntrezGene`, `BIOTYPE` is the `gene_biotype` of the gene (with `misc_RNA` for a non-coding transcript of a protein-coding gene), and `CANONICAL` marks the MANE Select transcript of a gene that has one and otherwise the transcript that spans the most of the genome.
A start codon the annotation reads as methionine (`transl_except ... aa:Met`) is `M` in `Amino_acids`.
The mature microRNA records that hang off a precursor get no row, as in VEP's.
Those rules were read off VEP's cache run and not off VEP's code, because the cache is built from a database; [VEP_DIVERGENCE.md, Against VEP's own cache](VEP_DIVERGENCE.md#against-veps-own-cache) has the evidence and the cells that are left.

That comparison is the case this page is about.
Chromosomes 21 and 22 of VEP 116.3's own RefSeq cache were run through VEP (`--cache --offline --refseq --exclude_predicted`) and through fastVEP over the matching GFF3: 203,125 rows of 21,052 changes at the edges of the transcripts, and the HG002 and ClinVar samples on those two chromosomes.
After the changes of this series the cells that differ are in these columns, and all of them have a cause: `HGVSp` (10,631 of the 203,125 rows of the edge set), `Consequence` and `IMPACT` (1,055 and 1,055 rows of it), which are what VEP writes at a start codon, at a stop codon or across an intron where its output contradicts the sequence (divergences 12 to 17 of the divergence document, which says which tool is right and why), `CANONICAL` where the span rule misses (150 rows of the edge set, and 15 of the 15,426 rows of HG002), and `HGVSc` and `HGVS_OFFSET` of an allele of a record with several ALTs (122 and 146 rows of HG002).
On HG002, 2 rows are only fastVEP's and 1 only VEP's.
On the ClinVar sample, 26,383 rows, 2 cells differ, both of them a frameshift that VEP's walk over the genome carries out of an exon, where it writes no protein change and fastVEP writes the one the transcript has.
The `BIOTYPE` rule was read off chr21 and chr22: a biotype that is not on them keeps the name VEP's `--gff` reader gives.

Two things are not covered:

- VEP's RefSeq cache carries transcripts whose sequence differs from the genome, repaired from a BAM (`--bam`, `--use_transcript_ref`; `_rna_edit` in the cache, `BAM_EDIT=OK` in the output): 59 of the edge set's 1,061 and 194 and 57 of the two samples'.
  fastVEP builds transcript sequence from the GFF3 and the FASTA and has no equivalent, so those were left out of the comparison and its agreement with VEP on them is unknown.
- The cache's own columns (`REFSEQ_MATCH`, `REFSEQ_OFFSET`, `GIVEN_REF`, `USED_REF` and `BAM_EDIT`) are not written.

## The output

The tab file and the VCF `CSQ` follow the cells of VEP 116.3, which was run to check them: `Uploaded_variation` (the VCF line itself for a record with no ID), `Location`, `-` for an absent value, `Amino_acids` as a single letter where the residue is unchanged, `HGVSp` with `=` as `%3D` (unless `--no_escape`), and no feature on an intergenic row.
The columns are the ones you ask for, in VEP's order.
[README.md](../README.md#tab-output) lists the rules.
The columns VEP computes per allele and per known variant (`AF` and its populations, `MAX_AF`, `MAX_AF_POPS`, `CLIN_SIG`, `SOMATIC`, `PHENO`, `PUBMED`) are written as VEP writes them when you give the flag that asks for them: 106,673 (variant, allele) pairs of HG002's chr21 and chr22 against VEP's real cache, with `--af_1kg --check_existing` and again with all seven flags, and not one of the 4.7 million cells differs.
The `AF`, `CLIN_SIG`, `SOMATIC`, `PHENO` and `PUBMED` fields of the VCF `CSQ` keep their place among its fields and hold VEP's per-allele value or nothing; fastVEP used to summarise the whole variant in them (its `AF` was the first frequency of any known variant of the site, which is not VEP's `AF`, and for a record with several alleles it put one allele's value on every allele's row) and no longer does.
Two distinct allele-specific `CLIN_SIG` values are joined in a fixed order, where VEP's order is that of a Perl hash and differs from one run to the next.
JSON carries none of these columns and says so for each flag given.

### What will differ

Compared with VEP 116.3 run on the same GFF3 (`--gff`), on a whole-genome callset (20,241 records of GIAB HG002, annotated by both with the same GRCh38 Ensembl 115 GFF3 and FASTA), 3,851 of the 5,023,485 cells of 27 columns of the tab files differ (0.08 %), in four columns (`HGVS_OFFSET`, `HGVSc`, `HGVSg` and, for 17 cells, `HGVSp`), and every other column is the same on every row.
The remaining three, `SYMBOL_SOURCE`, `HGNC_ID` and `FLAGS`, are the ones a GFF3 cannot carry and a `--gff` run writes as `-` (below).
Compared with VEP's own cache, which is what `vep --cache` writes and what this page is about, the rows of HG002 on chr21 and chr22 differ in `HGVSc` and `HGVS_OFFSET` (338 and 437 of 81,309 rows on the Ensembl cache, 122 and 146 of 15,426 on the RefSeq one), in `FLAGS` and `HGNC_ID` on the Ensembl one (35 and 533), and in `CANONICAL` on the RefSeq one (15).
[VEP_DIVERGENCE.md, Part 3](VEP_DIVERGENCE.md#part-3---the-tab-output-against-vep-1163) lists each class with its size and says which tool is right.
The ones you are likely to meet:

- **A record with several ALT alleles.**
  VEP 116 no longer splits it into one variant per allele, and then writes `HGVSc` of an insertion or deletion allele without the 3' shift and without `HGVS_OFFSET`, and `HGVSg` of an allele that trims to a substitution as a `delins` that does not reproduce the variant.
  fastVEP writes every allele as it would be alone.
  On HG002, 36 % of the rows of such records differ in `HGVSc`, and nothing else does because of it.
  A `*` ALT counts as another allele: `TC>T,*` is shifted by fastVEP and not by VEP.
- **`HGVSp` at a stop codon.**
  VEP writes `p.Ter583=` beside `stop_lost` on 9,166 rows of its own Ensembl cache run, because it reads an unresolved residue as a terminator; fastVEP writes what the edited protein says (`p.Ter583TrpextTer29`).
- **A change at the start codon.**
  fastVEP does not call the initiator lost where the protein still begins with its Met, or where the whole coding sequence is behind the edit (a deletion from the 5' UTR into the initiator inside a repeat, `c.1dup`): VEP writes `start_lost` beside `start_retained_variant`, HIGH, two terms that contradict each other, and where the change does not reach the UTR a `frameshift_variant` and `p.Met1?` as well, for a protein that has no change.
  Where the initiator is lost fastVEP writes `p.Met1?`, as HGVS has it, where VEP writes a range, several residue names, a later residue, or nothing.
  Of the 2,026,951 rows of the Ensembl cache's edge set (45,916 changes at the edges of the transcripts of chr21 and chr22), 9,385 differ in `Consequence` and 54,864 in `HGVSp`, and none of HG002's rows does; on the ClinVar sample 66 rows differ in `Consequence` from a `--gff` run (54 of them on a transcript the cache flags incomplete) and 23 cells in `HGVSp` from VEP's cache.
- **A change the 3'-rule carries over an intron.**
  VEP moves it along the genome and describes it where the walk ends, which is another change of the transcript or none, so it writes no protein change for 582 coding rows of the ClinVar sample and 17 of HG002's; fastVEP writes the protein change the change is in the transcript, and an independent check that applies each description back to the sequence finds every one right.
- **`SYMBOL_SOURCE`, `HGNC_ID`, `FLAGS` and the terms that depend on them.**
  fastVEP follows VEP's cache, which has them (`HGNC` with its accession, `cds_start_NF` and `cds_end_NF`), and a `vep --gff` run has `-` in all three, so a comparison against `--gff` output differs in them on most rows and should leave them out.
  A transcript the cache flags incomplete has no start or stop codon to lose, which a `--gff` run cannot know: on the ClinVar 2-star+ sample (3,136,295 rows) 66 `Consequence` rows differ from a `--gff` run, 54 of them on such a transcript (`start_lost` in VEP, `missense_variant` here), where against VEP's Ensembl cache the same sample differs in 1 `Consequence` cell, a mature microRNA.
  `HGNC_ID` also differs where VEP copies one gene's ID to every other gene of the same symbol that it loaded.
- **A transcript the cache does not hold.**
  An Ensembl GFF3 has readthrough transcripts and transcripts of biotype `artifact`, and VEP's cache does not: give `--exclude-readthrough --exclude-artifact` for the cache's transcript set.
  A RefSeq GFF3 has the predicted ones: give `--exclude_predicted`.
  Without them fastVEP writes rows VEP does not.
- **An unlisted GFF3 feature type.**
  VEP's `--gff` reader keeps only some feature types (`unconfirmed_transcript` and `scRNA` are not among them), so a variant near one is `intergenic_variant` in VEP's `--gff` run and an upstream or downstream variant in fastVEP: 121 rows only fastVEP writes on HG002, and 7 only VEP writes, besides the 186,055 rows both do.
- **A replacement across the first or last base of a transcript.**
  VEP writes an `HGVSc` range the wrong way round (`c.*349_*348insGT`); 11 cells of the 36,958 rows of the edge set.

The rest of the list that used to be here is gone.
A change that runs off the end of a transcript or an exon, an indel written with bases the alleles share at the back, and a `*` allele are written as VEP writes them, and so are the position, `EXON` and `INTRON` cells at the edge of an exon, and the consequence terms at a start or stop codon, which were fixed in this series and are measured in the divergence document.

## Before you switch

Run both tools on a sample of your own sites, with the flags you use, and compare row by row.
`validation/run_divergence.sh` and `validation/compare_rows.py` do that for the Ensembl GFF3 and write one file per disagreeing field, which is the quickest way to see whether a disagreement touches the fields you read.
`validation/compare_tab.py` compares two tab files cell by cell, `validation/check_hgvs_roundtrip.py` applies each `HGVSc`, `HGVSp` and `HGVSg` of a tab file back to the sequence and says which are not the change the VCF made (it needs no VEP run, and it is how the defects of VEP's own output that this series does not copy were found), `validation/run_cache_oracle.sh` runs VEP on its own cache beside fastVEP over the matching GFF3 and compares the two (`validation/compare_cache_run.py`), and `validation/run_known_columns.sh` does the same for the known-variant columns.
