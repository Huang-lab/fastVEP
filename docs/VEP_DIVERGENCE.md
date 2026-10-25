# Where fastVEP differs from Ensembl VEP

fastVEP's consequence and HGVS output is a port of Ensembl VEP's own model, and the goal is to
agree with it.
Agreement is the *evidence* that the port is faithful, not the objective: the output is a
prediction a clinician may act on, so where VEP is demonstrably wrong in a way that changes a
call, fastVEP is right instead and says so here.

This document is the complete list, in both directions, and every row in it is measured rather
than argued.

## What is measured

Against **real Ensembl VEP 115.1** (Docker `ensemblorg/ensembl-vep:release_115.1`, `--gff` mode,
the same GRCh38 Ensembl 115 GFF3 and FASTA fastVEP reads, `--hgvs --symbol --canonical
--allele_number`), on three inputs that exercise different code paths:

| Input | Variants | Matched (variant, allele, transcript) rows | Coding rows |
|---|---:|---:|---:|
| ClinVar 2-star+, stratified towards the hard shapes | 6,600 | 151,684 | 74,150 (48.9 %) |
| GIAB HG002 WGS, 1-in-200 systematic | 20,241 | 118,956 | 861 (0.7 %) |
| ClinVar in-frame deletions | 400 | 9,141 | 5,212 (57.0 %) |

The ClinVar sample is built to be hard: 54.5 % of its variants are not SNVs, against 7.2 % of the
673,660-variant ClinVar 2-star+ set it is drawn from, so its rates are the rates on the shapes
that disagree.
The HG002 sample is the other end of that range, and is what an ordinary WGS callset looks like.

Reproduce all three with

```bash
bash validation/run_divergence.sh            # samples, both tools, per-field counts
```

which cuts the samples with `validation/sample_variants.py` (deterministic - every k-th record of
a class, no seed to carry around), runs both tools, and compares row by row with
`validation/compare_rows.py`, which also writes one TSV per field holding every disagreeing row.
The counts below are those files.

Two things the comparison counts but does not rate.
**Rows only fastVEP has**: 9,579 on the ClinVar sample and 67,226 genome-wide, almost all of them
transcripts of `ncRNA_gene` records, which VEP's `--gff` reader discards with a warning and
fastVEP reads (#98).
**Rows only VEP has**: none on the ClinVar sample, and 9,870 genome-wide, every one of them an
`intergenic_variant` on a variant whose only neighbours are those same non-coding genes.
Neither is a disagreement about a variant both tools annotated.

## Current agreement

| Field | Scope | ClinVar sample | Genome-wide (HG002) |
|---|---|---:|---:|
| `Amino_acids` | coding rows | 0 rows, **100 %** | 0 rows, **100 %** |
| `Codons` | coding rows | 0 rows, **100 %** | 0 rows, **100 %** |
| Splice terms | all rows | 1 row, 99.999 % | 0 rows, **100 %** |
| Consequence terms | coding rows | 59 rows, 99.920 % | 0 rows, **100 %** |
| Whole consequence set | all rows | 60 rows, 99.960 % | 0 rows, **100 %** |
| `IMPACT` | all rows | 47 rows, 99.969 % | 0 rows, **100 %** |
| `HGVSc` | all rows | 41 rows, 99.973 % | 23 rows, 99.981 % |
| `HGVSp` | all rows | 778 rows, 99.487 % | 17 rows, 99.986 % |

Every field that carries a clinical call agrees completely on a genome-wide callset.
`HGVSc` and `HGVSp` are the exceptions, and both are down to the two gaps in Part 2: a description
one tool writes and the other declines to.
No row on any of the three inputs has both tools naming a `c.` change and naming it differently.

The 400-variant in-frame deletion set is the only input that exercises protein-level
3'-normalisation at a protein terminus: `IMPACT`, `Amino_acids`, `Codons` and the splice terms
agree on all of it, 10 consequence rows and 9 `HGVSc` rows are the gaps listed in Part 2, and all
43 of its `HGVSp` rows are divergence 9 below, on a single variant.
Those 43 are 99.17 % agreement counted over the 5,192 rows where both tools name a protein change,
which is the figure the manuscript reports, and 99.53 % counted over all 9,141 matched rows.

---

## Part 1 - deliberate divergences

Each is a case where Ensembl's output is either not valid HGVS, or contradicts what the sequence
says, and where matching it would degrade a clinical call.

### 1. A change that introduces a premature stop is HIGH, not a moderate in-frame insertion

**36 consequence rows and all 36 of them `IMPACT`** - the largest divergence, and the one that
moves an ACMG call.

`stop_retained` in `Utils/VariationEffect.pm` (release/115, l. 1284) defers to
`ref_eq_alt_sequence` in the same file (l. 1321), whose first clause is:

```perl
return 1 if ( ($ref_pep eq substr($alt_pep, 0, 1) && $alt_pep =~ /\*/) || ... );
```

That asks whether the replacement keeps the residue it starts on and introduces a terminator
*anywhere* - not whether the annotated terminator survived.
Both `frameshift` (l. 1435) and `stop_gained` (l. 1208) return 0 when `stop_retained` holds, so
the variant comes out `inframe_insertion,stop_retained_variant`, MODERATE.

| Variant | Transcript | ClinVar | Ensembl VEP 115.1 | fastVEP |
|---|---|---|---|---|
| BRCA1 `c.5030_5033dup` | ENST00000357654 | Pathogenic | `inframe_insertion,stop_retained_variant`, MODERATE | `stop_gained,frameshift_variant`, HIGH |
| BRCA1 `c.1499_1508dup` | ENST00000357654 | Pathogenic | same, MODERATE | HIGH |
| BRCA2 `c.3205_3206insAATTGCAGTCAATTAATAT` | ENST00000544455 | Pathogenic | same, MODERATE | HIGH |
| ITGB3 `c.122_125dup` | ENST00000559488 | Pathogenic | same, MODERATE | HIGH |
| LDLR `c.1309_1310insTCGCTCTGGACACGTAGGTGG` | ENST00000557933 | Pathogenic | same, MODERATE | `stop_gained,inframe_insertion`, HIGH |

VEP's own `Amino_acids` for the first of these is `N/N*X`: the reference codon translates to Asn,
the edited window to Asn then a terminator then an incomplete codon.
Nothing about that says the annotated stop at residue 1863 survived.

Across the whole ClinVar 2-star+ set the clause fires on **52 variants** where fastVEP reports
`stop_gained` - 35 of them frameshifts and 17 in-frame insertions - in ABCC9, BRCA1, BRCA2, BRIP1,
CLCN7, COL7A1, DNAH11, DNAI1, DSP, EXT1, FLCN, GATA2, ITGB3, KCNA2, LDLR, LZTR1, MANBA, MLH1,
MSH6, PCDH15, PMS2, PTS, RPGR, RUNX1, RYR2, SACS, SAMD9L, SDHB, SPG11, SPTAN1, STK11, TGFB2,
TMEM127, TP53, TREX1, TRNT1, TTN and USH2A.
**38 of the 52 are Pathogenic or Likely pathogenic in ClinVar.**
All 52 would lose PVS1 under Ensembl's reading, because neither `inframe_insertion` nor
`stop_retained_variant` is a null variant.

fastVEP reproduces Ensembl's other two `ref_eq_alt_sequence` clauses, which do test the
terminator: one asks whether it sits at the same residue on both sides, the other whether the
edited protein still matches the reference over the reference's own length and grows by fewer than
three residues past it.
Only the first clause is refused, in `terms_for_window`
(`crates/fastvep-consequence/src/predictor.rs`).

### 2. A start codon the variant deletes is lost, whichever spelling HGVS prefers

**10 consequence rows, all 10 of them `IMPACT`: HIGH here, LOW in VEP.**

`_ins_del_start_altered` (`VariationEffect.pm`, l. 1028) edits the transcript and asks whether the
initiator survived - but it reads `$bvfo->cdna_start`, the *shifted* coordinate, while
`_overlaps_start_codon` (l. 965) reads `cdna_start_unshifted`.
So a deletion the 3'-rule can rewrite as removing something else is tested at a position it does
not occupy, and the initiator "survives".

MLH1 `3:37014476 AAATG>A` on ENST00000536378 deletes the four bases `AATG`, which is the initiator
`ATG` at 37014478-80 plus the base before it.
The block repeats immediately (`…A AATG AATG G…`), so HGVS writes the deletion one copy along, as
`c.4_7del`, and both tools agree on that string.
VEP then calls the row `start_retained_variant`, LOW - while writing `p.Asn2ValfsTer10` for its
own HGVSp, which is a frameshift by any reading.
fastVEP calls it `start_lost`, HIGH.

Nine of the ten rows are that variant across MLH1's transcripts.

### 3. A synonymous multi-residue window names the whole span

**60 `HGVSp` rows.**

A change spanning two codons that leaves both residues unchanged is `p.Leu2672_Ile2673=`
(FLNC `7:128858463 CATT>TATC`, ENST00000950263).
Ensembl writes `p.LeuIle2672=` - two three-letter codes sharing one position - which is not a form
HGVS defines.
fastVEP writes the span.

### 4. A `*` coordinate keeps the exonic anchor it is counted from

**23 `HGVSc` rows on the ClinVar sample, 12 genome-wide.**

For an intronic position in the part of a transcript that follows the terminator, Ensembl writes
the offset as if it were the whole coordinate, dropping the exonic anchor it is counted from:

| Variant | Transcript | Ensembl VEP 115.1 | fastVEP |
|---|---|---|---|
| PMS2 `7:5992071 AAGG>A` | ENST00000699951 | `c.*-17_*-15del` | `c.804-17_804-15del` |
| SACS `13:23332502 G>GC` | ENST00000683680 | `c.*-2551dup` | `c.2319-2551dup` |
| COMMD7 `20:32738081 C>T` | ENST00000610160 | `c.*4395G>A` | `c.174+4395G>A` |
| CCDC88A `2:55298288 GA>G` | ENST00000644415 | `c.*1539del` | `c.898+1539del` |

HGVS `*N` positions count forward from the terminator in the transcript, so `*-944` is
self-contradictory and `*4395` names an mRNA base 4,395 residues past a stop the transcript does
not have; in both forms the `*` carries no exonic position at all.
fastVEP anchors the position to the exonic base it is offset from.
COMMD7 is the worked example: the transcript's last CDS base is `c.174` at genomic 32742476, the
variant sits at 32738081 in the intron that follows it, and 32742476 - 32738081 = 4,395.

The anchor numbers are not independently verified beyond that one; what is verified is that
Ensembl's form is malformed.

### 5. A frameshift's terminator distance counts to the stop the frameshift creates

**411 `HGVSp` rows** - the largest `HGVSp` divergence of all.

`p.Leu151ProfsTer39` says the new reading frame runs 39 residues from the first changed one before
it hits a stop.
That is a fact about the edited protein, so it can be computed rather than argued about.

Six rows were checked by rebuilding the transcript from the Ensembl 115 GFF3 and FASTA, applying
the variant, and translating - independently of either tool, and continuing past the annotated
terminator because that is where a frameshift's stop usually lands.
The reconstruction is corroborated by the two tools themselves: on every row it names the same
reference residue at the same position as both of them, so all three agree on where the frame
breaks and only the distance to the new stop is in dispute.

| Gene | Transcript | Variant | Computed | fastVEP | Ensembl VEP 115.1 |
|---|---|---|---:|---:|---:|
| MPV17 | ENST00000233545 | `c.451dup` | **39** | 39 | 41 |
| ARSB | ENST00000264914 | `5:78780421 AG>A` | **48** | 48 | 175 |
| VWF | ENST00000321023 | `12:6123146 C>CA` | **302** | 302 | 56 |
| CFI | ENST00000882820 | `4:109740925 C>CT` | **14** | 14 | 22 |
| PCDH15 | ENST00000373956 | `10:54079377 TCAGT…>T` | **1275** | 1275 | 23 |
| NTHL1 | ENST00000651522 | `16:2040005 GC>G` | **55** | 55 | `Ter?` |

fastVEP matched the computed answer on all six, over differences of 2 to 1,252 and in both
directions - Ensembl is shorter on some rows and longer on others, and it is wrong either way.
The NTHL1 row is counted under divergence 6 rather than here, because VEP declines to give a
distance at all; it is in this table because it was checked the same way.

In every one the new stop lies past the reference protein's own terminator, which is the whole
difficulty: the frameshift runs off the end of the annotated CDS and the stop that ends it is in
what was the 3' UTR.
fastVEP translates a CDS that continues downstream for exactly that reason.
Ensembl's `_get_alternate_cds` builds its edited CDS from 3'-shifted coordinates while taking
consequences from unshifted ones, and the two do not line up.

### 6. A frame with no stop in the transcript is reported as having none

**51 `HGVSp` rows**, all of them `Ter?` on one side and a distance on the other.

The NTHL1 row in the table above is one of them, checked the same way: the stop really is 55
residues along, and `Ter?` is not an answer about it.
Where fastVEP writes `Ter?` instead - ITGA2B `p.Ter1040Trpfs` on ENST00000262407, CLN3
`p.Gly37Valfs` on ENST00000561505, RUNX1 `p.Arg380Profs` on ENST00000399240 - the shifted frame
runs to the end of the transcript without another stop, and there is nothing to count.

### 7. A lost terminator's extension is counted, not written off

**20 `HGVSp` rows on the ClinVar sample, 1 genome-wide.**

`p.Ter486GluextTer36` says the terminator at 486 became Glu and the protein runs 36 residues
further before the next stop.
Both halves are in the sequence, so both are computed:

| Gene | Transcript | Variant | Computed | fastVEP | Ensembl VEP 115.1 |
|---|---|---|---:|---:|---:|
| TP53 | ENST00000620739 | `17:7669611 A>T` | **9** | `p.Ter355ArgextTer9` | `p.Ter355ArgextTer1` |
| MUTYH | ENST00000674679 | `1:45333589 TC>AA` | **33** | `p.Ter48PheextTer33` | `p.Ter48PheextTer19` |
| DNAH10 | ENST00000538983 | `12:123932637 A>G` | **3** | `p.Ter82TrpextTer3` | `p.Ter82TrpextTer3` |

VEP's distances here fail for the same reason as divergence 5, and by the same shifted-coordinate
mechanism: on the TP53 variant it returns 1, 2, 5, 11, 30 and 106 on six different transcripts of
the same gene, all of which reconstruct to 9.
The DNAH10 row is included because it is one VEP gets right, and fastVEP now agrees with it -
fastVEP used to write `p.Trp82ext*?` for every one of these, naming neither the residue that was
lost nor the protein the loss adds.

### 8. A frameshift's first changed residue is read off the edited protein

**3 `HGVSp` rows**, and the two variants behind them both reconstruct in fastVEP's favour:

| Gene | Transcript | Variant | Computed | fastVEP | Ensembl VEP 115.1 |
|---|---|---|---|---|---|
| PCDH15 | ENST00000373956 | `10:54079377 TCAGT…>T` | Gly650→Asp, stop 1,275 on | `p.Gly650AspfsTer1275` | `p.Gly650GlufsTer23` |
| FGFR2 | ENST00000429361 | `10:121479877 TTA>T` | Ter372→Asn, stop 7 on | `p.Ter372AsnfsTer7` | `p.Ter372ThrfsTer?` |

VEP names a residue the edited frame does not produce at that position.
**This entry was listed as a fastVEP defect before it was computed**; computing it reverses the
verdict, exactly as it did for divergence 5.

### 9. The smaller ones

| What | Rows | fastVEP | Ensembl VEP 115.1 | Why ours |
|---|---:|---|---|---|
| A two-residue deletion in a poly-Glu C-terminus. NT5C2 `10:103089678 TTCCTCC>T`, ENSP00000502205, protein 561 aa ending `EEEEEEE` | 43 of the in-frame deletion set, on one variant | `p.Glu560_Glu561del` | `p.Glu559_Glu560del` | The 3'-most placement is the one HGVS asks for. VEP's `_shift_3prime` halts *n*-1 residues before the terminus for an *n*-residue change ([#94](https://github.com/Huang-lab/fastVEP/issues/94)) |
| A transcript whose CDS is annotated 5' incomplete (first CDS record carries phase 1 or 2, so no initiator is annotated). TCIRG1 ENST00000698256, NF1 ENST00000696141, NFIA ENST00000496712 | 3 | no `start_lost` | adds `start_lost` | There is no `ATG` to lose. Declining is the lower-impact answer, which is worth saying plainly |
| `start_lost` and `start_retained_variant` on the same row, where a *substitution* replaces a non-canonical initiator with `ATG`. NBEA `13:35171273 C>T` (`aCg`→`aTg`), ENST00000629018 | 2 | `start_lost` | `start_lost,start_retained_variant` | The pair contradicts itself. `IMPACT` is HIGH either way. Ensembl reaches `start_retained_variant` for a substitution through `_snp_start_altered`, which this window does not model; for a length change both tools can report the pair, and no variant of the ClinVar 2-star+ set does |
| The polypyrimidine tract on a transcript carrying a frameshift intron. PTEN `10:87933000 C>CT`, ENST00000693560 | 1 | `splice_polypyrimidine_tract_variant` at -3 to -17 | the term stops at -13 | An intron of 12 bp or less makes VEP treat every exon as 12 bases wider (`_overlapped_exons`, `BaseTranscriptVariation.pm` l. 861). At the same 15 positions the other 14 PTEN transcripts get the full window from VEP, and this is the only one with a 1 bp intron |

### 10. `--pick` keeps every allele of the transcript it picks

VEP's `--pick` is documented as "pick one line or block of consequence data per variant", and it
means it literally.
On TP53 `17:7676154 G>A,C`, where both tools report 76 consequence entries unpicked, VEP's
`--pick` leaves **one** and fastVEP's leaves **two**.
The entry VEP drops is not a duplicate: it is the only annotation of the `C` allele anywhere in
the record.

| Option | fastVEP entries | VEP entries | Agree |
|---|---:|---:|:-:|
| `--pick` | 2 | 1 | no |
| `--pick-allele` | 2 | 2 | yes |
| `--pick-allele-gene` | 2 | 2 | yes |
| `--flag-pick` | 76, 2 flagged | 76, 1 flagged | no |
| `--flag-pick-allele` | 76, 2 flagged | 76, 2 flagged | yes |
| `--flag-pick-allele-gene` | 76, 2 flagged | 76, 2 flagged | yes |

A missing alt allele is indistinguishable from an alt that had no consequence, in a file whose
reader is deciding what to report on.
So `--pick` reduces to one transcript and keeps that transcript's annotation for every allele,
and the option for VEP's behaviour is `--pick-allele`, which agrees with it exactly.
`--flag-pick` flags the same set `--pick` would have kept, for the same reason.

The four allele-scoped options agree with VEP on which transcript wins as well as how many
entries survive: over five two-alt sites at overlapping-gene loci on chromosome 22 (HDHD5,
HIRA/C22orf39, GNB1L/RTL10, PI4KA/SERPIND1, CLTCL1/SLC25A1), `--pick-allele-gene` and
`--flag-pick-allele-gene` picked the same (allele, gene, transcript) set as VEP at 5 of 5, once
the comparison is restricted to the genes VEP reports at all - see the `ncRNA_gene` note under
[What is measured](#what-is-measured).

### 11. An allele of a record with several ALTs is described as it would be alone

**HG002: 35 `HGVSg` cells, 1,717 `HGVSc` cells, 2,070 `HGVS_OFFSET` cells.
The world: 709, 108 and 78.**

Ensembl 116 no longer splits a multi-allelic record into one variant per ALT, as 115 did (`InputBuffer.pm` has `$self->{minimal}=1 if $vf->{minimised}` commented out), so what it writes for an allele now depends on what else the record holds.
Splitting every such record into one record per ALT (`validation/split_alts.py`) shows how much.
The same allele, in the same transcript, differs between its record and a record of its own in **1,071 of 1,879 rows** of VEP's output, and in **280** of fastVEP's.
For `HGVSg`, `HGVSc` and `HGVS_OFFSET` it is 709, 107 and 80 cells for VEP, and 8, 0 and 0 for fastVEP.
(What is left of fastVEP's 280 is mostly the span of the untrimmed record - 140 `cDNA_position`, 131 `DISTANCE`, 29 `CDS_position` and 26 `Codons` cells, and 17 `Consequence`, `IMPACT`, `Amino_acids`, `HGVSp` and `EXON` cells that follow from it - which is what Part 2 says, under *Multi-allelic records, and what is still untrimmed on them*, and Ensembl 116 has it as well, in 158 `DISTANCE` and 154 `cDNA_position` cells.
The 8 `HGVSg` cells are an inversion inside shared bases, `inv` in the record and `delins` alone, as Ensembl writes it alone.)

With every ALT its own record, fastVEP and VEP 116.3 write the same `HGVSg` on all 7,213 rows of the world, and the same `HGVS_OFFSET` on every row on which both write an `HGVSc`.

**`HGVSg`.**
`hgvs_genomic` trims an allele of a multi-allelic record against the reference and takes its class from the trimmed pair, but a change that is then a substitution or a replacement is described by the *untrimmed* allele (`$check_allele` is copied before the trim, and only an insertion or a deletion has it replaced).
`5:100381015 CATAA>AATAA,C` is a substitution of one base, and Ensembl 116.3 writes `5:g.100381015delinsAATAA`: the span of the substitution with the five bases of the allele.
Applied to the genome, that inserts four bases the variant does not have.
fastVEP writes `5:g.100381015C>A`, which is also the change Ensembl's own `HGVSc` of that allele describes.

That is checked rather than argued, by applying each description to the genome (`validation/check_hgvsg_roundtrip.py`).
In the world, the 709 cells are **336 (record, allele) pairs, and for all 336 Ensembl's description does not reproduce the allele and fastVEP's does**.
On HG002 three records are affected, `5:100381015 CATAA>AATAA,C`, `8:86711835 CTATA>ATATA,C` and `10:90497339 AT>TT,A`, and 20,238 of 20,241 records have the same set of descriptions in both tools.
Nothing else differs, which is what [`hgvsg_vep`](../crates/fastvep-annotate/src/hgvs_normalize.rs) being a transcription of `hgvs_genomic` is for: it agrees with VEP 116.3 on all 1,472 alleles of `crates/fastvep-annotate/fixtures/vep_hgvsg` and on every row of the world that is not this one.
Both generators keep the records this affects out of the committed fixtures, and `vep_tab_world/gen_world.py --all` puts them back.

**`HGVSc` and `HGVS_OFFSET`.**
`_return_3prime` (`TranscriptVariationAllele.pm`) returns before it shifts anything `unless ($var_class eq 'insertion' || $var_class eq 'deletion')`, where `$var_class` is `$vf->var_class`: the class of the whole record, not of the allele.
A record whose ALTs include an indel and anything else, a second ALT or a `*`, is `indel` or `sequence_alteration`, and the shift is never tried.
In the world, none of the 1,839 rows of multi-allelic records of any class but `insertion` has an `HGVS_OFFSET` in VEP's output, and 5 of the 37 rows of class `insertion` have one.
So the same deletion is written `c.92-111634del` in its record and `c.92-111619del` alone, fifteen bases from where HGVS places it, in `1:14068799 CTT>CT,C` on ENST00000636203, with `HGVS_OFFSET` empty and 15; and an insertion that repeats its neighbours is `c.600-2233_600-2232insCACACACA` in `1:227162868 A>ATGTGTGTG,G` on ENST00000334218 and `c.600-2240_600-2233dup` alone.
Both pairs are Ensembl 116.3 run on those records, as they are in HG002 and as split by `split_alts.py`.
HGVS places a change at its most 3' position and writes a repeat as a duplication whatever else the VCF line holds, so fastVEP gives every allele the description it has alone, and counts `HGVS_OFFSET` from where the VCF anchors the allele (the base all its alleles share is the anchor), which is also the same alone and in the record.

For a user this is the part of #139 that does not go away on migrating.
On HG002, 1,717 of the 4,780 rows of records with several ALTs (36 %) differ in `HGVSc` because of it, which is 0.92 % of all 186,055 rows.
A comparison that has to match Ensembl 116 string for string on those rows needs Ensembl's unshifted form, which fastVEP does not write; one that wants the same string whatever the record looks like needs fastVEP's.

---

## Part 2 - gaps, where Ensembl is right

These are fastVEP's remaining defects, not disagreements.

### `--pick`'s last tie-break

| Gap | Rows | What it looks like |
|---|---:|---|
| The ninth tie-break, transcript length, is not available | 1 of 5 two-alt overlapping-gene sites on chromosome 22 | PI4KA/SERPIND1 at `22:20780916`. Both transcripts are MANE Select and canonical for their own gene, both TSL 1, both protein_coding, both carry a CCDS, and both are `intron_variant` - so all eight tiers of `--pick-order` tie. VEP then breaks the tie on transcript length and reports PI4KA; fastVEP breaks it on transcript ID and reports SERPIND1 |

`TranscriptVariation` does not carry a transcript length, which is why `--pick-order length` is
rejected with a message saying so rather than accepted and ignored.
The tie is only reachable where a pick compares transcripts of *different* genes, so it affects
`--pick`, `--pick-allele` and their flagging forms, and never `--pick-allele-gene`, which compares
within one gene.

### Consequence terms

| Gap | Rows | What it looks like |
|---|---:|---|
| `3_prime_UTR_variant` / `5_prime_UTR_variant` missing on a delins that spans the CDS boundary | 8 on the ClinVar sample, 10 on the in-frame deletion set | fastVEP `coding_sequence_variant`, VEP `3_prime_UTR_variant,coding_sequence_variant` (CHEK2 `22:28711907 ACT>A`, ENST00000439200) |

That is the whole list: the other 52 consequence rows on the ClinVar sample are the divergences
above.

### HGVSc

| Gap | Rows | What it looks like |
|---|---:|---|
| fastVEP names one where VEP names none | 18 on the ClinVar sample, 11 genome-wide, 9 on the in-frame deletion set | fastVEP `c.-66_-60del` for HNF4A `20:44401297 GGGAGGGC>G` on ENST00000415691, a deletion that straddles the transcript's first base; VEP writes nothing |

That is the whole list for `HGVSc`.
The 23 rows left genome-wide are 11 of those and 12 of divergence 4; the 41 on the ClinVar sample
are 18 of those and 23 where VEP writes a malformed `*` coordinate.
There is no row on any of the three inputs where both tools write a description and the two
disagree.

#### Multi-allelic records, and what is still untrimmed on them

A VCF record carries one position for the whole site, so a multi-allelic indel can only have the
one base *every* allele shares stripped from it.
An allele that is a clean deletion therefore reaches the annotator as a replacement.
TTC28 `22:28225368 AAAGAAG>AAAG,A` on ENST00000612946 is the shape, and the `HGVSc` column is what
this used to get wrong:

| Allele | Ensembl VEP 115.1 | fastVEP |
|---|---|---|
| `AAAG` (a 3 bp deletion) | `-`, `c.553-61772_553-61770del` | `AAG`, `c.553-61772_553-61770del` |
| `A` (a 6 bp deletion) | `-`, `c.553-61775_553-61770del` | `-`, `c.553-61775_553-61770del` |

Ensembl closes this twice over, and neither is gated on `--minimal`.
`Parser.pm`'s `post_process_vfs` sends any record whose alleles differ in length through
`minimise_alleles`, and `InputBuffer::split_variants` then makes it one VariationFeature per ALT,
each trimmed against the reference on its own and annotated separately before being rejoined for
output; `_clip_alleles` in `TranscriptVariationAllele.pm` clips again while building the notation.
fastVEP does the second of those, where the description is built, which is why the descriptions
now agree: **801 of the 824 genome-wide `HGVSc` rows**, and the 11 records where the two tools'
`HGVSg` sets differed are down to 10.

What is still untrimmed is the **reported allele** - `AAG` above where VEP prints `-` - and the
positional fields beside it, which describe the site's span.
The reach of that is wider than it sounds: a comparison keyed on the allele string does not line
those rows up at all, which is why `validation/compare_rows.py` pairs them by `ALLELE_NUM`.
Closing it means splitting the site into one variant per allele the way Ensembl does, which would
put at risk the thing this document's strongest result rests on - the consequence caller reads its
differing region out of the untrimmed pair and agrees with VEP on every row of all three inputs.
That is a change to make on its own evidence, not as a side effect of a nomenclature fix, and it
is not made here.

That is Ensembl 115.1.
Release 116 no longer splits a multi-allelic record into one variant per ALT, which changes what it
writes for them: [Part 3](#part-3---the-tab-output-against-vep-1163) has the measurement.

### HGVSp

| Gap | Rows | What it looks like |
|---|---:|---|
| fastVEP names one where VEP names none | 196 on the ClinVar sample, 17 genome-wide | 177 of the 196 are a frameshift that also earns `splice_region_variant`, where VEP declines to name the protein change at all: MSH6 `2:47806357 TG>T` on ENST00000936511, fastVEP `p.Ala1277HisfsTer4` |
| VEP names one where fastVEP does not | 27 | TRAF3IP1 `2:238332821 TAGTC>T` on ENST00000373327, VEP `p.Ser306GlnfsTer34`; 15 of the 27 are a delins spanning a splice acceptor |
| A frameshift in the first codons, which VEP writes as unknown | 5 | VEP `p.Ala2_?1` (NF1 ENST00000696141) and `p.Leu2?` (TCIRG1 ENST00000698256), fastVEP `p.Ala2GlnfsTer7` and `p.Leu2CysfsTer25`. Three of the five are the 5' incomplete transcripts of divergence 9; the other two (DYNC2H1 ENST00000528670, SGCB ENST00000514133) have a complete initiator and a frameshift that starts two residues into it |
| `Amino_acids` of `X` from an ambiguous reference codon | 3 | VEP `p.Ter157=`, fastVEP `p.Xaa157=`. Which is right depends on whether the position is the terminator or genuinely unknown; unresolved |
| A delins naming an unknown residue the window reaches | 2 | IL7R `5:35860925 GC>TT` on ENST00000515665, fastVEP `p.Gln52_Xaa53delinsHisXaa` where VEP names the resolved change alone, `p.Gln52His`; ATM `11:108310201 A>AT` on ENST00000529588 is the other |

The shape of that list is worth stating.
Of the 778 `HGVSp` rows on the ClinVar sample, **545 are places fastVEP diverges on purpose and
has been checked against the sequence**, and 233 are gaps - 196 of which are fastVEP naming a
change where Ensembl declines to.

---

## Part 3 - the tab output, against VEP 116.3

Everything above was measured against Ensembl VEP 115.1.
The tab output's optional columns, `HGVSg` and `HGVS_OFFSET` had no measurement at all, and the release most people run now is 116, so they have their own, against **Ensembl VEP 116.3** (Docker `ensemblorg/ensembl-vep:release_116.3`).
It reads the same GRCh38 Ensembl 115 GFF3 and FASTA as above and is run with `--tab --symbol --biotype --canonical --numbers --variant_class --show_ref_allele --hgvs --hgvsg`.
`validation/compare_tab.py` compares its file with fastVEP's cell by cell, pairing rows on (`Uploaded_variation`, `Allele`, `Feature`).

| Input | Records | Paired rows | Cells compared | Cells that differ |
|---|---:|---:|---:|---:|
| GIAB HG002 WGS, 1-in-200 systematic (the sample above) | 20,241 | 186,055 | 5,581,650 | 3,876 (0.069 %) |
| A generated world: both strands, coding and non-coding genes, 372 of its 3,000 records with several ALTs | 3,000 | 7,213 | 216,390 | 1,162 (0.54 %) |

The world is built to be hard and not to be representative: it has far more indels in repeats, at transcript ends and at exon edges than any genome, and more multi-allelic sites.
Its rates are the rates on those shapes, and HG002 is what an ordinary WGS callset looks like.
Reproduce the world with `bash validation/run_tab_world.sh` (a record count, a seed and an output directory are its optional arguments; the defaults are 3,000 and 139), which prints every number below that is about it, and HG002 with

```bash
vep --input_file hg002_sample.vcf --gff Homo_sapiens.GRCh38.115.sorted.gff3.gz \
    --fasta Homo_sapiens.GRCh38.dna.primary_assembly.fa --tab --force_overwrite --no_stats \
    --symbol --biotype --canonical --numbers --variant_class --show_ref_allele --hgvs --hgvsg \
    --output_file vep.tab
fastvep annotate -i hg002_sample.vcf --gff3 Homo_sapiens.GRCh38.115.gff3 \
    --fasta Homo_sapiens.GRCh38.dna.primary_assembly.fa --output-format tab -o fv.tab \
    --symbol --biotype --canonical --numbers --variant-class --show-ref-allele --hgvs --hgvsg
python3 validation/compare_tab.py vep.tab fv.tab \
    --ignore SOURCE,Homo_sapiens.GRCh38.115.sorted.gff3.gz
```

On HG002 the 3,876 cells are in five columns: `HGVS_OFFSET` 2,073, `HGVSc` 1,742, `HGVSp` 17, `HGVSg` 35 and `cDNA_position` 9.
`Consequence`, `IMPACT`, `Amino_acids`, `Codons`, `CDS_position`, `Protein_position`, `Location`, `REF_ALLELE`, `VARIANT_CLASS`, `DISTANCE`, `STRAND`, `FLAGS`, `EXON`, `INTRON`, `SYMBOL` and `BIOTYPE` agree on every paired row.

Rows that are not paired: 127 are only fastVEP's (121 on transcripts, 6 for a `*` allele) and 7 only VEP's.
VEP's `--gff` reader keeps only the feature types of `%INCLUDE_FEATURE_TYPES` (`GFF.pm`), its log says it skipped parents "due to invalid types", and `unconfirmed_transcript` and `scRNA`, the types of the three transcripts looked at, are not on that list.
The three of the seven VEP-only rows looked at are an `intergenic_variant` where fastVEP reports one of those transcripts upstream or downstream.
Not all 134 were checked, and none of it is new (see [What is measured](#what-is-measured)).

| What differs | HG002 cells | World cells | Who is right |
|---|---:|---:|---|
| `HGVSc` and `HGVS_OFFSET` of an allele in a record with several ALTs | 1,717 and 2,070 | 108 and 78 | fastVEP: [divergence 11](#11-an-allele-of-a-record-with-several-alts-is-described-as-it-would-be-alone) |
| `HGVSg` of an allele in a record with several ALTs | 35 | 709 | fastVEP: divergence 11 |
| `HGVSc` of an intronic base after the terminator | 12 | 0 | fastVEP: divergence 4 |
| `cDNA_position`, `CDS_position`, `Protein_position` written `N-?` or `?-N` | 9, 0, 0 | 68, 43, 43 | VEP: [a gap](#positions-and-hgvsc-at-the-edge-of-a-transcript-or-an-exon) |
| `HGVSc` of a deletion that runs off the end of the transcript | 9 | 16 | VEP: a gap, in the same section |
| `INTRON` of a change that spans an exon edge | 0 | 16 | VEP: a gap, in the same section |
| `HGVSc` that one tool writes and the other does not | 4 | 8 | the Part 2 `HGVSc` gap, in fastVEP's direction |
| `HGVSc` of an inversion padded with bases it shares | 0 | 15 | fastVEP: [`inv` is the form HGVS asks for](#an-inversion-inside-padding) |
| `Consequence` and `IMPACT` where the change is at the terminator | 0 and 0 | 10 and 5 | not decided: [see below](#a-change-at-the-terminator) |
| `HGVSp` | 17 | 32 | divergences 3, 5 and 7 and the Part 2 `HGVSp` rows |

Every line is the cells of one class, and a row can be in more than one.

### Positions, and `HGVSc`, at the edge of a transcript or an exon

**HG002: 9 `cDNA_position` cells, and the `HGVSc` of the same 9 rows.
The world: 154 position cells, 16 `INTRON` cells and 16 `HGVSc` cells.**

Where one end of a change is not in an exon, Ensembl writes `?` for it: `?-7` for `2:240646354 CTGTCTTCCTGTGACGT>C` on ENST00000847201 in HG002, a 16-base deletion of which Ensembl's description covers 7, and `199-?` for `chr10_2497_CAGCTGC/C` of the world.
fastVEP writes the one coordinate it has, `7` and `199`.
The `HGVSc` of such a change differs in the same rows: Ensembl names the part inside the transcript, `n.1_7del`, and fastVEP names the base at the edge, `n.7del`, or in 10 rows of the world nothing at all.
One base for a change that removed seven is the less accurate of the two, and it is the same defect as the positions.
Where a change overlaps an exon and the intron beside it, Ensembl also writes `INTRON` beside `EXON` (`1/2` for `rs544769` of the world), and fastVEP only `EXON`.
Not closed here.

### An inversion inside padding

**The world: 15 `HGVSc` cells.**

`AAC>AGT` replaces `AC` with `GT`, which is its reverse complement.
HGVS asks for `inv`, and fastVEP writes `c.198_199inv`.
Ensembl writes `c.198_199delinsGT`: `hgvs_variant_notation` tests for the inversion on the alleles as it holds them, `AAC` and `AGT`, which are not reverse complements, and `_clip_alleles` afterwards clips the shared `A` and leaves the type as `delins`.
Both describe the change.
A replacement that is wholly an inversion (`CA>TG`) is `inv` in both tools, as the 121 of the ClinVar sample are.

### A change at the terminator

**The world: 10 `Consequence` and 5 `IMPACT` cells.
HG002: none.**

Ten rows of the world, all within a few bases of the terminator codon of ENST00000000001, ENST00000000002 or ENST00000000003, get a different term in each tool.
`chr10_3245_A/AGCATG` is `frameshift_variant` in Ensembl, with `p.Ter218TrpextTer25`, and `inframe_insertion,stop_retained_variant` in fastVEP, whose five inserted bases are not a multiple of three.
`chr10_10095_TCAGT/GCGCA` is `coding_sequence_variant` (MODIFIER) in Ensembl and `stop_lost` (HIGH) in fastVEP.
It was not diagnosed, so it is listed and not argued either way: the in-frame reading of a 5-base insertion looks wrong, and nothing here shows that Ensembl's is right.

### What a `*` allele and a symbolic one get

fastVEP writes no `HGVSg` for them, as Ensembl does not (`next if $allele =~ m/[^ACGT\-]/ig` in `hgvs_genomic`), and used to write `chr10:g.1400?`, which is not HGVS.
Ensembl writes no row for a `*` allele at all, and none for a record whose only ALT is `*`; fastVEP writes one, with an empty `Consequence` and `MODIFIER`, and 6 rows of HG002 are those.
An ALT of `N` gets no `HGVSg` in either tool.
A VCF REF of `N` gets the genome's base in both tools' `HGVSg`, and in Ensembl's `HGVSc` (`c.201+300C>A`) but not in fastVEP's, which writes the REF it was given (`c.201+300N>A`).

### Lone indels that are not minimal

Ensembl minimises a record with one ALT that is an indel from both ends (`minimise_alleles`), where fastVEP strips only the base the alleles share at the front.
The two are the same for a record that is already minimal and padded by one base, which is every indel of HG002 (0 of its 2,818 are not), and they differ for one that shares bases at the back as well, or two or more at the front.
`Allele`, `Location`, `REF_ALLELE`, `VARIANT_CLASS` and `DISTANCE` then differ, `HGVSc` and `HGVSg` do not (both clip), and some rows get a different term.
The world has none of them as it is, and the same world with each ALT its own record has 197 rows of 7,213, in 93 records: for `1306 AG>ACTATAATG` Ensembl writes `Allele` `CTATAAT` at `1306-1307` as an `insertion`, and fastVEP `CTATAATG` at `1307` as an `indel`.
8 `Consequence` cells differ on those rows (4 of them `IMPACT`), `intron_variant` MODIFIER in Ensembl and `splice_donor_variant` HIGH in fastVEP for an insertion next to a donor that fastVEP reads as replacing it.
Not closed here: it changes what the whole annotation reads, and wants its own measurement.


---

## What is *not* a divergence

Four things that look like one and are not.

**Consequences are not 3'-shifted.**
`TranscriptVariation.pm` (l. 132) sets
`$self->{shifted} = (defined($args{'-no_shift'}) && !$args{'-no_shift'})`, so the flag is false
unless `--no_shift 0` is passed explicitly, and `get_all_OverlapConsequences`
(`BaseVariationFeatureOverlapAllele.pm`, l. 269) applies the shift offset to a predicate's
coordinates only when it is set.
HGVS *is* shifted, by a separate call to `_return_3prime`.
fastVEP does the same.
When this was tried the other way round, shifting the span before the splice predicates took that
run from 365 mismatched rows to 2,262.

**Splice sites are matched against the bases that differ, not the variant's span.**
`_get_differing_regions` (`VariationFeatureOverlapAllele.pm`, l. 405) XORs the two allele strings
position by position with the shorter padded, and groups the result into contiguous runs.
That is both narrower than the span - a matching interior is not tested - and wider, because the
padding puts every base past the reference allele's end into a region.
fastVEP reproduces it, which is why the splice terms disagree on exactly one row in 151,684.

**An HGVSc offset does not say where the variant is, in either tool.**
`c.` is a display form and 3'-shifted, and the shift runs away from the exon on a donor and toward
it on an acceptor, so the `+N` / `-N` token can be a long way from the change that earned the
consequence term beside it.
Real VEP 115.1 writes `ENST00000379370.7:c.4298+21_4298+55del` for an AGRN deletion it calls
`splice_donor_variant`, and `ENST00000676179.1:c.2043-9dup` for a KIF1B insertion that earns no
splice term at all.
fastVEP writes the same strings, because they are the right descriptions.

What follows is that no consumer of these fields may recover a position by parsing one, and that
"where the change is" and "where HGVS numbers it" are two different questions rather than one
question with a right and a wrong answer to it.
Both are computed off the transcript now, and the ACMG criteria take whichever their own rule
asks for.
BP7's deep-intronic extension reads where the change sits, so that it cannot call a variant
deep-intronic while the row beside it says `splice_donor_region_variant`.
PVS1, PS1's splice path, PM1, PM4, BP3 and BP4's splice path read where HGVS *places* it: the clipped alleles are shifted
3' over the reference, bounded by the transcript as the description's shift is, and the result is
classed as wholly intronic beyond `+-2`, exonic in or out of frame, or neither.
On a donor the 3' end of an indel's range is the furthest into the intron it can be written, so a
placement past `+2` *proves* the canonical `GT` survives the edit, and no exon base is edited
either; on an acceptor the shift runs toward the exon, so a deletion written across the acceptor
can arrive wholly exonic and in frame (MLH1 `TAGG>T`, or FOXP2 `CCAG>C` down a `CAG` repeat).
Neither case changes a consequence term, which is left as the predictor reports it: the placement is used only where a
criterion is asking about the edited sequence and not about the record.
No criterion parses an HGVSc for it.
[`analysis/acmg_benchmark/RUN_VERSIONS.md`](../analysis/acmg_benchmark/RUN_VERSIONS.md) carries
the worked examples and the benchmark runs on either side of the choice.

**`splice_region_variant` is decided by the last differing region, not the union.**
`_intron_effects` (`BaseTranscriptVariationAllele.pm`, l. 215) assigns rather than or-assigns it
inside its region loop, so a change whose first differing base is in the splice region and whose
last is not comes out with no splice term.
Reproducing the tidier rule instead put a `splice_region_variant` on 63 rows of that run which VEP
calls plain missense.

---

## Fixed since this document was last measured

Ten defects are gone, and the numbers above are what is left.
Each was measured the same way before and after, on the same samples: the row counts below are
what the fix removed.

| What was wrong | Rows it cost | What it was |
|---|---:|---|
| `HGVSg` was the plain description, where Ensembl's 3'-rule over the genome moves an indel and names a repeat a `dup` | 10 of 13,535 records of HG002, against 115.1 | An insertion sat where the VCF put it and never collapsed to a `dup` (`9:905488 C>CTGTGTGTG` was `g.905488_905489insTGTGTGTG` for Ensembl's `g.905507_905514dup`). `hgvsg_vep` transcribes `hgvs_genomic`; what is left is divergence 11 |
| The protein in front of every `HGVSp` was version 1: `ENSP00000363557.1:p.Thr266Met` for Ensembl's `.3`, and `NP_000001.2.1` for a RefSeq protein that already carries its version | 184 of the 844 `HGVSp` rows of HG002 against 116.3, and every RefSeq protein but the ones at version 1 | The GFF3's CDS records carry the version (37,293 of the 234,024 proteins in Ensembl 115's file have one other than 1) and VEP 116 reads it (`$ordered_cdss->[0]->{attributes}->{version} || 1`; 115 gave every translation version 1). `GffCds` parsed it and nothing read it. Read now, and appended to an id only if the id has none (`Transcript::versioned_protein_id`); the cache format moved on so that no old one is trusted |
| Each allele of a multi-allelic record kept the bases it shared with the reference | 801 of 824 genome-wide `HGVSc` rows, and 1 of 11 `HGVSg` records | A site is trimmed once for the whole record, so an allele that is a clean deletion arrives as a replacement and was written as one: `c.553-61775_553-61770delinsCTT` for TTC28's `c.553-61772_553-61770del`. Clipped now where the description is built, which is the second of the two places Ensembl clips. What the clip does *not* reach is in Part 2 |
| The HGVSc 3'-shift ran on the *spliced* sequence | 563 of 602 `HGVSc` rows | It could not follow a repeat out of an exon into the intron (`c.220del` for VEP's `c.220+1del`), and it would walk a deletion over a splice junction into the next exon, naming a block that is contiguous in the mRNA and not in the DNA a `c.` description is numbered against (BRCA1 `c.71_81del` for VEP's `c.70_80del`). The shift now runs on the genome, bounded by the transcript, and each end is written where it lands |
| A lost terminator was written `p.Glu486ext*?` | 30 `HGVSp` rows | The extension is a fact about the edited transcript, and `ext*?` is the form HGVS keeps for one nobody can measure. Now divergence 7 |
| The 3'-rule's rotation carried residues out from behind a terminator | 23 `HGVSp` rows | `insValAlaLeuAspThrTerVal` named a Val the protein never has |
| `ref_eq_alt_sequence`'s second clause was read as a position test | 18 consequence rows | It is a comparison of sequences: an insertion whose residues repeat what follows them holds the clause anywhere the repeat reaches the terminator, not only on the last residue (MSH6 `c.4106_4108dup`, ENST00000936511) |
| `start_lost` stopped at its first test | 4 rows of the in-frame deletion set, all four `IMPACT` | Ensembl falls through to the peptide test when the coordinate test declines because the edit is an in-frame indel, and an in-frame deletion that removes the initiator is a start loss there (KCNA2 `c.3_11del`, ENST00000639048, HIGH not MODERATE) |
| An unknown residue silenced `missense_variant` | 1 consequence row, `IMPACT` | `QX/HX` resolves one residue and not the next; Ensembl reports `missense_variant,coding_sequence_variant`, and reporting the second alone made a real missense MODIFIER (IL7R ENST00000515665) |
