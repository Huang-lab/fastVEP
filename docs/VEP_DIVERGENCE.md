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

That is release 115.1.
Against 116.3, which is what most people run, `HGVSp` has its own measurement in [Part 3](#hgvsp-at-a-start-or-a-stop-codon-in-a-repeat-and-at-the-end-of-a-protein).
None of the 186,055 rows of HG002 differs, and 1,818 cells of the 3,136,295 rows of the ClinVar sample do, of which 1,492 are the `--gff` reader's defect described in divergences 5 to 8 and the other 326 are listed there one by one.

---

## Part 1 - deliberate divergences

Each is a case where Ensembl's output is either not valid HGVS, or contradicts what the sequence
says, and where matching it would degrade a clinical call.

### 1. A change that introduces a premature stop is HIGH, not a moderate in-frame insertion

**Against release 116.3 this is no longer a divergence.**
`ref_eq_alt_sequence` was rewritten in release 116 and its first clause, quoted below, is gone.
VEP 116.3 writes `stop_gained,frameshift_variant`, HIGH, for all five variants of the table, and agrees with fastVEP on every `Consequence`, `IMPACT`, `Amino_acids` and `Codons` cell of the 6,021 ClinVar insertions and deletions of the five genes they are in (190,397 rows), and on a 1-in-5 sample of ClinVar 2-star+ (3,136,295 rows) but for the 17 rows of [Part 3](#coding-consequence-terms).
What follows is what release 115.1 did, which is why fastVEP refused the clause, and what the numbers below were measured against.

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

**Against release 116.3 this is no longer a divergence, because fastVEP no longer makes it.**
fastVEP asks the edited sequence now, as `_ins_del_start_altered` does, and the edited sequence of MLH1 `3:37014476 AAATG>A` reads `ATG` at the initiator again: the deleted `AATG` is one copy of a repeat (`A AATG AATG G`), and deleting either copy leaves the same transcript.
It is `start_retained_variant,5_prime_UTR_variant`, LOW, in VEP 116.3 and in fastVEP, where it was `start_lost`, HIGH.
The change is real and so is its cost: the variant is a frameshift after Met1 (`p.Asn2ValfsTer10`, in VEP's own `HGVSp`) and carries no HIGH term, because a deletion that starts in the 5' UTR has no CDS start and Ensembl's `frameshift` needs one.
Nine of the ten rows of the 6,600-variant sample are this variant across MLH1's transcripts.
What follows is what release 115.1 did, and what the numbers below were measured against.

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

**60 `HGVSp` rows against 115.1.
Against 116.3: 162 of the ClinVar sample, 5 of the generated world and 2 of the edge set.**

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

**Against release 116.3 this is not a difference with VEP, and the mechanism this section used to give was wrong.**
It is a defect of VEP's `--gff` reader on the reverse strand, and a cache does not have it (see [Part 3](#hgvsp-at-a-start-or-a-stop-codon-in-a-repeat-and-at-the-end-of-a-protein)).
On the ClinVar sample 1,492 `HGVSp` cells differ from stock VEP 116.3 in a count or a residue of a frameshift or an extension, every one of them on a reverse-strand transcript, and none differs once the reader is corrected.
What follows is what release 115.1 wrote and what the table was measured against.

**411 `HGVSp` rows** against 115.1, the largest `HGVSp` divergence of all.

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

The reason Ensembl's count is wrong is that it reads that UTR from the wrong place.
This section used to say that `_get_alternate_cds` builds its edited CDS from 3'-shifted coordinates
while taking consequences from unshifted ones, and that the two do not line up.
That was not checked against the code, and it is not what the code does.
`_get_alternate_cds` appends `$self->transcript_variation->_three_prime_utr`, which is
`BaseTranscriptVariation::_generic_utr`, which transfers the transcript to its feature slice and
asks it for its UTR, and for a transcript on the reverse strand built from a GFF3 that returns a
sequence of the right length from the wrong place in the genome ([Part 3](#coding-consequence-terms)
shows the printed sequence).
Every transcript the tables of divergences 5 to 8 name is on the reverse strand but one, DNAH10 in divergence 7, and that is the one row on which VEP agreed.

Four of the six rows were run again on 116.3, stock and with `_generic_utr` taking the UTR from `spliced_seq`, as the core API's own accessors do:

| Gene | Transcript | Stock VEP 116.3 | VEP 116.3 with the UTR corrected | fastVEP |
|---|---|---:|---:|---:|
| ARSB | ENST00000264914 | 175 | **48** | 48 |
| VWF | ENST00000321023 | 56 | **302** | 302 |
| CFI | ENST00000882820 | 22 | **14** | 14 |
| NTHL1 | ENST00000651522 | `Ter?` | **55** | 55 |

Stock 116.3 writes the numbers 115.1 did, and corrected it writes fastVEP's, on all four transcripts (and on the 15 others of CFI and 17 others of NTHL1 that carry the variant).
MPV17 and PCDH15 were not run again.

### 6. A frame with no stop in the transcript is reported as having none

**Against 116.3 this is divergence 5 again**: the NTHL1 row above is `Ter?` in stock 116.3 and 55 with the UTR corrected, and on the ClinVar sample the `Ter?` of stock VEP against a distance is 106 cells, every one on the reverse strand, and the converse (60 cells) the same.

**51 `HGVSp` rows** against 115.1, all of them `Ter?` on one side and a distance on the other.

The NTHL1 row in the table above is one of them, checked the same way: the stop really is 55
residues along, and `Ter?` is not an answer about it.
Where fastVEP writes `Ter?` instead - ITGA2B `p.Ter1040Trpfs` on ENST00000262407, CLN3
`p.Gly37Valfs` on ENST00000561505, RUNX1 `p.Arg380Profs` on ENST00000399240 - the shifted frame
runs to the end of the transcript without another stop, and there is nothing to count.

### 7. A lost terminator's extension is counted, not written off

**Against 116.3 this is divergence 5 again**: 205 cells of the ClinVar sample differ from stock VEP in the count of an extension and 38 are `Ter?` on one side, all of them on the reverse strand, and none differs with the UTR corrected.

**20 `HGVSp` rows on the ClinVar sample, 1 genome-wide** against 115.1.

`p.Ter486GluextTer36` says the terminator at 486 became Glu and the protein runs 36 residues
further before the next stop.
Both halves are in the sequence, so both are computed:

| Gene | Transcript | Variant | Computed | fastVEP | Ensembl VEP 115.1 |
|---|---|---|---:|---:|---:|
| TP53 | ENST00000620739 | `17:7669611 A>T` | **9** | `p.Ter355ArgextTer9` | `p.Ter355ArgextTer1` |
| MUTYH | ENST00000674679 | `1:45333589 TC>AA` | **33** | `p.Ter48PheextTer33` | `p.Ter48PheextTer19` |
| DNAH10 | ENST00000538983 | `12:123932637 A>G` | **3** | `p.Ter82TrpextTer3` | `p.Ter82TrpextTer3` |

VEP's distances here fail for the same reason as divergence 5, which is the UTR it reads and not
a shift: on the TP53 variant it returns 1, 2, 5, 11, 30 and 106 on six different transcripts of
the same gene, all of which reconstruct to 9, because each of the six has a different 3' UTR
and VEP reads each of them from somewhere else.
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

Both are on the reverse strand, and so is the cause: with the UTR corrected, no first residue of the ClinVar sample's frameshifts differs from VEP's but the two mitochondrial rows of divergence 9, where VEP translates the edited sequence with the standard code on chrM.
Against 116.3 this is divergence 5 as well.
There is one more way a first residue can differ, and it is not about the UTR: [divergence 12](#12-an-insertion-the-3-rule-carries-further-than-it-is-long-on-the-reverse-strand-is-described-as-another-insertion).

### 9. The smaller ones

| What | Rows | fastVEP | Ensembl VEP 115.1 | Why ours |
|---|---:|---|---|---|
| A two-residue deletion in a poly-Glu C-terminus. NT5C2 `10:103089678 TTCCTCC>T`, ENSP00000502205, protein 561 aa ending `EEEEEEE` | 43 of the in-frame deletion set, on one variant; 31 of the ClinVar sample against 116.3 | `p.Glu560_Glu561del` | `p.Glu559_Glu560del` | The 3'-most placement is the one HGVS asks for. Ensembl stops short of the end of the protein in two ways, both of which are in release 116 as they were in 115: `_shift_3prime` halts *n*-1 residues before the terminus for an *n*-residue change ([#94](https://github.com/Huang-lab/fastVEP/issues/94)), and the residues it looks at come from `_get_surrounding_peptides`, which returns nothing when `length <= position`, so a change that ends one residue short of the end cannot see the last one and is not shifted at all (NT5C2 is that case) |
| A transcript whose CDS is annotated 5' incomplete (first CDS record carries phase 1 or 2, so no initiator is annotated). TCIRG1 ENST00000698256, NF1 ENST00000696141, NFIA ENST00000496712 | 3 | no `start_lost` | adds `start_lost` | There is no `ATG` to lose. Declining is the lower-impact answer, which is worth saying plainly |
| `start_lost` and `start_retained_variant` on the same row, where a *substitution* replaces a non-canonical initiator with `ATG`. NBEA `13:35171273 C>T` (`aCg`→`aTg`), ENST00000629018 | 2 | `start_lost` | `start_lost,start_retained_variant` | **Gone against 116.3**: the window models `_snp_start_altered` now and both tools write the pair, and `p.Thr1?`. The pair contradicts itself, and `IMPACT` is HIGH either way |
| The polypyrimidine tract on a transcript carrying a frameshift intron. PTEN `10:87933000 C>CT`, ENST00000693560 | 1 | `splice_polypyrimidine_tract_variant` at -3 to -17 | the term stops at -13 | **Gone**: every exon of such a transcript is 12 bases wider for the `exon` test, as in `_overlapped_exons`, and an intron of 13 bases or fewer is one, which is not what this row said (see [Frameshift introns](#frameshift-introns)) |

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

### 12. An insertion the 3'-rule carries further than it is long, on the reverse strand, is described as another insertion

**Against 116.3: 49 `HGVSp` cells of the edge set, 213 of 1,884,993 rows of VEP's own Ensembl cache run and 24 of 203,125 of its RefSeq cache run.
None on HG002 or on the ClinVar sample.**

To describe the protein, Ensembl moves an insertion to its most 3' position (`_return_3prime`) and then turns the inserted string, with `shift_feature_seqs`, so that it names what is inserted *there*: one base round for each base the change moved.
The turning is a loop of `seq_length - shift_length` steps on the reverse strand.
That is the right number of turns when the change has moved no further than the insertion is long, and a negative number, which is none, when it has moved further: `AGA` that has moved four bases along a repeat is still `AGA`, where the insertion it stands for is `AAG`.
The edited CDS that the peptides and the counts are read from is then made with a string that is not the variant's, and the notation names a different insertion.

Printed from inside VEP 116.3 for `chr10_10108_A/AAGA` on a reverse-strand transcript: shift 4, allele `AGA`, and `_get_alternate_cds` is handed `TCT` (the reverse complement of `AGA`) at the moved position, where the string turned four times is `CTT`.
VEP writes `p.Leu199_Ala200insSer` (`TCT` is Ser) and fastVEP `p.Leu199dup` (`CTT` is Leu, which Leu199 repeats).
VEP's own `HGVSc` for that row says it is a duplication (`c.595_597dup`).
The same on `21:33465378 C>CACC` of ENST00000420455 in the Ensembl cache run, shift 4 of a three-base insertion: VEP `p.Asp83delinsGlyTyr`, fastVEP `p.Gly82dup`, and VEP's `HGVSc` `c.245_247dup` with a `HGVS_OFFSET` of -4.

Reproducing the loop (a throwaway change: no turning on the reverse strand when the shift is longer than the insertion) makes 49 of 49 cells of the edge set agree, which is how the cause was told from the others.
fastVEP turns the string the number of times it moved, on either strand.

### 13. A change at a stop codon that moves the frame, and an incomplete or ambiguous codon: an unresolved residue is not a terminator

**Against 116.3: 61 `HGVSp` cells of the edge set, 12,356 of 1,884,993 rows of VEP's Ensembl cache run, 5,169 of 203,125 of its RefSeq cache run, and 122 rows of the ClinVar sample.**

Ensembl rewrites the first `*` of each peptide to `X` before it tests them, and then rewrites every `Xaa` to `Ter`.
The `*` becomes a residue so that the lengths can be compared; the `X` that the window can contain, from a codon that is not complete or that holds an `N`, is then a `Ter` too.
Two things follow, and both are measured.

- A change at a stop codon that moves the frame (a deletion of one base, an insertion of two) leaves a window that ends part-way through a codon, so its last residue is `X`, and the notation reads it as a stop.
  Release 116 makes it worse: `frameshift` answers no for a change that starts on the terminator, so the notation is not even the frameshift's.
  `21:25945913 CT>C` of ENST00001003897 is `stop_lost`, HIGH, in VEP's own `Consequence` and `p.Ter583=` in its `HGVSp`: the protein loses its stop and the notation says nothing changed.
  VEP's `Consequence` says `stop_lost` on every one of the rows of this kind in the Ensembl cache run, and its `HGVSp` is `p.Ter583=` on 9,166 of them and `p.Ter583delinsTrpTer` or its like on 3,190.
- A transcript whose CDS ends part-way through its last codon (`incomplete_terminal_codon_variant`), or a codon with an `N` in it, is `Xaa` and not `Ter`: `p.Ter192=` for ClinVar's `17:65558046 G>A` on ENST00000585045, where the residue is the one the CDS does not complete.
  111 rows of the ClinVar sample are these, 54 of them an incomplete terminal codon and 57 a codon with an `N`, and 11 more are the first kind of change.

fastVEP does not read an unresolved residue as a stop.
`Ter` is a `*` and `Xaa` an `X`, throughout.
For a change that moves the frame (its length against its span in the CDS), at the terminator or not, and whose window reaches no stop, it describes the change from the edited translation, as it does a frameshift: `p.Ter583TrpextTer29`, where the stop became a Trp and the protein runs 29 residues on, and `p.Ter583=` only when the edited translation keeps the stop in the same place.
Where the window does reach a stop, nothing after the stop is protein, and the notation is Ensembl's own.

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

None that was measured against release 116.3: the rows of [Part 3](#coding-consequence-terms) are what is left.
The one gap this section listed against 115.1, `3_prime_UTR_variant` and `5_prime_UTR_variant` missing on a change that spans the end of a transcript whose coding region reaches it (CHEK2 `22:28711907 ACT>A`, ENST00000439200), was Ensembl's `_before_coding` and `_after_coding` being plain `overlap` calls against an interval that can be empty.
It is `stop_retained_variant,3_prime_UTR_variant` in both tools now.

### HGVSc

None that was measured against release 116.3.
The one gap this section listed against 115.1, fastVEP naming a change where VEP names none (18 rows of the ClinVar sample, 11 genome-wide and 9 of the in-frame deletion set: `c.-66_-60del` for HNF4A `20:44401297 GGGAGGGC>G` on ENST00000415691, a deletion that straddles the transcript's first base), is closed.
Ensembl clamps a change that runs past an end of a transcript to the transcript, works out the 3'-rule for the change as it was and adds that distance to the clamped coordinates, and writes nothing for one the shift then carries off the end ([Part 3](#hgvsc-at-the-edge-of-a-transcript)).
HNF4A is `c.-63_-60del` in both tools now.
What differs from 116.3 now is divergence 4 (12 cells of HG002, 308 of the ClinVar sample, 390 of VEP's Ensembl cache run), divergence 11, and the 11 cells of the edge set at the ends of a transcript that Part 3 lists.

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

None that was measured against release 116.3, in either direction.
The five kinds of gap this section listed against 115.1 are gone.
The variants that stood for four of them (MSH6, TRAF3IP1, IL7R and ATM, and HNF4A for `HGVSc`) agree with VEP 116.3 on every `HGVSc` and `HGVSp` cell of all 115 rows they have, and the ClinVar sample has no row left in which one tool names a protein change and the other does not.

| Gap against 115.1 | Rows | What closed it |
|---|---:|---|
| fastVEP names one where VEP names none | 196 on the ClinVar sample, 17 genome-wide | 177 of them a frameshift that also earns `splice_region_variant` (MSH6 `2:47806357 TG>T`, `c.3828+1del`). Ensembl describes the change after the 3'-rule has moved it, and a change the rule carries out of the coding sequence has no protein position |
| VEP names one where fastVEP does not | 27 | TRAF3IP1 `2:238332821 TAGTC>T`, `p.Ser306GlnfsTer34` in both. The same rule from the other side: the notation is written for the change where it is moved to, and fastVEP did not have a protein position for the one it was given |
| A frameshift in the first codons, which VEP writes as unknown | 5 | The notation of a change that loses the initiator is the unknown one (`p.Met1?`), whatever else the change is, and fastVEP writes it. What differs now is the `Consequence` of a transcript annotated 5' incomplete (divergence 9), which has `start_lost` in VEP and not here: 9 rows of the ClinVar sample |
| `Amino_acids` of `X` from an ambiguous reference codon | 3 | Not closed, and not a gap: `p.Xaa157=` is what the codon is, and [divergence 13](#13-a-change-at-a-stop-codon-that-moves-the-frame-and-an-incomplete-or-ambiguous-codon-an-unresolved-residue-is-not-a-terminator) says why |
| A delins naming an unknown residue the window reaches | 2 | IL7R `5:35860925 GC>TT` is `p.Gln52His` in both, and ATM `11:108310201 A>AT` `p.Ter1907LeuextTer1121`: the residue both windows share is clipped before the notation is written |

---

## Part 3 - the tab output, against VEP 116.3

Everything above was measured against Ensembl VEP 115.1.
The tab output's optional columns, `HGVSg` and `HGVS_OFFSET` had no measurement at all, and the release most people run now is 116, so they have their own, against **Ensembl VEP 116.3** (Docker `ensemblorg/ensembl-vep:release_116.3`).
It reads the same GRCh38 Ensembl 115 GFF3 and FASTA as above and is run with `--tab --symbol --biotype --canonical --numbers --variant_class --show_ref_allele --hgvs --hgvsg`.
`validation/compare_tab.py` compares its file with fastVEP's cell by cell, pairing rows on (`Uploaded_variation`, `Allele`, `Feature`).

| Input | Records | Paired rows | Cells compared | Cells that differ |
|---|---:|---:|---:|---:|
| GIAB HG002 WGS, 1-in-200 systematic (the sample above) | 20,241 | 186,055 | 5,023,485 | 3,834 (0.076 %); 3,867 before the `HGVSp` and the transcript-end `HGVSc` |
| A generated world: both strands, coding and non-coding genes, 372 of its 3,000 records with several ALTs | 3,000 | 7,213 | 216,390 | 936 (0.43 %); 973 before |
| The same world's edges: every shape of change at every position within ten bases of every exon, CDS and transcript boundary, each its own single-ALT record (`fixtures/vep_tab_edges`) | 17,672 | 36,958 | 1,108,740 | 263 (0.024 %); 2,053 before the HGVS notation, 2,689 before the coding terms and 7,828 before the position fix |
| The same world's lone indels written with shared bases (`fixtures/vep_tab_lone_indels`) | 1,500 | 3,176 | 95,280 | 9 (0.009 %); 52 before the HGVS notation, 59 before the coding terms and 14,526 before the minimisation |
| A 1-in-5 sample of ClinVar 2-star+, without the HGVS columns | 134,732 | 3,136,295 | 72,134,785 | 102 (0.00014 %): 20 once a frameshift intron is read as Ensembl reads it and 82 that follow from the CDS flags of VEP's cache (below); 43 and 275 before |
| The same sample with them (`--hgvs --hgvsg`) | 134,732 | 3,136,295 | 84,679,965 | 2,263 (0.0027 %), 760 once VEP's UTR is corrected (below); 3,354 before the HGVS notation |

Since the change that reads a transcript as VEP's cache holds it ([Against VEP's own cache](#against-veps-own-cache)), the comparisons of HG002 and of the ClinVar sample against a `--gff` run leave out three columns, `SYMBOL_SOURCE`, `HGNC_ID` and `FLAGS` (`--ignore`): a GFF3 carries none of them, VEP's `--gff` writes `-`, and fastVEP now follows the cache.
That is 27 columns of 30, and 23 of 26 without the HGVS columns, which is what the cell counts above are; the generated worlds have no such metadata and agree in all 30.

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

On HG002 the 3,867 cells are in four columns: `HGVS_OFFSET` 2,073, `HGVSc` 1,742, `HGVSp` 17 and `HGVSg` 35.
`Consequence`, `IMPACT`, `Amino_acids`, `Codons`, `cDNA_position`, `CDS_position`, `Protein_position`, `Location`, `REF_ALLELE`, `VARIANT_CLASS`, `DISTANCE`, `STRAND`, `FLAGS`, `EXON`, `INTRON`, `SYMBOL` and `BIOTYPE` agree on every paired row.

Rows that are not paired: 121 are only fastVEP's, all on transcripts (it also wrote 6 for a `*` allele, which VEP does not, and no longer does), and 7 are only VEP's.
VEP's `--gff` reader keeps only the feature types of `%INCLUDE_FEATURE_TYPES` (`GFF.pm`), its log says it skipped parents "due to invalid types", and `unconfirmed_transcript` and `scRNA`, the types of the three transcripts looked at, are not on that list.
The three of the seven VEP-only rows looked at are an `intergenic_variant` where fastVEP reports one of those transcripts upstream or downstream.
Not all 128 were checked, and none of it is new (see [What is measured](#what-is-measured)).

| What differs | HG002 cells | World cells | Who is right |
|---|---:|---:|---|
| `HGVSc` and `HGVS_OFFSET` of an allele in a record with several ALTs | 1,717 and 2,070 | 108 and 78 | fastVEP: [divergence 11](#11-an-allele-of-a-record-with-several-alts-is-described-as-it-would-be-alone) |
| `HGVSg` of an allele in a record with several ALTs | 35 | 709 | fastVEP: divergence 11 |
| `HGVSc` of an intronic base after the terminator | 12 | 0 | fastVEP: divergence 4 |
| `HGVSc` of a deletion that runs off the end of the transcript | 9, now none | 16, now none | was a gap, [closed](#hgvsc-at-the-edge-of-a-transcript) |
| `HGVSc` that one tool writes and the other does not | 4, now none | 8, now none | was the Part 2 `HGVSc` gap, in fastVEP's direction; [closed](#hgvsc-at-the-edge-of-a-transcript) |
| `HGVSc` of an inversion padded with bases it shares | 0 | 15 | fastVEP: [`inv` is the form HGVS asks for](#an-inversion-inside-padding) |
| `Consequence` and `IMPACT` of a change at a start or stop codon of a reverse-strand transcript | 0 and 0 | 0 and 0 (9 and 3 on the edge set) | VEP: its `--gff` reader builds the UTR from the wrong place, [see below](#coding-consequence-terms) |
| `HGVSp` | 17, now none | 32, now 23 | divergences 3, 12 and 13, and VEP's `--gff` reader: [below](#hgvsp-at-a-start-or-a-stop-codon-in-a-repeat-and-at-the-end-of-a-protein) |

Every line is the cells of one class, and a row can be in more than one.

### Positions and `EXON`/`INTRON` at the edge of an exon

Before: **HG002 9 `cDNA_position` cells.
The world 154 position cells and 16 `INTRON` cells.
The edge set (below) 1,784 `cDNA_position`, 1,268 `CDS_position`, 1,229 `Protein_position`, 255 `EXON` and 603 `INTRON` cells.
After: none.**

Where one end of a change is not in an exon, Ensembl writes `?` for it: `?-7` for `2:240646354 CTGTCTTCCTGTGACGT>C` on ENST00000847201 in HG002, a 16-base deletion of which the transcript has 7, and `199-?` for `chr10_2497_CAGCTGC/C` of the world.
fastVEP wrote the one coordinate it had, `7` and `199`, which names one base of a change that covers sixteen.

The rule is Ensembl's mapper's, read from `Mapper::map_coordinates`, `map_insert`, `TranscriptMapper::genomic2cds` and `BaseTranscriptVariation::cdna_start` of release 116 and checked against VEP 116.3 itself:

- The span is mapped through the exons in transcript order, which gives coordinates and gaps, and the cell is the start of the first piece and the end of the last (`format_coords`).
  A gap is any stretch that is not exon, so an end has a coordinate exactly when the 5'-most base (for the start) or the 3'-most (for the end) is in an exon.
  On the reverse strand the ends swap: the genomically lower end of a change is its 3' end, so the deletion of `2491-2502` on ENST00000000005, whose exon starts at 2500, is `200-?` and not `?-200`.
- The CDS and the protein are the same test against the coding sequence, with the UTR part of a piece made a gap (`genomic2cds`): a deletion that begins in the 5' UTR and ends in the initiator is `?-3` and `?-1`.
- An insertion is mapped by its two flanking bases with the gaps dropped, and a flank in an exon contributes its own position and the one across the boundary, whether that exists or not: an insertion in the intron right after exon 1's last base (cDNA 601) is `601-602`, and so is one right before exon 2's first base.
  For the CDS the piece `(c + 1, c)` is a gap when `c` is before the first coding base or `c + 1` after the last, so an insertion between the last UTR base and the initiator, or after the terminator, has no CDS or protein position.
- `EXON` and `INTRON` list every exon and intron that the change overlaps by `overlap($vf_start, $vf_end, ...)`, which for an insertion (`start = end + 1`) is true only when it sits between two bases of the feature.
  A change over two exons is `1-2/3` and the intron between them `1/2`; an insertion at an exon's edge is in neither, and VEP writes `-` for both.

`fixtures/vep_tab_edges` is the evidence: 17,672 single-ALT records at every position within ten bases of every exon and CDS boundary of the world, on both strands, in coding and non-coding transcripts (SNVs, MNVs, deletions of 1 to 30 bases, insertions, delins, and deletions of a whole exon and of a whole intron), and what VEP 116.3 wrote for them (36,958 rows).
After the change every cell of the six columns agrees on every row; the committed fixture is a sample of 2,271 of them, and `vep_tab_edges_golden.rs` holds fastVEP to it.
The JSON writer follows VEP's, which splits the cell and leaves out a `?`: `cdna_end` without a `cdna_start`.

### `HGVSc` at the edge of a transcript

**HG002: 9 rows, and 4 more that one tool wrote and the other did not.
The world: 16 and 8.
The edge set: 973 `HGVSc` cells and 230 `HGVS_OFFSET`.
After: none, none, and 11 and 3.**

The `HGVSc` of a change that has an end outside the transcript differed in the same rows as the positions did: Ensembl names the part inside the transcript, `n.1_7del`, and fastVEP named the base at the edge, `n.7del`, or in 10 rows of the world nothing at all.
Four rules of `hgvs_transcript`, each read from release 116 and checked on VEP 116.3 with its code instrumented, make up what Ensembl does, and fastVEP now does the same:

- HGVS is written only for a change that is `within_feature`, a plain `overlap` with the transcript.
  An insertion is `start = end + 1`, so it overlaps only when both of the bases it sits between are inside, and one in front of the first base or behind the last is an `upstream_gene_variant` or a `downstream_gene_variant` and has no `HGVSc`.
  The 3'-rule used to carry such an insertion into the transcript, and `n.3_4dup` was written for a change that VEP says nothing about.
- A change that overlaps and runs past an end has its coordinates clamped to the transcript (`_var2transcript_slice_coords`: `1` and the transcript's length).
  The reference allele is what the clamped span holds and the alternate keeps its whole length, so a replacement of the transcript's first base and the base before it by `GCG` is `n.1delinsGCG`.
- The 3'-rule is worked out for the change *as it was*, over the chromosome, and its length is added to the clamped coordinates.
  A deletion of twelve bases of which four are in the transcript, in a repeat that carries it three bases on, is `n.4_7del` with a `HGVS_OFFSET` of 3, where shifting the four that are inside gives `n.1_4del`.
  It is not a description an HGVS reader would choose, and it is the one a replacement of VEP's output has to reproduce.
  The walk reads 1,000 bases of the chromosome either side of the change whatever the transcript is.
- A description the shift carries past the last base of the transcript is not written (`return undef if $slice->length < $slice_end + $offset`): a deletion in the last bases of a transcript, in a repeat that runs on beyond it, has no `HGVSc` in VEP.

HNF4A `20:44401297 GGGAGGGC>G` is `c.-63_-60del` in both tools now.

What is left is 11 `HGVSc` cells and 3 `HGVS_OFFSET` of the edge set, of one shape: a replacement whose ends straddle the first or the last base of a transcript, which reaches Ensembl's arithmetic where it is least defined.
VEP writes the range of an insertion past the 3' end of a reverse-strand transcript the wrong way round, `c.*349_*348insGT` where fastVEP writes `c.*348_*349insGT`, and writes an insertion (`c.*349_*348insG`, `c.-400_-399insG`) where fastVEP finds a duplication (`c.*348dup`, `c.-397dup`), with no `HGVS_OFFSET` where fastVEP has one.
A range whose end is before its start is not HGVS, and the duplication is the form HGVS asks for, so these are not reproduced; none of them is on HG002, the ClinVar sample or the Ensembl and RefSeq cache runs.

### `HGVSp` at a start or a stop codon, in a repeat and at the end of a protein

**`HGVSp` cells that differ from stock VEP 116.3, before and after:**

| Input | Before | After | After, with VEP's UTR corrected |
|---|---:|---:|---:|
| HG002, 1-in-200 (186,055 rows) | 17 | 0 | 0 |
| The generated world (7,213 rows) | 32 | 23 | 10 |
| The edge set (36,958 rows) | 838 | 237 | 113 |
| The lone-indel set (3,176 rows) | 30 | 9 | 1 |
| The ClinVar sample (3,136,295 rows) | 2,837 | 1,818 | 326 |
| VEP's own Ensembl cache, chr21 and chr22: 45,916 changes at the start, the stop, every exon boundary and the ends of every transcript, the 1,884,993 rows of transcripts the cache calls complete | 99,670 | 12,604 | not needed: a cache holds its UTRs |
| VEP's own RefSeq cache, the same: 21,052 changes, the 203,125 rows of transcripts the cache has no sequence edits for | 22,680 | 5,307 | not needed |

The cache rows are VEP 116.3 run with `--cache --offline --assembly GRCh38` (and `--refseq --exclude_predicted` for RefSeq) over the cache of chr21 and chr22, and compared with `python3 validation/compare_tab.py vep.tab fv.tab --vep-where FLAGS=- --vep-where-not BAM_EDIT=OK`.
A GFF3 has no `cds_start_NF` or `cds_end_NF`, so the transcripts a cache calls incomplete at an end are left out, and so are the RefSeq transcripts it edits (`_rna_edit`, `BAM_EDIT=OK`), whose sequence no GFF3 holds.
`validation/classify_hgvsp.py` groups the cells that differ by the shape of what each tool wrote, which is how the causes below were told apart.
The corrected column is VEP with `BaseTranscriptVariation::_generic_utr` reading the UTR from `spliced_seq`, as the core API's own accessors do.

Each cell that is left has a cause, below.

**What it was.**
fastVEP wrote the notation from rules read off Ensembl's output: a chain of eight cases (substitution, in-frame indel with an anchor that had to be corroborated against the peptide, frameshift, stop lost to a substitution, and so on), in each of the two per-variant loops.
They agreed with VEP on nearly every row and failed at the places a rule has no case for.
Ensembl's notation is not derived from the protein the change makes.
It is `hgvs_protein`, a sequence of string operations on the peptides of the codon window the change edits, which is the `Amino_acids` pair, of the change *after* the 3'-rule has moved it, and what it writes at a start codon, a stop codon or the end of a protein is the sum of those operations and nothing else:

| Change | Ensembl VEP 116.3 | fastVEP before |
|---|---|---|
| ENST00000427445 `21:10538724 A>AA`, a frameshift that loses the initiator | `p.Met1?` | `p.Met1AsnfsTer3` |
| ENST00001003897 `21:25945916 A>AT`, an insertion between the last codon and the stop | `p.Ter583MetextTer43` | `p.Ter583MetfsTer44` |
| ENST00001003897 `21:25945913 CTCA>C`, the stop codon deleted in frame | `p.Ter583delextTer23` | `p.Ter583del` |
| ENST00000346798 `21:25975228 G>GGGT`, an insertion the 3'-rule carries over an exon boundary | `p.Gln433_His434insThr` | none |
| ENST00000286835 `21:31694203 C>CC`, a frameshift the 3'-rule carries out of the coding sequence | none | `p.Ser442ValfsTer8` |

**What it is.**
`fastvep_hgvs::protein` is a transcription of `hgvs_protein` and the routines it calls (release 116), one Rust function to each Perl sub and in the same order, so that a difference between the two can be found by reading them side by side.
`fastvep_annotate::hgvsp_for_allele`, which both per-variant loops call, works out what it is given: the 3'-shift over the genome (the one `HGVSc` uses, 1,000 bases either side), the change's protein position *where it was moved to*, the window of codons there, and the translation of the edited CDS with its 3' UTR when a frameshift or a lost terminator has to be counted.
Two details of Ensembl's order matter.
`$pre->{coding}` is asked of the change where it was, and `translation_start` of the change where it was moved to, so a change that the rule carries out of the CDS has no notation and one it carries onto the initiator has one.
And `frameshift` is the one predicate evaluated afresh at that point, on the CDS span of the moved change, where `start_lost` and `stop_lost` are cached from the first evaluation: a single base inserted in front of the initiator is not a `frameshift_variant`, and once the rule has carried it onto the first codon it is `p.Met1AsnfsTer?`.

It is checked three ways, none of them a reading of the source alone.
VEP 116.3's own output (above); an in-tree property test that applies each of 6,000 random in-frame descriptions back to the reference protein and requires the protein the variant makes; and the unit tests of each branch.

**What is left, and who is right.**

| Cells | What | Who is right |
|---|---|---|
| 1,492 of the ClinVar sample's 1,818, 124 of the edge set's 237 | A count, a residue or `Ter?` of a frameshift or an extension, on a reverse-strand transcript | VEP's `--gff` reader: [divergence 5](#5-a-frameshifts-terminator-distance-counts-to-the-stop-the-frameshift-creates). Not in VEP's caches, which agree with fastVEP on these |
| 162 of the ClinVar sample, 2 of the edge set | A synonymous change across two codons | fastVEP: [divergence 3](#3-a-synonymous-multi-residue-window-names-the-whole-span) |
| 111 of the ClinVar sample, 12,356 of VEP's Ensembl cache run, 5,169 of its RefSeq cache run, 61 of the edge set | An unresolved residue read as a terminator, and a change at a stop codon that moves the frame | fastVEP: [divergence 13](#13-a-change-at-a-stop-codon-that-moves-the-frame-and-an-incomplete-or-ambiguous-codon-an-unresolved-residue-is-not-a-terminator) |
| 49 of the edge set, 213 of the Ensembl cache run, 24 of the RefSeq | An insertion on the reverse strand that the 3'-rule carries further than it is long | fastVEP: [divergence 12](#12-an-insertion-the-3-rule-carries-further-than-it-is-long-on-the-reverse-strand-is-described-as-another-insertion) |
| 31 of the ClinVar sample, 3 of the Ensembl cache run, 1 of the edge set | A deletion or an insertion in a repeat at the end of a protein, that Ensembl does not shift as far as the end | fastVEP: [divergence 9](#9-the-smaller-ones), NT5C2 |
| 9 of the ClinVar sample | `start_lost` in VEP and not here, on a transcript annotated 5' incomplete | open: [divergence 9](#9-the-smaller-ones) |
| 2 of the ClinVar sample | A frameshift on chrM, where VEP translates the edited sequence with the standard code (`_get_alternate_cds` calls `translate()` without a table) and fastVEP with the mitochondrial one | fastVEP |
| 32 of the Ensembl cache run, 114 of the RefSeq cache run | A transcript whose initiator is not `ATG`: `p.Leu1?` here where VEP writes `p.Met1?` | was VEP's reading of the transcript's own annotation, and is read now: [against VEP's own cache](#against-veps-own-cache) |

**What Ensembl writes that is not HGVS, and fastVEP writes too.**
The notation of a change that loses the initiator is `p.Met1?`, and Ensembl's own cache run has it on 71,797 of the 96,350 `start_lost` rows of the complete transcripts.
On 11,308 it is a window of more than one residue, `p.MetAsn1_?2`, two residue names sharing one position, which is the form `p.LeuIle2672=` was.
On 2,028 (2.1 %) the range is inverted and names a residue that is not the first: for an insertion in front of the initiator the window has no residue before it, and Ensembl's `substr` reads one from the end of the protein, so `p.Asp1_?0` names the last residue of the protein as if it were residue 1.
fastVEP writes all three as Ensembl does, because they are what Ensembl writes and a migration compares them; each is `p.Met1?` as HGVS has it, and which of them fastVEP should write is a decision for the maintainers and not a side effect of this one.

### An inversion inside padding

**The world: 15 `HGVSc` cells.**

`AAC>AGT` replaces `AC` with `GT`, which is its reverse complement.
HGVS asks for `inv`, and fastVEP writes `c.198_199inv`.
Ensembl writes `c.198_199delinsGT`: `hgvs_variant_notation` tests for the inversion on the alleles as it holds them, `AAC` and `AGT`, which are not reverse complements, and `_clip_alleles` afterwards clips the shared `A` and leaves the type as `delins`.
Both describe the change.
A replacement that is wholly an inversion (`CA>TG`) is `inv` in both tools, as the 121 of the ClinVar sample are.

### Coding consequence terms

**The edge set: 447 `Consequence`, 176 `IMPACT`, 12 `Amino_acids` and 12 `Codons` cells.
After: 9, 3, none and none, and those are VEP's.
The generated world: 10 `Consequence` and 5 `IMPACT` cells, none after, and its 17 default columns agree on all 122,621 cells.
The lone-indel set: 4 and 3, none after.
The ClinVar sample: 275 cells of 81,543,670, 43 after (33 rows), and 20 once a frameshift intron is read as Ensembl reads it (17 rows, below, and [Frameshift introns](#frameshift-introns)), and 102 against a `--gff` run once the CDS flags of VEP's cache are read, which a GFF3 does not have ([Against VEP's own cache](#against-veps-own-cache)).
VEP's own RefSeq cache, 150 transcripts of chr21 and chr22 (21,052 changes at their start codons, their stop codons and one internal exon boundary each; 203,125 rows on the transcripts the cache has no sequence edits for): 14,109 rows differed in a term, an impact, a position or a residue, and 48 do, all of them three transcripts with a non-AUG start codon (below).**

The consequence terms were a port of release 115's `VariationEffect.pm`, and release 116 changed it.
The diff of the two files is 155 lines, and every predicate it touches is about the terminator or the frame (`inframe_insertion`, `inframe_deletion`, `stop_gained`, `stop_lost`, `stop_retained`, `ref_eq_alt_sequence`, `frameshift`, and two new helpers that count an insertion's length into the position):

| Predicate | Release 115 | Release 116 |
|---|---|---|
| `ref_eq_alt_sequence`, behind `stop_retained_variant` for a peptide that is whole | held if the replacement began with the residue it replaced and introduced a terminator *anywhere*, or if what spilled past the protein's own length was shorter than three residues, or if both peptides held a terminator at the same index | holds if the protein with the alternate residues put in is the reference protein and what spilled past its length *opens with a terminator*, or if both hold a terminator at the same index |
| `stop_retained` | `ref_eq_alt_sequence` for any alternate peptide that is not empty, and for none `_overlaps_stop_codon` and `_ins_del_stop_altered` | `ref_eq_alt_sequence` only for one that is not empty and has no `X`, and otherwise the `_cil` pair: the same two tests with the length of an inserted sequence counted into the position |
| `stop_lost` | the peptides, or the edited sequence when there is no peptide | declines on a partial codon; asks the edited sequence whenever the alternate peptide holds an `X` |
| `stop_gained` | declined when `stop_retained` held | also declines when `stop_lost` holds |
| `frameshift` | the length arithmetic on the CDS span | declines when the reference peptide begins with the terminator: a change that starts on it does not move the frame |
| `inframe_insertion`, `inframe_deletion` | | `inframe_insertion` declines when both peptides are `*`, and `inframe_deletion` when the reference peptide is `*` |

The port had also taken two decisions of its own that the table does not cover: a change that reached past an end of the CDS was handled by a separate path, and a transcript whose CDS began part-way through a codon had no initiator at all.
Both are gone: one set of tests, `crates/fastvep-consequence/src/edits.rs`, is a transcription of the five predicates that edit the transcript (`_ins_del_start_altered`, `_inv_start_altered`, `_snp_start_altered`, `_ins_del_stop_altered` and its `_cil` twin) and the three `overlaps` tests that guard them, and `terms_for_window` is the predicate graph of release 116 in the order its predicates call one another.
The order matters and is not a severity order: `stop_retained` gates both `frameshift` and `stop_gained`, `start_lost` gates `inframe_insertion` and `protein_altering_variant`, and `inframe_insertion` is evaluated by `start_lost` with `start_lost` still at its default of 0, so the pair is decided in two steps.

Six more things that the probes found, each checked against VEP 116.3:

- **Every term has an `include`** in `Config.pm` that decides whether its predicate is asked at all: `frameshift_variant` only of a change that moves the length (`snp => 0`), `inframe_insertion` and `inframe_deletion` only of one that adds or removes bases, `missense_variant` only of one that does neither.
  Predicates that call each other do not go through it.
  Without the gate a frameshift whose alternate peptide is as long as the reference (`K/X`) is a `missense_variant`.
- **A change that crosses an intron is measured in the CDS.**
  `codon` takes the CDS span (`cds_end - cds_start + 1`) as the length of what the whole alternate allele replaces, so a deletion of 17 coding bases and the 600 intron bases between them is a frameshift with the six codons around it named (`LLLPSI/X`, `p.Leu98HisfsTer3`).
  It was `coding_sequence_variant` with no residues.
- **`_before_coding` and `_after_coding` are plain `overlap` calls against an interval that can hold no base.**
  For a transcript whose coding region reaches its own end the interval is `[end + 1, end]`, and a change that spans that end is a `3_prime_UTR_variant`: `stop_retained_variant,3_prime_UTR_variant` for KRIT1-215's `7:92241100-92241103 ACTT>-` and for eight more on the ClinVar sample.
- **The alternate codon is read from the edited CDS with the 3' UTR appended** (`_get_alternate_cds`), and the reference codon from the CDS alone: a change in a last codon that the CDS leaves short is two bases on one side and three on the other (`X/A`, `gA/gCt`).
- **`non_coding_exon_variant` is a plain `overlap` of the change with an exon**, so an insertion needs both of the bases it sits between inside it.
  One on an exon's first base is a `non_coding_transcript_variant`: 28 rows of the ClinVar sample.
- **An insertion counts as the two bases it sits between** (`min_vf`, `max_vf`) for whether the change reaches the coding region, which is what makes one on the first coding base a coding variant whose predicates then decline for their own reasons.

Two deliberate divergences of Part 1 are gone: the first, because VEP rewrote the clause, and the second, because the faithful rule gives VEP's answer.
Their headings say so.

What is left on the ClinVar sample is 17 rows, and all of them have a cause (it was 33: the other 16 are on two PTEN transcripts with a frameshift intron, [below](#frameshift-introns); and against a `--gff` run there are 44 more, on transcripts the cache flags incomplete, [below](#against-veps-own-cache)):

| Rows | What | Who is right |
|---:|---|---|
| 8 | Reverse-strand transcripts, where a UTR sequence decides `start_retained_variant` or `stop_lost` | VEP's `--gff` reader, which is wrong (below) |
| 9 | A transcript whose CDS begins part-way through a codon (phase 1 or 2, no `cds_start_NF` tag in the GFF3) and a change in its second codon: VEP `frameshift_variant,start_lost`, fastVEP `frameshift_variant` | open: VEP's `--gff` reader never sets `cds_start_NF` (nothing in `GFF.pm` or `BaseGXF.pm` does), so it runs the start-codon predicates on a codon that is not there; a cache dumped from Ensembl's database keeps the attribute (`Dumper/Core.pm` lists it) and `_overlaps_start_codon` returns 0 for a transcript that has it. That the transcripts of these rows have it in the cache was not checked |

**VEP's `--gff` reader builds the UTR of a reverse-strand transcript from the wrong place.**
`_ins_del_start_altered` and the stop tests read the UTR sequence from `BaseTranscriptVariation::_generic_utr`, which transfers the transcript to its feature slice and asks it for `five_prime_utr` or `three_prime_utr`.
For a transcript on the reverse strand built from a GFF3 that returns a sequence of the right length and the wrong content: printed from inside VEP 116.3, the 400 bases of ENST00000000003's 5' UTR end in the reverse complement of genome 10001-10014, which is inside its 3' UTR, while its `spliced_seq` is right.
With `_generic_utr` taking the UTR from `spliced_seq`, as the core API's accessors do, VEP and fastVEP agree on every `Consequence`, `IMPACT`, `Amino_acids` and `Codons` cell of all 36,957 rows of the edge set, and on 164 of 164 real ClinVar rows that changed in this commit.
On VEP's RefSeq cache, which holds its sequences and does not go through that code, the check is direct: of 203,125 rows on unedited transcripts `--gff` mode disagrees with the cache on 897 (in a term, an impact, a position or a residue), every one of them on the reverse strand, and fastVEP agrees with the cache on 849 of the 897.
The 48 others are the three RefSeq transcripts of MRTFA, whose start codon is `GUG`: the cache reads the first residue as `M` (the CDS records say `transl_except=(pos:complement(40552344..40552346),aa:Met)`) and a codon table says `V`, which is what `--gff` mode and fastVEP write.
The committed edge fixture lists the one row of it that its sample holds, and asserts that it still differs.
The same accessor supplies the 3' UTR that a frameshift's or an extension's count is made against, and that is the whole of what divergences 5 to 8 turned out to be.

On the 134,732-record ClinVar sample 161 of the 164 rows that changed in this commit agree with stock VEP 116.3, and the three that do not are those reverse-strand rows; the build before it agreed with 2 of the 164.
The 164 rows are 71 variants, 0.005 % of the 3,137,188 rows, and the biggest groups of them are `inframe_insertion,stop_retained_variant` to `frameshift_variant` (24), to `inframe_insertion` (18) and to `stop_retained_variant` (10), `stop_gained,inframe_insertion` to `inframe_insertion,stop_retained_variant` (15), and `non_coding_transcript_exon_variant` to `non_coding_transcript_variant` (28).

### Frameshift introns

Ensembl calls an intron a frameshift intron when its last base is at most 12 bases past its first (`abs($intron_end - $intron_start) <= 12` in `BaseTranscriptVariation.pm`), which is an intron of 13 bases or fewer.
The comment above that test says 12 bases, and the test is what runs; this document said 12 as well, and so did fastVEP's code (`< 12`), which the 13-base intron of the world below showed to be wrong.
Such an intron is there so that a transcript agrees with an assembly that has an extra or a missing base, and it changes what Ensembl says about the whole transcript.
Four things follow in release 116, each read from the Perl and each checked by running VEP 116.3:

- **A change inside the intron earns no splice term and no `intron_variant`.**
  `_intron_effects` sets `within_frameshift_intron` and goes on to the next intron, in each of its two loops.
  fastVEP did this already, for the introns of 12 bases or fewer.
- **Inside the coding region it is a `coding_sequence_variant`.**
  `within_cds` has a second way to hold, for a change in a frameshift intron that overlaps the coding region, and `coding_unknown` follows from it when no peptide is made.
  fastVEP wrote the catch-all, `coding_transcript_variant`, which in Ensembl is a change that covers a whole transcript.
- **Every exon is 12 bases wider on each side for the "does it touch an exon" test** (`_overlapped_exons`, `$stretch`), in that transcript and in no other.
  `splice_polypyrimidine_tract_variant` has `exon => 0`, so it is gone from the 3 to 12 bases before every exon of the transcript, not only the one beside the short intron, and stays at 13 to 17.
  The test is asked of the change's two flanking bases sorted (`min_vf`, `max_vf`), so an insertion touches an exon if either of the bases it sits between does.
- **The same test decides whether a change is `coding` for `HGVSp`** (`_bvfo_preds`: the `exon` flag decides whether the mapper is asked, and the mapper answers a change in an intron with a gap, which is a list that is not empty).
  `hgvs_protein` then moves the change with the 3'-rule, and when that carries it into the exon Ensembl writes the protein change: `p.Glu56SerfsTer136` for the deletion of one base in the two-base intron of the SHANK3 transcript ENST00000673971 (`22:50674581`, `c.166del`).
  fastVEP wrote none.

One more difference came out of the same rows, and it is not about introns.
Ensembl writes the reference codons in uppercase over `length($self->feature_seq)` bases from the codon position (`display_codon`), and for the reference allele that is the length of the change in the genome: a deletion of five bases over a one-base intron removes four bases of the CDS, and `CACCGc` has five uppercase letters where fastVEP wrote four (`CACCgc`).
The codons are cut with the CDS span, as before; only the case follows the genome.

**Measured against VEP 116.3.**
A world of seven transcripts with introns of 1 to 60 bases, on both strands, coding and not, two of them controls whose introns are all 14 bases or more (`crates/fastvep-cli/fixtures/vep_tab_fsintron`): a substitution, deletions of 1, 3, 5 and 15 bases and insertions of 1 and 3 at every position within 18 bases of an exon edge and in every intron of 40 bases or fewer.
Over all 11,515 records (31,024 rows) 977 `Consequence`, 574 `IMPACT`, 44 `Codons` and 8 `HGVSp` cells differed from VEP, and none does now.
The committed fixture is a sample of it, 2,879 records and 7,756 rows (every substitution, every eighth of the rest), on which 234, 136, 6 and 2 cells differed and none does.
VEP's own caches, over a substitution, deletions of 1 and 3 bases and insertions of 1 and 3 at every position within 25 bases of every exon edge of the transcripts of chr21 and chr22 that have such an intron, and in each of those introns:

| Cache | Transcripts | Rows | Before | After |
|---|---|---:|---|---|
| Ensembl 116 | ENST00000673971 and ENST00000692848 (the two SHANK3 isoforms, each with a 2-base intron) and ENST00000624847 (a lncRNA with two) | 124,027 | `Consequence` 2,342, `IMPACT` 946, `HGVSp` 57 | none, and `HGVSp` 2 |
| RefSeq 116 | NM_001396027.1 (two 1-base introns) and NM_001372044.2 (two 2-base introns) | 15,084 | `Consequence` 1,208, `IMPACT` 484, `Codons` 2, `HGVSp` 39 | none, and `HGVSp` 5 |

The `HGVSp` cells that are left are the stop-codon family of [divergence 13](#13-a-change-at-a-stop-codon-that-moves-the-frame-and-an-incomplete-or-ambiguous-codon-an-unresolved-residue-is-not-a-terminator): on the RefSeq transcript `p.Ter232%3D` against `p.Ter232delextTer?` and `p.Ter232delinsTrpTer` against `p.Ter232TrpextTer?`, and on ENST00000673995 of the Ensembl probe, which has no frameshift intron, `p.Ter1%3D` against `p.Xaa1%3D` and `p.Ter1_Pro2insThr` against `p.Xaa1_Pro2insThr`.
Both probes are made with `validation/gen_short_intron_probe.py`, and the numbers are `validation/compare_tab.py`'s, restricted to those columns.
Nothing outside a transcript with such an intron changes: the tab files of HG002 (186,055 rows) and of the edge set (36,958 rows) are byte for byte what they were.
On the 1-in-5 sample of ClinVar 2-star+ (3,136,295 rows) 16 rows change, 13 on ENST00000693560 and 3 on ENST00000850883, two PTEN transcripts with a one-base intron (16 `Consequence` and 7 `IMPACT` cells), and all 16 are VEP's now.

### `--total_length`, `--no_escape` and `--force_overwrite`

The reporter's command line carries these three, and `fastvep annotate` refused each of them with an error, so a migrated command failed at once.
They are accepted now, in VEP's spelling (`--total_length`) and in the hyphenated one the rest of the CLI uses (`--total-length`), and each does what VEP 116.3 does.
That was read from `OutputFactory.pm` and `OutputFactory/VCF.pm` and measured by running VEP on the same input.

**`--total_length`** writes the length each position is counted in after the cell: `451-?/1403` for the cDNA, `51/654` for the coding sequence with its stop codon, `17/217` for the protein without it.
VEP writes it whether or not the cell has a value.
`OutputFactory.pm` appends it inside the two branches that test the cached `_bvfo_preds` flags, `exon` for the cDNA and `coding` for the other two, and `format_coords` of a span with no coordinate at either end is `-`, so the cell is `-/1403`.
That is 76 of the 4,753 rows of the edge fixture for the cDNA and 62 for the CDS and the protein: a change that covers an exon from intron to intron, or from outside the transcript, and an insertion immediately before the initiator or after the terminator.
In a transcript with a frameshift intron it is also every change within 12 bases of an exon (301 rows for the cDNA and 212 for the CDS in the frameshift-intron world below).
The two flags are the ones `HGVSp` asks for, `Transcript::pre_exon` and `Transcript::pre_coding`, defined once.
VEP's VCF writer blanks a cell that is a bare `-` and not `-/1403`, so a CSQ field has the dash over the total and an empty field where there is no total.
The structural variants fastVEP writes no position for are given no total either.

- **JSON.** VEP's JSON writer makes `cdna_start` and `cdna_end` by splitting the cell at its first `-`, and with a total that gives the wrong answer: 1,525 of the 1,633 ranges of the edge fixture have their start for an end, and the 90 single positions are the string `"5/1403"`, for both.
  fastVEP does not reproduce that and writes no totals in JSON, and says so on stderr.
  This is a decision for the maintainers; the conservative reading is that a wrong end is not worth matching.
- **Without `--fasta`** the lengths are worked out from the coding span and its phase, which is the length of the translateable sequence, and the protein length assumes the coding sequence ends in a stop codon.
  For a transcript whose coding sequence does not (`cds_end_NF`) that is one residue short.
  This was read from the code and not measured: the runs below have a FASTA.

**`--no_escape`** leaves the `=` of a synonymous `HGVSp` alone, in the tab file and in the VCF (`p.Val3=`, not `p.Val3%3D`), and changes the space of a gene symbol in the VCF.
VEP puts `%20` for the space of a `SYMBOL` when it builds the value for VCF output (`$hash->{SYMBOL} =~ s/\s/\%20/g if ... eq 'vcf' && !$self->{no_escape}`), and its VCF writer then turns whatever whitespace is left into `_`: `GENE%20A` by default and `GENE_A` with the flag, measured with a gene named `GENE A`.
fastVEP wrote `_` always, which is VEP's `--no_escape` and not what it does by default, so the default changed too.
The tab file never escaped a symbol, in either tool.
VEP's own JSON writes a raw `=` whatever the flag, as fastVEP's does.

- **GFF3 percent-escapes** came out of testing it and are not part of this change: VEP's GFF reader leaves `Name=GENE%20B` as `GENE%20B`, and fastVEP's decodes it (`GENE B`, which the VCF writes as `GENE%20B`, or `GENE_B` with the flag), as the GFF3 specification asks.
  A symbol with an escaped comma is `GENE%2CC` in VEP and `GENE&C` in a fastVEP CSQ.

**`--force_overwrite`** is accepted and does nothing: VEP refuses to write over an existing output file without it, and `fastvep annotate` always writes over its output.

Measured against VEP 116.3, cell by cell, with the flags on in both tools:

| Input | Rows | Position cells (cDNA, CDS, protein) | Other columns |
|---|---:|---|---|
| The edge fixture, tab (`vep_tab_edges`, 2,271 changes) | 4,753 | **all agree** | the 36 cells listed there |
| The same, VCF `CSQ` | 4,753 entries | **all 14,259 agree** | |
| The frameshift-intron world (`vep_tab_fsintron`, 2,879 changes), where 301 rows are `-/N` for the cDNA and 212 for the CDS | 7,756 | **all agree** | all 131,852 cells agree |
| A RefSeq-style world with a symbol that has a space and synonymous changes (`vep_flags_compat`, 760 records) | 1,418 | all agree | all 42,540 cells agree, with each flag and without |
| HG002 sample, Ensembl 115 GFF3, `--gff` (20,241 records) | 186,055 | all agree | the 3,834 cells of Part 3 and no others |
| VEP's Ensembl cache, chr21 and chr22, 45,916 changes, complete transcripts | 1,884,993 | all agree | the `HGVSp`, `HGVSc`, `Amino_acids` and `Consequence` cells of the `HGVSp` section, and the `SYMBOL_SOURCE` and `HGNC_ID` a GFF3 does not carry |
| VEP's RefSeq cache, the same, unedited transcripts | 203,125 | all agree | those of the `HGVSp` section, and the `SYMBOL_SOURCE`, `BIOTYPE` and `CANONICAL` of a RefSeq GFF3 |

The committed goldens are `crates/fastvep-cli/fixtures/vep_flags_compat/golden_total_length.tsv.gz` and `golden_no_escape.tsv.gz`, with a `regen.sh` that runs VEP in Docker, and `vep_flags_compat_golden.rs` holds fastVEP against them: both fail if the dash over a total is dropped.

### Alleles that are not plain bases

**HG002: the 6 rows written for a `*` allele.
The allele set (`fixtures/vep_tab_alleles`): 100 rows only fastVEP wrote and 434 of the 10,800 cells of the 360 rows both wrote, before; none, and 10 (divergence 11), now.**

Ensembl translates an allele only when it, and the reference, are runs of `A`, `C`, `G` and `T` (`seq_is_unambiguous_dna`, `^[ACGT-]+$`).
Anything else has no peptide there, so no `Amino_acids`, no `HGVSp`, no `HGVSc`, and no term that is about the protein.
Each kind was run through VEP 116.3 at nine kinds of place (both strands of a coding exon, a UTR, an intron, a non-coding exon, and an intergenic base; 171 records) and fastVEP now agrees on every row and every cell but the 10 of divergence 11, except where the last column says otherwise:

| ALT or REF | Ensembl VEP 116.3 | fastVEP before | Now |
|---|---|---|---|
| `*`, the spanning deletion | no row for it; a record whose only ALT is `*` has no row at all, in a VCF no `CSQ`, in JSON no consequences | a row for it on every transcript and an intergenic one, `Consequence` empty, `MODIFIER` | the same as VEP |
| `.` | the record is not annotated, and a VCF output does not have the line | the rows of a `*` allele | no row; the VCF keeps the line, unannotated |
| an ambiguity code in the ALT (`N`, `R`, `Y`) | in the coding sequence `coding_sequence_variant` MODIFIER, `Amino_acids`, `HGVSc` and `HGVSp` `-`, or `frameshift_variant` HIGH if its length moves the frame; elsewhere the term of the place, with no `HGVSc` | `missense_variant,coding_sequence_variant` MODERATE, `K/X`, `c.91A>N`, `p.Lys31Xaa`; `p.Gln28XaafsTer9` for an insertion | the same as VEP |
| an ambiguity code in the REF | `coding_sequence_variant`, no `Amino_acids` or `HGVSp`, but `HGVSc` against the genome's base: `c.101T>C` for `N>C` at a `T` | `missense_variant` `V/A`, `c.101N>C` | the same as VEP |
| `<NON_REF>` | the whole record is skipped, with a warning, even beside a real ALT | beside a real ALT a `frameshift_variant` row for it (it was read as a symbolic deletion) | **no row for it, and the real ALT beside it annotated**: VEP loses a real variant of a gVCF record, and this does not |
| `<*>` | a `coding_sequence_variant` row for it, which says nothing | beside a real ALT a `frameshift_variant` row | **no row for it**, for the same reason as `<NON_REF>` |

An ambiguity code is complemented on the reverse strand as Ensembl's `reverse_comp` has it (`R` is `Y`), which the `Codons` cell shows; `Allele::is_plain_bases` holds the rule.
fastVEP also writes no `HGVSg` for them, as Ensembl does not (`next if $allele =~ m/[^ACGT\-]/ig` in `hgvs_genomic`), and used to write `chr10:g.1400?`, which is not HGVS.

Two things are left as they were and are not claimed to be VEP's: an ambiguous allele that overlaps the stop codon (VEP says `stop_retained_variant` for a length change there, because its test of the edited codon declines to run on an unknown base), and the symbolic alleles of a structural variant, which have their own predictor.

### Lone indels that are not minimal

**The world with each ALT its own record: 197 rows of 7,213 had another `Allele`, `Location`, `REF_ALLELE` and `VARIANT_CLASS`, in 93 records, and 8 `Consequence` cells differed on them (4 of them `IMPACT`).
The lone-indel set: 14,526 of 95,280 cells differed.
After: none of those, and 59 cells, 52 once the coding terms were ported and 9 once the HGVS notation was: the `HGVSp` of eight changes at the stop codon of a reverse-strand transcript, which is VEP's `--gff` reader (divergence 5), and of one that moves the frame there (divergence 13).**

Ensembl minimises a record with one ALT that is an indel from both ends (`Parser::post_process_vfs`, `minimise_alleles`, `trim_sequences`): the bases the alleles share at the front, then at the back, until one of them is empty.
fastVEP stripped only the first base they share.
The two are the same for a record that is already minimal and padded by one base, which is every indel of HG002 (0 of its 2,818 are not), and they differed for one that shares bases at the back as well, or two or more at the front.
`CAGT>CT` is the deletion of `AG`, and `CCAT>CCT` the deletion of `A`, not a deletion of `AGT` and of `CAT`; `Allele`, `Location`, `REF_ALLELE`, `VARIANT_CLASS` and `DISTANCE` follow from it, and so does what the consequence step sees: for `1306 AG>ACTATAATG` Ensembl writes `Allele` `CTATAAT` at `1306-1307` as an `insertion`, and fastVEP wrote `CTATAATG` at `1307` as an `indel`, with `splice_donor_variant` HIGH for an insertion next to a donor that it read as replacing it.
`HGVSc` and `HGVSg` were not affected (both clip), and still are not.

A record with several ALTs is not minimised in release 116, and is not here: in 115 a minimised indel switched on `split_variants` for the whole input (`InputBuffer::next`), which split such a record into one variant per ALT, and in 116 those lines are commented out, so it is annotated whole.
`fixtures/vep_tab_lone_indels` is 1,500 records written the way a caller that does not trim writes them, a prefix the alleles share (the anchor and up to two more bases), the bases that change and a suffix they share (none to three), as a deletion, an insertion or a replacement, and what VEP 116.3 wrote for them; `vep_tab_lone_indels_golden.rs` holds fastVEP to it.
All the columns it compares, `Consequence`, `IMPACT`, `HGVSc` and `HGVS_OFFSET` among them, agree on all 3,176 rows, and `HGVSp` on all but the 9 cells `vep_tab_lone_indels_golden.rs` lists, each asserted to still differ.

### A RefSeq GFF3

Everything above is Ensembl's GFF3, and #139 runs RefSeq, so one RefSeq run was measured the same way: Ensembl VEP 116.3 and fastVEP over chr21 and chr22 of RefSeq's GFF3 for GRCh38.p14 (GCF_000001405.40, 134,543 records, the contigs renamed `21` and `22` so that the GRCh38 FASTA above answers for them) and the 527 variants of HG002 on those two chromosomes.
2,717 rows are paired and 135 of their 81,510 cells differ (0.17 %), all of them divergence 11: `HGVSc` 62 and `HGVS_OFFSET` 73, on the 7 of the 527 records that have several ALTs.
`Feature`, `Gene`, `BIOTYPE`, `SYMBOL`, `Consequence`, `IMPACT`, the positions, `Amino_acids` and `Codons` agree on every paired row.
fastVEP writes 8 rows VEP does not: those of the mature microRNA records (`miRNA`, `rna-MIR185`), which hang off a transcript, and which VEP writes nothing for.

That is since fastVEP names a RefSeq record as VEP's GFF reader does (`BaseGXF.pm`, `GFF.pm`), and it did not.
NCBI's GFF3 gives a transcript the ID `rna-NM_000001.3` and the attribute `transcript_id=NM_000001.3`, a gene the ID `gene-TP53` and its transcripts the Dbxref `GeneID:7157`, and states no biotype.
VEP names the transcript by its `transcript_id` (its ID as it is where it has none: an NCBI tRNA record has none, and is `rna-TRG-GCC1-5`), the gene by the `GeneID`, and guesses a biotype the record does not state: `protein_coding` for an mRNA, the `ncrna_class` of an ncRNA, `miRNA` for a transcript of a gene described as a microRNA, `IG_V_gene` for a `V_gene_segment`, and the `gbkey` for the rest (`ncRNA` for an lnc_RNA, `misc_RNA`, `rRNA`, `tRNA`).
fastVEP named them by the ID, by the gene's own ID and by the SO term, so that 175 of the 2,717 rows paired, none of them a transcript row but the tRNA's, and with `rna-` taken off by hand `BIOTYPE` differed on 2,542 of 2,716 rows (`mRNA` for `protein_coding` on 2,130) and `Gene` on every transcript row.
Where VEP can name no biotype it drops the transcript with a warning, and fastVEP keeps it, named by its SO term.

Two records can now have one name: RefSeq writes a transcript that aligns in two places twice (`rna-NM_x.1` and `rna-NM_x.1-2`, both `transcript_id=NM_x.1`; 6,111 of the 186,133 ids in the file), and VEP gives both the one name.

The case the reporter runs, VEP's RefSeq *cache* (`--refseq`), is measured next.
`fixtures/vep_refseq_gff` holds a RefSeq-style GFF3 with one transcript of each kind a human one has and what VEP 116.3 wrote for it (138 CSQ entries), and `refseq_names_vep.rs` holds fastVEP to it.

### Against VEP's own cache

A cache is not a GFF3.
Ensembl builds it from its database, so what it holds of a transcript is the database's: the symbol source, the canonical flag, the biotype, whether the CDS is complete at either end, and whether the first codon is read as Met.
A GFF3 carries some of that differently and some not at all, and the reporter of #139 runs `vep --cache --offline --refseq`.
Chromosomes 21 and 22 of VEP 116.3's two human caches (`homo_sapiens/116_GRCh38` and `homo_sapiens_refseq/116_GRCh38`, about 700 MB each, taken from the `indexed_vep_cache` archives) were run through VEP 116.3 itself with `--cache --offline` (and `--refseq --exclude_predicted` for RefSeq's) and the columns above, and through fastVEP over the matching GFF3: Ensembl 116's, and RefSeq's GRCh38.p14 with the contigs renamed and `NM_` and `NR_` records only, which is what `--exclude_predicted` leaves.
Rows are paired as above.
`validation/run_cache_oracle.sh` makes both files and `validation/compare_cache_run.py` compares them, by cell and by transcript.
A RefSeq cache also holds transcripts that VEP edits with an alignment (`_rna_edit`, `BAM_EDIT=OK`: 59 of the edge set's 1,061 and 194 and 57 of the two samples'), which no GFF3 can reproduce, and they are left out of both files.

Four inputs: the edge set of each cache (every shape of change within ten bases of every exon, CDS and transcript boundary: 45,916 variants on the Ensembl cache and 21,052 on the RefSeq one), every 20th variant of HG002 on the two chromosomes (5,275) and the ClinVar 2-star+ sample on them (4,712).
Cells that differ, by column, with the previous commit's fastVEP (`3082129`) first and this one's second; the Ensembl runs are with `--exclude-readthrough --exclude-artifact`, below.

| Input | Pairs | Rows one tool alone has | Cells that differ |
|---|---:|---|---|
| Ensembl cache, edge set | 2,026,951 | fastVEP 18,733 -> 0 | `SYMBOL_SOURCE` 1,994,540 -> 0, `HGNC_ID` 1,994,048 -> 3,264, `FLAGS` 141,958 -> 642, `Consequence` 1,694 -> 0, `IMPACT` 566 -> 0, `Amino_acids` 32 -> 0, `HGVSc` 417 -> 417, `HGVSp` 14,766 -> 13,842 |
| RefSeq cache, edge set | 203,125 | fastVEP 436 -> 0 | `SYMBOL_SOURCE` 203,125 -> 0, `BIOTYPE` 1,677 -> 0, `CANONICAL` 26,757 -> 150, `Consequence` 6 -> 0, `IMPACT` 3 -> 0, `Amino_acids` 48 -> 0, `HGVSp` 5,307 -> 5,193 |
| Ensembl cache, HG002 | 81,309 | fastVEP 559 -> 0, VEP 1 -> 0 | `FLAGS` 2,232 -> 35, `SYMBOL_SOURCE` 74,363 -> 0, `HGNC_ID` 74,311 -> 533, `HGVSc` 338 -> 338, `HGVS_OFFSET` 437 -> 437 |
| RefSeq cache, HG002 | 15,426 | fastVEP 66 -> 2, VEP 1 -> 1 | `SYMBOL_SOURCE` 13,222 -> 0, `BIOTYPE` 1,643 -> 0, `CANONICAL` 3,269 -> 15, `HGVSc` 122 -> 122, `HGVS_OFFSET` 146 -> 146 |
| Ensembl cache, ClinVar | 168,813 | fastVEP 872 -> 0 | `Consequence` 3 -> 1, `IMPACT` 2 -> 0, `FLAGS` 9,586 -> 384, `SYMBOL_SOURCE` 167,127 -> 0, `HGNC_ID` 167,117 -> 34, `HGVSc` 10 -> 10, `HGVSp` 4 -> 4 |
| RefSeq cache, ClinVar | 26,383 | fastVEP 240 -> 0 | `SYMBOL_SOURCE` 26,383 -> 0, `BIOTYPE` 471 -> 0, `CANONICAL` 5,903 -> 0: none differ |

What was different, and what a cache does:

- **`SYMBOL_SOURCE` and `HGNC_ID`.**
  A RefSeq gene with one NCBI gene id is `EntrezGene` on every row of the cache and has no `HGNC_ID`, though its Dbxref has one (`prefetch_gene_ids`).
  An Ensembl gene says where its symbol is from at the end of its description (`[Source:HGNC Symbol;Acc:HGNC:12023]`): of the 2,031 genes of the Ensembl cache run, the 1,153 that say `HGNC Symbol` are `HGNC` with that accession, the 27 that say `RFAM` are `RFAM`, the 13 that say `NCBI gene (formerly Entrezgene)` are `EntrezGene`, and the 838 that say nothing have neither.
- **`BIOTYPE`.**
  The cache gives a RefSeq transcript the `gene_biotype` of its gene, except `misc_RNA` for a non-coding transcript of a protein-coding gene: 503 of 503 `lnc_RNA` records of lncRNA genes are `lncRNA` (VEP's `--gff` reader and fastVEP before wrote their `gbkey`, `ncRNA`), the 67 `transcript` records of transcribed pseudogenes `transcribed_pseudogene`, and an `snRNA` or `snoRNA` its gene's.
  Only the pairs seen on chr21 and chr22 are changed; the biotypes RefSeq has that are not on them (`tRNA`, `scRNA`, `antisense_RNA`, `vault_RNA`) keep the name the `--gff` reader gives.
- **`CANONICAL`, `MANE` and `MANE_SELECT`.**
  RefSeq's GFF3 has no `Ensembl_canonical` tag, and writes `tag=MANE Select` with a space, where Ensembl writes `MANE_Select`; the record names Ensembl's transcript in its Dbxref, which is what the RefSeq cache writes as `MANE_SELECT`.
  The cache marks the MANE Select transcript of 640 genes of 640, and otherwise the transcript that spans the most of the genome, the first in the file among equals: 451 of the 455 genes that have no MANE Select transcript (366 have one transcript, and of the other 89 the cache's is the longest-spanning in 85, and the first of those in the file in all 85).
  Neither the length of a transcript's exons nor of its CDS does better than the span (63 and 52 of 89).
  The rule needs all of a gene's transcripts.
  A whole file has them.
  A tabix read has the ones that overlap its input, and the longest of those is not the gene's, so a gene across the edge of a region is read again, whole, for the flag (and for it alone: its other transcripts are not added to the run).
  A plain and a bgzipped copy of one sorted file then mark the same transcripts: on the sorted GFF3 of chromosomes 21 and 22 the two runs differed in `CANONICAL` on 4,266 of 212,240 rows before, and on none after.
- **The mature miRNA of a precursor.**
  A record that hangs off a transcript (`rna-MIR185`, `Parent=rna-NR_029706.1`) is where in the precursor the miRNA is, and neither the cache nor the `--gff` reader has a transcript of it: 436 rows of the edge set, and 66 and 240 of the samples, were fastVEP's alone.
- **`FLAGS`, and the first codon.**
  Ensembl's GFF3 export has no `cds_start_NF` or `cds_end_NF`, and the cache has them on 5 % of its coding transcripts: a change at the last codon of an incomplete transcript was `stop_lost`, a frameshift at its first was `start_lost` (1,694 `Consequence` and 566 `IMPACT` cells of the edge set).
  They are read off the sequence once a run has the FASTA, and agree with the cache on 12,859 (start) and 12,858 (end) of 12,867 coding transcripts: `cds_end_NF` when the CDS does not end in a stop codon in frame (all 621 the cache flags, and 9 it does not), `cds_start_NF` when the first CDS exon starts part-way through a codon (all 229) or its first codon is not `ATG` and there is no 5' UTR above it (a transcript that begins with another codon and has a UTR is complete in the cache, 18 of 18; one without is not, 128 of 130).
  A GFF3 that carries the flags (GENCODE's) keeps them.
  The first codon: where the annotation reads a start that is not `ATG` as Met, residue 1 is `M`: Ensembl's `initial_met` edit, applied to the reference peptide of the codon and to nothing else, which the cache has on all 14 complete transcripts that start with `CTG`, `GTG` or `ACG` and none of the six that start with something else, and RefSeq's `transl_except=(pos:...,aa:Met)` (MRTFA's `GUG`).
  `Transcript::translate` also writes a first codon that the table lists as a start (`TTG`, `CTG` and `ATG`, read from BioPerl in VEP's image) as `M`, and so does fastVEP, which is the protein a frameshift's `p.Met1?` is read against.
  `Amino_acids` reads `M/L` for a `CTG` start and `M/V` for MRTFA's, `p.Met1?` replaces `p.Leu1?` and `p.Val1?` (32 and 48 `Amino_acids` cells, 32 and 114 `HGVSp` cells).
- **Transcripts the cache does not hold.**
  VEP's Ensembl cache is built without readthrough transcripts and without those of biotype `artifact` (`AnnotationSource::Database::Transcript`): 109 of the 23,382 transcripts of chr21 and chr22 in Ensembl's GFF3, and the only ones of them the cache lacks.
  `--exclude-readthrough`, `--exclude-artifact` and `--exclude-predicted` (VEP's `--exclude_predicted`: `not stable_id match ^X._`) leave them out of a run on request, and with them the Ensembl runs have the cache's transcript set row for row.
  None is on by default.

Every cell that is left has a cause:

- **`HGVSc` and `HGVS_OFFSET`**: an allele of a record with several ALTs, which is [divergence 11](#11-an-allele-of-a-record-with-several-alts-is-described-as-it-would-be-alone) (HG002: 121 and 146 on the RefSeq cache, 328 and 437 on the Ensembl one), and the `*` coordinate of an intron beside the last coding base, which is [divergence 4](#4-a--coordinate-keeps-the-exonic-anchor-it-is-counted-from): Ensembl's cache has it too, on 417 cells of the Ensembl edge set (31 of 9,013 transcripts), 10 of the Ensembl HG002 run, 1 of the RefSeq one and 10 of the ClinVar sample.
- **`HGVSp`**: the classes of divergences 3, 9, 12 and 13 and the stop-codon family, listed in the `HGVSp` section above, and nothing else: of the Ensembl edge set's 13,842, 12,572 are on transcripts the cache calls complete (12,604 before the initiator, less its 32) and 1,270 on incomplete ones, where `p.Ter157%3D` is VEP's `Xaa` read as a terminator; of the RefSeq edge set's 5,193, the stop-codon family.
  4 of the ClinVar sample's are the same.
- **`CANONICAL`**: 4 of the 455 genes without a MANE Select transcript, and 4 genes the cache has no canonical transcript for at all (coding genes whose transcripts are all `RefSeq Select`), which nothing in the GFF3 explains and which come from the database's own `canonical_transcript_id`: 4 transcripts of HG002 (15 cells) and 1 of the edge set (150 rows).
- **`HGNC_ID`**: VEP copies one gene's HGNC ID to every other gene of the same symbol it loaded together (`merge_features`), so `BCRP3` of ENSG00000215481, whose symbol source is `EntrezGene`, has the `HGNC:1016` of ENSG00000284128.
  That depends on which transcripts are loaded in one call and is not reproduced: 34 transcripts of the edge set, 159 of HG002 and 34 of the ClinVar sample.
- **`FLAGS`**: 4 transcripts of the edge set, 14 of HG002 and 8 of the ClinVar sample, of 12,867 coding ones: the sequence rules above differ from the database's on 9 transcripts' `cds_end_NF` (a CDS of 698 bases that the database calls complete) and on 8 transcripts' `cds_start_NF` (6 the database flags that start with `ATG`, and 2 the rule flags that it does not).
- **`Consequence`**: 1 cell, `mature_miRNA_variant` on ENST00000290239: the cache's `miRNA` attribute holds the mature microRNA's coordinates, which Ensembl's GFF3 export does not carry.
- **Rows one tool alone has**: 2 of fastVEP's, the tRNA `rna-TRG-GCC1-5`, a record with no RefSeq accession, which the cache does not hold, and the intergenic row VEP writes beside it.

A `--gff` run of VEP has none of the first three: a GFF3 has no symbol source, no HGNC ID and no CDS flags, so VEP writes `-` for them.
Against the HG002 sample's `--gff` run (the first table of this part) fastVEP now differs in `SYMBOL_SOURCE` (153,728 of the 5,581,650 cells), `HGNC_ID` (152,539) and `FLAGS` (8,631), and in nothing else that it did not before: it follows the cache, which is what `vep --cache` writes and what the reporter runs.
The CDS flags reach the terms as well.
A transcript the cache flags `cds_start_NF` or `cds_end_NF` has no start or stop codon to lose, and a `--gff` run, which has no flags, calls a change in it `start_lost`, `stop_lost` or `stop_retained_variant`.
Against the ClinVar 2-star+ sample's `--gff` run (3,136,295 rows) fastVEP therefore differs in the `Consequence` of 61 rows (17 before the flags), the `IMPACT` of 41 (3 before) and the `HGVSp` of 1,853 cells (1,818 before): the 44 rows, 38 `IMPACT` cells and 35 `HGVSp` cells that the flags changed are all on a transcript fastVEP flags incomplete, and against VEP's Ensembl cache the same sample differs in 1 `Consequence` cell.
HG002 and the generated worlds are unchanged.

Not reproduced, and not in the columns above:
`APPRIS` (Ensembl's GFF3 has no tag for it: 734 of 3,650 rows of a 150-variant sample carry one in the cache), and the `MANE_SELECT` of an Ensembl transcript, which the Ensembl cache writes as RefSeq's accession (`NM_001754.5`) and the GFF3 does not name (fastVEP writes the transcript's own, as before).
`TSL`, `CCDS`, `ENSP`, `CANONICAL` and `MANE` of an Ensembl transcript agree on that sample.

### The columns of the known variants

`AF` and the population frequencies, `MAX_AF` and `MAX_AF_POPS`, `CLIN_SIG`, `SOMATIC`, `PHENO` and `PUBMED` come from the known variants of VEP's variation cache (`--cache-dir`).
The reporter of #139 runs `--af_1kg --check_existing`.
VEP fills them with two subs of `OutputFactory.pm`, `add_colocated_variant_info` and `add_colocated_frequency_data`, and `crates/fastvep-annotate/src/colocated.rs` is a transcription of them, one function to each.
They were measured against **Ensembl VEP 116.3** the way the rest of this part is, with `crates/fastvep-cache/fixtures/vep_known_columns/regen.sh` and `validation/compare_known_columns.py`.

What these columns say is a function of the input variant and its allele, not of the transcript.
The comparison therefore pairs the two files on (`Uploaded_variation`, `Allele`), whatever transcripts each wrote for the variant, and requires every row of a pair to say the same thing within each file (in no run below does a pair have rows that disagree, in either file).

| Input | Flags | (variant, allele) pairs | Cells compared | Not `-` in VEP | Cells that differ |
|---|---|---:|---:|---:|---:|
| A generated world over a cache with the 42 real columns: 50 hand-written corners of the Perl and 396 generated pairs with random payloads, on transcript rows and on intergenic ones (`fixtures/vep_known_columns`) | eight sets, one golden file each | 673 rows each | 63,262 | 7,489 | 0 |
| GIAB HG002, chr21 and chr22 of the whole callset (105,494 records), against VEP's real release-116 Ensembl cache | `--af_1kg --check_existing` | 106,673 | 1,066,730 | 704,150 | 0 |
| the same | `--af --af_1kg --af_gnomade --af_gnomadg --max_af --pubmed --check_existing` | 106,673 | 3,626,882 | 2,168,576 | 0 (2 before the PubMed fix below) |
| the same records, VEP in `--refseq --exclude_predicted` cache mode over the RefSeq cache | `--af_1kg --check_existing` | 106,673 | 1,066,730 | 704,150 | 0 |
| ClinVar 2-star+ on chr21 and chr22 (4,712 records) | all of them, and the reporter's | 4,704 | 159,936 and 47,040 | 102,363 and 23,925 | 0 and 0 |
| 25,346 records drawn from the rarer shapes of that cache (allele-specific `CLIN_SIG`, records with no alleles, failed and somatic ones, `pubmed`, three or more alleles with gnomAD frequencies, fifteen or more populations, indels with gnomAD frequencies) | all of them | 37,590 | 1,278,060 | 507,658 | 0 |

The RefSeq cache's variation files are byte-identical to the Ensembl cache's (`cmp` on chr21 and chr22), so the reporter's mode reads the same known variants; VEP's `--refseq` run is there to show that nothing in it changes these columns.
`Existing_variation` agrees in all of these too.
The real cache of chr21 and chr22 has no record on the minus strand (every one of its 29.7 million `strand` cells is `.`), so that path is measured on the generated world only.

Reproduce the world with `bash crates/fastvep-cache/fixtures/vep_known_columns/regen.sh` (it needs Docker, `bgzip` and `tabix`), and a real run with

```bash
vep --cache --offline --dir_cache CACHE --assembly GRCh38 --fasta GRCh38.fa --input_file calls.vcf \
    --tab --af_1kg --check_existing --output_file vep.tab
fastvep annotate -i calls.vcf --gff3 genes.gff3 --cache-dir CACHE/homo_sapiens/116_GRCh38 \
    --output-format tab --af-1kg --check-existing -o fv.tab
python3 validation/compare_known_columns.py vep.tab fv.tab
```

What VEP does is not what one would write from a description of the columns, and each of these was found in the Perl and then confirmed by a run:

- `AF` is formatted `%.4f` and no other population is: `AFR_AF` is the text the cache has (`1e-05`, `0.00123456789`).
  An exact tie in binary rounds to even, as `printf` does (`0.03125` is `0.0312`).
- The frequency of an allele the record does not list is interpolated, as `1 - total`, for `AF` of a record with two alleles only, and printed with 15 significant digits (`0.7`, not `0.7000000000000001`).
  A record of three alleles that lists one has no `AF` for the others, and no other population is ever interpolated.
- The frequencies are read from *every* known variant of the site for every allele, in the order the cache holds them, and not in the order `Existing_variation` is sorted into.
  That includes a known variant that matched only another allele of the site: a record on the other strand, `T/G` with `AF=G:0.3`, attached to the input `A>C,T` through its `G`, gives the input's other allele `T` an `AF` of `0.7` under `--af` and nothing under `--af_1kg` alone, because the unmatched record is read only when `af` is among the groups (or `--max_af` puts it there).
  The row has an `AF` and no `Existing_variation`.
- `MAX_AF` leaves out the combined populations (`AF`, `gnomADe`, `gnomADg`) and lists every population that ties for it, once for each known variant that had it: `EUR,gnomADg_NFE,EUR,SAS`.
  When every population is 0 it is `0` and lists them all.
- `CLIN_SIG` is what the record says for the allele when it has `clin_sig_allele`, and replaces the list of the variant that `clin_sig` gives; once any known variant of the site has a `clin_sig_allele`, the `clin_sig` of those after it is not read.
  When the record's reference allele is not the input's, the alleles of `clin_sig_allele` are reverse-complemented first, and a deletion whose `clin_sig_allele` names the deleted bases (`AAGT:pathogenic`) finds nothing for the allele `-`.
- `SOMATIC` and `PHENO` list one `0` or `1` for each known variant, and are kept whole if any is `1` (`0,0,1`).
  `PUBMED` is every id of every known variant, not deduplicated.
- Every `--af*` flag, `--max_af` and `--pubmed` switch `--check_existing` on, which adds `CLIN_SIG`, `SOMATIC` and `PHENO`; the columns come in the order of `@FLAG_FIELDS`, after `HGVSg`, with `AF` once where `--af` and `--af_1kg` both name it.
- A trailing comma in the cache's `pubmed` (two records of HG002's chr21 and chr22 meet one) is no id: Perl's `split` drops the empty field at the end.
  fastVEP kept it, and the first real-data run found it.

Where fastVEP differs on purpose:

- **`--cache-dir` is `--check_existing`.**
  VEP's `--cache` is also its transcript cache, so it needs `--check_existing` to look the known variants up.
  fastVEP's `--cache-dir` is the variation cache and has no other use (gene models come from a GFF3), so naming one fills `Existing_variation` and, as VEP's `--check_existing` does, `CLIN_SIG`, `SOMATIC` and `PHENO`: a flag that could do nothing without a second one would be a trap.
  The frequencies and the PubMed ids need their flags, as in VEP.
- **The CSQ has no summary of the whole variant.**
  `AF`, `CLIN_SIG`, `SOMATIC`, `PHENO` and `PUBMED` were fields of fastVEP's `CSQ` before these columns, and held a summary of the whole variant: `AF` the first frequency of any known variant of the site (gnomAD's where it had one, which is not VEP's `AF`), `CLIN_SIG` the first clinical significance.
  For a record with several ALTs that puts one allele's value on every allele's row, which is a plausible-looking wrong answer, so they are gone.
  The five fields keep their place among the 49 of the default list, so that a reader that finds a field by its position keeps finding it, and hold VEP's per-allele value or nothing.
  A reader of the old `AF` now gives `--af-1kg` (VEP's `AF`) or `--af-gnomade`; one of the old `CLIN_SIG`, `SOMATIC` or `PHENO` changes only for a record whose alleles differ.
- **`CLIN_SIG` with two distinct allele-specific values.**
  VEP joins them with `;` in the order of a Perl hash (`keys %clin_sigs`), which differs from one run to the next, so its own output is not reproducible.
  fastVEP joins them in the order of the known variants as `Existing_variation` lists them, which is the same every time.
  The golden tests and the comparison script compare the parts as a set.
- **A flag that reads the cache, without one,** stops the run with the flag's name and `--cache-dir`: VEP would write `-` for every variant.
  `--everything` has the frequency flags only when a `--cache-dir` is given.
- **A cache without the columns of a flag** says so once, where it is opened (`--af-gnomadg asks for the gnomAD genomes frequencies, and this cache's variation_cols has none of their columns`) and then writes empty columns, as VEP does; VEP says nothing, except that `--af_gnomad` against a cache with no `gnomADe` is an error unless `--everything` is given.
- **JSON** carries none of these columns and says so for each flag given (VEP's JSON has a `colocated_variants` block this does not write).

Not accepted, and rejected as unknown flags where the other tool has them: `--af_esp` and `--af_exac` (VEP 116.3 has neither: `Config.pm` lists `af`, `af_1kg`, `af_gnomade`, `af_gnomad`, `af_gnomadg` and `max_af`, and its cache has no ESP or ExAC column), `--failed`, `--no_check_alleles`, `--exclude_null_alleles`, `--var_synonyms`, `--clinvar_somatic_classification`, `--check_frequency`, `--filter_common`, `--freq_*` and `--old_maf`.
`--af_gnomad`, VEP's old name for `--af_gnomade`, is accepted, and every flag of fastVEP whose name VEP spells with underscores is accepted in that spelling too.

One difference was found on the way and is not about these columns: the `Allele` cell of an input with lower-case bases is written lower-case by fastVEP and upper-case by VEP, which upper-cases its input.

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

Eighteen defects are gone, and the numbers above are what is left.
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
| `cDNA_position`, `CDS_position`, `Protein_position`, `EXON` and `INTRON` of a change with an end outside an exon, and of an insertion at an exon's edge | 9 cDNA cells of HG002 against 116.3; on the edge set 1,784, 1,268, 1,229, 255 and 603 cells | The one coordinate a change had was written where Ensembl writes `?-7` or `601-?`, an insertion's flank where it writes the two positions it sits between, and `EXON` and `INTRON` named the exon and intron of the start alone. The rule is Ensembl's mapper's, in Part 3 |
| An ALT of `*`, `.`, `<NON_REF>` or `<*>`, and an ambiguity code (`N`) in an ALT or a REF | on the allele set 100 rows only fastVEP wrote and 434 cells of 360 rows; HG002 6 rows | `*` was written as a row with no consequence, `.` as a `*`, and `<NON_REF>` or `<*>` beside a real ALT as a `frameshift_variant`; an `N` was a `MODERATE` `missense_variant` with `p.Lys31Xaa` and an `HGVSc` where Ensembl, which cannot translate it, writes `coding_sequence_variant` and nothing else. Part 3 has the table |
| The coding consequence terms were release 115's | 275 cells of the ClinVar sample against 116.3 (191 `Consequence`, 82 `IMPACT`); 447 and 176 of the edge set; 14,109 of 203,125 rows of VEP's RefSeq cache on start and stop codons | Release 116 rewrote `ref_eq_alt_sequence` and changed six more predicates; the port also had a separate path for a change reaching past the CDS and no initiator for a CDS that begins mid-codon. One transcription of `VariationEffect.pm`'s graph and of its edit tests replaces them, with the `include` gates of `Config.pm`, a window measured in the CDS, the UTR interval of `_before_coding` and the 3' UTR on the alternate codon. Part 3 |
| A lone indel was trimmed of its first shared base only | the world with each ALT its own record 197 rows of 7,213; the lone-indel set 2,950 `Allele` and 2,811 `Location` cells of 3,176 rows | `CAGT>CT` was a deletion of `AGT` at the next base, where Ensembl, which trims every shared base from both ends, writes the deletion of `AG`; the `Allele`, `Location`, `REF_ALLELE`, `VARIANT_CLASS`, positions and some consequence terms followed. Part 3 |
| `HGVSp` was a chain of rules read off Ensembl's output, and failed where a rule had no case: at a start codon, at a stop codon, and where the 3'-rule carries a change over an exon boundary | 99,670 cells of VEP's Ensembl cache run, 22,680 of its RefSeq cache run, 838 of the edge set, 2,837 of the ClinVar sample, 17 of HG002; now 12,604, 5,307, 237, 1,818 and 0, and every one of those has a cause (Part 3) | `fastvep_hgvs::protein` is a transcription of `hgvs_protein` and the routines it calls, and `fastvep_annotate::hgvsp_for_allele` hands it what Ensembl hands it: the change where the 3'-rule moved it, its protein position and codon window there, and the edited translation. The notation of a change that loses the initiator is `p.Met1?`, of a lost terminator `p.Ter583delextTer23`, of an insertion that the rule carries out of the coding sequence none. It also fixed the cause behind divergences 5 to 8, which was not the one this document gave |
| `HGVSc` at an end of a transcript: a change that overlaps it and runs past an end was written for the part that is inside, and an insertion beside it, which is not `within_feature`, was written all the same | 973 cells of the edge set, 230 `HGVS_OFFSET`, 13 of HG002 (9 and 4), 24 of the world | Ensembl clamps the change to the transcript, works out the 3'-rule for the change as it was and adds its length to the clamped coordinates, writes nothing for a change that is not a plain `overlap` of the transcript, and nothing for one the shift carries off the end. Part 3 |
| A transcript's symbol source, canonical flag, MANE accession and biotype were a GFF3's and not the cache's, and Ensembl's `cds_start_NF` and `cds_end_NF`, which its GFF3 export drops, were not there | `SYMBOL_SOURCE` 203,125 of 203,125 rows of VEP's RefSeq cache run and 1,994,540 of its 2,026,951 Ensembl rows; `CANONICAL` 26,757; `BIOTYPE` 1,677; `FLAGS` 141,958; `Consequence` 1,694 and `IMPACT` 566 | The RefSeq cache writes `EntrezGene`, the MANE Select transcript as canonical and the `gene_biotype` of the gene, and an Ensembl one the database's flags, which a GFF3 does not have. Read from the records, and from the sequence where there are none; the cache format moved on. [Against VEP's own cache](#against-veps-own-cache) |
| The first residue of a transcript that the annotation reads a non-ATG initiator as Met was the code table's | 32 `Amino_acids` and 32 `HGVSp` cells of VEP's Ensembl cache run, 48 and 114 of its RefSeq one | `CTG`, `GTG` and `ACG` starts are `M` in the cache (Ensembl's `initial_met`, RefSeq's `transl_except ... aa:Met`) and `L` or `V` in the table, so `Amino_acids` was `L/X` for `M/X` and `p.Leu1?` for `p.Met1?` |
