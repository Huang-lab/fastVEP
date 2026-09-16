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

---

## Part 2 - gaps, where Ensembl is right

These are fastVEP's remaining defects, not disagreements.

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

### HGVSg

`HGVSg` is not in the CSQ field set, so `validation/compare_rows.py` does not rate it; it was
measured separately against the same HG002 sample, run with `--hgvsg` on VEP's side and JSON
output on fastVEP's, comparing the set of `g.` descriptions each tool writes per record.
They agree on **13,525 of 13,535 records (99.926 %)**.

Every one of the ten is the same gap: Ensembl applies the 3'-rule over the genome and fastVEP does
not, so an insertion sits where the VCF put it rather than where HGVS asks for it, and never
collapses to a `dup`.
`9:905488 C>CTGTGTGTG` is `9:g.905488_905489insTGTGTGTG` here and `9:g.905507_905514dup` there.
That shift is not the transcript-direction one `HGVSc` gets: it runs along the genome whichever
way the gene points, and it is a separate piece of work.

Two of the ten are a record where VEP's own `HGVSg` disagrees with VEP's own `Allele`: for
`5:100381015 CATAA>AATAA,C` it writes `Allele` `A` and `HGVSc` `n.188+196G>T` - the substitution
that it is - beside `HGVSg` `5:g.100381015delinsAATAA`.
fastVEP writes `5:g.100381015C>A`, which is the same change as its own `HGVSc`.

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
What that means is that no consumer of these fields may recover a position by parsing one: the
ACMG criteria that did are the subject of a run-versions entry
([`analysis/acmg_benchmark/RUN_VERSIONS.md`](../analysis/acmg_benchmark/RUN_VERSIONS.md)), and
they now take the offset from the transcript instead.

**`splice_region_variant` is decided by the last differing region, not the union.**
`_intron_effects` (`BaseTranscriptVariationAllele.pm`, l. 215) assigns rather than or-assigns it
inside its region loop, so a change whose first differing base is in the splice region and whose
last is not comes out with no splice term.
Reproducing the tidier rule instead put a `splice_region_variant` on 63 rows of that run which VEP
calls plain missense.

---

## Fixed since this document was last measured

Eight defects are gone, and the numbers above are what is left.
Each was measured the same way before and after, on the same samples: the row counts below are
what the fix removed.

| What was wrong | Rows it cost | What it was |
|---|---:|---|
| Each allele of a multi-allelic record kept the bases it shared with the reference | 801 of 824 genome-wide `HGVSc` rows, and 1 of 11 `HGVSg` records | A site is trimmed once for the whole record, so an allele that is a clean deletion arrives as a replacement and was written as one: `c.553-61775_553-61770delinsCTT` for TTC28's `c.553-61772_553-61770del`. Clipped now where the description is built, which is the second of the two places Ensembl clips. What the clip does *not* reach is in Part 2 |
| The HGVSc 3'-shift ran on the *spliced* sequence | 563 of 602 `HGVSc` rows | It could not follow a repeat out of an exon into the intron (`c.220del` for VEP's `c.220+1del`), and it would walk a deletion over a splice junction into the next exon, naming a block that is contiguous in the mRNA and not in the DNA a `c.` description is numbered against (BRCA1 `c.71_81del` for VEP's `c.70_80del`). The shift now runs on the genome, bounded by the transcript, and each end is written where it lands |
| A lost terminator was written `p.Glu486ext*?` | 30 `HGVSp` rows | The extension is a fact about the edited transcript, and `ext*?` is the form HGVS keeps for one nobody can measure. Now divergence 7 |
| The 3'-rule's rotation carried residues out from behind a terminator | 23 `HGVSp` rows | `insValAlaLeuAspThrTerVal` named a Val the protein never has |
| `ref_eq_alt_sequence`'s second clause was read as a position test | 18 consequence rows | It is a comparison of sequences: an insertion whose residues repeat what follows them holds the clause anywhere the repeat reaches the terminator, not only on the last residue (MSH6 `c.4106_4108dup`, ENST00000936511) |
| `start_lost` stopped at its first test | 4 rows of the in-frame deletion set, all four `IMPACT` | Ensembl falls through to the peptide test when the coordinate test declines because the edit is an in-frame indel, and an in-frame deletion that removes the initiator is a start loss there (KCNA2 `c.3_11del`, ENST00000639048, HIGH not MODERATE) |
| An unknown residue silenced `missense_variant` | 1 consequence row, `IMPACT` | `QX/HX` resolves one residue and not the next; Ensembl reports `missense_variant,coding_sequence_variant`, and reporting the second alone made a real missense MODIFIER (IL7R ENST00000515665) |
