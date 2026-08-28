# Round 6 - note for the medical geneticist

Standalone. Round 5's table went out and has not come back, so nothing here assumes you have
looked at `discordant_59_round5.xlsx`; that round's 59 rows are carried forward into the
attachment below, and everything since is added to them.

Two rounds of work have landed since the numbers you last commented on (v20). Neither touched a
criterion, a threshold or a database. Both were in the **annotation** layer - what consequence a
variant is given, and what HGVS name it gets - and in two places where the classifier was reading
that annotation wrongly. The criteria then fire differently because they are reading a different
consequence.

One attachment, **`discordant_64_round6.xlsx`** (built at
`data/benchmark/discordance_review/`, which is not tracked - the workbooks are sent, not
committed): only the calls still discordant, 64 rows where
fastVEP commits to a direction ClinVar 2-star+ contradicts. Everything you resolved has dropped
out. Your notes and rulings are carried onto the rows they belong to. **4 of the 64 you have never
been sent, 47 you were sent without ruling on, 13 carry your ruling.** The same six blank columns
on the right as last time.

---

## What changed since v20

Discordance with ClinVar is still not itself the target; faithfulness to the guideline is. What is
new this time is that there is now an **external check on the annotation itself**. Everything below
was measured against real Ensembl VEP 115.1 running on the same GFF3 and FASTA fastVEP reads, over
150,725 transcript rows from a 6,600-variant sample of your truth set.

### 1. The consequence a variant gets is now Ensembl's answer, not an approximation of it

fastVEP had four separate pieces of code deciding what a change inside a coding sequence was - one
for substitutions, one for insertions, one for deletions, one for changes that do both - and each
had a different, partial version of the rules. They are now one calculation: the whole codons the
variant touches, translated on both sides, with every term read off the resulting peptide pair the
way Ensembl reads it.

| | before (v20) | after (v26) |
|---|---:|---:|
| Amino-acid change agrees with Ensembl VEP | 77.7 % | **100 %** |
| Codon change agrees | 72.5 % | **100 %** |
| Splice terms agree | 99.0 % | **100 %** |
| Whole consequence set agrees | 98.0 % | **99.9 %** |

The practical effect on your rows: a deletion whose *second* base landed on a splice site used to
be missed entirely, and a change that both deletes and inserts bases was typed from the direction
of its length change rather than from the protein it produces.

### 2. Splice terms are no longer hidden behind each other

`splice_donor_5th_base_variant`, `splice_donor_region_variant` and
`splice_polypyrimidine_tract_variant` sit *beside* `splice_donor_variant` in Ensembl's output.
fastVEP was treating them as a fallback, reporting only the most severe. `intron_variant` likewise
sits beside an exonic term. This is cosmetic for the classification but it is what let defect (4)
below become visible.

### 3. HGVS follows the 3'-rule before it names anything

An insertion inside a repeat is described at its most 3' position, and only there does it read as a
duplication: `c.4177_4182dup` where fastVEP wrote `c.4172_4173insCACCAG`. A change spanning an exon
boundary is written over both of its ends. HGVSp agreement with Ensembl went from 83 % to 99 % on
coding rows, HGVSc from 88 % to 99 %.

This matters beyond cosmetics because PVS1 reads the intronic offset off the HGVSc string - see (4).

### 4. Two defects, both of which reported a truncating variant as something milder

**PVS1 was reading the wrong end of a span.** PVS1's splice track applies only to the canonical
±1/±2 bases, and the gate reads that offset out of the HGVSc. A change like OAT
`c.764_771+9delinsTTAGCTGTTTGTATCACACCA` starts in the exon and runs nine bases into the intron -
it removes +1 and +2 on the way - but the offset was taken from the *endpoint*, `+9`, so PVS1
concluded the variant was deep intronic and stood down. **Fourteen ClinVar-pathogenic donor and
acceptor deletions lost PVS1 to this**: OAT, CEP290, MSH6, ATM (×2), MLH1, MYO7A, PYGL, IFT140,
CLN5, HADHA, CYP4V2, SERAC1, BRCA1. The gate now reads the smallest offset the span *covers*.

It only started firing at all because HGVSc improved: before this round those spans had no HGVSc,
so there was no offset to misread.

**A duplication was being written across an exon boundary.** Offsets count outward from the exon
and do not run through zero, so stepping four bases back from `+1` does not reach `-2`. The
insertion-to-duplication conversion did exactly that, naming three bases in the previous intron.

### 5. Where we now deliberately disagree with Ensembl

One place, and it is the one I would most like your view on.

Ensembl has a clause that asks whether a replacement keeps the residue it starts on and introduces
a terminator *anywhere* - not whether the annotated stop codon survived. It fires on any frameshift
whose new stop lands in the codon just after the insertion point, and both `frameshift_variant` and
`stop_gained` defer to it. So real VEP 115.1 reports:

| Variant | ClinVar | Ensembl VEP 115.1 | fastVEP |
|---|---|---|---|
| BRCA1 `c.5030_5033dup` | Pathogenic (3-star) | `inframe_insertion,stop_retained_variant`, MODERATE | `stop_gained,frameshift_variant`, HIGH |
| BRCA1 `c.1499_1508dup` | Pathogenic (3-star) | same, MODERATE | HIGH |
| BRCA2 `c.3205_3206insAATTGCAGTCAATTAATAT` | Pathogenic (3-star) | same, MODERATE | HIGH |
| TP53 `c.895_919dup` | Pathogenic/Likely_pathogenic | same, MODERATE | HIGH |
| ITGB3 `c.122_125dup` | Pathogenic (3-star) | same, MODERATE | HIGH |

34 variants in your truth set are of this shape, also in SDHB, USH2A, ATM, MSH6, BRIP1, GATA2,
SPG11, TTN, EXT1, CLCN7, COL7A1, MANBA, TMEM127 and TRNT1. fastVEP does not reproduce it, so those
34 keep PVS1. Reversing that is a one-line change and would move all 34 out of P/LP.

Ensembl's own amino-acid string for the first of these is `N/N*X`: Asn becomes Asn, then a
terminator, then an incomplete codon. Nothing in it says the annotated stop at residue 1863
survived. Every other divergence between the two tools is listed in `docs/VEP_DIVERGENCE.md`; this
is the only one that changes a call.

---

## Where the numbers went

| | v20 | **v26** |
|---|---:|---:|
| Exact match | 425,145 | **425,322** |
| Same-direction | 522,101 | **522,473** |
| Exact / same-direction rate | 63.1 % / 77.5 % | 63.1 % / 77.6 % |
| P/LP recall | 62.03 % | **62.42 %** |
| ClinVar P/LP called B/LB | 30 | **30** |
| ClinVar B/LB called P/LP | 29 | **34** |
| Opposite-direction calls | 59 | **64** |

**370 more ClinVar-pathogenic variants now reach P/LP**, and nothing moved from P/LP into B/LB.

The opposite-direction count went the wrong way by five, and it is worth being precise about which
five. Exactly **one** row resolved and **six** are new, all six in the same direction - called
pathogenic where ClinVar says benign - and all six splice-adjacent:

| Gene | Variant | ClinVar | Ensembl VEP agrees the consequence is HIGH? |
|---|---|---|---|
| NFKB2 | `10:102397405` 24 bp del | Benign (2-star) | yes |
| FOXP2 | `7:114631524 CCAG>C` | Benign/Likely_benign | yes |
| LMX1B | `9:126693802` 33 bp del | Likely_benign | yes |
| MSH6 | `2:47806358 G>GGTAT` | Likely_benign | yes |
| DSP | `6:7565520 C>CGTAA` | Likely_benign | yes |
| MSH6 | `2:47806651 G>GGTAACTAACTAACTA` | Likely_benign | no - this is the divergence in (5) |

The first three genuinely remove part of a donor or acceptor dinucleotide, so the disagreement is
with ClinVar's classification rather than with the annotation. The last three are the same shape as
each other and are the question in the next section.

Separately, **ten** ClinVar-pathogenic variants moved from P/LP to **VUS** rather than to B/LB, and
those are in the round-5 list you have not yet seen: single-base duplications two bases into the
intron, `c.N+2dup`, in MECR, COL11A1, PCDH15, SLX4, BRIP1, SLC7A9, ERCC2, NEB, GLB1 and NPHP3.
Ensembl agrees with the new `splice_region_variant` call on all ten - duplicating the base after
`+2` leaves the `GT` intact. Whether a `+2` duplication should still reach PVS1 is a criteria
question, and it is the same reasoning as the next section.

---

## Three questions for you

### A. A pure insertion at the last base of an exon

MSH6 `c.3801+3_3801+6dup`, MSH6 `c.4001+1_4001+25dup` and DSP `c.939+3_939+6dup`, all ClinVar
Likely_benign, all now Likely_pathogenic on PVS1.

Read where the VCF puts it, each of these is an insertion at the exon's final base, which shifts
the reading frame. Read under the HGVS 3'-rule, it slides into the intron and becomes a small
intronic duplication with no coding effect. Ensembl's own HGVSc takes the second reading; its
consequence takes the first. ClinVar's submitters evidently took the second.

PVS1 already has a gate for exactly this ambiguity, from your round-2 review: a pure insertion at a
*canonical splice site* does not get PVS1 without positive SpliceAI support, because inserting a
base beside the dinucleotide need not destroy it. That gate lives on PVS1's splice track. These
three reach PVS1 on its **frameshift** track, where it does not apply.

**Should the same reasoning extend to the frameshift track when the insertion sits on an exon's
last base?** Same question, differently dressed, as the ten `c.N+2dup` variants above.

### B. An inversion over a splice acceptor

NBN `8:89958851`, an 18 bp inversion, `c.995-23_998inv`, ClinVar Pathogenic/Likely_pathogenic.
fastVEP now calls it `coding_sequence_variant,intron_variant`, MODIFIER, VUS - down from
`splice_acceptor_variant`, HIGH, LP.

Ensembl matches a variant against a splice site over the bases that actually *differ* between the
two alleles, not over the whole span, and an inversion leaves its palindromic positions unchanged.
Here the differing bases straddle the acceptor without covering it, so Ensembl calls it MODIFIER
and fastVEP now agrees. An 18 bp inversion spanning an acceptor is nonetheless very likely to
disrupt splicing.

**Should an inversion over a splice site keep PVS1?** That would be a criteria rule, or a
SpliceAI-driven one - the annotation is now doing what Ensembl does.

### C. The deliberate divergence in (5)

Keep it, or match Ensembl? Keeping it holds 34 ClinVar-pathogenic truncating variants in P/LP.

---

## Still open from round 5

Nothing below has moved, because none of it is computable from a variant and a database:

1. **The carve-out where in-frame indels in a splice region may still collect BP4** - 24 of 6,277.
   Your ruling settles it; removing it is a one-line change.
2. **A rule for common phenotypes, late onset and low penetrance.** Your three examples want three
   different answers (GJB2 a looser published bar, SERPINA1 a per-allele exception, MEFV possibly a
   gene-level exemption). The gene-disease attribute table is built - inheritance 96 %, onset 91 %,
   prevalence 85 %, penetrance 0 %, because no public source publishes it.
3. **The hypomorph list.** Three ship (GAA `c.-32-13T>G`, CFTR `c.1210-11T>G`, SPTA1 `c.4339-99C>T`),
   each naming which frequency criteria it suppresses. RBM8A `c.-21G>A` is waiting on your
   confirmation of the allele.
4. **Founder alleles.** gnomAD's grpmax deliberately excludes Finnish, Ashkenazi and "remaining" as
   bottlenecked, so those alleles are invisible to BA1/BS1 by construction. FMR1 `c.1126-1G>A`
   (Finnish 1.06 %) is still unresolved.
5. **The 47 rows in the attachment you were sent but have not ruled on**, plus the 4 that are new.

On your own rulings, the position is unchanged from v20 and I would rather say so than tune toward
it: of the 76 rows where you stated a call, we agree on 1, defer to VUS on 55, and remain opposite
on 11. None of this round's work moved any of them - the notes on those rows cite segregation,
phase and functional data, and none of that is derivable from a variant and a database.

**Most useful next:** the three questions above, then items 1-4, then the 51 unruled rows in the
attachment.
