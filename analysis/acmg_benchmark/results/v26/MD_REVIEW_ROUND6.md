# Round 6 - note for the medical geneticist

Standalone, covering everything since the run your last comments were written against (v9). The
round-5 table went out and has not come back, so nothing here assumes you have read
`discordant_59_round5.xlsx` or the note that went with it: both are folded in below.

One attachment, **`discordant_64_round6.xlsx`**: only the calls still discordant, 64 rows where
fastVEP commits to a direction ClinVar 2-star+ contradicts. Everything you resolved has dropped
out. Your notes and rulings are carried onto the rows they belong to. **4 of the 64 you have never
been sent, 47 you were sent without ruling on, 13 carry your ruling.** Six blank columns on the
right, two of them added for your points 2 and 6 below.

Discordance with ClinVar is not itself a defect. Faithfulness to the guideline is the target. This
round adds a second target that did not exist before - **faithfulness to Ensembl's own annotation**,
now measurable, because both tools can be run over the same GFF3 and FASTA and compared row by row.

---

## Part 1 - your six points, and where each one landed

### 1. No BP4 for in-frame indels

In force. Of 6,277 in-frame variants, 24 still collect BP4, and every one sits in a splice region -
the one deliberate exception, where SpliceAI is calibrated for the question being asked. **Your
ruling on that carve-out would settle it; removing it is a one-line change.**

### 2. Common phenotypes, late onset, low penetrance

Partly built. BS2 tests a maximum credible prevalence, ClinVar's low-penetrance terms block benign
frequency evidence, and 304 published VCEP bars across 117 genes are built but not yet default.

What does not exist is a switch turning frequency evidence off for a disorder, because your three
examples want three different answers: GJB2 a looser published bar, SERPINA1's Z allele a per-allele
exception, MEFV possibly a gene-level exemption (Ghosh lists two alleles; whether the gene is exempt
is a judgement I do not have). The gene-disease attribute table now covers 6,857 pairs - inheritance
96 %, onset 91 %, prevalence 85 %, **penetrance 0 %**, because no public source publishes it. That
last field is exactly what Richards' BS1/BS2 precondition asks about. **Still waiting on your rule,
if there is one.**

### 3. BS2 on one or two homozygotes

Fixed. BS2 no longer counts homozygotes; it asks whether there are more individuals with no working
copy than the disorder could account for. 777 homozygotes are needed at gnomAD v4's cohort size.
The declined message is your argument verbatim: "consistent with a late-onset, reduced-penetrance or
variably expressive disorder rather than tolerance."

### 4. Founder alleles

BA1 and BS1 test the filtering AF, which removes the thin-sampling half of the problem; curation is
the only route to the rest. Worth knowing: gnomAD's grpmax deliberately excludes Finnish, Ashkenazi
and "remaining" as bottlenecked, so those founder alleles are invisible to BA1/BS1 by construction.
That is why Ghosh hand-listed ACADS and BTD, and why FMR1 `c.1126-1G>A` stays unresolved - its
Finnish frequency is 1.06 %, and grpmax cannot see it. **Still waiting on the list.**

### 5. Hypomorphs in trans

Your sharpest point, implemented. If a variant causes disease only in trans with a rarer null, being
common is what the disease model predicts, and unaffected homozygotes are too. Each exception entry
names which frequency criteria it suppresses; three ship blocking BA1, BS1 and BS2, all found in
your own table:

| Variant | AF | Hom | Was | Now |
|---|---:|---:|---|---|
| GAA `c.-32-13T>G` (late-onset Pompe) | 5.4e-3 | 23 | Likely benign | VUS |
| CFTR `c.1210-11T>G` (5T) | 9.8e-3 | 27 | Likely benign | VUS |
| SPTA1 `c.4339-99C>T` (αLELY) | 6.6e-3 | 40 | Likely benign | VUS |

None reaches 5 %, which is why a BA1-only list never touched them. None reaches Likely pathogenic
either, and should not - the deciding evidence is PM3. **RBM8A `c.-21G>A` belongs here; I have not
added it without your confirmation of the allele.**

### 6. Functional data over BP7

Implemented from your OCA2 case: curated PS3/BS3 suppresses PP3, BP4 and BP7 by name. Your point
that predictors are often trained on that data is why suppression runs both ways. OCA2 `c.1080C>T`
goes from Likely benign to VUS on one line of curated PS3. We will not mine PubMed for PS3/BS3 - a
wrong PMID in a clinical report is worse than no evidence.

### The defect your comments uncovered

When gnomAD rejects a site, we withheld BA1, BS1, BS2 **and** PM2 - even-handed in intention,
one-sided in effect. The benign side loses up to Standalone plus two Strong; the pathogenic side
loses one Supporting; and PVS1, read from the alignments just declared untrustworthy, stayed at full
strength. RAI1 `c.840del` - 48,739 gnomAD homozygotes, ClinVar Benign - was Likely pathogenic on
PVS1 alone, with seventeen more like it. Now Uncertain. That fix alone took opposite-direction calls
from 76 to 59.

---

## Part 2 - what has happened since, and it is all in the annotation layer

Two further rounds landed. Neither touched a criterion, a threshold or a database. Both changed what
consequence a variant is given and what HGVS name it gets, and the criteria then fire differently
because they are reading a different consequence.

The reason this became worth doing: fastVEP and Ensembl VEP can now be run over the same GFF3 and
FASTA and compared row by row. Over a 6,600-variant sample of your truth set - 150,725 transcript
rows - fastVEP disagreed with Ensembl on a quarter of all amino-acid changes.

### 7. The consequence is now Ensembl's answer, not an approximation of it

Four separate pieces of code decided what a change inside a coding sequence was - one for
substitutions, one for insertions, one for deletions, one for changes that do both - and each had a
different, partial version of the rules. They are now one calculation: the whole codons the variant
touches, translated on both sides, with every term read off the resulting peptide pair the way
Ensembl reads it.

| | before | after |
|---|---:|---:|
| Amino-acid change agrees with Ensembl VEP | 77.7 % | **100 %** |
| Codon change agrees | 72.5 % | **100 %** |
| Splice terms agree | 99.0 % | **100 %** |
| Whole consequence set agrees | 98.0 % | **99.9 %** |
| HGVSp agrees (coding rows) | 83.4 % | **98.8 %** |
| HGVSc agrees (coding rows) | 87.7 % | **99.5 %** |

The practical effect on your rows: a deletion whose *second* base landed on a splice site used to be
missed entirely; a change that both deletes and inserts bases was typed from the direction of its
length change rather than from the protein it produces; and the extended splice terms
(`splice_donor_5th_base_variant` and the rest) were being hidden behind `splice_donor_variant`
instead of reported beside it, which is what Ensembl does.

### 8. Two defects, both of which reported a truncating variant as something milder

**PVS1 was reading the wrong end of a span.** PVS1's splice track applies only to the canonical
±1/±2 bases, and the gate reads that offset out of the HGVSc string. A change like OAT
`c.764_771+9delinsTTAGCTGTTTGTATCACACCA` starts in the exon and runs nine bases into the intron - it
removes +1 and +2 on the way - but the offset was taken from the *endpoint*, `+9`, so PVS1 concluded
the variant was deep intronic and stood down. **Fourteen ClinVar-pathogenic donor and acceptor
deletions lost PVS1 to this**: OAT, CEP290, MSH6, ATM (×2), MLH1, MYO7A, PYGL, IFT140, CLN5, HADHA,
CYP4V2, SERAC1, BRCA1. The gate now reads the smallest offset the span *covers*.

It only started firing at all because HGVSc improved. Before this round those spans had no HGVSc, so
there was no offset to misread - the same shape of problem as the gnomAD-rejection defect above:
withholding a field is not neutral.

**A duplication was being written across an exon boundary.** Offsets count outward from the exon and
do not run through zero, so stepping four bases back from `+1` does not reach `-2`. The
insertion-to-duplication conversion did exactly that, naming three bases in the previous intron.

### 9. One place where we now deliberately disagree with Ensembl

Ensembl has a clause that asks whether a replacement keeps the residue it starts on and introduces a
terminator *anywhere* - not whether the annotated stop codon survived. It fires on any frameshift
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
SPG11, TTN, EXT1, CLCN7, COL7A1, MANBA, TMEM127 and TRNT1. Ensembl's own amino-acid string for the
first is `N/N*X`: Asn becomes Asn, then a terminator, then an incomplete codon. Nothing in it says
the annotated stop at residue 1863 survived. fastVEP does not reproduce it, so those 34 keep PVS1.
Reversing that is a one-line change and would move all 34 out of P/LP.

Every other difference between the two tools is enumerated in `docs/VEP_DIVERGENCE.md`, with counts.
This is the only one that changes a call.

---

## Where the numbers went

| | v9 (your comments) | v20 (round 5) | **v26 (now)** |
|---|---:|---:|---:|
| Exact match | 403,507 | 425,145 | **425,322** |
| Same-direction | 477,965 | 522,101 | **522,473** |
| Exact-match rate | 59.9 % | 63.1 % | **63.1 %** |
| Same-direction rate | 71.0 % | 77.5 % | **77.6 %** |
| P/LP recall | 45.58 % | 62.03 % | **62.42 %** |
| B/LB recall | 51.93 % | 62.53 % | **62.53 %** |
| ClinVar P/LP called B/LB | 18 | 30 | **30** |
| ClinVar B/LB called P/LP | 28 | 29 | **34** |
| Opposite-direction calls | 46 | 59 | **64** |

**15,803 more ClinVar-pathogenic variants reach P/LP than at v9**, and 30,171 more benign ones reach
B/LB.

The opposite-direction count rose from 46 to 64, and that needs saying plainly rather than buried.
Most of it is the price of decisiveness: at v9 fastVEP committed to a direction on 29.5 % of the set
and gave an exact `Pathogenic` to 146 of your 79,823 ClinVar-Pathogenic variants. It now commits on
36.5 %. Per decisive call the opposite rate has barely moved - 0.023 % at v9, 0.024 % at v20,
0.026 % now - so the extra opposite calls are the tail of 47,000 additional decisive ones, not a
new class of error. That framing is not an excuse for the six new ones, which are itemised below.

Between v20 and v26 exactly **one** row resolved and **six** are new, all six in the same direction
and all six splice-adjacent:

| Gene | Variant | ClinVar | Does Ensembl agree the consequence is HIGH? |
|---|---|---|---|
| NFKB2 | `10:102397405` 24 bp del | Benign (2-star) | yes |
| FOXP2 | `7:114631524 CCAG>C` | Benign/Likely_benign | yes |
| LMX1B | `9:126693802` 33 bp del | Likely_benign | yes |
| MSH6 | `2:47806358 G>GGTAT` | Likely_benign | yes |
| DSP | `6:7565520 C>CGTAA` | Likely_benign | yes |
| MSH6 | `2:47806651 G>GGTAACTAACTAACTA` | Likely_benign | no - this is the divergence in (9) |

The first three genuinely remove part of a donor or acceptor dinucleotide, so the disagreement there
is with ClinVar's classification rather than with the annotation. The last three are one shape and
are question A below.

Separately, **ten** ClinVar-pathogenic variants moved out of P/LP into **VUS** - not into B/LB, so
they are not in the count above. All ten are single-base duplications two bases into the intron,
`c.N+2dup`, in MECR, COL11A1, PCDH15, SLX4, BRIP1, SLC7A9, ERCC2, NEB, GLB1 and NPHP3. Ensembl
agrees with the new `splice_region_variant` call on all ten: duplicating the base after `+2` leaves
the `GT` intact. This is the class you already reached in the *classifier* in round 2 (PTEN
`c.802-2dupA`, BRIP1 `c.2258-2dup`); what changed is that the consequence itself is now right rather
than PVS1's gate carrying the whole burden. Whether a `+2` duplication should still reach PVS1 is
question A too.

---

## Three questions for you

### A. A pure insertion or duplication at the very edge of an exon

Two faces of one question, thirteen variants.

- MSH6 `c.3801+3_3801+6dup`, MSH6 `c.4001+1_4001+25dup`, DSP `c.939+3_939+6dup` - all ClinVar
  Likely_benign, all now Likely_pathogenic on PVS1.
- The ten `c.N+2dup` above - all ClinVar Pathogenic or Likely_pathogenic, all now VUS.

Read where the VCF puts it, an insertion at an exon's final base shifts the reading frame. Read under
the HGVS 3'-rule, it slides into the intron and becomes a small intronic duplication with no coding
effect. Ensembl's own HGVSc takes the second reading while its consequence takes the first. ClinVar's
submitters took the second for the MSH6 and DSP alleles and the first for the `+2dup` ones.

PVS1 already has a gate for exactly this ambiguity, from your round-2 review: a pure insertion at a
canonical splice site does not get PVS1 without positive SpliceAI support, because inserting a base
beside the dinucleotide need not destroy it. That gate lives on PVS1's **splice** track. The MSH6 and
DSP variants reach PVS1 on its **frameshift** track, where it does not apply.

**Should the same reasoning extend to the frameshift track when the insertion sits on an exon's last
base? And should a `+2` duplication reach PVS1 in some genes anyway?** Both are criteria decisions,
and an exception list or a SpliceAI-positive path would implement either.

### B. An inversion over a splice acceptor

NBN `8:89958851`, 18 bp, `c.995-23_998inv`, ClinVar Pathogenic/Likely_pathogenic. fastVEP now calls
it `coding_sequence_variant,intron_variant`, MODIFIER, VUS - down from `splice_acceptor_variant`,
HIGH, LP.

Ensembl matches a variant against a splice site over the bases that actually *differ* between the
two alleles, not over the whole span, and an inversion leaves its palindromic positions unchanged.
Here the differing bases straddle the acceptor without covering it, so Ensembl calls it MODIFIER and
fastVEP now agrees. An 18 bp inversion spanning an acceptor is nonetheless very likely to disrupt
splicing.

**Should an inversion over a splice site keep PVS1?** A criteria rule, or a SpliceAI-driven one - the
annotation is doing what Ensembl does.

### C. The deliberate divergence in (9)

Keep it, or match Ensembl? Keeping it holds 34 ClinVar-pathogenic truncating variants in P/LP.

---

## Against your own rulings

Unchanged, and I would rather report the ceiling than tune toward it. Of the 76 rows where you stated
a call: we agree on 1, defer to VUS on 55, remain opposite on 11. **None of the annotation work moved
any of them.** The notes on those rows cite segregation, phase and functional data, and none of that
is derivable from a variant and a database.

**Most useful next:** the three questions above; then the hypomorph list, the founder alleles, your
rule for point 2 if there is one, and your ruling on the in-frame splice-region BP4 carve-out; then
the 51 unruled rows in the attachment.
