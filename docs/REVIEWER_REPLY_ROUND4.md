# Reply to the round-3 review: the BA1 exception list

Thank you for Table 1 of Ghosh 2018. Two things came of it, and the second is the more important.

**The list was already shipped. It had never once fired.**

All nine variants were in the defaults, cited to the right paper. The matcher compared the full HGVS string the pipeline emits, `ENST00000357618.10:c.187C>G`, against the bare `c.187C>G` the list holds. It never matched, on any variant, since the day it was written. Unit tests passed because they fed the matcher bare `c.` tokens the pipeline does not produce.

So four of the nine were being called **Benign on frequency alone**, which is precisely what the SVI list exists to prevent:

| Variant | Was | Now |
|---|---|---|
| HFE c.187C>G (p.His63Asp) | Benign | VUS |
| MEFV c.1105C>T (p.Pro369Ser) | Benign | VUS |
| PIBF1 c.1214G>A (p.Arg405Gln) | Benign | VUS |
| ACAD9 c.-44_-41dupTAAG | Benign | Likely benign |

Verified end to end on all nine, not in a unit test. HFE c.845G>A (p.Cys282Tyr) was already protected, but by a different mechanism - ClinVar's own `Pathogenic, low penetrance` term - not by the list.

## Two matching problems your table let us fix properly

Matching on HGVS at all was the underlying mistake, and your table has the columns that fix it. Entries are now keyed on **genomic coordinate**, with the HGVS kept as a fallback and your ClinVar variation IDs carried for provenance, because two of the nine cannot be matched on their published HGVS at all:

- **BTD** is `c.1330G>C` on NM_000060.4 in your table and `c.1270G>C` on the Ensembl canonical transcript in our output. A `c.` token only identifies a variant relative to the transcript it was written for.
- **ACAD9** is `c.-44_-41dupTAAG` in your table and `c.-45_-44insTAAG` in ours. Same variant, two valid spellings.

Both now match on coordinate. The GRCh37 positions in your table were lifted to GRCh38 through the ClinVar IDs you also supplied.

## Why our benchmark could not have caught this

We re-ran the full 673,660-variant ClinVar 2-star+ set with the fix. **Not one call changed.** That is not a disappointing result; it is the explanation for how a dead exception list survived this long.

Seven of the nine variants are `Conflicting classifications of pathogenicity` in ClinVar, so they carry no consensus label and are not in the truth set at all. Of the two that are scored, neither call moves. The exception list exists for variants where expert opinion is divided or where a variant is common *and* pathogenic - exactly the variants a consensus-truth benchmark discards. **Our headline metric is structurally blind to this whole class**, and we have written that into the limitations rather than leaving it as folklore.

This is the strongest argument yet for your kind of review over any amount of aggregate scoring.

## Your four rows in these genes are different alleles, and the list does not reach them

I want to be exact, because it would be easy to imply more than is true. Your review flagged GJB2, ACADS, MEFV and BTD - and in every case a **different variant** from the one Ghosh lists:

| Your row | On the Ghosh list? | Your note |
|---|---|---|
| GJB2 c.167del | No; the list has c.109G>A | "hearing loss prevalence is high, not a rare phenotype, GJB2 LoF variants are high in population frequency" |
| ACADS c.319C>T | No; the list has c.511C>T | "founder variant in AJ population; 9 hom in gnomAD; incomplete penetrance" |
| MEFV c.2177T>C | No; the list has c.1105C>T and c.1223G>A | "common phenotype, low penetrance, BS1 and BS2 do not apply" |
| BTD c.\*159G>A | No; the list has c.1330G>C | (no note) |

The mechanism now works and is empty of the entries your cases need. That converts your review into something very concrete: **name the alleles and they ship as defaults, with your attribution.** Gene, HGVS or coordinate, and one line of reasoning is enough. The same applies to the ACADS-in-Ashkenazi, RSPH9, C9 and OCA2 alleles you raised in round 2.

There is a second reason the Ghosh list is thinner than it looks for founder alleles, and it is worth your attention. BA1 and BS1 test gnomAD's **grpmax** filtering allele frequency, and gnomAD deliberately excludes the Finnish, Ashkenazi and "remaining" groups from grpmax because they are bottlenecked. So a Finnish or Ashkenazi founder allele is invisible to that statistic by construction - which is exactly why Ghosh had to list ACADS and BTD by hand, both annotated in your table as ">5% MAF only in Finnish". A curated entry is not a workaround here; it is the only mechanism available.

## Your FMR1 question, answered precisely

> "210 hemizygotes and 1 homozygote in gnomAD v4; not disease causing; I am not sure why the population frequencies were not applied."

They were applied. Both codes evaluated it and both declined, and the numbers say why:

```
FMR1 c.1126-1G>A
  BS2  1 homozygote and 173 hemizygotes among 722,456 individuals
       95% lower bound 2.12e-4, against the 1e-3 prevalence bar   -> declines
  BS1  grpmax filtering AF 5.24e-4, against the 5e-3 bar          -> declines
  PM2  AF 5.61e-4 is above the 4e-5 bar                           -> also declines
```

Hemizygote counting is working - BS2 names the 173 explicitly, where before this whole class was invisible. (173 rather than your 210 because we count `AC_XY`, gnomAD's XY-sample allele count.)

**This variant is the frequency dead zone, exactly.** At 5.24e-4 it is fifteen times too common for PM2 and ten times too rare for BS1, so it receives no frequency evidence in either direction and PVS1 stands alone at Likely pathogenic. It is the single cleanest example of the question we asked you last round, and it is still the question we would most like answered.

Its Finnish frequency is 1.06%, which would settle it - and grpmax excludes Finnish, per the paragraph above.

## What else moved, and the updated workbook

Two attachments, both your own annotated file with columns appended on the right. Your notes, highlighting and row order are untouched.

- **`discordant_122_round4.xlsx`** - current calls under the shipping defaults.
- **`discordant_122_round4_with_vcep_bars.xlsx`** - the same rows with published ClinGen VCEP frequency thresholds loaded, which is a candidate default we have not turned on. See below.

Of your 122 rows: **90 are now VUS**, 32 still commit to a direction that contradicts ClinVar. Against your own adjudication, on the 62 rows where your highlighting states a call, we now agree on 1, defer to VUS on 50, and remain opposite on 11. Deferring where you would commit is the honest automated answer when the deciding evidence is case-level or literature-level, but it is not a good answer, and the 11 are listed below.

Since round 3, one criterion moved a large number of calls, though none of yours: PS1's splicing path (Walker 2023) had never fired because nothing built the index it reads. It now fires on 4,776 variants, 98.8% of them pathogenic in ClinVar, and moved 3,793 calls from Likely pathogenic to Pathogenic - 3,363 of them correctly. The mechanism is arithmetic and worth knowing: PVS1 is 8 points and PM2_Supporting is 1, and Pathogenic needs 10.

## Published VCEP thresholds: 15 of your rows, 1 that moves

We now read BA1, BS1 and PM2 bars straight from the ClinGen CSpec Registry - 304 bars across 117 genes, each traceable to the sentence it was parsed from. Fifteen of your rows are in genes with a published bar: ENG ×2, MSH2, PTEN ×2, APC, HBA2 ×3, MYO15A, MSH6 ×2, ATM ×2, BRCA2.

Loading them moves exactly one, and in the right direction:

```
MSH2 c.944G>T   truth Likely benign   VUS -> Likely benign
                (Lynch VCEP BS1 = 1e-4; grpmax AF 1.76e-3)
```

None of the other fourteen move, and none break.

**We have not made this a default, and the reason is in your territory.** Across the whole benchmark it buys 2,028 more correct Benign calls but introduces four new false-benign ones, and two of those four are founder variants - ANO5 and CAPN3, both limb-girdle muscular dystrophy - that the panels excluded *in the prose of the specification*, which our parser cannot read. 32 of the 286 usable bars carry an exclusion of that shape. Applying a tightened bar without its exceptions is applying half a specification, and the half omitted is the half that protects founder alleles. This is the same problem as the Ghosh list, one layer up.

## The 11 rows where we still contradict you

| Gene | Variant | You | Us | Why it is still open |
|---|---|---|---|---|
| CBS | c.844_845ins…, c.832_833ins… ×2 | benign | LP | 68 bp insertions; PVS1 fires. Your note says mapping is unreliable for insertions this large and that ClinVar's Benign looks wrong, which contradicts the green highlight on those rows. **Which did you mean?** |
| CHD7 | c.5429G>A | benign | LP | 30 hets; PM1 has 3 pathogenic and 0 benign neighbours, so the "without benign variation" test passes. Your objection is gene-level - variants are spread out - and needs a local-versus-background density test no guideline specifies. |
| PTCH1 | c.3472A>T | benign | LP | 15 hets; same shape as CHD7. |
| PTEN | intronic splice region | benign | LP | You cite the PTEN VCEP (PMID 30311380). PVS1 is firing on a non-canonical intronic variant. This one looks like ours to fix and we are treating it as a bug. |
| CYC1 | c.454-7_454-2dup | benign | LP | A duplication in the acceptor region. We already fixed the case where an insertion re-inserts the identical base (PTEN c.802-2dupA); this is the more general form and is not fixed. AF 2.8e-3 with 10 homozygotes. |
| FMR1 | c.1126-1G>A | benign | LP | The dead zone, above. |
| OCA2 | c.1080C>T | pathogenic | LB | Needs the splicing functional data as a `--functional-evidence` entry; the mechanism is in place and the file is yours to supply. |
| PRR12 | c.3505C>T | pathogenic | LB | You judge BP1 wrongly called. BP1 now requires that missense is not an established mechanism, and PRR12 has too few ClinVar missense entries to trip that test. |
| IQSEC1 | c.962G>A | pathogenic | LB | Same as PRR12. |

## What we are asking for

In order of how much each would move:

1. **Founder and hypomorphic alleles to add to the BA1 exception list.** The mechanism is now real and verified; the list is Ghosh's nine. Gene plus HGVS or coordinate, plus a line of reasoning.
2. **A BS1 answer for dominant, early-onset disorders with no published VCEP bar.** Still the single highest-value decision, and FMR1 above shows why. Per-disease is fine, and preferable.
3. **CBS: your call.** Your highlight and your note disagree, and three rows turn on it.
4. **The eight rows you marked wrong in both directions.** These are excluded from scoring rather than counted. A final call on them would give us a small expert test set independent of ClinVar, which is worth more to us than the ClinVar comparison.
