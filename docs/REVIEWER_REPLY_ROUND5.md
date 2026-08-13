# Reply to the round-5 review: frequency evidence, hypomorphs, and functional data

Three attachments.

- **`strong_discordant_round5.xlsx`** - every current opposite-direction call, 59 rows, with your round-4 ruling carried across on the 13 you coloured.
- **`new_for_review_round5.xlsx`** - the 46 of those that carry no ruling from you yet. This is the ask.
- **`discordant_122_round5.xlsx`** - your own annotated file with the current call appended, so nothing you wrote has to be re-located.

Both new tables have six judgement columns, two of which exist because of your comments below: one for *why* a frequency criterion should not apply to a given disorder, and one for a PMID that settles the case.

One framing point before the substance, because it changes what we are optimising.
We agree that discordance with ClinVar is not itself a defect.
The target is a faithful application of the guideline, and where a faithful application disagrees with a ClinVar consensus we would rather be faithful and say so.
What we cannot accept is discordance we cannot explain, which is why every row in these tables carries the criteria that fired and the evidence behind them.

Your six points, in order.
Two were already implemented, in both cases because of your earlier rounds.
Two were partly implemented and are now closed.
Two are open and need you.

---

## 1. "No reliable predictors for in-frame del or dup, so BP4 should not be used for in-frame variants"

**Agreed and already in force.**
BP4 refuses in-frame insertions, in-frame deletions, stop-lost and protein-altering variants outright.
The reason recorded in the code is yours from round 2: an in-frame deletion with SpliceAI 0.00 was collecting BP4 Supporting and, with BS2, reaching Likely benign on AFG2A, RSPH9, AMHR2, OCA2 and CYP24A1.

Across the whole 673,660-variant benchmark, 6,277 variants have an in-frame consequence and **24 of them collect BP4**.
Every one of those 24 also overlaps a splice region, which is the single deliberate exception: where a variant sits at a splice site, SpliceAI is a calibrated predictor of the question actually being asked, and Walker 2023 supports using it.
A typical row is `PM2_Supporting & PM4 & BP4` on an `inframe_deletion & splice_region_variant`.

**We would like your ruling on that carve-out.**
If your view is that an in-frame indel gets no BP4 under any circumstances, including at a splice site, it is a one-line change and we will make it.

## 2. "Common Mendelian phenotypes where BS1 and BS2 may not apply; prevalence, late onset, low penetrance, variable expressivity"

**Agreed, partly built, and the remaining part is a question only you can answer.**

Four mechanisms exist today, and they cover different parts of what you are describing:

| Mechanism | State | What it handles |
|---|---|---|
| BS2 tested against a maximum credible prevalence, not a raw count | shipped, default `1e-3` | late onset, reduced penetrance, variable expressivity - see point 3 |
| ClinVar `low penetrance` / `risk allele` terms block all benign frequency evidence | shipped | variants the archive has already flagged |
| Per-gene BA1/BS1 bars from the ClinGen CSpec Registry | built, **not default** | 304 published bars across 117 genes |
| Per-allele exception list | shipped, now covers BA1, BS1 and BS2 | see point 5 |

What does not exist is a switch that says "for this gene-disease pair, frequency evidence does not apply at all".
We have deliberately not invented one, because the three examples you gave want three different answers and we would rather you choose than guess:

- **Hearing loss / GJB2.** The Hearing Loss VCEP has published bars for this gene. The right answer looks like a looser bar, not a disabled criterion.
- **Alpha-1 antitrypsin deficiency / SERPINA1.** The Z allele is common, penetrance is partial and environmental. The right answer looks like a per-allele exception, in the shape of point 5.
- **Familial Mediterranean fever / MEFV.** Ghosh 2018 already lists two MEFV alleles as BA1 exceptions on low-penetrance grounds. Whether the whole gene should be exempt is a judgement we do not have.

We also built the raw material for a general answer and it fell short in one specific place.
There is now a gene-disease attribute table covering 6,857 gene-disease pairs from Orphanet: inheritance for 96.2 %, age of onset for 91.4 %, a prevalence class for 84.8 %.
**Penetrance is populated for 0 %**, because no public structured source states it.
So an automated "is this disorder fully penetrant by early adulthood" test cannot be built from public data, and the Richards 2015 precondition on BS1 and BS2 is exactly that question.

**What we need:** the rule you would apply, at whatever granularity is natural to you. Per gene, per disorder, or a general principle we can encode. `BS1_BS2_NOT_APPLICABLE_WHY` in the attached tables is where individual cases can go.

## 3. "BS2 applied for AR diseases where gnomAD had a single or a few homozygotes"

**Fixed after your round-2 review, and worth showing the arithmetic because it answers the general case, not just the examples.**

BS2 no longer counts homozygotes. It asks whether there are *more individuals with no working copy of the gene than the disorder itself could account for*, by comparing the 95 % lower bound on their frequency against a maximum credible prevalence of 1 in 1,000.
At gnomAD v4's cohort size that means **777 homozygotes are required before BS2 can fire on a recessive disorder**. One or two cannot do it, at any cohort size, and there is a hard floor at two regardless.

The wording of the declined result is deliberately your argument:

```
23 homozygote(s) among 707,639 gnomAD individuals gives a 95% lower bound of 2.22e-5,
which does not exceed the 1e-3 prevalence bar; consistent with a late-onset,
reduced-penetrance or variably expressive disorder rather than tolerance
```

That is GAA `c.-32-13T>G` with 23 homozygotes, and CFTR `c.1210-11T>G` with 27 and SPTA1 `c.4339-99C>T` with 40 read the same way. All three decline BS2.

## 4. "Founder variants may have increased representation in gnomAD in that population; carefully apply BS1, BS2"

**Half addressed by construction, half only addressable by curation, and there is a wrinkle in gnomAD worth your attention.**

BA1 and BS1 do not test a point-estimate allele frequency. They test gnomAD's filtering allele frequency, the lower bound of the 95 % confidence interval maximised across ancestry groups (Whiffin 2017).
That removes the small-sample half of the founder problem: a frequency inflated by thin sampling in one group no longer reaches the bar.
It does nothing about the other half. A founder allele that really is common in a well-sampled population has a high filtering AF too, and only a curated exception keeps the benign criteria off it.

The wrinkle: gnomAD's `grpmax` **deliberately excludes the Finnish, Ashkenazi and "remaining" groups** because they are bottlenecked.
So a Finnish or Ashkenazi founder allele is invisible to BA1 and BS1 by construction, in the benign direction.
This is why Ghosh 2018 had to hand-list ACADS and BTD, both annotated in their table as "detected at >5 % MAF only in Finnish".
It also explains FMR1 `c.1126-1G>A` from your last round: its Finnish frequency is 1.06 %, which would settle it, and grpmax cannot see it.

**What we need:** founder alleles, by gene plus HGVS or coordinate, plus one line of reasoning. They ship as defaults with your attribution. The mechanism is verified end to end and currently holds twelve entries.

## 5. "AR diseases where a common hypomorphic variant in trans with a LoF variant causes disease (TAR syndrome and others)"

**This was the sharpest point and it is now implemented.**

You are describing a class the frequency criteria are structurally wrong about, not a threshold that needs tuning.
If a variant causes disease only in trans with a rarer null allele, then it being common is what the disease model *predicts*, and homozygotes being unaffected is also what the model predicts.
Reading either as evidence against pathogenicity inverts the mechanism.

The curated exception list, which previously covered BA1 only, now takes a `blocks` field naming which frequency criteria an entry suppresses.
Three hypomorphic alleles ship as defaults, all three drawn from the discordance table you have been reviewing:

| Variant | gnomAD AF | Hom | Was | Now |
|---|---:|---:|---|---|
| GAA `c.-32-13T>G` (late-onset Pompe) | 5.4e-3 | 23 | Likely benign | VUS |
| CFTR `c.1210-11T>G` (the 5T allele) | 9.8e-3 | 27 | Likely benign | VUS |
| SPTA1 `c.4339-99C>T` (alpha-LELY) | 6.6e-3 | 40 | Likely benign | VUS |

All three are ClinVar Pathogenic, all three were being called Likely benign on BS1, and none of them reaches 5 % - which is why a BA1-only exception list never touched them.
Note that none of the three reaches Likely pathogenic even now, and should not: the deciding evidence is PM3, the second allele in trans, which no annotator can compute from a single variant.
VUS is the honest automated answer for all three.

Your TAR example, RBM8A `c.-21G>A`, is not in this benchmark because it is not 2-star in ClinVar. It should be on the list and we have not added it without your confirmation of the allele.

**What we need:** the rest of the list. Every hypomorphic allele you would name gets an entry, and the entry says in prose why frequency does not apply to it.

## 6. "How can we incorporate literature-based information? BP7 was applied where experimental evidence for a splice defect is known"

**Already implemented, in direct response to your OCA2 case, and here is the demonstration.**

Curated functional evidence now outranks every computational prediction. When PS3 or BS3 is met, PP3, BP4 and BP7 are all suppressed with the reason named.
Your reasoning that the predictors were often *trained on* that same functional data is exactly right, and it is why suppression runs in both directions rather than letting the predictor cast a dissenting vote.

OCA2 `c.1080C>T`, your row, run twice on the shipping binary:

```
without functional evidence                     -> Likely benign
    BP4  met   SpliceAI max_ds=0.01
    BP7  met   Synonymous, no predicted splice impact, PhyloP 0.02

with one line of curated PS3                    -> VUS
    PS3  met   Curated functional evidence: assay shows a damaging effect at Strong
    BP4  not met   Superseded by functional evidence: PS3 rests on a curated assay,
                   and ClinGen SVI ranks experimental evidence above computational prediction
    BP7  not met   Superseded by functional evidence: (same)
```

The input is a four-column TSV keyed by coordinate, with an optional strength per entry because Brnich 2020 makes assay strength a judgement about the assay rather than a constant:

```
#chrom  pos       ref  alt  criterion  strength  pmid      note
15      27990612  G    A    PS3        Strong    <PMID>    minigene shows aberrant splicing
```

**On literature mining, our position is firm and we would like you to challenge it if you disagree.**
fastVEP will not search PubMed to populate PS3 or BS3.
A wrong PMID in a clinical report is worse than no evidence at all, and there is no way to validate an automated literature call at the scale this tool runs at.
What we will do is consume a curated file, which is how every VCEP pipeline works.

If you would endorse a public structured source - ClinGen's VCEP-approved functional assertions, or MaveDB for the multiplexed assays - we will ship a seed file built from it and cite it. That is a decision about clinical acceptability, so it is yours, not ours.

Note also what the OCA2 example shows about the ceiling: even with a sound Strong assay, the variant lands at VUS, because PS3 plus PM2_Supporting is 5 points and Likely pathogenic needs 6. The missing point is PM3 or PP4, both case-level. This is the same wall as the 55 rows below.

---

## What the two fixes cost and bought

Both changes above were measured on the full ClinVar 2-star-or-better set.

| | v18 (your last table) | v19 | **v20 (attached)** |
|---|---:|---:|---:|
| Exact match | 63.1 % | 63.1 % | **63.1 %** |
| Same-direction | 77.5 % | 77.5 % | **77.5 %** |
| Opposite-direction calls | 76 | 62 | **59** |
| Called pathogenic, ClinVar benign | 43 | 29 | **29** |
| Called benign, ClinVar pathogenic | 33 | 33 | **30** |

**v19** fixed a defect your comments led us to, which none of us had named.
When gnomAD rejects a site - its own FILTER, a low-complexity tract, a segmental duplication - fastVEP withheld BA1, BS1, BS2 *and* PM2, which is even-handed in intention and badly one-sided in effect.
The benign side loses up to Standalone plus two Strong criteria; the pathogenic side loses one Supporting; and PVS1, read from the very alignments just declared untrustworthy, is kept at full strength.
RAI1 `c.840del` - **48,739 gnomAD homozygotes**, ClinVar Benign - was reported Likely pathogenic on PVS1 alone, and seventeen more had the same shape.
Those calls are now Uncertain, with the veto named in the summary.
The rule fires only where the suppressed frequency would itself have met BA1 or BS1, which matters: without that condition it demoted 2,630 correct pathogenic calls to remove 18 wrong ones.
With it, 14 wrong calls go for 6 correct ones.

**v20** is point 5 above: three fewer false-benign calls, at zero cost to anything else.

We also fixed a smaller hole found on the way. BS1 stands down where BA1 already covers the frequency, but it was comparing against the global 5 % default rather than the gene's own published bar, so in a gene with a looser VCEP bar a common allele could collect no frequency evidence in either direction.

## Where we stand against your own rulings

Of your 122 rows: **93 are now VUS**, 20 Likely pathogenic, 9 Likely benign.
Scored against your colour coding on the 75 rows where you stated a call: 1 agree, **55 where we now say VUS and you committed**, 11 still opposite, and 8 where you rejected both calls.

The 55 is the real number and it has barely moved in three rounds.
Your notes on those rows cite segregation, phase and functional data. None of it is computable from a variant and an annotation database, and no amount of threshold tuning will produce it.
That is a statement about the ceiling of automated classification, not an excuse, and we would rather report it than tune toward your answers by other means.

The 11 still opposite are: PTEN (intronic PVS1, our bug, still open), OCA2 (needs the functional file from point 6), CBS x3 (68 bp insertions), CHD7 and PTCH1 (PM1 density), PRR12 and IQSEC1 (BP1), CYC1 (splice duplication), FMR1 (the frequency dead zone).

## What we are asking for

1. **The hypomorph list.** Point 5's mechanism is real and holds three entries, all found in your own table. Name the rest.
2. **The founder-allele list.** Same mechanism, same shape, and grpmax means curation is the only route.
3. **Your rule for point 2**, at whatever granularity is natural: which disorders take frequency evidence off the table entirely, and which just want a different bar.
4. **The 46 rows in `new_for_review_round5.xlsx`.** These have never been adjudicated. Even a direction and one line each would give us an expert test set independent of ClinVar, which is worth more to us than the ClinVar comparison.
5. **Two narrow rulings:** BP4 on an in-frame indel that overlaps a splice site (point 1), and whether a public functional-evidence source is clinically acceptable as a seed (point 6).
