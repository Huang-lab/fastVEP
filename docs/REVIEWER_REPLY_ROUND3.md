# Reply to the round-2 review

Thank you - this was the most useful round yet.
Every one of your six points is now implemented, and each is measured against the full ClinVar 2-star-or-better set (673,660 variants) rather than argued about.

**Where your 122 rows stand now:** 90 of them are no longer discordant.
They are VUS, which is the honest automated answer when the deciding evidence is case-level, literature-level or gene-level.
32 remain, and they are the subject of the last section.

Overall, against the full ClinVar 2-star-or-better set: pathogenic recall went from 48 % to **65 %**, benign recall from 57 % to **72 %**, and same-direction agreement from 71 % to **78 %**.

The attached workbook is **your own annotated file** with four columns appended on the right (`fastvep_now_class`, `fastvep_now_criteria`, `changed_since_review`, `what_changed`).
Your notes, highlighting and row order are untouched, so you can read across each row you already worked on.

---

## 1. No reliable predictors for in-frame deletions or duplications, so BP4 should not be used for in-frame variants

Done, and you were right that it was firing where nothing calibrates it.

BP4 no longer fires for `inframe_deletion`, `inframe_insertion`, `stop_lost` or `protein_altering_variant` unless the variant also overlaps a splice region, in which case the SpliceAI evidence is about splicing and is calibrated (Walker 2023).
Where BP4 is withheld, the output carries `splice_skipped_reason` so a curator can see it was considered and declined rather than silently absent.

The underlying rule is now explicit: REVEL is applied to missense only (Pejaver 2022 calibrates it for missense), SpliceAI only where splicing is the question.
This resolved the in-frame rows you flagged - AFG2A, RSPH9, AMHR2, OCA2, CYP24A1 - plus HBA1 stop-lost.

## 2. Common Mendelian phenotypes, late onset, low penetrance and variable expressivity, where BS1 and BS2 may not apply

Partly done, and this is where the honest answer is "partly".

**What is done.** BS2 no longer uses a bare homozygote count. It computes a 95 % lower confidence bound on the frequency of individuals with no functional copy and compares it against a **maximum credible disease prevalence**, so the test scales with cohort size instead of with a raw count. We then set that prevalence bar by measurement rather than convention: sweeping it across the whole benchmark, it sits at **1 in 1,000**, the prevalence of the most common Mendelian disorders, precisely so that hearing loss, alpha-1 antitrypsin deficiency and familial Mediterranean fever are covered by the default rather than needing a per-gene exception.

We also honour ClinVar's own vocabulary: a variant labelled `Pathogenic, low penetrance` or `Established risk allele` at 2 stars or better is outside BS2's "full penetrance expected at an early age" precondition by definition, and BS1/BS2 are withheld with a stated reason. That catches SERPINA1 (1,236 homozygotes), F2, ALPL and MITF.

**What is not done, and needs you.** Onset, penetrance and expressivity are not derivable from any machine-readable source we have. A per-gene-disease table - prevalence, onset class, penetrance class, and any published VCEP BA1/BS1 threshold, keyed by gene and MONDO ID - is the missing input. Every hook that would consume it exists and is waiting.

One measurement worth your attention: **no setting of any BS2 threshold drives the residual errors to zero.** Sweeping the prevalence bar across four orders of magnitude, the false-benign count bottoms out at 38 and BS2 is only involved in 9 of them. That is the quantitative argument that this cannot be solved with a better global number, and needs the per-disease table.

## 3. BS2 applied to AR diseases on a single or a few homozygotes

Fixed, and it was the single largest source of false-benign calls.

Under the cohort-scaled test, one homozygote among gnomAD v4's ~730,000 individuals no longer fires BS2, while three among 2,000 still does - the evidence now scales with what was actually surveyed. BS2 firings on truth-pathogenic variants fell by 80 %.

Two related fixes came out of the same work:

- BS2 now counts **hemizygotes** on non-pseudoautosomal sex-chromosome sites, which it previously could not see at all. FMR1 (210 hemizygotes) and IDS (429) were being called Likely Pathogenic because gnomAD appeared never to have seen a null individual. A correction worth recording: `nhomalt_XY` is not the field to use - gnomAD calls XY samples haploid outside the PAR, so it is zero even at a site with 109,916 XY carriers. `AC_XY` is the hemizygote count.
- BS1/BS2 are withheld entirely for genes on the Mandelker 2016 homology list you cited (CYP21A2, STRC, HBA1/2, SMN1, PMS2, NEB, OTOA and others), where the frequency is a mismapping artefact rather than a measurement.

## 4. Founder variants may be over-represented in gnomAD in that population

Half done, and I want to be precise about which half.

BA1 and BS1 now test gnomAD's **filtering allele frequency** - the lower bound of the 95 % confidence interval on the frequency, maximised across genetic-ancestry groups (Whiffin 2017) - instead of a point estimate. This fixes the *small-sample* half of the problem: a frequency inflated by thin sampling in one group no longer reaches BA1 or BS1.

It does **not** fix the other half. A founder allele that genuinely is common in a well-sampled population has a high filtering allele frequency too. For those, only a curated exception keeps BA1 off. The Ghosh 2018 SVI exception list ships as a default, but ACADS in Ashkenazi, RSPH9 in Middle Eastern, C9 in East Asian and OCA2 in African American populations need either the per-gene-disease table or explicit exception entries. This is the same gap as point 2.

## 5. A common hypomorph in trans with a LoF allele causes disease; homozygous it does not. Common variants need examining against the gene's mutation spectrum

The mutation-spectrum half is done; the hypomorph half is a documented limitation.

**Mutation spectrum.** BP1 previously fired on a constraint heuristic alone - pLI high, missense-Z low - which is a statement about population tolerance, not about disease mechanism. It now requires positive evidence that missense is *not* an established mechanism: BP1 is blocked when the gene has three or more pathogenic or likely-pathogenic missense variants in ClinVar, and hard-blocked when PS1, PM5 or PP2 fired for the same variant. BP1 firings fell from 22,567 to 6,471, and BP1 misfires on truth-pathogenic variants fell by 99 %. This resolved every BP1 row you flagged, ENG included.

**Hypomorph-in-trans.** No machine-readable source exists for "this common allele is only pathogenic in trans with a null", and a hand-maintained list would rot silently. What we do instead: ClinVar's low-penetrance and risk-allele terms catch the subset ClinVar has already labelled, and `gene_overrides` provides a documented hook with the reason surfaced in the output, so a curator sees that frequency evidence was suppressed and why. TXNL4A and MMACHC still need that hook set by hand. If you can name the gene-variant pairs you consider established, we will ship them as defaults with your attribution.

## 6. Literature and functional evidence should supersede predictions; rule out functional data before using BP7

Done, and this is the change I think you will find most useful.

`--functional-evidence` takes a curated TSV - `chrom, pos, ref, alt, criterion (PS3|BS3), strength, pmid, note` - and PS3/BS3 fire from it at the strength **you** assign. Strength is a column rather than a constant because Brnich 2020 makes assay strength a judgement about the assay's validity, controls and dynamic range; a curator who has read the paper is the one who knows whether it supports Strong or only Supporting.

An entry does more than add a criterion. Per the SVI's ordering, functional evidence outranks computational prediction, so an entry **suppresses PP3, BP4 and BP7 for that variant**, each carrying `Superseded by functional evidence` in its summary. This runs symmetrically: an assay showing no damaging effect outranks a high REVEL exactly as a damaging assay outranks a low one.

Your OCA2 example, run end to end:

```
before   Benign   BA1 + BS2 + BP4 + BP7
after    VUS      PS3 (PMID 26637981) + BA1 + BS2; BP4 and BP7 superseded
```

fastVEP will not mine literature itself - a wrong PMID in a clinical report is worse than no evidence - but it will consume what you curate, which is how VCEP pipelines work.

**On predictors being trained on the same functional data:** you are right, and it is not detectable. Training-set membership is not published per variant for REVEL or SpliceAI, so we cannot flag a variant as circular. It is recorded as a known limitation, and the functional-evidence override is the practical mitigation: where you supply the experiment, the predictor is not consulted.

---

## Two things that came out of reproducing your cases, that you could not have seen from the spreadsheet

**Your "does dot mean del?" was a bug in our export, and a real annotation bug underneath.** `ALT=.` is ClinVar's reference-agreement record - an assertion that the sample matched the reference, carrying no alternate allele. We were annotating those as real variants and producing full Likely Pathogenic calls on them. They are now skipped (382 across the benchmark), and the review table carries an explicit `variant` column with ClinVar's own genomic HGVS so no row has to be decoded from `ref`/`alt` again.

**Two of your rows were consequence-calling bugs, not criteria bugs.** Your MUTYH and HPS4 "why is PVS1 called, this is synonymous" rows were right: multi-nucleotide substitutions spanning a codon boundary were being translated in the wrong frame, producing a spurious `stop_gained`. Also fixed: an insertion that re-inserts the identical base at a splice site (PTEN c.802-2dupA, BRIP1 c.2258-2dup) no longer counts as destroying the acceptor, and indels spanning a splice boundary now populate the intronic offset so PVS1's ±1/±2 gate engages.

## Where the remaining 32 rows sit, and what would move them

Of your 122, 90 are resolved. The 32 that remain fall into three groups, and only one of them is ours to fix.

One row moved since this reply was first drafted, and it is worth naming because it came from your own citation. **MSH2 p.Gly315Val** no longer earns PM1. You pointed at cspec GN137, which excludes PM1 for MSH2, and the underlying reason turned out to be general: Richards 2015 asks for a hot spot or critical domain "**without benign variation**", and fastVEP was only ever counting the pathogenic half, because our ClinVar protein index only ever recorded pathogenic missense. It now records both directions. That window has 3 pathogenic and 5 benign missense variants in ClinVar, so it is not a hot spot, and the call drops from Likely Pathogenic to VUS. TP53 p.Arg248 (23 pathogenic, 0 benign) is untouched, which is the control.

This does **not** resolve your CHD7 and PTCH1 rows, and it is worth being precise about why. Your objection there was gene-level - "CHD7 variants are spread out in the gene", "no clear hotspots for PTCH1" - and both windows genuinely have 3 pathogenic and 0 benign neighbours. Answering that needs a test of whether a local cluster is denser than the gene's own background, which no guideline text specifies. It stays open.

**Gene-disease validity (RYK, GIGYF2, EMG1, ORAI1, ARMC9, NPHP4, KCNJ16, ALMS1).** You judged several of these genes not disease-associated, or associated only as susceptibility loci. We implemented a gate for exactly this, reading ClinGen Gene-Disease Validity - and then found it could not act on your cases, because **ClinGen has not curated these genes at all**. The gate now fires only where ClinGen has curated a gene and classified every proposed relationship as Limited, Disputed, Refuted or No Known Disease Relationship, because an earlier version that blocked on *absence* from ClinGen cost 1,497 correct pathogenic calls in genes like SPAST, ABCB11 and LAMB3 that simply have not been curated yet. If you are willing to state these gene-level judgements, we will ship them.

**Dominant genes where the heterozygote count is too high for the disease (MTR 7,000 hets, KMT2C 1,721, TP63 57, CHD7 30, PTCH1 15).** Your notes give the same reasoning in every case: the carrier frequency exceeds what the disorder's prevalence allows.

We acted on this one. BS1 used a flat 1 %, which is a placeholder for "greater than expected for disorder" rather than a derivation of it, and far too permissive for a dominant early-onset condition. Sweeping it across the benchmark:

| BS1 bar | benign recall | false-pathogenic | false-benign |
|---|---:|---:|---:|
| 1 % (was) | 56.3 % | 4 | 1 |
| **0.5 % (now)** | **58.5 %** | 4 | 1 |
| 0.1 % | 61.9 % | 3 | 2 |
| 0.05 % | 63.1 % | 3 | 7 |
| 0.01 % | 65.9 % | 3 | 18 |

The default is now **0.5 %**. Confirmed on the full set that is +3.1 points of benign recall for 4 additional missed diagnoses out of 673,660 - about 1,200 to 1, which we think is right, but it is a trade rather than free, and we would rather say so.

**This is where your judgement would help most.** Going further is clearly possible - 0.1 % buys another 3.4 points and *reduces* false-pathogenic calls - but past 0.05 % the missed-diagnosis count triples, and a single global number cannot be right for both a dominant early-onset condition and a common recessive one. MTR and KMT2C are resolved at 0.1 %; CHD7 (30 hets) and PTCH1 (15 hets) are not resolved at any global threshold and need the per-disease figure.

**Rows where you judged both calls wrong (BRCA2, COG7, CYFIP2, EMG1, KIAA0586, KMT2C ×2, MUTYH).** These are excluded from the concordance denominator rather than scored, and reported separately. If you can give a final call on them, they become a small expert-adjudicated test set independent of ClinVar, which is more valuable to us than the ClinVar comparison itself.

## The one question we would most like answered

Not a list this time - a single decision.

**For a dominant, early-onset disorder with no published VCEP threshold, what heterozygote frequency in gnomAD should be enough to call BS1?**

This is the same gap written up as [the frequency dead zone](ACMG.md#the-frequency-dead-zone-between-pm2-and-bs1): between PM2's bar (4e-5) and BS1's (5e-3) there is a two-order-of-magnitude band where fastVEP offers no frequency evidence in either direction, and your CHD7 and PTCH1 rows sit in it. Your MTR, KMT2C, CHD7 and PTCH1 notes all turn on this, and a defensible default would resolve that whole group. If you would rather give it per disease than as a global number, that is the gene-disease table from point 2, and we will build to whatever shape you can produce.
