# Round 7 - note for the medical geneticist

This answers `discordant_64_round6.xlsx` as you returned it, with a verdict on every one of the 64 rows and a longer note on 18 of them.
The same workbook with fastVEP's current call appended to every row is `discordant_64_round7.xlsx`.
This note covers every row, including the ones that needed no change, so it can be read without the workbook beside it.
The measurements behind each statement are in [`RUN_VERSIONS.md`](../../RUN_VERSIONS.md), section "v36 to v37".

Two things about how to read the numbers.
"ClinVar" below means ClinVar's own 2-star-and-better classification, which is a reference and not ground truth.
Discordance with it is not itself a defect.
Where a row is unchanged below, the reason is stated, and where that reason is "we need a PMID from you", it says so.

## What changed

Ten of the 64 rows are no longer opposite-direction.
Over the whole set of 673,660 ClinVar 2-star+ variants, opposite-direction calls fall from 59 to 49, and all 59 were rows of your workbook.
Five more rows (LMX1B, NFKB2, AGRN, LZTR1, PTEN) were already not opposite-direction on the run this round started from.
Your workbook was built from an earlier run, before the 3'-shift work of the previous rounds.

| Row | Was | Now | Why |
|---|---|---|---|
| MSH6 `c.3801+3_3801+6dup` | Likely pathogenic (PVS1) | VUS | The record anchors on the last base of exon 8, so its consequence is `frameshift_variant`. HGVS places the same edit four bases into the intron, where it changes no exon base and leaves the donor `GT` intact, so PVS1 no longer applies. Your note: "VUS: Generally if rare intronic dup then the assumption is it doesn't do much, utmost a VUS". |
| MSH6 `c.4001+1_4001+25dup` | Likely pathogenic (PVS1) | VUS | Same: 25 intronic bases duplicated, exon and donor intact. Your note: "LB: intronic in MANE". |
| DSP `c.939+3_939+6dup` | Likely pathogenic (PVS1) | VUS | Same. Your note: "DSP(NM_004415.4):c.939+3_939+6dup - 3bp dup in splice region ,,likely not affect splicing as the same bp were duplicated…so LB". |
| FOXP2 `c.618_620del` | Likely pathogenic (PVS1) | VUS | The change is one unit of a `CAG` repeat that runs across the acceptor, written by the record as the intron's final `CAG`, which includes the acceptor `AG`, so it is `splice_acceptor_variant`. Moved to the 3' end of the repeat, which is where HGVS places it, it is a `CAG` deleted inside the exon, in frame, with the acceptor `AG` intact. Your note: "B; AD, Mane S: NM_014491.4: Q209del? p.(Gln209del)?, common". |
| SNRNP200 p.Ser1087Leu | Likely benign (BS2, BP4) | VUS | BS2 for a dominant gene now asks whether gnomAD holds more carriers than the disorder could account for, and not whether it holds five. Nine carriers among 1.4 million alleles do not clear it. A splice BP4 (SpliceAI 0) no longer counts against a missense REVEL calls damaging, because the protein and splicing effects are assessed separately and the more deleterious stands (Walker 2023, which we quote from the preprint because we could not check the published text). Your note: "only 9 people in gnomadv4- BS2 wrongly called; Bp4 why?". |
| G6PD `c.563C>T` and `c.202G>A` | Benign (BA1/BS1/BS2) | VUS | On the curated frequency-exception list, from your notes. |
| SERPINC1 `c.439A>G` | Likely benign (BS1, BP4) | VUS | Same, from your note "founder variant in the African population". |
| CDH23 p.Thr2454=, COL7A1 p.Pro1699= | Likely benign (BP4, BP7) | VUS | BP7 must not be applied to a synonymous change at an exon's first base or last three (Walker 2023). The criterion read that rule and the pipeline never told it where the exon edge was. |

For COL7A1 your call depends on the inheritance ("likely LP for AR; for AD it will be LB"), and we report one call per gene, so VUS is the call for neither.

Nothing else about the annotation changed.
Across all 673,660 variants no consequence, HGVS or transcript differs, and on the GIAB HG002 genome the annotation is byte-identical to what it was; only the criteria and the calls differ.

## What it cost

The changes that fixed your rows move other calls too, and not all of them toward ClinVar.
This is the account of both directions.

- **Agreement with ClinVar's Benign and Likely benign classifications falls.**
  Net, 805 same-direction calls are lost (Likely benign 66,046 to 65,461, Benign 111,965 to 111,745), and almost all of them are VUS now.
  Two changes supply most of them: BP7 no longer applies at an exon's edge, which removes it from 330 variants ClinVar calls Benign or Likely benign, and BS2 no longer fires on a few carriers, which costs 474 benign-direction calls.
  The aggregate still rises, because ClinVar-VUS variants gain 1,085 exact matches and Pathogenic and Likely pathogenic gain 64 same-direction calls.
- **Six ClinVar Pathogenic or Likely pathogenic variants fall to VUS** because PVS1 no longer applies to them: ALPL `c.69_74del`, KCNQ1 `c.478_480del` and RYR1 `c.958_960del` (in-frame deletions at an acceptor), and LDLR `c.190+2_190+3dup`, DYSF `c.1267_1276+4dup` and `c.1272_1276+3dup` (duplications at a donor).
  In each, the edit leaves the exon and the canonical dinucleotide in place.
  ClinVar's classification rests on evidence we cannot compute from the variant.
  **Do you agree that VUS is the right call for these, with the stand-down recorded in PVS1's details?**
- **32 variants ClinVar calls VUS become Likely pathogenic** because BS2 no longer stands against them.
  22 carry a PVS1 that BS2 had been cancelling, with PM2 beside it or on its own (SCN11A x9, INF2 and TNPO3 x3 each, UMOD and ELN x2 each, and one each in HSPB8, TERT and DNAJB6).
  That is the same kind of call as the rows you reviewed, a PVS1 reaching Likely pathogenic.
  The other 10 are missense variants on PP3 or PM1 (CASR x5, FLNB x3, NOTCH2 and DNM2).
  The sweep that chose the new bar recovered 70 ClinVar Pathogenic or Likely pathogenic calls at that price.

## Rows that did not change, and why

### Where your own call is not ClinVar's (not a fastVEP defect)

| Row | Your call | Note |
|---|---|---|
| ARMC9 `c.2205dup` | LP for the AR condition | ClinVar is Benign; fastVEP is Likely pathogenic. |
| ORAI1 `c.142dup` | fastVEP correct | AD and AR mechanisms differ; "LP seems plausible". |
| RYR1 `c.14510del` | LP for the AR condition | ClinVar's Likely benign is for MH, a missense gain-of-function disorder. |
| RYK `c.10dup`, EMG1 `c.126dup` | VUS | You rule neither ClinVar's Benign nor our Likely pathogenic. |
| COL7A1 | LP for AR, LB for AD | See above. |

### Needs evidence we do not mine: your PMIDs would settle these

Curated PS3 or BS3 already suppresses PP3, BP4 and BP7 by name, and takes a PMID per variant through `--functional-evidence`.
The workbook's PMID column came back blank, and we do not search PubMed for functional evidence, because a wrong PMID in a clinical report is worse than none.

Twenty rows are Likely benign on BP4 and BP7 (some with PM2_Supporting): OCA2 `c.1080C>T`, POLR3A `c.1909+22G>A`, PNKP `c.1386+49_1387-33del`, MYO15A p.Gly3287=, NTRK1 `c.851-33T>A`, NPHP1 `c.771+58C>T` and `c.771+124C>T`, PGM3 `c.-2-185C>T`, TH `c.90+9dup`, `c.90+13G>A` and `c.90+13G>T`, PAH `c.1199+17G>A`, BBS4 `c.332+8T>C`, SLC4A11 `c.44-69del`, RTEL1 `c.396-26C>T`, SLC39A4 `c.192+19G>A`, LAMA3 `c.1273+51C>T`, NCF4 `c.758+50C>T`, MOCS1 `c.124-122_124-121insA`, and DST `c.1779-222C>T`.
As we read Walker 2023, that is the guideline's own outcome for a variant outside the splice region with a SpliceAI at or below 0.1, so these are not misapplied criteria; the evidence that overturns them is case data and functional data, such as your "genetic evidence and experimental evidence for affecting splicing" for PNKP.

### A null variant at a frequency you think is too high for the disease

KDM2B, PHKB, EDN3, RAD21, RUNX1, FMR1, CAMK2B, TP63, TBC1D8B and BRCA1 `c.594-2A>C` are all Likely pathogenic on PVS1 with a gnomAD frequency you rule too high for the disease.
Your CAMK2B note puts the principle: "PVS1 application should be nuanced", and that for loss-of-function AD genes the variant has to be rare in population databases.
PIEZO1 is in the same group, but the gnomAD record there failed gnomAD's own quality control (FILTER=AS_VQSR), so the frequency criteria are withheld and the call is PVS1 alone.

A frequency rule on PVS1 would address them, and the data say what it would cost.
Among PVS1 variants in genes whose established ClinGen tags are autosomal dominant only, with a population maximum between 3e-4 and 3e-3, 43 are ClinVar Pathogenic or Likely pathogenic and 5 are Benign or Likely benign.
Above 3e-3 it is 8 against 17.
Above 3e-3, then, a rule would correct 17 calls and miss 8 diagnoses.
**Is that a trade you would make, and if so is the line 3e-3?**

### Judgement calls that need a decision from you

| Row | Question |
|---|---|
| CHD7, PTCH1, FBN1 (PM1) | You say "PM1 is wrongly called ..CHD7 variants are spreadout in the gene.." and "no clear hotspots for PTCH1". PM1 fires on a window with three pathogenic missense. A test of each window against the gene's own ClinVar density fixes all three, and on an earlier build it moved about 2,700 other calls, 864 toward ClinVar and 1,865 away, so it is not in this change. Do you want it available as an opt-in? |
| DST `c.1779-222C>T` | You note "whenever MANE AND MANE PLUS CLINICAL ARE AVAILABLE MANE PLUS CLINICAL TAKES PRECEDENCE". As a default it lowers same-direction agreement by 402 calls and resolves this row. Would you accept the Plus Clinical call reported beside the MANE Select one when they differ, rather than instead of it? |
| G6PD (whole gene) | We exempted the two alleles you listed. Is every G6PD variant exempt from BA1, BS1 and BS2? |
| IQSEC1, PRR12 (BP1) | BP1 fires on 6,471 variants and only 10 are ClinVar Pathogenic or Likely pathogenic; two of those 10 are these rows. Do you want BP1 off for a gene with any ClinVar pathogenic missense? |
| CBS, three 68-bp insertions | You rule for ClinVar. Nothing in the frequency or consequence rules reaches it; it would need a curated entry with a source. |
| PADI3 `c.1669C>T`, SLC12A3 `c.539C>A` | BS1 fires on a population maximum of 6.0e-3 and 8.4e-3. You ask for a deeper dive on PADI3 and rule BS1 incorrect on SLC12A3. Are these founder alleles for the exception list? We need the population and a PMID. |
| NOTCH2 `c.17_18del` | You see it fail QC in gnomAD v2. We load v4, where the record holds 6 alleles, its FILTER is AS_VQSR (it failed gnomAD's own quality control), and the frequency criteria are withheld by the quality gate, so the call is PVS1 alone. |
| PRODH `c.554G>A` | gnomAD says the site "is covered in fewer than 50% of individuals". Our record holds 26,422 alleles where a well-covered site holds about 1.46 million, about 1.8%, and our only coverage guard is a floor of 2,000 alleles. Do you want a coverage-fraction rule for PM2 and PVS1? |
| TSC2 `c.3850C>T` | You note "TSC2 variants in DBs can be somatic as well". No field we carry says so; it would need a curated entry per variant. |
| BS2 and reduced penetrance | BS2 reads disease names, not ClinGen's inheritance tag, so it reaches about one gene in twenty. We measured reading the tag and have not made the change. It gains about 2,800 same-direction calls, and it also calls TTR p.Val142Ile and VWF p.Arg854Gln Benign and moves four loss-of-function alleles (NEK1, DDX41, SRCAP, PROKR2) from Likely pathogenic to VUS, all six of them ClinVar Pathogenic or Likely pathogenic. BS2 needs full penetrance at an early age, and the tag only says dominant or recessive. ClinGen labels NEK1 a susceptibility and DDX41 a predisposition syndrome, TTR p.Val142Ile belongs to a disorder usually described as incompletely penetrant (we have not checked this against a source), and ClinVar lists type 2N among the conditions of VWF p.Arg854Gln, a subtype Orphanet lists as recessive, while ClinGen's two VWF entries are both tagged dominant. Is there a penetrance source you would trust, so BS2 could be limited to disorders with full penetrance at an early age? |

## One thing we could not do from the review

**RUNX1 `c.-59-1G>A` and RAD21 `c.-32-1G>A`.**
For RUNX1 you cite the MM-VCEP: the variant is a canonical splice site in the 5' UTR whose effect is not the same across isoforms, so "PVS1 would not apply per the MM-VCEP".
For RAD21 you rule "LB for AD disease; canonical splice site in intron 1; but very common in gnomad".
Both are gene-specific reasoning, and neither is a rule fastVEP can apply to a gene it has no entry for.
A blanket rule that PVS1 does not apply to an intron wholly in the 5' UTR would remove PVS1 from 45 variants in the ClinVar 2-star+ set, 26 of them ClinVar Pathogenic or Likely pathogenic (BRCA2 five times, PAX6, GJB2 and GAA twice each, MMUT among them).
It would need a per-gene list, and a source for each entry.
