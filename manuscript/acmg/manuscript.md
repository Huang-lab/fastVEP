# Iterative incorporation of ACMG/AMP and ClinGen SVI guidance into an open-source variant classifier: expert review, how thresholds were chosen, and discordance with ClinVar 2-star+ classifications

**Kuan-lin Huang**^1,*^

^1^ Department of Genetics and Genomic Sciences, Precision Omics Laboratory, Icahn School of Medicine at Mount Sinai, New York, NY, USA

\* Corresponding author. Website: https://PrecisionOmics.org

---

## Abstract

**Background.**
Software that applies the ACMG/AMP standards must adopt each ClinGen refinement and set the thresholds that no recommendation states.
We describe how the classification engine of fastVEP, an open-source reimplementation of the Ensembl Variant Effect Predictor, was built, reviewed by a medical geneticist and measured against ClinVar's classifications of review status two stars or better, which are a reference and not ground truth.

**Methods.**
Thirty-seven numbered versions were scored on 673,660 ClinVar variants (release of 4 April 2026; 31 have stored scored output) over six rounds of review, three of which returned rulings.
For this paper 89 further runs re-measured all but four groups of the thresholds that no guideline states and tested departures from the sources.
The final version was applied unchanged to ClinVar's release of 4 October 2026, including 107,528 variants that were not in the April set.

**Results.**
On the April set, where every threshold was tuned, agreement in direction rose from 51.0% (v1, reconstructed) to 77.6% and was flat after the thirteenth version, while opposite-direction calls went from 25 (v10) to 94 and back to 49 (59 under a strict leave-one-out reading, three of the added calls Benign).
None of the 4,017 Pathogenic or 67,700 Benign calls is among the 49, and rules with no combining rule or safeguard reached 84.6% agreement with 1,730 contradictions.
Review of the discordant calls first identified 30 of 64 recorded changes (by our coding), and the largest class (21) was a condition of a guideline that a criterion had not applied.
The ClinGen file that the benchmark read lacked 104 curations, among them CFTR, all of whose 558 null variants were VUS.
BP7's conservation test, which the splicing recommendation cautions against without an empirical justification, costs 16,482 benign-side calls that agree with ClinVar to avoid two that contradict it.
On the 107,528 new variants agreement in direction was 70.1% and contradictions were 2.3 times as frequent (95% interval 1.2 to 4.0).

**Conclusions.**
The largest class of change was a condition that the standards state and a first implementation omitted.
Agreement in direction does not separate cautious from aggressive calling, so contradictions belong beside it.
The remaining discordance is consistent with evidence that a single-sample VCF does not carry, with design choices such as approximating PVS1's mechanism condition, and with an input that was incomplete for a reason we did not establish.

**Keywords:** variant classification, ACMG/AMP, ClinGen SVI, ClinVar, variant interpretation, benchmarking, fastVEP

---

## 1. Introduction

The interpretation of germline sequence variants rests on the 2015 standards of the American College of Medical Genetics and Genomics and the Association for Molecular Pathology (ACMG/AMP), which define 28 criteria of pathogenic and benign evidence and the rules that combine them into five classes: Pathogenic, Likely pathogenic, Uncertain significance (VUS), Likely benign and Benign (Richards et al., 2015).
Applied by different laboratories, the same standards can give different answers for the same variant, and the Clinical Genome Resource (ClinGen) has worked to narrow that variation, through general recommendations from its Sequence Variant Interpretation (SVI) Working Group and through gene-specific specifications written by Variant Curation Expert Panels (Rehm et al., 2015; Harrison et al., 2017).

The general recommendations are many, and they change criteria that looked settled.
The loss-of-function criterion PVS1 became a decision tree with graded strengths (Abou Tayoun et al., 2018).
The stand-alone benign frequency criterion BA1 acquired a list of exceptions (Ghosh et al., 2018).
The Bayesian reading of the framework turned the combining table into a point system (Tavtigian et al., 2018; Tavtigian et al., 2020).
Absence from population databases (PM2) was reduced to Supporting strength (ClinGen SVI, 2020), and computational evidence (PP3 and BP4) was calibrated tool by tool (Pejaver et al., 2022).
Splicing evidence received its own recommendations (Walker et al., 2023), and the SVI proposed that laboratories stop using the reputable-source criteria PP5 and BP6 (Biesecker and Harrison, 2018).

Software that applies the framework at scale must commit to each of these decisions, and to many that no recommendation settles.
It has to say which numeric frequency counts as "extremely rare", how many carriers make a variant "observed in a healthy adult", and how far into an intron a variant must lie before splicing can be set aside.
Several tools automate part of the framework (Li and Wang, 2017; Nykamp et al., 2017; Kopanos et al., 2019).
A common way to evaluate one is to compare its output with ClinVar (Landrum et al., 2018), whose classifications are a widely shared reference.

That comparison has limits that are easy to state and easy to forget.
ClinVar's review status counts submitters and expert panels, not truth.
Its classifications can rest on the same population-frequency, computational and consequence evidence that the software reads, so agreement can be partly circular.
A single concordance figure cannot say whether a disagreement is a defect in the software, an out-of-date or mistaken ClinVar record, or evidence that no annotator can see.
And the failure that matters most in a clinical setting is not a crash but a plausible wrong answer, such as a pathogenic variant called benign, or a criterion that silently never fires.
Neither need show up in an agreement figure.

Here we describe how the ACMG/AMP classification engine of fastVEP, an open-source Rust reimplementation of the Ensembl Variant Effect Predictor (McLaren et al., 2016), was built, reviewed and measured.
The engine runs inside the annotation pass, so a benchmark of 673,660 ClinVar variants of review status two stars or better takes about two minutes to annotate (a median of 109 seconds over the runs of this paper, made up to three at a time on one Apple M4 machine with 10 cores and 24 GB) and about ten seconds to score.
That speed allowed 37 numbered benchmark runs in 24 weeks (21 April to 6 October 2026), and a medical geneticist reviewed the calls on which the software and ClinVar disagreed in direction over six numbered rounds, in three of which they returned an annotated workbook.

We report the following.
First, how each guideline was incorporated and how agreement with ClinVar moved from one numbered run to the next (a run bundles several changes, so a change in agreement cannot be assigned to one guideline).
Second, what expert review found that the benchmark could not, and how the changes fall into a small number of classes.
Third, how each numeric threshold that no guideline states was chosen and how much it matters in the final classifier, re-measured on the final code and not recalled from the version in which it was set.
Fourth, how much of the agreement depends on ClinVar-derived inputs.
Fifth, the discordance of the final classifier with ClinVar's 2-star-and-better classifications, on the benchmark snapshot, on a ClinVar release six months later, and on the variants that reached two stars after all tuning was done.
Along the way we read each recommendation that the classifier cites for a number or a grade against its source, and report where the implementation departs from it.

---

## 2. Methods

### 2.1 The classifier and its evidence sources

The classification engine is the `fastvep-classification` crate of fastVEP (version 0.4.0): 14,857 lines of Rust with 340 unit tests at the commit analysed.
It is called from the annotation loop after consequence prediction and HGVS naming, once per allele and transcript, and with `--pick`, which was used throughout, one transcript is reported per variant, chosen by VEP's default order (MANE Select, MANE Plus Clinical, canonical, APPRIS, TSL, biotype, CCDS, consequence severity; MANE is described by Morales et al., 2022).
The Ensembl GFF3 that we used holds no APPRIS attribute, so that step of the order does nothing here.
Both per-variant loops of the code base, the shared annotation library and the command-line batch path, call the same classifier and the same helper for splice context, and each has end-to-end tests.

The engine implements all 28 ACMG/AMP criteria, and a run reports 26 of them, because PP5 and BP6 are off by default, following the SVI's proposal that laboratories stop using them (Biesecker and Harrison, 2018).
In the benchmark, which has one variant per record and no pedigree or phenotype, 15 of the 26 can fire: PVS1, PS1, PM1, PM2, PM4, PM5, PP2 and PP3 on the pathogenic side, and BA1, BS1, BS2, BP1, BP3, BP4 and BP7 on the benign side.
PM3 needs the genotype of a proband, which the sites-only benchmark file lacks, so it is never met (it is reported as evaluated in genes that are not read as recessive and as not evaluated in genes that are), and the other 10 need case-level data and are reported as not evaluated: PS2, PS3, PS4, PM6, PP1, PP4, BS3, BS4, BP2 and BP5.
Supplementary Table S12 gives how often each computable criterion fires in v37, at which strengths, and on which ClinVar classes.
Two gene-level gates are referred to in the text: BP1 is withheld (the veto) in a gene with at least three pathogenic missense variants in the ClinVar protein index, where missense is an established mechanism, and PVS1 is withheld by a validity gate and a mechanism gate, the validity gate also withholding PP2 and PM1 (Table 1).
PS3 and BS3 are read from a curated file when one is supplied, PS2 and PM6 from the genotypes of a proband's parents, and PM3 and BP2 from the proband's genotype together with its other variants in the gene; none of these routes is used in the benchmark.
Strength follows the guidelines where they grade a criterion: PVS1 by a decision tree after that of Abou Tayoun et al. (2018), whose departures from the published tree are listed in Section 3.2, PP3 and BP4 by the REVEL bands of Pejaver et al. (2022), and PM2 at Supporting (ClinGen SVI, 2020).
Table 1 lists the guidance incorporated, when it first landed, and what the classifier does with it.

**Table 1. Guidance incorporated into the classifier.**
"Written" is the date of the commit, or, for items built in the August sprint that was merged as one change, of the benchmark run that first contained them; "live" is the run from which the criterion influenced a call, which for several items is much later than the date it was written (Section 3.2).
A numbered run bundles every change made since the one before it, so the change in agreement between two runs cannot be assigned to one row.

| Source | Criteria | What the classifier does | Written (2026) | Live |
|---|---|---|---|---|
| ACMG/AMP standards (Richards et al., 2015) | all 28 | Evaluates every criterion; five classes; the 2015 combining table is kept on the benign side | 21 Apr | v1 |
| Bayesian framework and point system (Tavtigian et al., 2018, 2020) | combining rules | Pathogenic side scored 1, 2, 4 and 8 points for Supporting to Very Strong, 10 or more Pathogenic, 6 to 9 Likely pathogenic; a Benign Very Strong counts as two Strong | 12 Aug | v11 |
| PVS1 recommendations (Abou Tayoun et al., 2018) | PVS1 | A decision tree after the published one, over the prediction of nonsense-mediated decay (NMD), last exon, fraction of protein removed, critical region and start loss; it departs from the tree in two grades and in how it predicts NMD (Section 3.2) | 29 Apr | v8 (graded), v13 (inputs derived) |
| BA1 exception list (Ghosh et al., 2018) | BA1 | Nine alleles exempt from BA1, matched by genomic coordinate first and HGVS second; extended by three hypomorphic and three reviewer-named alleles, which also block BS1 and BS2 | 29 Apr | v18 |
| PM2 recommendation, version 1.0 (ClinGen SVI, 2020) | PM2, combining | PM2 at Supporting; one Very Strong plus one Supporting reaches Likely pathogenic. The recommendation gives no frequency bar, and the bars are not from it (Section 2.5) | 29 Apr | v6 (absence read as rarity) |
| PM3 recommendation, version 1.0 (ClinGen SVI, 2019) | PM3 | Points-based in-trans scoring; needs the genotype of a proband, which the benchmark lacks | 29 Apr | n/a |
| Calibration of in-silico tools (Pejaver et al., 2022) | PP3, BP4 | REVEL at 0.644, 0.773 and 0.932 for PP3 Supporting, Moderate and Strong; at 0.290, 0.183, 0.016 and 0.003 for BP4 Supporting to Very Strong; missense only; PP3 is dropped beside PS1 or PM5, and PM1 beside a Strong PP3 or beside PVS1 | 29 Apr | v1 |
| Splicing recommendations (Walker et al., 2023) | PP3, BP4, BP7, PS1, PVS1 | SpliceAI at least 0.2 for PP3 and at most 0.1 for BP4, both at Supporting strength; BP7 not for a synonymous change at an exon's first base or last three, and for an intronic change only at +7 or beyond on a donor side and -21 or beyond on an acceptor side, but with the 2015 conservation test kept and SpliceAI below 0.2, where the recommendation applies BP7 after BP4 (SpliceAI at most 0.1) and cautions against the conservation test (Sections 3.2 and 3.3); PS1 for a variant on the same canonical splice dinucleotide as a Pathogenic one; splice BP4 set aside beside a REVEL PP3 | 29 Apr | v4 (PP3, BP4, BP7), v16 (PS1), v37 (BP7 exon edge) |
| Filtering allele frequency (Whiffin et al., 2017; ClinGen VCEP Review Committee, 2024) | BA1, BS1 | Tested against gnomAD's grpmax filtering allele frequency (the highest, across ancestry groups, of the lower bound of the 95% interval of each group's frequency) rather than a point estimate | 12 Aug | v10 |
| Homologous genes (Mandelker et al., 2016) | BA1, BS1, BS2, PM2 | Not evaluated in 29 genes whose population frequencies are unreliable | 12 Aug | v9 |
| Functional evidence (Brnich et al., 2019) | PS3, BS3 | Curated file with a strength column; supersedes PP3, BP4 and BP7 for that variant | 12 Aug | available from v10; not used in the benchmark |
| ClinGen gene-disease validity (Strande et al., 2017; ClinGen, 2026) and curated disease mechanisms | PVS1, PP2, PM1 | Withheld where every ClinGen curation of the gene is Limited, Disputed, Refuted or No Known Disease Relationship (the validity gate), and where the curated mechanism excludes loss of function (the mechanism gate, which applies to PVS1 only) | 12 Aug | v10 (a first version blocked on absence), v11 (curated-negative only) |
| Published VCEP bars (ClinGen CSpec Registry) | BA1, BS1, PM2 | 286 bars across 117 genes, each with its source sentence kept; of 411 criterion entries read, 18 are withheld where two panels disagree and 107 were rejected, among them 37 that require absence, 27 conditioned on inheritance and 21 with no threshold (Section 3.3); opt-in, not in the benchmark | 13 Aug | opt-in |
| Genes with multiple disorders (ClinGen, 2024) | all | Per-disorder override scaffold only; not active | 29 Apr | n/a |

The criteria read nine evidence sources (Table 2).
All are public.
Two of them (ClinVar and ClinGen) come from the community that produces the reference, and we did not check how far the training data of REVEL and SpliceAI overlap the reference, so the agreement is partly circular in ways that Section 3.4 can measure only for the ClinVar-derived inputs.

**Table 2. Evidence sources read by the classifier in the benchmark.**

| Source | Version | Criteria fed |
|---|---|---|
| gnomAD exomes (gnomAD consortium, n.d.) | v4.1: about 730,000 exomes, allele number about 1.46 million at a well-covered autosomal site; allele and carrier counts, homozygotes, hemizygotes (AC_XY), per-ancestry frequency, filtering allele frequency, FILTER, segmental-duplication and low-complexity flags | PM2, BA1, BS1, BS2 |
| gnomAD gene constraint (Karczewski et al., 2020) | v4.1 metrics, 18,173 genes in the database built from them: PVS1 takes a gene as intolerant of loss of function at pLI of 0.9 or more or LOEUF of 0.35 or less, PP2 asks for a missense z above 3.09, and BP1 for pLI of 0.9 or more with a missense z below 2 | PVS1, PP2, BP1 |
| REVEL (Ioannidis et al., 2016) | v1.3 | PP3, BP4 (missense) |
| SpliceAI (Jaganathan et al., 2019) | as distributed in the gnomAD v4.1 site VCF | PP3, BP4 (splicing), BP7, and the check that a canonical splice call of PVS1 must pass |
| PhyloP (Pollard et al., 2010) | as distributed in the gnomAD v4.1 site VCF | BP7 |
| ClinVar (Landrum et al., 2018), allele level | release of 27 June 2026 | BS1 and BS2 (a gate for low-penetrance and risk alleles); PM2 (a backstop that reads ClinVar's own population frequencies where gnomAD has no record); PS1 and PM1 (leaving the variant's own record out of the counts) |
| ClinVar, protein and splice index | `variant_summary` downloaded 29 April 2026 | PS1, PM1, PM5, BP1, PVS1 critical region |
| ClinGen gene-disease validity (Strande et al., 2017) | file created 30 April 2026: 2,924 genes, 3,546 curations, 2,739 of them Definitive, Strong or Moderate and the rest Limited, Disputed, Refuted or No Known Disease Relationship; it lacked 104 curations that ClinGen dates before that day (Section 3.5) | PVS1, PP2, PM1 gate; inheritance for PM2, BS2 and PM3 |
| RepeatMasker (Smit et al., n.d.) | UCSC track | BP3 |
| Ensembl (McLaren et al., 2016) | release 115 GFF3 and GRCh38 primary-assembly FASTA | consequence, HGVS, PVS1 tree inputs |

### 2.2 Reference set, scoring and definitions

The reference is ClinVar's VCF for GRCh38 dated 4 April 2026.
We kept records whose review status is at least two stars (criteria provided, multiple submitters, no conflicts; reviewed by expert panel; practice guideline) and that are single-nucleotide variants or indels of at most 100 bases.
A record's class is read from its CLNSIG by substring match in a fixed order (Pathogenic/Likely pathogenic, then Benign/Likely benign, then Likely pathogenic, Likely benign, a leading Pathogenic or Benign, uncertain significance, and conflicting classifications), and a record that matches none, a drug response for example, is dropped.
That gives 673,660 records: 79,823 Pathogenic, 13,989 Likely pathogenic, 295,298 VUS, 128,038 Likely benign and 156,512 Benign, of which 652,400 are two-star, 21,237 are three-star and 23 are four-star.
ClinVar's combined labels are mapped to the first class (Pathogenic/Likely pathogenic to Pathogenic, Benign/Likely benign to Benign), so exact agreement is a strict measure and agreement in direction is the fair one.
The VUS class is the 295,298 variants that submitters classify as uncertain without conflict, except that two RYR1 records of expert-panel review whose classification is conflicting are labelled VUS by the label script, and fastVEP calls both Likely pathogenic.
Four records of CFH and CFI whose classifications are Likely pathogenic and Likely pathogenic with low penetrance are labelled Pathogenic by a pattern match of the same script.
We kept both labels so that the 37 runs stay comparable: relabelling the four changes exact agreement by one variant and agreement in direction by none, and leaving out the two RYR1 records lowers the number of ClinVar VUS called pathogenic-side by two.
A variant's gene is the first symbol of the GENEINFO field of its ClinVar record, and the population maximum is the highest allele frequency among the ten population groups that the annotation carries, "other" and "remaining" included, which is what the frequency criteria fall back on where there is no filtering allele frequency.
We call these classifications a reference and not ground truth throughout.

Pathogenic-side means Pathogenic or Likely pathogenic and benign-side means Likely benign or Benign (P/LP and B/LB in the tables), and the recall of a ClinVar group is the share of its variants that the classifier puts in the same direction.
A committed call is a call of one of those four classes, and a quiet variant is a synonymous, intronic (polypyrimidine-tract variants included), untranslated, flanking (upstream or downstream of a gene) or non-coding one.
"Missed diagnosis" and "false-pathogenic" name the direction in which a call contradicts ClinVar, and "correct" means in agreement with it; none of them says that ClinVar is right, and the reviewer supported ClinVar's class for 13 of the 24 missed diagnoses and 17 of the 25 false-pathogenic calls of the final classifier (Section 3.6).

The scorer takes the call on the picked transcript.
A call is exact if it equals ClinVar's class; same-direction if both are pathogenic-side (Pathogenic or Likely pathogenic), both benign-side, or both VUS; and opposite if ClinVar is pathogenic-side and the call benign-side, or the reverse.
An opposite call is a false-benign call when ClinVar is pathogenic-side (a missed diagnosis) and a false-pathogenic call when ClinVar is benign-side.
Rates use all 673,660 records as the denominator, which is the smaller figure: the 295,298 ClinVar VUS cannot be contradicted, so the 49 opposite-direction calls of the final classifier are 0.0073% of the records, 0.013% of the 378,362 with a pathogenic-side or benign-side class and 0.020% of the 244,093 committed calls.
Among them, 382 ClinVar reference-agreement records (ALT is `.`) produce no call and count as not concordant.

### 2.3 Iterative development

The unit of the process is a run: the full benchmark annotated with `fastvep annotate --acmg --pick` and scored by one script.
After every change to criteria, thresholds, data or annotation, the benchmark was re-run and compared with the run before it.
We refer to each numbered run as v1 to v37.
From v37 on, the harness records what produced an annotation (flags, the hash of the binary and of every configuration, which inputs were read, and the commit) and refuses to reuse an annotation that something else produced.
From v9 on, runs are scored in a leave-one-out mode in which a variant's own ClinVar record is excluded from the evidence that PS1 and PM1 read; a run that lets the record count is called ClinVar-informed.
PM5 cannot match its own record, because it asks for a different amino acid change at the same residue, and the critical region of PVS1 reads residues downstream of the variant.
The veto of BP1 and the allele-level ClinVar gates read the record as it is (Section 3.4).
Where a change touched the annotation layer (consequence or HGVS), every variant whose call moved was re-annotated with Ensembl VEP release 115.1 on the same GFF3 and FASTA with `--pick`, after confirming that both tools reported the same transcript, to establish which tool the change brought the software closer to; the written record keeps the outcome only as ledger items (3 were found by this comparison), and we do not report the comparison itself.

Not every numbered run is a step on the way.
Run v3 was a partial run that was overwritten before it was preserved, v21 and v22 repeat v20, v24 and v25 were intermediate runs that are not shipped states, and v17 is the opt-in run with the published VCEP bars, which v18 does not carry.
We mark v17 in the figures and leave it out of every range that the text quotes.
The raw output of v1 to v6 was overwritten, so v2 and v4 to v6 are quoted from the figures recorded at the time, and v1 from a matrix that was reconstructed from the documentation of the run: scored matrices are stored for 31 of the 37 runs, v7 to v37.
Every rate is over the 673,660 records of the benchmark; v1 annotated only 627,376 of them, and its recorded 54.7% is over those, where over all records it is 51.0%.
We could not check the denominators of v2 to v6.

We also kept a ledger of the changes that the project's written record (RUN_VERSIONS.md and the documents of the review rounds) describes as changing a criterion, a threshold, an input or the scoring: 64 items, each coded for its class and for the route by which it was first found.
The coding was done by one coder from the written record, was not repeated independently, and is not a record of every commit.

### 2.4 Expert review

A medical geneticist reviewed the calls on which the software and ClinVar disagreed in direction.
Each row of the first workbook carried 75 columns of evidence, and later workbooks carried more: ClinVar's record, the criteria that fired, REVEL, SpliceAI, PhyloP, gnomAD counts and frequencies by ancestry, gene constraint and the gene's ClinGen disease annotations.
The reviewer saw only the opposite-direction calls, with ClinVar's class and ours both visible, and not the pathogenic-side variants that we left as VUS, nor any concordant call.
They are a medical geneticist [qualifications and affiliation to be added by the author, with the reviewer's consent]; they were not blind to the software or to its calls, and they supplied two inputs to the classifier (the table of Ghosh et al. and the names of three alleles; ledger items 50 and 63), so they were a contributor as well as an evaluator, and their rulings are one expert's reading of rows on which both calls were visible and not independent evidence of which call is right.
The reviewer returned annotated workbooks in rounds 1, 2 and 6, with the verdict recorded as a cell colour in rounds 1 and 2 and as text in round 6, and a covering email in round 2; rounds 3 to 5 were tables and notes we sent after re-runs, which were not returned marked.
The workbook of round 1 carries a list headed "Discordance categories" that names six kinds of discordance, and a key of two colours (a correct call, and both calls discordant with the reviewer's call); the key of round 2 adds a query, a request for more evidence, and a marker for genes whose population frequencies cannot be relied on because of homology.
We decoded the colours against each workbook's key, and read any green on a call as that call being correct (in round 2 the green on the cells is a different shade from the key's, and in round 1 three greens on fastVEP's call are darker than the key's).
In round 6 the workbook asked for the reviewer's class and for who is right in two columns, which came back empty, and the reviewer wrote their ruling as free text in the verdict column of all 64 rows and added a note to 18.
We mapped the rulings of all three rounds onto five categories (ClinVar's call supported, fastVEP's call supported, neither, needs more evidence or depends on the disorder, and no verdict).
The round-6 text was coded by us from its content by five rules given in the supplement, and the reviewer did not assign these categories.
Between rounds the reviewer's notes were carried onto the rows they belong to, and each round was answered with a note that said what had changed, what it cost, and what was measured and not done.

### 2.5 How thresholds were chosen

Thresholds fall into three kinds.
Some are fixed by a guideline and were not tuned: the REVEL bands, the SpliceAI cut-offs, BA1 at 5%, and the points that map to each class.
Some are not stated by any guideline: the PM2 bars, the BS1 bar, the BS2 bars, BP7's far boundary into the intron, BP7's PhyloP bar and the count of pathogenic missense variants at which BP1's veto applies.
Most of these were set by sweeping one parameter over the benchmark with everything else at its default, and three (the recessive-only PM2 bar, BP7's PhyloP bar and the count of BP1's veto) were set without a sweep and were swept for the first time for this paper.
A few are structural, such as which combining scheme is used on which side, and were compared directly.
The rules applied when a sweep did not decide, written down as they were used, were these.
Prefer a value inside the range that published VCEP specifications use; we measured that range from the extracted CSpec table (Section 3.3).
Prefer the safer side of a trade, because a missed diagnosis costs more than a variant left as VUS.
Take the knee of a curve where there is one, which we took to be the point past which the exchange rate (the extra correct calls per extra opposite-direction call) falls from hundreds of calls per added error to about a hundred or fewer.
Where a reviewed variant showed a rule to be wrong, fix the rule; where it showed a variant to be an exception to a correct rule, put it on the curated exception list; and never move a threshold to repair one variant.
Both the exception list and the dominant-gene BS2 rule were written after the reviewer's rows, so the 49 residual calls of Section 3.6 are in-sample: without the exception list there are 55, and there are 59 when the exception list, BP1's own-record count and ClinVar's low-penetrance gate are removed together (Section 3.4).
Four more groups of thresholds that are ours (seven numbers) were not varied, because we made no sweep of them: the pLI and LOEUF bars of PVS1 (0.9 and 0.35), the missense z-score bars of PP2 and BP1 (3.09 and 2), PM1's window and count of pathogenic neighbours, and the floor of homozygotes for BS2.

### 2.6 Re-measuring on the final classifier

The sweeps behind the shipped values were run on earlier versions of the code, in some cases on a sample, and the combining rule changed after several of them.
For this paper we re-ran every sweep of a threshold that no guideline states on the final classifier, as full-benchmark runs that differ from v37 in one setting: 8 values of the PM2 bar for genes not read as recessive-only, 6 of the PM2 bar for genes read as recessive-only, 6 of the BS1 bar, 4 of the BS2 homozygote bar, 7 of BP7's far boundary, 5 of BP7's PhyloP bar and 5 of the count of BP1's veto (41 runs).
Four more values of the PM2 bar were run under the 2015 combining table, to see whether the bar's importance belongs to the bar or to the rule (the interaction of Section 3.3), and we used the 8 runs of the dominant-gene BS2 settings (5 of the prevalence bar and 3 of the allele-count floor) that had already been made on this code.
To measure what each design decision buys, we also ran 16 single-switch ablations and 9 leave-one-source-out ablations, in which one source is removed from the annotation directory so that the criteria that read it lose their input (for gnomAD the run also tells PM2 not to read an absent record as rarity, because PM2 would otherwise be met for every variant).
The two point cutoffs of the pathogenic side, Likely pathogenic from 6 points and Pathogenic from 10, are constants in the code and not settings, so for 9 further runs we built the binary with a patch that changes eight lines (six added, two removed) and reads them from the environment (`manuscript/acmg/replay_validation/combiner_cutoff_env.patch`).
Two departures from the sources that we found while reading them again were measured with a second measurement build, made with a patch that exposes through the environment BP1's veto counting a variant's own ClinVar record and BP7's SpliceAI bar (`bp_hooks_env.patch`): a control, BP1 without the own record, BP7 with SpliceAI at most 0.1, and BP7 with that and no conservation test (4 runs); a fifth run with the same build removes together the three adjustments that read the variant's own ClinVar record or followed the reviewer's rows (BP1's veto counting that record, ClinVar's low-penetrance gate and the curated frequency-exception list), each of which had been measured alone; two more runs with that build give BP1's veto a count of 1 and a count of 2 and leave the variant's own record out, the pair that Section 4.5 proposes together.
Three runs change an input to see what the gap in PVS1 of Section 3.5 is made of: the ClinGen file with the curations that its export lacked put back, PVS1 licensed for genes that Orphanet labels as a loss-of-function disease, and both (3 runs).
With nothing set, each measurement build is the shipped classifier: the control of each reproduces v37 in every cell of the concordance matrix and in the firing count of every criterion.
All 89 runs (Supplementary Table S1) were made by the same harness: 73 with the binary of v37 (SHA-256 `d12c0d8e...`), 9 with the build of the first patch and 7 with the build of the second.
The binary of v37 was built at commit 5d2521f; the commits above it up to the head of that branch add tests, comments and two guards that no benchmark input reaches, and a binary built from the head produced an annotation byte-identical to v37's (RUN_VERSIONS.md).
The provenance file of each run is in the repository, and records the flags, the hash of the binary and of every configuration, the inputs read (by name and size, not content), and the commit and number of uncommitted files of the checkout that ran the harness.
The harness refuses to reuse an annotation unless the binary, the flags, the configuration and the inputs are those that made it; the environment variable that sets a run of a patched build is not part of that check and is recorded in the run matrix, not in the provenance file.

### 2.7 Offline replay of the combining rules

The combining step is a pure function of the criteria a variant met, and the output records them.
We therefore re-derived calls under alternative schemes, and under alternative band edges, from the recorded criteria without re-annotating.
The replay of the shipped scheme reproduces v37 for 673,641 of 673,660 variants.
The remaining 19 are variants that the replay calls Likely pathogenic and the classifier holds at VUS by its site-level frequency veto: when gnomAD rejects a site (by its own FILTER, a low-complexity tract or a segmental duplication) and the frequency that was withheld would itself have met BA1 or BS1, a Likely pathogenic call becomes VUS.
The output does not record the detail the veto needs, so the replay holds those 19 at VUS wherever it would otherwise call them Likely pathogenic.
The replay was checked against real runs of the same configurations (Supplementary Table S4): the Richards 2015 table on both sides, and the point cutoffs of Table 4.
It reproduces the table run, the control and 6 of the 8 runs with a changed cutoff cell for cell, and the other two (Likely pathogenic from 4 and from 5 points) within 3 and 1 variants, each a variant that the replay calls Likely pathogenic and the run, through the veto, does not.
We use the real runs wherever one exists, and the replay only for the benign-side point scales, which have no switch in the shipped code and for which we made no real run (Section 3.3).

### 2.8 A later ClinVar release

To test the final classifier on variants that none of its tuning used, we applied it unchanged to ClinVar's release of 4 October 2026, extracted the same way (781,135 records at two stars or better).
The gnomAD databases the classifier reads were extracted for the regions of the April set, and a variant outside them has no gnomAD record, which the PM2 rule would read as absence.
We therefore kept the variants whose position lies inside those regions, the rule that decided which April variants had gnomAD data, applied to every variant alike, and we report how many it drops.
Variants are assigned to four strata by comparison with the April set: new (not in it), unchanged, class changed within its direction, and class changed across a direction.
A new variant is not one that ClinVar did not hold in April: 98.5% of them have a record in the `variant_summary` of 29 April from which the protein index was built, and for 92% of them it is a single-submitter record (Section 3.7).
The extract keeps only positions inside the April regions, so the new variants lie near April variants.
The ClinVar-derived indices the classifier reads were not rebuilt.
Seven runs of the changes that Section 4.5 proposes were made on the same variants, each differing from the unchanged classifier in one change: the completed ClinGen file, Orphanet's loss-of-function label, both, BP7 as Walker et al. describe it, BP1's veto without the variant's own record, BP1's veto at a count of 1, and no allele-number floor.
They were found from the April set and by reading the sources (Section 3.2), and the new variants were not scored by any sweep, review round or run, so they are the nearest thing to a test of them that the data allow.
To ask whether the higher rate of opposite-direction calls on the new variants is a difference in kind or in mix, we compared three things: the ratio of the rate on the new variants to the rate on the April set, with an interval from a bootstrap that resamples genes jointly across the two sets; the number of opposite-direction calls that the new variants would have at the April rate of their own stratum (ClinVar's direction or class, each alone or by consequence group), against the number observed; and the recall of each ClinVar class on the new variants against what the April recall of each band of population-maximum allele frequency would give them at their own mix of bands.

### 2.9 Statistics and software

Variants within a gene are not independent, so confidence intervals for the headline rates come from a gene-cluster bootstrap with 2,000 resamples (95% percentile interval), and proportions inside one table cell from the exact binomial (Clopper-Pearson) interval (Clopper and Pearson, 1934), which treats variants as independent.
The study is descriptive and we report no hypothesis tests; the interval of the ratio of the rates in Section 3.7 comes from a bootstrap that resamples genes jointly across the two sets, and so allows for the uncertainty of the April rate, which a Poisson comparison with it would not.
Analysis used Python 3.9.6 with pandas 2.3.3, NumPy 2.0.2, SciPy 1.13.1 and matplotlib 3.9.4; the binary was built with Rust 1.98.0.
Every number in the text and every table and figure is produced by `analysis/acmg_benchmark/scripts/11_manuscript_runs.py`, `12_per_variant_table.py`, `13_manuscript_stats.py`, `14_temporal_holdout.py` and `15_clingen_complete_export.py`, and by `manuscript/acmg/generate_figures.py`, from the scored runs.

## 3. Results

### 3.1 Thirty-seven runs: where agreement came from, and where it stopped

Same-direction agreement with ClinVar rose from 51.0% in v1 (a matrix reconstructed from the run's documentation) to 77.6% in v37 (Figure 1A; Supplementary Table S10 gives every run).
Almost all of that rise came in the first 13 runs, and the largest steps were mostly plumbing and not new guidance.
Agreement rose 8.8 points between v2 and v4, when SpliceAI and PhyloP results finally reached the criteria that read them (items 6 and 7 of the ledger of Section 3.2), and 5.2 points between v5 and v6, when a variant absent from gnomAD began to count as absent for PM2 and BP4's splice path was confined to splice questions (items 9 and 10 of the ledger), two conditions of the guidelines that the criteria had not applied.
The first review (v7 to v8) added 2.8 points.
The next run, v9, lost 2.7 points, and that loss was the correct result: PS1 and PM1 had been counting each variant's own ClinVar record, and 94% of PS1's firings were self-matches, so scoring them leave-one-out removed agreement that had never been evidence.
Of the 2.7 points, about 0.6 were the cost of the review's criteria fixes and about 2.0 were the leave-one-out scoring (RUN_VERSIONS.md).
Recall of ClinVar's pathogenic-side variants fell from 63.2% at v8 to 45.6% at v9 for the same reason, and the 62.5% of v37 is 0.7 points below the 63.2% that v8 reached by counting each variant's own record.
Agreement then rose by 3.2 points at v10, which bundled a rebuilt gnomAD database with its quality fields, the gene-level gates, and the first two measured thresholds (the PM2 bar and the BS2 prevalence bar), by 0.9 points each for the SVI point system (v11) and the BS1 bar (v12), and by 1.6 points at v13, when `--acmg` was made to enable HGVS and the criteria that read it came alive.

After v13 the headline stopped moving.
Across the next 23 runs (v14 to v37, without v17) same-direction agreement stayed between 77.50% and 77.62%, and no run moved it by more than 0.07 points.
The work did not stop, and its effect shows in the other panels.
Opposite-direction calls fell from 94 at v13 to 49 at v37, false-benign calls from 50 to 24, and false-pathogenic calls from 44 to 25 (Figure 1B).
The path was not monotonic.
Opposite-direction calls had reached a minimum of 25 at v10 and then rose to 56 at v11, when a lone PVS1 became Likely pathogenic under the point system, which gained 6,632 correct pathogenic-side calls for 31 additional false-pathogenic ones.
They rose again at v13, when BP7's intronic extension started to fire, which gained 11,339 benign-side calls for 36 additional missed diagnoses.
Both were trades made on purpose and recorded as trades, and Section 3.3 returns to how they were weighed.
The matrix stored for v1 holds 31 opposite-direction calls, fewer than v37, because v1 committed to only 74,355 variants where v37 commits to 244,093; per committed call, contradictions fell from 0.042% in v1 to 0.020% in v37, with a minimum of 0.011% at v10.

Recall of ClinVar's pathogenic-side variants in direction has stayed between 62.0 and 62.5% since v11, and recall of its benign-side variants between 62.3 and 62.7% since v13 (Figure 1C).
That flat line is the structural fact about the benchmark that Section 3.5 examines.
The agreement figures also have an uneven past: the raw output of v1 to v6 was not kept, and v2 and v4 to v6 are quoted from the figures recorded at the time with denominators that we could not check, so their values carry less weight than the stored runs from v7.

### 3.2 What review of the discordant calls found

**Review of the discordant calls first identified 30 of the 64 changes, by our coding, and the largest class was a guideline condition left out.**
The ledger of Section 2.3 holds 64 items, each coded by us for its class and for how it was first found (Figure 3, Supplementary Table S2).
Expert review first found 30 of the 64 (47%), code audit 8, our own reading of the discordant calls 7, the first real-data run 5, an aggregate figure of the benchmark 5, a threshold sweep 4, a comparison with Ensembl VEP 3, and a regression in a later run 2.
All of the review's input was the benchmark's own list of opposite-direction calls, so this compares the reviewer's reading of that list with other routes and not two independent ways of finding the same thing, and "first" is our reading of the written record.
Two of the seven that we found by reading the discordant calls ourselves, the self-matches of PS1 and PM1 and the two modes of reporting that followed, came from reproducing the reviewer's rows, which could not have shown them.
By class, 21 items were a condition of a guideline that the criterion had not applied: BS2's "full penetrance expected at an early age", PM1's "without benign variation", BP1's need for a truncating-only mechanism, PVS1's need for a gene-disease relationship and for loss of function as the mechanism, and BP4's confinement to variants for which a calibrated predictor exists.
Seventeen were wiring defects, in which a criterion or an input was inert, misread or fed the wrong field.
Eight concerned how far a population frequency can be trusted, five were annotation differences from Ensembl VEP, five were thresholds set by measurement, four concerned scoring and reporting, two were circularity with ClinVar, and two were combining rules.

**The reviewer supported ClinVar's call far more often than the software's.**
The medical geneticist returned annotated workbooks in three rounds (Figure 2A).
In round 1 the workbook held 101 opposite-direction calls from v7, and the reviewer judged ClinVar's call correct on 86, fastVEP's on 3, neither on 8, and left 4 with a colour that was not in the key.
Round 2 held 122 calls from v8, of which 93 were new and 29 carried over: ClinVar's call was judged correct on 65, fastVEP's on 2, neither on 8, more evidence was needed on 1, and 46 carried no recorded verdict.
Seven of the 65 (four in CYP21A2 and three in STRC) carry the marker for a gene with a homology problem, and on six of them it is on fastVEP's call as well as the green on ClinVar's; the project's note of that round counted those six as not marked, which gives 59 and 52 where we count 65 and 46.
Round 6 held 64 calls from v26: by our coding of the returned text, the reviewer supported ClinVar on 43, fastVEP on 4, neither on 9, and held 8 pending more evidence or by disorder.
So when the software and ClinVar disagreed in direction, the reviewer sided with ClinVar on 89% of the 97 rows of round 1 that carried a verdict in the key, 86% of the 76 rows of round 2 that carried one, and 67% of the 64 rows of round 6.
The three percentages come from different instruments, which are colour decoding against a two-colour key (round 1), colour decoding against a richer key (round 2) and our coding of free text (round 6, where 14 rows were coded from the evidence they cite and not from a stated class), so they are not one series and are by our coding.
The list of discordance categories that came with round 1 names six kinds of discordance: the use of PM1, PS1 or PM5 where it does not apply, a pathogenic call that rests on a founder or low-penetrance variant, and loss of function that is not the disease mechanism.
Three of the kinds match what the ledger later recorded: PM1 and PS1 applied where they did not apply (items 17 and 23), a founder or low-penetrance variant that is an exception to the benign criteria (items 26, 52 and 63), and loss of function that is not the mechanism (item 33).

**The first review saw about a third of the discordant calls.**
The workbook of round 1 listed 101 opposite-direction calls, but v7 had 313.
The scorer wrote its discrepancy log with a cap of 10,000 rows; the soft discrepancies (a ClinVar pathogenic call scored VUS) filled it before the end of the genome, and the export stopped at chromosome 5, silently dropping every opposite-direction call after it.
The cap was removed for opposite-direction calls in v8, whose export is complete, so 212 of v7's 313 had never been listed.
Nothing about the 101 rows was wrong, but the sample was the first five chromosomes and not a random one.
It was enough to expose three gaps in how criteria were executed (PVS1's grading, BP4's SpliceAI path and PM2's gnomAD lookup) and no more.

**Wiring defects left criteria inert for months, and for three of them the headline barely moved.**
Eight of the 17 wiring defects were found after v9.
Five of those eight were criteria that the documentation described and the code did not deliver.
`--acmg` did not enable HGVS, so PVS1's splice gate, BP7's intronic extension and the BA1 exception list all read an empty string, and no run from v9 to v12 passed the flag (found at v13).
BP3 had never fired, because its repeat database was never built (v15).
PS1's splice path had never fired, because the field it reads was never populated (v16).
The BA1 exception list had never matched a variant, because its matcher compared a full HGVS string with a bare `c.` token, and the unit tests passed because they fed it the bare token that the pipeline never produces (v18).
BP7's exon-edge exclusion was described and was read by the criterion, but the pipeline never set the field (v37).
The other three were of a different kind: two criteria read an HGVS display string and broke when that string changed (v28 to v35, found as regressions), and PVS1 read the end of a deletion's span and not the offset the span covers, which cost 14 ClinVar-pathogenic donor and acceptor deletions their PVS1 (v26, found by mining the discordant calls).

Whether the benchmark could see one of these when it was fixed depends on what the benchmark measures.
Building BP3's database moved same-direction agreement by 24 calls (0.004 points) and no opposite-direction call.
Switching on PS1's splice path left same-direction agreement unchanged and moved exact agreement by 2,952 calls (0.44 points), because it turned 3,793 Likely pathogenic calls into Pathogenic, a change the same-direction measure does not count.
Making the BA1 exception list match for the first time changed no cell of the concordance matrix at all.
The reason is structural: seven of the nine alleles of Ghosh et al. (2018), which are nine variants in seven genes, are not in the April reference set.
In the June and October releases six of the seven are "conflicting classifications", and the seventh (HFE p.Cys282Tyr) is a Pathogenic record with "low penetrance" and "risk factor" terms, which the October 2-star+ set does include.
The reference therefore discards most of the population that the exception list exists for, and no aggregate computed against it can validate that criterion, which is why an end-to-end test over all nine alleles covers it instead.

**Reading the sources again for this paper found departures that the project's documentation had described as following them.**
While preparing this paper we read each recommendation that we cite for a number or a grade against its source and compared the implementation with it, and checked the inputs against independent listings.
Eight findings came out of it, and the classifier was not changed for any of them.
First, PVS1 grades an NMD-escaping truncation outside a critical region as Moderate when at least 10% of the protein is removed and Supporting when less is, where the tree of Abou Tayoun et al. (2018) gives Strong and Moderate, and it grades a canonical splice variant in the last exon as Moderate whatever the fraction removed, where the tree has no last-exon branch and grades by critical region and fraction.
Regrading both by one step would change 7 of 673,660 calls (6 from VUS to Likely pathogenic and 1 from Benign to VUS), by a replay that is an upper bound, because the recorded code does not say which branch gave a PVS1_Moderate; so few change because a grade one step lower plus PM2 seldom reaches 6 points either way.
Second, PVS1 predicts NMD by whether the variant lies in the last exon, where the tree predicts it by the 50-nucleotide rule, which also lets a stop in the last 50 nucleotides of the penultimate exon escape; the exact rule is a setting, the code comment calls it the guideline-faithful one, and it is off by default because of what it costs against ClinVar: 603 calls in direction and 0.73 points of pathogenic-side recall, to avoid one false-pathogenic call (Section 3.4).
That is a default chosen on agreement with ClinVar, a tuning to the reference that is not the dependence on ClinVar-derived inputs that Section 3.4 measures.
Third, BP7 keeps the 2015 test of poor conservation (PhyloP below 2.0) and takes a SpliceAI score below 0.2 as no impact, where Walker et al. (2023) apply BP7 after BP4, which for SpliceAI is a score of 0.1 or below, and caution against a conservation test; Section 3.3 measures it.
Fourth, the allele-number floor has a narrower source than its use (Section 3.3).
Fifth, the recessive-gene PM2 bar of 7e-5 was attributed in the code to the SVI's PM2 recommendation, which gives no frequency at all; the number equals the recessive bar of the Hearing Loss expert panel, the project's record does not say where it came from, and Section 3.3 sweeps it.
Sixth, the veto of BP1 counts the variant's own ClinVar record, which PS1 and PM1 exclude (Section 3.4).
Seventh, the ClinGen file that the benchmark read lacked 104 curations that ClinGen dates before the day of the export, among them CFTR (Section 3.5).
Eighth, three citations were wrong: the PMID of Ghosh et al. (2018), which was wrong in three places, the journal of Tavtigian et al. (2020), which is Human Mutation and not Genetics in Medicine, and the date of the ClinGen guidance on genes with multiple disorders, which is January 2024 and was recorded as July 2025.
We corrected the citations and the comments in the repository, and Sections 3.3 to 3.5 measure what each departure would change, so that the decision to change the classifier can be made on a number (Section 4.5).

### 3.3 Where the thresholds sit

Each threshold that no guideline states was set once, in a particular run, against the combining rule and the other criteria as they then stood.
Re-measuring them on the final classifier changed what we can say about several (Figure 4, Table 3).

**Table 3. Thresholds that no guideline states: shipped values, and what the final classifier shows.**
Exchange rates are extra correct benign-side (or pathogenic-side) calls per extra missed diagnosis (or false-pathogenic call) for the step on either side of the shipped value; "free" means no added error.
Supplementary Table S3 gives every step of every sweep.

<!-- table:thresholds:begin -->
| Threshold: shipped value | How it was set | Last step before the shipped value | First step beyond it | On the final classifier |
|:---|:---|:---|:---|:---|
| PM2 bar, genes not read as recessive-only: 4e-5 | Swept in v10 under the 2015 table; the median of the published VCEP bars is 4e-5 (64 genes) | 2e-5 to 4e-5: +13 calls, +1 (13 per added error) | 4e-5 to 1e-4: +6 calls, +0 (no added error) | Recall 62.0% to 62.5% across strict absence to 5e-4 under the point system; 38.8% to 57.4% under the 2015 table |
| PM2 bar, genes read as recessive-only: 7e-5 | Not from the SVI PM2 recommendation, which the code cited and which gives no frequency; equal to the Hearing Loss expert panel's recessive bar. First swept here | 2e-5 to 7e-5: +0 calls, +0 (no added error) | 7e-5 to 2e-4: +6 calls, +0 (no added error) | Flat: recall 62.46% to 62.48% and agreement in direction within 10 calls from strict absence to 5e-4 |
| BS1 bar: 5e-3 | Swept on a 1-in-10 sample, where the move from 1e-2 looked free; on the full set it cost 4 more missed diagnoses for 2.3 points of benign-side recall (6,443 calls), and was kept (v12). The median of the published VCEP bars is 2.5e-4 | 1e-2 to 5e-3: +5,289 calls, +2 (2,645 per added error) | 5e-3 to 2e-3: +5,565 calls, +6 (928 per added error) | Knee at 1e-3 (112 per added error below it); 5e-3 is on the cautious side of the knee |
| BS2 homozygote bar, every gene not read as dominant-only: 1e-3 | Swept on the full benchmark (12 Aug): no knee, so set on the prevalence of the commonest recessive conditions | 1e-2 to 1e-3: +304 calls, +0 (no added error) | 1e-3 to 1e-4: +1,389 calls, +3 (463 per added error) | No knee: 463, 291 and 533 calls per added error below the shipped value |
| BS2 carrier bar, genes read as dominant-only: 3e-4, with an allele count of at least 5 | Added after review of one gene (SNRNP200) and swept on the full benchmark (v37) | - | - | Flat: 24 missed diagnoses and 25 false-pathogenic calls at every setting from off to 3e-3 |
| BP1 veto: 3 pathogenic missense variants in the gene | Added in v9 with the veto, after the round-2 review of BP1 firing in genes with a pathogenic missense spectrum; the count was never swept before this paper | 2 to 3: +17 calls, +1 (17 per added error) | 3 to 5: +79 calls, +3 (26 per added error) | At most 132 calls per added error at every step, so there is no knee and the shipped count is on the aggressive side of the stricter ones; a count of 1 gives +1,716 calls in direction and 149 fewer correct benign-side calls, and 2 fewer missed diagnoses (1 fewer when the veto leaves the variant's own record out) |
| BP7 far boundary: 300 nt | Set at the knee of a sweep (v14) | 200 to 300: +157 calls, +1 (157 per added error) | 300 to 500: +40 calls, +1 (40 per added error) | Knee confirmed: 157 to 849 calls per added error below 300, about 26 to 40 above |
| BP7 PhyloP bar: 2.0 | In the first commit, with no source recorded (the 2015 standards say "not highly conserved" and give no number). First swept here | 1.5 to 2: +4,571 calls, +3 (1,524 per added error) | 2 to 3: +6,539 calls, +0 (no added error) | No knee: at least 1,500 calls per added error at every step; switching the test off adds 16,482 benign-side calls and 2 missed diagnoses (Section 3.3) |
| Allele-number floor for BA1 and BS1: 2,000 | From the dataset-size condition of the BA1 recommendation, applied here to the site's allele number and to BS1 | - | - | Removing it adds 774 calls in direction and no opposite-direction call |
| Likely pathogenic from 6 points | Fixed by the point system (Tavtigian et al., 2020); varied here as a sensitivity analysis | - | - | 5 points: +10,523 pathogenic-side calls, +5,947 ClinVar-VUS and +19 benign-side called pathogenic-side (Table 4) |
| Pathogenic from 10 points | Fixed by the point system | - | - | 9 to 12 points move agreement in direction by at most 3 calls |
| BA1 at 5% | Fixed by the guideline; the median of the published VCEP bars is 1e-3 | - | - | Not varied |
<!-- table:thresholds:end -->

**The PM2 bar now matters little, because of the combining rule.**
The classifier uses one PM2 bar for genes that it does not read as recessive-only (dominant, both, or of unknown inheritance, where inheritance is read from the names of the diseases in the ClinGen file: of its 2,924 genes the classifier reads 67 as recessive-only, 50 as dominant-only, 12 as both and 2,795 as neither, Section 4.4) and another for genes that it reads as recessive-only.
When the first was moved from strict absence to 4e-5 in v10, the sweep of that time recorded a rise of 19 points in pathogenic recall, with a knee.
On the final classifier the same sweep is almost flat.
Across strict absence and every bar from 1e-6 to 5e-4, recall of ClinVar's pathogenic-side variants moves from 61.97% to 62.49%, which is 490 variants of 93,812.
It reaches 62.48% at the shipped 4e-5, and the two steps above that change it by 6 and 5 variants.
False-pathogenic calls rise from 22 at strict absence to 25 at 4e-5 and 26 at 2e-4, and agreement in direction is higher at strict absence (523,414 calls) than at the shipped bar (522,886), because ClinVar-VUS variants that carry PM2 plus other evidence stay VUS.
The bar for recessive-only genes, which the code attributed to the SVI's recommendation, which gives no frequency, and which equals the recessive bar of the Hearing Loss expert panel, is flat as well: from strict absence to 5e-4, recall of pathogenic-side variants is between 62.46% and 62.48%, agreement in direction varies by 10 calls, and the numbers of missed diagnoses and false-pathogenic calls do not change.
The two bars therefore act on nearly every gene and on 67 genes, and a flat sweep says what a bar does to the calls it reaches, not what a bar would do if inheritance were read from a correct field.
The reason the knee disappeared is the combining rule, not the bar.
Under the 2015 combining table, which is one switch away in the same classifier, the same sweep is steep (Figure 4B).
Recall of pathogenic-side variants is 38.8% at strict absence, 46.8% at 1e-6, 54.1% at 1e-5, 56.7% at 4e-5 and 57.4% at 1e-4.
Raising the bar from strict absence to 4e-5 gains 16,785 of 93,812 variants for 7 more false-pathogenic calls, which is the knee that v10 recorded: the sweep then gave 37.8% at strict absence and 56.8% at 4e-5, against 38.8% and 56.7% here.
Under that table a lone PVS1 matches no row, and a null variant whose only companion is PM2 reaches Likely pathogenic only through the rule that the SVI added with its PM2 recommendation, one Very Strong criterion plus one Supporting criterion (ClinGen SVI, 2020).
Of the 51,348 pathogenic-side calls that the table makes with PM2 and would not make without it, 45,855 (89%) rest on PVS1 and PM2 and nothing else.
The point system counts a lone PVS1 as 8 points, so under it only 5,493 calls depend on PM2, and none of them is a PVS1-and-PM2 pair.
A bar can move only the calls that need the criterion it gates, and the combining rule decides how many do.

**The BS1, BS2 and BP7 bars trade benign-side calls for missed diagnoses, and each has a different shape.**
All three pass through the shipped point (Figure 4A), and Table 3 gives what each step buys and costs.
BP7's far boundary reproduces the knee found in v14.
Moving it from 200 to 300 nucleotides buys 157 correct benign-side calls for one more missed diagnosis, and moving it from 300 to 500 buys 40 for one more, with 51 for two more beyond that and 457 for fourteen with no bound at all.
The shipped 300 is at the corner: below it a step is worth 157 to 849 calls per missed diagnosis, above it about 26 to 40.
Turning the extension off is not free either: it gives up 9,334 correct benign-side calls to avoid 11 missed diagnoses (the step from 0 to 50 nucleotides).
The BS2 homozygote bar, which applies to every gene that is not read as dominant-only, has no knee, as in the original sweep.
The steps from 1e-3 to 1e-6 are worth 463, 291 and 533 calls per missed diagnosis, so no setting is picked out by the data and the choice rests on how the two errors are weighed.
The BS1 bar is the one that moved.
Between 2e-2 and 5e-3 (shipped), lowering the bar buys 11,489 benign-side calls at the cost of 2 missed diagnoses.
Below the shipped value the exchange rate is 928 calls per missed diagnosis from 5e-3 to 2e-3, 625 from 2e-3 to 1e-3, and then 112 from 1e-3 to 5e-4 and 85 from 5e-4 to 2e-4.
By the knee rule that fixed BP7's bound, which looks for the point past which the exchange rate falls from hundreds of calls per missed diagnosis to about a hundred or fewer, the knee of the BS1 curve is at 1e-3, and the shipped 5e-3 sits on its cautious side: 3.1 points of benign-side recall (8,692 calls) and 11 missed diagnoses (24 against 35) separate them.
The published bars point the same way.
In the table of ClinGen VCEP specifications we extracted, the median BS1 bar is 2.5e-4 (107 genes with a usable bar, 102 of them stricter than our 5e-3), and the median BA1 bar is 1e-3 against 5e-2 (115 genes, 111 stricter).
The median PM2 bar is 4e-5 (64 genes, range 1.3e-6 to 5e-4), which is the value we ship for every gene that is not read as recessive-only.
That median leaves out the panels that give no single number: of 137 PM2 entries the extraction took 64 as bars, and set aside 37 that require absence, 14 with no threshold, 9 whose bar depends on inheritance (the Hearing Loss panel's recessive bar of 7e-5 among them), 6 in superscript notation, 4 where two panels disagree and 3 that are ambiguous or unspecified.
Across BA1, BS1 and PM2 together the extraction read 411 criterion entries and took 286 as bars (Table 1).
We have not changed the shipped value in this work.
A global frequency bar is a poor stand-in for a bar that depends on a disorder's prevalence and penetrance, which is what "greater than expected for the disorder" asks for, and that is a decision for a reviewer, not for a sweep.

**BP7's conservation test gives up more benign-side calls per missed diagnosis avoided than any other threshold in the ranges swept, and it had never been swept.**
BP7 asks that the nucleotide not be highly conserved, which the classifier reads as a PhyloP score below 2.0, a number from the first commit with no source recorded.
Raising the bar from 1 to 1.5 gains 5,520 correct benign-side calls with no added missed diagnosis, and from 1.5 to the shipped 2 gains 4,571 for 3 more (1,524 per added error).
Beyond the shipped value, raising it to 3 gains 6,539 for none, to 4 gains 3,520 for one, and switching the test off gains 6,423 for one more.
Every step buys at least 1,500 calls per added missed diagnosis, so unlike the BS1 bar the curve has no knee, and the rate for switching the test off, about 8,200 correct calls per added missed diagnosis, rests on two events.
Switching the test off, which follows the caution of Walker et al. (2023) against a conservation test without an empirical justification, adds 16,482 correct benign-side calls and 127 ClinVar-VUS variants called benign-side, raises agreement in direction by 16,355 calls (2.4 points) and raises missed diagnoses from 24 to 26.
These gains are measured against ClinVar's benign-side labels, which submitters may themselves have assigned by applying BP4 and BP7, so for this criterion the reference is not independent of the evidence.
Asking that SpliceAI be at most 0.1, as BP4 does, removes 313 correct benign-side calls from v37, and asking it on top of switching the test off removes 402, neither changing an error, so the two together, which is BP7 as the recommendation describes it, give 538,839 calls in direction (80.0%), 15,953 more than v37, for the same two additional missed diagnoses.
We have not changed the classifier: the two variants (NF1 and MAPT) are ClinVar-pathogenic and have PM2 and BP4 beside BP7, and it is for the maintainers to weigh two missed diagnoses against 16,000 correct calls (Section 4.5).
On the 107,528 variants that were not in the April set, BP7 as the recommendation describes it adds 4,933 calls in direction and no missed diagnosis (Section 3.7).

**The count of BP1's veto had never been swept, and a stricter count calls fewer ClinVar-VUS variants benign-side and avoids two missed diagnoses, one of which rests on the variant's own record.**
BP1 is withheld in a gene that has at least 3 pathogenic missense variants in the ClinVar protein index, a count that was set with the veto in v9 and never varied (the veto itself avoids 27 missed diagnoses relative to none).
At a count of 1, BP1 is withheld in every gene with one pathogenic missense record: agreement in direction rises by 1,716 calls, missed diagnoses fall from 24 to 22, false-pathogenic calls do not change, and 149 correct benign-side calls are given up, while the benign-side calls on ClinVar-VUS variants fall by 1,865.
The rise is those 1,865 variants staying VUS, which is where ClinVar leaves 44% of its records, less the 149 correct calls given up; it is not a finding that the benign-side calls on them were wrong (Section 3.5).
The two missed diagnoses that a count of 1 removes are the BP1 pair of Section 3.6 (IQSEC1 and PRR12): the protein index holds one pathogenic missense record for PRR12, the variant itself, so its removal rests on its own record, and two for IQSEC1, the variant and p.Thr343Met, so its removal rests on its own record at a count of 2 and not at a count of 1.
With the variant's own record left out of the count as well (Section 3.4), a count of 1 therefore removes IQSEC1 only (23 missed diagnoses, with the same 1,716 calls and 149 correct benign-side calls) and a count of 2 neither (24).
At a count of 2 without that change the figures are 849 more calls in direction, 23 missed diagnoses (IQSEC1 removed) and 17 correct benign-side calls given up.
Above the shipped value each step buys correct benign-side calls at the price of missed diagnoses: a count of 5 gives 3 more, 10 gives 5 more and 25 gives 12 more.
No step gives more than 132 calls per added missed diagnosis, so the exchange rate never falls from hundreds of calls to about a hundred, the curve has no knee in the sense of Section 2.5, and the shipped count of 3 lies on its aggressive side; a count of 1 or 2 is the cautious choice, though what it avoids is two missed diagnoses on two variants, one of which (PRR12) rests on its own record (Table 3).
Both are among the 49 and were found in the set on which the count was swept, and on the new variants a count of 1 avoids no missed diagnosis (Section 3.7).
We have not changed it, and the count bears on the finding of Section 3.4 that the veto counts the variant's own record.

**The BS2 carrier bar for genes read as dominant-only is flat.**
Switching the bar off, so that an allele count of at least five alone suffices, and setting it at 3e-3 change the false-benign and false-pathogenic counts by nothing (24 and 25 throughout), on the 50 genes that the classifier reads as dominant-only.
The shipped 3e-4 recovers 70 pathogenic-side calls relative to no bar (58,540 to 58,610) and gives up 474 benign-side calls (177,680 to 177,206), and nothing above it recovers another.

**On the pathogenic side the points that define the classes move the result more than any numeric bar that we swept on the shipped classifier, and the evidence sits in an awkward place on that scale.**
Figure 5A shows, for variants that met no benign evidence, where ClinVar's classes fall on the pathogenic points the classifier reached (Supplementary Table S11).
Variants with 0 or 1 point (27,027 and 105,989 of them) are 3.9% and 7.9% ClinVar pathogenic-side.
At 2 points the share is 10.7%, at 3 points 21.9%, and at 4 points 26.5%.
At 5 points it jumps to 64.0% (10,485 of 16,390), with 35.9% ClinVar VUS; at 6 points the share is 59.1%; from 7 points it is 94.8% or more, reaching 99.5% at 10.
The 5-point group is a cluster: PM1 with PM2 and PM5, PVS1 at Strong strength with PM2, and PP3 at Strong strength with PM2 each reach exactly 5, one point short of Likely pathogenic.
Moving the Likely pathogenic cutoff from 6 to 5 points, in a real run of the classifier with that one constant changed (Table 4), raises recall of ClinVar's pathogenic-side variants from 62.5% to 73.7% (10,523 more variants).
It does so at the price of 5,947 more ClinVar-VUS variants and 19 more ClinVar benign-side variants called pathogenic-side, so 1.8 pathogenic-side variants are gained for each variant added that ClinVar does not call pathogenic-side (whether the VUS among them are errors is the question of Section 3.5).
Agreement in direction rises by 4,576 calls, and exact agreement falls by 3,624.
A cutoff of 4 reaches 76.7% recall, but 17,452 ClinVar-VUS and 139 benign-side variants are then called pathogenic-side, and agreement in direction falls below the shipped value.
Raising the cutoff to 7 or 8 gives up 3,257 or 3,532 pathogenic-side calls and avoids 3 false-pathogenic calls, and at 9 a lone PVS1 (8 points) is no longer Likely pathogenic, recall falls to 52.8% and false-pathogenic calls to 11.
No cutoff changes the number of missed diagnoses, which is 24 in every run.
Tavtigian's framework places a posterior of 0.90 at 6 points and 0.81 at 5 under its prior of 0.10, and the benchmark's own prior among defined classes is 0.25, so the comparison with ClinVar cannot tell which cutoff is calibrated; it can say what each costs.
The Pathogenic cutoff of 10 points does not move agreement in direction: moving it between 9 and 12 changes that agreement by at most 3 calls, because the cutoff separates two tiers of the same direction.
It does move exact agreement, which is 462,778 at 9 points against 425,869 at 10.
Of the 58,230 Likely pathogenic calls, 46,718 (80%) are at exactly 9 points, 45,855 of them PVS1 with PM2 and nothing else, and 41,234 of those 46,718 are ClinVar Pathogenic.
The tier boundary therefore falls one point above where a lone null variant with PM2 lands, and we have kept it at the framework's 10, because what lifts such a variant to Pathogenic in laboratory practice is case-level evidence that a single-sample VCF does not carry.

**Table 4. The two point cutoffs of the pathogenic side, by real runs.**
Each row is a run of the final classifier with one of the two cutoffs changed, which the code holds as constants (Section 2.6); the offline replay of Section 2.7 reproduces every row cell for cell except the first two, where it differs by 3 and 1 variants.
Rates are over 673,660 records.

<!-- table:cutoffs:begin -->
| Setting | ClinVar P/LP called P/LP (recall) | ClinVar VUS called P/LP | ClinVar B/LB called P/LP | ClinVar P/LP called LB/B | Same-direction agreement | Exact agreement |
|:---|---:|---:|---:|---:|---:|---:|
| Likely pathogenic from 4 points | 71,956 (76.7%) | 17,452 | 139 | 24 | 522,381 (77.5%) | 414,959 (61.6%) |
| Likely pathogenic from 5 points | 69,133 (73.7%) | 9,559 | 44 | 24 | 527,462 (78.3%) | 422,245 (62.7%) |
| Likely pathogenic from 6 points (shipped) | 58,610 (62.5%) | 3,612 | 25 | 24 | 522,886 (77.6%) | 425,869 (63.2%) |
| Likely pathogenic from 7 points | 55,353 (59.0%) | 1,356 | 22 | 24 | 521,885 (77.5%) | 427,452 (63.5%) |
| Likely pathogenic from 8 points | 55,078 (58.7%) | 1,340 | 22 | 24 | 521,626 (77.4%) | 427,432 (63.4%) |
| Likely pathogenic from 9 points | 49,558 (52.8%) | 1,166 | 11 | 24 | 516,280 (76.6%) | 426,330 (63.3%) |
| Pathogenic from 9 points | 58,613 (62.5%) | 3,612 | 26 | 24 | 522,889 (77.6%) | 462,778 (68.7%) |
| Pathogenic from 10 points (shipped) | 58,610 (62.5%) | 3,612 | 25 | 24 | 522,886 (77.6%) | 425,869 (63.2%) |
| Pathogenic from 11 points | 58,610 (62.5%) | 3,612 | 25 | 24 | 522,886 (77.6%) | 422,795 (62.8%) |
| Pathogenic from 12 points | 58,610 (62.5%) | 3,612 | 25 | 24 | 522,886 (77.6%) | 422,757 (62.8%) |
<!-- table:cutoffs:end -->

The benign side shows that the choice of scheme matters more than the choice of floor.
The shipped scheme applies the 2015 table to benign evidence.
Scoring both sides on the point scale with Likely benign from -1 points, as in the published point system, raises recall of ClinVar's benign-side variants from 62.3% to 73.9%, but raises missed diagnoses from 24 to 380, of which 218 are a lone BP4, and lowers recall of ClinVar's VUS variants from 97.2% to 69.2%.
Tightening the floor to -2 gives 49.6% recall and 62 missed diagnoses, and to -3 gives 37.0% and 26.
The table therefore does better than the floors of -2 and -3 on recall and on missed diagnoses alike, and the floor of -1 is a trade of 356 more missed diagnoses for 11.6 points of recall (Figure 5C).
These three are replays of the recorded criteria and not real runs (Section 2.7).

**An allele-number floor whose source is narrower than its use bought nothing.**
BA1 and BS1 require an allele number of at least 2,000 at the variant's site in gnomAD.
The figure comes from the BA1 recommendation, which asks for an allele frequency above 5% in a general continental population dataset of at least 2,000 observed alleles (Ghosh et al., 2018).
The project's documentation attributed it to the March 2024 communication on gnomAD v4, which states no allele-number floor and recommends the grpmax filtering allele frequency (ClinGen VCEP Review Committee, 2024).
Here it is applied to the allele number at the variant's own site, and to BS1 as well, and neither use is stated by either source.
Removing it raises same-direction agreement by 774 calls and changes no opposite-direction call (Section 3.4).

### 3.4 What each decision and each source buys, and how much of the agreement is circular

Figure 8 shows, for every ablation of Section 2.6, the change from v37 in agreement in direction, in recall of ClinVar's pathogenic-side variants and in opposite-direction calls.
A setting that lowers agreement with ClinVar is not for that reason wrong, because ClinVar rewards calls and a safeguard withholds some.
We therefore read each safeguard by what it removes in errors and what it costs in calls, as in Section 3.3.

**Counting a variant's own ClinVar record would add 14.3 points of pathogenic-side recall, and removing every ClinVar-derived input costs 0.77 points of agreement and raises the opposite-direction calls from 49 to 74.**
The shipped classifier excludes a variant's own ClinVar record from the evidence that PS1 and PM1 read (PM5 cannot match it, Section 2.3).
Letting it count, as the runs before v9 did, raises agreement in direction by 13,440 calls (2.00 points) and recall of pathogenic-side variants from 62.5% to 76.8%, which is the figure a self-referential evaluation would have reported, and removes 7 of the 49 opposite-direction calls.
Removing the ClinVar protein and splice index, which holds the classifications of other variants that PS1, PM1, PM5 and BP1 read, lowers agreement by 5,159 calls (0.77 points) and pathogenic-side recall by 2.1 points.
Removing the allele-level ClinVar record as well changes one further call.
The dependence on ClinVar-derived inputs is therefore small in agreement and larger in errors: without the index, missed diagnoses rise from 24 to 51, and from 24 to 52 with the allele-level record as well (74 opposite-direction calls against 49), and the 27 are the same variants as appear when BP1's veto is switched off, because BP1 is then no longer withheld in genes where missense variants are known to be pathogenic.
The veto is protective, but it protects by reading ClinVar's other classifications in the gene, which in a benchmark against ClinVar is part of the circularity that this section sets out to measure.
The veto itself counts the variant's own record, which PS1 and PM1 do not: in a measurement build in which the veto leaves it out, missed diagnoses are 27 and not 24, so 3 of the 27 that the veto avoids rest on the variant's own record (FBXO28, PCGF2 and SBF1, each a ClinVar-pathogenic missense variant whose own entry is the third of the three that the veto needs).
A benchmark against ClinVar cannot separate that dependence from the classifier's own skill, because other variants' ClinVar classifications are evidence in the guideline's own terms (PS1 and PM5 are defined by them).
Two further adjustments also read the variant's own ClinVar record or followed the reviewer's rows: ClinVar's low-penetrance gate (1 missed diagnosis when removed) and the curated frequency-exception list (6).
Removed together with the veto's own-record count, the three raise missed diagnoses from 24 to 34, exactly the sum of the three alone (3, 1 and 6), and opposite-direction calls from 49 to 59, so 59 is the figure when the adjustments that read the own record or followed the reviewer's rows are removed together, and the veto itself and the other changes that followed the review remain.
Three of the ten added calls are Benign calls on ClinVar-pathogenic variants (one SERPINA1 variant once the low-penetrance gate is gone, and two G6PD variants once the exception list is), so the zero for Benign calls that Section 3.5 reports rests on two of the 15 curated alleles (both G6PD) and one ClinVar-derived gate.

**Three sources carry the benign side and one carries most of the pathogenic side.**
Removing SpliceAI lowers agreement in direction by 106,677 calls (15.8 points) and recall of benign-side variants from 62.3% to 24.6%.
Removing PhyloP lowers it by 79,399 calls (11.8 points; benign-side recall 34.2%), and removing gnomAD by 46,539 (6.9 points; 46.2%).
Without SpliceAI, 199,860 of the 248,731 BP4 firings are lost, and BP7, which needs SpliceAI and PhyloP together, is lost entirely (132,615 firings).
gnomAD feeds BA1, BS1, BS2 and PM2, and all four stop firing without it.
On the pathogenic side, removing the ClinGen gene-disease file lowers agreement by 43,193 calls (6.4 points) and recall of pathogenic-side variants from 62.5% to 15.7%, because for 48,545 of the 61,341 firings of PVS1 at any strength (79%) a ClinGen curation of Moderate strength or better is the only route by which the gene is taken to act by loss of function, the others being gnomAD constraint and a curated mechanism; the file also supplies the inheritance that PM2, BS2 and PM3 read, so the loss is not all PVS1's.
Removing gnomAD gene constraint takes out PP2 and BP1 entirely and 721 PVS1 firings and lowers pathogenic-side recall by 4.0 points, and it raises agreement in direction by 1,905 calls, because variants that lose their criteria become VUS and 44% of the reference is VUS.
Removing REVEL, the source specific to missense variants, lowers agreement by 0.78 points, and removing RepeatMasker (BP3) by 31 calls.

**The same two sources carry most of the missed diagnoses.**
Removing SpliceAI or PhyloP lowers missed diagnoses from 24 to 4, because BP7 is behind 20 of the 24.
BP7 is a criterion of high yield and high cost: 132,202 of its 132,615 firings (99.7%) are on ClinVar benign-side variants, and without it 79,399 calls in direction are lost to avoid those 20 missed diagnoses, about 4,000 calls for each.
Section 3.3 shows that much of that price is paid for the conservation test and not for BP7's other conditions, and Section 3.6 returns to the 20.

**Most design decisions pay for themselves, and three do not clearly do so.**
Four safeguards remove errors and cost no calls: the BP1 veto (keeping it gains 4,156 calls and avoids 27 missed diagnoses, and a stricter count would avoid 2 more, or 1 with the own record left out, Section 3.3), the gate on disease mechanism for PVS1 (147 calls and 13 false-pathogenic calls), the curated list of 15 frequency-exception alleles (6 missed diagnoses) and the ClinVar low-penetrance gate (1).
Reading absence from gnomAD as rarity for PM2 gains 1,499 calls and no error.
Three more remove errors at a price in calls: the gene-disease validity gate costs 193 calls to avoid 10 false-pathogenic calls (about 19 calls per error), the homologous-gene list 391 calls to avoid 2 missed diagnoses (about 196 per error), and the filtering allele frequency 2,263 calls to avoid 8 missed diagnoses (about 283 per missed diagnosis), while adding 4 false-pathogenic calls.
Their prices are all below the 928 calls per added error at which we kept the BS1 bar where it was, so by that standard each is worth paying.
The standard was not applied uniformly, which Section 4.4 notes.
Three cost more than they remove on this benchmark.
The gnomAD region flags (segmental duplications and low-complexity tracts) withhold the frequency criteria at a price of 3,305 calls to avoid one false-pathogenic call, the allele-number floor withholds 774 calls to avoid none, and PM1's requirement that the window hold no benign variation withholds 272 to avoid none.
We have not changed any of them, because safeguards whose targets are rare are poorly tested by a benchmark of consensus classifications.
The allele-number floor is the first candidate for removal: neither source it can be traced to states it in the form used here (Section 3.3), and the data do not support it.
Two opt-in settings change little: capping PP3 at Moderate raises agreement by 712 calls and lowers pathogenic-side recall by 0.34 points with no change in errors, and the exact 50-nucleotide NMD rule for PVS1 loses 603 calls to avoid one false-pathogenic call.

**Two settings would trade errors for recall in the other direction.**
PM2 at Moderate, as in the 2015 standards, raises recall of pathogenic-side variants by 11.0 points (10,280 more variants) and agreement in direction by 4,475 calls, for 14 more false-pathogenic calls (about 320 calls per error).
The 2015 table on both sides loses 5,255 calls and 5.8 points of recall to avoid 11 false-pathogenic calls, which is the combining-rule effect of Section 3.3.
Neither is a clear improvement by the standard of this paper, and the first contradicts the ClinGen recommendation that PM2 be used at Supporting strength.

### 3.5 Discordance of the final classifier with ClinVar 2-star+ classifications

On the benchmark snapshot, v37 reproduces ClinVar's class for 425,869 of 673,660 variants (63.2%, 95% interval 61.7 to 64.6) and its direction for 522,886 (77.6%, 76.8 to 78.4; gene-cluster bootstrap over 10,232 groups of variants, which are 10,231 genes and one group of the 4 records that name none).
Table 5 divides every record by how the call relates to ClinVar's class (Figure 6A).

**Table 5. Every ClinVar 2-star+ record, by how the v37 call relates to ClinVar's class.**

<!-- table:decomposition:begin -->
| Kind | Cells | Variants | Share of all records |
|:---|:---|---:|---:|
| Same class as ClinVar | the diagonal | 425,869 | 63.2% |
| Same direction, different tier | P and LP, LB and B | 97,017 | 14.4% |
| &nbsp;&nbsp;ClinVar Pathogenic called Likely pathogenic |  | 48,301 | 7.2% |
| &nbsp;&nbsp;ClinVar Benign called Likely benign |  | 46,282 | 6.9% |
| ClinVar P/LP called VUS | a pathogenic-side class, no call made | 35,178 | 5.2% |
| ClinVar LB/B called VUS | a benign-side class, no call made | 106,937 | 15.9% |
| ClinVar VUS called P/LP |  | 3,612 | 0.5% |
| ClinVar VUS called LB/B |  | 4,616 | 0.7% |
| Opposite direction |  | 49 | 0.007% |
| &nbsp;&nbsp;ClinVar P/LP called LB/B (missed diagnoses) |  | 24 | 0.004% |
| &nbsp;&nbsp;ClinVar LB/B called P/LP |  | 25 | 0.004% |
| No call | reference-agreement records (ALT is a dot) | 382 | 0.06% |
| All records |  | 673,660 | 100% |
<!-- table:decomposition:end -->

**The discordance is almost entirely one side being undecided, and contradictions are rare.**
Of the 150,774 records that are not in the same direction as ClinVar, 142,115 (94.3%) are a pathogenic-side or benign-side ClinVar class that v37 leaves as VUS, 8,228 (5.5%) are the reverse, 382 have no call, and 49 (0.03%) are opposite.
The opposite-direction rate is 0.0073% of records (95% interval 0.0051 to 0.0098% by gene-cluster bootstrap; 0.013% of the records whose ClinVar class is pathogenic-side or benign-side), made up of 24 missed diagnoses (0.026% of ClinVar's 93,812 pathogenic-side variants) and 25 false-pathogenic calls (0.009% of its 284,550 benign-side variants).
A further 97,017 records (14.4%) are in the right direction at the wrong tier, and 94,583 of them are ClinVar Pathogenic called Likely pathogenic (48,301) or ClinVar Benign called Likely benign (46,282): exact agreement measures the tier boundary more than the direction.
The agreement is also asymmetric in what it says about each call.
Of the 62,247 variants called pathogenic-side, 94.2% (58,610; 92.6 to 95.5) are ClinVar pathogenic-side and 5.8% ClinVar VUS, and of the 181,846 called benign-side, 97.4% (177,206; 96.5 to 98.2) are ClinVar benign-side (Figure 6C).
Recall is much lower than precision on both sides: 62.5% of ClinVar's pathogenic-side variants (58.2 to 66.6) and 62.3% of its benign-side variants (60.6 to 63.9) receive a call in the same direction, against 97.2% of its VUS variants (96.6 to 97.7).
The wide interval of the first is real: recall is a property of the gene.
Among the 1,273 genes with at least ten ClinVar pathogenic-side variants, 232 (18%) have none called pathogenic-side and 231 (18%) have 90% or more (Supplementary Table S9).

**Agreement in direction does not separate cautious from aggressive calling.**
Calling every variant VUS agrees with ClinVar in direction for 43.8%, and rules that read the inputs of some ACMG criteria but apply no combining rule, gene-level gate or safeguard (Supplementary Table S13) reach v37's agreement and more only by making many more contradictions (Figure 9).
Calling a variant benign-side when its population-maximum allele frequency is 5e-3 or more gives 61.6% with 159 contradictions, and at 1e-3 (the median published BA1 bar) 64.0% with 813; a higher bar has fewer contradictions and agrees less, 52 contradictions at 59.7% for 1e-2 and 22 at 57.1% for 2e-2, so with the number of contradictions that v37 has the rule agrees for about 18 points less.
Adding a pathogenic-side call for a stop-gained, frameshift, canonical splice or start-lost variant that is absent from gnomAD or below 1e-4 gives 72.8% with 867, and adding a benign-side call for a quiet variant with SpliceAI below 0.2 gives 84.6%, 6.9 points above v37, with 1,730 contradictions (1,676 missed diagnoses and 54 false-pathogenic calls), 35 times as many as v37.
Without the frequency rule, the last two rules give 77.9%, 0.3 points above v37, with 950 contradictions (896 missed diagnoses and 54 false-pathogenic calls), 19 times as many, and no bar of the first rule brings the three rules below 950 contradictions.
The last rule is not a bargain by the standard of Section 2.5 either: over the two-rule version at 1e-3 it adds 86,510 correct benign-side calls for 863 more missed diagnoses, 100 per error.
The population maximum used by these rules is a point estimate over ten groups, so it is more aggressive than the filtering allele frequency that BA1 tests (in v37 a point estimate in place of that frequency costs 8 more missed diagnoses, Section 3.4), but it is the bar that sets the contradictions of the first rule (22, 52, 159 and 813 at 2e-2, 1e-2, 5e-3 and 1e-3).
The rules are not worse than v37 at everything: the two-rule version has a higher recall of ClinVar's pathogenic-side variants (66.0% against 62.5%) and a slightly higher share of its pathogenic-side calls that ClinVar calls pathogenic-side (95.4% against 94.2%), at 54 false-pathogenic calls against 25, which is in part the PVS1 gap below seen from the other side: of the 9,967 pathogenic-side calls that the second rule makes correctly and v37 leaves as VUS, 5,749 carry no PVS1 in v37 and 4,218 carry it at a lower strength or held back, and 43 of its 54 false-pathogenic calls are VUS in v37.
So agreement alone cannot say which of two classifiers is the more cautious, and v37 is better read with its contradictions beside its agreement.
The tiers of v37 are informative in the same way, with their denominators: none of its 4,017 Pathogenic calls (upper end of the exact 95% interval 0.092%) and none of its 67,700 Benign calls (0.0054%) contradicts ClinVar, against 25 of 58,230 Likely pathogenic calls (0.043%) and 24 of 114,146 Likely benign calls (0.021%).
The zero for Benign is informative, because 14 would be expected at the Likely benign rate, and the zero for Pathogenic is not, because 1.7 would; both hold for this in-sample set only, and the zero for Benign rests on two of the 15 curated alleles (both G6PD) and one ClinVar-derived gate, since three Benign calls contradict ClinVar when they are removed (Section 3.4).
Its 4,017 Pathogenic calls are ClinVar pathogenic-side for 3,996 (99.5%) and VUS for 21, and its 67,700 Benign calls are ClinVar benign-side for 67,455 (99.6%) and VUS for 245; the 25 false-pathogenic calls are all Likely pathogenic and the 24 missed diagnoses all Likely benign.

**Review status separates the classes, and the highest tier is one gene.**
The 3-star records, which are expert-panel classifications, are called pathogenic-side more often than the 2-star records (67.8% against 61.7% of their pathogenic-side variants), and called benign-side much less often (36.6% against 62.8% of their benign-side variants, with 63.3% left as VUS).
Expert-panel benign-side classifications can rest on disease-specific frequency bars, and the published bars are stricter than ours (Section 3.3), which is a possible reason for the difference and one that this benchmark cannot establish.
All 23 four-star records (practice-guideline classifications) are variants of CFTR, and v37 calls 22 of them VUS (Figure 7A).

**Pathogenic-side recall is set by the consequence class.**
Variants called pathogenic-side are almost all null variants: 54,498 of the 58,610 correct pathogenic-side calls (93.0%) are frameshift, nonsense or canonical splice variants, and 4,109 are missense (Figure 7A, Supplementary Table S5).
Among ClinVar's pathogenic-side variants, 88.4% of canonical splice variants and 82.6% of frameshift and nonsense variants are called pathogenic-side, against 17.8% of missense variants, 0.3% of in-frame variants and none of the start-loss, stop-loss, splice-region, intronic, synonymous, UTR and other variants.
Among ClinVar's benign-side variants the picture is the reverse: 79.2% of intronic, 65.7% of synonymous, 62.5% of UTR and 49.8% of missense variants are called benign-side (Figure 7B).

**Where the 35,178 pathogenic-side variants called VUS come from.**
They are not an even sample.
By the evidence that was met, 10,541 (30.0%) reached 5 points, one short of Likely pathogenic, 10,597 (30.1%) reached 0 or 1 point, and 14,027 reached 2 to 4 points; 13 reached 6 points or more and were held at VUS by a veto or by benign evidence.
The 5-point group is the cluster of Section 3.3 (PM1 with PM2 and PM5, PVS1 at Strong strength with PM2, and PP3 at Strong strength with PM2).
By consequence, 18,992 (54.0%) are missense variants, 9,129 (26.0%) frameshift or nonsense, 2,041 (5.8%) in flanking or non-coding sequence, and the rest splice-region, canonical splice, in-frame, start-loss and stop-loss, intronic, UTR and synonymous variants.
The concentration in genes is high: the 10 genes with most of them hold 11.7%, the top 100 hold 40.5% and the top 500 hold 74.7%, of 2,716 genes with any (CFTR 778, APC 712, LDLR 517, PAH 453, GCK 344, ABCA4 298, GLA 264, BRCA1 255, NF1 253 and SACS 246).

**A gap in PVS1 accounts for more than half of the truncating variants that are missed, and an incomplete input explains a small part of it.**
Of the 52,581 ClinVar pathogenic-side stop-gained and frameshift variants, 47,397 carry a PVS1 criterion at some strength and 5,184, in 625 genes (620 by the symbol of the annotation, which differs from ClinVar's first gene symbol for 55 of the variants), carry none (Supplementary Table S6).
Those 5,184 are 57% of the 9,129 that are called VUS.
PVS1 asks three questions in order.
The validity gate asks whether every ClinGen curation of the gene is weak, the mechanism gate whether the curated mechanism excludes loss of function, and the gene must then be taken to act by loss of function, which the classifier reads from one of three routes: gnomAD constraint (pLI at least 0.9 or LOEUF at most 0.35), a curated loss-of-function mechanism, or a ClinGen curation of Moderate strength or better.
The second route is nearly empty: of the 26 genes in the shipped table of mechanisms, 20 are gain of function, one is dominant negative and five name both loss and gain of function, and those five are each licensed by constraint or by a curation already.
The validity gate has one leak that changes no call: of the 505 genes whose ClinGen curations are all weak, one (NSMCE3) escapes it because its label ends in a semicolon, which the loader reads as a separator, and it has no variant in the benchmark.
For 4,696 of the 5,184 variants (91%, in 575 genes) none of the three routes is met: the gene tolerates heterozygous loss in gnomAD, has no curated mechanism and has no curation of Moderate strength or better in the file.
For 463 (in 38 genes) the validity gate withholds PVS1, and for 5 (in 4 genes) the mechanism gate does, and 20 meet a route and are not graded as PVS1 for another reason.
Absence has the effect here that it had in the first version of the validity gate (v10, Table 1), which blocked PVS1 where a gene had no curation and which v11 replaced with a gate that blocks only on curated-negative evidence: a gene with no curation is not blocked now, but with no other route to loss of function it is not licensed either.
All 76 variants in genes that gnomAD constrains (RASA1 40, TGFB3 13, LRP2 11, VPS13D 5, SLC12A5 3, DMXL2 2 and PTPN11 2) are held by a gate, 74 by the validity gate and 2 by the mechanism gate.
A gate that works on the gene cannot see which disease a variant is pathogenic for.
The only ClinGen curation of RASA1 in the file is a Disputed one for Noonan syndrome, where ClinVar's records name capillary malformation-arteriovenous malformation for 39 of its 40 variants, and the only curation of LRP2 is "no known disease relationship" with congenital heart disease, where ClinVar's records name Donnai-Barrow syndrome for 10 of its 11.
The ClinGen guidance on genes with several disorders (ClinGen, 2024), whose per-disorder scaffold is in the classifier and is not active, is the answer to that.

The largest single gene in the gap is an input that was incomplete.
CFTR holds 439 of the 5,184 variants, and its 558 null variants in the benchmark (stop-gained, frameshift and canonical splice) are all called VUS and none carries PVS1; all 23 four-star records are variants of CFTR, and v37 calls 22 of them VUS.
ClinGen's export of 11 October 2026 lists CFTR as Definitive for cystic fibrosis (autosomal recessive, classified on 1 June 2022), but the export of 30 April 2026 that the benchmark read does not contain it, and nor do 103 other curations that ClinGen dates before that day: 104 curations in 100 genes, 76 of them Definitive, Strong or Moderate, in 73 genes.
The completed file adds the 112 rows of the October export that are dated on or before 30 April and whose gene and disease pair the April file lacks: the 104 curations, seven rows that re-label, under another disease identifier, an assertion that the April file lists, and one row (GATM, Moderate) dated 04:00 UTC on 30 April, which is after the April file was saved (02:19 UTC, by its modification time) and changes none of the 163 calls on GATM's variants.
Their dates are recent: 83 of the 104 (80%) were classified in 2025 or 2026 and 45 (43%) from January to April 2026, against 664 of the 3,546 curations that the file lists (19%) and 159 (4.5%), where an omission independent of the date would leave about 19 and 5, which points to an export that lagged behind ClinGen's curations; five are older than 2024, CFTR (2022) among them, which a lag does not explain, and we have not established the cause.
We found it only by asking why CFTR had no PVS1.
Putting all 112 back, in a run that differs from v37 in that file alone, gives PVS1 to 674 of the 5,184 variants (13%, in 23 genes: CFTR 439, ELP1 85, EVC 46 and SPG7 41 among them, which are 611 of the 674) and a pathogenic-side call to 653.
Of CFTR's 558 null variants, 556 then carry PVS1 and 542 are called pathogenic-side, and 12 of the 23 four-star records are called pathogenic-side where 1 was.
The run adds 852 correct pathogenic-side calls and raises recall of pathogenic-side variants from 62.5% to 63.4%, and it adds 19 ClinVar-VUS variants called pathogenic-side and 2 false-pathogenic calls, both frameshift variants of ALPK1 that ClinVar calls Benign or Likely benign at allele frequencies of 2.3e-4 and 6.4e-4; missed diagnoses do not change.

Orphanet's loss-of-function label reaches more of the gap than the completed file does, and neither closes it.
Licensing PVS1 for the 1,042 genes that Orphanet labels as a loss-of-function disease and that the shipped table does not name (the association type "disease-causing germline mutation(s) (loss of function)" of Orphadata product 6, version of 23 June 2026, where a gene that also has a gain-of-function association is taken as both; Orphanet, 2026) gives PVS1 to 1,506 of the variants and a pathogenic-side call to 1,421, adds 1,829 correct pathogenic-side calls and 31 ClinVar-VUS variants called pathogenic-side, and adds 5 false-pathogenic calls (MAP3K6, ALPK1 twice, VPS13C and LTBP4, two of them splice-acceptor variants with population-maximum frequencies of 0.8%).
With the complete ClinGen file as well, 1,714 variants gain PVS1 and 1,624 are called pathogenic-side, agreement in direction rises by 2,053 calls (0.3 points), recall of pathogenic-side variants to 64.7%, and false-pathogenic calls by 5, about 411 calls per added error.
Even so 3,470 of the 5,184 variants (67%) carry no PVS1, and closing the rest needs a source of loss-of-function mechanism that does not depend on constraint, such as the PVS1 applicability stated in VCEP specifications, which our extraction of the CSpec registry does not read; we have not changed the classifier here.
In effect the classifier implements the standards' condition that loss of function be a known mechanism as gene-disease validity or constraint, which are proxies: a Definitive curation says that a gene causes disease and not that it does so by loss of function, and the curated-mechanism route licenses nothing in this benchmark.
The other 3,945 truncating variants that are called VUS carry PVS1 at a lower strength (3,062 at Strong, 463 at Moderate, 416 at Supporting and 4 at full strength but held back by other evidence), and 1,434 of the 1,446 canonical splice variants that are called VUS carry no PVS1 either.

**Where the 106,937 benign-side variants called VUS come from.**
Of the ClinVar benign-side variants that are called VUS, 31,912 (29.8%) carry PM2 and nothing else (they are rare in gnomAD and no benign evidence fires), 20,604 (19.3%) carry BP4 with PM2, 17,467 (16.3%) a lone BP4, 14,431 (13.5%) no criterion at all, 5,586 (5.2%) a lone BS1, and the rest smaller combinations.
Under the 2015 table that the benign side keeps, one supporting benign criterion does not reach Likely benign, and Section 3.3 shows what the alternative would cost.
At least 52,516 of these variants (49%) carry PM2, which is rarity and so does not separate a benign variant from a pathogenic one; 17,467 carry a lone BP4 and 5,586 a lone BS1, and neither reaches Likely benign under the 2015 table; and 14,431 carry no criterion at all.

**The 3,612 and 4,616 variants where v37 commits and ClinVar's submitters see insufficient evidence are not obviously errors.**
The pathogenic-side calls are PP2 with PM2 and PP3 at Strong strength (1,195), PVS1 with PM2 (1,135), PM1 with PM2, PP2 and PP3 at Moderate strength (587), PM1 with PM2, PM5 and PP2 (448), and a lone PVS1 (163).
The benign-side calls are chiefly BP1 with BP4 at Supporting or Moderate strength, with or without PM2 (3,677 of the 4,616).
Each pathogenic-side call stands on 6 or more points of computable evidence, which counts correlated predictors (PP2, PP3 and PM1 for a missense variant) as if they were independent, as the point system does, and each benign-side call on two or more benign criteria or on one at Very Strong or Stand-alone strength (122 of the 4,616).
A ClinVar VUS can mean that submitters weighed the evidence and found it insufficient, or that they did not weigh it to a class.
The benchmark cannot say which is right, and Section 3.7 asks whether ClinVar later moved toward them.

### 3.6 The 49 opposite-direction calls that remain

The last review round (v26) held 64 opposite-direction calls, and every one of the 49 in v37 is among them, so none of the 49 is new since that table was sent (one run, v31, briefly had 65 opposite-direction calls; Table 6, Figure 2B).
The 15 that were repaired all became VUS and none became the other class, and the reviewer had supported ClinVar's call on 13 of them (of the other two, one was coded neither and one depends on the disorder).
Each of the 49 had already been in the table of an earlier round (17 since round 3, 2 since round 4 and 30 since round 5), so they are the calls that several rounds of changes left standing.

**Table 6. The 49 opposite-direction calls of v37, by the evidence behind the call and the ruling of the last review round.**
Rulings are our coding of the reviewer's text (Section 2.4), not categories the reviewer assigned.

<!-- table:residual:begin -->
| Evidence behind the call | Calls | ClinVar supported | fastVEP supported | Neither | Needs more evidence | Depends on the disorder |
|:---|---:|---:|---:|---:|---:|---:|
| PVS1 on a ClinVar benign or likely benign variant | 22 | 14 | 2 | 3 | 2 | 1 |
| BP4 and BP7 on a ClinVar pathogenic or likely pathogenic variant | 20 | 10 | 2 | 5 | 3 | 0 |
| PM1, PP2 and PP3 on a ClinVar benign or likely benign variant | 3 | 3 | 0 | 0 | 0 | 0 |
| BS1 on a ClinVar pathogenic or likely pathogenic variant | 2 | 1 | 0 | 0 | 1 | 0 |
| BP1 on a ClinVar pathogenic or likely pathogenic variant | 2 | 2 | 0 | 0 | 0 | 0 |
| All | 49 | 30 | 4 | 8 | 6 | 1 |
<!-- table:residual:end -->

For 30 of the 49 (61%) the reviewer supported ClinVar's call (13 of the 24 missed diagnoses and 17 of the 25 false-pathogenic calls), for 4 (8%) fastVEP's, for 8 the reviewer's own call was VUS, for 6 more evidence was needed and for 1 the answer depends on the disorder.
So for 30 of the 49 the reviewer judged ClinVar's call the better, which is not the same as the classifier being wrong in each of them, and for the other 19 the reviewer did not side with ClinVar.
Of the 30 rulings for ClinVar, 11 follow from a class that the reviewer stated, 9 from the reviewer's statement that ClinVar is correct, and 10 from evidence that the reviewer gave for ClinVar's class without naming one (Supplementary Table S7).
Five of the 49 (10%) are no longer in the 2-star+ set in the October release, where ClinVar records them as having conflicting classifications, against 1.3% of all April variants that left it.
The 49 fall into five groups, and for most of them the reviewer's remarks name what decides the call, which is mostly evidence that a single-sample VCF does not carry.

**PVS1 on a variant ClinVar calls benign or likely benign (22).**
Eleven of these are a lone PVS1 and 11 are PVS1 with PM2, and all are loss-of-function variants in a gene that the classifier takes to be sensitive to loss of function.
Six are seen in 611 to 4,892 gnomAD alleles (EDN3, FMR1, KDM2B, PHKB, PIEZO1 and RAD21), at allele frequencies of 5.6e-4 to 3.3e-3 across gnomAD and population-maximum frequencies of 0.23% to 1.5%.
BS1, which tests the filtering allele frequency against 5e-3, did not fire for any of the six, and their population-maximum frequencies are 9 to 59 times the median published VCEP bar of 2.5e-4; the reviewer supported ClinVar on all six, remarking for five of them that the variant is too common in gnomAD.
Two of the six (KDM2B and PHKB) are also truncations in the first coding exon, where the reviewer noted that NMD can be escaped; neither the published PVS1 tree nor the classifier has a rule for it, and the classifier predicts NMD for every truncation outside the last exon.
The reviewer's other reasons were a haplotype in which the variant is always in cis with a second variant and the pair is not pathogenic (BRCA1); a splice-site variant in the 5' UTR of one RUNX1 isoform; a gene without an OMIM disease (RYK); a site that failed quality control in gnomAD (NOTCH2); a site covered in under half of gnomAD exomes (PRODH); and, for two variants (RYR1 and ARMC9), a call that is correct for the recessive disorder.

**BP4 and BP7 on a variant ClinVar calls pathogenic or likely pathogenic (20).**
These are the 20 missed diagnoses of Section 3.4.
Seventeen are intronic, two synonymous and one in the splice region, with SpliceAI below 0.2 and PhyloP below 2.0, which is what BP7 asks.
The reviewer supported ClinVar for 10, and for several noted genetic or functional evidence of an effect on splicing (POLR3A c.1909+22G>A, at an allele frequency of 2.5e-3 in gnomAD, is one: the reviewer cites affected individuals and some functional data).
For 9 the reviewer's remark is that ClinVar's classification appears to predate the current MANE transcripts, including the two TH variants on which the reviewer supported fastVEP.
One of the 20 (DST) is a stop-gained variant on the MANE Plus Clinical transcript and a deep-intronic variant on the MANE Select transcript that the pick chose, and the reviewer asked that Plus Clinical take precedence where both exist, which is a question of transcript policy that we have not changed.

**Missense variants: PM1, PP2 and PP3 on a variant ClinVar calls benign or likely benign (3), BS1 on pathogenic missense variants (2) and BP1 on pathogenic missense variants (2).**
The first three (CHD7, FBN1 and PTCH1) reach 6 points through PM1, PM2, PP2 and PP3 at Moderate strength, and the reviewer supported ClinVar on all three.
The BS1 pair (PADI3 and SLC12A3) are seen in 328 and 338 gnomAD alleles, with population-maximum frequencies of 0.6% and 0.8%; the reviewer noted a homozygote in gnomAD and a disorder (Gitelman syndrome) for SLC12A3, supported ClinVar there, and asked for a deeper look at PADI3.
The BP1 pair (IQSEC1 and PRR12) are in genes with two and one pathogenic or likely pathogenic missense records in the ClinVar protein index, the variant's own among them, so BP1's veto does not apply (a count of 1 would withhold BP1 for both, for PRR12 only through its own record and for IQSEC1 also through p.Thr343Met, Section 3.3), and the reviewer supported ClinVar on both.

No change of one threshold would remove these calls as a class, except that a stricter count of BP1's veto acts on the BP1 pair (Section 3.3).
A disease-specific frequency bar could act on about eight of them (the six common variants of PVS1 and the two of BS1), a transcript policy on a handful, and for the rest of the others the reviewer's remarks point to phase, case-level or functional evidence.
We have not put them on the curated exception list, which holds 15 alleles, because a list long enough to hide the residual would no longer be a statement about the classifier.

### 3.7 A later ClinVar release, and the variants that were not in the April set

ClinVar's release of 4 October 2026 holds 781,135 records at two stars or better: 665,057 of the April set, which are still at two stars or better, and 116,078 that are new.
A further 8,603 of the April variants have left the set, among them five of the 49 residual calls (Section 3.6).
The gnomAD extract the classifier reads covers every April variant, and 107,528 of the new ones; the rule of Section 2.8 drops the other 8,550, which is 7.4% of the new records and 1.1% of the release.
Most of the new variants were already in ClinVar in April, at a lower review status: 105,920 of the 107,528 (98.5%) have a record in the `variant_summary` of 29 April from which the protein index was built, and for 98,623 (92% of all the new variants) it is a single-submitter record.
They are variants that reached two stars since, and 9,234 of them already had a pathogenic or likely pathogenic record in April.
Table 7 gives the final classifier, unchanged, on the April set and on each stratum of the covered release (Supplementary Table S8).

**Table 7. The final classifier on the April benchmark and on the later ClinVar release, by stratum.**
A variant is new if it is not in the April set, unchanged if it has the same class in both releases, and has moved within or across a direction if its class changed.
Opposite-direction calls are given as the total, with missed diagnoses and false-pathogenic calls in brackets.

<!-- table:temporal:begin -->
| Set | Variants | Same-direction agreement | ClinVar P/LP called P/LP | ClinVar LB/B called LB/B | ClinVar VUS called VUS | Opposite-direction calls (missed diagnoses and false-pathogenic) |
|:---|---:|---:|---:|---:|---:|---:|
| April release, the benchmark set | 673,660 | 77.6% | 62.5% | 62.3% | 97.2% | 49 (24 and 25) |
| October release, covered by the gnomAD extract | 772,585 | 76.6% | 62.7% | 61.2% | 97.3% | 64 (26 and 38) |
| &nbsp;&nbsp;unchanged since April | 662,668 | 77.7% | 63.1% | 62.5% | 97.2% | 44 (20 and 24) |
| &nbsp;&nbsp;class moved within its direction | 2,232 | 54.1% | 37.4% | 59.0% | - | 0 (0 and 0) |
| &nbsp;&nbsp;class moved across a direction | 157 | 19.1% | 16.2% | 7.1% | 77.3% | 2 (0 and 2) |
| &nbsp;&nbsp;new: not in the April set | 107,528 | 70.1% | 61.0% | 55.2% | 97.3% | 18 (6 and 12) |
<!-- table:temporal:end -->

**On variants that were never in the benchmark, agreement in direction is 70.1% against 77.6%, mostly on the benign side, and contradictions are 2.3 times as frequent.**
Across the 772,585 covered records, agreement in direction is 76.6% against 77.6% in April, with 64 opposite-direction calls (26 missed diagnoses and 38 false-pathogenic).
On the 662,668 variants whose class did not change the figures are those of April (77.7%, 44 opposite), as they should be.
On the 107,528 new variants agreement in direction is 70.1%, and recall of pathogenic-side variants is 61.0% (April 62.5%), of VUS 97.3% (97.2%) and of benign-side variants 55.2% (62.3%).
The new variants are a different mix of ClinVar classes (9% pathogenic-side, 34% VUS and 57% benign-side, against 14%, 44% and 42% in April), and at the April recall of each direction their agreement would be 74.2%, so taking the mix first 3.4 of the 7.5 points are the mix and 4.1 are lower recall, nearly all of it benign-side (the other order gives 4.4 and 3.2, splitting the same 7.5 points).
The fall on the benign side is mostly a change in composition (Figure 7C).
New Benign variants are much rarer than the old ones (median population-maximum allele frequency 0.28% against 3.6%, and 31% above 0.5% against 70%), and Benign is the class that BA1 can reach by frequency alone: recall of Benign variants falls from 71.4% to 52.3%, while recall of Likely benign variants, which frequency rarely decides, rises from 51.1% to 55.7%.
At the April recall of each of nine bands of population-maximum allele frequency, the new Benign variants would have a recall of 54.9%, so the mix of bands explains 86% of the fall (16.5 of 19.1 points), and for Likely benign variants it explains most of the rise (expected 54.9%).
We made the split into nine bands after a split into six, with one band from 1e-2 up and one from 1e-3 to 1e-2, which explained 70%, had shown where the rest lay, so the 86% describes the fall and does not test it.
What remains is mostly in the highest bands (the two at 2% or more carry 1.8 of the 2.7 points that the mix does not explain): the new Benign variants at 2% or more are called benign-side for 58% (326 of 564) against 91% of 83,527 April ones, and 157 of the 238 that are not carry BP4 alone, so neither BA1 nor BS1 fired on them.
Annotating those 238 again with the details of each criterion shows why.
For 189 (79%) BA1 and BS1 were withheld by design: 124 sit in a low-complexity region that gnomAD flags, 56 at a site whose gnomAD record fails its own filter, 5 in a segmental duplication and 4 in a gene of the homologous-gene list, each counted once, by the first of these that applies (Section 3.4).
For 17 (7%) the criteria were evaluated and the filtering allele frequency is below the bar, which is how one small group can inflate a population maximum.
Another 31 met BS1 and stay at VUS, 23 with BS1 alone, which the 2015 table on the benign side does not call, and 8 with pathogenic-side evidence beside it.
One more has an allele number below the floor.
The 189 are 34% of the 564 new Benign variants of that frequency, against 4.8% (4,035 of 83,527) of the April ones, so the new variants are far more often at sites where the classifier does not trust the gnomAD frequency, and we have not established why.
Recall of ClinVar Pathogenic variants is 64.7% against 65.0%, and of Likely pathogenic variants 51.4% against 48.3%.

**The rate of opposite-direction calls is higher on new variants, and the kinds resemble the old ones.**
Among the new variants it is 0.017% (18 of 107,528; 95% interval 0.009 to 0.025% by gene-cluster bootstrap), 2.3 times the 0.0073% of the April set (0.005 to 0.010%), and the interval of the ratio, from a bootstrap that resamples genes jointly across the two sets, is 1.2 to 4.0.
The mix of ClinVar class and consequence does not explain the excess, taking the April rates as known: at the April rate of its own stratum, the new set would have had 5.8 opposite-direction calls by ClinVar direction and consequence group, 6.1 by ClinVar class and consequence group, 7.1 by class alone and 7.9 by direction alone, against the 18 observed.
The April figure was reached after every contradiction in it had been reviewed and many repaired, and the new variants had had no such attention, which would raise their rate whatever else differs.
We cannot separate that from other differences between the sets, such as how the new variants reached two stars (nine in ten were single-submitter records in April).
The 18 are of three kinds that we had named on the April residual: 9 are PVS1 on a ClinVar benign-side variant (ALMS1, CACNA1G, CCNQ, MAP3K14, ORAI1, PRKDC, SMARCA2, TRPV4, UNC93B1), 6 are BP4 and BP7 on a ClinVar pathogenic-side variant (ABCA4, HEXA, IFT122, NUP107, PKD1, SLC7A9), and 3 are missense variants called pathogenic-side on PP2 and PP3 (COL1A1 twice and KCNH1).
The kinds are coarse and were defined after the April residual was seen, and 18 events give little power to find a new one, so this shows that the new contradictions resemble the old ones and not that no other kind exists.

**The changes proposed in Section 4.5 do not fail on the new variants where they can be tried.**
BP7 as the splicing recommendation describes it adds 4,933 calls in direction on the 107,528 new variants (4.6% of them, against 2.4% of the April set) and no missed diagnosis, where on the April set it added 15,953 and 2.
The completed ClinGen file adds 47 calls in direction and no opposite-direction call, where on the April set it added 834 and 2 false-pathogenic calls; Orphanet's label adds 272 for 2 false-pathogenic calls (1,795 and 5 on the April set), and the two together 308 and 2 (2,053 and 5); the rates for Orphanet's label and for both are 136 and 154 calls per added false-pathogenic call on the new variants against 359 and 411 on the April set, and on the new variants rest on 2 events.
BP1's veto without the variant's own record changes no call, and the test is weak: the veto changes for a variant only where the gene's count is exactly 3, which is 2 of the 94 new ClinVar-pathogenic missense variants in genes where BP1 can fire (NUP85 and RPS29, whose own records are in the index), and neither has BP4 beside BP1 to make a call from.
The measurement build discounts one record whenever ClinVar's allele-level record of 27 June is pathogenic, which for up to 17 of the 94 (those whose own record the April index does not hold as pathogenic) discounts a record that is not in their count, and the run changes no call for them either.
BP1's veto at a count of 1 gives 122 more calls in direction (157 fewer ClinVar-VUS variants called benign-side and 35 fewer correct benign-side calls) and changes no opposite-direction call, because none of the 6 missed diagnoses among the new variants rests on BP1; on the April set it gave 1,716 calls and removed 2 missed diagnoses.
Removing the allele-number floor changes no call, which tests it on at most 3 variants: 33 new variants have an allele number below 2,000, and 3 of them have a population maximum at or above the BS1 bar of 5e-3.
These are seven runs on one set of variants, and the new variants lie near April variants (Section 2.8) and were mostly in ClinVar in April, so they are evidence that the changes do not fail on variants that no tuning scored and not that they would be free.

**ClinVar moved few of the calls the classifier made against a VUS, and the 157 variants that changed direction were mostly ones the classifier had left as VUS.**
ClinVar changed the direction of 157 of the 665,057 April variants in six months, 36 from VUS to a pathogenic-side class, 98 from VUS to a benign-side class, 22 into VUS (14 from pathogenic-side and 8 from benign-side) and one from benign-side to pathogenic-side.
The classifier called VUS for 119 of the 134 that left VUS, and agreed with ClinVar's new class for 30 of the 157 and with its earlier class for 124.
For the 22 moved into VUS it agreed with the new class for 17.
Of the 8,228 April variants where the classifier committed against a ClinVar VUS, 13 later moved in the direction of its call and 2 moved against it, so a test of whether ClinVar later moved toward those calls cannot be made from six months of reclassification.
The 44 residual calls that are still in the set are unchanged in class.

**What the later release does not test.**
The new variants lie in 4,124 genes, all but one of which are in the April set, so this is a hold-out on variants and not on genes, and every gene-level setting (the exception list, the homologous-gene list, the ClinGen gene tags) was set with the April genes in view.
The ClinVar-derived indices the classifier reads were not rebuilt, so the protein index lacks every record added since April, and holds the April record of each new variant that was then a classified missense variant, which for BP1's veto is the variant's own record (for 77 of the 94 new ClinVar-pathogenic missense variants in genes where BP1 can fire, and for the other 17 the April index does not hold it as pathogenic), as it is for 788 of the 790 April variants of the same kind.
The records of other variants that were added or reclassified since lower PS1 and PM5 for a new variant, lower or raise PM1 (fewer pathogenic neighbours lower it and fewer benign ones raise it), and make the veto of BP1 less likely to apply, which could add missed diagnoses.
A new variant that already has a pathogenic record in the allele-level annotation of 27 June has that record left out of the index counts of PS1 and PM1 although the April index does not hold it, which lowers both further.
None of the 18 opposite-direction calls on the new variants rests on BP1, PS1 or PM5 and one rests on PM1, so the missing records can have produced at most one of them.
And six months is short: ClinVar reclassifies slowly, so the later release adds mostly second submissions to variants that ClinVar already held, and little new knowledge about old variants.

## 4. Discussion

### 4.1 Principal findings

Five results of this work seem to us worth testing on other classifiers, which we have not done: they come from one tool, one reviewer and one benchmark.

**The aggregate stopped rewarding the work.**
Agreement in direction with ClinVar reached 77.6% at v13 and stayed between 77.50% and 77.62% for the next 23 runs, while the work of those runs cut opposite-direction calls from 94 to 49 and the ledger records 20 further changes.
A project that watched only agreement would have concluded at v13 that the work was finished, and simple rules with no combining rule or safeguard reach higher agreement with 35 times as many contradictions (Section 3.5).
The two largest steps before v13 (8.8 and 5.2 points, from figures recorded at the time) were mostly plumbing: SpliceAI and PhyloP reaching the criteria that read them, and two conditions of the guidelines that criteria had not applied (an absent variant read as absent for PM2, and BP4's confinement to splice questions).
The measures that moved late were counts of rare, severe errors, which are small numbers, and they are the ones a clinician would act on.

**Review of the discordant calls first identified the most changes, and mostly by reading the conditions a criterion's text states.**
Review first identified 30 of the 64 ledger items (47%, by our coding), and the largest class (21 of 64) was a condition of a guideline that a criterion had not applied: BS2's full penetrance at an early age, PM1's absence of benign variation, BP1's truncating-only mechanism, PVS1's need for a gene-disease relationship and for loss of function as the mechanism.
The reviewer supported ClinVar's call on 89%, 86% and 67% of the rows with a verdict in the three returned rounds (by three different instruments, Section 3.2).
A benchmark cannot find a missing condition when the omission does not change an aggregate, and several did not: building BP3's database moved agreement by 24 calls, and making the BA1 exception list match for the first time moved no cell of the matrix, because seven of the nine alleles it exists for are not in the April reference set.
Review has a sampling problem of its own, and ours was a cap of 10,000 rows in a discrepancy log that silently limited the first round to the first five chromosomes.

**A threshold belongs to the combining rule beside it.**
The PM2 bar moved pathogenic-side recall by 19 points when it was set in v10, by 0.5 in v37, and by 17.9 again when the 2015 table is restored, because the table needs PM2 to lift a null variant to Likely pathogenic and the point system does not.
The cutoffs that define the classes move the pathogenic side more than any numeric bar that we swept on the shipped classifier (11 percentage points of recall for one point of cutoff, against 0.5 for the whole range of the PM2 bar), and they are fixed by the framework and not by the data.
We showed this for the PM2 bar, the only bar that we re-ran under the 2015 table.
Thresholds that were tuned should be re-measured whenever the rule beside them changes, and the exchange rate between calls gained and errors added is a better unit than the agreement figure.
Re-measuring also turned up two thresholds that had never been swept and mattered: BP7's conservation bar, which matters most on the benign side, and the count of BP1's veto, where a stricter value calls fewer ClinVar-VUS variants benign-side and avoids two missed diagnoses, one of which rests on the variant's own record.

**Agreement with ClinVar partly reflects the classifier's own inputs, and the dependence that remains is small in agreement and larger in errors.**
Counting a variant's own ClinVar record would add 14.3 points of pathogenic-side recall, and the classifier excludes it from PS1 and PM1.
Removing every ClinVar-derived input costs 0.77 points of agreement and adds 28 missed diagnoses (24 to 52) and takes the opposite-direction calls from 49 to 74, 27 of the missed diagnoses because it removes the BP1 veto.
The veto reads other variants' ClinVar classifications to apply BP1's condition on the mechanism of the gene, so in a benchmark against ClinVar its protection is itself part of the dependence.
It also counts the variant's own record, and 3 of the 27 rest on it: with that and the two other adjustments that read the own record or followed the reviewer's rows removed together, the 24 missed diagnoses are 34 and the 49 opposite-direction calls 59.
Reporting both figures is what makes the claim checkable.

**An input can be incomplete without anything showing, and a check against an independent listing can find it.**
The ClinGen file that the benchmark read lacked 104 curations, among them CFTR, and the classifier degraded in the way its design implies: a gene with no curation is not blocked, it only has no route to loss of function, so every null variant of CFTR was called VUS.
The design is the defect worth naming, because a gene that is absent from the curated file silently has no route, and an explicit flag for a gene that is not in the file would have listed it, among 7,391 of the benchmark's 10,231 genes.
Nothing in an aggregate showed it, and the review could not, because a VUS is not an opposite-direction call; even the benchmark's four-star records, all of them CFTR, were called VUS for 22 of 23.
We found it by asking why one gene had no PVS1.
A mechanical check would have found it sooner: list the genes in which ClinVar's pathogenic truncating variants have no route to PVS1, and compare them with an independent listing.
That check flags 625 genes (576 of them not in the file), of which the completed ClinGen file explains 23 and the file with Orphanet's label 170, so it is a prompt for review and not a detector, and here one input defect sat among many gaps of design.

### 4.2 Reading the discordance

A ClinVar agreement of 63% exact and 78% in direction is not a quality score, and the decomposition of Table 5 shows why.
Almost none of the discordance is contradiction: 49 of 673,660 records (0.0073%) are opposite, and 24 of them are missed diagnoses.
Almost all the rest is a defined ClinVar class that the classifier leaves as VUS (142,115 records).
Part of that is evidence that a single-sample VCF does not carry (segregation, de novo status, phenotype specificity, functional assays, case counts), and 13 of the 28 criteria cannot fire in this benchmark.
Part is design: 10,541 of the 35,178 pathogenic-side variants called VUS sit one point below the Likely pathogenic cutoff, 5,184 are the PVS1 gap, and of the 106,937 benign-side variants called VUS 23,053 carry a lone BP4 or BS1, which the 2015 table on the benign side does not call, and 31,912 carry PM2 only.
And part is an input that was incomplete (Section 3.5).
We cannot apportion the three from the aggregate or from the VCF.
A VUS from an annotator that has these inputs is a statement about the inputs, and it does not mean what a VUS from a laboratory means.
A further 97,017 records are in the right direction at the wrong tier, and for the pathogenic-side half the boundary falls one point above where a lone null variant with PM2 lands: 41,234 of the 48,301 ClinVar Pathogenic variants called Likely pathogenic are at 9 points.
The classifier's calls are therefore best read with their error profile beside them: of 62,247 pathogenic-side calls 94.2% are pathogenic-side in ClinVar, 5.8% are ClinVar VUS and 0.04% are benign-side in ClinVar, and none of the 4,017 Pathogenic and the 67,700 Benign calls contradicts ClinVar in this in-sample set.

The residual contradictions are informative because they are few and each was reviewed by one medical geneticist, whose rulings are one expert's reading and by our coding.
In the last round the reviewer supported ClinVar on 61% of the 49 and fastVEP on 8%, with a different class or more evidence needed for the rest, and the reasons they gave were rarely about the arithmetic of the criteria.
They were the frequency a disease model allows, the phase of a variant, a truncation that escapes NMD in the first exon, an isoform, a transcript that predates MANE, and a recessive disorder.
Each is evidence or a policy that a single-sample, single-transcript annotator can only approximate, and several can be approximated better (a disease-specific frequency bar, a first-exon rule, a transcript policy).
We have not done so here, because each is a decision about the classifier and not about this paper's measurements, and because the residual is the sample on which the classifier was reviewed: a fix tuned on it could not be tested on it; the hold-out of Section 3.7 is the test that the data allow, and the changes of Section 4.5 that could be run there do not fail on variants the review never saw.

ClinVar is a reference and not ground truth, and three observations in this work bear on how far to trust it.
Five of the 49 residual variants (10%) had left the 2-star+ set by October, all with conflicting classifications, against 1.3% of all April variants that left it for any reason.
ClinVar moved 157 of 665,057 variants across a direction in six months.
And the label of a variant can depend on the haplotype it lies on: for BRCA1 c.594-2A>C the reviewer's note records that the variant is always in cis with c.641A>G and that the haplotype, not the variant, has been shown not to be pathogenic, which no single-variant annotator could reproduce.
A discordance with ClinVar is a prompt to look, and the reviewer's reasons suggest that it is often a prompt about evidence and not about the arithmetic of the criteria, which one reviewer's reading of 49 rows cannot show for the other discordance.

### 4.3 Practical lessons

Several lessons are of practical use to anyone building a classifier of this kind.
Test a criterion with the pipeline's own output and not a hand-written input: the BA1 exception matcher passed its unit tests because they fed it the bare token that the pipeline never produces.
Check that a benchmark can see what is being tested, because it cannot validate a criterion whose population it excludes.
Keep a leave-one-out mode and report the ClinVar-informed figure beside it, so that a reader can see the size of the circularity, and apply it to every ClinVar-derived input, which we did not do for BP1's veto.
Report a baseline of simple rules, and the contradictions, beside the agreement figure, because agreement can be exceeded by calling aggressively.
Record what produced every annotation, and refuse to reuse one that something else produced, because a cache that was right for an older parser is intact and is not safe to read, as the repository's own transcript cache showed when a parser was fixed (issues 95 and 98).
Check each input against an independent listing, because a missing record degrades quietly.
Count severe errors as a primary measure beside agreement, and say which errors weigh more: by that count v10 (25 opposite-direction calls) was better than v37 (49), and it was given up for 6,632 correct pathogenic-side calls, which is a trade and not an improvement.
And give the reviewer the whole row: the first workbook carried 75 columns of evidence, and the list of discordance categories that the reviewer attached to round 1 anticipates three classes of the ledger.

### 4.4 Limitations

The study has limits that the reader should weigh.
A single medical geneticist reviewed the calls, was not blind to either call, supplied inputs as well as rulings, saw only the opposite-direction calls and none of the 35,178 pathogenic-side variants that we left as VUS or any concordant call, and returned a marked workbook in three of six rounds; the round-6 categories are our coding of their text, and we have not asked the reviewer to confirm our coding of rounds 2 and 6.
The ledger of 64 items and its classes are our own coding of the project's record and have not been independently coded.
All thresholds were tuned on the April benchmark and its discordances, so the April figures are in-sample (the 49 are 59 when the three adjustments that read the own record or followed the reviewer's rows are removed), and the October hold-out is a hold-out of variants and not of genes, of a release six months later, whose new variants were mostly in ClinVar in April at a lower review status.
The price of an error was not applied uniformly: the trades accepted at v11 and v13 cost 214 and 315 correct calls per added error (on 31 and 36 errors), and at the BS2 and BS1 bars we declined trades at 463 and 928 (on 3 and 6 errors), each decision resting on its own recorded reasons.
The ClinVar-derived inputs of the classifier were built from releases of 29 April and 27 June 2026 and the reference is the VCF of 4 April, so the three are not one snapshot.
The reference is a selection (two stars and above, SNVs and indels of at most 100 bases): it excludes records with conflicting classifications (six of the nine BA1-exception alleles of Ghosh et al. are among them in later releases, and a seventh, HFE, was absent from April and is in October) and records with a single submitter.
It also has no negative controls, and its mix of classes (25% pathogenic-side among the defined classes) fixes the precision figures: a different mix would give different precision for the same calls.
A variant's call is the call on one picked transcript, and other transcripts can differ.
The classifier evaluates 15 of 28 criteria on a single-sample VCF without phenotype, so its recall of ClinVar's classes has a ceiling that is not a fault of the thresholds.
The classifier reads the gnomAD exomes only, and PM2 reads absence from them as rarity, including at sites that the exomes do not cover; requiring a record costs 0.22 points of agreement and 2.9 points of pathogenic-side recall (Figure 8A), and PM2 fires on 64% of the variants.
The exchange rates of Section 3.3 are counts of differences between deterministic runs and carry no sampling interval: a rate computed from one to three added errors is as uncertain as a count of that size (the 95% Poisson interval of a count of 2 is 0.2 to 7.2).
The offline replay of the benign-side point scales was not checked by real runs, and the cutoff and measurement builds are patched builds of the shipped code.
We made no comparison with another tool, and the bulk outputs of the runs (annotated VCFs of 250 to 285 MB and per-variant tables of about 20 MB, compressed) are not archived with the repository; its scored matrices and provenance are.
The gene inheritance mode that BS2 and PM2 read is the weakest input: the classifier reads it from the names of the diseases, so that of the 2,924 genes of the ClinGen file it takes 67 as recessive-only, 50 as dominant-only and 12 as both, and of the 740 genes with a benchmark variant whose established ClinGen tags are dominant only, 39 register as dominant, and of the 1,285 whose established tags are recessive only, 60 register as recessive.
Making BS2 read the tag gains 2,812 calls in direction and moves six ClinVar pathogenic-side variants toward benign, so it was measured and not shipped (RUN_VERSIONS.md), and the flat sweeps of the recessive-only PM2 bar and the dominant-only BS2 bar (Section 3.3) say nothing about what those bars would do on a correct field.
Where the project's documentation attributed a number or a grade to a published recommendation that we could not find stated there, or where the implementation departs from the recommendation, we say so in Sections 3.2 and 3.3.
fastVEP's classification is a decision-support output for research use: it is not a laboratory-validated classification and does not replace expert review.

### 4.5 Next steps

The findings suggest an order of work, and the measurements of Sections 3.3 to 3.5 give each step a price; Section 3.7 tried the changes that could be run on the 107,528 new variants.
First, complete the ClinGen input, which adds 852 correct pathogenic-side calls (834 in direction) for 2 false-pathogenic calls on the April set and 51 (47 in direction) with no opposite-direction call on the new variants; the check of Section 4.1 and a flag for a gene that is not in the curated file were not measured, and the flag would list 7,391 of the benchmark's genes.
Second, a source of established loss-of-function mechanism that does not depend on constraint (the PVS1 applicability stated in VCEP specifications, or Orphanet's label, which adds 272 calls in direction for 2 false-pathogenic calls on the new variants), together with handling of genes that have several disorders, because 3,470 ClinVar truncating variants would still carry no PVS1 after both inputs are used.
Third, BP7 as the splicing recommendation describes it, which gains 15,953 calls in direction for 2 missed diagnoses on the April set and 4,933 for none on the new variants, and which only the maintainers can weigh.
Fourth, BP1's veto leaving out the variant's own record (3 more missed diagnoses on the April set and none on the new variants) with a stricter count: the pair at a count of 1 gives 1,716 more calls in direction and 149 fewer correct benign-side calls on the April set, with 23 missed diagnoses against v37's 24, and the count of 1 alone gives 122 more calls on the new variants and avoids no missed diagnosis there.
Fifth, a transcript policy for MANE Plus Clinical, a first-exon rule for PVS1, and disease-specific frequency bars for BA1 and BS1 where a VCEP has published one, which the classifier already supports as an opt-in table.
Sixth, review of the allele-number floor, whose source is narrower than its use and whose removal would add 774 calls in direction on the April set and none on the new variants, and of the two grades and the NMD default of PVS1 (Section 3.2).
Each of these can be measured by the harness of this paper before it is shipped, and none should be tuned to the 49, which is why the benefit of a stricter count of BP1's veto, which avoids two of them, is in-sample.

## 5. Conclusions

Building a software classifier for the ACMG/AMP standards was largely the work of applying conditions that the standards state and a first implementation omitted (21 of 64 changes), and review of the discordant calls first identified more of the changes (30, by our coding) than any other route.
On the 673,660 ClinVar 2-star+ variants on which every threshold was tuned, the final version agrees with ClinVar's direction for 77.6% and contradicts it for 49 (0.0073% of the records and 0.013% of those with a pathogenic-side or benign-side class), which are 59 under a strict leave-one-out reading, three of the added calls Benign; none of its 4,017 Pathogenic or 67,700 Benign calls is among the 49, and rules with no combining rule or safeguard reach higher agreement with 35 times as many contradictions.
On the 107,528 variants that were not in the April set, agreement in direction is 70.1%, mostly because of the benign side, and contradictions are 2.3 times as frequent (interval 1.2 to 4.0), more than their mix predicts at the April rates, and of the same coarse kinds.
The thresholds that no guideline states are less decisive than the combining rule around them for the PM2 bar, the one we tested under the 2015 table, and two that had never been swept, BP7's conservation test and the count of BP1's veto, each exchange calls and errors at a rate that the maintainers can now weigh, which we have not changed; on the pathogenic side the points that define the classes move the result more than the bars that we swept on the shipped classifier.
The remaining discordance is consistent with evidence that a single-sample VCF does not carry and with design choices of the classifier, among them PVS1's condition on mechanism, which is approximated by gene-disease validity or constraint.
The largest gap, in PVS1, is in part an input of unknown cause that was incomplete (CFTR and 103 other curations, which explain 13% of the gap) and in part a design gap, of which Orphanet's loss-of-function label closes another 20%.

---

## References

Abou Tayoun, A. N., Pesaran, T., DiStefano, M. T., et al. (2018). Recommendations for interpreting the loss of function PVS1 ACMG/AMP variant criterion. *Human Mutation*, 39(11), 1517-1524. https://doi.org/10.1002/humu.23626

Biesecker, L. G., & Harrison, S. M., for the ClinGen Sequence Variant Interpretation Working Group. (2018). The ACMG/AMP reputable source criteria for the interpretation of sequence variants. *Genetics in Medicine*, 20(12), 1687-1688. https://doi.org/10.1038/gim.2018.42

Brnich, S. E., Abou Tayoun, A. N., Couch, F. J., et al. (2019). Recommendations for application of the functional evidence PS3/BS3 criterion using the ACMG/AMP sequence variant interpretation framework. *Genome Medicine*, 12(1), 3. https://doi.org/10.1186/s13073-019-0690-2

ClinGen. (2024). ClinGen's guidance: classifying variants in genes associated with multiple disorders, version 1 (January 2024). https://clinicalgenome.org/docs/clingens-guidance-classifying-variants-in-genes-associated-with-multiple-disorders/

ClinGen. (2026). Gene-disease validity standard operating procedures, version 12 (2 February 2026). https://clinicalgenome.org/docs/gene-disease-validity-standard-operating-procedures-version-12/

ClinGen. (n.d.). Criteria Specification (CSpec) Registry. https://cspec.genome.network

ClinGen SVI (ClinGen Sequence Variant Interpretation Working Group). (2019). SVI recommendation for in trans criterion (PM3), version 1.0 (approved 2 May 2019). https://clinicalgenome.org/docs/pm3-recommendation-for-in-trans-criterion-pm3-version-1.0/

ClinGen SVI (ClinGen Sequence Variant Interpretation Working Group). (2020). SVI recommendation for absence/rarity (PM2), version 1.0 (approved 4 September 2020). https://clinicalgenome.org/docs/pm2-recommendation-for-absence-rarity/

ClinGen VCEP Review Committee. (2024). ClinGen SVI VCEP Review Committee communication to ClinGen VCEPs regarding gnomAD v4, version 1.0 (8 March 2024). Later versions (August 2024 and June 2025) keep the recommendation to use the grpmax filtering allele frequency for BA1 and BS1. https://clinicalgenome.org/docs/clingen-guidance-to-vceps-regarding-the-use-of-gnomad-v4/

Clopper, C. J., & Pearson, E. S. (1934). The use of confidence or fiducial limits illustrated in the case of the binomial. *Biometrika*, 26(4), 404-413. https://doi.org/10.1093/biomet/26.4.404

gnomAD consortium. (n.d.). gnomAD v4.1 (the Genome Aggregation Database). https://gnomad.broadinstitute.org

Ghosh, R., Harrison, S. M., Rehm, H. L., Plon, S. E., & Biesecker, L. G., for the ClinGen Sequence Variant Interpretation Working Group. (2018). Updated recommendation for the benign stand-alone ACMG/AMP criterion. *Human Mutation*, 39(11), 1525-1530. https://doi.org/10.1002/humu.23642

Harrison, S. M., Dolinsky, J. S., Knight Johnson, A. E., et al. (2017). Clinical laboratories collaborate to resolve differences in variant interpretations submitted to ClinVar. *Genetics in Medicine*, 19(10), 1096-1104. https://doi.org/10.1038/gim.2017.14

Ioannidis, N. M., Rothstein, J. H., Pejaver, V., et al. (2016). REVEL: an ensemble method for predicting the pathogenicity of rare missense variants. *American Journal of Human Genetics*, 99(4), 877-885. https://doi.org/10.1016/j.ajhg.2016.08.016

Jaganathan, K., Kyriazopoulou Panagiotopoulou, S., McRae, J. F., et al. (2019). Predicting splicing from primary sequence with deep learning. *Cell*, 176(3), 535-548. https://doi.org/10.1016/j.cell.2018.12.015

Karczewski, K. J., Francioli, L. C., Tiao, G., et al. (2020). The mutational constraint spectrum quantified from variation in 141,456 humans. *Nature*, 581(7809), 434-443. https://doi.org/10.1038/s41586-020-2308-7

Kopanos, C., Tsiolkas, V., Kouris, A., et al. (2019). VarSome: the human genomic variant search engine. *Bioinformatics*, 35(11), 1978-1980. https://doi.org/10.1093/bioinformatics/bty897

Landrum, M. J., Lee, J. M., Benson, M., et al. (2018). ClinVar: improving access to variant interpretations and supporting evidence. *Nucleic Acids Research*, 46(D1), D1062-D1067. https://doi.org/10.1093/nar/gkx1153

Li, Q., & Wang, K. (2017). InterVar: clinical interpretation of genetic variants by the 2015 ACMG-AMP guidelines. *American Journal of Human Genetics*, 100(2), 267-280. https://doi.org/10.1016/j.ajhg.2017.01.004

Mandelker, D., Schmidt, R. J., Ankala, A., et al. (2016). Navigating highly homologous genes in a molecular diagnostic setting: a resource for clinical next-generation sequencing. *Genetics in Medicine*, 18(12), 1282-1289. https://doi.org/10.1038/gim.2016.58

McLaren, W., Gil, L., Hunt, S. E., et al. (2016). The Ensembl Variant Effect Predictor. *Genome Biology*, 17(1), 122. https://doi.org/10.1186/s13059-016-0974-4

Morales, J., Pujar, S., Loveland, J. E., et al. (2022). A joint NCBI and EMBL-EBI transcript set for clinical genomics and research. *Nature*, 604(7905), 310-315. https://doi.org/10.1038/s41586-022-04558-8

Nykamp, K., Anderson, M., Powers, M., et al. (2017). Sherloc: a comprehensive refinement of the ACMG-AMP variant classification criteria. *Genetics in Medicine*, 19(10), 1105-1117. https://doi.org/10.1038/gim.2017.37

Orphanet. (2026). Orphadata, product 6: gene-disorder associations with their association types (release of 23 June 2026; CC BY 4.0). INSERM. https://www.orphadata.com

Pejaver, V., Byrne, A. B., Feng, B. J., et al. (2022). Calibration of computational tools for missense variant pathogenicity classification and ClinGen recommendations for PP3/BP4 criteria. *American Journal of Human Genetics*, 109(12), 2163-2177. https://doi.org/10.1016/j.ajhg.2022.10.013

Pollard, K. S., Hubisz, M. J., Rosenbloom, K. R., & Siepel, A. (2010). Detection of nonneutral substitution rates on mammalian phylogenies. *Genome Research*, 20(1), 110-121. https://doi.org/10.1101/gr.097857.109

Rehm, H. L., Berg, J. S., Brooks, L. D., et al. (2015). ClinGen - the Clinical Genome Resource. *New England Journal of Medicine*, 372(23), 2235-2242. https://doi.org/10.1056/NEJMsr1406261

Richards, S., Aziz, N., Bale, S., et al. (2015). Standards and guidelines for the interpretation of sequence variants: a joint consensus recommendation of the American College of Medical Genetics and Genomics and the Association for Molecular Pathology. *Genetics in Medicine*, 17(5), 405-424. https://doi.org/10.1038/gim.2015.30

Smit, A. F. A., Hubley, R., & Green, P. (n.d.). RepeatMasker. http://www.repeatmasker.org (software; the authors ask to be cited by name and web page).

Strande, N. T., Riggs, E. R., Buchanan, A. H., et al. (2017). Evaluating the clinical validity of gene-disease associations: an evidence-based framework developed by the Clinical Genome Resource. *American Journal of Human Genetics*, 100(6), 895-906. https://doi.org/10.1016/j.ajhg.2017.04.015

Tavtigian, S. V., Greenblatt, M. S., Harrison, S. M., et al. (2018). Modeling the ACMG/AMP variant classification guidelines as a Bayesian classification framework. *Genetics in Medicine*, 20(9), 1054-1060. https://doi.org/10.1038/gim.2017.210

Tavtigian, S. V., Harrison, S. M., Boucher, K. M., & Biesecker, L. G. (2020). Fitting a naturally scaled point system to the ACMG/AMP variant classification guidelines. *Human Mutation*, 41(10), 1734-1737. https://doi.org/10.1002/humu.24088

Walker, L. C., de la Hoya, M., Wiggins, G. A. R., et al. (2023). Using the ACMG/AMP framework to capture evidence related to predicted and observed impact on splicing: recommendations from the ClinGen SVI Splicing Subgroup. *American Journal of Human Genetics*, 110(7), 1046-1067. https://doi.org/10.1016/j.ajhg.2023.06.002

Whiffin, N., Minikel, E., Walsh, R., et al. (2017). Using high-resolution variant frequencies to empower clinical genome interpretation. *Genetics in Medicine*, 19(10), 1151-1158. https://doi.org/10.1038/gim.2017.26

---

## Data and code availability

fastVEP is open source under the Apache License 2.0 at https://github.com/Huang-lab/fastVEP.
The analysis of this paper is in the repository under `analysis/acmg_benchmark/` (the benchmark harness `run_benchmark.sh`, the run record `RUN_VERSIONS.md`, and the scripts `11_manuscript_runs.py`, `12_per_variant_table.py`, `13_manuscript_stats.py`, `14_temporal_holdout.py` and `15_clingen_complete_export.py`) and under `manuscript/acmg/` (this manuscript, the tables it reads in `data/`, `build_tables.py`, `verify_numbers.py`, `generate_figures.py`, the figures, and the patches `replay_validation/combiner_cutoff_env.patch` and `replay_validation/bp_hooks_env.patch` that expose the two point cutoffs, and BP1's own-record count and BP7's SpliceAI bar, for the runs of Table 4 and of Sections 3.3 and 3.4).
The scored concordance matrix, the criterion firing counts and the provenance file of each run made for this paper (89 on the April set and 7 on the October release) are in `analysis/acmg_benchmark/results/v37_manuscript/`; the provenance records the flags, the SHA-256 of the binary and of every configuration, the inputs read and the commit of the checkout that ran the harness (5d8004a for 63 of the 89 runs, and 1dfa488 for 26, the head of the branch that was squash-merged as 5d8004a, whose tree is the same).
The annotated VCFs (about 250 MB each), the per-variant tables, the annotation databases and the ClinVar, gnomAD, ClinGen and Ensembl inputs are not in the repository.
Table 2 and Section 2.2 give their sources and versions, and the ClinVar files are the weekly releases dated 4 April and 4 October 2026.
Orphanet's association table is the public Orphadata product 6 (`build_gene_disease_attributes.py` gives the command).
The ClinGen export that completes the file of 30 April is a public download of 11 October 2026 (`15_clingen_complete_export.py` gives the command), and the 112 rows it adds, 104 of them curations that the April file lacked, are listed with the status of each in `manuscript/acmg/data/clingen_missing_curations.csv`.
The reviewer's workbooks are not in the repository.

## Declarations

**Author contributions.** [To be completed by the author.]

**Use of AI assistance.** [To be completed by the author.
The analysis scripts, the figures and the first draft of this manuscript were produced with the assistance of Claude (Anthropic) working as a coding agent in the project repository, and the author has not yet reviewed them in full.]

**Funding.** [To be completed by the author.]

**Competing interests.** [To be completed by the author.]

**Acknowledgements.** The medical geneticist who reviewed the calls over six rounds [name and affiliation to be added, with their consent].

**Ethics.** The study uses public records of variants and no patient data, and no human participants were recruited.
The reviewer's workbooks are not public, and the manuscript paraphrases the reviewer's rulings; whether they may be quoted is for the author to settle with the reviewer.

---

## Figure legends

**Figure 1. Agreement with ClinVar 2-star+ classifications across 37 benchmark runs.**
(A) Same-direction and exact agreement as a percentage of 673,660 variants.
Open markers are v1, a matrix reconstructed from the documentation of the run, and v2 and v4 to v6, which were recorded at the time; no raw output of them was kept.
v1 annotated 627,376 of the records, and its 54.7% of those is 51.0% of all records.
The cross marks v17, the opt-in run with the published VCEP bars, which is not a step on the way and is left out of the lines.
(B) Opposite-direction calls on a log scale: ClinVar pathogenic-side variants called benign-side (missed diagnoses), ClinVar benign-side variants called pathogenic-side, and their total; open markers as in A, and v6 is the recorded total.
(C) Same-direction recall for ClinVar's pathogenic-side and benign-side variants; the values of v1 (10.9% and 19.0%) are below the axis.
The 13 numbered markers are the changes listed under the panels.
Run v3 was a partial run that was overwritten, v21 and v22 repeat v20, and v24 and v25 were intermediate runs.
ClinVar's classifications are a reference and not ground truth.

**Figure 2. What the reviewer ruled.**
(A) The rulings on the rows of the three workbooks the reviewer returned, in rounds 1, 2 and 6, in five categories: ClinVar's call supported, fastVEP's call supported, neither (the reviewer's call is VUS or another class), more evidence needed or dependent on the disorder, and no verdict recorded.
Rounds 1 and 2 are decoded from cell colours against each workbook's key, and round 6 is our coding of the reviewer's text (Supplementary Table S7).
(B) The 49 opposite-direction calls of v37, by the evidence behind the call, with the rulings of the last round.

**Figure 3. The 64 recorded changes by class and by the route that first found them.**
Each cell counts changes in the ledger of Supplementary Table S2.
The column on the right gives the total of each class, and the row beneath the cells the total of each route.
Classes and routes are our coding of the project's record.

**Figure 4. Thresholds re-measured on the final classifier.**
(A) Five thresholds on the benign side trade correct benign-side calls for missed diagnoses (the pink line is the count of BP1's veto, from 1 to 25 pathogenic missense variants); open markers are swept values and the filled black marker is the shipped one.
(B) Recall of ClinVar's pathogenic-side variants against the PM2 bar for genes not read as recessive-only, under the SVI point system (shipped) and under the 2015 combining table, and against the bar for genes read as recessive-only under the point system.
(C) The BS2 carrier bar for genes read as dominant-only: the change in recall of pathogenic-side (red) and benign-side (blue) variants from the shipped bar.

**Figure 5. The points that define the classes.**
(A) For variants with no benign evidence, the ClinVar classes of the variants at each number of pathogenic points; the lines mark Likely pathogenic from 6 points and Pathogenic from 10.
(B) Real runs with the Likely pathogenic cutoff changed from 4 to 9 points: recall of pathogenic-side variants against the number of ClinVar VUS and benign-side variants called pathogenic-side.
(C) The benign side under the 2015 table (shipped) and under a point scale with Likely benign from -1, -2 and -3 points, replayed from the recorded criteria and not real runs.
In B and C a filled black marker is the shipped setting.

**Figure 6. The anatomy of the discordance of v37 with ClinVar 2-star+ classifications.**
(A) The five-tier concordance matrix: counts, and the share of its ClinVar class (the row) that each cell holds; the outlined cells are opposite-direction calls.
(B) What v37 says of each ClinVar class, as P/LP, VUS and B/LB; 382 records with no call (0.13% of the ClinVar benign-side class) are not drawn.
(C) What ClinVar says behind each v37 call.

**Figure 7. Discordance by consequence class and review status, and the later release.**
(A, B) The share of ClinVar's pathogenic-side and benign-side variants that v37 calls P/LP, VUS and B/LB, by consequence class and by review status; an asterisk marks a row whose bar stops short of 100% because records have no call.
(C) Recall in direction for each ClinVar class on the April set, on the October variants whose class did not change, and on the new October variants that were not in the April set; whiskers are exact binomial 95% intervals, which treat variants as independent.

**Figure 8. What each design decision and each evidence source buys, and what each departure from the sources would change.**
(A) Sixteen settings changed one at a time (switched off, or on for the two opt-in settings): the change from v37 in agreement in direction (percentage points of all records), in recall of pathogenic-side variants, in missed diagnoses and in false-pathogenic calls.
(B) Nine evidence sources removed one at a time from the annotation directory.
(C) Eight runs that test a departure from the sources or an adjustment: BP1's veto leaving out the variant's own ClinVar record, BP7's SpliceAI bar and conservation test (alone and together), the ClinGen file with the curations that its export lacked, PVS1 licensed by Orphanet's loss-of-function label (alone and with the ClinGen file), and the own-record veto, the low-penetrance gate and the exception list removed together.
The BP1 run, the SpliceAI run, the run with both BP7 changes and the run with the three adjustments use the measurement build of Section 2.6.

**Figure 9. Agreement in direction against opposite-direction calls.**
Each point is a classifier, on 673,660 variants: the versions v7 to v37 with stored scored matrices (blue, joined in order), the runs of this paper and the earlier runs of the dominant-gene BS2 settings, each v37 with one change (grey), every variant called VUS (the diamond at zero), and three rules with no combining rule, gene-level gate or safeguard at nine frequency bars each, with the second rule alone and the second with the third also given with no frequency rule (the filled markers).
R1 calls a variant benign-side when its population-maximum allele frequency is at or above the bar (the bars at the two ends of each line are labelled, and the point at 0.001 of the last rule is annotated), R2 adds a pathogenic-side call for a stop-gained, frameshift, canonical splice or start-lost variant that is absent from gnomAD or below 1e-4, and R3 adds a benign-side call for a quiet variant (synonymous, intronic, untranslated, flanking or non-coding) with SpliceAI below 0.2.
The horizontal axis is on a log scale above 10.

---

## Supplementary materials

All supplementary tables are CSV or TSV files in `manuscript/acmg/data/` of the repository, and each is written by `analysis/acmg_benchmark/scripts/13_manuscript_stats.py`, `11_manuscript_runs.py` or `15_clingen_complete_export.py` from the scored runs.

**Supplementary Table S1.** The 89 runs made for this paper and the 8 earlier runs of the dominant-gene BS2 settings: concordance matrix, firing counts and settings (`runs_v37_variants.csv`; the run matrix, the matrices of the seven runs on the October release and the provenance are in `analysis/acmg_benchmark/results/v37_manuscript/`), the comparison of the controls of the two measurement builds with v37 (`controls_v37.csv`), the calls that became opposite-direction in the runs that test BP1's veto and BP7 (`hook_new_opposite_v37.csv`), the calls that a stricter count of BP1's veto removes, and the records of the protein index behind them (`bp1_veto_count_resolved_v37.csv`, `bp1_index_counts_v37.csv`), and the binary, commit and annotation time of each run (`provenance_v37.csv`).

**Supplementary Table S2.** The ledger of 64 changes: run, class, route of discovery and source in the project record (`change_ledger.csv`).

**Supplementary Table S3.** What each step of each threshold sweep buys and costs (`threshold_steps_v37.csv`), and the published VCEP bars (`vcep_bars_v37.csv`, `vcep_audit_v37.csv`, `vcep_audit_pm2_v37.csv`).

**Supplementary Table S4.** The offline replay of the combining rules under alternative schemes and cutoffs, its validation against real runs, and the number of calls that the replay cannot reproduce because the site-level frequency veto holds them at VUS (`combining_rules_replay_v37.csv`, `replay_validation_v37.csv`, `pm2_dependence_v37.csv`, `full_points_false_benign_signatures_v37.csv`, `replay_veto_v37.csv`).

**Supplementary Table S5.** Concordance of v37 by consequence class and by review status, with the 23 four-star records (`confusion_5x5_v37.csv`, `confusion_3x3_v37.csv`, `what_is_behind_each_call_v37.csv`, `by_consequence_v37.csv`, `by_review_status_v37.csv`, `four_star_variants_v37.csv`).

**Supplementary Table S6.** The PVS1 gap: truncating variants without PVS1, by gene and by the reason PVS1 is absent (`pvs1_gap_v37.csv`, `pvs1_gap_genes_v37.csv`, `pvs1_gap_reasons_v37.csv`, `pvs1_gap_constrained_v37.csv`); the 112 rows that complete the ClinGen export of 30 April, 104 of them curations that it lacked (`clingen_missing_curations.csv`, `clingen_missing_summary_v37.csv`); what the three runs that change an input do to the gap (`pvs1_gap_closure_v37.csv`, `pvs1_gap_closure_genes_v37.csv`, `pvs1_gap_closure_new_false_pathogenic_v37.csv`); and the composition of the pathogenic-side variants called VUS (`plp_to_vus_points_v37.csv`, `plp_to_vus_by_consequence_v37.csv`, `plp_to_vus_gene_concentration_v37.csv`, `signatures_offdiagonal_v37.csv`).

**Supplementary Table S7.** The 49 opposite-direction calls of v37 with their criteria, gnomAD counts and rulings (`opposite_direction_v37.csv`, `opposite_direction_by_group_v37.csv`), the rulings of the three returned workbooks (`review_rounds.csv`), the 64 rows of the last round and what became of them (`round6_rows_v37.csv`), and the rulings with our coding of them (`round6_rulings.tsv`, `round6_rules_v37.csv`).
We coded the reviewer's round-6 text by five rules, which a row can satisfy in more than one way; the column `rule` of `round6_rulings.tsv` records the one that decided it.
(1, class, 29 rows) The text states a class: it is compared with the two calls and coded ClinVar when its direction is ClinVar's, fastVEP when it is fastVEP's, and neither when it is VUS; a pair of classes that includes VUS (one row) was coded by its other member.
(2, stated, 13 rows) The text says who is right, or that both calls are discordant with the reviewer's: coded as stated.
(3, disorder, 2 rows) The text makes the class depend on the disorder or the inheritance, in a gene with both dominant and recessive disease: coded depends on the disorder, whatever class it also states.
(4, review, 6 rows) The text asks for more review, or gives a class only on a condition: coded needs more evidence.
(5, evidence, 14 rows) The text gives evidence for or against a call without a class or a statement of who is right: coded for the call that the evidence supports (13 ClinVar and 1 fastVEP).
Of the 49 rows that are still opposite, 11 of the 30 coded ClinVar come from rule 1, 9 from rule 2 and 10 from rule 5.
The reviewer's text is not reproduced here.

**Supplementary Table S8.** The later ClinVar release: strata, class-wise recall, coverage, opposite-direction rates, the ratio of the rates with its interval and the expected numbers at the April rates, recall standardised for the allele-frequency mix, the new Benign variants of 2% or more and why BA1 and BS1 did not fire on those that were not called benign-side, agreement in direction standardised for the mix of ClinVar classes, how the new variants stand in the April protein index, the seven proposed changes run on the new variants, the 157 variants that changed direction and what became of the 49 residual calls, and the number of genes that the new variants fall in, with the criteria of their opposite-direction calls (`temporal_v37.csv`, `temporal_by_class_v37.csv`, `temporal_coverage_v37.csv`, `temporal_opposite_rates_v37.csv`, `temporal_ratio_v37.csv`, `temporal_standardised_recall_v37.csv`, `temporal_recall_by_af_band_v37.csv`, `temporal_direction_changes_v37.csv`, `temporal_direction_change_variants_v37.csv`, `temporal_new_opposite_v37.csv`, `april_opposite_in_october_v37.csv`, `temporal_new_composition_v37.csv`, `temporal_new_top_genes_v37.csv`, `temporal_new_genes_v37.csv`, `temporal_high_af_benign_v37.csv`, `temporal_high_af_blockers_v37.csv`, `temporal_agreement_mix_v37.csv`, `temporal_new_in_index_v37.csv`, `temporal_proposals_v37.csv`).

**Supplementary Table S9.** Gene-level concordance of v37 (`gene_level_v37.csv`, `plp_recall_by_gene_v37.csv`) and the headline rates with their intervals (`headline_v37.csv`).

**Supplementary Table S10.** The runs v1 to v37: agreement, opposite-direction calls and recall of each, with the number of variants each run committed to, the denominator of each, and the marker of v17 (`trajectory.csv`; v1 is a reconstructed matrix and v2 and v4 to v6 are the figures recorded at the time).

**Supplementary Table S11.** The ClinVar classes of the variants at each number of pathogenic points, and the posterior probability of pathogenicity that the Bayesian framework assigns to each number of points (`points_by_clinvar_v37.csv`, `tavtigian_posterior.csv`), and the aggregates that the text quotes from the per-variant table (`facts_v37.csv`).

**Supplementary Table S12.** How often each of the 15 computable criteria fires in v37, at which strengths, and on variants of which ClinVar class (`criteria_firing_v37.csv`).

**Supplementary Table S13.** The rules with no combining rule, gene-level gate or safeguard and what they agree with ClinVar on, at one bar and at nine, beside v37, the second rule against the PVS1 gap, and what each tier of the v37 call is worth, with the rate of opposite-direction calls and its interval (`baselines_v37.csv`, `baselines_curve_v37.csv`, `baselines_r2_vs_v37.csv`, `precision_by_tier_v37.csv`).
