# Iterative incorporation of ACMG/AMP and ClinGen SVI guidance into an open-source variant classifier: expert review, how thresholds were chosen, and discordance with ClinVar 2-star+ classifications

**Mythily Ganapathi**^1^, **Wendy Chung**^2^, **Kuan-lin Huang**^3,*^

^1^ Department of Pathology and Cell Biology, Columbia University Irving Medical Center, New York, NY, USA

^2^ Department of Pediatrics, Boston Children's Hospital and Harvard Medical School, Boston, MA, USA

^3^ Department of Genetics and Genomic Sciences, Precision Omics Laboratory, Icahn School of Medicine at Mount Sinai, New York, NY, USA

\* Corresponding author: kuan-lin.huang@mssm.edu. Website: https://PrecisionOmics.org

Running title: ACMG/AMP classification in fastVEP

---

## Abstract

**Background.**
The ACMG/AMP criteria are underspecified (no guideline states many of the thresholds that software must set, and the evidence can be combined by a table or by points that do not give the same call for every combination) and fragmented (the standards of 2015, a long series of ClinGen recommendations that each change criteria that looked settled, and gene-specific specifications).
Our contribution is how to implement such criteria in software so that the calls agree closely with community opinion, taken here as ClinVar's classifications of review status two stars or better (a reference and not ground truth) and the review of one medical geneticist, and what each choice costs in contradictions, among the settings that we measured.
We describe the classification engine of fastVEP, an open-source reimplementation of the Ensembl Variant Effect Predictor, in which the choices that the guidelines leave open were measured against that opinion.

**Methods.**
Thirty-eight numbered versions were scored on 673,660 ClinVar variants (release of 4 April 2026; 32 have stored scored output) over six rounds of review, three of which returned rulings.
The input lists variants only, with no genotype, pedigree or phenotype of any individual (a sites-only VCF), and from it 15 of the 28 ACMG/AMP criteria can be evaluated.
The criteria that are met are combined by a point scale on the pathogenic side (1, 2, 4 and 8 points for Supporting, Moderate, Strong and Very Strong evidence; 6 to 9 points is Likely pathogenic and 10 or more is Pathogenic) and by the 2015 combining table on the benign side.
The latest version (v38) reads the mode of inheritance that ClinGen curates, does not count BS2's carriers beside BS1 or BA1 and classifies a gene once for each of its modes of inheritance, where v37, the version that the review rounds ended with, read the words of the disease name; both are versions of the classifier, and every run is on GRCh38.
For this paper 89 runs of v38 on the April set re-measured all but four groups of the thresholds that no guideline states and tested departures from the sources, and 7 more applied the changes that Section 4.5 proposes to a later release; the same 96 runs of v37 give the before of every difference that v38 makes.
v38 was applied unchanged to ClinVar's release of 4 October 2026, including 107,528 variants that were not in the April set.

**Results.**
On the April set, where every threshold was tuned, agreement in direction rose from 51.0% (v1, reconstructed) to 77.8% (v38), flat from the thirteenth version to v37 (77.50% to 77.62%) and 0.22 points higher at v38, while opposite-direction calls went from 25 (v10) to 94 and back to 48 (58 under a strict leave-one-out reading, three of the added calls Benign).
None of the 4,038 Pathogenic or 67,733 Benign calls is among the 48, and rules with no combining rule or safeguard reached 84.6% agreement with 1,728 contradictions.
Most of the disagreement was not contradiction but a variant that ClinVar classifies as pathogenic-side or benign-side and fastVEP leaves as a variant of uncertain significance (VUS): 140,438 of the 149,280 records that are not in ClinVar's direction, among them 35,182 pathogenic-side and 105,256 benign-side variants.
These VUS calls have identifiable causes.
Thirteen criteria are not applied to a sites-only input (PS2, PM6, PM3 and BP2 need the genotypes of a proband; PS3 and BS3 need functional assays; PS4, PP1, BS4, PP4 and BP5 need case-level data; PP5 and BP6 are off by default).
Of the pathogenic-side variants called VUS, 10,693 sit one point short of Likely pathogenic and 5,184 are truncating variants that carry no PVS1, and 21,614 of the benign-side variants called VUS carry a lone BP4 or BS1, which the 2015 combining rules do not call benign.
Review of the discordant calls first identified 30 of 64 recorded changes (by our coding), and the largest class (21) was a condition of a guideline that a criterion had not applied.
The reviewer's remarks on the residual 48 opposite-direction calls pointed to what a single-variant annotator does not read: the disease model (a mechanism other than loss of function, and the mode of inheritance), the frequency that a disease allows (a variant's population frequency judged against the prevalence of the disease), the choice of transcript, and the haplotype that a variant lies on.
The ClinGen file that the benchmark read lacked 104 curations, among them CFTR, all of whose 558 null variants, each null on the gene's MANE Select transcript, were VUS.
BP7's conservation test, which the splicing recommendation cautions against without an empirical justification, costs 16,116 benign-side calls that agree with ClinVar to avoid two that contradict it.
On the 107,528 new variants agreement in direction was 70.3% and contradictions were 2.3 times as frequent (95% interval 1.2 to 4.1).
Against v37, v38 reads the mode of inheritance that ClinGen curates, which makes 2,056 genes dominant-only or recessive-only where the words of the disease names made 117, and adds 1,494 calls in direction (0.22 points) with no added reference-discordant benign call and one fewer reference-discordant pathogenic call, but four ClinVar-Pathogenic loss-of-function variants move from Likely pathogenic to VUS.
Counting one gnomAD observation as two benign criteria would have added 1,292 more calls in direction and called two ClinVar-Pathogenic variants Benign, so v38 does not.
Reading gnomAD's genomes as well as its exomes, through its joint release, would add 9,321 calls in direction (1.38 points) on every chromosome for one more reference-discordant benign call and two more reference-discordant pathogenic calls, and the release lacks the flags that withhold the frequency criteria in repeats.

**Conclusions.**
What a reader can reuse is the set of implementation choices and the way they were weighed: the default scheme, the reading of each recommendation as the conditions that its text states, and the default of each threshold that no guideline states, each with the agreement and the contradictions that it measured, and the checks that recompute the numbers that the text quotes from the stored runs.
The largest class of change was a condition that the standards state and a first implementation omitted.
Agreement in direction does not separate cautious from aggressive calling, so contradictions belong beside it.
A flat sweep can mean that the setting reaches almost nothing, or that the rule beside it makes the setting irrelevant.
The agreement is bounded: the thresholds were tuned on the April set, so the figures on the new variants are the fair ones; ClinVar can be wrong or out of date, and the expert opinion is one reviewer's; and the remaining discordance is consistent with evidence that a single-sample VCF does not carry, with design choices such as approximating PVS1's mechanism condition, and with an input that was incomplete for a reason we did not establish.

**Keywords:** variant classification, ACMG/AMP, ClinGen SVI, ClinVar, variant interpretation, benchmarking, fastVEP

---

## 1. Introduction

The interpretation of germline sequence variants rests largely on the 2015 standards of the American College of Medical Genetics and Genomics and the Association for Molecular Pathology (ACMG/AMP), which define 28 criteria of pathogenic and benign evidence and the rules that combine them into five classes: Pathogenic, Likely pathogenic, Uncertain significance (VUS), Likely benign and Benign (Richards et al., 2015).
Applied by different laboratories, the same standards can give different interpretations of the same variant, and the Clinical Genome Resource (ClinGen) has worked to minimize that variation, through general recommendations from its Sequence Variant Interpretation (SVI) Working Group and through gene-specific specifications written by Variant Curation Expert Panels (Rehm et al., 2015; Harrison et al., 2017).

The general recommendations are many, and they change criteria that looked settled.
The loss-of-function criterion PVS1 became a decision tree with graded strengths (Abou Tayoun et al., 2018).
The stand-alone benign frequency criterion BA1 acquired a list of exceptions (Ghosh et al., 2018).
The Bayesian reading of the framework turned the combining table into a point system (Tavtigian et al., 2018; Tavtigian et al., 2020).
Absence from population databases (PM2) was reduced to Supporting strength (ClinGen SVI, 2020), and computational evidence (PP3 and BP4) was calibrated tool by tool (Pejaver et al., 2022).
Splicing evidence received its own recommendations (Walker et al., 2023), and the SVI proposed that laboratories stop using the reputable-source criteria PP5 and BP6 (Biesecker and Harrison, 2018).

Software that applies the framework at scale must commit to each of these decisions, and to many that no recommendation settles.
It has to say what frequency counts as "extremely rare", how many carriers make a variant "observed in a healthy adult", and how far into an intron a variant must lie before splicing effects can be discounted.
Several tools automate part of the framework (Li and Wang, 2017; Nykamp et al., 2017; Kopanos et al., 2019).
A common way to evaluate one is to compare its output with ClinVar (Landrum et al., 2018), whose classifications are a widely shared reference but not ground truth (Section 4.2).
A single concordance figure cannot say whether a disagreement is a defect in the software, an out-of-date or mistaken ClinVar record, or evidence that no annotator can see.
And the failure that matters most in a clinical setting is not a crash but a plausible wrong answer, such as a pathogenic variant called benign, or a criterion that silently never fires.
Neither need show up in an agreement figure.

In two words, the criteria are underspecified and fragmented.
They are underspecified because no guideline states many of the quantities that software must set (the examples above: what frequency is extremely rare, how many carriers make a variant observed in a healthy adult, how far into an intron a variant must lie), and because the evidence can be combined by the table of 2015 or by points, which do not give the same call for every combination of criteria.
They are fragmented because the rules are spread over the standards of 2015, the general recommendations listed above, each of which changes criteria that looked settled, and the gene-specific specifications of Variant Curation Expert Panels, so that an implementation must commit to a reading of each of them.

Our contribution is how to implement such criteria in software so that the calls agree closely with community opinion, and what each choice costs in contradictions.
We take community opinion to be two things: ClinVar's classifications of review status two stars or better, which are a reference and not ground truth, and the review of one medical geneticist, who read the calls on which the software and ClinVar disagreed in direction over six rounds.
The paper gives the implementation choices that follow (the reading of each recommendation, the default of each threshold that no guideline states, and the combining scheme, with the reasons), the agreement and the contradictions that each measured, and the checks that recompute the numbers that the text quotes from the stored runs.
Agreement is judged by four measures: exact agreement with ClinVar's class, agreement in direction, opposite-direction calls (a pathogenic-side call where ClinVar is benign-side, or the reverse) and the recall of each side, because agreement in direction alone does not separate cautious from aggressive calling.

The claim has four bounds.
First, "closely" is relative to the settings that we measured, and a default is not always the setting that agrees most: where a sweep did not decide, thresholds were set inside the range of published expert-panel specifications, on the safer side of a trade or at the knee of the exchange rate between correct calls gained and errors added (Section 2.5), and three were not swept at all until this paper.
Second, the thresholds were tuned on the April benchmark, so only the variants that ClinVar added in October are out of sample.
Third, ClinVar can be wrong or out of date, and the expert is one reviewer, who saw both calls.
Fourth, an input of variants alone, with no genotype, pedigree or phenotype of any individual, allows 15 of the 28 criteria, which caps what any annotator of such an input can reach.

Here we describe how the ACMG/AMP classification engine of fastVEP, an open-source Rust reimplementation of the Ensembl Variant Effect Predictor (McLaren et al., 2016), was built, reviewed and measured.
The engine runs inside the annotation pass, so a benchmark of 673,660 ClinVar variants of review status two stars or better takes about two minutes to annotate in the runs of v37 (a median of 127 seconds over its 89 runs on the April set, made on 4 November 2026 up to three at a time on one Apple M4 machine with 10 cores and 24 GB that another job was also using) and about ten seconds to score.
Release 0.5.0, which holds v38, two speed-ups of its scoring and the work on the speed of annotation that followed them, and which reproduced 42 of the runs of v38 in every scored file (RUN_VERSIONS.md), was timed alone on the same machine on 7 November 2026 and took 25.3 seconds on every core (the median of three runs), 74 seconds on one thread (one run) and 11 seconds to score; the runs of v38 made for this paper were made three at a time and were not timed on a quiet machine.
That speed allowed 37 numbered benchmark runs in 24 weeks (21 April to 6 October 2026), a 38th (v38) in November, and a medical geneticist reviewed the calls on which the software and ClinVar disagreed in direction over six numbered rounds, in three of which they returned an annotated workbook.

Two versions of the classifier appear throughout.
v37 is the version in which the review rounds ended: it took a gene's mode of inheritance from the words of its disease names.
v38 is the latest version (release 0.5.0): it reads the mode of inheritance that ClinGen curates, does not count BS2's carriers beside BS1 or BA1, and classifies a gene once for each of its modes of inheritance.
Both are versions of the classifier, every run of either is on GRCh38, and every number is v38's unless the text names v37.

We report six things, and each answers part of how such criteria should be implemented.
First, how each guideline was incorporated and how agreement with ClinVar moved from one numbered run to the next (a run bundles several changes, so a change in agreement cannot be assigned to one guideline), which gives the readings that the classifier implements.
Second, what expert review found that the benchmark could not, and how the changes fall into a small number of classes, which shows where a reading goes wrong.
Third, how each numeric threshold that no guideline states was chosen and how much it matters in v38, re-measured on that code and not recalled from the version in which it was set, which shows which defaults to take care over.
Fourth, how much of the agreement depends on ClinVar-derived inputs, which bounds how much of it is circular.
Fifth, the discordance of v38 with ClinVar's 2-star-and-better classifications, on the benchmark snapshot, on a ClinVar release six months later, and on the variants that reached two stars after all tuning was done, which gives how closely the result agrees, in sample and out of it, and at what cost in contradictions.
Sixth, what v38 changes from v37, the version that the review rounds ended with, measured by the same runs of both (Section 3.8), which shows what one further reading of a source is worth.
Along the way we read each recommendation that the classifier cites for a number or a grade against its source, and report where the implementation departs from it.

---

## 2. Methods

### 2.1 The classifier and its evidence sources

The classification engine is the `fastvep-classification` crate of fastVEP (version 0.5.0-rc.1): 16,746 lines of Rust with 398 unit tests at the tag that made the runs of this paper.
That tag (`v0.5.0-rc.1`, c12d9f9) is the first release candidate of release 0.5.0 and the only tag of it so far, so that "release 0.5.0" in this paper means that candidate, and the classifier in it is v38 (5293cc0), the latest version.
v37 (5d8004a) is the version in which the review rounds ended: it took the mode of inheritance from the words of the disease names, and its crate was 14,857 lines of Rust with 340 unit tests (version 0.4.0).
v38 reads the mode of inheritance that ClinGen curates, does not count BS2's carriers beside BS1 or BA1, and classifies a gene once for each of its modes of inheritance (Sections 2.10 and 3.8).
v37 and v38 are versions of the classifier and not builds of the genome: every run of this paper, of either version, is on GRCh38 (Sections 2.2 and 2.8).
Every number of this paper is v38's unless it names v37, which it does for the history of the numbered versions v1 to v37 and as the before of a difference between the two (Section 3.8).
It is called from the annotation loop after consequence prediction and HGVS naming, once per allele and transcript, and with `--pick`, which was used throughout, one transcript is reported per variant, chosen by VEP's default order (MANE Select, MANE Plus Clinical, canonical, APPRIS, TSL, biotype, CCDS, consequence severity; MANE is described by Morales et al., 2022).
The Ensembl GFF3 that we used holds no APPRIS attribute, so that step of the order does nothing here.
Both per-variant loops of the code base, the shared annotation library and the command-line batch path, call the same classifier and the same helper for splice context, and each has end-to-end tests.

The engine implements all 28 ACMG/AMP criteria, and a run reports 26 of them, because PP5 and BP6 are off by default, following the SVI's proposal that laboratories stop using them (Biesecker and Harrison, 2018).
The benchmark input is a sites-only VCF: for each variant it gives the position and the alleles, and it carries no genotype, pedigree or phenotype of any individual.
From it 15 of the 28 criteria can be applied: PVS1, PS1, PM1, PM2, PM4, PM5, PP2 and PP3 on the pathogenic side, and BA1, BS1, BS2, BP1, BP3, BP4 and BP7 on the benign side.
The other 13 are not applied in this paper, and they fall into four groups.
Four need the genotypes of a proband and of its parents, or the proband's other variants in the gene: PS2 and PM6 (de novo status), and PM3 and BP2 (in trans or in cis with another variant).
The engine reads them from a trio VCF (`--proband`, `--mother`, `--father`), which the benchmark does not have.
Two need functional assays, PS3 and BS3, which the engine reads from a curated file when one is supplied (`--functional-evidence`), and none was.
Five need case-level data that the engine does not read: PS4 (prevalence in affected individuals against controls), PP1 and BS4 (segregation in a family), PP4 (a phenotype specific to the gene's disease) and BP5 (an alternative cause of the disease in the patient).
PP5 and BP6 are off by default.
To check this against the engine and not against its documentation, we annotated every 40th record of the benchmark (16,835 records) to JSON, which reports for each criterion whether it could be evaluated and whether it was met (Supplementary Table S1).
PS2, PS3, PS4, PM6, PP1, PP4, BS3, BS4, BP2 and BP5 were not evaluated for any of them.
PM3 was evaluated for 9,918, in genes that are not read as recessive, where the criterion does not apply, and was met for none.
Each of the other 15 was evaluated for some records and met for some.
Supplementary Table S1 also gives how often each of the 15 fires in v38, at which strengths, and on which ClinVar classes.
Two gene-level gates are referred to in the text: BP1 is withheld (the veto) in a gene with at least three pathogenic missense variants in the ClinVar protein index, where missense is an established mechanism, and PVS1 is withheld by a validity gate, which reads ClinGen's gene-disease validity curations, and by a mechanism gate, the validity gate also withholding PP2 and PM1 (Table 1; the sources are in Table 2).
Strength follows the guidelines where they grade a criterion: PVS1 by a decision tree after that of Abou Tayoun et al. (2018), whose departures from the published tree are listed in Section 3.2, PP3 and BP4 by the REVEL bands of Pejaver et al. (2022), and PM2 at Supporting (ClinGen SVI, 2020).
Criteria that read the same signal are reconciled in one pass before the call is combined, and a criterion that the pass sets aside stays in the output with the reason.
PM1 is set aside beside PVS1 and beside PP3 at Strong strength; PP3 beside PVS1 at a canonical splice site, when PP3 was met from SpliceAI, and beside PS1 or PM5, when it was met from REVEL; BP4 from SpliceAI beside PP3 from REVEL; PP3, BP4 and BP7 beside curated functional evidence, which the benchmark does not use; and, from v38, BS2's carrier count beside BS1 or BA1 (Section 2.10).
Criteria that are met are combined into a class by the default scheme, which is points on the pathogenic side and the 2015 combining table on the benign side (Table 1).
On the pathogenic side the point scale is that of Tavtigian et al. (2020): a criterion of Supporting, Moderate, Strong or Very Strong strength adds 1, 2, 4 or 8 points, and a variant with 6 to 9 points is Likely pathogenic and one with 10 or more is Pathogenic.
On the benign side the table is that of Richards et al. (2015): a Stand-alone criterion alone, or two Strong criteria, make a variant Benign, and one Strong and one Supporting criterion, or two Supporting criteria, make it Likely benign.
A call on one side stands unless the other side also reaches a call, or holds evidence of Strong strength or above without reaching one, and the variant is then a VUS; the site-level frequency veto of Section 2.7 also holds some Likely pathogenic calls at VUS.
Two alternatives were measured: the 2015 table on both sides, which is a setting of the classifier (`use_point_system = false`) and was run for this paper, and points on both sides, which has no setting and exists only as a replay of the recorded criteria (Section 2.7).
The default scheme is split because the table has no row for a lone Very Strong criterion, which the points call Likely pathogenic, and because points on the benign side would call a lone Supporting criterion Likely benign (from -1 points), which the table does not (Section 3.3).
Where the text says points or the point scale and names no side, it means the pathogenic half of the default scheme.
PP2, PM1 and PP3 on one missense variant are otherwise counted as three independent criteria, as the point scale counts them, although they can be correlated, and a combined cap on the three was measured and not adopted (Section 3.8).
Table 1 lists the guidance incorporated, when it first landed, and what the classifier does with it.

**Table 1. Guidance incorporated into the classifier.**
"Written" is the date of the commit, or, for items built in the August sprint that was merged as one change, of the benchmark run that first contained them; "live" is the run from which the criterion influenced a call, which for several items is much later than the date it was written (Section 3.2).
A numbered run bundles every change made since the one before it, so the change in agreement between two runs cannot be assigned to one row.

| Source | Criteria | What the classifier does | Written (2026) | Live |
|---|---|---|---|---|
| ACMG/AMP standards (Richards et al., 2015) | all 28 | Implements every criterion (a run reports 26, and 15 can be applied to a sites-only input, Section 2.1); five classes; the 2015 combining table is kept on the benign side | 21 Apr | v1 |
| Bayesian framework and point system (Tavtigian et al., 2018, 2020) | combining rules | Pathogenic side scored 1, 2, 4 and 8 points for Supporting to Very Strong, 10 or more Pathogenic, 6 to 9 Likely pathogenic; a Benign Very Strong counts as two Strong | 12 Aug | v11 |
| PVS1 recommendations (Abou Tayoun et al., 2018) | PVS1 | A decision tree after the published one, over the prediction of nonsense-mediated decay (NMD), last exon, fraction of protein removed, critical region and start loss; it departs from the tree in two grades and in how it predicts NMD (Section 3.2) | 29 Apr | v8 (graded), v13 (inputs derived) |
| BA1 exception list (Ghosh et al., 2018) | BA1 | Nine alleles exempt from BA1, matched by genomic coordinate first and HGVS second; extended by three hypomorphic and three reviewer-named alleles, which also block BS1 and BS2 | 29 Apr | v18 |
| PM2 recommendation, version 1.0 (ClinGen SVI, 2020) | PM2, combining | PM2 at Supporting; one Very Strong plus one Supporting reaches Likely pathogenic. The recommendation gives no frequency threshold, and the thresholds are not from it (Section 2.5) | 29 Apr | v6 (absence read as rarity) |
| PM3 recommendation, version 1.0 (ClinGen SVI, 2019) | PM3 | Points-based in-trans scoring; needs the genotype of a proband, which the benchmark lacks | 29 Apr | n/a |
| Calibration of in-silico tools (Pejaver et al., 2022) | PP3, BP4 | REVEL at 0.644, 0.773 and 0.932 for PP3 Supporting, Moderate and Strong; at 0.290, 0.183, 0.016 and 0.003 for BP4 Supporting to Very Strong; missense only; PP3 is dropped beside PS1 or PM5, and PM1 beside a Strong PP3 or beside PVS1 | 29 Apr | v1 |
| Splicing recommendations (Walker et al., 2023) | PP3, BP4, BP7, PS1, PVS1 | SpliceAI at least 0.2 for PP3 and at most 0.1 for BP4, both at Supporting strength; BP7 not for a synonymous change at an exon's first base or last three, and for an intronic change only at +7 or beyond on a donor side and -21 or beyond on an acceptor side, but with the 2015 conservation test kept and SpliceAI below 0.2, where the recommendation applies BP7 after BP4 (SpliceAI at most 0.1) and cautions against the conservation test (Sections 3.2 and 3.3); PS1 for a variant on the same canonical splice dinucleotide as a Pathogenic one; splice BP4 set aside beside a REVEL PP3 | 29 Apr | v4 (PP3, BP4, BP7), v16 (PS1), v37 (BP7 exon edge) |
| Filtering allele frequency (Whiffin et al., 2017; ClinGen VCEP Review Committee, 2024) | BA1, BS1 | Tested against gnomAD's grpmax filtering allele frequency (the highest, across ancestry groups, of the lower bound of the 95% interval of each group's frequency) rather than a point estimate | 12 Aug | v10 |
| Homologous genes (Mandelker et al., 2016) | BA1, BS1, BS2, PM2 | Not evaluated in 29 genes whose population frequencies are unreliable | 12 Aug | v9 |
| Functional evidence (Brnich et al., 2019) | PS3, BS3 | Curated file with a strength column; supersedes PP3, BP4 and BP7 for that variant | 12 Aug | available from v10; not used in the benchmark |
| ClinGen gene-disease validity (Strande et al., 2017; ClinGen, 2026) and curated disease mechanisms | PVS1, PP2, PM1 | Withheld where every ClinGen curation of the gene is Limited, Disputed, Refuted or No Known Disease Relationship (the validity gate), and where the curated mechanism excludes loss of function (the mechanism gate, which applies to PVS1 only); from v38 the mode of inheritance that ClinGen curates for each disease is read as well, in place of the words of the disease name (Section 3.8) | 12 Aug | v10 (a first version blocked on absence), v11 (curated-negative only), v38 (inheritance) |
| Published VCEP thresholds (ClinGen CSpec Registry) | BA1, BS1, PM2 | 286 thresholds across 117 genes, each with its source sentence kept; of 411 criterion entries read, 18 are withheld where two panels disagree and 107 were rejected, among them 37 that require absence, 27 conditioned on inheritance and 21 with no threshold (Section 3.3); opt-in, not in the benchmark | 13 Aug | opt-in |
| Genes with multiple disorders (ClinGen, 2024) | all | A gene is classified once for each mode of inheritance of its established diseases, besides the call for the gene as a whole, which is unchanged (Section 3.8); the per-disorder override scaffold of 29 April is still not read | 29 Apr (scaffold); after v37 (the calls) | n/a, the call for the gene as a whole is unchanged |

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
| ClinGen gene-disease validity (Strande et al., 2017) | file created 30 April 2026: 2,924 genes, 3,546 curations, 2,739 of them Definitive, Strong or Moderate and the rest Limited, Disputed, Refuted or No Known Disease Relationship; it lacked 104 curations that ClinGen dates before that day (Section 3.5) | PVS1, PP2, PM1 gate; inheritance for PM2, BS2, PM3 and BP2 (v37 read it from the words of the disease names, v38 from the curated field of each curation) |
| RepeatMasker (Smit et al., n.d.) | UCSC track | BP3 |
| Ensembl (McLaren et al., 2016) | release 115 GFF3 and GRCh38 primary-assembly FASTA | consequence, HGVS, PVS1 tree inputs |

### 2.2 Reference set, scoring and definitions

The reference is ClinVar's VCF for GRCh38 dated 4 April 2026.
We kept records whose review status is at least two stars (criteria provided, multiple submitters, no conflicts; reviewed by expert panel; practice guideline) and that are single-nucleotide variants or indels of at most 100 bases.
A record's class is read from its CLNSIG by substring match in a fixed order (Pathogenic/Likely pathogenic, then Benign/Likely benign, then Likely pathogenic, Likely benign, a leading Pathogenic or Benign, uncertain significance, and conflicting classifications), and a record that matches none, a drug response for example, is dropped.
That gives 673,660 records: 79,823 Pathogenic, 13,989 Likely pathogenic, 295,298 VUS, 128,038 Likely benign and 156,512 Benign, of which 652,400 are two-star, 21,237 are three-star and 23 are four-star.
ClinVar's combined labels are mapped to the first class (Pathogenic/Likely pathogenic to Pathogenic, Benign/Likely benign to Benign), so exact agreement is a strict measure and agreement in direction is the fair one.
The VUS class is the 295,298 variants that submitters classify as uncertain without conflict, except that two RYR1 records of expert-panel review whose classification is conflicting are labeled VUS by the label script, and fastVEP calls both Likely pathogenic.
Four records of CFH and CFI whose classifications are Likely pathogenic and Likely pathogenic with low penetrance are labeled Pathogenic by a pattern match of the same script.
We kept both labels so that the 38 runs stay comparable: relabeling the four changes exact agreement by one variant and agreement in direction by none, and leaving out the two RYR1 records lowers the number of ClinVar VUS called pathogenic-side by two.
A variant's gene is the first symbol of the GENEINFO field of its ClinVar record, and the population maximum is the highest allele frequency among the ten population groups that the annotation carries, "other" and "remaining" included, which is what the frequency criteria fall back on where there is no filtering allele frequency.
We call these classifications a reference and not ground truth throughout.
Together with the expert review of Section 2.4 they are what this paper means by community opinion: the classifications on which ClinVar's submitters or an expert panel agree, and one reviewer's reading, and not a poll of the field.

Pathogenic-side means Pathogenic or Likely pathogenic and benign-side means Likely benign or Benign (P/LP and B/LB in the tables), and the recall of a ClinVar group is the share of its variants that the classifier puts in the same direction.
A committed call is a call of one of those four classes, and a quiet variant is a synonymous, intronic (polypyrimidine-tract variants included), untranslated, flanking (upstream or downstream of a gene) or non-coding one.
"Reference-discordant benign" and "reference-discordant pathogenic" name the direction in which a call contradicts ClinVar, and "correct" means in agreement with it; none of them says that ClinVar is right or that a call is a diagnostic error, and the reviewer supported ClinVar's class for 13 of the 24 reference-discordant benign calls and 16 of the 24 reference-discordant pathogenic calls of v38 (Section 3.6).

The scorer takes the call on the picked transcript.
A call is exact if it equals ClinVar's class; same-direction if both are pathogenic-side (Pathogenic or Likely pathogenic), both benign-side, or both VUS; and opposite if ClinVar is pathogenic-side and the call benign-side, or the reverse.
An opposite call is a reference-discordant benign call when ClinVar is pathogenic-side and a reference-discordant pathogenic call when ClinVar is benign-side.
Rates use all 673,660 records as the denominator, which is the smaller figure: the 295,298 ClinVar VUS cannot be contradicted, so the 48 opposite-direction calls of v38 are 0.0071% of the records, 0.013% of the 378,362 with a pathogenic-side or benign-side class and 0.020% of the 245,954 committed calls.
Among the 673,660 records, 382 have `.` as their alternative allele.
Each of them names the reference sequence itself in its HGVS expression (for example `NC_000001.11:g.12358644=`), and ClinVar classifies every one of them as benign-side (377 Benign and 5 Likely benign in the labels of the benchmark).
There is no alternative allele to annotate, so they receive no call, and the scorer counts them as not concordant, which lowers the recall of benign-side variants by at most 382 of 284,550.

### 2.3 Iterative development

The unit of the process is a run: the full benchmark annotated with `fastvep annotate --acmg --pick` and scored by one script.
After every change to criteria, thresholds, data or annotation, the benchmark was re-run and compared with the run before it.
We refer to each numbered run as v1 to v38.
From v37 on, the harness records what produced an annotation (flags, the hash of the binary and of every configuration, which inputs were read, and the commit) and refuses to reuse an annotation that something else produced.
From v9 on, runs are scored in a leave-one-out mode in which a variant's own ClinVar record is excluded from the evidence that PS1 and PM1 read; a run that lets the record count is called ClinVar-informed.
PM5 cannot match its own record, because it asks for a different amino acid change at the same residue, and the critical region of PVS1 reads residues downstream of the variant.
The veto of BP1 and the allele-level ClinVar gates read the record as it is (Section 3.4).
Where a change touched the annotation layer (consequence or HGVS), every variant whose call moved was re-annotated with Ensembl VEP release 115.1 on the same GFF3 and FASTA with `--pick`, after confirming that both tools reported the same transcript, to establish which tool the change brought the software closer to; the written record keeps the outcome only as ledger items (3 were found by this comparison), and we do not report the comparison itself.

Not every numbered run is a step on the way.
Run v3 was a partial run that was overwritten before it was preserved, v21 and v22 repeat v20, v24 and v25 were intermediate runs that are not versions that the project kept, and v17 is the opt-in run with the published VCEP thresholds, which v18 does not carry.
We mark v17 in the figures and leave it out of every range that the text quotes.
The raw output of v1 to v6 was overwritten, so v2 and v4 to v6 are quoted from the figures recorded at the time, and v1 from a matrix that was reconstructed from the documentation of the run: scored matrices are stored for 32 of the 38 runs, v7 to v38.
Every rate is over the 673,660 records of the benchmark; v1 annotated only 627,376 of them, and its recorded 54.7% is over those, where over all records it is 51.0%.
We could not check the denominators of v2 to v6.

We also kept a ledger of the changes that the project's written record (RUN_VERSIONS.md and the documents of the review rounds) describes as changing a criterion, a threshold, an input or the scoring: 64 items, each coded for its class and for the route by which it was first found.
The coding was done by one coder from the written record, was not repeated independently, and is not a record of every commit.

### 2.4 Expert review

A clinical molecular geneticist (M.G., an author of this paper) reviewed the calls on which the software and ClinVar disagreed in direction.
Each row of the first workbook carried 75 columns of evidence, and later workbooks carried more: ClinVar's record, the criteria that fired, REVEL, SpliceAI, PhyloP, gnomAD counts and frequencies by ancestry, gene constraint and the gene's ClinGen disease annotations.
The reviewer saw only the opposite-direction calls, with ClinVar's class and ours both visible, and not the pathogenic-side variants that we left as VUS, nor any concordant call.
The reviewer was not blind to the software or to its calls, and supplied two inputs to the classifier (the table of Ghosh et al. and the names of three alleles; ledger items 50 and 63), so they were a contributor as well as an evaluator, and their rulings are one expert's reading of rows on which both calls were visible and not independent evidence of which call is right.
A second clinical geneticist (W.C., an author) advised on a subset of the cases; the rulings that this paper counts are those of the systematic review by M.G.
The reviewer returned annotated workbooks in rounds 1, 2 and 6, with the verdict recorded as a cell color in rounds 1 and 2 and as text in round 6, and a covering email in round 2; rounds 3 to 5 were tables and notes we sent after re-runs, which were not returned marked.
The workbook of round 1 carries a list headed "Discordance categories" that names six kinds of discordance, and a key of two colors (a correct call, and both calls discordant with the reviewer's call); the key of round 2 adds a query, a request for more evidence, and a marker for genes whose population frequencies cannot be relied on because of homology.
We decoded the colors against each workbook's key, and read any green on a call as that call being correct (in round 2 the green on the cells is a different shade from the key's, and in round 1 three greens on fastVEP's call are darker than the key's).
In round 6 the workbook asked for the reviewer's class and for who is right in two columns, which came back empty, and the reviewer wrote their ruling as free text in the verdict column of all 64 rows and added a note to 18.
We mapped the rulings of all three rounds onto five categories (ClinVar's call supported, fastVEP's call supported, neither, needs more evidence or depends on the disorder, and no verdict).
The round-6 text was coded by us from its content by five rules given in the supplement, and the reviewer did not assign these categories.
Between rounds the reviewer's notes were carried onto the rows they belong to, and each round was answered with a note that said what had changed, what it cost, and what was measured and not done.

### 2.5 How thresholds were chosen

Thresholds fall into three kinds.
Some are fixed by a guideline and were not tuned: the REVEL bands, the SpliceAI cut-offs, BA1 at 5%, and the points that map to each class.
Some are not stated by any guideline: the PM2 thresholds, the BS1 threshold, the BS2 thresholds, BP7's intronic distance, BP7's PhyloP threshold and the count of pathogenic missense variants at which BP1's veto applies.
Most of these were set by sweeping one parameter over the benchmark with everything else at its default, and three (the recessive-only PM2 threshold, BP7's PhyloP threshold and the count of BP1's veto) were set without a sweep and were swept for the first time for this paper.
A few are structural, such as which combining scheme is used on which side, and were compared directly.
The rules applied when a sweep did not decide, written down as they were used, were these.
Prefer a value inside the range that published VCEP specifications use; we measured that range from the extracted CSpec table (Section 3.3).
Prefer the safer side of a trade, because, if ClinVar is right, a benign call on a pathogenic variant costs more than a variant left as VUS.
Take the knee of a curve where there is one, which we took to be the point past which the exchange rate (the extra correct calls per extra opposite-direction call) falls from hundreds of calls per added error to about a hundred or fewer.
Where a reviewed variant showed a rule to be wrong, fix the rule; where it showed a variant to be an exception to a correct rule, put it on the curated exception list; and never move a threshold to repair one variant.
Both the exception list and the dominant-gene BS2 rule were written after the reviewer's rows, so the 48 residual calls of Section 3.6 are in-sample: without the exception list there are 54, and there are 58 when the exception list, BP1's own-record count and ClinVar's low-penetrance gate are removed together (Section 3.4).
Four more groups of thresholds that are ours (seven numbers) were not varied, because we made no sweep of them: the pLI and LOEUF thresholds of PVS1 (0.9 and 0.35), the missense z-score thresholds of PP2 and BP1 (3.09 and 2), PM1's window and count of pathogenic neighbours, and the floor of homozygotes for BS2.

### 2.6 Re-measuring on v38

The sweeps behind the default values were run on earlier versions of the code, in some cases on a sample, and the combining rule changed after several of them.
For this paper we re-ran every sweep of a threshold that no guideline states on v38, as full-benchmark runs that differ from v38 in one setting: 8 values of the PM2 threshold for genes not read as recessive-only, 6 of the PM2 threshold for genes read as recessive-only, 6 of the BS1 threshold, 4 of the BS2 homozygote threshold, 7 of BP7's intronic distance, 5 of BP7's PhyloP threshold and 5 of the count of BP1's veto (41 runs).
Four more values of the PM2 threshold were run under the 2015 combining table, to see whether the threshold's importance belongs to the threshold or to the rule (the interaction of Section 3.3), and we used the 8 runs of the dominant-gene BS2 settings (5 of the prevalence threshold and 3 of the allele-count floor), which Section 3.8 had already made on v38.
To measure what each design decision buys, we also ran 16 single-switch ablations and 9 leave-one-source-out ablations, in which one source is removed from the annotation directory so that the criteria that read it lose their input (for gnomAD the run also tells PM2 not to read an absent record as rarity, because PM2 would otherwise be met for every variant).
The two point cutoffs of the pathogenic side, Likely pathogenic from 6 points and Pathogenic from 10, are constants in the code and not settings, so for 9 further runs we built the binary with a patch that changes eight lines (six added, two removed) and reads them from the environment (`manuscript/acmg/replay_validation/combiner_cutoff_env.patch`).
Two departures from the sources that we found while reading them again were measured with a second measurement build, made with a patch that exposes through the environment BP1's veto counting a variant's own ClinVar record and BP7's SpliceAI threshold (`bp_hooks_env.patch`): a control, BP1 without the own record, BP7 with SpliceAI at most 0.1, and BP7 with that and no conservation test (4 runs); a fifth run with the same build removes together the three adjustments that read the variant's own ClinVar record or followed the reviewer's rows (BP1's veto counting that record, ClinVar's low-penetrance gate and the curated frequency-exception list), each of which had been measured alone; two more runs with that build give BP1's veto a count of 1 and a count of 2 and leave the variant's own record out, the pair that Section 4.5 proposes together.
Three runs change an input to see what the gap in PVS1 of Section 3.5 is made of: the ClinGen file with the curations that its export lacked put back, PVS1 licensed for genes that Orphanet labels as a loss-of-function disease, and both (3 runs).
With nothing set, each measurement build is the default classifier: the control of each reproduces v38 in every cell of the concordance matrix and in the firing count of every criterion (Appendix Table A1).
All 89 runs (Supplementary Table S2) were made by the same harness: 73 with the plain build of v38 (SHA-256 `70b7c684...`), 9 with the build of the first patch and 7 with the build of the second.
The three binaries are builds of the tag `v0.5.0-rc.1` (c12d9f9), the release that holds v38, the second and the third with one patch each, and each was made from a clean export of the tag in a target directory of its own (`build_manuscript_binaries.sh`).
Eighteen of the plain runs are sweeps of the PM2 thresholds and of the BS2 homozygote threshold that Section 3.8 had made first, with the binary of the working tree that became v38 (SHA-256 `f184277d...`); each was made again on 7 November 2026 with the build above, with every scored file identical, and the other 71 runs were made once, afterwards, 55 with the build above and 16 with the builds of the two patches.
The same 96 runs were made for v37, with binaries built from commit 5d8004a (plain build SHA-256 `1b8e7c46...`), twice: first in October 2026 with binaries built at earlier commits of the same branch (the plain build at commit 5d2521f, SHA-256 `d12c0d8e...`), and again on 4 November 2026 with the builds of commit 5d8004a, the main run `v37` and the earlier runs of the BS2 threshold included.
The two sets agree in every scored file of every run, and in every record of the annotated VCFs that both keep; the only difference is that five CSV files of each run made on 6 October end their lines with CR LF, which the scoring script no longer writes (RUN_VERSIONS.md).
The provenance file of each run of either version is in the repository, and records the flags, the hash of the binary and of every configuration, the inputs read (by name and size, not content), and the commit and number of uncommitted files of the checkout that ran the harness.
The harness refuses to reuse an annotation unless the binary, the flags, the configuration and the inputs are those that made it; the environment variable that sets a run of a patched build is not part of that check and is recorded in the run matrix, not in the provenance file.

### 2.7 Offline replay of the combining rules

The combining step is a pure function of the criteria a variant met, and the output records them.
We therefore re-derived calls under alternative schemes, and under alternative band edges, from the recorded criteria without re-annotating.
The replay of the default scheme reproduces v38 for 673,641 of 673,660 variants.
The remaining 19 are variants that the replay calls Likely pathogenic and the classifier holds at VUS by its site-level frequency veto: when gnomAD rejects a site (by its own FILTER, a low-complexity tract or a segmental duplication) and the frequency that was withheld would itself have met BA1 or BS1, a Likely pathogenic call becomes VUS.
The output does not record the detail the veto needs, so the replay holds those 19 at VUS wherever it would otherwise call them Likely pathogenic.
The replay was checked against real runs of the same configurations (Appendix Table A2): the Richards 2015 table on both sides, and the point cutoffs of Section 3.3.
It reproduces the table run, the control and 6 of the 8 runs with a changed cutoff cell for cell, and the other two (Likely pathogenic from 4 and from 5 points) within 3 and 1 variants, each a variant that the replay calls Likely pathogenic and the run, through the veto, does not.
We use the real runs wherever one exists, and the replay only for the benign-side point scales, which have no switch in the default code and for which we made no real run (Section 3.3).

### 2.8 A later ClinVar release

To test v38 on variants that none of its tuning used, we applied it unchanged to ClinVar's release of 4 October 2026, extracted the same way (781,135 records at two stars or better).
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
Every number in the text and every table and figure is produced by `analysis/acmg_benchmark/scripts/11_manuscript_runs.py`, `12_per_variant_table.py`, `13_manuscript_stats.py`, `14_temporal_holdout.py` and `15_clingen_complete_export.py`, and by `manuscript/acmg/generate_figures.py`, from the scored runs; the aggregates that no other table holds are in Appendix Table A3.

### 2.10 The changes of v38 and how they were measured

v38 differs from v37 in four ways, and each change to a call was measured with the harness of Section 2.6 against the commit before the change.
That commit has the classification code of v37 and a newer annotation layer (Section 3.8).
The derivation of the pseudoautosomal flag and the report of a duplicated source change no call of the benchmark and are covered by tests only (`gnomad.rs` and `tests/sa_duplicate_source.rs`).
Four smaller changes touch no call of the benchmark and are covered by tests too.
BP2's trans arm needs a gene whose established diseases are dominant and none recessive, as BS2's carrier arm does, because a gene with a recessive disease explains the same observation as the evidence that PM3 counts (BP2 needs a proband, which the benchmark lacks).
The annotation library behind both web servers looks an indel up in the databases in the form the VCF wrote it, as the command line does, where it had asked in the form the parser trims it to and missed it: on 120 ClinVar indels of chromosomes 21 and 22 it found gnomAD for none and the command line for 67, and now the two agree on all 120 (RUN_VERSIONS.md).
The converter of the ClinGen export finds its columns by their names and stops if one is missing.
And a gene database is written in the order of its symbols, so that building it twice gives the same file.

First, the mode of inheritance is read from the field that ClinGen curates.
A ClinGen curation reaches the classifier as one phenotype string that ends in a tag, `<disease label> (ClinGen <classification>/<mode of inheritance>, <MONDO identifier>)`, which the converter of the ClinGen export writes from the curation's own columns.
Up to v37 the criteria that branch on a gene's inheritance (PM2, BS2, PM3 and BP2) searched the phenotype string for the words "dominant" and "recessive", and the tag was read for its classification alone.
From v38 one module reads the tag: AD is dominant, AR is recessive, SD (semidominant) is both, and XL, MT and UD are neither, which is the arm that each criterion takes where inheritance is unknown.
A mode that the module does not know is read as no information and is not guessed at.
A phenotype with no tag, which is all that the legacy OMIM source has, is read by its words as before, and where a gene has tagged phenotypes the untagged pieces that the database builder makes by splitting a label at a semicolon are not diseases of their own.
A test builds the gene database with the command-line builder from 263 real curations in 185 genes of ClinGen's export of 11 October 2026, and checks what each gene reads as against the curations' own columns.
The genes cover every mode of inheritance that the export uses, the three labels that contain a semicolon and the three X-linked disorders whose label says "dominant".

Second, BS2's carrier count is set aside beside BS1 or BA1.
In a dominant gene BS2 counts carriers, which is the allele count less the homozygotes and so the count that BS1 and BA1 turn into a frequency, and the threshold of its carrier arm is lower than theirs.
BS1 and BA1 test the filtering allele frequency of the population with the highest one and the carrier arm tests the global carrier count, so a lower threshold does not make the arm fire wherever they do (it fired on 7,234 of the 7,328 variants that meet BS1, Section 3.8), but where both fire they summarize one gnomAD record twice.
We read that as one observation counted twice, which is our reading and not a statement of the ClinGen SVI, and the 2015 table that the benign side keeps reads two Strong criteria as Benign.
The reconciliation pass of Section 2.1 now sets the carrier count aside beside either, and the arm acts alone only below BS1's threshold.
The arm for homozygotes and hemizygotes is not set aside: a homozygote is a genotype that an allele frequency does not give, and the hemizygotes of the non-pseudoautosomal X are part of the allele count, so for them the doubling remains (Section 3.8 counts it).

Third, a gene with more than one mode of inheritance is classified once for each mode.
Where the Definitive, Strong or Moderate curations of a gene carry more than one mode, the classification is repeated for each mode with the gene's disease annotation restricted to the curations of that mode.
The calls are reported beside the call of the gene as a whole, in the JSON output as `per_disorder`, in the CSQ field `ACMG_DISORDERS` and in the web interface.
The call of the gene as a whole is the call it would have without them, and it is the call that every table of this paper scores, because a ClinVar record of the benchmark has no disease-specific class to score a call against.

Fourth, the allele counts that the frequency criteria read can come from gnomAD's joint release of exomes and genomes, and the loader guards what that exposed.
In the headers read on 4 November 2026 the joint release declares none of the three flags that the exomes carry for regions where gnomAD's counts are unreliable and for the pseudoautosomal regions (`lcr`, `segdup` and `non_par`), so a database built from it cannot withhold BA1, BS1, BS2 and PM2 in a low-complexity region or a segmental duplication, and cannot tell BS2 which hemizygotes to count.
The database builder now derives the pseudoautosomal flag from the GRCh38 coordinates, which agree with the flag in gnomAD's genomes file at the edges that could be read (exactly at the end of the first region of the X chromosome, and within the gaps between the records at the others), and prints a warning that the other two are missing.
The builder also reads the names of the filters that failed in each callset from two keys of the joint release (`exomes_filters` and `genomes_filters`), because the FILTER column of the release says only which side failed, and it reports a filter that fired in either callset for the record.
The annotation loops of the library and of the command line record a source that two databases both answer, which until v38 the later path of the two settled without a message, and the command line says so at the end of the run; the case that exposed it was a second gnomAD database in the same directory.
A script extracts the benchmark's regions of the exomes, the genomes or the joint release from gnomAD's public files and builds the databases (`analysis/acmg_benchmark/scripts/sa_sources/build_gnomad_release_per_chrom.sh`).

Each alternative to a change is a run that differs from v38 in one setting, made with a build of a patch that reads the setting from the environment (`replay_validation/inheritance_alternatives_env.patch`).
The control of that build, with nothing set, reproduces v38 in all 30 cells of the concordance matrix and in the firing count of all 24 criterion codes (one for each strength at which a criterion fired, Section 3.8).
The alternatives are counting BS2's carriers beside BS1 and BA1, as the reading first ran; withholding the carrier arm for a disorder whose label says susceptibility or predisposition; reading the field for PM2, PM3 and BP2 and leaving BS2 as in v37; and capping PP2, PM1 and PP3 together at 4 points (Strong), the cap that the reconciliation pass already applies to PP3 and PM1.
The four thresholds that read inheritance (the PM2 threshold for genes read as recessive-only and the threshold for the rest, the BS2 threshold for individuals with no functional copy, and the BS2 carrier threshold with its floor on the allele count) were swept again on v38 in 26 runs.
The adopted change, the first and the last alternative and the commit before the change were applied unchanged to the release of 4 October 2026 (Section 2.8) in 4 runs.
Three runs on every chromosome compare the gnomAD exomes with the joint release and with the exomes when the region flags withhold nothing, and a pilot of the same three runs on chromosomes 21 and 22, made before the other chromosomes had been extracted, was kept as a check.
In all there are 38 runs on the April set and 4 on the October release, each with the provenance record of Section 2.6.
All but one of the 42 were made twice, as the runs of v37 were (Section 2.6): first with binaries built from earlier states of the same working tree, before its last edits, and again on 4 November 2026 with binaries built from the final tree by `build_manuscript_binaries.sh` (the commit before the change, the adopted change, and the adopted change with the patch).
The two sets agree in every scored file of every run and in every record of the annotated VCFs that both keep (RUN_VERSIONS.md); the run on the joint release on every chromosome was made once, on 5 November 2026, with the binary of the second set.
All 42 were made once more on 7 November 2026 with binaries built from the tag of release 0.5.0, and agree with the first runs in every scored file (Section 2.6 and Data and code availability).

## 3. Results

### 3.1 Thirty-eight runs: where agreement came from, and where it stopped

Same-direction agreement with ClinVar rose from 51.0% in v1 (a matrix reconstructed from the run's documentation) to 77.8% in v38 (Figure S1A; Supplementary Table S2 gives every run).
Almost all of that rise came in the first 13 runs, and the largest steps were mostly plumbing and not new guidance.
Agreement rose 8.8 points between v2 and v4, when SpliceAI and PhyloP results finally reached the criteria that read them (items 6 and 7 of the ledger of Section 3.2), and 5.2 points between v5 and v6, when a variant absent from gnomAD began to count as absent for PM2 and BP4's splice path was confined to splice questions (items 9 and 10 of the ledger), two conditions of the guidelines that the criteria had not applied.
The first review (v7 to v8) added 2.8 points.
The next run, v9, lost 2.7 points, and that loss was the correct result: PS1 and PM1 had been counting each variant's own ClinVar record, and 94% of PS1's firings were self-matches, so scoring them leave-one-out removed agreement that had never been evidence.
Of the 2.7 points, about 0.6 were the cost of the review's criteria fixes and about 2.0 were the leave-one-out scoring (RUN_VERSIONS.md).
Recall of ClinVar's pathogenic-side variants fell from 63.2% at v8 to 45.6% at v9 for the same reason, and the 62.5% of v38 is 0.7 points below the 63.2% that v8 reached by counting each variant's own record.
Agreement then rose by 3.2 points at v10, which bundled a rebuilt gnomAD database with its quality fields, the gene-level gates, and the first two measured thresholds (the PM2 threshold and the BS2 prevalence threshold), by 0.9 points each for the SVI point scale on the pathogenic side (v11) and the BS1 threshold (v12), and by 1.6 points at v13, when `--acmg` was made to enable HGVS and the criteria that read it came alive.

After v13 the headline stopped moving.
Across the next 23 runs (v14 to v37, without v17) same-direction agreement stayed between 77.50% and 77.62%, and no run moved it by more than 0.07 points.
The 38th run, v38, which reads the mode of inheritance that ClinGen curates (Section 3.8), moved it by 0.22 points to 77.84%, the one step after v13 of more than 0.07 points.
The work did not stop, and its effect shows in the other panels.
Opposite-direction calls fell from 94 at v13 to 49 at v37 and 48 at v38, reference-discordant benign calls from 50 to 24 at both, and reference-discordant pathogenic calls from 44 to 25 at v37 and 24 at v38 (Figure S1B).
The path was not monotonic.
Opposite-direction calls had reached a minimum of 25 at v10 and then rose to 56 at v11, when a lone PVS1 became Likely pathogenic under the point scale of the pathogenic side, which gained 6,632 correct pathogenic-side calls for 31 additional reference-discordant pathogenic calls.
They rose again at v13, when BP7's intronic extension started to fire, which gained 11,339 benign-side calls for 36 additional reference-discordant benign calls.
Both were trades made on purpose and recorded as trades, and Section 3.3 returns to how they were weighed.
The matrix stored for v1 holds 31 opposite-direction calls, fewer than v38, because v1 committed to only 74,355 variants where v38 commits to 245,954; per committed call, contradictions fell from 0.042% in v1 to 0.020% in v37 and in v38, with a minimum of 0.011% at v10.

Recall of ClinVar's pathogenic-side variants in direction has stayed between 62.0 and 62.5% since v11, and recall of its benign-side variants between 62.3 and 62.7% from v13 to v37 and 62.9% at v38 (Figure S1C).
That flat line is the structural fact about the benchmark that Section 3.5 examines.
The agreement figures also have an uneven past: the raw output of v1 to v6 was not kept, and v2 and v4 to v6 are quoted from the figures recorded at the time with denominators that we could not check, so their values carry less weight than the stored runs from v7.

### 3.2 What review of the discordant calls found

**Review of the discordant calls first identified 30 of the 64 changes, by our coding, and the largest class was a guideline condition left out.**
The ledger of Section 2.3 holds 64 items, each coded by us for its class and for how it was first found (Figure S2, Supplementary Table S3).
Expert review first found 30 of the 64 (47%), code audit 8, our own reading of the discordant calls 7, the first real-data run 5, an aggregate figure of the benchmark 5, a threshold sweep 4, a comparison with Ensembl VEP 3, and a regression in a later run 2.
All of the review's input was the benchmark's own list of opposite-direction calls, so this compares the reviewer's reading of that list with other routes and not two independent ways of finding the same thing, and "first" is our reading of the written record.
Two of the seven that we found by reading the discordant calls ourselves, the self-matches of PS1 and PM1 and the two modes of reporting that followed, came from reproducing the reviewer's rows, which could not have shown them.
By class, 21 items were a condition of a guideline that the criterion had not applied: BS2's "full penetrance expected at an early age", PM1's "without benign variation", BP1's need for a truncating-only mechanism, PVS1's need for a gene-disease relationship and for loss of function as the mechanism, and BP4's confinement to variants for which a calibrated predictor exists.
Seventeen were wiring defects, in which a criterion or an input was inert, misread or fed the wrong field.
Eight concerned how far a population frequency can be trusted, five were annotation differences from Ensembl VEP, five were thresholds set by measurement, four concerned scoring and reporting, two were circularity with ClinVar, and two were combining rules.

**The reviewer supported ClinVar's call far more often than the software's.**
The medical geneticist returned annotated workbooks in three rounds (Figure S3A, Appendix Table A4).
In round 1 the workbook held 101 opposite-direction calls from v7, and the reviewer judged ClinVar's call correct on 86, fastVEP's on 3, neither on 8, and left 4 with a color that was not in the key.
Round 2 held 122 calls from v8, of which 93 were new and 29 carried over: ClinVar's call was judged correct on 65, fastVEP's on 2, neither on 8, more evidence was needed on 1, and 46 carried no recorded verdict.
Seven of the 65 (four in CYP21A2 and three in STRC) carry the marker for a gene with a homology problem, and on six of them it is on fastVEP's call as well as the green on ClinVar's; the project's note of that round counted those six as not marked, which gives 59 and 52 where we count 65 and 46.
Round 6 held 64 calls from v26: by our coding of the returned text, the reviewer supported ClinVar on 43, fastVEP on 4, neither on 9, and held 8 pending more evidence or by disorder.
So when the software and ClinVar disagreed in direction, the reviewer sided with ClinVar on 89% of the 97 rows of round 1 that carried a verdict in the key, 86% of the 76 rows of round 2 that carried one, and 67% of the 64 rows of round 6.
The three percentages come from different instruments, which are color decoding against a two-color key (round 1), color decoding against a richer key (round 2) and our coding of free text (round 6, where 14 rows were coded from the evidence they cite and not from a stated class), so they are not one series and are by our coding.
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
Regrading both by one step would change 2 of 673,660 calls (1 from VUS to Likely pathogenic and 1 from Benign to VUS), by a replay that is an upper bound, because the recorded code does not say which branch gave a PVS1_Moderate; so few change because a grade one step lower plus PM2 seldom reaches 6 points either way.
The same replay on the recorded criteria of v37 changed 7 (6 from VUS to Likely pathogenic and 1 from Benign to VUS).
Second, PVS1 predicts NMD by whether the variant lies in the last exon, where the tree predicts it by the 50-nucleotide rule, which also lets a stop in the last 50 nucleotides of the penultimate exon escape; the exact rule is a setting, the code comment calls it the guideline-faithful one, and it is off by default because of what it costs against ClinVar: 603 calls in direction and 0.73 points of pathogenic-side recall, to avoid one reference-discordant pathogenic call (Section 3.4).
That is a default chosen on agreement with ClinVar, a tuning to the reference that is not the dependence on ClinVar-derived inputs that Section 3.4 measures.
Third, BP7 keeps the 2015 test of poor conservation (PhyloP below 2.0) and takes a SpliceAI score below 0.2 as no impact, where Walker et al. (2023) apply BP7 after BP4, which for SpliceAI is a score of 0.1 or below, and caution against a conservation test; Section 3.3 measures it.
Fourth, the allele-number floor has a narrower source than its use (Section 3.3).
Fifth, the recessive-gene PM2 threshold of 7e-5 was attributed in the code to the SVI's PM2 recommendation, which gives no frequency at all; the number equals the recessive threshold of the Hearing Loss expert panel, the project's record does not say where it came from, and Section 3.3 sweeps it.
Sixth, the veto of BP1 counts the variant's own ClinVar record, which PS1 and PM1 exclude (Section 3.4).
Seventh, the ClinGen file that the benchmark read lacked 104 curations that ClinGen dates before the day of the export, among them CFTR (Section 3.5).
Eighth, three citations were wrong: the PMID of Ghosh et al. (2018), which was wrong in three places, the journal of Tavtigian et al. (2020), which is Human Mutation and not Genetics in Medicine, and the date of the ClinGen guidance on genes with multiple disorders, which is January 2024 and was recorded as July 2025.
We corrected the citations and the comments in the repository, and Sections 3.3 to 3.5 measure what each departure would change, so that the decision to change the classifier can be made on a number (Section 4.5).

### 3.3 Where the thresholds sit

Each threshold that no guideline states was set once, in a particular run, against the combining rule and the other criteria as they then stood.
Re-measuring them on v38, the latest classifier, changed what we can say about several (Figure 1, Table 3).

**Table 3. Thresholds that no guideline states: default values, and what v38 shows.**
Exchange rates are extra correct benign-side (or pathogenic-side) calls per extra reference-discordant benign call (or reference-discordant pathogenic call) for the step on either side of the default value; "free" means no added error.
Correct calls, reference-discordant benign calls and reference-discordant pathogenic calls are counted against ClinVar 2-star+ classifications (Section 2.2).
Supplementary Table S4 gives every step of every sweep.

<!-- table:thresholds:begin -->
| Threshold: default value | How it was set | Last step before the default value | First step beyond it | On v38 |
|:---|:---|:---|:---|:---|
| PM2 threshold, genes not read as recessive-only: 4e-5 | Swept in v10 under the 2015 table; the median of the published VCEP thresholds is 4e-5 (64 genes) | 2e-5 to 4e-5: +9 calls, +1 (9 per added error) | 4e-5 to 1e-4: +5 calls, +0 (no added error) | Recall 62.0% to 62.5% across strict absence to 5e-4 under the default scheme; 51.0% to 57.2% under the 2015 table on both sides |
| PM2 threshold, genes read as recessive-only: 7e-5 | Not from the SVI PM2 recommendation, which the code cited and which gives no frequency; equal to the Hearing Loss expert panel's recessive threshold. First swept here | 2e-5 to 7e-5: +4 calls, +0 (no added error) | 7e-5 to 2e-4: +7 calls, +0 (no added error) | Flat: recall 62.43% to 62.48% and agreement in direction within 27 calls from strict absence to 5e-4 |
| BS1 threshold: 5e-3 | Swept on a 1-in-10 sample, where the move from 1e-2 looked free; on the full set it cost 4 more reference-discordant benign calls for 2.3 points of benign-side recall (6,443 calls), and was kept (v12). The median of the published VCEP thresholds is 2.5e-4 | 1e-2 to 5e-3: +4,493 calls, +2 (2,247 per added error) | 5e-3 to 2e-3: +4,925 calls, +6 (821 per added error) | Knee at 1e-3 (99 per added error below it); 5e-3 is on the cautious side of the knee |
| BS2 homozygote threshold, every gene not read as dominant-only: 1e-3 | Swept on the full benchmark (12 Aug): no knee, so set on the prevalence of the commonest recessive conditions | 1e-2 to 1e-3: +256 calls, +0 (no added error) | 1e-3 to 1e-4: +1,203 calls, +3 (401 per added error) | No knee: 401, 272 and 487 calls per added error below the default value |
| BS2 carrier threshold, genes read as dominant-only: 3e-4, with an allele count of at least 5 | Added after review of one gene (SNRNP200) and swept on the full benchmark (v38) | - | - | Not flat: 24 to 47 reference-discordant benign calls and 18 to 24 reference-discordant pathogenic calls from off to 3e-3 |
| BP1 veto: 3 pathogenic missense variants in the gene | Added in v9 with the veto, after the round-2 review of BP1 firing in genes with a pathogenic missense spectrum; the count was never swept before this paper | 2 to 3: +16 calls, +1 (16 per added error) | 3 to 5: +81 calls, +3 (27 per added error) | At most 111 calls per added error at every step, so there is no knee and the default count is on the aggressive side of the stricter ones; a count of 1 gives +1,735 calls in direction and 127 fewer correct benign-side calls, and 2 fewer reference-discordant benign calls (1 fewer when the veto leaves the variant's own record out) |
| BP7 intronic distance: 300 nt | Set at the knee of a sweep (v14) | 200 to 300: +153 calls, +1 (153 per added error) | 300 to 500: +37 calls, +1 (37 per added error) | Knee confirmed: 153 to 819 calls per added error below 300, about 24 to 37 above |
| BP7 PhyloP threshold: 2.0 | In the first commit, with no source recorded (the 2015 standards say "not highly conserved" and give no number). First swept here | 1.5 to 2: +4,505 calls, +3 (1,502 per added error) | 2 to 3: +6,422 calls, +0 (no added error) | No knee: at least 1,500 calls per added error at every step; switching the test off adds 16,116 benign-side calls and 2 reference-discordant benign calls (Section 3.3) |
| Allele-number floor for BA1 and BS1: 2,000 | From the dataset-size condition of the BA1 recommendation, applied here to the site's allele number and to BS1 | - | - | Removing it adds 728 calls in direction and no opposite-direction call |
| Likely pathogenic from 6 points | Fixed by the point system (Tavtigian et al., 2020); varied here as a sensitivity analysis | - | - | 5 points: +10,675 pathogenic-side calls, +6,003 ClinVar-VUS and +17 benign-side called pathogenic-side (Table 4) |
| Pathogenic from 10 points | Fixed by the point system | - | - | 9 to 12 points move agreement in direction by at most 3 calls |
| BA1 at 5% | Fixed by the guideline; the median of the published VCEP thresholds is 1e-3 | - | - | Not varied |
<!-- table:thresholds:end -->

**The PM2 threshold now matters little, because of the combining rule.**
The classifier uses one PM2 threshold for genes that it does not read as recessive-only (dominant, both, or of unknown inheritance) and another for genes that it reads as recessive-only, and v38 takes that reading from the mode of inheritance that ClinGen curates: of the 2,924 genes of the ClinGen file it reads 1,310 as recessive-only, 746 as dominant-only, 169 as both and 699 as neither (Section 3.8), where v37, which took it from the names of the diseases, read 67, 50, 12 and 2,795.
When the first was moved from strict absence to 4e-5 in v10, the sweep of that time recorded a rise of 19 points in pathogenic recall, with a knee.
On v38 the same sweep is almost flat, and on v37 it was almost flat as well.
Across strict absence and every threshold from 1e-6 to 5e-4, recall of ClinVar's pathogenic-side variants moves from 61.98% to 62.48%, which is 470 variants of 93,812.
It reaches 62.47% at the default 4e-5, and the two steps above that change it by 5 and 5 variants.
Reference-discordant pathogenic calls rise from 21 at strict absence to 24 at 4e-5 and 25 at 2e-4, and agreement in direction is higher at strict absence (524,878 calls) than at the default threshold (524,380), because ClinVar-VUS variants that carry PM2 plus other evidence stay VUS.
The threshold for recessive-only genes, which the code attributed to the SVI's recommendation, which gives no frequency, and which equals the recessive threshold of the Hearing Loss expert panel, is flat as well: from strict absence to 5e-4, recall of pathogenic-side variants is between 62.43% and 62.48%, agreement in direction varies by 27 calls, and the numbers of reference-discordant benign calls and reference-discordant pathogenic calls do not change.
The two thresholds therefore act on nearly every gene and on the 1,310 genes read as recessive-only, and both sweeps are flat.
On v37 the second threshold reached 67 genes, and a flat sweep there said what a threshold does to the calls it reaches and not what it would do if inheritance were read from the curated field; v38 reads that field, and the sweep stays flat on the genes it now reaches (Section 3.8).
The reason the knee disappeared is the combining rule, not the threshold.
Under the 2015 combining table, which is one switch away in the same classifier, the same sweep is steep (Figure 1B).
Recall of pathogenic-side variants is 51.0% at strict absence, 54.8% at 1e-6, 56.6% at 1e-5, 57.1% at 4e-5 and 57.2% at 1e-4.
Raising the threshold from strict absence to 4e-5 gains 5,700 of 93,812 variants (6.1 points) for 7 more reference-discordant pathogenic calls.
The sweep that v10 recorded under this table, 37.8% at strict absence and 56.8% at 4e-5, a rise that it took for a knee, was reproduced by v37 (38.8% and 56.7%, a gain of 16,785 variants, 17.9 points); v38 rises at every step as v37 does and by about a third as much (51.0% and 57.1%, a gain of 5,700 variants, 6.1 points).
It is still the steeper sweep: it moves recall 12 times as far as the same sweep under the default scheme (5,700 against 470 variants).
Why v38 flattens it is a reading that the firing counts support and that no run has tested, because no run sweeps both PM2 thresholds under the 2015 table: at strict absence PM2 still fires on 229,154 variants under v38 against 133,152 under v37, and the 1,310 genes that v38 reads as recessive-only keep their own threshold of 7e-5, which this sweep does not change.
Under that table a lone PVS1 matches no row, and a null variant whose only companion is PM2 reaches Likely pathogenic only through the rule that the SVI added with its PM2 recommendation, one Very Strong criterion plus one Supporting criterion (ClinGen SVI, 2020).
Of the 51,705 pathogenic-side calls that the table makes with PM2 and would not make without it, 46,209 (89%) rest on PVS1 and PM2 and nothing else.
The point scale of the default scheme counts a lone PVS1 as 8 points, so under it only 5,496 calls depend on PM2, and none of them is a PVS1-and-PM2 pair.
A threshold can move only the calls that need the criterion it gates, and the combining rule decides how many do.

**The BS1, BS2 and BP7 thresholds trade benign-side calls for reference-discordant benign calls, and each has a different shape.**
All three pass through the default point (Figure 1A), and Table 3 gives what each step buys and costs.
BP7's intronic distance reproduces the knee found in v14.
Moving it from 200 to 300 nucleotides buys 153 correct benign-side calls for one more reference-discordant benign call, and moving it from 300 to 500 buys 37 for one more, with 48 for two more beyond that and 440 for fourteen with no bound at all.
The default 300 is at the corner: below it a step is worth 153 to 819 calls per reference-discordant benign call, above it about 24 to 37.
Turning the extension off is not free either: it gives up 9,004 correct benign-side calls to avoid 11 reference-discordant benign calls (the step from 0 to 50 nucleotides).
The BS2 homozygote threshold, which applies to every gene that is not read as dominant-only, has no knee, as in the original sweep.
The steps from 1e-3 to 1e-6 are worth 401, 272 and 487 calls per reference-discordant benign call, so no setting is picked out by the data and the choice rests on how the two errors are weighed.
The BS1 threshold is the one that moved.
Between 2e-2 and 5e-3 (default), lowering the threshold buys 9,910 benign-side calls at the cost of 2 reference-discordant benign calls.
Below the default value the exchange rate is 821 calls per reference-discordant benign call from 5e-3 to 2e-3, 562 from 2e-3 to 1e-3, and then 99 from 1e-3 to 5e-4 and 79 from 5e-4 to 2e-4.
By the knee rule that fixed BP7's bound, which looks for the point past which the exchange rate falls from hundreds of calls per reference-discordant benign call to about a hundred or fewer, the knee of the BS1 curve is at 1e-3, and the default 5e-3 sits on its cautious side: 2.7 points of benign-side recall (7,735 calls) and 11 reference-discordant benign calls (24 against 35) separate them.
The published thresholds point the same way.
In the table of ClinGen VCEP specifications we extracted, the median BS1 threshold is 2.5e-4 (107 genes with a usable threshold, 102 of them stricter than our 5e-3), and the median BA1 threshold is 1e-3 against 5e-2 (115 genes, 111 stricter).
The median PM2 threshold is 4e-5 (64 genes, range 1.3e-6 to 5e-4), which is the value we ship for every gene that is not read as recessive-only.
That median leaves out the panels that give no single number: of 137 PM2 entries the extraction took 64 as thresholds, and set aside 37 that require absence, 14 with no threshold, 9 whose threshold depends on inheritance (the Hearing Loss panel's recessive threshold of 7e-5 among them), 6 in superscript notation, 4 where two panels disagree and 3 that are ambiguous or unspecified (Appendix Table A5).
Across BA1, BS1 and PM2 together the extraction read 411 criterion entries and took 286 as thresholds (Table 1).
We have not changed the default value in this work.
A global frequency threshold is a poor stand-in for a threshold that depends on a disorder's prevalence and penetrance, which is what "greater than expected for the disorder" asks for, and that is a decision for a reviewer, not for a sweep.

**BP7's conservation test gives up more benign-side calls per reference-discordant benign call avoided than any other threshold in the ranges swept, and it had never been swept.**
BP7 asks that the nucleotide not be highly conserved, which the classifier reads as a PhyloP score below 2.0, a number from the first commit with no source recorded.
Raising the threshold from 1 to 1.5 gains 5,423 correct benign-side calls with no added reference-discordant benign call, and from 1.5 to the default 2 gains 4,505 for 3 more (1,502 per added error).
Beyond the default value, raising it to 3 gains 6,422 for none, to 4 gains 3,440 for one, and switching the test off gains 6,254 for one more.
Every step buys at least 1,500 calls per added reference-discordant benign call, so unlike the BS1 threshold the curve has no knee, and the rate for switching the test off, about 8,100 correct calls per added reference-discordant benign call, rests on two events.
Switching the test off, which follows the caution of Walker et al. (2023) against a conservation test without an empirical justification, adds 16,116 correct benign-side calls and 127 ClinVar-VUS variants called benign-side, raises agreement in direction by 15,989 calls (2.4 points) and raises reference-discordant benign calls from 24 to 26.
These gains are measured against ClinVar's benign-side labels, which submitters may themselves have assigned by applying BP4 and BP7, so for this criterion the reference is not independent of the evidence.
Asking that SpliceAI be at most 0.1, as BP4 does, removes 338 correct benign-side calls from v38, and asking it on top of switching the test off removes 437, neither changing an error, so the two together, which is BP7 as the recommendation describes it, give 539,932 calls in direction (80.1%), 15,552 more than v38, for the same two additional reference-discordant benign calls.
We have not changed the classifier: the two variants (NF1 and MAPT) are ClinVar-pathogenic and have PM2 and BP4 beside BP7, and it is for the maintainers to weigh two reference-discordant benign calls against 16,000 correct calls (Section 4.5).
On the 107,528 variants that were not in the April set, BP7 as the recommendation describes it adds 4,859 calls in direction and no reference-discordant benign call (Section 3.7).

**The count of BP1's veto had never been swept, and a stricter count calls fewer ClinVar-VUS variants benign-side and avoids two reference-discordant benign calls, one of which rests on the variant's own record.**
BP1 is withheld in a gene that has at least 3 pathogenic missense variants in the ClinVar protein index, a count that was set with the veto in v9 and never varied (the veto itself avoids 28 reference-discordant benign calls relative to none).
At a count of 1, BP1 is withheld in every gene with one pathogenic missense record: agreement in direction rises by 1,735 calls, reference-discordant benign calls fall from 24 to 22, reference-discordant pathogenic calls do not change, and 127 correct benign-side calls are given up, while the benign-side calls on ClinVar-VUS variants fall by 1,862.
The rise is those 1,862 variants staying VUS, which is where ClinVar leaves 44% of its records, less the 127 correct calls given up; it is not a finding that the benign-side calls on them were wrong (Section 3.5).
The two reference-discordant benign calls that a count of 1 removes are the BP1 pair of Section 3.6 (IQSEC1 and PRR12): the protein index holds one pathogenic missense record for PRR12, the variant itself, so its removal rests on its own record, and two for IQSEC1, the variant and p.Thr343Met, so its removal rests on its own record at a count of 2 and not at a count of 1.
With the variant's own record left out of the count as well (Section 3.4), a count of 1 therefore removes IQSEC1 only (23 reference-discordant benign calls, with the same 1,735 calls and 127 correct benign-side calls) and a count of 2 neither (24).
At a count of 2 without that change the figures are 850 more calls in direction, 23 reference-discordant benign calls (IQSEC1 removed) and 16 correct benign-side calls given up.
Above the default value each step buys correct benign-side calls at the price of reference-discordant benign calls: a count of 5 gives 3 more, 10 gives 5 more and 25 gives 12 more.
No step gives more than 111 calls per added reference-discordant benign call, so the exchange rate never falls from hundreds of calls to about a hundred, the curve has no knee in the sense of Section 2.5, and the default count of 3 lies on its aggressive side; a count of 1 or 2 is the cautious choice, though what it avoids is two reference-discordant benign calls on two variants, one of which (PRR12) rests on its own record (Table 3).
Both are among the 48 and were found in the set on which the count was swept, and on the new variants a count of 1 avoids no reference-discordant benign call (Section 3.7).
We have not changed it, and the count bears on the finding of Section 3.4 that the veto counts the variant's own record.

**The BS2 carrier threshold for genes read as dominant-only is a trade, and the default is the last step that adds no error.**
On the 746 genes that v38 reads as dominant-only the threshold is not flat.
Switching it off, so that an allele count of at least five alone suffices, gives 47 reference-discordant benign calls and 18 reference-discordant pathogenic calls, and 3e-5 and 1e-4 give 31 and 27 reference-discordant benign calls, against 24 at the default 3e-4 and at 1e-3 and 3e-3.
From 3e-3 to the default each step adds benign-side calls and no error (448 and 1,012), and below it a step buys 715, 678 and 315 calls per added reference-discordant benign call.
The default 3e-4 recovers 439 pathogenic-side calls relative to no threshold (58,167 to 58,606) and gives up 9,897 benign-side calls (188,785 to 178,888), and 1e-3 recovers 4 more pathogenic-side calls (58,610).
With the threshold off, a floor of 20 alleles gives 33 reference-discordant benign calls, a floor of 100 gives 26 and a floor of 1,000 gives 24.
We have not changed the default.
On v37 the same sweep was flat, because the threshold reached the 50 genes that the names of the diseases made dominant-only, and its reference-discordant benign and reference-discordant pathogenic calls were 24 and 25 at every setting (Section 3.8).

**On the pathogenic side the points that define the classes move the result more than any numeric threshold that we swept on the default classifier, and the evidence sits in an awkward place on that scale.**
Figure 2A shows, for variants that met no benign evidence, where ClinVar's classes fall on the pathogenic points the classifier reached (Appendix Table A6).
Variants with 0 or 1 point (24,535 and 106,915 of them) are 4.1% and 7.8% ClinVar pathogenic-side.
At 2 points the share is 10.6%, at 3 points 21.9%, and at 4 points 25.7%.
At 5 points it jumps to 64.1% (10,637 of 16,594), with 35.8% ClinVar VUS; at 6 points the share is 59.0%; from 7 points it is 94.9% or more, reaching 99.5% at 10.
The 5-point group is a cluster: PM1 with PM2 and PM5, PVS1 at Strong strength with PM2, and PP3 at Strong strength with PM2 each reach exactly 5, one point short of Likely pathogenic.
Moving the Likely pathogenic cutoff from 6 to 5 points, in a real run of the classifier with that one constant changed (Table 4), raises recall of ClinVar's pathogenic-side variants from 62.5% to 73.9% (10,675 more variants).
It does so at the price of 6,003 more ClinVar-VUS variants and 17 more ClinVar benign-side variants called pathogenic-side, so 1.8 pathogenic-side variants are gained for each variant added that ClinVar does not call pathogenic-side (whether the VUS among them are errors is the question of Section 3.5).
Agreement in direction rises by 4,672 calls, and exact agreement falls by 3,670.
A cutoff of 4 reaches 76.7% recall, but 17,464 ClinVar-VUS and 129 benign-side variants are then called pathogenic-side, and agreement in direction falls below the default value.
Raising the cutoff to 7 or 8 gives up 3,254 or 3,532 pathogenic-side calls and avoids 3 reference-discordant pathogenic calls, and at 9 a lone PVS1 (8 points) is no longer Likely pathogenic, recall falls to 53.2% and reference-discordant pathogenic calls to 11.
No cutoff changes the number of reference-discordant benign calls, which is 24 in every run.
Tavtigian's framework places a posterior of 0.90 at 6 points and 0.81 at 5 under its prior of 0.10, and the benchmark's own prior among defined classes is 0.25, so the comparison with ClinVar cannot tell which cutoff is calibrated; it can say what each costs.
The Pathogenic cutoff of 10 points does not move agreement in direction: moving it between 9 and 12 changes that agreement by at most 3 calls, because the cutoff separates two tiers of the same direction.
It does move exact agreement, which is 463,222 at 9 points against 426,022 at 10.
Of the 58,200 Likely pathogenic calls, 47,055 (81%) are at exactly 9 points, 46,209 of them PVS1 with PM2 and nothing else, and 41,549 of those 47,055 are ClinVar Pathogenic.
The tier boundary therefore falls one point above where a lone null variant with PM2 lands, and we have kept it at the framework's 10, because what lifts such a variant to Pathogenic in laboratory practice is case-level evidence that a single-sample VCF does not carry.

**Table 4. The two point cutoffs of the pathogenic side, by real runs.**
A point is the unit of the scale that combines pathogenic criteria (Section 2.1): evidence of Supporting, Moderate, Strong or Very Strong strength counts 1, 2, 4 or 8 points, and a variant with 6 to 9 points is Likely pathogenic and one with 10 or more is Pathogenic.
Each row is a run of v38 with one of the two cutoffs changed, which the code holds as constants (Section 2.6); the offline replay of Section 2.7 reproduces every row cell for cell except the first two, where it differs by 3 and 1 variants.
Rates are over 673,660 records.
ClinVar P/LP, VUS and B/LB are classes of ClinVar 2-star+ records of the release of 4 April 2026 (Section 2.2).

<!-- table:cutoffs:begin -->
| Setting | ClinVar P/LP called P/LP (recall) | ClinVar VUS called P/LP | ClinVar B/LB called P/LP | ClinVar P/LP called LB/B | Same-direction agreement | Exact agreement |
|:---|---:|---:|---:|---:|---:|---:|
| Likely pathogenic from 4 points | 71,964 (76.7%) | 17,464 | 129 | 24 | 523,871 (77.8%) | 415,098 (61.6%) |
| Likely pathogenic from 5 points | 69,281 (73.9%) | 9,611 | 41 | 24 | 529,052 (78.5%) | 422,352 (62.7%) |
| Likely pathogenic from 6 points (default) | 58,606 (62.5%) | 3,608 | 24 | 24 | 524,380 (77.8%) | 426,022 (63.2%) |
| Likely pathogenic from 7 points | 55,352 (59.0%) | 1,349 | 21 | 24 | 523,385 (77.7%) | 427,608 (63.5%) |
| Likely pathogenic from 8 points | 55,074 (58.7%) | 1,333 | 21 | 24 | 523,123 (77.7%) | 427,588 (63.5%) |
| Likely pathogenic from 9 points | 49,918 (53.2%) | 1,164 | 11 | 24 | 518,136 (76.9%) | 426,506 (63.3%) |
| Pathogenic from 9 points | 58,609 (62.5%) | 3,608 | 25 | 24 | 524,383 (77.8%) | 463,222 (68.8%) |
| Pathogenic from 10 points (default) | 58,606 (62.5%) | 3,608 | 24 | 24 | 524,380 (77.8%) | 426,022 (63.2%) |
| Pathogenic from 11 points | 58,606 (62.5%) | 3,608 | 24 | 24 | 524,380 (77.8%) | 422,929 (62.8%) |
| Pathogenic from 12 points | 58,606 (62.5%) | 3,608 | 24 | 24 | 524,380 (77.8%) | 422,891 (62.8%) |
<!-- table:cutoffs:end -->

**The default scheme is points on the pathogenic side and the 2015 table on the benign side, and each half was chosen over a measured alternative.**
Three schemes can be set side by side on the recorded criteria of v38 (Appendix Table A2): the default scheme, the 2015 table on both sides, and the point scale on both sides.
The first two are real runs, the second through the setting `use_point_system = false`, and the third has no setting in the classifier and is a replay of the recorded criteria (Section 2.7).
The default scheme has 426,022 exact calls and 524,380 in direction, recall of ClinVar's pathogenic-side variants of 62.5% and of its benign-side variants of 62.9%, 24 reference-discordant benign calls and 24 reference-discordant pathogenic calls.
The 2015 table on both sides has 424,946 exact calls and 519,480 in direction, recall of 57.1% of ClinVar's pathogenic-side variants and 62.9% of its benign-side variants, and 24 reference-discordant benign calls and 14 reference-discordant pathogenic calls, so it differs from the default scheme on the pathogenic side only, where it gives up 4,900 calls in direction and 5.4 points of recall to avoid 10 reference-discordant pathogenic calls (Section 3.4).
On the benign side the choice of scheme matters more than the choice of floor.
Scoring both sides on the point scale with Likely benign from -1 points, as in the published point system, raises recall of ClinVar's benign-side variants from 62.9% to 74.1%, but raises reference-discordant benign calls from 24 to 371, of which 205 are a lone BP4, and lowers recall of ClinVar's VUS variants from 97.2% to 69.6%.
Tightening the floor to -2 gives 49.9% recall and 64 reference-discordant benign calls, and to -3 gives 38.4% and 30.
The table therefore does better than the floors of -2 and -3 on recall and on reference-discordant benign calls alike, and the floor of -1 is a trade of 347 more reference-discordant benign calls for 11.2 points of recall (Figure 2C).
These three point scales are replays of the recorded criteria and not real runs (Section 2.7), and each has fewer exact calls (350,142, 367,909 and 364,989) and fewer calls in direction (474,950, 474,029 and 455,333) than the default scheme.

**An allele-number floor whose source is narrower than its use bought nothing.**
BA1 and BS1 require an allele number of at least 2,000 at the variant's site in gnomAD.
The figure comes from the BA1 recommendation, which asks for an allele frequency above 5% in a general continental population dataset of at least 2,000 observed alleles (Ghosh et al., 2018).
The project's documentation attributed it to the March 2024 communication on gnomAD v4, which states no allele-number floor and recommends the grpmax filtering allele frequency (ClinGen VCEP Review Committee, 2024).
Here it is applied to the allele number at the variant's own site, and to BS1 as well, and neither use is stated by either source.
Removing it raises same-direction agreement by 728 calls and changes no opposite-direction call (Section 3.4).

### 3.4 What each decision and each source buys, and how much of the agreement is circular

Figure 3 shows, for every ablation of Section 2.6, the change from v38 in agreement in direction, in recall of ClinVar's pathogenic-side variants and in opposite-direction calls.
A setting that lowers agreement with ClinVar is not for that reason wrong, because ClinVar rewards calls and a safeguard withholds some.
We therefore read each safeguard by what it removes in errors and what it costs in calls, as in Section 3.3.

**Counting a variant's own ClinVar record would add 14.3 points of pathogenic-side recall, and removing every ClinVar-derived input costs 0.78 points of agreement and raises the opposite-direction calls from 48 to 75.**
The default classifier excludes a variant's own ClinVar record from the evidence that PS1 and PM1 read (PM5 cannot match it, Section 2.3).
Letting it count, as the runs before v9 did, raises agreement in direction by 13,445 calls (2.00 points) and recall of pathogenic-side variants from 62.5% to 76.8%, which is the figure a self-referential evaluation would have reported, and removes 7 of the 48 opposite-direction calls.
Removing the ClinVar protein and splice index, which holds the classifications of other variants that PS1, PM1, PM5 and BP1 read, lowers agreement by 5,236 calls (0.78 points) and pathogenic-side recall by 2.1 points.
Removing the allele-level ClinVar record as well changes one further call.
The dependence on ClinVar-derived inputs is therefore small in agreement and larger in errors: without the index, reference-discordant benign calls rise from 24 to 52, and from 24 to 54 with the allele-level record as well (75 opposite-direction calls against 48), and the 28 are the same variants as appear when BP1's veto is switched off, because BP1 is then no longer withheld in genes where missense variants are known to be pathogenic.
The veto is protective, but it protects by reading ClinVar's other classifications in the gene, which in a benchmark against ClinVar is part of the circularity that this section sets out to measure.
The veto itself counts the variant's own record, which PS1 and PM1 do not: in a measurement build in which the veto leaves it out, reference-discordant benign calls are 27 and not 24, so 3 of the 28 that the veto avoids rest on the variant's own record (FBXO28, PCGF2 and SBF1, each a ClinVar-pathogenic missense variant whose own entry is the third of the three that the veto needs).
A benchmark against ClinVar cannot separate that dependence from the classifier's own skill, because other variants' ClinVar classifications are evidence in the guideline's own terms (PS1 and PM5 are defined by them).
Two further adjustments also read the variant's own ClinVar record or followed the reviewer's rows: ClinVar's low-penetrance gate (1 reference-discordant benign call when removed) and the curated frequency-exception list (6).
Removed together with the veto's own-record count, the three raise reference-discordant benign calls from 24 to 34, exactly the sum of the three alone (3, 1 and 6), and opposite-direction calls from 48 to 58, so 58 is the figure when the adjustments that read the own record or followed the reviewer's rows are removed together, and the veto itself and the other changes that followed the review remain.
Three of the ten added calls are Benign calls on ClinVar-pathogenic variants (one SERPINA1 variant once the low-penetrance gate is gone, and two G6PD variants once the exception list is), so the zero for Benign calls that Section 3.5 reports rests on two of the 15 curated alleles (both G6PD) and one ClinVar-derived gate.

**Three sources carry the benign side and one carries most of the pathogenic side.**
Removing SpliceAI lowers agreement in direction by 108,084 calls (16.0 points) and recall of benign-side variants from 62.9% to 24.7%.
Removing PhyloP lowers it by 77,774 calls (11.5 points; benign-side recall 35.4%), and removing gnomAD by 48,029 (7.1 points; 46.2%).
Without SpliceAI, 199,869 of the 248,740 BP4 firings are lost, and BP7, which needs SpliceAI and PhyloP together, is lost entirely (132,615 firings).
gnomAD feeds BA1, BS1, BS2 and PM2, and all four stop firing without it.
On the pathogenic side, removing the ClinGen gene-disease file lowers agreement by 44,685 calls (6.6 points) and recall of pathogenic-side variants from 62.5% to 15.7%, because for 48,540 of the 61,334 firings of PVS1 at any strength (79%) a ClinGen curation of Moderate strength or better is the only route by which the gene is taken to act by loss of function, the others being gnomAD constraint and a curated mechanism; the file also supplies the inheritance that PM2, BS2 and PM3 read, so the loss is not all PVS1's.
Removing gnomAD gene constraint takes out PP2 and BP1 entirely and 721 PVS1 firings and lowers pathogenic-side recall by 4.0 points, and it raises agreement in direction by 1,933 calls, because variants that lose their criteria become VUS and 44% of the reference is VUS.
Removing REVEL, the source specific to missense variants, lowers agreement by 0.88 points, and removing RepeatMasker (BP3) by 52 calls.

**The same two sources carry most of the reference-discordant benign calls.**
Removing SpliceAI or PhyloP lowers reference-discordant benign calls from 24 to 4, because BP7 is behind 20 of the 24.
BP7 is a criterion of high yield and high cost: 132,202 of its 132,615 firings (99.7%) are on ClinVar benign-side variants, and without it 77,774 calls in direction are lost to avoid those 20 reference-discordant benign calls, about 4,000 calls for each.
Section 3.3 shows that much of that price is paid for the conservation test and not for BP7's other conditions, and Section 3.6 returns to the 20.

**Most design decisions pay for themselves, and three do not clearly do so.**
Four safeguards remove errors and cost no calls: the BP1 veto (keeping it gains 4,233 calls and avoids 28 reference-discordant benign calls, and a stricter count would avoid 2 more, or 1 with the own record left out, Section 3.3), the gate on disease mechanism for PVS1 (146 calls and 12 reference-discordant pathogenic calls), the curated list of 15 frequency-exception alleles (6 reference-discordant benign calls) and the ClinVar low-penetrance gate (1).
Reading absence from gnomAD as rarity for PM2 gains 1,499 calls and no error.
Three more remove errors at a price in calls: the gene-disease validity gate costs 193 calls to avoid 10 reference-discordant pathogenic calls (about 19 calls per error), the homologous-gene list 376 calls to avoid 1 reference-discordant benign call (about 376 per error), and the filtering allele frequency 1,917 calls to avoid 8 reference-discordant benign calls (about 240 per reference-discordant benign call), while adding 4 reference-discordant pathogenic calls.
Their prices are all below the 821 calls per added error at which we kept the BS1 threshold where it was, so by that standard each is worth paying.
The standard was not applied uniformly, which Section 4.4 notes.
Three cost more than they remove on this benchmark.
The gnomAD region flags (segmental duplications and low-complexity tracts) withhold the frequency criteria at a price of 3,473 calls to avoid one reference-discordant pathogenic call, the allele-number floor withholds 728 calls to avoid none, and PM1's requirement that the window hold no benign variation withholds 272 to avoid none.
We have not changed any of them, because safeguards whose targets are rare are poorly tested by a benchmark of consensus classifications.
The allele-number floor is the first candidate for removal: neither source it can be traced to states it in the form used here (Section 3.3), and the data do not support it.
Two opt-in settings change little: capping PP3 at Moderate raises agreement by 715 calls and lowers pathogenic-side recall by 0.34 points with no change in errors, and the exact 50-nucleotide NMD rule for PVS1 loses 603 calls to avoid one reference-discordant pathogenic call.

**Two settings would trade errors for recall in the other direction.**
PM2 at Moderate, as in the 2015 standards, raises recall of pathogenic-side variants by 11.1 points (10,432 more variants) and agreement in direction by 4,567 calls, for 14 more reference-discordant pathogenic calls (about 326 calls per error).
The 2015 table on both sides loses 4,900 calls and 5.4 points of recall to avoid 10 reference-discordant pathogenic calls, which is the combining-rule effect of Section 3.3.
Neither is a clear improvement by the standard of this paper, and the first contradicts the ClinGen recommendation that PM2 be used at Supporting strength.

### 3.5 Discordance of v38 with ClinVar 2-star+ classifications

On the benchmark snapshot, v38 reproduces ClinVar's class for 426,022 of 673,660 variants (63.2%, 95% interval 61.7 to 64.6) and its direction for 524,380 (77.8%, 77.0 to 78.6; gene-cluster bootstrap over 10,232 groups of variants, which are 10,231 genes and one group of the 4 records that name none).
Table 5 divides every record by how the call relates to ClinVar's class (Figure 4A; Supplementary Table S5 holds the matrices and the headline rates).

**Table 5. Every ClinVar 2-star+ record, by how the v38 call relates to ClinVar's class.**

<!-- table:decomposition:begin -->
| Kind | Cells | Variants | Share of all records |
|:---|:---|---:|---:|
| Same class as ClinVar | the diagonal | 426,022 | 63.2% |
| Same direction, different tier | P and LP, LB and B | 98,358 | 14.6% |
| &nbsp;&nbsp;ClinVar Pathogenic called Likely pathogenic |  | 48,277 | 7.2% |
| &nbsp;&nbsp;ClinVar Benign called Likely benign |  | 47,608 | 7.1% |
| ClinVar P/LP called VUS | a pathogenic-side class, no call made | 35,182 | 5.2% |
| ClinVar LB/B called VUS | a benign-side class, no call made | 105,256 | 15.6% |
| ClinVar VUS called P/LP |  | 3,608 | 0.5% |
| ClinVar VUS called LB/B |  | 4,804 | 0.7% |
| Opposite direction |  | 48 | 0.007% |
| &nbsp;&nbsp;ClinVar P/LP called LB/B (reference-discordant benign calls) |  | 24 | 0.004% |
| &nbsp;&nbsp;ClinVar LB/B called P/LP |  | 24 | 0.004% |
| No call | records with no alternative allele (ALT is a dot) | 382 | 0.06% |
| All records |  | 673,660 | 100% |
<!-- table:decomposition:end -->

**The discordance is almost entirely one side being undecided, and contradictions are rare.**
Of the 149,280 records that are not in the same direction as ClinVar, 140,438 (94.1%) are a pathogenic-side or benign-side ClinVar class that v38 leaves as VUS, 8,412 (5.6%) are the reverse, 382 (the records with no alternative allele, Section 2.2) have no call, and 48 (0.03%) are opposite.
The opposite-direction rate is 0.0071% of records (95% interval 0.0050 to 0.0096% by gene-cluster bootstrap; 0.013% of the records whose ClinVar class is pathogenic-side or benign-side), made up of 24 reference-discordant benign calls (0.026% of ClinVar's 93,812 pathogenic-side variants) and 24 reference-discordant pathogenic calls (0.008% of its 284,550 benign-side variants).
A further 98,358 records (14.6%) are in the right direction at the wrong tier, and 95,885 of them are ClinVar Pathogenic called Likely pathogenic (48,277) or ClinVar Benign called Likely benign (47,608): exact agreement measures the tier boundary more than the direction.
The agreement is also asymmetric in what it says about each call.
Of the 62,238 variants called pathogenic-side, 94.2% (58,606; 92.6 to 95.5) are ClinVar pathogenic-side and 5.8% ClinVar VUS, and of the 183,716 called benign-side, 97.4% (178,888; 96.5 to 98.1) are ClinVar benign-side (Figure 4C).
Recall is much lower than precision on both sides: 62.5% of ClinVar's pathogenic-side variants (58.2 to 66.6) and 62.9% of its benign-side variants (61.2 to 64.5) receive a call in the same direction, against 97.2% of its VUS variants (96.5 to 97.7).
The wide interval of the first is real: recall is a property of the gene.
Among the 1,273 genes with at least ten ClinVar pathogenic-side variants, 233 (18%) have none called pathogenic-side and 231 (18%) have 90% or more (Appendix Table A7).

**Agreement in direction does not separate cautious from aggressive calling.**
Calling every variant VUS agrees with ClinVar in direction for 43.8%, and rules that read the inputs of some ACMG criteria but apply no combining rule, gene-level gate or safeguard (Supplementary Table S6) reach v38's agreement and more only by making many more contradictions (Figure S4).
Calling a variant benign-side when its population-maximum allele frequency is 5e-3 or more gives 61.6% with 159 contradictions, and at 1e-3 (the median published BA1 threshold) 64.0% with 813; a higher threshold has fewer contradictions and agrees less, 52 contradictions at 59.7% for 1e-2 and 22 at 57.1% for 2e-2, so with the number of contradictions that v38 has the rule agrees for about 18 points less.
Adding a pathogenic-side call for a stop-gained, frameshift, canonical splice or start-lost variant that is absent from gnomAD or below 1e-4 gives 72.8% with 865, and adding a benign-side call for a quiet variant with SpliceAI below 0.2 gives 84.6%, 6.7 points above v38, with 1,728 contradictions (1,676 reference-discordant benign calls and 52 reference-discordant pathogenic calls), 36 times as many as v38.
Without the frequency rule, the last two rules give 77.9%, 0.1 points above v38, with 948 contradictions (896 reference-discordant benign calls and 52 reference-discordant pathogenic calls), 20 times as many, and no threshold of the first rule brings the three rules below 948 contradictions.
The last rule is not a bargain by the standard of Section 2.5 either: over the two-rule version at 1e-3 it adds 86,510 correct benign-side calls for 863 more reference-discordant benign calls, 100 per error.
The population maximum used by these rules is a point estimate over ten groups, so it is more aggressive than the filtering allele frequency that BA1 tests (in v38 a point estimate in place of that frequency costs 8 more reference-discordant benign calls, Section 3.4), but it is the threshold that sets the contradictions of the first rule (22, 52, 159 and 813 at 2e-2, 1e-2, 5e-3 and 1e-3).
The rules are not worse than v38 at everything: the two-rule version has a higher recall of ClinVar's pathogenic-side variants (66.0% against 62.5%) and a slightly higher share of its pathogenic-side calls that ClinVar calls pathogenic-side (95.4% against 94.2%), at 52 reference-discordant pathogenic calls against 24, which is in part the PVS1 gap below seen from the other side: of the 9,970 pathogenic-side calls that the second rule makes correctly and v38 leaves as VUS, 5,749 carry no PVS1 in v38 and 4,221 carry it at a lower strength or held back, and 41 of its 52 reference-discordant pathogenic calls are VUS in v38.
So agreement alone cannot say which of two classifiers is the more cautious, and v38 is better read with its contradictions beside its agreement.
The tiers of v38 are informative in the same way, with their denominators: none of its 4,038 Pathogenic calls (upper end of the exact 95% interval 0.091%) and none of its 67,733 Benign calls (0.0054%) contradicts ClinVar, against 24 of 58,200 Likely pathogenic calls (0.041%) and 24 of 115,983 Likely benign calls (0.021%).
The zero for Benign is informative, because 14 would be expected at the Likely benign rate, and the zero for Pathogenic is not, because 1.7 would; both hold for this in-sample set only, and the zero for Benign rests on two of the 15 curated alleles (both G6PD) and one ClinVar-derived gate, since three Benign calls contradict ClinVar when they are removed (Section 3.4).
Its 4,038 Pathogenic calls are ClinVar pathogenic-side for 4,017 (99.5%) and VUS for 21, and its 67,733 Benign calls are ClinVar benign-side for 67,408 (99.5%) and VUS for 325; the 24 reference-discordant pathogenic calls are all Likely pathogenic and the 24 reference-discordant benign calls all Likely benign.

**Review status separates the classes, and the highest tier is one gene.**
The 3-star records, which are expert-panel classifications, are called pathogenic-side more often than the 2-star records (67.8% against 61.7% of their pathogenic-side variants), and called benign-side much less often (37.6% against 63.4% of their benign-side variants, with 62.2% left as VUS).
Expert-panel benign-side classifications can rest on disease-specific frequency thresholds, and the published thresholds are stricter than ours (Section 3.3), which is a possible reason for the difference and one that this benchmark cannot establish.
All 23 four-star records (practice-guideline classifications) are variants of CFTR, and v38 calls 22 of them VUS (Figure 5A).

**Pathogenic-side recall is set by the consequence class.**
Variants called pathogenic-side are almost all null variants: 54,494 of the 58,606 correct pathogenic-side calls (93.0%) are frameshift, nonsense or canonical splice variants, and 4,109 are missense (Figure 5A, Supplementary Table S5).
Among ClinVar's pathogenic-side variants, 88.4% of canonical splice variants and 82.6% of frameshift and nonsense variants are called pathogenic-side, against 17.8% of missense variants, 0.4% of in-frame variants and none of the start-loss, stop-loss, splice-region, intronic, synonymous, UTR and other variants.
Among ClinVar's benign-side variants the picture is the reverse: 79.3% of intronic, 66.0% of synonymous, 63.8% of UTR and 51.8% of missense variants are called benign-side (Figure 5B).

**Where the 35,182 pathogenic-side variants called VUS come from.**
They are not an even sample.
By the evidence that was met, 10,693 (30.4%) reached 5 points, one short of Likely pathogenic, 10,533 (29.9%) reached 0 or 1 point, and 13,939 reached 2 to 4 points; 17 reached 6 points or more and were held at VUS by a veto or by benign evidence.
The 5-point group is the cluster of Section 3.3 (PM1 with PM2 and PM5, PVS1 at Strong strength with PM2, and PP3 at Strong strength with PM2).
By consequence, 18,992 (54.0%) are missense variants, 9,136 (26.0%) frameshift or nonsense, 2,041 (5.8%) in flanking or non-coding sequence, and the rest splice-region, canonical splice, in-frame, start-loss and stop-loss, intronic, UTR and synonymous variants.
The concentration in genes is high: the 10 genes with most of them hold 11.7%, the top 100 hold 40.5% and the top 500 hold 74.7%, of 2,717 genes with any (CFTR 778, APC 712, LDLR 517, PAH 453, GCK 344, ABCA4 298, GLA 264, BRCA1 255, NF1 253 and SACS 246).

**A gap in PVS1 accounts for more than half of the truncating variants that are missed, and an incomplete input explains a small part of it.**
Of the 52,584 ClinVar pathogenic-side stop-gained and frameshift variants, 47,400 carry a PVS1 criterion at some strength and 5,184, in 625 genes (620 by the symbol of the annotation, which differs from ClinVar's first gene symbol for 55 of the variants), carry none (Supplementary Table S7).
Those 5,184 are 57% of the 9,136 that are called VUS.
PVS1 asks three questions in order.
The validity gate asks whether every ClinGen curation of the gene is weak, the mechanism gate whether the curated mechanism excludes loss of function, and the gene must then be taken to act by loss of function, which the classifier reads from one of three routes: gnomAD constraint (pLI at least 0.9 or LOEUF at most 0.35), a curated loss-of-function mechanism, or a ClinGen curation of Moderate strength or better.
The second route is nearly empty: of the 26 genes in the default table of mechanisms, 20 are gain of function, one is dominant negative and five name both loss and gain of function, and those five are each licensed by constraint or by a curation already.
The validity gate had one leak that changed no call in v37: of the 505 genes whose ClinGen curations are all weak, one (NSMCE3) escaped it because its label ends in a semicolon, which the database builder reads as a separator, so that the label and the tag landed in different phenotype strings and the label, with no tag, counted as a disease of its own, and it has no variant in the benchmark.
v38 closes it, because it drops an untagged piece where a tagged one exists (Section 2.10).
For 4,696 of the 5,184 variants (91%, in 575 genes) none of the three routes is met: the gene tolerates heterozygous loss in gnomAD, has no curated mechanism and has no curation of Moderate strength or better in the file.
For 463 (in 38 genes) the validity gate withholds PVS1, and for 5 (in 4 genes) the mechanism gate does, and 20 meet a route and are not graded as PVS1 for another reason.
Absence has the effect here that it had in the first version of the validity gate (v10, Table 1), which blocked PVS1 where a gene had no curation and which v11 replaced with a gate that blocks only on curated-negative evidence: a gene with no curation is not blocked now, but with no other route to loss of function it is not licensed either.
All 76 variants in genes that gnomAD constrains (RASA1 40, TGFB3 13, LRP2 11, VPS13D 5, SLC12A5 3, DMXL2 2 and PTPN11 2) are held by a gate, 74 by the validity gate and 2 by the mechanism gate (Appendix Table A8).
A gate that works on the gene cannot see which disease a variant is pathogenic for.
The only ClinGen curation of RASA1 in the file is a Disputed one for Noonan syndrome, where ClinVar's records name capillary malformation-arteriovenous malformation for 39 of its 40 variants, and the only curation of LRP2 is "no known disease relationship" with congenital heart disease, where ClinVar's records name Donnai-Barrow syndrome for 10 of its 11.
The ClinGen guidance on genes with several disorders (ClinGen, 2024) is the answer to that, and v38 follows it as far as the mode of inheritance: a gene is classified once for each mode of inheritance of its established diseases (Section 3.8), which does not tell two diseases of one mode apart, and neither RASA1 nor LRP2 has an established curation to classify against.

The largest single gene in the gap is an input that was incomplete.
CFTR holds 439 of the 5,184 variants, and its 558 null variants in the benchmark (stop-gained, frameshift and canonical splice) are all called VUS and none carries PVS1; all 23 four-star records are variants of CFTR, and v38 calls 22 of them VUS.
All 558 are null on CFTR's MANE Select transcript (ENST00000003084.11), which is the transcript that the pick chose for each of them, so the failure does not come from the choice of transcript.
ClinVar classifies 556 of the 558 as pathogenic-side (528 Pathogenic and 28 Likely pathogenic) and 2 as VUS, and Supplementary Table S7 lists every one of them with the call that the completed ClinGen file gives.
ClinGen's export of 11 October 2026 lists CFTR as Definitive for cystic fibrosis (autosomal recessive, classified on 1 June 2022), but the export of 30 April 2026 that the benchmark read does not contain it, and nor do 103 other curations that ClinGen dates before that day: 104 curations in 100 genes, 76 of them Definitive, Strong or Moderate, in 73 genes.
The completed file adds the 112 rows of the October export that are dated on or before 30 April and whose gene and disease pair the April file lacks: the 104 curations, seven rows that re-label, under another disease identifier, an assertion that the April file lists, and one row (GATM, Moderate) dated 04:00 UTC on 30 April, which is after the April file was saved (02:19 UTC, by its modification time) and changes none of the 163 calls on GATM's variants.
It does change the criteria of 6 of them under v38, which reads the curated tag: the established tag of the gene is AR alone in the April file and AD and AR in the completed one, so the gene stops being read as recessive-only, PM2's threshold for it falls from 7e-5 to 4e-5, and the 6 variants, seen at allele frequencies of 4.2e-5 to 6.1e-5, lose PM2 without changing a call.
Their dates are recent: 83 of the 104 (80%) were classified in 2025 or 2026 and 45 (43%) from January to April 2026, against 664 of the 3,546 curations that the file lists (19%) and 159 (4.5%), where an omission independent of the date would leave about 19 and 5, which points to an export that lagged behind ClinGen's curations; five are older than 2024, CFTR (2022) among them, which a lag does not explain, and we have not established the cause.
We found it only by asking why CFTR had no PVS1.
Putting all 112 back, in a run that differs from v38 in that file alone, gives PVS1 to 674 of the 5,184 variants (13%, in 23 genes: CFTR 439, ELP1 85, EVC 46 and SPG7 41 among them, which are 611 of the 674) and a pathogenic-side call to 653.
Of CFTR's 558 null variants, 556 then carry PVS1 and 542 are called pathogenic-side, and 12 of the 23 four-star records are called pathogenic-side where 1 was.
The run adds 852 correct pathogenic-side calls and raises recall of pathogenic-side variants from 62.5% to 63.4%, and it adds 19 ClinVar-VUS variants called pathogenic-side and 1 reference-discordant pathogenic call, a frameshift variant of ALPK1 that ClinVar calls Benign, at an allele frequency of 2.3e-4; reference-discordant benign calls do not change.

Orphanet's loss-of-function label reaches more of the gap than the completed file does, and neither closes it.
Licensing PVS1 for the 1,042 genes that Orphanet labels as a loss-of-function disease and that the default table does not name (the association type "disease-causing germline mutation(s) (loss of function)" of Orphadata product 6, version of 23 June 2026, where a gene that also has a gain-of-function association is taken as both; Orphanet, 2026) gives PVS1 to 1,506 of the variants and a pathogenic-side call to 1,421, adds 1,829 correct pathogenic-side calls and 31 ClinVar-VUS variants called pathogenic-side, and adds 5 reference-discordant pathogenic calls (MAP3K6, ALPK1 twice, VPS13C and LTBP4, two of them splice-acceptor variants with population-maximum frequencies of 0.8%).
With the complete ClinGen file as well, 1,714 variants gain PVS1 and 1,624 are called pathogenic-side, agreement in direction rises by 2,067 calls (0.3 points), recall of pathogenic-side variants to 64.7%, and reference-discordant pathogenic calls by 4, about 517 calls per added error.
Even so 3,470 of the 5,184 variants (67%) carry no PVS1, and closing the rest needs a source of loss-of-function mechanism that does not depend on constraint, such as the PVS1 applicability stated in VCEP specifications, which our extraction of the CSpec registry does not read; we have not changed the classifier here.
In effect the classifier implements the standards' condition that loss of function be a known mechanism as gene-disease validity or constraint, which are proxies: a Definitive curation says that a gene causes disease and not that it does so by loss of function, and the curated-mechanism route licenses nothing in this benchmark.
The other 3,952 truncating variants that are called VUS carry PVS1 at a lower strength (3,066 at Strong, 462 at Moderate, 416 at Supporting and 8 at full strength but held back by other evidence), and 1,434 of the 1,446 canonical splice variants that are called VUS carry no PVS1 either.

**Where the 105,256 benign-side variants called VUS come from.**
Of the ClinVar benign-side variants that are called VUS, 32,168 (30.6%) carry PM2 and nothing else (they are rare in gnomAD and no benign evidence fires), 21,073 (20.0%) carry BP4 with PM2, 15,943 (15.1%) a lone BP4, 13,589 (12.9%) no criterion at all, 5,671 (5.4%) a lone BS1, and the rest smaller combinations.
Under the 2015 table that the benign side keeps, one supporting benign criterion does not reach Likely benign, and Section 3.3 shows what the alternative would cost.
At least 53,241 of these variants (51%) carry PM2, which is rarity and so does not separate a benign variant from a pathogenic one; 15,943 carry a lone BP4 and 5,671 a lone BS1, and neither reaches Likely benign under the 2015 table; and 13,589 carry no criterion at all.

**The 3,608 and 4,804 variants where v38 commits and ClinVar's submitters see insufficient evidence are not obviously errors.**
The pathogenic-side calls are PP2 with PM2 and PP3 at Strong strength (1,198), PVS1 with PM2 (1,133), PM1 with PM2, PP2 and PP3 at Moderate strength (587), PM1 with PM2, PM5 and PP2 (448), and a lone PVS1 (158).
The benign-side calls are chiefly BP1 with BP4 at Supporting or Moderate strength, with or without PM2 (3,673 of the 4,804).
Each pathogenic-side call stands on 6 or more points of computable evidence, which counts correlated predictors (PP2, PP3 and PM1 for a missense variant) as if they were independent, as the point system does, and each benign-side call on two or more benign criteria or on one at Very Strong or Stand-alone strength (122 of the 4,804).
A ClinVar VUS can mean that submitters weighed the evidence and found it insufficient, or that they did not weigh it to a class.
The benchmark cannot say which is right, and Section 3.7 asks whether ClinVar later moved toward them.

### 3.6 The 48 opposite-direction calls that remain

The last review round (v26) held 64 opposite-direction calls, and every one of the 48 in v38 is among them, so none of the 48 is new since that table was sent (one run, v31, briefly had 65 opposite-direction calls; Table 6, Figure S3B).
The 16 that were repaired all became VUS and none became the other class, and the reviewer had supported ClinVar's call on 14 of them (of the other two, one was coded neither and one depends on the disorder).
Each of the 48 had already been in the table of an earlier round (17 since round 3, 2 since round 4 and 29 since round 5), so they are the calls that several rounds of changes left standing.

**Table 6. The 48 opposite-direction calls of v38, by the evidence behind the call and the ruling of the last review round.**
Rulings are our coding of the reviewer's text (Section 2.4), not categories the reviewer assigned, and Supplementary Table S8 lists the 48 calls with their criteria, gnomAD counts and rulings.
ClinVar here is ClinVar 2-star+ of the release of 4 April 2026 (Section 2.2).

<!-- table:residual:begin -->
| Evidence behind the call | Calls | ClinVar supported | fastVEP supported | Neither | Needs more evidence | Depends on the disorder |
|:---|---:|---:|---:|---:|---:|---:|
| PVS1 on a ClinVar benign or likely benign variant | 21 | 13 | 2 | 3 | 2 | 1 |
| BP4 and BP7 on a ClinVar pathogenic or likely pathogenic variant | 20 | 10 | 2 | 5 | 3 | 0 |
| PM1, PP2 and PP3 on a ClinVar benign or likely benign variant | 3 | 3 | 0 | 0 | 0 | 0 |
| BS1 on a ClinVar pathogenic or likely pathogenic variant | 2 | 1 | 0 | 0 | 1 | 0 |
| BP1 on a ClinVar pathogenic or likely pathogenic variant | 2 | 2 | 0 | 0 | 0 | 0 |
| All | 48 | 29 | 4 | 8 | 6 | 1 |
<!-- table:residual:end -->

For 29 of the 48 (60%) the reviewer supported ClinVar's call (13 of the 24 reference-discordant benign calls and 16 of the 24 reference-discordant pathogenic calls), for 4 (8%) fastVEP's, for 8 the reviewer's own call was VUS, for 6 more evidence was needed and for 1 the answer depends on the disorder.
So for 29 of the 48 the reviewer judged ClinVar's call the better, which is not the same as the classifier being wrong in each of them, and for the other 19 the reviewer did not side with ClinVar.
Of the 29 rulings for ClinVar, 10 follow from a class that the reviewer stated, 9 from the reviewer's statement that ClinVar is correct, and 10 from evidence that the reviewer gave for ClinVar's class without naming one (Appendix Table A4).
Five of the 48 (10%) are no longer in the 2-star+ set in the October release, where ClinVar records them as having conflicting classifications, against 1.3% of all April variants that left it.
The 48 fall into five groups, and for most of them the reviewer's remarks name what decides the call, which is mostly evidence that a single-sample VCF does not carry.

**PVS1 on a variant ClinVar calls benign or likely benign (21).**
Ten of these are a lone PVS1 and 11 are PVS1 with PM2, and all are loss-of-function variants in a gene that the classifier takes to be sensitive to loss of function.
Five are seen in 611 to 4,892 gnomAD alleles (EDN3, FMR1, KDM2B, PHKB and PIEZO1), at allele frequencies of 5.6e-4 to 3.3e-3 across gnomAD and population-maximum frequencies of 0.37% to 1.5%.
BS1, which tests the filtering allele frequency against 5e-3, did not fire for any of the five, and their population-maximum frequencies are 15 to 59 times the median published VCEP threshold of 2.5e-4; the reviewer supported ClinVar on all five and, as we read the returned workbook of round 6, which is not in the repository, remarked that the variant is too common in gnomAD for KDM2B and PHKB, that it has more than 1,000 carriers for EDN3 and PIEZO1, and, for FMR1, gave its 210 hemizygotes and one homozygote and asked why the population frequencies were not applied.
Two of the five (KDM2B and PHKB) are also truncations in the first coding exon, where the reviewer noted that NMD can be escaped; neither the published PVS1 tree nor the classifier has a rule for it, and the classifier predicts NMD for every truncation outside the last exon.
The reviewer's other reasons were a haplotype in which the variant is always in cis with a second variant and the pair is not pathogenic (BRCA1); a splice-site variant in the 5' UTR of one RUNX1 isoform; a gene without an OMIM disease (RYK); a site that failed quality control in gnomAD (NOTCH2); a site covered in under half of gnomAD exomes (PRODH); and, for two variants (RYR1 and ARMC9), a call that is correct for the recessive disorder.

**BP4 and BP7 on a variant ClinVar calls pathogenic or likely pathogenic (20).**
These are the 20 reference-discordant benign calls of Section 3.4.
Seventeen are intronic, two synonymous and one in the splice region, with SpliceAI below 0.2 and PhyloP below 2.0, which is what BP7 asks.
The reviewer supported ClinVar for 10, and for several noted genetic or functional evidence of an effect on splicing (POLR3A c.1909+22G>A, at an allele frequency of 2.5e-3 in gnomAD, is one: the reviewer cites affected individuals and some functional data).
For 9 the reviewer's remark is that ClinVar's classification appears to predate the current MANE transcripts, including the two TH variants on which the reviewer supported fastVEP.
One of the 20 (DST) is a stop-gained variant on the MANE Plus Clinical transcript and a deep-intronic variant on the MANE Select transcript that the pick chose, and the reviewer asked that Plus Clinical take precedence where both exist, which is a question of transcript policy that we have not changed.

**Missense variants: PM1, PP2 and PP3 on a variant ClinVar calls benign or likely benign (3), BS1 on pathogenic missense variants (2) and BP1 on pathogenic missense variants (2).**
The first three (CHD7, FBN1 and PTCH1) reach 6 points through PM1, PM2, PP2 and PP3 at Moderate strength, and the reviewer supported ClinVar on all three.
The BS1 pair (PADI3 and SLC12A3) are seen in 328 and 338 gnomAD alleles, with population-maximum frequencies of 0.6% and 0.8%; the reviewer noted a homozygote in gnomAD and a disorder (Gitelman syndrome) for SLC12A3, supported ClinVar there, and asked for a deeper look at PADI3.
The BP1 pair (IQSEC1 and PRR12) are in genes with two and one pathogenic or likely pathogenic missense records in the ClinVar protein index, the variant's own among them, so BP1's veto does not apply (a count of 1 would withhold BP1 for both, for PRR12 only through its own record and for IQSEC1 also through p.Thr343Met, Section 3.3), and the reviewer supported ClinVar on both.

No change of one threshold would remove these calls as a class, except that a stricter count of BP1's veto acts on the BP1 pair (Section 3.3).
A disease-specific frequency threshold could act on about seven of them (the five common variants of PVS1 and the two of BS1), a transcript policy on a handful, and for the rest of the others the reviewer's remarks point to phase, case-level or functional evidence.
We have not put them on the curated exception list, which holds 15 alleles, because a list long enough to hide the residual would no longer be a statement about the classifier.

### 3.7 A later ClinVar release, and the variants that were not in the April set

ClinVar's release of 4 October 2026 holds 781,135 records at two stars or better: 665,057 of the April set, which are still at two stars or better, and 116,078 that are new.
A further 8,603 of the April variants have left the set, among them five of the 48 residual calls (Section 3.6).
The gnomAD extract the classifier reads covers every April variant, and 107,528 of the new ones; the rule of Section 2.8 drops the other 8,550, which is 7.4% of the new records and 1.1% of the release.
Most of the new variants were already in ClinVar in April, at a lower review status: 105,920 of the 107,528 (98.5%) have a record in the `variant_summary` of 29 April from which the protein index was built, and for 98,623 (92% of all the new variants) it is a single-submitter record.
They are variants that reached two stars since, and 9,234 of them already had a pathogenic or likely pathogenic record in April.
Table 7 gives v38, unchanged, on the April set and on each stratum of the covered release (Supplementary Table S9; Appendix Table A9 holds the rest of this section).

**Table 7. v38 on the April benchmark and on the later ClinVar release, by stratum.**
A variant is new if it is not in the April set, unchanged if it has the same class in both releases and has moved within or across a direction if its class changed.
Opposite-direction calls are given as the total, with reference-discordant benign calls and reference-discordant pathogenic calls in brackets.
ClinVar P/LP, LB/B and VUS are classes of ClinVar 2-star+ records of the release of 4 April 2026 for the April set and of 4 October 2026 for the later release.

<!-- table:temporal:begin -->
| Set | Variants | Same-direction agreement | ClinVar P/LP called P/LP | ClinVar LB/B called LB/B | ClinVar VUS called VUS | Opposite-direction calls (reference-discordant benign and reference-discordant pathogenic) |
|:---|---:|---:|---:|---:|---:|---:|
| April release, the benchmark set | 673,660 | 77.8% | 62.5% | 62.9% | 97.2% | 48 (24 and 24) |
| October release, covered by the gnomAD extract | 772,585 | 76.8% | 62.7% | 61.7% | 97.2% | 63 (26 and 37) |
| &nbsp;&nbsp;unchanged since April | 662,668 | 77.9% | 63.1% | 63.1% | 97.2% | 43 (20 and 23) |
| &nbsp;&nbsp;class moved within its direction | 2,232 | 54.9% | 37.4% | 60.1% | - | 0 (0 and 0) |
| &nbsp;&nbsp;class moved across a direction | 157 | 19.1% | 16.2% | 7.1% | 77.3% | 2 (0 and 2) |
| &nbsp;&nbsp;new: not in the April set | 107,528 | 70.3% | 61.0% | 55.6% | 97.3% | 18 (6 and 12) |
<!-- table:temporal:end -->

**On variants that were never in the benchmark, agreement in direction is 70.3% against 77.8%, mostly on the benign side, and contradictions are 2.3 times as frequent.**
Across the 772,585 covered records, agreement in direction is 76.8% against 77.8% in April, with 63 opposite-direction calls (26 reference-discordant benign calls and 37 reference-discordant pathogenic calls).
On the 662,668 variants whose class did not change, the figures are those of April (77.9%, 43 opposite), as they should be.
On the 107,528 new variants agreement in direction is 70.3%, and recall of pathogenic-side variants is 61.0% (April 62.5%), of VUS 97.3% (97.2%) and of benign-side variants 55.6% (62.9%).
The new variants are a different mix of ClinVar classes (9% pathogenic-side, 34% VUS and 57% benign-side, against 14%, 44% and 42% in April), and at the April recall of each direction their agreement would be 74.5%, so taking the mix first 3.3 of the 7.5 points are the mix and 4.2 are lower recall, nearly all of it benign-side (the other order gives 4.3 and 3.2, splitting the same 7.5 points).
The fall on the benign side is mostly a change in composition (Figure 5C).
New Benign variants are much rarer than the old ones (median population-maximum allele frequency 0.28% against 3.6%, and 31% above 0.5% against 70%), and Benign is the class that BA1 can reach by frequency alone: recall of Benign variants falls from 72.2% to 53.5%, while recall of Likely benign variants, which frequency rarely decides, rises from 51.5% to 56.0%.
At the April recall of each of nine bands of population-maximum allele frequency, the new Benign variants would have a recall of 57.6%, so the mix of bands explains 78% of the fall (14.6 of 18.7 points), and for Likely benign variants it explains most of the rise (expected 55.3%).
We made the split into nine bands after a split into six, with one band from 1e-2 up and one from 1e-3 to 1e-2, which explained 63%, had shown where the rest lay, so the 78% describes the fall and does not test it.
What remains is spread over the bands: the two at 2% or more carry 1.8 of the 4.1 points that the mix does not explain (43%), the two from 5e-3 to 2e-2 carry 1.2, and the band from 1e-5 to 1e-4 carries 0.8.
The two at 2% or more have the widest gaps between the April recall and the new recall (23 and 35 points), and there the new Benign variants are called benign-side for 58% (326 of 564) against 91% of 83,527 April ones, and 157 of the 238 that are not carry BP4 alone, so neither BA1 nor BS1 fired on them.
Annotating those 238 again with the details of each criterion shows why.
For 189 (79%) BA1 and BS1 were withheld by design: 124 sit in a low-complexity region that gnomAD flags, 56 at a site whose gnomAD record fails its own filter, 5 in a segmental duplication and 4 in a gene of the homologous-gene list, each counted once, by the first of these that applies (Section 3.4).
For 16 (7%) the criteria were evaluated and the filtering allele frequency is below the threshold, which is how one small group can inflate a population maximum.
Another 32 met BS1 and stay at VUS, 24 with BS1 alone, which the 2015 table on the benign side does not call, and 8 with pathogenic-side evidence beside it.
One more has an allele number below the floor.
The 189 are 34% of the 564 new Benign variants of that frequency, against 4.8% (4,035 of 83,527) of the April ones, so the new variants are far more often at sites where the classifier does not trust the gnomAD frequency, and we have not established why.
Recall of ClinVar Pathogenic variants is 64.7% against 65.0%, and of Likely pathogenic variants 51.4% against 48.3%.

**The rate of opposite-direction calls is higher on new variants, and the kinds resemble the old ones.**
Among the new variants it is 0.017% (18 of 107,528; 95% interval 0.009 to 0.025% by gene-cluster bootstrap), 2.3 times the 0.0071% of the April set (0.005 to 0.010%), and the interval of the ratio, from a bootstrap that resamples genes jointly across the two sets, is 1.2 to 4.1.
The mix of ClinVar class and consequence does not explain the excess, taking the April rates as known: at the April rate of its own stratum, the new set would have had 5.7 opposite-direction calls by ClinVar direction and consequence group, 5.9 by ClinVar class and consequence group, 6.7 by class alone and 7.7 by direction alone, against the 18 observed.
The April figure was reached after every contradiction in it had been reviewed and many repaired, and the new variants had had no such attention, which would raise their rate whatever else differs.
We cannot separate that from other differences between the sets, such as how the new variants reached two stars (nine in ten were single-submitter records in April).
The 18 are of three kinds that we had named on the April residual: 9 are PVS1 on a ClinVar benign-side variant (ALMS1, CACNA1G, CCNQ, MAP3K14, ORAI1, PRKDC, SMARCA2, TRPV4, UNC93B1), 6 are BP4 and BP7 on a ClinVar pathogenic-side variant (ABCA4, HEXA, IFT122, NUP107, PKD1, SLC7A9), and 3 are missense variants called pathogenic-side on PP2 and PP3 (COL1A1 twice and KCNH1).
The kinds are coarse and were defined after the April residual was seen, and 18 events give little power to find a new one, so this shows that the new contradictions resemble the old ones and not that no other kind exists.

**The changes proposed in Section 4.5 do not fail on the new variants where they can be tried.**
BP7 as the splicing recommendation describes it adds 4,859 calls in direction on the 107,528 new variants (4.5% of them, against 2.3% of the April set) and no reference-discordant benign call, where on the April set it added 15,552 and 2.
The completed ClinGen file adds 50 calls in direction and no opposite-direction call, where on the April set it added 848 and 1 reference-discordant pathogenic call; Orphanet's label adds 272 for 2 reference-discordant pathogenic calls (1,795 and 5 on the April set), and the two together 311 and 2 (2,067 and 4); the rates for Orphanet's label and for both are 136 and 156 calls per added reference-discordant pathogenic call on the new variants against 359 and 517 on the April set, and on the new variants rest on 2 events.
BP1's veto without the variant's own record changes no call, and the test is weak: the veto changes for a variant only where the gene's count is exactly 3, which is 2 of the 94 new ClinVar-pathogenic missense variants in genes where BP1 can fire (NUP85 and RPS29, whose own records are in the index), and neither has BP4 beside BP1 to make a call from.
The measurement build discounts one record whenever ClinVar's allele-level record of 27 June is pathogenic, which for up to 17 of the 94 (those whose own record the April index does not hold as pathogenic) discounts a record that is not in their count, and the run changes no call for them either.
BP1's veto at a count of 1 gives 124 more calls in direction (157 fewer ClinVar-VUS variants called benign-side and 33 fewer correct benign-side calls) and changes no opposite-direction call, because none of the 6 reference-discordant benign calls among the new variants rests on BP1; on the April set it gave 1,735 calls and removed 2 reference-discordant benign calls.
Removing the allele-number floor changes no call, which tests it on at most 3 variants: 33 new variants have an allele number below 2,000, and 3 of them have a population maximum at or above the BS1 threshold of 5e-3.
These are seven runs on one set of variants, and the new variants lie near April variants (Section 2.8) and were mostly in ClinVar in April, so they are evidence that the changes do not fail on variants that no tuning scored and not that they would be free.

**ClinVar moved few of the calls the classifier made against a VUS, and the 157 variants that changed direction were mostly ones the classifier had left as VUS.**
ClinVar changed the direction of 157 of the 665,057 April variants in six months, 36 from VUS to a pathogenic-side class, 98 from VUS to a benign-side class, 22 into VUS (14 from pathogenic-side and 8 from benign-side) and one from benign-side to pathogenic-side.
The classifier called VUS for 119 of the 134 that left VUS, and agreed with ClinVar's new class for 30 of the 157 and with its earlier class for 124.
For the 22 moved into VUS it agreed with the new class for 17.
Of the 8,412 April variants where the classifier committed against a ClinVar VUS, 13 later moved in the direction of its call and 2 moved against it, so a test of whether ClinVar later moved toward those calls cannot be made from six months of reclassification.
The 43 residual calls that are still in the set are unchanged in class.

**What the later release does not test.**
The new variants lie in 4,124 genes, all but one of which are in the April set, so this is a hold-out on variants and not on genes, and every gene-level setting (the exception list, the homologous-gene list, the ClinGen gene tags) was set with the April genes in view.
The ClinVar-derived indices the classifier reads were not rebuilt, so the protein index lacks every record added since April, and holds the April record of each new variant that was then a classified missense variant, which for BP1's veto is the variant's own record (for 77 of the 94 new ClinVar-pathogenic missense variants in genes where BP1 can fire, and for the other 17 the April index does not hold it as pathogenic), as it is for 788 of the 790 April variants of the same kind.
The records of other variants that were added or reclassified since lower PS1 and PM5 for a new variant, lower or raise PM1 (fewer pathogenic neighbours lower it and fewer benign ones raise it), and make the veto of BP1 less likely to apply, which could add reference-discordant benign calls.
A new variant that already has a pathogenic record in the allele-level annotation of 27 June has that record left out of the index counts of PS1 and PM1 although the April index does not hold it, which lowers both further.
None of the 18 opposite-direction calls on the new variants rests on BP1, PS1 or PM5 and one rests on PM1, so the missing records can have produced at most one of them.
And six months is short: ClinVar reclassifies slowly, so the later release adds mostly second submissions to variants that ClinVar already held, and little new knowledge about old variants.

### 3.8 What v38 changed: reading the mode of inheritance that ClinGen curates

Sections 3.1 to 3.7 are about v38, the latest classifier, and this section says what it changed from v37, the version in which the review rounds ended, and why.
Up to v37 the criteria that branch on a gene's inheritance took it from the words of the disease names in the ClinGen file, and the mode of inheritance that ClinGen curates for each disease, which is in the same string, was not read (Section 2.10).
v37 left the field unread: limitation 20 of `docs/ACMG.md` recorded it with the counts below, and the record of v37 (RUN_VERSIONS.md) measured a build that made BS2 read the tag (`v37_adtag`), which gained 2,812 calls in direction and moved six ClinVar pathogenic-side variants toward benign: TTR p.Val142Ile and VWF p.Arg854Gln to Benign, and loss-of-function alleles of NEK1, DDX41, SRCAP and PROKR2 to VUS.
The sweeps of v37 did not show what the field would change: the two that it should have moved, the PM2 threshold for genes read as recessive-only and the BS2 carrier threshold for genes read as dominant-only, were flat, and they reached 67 and 50 genes (Section 3.3).
v38 reads the field.
Every run of Sections 2.6 to 2.8 was made with both versions, the 96 runs of v38 and the same 96 of v37, by the same harness on the same inputs, so that a pair of runs differs in the build that made it, the classifier and the annotation layer that it was built with (below), and in nothing else (`analysis/acmg_benchmark/results/v38_family_vs_v37.md` sets them side by side).
The tables of this section set the change against the commit before it, which has the classification code of v37 (the two differ in comments only) on a newer annotation layer.
Its run differs from v37's in the criteria of 25 indels and in the call of 4 of them, ClinVar VUS variants that v37 called Likely pathogenic, so it has 425,873 exact calls and 522,890 in direction where v37 has 425,869 and 522,886.
Rates are over the 673,660 records, as in the rest of the paper, and 673,278 of them have a call.
The runs and their controls are listed in Appendix Table A10.

**The field reads 2,056 genes as dominant-only or recessive-only where the names read 117.**
Of the 2,924 genes of the ClinGen file, the names gave 50 dominant-only, 67 recessive-only, 12 both and 2,795 neither, and the field gives 746, 1,310, 169 and 699 (Appendix Table A11).
Among the 2,840 genes with a benchmark variant, the names read as dominant 39 of the 740 whose established curations are all autosomal dominant, and as recessive 60 of the 1,285 whose established curations are all autosomal recessive, and of the 114 genes with both kinds of curation 15 read as dominant and 13 as recessive.
Of the other 701 genes, 454 have no established curation, 162 are X-linked only, 32 are mitochondrial only and 53 carry a semidominant curation, which is read as both dominant and recessive.
In the file 129 genes, 127 of them with a benchmark variant, have established curations of more than one mode of inheritance.

**On the benchmark v38 gives 1,494 more calls in direction than v37 and one opposite-direction call fewer, and the gain is BS2's carrier arm reaching the dominant genes.**
v38 agrees with ClinVar's class for 426,022 of the 673,660 records and with its direction for 524,380 (77.84%), where v37 agrees for 425,869 and 522,886 (77.62%): 153 and 1,494 more calls, which is 0.22 points of the records.
Of the 1,494, 4 are the commit before the change and 1,490 are the change itself, so that against the commit before it v38 gives 149 more exact calls and 1,490 more calls in direction (Table 8).
The gain is on the benign side: 1,682 more ClinVar benign-side variants are called benign-side (recall 62.28% to 62.87%) and 4 fewer pathogenic-side variants are called pathogenic-side (62.48% to 62.47%).
Opposite-direction calls go from 49 to 48 (24 reference-discordant benign calls and 25 reference-discordant pathogenic calls under v37, 24 and 24 under v38): one call leaves the set and none enters it.
The call that leaves is RAD21 8:116866762 C>T, a splice-acceptor variant that ClinVar classifies as Likely benign at two stars: v37 called it Likely pathogenic on PVS1 alone, and v38 calls it VUS on PVS1 and BS2, because the curated tag of the gene is AD and gnomAD holds the variant in 2,630 of 1,381,576 alleles with 3 homozygotes, so BS2's carrier arm applies, and Strong benign evidence that does not reach a call holds a pathogenic call at VUS (Section 2.1).
In round 6 the reviewer supported ClinVar's call on this variant.

**Table 8. The classifier before and after reading ClinGen's curated mode of inheritance, and the alternatives that were measured.**
Each row is a run of the harness on the April set (Section 2.10), and the first row is the commit before the change.
Reference-discordant benign calls are ClinVar pathogenic-side variants called benign-side, reference-discordant pathogenic calls are ClinVar benign-side variants called pathogenic-side, and rates are over 673,660 records.
The recall of ClinVar's benign-side variants is over the 284,168 that have a call, because the binaries that made these runs write no record for the 382 reference-agreement records (ALT is a dot) that Sections 3.1 to 3.7 count as no-calls, which makes it 0.1 points higher than the recall over the 284,550 of those sections.
Appendix Table A10 holds every run.

<!-- table:inheritance:begin -->
| Setting | ClinVar P/LP called P/LP (recall) | ClinVar B/LB called LB/B (recall) | Reference-discordant benign calls | Reference-discordant pathogenic calls | Same-direction agreement | Exact agreement |
|:---|---:|---:|---:|---:|---:|---:|
| Before: inheritance read from the disease names | 58,610 (62.5%) | 177,206 (62.4%) | 24 | 25 | 522,890 (77.6%) | 425,873 (63.2%) |
| **After (default)**: curated inheritance; BS2's carriers not counted beside BS1 or BA1 | 58,606 (62.5%) | 178,888 (63.0%) | 24 | 24 | 524,380 (77.8%) | 426,022 (63.2%) |
| Curated inheritance; BS2's carriers counted beside BS1 and BA1 | 58,606 (62.5%) | 180,184 (63.4%) | 26 | 24 | 525,672 (78.0%) | 429,637 (63.8%) |
| Default, and no BS2 carrier arm for a susceptibility or predisposition disorder | 58,608 (62.5%) | 178,823 (62.9%) | 24 | 24 | 524,327 (77.8%) | 426,009 (63.2%) |
| Curated inheritance for PM2, PM3 and BP2; BS2 as before | 58,610 (62.5%) | 177,206 (62.4%) | 24 | 25 | 522,887 (77.6%) | 425,889 (63.2%) |
| Default, with PP2, PM1 and PP3 capped together at 4 points | 58,286 (62.1%) | 178,888 (63.0%) | 24 | 21 | 524,647 (77.9%) | 426,534 (63.3%) |
<!-- table:inheritance:end -->

v38 gives a benign-side call to 1,809 ClinVar benign-side variants that were VUS and takes it from 127, so that recall of ClinVar's benign-side variants rises from 62.4% to 63.0% (over the 284,168 that have a call; 62.3% to 62.9% over the 284,550 of Sections 3.1 to 3.7), and it gives a benign-side call to 188 more ClinVar-VUS variants (4,804 against 4,616) and changes the number called pathogenic-side by none (3,608).
Reference-discordant benign calls stay at 24, and reference-discordant pathogenic calls fall from 25 to 24.
Four ClinVar pathogenic-side variants lose their pathogenic-side call, which is taken up below.
PM2 is met by 9,335 more variants (441,577 against 432,242), because the 1,310 genes now read as recessive-only use PM2's higher threshold (7e-5 against 4e-5), but that moves little: the run that reads the field for PM2, PM3 and BP2 and leaves BS2 as it was adds 16 exact calls and takes 3 from agreement in direction (Table 8).
The gain is BS2's.
Its carrier arm now reaches the 746 genes read as dominant-only and not 50, and BS2 is met by 42,900 variants against 45,842 before, and by 59,003 when it is also counted beside BS1 and BA1 (Appendix Table A10).
It is met by fewer variants than before because the rule that sets the carrier count aside withholds 16,103 carrier firings, the difference between the 59,003 and the 42,900, and because the genes that are newly dominant-only no longer take the arm for homozygotes.

**Counting BS2 beside BS1 would have called two ClinVar-Pathogenic variants Benign, and v38 does not.**
The reading as it first ran, with BS2's carriers counted beside BS1 and BA1, gives 3,615 more exact calls and 1,292 more calls in direction than v38 (Table 8).
It also calls two ClinVar-Pathogenic variants Benign that are VUS in v38: TTR p.Val142Ile (722 gnomAD alleles, 3 homozygotes; PM5, BS1 and BS2) and VWF p.Arg854Gln (7,432 alleles, 23 homozygotes; PM1, PM5, PP2, BS1 and BS2).
In genes whose established curations are all autosomal dominant the carrier arm is met on 7,234 of the 7,328 variants that meet BS1 and on all 8,869 that meet BA1, which is 16,103 firings that add a second Strong criterion for one gnomAD record, which we read as one observation counted twice: the carrier arm reads the global counts and BS1 and BA1 the filtering allele frequency of the population with the highest one, so the two are summaries of the same record, and that reading is ours and not a statement of the ClinGen SVI.
A stand-alone BA1 is Benign already, so the doubling changes a call only beside BS1: it changes 6,494 calls, 5,192 from Likely benign and 1,302 from VUS to Benign, of which 6,478 are ClinVar benign-side variants, 14 are ClinVar VUS variants and 2 are TTR and VWF (Appendix Table A12).
The guard is not free.
Against counting both it gives up 3,615 exact calls and 1,292 calls in direction on the April set, almost all on ClinVar benign-side variants that a frequency makes Benign when it is counted twice.
On the 107,525 new variants of the October release it gives up 150 exact calls and 46 in direction and avoids no contradiction, because neither run has one among them.
So it rests on the double counting and on two variants, and not on a gain that the benchmark shows on variants that were not seen when it was written.
We kept it, because a call that counts one observation twice reads as more certain than the evidence is.

**The arm for homozygotes and hemizygotes is not set aside, and the same doubling remains in it for 1,301 variants.**
A homozygote is a genotype that an allele frequency does not give, but the hemizygotes of the non-pseudoautosomal X are part of the allele count that BS1 and BA1 read, so for them BS2 repeats the frequency.
BS2 is met beside BS1 or BA1 through that arm on 36,818 variants of the benchmark, 1,598 of them on the X chromosome.
BA1 alone is Benign, so the doubling changes a call only beside BS1 without BA1, which is 1,301 variants (443 on the X chromosome), all called Benign, of which 1,298 are ClinVar benign-side variants and 3 are ClinVar VUS variants.
No ClinVar pathogenic-side variant is among them, so the benchmark shows no harm from it, and setting the arm aside would only lower the call of variants that ClinVar also calls benign-side; we left it as it was and list it among the limitations (Section 4.4).

**Four ClinVar-Pathogenic variants lose their call to BS2's carrier arm.**
In four loss-of-function variants, NEK1 p.Ser1036Ter, DDX41 p.Asp140GlyfsTer2, SRCAP p.Pro1879ThrfsTer21 and PROKR2 p.His20MetfsTer24, the carrier arm is met beside PVS1, which takes each from Likely pathogenic to VUS (Appendix Table A12).
Each is seen in 295 to 324 gnomAD alleles, and none as a homozygote.
BS2 asks for full penetrance expected at an early age, and the classifier has no field for penetrance or age at onset, so it applies the arm to every gene read as dominant.
Withholding the arm for a disorder whose label says susceptibility or predisposition restores two of the four and gives up 65 correct benign-side calls (53 calls in direction on balance), and it reads the name of a disease again, which is what v38 stops doing, so we did not adopt it and the four are a limitation (Section 4.4).
The same four variants lose the call on the October release.

**v38 reads the field although v37 declined to, and the reason is not the aggregate.**
v37 left the tag unread because reading it for BS2 moved six ClinVar pathogenic-side variants toward benign, and the guards that removed them were a pair chosen after the six were found (`docs/ACMG.md`, limitation 20 as of v37).
The unread field makes the classifier answer a question about the name of a disease where the source curates the answer, for 2,056 genes where the names answered for 117, and that is a defect in what is read and not a threshold to set.
What the measurement adds is that two of the six have a cause that is not the two variants: the carrier count is the frequency that BS1 and BA1 have already counted, so the rule that sets it aside is stated by what it counts and not by a list of variants, although it was written after TTR and VWF were seen (below), and it changes 6,494 calls.
The other four are not removed.
They move from Likely pathogenic to VUS and not to a benign-side call, and we accept them as the price of an arm that applies to every gene read as dominant because the classifier has no penetrance to read, which is what v37 said a wider reading would need.
The label guard that restores two of the four gives up about 32 correct benign-side calls for each variant restored, which is below the hundred that Section 2.5 takes for a knee, so by the exchange rate alone we would adopt it.
We did not, because it is a filter on the name of a disease that was written to repair two variants seen in the results, which Section 2.5 rules out ("never move a threshold to repair one variant"), and because it reads the name of a disease, which v38 stops doing.
That is a decision and not a measurement, and the run that prices the other choice is `v38_alt_label_guard`.

**A combined cap on PP2, PM1 and PP3 trades 320 correct pathogenic-side calls for 587 VUS kept at VUS and 3 contradictions avoided.**
PP2, PM1 and PP3 on one missense variant are three criteria that the point system counts as independent, except that PM1 is set aside beside a Strong PP3 (Section 2.1).
All three are met for 1,591 variants, 913 of them called Likely pathogenic.
Capping the three together at 4 points changes 910 calls, all from Likely pathogenic to VUS (Table 8, Appendix Table A12): 320 of them are ClinVar pathogenic-side variants, 587 are ClinVar VUS variants and 3 are ClinVar benign-side variants that were called pathogenic-side.
It raises agreement in direction by 267 calls and exact agreement by 512, lowers recall of ClinVar's pathogenic-side variants from 62.5% to 62.1%, and removes 3 of the 24 reference-discordant pathogenic calls, at 107 correct pathogenic-side calls given up for each.
On the October release the same cap changes 1,018 calls and gives up 354 correct pathogenic-side calls for the 4 contradictions it avoids, and on the new variants it changes 123 calls and gives up 37 for the one it avoids.
We did not adopt it, because it changes the combining rule and does not correct a defect: the point scale that the classifier follows on the pathogenic side counts the criteria as independent, the benchmark does not separate the two readings (the 587 ClinVar VUS variants that it keeps at VUS are not shown to be errors, Section 3.5), and an exchange of 107 calls per contradiction is at the level that Section 2.5 takes for a knee, so it does not decide.

**A gene with more than one mode of inheritance is classified once for each, and the calls differ for 0.55% of the variants of such genes.**
Of the 673,278 variants with a call, 87,599 (13.0%) are in 128 genes whose established curations carry more than one mode of inheritance, by the symbol of the annotation, and each is classified against each mode besides the gene as a whole (Appendix Table A13).
The modes meet different criteria for 4,437 of the variants, always BS2 (2,976) or PM2 (1,461), the two criteria that branch on inheritance and can be met in the benchmark.
They give different classes for 482 variants (0.55%) and different directions for 435 (0.50%).
Among the ClinVar classes of such genes the modes give different classes for 15 of 20,754 pathogenic-side variants, 88 of 34,600 VUS variants and 379 of 32,245 benign-side variants.
In no gene do they differ for more than 30 variants (PKD1 30, BRCA2 24 and RYR1 19).
The benchmark cannot say which mode's call is right, because a ClinVar record of it has no disease-specific class: the calls of the modes are not scored, and what the benchmark shows is that they seldom disagree.

**Against v37 the sweeps of Section 3.3 keep their knees, and the two thresholds that reach genes by inheritance change.**
Of the 12 sweeps, 4 differ from v37 in the sign of a step, the position of a knee or the order of their settings, and none differs in the knee: the rule of Section 2.5 finds the knee of the BS1 threshold at 1e-3 and of BP7's distance at 300 nucleotides under both versions (Appendix Table A14 holds the sweeps of the four thresholds that read inheritance).
The BS2 carrier threshold and its floor stop being flat: with the threshold off, reference-discordant benign calls are 47 under v38 and 24 under v37, and a floor of 20 alleles gives 33 and 24 (Section 3.3).
The PM2 sweep under the 2015 table loses most of its range, from a gain of 16,785 variants to one of 5,700.
The PM2 threshold for genes read as recessive-only reaches 1,310 genes where it reached 67, and its sweep stays flat: recall moves across 62.43% to 62.48% and agreement in direction by 27 calls (10 under v37), neither kind of reference-discordant call changes, and its settings change order by a few calls, which is not a change of shape in the sense that matters.
The BS2 threshold for individuals with no functional copy keeps the shape it had: below the default 1e-3 each step buys 401, 272 and 487 correct benign-side calls per added reference-discordant benign call, against 463, 291 and 533 under v37.

**Against v37 the ablations of Section 3.4 keep their effects, except for the removal of the ClinGen file.**
Of the 245 effects of the single switches, the removals of one source, the hook runs and the runs that fill the PVS1 gap, on seven measures, 11 changed sign, appeared or vanished, or changed by more than 20%, and 5 of the 11 are of fewer than 10 calls, which a single variant makes.
The one that is large is the removal of the ClinGen file: it costs 1,845 correct benign-side calls under v38 where it cost 163 under v37.
The ClinVar benign-side variants that the removal turns into VUS are 1,858 under v38 against 177 under v37, PM2 at Supporting strength fires on 10,011 fewer variants under v38 and on 676 fewer under v37, and BS2 fires on 1,984 more variants under v38 and on 958 fewer under v37.
PM2 and BS2 are the two criteria that branch on the mode of inheritance, and the removal takes away the field that v38 reads for them and the words that v37 read; which of the two accounts for the loss of benign-side calls we have not traced, because no per-variant table was made for this run.
The completed ClinGen file changes the correct benign-side calls by 20 under v38 where it changed them by 1 under v37, and PM2 at Supporting strength by 172 variants and BS2 by 47, where under v37 it changed them by 0 and 8; that too we have not traced.

**On the new variants of the October release the effect is 258 more calls in direction, and the guard costs the same.**
On the 772,585 covered records v38 gives 261 more exact calls and 1,762 more calls in direction than the commit before it, with 26 reference-discordant benign calls before and after and 38 reference-discordant pathogenic calls before and 37 after (Appendix Table A15), 1,504 of the 1,762 on variants of the April set.
On the 107,525 new variants with a call it gives 101 more exact calls and 258 more in direction, with 6 reference-discordant benign calls and 12 reference-discordant pathogenic calls before and after.
Against v37, v38 agrees in direction for 1,767 more of the 772,585 covered records (76.6% to 76.8%), with 63 opposite-direction calls where v37 has 64 (the RAD21 variant again), and on the 107,528 new variants agreement in direction goes from 70.1% to 70.3% with the same 18 opposite-direction calls, which are 2.3 times the April rate under both versions.
Counting BS2 beside BS1 and BA1 adds 4,034 exact calls and 3,096 in direction over the commit before the change, and two reference-discordant benign calls, TTR p.Val142Ile and VWF p.Arg854Gln again, both of them old variants.

**gnomAD's genomes were not read by the frequency criteria, and reading them through the joint release gives more calls in direction than any other change of this section, for three more opposite-direction calls.**
The classifier read the exomes only, in the benchmark (Section 4.4).
The joint v4.1 release adds the genomes' counts to the exomes', and it declares no `lcr`, `segdup` or `non_par` flag (Section 2.10).
On every chromosome (673,278 variants) a database built from it holds a record for every variant that the exomes hold, with an allele count, an allele number and a count of homozygotes that are never lower (537,038 variants), and one for 17,382 more (2.6%) that the exomes lack (Appendix Table A16).
It changes 13,302 calls (2.0%), 12,983 of them on ClinVar benign-side variants, 209 on VUS and 110 on pathogenic-side variants, and gives 10,190 more exact calls (1.51 points) and 9,321 more in direction (1.38 points) than the exomes, for one more reference-discordant benign call (24 to 25) and two more reference-discordant pathogenic calls (24 to 26) (Table 9).
ClinVar benign-side variants called benign-side rise from 178,888 to 188,310, and pathogenic-side variants called pathogenic-side from 58,606 to 58,632.
Of the 13,302 calls that move, 4,485 are on variants that the exomes lack, among them 4,242 of the 6,476 intronic variants whose call moves, and the exomes' records of the 2,171 UTR variants whose call moves rest on a median of 394 alleles, against 152,416 in the joint release (Appendix Table A16).
BA1 is gained on 10,115 of the moved calls and lost on 339.
Each of the 22 autosomes and the X gains calls, with 1.3% to 2.7% of its variants changing, and nothing moves on the Y or the mitochondrial contigs, which have 4 and 550 variants.

**Table 9. gnomAD's joint release against its exomes, with the region flags separated out.**
The joint release declares no `lcr` or `segdup` flag, so its run differs from the exomes in two ways, the genomes' counts and the missing flags, and the run of the exomes with the flags withholding nothing shows the second alone.
Changes are from the run of the exomes.

<!-- table:joint:begin -->
| Source of the counts | Variants | ClinVar P/LP called P/LP | ClinVar B/LB called LB/B | Reference-discordant benign calls | Reference-discordant pathogenic calls | Same-direction agreement (change) | Exact agreement (change) |
|:---|---:|---:|---:|---:|---:|---:|---:|
| gnomAD exomes | 673,278 | 58,606 | 178,888 | 24 | 24 | 524,380 | 426,022 |
| gnomAD exomes, region flags not withholding | 673,278 | 58,634 | 182,437 | 24 | 25 | 527,853 (+3,473) | 430,097 (+4,075) |
| gnomAD joint exomes and genomes | 673,278 | 58,632 | 188,310 | 25 | 26 | 533,701 (+9,321) | 436,212 (+10,190) |
<!-- table:joint:end -->

Part of the gain is the flags that the joint release lacks and not the genomes: the exomes with the flags withholding nothing give 4,075 more exact calls and 3,473 more in direction, with one more reference-discordant pathogenic call and no more reference-discordant benign call, so the genomes add 6,115 and 5,848, with one more of each.
A pilot of the same three runs on chromosomes 21 and 22 alone, made before the other chromosomes had been extracted, gives the same call and the same criteria for each of its 23,538 variants as the three runs of every chromosome do (Appendix Table A16).

Five variants change their standing as an opposite-direction call, four added and one removed (Appendix Table A16).
Four of the five are at the threshold of BS1, a filtering allele frequency of 0.5%, and the fifth belongs to the missing flags.
KMT2C p.Cys391Ter (ClinVar Likely benign) and RAI1 p.Gln280HisfsTer83 (Benign) are null variants at sites that failed gnomAD's variant-quality filter, which withholds the frequency criteria, and the site-level frequency veto (Section 2.7) holds a call that rests on PVS1 alone at VUS where the withheld frequency would have met BS1.
Their filtering allele frequencies are 0.522% and 0.537% in the exomes and 0.495% and 0.498% in the joint release, so the veto does not act there and both are called Likely pathogenic.
PIEZO1's splice-acceptor deletion (ClinVar Benign) is the same case in the other direction, with 0.198% in the exomes and 0.947% in the joint release, so the veto holds it at VUS where the exomes called it Likely pathogenic.
TXNL4A's deletion upstream of the gene (ClinVar Pathogenic) is covered by 364,572 alleles in the exomes, with a filtering allele frequency of 0.384%, and by 516,266 in the joint release, with 0.560%, so BS1 fires and the call is Likely benign.
NF1 p.Met895Ile (ClinVar Likely benign) is called Likely pathogenic because PM2, which the segmental-duplication flag of the exomes had withheld, fires.
Thirteen ClinVar pathogenic-side variants lose a pathogenic-side call and 39 gain one.
Ten of the 13 are missense variants whose record in the joint release carries a filter of gnomAD's quality control (AC0 for eight, AS_VQSR for one and both for one), which withholds PM2: five have a clean record in the exomes, with one to four alleles, and the filter is the genomes', which the database reports because it reports a filter that fired in either callset (Section 2.10), and five are not in the exomes at all.
Three are frameshifts, MYO6 p.Gln918ThrfsTer24 and PRRT2 p.Arg217ProfsTer8 and p.Arg217GlufsTer12, at low-complexity sites, where the flag of the exomes withheld BS2 and the joint release has no flag.
Against the exomes without the flags, 741 variants with a positive joint count lose PM2 to an AC0 in one callset, and 11 calls move, 9 of them on ClinVar pathogenic-side variants.

The change costs about 3,100 calls in direction for each added opposite-direction call (9,321 for three), where Section 2.5 puts the knee of a curve at about a hundred, and the genomes' share alone costs about 2,900 (5,848 for two).
The benchmark of this paper still reads the exomes, because every other number in it was measured on them and not because the joint release measures worse.
What the joint release gives up is the flags, and the benchmark prices only part of what they do: withholding nothing gains 3,473 calls in direction and costs one reference-discordant pathogenic call, and the flags keep three pathogenic-side frameshifts at Likely pathogenic.
It cannot show what they prevent at sites that nobody has classified, and ClinVar's classified variants are all that it measures.
Which release to build is therefore the choice of the user of a database, and `docs/ACMG_SETUP.md` gives the numbers for it.

**What these measurements do not show.**
v38 was developed on the April set, and only the 107,525 new variants of the October release were not seen when it was written: the effect of the change on them, against the commit before it, is 258 calls in direction, against 1,762 on the whole release (259 and 1,767 against v37), and they lie near April variants and were mostly in ClinVar in April (Section 2.8).
The gains are measured against ClinVar's benign-side classes, which submitters may have assigned from the same gnomAD frequencies that BA1, BS1 and BS2 read, so for those criteria the reference is not independent of the evidence (Section 3.4).
The guard against counting BS2 beside BS1 was written after TTR and VWF were seen, and both are old variants, so no variant of the benchmark tests it.
The calls for each mode of inheritance are not scored, and the benchmark has no class that could score them.
The reading is as good as the field it reads: ClinGen curates one mode for each curation, and the file of the benchmark is the export of April, which lacks 104 curations (Section 3.5).

## 4. Discussion

**What we recommend implementing, and how closely it agrees with community opinion.**
Given criteria that are underspecified and fragmented, we recommend implementing them as follows, and each recommendation points to the section that holds its numbers.
First, combine the criteria that are met by points on the pathogenic side and by the 2015 table on the benign side.
The table has no row for a lone Very Strong criterion, which the points call Likely pathogenic, and points on the benign side would call a lone Supporting criterion Likely benign, which the table does not.
Scoring both sides by points raises benign-side recall and adds hundreds of reference-discordant benign calls, and the table on both sides gives up pathogenic-side recall to avoid ten reference-discordant pathogenic calls (Section 3.3).
Second, implement each criterion with the conditions that its text states, which was the largest class of change (Section 3.2), and read the field that a source curates and not the prose beside it, because v38, which does, gives more calls in direction than v37 (Section 3.8).
Third, treat each threshold that no guideline states as a measured choice.
The points that define the classes, BP7's conservation test and the BS1 and BS2 thresholds matter, each at its own exchange rate between correct calls gained and errors added, and the PM2 thresholds under the default scheme, the cutoff for Pathogenic and the allele-number floor each move agreement in direction by fewer than a thousand calls of 673,660 (Section 3.3).
Fourth, check each input against an independent listing and keep a leave-one-out figure beside every ClinVar-derived input, because the ClinGen file that the benchmark read silently lacked CFTR, and counting a variant's own ClinVar record would add 14.3 points of pathogenic-side recall that comes from the answer itself (Sections 3.4 and 3.5).
The traps that we measured are the ones to avoid.
Counting one gnomAD observation as two benign criteria adds 1,292 calls in direction and calls two ClinVar-Pathogenic variants Benign (Section 3.8), and rules with no combining rule or safeguard reach 84.6% agreement in direction with 1,728 contradictions where v38 has 48 (Section 3.5).
With these choices v38 agrees with ClinVar's class for 63% and with its direction for 78% of the 673,660 records on which every threshold was tuned, and contradicts its direction for 48 (0.0071%).
On the 107,528 variants that were not in that set, agreement in direction is 70.3% and contradictions are 2.3 times as frequent (Sections 3.1, 3.5 and 3.7).
That is how closely the implementation agrees, and it is not a quality score: most of the disagreement is a variant that ClinVar classifies and the classifier leaves as VUS, for reasons that the input, the design and an incomplete source explain in part (Section 4.2).
The recommendation has the bounds of Section 4.4.
"Closely" is relative to the settings that we measured, and a default is not always the setting that agrees most: agreement in direction is higher with the PM2 threshold at strict absence than at the default, and BP7's conservation test, which the default keeps, gives up many correct benign-side calls to avoid two contradictions, because where a sweep did not decide the defaults were set by the rules of Section 2.5 and not by maximizing agreement (Section 3.3).
The thresholds were tuned on the April set, so the figures on the new variants are the fair ones, and they are a hold-out of variants and not of genes (Section 3.7).
ClinVar can be wrong or out of date, and the expert opinion is one reviewer's, who saw both calls (Sections 4.2 and 4.4).
A sites-only input evaluates 15 of the 28 criteria, which caps what any annotator of such an input can reach, and all of this comes from one tool, one reviewer and one benchmark.

### 4.1 Principal findings

Six results of this work seem to us worth testing on other classifiers, which we have not done: they come from one tool, one reviewer and one benchmark.

**The aggregate stopped rewarding the work.**
Agreement in direction with ClinVar reached 77.6% at v13 and stayed between 77.50% and 77.62% for the next 23 runs (v14 to v37), while the work of those runs cut opposite-direction calls from 94 to 49 and the ledger records 20 further changes; the one change since, v38, moved it by 0.22 points, to 77.84%, and cut the opposite-direction calls to 48.
A project that watched only agreement would have concluded at v13 that the work was finished, and simple rules with no combining rule or safeguard reach higher agreement with 36 times as many contradictions (Section 3.5).
The two largest steps before v13 (8.8 and 5.2 points, from figures recorded at the time) were mostly plumbing: SpliceAI and PhyloP reaching the criteria that read them, and two conditions of the guidelines that criteria had not applied (an absent variant read as absent for PM2, and BP4's confinement to splice questions).
The measures that moved late were counts of rare, severe errors, which are small numbers, and they are the ones a clinician would act on.

**Review of the discordant calls first identified the most changes, and mostly by reading the conditions a criterion's text states.**
Review first identified 30 of the 64 ledger items (47%, by our coding), and the largest class (21 of 64) was a condition of a guideline that a criterion had not applied: BS2's full penetrance at an early age, PM1's absence of benign variation, BP1's truncating-only mechanism, PVS1's need for a gene-disease relationship and for loss of function as the mechanism.
The reviewer supported ClinVar's call on 89%, 86% and 67% of the rows with a verdict in the three returned rounds (by three different instruments, Section 3.2).
A benchmark cannot find a missing condition when the omission does not change an aggregate, and several did not: building BP3's database moved agreement by 24 calls, and making the BA1 exception list match for the first time moved no cell of the matrix, because seven of the nine alleles it exists for are not in the April reference set.
Review has a sampling problem of its own, and ours was a cap of 10,000 rows in a discrepancy log that silently limited the first round to the first five chromosomes.

**A threshold belongs to the combining rule beside it.**
The PM2 threshold moved pathogenic-side recall by 19 points when it was set in v10, by 0.5 in v38 (as in v37), and by 6.1 when the 2015 table is restored on the pathogenic side (17.9 in v37, which reproduced the sweep of v10), because the table needs PM2 to lift a null variant to Likely pathogenic and the point scale does not.
The cutoffs that define the classes move the pathogenic side more than any numeric threshold that we swept on the default classifier (11 percentage points of recall for one point of cutoff, against 0.5 for the whole range of the PM2 threshold), and they are fixed by the framework and not by the data.
We showed this for the PM2 threshold, the only threshold that we re-ran under the 2015 table.
Thresholds that were tuned should be re-measured whenever the rule beside them changes, and the exchange rate between calls gained and errors added is a better unit than the agreement figure.
Re-measuring also turned up two thresholds that had never been swept and mattered: BP7's conservation threshold, which matters most on the benign side, and the count of BP1's veto, where a stricter value calls fewer ClinVar-VUS variants benign-side and avoids two reference-discordant benign calls, one of which rests on the variant's own record.

**Agreement with ClinVar partly reflects the classifier's own inputs, and the dependence that remains is small in agreement and larger in errors.**
Counting a variant's own ClinVar record would add 14.3 points of pathogenic-side recall, and the classifier excludes it from PS1 and PM1.
Removing every ClinVar-derived input costs 0.78 points of agreement and adds 30 reference-discordant benign calls (24 to 54) and takes the opposite-direction calls from 48 to 75, 28 of the reference-discordant benign calls because it removes the BP1 veto.
The veto reads other variants' ClinVar classifications to apply BP1's condition on the mechanism of the gene, so in a benchmark against ClinVar its protection is itself part of the dependence.
It also counts the variant's own record, and 3 of the 28 rest on it: with that and the two other adjustments that read the own record or followed the reviewer's rows removed together, the 24 reference-discordant benign calls are 34 and the 48 opposite-direction calls 58.
Reporting both figures is what makes the claim checkable.

**An input can be incomplete without anything showing, and a check against an independent listing can find it.**
The ClinGen file that the benchmark read lacked 104 curations, among them CFTR, and the classifier degraded in the way its design implies: a gene with no curation is not blocked, it only has no route to loss of function, so every null variant of CFTR was called VUS (Supplementary Table S7 lists them).
The design is the defect worth naming, because a gene that is absent from the curated file silently has no route, and an explicit flag for a gene that is not in the file would have listed it, among 7,391 of the benchmark's 10,231 genes.
Nothing in an aggregate showed it, and the review could not, because a VUS is not an opposite-direction call; even the benchmark's four-star records, all of them CFTR, were called VUS for 22 of 23.
We found it by asking why one gene had no PVS1.
A mechanical check would have found it sooner: list the genes in which ClinVar's pathogenic truncating variants have no route to PVS1, and compare them with an independent listing.
That check flags 625 genes (576 of them not in the file), of which the completed ClinGen file explains 23 and the file with Orphanet's label 170, so it is a prompt for review and not a detector, and here one input defect sat among many gaps of design.

**A flat sweep can mean that the setting reaches almost nothing, or that the rule beside it makes the setting irrelevant, and an input read from prose can hide the first.**
The sweep of the BS2 carrier threshold was flat in v37 because the classifier took inheritance from the words of the disease names, which made 50 genes dominant-only where the curated field in the same string makes 746; read from the field the threshold is a trade (Sections 3.3 and 3.8).
The sweeps of the PM2 thresholds were flat for the other reason, the combining rule (Section 3.3), and they stayed flat when the field raised the reach of the recessive-only threshold from 67 genes to 1,310.
Across the two thresholds the words made 117 of the 2,924 genes of the ClinGen file dominant-only or recessive-only, and the curated field makes 2,056.
Reading the field exposed a cost that no sweep of v37 could show, one gnomAD observation counted as two Strong benign criteria: when BS2 first read the field, on a patched build of v37, two ClinVar-Pathogenic variants became Benign, and they did again on the reading as it first ran in v38, until BS2's carrier count was set aside beside BS1 and BA1 (Section 3.8).
A flat curve is an observation about the calls that a setting reaches, and it should be read with the number of calls it reaches and the rule beside it.

### 4.2 Reading the discordance

A ClinVar agreement of 63% exact and 78% in direction is not a quality score, and the decomposition of Table 5 shows why.
Almost none of the discordance is contradiction: 48 of 673,660 records (0.0071%) are opposite, and 24 of them are reference-discordant benign calls.
Almost all the rest is a defined ClinVar class that the classifier leaves as VUS (140,438 records).
Part of that is evidence that a single-sample VCF does not carry (segregation, de novo status, phenotype specificity, functional assays, case counts), and 13 of the 28 criteria cannot fire in this benchmark.
Part is design: 10,693 of the 35,182 pathogenic-side variants called VUS sit one point below the Likely pathogenic cutoff, 5,184 are the PVS1 gap, and of the 105,256 benign-side variants called VUS 21,614 carry a lone BP4 or BS1, which the 2015 table on the benign side does not call, and 32,168 carry PM2 only.
And part is an input that was incomplete (Section 3.5).
We cannot apportion the three from the aggregate or from the VCF.
A VUS from an annotator that has these inputs is a statement about the inputs, and it does not mean what a VUS from a laboratory means.
A further 98,358 records are in the right direction at the wrong tier, and for the pathogenic-side half the boundary falls one point above where a lone null variant with PM2 lands: 41,549 of the 48,277 ClinVar Pathogenic variants called Likely pathogenic are at 9 points.
The classifier's calls are therefore best read with their error profile beside them: of 62,238 pathogenic-side calls 94.2% are pathogenic-side in ClinVar, 5.8% are ClinVar VUS and 0.04% are benign-side in ClinVar, and none of the 4,038 Pathogenic and the 67,733 Benign calls contradicts ClinVar in this in-sample set.

The residual contradictions are informative because they are few and each was reviewed by one medical geneticist, whose rulings are one expert's reading and by our coding.
In the last round the reviewer supported ClinVar on 60% of the 48 and fastVEP on 8%, with a different class or more evidence needed for the rest, and the reasons they gave were rarely about the arithmetic of the criteria.
They were the frequency a disease model allows, the phase of a variant, a truncation that escapes NMD in the first exon, an isoform, a transcript that predates MANE, and a recessive disorder.
Each is evidence or a policy that a single-sample, single-transcript annotator can only approximate, and several can be approximated better (a disease-specific frequency threshold, a first-exon rule, a transcript policy).
We have not done so here, because each is a decision about the classifier and not about this paper's measurements, and because the residual is the sample on which the classifier was reviewed: a fix tuned on it could not be tested on it; the hold-out of Section 3.7 is the test that the data allow, and the changes of Section 4.5 that could be run there do not fail on variants the review never saw.

ClinVar is a reference and not ground truth.
Its review status records how many submitters agree, or that an expert panel reviewed the record, and is not an independent confirmation that a variant is pathogenic or benign.
Its classifications can rest on the same population-frequency, computational and consequence evidence that the software reads, so agreement can be partly circular, which Section 3.4 measures for the ClinVar-derived inputs and not for REVEL or SpliceAI.
Three observations in this work bear on how far to trust it.
Five of the 48 residual variants (10%) had left the 2-star+ set by October, all with conflicting classifications, against 1.3% of all April variants that left it for any reason.
ClinVar moved 157 of 665,057 variants across a direction in six months.
And the label of a variant can depend on the haplotype it lies on: for BRCA1 c.594-2A>C the reviewer's note records that the variant is always in cis with c.641A>G and that the haplotype, not the variant, has been shown not to be pathogenic, which no single-variant annotator could reproduce.
A discordance with ClinVar is a prompt to look, and the reviewer's reasons suggest that it is often a prompt about evidence and not about the arithmetic of the criteria, which one reviewer's reading of 48 rows cannot show for the other discordance.

### 4.3 Practical lessons

Several lessons are of practical use to anyone building a classifier of this kind.
Test a criterion with the pipeline's own output and not a hand-written input: the BA1 exception matcher passed its unit tests because they fed it the bare token that the pipeline never produces.
Check that a benchmark can see what is being tested, because it cannot validate a criterion whose population it excludes.
Keep a leave-one-out mode and report the ClinVar-informed figure beside it, so that a reader can see the size of the circularity, and apply it to every ClinVar-derived input, which we did not do for BP1's veto.
Report a baseline of simple rules, and the contradictions, beside the agreement figure, because agreement can be exceeded by calling aggressively.
Record what produced every annotation, and refuse to reuse one that something else produced, because a cache that was right for an older parser is intact and is not safe to read, as the repository's own transcript cache showed when a parser was fixed (issues 95 and 98).
Check each input against an independent listing, because a missing record degrades quietly.
Read the field that a source curates and not the prose beside it, and test the reading on the source's own rows: the words of the disease names made 117 genes dominant-only or recessive-only where the curated field makes 2,056, the sweep of the BS2 carrier threshold that reached 50 of them was flat, the effect of the field on a call was seen only in a measurement of a patched build, and a test over 185 genes of real ClinGen rows now fails if the reading moves.
Check what a second source of the same kind does to the first: a second gnomAD database in the same directory answered in place of the first without a message until the command line reported it.
Count severe errors as a primary measure beside agreement, and say which errors weigh more: by that count v10 (25 opposite-direction calls) was better than v38 (48), and it was given up for 6,632 correct pathogenic-side calls, which is a trade and not an improvement.
And give the reviewer the whole row: the first workbook carried 75 columns of evidence, and the list of discordance categories that the reviewer attached to round 1 anticipates three classes of the ledger.

### 4.4 Limitations

The study has limits that the reader should weigh.
A single medical geneticist reviewed the calls, was not blind to either call, supplied inputs as well as rulings, saw only the opposite-direction calls and none of the 35,182 pathogenic-side variants that we left as VUS or any concordant call, and returned a marked workbook in three of six rounds; the round-6 categories are our coding of their text, and we have not asked the reviewer to confirm our coding of rounds 2 and 6.
The ledger of 64 items and its classes are our own coding of the project's record and have not been independently coded.
All thresholds were tuned on the April benchmark and its discordances, so the April figures are in-sample (the 48 are 58 when the three adjustments that read the own record or followed the reviewer's rows are removed), and the October hold-out is a hold-out of variants and not of genes, of a release six months later, whose new variants were mostly in ClinVar in April at a lower review status.
The price of an error was not applied uniformly: the trades accepted at v11 and v13 cost 214 and 315 correct calls per added error (on 31 and 36 errors), and at the BS2 and BS1 thresholds we declined trades at 401 and 821 (on 3 and 6 errors), and in v38 we declined the label guard that restores two ClinVar-Pathogenic variants at about 32 correct benign-side calls for each, below the knee of Section 2.5, because it reads the name of a disease (Section 3.8); each decision rests on its own recorded reasons.
The ClinVar-derived inputs of the classifier were built from releases of 29 April and 27 June 2026 and the reference is the VCF of 4 April, so the three are not one snapshot.
The reference is a selection (two stars and above, SNVs and indels of at most 100 bases): it excludes records with conflicting classifications (six of the nine BA1-exception alleles of Ghosh et al. are among them in later releases, and a seventh, HFE, was absent from April and is in October) and records with a single submitter.
It also has no negative controls, and its mix of classes (25% pathogenic-side among the defined classes) fixes the precision figures: a different mix would give different precision for the same calls.
A variant's call is the call on one picked transcript, and other transcripts can differ.
The classifier evaluates 15 of 28 criteria on a single-sample VCF without phenotype, so its recall of ClinVar's classes has a ceiling that is not a fault of the thresholds.
The benchmark's classifier reads the gnomAD exomes only, and PM2 reads absence from them as rarity, including at sites that the exomes do not cover; requiring a record costs 0.22 points of agreement and 2.9 points of pathogenic-side recall (Figure 3A), and PM2 fires on 66% of the variants.
gnomAD's genomes are not read by the benchmark; read through the joint release on every chromosome they give 9,321 more calls in direction and three more opposite-direction calls, and the release lacks the flags of low-complexity regions and segmental duplications, whose protection a benchmark of ClinVar's classified variants cannot measure (Section 3.8).
The exchange rates of Section 3.3 are counts of differences between deterministic runs and carry no sampling interval: a rate computed from one to three added errors is as uncertain as a count of that size (the 95% Poisson interval of a count of 2 is 0.2 to 7.2).
The offline replay of the benign-side point scales was not checked by real runs, and the cutoff and measurement builds are patched builds of the default code.
We made no comparison with another tool, and the bulk outputs of the runs (annotated VCFs of 250 to 285 MB and per-variant tables of about 20 MB, compressed) are not archived with the repository; its scored matrices and provenance are.
The mode of inheritance that BS2 and PM2 read was the weakest input of v37, which took it from the names of the diseases: of the 2,924 genes of the ClinGen file it read 67 as recessive-only, 50 as dominant-only and 12 as both, and of the 740 genes with a benchmark variant whose established ClinGen tags are dominant only, 39 registered as dominant, and of the 1,285 whose established tags are recessive only, 60 registered as recessive.
v38 reads the curated field, and every number of Sections 3.1 to 3.7 rests on it, the sweeps of the recessive-only PM2 threshold and the dominant-only BS2 threshold of Section 3.3 among them; Section 3.8 sets them against v37's.
v38 has limits of its own.
The field gives a mode and not a penetrance, so BS2's carrier arm applies to every gene read as dominant, and four ClinVar-Pathogenic loss-of-function variants lose their call to it.
The calls for each mode of inheritance of a gene are not scored, because a ClinVar record of the benchmark has no disease-specific class.
BS2's arm for homozygotes and hemizygotes still counts, beside BS1 and without BA1, 1,301 variants that BS1 and BS2 call Benign on one frequency, none of them a ClinVar pathogenic-side variant (Section 3.8).
And the guard against counting BS2 beside BS1 was written after two of the variants it protects had been seen, so no variant of the benchmark tests it.
Three limits of the implementation remain in v38, the release that this paper describes, each to a different degree.
Inheritance, which v37 took from the words of the disease names, is read from the mode that ClinGen curates for each of the 2,924 genes of the ClinGen file; a gene that is known only from an OMIM phenotype has no tag and is still read by its words, and so is a gene that the file does not hold (Section 3.8).
PM2 still reads a variant with no gnomAD record as absent, including at sites that gnomAD's exomes do not cover, because no coverage is read (the setting `pm2_absent_when_no_record` turns that off, at the cost given above), and reading the genomes through the joint release does not change that.
And the ClinGen file is an input that the release does not ship: the converter checks the format of the export that it is given and cannot tell that curations are missing from it, so an incomplete export is read as a complete one, which is how the benchmark came to read one that lacked CFTR (Section 3.5); the setup guide tells the user to download ClinGen's current export, and the mechanical check that would have found the gap sooner (Section 4.1) is an analysis of the harness and not a part of the release.
Where the project's documentation attributed a number or a grade to a published recommendation that we could not find stated there, or where the implementation departs from the recommendation, we say so in Sections 3.2 and 3.3.
fastVEP's classification is a decision-support output for research use: it is not a laboratory-validated classification and does not replace expert review.

### 4.5 Next steps

The findings suggest an order of work, and the measurements of Sections 3.3 to 3.5 give each step a price; Section 3.7 tried the changes that could be run on the 107,528 new variants.
First, complete the ClinGen input, which adds 852 correct pathogenic-side calls (848 in direction) for 1 reference-discordant pathogenic call on the April set and 51 (50 in direction) with no opposite-direction call on the new variants; the check of Section 4.1 and a flag for a gene that is not in the curated file were not measured, and the flag would list 7,391 of the benchmark's genes.
Second, a source of established loss-of-function mechanism that does not depend on constraint (the PVS1 applicability stated in VCEP specifications, or Orphanet's label, which adds 272 calls in direction for 2 reference-discordant pathogenic calls on the new variants), together with a way to tell which disease a variant is pathogenic for (v38 classifies a gene once for each mode of inheritance, which does not tell two diseases of one mode apart), because 3,470 ClinVar truncating variants would still carry no PVS1 after both inputs are used.
Third, BP7 as the splicing recommendation describes it, which gains 15,552 calls in direction for 2 reference-discordant benign calls on the April set and 4,859 for none on the new variants, and which only the maintainers can weigh.
Fourth, BP1's veto leaving out the variant's own record (3 more reference-discordant benign calls on the April set and none on the new variants) with a stricter count: the pair at a count of 1 gives 1,735 more calls in direction and 127 fewer correct benign-side calls on the April set, with 23 reference-discordant benign calls against the 24 of v38, and the count of 1 alone gives 124 more calls on the new variants and avoids no reference-discordant benign call there.
Fifth, a transcript policy for MANE Plus Clinical, a first-exon rule for PVS1, and disease-specific frequency thresholds for BA1 and BS1 where a VCEP has published one, which the classifier already supports as an opt-in table.
Sixth, review of the allele-number floor, whose source is narrower than its use and whose removal would add 728 calls in direction on the April set and none on the new variants, and of the two grades and the NMD default of PVS1 (Section 3.2).
Seventh, gnomAD's genomes, read through the joint release, which gives 9,321 more calls in direction for three more opposite-direction calls, about 3,100 for each, and gives up the flags of low-complexity regions and segmental duplications (Section 3.8).
The prices of all seven steps are measured on v38; a completed ClinGen file also changes the inheritance that v38 reads for the genes it adds (CFTR is autosomal recessive), which changes the criteria of some variants (Section 3.8).
Each of these can be measured by the harness of this paper before it is adopted, and none should be tuned to the 48, which is why the benefit of a stricter count of BP1's veto, which avoids two of them, is in-sample.

## 5. Conclusions

The contribution of this work is how to implement underspecified and fragmented ACMG/AMP criteria in software so that the calls agree closely with community opinion, here ClinVar's classifications of review status two stars or better and one reviewer's rulings, and what each choice costs in contradictions.
The choices, the agreement and the contradictions that each measured, and the checks that recompute the numbers that the text quotes are in this paper for reuse, within the bounds of Section 4.4.

Building a software classifier for the ACMG/AMP standards was largely the work of applying conditions that the standards state and a first implementation omitted (21 of 64 changes), and review of the discordant calls first identified more of the changes (30, by our coding) than any other route.
On the 673,660 ClinVar 2-star+ variants on which every threshold was tuned, v38 agrees with ClinVar's direction for 77.8% and contradicts it for 48 (0.0071% of the records and 0.013% of those with a pathogenic-side or benign-side class), which are 58 under a strict leave-one-out reading, three of the added calls Benign; none of its 4,038 Pathogenic or 67,733 Benign calls is among the 48, and rules with no combining rule or safeguard reach higher agreement with 36 times as many contradictions.
On the 107,528 variants that were not in the April set, agreement in direction is 70.3%, mostly because of the benign side, and contradictions are 2.3 times as frequent (interval 1.2 to 4.1), more than their mix predicts at the April rates, and of the same coarse kinds.
The thresholds that no guideline states are less decisive than the combining rule around them for the PM2 threshold, the one we tested under the 2015 table, and two that had never been swept, BP7's conservation test and the count of BP1's veto, each exchange calls and errors at a rate that the maintainers can now weigh, which we have not changed; on the pathogenic side the points that define the classes move the result more than the thresholds that we swept on the default classifier.
Against v37, v38 reads the mode of inheritance that ClinGen curates, which v37 took from the words of the disease names, and gives 1,494 more calls in direction on the April set and 1,767 on the covered October release, 259 of them on the new variants, with no added reference-discordant benign call and one fewer reference-discordant pathogenic call, though four ClinVar-Pathogenic loss-of-function variants move from Likely pathogenic to VUS; it does not count BS2 beside BS1 or BA1, which would have added 1,292 calls in direction and called two ClinVar-Pathogenic variants Benign.
Reading gnomAD's genomes through its joint release adds 9,321 calls in direction on every chromosome and three opposite-direction calls, and gives up the flags that withhold the frequency criteria in repeats.
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
The analysis of this paper is in the repository under `analysis/acmg_benchmark/` (the benchmark harness `run_benchmark.sh`, the run record `RUN_VERSIONS.md`, and the scripts `11_manuscript_runs.py`, `12_per_variant_table.py`, `13_manuscript_stats.py`, `14_temporal_holdout.py`, `15_clingen_complete_export.py` and `16_compare_runs.py`, which compares a rerun with the runs it repeats, `20_compare_families.py`, which sets the runs of v38 against the same runs of v37, and `build_manuscript_binaries.sh`, which builds the binaries) and under `manuscript/acmg/` (this manuscript, the tables it reads in `data/` and the workbook `supplementary_tables.xlsx` that gathers them, `build_tables.py`, `build_supplement.py`, `verify_numbers.py`, `generate_figures.py`, the figures, and the patches `replay_validation/combiner_cutoff_env.patch` and `replay_validation/bp_hooks_env.patch` that expose the two point cutoffs, and BP1's own-record count and BP7's SpliceAI threshold, for the runs of Table 4 and of Sections 3.3 and 3.4).
The scored concordance matrix, the criterion firing counts and the provenance file of each run made for this paper (96 runs of v38, 89 on the April set and 7 on the October release, and the same 96 runs of v37) are in `analysis/acmg_benchmark/results/v37_manuscript/`, a directory that keeps the name it had when it held the runs of v37 alone; the provenance records the flags, the SHA-256 of the binary and of every configuration, the inputs read, and the commit and the number of uncommitted files of the checkout that ran the harness, which is not the commit of the binary: `binaries_v38.tsv` in the same directory says what each of the three binaries of v38 was built from, and `binaries.tsv` what each of the three of v37 and of the three builds of the working tree that made the runs of Section 3.8 first was built from.
The binaries of v38 are three builds of the tag `v0.5.0-rc.1` (`c12d9f9`, release 0.5.0): the plain one, and the two with the patches that expose the point cutoffs and the BP1 and BP7 settings (`build_manuscript_binaries.sh`).
The annotated VCFs (about 250 MB each), the per-variant tables, the annotation databases and the ClinVar, gnomAD, ClinGen and Ensembl inputs are not in the repository.
Table 2 and Section 2.2 give their sources and versions, and the ClinVar files are the weekly releases dated 4 April and 4 October 2026.
Orphanet's association table is the public Orphadata product 6 (`build_gene_disease_attributes.py` gives the command).
The ClinGen export that completes the file of 30 April is a public download of 11 October 2026 (`15_clingen_complete_export.py` gives the command), and the 112 rows it adds, 104 of them curations that the April file lacked, are listed with the status of each in `manuscript/acmg/data/clingen_missing_curations.csv`.
The runs that Section 3.8 made on its own, 38 on the April set and 4 on the October release, are in the same files (`runs_v38.tsv`, `holdout_runs_v38.tsv` and `run_matrix_v38.tsv`, with their provenance files), and 18 of the 38 are sweeps that belong to the family of 96 as well, so the family adds 71 runs on the April set and 7 on the October release to them.
The runs of Section 3.8 were made with builds of the working tree that became v38, so the source of those binaries is that tree and not a commit.
The classifier of v38 is the commit `5293cc0`, merged with its documentation and two speed-ups in pull request #167 (`6504819`), and the 38 runs on the April set and the 4 on the October release were made again on 7 November 2026 by binaries built from the tree of release 0.5.0 (the tag `v0.5.0-rc.1`, `c12d9f9`): all 42 agree with these runs in every scored file, and the 15 annotated VCFs that both keep agree record by record (`analysis/acmg_benchmark/results/v38_v0.5.0/rerun_comparison.csv`; `analysis/acmg_benchmark/RUN_VERSIONS.md`, "Release 0.5.0 reproduces v38").
Two later changes that make the annotation faster (`RUN_VERSIONS.md`, "After v38") reproduce the annotated VCF of the `v38` run record for record, so the code that is merged gives the calls reported in Section 3.8.
The alternatives of Section 3.8 are read from the environment by a build of `manuscript/acmg/replay_validation/inheritance_alternatives_env.patch`.
The real ClinGen rows that test the reading of the field are in `crates/fastvep-cli/fixtures/clingen/` and are made by `analysis/acmg_benchmark/scripts/sa_sources/make_clingen_parse_fixture.py`, and the databases of gnomAD's exomes, genomes and joint release are built by `build_gnomad_release_per_chrom.sh` in the same directory.
`17_joint_boundary_cases.py` annotates the variants that the joint release turns into errors, or out of a pathogenic-side call, through both databases.
`18_criteria_status.py` annotates every 40th record of the benchmark to JSON and writes the status of the 28 criteria (Supplementary Table S1).
The reviewer's workbooks are not in the repository, and the remarks of the reviewer that Sections 3.2 and 3.6 report are read from them and are not checked by `verify_numbers.py`.

The results of the two versions are in these places.
The results of v38, which is the version of every number of this paper unless it names v37, are the 96 runs of v38 in `analysis/acmg_benchmark/results/v37_manuscript/` (`runs_v38.tsv`, `holdout_runs_v38.tsv`, `run_matrix_v38.tsv` and `binaries_v38.tsv`, with the provenance and the configuration of each run in the subdirectories provenance and configs, where the files of these runs are named v38_ and the name of the run), their repetition by release 0.5.0 in `analysis/acmg_benchmark/results/v38_v0.5.0/` (the 42 runs of Section 3.8, 18 of which are the sweeps that they share with the family, all identical to the first runs), the 67 tables of `manuscript/acmg/data/` that have a twin for each version (`headline_v38.csv` and the others, `_v38` in the name), the 23 files whose names begin with v38_ (Appendix Tables A10 to A16, which are sheets A10 to A16 of `supplementary_tables.xlsx`), Tables 3 to 9, Figures 1 to 5 and S1 to S4, Supplementary Tables S1 to S9 and Appendix Tables A1 to A16, and the sections "v38", "The full rerun of 4 November 2026", "After v38", "Release 0.5.0 reproduces v38" and "The family v38" of `RUN_VERSIONS.md`.
The results of v37, the before of every difference that v38 makes, are the same 96 runs of v37 in the same directory (`runs.tsv`, `holdout_runs.tsv`, `run_matrix.tsv` and `binaries.tsv`) and the 67 twin tables with `_v37` in the name, which `13_manuscript_stats.py --base v37` writes and `build_tables.py --base v37`, `build_supplement.py --base v37 --out` and `generate_figures.py --base v37 --out` turn into the tables, the workbook and the figures of v37; `verify_numbers.py --base v37` lists, sentence by sentence, the numbers of the text that v38 changed.
`analysis/acmg_benchmark/results/v38_family_vs_v37.md`, with its five CSV files, sets the two versions against each other run for run.
Seven files of `manuscript/acmg/data/` belong to neither version: `change_ledger.csv`, `review_rounds.csv` and `trajectory.csv` hold the history of the numbered versions, `round6_rulings.tsv` is the reviewer's own ruling, `clingen_missing_curations.csv` lists the rows of the ClinGen export, `criteria_status.csv` is read from the engine of v38, and `tavtigian_posterior.csv` is arithmetic of the framework.
The annotated VCFs and the per-variant tables of the runs are in `data/benchmark/manuscript_runs/` on the machine that made them, and are not in the repository.
The analyses of Sections 3.1 to 3.7 were made on v38 and on v37; what was not repeated on v38 is what needs a reviewer or a run that the repository does not have, which is the review rounds, the ledger and the runs of round 6, which need builds that the repository does not hold.

## Declarations

**Author contributions.** M.G. carried out the clinical review of discordant calls over six rounds and contributed curated variant inputs.
W.C. advised on the clinical interpretation of a subset of cases.
K.-L.H. designed and implemented the classifier, ran the benchmark and the measurements, and wrote the manuscript.
All authors read and approved the final version.
[To be confirmed by all authors.]

**Use of AI assistance.** [To be completed by the author.
The analysis scripts, the figures and the first draft of this manuscript were produced with the assistance of Claude (Anthropic) working as a coding agent in the project repository, and the author has not yet reviewed them in full.]

**Funding.** [To be completed by the author.]

**Competing interests.** [To be completed by the author.]

**Ethics.** The study uses public records of variants and no patient data, and no human participants were recruited.
The reviewer is an author (M.G.), the review workbooks are not public, and the manuscript paraphrases rather than quotes the rulings they record.

---

## Figure legends

In every figure, a ClinVar class (P/LP, VUS and B/LB, or pathogenic-side and benign-side) is the class of a ClinVar record of review status two stars or better, which the labels call "ClinVar 2-star+": criteria provided, multiple submitters, no conflicts; reviewed by expert panel; practice guideline.
The release is that of 4 April 2026, and panel C of Figure 5 also uses the release of 4 October 2026.
P/LP is Pathogenic or Likely pathogenic, B/LB is Likely benign or Benign, and VUS is uncertain significance.
Each figure repeats this in a footnote.
ClinVar's classifications are a reference and not ground truth.

**Figure 1. Thresholds that no guideline states trade calls for errors at different rates, and the PM2 threshold matters far more when the pathogenic side is scored by the 2015 combining table than when it is scored by points.**
The thresholds of Section 3.3, re-measured on v38.
(A) Five thresholds on the benign side trade correct benign-side calls for reference-discordant benign calls (the pink line is the count of BP1's veto, from 1 to 25 pathogenic missense variants); open markers are swept values and the filled black marker is the default one.
(B) Recall of ClinVar 2-star+ pathogenic-side variants against the PM2 threshold for genes not read as recessive-only, under the default scheme, which scores the pathogenic side by points ("Points"), and under the 2015 combining table on both sides ("2015 table"), and against the threshold for genes read as recessive-only under the default scheme ("Points, recessive-only").
The benign side is the 2015 table in all three, so only the pathogenic side differs.
(C) The BS2 carrier threshold for genes read as dominant-only (v38, which reads 746 genes as dominant-only, where v37 read 50): the change in recall of ClinVar 2-star+ pathogenic-side (red) and benign-side (blue) variants from the default threshold.

**Figure 2. The points that define the classes move the pathogenic side more than any numeric threshold that we swept on the default classifier.**
(A) For variants with no benign evidence, the ClinVar 2-star+ classes of the variants at each number of pathogenic points; the lines mark Likely pathogenic from 6 points and Pathogenic from 10.
(B) Real runs with the Likely pathogenic cutoff changed from 4 to 9 points: recall of ClinVar 2-star+ pathogenic-side variants against the number of ClinVar 2-star+ VUS and benign-side variants called pathogenic-side.
(C) The benign side under the 2015 table (default) and under a point scale with Likely benign from -1, -2 and -3 points, replayed from the recorded criteria and not real runs.
In B and C a filled black marker is the default setting.

**Figure 3. Three evidence sources carry the benign side and one carries most of the pathogenic side, most design decisions pay for themselves, and each departure from the sources has a price.**
What each design decision and each evidence source buys, and what each departure from the sources would change.
(A) Sixteen settings changed one at a time (switched off, or on for the two opt-in settings): the change from v38 in agreement in direction (percentage points of all records), in recall of ClinVar 2-star+ pathogenic-side variants, in reference-discordant benign calls and in reference-discordant pathogenic calls.
(B) Nine evidence sources removed one at a time from the annotation directory.
(C) Eight runs that test a departure from the sources or an adjustment: BP1's veto leaving out the variant's own ClinVar record, BP7's SpliceAI threshold and conservation test (alone and together), the ClinGen file with the curations that its export lacked, PVS1 licensed by Orphanet's loss-of-function label (alone and with the ClinGen file), and the own-record veto, the low-penetrance gate and the exception list removed together.
The BP1 run, the SpliceAI run, the run with both BP7 changes and the run with the three adjustments use the measurement build of Section 2.6.

**Figure 4. Almost all of the discordance of v38 with ClinVar 2-star+ classifications is one side left undecided, and contradictions are rare.**
(A) The five-tier concordance matrix: counts, and the share of its ClinVar 2-star+ class (the row) that each cell holds; the outlined cells are opposite-direction calls.
(B) What v38 says of each ClinVar 2-star+ class, as P/LP, VUS and B/LB; 382 records with no call (0.13% of the ClinVar 2-star+ benign-side class) are not drawn.
(C) What ClinVar 2-star+ says behind each v38 call.

**Figure 5. Recall varies with the consequence class and the review status, and on the new variants it falls mostly on the benign side.**
(A, B) The share of ClinVar 2-star+ pathogenic-side and benign-side variants that v38 calls P/LP, VUS and B/LB, by consequence class and by review status; an asterisk marks a row whose bar stops short of 100% because records have no call.
(C) Recall in direction for each ClinVar 2-star+ class on the April 2026 release, on the October 2026 variants whose class did not change, and on the new October 2026 variants that were not in the April set; whiskers are exact binomial 95% intervals, which treat variants as independent.

---

## Supplementary materials

The supplementary figures show how the classifier was developed and reviewed, and the nine supplementary tables are the tables that summarize the data best; every other table is in the appendix.
All of them are in one workbook, `manuscript/acmg/supplementary_tables.xlsx`, which has a sheet of contents and a sheet for each table (S1 to S9 and A1 to A16).
Each table is also a set of CSV or TSV files in `manuscript/acmg/data/` of the repository, which the entries below name, and each of those files is written by `analysis/acmg_benchmark/scripts/13_manuscript_stats.py`, `11_manuscript_runs.py`, `14_temporal_holdout.py` or `15_clingen_complete_export.py` from the scored runs.
The files keep the names that the scripts read: a column `false_benign` or `false_pathogenic` counts the reference-discordant benign or reference-discordant pathogenic calls of the text, `shipped` marks the default setting, and `bar` and `far boundary` in a column or a setting name are a threshold and BP7's intronic distance.

### Supplementary figures

In Figures S1, S3 and S4 a ClinVar class is that of a ClinVar 2-star+ record of the release of 4 April 2026, as in the main figures, and each figure says so in a footnote.
Figure S2 shows no ClinVar class.

**Figure S1. Agreement with ClinVar 2-star+ classifications across the benchmark runs v1 to v38.**
(A) Same-direction and exact agreement as a percentage of 673,660 variants.
Open markers are v1, a matrix reconstructed from the documentation of the run, and v2 and v4 to v6, which were recorded at the time; no raw output of them was kept.
v1 annotated 627,376 of the records, and its 54.7% of those is 51.0% of all records.
The cross marks v17, the opt-in run with the published VCEP thresholds, which is not a step on the way and is left out of the lines.
(B) Opposite-direction calls on a log scale: ClinVar 2-star+ pathogenic-side variants called benign-side (reference-discordant benign calls), ClinVar 2-star+ benign-side variants called pathogenic-side, and their total; open markers as in A, and v6 is the recorded total.
(C) Same-direction recall for ClinVar 2-star+ pathogenic-side and benign-side variants; the values of v1 (10.9% and 19.0%) are below the axis.
The 13 numbered markers are the changes listed under the panels.
Run v3 was a partial run that was overwritten, v21 and v22 repeat v20, and v24 and v25 were intermediate runs.
ClinVar's classifications are a reference and not ground truth.

**Figure S2. The 64 recorded changes by class and by the route that first found them.**
Each cell counts changes in the ledger of Supplementary Table S3.
The column on the right gives the total of each class, and the row beneath the cells the total of each route.
Classes and routes are our coding of the project's record.

**Figure S3. What the reviewer ruled.**
(A) The rulings on the rows of the three workbooks the reviewer returned, in rounds 1, 2 and 6, in five categories: the ClinVar 2-star+ call supported, fastVEP's call supported, neither (the reviewer's call is VUS or another class), more evidence needed or dependent on the disorder, and no verdict recorded.
Rounds 1 and 2 are decoded from cell colors against each workbook's key, and round 6 is our coding of the reviewer's text (Appendix Table A4).
(B) The 48 opposite-direction calls of v38, by the evidence behind the call, with the rulings of the last round.

**Figure S4. Agreement in direction against opposite-direction calls.**
Each point is a classifier, on 673,660 variants: the versions v7 to v38 with stored scored matrices (blue, joined in order), the runs of this paper and the earlier runs of the dominant-gene BS2 settings, each v38 with one change (gray), every variant called VUS (the diamond at zero), and three rules with no combining rule, gene-level gate or safeguard at nine frequency thresholds each, with the second rule alone and the second with the third also given with no frequency rule (the filled markers).
R1 calls a variant benign-side when its population-maximum allele frequency is at or above the threshold (the thresholds at the two ends of each line are labeled, and the point at 0.001 of the last rule is annotated), R2 adds a pathogenic-side call for a stop-gained, frameshift, canonical splice or start-lost variant that is absent from gnomAD or below 1e-4, and R3 adds a benign-side call for a quiet variant (synonymous, intronic, untranslated, flanking or non-coding) with SpliceAI below 0.2.
The horizontal axis is on a log scale above 10.

### Supplementary tables

A ClinVar class in these tables is the class of a ClinVar 2-star+ record of the release of 4 April 2026, and of 4 October 2026 where a table says that it is about the later release.

**Supplementary Table S1.** The 28 criteria and the benchmark: what each needs, whether the benchmark can apply it and, on every 40th record of the benchmark annotated to JSON, how often the engine could evaluate it and how often it met it (`criteria_status.csv`, made by `18_criteria_status.py`); and how often each of the 15 computable criteria fires in v38, at which strengths, and on variants of which ClinVar class (`criteria_firing_v38.csv`).

**Supplementary Table S2.** The versions and the runs of this paper.
For 37 of the 38 versions (v3 was a partial run that was overwritten before it was preserved), the agreement, the opposite-direction calls and the recall of each, with the number of variants each committed to, the denominator of each, and the marker of v17 (`trajectory.csv`; v1 is a reconstructed matrix and v2 and v4 to v6 are the figures recorded at the time).
For v38 itself, the 89 runs made for this paper and the 8 earlier runs of the dominant-gene BS2 settings, with the main run (98 rows), the concordance matrix, the firing counts and the settings (`runs_v38_variants.csv`; the run matrix, the matrices of the seven runs on the October release and the provenance files are in `analysis/acmg_benchmark/results/v37_manuscript/`).

**Supplementary Table S3.** The ledger of 64 changes: run, class, route of discovery and source in the project record (`change_ledger.csv`).

**Supplementary Table S4.** What each step of each threshold sweep buys and costs (`threshold_steps_v38.csv`), and the median, quartiles and range of the published VCEP thresholds (`vcep_bars_v38.csv`).

**Supplementary Table S5.** Concordance of v38 with ClinVar 2-star+ classifications: the five-tier and three-tier matrices, what ClinVar says behind each call, the concordance by consequence class and by review status, the 23 four-star records, and the headline rates with their intervals (`confusion_5x5_v38.csv`, `confusion_3x3_v38.csv`, `what_is_behind_each_call_v38.csv`, `by_consequence_v38.csv`, `by_review_status_v38.csv`, `four_star_variants_v38.csv`, `headline_v38.csv`).

**Supplementary Table S6.** The rules with no combining rule, gene-level gate or safeguard and what they agree with ClinVar on, at one threshold and at nine, beside v38; the second rule against the PVS1 gap; and what each tier of the v38 call is worth, with the rate of opposite-direction calls and its interval (`baselines_v38.csv`, `baselines_curve_v38.csv`, `baselines_r2_vs_v38.csv`, `precision_by_tier_v38.csv`).

**Supplementary Table S7.** The PVS1 gap and the ClinGen file: truncating variants without PVS1, by gene and by the reason PVS1 is absent (`pvs1_gap_v38.csv`, `pvs1_gap_genes_v38.csv`, `pvs1_gap_reasons_v38.csv`); the 112 rows that complete the ClinGen export of 30 April, 104 of them curations that it lacked (`clingen_missing_curations.csv`, `clingen_missing_summary_v38.csv`); the 558 null variants of CFTR, each with its transcript, its call in v38 and the call that the completed file gives (`cftr_null_variants_v38.csv`); and what the three runs that change an input do to the gap (`pvs1_gap_closure_v38.csv`).

**Supplementary Table S8.** The 48 opposite-direction calls of v38 with their criteria, gnomAD counts and rulings, and the rulings by kind of call (`opposite_direction_v38.csv`, `opposite_direction_by_group_v38.csv`).

**Supplementary Table S9.** The later ClinVar release: the strata, recall by class, coverage, the opposite-direction rates and their ratio with its interval, the 18 opposite-direction calls on the new variants, and the seven proposed changes run on the new variants (`temporal_v38.csv`, `temporal_by_class_v38.csv`, `temporal_coverage_v38.csv`, `temporal_opposite_rates_v38.csv`, `temporal_ratio_v38.csv`, `temporal_new_opposite_v38.csv`, `temporal_proposals_v38.csv`).

## Appendix

Further tables, in the order in which the text cites them.
A ClinVar class in these tables is, as in the supplementary tables, the class of a ClinVar 2-star+ record.

**Appendix Table A1.** The measurement builds and the provenance of every run: the comparison of the controls of the two builds with v38 (`controls_v38.csv`), the calls that became opposite-direction in the runs that test BP1's veto and BP7 (`hook_new_opposite_v38.csv`), the calls that a stricter count of BP1's veto removes and the records of the protein index behind them (`bp1_veto_count_resolved_v38.csv`, `bp1_index_counts_v38.csv`), and the binary, commit and annotation time of each run (`provenance_v38.csv`).

**Appendix Table A2.** The offline replay of the combining rules under alternative schemes and cutoffs, its validation against real runs, and the number of calls that the replay cannot reproduce because the site-level frequency veto holds them at VUS (`combining_rules_replay_v38.csv`, `replay_validation_v38.csv`, `replay_veto_v38.csv`, `pm2_dependence_v38.csv`, `full_points_false_benign_signatures_v38.csv`).

**Appendix Table A3.** The aggregates that the text quotes from the per-variant table and that no other table holds (`facts_v38.csv`).

**Appendix Table A4.** The rulings of the three returned workbooks (`review_rounds.csv`), the 64 rows of the last round and what became of them (`round6_rows_v38.csv`), and the rulings with our coding of them (`round6_rulings.tsv`, `round6_rules_v38.csv`).
We coded the reviewer's round-6 text by five rules, which a row can satisfy in more than one way; the column `rule` of `round6_rulings.tsv` records the one that decided it.
(1, class, 29 rows) The text states a class: it is compared with the two calls and coded ClinVar when its direction is ClinVar's, fastVEP when it is fastVEP's, and neither when it is VUS; a pair of classes that includes VUS (one row) was coded by its other member.
(2, stated, 13 rows) The text says who is right, or that both calls are discordant with the reviewer's: coded as stated.
(3, disorder, 2 rows) The text makes the class depend on the disorder or the inheritance, in a gene with both dominant and recessive disease: coded depends on the disorder, whatever class it also states.
(4, review, 6 rows) The text asks for more review, or gives a class only on a condition: coded needs more evidence.
(5, evidence, 14 rows) The text gives evidence for or against a call without a class or a statement of who is right: coded for the call that the evidence supports (13 ClinVar and 1 fastVEP).
Of the 48 rows that are still opposite, 10 of the 29 coded ClinVar come from rule 1, 9 from rule 2 and 10 from rule 5.
The reviewer's text is not reproduced here.

**Appendix Table A5.** What the extraction of the published VCEP thresholds did with each criterion entry it read, for BA1, BS1 and PM2 together and for PM2 alone (`vcep_audit_v38.csv`, `vcep_audit_pm2_v38.csv`).

**Appendix Table A6.** The points scale: the ClinVar classes of the variants at each number of pathogenic points, the posterior probability of pathogenicity that the Bayesian framework assigns to each number of points, and the composition of the pathogenic-side variants called VUS (`points_by_clinvar_v38.csv`, `tavtigian_posterior.csv`, `plp_to_vus_points_v38.csv`, `plp_to_vus_by_consequence_v38.csv`, `plp_to_vus_gene_concentration_v38.csv`, `signatures_offdiagonal_v38.csv`).

**Appendix Table A7.** Gene-level concordance of v38 (`gene_level_v38.csv`, `plp_recall_by_gene_v38.csv`).

**Appendix Table A8.** The PVS1 gap in detail: the genes that gnomAD constrains, the genes that gain PVS1 in each of the three runs that change an input, and the reference-discordant pathogenic calls that those runs add (`pvs1_gap_constrained_v38.csv`, `pvs1_gap_closure_genes_v38.csv`, `pvs1_gap_closure_new_false_pathogenic_v38.csv`).

**Appendix Table A9.** The later ClinVar release in detail: recall by band of allele frequency and standardized for the mix of bands, agreement in direction standardized for the mix of ClinVar classes, the new Benign variants of 2% or more and why BA1 and BS1 did not fire on those that were not called benign-side, the composition of the new variants, their genes and how they stand in the April protein index, the 157 variants that changed direction, and what became of the 48 residual calls (`temporal_recall_by_af_band_v38.csv`, `temporal_standardised_recall_v38.csv`, `temporal_agreement_mix_v38.csv`, `temporal_high_af_benign_v38.csv`, `temporal_high_af_blockers_v38.csv`, `temporal_new_composition_v38.csv`, `temporal_new_top_genes_v38.csv`, `temporal_new_genes_v38.csv`, `temporal_new_in_index_v38.csv`, `temporal_direction_changes_v38.csv`, `temporal_direction_change_variants_v38.csv`, `april_opposite_in_october_v38.csv`).

**Appendix Table A10.** The runs of Section 3.8 and their control: the concordance of the adopted change, of the commit before it and of each alternative, the control of the measurement build against v38, how often each criterion fires, and the facts that the text quotes from the per-variant tables (`v38_runs.csv`, `v38_controls.csv`, `v38_criteria_firing.csv`, `v38_facts.csv`).

**Appendix Table A11.** The reading of the mode of inheritance: the genes of the ClinGen file, and those with a benchmark variant, as dominant-only, recessive-only, both or neither by the words of the disease names (v37) and by the curated field (v38), by the established tags of the gene, with the genes that have no established curation or another mode, and the genes with more than one mode (`v38_inheritance_reading.csv`).

**Appendix Table A12.** What moved: the calls that changed between v38 and the commit before it and between v38 and each alternative, by ClinVar class, the ClinVar pathogenic-side variants among them with the criteria before and after, and the cap on PP2, PM1 and PP3 on the April set and the October release (`v38_transitions.csv`, `v38_clinvar_plp_moves.csv`, `v38_triple.csv`).

**Appendix Table A13.** The calls for each mode of inheritance: how many variants and genes have them, how often they differ and by which criteria, the pairs of calls by ClinVar class, and the genes (`v38_per_disorder.csv`, `v38_per_disorder_pairs.csv`, `v38_per_disorder_genes.csv`).

**Appendix Table A14.** The four thresholds that read inheritance, swept again on v38: each run, and what each step buys and costs (`v38_sweeps.csv`, `v38_threshold_steps.csv`).

**Appendix Table A15.** The October release: the commit before the change, the adopted change and two alternatives on the covered variants and on the new ones, and the ClinVar pathogenic-side variants whose call moved (`v38_holdout.csv`, `v38_holdout_plp_moves.csv`).

**Appendix Table A16.** gnomAD's joint release against its exomes: the three runs on every chromosome and the pilot of the same three on chromosomes 21 and 22, the calls and criteria that moved, the check that the joint release holds the records and counts of the exomes, the check that the pilot's calls are those of the runs of every chromosome, the change on each chromosome, the coverage behind the calls that move, the variants whose standing as an error or as a pathogenic-side call changed, annotated through both databases, and the variants that lose PM2 against the exomes without the flags, by the reason (`v38_joint.csv`, `v38_joint_changes.csv`, `v38_joint_records.csv`, `v38_joint_pilot_check.csv`, `v38_joint_by_chrom.csv`, `v38_joint_coverage.csv`, `v38_joint_errors.csv`, `v38_joint_pm2_lost.csv`).
