# ACMG classification manuscript

A manuscript on how the ACMG/AMP classifier of fastVEP was built, reviewed by a medical geneticist and measured, how its thresholds were chosen, and what its discordance with ClinVar 2-star+ classifications looks like in the final version (v37) and on a later ClinVar release, and what changed after v37, when the classifier began to read the mode of inheritance that ClinGen curates (v38, Section 3.8).
ClinVar's 2-star+ classifications are a reference and not ground truth, throughout.

## Files

| Path | What |
|---|---|
| `manuscript.md` | The manuscript: abstract, methods, results, discussion, references, the legends of the five main figures, and the supplement (four supplementary figures, nine supplementary tables) and appendix (sixteen further tables) |
| `figures/` | Figures 1 to 5 (the main figures, each stating a conclusion) and S1 to S4 (the process figures: the runs, the ledger of changes, the reviewer's rulings and the frontier of agreement against contradictions) as PNG (300 dpi), written by `generate_figures.py` together with a PDF of each (the repository ignores `*pdf`, so only the PNGs are tracked); every label that names a ClinVar class says "ClinVar 2-star+", and each figure has a footnote that gives the review status, the release and the class abbreviations |
| `generate_figures.py` | Draws every figure from the tables in `data/` and nothing else |
| `build_tables.py` | Builds the tables that carry numbers (Tables 3 to 9) from `data/` and keeps `manuscript.md` in step with them |
| `supplementary_tables.xlsx` | Every supplementary table (S1 to S9) and appendix table (A1 to A16) in one workbook: a sheet of contents and a sheet for each table, which holds the CSV or TSV files of `data/` that make it up, one under another, each under its file name |
| `build_supplement.py` | Writes `supplementary_tables.xlsx` from `data/`; which file is in which table is read from the table entries of `manuscript.md`, so the workbook cannot disagree with it |
| `verify_numbers.py` | Recomputes the numbers quoted in the text from `data/` and fails when a sentence does not contain them; it also checks that every figure and table that the text cites exists and is cited in order, that there are fewer than ten supplementary tables, that every file of `data/` is in exactly one table, and that `supplementary_tables.xlsx` holds what those files say |
| `data/` | The supplementary and appendix tables, 97 small CSV and TSV files (each is named by exactly one table of the manuscript, which `verify_numbers.py` checks), written by `analysis/acmg_benchmark/scripts/13_manuscript_stats.py`, `11_manuscript_runs.py`, `14_temporal_holdout.py`, `15_clingen_complete_export.py`, `17_joint_boundary_cases.py` and `18_criteria_status.py` |
| `replay_validation/combiner_cutoff_env.patch` | A patch of eight lines that makes the two point cutoffs of the pathogenic side readable from the environment, for the nine cutoff runs of Table 4 |
| `replay_validation/bp_hooks_env.patch` | A patch that makes BP1's own-record count and BP7's SpliceAI bar readable from the environment, for the runs of Sections 3.3 and 3.4 |
| `replay_validation/inheritance_alternatives_env.patch` | A patch, against the tree of v38, that makes the alternatives of Section 3.8 readable from the environment: BS2's carriers counted beside BS1 and BA1, no carrier arm for a susceptibility or predisposition disorder, BS2 read from the disease names as in v37, and a combined cap on PP2, PM1 and PP3 |

The scored matrix, firing counts and provenance of every run made for the paper are in `analysis/acmg_benchmark/results/v37_manuscript/`, with `binaries.tsv` (what each binary was built from) and `rerun_comparison.csv` (each run and table against the first run it repeats).
The commands that regenerate everything are in `analysis/acmg_benchmark/README.md`.

## Checking the manuscript

```bash
python3 manuscript/acmg/build_tables.py --check    # the tables in the text are what the data say
python3 manuscript/acmg/build_supplement.py --check    # the workbook of supplementary tables is what the data say
python3 manuscript/acmg/verify_numbers.py          # the numbers in the sentences are what the data say
```

`verify_numbers.py` finds a sentence by a phrase that occurs in no other line and checks that every number recomputed from the data is in it.
It does not check numbers that are not a function of the data.
Those are listed under "What is not checked" below.

## What was verified against sources, and by whom

- Every reference and its PMID was checked against PubMed and Crossref by a subagent, which found one year to correct (Brnich et al. is 2019 in PubMed) and could not verify the form of citation for RepeatMasker, whose site was down.
- The statements about what ClinGen and ACMG/AMP recommendations say (Tavtigian points and bands, Pejaver REVEL intervals, Walker splicing rules, the Abou Tayoun PVS1 tree, the Ghosh exception list, the SVI PM2 and PM3 documents, the VCEP Review Committee communication on gnomAD v4, the multi-disorder guidance, Whiffin's filtering allele frequency, Biesecker and Harrison on PP5 and BP6) were checked by a second subagent against the full text of each source.
  The PVS1 result, which changed the manuscript most, was read a second time against the saved full text and Figure 1 of Abou Tayoun et al. (2018).
  The statement about BP7 in Walker et al. (2023) was read a third time in the course of the final revision.
- The manuscript was then audited by three further subagents (its numbers, its claims against the code and the project's documents, and its argument as a peer reviewer would read it), and the audit's findings were checked against the code and the data before they were acted on.
  The audits found, and the manuscript and documents now correct: the omitted ClinGen curations are 104 and not 112 (seven of the 112 rows re-label a curation that the April file lists and one is dated after that file was downloaded), PM3 needs a proband's genotype and not a trio, a stale protein index can raise PM1 as well as lower it, the leave-one-source-out runs do not all report not evaluated (and the gnomAD run also switches a PM2 setting), the environment variable of a patched-build run and the commit of 26 runs were not what the text said, PhyloP feeds only BP7 and SpliceAI also feeds PVS1, one author's name, and two quirks of the label script.
  The peer review also asked for the in-sample figures to be separated from the hold-out, for the baselines to be shown as a curve, and for the changes that the manuscript proposes to be tried on the variants that none of them was tuned on.
  We did: a strict leave-one-out run (59 opposite-direction calls), the baseline curve (Figure S4), a finer allele-frequency standardisation of the benign-side fall, a sweep of BP1's veto count that the audit found had never been made, and runs of the proposed changes on the new variants (Section 3.7).
  Not done: asking the reviewer to confirm our coding of rounds 2 and 6, a comparison with another tool, and the filtering allele frequency in the per-variant tables.
- A fourth subagent then audited only the sentences that the revision had added or changed, against the data and the code.
  Every one of its thirteen main points held when we recomputed it, and the manuscript now corrects them: one of the two missed diagnoses that a stricter count of BP1's veto avoids (PRR12) rests on the variant's own record, and the other (IQSEC1) does at a count of 2 only (a run of the pair, with the own record left out, avoids one); the strict reading adds three Benign calls, so the zero for Benign calls does not survive it; a quiet variant, as the baseline rule defines it, includes the flanking variants and the polypyrimidine tract, which the text left out (without them the rule has 1,036 contradictions and not 1,730); 98.5% of the new variants are in the April protein index with their April record, because they were in ClinVar at a lower review status, so the hold-out of BP1's own-record count is a weak test and the new variants are not variants that ClinVar did not hold; the date skew of the omitted ClinGen curations does not fit a random omission; the baseline family has fewer contradictions at higher bars and without a frequency rule than the text said; only one further source (Orphanet) was measured and not two; and several sentences compared numbers of different kinds.
  Three of the points needed new runs or tables (BP1's count with the own record left out, BP1's count on the new variants, and why BA1 and BS1 did not fire on the high-frequency new Benign variants, which is gnomAD's low-complexity, filter and segmental-duplication flags and the classifier's homologous-gene list for 79% of them).
- A fifth subagent audited those corrections.
  Its five points of meaning held when we recomputed them and are corrected: the removal of IQSEC1 by a stricter count of BP1's veto rests on its own record at a count of 2 and not at a count of 1, the protein index holds the April record only of the new variants that were then classified missense variants (77 of the 94 that BP1 can act on), the rules without a frequency rule agree for 0.3 points more than v37 and not the same, 23 of the 31 new Benign variants that met BS1 and stay at VUS have BS1 alone, and the zero for Benign calls rests on two of the 15 curated alleles and one gate.
- The ClinGen export of 30 April 2026 was compared with an export downloaded on 11 October 2026, which is how the 104 omitted curations were found.
- Sections 2.10 and 3.8 (v38) were written after those audits, and were then audited by two further subagents on 4 November 2026: one read the diff of the code for defects, the other read the manuscript and the documents for logic, consistency and honesty.
  Their findings were checked against the code, the data and the repository's own record before they were acted on.
  They found, and the code, the documents and the manuscript now correct: the library behind both web servers looked an indel up in the databases in a form that never matches (confirmed end to end on 120 real indels), a warning about missing region flags that the command line never shows, BP2's trans arm opened by a gene's recessive disease, the history of the unread tag (v37 had measured a build that read it and declined), the replication on the October release (258 of the 1,762 calls are on new variants), and a doubling that BS2's arm for homozygotes and hemizygotes still has (measured and documented, not changed).
  Every number in them is recomputed from `data/v38_*.csv` by `verify_numbers.py`; the statements about the code were checked against it (the classification crate of the commit before v38 differs from v37's in comments only, by `git diff`), the ones about the sources against the ClinGen exports of 30 April and 11 October 2026 (the three labels with a semicolon, the three X-linked disorders whose label says "dominant", the review's counts of modes of inheritance, which are the October export's), and the reading of the field by a test over 185 genes of real ClinGen rows.
  The runs were made twice, and the binaries of the second set were built from the final tree: the two sets agree in every scored file of every run and in every record of the annotated VCFs that both keep (`analysis/acmg_benchmark/RUN_VERSIONS.md`, "The full rerun of 4 November 2026").
- Not re-verified: the exact strength grades of PS1 for a splice variant (Walker et al., 2023, Table 3) against our implementation, and the claim that Ensembl VEP and fastVEP agree on each annotation difference that the ledger mentions.
- The ledger of 64 changes, its classes and its routes of discovery are our coding of `RUN_VERSIONS.md` and the round documents.
  They have not been coded independently.

## Findings of the verification that are not in the classifier

The manuscript reports these in Sections 3.2 to 3.5, with what each would change, and the repository's comments and documentation were corrected for the citations and for the attributions.
The classifier was not changed for these, because each is a decision that deserves its own review, except for two that Section 3.8 describes (v38): the mode of inheritance is read from the field that ClinGen curates, and NSMCE3's label no longer escapes the validity gate.

- PVS1 grades an NMD-escaping truncation outside a critical region one step below the published tree, and a canonical splice variant in the last exon as Moderate whatever the fraction removed (7 calls, by replay).
- BP7 keeps the 2015 conservation test and uses SpliceAI below 0.2, where Walker et al. apply BP7 after BP4 and caution against a conservation test (the test costs 16,482 correct benign-side calls to avoid 2 missed diagnoses; with SpliceAI at most 0.1 as well, 15,953 more calls in direction for the same 2).
- The allele-number floor of 2,000 comes from the dataset-size condition of the BA1 recommendation and is applied here to the site's allele number and to BS1 (it buys nothing: 774 calls, no error).
- The recessive-gene PM2 bar of 7e-5 was attributed to the SVI PM2 recommendation, which gives no frequency; it equals the Hearing Loss VCEP's recessive bar, and the project's record does not say where it came from (swept: flat).
- BP1's veto counts the variant's own ClinVar record, which PS1 and PM1 exclude (3 missed diagnoses rest on it).
  With that, ClinVar's low-penetrance gate and the 15-allele exception list removed together, there are 59 opposite-direction calls, and three of the 10 added are Benign calls on ClinVar-pathogenic variants (two G6PD and one SERPINA1), so the zero for Benign calls does not survive that reading.
- The ClinGen file of the benchmark lacked 104 curations that ClinGen dates before the day of its export, CFTR among them (with them restored, 674 more truncating variants carry PVS1 and 12 of the 23 four-star records are called pathogenic-side where 1 was).
- PVS1's gene-level gates cannot tell which disease a variant is pathogenic for (RASA1, LRP2 and PTPN11 are examples).
- BP1's veto count of 3 was never swept before this paper, and a stricter count calls fewer ClinVar-VUS variants benign-side: a count of 1 gains 1,716 calls in direction (1,865 ClinVar-VUS variants no longer called benign-side, less 149 correct benign-side calls given up) and removes 2 missed diagnoses (IQSEC1 and PRR12).
  PRR12 owes the veto to the variant's own record, and IQSEC1 does at a count of 2 only (the index holds a second record for it, p.Thr343Met), so with the own record left out a count of 1 removes IQSEC1 only (a real run) and a count of 2 neither.
- PVS1 predicts NMD by the last-exon proxy where Abou Tayoun et al. use the 50-nucleotide rule, a default chosen on agreement with ClinVar (603 calls in direction and 0.73 points of pathogenic-side recall to avoid one false-pathogenic call).
- The validity gate missed NSMCE3 in v37, a gene whose ClinGen curations are all weak, because its label ends in a semicolon (one of 505 weak-only genes, with no benchmark variant); v38 drops an untagged piece of a split label where a tagged one exists.
- The label script maps two RYR1 expert-panel records with a conflicting classification to VUS and four CFH and CFI records classified Likely pathogenic with low penetrance to Pathogenic.
  We kept the labels so that the 37 runs stay comparable; relabelling changes exact agreement by one variant, and leaving out the two RYR1 records lowers the ClinVar VUS called pathogenic-side by two.
- Two comments in the classifier were wrong and are corrected: the one that said the `LOF_and_GOF` mechanisms are inert (they license PVS1 through the curated-mechanism route, which changes no call in the benchmark), and the one that said the ClinGen phenotypes are Moderate or better only (they are every classification, and the validity gate is what withholds the weak ones).

## Decisions that are the author's

These are marked in the manuscript with square brackets or are not in it at all.

1. **The reviewer.**
   The reviewer is M.G., an author since the co-authors' version of 4 November 2026, and W.C., who advised on a subset of the cases, is an author too; the manuscript names them in the author line and in Sections 2.4 and the Declarations.
   The manuscript describes the reviewer's rulings and paraphrases their reasons.
   It quotes only the heading of the list that came with the round-1 workbook ("Discordance categories") and does not reproduce their notes, and `data/round6_rulings.tsv` holds our coding of the text and not the text.
   Whether the text may be quoted or published as a supplementary column is for the authors to settle.
   The workbooks they returned are not in the repository.
2. **Disclosure of AI assistance.**
   The scripts, the figures and the first draft were produced with Claude working as a coding agent in this repository.
   The Declarations section holds a placeholder for the author's statement.
3. **Author contributions** are the text of the draft that the co-authors read, marked for all authors to confirm, and **funding and competing interests** are placeholders.
4. **The title and the abstract** state the main findings and should be read by the author before submission.
   The abstract is about 800 words, the length of a long-form paper, and has to be cut to the limit of the journal chosen; the language-edited version that the co-authors read (a .docx, not in the repository, of v37 only) has a shorter one and another title.
5. **Pronouns.**
   The manuscript refers to the reviewer as "the reviewer" and "they", because their pronouns are not recorded.
6. **Whether to change the classifier** on the findings above, and in what order (Section 4.5).
   Each has a price measured in this work, and none should be tuned to the 49 residual calls.
7. **The commit of v38.**
   The runs of Section 3.8 were made by a binary built from the working tree that became v38, so their provenance records the commit of the checkout that ran the harness (the commit before the change, and for `v38_joint`, made on 5 November 2026 after the change was committed on a local branch, that commit) and a count of uncommitted files, and `analysis/acmg_benchmark/results/v37_manuscript/binaries.tsv` says what each binary was built from.
   The manuscript marks the place for the commit in square brackets, and the author should fill it in when the change is merged.
8. **The decisions that v38 took** (Section 3.8).
   v37 had measured a build that made BS2 read the tag (`v37_adtag`) and left the tag unread, because six ClinVar pathogenic-side variants moved toward benign; v38 reads it, for the reason given in Section 3.8 (the unread field answers a question about the name of a disease, and the guard is stated by what it counts), and accepts that four of the six move from Likely pathogenic to VUS.
   The guard against counting BS2 beside BS1 or BA1 was adopted at a cost of 3,615 exact calls on the April set.
   A combined cap on PP2, PM1 and PP3 was not adopted (it raises agreement and lowers recall of pathogenic-side variants).
   The label guard that restores two of the four was not adopted although its exchange rate, about 32 correct benign-side calls for each variant restored, is below the knee of Section 2.5, because it filters on the name of a disease and was written to repair two variants seen in the results.
   BS2's arm for homozygotes and hemizygotes was left as it was, with the same doubling on the X chromosome for 1,301 variants none of which is a ClinVar pathogenic-side variant.
   All are reasoned in the text and measured on the April set and the October release, and each can be reversed with the numbers of Table 8 and Appendix Tables A10 and A12.
9. **The gnomAD source of the benchmark.**
   Section 3.8 measures gnomAD's joint release of exomes and genomes on every chromosome: 10,190 more exact calls and 9,321 more in direction than the exomes, for one more missed diagnosis and two more false-pathogenic calls, and without the region flags that withhold the frequency criteria.
   `docs/ACMG_SETUP.md` recommends it for classification where no policy requires the region veto, and the benchmark of the paper still reads the exomes, because every other number in the paper was measured on them.
   Moving the benchmark to the joint release would change all of those numbers, so it is left to the author.
10. **Archiving.**
   The bulk outputs (annotated VCFs, per-variant tables of about 20 MB each) are not in the repository.
   A tagged release with a DOI, and a compressed deposit of the per-variant tables, would let every table in `data/` be re-derived without a re-run.

## What is not checked

`verify_numbers.py` guards the numbers, not the wording.
A sentence can pass and still say something the data do not support.
It does not check the following, which are not functions of `data/`.

- Dates, version numbers and constants of the code (the 2,000-allele floor, the pLI and LOEUF bars, the missense z-score bars, the 15 alleles of the exception list, the 26 criteria that a run reports; `18_criteria_status.py` checks the criteria that are not applied to a sites-only input against the engine's JSON on a sample of the benchmark, which is not the whole benchmark, and stops when the engine disagrees with the groups it states).
- The counts read from the project's own record: the split of the 2.7 points that v9 lost (0.6 and 2.0), the 10,000-row cap and the chromosome at which it stopped the first export, the 3,793 Likely pathogenic calls that PS1's splice path turned into Pathogenic, the 14 deletions that lost their PVS1, and the 19-point knee that v10 recorded.
- What each recommendation says (Section 3.2 and Table 1), which was checked by subagents and by reading, and the release dates of ClinVar and ClinGen inputs.
- The reviewer's remarks, which are paraphrased in Section 3.6 and coded, not stored.
- The statements of Sections 2.10 and 3.8 that are facts of the code or of a source file and not of `data/`: that the classification code of the commit before v38 differs from v37's in comments only (`git diff 5d8004a 9d84be0 -- crates/fastvep-classification`), that the three labels with a semicolon and the three X-linked disorders whose label says "dominant" (GJB1, PDK3 and PHEX) are those of the exports of 30 April and 11 October 2026, that the joint release declares no `lcr`, `segdup` or `non_par` flag in the headers read on 4 November 2026, that its FILTER column says only which side failed (the definitions of `EXOMES_FILTERED`, `GENOMES_FILTERED`, `BOTH_FILTERED` and `PASS` in its header, read on 5 November 2026) while the builder reports a filter that fired in either callset (`crates/fastvep-sa/src/sources/gnomad.rs`), and that the pseudoautosomal coordinates agree with the flag in gnomAD's genomes file at the edges that could be read.
- What the record of v37 says about the build that made BS2 read the tag (`v37_adtag`: 2,812 more calls in direction and 3,809 more exact, six ClinVar pathogenic-side variants moved toward benign), which is in `docs/ACMG.md` (limitation 20) and `RUN_VERSIONS.md` and was recomputed from its stored matrix when Section 3.8 was written.
- The end-to-end comparison of the library with the command line on 120 ClinVar indels of chromosomes 21 and 22 (gnomAD found for none and for 67 before the change, and the two agree on all 120 after it), which was a one-off with the benchmark's databases and is recorded in `RUN_VERSIONS.md`.
- That a build of the final tree reproduces the annotated VCF of the `v38` run byte for byte, after the two speed-ups of 5 November 2026 (a one-off `cmp` against the stored VCF, recorded in `RUN_VERSIONS.md`, "After v38").
- That no multi-mode gene of the ClinGen exports of 30 April and 11 October 2026 has an X-linked, mitochondrial or undetermined curation among its established ones: their modes are AD with AR, AD with SD, AR with SD, and all three.

Other limits:

- The hold-out of Section 3.7 is a hold-out of variants and not of genes, and the ClinVar-derived indices of the classifier were not rebuilt for it.
  Its new variants are mostly variants that ClinVar held in April at a lower review status (98.5% have a record in the `variant_summary` that the protein index was built from), so they are variants that no sweep or review scored, and the protein index holds the April record of those that were then classified missense variants.
- The cutoff and measurement builds are patched builds of the shipped code, and the replay of the benign-side point scales has no real run to be checked against.
- The bulk outputs of the runs (annotated VCFs, per-variant tables) are under `data/benchmark/`, which is not in the repository, so a reader cannot re-derive `data/` without rerunning the benchmark.
