# ACMG classification manuscript

A manuscript on how the ACMG/AMP classifier of fastVEP was built, reviewed by a medical geneticist and measured, how its thresholds were chosen, and what its discordance with ClinVar 2-star+ classifications looks like in the latest version (v38) and on a later ClinVar release, and what v38 changed from v37, the version in which the review rounds ended, when the classifier began to read the mode of inheritance that ClinGen curates (Section 3.8).
v37 and v38 are versions of the classifier, and every run of either is on GRCh38.
Every number of the manuscript is v38's unless the text names v37, which it does for the history of the numbered versions and as the before of a difference.
ClinVar's 2-star+ classifications are a reference and not ground truth, throughout.

## Vocabulary

The manuscript, its figures and its tables say "default", "threshold", "reference-discordant benign call", "reference-discordant pathogenic call" and BP7's "intronic distance", as the language-edited draft that the co-authors read does.
They no longer say "shipped", "bar", "missed diagnosis", "false-benign" or "false-pathogenic call", or BP7's "far boundary".
The two discordant classes are named for what was measured, a call that disagrees with a ClinVar 2-star+ class, and not for a verdict that the call is wrong or that a diagnosis was missed.
"Default" is the setting of the version of the classifier that the sentence is about (v38, unless the text names v37).
"Default scheme" is the combining scheme of the default, which is the same in v37 and v38: the point scale of Tavtigian et al. (2020) on the pathogenic side and the 2015 combining table on the benign side.
"The point scale", or "points" with the side, is its pathogenic half, and "the 2015 table on the benign side" is its other half.
"The 2015 table on both sides", or "under the 2015 table" of a threshold's sweep, is the setting `use_point_system = false`, which was run.
"Points on both sides" has no setting and exists only as an offline replay of the recorded criteria (Section 2.7).
The text does not call the default "the point system", which reads as points on both sides: "the point system" is the scheme of Tavtigian et al. itself, and "the point scale" is the pathogenic half of the default.
The older words remain only where renaming them would break something that reads them, or would change what the classifier prints: the columns `false_benign` and `false_pathogenic` and the words `shipped` and `bar` in the files of `data/`, the setting names in their `knob` column (for example "BS1 bar" and "BP7 far boundary"), the run names (`v38_oct_shipped`), the scripts and the run descriptions that they write into `analysis/acmg_benchmark/results/`, and the classifier's own output text (the evidence summaries of BS2, BP7 and PM2 and one configuration error, which say "prevalence bar", "far boundary" and "rarity bar").
Read `false_benign` as the reference-discordant benign calls of the text, and `false_pathogenic` as its reference-discordant pathogenic calls.

## Files

| Path | What |
|---|---|
| `manuscript.md` | The manuscript: abstract, methods, results, discussion, references, the legends of the five main figures, and the supplement (four supplementary figures, nine supplementary tables) and appendix (sixteen further tables) |
| `figures/` | Figures 1 to 5 (the main figures, each stating a conclusion) and S1 to S4 (the process figures: the runs, the ledger of changes, the reviewer's rulings and the frontier of agreement against contradictions) as PNG (300 dpi), written by `generate_figures.py` together with a PDF of each (the repository ignores `*pdf`, so only the PNGs are tracked); every label that names a ClinVar class says "ClinVar 2-star+", and each figure has a footnote that gives the review status, the release and the class abbreviations |
| `generate_figures.py` | Draws every figure from the tables in `data/` and nothing else; `--base v37 --out DIR` draws the figures of v37 into a directory of your own |
| `build_tables.py` | Builds the tables that carry numbers (Tables 3 to 9) from `data/` and keeps `manuscript.md` in step with them; `--base v37` checks them against the tables of v37 (the stale ones are what v38 changed) |
| `supplementary_tables.xlsx` | Every supplementary table (S1 to S9) and appendix table (A1 to A16) in one workbook: a sheet of contents and a sheet for each table, which holds the CSV or TSV files of `data/` that make it up, one under another, each under its file name |
| `build_supplement.py` | Writes `supplementary_tables.xlsx` from `data/`; which file is in which table is read from the table entries of `manuscript.md`, so the workbook cannot disagree with it; `--base v37 --write --out PATH` writes the workbook of v37 to a file of your own |
| `verify_numbers.py` | Recomputes the numbers quoted in the text from `data/` (the files of v38 by default) and fails when a sentence does not contain them, a number being found only where it is a number of its own (not inside a longer number, a name such as `PM2` or `v38`, a hash or a power such as `4e-5`); `--base v37` recomputes them from the files of v37, and its failures are the sentences that v38 changed; it also checks that every figure and table that the text cites exists and is cited in order, that there are fewer than ten supplementary tables, that every file of `data/` is in exactly one table, and that `supplementary_tables.xlsx` holds what those files say |
| `data/` | The supplementary and appendix tables, 164 small CSV and TSV files (each is named by exactly one table of the manuscript, which `verify_numbers.py` checks): for every table that measures the final classifier a file for v38, which the manuscript reads, and a twin for v37 (`_v38` and `_v37` in the names, 67 pairs), the 23 files of Section 3.8 (`v38_` in the names), and seven that belong to neither version (the history of the numbered versions, the reviewer's ruling and three inputs); written by `analysis/acmg_benchmark/scripts/13_manuscript_stats.py`, `11_manuscript_runs.py`, `14_temporal_holdout.py`, `15_clingen_complete_export.py`, `17_joint_boundary_cases.py` and `18_criteria_status.py` |
| `replay_validation/combiner_cutoff_env.patch` | A patch of eight lines that makes the two point cutoffs of the pathogenic side readable from the environment, for the nine cutoff runs of Table 4 |
| `replay_validation/bp_hooks_env.patch` | A patch that makes BP1's own-record count and BP7's SpliceAI threshold readable from the environment, for the runs of Sections 3.3 and 3.4 |
| `replay_validation/inheritance_alternatives_env.patch` | A patch, against the tree of v38, that makes the alternatives of Section 3.8 readable from the environment: BS2's carriers counted beside BS1 and BA1, no carrier arm for a susceptibility or predisposition disorder, BS2 read from the disease names as in v37, and a combined cap on PP2, PM1 and PP3 |

The scored matrix, firing counts and provenance of every run made for the paper (the 96 runs of v38 and the same 96 of v37) are in `analysis/acmg_benchmark/results/v37_manuscript/` (a directory that keeps the name it had when it held the runs of v37 alone), with `binaries_v38.tsv` and `binaries.tsv` (what each binary of v38 and of v37 was built from) and `rerun_comparison.csv` of `results/v38_v0.5.0/` (each run and table against the first run it repeats).
`analysis/acmg_benchmark/results/v38_family_vs_v37.md` sets the two versions against each other run for run, and says in its Section 9 which statements and numbers of the manuscript would read differently under v37.
The commands that regenerate everything are in `analysis/acmg_benchmark/README.md`.

## Checking the manuscript

```bash
python3 manuscript/acmg/build_tables.py --check    # the tables in the text are what the data say
python3 manuscript/acmg/build_supplement.py --check    # the workbook of supplementary tables is what the data say
python3 manuscript/acmg/verify_numbers.py          # the numbers in the sentences are what the data say
python3 manuscript/acmg/verify_numbers.py --base v37    # not a check of the text: the sentences whose numbers v38 changed from v37
```

`verify_numbers.py` finds a sentence by a phrase that occurs in no other line and checks that every number recomputed from the data is in it.
It does not check numbers that are not a function of the data.
Those are listed under "What is not checked" below.

## What was verified against sources, and by whom

The audits below are those of the manuscript as it was written for v37, and the numbers that they quote are v37's unless the text says v38.
The manuscript was then re-based on v38, the latest version, which `verify_numbers.py` and `20_compare_families.py` check (see "The re-basing on v38" below).

- Every reference and its PMID was checked against PubMed and Crossref by a subagent, which found one year to correct (Brnich et al. is 2019 in PubMed) and could not verify the form of citation for RepeatMasker, whose site was down.
- The statements about what ClinGen and ACMG/AMP recommendations say (Tavtigian points and bands, Pejaver REVEL intervals, Walker splicing rules, the Abou Tayoun PVS1 tree, the Ghosh exception list, the SVI PM2 and PM3 documents, the VCEP Review Committee communication on gnomAD v4, the multi-disorder guidance, Whiffin's filtering allele frequency, Biesecker and Harrison on PP5 and BP6) were checked by a second subagent against the full text of each source.
  The PVS1 result, which changed the manuscript most, was read a second time against the saved full text and Figure 1 of Abou Tayoun et al. (2018).
  The statement about BP7 in Walker et al. (2023) was read a third time in the course of the final revision.
- The manuscript was then audited by three further subagents (its numbers, its claims against the code and the project's documents, and its argument as a peer reviewer would read it), and the audit's findings were checked against the code and the data before they were acted on.
  The audits found, and the manuscript and documents now correct: the omitted ClinGen curations are 104 and not 112 (seven of the 112 rows re-label a curation that the April file lists and one is dated after that file was downloaded), PM3 needs a proband's genotype and not a trio, a stale protein index can raise PM1 as well as lower it, the leave-one-source-out runs do not all report not evaluated (and the gnomAD run also switches a PM2 setting), the environment variable of a patched-build run and the commit of 26 runs were not what the text said, PhyloP feeds only BP7 and SpliceAI also feeds PVS1, one author's name, and two quirks of the label script.
  The peer review also asked for the in-sample figures to be separated from the hold-out, for the baselines to be shown as a curve, and for the changes that the manuscript proposes to be tried on the variants that none of them was tuned on.
  We did: a strict leave-one-out run (59 opposite-direction calls under v37, 58 under v38), the baseline curve (Figure S4), a finer allele-frequency standardisation of the benign-side fall, a sweep of BP1's veto count that the audit found had never been made, and runs of the proposed changes on the new variants (Section 3.7).
  Not done: asking the reviewer to confirm our coding of rounds 2 and 6, a comparison with another tool, and the filtering allele frequency in the per-variant tables.
- A fourth subagent then audited only the sentences that the revision had added or changed, against the data and the code.
  Every one of its thirteen main points held when we recomputed it, and the manuscript now corrects them: one of the two reference-discordant benign calls that a stricter count of BP1's veto avoids (PRR12) rests on the variant's own record, and the other (IQSEC1) does at a count of 2 only (a run of the pair, with the own record left out, avoids one); the strict reading adds three Benign calls, so the zero for Benign calls does not survive it; a quiet variant, as the baseline rule defines it, includes the flanking variants and the polypyrimidine tract, which the text left out (without them the rule had, under v37, 1,036 contradictions and not 1,730); 98.5% of the new variants are in the April protein index with their April record, because they were in ClinVar at a lower review status, so the hold-out of BP1's own-record count is a weak test and the new variants are not variants that ClinVar did not hold; the date skew of the omitted ClinGen curations does not fit a random omission; the baseline family has fewer contradictions at higher thresholds and without a frequency rule than the text said; only one further source (Orphanet) was measured and not two; and several sentences compared numbers of different kinds.
  Three of the points needed new runs or tables (BP1's count with the own record left out, BP1's count on the new variants, and why BA1 and BS1 did not fire on the high-frequency new Benign variants, which is gnomAD's low-complexity, filter and segmental-duplication flags and the classifier's homologous-gene list for 79% of them).
- A fifth subagent audited those corrections.
  Its five points of meaning held when we recomputed them and are corrected: the removal of IQSEC1 by a stricter count of BP1's veto rests on its own record at a count of 2 and not at a count of 1, the protein index holds the April record only of the new variants that were then classified missense variants (77 of the 94 that BP1 can act on), the rules without a frequency rule agree for 0.3 points more than v37 (0.1 more than v38) and not the same, 23 of the 31 new Benign variants that met BS1 and stay at VUS have BS1 alone, and the zero for Benign calls rests on two of the 15 curated alleles and one gate.
- The ClinGen export of 30 April 2026 was compared with an export downloaded on 11 October 2026, which is how the 104 omitted curations were found.
- Sections 2.10 and 3.8 (v38) were written after those audits, and were then audited by two further subagents on 4 November 2026: one read the diff of the code for defects, the other read the manuscript and the documents for logic, consistency and honesty.
  Their findings were checked against the code, the data and the repository's own record before they were acted on.
  They found, and the code, the documents and the manuscript now correct: the library behind both web servers looked an indel up in the databases in a form that never matches (confirmed end to end on 120 real indels), a warning about missing region flags that the command line never shows, BP2's trans arm opened by a gene's recessive disease, the history of the unread tag (v37 had measured a build that read it and declined), the replication on the October release (258 of the 1,762 calls are on new variants), and a doubling that BS2's arm for homozygotes and hemizygotes still has (measured and documented, not changed).
  Every number in them is recomputed from `data/v38_*.csv` by `verify_numbers.py`; the statements about the code were checked against it (the classification crate of the commit before v38 differs from v37's in comments only, by `git diff`), the ones about the sources against the ClinGen exports of 30 April and 11 October 2026 (the three labels with a semicolon, the three X-linked disorders whose label says "dominant", the review's counts of modes of inheritance, which are the October export's), and the reading of the field by a test over 185 genes of real ClinGen rows.
  The runs were made twice, and the binaries of the second set were built from the final tree: the two sets agree in every scored file of every run and in every record of the annotated VCFs that both keep (`analysis/acmg_benchmark/RUN_VERSIONS.md`, "The full rerun of 4 November 2026").
- The re-basing on v38 (below) was audited on 8 November 2026 by two further subagents: one read the manuscript sentence by sentence against the data, and the other read the diff of the guard and of the scripts and tried to break the guard by changing each number of a checked sentence (a mutation test).
  Their findings were recomputed before they were acted on.
  They found, and the manuscript, the guard and the scripts now correct: a count of truncating variants at full strength, the denominator of the benign-side calls of ClinVar's VUS variants (the guard held it as a typed 4,616, which is v37's), a span of Table 3 that left the default out where the text did not, a sentence about the added reference-discordant pathogenic calls of one run that still described two variants where v38 has one, a leak of the validity gate that the text assigned to v38 and that v38 closes, a statement that a sweep's knee is not reproduced where its rise is not, the description of release 0.5.0 (it holds the work on speed that followed v38, which the repetition of 42 runs shows to change no call), and several statements of how the runs were made.
  In the guard, a number was found inside a name, a hash, a power or a percentage (so the 1 of `PVS1` passed for a count spelled "one"), and `has()` now refuses that; a collection step of `11_manuscript_runs.py` would have written the two run tables of v37 with a header and no row on a machine that holds only the family v38, which is the family it now runs by default; and the comparison document dated the runs of the family v38 in UTC where the paper dates in New York.
- Not re-verified: the exact strength grades of PS1 for a splice variant (Walker et al., 2023, Table 3) against our implementation, and the claim that Ensembl VEP and fastVEP agree on each annotation difference that the ledger mentions.
- The ledger of 64 changes, its classes and its routes of discovery are our coding of `RUN_VERSIONS.md` and the round documents.
  They have not been coded independently.

## The re-basing on v38

The manuscript was written for v37, the version in which the review rounds ended, and a co-author read its numbers as v37's.
It now says v38, the latest classifier (release 0.5.0), everywhere that it quotes a result, and v37 only as the history of the numbered versions and as the before of a difference (Section 3.8).
v37 and v38 are versions of the classifier and not builds of the genome: every run of either is on GRCh38.
The default scheme is the same in both: points on the pathogenic side and the 2015 table on the benign side.

- What v38 adds to v37 is the mode of inheritance read from the field that ClinGen curates (and not from the words of the disease names), BS2's carrier count set aside beside BS1 or BA1, and a gene classified once for each of its modes of inheritance.
- Every table, figure and number is v38's.
  The twin of each table for v37 is in `data/` (`_v37` in the name), and `build_tables.py`, `build_supplement.py`, `generate_figures.py` and `verify_numbers.py` take `--base v37` for it.
- Section 3.8 is the account of what v38 changed, from `analysis/acmg_benchmark/results/v38_family_vs_v37.md`, which sets the 96 runs of v38 against the same 96 of v37, and its Section 9 lists what would read differently under v37.
- The benign-side recall of Sections 3.1 to 3.7 is over the 284,550 ClinVar benign-side variants, with the 382 records that have no alternative allele counted as no-calls, as in v37's tables.
  The tables of Section 3.8 are over the 284,168 that have a call, because the binaries of v38 write no record for those 382, and Table 8's caption says so.
- Not re-measured with the complete databases: the time of annotation, which the manuscript quotes for release 0.5.0, measured alone on 7 November 2026 with the first databases (`results/v38_v0.5.0/annotation_time.tsv`), and for the runs of v37 made on 8 November 2026 (a median of the provenance of the 89 runs, up to seven at a time with other work on the machine).
  The runs of the paper are not timings, and `scripts/19_time_annotation.py` now reads the complete stack, so that a new measurement is one.
- Not verified, and said so where it is stated: the reading that the firing counts support for why the PM2 sweep under the 2015 table flattens on v38 (no run sweeps both PM2 thresholds under the 2015 table), and which of the curated tag and the disease names accounts for the larger loss of benign-side calls when the ClinGen file is removed (no per-variant table was made for that run).
- The language-edited .docx that the co-authors read still says v37 throughout.
  It has to be brought to the text of `manuscript.md` by the corresponding author, and nothing in the repository edits it.

## Findings of the verification that are not in the classifier

The manuscript reports these in Sections 3.2 to 3.5, with what each would change, and the repository's comments and documentation were corrected for the citations and for the attributions.
The classifier was not changed for these, because each is a decision that deserves its own review, except for two that Section 3.8 describes (v38): the mode of inheritance is read from the field that ClinGen curates, and NSMCE3's label no longer escapes the validity gate.

- PVS1 grades an NMD-escaping truncation outside a critical region one step below the published tree, and a canonical splice variant in the last exon as Moderate whatever the fraction removed (2 calls on v38 and 7 on v37, by replay).
- BP7 keeps the 2015 conservation test and uses SpliceAI below 0.2, where Walker et al. apply BP7 after BP4 and caution against a conservation test (the test costs 16,116 correct benign-side calls to avoid 2 reference-discordant benign calls; with SpliceAI at most 0.1 as well, 15,552 more calls in direction for the same 2).
- The allele-number floor of 2,000 comes from the dataset-size condition of the BA1 recommendation and is applied here to the site's allele number and to BS1 (it buys nothing: 728 calls, no error).
- The recessive-gene PM2 threshold of 7e-5 was attributed to the SVI PM2 recommendation, which gives no frequency; it equals the Hearing Loss VCEP's recessive threshold, and the project's record does not say where it came from (swept: flat).
- BP1's veto counts the variant's own ClinVar record, which PS1 and PM1 exclude (3 reference-discordant benign calls rest on it).
  With that, ClinVar's low-penetrance gate and the 15-allele exception list removed together, there are 58 opposite-direction calls, and three of the 10 added are Benign calls on ClinVar-pathogenic variants (two G6PD and one SERPINA1), so the zero for Benign calls does not survive that reading.
- The ClinGen file of the benchmark lacked 104 curations that ClinGen dates before the day of its export, CFTR among them (with them restored, 674 more truncating variants carry PVS1 and 12 of the 23 four-star records are called pathogenic-side where 1 was).
- PVS1's gene-level gates cannot tell which disease a variant is pathogenic for (RASA1, LRP2 and PTPN11 are examples).
- BP1's veto count of 3 was never swept before this paper, and a stricter count calls fewer ClinVar-VUS variants benign-side: a count of 1 gains 1,735 calls in direction (1,862 ClinVar-VUS variants no longer called benign-side, less 127 correct benign-side calls given up) and removes 2 reference-discordant benign calls (IQSEC1 and PRR12).
  PRR12 owes the veto to the variant's own record, and IQSEC1 does at a count of 2 only (the index holds a second record for it, p.Thr343Met), so with the own record left out a count of 1 removes IQSEC1 only (a real run) and a count of 2 neither.
- PVS1 predicts NMD by the last-exon proxy where Abou Tayoun et al. use the 50-nucleotide rule, a default chosen on agreement with ClinVar (603 calls in direction and 0.73 points of pathogenic-side recall to avoid one reference-discordant pathogenic call).
- The validity gate missed NSMCE3 in v37, a gene whose ClinGen curations are all weak, because its label ends in a semicolon (one of 505 weak-only genes, with no benchmark variant); v38 drops an untagged piece of a split label where a tagged one exists.
- The label script maps two RYR1 expert-panel records with a conflicting classification to VUS and four CFH and CFI records classified Likely pathogenic with low penetrance to Pathogenic.
  We kept the labels so that the 38 runs stay comparable; relabelling changes exact agreement by one variant, and leaving out the two RYR1 records lowers the ClinVar VUS called pathogenic-side by two.
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
   The abstract is about 1,215 words (it was about 860 before the sentences that set v38 against v37 and the statement of the contribution were added), the length of a long-form paper, and has to be cut to the limit of the journal chosen; the Background and the Conclusions are the part that states the contribution, and the Methods and the Results hold the detail that can go; the language-edited version that the co-authors read (a .docx, not in the repository, which says v37 throughout) has a shorter one and another title.
   The title was not changed: it describes the process (incorporation, review, thresholds, discordance) and does not state the contribution, and two titles that do were proposed in the hand-off for the corresponding author to choose from or to reject.
5. **"Community opinion".**
   The abstract, the introduction, the discussion and the conclusions say that the contribution is how to implement underspecified and fragmented criteria so that the calls agree closely with community opinion.
   Section 2.2 defines the term as ClinVar's classifications of review status two stars or better and one reviewer's reading, and says that it is not a poll of the field.
   Whether the term is acceptable for an expert component of one reviewer who saw both calls is for the authors to settle, and so is whether the sentence that says what the paper recommends (the opening of Section 4) belongs in the Discussion or in the Introduction.
6. **Pronouns.**
   The manuscript refers to the reviewer as "the reviewer" and "they", because their pronouns are not recorded.
7. **Whether to change the classifier** on the findings above, and in what order (Section 4.5).
   Each has a price measured in this work, and none should be tuned to the 48 residual calls.
8. **The commit of v38.**
   The runs of Section 3.8 were made by a binary built from the working tree that became v38, so their provenance records the commit of the checkout that ran the harness (the commit before the change, and for `v38_joint`, made on 5 November 2026 after the change was committed on a local branch, that commit) and a count of uncommitted files, and `analysis/acmg_benchmark/results/v37_manuscript/binaries.tsv` says what each binary was built from.
   The change is the commit `5293cc0`, merged in pull request #167 (`6504819`), and the manuscript names both; binaries built from the tree of release 0.5.0 (the tag `v0.5.0-rc.1`) made the 42 runs again on 7 November 2026 and reproduce them (`analysis/acmg_benchmark/RUN_VERSIONS.md`, "Release 0.5.0 reproduces v38").
9. **The decisions that v38 took** (Section 3.8).
   v37 had measured a build that made BS2 read the tag (`v37_adtag`) and left the tag unread, because six ClinVar pathogenic-side variants moved toward benign; v38 reads it, for the reason given in Section 3.8 (the unread field answers a question about the name of a disease, and the guard is stated by what it counts), and accepts that four of the six move from Likely pathogenic to VUS.
   The guard against counting BS2 beside BS1 or BA1 was adopted at a cost of 3,615 exact calls on the April set.
   A combined cap on PP2, PM1 and PP3 was not adopted (it raises agreement and lowers recall of pathogenic-side variants).
   The label guard that restores two of the four was not adopted although its exchange rate, about 32 correct benign-side calls for each variant restored, is below the knee of Section 2.5, because it filters on the name of a disease and was written to repair two variants seen in the results.
   BS2's arm for homozygotes and hemizygotes was left as it was, with the same doubling on the X chromosome for 1,301 variants none of which is a ClinVar pathogenic-side variant.
   All are reasoned in the text and measured on the April set and the October release, and each can be reversed with the numbers of Table 8 and Appendix Tables A10 and A12.
10. **The tag that the manuscript names.**
   The manuscript names the first release candidate of release 0.5.0 (`v0.5.0-rc.1`, `c12d9f9`) as the code of v38 and says that "release 0.5.0" means that candidate, because no other tag of it exists yet.
   When v0.5.0 is tagged, `verify_numbers.py` will still read the candidate.
   Whether the final tag changes any call has to be shown by the repetition of the runs (`scripts/16_compare_runs.py`, as for the 42 runs of `results/v38_v0.5.0/`), and the text has to name the final tag if the author wants it named.
11. **The gnomAD source of the benchmark.**
   Section 3.8 measures gnomAD's joint release of exomes and genomes on every chromosome: 10,193 more exact calls and 9,187 more in direction than the exomes, for no more reference-discordant benign call and two more reference-discordant pathogenic calls, and without the region flags that withhold the frequency criteria.
   `docs/ACMG_SETUP.md` recommends it for classification where no policy requires the region veto, and the benchmark of the paper still reads the exomes, because every other number in the paper was measured on them.
   Moving the benchmark to the joint release would change all of those numbers, so it is left to the author.
12. **Archiving.**
   The bulk outputs (annotated VCFs, per-variant tables of about 20 MB each) are not in the repository.
   A tagged release with a DOI, and a compressed deposit of the per-variant tables, would let every table in `data/` be re-derived without a re-run.
13. **The SpliceAI and PhyloP databases of the benchmark lacked records that the gnomAD file carries, and the benchmark was re-based on complete ones (8 November 2026).**
   The first databases were distilled in April 2026 from regions of the gnomAD v4.1 exomes file and lacked 8.2% (SpliceAI) and 8.3% (PhyloP) of the values that exist for the benchmark's variants, in 14 batches of the extraction that end early and, for PhyloP, in 457 indels (`analysis/acmg_benchmark/results/complete_databases/README.md`).
   Every run was made again with databases built from the whole file, the October release is scored whole, and the manuscript, its tables, workbook and figures were regenerated: Section 2.1 says what happened, Section 3.9 and Table 10 set the first databases against the complete ones, and the first drafts, including the version that the co-authors read, are the numbers of the first databases.
   v38 on the April set goes from 77.84% to 79.32% of agreement in direction and from 48 to 51 opposite-direction calls; the three added are ClinVar-Pathogenic variants (TMEM216, VWF, GTPBP3) that the reviewer has not seen.
   The decisions that are the authors': whether the history of the first drafts stays in the main text (Sections 2.1, 3.9 and 4.4) or goes to the supplement; whether the reviewer rules on the three added calls; and whether v1 to v36, which were not repeated, are rebuilt (their binaries are not kept).

## What is not checked

`verify_numbers.py` guards the numbers, not the wording.
A sentence can pass and still say something the data do not support.
It does not check the following, which are not functions of `data/`.

- Dates, version numbers and constants of the code (the 2,000-allele floor, the pLI and LOEUF thresholds, the missense z-score thresholds, the 15 alleles of the exception list, the 26 criteria that a run reports; `18_criteria_status.py` checks the criteria that are not applied to a sites-only input against the engine's JSON on a sample of the benchmark, which is not the whole benchmark, and stops when the engine disagrees with the groups it states).
- The counts read from the project's own record: the split of the 2.7 points that v9 lost (0.6 and 2.0), the 10,000-row cap and the chromosome at which it stopped the first export, the 3,793 Likely pathogenic calls that PS1's splice path turned into Pathogenic, the 14 deletions that lost their PVS1, and the 19-point knee that v10 recorded.
- What each recommendation says (Section 3.2 and Table 1), which was checked by subagents and by reading, and the release dates of ClinVar and ClinGen inputs.
- The reviewer's remarks, which are paraphrased in Section 3.6 and coded, not stored; the remark on how common each of the five PVS1 variants is in gnomAD ("too common" for KDM2B and PHKB, more than 1,000 carriers for EDN3 and PIEZO1, the hemizygote count for FMR1) was read from the returned workbook of round 6 when the manuscript was re-based on v38, and `verify_numbers.py` does not check it.
- The statements of Sections 2.10 and 3.8 that are facts of the code or of a source file and not of `data/`: that the classification code of the commit before v38 differs from v37's in comments only (`git diff 5d8004a 9d84be0 -- crates/fastvep-classification`), that the three labels with a semicolon and the three X-linked disorders whose label says "dominant" (GJB1, PDK3 and PHEX) are those of the exports of 30 April and 11 October 2026, that the joint release declares no `lcr`, `segdup` or `non_par` flag in the headers read on 4 November 2026, that its FILTER column says only which side failed (the definitions of `EXOMES_FILTERED`, `GENOMES_FILTERED`, `BOTH_FILTERED` and `PASS` in its header, read on 5 November 2026) while the builder reports a filter that fired in either callset (`crates/fastvep-sa/src/sources/gnomad.rs`), and that the pseudoautosomal coordinates agree with the flag in gnomAD's genomes file at the edges that could be read.
- What the record of v37 says about the build that made BS2 read the tag (`v37_adtag`: 2,812 more calls in direction and 3,809 more exact, six ClinVar pathogenic-side variants moved toward benign), which is in `docs/ACMG.md` (limitation 20) and `RUN_VERSIONS.md` and was recomputed from its stored matrix when Section 3.8 was written.
- The end-to-end comparison of the library with the command line on 120 ClinVar indels of chromosomes 21 and 22 (gnomAD found for none and for 67 before the change, and the two agree on all 120 after it), which was a one-off with the benchmark's databases and is recorded in `RUN_VERSIONS.md`.
- That a build of the final tree reproduces the annotated VCF of the `v38` run byte for byte, after the two speed-ups of 5 November 2026 (a one-off `cmp` against the stored VCF, recorded in `RUN_VERSIONS.md`, "After v38").
- That no multi-mode gene of the ClinGen exports of 30 April and 11 October 2026 has an X-linked, mitochondrial or undetermined curation among its established ones: their modes are AD with AR, AD with SD, AR with SD, and all three.

Other limits:

- The hold-out of Section 3.7 is a hold-out of variants and not of genes, and the ClinVar-derived indices of the classifier were not rebuilt for it.
  Its new variants are mostly variants that ClinVar held in April at a lower review status (98.5% have a record in the `variant_summary` that the protein index was built from), so they are variants that no sweep or review scored, and the protein index holds the April record of those that were then classified missense variants.
- The cutoff and measurement builds are patched builds of the default code, and the replay of the benign-side point scales has no real run to be checked against.
- The bulk outputs of the runs (annotated VCFs, per-variant tables) are under `data/benchmark/`, which is not in the repository, so a reader cannot re-derive `data/` without rerunning the benchmark.
