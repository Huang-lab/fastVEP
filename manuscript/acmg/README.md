# ACMG classification manuscript

A manuscript on how the ACMG/AMP classifier of fastVEP was built, reviewed by a medical geneticist and measured, how its thresholds were chosen, and what its discordance with ClinVar 2-star+ classifications looks like in the final version (v37) and on a later ClinVar release.
ClinVar's 2-star+ classifications are a reference and not ground truth, throughout.

## Files

| Path | What |
|---|---|
| `manuscript.md` | The manuscript: abstract, methods, results, discussion, references, figure legends and supplementary list |
| `figures/` | Figures 1 to 9 as PNG (300 dpi), written by `generate_figures.py` together with a PDF of each (the repository ignores `*pdf`, so only the PNGs are tracked) |
| `generate_figures.py` | Draws every figure from the tables in `data/` and nothing else |
| `build_tables.py` | Builds the tables that carry numbers (Tables 3 to 7) from `data/` and keeps `manuscript.md` in step with them |
| `verify_numbers.py` | Recomputes the numbers quoted in the text from `data/` and fails when a sentence does not contain them |
| `data/` | The supplementary tables: small CSV and TSV files written by `analysis/acmg_benchmark/scripts/13_manuscript_stats.py`, `11_manuscript_runs.py`, `14_temporal_holdout.py` and `15_clingen_complete_export.py` |
| `replay_validation/combiner_cutoff_env.patch` | A patch of eight lines that makes the two point cutoffs of the pathogenic side readable from the environment, for the nine cutoff runs of Table 4 |
| `replay_validation/bp_hooks_env.patch` | A patch that makes BP1's own-record count and BP7's SpliceAI bar readable from the environment, for the runs of Sections 3.3 and 3.4 |

The scored matrix, firing counts and provenance of every run made for the paper are in `analysis/acmg_benchmark/results/v37_manuscript/`.
The commands that regenerate everything are in `analysis/acmg_benchmark/README.md`.

## Checking the manuscript

```bash
python3 manuscript/acmg/build_tables.py --check    # the tables in the text are what the data say
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
  We did: a strict leave-one-out run (59 opposite-direction calls), the baseline curve (Figure 9), a finer allele-frequency standardisation of the benign-side fall, a sweep of BP1's veto count that the audit found had never been made, and runs of the proposed changes on the new variants (Section 3.7).
  Not done: asking the reviewer to confirm our coding of rounds 2 and 6, a comparison with another tool, and the filtering allele frequency in the per-variant tables.
- A fourth subagent then audited only the sentences that the revision had added or changed, against the data and the code.
  Every one of its thirteen main points held when we recomputed it, and the manuscript now corrects them: one of the two missed diagnoses that a stricter count of BP1's veto avoids (PRR12) rests on the variant's own record, and the other (IQSEC1) does at a count of 2 only (a run of the pair, with the own record left out, avoids one); the strict reading adds three Benign calls, so the zero for Benign calls does not survive it; a quiet variant, as the baseline rule defines it, includes the flanking variants and the polypyrimidine tract, which the text left out (without them the rule has 1,036 contradictions and not 1,730); 98.5% of the new variants are in the April protein index with their April record, because they were in ClinVar at a lower review status, so the hold-out of BP1's own-record count is a weak test and the new variants are not variants that ClinVar did not hold; the date skew of the omitted ClinGen curations does not fit a random omission; the baseline family has fewer contradictions at higher bars and without a frequency rule than the text said; only one further source (Orphanet) was measured and not two; and several sentences compared numbers of different kinds.
  Three of the points needed new runs or tables (BP1's count with the own record left out, BP1's count on the new variants, and why BA1 and BS1 did not fire on the high-frequency new Benign variants, which is gnomAD's low-complexity, filter and segmental-duplication flags and the classifier's homologous-gene list for 79% of them).
- A fifth subagent audited those corrections.
  Its five points of meaning held when we recomputed them and are corrected: the removal of IQSEC1 by a stricter count of BP1's veto rests on its own record at a count of 2 and not at a count of 1, the protein index holds the April record only of the new variants that were then classified missense variants (77 of the 94 that BP1 can act on), the rules without a frequency rule agree for 0.3 points more than v37 and not the same, 23 of the 31 new Benign variants that met BS1 and stay at VUS have BS1 alone, and the zero for Benign calls rests on two of the 15 curated alleles and one gate.
- The ClinGen export of 30 April 2026 was compared with an export downloaded on 11 October 2026, which is how the 104 omitted curations were found.
- Not re-verified: the exact strength grades of PS1 for a splice variant (Walker et al., 2023, Table 3) against our implementation, and the claim that Ensembl VEP and fastVEP agree on each annotation difference that the ledger mentions.
- The ledger of 64 changes, its classes and its routes of discovery are our coding of `RUN_VERSIONS.md` and the round documents.
  They have not been coded independently.

## Findings of the verification that are not in the classifier

The manuscript reports these in Sections 3.2 to 3.5, with what each would change, and the repository's comments and documentation were corrected for the citations and for the attributions.
The classifier was not changed, because each is a decision that deserves its own review.

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
- The validity gate misses NSMCE3, a gene whose ClinGen curations are all weak, because its label ends in a semicolon (one of 505 weak-only genes, with no benchmark variant).
- The label script maps two RYR1 expert-panel records with a conflicting classification to VUS and four CFH and CFI records classified Likely pathogenic with low penetrance to Pathogenic.
  We kept the labels so that the 37 runs stay comparable; relabelling changes exact agreement by one variant, and leaving out the two RYR1 records lowers the ClinVar VUS called pathogenic-side by two.
- Two comments in the classifier were wrong and are corrected: the one that said the `LOF_and_GOF` mechanisms are inert (they license PVS1 through the curated-mechanism route, which changes no call in the benchmark), and the one that said the ClinGen phenotypes are Moderate or better only (they are every classification, and the validity gate is what withholds the weak ones).

## Decisions that are the author's

These are marked in the manuscript with square brackets or are not in it at all.

1. **The reviewer.**
   The manuscript describes the reviewer's rulings and paraphrases their reasons.
   It quotes only the heading of the list that came with the round-1 workbook ("Discordance categories") and does not reproduce their notes, and `data/round6_rulings.tsv` holds our coding of the text and not the text.
   Whether to name them, in the acknowledgements or as an author, and whether the text may be quoted or published as a supplementary column, is for the author to settle with them.
   The workbooks they returned are not in the repository.
2. **Disclosure of AI assistance.**
   The scripts, the figures and the first draft were produced with Claude working as a coding agent in this repository.
   The Declarations section holds a placeholder for the author's statement.
3. **Author contributions, funding and competing interests** are placeholders.
4. **The title and the abstract** state the main findings and should be read by the author before submission.
5. **Pronouns.**
   The manuscript refers to the reviewer as "the reviewer" and "they", because their pronouns are not recorded.
6. **Whether to change the classifier** on the findings above, and in what order (Section 4.5).
   Each has a price measured in this work, and none should be tuned to the 49 residual calls.
7. **Archiving.**
   The bulk outputs (annotated VCFs, per-variant tables of about 20 MB each) are not in the repository.
   A tagged release with a DOI, and a compressed deposit of the per-variant tables, would let every table in `data/` be re-derived without a re-run.

## What is not checked

`verify_numbers.py` guards the numbers, not the wording.
A sentence can pass and still say something the data do not support.
It does not check the following, which are not functions of `data/`.

- Dates, version numbers and constants of the code (the 2,000-allele floor, the pLI and LOEUF bars, the missense z-score bars, the 15 alleles of the exception list, the 26 criteria that a run reports and the 10 that are not evaluated).
- The counts read from the project's own record: the split of the 2.7 points that v9 lost (0.6 and 2.0), the 10,000-row cap and the chromosome at which it stopped the first export, the 3,793 Likely pathogenic calls that PS1's splice path turned into Pathogenic, the 14 deletions that lost their PVS1, and the 19-point knee that v10 recorded.
- What each recommendation says (Section 3.2 and Table 1), which was checked by subagents and by reading, and the release dates of ClinVar and ClinGen inputs.
- The reviewer's remarks, which are paraphrased in Section 3.6 and coded, not stored.

Other limits:

- The hold-out of Section 3.7 is a hold-out of variants and not of genes, and the ClinVar-derived indices of the classifier were not rebuilt for it.
  Its new variants are mostly variants that ClinVar held in April at a lower review status (98.5% have a record in the `variant_summary` that the protein index was built from), so they are variants that no sweep or review scored, and the protein index holds the April record of those that were then classified missense variants.
- The cutoff and measurement builds are patched builds of the shipped code, and the replay of the benign-side point scales has no real run to be checked against.
- The bulk outputs of the runs (annotated VCFs, per-variant tables) are under `data/benchmark/`, which is not in the repository, so a reader cannot re-derive `data/` without rerunning the benchmark.
