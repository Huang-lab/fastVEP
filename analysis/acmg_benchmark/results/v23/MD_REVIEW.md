# Medical-geneticist review queue (v23 benchmark)

This round has no criteria changes.
Every difference below comes from the **annotation** layer: three consequence-calling bugs were fixed, and the ACMG criteria then fired differently because they were reading a different consequence.
So the question for this review is not "is the rule right" but "is the new consequence right, and does the resulting call look right to you".

For the first part of that question there is now an external answer.
Every variant whose call changed was re-run through **real Ensembl VEP 115.1** (Docker `ensemblorg/ensembl-vep:release_115.1`, same GFF3 and FASTA fastVEP uses, `--pick`), and VEP picked the same transcript in all 243 cases.

| | before (v22) | after (v23) |
|---|---:|---:|
| Changed calls whose IMPACT tier matches real VEP | **5 / 243** | **243 / 243** |

That is the headline.
The calls that moved were, with five exceptions, calls where fastVEP's own consequence disagreed with VEP, and they now agree.

## What changed in the annotation layer

**1. A splice site is matched against the whole variant, not its first base.**
Every splice predicate tested only the variant's leftmost base, so a deletion or delins whose *second* base landed on the donor or acceptor dinucleotide was missed.
Ensembl tests each site by overlap with the whole span.
This is the source of nearly all the movement below.

**2. A `delins` is read from its peptides.**
A change that both deletes and inserts bases was previously reported as `inframe_deletion` or `inframe_insertion` purely from the direction of its length change, with `Amino_acids` of the form `E/X` and an empty `Codons`.
It now translates the affected codon window on both sides.
VEP calls most of these `protein_altering_variant`, some `stop_gained`, some `start_lost`.

**3. A codon window is refused when the variant is not contiguous in the CDS.**
A change straddling a splice junction was translated as if its bases were adjacent in the coding sequence, inventing a residue change.
Ten such rows in a 6,600-variant sample were a false `stop_gained`.

## The 243 changed calls

`call_changes_v22_to_v23.tsv` carries every one, with the ClinVar 2-star+ truth, both fastVEP calls and criteria, both consequences, and **what VEP 115.1 says**.
The `review_bucket` column sorts them into what needs your attention:

| Bucket | n | What it is | Needs review? |
|---|---:|---|---|
| `A_gained_PLP_agrees_ClinVar` | 193 | Now P/LP, and ClinVar says P/LP | No - these are recovered true positives |
| `B_gained_PLP_opposes_ClinVar` | 3 | Now P/LP, ClinVar says B/LB | **Yes** |
| `C_gained_PLP_ClinVar_VUS` | 9 | Now P/LP, ClinVar says VUS | Worth a look |
| `D_lost_PLP_ClinVar_says_PLP` | 10 | No longer P/LP, ClinVar says P/LP | **Yes** |
| `E_lost_PLP_agrees_ClinVar` | 1 | No longer P/LP, ClinVar says B/LB | No |
| `F_lost_PLP_ClinVar_VUS` | 21 | No longer P/LP, ClinVar says VUS | Worth a look |
| `G_strengthened_within_PLP` | 3 | LP to P | No |
| `H_other` | 3 | VUS to LB, and similar | No |

Net effect on the whole 673,660-variant benchmark: **+90 exact matches, +196 same-direction**, and opposite-direction calls **59 to 61**.

## The two buckets that need you

### D. Ten pathogenic splice duplications that are no longer PVS1 (10 variants)

All ten are the same shape: a single-base duplication two bases into the intron, `c.N+2dup`.

| Gene | Variant | ClinVar | v22 | v23 | VEP 115.1 |
|---|---|---|---|---|---|
| MECR | `c.830+2dup` | P/LP | splice_donor (HIGH), **LP** | splice_region (LOW), VUS | splice_region_variant, LOW |
| COL11A1 | `c.3816+2dup` | P/LP | splice_donor (HIGH), **P** | splice_region (LOW), VUS | splice_region_variant, LOW |
| PCDH15 | `c.3717+2dup` | P/LP | splice_donor (HIGH), **LP** | splice_region (LOW), VUS | splice_region_variant, LOW |
| SLX4 | `c.1163+2dup` | P/LP | splice_donor (HIGH), **P** | splice_region (LOW), VUS | splice_region_variant, LOW |
| BRIP1 | `c.2492+2dup` | P/LP | splice_donor (HIGH), **LP** | splice_region (LOW), VUS | splice_region_variant, LOW |
| SLC7A9 | `c.1399+2dup` | P | splice_donor (HIGH), **LP** | splice_region (LOW), VUS | splice_region_variant, LOW |
| ERCC2 | `c.1479+2dup` | P/LP | splice_donor (HIGH), **LP** | splice_region (LOW), VUS | splice_region_variant, LOW |
| NEB | `c.1470+2dup` | LP | splice_donor (HIGH), **LP** | splice_region (LOW), VUS | splice_region_variant, LOW |
| GLB1 | `c.75+2dup` | P | splice_donor (HIGH), **LP** | splice_region (LOW), VUS | splice_region_variant, LOW |
| NPHP3 | `c.3812+2dup` | P/LP | splice_donor (HIGH), **LP** | splice_region (LOW), VUS | splice_region_variant, LOW |

**Ensembl agrees with the new call on all ten**: inserting a base after the `+2` position does not change the `GT` dinucleotide itself, so it is not a canonical-splice variant.
This is the same class you flagged in round 2 (`PTEN c.802-2dupA`, `BRIP1 c.2258-2dup`), where the fix was made in the *classifier* - PVS1 stands down for an insertion at the site without SpliceAI support.
The annotation layer was still calling them `splice_donor_variant`, and PVS1's gate was carrying the whole burden; now the consequence itself is right.

**The question for you**: ClinVar 2-star+ calls all ten P/LP, and the literature on `+2dup` alleles is not the same as the literature on `+2` substitutions.
If a `+2` duplication should still reach PVS1 in some genes, that is a criteria decision (an exception list, or a SpliceAI-positive path), not an annotation one, and we would want your view before adding it.

### B. Three variants that became P/LP against a ClinVar benign call (3 variants)

| Gene | Variant | ClinVar | v23 | VEP 115.1 |
|---|---|---|---|---|
| NFKB2 | `10:102397405 ACGGGT...GGGT>A` (24 bp del) | Benign (2*) | splice_donor + inframe_deletion, HIGH, **P** `[PVS1&PM4]` | splice_donor_variant, HIGH |
| FOXP2 | `7:114631524 CCAG>C` | Benign/Likely_benign (2*) | splice_acceptor, HIGH, **LP** `[PVS1]` | splice_acceptor_variant, HIGH |
| LMX1B | `9:126693802 GCTGGG...CAGGGC>G` (33 bp del) | Likely_benign (2*) | splice_donor + inframe_deletion, HIGH, **P** `[PVS1&PM2_Supporting&PM4&BP3]` | splice_donor_variant, HIGH |

All three are deletions that genuinely remove part of a donor or acceptor dinucleotide, and Ensembl calls all three HIGH.
The disagreement is therefore with ClinVar's classification, not with the annotation.
Two are in repeat-rich introns where the deletion 3'-shifts well away from the site in HGVS (`NFKB2 c.502+22_502+45del`, `LMX1B c.886+7_886+39del`), which is a plausible reason a submitter judged them benign.

## What did not change

No ACMG criterion logic, no thresholds, no supplementary-annotation data.
The SA stack is identical to v22.
Criterion firing counts move only where the consequence moved: PVS1 gains the multi-base variants that reach a splice site and loses the `+2dup` insertions.

## Reproducing

```bash
bash analysis/acmg_benchmark/scripts/run_benchmark.sh v23
```

`concordance_summary.txt`, `concordance_matrix.csv`, `criterion_firing_rates.csv` and `rule_distribution.csv` in this directory are that run's output.
