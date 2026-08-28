# Medical-geneticist review queue (v26 benchmark)

Like v23, this round changes no criterion logic, no threshold and no supplementary-annotation data.
The differences come from the annotation layer, plus two corrections to how the classifier *reads* that annotation.
So the question is again "is the new consequence right, and does the resulting call look right to you" - and for the first part there is an external answer.

Every variant whose call changed was re-run through **real Ensembl VEP 115.1** (Docker `ensemblorg/ensembl-vep:release_115.1`, the same GFF3 and FASTA fastVEP uses, `--pick`).
VEP picked the same transcript in all 290 cases, so the comparison is like for like.

| | before (v23) | after (v26) |
|---|---:|---:|
| Changed calls whose IMPACT tier matches real VEP | 232 / 290 | **276 / 290** |

The 14 that still do not match are all one class, and they are deliberate.
[See below](#where-we-now-disagree-with-vep-on-purpose) - fastVEP reports HIGH for a truncating variant that VEP reports as a moderate in-frame insertion.

## What changed in the annotation layer

**1. Every consequence term now comes from Ensembl's own model rather than an approximation of it.**
A change inside the CDS used to be described by one of four code paths - substitution, pure insertion, pure deletion, delins - and each had a different, partial version of the rules.
They are now one codon window, translated on both sides, with the terms read off the peptide and codon pair the way `VariationEffect.pm` reads them.
Over a 6,600-variant ClinVar sample of 150,725 transcript rows, `Amino_acids` and `Codons` went from 77.7 % and 72.5 % agreement with real VEP to **100 %**, and the coding terms from 98.0 % to **100 %**.

**2. Splice terms come from one pass over the introns, and the extended terms are not a fallback.**
`splice_donor_5th_base_variant`, `splice_donor_region_variant` and `splice_polypyrimidine_tract_variant` sit *beside* `splice_donor_variant` in Ensembl's output; fastVEP was hiding them behind it.
`intron_variant` likewise sits beside an exonic term.
Splice terms now match real VEP on all 150,725 rows.

**3. HGVS follows the 3'-rule before it names anything.**
An insertion in a repeat is written at its most 3' position, and only there does it read as a duplication: `c.4177_4182dup` where fastVEP wrote `c.4172_4173insCACCAG`.
A span with one end in an exon and the other in an intron is written over both ends.
HGVSc agreement went from 87.4 % to **99.5 %** on coding rows, HGVSp from 83.4 % to **98.8 %**.

## What changed in how the classifier reads it

**4. PVS1 no longer stands down on a splice deletion that reaches the dinucleotide.**
PVS1's splice track applies only to the canonical ±1/±2 bases, and the gate reads that offset off the HGVSc string.
A span like OAT `c.764_771+9delinsTTAGCTGTTTGTATCACACCA` starts in the exon and runs nine bases into the intron - it covers +1 and +2 on the way - but the offset was read from the *endpoint*.
This only became visible once those spans had an HGVSc at all.
Fourteen ClinVar-pathogenic donor and acceptor deletions get PVS1 back: OAT, CEP290, MSH6, ATM (×2), MLH1, MYO7A, PYGL, IFT140, CLN5, HADHA, CYP4V2, SERAC1, BRCA1.

**5. A duplication is not written across an exon boundary.**
Offsets count from the exon and do not run through zero, so stepping four bases back from `+1` does not reach `-2`.
The `ins` to `dup` conversion did exactly that, naming three bases in the previous intron.

## The 290 changed calls

`call_changes_v23_to_v26.tsv` carries every one, with the ClinVar 2-star+ truth, both fastVEP calls and criteria, both consequences, and what VEP 115.1 says.
The `review_bucket` column sorts them:

| Bucket | n | What it is | Needs review? |
|---|---:|---|---|
| `A_gained_PLP_agrees_ClinVar` | 188 | Now P/LP, and ClinVar says P/LP | No - recovered true positives |
| `B_gained_PLP_opposes_ClinVar` | 3 | Now P/LP, ClinVar says B/LB | **Yes** |
| `C_gained_PLP_ClinVar_VUS` | 18 | Now P/LP, ClinVar says VUS | Worth a look |
| `D_lost_PLP_ClinVar_says_PLP` | 1 | No longer P/LP, ClinVar says P/LP | **Yes** |
| `G_reweighted_within_PLP` | 73 | P to LP or LP to P | No |
| `H_other` | 7 | VUS to LB and similar | No |

Against the whole 673,660-variant benchmark:

| | v23 | **v26** | change |
|---|---:|---:|---:|
| Exact match | 425,235 | **425,322** | +87 |
| Same-direction | 522,297 | **522,473** | +176 |
| Opposite-direction | 61 | 64 | +3 |

Collapsed to the three directions a clinician reads:

| ClinVar / fastVEP | P/LP | VUS | B/LB | NoCall | total |
|---|---:|---:|---:|---:|---:|
| **P/LP** | 58,374 → **58,561** | 35,408 → 35,221 | 30 | 0 | 93,812 |
| **VUS** | 3,601 → 3,619 | 285,997 → 285,979 | 5,700 | 0 | 295,298 |
| **B/LB** | 31 → **34** | 106,211 → 106,201 | 177,926 → 177,933 | 382 | 284,550 |

P/LP recall 62.22 % → **62.42 %** (+187 recovered pathogenic calls).
False-benign (ClinVar P/LP called B/LB) unchanged at 30.
False-pathogenic (ClinVar B/LB called P/LP) 31 → 34, and all three new ones are in bucket B below.

## The two buckets that need you

### B. Three variants that became P/LP against a ClinVar benign call

| Gene | Variant | ClinVar | v26 | VEP 115.1 |
|---|---|---|---|---|
| MSH6 | `2:47806358 G>GGTAT` | Likely_benign (2*) | frameshift + splice_region, HIGH, **LP** `[PVS1&PM2_Supporting]` | `frameshift_variant,splice_region_variant`, HIGH |
| MSH6 | `2:47806651 G>GGTAACTAACTAACTA` | Likely_benign (2*) | stop_gained + frameshift + splice_region, HIGH, **LP** | `inframe_insertion,splice_region_variant,stop_retained_variant`, MODERATE |
| DSP | `6:7565520 C>CGTAA` | Likely_benign (2*) | frameshift + splice_region, HIGH, **LP** | `frameshift_variant,splice_region_variant`, HIGH |

All three are the same shape: **a pure insertion at an exon's last base**.
Read where the VCF puts it, the insertion is in the exon and shifts the frame, which is what both fastVEP and (on two of the three) VEP report.
Read under the HGVS 3'-rule, it slides into the intron - VEP's own HGVSc for the MSH6 variant is `c.3801+3_3801+6dup`, entirely intronic - and then it is a small intronic duplication with no coding effect.
ClinVar's submitters evidently took the second reading.

PVS1 already has a gate for exactly this ambiguity, added in the round-2 review: a pure insertion at a *canonical splice site* does not get PVS1 without positive SpliceAI support, because inserting a base beside the dinucleotide need not destroy it.
That gate lives on PVS1's splice track, and these three reach PVS1 on its **frameshift** track instead, where it does not apply.

**The question for you**: should the same "a pure insertion may not do what its unshifted position suggests" reasoning extend to the frameshift track when the insertion sits on the last base of an exon?
That is a criteria decision, and we would want your view before making it.

### D. One pathogenic variant that is no longer P/LP

| Gene | Variant | ClinVar | v23 | v26 | VEP 115.1 |
|---|---|---|---|---|---|
| NBN | `8:89958851 AATCCTGTAAATCACACA>CTTTCTACTTGTGTGATT` | Pathogenic/Likely_pathogenic (2*) | splice_acceptor, HIGH, **LP** `[PVS1&PM2_Supporting]` | coding_sequence + intron, MODIFIER, VUS | `coding_sequence_variant,intron_variant`, MODIFIER |

This is an **inversion**: the alternate allele is the reverse complement of the reference, `c.995-23_998inv`.
Ensembl matches a variant against a splice site over the bases that actually *differ* between the two alleles, not over the whole span, and an inversion leaves its palindromic positions unchanged.
Here the differing bases straddle the acceptor without covering it, so VEP calls it MODIFIER and fastVEP now agrees.

**The question for you**: an 18-base inversion spanning an acceptor site is very likely to disrupt splicing whatever the base-by-base comparison says.
If inversions over a splice site should keep PVS1, that is a criteria rule (or a SpliceAI-driven one), not an annotation change - the annotation is now doing what Ensembl does.

## Where we now disagree with VEP, on purpose

Ensembl's `ref_eq_alt_sequence` has a clause that asks whether the replacement keeps the residue it starts on and introduces a terminator *anywhere*.
It does not ask whether the annotated terminator survived, so it fires on any frameshift whose new stop lands in the codon just after the insertion point, and `frameshift_variant` and `stop_gained` both defer to it.

Real VEP 115.1 therefore reports:

| Variant | ClinVar | VEP 115.1 | fastVEP v26 |
|---|---|---|---|
| BRCA1 `c.5030_5033dup` | Pathogenic (3*) | `inframe_insertion,stop_retained_variant`, MODERATE | `stop_gained,frameshift_variant`, HIGH |
| BRCA1 `c.1499_1508dup` | Pathogenic (3*) | same, MODERATE | HIGH |
| BRCA2 `c.3205_3206insAATTGCAGTCAATTAATAT` | Pathogenic (3*) | same, MODERATE | HIGH |
| TP53 `c.895_919dup` | Pathogenic/Likely_pathogenic | same, MODERATE | HIGH |

and 30 more in the ClinVar 2-star set, including SDHB, USH2A, ATM, MSH6, BRIP1, ITGB3, GATA2, SPG11 and TTN.
Reproducing this would have reported a known truncating variant as a moderate in-frame insertion, so fastVEP does not.

That single divergence is the whole of what is left between the two tools: **101 of the 117 consequence-set differences** over a 6,600-variant sample of 150,725 rows, and **369 of the 370 IMPACT differences** on a second sample chosen adversarially from the rows this work changes.
Every one of them is fastVEP saying HIGH for a variant that truncates the protein.

If you would rather fastVEP matched VEP here, that is a one-line change and we should discuss it - but it would move 34 ClinVar-pathogenic calls in the benchmark out of P/LP.

## What did not change

No ACMG criterion logic, no thresholds, no supplementary-annotation data.
The SA stack is identical to v23.
Criterion firing counts move only where the consequence or the HGVSc moved.

## Reproducing

```bash
bash analysis/acmg_benchmark/scripts/run_benchmark.sh v26
```

`concordance_summary.txt`, `concordance_matrix.csv`, `criterion_firing_rates.csv` and `rule_distribution.csv` in this directory are that run's output.
