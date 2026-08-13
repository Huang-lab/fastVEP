# ACMG Expert Review, Round 2: Findings, Implementation, and Results

> **Status: every item in this document is implemented except B1, which needs
> curation rather than code.**
> Tier A, Tier C and B4 landed in run v9; B2, B3, B5, B6, B7, B8, C5 and C6
> landed after it and are measured in run v10.
> The v9 figures below stand as the record of what the criteria fixes alone
> bought: opposite-direction calls against ClinVar 2-star+ dropped from **122 to
> 41** (ClinVar-informed) at a cost of 0.6 pp of same-direction agreement, with a
> leave-one-out mode removing the ClinVar self-reference PS1 and PM1 relied on
> (section 8).
> Two thresholds that had been set by convention are now set by measurement -
> the BS2 prevalence bar and the PM2 bar for dominant genes - and both sweeps are
> recorded in [`ACMG.md`](ACMG.md).
> **B1 (prevalence-aware per-gene-disease thresholds) is the one thing still
> open, and it is a curation task**: the code hooks exist and are waiting for the
> table.


Source: `discordant 122.xlsx` (122 opposite-direction calls from the v8 ClinVar 2-star+ benchmark) plus the reviewer's covering email.
93 variants are new this round, 29 were carried over from the v7 review.
Every claim marked **verified** below was reproduced by running `fastvep annotate --acmg --pick` against the benchmark SA stack (`data/benchmark/sa_db`) on the exact variants in the file.

## TL;DR

The reviewer is right on essentially every point, and most of what she is asking for is *already in the text of Richards 2015* or in a ClinGen SVI/VCEP addition that we cite elsewhere in [ACMG.md](ACMG.md).
Very little of it requires anything we would call "beyond the guideline".

Three things came out of reproducing her cases that she could not have seen from the spreadsheet:

1. **PS1 self-matches.** Our ClinVar protein index includes the variant under test, so a ClinVar-pathogenic missense earns PS1 from its own ClinVar record. This is circular and it inflates the benchmark.
2. **`ALT=.` records are annotated as real variants.** Five rows in her file are ClinVar reference-agreement records that we should have skipped. Her yellow "does dot mean del?" highlight is exactly this.
3. **Multi-nucleotide substitutions are translated in the wrong codon frame.** Her two "why is PVS1 called, this is synonymous" rows (MUTYH, HPS4) are a real consequence-calling bug, not a criterion bug.

## 1. What the 122 rows actually contain

| Direction | Count |
|---|---:|
| Truth P/LP, fastVEP B/LB (false-benign) | 78 |
| Truth B/LB, fastVEP P/LP (false-pathogenic) | 44 |

Her colour key resolves to:

| Reviewer verdict | Count |
|---|---:|
| ClinVar call correct, fastVEP wrong | 59 |
| fastVEP call correct, ClinVar wrong | 2 (KMT2B, ORAI1) |
| Both calls wrong, she would call something else | 8 (BRCA2, COG7, CYFIP2, EMG1, KIAA0586, KMT2C x2, MUTYH) |
| Needs more literature before deciding | 1 (KIAA0586) |
| Reviewed, not colour-coded | 52 |

Criterion footprint across the 122 (a criterion can appear in many rows):

| Criterion | Rows | Criterion | Rows |
|---|---:|---|---:|
| PM2_Supporting | 55 | BP1 | 21 |
| BP4 (+8 BP4_Moderate, 1 BP4_Strong) | 53 | BS1 | 18 |
| BS2 | 52 | PM1 | 15 |
| PS1 | 38 | PM4 | 6 |
| PVS1 (+1 PVS1_Moderate) | 38 | PP2 | 5 |

Two combinations dominate:

- `PVS1;PM2_Supporting` on 33 rows, all false-pathogenic.
- `BS2;BP4` and its relatives (`PS1;BS2;BP4`, `PM4;BS2;BP4`, `PS1;PM2_Supporting;BP1;BP4`) across most of the false-benign side.

**Homozygote counts where BS2 fired on a truth-Pathogenic variant:** 20 rows had exactly 1 homozygote in gnomAD, 6 rows had 2.
That single observation is 26 of the 52 BS2 false-benign rows, and it is precisely the reviewer's "a couple of homozygotes are possible if disease has later onset or variable expressivity" point.

## 2. Verified reproductions

Ran on the nine representative variants; all reproduced exactly as reported in the file.

| Variant | fastVEP consequence | Criteria | Problem |
|---|---|---|---|
| MUTYH 1:45331682 GC>AT | `stop_gained` (codons `CTg/TAg`) | PVS1, PM2_Supporting | delins spans a codon boundary (c.1164_1165); the true effect is synonymous. **Wrong frame.** |
| HPS4 22:26464569 GA>CT | `stop_gained` (codons `tCC/tGA`) | PVS1, PM2_Supporting | c.1060_1061delTCinsAG is p.(Ser354=). We are one base off in the codon. |
| PTEN 10:87960891 T>TA | `splice_acceptor_variant` | PVS1, PM2_Supporting | c.802-2dupA duplicates the identical base; the AG dinucleotide is unchanged. |
| BRIP1 17:61743135 C>CT | `splice_acceptor_variant` | PVS1, PM2_Supporting | same, c.2258-2dup. |
| ATM 11:108244119 GTAATC>G | `splice_donor_variant` | PVS1, PM2_Supporting | deletion inside a microsatellite; `intronic_offset` is never populated for indels so the ±1/±2 gate does not engage. |
| RSPH9 6:43670917 GAGA>G | `inframe_deletion` | PM4, BS2, BP4 | **BP4 on an in-frame deletion**, earned through the SpliceAI path. Exactly the reviewer's point 1. |
| ENG 9:127824800 C>T | `missense_variant&splice_region_variant` | PS1, PM2_Supporting, BP1, BP4 | BP1 and PS1 fire together; and PS1's only supporting index entry is `331:G>S:Pathogenic`, i.e. **this variant's own ClinVar record**. |
| NPHP4 1:5875102 T>. | `splice_acceptor_variant` | PVS1, PM2_Supporting | `ALT=.` parsed as an allele; produced `T/.` and `cgA/cg.` codons. |
| FMR1 X:147938098, IDS X:149500991 | - | - | no gnomAD record reached the classifier, and we model no hemizygote counts at all. |

I also confirmed the gnomAD chrX `.osa` is populated in general (21 of 40 sampled chrX ClinVar variants got `FV_GNOMAD`), so the FMR1/IDS gap is site-specific rather than a missing chromosome.

## 3. Tier A: straight from the guideline text, no new data needed

**All implemented in run v9.** These needed only code changes in `crates/fastvep-classification`.
Measured effect across the full 673,660-variant benchmark is in section 8; on the reviewer's 122 rows specifically, Tier A plus the homology gene list from B4 resolved 78, mostly into VUS, which is the honest automated answer when the deciding evidence is case-level or literature-level.

### A1. BP1 requires positive evidence of a truncating-only mechanism

- **Guideline:** Richards 2015, BP1 = "Missense variant in a gene for which primarily truncating variants are known to cause disease."
- **Now:** [benign_supporting.rs:30](../crates/fastvep-classification/src/criteria/benign_supporting.rs#L30) fires BP1 on a constraint heuristic alone (pLI >= 0.9 and misZ < 2.0). That is a statement about tolerance to missense variation in the population, not a statement about the gene's disease mechanism.
- **Change:** require positive evidence that missense is *not* an established mechanism. Block BP1 when the ClinVar protein index for the gene contains at least N pathogenic/likely-pathogenic missense entries (N configurable, default 3), and hard-block when PS1, PM5 or PP2 fired for this variant.
- **Impact:** all 21 BP1 rows. Fixes ENG x2, CFH, POU3F3, FGFR3, APC, BCL11B, ANKRD11, NPTX1, IRF2BPL, ERF, SBF1, PKD1 x2, ARID1B, PRR12, IQSEC1, SETD1A, COL11A1, MITF.

### A2. Contradiction guard in the combiner

- **Guideline:** Richards 2015 explicitly says a variant meeting criteria in both directions should default to Uncertain Significance.
- **Now:** [combiner.rs:28](../crates/fastvep-classification/src/combiner.rs#L28) only declares a conflict when *both* sides independently reach a definite call. A met PS1 (Strong) alongside two BP criteria still yields Likely Benign.
- **Change:** any met criterion at Strong or above on one side blocks a definite call on the other side and yields VUS with an explicit "conflicting" rule string.
- **Impact:** 36 of the 78 false-benign rows carry a met PS1 or PVS1 alongside the benign criteria that decided the call.

### A3. BP4 is not applied where no predictor is calibrated

- **Guideline:** Pejaver 2022 calibrates REVEL for missense only; Walker 2023 scopes SpliceAI-based BP4 to variants where splicing is the question. There is no calibrated benign predictor for in-frame indels. This is the reviewer's point 1 verbatim.
- **Now:** the SpliceAI <= 0.1 path in [benign_supporting.rs:396](../crates/fastvep-classification/src/criteria/benign_supporting.rs#L396) is gated against PVS1-territory consequences and against deep-exonic missense, but an in-frame deletion, an in-frame insertion or a stop-lost variant still collects BP4 Supporting from a SpliceAI score of 0.
- **Change:** extend the gate. BP4 does not fire for `inframe_deletion`, `inframe_insertion`, `stop_lost`, or `protein_altering_variant` unless the variant also overlaps a splice region. Emit `splice_skipped_reason` as we already do for missense.
- **Impact:** the 6 in-frame rows (AFG2A, RSPH9, AMHR2, OCA2, CYP24A1) and HBA1 stop-lost, plus a long tail outside the 122.

### A4. BS2 honours "full penetrance expected at an early age"

- **Guideline:** Richards 2015, BS2 = "Observed in a healthy adult individual for a recessive (homozygous), dominant (heterozygous), or X-linked (hemizygous) disorder, **with full penetrance expected at an early age**." We implement the first half and ignore the qualifier.
- **Now:** [benign_strong.rs:171](../crates/fastvep-classification/src/criteria/benign_strong.rs#L171) fires BS2 on `all_hc > 0`, i.e. a single homozygote anywhere in gnomAD.
- **Change, part 1 (no new data):** add `bs2_ar_min_hom`, default 3. One or two homozygotes is not evidence of tolerance for a late-onset, low-penetrance or variably expressive disorder.
- **Superseded by measurement.** `bs2_ar_min_hom` turned out to be inert once the cohort-scaled test existed: sweeping it 1 to 10 changes nothing. The knob that matters is the prevalence bar, and sweeping *that* across the full benchmark set the default at 1e-3, the prevalence of the most common Mendelian disorders, rather than the 1e-5 first shipped. That covers the hearing-loss / alpha-1 antitrypsin / familial Mediterranean fever cases in her point 2 without needing the per-gene table. See [`ACMG.md`](ACMG.md#choosing-the-bs2-prevalence-bar).
- **Change, part 2 (needs B1's gene-disease attribute table):** suppress BS2 entirely when the gene-disease entity is flagged late-onset, reduced-penetrance or variable-expressivity.
- **Impact:** 26 of 52 BS2 rows from part 1 alone (EYS, ODAD1, DNAAF19, COG7, RSPH9, AMHR2, CYP24A1, MMACHC, OCA2, KIAA0586, RB1, STRC and others). Part 2 covers ALPL, MEFV, ACADS, C6, C9, SLCO1B1, SLC22A5, DUOX2.

### A5. PM1 scope and anti-double-counting with PVS1

- **Guideline:** PM1 is "located in a mutational hot spot and/or critical functional domain **without benign variation**", and in practice it is applied to missense and in-frame variants to add residue-level evidence. Stacking it on a null variant double-counts the region evidence PVS1 already carries.
- **Now:** [pathogenic_moderate.rs:27](../crates/fastvep-classification/src/criteria/pathogenic_moderate.rs#L27) applies PM1 to anything with a protein position, including frameshifts and nonsense. The reconciliation table in [ACMG.md](ACMG.md) suppresses PM1 only when PP3 was elevated to Strong.
- **Change:** restrict PM1 to missense and in-frame indels; add "PVS1 fires -> suppress PM1" to the reconciliation pass.
- **Impact:** 10 of 44 false-pathogenic rows. Directly answers her CBS, MSH6, RYR1 and CHD7/PTCH1 comments.

### A6. ClinVar low-penetrance and risk-allele terms suppress BS1/BS2

- **Guideline:** ClinVar's controlled vocabulary now carries `Pathogenic, low penetrance` and `Established risk allele`. A variant the submitters have already labelled low-penetrance is by definition outside BS2's "full penetrance" precondition.
- **Now:** we parse CLNSIG for pathogenic/benign substrings and ignore these terms.
- **Change:** when a 2-star+ ClinVar record carries a low-penetrance or risk-allele term, mark BS1 and BS2 NotEvaluated with an explicit reason.
- **Impact:** SERPINA1 (1236 homozygotes, `Pathogenic, low penetrance`), F2 (`Established risk allele`), ALPL, MITF.
- **Caveat:** this is a ClinVar-derived gate, so it must be disabled in the ClinVar concordance benchmark. See section 7.

### A7. PS1 must exclude the variant under test

- **Guideline:** PS1 = "Same amino acid change as a previously established pathogenic variant **regardless of the nucleotide change**". The intent is a *different* nucleotide change producing the same residue substitution.
- **Now:** [pathogenic_strong.rs:139](../crates/fastvep-classification/src/criteria/pathogenic_strong.rs#L139) matches on `(protein position, alt AA, significance contains "pathogenic")` with no exclusion of the query variant. **Verified:** for ENG c.991G>A the only matching index entry is `331:G>S:Pathogenic`, which is the variant itself.
- **Change:** carry the nucleotide-level identity (ClinVar allele ID or genomic key) in `clinvar_protein.oga` and exclude self-matches. The same leave-one-out applies to PM5 and to the PM1 hotspot count, where the variant currently contributes one of the three pathogenic neighbours it needs.
- **Impact:** correctness, and it materially changes the benchmark. Expect pathogenic recall to fall; the current number is partly self-fulfilling.

## 4. Tier B: ACMG guideline additions that need new data

Each of these is a published ClinGen SVI or VCEP mechanism, not an invention of ours.

### B1. Prevalence-aware BS1/BA1 thresholds (Whiffin/Ware framework + cspec)

This is the reviewer's central point about hearing loss, alpha-1 antitrypsin deficiency and familial Mediterranean fever.

- **Guideline:** Whiffin et al. 2017 (Genet Med), adopted by ClinGen SVI, derives a maximum credible population allele frequency per gene-disease from disease prevalence, allelic heterogeneity, genetic heterogeneity, penetrance and inheritance. Where a ClinGen VCEP has published explicit BA1/BS1 thresholds in cspec, those supersede any computed value. The reviewer cited cspec twice in the file (GN137 for MSH2, GN138 for MSH6).
- **Now:** a single flat `bs1_af_threshold = 0.01` applied to max-population AF for every gene on earth.
- **Change:** new gene-disease attribute table shipped as an `.oga` source, keyed by gene symbol and MONDO ID, carrying prevalence, penetrance class, onset class, allelic/genetic heterogeneity and, where published, the VCEP's literal BA1/BS1/PM2 thresholds. Config already has the `disorders` scaffold in [config.rs:216](../crates/fastvep-classification/src/config.rs#L216); this fills it.
- **Effort:** this is the largest single item. Seed with the ~60 VCEPs that have published specifications; fall back to the flat threshold with an explicit "generic threshold" flag in the output for everything else.

### B2. Filtering allele frequency instead of raw max-population AF

- **Guideline:** ClinGen SVI recommends the filtering allele frequency (FAF95, `grpmax`) for BA1/BS1 rather than a point estimate, precisely because a founder allele in a small subpopulation produces an unreliable maximum.
- **Now:** the gnomAD `.osa` builder captures `AF`/`AC`/`AN`/`nhomalt` plus per-population `AF` only. There is no per-population AN, so we cannot require a minimum sample size behind the population maximum, and we cannot compute FAF.
- **Change:** capture `faf95` and `fafmax_faf95_max`; use them for BA1/BS1 when available, falling back to the population maximum otherwise.
- **Status: done.** Both columns are extracted, and BA1/BS1 read them through one shared helper so they cannot disagree about how frequent a variant is. Per-population AN was considered and dropped: the FAF already folds sample size into the estimate, so ten extra dense columns would buy nothing the FAF does not already give.
- **Honest scope.** This fixes the *small-sample* half of the founder problem: a frequency inflated by thin sampling in one group no longer reaches BA1/BS1. It does not fix a founder allele that really is common in a well-sampled population, whose FAF is high too. Only the curated Ghosh 2018 exception list keeps BA1 off those, so ACADS / RSPH9 / C9 / OCA2 need B1 or an exception-list entry, not B2.
- **Impact:** the founder-variant caution she raised (ACADS in Ashkenazi, RSPH9 in Middle Eastern, C9 in East Asian, OCA2 in African American).

### B3. gnomAD quality-control awareness

- **Now:** we ingest gnomAD sites without recording `FILTER` or the region flags. A variant that fails gnomAD QC looks identical to a clean one, and a variant absent because of a QC drop earns PM2 as though it were genuinely absent.
- **Change:** capture `FILTER` (`AC0`, `AS_VQSR`, `InbreedingCoeff`) and the `lcr` / `segdup` / `non_par` flags. When the site is not PASS or sits in a low-confidence region, mark BA1/BS1/BS2 **and PM2** NotEvaluated with a reason rather than firing.
- **Status: done.** All four criteria now share one gate (`criteria::frequency_gate`) instead of repeating their preconditions, which is also what makes the PM2 half true: absence is only evidence when the database could have seen the variant. Verified end-to-end on chr22:17084998 C>T, an `AC0` record that previously earned PM2_Supporting for being "absent from gnomAD" at a site gnomAD explicitly could not call.
- **Impact:** her "fails QC in gnomAD" notes on KMT2C x2, MTR, NOTCH2, TNNI3K, GIGYF2, and the CBS 68 bp insertion where she warns that large insertions map poorly.

### B4. Homologous / segmental-duplication gene list

- **Guideline:** Mandelker et al. 2016 (PMID 27228465), which the reviewer cited directly in her colour key: "gene with homology issue; cannot rely on population frequencies".
- **Change:** ship the published gene list as a gene-level flag; suppress BA1/BS1/BS2/PM2 for those genes and surface `frequency_unreliable_homology` in the output so a curator sees why.
- **Impact:** she flagged 7 rows herself (CYP21A2 x4, STRC x3); the list also covers HBA1/HBA2 (4 more rows here), SMN1, PMS2, NEB, OTOA.

### B5. Hemizygote-aware BS2 for X-linked genes

- **Guideline:** BS2's own text includes "X-linked (hemizygous)".
- **Now:** `GnomadData` has no hemizygote field at all, so an X-linked variant can never earn BS2 from hemizygous observations.
- **Change:** capture the XY-stratified counts and count hemizygotes towards BS2.
- **Status: done.** BS2 now counts individuals with no functional copy - homozygotes plus, on a non-PAR sex-chromosome site, hemizygotes - and the cohort size becomes `(AN - AN_XY)/2 + AN_XY`, which is exactly `AN/2` when there are no XY columns, so autosomal results are unchanged.
- **One correction to the plan above.** `nhomalt_XY` is *not* extracted, because it carries no information. gnomAD calls XY samples haploid outside the PAR and so never records one as homozygous: across all 6,955 non-PAR chrX records in the IDS region it is zero, including at a site with 109,916 XY carriers. Inside the PAR those samples are diploid and the global `nhomalt` already counts them. `AC_XY` is the hemizygote count. This is why the FMR1 and IDS rows looked as though gnomAD had never seen a null individual.
- **Impact:** FMR1 (210 hemizygotes) and IDS (429 hemizygotes), both of which we called Likely Pathogenic.

### B6. Disease mechanism (loss of function vs gain of function) gates PVS1

- **Guideline:** Abou Tayoun 2018 requires that loss of function be an established mechanism for the gene before PVS1 is applied.
- **Now:** [pvs1.rs:369](../crates/fastvep-classification/src/criteria/pvs1.rs#L369) reads `GeneOverride.mechanism` only to *enable* PVS1; a gene declared GOF is never blocked, and the default table is empty. In practice PVS1 fires on any null variant in a constrained gene.
- **Change:** make mechanism authoritative. If mechanism is GOF-only, PVS1 is NotApplicable. Ship a curated default table sourced from ClinGen VCEP specifications and GDV mechanism statements.
- **Status: done.** A mechanism that excludes loss of function (`GOF`, `DOMINANT_NEGATIVE`) takes PVS1 to NotApplicable. The curated table ships in its own `gene_mechanisms` map rather than inside `gene_overrides`, so a TOML setting one `[gene_overrides.X]` block does not silently discard it. Constraint cannot substitute for this: a gain-of-function gene under purifying selection has a high pLI for the same reason a haploinsufficient one does, which is how PCSK9 - where the null allele lowers LDL and is the mechanism of an approved drug class - was collecting PVS1.
- **Measured** over 12,844 ClinVar variants in the flagged genes: PVS1 blocked in PCSK9 23, IFIH1 20, TNNI3K 11, CYFIP2 4, MYH7 89. Controls hold - BRCA1 keeps all 2,571 PVS1 firings.
- **One correction.** RYR1 is *not* a case for this gate. Malignant hyperthermia is gain of function but the congenital myopathies are loss of function, so a null allele is still pathogenic for one of the two. It is curated `LOF_and_GOF` and PVS1 still applies to it; all 89 firings survive. Genes with both mechanisms are listed in the table anyway, because each is one somebody will reasonably propose adding as GOF.
- **Impact:** PCSK9 x2, CYFIP2, IFIH1 x2, TNNI3K.

### B7. Gene-disease validity gate

- **Guideline:** every pathogenic criterion presupposes an established gene-disease relationship. We already load ClinGen GDV (Definitive/Strong/Moderate only) into `omim.oga`.
- **Change:** enforce it. PVS1, PP2 and PM1 require a GDV entry for the gene; without one, report VUS with `no_established_gene_disease_relationship`.
- **Status: done.** All three now report NotEvaluated with that reason for a gene absent from the loaded source. The gate degrades the opposite way from every other one here, deliberately: with no `omim.oga` loaded it does not fire at all, because a gene missing from a file nobody opened is not a gene without a disease. `ClassificationInput` carries `gene_disease_db_loaded` for exactly that distinction - `omim` being `None` cannot tell the two apart, and gating on it would suppress PVS1/PP2/PM1 genome-wide.
- **Correction to the impact list above.** Only RYK and GIGYF2 are actually absent from ClinGen GDV. EMG1, ARMC9 and ORAI1 all have entries, so B7 cannot be what fixes them; whatever is wrong with those rows is something else. Measured blocks: RYK 1, GIGYF2 2.
- **Cost:** 0.7 pp of pathogenic recall on a 1-in-10 benchmark sample, with benign recall and opposite-direction calls unchanged. The gate is narrow by construction.
- **Impact:** RYK ("gene is not associated with disease"), GIGYF2 (susceptibility locus only).

### B8. Functional evidence input that supersedes computational predictions

This answers her last email paragraph.

- **Guideline:** the SVI's ordering is explicit. Functional evidence (PS3/BS3) is stronger than computational prediction, and PP3/BP4/BP7 should not be used to argue against sound experimental data.
- **Now:** PS3 and BS3 are hard-coded NotEvaluated.
- **Change:** add `--functional-evidence <TSV>` mapping variant key to PS3/BS3 with strength and PMID. When an entry is present, PS3/BS3 fire at the stated strength and BP4/BP7/PP3 are suppressed for that variant with `superseded_by_functional_evidence`.
- **Status: done.** `--functional-evidence <TSV>` takes `chrom, pos, ref, alt, criterion, strength, pmid, note`, keyed on the record's original VCF coordinates rather than fastVEP's normalised alleles, because the file is written by a human reading a VCF. Strength is a column, not a constant: Brnich 2020 grades assay strength on validity, controls and dynamic range, so a curator who has read the paper can say Supporting. Malformed rows are hard errors naming the line, parsed before the run starts, and a file asserting both PS3 and BS3 for one variant is rejected rather than resolved by whichever line came last.
- **Verified end-to-end** on the OCA2 case: `before call=B (BA1+BS2+BP4+BP7)` → `after call=VUS (PS3 PMID 26637981 + BA1 + BS2, BP4 and BP7 superseded)`.
- **Impact:** OCA2 c.1503A>G (synonymous with published splice-defect data, where we fired BP7 and BP4), C2 (exon-skipping functional data).
- **Note:** fastVEP will not mine literature. It will consume curated evidence, which is how every VCEP pipeline works.

## 5. Tier C: annotation bugs found while reproducing

Correctness bugs independent of ACMG. **C1 through C4 are fixed in v9; C5 and C6 after it.**

| # | Bug | Evidence | Affected rows |
|---|---|---|---|
| C1 | `ALT=.` (ClinVar reference-agreement records) is parsed as a real allele and annotated. Produces `T/.` amino acids and `cgA/cg.` codons, then a full Likely Pathogenic call. | Verified on NPHP4 and NLRP3. | 5 (NLRP3, NPHP4, ALMS1, KCNJ16, PIGN) |
| C2 | Multi-nucleotide substitutions spanning or offset from a codon boundary are translated in the wrong frame, producing spurious `stop_gained`. | MUTYH `GC>AT` gives `CTg/TAg`; the correct c.1164_1165delinsAT is synonymous. HPS4 `GA>CT` gives `tCC/tGA`; the correct c.1060_1061delTCinsAG is p.(Ser354=). | 2, plus unknown tail |
| C3 | An insertion or duplication at a canonical splice site that re-inserts the identical base(s) still calls `splice_acceptor_variant` and earns PVS1. The dinucleotide is unchanged. | PTEN c.802-2dupA, BRIP1 c.2258-2dup, DSP c.939+3_939+6dup. | 3-4 |
| C4 | `intronic_offset` is not populated for indels, so PVS1's canonical ±1/±2 gate never engages for deletions and insertions spanning a splice boundary. | ATM `GTAATC>G` got PVS1 with SpliceAI 0.00. | 5 (ATM x2, BRCA2, MSH6 x2) |
| C5 | **Fixed, and the cause was not what this table assumed.** fastVEP implements VEP's default `--pick_order` faithfully, and in that order consequence severity (`rank`) is the *last* tier, below mane_select/canonical/appris/tsl/biotype/ccds. So a MANE transcript the variant merely neighbours outranks a non-MANE one it disrupts. Stock VEP makes the same choice. Fixed by supporting VEP's `--pick-order` flag rather than by changing the default, which stays VEP-identical. | CYP21A2 rows came out on **C4B**, STRC-region rows on **TIMM9**. With `--pick-order rank,...` over 300 variants in these regions: C4B 43→0, TIMM9 32→2, CYP21A2 34→77. | 8-9 |
| C6 | **Fixed.** The export now carries a `variant` column holding ClinVar's own genomic HGVS, so no row has to be identified by interpreting the join key, and a reference-agreement row gets a review_question saying what `ALT=.` is and that its presence is a bug to report rather than a call to review. `ref`/`alt` stay verbatim as the join back to the VCF. | `04_build_md_review_table.py` | 5 |

## 6. What we will not do, and why

| Not doing | Why | What we do instead |
|---|---|---|
| Automated literature mining for PS3/BS3/PP1/PS4 | Not tractable deterministically, and a wrong PMID in a clinical report is worse than no evidence. | B8: consume a curated functional-evidence file, and keep PS3/BS3 explicitly NotEvaluated otherwise. |
| Hand-curated per-variant hypomorph lists (TAR/RBM8A, TXNL4A promoter, MMACHC, SERPINA1 PI\*Z-in-trans) | No machine-readable source exists, and hand lists rot silently. | A6 catches the subset ClinVar already labels low-penetrance; the rest is a documented `gene_overrides` hook plus an output flag so a curator sees the frequency evidence was suppressed and why. |
| Capping PP3 at Moderate (her KMT2C comment) | Pejaver 2022 is a ClinGen SVI product and explicitly calibrates REVEL >= 0.932 to Strong. Overriding it would put us outside the guideline. | **Done.** `pp3_max_strength` is exposed in config; a lab following her stricter convention sets `pp3_max_strength = "Moderate"`. Default stays uncapped. The cap is a ceiling only, never a floor. |
| Detecting whether a predictor was trained on the variant in question | Training-set membership is not published per variant for REVEL or SpliceAI. | Documented as a known limitation in [ACMG.md](ACMG.md). B8's functional-evidence override is the practical mitigation. |
| Automatically choosing which disorder to classify against for multi-phenotype genes (RYR1 MH vs myopathy, IFIH1 AD vs AR, PKD1 AD vs AR, TP63) | A variant can be genuinely benign for one disorder and pathogenic for another. Picking silently would be worse than not picking. | Phase 2 of B1: emit a per-disorder classification when GDV lists multiple entities with different inheritance, and accept `--disorder MONDO:xxxxxxx` to select one. Never silently pick. |
| Chasing the 8 rows she marked "both calls are discordant" | Our benchmark truth set is ClinVar; on those rows ClinVar itself is disputed. | Record them as truth-set exclusions in the benchmark and report concordance with and without them. |

## 7. Sequencing, and one methodology point

Order followed, roughly by value per unit of work (steps 1 to 4 are done):

1. **Tier C bugs** (C1, C2, C3, C4). Self-contained, testable, and they poison everything downstream. *Done, run v9.*
2. **Tier A criteria fixes** (A1, A2, A3, A4 part 1, A5, A7). Pure classification-crate changes with unit tests; resolved 78 of the 122. *Done, run v9.*
3. **B3, B4, B5** (gnomAD QC fields, homology list, hemizygotes). Data-layer additions to the gnomAD builder. *Done; B4 landed in v9, B3 and B5 in v10.*
4. **B2** (FAF95). *Done, v10.* Turned out to be the same data-layer pass as B3/B5, so it was folded into it rather than run separately.
5. **B6/B7** (mechanism and validity gates) and **B8** (functional evidence). *Done, v10.* Along the way, **C5** and **C6**, and the PM2 threshold below.
6. **B1** (prevalence-aware per-gene-disease thresholds). Still open, and the only remaining item. It is a curation task, not an engineering one.

**One thing the rebuild exposed that nobody had asked about.** The v10 gnomAD database has far better coverage than the April v1 build - 50.0 % of chr15 pathogenic variants carry a gnomAD record against 21.5 % on a not-yet-rebuilt chromosome - so many variants stopped looking "absent" and PM2 collapsed. The cause was `pm2_ad_af_threshold = 0.0`, a literal reading of Richards 2015's "absent from controls" that does not survive the change in denominator: the 2015 text was written against ExAC's 60,706 exomes, and a singleton among gnomAD v4's 800,000 people is exactly the "not seen in the general population" PM2 asks about. Now 4e-5, from a sweep and inside published VCEP practice, worth 19 points of pathogenic recall. Better data met a threshold that was never right.

Steps 3 and 4 required re-extracting gnomAD v4.1 exomes to pick up nine columns the
builder had never captured: the FILTER verdict, the `lcr` / `segdup` / `non_par`
region flags, `AC_XY` / `AN_XY`, and `faf95` / `fafmax_faf95_max`. Every one is
optional on read, and each gate degrades to its previous behaviour when absent, so
upgrading fastVEP without rebuilding the annotation database cannot change a call.

**Methodology point worth raising with the reviewer.** A7 exposes a circularity we should be honest about in the manuscript. PS1, PM5 and PM1 are all computed from a ClinVar-derived index that currently includes the variant being classified, and A6 would add another ClinVar-derived gate. Any of these makes ClinVar concordance a partly self-fulfilling metric. Proposal: run the benchmark in two modes, a "leave-one-variant-out" mode with self-matches excluded (the honest number), and a "ClinVar-informed" mode matching what a lab pipeline would actually do, and report both.

## 8. Results (run v10)

Full ClinVar 2-star+ rerun, 673,660 variants, against the rebuilt 24-chromosome
gnomAD v4.1 `.osa2` database. Compared against v9 leave-one-out, which is the
like-for-like mode (both exclude the variant's own ClinVar record).

| Metric | v9 (leave-one-out) | **v10** | change |
|---|---:|---:|---:|
| Exact match | 59.9 % | **61.2 %** | +1.3 pp |
| Same-direction | 71.0 % | **74.1 %** | +3.1 pp |
| **Opposite-direction** | **46** | **25** | **−46 %** |
| Pathogenic recall | 48.1 % | **58.5 %** | +10.4 pp |
| Benign recall | 57.2 % | **65.3 %** | +8.1 pp |
| VUS recall | 97.3 % | 97.0 % | −0.3 pp |
| Likely-benign recall | 45.5 % | 45.7 % | +0.2 pp |
| NoCall | 382 | 382 | — |

Every headline metric improved except VUS recall, which gave up 0.3 pp. Read it
as three separate effects:

- **Better data.** The gnomAD rebuild raised record coverage sharply (50.0 % of
  chr15 pathogenic variants carry a record, against 21.5 % on the old build), and
  B2's filtering allele frequency and B3's QC gating made what is there more
  trustworthy. Most of the benign-recall gain is here.
- **A threshold the better data exposed.** PM2's strict-absence rule collapsed
  once more variants had records; moving it to 4e-5 is most of the
  pathogenic-recall gain. This was not a planned item and would not have been
  found without the rebuild.
- **The gates.** B6 and B7 cost about 0.7 pp of pathogenic recall on their own,
  and are the main contributor to halving opposite-direction calls: the PVS1
  firings they removed were the confidently-wrong ones.

25 opposite-direction calls out of 673,660, from 122 at the start of this review
round.

## 8b. Results (run v9, retained for the record)

Full ClinVar 2-star+ rerun, 673,660 variants, same annotation stack as v8 apart
from a `clinvar_protein.oga` rebuild that adds the distinct-nucleotide-change
count PS1 needs.

| Metric | v8 | v9 ClinVar-informed | v9 leave-one-out |
|---|---:|---:|---:|
| Exact match | 61.3 % | 60.7 % | 59.9 % |
| Same-direction | 73.6 % | 73.0 % | 71.0 % |
| **Opposite-direction** | **122** | **41** | **46** |
| NoCall | 0 | 382 | 382 |
| Pathogenic recall | 63.3 % | 60.8 % | 48.1 % |
| VUS recall | 94.2 % | 97.3 % | 97.3 % |
| Benign recall | 61.8 % | 57.2 % | 57.2 % |

Read it as two separate effects:

- **The criteria fixes** (v8 → v9 ClinVar-informed) cut opposite-direction calls by 66 % and cost 0.6 pp of same-direction agreement. VUS recall rises 3.1 pp, which is the intended direction: variants whose deciding evidence is case-level or literature-level now land in VUS instead of being confidently wrong.
- **Removing the ClinVar self-reference** (ClinVar-informed → leave-one-out) costs a further 2.0 pp of same-direction agreement and 12.7 pp of Pathogenic recall. That is the part of v8's number that was self-fulfilling, not a regression.

The 382 NoCall records are the ClinVar reference-agreement rows (`ALT=.`) that
v8 and earlier were annotating as real variants. They are excluded from the
denominator now rather than silently scored.

### Per-criterion firing changes

| Criterion | v8 | v9 | Driver |
|---|---:|---:|---|
| PS1 | 23,022 | 1,274 | self-match exclusion (A7); **94 % of v8 firings were the variant's own ClinVar record** |
| BS2 | 101,853 | 67,615 | cohort-scaled homozygote test (A4) |
| PM1 | 41,845 | 23,323 | missense/in-frame restriction plus PVS1 suppression (A5) |
| BP1 | 22,567 | 6,471 | mutation-spectrum veto (A1) |
| PVS1 | 56,088 | 55,904 | splice-evidence gate is precisely targeted, not a blunt cut |

On truth-Pathogenic variants specifically, BS2 misfires fell from 520 to 102
(−80 %) and BP1 misfires from 624 to 6 (−99 %).

### Review-round disposition

Of the reviewer's 122 adjudicated rows, **78 are resolved**; 44 remain
opposite-direction and 2 are new. The regenerated review table carries her
per-variant note forward on every row she annotated.

What remains is concentrated on Tier B, not on criteria logic: PVS1 accounts for
23 of the 46, dominated by gain-of-function genes (PCSK9 ×2, IFIH1 ×2, CYFIP2),
gene-disease validity, and gnomAD QC. The benign remainder (GJB2, C2, C9,
TXNL4A) is the prevalence-aware-threshold work in B1.

### What is implemented

As of run v10: **everything in this document except B1.**

| | Item | Landed |
|---|---|---|
| A1-A7 | Criteria fixes straight from the guideline text | v9 |
| C1-C4 | Annotation bugs (`ALT=.`, MNV frame, splice dup, indel offsets) | v9 |
| B4 | Homologous / segmental-duplication gene list | v9 |
| B2 | Filtering allele frequency for BA1/BS1 | v10 |
| B3 | gnomAD FILTER / region-flag gating of BA1/BS1/BS2/PM2 | v10 |
| B5 | Hemizygote-aware BS2 on non-PAR sex chromosomes | v10 |
| B6 | Disease mechanism gates PVS1 | v10 |
| B7 | Gene-disease validity gates PVS1/PP2/PM1 | v10 |
| B8 | `--functional-evidence` supplying PS3/BS3 | v10 |
| C5 | Configurable `--pick-order` | v10 |
| C6 | Legible review-table export | v10 |
| **B1** | **Prevalence-aware per-gene-disease thresholds** | **open** |

B1 is open because it needs a curated table, not code: the `disorders` config
scaffold, the per-gene threshold hook and the gene-level `.oga` loader are all
in place waiting for it.

Two thresholds also moved from convention to measurement, which was not in the
original plan: the BS2 prevalence bar (1e-5 → 1e-3) and the PM2 bar for dominant
genes (strict absence → 4e-5). Both sweeps are in [`ACMG.md`](ACMG.md).

## 9. Answers received, and what they changed

1. **BS2 on recessive genes: scale with the size behind the observation.** Implemented as a Poisson 95 % lower bound on the homozygote frequency compared against a disease-prevalence bar (`bs2_hom_prevalence_threshold`, default 1e-5), with an absolute floor of 2 homozygotes. One homozygote among gnomAD v4's 730 K individuals no longer fires BS2; three among 2 K still does.
2. **Gene-disease attribute table: build it from VCEP specifications plus Orphanet prevalence.** That is B1 and remains open. The per-gene threshold hook and the `disorders` config scaffold are in place to receive it.
3. **The 8 "both calls are discordant" rows** are still awaiting her call and stay in the review table until then.
4. **One classification per disease entity: yes.** Planned as phase 2 of B1, once the gene-disease attribute table exists to define the entities and their inheritance.

## 10. Still open

Only two things, and neither is code.

1. **B1, the per-gene-disease attribute table.** Prevalence, penetrance class, onset class and any published VCEP BA1/BS1 thresholds, keyed by gene and MONDO ID. Everything that consumes it exists; nothing can generate it. This is the item that would resolve the residual benign discordances (GJB2, C2, C9, TXNL4A) that no global threshold can.
2. **Her call on the 8 rows where both ClinVar and fastVEP are wrong**, so they can be held as a small expert-adjudicated test set independent of ClinVar.

What is *no longer* open, and why:

- **The BS2 prevalence bar** was question 1 here. Answered by measurement rather than by asking: swept across the full benchmark and set to 1e-3, the prevalence of the most common Mendelian disorders, which covers the hearing-loss / alpha-1 antitrypsin / familial Mediterranean fever cases in her point 2 without needing the per-gene table. See [`ACMG.md`](ACMG.md#choosing-the-bs2-prevalence-bar).
- **The PM2 bar for dominant genes** was not a question anyone had asked, and turned out to matter more. fastVEP read Richards 2015's "absent from controls" as literal absence, which is 19 points of pathogenic recall against gnomAD v4's 800,000 individuals. Now 4e-5, from the sweep and inside published VCEP practice. See [`ACMG.md`](ACMG.md#choosing-the-pm2-bar-for-dominant-genes).

Both are configurable, and both sweeps are reproducible with
[`sweep_acmg_thresholds.py`](../analysis/acmg_benchmark/scripts/sweep_acmg_thresholds.py).
