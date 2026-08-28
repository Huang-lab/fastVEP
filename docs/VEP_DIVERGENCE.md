# Where fastVEP differs from Ensembl VEP

fastVEP's consequence and HGVS output is a port of Ensembl VEP's own model, and the goal is to
agree with it. Agreement is the *evidence* that the port is faithful, not the objective: the
output is a prediction a clinician may act on, so where VEP is demonstrably wrong in a way that
changes a call, fastVEP is right instead and says so here.

This document is the complete list, in both directions. Every row count is measured against
**real Ensembl VEP 115.1** (Docker `ensemblorg/ensembl-vep:release_115.1`, `--gff` mode, the same
GRCh38 Ensembl 115 GFF3 and FASTA fastVEP reads, `--hgvs --symbol --canonical`) over a
6,600-variant stratified sample of the ClinVar 2-star+ set: **150,725 matched (variant, allele,
transcript) rows**, 72,632 of them with a coding consequence.

That sample is stratified to be hard - 54.5 % of its variants are not SNVs, against 7.2 % of the
ClinVar 2-star+ set it is drawn from - so the rates below are the rates on the shapes that
disagree, not the rates a typical callset sees. The genome-wide section at the end gives the other
end of that range.

## Current agreement

| Field | Scope | Rows disagreeing | Agreement |
|---|---|---:|---:|
| Consequence terms | coding rows | 105 | 99.86 % |
| `Amino_acids` | coding rows | 0 | **100 %** |
| `Codons` | coding rows | 0 | **100 %** |
| Splice terms | all rows | 0 | **100 %** |
| Whole consequence set | all rows | 117 | 99.92 % |
| `IMPACT` | all rows | 104 | 99.93 % |
| `HGVSc` | all rows | 597 | 99.60 % |
| `HGVSp` | all rows | 856 | 99.43 % |

Reproduce with `analysis/acmg_benchmark/` for the variant set and the harness described in
[Performance Benchmarks](../README.md#performance-benchmarks) for the VEP invocation.

---

## Part 1 - deliberate divergences

Six. Each is a case where Ensembl's output is either not valid HGVS, or contradicts what the
sequence says, and where matching it would degrade a clinical call.

### 1. A frameshift that introduces a premature stop is HIGH, not a moderate in-frame insertion

**101 consequence rows, 101 `IMPACT` rows and 101 `HGVSp` rows** - the largest divergence, and the
only one that moves an ACMG call.

`stop_retained` in `Utils/VariationEffect.pm` defers to `ref_eq_alt_sequence`
(`TranscriptVariationAllele.pm` release/115, l. 1321), whose first clause is:

```perl
return 1 if ( ($ref_pep eq substr($alt_pep, 0, 1) && $alt_pep =~ /\*/) || ... );
```

That asks whether the replacement keeps the residue it starts on and introduces a terminator
*anywhere* - not whether the annotated terminator survived. It therefore fires on any frameshift
whose new stop lands in the codon just after the insertion point. Both `frameshift` (l. 1435) and
`stop_gained` (l. 1208) return 0 when `stop_retained` holds, so the variant comes out
`inframe_insertion,stop_retained_variant`, MODERATE.

| Variant | ClinVar 2-star+ | Ensembl VEP 115.1 | fastVEP |
|---|---|---|---|
| BRCA1 `c.5030_5033dup` | Pathogenic (3-star) | `inframe_insertion,stop_retained_variant`, MODERATE | `stop_gained,frameshift_variant`, HIGH |
| BRCA1 `c.1499_1508dup` | Pathogenic (3-star) | same, MODERATE | HIGH |
| BRCA2 `c.3205_3206insAATTGCAGTCAATTAATAT` | Pathogenic (3-star) | same, MODERATE | HIGH |
| TP53 `c.895_919dup` | Pathogenic/Likely_pathogenic | same, MODERATE | HIGH |
| ITGB3 `c.122_125dup` | Pathogenic (3-star) | same, MODERATE | HIGH |

VEP's own `Amino_acids` for the first of these is `N/N*X`: the reference codon translates to Asn,
the edited window to Asn then a terminator then an incomplete codon. Nothing about that says the
annotated stop at residue 1863 survived. There are 34 such variants in the ClinVar 2-star+ set,
also in SDHB, USH2A, ATM, MSH6, BRIP1, GATA2, SPG11, TTN, EXT1, CLCN7, COL7A1, MANBA, TMEM127 and
TRNT1.

fastVEP reproduces Ensembl's other two `ref_eq_alt_sequence` clauses, which do test the
terminator: one asks whether it sits at the same residue on both sides, the other whether the
edited protein still matches the reference over the reference's own length and grows by fewer than
three residues past it. Only the first clause is refused.

**Consequence for ACMG:** 34 variants keep PVS1 that would otherwise lose it. Reversing this
divergence is a one-line change in `terms_for_window`
(`crates/fastvep-consequence/src/predictor.rs`) and would move those 34 out of P/LP.

### 2. A synonymous multi-residue window names the whole span

**137 `HGVSp` rows.**

A change spanning two codons that leaves both residues unchanged is `p.Leu346_Leu347=`. Ensembl
writes `p.LeuLeu346=` - two three-letter codes sharing one position - which is not a form HGVS
defines. fastVEP writes the span.

(fastVEP previously named one residue of the pair, `p.Leu347=`, which was well formed but arbitrary
and picked the second residue on a reverse-strand transcript. Neither tool agreed with the other
before this change either; the count is unchanged and the description is now correct.)

### 3. A delins names every reference residue it replaces

**18 `HGVSp` rows.**

For `Amino_acids` of `STHYHSLV/STNEW**V` at residues 1248-1255, translation stops at the first
terminator the change introduces, so residues 1250 through 1255 are all gone. fastVEP writes
`p.His1250_Val1255delinsAsnGluTrpTer`. Ensembl writes `p.His1250_Leu1254delinsAsnGluTrpTer`,
dropping the last reference residue while keeping the truncated replacement - five residues
replaced by four, where the window replaced six.

### 4. A `*` coordinate is not given a negative offset

**17 `HGVSc` rows.**

Ensembl writes `c.*-944del`. HGVS `*N` positions count *forward* from the terminator, so `*-944`
is self-contradictory, and the `*` here carries no number at all. fastVEP anchors the position to
the exonic base that follows it instead: `c.2319-944del`.

The anchor number itself is not independently verified - only that Ensembl's form is malformed.

### 5. A frameshift's terminator distance counts to the stop the frameshift creates

**373 `HGVSp` rows** - the largest `HGVSp` divergence of all.

`p.Leu151ProfsTer39` says the new reading frame runs 39 residues from the first changed one before
it hits a stop. That is a fact about the edited protein, so it can be computed rather than argued
about.

Ten of these rows were checked by rebuilding the CDS from the Ensembl 115 GFF3 and FASTA, applying
the variant and translating, independently of either tool. The reconstruction is corroborated by the
reference protein lengths it produced - TP53 393, FLCN 579, SLC17A5 495, MPV17 176, NRIP1 1158, all
matching UniProt.

| Gene | Variant | Computed | fastVEP | Ensembl VEP 115.1 |
|---|---|---:|---:|---:|
| MPV17 | `c.451dup` | **39** | 39 | 50 |
| NRIP1 | `c.3465_3468del` | **15** | 15 | 6 |
| GSS | `c.1295_1296del` | **62** | 62 | 49 |
| CNGB3 | `c.2085del` | **134** | 134 | 118 |
| DYM | `c.2043del` | **94** | 94 | 81 |
| TP53 | `p.Asp393Thrfs` | **29** | 29 | 89 |
| ARSA | `p.Arg498Profs` | **76** | 76 | 21 |
| SLC17A5 | `p.Thr448Profs` | **54** | 54 | 114 |
| FLCN | `p.Ala541Cysfs` | **61** | 61 | 60 |
| RHCE region | `p.Pro573Leufs` | **108** | 108 | 197 |

fastVEP matched on all ten, over differences of 1 to 89 and in both directions - Ensembl is shorter
on roughly two thirds of the rows and longer on the rest, and it is wrong either way.

In every one of the ten the new stop lies **past the reference protein's own terminator**, which is
the whole difficulty: the frameshift runs off the end of the annotated CDS and the stop that ends it
is in what was the 3' UTR. fastVEP translates a CDS that continues downstream for exactly that
reason. Ensembl's `_get_alternate_cds` builds its edited CDS from 3'-shifted coordinates while
taking consequences from unshifted ones, and the two do not line up.

**This entry previously sat in Part 2, as a fastVEP defect.** It was recorded there from reading
Ensembl's source rather than from computing the answer, and computing it reverses the verdict.

### 6. A frame with no stop in the transcript is reported as having none

**62 `HGVSp` rows.** Ensembl reports no stop on 43 of them and fastVEP on 19, so this is not one
tool being systematically more cautious.

Checked the same way as divergence 5 - rebuild the CDS, apply the variant, translate - on seven
rows, spanning both directions:

| Gene | Ensembl VEP 115.1 | fastVEP | Computed |
|---|---|---|---|
| ITGA2B `p.Ter1040Trpfs` | `Ter13` | `Ter?` | **no stop in 62 UTR codons** |
| CLN3 `p.Gly37Valfs` | `Ter43` | `Ter?` | **no stop** |
| BRCA1 `p.Cys44Leufs` | `Ter26` | `Ter?` | **no stop** |
| RUNX1 `p.Arg380Profs` | `Ter62` | `Ter?` | **no stop** |
| FBN1 `p.Arg283Serfs` | `Ter71` | `Ter?` | **no stop** |
| CARS2 `p.Val497Glyfs` | `Ter?` | `Ter101` | **101** |
| NTHL1 `p.Ter305Metfs` | `Ter?` | `Ter16` | **16** |

fastVEP was right on all seven. Where it writes `Ter?` the shifted frame really does run to the end
of the transcript without a stop, and where it names a distance the stop is there.

The ITGA2B reconstruction is worth stating because it is the strongest of the five: the reference
protein comes out at 1,039 residues starting `MARALCPL`, residue 1040 translates to Trp under the
edit - which is what *both* tools call it - and the remaining 62 codons of 3' UTR contain no stop in
that frame.

NTHL1 also settles a second disagreement in the same string: the terminator becomes Met, as fastVEP
writes, not Lys.

---

## Part 2 - gaps, where Ensembl is right

These are fastVEP's remaining defects, not disagreements. Listed with counts so the size of each
is on the record.

### HGVSc

| Gap | Rows | What it looks like |
|---|---:|---|
| The shift crosses the exon boundary and the two tools land on opposite sides of it | 344 | VEP `c.732dup`, fastVEP `c.731-1_731insA`; VEP `n.1134+1del`, fastVEP `n.1135del` |
| Both ends intronic, but a span reaching the boundary is still placed differently | 151 | VEP `c.956-5_957del`, fastVEP `c.956-6_956del` |
| Exonic 3'-shift off by one | 66 | VEP `c.2403_2406del`, fastVEP `c.2404_2407del` |
| fastVEP emits an HGVSc where VEP emits none | 19 | |

These four and divergence 4 account for all 597 `HGVSc` mismatches exactly:
344 + 151 + 66 + 19 + 17.

495 of the 597 - the first two rows - are the same unfinished piece of work: the 3'-shift is
bounded by the intron it starts in, so a span that should travel across the splice site stops at
it. The exonic shift and the intronic shift are each correct within their own territory.

**Fixed since this document was first written.** The largest gap it listed, *"intronic 3'-shift off
by one, and the anchor side chosen for a duplication deep in a long intron"* (1,238 rows), was
neither off by one nor about the anchor side. The duplication's position was derived by walking the
*unshifted* insertion point one base at a time while the next base repeated the current one, which
slides through a homopolymer and nothing else, so a `TG` insertion in a `TGTG…` repeat never moved
at all - up to 48 bases from where HGVS puts it, and only 16 % of the disagreements were a single
base. It also explains a pattern the row counts hid: because that walk stays where the VCF put the
variant, fastVEP's anchor was the genomically-left one on both strands, so it read as
under-shifting on forward-strand transcripts and over-shifting on reverse-strand ones.

### HGVSp

| Gap | Rows | What it looks like |
|---|---:|---|
| fastVEP emits an HGVSp where VEP emits none | 113 | mostly spans with one end outside the CDS |
| VEP emits one where fastVEP does not | 23 | VEP `p.Glu480del` |
| A delins naming the same replacement one residue further along | 15 | VEP `p.Thr454_Thr455delinsHisPro`, fastVEP `p.Thr455_Thr456delinsHisPro`. The protein-side twin of the exonic 3'-shift off by one above |
| Other | 10 | VEP `p.Ter44GlnextTer20`, fastVEP `p.Gln44ext*?`; VEP `p.Met4_?1` |
| `Amino_acids` of `X` from an ambiguous reference codon | 3 | VEP `p.Ter157=`, fastVEP `p.Xaa157=`. Which is right depends on whether the position is the terminator or genuinely unknown; unresolved |
| A frameshift differing in its first changed residue | 1 | |

These six and the five `HGVSp`-bearing divergences account for all 856 mismatches exactly:
373 + 137 + 113 + 101 + 62 + 23 + 18 + 15 + 10 + 3 + 1.

Note the shape of that list. 691 of the 856 are places fastVEP diverges on purpose and has been
checked against the sequence; 162 are gaps, and 113 of those are simply fastVEP naming a change
where Ensembl declines to. The `HGVSp` column is not 856 defects.

### Consequence terms

| Gap | Rows |
|---|---:|
| `3_prime_UTR_variant` / `5_prime_UTR_variant` missing on a delins that spans the CDS boundary | 6 |
| `non_coding_transcript_variant` where fastVEP says `non_coding_transcript_exon_variant`, for an insertion at an exon edge | 6 |
| `incomplete_terminal_codon_variant` missing on one boundary-spanning delins | 1 |

Three rows that were in this table are not gaps, and two of them are `IMPACT` differences. Each was
checked by rebuilding the transcript and translating it:

| Variant | Ensembl VEP 115.1 | fastVEP | What the sequence says |
|---|---|---|---|
| TSC1 `9:132923438 GCATGGTTATCAA>AC` | `stop_retained_variant`, LOW | `stop_lost`, **HIGH** | The change covers the last four bases of the CDS, the annotated `TGA` among them. The edited frame has no stop at residue 126 and runs to 353. The stop is lost. |
| ALDH3A2 `17:19663333 CCC>GGGCTAAAAGTACT` | `start_lost`, HIGH | `protein_altering_variant`, MODERATE | The transcript's first CDS record carries **phase 2**: the CDS begins mid-codon, 5' incomplete, and no initiator is annotated. There is no `ATG` to lose. |
| ST3GAL5 `2:85861216 TC>T` | `frameshift_variant,start_lost` | `frameshift_variant` | Same shape - the first CDS record in transcript order carries **phase 1**. |

The two `start_lost` rows are the same mistake: Ensembl treats the first codon of the CDS as the
initiator whether or not it is one, and on a transcript whose CDS is annotated as 5' incomplete it
is not. Declining to call `start_lost` there is the lower-impact answer, which is worth saying
plainly - it is not a case where the cautious reading and the correct one agree.

---

## What is *not* a divergence

Three things that look like one and are not.

**Consequences are not 3'-shifted.** `TranscriptVariation.pm` sets
`$self->{shifted} = (defined($args{-no_shift}) && !$args{-no_shift})`, so the flag is false unless
`--no_shift 0` is passed explicitly, and `get_all_OverlapConsequences` applies the shift offset to
a predicate's coordinates only when it is set. HGVS *is* shifted, by a separate call to
`_return_3prime`. fastVEP does the same. Shifting the span before the splice predicates makes
agreement much worse (365 mismatched rows becomes 2,262).

**Splice sites are matched against the bases that differ, not the variant's span.**
`_get_differing_regions` (`VariationFeatureOverlapAllele.pm`) XORs the two allele strings position
by position with the shorter padded, and groups the result into contiguous runs. That is both
narrower than the span - a matching interior is not tested - and wider, because the padding puts
every base past the reference allele's end into a region. fastVEP reproduces it, which is why the
splice terms agree on all 150,725 rows.

**`splice_region_variant` is decided by the last differing region, not the union.**
`_intron_effects` assigns rather than or-assigns it inside its region loop, so a change whose first
differing base is in the splice region and whose last is not comes out with no splice term.
Reproducing the tidier rule instead put a `splice_region_variant` on 63 rows that VEP calls plain
missense.

---

## Genome-wide, on an ordinary callset

The table at the top is measured on a sample stratified towards the hard shapes. The same
comparison over a systematic 1-in-200 sample of the GIAB HG002 WGS callset - 20,241 variants,
**122,317 matched rows**, real Ensembl VEP 115.1 through the same harness - is the other end of the
range:

| Field | ClinVar sample | HG002 genome-wide |
|---|---:|---:|
| Whole consequence set | 99.92 % | **100.000 %** (0 rows) |
| `IMPACT` | 99.93 % | **100.000 %** (0 rows) |
| `HGVSp` | 99.43 % | **99.985 %** (18 rows) |
| `HGVSc` | 99.60 % | **99.879 %** (148 rows) |

The three fields that carry a clinical call agree completely on a genome-wide callset, because the
divergences that remain need a coding indel and only 0.5 % of genome-wide transcript rows have a
coding consequence at all, against 45.2 % in the stratified sample.

`HGVSc` is the exception and does not dilute: before the intronic duplication fix it was **98.74 %**
genome-wide against 98.82 % on the ClinVar sample, and **87.0 %** on insertion rows, because deep
intronic repeats are where a WGS callset's indels live and that was exactly the failing case.

Most of the 148 rows left come from **multi-allelic VCF records** (`TAA>TA,T`), which are 1.18 % of
this callset. The shared first base is stripped across all alleles, but each allele is not then
trimmed against the reference on its own, so `ACAC>AC` stays a four-base replacement where it is a
two-base deletion.

That is a variant-parsing gap rather than an HGVS one, and it reaches further than the name. Written
three ways at the same site, the same 5-base deletion comes out:

| VCF | Ensembl VEP 115.1 | fastVEP |
|---|---|---|
| `GAAGAA>G` | `-`, `frameshift_variant`, `c.1364_1368del` | `-`, `frameshift_variant`, `c.1365_1369del` |
| `GAAGAAA>GA` | `-`, `frameshift_variant,splice_region_variant`, `c.1364_1368del` | `A`, `frameshift_variant,splice_region_variant`, `c.1361_1366delinsA` |

VEP trims the shared suffix and reports the allele as a clean deletion; fastVEP carries the extra
base, which widens the span by one, turns the `del` into a `delins`, and reports `A` where VEP
reports `-`. A comparison keyed on the allele does not even line those rows up.

Neither dataset here contains a *single*-alt record needing that trim - both are already
parsimonious - so the gap is reachable only through multi-allelic sites and through callers that
emit non-minimal records. Fixing it means carrying a position per allele rather than one per site,
which the variant representation does not do today. It is not fixed here.

## The earlier divergence, still standing

Before this list existed, one divergence was documented in the README: a two-residue deletion at a
protein C-terminus where fastVEP writes `p.Glu560_Glu561del` and VEP writes `p.Glu559_Glu560del`.
The cause is the scan bound in VEP's `_shift_3prime`, which halts *n*-1 residues before the
terminus for an *n*-residue change. fastVEP keeps the 3'-maximal answer HGVS specifies. See
[issue #94](https://github.com/Huang-lab/fastVEP/issues/94).
