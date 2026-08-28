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

Four. Each is a case where Ensembl's output is either not valid HGVS or contradicts what the
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
| Frameshift terminator distance | 415 | VEP `p.Ser318PhefsTer181`, fastVEP `p.Ser318PhefsTer282`. Ensembl builds the edited CDS from *3'-shifted* coordinates for HGVS (`_get_alternate_cds`) while computing consequences from unshifted ones; fastVEP uses unshifted for both |
| fastVEP emits an HGVSp where VEP emits none | 113 | mostly spans with one end outside the CDS |
| VEP emits one where fastVEP does not | 23 | |
| First changed residue of a frameshift at the annotated terminator | 46 | VEP `p.Ter309LysfsTer14`, fastVEP `p.Ter309MetfsTer16` |
| `Amino_acids` of `X` from an ambiguous reference codon | 3 | VEP `p.Ter157=`, fastVEP `p.Xaa157=`. Which is right depends on whether the position is the terminator or genuinely unknown; unresolved |

The eight buckets above and the four divergences account for all 856 `HGVSp` mismatches exactly:
415 + 137 + 113 + 101 + 46 + 23 + 18 + 3.

### Consequence terms

| Gap | Rows |
|---|---:|
| `3_prime_UTR_variant` / `5_prime_UTR_variant` missing on a delins that spans the CDS boundary | 6 |
| `non_coding_transcript_variant` where fastVEP says `non_coding_transcript_exon_variant`, for an insertion at an exon edge | 6 |
| `start_lost` missing, or reported as `protein_altering_variant` | 2 |
| `stop_retained_variant` where fastVEP says `stop_lost` | 1 |
| `incomplete_terminal_codon_variant` missing on one boundary-spanning delins | 1 |

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

Of the 148 rows left, **108 come from multi-allelic VCF records** (`TAA>TA,T`), where the shared
first base is stripped across all alleles but each allele is not then trimmed on its own, so
`ACAC>AC` is described as a four-base deletion rather than a two-base one. That is a variant-parsing
gap rather than an HGVS one, it changes the consequence as well as the name, and it is not fixed
here.

## The earlier divergence, still standing

Before this list existed, one divergence was documented in the README: a two-residue deletion at a
protein C-terminus where fastVEP writes `p.Glu560_Glu561del` and VEP writes `p.Glu559_Glu560del`.
The cause is the scan bound in VEP's `_shift_3prime`, which halts *n*-1 residues before the
terminus for an *n*-residue change. fastVEP keeps the 3'-maximal answer HGVS specifies. See
[issue #94](https://github.com/Huang-lab/fastVEP/issues/94).
