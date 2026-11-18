# Structural variants

fastVEP annotates a structural variant (SV) the way Ensembl VEP 116.3 does.
The same records are read as SVs, with the same class, the same interval, the same transcripts paired with them and the same terms for each pair.
`crates/fastvep-cli/tests/vep_tab_sv_golden.rs` holds that to VEP cell by cell, over the 17,284 rows that VEP writes for the 6,692 records of `crates/fastvep-cli/fixtures/vep_tab_sv` (`regen.sh` there runs VEP itself, from the pinned 116.3 image, in `--gff` mode).
This page is the rules, read from the Perl source of that release and checked against what it writes, and the places where fastVEP is not VEP.

Where the code is:

| Crate | Module | What it holds |
|---|---|---|
| `fastvep-core` | `sv` | the 18 classes (`SvClass`), `get_SO_term`, the alleles of a variant and the sites they sit at, the expansion of a tandem repeat into bases |
| `fastvep-io` | `sv_record` | how a VCF record becomes an SV: routing, the class, the interval, the skips |
| `fastvep-consequence` | `sv_predictor` | the terms of one variant against one transcript |
| `fastvep-annotate` | `sv` | which transcripts a variant is paired with, the rows of one that has none, and tandem repeats read as bases |

Both per-variant loops (`fastvep annotate` and the library that the web server runs) call `fastvep_annotate::sv`, so a structural variant has one answer whichever way it is asked.
`crates/fastvep-cli/tests/loop_parity.rs` holds them to it.

## 1. Which records are structural variants

A record is an SV when its REF and ALTs, joined, are not just `A`, `C`, `G` and `T`, and either `INFO/SVTYPE` is set or the ALTs joined with `,` hold a symbolic allele or a breakend (`/[<\[\]][^\*]+[>\]\[]|^\.\w+|\w+\.$/`).
So `<DEL>`, `A[chr1:2500[`, `.A` and `A.` are SVs, and `<*>` is not: the `[^\*]` needs a character other than `*` after the `<`.
A `T` beside a `<DEL>` in one record is an SV record, and its `T` is an allele of it (section 6).

## 2. The class

The class is read from the ALTs joined with `/`.
When `SVTYPE` is set and the ALTs do not start with `DEL`, `INS`, `DUP`, `INV`, `CNV` or `CN` and a number (case-sensitive, with or without the `<`), the class is read from `SVTYPE` instead.
`SvClass::from_type` is the port of `get_SO_term`, and the order of its steps decides:

1. a type repeated (`<INS>,<INS>`) is the type once;
2. `INS:ME` and `DEL:ME`, with `:ALU`, `:HERV`, `:LINE1` (or `:L1`) or `:SVA` after them, are the mobile element classes, and an element it does not know is the plain mobile element class;
3. `DEL` together with `DUP` is a copy number variation;
4. `DUP:TANDEM` is a tandem duplication;
5. `CNV:TR` is a tandem repeat;
6. `CN` and a digit: `<CN0>` is a deletion, `<CN2>` a duplication, any other number a copy number variation;
7. `CNV` is a copy number variation;
8. a bracket anywhere, or a `.` at either end, is a breakend (`chromosome_breakpoint`);
9. `<...>` is the name between the brackets up to its first `:`;
10. otherwise the type as it stands.

The name is then looked up, case-sensitively, among the 18 classes of Ensembl's `%SO_TERMS`: `insertion`, `mobile_element_insertion`, `Alu_insertion`, `HERV_insertion`, `LINE1_insertion`, `SVA_insertion`, `deletion`, `mobile_element_deletion`, `Alu_deletion`, `HERV_deletion`, `LINE1_deletion`, `SVA_deletion`, `tandem_repeat`, `tandem_duplication`, `duplication`, `copy_number_variation`, `inversion` and `chromosome_breakpoint`.
A name that is not among them (`<CPX>`, `<CTX>`, `<STR>`) is not supported, and the record is skipped (section 4).
`<del>` is not a deletion.

The class is `VARIANT_CLASS`, and it is the `Allele` column of every row but a breakend's.

## 3. The interval

The start is `POS + 1` for every symbolic class: the base at `POS` is the padding base.
The end is `start + |SVLEN| - 1` when `SVLEN` is set (the first value of several, and one that is not a number is zero), else `END`, else `start`.
`SVLEN` wins over `END`.
For a breakend `END` is dropped first, since a caller may put the mate's position there.

So an `<INS>` with `SVLEN=300` covers the 300 bases after `POS` (`--sv-insertion-as-point` makes every insertion the point instead; see section 11), an `<INS>` with `END=POS` is a zero-length interval between `POS` and `POS + 1` (Ensembl's insertion convention, `start = end + 1`), and an `<INS>` with neither is the one base `POS + 1`.
A record with no `SVLEN` or `END` is one base for every class but a deletion (section 4).

## 4. What is skipped

VEP skips these records with a warning, and so does fastVEP (the command line prints the first few and a total, the web server returns each as a warning):

- a class that is not supported, named in the warning;
- a deletion (any class named `deletion`) whose start is after its end, or that has neither `SVLEN` nor `END`: "deletion looks incomplete";
- a record whose start is more than one past its end;
- a variant, of any kind, whose end less its start is more than `--max-sv-size` (`--max_sv_size` is accepted too; 10,000,000 by default, `-1` for no limit).

## 5. Tandem repeats

A record whose class is a tandem one (`tandem_repeat` or `tandem_duplication`) and that has `RUS` and `RUC` or `RB` in its INFO is not annotated as a structural variant, when there is a reference sequence to read.
VEP reads the bases of its interval (the end is `start + max(SVLEN) - 1`), builds each allele by repeating the repeat unit, and annotates a small variant with the allele string `REF/R1/R2`.
`expand_tandem_repeat` does the same:

- `RN` is how many repeats each allele has (`.` for none), and when it is absent the alleles are counted by the `CNV:TR` in the type;
- a count is `RUC`, or `RB` divided by the unit's length; a unit of `.` is `N`;
- a list that runs out counts as none: `RN=2` over one unit and count is the first repeat alone;
- an allele longer than `--max-sv-size`, or a repeat that lists more than 1,024 alleles, leaves the record a structural variant (VEP has no limit on the number, and 2 MB of `RN` was 1.1 GB of rows), and `-1` is read as no repeat being short enough, as in VEP (`-1 || 5000` is `-1`, and nothing is not longer than it; the five tandem repeats of the golden's INFO shapes, run with `--max_sv_size -1`, give the same rows in both programs);
- an allele that is the reference, or one that repeats another, has no row;
- a record with no allele (`<DUP:TANDEM>` with the INFO and no `RN`, which counts only `CNV:TR`) has no row at all;
- one allele of another length than the reference is trimmed of the bases it shares with it at either end, and several alleles are not trimmed.

Without a reference, or without `RUS`, the record stays a structural variant of its class (`tandem_repeat` or `tandem_duplication`).

## 6. Alleles and breakends

A record that is not a breakend has one allele for each of its ALTs, numbered from 1, and each is written with the class as its `Allele`.
A record with two ALTs of one class has two identical rows, as in VEP.

A breakend is a bracket form (`A[chr1:2500[`, `]chr1:100]T`), or a symbolic `<BND>` with `CHR2` and `END2` (or Manta's `END`), which VEP rewrites to the bracket form.
It has the alleles of its own position and of each mate.
Its own is written `A.`, or `.A` when any of the mates' ALT pieces starts with a bracket, and numbered 0, and a mate is written as the ALT piece it came from.
A single breakend (`A.`, `.A`) has no mate and one allele, written `chromosome_breakpoint`.
A bare `<BND>`, with no `CHR2`, has no mate either, but its allele string is `A/<BND>` (the reference and the ALT), which VEP splits into two alleles, so it has two identical rows of that term for each transcript (the golden holds both).
A mate is on the contig it names, which is matched by the same names as the reference's (`chr1` and `1` are one).

## 7. Which transcripts a variant is paired with

A transcript is paired with a variant when the variant's interval overlaps the transcript widened by `--distance` on each side (the upstream side is its 5' side), or when a mate's position does, on the mate's contig.
Rows come out in the order of the transcripts' stable IDs.

A transcript has a row for each allele that is close to it: an allele that has no position of its own always is, and a breakend's own position or a mate's is only when it is on the transcript's contig and within 5,000 bases of it, which is a constant of Ensembl's (`MAX_DISTANCE_FROM_TRANSCRIPT`) and not `--distance`.
A transcript that is paired and has no allele close to it has no row (a large `--distance` pairs a transcript that no breakend is within 5,000 bases of).
A row that no term holds for is the `intergenic_variant` that nothing else displaced, with the transcript as its Feature: a mate within `--distance` of a transcript that is not beside it in the sense of upstream and downstream (which are asked of the variant's own interval) is one.

A variant with no transcript paired with it has one `intergenic_variant` row for each allele, a breakend's own and each mate's, with no feature.

## 8. The terms

For each (transcript, allele) the terms are Ensembl's `get_all_OverlapConsequences`: the definitions are tried in order of tier, the ones that a flag of the pair rules out are skipped, and the loop stops at the first tier above the lowest of tier 1 or 2 that matched.
A variant that completely covers a transcript has the flags of a complete overlap and nothing else, so only the terms of tier 1, the feature terms, and `coding_transcript_variant`, `non_coding_transcript_variant` and `NMD_transcript_variant` can apply to it.
A pair that no definition holds for is `intergenic_variant`.
The terms that are of a `VariationFeature` (splice, `stop_gained`, `missense_variant`, `synonymous_variant` and the like) never apply to an SV.

In rank order, with the rule for an SV as `sv_predictor.rs` has it:

| Rank | Term | Rule |
|---|---|---|
| 1 | `transcript_ablation` | covers the transcript completely, and is a deletion class |
| 5 | `frameshift_variant` | a deletion inside one exon, whose length (end less start, plus one) is not a multiple of 3; needs the CDS to be hit |
| 6 | `stop_lost` | a deletion over the stop codon's bases, unless it starts in a last codon that is partial |
| 7 | `start_lost` | both ends of the variant are in an exon, the span between their cDNA coordinates has a base of the start codon, and so has its genomic interval; the transcript's start must be known (no `cds_start_NF`); there is no test of the class, so an inversion, a duplication or an insertion has it too |
| 8 | `transcript_amplification` | covers the transcript completely, and is a duplication, a tandem duplication or a tandem repeat |
| 9 | `feature_elongation` | a duplication, tandem repeat or insertion class, inside the transcript, with a base in the cDNA |
| 10 | `feature_truncation` | a breakend: its own position (or the mate's, for a mate's allele) is in the transcript; a deletion: in the cDNA, and either a partial overlap that reaches past an end or inside the transcript |
| 12 | `inframe_deletion` | a deletion inside the coding part of one exon, a multiple of 3 long |
| 18 | `splice_polypyrimidine_tract_variant` | the intron effects of its whole span |
| 20 | `start_retained_variant` | the cDNA test of `start_lost` alone, so it is written beside it |
| 23 | `coding_sequence_variant` | not a complete overlap, a base in the CDS, and neither `frameshift_variant` nor `inframe_deletion` |
| 25, 26 | `5_prime_UTR_variant`, `3_prime_UTR_variant` | a base of the UTR, in the cDNA |
| 27 | `non_coding_transcript_exon_variant` | a transcript without a protein that the variant does not completely cover, a base in an exon |
| 28 | `intron_variant` | the intron effects of its whole span |
| 29 | `NMD_transcript_variant` | the transcript's biotype is `nonsense_mediated_decay`, and the variant overlaps it |
| 30 | `non_coding_transcript_variant` | a transcript with no protein that the variant overlaps, and no rank 27 (it has no base in an exon, or it covers the whole transcript) |
| 31 | `coding_transcript_variant` | a protein-coding transcript that the variant completely covers |
| 32, 33 | `upstream_gene_variant`, `downstream_gene_variant` | beside the transcript, within `--distance`, and not overlapping it |

`inframe_insertion` is never given to an SV (Ensembl's code says it cannot know what is inserted).

"In the cDNA" is Ensembl's `within_cdna`: the exon mapper has a coordinate for the variant that is inside the transcript's cDNA.
An insertion is the pair of bases it sits between, and the mapper makes a coordinate of the one that is in an exon, `(c + 1, c)`, so an insertion on the outside of the first or last base of a transcript has `(1, 0)` or `(L + 1, L)` and is not in the cDNA.
It has `upstream_gene_variant` or `downstream_gene_variant` at a distance of 0, and no `feature_elongation`.

The `IMPACT` is that of the term of lowest rank.
In release 116.3 `feature_elongation` and `feature_truncation` are `HIGH`, as are `transcript_ablation` and `transcript_amplification`.

## 9. The columns of a row

- `Allele` is the class, or the breakend's string.
- `cDNA_position`, `CDS_position`, `Protein_position`, `EXON` and `INTRON` are the cells of the exon mapper for the whole interval: a position when the variant has an exon (and the CDS and protein when it has a coding part as well), a `?` for an end that is not in an exon, and nothing for a variant that covers the transcript.
- `DISTANCE` is written only for `upstream_gene_variant` and `downstream_gene_variant`.
- `REF_ALLELE`, `UPLOADED_ALLELE`, `Existing_variation` and the HGVS columns are `-`.
- The allele of a tandem repeat that was read as bases is a small variant's, and its columns are a small variant's.

## 10. Where fastVEP is not VEP

**A breakend without a mate inside a transcript is `feature_truncation`.**
Run as it is by default, VEP writes the term for `A.`, `.A` and a bare `<BND>` that lies in a transcript.
Run with `--no_stats` (or with a `--pick` that has to compare the rows), it does not, and the same record has `coding_sequence_variant` alone.
The cause, as the source reads: `OutputFactory.pm` reads `$vfoa->{breakend}->{string}`, which makes an empty `{breakend}` on an allele that has none, and `feature_truncation` asks `within_feature` about that breakend; where the terms were worked out before that line (the statistics do it, and so does a pick) the question was asked about the variant's own position, which is in the transcript.
The default run is the one that is the code's intent, so it is what fastVEP writes.
Of the 104 rows of 35 breakend records (101 distinct tuples of the strict key), 49 rows (47 tuples) differ between the two runs, all of them the addition of this term (and the `HIGH` that comes with it).
`regen.sh` runs VEP without `--no_stats` for this reason.

**`<*>` and `<NON_REF>` give no row.**
VEP annotates `<*>` as a small variant whose allele is called `<*>`, and skips the whole record of `<NON_REF>`, a real ALT beside it included (`A G,<NON_REF>` has no row for `G`).
fastVEP gives neither of them a row, and annotates the real ALT beside them.
This is the same as for a gVCF's reference block.
A record whose only ALT is a `.` or a `<*>` beside an `SVTYPE` has no row either, where VEP writes the rows of the class that the `SVTYPE` names (`G . SVTYPE=DEL;END=1300` gave VEP's deletion rows and none here, and `G A,<NON_REF> SVTYPE=DEL;END=1300` gave its deletion rows where fastVEP annotates the SNV; both run against VEP 116.3 in the review of this branch, not repeated here).

**Rows that VEP does not write the same way twice are written one way.**
For a breakend whose mate is on another contig, which of its rows VEP writes varies from run to run: three runs over the 1,400 records of the breakend file of the vep-rs set gave 68,938, 68,846 and 68,880 rows, the 6 records that differ are all breakends of this kind, and fastVEP scores 0.9816, 0.9810 and 0.9812 against them.
So that variation is not what keeps the file from 1.0000.
Of the 2,575 tuples that differ from the run of the table, 903 are on 40 single breakends, where VEP pairs a transcript of another chromosome (the next paragraph), and 1,672 are rows for the mates of breakends on another chromosome that only fastVEP writes, against each of the three runs.
The contigs of that file are spelled one way throughout, so it is not the spelling below, and no cause was found.
A lead, not shown for this file: the transcripts that VEP reads for a stretch of records come from the span of the records read with them (`filter_features_by_min_max` of `AnnotationSource.pm`), and a differential run in the review of this branch gave a mate on another contig two `intergenic_variant` rows when its record was annotated alone and none beside a neighbour on the other contig (not repeated here).
fastVEP writes a fixed set, with the mate's contig matched by name.
The fixture leaves such records out.

**A transcript of another contig is never paired.**
In `--gff` mode a transcript is a bare record when `AnnotationType/Transcript.pm` asks which variants overlap it, and has no slice until `lazy_load_transcript` makes one a few lines later, so the contig it asks about is undefined, and `_chr_matches_for_vf_overlap` of `InputBuffer.pm` answers yes to anything: a variant can be paired with a transcript on another contig by the numbers of its interval.
It happens: of the 8,158 records of a gnomAD-SV v4.1 chromosome 21 sample, 905 have rows of VEP on a transcript of another chromosome (11,930 of its 109,318 rows).
fastVEP matches the contig first.

**A contig that the reference does not have.**
VEP skips a record on a contig that its annotation sources and synonyms do not name, and says so (`Parser.pm`: "Chromosome ... not found in annotation sources or synonyms"), for a structural variant as for a small one.
fastVEP warns once for each such contig and goes on with the record, which is an `intergenic_variant` row, because a database of the run can still say something about it (`annotate.rs`).

**A mate's contig is matched with or without `chr`.**
VEP decides which alleles of a breakend a transcript gets a row for with `_close_to_feature`, which compares contig names with `_compare_seq_region_names` (`VariationEffect.pm`): the names are equal when they are the same letter for letter, case aside, or when a `--synonyms` file pairs them.
So `chr21` and `21` are two contigs to it, and fastVEP takes them for one, as it does for its reference (`chr1` and `1`).
gnomAD-SV spells the contig `21` in its contig column and `chr21` in `CHR2`, and for such a breakend VEP's rows depend on how the GFF3 spells the contig: against a GFF3 that names it `21` they are those of the record's own position and none of its mate's, and against one that names it `chr21` the rows on transcripts are those of its mate and none of its own position's.
fastVEP writes the rows of both.
Of the 8,158 records of the gnomAD-SV sample, 572 have mate rows that only fastVEP writes, and every one is a symbolic `<BND>` whose `CHR2` is spelled with `chr` in a file whose contig column says `21`: 463 with `CHR2=chr21`, where the mate is on the same chromosome, and 109 with a mate on `chrX`, `chr22` or `chrY`.

To see that this is the whole of it, the sample is annotated with the contig written one way throughout: the contig column of the VCF `chr21`, and a GFF3 and a FASTA of chromosome 21 that name it `chr21` (`validation/sv_one_spelling.sh` does this, and prints the score).
VEP's rows (104,823, the same in three runs) and fastVEP's (104,966) then have a strict F1 of 0.9993, 104,816 tuples in common.
The 150 only fastVEP has are 136 on `unconfirmed_transcript` features, a type VEP's reader ignores (its log says so: "Ignoring 'unconfirmed_transcript' feature_type ... This feature_type is not supported in VEP", and see [VEP_DIVERGENCE.md](VEP_DIVERGENCE.md)), and 14 intergenic rows of 7 breakends whose mate is on a chromosome that this reference does not have (`chr22` and `chrY`), for which VEP writes no row at all (not traced to a line of its code); the 7 only VEP has are the `intergenic_variant` rows of 7 deletions and duplications whose only neighbours are such transcripts.

**Not written: `--overlaps`.**
VEP's `--overlaps` (`OverlapBP` and `OverlapPC`, and `bp_overlap` and `percentage_overlap` in JSON) is not implemented for any variant.
Regulatory features have no rows either, since fastVEP has no regulatory build.

**`--acmg` gives a structural variant no classification.**
The criteria are Richards 2015's, for a sequence change, and a structural variant has the terms of one in its rows (`frameshift_variant` for a deletion inside an exon, `start_lost` for an inversion over a start codon), which PVS1 would read as what they are for a base.
The criteria for copy number variants are another document's (Riggs 2020), and are not implemented.

**A tabix GFF3 is read for what a variant reaches.**
A run over a bgzipped, indexed GFF3 loads only the regions of its VCF's variants, and for a structural variant that is the interval it covers and the position of each mate (on the mate's contig), and not `POS +- distance` alone.

**JSON positions.**
VEP's JSON writer makes `cdna_start` and `cdna_end` (and the CDS and protein pairs) by splitting the cell at its first `-` and takes the end only when it is digits or `?` (`OutputFactory/JSON.pm`), so a range that has its total after it (`32-51/542`, with `--total_length`) has its start for its end.
That is so for every variant, and a structural variant's rows have totals like any other (`VEP_DIVERGENCE.md`, Part 3, `--total_length`).
fastVEP writes no totals in JSON and the end that the tab has; without `--total_length` the JSON of VEP agrees with its tab.
The VCF `CSQ` and the tab agree with VEP's.

**Five cells of the fixture.**
A tandem repeat with more than one allele is read as a small variant with several ALTs, and fastVEP describes each allele of such a record as it would be alone, where VEP 116 does not ([VEP_DIVERGENCE.md](VEP_DIVERGENCE.md), Part 1, section 11).
So `HGVSc` of an insertion into the repeat beside a deletion allele is the `dup` of the 3' rule where VEP writes the `ins` where it was put, VEP writes no `HGVS_OFFSET` for the three alleles that fastVEP moves, and VEP's `HGVSg` of one `delins` starts six bases after where the change does.
`KNOWN` in the golden test names the five cells, and asserts they still differ.

**`--pick` and `--flag-pick`.**
The transcript a pick chooses for an SV is VEP's, except where all eight criteria tie, where VEP goes on to the length of the transcript and fastVEP to its ID (the gap of [VEP_DIVERGENCE.md](VEP_DIVERGENCE.md), Part 2).
Of the 344 records of the fixture whose `PICK` column differs under `--flag-pick`, 233 have another transcript flagged, and every one of them ties on canonical, biotype and the rank of its terms in a world that has no MANE, TSL, APPRIS or CCDS annotation; 85 have the same transcript with more of its alleles flagged (divergence 10), 21 are variants that no transcript is paired with, whose `intergenic_variant` rows are flagged by every flagging pick (below), and 5 have rows that a key cannot tell apart (a record with several ALTs of one class, the bare `<BND>`).
A variant that no transcript is paired with has its `intergenic_variant` rows flagged by every flagging pick.
VEP flags all of them for `--flag_pick_allele_gene` (its per-gene pick returns every row that is not a transcript's), as fastVEP does, and only the first for `--flag_pick`, and for an SV for `--flag_pick_allele` as well (a two-allele deletion/duplication record and a breakend with a mate, run against VEP 116.3).
fastVEP flags the row of every allele in those two cases, for the reasons of divergence 10 and of the paragraph below.
`--pick` keeps every allele of the transcript it picks (divergence 10), so a breakend has a row for its own allele and for each mate's, where VEP has one row.
For an SV VEP makes `--pick_allele` the same as `--pick`, `--pick_allele_gene` the same as `--per_gene`, and `--flag_pick_allele` the same as `--flag_pick` (`filter_StructuralVariationOverlapAlleles` of `OutputFactory.pm`): one row for the variant, or one for each gene, whatever its alleles.
fastVEP picks for each ALT of the record as it does for a small variant, so a record with several ALTs keeps a row for each of them, and a breakend with one ALT keeps one under `--pick-allele`, whose own position and mates share that ALT, as VEP does.
Which of the rows of one transcript is flagged by `--flag-pick-allele-gene` can differ: for one breakend with a mate, VEP flagged the row of the mate's allele and fastVEP the row of the breakend's own (found by the review of this branch, not repeated here).
The `FLAGS` column is fastVEP's (`cds_start_NF` and `cds_end_NF` read from the sequence), which a `--gff` run of VEP does not write, and so `start_lost` and `start_retained_variant` are not written for a transcript that has `cds_start_NF`, where VEP writes them ([VEP_DIVERGENCE.md](VEP_DIVERGENCE.md), Part 2).

## 11. VEP rules that are kept

These look wrong, or surprising, and fastVEP does what VEP does, because the rows are to be VEP's.

- An `<INS>` with an `SVLEN` is an interval of that many bases, so a long insertion can reach a UTR, an intron or a neighbouring gene that its insertion point does not.
  On the 1,030 insertions of a gnomAD-SV v4.1 chr21 sample (mobile elements included), reading the same records as the point between `POS` and `POS + 1` changes the rows of 92 of them (8.9 per cent): 47 (4.6 per cent) have no exon, UTR or coding term as a point (`coding_sequence_variant` of 25 records becomes 6, `non_coding_transcript_exon_variant` 67 becomes 32, and none gains one), and `feature_elongation`, which is `HIGH`, moves both ways: 29 records lose it and 7 gain it (72 have it as an interval and 50 as a point), since it is asked of a variant inside the transcript (section 8) and an interval that runs past an end of a transcript is not inside it, where its point can be.
  The point is what the VCF specification describes, and it is not VEP's reading, so `--sv-insertion-as-point` asks for it and it is not the default: without the option the rows are VEP's.
- `start_lost` has no test for the class, `start_retained_variant` is written beside it, and `start_lost` needs both ends of the variant to be in an exon.
- `frameshift_variant` and `inframe_deletion` are for a deletion inside one exon, by its length and not by the bases of the CDS it takes, so a deletion across two exons has neither.
- `feature_elongation` and `feature_truncation` are `HIGH`, and an insertion into an intron or a deletion inside one has neither.
- A breakend's mate can have `feature_truncation` and a flank term at once, since the first is asked of the mate's position and the flank of the variant's.
- A record with several ALTs of one class has identical rows, so its output grows with the square of the number of its ALTs (10,001 `<DEL>` ALTs gave 1.8 GB of rows in a fuzz of the program); the web server holds the ALT list of a record to 256.

## 12. How this is tested

- `tests/vep_tab_sv_golden.rs`: every row of the fixture, 23 columns cell by cell, and the order of the rows, against VEP.
  The fixture is 10 transcripts of every shape (both strands, with UTRs, with a coding sequence that fills its exons, with a partial last codon, with a frameshift intron, with a skipped exon, non-coding, NMD, miRNA, antisense), and the variants are every class laid over every place on every transcript where one end can be on another kind of sequence than the other, with the records that are about reading (`END` and `SVLEN`, the type from `SVTYPE`, several ALTs, breakends of every notation, tandem repeats, the skips).
- `validation/run_sv_world.py`: the same comparison on a world that nothing was fitted to.
  It runs `gen_world.py` for each seed it is given, runs the pinned VEP image and the binary over the world, and compares the rows, their order and the 23 columns as the golden does.
  Eight seeds besides the golden's own (11, 23, 37, 41, 59, 73, 97 and 131) are 139,860 rows and 3,216,780 cells of VEP 116.3: every row is in both tools and in the same order, and every cell is VEP's but for the ones of divergence 11 (5 to 7 for each seed, all of a tandem repeat with several alleles), and the golden's seed run again gives its 17,284 rows and the same 5 cells.
- `tests/loop_parity.rs`: the two loops write the same rows for it.
- The unit tests of `fastvep-core::sv`, `fastvep-io::sv_record` and `fastvep-consequence::sv_predictor`: the class of each notation, the interval of each record, the INFO values that Perl reads in its own way (`KEY=` is a value, a count of 0 is false, a length past the largest number is that number), the tandem repeats, and each term's rule on a transcript made for it.
- `tests/gff3_tabix_region_cache.rs`: a structural variant is read from a tabix GFF3 for the span it covers and for its mates, and a record that reaches past what the index addresses does not empty its contig; `crates/fastvep-web/tests/bounds.rs`: a record that lists more repeat alleles than a request may hold is refused.
- `the_ranks_and_impacts_are_those_of_release_116_3` in `fastvep-core`: every rank and impact of the 41 terms, from `Utils/Constants.pm`.

To regenerate the fixture after a change to the generator: `crates/fastvep-cli/fixtures/vep_tab_sv/regen.sh` (needs Docker, `bgzip` and `tabix`).

## 13. Speed

Measured on 2026-11-17 on an Apple M4, over the GRCh38 Ensembl 115 GFF3 and FASTA, with tab output: the median of 5 runs of fastVEP and 3 of Ensembl VEP 116.3 (Docker, one process, `--gff`, `--no_stats`), each after a run that is not timed (`benchmarks/run_sv_head_to_head.sh`, results in `benchmarks/results/sv_head_to_head/`).

| Input | Records | Ensembl VEP | fastVEP, 1 thread | fastVEP, 10 cores |
|---|---:|---:|---:|---:|
| gnomAD-SV v4.1, chromosome 21, a 1-in-4 sample | 8,158 | 20.9s | 1.25s (16.7x) | 0.70s (30.0x) |
| ClinVar structural variants, chromosome 21 | 3,060 | 14.5s | 0.94s (15.4x) | 0.39s (36.8x) |

These inputs are small, so most of a run of fastVEP is its transcript cache (0.86s of the 0.94s of the ClinVar sample on one thread), while the start-up of Ensembl VEP is 0.36s of its 14.5s and 20.9s and the rest is its records; the ratio is larger for a larger input, because the time of fastVEP grows slowly with it (the chr22 table of the README).
The rules change what a run over small variants writes in two places only, and what it costs by a little.
On every benchmark input of the repository (`benchmarks/run_identity_sweep.sh`) its output is the same bytes as master's and as a rebuild of the release candidate's but for two changes that were made on purpose: `--flag-pick*` flags the `intergenic_variant` rows that it left without a flag (7,001 of the 50,284 variants of chromosome 22 of HG002), as Ensembl VEP does, and the ranks of `splice_donor_5th_base_variant` and `splice_region_variant` are release 116.3's, which changes the `most_severe_consequence` of 3 of the 4,048,342 JSON records of HG002 (each a variant with the one term on a transcript and the other on another).
The fastest run of each build, over the rounds of 16 workloads, differs by -1.0% to +1.4% from one build to the other, and the medians by -7.0% to +2.9%; a laptop moves a round by about 3%, and by tens of percent when it stalls (`benchmarks/run_small_variant_regression.sh`, `benchmarks/results/small_variant_regression/`).
The work that a run does is told by the instructions that it retires, which no other load on the machine moves: on one thread the two builds differ by -0.39% to +0.23% in 6 of the 8 workloads, by +0.67% for JSON and by +1.86% for tab output, which also takes 1.4% longer in its fastest run; the cost of tab output is not in the lines of the tab writer that changed (restoring them as master has them leaves it; `tab_output_cost.tsv` of the results has the instructions of each build, on one thread over the 1-in-5 sample, the median of three runs).
A row of a variant carries the number of the structural variant's allele that it answers for (four bytes) and is 616 bytes, within the 632 that the loop was measured with, because the ACMG classification that only a run with `--acmg` fills is a pointer to it and not the value; `a_row_stays_within_the_size_the_loop_was_measured_with` holds it to the 632.
The release binary is 6.02 MB with the rules, 3.1% more than the 5.83 MB of the release candidate that the tables of the README were measured with.
