#!/usr/bin/env python3
"""Count the HGVSp cells that differ between a VEP tab file and a fastVEP tab file by the rule that explains them.

usage: classify_hgvsp_by_rule.py VEP.tab FASTVEP.tab [--scope all|complete|refseq|incomplete]

Rows are paired on (Uploaded_variation, Allele, Feature), and a cell that differs is put in the first
class below whose test it passes, so the classes are exclusive and add up to the cells that differ.
Each class is a section of docs/VEP_DIVERGENCE.md; "other" is what is left, and is the stop-codon
family of divergence 13 and the cells that nothing here explains. `classify_hgvsp.py` groups the same
cells by the shape of the two strings and is the way to look into "other".

--scope:
  all         every paired row
  complete    the rows of VEP's output whose FLAGS is `-`, which is the transcripts a cache calls complete
              (a GFF3 has no `cds_start_NF` or `cds_end_NF`, so a run over one is compared on these)
  refseq      `complete` and BAM_EDIT is not `OK` (a transcript VEP edits from an alignment has a sequence
              that no GFF3 holds)
  incomplete  the rows whose FLAGS is not `-`
"""
import argparse
import collections
import re

UTR = {"5_prime_UTR_variant", "3_prime_UTR_variant"}
CODING = {
    "stop_gained", "frameshift_variant", "stop_lost", "start_lost", "inframe_insertion", "inframe_deletion",
    "missense_variant", "protein_altering_variant", "incomplete_terminal_codon_variant", "start_retained_variant",
    "stop_retained_variant", "synonymous_variant", "coding_sequence_variant",
}


def rows(path):
    index = None
    with open(path) as f:
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                index = {h: i for i, h in enumerate(line.rstrip("\n").lstrip("#").split("\t"))}
                continue
            yield line.rstrip("\n").split("\t"), index


def description(cell):
    return cell.split(":", 1)[1] if ":" in cell else cell


def classify(vep_cell, ours_cell, vep_terms, our_terms, vep_row, index):
    v, o = description(vep_cell), description(ours_cell)
    if "start_lost" in our_terms and o.endswith("1?"):
        if vep_cell == "-":
            return "15: p.Met1? where VEP writes nothing"
        if "_?" in v:
            return "15: p.Met1? where VEP writes a range or several residue names"
        if re.fullmatch(r"p\.[A-Za-z]{3}\d+\?", v):
            return "15: p.Met1? where VEP names a later residue"
    if "start_lost" in vep_terms and "start_lost" not in our_terms:
        return "14: no start_lost where the protein keeps its Met, or the whole coding sequence is behind the edit"
    if "start_lost" in our_terms and "start_lost" not in vep_terms and o.endswith("1?"):
        return "14: start_lost and p.Met1? where VEP writes p.Met1del, for a whole ATG deleted on a transcript with no 5' UTR"
    if ours_cell == "-" and (our_terms & UTR) and not (our_terms & CODING) and vep_cell != "-":
        return "16: a UTR-only row has no protein notation"
    if vep_cell == "-" and ours_cell != "-" and (our_terms & CODING):
        return "17: a protein change where VEP writes none (the walk over the transcript)"
    if "HGVS_OFFSET" in index and vep_row[index["STRAND"]] == "-1" and vep_row[index["VARIANT_CLASS"]] == "insertion":
        shift = vep_row[index["HGVS_OFFSET"]]
        if shift not in ("-", "") and abs(int(shift)) > len(vep_row[index["Allele"]]):
            return "12: a reverse-strand insertion moved further than it is long"
    return "other"


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("vep")
    parser.add_argument("fastvep")
    parser.add_argument("--scope", choices=["all", "complete", "refseq", "incomplete"], default="all")
    args = parser.parse_args()

    ours = {}
    for row, index in rows(args.fastvep):
        key = (row[index["Uploaded_variation"]], row[index["Allele"]], row[index["Feature"]])
        ours[key] = (row[index["HGVSp"]], row[index["Consequence"]])

    classes, example = collections.Counter(), {}
    paired = differing = 0
    for row, index in rows(args.vep):
        flags = row[index["FLAGS"]]
        if args.scope in ("complete", "refseq") and flags != "-":
            continue
        if args.scope == "refseq" and row[index["BAM_EDIT"]] == "OK":
            continue
        if args.scope == "incomplete" and flags == "-":
            continue
        key = (row[index["Uploaded_variation"]], row[index["Allele"]], row[index["Feature"]])
        if key not in ours:
            continue
        paired += 1
        vep_cell, row_terms = row[index["HGVSp"]], row[index["Consequence"]]
        ours_cell, our_consequence = ours[key]
        if vep_cell == ours_cell:
            continue
        differing += 1
        name = classify(vep_cell, ours_cell, set(row_terms.split(",")), set(our_consequence.split(",")), row, index)
        classes[name] += 1
        example.setdefault(name, (key, vep_cell, ours_cell))

    print(f"{args.vep}: {paired} paired rows, {differing} HGVSp cells differ")
    for name, count in classes.most_common():
        key, vep_cell, ours_cell = example[name]
        print(f"{count:8d}  {name}")
        print(f"          e.g. {key[0]} {key[2]}: VEP {vep_cell}, fastVEP {ours_cell}")


if __name__ == "__main__":
    main()
