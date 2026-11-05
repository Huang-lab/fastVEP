#!/usr/bin/env python3
"""Convert ClinGen Gene-Disease Validity CSV → genemap2.txt-style TSV
for the `fastvep sa-build --source omim` pipeline.

ClinGen Gene-Disease Validity (GDV) is the **preferred** disease-gene
source per ClinGen SVI / Abou Tayoun 2018: a multi-curator scored
rubric (Strehlow et al. 2024) that, unlike OMIM, explicitly excludes
disputed, refuted, and limited associations. fastvep's PVS1 evaluator
consults `omim.phenotypes` non-emptiness in `is_lof_intolerant_gene`
(`PVS1.rs:311`) — the json_key is named `omim` for historical reasons,
but the canonical content is ClinGen GDV.

OMIM `genemap2.txt` (registration-gated at omim.org) remains supported
as a legacy alternative. Both populate the same `.oga` schema.

Emits **every** ClinGen classification, including the weak and negative ones
(`Limited`, `Disputed`, `Refuted`, `No Known Disease Relationship`).

That is a change from the original behaviour, which kept only Definitive,
Strong and Moderate. Dropping the rest made the file unable to distinguish two
very different states that the classifier needs to tell apart:

  * a gene ClinGen curated and concluded has little or no evidence behind it -
    positive evidence against a gene-disease relationship, and a good reason
    not to apply PVS1/PP2/PM1; from
  * a gene ClinGen has simply not reached yet - which says nothing at all.
    ClinGen has curated ~2,400 genes to Definitive/Strong/Moderate out of a
    genome, so absence is overwhelmingly the second case. SPAST, ABCB11 and
    LAMB3 are all absent, and all three cause disease.

The classifier reads the classification back out of the phenotype string and
decides; see `OmimData::has_established_relationship`.

Output is tab-separated with 13 columns matching the subset
`parse_omim_genemap` reads (col 5: MIM number / identifier, col 8:
Approved Gene Symbol, col 12: phenotype list joined with ";") — this
must match the real genemap2.txt column layout, since that's what
`parse_omim_genemap` parses (see issue #64: these two columns were
previously swapped in both the parser and this script).
"""

import csv
import sys
from collections import defaultdict

# Classifications that count as an established gene-disease relationship, per
# the SVI-recommended evidence threshold. Everything else is emitted too, but
# marked, so the classifier can tell "curated as weak" from "not curated".
ESTABLISHED = {"Definitive", "Strong", "Moderate"}

# The columns this reads, by the names the export gives them in its header row.
COLUMNS = ("GENE SYMBOL", "DISEASE LABEL", "DISEASE ID (MONDO)", "MOI", "CLASSIFICATION")

# The modes of inheritance the classifier knows (`ModeOfInheritance::from_code`). Another code, or none, is
# written as it stands and read as no inheritance information, so it is said here, at the build.
KNOWN_MOI = {"AD", "AR", "XL", "SD", "MT", "UD"}


def main():
    if len(sys.argv) != 3:
        sys.exit("usage: clingen_gdv_to_oga.py <clingen_csv> <out_tsv>")
    src, dst = sys.argv[1], sys.argv[2]

    by_gene = defaultdict(list)
    unknown_moi = defaultdict(int)
    columns = None  # name -> index, from the header row
    with open(src) as f:
        for row in csv.reader(f):
            if not row:
                continue
            if row[0] == "GENE SYMBOL":
                # By name and not by position: a column that moved would put another column's value in the
                # tag, and a tag whose mode of inheritance is wrong is read without a word as a gene with none.
                missing = [c for c in COLUMNS if c not in row]
                if missing:
                    sys.exit(f"{src}: the header row has no column {missing}: {row}")
                columns = {c: row.index(c) for c in COLUMNS}
                continue
            gene = row[0].strip('"')
            if not gene or gene == "CLINGEN GENE DISEASE VALIDITY CURATIONS":
                continue
            if gene.startswith("+") or gene.startswith("FILE") or gene.startswith("WEBPAGE"):
                continue
            if columns is None:
                sys.exit(f"{src}: a curation ({gene}) comes before the header row, so its columns are not known")
            if len(row) <= max(columns.values()):
                continue
            disease = row[columns["DISEASE LABEL"]].strip('"')
            moi = row[columns["MOI"]].strip('"')
            cls = row[columns["CLASSIFICATION"]].strip('"')
            if not cls:
                continue
            if moi not in KNOWN_MOI:
                unknown_moi[moi] += 1
            # genemap2 phenotype format: "<disease label>, <MIM>, <inheritance>"
            # We don't have OMIM MIM numbers, so use MONDO from the MONDO column.
            mondo = row[columns["DISEASE ID (MONDO)"]].strip('"') or "."
            phen = f"{disease} (ClinGen {cls}/{moi}, {mondo})"
            by_gene[gene].append(phen)
    if columns is None:
        sys.exit(f"{src}: no header row with GENE SYMBOL in it, so this is not a ClinGen Gene-Disease Validity export")
    if unknown_moi:
        print(
            f"warning: curations with a mode of inheritance the classifier does not know are read as having none: "
            f"{dict(sorted(unknown_moi.items(), key=lambda kv: kv[0]))}",
            file=sys.stderr,
        )

    with open(dst, "w") as out:
        # genemap2 has a comment header; parse_omim_genemap skips lines
        # starting with '#'.
        out.write(
            "# Synthetic genemap2 from ClinGen Gene-Disease Validity\n"
            "# col5=mim_proxy (always 0), col8=gene_symbol, col12=phenotypes\n"
        )
        for gene, phenotypes in sorted(by_gene.items()):
            cols = [""] * 13
            cols[5] = "0"  # MIM proxy; classifier only checks col 12 non-emptiness
            cols[8] = gene
            cols[12] = "; ".join(phenotypes)
            out.write("\t".join(cols) + "\n")

    print(f"wrote {dst}: {len(by_gene)} genes", file=sys.stderr)


if __name__ == "__main__":
    main()
