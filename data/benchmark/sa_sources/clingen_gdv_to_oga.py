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

Includes only `Definitive`, `Strong`, `Moderate` ClinGen classifications
(excludes `Limited`, `Disputed`, `Refuted`, `No Known Disease
Relationship`) to align with the SVI-recommended evidence threshold.

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

INCLUDED = {"Definitive", "Strong", "Moderate"}


def main():
    if len(sys.argv) != 3:
        sys.exit("usage: clingen_gdv_to_oga.py <clingen_csv> <out_tsv>")
    src, dst = sys.argv[1], sys.argv[2]

    by_gene = defaultdict(list)
    with open(src) as f:
        for row in csv.reader(f):
            if not row or len(row) < 7:
                continue
            gene = row[0].strip('"')
            if not gene or gene in (
                "GENE SYMBOL",
                "CLINGEN GENE DISEASE VALIDITY CURATIONS",
            ):
                continue
            if gene.startswith("+") or gene.startswith("FILE") or gene.startswith("WEBPAGE"):
                continue
            disease = row[2].strip('"')
            moi = row[4].strip('"')
            cls = row[6].strip('"')
            if cls not in INCLUDED:
                continue
            # genemap2 phenotype format: "<disease label>, <MIM>, <inheritance>"
            # We don't have OMIM MIM numbers, so use MONDO from row[3].
            mondo = row[3].strip('"') or "."
            phen = f"{disease} (ClinGen {cls}/{moi}, {mondo})"
            by_gene[gene].append(phen)

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
