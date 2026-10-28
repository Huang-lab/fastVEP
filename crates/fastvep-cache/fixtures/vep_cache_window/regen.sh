#!/usr/bin/env bash
# Regenerates this directory: a variation cache of the shape VEP's own converter writes
# (`7/all_vars.gz` and `.csi`, bgzip then `tabix -C -s 1 -b 5 -e 5`) holding about ten thousand
# records in the densities a real one has, for the test of the reader's sliding window.
#
#   ./regen.sh
#
# Needs python3, bgzip and tabix (htslib).
set -euo pipefail
export LC_ALL=C
here="$(cd "$(dirname "$0")" && pwd)"
mkdir -p "$here/7"
python3 - "$here" <<'PY'
import random, sys
here = sys.argv[1]
rng = random.Random(139)
cols = "chr,variation_name,failed,somatic,start,end,allele_string,strand,minor_allele,minor_allele_freq,clin_sig,phenotype_or_disease,clin_sig_allele,pubmed,var_synonyms,AFR,AMR,EAS,EUR,SAS".split(",")
starts = []
# a dense region, a record every few bases, like a real chromosome
p = 100_000
for _ in range(6_700):
    p += rng.choice([1, 2, 3, 3, 4, 5, 7])
    starts.append(p)
# a hot site: three hundred known variants at one position
starts += [140_000] * 300
# a sparse region
starts += [200_000 + 5_000 * i for i in range(40)]
# a gap of a megabase, then a dense region again, with indels whose end is before their start
p = 1_500_000
for _ in range(3_000):
    p += rng.choice([1, 2, 3, 4])
    starts.append(p)
starts.sort()
lines = []
for i, s in enumerate(starts):
    kind = rng.random()
    if kind < 0.75:
        alleles, end = "C/T", s
    elif kind < 0.9:
        alleles, end = "-/A", s - 1          # an insertion: end before start
    else:
        alleles, end = "AGT/-", s + 2
    lines.append("\t".join([
        "7", f"rs{i}", "0", "0", str(s), str(end), alleles, "1", ".", ".", ".", "0", ".", ".", ".",
        ".", ".", ".", ".", ".",
    ]))
open(f"{here}/7/all_vars", "w").write("\n".join(lines) + "\n")
open(f"{here}/info.txt", "w").write(
    "# CACHE UPDATED 2026-04-08 00:00:00\nspecies\thomo_sapiens\nassembly\tGRCh38\n"
    "cache_region_size\t1000000\nvariation_cols\t" + ",".join(cols) + "\nvar_type\ttabix\n")
print(f"{len(lines)} records")
PY
bgzip -f "$here/7/all_vars"
tabix -f -C -s 1 -b 5 -e 5 "$here/7/all_vars.gz"
