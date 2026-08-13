#!/usr/bin/env python3
"""Convert the UCSC RepeatMasker track to a BED interval file for BP3.

BP3 is "in-frame deletions/insertions in a repetitive region without a known
function" (Richards 2015). fastVEP evaluates it from a positional `.osi`
interval database, and without one BP3 reports `evaluated: false` rather than
guessing - so the criterion has been inert for want of a data file.

Source
------
UCSC `rmsk.txt.gz` for hg38, the MySQL dump of the RepeatMasker track:

    https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz

Columns, 0-based, no header:

    0 bin          5 genoName    10 repName    13 repStart
    1 swScore      6 genoStart   11 repClass   14 repEnd
    2 milliDiv     7 genoEnd     12 repFamily  15 repLeft
    3 milliDel     8 genoLeft                  16 id
    4 milliIns     9 strand

`genoStart`/`genoEnd` are already half-open 0-based, which is BED's own
convention, so no coordinate arithmetic is needed - a fencepost error here
would silently shift every interval by one base.

Usage
-----
    python3 repeatmasker_to_bed.py rmsk.txt.gz > repeatmasker.bed
    fastvep sa-build --source custom_bed --name repeatmasker \\
        -i repeatmasker.bed -o data/benchmark/sa_db/repeatmasker

The `--name repeatmasker` matters: the classifier finds the track by looking
for a supplementary key containing "repeat", so a differently-named database is
loaded and then ignored.

Scope
-----
All repeat classes are emitted, which is the literal reading of "repetitive
region". The class is carried in the BED name column so a curator can see *what
kind* of repeat drove BP3 - an in-frame indel inside a simple tandem repeat and
one inside an exonized Alu are not equally good arguments for benignity, and
the distinction is invisible if the name is dropped.

Only primary-assembly contigs are kept. Alt haplotypes, unplaced and unlocalised
scaffolds carry repeat annotations too, but no variant callset this is used with
reports on them, and they triple the interval count for nothing.
"""

from __future__ import annotations

import argparse
import gzip
import re
import sys

# chr1-22, chrX, chrY, chrM. Anything with an underscore is an alt/random/Un
# scaffold; `_alt`, `_random` and `chrUn_*` all match that.
PRIMARY = re.compile(r"^chr(\d{1,2}|X|Y|M)$")

GENO_NAME, GENO_START, GENO_END = 5, 6, 7
REP_NAME, REP_CLASS, REP_FAMILY = 10, 11, 12
MIN_FIELDS = 13


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("rmsk", help="UCSC rmsk.txt.gz (or uncompressed rmsk.txt)")
    ap.add_argument("--keep-alt-contigs", action="store_true",
                    help="also emit alt/random/unplaced scaffolds")
    ap.add_argument("--strip-chr", action="store_true",
                    help="emit '1' rather than 'chr1'. Not normally needed - the "
                         ".osi reader resolves chr-prefixed and bare names to "
                         "each other - but available if a downstream tool does not")
    ap.add_argument("--class", dest="classes", default=None,
                    help="comma-separated repClass allowlist, e.g. "
                         "'Simple_repeat,Low_complexity'. Default: every class")
    args = ap.parse_args()

    wanted = None
    if args.classes:
        wanted = {c.strip() for c in args.classes.split(",") if c.strip()}

    opener = gzip.open if args.rmsk.endswith(".gz") else open
    kept = skipped_contig = skipped_class = malformed = 0

    with opener(args.rmsk, "rt") as fh:
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) <= MIN_FIELDS:
                malformed += 1
                continue
            chrom = cols[GENO_NAME]
            if not args.keep_alt_contigs and not PRIMARY.match(chrom):
                skipped_contig += 1
                continue
            rep_class = cols[REP_CLASS]
            if wanted is not None and rep_class not in wanted:
                skipped_class += 1
                continue
            try:
                start = int(cols[GENO_START])
                end = int(cols[GENO_END])
            except ValueError:
                malformed += 1
                continue
            # Degenerate intervals would index as empty and can only confuse a
            # downstream overlap test.
            if end <= start:
                malformed += 1
                continue
            if args.strip_chr:
                chrom = chrom[3:] if chrom.startswith("chr") else chrom
            name = f"{rep_class}/{cols[REP_FAMILY]}/{cols[REP_NAME]}"
            print(f"{chrom}\t{start}\t{end}\t{name}")
            kept += 1

    print(f"kept {kept:,} intervals; skipped {skipped_contig:,} on non-primary "
          f"contigs, {skipped_class:,} by class filter, {malformed:,} malformed",
          file=sys.stderr)
    if kept == 0:
        print("no intervals emitted - is this really a UCSC rmsk dump?", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
