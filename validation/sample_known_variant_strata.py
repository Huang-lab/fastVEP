#!/usr/bin/env python3
"""sample_known_variant_strata.py CACHE_DIR FASTA OUT.vcf [K]

A VCF whose records are drawn from the rarer shapes of the real VEP variation cache of chr21 and
chr22, K per stratum per chromosome (default 1500), so that the known-variant columns are exercised
on payloads a whole-genome callset seldom meets: allele-specific CLIN_SIG (with and without a deleted
base for its key), records with no alleles, failed and somatic ones, PubMed ids, three or more
alleles with gnomAD frequencies, fifteen or more populations, indels with gnomAD frequencies.

CACHE_DIR is a cache's version directory (`.../homo_sapiens/116_GRCh38`) with `{21,22}/all_vars.gz`
(read with bgzip), FASTA is the
primary assembly with a .fai (Ensembl contig names: 21, 22). Each record is turned into the VCF line a
caller would write: SNVs as they are, an insertion or deletion anchored on the base before it, a
record with no alleles as an SNV at its start. Deterministic: seed 139. Compare the two tools over
the result with validation/run_known_columns.sh."""
import random
import subprocess
import sys

cache, fasta, out = sys.argv[1:4]
K = int(sys.argv[4]) if len(sys.argv) > 4 else 1500
rng = random.Random(139)

COLS = ("chr,variation_name,failed,somatic,start,end,allele_string,strand,clin_sig,phenotype_or_disease,"
        "clin_sig_allele,clin_sig_ref_allele,clinical_impact,pubmed,var_synonyms,AF,AFR,AMR,EAS,EUR,SAS,"
        "gnomADe,gnomADe_AFR,gnomADe_AMR,gnomADe_ASJ,gnomADe_EAS,gnomADe_FIN,gnomADe_NFE,gnomADe_REMAINING,"
        "gnomADe_SAS,gnomADe_MID,gnomADg,gnomADg_AFR,gnomADg_AMI,gnomADg_AMR,gnomADg_ASJ,gnomADg_EAS,"
        "gnomADg_FIN,gnomADg_MID,gnomADg_NFE,gnomADg_REMAINING,gnomADg_SAS").split(",")
I = {c: i for i, c in enumerate(COLS)}

# the FASTA, by its index
fai = {}
for line in open(fasta + ".fai"):
    name, length, offset, bases, width = line.split()[:5]
    fai[name] = (int(length), int(offset), int(bases), int(width))
fh = open(fasta, "rb")


def fetch(chrom, start, end):
    """1-based inclusive."""
    length, offset, bases, width = fai[chrom]
    start0, end0 = start - 1, end
    first = offset + (start0 // bases) * width + start0 % bases
    last = offset + ((end0 - 1) // bases) * width + (end0 - 1) % bases
    fh.seek(first)
    return fh.read(last - first + 1).replace(b"\n", b"").decode().upper()


strata = {
    "clin_sig_allele": lambda f, a: f[I["clin_sig_allele"]] != ".",
    "clin_sig_allele_indel": lambda f, a: f[I["clin_sig_allele"]] != "." and "-" in a,
    "no_alleles": lambda f, a: "/" not in a,
    "failed": lambda f, a: f[I["failed"]] != ".",
    "somatic": lambda f, a: f[I["somatic"]] != "." and "/" in a,
    "pubmed": lambda f, a: f[I["pubmed"]] != ".",
    "multi_gnomadg": lambda f, a: a.count("/") >= 2 and f[I["gnomADg"]] != ".",
    "indel_gnomad": lambda f, a: "-" in a and (f[I["gnomADg"]] != "." or f[I["gnomADe"]] != "."),
    "af_multi": lambda f, a: a.count("/") >= 2 and f[I["AF"]] != ".",
    "many_pops": lambda f, a: sum(1 for c in COLS[I["AF"]:] if f[I[c]] != ".") >= 15,
    "phenotype": lambda f, a: f[I["phenotype_or_disease"]] != ".",
    "plain_snv_af": lambda f, a: len(a) == 3 and "/" in a and f[I["AF"]] != ".",
}

chosen = {}
for chrom in ("21", "22"):
    reservoirs = {s: [] for s in strata}
    seen = {s: 0 for s in strata}
    proc = subprocess.Popen(["bgzip", "-dc", f"{cache}/{chrom}/all_vars.gz"], stdout=subprocess.PIPE, text=True)
    for line in proc.stdout:
        f = line.rstrip("\n").split("\t")
        if len(f) != 42:
            continue
        a = f[I["allele_string"]]
        for s, pred in strata.items():
            if pred(f, a):
                seen[s] += 1
                r = reservoirs[s]
                if len(r) < K:
                    r.append(f)
                else:
                    j = rng.randrange(seen[s])
                    if j < K:
                        r[j] = f
    for s, r in reservoirs.items():
        print(chrom, s, "seen", seen[s], "kept", len(r), file=sys.stderr)
        for f in r:
            chosen[(chrom, int(f[I["start"]]), f[I["variation_name"]])] = f

other = {"A": "C", "C": "G", "G": "T", "T": "A"}
vcf = []
skipped = 0
for (chrom, start, name), f in sorted(chosen.items()):
    a = f[I["allele_string"]]
    end = int(f[I["end"]]) if f[I["end"]] != "." else start
    if "/" not in a:
        ref = fetch(chrom, start, start)
        alts = [other.get(ref, "A")]
        pos = start
    else:
        alleles = a.split("/")
        ref, alts = alleles[0], alleles[1:]
        if any(not x.replace("-", "").isalpha() for x in alleles):
            skipped += 1
            continue
        if "-" in alleles:
            if len(alleles) != 2:
                skipped += 1
                continue
            anchor = fetch(chrom, start - 1, start - 1)
            pos = start - 1
            ref = anchor + ("" if ref == "-" else ref)
            alts = [anchor + ("" if alts[0] == "-" else alts[0])]
        else:
            pos = start
            ref = fetch(chrom, start, start + len(ref) - 1) if len(ref) == 1 else ref
    vcf.append((chrom, pos, ref, alts))

vcf.sort(key=lambda r: (int(r[0]), r[1], r[2], r[3]))
with open(out, "w") as o:
    o.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    for n, (chrom, pos, ref, alts) in enumerate(vcf):
        o.write(f"{chrom}\t{pos}\t.\t{ref}\t{','.join(alts)}\t.\t.\t.\n")
print(len(chosen), "records chosen,", len(vcf), "written,", skipped, "skipped", file=sys.stderr)
