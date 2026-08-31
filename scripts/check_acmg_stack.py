#!/usr/bin/env python3
"""Prove that every supplementary annotation the ACMG classifier reads is
actually reaching it.

Why this exists: an incomplete `--sa-dir` does not fail. `fastvep annotate
--acmg` runs, emits a classification for every variant, and silently drops the
evidence it has no data for. A missing REVEL database does not produce an error;
it produces VUS. So the useful question after a build is not "did sa-build
succeed" but "does each source answer a query", and that can only be measured by
annotating something and looking at what came back.

The probe variants are eight real ClinVar 2-star+ records on GRCh38, chosen to
span the consequence classes the criteria branch on: nonsense, frameshift-free
in-frame deletion, intronic, synonymous, missense, and a second nonsense. They
are hard-coded rather than sampled so that two runs of this script are
comparable.

One source is checked indirectly. An interval database only yields an
annotation when the position falls inside an interval, so no RepeatMasker key
appears for a variant that simply is not in a repeat, and its absence cannot be
distinguished from an unloaded file by looking at the annotation alone. BP3
carries that distinction instead: it reports `evaluated: false` when no repeat
source was loaded and `evaluated: true` when one was, so BP3's flag is the
RepeatMasker health check.

Usage:
    python3 scripts/check_acmg_stack.py \
        --gff3 Homo_sapiens.GRCh38.115.gff3 \
        --fasta Homo_sapiens.GRCh38.dna.primary_assembly.fa \
        --sa-dir sa_databases/

Exits non-zero when a source in the reference stack answered nothing, so this
can gate a deployment.
"""

import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile

# chrom, pos, id, ref, alt. GRCh38, no `chr` prefix so the IDs match an Ensembl
# GFF3 without renaming; fastVEP normalises contig names on the SA side.
PROBES = [
    ("1", "1050461", "AGRN_nonsense", "C", "T"),
    ("1", "2304066", "SKI_inframe_del", "TCGGAGG", "T"),
    ("2", "178527324", "TTN_intronic", "C", "CA"),
    ("7", "117480106", "CFTR_synonymous", "G", "T"),
    ("10", "87862043", "PTEN_missense", "A", "T"),
    ("14", "23412862", "MYH7_missense", "C", "T"),
    ("17", "43045685", "BRCA1_missense", "T", "A"),
    ("17", "43045709", "BRCA1_nonsense", "AG", "A"),
]

# (display name, where to look, key). "allele" keys sit on each transcript
# consequence; "gene" keys sit under the variant's `genes` map.
ALLELE_SOURCES = [
    ("gnomAD", "gnomad", "BA1 BS1 BS2 PM2"),
    ("ClinVar", "clinvar", "PM2 backstop, PP5/BP6, PM3/BP2"),
    ("REVEL", "revel", "PP3 BP4 (missense)"),
    ("SpliceAI", "spliceAI", "PVS1 splice, PP3 BP4 (splice), BP7"),
    ("PhyloP", "phylop", "BP7 conservation tier"),
]
GENE_SOURCES = [
    ("ClinGen GDV / OMIM", "omim", "PVS1 gate, PP2/PM1 gate"),
    ("gnomAD constraint", "gnomad_genes", "PVS1 pLI/LOEUF, PP2, BP1"),
    ("ClinVar protein index", "clinvar_protein", "PS1 PM1 PM5 BP1"),
]

# Criteria that never fire from a supplementary annotation, so a run without
# trio samples, curated functional evidence or case-control statistics leaves
# them unevaluated by design. Listing them here keeps them out of the report's
# "missing data" section, where they would look like a broken setup.
NON_SA_CRITERIA = {
    "PS2": "de novo confirmed; needs --proband/--mother/--father",
    "PM6": "de novo assumed; needs --proband/--mother/--father",
    "PP1": "segregation; needs family data",
    "BS4": "lack of segregation; needs family data",
    "PS3": "functional; needs --functional-evidence",
    "BS3": "functional; needs --functional-evidence",
    "PS4": "case-control; NotEvaluated by default per SVI",
    "PP4": "phenotype specificity; needs patient phenotype",
    "BP5": "alternate molecular cause; needs a case context",
    "BP2": "in-cis/in-trans; needs a companion variant at the same site",
    "PM3": "in-trans with pathogenic; needs genotypes or a companion variant",
}


def run_annotate(args, vcf_path, json_path):
    cmd = [
        args.fastvep, "annotate",
        "-i", vcf_path,
        "-o", json_path,
        "--gff3", args.gff3,
        "--sa-dir", args.sa_dir,
        "--acmg", "--pick",
        "--output-format", "json",
    ]
    if args.fasta:
        cmd += ["--fasta", args.fasta]
    if args.acmg_config:
        cmd += ["--acmg-config", args.acmg_config]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        sys.stderr.write(proc.stdout)
        sys.stderr.write(proc.stderr)
        sys.exit(f"fastvep annotate failed with exit code {proc.returncode}")
    return proc.stderr


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--gff3", required=True)
    p.add_argument("--fasta", help="strongly recommended: without it there are "
                                   "no coding consequences, so every "
                                   "missense-dependent criterion goes quiet")
    p.add_argument("--sa-dir", required=True)
    p.add_argument("--acmg-config")
    p.add_argument("--fastvep", default=shutil.which("fastvep") or "fastvep")
    p.add_argument("--keep-json", help="write the probe's JSON output here "
                                       "instead of a temporary file")
    args = p.parse_args()

    tmpdir = tempfile.mkdtemp(prefix="fastvep-acmg-check-")
    vcf_path = os.path.join(tmpdir, "probe.vcf")
    json_path = args.keep_json or os.path.join(tmpdir, "probe.json")
    with open(vcf_path, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for chrom, pos, vid, ref, alt in PROBES:
            f.write(f"{chrom}\t{pos}\t{vid}\t{ref}\t{alt}\t.\t.\t.\n")

    stderr = run_annotate(args, vcf_path, json_path)
    loaded = [ln for ln in stderr.splitlines() if "annotation source(s)" in ln]

    with open(json_path) as f:
        variants = json.load(f)

    allele_hits = {key: 0 for _, key, _ in ALLELE_SOURCES}
    gene_hits = {key: set() for _, key, _ in GENE_SOURCES}
    criteria = {}
    n_csq = 0
    n_missense = 0
    for v in variants:
        for _, ga in (v.get("genes") or {}).items():
            for key in gene_hits:
                if key in ga:
                    gene_hits[key].add(v.get("id", "?"))
        for tc in v.get("transcript_consequences") or []:
            n_csq += 1
            if "missense_variant" in (tc.get("consequence_terms") or []):
                n_missense += 1
            for key in allele_hits:
                if tc.get(key) is not None:
                    allele_hits[key] += 1
            acmg = tc.get("acmg") or {}
            for c in acmg.get("criteria") or []:
                seen, total = criteria.get(c["code"], (0, 0))
                criteria[c["code"]] = (seen + (1 if c.get("evaluated") else 0),
                                       total + 1)

    print(f"Probed {len(variants)} variants, {n_csq} picked consequences "
          f"({n_missense} missense)")
    for ln in loaded:
        print(f"  {ln}")
    print()

    missing = []
    rows = []
    for name, key, used_by in ALLELE_SOURCES:
        n = allele_hits[key]
        ok = n > 0
        rows.append((name, key, "OK" if ok else "MISSING",
                     f"{n}/{n_csq} annotated", used_by))
        if not ok:
            missing.append(name)
    for name, key, used_by in GENE_SOURCES:
        n = len(gene_hits[key])
        ok = n > 0
        rows.append((name, key, "OK" if ok else "MISSING",
                     f"{n}/{len(variants)} variants' genes", used_by))
        if not ok:
            missing.append(name)
    # RepeatMasker, via BP3's evaluated flag (see the module docstring).
    bp3_seen, bp3_total = criteria.get("BP3", (0, 0))
    ok = bp3_total > 0 and bp3_seen == bp3_total
    rows.append(("RepeatMasker", "(BP3 evaluated)", "OK" if ok else "MISSING",
                 f"BP3 evaluated {bp3_seen}/{bp3_total}", "BP3"))
    if not ok:
        missing.append("RepeatMasker")

    w = max(len(r[0]) for r in rows)
    print(f"{'source'.ljust(w)}  {'key':21}  {'status':8}  evidence")
    for name, key, status, detail, used_by in rows:
        print(f"{name.ljust(w)}  {key:21}  {status:8}  {detail:28}  {used_by}")
    print()

    unevaluated = sorted(c for c, (seen, _) in criteria.items() if seen == 0)
    expected = [c for c in unevaluated if c in NON_SA_CRITERIA]
    unexpected = [c for c in unevaluated if c not in NON_SA_CRITERIA]
    if expected:
        print("Never evaluated, by design (no supplementary annotation supplies these):")
        for c in expected:
            print(f"  {c:6} {NON_SA_CRITERIA[c]}")
        print()
    if unexpected:
        print("Never evaluated, NOT expected. Chase these:")
        for c in unexpected:
            print(f"  {c}")
        print()

    if missing:
        print(f"INCOMPLETE: {len(missing)} source(s) answered nothing: "
              f"{', '.join(missing)}")
        print("See docs/ACMG_SETUP.md for what each one unlocks.")
        return 1
    print("Complete: every source in the reference stack answered at least one query.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
