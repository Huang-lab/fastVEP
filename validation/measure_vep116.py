#!/usr/bin/env python3
"""fastVEP against Ensembl VEP 116.3, cell by cell: the counts of Section 3.3 of the manuscript and of Part 3 of docs/VEP_DIVERGENCE.md, in one command.

    validation/measure_vep116.py --bin target/release/fastvep --out benchmarks/results/v0.5.0 [--only hg002,clinvar,example,inframe,world,edges,lone]

What each input is, and how it is run (the commands of Part 3 of docs/VEP_DIVERGENCE.md):

    example   the 173 variants of the VEP example set (validation/human/vep_example_GRCh38.vcf), over the human gene model
    inframe   the 400 ClinVar in-frame deletions (validation/human/clinvar_inframe_deletions.vcf)
    hg002     a systematic 1-in-200 sample of the GIAB HG002 call set (20,241 variants; validation/sample_variants.py --every 200)
    clinvar   a systematic 1-in-5 sample of the ClinVar 2-star+ call set (134,732), twice: without the HGVS columns and with them
    world     a generated world of 3,000 records on both strands, 372 of them with several ALTs (the fixtures' gen_world.py --all, seed 139)
    edges     every shape of change within ten bases of every exon, CDS and transcript boundary of that world (gen_edges.py, every record)
    lone      its lone indels written with shared bases (gen_lone.py, 1,500 records)
    refseq    chromosomes 21 and 22 of NCBI's RefSeq GFF3 for GRCh38.p14, the contigs renamed 21 and 22 so that the GRCh38 FASTA answers for them, and the
              records of the HG002 sample on those two chromosomes; it needs `--refseq-gff3 GCF_000001405.40_GRCh38.p14_genomic.gff.gz` (NCBI's file, which
              no script of the repository fetches) and is run only when that is given. Its cells that differ are counted by column in
              vep116_refseq_columns.csv: a RefSeq GFF3 carries no biotype and no canonical flag, so on those two columns VEP's `--gff` reader writes its own
              guess and fastVEP, which reads a transcript as VEP's cache holds it, writes the cache's

Ensembl VEP runs in Docker from the image of release 116.3 (pinned by digest, as the fixtures are made), in `--gff` mode over a bgzip and tabix
copy of the same GFF3 that fastVEP reads, with `--tab --symbol --biotype --canonical --numbers --variant_class --show_ref_allele --hgvs --hgvsg`
(`--hgvs --hgvsg` left out of the second ClinVar run). `compare_tab.py` pairs the rows on (variant, allele, transcript) and compares every cell of the
columns both write; against the human gene model it leaves out `SYMBOL_SOURCE`, `HGNC_ID` and `FLAGS` (a GFF3 carries none of them, and fastVEP follows
VEP's own cache for them) and the two columns that name the annotation file, which is 27 columns with the HGVS columns and 23 without.

For the human inputs the cells that differ are then classified, so that "each has a documented cause" is something this script computes: an allele
of a record with several ALTs (divergence 11 of docs/VEP_DIVERGENCE.md), the `*` offset of an intronic base after the terminator (divergence 4), a
protein change that Ensembl writes none for (divergence 17), and anything else, which is reported as `unclassified` and should be none.

Writes vep116_agreement.csv (one row per input), vep116_hg002_columns.csv (the HG002 sample by column and by cause), vep116_run.tsv (the image, the
binary, the date) and the differing rows of the human inputs (vep116_diff/) into --out. A run of some of the inputs replaces their rows of
vep116_agreement.csv and leaves the others, and adds a block to vep116_run.tsv. Needs Docker, bgzip and tabix; standard library only.
"""

from __future__ import annotations

import argparse
import collections
import csv
import gzip
import hashlib
import re
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
HUMAN = (ROOT / "test_data" / "organisms" / "human").resolve()
CLINVAR = (ROOT / "data" / "benchmark" / "clinvar_2star.vcf").resolve()
FIXTURES = ROOT / "crates" / "fastvep-cli" / "fixtures"
# release_116.3, by digest, as the fixtures are made (docs/VEP_DIVERGENCE.md, Part 3)
IMAGE = "ensemblorg/ensembl-vep@sha256:2a46b0ba4ac544a87136aad92af92050ccbbf09098c1ce306d7e88ef7ed19aa3"
VEP_HGVS = ["--hgvs", "--hgvsg"]
VEP_OTHER = ["--symbol", "--biotype", "--canonical", "--numbers", "--variant_class", "--show_ref_allele"]
OUR_HGVS = ["--hgvs", "--hgvsg"]
OUR_OTHER = ["--symbol", "--biotype", "--canonical", "--numbers", "--variant-class", "--show-ref-allele"]
HUMAN_IGNORE = "SOURCE,Homo_sapiens.GRCh38.115.sorted.gff3.gz,SYMBOL_SOURCE,HGNC_ID,FLAGS"
WORLD_IGNORE = "SOURCE,genes.sorted.gff3.gz"
REFSEQ_IGNORE = "SOURCE,refseq.sorted.gff3.gz,SYMBOL_SOURCE,HGNC_ID,FLAGS"  # a GFF3 has none of the three, RefSeq's no more than Ensembl's
REFSEQ_CONTIGS = {"NC_000021.9": "21", "NC_000022.11": "22"}  # chromosomes 21 and 22 of GRCh38.p14, as the GRCh38 FASTA names them
ALL = ("example", "inframe", "hg002", "clinvar", "world", "edges", "lone", "refseq")
ORDER = ["example", "inframe_deletions", "hg002_1in200", "clinvar_1in5", "clinvar_1in5_hgvs", "world", "edges", "lone_indels", "refseq_gff3"]


def run(cmd: list[str], **kw) -> subprocess.CompletedProcess:
    r = subprocess.run(cmd, capture_output=True, text=True, **kw)
    if r.returncode != 0:
        sys.exit(f"FAILED (exit {r.returncode}): {' '.join(cmd)[:300]}\n{r.stderr[-800:]}")
    return r


def vep(work: Path, refs: dict[str, Path], vcf: str, gff: str, fasta: str, out: str, with_hgvs: bool) -> float:
    """One run of VEP in a fresh container: the wall-clock seconds."""
    mounts = ["-v", f"{work}:/data"] + [x for name, p in refs.items() for x in ("-v", f"{p}:/{name}:ro")]
    cmd = ["docker", "run", "--rm", "-e", "VEP_DATABASE=0", *mounts, IMAGE, "vep", "--input_file", f"/data/{vcf}", "--gff", gff, "--fasta", fasta,
           "--output_file", f"/data/{out}", "--tab", "--force_overwrite", "--no_stats", *VEP_OTHER, *(VEP_HGVS if with_hgvs else [])]
    t0 = time.perf_counter()
    r = subprocess.run(cmd, capture_output=True, text=True)
    wall = time.perf_counter() - t0
    if r.returncode != 0 or not (work / out).exists():
        sys.exit(f"VEP failed (exit {r.returncode}): {' '.join(cmd)[:300]}\n{r.stderr[-800:]}")
    return wall


def fastvep(binary: Path, work: Path, vcf: Path, gff3: Path, fasta: Path, out: str, with_hgvs: bool) -> float:
    t0 = time.perf_counter()
    run([str(binary), "annotate", "-i", str(vcf), "--gff3", str(gff3), "--fasta", str(fasta), "--output-format", "tab", "-o", str(work / out),
         *OUR_OTHER, *(OUR_HGVS if with_hgvs else [])])
    return time.perf_counter() - t0


def compare(work: Path, vep_tab: str, our_tab: str, ignore: str, diff_dir: Path | None) -> dict:
    """compare_tab.py over the two files, and what it prints, as numbers."""
    cmd = [sys.executable, str(HERE / "compare_tab.py"), str(work / vep_tab), str(work / our_tab), "--ignore", ignore]
    if diff_dir is not None:
        diff_dir.mkdir(parents=True, exist_ok=True)
        cmd += ["--out", str(diff_dir)]
    text = run(cmd).stdout
    m = re.search(r"== cells: (\d+) paired rows, (\d+) cells compared, (\d+) differ ==", text)
    k = re.search(r"paired keys: (\d+), VEP only: (\d+), fastVEP only: (\d+)", text)
    if not (m and k):
        sys.exit(f"could not read the output of compare_tab.py:\n{text[-800:]}")
    per_column = {c: int(n) for c, n in re.findall(r"^\[(\w+)\] (\d+) differ", text, re.M)}
    paired, cells, differ = (int(x) for x in m.groups())
    return dict(paired_rows=paired, cells_compared=cells, cells_differ=differ, columns_compared=cells // paired if paired else 0,
                only_vep=int(k.group(2)), only_fastvep=int(k.group(3)), per_column=per_column)


def multi_alt_keys(vcf: Path) -> set[str]:
    """The names VEP gives the records that have several ALTs: `chrom_pos_ref/alt1/alt2`."""
    keys = set()
    with vcf.open() as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.split("\t", 5)
            if "," in c[4]:
                keys.add(f"{c[0]}_{c[1]}_{c[3]}/{c[4].replace(',', '/')}")
    return keys


def classify(diff_dir: Path, vcf: Path) -> dict[str, dict[str, int]]:
    """{column: {cause: cells}} for the cells that differ: the causes are those of docs/VEP_DIVERGENCE.md."""
    multi = multi_alt_keys(vcf)
    out: dict[str, dict[str, int]] = {}
    for path in sorted(diff_dir.glob("diff_*.tsv")):
        col = path.stem[len("diff_"):]
        tally: collections.Counter = collections.Counter()
        with path.open() as f:
            next(f)
            for line in f:
                variant, _allele, _feature, v, o = line.rstrip("\n").split("\t")
                if variant in multi and col in ("HGVSc", "HGVS_OFFSET", "HGVSg"):
                    tally["divergence_11_several_alts"] += 1
                elif col == "HGVSc" and ":c.*" in v and ":c.*" not in o:
                    tally["divergence_4_terminator_intron"] += 1
                elif col == "HGVSp" and v == "-" and o != "-":
                    tally["divergence_17_three_prime_rule"] += 1
                else:
                    tally["unclassified"] += 1
        out[col] = dict(tally)
    return out


def sorted_gff3(src: Path, dst: Path) -> Path:
    """A bgzip and tabix copy of a GFF3, sorted, which is what VEP reads in --gff mode."""
    if not (dst.exists() and dst.with_name(dst.name + ".tbi").exists()):
        text = subprocess.run(["grep", "-v", "^#", str(src)], capture_output=True, text=True).stdout.split("\n")
        rows = sorted((l for l in text if l), key=lambda l: (l.split("\t")[0], int(l.split("\t")[3]), int(l.split("\t")[4])))
        header = [l for l in src.read_text().split("\n") if l.startswith("#")]
        body = "\n".join(header + rows) + "\n"
        with dst.open("wb") as f:
            p = subprocess.run(["bgzip", "-c"], input=body.encode(), stdout=f)
            if p.returncode != 0:
                sys.exit("bgzip failed")
        run(["tabix", "-f", "-p", "gff", str(dst)])
    return dst


def refseq_chr2122(gff_gz: Path, dst: Path) -> int:
    """The records of RefSeq's GFF3 on chromosomes 21 and 22, with the contigs renamed; the number of records."""
    n = 0
    with gzip.open(gff_gz, "rt") as fin, dst.open("w") as fout:
        fout.write("##gff-version 3\n")
        for line in fin:
            if line.startswith("#"):
                continue
            c = line.split("\t", 1)
            if c[0] in REFSEQ_CONTIGS:
                fout.write(REFSEQ_CONTIGS[c[0]] + "\t" + c[1])
                n += 1
    return n


def keep_chromosomes(src: Path, dst: Path, chroms: set[str]) -> int:
    n = 0
    with src.open() as fin, dst.open("w") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
            elif line.split("\t", 1)[0] in chroms:
                fout.write(line)
                n += 1
    return n


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bin", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--only", default="", help=f"comma-separated: {', '.join(ALL)}")
    ap.add_argument("--work", type=Path, default=ROOT / "benchmarks" / "output" / "vep116")
    ap.add_argument("--refseq-gff3", type=Path, help="NCBI's GCF_000001405.40_GRCh38.p14_genomic.gff.gz, for the refseq input")
    args = ap.parse_args()
    parts = [p for p in (args.only.split(",") if args.only else ALL) if p]
    if not args.only and args.refseq_gff3 is None:
        parts.remove("refseq")  # the one input whose source file no script of the repository fetches
    for p in parts:
        if p not in ALL:
            sys.exit(f"unknown input {p!r}")
    if "refseq" in parts and (args.refseq_gff3 is None or not args.refseq_gff3.exists()):
        sys.exit("the refseq input needs --refseq-gff3 (NCBI's GCF_000001405.40_GRCh38.p14_genomic.gff.gz)")
    info = subprocess.run(["docker", "info", "--format", "{{.ServerVersion}} cpus={{.NCPU}} mem={{.MemTotal}}"], capture_output=True, text=True).stdout.strip()
    if not info or "cpus=0" in info:
        sys.exit("Docker's engine is not answering")
    args.out.mkdir(parents=True, exist_ok=True)
    work = args.work.resolve()
    work.mkdir(parents=True, exist_ok=True)
    binary = args.bin.resolve()
    started = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    human_gff3, human_fa = HUMAN / "Homo_sapiens.GRCh38.115.gff3", HUMAN / "Homo_sapiens.GRCh38.dna.primary_assembly.fa"
    human_sorted = HUMAN / "Homo_sapiens.GRCh38.115.sorted.gff3.gz"
    if any(p in parts for p in ("example", "inframe", "hg002", "clinvar")):
        if not human_sorted.exists():
            sys.exit(f"missing {human_sorted}: validation/run_divergence.sh makes it")
    rows, hg002_cols = [], None
    timings = []

    def record(name: str, label: str, vcf: Path, r: dict, t_vep: float, t_our: float) -> None:
        n_rec = sum(1 for l in vcf.open() if not l.startswith("#"))
        print(f"{label:22s} {n_rec:>8,} records {r['paired_rows']:>10,} rows {r['cells_compared']:>12,} cells {r['cells_differ']:>7,} differ   "
              f"(VEP {t_vep:.0f} s, fastVEP {t_our:.1f} s)", flush=True)
        timings.append((label, t_vep, t_our))
        rows.append(dict(input=name, records=n_rec, **{k: v for k, v in r.items() if k != "per_column"}))

    def human(name: str, vcf: Path, with_hgvs: bool, label: str, classify_causes: bool = False) -> None:
        """`vcf` is in the work directory, which the container mounts as /data."""
        nonlocal hg002_cols
        vep_out, our_out = f"{label}.vep.tab", f"{label}.fastvep.tab"
        t_vep = vep(work, {"human": HUMAN}, vcf.name, "/human/" + human_sorted.name, "/human/" + human_fa.name, vep_out, with_hgvs)
        t_our = fastvep(binary, work, vcf, human_gff3, human_fa, our_out, with_hgvs)
        diff_dir = args.out / "vep116_diff" / label
        r = compare(work, vep_out, our_out, HUMAN_IGNORE, diff_dir)
        record(name, label, vcf, r, t_vep, t_our)
        if classify_causes:
            hg002_cols = (r, classify(diff_dir, vcf))
        for f in (work / vep_out, work / our_out):
            f.unlink()

    def small(name: str, d: Path) -> None:
        """A generated world: `d` holds in.vcf, genes.gff3 and genome.fa, and the container mounts it as /data."""
        sorted_gff3(d / "genes.gff3", d / "genes.sorted.gff3.gz")
        t_vep = vep(d, {}, "in.vcf", "/data/genes.sorted.gff3.gz", "/data/genome.fa", "vep.tab", True)
        t_our = fastvep(binary, d, d / "in.vcf", d / "genes.gff3", d / "genome.fa", "fastvep.tab", True)
        r = compare(d, "vep.tab", "fastvep.tab", WORLD_IGNORE, None)
        record(name, d.name, d / "in.vcf", r, t_vep, t_our)

    # ---- the human inputs -----------------------------------------------------------------------------------------------------
    if "example" in parts:
        shutil_copy(HERE / "human" / "vep_example_GRCh38.vcf", work / "example.vcf")
        human("example", work / "example.vcf", True, "example")
    if "inframe" in parts:
        shutil_copy(HERE / "human" / "clinvar_inframe_deletions.vcf", work / "inframe.vcf")
        human("inframe_deletions", work / "inframe.vcf", True, "inframe")
    if "hg002" in parts:
        run([sys.executable, str(HERE / "sample_variants.py"), str(HUMAN / "human_giab_hg002_full.vcf"), str(work / "hg002.vcf"), "--every", "200"])
        human("hg002_1in200", work / "hg002.vcf", True, "hg002", classify_causes=True)
    if "clinvar" in parts:
        run([sys.executable, str(HERE / "sample_variants.py"), str(CLINVAR), str(work / "clinvar.vcf"), "--every", "5"])
        human("clinvar_1in5", work / "clinvar.vcf", False, "clinvar_nohgvs")
        human("clinvar_1in5_hgvs", work / "clinvar.vcf", True, "clinvar_hgvs")

    # ---- the generated worlds -----------------------------------------------------------------------------------------------------
    world_dir = work / "world_gen"
    if "world" in parts:
        world_dir.mkdir(parents=True, exist_ok=True)
        run([sys.executable, str(FIXTURES / "vep_tab_world" / "gen_world.py"), str(world_dir), "3000", "139", "--all"])
        small("world", world_dir)
    fx_world = FIXTURES / "vep_tab_world"
    for name, gen, extra in (("edges", FIXTURES / "vep_tab_edges" / "gen_edges.py", ["1"]), ("lone", FIXTURES / "vep_tab_lone_indels" / "gen_lone.py", ["1500", "139"])):
        if name in parts:
            d = work / f"{name}_gen"
            d.mkdir(parents=True, exist_ok=True)
            run([sys.executable, str(gen), str(d / "in.vcf"), *extra])
            for f in ("genes.gff3", "genome.fa"):
                shutil_copy(fx_world / f, d / f)
            small("lone_indels" if name == "lone" else name, d)
    if "refseq" in parts:
        d = work / "refseq_gen"
        d.mkdir(parents=True, exist_ok=True)
        n_gff = refseq_chr2122(args.refseq_gff3, d / "refseq.gff3")
        if not (work / "hg002.vcf").exists():
            run([sys.executable, str(HERE / "sample_variants.py"), str(HUMAN / "human_giab_hg002_full.vcf"), str(work / "hg002.vcf"), "--every", "200"])
        keep_chromosomes(work / "hg002.vcf", d / "in.vcf", set(REFSEQ_CONTIGS.values()))
        print(f"refseq: {n_gff:,} records of the GFF3 on chromosomes 21 and 22", flush=True)
        sorted_gff3(d / "refseq.gff3", d / "refseq.sorted.gff3.gz")
        t_vep = vep(d, {"human": HUMAN}, "in.vcf", "/data/refseq.sorted.gff3.gz", "/human/" + human_fa.name, "vep.tab", True)
        t_our = fastvep(binary, d, d / "in.vcf", d / "refseq.gff3", human_fa, "fastvep.tab", True)
        r = compare(d, "vep.tab", "fastvep.tab", REFSEQ_IGNORE, args.out / "vep116_diff" / "refseq")
        record("refseq_gff3", "refseq", d / "in.vcf", r, t_vep, t_our)
        with (args.out / "vep116_refseq_columns.csv").open("w", newline="") as f:
            w = csv.writer(f, lineterminator="\n")
            w.writerow(["column", "paired_rows", "cells_differ"])
            for c, n in sorted(r["per_column"].items(), key=lambda kv: -kv[1]):
                w.writerow([c, r["paired_rows"], n])
    # ---- what is written ------------------------------------------------------------------------------------------------------------
    # The rows of an earlier run of other inputs stay.
    csv_path = args.out / "vep116_agreement.csv"
    header = ["input", "records", "paired_rows", "cells_compared", "cells_differ", "columns_compared", "only_vep", "only_fastvep"]
    kept = {}
    if csv_path.exists():
        with csv_path.open(newline="") as f:
            kept = {r["input"]: [r[c] for c in header] for r in csv.DictReader(f)}
    for r in rows:
        kept[r["input"]] = [r[c] for c in header]
    with csv_path.open("w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(header)
        for name in sorted(kept, key=ORDER.index):
            w.writerow(kept[name])
    if hg002_cols:
        r, causes = hg002_cols
        shared = ["Consequence", "IMPACT", "Amino_acids", "Codons", "cDNA_position", "CDS_position", "Protein_position", "Location", "REF_ALLELE", "VARIANT_CLASS",
                  "DISTANCE", "STRAND", "EXON", "INTRON", "SYMBOL", "BIOTYPE", "HGVSc", "HGVS_OFFSET", "HGVSg", "HGVSp"]
        with (args.out / "vep116_hg002_columns.csv").open("w", newline="") as f:
            w = csv.writer(f, lineterminator="\n")
            w.writerow(["column", "paired_rows", "cells_differ", "divergence_11_several_alts", "divergence_4_terminator_intron", "divergence_17_three_prime_rule", "unclassified"])
            for c in shared:
                ca = causes.get(c, {})
                w.writerow([c, r["paired_rows"], r["per_column"].get(c, 0), ca.get("divergence_11_several_alts", 0), ca.get("divergence_4_terminator_intron", 0),
                            ca.get("divergence_17_three_prime_rule", 0), ca.get("unclassified", 0)])
            other = {c: n for c, n in r["per_column"].items() if c not in shared}
            for c, n in other.items():
                w.writerow([c, r["paired_rows"], n, 0, 0, 0, n])
    version = run(["docker", "run", "--rm", IMAGE, "vep", "--help"]).stdout
    vline = re.search(r"ensembl-vep\s*:\s*(\S+)", version)
    with (args.out / "vep116_run.tsv").open("a") as f:  # one block for each run of the script
        f.write(f"started_utc\t{started}\nimage\t{IMAGE}\nvep\t{vline.group(1) if vline else '?'}\ndocker\t{info}\n")
        f.write(f"fastvep\t{binary.name}\t{run([str(binary), '--version']).stdout.strip()}\tsha256 {hashlib.sha256(binary.read_bytes()).hexdigest()}\n")
        f.write("inputs\t" + ",".join(parts) + "\n")
        for label, tv, to in timings:
            f.write(f"seconds_{label}\tVEP {tv:.0f}\tfastVEP {to:.1f}\n")
        f.write("\n")
    print("wrote vep116_agreement.csv" + (", vep116_hg002_columns.csv" if hg002_cols else "") + ", vep116_run.tsv")


def shutil_copy(src: Path, dst: Path) -> None:
    dst.write_bytes(src.read_bytes())


if __name__ == "__main__":
    main()
