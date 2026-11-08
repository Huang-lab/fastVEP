#!/usr/bin/env python3
"""Count what Table 1, Table 9 and Figure 1 of the manuscript quote about the code, from the source.

    manuscript/code_metrics.py --out benchmarks/results/v0.5.0        # writes code_metrics.csv, so_terms.csv, crate_dependencies.csv
    manuscript/code_metrics.py --root /path/to/another/checkout --print

The rules are the ones that reproduce the numbers of `code_metrics.csv` as they stood on 29 July 2026 (the lines of every crate, to the
line, at the commit that wrote them), so a change in a number is a change in the code and not in the way it was counted:

    loc_<crate>       the physical lines (newline characters, as `wc -l` counts) of every .rs file under crates/<crate>/src, blank
                      lines, comments and the `#[cfg(test)]` modules inside those files included; `tests/` and `examples/` are not
    tests_<crate>     the `#[test]` and `#[tokio::test]` attributes under crates/<crate>/src and crates/<crate>/tests (an `#[ignore]`d
                      test is counted, so this is what `cargo test` runs plus what it ignores; `tests_ignored` is how many)
    so_consequence_types           the variants of `enum Consequence` (crates/fastvep-core/src/consequence.rs)
    so_consequence_types_emitted   those that crates/fastvep-consequence/src names outside its test modules, which is to say the ones
                                   a predictor can produce: the others are defined and never emitted (`so_terms.csv` says which)
    csq_output_fields              the entries of DEFAULT_CSQ_FIELDS (crates/fastvep-io/src/output.rs), of which the five that
                                   output.rs lists as `NO_WRITER` are written empty
    web_gui_lines                  the lines of web/index.html
    crates                         the members of the workspace

`so_terms.csv` is the rank, the term, the impact and whether a predictor emits it, for each of the 49 variants: Supplementary Table S1 and the
colours of Figure 4 are drawn from it. `crate_dependencies.csv` is the edges between the workspace crates as `cargo metadata` gives them
(a normal dependency, not a dev-dependency): the arrows of Figure 1.

Standard library only; `cargo metadata` is run for the edges and skipped if cargo is not there.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
CRATES = ["fastvep-core", "fastvep-genome", "fastvep-cache", "fastvep-consequence", "fastvep-hgvs", "fastvep-io", "fastvep-filter",
          "fastvep-sa", "fastvep-annotate", "fastvep-classification", "fastvep-cli", "fastvep-web"]
TEST_ATTR = re.compile(r"^\s*#\[(?:tokio::)?test\b", re.M)
IGNORE_ATTR = re.compile(r"^\s*#\[ignore\b", re.M)


def newlines(path: Path) -> int:
    return path.read_bytes().count(b"\n")


def rs_files(directory: Path) -> list[Path]:
    return sorted(directory.rglob("*.rs")) if directory.is_dir() else []


def before_tests(path: Path) -> str:
    """The text of a file before its first `#[cfg(test)]`: what a predictor does, and not what its tests name."""
    text = path.read_text(errors="replace")
    i = text.find("#[cfg(test)]")
    return text if i < 0 else text[:i]


def enum_variants(consequence_rs: Path) -> list[str]:
    text = consequence_rs.read_text()
    body = text[text.index("pub enum Consequence"):]
    body = body[:body.index("\n}\n")]
    return re.findall(r"^\s{4}([A-Z][A-Za-z0-9]+),", body, re.M)


def match_arms(text: str, function: str) -> str:
    start = text.index(f"pub fn {function}(")
    return text[start:text.index("\n    }\n", start)]


def so_terms(root: Path) -> list[list[str]]:
    consequence_rs = root / "crates/fastvep-core/src/consequence.rs"
    text = consequence_rs.read_text()
    variants = enum_variants(consequence_rs)
    rank = {v: int(r) for v, r in re.findall(r"Self::([A-Za-z0-9]+) => (\d+),", match_arms(text, "rank"))}
    term = dict(re.findall(r'Self::([A-Za-z0-9]+) => "([a-zA-Z0-9_]+)"', match_arms(text, "so_term")))
    impact: dict[str, str] = {}
    for arm in re.split(r"\n\s*\n", match_arms(text, "impact")):
        names = re.findall(r"Self::([A-Za-z0-9]+)", arm)
        level = re.search(r"Impact::(\w+)", arm)
        if level:
            for n in names:
                impact[n] = level.group(1).upper()
    predictors = "\n".join(before_tests(p) for p in rs_files(root / "crates/fastvep-consequence/src"))
    rows = [["rank", "so_term", "impact", "emitted_by_a_predictor"]]
    for v in sorted(variants, key=lambda v: rank[v]):
        rows.append([str(rank[v]), term[v], impact[v], "yes" if re.search(rf"Consequence::{v}\b", predictors) else "no"])
    return rows


def csq_fields(root: Path) -> int:
    text = (root / "crates/fastvep-io/src/output.rs").read_text()
    start = text.index("pub const DEFAULT_CSQ_FIELDS")
    return len(re.findall(r'^\s+"[^"]+",$', text[start:text.index("];", start)], re.M))


def dependencies(root: Path) -> list[list[str]]:
    try:
        meta = json.loads(subprocess.run(["cargo", "metadata", "--offline", "--no-deps", "--format-version", "1"], cwd=root,
                                         capture_output=True, text=True, check=True).stdout)
    except (OSError, subprocess.CalledProcessError):
        return []
    names = {p["name"] for p in meta["packages"]}
    rows = [["crate", "depends_on"]]
    for p in sorted(meta["packages"], key=lambda p: CRATES.index(p["name"]) if p["name"] in CRATES else 99):
        for d in sorted({d["name"] for d in p["dependencies"] if d["name"] in names and d["kind"] in (None, "normal")}):
            rows.append([p["name"], d])
    return rows


def metrics(root: Path) -> list[list[str]]:
    rows = [["metric", "value"]]
    loc = {c: sum(newlines(f) for f in rs_files(root / "crates" / c / "src")) for c in CRATES}
    tests = {c: sum(len(TEST_ATTR.findall(f.read_text(errors="replace"))) for d in ("src", "tests") for f in rs_files(root / "crates" / c / d)) for c in CRATES}
    rows += [[f"loc_{c}", str(loc[c])] for c in CRATES]
    rows.append(["total_source_lines", str(sum(loc.values()))])
    rows += [[f"tests_{c}", str(tests[c])] for c in CRATES]
    rows.append(["tests_total", str(sum(tests.values()))])
    rows.append(["tests_ignored", str(sum(len(IGNORE_ATTR.findall(f.read_text(errors="replace"))) for c in CRATES for d in ("src", "tests") for f in rs_files(root / "crates" / c / d)))])
    terms = so_terms(root)
    rows.append(["so_consequence_types", str(len(terms) - 1)])
    rows.append(["so_consequence_types_emitted", str(sum(1 for t in terms[1:] if t[3] == "yes"))])
    rows.append(["web_gui_lines", str(newlines(root / "web/index.html"))])
    rows.append(["csq_output_fields", str(csq_fields(root))])
    rows.append(["crates", str(len(CRATES))])
    return rows


def write(path: Path, rows: list[list[str]]) -> None:
    with path.open("w", newline="") as f:
        csv.writer(f, lineterminator="\n").writerows(rows)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", type=Path, default=HERE.parent, help="the checkout to count (default: this one)")
    ap.add_argument("--out", type=Path, help="directory for code_metrics.csv, so_terms.csv and crate_dependencies.csv")
    ap.add_argument("--print", action="store_true", help="print code_metrics.csv")
    args = ap.parse_args()
    root = args.root.resolve()
    missing = [c for c in CRATES if not (root / "crates" / c).is_dir()]
    if missing:
        sys.exit(f"not a fastVEP checkout, or the crates have changed: missing {missing}")
    m, terms, deps = metrics(root), so_terms(root), dependencies(root)
    if args.print or not args.out:
        for k, v in m[1:]:
            print(f"{k:34s} {v}")
    if args.out:
        args.out.mkdir(parents=True, exist_ok=True)
        write(args.out / "code_metrics.csv", m)
        write(args.out / "so_terms.csv", terms)
        if deps:
            write(args.out / "crate_dependencies.csv", deps)


if __name__ == "__main__":
    main()
