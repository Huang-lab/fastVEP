#!/usr/bin/env python3
"""Which JSON records two runs of `fastvep annotate --output-format json` over one input write differently.

    benchmarks/json_line_diff.py A.json B.json

fastVEP pretty-prints its JSON, a key to a line and the keys of an object in order, so two outputs of one input line up line for
line for as long as their records have the same shape. This reads both in step and prints a line for every line that differs:
the record (contig, position, alleles), the key the line is on, and the two values. At the end it says how many lines were
compared, how many differ and how many records they are in.

It stops at the first line that one file has and the other has not, since what follows no longer lines up, and says so.
Exit status: 0 when the files are the same, 1 when they differ, 2 when they stop lining up.
"""
import re
import sys

KEY = re.compile(r'^\s*"([^"]+)"\s*:\s*(.*?),?\s*$')
RECORD_FIELD = re.compile(r'^  "(allele_string|seq_region_name|start)"\s*:\s*(.*?),?\s*$')


def main(a_path: str, b_path: str) -> int:
    compared = differing = 0
    records = set()
    record_no = 0
    current: dict[str, str] = {}
    # the differing lines of this record: its position (`start`) comes after its `most_severe_consequence`, so they are printed
    # when the record is over
    waiting: list[tuple[str, str, str, dict[str, str]]] = []
    status = 0

    def flush() -> None:
        for key, va, vb, record in waiting:
            print(f"{record.get('seq_region_name', '?')}\t{record.get('start', '?')}\t{record.get('allele_string', '?')}\t{key}\t{va}\t{vb}")
        waiting.clear()

    with open(a_path, "rb") as fa, open(b_path, "rb") as fb:
        while True:
            la, lb = fa.readline(), fb.readline()
            if not la and not lb:
                break
            compared += 1
            if not la or not lb:
                flush()
                print(f"line {compared}: one file ends and the other goes on; the rest does not line up", file=sys.stderr)
                status = 2
                break
            sa = la.decode("utf-8", "replace").rstrip("\n")
            if sa == "{":
                flush()
                record_no += 1
                current = {}
            m = RECORD_FIELD.match(sa)
            if m:
                current[m.group(1)] = m.group(2).strip('"')
            if la != lb:
                differing += 1
                records.add(record_no)
                sb = lb.decode("utf-8", "replace").rstrip("\n")
                ka, kb = KEY.match(sa), KEY.match(sb)
                key = (ka or kb).group(1) if (ka or kb) else "(no key)"
                va = ka.group(2) if ka else sa.strip()
                vb = kb.group(2) if kb else sb.strip()
                waiting.append((key, va, vb, current))
                status = status or 1
    flush()
    print(f"lines compared: {compared}; lines that differ: {differing}; records they are in: {len(records)}")
    return status


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    sys.exit(main(sys.argv[1], sys.argv[2]))
