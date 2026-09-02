#!/usr/bin/env python3
"""Step 2: Merge transvar.input from .avinput into _multianno.txt.

Matches rows by (Chr, Start, End, Ref, Alt) and appends the transvar.input
value as a new column at the end of _multianno.txt.
"""

import sys
import argparse
from pathlib import Path
from typing import Optional


def load_avinput(avinput_path: str) -> dict[tuple[str, str, str, str, str], str]:
    """Load .avinput and return a dict keyed by (chr,start,end,ref,alt) -> transvar.input."""
    lookup: dict[tuple[str, str, str, str, str], str] = {}
    for line in Path(avinput_path).read_text().splitlines():
        stripped = line.rstrip()
        if not stripped or stripped.startswith("#"):
            continue
        fields = stripped.split("\t")
        if len(fields) < 6:
            continue
        key = (fields[0], fields[1], fields[2], fields[3], fields[4])
        lookup[key] = fields[5]
    return lookup


def process(
    avinput_path: str,
    multianno_path: str,
    output_path: Optional[str] = None,
) -> None:
    if not Path(avinput_path).is_file():
        sys.exit(f"Error: {avinput_path} not found")
    if not Path(multianno_path).is_file():
        sys.exit(f"Error: {multianno_path} not found")
    if output_path is None:
        output_path = multianno_path

    out_dir = Path(output_path).parent
    if not out_dir.is_dir():
        sys.exit(f"Error: output directory {out_dir} does not exist")

    lookup = load_avinput(avinput_path)
    lines = Path(multianno_path).read_text().splitlines()
    out_lines: list[str] = []

    header_added = False
    matched = 0
    unmatched = 0

    for line in lines:
        stripped = line.rstrip()
        if not stripped:
            out_lines.append(stripped)
            continue

        # Handle comment/header lines — pass through
        if stripped.startswith("#"):
            # Find the last header line (column names)
            if header_added:
                out_lines.append(stripped)
                continue

            # This is likely the column header line
            fields = stripped.split("\t")
            # Detect if this looks like a multianno header
            if "Chr" in fields or "chr" in fields or fields[0].lower().startswith("chr"):
                out_lines.append(f"{stripped}\ttransvar.input")
                header_added = True
            else:
                out_lines.append(stripped)
            continue

        fields = stripped.split("\t")
        # Multianno columns: Chr(0), Start(1), End(2), Ref(3), Alt(4), ...
        if len(fields) >= 5:
            key = (fields[0], fields[1], fields[2], fields[3], fields[4])
            tv = lookup.get(key, "-")
            if tv != "-":
                matched += 1
            else:
                unmatched += 1
            out_lines.append(f"{stripped}\t{tv}")
        else:
            out_lines.append(f"{stripped}\t-")
            unmatched += 1

    Path(output_path).write_text("\n".join(out_lines) + "\n")
    print(
        f"Done. {matched} matched, {unmatched} unmatched. "
        f"Output written to {output_path}"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Merge transvar.input from .avinput into _multianno.txt"
    )
    parser.add_argument("avinput", help="Input .avinput file (with transvar.input column)")
    parser.add_argument("multianno", help="Input _multianno.txt file")
    parser.add_argument("-o", "--output", help="Output file (default: overwrite multianno)")
    args = parser.parse_args()
    process(args.avinput, args.multianno, args.output)
