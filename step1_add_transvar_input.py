#!/usr/bin/env python3
"""Step 1: Add transvar.input column to .avinput file.

Reads .avinput (chr, start, end, ref, alt) and appends a transvar.input column
based on variant type classification.
"""

import re
import sys
import argparse
from pathlib import Path
from typing import Optional


_SEQ_RE = re.compile(r"^[ACTGactg]+$")


def build_transvar_input(chr_: str, start: str, end: str, ref: str, alt: str) -> str:
    """Generate transvar.input string from variant fields."""
    if ref == "0" or alt == "0":
        return "-"

    if ref == "-":
        return f"{chr_}:g.{start}_{end}ins{alt}"

    if alt == "-":
        return f"{chr_}:g.{start}_{end}del"

    if _SEQ_RE.match(ref) and _SEQ_RE.match(alt):
        return f"{chr_}:g.{start}_{end}delins{alt}"

    return "-"


def process(input_path: str, output_path: Optional[str] = None) -> None:
    if not Path(input_path).is_file():
        sys.exit(f"Error: {input_path} not found")

    if output_path is None:
        output_path = input_path

    out_dir = Path(output_path).parent
    if not out_dir.is_dir():
        sys.exit(f"Error: output directory {out_dir} does not exist")

    lines = Path(input_path).read_text().splitlines()
    out_lines: list[str] = []

    for line in lines:
        stripped = line.rstrip()
        if not stripped or stripped.startswith("#"):
            out_lines.append(stripped)
            continue

        fields = stripped.split("\t")
        if len(fields) < 5:
            out_lines.append(f"{stripped}\t-")
            continue

        chr_, start, end, ref, alt = fields[0], fields[1], fields[2], fields[3], fields[4]
        tv = build_transvar_input(chr_, start, end, ref, alt)
        out_lines.append(f"{stripped}\t{tv}")

    Path(output_path).write_text("\n".join(out_lines) + "\n")

    variant_count = len([l for l in out_lines if l and not l.startswith("#")])
    print(f"Done. {variant_count} variants processed, output written to {output_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Add transvar.input column to .avinput")
    parser.add_argument("input", help="Input .avinput file")
    parser.add_argument("-o", "--output", help="Output file (default: overwrite input)")
    args = parser.parse_args()
    process(args.input, args.output)
