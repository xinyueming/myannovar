#!/usr/bin/env python3
"""Step 1: Add transvar.input column to .avinput file.

Reads .avinput (chr, start, end, ref, alt) and appends a transvar.input column
based on variant type classification.
"""

import re
import argparse
from pathlib import Path


def build_transvar_input(chr_: str, start: str, end: str, ref: str, alt: str) -> str:
    """Generate transvar.input string from variant fields."""
    if ref == "0" or alt == "0":
        return "-"

    if ref == "-":
        return f"{chr_}:g.{start}_{end}ins{alt}"

    if alt == "-":
        return f"{chr_}:g.{start}_{end}del"

    if re.match(r"[ACTGactg]+", ref) and re.match(r"[ACTGactg]+", alt):
        return f"{chr_}:g.{start}_{end}delins{alt}"

    return "-"


def process(input_path: str, output_path: str | None = None) -> None:
    if output_path is None:
        output_path = input_path

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

    # Add header column name if first non-empty/non-comment line has 5+ fields
    # and the original first line looks like a header
    first_data = out_lines[0] if out_lines else ""
    if first_data and not first_data.startswith("#"):
        header_fields = first_data.split("\t")
        # If the first line itself has data fields, it might not have a separate header
        # Only prepend header if the file has a header-like first line
        pass  # No change needed for headerless files

    Path(output_path).write_text("\n".join(out_lines) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Add transvar.input column to .avinput")
    parser.add_argument("input", help="Input .avinput file")
    parser.add_argument("-o", "--output", help="Output file (default: overwrite input)")
    args = parser.parse_args()
    process(args.input, args.output)
