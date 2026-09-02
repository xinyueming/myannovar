#!/usr/bin/env python3
"""Step 3: Extract transvar.input values from _multianno.txt.

Reads _multianno.txt, extracts the transvar.input column, filters out `-`,
deduplicates, and writes a clean transvar.input file for TransVar annotation.
"""

import sys
import argparse
from pathlib import Path
from typing import Optional


def process(input_path: str, output_path: Optional[str] = None) -> None:
    if not Path(input_path).is_file():
        sys.exit(f"Error: {input_path} not found")
    if output_path is None:
        output_path = Path(input_path).parent / "transvar.input"

    out_dir = Path(output_path).parent
    if not out_dir.is_dir():
        sys.exit(f"Error: output directory {out_dir} does not exist")

    lines = Path(input_path).read_text().splitlines()
    tv_col: Optional[int] = None
    values: set[str] = set()

    for line in lines:
        stripped = line.rstrip()
        if not stripped or stripped.startswith("#"):
            # Check if this is the header line with column names
            fields = stripped.split("\t")
            for i, f in enumerate(fields):
                if f.strip() == "transvar.input":
                    tv_col = i
                    break
            continue

        if tv_col is None:
            continue

        fields = stripped.split("\t")
        if tv_col < len(fields):
            val = fields[tv_col]
            if val != "-":
                values.add(val)

    out_lines = sorted(values)
    Path(output_path).write_text("\n".join(out_lines) + "\n")
    print(f"Done. {len(out_lines)} unique transvar.input entries written to {output_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract transvar.input from _multianno.txt"
    )
    parser.add_argument("input", help="Input _multianno.txt file")
    parser.add_argument("-o", "--output", help="Output file (default: transvar.input next to input)")
    args = parser.parse_args()
    process(args.input, args.output)
