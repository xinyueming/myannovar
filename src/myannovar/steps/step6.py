#!/usr/bin/env python3
"""Step 6: Replace AAChange.refGeneWithVer in _multianno.txt with transvar annotations.

Matches rows by transvar.input column, and replaces the AAChange.refGeneWithVer
column with the AAChange.transvar value from transvar.multianno.
"""

import sys
import argparse
from pathlib import Path
from typing import Optional, Dict, List


def load_multianno(path: str) -> Dict[str, str]:
    """Load transvar.multianno and return dict keyed by transvar.input -> AAChange.transvar."""
    lookup: Dict[str, str] = {}
    for line in Path(path).read_text().splitlines():
        stripped = line.rstrip()
        if not stripped or stripped.startswith("#"):
            continue
        fields = stripped.split("\t")
        if len(fields) >= 2:
            lookup[fields[0]] = fields[1]
    return lookup


def find_column_index(header_fields: List[str], *names: str) -> Optional[int]:
    """Find column index by trying multiple candidate names."""
    for i, f in enumerate(header_fields):
        if f in names:
            return i
    return None


def process(
    multianno_path: str,
    transvar_multianno_path: str,
    output_path: Optional[str] = None,
) -> None:
    if not Path(multianno_path).is_file():
        sys.exit(f"Error: {multianno_path} not found")
    if not Path(transvar_multianno_path).is_file():
        sys.exit(f"Error: {transvar_multianno_path} not found")
    if output_path is None:
        output_path = multianno_path
    else:
        output_path = Path(output_path)

    out_dir = Path(output_path).parent
    if not out_dir.is_dir():
        sys.exit(f"Error: output directory {out_dir} does not exist")

    tv_lookup = load_multianno(transvar_multianno_path)

    lines = Path(multianno_path).read_text().splitlines()
    out_lines: List[str] = []

    tv_col: Optional[int] = None
    aachange_col: Optional[int] = None
    replaced = 0
    unchanged = 0

    for line in lines:
        stripped = line.rstrip()
        if not stripped:
            out_lines.append(stripped)
            continue

        # Detect header line (may not have # prefix in ANNOVAR multianno)
        if tv_col is None or aachange_col is None:
            fields = stripped.split("\t")
            first = fields[0].lstrip("#").strip()
            if first == "Chr" or first.lower() == "chr" or stripped.startswith("#"):
                if tv_col is None:
                    tv_col = find_column_index(fields, "transvar.input")
                if aachange_col is None:
                    aachange_col = find_column_index(
                        fields, "AAChange.refGeneWithVer", "AAChange.refGene"
                    )
                # Remove transvar.input column from header
                if tv_col is not None:
                    fields.pop(tv_col)
                    # After removal, aachange_col may shift
                    if aachange_col > tv_col:
                        aachange_col -= 1
                    out_lines.append("\t".join(fields))
                else:
                    out_lines.append(stripped)
                continue

        fields = stripped.split("\t")
        if tv_col is None or aachange_col is None:
            out_lines.append(stripped)
            unchanged += 1
            continue

        # Replace AAChange if transvar annotation exists
        if tv_col < len(fields):
            tv_input = fields[tv_col]
            replacement = tv_lookup.get(tv_input)
            if replacement is not None and aachange_col < len(fields):
                fields[aachange_col] = replacement
                replaced += 1
            else:
                unchanged += 1
        else:
            unchanged += 1

        # Remove transvar.input column
        if tv_col is not None and tv_col < len(fields):
            fields.pop(tv_col)

        out_lines.append("\t".join(fields))

    Path(output_path).write_text("\n".join(out_lines) + "\n")
    print(
        f"Done. {replaced} replaced, {unchanged} unchanged. "
        f"Output written to {output_path}"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Replace AAChange in _multianno.txt with transvar annotations"
    )
    parser.add_argument("multianno", help="Input _multianno.txt (with transvar.input column)")
    parser.add_argument("transvar_multianno", help="Input transvar.multianno file")
    parser.add_argument("-o", "--output", help="Output file (default: overwrite multianno)")
    args = parser.parse_args()
    process(args.multianno, args.transvar_multianno, args.output)
