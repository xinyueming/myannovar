#!/usr/bin/env python3
"""Step 6: Replace AAChange.refGeneWithVer in _multianno.txt with transvar annotations.

Matches rows by transvar.input column, replaces the AAChange.refGeneWithVer
column with the AAChange.transvar value from transvar.multianno, and
replaces Otherinfo4/5/7/8 (CHROM/POS/REF/ALT) with transvar's --gseq values.
"""

import argparse
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional


@dataclass(frozen=True)
class TransVarEntry:
    aachange: str
    chrom: str = ""
    pos: str = ""
    ref: str = ""
    alt: str = ""


def load_multianno(path: str) -> Dict[str, TransVarEntry]:
    """Load transvar.multianno and return dict keyed by transvar.input."""
    lookup: Dict[str, TransVarEntry] = {}
    lines = Path(path).read_text().splitlines()
    if not lines:
        return lookup

    # Detect header format
    header = lines[0].split("\t")
    has_gseq = len(header) >= 6 and header[2] == "CHROM"

    for line in lines[1:]:
        stripped = line.rstrip()
        if not stripped or stripped.startswith("#"):
            continue
        fields = stripped.split("\t")
        if len(fields) < 2:
            continue
        if has_gseq and len(fields) >= 6:
            lookup[fields[0]] = TransVarEntry(
                aachange=fields[1],
                chrom=fields[2],
                pos=fields[3],
                ref=fields[4],
                alt=fields[5],
            )
        else:
            lookup[fields[0]] = TransVarEntry(aachange=fields[1])

    return lookup


def find_column_index(header_fields: List[str], *names: str) -> Optional[int]:
    """Find column index by trying multiple candidate names."""
    for i, f in enumerate(header_fields):
        if f in names:
            return i
    return None


def _find_otherinfo_indices(header_fields: List[str]) -> Dict[int, str]:
    """Return {index: name} for Otherinfo4/5/7/8 columns."""
    result = {}
    for i, f in enumerate(header_fields):
        if f == "Otherinfo4":
            result[i] = "CHROM"
        elif f == "Otherinfo5":
            result[i] = "POS"
        elif f == "Otherinfo7":
            result[i] = "REF"
        elif f == "Otherinfo8":
            result[i] = "ALT"
    return result


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
    otherinfo_cols: Dict[int, str] = {}
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
                otherinfo_cols = _find_otherinfo_indices(fields)

                # Remove transvar.input column from header
                if tv_col is not None:
                    fields.pop(tv_col)
                    # After removal, adjust indices
                    if aachange_col is not None and aachange_col > tv_col:
                        aachange_col -= 1
                    otherinfo_cols = {
                        (i - 1 if i > tv_col else i): name
                        for i, name in otherinfo_cols.items()
                    }
                    out_lines.append("\t".join(fields))
                else:
                    out_lines.append(stripped)
                continue

        fields = stripped.split("\t")
        if tv_col is None or aachange_col is None:
            out_lines.append(stripped)
            unchanged += 1
            continue

        # Replace AAChange and Otherinfo if transvar annotation exists
        if tv_col < len(fields):
            tv_input = fields[tv_col]
            entry = tv_lookup.get(tv_input)
            if entry is not None:
                if aachange_col is not None and aachange_col < len(fields):
                    fields[aachange_col] = entry.aachange
                # Replace Otherinfo4/5/7/8 with transvar gseq values
                for idx, col_name in otherinfo_cols.items():
                    if idx < len(fields):
                        value = getattr(entry, col_name.lower(), "")
                        if value:
                            fields[idx] = value
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
