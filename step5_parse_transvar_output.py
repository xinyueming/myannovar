#!/usr/bin/env python3
"""Step 5: Parse TransVar output into transvar.multianno format.

Parses transvar.output, extracts gene/transcript/region/cDNA/protein info,
groups by input variant, and writes a two-column transvar.multianno file.
"""

import re
import sys
import argparse
from collections import OrderedDict
from pathlib import Path
from typing import List, Optional


_REGION_EXON_RE = re.compile(r"exon_(\d+)")
_REGION_INTRON_RE = re.compile(r"intron_between_exon_(\d+)_and_(\d+)")
_COORD_SEP = "/"


def parse_region(region: str) -> str:
    """Convert region field to exon/intron label."""
    m = _REGION_INTRON_RE.search(region)
    if m:
        return f"intron{m.group(1)}"
    m = _REGION_EXON_RE.search(region)
    if m:
        return f"exon{m.group(1)}"
    return region.strip("[] ").split("_")[-1] if region else "-"


def parse_coordinates(coords: str):
    """Split gDNA/cDNA/protein, return (cDNA, protein) or ('-', '-')."""
    parts = coords.split(_COORD_SEP)
    if len(parts) >= 3:
        return parts[1], parts[2]
    if len(parts) == 2:
        return parts[0], "-"
    return "-", "-"


def parse_line(line: str):
    """Parse a single transvar.output line. Returns dict or None."""
    fields = line.rstrip().split("\t")
    if len(fields) < 7:
        return None
    return {
        "input": fields[0],
        "transcript": fields[1],
        "gene": fields[2],
        "strand": fields[3],
        "coordinates": fields[4],
        "region": fields[5],
        "info": fields[6],
    }


def process(input_path: str, output_path: Optional[str] = None) -> None:
    if not Path(input_path).is_file():
        sys.exit(f"Error: {input_path} not found")
    if output_path is None:
        output_path = Path(input_path).parent / "transvar.multianno"
    else:
        output_path = Path(output_path)

    out_dir = output_path.parent
    if not out_dir.is_dir():
        sys.exit(f"Error: output directory {out_dir} does not exist")

    # Group by input, preserving order
    groups: OrderedDict[str, List[str]] = OrderedDict()

    for line in Path(input_path).read_text().splitlines():
        stripped = line.rstrip()
        if not stripped or stripped.startswith("#"):
            continue
        parsed = parse_line(stripped)
        if parsed is None:
            continue

        region = parse_region(parsed["region"])
        cDNA, protein = parse_coordinates(parsed["coordinates"])
        entry = f"{parsed['gene']}:{parsed['transcript']}:{region}:{cDNA}:{protein}"

        key = parsed["input"]
        if key not in groups:
            groups[key] = []
        groups[key].append(entry)

    # Write output
    out_lines = ["transvar.input\tAAChange.transvar"]
    for key, entries in groups.items():
        out_lines.append(f"{key}\t{','.join(entries)}")

    Path(output_path).write_text("\n".join(out_lines) + "\n")
    print(f"Done. {len(groups)} variants written to {output_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Parse TransVar output into transvar.multianno format"
    )
    parser.add_argument("input", help="Input transvar.output file")
    parser.add_argument("-o", "--output", help="Output file (default: transvar.multianno next to input)")
    args = parser.parse_args()
    process(args.input, args.output)
