#!/usr/bin/env python3
"""Step 5: Parse TransVar output into transvar.multianno format.

Parses transvar.output, extracts gene/transcript/region/cDNA/protein info,
groups by input variant, and writes transvar.multianno with optional
CHROM/POS/REF/ALT from --gseq output.
"""

import argparse
import re
import sys
from collections import OrderedDict
from pathlib import Path
from typing import Dict, Optional

_REGION_EXON_RE = re.compile(r"exon_(\d+)")
_REGION_INTRON_RE = re.compile(r"intron_between_exon_(\d+)_and_(\d+)")
_TRANSCRIPT_CLEAN = re.compile(r"\s*\([^)]*\)$")
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
    """Parse a single transvar.output line. Returns dict or None.

    Without --gseq: 7 columns (input..info)
    With --gseq:    11 columns (adds CHROM, POS, REF, ALT)
    """
    fields = line.rstrip().split("\t")
    if len(fields) < 7:
        return None
    result: Dict[str, str] = {
        "input": fields[0],
        "transcript": fields[1],
        "gene": fields[2],
        "strand": fields[3],
        "coordinates": fields[4],
        "region": fields[5],
        "info": fields[6],
        "chrom": fields[7] if len(fields) > 7 else "",
        "pos": fields[8] if len(fields) > 8 else "",
        "ref": fields[9] if len(fields) > 9 else "",
        "alt": fields[10] if len(fields) > 10 else "",
    }
    return result


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
    # groups[key] = {"aachange": [...], "chrom": str, "pos": str, ...}
    groups: OrderedDict[str, dict] = OrderedDict()

    for line in Path(input_path).read_text().splitlines():
        stripped = line.rstrip()
        if not stripped or stripped.startswith("#"):
            continue
        # Skip header line
        if stripped.startswith("input\t"):
            continue
        parsed = parse_line(stripped)
        if parsed is None:
            continue

        region = parse_region(parsed["region"])
        cDNA, protein = parse_coordinates(parsed["coordinates"])
        transcript = _TRANSCRIPT_CLEAN.sub("", parsed["transcript"])
        entry = f"{parsed['gene']}:{transcript}:{region}:{cDNA}:{protein}"

        key = parsed["input"]
        if key not in groups:
            groups[key] = {
                "aachange": [],
                "chrom": parsed.get("chrom", ""),
                "pos": parsed.get("pos", ""),
                "ref": parsed.get("ref", ""),
                "alt": parsed.get("alt", ""),
            }
        groups[key]["aachange"].append(entry)

    # Detect whether we have gseq data
    has_gseq = any(
        g["chrom"] for g in groups.values()
    )

    # Write output
    if has_gseq:
        out_lines = [
            "transvar.input\tAAChange.transvar\tCHROM\tPOS\tREF\tALT"
        ]
        for key, g in groups.items():
            out_lines.append(
                f"{key}\t{','.join(g['aachange'])}\t{g['chrom']}\t{g['pos']}\t{g['ref']}\t{g['alt']}"
            )
    else:
        out_lines = ["transvar.input\tAAChange.transvar"]
        for key, g in groups.items():
            out_lines.append(f"{key}\t{','.join(g['aachange'])}")

    Path(output_path).write_text("\n".join(out_lines) + "\n")
    print(f"Done. {len(groups)} variants written to {output_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Parse TransVar output into transvar.multianno format"
    )
    parser.add_argument("input", help="Input transvar.output file")
    parser.add_argument(
        "-o", "--output",
        help="Output file (default: transvar.multianno next to input)",
    )
    args = parser.parse_args()
    process(args.input, args.output)
