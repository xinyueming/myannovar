#!/usr/bin/env python3
"""Step 7: Convert _multianno.txt to VCF format.

Converts ANNOVAR multianno output to standard VCF 4.2 format,
preserving AAChange.transvar annotation in the INFO field.
"""

import sys
import argparse
from pathlib import Path
from typing import Optional, List


def build_info_fields(fields: List[str], header_fields: List[str]) -> str:
    """Build VCF INFO string from multianno columns."""
    parts = []
    for i, name in enumerate(header_fields):
        if i == 0:
            continue  # skip Chr
        if name in ("Chr", "Start", "End", "Ref", "Alt", "transvar.input"):
            continue
        value = fields[i] if i < len(fields) else "."
        # Escape commas and semicolons in values
        value = value.replace(";", "%3B").replace("=", "%3D")
        key = name.replace(".", "_")
        parts.append(f"{key}={value}")
    return ";".join(parts) if parts else "."


def process(input_path: str, output_path: Optional[str] = None) -> None:
    if not Path(input_path).is_file():
        sys.exit(f"Error: {input_path} not found")
    if output_path is None:
        output_path = Path(input_path).with_suffix(".vcf")
    else:
        output_path = Path(output_path)

    out_dir = Path(output_path).parent
    if not out_dir.is_dir():
        sys.exit(f"Error: output directory {out_dir} does not exist")

    lines = Path(input_path).read_text().splitlines()
    out_lines: List[str] = []

    # VCF header
    out_lines.append("##fileformat=VCFv4.2")
    out_lines.append('##INFO=<ID=ANN,Number=.,Type=String,Description="ANNOVAR annotations">')
    out_lines.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")

    header_fields: List[str] = []
    variant_count = 0

    for line in lines:
        stripped = line.rstrip()
        if not stripped:
            continue

        if stripped.startswith("#"):
            # Parse header for column names
            fields = stripped.lstrip("# ").split("\t")
            if fields[0].lower() in ("chr", "chrom", "chromosome"):
                header_fields = fields
            continue

        fields = stripped.split("\t")
        if len(fields) < 5:
            continue

        chrom = fields[0]
        pos = fields[1]
        ref = fields[3]
        alt = fields[4]

        # Normalize indels
        if ref == "-":
            # Insertion: use empty ref, alt as-is
            ref = ""
        elif alt == "-":
            # Deletion: use empty alt
            alt = ""

        info = build_info_fields(fields, header_fields)

        out_lines.append(f"{chrom}\t{pos}\t.\t{ref if ref else '.'}\t{alt if alt else '.'}\t.\t.\t{info}")
        variant_count += 1

    Path(output_path).write_text("\n".join(out_lines) + "\n")
    print(f"Done. {variant_count} variants written to {output_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Convert _multianno.txt to VCF format"
    )
    parser.add_argument("input", help="Input _multianno.txt file")
    parser.add_argument("-o", "--output", help="Output VCF file (default: same name with .vcf extension)")
    args = parser.parse_args()
    process(args.input, args.output)
