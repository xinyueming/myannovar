#!/usr/bin/env python3
"""Step 7: Convert _multianno.txt to VCF format.

Reconstructs a polished VCF by merging multianno annotations into the
original VCF's INFO field, preserving FORMAT and sample columns.
"""

import sys
import argparse
from pathlib import Path
from typing import List, Optional


def _vcf_escape_value(value: str) -> str:
    """Escape a value for VCF INFO field. Spaces become underscores."""
    if not value or value == ".":
        return "."
    # Spaces → underscores (ANNOVAR convention)
    value = value.replace(" ", "_")
    # VCF spec: escape ; and = in values
    value = value.replace(";", "%3B").replace("=", "%3D")
    return value


def process(
    multianno_path: str,
    original_vcf: str,
    output_path: Optional[str] = None,
    build: Optional[str] = None,
) -> None:
    if not Path(multianno_path).is_file():
        sys.exit(f"Error: {multianno_path} not found")
    if not Path(original_vcf).is_file():
        sys.exit(f"Error: {original_vcf} not found")
    if output_path is None:
        output_path = Path(multianno_path).with_suffix(".vcf")
    else:
        output_path = Path(output_path)

    out_dir = Path(output_path).parent
    if not out_dir.is_dir():
        sys.exit(f"Error: output directory {out_dir} does not exist")

    # 1. Read original VCF header
    vcf_lines = Path(original_vcf).read_text().splitlines()
    vcf_header = [line for line in vcf_lines if line.startswith("##")]
    vcf_header_line = None
    for line in vcf_lines:
        if line.startswith("#CHROM"):
            vcf_header_line = line
            break

    if vcf_header_line is None:
        sys.exit(f"Error: no #CHROM line found in {original_vcf}")

    # Parse VCF header columns
    vcf_cols = vcf_header_line.lstrip("#").split("\t")

    # 2. Read multianno
    ma_lines = Path(multianno_path).read_text().splitlines()
    ma_header = None
    ma_data_start = 0
    for i, line in enumerate(ma_lines):
        stripped = line.rstrip()
        if stripped:
            ma_header = stripped
            ma_data_start = i + 1
            break

    if ma_header is None:
        sys.exit(f"Error: no header found in {multianno_path}")

    ma_cols = ma_header.split("\t")
    n_ma_cols = len(ma_cols)

    # 3. Identify Otherinfo columns and annotation columns
    otherinfo_indices: List[int] = []
    for i, col in enumerate(ma_cols):
        if col.startswith("Otherinfo"):
            otherinfo_indices.append(i)

    if not otherinfo_indices:
        sys.exit(f"Error: no Otherinfo columns found in {multianno_path}")

    oi = lambda n: otherinfo_indices[n - 1] if n <= len(otherinfo_indices) else None

    info_col = oi(11)
    format_col = oi(12)
    sample_cols = [oi(n) for n in range(13, len(otherinfo_indices) + 1) if oi(n) is not None]

    # Annotation columns: from col 5 (Func...) to just before Otherinfo1
    annot_start = 5
    annot_end = otherinfo_indices[0]
    annot_col_names = ma_cols[annot_start:annot_end]

    # 4. Build output
    out_lines: List[str] = []
    out_lines.extend(vcf_header)
    out_lines.append(vcf_header_line)

    # Build ANNOVAR_DATE from build or use current date
    import datetime
    if build:
        annovar_date = build  # e.g. "hg19" → use as-is? No, expected is "2020-06-08"
        # The expected ANNOVAR_DATE is the date ANNOVAR database was built, not the build name
        # We'll just add a placeholder
        annovar_date = datetime.date.today().isoformat()
    else:
        annovar_date = datetime.date.today().isoformat()

    variant_count = 0
    for line in ma_lines[ma_data_start:]:
        stripped = line.rstrip()
        if not stripped:
            continue

        fields = stripped.split("\t")
        if len(fields) < n_ma_cols:
            fields.extend(["."] * (n_ma_cols - len(fields)))

        chrom = fields[0]
        pos = fields[1]
        ref = fields[3]
        alt = fields[4]

        ref_vcf = "" if ref == "-" else ref
        alt_vcf = "" if alt == "-" else alt

        qual = fields[oi(9)] if oi(9) is not None else "."
        filt = fields[oi(10)] if oi(10) is not None else "."
        id_val = fields[oi(6)] if oi(6) is not None else "."
        info_orig = fields[info_col] if info_col is not None else "."
        format_val = fields[format_col] if format_col is not None else "."

        # Build annotation key=value pairs (keep dots in column names)
        annot_parts = [f"ANNOVAR_DATE={annovar_date}"]
        for j, col_name in enumerate(annot_col_names):
            col_idx = annot_start + j
            value = fields[col_idx] if col_idx < len(fields) else "."
            escaped_value = _vcf_escape_value(value)
            annot_parts.append(f"{col_name}={escaped_value}")

        # Merge annotations into INFO
        if info_orig and info_orig != ".":
            info_merged = f"{info_orig};{';'.join(annot_parts)}"
        else:
            info_merged = ";".join(annot_parts)

        # Build output line
        out_parts = [chrom, pos, id_val, ref_vcf if ref_vcf else ".",
                     alt_vcf if alt_vcf else ".", qual, filt, info_merged]
        if format_val and format_val != ".":
            out_parts.append(format_val)
            for si in sample_cols:
                out_parts.append(fields[si] if si < len(fields) else ".")

        out_lines.append("\t".join(out_parts))
        variant_count += 1

    Path(output_path).write_text("\n".join(out_lines) + "\n")
    print(f"Done. {variant_count} variants written to {output_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Convert _multianno.txt to VCF format"
    )
    parser.add_argument("multianno", help="Input _multianno.txt file (with transvar.input removed)")
    parser.add_argument("original_vcf", help="Original input VCF file (for header and sample data)")
    parser.add_argument("-o", "--output", help="Output VCF file (default: same name with .vcf extension)")
    args = parser.parse_args()
    process(args.multianno, args.original_vcf, args.output)
