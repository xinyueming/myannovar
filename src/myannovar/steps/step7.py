"""Step 7: Convert _multianno.txt to VCF format.

Port of ANNOVAR's backConvertVCF (table_annovar.pl:113-202).
Reconstructs a polished VCF by merging multianno annotations into the
original VCF's INFO field, preserving FORMAT and sample columns.
"""

from __future__ import annotations

import argparse
import datetime
import gzip
import re
import sys
from pathlib import Path

# Frequency column prefixes that get Type=Float in ##INFO header
_FREQ_PREFIX = re.compile(
    r"^(1000g\d+|esp\d+|cg\d+|popfreq|nci\d+|ExAC|gnomAD)"
)


def _anno_escape(value: str) -> str:
    """Escape a value for VCF INFO field, matching ANNOVAR backConvertVCF.

    ANNOVAR does:
      $field =~ s/\\s/_/g;       # spaces -> underscores
      $field =~ s/;/\\\\x3b/g;   # semicolons -> literal \\x3b
      $field =~ s/=/\\\\x3d/g;   # equals -> literal \\x3d
    """
    if not value or value == ".":
        return "."
    value = value.replace(" ", "_")
    value = value.replace(";", "\\x3b")
    value = value.replace("=", "\\x3d")
    return value


def _is_otherinfo(col: str) -> bool:
    return col.startswith("Otherinfo")


def _skip_for_info_header(col: str) -> bool:
    """Columns that should NOT get a ##INFO header line."""
    return col in ("Chr", "Start", "End", "Ref", "Alt") or _is_otherinfo(col)


def _anno_date() -> str:
    return datetime.datetime.now(tz=datetime.timezone.utc).date().isoformat()


def _find_otherinfo_indices(cols: list[str]) -> list[int]:
    """Return sorted list of indices for columns starting with Otherinfo."""
    return sorted(i for i, c in enumerate(cols) if c.startswith("Otherinfo"))


def _build_anno_string(
    fields: list[str],
    col_names: list[str],
    annot_end_idx: int,
    date: str,
) -> str:
    """Build the ANNOVAR annotation string for one variant.

    Mirrors Perl: for my $i (5 .. @name-2)
      anno_string .= ";$name[$i]=$field[$i]"
    Then appends ;ALLELE_END
    """
    parts = [f"ANNOVAR_DATE={date}"]
    for i in range(5, annot_end_idx):
        name = col_names[i]
        val = fields[i] if i < len(fields) else "."
        parts.append(f"{name}={_anno_escape(val)}")
    return ";" + ";".join(parts) + ";ALLELE_END"


def _vcf_columns_from_otherinfo(
    fields: list[str],
    oi: list[int],
    annot_string: str,
) -> list[str]:
    """Extract VCF columns from Otherinfo, appending annot_string to INFO.

    Mapping (1-indexed in VCF, 0-indexed in oi list):
      Otherinfo4 -> CHROM
      Otherinfo5 -> POS
      Otherinfo6 -> ID
      Otherinfo7 -> REF
      Otherinfo8 -> ALT
      Otherinfo9 -> QUAL
      Otherinfo10 -> FILTER
      Otherinfo11 -> INFO  (original VCF INFO, annot_string appended)
      Otherinfo12 -> FORMAT
      Otherinfo13+ -> SAMPLE columns
    """
    if len(oi) < 11:
        chrom = fields[0]
        pos = fields[1]
        id_val = fields[2] if len(fields) > 2 else "."
        ref = fields[3] if len(fields) > 3 else "."
        alt = fields[4] if len(fields) > 4 else "."
        qual = "."
        filt = "."
        info_orig = "."
        format_val = "."
        samples: list[str] = []
    else:
        chrom = fields[oi[3]] if oi[3] < len(fields) else "."
        pos = fields[oi[4]] if oi[4] < len(fields) else "."
        id_val = fields[oi[5]] if oi[5] < len(fields) else "."
        ref = fields[oi[6]] if oi[6] < len(fields) else "."
        alt = fields[oi[7]] if oi[7] < len(fields) else "."
        qual = fields[oi[8]] if oi[8] < len(fields) else "."
        filt = fields[oi[9]] if oi[9] < len(fields) else "."
        info_orig = fields[oi[10]] if oi[10] < len(fields) else "."
        format_val = fields[oi[11]] if oi[11] < len(fields) else "."
        samples = [
            fields[oi[n]] if oi[n] < len(fields) else "."
            for n in range(12, len(oi))
        ]

    # Build INFO: original VCF INFO + annotation string
    if info_orig and info_orig != ".":
        info_merged = f"{info_orig}{annot_string}"
    else:
        info_merged = annot_string.lstrip(";")

    out = [chrom, pos, id_val if id_val else ".",
           ref if ref else ".", alt if alt else ".",
           qual if qual else ".", filt if filt else ".",
           info_merged]
    if format_val and format_val != ".":
        out.append(format_val)
        out.extend(samples)

    return out


def _sample_key(fields: list[str], oi: list[int]) -> tuple:
    """Key for multi-allelic merging: sample columns from Otherinfo13+."""
    if len(oi) < 13:
        return ()
    return tuple(
        fields[oi[n]] if oi[n] < len(fields) else "."
        for n in range(12, len(oi))
    )


def process(
    multianno_path: str,
    original_vcf: str,
    output_path: str | None = None,
    build: str | None = None,
) -> None:
    if not Path(multianno_path).is_file():
        sys.exit(f"Error: {multianno_path} not found")
    if not Path(original_vcf).is_file():
        sys.exit(f"Error: {original_vcf} not found")
    if output_path is None:
        output_path = str(Path(multianno_path).with_suffix(".vcf"))
    else:
        output_path = str(output_path)

    out_dir = Path(output_path).parent
    if not out_dir.is_dir():
        sys.exit(f"Error: output directory {out_dir} does not exist")

    # 1. Read original VCF headers
    if str(original_vcf).endswith(".gz"):
        with gzip.open(original_vcf, "rt", encoding="utf-8") as f:
            vcf_lines = f.read().splitlines()
    else:
        vcf_lines = Path(original_vcf).read_text().splitlines()
    vcf_header = [line for line in vcf_lines if line.startswith("##")]
    vcf_header_line = None
    for line in vcf_lines:
        if line.startswith("#CHROM"):
            vcf_header_line = line
            break
    if vcf_header_line is None:
        sys.exit(f"Error: no #CHROM line found in {original_vcf}")

    # 2. Read multianno header
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

    # 3. Identify Otherinfo column indices
    oi_indices = _find_otherinfo_indices(ma_cols)
    if not oi_indices:
        sys.exit(f"Error: no Otherinfo columns found in {multianno_path}")

    # 4. Build output headers
    out_lines: list[str] = []
    out_lines.extend(vcf_header)

    # Insert ##INFO headers for annotation columns
    out_lines.append(
        '##INFO=<ID=ANNOVAR_DATE,Number=1,Type=String,'
        'Description="Flag the start of ANNOVAR annotation for one '
        'alternative allele">'
    )
    for col in ma_cols:
        if _skip_for_info_header(col):
            continue
        if _FREQ_PREFIX.match(col):
            typ = "Float"
            desc = f"{col} annotation provided by ANNOVAR"
        else:
            typ = "String"
            desc = f"{col} annotation provided by ANNOVAR"
        out_lines.append(
            f'##INFO=<ID={col},Number=.,Type={typ},'
            f'Description="{desc}">'
        )
    out_lines.append(
        '##INFO=<ID=ALLELE_END,Number=0,Type=Flag,'
        'Description="Flag the end of ANNOVAR annotation for one '
        'alternative allele">'
    )
    out_lines.append(vcf_header_line)

    # 5. Process data rows
    # anno_string loop: columns 5 to @name-2 (i.e. index n_ma_cols - 2)
    annot_end_idx = n_ma_cols - 1  # exclusive, so range(5, annot_end_idx)
    date = _anno_date()

    variant_count = 0
    pre_sample_key: tuple | None = None
    pre_anno_string: str = ""
    pre_vcf_cols: list[str] | None = None

    def _flush() -> None:
        nonlocal variant_count
        if pre_vcf_cols is None:
            return
        # Append accumulated annotation string to INFO (index 7)
        pre_vcf_cols[7] += pre_anno_string
        out_lines.append("\t".join(pre_vcf_cols))
        variant_count += 1

    for line in ma_lines[ma_data_start:]:
        stripped = line.rstrip()
        if not stripped:
            continue
        fields = stripped.split("\t")
        if len(fields) < n_ma_cols:
            fields.extend(["."] * (n_ma_cols - len(fields)))

        # Build annotation string (cols 5 .. @name-2)
        anno = _build_anno_string(fields, ma_cols, annot_end_idx, date)

        # Multi-allelic merging: compare sample columns
        cur_key = _sample_key(fields, oi_indices)
        if pre_sample_key is not None and cur_key == pre_sample_key:
            # Same locus, merge annotation string
            pre_anno_string += anno
        else:
            # Flush previous
            _flush()
            # Start new
            pre_vcf_cols = _vcf_columns_from_otherinfo(
                fields, oi_indices, ""
            )
            pre_anno_string = anno
            pre_sample_key = cur_key

    # Flush last row
    _flush()

    Path(output_path).write_text("\n".join(out_lines) + "\n")
    print(f"Done. {variant_count} variants written to {output_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Convert _multianno.txt to VCF format"
    )
    parser.add_argument(
        "multianno",
        help="Input _multianno.txt file (with transvar.input removed)",
    )
    parser.add_argument(
        "original_vcf",
        help="Original input VCF file (for header and sample data)",
    )
    parser.add_argument(
        "-o", "--output",
        help="Output VCF file (default: same name with .vcf extension)",
    )
    args = parser.parse_args()
    process(args.multianno, args.original_vcf, args.output)
