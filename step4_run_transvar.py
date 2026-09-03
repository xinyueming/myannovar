#!/usr/bin/env python3
"""Step 4: Run TransVar genomic annotation.

Calls: transvar ganno -l transvar.input --refseq > transvar.output
"""

import sys
import shutil
import argparse
import subprocess
from pathlib import Path
from typing import Optional


def process(input_path: str, output_path: Optional[str] = None, refseq: bool = True) -> None:
    if not shutil.which("transvar"):
        sys.exit("Error: transvar not found. Install: pip install transvar")

    if not Path(input_path).is_file():
        sys.exit(f"Error: {input_path} not found")

    if output_path is None:
        output_path = Path(input_path).parent / "transvar.output"
    else:
        output_path = Path(output_path)

    out_dir = output_path.parent
    if not out_dir.is_dir():
        sys.exit(f"Error: output directory {out_dir} does not exist")

    cmd = ["transvar", "ganno", "-l", input_path]
    if refseq:
        cmd.append("--refseq")

    print(f"Running: {' '.join(cmd)}")

    with Path(output_path).open("w") as out:
        result = subprocess.run(cmd, stdout=out, stderr=sys.stderr)

    if result.returncode != 0:
        sys.exit(f"Error: transvar exited with code {result.returncode}")

    lines = Path(output_path).read_text().strip().splitlines()
    print(f"Done. {len(lines)} annotations written to {output_path}")


def main():
    parser = argparse.ArgumentParser(description="Run TransVar genomic annotation")
    parser.add_argument("input", help="Input transvar.input file")
    parser.add_argument("-o", "--output", help="Output file (default: transvar.output next to input)")
    parser.add_argument("--no-refseq", action="store_true", help="Skip --refseq flag")
    args = parser.parse_args()
    process(args.input, args.output, refseq=not args.no_refseq)


if __name__ == "__main__":
    main()
