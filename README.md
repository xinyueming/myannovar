# MyANNOVAR

ANNOVAR + TransVar variant annotation pipeline.

## Installation

```bash
pip install -e .
```

### Requirements

- Python >= 3.8
- Perl (for ANNOVAR)
- TransVar (`pip install transvar`)

## Quick Start

```bash
# One-command annotation (ANNOVAR + TransVar)
myannovar run -i input.vcf -o result.vcf \
  --humandb /path/to/humandb \
  --protocol "refGeneWithVer" \
  --operation "g" \
  --argument "'-hgvs'"

# Step-by-step: ANNOVAR only
myannovar annovar -i input.vcf -o sample.var \
  --humandb /path/to/humandb \
  --protocol "refGeneWithVer" \
  --operation "g" \
  --argument "'-hgvs'"

# Step-by-step: TransVar annotation
myannovar annotate -i sample.avinput -m sample_multianno.txt -o result.vcf
```

## Commands

### `myannovar run`

Run the complete ANNOVAR + TransVar pipeline in one command.

```bash
myannovar run -i input.vcf -o result.vcf \
  --humandb /path/to/humandb \
  --protocol "refGeneWithVer,cytoBand,clinvar_20220320" \
  --operation "g,r,f,f,f" \
  --argument "'-hgvs',,,"
```

| Option | Description |
|--------|-------------|
| `-i, --input` | Input VCF file |
| `-o, --output` | Output VCF file |
| `--humandb` | ANNOVAR human database directory |
| `-b, --build` | Genome build (default: hg38) |
| `--protocol` | Protocol string (e.g. refGeneWithVer,cytoBand) |
| `--operation` | Operation string (e.g. g,r,f) |
| `--argument` | Argument string (e.g. '-hgvs',,) |
| `--annovar-dir` | ANNOVAR scripts directory (optional, auto-detected) |
| `--keep-temp` | Keep intermediate files |
| `--refseq / --no-refseq` | Enable/disable RefSeq annotations |
| `-v, --verbose` | Verbose output |

### `myannovar annovar`

Run ANNOVAR table_annovar.pl only. Produces `.avinput` and `_multianno.txt`.

```bash
myannovar annovar -i input.vcf -o sample.var \
  --humandb /path/to/humandb \
  --protocol "refGeneWithVer,cytoBand,clinvar_20220320" \
  --operation "g,r,f,f,f" \
  --argument "'-hgvs',,,"
```

| Option | Description |
|--------|-------------|
| `-i, --input` | Input VCF file |
| `-o, --output` | Output prefix (e.g. sample.var) |
| `--humandb` | ANNOVAR human database directory |
| `-b, --build` | Genome build (default: hg38) |
| `--protocol` | Protocol string (e.g. refGeneWithVer,cytoBand) |
| `--operation` | Operation string (e.g. g,r,f) |
| `--argument` | Argument string (e.g. '-hgvs',,) |
| `--annovar-dir` | ANNOVAR scripts directory (optional) |
| `-v, --verbose` | Verbose output |

### `myannovar annotate`

Run the TransVar 7-step annotation pipeline on existing ANNOVAR output.

```bash
myannovar annotate -i sample.avinput -m sample_multianno.txt -o result.vcf
```

| Option | Description |
|--------|-------------|
| `-i, --input` | Input .avinput file |
| `-m, --multianno` | Input _multianno.txt file |
| `-o, --output` | Output VCF file |
| `--keep-temp` | Keep intermediate files |
| `--refseq / --no-refseq` | Enable/disable RefSeq annotations |
| `-v, --verbose` | Verbose output |

## Pipeline Steps

The `annotate` and `run` commands execute these steps internally:

| Step | Description |
|------|-------------|
| 1 | Add transvar.input column to .avinput |
| 2 | Merge transvar.input into _multianno.txt |
| 3 | Extract unique transvar.input values |
| 4 | Run TransVar ganno |
| 5 | Parse TransVar output into transvar.multianno |
| 6 | Replace AAChange.refGeneWithVer with TransVar annotations |
| 7 | Convert _multianno.txt to VCF format |

## Input/Output Format

### Input

Standard VCF 4.2 format.

### Output

VCF 4.2 format with ANNOVAR annotations in the INFO field.

```vcf
##fileformat=VCFv4.2
##INFO=<ID=ANN,Number=.,Type=String,Description="ANNOVAR annotations">
#CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO
chr1    100  .   A    T    .     .       Start=100;End=100;...
```

## Configuration

### ANNOVAR Scripts Directory

The pipeline auto-detects ANNOVAR scripts in this order:

1. `--annovar-dir` CLI option
2. `ANNOVAR_DIR` environment variable
3. Built-in `annovar_scripts/` (bundled with the package)

### Environment Variables

| Variable | Description |
|----------|-------------|
| `ANNOVAR_DIR` | Path to ANNOVAR scripts directory |

## Examples

```bash
# Multiple databases
myannovar run -i sample.vcf -o result.vcf \
  --humandb /path/to/humandb \
  --protocol "refGeneWithVer,cytoBand,gnomad_genome" \
  --operation "g,r,f" \
  --argument "'-hgvs',,"

# Use hg19 build
myannovar run -i sample.vcf -o result.vcf \
  --humandb /path/to/humandb \
  --protocol "refGeneWithVer" \
  --operation "g" \
  --argument "'-hgvs'" \
  -b hg19

# Keep intermediate files for debugging
myannovar run -i sample.vcf -o result.vcf \
  --humandb /path/to/humandb \
  --protocol "refGeneWithVer" \
  --operation "g" \
  --argument "'-hgvs'" \
  --keep-temp

# Verbose logging
myannovar run -i sample.vcf -o result.vcf \
  --humandb /path/to/humandb \
  --protocol "refGeneWithVer" \
  --operation "g" \
  --argument "'-hgvs'" \
  -v
```

## Troubleshooting

**`transvar not found`**: Install TransVar first: `pip install transvar`

**`perl not found`**: Install Perl: `apt install perl` (Ubuntu) or `brew install perl` (macOS)

**`table_annovar.pl failed`**: Check that your input VCF is valid and ANNOVAR databases are downloaded

**`ANNOVAR scripts not found`**: Set `ANNOVAR_DIR` environment variable or use `--annovar-dir`

## Programmatic Usage

```python
from myannovar.annovar import AnnovarRunner
from myannovar.pipeline import Pipeline
from myannovar.steps.step1 import process as step1_process

# Use AnnovarRunner directly
runner = AnnovarRunner()
result = runner.run(
    input_file="input.vcf",
    output_prefix="sample",
    build="hg38",
    humandb="/path/to/humandb",
    protocol="refGeneWithVer",
    operation="g",
    argument="'-hgvs'",
)

# Use Pipeline for full orchestration
pipeline = Pipeline(keep_temp=False)
pipeline.run_full(
    input_vcf="input.vcf",
    output_vcf="result.vcf",
    build="hg38",
    humandb="/path/to/humandb",
    protocol="refGeneWithVer",
    operation="g",
    argument="'-hgvs'",
)

# Use individual steps
step1_process("input.avinput", "output.avinput")
```

## Project Structure

```
myannovar/
├── pyproject.toml
├── README.md
├── src/
│   └── myannovar/
│       ├── __init__.py
│       ├── annovar.py          # AnnovarRunner
│       ├── cli.py              # CLI entry point
│       ├── pipeline.py         # Pipeline orchestration
│       ├── logging_config.py
│       ├── steps/              # TransVar steps
│       │   ├── step1.py
│       │   ├── step2.py
│       │   └── ...
│       └── annovar_scripts/    # ANNOVAR Perl scripts
├── tests/
└── test_data/
```

## License

MIT