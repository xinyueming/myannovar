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
myannovar run -i input.vcf -o result.vcf --db refGene --build hg38

# Step-by-step: ANNOVAR only
myannovar annovar -i input.vcf -o sample --db refGene --build hg38

# Step-by-step: TransVar annotation
myannovar annotate -i sample.avinput -m sample_multianno.txt -o result.vcf
```

## Commands

### `myannovar run`

Run the complete ANNOVAR + TransVar pipeline in one command.

```bash
myannovar run -i input.vcf -o result.vcf --db refGene --db avsnp150 --build hg38
```

| Option | Description |
|--------|-------------|
| `-i, --input` | Input VCF file |
| `-o, --output` | Output VCF file |
| `-d, --db` | Annotation database (can be specified multiple times) |
| `-b, --build` | Genome build (default: hg38) |
| `--annovar-dir` | ANNOVAR scripts directory (optional, auto-detected) |
| `--keep-temp` | Keep intermediate files |
| `--refseq / --no-refseq` | Enable/disable RefSeq annotations |
| `-v, --verbose` | Verbose output |

### `myannovar annovar`

Run ANNOVAR table_annovar.pl only. Produces `.avinput` and `_multianno.txt`.

```bash
myannovar annovar -i input.vcf -o sample --db refGene --build hg38
```

### `myannovar annotate`

Run the TransVar 7-step annotation pipeline on existing ANNOVAR output.

```bash
myannovar annotate -i sample.avinput -m sample_multianno.txt -o result.vcf
```

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
myannovar run -i sample.vcf -o result.vcf -d refGene -d avsnp150 -d gnomad_genome

# Use hg19 build
myannovar run -i sample.vcf -o result.vcf -d refGene -b hg19

# Keep intermediate files for debugging
myannovar run -i sample.vcf -o result.vcf -d refGene --keep-temp

# Verbose logging
myannovar run -i sample.vcf -o result.vcf -d refGene -v
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
from myannovar.steps import step1_process, step7_process

# Use AnnovarRunner directly
runner = AnnovarRunner()
result = runner.run("input.vcf", "sample", ["refGene"], "hg38")

# Use Pipeline for full orchestration
pipeline = Pipeline(keep_temp=False)
pipeline.run_full("input.vcf", "result.vcf", ["refGene"], "hg38")

# Use individual steps
step1_process("input.avinput", "output.avinput")
```
