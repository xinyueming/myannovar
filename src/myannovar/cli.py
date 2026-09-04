"""CLI entry point for myannovar."""

import argparse
import logging
from pathlib import Path
from typing import List, Optional

from myannovar.annovar import AnnovarRunner
from myannovar.logging_config import setup_logging
from myannovar.pipeline import Pipeline


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="myannovar",
        description="ANNOVAR + TransVar variant annotation pipeline",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    # -- annovar -------------------------------------------------------
    p_annovar = sub.add_parser(
        "annovar",
        help="Run ANNOVAR table_annovar.pl",
    )
    p_annovar.add_argument("-i", "--input", required=True, help="Input VCF file")
    p_annovar.add_argument("-o", "--output", required=True, help="Output prefix")
    p_annovar.add_argument(
        "-d", "--db", required=True, action="append",
        help="Annotation database (can be specified multiple times)",
    )
    p_annovar.add_argument("-b", "--build", default="hg38", help="Genome build (default: hg38)")
    p_annovar.add_argument("--annovar-dir", help="ANNOVAR scripts directory")
    p_annovar.add_argument("-v", "--verbose", action="store_true", help="Verbose output")
    p_annovar.set_defaults(func=_cmd_annovar)

    # -- annotate ------------------------------------------------------
    p_annotate = sub.add_parser(
        "annotate",
        help="Run TransVar 7-step annotation pipeline",
    )
    p_annotate.add_argument("-i", "--input", required=True, help="Input .avinput file")
    p_annotate.add_argument("-m", "--multianno", required=True, help="Input _multianno.txt file")
    p_annotate.add_argument("-o", "--output", required=True, help="Output VCF file")
    p_annotate.add_argument(
        "--keep-temp", action="store_true",
        help="Keep intermediate files",
    )
    p_annotate.add_argument(
        "--refseq", action="store_true", default=True,
        help="Use RefSeq annotations (default: enabled)",
    )
    p_annotate.add_argument(
        "--no-refseq", action="store_true",
        help="Disable RefSeq annotations",
    )
    p_annotate.add_argument("-v", "--verbose", action="store_true", help="Verbose output")
    p_annotate.set_defaults(func=_cmd_annotate)

    # -- run -----------------------------------------------------------
    p_run = sub.add_parser(
        "run",
        help="Run full ANNOVAR + TransVar pipeline",
    )
    p_run.add_argument("-i", "--input", required=True, help="Input VCF file")
    p_run.add_argument("-o", "--output", required=True, help="Output VCF file")
    p_run.add_argument(
        "-d", "--db", required=True, action="append",
        help="Annotation database (can be specified multiple times)",
    )
    p_run.add_argument("-b", "--build", default="hg38", help="Genome build (default: hg38)")
    p_run.add_argument("--annovar-dir", help="ANNOVAR scripts directory")
    p_run.add_argument(
        "--keep-temp", action="store_true",
        help="Keep intermediate files",
    )
    p_run.add_argument(
        "--refseq", action="store_true", default=True,
        help="Use RefSeq annotations (default: enabled)",
    )
    p_run.add_argument(
        "--no-refseq", action="store_true",
        help="Disable RefSeq annotations",
    )
    p_run.add_argument("-v", "--verbose", action="store_true", help="Verbose output")
    p_run.set_defaults(func=_cmd_run)

    return parser


# -- command handlers --------------------------------------------------

def _cmd_annovar(args: argparse.Namespace) -> None:
    setup_logging(args.verbose)
    logger = logging.getLogger(__name__)

    runner = AnnovarRunner(annovar_dir=args.annovar_dir)
    logger.info("Running ANNOVAR: %s -> %s", args.input, args.output)
    result = runner.run(
        input_file=args.input,
        output_prefix=args.output,
        db_names=args.db,
        build=args.build,
    )
    logger.info("avinput:  %s", result["avinput"])
    logger.info("multianno: %s", result["multianno"])


def _cmd_annotate(args: argparse.Namespace) -> None:
    setup_logging(args.verbose)
    logger = logging.getLogger(__name__)

    refseq = not args.no_refseq
    pipeline = Pipeline(
        workdir=Path(args.input).parent,
        keep_temp=args.keep_temp,
    )
    output = pipeline.run_transvar(
        avinput=args.input,
        multianno=args.multianno,
        output=args.output,
        refseq=refseq,
    )
    logger.info("Annotation complete: %s", output)


def _cmd_run(args: argparse.Namespace) -> None:
    setup_logging(args.verbose)
    logger = logging.getLogger(__name__)

    refseq = not args.no_refseq
    pipeline = Pipeline(
        workdir=Path(args.input).parent,
        annovar_dir=args.annovar_dir,
        keep_temp=args.keep_temp,
    )
    output = pipeline.run_full(
        input_vcf=args.input,
        output_vcf=args.output,
        db=args.db,
        build=args.build,
        refseq=refseq,
    )
    logger.info("Pipeline complete: %s", output)


# -- main ------------------------------------------------------------

def main(argv: Optional[List[str]] = None) -> None:
    parser = _build_parser()
    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
