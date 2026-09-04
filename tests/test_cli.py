"""Tests for CLI entry point."""

import pytest
from myannovar.cli import _build_parser


class TestParser:
    """Test CLI argument parsing."""

    def setup_method(self):
        self.parser = _build_parser()

    # -- annovar subcommand --------------------------------------------

    def test_annovar_required_args(self):
        args = self.parser.parse_args([
            "annovar", "-i", "input.vcf", "-o", "sample",
            "-d", "refGene",
        ])
        assert args.command == "annovar"
        assert args.input == "input.vcf"
        assert args.output == "sample"
        assert args.db == ["refGene"]
        assert args.build == "hg38"
        assert args.annovar_dir is None
        assert not args.verbose

    def test_annovar_multiple_dbs(self):
        args = self.parser.parse_args([
            "annovar", "-i", "input.vcf", "-o", "sample",
            "-d", "refGene", "-d", "avsnp150",
        ])
        assert args.db == ["refGene", "avsnp150"]

    def test_annovar_build(self):
        args = self.parser.parse_args([
            "annovar", "-i", "input.vcf", "-o", "sample",
            "-d", "refGene", "-b", "hg19",
        ])
        assert args.build == "hg19"

    def test_annovar_verbose(self):
        args = self.parser.parse_args([
            "annovar", "-i", "input.vcf", "-o", "sample",
            "-d", "refGene", "-v",
        ])
        assert args.verbose

    # -- annotate subcommand -------------------------------------------

    def test_annotate_required_args(self):
        args = self.parser.parse_args([
            "annotate", "-i", "sample.avinput", "-m", "sample_multianno.txt",
            "-o", "result.vcf",
        ])
        assert args.command == "annotate"
        assert args.input == "sample.avinput"
        assert args.multianno == "sample_multianno.txt"
        assert args.output == "result.vcf"
        assert not args.keep_temp
        assert args.refseq is True

    def test_annotate_keep_temp(self):
        args = self.parser.parse_args([
            "annotate", "-i", "sample.avinput", "-m", "sample_multianno.txt",
            "-o", "result.vcf", "--keep-temp",
        ])
        assert args.keep_temp

    def test_annotate_no_refseq(self):
        args = self.parser.parse_args([
            "annotate", "-i", "sample.avinput", "-m", "sample_multianno.txt",
            "-o", "result.vcf", "--no-refseq",
        ])
        assert args.no_refseq

    # -- run subcommand ------------------------------------------------

    def test_run_required_args(self):
        args = self.parser.parse_args([
            "run", "-i", "input.vcf", "-o", "result.vcf",
            "-d", "refGene",
        ])
        assert args.command == "run"
        assert args.input == "input.vcf"
        assert args.output == "result.vcf"
        assert args.db == ["refGene"]
        assert args.build == "hg38"
        assert not args.keep_temp
        assert args.refseq is True

    def test_run_all_args(self):
        args = self.parser.parse_args([
            "run", "-i", "input.vcf", "-o", "result.vcf",
            "-d", "refGene", "-d", "avsnp150",
            "-b", "hg19", "--annovar-dir", "/opt/annovar",
            "--keep-temp", "--no-refseq", "-v",
        ])
        assert args.db == ["refGene", "avsnp150"]
        assert args.build == "hg19"
        assert args.annovar_dir == "/opt/annovar"
        assert args.keep_temp
        assert args.no_refseq
        assert args.verbose

    # -- help / error --------------------------------------------------

    def test_no_command_raises(self):
        with pytest.raises(SystemExit):
            self.parser.parse_args([])
