"""Tests for Pipeline module."""

import pytest
from pathlib import Path
from unittest.mock import patch, MagicMock

from myannovar.pipeline import Pipeline

_ANNOVAR_KW = {
    "build": "hg38",
    "humandb": "/data/humandb",
    "protocol": "refGene,avsnp150",
    "operation": "g,f",
    "argument": "'-hgvs',",
}


class TestPipelineInit:
    """Test Pipeline initialization."""

    def test_default_workdir(self):
        p = Pipeline()
        assert p.workdir == Path.cwd()

    def test_custom_workdir(self, tmp_path: Path):
        p = Pipeline(workdir=tmp_path)
        assert p.workdir == tmp_path

    def test_keep_temp(self):
        p = Pipeline(keep_temp=True)
        assert p.keep_temp is True


class TestPipelineCleanup:
    """Test intermediate file cleanup."""

    def test_cleanup_removes_files(self, tmp_path: Path):
        f = tmp_path / "temp.txt"
        f.write_text("test")
        p = Pipeline(workdir=tmp_path)
        p._track(str(f))
        p._cleanup()
        assert not f.exists()
        assert len(p._temp_files) == 0

    def test_cleanup_ignores_missing(self, tmp_path: Path):
        p = Pipeline(workdir=tmp_path)
        p._track(str(tmp_path / "nonexistent.txt"))
        p._cleanup()  # should not raise

    def test_keep_temp_preserves_files(self, tmp_path: Path):
        f = tmp_path / "temp.txt"
        f.write_text("test")
        p = Pipeline(workdir=tmp_path, keep_temp=True)
        p._track(str(f))
        # _cleanup is not called when keep_temp=True
        assert f.exists()


class TestPipelineRunTransvar:
    """Test run_transvar with mocked step functions."""

    @patch("myannovar.pipeline.step1_process")
    @patch("myannovar.pipeline.step2_process")
    @patch("myannovar.pipeline.step3_process")
    @patch("myannovar.pipeline.step4_process")
    @patch("myannovar.pipeline.step5_process")
    @patch("myannovar.pipeline.step6_process")
    @patch("myannovar.pipeline.step7_process")
    def test_run_transvar_calls_all_steps(
        self, s7, s6, s5, s4, s3, s2, s1, tmp_path: Path
    ):
        avinput = str(tmp_path / "sample.avinput")
        multianno = str(tmp_path / "sample_multianno.txt")
        output = str(tmp_path / "result.vcf")

        # Create minimal input files
        Path(avinput).write_text("chr1\t100\t100\tA\tT\n")
        Path(multianno).write_text("#Chr\tStart\tEnd\tRef\tAlt\nchr1\t100\t100\tA\tT\n")

        p = Pipeline(workdir=tmp_path, keep_temp=True)
        result = p.run_transvar(avinput, multianno, output)

        assert result == output
        s1.assert_called_once()
        s2.assert_called_once()
        s3.assert_called_once()
        s4.assert_called_once()
        s5.assert_called_once()
        s6.assert_called_once()
        s7.assert_called_once()


class TestPipelineRunFull:
    """Test run_full with mocked AnnovarRunner and steps."""

    @patch("myannovar.pipeline.step1_process")
    @patch("myannovar.pipeline.step2_process")
    @patch("myannovar.pipeline.step3_process")
    @patch("myannovar.pipeline.step4_process")
    @patch("myannovar.pipeline.step5_process")
    @patch("myannovar.pipeline.step6_process")
    @patch("myannovar.pipeline.step7_process")
    @patch("myannovar.pipeline.AnnovarRunner")
    def test_run_full_calls_annovar_and_steps(
        self, MockRunner, s7, s6, s5, s4, s3, s2, s1, tmp_path: Path
    ):
        input_vcf = str(tmp_path / "input.vcf")
        output_vcf = str(tmp_path / "result.vcf")

        # Create input VCF
        Path(input_vcf).write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\nchr1\t100\t.\tA\tT\n")

        # Mock AnnovarRunner to create expected output files
        mock_runner = MagicMock()
        MockRunner.return_value = mock_runner
        mock_runner.run.return_value = {
            "avinput": str(tmp_path / "result.avinput"),
            "multianno": str(tmp_path / "result_multianno.txt"),
        }

        # Create ANNOVAR output files
        (tmp_path / "result.avinput").write_text("chr1\t100\t100\tA\tT\t-\n")
        (tmp_path / "result_multianno.txt").write_text(
            "#Chr\tStart\tEnd\tRef\tAlt\ttransvar.input\nchr1\t100\t100\tA\tT\t-\n"
        )

        p = Pipeline(workdir=tmp_path, keep_temp=True)
        result = p.run_full(
            input_vcf, output_vcf, **_ANNOVAR_KW,
        )

        assert result == output_vcf
        mock_runner.run.assert_called_once()
