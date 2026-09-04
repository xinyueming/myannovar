"""Tests for AnnovarRunner module."""

import pytest
from pathlib import Path
from unittest.mock import patch, MagicMock

from myannovar.annovar import AnnovarRunner


class TestAnnovarRunnerDirResolution:
    """Test ANNOVAR directory resolution."""

    def test_explicit_dir(self, tmp_path: Path):
        script = tmp_path / "table_annovar.pl"
        script.touch()
        runner = AnnovarRunner(annovar_dir=str(tmp_path))
        assert runner.annovar_dir == tmp_path.resolve()

    def test_explicit_dir_not_found(self, tmp_path: Path):
        with pytest.raises(FileNotFoundError, match="table_annovar.pl not found"):
            AnnovarRunner(annovar_dir=str(tmp_path))

    def test_builtin_dir(self):
        """Built-in annovar_scripts/ should be found."""
        runner = AnnovarRunner()
        assert runner.annovar_dir.name == "annovar_scripts"


class TestAnnovarRunnerRun:
    """Test AnnovarRunner.run()."""

    _KW = {
        "build": "hg38",
        "humandb": "/data/humandb",
        "protocol": "refGene,avsnp150",
        "operation": "g,f",
        "argument": "'-hgvs',",
    }

    def test_perl_not_found(self, tmp_path: Path):
        script = tmp_path / "table_annovar.pl"
        script.touch()
        runner = AnnovarRunner(annovar_dir=str(tmp_path))

        with patch("myannovar.annovar.shutil.which", return_value=None):
            with pytest.raises(FileNotFoundError, match="perl not found"):
                runner.run("input.vcf", "sample.var", **self._KW)

    def test_command_construction(self, tmp_path: Path):
        """Verify the command line matches real table_annovar.pl usage."""
        script = tmp_path / "table_annovar.pl"
        script.touch()
        runner = AnnovarRunner(annovar_dir=str(tmp_path))

        with patch("myannovar.annovar.shutil.which", return_value="/usr/bin/perl"):
            with patch("myannovar.annovar.subprocess.run") as mock_run:
                mock_run.return_value = MagicMock(returncode=0)
                runner.run("input.vcf", "sample.var", **self._KW)

                call_args = mock_run.call_args[0][0]
                assert call_args[0] == "perl"
                assert str(script) in call_args
                assert "input.vcf" in call_args
                # humandb as positional after input
                assert "/data/humandb" in call_args
                assert "-buildver" in call_args
                assert "hg38" in call_args
                assert "-out" in call_args
                assert "sample.var" in call_args
                # Required flags
                assert "-remove" in call_args
                assert "-vcfinput" in call_args
                assert "--polish" in call_args
                # No removed flags
                assert "-csvlog" not in call_args
                assert "-xref" not in call_args
                # Custom values passed through
                assert "-protocol" in call_args
                assert "refGene,avsnp150" in call_args
                assert "g,f" in call_args

    def test_argument_passed_verbatim(self, tmp_path: Path):
        """--argument value should be passed exactly as provided."""
        script = tmp_path / "table_annovar.pl"
        script.touch()
        runner = AnnovarRunner(annovar_dir=str(tmp_path))

        with patch("myannovar.annovar.shutil.which", return_value="/usr/bin/perl"):
            with patch("myannovar.annovar.subprocess.run") as mock_run:
                mock_run.return_value = MagicMock(returncode=0)
                runner.run(
                    "input.vcf", "sample.var",
                    build="hg38", humandb="/db",
                    protocol="db1,db2,db3", operation="g,r,f",
                    argument="'-hgvs',,",
                )

                call_args = mock_run.call_args[0][0]
                idx = call_args.index("--argument")
                assert call_args[idx + 1] == "'-hgvs',,"

    def test_run_failure(self, tmp_path: Path):
        """Should raise on non-zero return code."""
        script = tmp_path / "table_annovar.pl"
        script.touch()
        runner = AnnovarRunner(annovar_dir=str(tmp_path))

        with patch("myannovar.annovar.shutil.which", return_value="/usr/bin/perl"):
            with patch("myannovar.annovar.subprocess.run") as mock_run:
                mock_run.return_value = MagicMock(returncode=1, stderr="error")
                with pytest.raises(RuntimeError, match="table_annovar.pl failed"):
                    runner.run("input.vcf", "sample.var", **self._KW)
