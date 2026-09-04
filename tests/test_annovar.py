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
        # Should resolve to the installed annovar_scripts/
        assert runner.annovar_dir.name == "annovar_scripts"


class TestAnnovarRunnerRun:
    """Test AnnovarRunner.run()."""

    def test_perl_not_found(self, tmp_path: Path):
        """Should raise if perl is not in PATH."""
        script = tmp_path / "table_annovar.pl"
        script.touch()
        runner = AnnovarRunner(annovar_dir=str(tmp_path))

        with patch("myannovar.annovar.shutil.which", return_value=None):
            with pytest.raises(FileNotFoundError, match="perl not found"):
                runner.run("input.vcf", "sample", ["refGene"])

    def test_command_construction(self, tmp_path: Path):
        """Verify the command line is built correctly."""
        script = tmp_path / "table_annovar.pl"
        script.touch()
        runner = AnnovarRunner(annovar_dir=str(tmp_path))

        with patch("myannovar.annovar.shutil.which", return_value="/usr/bin/perl"):
            with patch("myannovar.annovar.subprocess.run") as mock_run:
                mock_run.return_value = MagicMock(returncode=0)
                runner.run("input.vcf", "sample", ["refGene", "avsnp150"], "hg19")

                call_args = mock_run.call_args[0][0]
                assert "perl" in call_args
                assert str(script) in call_args
                assert "input.vcf" in call_args
                assert "-buildver" in call_args
                assert "hg19" in call_args
                assert "-protocol" in call_args
                assert "refGene,avsnp150" in call_args

    def test_run_failure(self, tmp_path: Path):
        """Should raise on non-zero return code."""
        script = tmp_path / "table_annovar.pl"
        script.touch()
        runner = AnnovarRunner(annovar_dir=str(tmp_path))

        with patch("myannovar.annovar.shutil.which", return_value="/usr/bin/perl"):
            with patch("myannovar.annovar.subprocess.run") as mock_run:
                mock_run.return_value = MagicMock(returncode=1, stderr="error")
                with pytest.raises(RuntimeError, match="table_annovar.pl failed"):
                    runner.run("input.vcf", "sample", ["refGene"])
