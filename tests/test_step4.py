"""Tests for step4_run_transvar.py."""

from unittest import mock
from pathlib import Path
import sys
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from step4_run_transvar import process


class TestProcess:
    def test_transvar_not_installed(self, tmp_path: Path):
        inp = tmp_path / "transvar.input"
        inp.touch()
        with mock.patch("step4_run_transvar.shutil.which", return_value=None):
            with pytest.raises(SystemExit, match="transvar not found"):
                process(str(inp))

    def test_file_not_found(self):
        with mock.patch("step4_run_transvar.shutil.which", return_value="/usr/bin/transvar"):
            with pytest.raises(SystemExit, match="not found"):
                process("nonexistent.input")

    def test_command_construction(self, tmp_path: Path):
        inp = tmp_path / "transvar.input"
        inp.write_text("chr1:g.100_100insA\n")
        out = tmp_path / "transvar.output"

        with mock.patch("step4_run_transvar.shutil.which", return_value="/usr/bin/transvar"):
            with mock.patch("subprocess.run") as mock_run:
                mock_run.return_value = mock.Mock(returncode=0)
                process(str(inp), str(out))
                mock_run.assert_called_once()
                cmd = mock_run.call_args[0][0]
                assert cmd[0] == "transvar"
                assert cmd[1] == "ganno"
                assert cmd[2] == "-l"
                assert cmd[3] == str(inp)
                assert "--refseq" in cmd

    def test_command_without_refseq(self, tmp_path: Path):
        inp = tmp_path / "transvar.input"
        inp.write_text("chr1:g.100_100insA\n")
        out = tmp_path / "transvar.output"

        with mock.patch("step4_run_transvar.shutil.which", return_value="/usr/bin/transvar"):
            with mock.patch("subprocess.run") as mock_run:
                mock_run.return_value = mock.Mock(returncode=0)
                process(str(inp), str(out), refseq=False)
                cmd = mock_run.call_args[0][0]
                assert "--refseq" not in cmd
