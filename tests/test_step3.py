"""Tests for step3_extract_transvar_input.py."""

import pytest
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from step3_extract_transvar_input import process


DATA = Path(__file__).resolve().parent.parent / "test_data"


class TestProcess:
    def test_normal_extraction(self, tmp_path: Path):
        # Use step2 output (with transvar.input column) as input
        out = tmp_path / "output.txt"
        process(str(DATA / "step2_output.txt"), str(out))
        lines = Path(out).read_text().strip().splitlines()
        assert len(lines) == 3

    def test_no_valid_entries(self, tmp_path: Path):
        src = tmp_path / "all_dash.txt"
        src.write_text("# Chr\tStart\tEnd\tRef\tAlt\ttransvar.input\nchr1\t100\t100\tA\tT\t-\n")
        out = tmp_path / "output.txt"
        process(str(src), str(out))
        lines = Path(out).read_text().strip().splitlines()
        assert len(lines) == 0

    def test_dedup(self, tmp_path: Path):
        src = tmp_path / "dup.txt"
        src.write_text(
            "# Chr\tStart\tEnd\tRef\tAlt\ttransvar.input\n"
            "chr1\t100\t100\t-\tA\tchr1:g.100_100insA\n"
            "chr1\t100\t100\t-\tA\tchr1:g.100_100insA\n"
            "chr2\t200\t200\tG\t-\tchr2:g.200_200del\n"
        )
        out = tmp_path / "output.txt"
        process(str(src), str(out))
        lines = Path(out).read_text().strip().splitlines()
        assert len(lines) == 2

    def test_file_not_found(self):
        with pytest.raises(SystemExit):
            process("nonexistent.txt")
