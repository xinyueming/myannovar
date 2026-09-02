"""Tests for step2_merge_to_multianno.py."""

import pytest
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from step2_merge_to_multianno import load_avinput, process


DATA = Path(__file__).resolve().parent.parent / "test_data"


class TestLoadAvinput:
    def test_load_valid_avinput(self):
        lookup = load_avinput(str(DATA / "step2.avinput"))
        assert len(lookup) == 3
        assert lookup[("chr1", "100", "100", "-", "A")] == "chr1:g.100_100insA"
        assert lookup[("chr1", "100", "102", "G", "-")] == "chr1:g.100_102del"
        assert lookup[("chr1", "100", "101", "A", "T")] == "chr1:g.100_101delinsT"


class TestProcess:
    def test_normal_match(self, tmp_path: Path):
        out = tmp_path / "out.txt"
        process(str(DATA / "step2.avinput"), str(DATA / "test_multianno.txt"), str(out))

        lines = out.read_text().strip().splitlines()
        # 1 header + 4 data lines
        assert len(lines) == 5
        # header has transvar.input column
        assert lines[0].endswith("transvar.input")
        # matched rows have transvar.input values
        assert lines[1].endswith("chr1:g.100_100insA")
        assert lines[2].endswith("chr1:g.100_102del")
        # unmatched row has "-"
        assert lines[3].endswith("-")
        assert lines[4].endswith("chr1:g.100_101delinsT")

    def test_file_not_found_avinput(self):
        with pytest.raises(SystemExit):
            process("nonexistent.avinput", str(DATA / "test_multianno.txt"))

    def test_file_not_found_multianno(self):
        with pytest.raises(SystemExit):
            process(str(DATA / "step2.avinput"), "nonexistent.txt")
