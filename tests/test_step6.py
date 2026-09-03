"""Tests for step6_replace_aachange.py."""

import pytest
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from step6_replace_aachange import load_multianno, find_column_index, process


DATA = Path(__file__).resolve().parent.parent / "test_data"


class TestLoadMultianno:
    def test_load(self):
        lookup = load_multianno(str(DATA / "step6_transvar_multianno.txt"))
        # Note: header line "transvar.input" is also loaded as a key
        assert len(lookup) == 3
        assert "PIK3CA:NM_006218.2:exon10:c.1633G>A:p.E545K" in lookup["chr3:g.178936091_178936091delinsA"]


class TestFindColumnIndex:
    def test_found(self):
        assert find_column_index(["Chr", "Start", "AAChange.refGeneWithVer"], "AAChange.refGeneWithVer") == 2

    def test_not_found(self):
        assert find_column_index(["Chr", "Start"], "Missing") is None

    def test_multiple_names(self):
        assert find_column_index(["Chr", "AAChange.refGene"], "AAChange.refGeneWithVer", "AAChange.refGene") == 1


class TestProcess:
    def test_normal_replacement(self, tmp_path: Path):
        out = tmp_path / "output.txt"
        process(
            str(DATA / "step6_multianno.txt"),
            str(DATA / "step6_transvar_multianno.txt"),
            str(out),
        )
        lines = out.read_text().strip().splitlines()
        assert len(lines) == 4  # header + 3 data
        # First two rows should have transvar annotations
        assert "exon10" in lines[1]
        assert "exon2" in lines[2]
        # Third row should be unchanged (no match)
        assert "R167*" in lines[3]

    def test_file_not_found_multianno(self):
        with pytest.raises(SystemExit):
            process("nonexistent.txt", str(DATA / "step6_transvar_multianno.txt"))

    def test_file_not_found_transvar(self):
        with pytest.raises(SystemExit):
            process(str(DATA / "step6_multianno.txt"), "nonexistent.txt")
