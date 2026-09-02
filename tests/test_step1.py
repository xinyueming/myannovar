"""Tests for step1_add_transvar_input.py."""

import pytest
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from step1_add_transvar_input import build_transvar_input, process


class TestBuildTransvarInput:
    """Test build_transvar_input() for each rule branch."""

    def test_insertion(self):
        assert build_transvar_input("chr1", "100", "100", "-", "A") == "chr1:g.100_100insA"

    def test_deletion(self):
        assert build_transvar_input("chr1", "100", "102", "G", "-") == "chr1:g.100_102del"

    def test_ref_zero(self):
        assert build_transvar_input("chr1", "100", "100", "0", "T") == "-"

    def test_alt_zero(self):
        assert build_transvar_input("chr1", "100", "100", "A", "0") == "-"

    def test_snv(self):
        assert build_transvar_input("chr1", "100", "101", "A", "T") == "chr1:g.100_101delinsT"

    def test_multinuc_delins(self):
        assert build_transvar_input("chr2", "200", "205", "ACTG", "T") == "chr2:g.200_205delinsT"

    def test_lowercase_seq(self):
        assert build_transvar_input("chrX", "10", "10", "-", "atcg") == "chrX:g.10_10insatcg"

    def test_non_seq_ref(self):
        assert build_transvar_input("chr1", "100", "100", "N", "T") == "-"

    def test_other(self):
        assert build_transvar_input("chr1", "100", "100", ".", ".") == "-"


class TestProcess:
    """Test end-to-end process() function."""

    def test_sample_avinput(self, tmp_path: Path):
        src = Path(__file__).resolve().parent.parent / "test_data" / "sample.avinput"
        out = tmp_path / "output.avinput"
        process(str(src), str(out))

        lines = out.read_text().strip().splitlines()
        # 1 comment + 5 data lines
        assert len(lines) == 6
        assert lines[0].startswith("#")
        assert lines[1].endswith("chr1:g.100_100insA")
        assert lines[2].endswith("chr1:g.100_102del")
        assert lines[3].endswith("-")
        assert lines[4].endswith("chr1:g.100_101delinsT")
        assert lines[5].endswith("chr2:g.200_205delinsT")

    def test_edge_avinput(self, tmp_path: Path):
        src = Path(__file__).resolve().parent.parent / "test_data" / "edge.avinput"
        out = tmp_path / "output.avinput"
        process(str(src), str(out))

        lines = out.read_text().strip().splitlines()
        assert len(lines) == 6
        assert lines[1].endswith("-")  # ref=0
        assert lines[2].endswith("chr2:g.100_100inst")  # lowercase
        assert lines[3].endswith("-")  # only 4 columns
        assert lines[4].endswith("-")  # non-seq ref N
        assert lines[5].endswith("chr5:g.400_400del")  # alt=-

    def test_file_not_found(self):
        with pytest.raises(SystemExit):
            process("nonexistent.avinput")
