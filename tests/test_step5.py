"""Tests for step5_parse_transvar_output.py."""

import pytest
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from step5_parse_transvar_output import parse_region, parse_coordinates, process


DATA = Path(__file__).resolve().parent.parent / "test_data"


class TestParseRegion:
    def test_exon(self):
        assert parse_region("inside_[cds_in_exon_10]") == "exon10"

    def test_intron(self):
        assert parse_region("inside_[intron_between_exon_2_and_3]") == "intron2"

    def test_unknown(self):
        assert parse_region("upstream") == "upstream"


class TestParseCoordinates:
    def test_full_gdna_cdna_protein(self):
        assert parse_coordinates("chr3:g.178936091G>A/c.1633G>A/p.E545K") == ("c.1633G>A", "p.E545K")

    def test_two_parts(self):
        # Two parts: gDNA/cDNA → returns (gDNA, "-")
        assert parse_coordinates("chr1:g.100_101insA/c.100_101insA") == ("chr1:g.100_101insA", "-")

    def test_single_part(self):
        assert parse_coordinates("c.100A>T") == ("-", "-")


class TestProcess:
    def test_normal(self, tmp_path: Path):
        out = tmp_path / "output.txt"
        process(str(DATA / "test_transvar.output"), str(out))

        lines = out.read_text().strip().splitlines()
        assert len(lines) == 4  # header + 3 variants
        assert lines[0] == "transvar.input\tAAChange.transvar"
        # Merged transcripts
        assert "PIK3CA:NM_006218.2:exon10:c.1633G>A:p.E545K" in lines[1]
        assert "PIK3CA:NM_006218.3:exon10:c.1633G>A:p.E545K" in lines[1]
        # Intron
        assert "intron2" in lines[3]

    def test_file_not_found(self):
        with pytest.raises(SystemExit):
            process("nonexistent.output")
