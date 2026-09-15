"""Tests for step5_parse_transvar_output.py."""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from myannovar.steps.step5 import (
    parse_coordinates,
    parse_line,
    parse_region,
    process,
)


class TestParseRegion:
    def test_exon(self):
        assert parse_region("inside_[cds_in_exon_10]") == "exon10"

    def test_intron(self):
        assert parse_region("inside_[intron_between_exon_2_and_3]") == "intron2"

    def test_unknown(self):
        assert parse_region("upstream") == "upstream"


class TestParseCoordinates:
    def test_full_gdna_cdna_protein(self):
        assert parse_coordinates("chr3:g.178936091G>A/c.1633G>A/p.E545K") == (
            "c.1633G>A",
            "p.E545K",
        )

    def test_two_parts(self):
        assert parse_coordinates("chr1:g.100_101insA/c.100_101insA") == (
            "chr1:g.100_101insA",
            "-",
        )

    def test_single_part(self):
        assert parse_coordinates("c.100A>T") == ("-", "-")


class TestParseLine:
    def test_basic_7_columns(self):
        line = "chr10:g.100T>C\tNM_001\tGENE1\t-\tchr10:g.100T>C/c.100T>C/p.X100Y\tinside_[exon_1]\tCSQN=SNV"
        result = parse_line(line)
        assert result is not None
        assert result["input"] == "chr10:g.100T>C"
        assert result["gene"] == "GENE1"
        assert result["chrom"] == ""

    def test_with_gseq_11_columns(self):
        line = "chr10:g.100T>C\tNM_001\tGENE1\t-\tchr10:g.100T>C/c.100T>C/p.X100Y\tinside_[exon_1]\tCSQN=SNV\tchr10\t100\tT\tC"
        result = parse_line(line)
        assert result is not None
        assert result["chrom"] == "chr10"
        assert result["pos"] == "100"
        assert result["ref"] == "T"
        assert result["alt"] == "C"


class TestProcess:
    def _write_transvar_output(self, path: Path, with_gseq: bool = False) -> Path:
        if with_gseq:
            path.write_text(
                "input\ttranscript\tgene\tstrand\tcoordinates\tregion\tinfo\tCHROM\tPOS\tREF\tALT\n"
                "chr10:g.100T>C\tNM_001 (protein_coding)\tGENE1\t-\tchr10:g.100T>C/c.100T>C/p.X100Y\tinside_[exon_1]\tCSQN=SNV\tchr10\t100\tT\tC\n"
                "chr10:g.100T>C\tNM_002 (protein_coding)\tGENE1\t-\tchr10:g.100T>C/c.100T>C/p.X100Y\tinside_[exon_2]\tCSQN=SNV\tchr10\t100\tT\tC\n"
                "chr3:g.200A>G\tNM_003 (protein_coding)\tGENE2\t+\tchr3:g.200A>G/c.200A>G/p.K200R\tinside_[intron_between_exon_2_and_3]\tCSQN=Intronic\tchr3\t200\tA\tG\n"
            )
        else:
            path.write_text(
                "input\ttranscript\tgene\tstrand\tcoordinates\tregion\tinfo\n"
                "chr10:g.100T>C\tNM_001 (protein_coding)\tGENE1\t-\tchr10:g.100T>C/c.100T>C/p.X100Y\tinside_[exon_1]\tCSQN=SNV\n"
                "chr10:g.100T>C\tNM_002 (protein_coding)\tGENE1\t-\tchr10:g.100T>C/c.100T>C/p.X100Y\tinside_[exon_2]\tCSQN=SNV\n"
                "chr3:g.200A>G\tNM_003 (protein_coding)\tGENE2\t+\tchr3:g.200A>G/c.200A>G/p.K200R\tinside_[intron_between_exon_2_and_3]\tCSQN=Intronic\n"
            )
        return path

    def test_without_gseq(self, tmp_path: Path):
        inp = self._write_transvar_output(tmp_path / "tv.output", with_gseq=False)
        out = tmp_path / "transvar.multianno"
        process(str(inp), str(out))
        lines = out.read_text().strip().splitlines()
        assert lines[0] == "transvar.input\tAAChange.transvar"
        assert len(lines) == 3  # header + 2 unique inputs
        assert "GENE1:NM_001:exon1" in lines[1]
        assert "GENE1:NM_002:exon2" in lines[1]
        assert "intron2" in lines[2]

    def test_with_gseq(self, tmp_path: Path):
        inp = self._write_transvar_output(tmp_path / "tv.output", with_gseq=True)
        out = tmp_path / "transvar.multianno"
        process(str(inp), str(out))
        lines = out.read_text().strip().splitlines()
        assert lines[0] == "transvar.input\tAAChange.transvar\tCHROM\tPOS\tREF\tALT"
        assert len(lines) == 3
        # Check gseq values
        parts = lines[1].split("\t")
        assert parts[2] == "chr10"
        assert parts[3] == "100"
        assert parts[4] == "T"
        assert parts[5] == "C"

    def test_file_not_found(self):
        with pytest.raises(SystemExit):
            process("nonexistent.txt")
