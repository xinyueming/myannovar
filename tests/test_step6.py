"""Tests for step6_replace_aachange.py."""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from myannovar.steps.step6 import (
    find_column_index,
    load_multianno,
    process,
)


class TestLoadMultianno:
    def test_load_without_gseq(self, tmp_path: Path):
        ma = tmp_path / "tv.ma"
        ma.write_text(
            "transvar.input\tAAChange.transvar\n"
            "chr3:g.100T>C\tGENE1:NM_001:exon1:c.100T>C:p.X100Y\n"
            "chr5:g.200A>G\tGENE2:NM_002:exon3:c.200A>G:p.K200R\n"
        )
        lookup = load_multianno(str(ma))
        assert len(lookup) == 2
        assert lookup["chr3:g.100T>C"].aachange == "GENE1:NM_001:exon1:c.100T>C:p.X100Y"
        assert lookup["chr3:g.100T>C"].chrom == ""

    def test_load_with_gseq(self, tmp_path: Path):
        ma = tmp_path / "tv.ma"
        ma.write_text(
            "transvar.input\tAAChange.transvar\tCHROM\tPOS\tREF\tALT\n"
            "chr3:g.100T>C\tGENE1:NM_001:exon1:c.100T>C:p.X100Y\tchr3\t100\tT\tC\n"
            "chr5:g.200A>G\tGENE2:NM_002:exon3:c.200A>G:p.K200R\tchr5\t200\tA\tG\n"
        )
        lookup = load_multianno(str(ma))
        assert len(lookup) == 2
        entry = lookup["chr3:g.100T>C"]
        assert entry.chrom == "chr3"
        assert entry.pos == "100"
        assert entry.ref == "T"
        assert entry.alt == "C"


class TestFindColumnIndex:
    def test_found(self):
        assert find_column_index(["Chr", "Start", "AAChange.refGeneWithVer"], "AAChange.refGeneWithVer") == 2

    def test_not_found(self):
        assert find_column_index(["Chr", "Start"], "Missing") is None

    def test_multiple_names(self):
        assert (
            find_column_index(
                ["Chr", "AAChange.refGene"],
                "AAChange.refGeneWithVer",
                "AAChange.refGene",
            )
            == 1
        )


class TestProcess:
    def _write_multianno(self, path: Path) -> Path:
        """Write a step2_multianno with transvar.input column."""
        path.write_text(
            "Chr\tStart\tEnd\tRef\tAlt\tFunc.refGeneWithVer\tGene.refGeneWithVer\t"
            "AAChange.refGeneWithVer\ttransvar.input\t"
            "Otherinfo1\tOtherinfo2\tOtherinfo3\tOtherinfo4\tOtherinfo5\t"
            "Otherinfo6\tOtherinfo7\tOtherinfo8\tOtherinfo9\tOtherinfo10\t"
            "Otherinfo11\tOtherinfo12\n"
            "chr3\t100\t100\tT\tC\texonic\tGENE1\t"
            "old_value\tchr3:g.100T>C\t"
            ".\t.\t.\told_chr\t100\t.\told_T\told_C\t30\tPASS\tSTATUS=ok\tGT\n"
            "chr5\t200\t200\tA\tG\texonic\tGENE2\t"
            "old_value2\tchr5:g.200A>G\t"
            ".\t.\t.\told_chr5\t200\t.\told_A\told_G\t25\tPASS\tSTATUS=ok2\tGT\n"
            "chr7\t300\t300\tC\tT\tintronic\tGENE3\t"
            "R167*\tno_match_key\t"
            ".\t.\t.\told_chr7\t300\t.\told_C\told_T\t20\tPASS\tSTATUS=ok3\tGT\n"
        )
        return path

    def _write_transvar_multianno(self, path: Path, with_gseq: bool = True) -> Path:
        if with_gseq:
            path.write_text(
                "transvar.input\tAAChange.transvar\tCHROM\tPOS\tREF\tALT\n"
                "chr3:g.100T>C\tGENE1:NM_001:exon1:c.100T>C:p.X100Y\tchr3\t100\tT\tC\n"
                "chr5:g.200A>G\tGENE2:NM_002:exon3:c.200A>G:p.K200R\tchr5\t200\tA\tG\n"
            )
        else:
            path.write_text(
                "transvar.input\tAAChange.transvar\n"
                "chr3:g.100T>C\tGENE1:NM_001:exon1:c.100T>C:p.X100Y\n"
                "chr5:g.200A>G\tGENE2:NM_002:exon3:c.200A>G:p.K200R\n"
            )
        return path

    def test_normal_replacement_with_gseq(self, tmp_path: Path):
        ma = self._write_multianno(tmp_path / "multianno.txt")
        tv = self._write_transvar_multianno(tmp_path / "tv.ma", with_gseq=True)
        out = tmp_path / "output.txt"
        process(str(ma), str(tv), str(out))
        lines = out.read_text().strip().splitlines()
        assert len(lines) == 4  # header + 3 data

        # Check AAChange replaced
        assert "GENE1:NM_001:exon1" in lines[1]
        assert "GENE2:NM_002:exon3" in lines[2]
        # Third row unchanged (no match)
        assert "R167*" in lines[3]

        # Check Otherinfo replaced by value presence in the row
        parts1 = lines[1].split("\t")
        parts2 = lines[2].split("\t")
        # Row 1: should have chr3, 100, T, C from transvar
        assert "chr3" in parts1
        assert "100" in parts1
        # Row 2: should have chr5, 200, A, G from transvar
        assert "chr5" in parts2
        assert "200" in parts2
        # Old values should NOT be present in replaced rows
        assert "old_chr" not in parts1
        assert "old_T" not in parts1
        # Third row (no match) should still have old values
        parts3 = lines[3].split("\t")
        assert "old_chr7" in parts3

    def test_without_gseq(self, tmp_path: Path):
        ma = self._write_multianno(tmp_path / "multianno.txt")
        tv = self._write_transvar_multianno(tmp_path / "tv.ma", with_gseq=False)
        out = tmp_path / "output.txt"
        process(str(ma), str(tv), str(out))
        lines = out.read_text().strip().splitlines()
        # AAChange replaced but Otherinfo unchanged
        assert "GENE1:NM_001:exon1" in lines[1]
        header = lines[0].split("\t")
        oi4_idx = header.index("Otherinfo4")
        parts1 = lines[1].split("\t")
        # Without gseq, Otherinfo should NOT be replaced
        assert parts1[oi4_idx] == "old_chr"

    def test_file_not_found_multianno(self):
        with pytest.raises(SystemExit):
            process("nonexistent.txt", "also_nonexistent.txt")

    def test_file_not_found_transvar(self, tmp_path: Path):
        ma = self._write_multianno(tmp_path / "multianno.txt")
        with pytest.raises(SystemExit):
            process(str(ma), "nonexistent.txt")
