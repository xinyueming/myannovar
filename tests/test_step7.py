"""Tests for step7_to_vcf.py."""

import pytest
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from step7_to_vcf import build_info_fields, process


DATA = Path(__file__).resolve().parent.parent / "test_data"


class TestBuildInfoFields:
    def test_basic(self):
        fields = ["chr1", "100", "100", "G", "A", "exonic", "TP53"]
        header = ["Chr", "Start", "End", "Ref", "Alt", "Func.refGene", "Gene.refGene"]
        info = build_info_fields(fields, header)
        assert "Func_refGene=exonic" in info
        assert "Gene_refGene=TP53" in info

    def test_escaping(self):
        fields = ["chr1", "100", "100", "G", "A", "a;b"]
        header = ["Chr", "Start", "End", "Ref", "Alt", "Test"]
        info = build_info_fields(fields, header)
        assert "%3B" in info

    def test_skips_position_columns(self):
        fields = ["chr1", "100", "100", "G", "A"]
        header = ["Chr", "Start", "End", "Ref", "Alt"]
        info = build_info_fields(fields, header)
        assert info == "."


class TestProcess:
    def test_vcf_header(self, tmp_path: Path):
        out = tmp_path / "output.vcf"
        process(str(DATA / "step7_multianno.txt"), str(out))
        lines = out.read_text().strip().splitlines()
        assert lines[0] == "##fileformat=VCFv4.2"
        assert lines[2].startswith("#CHROM")

    def test_snv(self, tmp_path: Path):
        out = tmp_path / "output.vcf"
        process(str(DATA / "step7_multianno.txt"), str(out))
        lines = out.read_text().strip().splitlines()
        # SNV row (first data line after headers)
        assert "chr3" in lines[3]
        assert "\tG\tA\t" in lines[3]

    def test_insertion(self, tmp_path: Path):
        out = tmp_path / "output.vcf"
        process(str(DATA / "step7_multianno.txt"), str(out))
        lines = out.read_text().strip().splitlines()
        # Insertion row: ref=-, alt=A
        for line in lines[3:]:
            if line.startswith("chr1"):
                parts = line.split("\t")
                assert parts[3] == "."  # ref is empty → .
                assert parts[4] == "A"  # alt is A
                break

    def test_deletion(self, tmp_path: Path):
        out = tmp_path / "output.vcf"
        process(str(DATA / "step7_multianno.txt"), str(out))
        lines = out.read_text().strip().splitlines()
        # Deletion row: ref=AGT, alt=-
        for line in lines[3:]:
            if line.startswith("chr2"):
                parts = line.split("\t")
                assert parts[3] == "AGT"  # ref is AGT
                assert parts[4] == "."  # alt is empty → .
                break

    def test_file_not_found(self):
        with pytest.raises(SystemExit):
            process("nonexistent.txt")
