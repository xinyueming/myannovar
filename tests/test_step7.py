"""Tests for step7 (backConvertVCF port)."""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from myannovar.steps.step7 import (
    _anno_escape,
    _build_anno_string,
    _find_otherinfo_indices,
    _is_otherinfo,
    process,
)


class TestAnnoEscape:
    def test_spaces_to_underscores(self):
        assert _anno_escape("nonsynonymous SNV") == "nonsynonymous_SNV"

    def test_semicolon_escape(self):
        assert _anno_escape("a;b") == "a\\x3bb"

    def test_equals_escape(self):
        assert _anno_escape("a=b") == "a\\x3db"

    def test_combined_escape(self):
        result = _anno_escape("status=ok;val=1")
        assert "\\x3d" in result
        assert "\\x3b" in result

    def test_dot_passthrough(self):
        assert _anno_escape(".") == "."

    def test_empty_passthrough(self):
        assert _anno_escape("") == "."


class TestHelpers:
    def test_is_otherinfo(self):
        assert _is_otherinfo("Otherinfo1") is True
        assert _is_otherinfo("Otherinfo14") is True
        assert _is_otherinfo("Func.refGene") is False

    def test_find_otherinfo_indices(self):
        cols = ["Chr", "Otherinfo1", "Func", "Otherinfo2"]
        assert _find_otherinfo_indices(cols) == [1, 3]


class TestBuildAnnoString:
    def test_basic(self):
        fields = ["chr1", "100", "100", "G", "A", "exonic", "TP53", ".", ".", "v1"]
        names = ["Chr", "Start", "End", "Ref", "Alt", "Func", "Gene", "Oi1", "Oi2", "Oi3"]
        # annot_end_idx = len(names) - 1 = 9, so range(5, 9) = cols 5,6,7,8
        result = _build_anno_string(fields, names, annot_end_idx=9, date="2024-01-01")
        assert "ANNOVAR_DATE=2024-01-01" in result
        assert "Func=exonic" in result
        assert "Gene=TP53" in result
        assert "ALLELE_END" in result

    def test_escapes_semicolons(self):
        fields = ["chr1", "100", "100", "G", "A", "a;b", "."]
        names = ["Chr", "Start", "End", "Ref", "Alt", "Test", "Oi1"]
        # annot_end_idx = len(names)-1 = 6, range(5,6) = col 5 (Test)
        result = _build_anno_string(fields, names, annot_end_idx=6, date="2024-01-01")
        assert "\\x3b" in result
        assert "Test=a\\x3bb" in result


class TestProcess:
    def _write_test_vcf(self, path: Path) -> Path:
        """Write a minimal VCF for use as original_vcf input."""
        path.write_text(
            "##fileformat=VCFv4.2\n"
            "##INFO=<ID=STATUS,Number=1,Type=String,Description=\"Status\">\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\n"
        )
        return path

    def _write_test_multianno(self, path: Path) -> Path:
        """Write a minimal multianno with Otherinfo columns."""
        cols = [
            "Chr", "Start", "End", "Ref", "Alt",
            "Func.refGene", "Gene.refGene",
            "Otherinfo1", "Otherinfo2", "Otherinfo3",
            "Otherinfo4", "Otherinfo5", "Otherinfo6",
            "Otherinfo7", "Otherinfo8", "Otherinfo9",
            "Otherinfo10", "Otherinfo11", "Otherinfo12",
            "Otherinfo13",  # SAMPLE column (needed for unique merge key)
        ]
        header = "\t".join(cols)
        # One data row: chr3 100 G>A with sample data
        data = "\t".join([
            "chr3", "100", "100", "G", "A",
            "exonic", "TP53",
            ".", ".", ".",  # Otherinfo1-3
            "chr3", "100", ".",  # Otherinfo4-6 (CHROM, POS, ID)
            "G", "A", "30",  # Otherinfo7-9 (REF, ALT, QUAL)
            "PASS", "STATUS=test",  # Otherinfo10-11 (FILTER, INFO)
            "GT", "0/1",  # Otherinfo12-13 (FORMAT, SAMPLE)
        ])
        path.write_text(header + "\n" + data + "\n")
        return path

    def test_vcf_header(self, tmp_path: Path):
        vcf = self._write_test_vcf(tmp_path / "orig.vcf")
        ma = self._write_test_multianno(tmp_path / "ma.txt")
        out = tmp_path / "output.vcf"
        process(str(ma), str(vcf), str(out))
        lines = out.read_text().strip().splitlines()
        assert lines[0] == "##fileformat=VCFv4.2"
        # ANNOVAR_DATE header present
        assert any("ANNOVAR_DATE" in l for l in lines)
        # ALLELE_END header present
        assert any("ALLELE_END" in l for l in lines)
        # Func.refGene header present
        assert any("Func.refGene" in l for l in lines)
        # CHROM line
        assert any(l.startswith("#CHROM") for l in lines)

    def test_data_row(self, tmp_path: Path):
        vcf = self._write_test_vcf(tmp_path / "orig.vcf")
        ma = self._write_test_multianno(tmp_path / "ma.txt")
        out = tmp_path / "output.vcf"
        process(str(ma), str(vcf), str(out))
        lines = out.read_text().strip().splitlines()
        # Find first data line
        data_lines = [l for l in lines if not l.startswith("#")]
        assert len(data_lines) == 1
        parts = data_lines[0].split("\t")
        assert parts[0] == "chr3"
        assert parts[1] == "100"
        assert parts[3] == "G"
        assert parts[4] == "A"
        assert "ANNOVAR_DATE" in parts[7]
        assert "Func.refGene=exonic" in parts[7]
        assert "Gene.refGene=TP53" in parts[7]
        assert "ALLELE_END" in parts[7]

    def test_insertion_ref_dash(self, tmp_path: Path):
        vcf = self._write_test_vcf(tmp_path / "orig.vcf")
        ma = self._write_test_multianno(tmp_path / "ma.txt")
        # Add an insertion row: multianno Ref=-, Otherinfo7=.(no original VCF REF)
        # In real VCF, insertions have a reference base (e.g. "T" for T>TA).
        # When multianno has Ref=-, Otherinfo7 may be "." if original VCF had no REF.
        with open(ma, "a") as f:
            f.write("\t".join([
                "chr1", "200", "200", "-", "A",
                "intronic", "BRCA1",
                ".", ".", ".",
                "chr1", "200", ".",
                ".", "A", "25",
                "PASS", "TYPE=INS",
                "GT", "0/0",
            ]) + "\n")
        out = tmp_path / "output.vcf"
        process(str(ma), str(vcf), str(out))
        lines = out.read_text().strip().splitlines()
        data_lines = [l for l in lines if not l.startswith("#")]
        ins = [l for l in data_lines if "intronic" in l and "BRCA1" in l]
        assert len(ins) == 1
        parts = ins[0].split("\t")
        assert parts[3] == "."  # Otherinfo7=. becomes .
        assert parts[4] == "A"

    def test_deletion_alt_dash(self, tmp_path: Path):
        vcf = self._write_test_vcf(tmp_path / "orig.vcf")
        ma = self._write_test_multianno(tmp_path / "ma.txt")
        # Deletion: multianno Alt=-, Otherinfo8=.(no original VCF ALT)
        with open(ma, "a") as f:
            f.write("\t".join([
                "chr2", "300", "300", "AGT", "-",
                "exonic", "EGFR",
                ".", ".", ".",
                "chr2", "300", ".",
                "AGT", ".", "20",
                "PASS", "TYPE=DEL",
                "GT", "1/0",
            ]) + "\n")
        out = tmp_path / "output.vcf"
        process(str(ma), str(vcf), str(out))
        lines = out.read_text().strip().splitlines()
        data_lines = [l for l in lines if not l.startswith("#")]
        dels = [l for l in data_lines if "exonic" in l and "EGFR" in l]
        assert len(dels) == 1
        parts = dels[0].split("\t")
        assert parts[3] == "AGT"
        assert parts[4] == "."  # Otherinfo8=. becomes .

    def test_file_not_found(self):
        with pytest.raises(SystemExit):
            process("nonexistent.txt", "also_nonexistent.vcf")

    def test_uses_original_vcf_header(self, tmp_path: Path):
        vcf = tmp_path / "orig.vcf"
        vcf.write_text(
            "##fileformat=VCFv4.2\n"
            "##custom=header\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        )
        ma = self._write_test_multianno(tmp_path / "ma.txt")
        out = tmp_path / "output.vcf"
        process(str(ma), str(vcf), str(out))
        lines = out.read_text().strip().splitlines()
        assert lines[0] == "##fileformat=VCFv4.2"
        assert lines[1] == "##custom=header"
        # ANNOVAR headers inserted before #CHROM
        chrom_idx = next(i for i, l in enumerate(lines) if l.startswith("#CHROM"))
        assert chrom_idx > 2  # At least ANNOVAR_DATE + ALLELE_END inserted
