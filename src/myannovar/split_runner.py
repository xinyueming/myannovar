"""SplitRunner - split VCF by chromosome, annotate in parallel, merge results."""

from __future__ import annotations

import logging
import shutil
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)

# Human chromosome sort order: 1-22, X, Y, M/MT, then others alphabetically
_CHROM_ORDER = {str(i): i for i in range(1, 23)}
_CHROM_ORDER.update({"X": 23, "Y": 24, "M": 25, "MT": 25})


def _sort_key(chrom: str) -> tuple:
    """Sort key for chromosome names."""
    c = chrom.replace("chr", "").replace("Chr", "").replace("CHR", "")
    return (_CHROM_ORDER.get(c, 99), c)


def _check_bcftools() -> str:
    """Return bcftools path or raise."""
    path = shutil.which("bcftools")
    if path is None:
        raise FileNotFoundError(
            "bcftools not found in PATH. Install bcftools >= 1.10."
        )
    return path


def _run(cmd: list[str], *, check: bool = True) -> subprocess.CompletedProcess:
    logger.debug("Running: %s", " ".join(cmd))
    return subprocess.run(cmd, capture_output=True, text=True, check=check)


def get_chromosomes(vcf_path: str) -> list[str]:
    """Extract unique chromosome names from VCF data rows.

    Skips header lines (## and #CHROM), reads column 1, returns
    sorted list using biological chromosome order.
    """
    bcftools = _check_bcftools()
    result = _run(
        [bcftools, "query", "-l", vcf_path],
        check=False,
    )
    # bcftools query -l lists samples, not chromosomes.
    # Use bcftools view -H to get data rows only (no header).
    result = _run(
        [bcftools, "view", "-H", vcf_path],
    )
    chroms = set()
    for line in result.stdout.splitlines():
        parts = line.split("\t", 1)
        if parts:
            chroms.add(parts[0])
    if not chroms:
        raise ValueError(f"No data rows found in {vcf_path}")
    return sorted(chroms, key=_sort_key)


def _ensure_indexed(vcf_path: str) -> str:
    """Ensure VCF is tabix-indexed. Returns (possibly compressed and indexed) VCF path."""
    bcftools = _check_bcftools()
    path = Path(vcf_path)

    # Check for existing index (.csi or .tbi)
    for ext in (".csi", ".tbi"):
        if path.with_suffix(path.suffix + ext).exists() or path.with_name(path.name + ext).exists():
            return vcf_path

    # If uncompressed VCF, compress first
    if not vcf_path.endswith(".gz"):
        compressed = path.with_suffix(path.suffix + ".gz")
        _run([bcftools, "view", vcf_path, "-Oz", "-o", str(compressed)])
        vcf_path = str(compressed)

    # Create tabix index
    _run([bcftools, "index", "-t", vcf_path])
    return vcf_path


def split_vcf(
    vcf_path: str,
    tmp_dir: Path,
    chroms: list[str],
) -> dict[str, Path]:
    """Split VCF by chromosome using bcftools view.

    Creates {tmp_dir}/{chrom}/input.vcf.gz for each chromosome.
    Returns {chrom: vcf_path} mapping.
    """
    bcftools = _check_bcftools()

    # Ensure VCF is indexed (required for -r region queries)
    vcf_path = _ensure_indexed(vcf_path)

    result = {}
    for chrom in chroms:
        chrom_dir = tmp_dir / chrom
        chrom_dir.mkdir(parents=True, exist_ok=True)
        out_vcf = chrom_dir / "input.vcf.gz"
        _run(
            [bcftools, "view", "-r", chrom, "-Oz", "-o", str(out_vcf), vcf_path],
        )
        # Skip chromosomes with no variants (empty output file)
        if out_vcf.stat().st_size == 0:
            logger.warning("Chromosome %s has no variants, skipping", chrom)
            out_vcf.unlink()
            continue
        result[chrom] = out_vcf
        logger.info("Split %s -> %s", chrom, out_vcf)
    return result


def merge_vcfs(
    input_vcf: str,
    output_vcf: str,
    tmp_dir: Path,
    chroms: list[str],
) -> None:
    """Merge annotated VCFs by chromosome order using bcftools concat.

    The annotated output files are named {output_prefix}.{build}_multianno.vcf
    (or .vcf.gz if the original was gzipped).
    """
    bcftools = _check_bcftools()
    annotated_files: list[Path] = []
    for chrom in chroms:
        chrom_dir = tmp_dir / chrom
        # Find the annotated VCF - it should be the output_vcf basename
        # Since each chromosome writes to the same output_vcf path in its own dir,
        # we need to use a chromosome-specific output name.
        # The convention is: {chrom}_annotated.vcf in each chrom_dir.
        annotated = chrom_dir / f"{chrom}_annotated.vcf"
        if not annotated.is_file():
            raise FileNotFoundError(f"Annotated VCF not found: {annotated}")
        annotated_files.append(annotated)

    out_path = Path(output_vcf)
    if str(output_vcf).endswith(".gz"):
        _run(
            [
                bcftools,
                "concat",
                "-Oz",
                "-o",
                str(out_path),
                *[str(f) for f in annotated_files],
            ],
        )
        _run([bcftools, "index", str(out_path)], check=False)
    else:
        _run(
            [
                bcftools,
                "concat",
                "-Ov",
                "-o",
                str(out_path),
                *[str(f) for f in annotated_files],
            ],
        )
    logger.info("Merged %d chromosomes -> %s", len(chroms), output_vcf)


def merge_multianno(tmp_dir: Path, chroms: list[str], output_path: str) -> None:
    """Merge step6_multianno.txt from each chromosome directory.

    Takes the header from the first file only, then concatenates
    data rows from all chromosomes in order.
    """
    header_written = False
    with open(output_path, "w") as out:
        for chrom in chroms:
            chrom_step6 = tmp_dir / chrom / f"{chrom}_annotated.txt"
            if not chrom_step6.is_file():
                logger.warning("step6_multianno.txt not found for %s, skipping", chrom)
                continue
            lines = chrom_step6.read_text().splitlines()
            if not lines:
                continue
            if not header_written:
                out.write(lines[0] + "\n")
                header_written = True
            for line in lines[1:]:
                if line.strip():
                    out.write(line + "\n")
    logger.info("Merged step6_multianno.txt -> %s", output_path)


class SplitRunner:
    """Split VCF by chromosome, annotate in parallel, merge results."""

    def __init__(self, jobs: int = 1, keep_temp: bool = False):
        self.jobs = max(1, jobs)
        self.keep_temp = keep_temp
        self._tmp_dir: Optional[Path] = None

    def _run_single_chrom(
        self,
        chrom: str,
        chrom_vcf: Path,
        output_vcf: str,
        workdir: Path,
        annovar_dir: Optional[str],
        build: str,
        humandb: str,
        protocol: str,
        operation: str,
        argument: str,
        refseq: bool,
        polish: bool,
    ) -> str:
        """Run full annotation pipeline for a single chromosome."""
        from myannovar.pipeline import Pipeline

        chrom_dir = chrom_vcf.parent
        logger.info("Starting annotation for %s (dir=%s)", chrom, chrom_dir)

        t0 = time.time()

        # Chromosome-specific output filename to avoid collisions
        chrom_output = str(chrom_dir / f"{chrom}_annotated.vcf")

        pipeline = Pipeline(
            workdir=chrom_dir,
            annovar_dir=annovar_dir,
            keep_temp=True,  # keep intermediate files for debugging if one fails
        )

        # output_prefix for ANNOVAR: use chrom name
        output_prefix = chrom

        pipeline.run_full(
            input_vcf=str(chrom_vcf),
            output_vcf=chrom_output,
            build=build,
            humandb=humandb,
            protocol=protocol,
            operation=operation,
            argument=argument,
            refseq=refseq,
            polish=polish,
        )

        elapsed = time.time() - t0
        logger.info("Chromosome %s complete (%.1fs)", chrom, elapsed)

        return chrom_output

    def run(
        self,
        input_vcf: str,
        output_vcf: str,
        annovar_dir: Optional[str],
        build: str,
        humandb: str,
        protocol: str,
        operation: str,
        argument: str,
        refseq: bool = True,
        polish: bool = True,
    ) -> str:
        """Full split-annotate-merge pipeline.

        Args:
            input_vcf: input VCF file path
            output_vcf: final output VCF path
            annovar_dir: ANNOVAR scripts directory
            build: genome build (hg19/hg38)
            humandb: human database directory
            protocol: comma-separated protocol
            operation: comma-separated operation
            argument: comma-separated argument
            refseq: use RefSeq annotations
            polish: polish protein notation

        Returns:
            Path to the merged output VCF.
        """
        _check_bcftools()
        workdir = Path(input_vcf).resolve().parent
        t_start = time.time()

        # Step 1: Get chromosomes
        logger.info("Scanning chromosomes from %s", input_vcf)
        chroms = get_chromosomes(input_vcf)
        logger.info("Found %d chromosomes: %s", len(chroms), chroms)

        if len(chroms) == 1:
            logger.info("Only 1 chromosome, running single-threaded")
            self.jobs = 1

        # Step 2: Create temp directory
        timestamp = int(time.time())
        self._tmp_dir = workdir / f"split_tmp_{timestamp}"
        self._tmp_dir.mkdir(parents=True, exist_ok=True)
        logger.info("Temp directory: %s", self._tmp_dir)

        # Step 3: Split VCF
        logger.info("Splitting VCF into %d chromosomes", len(chroms))
        split_files = split_vcf(input_vcf, self._tmp_dir, chroms)

        # Step 4: Run in parallel
        logger.info("Running annotation with %d jobs", self.jobs)
        annotated: dict[str, str] = {}
        errors: list[tuple[str, Exception]] = []

        with ThreadPoolExecutor(max_workers=self.jobs) as executor:
            futures = {
                executor.submit(
                    self._run_single_chrom,
                    chrom=chrom,
                    chrom_vcf=vcf_path,
                    output_vcf=output_vcf,
                    workdir=workdir,
                    annovar_dir=annovar_dir,
                    build=build,
                    humandb=humandb,
                    protocol=protocol,
                    operation=operation,
                    argument=argument,
                    refseq=refseq,
                    polish=polish,
                ): chrom
                for chrom, vcf_path in split_files.items()
            }
            for future in as_completed(futures):
                chrom = futures[future]
                try:
                    result = future.result()
                    annotated[chrom] = result
                    logger.info("Chromosome %s complete", chrom)
                except Exception as e:
                    errors.append((chrom, e))
                    logger.error("Chromosome %s failed: %s", chrom, e)

        if errors:
            logger.error(
                "%d chromosome(s) failed. Temp files kept at %s for debugging.",
                len(errors),
                self._tmp_dir,
            )
            details = "; ".join(f"{c}: {e}" for c, e in errors)
            raise RuntimeError(
                f"Annotation failed for: {', '.join(c for c, _ in errors)}. "
                f"Details: {details}"
            )

        # Step 5: Merge results
        logger.info("Merging %d annotated chromosomes", len(annotated))
        merge_vcfs(input_vcf, output_vcf, self._tmp_dir, chroms)

        # Step 5b: Merge step6_multianno.txt from each chromosome
        txt_output = str(Path(output_vcf).with_suffix(".txt"))
        merge_multianno(self._tmp_dir, chroms, txt_output)

        # Step 6: Cleanup
        if not self.keep_temp and self._tmp_dir:
            logger.info("Cleaning up temp directory: %s", self._tmp_dir)
            shutil.rmtree(self._tmp_dir, ignore_errors=True)

        elapsed_total = time.time() - t_start
        logger.info("Split pipeline complete: %s (%.1fs total)", output_vcf, elapsed_total)
        return output_vcf
