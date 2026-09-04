"""Pipeline - orchestrates ANNOVAR + TransVar annotation."""

import logging
from pathlib import Path
from typing import List, Optional

from myannovar.annovar import AnnovarRunner
from myannovar.steps.step1 import process as step1_process
from myannovar.steps.step2 import process as step2_process
from myannovar.steps.step3 import process as step3_process
from myannovar.steps.step4 import process as step4_process
from myannovar.steps.step5 import process as step5_process
from myannovar.steps.step6 import process as step6_process
from myannovar.steps.step7 import process as step7_process

logger = logging.getLogger(__name__)

_STEP_NAMES = {
    1: "Add transvar.input column",
    2: "Merge transvar.input into multianno",
    3: "Extract unique transvar.input values",
    4: "Run TransVar ganno",
    5: "Parse TransVar output",
    6: "Replace AAChange with TransVar annotations",
    7: "Convert to VCF",
}


class Pipeline:
    """Orchestrate ANNOVAR + TransVar annotation pipeline."""

    def __init__(
        self,
        workdir: Optional[Path] = None,
        annovar_dir: Optional[str] = None,
        keep_temp: bool = False,
    ):
        self.workdir = workdir or Path.cwd()
        self.keep_temp = keep_temp
        self._annovar_dir = annovar_dir
        self._temp_files: List[Path] = []

    # -- high-level entry points ---------------------------------------

    def run_annovar(
        self,
        input_vcf: str,
        output_prefix: str,
        build: str,
        humandb: str,
        protocol: str,
        operation: str,
        argument: str,
    ) -> dict:
        """Run ANNOVAR annotation.

        Returns dict with ``avinput`` and ``multianno`` paths.
        """
        logger.info("Running ANNOVAR: %s (build=%s, protocol=%s)", input_vcf, build, protocol)
        runner = AnnovarRunner(annovar_dir=self._annovar_dir)
        result = runner.run(
            input_file=input_vcf,
            output_prefix=output_prefix,
            build=build,
            humandb=humandb,
            protocol=protocol,
            operation=operation,
            argument=argument,
        )
        logger.info("ANNOVAR output: avinput=%s, multianno=%s", result["avinput"], result["multianno"])
        return result

    def run_transvar(
        self,
        avinput: str,
        multianno: str,
        output: str,
        refseq: bool = True,
    ) -> str:
        """Run TransVar 7-step annotation pipeline.

        Args:
            avinput: .avinput file (from ANNOVAR)
            multianno: _multianno.txt file (from ANNOVAR)
            output: output VCF path
            refseq: whether to use --refseq for TransVar

        Returns:
            Path to the output VCF file.
        """
        logger.info("Starting TransVar pipeline (%d steps)", len(_STEP_NAMES))
        workdir = Path(avinput).resolve().parent

        for step_num, (step_fn, step_args) in enumerate(
            [
                (step1_process, (avinput, str(workdir / "step1.avinput"))),
                (step2_process, (str(workdir / "step1.avinput"), multianno, str(workdir / "step2_multianno.txt"))),
                (step3_process, (str(workdir / "step2_multianno.txt"), str(workdir / "transvar.input"))),
                (step4_process, (str(workdir / "transvar.input"), str(workdir / "transvar.output"), refseq)),
                (step5_process, (str(workdir / "transvar.output"), str(workdir / "transvar.multianno"))),
                (step6_process, (str(workdir / "step2_multianno.txt"), str(workdir / "transvar.multianno"), str(workdir / "step6_multianno.txt"))),
                (step7_process, (str(workdir / "step6_multianno.txt"), output)),
            ],
            start=1,
        ):
            name = _STEP_NAMES[step_num]
            logger.info("Step %d: %s", step_num, name)
            out_path = step_args[-1] if isinstance(step_args[-1], str) else str(step_args[-1])
            self._track(out_path)
            step_fn(*step_args)

        if not self.keep_temp:
            logger.debug("Cleaning up %d intermediate files", len(self._temp_files))
            self._cleanup()

        logger.info("TransVar pipeline complete: %s", output)
        return output

    def run_full(
        self,
        input_vcf: str,
        output_vcf: str,
        build: str,
        humandb: str,
        protocol: str,
        operation: str,
        argument: str,
        refseq: bool = True,
    ) -> str:
        """Run complete ANNOVAR + TransVar pipeline.

        Args:
            input_vcf: input VCF file
            output_vcf: output VCF file
            build: genome build
            humandb: human database directory
            protocol: comma-separated protocol string
            operation: comma-separated operation string
            argument: comma-separated argument string
            refseq: whether to use RefSeq

        Returns:
            Path to the output VCF file.
        """
        logger.info("Starting full pipeline: %s -> %s", input_vcf, output_vcf)
        output_prefix = Path(output_vcf).stem
        avinput_path = str(Path(input_vcf).parent / f"{output_prefix}.avinput")
        multianno_path = str(Path(input_vcf).parent / f"{output_prefix}_multianno.txt")

        # Phase 1: ANNOVAR
        self.run_annovar(
            input_vcf, output_prefix, build, humandb,
            protocol, operation, argument,
        )

        # Phase 2: TransVar
        self.run_transvar(avinput_path, multianno_path, output_vcf, refseq)

        if not self.keep_temp:
            self._cleanup()

        logger.info("Full pipeline complete: %s", output_vcf)
        return output_vcf

    # -- internal ------------------------------------------------------

    def _track(self, path: str) -> None:
        self._temp_files.append(Path(path))

    def _cleanup(self) -> None:
        """Remove intermediate files."""
        for f in self._temp_files:
            if f.is_file():
                try:
                    f.unlink()
                except OSError:
                    pass
        self._temp_files.clear()
