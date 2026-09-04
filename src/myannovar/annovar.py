"""AnnovarRunner - wraps ANNOVAR table_annovar.pl for Python."""

import logging
import os
import shutil
import subprocess
from pathlib import Path
from typing import Dict, List, Optional

logger = logging.getLogger(__name__)


class AnnovarRunner:
    """Execute ANNOVAR table_annovar.pl and return output file paths."""

    SCRIPT_NAME = "table_annovar.pl"

    def __init__(self, annovar_dir: Optional[str] = None):
        self._annovar_dir = self._resolve_dir(annovar_dir)

    # -- directory resolution ------------------------------------------

    def _resolve_dir(self, user_dir: Optional[str]) -> Path:
        """Find ANNOVAR script directory."""
        if user_dir is not None:
            p = Path(user_dir)
            if not (p / self.SCRIPT_NAME).is_file():
                raise FileNotFoundError(
                    f"table_annovar.pl not found in {p}"
                )
            return p.resolve()

        # 1. Environment variable
        env_dir = os.environ.get("ANNOVAR_DIR")
        if env_dir:
            p = Path(env_dir)
            if (p / self.SCRIPT_NAME).is_file():
                return p.resolve()

        # 2. Built-in annovar_scripts/ (installed with package)
        builtin = Path(__file__).parent / "annovar_scripts"
        if (builtin / self.SCRIPT_NAME).is_file():
            return builtin.resolve()

        raise FileNotFoundError(
            "ANNOVAR scripts not found. Set ANNOVAR_DIR or provide annovar_dir."
        )

    @property
    def annovar_dir(self) -> Path:
        return self._annovar_dir

    # -- run -----------------------------------------------------------

    def run(
        self,
        input_file: str,
        output_prefix: str,
        build: str,
        humandb: str,
        protocol: str,
        operation: str,
        argument: str,
        other_args: Optional[List[str]] = None,
    ) -> Dict[str, str]:
        """Run table_annovar.pl and return paths to output files.

        Args:
            input_file: input VCF file path
            output_prefix: output file prefix (e.g. "sample.var")
            build: genome build version (hg19/hg38)
            humandb: human database directory path
            protocol: comma-separated protocol string
                      (e.g. "refGeneWithVer,cytoBand,clinvar_20220320")
            operation: comma-separated operation string
                       (e.g. "g,r,f,f,f")
            argument: comma-separated argument string
                      (e.g. "'-hgvs',,,")
            other_args: additional arguments to append

        Returns:
            {"avinput": str, "multianno": str}
        """
        logger.debug("ANNOVAR dir: %s", self._annovar_dir)
        if not shutil.which("perl"):
            raise FileNotFoundError("perl not found in PATH")

        script = self._annovar_dir / self.SCRIPT_NAME
        if not script.is_file():
            raise FileNotFoundError(f"{script} not found")

        cmd = [
            "perl",
            str(script),
            str(input_file),
            str(humandb),
            "-buildver", build,
            "-out", output_prefix,
            "-remove",
            "-protocol", protocol,
            "-operation", operation,
            "-nastring", ".",
            "-vcfinput",
            "--polish",
            "--argument", argument,
        ]

        if other_args:
            cmd.extend(other_args)

        logger.debug("Running: %s", " ".join(cmd))

        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            logger.error("ANNOVAR failed: %s", result.stderr.strip())
            raise RuntimeError(
                f"table_annovar.pl failed (rc={result.returncode}):\n{result.stderr}"
            )

        workdir = Path(input_file).resolve().parent
        avinput = workdir / f"{output_prefix}.avinput"
        multianno = workdir / f"{output_prefix}_multianno.txt"

        logger.info("ANNOVAR complete: avinput=%s, multianno=%s", avinput, multianno)

        return {
            "avinput": str(avinput),
            "multianno": str(multianno),
        }
