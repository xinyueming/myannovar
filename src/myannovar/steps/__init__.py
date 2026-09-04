"""Step modules for TransVar annotation pipeline.

Each module exports a ``process()`` function::

    from myannovar.steps.step1 import process
    process("input.avinput", "output.avinput")
"""

from myannovar.steps.step1 import process as step1_process
from myannovar.steps.step2 import process as step2_process
from myannovar.steps.step3 import process as step3_process
from myannovar.steps.step4 import process as step4_process
from myannovar.steps.step5 import process as step5_process
from myannovar.steps.step6 import process as step6_process
from myannovar.steps.step7 import process as step7_process

__all__ = [
    "step1_process",
    "step2_process",
    "step3_process",
    "step4_process",
    "step5_process",
    "step6_process",
    "step7_process",
]
