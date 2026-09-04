"""Logging configuration for myannovar."""

import logging


def setup_logging(verbose: bool = False) -> None:
    """Configure logging for the application.

    Args:
        verbose: if True, set level to DEBUG; otherwise INFO.
    """
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S",
    )
