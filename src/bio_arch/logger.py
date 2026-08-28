"""Structured logging setup for bio_arch.

Provides consistent, amateur-friendly console logs and optional detailed file logging.
"""

from __future__ import annotations

import logging
from pathlib import Path
import sys


def setup_logger(
    name: str = "bio_arch",
    level: int = logging.INFO,
    log_file: Path | str | None = None,
) -> logging.Logger:
    """Configure and return a structured logger.

    Args:
        name: Logger name.
        level: Logging level (default INFO).
        log_file: Optional file path to persist detailed execution logs.

    Returns:
        Configured logging.Logger instance.
    """
    logger = logging.getLogger(name)
    logger.setLevel(level)

    # Avoid duplicate handlers if setup is called multiple times
    if logger.hasHandlers():
        logger.handlers.clear()

    # Formatter for console: clean and concise
    console_format = logging.Formatter(
        fmt="[%(levelname)s] %(asctime)s - %(name)s: %(message)s",
        datefmt="%H:%M:%S",
    )
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(level)
    console_handler.setFormatter(console_format)
    logger.addHandler(console_handler)

    if log_file:
        file_path = Path(log_file)
        file_path.parent.mkdir(parents=True, exist_ok=True)
        file_format = logging.Formatter(
            fmt="%(asctime)s [%(levelname)s] %(name)s (%(filename)s:%(lineno)d): %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )
        file_handler = logging.FileHandler(file_path, encoding="utf-8")
        file_handler.setLevel(logging.DEBUG)  # Capture all debug info in file
        file_handler.setFormatter(file_format)
        logger.addHandler(file_handler)

    return logger
