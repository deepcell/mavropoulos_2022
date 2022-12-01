"""Module containing common logging functionality."""

from __future__ import annotations

import argparse
import logging
import os

from . import Path, PathLike, parse_path


def init_logging(
    *,
    level: str | int = logging.INFO,
    log_filename: PathLike | None = None,
    log_filemode: str = "w",
    testing: bool = False,
) -> None:
    """Initialize logging.

    Call this once and only once per script.

    Args:
        level: logging level
        log_filename: log output file
        log_filemode: log file mode (default=w)
        testing: if testing is True then do not add timestamps to log entries

    """
    # Initialize logging
    log_format = "%(levelname)s\t%(name)s\t%(message)s"
    if not testing:
        log_format = "%(asctime)s\t" + log_format

    if log_filename:
        parse_path(log_filename).parent.mkdir(exist_ok=True, parents=True)

    logging.basicConfig(
        level=level,
        format=log_format,
        datefmt="%y-%m-%d %H:%M:%S",
        filename=log_filename,
        filemode=log_filemode,
    )


def parse_loglevel(value: str) -> str | int:
    """Parse a string value as a log level.

    Valid log levels are integers or a string: "CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG".  Comparison is done in
    a case-insensitive manner.

    Args:
        value: string input

    Returns:
        valid string or integer logging level

    Examples:
        >>> parse_loglevel('INFO')
        'INFO'
        >>> parse_loglevel('debug')
        'DEBUG'
        >>> parse_loglevel('10')
        10
        >>> parse_loglevel('foo')
        Traceback (most recent call last):
        ...
        ValueError: Invalid log level: foo

    """
    if value.upper() in ["CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"]:
        return value.upper()

    try:
        return int(value)
    except ValueError:
        pass

    raise ValueError(f"Invalid log level: {value}")


def add_logging_args(parser: argparse.ArgumentParser | argparse._ArgumentGroup) -> None:
    """Add arguments to configure logging.

    Args:
        parser: argument parser

    """
    parser.add_argument(
        "--testing",
        action="store_true",
        help="Set when testing to remove timestamps from logging output",
    )

    # Specify file in which to store logging output
    parser.add_argument("--log-file", type=Path, metavar="FILE", help="Write logging messages to FILE")

    # Specify logging verbosity
    parser.add_argument(
        "--log-level",
        metavar="LEVEL",
        type=parse_loglevel,
        default=os.environ.get("LOGLEVEL", "INFO"),
        help="Set logging verbosity to an integer or one of these values: DEBUG, INFO, WARNING, ERROR, CRITICAL."
        " Default=INFO unless overridden by setting the LOGLEVEL environment variable.)",
    )
