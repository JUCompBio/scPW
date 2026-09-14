"""
Structured logging: console and file, one line per event.

The problem this solves. The parallel grid ran workers as subprocesses that returned a
whole chunk of results at once, so nothing appeared until six trainings had finished
and the CSV stayed empty for tens of minutes. Progress was invisible and an interrupted
run lost everything in flight.

So every worker logs each training as it completes, to stdout and to a JSONL file. The
console output is human-readable; the JSONL is machine-readable so a run can be
analysed or resumed without parsing prose.

Log files go to logs/<name>_<timestamp>.jsonl. They are gitignored: they are a record
of a run, not an input to one.
"""

from __future__ import annotations

import logging
import os
import sys
import time
from pathlib import Path

import structlog

ROOT = Path(__file__).resolve().parents[1]
LOG_DIR = ROOT / "logs"


def setup(name: str = "run", to_file: bool = True) -> structlog.BoundLogger:
    """Configure structlog for console plus optional JSONL file output.

    Args:
        name: stem for the log file.
        to_file: write JSONL alongside the console output. Off for short-lived
            subprocesses that would otherwise each create their own file.

    Returns:
        A bound logger. Callers should pass structured keys rather than formatted
        strings, so the JSONL stays queryable.
    """
    LOG_DIR.mkdir(exist_ok=True)

    processors = [
        structlog.contextvars.merge_contextvars,
        structlog.processors.add_log_level,
        structlog.processors.TimeStamper(fmt="%H:%M:%S"),
    ]

    if to_file:
        # A second sink via the stdlib, so the JSONL gets every event the console
        # does without the caller having to log twice.
        stamp = time.strftime("%Y%m%d_%H%M%S")
        path = LOG_DIR / f"{name}_{stamp}.jsonl"
        handler = logging.FileHandler(path, encoding="utf-8")
        handler.setFormatter(logging.Formatter("%(message)s"))
        root = logging.getLogger("scpw")
        root.setLevel(logging.INFO)
        root.handlers = [handler]
        root.propagate = False
        os.environ["SCPW_LOG_FILE"] = str(path)

    structlog.configure(
        processors=processors + [
            structlog.dev.ConsoleRenderer(colors=sys.stdout.isatty()),
        ],
        wrapper_class=structlog.make_filtering_bound_logger(logging.INFO),
        cache_logger_on_first_use=True,
    )
    return structlog.get_logger()


def log_jsonl(**fields) -> None:
    """Append one JSON record to the run's log file, if one is configured.

    Kept separate from the console logger because subprocess workers write their
    progress to the parent's file by path rather than by inheriting a handler.
    """
    path = os.environ.get("SCPW_LOG_FILE")
    if not path:
        return
    import json

    fields.setdefault("ts", time.strftime("%H:%M:%S"))
    with open(path, "a", encoding="utf-8") as fh:
        fh.write(json.dumps(fields) + "\n")
