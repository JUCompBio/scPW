"""
Structured experiment logging for scPathWalk.

Every pipeline invocation records exactly one append-only JSON line to
``<root>/runs.jsonl`` capturing the run identity, code version, full
hyperparameter set, metrics, and output paths. Per-run artefacts (figures,
rank CSVs) are written under ``<root>/<run_id>/`` so that repeated runs never
clobber one another.

Usage:
    from scpathwalk.logging_utils import RunLogger

    logger = RunLogger(method="scPathWalk", dataset="my_sample", params=p)
    ...
    logger.log_metrics("scPathWalk", {"silhouette": 0.62, ...})
    logger.log_artifact("umap_clusters", path)
    logger.finalize(status="ok")
"""

from __future__ import annotations

import json
import platform
import subprocess
import time
import uuid
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import structlog
from rich.console import Console
from rich.table import Table

RESULTS_ROOT = Path("results")
RUNS_LOG = RESULTS_ROOT / "runs.jsonl"

console = Console()
_LOGGING_CONFIGURED = False


def configure_logging(level: str = "info") -> structlog.stdlib.BoundLogger:
    """Configure structlog with a rich-friendly console renderer (idempotent).

    Returns a bound logger. Safe to call multiple times; only configures once.
    """
    global _LOGGING_CONFIGURED
    if not _LOGGING_CONFIGURED:
        structlog.configure(
            processors=[
                structlog.contextvars.merge_contextvars,
                structlog.processors.add_log_level,
                structlog.processors.TimeStamper(fmt="%H:%M:%S"),
                structlog.dev.ConsoleRenderer(colors=True),
            ],
            wrapper_class=structlog.make_filtering_bound_logger(
                {"debug": 10, "info": 20, "warning": 30, "error": 40}.get(level, 20)
            ),
            logger_factory=structlog.PrintLoggerFactory(),
            cache_logger_on_first_use=True,
        )
        _LOGGING_CONFIGURED = True
    return structlog.get_logger()



def _git_commit() -> str | None:
    """Return the current git commit hash, or None if unavailable."""
    try:
        out = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            capture_output=True, text=True, timeout=10,
        )
        if out.returncode == 0:
            return out.stdout.strip()
    except Exception:
        pass
    return None


def _git_dirty() -> bool | None:
    """Return True if the working tree has uncommitted changes."""
    try:
        out = subprocess.run(
            ["git", "status", "--porcelain"],
            capture_output=True, text=True, timeout=10,
        )
        if out.returncode == 0:
            return bool(out.stdout.strip())
    except Exception:
        pass
    return None


def _json_safe(obj: Any) -> Any:
    """Best-effort conversion of values into JSON-serializable forms."""
    if isinstance(obj, dict):
        return {k: _json_safe(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_json_safe(v) for v in obj]
    if isinstance(obj, (str, int, float, bool)) or obj is None:
        return obj
    # numpy scalars, Path, etc.
    try:
        import numpy as np
        if isinstance(obj, np.generic):
            return obj.item()
    except Exception:
        pass
    return str(obj)


class RunLogger:
    """Collects a single structured run record and appends it to the run log."""

    def __init__(self, method: str, dataset: str, params: dict, tag: str | None = None):
        self.run_id = datetime.now().strftime("%Y%m%d_%H%M%S") + "_" + uuid.uuid4().hex[:6]
        self.method = method
        self.dataset = dataset
        self.tag = tag
        self._start = time.time()
        self.run_dir = RESULTS_ROOT / self.run_id
        self.run_dir.mkdir(parents=True, exist_ok=True)

        configure_logging()
        self.log = structlog.get_logger().bind(
            run_id=self.run_id, method=method, dataset=dataset,
        )

        self.record: dict[str, Any] = {
            "run_id": self.run_id,
            "timestamp_utc": datetime.now(timezone.utc).isoformat(),
            "method": method,
            "dataset": dataset,
            "tag": tag,
            "git_commit": _git_commit(),
            "git_dirty": _git_dirty(),
            "platform": platform.platform(),
            "numerics": self.numerics_versions(),
            "params": _json_safe(params),
            "metrics": {},        # {method_label: {metric: value}}
            "artifacts": {},      # {name: path}
            "status": "running",
            "duration_sec": None,
            "error": None,
        }
        self.log.info("run_started", git_commit=self.record["git_commit"],
                      git_dirty=self.record["git_dirty"])

    def stage(self, name: str, **fields: Any) -> None:
        """Emit a structured log event for a pipeline stage."""
        self.log.info(name, **fields)

    def log_metrics(self, label: str, metrics: dict) -> None:
        """Record a metric dict under a label (e.g., 'scPathWalk', 'scVI_baseline')."""
        self.record["metrics"][label] = _json_safe(metrics)
        self.log.info("metrics_recorded", label=label,
                      **{k: (round(float(v), 4) if _is_number(v) else v)
                         for k, v in metrics.items()})
        print_metrics_table(metrics, title=f"{label} - {self.dataset}")


    @staticmethod
    def numerics_versions() -> dict[str, str]:
        """Versions of the libraries that determine the numbers.

        Recorded per run because these are not incidental. Upgrading torch from
        2.10 to 2.11 shifted one dataset's cluster count from 14 to 18 at a fixed seed
        and configuration, which is larger than most effects being measured. A run
        is only comparable to another run with the same stack.
        """
        import importlib.metadata as md

        out: dict[str, str] = {}
        for pkg in (
            "torch", "scvi-tools", "scanpy", "anndata", "numpy", "scipy",
            "scikit-learn", "networkx", "igraph", "leidenalg", "umap-learn",
            "pandas",
        ):
            try:
                out[pkg] = md.version(pkg)
            except md.PackageNotFoundError:
                out[pkg] = "absent"
        return out

    def log_wins(
        self,
        method: str = "scPathWalk",
        reference: str = "scVI_baseline",
    ) -> dict | None:
        """Score ``method`` against ``reference`` on the geometric metrics.

        Both metric dicts must already be recorded. Writes a ``wins`` block into
        the run record and prints it, so the win count for a run is visible in the
        console and in the ledger without post-hoc computation.

        ``n_clusters`` is reported alongside deliberately: Silhouette and CH both
        respond to partition granularity, and across the first sweep the metric
        gap correlated with the cluster-count gap (corr = -0.79 for CH). A win
        accompanied by a large cluster-count difference deserves less weight.
        """
        higher_is_better = {
            "silhouette": True,
            "calinski_harabasz": True,
            "davies_bouldin": False,
        }

        got = self.record["metrics"].get(method)
        ref = self.record["metrics"].get(reference)
        if not got or not ref:
            return None

        detail: dict[str, dict] = {}
        won = 0
        for metric, higher in higher_is_better.items():
            a, b = got.get(metric), ref.get(metric)
            if a is None or b is None:
                continue
            win = (a > b) if higher else (a < b)
            won += int(win)
            detail[metric] = {
                "value": round(float(a), 6),
                "reference": round(float(b), 6),
                "delta": round(float(a) - float(b), 6),
                "win": bool(win),
            }

        summary = {
            "method": method,
            "reference": reference,
            "wins": won,
            "of": len(detail),
            "metrics": detail,
            "n_clusters": got.get("n_clusters"),
            "n_clusters_reference": ref.get("n_clusters"),
        }
        self.record["wins"] = summary

        self.log.info(
            "wins_recorded",
            wins=f"{won}/{len(detail)}",
            won_on=",".join(k for k, v in detail.items() if v["win"]) or "none",
            delta_k=(
                (got.get("n_clusters") or 0) - (ref.get("n_clusters") or 0)
            ),
        )
        print_wins_table(summary, title=f"{method} vs {reference} - {self.dataset}")
        return summary

    def log_artifact(self, name: str, path: str | Path) -> None:
        """Record the path of a produced artefact (figure, CSV, etc.)."""
        self.record["artifacts"][name] = str(path)

    def log_extra(self, key: str, value: Any) -> None:
        """Record an arbitrary scalar/structure into the run record."""
        self.record[key] = _json_safe(value)

    def finalize(self, status: str = "ok", error: str | None = None) -> Path:
        """Write the run record as one JSON line to the ledger."""
        self.record["status"] = status
        self.record["error"] = error
        self.record["duration_sec"] = round(time.time() - self._start, 2)

        RESULTS_ROOT.mkdir(parents=True, exist_ok=True)
        with open(RUNS_LOG, "a", encoding="utf-8") as fh:
            fh.write(json.dumps(self.record) + "\n")

        # Also drop a standalone copy in the run dir for convenience.
        with open(self.run_dir / "run.json", "w", encoding="utf-8") as fh:
            json.dump(self.record, fh, indent=2)

        log_fn = self.log.info if status == "ok" else self.log.error
        log_fn("run_finalized", status=status,
               duration_sec=self.record["duration_sec"], error=error)
        return self.run_dir



def load_runs(path: str | Path = RUNS_LOG):
    """Load the run log into a pandas DataFrame (one row per run)."""
    import pandas as pd
    rows = []
    p = Path(path)
    if not p.exists():
        return pd.DataFrame()
    with open(p, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if line:
                rows.append(json.loads(line))
    return pd.DataFrame(rows)


def _is_number(v: Any) -> bool:
    try:
        float(v)
        return True
    except (TypeError, ValueError):
        return False


def print_metrics_table(metrics: dict, title: str = "Metrics") -> None:
    """Pretty-print a single metrics dict with rich."""
    table = Table(title=title, show_header=True, header_style="bold cyan")
    table.add_column("Metric")
    table.add_column("Value", justify="right")
    for k, v in metrics.items():
        table.add_row(k, f"{float(v):.4f}" if _is_number(v) else str(v))
    console.print(table)


def print_wins_table(summary: dict, title: str = "Wins") -> None:
    """Pretty-print a win summary produced by ``RunLogger.log_wins``."""
    won, of = summary["wins"], summary["of"]
    colour = "green" if won == of else ("yellow" if won else "red")
    table = Table(
        title=f"{title}   [{colour}]{won}/{of}[/{colour}]",
        show_header=True,
        header_style="bold cyan",
    )
    table.add_column("Metric")
    table.add_column(summary["method"], justify="right")
    table.add_column(summary["reference"], justify="right")
    # ASCII only: the console encoding is cp1252 when output is piped on Windows,
    # and a non-encodable character here aborted an otherwise successful run.
    table.add_column("delta", justify="right")
    table.add_column("win", justify="center")

    for metric, d in summary["metrics"].items():
        fmt = ".1f" if metric == "calinski_harabasz" else ".4f"
        table.add_row(
            metric,
            f"{d['value']:{fmt}}",
            f"{d['reference']:{fmt}}",
            f"{d['delta']:+{fmt}}",
            "[green]yes[/green]" if d["win"] else "[red]no[/red]",
        )

    k, k_ref = summary.get("n_clusters"), summary.get("n_clusters_reference")
    if k is not None and k_ref is not None:
        table.add_row("n_clusters", str(k), str(k_ref), f"{k - k_ref:+d}", "-")
    console.print(table)


def print_runs_table(path: str | Path = RUNS_LOG,
                      metric_label: str = "scPathWalk",
                      metric_keys: list[str] | None = None) -> None:
    """Pretty-print a comparison table across all logged runs with rich.

    ``metric_label`` selects which method label inside each run's ``metrics``
    dict to display (e.g., 'scPathWalk' or 'scVI_baseline').
    """
    records = []
    p = Path(path)
    if p.exists():
        for line in p.read_text(encoding="utf-8").splitlines():
            line = line.strip()
            if line:
                records.append(json.loads(line))

    if not records:
        console.print("[yellow]No runs recorded yet.[/yellow]")
        return

    if metric_keys is None:
        metric_keys = sorted({
            k for r in records for k in r.get("metrics", {}).get(metric_label, {})
        })

    table = Table(title=f"scPathWalk runs ({metric_label})",
                  show_header=True, header_style="bold magenta")
    table.add_column("dataset")
    table.add_column("status")
    for mk in metric_keys:
        table.add_column(mk, justify="right")
    table.add_column("tag", style="dim")
    table.add_column("run_id", style="dim")

    for r in records:
        m = r.get("metrics", {}).get(metric_label, {})
        row = [r.get("dataset", "?"), r.get("status", "?")]
        for mk in metric_keys:
            v = m.get(mk)
            row.append(f"{float(v):.4f}" if _is_number(v) else "-")
        row.append(str(r.get("tag") or ""))
        row.append(r.get("run_id", "?"))
        table.add_row(*row)

    console.print(table)

