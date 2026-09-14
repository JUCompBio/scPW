"""
Stage every external interactome source into the local cache.

Run this once before the first run. Everything else reads only the cache, so runs
are reproducible and cannot be perturbed by a transient API failure or an upstream
database release. Each source is written to ``data/interactome/<name>.parquet`` with
its row count and SHA-256 recorded in ``data/interactome_manifest.json``.

The training sources shape the gene network. The hold-out sources are never read
during training and exist only to evaluate results afterwards; the contamination
check below reports how much they overlap, since a hold-out that merely restates a
training source would not be independent evidence.

Usage:
    uv run fetch_interactomes.py
    uv run fetch_interactomes.py --refresh                     # re-fetch all
    uv run fetch_interactomes.py --only string_ppi             # one source
    uv run fetch_interactomes.py --summarize data/sample.h5ad  # layer report
"""

from __future__ import annotations

import argparse
import json
import traceback
from pathlib import Path

import anndata as ad  # noqa: E402

from scpathwalk.interactome import (  # noqa: E402
    FETCHERS,
    HOLDOUT_SOURCES,
    SPECS,
    TRAINING_SOURCES,
    ensure_source,
    manifest_hash,
    summarize_layers,
    verify_disjoint,
)

ORDER = list(TRAINING_SOURCES) + list(HOLDOUT_SOURCES)


def stage(names: list[str], refresh: bool) -> dict[str, str]:
    status: dict[str, str] = {}
    for name in names:
        spec = SPECS[name]
        print(f"\n=== {name}  [{spec.role}/{spec.kind}] ===")
        print(f"    {spec.primary_citation}")
        try:
            df = ensure_source(name, refresh=refresh)
            print(f"    staged {len(df):,} rows -> data/interactome/{name}.parquet")
            status[name] = "ok"
        except Exception as exc:  # noqa: BLE001
            print(f"    FAILED: {type(exc).__name__}: {exc}")
            traceback.print_exc(limit=2)
            status[name] = f"failed: {type(exc).__name__}"
    return status


def main() -> int:
    ap = argparse.ArgumentParser(description="Stage interactome sources")
    ap.add_argument("--refresh", action="store_true", help="re-fetch even if cached")
    ap.add_argument("--only", nargs="*", choices=sorted(FETCHERS), help="subset to fetch")
    ap.add_argument("--summarize", metavar="H5AD", type=Path,
                    help="report per-layer coverage on this dataset's genes")
    args = ap.parse_args()

    names = args.only if args.only else ORDER
    status = stage(names, refresh=args.refresh)

    print("\n=== contamination check (training vs hold-out) ===")
    report = verify_disjoint(strict=False)
    if not report["pairwise"]:
        print("    nothing to compare yet")
    for key, vals in sorted(report["pairwise"].items()):
        print(
            f"    {key:38s} inter={vals['intersection']:>9,}  "
            f"jaccard={vals['jaccard']:.4f}  frac_holdout={vals['frac_of_holdout']:.4f}"
        )
    Path("data").mkdir(exist_ok=True)
    Path("data/contamination_report.json").write_text(json.dumps(report, indent=2))

    if args.summarize:
        if not args.summarize.exists():
            print(f"\n{args.summarize} not found")
        else:
            adata = ad.read_h5ad(args.summarize, backed="r")
            universe = {g.upper() for g in adata.var_names}
            print(f"\n=== layers on {args.summarize.name} "
                  f"({len(universe):,} genes) ===")
            print(summarize_layers(universe).to_string(index=False))

    print(f"\nmanifest hash: {manifest_hash()}")
    failed = {k: v for k, v in status.items() if v != "ok"}
    if failed:
        print(f"FAILED: {failed}")
        return 1
    print("all sources staged")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
