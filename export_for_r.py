"""
Export a QC-filtered matrix to MTX so the R comparison methods see the same input.

This is the detail that decides whether the R comparison means anything. If Seurat
and Pagoda2 load raw files and apply their own filtering, any difference in the
results mixes the clustering method with the preprocessing.

So the export is the output of ``util.process_adata`` -- the identical cell and gene
set the Python arm trains on -- written as counts, because Seurat and Pagoda2 both
expect raw counts and do their own normalisation. Writing normalised values would
double-normalise inside R.

Counts come from ``adata.layers["counts"]`` where present, since ``adata.X`` has
already been normalised and log-transformed by the time ``process_adata`` returns.

Usage:
    uv run export_for_r.py --h5ad data/my_sample.h5ad
    uv run export_for_r.py --h5ad data/my_sample.h5ad --name my_sample
"""

from __future__ import annotations

import argparse
import logging
import os
import sys
import warnings
from pathlib import Path

os.environ.setdefault("PYTHONHASHSEED", "0")

import numpy as np  # noqa: E402
import scanpy as sc  # noqa: E402
import scipy.io as sio  # noqa: E402
import scipy.sparse as sp  # noqa: E402

from util import process_adata  # noqa: E402

warnings.filterwarnings("ignore")
sc.settings.verbosity = 0
logging.getLogger("scvi").setLevel(logging.ERROR)


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Export a QC-filtered matrix to MTX for the R comparators.")
    ap.add_argument("--h5ad", required=True, type=Path,
                    help="Input AnnData file of raw counts.")
    ap.add_argument("--name", default=None,
                    help="Export directory name. Defaults to the input file stem.")
    ap.add_argument("--out", type=Path, default=Path("data/r_export"),
                    help="Root directory for exports.")
    ap.add_argument("--n-top-genes", type=int, default=4000)
    args = ap.parse_args()

    if not args.h5ad.exists():
        print(f"error: {args.h5ad} not found", file=sys.stderr)
        return 1

    name = args.name or args.h5ad.stem
    proc = process_adata(sc.read_h5ad(args.h5ad), n_top_genes=args.n_top_genes)

    # Raw counts, not the normalised X. Seurat's NormalizeData and Pagoda2 both
    # expect counts and would otherwise normalise an already-normalised matrix.
    if "counts" in proc.layers:
        X = proc.layers["counts"]
        source = "layers['counts']"
    else:
        X = proc.X
        source = "X (no counts layer found)"

    X = sp.csr_matrix(X) if not sp.issparse(X) else X.tocsr()
    # MTX convention for 10x-style data is genes x cells.
    mat = X.T.tocoo()

    d = args.out / name
    d.mkdir(parents=True, exist_ok=True)
    integer = np.allclose(mat.data, np.round(mat.data))
    sio.mmwrite(str(d / "matrix.mtx"), mat,
                field="integer" if integer else "real")
    (d / "features.tsv").write_text("\n".join(proc.var_names) + "\n",
                                    encoding="utf-8")
    (d / "barcodes.tsv").write_text("\n".join(proc.obs_names) + "\n",
                                    encoding="utf-8")

    print(f"  {name:<20} {proc.n_obs:>6} cells x {proc.n_vars:>5} genes  "
          f"from {source}  {'integer counts' if integer else 'NON-INTEGER'}")
    if not integer:
        print("    WARNING: values are not integers, so R-side normalisation will "
              "be applied to already-transformed data")

    print(f"\nwrote {d}/{{matrix.mtx,features.tsv,barcodes.tsv}}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
