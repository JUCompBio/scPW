"""
Comparison methods.

Baselines beyond Pagoda2, scGSVA and scVI, specifically
simpler and widely adopted workflows. The two here are the standard
Python routes to a clustering, and every parameter is explicit so the comparison
is reproducible.

All methods take the same QC-filtered AnnData the pipeline uses and return
``(adata, metrics)`` with a ``leiden_scVI`` label column, so downstream metric and
coherence code needs no special cases. The column name is a historical artefact of
the scVI-only pipeline; it means "the clustering to evaluate" rather than anything
scVI-specific.
"""

from __future__ import annotations

import anndata as ad
import numpy as np
import scanpy as sc

CLUSTER_KEY = "leiden_scVI"


def scanpy_pca_leiden(
    adata: ad.AnnData,
    n_comps: int = 50,
    n_neighbors: int = 15,
    resolution: float = 0.55,
    seed: int = 420,
    scale: bool = True,
) -> ad.AnnData:
    """The conventional Scanpy workflow: scale, PCA, neighbours, Leiden.

    This is the baseline most readers have in mind when asking
    whether a deep generative model earns its complexity. Deliberately plain: no
    batch correction, no denoising, default Scanpy parameters other than the
    resolution, which is matched to every other method.
    """
    out = adata.copy()
    if scale:
        # zero_center=False keeps the matrix sparse-friendly and avoids the memory
        # blow-up on the larger datasets; max_value clips the long tail that
        # otherwise dominates the leading components.
        sc.pp.scale(out, zero_center=False, max_value=10)

    n_comps = int(min(n_comps, min(out.shape) - 1))
    sc.tl.pca(out, n_comps=n_comps, svd_solver="arpack", random_state=seed)
    out.obsm["X_pca_rep"] = out.obsm["X_pca"]

    sc.pp.neighbors(out, n_neighbors=n_neighbors, use_rep="X_pca_rep",
                    random_state=seed)
    sc.tl.leiden(out, key_added=CLUSTER_KEY, resolution=resolution,
                 random_state=seed)
    return out


def scanpy_hvg_leiden(
    adata: ad.AnnData,
    resolution: float = 0.55,
    seed: int = 420,
    n_neighbors: int = 15,
) -> ad.AnnData:
    """Leiden directly on the expression matrix, without dimensionality reduction.

    A floor rather than a competitor: it shows how much of any method's advantage
    comes from having a learned representation at all, as opposed to from the
    particular representation being proposed.
    """
    out = adata.copy()
    out.obsm["X_expr"] = (
        out.X.toarray() if hasattr(out.X, "toarray") else np.asarray(out.X)
    )
    sc.pp.neighbors(out, n_neighbors=n_neighbors, use_rep="X_expr",
                    random_state=seed)
    sc.tl.leiden(out, key_added=CLUSTER_KEY, resolution=resolution,
                 random_state=seed)
    return out


# Latent-space key per method, for metric computation.
REPRESENTATION = {
    "scanpy_pca_leiden": "X_pca_rep",
    "scanpy_hvg_leiden": "X_expr",
}

BASELINES = {
    "scanpy_pca_leiden": scanpy_pca_leiden,
    "scanpy_hvg_leiden": scanpy_hvg_leiden,
}


def run_baseline(name: str, adata: ad.AnnData, **kwargs) -> tuple[ad.AnnData, str]:
    """Run a named baseline, returning the result and its representation key."""
    if name not in BASELINES:
        raise KeyError(f"unknown baseline '{name}'; known: {sorted(BASELINES)}")
    return BASELINES[name](adata, **kwargs), REPRESENTATION[name]
