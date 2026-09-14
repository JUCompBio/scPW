"""
Gene adjacency for the network-informed dispersion prior.

scVI learns one inverse dispersion parameter per gene, freely and independently. A
gene detected in a small fraction of cells gets a badly determined value, and on a
typical 4,000-gene panel the median cell carries readings for a few dozen genes, so
most genes are in that regime.

Classical count modelling handles this by partial pooling: edgeR and DESeq2 shrink
each gene's dispersion toward a common trend, and single-cell versions pool across
genes of similar expression. This module supplies the ingredient for pooling by
**network adjacency** instead, so a poorly measured gene borrows from its
interactors rather than from unrelated genes that happen to be equally expressed.

Two decisions here are load-bearing.

**Which layers to include.** A data-derived co-expression layer at a low correlation
threshold reaches a density around 0.74, so nearly every gene pair is a neighbour.
Averaging over that many neighbours approximates a global mean, which turns the
penalty into ordinary uniform shrinkage: the pooling is no longer *network*-specific.
The default is therefore the curated layers only. It is also far cheaper, dropping
millions of edges to order 10^5.

**Gene order.** The matrix multiplies the model's per-gene parameter vector, so it
must follow ``adata.var_names`` exactly. A misaligned matrix would still train and
still produce a falling penalty while pairing every gene with the wrong neighbours,
which is the sort of failure that yields confident nonsense.
"""

from __future__ import annotations

import networkx as nx
import numpy as np
import scipy.sparse as sp

# Curated biological relationships. Co-expression and anti-correlation are excluded
# by default: they are data-derived, and at low thresholds dense enough to wash out
# the curated structure.
CURATED_LAYERS = ("PPI", "Pathways", "TFRegulons", "Complexes")


def gene_adjacency(
    graph: nx.Graph,
    var_names: list[str],
    normalize: bool = True,
    layers: tuple[str, ...] | None = CURATED_LAYERS,
    max_density: float = 0.05,
) -> sp.csr_matrix:
    """Sparse gene-gene adjacency aligned to ``var_names``.

    Args:
        graph: fused multiplex. Non-gene nodes carry ``node_kind != "gene"`` and are
            dropped, since the prior is over genes.
        var_names: gene order to align to, normally ``adata.var_names``.
        normalize: row-normalise so each gene's neighbour weights sum to 1, making
            the penalty a comparison against a weighted *mean*. Without it, hub genes
            would be penalised simply for having many partners.
        layers: keep only edges whose ``layer`` attribute is in this tuple. ``None``
            keeps everything. Edges lacking a ``layer`` attribute are kept, so a
            graph built without layer bookkeeping still works.
        max_density: warn above this density, since a dense adjacency makes the
            penalty behave like uniform shrinkage rather than network shrinkage.

    Returns:
        ``(n_genes, n_genes)`` CSR matrix in ``var_names`` order.
    """
    index = {g.upper(): i for i, g in enumerate(var_names)}
    n = len(var_names)

    # Collect with a comprehension and convert via numpy. The first version looped
    # in Python over every edge, which dominated runtime at 5.9M edges.
    #
    # fuse_multiplex records the contributing layers per edge as a LIST under the key
    # "layers", not a scalar "layer". Reading the wrong key would make the filter a
    # no-op and silently admit all 5.9M co-expression edges, so an edge is kept only
    # when at least one of its recorded layers is wanted, and edges with no record at
    # all are kept so a graph built without layer bookkeeping still works.
    if layers is not None:
        keep = set(layers)
        edges = [
            (u, v, d.get("weight", 1.0))
            for u, v, d in graph.edges(data=True)
            if not d.get("layers") or (keep & set(d["layers"]))
        ]
    else:
        edges = [(u, v, d.get("weight", 1.0)) for u, v, d in graph.edges(data=True)]

    if not edges:
        return sp.csr_matrix((n, n), dtype=np.float32)

    m = len(edges)
    u_arr = np.fromiter((index.get(str(e[0]).upper(), -1) for e in edges),
                        dtype=np.int64, count=m)
    v_arr = np.fromiter((index.get(str(e[1]).upper(), -1) for e in edges),
                        dtype=np.int64, count=m)
    w_arr = np.fromiter((float(e[2]) for e in edges), dtype=np.float32, count=m)

    # Drop edges touching anything absent from the expression matrix (pathway set
    # nodes, genes cut by HVG selection), self-loops, and non-positive weights.
    ok = (u_arr >= 0) & (v_arr >= 0) & (u_arr != v_arr) & (w_arr > 0)
    u_arr, v_arr, w_arr = u_arr[ok], v_arr[ok], w_arr[ok]
    if u_arr.size == 0:
        return sp.csr_matrix((n, n), dtype=np.float32)

    # Symmetrise: the penalty is undirected.
    rows = np.concatenate([u_arr, v_arr])
    cols = np.concatenate([v_arr, u_arr])
    vals = np.concatenate([w_arr, w_arr])

    A = sp.coo_matrix((vals, (rows, cols)), shape=(n, n), dtype=np.float32).tocsr()
    # Several layers may connect the same pair; summing is intended, since
    # fuse_multiplex has already balanced the layer masses.
    A.sum_duplicates()

    density = A.nnz / max(n * n, 1)
    if density > max_density:
        print(f"  WARNING: adjacency density {density:.3f} exceeds {max_density}. "
              f"At this many neighbours per gene the penalty approximates a global "
              f"mean, i.e. uniform shrinkage rather than network shrinkage.")

    if normalize:
        deg = np.asarray(A.sum(axis=1)).ravel()
        inv = np.zeros_like(deg)
        nz = deg > 0
        inv[nz] = 1.0 / deg[nz]
        # Isolated genes keep an all-zero row and contribute nothing to the penalty,
        # which is correct: there are no neighbours to borrow from.
        A = sp.diags(inv, dtype=np.float32) @ A

    return A.tocsr()


def adjacency_summary(A: sp.csr_matrix, var_names: list[str]) -> dict:
    """Diagnostics to log before training.

    ``connected_frac`` is the number to read first: genes with no neighbours are
    untouched by the prior, so if it is small the mechanism cannot do much at any
    tau.
    """
    deg = np.asarray((A > 0).sum(axis=1)).ravel()
    connected = int((deg > 0).sum())
    n = len(var_names)
    return {
        "n_genes": n,
        "n_edges": int(A.nnz // 2),
        "density": float(A.nnz / max(n * n, 1)),
        "connected_genes": connected,
        "connected_frac": float(connected / max(n, 1)),
        "mean_degree": float(deg.mean()),
        "median_degree": float(np.median(deg)),
        "max_degree": int(deg.max()) if deg.size else 0,
        "isolated_genes": int((deg == 0).sum()),
    }
