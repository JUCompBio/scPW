"""PageRank-based gene reweighting of the expression matrix."""

import numpy as np
import pandas as pd
import networkx as nx
from scipy.sparse import issparse
from sklearn.preprocessing import MinMaxScaler


def compute_pagerank(G: nx.Graph, alpha: float = 0.9) -> pd.DataFrame:
    """Run PageRank on the union graph and return a DataFrame of gene weights.

    Long node names (>15 chars, typically pathway names) are dropped.
    """
    try:
        pr = nx.pagerank(G, alpha=alpha, max_iter=1000)
    except nx.PowerIterationFailedConvergence:
        print("PageRank convergence failed -- reducing tolerance and retrying...")
        pr = nx.pagerank(G, alpha=alpha, max_iter=10000, tol=1e-3)

    gene_weights = pd.DataFrame(pr.items(), columns=["term", "value"])
    gene_weights = gene_weights[gene_weights["term"].str.len() <= 15].reset_index(
        drop=True
    )
    return gene_weights


def reweight_expression(adata, gene_weights: pd.DataFrame) -> pd.DataFrame:
    """Scale the original expression matrix by MinMax-scaled PageRank weights.

    Returns a dense DataFrame of reweighted expression values.
    """
    X = adata.X.todense() if issparse(adata.X) else adata.X
    gexp = pd.DataFrame(X, columns=adata.var_names.tolist())

    shared = [c for c in gexp.columns if c in gene_weights["term"].values]
    gexp = gexp[shared]

    scaler = MinMaxScaler()
    gene_weights = gene_weights.copy()
    gene_weights["value"] = scaler.fit_transform(gene_weights[["value"]])

    weight_map = dict(zip(gene_weights["term"], gene_weights["value"]))
    for col in gexp.columns:
        if col in weight_map:
            gexp[col] = gexp[col] * weight_map[col]

    return gexp
