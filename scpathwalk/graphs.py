"""Graph construction: correlation, PPI (STRING DB), and KEGG pathway graphs."""

import numpy as np
import networkx as nx
import pandas as pd
import requests
from scipy.sparse import issparse


def build_correlation_graph(adata, cutoff: float = 0.01) -> nx.Graph:
    """Build a gene-gene correlation graph from expression data.

    Edges connect genes whose Pearson correlation exceeds *cutoff*.
    """
    matrix = adata.X.todense() if issparse(adata.X) else adata.X
    corr_matrix = np.corrcoef(matrix.T)
    gene_names = adata.var_names.tolist()

    G = nx.Graph()
    n = len(gene_names)
    for i in range(n):
        for j in range(i + 1, n):
            if abs(corr_matrix[i, j]) > cutoff:
                G.add_edge(gene_names[i], gene_names[j], weight=corr_matrix[i, j])
    return G


def build_ppi_graph(
    gene_list: list[str],
    species: int = 9606,
    required_score: int = 400,
) -> nx.Graph:
    """Fetch a PPI network from STRING DB for the given genes."""
    url = "https://version-11-5.string-db.org/api/tsv-no-header/network"
    params = {
        "identifiers": "%0d".join(gene_list),
        "species": species,
        "required_score": required_score,
        "network_type": "functional",
        "show_query_node_labels": 1,
        "caller_identity": "scpathwalk",
    }
    resp = requests.post(url, data=params, timeout=120)
    resp.raise_for_status()

    edges = []
    for line in resp.text.strip().split("\n"):
        parts = line.strip().split("\t")
        p1, p2 = parts[2], parts[3]
        exp_score = float(parts[10])
        if exp_score > 0:
            edges.append({"src": p1, "dst": p2, "value": exp_score})

    if not edges:
        return nx.Graph()

    df = pd.DataFrame(edges)
    return nx.from_pandas_edgelist(df, "src", "dst", "value")


def build_pathway_graph(
    gene_list: list[str], gene_sets: str = "KEGG_2021_Human"
) -> tuple[nx.Graph, pd.DataFrame]:
    """Query Enrichr for KEGG pathways and build a clique graph.

    Returns the pathway graph and the enrichment results DataFrame.
    """
    import gseapy

    enr = gseapy.enrichr(
        gene_list=gene_list, gene_sets=gene_sets, outdir=None, no_plot=True
    )
    path_df = enr.results[["Term", "P-value", "Genes"]].copy()

    G = nx.Graph()
    for _, row in path_df.iterrows():
        genes = row["Genes"].split(";")
        clique = nx.complete_graph(genes)
        G = nx.compose(G, clique)

    return G, path_df


def merge_graphs(G_corr: nx.Graph, G_ppi: nx.Graph, G_path: nx.Graph) -> nx.Graph:
    """Merge the three gene-gene graphs into a single union graph, labelling layers."""
    for u, v in G_corr.edges():
        G_corr[u][v]["layer"] = "correlation"
    for u, v in G_ppi.edges():
        G_ppi[u][v]["layer"] = "PPI"
    for u, v in G_path.edges():
        G_path[u][v]["layer"] = "Pathways"

    G = nx.compose_all([G_corr, G_ppi, G_path])
    return G
