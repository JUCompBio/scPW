"""
Graph construction: co-expression, PPI, pathway, TF-regulon and complex layers.

Curated layers are read from the staged cache in ``data/interactome`` (see
``scpathwalk.interactome``) rather than fetched live, so a benchmarking grid is
reproducible and does not depend on API availability. The live builders are kept
for provenance but are no longer on the default path.
"""

import time

import numpy as np
import networkx as nx
import pandas as pd
import requests
from scipy.sparse import issparse

from scpathwalk.interactome import (
    edges_to_graph,
    load_source,
    sets_to_bipartite,
    sets_to_weighted_cliques,
)


def build_correlation_graph(
    adata, cutoff: float = 0.01, signed: bool = False
) -> nx.Graph | dict[str, nx.Graph]:
    """Build gene-gene co-expression graph(s) from expression data.

    Vectorised over the upper triangle rather than looping in Python: at 4000
    genes the loop version evaluates ~8M iterations per run.

    Args:
        adata: expression object; ``adata.X`` is used.
        cutoff: absolute Pearson correlation threshold.
        signed: if True, return ``{"correlation": G_pos, "anticorrelation": G_neg}``
            so that repression is represented as its own relation instead of being
            merged into co-expression by an absolute value. Negative weights are
            meaningless in a random walk, so the magnitude is stored and the sign
            is carried by the layer.

    Returns:
        A single graph, or a dict of graphs when ``signed=True``.
    """
    matrix = adata.X.todense() if issparse(adata.X) else adata.X
    corr = np.corrcoef(np.asarray(matrix).T)
    genes = np.asarray(adata.var_names.tolist())

    iu = np.triu_indices(corr.shape[0], k=1)
    r = corr[iu]
    finite = np.isfinite(r)
    src, dst, r = genes[iu[0]][finite], genes[iu[1]][finite], r[finite]

    keep = np.abs(r) > cutoff
    src, dst, r = src[keep], dst[keep], r[keep]

    def _mk(mask: np.ndarray) -> nx.Graph:
        df = pd.DataFrame(
            {"src": src[mask], "dst": dst[mask], "weight": np.abs(r[mask])}
        )
        G = nx.Graph()
        # var_names order is already deterministic, but sort anyway so every layer
        # in the multiplex shares one node ordering.
        G.add_nodes_from(sorted(genes.tolist()))
        if not df.empty:
            G.add_weighted_edges_from(df.itertuples(index=False, name=None))
        return G

    if not signed:
        return _mk(np.ones(len(r), dtype=bool))
    return {"correlation": _mk(r > 0), "anticorrelation": _mk(r < 0)}


def _resolve_string_ids(gene_list: list[str], species: int = 9606) -> pd.DataFrame:
    """Map gene names to STRING DB identifiers."""
    url = "https://version-11-5.string-db.org/api/tsv-no-header/get_string_ids"
    params = {
        "identifiers": "\r".join(gene_list),
        "species": species,
        "limit": 1,
        "echo_query": 1,
        "caller_identity": "scpathwalk",
    }
    resp = requests.post(url, data=params, timeout=120)
    resp.raise_for_status()

    rows = []
    for line in resp.text.strip().split("\n"):
        parts = line.split("\t")
        rows.append({"GeneName": parts[0], "STRING": parts[2]})
    return pd.DataFrame(rows)


def build_ppi_graph(
    gene_list: list[str],
    species: int = 9606,
    required_score: int = 400,
    batch_size: int = 1500,
    max_retries: int = 4,
) -> nx.Graph:
    """Fetch a PPI network from STRING DB for the given genes.

    Genes are queried in batches (STRING limits the number of identifiers per
    request) and each batch is retried with exponential backoff to tolerate
    transient DNS / network failures.
    """
    url = "https://version-11-5.string-db.org/api/tsv-no-header/network"

    def _post_batch(batch: list[str]) -> str:
        params = {
            "identifiers": "%0d".join(batch),
            "species": species,
            "required_score": required_score,
            "network_type": "functional",
            "show_query_node_labels": 1,
            "caller_identity": "scpathwalk",
        }
        last_exc = None
        for attempt in range(max_retries):
            try:
                resp = requests.post(url, data=params, timeout=180)
                resp.raise_for_status()
                return resp.text
            except requests.exceptions.RequestException as exc:
                last_exc = exc
                wait = 2 ** attempt
                print(f"  STRING request failed (attempt {attempt + 1}/{max_retries}); "
                      f"retrying in {wait}s...")
                time.sleep(wait)
        raise RuntimeError(f"STRING request failed after {max_retries} attempts") from last_exc

    edges = []
    for start in range(0, len(gene_list), batch_size):
        batch = gene_list[start:start + batch_size]
        text = _post_batch(batch)
        for line in text.strip().split("\n"):
            if not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) < 11:
                continue
            p1, p2 = parts[2], parts[3]
            exp_score = float(parts[10])
            if exp_score > 0:
                edges.append({"src": p1, "dst": p2, "value": exp_score})

    if not edges:
        return nx.Graph()

    df = pd.DataFrame(edges)
    return nx.from_pandas_edgelist(df, "src", "dst", "value")


def build_pathway_graph(gene_list: list[str], gene_sets: str = "KEGG_2021_Human") -> tuple[nx.Graph, pd.DataFrame]:
    """Query Enrichr for KEGG pathways and build a clique graph.

    Returns the pathway graph and the enrichment results DataFrame.
    """
    import gseapy

    enr = gseapy.enrichr(gene_list=gene_list, gene_sets=gene_sets, outdir=None, no_plot=True)
    path_df = enr.results[["Term", "P-value", "Genes"]].copy()

    G = nx.Graph()
    for _, row in path_df.iterrows():
        genes = row["Genes"].split(";")
        clique = nx.complete_graph(genes)
        G = nx.compose(G, clique)

    return G, path_df


# ── Cached curated layers ─────────────────────────────────────────────────────
CURATED_LAYERS = {
    "PPI": ("string_ppi", "edges", None),
    "Pathways": ("kegg_pathways", "sets", "term"),
    "TFRegulons": ("collectri_tf", "edges", None),
    "Complexes": ("corum_complexes", "sets", "complex"),
}


def build_curated_layers(
    gene_list: list[str],
    include: tuple[str, ...] = ("PPI", "Pathways", "TFRegulons", "Complexes"),
    pathway_repr: str = "bipartite",
) -> dict[str, nx.Graph]:
    """Build curated layers from the staged cache, restricted to ``gene_list``.

    Args:
        gene_list: gene universe (HVGs). Matched case-insensitively.
        include: which curated layers to build.
        pathway_repr: ``"bipartite"`` (gene<->set, |S| edges per set) or
            ``"clique"`` (1/|S|-weighted clique expansion). Bipartite is the
            default because unweighted cliques make a gene's importance scale with
            the size of the sets it belongs to: a 200-gene pathway contributes
            19,900 edges against 45 for a 10-gene pathway.

    Returns:
        Mapping of layer name to graph, ready for :func:`fuse_multiplex`.
    """
    universe = {g.upper() for g in gene_list}
    out: dict[str, nx.Graph] = {}

    for layer in include:
        if layer not in CURATED_LAYERS:
            raise KeyError(f"unknown curated layer '{layer}'")
        source, kind, id_col = CURATED_LAYERS[layer]
        df = load_source(source)
        if kind == "edges":
            out[layer] = edges_to_graph(df, universe)
        elif pathway_repr == "clique":
            out[layer] = sets_to_weighted_cliques(df, universe, id_col)
        else:
            out[layer] = sets_to_bipartite(df, universe, id_col)
    return out


def pathway_membership_table(gene_list: list[str]) -> pd.DataFrame:
    """KEGG pathway membership from the cache, shaped like the Enrichr output.

    Provides the ``Term`` / ``Genes`` columns that the downstream pathway
    annotation and plotting code expects, without a live Enrichr call.
    """
    universe = {g.upper() for g in gene_list}
    kegg = load_source("kegg_pathways")
    sub = kegg[kegg["gene"].isin(universe)]
    rows = [
        {"Term": term, "P-value": np.nan, "Genes": ";".join(sorted(set(genes))),
         "n_genes": len(set(genes))}
        for term, genes in sub.groupby("term")["gene"]
        if len(set(genes)) >= 2
    ]
    return pd.DataFrame(rows).sort_values("n_genes", ascending=False, ignore_index=True)


def fuse_multiplex(
    layers: dict[str, nx.Graph],
    layer_weights: dict[str, float] | None = None,
    normalize: str = "mass",
) -> tuple[nx.Graph, pd.DataFrame]:
    """Fuse named layers into one graph with per-layer normalisation.

    Summing adjacency matrices lets the densest layer dominate. Measured on a
    whole-transcriptome panel, a co-expression layer at a low correlation threshold
    contributed 98.16% of all edges against 0.13% for PPI: a 780-fold imbalance that
    per-edge weights of order 1 cannot offset, leaving the curated layers
    effectively inert.

    Following the multiplex formulation of Valdeolivas et al. (Bioinformatics
    2019), each layer is normalised *before* fusion so its total contribution is
    comparable, and the requested ``layer_weights`` then act on already-comparable
    quantities.

    Args:
        layers: mapping of layer name to graph.
        layer_weights: relative contribution per layer after normalisation.
        normalize:
            ``"mass"``  — scale each layer so its edge weights sum to 1, making
                          the layer's total contribution equal to its weight.
            ``"degree"`` — symmetric degree normalisation within each layer
                          (w_uv / sqrt(d_u * d_v)), which additionally damps hubs.
            ``"none"``  — legacy behaviour, retained for the ablation that
                          demonstrates the dominance problem.

    Returns:
        ``(G_union, contribution_report)`` where the report gives, per layer, the
        edge count and the share of total fused weight it ends up carrying. That
        share is the number to check: if one layer still holds ~98%, fusion has
        not worked.
    """
    if layer_weights is None:
        layer_weights = {}

    scaled: dict[str, dict[tuple[str, str], float]] = {}
    for name, G in layers.items():
        if G is None or G.number_of_edges() == 0:
            scaled[name] = {}
            continue

        raw = {
            (u, v): float(d.get("weight", 1.0) or 1.0) for u, v, d in G.edges(data=True)
        }

        if normalize == "mass":
            total = sum(abs(w) for w in raw.values()) or 1.0
            raw = {e: w / total for e, w in raw.items()}
        elif normalize == "degree":
            deg: dict[str, float] = {}
            for (u, v), w in raw.items():
                deg[u] = deg.get(u, 0.0) + abs(w)
                deg[v] = deg.get(v, 0.0) + abs(w)
            raw = {
                (u, v): w / np.sqrt(max(deg[u], 1e-12) * max(deg[v], 1e-12))
                for (u, v), w in raw.items()
            }
            total = sum(abs(w) for w in raw.values()) or 1.0
            raw = {e: w / total for e, w in raw.items()}
        elif normalize != "none":
            raise ValueError(f"unknown normalize mode '{normalize}'")

        lw = layer_weights.get(name, 1.0)
        scaled[name] = {e: w * lw for e, w in raw.items()}

    # Sorted node insertion so the fused graph's ordering does not depend on
    # Python's per-process string hash seed; igraph derives vertex indices from
    # insertion order and Leiden breaks ties on them.
    G_union = nx.Graph()
    all_nodes: dict[str, dict] = {}
    for name, G in layers.items():
        if G is not None:
            for node, data in G.nodes(data=True):
                all_nodes.setdefault(node, dict(data))
    for node in sorted(all_nodes):
        G_union.add_node(node, **all_nodes[node])

    for name, edges in scaled.items():
        for (u, v), w in sorted(edges.items()):
            if G_union.has_edge(u, v):
                G_union[u][v]["weight"] += w
                G_union[u][v]["layers"].append(name)
            else:
                G_union.add_edge(u, v, weight=w, layers=[name])

    grand_total = sum(abs(w) for edges in scaled.values() for w in edges.values()) or 1.0
    report = pd.DataFrame(
        [
            {
                "layer": name,
                "edges": len(edges),
                "weight_sum": sum(abs(w) for w in edges.values()),
                "weight_share": sum(abs(w) for w in edges.values()) / grand_total,
                "edge_share": len(edges)
                / max(sum(len(e) for e in scaled.values()), 1),
            }
            for name, edges in scaled.items()
        ]
    ).sort_values("weight_share", ascending=False, ignore_index=True)

    return G_union, report
