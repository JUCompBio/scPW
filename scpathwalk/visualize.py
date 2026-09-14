"""
Visualization: UMAP plots, ranked genes, pathway overlays.

The Agg backend is selected before pyplot is imported so figures are rendered to
file and never to a window. Pipeline runs are batch jobs; popping up windows
interrupts whatever else the machine is doing and blocks on some backends.
"""

import math
from collections import Counter
from pathlib import Path

import matplotlib

matplotlib.use("Agg", force=True)

import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402
import scanpy as sc  # noqa: E402

plt.ioff()


def plot_clusters(adata, friendly_name: str) -> None:
    """UMAP coloured by Leiden clusters."""
    sc.pl.umap(
        adata,
        color=["leiden_scVI"],
        frameon=False,
        show=False,
        save=f"_{friendly_name}_pathwaycluster",
    )
    plt.close("all")


def rank_and_save_genes(
    adata,
    friendly_name: str,
    rank_dir: str = "rank",
    n_genes: int = 50,
    dotplot: bool = True,
) -> None:
    """Rank genes per cluster (Wilcoxon) and save CSV files.

    Args:
        adata: clustered object with ``leiden_scVI`` in ``obs``.
        friendly_name: filename stem used to name the output files.
        rank_dir: output directory.
        n_genes: markers retained per cluster.
        dotplot: whether to render the ranking dot plot. Turn it off for a run
            whose markers are needed for the coherence metrics but whose figure
            is not wanted.
    """
    Path(rank_dir).mkdir(parents=True, exist_ok=True)

    adata.var_names_make_unique()
    sc.tl.rank_genes_groups(adata, groupby="leiden_scVI", method="wilcoxon")

    if dotplot:
        sc.pl.rank_genes_groups_dotplot(
            adata,
            groupby="leiden_scVI",
            standard_scale="var",
            n_genes=5,
            show=False,
            save=f"{friendly_name}_ranking",
        )
        plt.close("all")

    cluster_labels = sorted(set(adata.obs["leiden_scVI"].tolist()))
    frames = []
    for label in cluster_labels:
        df = sc.get.rank_genes_groups_df(adata, group=label).head(n_genes)
        df["label"] = str(label)
        frames.append(df)
    pd.concat(frames).to_csv(f"{rank_dir}/{friendly_name}_ranks.csv", index=False)

    # Full ranked-gene export, computed for the CSV only. The t-test ranking is
    # written out but not plotted; drawing it without a `save` target opens an
    # interactive window, so it is closed immediately.
    sc.tl.rank_genes_groups(adata, "leiden_scVI", method="t-test", key_added="ranked_genes")
    ranked = sc.get.rank_genes_groups_df(adata, None, key="ranked_genes")
    ranked.to_csv(f"{rank_dir}/{friendly_name}_full_clustering_ranked.csv", index=False)
    plt.close("all")


def plot_top_gene_signals(adata, friendly_name: str, n_per_cluster: int = 1) -> None:
    """UMAP overlaid with cluster-discriminative top marker genes.

    Instead of showing the top genes of a single cluster, this selects the top
    ``n_per_cluster`` differentially expressed genes from *each* Leiden cluster,
    so the overlay illustrates cluster-specific expression structure.
    """
    cluster_labels = sorted(set(adata.obs["leiden_scVI"].astype(str)))
    genes: list[str] = []
    for label in cluster_labels:
        df = sc.get.rank_genes_groups_df(adata, group=label).head(n_per_cluster)
        for g in df["names"].tolist():
            if g not in genes:
                genes.append(g)
    sc.pl.umap(
        adata,
        color=[*genes, "leiden_scVI"],
        legend_loc="on data",
        frameon=False,
        ncols=3,
        show=False,
        save=f"_{friendly_name}_signal",
    )
    plt.close("all")


def annotate_pathway_terms(adata, path_df: pd.DataFrame) -> None:
    """Assign the most-frequent pathway term to each gene in adata.var."""
    cluster_labels = sorted(set(adata.obs["leiden_scVI"].tolist()))
    path_expl = path_df.assign(Genes=path_df["Genes"].str.split(";")).explode("Genes")

    frames = []
    for label in cluster_labels:
        top_genes = sc.get.rank_genes_groups_df(adata, group=label)["names"].tolist()
        filtered = path_expl[path_expl["Genes"].isin(top_genes)].copy()
        filtered["cluster"] = str(label)
        frames.append(filtered)

    top_pw = pd.concat(frames)

    merged = pd.merge(adata.var, top_pw, left_on="gene_ids", right_on="Genes", how="left")

    # For each gene keep the pathway with highest count
    pw_counts = merged["Term"].value_counts().reset_index()
    pw_counts.columns = ["Term", "Pathway_Count"]
    merged = merged.merge(pw_counts, on="Term", how="left")

    # Drop genes with no pathway match before groupby, then rejoin
    has_pw = merged.dropna(subset=["Pathway_Count"])
    if not has_pw.empty:
        best = has_pw.loc[has_pw.groupby("gene_ids")["Pathway_Count"].idxmax()]
        best = best.drop(columns=["Pathway_Count"])
    else:
        best = merged.drop(columns=["Pathway_Count"]).head(0)  # empty frame

    final = pd.merge(adata.var[["gene_ids"]], best, on="gene_ids", how="left")
    # Reset index to match adata shape
    final.index = adata.var.index
    adata.var = final


def plot_pathway_umap(adata, friendly_name: str, n_terms: int = 11) -> None:
    """UMAP coloured by top pathway terms."""
    terms = adata.var["Term"].dropna().tolist()
    top_terms = [t for t, _ in Counter(terms).most_common(n_terms)]

    if not top_terms:
        print("No pathway terms found — skipping pathway UMAP.")
        return

    sc.pl.umap(
        adata,
        color=[*top_terms, "leiden_scVI"],
        gene_symbols="Term",
        legend_loc="on data",
        frameon=False,
        ncols=3,
        show=False,
        save=f"_{friendly_name}_pathway",
    )
    plt.close("all")
