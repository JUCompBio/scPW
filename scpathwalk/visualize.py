"""Visualization: UMAP plots, ranked genes, pathway overlays."""

import math
from collections import Counter
from pathlib import Path

import pandas as pd
import scanpy as sc


def plot_clusters(adata, friendly_name: str) -> None:
    """UMAP coloured by Leiden clusters."""
    sc.pl.umap(
        adata,
        color=["leiden_scVI"],
        frameon=False,
        save=f"_{friendly_name}_pathwaycluster",
    )


def rank_and_save_genes(adata, friendly_name: str, rank_dir: str = "rank") -> None:
    """Rank genes per cluster (Wilcoxon) and save CSV files."""
    Path(rank_dir).mkdir(exist_ok=True)

    adata.var_names_make_unique()
    sc.tl.rank_genes_groups(adata, groupby="leiden_scVI", method="wilcoxon")

    sc.pl.rank_genes_groups_dotplot(
        adata,
        groupby="leiden_scVI",
        standard_scale="var",
        n_genes=5,
        save=f"{friendly_name}_ranking",
    )

    cluster_labels = sorted(set(adata.obs["leiden_scVI"].tolist()))
    frames = []
    for label in cluster_labels:
        df = sc.get.rank_genes_groups_df(adata, group=label).head(50)
        df["label"] = str(label)
        frames.append(df)
    pd.concat(frames).to_csv(f"{rank_dir}/{friendly_name}_ranks.csv", index=False)

    sc.tl.rank_genes_groups(
        adata, "leiden_scVI", method="t-test", key_added="ranked_genes"
    )
    sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, key="ranked_genes")
    ranked = sc.get.rank_genes_groups_df(adata, None, key="ranked_genes")
    ranked.to_csv(
        f"{rank_dir}/{friendly_name}_full_clustering_ranked.csv", index=False
    )


def plot_top_gene_signals(adata, friendly_name: str) -> None:
    """UMAP overlaid with top genes from cluster 0."""
    genes = sc.get.rank_genes_groups_df(adata, group="0").head(5)["names"]
    sc.pl.umap(
        adata,
        color=[*genes, "leiden_scVI"],
        legend_loc="on data",
        frameon=False,
        ncols=3,
        save=f"_{friendly_name}_signal",
    )


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

    merged = pd.merge(
        adata.var, top_pw, left_on="gene_ids", right_on="Genes", how="left"
    )

    pw_counts = merged["Term"].value_counts().reset_index()
    pw_counts.columns = ["Term", "Pathway_Count"]
    merged = merged.merge(pw_counts, on="Term", how="left")

    has_pw = merged.dropna(subset=["Pathway_Count"])
    if not has_pw.empty:
        best = has_pw.loc[has_pw.groupby("gene_ids")["Pathway_Count"].idxmax()]
        best = best.drop(columns=["Pathway_Count"])
    else:
        best = merged.drop(columns=["Pathway_Count"]).head(0)

    final = pd.merge(adata.var[["gene_ids"]], best, on="gene_ids", how="left")
    final.index = adata.var.index
    adata.var = final


def plot_pathway_umap(adata, friendly_name: str, n_terms: int = 11) -> None:
    """UMAP coloured by top pathway terms."""
    terms = adata.var["Term"].dropna().tolist()
    top_terms = [t for t, _ in Counter(terms).most_common(n_terms)]

    if not top_terms:
        print("No pathway terms found -- skipping pathway UMAP.")
        return

    sc.pl.umap(
        adata,
        color=[*top_terms, "leiden_scVI"],
        gene_symbols="Term",
        legend_loc="on data",
        frameon=False,
        ncols=3,
        save=f"_{friendly_name}_pathway",
    )
