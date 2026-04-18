"""
run_pipeline.py -- Main scPathWalk pipeline.

Runs the full analysis for a single dataset:
  1. Load preprocessed .h5ad and select HVGs
  2. Build gene-gene graphs (correlation, PPI, KEGG pathway)
  3. Leiden-partition each graph, then build weighted consensus
  4. PageRank on the merged union graph -> gene importance
  5. Reweight expression matrix
  6. scVI embedding -> Leiden clustering
  7. Evaluate (Silhouette / Calinski-Harabasz / Davies-Bouldin)
  8. Visualize (UMAP, ranked genes, pathway overlay)
  9. Optionally run a vanilla scVI baseline for comparison

Usage:
    uv run run_pipeline.py --dataset ovary
    uv run run_pipeline.py --dataset 9k_cervical_cancer --no-baseline
"""

import argparse
import warnings
from pathlib import Path

import anndata as ad
import networkx as nx
import scanpy as sc

from util import process_adata
from scpathwalk.graphs import (
    build_correlation_graph,
    build_ppi_graph,
    build_pathway_graph,
    merge_graphs,
)
from scpathwalk.partition import (
    leiden_partition,
    map_partitions_to_genes,
    weighted_consensus_clustering,
)
from scpathwalk.reweight import compute_pagerank, reweight_expression
from scpathwalk.embed import setup_and_train, embed_and_cluster
from scpathwalk.evaluate import compute_metrics
from scpathwalk.visualize import (
    plot_clusters,
    rank_and_save_genes,
    plot_top_gene_signals,
    annotate_pathway_terms,
    plot_pathway_umap,
)

warnings.filterwarnings("ignore")

# -- Default parameters --------------------------------------------------------
DEFAULTS = dict(
    n_top_genes=1200,
    corr_cutoff=0.01,
    string_score=400,
    pagerank_alpha=0.9,
    consensus_weights=[0.7, 0.7, 1.0],
    n_latent=5,
    n_epochs=30,
    leiden_resolution=0.55,
    baseline_epochs=10,
    baseline_resolution=0.4,
)

LABEL_COLUMN_CANDIDATES = ["labels", "celltype.l2", "cell_types"]


def _assign_labels(adata_src, adata_dst):
    """Try to copy cell-type labels from the source to destination AnnData."""
    for col in LABEL_COLUMN_CANDIDATES:
        if col in adata_src.obs.columns:
            adata_dst.obs["labels"] = adata_src.obs[col]
            if col != "labels":
                adata_src.obs["labels"] = adata_src.obs[col]
            return
    print("  (no cell-type labels found -- skipping label transfer)")


def run(dataset: str, run_baseline: bool = True, **params):
    """Execute the full scPathWalk pipeline for one dataset."""
    p = {**DEFAULTS, **params}
    processed_dir = Path("processed")
    h5ad_path = processed_dir / f"{dataset}.h5ad"

    sc.settings.set_figure_params(dpi=180, dpi_save=600, format="pdf")

    # 1. Load & HVG subset
    print(f"\n{'='*70}")
    print(f"  scPathWalk pipeline -- {dataset}")
    print(f"{'='*70}")

    adata = sc.read_h5ad(h5ad_path)
    print(f"Loaded {h5ad_path}: {adata.n_obs} cells x {adata.n_vars} genes")

    adata_proc = process_adata(adata, n_top_genes=p["n_top_genes"])
    gene_list = adata_proc.var_names.tolist()
    print(f"After HVG selection: {len(gene_list)} genes")

    # 2. Correlation graph
    print("\n[CORR] Building correlation graph...")
    G_corr = build_correlation_graph(adata_proc, cutoff=p["corr_cutoff"])
    print(f"  {G_corr.number_of_nodes()} nodes, {G_corr.number_of_edges()} edges")

    _, corr_map = leiden_partition(G_corr)
    adata_proc.var["corr_partitions"] = map_partitions_to_genes(gene_list, corr_map)

    # 3. PPI graph
    print("\n[PPI] Fetching STRING DB network...")
    corr_genes = list(G_corr.nodes())
    G_ppi = build_ppi_graph(corr_genes, required_score=p["string_score"])
    print(f"  {G_ppi.number_of_nodes()} nodes, {G_ppi.number_of_edges()} edges")

    _, ppi_map = leiden_partition(G_ppi)
    adata_proc.var["ppi_partitions"] = map_partitions_to_genes(gene_list, ppi_map)

    # 4. Pathway graph
    print("\n[PW] Querying KEGG pathways via Enrichr...")
    G_path, path_df = build_pathway_graph(corr_genes)
    print(f"  {G_path.number_of_nodes()} nodes, {G_path.number_of_edges()} edges")

    _, pw_map = leiden_partition(G_path)
    adata_proc.var["pathway_partitions"] = map_partitions_to_genes(gene_list, pw_map)

    # 5. Consensus partition
    print("\n[CONSENSUS] Building weighted consensus partition...")
    partition_lists = [
        adata_proc.var["corr_partitions"].tolist(),
        adata_proc.var["ppi_partitions"].tolist(),
        adata_proc.var["pathway_partitions"].tolist(),
    ]
    consensus = weighted_consensus_clustering(partition_lists, p["consensus_weights"])
    adata_proc.var["weighted_consensus_partitions"] = consensus

    # 6. Merge graphs & PageRank
    print("\n[RWR] Merging graphs and running PageRank...")
    G_union = merge_graphs(G_corr, G_ppi, G_path)
    gene_weights = compute_pagerank(G_union, alpha=p["pagerank_alpha"])
    print(f"  {len(gene_weights)} gene weights computed")

    # 7. Reweight expression
    print("\n[REWEIGHT] Scaling expression by PageRank scores...")
    gexp = reweight_expression(adata, gene_weights)

    # 8. Build rescaled AnnData
    adata_rescaled = ad.AnnData(gexp)
    _assign_labels(adata, adata_rescaled)
    adata_rescaled.var["gene_ids"] = gexp.columns.tolist()

    # 9. scVI embedding + clustering
    print("\n[EMBED] Training scVI autoencoder...")
    model = setup_and_train(
        adata_rescaled, n_latent=p["n_latent"], n_epochs=p["n_epochs"]
    )
    adata_rescaled = embed_and_cluster(
        adata_rescaled, model, resolution=p["leiden_resolution"]
    )

    # 10. Evaluate
    metrics = compute_metrics(adata_rescaled)
    print(f"\n[METRICS] scPathWalk ({dataset}):")
    for k, v in metrics.items():
        print(f"  {k:>20s}: {v:.4f}")

    # 11. Visualize
    print("\n[VIZ] Generating plots...")
    plot_clusters(adata_rescaled, dataset)
    rank_and_save_genes(adata_rescaled, dataset)
    plot_top_gene_signals(adata_rescaled, dataset)

    print("[VIZ] Annotating pathway terms...")
    annotate_pathway_terms(adata_rescaled, path_df)
    plot_pathway_umap(adata_rescaled, dataset)

    # 12. Vanilla baseline (optional)
    if run_baseline:
        print(f"\n{'~'*70}")
        print(f"  Vanilla scVI baseline -- {dataset}")
        print(f"{'~'*70}")
        model_base = setup_and_train(
            adata, n_latent=p["n_latent"], n_epochs=p["baseline_epochs"]
        )
        adata = embed_and_cluster(
            adata, model_base, resolution=p["baseline_resolution"]
        )
        base_metrics = compute_metrics(adata)
        print(f"\n[METRICS] Vanilla scVI ({dataset}):")
        for k, v in base_metrics.items():
            print(f"  {k:>20s}: {v:.4f}")
        sc.pl.umap(
            adata,
            color=["leiden_scVI"],
            frameon=False,
            save=f"_{dataset}_defaultcluster",
        )

    print(f"\nDone -- {dataset}")
    return adata_rescaled


def main():
    all_datasets = [
        "3praw",
        "4plex",
        "9k_cervical_cancer",
        "GSE161529",
        "GSE228499",
        "GSE266351",
        "ovary",
    ]
    parser = argparse.ArgumentParser(description="Run the scPathWalk pipeline")
    parser.add_argument(
        "--dataset",
        required=True,
        choices=all_datasets,
        help="Name of the preprocessed dataset to analyse",
    )
    parser.add_argument(
        "--no-baseline",
        action="store_true",
        help="Skip the vanilla scVI baseline comparison",
    )
    args = parser.parse_args()
    run(args.dataset, run_baseline=not args.no_baseline)


if __name__ == "__main__":
    main()
