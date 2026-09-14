"""
Run scPathWalk on one dataset, with the matched control.

What this does, in order:

  1. Preprocess the counts (normalise, log1p, doublet removal, HVG subset).
  2. Build the curated gene-gene adjacency from the staged interactome cache.
  3. Train the network-regularised model and, unless disabled, the tau=0 control.
  4. Cluster the learned representation with Leiden and score both arms.

The adjacency is built from the pairwise curated layers: STRING protein
interactions and CollecTRI transcription-factor regulons. KEGG and CORUM are
set-valued -- they name groups of genes rather than gene pairs -- and under the
default bipartite representation they contribute no gene-gene edges, so they do
not enter the penalty. Pass --pathway-repr clique to expand each set into all
pairs instead; see scpathwalk.graphs.build_curated_layers.

Stage the interactome sources once before the first run:

    uv run fetch_interactomes.py

Usage:
    uv run run_scpathwalk.py --h5ad data/my_sample.h5ad
    uv run run_scpathwalk.py --h5ad data/my_sample.h5ad --tau 300 --n-latent 2
    uv run run_scpathwalk.py --h5ad data/my_sample.h5ad --no-control
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import sys
import warnings
from pathlib import Path

# Set before torch and scanpy are imported, so hash-dependent orderings are fixed.
os.environ.setdefault("PYTHONHASHSEED", "0")

import numpy as np  # noqa: E402
import scanpy as sc  # noqa: E402
import scvi  # noqa: E402

from scpathwalk.dispersion import CURATED_LAYERS, adjacency_summary, gene_adjacency  # noqa: E402
from scpathwalk.evaluate import compute_metrics  # noqa: E402
from scpathwalk.graphs import build_curated_layers, fuse_multiplex  # noqa: E402
from scpathwalk.network_vae import SCVINetworkPrior, theta_alignment  # noqa: E402
from util import process_adata  # noqa: E402

warnings.filterwarnings("ignore")
sc.settings.verbosity = 0
for _name in ("lightning.pytorch", "lightning_fabric", "pytorch_lightning", "scvi"):
    logging.getLogger(_name).setLevel(logging.ERROR)

# Equal and untuned across the pairwise layers. Complexes carry a lower weight
# only because their co-membership claim is weaker than a measured interaction;
# under the default bipartite representation both set-valued layers contribute
# nothing to the fitted adjacency in any case.
LAYER_WEIGHTS = {
    "correlation": 1.0,
    "anticorrelation": 0.5,
    "PPI": 1.0,
    "Pathways": 1.0,
    "TFRegulons": 1.0,
    "Complexes": 0.7,
}


def build_adjacency(adata, layers: tuple[str, ...], pathway_repr: str):
    """Curated gene-gene adjacency aligned to ``adata.var_names``.

    Nothing here consumes a random seed, so the result is a pure function of the
    gene panel and the staged sources.
    """
    genes = adata.var_names.tolist()
    curated = build_curated_layers(genes, include=layers, pathway_repr=pathway_repr)
    fused, report = fuse_multiplex(curated, LAYER_WEIGHTS, normalize="mass")
    A = gene_adjacency(fused, genes, layers=layers)
    return A, report


def train_arm(adata, A, tau: float, seed: int, n_latent: int, n_epochs: int,
              resolution: float, batch_key: str | None):
    """Train one arm and cluster it. ``tau=0`` gives the unpenalised control."""
    scvi.settings.seed = seed
    work = adata.copy()

    # Raw counts, not the normalised matrix: the ZINB likelihood is a count model.
    layer = "counts" if "counts" in work.layers else None
    SCVINetworkPrior.setup_anndata(work, layer=layer, batch_key=batch_key)

    model = SCVINetworkPrior(work, n_latent=n_latent, adjacency=A, tau=tau)
    model.train(n_epochs, enable_progress_bar=False)

    work.obsm["X_scVI"] = model.get_latent_representation()
    sc.pp.neighbors(work, use_rep="X_scVI", random_state=seed)
    sc.tl.leiden(work, key_added="leiden_scVI", resolution=resolution,
                 random_state=seed)

    metrics = compute_metrics(work)
    metrics["n_clusters"] = int(work.obs["leiden_scVI"].nunique())
    metrics["tau"] = tau
    metrics.update(theta_alignment(model.module))
    return work, metrics


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Run scPathWalk and its matched tau=0 control on one dataset.")
    ap.add_argument("--h5ad", required=True, type=Path,
                    help="Input AnnData file of raw counts.")
    ap.add_argument("--out", type=Path, default=Path("out"),
                    help="Directory for metrics, cluster assignments and markers.")
    ap.add_argument("--tau", type=float, default=30.0,
                    help="Penalty strength. 0 disables the network term entirely.")
    ap.add_argument("--n-latent", type=int, default=3)
    ap.add_argument("--n-epochs", type=int, default=30)
    ap.add_argument("--resolution", type=float, default=0.30,
                    help="Leiden resolution, shared by both arms.")
    ap.add_argument("--n-top-genes", type=int, default=4000)
    ap.add_argument("--seed", type=int, default=420)
    ap.add_argument("--batch-key", default=None,
                    help="Column in obs to pass as a batch covariate, e.g. sample.")
    ap.add_argument("--layers", nargs="*", default=list(CURATED_LAYERS),
                    choices=list(CURATED_LAYERS),
                    help="Curated layers to include in the adjacency.")
    ap.add_argument("--pathway-repr", default="bipartite",
                    choices=["bipartite", "clique"],
                    help="How to turn gene sets into edges. Bipartite contributes "
                         "no gene-gene edges; clique expands each set into all pairs.")
    ap.add_argument("--no-control", action="store_true",
                    help="Skip the tau=0 arm.")
    args = ap.parse_args()

    if not args.h5ad.exists():
        print(f"error: {args.h5ad} not found", file=sys.stderr)
        return 1

    args.out.mkdir(parents=True, exist_ok=True)
    layers = tuple(args.layers)

    print(f"reading {args.h5ad}")
    adata = process_adata(sc.read_h5ad(args.h5ad), n_top_genes=args.n_top_genes)
    print(f"  {adata.n_obs} cells x {adata.n_vars} genes after preprocessing")

    print("building the curated adjacency")
    A, report = build_adjacency(adata, layers, args.pathway_repr)
    summary = adjacency_summary(A, adata.var_names.tolist())
    print(report.to_string(index=False))
    print(f"  fitted adjacency: {summary['n_edges']} edges, "
          f"density {summary['density']:.4f}, "
          f"{summary['connected_frac']:.1%} of genes connected")
    if summary["n_edges"] == 0:
        print("error: the adjacency is empty, so the penalty cannot act. Stage the "
              "sources with 'uv run fetch_interactomes.py'.", file=sys.stderr)
        return 1

    rows = []
    arms = [("scpathwalk", args.tau)]
    if not args.no_control and args.tau != 0.0:
        arms.append(("control", 0.0))

    for name, tau in arms:
        print(f"\ntraining {name} (tau={tau:g})")
        work, metrics = train_arm(
            adata, A, tau, args.seed, args.n_latent, args.n_epochs,
            args.resolution, args.batch_key)
        metrics["arm"] = name
        rows.append(metrics)
        print(f"  silhouette {metrics['silhouette']:.4f}  "
              f"CH {metrics['calinski_harabasz']:.0f}  "
              f"DB {metrics['davies_bouldin']:.4f}  "
              f"k {metrics['n_clusters']}")

        stem = args.out / f"{args.h5ad.stem}_{name}"
        work.obs[["leiden_scVI"]].to_csv(f"{stem}_clusters.csv")
        np.save(f"{stem}_latent.npy", work.obsm["X_scVI"])
        sc.tl.rank_genes_groups(work, "leiden_scVI", method="wilcoxon")
        sc.get.rank_genes_groups_df(work, None).to_csv(
            f"{stem}_markers.csv", index=False)

    metrics_path = args.out / f"{args.h5ad.stem}_metrics.json"
    metrics_path.write_text(json.dumps(
        {"adjacency": summary, "config": {
            "tau": args.tau, "n_latent": args.n_latent, "n_epochs": args.n_epochs,
            "resolution": args.resolution, "seed": args.seed,
            "layers": list(layers), "pathway_repr": args.pathway_repr,
            "batch_key": args.batch_key, "n_top_genes": args.n_top_genes,
        }, "arms": rows}, indent=2, default=float), encoding="utf-8")

    if len(rows) == 2:
        delta = rows[0]["silhouette"] - rows[1]["silhouette"]
        print(f"\npaired silhouette difference: {delta:+.4f}")

    print(f"wrote {metrics_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
