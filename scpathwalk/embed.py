"""
scVI embedding, Leiden clustering, and UMAP computation.
"""

import anndata as ad
import scanpy as sc
import scvi


def setup_and_train(
    adata: ad.AnnData,
    n_latent: int = 5,
    n_epochs: int = 30,
    seed: int = 420,
    global_seed: int = 69,
    progress_bar: bool = False,
) -> scvi.model.SCVI:
    """Configure seeds, set up scVI, train, and return the model.

    ``progress_bar`` is off by default: across a benchmarking grid the Lightning
    bar emits hundreds of lines per dataset and obscures the structured log.
    """
    scvi.settings.seed = seed
    scvi.seed = global_seed
    scvi.settings.progress_bar_style = "tqdm"

    scvi.model.SCVI.setup_anndata(adata)
    model = scvi.model.SCVI(adata, n_latent=n_latent)
    model.train(n_epochs, enable_progress_bar=progress_bar)
    return model


def embed_and_cluster(
    adata: ad.AnnData,
    model: scvi.model.SCVI,
    resolution: float = 0.55,
    min_dist: float = 0.3,
    seed: int = 420,
) -> ad.AnnData:
    """Get latent representation, build UMAP, run Leiden clustering.

    ``seed`` is threaded into both UMAP and Leiden so a run is reproducible and so
    seed-to-seed variance can be measured.
    """
    latent_key = "X_scVI"
    cluster_key = "leiden_scVI"

    adata.obsm[latent_key] = model.get_latent_representation()
    sc.pp.neighbors(adata, use_rep=latent_key, random_state=seed)
    sc.tl.umap(adata, min_dist=min_dist, random_state=seed)
    sc.tl.leiden(adata, key_added=cluster_key, resolution=resolution, random_state=seed)

    return adata
