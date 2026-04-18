"""scVI embedding, Leiden clustering, and UMAP computation."""

import anndata as ad
import scanpy as sc
import scvi


def setup_and_train(
    adata: ad.AnnData,
    n_latent: int = 5,
    n_epochs: int = 30,
    seed: int = 420,
    global_seed: int = 69,
) -> scvi.model.SCVI:
    """Configure seeds, set up scVI, train, and return the model."""
    scvi.settings.seed = seed
    scvi.seed = global_seed

    scvi.model.SCVI.setup_anndata(adata)
    model = scvi.model.SCVI(adata, n_latent=n_latent)
    model.train(n_epochs)
    return model


def embed_and_cluster(
    adata: ad.AnnData,
    model: scvi.model.SCVI,
    resolution: float = 0.55,
    min_dist: float = 0.3,
) -> ad.AnnData:
    """Get latent representation, build UMAP, run Leiden clustering."""
    latent_key = "X_scVI"
    cluster_key = "leiden_scVI"

    adata.obsm[latent_key] = model.get_latent_representation()
    sc.pp.neighbors(adata, use_rep=latent_key)
    sc.tl.umap(adata, min_dist=min_dist)
    sc.tl.leiden(adata, key_added=cluster_key, resolution=resolution)

    return adata
