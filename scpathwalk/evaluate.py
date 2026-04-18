"""Clustering evaluation metrics."""

from sklearn.metrics import (
    silhouette_score,
    calinski_harabasz_score,
    davies_bouldin_score,
)


def compute_metrics(
    adata, latent_key: str = "X_scVI", cluster_key: str = "leiden_scVI"
) -> dict:
    """Compute Silhouette, Calinski-Harabasz, and Davies-Bouldin scores."""
    X = adata.obsm[latent_key]
    labels = adata.obs[cluster_key]

    return {
        "silhouette": silhouette_score(X, labels),
        "calinski_harabasz": calinski_harabasz_score(X, labels),
        "davies_bouldin": davies_bouldin_score(X, labels),
    }
