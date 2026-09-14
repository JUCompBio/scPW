#%%
import numpy as np
import scanpy as sc
import anndata as ad

#%%
def _looks_preprocessed(adata, n_top_genes):
    """Has this object already been through process_adata?

    Two independent signs, both required, so a genuinely raw small panel is not
    mistaken for a processed object:

    * ``X`` is non-integer, i.e. already normalised and log-transformed;
    * the gene count is at or below the HVG target, i.e. already subset.

    This matters because a caller that saves a processed ``.h5ad`` and then passes
    it back through would otherwise re-normalise already-normalised values and,
    worse, re-apply ``min_genes=100``. That threshold is calibrated for a full
    transcriptome of ~20,000 genes; on a 4,000-gene HVG subset most genes are zero
    in any given cell, so a second pass discards the majority of cells -- not poor
    quality cells, just cells whose HVGs happen to be sparse.
    """
    X = adata.X
    sample = X[:200].toarray() if hasattr(X, "toarray") else np.asarray(X[:200])
    non_integer = not np.allclose(sample, np.round(sample))
    already_subset = adata.n_vars <= n_top_genes
    return bool(non_integer and already_subset)


def process_adata(raw_adata, n_top_genes=4000, filter=True, force=False):
    """Normalise, log-transform, QC-filter and HVG-subset a raw AnnData.

    Idempotent: an object that has already been through this is returned
    unchanged (see :func:`_looks_preprocessed`). Pass ``force=True`` to override.
    """
    adata = raw_adata.copy()
    adata.obs_names_make_unique()

    if not force and _looks_preprocessed(adata, n_top_genes):
        # Already processed. Re-running would double-normalise and re-apply a
        # cell filter calibrated for the full transcriptome to an HVG subset.
        if "counts" not in adata.layers:
            adata.layers["counts"] = adata.X.copy()
        return adata

    adata.layers["counts"] = adata.X.copy()  # preserve counts before transforming
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata  # freeze the state in `.raw`

    if filter:
        sc.pp.filter_cells(adata, min_genes=100)
        sc.pp.filter_genes(adata, min_cells=3)

    try:  # with a batch key
        sc.pp.scrublet(adata, batch_key="sample")
        sc.pp.highly_variable_genes(
            adata,
            n_top_genes=n_top_genes,
            subset=True,
            layer="counts",
            flavor="seurat_v3",
            batch_key="sample"
            )
    except Exception:  # else without the batch key
        sc.pp.highly_variable_genes(
            adata,
            n_top_genes=n_top_genes,
            subset=True,
            layer="counts",
            flavor="seurat_v3",
            )

    return adata

#%%
def combine_pex(adata, pex_label = 'protein_expression'):
    pex = adata.obsm[pex_label]
    try:
        pex_adata = ad.AnnData(X=pex, obs=adata.obs, dtype=pex.dtype)
    except:
        pex = pex.to_numpy()
        pex_adata = ad.AnnData(X=pex, obs=adata.obs, dtype=pex.dtype)
    sc.pp.normalize_total(pex_adata, target_sum=1e4)
    sc.pp.log1p(pex_adata)
    full_adata = ad.concat([adata, pex_adata], axis=1)
    full_adata.obs = adata.obs
    full_adata.uns = adata.uns
    full_adata.obsm = adata.obsm
    full_adata.var_names_make_unique()

    pex_idx = pex_adata.X.shape[1]
    return (full_adata, pex_idx)