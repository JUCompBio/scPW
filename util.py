import scanpy as sc


def process_adata(raw_adata, n_top_genes=4000, filter=True):
    """Preprocess AnnData: normalize, log1p, Scrublet, Seurat v3 HVG selection."""
    adata = raw_adata.copy()
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.obs_names_make_unique()
    adata.raw = adata

    if filter:
        sc.pp.filter_cells(adata, min_genes=100)
        sc.pp.filter_genes(adata, min_cells=3)

    try:
        sc.pp.scrublet(adata, batch_key="sample")
        sc.pp.highly_variable_genes(
            adata,
            n_top_genes=n_top_genes,
            subset=True,
            layer="counts",
            flavor="seurat_v3",
            batch_key="sample",
        )
    except Exception:
        sc.pp.highly_variable_genes(
            adata,
            n_top_genes=n_top_genes,
            subset=True,
            layer="counts",
            flavor="seurat_v3",
        )

    return adata
