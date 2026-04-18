# scPathWalk

**scPathWalk: Improving Cell Subpopulation Discovery through Interactome Integration**

A pathway-guided single-cell RNA-seq clustering framework that integrates gene co-expression, protein-protein interactions (STRING DB), and KEGG pathway annotations into a unified network. Gene importance is estimated via Random Walk with Restart (RWR) and used to reweight the expression matrix before variational autoencoder embedding (scVI) and Leiden clustering.

## Requirements

- Python >= 3.10
- [uv](https://docs.astral.sh/uv/) (recommended package manager)
- Internet access (STRING DB and Enrichr API calls)
- GPU recommended for scVI training (CPU works but is slower)

## Setup

```bash
git clone https://github.com/JUCompBio/scPW.git
cd scPW
uv sync
```

## Usage

### Step 1: Download preprocessed data

Preprocessed `.h5ad` files for all seven datasets are available as [release assets](https://github.com/JUCompBio/scPW/releases/tag/v0.1.0):

```bash
mkdir processed
gh release download v0.1.0 --repo JUCompBio/scPW --dir processed
```

Alternatively, to preprocess from raw 10x Genomics files yourself, place them under `raw/<dataset>/` and run:

```bash
uv run process_data.py --dataset ovary
```

### Step 2: Run the scPathWalk pipeline

```bash
uv run run_pipeline.py --dataset ovary
```

Add `--no-baseline` to skip the vanilla scVI comparison.

## Pipeline Overview

1. **Preprocessing** -- Normalize, log-transform, Scrublet doublet detection, Seurat v3 HVG selection
2. **Correlation graph** -- Pearson correlation between genes (threshold tau=0.01)
3. **PPI graph** -- STRING DB protein-protein interactions (score >= 400)
4. **Pathway graph** -- KEGG pathway cliques via Enrichr
5. **Consensus partition** -- Weighted Leiden consensus across the three graphs
6. **PageRank (RWR)** -- Gene importance scores on the merged union graph (alpha=0.9)
7. **Reweight** -- MinMaxScale + expression rescaling by RWR scores
8. **Embed + Cluster** -- scVI (latent dim=5, 30 epochs) + Leiden clustering
9. **Evaluate** -- Silhouette, Calinski-Harabasz, Davies-Bouldin metrics
10. **Visualize** -- UMAP plots, ranked gene dot plots, pathway annotation overlays

## Datasets

| Dataset   | Tissue / Disease         | Source       |
|-----------|--------------------------|--------------|
| 3praw     | Mammary gland            | 10x Genomics |
| 4plex     | Breast tissue            | 10x Genomics |
| Cervical  | Cervical cancer          | 10x Genomics |
| Ovary     | Ovarian cancer           | 10x Genomics |
| GSE161529 | Breast cancer atlas      | GEO          |
| GSE228499 | HR+/HER2- breast cancer  | GEO          |
| GSE266351 | TNBC (targeted panel)    | GEO          |

## Project Structure

```
scPW/
  process_data.py          # Step 1: raw 10x -> .h5ad
  run_pipeline.py          # Step 2: full scPathWalk pipeline
  util.py                  # Shared preprocessing (normalize, HVG)
  scpathwalk/
    __init__.py
    graphs.py              # Correlation, PPI, pathway graph construction
    partition.py           # Leiden partitioning, consensus clustering
    reweight.py            # PageRank + expression reweighting
    embed.py               # scVI training, UMAP, Leiden clustering
    evaluate.py            # Clustering quality metrics
    visualize.py           # UMAP plots, gene ranking, pathway overlay
```

## Outputs

- `figures/` -- UMAP and dot plot PDFs
- `rank/` -- Per-cluster ranked gene CSVs

## Citation

If you use scPathWalk in your research, please cite:

```bibtex
@article{dey2026scpathwalk,
  title={scPathWalk: Improving Cell Subpopulation Discovery through Interactome Integration},
  author={Dey, Ashmita and Das, Rangan and Maulik, Ujjwal and Bandyopadhyay, Sanghamitra},
  journal={IEEE Access},
  year={2026}
}
```

## License

MIT
