# scPathWalk

Network-regularised gene dispersion for single-cell RNA-seq clustering.

A variational autoencoder over scRNA-seq counts learns one inverse dispersion per
gene, freely and independently. Most genes are detected in few cells, so those
estimates are poorly determined — and because dispersion controls how heavily a
gene's residuals steer training, a badly estimated value quietly distorts the
learned representation.

scPathWalk pools those estimates along curated gene relationships. A penalty on the
training objective pulls each gene's log-dispersion toward the weighted average of
its network neighbours, so a poorly measured gene borrows from its interaction
partners rather than from unrelated genes that happen to be equally expressed.

The count matrix is never modified: no imputation, smoothing, or reweighting. Only
the model's own assessment of per-gene reliability changes. Setting `tau=0` removes
the penalty exactly, which gives a matched control differing in one coefficient.

## Install

Requires Python 3.11+ and [uv](https://docs.astral.sh/uv/).

```bash
git clone https://github.com/JUCompBio/scPW.git
cd scPW
uv sync
```

`pyproject.toml` pins the CUDA 12.8 wheel channel for torch. Remove the
`[tool.uv.sources]` and `[[tool.uv.index]]` blocks for a CPU-only install; training
works on CPU, more slowly.

## Stage the interactome

Run once. Everything afterwards reads the local cache, so results cannot shift
because an upstream database was updated or an API was briefly unavailable.

```bash
uv run fetch_interactomes.py
```

This writes `data/interactome/*.parquet` plus a manifest recording each source's
row count and SHA-256. Requires network access.

| Source | Content | Role |
|---|---|---|
| STRING v12.0 | protein interactions, pairwise | training |
| CollecTRI | TF→target regulons, pairwise | training |
| KEGG 2021 | pathway membership, set-valued | training, see note |
| CORUM | complex membership, set-valued | training, see note |
| BioGRID | physical interactions | hold-out evaluation |
| Reactome 2022 | pathway membership | hold-out evaluation |
| GO BP 2023 | biological process terms | hold-out evaluation |

**Note on the set-valued sources.** STRING and CollecTRI name gene *pairs*, so their
records become network edges directly. KEGG and CORUM name *groups* — "this pathway
contains these 124 genes" asserts nothing about any particular pair. Under the
default `bipartite` representation each group becomes its own node with member genes
attached to it, so those layers contribute **no gene–gene edges** and do not enter
the penalty. The fitted network is STRING plus CollecTRI.

Pass `--pathway-repr clique` to expand each group into all pairs instead. That is a
modelling assumption, not an observation, and it makes a gene's influence scale with
the size of the groups it belongs to.

The hold-out sources are never read during training. `fetch_interactomes.py` reports
the Jaccard overlap between training and hold-out sources, since a hold-out that
merely restates a training source is not independent evidence.

## Run

Input is an AnnData file of raw counts.

```bash
uv run run_scpathwalk.py --h5ad data/my_sample.h5ad
```

This preprocesses, builds the adjacency, trains both the penalised arm and the
`tau=0` control, clusters each with Leiden, and writes metrics, cluster assignments,
latent representations and marker tables to `out/`.

```bash
# stronger penalty, two latent dimensions
uv run run_scpathwalk.py --h5ad data/my_sample.h5ad --tau 300 --n-latent 2

# multi-sample data: pass the sample column as a batch covariate
uv run run_scpathwalk.py --h5ad data/my_sample.h5ad --batch-key sample

# skip the control arm
uv run run_scpathwalk.py --h5ad data/my_sample.h5ad --no-control
```

Key arguments:

| Argument | Default | Meaning |
|---|---|---|
| `--tau` | 30 | penalty strength; 0 disables the network term |
| `--n-latent` | 3 | latent dimensionality |
| `--n-epochs` | 30 | training epochs |
| `--resolution` | 0.30 | Leiden resolution, shared by both arms |
| `--n-top-genes` | 4000 | highly variable genes retained |
| `--seed` | 420 | initialisation, posterior samples, Leiden tie-breaking |
| `--layers` | all four | curated layers to include |
| `--pathway-repr` | bipartite | how set-valued sources become edges |

`run_scpathwalk.py --help` lists the rest.

## Comparison methods

Seurat, Pagoda2 and scGSVA run in Docker containers, on the same QC-filtered matrix
the Python arm trains on. Exporting that matrix is what makes the comparison
meaningful: if each method applied its own filtering, differences would mix the
clustering method with the preprocessing.

```bash
# 1. export the same cells and genes the Python arm uses
uv run export_for_r.py --h5ad data/my_sample.h5ad

# 2. build the scGSVA image once (no official image exists)
docker build -t scpw/scgsva:pinned -f docker/scgsva.Dockerfile .

# 3. run a method
uv run run_comparators.py --method seurat --datasets my_sample
uv run run_comparators.py --method pagoda2 --datasets my_sample
uv run run_comparators.py --method scgsva  --datasets my_sample
```

Add `--dry-run` to see what would run. `--workers` and `--threads` bound container
parallelism; their product should stay within the physical core count, or the
containers contend for cores and run slower than serial. `--memory` caps each
container so one large dataset cannot exhaust the host.

Each method is scored in the representation it clusters on, as its own authors
evaluate it. Those spaces differ in dimensionality, so cross-method index values
describe the spaces each method delivers rather than being strictly comparable. The
`tau=0` comparison shares the representation family and the clustering procedure.

Both euclidean and cosine Silhouette are written out. The same partition in the same
embedding can differ by more than an order of magnitude between the two, and methods
differ in which distance they use internally, so reporting one would privilege
whichever method happened to match it.

## Layout

```
run_scpathwalk.py       entry point: train, cluster, score
fetch_interactomes.py   stage the external sources into data/interactome/
export_for_r.py         export a QC-filtered matrix to MTX for the R methods
run_comparators.py      dispatch the R methods in Docker containers
util.py                 shared preprocessing (process_adata)
scpathwalk/
  network_vae.py        the model and the dispersion penalty
  dispersion.py         build the gene adjacency, apply the layer filter
  graphs.py             layer construction and mass-normalised fusion
  interactome.py        source specs, download, caching, contamination check
  embed.py              latent representation, UMAP, Leiden
  evaluate.py           Silhouette, Calinski-Harabasz, Davies-Bouldin
  coherence.py          marker coherence against the hold-out sources
  baselines.py          PCA and Scanpy comparison pipelines
  visualize.py          UMAP plots, marker rankings, pathway overlays
  logging_setup.py      structured console and file logging
  logging_utils.py      run records and environment capture
r/                      the R clustering scripts and shared metric definitions
docker/                 Dockerfile for the scGSVA image
```

## Notes

- `fetch_interactomes.py` needs network access. Everything afterwards reads the cache.
- The ZINB likelihood is a count model, so `run_scpathwalk.py` registers
  `layers["counts"]` when present. `util.process_adata` preserves raw counts there
  before normalising `X`.
- Genes with no neighbours in the network receive no penalty term and keep the
  estimate the model produces on its own.
- Cluster-quality indices respond to the number of clusters as well as to structure,
  so cluster counts are reported alongside every metric.

## Licence

MIT. See `LICENSE`.
