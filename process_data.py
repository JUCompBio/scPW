"""
process_data.py -- Convert raw 10x Genomics data to preprocessed .h5ad files.

Usage:
    uv run process_data.py                    # Process all datasets
    uv run process_data.py --dataset 3praw    # Process a single dataset
"""

import argparse
import glob
import warnings
from pathlib import Path

import scanpy as sc
from tqdm.auto import tqdm

from util import process_adata

warnings.filterwarnings("ignore")

SINGLE_DATASETS = ["3praw", "4plex", "9k_cervical_cancer", "GSE161529", "ovary"]
MULTI_DATASETS = ["GSE228499", "GSE266351"]

RAW_DIR = Path("raw")
PROCESSED_DIR = Path("processed")


def process_single(name: str) -> None:
    """Load a single-sample 10x dataset, preprocess, and save."""
    print(f"\n{'='*60}\nProcessing {name}\n{'='*60}")
    adata = sc.read_10x_mtx(RAW_DIR / name)
    print(adata)
    adata_pp = process_adata(adata)
    out = PROCESSED_DIR / f"{name}.h5ad"
    adata_pp.write_h5ad(out, compression="gzip")
    print(f"Saved -> {out}  ({adata_pp.n_obs} cells x {adata_pp.n_vars} genes)")


def process_multi(name: str) -> None:
    """Load a multi-sample 10x dataset, merge batches, preprocess, and save."""
    print(f"\n{'='*60}\nProcessing {name} (multi-sample)\n{'='*60}")

    files = glob.glob(str(RAW_DIR / name / "*.gz"))
    file_names = [Path(f).name for f in files]
    prefixes = sorted(set(fn.rsplit("_", 1)[0] + "_" for fn in file_names))

    print(f"Found {len(prefixes)} sample prefixes")

    adata = sc.read_10x_mtx(RAW_DIR / name, prefix=prefixes[0])
    print(f"  Loaded {prefixes[0]}: {adata.n_obs} cells")

    for pfx in tqdm(prefixes[1:], desc="Loading batches"):
        batch = sc.read_10x_mtx(RAW_DIR / name, prefix=pfx)
        print(f"  Loaded {pfx}: {batch.n_obs} cells")
        adata = adata.concatenate(batch, batch_key="sample")

    print(f"Merged: {adata.n_obs} cells x {adata.n_vars} genes")

    adata_pp = process_adata(adata)
    out = PROCESSED_DIR / f"{name}.h5ad"
    adata_pp.write_h5ad(out, compression="gzip")
    print(f"Saved -> {out}  ({adata_pp.n_obs} cells x {adata_pp.n_vars} genes)")


def main():
    parser = argparse.ArgumentParser(description="Preprocess raw 10x data -> .h5ad")
    parser.add_argument(
        "--dataset",
        type=str,
        default=None,
        choices=SINGLE_DATASETS + MULTI_DATASETS,
        help="Process a single dataset (default: all)",
    )
    args = parser.parse_args()

    PROCESSED_DIR.mkdir(exist_ok=True)

    if args.dataset:
        datasets = [args.dataset]
    else:
        datasets = SINGLE_DATASETS + MULTI_DATASETS

    for name in datasets:
        if name in MULTI_DATASETS:
            process_multi(name)
        else:
            process_single(name)

    print("\nDone.")


if __name__ == "__main__":
    main()
