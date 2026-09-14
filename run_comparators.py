"""
Run the R comparison methods in parallel Docker containers.

Three methods, each in its own image, each reading the same QC-filtered matrices
that the Python arm trains on. Export those first with ``export_for_r.py``, so a
difference between arms cannot come from different preprocessing.

Two things have to be constrained or parallel containers are slower than serial.

**Threads.** R inside the container sees every logical CPU and ``OMP_NUM_THREADS``
is unset, so its BLAS spawns a thread per core by default. Four such containers
would oversubscribe the machine several times over and spend most of their time
context-switching. Each container is therefore given ``--cpus`` and a matching
``OMP_NUM_THREADS``, so workers x threads stays within the physical core count.

**Memory.** Large datasets make Seurat's ``FindAllMarkers`` and a full distance
matrix both memory-hungry. Containers get an explicit ``--memory`` limit so one
dataset cannot exhaust the host and take the other workers with it; a container
that exceeds its limit is killed with a clear OOM rather than degrading everything.

Datasets are dispatched largest-first, because the tail of a parallel run is set by
the longest single job.

Usage:
    uv run run_comparators.py --datasets my_sample --dry-run
    uv run run_comparators.py --method seurat --datasets my_sample other_sample
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parent

IMAGES = {
    "seurat": "satijalab/seurat:5.0.0",
    "pagoda2": "pkharchenkolab/pagoda2:latest",
    # No official scGSVA image exists. Build it once from docker/scgsva.Dockerfile:
    #   docker build -t scpw/scgsva:pinned -f docker/scgsva.Dockerfile .
    # Installing scGSVA from GitHub inside every run would make each benchmark
    # depend on the network and on whatever the default branch happens to be.
    "scgsva": "scpw/scgsva:pinned",
}

SCRIPTS = {
    "seurat": "r/seurat_cluster.R",
    "pagoda2": "r/pagoda2_cluster.R",
    "scgsva": "r/scgsva_cluster.R",
}

# Windows Docker Desktop does not always put docker on PATH.
WINDOWS_DOCKER = r"C:\Program Files\Docker\Docker\resources\bin\docker.exe"


def docker_bin() -> str:
    """Resolve the docker executable, falling back to the Windows install path."""
    found = shutil.which("docker")
    if found:
        return found
    if Path(WINDOWS_DOCKER).exists():
        return WINDOWS_DOCKER
    return "docker"


def run_one(method: str, dataset: str, seed: int, resolution: float,
            threads: int, memory: str, in_root: str, out_root: str
            ) -> tuple[str, bool, str]:
    """One container for one (method, dataset). Returns (dataset, ok, message)."""
    cmd = [
        docker_bin(), "run", "--rm",
        "-v", f"{ROOT}:/work",
        "-w", "/work",
        # Hard CPU and memory ceilings. Without --cpus the container's BLAS sees
        # every host core regardless of how many containers are running.
        "--cpus", str(threads),
        "--memory", memory,
        "-e", f"OMP_NUM_THREADS={threads}",
        "-e", f"MKL_NUM_THREADS={threads}",
        "-e", f"OPENBLAS_NUM_THREADS={threads}",
        IMAGES[method],
        "Rscript", f"/work/{SCRIPTS[method]}",
        f"/work/{in_root}/{dataset}", f"/work/{out_root}/{method}",
        dataset, str(seed), str(resolution), str(threads),
    ]
    t0 = time.time()
    p = subprocess.run(cmd, capture_output=True, text=True)
    dt = time.time() - t0

    if p.returncode != 0:
        tail = (p.stderr or p.stdout or "").strip().splitlines()[-3:]
        return dataset, False, f"exit {p.returncode} after {dt:.0f}s: {' | '.join(tail)}"
    done = [line for line in p.stdout.splitlines() if line.startswith("DONE")]
    return dataset, True, (done[-1] if done else f"completed in {dt:.0f}s") + f" [{dt:.0f}s]"


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Run the R comparison methods in Docker containers.")
    ap.add_argument("--method", default="seurat", choices=sorted(IMAGES))
    ap.add_argument("--datasets", nargs="+", required=True,
                    help="Dataset names, matching the export directories.")
    ap.add_argument("--in-root", default="data/r_export",
                    help="Directory holding one subdirectory per exported dataset.")
    ap.add_argument("--out-root", default="out",
                    help="Results are written to <out-root>/<method>/.")
    ap.add_argument("--seed", type=int, default=420)
    ap.add_argument("--resolution", type=float, default=0.55,
                    help="Clustering resolution. Pagoda2 ignores this; it exposes "
                         "no resolution parameter and uses its own default.")
    ap.add_argument("--workers", type=int, default=4,
                    help="Containers in parallel.")
    ap.add_argument("--threads", type=int, default=2,
                    help="CPUs per container. workers x threads should not exceed "
                         "the physical core count, or the containers contend.")
    ap.add_argument("--memory", default="12g",
                    help="Per-container limit. Large datasets need headroom for "
                         "FindAllMarkers and the distance matrix.")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    missing = [d for d in args.datasets
               if not (ROOT / args.in_root / d / "matrix.mtx").exists()]
    if missing:
        print(f"no exported matrices for: {', '.join(missing)}", file=sys.stderr)
        print(f"run 'uv run export_for_r.py --h5ad <file>' first, so R sees the "
              f"same QC-filtered cells and genes as the Python arm", file=sys.stderr)
        return 1

    # Largest first: the tail of a parallel run is set by the longest single job.
    order = sorted(args.datasets, reverse=True, key=lambda d: (
        (ROOT / args.in_root / d / "matrix.mtx").stat().st_size))

    print(f"{args.method}: {len(order)} datasets, {args.workers} containers "
          f"x {args.threads} cpus, {args.memory} each")
    print(f"image {IMAGES[args.method]}, seed {args.seed}, "
          f"resolution {args.resolution}")
    print(f"output -> {args.out_root}/{args.method}/")
    if args.dry_run:
        for d in order:
            print(f"  {d}")
        return 0

    (ROOT / args.out_root / args.method).mkdir(parents=True, exist_ok=True)
    started = time.time()
    ok = 0
    # Threads, not processes: each task only waits on a subprocess.
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futs = {pool.submit(run_one, args.method, d, args.seed, args.resolution,
                            args.threads, args.memory, args.in_root,
                            args.out_root): d for d in order}
        for f in as_completed(futs):
            dataset, good, msg = f.result()
            ok += int(good)
            flag = "" if good else "FAILED  "
            print(f"  {flag}{dataset:<20} {msg}", flush=True)

    print(f"\n{ok}/{len(order)} succeeded in {(time.time() - started) / 60:.1f} min")
    return 0 if ok == len(order) else 1


if __name__ == "__main__":
    raise SystemExit(main())
