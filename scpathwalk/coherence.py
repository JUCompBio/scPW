"""
Biological coherence of cluster marker sets, measured on hold-out resources.

The geometric metrics (Silhouette, Calinski-Harabasz, Davies-Bouldin) ask whether
clusters are compact and separated in the embedding. They do not ask whether the
clusters correspond to anything biological. scPathWalk pools dispersion by
interactome centrality, so what it optimises is functional coherence of the genes
that define each cluster, and that is what this module measures.

Two families of measure, both against resources the method never consumes:

**Marker Network Coherence (MNC).** For each cluster, the density of the subgraph
induced by its top-k markers in a hold-out interactome (BioGRID physical),
z-scored against random gene sets of the same size.

**Hold-out pathway enrichment.** Per-cluster over-representation in Reactome and
GO Biological Process. Less novel than MNC and correspondingly easier for a
reader to accept, since per-cluster enrichment is standard practice.

A caution that shapes the design. MNC rewards markers being interactome-connected,
and scPathWalk up-weights interactome-central genes, so there is a route by which
it could win for a partly circular reason: hub proteins are hubs in *both* STRING
and BioGRID, because both databases reflect which proteins have been studied most.
``degree_controlled_mnc`` addresses this by comparing against a null that matches
the observed marker set on BioGRID degree, so a gain that is purely "picked
better-connected genes" does not register as coherence. Report both.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

from scpathwalk.interactome import load_source

# Databases that must never be used to build a training layer.
HOLDOUT_NETWORK = "biogrid_physical"
HOLDOUT_PATHWAYS = ("reactome_pathways", "gobp_pathways")


# ── Network coherence ─────────────────────────────────────────────────────────
def load_holdout_adjacency(source: str = HOLDOUT_NETWORK) -> dict[str, set[str]]:
    """Load a hold-out interactome as an adjacency mapping keyed by gene symbol."""
    edges = load_source(source)
    adjacency: dict[str, set[str]] = {}
    for a, b in zip(edges["key_a"].to_numpy(), edges["key_b"].to_numpy()):
        adjacency.setdefault(a, set()).add(b)
        adjacency.setdefault(b, set()).add(a)
    return adjacency


def induced_subgraph_density(adjacency: dict[str, set[str]], genes) -> float:
    """Density of the subgraph induced by ``genes``: edges present / possible.

    Returns 0.0 for fewer than two genes present in the network.
    """
    present = [g for g in dict.fromkeys(genes) if g in adjacency]
    n = len(present)
    if n < 2:
        return 0.0
    member = set(present)
    edges = sum(len(adjacency[g] & member) for g in present) / 2.0
    return float(edges / (n * (n - 1) / 2.0))


def _null_densities(
    adjacency: dict[str, set[str]],
    universe: list[str],
    size: int,
    n_perm: int,
    rng: np.random.Generator,
) -> np.ndarray:
    """Densities of ``n_perm`` random gene sets of the given size."""
    universe_array = np.asarray(universe)
    if size < 2 or len(universe_array) < size:
        return np.zeros(n_perm)
    return np.array([
        induced_subgraph_density(adjacency, rng.choice(universe_array, size, replace=False))
        for _ in range(n_perm)
    ])


def _degree_matched_null(
    adjacency: dict[str, set[str]],
    universe: list[str],
    observed: list[str],
    n_perm: int,
    rng: np.random.Generator,
    n_bins: int = 10,
) -> np.ndarray:
    """Null that matches the observed set's BioGRID degree distribution.

    Genes are binned by degree; each draw samples the same number from each bin as
    the observed set contains. A method that merely selects better-connected genes
    scores near zero against this null, so a positive result here is evidence of
    coherence beyond hubness.
    """
    degrees = {g: len(adjacency.get(g, ())) for g in universe}
    if not degrees:
        return np.zeros(n_perm)

    values = np.array([degrees[g] for g in universe], dtype=float)
    # Rank-based bin edges: degree is heavy-tailed, so equal-width bins would put
    # almost every gene in the first bin.
    quantiles = np.quantile(values, np.linspace(0, 1, n_bins + 1))
    quantiles[-1] += 1e-9
    bin_of = np.digitize(values, quantiles[1:-1])
    pools: dict[int, np.ndarray] = {
        b: np.asarray([g for g, bb in zip(universe, bin_of) if bb == b])
        for b in range(n_bins)
    }

    observed_present = [g for g in observed if g in degrees]
    if len(observed_present) < 2:
        return np.zeros(n_perm)
    obs_bins = np.digitize(
        np.array([degrees[g] for g in observed_present], dtype=float), quantiles[1:-1]
    )
    wanted = pd.Series(obs_bins).value_counts().to_dict()

    out = []
    for _ in range(n_perm):
        sample: list[str] = []
        for b, count in wanted.items():
            pool = pools.get(int(b), np.array([]))
            if len(pool) == 0:
                continue
            take = min(count, len(pool))
            sample.extend(rng.choice(pool, take, replace=False).tolist())
        out.append(induced_subgraph_density(adjacency, sample))
    return np.asarray(out)


def mnc_per_cluster(
    markers: pd.DataFrame,
    universe: list[str],
    adjacency: dict[str, set[str]] | None = None,
    top_k: int = 50,
    n_perm: int = 1000,
    seed: int = 0,
    degree_control: bool = True,
) -> pd.DataFrame:
    """Marker Network Coherence per cluster.

    Args:
        markers: long table with ``names`` (gene) and ``label`` (cluster) columns,
            ordered so that the most significant markers come first, as written by
            ``visualize.rank_and_save_genes``.
        universe: genes the clustering could have chosen from (the HVG set). The
            null must be drawn from this, not from all of BioGRID, or the
            comparison is against a different gene population.
        adjacency: hold-out network; loaded from cache when omitted.
        top_k: markers per cluster.
        n_perm: permutations for each null.
        seed: RNG seed.
        degree_control: also compute the degree-matched z-score.

    Returns:
        One row per cluster with the observed density, both z-scores, and an
        empirical right-tail p-value against the size-matched null.
    """
    if adjacency is None:
        adjacency = load_holdout_adjacency()

    rng = np.random.default_rng(seed)
    universe_upper = [g.upper() for g in dict.fromkeys(universe)]
    in_network = [g for g in universe_upper if g in adjacency]

    rows = []
    null_cache: dict[int, np.ndarray] = {}

    for label, group in markers.groupby("label", sort=True):
        genes = [str(g).upper() for g in group["names"].head(top_k)]
        present = [g for g in dict.fromkeys(genes) if g in adjacency]
        observed = induced_subgraph_density(adjacency, present)

        size = len(present)
        if size not in null_cache:
            null_cache[size] = _null_densities(
                adjacency, in_network, size, n_perm, rng
            )
        null = null_cache[size]
        sd = null.std(ddof=1) if null.size > 1 else 0.0
        z = float((observed - null.mean()) / sd) if sd > 0 else 0.0
        pvalue = float((null >= observed).mean()) if null.size else 1.0

        record = {
            "cluster": str(label),
            "n_markers": len(genes),
            "n_in_network": size,
            "density": observed,
            "null_mean": float(null.mean()) if null.size else 0.0,
            "null_sd": float(sd),
            "mnc_z": z,
            "mnc_p": pvalue,
            "mean_degree": float(
                np.mean([len(adjacency[g]) for g in present]) if present else 0.0
            ),
        }

        if degree_control:
            dnull = _degree_matched_null(adjacency, in_network, present, n_perm, rng)
            dsd = dnull.std(ddof=1) if dnull.size > 1 else 0.0
            record["mnc_z_degree_matched"] = (
                float((observed - dnull.mean()) / dsd) if dsd > 0 else 0.0
            )
            record["degree_null_mean"] = float(dnull.mean()) if dnull.size else 0.0

        rows.append(record)

    return pd.DataFrame(rows)


def aggregate_mnc(per_cluster: pd.DataFrame) -> dict[str, float]:
    """Dataset-level summary. Median across clusters, which is robust to the one
    or two very small clusters most partitions contain."""
    if per_cluster.empty:
        return {}
    out = {
        "mnc_median": float(per_cluster["mnc_z"].median()),
        "mnc_mean": float(per_cluster["mnc_z"].mean()),
        "n_clusters": int(len(per_cluster)),
        "frac_significant": float((per_cluster["mnc_p"] < 0.05).mean()),
        "mean_degree": float(per_cluster["mean_degree"].mean()),
    }
    if "mnc_z_degree_matched" in per_cluster:
        out["mnc_median_degree_matched"] = float(
            per_cluster["mnc_z_degree_matched"].median()
        )
    return out


# ── Hold-out pathway enrichment ───────────────────────────────────────────────
def enrichment_per_cluster(
    markers: pd.DataFrame,
    universe: list[str],
    source: str,
    top_k: int = 50,
    min_set_size: int = 5,
    max_set_size: int = 500,
) -> pd.DataFrame:
    """Per-cluster over-representation against a hold-out gene-set collection.

    One-sided Fisher exact test per (cluster, gene set), Benjamini-Hochberg
    corrected within each cluster. Reports the strongest term per cluster.

    The universe is restricted to the HVG set so the test asks "enriched among the
    genes this clustering could have selected", not "enriched among all genes",
    which would inflate significance for any HVG-derived marker list.
    """
    sets = load_source(source)
    universe_upper = {g.upper() for g in universe}
    sets = sets[sets["gene"].isin(universe_upper)]

    id_col = "complex" if "complex" in sets.columns else "term"
    grouped = sets.groupby(id_col)["gene"].apply(set)
    grouped = grouped[grouped.map(len).between(min_set_size, max_set_size)]
    if grouped.empty:
        return pd.DataFrame()

    n_universe = len(universe_upper)
    rows = []

    for label, group in markers.groupby("label", sort=True):
        genes = {str(g).upper() for g in group["names"].head(top_k)}
        genes &= universe_upper
        if len(genes) < 2:
            continue

        records = []
        for term, members in grouped.items():
            overlap = len(genes & members)
            if overlap == 0:
                continue
            table = [
                [overlap, len(genes) - overlap],
                [len(members) - overlap,
                 n_universe - len(genes) - len(members) + overlap],
            ]
            _, pvalue = stats.fisher_exact(table, alternative="greater")
            records.append((term, overlap, len(members), pvalue))

        if not records:
            continue

        frame = pd.DataFrame(records, columns=["term", "overlap", "set_size", "pval"])
        frame = frame.sort_values("pval", ignore_index=True)
        m = len(frame)
        frame["padj"] = np.minimum.accumulate(
            (frame["pval"] * m / (frame.index.to_numpy() + 1))[::-1]
        )[::-1].clip(0, 1)

        best = frame.iloc[0]
        rows.append({
            "cluster": str(label),
            "n_markers": len(genes),
            "top_term": best["term"],
            "overlap": int(best["overlap"]),
            "set_size": int(best["set_size"]),
            "padj": float(best["padj"]),
            "neg_log10_padj": float(-np.log10(max(best["padj"], 1e-300))),
            "n_significant": int((frame["padj"] < 0.05).sum()),
        })

    return pd.DataFrame(rows)


def aggregate_enrichment(per_cluster: pd.DataFrame) -> dict[str, float]:
    """Dataset-level enrichment summary."""
    if per_cluster.empty:
        return {}
    return {
        "median_neg_log10_padj": float(per_cluster["neg_log10_padj"].median()),
        "frac_clusters_enriched": float((per_cluster["padj"] < 0.05).mean()),
        "mean_n_significant": float(per_cluster["n_significant"].mean()),
        "n_clusters": int(len(per_cluster)),
    }


def pathway_uniqueness(
    markers: pd.DataFrame,
    universe: list[str],
    source: str,
    top_k: int = 50,
) -> float:
    """Mean pairwise Jaccard *distance* between clusters' top marker sets.

    High values mean clusters are characterised by distinct gene programmes rather
    than by overlapping variations on the same one, which is the question behind
    asking whether reported markers are uniquely enriched in specific clusters.
    """
    universe_upper = {g.upper() for g in universe}
    per_cluster = [
        {str(g).upper() for g in group["names"].head(top_k)} & universe_upper
        for _, group in markers.groupby("label", sort=True)
    ]
    per_cluster = [s for s in per_cluster if s]
    if len(per_cluster) < 2:
        return float("nan")

    distances = []
    for i in range(len(per_cluster)):
        for j in range(i + 1, len(per_cluster)):
            a, b = per_cluster[i], per_cluster[j]
            union = len(a | b)
            distances.append(1.0 - (len(a & b) / union if union else 0.0))
    return float(np.mean(distances))
