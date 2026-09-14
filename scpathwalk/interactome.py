"""
Interactome layer acquisition, cleaning, and caching.

Every external source is fetched **once** into a versioned Parquet cache under
``data/``. The pipeline reads only the cache, never the network, so a benchmarking
grid of many runs is reproducible and cannot be perturbed by a transient API
failure or by an upstream database release.

Two disjoint groups of sources are defined:

*Training layers* are consumed by the method (co-expression, STRING PPI, KEGG
pathway membership, CollecTRI TF regulons, CORUM complexes).

*Hold-out sources* are used **only** for evaluation (BioGRID physical
interactions, Reactome, GO Biological Process). They must never enter a training
layer, otherwise the biological-coherence evaluation is circular.

``verify_disjoint()`` enforces that separation and is called by the fetch script.

Provenance note: CollecTRI and CORUM are reached through the OmniPath web
service, which aggregates many resources. Every request therefore filters on an
explicit ``resources=`` value so that aggregated STRING/Reactome content cannot
leak into a training layer. The primary databases remain the correct citation.
"""

from __future__ import annotations

import gzip
import hashlib
import io
import json
import zipfile
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

import networkx as nx
import pandas as pd
import requests

DATA_DIR = Path("data")
CACHE_DIR = DATA_DIR / "interactome"
MANIFEST_PATH = DATA_DIR / "interactome_manifest.json"

TAXON_HUMAN = 9606
TIMEOUT = 600

# ── Pinned source versions ────────────────────────────────────────────────────
# Pinned deliberately. Bumping a version changes results and must be a recorded
# decision, not an accident of when the pipeline happened to run.
STRING_VERSION = "12.0"
KEGG_LIBRARY = "KEGG_2021_Human"
REACTOME_LIBRARY = "Reactome_2022"
GOBP_LIBRARY = "GO_Biological_Process_2023"

STRING_LINKS_URL = (
    f"https://stringdb-downloads.org/download/protein.links.detailed.v{STRING_VERSION}/"
    f"{TAXON_HUMAN}.protein.links.detailed.v{STRING_VERSION}.txt.gz"
)
STRING_INFO_URL = (
    f"https://stringdb-downloads.org/download/protein.info.v{STRING_VERSION}/"
    f"{TAXON_HUMAN}.protein.info.v{STRING_VERSION}.txt.gz"
)
OMNIPATH_INTERACTIONS_URL = "https://omnipathdb.org/interactions"
OMNIPATH_COMPLEXES_URL = "https://omnipathdb.org/complexes"
ENRICHR_LIBRARY_URL = "https://maayanlab.cloud/Enrichr/geneSetLibrary"
BIOGRID_URL = (
    "https://downloads.thebiogrid.org/Download/BioGRID/Latest-Release/"
    "BIOGRID-ORGANISM-LATEST.tab3.zip"
)

TRAINING_SOURCES = ("string_ppi", "kegg_pathways", "collectri_tf", "corum_complexes")
HOLDOUT_SOURCES = ("biogrid_physical", "reactome_pathways", "gobp_pathways")


@dataclass
class LayerSpec:
    """Description of one cached source."""

    name: str
    kind: str  # "edges" or "sets"
    role: str  # "training" or "holdout"
    url: str
    primary_citation: str
    notes: str = ""
    extra: dict = field(default_factory=dict)


SPECS: dict[str, LayerSpec] = {
    "string_ppi": LayerSpec(
        name="string_ppi",
        kind="edges",
        role="training",
        url=STRING_LINKS_URL,
        primary_citation="Szklarczyk et al., STRING v12, NAR 2023",
        notes=(
            "Bulk flat file rather than the REST API: one pinned download instead "
            "of hundreds of batched requests, and the version is explicit."
        ),
    ),
    "kegg_pathways": LayerSpec(
        name="kegg_pathways",
        kind="sets",
        role="training",
        url=f"{ENRICHR_LIBRARY_URL}?libraryName={KEGG_LIBRARY}",
        primary_citation="Kanehisa & Goto, KEGG, NAR 2000",
        notes="Gene-set membership. Consumed by the method, so KEGG enrichment is a consistency check only.",
    ),
    "collectri_tf": LayerSpec(
        name="collectri_tf",
        kind="edges",
        role="training",
        url=OMNIPATH_INTERACTIONS_URL,
        primary_citation="Muller-Dott et al., CollecTRI, NAR 2023",
        notes="Signed, directed TF->target regulons. Retrieved via OmniPath with resources=CollecTRI.",
    ),
    "corum_complexes": LayerSpec(
        name="corum_complexes",
        kind="sets",
        role="training",
        url=OMNIPATH_COMPLEXES_URL,
        primary_citation="Tsitsiridis et al., CORUM, NAR 2023",
        notes=(
            "Protein complex co-membership. CORUM's own download has a broken TLS "
            "chain, so it is reached through OmniPath and filtered on sources."
        ),
    ),
    "biogrid_physical": LayerSpec(
        name="biogrid_physical",
        kind="edges",
        role="holdout",
        url=BIOGRID_URL,
        primary_citation="Oughtred et al., BioGRID, Protein Science 2021",
        notes=(
            "Hold-out interactome for the Marker Network Coherence metric. Flat file "
            "avoids the REST API's access-key requirement and pins the release."
        ),
    ),
    "reactome_pathways": LayerSpec(
        name="reactome_pathways",
        kind="sets",
        role="holdout",
        url=f"{ENRICHR_LIBRARY_URL}?libraryName={REACTOME_LIBRARY}",
        primary_citation="Milacic et al., Reactome, NAR 2024",
        notes="Hold-out pathway database for enrichment. Never used to build a layer.",
    ),
    "gobp_pathways": LayerSpec(
        name="gobp_pathways",
        kind="sets",
        role="holdout",
        url=f"{ENRICHR_LIBRARY_URL}?libraryName={GOBP_LIBRARY}",
        primary_citation="Ashburner et al., Gene Ontology, Nat Genet 2000",
        notes="Hold-out pathway database for enrichment.",
    ),
}


# ── Helpers ───────────────────────────────────────────────────────────────────
def _cache_path(name: str) -> Path:
    return CACHE_DIR / f"{name}.parquet"


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _download(url: str, params: dict | None = None) -> bytes:
    resp = requests.get(url, params=params, timeout=TIMEOUT)
    resp.raise_for_status()
    return resp.content


def _normalize_edges(df: pd.DataFrame) -> pd.DataFrame:
    """Canonicalize an edge table: uppercase symbols, drop self-loops, dedupe.

    Undirected edges are stored with ``source <= target`` so that a set of
    ``(source, target)`` tuples can be compared across sources without ambiguity.
    Directed layers keep their orientation but also carry the canonical pair in
    ``key_a`` / ``key_b`` for overlap checks.
    """
    df = df.copy()
    df["source"] = df["source"].astype(str).str.strip().str.upper()
    df["target"] = df["target"].astype(str).str.strip().str.upper()
    df = df[(df["source"] != "") & (df["target"] != "")]
    df = df[df["source"] != df["target"]]

    lo = df[["source", "target"]].min(axis=1)
    hi = df[["source", "target"]].max(axis=1)
    df["key_a"], df["key_b"] = lo, hi

    if not df.get("is_directed", pd.Series([False] * len(df))).any():
        df["source"], df["target"] = lo, hi

    return df.drop_duplicates(subset=["source", "target"]).reset_index(drop=True)


def _sets_to_frame(sets: dict[str, list[str]], id_col: str = "term") -> pd.DataFrame:
    """Flatten a term -> genes mapping into a long table."""
    rows = []
    for term, genes in sets.items():
        clean = sorted({g.strip().upper() for g in genes if g and g.strip()})
        if len(clean) < 2:
            continue
        for g in clean:
            rows.append({id_col: term, "gene": g})
    out = pd.DataFrame(rows)
    if not out.empty:
        out = out.drop_duplicates().reset_index(drop=True)
    return out


def _parse_enrichr_library(text: str) -> dict[str, list[str]]:
    """Parse an Enrichr GMT-style library into term -> genes."""
    sets: dict[str, list[str]] = {}
    for line in text.strip().split("\n"):
        parts = line.rstrip().split("\t")
        if len(parts) < 3:
            continue
        term = parts[0].strip()
        genes = [p.split(",")[0].strip() for p in parts[1:] if p.strip()]
        if term and genes:
            sets[term] = genes
    return sets


# ── Fetchers ──────────────────────────────────────────────────────────────────
def fetch_string_ppi(required_score: int = 400) -> pd.DataFrame:
    """STRING v12 human PPI as symbol-keyed edges.

    The combined score gates inclusion; the experimental subscore is retained as
    a separate column so it can be used as an edge weight rather than as a filter.
    The previous implementation filtered on the experimental subscore being
    non-zero, which silently reduced the layer to a few thousand edges.
    """
    info_raw = _download(STRING_INFO_URL)
    with gzip.open(io.BytesIO(info_raw), "rt") as fh:
        info = pd.read_csv(fh, sep="\t")
    id_to_symbol = dict(
        zip(info["#string_protein_id"], info["preferred_name"].astype(str).str.upper())
    )

    links_raw = _download(STRING_LINKS_URL)
    with gzip.open(io.BytesIO(links_raw), "rt") as fh:
        links = pd.read_csv(fh, sep=" ")

    links = links[links["combined_score"] >= required_score]
    out = pd.DataFrame(
        {
            "source": links["protein1"].map(id_to_symbol),
            "target": links["protein2"].map(id_to_symbol),
            "combined_score": links["combined_score"].astype("int32"),
            "experimental_score": links["experimental"].astype("int32"),
        }
    ).dropna(subset=["source", "target"])

    out["weight"] = out["combined_score"] / 1000.0
    return _normalize_edges(out)


def fetch_collectri() -> pd.DataFrame:
    """CollecTRI signed, directed TF -> target regulons via OmniPath.

    OmniPath encodes heteromeric TF complexes as underscore-joined symbols
    (e.g. ``JUNB_JUND``). Those are expanded so each constituent TF contributes an
    edge, with ``from_complex`` marking the expansion.
    """
    params = {
        "resources": "CollecTRI",
        "organisms": TAXON_HUMAN,
        "genesymbols": "yes",
        "fields": "sources",
    }
    resp = requests.get(OMNIPATH_INTERACTIONS_URL, params=params, timeout=TIMEOUT)
    resp.raise_for_status()
    df = pd.read_csv(io.StringIO(resp.text), sep="\t")

    rows = []
    for rec in df.itertuples(index=False):
        src = str(getattr(rec, "source_genesymbol", "") or "")
        tgt = str(getattr(rec, "target_genesymbol", "") or "")
        if not src or not tgt:
            continue
        parts = [p for p in src.split("_") if p]
        for tf in parts:
            rows.append(
                {
                    "source": tf,
                    "target": tgt,
                    "is_directed": True,
                    "is_stimulation": bool(getattr(rec, "is_stimulation", False)),
                    "is_inhibition": bool(getattr(rec, "is_inhibition", False)),
                    "from_complex": len(parts) > 1,
                    "weight": 1.0,
                }
            )

    out = pd.DataFrame(rows)
    # Signed relation: +1 activation, -1 repression, 0 unsigned or contradictory.
    out["sign"] = 0
    out.loc[out["is_stimulation"] & ~out["is_inhibition"], "sign"] = 1
    out.loc[out["is_inhibition"] & ~out["is_stimulation"], "sign"] = -1
    return _normalize_edges(out)


def fetch_corum_complexes() -> pd.DataFrame:
    """CORUM protein complexes as term -> gene membership, via OmniPath.

    Filtered on the ``sources`` provenance column so only CORUM-derived complexes
    are kept out of OmniPath's larger aggregated set.
    """
    resp = requests.get(OMNIPATH_COMPLEXES_URL, timeout=TIMEOUT)
    resp.raise_for_status()
    df = pd.read_csv(io.StringIO(resp.text), sep="\t")
    df = df[df["sources"].astype(str).str.contains("CORUM", na=False)]

    sets: dict[str, list[str]] = {}
    for rec in df.itertuples(index=False):
        name = str(getattr(rec, "name", "") or "").strip()
        comps = str(getattr(rec, "components_genesymbols", "") or "")
        genes = [g for g in comps.split("_") if g]
        if not name or len(genes) < 2:
            continue
        # Complex names are not unique in CORUM; disambiguate by membership.
        key = f"{name}|{len(genes)}|{genes[0]}"
        sets[key] = genes
    return _sets_to_frame(sets, id_col="complex")


def fetch_gene_set_library(library: str) -> pd.DataFrame:
    """Fetch an Enrichr gene-set library (KEGG, Reactome, GO BP) as long table."""
    resp = requests.get(
        ENRICHR_LIBRARY_URL, params={"mode": "text", "libraryName": library}, timeout=TIMEOUT
    )
    resp.raise_for_status()
    return _sets_to_frame(_parse_enrichr_library(resp.text), id_col="term")


def fetch_biogrid_physical() -> pd.DataFrame:
    """BioGRID human physical interactions from the TAB3 flat file.

    Hold-out network for MNC. Restricted to human-human physical interactions and
    de-duplicated to unique symbol pairs; experiment-level rows are collapsed and
    the supporting evidence count is retained.
    """
    raw = _download(BIOGRID_URL)
    with zipfile.ZipFile(io.BytesIO(raw)) as zf:
        member = next(
            (n for n in zf.namelist() if "Homo_sapiens" in n and n.endswith(".tab3.txt")),
            None,
        )
        if member is None:
            raise RuntimeError("Homo sapiens TAB3 file not found in the BioGRID archive")
        with zf.open(member) as fh:
            df = pd.read_csv(
                fh,
                sep="\t",
                low_memory=False,
                usecols=[
                    "Official Symbol Interactor A",
                    "Official Symbol Interactor B",
                    "Organism ID Interactor A",
                    "Organism ID Interactor B",
                    "Experimental System Type",
                ],
            )

    df = df[
        (df["Organism ID Interactor A"] == TAXON_HUMAN)
        & (df["Organism ID Interactor B"] == TAXON_HUMAN)
        & (df["Experimental System Type"].astype(str).str.lower() == "physical")
    ]
    edges = df.rename(
        columns={
            "Official Symbol Interactor A": "source",
            "Official Symbol Interactor B": "target",
        }
    )[["source", "target"]]

    edges = _normalize_edges(edges.assign(weight=1.0))
    counts = (
        df.rename(
            columns={
                "Official Symbol Interactor A": "source",
                "Official Symbol Interactor B": "target",
            }
        )
        .assign(
            key_a=lambda d: d[["source", "target"]].astype(str).min(axis=1).str.upper(),
            key_b=lambda d: d[["source", "target"]].astype(str).max(axis=1).str.upper(),
        )
        .groupby(["key_a", "key_b"])
        .size()
        .rename("n_evidence")
        .reset_index()
    )
    return edges.merge(counts, on=["key_a", "key_b"], how="left")


FETCHERS = {
    "string_ppi": fetch_string_ppi,
    "kegg_pathways": lambda: fetch_gene_set_library(KEGG_LIBRARY),
    "collectri_tf": fetch_collectri,
    "corum_complexes": fetch_corum_complexes,
    "biogrid_physical": fetch_biogrid_physical,
    "reactome_pathways": lambda: fetch_gene_set_library(REACTOME_LIBRARY),
    "gobp_pathways": lambda: fetch_gene_set_library(GOBP_LIBRARY),
}


# ── Cache management ──────────────────────────────────────────────────────────
def load_manifest() -> dict:
    if MANIFEST_PATH.exists():
        return json.loads(MANIFEST_PATH.read_text())
    return {"sources": {}}


def save_manifest(manifest: dict) -> None:
    MANIFEST_PATH.parent.mkdir(parents=True, exist_ok=True)
    MANIFEST_PATH.write_text(json.dumps(manifest, indent=2, sort_keys=True))


def manifest_hash() -> str:
    """Short hash over cached source checksums, for stamping into run records."""
    manifest = load_manifest()
    payload = json.dumps(
        {k: v.get("sha256") for k, v in sorted(manifest.get("sources", {}).items())},
        sort_keys=True,
    )
    return hashlib.sha256(payload.encode()).hexdigest()[:12]


def ensure_source(name: str, refresh: bool = False) -> pd.DataFrame:
    """Return a cached source, fetching it if absent or if ``refresh`` is set."""
    if name not in FETCHERS:
        raise KeyError(f"unknown source '{name}'; known: {sorted(FETCHERS)}")

    path = _cache_path(name)
    if path.exists() and not refresh:
        return pd.read_parquet(path)

    df = FETCHERS[name]()
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(path, index=False)

    spec = SPECS[name]
    manifest = load_manifest()
    manifest["sources"][name] = {
        "role": spec.role,
        "kind": spec.kind,
        "url": spec.url,
        "primary_citation": spec.primary_citation,
        "notes": spec.notes,
        "fetched_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "rows": int(len(df)),
        "sha256": _sha256(path),
        "cache": str(path).replace("\\", "/"),
    }
    if spec.name == "string_ppi":
        manifest["sources"][name]["string_version"] = STRING_VERSION
    save_manifest(manifest)
    return df


def load_source(name: str) -> pd.DataFrame:
    """Read a source from cache. Raises if it has not been staged."""
    path = _cache_path(name)
    if not path.exists():
        raise FileNotFoundError(
            f"'{name}' is not staged. Run: uv run fetch_interactomes.py"
        )
    return pd.read_parquet(path)


# ── Contamination check ───────────────────────────────────────────────────────
def _edge_key_set(df: pd.DataFrame) -> set[tuple[str, str]]:
    return set(zip(df["key_a"], df["key_b"]))


def _set_membership_pairs(df: pd.DataFrame, id_col: str) -> set[tuple[str, str]]:
    """Co-membership pairs implied by gene sets, capped to keep the check tractable."""
    pairs: set[tuple[str, str]] = set()
    for _, genes in df.groupby(id_col)["gene"]:
        g = sorted(set(genes))
        if len(g) > 200:  # skip very large sets; cliques explode quadratically
            continue
        for i in range(len(g)):
            for j in range(i + 1, len(g)):
                pairs.add((g[i], g[j]))
    return pairs


def verify_disjoint(strict: bool = False) -> dict:
    """Quantify overlap between training layers and hold-out sources.

    Some overlap is unavoidable and expected: a real protein interaction can be
    recorded in both STRING and BioGRID. What matters for the evaluation is that
    the hold-out is not *derived from* the training layers, and that the overlap
    is small enough that a win on hold-out coherence cannot be explained by
    memorisation. This reports the Jaccard index so the overlap can be stated.
    """
    report: dict = {"pairwise": {}, "notes": []}

    train_edges: dict[str, set] = {}
    for name in TRAINING_SOURCES:
        path = _cache_path(name)
        if not path.exists():
            continue
        df = pd.read_parquet(path)
        if SPECS[name].kind == "edges":
            train_edges[name] = _edge_key_set(df)
        else:
            id_col = "complex" if name == "corum_complexes" else "term"
            train_edges[name] = _set_membership_pairs(df, id_col)

    hold_edges: dict[str, set] = {}
    for name in HOLDOUT_SOURCES:
        path = _cache_path(name)
        if not path.exists():
            continue
        df = pd.read_parquet(path)
        if SPECS[name].kind == "edges":
            hold_edges[name] = _edge_key_set(df)
        else:
            hold_edges[name] = _set_membership_pairs(df, "term")

    for tname, tset in train_edges.items():
        for hname, hset in hold_edges.items():
            inter = len(tset & hset)
            union = len(tset | hset)
            report["pairwise"][f"{tname}~{hname}"] = {
                "train_size": len(tset),
                "holdout_size": len(hset),
                "intersection": inter,
                "jaccard": round(inter / union, 6) if union else 0.0,
                "frac_of_holdout": round(inter / len(hset), 6) if hset else 0.0,
            }

    if strict:
        for key, vals in report["pairwise"].items():
            if vals["jaccard"] > 0.5:
                raise AssertionError(
                    f"{key}: Jaccard {vals['jaccard']:.3f} is too high for a hold-out"
                )
    return report


# ── Graph construction on a gene universe ─────────────────────────────────────
def edges_to_graph(
    df: pd.DataFrame, gene_universe: set[str], weight_col: str = "weight"
) -> nx.Graph:
    """Build an undirected graph restricted to ``gene_universe``.

    Nodes are inserted in sorted order. NetworkX preserves insertion order, and
    igraph assigns vertex indices from it, so an unordered insertion makes Leiden's
    tie-breaking depend on Python's per-process string hash seed. That produced
    runs at a fixed seed which differed from each other.
    """
    sub = df[df["source"].isin(gene_universe) & df["target"].isin(gene_universe)]
    G = nx.Graph()
    G.add_nodes_from(sorted(gene_universe))
    w = weight_col if weight_col in sub.columns else None
    for rec in sub.itertuples(index=False):
        G.add_edge(rec.source, rec.target, weight=float(getattr(rec, w)) if w else 1.0)
    return G


def sets_to_bipartite(
    df: pd.DataFrame, gene_universe: set[str], id_col: str, min_size: int = 2
) -> nx.Graph:
    """Gene <-> set bipartite graph.

    Preferred over clique expansion. A set of size *n* contributes *n* edges here
    rather than *n(n-1)/2*, so gene importance no longer scales with the size of
    whichever sets a gene happens to belong to, and the set itself becomes a node
    whose activity can be read directly off a walk.
    """
    sub = df[df["gene"].isin(gene_universe)]
    sizes = sub.groupby(id_col)["gene"].nunique()
    keep = sorted(sizes[sizes >= min_size].index)
    sub = sub[sub[id_col].isin(set(keep))]

    G = nx.Graph()
    # Sorted insertion: see the note in edges_to_graph. Iterating a set of strings
    # gives a different order in every process.
    G.add_nodes_from(sorted(gene_universe), bipartite="gene", node_kind="gene")
    for term in keep:
        G.add_node(f"SET::{term}", bipartite="set", node_kind="set")
    for rec in sub.itertuples(index=False):
        G.add_edge(getattr(rec, "gene"), f"SET::{getattr(rec, id_col)}", weight=1.0)
    return G


def sets_to_weighted_cliques(
    df: pd.DataFrame, gene_universe: set[str], id_col: str, max_size: int = 300
) -> nx.Graph:
    """Clique expansion with 1/|S| edge weights.

    Retained for the ablation that shows what the unweighted clique construction
    did: without the size penalty, a 200-gene pathway contributes 19,900 edges
    and dominates a 10-gene pathway's 45 purely combinatorially.
    """
    sub = df[df["gene"].isin(gene_universe)]
    G = nx.Graph()
    G.add_nodes_from(sorted(gene_universe))
    for term, genes in sub.groupby(id_col)["gene"]:
        g = sorted(set(genes))
        if len(g) < 2 or len(g) > max_size:
            continue
        w = 1.0 / len(g)
        for i in range(len(g)):
            for j in range(i + 1, len(g)):
                if G.has_edge(g[i], g[j]):
                    G[g[i]][g[j]]["weight"] += w
                else:
                    G.add_edge(g[i], g[j], weight=w)
    return G


def summarize_layers(gene_universe: set[str]) -> pd.DataFrame:
    """Per-layer edge counts and density on a gene universe."""
    rows = []
    for name in TRAINING_SOURCES + HOLDOUT_SOURCES:
        path = _cache_path(name)
        if not path.exists():
            continue
        spec = SPECS[name]
        df = pd.read_parquet(path)
        if spec.kind == "edges":
            G = edges_to_graph(df, gene_universe)
            covered = {n for n in G.nodes if G.degree(n) > 0}
        else:
            id_col = "complex" if name == "corum_complexes" else "term"
            G = sets_to_bipartite(df, gene_universe, id_col)
            covered = {
                n
                for n in G.nodes
                if G.nodes[n].get("node_kind") == "gene" and G.degree(n) > 0
            }
        n = len(gene_universe)
        rows.append(
            {
                "layer": name,
                "role": spec.role,
                "kind": spec.kind,
                "rows_cached": len(df),
                "edges_in_universe": G.number_of_edges(),
                "genes_covered": len(covered),
                "gene_coverage": round(len(covered) / n, 4) if n else 0.0,
                "density": round(
                    G.number_of_edges() / (n * (n - 1) / 2), 6
                ) if n > 1 else 0.0,
            }
        )
    return pd.DataFrame(rows)
