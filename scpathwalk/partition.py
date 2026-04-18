"""Leiden partitioning and weighted consensus clustering."""

import igraph as ig
import leidenalg
import networkx as nx
import numpy as np


def leiden_partition(
    G: nx.Graph, seed: int = 0
) -> tuple[list[list[int]], dict[str, int]]:
    """Run Leiden community detection on a NetworkX graph.

    Returns:
        partitions: list of node-index lists per community
        node_to_partition: dict mapping node name -> partition id
    """
    g = ig.Graph.from_networkx(G)
    partitions = list(
        leidenalg.find_partition(
            g, partition_type=leidenalg.ModularityVertexPartition, seed=seed
        )
    )

    node_to_partition = {}
    for part_id, members in enumerate(partitions):
        for node_idx in members:
            node_to_partition[g.vs[node_idx]["_nx_name"]] = part_id

    return partitions, node_to_partition


def map_partitions_to_genes(
    gene_list: list[str], node_to_partition: dict[str, int]
) -> list[int]:
    """Map partition IDs onto a gene list, using -1 for genes not in the graph."""
    return [node_to_partition.get(g, -1) for g in gene_list]


def weighted_consensus_clustering(
    partition_lists: list[list[int]],
    weights: list[float],
) -> list[int]:
    """Build a weighted consensus partition from multiple partition vectors.

    Each gene is assigned to the cluster that received the highest weighted vote.
    """
    num_samples = len(partition_lists[0])
    num_clusters = max(max(pl) for pl in partition_lists) + 1
    consensus_matrix = np.zeros((num_samples, num_clusters), dtype=float)

    for weight, pl in zip(weights, partition_lists):
        for j, cluster_id in enumerate(pl):
            if cluster_id >= 0:
                consensus_matrix[j, cluster_id] += weight

    return np.argmax(consensus_matrix, axis=1).tolist()
