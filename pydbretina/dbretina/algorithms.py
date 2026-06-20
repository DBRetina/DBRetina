"""Graph clustering algorithms for the REST API, backed by igraph.

This module is the registry + dispatcher used by the ``/api/v1/graph/cluster``,
``/api/v1/graph/components`` and ``/api/v1/algorithms/clustering`` endpoints. It wraps
igraph's community-detection and connected-components implementations — no custom graph
algorithms are implemented here.

``run_clustering`` operates on an :class:`igraph.Graph` and returns membership keyed by
**node index** (the integer vertex id used to build the graph), matching how the REST
endpoints map back to group names.
"""

import math
from collections import Counter
from dataclasses import dataclass


@dataclass
class ClusteringResult:
    """Result of a clustering run.

    membership: {node_index: cluster_id}; a cluster_id of ``-1`` means the node was
                filtered out (e.g. a connected component smaller than ``min_size``).
    """

    membership: dict
    num_clusters: int
    cluster_sizes: dict
    modularity: float
    algorithm: str
    parameters: dict


ALGORITHMS = {
    "leiden": {
        "name": "leiden",
        "display_name": "Leiden",
        "description": "Leiden community detection (modularity objective). Recommended.",
        "parameters": {
            "resolution": "float, optional (default 1.0)",
            "n_iterations": "int, optional",
        },
    },
    "louvain": {
        "name": "louvain",
        "display_name": "Louvain (multilevel)",
        "description": "Louvain (multilevel) modularity-based community detection.",
        "parameters": {},
    },
    "connected_components": {
        "name": "connected_components",
        "display_name": "Connected Components",
        "description": (
            "Weakly connected components. Components smaller than min_size are "
            "dropped (assigned cluster id -1)."
        ),
        "parameters": {"min_size": "int, optional (default 1)"},
    },
}


def list_algorithms() -> list:
    """Return the available clustering algorithms with their metadata."""
    return list(ALGORITHMS.values())


def get_algorithm(name: str) -> dict:
    """Return the metadata for a single algorithm, or raise ValueError if unknown."""
    if name not in ALGORITHMS:
        raise ValueError(
            f"Unknown algorithm '{name}'. Valid: {', '.join(sorted(ALGORITHMS))}"
        )
    return ALGORITHMS[name]


def _resolve_weights(graph, weights):
    """Prefer explicit weights (if the right length), else the graph's edge weights."""
    if weights is not None and len(weights) == graph.ecount():
        return list(weights)
    if "weight" in graph.es.attributes():
        return graph.es["weight"]
    return None


def run_clustering(graph, algorithm: str = "leiden", weights=None, min_size: int = 1, **params):
    """Cluster an :class:`igraph.Graph`; returns a :class:`ClusteringResult`.

    Args:
        graph: an undirected igraph.Graph.
        algorithm: one of ``ALGORITHMS`` ("leiden", "louvain", "connected_components").
        weights: optional per-edge weights (list aligned to graph.es).
        min_size: connected_components only — drop components smaller than this.
        **params: forwarded to the underlying igraph community method.
    """
    if algorithm not in ALGORITHMS:
        raise ValueError(
            f"Unknown algorithm '{algorithm}'. Valid: {', '.join(sorted(ALGORITHMS))}"
        )

    if graph.vcount() == 0:
        return ClusteringResult({}, 0, {}, 0.0, algorithm, dict(params))

    if algorithm == "connected_components":
        raw = list(graph.connected_components(mode="weak").membership)
        sizes = Counter(raw)
        kept = sorted(cid for cid, sz in sizes.items() if sz >= max(1, int(min_size)))
        relabel = {cid: i for i, cid in enumerate(kept)}
        membership = {idx: relabel.get(cid, -1) for idx, cid in enumerate(raw)}
        cluster_sizes = {relabel[cid]: sizes[cid] for cid in kept}
        return ClusteringResult(
            membership, len(kept), cluster_sizes, 0.0, algorithm, dict(params)
        )

    # Community detection (leiden / louvain)
    w = _resolve_weights(graph, weights)
    if algorithm == "leiden":
        vc = graph.community_leiden(objective_function="modularity", weights=w, **params)
    else:  # louvain
        vc = graph.community_multilevel(weights=w, **params)

    membership = {idx: cid for idx, cid in enumerate(vc.membership)}
    cluster_sizes = dict(Counter(vc.membership))
    try:
        modularity = float(vc.modularity)
        if not math.isfinite(modularity):  # no-edge graphs -> NaN, not JSON-serializable
            modularity = 0.0
    except Exception:
        modularity = 0.0
    return ClusteringResult(
        membership, len(cluster_sizes), cluster_sizes, modularity, algorithm, dict(params)
    )
