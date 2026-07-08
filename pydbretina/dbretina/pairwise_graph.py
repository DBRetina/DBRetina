"""In-memory igraph-backed graph engine for DBRetina pairwise data.

``PairwiseGraph`` builds an undirected :class:`igraph.Graph` from filtered pairwise
results and exposes neighborhood traversal, shortest paths, connected components,
community detection and PageRank — all delegated to igraph (no custom graph
algorithms).

Vertices are the full node universe from ``store.get_names_map()`` (so isolated
groups exist as degree-0 vertices). Edges come from
``store.to_pandas(metric=, cutoff=)``; each edge carries the filter ``weight`` and
``shared_features`` as edge attributes. Group ids (``gid``) are stable across
subgraphs; the ``gid <-> vertex index`` mapping is kept internally.

Usage::

    store = PairwiseStore("/path/to/experiment_pairwise")
    graph = PairwiseGraph(store, metric="jaccard", cutoff=30)

    components = graph.connected_components()   # {name: component_id}
    pr = graph.pagerank()                        # {name: score}
    sub = graph.subgraph_around("kegg_alzheimer", hops=2)
"""

import warnings
from typing import Optional

import pandas as pd

from .pairwise_store import PairwiseStore


class PairwiseGraph:
    """Graph query engine over DBRetina pairwise data using in-memory igraph."""

    def __init__(
        self,
        store: PairwiseStore,
        metric: str,
        cutoff: float = 0.0,
        db_path: Optional[str] = None,  # accepted for API compat; unused
    ):
        """Build a graph from filtered pairwise data.

        Args:
            store: PairwiseStore instance to read data from.
            metric: Metric to filter edges on (e.g. "jaccard"); also the edge weight.
            cutoff: Minimum metric value to include an edge.
            db_path: Ignored (kept for backwards-compatible signature).
        """
        self._store = store
        self._metric = metric
        self._cutoff = cutoff
        self._names_map = store.get_names_map()
        self._name_to_id = {v: k for k, v in self._names_map.items()}

        # Carry *every* metric the dataset has as an edge attribute so the
        # client can filter on any metric without re-querying. ``weight`` stays
        # the active filter metric (the graph algorithms use it); the active
        # metric is ALSO kept as its own named attribute so the client can
        # filter on it like any other.
        self._metric_cols = list(store.available_metrics)
        if metric not in self._metric_cols:
            self._metric_cols.append(metric)
        sel_cols = ["group_1_id", "group_2_id", "shared_features"] + self._metric_cols
        df = store.to_pandas(metric=metric, cutoff=cutoff, columns=sel_cols)
        self._build_from_frame(df)

    def _build_from_frame(self, df: pd.DataFrame):
        """Construct the igraph.Graph from an edges DataFrame (gid columns)."""
        import igraph as ig

        # Stable, sorted vertex ordering keyed by group id.
        gids = sorted(self._names_map.keys())
        self._gid_to_idx = {gid: i for i, gid in enumerate(gids)}
        self._idx_to_gid = gids

        g = ig.Graph(n=len(gids), directed=False)
        g.vs["gid"] = gids
        g.vs["name"] = [self._names_map[gid] for gid in gids]

        # Which metric columns are actually present in this frame? (subgraph
        # rebuilds carry whatever edges_dataframe() produced.)
        metric = self._metric
        cols = set(df.columns)
        metric_cols = [m for m in getattr(self, "_metric_cols", []) if m in cols]
        # The active filter metric must always be readable for ``weight``.
        if metric not in cols and metric in metric_cols:
            metric_cols.remove(metric)
        self._metric_cols = metric_cols

        edges, weights, shared = [], [], []
        metric_vals = {m: [] for m in metric_cols}
        idx = self._gid_to_idx
        for row in df.itertuples(index=False):
            s = int(row.group_1_id)
            d = int(row.group_2_id)
            if s in idx and d in idx:
                edges.append((idx[s], idx[d]))
                weights.append(float(getattr(row, metric)) if metric in cols else 0.0)
                shared.append(int(row.shared_features))
                for m in metric_cols:
                    metric_vals[m].append(float(getattr(row, m)))

        if edges:
            g.add_edges(edges)
            g.es["weight"] = weights
            g.es["shared_features"] = shared
            for m in metric_cols:
                g.es[m] = metric_vals[m]

        self._g = g
        self._num_nodes = g.vcount()
        self._num_edges = g.ecount()

    # ── Lifecycle ─────────────────────────────────────────────────────

    def close(self):
        """Release the in-memory graph (no external resources to free)."""
        self._g = None

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()

    # ── Properties ────────────────────────────────────────────────────

    @property
    def num_nodes(self) -> int:
        return self._num_nodes

    @property
    def num_edges(self) -> int:
        return self._num_edges

    @property
    def metric(self) -> str:
        return self._metric

    @property
    def cutoff(self) -> float:
        return self._cutoff

    # ── Edge access ───────────────────────────────────────────────────

    def edges_dataframe(self) -> pd.DataFrame:
        """Return this graph's edges as a DataFrame keyed by group id.

        Columns: ``src``/``dst`` (int group ids), ``weight`` (float, the filter
        metric value), ``shared_features`` (int), plus one float column per
        available metric (``ochiai``, ``jaccard``, …). One row per undirected
        edge. The extra metric columns let ``subgraph`` round-trip every metric.
        """
        g = self._g
        idx_to_gid = self._idx_to_gid
        metric_cols = list(getattr(self, "_metric_cols", []))
        if g.ecount() == 0:
            data = {
                "src": pd.Series([], dtype="int64"),
                "dst": pd.Series([], dtype="int64"),
                "weight": pd.Series([], dtype="float64"),
                "shared_features": pd.Series([], dtype="int64"),
            }
            for m in metric_cols:
                data[m] = pd.Series([], dtype="float64")
            return pd.DataFrame(data)
        weights = g.es["weight"]
        shared = g.es["shared_features"]
        src, dst = [], []
        for e in g.es:
            src.append(idx_to_gid[e.source])
            dst.append(idx_to_gid[e.target])
        data = {
            "src": pd.array(src, dtype="int64"),
            "dst": pd.array(dst, dtype="int64"),
            "weight": pd.array([float(w) for w in weights], dtype="float64"),
            "shared_features": pd.array([int(s) for s in shared], dtype="int64"),
        }
        for m in metric_cols:
            data[m] = pd.array([float(v) for v in g.es[m]], dtype="float64")
        return pd.DataFrame(data)

    def to_igraph(self):
        """Return the underlying :class:`igraph.Graph` (vertices carry ``gid``/``name``)."""
        return self._g

    # ── Graph Operations ──────────────────────────────────────────────

    def subgraph(self, group_names: list[str]) -> "PairwiseGraph":
        """Extract a subgraph containing only the specified groups.

        Args:
            group_names: List of group names to include.

        Returns:
            A new PairwiseGraph with those nodes and the induced edges.
        """
        gids = []
        for name in group_names:
            gid = self._name_to_id.get(name.lower())
            if gid is None:
                raise KeyError(f"Group not found: {name}")
            gids.append(gid)
        gid_set = set(gids)
        sub_names = {gid: self._names_map[gid] for gid in gid_set}

        df = self.edges_dataframe()
        if not df.empty:
            mask = df["src"].isin(gid_set) & df["dst"].isin(gid_set)
            df = df.loc[mask]

        return PairwiseGraph._from_edge_frame(
            df, sub_names, self._metric, self._cutoff, self._store
        )

    def subgraph_around(self, center: str, hops: int = 1) -> "PairwiseGraph":
        """Extract a subgraph around a center node within N hops.

        Args:
            center: Name of the center group.
            hops: Number of hops to include (default 1 = direct neighbors).

        Returns:
            A new PairwiseGraph with the center, all nodes within ``hops``, and
            the edges induced among them.
        """
        gid = self._resolve_id(center)
        center_idx = self._gid_to_idx[gid]
        # igraph neighborhood includes the center vertex itself.
        idxs = self._g.neighborhood(vertices=center_idx, order=hops)
        idx_to_gid = self._idx_to_gid
        group_names = [self._names_map[idx_to_gid[i]] for i in idxs]
        return self.subgraph(group_names)

    def connected_components(
        self, group_names: Optional[list[str]] = None
    ) -> dict[str, int]:
        """Find weakly connected components.

        Args:
            group_names: Optional subset of group names to restrict analysis to.
                         If None, uses all nodes in the graph.

        Returns:
            Dict mapping group_name -> component_id.
        """
        if group_names is not None:
            return self.subgraph(group_names).connected_components()

        membership = self._g.connected_components(mode="weak").membership
        names = self._g.vs["name"]
        return {names[i]: membership[i] for i in range(self._g.vcount())}

    def community_detection(
        self, method: str = "leiden", group_names: Optional[list[str]] = None
    ) -> dict[str, int]:
        """Detect communities using igraph.

        Args:
            method: "leiden" (modularity objective) or "louvain" (multilevel).
            group_names: Optional subset of group names to restrict analysis to.

        Returns:
            Dict mapping group_name -> community_id.
        """
        if group_names is not None:
            return self.subgraph(group_names).community_detection(method=method)

        g = self._g
        weights = "weight" if g.ecount() > 0 and "weight" in g.es.attributes() else None
        if method == "leiden":
            membership = g.community_leiden(
                objective_function="modularity", weights=weights
            ).membership
        elif method == "louvain":
            membership = g.community_multilevel(weights=weights).membership
        else:
            raise ValueError(f"Unknown method '{method}'. Use 'leiden' or 'louvain'.")

        names = g.vs["name"]
        return {names[i]: membership[i] for i in range(g.vcount())}

    def pagerank(
        self, damping: float = 0.85, group_names: Optional[list[str]] = None
    ) -> dict[str, float]:
        """Compute PageRank scores.

        Args:
            damping: Damping factor (default 0.85).
            group_names: Optional subset of group names to restrict PageRank to.

        Returns:
            Dict mapping group_name -> PageRank score.
        """
        if group_names is not None:
            return self.subgraph(group_names).pagerank(damping=damping)

        g = self._g
        weights = "weight" if g.ecount() > 0 and "weight" in g.es.attributes() else None
        pr = g.pagerank(damping=damping, weights=weights)
        names = g.vs["name"]
        return {names[i]: pr[i] for i in range(g.vcount())}

    def shortest_path_full(self, source: str, target: str) -> dict:
        """Find the shortest path between two groups (by hop count).

        Args:
            source: Source group name.
            target: Target group name.

        Returns:
            Dict with ``source``, ``target``, ``path_length`` (hops; -1 if
            disconnected), ``path_nodes`` (list of group names) and ``connected``.
        """
        src_id = self._resolve_id(source)
        tgt_id = self._resolve_id(target)
        src_idx = self._gid_to_idx[src_id]
        tgt_idx = self._gid_to_idx[tgt_id]

        # Unweighted shortest path (by number of edges), matching the previous
        # hop-count semantics. igraph warns to stderr when the target is unreachable;
        # that's an expected case (we report connected=False below), so silence it.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            paths = self._g.get_shortest_paths(src_idx, to=tgt_idx, output="vpath")
        path = paths[0] if paths else []

        if not path:
            return {
                "source": source,
                "target": target,
                "path_length": -1,
                "path_nodes": [],
                "connected": False,
            }

        names = self._g.vs["name"]
        path_nodes = [names[i] for i in path]
        return {
            "source": source,
            "target": target,
            "path_length": len(path) - 1,
            "path_nodes": path_nodes,
            "connected": True,
        }

    # ── Internal ──────────────────────────────────────────────────────

    def _resolve_id(self, group_name: str) -> int:
        gid = self._name_to_id.get(group_name.lower())
        if gid is None:
            raise KeyError(f"Group not found: {group_name}")
        return gid

    @classmethod
    def _from_edge_frame(
        cls,
        edges_df: pd.DataFrame,
        names_map: dict[int, str],
        metric: str,
        cutoff: float,
        store: PairwiseStore,
    ) -> "PairwiseGraph":
        """Build a PairwiseGraph directly from an edges DataFrame (internal).

        ``edges_df`` must have columns ``src``, ``dst``, ``weight``,
        ``shared_features`` keyed by group id, plus one float column per available
        metric (e.g. ``ochiai``, ``jaccard``, …) which are carried onto the edges.
        Used by ``subgraph``.
        """
        obj = object.__new__(cls)
        obj._store = store
        obj._metric = metric
        obj._cutoff = cutoff
        obj._names_map = dict(names_map)
        obj._name_to_id = {v: k for k, v in names_map.items()}

        # Per-metric columns carried by edges_dataframe() (anything that isn't a
        # structural column). _build_from_frame reads these back as edge attrs.
        structural = {"src", "dst", "weight", "shared_features"}
        metric_cols = [c for c in edges_df.columns if c not in structural]
        obj._metric_cols = metric_cols

        # Reshape to the gid-column layout _build_from_frame expects.
        if edges_df.empty:
            frame = pd.DataFrame(
                columns=["group_1_id", "group_2_id", "shared_features", metric]
                + [m for m in metric_cols if m != metric]
            )
        else:
            cols = {
                "group_1_id": edges_df["src"].to_numpy(),
                "group_2_id": edges_df["dst"].to_numpy(),
                "shared_features": edges_df["shared_features"].to_numpy(),
                metric: edges_df["weight"].to_numpy(),
            }
            for m in metric_cols:
                cols[m] = edges_df[m].to_numpy()
            frame = pd.DataFrame(cols)
        obj._build_from_frame(frame)
        return obj

    def __repr__(self):
        return (
            f"PairwiseGraph(metric='{self._metric}', cutoff={self._cutoff}, "
            f"nodes={self._num_nodes}, edges={self._num_edges:,})"
        )

    def _repr_html_(self):
        """Rich HTML display for Jupyter notebooks."""
        try:
            n_components = len(self._g.connected_components(mode="weak"))
        except Exception:
            n_components = "?"

        try:
            degrees = self._g.degree()
            avg_deg = sum(degrees) / len(degrees) if degrees else 0.0
            max_deg = max(degrees) if degrees else 0
            avg_deg_str = f"{avg_deg:.1f}"
        except Exception:
            avg_deg_str = "?"
            max_deg = "?"

        return (
            f"<div style='border:1px solid #ddd;padding:12px;border-radius:6px;"
            f"max-width:400px;font-family:sans-serif;font-size:13px'>"
            f"<b>PairwiseGraph</b><br>"
            f"<table style='margin:8px 0;border-collapse:collapse'>"
            f"<tr><td>Nodes</td><td><b>{self._num_nodes:,}</b></td></tr>"
            f"<tr><td>Edges</td><td><b>{self._num_edges:,}</b></td></tr>"
            f"<tr><td>Filter</td><td>{self._metric} &ge; {self._cutoff}</td></tr>"
            f"<tr><td>Components</td><td>{n_components}</td></tr>"
            f"<tr><td>Avg degree</td><td>{avg_deg_str}</td></tr>"
            f"<tr><td>Max degree</td><td>{max_deg}</td></tr>"
            f"</table>"
            f"</div>"
        )
