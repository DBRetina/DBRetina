"""KuzuDB-backed graph query engine for DBRetina pairwise data."""

import pathlib
import shutil
import tempfile
from typing import Optional

try:
    import kuzu
except ImportError:
    kuzu = None  # type: ignore


def _require_kuzu():
    if kuzu is None:
        raise ImportError(
            "KuzuDB is required for graph features. "
            "Install with: pip install 'DBRetina[server]'"
        )


import pandas as pd
import pyarrow as pa

from .pairwise_store import PairwiseStore


class PairwiseGraph:
    """Graph query engine over DBRetina pairwise data using KuzuDB.

    Builds an embedded graph database from filtered pairwise results,
    enabling Cypher queries, neighborhood traversal, shortest paths,
    and graph algorithm wrappers.

    Usage::

        store = PairwiseStore("/path/to/experiment_pairwise")
        graph = PairwiseGraph(store, metric="jaccard", cutoff=30)

        # Cypher queries
        df = graph.cypher(
            "MATCH (a:Group)-[r:SIMILAR_TO]->(b:Group) "
            "RETURN a.name, b.name, r.jaccard "
            "ORDER BY r.jaccard DESC LIMIT 10"
        )

        # Neighborhood
        neighbors = graph.neighbors("kegg_alzheimer", hops=2)

        # Graph algorithms
        components = graph.connected_components()
        pr = graph.pagerank()
    """

    def __init__(
        self,
        store: PairwiseStore,
        metric: str,
        cutoff: float = 0.0,
        db_path: Optional[str] = None,
    ):
        """Create a graph from filtered pairwise data.

        Args:
            store: PairwiseStore instance to read data from.
            metric: Metric to filter edges on (e.g. "jaccard").
            cutoff: Minimum metric value to include an edge.
            db_path: Optional path for the KuzuDB database.
                     If None, uses a temporary directory.
        """
        _require_kuzu()
        self._store = store
        self._metric = metric
        self._cutoff = cutoff
        self._names_map = store.get_names_map()
        self._name_to_id = {v: k for k, v in self._names_map.items()}

        # Create KuzuDB database
        if db_path:
            self._db_dir = pathlib.Path(db_path)
            self._owns_dir = False
        else:
            # KuzuDB creates the directory itself; give it a non-existent path
            tmp_parent = pathlib.Path(tempfile.mkdtemp(prefix="dbretina_graph_"))
            self._db_dir = tmp_parent / "db"
            self._owns_dir = True
            self._tmp_parent = tmp_parent

        self._db = kuzu.Database(str(self._db_dir))
        self._conn = kuzu.Connection(self._db)

        # Build the graph
        self._build_graph()

    def _build_graph(self):
        """Populate KuzuDB with nodes and edges from the store."""
        conn = self._conn

        # Create node table
        conn.execute(
            "CREATE NODE TABLE IF NOT EXISTS `Group`("
            "id UINT32, name STRING, PRIMARY KEY(id))"
        )

        # Create relationship table with all metrics
        rel_cols = (
            "shared_features UINT64, "
            "containment FLOAT, "
            "ochiai FLOAT, "
            "jaccard FLOAT, "
            "csi FLOAT, "
            "dice FLOAT, "
            "odds_ratio FLOAT"
        )
        if self._store.has_pvalue:
            rel_cols += ", pvalue DOUBLE"

        conn.execute(
            f'CREATE REL TABLE IF NOT EXISTS SIMILAR_TO('
            f'FROM `Group` TO `Group`, {rel_cols})'
        )

        # Insert nodes via CSV COPY
        names_map = self._names_map
        tmp_parent = self._db_dir.parent
        if names_map:
            node_csv = tmp_parent / "_nodes.csv"
            with open(node_csv, "w", newline="") as f:
                f.write("id,name\n")
                for gid, name in names_map.items():
                    # Escape quotes in names for CSV
                    escaped = name.replace('"', '""')
                    f.write(f'{gid},"{escaped}"\n')
            conn.execute(
                f'COPY `Group` FROM "{node_csv}" (HEADER=TRUE)'
            )
            node_csv.unlink(missing_ok=True)

        # Insert edges from filtered pairs via CSV COPY
        edge_cols = [
            "group_1_id", "group_2_id", "shared_features",
            "containment", "ochiai", "jaccard", "csi", "dice", "odds_ratio",
        ]
        if self._store.has_pvalue:
            edge_cols.append("pvalue")

        df = self._store.to_pandas(
            metric=self._metric, cutoff=self._cutoff, columns=edge_cols
        )

        if len(df) > 0:
            edge_csv = tmp_parent / "_edges.csv"
            df.to_csv(edge_csv, index=False, header=True)
            conn.execute(
                f'COPY SIMILAR_TO FROM "{edge_csv}" '
                f'(HEADER=TRUE)'
            )
            edge_csv.unlink(missing_ok=True)

        self._num_nodes = len(names_map)
        self._num_edges = len(df)

    def close(self):
        """Close the graph database."""
        self._conn = None
        self._db = None
        if self._owns_dir:
            parent = getattr(self, "_tmp_parent", self._db_dir)
            if parent.exists():
                shutil.rmtree(parent, ignore_errors=True)

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

    # ── Cypher Queries ────────────────────────────────────────────────

    def cypher(self, query: str, parameters: Optional[dict] = None) -> pd.DataFrame:
        """Execute a Cypher query and return results as a DataFrame.

        The graph has:
        - Node table ``Group`` with properties: ``id`` (UINT32), ``name`` (STRING)
        - Relationship table ``SIMILAR_TO`` with properties:
          ``shared_features``, ``containment``, ``ochiai``, ``jaccard``,
          ``csi``, ``dice``, ``odds_ratio``, and optionally ``pvalue``.

        Args:
            query: Cypher query string.
            parameters: Optional query parameters.

        Returns:
            pandas DataFrame with query results.
        """
        if parameters:
            result = self._conn.execute(query, parameters=parameters)
        else:
            result = self._conn.execute(query)
        return result.get_as_df()

    # ── Graph Operations ──────────────────────────────────────────────

    def neighbors(self, group_name: str, hops: int = 1) -> pd.DataFrame:
        """Find neighbors of a group within N hops.

        Args:
            group_name: Name of the source group.
            hops: Number of hops (1 = direct neighbors).

        Returns:
            DataFrame with neighbor names and edge metrics.
        """
        gid = self._resolve_id(group_name)
        if hops == 1:
            return self.cypher(
                f'MATCH (a:`Group`)-[r:SIMILAR_TO]-(b:`Group`) '
                f"WHERE a.id = $gid "
                f"RETURN b.name AS neighbor, b.id AS neighbor_id, "
                f"r.{self._metric} AS {self._metric}, "
                f"r.shared_features AS shared_features "
                f"ORDER BY r.{self._metric} DESC",
                parameters={"gid": gid},
            )
        else:
            return self.cypher(
                f'MATCH (a:`Group`)-[r:SIMILAR_TO* 1..{hops}]-(b:`Group`) '
                f"WHERE a.id = $gid AND a.id <> b.id "
                f"RETURN DISTINCT b.name AS neighbor, b.id AS neighbor_id",
                parameters={"gid": gid},
            )

    def shortest_path(self, source: str, target: str) -> pd.DataFrame:
        """Find the shortest path between two groups.

        Args:
            source: Source group name.
            target: Target group name.

        Returns:
            DataFrame with path information.
        """
        src_id = self._resolve_id(source)
        tgt_id = self._resolve_id(target)
        return self.cypher(
            f'MATCH (a:`Group`)-[r:SIMILAR_TO* SHORTEST 1..30]-(b:`Group`) '
            f"WHERE a.id = $src AND b.id = $tgt "
            f"RETURN length(r) AS path_length",
            parameters={"src": src_id, "tgt": tgt_id},
        )

    def degree_distribution(self) -> pd.DataFrame:
        """Compute degree distribution of the graph."""
        return self.cypher(
            'MATCH (a:`Group`)-[r:SIMILAR_TO]-() '
            "RETURN a.name AS group_name, a.id AS group_id, COUNT(r) AS degree "
            "ORDER BY degree DESC"
        )

    def subgraph(self, group_names: list[str]) -> "PairwiseGraph":
        """Extract a subgraph containing only the specified groups.

        Creates a new PairwiseGraph with only the specified nodes and
        edges between them.

        Args:
            group_names: List of group names to include.

        Returns:
            A new PairwiseGraph instance containing only the specified groups.
        """
        gids = []
        for name in group_names:
            gid = self._name_to_id.get(name.lower())
            if gid is None:
                raise KeyError(f"Group not found: {name}")
            gids.append(gid)
        gid_set = set(gids)

        placeholders = ", ".join(str(g) for g in gid_set)
        pvalue_col = ", r.pvalue AS pvalue" if self._store.has_pvalue else ""
        edges_df = self.cypher(
            f'MATCH (a:`Group`)-[r:SIMILAR_TO]->(b:`Group`) '
            f'WHERE a.id IN [{placeholders}] AND b.id IN [{placeholders}] '
            f'RETURN a.id AS group_1_id, b.id AS group_2_id, '
            f'r.shared_features AS shared_features, '
            f'r.containment AS containment, r.ochiai AS ochiai, '
            f'r.jaccard AS jaccard, r.csi AS csi, r.dice AS dice, '
            f'r.odds_ratio AS odds_ratio'
            + pvalue_col
        )

        return PairwiseGraph._from_edges(
            edges_df=edges_df,
            names_map={gid: self._names_map[gid] for gid in gid_set},
            metric=self._metric,
            cutoff=self._cutoff,
            has_pvalue=self._store.has_pvalue,
            store=self._store,
        )

    def subgraph_around(self, center: str, hops: int = 1) -> "PairwiseGraph":
        """Extract a subgraph around a center node within N hops.

        Args:
            center: Name of the center group.
            hops: Number of hops to include (default 1 = direct neighbors).

        Returns:
            A new PairwiseGraph containing the center and all nodes
            within the specified number of hops, plus edges between them.
        """
        gid = self._resolve_id(center)
        neighbors_df = self.cypher(
            f'MATCH (a:`Group`)-[r:SIMILAR_TO* 1..{hops}]-(b:`Group`) '
            f'WHERE a.id = $gid '
            f'RETURN DISTINCT b.id AS neighbor_id',
            parameters={"gid": gid},
        )

        neighbor_ids = set(neighbors_df["neighbor_id"].tolist())
        neighbor_ids.add(gid)

        group_names = [
            self._names_map[nid] for nid in neighbor_ids
            if nid in self._names_map
        ]
        return self.subgraph(group_names)

    def connected_components(
        self, group_names: Optional[list[str]] = None
    ) -> dict[str, int]:
        """Find weakly connected components.

        Args:
            group_names: Optional list of group names to restrict analysis to.
                         If None, uses all nodes in the graph.

        Returns:
            Dict mapping group_name -> component_id.
        """
        try:
            import igraph as ig
        except ImportError:
            raise ImportError(
                "igraph is required for connected_components(). "
                "Install with: pip install igraph"
            )

        edges_df = self.cypher(
            'MATCH (a:`Group`)-[r:SIMILAR_TO]->(b:`Group`) '
            "RETURN a.id AS src, b.id AS dst"
        )

        all_ids = self._resolve_group_ids(group_names)
        if edges_df.empty:
            names = self._names_map if group_names is None else {
                gid: self._names_map[gid] for gid in all_ids
            }
            return {name: i for i, name in enumerate(names.values())}

        id_idx = {gid: i for i, gid in enumerate(all_ids)}
        g = ig.Graph(n=len(all_ids), directed=False)
        edge_list = [
            (id_idx[row["src"]], id_idx[row["dst"]])
            for _, row in edges_df.iterrows()
            if row["src"] in id_idx and row["dst"] in id_idx
        ]
        g.add_edges(edge_list)
        membership = g.connected_components().membership
        return {
            self._names_map[gid]: membership[i]
            for i, gid in enumerate(all_ids)
            if gid in self._names_map
        }

    def community_detection(
        self, method: str = "leiden", group_names: Optional[list[str]] = None
    ) -> dict[str, int]:
        """Detect communities using igraph algorithms.

        Args:
            method: "leiden" or "louvain".
            group_names: Optional list of group names to restrict analysis to.
                         If None, uses all nodes in the graph.

        Returns:
            Dict mapping group_name -> community_id.
        """
        try:
            import igraph as ig
        except ImportError:
            raise ImportError(
                "igraph is required for community_detection(). "
                "Install with: pip install igraph"
            )

        edges_df = self.cypher(
            f'MATCH (a:`Group`)-[r:SIMILAR_TO]->(b:`Group`) '
            f"RETURN a.id AS src, b.id AS dst, r.{self._metric} AS weight"
        )

        all_ids = self._resolve_group_ids(group_names)
        if edges_df.empty:
            names = self._names_map if group_names is None else {
                gid: self._names_map[gid] for gid in all_ids
            }
            return {name: 0 for name in names.values()}

        id_idx = {gid: i for i, gid in enumerate(all_ids)}
        g = ig.Graph(n=len(all_ids), directed=False)
        edge_list = []
        weights = []
        for _, row in edges_df.iterrows():
            if row["src"] in id_idx and row["dst"] in id_idx:
                edge_list.append((id_idx[row["src"]], id_idx[row["dst"]]))
                weights.append(float(row["weight"]))
        g.add_edges(edge_list)
        if weights:
            g.es["weight"] = weights

        if method == "leiden":
            membership = g.community_leiden(weights="weight" if weights else None).membership
        elif method == "louvain":
            membership = g.community_multilevel(weights="weight" if weights else None).membership
        else:
            raise ValueError(f"Unknown method '{method}'. Use 'leiden' or 'louvain'.")

        return {
            self._names_map[gid]: membership[i]
            for i, gid in enumerate(all_ids)
            if gid in self._names_map
        }

    def pagerank(
        self, damping: float = 0.85, group_names: Optional[list[str]] = None
    ) -> dict[str, float]:
        """Compute PageRank scores.

        Args:
            damping: Damping factor (default 0.85).
            group_names: Optional list of group names to restrict PageRank to.
                         If None, uses all nodes in the graph.

        Returns:
            Dict mapping group_name -> PageRank score.
        """
        try:
            import igraph as ig
        except ImportError:
            raise ImportError(
                "igraph is required for pagerank(). "
                "Install with: pip install igraph"
            )

        edges_df = self.cypher(
            f'MATCH (a:`Group`)-[r:SIMILAR_TO]->(b:`Group`) '
            f"RETURN a.id AS src, b.id AS dst, r.{self._metric} AS weight"
        )

        all_ids = self._resolve_group_ids(group_names)
        id_idx = {gid: i for i, gid in enumerate(all_ids)}
        g = ig.Graph(n=len(all_ids), directed=False)
        edge_list = []
        weights = []
        for _, row in edges_df.iterrows():
            if row["src"] in id_idx and row["dst"] in id_idx:
                edge_list.append((id_idx[row["src"]], id_idx[row["dst"]]))
                weights.append(float(row["weight"]))
        g.add_edges(edge_list)
        if weights:
            g.es["weight"] = weights

        pr = g.pagerank(damping=damping, weights="weight" if weights else None)
        return {
            self._names_map[gid]: pr[i]
            for i, gid in enumerate(all_ids)
            if gid in self._names_map
        }

    # ── Name Resolution ──────────────────────────────────────────────

    def search_groups(self, pattern: str) -> dict[int, str]:
        """Search for groups whose name contains the pattern (case-insensitive).

        Delegates to the underlying PairwiseStore.
        """
        return self._store.search_groups(pattern)

    # ── Internal ──────────────────────────────────────────────────────

    def _resolve_id(self, group_name: str) -> int:
        gid = self._name_to_id.get(group_name.lower())
        if gid is None:
            raise KeyError(f"Group not found: {group_name}")
        return gid

    def _resolve_group_ids(self, group_names: Optional[list[str]] = None) -> list[int]:
        """Resolve group names to sorted list of IDs for igraph construction.

        Args:
            group_names: Optional list of group names. If None, returns all IDs.

        Returns:
            Sorted list of group IDs.
        """
        if group_names is None:
            return sorted(self._names_map.keys())
        target_ids = set()
        for name in group_names:
            gid = self._name_to_id.get(name.lower())
            if gid is not None:
                target_ids.add(gid)
        return sorted(target_ids)

    @classmethod
    def _from_edges(
        cls,
        edges_df: pd.DataFrame,
        names_map: dict[int, str],
        metric: str,
        cutoff: float,
        has_pvalue: bool,
        store: PairwiseStore,
        db_path: Optional[str] = None,
    ) -> "PairwiseGraph":
        """Build a PairwiseGraph directly from edge data (internal constructor)."""
        _require_kuzu()
        obj = object.__new__(cls)
        obj._store = store
        obj._metric = metric
        obj._cutoff = cutoff
        obj._names_map = names_map
        obj._name_to_id = {v: k for k, v in names_map.items()}

        if db_path:
            obj._db_dir = pathlib.Path(db_path)
            obj._owns_dir = False
        else:
            tmp_parent = pathlib.Path(tempfile.mkdtemp(prefix="dbretina_subgraph_"))
            obj._db_dir = tmp_parent / "db"
            obj._owns_dir = True
            obj._tmp_parent = tmp_parent

        obj._db = kuzu.Database(str(obj._db_dir))
        obj._conn = kuzu.Connection(obj._db)

        conn = obj._conn
        conn.execute(
            "CREATE NODE TABLE IF NOT EXISTS `Group`("
            "id UINT32, name STRING, PRIMARY KEY(id))"
        )
        rel_cols = (
            "shared_features UINT64, containment FLOAT, ochiai FLOAT, "
            "jaccard FLOAT, csi FLOAT, dice FLOAT, odds_ratio FLOAT"
        )
        if has_pvalue:
            rel_cols += ", pvalue DOUBLE"
        conn.execute(
            f'CREATE REL TABLE IF NOT EXISTS SIMILAR_TO('
            f'FROM `Group` TO `Group`, {rel_cols})'
        )

        tmp_parent_path = obj._db_dir.parent
        if names_map:
            node_csv = tmp_parent_path / "_nodes.csv"
            with open(node_csv, "w", newline="") as f:
                f.write("id,name\n")
                for gid, name in names_map.items():
                    escaped = name.replace('"', '""')
                    f.write(f'{gid},"{escaped}"\n')
            conn.execute(f'COPY `Group` FROM "{node_csv}" (HEADER=TRUE)')
            node_csv.unlink(missing_ok=True)

        if len(edges_df) > 0:
            edge_csv = tmp_parent_path / "_edges.csv"
            edges_df.to_csv(edge_csv, index=False, header=True)
            conn.execute(f'COPY SIMILAR_TO FROM "{edge_csv}" (HEADER=TRUE)')
            edge_csv.unlink(missing_ok=True)

        obj._num_nodes = len(names_map)
        obj._num_edges = len(edges_df)
        return obj

    def __repr__(self):
        return (
            f"PairwiseGraph(metric='{self._metric}', cutoff={self._cutoff}, "
            f"nodes={self._num_nodes}, edges={self._num_edges:,})"
        )

    def _repr_html_(self):
        """Rich HTML display for Jupyter notebooks."""
        try:
            cc = self.connected_components()
            n_components = len(set(cc.values()))
        except Exception:
            n_components = "?"

        try:
            deg = self.degree_distribution()
            avg_deg = sum(deg.values()) / len(deg) if deg else 0.0
            max_deg = max(deg.values()) if deg else 0
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
