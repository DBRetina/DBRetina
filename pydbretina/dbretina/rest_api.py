"""FastAPI REST server for DBRetina pairwise data + graph dashboard."""

import pathlib
import threading
from collections import OrderedDict
from typing import Optional

from fastapi import FastAPI, HTTPException, Query, Request
from fastapi.responses import JSONResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel

from .pairwise_store import PairwiseStore


# ── Graph cache ─────────────────────────────────────────────────

class _GraphCache:
    """LRU cache for PairwiseGraph instances keyed by (metric, cutoff)."""

    def __init__(self, store: PairwiseStore, max_size: int = 3):
        self._store = store
        self._max_size = max_size
        self._cache: OrderedDict = OrderedDict()
        self._lock = threading.Lock()

    def get(self, metric: str, cutoff: float):
        from .pairwise_graph import PairwiseGraph

        key = (metric, cutoff)
        with self._lock:
            if key in self._cache:
                self._cache.move_to_end(key)
                return self._cache[key]

        # Build outside lock (can be slow)
        graph = PairwiseGraph(self._store, metric, cutoff)

        with self._lock:
            if key in self._cache:
                graph.close()
                self._cache.move_to_end(key)
                return self._cache[key]
            self._cache[key] = graph
            while len(self._cache) > self._max_size:
                _, old = self._cache.popitem(last=False)
                old.close()
        return graph

    def close_all(self):
        with self._lock:
            for g in self._cache.values():
                g.close()
            self._cache.clear()


# ── App factory ─────────────────────────────────────────────────

def create_app(
    store: PairwiseStore,
    api_key: Optional[str] = None,
    dbri_path: Optional[str] = None,
    graph_metric: str = "ochiai",
    graph_cutoff: float = 20.0,
) -> FastAPI:
    """Create a FastAPI app wrapping a PairwiseStore.

    Args:
        store: Open PairwiseStore instance.
        api_key: Optional API key for authentication.
        dbri_path: Optional path to .dbri file for gene queries.
        graph_metric: Default metric for graph endpoints.
        graph_cutoff: Default cutoff for graph endpoints.
    """
    app = FastAPI(
        title="DBRetina API",
        version="2.0",
        description="REST API + interactive dashboard for DBRetina pairwise similarity data.",
    )
    app.state.store = store
    app.state.graph_cache = _GraphCache(store)
    app.state.dbri_path = dbri_path
    app.state.default_metric = graph_metric
    app.state.default_cutoff = graph_cutoff

    # ── CORS for dev mode ───────────────────────────────────────

    from fastapi.middleware.cors import CORSMiddleware
    app.add_middleware(
        CORSMiddleware,
        allow_origins=["*"],
        allow_methods=["*"],
        allow_headers=["*"],
    )

    # ── Auth middleware ──────────────────────────────────────────

    if api_key:
        @app.middleware("http")
        async def check_api_key(request: Request, call_next):
            path = request.url.path
            if path in ("/", "/docs", "/openapi.json", "/redoc") or not path.startswith("/api/"):
                return await call_next(request)
            key = request.headers.get("x-api-key") or request.query_params.get("api_key")
            if key != api_key:
                return JSONResponse(status_code=401, content={"detail": "Invalid API key"})
            return await call_next(request)

    # ── Helpers ──────────────────────────────────────────────────

    VALID_METRICS = ("containment", "ochiai", "jaccard", "csi", "dice", "odds_ratio", "pvalue")

    def validate_metric(metric: str):
        if metric not in VALID_METRICS:
            raise HTTPException(400, f"Unknown metric '{metric}'. Valid: {VALID_METRICS}")

    def get_graph(metric: Optional[str] = None, cutoff: Optional[float] = None):
        m = metric or app.state.default_metric
        c = cutoff if cutoff is not None else app.state.default_cutoff
        validate_metric(m)
        return app.state.graph_cache.get(m, c), m, c

    # ── Existing endpoints ──────────────────────────────────────

    @app.get("/api/v1/info")
    def get_info():
        """Dataset metadata."""
        return {
            "num_pairs": store.num_pairs,
            "num_groups": store.num_groups,
            "has_pvalue": store.has_pvalue,
            "has_genes": dbri_path is not None,
            "default_metric": app.state.default_metric,
            "default_cutoff": app.state.default_cutoff,
            "valid_metrics": list(VALID_METRICS),
            "manifest": store.manifest,
        }

    @app.get("/api/v1/groups")
    def list_groups(
        limit: int = Query(100, ge=1, le=10000),
        offset: int = Query(0, ge=0),
    ):
        """List all groups with IDs and names."""
        names = store.get_names_map()
        items = [{"id": gid, "name": name} for gid, name in sorted(names.items())]
        return {
            "total": len(items),
            "offset": offset,
            "limit": limit,
            "groups": items[offset : offset + limit],
        }

    @app.get("/api/v1/pairs")
    def get_pairs(
        metric: str = Query(..., description="Metric to filter on"),
        cutoff: float = Query(0.0, ge=0),
        limit: int = Query(1000, ge=1, le=100000),
        offset: int = Query(0, ge=0),
    ):
        """Query pairs filtered by metric and cutoff, with pagination."""
        validate_metric(metric)
        df = store.to_pandas(metric=metric, cutoff=cutoff, limit=offset + limit)
        if offset > 0:
            df = df.iloc[offset:]
        names = store.get_names_map()
        df["group_1_name"] = df["group_1_id"].map(names)
        df["group_2_name"] = df["group_2_id"].map(names)
        records = df.to_dict(orient="records")
        return {
            "offset": offset,
            "limit": limit,
            "count": len(records),
            "pairs": records,
        }

    @app.get("/api/v1/groups/{group_name}/pairs")
    def get_group_pairs(
        group_name: str,
        metric: Optional[str] = Query(None),
        cutoff: float = Query(0.0, ge=0),
        limit: int = Query(1000, ge=1, le=100000),
    ):
        """Get all pairs involving a specific group."""
        if metric:
            validate_metric(metric)
        try:
            reader = store.query_group(group_name, metric=metric, cutoff=cutoff)
        except KeyError:
            raise HTTPException(404, f"Group not found: {group_name}")

        rows = []
        names = store.get_names_map()
        for batch in reader:
            chunk = batch.to_pydict()
            for i in range(len(batch)):
                row = {col: chunk[col][i] for col in chunk}
                row["group_1_name"] = names.get(row.get("group_1_id"), "")
                row["group_2_name"] = names.get(row.get("group_2_id"), "")
                rows.append(row)
                if len(rows) >= limit:
                    break
            if len(rows) >= limit:
                break

        return {"group": group_name, "count": len(rows), "pairs": rows}

    @app.get("/api/v1/statistics")
    def get_statistics():
        """Metric distribution statistics."""
        return store.get_statistics()

    @app.get("/api/v1/statistics/{metric}")
    def get_metric_summary(metric: str):
        """Summary statistics for a single metric."""
        validate_metric(metric)
        return store.metric_summary(metric)

    class SQLQuery(BaseModel):
        query: str

    @app.post("/api/v1/sql")
    def execute_sql(body: SQLQuery):
        """Execute SQL against the pairs view."""
        try:
            result = store.sql(body.query)
            df = result.fetchdf()
            return {
                "columns": list(df.columns),
                "row_count": len(df),
                "rows": df.to_dict(orient="records"),
            }
        except Exception as e:
            raise HTTPException(400, f"SQL error: {e}")

    @app.get("/api/v1/top")
    def top_pairs(
        metric: str = Query(...),
        n: int = Query(20, ge=1, le=1000),
    ):
        """Top N pairs by a given metric."""
        validate_metric(metric)
        df = store.top_pairs(metric, n=n)
        return {"metric": metric, "count": len(df), "pairs": df.to_dict(orient="records")}

    # ── New: Search endpoints ───────────────────────────────────

    @app.get("/api/v1/groups/search")
    def search_groups(
        q: str = Query(..., min_length=1, description="Search pattern"),
        limit: int = Query(20, ge=1, le=100),
    ):
        """Search groups by name substring (case-insensitive)."""
        matches = store.search_groups(q)
        items = [{"id": gid, "name": name} for gid, name in sorted(matches.items())]
        return {"query": q, "count": len(items[:limit]), "matches": items[:limit]}

    @app.get("/api/v1/features/search")
    def search_features(
        q: str = Query(..., min_length=1, description="Gene/feature name"),
        limit: int = Query(20, ge=1, le=200),
    ):
        """Find groups containing a gene/feature."""
        if not dbri_path:
            raise HTTPException(400, "Gene index not available (start server with -i flag)")
        try:
            groups = store.search_by_feature(q)
        except Exception as e:
            raise HTTPException(400, str(e))
        return {"feature": q, "count": len(groups[:limit]), "groups": groups[:limit]}

    @app.get("/api/v1/shared-features")
    def shared_features(
        group_a: str = Query(...),
        group_b: str = Query(...),
    ):
        """Get shared genes/features between two groups."""
        if not dbri_path:
            raise HTTPException(400, "Gene index not available (start server with -i flag)")
        try:
            features = store.shared_features(group_a, group_b)
        except Exception as e:
            raise HTTPException(400, str(e))
        return {
            "group_a": group_a,
            "group_b": group_b,
            "count": len(features),
            "features": sorted(features),
        }

    # ── New: Graph endpoints ────────────────────────────────────

    def _graph_data_response(graph, metric_name, limit=5000):
        """Convert a PairwiseGraph to the standard node/edge JSON format."""
        import igraph as ig

        # Get edges
        edges_df = graph.cypher(
            f'MATCH (a:`Group`)-[r:SIMILAR_TO]->(b:`Group`) '
            f'RETURN a.id AS src, b.id AS dst, r.{metric_name} AS weight, '
            f'r.shared_features AS shared_features'
        )
        names_map = graph._names_map

        # Build igraph for degree + pagerank + communities
        all_ids = sorted(names_map.keys())
        id_idx = {gid: i for i, gid in enumerate(all_ids)}
        g = ig.Graph(n=len(all_ids), directed=False)

        edge_list = []
        weights = []
        edge_data = []
        for _, row in edges_df.iterrows():
            src, dst = int(row["src"]), int(row["dst"])
            if src in id_idx and dst in id_idx:
                edge_list.append((id_idx[src], id_idx[dst]))
                w = float(row["weight"])
                weights.append(w)
                edge_data.append({
                    "source": str(src), "target": str(dst),
                    "weight": round(w, 2),
                    "shared_features": int(row["shared_features"]),
                })

        if edge_list:
            g.add_edges(edge_list)
            g.es["weight"] = weights

        degrees = g.degree()
        try:
            pr = g.pagerank(weights="weight")
        except Exception:
            pr = [0.0] * g.vcount()

        try:
            communities = g.community_leiden(weights="weight", objective_function="modularity")
            membership = communities.membership
        except Exception:
            membership = [0] * g.vcount()

        # Filter to connected nodes only
        connected_ids = set()
        for ed in edge_data:
            connected_ids.add(int(ed["source"]))
            connected_ids.add(int(ed["target"]))

        nodes = []
        for gid in all_ids:
            if gid not in connected_ids:
                continue
            idx = id_idx[gid]
            nodes.append({
                "id": str(gid),
                "label": names_map[gid],
                "degree": degrees[idx],
                "community": membership[idx],
                "pagerank": round(pr[idx], 6),
            })

        truncated = False
        if len(nodes) > limit:
            # Keep top-degree nodes
            nodes.sort(key=lambda n: n["degree"], reverse=True)
            keep_ids = {n["id"] for n in nodes[:limit]}
            nodes = nodes[:limit]
            edge_data = [e for e in edge_data if e["source"] in keep_ids and e["target"] in keep_ids]
            truncated = True

        return {
            "nodes": nodes,
            "edges": edge_data,
            "meta": {
                "total_nodes": len(connected_ids),
                "total_edges": len(edges_df),
                "returned_nodes": len(nodes),
                "returned_edges": len(edge_data),
                "truncated": truncated,
                "metric": metric_name,
                "cutoff": graph._cutoff,
            },
        }

    @app.get("/api/v1/graph/data")
    def get_graph_data(
        metric: Optional[str] = Query(None),
        cutoff: Optional[float] = Query(None),
        limit: int = Query(5000, ge=1, le=50000),
    ):
        """Get graph nodes and edges for visualization."""
        graph, m, c = get_graph(metric, cutoff)
        return _graph_data_response(graph, m, limit=limit)

    @app.get("/api/v1/graph/neighborhood")
    def get_neighborhood(
        group: str = Query(..., description="Target group name"),
        metric: Optional[str] = Query(None),
        cutoff: Optional[float] = Query(None),
        hops: int = Query(2, ge=1, le=5),
    ):
        """Get subgraph around a target group."""
        graph, m, c = get_graph(metric, cutoff)
        try:
            sub = graph.subgraph_around(group, hops=hops)
        except KeyError:
            raise HTTPException(404, f"Group not found: {group}")
        try:
            return _graph_data_response(sub, m, limit=10000)
        finally:
            sub.close()

    @app.get("/api/v1/graph/communities")
    def get_communities(
        metric: Optional[str] = Query(None),
        cutoff: Optional[float] = Query(None),
        method: str = Query("leiden"),
    ):
        """Get community assignments for all nodes."""
        graph, m, c = get_graph(metric, cutoff)
        communities = graph.community_detection(method=method)
        from collections import Counter
        sizes = dict(Counter(communities.values()).most_common())
        return {
            "communities": communities,
            "num_communities": len(sizes),
            "sizes": sizes,
        }

    @app.get("/api/v1/graph/layout")
    def get_layout(
        metric: Optional[str] = Query(None),
        cutoff: Optional[float] = Query(None),
        algorithm: str = Query("fr", description="Layout algorithm: fr, drl, kk, circle"),
    ):
        """Get pre-computed node positions for graph layout."""
        import igraph as ig

        graph, m, c = get_graph(metric, cutoff)

        edges_df = graph.cypher(
            f'MATCH (a:`Group`)-[r:SIMILAR_TO]->(b:`Group`) '
            f'RETURN a.id AS src, b.id AS dst, r.{m} AS weight'
        )
        names_map = graph._names_map
        all_ids = sorted(names_map.keys())
        id_idx = {gid: i for i, gid in enumerate(all_ids)}

        g = ig.Graph(n=len(all_ids), directed=False)
        edge_list = []
        weights = []
        for _, row in edges_df.iterrows():
            src, dst = int(row["src"]), int(row["dst"])
            if src in id_idx and dst in id_idx:
                edge_list.append((id_idx[src], id_idx[dst]))
                weights.append(float(row["weight"]))
        if edge_list:
            g.add_edges(edge_list)
            g.es["weight"] = weights

        algo_map = {
            "fr": "fruchterman_reingold",
            "drl": "drl",
            "kk": "kamada_kawai",
            "circle": "circle",
        }
        layout_name = algo_map.get(algorithm, "fruchterman_reingold")

        if g.vcount() > 500 and layout_name == "fruchterman_reingold":
            layout_name = "drl"

        layout = g.layout(layout_name)

        positions = {}
        for gid in all_ids:
            idx = id_idx[gid]
            if g.degree(idx) > 0:
                positions[str(gid)] = [round(layout[idx][0], 4), round(layout[idx][1], 4)]

        return {"algorithm": layout_name, "positions": positions}

    class CypherQuery(BaseModel):
        query: str
        metric: Optional[str] = None
        cutoff: Optional[float] = None

    @app.post("/api/v1/cypher")
    def execute_cypher(body: CypherQuery):
        """Execute a Cypher query on the graph database."""
        graph, m, c = get_graph(body.metric, body.cutoff)
        try:
            df = graph.cypher(body.query)
            return {
                "columns": list(df.columns),
                "row_count": len(df),
                "rows": df.to_dict(orient="records"),
            }
        except Exception as e:
            raise HTTPException(400, f"Cypher error: {e}")

    # ── Static file serving for dashboard ───────────────────────

    dashboard_dir = pathlib.Path(__file__).parent / "dashboard_dist"
    if dashboard_dir.is_dir():
        app.mount("/", StaticFiles(directory=str(dashboard_dir), html=True), name="dashboard")

    # ── Shutdown hook ───────────────────────────────────────────

    @app.on_event("shutdown")
    def shutdown():
        app.state.graph_cache.close_all()

    return app
