"""FastAPI REST server for DBRetina pairwise data + graph dashboard."""

import asyncio
import math
import pathlib
import threading
import time
from collections import OrderedDict
from functools import wraps
from typing import Any, Callable, Literal, Optional

from fastapi import FastAPI, HTTPException, Query, Request
from fastapi.responses import JSONResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel, Field

from .pairwise_store import PairwiseStore
from .api_errors import (
    DBRetinaAPIError,
    DataNotFoundError,
    ValidationError,
    QueryTimeoutError,
    DataTooLargeError,
    QuerySyntaxError,
    UnsafeQueryError,
    FeatureNotAvailableError,
    AlgorithmError,
    validate_sql_safety,
    validate_metric as _validate_metric,
    validate_cutoff as _validate_cutoff,
)


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

    # ── Sandbox the SQL surface: lock DuckDB to in-memory data, no filesystem ──
    # Untrusted queries (POST /api/v1/sql) run against this connection; disabling
    # external access blocks read_text/read_csv/COPY/ATTACH etc. — defense that does
    # not depend on a keyword denylist.
    store.harden_readonly()

    # ── CORS (restrict to the dashboard dev origin; never wildcard) ──────
    # The bundled dashboard is served same-origin (no CORS needed); only the Vite
    # dev server is cross-origin.
    from fastapi.middleware.cors import CORSMiddleware
    app.add_middleware(
        CORSMiddleware,
        allow_origins=["http://localhost:5173", "http://127.0.0.1:5173"],
        allow_credentials=False,
        allow_methods=["GET", "POST"],
        allow_headers=["x-api-key", "content-type"],
    )

    # ── Exception handlers ────────────────────────────────────────

    @app.exception_handler(DBRetinaAPIError)
    async def dbretina_error_handler(request: Request, exc: DBRetinaAPIError):
        """Handle all custom DBRetina exceptions."""
        return JSONResponse(
            status_code=exc.status_code,
            content=exc.to_dict(),
        )

    @app.exception_handler(Exception)
    async def general_error_handler(request: Request, exc: Exception):
        """Handle unexpected exceptions with a generic error response."""
        # Log the error for debugging (in production, use proper logging)
        import traceback
        traceback.print_exc()
        return JSONResponse(
            status_code=500,
            content={
                "error": True,
                "error_code": "INTERNAL_ERROR",
                "detail": "An unexpected error occurred. Please try again.",
            },
        )

    # ── Auth middleware ──────────────────────────────────────────

    if api_key:
        @app.middleware("http")
        async def check_api_key(request: Request, call_next):
            path = request.url.path
            if path in ("/", "/docs", "/openapi.json", "/redoc") or not path.startswith("/api/"):
                return await call_next(request)
            key = request.headers.get("x-api-key")
            if key != api_key:
                return JSONResponse(status_code=401, content={"detail": "Invalid API key"})
            return await call_next(request)

    # ── Helpers ──────────────────────────────────────────────────

    DEFAULT_TIMEOUT = 30.0  # seconds
    MAX_PAIRS_RESPONSE = 100000
    MAX_GRAPH_NODES = 50000

    # Metric-specific cutoff ranges
    METRIC_CUTOFF_RANGES = {
        "containment": (0.0, 100.0),
        "ochiai": (0.0, 100.0),
        "jaccard": (0.0, 100.0),
        "csi": (0.0, 100.0),
        "dice": (0.0, 100.0),
        "odds_ratio": (0.0, float("inf")),  # odds ratio can be very large
        "pvalue": (0.0, 1.0),
    }

    def validate_metric(metric: str):
        """Validate metric name (respects this dataset's available metrics, e.g. pvalue)."""
        _validate_metric(metric, store.available_metrics)

    def validate_cutoff(cutoff: float, metric: str = "ochiai"):
        """Validate cutoff value based on metric type."""
        # Reject non-finite cutoffs (inf / -inf / NaN) for EVERY metric: even
        # metrics with an infinite range max would otherwise pass the bounds
        # check below and be interpolated into SQL as the bare token `inf`/`nan`,
        # which DuckDB parses as a (missing) column -> uncaught 500.
        if not math.isfinite(cutoff):
            raise ValidationError(
                detail=f"Cutoff must be a finite number, got {cutoff}",
                field="cutoff",
                value=cutoff,
            )
        min_val, max_val = METRIC_CUTOFF_RANGES.get(metric, (0.0, 100.0))
        # For infinity max, just check minimum
        if max_val == float("inf"):
            if cutoff < min_val:
                raise ValidationError(
                    detail=f"Cutoff for {metric} must be >= {min_val}",
                    field="cutoff",
                    value=cutoff,
                )
        else:
            _validate_cutoff(cutoff, min_val, max_val)

    def check_data_size(size: int, max_size: int, resource: str = "data"):
        """Check if requested data size is within limits."""
        if size > max_size:
            raise DataTooLargeError(
                detail=f"Requested {resource} ({size:,} items) exceeds limit ({max_size:,})",
                requested_size=size,
                max_size=max_size,
                suggestion=f"Use pagination or filters to reduce result size",
            )

    def get_graph(metric: Optional[str] = None, cutoff: Optional[float] = None):
        m = metric or app.state.default_metric
        c = cutoff if cutoff is not None else app.state.default_cutoff
        validate_metric(m)
        validate_cutoff(c, m)
        return app.state.graph_cache.get(m, c), m, c

    # ── Existing endpoints ──────────────────────────────────────

    # ── Health & Info endpoints ────────────────────────────────

    @app.get("/api/v1/health")
    def health_check():
        """Health check endpoint for monitoring."""
        try:
            # Quick connectivity check
            _ = store.num_pairs
            status = "healthy"
            details = {
                "store": "connected",
                "cache_size": len(app.state.graph_cache._cache),
            }
        except Exception as e:
            status = "unhealthy"
            details = {"error": str(e)}

        return {
            "status": status,
            "timestamp": time.time(),
            "version": "2.0",
            "details": details,
        }

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
            "valid_metrics": list(store.available_metrics),
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
        validate_cutoff(cutoff, metric)
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
        # cutoff is interpolated into SQL by query_group; reject non-finite /
        # out-of-range values here (this endpoint never validated it before, so a
        # cutoff of inf reached DuckDB as a bare token -> 500).
        validate_cutoff(cutoff, metric or app.state.default_metric)
        try:
            reader = store.query_group(group_name, metric=metric, cutoff=cutoff)
        except KeyError:
            raise DataNotFoundError(
                detail=f"Group not found: {group_name}",
                resource_type="group",
                resource_id=group_name,
            )

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

    @app.get("/api/v1/groups/{group_name}/genes")
    def get_group_genes(
        group_name: str,
        limit: int = Query(50, ge=1, le=500),
    ):
        """Get all genes/features associated with a specific group."""
        if not store._dbri_path:
            raise HTTPException(
                status_code=400,
                detail="Gene data not available. Server was started without -i/--index-prefix.",
            )
        try:
            genes = store.get_group_genes(group_name)
            return {
                "group": group_name,
                "count": len(genes),
                "genes": sorted(genes)[:limit],
            }
        except KeyError:
            raise DataNotFoundError(
                detail=f"Group not found: {group_name}",
                resource_type="group",
                resource_id=group_name,
            )

    @app.get("/api/v1/groups/{group_name}/metric-profile")
    def get_group_metric_profile(group_name: str):
        """Get aggregated metric profile (AVG, MAX, COUNT) for a group."""
        try:
            metrics = store.group_metric_profile(group_name)
        except KeyError:
            raise DataNotFoundError(
                detail=f"Group not found: {group_name}",
                resource_type="group",
                resource_id=group_name,
            )
        except Exception as e:
            raise DBRetinaAPIError(detail=f"Metric profile failed: {e}")
        return {"group": group_name, "metrics": metrics}

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

    # ── Advanced Filtering ─────────────────────────────────────────

    class FilterCondition(BaseModel):
        """A single filter condition."""
        metric: str
        operator: str  # ">=", "<=", ">", "<", "==", "!=", "between"
        value: float | list[float]  # single value or [min, max] for "between"

    class FilterRequest(BaseModel):
        """Multi-metric filter request."""
        filters: list[FilterCondition]
        logic: str = "AND"  # "AND" or "OR"
        limit: int = 5000
        offset: int = 0

    @app.post("/api/v1/pairs/filter")
    def filter_pairs(body: FilterRequest):
        """Filter pairs using multiple metric conditions.

        Allows combining multiple filter conditions with AND/OR logic.
        This is more powerful than the simple GET /pairs endpoint which
        only supports a single metric/cutoff filter.

        Example:
            {
                "filters": [
                    {"metric": "ochiai", "operator": ">=", "value": 50},
                    {"metric": "pvalue", "operator": "<=", "value": 0.05}
                ],
                "logic": "AND",
                "limit": 1000
            }
        """
        if not body.filters:
            raise ValidationError(detail="At least one filter condition is required")

        if body.logic not in ("AND", "OR"):
            raise ValidationError(
                detail="Logic must be 'AND' or 'OR'",
                field="logic",
                value=body.logic,
                allowed_values=["AND", "OR"],
            )

        # Validate all metrics and operators
        valid_operators = {">=", "<=", ">", "<", "==", "!=", "between"}
        for filt in body.filters:
            validate_metric(filt.metric)
            if filt.operator not in valid_operators:
                raise ValidationError(
                    detail=f"Invalid operator '{filt.operator}'",
                    field="operator",
                    value=filt.operator,
                    allowed_values=list(valid_operators),
                )
            if filt.operator == "between":
                if not isinstance(filt.value, list) or len(filt.value) != 2:
                    raise ValidationError(
                        detail="'between' operator requires [min, max] value",
                        field="value",
                    )
            # Reject non-finite numeric values: they are interpolated into the SQL
            # WHERE clause below and would reach DuckDB as the bare token inf/nan
            # (parsed as a missing column -> uncaught 500).
            _vals = filt.value if isinstance(filt.value, list) else [filt.value]
            for _v in _vals:
                if isinstance(_v, float) and not math.isfinite(_v):
                    raise ValidationError(
                        detail=f"Filter value must be finite, got {_v}",
                        field="value",
                        value=_v,
                    )

        # Build SQL WHERE clause
        conditions = []
        for filt in body.filters:
            if filt.operator == "between":
                min_val, max_val = filt.value
                conditions.append(f"{filt.metric} >= {min_val} AND {filt.metric} <= {max_val}")
            else:
                conditions.append(f"{filt.metric} {filt.operator} {filt.value}")

        connector = f" {body.logic} "
        where_clause = connector.join(f"({c})" for c in conditions)

        # Execute query
        query = f"SELECT * FROM pairs WHERE {where_clause} LIMIT {body.limit} OFFSET {body.offset}"
        try:
            result = store.sql(query)
            df = result.fetchdf()
        except Exception as e:
            raise QuerySyntaxError(
                detail=f"Filter query failed: {e}",
                query_type="sql",
            )

        # Add names
        names = store.get_names_map()
        if "group_1_id" in df.columns:
            df["group_1_name"] = df["group_1_id"].map(names)
        if "group_2_id" in df.columns:
            df["group_2_name"] = df["group_2_id"].map(names)

        return {
            "count": len(df),
            "filters": [f.model_dump() for f in body.filters],
            "logic": body.logic,
            "limit": body.limit,
            "offset": body.offset,
            "pairs": df.to_dict(orient="records"),
        }

    @app.post("/api/v1/sql")
    def execute_sql(body: SQLQuery):
        """Execute SQL against the pairs view (read-only)."""
        # Validate query safety
        validate_sql_safety(body.query)

        try:
            result = store.sql(body.query)
            df = result.fetchdf()

            # Check result size
            if len(df) > MAX_PAIRS_RESPONSE:
                raise DataTooLargeError(
                    detail=f"Query returned too many rows ({len(df):,})",
                    requested_size=len(df),
                    max_size=MAX_PAIRS_RESPONSE,
                    suggestion="Add LIMIT clause to your query",
                )

            return {
                "columns": list(df.columns),
                "row_count": len(df),
                "rows": df.to_dict(orient="records"),
            }
        except (DataTooLargeError, UnsafeQueryError):
            raise
        except Exception as e:
            raise QuerySyntaxError(
                detail=f"SQL error: {e}",
                query_type="sql",
                query=body.query,
            )

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

    @app.get("/api/v1/genes/search")
    def search_genes(
        q: str = Query(..., min_length=1, description="Gene name query"),
        limit: int = Query(20, ge=1, le=100),
    ):
        """Search gene names by prefix/substring match for autocomplete."""
        if not dbri_path:
            raise FeatureNotAvailableError(
                detail="Gene index not available",
                feature="gene_search",
                requirement="Start server with -i flag to specify index file",
            )
        try:
            matches = store.search_gene_names(q, limit=limit)
        except Exception as e:
            raise DBRetinaAPIError(detail=f"Gene search failed: {e}")
        return {"query": q, "genes": matches}

    @app.get("/api/v1/features/search")
    def search_features(
        q: str = Query(..., min_length=1, description="Gene/feature name"),
        limit: int = Query(20, ge=1, le=200),
    ):
        """Find groups containing a gene/feature."""
        if not dbri_path:
            raise FeatureNotAvailableError(
                detail="Gene index not available",
                feature="gene_search",
                requirement="Start server with -i flag to specify index file",
            )
        try:
            groups = store.search_by_feature(q)
        except Exception as e:
            raise DataNotFoundError(
                detail=f"Feature search failed: {e}",
                resource_type="feature",
                resource_id=q,
            )
        return {"feature": q, "count": len(groups[:limit]), "groups": groups[:limit]}

    @app.get("/api/v1/shared-features")
    def shared_features(
        group_a: str = Query(...),
        group_b: str = Query(...),
    ):
        """Get shared genes/features between two groups."""
        if not dbri_path:
            raise FeatureNotAvailableError(
                detail="Gene index not available",
                feature="shared_features",
                requirement="Start server with -i flag to specify index file",
            )
        try:
            features = store.shared_features(group_a, group_b)
        except KeyError as e:
            raise DataNotFoundError(
                detail=f"Group not found: {e}",
                resource_type="group",
            )
        except Exception as e:
            raise DBRetinaAPIError(detail=f"Failed to get shared features: {e}")
        return {
            "group_a": group_a,
            "group_b": group_b,
            "count": len(features),
            "features": sorted(features),
        }

    # ── Gene Importance & Analysis ─────────────────────────────────

    # Allowed gene-importance scoring methods (mirrors gene_importance.hub_genes
    # dispatch). Declared as a Literal so an unknown method is a clean 422 at the
    # request boundary instead of silently falling through to projection.
    GeneMethod = Literal["hypergraph", "edge_weighted", "projection"]

    class HubGenesRequest(BaseModel):
        """Request for hub genes analysis."""
        group_name: str
        method: GeneMethod = "hypergraph"
        # hops feeds igraph neighborhood order (must be non-negative); bound it
        # like the sibling neighborhood endpoint (ge=1, le=5).
        hops: int = Field(2, ge=1, le=5)
        # top_n slices the ranked list; a negative value silently truncates.
        top_n: int = Field(30, ge=1)
        metric: str = "ochiai"
        cutoff: float = 20.0

    class ExplainPairRequest(BaseModel):
        """Request to explain gene connection between two groups."""
        group_a: str
        group_b: str
        method: GeneMethod = "hypergraph"

    class ClusterGenesRequest(BaseModel):
        """Request for cluster gene analysis."""
        node_names: list[str]
        method: GeneMethod = "hypergraph"
        # top_n slices the ranked list; a negative value silently truncates.
        top_n: int = Field(50, ge=1)

    def _get_gene_importance():
        """Get or create GeneImportance instance."""
        if not dbri_path:
            raise FeatureNotAvailableError(
                detail="Gene index not available",
                feature="gene_importance",
                requirement="Start server with -i flag to specify index file",
            )
        from .gene_importance import GeneImportance
        return GeneImportance(store)

    @app.post("/api/v1/genes/hub-genes")
    def get_hub_genes(body: HubGenesRequest):
        """Find the most important genes for a disease/group.

        Uses three scoring methods:
        - **hypergraph**: TF-IDF style (local enrichment × inverse global frequency)
        - **edge_weighted**: Sum of edge weights per gene
        - **projection**: PageRank on gene co-occurrence graph
        """
        gi = _get_gene_importance()
        # The edge_weighted method is the only one that consumes ``metric`` and
        # ``cutoff``; reject an unknown metric or a non-finite / out-of-range
        # cutoff up front as a 400. (cutoff is interpolated into SQL by
        # edge_weighted_scores, so inf/NaN would otherwise 500 and a negative
        # value would be a wrong-but-200.)
        if body.method == "edge_weighted":
            _validate_metric(body.metric, store.available_metrics)
            validate_cutoff(body.cutoff, body.metric)
        try:
            df = gi.hub_genes(
                group_name=body.group_name,
                hops=body.hops,
                method=body.method,
                top_n=body.top_n,
                metric=body.metric,
                cutoff=body.cutoff,
            )
            genes = df.to_dict(orient="records")
            return {
                "group": body.group_name,
                "method": body.method,
                "hops": body.hops,
                "genes": genes,
            }
        except KeyError as e:
            raise DataNotFoundError(
                detail=f"Group not found: {e}",
                resource_type="group",
                resource_id=body.group_name,
            )
        except ValueError as e:
            # Bad user input the model couldn't catch (e.g. an unknown metric
            # reaching the scorer) is a client error, not a server fault.
            raise ValidationError(detail=str(e))
        except Exception as e:
            raise DBRetinaAPIError(detail=f"Hub genes analysis failed: {e}")

    @app.post("/api/v1/genes/explain-pair")
    def explain_pair(body: ExplainPairRequest):
        """Rank shared genes between two diseases by importance.

        Returns genes sorted by their contribution to the similarity
        between the two groups, with global frequency and specificity info.
        """
        gi = _get_gene_importance()
        try:
            df = gi.explain_pair(
                group_a=body.group_a,
                group_b=body.group_b,
                method=body.method,
            )
            genes = df.to_dict(orient="records")
            return {
                "group_a": body.group_a,
                "group_b": body.group_b,
                "method": body.method,
                "gene_count": len(genes),
                "genes": genes,
            }
        except KeyError as e:
            raise DataNotFoundError(
                detail=f"Group not found: {e}",
                resource_type="group",
            )
        except Exception as e:
            raise DBRetinaAPIError(detail=f"Explain pair failed: {e}")

    @app.get("/api/v1/genes/{gene_name}/statistics")
    def get_gene_statistics(gene_name: str):
        """Get prevalence statistics for a specific gene.

        Returns how many groups contain this gene and what percentage
        of the total groups that represents.
        """
        if not dbri_path:
            raise FeatureNotAvailableError(
                detail="Gene index not available",
                feature="gene_statistics",
                requirement="Start server with -i flag to specify index file",
            )
        try:
            gene_index = store._get_gene_index()
            groups_with_gene = gene_index.get(gene_name, set())

            if not groups_with_gene:
                # Try case-insensitive search
                gene_lower = gene_name.lower()
                for g, group_set in gene_index.items():
                    if g.lower() == gene_lower:
                        groups_with_gene = group_set
                        gene_name = g  # Use the actual case
                        break

            if not groups_with_gene:
                raise DataNotFoundError(
                    detail=f"Gene not found: {gene_name}",
                    resource_type="gene",
                    resource_id=gene_name,
                )

            total_groups = store.num_groups
            group_count = len(groups_with_gene)
            prevalence = (group_count / total_groups * 100) if total_groups > 0 else 0

            return {
                "gene": gene_name,
                "group_count": group_count,
                "total_groups": total_groups,
                "prevalence_percent": round(prevalence, 2),
            }
        except DataNotFoundError:
            raise
        except Exception as e:
            raise DBRetinaAPIError(detail=f"Gene statistics failed: {e}")

    @app.get("/api/v1/genes/{gene_name}/groups")
    def get_gene_groups(
        gene_name: str,
        limit: int = Query(50, ge=1, le=500),
    ):
        """Get all groups containing a specific gene.

        Returns groups sorted by their degree (connectivity) in the
        current graph, helping identify which diseases are most
        connected that involve this gene.
        """
        if not dbri_path:
            raise FeatureNotAvailableError(
                detail="Gene index not available",
                feature="gene_groups",
                requirement="Start server with -i flag to specify index file",
            )
        try:
            gene_index = store._get_gene_index()
            groups_with_gene = gene_index.get(gene_name, set())

            if not groups_with_gene:
                # Try case-insensitive search
                gene_lower = gene_name.lower()
                for g, group_set in gene_index.items():
                    if g.lower() == gene_lower:
                        groups_with_gene = group_set
                        gene_name = g
                        break

            if not groups_with_gene:
                raise DataNotFoundError(
                    detail=f"Gene not found: {gene_name}",
                    resource_type="gene",
                    resource_id=gene_name,
                )

            # Get group details
            names_map = store.get_names_map()
            groups = []
            for group_name_lower in groups_with_gene:
                # Find the group ID
                for gid, name in names_map.items():
                    if name.lower() == group_name_lower.lower():
                        groups.append({
                            "id": gid,
                            "name": name,
                        })
                        break

            # Sort by name for now (could add degree if graph is available)
            groups.sort(key=lambda x: x["name"])

            return {
                "gene": gene_name,
                "group_count": len(groups),
                "groups": groups[:limit],
            }
        except DataNotFoundError:
            raise
        except Exception as e:
            raise DBRetinaAPIError(detail=f"Gene groups failed: {e}")

    @app.post("/api/v1/genes/cluster-analysis")
    def analyze_cluster_genes(body: ClusterGenesRequest):
        """Analyze which genes define a cluster of diseases.

        Given a list of disease/group names, finds genes that are
        enriched within this cluster compared to the global dataset.
        """
        if not dbri_path:
            raise FeatureNotAvailableError(
                detail="Gene index not available",
                feature="cluster_analysis",
                requirement="Start server with -i flag to specify index file",
            )
        try:
            gene_sets = store._load_gene_sets()
            gene_index = store._get_gene_index()
            total_groups = store.num_groups

            # Collect genes from cluster
            cluster_genes: dict[str, int] = {}
            valid_nodes = 0
            for node_name in body.node_names:
                node_key = node_name.lower()
                if node_key in gene_sets:
                    valid_nodes += 1
                    for gene in gene_sets[node_key]:
                        cluster_genes[gene] = cluster_genes.get(gene, 0) + 1

            if valid_nodes == 0:
                raise DataNotFoundError(
                    detail="No valid groups found in cluster",
                    resource_type="cluster",
                )

            # Score genes by enrichment
            import math
            scored_genes = []
            for gene, in_cluster_count in cluster_genes.items():
                global_count = len(gene_index.get(gene, set()))
                if global_count > 0:
                    # TF-IDF style scoring
                    tf = in_cluster_count / valid_nodes
                    idf = math.log(total_groups / global_count) if global_count < total_groups else 0
                    score = tf * idf
                    scored_genes.append({
                        "gene": gene,
                        "score": round(score, 6),
                        "in_cluster_count": in_cluster_count,
                        "cluster_size": valid_nodes,
                        "global_count": global_count,
                    })

            # Sort by score descending
            scored_genes.sort(key=lambda x: -x["score"])

            return {
                "cluster_size": valid_nodes,
                "total_genes": len(scored_genes),
                "method": body.method,
                "genes": scored_genes[:body.top_n],
            }
        except DataNotFoundError:
            raise
        except Exception as e:
            raise DBRetinaAPIError(detail=f"Cluster analysis failed: {e}")

    # ── New: Graph endpoints ────────────────────────────────────

    def _graph_data_response(graph, metric_name, limit=5000):
        """Convert a PairwiseGraph to the standard node/edge JSON format."""
        # Reuse the graph's already-built igraph (vertex i has gid == idx_to_gid[i]).
        g = graph.to_igraph()
        names_map = graph._names_map
        id_idx = graph._gid_to_idx
        all_ids = graph._idx_to_gid

        # The graph's edge weights ARE this metric; read them back per edge.
        weights = g.es["weight"] if g.ecount() > 0 else None
        # Every metric the dataset carries is also stored as an edge attribute
        # (see PairwiseGraph). Expose each one so the client can filter on any
        # metric without re-querying. ``available_metrics`` excludes pvalue when
        # it wasn't computed, so no pvalue field appears on those datasets.
        edge_attrs = set(g.es.attributes()) if g.ecount() > 0 else set()
        metric_cols = [m for m in store.available_metrics if m in edge_attrs]
        edge_data = []
        for e in g.es:
            src = all_ids[e.source]
            dst = all_ids[e.target]
            ed = {
                "source": str(src), "target": str(dst),
                "weight": round(float(e["weight"]), 2),
                "shared_features": int(e["shared_features"]),
            }
            for m in metric_cols:
                ed[m] = round(float(e[m]), 6)
            edge_data.append(ed)

        degrees = g.degree()
        try:
            pr = g.pagerank(weights=weights)
        except Exception:
            pr = [0.0] * g.vcount()

        try:
            communities = g.community_leiden(weights=weights, objective_function="modularity")
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
                "total_edges": g.ecount(),
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
            raise DataNotFoundError(
                detail=f"Group not found: {group}",
                resource_type="group",
                resource_id=group,
            )
        try:
            return _graph_data_response(sub, m, limit=10000)
        finally:
            sub.close()

    @app.get("/api/v1/graph/shortest-path")
    def get_shortest_path(
        source: str = Query(..., description="Source group name"),
        target: str = Query(..., description="Target group name"),
        metric: Optional[str] = Query(None),
        cutoff: Optional[float] = Query(None),
    ):
        """Find the shortest path between two groups in the graph.

        Returns the path as a list of group names, along with path length
        and connectivity status. If no path exists, returns connected=False.
        """
        graph, m, c = get_graph(metric, cutoff)
        try:
            result = graph.shortest_path_full(source, target)

            # If connected and we have gene data, get shared genes along path
            if result.get("connected") and dbri_path and len(result.get("path_nodes", [])) >= 2:
                path_nodes = result["path_nodes"]
                shared_genes_along_path = []

                # Get shared genes between consecutive pairs in path
                for i in range(len(path_nodes) - 1):
                    try:
                        shared = store.shared_features(path_nodes[i], path_nodes[i + 1])
                        shared_genes_along_path.append({
                            "from": path_nodes[i],
                            "to": path_nodes[i + 1],
                            "shared_count": len(shared),
                            "genes": sorted(shared)[:20],  # Limit to top 20
                        })
                    except Exception:
                        shared_genes_along_path.append({
                            "from": path_nodes[i],
                            "to": path_nodes[i + 1],
                            "shared_count": 0,
                            "genes": [],
                        })

                result["shared_genes_along_path"] = shared_genes_along_path

            return result

        except KeyError as e:
            raise DataNotFoundError(
                detail=f"Group not found: {e}",
                resource_type="group",
            )
        except Exception as e:
            raise DBRetinaAPIError(detail=f"Path finding failed: {e}")

    @app.get("/api/v1/graph/communities")
    def get_communities(
        metric: Optional[str] = Query(None),
        cutoff: Optional[float] = Query(None),
        method: str = Query("leiden"),
    ):
        """Get community assignments for all nodes (legacy endpoint).

        For more control over clustering parameters, use POST /api/v1/graph/cluster.
        """
        graph, m, c = get_graph(metric, cutoff)
        try:
            communities = graph.community_detection(method=method)
        except ValueError:
            # community_detection raises a clean ValueError for an unknown method;
            # surface it as a 400 (mirrors /graph/cluster's algorithm validation)
            # instead of letting it fall through to the catch-all 500.
            raise ValidationError(
                detail=f"Unknown community detection method '{method}'",
                field="method",
                value=method,
                allowed_values=["leiden", "louvain"],
            )
        from collections import Counter
        sizes = dict(Counter(communities.values()).most_common())
        return {
            "communities": communities,
            "num_communities": len(sizes),
            "sizes": sizes,
        }

    # ── Advanced Clustering Endpoints ──────────────────────────────

    @app.get("/api/v1/algorithms/clustering")
    def list_clustering_algorithms():
        """List all available clustering algorithms with their parameters.

        Returns information about each algorithm including:
        - name: Algorithm identifier for API calls
        - display_name: Human-readable name
        - description: Brief explanation of the algorithm
        - parameters: List of configurable parameters with types and defaults
        """
        from .algorithms import list_algorithms
        return {"algorithms": list_algorithms()}

    class ClusteringRequest(BaseModel):
        """Request body for clustering endpoint."""
        algorithm: str = "leiden"
        parameters: dict = {}
        metric: Optional[str] = None
        cutoff: Optional[float] = None

    @app.post("/api/v1/graph/cluster")
    def run_clustering_endpoint(body: ClusteringRequest):
        """Run a clustering algorithm with custom parameters.

        Provides full control over clustering algorithm selection and parameters.
        Use GET /api/v1/algorithms/clustering to see available algorithms and
        their configurable parameters.

        Returns:
        - membership: dict mapping group names to cluster IDs
        - num_clusters: number of clusters found
        - cluster_sizes: dict mapping cluster ID to number of members
        - modularity: quality score (higher is better)
        - algorithm: name of algorithm used
        - parameters: parameters used for clustering
        """
        from .algorithms import run_clustering, get_algorithm, ALGORITHMS

        # Validate algorithm name
        if body.algorithm not in ALGORITHMS:
            valid = list(ALGORITHMS.keys())
            raise ValidationError(
                detail=f"Unknown algorithm '{body.algorithm}'",
                field="algorithm",
                value=body.algorithm,
                allowed_values=valid,
            )

        graph, m, c = get_graph(body.metric, body.cutoff)

        # Reuse the graph's already-built igraph (read-only; run_clustering does
        # not mutate it). Its edge weights ARE this metric.
        g = graph.to_igraph()
        names_map = graph._names_map
        all_ids = graph._idx_to_gid
        id_idx = graph._gid_to_idx
        weights = g.es["weight"] if g.ecount() > 0 else None

        try:
            result = run_clustering(
                g,
                algorithm=body.algorithm,
                weights=weights,
                **body.parameters,
            )
        except (TypeError, ValueError) as e:
            # ``parameters`` is a free-form, user-supplied dict forwarded into
            # igraph. A bad value (e.g. resolution="notanum") or unexpected key
            # is a client error, so report 400 instead of a 500 ALGORITHM_ERROR.
            raise ValidationError(
                detail=f"Invalid clustering parameter: {e}",
                field="parameters",
                value=body.parameters,
            )
        except Exception as e:
            raise AlgorithmError(
                detail=f"Clustering failed: {e}",
                algorithm=body.algorithm,
                reason=str(e),
            )

        # Convert membership from node indices to group names
        membership_by_name = {}
        for gid in all_ids:
            idx = id_idx[gid]
            if idx in result.membership:
                membership_by_name[names_map[gid]] = result.membership[idx]

        return {
            "membership": membership_by_name,
            "num_clusters": result.num_clusters,
            "cluster_sizes": result.cluster_sizes,
            "modularity": result.modularity,
            "algorithm": result.algorithm,
            "parameters": result.parameters,
            "metric": m,
            "cutoff": c,
        }

    @app.get("/api/v1/graph/components")
    def get_connected_components(
        metric: Optional[str] = Query(None),
        cutoff: Optional[float] = Query(None),
        min_size: int = Query(1, ge=1, description="Minimum component size"),
    ):
        """Get connected components of the graph.

        Returns groups organized by their connected component. Useful for
        identifying isolated subgraphs that have no connections to the rest
        of the network.
        """
        from .algorithms import run_clustering

        graph, m, c = get_graph(metric, cutoff)

        # Reuse the graph's already-built igraph (connected_components is read-only).
        g = graph.to_igraph()
        names_map = graph._names_map
        all_ids = graph._idx_to_gid
        id_idx = graph._gid_to_idx

        result = run_clustering(
            g,
            algorithm="connected_components",
            min_size=min_size,
        )

        # Convert membership from node indices to group names
        membership_by_name = {}
        for gid in all_ids:
            idx = id_idx[gid]
            if idx in result.membership and result.membership[idx] != -1:
                membership_by_name[names_map[gid]] = result.membership[idx]

        # Group names by component
        components: dict[int, list[str]] = {}
        for name, comp_id in membership_by_name.items():
            if comp_id not in components:
                components[comp_id] = []
            components[comp_id].append(name)

        return {
            "num_components": len(components),
            "components": components,
            "membership": membership_by_name,
            "min_size_filter": min_size,
        }

    @app.get("/api/v1/graph/layout")
    def get_layout(
        metric: Optional[str] = Query(None),
        cutoff: Optional[float] = Query(None),
        algorithm: str = Query("fr", description="Layout algorithm: fr, drl, kk, circle"),
    ):
        """Get pre-computed node positions for graph layout."""
        graph, m, c = get_graph(metric, cutoff)

        # Reuse the graph's already-built igraph (g.layout is read-only).
        g = graph.to_igraph()
        names_map = graph._names_map
        all_ids = graph._idx_to_gid
        id_idx = graph._gid_to_idx
        weights = list(g.es["weight"]) if g.ecount() > 0 else []

        algo_map = {
            "fr": "fruchterman_reingold",
            "drl": "drl",
            "kk": "kamada_kawai",
            "circle": "circle",
        }
        layout_name = algo_map.get(algorithm, "fruchterman_reingold")

        if g.vcount() > 500 and layout_name == "fruchterman_reingold":
            layout_name = "drl"

        # Pass weight values directly (not attribute name) to avoid C-level crash
        if weights and layout_name in ("fruchterman_reingold", "kamada_kawai", "drl"):
            layout = g.layout(layout_name, weights=weights)
        else:
            layout = g.layout(layout_name)

        positions = {}
        for gid in all_ids:
            idx = id_idx[gid]
            if g.degree(idx) > 0:
                positions[str(gid)] = [round(layout[idx][0], 4), round(layout[idx][1], 4)]

        return {"algorithm": layout_name, "positions": positions}

    # The raw /api/v1/cypher endpoint was REMOVED for security: untrusted Cypher
    # allowed arbitrary file read/write plus function-style readers that crash the
    # process (segfault DoS) — none of which a keyword denylist or read-only mode
    # can contain. Graph features now run entirely on in-memory igraph via the typed
    # /api/v1/graph/* endpoints (no embedded graph database).

    # ── Export Endpoints ───────────────────────────────────────────

    from fastapi.responses import StreamingResponse
    import io

    @app.get("/api/v1/export/data/{format}")
    def export_data(
        format: str,
        metric: str = Query(..., description="Metric to filter on"),
        cutoff: float = Query(0.0, ge=0),
        limit: int = Query(100000, ge=1, le=1000000),
    ):
        """Export filtered pairwise data in various formats.

        Supported formats: csv, tsv, json

        Returns the data as a downloadable file.
        """
        validate_metric(metric)
        validate_cutoff(cutoff, metric)

        if format not in ("csv", "tsv", "json"):
            raise ValidationError(
                detail=f"Unsupported export format: {format}",
                field="format",
                value=format,
                allowed_values=["csv", "tsv", "json"],
            )

        # Get data
        df = store.to_pandas(metric=metric, cutoff=cutoff, limit=limit)
        names = store.get_names_map()
        if "group_1_id" in df.columns:
            df["group_1_name"] = df["group_1_id"].map(names)
        if "group_2_id" in df.columns:
            df["group_2_name"] = df["group_2_id"].map(names)

        # Reorder columns for better readability
        cols = df.columns.tolist()
        priority = ["group_1_name", "group_2_name", "shared_features"]
        ordered = [c for c in priority if c in cols] + [c for c in cols if c not in priority]
        df = df[ordered]

        # Generate output
        output = io.BytesIO()
        filename = f"dbretina_pairs_{metric}_{int(cutoff)}"

        if format == "csv":
            df.to_csv(output, index=False)
            media_type = "text/csv"
            filename += ".csv"
        elif format == "tsv":
            df.to_csv(output, index=False, sep="\t")
            media_type = "text/tab-separated-values"
            filename += ".tsv"
        else:  # json
            output.write(df.to_json(orient="records", indent=2).encode("utf-8"))
            media_type = "application/json"
            filename += ".json"

        output.seek(0)
        return StreamingResponse(
            output,
            media_type=media_type,
            headers={"Content-Disposition": f"attachment; filename={filename}"},
        )

    @app.get("/api/v1/export/graph/{format}")
    def export_graph(
        format: str,
        metric: Optional[str] = Query(None),
        cutoff: Optional[float] = Query(None),
        include_layout: bool = Query(False, description="Include node positions"),
        include_communities: bool = Query(True, description="Include community assignments"),
    ):
        """Export graph in various network formats.

        Supported formats: graphml, gexf, json, cytoscape

        Returns the graph as a downloadable file.
        """
        if format not in ("graphml", "gexf", "json", "cytoscape"):
            raise ValidationError(
                detail=f"Unsupported graph format: {format}",
                field="format",
                value=format,
                allowed_values=["graphml", "gexf", "json", "cytoscape"],
            )

        graph, m, c = get_graph(metric, cutoff)

        # This endpoint writes vertex attributes (id/community/x/y), so it must
        # operate on a COPY of the shared graph — never mutate the cached one.
        g = graph.to_igraph().copy()
        all_ids = graph._idx_to_gid
        weights = list(g.es["weight"]) if g.ecount() > 0 else []

        # Match the old export schema: vertices carry string "id" + "name" only;
        # drop the internal int "gid" so it doesn't leak into the export.
        g.vs["id"] = [str(gid) for gid in all_ids]
        if "gid" in g.vs.attributes():
            del g.vs["gid"]

        # Add communities if requested
        if include_communities and g.ecount() > 0:
            try:
                communities = g.community_leiden(weights=weights if weights else None)
                g.vs["community"] = communities.membership
            except Exception:
                g.vs["community"] = [0] * g.vcount()

        # Add layout if requested
        if include_layout:
            try:
                layout_algo = "fr" if g.vcount() < 500 else "drl"
                layout = g.layout(layout_algo, weights=weights if weights else None)
                g.vs["x"] = [coord[0] for coord in layout]
                g.vs["y"] = [coord[1] for coord in layout]
            except Exception:
                pass

        # Generate output
        output = io.BytesIO()
        filename = f"dbretina_graph_{m}_{int(c)}"

        def _graphml_bytes(graph):
            # igraph's write_graphml needs a real file descriptor (a BytesIO
            # raises io.UnsupportedOperation: fileno), so write to a temp file.
            import os
            import tempfile

            fd, tmp = tempfile.mkstemp(suffix=".graphml")
            os.close(fd)
            try:
                graph.write_graphml(tmp)
                with open(tmp, "rb") as fh:
                    return fh.read()
            finally:
                os.unlink(tmp)

        if format == "graphml":
            output.write(_graphml_bytes(g))
            media_type = "application/xml"
            filename += ".graphml"
        elif format == "gexf":
            # igraph doesn't have native GEXF, use GraphML
            output.write(_graphml_bytes(g))
            media_type = "application/xml"
            filename += ".gexf"
        elif format in ("json", "cytoscape"):
            # Cytoscape.js JSON format
            nodes = []
            for v in g.vs:
                node_data = {"id": v["id"], "name": v["name"]}
                if "community" in v.attributes():
                    node_data["community"] = v["community"]
                if "x" in v.attributes() and "y" in v.attributes():
                    node_data["x"] = v["x"]
                    node_data["y"] = v["y"]
                nodes.append({"data": node_data})

            edges = []
            for e in g.es:
                edge_data = {
                    "id": f"e{e.index}",
                    "source": g.vs[e.source]["id"],
                    "target": g.vs[e.target]["id"],
                    "weight": e["weight"],
                    "shared_features": e["shared_features"],
                }
                edges.append({"data": edge_data})

            import json
            graph_json = {
                "elements": {"nodes": nodes, "edges": edges},
                "meta": {
                    "metric": m,
                    "cutoff": c,
                    "node_count": g.vcount(),
                    "edge_count": g.ecount(),
                },
            }
            output.write(json.dumps(graph_json, indent=2).encode("utf-8"))
            media_type = "application/json"
            filename += ".json"

        output.seek(0)
        return StreamingResponse(
            output,
            media_type=media_type,
            headers={"Content-Disposition": f"attachment; filename={filename}"},
        )

    # ── Static file serving for dashboard ───────────────────────

    dashboard_dir = pathlib.Path(__file__).parent / "dashboard_dist"
    if dashboard_dir.is_dir():
        app.mount("/", StaticFiles(directory=str(dashboard_dir), html=True), name="dashboard")

    # ── Shutdown hook ───────────────────────────────────────────

    @app.on_event("shutdown")
    def shutdown():
        app.state.graph_cache.close_all()

    return app
