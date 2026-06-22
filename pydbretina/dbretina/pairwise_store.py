"""DuckDB-backed tabular query engine for DBRetina pairwise Parquet output."""

import json
import pathlib
import threading
from typing import Optional

try:
    import duckdb
except ImportError:
    duckdb = None  # type: ignore


def _require_duckdb():
    if duckdb is None:
        raise ImportError(
            "DuckDB is required for this feature. "
            "Install with: pip install 'DBRetina[server]'"
        )


import pyarrow as pa
import pyarrow.parquet as pq


class PairwiseStore:
    """Query DBRetina pairwise results stored as partitioned Parquet.

    Wraps DuckDB for efficient SQL-level access with predicate pushdown,
    column projection, and streaming reads over partitioned Parquet files.

    Usage::

        store = PairwiseStore("/path/to/experiment_pairwise")
        print(store.num_pairs, store.num_groups)

        # Streaming filtered scan
        for batch in store.filter_pairs("jaccard", 50):
            df = batch.to_pandas()

        # SQL escape hatch
        rel = store.sql("SELECT * FROM pairs WHERE ochiai > 80 LIMIT 10")
        print(rel.fetchdf())

        # Materialize to pandas
        df = store.to_pandas(metric="jaccard", cutoff=30)
    """

    METRICS = ("containment", "ochiai", "jaccard", "csi", "dice", "odds_ratio", "pvalue")

    def __init__(self, pairwise_dir: str, dbri_path: Optional[str] = None):
        _require_duckdb()
        self._dir = pathlib.Path(pairwise_dir)
        if not self._dir.is_dir():
            raise FileNotFoundError(f"Pairwise directory not found: {self._dir}")

        data_dir = self._dir / "data"
        if not data_dir.is_dir():
            raise FileNotFoundError(f"Data directory not found: {data_dir}")

        # Load manifest
        manifest_path = self._dir / "manifest.json"
        if manifest_path.exists():
            with open(manifest_path) as f:
                self._manifest = json.load(f)
        else:
            self._manifest = {}

        # Load names map
        names_path = self._dir / "names.parquet"
        if names_path.exists():
            names_table = pq.read_table(names_path)
            names_df = names_table.to_pandas()
            self._id_to_name = dict(zip(names_df["group_id"], names_df["group_name"]))
            self._name_to_id = {v: k for k, v in self._id_to_name.items()}
        else:
            self._id_to_name = {}
            self._name_to_id = {}

        # Load statistics
        stats_path = self._dir / "statistics.json"
        if stats_path.exists():
            with open(stats_path) as f:
                self._statistics = json.load(f)
        else:
            self._statistics = {}

        # Load group index for partition pruning
        group_index_path = self._dir / "group_index.parquet"
        if group_index_path.exists():
            gi_table = pq.read_table(group_index_path)
            gi_df = gi_table.to_pandas()
            self._group_partitions = dict(
                zip(gi_df["group_id"], gi_df["partition_ids"])
            )
        else:
            self._group_partitions = {}

        # DuckDB connection.
        #
        # A single PairwiseStore is shared across every ``serve`` REST endpoint,
        # and FastAPI runs the sync route handlers in a threadpool. A DuckDB
        # *connection* keeps the most recent query's result as mutable state, so
        # two threads interleaving execute()/fetch on the same connection corrupt
        # each other (fetchone() returns None -> 500; issue 055). DuckDB *cursors*
        # off one connection each carry independent result state, so per-query
        # cursors give real concurrency. We still take a short lock around cursor
        # creation (cursor() itself is not documented thread-safe) and release it
        # immediately; the returned cursor — including any streaming reader built
        # from it — is then safe to fetch from without holding the lock, which is
        # exactly what the ``sql()`` / ``to_arrow_reader()`` callers do. The
        # connection-level state set in harden_readonly() (the materialized
        # ``pairs`` table and ``enable_external_access=false``) is inherited by
        # every cursor, so the /sql sandbox is preserved.
        self._con = duckdb.connect()
        self._con_lock = threading.Lock()
        self._parquet_glob = str(data_dir / "*.parquet")

        # Lazy 'pairs' view over the parquet (ids + metrics only; cheap, no joins). The
        # served store re-materializes pairs WITH group names in harden_readonly() so the
        # /sql endpoint can query by name; CLI paths add names post-query via resolve_names().
        self._con.execute(
            f"CREATE VIEW pairs AS SELECT * FROM read_parquet('{self._parquet_glob}')"
        )

        # Gene set data (lazy-loaded from .dbri)
        self._dbri_path: Optional[str] = dbri_path
        self._gene_sets: Optional[dict[str, list[str]]] = None
        self._gene_index: Optional[dict[str, set[str]]] = None

    def close(self):
        """Close the DuckDB connection."""
        self._con.close()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()

    def _cursor(self):
        """Return a fresh DuckDB cursor with its own result state (thread-safe).

        Every concurrent query path goes through here so interleaved requests no
        longer clobber a shared connection's result (issue 055). The lock guards
        only the brief ``cursor()`` creation; callers fetch from the returned
        cursor without holding it, so streaming readers (``to_arrow_reader``)
        and the ``sql()`` escape hatch stay correct under concurrency.
        """
        with self._con_lock:
            return self._con.cursor()

    def _pairs_select_sql(self) -> str:
        """SELECT defining the served ``pairs`` table (built once in harden_readonly()).
        Joins the names map (when present) so ``pairs`` exposes ``group_1_name``/``group_2_name``
        next to the raw id/metric columns, for the /sql endpoint. The CLI lazy view stays raw."""
        names_path = self._dir / "names.parquet"
        if names_path.exists():
            return (
                f"SELECT p.*, "
                f"n1.group_name AS group_1_name, "
                f"n2.group_name AS group_2_name "
                f"FROM read_parquet('{self._parquet_glob}') p "
                f"LEFT JOIN read_parquet('{names_path}') n1 ON p.group_1_id = n1.group_id "
                f"LEFT JOIN read_parquet('{names_path}') n2 ON p.group_2_id = n2.group_id"
            )
        return f"SELECT * FROM read_parquet('{self._parquet_glob}')"

    def harden_readonly(self):
        """Sandbox the DuckDB connection for untrusted SQL (the ``serve`` /sql endpoint).

        Materializes the lazy ``pairs`` parquet view into an in-memory table, then disables
        DuckDB filesystem access, so file-reading SQL (read_text/read_csv/read_parquet/COPY/
        ATTACH/...) is rejected while queries over ``pairs`` keep working. Trade-off: loads the
        pairwise data into memory; intended for the served store only. Idempotent.
        """
        if getattr(self, "_hardened", False):
            return
        self._con.execute("DROP VIEW IF EXISTS pairs")
        self._con.execute(f"CREATE TABLE pairs AS {self._pairs_select_sql()}")
        self._con.execute("SET enable_external_access=false")
        self._hardened = True

    # ── Metadata ──────────────────────────────────────────────────────

    @property
    def directory(self) -> pathlib.Path:
        return self._dir

    @property
    def num_pairs(self) -> int:
        if "num_pairs" in self._manifest:
            return self._manifest["num_pairs"]
        # COUNT(*) on a fresh cursor always yields one row, but mirror
        # metric_summary's defensive None-handling so no fetchone() result is
        # ever unconditionally subscripted (issue 055).
        row = self._cursor().execute("SELECT COUNT(*) FROM pairs").fetchone()
        return row[0] if row is not None else 0

    @property
    def num_groups(self) -> int:
        if "num_groups" in self._manifest:
            return self._manifest["num_groups"]
        return len(self._id_to_name)

    @property
    def has_pvalue(self) -> bool:
        return self._manifest.get("has_pvalue", False)

    @property
    def available_metrics(self) -> tuple:
        """Metrics actually present in this dataset (pvalue only if it was computed)."""
        return tuple(m for m in self.METRICS if m != "pvalue" or self.has_pvalue)

    @property
    def manifest(self) -> dict:
        return dict(self._manifest)

    def get_names_map(self) -> dict[int, str]:
        return dict(self._id_to_name)

    def get_statistics(self) -> dict:
        return dict(self._statistics)

    def group_name(self, group_id: int) -> str:
        return self._id_to_name.get(group_id, f"unknown_{group_id}")

    def group_id(self, group_name: str) -> Optional[int]:
        return self._name_to_id.get(group_name.lower())

    # ── Query Methods ─────────────────────────────────────────────────

    def filter_pairs(
        self,
        metric: str,
        cutoff: float,
        columns: Optional[list[str]] = None,
    ) -> pa.RecordBatchReader:
        """Streaming filtered scan with predicate pushdown.

        Args:
            metric: Column name to filter on (e.g. "jaccard", "ochiai").
            cutoff: Minimum value for the metric (0-100 for similarity metrics).
            columns: Optional list of columns to project. None = all columns.

        Returns:
            Arrow RecordBatchReader for streaming consumption.
        """
        self._validate_metric(metric)
        col_clause = ", ".join(columns) if columns else "*"
        query = (
            f"SELECT {col_clause} FROM pairs "
            f"WHERE {metric} >= {cutoff}"
        )
        return self._cursor().execute(query).to_arrow_reader()

    def query_group(
        self,
        group_name: str,
        metric: Optional[str] = None,
        cutoff: float = 0.0,
        columns: Optional[list[str]] = None,
    ) -> pa.RecordBatchReader:
        """Query all pairs involving a specific group.

        Args:
            group_name: Name of the group to query.
            metric: Optional metric to filter on.
            cutoff: Minimum value for the metric.
            columns: Optional column projection.

        Returns:
            Arrow RecordBatchReader.
        """
        gid = self._name_to_id.get(group_name.lower())
        if gid is None:
            raise KeyError(f"Group not found: {group_name}")

        col_clause = ", ".join(columns) if columns else "*"
        where = f"WHERE (group_1_id = {gid} OR group_2_id = {gid})"
        if metric:
            self._validate_metric(metric)
            where += f" AND {metric} >= {cutoff}"

        query = f"SELECT {col_clause} FROM pairs {where}"
        return self._cursor().execute(query).to_arrow_reader()

    def shared_features(self, group_a: str, group_b: str) -> set[str]:
        """Return the set of features (genes) shared between two groups.

        Requires ``dbri_path`` to have been set at construction time.

        Args:
            group_a: Name of the first group.
            group_b: Name of the second group.

        Returns:
            Set of gene/feature names present in both groups.
        """
        gene_sets = self._load_gene_sets()
        key_a, key_b = group_a.lower(), group_b.lower()
        if key_a not in gene_sets:
            raise KeyError(f"Group not found in gene sets: {group_a}")
        if key_b not in gene_sets:
            raise KeyError(f"Group not found in gene sets: {group_b}")
        return set(gene_sets[key_a]) & set(gene_sets[key_b])

    def search_by_feature(self, feature_name: str) -> list[str]:
        """Return all group names that contain the given feature (gene).

        Requires ``dbri_path`` to have been set at construction time.

        Args:
            feature_name: Gene or feature name to search for (case-insensitive).

        Returns:
            Sorted list of group names containing this feature.
        """
        gene_index = self._get_gene_index()
        return sorted(gene_index.get(feature_name.lower(), set()))

    def search_gene_names(self, query: str, limit: int = 20) -> list[dict]:
        """Search gene names by prefix/substring match.

        Args:
            query: Search string (case-insensitive).
            limit: Max results to return.

        Returns:
            List of dicts with gene name and group count, sorted by
            exact-prefix first then substring, each sub-sorted by group count.
        """
        gene_index = self._get_gene_index()
        q = query.lower()
        prefix_matches = []
        substring_matches = []
        for gene_name, groups in gene_index.items():
            if gene_name.startswith(q):
                prefix_matches.append({"gene": gene_name, "group_count": len(groups)})
            elif q in gene_name:
                substring_matches.append({"gene": gene_name, "group_count": len(groups)})
        # Sort each bucket by group count descending (most prevalent first)
        prefix_matches.sort(key=lambda x: -x["group_count"])
        substring_matches.sort(key=lambda x: -x["group_count"])
        return (prefix_matches + substring_matches)[:limit]

    def get_group_genes(self, group_name: str) -> set[str]:
        """Return all genes/features associated with a specific group.

        Requires ``dbri_path`` to have been set at construction time.

        Args:
            group_name: Name of the group (case-insensitive).

        Returns:
            Set of gene/feature names belonging to this group.
        """
        gene_sets = self._load_gene_sets()
        key = group_name.lower()
        if key not in gene_sets:
            raise KeyError(f"Group not found in gene sets: {group_name}")
        return set(gene_sets[key])

    def iterate_all(
        self, columns: Optional[list[str]] = None
    ) -> pa.RecordBatchReader:
        """Full streaming iteration with optional column projection."""
        col_clause = ", ".join(columns) if columns else "*"
        return self._cursor().execute(f"SELECT {col_clause} FROM pairs").to_arrow_reader()

    def sql(self, query: str) -> duckdb.DuckDBPyRelation:
        """Execute arbitrary SQL against the 'pairs' view.

        The view 'pairs' contains all pairwise data. You can also use
        read_parquet() directly for advanced queries.

        Example::

            store.sql("SELECT group_1_id, AVG(jaccard) FROM pairs GROUP BY group_1_id")
        """
        # Fresh cursor: the returned relation is fetched by the caller (e.g. the
        # /sql endpoint, gene_importance), possibly after another request has run,
        # so it must own its result state rather than share the connection's.
        return self._cursor().execute(query)

    # ── Materializers ─────────────────────────────────────────────────

    def to_pandas(
        self,
        metric: Optional[str] = None,
        cutoff: float = 0.0,
        columns: Optional[list[str]] = None,
        limit: Optional[int] = None,
    ):
        """Materialize results to a pandas DataFrame.

        Args:
            metric: Optional metric to filter on.
            cutoff: Minimum value for the metric.
            columns: Optional column projection.
            limit: Optional row limit.

        Returns:
            pandas DataFrame.
        """
        import pandas as pd

        col_clause = ", ".join(columns) if columns else "*"
        where = ""
        if metric:
            self._validate_metric(metric)
            where = f"WHERE {metric} >= {cutoff}"
        limit_clause = f"LIMIT {limit}" if limit else ""
        query = f"SELECT {col_clause} FROM pairs {where} {limit_clause}"
        return self._cursor().execute(query).fetchdf()

    def to_arrow(
        self,
        metric: Optional[str] = None,
        cutoff: float = 0.0,
        columns: Optional[list[str]] = None,
    ) -> pa.Table:
        """Materialize results to an Arrow Table.

        Args:
            metric: Optional metric to filter on.
            cutoff: Minimum value for the metric.
            columns: Optional column projection.

        Returns:
            Arrow Table.
        """
        col_clause = ", ".join(columns) if columns else "*"
        where = ""
        if metric:
            self._validate_metric(metric)
            where = f"WHERE {metric} >= {cutoff}"
        query = f"SELECT {col_clause} FROM pairs {where}"
        return self._cursor().execute(query).fetch_arrow_table()

    # ── Name Resolution Helpers ───────────────────────────────────────

    def resolve_names(self, df):
        """Add group_1_name and group_2_name columns to a DataFrame.

        Modifies the DataFrame in-place and returns it.
        """
        if "group_1_id" in df.columns:
            df["group_1_name"] = df["group_1_id"].map(self._id_to_name)
        if "group_2_id" in df.columns:
            df["group_2_name"] = df["group_2_id"].map(self._id_to_name)
        return df

    def search_groups(self, pattern: str) -> dict[int, str]:
        """Search for groups whose name contains the pattern (case-insensitive).

        Args:
            pattern: Substring to search for in group names.

        Returns:
            Dict mapping group_id -> group_name for all matches.
        """
        needle = pattern.lower()
        return {
            gid: name
            for gid, name in self._id_to_name.items()
            if needle in name.lower()
        }

    # ── Aggregation Helpers ───────────────────────────────────────────

    def top_pairs(self, metric: str, n: int = 20):
        """Return top N pairs by a given metric as a named DataFrame."""
        self._validate_metric(metric)
        df = self._cursor().execute(
            f"SELECT * FROM pairs ORDER BY {metric} DESC LIMIT {n}"
        ).fetchdf()
        return self.resolve_names(df)

    def metric_summary(self, metric: str) -> dict:
        """Return summary statistics for a metric."""
        self._validate_metric(metric)
        row = self._cursor().execute(
            f"SELECT COUNT(*) as count, "
            f"MIN({metric}) as min, "
            f"MAX({metric}) as max, "
            f"AVG({metric}) as mean, "
            f"MEDIAN({metric}) as median, "
            f"STDDEV({metric}) as stddev "
            f"FROM pairs"
        ).fetchone()
        # Defensive: an aggregate over a (possibly empty) table should always yield
        # one row, but never unconditionally subscript a None fetchone() result —
        # that was the 500 surface for the connection race (issue 055).
        if row is None:
            return {"count": 0, "min": None, "max": None,
                    "mean": None, "median": None, "stddev": None}
        return {
            "count": row[0],
            "min": row[1],
            "max": row[2],
            "mean": row[3],
            "median": row[4],
            "stddev": row[5],
        }

    def group_metric_profile(self, group_name: str) -> list[dict]:
        """Compute AVG, MAX, COUNT for each metric column for a given group.

        Aggregates across all pairs where the group appears as either
        group_1 or group_2.

        Args:
            group_name: Name of the group (case-insensitive).

        Returns:
            List of dicts: [{"metric": str, "avg": float, "max": float, "count": int}, ...]

        Raises:
            KeyError: If the group is not found in the names map.
        """
        gid = self._name_to_id.get(group_name.lower())
        if gid is None:
            raise KeyError(f"Group not found: {group_name}")

        parts = []
        for m in self.available_metrics:
            parts.append(
                f"SELECT '{m}' AS metric, "
                f"AVG({m}) AS avg, "
                f"MAX({m}) AS max, "
                f"COUNT(*) AS count "
                f"FROM pairs "
                f"WHERE group_1_id = {gid} OR group_2_id = {gid}"
            )
        query = " UNION ALL ".join(parts)
        rows = self._cursor().execute(query).fetchall()
        return [
            {"metric": r[0], "avg": r[1], "max": r[2], "count": r[3]}
            for r in rows
        ]

    def group_pair_counts(self) -> dict[int, int]:
        """Count how many pairs each group participates in."""
        rows = self._cursor().execute(
            "SELECT gid, COUNT(*) as cnt FROM ("
            "  SELECT group_1_id as gid FROM pairs"
            "  UNION ALL"
            "  SELECT group_2_id as gid FROM pairs"
            ") GROUP BY gid ORDER BY cnt DESC"
        ).fetchall()
        return {r[0]: r[1] for r in rows}

    # ── Internal ──────────────────────────────────────────────────────

    def _load_gene_sets(self) -> dict[str, list[str]]:
        """Lazy-load raw gene sets from the .dbri index file."""
        if self._gene_sets is not None:
            return self._gene_sets
        if self._dbri_path is None:
            raise ValueError(
                "No dbri_path provided. Pass dbri_path= to PairwiseStore() "
                "to enable gene-level queries."
            )
        import _dbretina_internal as dbi

        raw = dbi.dbri_load_raw_gene_sets(self._dbri_path)
        if not raw:
            raise RuntimeError(
                f".dbri file lacks RAW_GENE_SETS section: {self._dbri_path}"
            )
        self._gene_sets = json.loads(raw).get("data", {})
        return self._gene_sets

    def _get_gene_index(self) -> dict[str, set[str]]:
        """Build or return cached inverted index: gene_name -> {group_names}."""
        if self._gene_index is not None:
            return self._gene_index
        gene_sets = self._load_gene_sets()
        index: dict[str, set[str]] = {}
        for group_name, genes in gene_sets.items():
            for gene in genes:
                index.setdefault(gene, set()).add(group_name)
        self._gene_index = index
        return self._gene_index

    def _validate_metric(self, metric: str):
        if metric not in self.available_metrics:
            raise ValueError(
                f"Unknown metric '{metric}'. Valid: {', '.join(self.available_metrics)}"
            )

    def __repr__(self):
        return (
            f"PairwiseStore('{self._dir}', "
            f"pairs={self.num_pairs:,}, groups={self.num_groups})"
        )

    def _repr_html_(self):
        """Rich HTML display for Jupyter notebooks."""
        manifest = self.manifest
        stats = self.get_statistics()
        metric_rows = ""
        for m in ("containment", "ochiai", "jaccard", "csi", "dice"):
            if m in stats.get("metrics", {}):
                s = stats["metrics"][m]
                metric_rows += (
                    f"<tr><td>{m}</td>"
                    f"<td>{s.get('min', ''):.1f}</td>"
                    f"<td>{s.get('mean', ''):.1f}</td>"
                    f"<td>{s.get('median', ''):.1f}</td>"
                    f"<td>{s.get('max', ''):.1f}</td></tr>"
                )
        if not metric_rows:
            # Fall back to live summary if no pre-computed stats
            for m in ("containment", "ochiai", "jaccard", "csi", "dice"):
                try:
                    s = self.metric_summary(m)
                    metric_rows += (
                        f"<tr><td>{m}</td>"
                        f"<td>{s['min']:.1f}</td>"
                        f"<td>{s['mean']:.1f}</td>"
                        f"<td>{s['median']:.1f}</td>"
                        f"<td>{s['max']:.1f}</td></tr>"
                    )
                except Exception:
                    pass

        cutoff_info = ""
        if manifest.get("cutoff_metric"):
            cutoff_info = (
                f"<tr><td>Cutoff</td>"
                f"<td colspan='4'>{manifest['cutoff_metric']} "
                f"&ge; {manifest.get('cutoff_threshold', 0)}</td></tr>"
            )

        return (
            f"<div style='border:1px solid #ddd;padding:12px;border-radius:6px;"
            f"max-width:500px;font-family:sans-serif;font-size:13px'>"
            f"<b>PairwiseStore</b><br>"
            f"<table style='margin:8px 0;border-collapse:collapse'>"
            f"<tr><td>Pairs</td><td colspan='4'><b>{self.num_pairs:,}</b></td></tr>"
            f"<tr><td>Groups</td><td colspan='4'><b>{self.num_groups:,}</b></td></tr>"
            f"<tr><td>Population</td><td colspan='4'>{manifest.get('population_size', '?'):,}</td></tr>"
            f"{cutoff_info}"
            f"<tr><td colspan='5' style='padding-top:6px'><b>Metric Summary</b></td></tr>"
            f"<tr style='font-size:11px;color:#666'>"
            f"<td></td><td>min</td><td>mean</td><td>median</td><td>max</td></tr>"
            f"{metric_rows}"
            f"</table>"
            f"<span style='color:#888;font-size:11px'>{self._dir}</span>"
            f"</div>"
        )
