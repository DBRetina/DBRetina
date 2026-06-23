"""Export DBRetina pairwise data to Neo4j or Memgraph via Bolt protocol."""

import os
import sys
import click
from dbretina.click_context import cli
import dbretina.dbretina_doc_url as dbretina_doc


class _TsvPairwiseStore:
    """Minimal PairwiseStore-shaped adapter over a standalone pairwise TSV.

    ``-d`` advertises a "Pairwise Parquet directory or TSV file", and the sibling
    ``export`` command ingests a bare pairwise TSV via a pandas read_csv fallback.
    export-neo4j, however, only ever opened a parquet store, so a genuine standalone
    TSV (no sibling parquet dir / .dbrp) was rejected with "Pairwise directory not
    found" (issue 070). This adapter exposes exactly the surface ``export_to_neo4j``
    consumes -- ``get_names_map()``, ``has_pvalue``, ``filter_pairs(metric, cutoff)``,
    ``iterate_all()``, ``close()`` -- reading the TSV with pandas (same dependency as
    ``export``) and yielding pyarrow RecordBatches so the existing
    ``record_batch.to_pandas()`` streaming loop is unchanged.
    """

    # TSV header -> the lowercase column names export_to_neo4j reads.
    _RENAME = {"group_1_ID": "group_1_id", "group_2_ID": "group_2_id"}

    def __init__(self, tsv_path: str):
        import pandas as pd

        self._path = tsv_path
        df = pd.read_csv(tsv_path, sep="\t", comment="#")
        df = df.rename(columns=self._RENAME)
        required = {"group_1_id", "group_2_id", "shared_features"}
        missing = required - set(df.columns)
        if missing:
            raise ValueError(
                f"'{tsv_path}' is not a DBRetina pairwise TSV "
                f"(missing columns: {', '.join(sorted(missing))})"
            )
        self._df = df
        self._has_pvalue = "pvalue" in df.columns

        # id -> name map from the two id/name column pairs (ids are unique per group).
        names: dict[int, str] = {}
        if "group_1_name" in df.columns:
            for gid, name in zip(df["group_1_id"], df["group_1_name"]):
                names[int(gid)] = name
        if "group_2_name" in df.columns:
            for gid, name in zip(df["group_2_id"], df["group_2_name"]):
                names[int(gid)] = name
        self._names = names

    @property
    def has_pvalue(self) -> bool:
        return self._has_pvalue

    def get_names_map(self) -> dict:
        return dict(self._names)

    def _batches(self, df):
        """Yield the frame as pyarrow RecordBatches (matches the parquet reader API)."""
        import pyarrow as pa

        return pa.Table.from_pandas(df, preserve_index=False).to_batches()

    def iterate_all(self):
        return self._batches(self._df)

    def filter_pairs(self, metric: str, cutoff: float):
        from dbretina.pairwise_store import passes_cutoff

        if metric not in self._df.columns:
            raise ValueError(
                f"Unknown metric '{metric}'. Valid: "
                f"{', '.join(c for c in self._df.columns if c not in ('group_1_id', 'group_2_id', 'group_1_name', 'group_2_name', 'shared_features'))}"
            )
        # Metric-aware filter: p-value keeps value <= cutoff, similarity metrics
        # keep value >= cutoff (issue 068 semantics), via the shared helper.
        mask = self._df[metric].apply(lambda v: passes_cutoff(float(v), cutoff, metric))
        return self._batches(self._df[mask])

    def close(self):
        pass

    def __repr__(self):
        return (
            f"TsvPairwiseStore('{self._path}', "
            f"pairs={len(self._df):,}, groups={len(self._names)})"
        )


def export_to_neo4j(
    store,
    uri: str,
    user: str = "neo4j",
    password: str = "neo4j",
    database: str = "neo4j",
    metric: str = None,
    cutoff: float = 0.0,
    batch_size: int = 5000,
    clear: bool = False,
    on_progress=None,
):
    """Stream pairwise data from a PairwiseStore into Neo4j/Memgraph.

    Creates ``(:Group {id, name})`` nodes and
    ``[:SIMILAR_TO {shared_features, containment, ochiai, jaccard, ...}]``
    relationships via batched UNWIND Cypher inserts.

    Args:
        store: Open PairwiseStore instance.
        uri: Bolt URI (e.g. ``bolt://localhost:7687``).
        user: Database username.
        password: Database password.
        database: Database name (Neo4j only, ignored by Memgraph).
        metric: Optional metric to filter edges.
        cutoff: Minimum metric value for edge inclusion.
        batch_size: Records per UNWIND transaction.
        clear: If True, delete all nodes/relationships first.
        on_progress: Optional callback ``(nodes_done, edges_done) -> None``.

    Returns:
        Tuple of (num_nodes, num_edges) inserted.
    """
    try:
        from neo4j import GraphDatabase
    except ImportError:
        raise ImportError(
            "neo4j driver required. Install with: pip install neo4j"
        )

    driver = GraphDatabase.driver(uri, auth=(user, password))

    try:
        with driver.session(database=database) as session:
            # Optionally clear existing data
            if clear:
                session.run("MATCH (n) DETACH DELETE n")

            # Create constraints/indexes for performance
            try:
                session.run(
                    "CREATE CONSTRAINT group_id IF NOT EXISTS "
                    "FOR (g:Group) REQUIRE g.id IS UNIQUE"
                )
            except Exception:
                # Memgraph or older Neo4j may not support this syntax
                try:
                    session.run("CREATE INDEX ON :Group(id)")
                except Exception:
                    pass

            # Insert nodes in batches
            names_map = store.get_names_map()
            node_batches = []
            batch = []
            for gid, name in names_map.items():
                batch.append({"id": int(gid), "name": name})
                if len(batch) >= batch_size:
                    node_batches.append(batch)
                    batch = []
            if batch:
                node_batches.append(batch)

            nodes_done = 0
            for b in node_batches:
                session.run(
                    "UNWIND $rows AS row "
                    "MERGE (g:Group {id: row.id}) "
                    "SET g.name = row.name",
                    rows=b,
                )
                nodes_done += len(b)
                if on_progress:
                    on_progress(nodes_done, 0)

            # Insert edges in batches via streaming
            edge_cols = [
                "group_1_id", "group_2_id", "shared_features",
                "containment", "ochiai", "jaccard", "csi", "dice", "odds_ratio",
            ]
            if store.has_pvalue:
                edge_cols.append("pvalue")

            reader = store.filter_pairs(metric, cutoff) if metric else store.iterate_all()

            edges_done = 0
            batch = []
            for record_batch in reader:
                df = record_batch.to_pandas()
                for _, row in df.iterrows():
                    edge = {
                        "src": int(row["group_1_id"]),
                        "dst": int(row["group_2_id"]),
                        "shared_features": int(row["shared_features"]),
                        "containment": float(row["containment"]),
                        "ochiai": float(row["ochiai"]),
                        "jaccard": float(row["jaccard"]),
                        "csi": float(row["csi"]),
                        "dice": float(row["dice"]),
                        "odds_ratio": float(row["odds_ratio"]),
                    }
                    if store.has_pvalue and "pvalue" in row:
                        edge["pvalue"] = float(row["pvalue"])
                    batch.append(edge)

                    if len(batch) >= batch_size:
                        _insert_edges(session, batch, store.has_pvalue)
                        edges_done += len(batch)
                        batch = []
                        if on_progress:
                            on_progress(nodes_done, edges_done)

            if batch:
                _insert_edges(session, batch, store.has_pvalue)
                edges_done += len(batch)
                if on_progress:
                    on_progress(nodes_done, edges_done)

    finally:
        driver.close()

    return nodes_done, edges_done


def _insert_edges(session, batch, has_pvalue):
    """UNWIND a batch of edge dicts into SIMILAR_TO relationships."""
    props = (
        "shared_features: row.shared_features, "
        "containment: row.containment, "
        "ochiai: row.ochiai, "
        "jaccard: row.jaccard, "
        "csi: row.csi, "
        "dice: row.dice, "
        "odds_ratio: row.odds_ratio"
    )
    if has_pvalue:
        props += ", pvalue: row.pvalue"

    session.run(
        "UNWIND $rows AS row "
        "MATCH (a:Group {id: row.src}) "
        "MATCH (b:Group {id: row.dst}) "
        f"CREATE (a)-[:SIMILAR_TO {{{props}}}]->(b)",
        rows=batch,
    )


# ── CLI Command ────────────────────────────────────────────────────


@cli.command(name="export-neo4j", epilog=dbretina_doc.doc_url("export-neo4j"), help_priority=13)
@click.option("-d", "--data", "data_path", required=True, type=click.Path(exists=True),
              help="Pairwise Parquet directory or TSV file")
@click.option("--uri", required=True, type=str, help="Bolt URI (e.g. bolt://localhost:7687)")
@click.option("--user", default="neo4j", show_default=True, help="Database username")
@click.option("--password", required=True, type=str, help="Database password")
@click.option("--database", default="neo4j", show_default=True, help="Database name (Neo4j only)")
@click.option("-m", "--metric", default=None, type=str, help="Filter edges by metric")
@click.option("-c", "--cutoff", default=0.0, type=float, help="Minimum metric value")
@click.option("--batch-size", default=5000, show_default=True, type=int, help="Records per transaction")
@click.option("--clear", is_flag=True, default=False, help="Delete all existing data first")
@click.pass_context
def main(ctx, data_path, uri, user, password, database, metric, cutoff, batch_size, clear):
    """Export pairwise data to Neo4j or Memgraph.

    \b
    Streams groups as (:Group) nodes and similarities as
    [:SIMILAR_TO] relationships via batched Cypher UNWIND inserts.

    \b
    Example:
        DBRetina export-neo4j -d experiment_pairwise \\
            --uri bolt://localhost:7687 --password secret \\
            -m jaccard -c 50
    """
    LOGGER = ctx.obj

    try:
        import neo4j  # noqa: F401
    except ImportError:
        LOGGER.ERROR(
            "neo4j driver required. Install with: pip install neo4j"
        )
        return

    data_path = os.path.abspath(data_path)

    from dbretina.compat import open_pairwise
    from dbretina.pairwise_store import PairwiseStore, cutoff_operator

    # Resolve -d to a store, mirroring how the sibling `export` command accepts
    # the same inputs (issue 070):
    #   1. open_pairwise: a parquet directory, or a .tsv WITH a sibling parquet dir.
    #   2. PairwiseStore(path): a parquet directory passed directly.
    #   3. a standalone pairwise .tsv (no parquet sibling) -> the TSV adapter, just
    #      like export's read_csv fallback. Previously this third case was rejected
    #      with "Pairwise directory not found" despite -d's help promising a TSV file.
    store = open_pairwise(data_path)
    if store is None and os.path.isdir(data_path):
        try:
            store = PairwiseStore(data_path)
        except Exception as e:
            LOGGER.ERROR(f"Could not open pairwise data: {e}")
            return
    if store is None:
        try:
            store = _TsvPairwiseStore(data_path)
        except Exception as e:
            LOGGER.ERROR(f"Could not open pairwise data: {e}")
            return

    LOGGER.INFO(f"Loaded: {store}")
    LOGGER.INFO(f"Exporting to {uri} (database={database})")
    if metric:
        # Metric-aware operator in the log: p-value filters keep value <= cutoff
        # (lower is better), similarity metrics keep value >= cutoff. The actual
        # PairwiseStore/adapter filter already uses this; the log hardcoded '>='.
        LOGGER.INFO(f"Filtering: {metric} {cutoff_operator(metric)} {cutoff}")

    def progress(nodes, edges):
        print(f"\r  nodes={nodes:,}  edges={edges:,}", end="", flush=True)

    try:
        n_nodes, n_edges = export_to_neo4j(
            store, uri=uri, user=user, password=password,
            database=database, metric=metric, cutoff=cutoff,
            batch_size=batch_size, clear=clear, on_progress=progress,
        )
        print()
        LOGGER.SUCCESS(f"Exported {n_nodes:,} nodes and {n_edges:,} edges to {uri}")
    except Exception as e:
        print()
        LOGGER.ERROR(f"Export failed: {e}")
    finally:
        store.close()
