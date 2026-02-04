"""Export DBRetina pairwise data to Neo4j or Memgraph via Bolt protocol."""

import os
import sys
import click
from dbretina.click_context import cli
import dbretina.dbretina_doc_url as dbretina_doc


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
    from dbretina.pairwise_store import PairwiseStore

    store = open_pairwise(data_path)
    if store is None:
        try:
            store = PairwiseStore(data_path)
        except Exception as e:
            LOGGER.ERROR(f"Could not open pairwise data: {e}")
            return

    LOGGER.INFO(f"Loaded: {store}")
    LOGGER.INFO(f"Exporting to {uri} (database={database})")
    if metric:
        LOGGER.INFO(f"Filtering: {metric} >= {cutoff}")

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
