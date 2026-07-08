"""DBRetina serve command — starts the REST API + dashboard server."""

import os
import click
from dbretina.click_context import cli
import dbretina.dbretina_doc_url as dbretina_doc


@cli.command(name="serve", epilog=dbretina_doc.doc_url("serve"), help_priority=12)
@click.option(
    "-d", "--data", "data_path", required=True,
    type=click.Path(exists=True),
    help="Path to pairwise Parquet directory or TSV file",
)
@click.option(
    "-i", "--index-prefix", default=None, type=str,
    help="Index file prefix (for gene overlay / shared features)",
)
@click.option(
    "-m", "--metric", default="ochiai", show_default=True,
    type=str, help="Default similarity metric for graph endpoints",
)
@click.option(
    "-c", "--cutoff", default=20.0, show_default=True,
    type=float, help="Default minimum metric value for graph endpoints",
)
@click.option(
    "--port", default=8080, show_default=True,
    type=int, help="HTTP port",
)
@click.option(
    "--host", default="127.0.0.1", show_default=True,
    type=str, help="Bind address",
)
@click.option(
    "--api-key", default=None, type=str,
    help="Optional API key for authentication (header: x-api-key)",
)
@click.pass_context
def main(ctx, data_path, index_prefix, metric, cutoff, port, host, api_key):
    """Start a REST API server with interactive dashboard.

    Serves the pairwise dataset over HTTP with JSON endpoints for
    filtering, group queries, statistics, SQL, and an interactive
    network graph dashboard.

    \b
    Example:
        DBRetina serve -d experiment_pairwise --port 8080
        DBRetina serve -d experiment_pairwise -i myindex -m ochiai -c 30
        curl http://localhost:8080/api/v1/info
    """
    LOGGER = ctx.obj

    # Try importing server deps
    try:
        import uvicorn
    except ImportError:
        LOGGER.ERROR(
            "uvicorn is required for the serve command. "
            "Install with: pip install 'DBRetina[server]' "
            "or: pip install fastapi uvicorn"
        )
        return

    try:
        from fastapi import FastAPI
    except ImportError:
        LOGGER.ERROR(
            "fastapi is required for the serve command. "
            "Install with: pip install 'DBRetina[server]' "
            "or: pip install fastapi uvicorn"
        )
        return

    # Open the data
    from dbretina.compat import open_pairwise
    from dbretina.pairwise_store import PairwiseStore

    data_path = os.path.abspath(data_path)
    dbri_path = None
    if index_prefix:
        dbri_path = os.path.abspath(f"{index_prefix}.dbri")
        if not os.path.exists(dbri_path):
            LOGGER.ERROR(f"Index file not found: {dbri_path}")
            return

    store = open_pairwise(data_path)
    if store is None:
        try:
            store = PairwiseStore(data_path, dbri_path=dbri_path)
        except Exception as e:
            LOGGER.ERROR(f"Could not open pairwise data at {data_path}: {e}")
            return
    elif dbri_path:
        store._dbri_path = dbri_path

    LOGGER.INFO(f"Loaded: {store}")
    if dbri_path:
        LOGGER.INFO(f"Gene index: {dbri_path}")
    LOGGER.INFO(f"Starting server on {host}:{port}")
    LOGGER.INFO(f"Dashboard: http://{host}:{port}/")
    if api_key:
        LOGGER.INFO("API key authentication enabled")

    from dbretina.rest_api import create_app

    app = create_app(
        store,
        api_key=api_key,
        dbri_path=dbri_path,
        graph_metric=metric,
        graph_cutoff=cutoff,
    )

    try:
        uvicorn.run(app, host=host, port=port, log_level="info")
    finally:
        store.close()
