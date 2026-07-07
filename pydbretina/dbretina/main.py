import sys

import click
from dbretina.click_context import cli
from dbretina.pairwise import main as pairwise_main
from dbretina.dataset_indexing import main as index_datasets
from dbretina.clustering import main as clustering
from dbretina.export import main as export
from dbretina.geneinfo import main as dbretina_geneinfo
from dbretina.query import main as dbretina_query
from dbretina.dedup import main as dbretina_dedup
from dbretina.bipartite_pairwise import main as bipartite_main
from dbretina.genenet import main as interactome_main
from dbretina.modularity_scan import main as modularity_main
from dbretina.setcov import main as setcov_main
from dbretina.bipartite_graph import main as grph_main
from dbretina.append import main as append_main
from dbretina.merge import main as merge_main


def _safe_add_command(module_path: str, command_attr: str, cli_name: str):
    """Safely import and register a CLI command with helpful error messages."""
    try:
        module = __import__(f"dbretina.{module_path}", fromlist=[command_attr])
        cli.add_command(getattr(module, command_attr), name=cli_name)
    except ImportError as e:
        @cli.command(name=cli_name)
        def _missing_cmd(err=str(e)):
            f"""Command '{cli_name}' requires additional dependencies."""
            raise click.ClickException(
                f"Command '{cli_name}' is unavailable: {err}\n"
                f"Try: pip install 'DBRetina[server]' or pip install 'DBRetina[all]'"
            )

# Core commands (always available)
cli.add_command(index_datasets, name="index")
cli.add_command(pairwise_main, name="pairwise")
cli.add_command(dbretina_query, name="query")
cli.add_command(clustering, name="cluster")
cli.add_command(export, name="export")
cli.add_command(dbretina_geneinfo, name="geneinfo")
cli.add_command(dbretina_dedup, name="dedup")
cli.add_command(bipartite_main, name="bipartite")
cli.add_command(interactome_main, name="interactome")
cli.add_command(modularity_main, name="modularity")
cli.add_command(setcov_main, name="setcov")
cli.add_command(grph_main, name="graph")
cli.add_command(append_main, name="append")
cli.add_command(merge_main, name="merge")

# Optional commands with external dependencies
_safe_add_command("serve", "main", "serve")
_safe_add_command("export_neo4j", "main", "export-neo4j")
_safe_add_command("cmd_search", "main", "search")
_safe_add_command("cmd_neighbors", "main", "neighbors")
_safe_add_command("cmd_shared_genes", "main", "shared-genes")
_safe_add_command("cmd_gene_search", "main", "gene-search")
_safe_add_command("cmd_genescore", "main", "genescore")
_safe_add_command("cmd_connect", "main", "connect")
_safe_add_command("cmd_module", "main", "module")
_safe_add_command("cmd_enrich", "main", "enrich")


if __name__ == '__main__':
    cli()
