import sys

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

if __name__ == '__main__':
    cli()
