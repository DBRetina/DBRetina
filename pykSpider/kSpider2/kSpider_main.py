import sys

from kSpider2.click_context import cli
from kSpider2.ks_pairwise import main as pairwise_main   # pylint: disable=relative-beyond-top-level
from kSpider2.ks_sketch_dbretina import main as dbretina_sketch
from kSpider2.ks_dataset_indexing import main as index_datasets
from kSpider2.ks_clustering import main as clustering
from kSpider2.ks_export import main as export
from kSpider2.ks_geneinfo import main as dbretina_geneinfo
from kSpider2.ks_query import main as dbretina_query
from kSpider2.ks_dedup import main as dbretina_dedup
from kSpider2.ks_bipartite_pairwise import main as bipartite_main
from kSpider2.ks_genenet import main as interactome_main
from kSpider2.ks_modularity_scan import main as modularity_main
from kSpider2.ks_setcov import main as setcov_main
from kSpider2.ks_bipartite_graph import main as grph_main
# from kSpider2.ks_clustmap import main as clustmap_main

# cli.add_command(dbretina_sketch, name="sketch")
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
# cli.add_command(clustmap_main, name="clustmap")

if __name__ == '__main__':
    cli()
