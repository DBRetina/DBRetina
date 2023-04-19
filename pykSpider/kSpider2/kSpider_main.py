import sys

from kSpider2.click_context import cli
from kSpider2.ks_pairwise import main as pairwise_main   # pylint: disable=relative-beyond-top-level
from kSpider2.ks_sketch_dbretina import main as dbretina_sketch
from kSpider2.ks_dataset_indexing import main as index_datasets
from kSpider2.ks_clustering import main as clustering
from kSpider2.ks_export import main as export
from kSpider2.ks_query_dbretina import main as dbretina_query

cli.add_command(dbretina_sketch, name="sketch")
cli.add_command(index_datasets, name="index")
cli.add_command(pairwise_main, name="pairwise")
cli.add_command(clustering, name="cluster")
cli.add_command(export, name="export")
cli.add_command(dbretina_query, name="query")



if __name__ == '__main__':
    cli()
