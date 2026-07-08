"""DBRetina gene-search command — find groups containing a gene/feature."""

import os
import sys

import click
from dbretina.click_context import cli
import dbretina.dbretina_doc_url as dbretina_doc


@cli.command(name="gene-search", epilog=dbretina_doc.doc_url("gene-search"), help_priority=16)
@click.option(
    "-d", "--data", "data_path", required=True,
    type=click.Path(exists=True),
    help="Path to pairwise Parquet directory or TSV file",
)
@click.option(
    "-i", "--index-prefix", "index_prefix", required=True,
    type=str, help="Index file prefix (for .dbri gene data)",
)
@click.argument("feature")
@click.option(
    "-o", "--output", default=None, type=str,
    help="Output file path (default: print to stdout)",
)
@click.pass_context
def main(ctx, data_path, index_prefix, feature, output):
    """Find all groups containing a given gene/feature.

    The search is case-insensitive.

    \b
    Example:
        DBRetina gene-search -d experiment_pairwise -i myindex "SHANK3"
        DBRetina gene-search -d experiment_pairwise -i myindex "PTEN" -o pten_diseases.tsv
    """
    LOGGER = ctx.obj

    from dbretina.compat import open_pairwise
    from dbretina.pairwise_store import PairwiseStore

    data_path = os.path.abspath(data_path)
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
    else:
        store._dbri_path = dbri_path

    try:
        groups = store.search_by_feature(feature)

        if not groups:
            LOGGER.INFO(f'{feature.upper()} not found in any group')
            return

        lines = ["group_name"] + groups
        text = "\n".join(lines)

        if output:
            with open(output, "w") as f:
                f.write(text + "\n")
            LOGGER.INFO(f'{feature.upper()} found in {len(groups)} groups -> {output}')
        else:
            print(text)
            LOGGER.INFO(f'{feature.upper()} found in {len(groups)} groups')
    finally:
        store.close()
