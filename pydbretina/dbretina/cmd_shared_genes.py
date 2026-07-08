"""DBRetina shared-genes command — show shared features between two groups."""

import os
import sys

import click
from dbretina.click_context import cli
import dbretina.dbretina_doc_url as dbretina_doc


@cli.command(name="shared-genes", epilog=dbretina_doc.doc_url("shared-genes"), help_priority=15)
@click.option(
    "-d", "--data", "data_path", required=True,
    type=click.Path(exists=True),
    help="Path to pairwise Parquet directory or TSV file",
)
@click.option(
    "-i", "--index-prefix", "index_prefix", required=True,
    type=str, help="Index file prefix (for .dbri gene data)",
)
@click.argument("group_a")
@click.argument("group_b")
@click.option(
    "-o", "--output", default=None, type=str,
    help="Output file path (default: print to stdout)",
)
@click.pass_context
def main(ctx, data_path, index_prefix, group_a, group_b, output):
    """Show shared features (genes) between two groups.

    \b
    Example:
        DBRetina shared-genes -d experiment_pairwise -i myindex \\
            "autism spectrum disorders" "schizophrenia"
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
        shared = store.shared_features(group_a, group_b)
        genes = sorted(shared)

        header = f'feature'
        lines = [header] + genes

        text = "\n".join(lines)

        if output:
            with open(output, "w") as f:
                f.write(text + "\n")
            LOGGER.INFO(
                f'{len(genes)} shared features between '
                f'"{group_a}" and "{group_b}" -> {output}'
            )
        else:
            print(text)
            LOGGER.INFO(
                f'{len(genes)} shared features between '
                f'"{group_a}" and "{group_b}"'
            )
    except KeyError as e:
        LOGGER.ERROR(str(e))
    finally:
        store.close()
