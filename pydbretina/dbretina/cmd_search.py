"""DBRetina search command — find groups by name pattern."""

import os
import sys

import click
from dbretina.click_context import cli
import dbretina.dbretina_doc_url as dbretina_doc


@cli.command(name="search", epilog=dbretina_doc.doc_url("search"), help_priority=13)
@click.option(
    "-d", "--data", "data_path", required=True,
    type=click.Path(exists=True),
    help="Path to pairwise Parquet directory or TSV file",
)
@click.argument("pattern")
@click.option(
    "-o", "--output", default=None, type=str,
    help="Output file path (default: print to stdout)",
)
@click.pass_context
def main(ctx, data_path, pattern, output):
    """Search for groups whose name contains a pattern.

    \b
    Example:
        DBRetina search -d experiment_pairwise "autism"
        DBRetina search -d experiment_pairwise "breast" -o breast_groups.tsv
    """
    LOGGER = ctx.obj

    from dbretina.compat import open_pairwise
    from dbretina.pairwise_store import PairwiseStore

    data_path = os.path.abspath(data_path)
    store = open_pairwise(data_path)
    if store is None:
        try:
            store = PairwiseStore(data_path)
        except Exception as e:
            LOGGER.ERROR(f"Could not open pairwise data at {data_path}: {e}")
            return

    try:
        matches = store.search_groups(pattern)

        if not matches:
            LOGGER.INFO(f'No groups found matching "{pattern}"')
            return

        lines = ["group_id\tgroup_name"]
        for gid, name in sorted(matches.items(), key=lambda x: x[1]):
            lines.append(f"{gid}\t{name}")

        text = "\n".join(lines)

        if output:
            with open(output, "w") as f:
                f.write(text + "\n")
            LOGGER.INFO(f'Found {len(matches)} groups matching "{pattern}" -> {output}')
        else:
            print(text)
            LOGGER.INFO(f'Found {len(matches)} groups matching "{pattern}"')
    finally:
        store.close()
