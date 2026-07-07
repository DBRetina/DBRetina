"""DBRetina module command — a group's characteristic module (its top neighbors)."""

import os

import click
from dbretina.click_context import cli
import dbretina.dbretina_doc_url as dbretina_doc


@cli.command(name="module", epilog=dbretina_doc.doc_url("module"), help_priority=19)
@click.option(
    "-d", "--data", "data_path", required=True,
    type=click.Path(exists=True),
    help="Path to pairwise Parquet directory or TSV file",
)
@click.argument("group_name")
@click.option(
    "-m", "--metric", required=True, type=str,
    help="Similarity metric to rank neighbors by",
)
@click.option(
    "-c", "--cutoff", required=True, type=float,
    help="Minimum metric value (0-100)",
)
@click.option(
    "--min-shared", "min_shared", default=0, show_default=True,
    type=click.IntRange(min=0),
    help="Drop neighbors sharing fewer than this many features",
)
@click.option(
    "-n", "--top", default=0, show_default=True,
    type=click.IntRange(min=0),
    help="Keep only the top N neighbors (0 = all)",
)
@click.option(
    "-o", "--output", default=None, type=str,
    help="Output file path (default: print to stdout)",
)
@click.pass_context
def main(ctx, data_path, group_name, metric, cutoff, min_shared, top, output):
    """Extract a group's characteristic module: the groups most similar to it.

    Prints one group name per line (the seed's neighbors passing the cutoff and
    --min-shared), a list you can feed to ``query -g``, ``geneinfo -g``, or
    ``enrich --module``.

    \b
    Example:
        DBRetina module -d experiment_pairwise "rheumatoid arthritis" -m containment -c 50 --min-shared 10 -o ra_module.txt
    """
    LOGGER = ctx.obj

    from dbretina.compat import open_pairwise
    from dbretina.pairwise_store import PairwiseStore, LOWER_IS_BETTER, cutoff_operator

    data_path = os.path.abspath(data_path)
    store = open_pairwise(data_path)
    if store is None:
        try:
            store = PairwiseStore(data_path)
        except Exception as e:
            LOGGER.ERROR(f"Could not open pairwise data at {data_path}: {e}")
            return

    try:
        gid = store.group_id(group_name)
        if gid is None:
            LOGGER.ERROR(f'Group not found: "{group_name}"')
            return
        store._validate_metric(metric)

        order_dir = "ASC" if metric in LOWER_IS_BETTER else "DESC"
        limit = f"LIMIT {top}" if top else ""
        df = store.sql(
            f"SELECT group_1_id, group_2_id, {metric} FROM pairs "
            f"WHERE (group_1_id = {gid} OR group_2_id = {gid}) "
            f"AND {metric} {cutoff_operator(metric)} {cutoff} "
            f"AND shared_features >= {min_shared} "
            f"ORDER BY {metric} {order_dir} {limit}"
        ).fetchdf()

        if len(df) == 0:
            LOGGER.INFO(f'No module members found for "{group_name}"')
            return

        df = store.resolve_names(df)
        seed = group_name.lower()
        members = [
            row["group_2_name"] if row["group_1_name"] == seed else row["group_1_name"]
            for _, row in df.iterrows()
        ]

        text = "\n".join(members)
        if output:
            with open(output, "w") as f:
                f.write(text + "\n")
            LOGGER.INFO(f'Module of "{group_name}": {len(members)} groups -> {output}')
        else:
            print(text)
            LOGGER.INFO(f'Module of "{group_name}": {len(members)} groups')
    except (ValueError, KeyError) as e:
        LOGGER.ERROR(str(e))
    finally:
        store.close()
