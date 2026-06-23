"""DBRetina neighbors command — show top co-occurring groups."""

import os
import sys

import click
from dbretina.click_context import cli
import dbretina.dbretina_doc_url as dbretina_doc


@cli.command(name="neighbors", epilog=dbretina_doc.doc_url("neighbors"), help_priority=14)
@click.option(
    "-d", "--data", "data_path", required=True,
    type=click.Path(exists=True),
    help="Path to pairwise Parquet directory or TSV file",
)
@click.argument("group_name")
@click.option(
    "-m", "--metric", required=True, type=str,
    help="Similarity metric to sort by (e.g. ochiai, jaccard)",
)
@click.option(
    "-c", "--cutoff", required=True, type=float,
    help="Minimum metric value (0-100)",
)
@click.option(
    "-n", "--top", default=20, show_default=True,
    type=click.IntRange(min=0),
    help="Number of top results to show",
)
@click.option(
    "-o", "--output", default=None, type=str,
    help="Output file path (default: print to stdout)",
)
@click.pass_context
def main(ctx, data_path, group_name, metric, cutoff, top, output):
    """Show top co-occurring groups for a target group.

    \b
    Example:
        DBRetina neighbors -d experiment_pairwise "autism spectrum disorders" -m ochiai -c 20
        DBRetina neighbors -d experiment_pairwise "schizophrenia" -m jaccard -c 10 -n 50
    """
    LOGGER = ctx.obj

    from dbretina.compat import open_pairwise
    from dbretina.pairwise_store import (
        PairwiseStore, LOWER_IS_BETTER, cutoff_operator, format_metric_value,
    )

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

        # p-value is "lower is better": keep pairs with pvalue <= cutoff and
        # rank the most-significant (smallest p) first; similarity metrics keep
        # >= cutoff and rank largest first (issue 068).
        order_dir = "ASC" if metric in LOWER_IS_BETTER else "DESC"
        df = store.sql(
            f"SELECT group_1_id, group_2_id, {metric}, jaccard, shared_features "
            f"FROM pairs "
            f"WHERE (group_1_id = {gid} OR group_2_id = {gid}) "
            f"AND {metric} {cutoff_operator(metric)} {cutoff} "
            f"ORDER BY {metric} {order_dir} LIMIT {top}"
        ).fetchdf()

        if len(df) == 0:
            LOGGER.INFO(
                f'No neighbors found for "{group_name}" '
                f"with {metric} {cutoff_operator(metric)} {cutoff}"
            )
            return

        df = store.resolve_names(df)

        lines = [f"neighbor\t{metric}\tjaccard\tshared_features"]
        for _, row in df.iterrows():
            if row["group_1_name"] == group_name.lower() or row["group_1_name"] == group_name:
                other = row.get("group_2_name", "")
            else:
                other = row.get("group_1_name", "")
            lines.append(
                f"{other}\t{format_metric_value(row[metric], metric)}\t"
                f"{row['jaccard']:.1f}\t{row['shared_features']}"
            )

        text = "\n".join(lines)

        if output:
            with open(output, "w") as f:
                f.write(text + "\n")
            LOGGER.INFO(
                f'Top {len(df)} neighbors of "{group_name}" '
                f"({metric} {cutoff_operator(metric)} {cutoff}) -> {output}"
            )
        else:
            print(text)
            LOGGER.INFO(
                f'Showing top {len(df)} neighbors of "{group_name}" '
                f"({metric} {cutoff_operator(metric)} {cutoff})"
            )
    except (ValueError, KeyError) as e:
        # Expected user-input errors (e.g. unknown/unavailable metric from
        # store._validate_metric) -> clean [ERROR] line, no traceback.
        # LOGGER.ERROR exits nonzero; finally below still runs.
        LOGGER.ERROR(str(e))
    finally:
        store.close()
