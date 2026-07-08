"""DBRetina connect command — shortest path between two groups in the similarity graph."""

import os

import click
from dbretina.click_context import cli
import dbretina.dbretina_doc_url as dbretina_doc


@cli.command(name="connect", epilog=dbretina_doc.doc_url("connect"), help_priority=18)
@click.option(
    "-d", "--data", "data_path", required=True,
    type=click.Path(exists=True),
    help="Path to pairwise Parquet directory or TSV file",
)
@click.argument("group_a")
@click.argument("group_b")
@click.option(
    "-m", "--metric", required=True, type=str,
    help="Similarity metric used for the edges",
)
@click.option(
    "-c", "--cutoff", required=True, type=float,
    help="Minimum metric value for an edge (0-100)",
)
@click.option(
    "-o", "--output", default=None, type=str,
    help="Output file path (default: print to stdout)",
)
@click.pass_context
def main(ctx, data_path, group_a, group_b, metric, cutoff, output):
    """Find the shortest path between two groups through the similarity network.

    Reports each group along the path and the edge metric at each step, which shows
    how two groups that do not overlap directly are linked through intermediaries.

    \b
    Example:
        DBRetina connect -d experiment_pairwise "rheumatoid arthritis" "atorvastatin" -m containment -c 40
    """
    LOGGER = ctx.obj

    from dbretina.compat import open_pairwise
    from dbretina.pairwise_store import PairwiseStore, cutoff_operator, format_metric_value
    from dbretina.pairwise_graph import PairwiseGraph

    data_path = os.path.abspath(data_path)
    store = open_pairwise(data_path)
    if store is None:
        try:
            store = PairwiseStore(data_path)
        except Exception as e:
            LOGGER.ERROR(f"Could not open pairwise data at {data_path}: {e}")
            return

    try:
        for name in (group_a, group_b):
            if store.group_id(name) is None:
                LOGGER.ERROR(f'Group not found: "{name}"')
                return
        store._validate_metric(metric)

        graph = PairwiseGraph(store, metric=metric, cutoff=cutoff)
        result = graph.shortest_path_full(group_a, group_b)

        if not result["connected"]:
            LOGGER.INFO(
                f'"{group_a}" and "{group_b}" are not connected at '
                f"{metric} {cutoff_operator(metric)} {cutoff}"
            )
            return

        path = result["path_nodes"]
        ig = graph.to_igraph()
        lines = [f"step\tgroup\t{metric}"]
        lines.append(f"0\t{path[0]}\t-")
        for i in range(1, len(path)):
            va = ig.vs.find(name=path[i - 1]).index
            vb = ig.vs.find(name=path[i]).index
            weight = ig.es[ig.get_eid(va, vb)]["weight"]
            lines.append(f"{i}\t{path[i]}\t{format_metric_value(weight, metric)}")

        text = "\n".join(lines)
        if output:
            with open(output, "w") as f:
                f.write(text + "\n")
            LOGGER.INFO(
                f'{result["path_length"]}-hop path from "{group_a}" to "{group_b}" -> {output}'
            )
        else:
            print(text)
            LOGGER.INFO(
                f'{result["path_length"]}-hop path from "{group_a}" to "{group_b}"'
            )
    except (ValueError, KeyError) as e:
        LOGGER.ERROR(str(e))
    finally:
        store.close()
