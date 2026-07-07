"""DBRetina enrich command — rank groups by enrichment for a reference module."""

import os
from collections import defaultdict

import click
from dbretina.click_context import cli
import dbretina.dbretina_doc_url as dbretina_doc


@cli.command(name="enrich", epilog=dbretina_doc.doc_url("enrich"), help_priority=20)
@click.option(
    "-d", "--data", "data_path", required=True,
    type=click.Path(exists=True),
    help="Path to pairwise Parquet directory or TSV file",
)
@click.option(
    "--module", "module_file", required=True,
    type=click.Path(exists=True),
    help="Single-column file of reference group names (the module)",
)
@click.option(
    "-m", "--metric", required=True, type=str,
    help="Similarity metric used for the neighbor edges",
)
@click.option(
    "-c", "--cutoff", required=True, type=float,
    help="Minimum metric value for a neighbor edge (0-100)",
)
@click.option(
    "--min-shared", "min_shared", default=0, show_default=True,
    type=click.IntRange(min=0),
    help="Drop neighbor edges sharing fewer than this many features",
)
@click.option(
    "--exclude-module", "exclude_module", is_flag=True, default=False,
    help="Drop the module's own members from the ranking (keep only new candidates)",
)
@click.option(
    "-n", "--top", default=25, show_default=True,
    type=click.IntRange(min=0),
    help="Number of top groups to show (0 = all)",
)
@click.option(
    "-o", "--output", default=None, type=str,
    help="Output file path (default: print to stdout)",
)
@click.pass_context
def main(ctx, data_path, module_file, metric, cutoff, min_shared, exclude_module, top, output):
    """Rank every group by enrichment of a reference module among its neighbors.

    For each group, counts how many of its neighbors (edges passing the cutoff and
    --min-shared) belong to the reference module, and scores that with a
    hypergeometric test over the group universe. Groups whose neighborhood is made
    of module members rank first. The reported p-value is uncorrected. Pass
    --exclude-module to drop the module's own members and keep only new candidates.

    \b
    Example:
        DBRetina enrich -d experiment_pairwise --module ra_module.txt -m containment -c 40
    """
    LOGGER = ctx.obj

    from scipy.stats import hypergeom

    from dbretina.compat import open_pairwise
    from dbretina.pairwise_store import PairwiseStore, cutoff_operator

    data_path = os.path.abspath(data_path)
    store = open_pairwise(data_path)
    if store is None:
        try:
            store = PairwiseStore(data_path)
        except Exception as e:
            LOGGER.ERROR(f"Could not open pairwise data at {data_path}: {e}")
            return

    try:
        store._validate_metric(metric)

        with open(module_file) as f:
            module_names = {line.strip().lower() for line in f if line.strip()}
        module_ids = {store.group_id(n) for n in module_names}
        module_ids.discard(None)
        if not module_ids:
            LOGGER.ERROR("None of the module groups were found in the data")
            return

        df = store.sql(
            f"SELECT group_1_id, group_2_id FROM pairs "
            f"WHERE {metric} {cutoff_operator(metric)} {cutoff} "
            f"AND shared_features >= {min_shared}"
        ).fetchdf()

        neighbors = defaultdict(set)
        for a, b in zip(df["group_1_id"], df["group_2_id"]):
            neighbors[a].add(b)
            neighbors[b].add(a)

        universe = store.num_groups
        module_size = len(module_ids)
        rows = []
        for group, nbrs in neighbors.items():
            if exclude_module and group in module_ids:
                continue
            hits = len(nbrs & module_ids)
            if hits == 0:
                continue
            p = float(hypergeom.sf(hits - 1, universe, module_size, len(nbrs)))
            rows.append((p, hits, len(nbrs), group))
        rows.sort(key=lambda r: (r[0], -r[1]))
        if top:
            rows = rows[:top]

        names = store.get_names_map()
        lines = ["group\thits\tneighbors\tenrichment_p\tin_module"]
        for p, hits, deg, group in rows:
            in_module = "yes" if group in module_ids else "no"
            lines.append(f"{names.get(group, group)}\t{hits}\t{deg}\t{p:.3g}\t{in_module}")

        text = "\n".join(lines)
        if output:
            with open(output, "w") as f:
                f.write(text + "\n")
            LOGGER.INFO(f"Enrichment against {module_size} module groups -> {output}")
        else:
            print(text)
            LOGGER.INFO(
                f"Ranked {len(rows)} groups by enrichment for {module_size} module groups"
            )
    except (ValueError, KeyError) as e:
        LOGGER.ERROR(str(e))
    finally:
        store.close()
