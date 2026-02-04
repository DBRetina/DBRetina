"""DBRetina genescore command — score genes by importance."""

import os
import sys

import click
from dbretina.click_context import cli
import dbretina.dbretina_doc_url as dbretina_doc


@cli.command(name="genescore", epilog=dbretina_doc.doc_url("genescore"), help_priority=17)
@click.option(
    "-d", "--data", "data_path", required=True,
    type=click.Path(exists=True),
    help="Path to pairwise Parquet directory or TSV file",
)
@click.option(
    "-i", "--index-prefix", "index_prefix", required=True,
    type=str, help="Index file prefix (for .dbri gene data)",
)
@click.argument("group_name")
@click.option(
    "--compare", "group_b", default=None, type=str,
    help="Second group name for explain-pair mode",
)
@click.option(
    "--method", default="projection", show_default=True,
    type=click.Choice(["edge_weighted", "hypergraph", "projection"],
                      case_sensitive=False),
    help="Scoring method",
)
@click.option(
    "--hops", default=2, show_default=True, type=int,
    help="Neighborhood hops (hypergraph/projection methods)",
)
@click.option(
    "-m", "--metric", default="ochiai", show_default=True, type=str,
    help="Metric for edge_weighted method",
)
@click.option(
    "-c", "--cutoff", default=20.0, show_default=True, type=float,
    help="Metric cutoff for edge_weighted method",
)
@click.option(
    "-n", "--top", default=30, show_default=True, type=int,
    help="Number of top genes to show",
)
@click.option(
    "-o", "--output", default=None, type=str,
    help="Output file path (default: print to stdout)",
)
@click.pass_context
def main(ctx, data_path, index_prefix, group_name, group_b, method,
         hops, metric, cutoff, top, output):
    """Score genes by importance for a disease/group.

    Three scoring methods are available:

    \b
      edge_weighted  — Sum of edge weights per gene
      hypergraph     — TF-IDF: local enrichment x inverse global frequency
      projection     — PageRank on gene co-occurrence graph

    \b
    Hub genes for a disease:
        DBRetina genescore -d experiment_pairwise -i myindex \\
            "autism spectrum disorders"

    \b
    Explain shared genes between two diseases:
        DBRetina genescore -d experiment_pairwise -i myindex \\
            "autism spectrum disorders" --compare "schizophrenia"

    \b
    Use hypergraph method:
        DBRetina genescore -d experiment_pairwise -i myindex \\
            "autism spectrum disorders" --method hypergraph
    """
    LOGGER = ctx.obj

    from dbretina.compat import open_pairwise
    from dbretina.pairwise_store import PairwiseStore
    from dbretina.gene_importance import GeneImportance

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
        gi = GeneImportance(store)

        if group_b:
            # Explain-pair mode
            LOGGER.INFO(
                f'Scoring shared genes between "{group_name}" and '
                f'"{group_b}" (method={method})'
            )
            df = gi.explain_pair(group_name, group_b, method=method)

            lines = ["gene\tscore\tglobal_freq\tspecificity"]
            for _, row in df.head(top).iterrows():
                lines.append(
                    f"{row['gene']}\t{row['score']:.6f}\t"
                    f"{row['global_freq']}\t{row['specificity']:.6f}"
                )

            text = "\n".join(lines)

            if output:
                # Write full results, not just top N
                full_lines = ["gene\tscore\tglobal_freq\tspecificity"]
                for _, row in df.iterrows():
                    full_lines.append(
                        f"{row['gene']}\t{row['score']:.6f}\t"
                        f"{row['global_freq']}\t{row['specificity']:.6f}"
                    )
                with open(output, "w") as f:
                    f.write("\n".join(full_lines) + "\n")
                LOGGER.INFO(f"{len(df)} shared genes scored -> {output}")
            else:
                print(text)
                LOGGER.INFO(f"{len(df)} shared genes scored")
        else:
            # Hub genes mode
            LOGGER.INFO(
                f'Scoring genes for "{group_name}" '
                f"(method={method}, hops={hops})"
            )
            df = gi.hub_genes(
                group_name, hops=hops, method=method, top_n=top,
            )

            lines = ["gene\tscore\tnum_diseases"]
            for _, row in df.iterrows():
                lines.append(
                    f"{row['gene']}\t{row['score']:.6f}\t{row['num_diseases']}"
                )

            text = "\n".join(lines)

            if output:
                with open(output, "w") as f:
                    f.write(text + "\n")
                LOGGER.INFO(f"Top {len(df)} genes -> {output}")
            else:
                print(text)
                LOGGER.INFO(f"Top {len(df)} genes scored")
    except KeyError as e:
        LOGGER.ERROR(str(e))
    except ImportError as e:
        LOGGER.ERROR(str(e))
    finally:
        store.close()
