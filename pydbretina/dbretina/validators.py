"""Shared validation callbacks for CLI options."""
import click

VALID_METRICS = ["containment", "ochiai", "jaccard", "csi", "dice", "odds_ratio", "pvalue"]
# Metrics accepted as the pairwise cutoff-FILTER metric (`pairwise -m`). This must
# match the C++ allowlist in src/pairwise.cpp (allowed_distances = {containment,
# ochiai, jaccard}); csi/dice ARE computed and emitted as columns, but the C++
# pairwise() rejects them as the cutoff metric with a raw std::invalid_argument.
# Keep this list == the C++ set so `pairwise -m csi/dice` fails with a clean Click
# error (exit 2) instead of an opaque C++ ValueError traceback. Other commands
# (query/cluster/export/graph) use VALID_METRICS, which still allows csi/dice.
SIMILARITY_METRICS = ["containment", "ochiai", "jaccard"]
# Metrics the `bipartite` command genuinely emits. Its output only ever carries
# containment/ochiai/jaccard columns (+ pvalue when the input has it); csi/dice/
# odds_ratio are computed by pairwise but bipartite never writes them. The
# permissive VALID_METRICS let `bipartite -m csi/dice/odds_ratio` exit 0 while
# silently producing output WITHOUT those columns (issue 048). Restrict bipartite
# -m to this set so an unsupported metric gives a clean Click error instead. This
# matches the bipartite --help text ['containment','ochiai','jaccard','pvalue'].
BIPARTITE_METRICS = ["containment", "ochiai", "jaccard", "pvalue"]
# scipy.cluster.hierarchy.linkage() accepted methods.
VALID_LINKAGE_METHODS = ["single", "complete", "average", "weighted", "centroid", "median", "ward"]

def validate_metric(ctx, param, value):
    """Validate metric is one of the known metrics."""
    if not value or value == "NA":
        return value
    if value.lower() not in VALID_METRICS:
        raise click.BadParameter(
            f"Invalid metric '{value}'. Choose from: {', '.join(VALID_METRICS)}"
        )
    return value.lower()

def validate_similarity_metric(ctx, param, value):
    """Validate metric is a similarity metric (not pvalue/odds_ratio)."""
    if value and value.lower() not in SIMILARITY_METRICS:
        raise click.BadParameter(
            f"Invalid similarity metric '{value}'. Choose from: {', '.join(SIMILARITY_METRICS)}"
        )
    return value.lower() if value else value

def validate_bipartite_metric(ctx, param, value):
    """Validate metric is one bipartite actually emits (issue 048).

    bipartite only writes containment/ochiai/jaccard (+ pvalue when present), so
    csi/dice/odds_ratio -- though valid pairwise metrics -- must be rejected here
    with a clean error rather than silently producing columnless output.
    """
    if not value or value == "NA":
        return value
    if value.lower() not in BIPARTITE_METRICS:
        raise click.BadParameter(
            f"Invalid metric '{value}'. Choose from: {', '.join(BIPARTITE_METRICS)}"
        )
    return value.lower()

def validate_linkage(ctx, param, value):
    """Validate linkage method is one accepted by scipy's linkage()."""
    if value and value.lower() not in VALID_LINKAGE_METHODS:
        raise click.BadParameter(
            f"Invalid linkage '{value}'. Choose from: {', '.join(VALID_LINKAGE_METHODS)}"
        )
    return value.lower() if value else value
