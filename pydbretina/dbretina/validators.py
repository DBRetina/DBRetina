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

def validate_linkage(ctx, param, value):
    """Validate linkage method is one accepted by scipy's linkage()."""
    if value and value.lower() not in VALID_LINKAGE_METHODS:
        raise click.BadParameter(
            f"Invalid linkage '{value}'. Choose from: {', '.join(VALID_LINKAGE_METHODS)}"
        )
    return value.lower() if value else value
