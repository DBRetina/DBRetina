#!/usr/bin/python
# -*- coding: utf-8 -*-
from __future__ import division
import sys
import _dbretina_internal as dbretina_internal
import click
from dbretina.click_context import cli
from dbretina.validators import validate_similarity_metric
import json
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import os
import dbretina.dbretina_doc_url as dbretina_doc

# Optional import for FDR correction
try:
    from statsmodels.stats.multitest import multipletests
    HAS_STATSMODELS = True
except ImportError:
    HAS_STATSMODELS = False

def plot_histogram(json_path, outout_file_path, use_log = False):
    # Load data from JSON file

    with open(json_path) as f:
        data = json.load(f)
        
    # Set up plot style and figure size
    sns.set(style="whitegrid")
    plt.rcParams.update({"font.size": 14})
    fig, ax = plt.subplots(figsize=(14, 8))

    # Colors for each metric
    palette = sns.color_palette("husl", 7)
    colors = {
        "ochiai": palette[0],
        "containment": palette[1],
        "jaccard": palette[2],
        "csi": palette[3],
        "dice": palette[4],
    }

    # Number of metrics and similarity ranges
    num_metrics = len(data)
    num_ranges = len(data[next(iter(data))])
    bar_width = 0.8 / num_metrics

    # Sort x-axis values
    x_labels = sorted(data[next(iter(data))].keys(), key=lambda x: (int(x.split('-')[0]), int(x.split('-')[1])))

    # Plot each metric
    x = np.arange(len(x_labels))
    for idx, (metric, values) in enumerate(data.items()):
        y = [values[label] for label in x_labels]
        ax.bar(x + idx * bar_width, y, color=colors[metric], width=bar_width, alpha=0.8, label=metric)

    # Set x-axis tick labels and legend
    ax.set_xticks(x + (num_metrics - 1) * bar_width / 2)
    ax.set_xticklabels(x_labels, rotation=45)
    ax.legend(loc = "upper right")

    # Set axis labels and title
    ax.set_xlabel("Similarity Range")
    ax.set_ylabel("Frequency")
    ax.set_title("Similarity Metrics Frequency Distribution", fontsize=18)

    # Set y-axis to log scale and adjust limits
    if use_log:
        ax.set_ylabel("Log Frequency")
        ax.set_yscale('log')

        # ax.set_ylim(bottom=0.5)
    # ax.set_yscale('function', functions=(lambda x: np.log2(x+0.01), lambda x: 2**x))
    # ax.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: f"{np.log2(x):.0f}" if x > 0 else "0"))
    # ax.set_ylim(bottom=0.5)

    # Save plot to a file
    plt.tight_layout()
    plt.savefig(outout_file_path, dpi=600)

def inject_index_command(index_prefix):
    return ""


def get_command():
    _sys_argv = sys.argv
    for i in range(len(_sys_argv)):
        if _sys_argv[i] == "-i":
            _sys_argv[i+1] = os.path.abspath(_sys_argv[i+1])
    return "#command: DBRetina " + " ".join(_sys_argv[1:])

def apply_fdr_correction(tsv_path, alpha=0.05, method='fdr_bh'):
    """
    Apply FDR correction to p-values in pairwise output TSV.

    Args:
        tsv_path: Path to the pairwise TSV file
        alpha: Significance threshold (default: 0.05)
        method: Correction method ('fdr_bh' for Benjamini-Hochberg)

    Returns:
        DataFrame with added qvalue and fdr_significant columns
    """
    if not HAS_STATSMODELS:
        raise ImportError(
            "statsmodels is required for FDR correction. "
            "Install with: pip install 'DBRetina[stats]' or pip install 'DBRetina[all]'"
        )

    # Read TSV, skipping comment lines
    df = pd.read_csv(tsv_path, sep='\t', comment='#')

    if 'pvalue' not in df.columns:
        raise ValueError("No 'pvalue' column found. Run pairwise with --pvalue flag first.")

    # Apply Benjamini-Hochberg correction
    _, qvalues, _, _ = multipletests(df['pvalue'].values, alpha=alpha, method=method)

    df['qvalue'] = qvalues
    df['fdr_significant'] = qvalues < alpha

    return df


@cli.command(name="pairwise", epilog=dbretina_doc.doc_url("pairwise"), help_priority=2)
@click.option('-i', '--index-prefix', required=True, type=click.STRING, help="Index file prefix")
@click.option('-t', '--threads', "user_threads", default=1, required=False, show_default=True, type=click.IntRange(min=1), help="number of cores")
@click.option('-m', '--metric', "similarity_type", required=False, default="containment", type=click.STRING, callback=validate_similarity_metric, help="select from ['containment', 'jaccard', 'ochiai']")
@click.option('-c', '--cutoff', required=False, type=click.FloatRange(0, 100, clamp=False), default=0.0, show_default=True, help="filter out similarities < cutoff")
@click.option('--pvalue', 'calculate_pvalue', is_flag=True, required = False, default = False, help="calculate Hypergeometric p-value")
@click.option('--fdr', 'apply_fdr', is_flag=True, default=False, help="Apply Benjamini-Hochberg FDR correction to p-values (requires --pvalue and statsmodels)")
@click.option('--fdr-alpha', 'fdr_alpha', default=0.05, type=float, show_default=True, help="FDR significance threshold")
@click.pass_context
def main(ctx, index_prefix, user_threads, similarity_type, cutoff, calculate_pvalue, apply_fdr, fdr_alpha):
    """
    Calculate pairwise similarities.

    STATISTICAL NOTES:

    P-VALUE: The hypergeometric p-value tests whether the observed gene overlap
    is greater than expected by chance. Low p-values indicate significant
    association. The p-value depends critically on the universe size (N), which
    is the total number of unique features in the index.

    MULTIPLE TESTING: When comparing many pairs, use --fdr to apply Benjamini-
    Hochberg correction. Without correction, at p<0.05, expect ~5% false
    positives across all comparisons.

    ODDS RATIO INTERPRETATION:
      - OR > 1: Positive association (more overlap than expected by chance)
      - OR = 1: No association (overlap equals random expectation)
      - OR < 1: Negative association (less overlap than expected)
      - OR = -1: Undefined (filtered from output)

    CONTAINMENT: Reported in the symmetric `containment` column as
    shared/min(s1,s2) (a single value per pair); directionality is not
    preserved.
    """
    
    # -i is a plain STRING (not click.Path), so existence isn't validated by
    # Click. Guard here so a missing/typo'd prefix gives a clean [ERROR] instead
    # of a raw RuntimeError traceback from the unwrapped dbretina_internal.pairwise
    # call below (issue 073).
    if not os.path.exists(f"{index_prefix}.dbri") and \
            not os.path.exists(f"{index_prefix}_raw.json"):
        ctx.obj.ERROR(f"index prefix '{index_prefix}' (.dbri / _raw.json) not found")
        sys.exit(1)

    commands = inject_index_command(index_prefix) + '\n' + get_command()

    # Validate FDR options
    if apply_fdr and not calculate_pvalue:
        ctx.obj.ERROR("--fdr requires --pvalue flag to be set")
        sys.exit(1)

    if apply_fdr and not HAS_STATSMODELS:
        ctx.obj.ERROR("FDR correction requires statsmodels. Install with: pip install 'DBRetina[stats]' or 'DBRetina[all]'")
        sys.exit(1)

    if calculate_pvalue:
        ctx.obj.INFO("Please wait for a while, calculating p-value may take a long time.")

    ctx.obj.INFO(
        f"Constructing the pairwise matrix using {user_threads} cores.")
    dbretina_internal.pairwise(index_prefix, user_threads, similarity_type, cutoff, commands, calculate_pvalue)

    # Warn users about multiple testing when p-values are computed without FDR
    if calculate_pvalue and not apply_fdr:
        ctx.obj.WARNING(
            "Computing raw p-values without FDR correction. "
            "For multiple testing correction, consider using --fdr flag."
        )

    stats_json_path = f"{index_prefix}_DBRetina_pairwise_stats.json"
    dbrp_path = f"{index_prefix}_DBRetina_pairwise.dbrp"
    tsv_path = f"{index_prefix}_DBRetina_pairwise.tsv"
    linear_histo = f"{index_prefix}_DBRetina_similarity_metrics_plot_linear.png"
    log_histo = f"{index_prefix}_DBRetina_similarity_metrics_plot_log.png"

    if os.path.exists(dbrp_path) and not os.path.exists(stats_json_path):
        # Read statistics from .dbrp binary file
        stats_json_str = dbretina_internal.dbrp_load_statistics(dbrp_path)
        stats_data = json.loads(stats_json_str)
        with open(stats_json_path, 'w') as f:
            json.dump(stats_data, f, indent=2)

    # Apply FDR correction if requested
    if apply_fdr and os.path.exists(tsv_path):
        ctx.obj.INFO(f"Applying Benjamini-Hochberg FDR correction (alpha={fdr_alpha})...")
        try:
            df = apply_fdr_correction(tsv_path, alpha=fdr_alpha)

            # Count significant pairs
            raw_sig = (df['pvalue'] < fdr_alpha).sum()
            fdr_sig = df['fdr_significant'].sum()

            ctx.obj.INFO(f"Raw significant pairs (p < {fdr_alpha}): {raw_sig}")
            ctx.obj.INFO(f"FDR-corrected significant pairs (q < {fdr_alpha}): {fdr_sig}")

            # Save corrected results
            fdr_tsv_path = f"{index_prefix}_DBRetina_pairwise_fdr.tsv"
            df.to_csv(fdr_tsv_path, sep='\t', index=False)
            ctx.obj.INFO(f"FDR-corrected results saved to: {fdr_tsv_path}")

        except Exception as e:
            ctx.obj.ERROR(f"FDR correction failed: {e}")

    ctx.obj.INFO(f"Plotting similarity metrics distribution to {linear_histo} and {log_histo}")
    plot_histogram(stats_json_path, linear_histo, use_log=False)
    plot_histogram(stats_json_path, log_histo, use_log=True)

    ctx.obj.SUCCESS("Done.")
