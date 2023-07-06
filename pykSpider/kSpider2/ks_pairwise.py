#!/usr/bin/python
# -*- coding: utf-8 -*-
from __future__ import division
import sys
import _kSpider_internal as kSpider_internal
import click
from kSpider2.click_context import cli
import json
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os
import kSpider2.dbretina_doc_url as dbretina_doc

def plot_histogram(json_path, outout_file_path, use_log = False):
    # Load data from JSON file

    with open(json_path) as f:
        data = json.load(f)
        
    # Set up plot style and figure size
    sns.set(style="whitegrid")
    plt.rcParams.update({"font.size": 14})
    fig, ax = plt.subplots(figsize=(14, 8))

    # Colors for each metric
    colors = {
        "ochiai": sns.color_palette("husl", 5)[1],
        "containment": sns.color_palette("husl", 5)[2],
        "jaccard": sns.color_palette("husl", 5)[3],
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
    ax.legend()

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
    extra_file = f"{index_prefix}.extra"
    if not os.path.exists(extra_file):
        return ""
    with open(extra_file, "r") as f:
        for line in f:
            line = line.strip().split(":")
            if line[0] == "command":
                return line[1]
        return ""


def get_command():
    _sys_argv = sys.argv
    for i in range(len(_sys_argv)):
        if _sys_argv[i] == "-i":
            _sys_argv[i+1] = os.path.abspath(_sys_argv[i+1])
    return "#command: DBRetina " + " ".join(_sys_argv[1:])

@cli.command(name="pairwise", epilog=dbretina_doc.doc_url("pairwise"), help_priority=2)
@click.option('-i', '--index-prefix', required=True, type=click.STRING, help="Index file prefix")
@click.option('-t', '--threads', "user_threads", default=1, required=False, type=int, help="number of cores")
@click.option('-m', '--metric', "similarity_type", required=False, default="containment", type=click.STRING, help="select from ['containment', 'jaccard', 'ochiai']")
@click.option('-c', '--cutoff', required=False, type=click.FloatRange(0, 100, clamp=False), default=0.0, show_default=True, help="filter out similarities < cutoff")
@click.option('--pvalue', 'calculate_pvalue', is_flag=True, required = False, default = False, help="calculate Hypergeometric p-value")
@click.pass_context
def main(ctx, index_prefix, user_threads, similarity_type, cutoff, calculate_pvalue):
    """
    Calculate pairwise similarities.
    """
    
    commands = inject_index_command(index_prefix) + '\n' + get_command()
    
    if calculate_pvalue:
        ctx.obj.INFO("Please wait for a while, calculating p-value may take a long time.")
    
    ctx.obj.INFO(
        f"Constructing the pairwise matrix using {user_threads} cores.")
    kSpider_internal.pairwise(index_prefix, user_threads, similarity_type, cutoff, commands, calculate_pvalue)
    stats_json_path = f"{index_prefix}_DBRetina_pairwise_stats.json"
    linear_histo = f"{index_prefix}_DBRetina_similarity_metrics_plot_linear.png"
    log_histo = f"{index_prefix}_DBRetina_similarity_metrics_plot_log.png"
    ctx.obj.INFO(f"Plotting similarity metrics distribution to {linear_histo} and {log_histo}")
    
    plot_histogram(stats_json_path, linear_histo, use_log=False)
    plot_histogram(stats_json_path, log_histo, use_log=True)

    ctx.obj.SUCCESS("Done.")
