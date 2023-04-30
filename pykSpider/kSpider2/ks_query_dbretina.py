#!/usr/bin/python
# -*- coding: utf-8 -*-
from __future__ import division

from click.decorators import option
import _kSpider_internal as kSpider_internal
import click
from kSpider2.click_context import cli
import matplotlib.pyplot as plt
import seaborn as sns


def plot_histogram(features_counts, output_file):
    # Set style and context to make a nicer plot
    sns.set_style("white")
    sns.set_context("talk")

    plt.figure(figsize=(10, 6))  # Set the figure size
    plot = sns.histplot(features_counts, color='skyblue', edgecolor='black',
                        stat='count', bins=50)  # Generate histogram with KDE

    plt.title('Histogram of features frequencies')  # Set the title
    plt.xlabel('Feature frequency')  # Set the x-label
    plt.ylabel('Count (log scale)')  # Set the y-label
    plt.yscale('log')

    # Add a legend
    # plot.legend(labels=['Cluster Sizes'])
    plt.savefig(output_file, dpi=500)


@cli.command(name="query", help_priority=6)
@click.option('-i', '--index-prefix', "index_prefix", required=True, type=click.STRING, help="index file prefix")
@click.option('-q', '--query', "query_file", required=True, type=click.Path(exists=True), help="line-separated file of supergroups")
@click.option('-o', '--output', "output_prefix", required=True, default=None, help="output file prefix")
@click.pass_context
def main(ctx, query_file, index_prefix, output_prefix):
    """
    Query DBRetina index.
    """

    key_val_suffix = "_feature_to_groups.tsv"
    val_key_suffix = "_features_count_per_group.tsv"
    features_to_groups_file = output_prefix + key_val_suffix
    counts_file = output_prefix + val_key_suffix
    features_counts = []
    with open(counts_file) as f:
        features_counts.extend(int(line.strip().split('\t')[1]) for line in f)

    output_file = f"{output_prefix}_features_count_per_group_histogram.png"
    ctx.obj.INFO(f"Plotting histogram of features frequencies to {output_file}")
    plot_histogram(features_counts, output_file)

    kSpider_internal.query(index_prefix, query_file, output_prefix)
    ctx.obj.INFO(f"writing query results to {features_to_groups_file}, and {counts_file}")
    ctx.obj.SUCCESS("Query done!")
