from __future__ import division
from collections import defaultdict
import os
import click
from kSpider2.click_context import cli
import rustworkx as rx
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import math
import networkx as nx
import plotly.graph_objects as go
import sys


def get_command():
    _sys_argv = sys.argv
    for i in range(len(_sys_argv)):
        if os.path.isfile(_sys_argv[i]):
            _sys_argv[i] = os.path.abspath(_sys_argv[i])
        if _sys_argv[i] == '-o':
            _sys_argv[i+1] = os.path.abspath(_sys_argv[i+1])
    return "DBRetina " + " ".join(_sys_argv[1:])


class Clusters:

    distance_to_col = {
        "min_cont": 5,
        "avg_cont": 6,
        "max_cont": 7,
        "ochiai": 8,
        "jaccard": 9,
    }

    seq_to_kmers = dict()
    names_map = dict()

    def __init__(self, logger_obj, pairwise_file, cut_off_threshold, dist_type, output_prefix):
        self.output_prefix = output_prefix
        self.Logger = logger_obj
        self.edges_batch_number = 10_000_000
        self.dist_type = dist_type
        self.cut_off_threshold = cut_off_threshold
        self.pairwise_file = pairwise_file
        self.shared_kmers_threshold = 200
        self.original_nodes = {}
        self.metadata = []
        self.Logger.INFO("Loading TSV pairwise file")
        if dist_type not in self.distance_to_col:
            logger_obj.ERROR("unknown distance!")
        self.dist_col = self.distance_to_col[dist_type]

        self.graph = rx.PyGraph()

        total_nodes_no = int(
            next(open(pairwise_file, 'r')).strip().split(':')[-1])
        nodes_range = range(1, total_nodes_no + 1)
        self.nodes_indeces = self.graph.add_nodes_from(list(nodes_range))

    def construct_graph(self):
        batch_counter = 0
        edges_tuples = []

        print("[i] constructing graph")
        with open(self.pairwise_file, 'r') as pairwise_tsv:

            # skip comments
            while True:
                pos = pairwise_tsv.tell()
                line = pairwise_tsv.readline()
                if not line.startswith('#'):
                    pairwise_tsv.seek(pos)
                    break
                else:
                    self.metadata.append(line)            
            self.metadata.append(f"#command: {get_command()}\n")

            next(pairwise_tsv)  # Skip header
            for row in pairwise_tsv:
                row = row.strip().split('\t')
                distance = float(row[self.dist_col])

                # don't make graph edge
                if distance < self.cut_off_threshold:
                    continue

                if batch_counter < self.edges_batch_number:
                    batch_counter += 1
                    seq1 = int(row[0]) - 1
                    seq2 = int(row[1]) - 1
                    self.original_nodes[int(row[0])] = row[2]
                    self.original_nodes[int(row[1])] = row[3]

                    edges_tuples.append((seq1, seq2, distance))
                else:
                    self.graph.add_edges_from(edges_tuples)
                    batch_counter = 0
                    edges_tuples.clear()

            if len(edges_tuples):
                self.graph.add_edges_from(edges_tuples)

    def plot_histogram2(self, cluster_sizes):
        plt.figure(figsize=(10, 6))
        sns.set_style('whitegrid')
        sns.histplot(cluster_sizes, kde=True, color='darkblue', bins=math.ceil(max(cluster_sizes)/10))
        plt.title('Histogram of Cluster Sizes', fontsize=20)
        plt.xlabel('Cluster Sizes', fontsize=15)
        plt.ylabel('Count', fontsize=15)
        plt.yscale('log')
        plt.savefig(f"{self.output_prefix}_DBRetina_clusters_{self.cut_off_threshold}%.png", dpi=500)
    
    def plot_histogram(self, cluster_sizes):
        # Set style and context to make a nicer plot
        sns.set_style("white")
        sns.set_context("talk")

        plt.figure(figsize=(10,6))  # Set the figure size
        plot = sns.histplot(cluster_sizes, kde=True, color='skyblue', binwidth=1, edgecolor='black')  # Generate histogram with KDE
        
        # Remove top and right axes spines
        sns.despine()

        plt.title('Histogram of Cluster Sizes with KDE')  # Set the title
        plt.xlabel('Cluster Sizes')  # Set the x-label
        plt.ylabel('Count')  # Set the y-label

        # Set xticks to integer values
        plt.xticks(np.arange(min(cluster_sizes), max(cluster_sizes)+1, 1.0))
        
        # Add a legend
        plot.legend(labels=['KDE', 'Cluster Sizes'])
        plt.savefig(f"{self.output_prefix}_DBRetina_clusters_{self.cut_off_threshold}%.png", dpi=500)

    def plot_bubbles(self, cluster_sizes):
         # Create a new figure
        fig, ax = plt.subplots(figsize=(10, 6))
        
        # Set the style using seaborn
        sns.set_style("whitegrid")
        
        # Create arrays with the same size as cluster_sizes with random values for the x and y axes
        x = np.random.rand(len(cluster_sizes))
        y = np.random.rand(len(cluster_sizes))

        # Scale the cluster sizes so they're more visually appealing
        # scaled_sizes = [i**2 for i in cluster_sizes]
        
        # Create a scatter plot where the size of each point is proportional to the cluster size
        scaled_sizes = [i**0.5 * 20 for i in cluster_sizes]  
        # scaled_sizes = [i / max(cluster_sizes) for i in cluster_sizes]
    
        # Normalize the cluster sizes for the color mapping
        normalized_sizes = [(i - min(cluster_sizes)) / (max(cluster_sizes) - min(cluster_sizes)) for i in cluster_sizes]



        scatter = ax.scatter(x, y, s=scaled_sizes, c=cluster_sizes, cmap='viridis', alpha=0.6, edgecolors="w", linewidth=2)

        # Add a colorbar
        cbar = plt.colorbar(scatter)
        cbar.set_label('Cluster Sizes', fontsize=15)
        
        ax.set_title('Bubble Plot of Cluster Sizes', fontsize=20)
        ax.set_xlabel('Random X', fontsize=15)
        ax.set_ylabel('Random Y', fontsize=15)
        plt.savefig(f"{self.output_prefix}_DBRetina_clusters_{self.cut_off_threshold}%_bubbles.png", dpi=500)

    def cluster_graph(self):

        cluster_sizes = []

        self.connected_components = rx.connected_components(self.graph)
        single_components = 0
        retworkx_export = f"{self.output_prefix}_DBRetina_clusters_{self.cut_off_threshold}%.tsv"
        # and {self.output} ...")
        self.Logger.INFO(f"writing {retworkx_export}")
        # rx.node_link_json(self.graph, path = retworkx_export)
        cluster_id = 1
        with open(retworkx_export, 'w') as CLUSTERS:
            for metadata_line in self.metadata:
                CLUSTERS.write(metadata_line)

            CLUSTERS.write(f"cluster_id\tcluster_size\tcluster_members\n")
            for component in self.connected_components:
                # uncomment to exclude single genome clusters from exporting
                if len(component) == 1 and list(component)[0] + 1 not in self.original_nodes:
                    continue
                named_component = [self.original_nodes[node + 1]
                                   for node in component]
                # CLUSTERS.write(cluster_id + '\t' + len(component) + '\t' + '|'.join(named_component) + '\n')
                CLUSTERS.write(
                    f"{cluster_id}\t{len(component)}\t{'|'.join(named_component)}\n")
                cluster_sizes.append(len(component))
                cluster_id += 1

        self.Logger.INFO("plotting cluster sizes histogram and bubble plot")
        self.plot_histogram(cluster_sizes)
        self.plot_bubbles(cluster_sizes)
        self.Logger.INFO(f"number of clusters: {cluster_id - 1}")


"""
TODO:
New help messages

1. containment cutoff (sim_cutoff): cluster sequences with (containment > cutoff) where containment = shared kmers % to the total kmers in the smallest node.
2. connectivity cutoff (con_cutoff): cluster sequences with (connectivity > cutoff) where connectivity = shared kmers % to the total kmers in the largest node.
3. min count cutoff (min_count): the min kmers count of a node to connect two clusters, otherwise the node will be reported twice in both clusters.
"""


@cli.command(name="cluster", help_priority=4)
@click.option('-c', '--cutoff', required=False, type=click.FloatRange(0, 100, clamp=False), default=0.0, show_default=True, help="cluster the supergroups with (distance > cutoff)")
@click.option('-o', '--output-prefix', "output_prefix", required=True, type=click.STRING, help="output file prefix")
@click.option('-p', '--pairwise', 'pairwise_file', required=False, type=click.Path(exists=True), help="filtered pairwise TSV file")
@click.option('-d', '--dist-type', "distance_type", required=True, show_default=True, type=click.STRING, help="select from ['min_cont', 'avg_cont', 'max_cont', 'ochiai', 'jaccard']")
@click.pass_context
def main(ctx, pairwise_file, cutoff, distance_type, output_prefix):
    """Graph-based clustering of the pairwise TSV file."""

    cutoff = float(cutoff)

    kCl = Clusters(logger_obj=ctx.obj, pairwise_file=pairwise_file,
                   cut_off_threshold=cutoff, dist_type=distance_type, output_prefix=output_prefix)
    ctx.obj.INFO("Building the main graph...")
    kCl.construct_graph()
    ctx.obj.INFO("Clustering...")
    kCl.cluster_graph()
