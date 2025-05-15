#!/usr/bin/python
# -*- coding: utf-8 -*-
from __future__ import division

from click.decorators import option
import _kSpider_internal as kSpider_internal
import click
from kSpider2.click_context import cli
import matplotlib.pyplot as plt
import seaborn as sns
import os
from collections import defaultdict
import csv
import kSpider2.dbretina_doc_url as dbretina_doc
import json

class Graph:
    
    # node_to_edges dict with default = 0
    node_to_edges = defaultdict(lambda: 0)
    node_to_size = {}
    all_groups = set()
    
    def __init__(self, index_prefix):
        self.parent = {}
        self.components = None
        self.index_prefix = index_prefix
        self.load_raw_json()
        

    def load_raw_json(self):
        raw_json_file = f"{self.index_prefix}_raw.json"
        with open(raw_json_file) as f:
            group_to_items = json.load(f)["data"]
        
        for group, items in group_to_items.items():
            self.node_to_size[group] = len(items)
            self.all_groups.add(group)
    
    def get_all_groups_set(self):
        return self.all_groups


    def find(self, x):
        if x not in self.parent:
            self.parent[x] = x
        if self.parent[x] != x:
            self.parent[x] = self.find(self.parent[x])
        return self.parent[x]

    def union(self, x, y):
        self.node_to_edges[x] += 1
        self.node_to_edges[y] += 1
        self.parent[self.find(x)] = self.find(y)

    def add_edge(self, x, y):
        self.union(x, y)

    def _compute_components(self):
        self.components = {}
        for node in self.parent:
            parent = self.find(node)
            if parent not in self.components:
                self.components[parent] = []
            self.components[parent].append(node)

    def get_connected_components(self):
        if self.components is None:
            self._compute_components()
        return self.components.values()
    
    def select_node_from_cluster(self, cluster):
        if len(cluster) == 1:
            return cluster[0]

        # select the node with the most edges, if tie, select the largest node
        return max(cluster, key=lambda x: (self.node_to_edges[x], self.node_to_size[x]))


def path_to_absolute_path(ctx, param, value):
    return value if value == "NA" else os.path.abspath(value)


@cli.command(name="dedup", epilog = dbretina_doc.doc_url("dedup"), help_priority=7)
@click.option('-i', '--index-prefix', 'index_prefix', required=True, type=click.STRING, help="Index file prefix")
@click.option('-p', '--pairwise', 'pairwise_file', callback=path_to_absolute_path, required=True, type=click.Path(exists=True), help="the pairwise TSV file")
@click.option('-c', '--cutoff', 'cutoff', required=True, type=click.FloatRange(0, 100, clamp=False), help="ochiai similarity cutoff")
@click.option('-o', '--output', "output_prefix", required=True, default=None, help="output file prefix")
@click.pass_context
def main(ctx, pairwise_file, cutoff, output_prefix, index_prefix):
    """
        Deduplicate the pairwise similarity file using ochiai similarity
    """

    LOGGER = ctx.obj
    ochiai_graph = Graph(index_prefix)


    #################################
    # 1. Extract all gene sets
    #################################
    all_groups = ochiai_graph.get_all_groups_set()
    original_groups_count = len(all_groups)
    
    #################################
    # 2. Pairwise file parsing
    #################################

    metric_to_col = {
        "containment": 5,
        "ochiai": 6,
        "jaccard": 7,
        "odds_ratio": 8,
        "pvalue": 9,
    }

    METRIC_COL = metric_to_col["ochiai"]


    LOGGER.INFO(f"parsing the pairwise file: {pairwise_file}")
    with open(pairwise_file, 'r') as f:
        # Skip comment lines
        while True:
            pos = f.tell()  # remember the position
            line = f.readline()
            if not line.startswith("#"):
                f.seek(pos)  # rewind to the position before the line
                break

        # skip the header
        next(f)

        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            ochiai_similarity = float(row[METRIC_COL])
            if ochiai_similarity >= cutoff:
                # Extract the gene_set names from columns 3 and 4
                gene_set1 = row[2]
                gene_set2 = row[3]

                ochiai_graph.add_edge(gene_set1, gene_set2)

    # Get the connected components
    LOGGER.INFO("extracting the connected components")
    connected_components = ochiai_graph.get_connected_components()

    if len(connected_components) == 0:
        LOGGER.ERROR(f"There is no at least one pair of groups with similarity >= {cutoff}")

    LOGGER.INFO(f"number of connected components: {len(connected_components)}")

    for cluster in connected_components:
        selected_node = ochiai_graph.select_node_from_cluster(cluster)
        unselected_nodes = set(cluster) - {selected_node}
        all_groups -= unselected_nodes


    LOGGER.INFO(f"original number of gene sets: {original_groups_count}")
    LOGGER.INFO(f"Number of gene sets after deduplication: {len(all_groups)}")
    LOGGER.INFO(f"Number of removed gene sets due to deduplication: {original_groups_count - len(all_groups)}")

    final_output = f"{output_prefix}_deduplicated_groups.txt"

    LOGGER.INFO(f"writing deduplicated gene sets to {final_output}")
    with open(final_output, 'w') as f:
        f.write('\n'.join(all_groups))
