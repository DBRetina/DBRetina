#!/usr/bin/python
# -*- coding: utf-8 -*-
from __future__ import division
import sys
import click
from kSpider2.click_context import cli
import matplotlib.pyplot as plt
import json
import os
import seaborn as sns
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm
import pandas as pd

class Interactome:
    def __init__(self, output_prefix):
        self.graph = {}
        self.output_prefix = output_prefix

    def add_edge(self, node1, node2):
        if node1 == node2:
            return
        if node1 not in self.graph:
            self.graph[node1] = {}
        if node2 not in self.graph:
            self.graph[node2] = {}
        if node2 in self.graph[node1]:
            self.graph[node1][node2] += 1
            self.graph[node2][node1] += 1
        else:
            self.graph[node1][node2] = 1
            self.graph[node2][node1] = 1

    def export(self):
        rows = ["feature_1\tfeature_2\tshared_groups"]
        for node1, edges in self.graph.items():
            rows.extend(
                f"{node1}\t{node2}\t{weight}"
                for node2, weight in edges.items()
                if node1 < node2
            )
        with open(f"{self.output_prefix}_interactome.tsv", 'w') as file:
            file.write("\n".join(rows))
    
        
        
    def plot_statistics(self):
        # Create a DataFrame from the graph data
        data = []
        for node1, edges in self.graph.items():
            data.extend(
                {'Node1': node1, 'Node2': node2, 'Weight': weight}
                for node2, weight in edges.items()
                if node1 < node2
            )
        df = pd.DataFrame(data)

        # Create a summary statistics table
        print("Summary statistics of the interactome weights:")
        summary_stats = df['Weight'].describe()
        print(pd.DataFrame.from_dict(dict(summary_stats), orient='index').to_string())

        # Plot a histogram of the weights
        plt.figure(figsize=(10, 6))
        sns.histplot(df['Weight'], kde=False, bins=30)
        plt.title('Histogram of Edge Weights')
        plt.xlabel('Edge Weight')
        plt.ylabel('Frequency')
        plt.savefig(f"{self.output_prefix}_interactome_histogram.png", dpi=500)

        # Boxplot to show spread of weights
        plt.figure(figsize=(10, 6))
        sns.boxplot(x=df['Weight'])
        plt.title('Boxplot of Edge Weights')
        plt.savefig(f"{self.output_prefix}_interactome_boxplot.png", dpi=500)

        # Scatter plot to show any potential relationship between nodes and weights
        df['NodePair'] = df.apply(lambda row: f"{row['Node1']} - {row['Node2']}", axis=1)
        plt.figure(figsize=(10, 6))
        sns.scatterplot(x='NodePair', y='Weight', data=df)
        plt.title('Scatter Plot of Node Pairs and Weights')
        plt.xticks(rotation=90)
        plt.savefig(f"{self.output_prefix}_interactome_scatter.png", dpi=500)
    
    def graph_export(self):
        self.G = nx.Graph()
        for node1, edges in self.graph.items():
            for node2, weight in edges.items():
                self.G.add_edge(node1, node2, weight=weight)

        nx.write_graphml(self.G, f"{self.output_prefix}_interactome.graphml")
        nx.write_gexf(self.G, f"{self.output_prefix}_interactome.gexf")
        nx.write_pajek(self.G, f"{self.output_prefix}_interactome.net")


def get_command():
    _sys_argv = sys.argv
    for i in range(len(_sys_argv)):
        if _sys_argv[i] == "-i":
            _sys_argv[i+1] = os.path.abspath(_sys_argv[i+1])
    return "#command: DBRetina " + " ".join(_sys_argv[1:])


@cli.command(name="interactome", help_priority=9)
# @click.option('-t', '--threads', "user_threads", default=1, required=False, type=int, help="number of cores")
@click.option('-i', '--index-prefix', "index_prefix", required=True, type=click.STRING, help="index file prefix")
@click.option('-p', '--pairwise', 'pairwise_file', required=True, type=click.Path(exists=True), help="pairwise TSV file")
@click.option('-o', '--output-prefix', "output_prefix", required=True, type=click.STRING, help="output file prefix")
@click.pass_context
def main(ctx, index_prefix, pairwise_file, output_prefix):

    ###################################
    # 1. Extract gene set pairs
    ###################################

    geneSet_pairs = set()
    metadata = []
    # iterate over the pairwise file to construct geneSet_pairs
    with open(pairwise_file, 'r') as pairwise_tsv:
        while True:
            pos = pairwise_tsv.tell()
            line = pairwise_tsv.readline()
            if not line.startswith('#'):
                pairwise_tsv.seek(pos)
                break
            else:
                metadata.append(line)            
        metadata.append(f"#command: {get_command()}\n")

        next(pairwise_tsv)
        for row in pairwise_tsv:
            row = row.strip().split('\t')
            geneSet_pairs.add((row[2], row[3]))
            
    ##############################################
    # 2. Map gene set to genes (group to features)
    ##############################################

    raw_json_file = f"{index_prefix}_raw.json"
    with open(raw_json_file, "r") as f:
        supergroups_to_features = json.loads(f.read())["data"]
    
    ##############################################
    #3. Build the interactome
    ##############################################
    
    interactome = Interactome(output_prefix)
    
    for geneSet_pair in tqdm(geneSet_pairs):
        geneSet1_features = supergroups_to_features[geneSet_pair[0]]
        geneSet2_features = supergroups_to_features[geneSet_pair[1]]
        
        # create edges between all pairs of features
        for feature1 in geneSet1_features:
            for feature2 in geneSet2_features:
                interactome.add_edge(feature1, feature2)
    
    interactome.export()
    # interactome.graph_export()
    # interactome.plot_statistics()