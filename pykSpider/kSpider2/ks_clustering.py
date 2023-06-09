from __future__ import division
import os
import click
from kSpider2.click_context import cli
import rustworkx as rx
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import math
import sys
import igraph as ig
import leidenalg as la


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
        "odds": 10,
    }

    seq_to_kmers = dict()
    names_map = dict()
    
    def add_edges(self, edges_tuples):
        pass
    
    def add_nodes(self, nodes):
        pass

    def _add_igraph_edges(self, edges_tuples):
        weights = []
        edges = []
        for edge in edges_tuples:
            edges.append((edge[0], edge[1]))
            weights.append(edge[2])
        
        self.graph.add_edges(edges)        
        self.graph.es["weight"] = weights
    
    def _add_rx_edges(self, edges_tuples):
        self.graph.add_edges_from(edges_tuples)
        
    def _add_igraph_nodes(self, nodes):
        
        node_ids = []
        node_sizes = []
        
        featuresCount_file = self.pairwise_file.replace("pairwise.tsv", "featuresNo.tsv")
        with open(featuresCount_file) as F:
            next(F)
            for line in F:
                line = line.strip().split("\t")
                node_ids.append(line[1])
                node_sizes.append(math.log2(int(line[2])))
                # node_sizes.append(int(line[2]))
                
        self.graph.add_vertices(node_ids)
        self.graph.vs["size"] = node_sizes
    
    def _add_rx_nodes(self, nodes):
        self.graph.add_nodes_from(nodes)

    def __init__(self, logger_obj, pairwise_file, cut_off_threshold, dist_type, output_prefix, commuinty):
        self.output_prefix = output_prefix
        self.Logger = logger_obj
        self.edges_batch_number = 10_000_000
        self.dist_type = dist_type
        self.cut_off_threshold = cut_off_threshold
        self.pairwise_file = pairwise_file
        self.shared_kmers_threshold = 200
        self.original_nodes = {}
        self.metadata = []
        self.community = commuinty
        self.Logger.INFO("Loading TSV pairwise file")
        if dist_type not in self.distance_to_col:
            logger_obj.ERROR("unknown distance!")
        self.dist_col = self.distance_to_col[dist_type]

        self.graph = ig.Graph() if commuinty else rx.PyGraph()
        self.add_edges = self._add_igraph_edges if commuinty else self._add_rx_edges
        self.add_nodes = self._add_igraph_nodes if commuinty else self._add_rx_nodes
        
        total_nodes_no = int(
            next(open(pairwise_file, 'r')).strip().split(':')[-1])
        nodes_range = range(1, total_nodes_no + 1)
        self.nodes_indeces = self.add_nodes(list(nodes_range))
    
    

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
                self.original_nodes[int(row[0])] = row[2]
                self.original_nodes[int(row[1])] = row[3]

                # don't make graph edge
                if distance < self.cut_off_threshold:
                    continue

                if batch_counter < self.edges_batch_number:
                    batch_counter += 1
                    seq1 = int(row[0]) - 1
                    seq2 = int(row[1]) - 1

                    edges_tuples.append((seq1, seq2, distance))
                else:
                    self.add_edges(edges_tuples)
                    batch_counter = 0
                    edges_tuples.clear()

            if len(edges_tuples):
                self.add_edges(edges_tuples)


    def construct_igraph(self):
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
                self.original_nodes[int(row[0])] = row[2]
                self.original_nodes[int(row[1])] = row[3]

                # don't make graph edge
                if distance < self.cut_off_threshold:
                    continue

                if batch_counter < self.edges_batch_number:
                    batch_counter += 1
                    seq1 = int(row[0]) - 1
                    seq2 = int(row[1]) - 1

                    edges_tuples.append((seq1, seq2, distance))
                else:
                    self.add_edges(edges_tuples)
                    batch_counter = 0
                    edges_tuples.clear()

            if len(edges_tuples):
                self.add_edges(edges_tuples)        
    

    def plot_histogram(self, cluster_sizes):
        # Set style and context to make a nicer plot
        sns.set_style("whitegrid")
        # sns.set_context("talk")

        plt.figure()  # Set the figure size
        plot = sns.histplot(cluster_sizes, color='skyblue', edgecolor='black', stat='count', bins=10, discrete=False)  # Generate histogram with KDE
        # plot = sns.(cluster_sizes, color='skyblue', edgecolor='black', stat='count', bins=50, hue=False)  # Generate histogram with KDE

        plt.title('Histogram of Cluster Sizes')  # Set the title
        plt.xlabel('Cluster Sizes')  # Set the x-label
        plt.ylabel('Count (log scale)')  # Set the y-label
        plt.yscale('log')
        # plt.xticks(np.arange(min(cluster_sizes), max(cluster_sizes)+1, 1))

        
        # Add a legend
        # plot.legend(labels=['Cluster Sizes'])
        # plt.show()
        plt.savefig(f"{self.output_prefix}_clusters_histogram.png", dpi=500)


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
        # normalized_sizes = [(i - min(cluster_sizes)) / (max(cluster_sizes) - min(cluster_sizes)) for i in cluster_sizes]



        scatter = ax.scatter(x, y, s=scaled_sizes, c=cluster_sizes, cmap='viridis', alpha=0.6, edgecolors="w", linewidth=2)

        # Add a colorbar
        cbar = plt.colorbar(scatter)
        cbar.set_label('Cluster Sizes', fontsize=15)
        
        ax.set_title('Bubble Plot of Cluster Sizes', fontsize=20)
        ax.set_xlabel('Random X', fontsize=15)
        ax.set_ylabel('Random Y', fontsize=15)
        plt.savefig(f"{self.output_prefix}_clusters_bubbles.png", dpi=500)

    def cluster_graph(self):

        cluster_sizes = []
        
        if self.community:
            self.connected_components = la.find_partition(self.graph, 
                                                          la.CPMVertexPartition, 
                                                          weights= self.graph.es['weight'],
                                                          node_sizes = self.graph.vs['size'],
                                                          seed = 42,
                                                          )

            # self.connected_components = la.CPMVertexPartition(
            #     self.graph, 
            #     # resolution_parameter = 0.6,
            #     )

            # optimiser = la.Optimiser()
            

            # optimiser.optimise_partition(self.connected_components)
            # refine_partition = la.ModularityVertexPartition(self.graph)
            # optimiser.move_nodes_constrained(refine_partition, self.connected_components)
            ig.plot(self.connected_components, 
                f"{self.output_prefix}_communiteis.png", 
                bbox=(1500, 1500), 
                vertex_label_size=5, 
                edge_arrow_size=0.5, 
                edge_arrow_width=0.5, 
                edge_width=0.5,
                edge_label_size=5, 
                edge_curved=False, 
                layout = 'auto')

            # la.find_partition(self.graph, la.ModularityVertexPartition)
        else:
            self.connected_components = rx.connected_components(self.graph)
        
        single_components = 0
        retworkx_export = f"{self.output_prefix}_clusters.tsv"
        # and {self.output} ...")
        self.Logger.INFO(f"writing {retworkx_export}")
        # rx.node_link_json(self.graph, path = retworkx_export)
        cluster_id = 1
        total_clustered_nodes = 0
        with open(retworkx_export, 'w') as CLUSTERS:
            for metadata_line in self.metadata:
                CLUSTERS.write(metadata_line)

            CLUSTERS.write(f"cluster_id\tcluster_size\tcluster_members\n")
            for component in self.connected_components:
                # uncomment to exclude single genome clusters from exporting
                if len(component) == 1 and list(component)[0] + 1 not in self.original_nodes:
                    continue                
                # if len(component) == 1: continue
                
                named_component = [self.original_nodes[node + 1] for node in component]
                # CLUSTERS.write(cluster_id + '\t' + len(component) + '\t' + '|'.join(named_component) + '\n')
                CLUSTERS.write(f"{cluster_id}\t{len(component)}\t{'|'.join(named_component)}\n")
                cluster_sizes.append(len(component))
                total_clustered_nodes += len(component)
                cluster_id += 1

        if self.community:
            ig.write(self.graph, f"{self.output_prefix}_clusters.gml", format="gml")
            ig.write(self.graph, f"{self.output_prefix}_clusters.graphml", format="graphml")
        

        self.Logger.INFO("plotting cluster sizes histogram and bubble plot")
        self.Logger.INFO(f"Total number of clustered supergroups: {total_clustered_nodes}")
        self.Logger.INFO(f"Average cluser size: {total_clustered_nodes / len(self.connected_components)}")
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
@click.option('-p', '--pairwise', 'pairwise_file', required=False, type=click.Path(exists=True), help="filtered pairwise TSV file")
@click.option('-d', '--dist-type', "distance_type", required=True, show_default=True, type=click.STRING, help="select from ['min_cont', 'avg_cont', 'max_cont', 'ochiai', 'jaccard', 'odds_ratio]")
@click.option("--community", "community", is_flag=True, help="clusters as communities", default=False)
@click.option('-c', '--cutoff', required=False, type=click.FloatRange(0, 100, clamp=False), default=0.0, show_default=True, help="cluster the supergroups with (distance > cutoff)")
@click.option('-o', '--output-prefix', "output_prefix", required=True, type=click.STRING, help="output file prefix")
@click.pass_context
def main(ctx, pairwise_file, cutoff, distance_type, output_prefix, community):
    """Graph-based clustering of the pairwise TSV file."""
    

    cutoff = float(cutoff)

    kCl = Clusters(logger_obj=ctx.obj, pairwise_file=pairwise_file,
                   cut_off_threshold=cutoff, dist_type=distance_type, output_prefix=output_prefix, commuinty=community)
    ctx.obj.INFO("Building the main graph...")
    kCl.construct_graph()
    ctx.obj.INFO("Clustering...")
    kCl.cluster_graph()
