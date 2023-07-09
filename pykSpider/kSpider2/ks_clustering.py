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
import kSpider2.dbretina_doc_url as dbretina_doc

def check_if_there_is_a_pvalue(pairwise_file):
    with open(pairwise_file) as F:
        for line in F:
            if not line.startswith("#"):
                return "pvalue" in line
            else:
                continue

def get_command():
    _sys_argv = sys.argv
    for i in range(len(_sys_argv)):
        if os.path.isfile(_sys_argv[i]):
            _sys_argv[i] = os.path.abspath(_sys_argv[i])
        if _sys_argv[i] == '-o':
            _sys_argv[i+1] = os.path.abspath(_sys_argv[i+1])
    return "DBRetina " + " ".join(_sys_argv[1:])

class Clusters:

    metric_to_col = {
        "containment": 5,
        "ochiai": 6,
        "jaccard": 7,
        "odds_ratio": 8,
        "pvalue": 9,
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
        
        # replace part of the file name that shares the same prefix
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

    def __init__(self, logger_obj, pairwise_file, cut_off_threshold, metric, output_prefix, community):
        self.output_prefix = output_prefix
        self.Logger = logger_obj
        self.edges_batch_number = 10_000_000
        self.metric = metric
        self.cut_off_threshold = cut_off_threshold
        self.pairwise_file = pairwise_file
        self.shared_kmers_threshold = 200
        self.original_nodes = {}
        self.metadata = []
        self.community = community
        self.Logger.INFO("Loading TSV pairwise file")
        if metric not in self.metric_to_col:
            logger_obj.ERROR("unknown metric!")
        self.metric_col = self.metric_to_col[metric]
        
        # check if pvalue
        if metric == "pvalue" and not check_if_there_is_a_pvalue(pairwise_file):
            logger_obj.ERROR("pvalue not found in pairwise file!")

        self.graph = ig.Graph() if community else rx.PyGraph()
        self.add_edges = self._add_igraph_edges if community else self._add_rx_edges
        self.add_nodes = self._add_igraph_nodes if community else self._add_rx_nodes
        
        total_nodes_no = int(
            next(open(pairwise_file, 'r')).strip().split(':')[-1])
        nodes_range = range(1, total_nodes_no + 1)
        self.nodes_indeces = self.add_nodes(list(nodes_range))
    
    

    def construct_graph(self):
        batch_counter = 0
        edges_tuples = []

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
                similarity = float(row[self.metric_col])
                self.original_nodes[int(row[0])] = row[2]
                self.original_nodes[int(row[1])] = row[3]

                # don't make graph edge
                if similarity < self.cut_off_threshold:
                    continue

                if batch_counter < self.edges_batch_number:
                    batch_counter += 1
                    seq1 = int(row[0]) - 1
                    seq2 = int(row[1]) - 1

                    edges_tuples.append((seq1, seq2, similarity))
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
                similarity = float(row[self.metric_col])
                self.original_nodes[int(row[0])] = row[2]
                self.original_nodes[int(row[1])] = row[3]

                # don't make graph edge
                if similarity < self.cut_off_threshold:
                    continue

                if batch_counter < self.edges_batch_number:
                    batch_counter += 1
                    seq1 = int(row[0]) - 1
                    seq2 = int(row[1]) - 1

                    edges_tuples.append((seq1, seq2, similarity))
                else:
                    self.add_edges(edges_tuples)
                    batch_counter = 0
                    edges_tuples.clear()

            if len(edges_tuples):
                self.add_edges(edges_tuples)        
    
    def plot_histogram(self, cluster_sizes, filename):
        # Create a figure and a set of subplots
        fig, ax = plt.subplots()

        # Set the style of Seaborn
        sns.set(style="whitegrid")

        # Find the range of the data
        data_range = [np.min(cluster_sizes), np.max(cluster_sizes)]

        # Create the bins - one for each possible value in the range of data
        bins = np.arange(data_range[0], data_range[1] + 2) - 0.5

        # Plot histogram with Seaborn
        sns.histplot(cluster_sizes, bins=bins, color='navy', kde=False, ax=ax)

        # Add labels and title
        ax.set(xlabel='Cluster Size', ylabel='Frequency (Log Scale)', 
            title='Histogram of Cluster Sizes')

        # Use a logarithmic scale on y-axis to handle outliers
        ax.set_yscale('log')

        # Let matplotlib handle the x-axis ticks automatically
        plt.xticks(rotation=90)

        # Remove top and right borders
        sns.despine()

        # Save the plot as a png file
        fig.savefig(filename)
        plt.close()
    

    def _plot_histogram(self, cluster_sizes):
        # Set style and context to make a nicer plot
        sns.set_style("whitegrid")
        # sns.set_context("talk")

        plt.figure()  # Set the figure size
        plot = sns.histplot(
            cluster_sizes, 
            color='skyblue', 
            edgecolor='black', 
            stat='count', 
            # bins=10, 
            # discrete=False
            )

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
            # TODO either plot both igraph and rustworkx or disable both
            # ig.plot(self.connected_components, 
            #     f"{self.output_prefix}_graph_plot.png", 
            #     bbox=(1500, 1500), 
            #     vertex_label_size=5, 
            #     edge_arrow_size=0.5, 
            #     edge_arrow_width=0.5, 
            #     edge_width=0.5,
            #     edge_label_size=5,
            #     edge_curved=False,
            #     layout = 'auto')

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
        # TODO: activate only if also generated from connected components mode.
        # if self.community:
        #     ig.write(self.graph, f"{self.output_prefix}_clusters.gml", format="gml")
        #     ig.write(self.graph, f"{self.output_prefix}_clusters.graphml", format="graphml")
        

        self.Logger.INFO("plotting cluster sizes histogram and bubble plot")
        self.plot_histogram(cluster_sizes, f"{self.output_prefix}_clusters_histogram.png")
        self.Logger.INFO(f"Total number of clustered supergroups: {total_clustered_nodes}")
        self.Logger.INFO(f"Average cluser size: {total_clustered_nodes / len(self.connected_components)}")
        self.plot_bubbles(cluster_sizes)
        self.Logger.INFO(f"number of clusters: {cluster_id - 1}")


@cli.command(name="cluster", epilog = dbretina_doc.doc_url("cluster"), help_priority=4)
@click.option('-p', '--pairwise', 'pairwise_file', required=True, type=click.Path(exists=True), help="pairwise TSV file")
@click.option('-m', '--metric', "metric", required=True, type=click.STRING, help="select from ['containment', 'ochiai', 'jaccard', 'pvalue']")
@click.option("--community", "community", is_flag=True, help="clusters as communities", default=False)
@click.option('-c', '--cutoff', required=True, type=click.FloatRange(0, 100, clamp=False), default=0.0, help="cluster the supergroups with (similarity > cutoff)")
@click.option('-o', '--output-prefix', "output_prefix", required=True, type=click.STRING, help="output file prefix")
@click.pass_context
def main(ctx, pairwise_file, cutoff, metric, output_prefix, community):
    """Graph-based clustering of the pairwise TSV file."""
    
    cutoff = float(cutoff)
    kCl = Clusters(logger_obj=ctx.obj, pairwise_file=pairwise_file,
                   cut_off_threshold=cutoff, metric=metric, output_prefix=output_prefix, community=community)
    ctx.obj.INFO("Building the main graph...")
    kCl.construct_graph()
    ctx.obj.INFO("Clustering...")
    kCl.cluster_graph()
