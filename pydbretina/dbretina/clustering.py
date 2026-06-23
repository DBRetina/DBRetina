from __future__ import division
import os
import click
from dbretina.click_context import cli
from dbretina.validators import validate_metric
from dbretina.pairwise_store import passes_cutoff
import rustworkx as rx
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import math
import sys
import igraph as ig
import leidenalg as la
import dbretina.dbretina_doc_url as dbretina_doc
import _dbretina_internal as dbretina_internal

def get_command():
    _sys_argv = sys.argv
    for i in range(len(_sys_argv)):
        if os.path.isfile(_sys_argv[i]):
            _sys_argv[i] = os.path.abspath(_sys_argv[i])
        if _sys_argv[i] == '-o':
            _sys_argv[i+1] = os.path.abspath(_sys_argv[i+1])
    return "DBRetina " + " ".join(_sys_argv[1:])


def resolve_features_no_file(pairwise_file):
    """Resolve the sibling ``*_DBRetina_featuresNo.tsv`` for a pairwise input.

    -p may be the pairwise TSV, the parquet directory, or the .dbrp binary;
    pairwise always emits ``<prefix>_DBRetina_featuresNo.tsv`` next to
    ``<prefix>_DBRetina_pairwise{,.tsv,.dbrp}``. The community node loader used a
    brittle ``replace("pairwise.tsv", "featuresNo.tsv")`` that no-ops for the
    directory/.dbrp forms and then open()s the directory (issue 016). Strip the
    pairwise extension and the ``_DBRetina_pairwise`` suffix to find the sibling.
    """
    base = pairwise_file
    for ext in (".tsv", ".dbrp"):
        if base.endswith(ext):
            base = base[:-len(ext)]
            break
    if base.endswith("_DBRetina_pairwise"):
        prefix = base[:-len("_DBRetina_pairwise")]
        return prefix + "_DBRetina_featuresNo.tsv"
    # Legacy / unrecognized naming: fall back to the original .tsv substring swap.
    return pairwise_file.replace("pairwise.tsv", "featuresNo.tsv")

class Clusters:

    metric_to_col = {
        "containment": 5,
        "ochiai": 6,
        "jaccard": 7,
        "csi": 8,
        "dice": 9,
        "odds_ratio": 10,
        "pvalue": 11,
    }

    metric_name_to_id = {
        "containment": 0,
        "ochiai": 1,
        "jaccard": 2,
        "csi": 3,
        "dice": 4,
        "odds_ratio": 5,
        "pvalue": 6,
    }

    # NOTE: per-run *mutable* attributes (group_to_features, names_map) are
    # initialized as INSTANCE attributes in __init__, NOT here. Declaring them at
    # class scope made them shared across instances, so a second Clusters() in the
    # same process inherited the first run's entries (issue 065, sibling of 045).
    # Only the read-only metric lookups above live at class scope.

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

        # Sibling featuresNo.tsv, resolved for the .tsv / parquet-dir / .dbrp forms
        # (issue 016). The file holds one line per group ordered by sequential ID,
        # so node sizes align with the gid-1 edge indices used in construct_graph.
        featuresCount_file = resolve_features_no_file(self.pairwise_file)
        if not os.path.isfile(featuresCount_file):
            self.Logger.ERROR(
                f"could not find the per-group feature counts file "
                f"'{featuresCount_file}' required for --community clustering.")
        with open(featuresCount_file) as F:
            next(F)
            for line in F:
                line = line.strip().split("\t")
                node_ids.append(line[1])
                size = int(line[2])
                # Apply configurable node weight transform
                if self.node_weight_transform == 'log2':
                    weight = math.log2(size) if size > 1 else 1.0
                elif self.node_weight_transform == 'sqrt':
                    weight = math.sqrt(size)
                elif self.node_weight_transform == 'linear':
                    weight = float(size)
                else:
                    weight = math.log2(size) if size > 1 else 1.0
                node_sizes.append(weight)

        self.graph.add_vertices(node_ids)
        self.graph.vs["size"] = node_sizes
    
    def _add_rx_nodes(self, nodes):
        self.graph.add_nodes_from(nodes)

    def __init__(self, logger_obj, pairwise_file, cut_off_threshold, metric, output_prefix, community, resolution=1.0, node_weight_transform='log2'):
        # Per-run mutable state, instance-scoped so repeated in-process
        # instantiations don't share or leak entries (issue 065).
        self.group_to_features = dict()
        self.names_map = dict()
        self.output_prefix = output_prefix
        self.Logger = logger_obj
        self.edges_batch_number = 10_000_000
        self.metric = metric
        self.cut_off_threshold = cut_off_threshold
        self.pairwise_file = pairwise_file
        self.shared_features_threshold = 200
        self.original_nodes = {}
        self.metadata = []
        self.community = community
        self.resolution = resolution
        self.node_weight_transform = node_weight_transform
        self.Logger.INFO("Loading TSV pairwise file")
        if metric not in self.metric_to_col:
            logger_obj.ERROR("unknown metric!")
        self.metric_col = self.metric_to_col[metric]
        
        # check if pvalue (format-aware: handles the .tsv / parquet-dir / .dbrp
        # forms; the old text open() crashed with IsADirectoryError on a parquet
        # directory, issue 053).
        from dbretina.compat import pairwise_has_pvalue
        if metric == "pvalue" and not pairwise_has_pvalue(pairwise_file):
            logger_obj.ERROR("pvalue not found in pairwise file!")

        self.graph = ig.Graph() if community else rx.PyGraph()
        self.add_edges = self._add_igraph_edges if community else self._add_rx_edges
        self.add_nodes = self._add_igraph_nodes if community else self._add_rx_nodes

        # Resolve the canonical .dbrp sibling from any -p form (tsv / parquet
        # dir / .dbrp); the old str.replace no-op'd for the dir/.dbrp forms and
        # fed the directory path into the binary reader -> "Invalid .dbrp file
        # (bad magic bytes)" (issue 051, same root cause as 046). Used here for
        # the node count and again in construct_graph for the edges.
        from dbretina.compat import resolve_dbrp_path
        self.dbrp_path = resolve_dbrp_path(pairwise_file)
        # Try PairwiseStore for node count
        try:
            from dbretina.compat import open_pairwise
            _store = open_pairwise(pairwise_file)
            if _store is not None:
                total_nodes_no = _store.num_groups
                _store.close()
            else:
                raise ValueError("no store")
        except Exception:
            if self.dbrp_path is not None:
                total_nodes_no = dbretina_internal.dbrp_get_num_groups(self.dbrp_path)
            elif os.path.isdir(pairwise_file):
                # No usable parquet store (no manifest.json) and no sibling
                # .dbrp: the TSV header read below would IsADirectoryError.
                logger_obj.ERROR(
                    f"'{pairwise_file}' is a pairwise directory without a usable "
                    f"parquet store (no manifest.json) and no sibling .dbrp; pass "
                    f"the pairwise TSV, its parquet directory, or the .dbrp to -p.")
            else:
                # Bare-TSV fallback: scan the comment header for the '#nodes:N'
                # line. It used to be the first line, but the .dbri-format commit
                # moved it below '# DBRetina pairwise output' (now line 5), so the
                # old `int(next(open(...))...)` parsed the wrong line and crashed
                # with "invalid literal for int()" (issue 044). Scan all comment
                # lines instead of assuming a position.
                total_nodes_no = None
                with open(pairwise_file, 'r') as _pw:
                    for _line in _pw:
                        if not _line.startswith('#'):
                            break  # past the header; stop before data rows
                        if _line.startswith('#nodes:'):
                            total_nodes_no = int(_line.strip().split(':')[-1])
                            break
                if total_nodes_no is None:
                    logger_obj.ERROR(
                        f"could not find the '#nodes:N' header line in "
                        f"'{pairwise_file}'; is it a valid pairwise TSV?")
        nodes_range = range(1, total_nodes_no + 1)
        self.nodes_indeces = self.add_nodes(list(nodes_range))
    
    

    def construct_graph(self):
        batch_counter = 0
        edges_tuples = []

        # Try Parquet/PairwiseStore first
        try:
            from dbretina.compat import open_pairwise
            store = open_pairwise(self.pairwise_file)
        except Exception:
            store = None

        if store is not None:
            self.metadata.append(f"#command: {get_command()}\n")
            names_map = store.get_names_map()
            df = store.to_pandas(
                metric=self.metric, cutoff=self.cut_off_threshold,
                columns=["group_1_id", "group_2_id", self.metric]
            )
            for _, row in df.iterrows():
                gid1 = int(row["group_1_id"])
                gid2 = int(row["group_2_id"])
                self.original_nodes[gid1] = names_map.get(gid1, "")
                self.original_nodes[gid2] = names_map.get(gid2, "")
                similarity = float(row[self.metric])
                seq1 = gid1 - 1
                seq2 = gid2 - 1
                if batch_counter < self.edges_batch_number:
                    batch_counter += 1
                    edges_tuples.append((seq1, seq2, similarity))
                else:
                    self.add_edges(edges_tuples)
                    batch_counter = 0
                    edges_tuples.clear()
            if len(edges_tuples):
                self.add_edges(edges_tuples)
            store.close()
            return
        # Fallback: existing .dbrp / TSV code. self.dbrp_path is resolve_dbrp_path's
        # result (an existing file or None), so a directory never reaches the
        # binary reader here (issue 051).
        elif self.dbrp_path is not None:
            self.metadata.append(f"#command: {get_command()}\n")
            metric_id = self.metric_name_to_id[self.metric]
            records = dbretina_internal.dbrp_filter_pairs(self.dbrp_path, metric_id, self.cut_off_threshold)
            for record in records:
                self.original_nodes[record["group_1_id"]] = record["group_1_name"]
                self.original_nodes[record["group_2_id"]] = record["group_2_name"]
                similarity = record[self.metric]
                seq1 = record["group_1_id"] - 1
                seq2 = record["group_2_id"] - 1

                if batch_counter < self.edges_batch_number:
                    batch_counter += 1
                    edges_tuples.append((seq1, seq2, similarity))
                else:
                    self.add_edges(edges_tuples)
                    batch_counter = 0
                    edges_tuples.clear()

            if len(edges_tuples):
                self.add_edges(edges_tuples)
        else:
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

                    # don't make graph edge. Metric-aware: pvalue keeps
                    # value <= cutoff, similarity metrics keep >= cutoff (068).
                    if not passes_cutoff(similarity, self.cut_off_threshold, self.metric):
                        continue

                    # Register nodes only for edges that PASS the cutoff, matching
                    # the store/.dbrp paths above (which iterate cutoff-filtered
                    # pairs). Populating original_nodes from every pre-cutoff row
                    # made the bare-TSV path emit singleton clusters for isolated
                    # nodes that the store/.dbrp paths omit (issue 063).
                    self.original_nodes[int(row[0])] = row[2]
                    self.original_nodes[int(row[1])] = row[3]

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

        # Try Parquet/PairwiseStore first
        try:
            from dbretina.compat import open_pairwise
            store = open_pairwise(self.pairwise_file)
        except Exception:
            store = None

        if store is not None:
            self.metadata.append(f"#command: {get_command()}\n")
            names_map = store.get_names_map()
            df = store.to_pandas(
                metric=self.metric, cutoff=self.cut_off_threshold,
                columns=["group_1_id", "group_2_id", self.metric]
            )
            for _, row in df.iterrows():
                gid1 = int(row["group_1_id"])
                gid2 = int(row["group_2_id"])
                self.original_nodes[gid1] = names_map.get(gid1, "")
                self.original_nodes[gid2] = names_map.get(gid2, "")
                similarity = float(row[self.metric])
                seq1 = gid1 - 1
                seq2 = gid2 - 1
                if batch_counter < self.edges_batch_number:
                    batch_counter += 1
                    edges_tuples.append((seq1, seq2, similarity))
                else:
                    self.add_edges(edges_tuples)
                    batch_counter = 0
                    edges_tuples.clear()
            if len(edges_tuples):
                self.add_edges(edges_tuples)
            store.close()
            return
        # Fallback: existing .dbrp / TSV code. self.dbrp_path is resolve_dbrp_path's
        # result (an existing file or None), so a directory never reaches the
        # binary reader here (issue 051).
        elif self.dbrp_path is not None:
            self.metadata.append(f"#command: {get_command()}\n")
            metric_id = self.metric_name_to_id[self.metric]
            records = dbretina_internal.dbrp_filter_pairs(self.dbrp_path, metric_id, self.cut_off_threshold)
            for record in records:
                self.original_nodes[record["group_1_id"]] = record["group_1_name"]
                self.original_nodes[record["group_2_id"]] = record["group_2_name"]
                similarity = record[self.metric]
                seq1 = record["group_1_id"] - 1
                seq2 = record["group_2_id"] - 1

                if batch_counter < self.edges_batch_number:
                    batch_counter += 1
                    edges_tuples.append((seq1, seq2, similarity))
                else:
                    self.add_edges(edges_tuples)
                    batch_counter = 0
                    edges_tuples.clear()

            if len(edges_tuples):
                self.add_edges(edges_tuples)
        else:
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

                    # don't make graph edge. Metric-aware: pvalue keeps
                    # value <= cutoff, similarity metrics keep >= cutoff (068).
                    if not passes_cutoff(similarity, self.cut_off_threshold, self.metric):
                        continue

                    # Register nodes only for edges that PASS the cutoff, matching
                    # the store/.dbrp paths above (which iterate cutoff-filtered
                    # pairs). Populating original_nodes from every pre-cutoff row
                    # made the bare-TSV path emit singleton clusters for isolated
                    # nodes that the store/.dbrp paths omit (issue 063).
                    self.original_nodes[int(row[0])] = row[2]
                    self.original_nodes[int(row[1])] = row[3]

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
        # np.min/np.max have no identity on an empty array; skip the plot when
        # there are no clusters rather than crashing.
        if len(cluster_sizes) == 0:
            self.Logger.WARNING("no clusters found; skipping cluster-size histogram")
            return

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

        # No edge passed the cutoff -> edge-less graph. Both the community
        # (Leiden) path (the 'weight' edge attribute is never created) and the
        # connected-components path (empty cluster_sizes -> np.min([])) crash
        # here. Degrade gracefully: emit a clean warning and write a valid
        # header-only clusters file with no clusters.
        # igraph (community) exposes ecount(); rustworkx (cc) exposes num_edges().
        edge_count = self.graph.ecount() if self.community else self.graph.num_edges()
        if edge_count == 0:
            self.Logger.WARNING(
                f"no pairs passed the cutoff ({self.cut_off_threshold}) for "
                f"metric '{self.metric}'; no clusters"
            )
            self.connected_components = []
            retworkx_export = f"{self.output_prefix}_clusters.tsv"
            self.Logger.INFO(f"writing {retworkx_export}")
            with open(retworkx_export, 'w') as CLUSTERS:
                for metadata_line in self.metadata:
                    CLUSTERS.write(metadata_line)
                CLUSTERS.write("cluster_id\tcluster_size\tcluster_members\n")
            self.Logger.INFO("Total number of clustered supergroups: 0")
            self.Logger.INFO("number of clusters: 0")
            return

        if self.community:
            # 'linear' feeds raw gene-set sizes into the CPM node_sizes; the size penalty
            # then dominates and typically blocks all merging at the default --resolution
            # (log2/sqrt compress the sizes). Warn rather than silently return all-singletons.
            if self.node_weight_transform == 'linear':
                self.Logger.WARNING(
                    "--node-weight-transform 'linear' uses raw gene-set sizes as CPM node "
                    "weights, which often yields one cluster per node at the default "
                    "--resolution; tune --resolution (lower) or use 'log2'/'sqrt'."
                )
            # Use Leiden algorithm with CPM (Constant Potts Model) partition
            # Resolution parameter controls cluster granularity:
            #   - Higher resolution = more, smaller clusters
            #   - Lower resolution = fewer, larger clusters
            self.connected_components = la.find_partition(
                self.graph,
                la.CPMVertexPartition,
                resolution_parameter=self.resolution,
                weights=self.graph.es['weight'],
                node_sizes=self.graph.vs['size'],
                seed=42,
            )
        else:
            self.connected_components = rx.connected_components(self.graph)
        
        single_components = 0
        retworkx_export = f"{self.output_prefix}_clusters.tsv"
        self.Logger.INFO(f"writing {retworkx_export}")
        cluster_id = 1
        total_clustered_nodes = 0
        with open(retworkx_export, 'w') as CLUSTERS:
            for metadata_line in self.metadata:
                CLUSTERS.write(metadata_line)

            CLUSTERS.write(f"cluster_id\tcluster_size\tcluster_members\n")
            for component in self.connected_components:
                # Skip isolated nodes that aren't in original data
                if len(component) == 1 and list(component)[0] + 1 not in self.original_nodes:
                    continue

                named_component = [self.original_nodes[node + 1] for node in component]
                CLUSTERS.write(f"{cluster_id}\t{len(component)}\t{'|'.join(named_component)}\n")
                cluster_sizes.append(len(component))
                total_clustered_nodes += len(component)
                cluster_id += 1

        self.Logger.INFO("plotting cluster sizes histogram and bubble plot")
        self.plot_histogram(cluster_sizes, f"{self.output_prefix}_clusters_histogram.png")
        self.Logger.INFO(f"Total number of clustered supergroups: {total_clustered_nodes}")
        self.Logger.INFO(f"Average cluser size: {total_clustered_nodes / len(self.connected_components)}")
        self.plot_bubbles(cluster_sizes)
        self.Logger.INFO(f"number of clusters: {cluster_id - 1}")


@cli.command(name="cluster", epilog = dbretina_doc.doc_url("cluster"), help_priority=4)
@click.option('-p', '--pairwise', 'pairwise_file', required=True, type=click.Path(exists=True), help="pairwise TSV file")
@click.option('-m', '--metric', "metric", required=True, type=click.STRING, callback=validate_metric, help="select from ['containment', 'ochiai', 'jaccard', 'pvalue']")
@click.option("--community", "community", is_flag=True, help="clusters as communities", default=False)
@click.option('-c', '--cutoff', required=True, type=click.FloatRange(0, 100, clamp=False), default=0.0, help="cluster the supergroups with (similarity > cutoff)")
@click.option('-o', '--output-prefix', "output_prefix", required=True, type=click.STRING, help="output file prefix")
@click.option('--resolution', 'resolution', default=1.0, type=float, show_default=True,
              help="Resolution parameter for Leiden CPM partition (higher = more clusters)")
@click.option('--node-weight-transform', 'node_weight_transform',
              type=click.Choice(['log2', 'linear', 'sqrt']), default='log2', show_default=True,
              help="Transform for node weights based on gene set size. 'linear' (raw size) "
                   "usually needs a tuned --resolution or it leaves every node its own cluster.")
@click.pass_context
def main(ctx, pairwise_file, cutoff, metric, output_prefix, community, resolution, node_weight_transform):
    """Graph-based clustering of the pairwise TSV file.

    RESOLUTION PARAMETER (--resolution):
      Controls cluster granularity when using --community flag.
      - Higher values (e.g., 2.0) produce more, smaller clusters.
      - Lower values (e.g., 0.5) produce fewer, larger clusters.
      - Default 1.0 is a reasonable starting point.

    NODE WEIGHT TRANSFORM (--node-weight-transform):
      Controls how gene set size affects node importance in community detection.
      - log2: Logarithmic transform (default, reduces impact of very large sets)
      - sqrt: Square root transform (moderate reduction)
      - linear: No transform (larger sets have proportionally more influence)
    """

    cutoff = float(cutoff)
    clusters = Clusters(
        logger_obj=ctx.obj,
        pairwise_file=pairwise_file,
        cut_off_threshold=cutoff,
        metric=metric,
        output_prefix=output_prefix,
        community=community,
        resolution=resolution,
        node_weight_transform=node_weight_transform
    )
    ctx.obj.INFO("Building the main graph...")
    clusters.construct_graph()
    ctx.obj.INFO(f"Clustering (resolution={resolution}, node_weight_transform={node_weight_transform})...")
    clusters.cluster_graph()
