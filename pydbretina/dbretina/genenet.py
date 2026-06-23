#!/usr/bin/python
# -*- coding: utf-8 -*-
from __future__ import division
import sys
import click
from dbretina.click_context import cli
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
import dbretina.dbretina_doc_url as dbretina_doc

class GeneNet:
    def __init__(self, output_prefix, net_kind="genenet"):
        self.graph = {}
        self.output_prefix = output_prefix
        # net_kind drives output-file suffix and log strings so the same code
        # can serve both the 'genenet' and 'interactome' commands without
        # mislabeling one as the other.
        self.net_kind = net_kind

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

    def export(self, ctx):
        rows = ["feature_1\tfeature_2\tshared_groups"]
        for node1, edges in self.graph.items():
            rows.extend(
                f"{node1}\t{node2}\t{weight}"
                for node2, weight in edges.items()
                if node1 < node2
            )
            
        output_file_name = f"{self.output_prefix}_{self.net_kind}.tsv"
        ctx.obj.INFO(f"Exporting {self.net_kind} to {output_file_name}")
        with open(output_file_name, 'w') as file:
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
        print(f"Summary statistics of the {self.net_kind} weights:")
        summary_stats = df['Weight'].describe()
        print(pd.DataFrame.from_dict(dict(summary_stats), orient='index').to_string())

        # Plot a histogram of the weights
        plt.figure(figsize=(10, 6))
        sns.histplot(df['Weight'], kde=False, bins=30)
        plt.title('Histogram of Edge Weights')
        plt.xlabel('Edge Weight')
        plt.ylabel('Frequency')
        plt.savefig(f"{self.output_prefix}_{self.net_kind}_histogram.png", dpi=500)

        # Boxplot to show spread of weights
        plt.figure(figsize=(10, 6))
        sns.boxplot(x=df['Weight'])
        plt.title('Boxplot of Edge Weights')
        plt.savefig(f"{self.output_prefix}_{self.net_kind}_boxplot.png", dpi=500)

        # Scatter plot to show any potential relationship between nodes and weights
        df['NodePair'] = df.apply(lambda row: f"{row['Node1']} - {row['Node2']}", axis=1)
        plt.figure(figsize=(10, 6))
        sns.scatterplot(x='NodePair', y='Weight', data=df)
        plt.title('Scatter Plot of Node Pairs and Weights')
        plt.xticks(rotation=90)
        plt.savefig(f"{self.output_prefix}_{self.net_kind}_scatter.png", dpi=500)
    
    def graph_export(self, graphml, gexf, ctx):
        if not graphml and not gexf:
            return

        self.G = nx.Graph()
        for node1, edges in self.graph.items():
            for node2, weight in edges.items():
                self.G.add_edge(node1, node2, weight=weight)

        # uppercase graph nodes
        self.G = nx.relabel_nodes(self.G, lambda x: x.upper())

        if graphml:
            ctx.obj.INFO(f"Exporting {self.net_kind} as graphml file")
            nx.write_graphml(self.G, f"{self.output_prefix}_{self.net_kind}.graphml")
            nx.write_gml(self.G, f"{self.output_prefix}_{self.net_kind}.gml")
            nx.write_weighted_edgelist(self.G, f"{self.output_prefix}_{self.net_kind}_weighted.edgelist")
            nx.write_graphml_xml(self.G, f"{self.output_prefix}_{self.net_kind}.graphml.xml")
        if gexf:
            ctx.obj.INFO(f"Exporting {self.net_kind} as gexf file")
            nx.write_gexf(self.G, f"{self.output_prefix}_{self.net_kind}.gexf")


def get_command():
    _sys_argv = sys.argv
    for i in range(len(_sys_argv)):
        if _sys_argv[i] == "-i":
            _sys_argv[i+1] = os.path.abspath(_sys_argv[i+1])
    return "#command: DBRetina " + " ".join(_sys_argv[1:])


@cli.command(name="genenet", epilog = dbretina_doc.doc_url("genenet"), help_priority=9)
# @click.option('-t', '--threads', "user_threads", default=1, required=False, type=int, help="number of cores") # TODO later in C++ version
@click.option('-i', '--index-prefix', "index_prefix", required=True, type=click.STRING, help="index file prefix")
@click.option('-p', '--pairwise', 'pairwise_file', required=True, type=click.Path(exists=True), help="pairwise TSV file")
@click.option('--graphml', 'graphml', is_flag=True, default = False, help="export genenet as graphml file")
@click.option('--gexf', 'gexf', is_flag=True, default = False, help="export genenet as gexf file")
@click.option('-o', '--output-prefix', "output_prefix", required=True, type=click.STRING, help="output file prefix")
@click.pass_context
def main(ctx, index_prefix, pairwise_file, output_prefix, graphml, gexf):
    """Construct a Gene Network.
 

    Detailed description:


    For a groups pairwise file, construct an gene network between the features of each group and all other features in the pairwise file.
    """

    ###################################
    # 1. Extract gene set pairs
    ###################################
    
    ctx.obj.INFO(f"Extracting gene set pairs from {pairwise_file}")

    geneSet_pairs = set()
    metadata = []

    # Try Parquet/PairwiseStore first
    try:
        from dbretina.compat import open_pairwise
        store = open_pairwise(pairwise_file)
    except Exception:
        store = None

    if store is not None:
        ctx.obj.INFO("using Parquet pairwise data via PairwiseStore")
        metadata.append(f"#command: {get_command()}\n")
        names_map = store.get_names_map()
        df = store.to_pandas(columns=["group_1_id", "group_2_id"])
        for _, row in df.iterrows():
            name1 = names_map.get(int(row["group_1_id"]), "")
            name2 = names_map.get(int(row["group_2_id"]), "")
            geneSet_pairs.add((name1, name2))
        store.close()
    else:
        # No usable parquet store. Resolve the canonical .dbrp sibling from any
        # -p form (tsv / parquet dir / .dbrp); resolve_dbrp_path is isfile-guarded
        # so a directory never reaches the C++ .dbrp reader. genenet needs only the
        # group name pairs (no metric/cutoff filtering), so iterate every pair.
        # Mirrors graph (bipartite_graph.py) / bipartite (bipartite_pairwise.py);
        # without this the binary .dbrp was text-read and crashed with a raw
        # UnicodeDecodeError, and a storeless parquet directory crashed with a raw
        # IsADirectoryError (issues 066/067). Covers genenet AND interactome
        # (shared callback).
        from dbretina.compat import resolve_dbrp_path
        dbrp_path = resolve_dbrp_path(pairwise_file)
        if dbrp_path is not None:
            import _dbretina_internal as dbretina_internal
            ctx.obj.INFO(f"found .dbrp file: {dbrp_path}, using binary pairwise reader")
            metadata.append(f"#command: {get_command()}\n")
            for rec in dbretina_internal.dbrp_iterate_all(dbrp_path):
                geneSet_pairs.add((rec['group_1_name'], rec['group_2_name']))
        elif os.path.isdir(pairwise_file):
            # A directory with no usable parquet store (open_pairwise returned
            # None: no manifest.json) and no sibling .dbrp; the TSV open() below
            # would IsADirectoryError. Clean, actionable error instead.
            ctx.obj.ERROR(
                f"'{pairwise_file}' is a pairwise directory without a usable "
                f"parquet store (no manifest.json) and no sibling .dbrp; pass the "
                f"pairwise TSV, its parquet directory, or the .dbrp to -p.")
        else:
            # Fallback: existing TSV code
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

                header = next(pairwise_tsv)
                header = header.strip().split('\t')
                try:
                    group1_index = header.index('group_1_name')
                    group2_index = header.index('group_2_name')
                except:
                    group1_index = header.index('group_1')
                    group2_index = header.index('group_2')

                for row in pairwise_tsv:
                    row = row.strip().split('\t')
                    geneSet_pairs.add((row[group1_index], row[group2_index]))
            
    ##############################################
    # 2. Map gene set to genes (group to features)
    ##############################################
    ctx.obj.INFO(f"Mapping gene sets to features")

    dbri_path = f"{index_prefix}.dbri"
    # -i is a plain STRING (not click.Path); guard so a missing prefix gives a
    # clean [ERROR] instead of a raw FileNotFoundError from open(_raw.json)
    # (issue 029). Covers both 'genenet' and 'interactome' (shared callback).
    if not os.path.exists(dbri_path) and \
            not os.path.exists(f"{index_prefix}_raw.json"):
        ctx.obj.ERROR(f"index prefix '{index_prefix}' (.dbri / _raw.json) not found")
    if os.path.exists(dbri_path):
        import _dbretina_internal as dbretina_internal
        raw_json_str = dbretina_internal.dbri_load_raw_gene_sets(dbri_path)
        supergroups_to_features = json.loads(raw_json_str)["data"]
    else:
        raw_json_file = f"{index_prefix}_raw.json"
        with open(raw_json_file, "r") as f:
            supergroups_to_features = json.loads(f.read())["data"]
    
    ##############################################
    #3. Build the genenet
    ##############################################
    
    # 'genenet' and 'interactome' are the same callback; distinguish by the
    # invoked command name so outputs/log strings aren't mislabeled.
    net_kind = ctx.info_name if ctx.info_name in ("genenet", "interactome") else "genenet"

    ctx.obj.INFO(f"Building the {net_kind}. Please wait...")
    genenet = GeneNet(output_prefix, net_kind=net_kind)
    
    for geneSet_pair in tqdm(geneSet_pairs):
        geneSet1_features = supergroups_to_features[geneSet_pair[0]]
        geneSet2_features = supergroups_to_features[geneSet_pair[1]]
        
        # create edges between all pairs of features
        for feature1 in geneSet1_features:
            for feature2 in geneSet2_features:
                genenet.add_edge(feature1, feature2)
                
    genenet.export(ctx)
    genenet.graph_export(graphml, gexf, ctx)
    label = "Interactome" if net_kind == "interactome" else "Gene Network"
    ctx.obj.SUCCESS(f"{label} has been constructed successfully.")
    # genenet.plot_statistics()