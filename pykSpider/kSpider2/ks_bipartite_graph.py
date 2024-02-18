#!/usr/bin/python
# -*- coding: utf-8 -*-
from __future__ import division
import sys
import _kSpider_internal as kSpider_internal
import click
import contextlib
from kSpider2.click_context import cli
import subprocess
import os
import pandas as pd
import networkx as nx
import plotly.graph_objects as go
import plotly.express as px
import kSpider2.dbretina_doc_url as dbretina_doc
from kSpider2.ks_setcov import main as setcov_main
import networkx as nx
from collections import defaultdict

def execute_bash_command(command):
    try:
        process = subprocess.Popen(
            command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, text=True
        )
        output, error = process.communicate()

        if process.returncode == 0:
            return True
        print("Command execution failed.", file=sys.stderr)
        print(f"Error:\n{error}", file=sys.stderr)
        return False

    except Exception as e:
        print(f"An error occurred: {e}", file=sys.stderr)
        return False

def get_command():
    _sys_argv = sys.argv
    for i in range(len(_sys_argv)):
        if os.path.isfile(_sys_argv[i]):
            _sys_argv[i] = os.path.abspath(_sys_argv[i])
    return "DBRetina " + " ".join(_sys_argv[1:])


def path_to_absolute_path(ctx, param, value):
    with contextlib.suppress(Exception):
        return os.path.abspath(value) if value is not None else None

def check_if_there_is_a_pvalue(pairwise_file):
    with open(pairwise_file) as F:
        for line in F:
            if not line.startswith("#"):
                return "pvalue" in line
            else:
                continue

def validate_all_files_exist(ctx, param, value):
    if value is None:
        return None
    for path in value:
        if not os.path.exists(path):
            raise click.BadParameter(f"File '{path}' doesn't exist")
    return value


class DBRetinaGraph:
    
    # graph with type of nx.Graph()
    graph : nx.Graph = None
    node_attributes = None    
    pairwise_file = None
    index_prefix = None
    metric = None
    cutoff : float = None
    metadata = []
    dbretina_str_escape = lambda x: x.lower().replace('"', '')
    target_to_gene_sets = defaultdict(set)

    node_to_fragmentation = defaultdict(int)
    node_to_heterogeneity = defaultdict(int)
    node_to_modularity = defaultdict(int)
    node_to_target_name = {}
    target_groups = set()
    
    target_to_targetGroupID = {}
    geneSetToTargetsArgumentID = {}
    
    gene_set_to_targetID = {}
    
    pairwise_df = None

    
    metric_to_col = {
        "containment": 5,
        "ochiai": 6,
        "jaccard": 7,
        "odds_ratio": 8,
        "pvalue": 9,
    }
    
    metric_col : int = None
    
    
    def __init__(self, pairwise_file, index_prefix, metric, cutoff, LOGGER, inter_targets, intra_targets, output_prefix):
        self.metric_col = self.metric_to_col[metric]
        self.pairwise_file = pairwise_file
        self.index_prefix = index_prefix
        self.inter_targets = inter_targets
        self.intra_targets = intra_targets
        self.cutoff = cutoff
        self.metric = metric
        self.LOGGER = LOGGER
        self.load_all_targets()
        self.parse_node_size(index_prefix)
        self.output_prefix = output_prefix

        
    def load_all_pairwise(self):
        self.pairwise_df = pd.read_csv(self.pairwise_file, sep='\t', comment='#', usecols=['group_1_name','group_2_name', self.metric])

    def filter_by_cutoff(self):
        self.pairwise_df = self.pairwise_df[self.pairwise_df[self.metric] >= self.cutoff]

        
    def load_all_targets(self):
        self.LOGGER.INFO("Loading all targets")

        # lambda function to get basename without extension (might contain dots)
        get_basename = lambda x: os.path.splitext(os.path.basename(x))[0]

        self.LOGGER.INFO("Loading inter-targets")        
        for _targetdb_id, targets_set in enumerate(self.inter_targets, start=1):
            for target_file in targets_set:
                gene_sets = self.load_groups(target_file)
                file_basename = get_basename(target_file)
                self.target_to_targetGroupID[file_basename] = f"inter_{_targetdb_id}"
                self.target_to_gene_sets[file_basename] = gene_sets
                for gene_set_name in gene_sets:
                    self.gene_set_to_targetID[gene_set_name] = file_basename
                    self.geneSetToTargetsArgumentID[gene_set_name] = f"inter_{_targetdb_id}"


        self.LOGGER.INFO("Loading intra-targets")
        for _targetdb_id, targets_set in enumerate(self.intra_targets, start=1):
            for target_file in targets_set:
                gene_sets = self.load_groups(target_file)
                file_basename = get_basename(target_file)
                self.target_to_targetGroupID[file_basename] = f"intra_{_targetdb_id}"
                self.target_to_gene_sets[file_basename] = gene_sets
                for gene_set_name in gene_sets:
                    self.gene_set_to_targetID[gene_set_name] = file_basename
                    self.geneSetToTargetsArgumentID[gene_set_name] = f"intra_{_targetdb_id}"
        
        
        # now we have self.target_to_targetGroupID to map interTarget to group ID (group might contain one or more targets))
        # and we have self.target_to_gene_sets to map target to gene sets
        # and geneSetToTargetsArgumentID to map gene set to target ID (intra or inter) and the serial ID
        
        # make sure there is no overlap between values of self.target_to_gene_sets
        # report common groups between targets to a file
        # i.e. no gene set is in two targets
        common_groups = set()
        for target_name, groups in self.target_to_gene_sets.items():
            for other_target_file, other_groups in self.target_to_gene_sets.items():
                if target_name == other_target_file:
                    continue
                if len(groups.intersection(other_groups)) > 0:
                    common_groups.update(groups.intersection(other_groups))
                    with open(f"{self.output_prefix}_ERROR_common_groups.txt", 'w') as f:
                        for group in common_groups:
                            f.write(f"{group}\n")
                    self.LOGGER.ERROR(f"Target files {target_name} and {other_target_file} have common groups.\nCheck {self.output_prefix}_ERROR_common_groups.txt")


        self.LOGGER.INFO(f"Total number of gene sets: {len(self.geneSetToTargetsArgumentID)}")
        self.LOGGER.INFO(f"Total number of intra-targets groups: {len(self.intra_targets)}")
        self.LOGGER.INFO(f"Total number of inter-targets groups: {len(self.inter_targets)}")
        self.LOGGER.INFO(f"Total number of targets: {len(self.target_to_gene_sets)}")
        


    def load_groups(self, groups_file):
        groups = set()
        with open(groups_file) as F:
            for line in F:
                escaped_group = line.strip().split('\t')[0].lower().replace('"', '')
                groups.add(escaped_group)
            
        return groups
    
    def parse_node_size(self, index_prefix):
        self.node_to_size = {}
        with open(f"{index_prefix}_groupID_to_featureCount.tsv") as F:
            next(F)
            for line in F:
                node, size = line.strip().split('\t')
                self.node_to_size[node] = int(size)


    def export_node_attributes(self, include_isolates):
        df_nodes = pd.DataFrame.from_dict(self.geneSetToTargetsArgumentID, orient='index', columns=['targetGroup'])
        df_nodes["target_name"] = df_nodes.index.map(self.gene_set_to_targetID)
        df_nodes['heterogeneity'] = df_nodes.index.map(self.node_to_heterogeneity)
        df_nodes['geneSet_size'] = df_nodes.index.map(self.node_to_size)
        df_nodes['fragmentation'] = df_nodes.index.map(self.node_to_fragmentation)        
        df_nodes['modularity'] = abs(df_nodes['fragmentation'] + df_nodes['heterogeneity'])
        df_nodes.index.name = 'id'
        self.LOGGER.INFO(f"total number of nodes: {len(df_nodes)}")
        # keep only nodes in self.nodes_with_edges
        if not include_isolates:
            df_nodes = df_nodes.loc[list(self.nodes_with_edges)]
            self.LOGGER.INFO(f"remaining nodes after removing isolates: {len(df_nodes)}")
    
        df_nodes.to_csv(f"{self.output_prefix}_nodes.tsv", sep='\t')


    def pairwise_file_iterator(self, output_prefix):
        with open(self.pairwise_file) as pairwise_tsv:
            while True:
                pos = pairwise_tsv.tell()
                line = pairwise_tsv.readline()
                if line.startswith('#'):
                    continue
                pairwise_tsv.seek(pos)
                break
            
            next(pairwise_tsv)  # Skip header
            for row in pairwise_tsv:
                row = row.strip().split('\t')
                similarity = float(row[self.metric_col])

                # first skip similarity
                if similarity < self.cutoff:
                    continue

                node_1 = row[2]
                node_2 = row[3]

                yield node_1, node_2, similarity
                

    def build_graph(self):
        total_edges = 0
        pairwise_iter = self.pairwise_file_iterator(self.output_prefix)
        self.nodes_with_edges = set()
        with open(f"{self.output_prefix}_edges.tsv", 'w') as f_edges:
            f_edges.write(f"from\tto\t{self.metric}\n")
            for gene_set_1, gene_set_2, similarity in pairwise_iter:
                gene_set_1_targetArgumentID = self.geneSetToTargetsArgumentID[gene_set_1]
                gene_set_2_targetArgumentID = self.geneSetToTargetsArgumentID[gene_set_2]
                
                group1_len = self.node_to_size[gene_set_1]
                group2_len = self.node_to_size[gene_set_2]

                if group1_len < group2_len:
                    self.node_to_fragmentation[gene_set_1] -= 1
                    self.node_to_heterogeneity[gene_set_2] += 1
                elif group1_len > group2_len:
                    self.node_to_fragmentation[gene_set_2] -= 1
                    self.node_to_heterogeneity[gene_set_1] += 1

                both_same_intra = gene_set_1_targetArgumentID == gene_set_2_targetArgumentID and gene_set_1_targetArgumentID.startswith('intra') and gene_set_2_targetArgumentID.startswith('intra')
                coming_from_differnt_arguments = gene_set_1_targetArgumentID != gene_set_2_targetArgumentID

                if both_same_intra or coming_from_differnt_arguments:
                    f_edges.write(f"{gene_set_1}\t{gene_set_2}\t{similarity}\n")
                    self.nodes_with_edges.add(gene_set_1)
                    self.nodes_with_edges.add(gene_set_2)
                    total_edges += 1
                else:
                    continue
        
        self.LOGGER.INFO(f"Total number of edges: {total_edges}")
                

def process_targets_option(ctx, param, value):
    file_groups = []
    for file_group in value:
        file_list = [f.strip() for f in file_group.split(',')]
        for file in file_list:
            if not os.path.exists(file):
                raise click.BadParameter(f"File does not exist: {file}")
        # file_list = [os.path.abspath(file) for file in file_list] # for absolute paths
        file_groups.append(file_list)
    return file_groups

### -----> As a start, only create edges and include node size, and edge weights. <----- ###

@cli.command(name="graph", epilog = dbretina_doc.doc_url("graph"), help_priority=9)
@click.option('-i', '--index-prefix', 'index_prefix', required=True, type=click.STRING, help="Index file prefix")
@click.option('-p', '--pairwise', 'pairwise_file', callback=path_to_absolute_path, required=True, type=click.Path(exists=True), help="the pairwise TSV file")
@click.option('--intra-targets', "intra_targets", multiple=True, callback=process_targets_option, required=False, help="comma separated list of TSV files with first column as gene sets")
@click.option('--inter-targets', "inter_targets", multiple=True, callback=process_targets_option, required=False, help="comma separated list of TSV files with first column as gene sets")
@click.option('-m', '--metric', "metric", required=True, type=click.STRING, help="Similarity metric ['containment', 'ochiai', 'jaccard', 'pvalue']")
@click.option('-c', '--cutoff', 'cutoff', required=False, type=click.FloatRange(0, 100, clamp=False), default=0.0, show_default = True, help="Include comparisons (similarity > cutoff)")
@click.option('-o', '--output', "output_prefix", required=True, type=click.STRING, help="output file prefix")
@click.option('--include-isolates', "include_isolates", is_flag=True, default=False, show_default = True, help="Include isolate nodes")
@click.option('--visualize', "visualize", is_flag=True, default=False, show_default = True, help="Visualize the graph")
@click.pass_context
def main(ctx, index_prefix, pairwise_file, intra_targets, inter_targets, metric, cutoff, output_prefix, include_isolates, visualize):
    """
        Export edges, nodes graph files for visualization.
        Optionally visualize the DBRetina's graph.
    """
    LOGGER = ctx.obj

    LOGGER.INFO(
        f"""Running DBRetina graph with parameters: \n\t\t 
        pairwise_file: {pairwise_file} \n\t\t 
        intra-targets: {intra_targets} \n\t\t 
        inter-targets: {inter_targets} \n\t\t 
        metric: {metric} \n\t\t 
        cutoff: {cutoff} \n\t\t 
        output_prefix: {output_prefix}
        """
    )
    
    
    db_graph = DBRetinaGraph(
        index_prefix=index_prefix,
        cutoff=cutoff,
        intra_targets=intra_targets,
        inter_targets=inter_targets,
        LOGGER=LOGGER,
        metric=metric,
        output_prefix=output_prefix,
        pairwise_file=pairwise_file
    )
    
    db_graph.build_graph()
    db_graph.export_node_attributes(include_isolates)
    
    if visualize:
        from kSpider2.dbretina_viz import DBRetinaViz

        edge_df = pd.read_csv(f"{output_prefix}_edges.tsv", sep='\t')
        node_df = pd.read_csv(f"{output_prefix}_nodes.tsv", sep='\t')
        viz = DBRetinaViz(edge_df, node_df)
        viz.plot(debug=False)
    
    LOGGER.SUCCESS("Done!")
    

