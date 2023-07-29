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
    target_to_groups = defaultdict(set)

    node_to_fragmentation = defaultdict(int)
    node_to_heterogeneity = defaultdict(int)
    node_to_modularity = defaultdict(int)
    node_to_target_name = {}
    target_groups = set()

    
    metric_to_col = {
        "containment": 5,
        "ochiai": 6,
        "jaccard": 7,
        "odds_ratio": 8,
        "pvalue": 9,
    }
    
    metric_col : int = None
    
    
    def __init__(self, pairwise_file, index_prefix, metric, cutoff, LOGGER, all_targets, output_prefix):
        self.metric_col = self.metric_to_col[metric]
        self.pairwise_file = pairwise_file
        self.index_prefix = index_prefix
        self.cutoff = cutoff
        self.metric = metric
        self.LOGGER = LOGGER
        self.load_all_targets(all_targets)
        self.parse_node_size(index_prefix)
        self.output_prefix = output_prefix
        
    def load_all_targets(self, all_targets):
        self.LOGGER.INFO("Loading all targets")

        # lambda function to get basename without extension (might contain dots)
        get_basename = lambda x: os.path.splitext(os.path.basename(x))[0]

        self.target_to_groups = defaultdict(set)
        for target_file in all_targets:
            groups = self.load_groups(target_file)
            target_name = get_basename(target_file)
            self.target_to_groups[target_name] = groups

        for target_name, groups in self.target_to_groups.items():
            for other_target_file, other_groups in self.target_to_groups.items():
                if target_name == other_target_file:
                    continue
                if len(groups.intersection(other_groups)) > 0:
                    print(f"common groups between {target_name} and {other_target_file}: {groups.intersection(other_groups)}")
                    self.LOGGER.ERROR(f"Target files {target_file} and {other_target_file} have common groups")

        # construct node_to_target_name
        for target_name, groups in self.target_to_groups.items():
            for group in groups:
                self.node_to_target_name[group] = target_name

        self.target_groups = set(self.node_to_target_name.keys())
        print(f"size of target groups: {len(self.target_groups)}")

    def load_query_file(self, query_file):
        query_groups = self.load_groups(query_file)
        for target_file, groups in self.target_to_groups.items():
            if len(groups.intersection(query_groups)) > 0:
                self.LOGGER.ERROR(f"Query file {query_file} and target file {target_file} have common groups")
        return query_groups
        
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
        
    def export_node_attributes(self, output_prefix):
        df_nodes = pd.DataFrame.from_dict(self.node_to_fragmentation, orient='index', columns=['fragmentation'])
        df_nodes["target_name"] = df_nodes.index.map(self.node_to_target_name)
        df_nodes['heterogeneity'] = df_nodes.index.map(self.node_to_heterogeneity)
        df_nodes['size'] = df_nodes.index.map(self.node_to_size)
        df_nodes['modularity'] = abs(df_nodes['fragmentation'] + df_nodes['heterogeneity'])
        df_nodes.index.name = 'id'
        df_nodes.to_csv(f"{output_prefix}_nodes.tsv", sep='\t')
    
    
    def build_graph():
        pass
     
    
class DBRetinaGraphWithQueryNoInternal(DBRetinaGraph):
    
    def build_graph(self, query_file):
  
        
        self.query_groups = self.load_query_file(query_file)

        with open(f"{self.output_prefix}_edges.tsv", 'w') as edges_file:
            edges_file.write(f"from\tto\t{self.metric}\n")

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
                    

                    # first skip similarity
                    if similarity < self.cutoff:
                        continue

                    node_1 = row[2]
                    node_2 = row[3]

                    group1_len = self.node_to_size[node_1]
                    group2_len = self.node_to_size[node_2]

                    if group1_len < group2_len:
                        self.node_to_fragmentation[node_1] -= 1
                        self.node_to_heterogeneity[node_2] += 1
                    elif group1_len > group2_len:
                        self.node_to_fragmentation[node_2] -= 1
                        self.node_to_heterogeneity[node_1] += 1

                    # second make sure one of source or target is in query groups and the other is in target groups
                                        
                    if node_1 in self.query_groups and node_2 in self.target_groups:
                        edges_file.write(f"{node_1}\t{node_2}\t{similarity}\n")
                    elif node_2 in self.query_groups and node_1 in self.target_groups:
                        edges_file.write(f"{node_2}\t{node_1}\t{similarity}\n")
                    else:
                        continue

            # export node attributes
            self.export_node_attributes(self.output_prefix)

class DBRetinaGraphWithQueryInternal(DBRetinaGraph):
    def build_graph(self, query_file = None):
        pass

class DBRetinaGraphTargetsOnly(DBRetinaGraph):
    def build_graph(self, query_file = None):
        pass
    
class DBRetinaGraphTargetsOnlyInternal(DBRetinaGraph):
    def build_graph(self, query_file = None):
        pass
    
# This should take multiple group/gmt files. A query group is optional.
# If query group is provided, then connections are made from query group to all other groups.
# if query group is not provided, then connections are made between all groups.
# if query group is provided, and groups internal edges flag is set, then connections are made between all groups except query group.
# we can use the setcov command get some node attributes (size, fragmentation, heterogeinety, modularity, status (can be representative or not))

### -----> As a start, only create edges and include node size, and edge weights. <----- ###

@cli.command(name="graph", epilog = dbretina_doc.doc_url("graph"), help_priority=9)
@click.option('-i', '--index-prefix', 'index_prefix', required=True, type=click.STRING, help="Index file prefix")
@click.option('-p', '--pairwise', 'pairwise_file', callback=path_to_absolute_path, required=True, type=click.Path(exists=True), help="the pairwise TSV file")
@click.option('-q', '--query-file', "query_file", callback=path_to_absolute_path, required=False, type=click.Path(exists=True), help="TSV file with first column as gene sets")
@click.option('--target', "targets", multiple=True, callback=validate_all_files_exist, required=True, help="multiple TSV files with first column as gene sets")
@click.option('-m', '--metric', "metric", required=True, type=click.STRING, help="Similarity metric ['containment', 'ochiai', 'jaccard', 'pvalue']")
@click.option('-c', '--cutoff', 'cutoff', required=False, type=click.FloatRange(0, 100, clamp=False), default=0.0, show_default = True, help="Include comparisons (similarity > cutoff)")
@click.option('--internal', "internal_edges", is_flag=True, default=False, help="include edges within same target gene sets")
@click.option('-o', '--output', "output_prefix", required=True, type=click.STRING, help="output file prefix")
@click.pass_context
def main(ctx, index_prefix, pairwise_file, query_file, targets, metric, cutoff, internal_edges, output_prefix):
    """
        Visualize.
    """
    LOGGER = ctx.obj
    
    LOGGER.INFO(
        f"""Running DBRetina graph with parameters: \n\t\t 
        pairwise_file: {pairwise_file} \n\t\t 
        query_file: {query_file} \n\t\t 
        targets: {targets} \n\t\t 
        metric: {metric} \n\t\t 
        cutoff: {cutoff} \n\t\t 
        internal_edges: {internal_edges} \n\t\t 
        output_prefix: {output_prefix}
        """
    )

    # Modes:
    # 1. query, with target files | mode_name = 'query_targets'
    # 2. query, with target files and internal edges | mode_name = 'query_targets_internal'
    # 3. target files only | mode_name = 'targets'
    # 4. target files only and internal edges | mode_name = 'targets_internal'
    # 5. query only (ERROR) | mode_name = 'query_only'
    
    if query_file and not targets:
        LOGGER.ERROR("Query file provided without target files. Please provide target files.")
    elif query_file:
        mode_name = 'query_targets_internal' if internal_edges else 'query_targets'
    elif targets:
        mode_name = 'targets_internal' if internal_edges else 'targets'
    else:
        LOGGER.ERROR("Please provide target files with or without query.")
        
        
    
    
    if mode_name == 'query_targets':
        graph = DBRetinaGraphWithQueryNoInternal(pairwise_file, index_prefix, metric, cutoff, LOGGER, targets, output_prefix)
        graph.build_graph(query_file)
        
        
    

