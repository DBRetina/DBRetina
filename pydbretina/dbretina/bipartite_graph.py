#!/usr/bin/python
# -*- coding: utf-8 -*-
from __future__ import division
import sys
import _dbretina_internal as dbretina_internal
import click
import contextlib
from dbretina.click_context import cli
import subprocess
import os
import pandas as pd
import networkx as nx
import plotly.graph_objects as go
import plotly.express as px
import dbretina.dbretina_doc_url as dbretina_doc
from dbretina.setcov import main as setcov_main
from dbretina.validators import validate_metric
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

def validate_all_files_exist(ctx, param, value):
    if value is None:
        return None
    for path in value:
        if not os.path.exists(path):
            raise click.BadParameter(f"File '{path}' doesn't exist")
    return value


class DBRetinaGraph:

    # NOTE: every per-run *mutable* attribute (dicts/sets/lists) is initialized as an
    # INSTANCE attribute in __init__, NOT here. Declaring them at class scope made them
    # shared across instances, so a second DBRetinaGraph() in the same process inherited
    # the first run's entries (issue 045). Only read-only constants live at class scope.
    dbretina_str_escape = lambda x: x.lower().replace('"', '')

    metric_to_col = {
        "containment": 5,
        "ochiai": 6,
        "jaccard": 7,
        "csi": 8,
        "dice": 9,
        "odds_ratio": 10,
        "pvalue": 11,
    }


    def __init__(self, pairwise_file, index_prefix, metric, cutoff, LOGGER, inter_targets, intra_targets, output_prefix):
        # Per-run mutable state, instance-scoped so concurrent/repeated in-process
        # instantiations don't share or leak entries (issue 045).
        self.graph = None
        self.node_attributes = None
        self.metadata = []
        self.pairwise_df = None
        self.metric_col = None
        self.target_to_gene_sets = defaultdict(set)
        self.node_to_fragmentation = defaultdict(int)
        self.node_to_heterogeneity = defaultdict(int)
        self.node_to_modularity = defaultdict(int)
        self.node_to_target_name = {}
        self.target_groups = set()
        self.target_to_targetGroupID = {}
        self.geneSetToTargetsArgumentID = {}
        self.gene_set_to_targetID = {}
        # Guard the lookup: the -m callback validates real metrics but lets the "NA"
        # sentinel through, and a bare metric_to_col[...] would raise a raw KeyError.
        if metric not in self.metric_to_col:
            LOGGER.ERROR(f"Invalid metric '{metric}'. Choose from: {', '.join(self.metric_to_col)}")
            sys.exit(1)
        # Reject -m pvalue on a dataset that lacks it BEFORE iterating, on every
        # input form. The store/.dbrp routes otherwise crashed uncaught: the
        # parquet route raised a raw ValueError ("Unknown metric 'pvalue'") from
        # PairwiseStore._validate_metric and the .dbrp route silently yielded 0
        # edges (issue 043). pairwise_has_pvalue is format-aware (tsv/dir/.dbrp).
        from dbretina.compat import pairwise_has_pvalue
        if metric == "pvalue" and not pairwise_has_pvalue(pairwise_file):
            LOGGER.ERROR("pvalue not found in pairwise file!")
        self.metric_col = self.metric_to_col[metric]
        self.pairwise_file = pairwise_file
        self.index_prefix = index_prefix
        self.inter_targets = inter_targets
        self.intra_targets = intra_targets
        self.cutoff = cutoff
        self.metric = metric
        self.LOGGER = LOGGER
        # Set before load_all_targets(): its common-groups error path writes to
        # f"{self.output_prefix}_ERROR_common_groups.txt", so output_prefix must
        # already exist (it previously was assigned only after, raising a masked
        # AttributeError whenever that branch fired).
        self.output_prefix = output_prefix
        self.load_all_targets()
        self.parse_node_size(index_prefix)

        
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
        dbri_path = f"{index_prefix}.dbri"
        if os.path.exists(dbri_path):
            import _dbretina_internal as dbretina_internal
            self.node_to_size = dbretina_internal.dbri_load_group_feature_counts(dbri_path)
        else:
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
        metric_name_to_id = {
            "containment": 0, "ochiai": 1, "jaccard": 2, "csi": 3,
            "dice": 4, "odds_ratio": 5, "pvalue": 6,
        }

        # 1. Prefer Parquet/PairwiseStore (handles the parquet pairwise DIRECTORY
        #    form and a .tsv with a sibling Parquet dir), matching how genenet /
        #    bipartite / query route their pairwise input.
        try:
            from dbretina.compat import open_pairwise
            store = open_pairwise(self.pairwise_file)
        except Exception:
            store = None

        if store is not None:
            self.LOGGER.INFO("using Parquet pairwise data via PairwiseStore")
            names_map = store.get_names_map()
            try:
                reader = store.filter_pairs(
                    self.metric, self.cutoff,
                    columns=["group_1_id", "group_2_id", self.metric],
                )
                for batch in reader:
                    d = batch.to_pydict()
                    ids1 = d["group_1_id"]
                    ids2 = d["group_2_id"]
                    sims = d[self.metric]
                    for gid1, gid2, sim in zip(ids1, ids2, sims):
                        yield names_map[gid1], names_map[gid2], float(sim)
            finally:
                store.close()
            return

        # 2. .dbrp fast-path. Resolve the canonical .dbrp sibling from any -p
        #    form (tsv / parquet dir / .dbrp); resolve_dbrp_path is isfile-guarded
        #    so a directory never reaches the C++ .dbrp reader, and it finds the
        #    sibling even when -p is the parquet-directory form (issue 051).
        from dbretina.compat import resolve_dbrp_path
        dbrp_path = resolve_dbrp_path(self.pairwise_file)
        if dbrp_path is not None and self.metric in metric_name_to_id:
            mid = metric_name_to_id[self.metric]
            records = dbretina_internal.dbrp_filter_pairs(dbrp_path, mid, self.cutoff)
            for rec in records:
                yield rec['group_1_name'], rec['group_2_name'], float(rec[self.metric])
        elif os.path.isdir(self.pairwise_file):
            # No usable parquet store (no manifest.json) and no sibling .dbrp:
            # the TSV open() below would IsADirectoryError. Clean error instead.
            self.LOGGER.ERROR(
                f"'{self.pairwise_file}' is a pairwise directory without a usable "
                f"parquet store (no manifest.json) and no sibling .dbrp; pass the "
                f"pairwise TSV, its parquet directory, or the .dbrp to -p.")
        else:
            # 3. TSV fallback.
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
                # Groups absent from every target file (no targets at all -> the
                # documented default, or a partial target list) are "ungrouped".
                # A bare dict subscript here would raise a raw KeyError on the first
                # such edge; default them instead so the run stays graceful.
                gene_set_1_targetArgumentID = self.geneSetToTargetsArgumentID.get(gene_set_1, 'ungrouped')
                gene_set_2_targetArgumentID = self.geneSetToTargetsArgumentID.get(gene_set_2, 'ungrouped')
                # Record ungrouped endpoints so they still appear as nodes.
                self.geneSetToTargetsArgumentID.setdefault(gene_set_1, gene_set_1_targetArgumentID)
                self.geneSetToTargetsArgumentID.setdefault(gene_set_2, gene_set_2_targetArgumentID)

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
                # Ungrouped nodes (no targets / partial targets) are not subject to
                # the inter-group suppression: keep any edge that touches one so the
                # no-targets default yields a full graph instead of an empty one.
                touches_ungrouped = 'ungrouped' in (gene_set_1_targetArgumentID, gene_set_2_targetArgumentID)

                if both_same_intra or coming_from_differnt_arguments or touches_ungrouped:
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
@click.option('-m', '--metric', "metric", required=True, type=click.STRING, callback=validate_metric, help="Similarity metric ['containment', 'ochiai', 'jaccard', 'csi', 'dice', 'odds_ratio', 'pvalue']")
@click.option('-c', '--cutoff', 'cutoff', required=False, type=click.FloatRange(0, 100, clamp=False), default=0.0, show_default = True, help="Include comparisons (similarity >= cutoff)")
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
        # dash/visdcc are optional, undeclared deps used only by --visualize.
        # Guard the import so a missing dep is a clean actionable error, not a
        # raw ModuleNotFoundError traceback. The edges/nodes TSVs are already
        # written above, so they are unaffected.
        try:
            from dbretina.dbretina_viz import DBRetinaViz
        except ImportError as e:
            LOGGER.ERROR(
                f"--visualize requires the optional 'dash' and 'visdcc' packages "
                f"(missing: {e.name}). Install them with: "
                f"pip install dash visdcc dash-bootstrap-components"
            )
            sys.exit(1)

        edge_df = pd.read_csv(f"{output_prefix}_edges.tsv", sep='\t')
        node_df = pd.read_csv(f"{output_prefix}_nodes.tsv", sep='\t')
        viz = DBRetinaViz(edge_df, node_df)
        viz.plot(debug=False)
    
    LOGGER.SUCCESS("Done!")
    

