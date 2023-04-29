#!/usr/bin/python
# -*- coding: utf-8 -*-
from __future__ import division

import os
import sys
import click
import pandas as pd
from kSpider2.click_context import cli
from scipy.cluster.hierarchy import linkage, to_tree, ClusterWarning
import numpy as np
from warnings import simplefilter
simplefilter("ignore", ClusterWarning)
from Bio import Phylo
from Bio.Phylo import BaseTree
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram
import plotly.figure_factory as ff
from scipy.spatial.distance import squareform
import seaborn as sns

# update recursive limit
sys.setrecursionlimit(15000)


# Recursive function to convert linkage tree to BioPython Tree
def convert_tree(tree):
    if tree.is_leaf():
        return BaseTree.Clade(branch_length=tree.dist, name=str(tree.id))
    else:
        left = convert_tree(tree.get_left())
        right = convert_tree(tree.get_right())
        return BaseTree.Clade(branch_length=tree.dist, clades=[left, right])


# Thanks to https://stackoverflow.com/a/31878514/3371177
def get_newick(node, parent_dist, leaf_names, newick='') -> str:
    """
    Convert sciply.cluster.hierarchy.to_tree()-output to Newick format.

    :param node: output of sciply.cluster.hierarchy.to_tree()
    :param parent_dist: output of sciply.cluster.hierarchy.to_tree().dist
    :param leaf_names: list of leaf names
    :param newick: leave empty, this variable is used in recursion.
    :returns: tree in Newick format
    """
    if node.is_leaf():
        return "%s:%.2f%s" % (leaf_names[node.id], parent_dist - node.dist, newick)
    else:
        if len(newick) > 0:
            newick = "):%.2f%s" % (parent_dist - node.dist, newick)
        else:
            newick = ");"
        newick = get_newick(node.get_left(), node.dist,
                            leaf_names, newick=newick)
        newick = get_newick(node.get_right(), node.dist,
                            leaf_names, newick=",%s" % (newick))
        newick = "(%s" % (newick)
        return newick


@cli.command(name="export", help_priority=5)
@click.option('-i', '--index-prefix', "index_prefix", required=True, type=click.STRING, help="Index file prefix")
@click.option('-p', '--pairwise', 'pairwise_file', required=False, type=click.Path(exists=True), help="filtered pairwise TSV file")
@click.option('--newick', "newick", is_flag=True, help="Convert pairwise (containment) matrix to newick format", default=False)
@click.option('-d', '--dist-type', "distance_type", required=False, default="max_cont", show_default=True, type=click.STRING, help="select from ['min_cont', 'avg_cont', 'max_cont', 'ochiai', 'jaccard']")
@click.option('-o', "overwritten_output", default="na", required=False, type=click.STRING, help="custom output file name prefix")
@click.pass_context
def main(ctx, index_prefix, pairwise_file, newick, distance_type, overwritten_output):
    """
    Export the pairwise TSV to dissimilarity matrix or newick format.
    """
    if(not pairwise_file):
        ctx.obj.INFO("No pairwise file provided, using the default one")
        pairwise_file = f"{index_prefix}_DBRetina_pairwise.tsv"

    LOGGER = ctx.obj

    distance_to_col = {
        "min_cont": 5,
        "avg_cont": 6,
        "max_cont": 7,
        "ochiai": 8,
        "jaccard": 9,
    }

    if distance_type not in distance_to_col:
        LOGGER.ERROR("unknown distance!")

    dist_col = distance_to_col[distance_type]    

    # Check for existing pairwise file
    for _file in [pairwise_file]:
        if not os.path.exists(_file):
            LOGGER.ERROR(f"File {_file} is not found.")

    """Parse DBRetina's pairwise
    """

    distances = {}    
    newick_out = pairwise_file.replace(".tsv", f"_{distance_type}.newick")
    distmatrix_out = pairwise_file.replace(".tsv", f"_{distance_type}_distmat.tsv")

    if overwritten_output != "na":
        distmatrix_out = f"{overwritten_output}_{distance_type}_distmat.tsv"
        newick_out = f"{overwritten_output}_{distance_type}.newick"

    with open(pairwise_file) as PAIRWISE:
        next(PAIRWISE) # Skip header
        for line in PAIRWISE:
            line = (line.strip().split('\t'))
            # add double quotes to group names
            grp1 = f"\"{line[2]}\""
            grp2 = f"\"{line[3]}\""
            dist_metric = float(line[dist_col])
            # convert xx% to 0.xx for the distance matrix
            distances[(grp1, grp2)] = dist_metric / 100


    elements = set()
    # for pair in distances.keys():
    #     elements.update(pair[:2])
    # index_map = {element: i for i, element in enumerate(sorted(elements))}
    # n = len(elements)
    # dist_matrix = np.zeros((n, n))
    # # Fill in the upper triangle of the distance matrix with the pairwise distances
    # for (src1, src2), dist in distances.items():
    #     i = index_map[src1]
    #     j = index_map[src2]
    #     dist_matrix[i, j] = dist_matrix[j, i] = dist
    # src_names = sorted(elements)
    # dist_df = pd.DataFrame(dist_matrix, index=src_names, columns=src_names)
    # dist_df.to_csv(distmatrix_out + ".new.tsv", sep='\t')


    unique_ids = sorted({x for y in distances for x in y})
    df = pd.DataFrame(index=unique_ids, columns=unique_ids)
    for k, v in distances.items():
        # 1-v for dissimilarity
        df.loc[k[0], k[1]] = 1-v
        df.loc[k[1], k[0]] = 1-v

    df = df.fillna(0)
    LOGGER.INFO(f"Writing distance matrix to {distmatrix_out}")
    df.to_csv(distmatrix_out, sep='\t')

    if newick:
        loaded_df = pd.read_csv(distmatrix_out, sep='\t')
        LOGGER.INFO(f"Writing newick to {newick_out}.")
        names = list(loaded_df.columns[1:])
        dist = loaded_df[loaded_df.columns[1:]].to_numpy()
        Z = linkage(dist, 'average')
        tree = to_tree(Z, False)

        newick = get_newick(tree, tree.dist, names)
        with open(newick_out, 'w') as NW:
            NW.write(newick)

    LOGGER.SUCCESS("Done.")
