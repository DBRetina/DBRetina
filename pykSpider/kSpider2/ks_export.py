#!/usr/bin/python
# -*- coding: utf-8 -*-
from __future__ import division

import os
import click
import pandas as pd
from kSpider2.click_context import cli
from scipy.cluster.hierarchy import linkage, to_tree, ClusterWarning
from warnings import simplefilter
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from scipy.cluster.hierarchy import linkage, to_tree
from scipy.spatial.distance import squareform, pdist
simplefilter("ignore", ClusterWarning)
import re
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
import plotly.io as pio
import dash_bio
from pycirclize import Circos
from io import StringIO
from Bio import Phylo
import kSpider2.dbretina_doc_url as dbretina_doc

def newick_str_escape(name):
    # Remove special characters
    name = re.sub(r'[:;(),\[\]\-\' ]', '_', name)
    return name

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

def generate_heatmap(df):
    # Ensure DataFrame is numeric
    # Set up figure and axes
    fig, ax = plt.subplots(figsize=(10, 8))

    # Generate a mask for the upper triangle
    mask = np.triu(np.ones_like(df, dtype=bool))

    # Generate a custom diverging colormap
    cmap = sns.diverging_palette(230, 20, as_cmap=True)

    # Draw the heatmap
    sns.heatmap(df, mask=mask, cmap=cmap, center=0, square=True, linewidths=.5, cbar_kws={"shrink": .5}, ax=ax, annot=False, fmt=".2f")

    # Set the title
    ax.set_title('Heatmap')

    # Rotate the x-axis labels
    # ax.set_xticklabels(ax.get_xticklabels(), rotation=45, horizontalalignment='right')
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, horizontalalignment='right')
    
    # Rotate the y-axis labels
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0)

    # Show the plot
    plt.tight_layout()
    plt.savefig("heatmap.png", dpi=500)


# click callback to make sure it's selected from ['ids', 'names']
def validate_labels(ctx, param, value):
    if value not in ['ids', 'names']:
        raise click.BadParameter('labels must be one of ids, names')
    return value

def dataframe_to_newick(df):
    np.fill_diagonal(df.values, 0)
    condensed_matrix = squareform(df)
    Z = linkage(condensed_matrix, 'single')
    tree = to_tree(Z, rd=False)
    return tree_to_newick(tree, list(df.columns))

def tree_to_newick(node, leaf_names):
    if node.is_leaf():
        return leaf_names[node.id]
    left_child = tree_to_newick(node.get_left(), leaf_names)
    right_child = tree_to_newick(node.get_right(), leaf_names)
    return f"({left_child}:{node.dist:.2f},{right_child}:{node.dist:.2f})"


def similarity_df_to_newick(similarity_df, method):
    # Convert similarity matrix to distance matrix
    dissimilarity_df = 100 - similarity_df

    # Ensure the diagonal is zero
    np.fill_diagonal(dissimilarity_df.values, 0)

    # Convert the distance matrix DataFrame to a condensed distance matrix for linkage
    # condensed_distance_matrix = squareform(dissimilarity_df)
    condensed_distance_matrix = pdist(dissimilarity_df, metric = 'euclidean')

    # Perform hierarchical/agglomerative clustering
    Z = linkage(condensed_distance_matrix, method)

    # Convert linkage matrix to a tree
    tree = to_tree(Z, rd=False)

    return tree_to_newick(tree, list(similarity_df.columns))

def tree_to_newick(node, leaf_names):
    if node.is_leaf():
        return leaf_names[node.id]
    left_child = tree_to_newick(node.get_left(), leaf_names)
    right_child = tree_to_newick(node.get_right(), leaf_names)
    return f"({left_child}:{node.dist:.2f},{right_child}:{node.dist:.2f})"


def check_if_there_is_a_pvalue(pairwise_file):
    with open(pairwise_file) as F:
        for line in F:
            if not line.startswith("#"):
                return "pvalue" in line
            else:
                continue

def export_heatmap(df, filename):
    plt.figure(figsize=(10,8)) # Adjust size as needed
    
    # Create a heatmap with seaborn
    cmap = sns.diverging_palette(220, 10, as_cmap=True) 
    heatmap = sns.heatmap(df, cmap=cmap) 

    # Adding padding to prevent cutting of the heatmap
    b, t = plt.ylim()
    b += 0.5
    t -= 0.5
    plt.ylim(b, t)

    # Customize the heatmap
    heatmap.set_title('Similarity Heatmap', fontdict={'fontsize':18}, pad=16)
    # plt.xlabel('Column Name', fontsize=14) 
    # plt.ylabel('Column Name', fontsize=14) 
    # legend title
    # heatmap.legend(title="Similarity", bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.yticks(rotation=0)

    # Save the heatmap to a file
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()  # Close the figure after saving it to a file

@cli.command(name="export", epilog = dbretina_doc.doc_url("export"), help_priority=5)
@click.option('-p', '--pairwise', 'pairwise_file', required=True, type=click.Path(exists=True), help="pairwise TSV file")
@click.option('-m', '--metric', "metric", required=True, type=click.STRING, help="select from ['containment', 'ochiai', 'jaccard', 'pvalue']")
@click.option('--newick', "newick", is_flag=True, help="Convert the distance matrix to newick tree format", default=False)
@click.option('-l', '--labels', "labels_selection", callback = validate_labels, required=False, default="names", show_default=True, type=click.STRING, help="select from ['ids', 'names']")
@click.option('--linkage', "linkage_method", required=False, default="ward", show_default=True, type=click.STRING, help="select from ['single', 'complete', 'average', 'weighted', 'centroid', 'median', 'ward']")
@click.option('-o', "--output", "output_prefix", required=True, type=click.STRING, help="output prefix")
@click.pass_context
def main(ctx, pairwise_file, newick, metric, output_prefix, labels_selection, linkage_method):
    """Export to distance matrix and tree formats.
    
    Export a pairwise TSV file into a distance matrix, newick-format file and circular dendrogram.
    
    """
 
    LOGGER = ctx.obj

    metric_to_col = {
        "containment": 5,
        "ochiai": 6,
        "jaccard": 7,
        "odds_ratio": 8,
        "pvalue": 9,
    }

    if metric not in metric_to_col:
        LOGGER.ERROR("unknown metric!")
        
    
    # check if pvalue
    if metric == "pvalue" and not check_if_there_is_a_pvalue(pairwise_file):
        LOGGER.ERROR("pvalue not found in pairwise file!")

    dist_col = metric_to_col[metric]    

    # Check for existing pairwise file
    for _file in [pairwise_file]:
        if not os.path.exists(_file):
            LOGGER.ERROR(f"File {_file} is not found.")


    # Parse DBRetina's pairwise
    
    if labels_selection == "ids": src1_label_col = 0; src2_label_col = 1
    else: src1_label_col = 2; src2_label_col = 3
    
    df = pd.read_csv(pairwise_file, sep='\t', usecols=[src1_label_col, src2_label_col, dist_col], comment='#')
    
    if labels_selection == "ids":    
        df[df.columns[0]] = df[df.columns[0]].astype(str)
        df[df.columns[1]] = df[df.columns[1]].astype(str)
    else: # escape newick_str_escape
        LOGGER.WARNING("Escaping newick characters in labels if any.")
        df[df.columns[0]] = df[df.columns[0]].apply(lambda x: newick_str_escape(x))
        df[df.columns[1]] = df[df.columns[1]].apply(lambda x: newick_str_escape(x))
    
    # df[df.columns[2]] = df[df.columns[2]].apply(lambda x: math.log2(x) if x != 0 else 0)
    
    similarity_df = df.pivot(index=df.columns[0], columns=df.columns[1], values=df.columns[2])

    similarity_df = similarity_df.combine_first(similarity_df.T).fillna(0)
    # np.fill_diagonal(similarity_df.values, math.log2(100))
    np.fill_diagonal(similarity_df.values, 100)
    
    # convert to distance matrix
    distance_matrix_df = 100 - similarity_df
    
    export_heatmap(similarity_df, f"{output_prefix}_heatmap.png")
    
    LOGGER.INFO(f"Writing heatmap to {output_prefix}_heatmap.html")
    # pio.write_html(fig, f"{output_prefix}_heatmap.html")
    # fig.write_image(f"{output_prefix}_heatmap.png", width=2700, height=2800, scale=1)

    ##################### PLOTLY END #####################
    
    
    # serialize distance matrix to binary format
    LOGGER.INFO(f"serializing the distance matrix to {output_prefix}_distmat.pkl")
    similarity_df.to_pickle(f"{output_prefix}_distmat.pkl")
    LOGGER.INFO(f"Writing distance matrix to {output_prefix}_distmat.tsv")
    similarity_df.to_csv(f"{output_prefix}_distmat.tsv", sep='\t')
    newick_out = f"{output_prefix}.newick"

    if newick:
        # Call the function with your similarity DataFrame and 'single' linkage method
        try:
            newick_string = similarity_df_to_newick(similarity_df, linkage_method)
            with open(newick_out, 'w') as NW:
                    NW.write(newick_string)
            
            # TODO: add a flag to visualize the tree
            # Beta: visualize the newick tree
            tree = Phylo.read(StringIO(newick_string), "newick")
            # Initialize circos sector with tree size
            circos = Circos(sectors={"Tree": tree.count_terminals()})
            sector = circos.sectors[0]
            track = sector.add_track((30, 100))
            track.tree(tree, leaf_label_size=6)
            LOGGER.INFO(f"Writing dendrogram tree to {output_prefix}_dendrogram.png")
            fig = circos.savefig(f"{output_prefix}_dendrogram.png", dpi=600)
            
        except RecursionError as err:
            LOGGER.ERROR(f"Couldn't handle the tree depth | {err}")

    LOGGER.SUCCESS("Done.")
