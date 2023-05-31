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
import math
import re
import plotly.graph_objects as go
import plotly.figure_factory as ff
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
import random
import string


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


def similarity_df_to_newick(similarity_df, method='average'):
    # Convert similarity matrix to dissimilarity matrix
    dissimilarity_df = 100 - similarity_df

    # Ensure the diagonal is zero
    np.fill_diagonal(dissimilarity_df.values, 0)

    # Convert the dissimilarity matrix DataFrame to a condensed distance matrix for linkage
    # condensed_distance_matrix = squareform(dissimilarity_df)
    condensed_distance_matrix = pdist(dissimilarity_df, 'euclidean')

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


@cli.command(name="export", help_priority=5)
@click.option('-p', '--pairwise', 'pairwise_file', required=True, type=click.Path(exists=True), help="filtered pairwise TSV file")
@click.option('-d', '--dist-type', "distance_type", required=False, default="max_cont", show_default=True, type=click.STRING, help="select from ['min_cont', 'avg_cont', 'max_cont', 'ochiai', 'jaccard']")
@click.option('--newick', "newick", is_flag=True, help="Convert the similarity matrix to newick tree format", default=False)
@click.option('-l', '--labels', "labels_selection", callback = validate_labels, required=False, default="ids", show_default=True, type=click.STRING, help="select from ['ids', 'names']")
@click.option('-o', "output_prefix", required=True, type=click.STRING, help="output prefix")
@click.pass_context
def main(ctx, pairwise_file, newick, distance_type, output_prefix, labels_selection):
    """Export to dissimilarity matrix and newick format.
    
    Export a pairwise TSV file to a dissimilarity matrix and (optionally) a newick-format file.
    
    """
 
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


    # Parse DBRetina's pairwise
    
    if labels_selection == "ids": src1_label_col = 0; src2_label_col = 1
    else: src1_label_col = 2; src2_label_col = 3
    
    df = pd.read_csv(pairwise_file, sep='\t', usecols=[src1_label_col, src2_label_col, dist_col], comment='#')
    
    if labels_selection == "ids":    
        df[df.columns[0]] = df[df.columns[0]].astype(str)
        df[df.columns[1]] = df[df.columns[1]].astype(str)
    else: # escape newick_str_escape
        df[df.columns[0]] = df[df.columns[0]].apply(lambda x: newick_str_escape(x))
        df[df.columns[1]] = df[df.columns[1]].apply(lambda x: newick_str_escape(x))
    
    # df[df.columns[2]] = df[df.columns[2]].apply(lambda x: math.log2(x) if x != 0 else 0)
    
    similarity_df = df.pivot(index=df.columns[0], columns=df.columns[1], values=df.columns[2])

    similarity_df = similarity_df.combine_first(similarity_df.T).fillna(0)
    # np.fill_diagonal(similarity_df.values, math.log2(100))
    np.fill_diagonal(similarity_df.values, 100)
    

    cmap = sns.color_palette("Spectral", as_cmap=True)
    g = sns.clustermap(
        similarity_df, 
        cmap=cmap, 
        center=0, 
        linewidths=.5, 
        figsize=(10, 10),
        row_cluster=True, 
        col_cluster=True, 
        vmin=0, 
        # vmax=math.log2(100),
        vmax=100,
        dendrogram_ratio=(0.1, 0.2),
        # method = 'single',
        # metric='braycurtis',
        )
    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=45, horizontalalignment='right')
    plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.cax.set_title('Similarity', loc='left', fontsize=12, fontweight='bold')
    
    
    
    ### PLOTLY
    
    import plotly.graph_objects as go
    import plotly.io as pio
    import dash_bio
    
    ##### DASH
    
    fig = dash_bio.Clustergram(
        data=similarity_df,
        column_labels=list(similarity_df.columns.values),
        row_labels=list(similarity_df.index),
        height=2800,
        width=2700,
        
        # display_ratio=[0.1, 0.7]
    )
    
    pio.write_html(fig, 'dash-heatmap.html')
    
        # Your precomputed data
    z_data = similarity_df.values

    fig = go.Figure(data=go.Heatmap(
        z=z_data,
        x=similarity_df.columns,
        y=similarity_df.index,
        colorscale='YlGnBu',  # Adjust to a more distinguishable color scale
        hoverongaps = False))

    fig.update_layout(
        title='Similarity Heatmap',
        autosize=True,
        width=2000,  # Adjust to control block size
        height=2000,  # Adjust to control block size
        xaxis_nticks=len(similarity_df.columns),  # try to show all x labels
        # yaxis_nticks=len(similarity_df.index)  # try to show all y labels
    )
    pio.write_html(fig, 'heatmap.html')
    # fig.show()
    
    #############################
    
    # PLOTLY 2
    from scipy.spatial import distance
    from scipy.cluster import hierarchy
    
    # Preprocessing the data
    # Check for any missing data
    if similarity_df.isnull().any().any():
        print('Missing data found. Filling with zeros.')
        similarity_df.fillna(0, inplace=True)

    # Preprocessing the data
    labels = similarity_df.columns
    data_array = similarity_df.values
    # Convert the similarity matrix to a condensed distance matrix
    np.fill_diagonal(similarity_df.values, 0)
    condensed_dist_matrix = distance.squareform(similarity_df)

    # Compute linkage
    linkage = hierarchy.linkage(condensed_dist_matrix, 'single')

    # Create a dendrogram figure
    dendro = ff.create_dendrogram(linkage, labels=labels, orientation='bottom')
    
    print(dendro)
    
    for i in range(len(dendro['data'])):
        dendro['data'][i]['yaxis'] = 'y2'

    # Create a side dendrogram
    dendro_side = ff.create_dendrogram(linkage, orientation='right')
    for i in range(len(dendro_side['data'])):
        dendro_side['data'][i]['xaxis'] = 'x2'

    # Add the side dendrogram to the dendrogram figure
    for data in dendro_side['data']:
        dendro.add_trace(data)

    # Create a heatmap
    dendro_leaves = [int(i) for i in dendro_side['layout']['yaxis']['ticktext']]
    heatmap_data = data_array[dendro_leaves, :]
    heatmap_data = heatmap_data[:, dendro_leaves]
    heatmap = go.Heatmap(x=dendro['layout']['xaxis']['tickvals'],
                        y=dendro_side['layout']['yaxis']['tickvals'],
                        z=heatmap_data, colorscale='Blues')
    dendro.add_trace(heatmap)

    # Configure layout
    dendro.update_layout({'width':800, 'height':800, 'showlegend':False, 'hovermode': 'closest',})
    dendro.update_layout(xaxis={'domain': [.15, 1], 'mirror': False, 'showgrid': False, 'showline': False, 'zeroline': False, 'ticks':""})
    dendro.update_layout(xaxis2={'domain': [0, .15], 'mirror': False, 'showgrid': False, 'showline': False, 'zeroline': False, 'showticklabels': False, 'ticks':""})
    dendro.update_layout(yaxis={'domain': [0, .85], 'mirror': False, 'showgrid': False, 'showline': False, 'zeroline': False, 'showticklabels': False, 'ticks': ""})
    dendro.update_layout(yaxis2={'domain':[.825, .975], 'mirror': False, 'showgrid': False, 'showline': False, 'zeroline': False, 'showticklabels': False, 'ticks':""})

    # Plot
    # dendro.show()

    # Plot
    pio.write_html(dendro, 'dendro.html')
    
    
    




    ##################### PLOTLY END #####################
    
    
    # serialize distance matrix to binary format
    LOGGER.INFO(f"serializing the distance matrix to {output_prefix}_distmat.pkl")
    similarity_df.to_pickle(f"{output_prefix}_distmat.pkl")
    LOGGER.INFO(f"Writing distance matrix to {output_prefix}_distmat.tsv")
    similarity_df.to_csv(f"{output_prefix}_distmat.tsv", sep='\t')
    

    newick_out = f"{output_prefix}.newick"

    LOGGER.INFO(f"Writing clustermap plot to {output_prefix}_clustermap.png")
    plt.savefig(f"{output_prefix}_clustermap.png", dpi=600)
    

    if newick:
        # Call the function with your similarity DataFrame and 'single' linkage method
        try:
            newick_string = similarity_df_to_newick(similarity_df, 'average')
            with open(newick_out, 'w') as NW:
                    NW.write(newick_string)

        except RecursionError as err:
            LOGGER.ERROR(f"Couldn't handle the tree depth | {err}")

    LOGGER.SUCCESS("Done.")
