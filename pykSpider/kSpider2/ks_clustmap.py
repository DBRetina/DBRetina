#!/usr/bin/python
# -*- coding: utf-8 -*-

# TODO: fix later, or remove.

from __future__ import division
from click.decorators import option
import _kSpider_internal as kSpider_internal
import click
from kSpider2.click_context import cli
import matplotlib.pyplot as plt
import seaborn as sns
import os
from collections import defaultdict
import csv
import pandas as pd
import json
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

def path_to_absolute_path(ctx, param, value):
    return value if value == "NA" else os.path.abspath(value)


@cli.command(name="clustmap", help_priority=10)
@click.option('-i', '--index-prefix', 'index_prefix', required=True, type=click.STRING, help="Index file prefix")
@click.option('-p', '--pairwise', 'pairwise_file', callback=path_to_absolute_path, required=True, type=click.Path(exists=True), help="the pairwise TSV file")
@click.option('-o', '--output', "output_prefix", required=True, default=None, help="output file prefix")
@click.pass_context
def main(ctx, pairwise_file, index_prefix, output_prefix):
    """
    Compute the modularity of gene sets
    """

    LOGGER = ctx.obj
    
    #################################
    # 1. Pairwise file parsing
    #################################
    
    gene_sets = set()    

    LOGGER.INFO(f"parsing the pairwise file: {pairwise_file}")
    with open(pairwise_file, 'r') as f:
        # Skip comment lines
        while True:
            pos = f.tell()  # remember the position
            line = f.readline()
            if not line.startswith("#"):
                f.seek(pos)  # rewind to the position before the line
                break

        # skip the header
        next(f)

        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            gene_set1 = row[2]
            gene_set2 = row[3]
            gene_sets.add(gene_set1)
            gene_sets.add(gene_set2)


   ##########################################
   # 2. Extract all gene sets and their genes
   ##########################################
   
    gene_set_to_genes = {}
    raw_json_file = f"{index_prefix}_raw.json"
    # load json to dict
    with open(raw_json_file, 'r') as f:
        gene_set_to_genes = json.load(f)["data"]
    
    selected_gene_set_to_genes = {
        gene_set: genes for gene_set, genes in gene_set_to_genes.items() if gene_set in gene_sets
    }
    
    
    

    df = pd.DataFrame({
        gene_set_name : pd.Series(genes) for gene_set_name, genes in selected_gene_set_to_genes.items()
    })
    
    df.to_csv("mock_data.csv", index=False)
    
    all_genes = []
    for col in df.columns:
        all_genes.extend(df[col].dropna().tolist())

    # Get unique genes
    all_genes = list(set(all_genes))

    # Create binary matrix with genes as index and columns as lists
    binary_matrix = pd.DataFrame(0, index=all_genes, columns=df.columns)

    # Fill in the matrix with 1 if gene is present in the list
    for col in df.columns:
        binary_matrix.loc[df[col].dropna(), col] = 1

    # Sort the columns by the first two letters of their names
    sorted_columns = sorted(binary_matrix.columns, key=lambda x: x[:2])
    binary_matrix = binary_matrix[sorted_columns]
    

    g = sns.clustermap(binary_matrix, cmap="viridis", linewidths=.5, figsize=(12, 12), col_cluster=False, row_cluster=True)
    g.fig.suptitle('Mock data genes distibution', fontsize=16)
    # turnoff legend
    g.cax.set_visible(False)
    plt.savefig('simulated_data_clustermap.png', dpi=400)
    plt.show()