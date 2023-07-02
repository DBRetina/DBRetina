#!/usr/bin/python
# -*- coding: utf-8 -*-
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



def path_to_absolute_path(ctx, param, value):
    return value if value == "NA" else os.path.abspath(value)


@cli.command(name="modularity", help_priority=7)
@click.option('-i', '--index-prefix', 'index_prefix', required=True, type=click.STRING, help="Index file prefix")
@click.option('-p', '--pairwise', 'pairwise_file', callback=path_to_absolute_path, required=True, type=click.Path(exists=True), help="the pairwise TSV file")
@click.option('-c', '--cutoff', 'cutoff', required=True, type=click.FloatRange(0, 100, clamp=False), help="containment cutoff")
@click.option('-o', '--output', "output_prefix", required=True, default=None, help="output file prefix")
@click.pass_context
def main(ctx, pairwise_file, cutoff, output_prefix, index_prefix):
    """
    Compute the modularity of gene sets
    """

    LOGGER = ctx.obj


    #################################
    # 1. Extract all gene sets and their lengths
    #################################
    all_groups = set()
    namesMap_file = f"{index_prefix}.namesMap"
    with open(namesMap_file, 'r') as f:
        next(f)
        for line in f:
            gene_set_name = line.strip().split('|')[1]
            all_groups.add(gene_set_name)

    original_groups_count = len(all_groups)
    
    
    featuresCount_file = f"{index_prefix}_groupID_to_featureCount.tsv"
    gene_set_to_length = {}
    with open (featuresCount_file, 'r') as f:
        next(f)
        for line in f:
            gene_set_name, length = line.strip().split('\t')
            gene_set_to_length[gene_set_name] = int(length)    
    
    #################################
    # 2. Pairwise file parsing
    #################################

    distance_to_col = {
        "containment": 5,
        "ochiai": 6,
        "jaccard": 7,
        "odds_ratio": 8,
        "pvalue": 9,
    }

    DISTANCE_COL = distance_to_col["containment"]

    gene_sets_nodes_data = defaultdict(lambda: {'fragmentation': 0, 'heterogeneity': 0})

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
            _containment = float(row[DISTANCE_COL])
            if _containment >= cutoff:
                # Extract the gene_set names from columns 3 and 4
                gene_set1 = row[2]
                gene_set2 = row[3]
                gene_set1_len = gene_set_to_length[gene_set1]
                gene_set2_len = gene_set_to_length[gene_set2]
                
                if gene_set1_len < gene_set2_len:
                    gene_sets_nodes_data[gene_set1]['fragmentation'] -= 1
                    gene_sets_nodes_data[gene_set2]['heterogeneity'] += 1
                elif gene_set1_len > gene_set2_len:
                    gene_sets_nodes_data[gene_set1]['heterogeneity'] += 1
                    gene_sets_nodes_data[gene_set2]['fragmentation'] -= 1


    # convert to dataframe
    df = pd.DataFrame.from_dict(gene_sets_nodes_data, orient='index')
    df.index.name = 'gene_set'
    # add modularity
    df['modularity'] = abs(df['fragmentation'] + df['heterogeneity'])
    # add unincluded gene sets with modularity, fragmentation and heterogeneity of 0
    unincluded_gene_sets = all_groups - set(df.index)
    
    for gene_set in unincluded_gene_sets:
        df.loc[gene_set] = [0, 0, 0]
    
    LOGGER.INFO(f"Writing the modularity file: {output_prefix}_modularity.tsv")
    df.to_csv(f"{output_prefix}_modularity.tsv", sep='\t', index=True, header=True)
