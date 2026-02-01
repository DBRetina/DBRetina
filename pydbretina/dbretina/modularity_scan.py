#!/usr/bin/python
# -*- coding: utf-8 -*-
from __future__ import division

from click.decorators import option
import _dbretina_internal as dbretina_internal
import click
from dbretina.click_context import cli
import matplotlib.pyplot as plt
import seaborn as sns
import os
from collections import defaultdict
import csv
import pandas as pd
import dbretina.dbretina_doc_url as dbretina_doc


def path_to_absolute_path(ctx, param, value):
    return value if value == "NA" else os.path.abspath(value)


@cli.command(name="modularity", epilog=dbretina_doc.doc_url("modularity"), help_priority=7)
@click.option('-i', '--index-prefix', 'index_prefix', required=True, type=click.STRING, help="Index file prefix")
@click.option('-p', '--pairwise', 'pairwise_file', callback=path_to_absolute_path, required=True, type=click.Path(exists=True), help="pairwise TSV file")
@click.option('-c', '--cutoff', 'cutoff', required=True, type=click.FloatRange(0, 100, clamp=False), help="containment cutoff")
@click.option('-o', '--output', "output_prefix", required=True, help="output file prefix")
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
    gene_set_to_length = {}
    dbri_path = f"{index_prefix}.dbri"
    if os.path.exists(dbri_path):
        import _dbretina_internal as dbretina_internal
        all_groups = set(dbretina_internal.dbri_load_names_list(dbri_path))
        gene_set_to_length = dbretina_internal.dbri_load_group_feature_counts(dbri_path)
    else:
        namesMap_file = f"{index_prefix}.namesMap"
        with open(namesMap_file, 'r') as f:
            next(f)
            for line in f:
                gene_set_name = line.strip().split('|')[1]
                all_groups.add(gene_set_name)
        featuresCount_file = f"{index_prefix}_groupID_to_featureCount.tsv"
        with open(featuresCount_file, 'r') as f:
            next(f)
            for line in f:
                gene_set_name, length = line.strip().split('\t')
                gene_set_to_length[gene_set_name] = int(length)

    original_groups_count = len(all_groups)    
    
    #################################
    # 2. Pairwise file parsing
    #################################

    metric_to_col = {
        "containment": 5,
        "ochiai": 6,
        "jaccard": 7,
        "csi": 8,
        "dice": 9,
        "odds_ratio": 10,
        "pvalue": 11,
    }

    DISTANCE_COL = metric_to_col["containment"]

    gene_sets_nodes_data = defaultdict(lambda: {'fragmentation': 0, 'heterogeneity': 0})

    LOGGER.INFO(f"parsing the pairwise file: {pairwise_file}")

    # Check for .dbrp binary pairwise file first
    dbrp_path = os.path.splitext(pairwise_file)[0] + ".dbrp"
    if not os.path.exists(dbrp_path):
        # Also try replacing _pairwise.tsv with _pairwise.dbrp
        dbrp_path = pairwise_file.replace("_pairwise.tsv", "_pairwise.dbrp")

    if os.path.exists(dbrp_path):
        LOGGER.INFO(f"found .dbrp file: {dbrp_path}, using binary pairwise reader")
        # containment metric_id = 0
        records = dbretina_internal.dbrp_filter_pairs(dbrp_path, 0, cutoff)
        for rec in records:
            gene_set1 = rec['group_1_name']
            gene_set2 = rec['group_2_name']
            gene_set1_len = gene_set_to_length[gene_set1]
            gene_set2_len = gene_set_to_length[gene_set2]

            if gene_set1_len < gene_set2_len:
                gene_sets_nodes_data[gene_set1]['fragmentation'] -= 1
                gene_sets_nodes_data[gene_set2]['heterogeneity'] += 1
            elif gene_set1_len > gene_set2_len:
                gene_sets_nodes_data[gene_set1]['heterogeneity'] += 1
                gene_sets_nodes_data[gene_set2]['fragmentation'] -= 1
    else:
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
