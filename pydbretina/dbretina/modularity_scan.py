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

    # Try Parquet/PairwiseStore first
    try:
        from dbretina.compat import open_pairwise
        store = open_pairwise(pairwise_file)
    except Exception:
        store = None

    if store is not None:
        LOGGER.INFO("using Parquet pairwise data via PairwiseStore")
        names_map = store.get_names_map()
        df = store.to_pandas(
            metric="containment", cutoff=cutoff,
            columns=["group_1_id", "group_2_id", "containment"]
        )
        for _, row in df.iterrows():
            gene_set1 = names_map.get(int(row["group_1_id"]), "")
            gene_set2 = names_map.get(int(row["group_2_id"]), "")
            gene_set1_len = gene_set_to_length.get(gene_set1, 0)
            gene_set2_len = gene_set_to_length.get(gene_set2, 0)

            if gene_set1_len < gene_set2_len:
                gene_sets_nodes_data[gene_set1]['fragmentation'] -= 1
                gene_sets_nodes_data[gene_set2]['heterogeneity'] += 1
            elif gene_set1_len > gene_set2_len:
                gene_sets_nodes_data[gene_set1]['heterogeneity'] += 1
                gene_sets_nodes_data[gene_set2]['fragmentation'] -= 1
        store.close()
    else:
        # Fallback: existing .dbrp / TSV code. Resolve the canonical .dbrp
        # sibling from any -p form (tsv / parquet dir / .dbrp); the old
        # str.replace no-op'd for the dir/.dbrp forms and fed the directory path
        # into the binary reader (issue 046).
        from dbretina.compat import resolve_dbrp_path
        dbrp_path = resolve_dbrp_path(pairwise_file)

        # No store and no .dbrp: a directory input is unreadable here (the TSV
        # open() below would IsADirectoryError). Emit a clean error instead.
        if dbrp_path is None and os.path.isdir(pairwise_file):
            LOGGER.ERROR(
                f"'{pairwise_file}' is a pairwise directory without a usable "
                f"parquet store (no manifest.json) and no sibling .dbrp; pass the "
                f"pairwise TSV, its parquet directory, or the .dbrp to -p.")

        if dbrp_path is not None:
            LOGGER.INFO(f"found .dbrp file: {dbrp_path}, using binary pairwise reader")
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
                while True:
                    pos = f.tell()
                    line = f.readline()
                    if not line.startswith("#"):
                        f.seek(pos)
                        break
                next(f)
                reader = csv.reader(f, delimiter="\t")
                for row in reader:
                    _containment = float(row[DISTANCE_COL])
                    if _containment >= cutoff:
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
    # When no pair passes the cutoff, gene_sets_nodes_data is empty and the
    # DataFrame has no columns. Seed the expected columns so every group is
    # still written with modularity/fragmentation/heterogeneity = 0 (the
    # intent of the unincluded-groups loop below) instead of crashing.
    if df.empty:
        df = pd.DataFrame(columns=['fragmentation', 'heterogeneity'])
        df.index.name = 'gene_set'
        LOGGER.WARNING(f"no pairs passed the cutoff ({cutoff}); all groups get modularity 0")
    # add modularity
    df['modularity'] = abs(df['fragmentation'] + df['heterogeneity'])
    # add unincluded gene sets with modularity, fragmentation and heterogeneity of 0
    unincluded_gene_sets = all_groups - set(df.index)
    
    for gene_set in unincluded_gene_sets:
        df.loc[gene_set] = [0, 0, 0]
    
    LOGGER.INFO(f"Writing the modularity file: {output_prefix}_modularity.tsv")
    df.to_csv(f"{output_prefix}_modularity.tsv", sep='\t', index=True, header=True)
