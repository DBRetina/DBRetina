#!/usr/bin/python
# -*- coding: utf-8 -*-
from __future__ import division
import sys
import kSpider_internal
import click
from kSpider2.click_context import cli
import matplotlib.pyplot as plt
import seaborn as sns
from kSpider2.ks_clustering import main as ks_clustering
from kSpider2.ks_filter import main as ks_filter
from tqdm import tqdm
import numpy as np
import pandas as pd
import os
import json
import random
from kSpider2.customLogger import Logger

def path_to_absolute_path(ctx, param, value):
    return value if value == "NA" else os.path.abspath(value)

def get_command():
    _sys_argv = sys.argv
    for i in range(len(_sys_argv)):
        if _sys_argv[i] == "-i":
            _sys_argv[i+1] = os.path.abspath(_sys_argv[i+1])
    return "#command: DBRetina " + " ".join(_sys_argv[1:])

def cluster_to_pathways(clusters_file):
    clusters = {}
    header_flag = True
    with open(clusters_file, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            if header_flag:
                header_flag = False
                continue
            cluster_id, cluster_size, cluster_members = line.strip().split("\t")
            clusters[int(cluster_id)] = cluster_members.split("|")

    return clusters


def plot_psi_vs_ppi(pathways_PSI, pathways_PPI):
    # we convert the keys and values to lists so we can plot them
    dict1 = pathways_PSI
    dict2 = pathways_PPI
    keys1 = list(dict1.keys())
    values1 = list(dict1.values())

    keys2 = list(dict2.keys())
    values2 = list(dict2.values())

    # create a figure and axes with plt.subplots
    fig, ax = plt.subplots()

    # set the style to a seaborn one for nicer plots
    sns.set_style("whitegrid")

    # plot both sets of values
    sns.scatterplot(x=keys1, y=values1, ax=ax, color='blue', s=100, label='PSI')
    sns.scatterplot(x=keys2, y=values2, ax=ax, color='red', s=100, label='PPI')

    # add a title and labels
    ax.set_title('Scatter plot of two dictionaries', fontweight='bold', fontsize=14)
    ax.set_xlabel('Keys', fontweight='bold', fontsize=12)
    ax.set_ylabel('Values', fontweight='bold', fontsize=12)

    # add legend
    ax.legend()

    # show the plot
    plt.show()
    
    

    

@cli.command(name="dedup", help_priority=7)
@click.option('-i', '--index-prefix', required=True, type=click.STRING, help="Index file prefix")
@click.option('-a', '--asc', "asc_file", required=True, type=click.Path(exists=True), help="associations file")
@click.option('--clusters-file', "clusters_file", callback=path_to_absolute_path, required=False, default="NA", type=click.Path(exists=False), help="DBRetina clusters file")
@click.pass_context
def main(ctx, index_prefix, clusters_file, asc_file):
    """
    Deduplicate a given clusters file.
    """
    
    Logger = ctx.obj
        
    """
        1. PSI & PPI
            - calculate Pathway Pleiotropy Index PPI for clusters' pathways
            - Calculate Pathway Specificity Index PSI for pathways from the whole dataset
    """

    hashes_json_file = f"{index_prefix}_hashes.json"
    geneSet = kSpider_internal.GeneSets(hashes_json_file)
    geneSet.build_from_index(index_prefix)
    geneSet.build_from_clusters_file(clusters_file)
    pathways_PSI = geneSet.get_pathways_psi()
    pathways_PPI = geneSet.get_pathways_ppi()

    
    # pathways_PPI to TSV
    pathways_PPI_file = f"{index_prefix}_pathways_PPI.tsv"
    with open(pathways_PPI_file, "w") as f:
        f.write("pathway\tPPI\n")
        for pathway, ppi in pathways_PPI.items():
            f.write(f"{pathway}\t{ppi}\n")
    
    pathway_to_length = geneSet.get_pathway_lengths()
    
    with open(f"{index_prefix}_raw.json", "r") as f:
        pathway_to_genes = json.loads(f.read())["data"]


    
    # plot pathways_PSI vs pathways_PPI scatter plot seaborn
    sns.set_theme(style="whitegrid")


    # Plot if needed
    # plot_psi_vs_ppi(pathways_PSI, pathways_PPI)
    
    
    """
    2. Create 80% Community Ochiai
    """
    
    pairwise_file = index_prefix + "_DBRetina_pairwise.tsv"
    cutoff = 80
    distance_type = "ochiai"
    output_prefix = index_prefix + "second_step_80_ochiai"
    community_detection = True
    
    Logger.ACTIVE = False
    ctx.invoke(ks_clustering, 
               pairwise_file = pairwise_file, 
               cutoff = cutoff, 
               distance_type = distance_type, 
               output_prefix = output_prefix, 
               community = community_detection,
               )
    Logger.ACTIVE = True
    
    """
    3. Delete duplicates in the clustering file
    """
    
    second_step_clusters_file = f"{output_prefix}_clusters.tsv"
    second_step_clusters = cluster_to_pathways(second_step_clusters_file)
    cluster_sizes = [len(cluster) for cluster in second_step_clusters.values()]
    plt.yscale('log')
    plt.hist(cluster_sizes, bins=10)
    plt.savefig(f"{output_prefix}_clusters_sizes.png")
    
    pathway_lengths_before_deduplication = 0
    pathway_lengths_after_deduplication = 0

    second_step_groups_file_for_filter = f"{output_prefix}_groups.tsv"
    with open(second_step_groups_file_for_filter, 'w') as CSV:
        for cluster_id, cluster in second_step_clusters.items():
            pathway_lengths_before_deduplication += len(cluster)

            # cluster_pathway_to_ppi = {pathway: pathways_PPI[pathway] for pathway in cluster}
            cluster_pathway_to_ppi = {}
            for pathway in cluster:
                if pathway in pathways_PPI:
                    cluster_pathway_to_ppi[pathway] = pathways_PPI[pathway]
                else:
                    print(f"Pathway {pathway} not in pathways_PPI")

            pathways = list(cluster_pathway_to_ppi.keys())
            pathways.sort(key=lambda x: (cluster_pathway_to_ppi[x], -pathway_to_length[x]))
            min_ppi_len = (cluster_pathway_to_ppi[pathways[0]], pathway_to_length[pathways[0]])
            candidates = [p for p in pathways if (cluster_pathway_to_ppi[p], pathway_to_length[p]) == min_ppi_len]
            chosen_pathway = random.choice(candidates)

            pathway_lengths_after_deduplication += 1
            CSV.write(f"{chosen_pathway}\n")


    deduplication_precentage = (pathway_lengths_before_deduplication - pathway_lengths_after_deduplication) / pathway_lengths_before_deduplication
    deduplication_precentage = round(deduplication_precentage * 100, 2)
    Logger.INFO(f"Step 2 deduplication summary: {pathway_lengths_before_deduplication}-{pathway_lengths_after_deduplication} = {deduplication_precentage}%")

    ## Keep only the remaining pathways in the clusters file
    geneSet.keep_only_these_pathways(second_step_groups_file_for_filter)
    
    
    """
    4. Filter the selected groups file max-containment 80%
    """
    
    
    Logger.ACTIVE = True
    max_cont_cutoff = 80
    step_two_filter_pairwise_file_prefix = f"{output_prefix}_max_cont_{max_cont_cutoff}_pairwise"
    ctx.invoke(
        ks_filter,
        pairwise_file = pairwise_file, 
        groups_file = second_step_groups_file_for_filter,
        distance_type = "max_cont",
        cutoff = max_cont_cutoff, 
        output_file = step_two_filter_pairwise_file_prefix,
        )
    Logger.ACTIVE = True
    
    
    """
    4. Calculate heterogeneity and fragmentation from the max_cont 80% edges
    """
    
    Logger.INFO("Calculating heterogeneity and fragmentation from the max_cont 80% edges")
    
    # geneSet.calculate_heterogeneity_and_fragmentation_from_pairwise(pairwise_file)
    geneSet.calculate_heterogeneity_and_fragmentation_from_pairwise(f"{step_two_filter_pairwise_file_prefix}.tsv")
    pathway_to_fragmentation = geneSet.get_pathways_fragmentation()


    """
    5. Set cover
    """
    
    Logger.INFO("Step 3: Set cover")
    
    GC = 100
    original_community_ids = list(cluster_to_pathways(clusters_file).keys())
    with open(f"{output_prefix}_best_pathways_cluster_final_deduplicates.tsv", 'w') as TSV:
        TSV.write("Pathway\tcoverage%\tLength\tPSI\tPPI\tFragmentation\tclust_id\tgenes\n")
        for cluster_id in tqdm(original_community_ids):
            selected_pathways = geneSet.non_iterative_set_cover(cluster_id, GC)
            for pathway, coverage in selected_pathways.items():
                genes = ','.join(pathway_to_genes[pathway])
                TSV.write(f"{pathway}\t{coverage}\t{pathway_to_length[pathway]}\t{pathways_PSI[pathway]}\t{pathways_PPI[pathway]}\t{pathway_to_fragmentation[pathway]}\t{cluster_id}\t{genes}\n")

    
    """
    # Single cluster test
    cluster_id = 1
    GC = 100
    best_pathways_cluster_1 = geneSet.non_iterative_set_cover(cluster_id, GC);
    with open(f"{output_prefix}_best_pathways_cluster_{cluster_id}.tsv", 'w') as TSV:
        TSV.write("Pathway\tcoverage%\tLength\tPSI\tPPI\tFragmentation\tclust_id\tgenes\n")

        for pathway, coverage in best_pathways_cluster_1.items():
            genes = ','.join(pathway_to_genes[pathway])
            TSV.write(f"{pathway}\t{coverage}\t{pathway_to_length[pathway]}\t{pathways_PSI[pathway]}\t{pathways_PPI[pathway]}\t{pathway_to_fragmentation[pathway]}\t{cluster_id}\t{genes}\n")
    """
    

    
    #################################################################
    exit()

    cluster_id = 13
    GC = 100
    best_pathways_cluster_1 = geneSet.proportionalSetCover(cluster_id, GC);
    # print("best_pathways_cluster_1: ", best_pathways_cluster_1)

    # sort best_pathways_cluster_1 dict by value and print the order
    best_pathways_cluster_sorted = dict(
        sorted(
            best_pathways_cluster_1.items(),
            key=lambda item: item[1],
            reverse=True,
        )
    )
    print("Pathway to score")
    print("best_pathways_cluster_sorted: ", json.dumps(best_pathways_cluster_sorted, indent=4))
    

    # order by rank
    best_pathways_cluster_sorted_rank = {
        key: rank
        for rank, key in enumerate(best_pathways_cluster_sorted, 1)
    }

    print("Pathway to rank")
    print("best_pathways_cluster_sorted_rank: ", json.dumps(best_pathways_cluster_sorted_rank, indent=4))

    # clusters_id_to_pathways = cluster_to_pathways(clusters_file)
    # with open(clusters_file.replace('tsv','_dedup.tsv'), 'w') as F:
    #     for cluster_id, cluster_pathways in clusters_id_to_pathways.items():
    #         print("[DEBUG] cluster_id: ", cluster_id)
    #         deduplicated_pathways = geneSet.proportionalSetCover(cluster_id, GC).keys()
    #         F.write(f"{cluster_id}\t{len(deduplicated_pathways)}\t-{len(cluster_pathways)-len(deduplicated_pathways)}\t{','.join(deduplicated_pathways)}\n")


    # DEBUG -----------------------------------------------
    # pathways_PSI dict to json file
    with open(f"{index_prefix}_pathways_PSI.json", "w") as f:
        f.write(json.dumps(pathways_PSI))

    # pathways_PPI dict to json file
    with open(f"{index_prefix}_pathways_PPI.json", "w") as f:
        f.write(json.dumps(pathways_PPI))
    # END DEBUG -------------------------------------------

    ctx.obj.SUCCESS("Done.")
