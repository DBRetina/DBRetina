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
from kSpider2.ks_pairwise import main as ks_pairwise
from kSpider2.ks_dataset_indexing import main as ks_index
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


def plot_psi_vs_ppi(pathways_PSI, pathways_PPI, output_file="PSI_vs_PPI.png"):
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
    plt.savefig(output_file, dpi=500)
    


@cli.command(name="dedup", help_priority=7)
@click.option('-i', '--index-prefix', required=True, type=click.STRING, help="Index file prefix")
@click.option('--clusters-file', "clusters_file", callback=path_to_absolute_path, required=False, default="NA", type=click.Path(exists=False), help="DBRetina clusters file")
@click.pass_context
def main(ctx, index_prefix, clusters_file):
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
    
    
    genes_ppi_psi_file = f"original_{index_prefix}_genes_ppi_psi.tsv"
    geneSet.export_genes_to_ppi_psi_tsv(genes_ppi_psi_file)
    

    # export genes ppi & psi to TSV
    pathways_PSI_file = f"original_{index_prefix}_pathways_PSI.tsv"
    geneSet.export_genes_to_ppi_psi_tsv(pathways_PSI_file)    
    pathway_to_length = geneSet.get_pathway_lengths()
    
    with open(f"{index_prefix}_raw.json", "r") as f:
        pathway_to_genes = json.loads(f.read())["data"]



    # Plot if needed
    plot_psi_vs_ppi(pathways_PSI, pathways_PPI, f"original_{index_prefix}_PSI_vs_PPI.png")
    
    
    """
    2. Create 80% Community Ochiai
    """
    
    pairwise_file = index_prefix + "_DBRetina_pairwise.tsv"
    cutoff = 80
    distance_type = "ochiai"
    output_prefix = index_prefix + "second_step_80_ochiai"
    community_detection = False
    
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
    
    original_full_pathways = dict()

    second_step_groups_file_for_filter = f"{output_prefix}_groups.tsv"
    with open(second_step_groups_file_for_filter, 'w') as CSV:
        for cluster_id, cluster in second_step_clusters.items():
            pathway_lengths_before_deduplication += len(cluster)

            # cluster_pathway_to_ppi = {pathway: pathways_PPI[pathway] for pathway in cluster}
            cluster_pathway_to_ppi = {}
            for pathway in cluster:
                # to keep tracking of the original pathways then duplicates
                original_full_pathways[pathway] = None
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
    5. Calculate heterogeneity and fragmentation from the max_cont 80% edges
    """
    
    Logger.INFO("Calculating heterogeneity and fragmentation from the max_cont 80% edges")
    
    # geneSet.calculate_heterogeneity_and_fragmentation_from_pairwise(pairwise_file)
    geneSet.calculate_heterogeneity_and_fragmentation_from_pairwise(f"{step_two_filter_pairwise_file_prefix}.tsv")
    pathway_to_fragmentation = geneSet.get_pathway_to_fragmentation()
    pathway_to_heterogeneity = geneSet.get_pathway_to_heterogeneity()
    pathway_to_modularity = geneSet.get_pathway_to_modularity()
    pathway_to_pcsi = geneSet.get_pathways_pcsi()
    
    
    
    # pathways_PPI_PSI to TSV
    pathways_metadata_file = f"{index_prefix}_pathways_metadata.tsv"
    with open(pathways_metadata_file, "w") as f:
        f.write("pathway\tPPI\tPSI\tPCSI\tfragmentation\theterogeneity\tmodularity\tlength\n")
        for pathway, ppi in pathways_PPI.items():
            psi = pathways_PSI[pathway]
            ppi = pathways_PPI[pathway]
            pcsi = pathway_to_pcsi[pathway]
            fragmentation = pathway_to_fragmentation.get(pathway, 'NA')
            heterogeneity = pathway_to_heterogeneity.get(pathway, 'NA')
            modularity = pathway_to_modularity.get(pathway, 'NA')
            length = pathway_to_length[pathway]
            f.write(f"{pathway}\t{ppi}\t{psi}\t{pcsi}\t{fragmentation}\t{heterogeneity}\t{modularity}\t{length}\n")


    """
    6. Set cover
    """
    
    Logger.INFO("Step 3: Set cover")
    
    GC = 100
    cluster_id_to_pathways = cluster_to_pathways(second_step_clusters_file)
    original_community_ids = list(cluster_id_to_pathways)
    new_associations_file = f"{output_prefix}_deduplicated_associations.tsv"
    deduplicated_cluster_info_file = f"{output_prefix}_best_pathways_cluster_final_deduplicates.tsv"
    
    # perform deduplication
    selected_pathways_to_coverage = dict()
    for cluster_id in tqdm(original_community_ids):
        _selected_pathways = geneSet.non_iterative_set_cover(cluster_id, GC)
        selected_pathways_to_coverage.update(_selected_pathways)
 
    with open(deduplicated_cluster_info_file, 'w') as TSV, open (new_associations_file, 'w') as new_Assciations:
        TSV_header = '\t'.join(["Pathway","coverage%","Length","PSI","PPI","fragmentation","heterogeneity","modularity","clust_id","passed","genes"])
        TSV.write(TSV_header)
        new_Assciations.write("Pathway\tGene\n")
        for cluster_id, cluster_pathways in tqdm(cluster_id_to_pathways.items()):
            for pathway in cluster_pathways:
                pathway_coverage = selected_pathways_to_coverage.get(pathway, 0.0)
                pathway_length = pathway_to_length[pathway]
                pathway_psi = pathways_PSI[pathway]
                pathway_ppi = pathways_PPI[pathway]
                pathway_fragmentation = pathway_to_fragmentation.get(pathway, 'NA')
                pathway_heterogeneity = pathway_to_heterogeneity.get(pathway, 'NA')
                pathway_modularity = pathway_to_modularity.get(pathway, 'NA')
                pathway_passed = pathway_coverage > 0
                genes = pathway_to_genes[pathway]
                
                TSV.write(f"{pathway}\t{pathway_coverage}\t{pathway_length}\t{pathway_psi}\t{pathway_ppi}\t{pathway_fragmentation}\t{pathway_heterogeneity}\t{pathway_modularity}\t{cluster_id}\t{pathway_passed}\t{','.join(genes)}\n")
                
                # for the new associations file
                if pathway_passed:
                    for gene in genes:
                        new_Assciations.write(f"{pathway}\t{gene}\n")
    
    """
    Deduplication assessement
        Export genes PPI and PSI after deduplication
    """
    
    # Index the new associations file
    dedup_index_prefix = f"idx_dedup_{index_prefix}"
    dedup_hashes_json_file = f"{dedup_index_prefix}_hashes.json"
    ctx.invoke(ks_index, asc_file = new_associations_file, output_prefix = dedup_index_prefix)
    dedup_geneSet = kSpider_internal.GeneSets(dedup_hashes_json_file)
    dedup_geneSet.build_from_index(dedup_index_prefix)
    
    
    # cluster at ochiai 30%
    ctx.invoke(ks_pairwise, 
               index_prefix = dedup_index_prefix,
               user_threads = 16,
    )

    dedup_pairwise_file = dedup_index_prefix + "_DBRetina_pairwise.tsv"

    
    dedup_ochiai_30_clusters_file_prefix = f"{dedup_index_prefix}_ochiai_30"
    ctx.invoke(ks_clustering, 
               pairwise_file = dedup_pairwise_file, 
               output_prefix = dedup_ochiai_30_clusters_file_prefix,
               distance_type = "ochiai",
               cutoff = 30,
               community = True,
               )

    

    dedup_geneSet.build_from_clusters_file(dedup_ochiai_30_clusters_file_prefix + "_clusters.tsv")
    dedup_geneSet.calculate_heterogeneity_and_fragmentation_from_pairwise(dedup_pairwise_file)
    dedup_pathway_to_fragmentation = dedup_geneSet.get_pathway_to_fragmentation()
    dedup_pathway_to_heterogeneity = dedup_geneSet.get_pathway_to_heterogeneity()
    dedup_pathway_to_modularity = dedup_geneSet.get_pathway_to_modularity()
    dedup_pathways_PSI = dedup_geneSet.get_pathways_psi()
    dedup_pathways_PPI = dedup_geneSet.get_pathways_ppi()
    dedup_pathways_PCSI = dedup_geneSet.get_pathways_pcsi()
    dedup_pathway_to_length = dedup_geneSet.get_pathway_lengths()
    dedup_pathways_metadata_file = f"dedup_{index_prefix}_pathways_metadata.tsv"
    
    with open(dedup_pathways_metadata_file, "w") as f:
        f.write("pathway\tPPI\tPSI\tPCSI\tfragmentation\theterogeneity\tmodularity\tlength\n")
        for pathway, ppi in dedup_pathways_PPI.items():
            psi = dedup_pathways_PSI[pathway]
            ppi = dedup_pathways_PPI[pathway]
            pcsi = dedup_pathways_PCSI[pathway]
            fragmentation = dedup_pathway_to_fragmentation.get(pathway, 'NA')
            heterogeneity = dedup_pathway_to_heterogeneity.get(pathway, 'NA')
            modularity = dedup_pathway_to_modularity.get(pathway, 'NA')
            length = dedup_pathway_to_length[pathway]
            f.write(f"{pathway}\t{ppi}\t{psi}\t{fragmentation}\t{heterogeneity}\t{modularity}\t{length}\n")
    

    # export genes ppi & psi to TSV
    dedup_genes_ppi_psi_file = f"dedup_{dedup_index_prefix}_genes_ppi_psi.tsv"
    dedup_geneSet.export_genes_to_ppi_psi_tsv(dedup_genes_ppi_psi_file)

    
