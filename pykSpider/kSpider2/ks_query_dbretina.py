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
import sys
import json
import hashlib

def hash_string_to_unsigned_int64(input_string):
    hash_object = hashlib.sha256(input_string.encode())
    hex_dig = hash_object.hexdigest()
    # Get the first 16 hexadecimal digits (64 bits) and convert to an integer
    return int(hex_dig[:16], 16)

def plot_histogram(features_counts, output_file):
    # Set style and context to make a nicer plot
    sns.set_style("whitegrid")
    # sns.set_context("talk")

    plt.figure(figsize=(10, 6))  # Set the figure size
    plot = sns.histplot(features_counts, color='skyblue', edgecolor='black',
                        stat='count', bins=50, discrete=True)  # Generate histogram with KDE

    plt.title('Histogram of features frequencies')  # Set the title
    plt.xlabel('Feature frequency')  # Set the x-label
    plt.ylabel('Count (log scale)')  # Set the y-label
    plt.yscale('log')

# Add a legend
    # plot.legend(labels=['Cluster Sizes'])
    # plt.show()
    plt.savefig(output_file, dpi=500)


def validate_numbers(ctx, param, value):
    if not len(value):
        return []
    try:
        return [int(num) for num in value.split(',')]
    except ValueError as e:
        raise click.BadParameter(
            'Numbers must be a comma-separated list of integers') from e


def path_to_absolute_path(ctx, param, value):
    return value if value == "NA" else os.path.abspath(value)


def inject_index_command(index_prefix):
    extra_file = f"{index_prefix}.extra"
    if not os.path.exists(extra_file):
        return ""
    with open(extra_file, "r") as f:
        for line in f:
            line = line.strip().split(":")
            if line[0] == "command":
                return line[1]
        return ""


def get_command():
    _sys_argv = sys.argv
    for i in range(len(_sys_argv)):
        if _sys_argv[i] == "-i":
            _sys_argv[i+1] = os.path.abspath(_sys_argv[i+1])
    return "#command: DBRetina " + " ".join(_sys_argv[1:])


def invert_json(json_file):
    inverted_json = {}
    for key, values in json_file.items():
        for value in values:            
            if value not in inverted_json:
                inverted_json[value] = [key]
            else:
                inverted_json[value].append(key)
    return inverted_json


@cli.command(name="query", help_priority=6)
@click.option('-i', '--index-prefix', "index_prefix", required=True, type=click.STRING, help="index file prefix")
@click.option('-g', '--groups-file', "groups_file", callback=path_to_absolute_path, required=False, default="NA", type=click.Path(exists=False), help="single-column supergroups file")
@click.option('--clusters-file', "clusters_file", callback=path_to_absolute_path, required=False, default="NA", type=click.Path(exists=False), help="DBRetina clusters file")
@click.option('--cluster-ids', "cluster_ids", callback=validate_numbers, required=False, default="", help="comma-separated list of cluster IDs")
@click.option('-o', '--output', "output_prefix", required=True, default=None, help="output file prefix")
@click.pass_context
def main(ctx, groups_file, clusters_file, cluster_ids, index_prefix, output_prefix):
    """Query DBRetina index.

Detailed description:


    Query a DBRetina index with a set of groups (provided as a single-column file or cluster IDs in a DBRetina cluster file). Output each feature and the associated supergroups.
    

Examples:

    1- groups file                    | DBRetina query -i index_prefix -g groups_file -o output_prefix

    
    2- clusters file with cluster IDs | DBRetina query -i index_prefix --clusters-file clusters_file --cluster-ids 1,2,3 -o output_prefix

    
    """
    
    
    # Inverting the index
    inverted_index_prefix = f"inverted_{index_prefix}"
    phmap_file = f"{inverted_index_prefix}.phmap"
    if not os.path.exists(phmap_file):
        # inverted index not found
        raw_json_file = f"{index_prefix}_raw.json"
        # load json in a dictionary:
        with open(raw_json_file, "r") as f:
            supergroups_to_features = json.loads(f.read())["data"]

        inverted_supergroups_to_features = invert_json(supergroups_to_features)

        new_json_raw_dict = {"data": {}, "metadata": {"filetype": "private"}}
        new_json_hashes_dict = {"data": {}, "metadata": {"filetype": "public"}}        

        for key, values in inverted_supergroups_to_features.items():
            if key not in new_json_raw_dict["data"]:
                new_json_raw_dict["data"][key] = []
                new_json_hashes_dict["data"][key] = []
            for value in values:                
                new_json_raw_dict["data"][key].append(value)
                new_json_hashes_dict["data"][key].append(str(hash_string_to_unsigned_int64(value)))

        new_json_raw_file = f"{inverted_index_prefix}_raw.json"
        new_json_hashes_file = f"{inverted_index_prefix}_hashes.json"
        
        with open(new_json_raw_file, "w") as f:
            json.dump(new_json_raw_dict, f)
        with open(new_json_hashes_file, "w") as f:
            json.dump(new_json_hashes_dict, f)
            
        # Create the inverted index
        kSpider_internal.dbretina_indexing(new_json_raw_file, inverted_index_prefix)
        # index_prefix = inverted_index_prefix

    # if all are NA
    if groups_file == "NA" and clusters_file == "NA" and not len(cluster_ids):
        ctx.obj.ERROR("DBRetina's query command requires a groups_file or (clusters_file and cluster_ids).")

    # if clusters_file then must be cluster_id
    if clusters_file != "NA" and not len(cluster_ids):
        ctx.obj.ERROR(
            "DBRetina's filter command requires cluster_id(s) if clusters_file is provided.")
    elif clusters_file == "NA" and len(cluster_ids):
        ctx.obj.ERROR(
            "DBRetina's filter command requires a clusters_file if cluster_id(s) is provided.")

    # can't filter by groups_file and clusters_file at the same time
    if groups_file != "NA" and clusters_file != "NA":
        ctx.obj.ERROR(
            "DBRetina's filter command can't filter by groups_file and clusters_file at the same time.")

    # commands = inject_index_command(index_prefix) + '\n' + get_command()
    commands = get_command()


    _tmp_file = ".DBRetina.tmp.group"
    query_file = ""
    all_ids = set()
    if clusters_file != "NA":
        query_file = _tmp_file
        with (open(clusters_file) as f, open(query_file, 'w') as W):            
            # skip comments
            while True:
                pos = f.tell()
                line = f.readline()
                if not line.startswith('#'):
                    f.seek(pos)
                    break
            # skip header
            next(f)
            for line in f:
                line = line.strip().split('\t')
                cluster_id = int(line[0])
                if cluster_id in cluster_ids:
                    all_ids.add(cluster_id)
                    W.write(line[2].replace('|', '\n') + '\n')

        unfound_ids = set(cluster_ids).difference(all_ids)
        if len(unfound_ids):
            ctx.obj.WARNING(
                f"Couldn't find the following cluster IDs: {unfound_ids}")

    else:
        query_file = groups_file


    features_to_groups_file = f"{output_prefix}_feature_to_groups.tsv"
    counts_file = f"{output_prefix}_features_count_per_group.tsv"
    kSpider_internal.query(index_prefix, inverted_index_prefix, query_file, output_prefix, commands)
    ctx.obj.INFO(
        f"writing query results to {features_to_groups_file}, and {counts_file}")

    # if _tmp_file exists, remove it
    if os.path.exists(_tmp_file):
        os.remove(_tmp_file)

    features_counts = []
    with open(counts_file) as f:
        for line in f:
            if not line.startswith('#'):
                break
                        
        features_counts.extend(int(line.strip().split('\t')[1]) for line in f)

    output_file = f"{output_prefix}_features_count_per_group_histogram.png"
    ctx.obj.INFO(
        f"Plotting histogram of features frequencies to {output_file}")
    plot_histogram(features_counts, output_file)


    ctx.obj.SUCCESS("Query done!")

