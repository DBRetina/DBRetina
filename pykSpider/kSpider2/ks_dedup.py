#!/usr/bin/python
# -*- coding: utf-8 -*-
from __future__ import division
import sys
import kSpider_internal
import click
from kSpider2.click_context import cli
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os
import json

def path_to_absolute_path(ctx, param, value):
    return value if value == "NA" else os.path.abspath(value)

def get_command():
    _sys_argv = sys.argv
    for i in range(len(_sys_argv)):
        if _sys_argv[i] == "-i":
            _sys_argv[i+1] = os.path.abspath(_sys_argv[i+1])
    return "#command: DBRetina " + " ".join(_sys_argv[1:])


@cli.command(name="dedup", help_priority=7)
@click.option('-i', '--index-prefix', required=True, type=click.STRING, help="Index file prefix")
@click.option('-a', '--asc', "asc_file", required=True, type=click.Path(exists=True), help="associations file")
@click.option('--clusters-file', "clusters_file", callback=path_to_absolute_path, required=False, default="NA", type=click.Path(exists=False), help="DBRetina clusters file")
@click.option('-d', '--dist-type', "distance_type", required=False, default="max_cont", show_default=True, type=click.STRING, help="select from ['min_cont', 'avg_cont', 'max_cont', 'ochiai', 'jaccard', 'pval']")
@click.pass_context
def main(ctx, index_prefix, clusters_file, distance_type, asc_file):
    """
    Deduplicate a given clusters file.
    """
    
    
    """
        1. Read the clusters file
    """
    hashes_json_file = f"{index_prefix}_hashes.json"
    geneSet = kSpider_internal.GeneSets(hashes_json_file)
    geneSet.build_from_index(index_prefix)
    geneSet.build_from_clusters_file(clusters_file)
    pathways_PSI = geneSet.get_pathways_psi()
    pathways_PPI = geneSet.get_pathways_ppi()
    
    
    
    """
        2. Calculate Pathway Pleiotropy Index PPI for clusters' pathways
    """
    
    
    """
        3. Calculate Pathway Specificity Index PSI for pathways from the whole dataset
    """
    
    
    # DEBUG -----------------------------------------------
    # pathways_PSI dict to json file
    with open(f"{index_prefix}_pathways_PSI.json", "w") as f:
        f.write(json.dumps(pathways_PSI))
    
    # pathways_PPI dict to json file
    with open(f"{index_prefix}_pathways_PPI.json", "w") as f:
        f.write(json.dumps(pathways_PPI))
    # END DEBUG -------------------------------------------
    
    ctx.obj.SUCCESS("Done.")
