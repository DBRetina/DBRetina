#!/usr/bin/python
# -*- coding: utf-8 -*-
from __future__ import division

import _kSpider_internal as kSpider_internal
import click
from kSpider2.click_context import cli
import os
import sys
import pandas as pd
import csv
import json

class StringHasher:
    FNV_prime = 1099511628211
    offset_basis = 14695981039346656037

    def __init__(self):
        pass

    def strHasher(self, s):
        h = self.offset_basis
        for byte in s.encode():
            h = h ^ byte
            # Make sure it's a 64-bit number
            h = (h * self.FNV_prime) & 0xFFFFFFFFFFFFFFFF 
        return h



def append_line_to_file(line, file):
    with open(file, "a") as f:
        f.write(line + "\n")


def get_command():
    _sys_argv = sys.argv
    for i in range(len(_sys_argv)):
        if os.path.isfile(_sys_argv[i]):
            _sys_argv[i] = os.path.abspath(_sys_argv[i])
        if _sys_argv[i] == '-o':
            _sys_argv[i+1] = os.path.abspath(_sys_argv[i+1])
    return "DBRetina " + " ".join(_sys_argv[1:])

def gmts_to_association(gmt_paths, tsv_path):
    with open(tsv_path, 'w', encoding="utf-8") as writer:
        writer.write(f"gene_set\tgene\n")
        # Process each GMT file
        for gmt_path in gmt_paths:
            with open(gmt_path, 'r') as f:
                for line in f:
                    split_line = line.strip().split('\t')
                    if len(split_line) < 3:
                        raise ValueError(f"Line '{line.strip()}' in file '{gmt_path}' doesn't adhere to GMT format")
                    
                    set_name = split_line[0]
                    genes = split_line[2:]
                    for gene in genes:
                        writer.write(f"{set_name}\t{gene}\n")

def sketch(association_file, output_prefix):
    def fnv1a_64(s: str) -> str:
        FNV_prime = 1099511628211
        offset_basis = 14695981039346656037

        h = offset_basis
        for byte in s.encode():
            h = h ^ byte
            h = (h * FNV_prime) & 0xFFFFFFFFFFFFFFFF  # Make sure it's a 64-bit number
        return str(h)

    
    
    # load association file without the header
    association_df = pd.read_csv(association_file, sep="\t", header=None, skiprows=1)
    
    association_df.columns = ["gene_set", "gene"]
    
    # everything to lowercase
    association_df["gene_set"] = association_df["gene_set"].str.lower()
    association_df["gene"] = association_df["gene"].str.lower()
    
    raw_json_dict = {"metadata": {}}
    raw_json_dict["metadata"]["filetype"] = "private"
    raw_json_dict["data"] = {}
    # raw_json_dict["data"] to dict groupby gene_set
    raw_json_dict["data"] = association_df.groupby("gene_set")["gene"].apply(list).to_dict()
    
    # write raw_json_dict to json
    raw_json_path = output_prefix + "_raw.json"
    with open(raw_json_path, "w") as f:
        json.dump(raw_json_dict, f)
    
    private_json_dict = {"metadata": {}}
    private_json_dict["metadata"]["filetype"] = "public"
    private_json_dict["data"] = {}
    # hash the gene set names
    association_df["gene"] = association_df["gene"].apply(lambda x: fnv1a_64(x))    
    association_df.to_csv(output_prefix + "FUCK.tsv", sep="\t", index=False)
    private_json_dict["data"] = association_df.groupby("gene_set")["gene"].apply(list).to_dict()
    
    private_json_path = output_prefix + "_hashes.json"
    with open(private_json_path, "w") as f:
        json.dump(private_json_dict, f)
    
    


@cli.command(name="index", help_priority=1)
@click.option('-a', '--asc', "asc_file", required=False, type=click.Path(exists=True), help="associations file col1: gene_set, col2: single gene. 1st line is header.")
@click.option('-g', '--gmt', "gmt_file", required=False, type=click.Path(exists=True), help="GMT file")
# @click.option('-n', '--names', "names_file", required=False, type=click.Path(exists=True), help="names file")
@click.option('-o', '--output', "output_prefix", required=True, help="output file prefix")
@click.pass_context
def main(ctx, asc_file, output_prefix, gmt_file):
    """
    Index the input data files.
    """
    
    # at least one of asc_file or gmt_file must be provided
    if not asc_file and not gmt_file:
        ctx.obj.ERROR("At least one of asc_file or gmt_file must be provided")
    
    # can't provide both asc_file and gmt_file
    if asc_file and gmt_file:
        ctx.obj.ERROR("Can't provide both asc_file and gmt_file")


    # if not names_file:
    #     names_file = "NA"
    names_file = "NA"


    if gmt_file:
        asc_file = f"generated_{output_prefix}_gmt_to_asc.tsv"
        gmts_to_association([gmt_file], asc_file)


    ctx.obj.INFO("Sketching in progress, please wait...")
    # kSpider_internal.sketch_dbretina(asc_file, names_file, output_prefix)
    sketch(asc_file, output_prefix)
    
    if "generated_" in asc_file and gmt_file:
        os.remove(asc_file)

    json_file = f"{output_prefix}_hashes.json"
    ctx.obj.SUCCESS("File(s) has been sketched.")
    ctx.obj.INFO("Indexing in progress, please wait...")
    kSpider_internal.dbretina_indexing(json_file, output_prefix)
    append_line_to_file(f"command:{get_command()}", f"{output_prefix}.extra")
    ctx.obj.SUCCESS("DONE!")
