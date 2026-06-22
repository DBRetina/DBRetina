#!/usr/bin/python
# -*- coding: utf-8 -*-
from __future__ import division

import _dbretina_internal as dbretina_internal
import click
from dbretina.click_context import cli
import os
import sys
import gzip
import pandas as pd
from collections import defaultdict
import json
import dbretina.dbretina_doc_url as dbretina_doc

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

def _is_gzipped(path):
    """Detect gzip by the 0x1f 0x8b magic bytes (authoritative), falling back to the
    .gz extension only when the file can't be sniffed. Magic-first so a file mislabeled
    .gz (plain text) is read as text, and gzip content without a .gz extension still works."""
    try:
        with open(path, "rb") as fh:
            return fh.read(2) == b"\x1f\x8b"
    except OSError:
        return path.endswith(".gz")


def _open_gmt(path):
    """Open a (possibly gzipped) GMT file for reading text. Supports .gmt.gz."""
    if _is_gzipped(path):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "r", encoding="utf-8")


def gmts_to_association(ctx, gmt_paths, tsv_path):
    # Track group names already seen so we can WARN (not fail) when a duplicate
    # group name is encountered: its genes are unioned into the same supergroup,
    # which is easy to miss without a notice (issue 036).
    seen_set_names = set()
    duplicate_set_names = []
    with open(tsv_path, 'w', encoding="utf-8") as writer:
        writer.write(f"gene_set\tgene\n")
        for gmt_path in gmt_paths:
            ctx.obj.INFO(f"Processing {gmt_path}")
            with _open_gmt(gmt_path) as f:
                for line in f:
                    if not line.strip():
                        continue  # tolerate blank lines (common in real GMTs)
                    split_line = line.strip().split('\t')
                    if len(split_line) < 3:
                        ctx.obj.ERROR(
                            f"Line '{line.strip()}' in file '{gmt_path}' doesn't adhere to GMT format"
                        )

                    set_name = split_line[0].replace('"', '')
                    if '|' in set_name:
                        ctx.obj.ERROR(
                            f"Pipe character '|' is not allowed in group names: '{split_line[0]}'"
                        )
                    # Names are lowercased downstream before the union, so detect
                    # duplicates case-insensitively to match the real merge.
                    set_name_key = set_name.lower()
                    if set_name_key in seen_set_names:
                        duplicate_set_names.append(set_name_key)
                    else:
                        seen_set_names.add(set_name_key)
                    genes = split_line[2:]
                    for gene in genes:
                        gene = gene.replace('"', '')
                        if '|' in gene:
                            ctx.obj.ERROR(
                                f"Pipe character '|' is not allowed in gene names: '{gene}'"
                            )
                        writer.write(f"{set_name}\t{gene}\n")
    if duplicate_set_names:
        unique_dupes = sorted(set(duplicate_set_names))
        ctx.obj.WARNING(
            f"Duplicate group name(s) encountered; their features were unioned "
            f"into a single supergroup: {', '.join(unique_dupes)}"
        )
                        
def fnv1a_64(s: str) -> str:
    FNV_prime = 1099511628211
    offset_basis = 14695981039346656037

    h = offset_basis
    for byte in s.encode():
        h = h ^ byte
        h = (h * FNV_prime) & 0xFFFFFFFFFFFFFFFF  # Make sure it's a 64-bit number
    return str(h)


def build_gene_set_json(ctx, association_files, output_prefix):
    # default dictionary string to list of 
    gene_set_to_genes = defaultdict(list)
    for asc in association_files:
        with open(asc) as asc_reader:
            next(asc_reader)
            for line in asc_reader:
                line = line.strip().lower().split('\t')
                group = line[0].replace('"', '')
                gene = line[1].replace('"', '')
                if '|' in group or '|' in gene:
                    ctx.obj.ERROR(
                        f"Pipe character '|' is not allowed in group or gene names: '{line[0]}\t{line[1]}'"
                    )
                gene_set_to_genes[group].append(gene)
                
    # create dataframe
    association_df = pd.DataFrame(dict([(k,pd.Series(v)) for k,v in gene_set_to_genes.items()]))
    association_df = association_df.melt(var_name='gene_set', value_name='gene')
    association_df = association_df.dropna()
    
    # --------------------- raw
    raw_json_dict = {"metadata": {}}
    raw_json_dict["metadata"]["filetype"] = "private"
    raw_json_dict["data"] = {}
    raw_json_dict["data"] = association_df.groupby("gene_set")["gene"].apply(list).to_dict()
    
    # to json
    raw_json_path = output_prefix + "_raw.json"
    with open(raw_json_path, 'w') as JSON_WRITER:
        JSON_WRITER.write(json.dumps(raw_json_dict))
    
    # --------------------- hashes
    hashes_json_dict = {"metadata": {}}
    hashes_json_dict["metadata"]["filetype"] = "public"
    association_df["gene"] = association_df["gene"].apply(lambda x: fnv1a_64(x))
    hashes_json_dict["data"] = association_df.groupby("gene_set")["gene"].apply(list).to_dict()
    
    hashes_json_path = output_prefix + "_hashes.json"
    with open(hashes_json_path, "w") as f:
        json.dump(hashes_json_dict, f)
    

    
def validate_all_files_exist(ctx, param, value):
    if value is None:
        return None
    for path in value:
        if not os.path.exists(path):
            raise click.BadParameter(f"File '{path}' doesn't exist")
    return value


@cli.command(name="index", help_priority=1, epilog=dbretina_doc.doc_url("index"))
@click.option('-a', '--asc', "asc_file", multiple=True, required=False, callback = validate_all_files_exist , help="associations file col1: supergroup, col2: single feature. 1st line is header.")
@click.option('-g', '--gmt', "gmt_file", multiple=True, required=False, callback = validate_all_files_exist, help="GMT file(s)")
# @click.option('-n', '--names', "names_file", required=False, type=click.Path(exists=True), help="names file")
@click.option('-o', '--output', "output_prefix", required=True, help="output file prefix")
@click.pass_context
def main(ctx, asc_file, output_prefix, gmt_file):
    """
    Index the input data files.
    """
    
    # create output_path directories if they don't exist
    parent_directories = os.path.dirname(os.path.abspath(output_prefix))    
    if not os.path.exists(parent_directories):
        print(f"Creating output directory {parent_directories}")
        os.makedirs(parent_directories)
        
    
    # at least one of asc_file or gmt_file must be provided
    if not asc_file and not gmt_file:
        ctx.obj.ERROR("At least one of asc_file or gmt_file must be provided")
    
    # can't provide both asc_file and gmt_file
    if asc_file and gmt_file:
        ctx.obj.ERROR("Can't provide both association and GMT files")


    # if not names_file:
    #     names_file = "NA"
    names_file = "NA"

    asc_from_gmt = False
    if gmt_file:
        output_dir = os.path.dirname(os.path.abspath(output_prefix))
        base = os.path.basename(output_prefix)
        asc_file = os.path.join(output_dir, f"generated_{base}_gmt_to_asc.tsv")
        gmts_to_association(ctx, list(gmt_file), asc_file)
        asc_file = [asc_file]
        asc_from_gmt = True


    ctx.obj.INFO("Indexing in progress, please wait...")
    # dbretina_internal.sketch_dbretina(asc_file, names_file, output_prefix)
    # sketch(asc_file, output_prefix)
    build_gene_set_json(ctx, asc_file, output_prefix)
    
    if asc_from_gmt:
        os.remove(asc_file[0])

    json_file = f"{output_prefix}_hashes.json"
    ctx.obj.SUCCESS("File(s) have been indexed.")
    ctx.obj.INFO("Indexing in progress, please wait...")
    dbretina_internal.dbretina_indexing(json_file, output_prefix)
    ctx.obj.SUCCESS(f"Index written to {output_prefix}.dbri")
