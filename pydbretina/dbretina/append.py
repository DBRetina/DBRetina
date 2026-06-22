#!/usr/bin/python
# -*- coding: utf-8 -*-
from __future__ import division

import _dbretina_internal as dbretina_internal
import click
from dbretina.click_context import cli
import os
from dbretina.dataset_indexing import build_gene_set_json, gmts_to_association, validate_all_files_exist
import dbretina.dbretina_doc_url as dbretina_doc


@cli.command(name="append", help_priority=1)
@click.option('-i', '--index', "existing_index", required=True, type=click.Path(exists=True), help="Existing .dbri index file")
@click.option('-a', '--asc', "asc_file", multiple=True, required=False, callback=validate_all_files_exist, help="associations file col1: supergroup, col2: single feature. 1st line is header.")
@click.option('-g', '--gmt', "gmt_file", multiple=True, required=False, callback=validate_all_files_exist, help="GMT file(s)")
@click.option('-o', '--output', "output_path", required=True, help="output .dbri file path")
@click.pass_context
def main(ctx, existing_index, asc_file, gmt_file, output_path):
    """
    Append new gene sets to an existing .dbri index.
    """

    if not asc_file and not gmt_file:
        ctx.obj.ERROR("At least one of asc_file or gmt_file must be provided")

    if asc_file and gmt_file:
        ctx.obj.ERROR("Can't provide both association and GMT files")

    # Ensure output ends with .dbri
    if not output_path.endswith(".dbri"):
        output_path = output_path + ".dbri"

    # Create output directories if needed
    parent_directories = os.path.dirname(os.path.abspath(output_path))
    if parent_directories and not os.path.exists(parent_directories):
        os.makedirs(parent_directories)

    # Derive a temp prefix for indexing
    output_prefix = output_path.rsplit(".dbri", 1)[0]

    asc_from_gmt = False
    if gmt_file:
        output_dir = os.path.dirname(os.path.abspath(output_prefix))
        asc_file = os.path.join(output_dir, f"generated_{os.path.basename(output_prefix)}_gmt_to_asc.tsv")
        gmts_to_association(ctx, list(gmt_file), asc_file)
        asc_file = [asc_file]
        asc_from_gmt = True

    ctx.obj.INFO("Indexing new gene sets...")
    build_gene_set_json(ctx, asc_file, output_prefix)

    if asc_from_gmt:
        os.remove(asc_file[0])

    json_file = f"{output_prefix}_hashes.json"
    raw_json_file = f"{output_prefix}_raw.json"
    ctx.obj.INFO(f"Appending to existing index {existing_index}...")
    dbretina_internal.dbretina_append(existing_index, json_file, raw_json_file, output_path)
    ctx.obj.SUCCESS(f"Updated index written to {output_path}")
