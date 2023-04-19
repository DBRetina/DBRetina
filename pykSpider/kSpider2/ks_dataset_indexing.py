#!/usr/bin/python
# -*- coding: utf-8 -*-
from __future__ import division

import _kSpider_internal as kSpider_internal
import click
from kSpider2.click_context import cli
import os
from glob import glob


@cli.command(name="index", help_priority=1)
@click.option('-a', '--asc', "asc_file", required=True, type=click.Path(exists=True), help="ASC file")
@click.option('-n', '--names', "names_file", required=False, type=click.Path(exists=True), help="Optional names file")
@click.option('-o', '--output', "output_prefix", required=False, default=None, help="index output file prefix")
@click.option('--inverted', is_flag=True, help="Inverted index (use this option if you want to query group names and get their genes.)")
@click.pass_context
def main(ctx, asc_file, names_file, output_prefix, inverted):
    """
    Index hashes JSON file.
    """
    
    if not names_file: names_file = "NA"

    if not output_prefix:
        output_prefix = os.path.basename(asc_file)
        output_prefix = os.path.splitext(output_prefix)[0]
        output_prefix = "idx" + "_" + output_prefix

    if inverted:
        output_prefix = f"{output_prefix}_group_to_genes"
    else:
        output_prefix = f"{output_prefix}_gene_to_groups"

    kSpider_internal.sketch_dbretina(asc_file, names_file, inverted)
    asc_basename = os.path.basename(asc_file)
    asc_basename_without_extension = os.path.splitext(asc_basename)[0]
    json_file = f"{asc_basename_without_extension}_public.json"

    ctx.obj.SUCCESS("File(s) has been sketched.")
    
    ctx.obj.INFO(f"Indexing {asc_file}")
    kSpider_internal.dbretina_indexing(json_file, output_prefix)
    ctx.obj.SUCCESS("DONE!")
