#!/usr/bin/python
# -*- coding: utf-8 -*-
from __future__ import division

from click.decorators import option
import _kSpider_internal as kSpider_internal
import click
from kSpider2.click_context import cli
import os
import json


@cli.command(name="sketch", help_priority=3, hidden=True)
@click.option('-q', '--query', "query_file", required=True, type=click.Path(exists=True), help="query line separated file")
@click.option('-i', '--index-prefix', "index_prefix", required=True, type=click.STRING, help="Index file prefix")
@click.option('-o', '--output', "output_prefix", required=False, default=None, help="index output file prefix")
@click.pass_context
def main(ctx, query_file, index_prefix, output_prefix):
    """
    Query DBRetina index.
    """

    inverted = False
    if "_group_to_genes" in index_prefix:
        inverted = True
        key_val_suffix = "_group_to_genes.tsv"
        val_key_suffix = "_gene_to_groupsCount.tsv"  
        ctx.obj.INFO("Inverted index detected. Group names will be expected in the query file.")
        
    else:
        key_val_suffix = "_gene_to_groups.tsv"
        val_key_suffix = "_group_to_genesCount.tsv"
        ctx.obj.INFO("Gene names will be expected in the query file.")
        
    kSpider_internal.query(index_prefix, query_file, output_prefix, inverted)
    ctx.obj.SUCCESS(f"Query results written in {output_prefix}{key_val_suffix} and {output_prefix}{val_key_suffix}.")
