#!/usr/bin/python
# -*- coding: utf-8 -*-
from __future__ import division

from click.decorators import option
import _kSpider_internal as kSpider_internal
import click
from kSpider2.click_context import cli
import os
import json


@cli.command(name="query", help_priority=6)
@click.option('-q', '--query', "query_file", required=True, type=click.Path(exists=True), help="query line separated file")
@click.option('-i', '--index-prefix', "index_prefix", required=True, type=click.STRING, help="index file prefix")
@click.option('-o', '--output', "output_prefix", required=True, default=None, help="output file prefix")
@click.pass_context
def main(ctx, query_file, index_prefix, output_prefix):
    """
    Query DBRetina index.
    """ 
  
    key_val_suffix = "_feature_to_groups.tsv"
    val_key_suffix = "_features_count_per_group.tsv"

    kSpider_internal.query(index_prefix, query_file, output_prefix)
    ctx.obj.SUCCESS(f"Query results written in {output_prefix}{key_val_suffix} and {output_prefix}{val_key_suffix}.")
