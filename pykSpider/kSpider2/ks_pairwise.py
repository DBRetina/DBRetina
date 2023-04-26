#!/usr/bin/python
# -*- coding: utf-8 -*-
from __future__ import division
import sys
import _kSpider_internal as kSpider_internal
import click
from kSpider2.click_context import cli


@cli.command(name="pairwise", help_priority=2)
@click.option('-i', '--index-prefix', required=True, type=click.STRING, help="Index file prefix")
@click.option('-t', '--threads', "user_threads", default=1, required=False, type=int, help="number of cores")
@click.option('-d', '--dist-type', "distance_type", required=False, default="max_cont", show_default=True, type=click.STRING, help="select from ['min_cont', 'avg_cont', 'max_cont', 'ochiai', 'jaccard']")
@click.option('-c', '--cutoff', required=False, type=click.FloatRange(0, 1, clamp=False), default=-1, show_default=True, help="filter out distances < cutoff")
@click.pass_context
def main(ctx, index_prefix, user_threads, distance_type, cutoff):
    """
    Generate pairwise TSV.
    """
    ctx.obj.INFO(
        f"Constructing the pairwise matrix using {user_threads} cores.")
    kSpider_internal.pairwise(
        index_prefix, user_threads, distance_type, cutoff)
    ctx.obj.SUCCESS("Done.")
