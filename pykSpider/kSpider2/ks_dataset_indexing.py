#!/usr/bin/python
# -*- coding: utf-8 -*-
from __future__ import division

import _kSpider_internal as kSpider_internal
import click
from kSpider2.click_context import cli
import os
from glob import glob


@cli.command(name="index", help_priority=2)
@click.option('-j', '--json', "json_file", required=True, help="hashes json file)")
@click.pass_context
def main(ctx, json_file):
    """
    Index hashes JSON file.
    """
    if not os.path.exists(json_file):
        ctx.obj.ERROR(f"{json_file} does not exist!")

    ctx.obj.INFO(
        f"Indexing hashes json file {json_file}")
    kSpider_internal.dbretina_indexing(json_file)
    ctx.obj.SUCCESS("DONE!")
