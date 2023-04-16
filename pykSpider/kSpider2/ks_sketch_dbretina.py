#!/usr/bin/python
# -*- coding: utf-8 -*-
from __future__ import division

from click.decorators import option
import _kSpider_internal as kSpider_internal
import click
from kSpider2.click_context import cli
import os
import json


@cli.command(name="sketch", help_priority=1, hidden=True)
@click.option('-a', '--asc', "asc_file", required=True, type=click.Path(exists=True), help="ASC file")
@click.option('-n', '--names', "names_file", required=False, default="none", type=click.Path(exists=True), help="Optional names file")
@click.pass_context
def main(ctx, asc_file, names_file):
    """
    Sketch a DBRetina files.
    """

    kSpider_internal.sketch_dbretina(asc_file, names_file)

    ctx.obj.SUCCESS("File(s) has been sketched.")
