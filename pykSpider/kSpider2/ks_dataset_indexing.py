#!/usr/bin/python
# -*- coding: utf-8 -*-
from __future__ import division

import _kSpider_internal as kSpider_internal
import click
from kSpider2.click_context import cli
import os
from glob import glob
import sys


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


@cli.command(name="index", help_priority=1)
@click.option('-a', '--asc', "asc_file", required=True, type=click.Path(exists=True), help="associations file")
@click.option('-n', '--names', "names_file", required=False, type=click.Path(exists=True), help="names file")
@click.option('-o', '--output', "output_prefix", required=True, help="output file prefix")
@click.pass_context
def main(ctx, asc_file, names_file, output_prefix):
    """
    Index the input data files.
    """

    if not names_file:
        names_file = "NA"

    if not output_prefix:
        output_prefix = os.path.basename(asc_file)
        output_prefix = os.path.splitext(output_prefix)[0]
        output_prefix = "idx" + "_" + output_prefix

    kSpider_internal.sketch_dbretina(asc_file, names_file, output_prefix)
    asc_basename = os.path.basename(asc_file)
    asc_basename_without_extension = os.path.splitext(asc_basename)[0]
    json_file = f"{output_prefix}_hashes.json"

    ctx.obj.SUCCESS("File(s) has been sketched.")

    ctx.obj.INFO(f"Indexing {asc_file}")
    kSpider_internal.dbretina_indexing(json_file, output_prefix)
    append_line_to_file(f"command:{get_command()}", f"{output_prefix}.extra")
    ctx.obj.SUCCESS("DONE!")
