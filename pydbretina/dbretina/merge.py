#!/usr/bin/python
# -*- coding: utf-8 -*-

import _dbretina_internal as dbretina_internal
import click
from dbretina.click_context import cli
import os
import dbretina.dbretina_doc_url as dbretina_doc


@cli.command(name="merge", help_priority=4, epilog=dbretina_doc.doc_url("merge"))
@click.option('-a', '--index-a', "index_a", required=True, type=click.Path(exists=True), help="First (base) .dbri index file")
@click.option('-b', '--index-b', "index_b", required=True, type=click.Path(exists=True), help="Second .dbri index file to merge in")
@click.option('-o', '--output', "output_path", required=True, help="Output .dbri file path")
@click.pass_context
def main(ctx, index_a, index_b, output_path):
    """
    Merge two .dbri indexes into one.

    Uses index A as the base and incorporates all gene sets from index B.
    Both indexes must contain GROUP_TO_FEATURE_SET sections.
    Gene set names must be unique across both indexes.
    """

    if not output_path.endswith('.dbri'):
        output_path += '.dbri'

    ctx.obj.INFO(f"Merging {index_a} (base) + {index_b}")
    ctx.obj.INFO("Merging in progress, please wait...")

    try:
        dbretina_internal.dbretina_merge(
            os.path.abspath(index_a),
            os.path.abspath(index_b),
            os.path.abspath(output_path)
        )
    except RuntimeError as e:
        # Expected merge-time errors from the C++ core (e.g. duplicate group
        # names across indexes) -> clean [ERROR] line, no Python traceback.
        # ctx.obj.ERROR exits nonzero.
        ctx.obj.ERROR(str(e))

    ctx.obj.SUCCESS(f"Merged index written to {output_path}")
