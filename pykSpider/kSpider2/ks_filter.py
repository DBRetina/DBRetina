#!/usr/bin/python
# -*- coding: utf-8 -*-
from __future__ import division
import sys
import _kSpider_internal as kSpider_internal
import click
from kSpider2.click_context import cli
import subprocess
import os


def execute_bash_command(command):
    try:
        process = subprocess.Popen(
            command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, text=True
        )
        output, error = process.communicate()

        if process.returncode == 0:
            return True
        print("Command execution failed.", file=sys.stderr)
        print(f"Error:\n{error}", file=sys.stderr)
        return False

    except Exception as e:
        print(f"An error occurred: {e}", file=sys.stderr)
        return False


@cli.command(name="filter", help_priority=3)
@click.option('-p', '--pairwise', 'pairwise_file', required=True, type=click.STRING, help="the pairwise TSV file")
@click.option('-g', '--groups-file', "groups_file", required=False, default="NA", type=click.Path(exists=True), help="single-column supergroups file")
@click.option('-d', '--dist-type', "distance_type", required=False, default="max_cont", show_default=True, type=click.STRING, help="select from ['min_cont', 'avg_cont', 'max_cont', 'ochiai', 'jaccard']")
@click.option('-c', '--cutoff', required=False, type=click.FloatRange(0, 100, clamp=False), default=0.0, show_default=True, help="filter out distances < cutoff")
@click.pass_context
def main(ctx, pairwise_file, groups_file, distance_type, cutoff):
    """
    Generate pairwise TSV.
    """

    ctx.obj.INFO(
        f"Filtering the pairwise matrix on the {distance_type} column with a cutoff of {cutoff}.")

    original_file = pairwise_file

    distance_to_col = {
        "min_cont": 5,
        "avg_cont": 6,
        "max_cont": 7,
        "ochiai": 8,
        "jaccard": 9,
    }
    awk_column = distance_to_col[distance_type] + \
        1  # +1 because awk is 1-indexed

    if cutoff > 0.0 and groups_file != "NA":
        ctx.obj.ERROR("Cannot filter by cutoff and groups at the same time.")

    if cutoff != 0.0:
        filtered_file = pairwise_file.replace(
            '.tsv', f"_{distance_type}_{cutoff}%.tsv")

        ctx.obj.INFO(
            f"Filtering by cutoff({cutoff}%) on {distance_type}\nPlease wait...")
        command = f"cat {original_file} | LC_ALL=C awk awk -F'\t' '{{if (${awk_column} >= {cutoff}) print $0}}' > {filtered_file}"
        result = execute_bash_command(command)
    elif groups_file != "NA":
        # filtered by : basename without the unknown extension
        filtered_by = "".join(os.path.basename(groups_file).split('.')[:-1])
        filtered_file = pairwise_file.replace(
            '.tsv', f"_filtered_by_{filtered_by}.tsv")
        ctx.obj.INFO(f"Filtering by groups file {groups_file}\nPlease wait...")
        # awk_script = f"""awk -v a="$(awk '{{a[$1]=1}} END{{for (k in a) printf "{{%s}}\\n", k}}' {groups_file})" 'BEGIN {{ while ( getline < a ) id_map[$1]=1 }} {{ if ( ($3 in id_map) && ($4 in id_map) ) {{ print $0 }} }}' {original_file} > {filtered_file}"""
        # awk_script = f"""awk -v a="$(awk '{{a[tolower($1)]=1}} END{{for (k in a) printf "{{%s}}\\n", k}}' {groups_file})" 'BEGIN {{ while ( getline < a ) id_map[tolower($1)]=1 }} {{ if ( (tolower($3) in id_map) && (tolower($4) in id_map) ) {{ gsub(/"/, "", $3); gsub(/"/, "", $4); print $0 }} }}' {original_file} > {filtered_file}"""
        # awk_script = f"""awk '{{ if ( (tolower($3) in id_map) && (tolower($4) in id_map) ) {{ gsub(/"/, "", $3); gsub(/"/, "", $4); print $0 }} }}' id_map="$(awk '{{a[tolower($1)]=1}} END{{for (k in a) printf "{{%s}}\\n", k}}' {groups_file})" {original_file} > {filtered_file}"""
        awk_script = f"""head -n1 {original_file} > {filtered_file} && awk 'BEGIN {{ while ( getline < "{groups_file}" ) id_map[tolower($1)]=1 }} {{ if ( (tolower($3) in id_map) && (tolower($4) in id_map) ) {{ gsub(/"/, "", $3); gsub(/"/, "", $4); print $0 }} }}' {original_file} >> {filtered_file}"""
        result = execute_bash_command(awk_script)

    ctx.obj.SUCCESS("Done.")
