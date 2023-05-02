#!/usr/bin/python
# -*- coding: utf-8 -*-
from __future__ import division
import sys
import _kSpider_internal as kSpider_internal
import click
from kSpider2.click_context import cli
import subprocess
import os


def is_awk_available():
    try:
        subprocess.run(["awk"], stdin=subprocess.DEVNULL,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        return True
    except FileNotFoundError:
        return False


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


def get_command():
    _sys_argv = sys.argv
    for i in range(len(_sys_argv)):
        if os.path.isfile(_sys_argv[i]):
            _sys_argv[i] = os.path.abspath(_sys_argv[i])
    return "DBRetina " + " ".join(_sys_argv[1:])


def increment_version(output):
    version = 1
    base_output = output.rsplit('.', 1)[0]
    if not os.path.isfile(output):
        return output
    while True:
        output_version = f"{base_output}_{version}.tsv"
        if not os.path.isfile(output_version):
            return output_version
        version += 1


def validate_numbers(ctx, param, value):
    if not len(value):
        return []
    try:
        return [int(num) for num in value.split(',')]
    except ValueError as e:
        raise click.BadParameter(
            'Numbers must be a comma-separated list of integers') from e


def path_to_absolute_path(ctx, param, value):
    return value if value == "NA" else os.path.abspath(value)



def get_extended_nodes(nodes, kmer_size):
    
    return extended_nodes



@cli.command(name="filter", help_priority=3)
@click.option('-p', '--pairwise', 'pairwise_file', callback=path_to_absolute_path, required=True, type=click.Path(exists=True), help="the pairwise TSV file")
@click.option('-g', '--groups-file', "groups_file", callback=path_to_absolute_path, required=False, default="NA", type=click.Path(exists=False), help="single-column supergroups file")
@click.option('--clusters-file', "clusters_file", callback=path_to_absolute_path, required=False, default="NA", type=click.Path(exists=False), help="DBRetina clusters file")
@click.option('--cluster-ids', "cluster_ids", callback=validate_numbers, required=False, default="", help="comma-separated list of cluster IDs")
@click.option('-d', '--dist-type', "distance_type", required=False, default="NA", show_default=True, type=click.STRING, help="select from ['min_cont', 'avg_cont', 'max_cont', 'ochiai', 'jaccard']")
@click.option('-c', '--cutoff', required=False, type=click.FloatRange(0, 100, clamp=False), default=0.0, show_default=True, help="filter out distances < cutoff")
@click.option('--extend', "extend", is_flag=True, default=False, show_default=True, help="include all supergroups that are linked to the given supergroups.")
@click.option('-o', '--output', "output_file", required=True, type=click.STRING, help="output file prefix")
@click.pass_context
def main(ctx, pairwise_file, groups_file, distance_type, cutoff, output_file, clusters_file, cluster_ids, extend):
    """Filter a pairwise file.


Detailed description:

    Filter a pairwise file by distance cutoff and/or a set of groups (provided as a single-column file or cluster IDs in a DBRetina cluster file).

Examples:

    1- distance cutoff only              | dbretina filter -p pairwise.tsv -d ochiai -c 60 -o filtered.tsv

    2- distance cutoff and groups file   | dbretina filter -p pairwise.tsv -d min_cont -c 97 -g groups.tsv -o filtered.tsv

    3- distance cutoff and a cluster IDs | dbretina filter -p pairwise.tsv -d max_cont -c 77 --clusters-file clusters.tsv --clusters-id 8 -o filtered.tsv

    4- groups file only                  | dbretina filter -p pairwise.tsv -g groups.tsv -o filtered.tsv

    5- cluster file with cluster IDs     | dbretina filter -p pairwise.tsv --clusters-file clusters.tsv --clusters-id 8 -o filtered.tsv 
    """
    
    # Extend must be used only when clusters file or groups file is provided
    if extend and groups_file == "NA" and clusters_file == "NA":
        ctx.obj.ERROR("DBRetina's filter command requires a groups_file or clusters_file if --extend is provided.")
        

    # check if not any option is provided for filteration
    if distance_type == "NA" and cutoff == 0.0 and groups_file == "NA" and clusters_file == "NA":
        ctx.obj.ERROR(
            "DBRetina's filter command requires at least one option to filter the pairwise file.")

    # if clusters_file then must be cluster_id
    if clusters_file != "NA" and not len(cluster_ids):
        ctx.obj.ERROR(
            "DBRetina's filter command requires cluster_id(s) if clusters_file is provided.")
    elif clusters_file == "NA" and len(cluster_ids):
        ctx.obj.ERROR(
            "DBRetina's filter command requires a clusters_file if cluster_id(s) is provided.")

    # can't filter by groups_file and clusters_file at the same time
    if groups_file != "NA" and clusters_file != "NA":
        ctx.obj.ERROR(
            "DBRetina's filter command can't filter by groups_file and clusters_file at the same time.")

    # if distance_type is provided then cutoff must be provided
    if distance_type != "NA" and cutoff == 0.0:
        ctx.obj.ERROR(
            "DBRetina's filter command requires a cutoff if distance_type is provided.")
    elif distance_type == "NA" and cutoff != 0.0:
        ctx.obj.ERROR(
            "DBRetina's filter command requires a distance_type if cutoff is provided.")

    if not is_awk_available():
        ctx.obj.ERROR(
            "DBRetina's filter command requires awk to be installed and available in the PATH.")

    distance_to_col = {
        "min_cont": 5,
        "avg_cont": 6,
        "max_cont": 7,
        "ochiai": 8,
        "jaccard": 9,
    }
    if distance_type in distance_to_col:
        # +1 because awk is 1-indexed
        awk_column = distance_to_col[distance_type] + 1
    elif distance_type != "NA":
        ctx.obj.ERROR(f"DBRetina's filter command doesn't support the distance_type {distance_type}.")

    # check if output_file already exist
    output_file += ".tsv"
    if os.path.exists(output_file):
        ctx.obj.WARNING(f"Output file {output_file} already exists, overwriting ...")

    if groups_file != "NA" and not os.path.exists(groups_file):
        ctx.obj.ERROR(f"Groups file {groups_file} doesn't exist.")

    comment_lines = 0
    with (open(pairwise_file) as f, open(output_file, 'w') as w):
        for line in f:
            comment_lines += 1
            if line.startswith("#"):
                w.write(line)
            else:
                w.write(f"#command: {get_command()}\n")
                w.write(line)
                break

    ctx.obj.INFO(
        f"Filtering the pairwise matrix on the {distance_type} column with a cutoff of {cutoff} and groups file {groups_file}."
    )

    _tmp_file = ".DBRetina.tmp.group"

    all_ids = set()
    if clusters_file != "NA":
        groups_file = _tmp_file
        with open(clusters_file) as f, open(_tmp_file, 'w') as W:

            # skip comments
            while True:
                pos = f.tell()
                line = f.readline()
                if not line.startswith('#'):
                    f.seek(pos)
                    break

            next(f)
            for line in f:
                line = line.strip().split('\t')
                cluster_id = int(line[0])
                if cluster_id in cluster_ids:
                    all_ids.add(cluster_id)
                    W.write(line[2].replace('|', '\n') + '\n')

    unfound_ids = set(cluster_ids).difference(all_ids)
    if len(unfound_ids):
        ctx.obj.WARNING(
            f"Couldn't find the following cluster IDs: {unfound_ids}")

    
    extended_ids_list = ".DBRetina_extended_ids_list"
    
    with open(extended_ids_list, 'w') as f:
        f.write("")
    
    if extend:
        # awk_script = f"""grep '^[^#;]' {pairwise_file} | tail -n+2 | LC_ALL=C awk -F'\t' 'BEGIN {{ while ( getline < "{groups_file}" ) {{ gsub(/"/, "", $1); id_map[tolower($1)]=1 }} }} {{ if ( ($3 in id_map) || ($4 in id_map) ) {{ print $3 >> "{extended_ids_list}"; print $4 >> "{extended_ids_list}"}} }}'"""
        
        awk_script = f"""grep '^[^#;]' {pairwise_file} | tail -n+2 | LC_ALL=C awk -F'\t' 'BEGIN {{ while ( getline < "{groups_file}" ) {{ gsub(/"/, "", $1); id_map[tolower($1)]=1 }} }} {{ if ( (tolower($3) in id_map) || (tolower($4) in id_map) ) {{ print $0 }} }}' | awk -F'\t' '{{if (${awk_column} >= {cutoff}) {{ print $3 >> "{extended_ids_list}"; print $4 >> "{extended_ids_list}"}}}}'"""
        
        result = execute_bash_command(awk_script)
        bash_script = f"""sort -u {extended_ids_list} -o {extended_ids_list}.uniq"""
        result = execute_bash_command(bash_script)
        groups_file = extended_ids_list + '.uniq'

   
    # filter by both cutoff and groups
    if cutoff != 0.0 and groups_file != "NA":
        awk_script = f"""grep '^[^#;]' {pairwise_file} | tail -n+2 | LC_ALL=C awk -F'\t' 'BEGIN {{ while ( getline < "{groups_file}" ) {{ gsub(/"/, "", $1); id_map[tolower($1)]=1 }} }} {{ if ( (tolower($3) in id_map) && (tolower($4) in id_map) ) {{ print $0 }} }}' | awk -F'\t' '{{if (${awk_column} >= {cutoff}) print $0}}' >> {output_file}"""
        result = execute_bash_command(awk_script)

    elif cutoff != 0.0:
        ctx.obj.INFO(
            f"Filtering the pairwise matrix on the {distance_type} column with a cutoff of {cutoff}.")
        command = f"grep '^[^#;]' {pairwise_file} | tail -n+2 | LC_ALL=C awk -F'\t' '{{if (${awk_column} >= {cutoff}) print $0}}' >> {output_file}"
        result = execute_bash_command(command)

    elif groups_file != "NA":
        ctx.obj.INFO(f"Filtering by groups file {groups_file}\nPlease wait...")
        awk_script = f"""grep '^[^#;]' {pairwise_file} | tail -n+2 |LC_ALL=C awk -F'\t' 'BEGIN {{ while ( getline < "{groups_file}" ) {{ gsub(/"/, "", $1); id_map[tolower($1)]=1 }} }} {{ if ( (tolower($3) in id_map) && (tolower($4) in id_map) ) {{ print $0 }} }}' >> {output_file}"""
        result = execute_bash_command(awk_script)


    # if _tmp_file exists, remove it
    if os.path.exists(_tmp_file):
        os.remove(_tmp_file)

    ctx.obj.SUCCESS("Done.")
