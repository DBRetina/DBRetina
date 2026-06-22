#!/usr/bin/python
# -*- coding: utf-8 -*-
from __future__ import division
import sys
import _dbretina_internal as dbretina_internal
import click
import contextlib
from dbretina.click_context import cli
from dbretina.validators import validate_metric
import subprocess
import os
import pandas as pd
import plotly.express as px
import dbretina.dbretina_doc_url as dbretina_doc
from plotly.graph_objects import Figure, Parcats
import seaborn as sns
import matplotlib.pyplot as plt
import json

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


def path_to_absolute_path(ctx, param, value):
    with contextlib.suppress(Exception):
        return os.path.abspath(value) if value is not None else None

def plot_bipartite(df_bipartite, color_metric, output_prefix):
    df_bipartite['color'] = df_bipartite[color_metric]
    lowest_color = df_bipartite['color'].min()
    largest_color = df_bipartite['color'].max()
    
    # Create a new column holding the original color metric data as categorical values
    df_bipartite[color_metric+'_cat'] = df_bipartite[color_metric].apply(lambda x: f"{color_metric}: {x:.3f}")

    fig = Figure(data=[Parcats(
        dimensions=[{'label': 'Group 1', 'values': df_bipartite['group_1']},
                    {'label': 'Group 2', 'values': df_bipartite['group_2']},
                    {'label': color_metric, 'values': df_bipartite[color_metric+'_cat'], 'visible': False}],
        line={'color': df_bipartite['color'],
              'colorscale': 'Solar', 'cmin': lowest_color, 'cmax': largest_color,
              'colorbar': {'title': color_metric, 'thickness': 10, 'orientation': 'h'}
              },
        labelfont={'size': 18, 'family': 'Times'},
        tickfont={'size': 16, 'family': 'Times'},
        arrangement='freeform'
    )])
    
    # Compute the number of unique categories to set the height
    n_categories = len(pd.concat([df_bipartite['group_1'], df_bipartite['group_2']]).unique())
    height_per_category = 20  # adjust this as necessary
    minimum_height = 1300
    fig_height = max(minimum_height, n_categories * height_per_category)
    fig.update_layout(height=fig_height)
    
    
    
    # Compute margin size based on longest label
    labels = pd.concat([df_bipartite['group_1'], df_bipartite['group_2']])
    left_margin = df_bipartite['group_1'].str.len().max() * 4
    right_margin = df_bipartite['group_2'].str.len().max() * 4
    
    fig.update_layout(margin={"r":right_margin,"l":left_margin, 'pad' : 5})

    
    fig.update_layout(coloraxis_colorbar=dict(orientation="v"))
    fig.write_html(f"{output_prefix}.html")
    fig.write_image(f"{output_prefix}.png", width=1920, height=1080, scale=5)


def similarities_distribution_histogram(df_bipartite, filename, json_output_file = None, log_scale = False):
    # Function to map values to ranges
    def map_value_to_range(value):
        lower = (int(value) // 5) * 5
        upper = lower + 5
        if upper > 100: 
            upper = 100
        return f"{lower}-{upper}"

    # List of metrics to consider
    metrics = ["containment", "ochiai", "jaccard"]

    # Initialize a dict to store data
    data = {}

    # Create all possible ranges with zero counts
    all_ranges = [f"{i}-{i+5}" for i in range(0, 105, 5)]
    for metric in metrics:
        data[metric] = dict.fromkeys(all_ranges, 0)

    # Iterate over metrics
    for metric in metrics:
        # Bin the similarity scores into ranges
        df_bipartite[metric+'_range'] = df_bipartite[metric].apply(map_value_to_range)
        
        # Count the number of pairs in each range
        counts = df_bipartite[metric+'_range'].value_counts().to_dict()

        # Update counts in data
        data[metric].update(counts)

    # Prepare data for plotting
    df = pd.DataFrame(data)

    # Convert index to integer for sorting
    df.index = df.index.map(lambda x: int(x.split('-')[0]))
    df = df.sort_index()

    # Reset index to string for correct x-axis labels
    df.index = df.index.map(lambda x: f"{x}-{x+5}")

    # Melt dataframe to have format suitable for seaborn
    df = df.reset_index().melt('index', var_name='Metric', value_name='Count')
    
    df['index'] = df['index'].replace('100-105', '100-100')
    
    if json_output_file:
        # Dictionary to store the results
        results = {}

        # Iterate over the metrics
        for metric in ['containment', 'ochiai', 'jaccard']:
            metric_df = df[df['Metric'] == metric]  # Filter the dataframe by the current metric
            metric_dict = dict(zip(metric_df['index'], metric_df['Count']))  # Create a dictionary: {range: count}
            results[metric] = metric_dict  # Add the dictionary to the results

        # Write the results to a JSON file
        with open(json_output_file, 'w') as f:
            json.dump(results, f, indent=4)


    # Create the plot
    plt.figure(figsize=(10, 6))
    sns.set(style="whitegrid")
    sns.set_color_codes("pastel")
    sns.barplot(x='index', y='Count', hue='Metric', data=df, palette="viridis")


    # Add title and labels
    plt.title('Histogram of Similarity Scores', fontsize=18)
    plt.xlabel('Similarity Score Range', fontsize=14)
    plt.ylabel('Count', fontsize=14)

    # Increase the size of the legend and x-axis ticks labels
    plt.legend(title='Metrics', title_fontsize='13', fontsize='12')
    plt.legend(loc='upper right')
    plt.xticks(fontsize=10, rotation=90)



    if log_scale:
        plt.ylabel('Count (log scale)')
    else:
        plt.ylabel('Count')

    # Rotate x-axis labels for better readability
    plt.xticks(rotation=90)
    
    # y-axis in log scale
    if log_scale:
        plt.yscale('log')

    # Show the plot
    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    # plt.show()


def plot_pivot_table(df_bipartite, metric, output_prefix):
    pivot_table = df_bipartite.pivot(index='group_1', columns='group_2', values=metric)
    fig = px.imshow(pivot_table)

    # add labels
    fig.update_layout(
        title="Pivot table of group1 vs group2",
        xaxis_title="Group 2",
        yaxis_title="Group 1",
        font=dict(
            family="Times",
            size=16,
            color="#7f7f7f"
        )
    )
    
    # Dynamic plot size based on pivot table dimensions
    cell_size = 50  # adjust this value as necessary
    width = pivot_table.shape[1] * cell_size
    height = pivot_table.shape[0] * cell_size
    fig.update_layout(width=width, height=height)

    # make hover label = metric
    fig.update_traces(hovertemplate = f"{metric}: %{{z:.2f}}")

    fig.write_html(f"{output_prefix}.html")
    fig.write_image(f"{output_prefix}.png")


@cli.command(name="bipartite", epilog = dbretina_doc.doc_url("bipartite"), help_priority=8)
@click.option('-p', '--pairwise', 'pairwise_file', callback=path_to_absolute_path, required=True, type=click.Path(exists=True), help="the pairwise TSV file")
@click.option('--group1', "group_1_file", callback=path_to_absolute_path, required=False, type=click.Path(exists=True), help="group1 single-column supergroups file")
@click.option('--group2', "group_2_file", callback=path_to_absolute_path, required=False, type=click.Path(exists=True), help="group2 single-column supergroups file")
@click.option('--gmt1', "gmt_1_file", callback=path_to_absolute_path, required=False, type=click.Path(exists=True), help="GMT file 1")
@click.option('--gmt2', "gmt_2_file", callback=path_to_absolute_path, required=False, type=click.Path(exists=True), help="GMT file 2")
@click.option('-m', '--metric', "metric", required=True, type=click.STRING, callback=validate_metric, help="Bipartite coloring based on ['containment', 'ochiai', 'jaccard', 'pvalue']")
@click.option('-c', '--cutoff', 'cutoff', required=False, type=click.FloatRange(0, 100, clamp=False), default=0.0, show_default = True, help="Include comparisons (similarity > cutoff)")
@click.option('--no-plot', "no_plot", is_flag=True, default=False, help="do not plot the bipartite graph")
@click.option('--no-1-1', "no_1_1", is_flag=True, default=False, help="do not include 1-1 mapping")
@click.option('-o', '--output', "output_prefix", required=True, type=click.STRING, help="output file prefix")
@click.pass_context
def main(ctx, pairwise_file, group_1_file, group_2_file, gmt_1_file, gmt_2_file, output_prefix, metric, no_plot, cutoff, no_1_1):
    """
        Create a bipartite connections between two group files.
    """
    LOGGER = ctx.obj

    # check if pvalue (format-aware: handles the .tsv / parquet-dir / .dbrp
    # forms; the old text open() crashed with IsADirectoryError on a parquet
    # directory, issue 053). Resolve once and reuse for the schema below.
    from dbretina.compat import pairwise_has_pvalue
    input_has_pvalue = pairwise_has_pvalue(pairwise_file)
    if metric == "pvalue" and not input_has_pvalue:
        LOGGER.ERROR("pvalue not found in pairwise file!")


    ###########################################################
    # 1. parse the two group files to dictionary for O(1) access
    ########################################################### 

    LOGGER.INFO("Parsing the two group files...")

    group1_dict = {}
    group2_dict = {}
    deleted_groups = []

    dbretina_str_escape = lambda x: x.lower().replace('"', '')

    # if gmt files are provided, convert them to group files
    if gmt_1_file and gmt_2_file:
        with open(gmt_1_file) as IN_GMT:
            for line in IN_GMT:
                group1_dict[dbretina_str_escape(line.strip().split("\t")[0])] = {}
        with open(gmt_2_file) as IN_GMT:
            for line in IN_GMT:
                group2_dict[dbretina_str_escape(line.strip().split("\t")[0])] = {}
    elif group_1_file and group_2_file:
        with open(group_1_file) as IN_GROUP:
            for line in IN_GROUP:
                group1_dict[dbretina_str_escape(line.strip())] = {}

        with open(group_2_file) as IN_GROUP:
            for line in IN_GROUP:
                group2_dict[dbretina_str_escape(line.strip())] = {}

    else:
        LOGGER.ERROR("Please provide either two GMT files or two group files.")


    # make sure there is no overlap between the two groups
    if set(group1_dict.keys()).intersection(set(group2_dict.keys())):
        # delete the overlapping groups
        LOGGER.WARNING("There is an overlap between the two groups. The overlapping groups will be removed from the second groups file.")
        for group in set(group1_dict.keys()).intersection(set(group2_dict.keys())):
            # del group1_dict[group]
            del group2_dict[group]
            deleted_groups.append(group)

    if len(deleted_groups):
        LOGGER.WARNING(
            f'The following groups were removed from the second groups file: {", ".join(deleted_groups)}'
        )
                    

    ###########################################################
    # 2. parse the pairwise file and create the bipartite graph
    ###########################################################

    LOGGER.INFO("Parsing the pairwise file...")

    metric_to_col = {
        "containment": 5,
        "ochiai": 6,
        "jaccard": 7,
        "csi": 8,
        "dice": 9,
        "odds_ratio": 10,
        "pvalue": 11,
    }

    # Guard the lookup: the -m callback validates real metrics but lets the "NA"
    # sentinel through, and a bare metric_to_col[...]/store.to_pandas(metric=...)
    # would surface a raw KeyError/ValueError traceback. Mirror the graph command.
    if metric not in metric_to_col:
        LOGGER.ERROR(f"Invalid metric '{metric}'. Choose from: {', '.join(metric_to_col)}")
        sys.exit(1)

    df_bipartite = pd.DataFrame(columns=["group_1", "group_2", "containment", "ochiai", "jaccard"])

    df_rows = []

    if input_has_pvalue:
        df_bipartite["pvalue"] = None

    # Try Parquet/PairwiseStore first
    try:
        from dbretina.compat import open_pairwise
        store = open_pairwise(pairwise_file)
    except Exception:
        store = None

    if store is not None:
        LOGGER.INFO("using Parquet pairwise data via PairwiseStore")
        metadata = [f"#command: {get_command()}\n"]
        names_map = store.get_names_map()
        cols = ["group_1_id", "group_2_id", "containment", "ochiai", "jaccard"]
        if store.has_pvalue:
            cols.append("pvalue")
            df_bipartite["pvalue"] = None
        pdf = store.to_pandas(metric=metric, cutoff=cutoff, columns=cols)
        for _, row in pdf.iterrows():
            _source_1 = names_map.get(int(row["group_1_id"]), "")
            _source_2 = names_map.get(int(row["group_2_id"]), "")

            if _source_1 in group1_dict and _source_2 in group2_dict:
                group1 = _source_1
                group2 = _source_2
            elif _source_1 in group2_dict and _source_2 in group1_dict:
                group1 = _source_2
                group2 = _source_1
            else:
                continue

            df_row = {
                "group_1": group1,
                "group_2": group2,
                "containment": float(row["containment"]),
                "ochiai": float(row["ochiai"]),
                "jaccard": float(row["jaccard"]),
            }
            if store.has_pvalue:
                df_row["pvalue"] = float(row["pvalue"])
            df_rows.append(df_row)
        store.close()
    else:
        metadata = []
        dbrp_path = pairwise_file.replace(".tsv", ".dbrp")
        metric_name_to_id = {
            "containment": 0, "ochiai": 1, "jaccard": 2, "csi": 3,
            "dice": 4, "odds_ratio": 5, "pvalue": 6,
        }

        if os.path.exists(dbrp_path):
            LOGGER.INFO(f"found .dbrp file: {dbrp_path}, using binary pairwise reader")
            mid = metric_name_to_id.get(metric, 0)
            records = dbretina_internal.dbrp_filter_pairs(dbrp_path, mid, cutoff)
            has_pvalue = len(records) > 0 and 'pvalue' in records[0]
            if has_pvalue:
                df_bipartite["pvalue"] = None
            metadata.append(f"#command: {get_command()}\n")
            for rec in records:
                _source_1 = rec['group_1_name']
                _source_2 = rec['group_2_name']

                if _source_1 in group1_dict and _source_2 in group2_dict:
                    group1 = _source_1
                    group2 = _source_2
                elif _source_1 in group2_dict and _source_2 in group1_dict:
                    group1 = _source_2
                    group2 = _source_1
                else:
                    continue

                df_row = {
                    "group_1": group1,
                    "group_2": group2,
                    "containment": float(rec['containment']),
                    "ochiai": float(rec['ochiai']),
                    "jaccard": float(rec['jaccard']),
                }
                if has_pvalue:
                    df_row["pvalue"] = float(rec['pvalue'])
                df_rows.append(df_row)
        else:
            with open(pairwise_file, 'r') as pairwise_tsv:
                while True:
                    pos = pairwise_tsv.tell()
                    line = pairwise_tsv.readline()
                    if not line.startswith('#'):
                        pairwise_tsv.seek(pos)
                        break
                    else:
                        metadata.append(line)
                metadata.append(f"#command: {get_command()}\n")

                next(pairwise_tsv)
                for row in pairwise_tsv:
                    row = row.strip().split('\t')

                    _source_1 = row[2]
                    _source_2 = row[3]

                    if float(row[metric_to_col[metric]]) < cutoff:
                        continue

                    if _source_1 in group1_dict and _source_2 in group2_dict:
                        group1 = _source_1
                        group2 = _source_2
                    elif _source_1 in group2_dict and _source_2 in group1_dict:
                        group1 = _source_2
                        group2 = _source_1
                    else:
                        continue

                    if "pvalue" in df_bipartite.columns:
                        df_row = {
                            "group_1": group1,
                            "group_2": group2,
                            "containment": float(row[metric_to_col["containment"]]),
                            "ochiai": float(row[metric_to_col["ochiai"]]),
                            "jaccard": float(row[metric_to_col["jaccard"]]),
                            "pvalue": float(row[metric_to_col["pvalue"]]),
                        }
                    else:
                        df_row = {
                            "group_1": group1,
                            "group_2": group2,
                            "containment": float(row[metric_to_col["containment"]]),
                            "ochiai": float(row[metric_to_col["ochiai"]]),
                            "jaccard": float(row[metric_to_col["jaccard"]]),
                        }

                    df_rows.append(df_row)

    LOGGER.INFO(f"Writing the bipartite TSV file to {output_prefix}_bipartite_pairwise.tsv")
    df_bipartite = pd.DataFrame(df_rows)
    
    # check if there is no rows 
    if df_bipartite.shape[0] == 0:
        LOGGER.ERROR("There is no overlap between the two groups. Please check the input files.")
    
    if no_1_1:
        # remove one-to-one overlaps between group_1 and group_2
        LOGGER.INFO("Removing one-to-one overlaps between group_1 and group_2")    # Count the number of connections for each member in both groups
        group_1_counts = df_bipartite['group_1'].value_counts()
        group_2_counts = df_bipartite['group_2'].value_counts()
        
        # Filter the DataFrame to only keep rows where the count is more than 1 for both groups
        df_bipartite = df_bipartite[(group_1_counts[df_bipartite['group_1']].values > 1) & 
                                    (group_2_counts[df_bipartite['group_2']].values > 1)]
        
        
    df_bipartite.to_csv(f"{output_prefix}_bipartite_pairwise.tsv", sep='\t', index=False)

    histogram_plot_file = f"{output_prefix}_similarity_metrics_histogram.png"
    LOGGER.INFO(f"Plotting the similarity metrics histogram to {histogram_plot_file}")
    similarities_distribution_histogram(df_bipartite, histogram_plot_file, log_scale = False)

    histogram_plot_file = f"{output_prefix}_similarity_metrics_histogram_log.png"
    json_stats_file = f"{output_prefix}_similarity_metrics_histogram.json"
    LOGGER.INFO(f"Writing the similarity metrics histogram to {json_stats_file}")
    LOGGER.INFO(f"Plotting the similarity metrics histogram (log-scale) to {histogram_plot_file}")
    similarities_distribution_histogram(df_bipartite, histogram_plot_file, json_output_file=json_stats_file, log_scale = True)

    # report if there are unmatched groups
    unique_matched_group1 = set(df_bipartite['group_1'].unique())
    unique_matched_group2 = set(df_bipartite['group_2'].unique())
    unfound_group1 = set(group1_dict.keys()).difference(unique_matched_group1)
    unfound_group2 = set(group2_dict.keys()).difference(unique_matched_group2)
    if unfound_group1 or unfound_group2:
        LOGGER.WARNING(f"Missing gene set names detected, reporting to {output_prefix}_missing_groups.txt")
        with open(f"{output_prefix}_missing_groups.txt", 'w') as OUT:
            if unfound_group1:
                OUT.write(f"Missing group1 names: {','.join(list(unfound_group1))}\n")
            if unfound_group2:
                OUT.write(f"Missing group2 names: {','.join(list(unfound_group2))}\n")


    if not no_plot:
        # Static image (PNG) export needs the optional 'kaleido' package; if plotting fails
        # (kaleido missing, or empty/degenerate data) warn and continue — the bipartite data
        # outputs are already written above, so a plotting failure must not crash the command.
        try:
            LOGGER.INFO(f"Writing the bipartite graph to {output_prefix}_bipartite.html / .png")
            plot_bipartite(df_bipartite, metric, f"{output_prefix}_bipartite")
            LOGGER.INFO(f"Writing the pivot table to {output_prefix}_pivot_table.html / .png")
            plot_pivot_table(df_bipartite, metric, f"{output_prefix}_pivot_table")
        except Exception as e:
            LOGGER.WARNING(
                f"plotting skipped (image export failed: {e}). "
                f"Install 'kaleido' for PNG export, or pass --no-plot to skip plotting."
            )

    LOGGER.SUCCESS("Done!")
