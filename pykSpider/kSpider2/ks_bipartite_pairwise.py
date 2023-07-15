#!/usr/bin/python
# -*- coding: utf-8 -*-
from __future__ import division
import sys
import _kSpider_internal as kSpider_internal
import click
import contextlib
from kSpider2.click_context import cli
import subprocess
import os
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import networkx as nx
import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output
import plotly.graph_objects as go
import plotly.express as px
import kSpider2.dbretina_doc_url as dbretina_doc
from sklearn.preprocessing import MinMaxScaler
from plotly.graph_objects import Figure, Parcats
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import plotly.subplots as sp
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

# TODO: Remove later
def working2_plot_bipartite(df_bipartite, output_prefix):
    from sklearn.preprocessing import LabelEncoder
    le = LabelEncoder()
    df_bipartite['group_1'] = le.fit_transform(df_bipartite['group_1'])
    fig = px.parallel_categories(df_bipartite[['group_1', 'group_2']], color="group_1", color_continuous_scale=px.colors.sequential.Inferno)
    fig.write_html(output_prefix + ".html")


def plot_bipartite(df_bipartite, color_metric, output_prefix):
    # Min-Max normalization for the color metric to get values in range [0,1]
    scaler = MinMaxScaler()
    df_bipartite['color'] = scaler.fit_transform(df_bipartite[[color_metric]])
    
    # Create a new column holding the original color metric data as categorical values
    df_bipartite[color_metric+'_cat'] = df_bipartite[color_metric].apply(lambda x: f"{color_metric}: {x:.3f}")

    fig = Figure(data=[Parcats(
        dimensions=[{'label': 'Group 1', 'values': df_bipartite['group_1']},
                    {'label': 'Group 2', 'values': df_bipartite['group_2']},
                    {'label': color_metric, 'values': df_bipartite[color_metric+'_cat'], 'visible': False}],
        line={'color': df_bipartite['color'],
              'colorscale': 'Jet', 'cmin': 0, 'cmax': 1,
              'colorbar': {'title': color_metric, 'thickness': 10, 'orientation': 'h'}
              },
        labelfont={'size': 18, 'family': 'Times'},
        tickfont={'size': 16, 'family': 'Times'},
        arrangement='freeform'
    )])
    
    fig.update_layout(coloraxis_colorbar=dict(orientation="v"))
    fig.write_html(f"{output_prefix}.html")
    fig.write_image(f"{output_prefix}.png")


# TODO: Remove later
def working_plot_bipartite(df_bipartite, output_prefix):
    B = nx.Graph()
    B.add_nodes_from(df_bipartite['group_1'], bipartite=0)
    B.add_nodes_from(df_bipartite['group_2'], bipartite=1)
    B.add_edges_from([(row['group_1'], row['group_2']) for _, row in df_bipartite.iterrows()])

    # separate by group
    l, r = nx.bipartite.sets(B)
    pos = {}

    # update position for node from each group
    pos |= ((node, (1, index)) for index, node in enumerate(l))
    pos.update((node, (2, index)) for index, node in enumerate(r))

    # Create node trace
    node_trace = go.Scatter(
        x=[pos[i][0] for i in B.nodes()],
        y=[pos[i][1] for i in B.nodes()],
        mode='markers',
        text=[f'Name: {str(n)}<br># of connections: {len(list(B.neighbors(n)))}' for n in B.nodes()],
        hoverinfo='text',
        marker=dict(
            showscale=True,
            colorscale='Viridis',
            reversescale=True,
            color=[len(list(B.neighbors(n))) for n in B.nodes()],
            size=10,
            colorbar=dict(
                thickness=15,
                title='Node Connections',
                xanchor='left',
                titleside='right'
            ),
            line=dict(width=2)))

    # Create edge trace
    edge_trace = go.Scatter(
        x=[],
        y=[],
        line=dict(width=0.5),
        hoverinfo='none',
        mode='lines')

    for edge in B.edges():
        x0, y0 = pos[edge[0]]
        x1, y1 = pos[edge[1]]
        edge_trace['x'] += (x0, x1, None)
        edge_trace['y'] += (y0, y1, None)
        # edge color based on the number of connections
        edge_trace['line']['color'] = "blue"

    # Create figure
    fig = go.Figure(data=[edge_trace, node_trace],
                    layout=go.Layout(
                        title='Network graph',
                        titlefont=dict(size=16),
                        showlegend=False,
                        hovermode='closest',
                        margin=dict(b=20, l=5, r=5, t=40),
                        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False)))

    fig.write_html(output_prefix + ".html")
    # export high quality image
    fig.write_image(output_prefix + ".png", width=1920, height=1080)


# TODO: Experminetal
def interactive_dashboard(df_bipartite):
    # create dash application
    app = dash.Dash(__name__)

    # get column names of the dataframe
    col_options = [dict(label=x, value=x) for x in df_bipartite.columns]

    # layout of the dashboard
    app.layout = html.Div([
        html.H1('Bipartite pairwise study', style={'textAlign': 'center', 'color': '#7FDBFF'}),


        html.Div(id='main-div', children=[
            html.Div([
                html.Div([
                    dcc.Dropdown(
                        id='xaxis-column',
                        options=col_options,
                        value='containment'
                    )], style={'width': '48%', 'display': 'inline-block'}),
                
                html.Div([
                    dcc.Dropdown(
                        id='yaxis-column',
                        options=col_options,
                        value='ochiai'
                    )], style={'width': '48%', 'float': 'right', 'display': 'inline-block'})
            ]),

            dcc.Dropdown(id='plot_type', options=[
                {'label': 'Scatter', 'value': 'scatter'},
                {'label': 'Bar', 'value': 'bar'},
                {'label': 'Box', 'value': 'box'},
                {'label': 'Heatmap', 'value': 'heatmap'},
                {'label': 'Parallel Categories', 'value': 'parcats'},
                {'label': 'Bipartite Network', 'value': 'bipartite'}],
                value='scatter'
            ),

            dcc.RangeSlider(
                id='pvalue-slider',
                min=df_bipartite['pvalue'].min(),
                max=df_bipartite['pvalue'].max(),
                value=[df_bipartite['pvalue'].min(), df_bipartite['pvalue'].max()],
                marks={str(pvalue): str(pvalue) for pvalue in df_bipartite['pvalue'].unique()},
                step=None
            ),

            # dcc.Graph(id='indicator-graphic')
            dcc.Graph(id='indicator-graphic', style={'height': '70vh'}, responsive=True)

        ])
    ])


    # callback to update graph when dropdown value is selected
    @app.callback(
        Output('indicator-graphic', 'figure'),
        [Input('xaxis-column', 'value'),
        Input('yaxis-column', 'value'),
        Input('plot_type', 'value'),
        Input('pvalue-slider', 'value')])
    
    def update_graph(xaxis_column_name, yaxis_column_name, plot_type, pvalue_range):
        dff = df_bipartite[(df_bipartite['pvalue'] >= pvalue_range[0]) & (df_bipartite['pvalue'] <= pvalue_range[1])]

        if plot_type == 'scatter':
            fig = px.scatter(dff, x=xaxis_column_name, y=yaxis_column_name)
        elif plot_type == 'bar':
            fig = px.bar(dff, x=xaxis_column_name, y=yaxis_column_name)
        elif plot_type == 'box':
            fig = px.box(dff, x=xaxis_column_name, y=yaxis_column_name)
        elif plot_type == 'heatmap':
            pivot_table = dff.pivot(index='group_1', columns='group_2', values='pvalue')
            fig = px.imshow(pivot_table)
        elif plot_type == 'parcats':
            fig = go.Figure(data=go.Parcats(
                dimensions=[
                    {'label': 'Group 1', 'values': dff['group_1']},
                    {'label': 'Group 2', 'values': dff['group_2']},
                    {'label': 'P-value', 'values': dff['pvalue']}]
            ))
        elif plot_type == 'bipartite':
            # Here you need to include the code to create bipartite graph. 
            # Creating a bipartite graph in Plotly can be complex, so you need to decide how to implement it according to your needs.
            fig = go.Figure()

        fig.update_layout(margin={'l': 40, 'b': 40, 't': 10, 'r': 0}, hovermode='closest')
        return fig
    
    app.run_server(debug=True)

def check_if_there_is_a_pvalue(pairwise_file):
    with open(pairwise_file) as F:
        for line in F:
            if not line.startswith("#"):
                return "pvalue" in line
            else:
                continue


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




@cli.command(name="bipartite", epilog = dbretina_doc.doc_url("bipartite"), help_priority=8)
@click.option('-p', '--pairwise', 'pairwise_file', callback=path_to_absolute_path, required=True, type=click.Path(exists=True), help="the pairwise TSV file")
@click.option('--group1', "group_1_file", callback=path_to_absolute_path, required=False, type=click.Path(exists=True), help="group1 single-column supergroups file")
@click.option('--group2', "group_2_file", callback=path_to_absolute_path, required=False, type=click.Path(exists=True), help="group2 single-column supergroups file")
@click.option('--gmt1', "gmt_1_file", callback=path_to_absolute_path, required=False, type=click.Path(exists=True), help="GMT file 1")
@click.option('--gmt2', "gmt_2_file", callback=path_to_absolute_path, required=False, type=click.Path(exists=True), help="GMT file 2")
# @click.option('-m', '--metric', "metric", required=True, type=click.STRING, help="select from ['containment', 'ochiai', 'jaccard', 'pvalue']")
@click.option('-o', '--output', "output_prefix", required=True, type=click.STRING, help="output file prefix")
@click.pass_context
def main(ctx, pairwise_file, group_1_file, group_2_file, gmt_1_file, gmt_2_file, output_prefix):
    """
        Create a bipartite connections between two group files.
    """
    LOGGER = ctx.obj
    
    FILE_HAS_PVALUE = check_if_there_is_a_pvalue(pairwise_file)
                    
    # must be two group files or two gmt files
    if (not gmt_1_file and not gmt_2_file) and (not group_1_file and not group_2_file):
        LOGGER.ERROR("Please provide either two GMT files or two group files.")
            

    ###########################################################
    # 1. parse the two group files to dictionary for O(1) access
    ########################################################### 

    LOGGER.INFO("Parsing the two group files...")

    group1_dict = {}
    group2_dict = {}
    unmatched_groups = []
    
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
        LOGGER.ERROR("There is an overlap between the two groups. Please make sure there is no overlap between the two groups.")

    ###########################################################
    # 2. parse the pairwise file and create the bipartite graph
    ###########################################################

    LOGGER.INFO("Parsing the pairwise file...")

    metric_to_col = {
        "containment": 5,
        "ochiai": 6,
        "jaccard": 7,
        "odds_ratio": 8,
        "pvalue": 9,
    }

    df_bipartite = pd.DataFrame(columns=["group_1", "group_2", "containment", "ochiai", "jaccard"])

    df_rows = []

    if check_if_there_is_a_pvalue(pairwise_file):
        df_bipartite["pvalue"] = None

    metadata = []
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
                

    ###########################################################
    # 3. create the bipartite graph
    ###########################################################    

    one_to_one = df_bipartite.drop_duplicates(subset=['group_1', 'group_2'])
    one_to_many = df_bipartite.groupby('group_1').filter(lambda x: len(x['group_2']) > 1)
    many_to_one = df_bipartite.groupby('group_2').filter(lambda x: len(x['group_1']) > 1)
    statistics = df_bipartite.describe()


    ###########################################################
    # draw the graph
    ## TODO: Nice scatter plot, but let's do it later 
    """ # << ---- ### DISABLED FOR NOW ### ---- 
    LOGGER.INFO(f"Writing the scatter graph to {output_prefix}_scatter.html")
    fig = px.scatter(df_bipartite, x="containment", y="ochiai", color="pvalue", size='jaccard', hover_data=['group_1','group_2'])
    fig.write_html(f"{output_prefix}_scatter.html")
    """
    
    ###########################################################

    LOGGER.INFO(f"Writing the bipartite graph to {output_prefix}_bipartite.html")
    LOGGER.INFO(f"Writing the bipartite graph to {output_prefix}_bipartite.png")
    plot_bipartite(df_bipartite, 'containment', f"{output_prefix}_bipartite")

    ###########################################################

    # Create a heatmap
    ######### HEATMAP | DSABLED FOR NOW #######################
    """ <<< ----- DISABLED FOR NOW -----
    LOGGER.INFO(f"Writing the heatmap to {output_prefix}_heatmap.html")
    pivot_table = df_bipartite.pivot(index='group_1', columns='group_2', values='pvalue')
    fig = px.imshow(pivot_table)

    # add labels
    fig.update_layout(
        title="Pivot table of p-values",
        xaxis_title="Group 2",
        yaxis_title="Group 1",
        font=dict(
            family="Courier New, monospace",
            size=18,
            color="#7f7f7f"
        )
    )
    # save to disk
    fig.write_html(f"{output_prefix}_pivot_heatmap.html")
    --------------------------------- >>> """    

    # TODO: DISABLED for now. No need to plot 3 categories 
    """
    ###########################################################
    ########### PARALLEL CATEGORIES PLOT ######################
    ########################################################### 
    LOGGER.INFO(f"Writing the parcats graph to {output_prefix}_parcats.html")
    df_bipartite['color_group'] = pd.Categorical(df_bipartite['group_1']).codes
    scaler = MinMaxScaler()
    df_bipartite['color'] = scaler.fit_transform(df_bipartite[[metric]])
    fig = go.Figure(data=
        go.Parcats(
            dimensions=[
                {'label': 'Group 1', 'values': df_bipartite['group_1']},
                {'label': 'Group 2', 'values': df_bipartite['group_2']},
                {'label': 'P-value', 'values': df_bipartite['pvalue']}],
            line={
                'color': df_bipartite['color_group'], 
                'colorscale': 'Jet',
                'colorbar': {'title': metric, 'thickness': 10, 'orientation': 'h'}
                },
            hoveron='color', 
            hoverinfo='count+probability',
            )
        )

    fig.write_html(f"{output_prefix}_parcats.html")
    fig.write_image(f"{output_prefix}_parcats.png", width=1920, height=1080, scale=2)
    """

    ###########################################################

    # TODO: [DEV] BETA
    # interactive_dashboard(df_bipartite)

    
    ## TODO: RESEARCH CODE STUDY PVALUE-CUTOFF #1
    #################################################################
    # 4. study the different statistics with different pvalue cutoffs
    #################################################################

    """ << --------- RESEARCH CODE ---------
    def calculate_connections(df, cutoff):
        # Filter the dataframe by the cutoff
        df_filtered = df[df['pvalue'] <= cutoff]

        # Count the number of unique group2 values for each group1 value
        counts_1 = df_filtered.groupby('group_1')['group_2'].nunique()
        # Count the number of unique group1 values for each group2 value
        counts_2 = df_filtered.groupby('group_2')['group_1'].nunique()

        # Calculate the number of 1-to-1, group1-to-many, and group2-to-many connections
        one_to_one = sum(counts_1 == 1) + sum(counts_2 == 1)
        group1_to_many = sum(counts_1 > 1)
        group2_to_many = sum(counts_2 > 1)

        # Total 1-to-many is the sum of group1-to-many and group2-to-many
        total_one_to_many = group1_to_many + group2_to_many

        return one_to_one, group1_to_many, group2_to_many, total_one_to_many


    LOGGER.INFO("Calculating the connections for different pvalue cutoffs (Please Wait!)")

    cutoffs = np.arange(df_bipartite['pvalue'].min(), df_bipartite['pvalue'].max(), 0.001)
    results_cutoffs = []
    for cutoff in cutoffs:
        one_to_one, group1_to_many, group2_to_many, total_one_to_many = calculate_connections(df_bipartite, cutoff)

        # Append the results for this cutoff to the dataframe
        results_cutoffs.append({
            'pvalue_cutoff': cutoff,
            '1-to-1': one_to_one,
            'group1-to-many': group1_to_many,
            'group2-to-many': group2_to_many,
            'total_1-to-many': total_one_to_many
        })

    results = pd.DataFrame(results_cutoffs)

    LOGGER.INFO(f"Writing the results to {output_prefix}_pvalues_cutoffs.tsv")
    results.to_csv(f"{output_prefix}_pvalues_cutoffs.tsv", sep='\t', index=False)

    #################################################################
    # 4.1 Plotting the analysis results for cutoff vs. connections

    #### 4.1.1 Correlation Heatmap of connections
    LOGGER.INFO(f"Writing the correlation heatmap to {output_prefix}_pvalues_correlation_heatmap.png")
    plt.figure(figsize=(10, 8))
    sns.set(style="white")
    # Calculate correlation matrix
    corr = results.iloc[:, 1:].corr()
    # Generate a mask for the upper triangle
    mask = np.triu(np.ones_like(corr, dtype=float))
    # Draw the heatmap
    sns.heatmap(corr, mask=mask, annot=True, fmt=".002f", linewidths=.5, cmap='coolwarm')
    plt.title('Correlation Heatmap of connections')
    plt.savefig(f"{output_prefix}_pvalues_correlation_heatmap.png")

    #### 4.1.2 Pairplot of connections    
    LOGGER.INFO(f"Writing the pairplot to {output_prefix}_pvalues_pairplot.png")
    sns.set(style="ticks", color_codes=True)
    # Exclude 'pvalue_cutoff' from the pairplot
    pairplot_data = results.iloc[:, :] # all rows, all columns except the first column
    # Draw the pairplot
    sns.pairplot(pairplot_data)
    plt.savefig(f"{output_prefix}_pvalues_pairplot.png", dpi=400)


    #### 4.1.3 Scatterplot of connections
    LOGGER.INFO(f"Writing the scatterplot to {output_prefix}_pvalues_cutoffs_scatterplot.html")
    fig = sp.make_subplots(rows=4, cols=1)
    # Create scatter plots for each metric against pvalue_cutoff
    for i, metric in enumerate(results.columns[1:], start=1):
        fig.add_trace(
            go.Scatter(x=results['pvalue_cutoff'], y=results[metric], mode='lines+markers', name=metric),
            row=i, col=1
        )

    # Update layout
    fig.update_layout(height=800, width=800, title_text="Metrics vs P-Value Cutoff")
    fig.write_html(f"{output_prefix}_pvalues_cutoffs_scatterplot.html")


    #### 4.1.4 Parallel plot
    LOGGER.INFO(f"Writing the parallel plot to {output_prefix}_pvalues_cutoffs_parallelplot.html")
    # Create a parallel coordinates plot
    fig = px.parallel_coordinates(results, color="pvalue_cutoff", 
                                labels={"1-to-1": "One-to-One", 
                                        "group1-to-many": "Group1-to-Many", 
                                        "group2-to-many": "Group2-to-Many", 
                                        "total_1-to-many": "Total One-to-Many"},
                                color_continuous_scale=px.colors.diverging.Tealrose,
                                color_continuous_midpoint=results['pvalue_cutoff'].median())

    fig.update_layout(
        title='Parallel Coordinates Plot for Different connections',
        autosize=True,
        # width=1000,
        # height=600,
    )

    fig.write_html(f"{output_prefix}_pvalues_cutoffs_parallelplot.html")
    
    #################################################################

    # Interactive Dashboard
    # interactive_dashboard(df_bipartite)
    
    >> --------- RESEARCH CODE --------- """
