#!/usr/bin/python
# -*- coding: utf-8 -*-
from __future__ import division
import sys
import _kSpider_internal as kSpider_internal
import click
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
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import plotly.subplots as sp
import kSpider2.dbretina_doc_url as dbretina_doc

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
    return value if value == "NA" else os.path.abspath(value)

def check_if_there_is_a_pvalue(pairwise_file):
    with open(pairwise_file) as F:
        for line in F:
            if not line.startswith("#"):
                return "pvalue" in line
            else:
                continue


def plot_bipartite(df_bipartite, output_prefix):
    B = nx.Graph()
    B.add_nodes_from(df_bipartite['group_1'], bipartite=0)
    B.add_nodes_from(df_bipartite['group_2'], bipartite=1)
    edges = [(row['group_1'], row['group_2'], row['pvalue']) for _, row in df_bipartite.iterrows()]
    B.add_weighted_edges_from(edges)

    # separate by group
    l, r = nx.bipartite.sets(B)
    pos = {}

    # update position for node from each group
    pos |= ((node, (1, index)) for index, node in enumerate(l))
    pos.update((node, (2, index)) for index, node in enumerate(r))

    # add edges as disconnected lines in a single trace
    edge_trace = go.Scatter(
        x=[],
        y=[],
        line=dict(width=0.5, color='#888'),
        hoverinfo='none',
        mode='lines')

    for edge in B.edges(data=True):
        x0, y0 = pos[edge[0]]
        x1, y1 = pos[edge[1]]
        edge_trace['x'] += (x0, x1, None)
        edge_trace['y'] += (y0, y1, None)

    # add nodes as scatter trace
    node_trace = go.Scatter(
        x=[],
        y=[],
        text=[],
        mode='markers',
        hoverinfo='text',
        marker=dict(
            showscale=True,
            colorscale='YlGnBu',
            reversescale=True,
            color=[],
            size=10,
            colorbar=dict(
                thickness=15,
                title='Node Connections',
                xanchor='left',
                titleside='right'
            ),
            line=dict(width=2)))

    for node in B.nodes():
        x, y = pos[node]
        node_trace['x'] += (x, )
        node_trace['y'] += (y, )

    # color node points by the number of connections
    for adjacencies in B.adjacency():
        node_trace['marker']['color'] += (len(adjacencies[1]), )
        node_info = f'Name: {str(adjacencies[0])}<br># of connections: {len(adjacencies[1])}'
        node_trace['text'] += (node_info, )

    # create plotly figure and add traces
    fig = go.Figure(data=[edge_trace, node_trace],
                    layout=go.Layout(
                        title='<br>Network graph',
                        titlefont=dict(size=16),
                        showlegend=False,
                        hovermode='closest',
                        margin=dict(b=20, l=5, r=5, t=40),
                        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False)))

    fig.write_html(output_prefix)


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


@cli.command(name="bipartite", help_priority=8)
@click.option('-p', '--pairwise', 'pairwise_file', callback=path_to_absolute_path, required=True, type=click.Path(exists=True), help="the pairwise TSV file")
@click.option('--group1', "group_1_file", callback=path_to_absolute_path, required=True, type=click.Path(exists=True), help="group1 single-column supergroups file")
@click.option('--group2', "group_2_file", callback=path_to_absolute_path, required=True, type=click.Path(exists=True), help="group2 single-column supergroups file")
@click.option('-o', '--output', "output_prefix", required=True, type=click.STRING, help="output file prefix")
@click.pass_context
def main(ctx, pairwise_file, group_1_file, group_2_file, output_prefix):
    """
        Create a bipartite relationships between two groups file
    """
    LOGGER = ctx.obj

    ###########################################################
    # 1. parse the two group files to dictionary for O(1) access
    ########################################################### 

    LOGGER.INFO("Parsing the two group files...")

    group1_dict = {}
    group2_dict = {}
    unmatched_groups = []

    with open(group_1_file) as IN_GROUP:
        for line in IN_GROUP:
            group1_dict[line.strip()] = {}

    with open(group_2_file) as IN_GROUP:
        for line in IN_GROUP:
            group2_dict[line.strip()] = {}

    # make sure there is no overlap between the two groups
    if set(group1_dict.keys()).intersection(set(group2_dict.keys())):
        LOGGER.ERROR("There is an overlap between the two groups. Please make sure there is no overlap between the two groups.")

    ###########################################################
    # 2. parse the pairwise file and create the bipartite graph
    ###########################################################

    LOGGER.INFO("Parsing the pairwise file...")

    distance_to_col = {
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
                    "containment": float(row[distance_to_col["containment"]]),
                    "ochiai": float(row[distance_to_col["ochiai"]]),
                    "jaccard": float(row[distance_to_col["jaccard"]]),
                    "pvalue": float(row[distance_to_col["pvalue"]]),
                }
            else:
                df_row = {
                    "group_1": group1,
                    "group_2": group2,
                    "containment": float(row[distance_to_col["containment"]]),
                    "ochiai": float(row[distance_to_col["ochiai"]]),
                    "jaccard": float(row[distance_to_col["jaccard"]]),
                }

            df_rows.append(df_row)

    LOGGER.INFO(f"Writing the bipartite TSV file to {output_prefix}_bipartite_full_relationships.tsv")
    df_bipartite = pd.DataFrame(df_rows)
    df_bipartite.to_csv(f"{output_prefix}_bipartite_full_relationships.tsv", sep='\t', index=False)
    
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
    LOGGER.INFO(f"Writing the scatter graph to {output_prefix}_bipartite.html")
    fig = px.scatter(df_bipartite, x="containment", y="ochiai", color="pvalue",
                    size='jaccard', hover_data=['group_1','group_2'])
    fig.write_html(f"{output_prefix}.html")
    ###########################################################

    LOGGER.INFO(f"Writing the bipartite graph to {output_prefix}_bipartite.html")
    plot_bipartite(df_bipartite, f"{output_prefix}_bipartite.html")

    ###########################################################

    # Create a heatmap
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

    ###########################################################

    LOGGER.INFO(f"Writing the parcats graph to {output_prefix}_parcats.html")
    fig = go.Figure(data=
    go.Parcats(
        dimensions=[
            {'label': 'Group 1',
             'values': df_bipartite['group_1']},
            {'label': 'Group 2',
             'values': df_bipartite['group_2']},
            {'label': 'P-value',
             'values': df_bipartite['pvalue']}]
        )
    )
    fig.write_html(f"{output_prefix}_parcats.html")

    ###########################################################

    # TODO: BETA
    # interactive_dashboard(df_bipartite)


    #################################################################
    # 4. study the different statistics with different pvalue cutoffs
    #################################################################


    def calculate_relationships(df, cutoff):
        # Filter the dataframe by the cutoff
        df_filtered = df[df['pvalue'] <= cutoff]

        # Count the number of unique group2 values for each group1 value
        counts_1 = df_filtered.groupby('group_1')['group_2'].nunique()
        # Count the number of unique group1 values for each group2 value
        counts_2 = df_filtered.groupby('group_2')['group_1'].nunique()

        # Calculate the number of 1-to-1, group1-to-many, and group2-to-many relationships
        one_to_one = sum(counts_1 == 1) + sum(counts_2 == 1)
        group1_to_many = sum(counts_1 > 1)
        group2_to_many = sum(counts_2 > 1)

        # Total 1-to-many is the sum of group1-to-many and group2-to-many
        total_one_to_many = group1_to_many + group2_to_many

        return one_to_one, group1_to_many, group2_to_many, total_one_to_many


    LOGGER.INFO("Calculating the relationships for different pvalue cutoffs (Please Wait!)")

    cutoffs = np.arange(df_bipartite['pvalue'].min(), df_bipartite['pvalue'].max(), 0.001)
    results_cutoffs = []
    for cutoff in cutoffs:
        one_to_one, group1_to_many, group2_to_many, total_one_to_many = calculate_relationships(df_bipartite, cutoff)

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
    # 4.1 Plotting the analysis results for cutoff vs. relationships

    #### 4.1.1 Correlation Heatmap of Relationships
    LOGGER.INFO(f"Writing the correlation heatmap to {output_prefix}_pvalues_correlation_heatmap.png")
    plt.figure(figsize=(10, 8))
    sns.set(style="white")
    # Calculate correlation matrix
    corr = results.iloc[:, 1:].corr()
    # Generate a mask for the upper triangle
    mask = np.triu(np.ones_like(corr, dtype=float))
    # Draw the heatmap
    sns.heatmap(corr, mask=mask, annot=True, fmt=".002f", linewidths=.5, cmap='coolwarm')
    plt.title('Correlation Heatmap of Relationships')
    plt.savefig(f"{output_prefix}_pvalues_correlation_heatmap.png")

    #### 4.1.2 Pairplot of Relationships    
    LOGGER.INFO(f"Writing the pairplot to {output_prefix}_pvalues_pairplot.png")
    sns.set(style="ticks", color_codes=True)
    # Exclude 'pvalue_cutoff' from the pairplot
    pairplot_data = results.iloc[:, :] # all rows, all columns except the first column
    # Draw the pairplot
    sns.pairplot(pairplot_data)
    plt.savefig(f"{output_prefix}_pvalues_pairplot.png", dpi=400)


    #### 4.1.3 Scatterplot of Relationships
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
        title='Parallel Coordinates Plot for Different Relationships',
        autosize=True,
        # width=1000,
        # height=600,
    )

    fig.write_html(f"{output_prefix}_pvalues_cutoffs_parallelplot.html")
    
    #################################################################

    # Interactive Dashboard
    # interactive_dashboard(df_bipartite)