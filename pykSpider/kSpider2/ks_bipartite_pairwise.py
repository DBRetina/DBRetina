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


def path_to_absolute_path(ctx, param, value):
    return value if value == "NA" else os.path.abspath(value)

def check_if_there_is_a_pvalue(pairwise_file):
    with open(pairwise_file) as F:
        for line in F:
            if not line.startswith("#"):
                return "pvalue" in line
            else:
                continue


def plot_bipartite(df_bipartite, output_file):
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

    fig.write_html(output_file)


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
@click.option('-o', '--output', "output_file", required=True, type=click.STRING, help="output file prefix")
@click.pass_context
def main(ctx, pairwise_file, group_1_file, group_2_file, output_file):
    """
        Create a bipartite relationships between two groups file
    """
    LOGGER = ctx.obj

    ###########################################################
    # 1. parse the two group files to dictionary for O(1) access
    ###########################################################    

    group1_dict = {}
    group2_dict = {}

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


    df_bipartite = pd.DataFrame(df_rows)
    df_bipartite.to_csv(output_file, sep='\t', index=False)

    ###########################################################
    # 3. create the bipartite graph
    ###########################################################    

    one_to_one = df_bipartite.drop_duplicates(subset=['group_1', 'group_2'])
    one_to_many = df_bipartite.groupby('group_1').filter(lambda x: len(x['group_2']) > 1)
    many_to_one = df_bipartite.groupby('group_2').filter(lambda x: len(x['group_1']) > 1)
    statistics = df_bipartite.describe()


    ###########################################################
    # draw the graph
    fig = px.scatter(df_bipartite, x="containment", y="ochiai", color="pvalue",
                    size='jaccard', hover_data=['group_1','group_2'])
    fig.write_html(f"{output_file}.html")
    ###########################################################

    plot_bipartite(df_bipartite, f"{output_file}_bipartite.html")
    
    ###########################################################
    
    # Create a heatmap
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
    fig.write_html(f"{output_file}_pivot_heatmap.html")
    
    ###########################################################
    
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
    fig.write_html(f"{output_file}_parcats.html")
    
    ###########################################################

    interactive_dashboard(df_bipartite)

    

        
    
    
        