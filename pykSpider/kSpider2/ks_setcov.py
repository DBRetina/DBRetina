#!/usr/bin/env python
import pandas as pd
import numpy as np
import csv
from collections import defaultdict
import argparse
import matplotlib.pyplot as plt
import networkx as nx
import math
from click.decorators import option
import click
from kSpider2.click_context import cli
import kSpider2.dbretina_doc_url as dbretina_doc
import json
import os
from kSpider2.ks_clustering import main as ks_clustering


class Graph:
    def __init__(self):
        self.parent = {}
        self.components = None

    def find(self, x):
        if x not in self.parent:
            self.parent[x] = x
        if self.parent[x] != x:
            self.parent[x] = self.find(self.parent[x])
        return self.parent[x]

    def union(self, x, y):
        self.parent[self.find(x)] = self.find(y)

    def add_edge(self, x, y):
        self.union(x, y)

    def _compute_components(self):
        self.components = {}
        for node in self.parent:
            parent = self.find(node)
            if parent not in self.components:
                self.components[parent] = []
            self.components[parent].append(node)

    def get_connected_components(self):
        if self.components is None:
            self._compute_components()
        return self.components.values()


class DeduplicateGroups():
    """
    A class for deduplicating groups.

    This class provides methods for deduplicating groups based on various criteria, such as item coverage, GPI, and CSI.
    """
    
    communities_clusters_file = ''
    associations_file = ''
    cluster_id_to_groups = dict()
    groups_to_items_no = dict()
    total_number_of_groups = 0
    main_pairwise_file = ''
    index_prefix = ''
    exact_ochiai_cutoff = 80
    ochiai_community_cutoff = 30
    GC = 100
    final_remaining_groups_count = 0
    final_remaining_groups = set()
    removed_exact_ochiai_groups = set()
    
    # Dataframes
    df_items_metadata = None
    df_groups_metadata = None
    df_associations = None
    df_item_to_CSI = None
    df_group_to_avg_CSI = None
    df_group_to_avg_gpi = None
    df_item_to_gpi = None
    df_group_to_modularity = None
    df_logging = None
    df_groups_per_item = None
    
    LOGGER = None # click custom logger
    ctx = None # click context
    
    metric_to_col = {
        "containment": 5,
        "ochiai": 6,
        "jaccard": 7,
        "odds_ratio": 8,
        "pvalue": 9,
    }
    
    def __init__(
        self, 
        associations_file, 
        index_prefix,
        GC = 100,
        containment_cutoff = 80,
        exact_ochiai_cutoff = 80,
        ochiai_community_cutoff = 30,
        ctx = None
        ):
        self.associations_file = associations_file
        self.ochiai_community_cutoff = ochiai_community_cutoff
        # Initialize a logging DataFrame
        self.df_logging = pd.DataFrame(columns=[
            'group', 'no_of_items', 'coverage %', 'gpi', 'CSI', 'fragmentation', 'heteroitemity', 'modularity', 'selected', 'cluster_id'
        ])
        self.main_pairwise_file = index_prefix + "_DBRetina_pairwise.tsv"
        self.GC = GC
        self.containment_cutoff = containment_cutoff
        self.index_prefix = index_prefix
        self.exact_ochiai_cutoff = exact_ochiai_cutoff
        self.LOGGER = ctx.obj
        self.ctx = ctx
        

    # returns a dictionary cluster_id -> groups
    def cluster_to_groups(self, clusters_file):
        cluster_id_to_groups = {}
        header_flag = True
        with open(clusters_file, 'r') as f:
            for line in f:
                if line.startswith("#"): continue
                if header_flag:
                    header_flag = False
                    continue
                cluster_id, cluster_size, cluster_members = line.strip().split("\t")
                cluster_id_to_groups[int(cluster_id)] = cluster_members.split("|")
        return cluster_id_to_groups


    def process_associations(self):
        self.df_associations = pd.read_csv(self.associations_file, sep="\t")
        # rename 'hgnc_symbol_ids to 'item'
        # rename second column to item
        self.df_associations.rename(columns={self.df_associations.columns[1]: 'item'}, inplace=True)
        self.df_associations.rename(columns={self.df_associations.columns[0]: 'group'}, inplace=True)
        # all to lower
        self.df_associations['item'] = self.df_associations['item'].str.lower()
        self.df_associations['group'] = self.df_associations['group'].str.lower()
        # calculate group to size
        self.groups_to_items_no = self.df_associations.groupby('group')['item'].nunique().to_dict()
        # constrcut dictionary of group to items list
        self.groups_to_items = self.df_associations.groupby('group')['item'].apply(set).to_dict()


    def build_item_to_CSI(self):
        # we have can use df_item_to_gpi to get the number of groups and clusters per item
        
        total_number_of_clusters = self.cluster_id_to_groups.keys().__len__()
        self.total_number_of_groups = self.df_associations['group'].nunique()
        self.df_groups_per_item = self.df_associations.groupby('item')['group'].nunique()
        
        self.df_item_to_CSI = pd.DataFrame({
            'item' : self.df_item_to_gpi['item'].values,
            'group_count': self.df_item_to_gpi['group_count'].values,
            'cluster_count': self.df_item_to_gpi['cluster_count'].values,
            'CSI': 100 * 
                    np.log2(self.df_item_to_gpi['cluster_count'].values / total_number_of_clusters) / 
                    np.log2(1.0 / total_number_of_clusters)
        })
        
    
    def build_group_to_CSI(self):
        #1 merge our original associations with item_to_CSI
        merged_df = pd.merge(self.df_associations, self.df_item_to_CSI, left_on='item', right_on='item', how='left')
        #2 Calculate the sum of CSI for each group and the number of items per group
        group_group = merged_df.groupby('group')['CSI'].agg(['sum', 'count'])
        #3 Calculate average CSI
        self.df_group_to_avg_CSI = pd.DataFrame({
            'average_CSI': group_group['sum'] / group_group['count']
        })
        
    def build_group_to_gpi(self):
        #1 Get the number of unique clusters and groups each item is found in
        self.cluster_id_to_groups = self.cluster_to_groups(self.communities_clusters_file)
        # Convert dictionary to dataframe
        cluster_to_groups_df = pd.DataFrame([(k, v) for k, vv in self.cluster_id_to_groups.items() for v in vv], 
                                            columns=['cluster_id', 'group'])

        #2 map cluster_id to items
        cluster_to_items = pd.merge(cluster_to_groups_df, self.df_associations, on='group', how='left')
        cluster_to_items = pd.merge(cluster_to_groups_df, self.df_associations, on='group', how='left')


        #3. Get the number of unique clusters and groups each item is found in
        item_counts = cluster_to_items.groupby('item').agg(
            cluster_count=pd.NamedAgg(column='cluster_id', aggfunc='nunique'),
            group_count=pd.NamedAgg(column='group', aggfunc='nunique')
        )
        
        #4. Calculate the Group Pleiotropy Index (GPI) for each item (Cd/Nd)
        item_counts['gpi'] = item_counts['cluster_count'] / item_counts['group_count']

        #5 All information in one dataframe
        self.df_item_to_gpi = pd.DataFrame({
            'item': item_counts.index,
            'cluster_count': item_counts['cluster_count'].values,
            'group_count': item_counts['group_count'].values,
            'gpi': item_counts['gpi'].values
        })

        #5. Merge cluster_to_items with item_to_gpi to get gpi for each item in each group
        # ['cluster_id', 'group', 'item', 'cluster_count', 'group_count', 'gpi']
        merged_df = pd.merge(cluster_to_items, self.df_item_to_gpi, on='item', how='inner')

        # Calculate the average GPI for each group
        group_group = merged_df.groupby('group')['gpi'].mean()

        # Create the DataFrame
        self.df_group_to_avg_gpi = pd.DataFrame({
            'group': group_group.index,
            'average_gpi': 100 * group_group.values
        })
             
        
    def export_item_to_gpi_CSI(self, file_name):
        # Merge item_to_gpi with item_to_CSI
        merged_df = pd.merge(self.df_item_to_CSI, self.df_item_to_gpi, on='item', how='left')
        # Export to CSV
        merged_df.to_csv(file_name, index=False, sep='\t')
        return file_name
    
    def export_group_to_gpi_CSI(self, file_name):
        # Merge group_to_gpi with group_to_CSI
        merged_df = pd.merge(self.df_group_to_avg_CSI, self.df_group_to_avg_gpi, on='group', how='left')
        # Export to CSV
        merged_df.to_csv(file_name, index=False, sep='\t')
        return file_name
        
    
    def process_pairwise_file(self, tsv_file, max_cont_threshold):
        # Initialize a dictionary to store group metrics
        group_metrics = defaultdict(lambda: {'fragmentation': 0, 'heteroitemity': 0})
        
        self.ochiai_graph = Graph()

        with open(tsv_file, 'r') as f:
            # Skip comment lines
            while True:
                pos = f.tell()  # remember the position
                line = f.readline()
                if not line.startswith("#"):
                    f.seek(pos)  # rewind to the position before the line
                    break

            # skip the header
            next(f)

            reader = csv.reader(f, delimiter="\t")
            for row in reader:
                # Extract the group names from columns 3 and 4
                group1 = row[2]
                group2 = row[3]
                
                containment = float(row[self.metric_to_col["containment"]])
                ochiai_similarity = float(row[self.metric_to_col["ochiai"]])
                
                # instead of re-iterating over the file, do it here
                if ochiai_similarity >= self.exact_ochiai_cutoff:
                    self.ochiai_graph.add_edge(group1, group2)
                
                # this for calculating modularity
                if containment < max_cont_threshold:
                    continue      

                # Get the lengths of the groups
                group1_len = self.groups_to_items_no[group1]
                group2_len = self.groups_to_items_no[group2]

                # Compute the metrics
                # Large node fragmented to small nodes
                # Large node heterogenous to small nodes
                if group1_len < group2_len:
                    group_metrics[group1]['fragmentation'] -= 1
                    group_metrics[group2]['heteroitemity'] += 1
                elif group1_len > group2_len:
                    group_metrics[group1]['heteroitemity'] += 1
                    group_metrics[group2]['fragmentation'] -= 1
                # else:  # if they are equal, update both
                #     print("ELSEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE")
                #     group_metrics[group1]['fragmentation'] -= 1
                #     group_metrics[group1]['heteroitemity'] += 1
                #     group_metrics[group2]['fragmentation'] -= 1
                #     group_metrics[group2]['heteroitemity'] += 1

        # Compute the modularity index for each group
        for metrics in group_metrics.values():
            fragmentation = metrics['fragmentation']
            heteroitemity = metrics['heteroitemity']
            metrics['modularity'] = abs(fragmentation + heteroitemity)

        # Convert the dictionary to a pandas DataFrame
        self.df_group_to_modularity = pd.DataFrame.from_dict(group_metrics, orient='index').reset_index()
        self.df_group_to_modularity.columns = ['group', 'fragmentation', 'heteroitemity', 'modularity']
        
    def build_groups_metadata(self):
        # This builds the groups metadata (gpi, CSI, frag, het, modularity, length)
        #1 merge df_group_to_avg_gpi and df_group_to_CSI
        self.df_groups_metadata = pd.merge(self.df_group_to_avg_gpi, self.df_group_to_avg_CSI, on='group', how='outer')
        # add deduplication status, default = remain

        #2 merge with df_group_indices
        self.df_groups_metadata = pd.merge(self.df_groups_metadata, self.df_group_to_modularity, on='group', how='outer')

        #3 Add no_of_items column from groups_to_items_no dictionary
        self.df_groups_metadata['no_of_items'] = self.df_groups_metadata['group'].apply(lambda x: self.groups_to_items_no[x])

        #4 Rename the columns
        self.df_groups_metadata['dedup'] = 'remained'
        self.df_groups_metadata.columns = ['group', 'average_gpi', 'average_CSI', 'fragmentation', 'heteroitemity', 'modularity', 'no_of_items', 'dedup']


    def remove_exact_ochiai_matches(self):
        connected_components = self.ochiai_graph.get_connected_components()
        unselected_groups = set()
        selected_groups = set()
        
        # create a temporary column for number of edges
        self.df_groups_metadata['total_edges'] = abs(self.df_groups_metadata["heteroitemity"]) + abs(self.df_groups_metadata["fragmentation"])

        for component_groups in connected_components:
            if len(component_groups) == 1:
                selected_groups.add(component_groups[0])
                continue
            
            # sort the groups by lowest average_CSI and highest no_of_items
            # modified to reflect `DBRetina dedup`
            sorted_groups = self.df_groups_metadata[self.df_groups_metadata['group'].isin(component_groups)].sort_values(by=['total_edges', 'no_of_items'], ascending=[False, False])
            selected_groups.add(sorted_groups.iloc[0]['group'])
            unselected_groups.update(sorted_groups.iloc[1:]['group'].tolist())
                
        # drop the temporary column
        self.df_groups_metadata.drop('total_edges', axis=1, inplace=True)
        
        # set deduplication status to exact_ochiai
        self.removed_exact_ochiai_groups = unselected_groups
        self.filtered_exact_ochiai_groups_count = len(unselected_groups)
        
        
        # set df_groups_metadata dedup to exact_ochiai if in self.removed_exact_ochiai_groups
        self.df_groups_metadata.loc[self.df_groups_metadata['group'].isin(self.removed_exact_ochiai_groups), 'dedup'] = 'exact_ochiai'

    
    def export_groups_metadata(self, file_name):
        # Save the dataframe to a CSV file
        # fill NA values with 0
        self.df_groups_metadata = self.df_groups_metadata.fillna(0)
        self.df_groups_metadata.to_csv(file_name, index=False, sep='\t', columns=['group', 'no_of_items', 'average_gpi', 'average_CSI', 'fragmentation', 'heteroitemity', 'modularity', 'dedup'])
        return file_name

    def cluster_to_universe_set(self, cluster_id):
        # Get groups for the given cluster ID
        cluster_groups = self.cluster_id_to_groups.get(cluster_id, [])

        # Filter df_associations to only include rows with groups in cluster_groups
        cluster_df = self.df_associations[self.df_associations['group'].isin(cluster_groups)]

        # Get the unique set of items associated with these groups
        return set(cluster_df['item'].unique())

    def unique_set_of_items(self):
        return set(self.df_associations['item'].unique())
    
        

    def modularity_based_set_cover(self, GC = 100):
        
        def _deduplicate_all_groups_instead_of_communities(GC = 100):
            selected_groups = {}
            uncovered_items = self.unique_set_of_items()
            total_items = len(uncovered_items)
            items_covered = 0
            sorted_groups_df = self.df_groups_metadata[
                self.df_groups_metadata['dedup'] != 'exact_ochiai'].sort_values(
                    by=['modularity', 'average_CSI', 'no_of_items'], ascending=[True, False, False]
                    )
            
            # Iterate through sorted groups
            for _, group_row in sorted_groups_df.iterrows():
                group = group_row['group']
                # Calculate the intersection of the current group items with the uncovered items
                group_items = self.groups_to_items[group]
                if common_items := group_items.intersection(uncovered_items):
                    # Add the group to the set cover
                    coverage = len(common_items) / total_items * 100
                    selected_groups[group] = coverage
                    # Remove the covered items from the uncovered set
                    uncovered_items -= common_items
                    # Update the count of covered items
                    items_covered += len(common_items)
                    # If the required item coverage is met, break the loop
                    if items_covered / total_items * 100 >= GC:
                        break

            return selected_groups
        
        selected_groups = list(_deduplicate_all_groups_instead_of_communities(GC).keys())
        self.final_remaining_groups_count = len(selected_groups)
        self.final_remaining_groups.update(selected_groups)
        
        self.df_groups_metadata.loc[
            ~self.df_groups_metadata['group'].isin(selected_groups) &
            ~self.df_groups_metadata['group'].isin(self.removed_exact_ochiai_groups), 'dedup'] = 'set-cov'

    def write_new_associations_file(self, new_associations_file):
        with open(new_associations_file, 'w') as NEW_ASSOCIATIONS_FILE:
            NEW_ASSOCIATIONS_FILE.write('group\titem\n')
            for group in self.final_remaining_groups:
                for item in self.groups_to_items[group]:
                    NEW_ASSOCIATIONS_FILE.write(f'{group}\t{item}\n')    


    def export_split_group_metadata(self, file_name):
        remaining_groups_df = self.df_groups_metadata[self.df_groups_metadata['group'].isin(self.final_remaining_groups)]
        removed_groups_df = self.df_groups_metadata[~self.df_groups_metadata['group'].isin(self.final_remaining_groups)]
        remaining_groups_df.to_csv(f'{file_name}_remaining_groups_metadata.tsv', index=False, sep='\t', columns=['group', 'no_of_items', 'average_gpi', 'average_CSI', 'fragmentation', 'heteroitemity', 'modularity'])
        removed_groups_df.to_csv(f'{file_name}_removed_groups_metadata.tsv', index=False, sep='\t', columns=['group', 'no_of_items', 'average_gpi', 'average_CSI', 'fragmentation', 'heteroitemity', 'modularity'])


    def export_deduplicated_gmt(self, file_name):
        file_name = f'{file_name}'
        groups_of_interest_df = self.df_associations[
            self.df_associations['group'].isin(self.final_remaining_groups)
        ]
        groups_to_items_df = groups_of_interest_df.groupby('group')['item'].apply(list).reset_index()
        with open(file_name, 'w') as f:
            fake_url = 'PLACEHOLDER_DESCRIPTION'
            for idx, row in groups_to_items_df.iterrows():
                row['item'] = [item.upper() for item in row['item']]
                f.write('\t'.join([row['group'], fake_url] + row['item']))
                # f.write('\t'.join([row['group'], row['group']] + row['item']))
                f.write('\n')

    def export_original_gmt(self, file_name):
        file_name = f'{file_name}'
        groups_to_items_df = self.df_associations.groupby('group')['item'].apply(list).reset_index()
        with open(file_name, 'w') as f:
            for idx, row in groups_to_items_df.iterrows():
                # upper case row['item']
                row['item'] = [item.upper() for item in row['item']]
                f.write('\t'.join([row['group'], row['group']] + row['item']))
                f.write('\n')

    def calculate_items_stats(self):
        # calculate overlap score for the original data
        original_overlap_score = self.df_groups_per_item.mean()
        final_stats = {
            "original_overlap_score": 0,
            "remaining_overlap_score": 0,
            "original_no_of_items": 0,
            "remaining_no_of_items": 0,
            "percentage_of_items_remaining": 0,
            'original_overlap_score': original_overlap_score,
            'original_no_of_items': len(self.df_groups_per_item),
        }
        # create a new dataframe removing groups not in final_remaining_groups
        df_remaining_associations = self.df_associations[
            self.df_associations['group'].isin(self.final_remaining_groups)]

        # calculate number of remaining groups per item
        df_remaining_groups_per_item = df_remaining_associations.groupby('item')['group'].nunique()

        final_stats['remaining_overlap_score'] = df_remaining_groups_per_item.mean()
        final_stats['remaining_no_of_items'] = len(df_remaining_groups_per_item)

        final_stats['percentage_of_items_remaining'] = final_stats['remaining_no_of_items'] / final_stats['original_no_of_items'] * 100

        return final_stats
    
    def perform_cli_community_detection(self, pairwise_file, cutoff, output_prefix):
        self.LOGGER.ACTIVE = False
        self.ctx.invoke(ks_clustering, 
                pairwise_file = pairwise_file, 
                cutoff = cutoff,
                metric = 'ochiai',
                output_prefix = output_prefix, 
                community = True,
                )
        self.LOGGER.ACTIVE = True
        
        # delete  if exist
        _histo_plot = f"{output_prefix}_clusters_histogram.png"
        _bubbles_plot = f"{output_prefix}_clusters_bubbles.png"
        if os.path.exists(_histo_plot):
            os.remove(_histo_plot)
        if os.path.exists(_bubbles_plot):
            os.remove(_bubbles_plot)
        
        return f"{output_prefix}_clusters.tsv"

    def build_all(self, output_prefix):
        
        self.LOGGER.INFO(f"Detecting communities with ochiai cutoff {self.ochiai_community_cutoff}")
        _communities_file_prefix = f"{output_prefix}_communities"
        self.communities_clusters_file = self.perform_cli_community_detection(
            self.main_pairwise_file, 
            self.containment_cutoff, 
            _communities_file_prefix)
        
        self.process_associations()
        self.build_group_to_gpi()
        self.build_item_to_CSI()
        self.build_group_to_CSI()
        self.LOGGER.SUCCESS("Calculated GPI and CSI")

        self.process_pairwise_file(self.main_pairwise_file, self.containment_cutoff)
        self.LOGGER.SUCCESS(f"Group modularities computed with containment cutoff {self.containment_cutoff}")

        self.build_groups_metadata()
        self.LOGGER.SUCCESS("Groups metadata built")

        _file_item_gpi_CSI = f"{output_prefix}_item_to_gpi_CSI.tsv"
        self.export_item_to_gpi_CSI(_file_item_gpi_CSI)
        self.LOGGER.SUCCESS(f"Exported Item to GPI CSI to {_file_item_gpi_CSI}")

        self.remove_exact_ochiai_matches()
        self.LOGGER.SUCCESS("Deduplication completed")

        self.modularity_based_set_cover(self.GC)
        self.LOGGER.SUCCESS("Set coverage process completed")

        self.export_groups_metadata(f"{output_prefix}_groups_metadata.tsv")
        self.LOGGER.SUCCESS(f"Exported groups metadata to {output_prefix}_groups_metadata.tsv")

        self.export_split_group_metadata(output_prefix)
        self.LOGGER.SUCCESS(f"Remaining groups metadata exported to {output_prefix}_remaining_groups_metadata.tsv") 
        self.LOGGER.SUCCESS(f"Removed groups metadata exported to {output_prefix}_removed_groups_metadata.tsv")

        self.write_new_associations_file(f"{output_prefix}_associations.tsv")
        self.LOGGER.SUCCESS(f"The new association file exported to {output_prefix}_associations.tsv")

        DEDUP_GMT = output_prefix + "_DEDUP_GMT.gmt"
        ORIGINAL_GMT = output_prefix + "_ORIGINAL_GMT.gmt"
        self.export_deduplicated_gmt(DEDUP_GMT)
        self.export_original_gmt(ORIGINAL_GMT)
        self.LOGGER.SUCCESS(f"Original and setcov GMT files exported to {DEDUP_GMT} and {ORIGINAL_GMT}")

        print("-------------------------")      
        self.LOGGER.INFO(f"original number of groups {self.total_number_of_groups}")
        number_of_groups_removed_from_set_cover = self.total_number_of_groups - self.filtered_exact_ochiai_groups_count - self.final_remaining_groups_count
        percentage_number_of_groups_removed_from_set_cover = 100 * number_of_groups_removed_from_set_cover / self.total_number_of_groups
        number_of_groups_removed_from_exact_duplicates = self.filtered_exact_ochiai_groups_count
        percentage_number_of_groups_removed_from_exact_duplicates = 100 * number_of_groups_removed_from_exact_duplicates / self.total_number_of_groups
        self.LOGGER.INFO(f"groups removed due to exact duplication: {self.filtered_exact_ochiai_groups_count} = {percentage_number_of_groups_removed_from_exact_duplicates}%")
        self.LOGGER.INFO(f"number of groups removed from set-cover only: {number_of_groups_removed_from_set_cover} = {percentage_number_of_groups_removed_from_set_cover}%")
        self.LOGGER.INFO(f"final number of groups {self.final_remaining_groups_count} = {100 * self.final_remaining_groups_count / self.total_number_of_groups}%")
        total_percentage = percentage_number_of_groups_removed_from_exact_duplicates + percentage_number_of_groups_removed_from_set_cover + (100 * self.final_remaining_groups_count / self.total_number_of_groups)
        self.LOGGER.INFO(f"[DEBUG] total percentage of groups: {total_percentage}%")

        items_stats = self.calculate_items_stats()
        self.LOGGER.INFO(f"original number of items {items_stats['original_no_of_items']}")
        self.LOGGER.INFO(f"remaining number of items {items_stats['remaining_no_of_items']}")
        self.LOGGER.INFO(f"percentage of items remaining {items_stats['percentage_of_items_remaining']}%")
        self.LOGGER.INFO(f"original overlap score {items_stats['original_overlap_score']}")
        self.LOGGER.INFO(f"remaining overlap score {items_stats['remaining_overlap_score']}")
    

class GraphBasedDeduplication(DeduplicateGroups):
    def graph_based_setcover(self, GC = 100):
        def plot_and_highligh_specific_node(G, node, filename, node_color='r'):
            # clear plt
            plt.clf()
            # fixed layout
            nx.draw(G, G_pos)
            # draw nodes with small labels
            nx.draw_networkx_nodes(G, G_pos, nodelist=[node], node_color=node_color)
            # small node label
            nx.draw_networkx_labels(G, G_pos, font_size=5)
            # show edge weights
            # nx.draw_networkx_edge_labels(G, pos, font_size=3)
            plt.savefig(filename, dpi=400)



        selected_groups = {}
        uncovered_items = self.unique_set_of_items()
        total_items = len(uncovered_items)
        items_covered = 0
        iteration = 1        

        # get initial sorting of groups
        sorted_groups_df = self.df_groups_metadata[
            self.df_groups_metadata['dedup'] != 'exact_ochiai'].sort_values(
                by=['modularity', 'average_CSI', 'no_of_items'], ascending=[True, False, False]
                )

        # drop rows having dedup = 'exact_ochiai'
        sorted_groups_df = sorted_groups_df[sorted_groups_df['dedup'] != 'exact_ochiai']

        # add similarity_bin column to sorted_groups_df
        sorted_groups_df['similarity_bin'] = 0
        sorted_groups_df['similarity'] = 0
        
        sorted_groups_df.to_csv(f"iteration_{iteration}_loop_sorted_groups_df.tsv", sep="\t")


        # index sorted_groups_df by group
        sorted_groups_df = sorted_groups_df.set_index('group')


        graph_remaining_groups_df = pd.DataFrame(columns = sorted_groups_df.columns)
        first_group_name = sorted_groups_df.iloc[0].name
        print(f"[DEBUG] first group name: {first_group_name}")

        # load pairwise tsv file
        df_edges = pd.read_csv(self.main_pairwise_file, sep='\t', comment='#', usecols=['group_1_name', 'group_2_name', 'ochiai'])
        

        # keep edges that have group_1_name and group_2_name in sorted_groups_df
        df_edges = df_edges[df_edges['group_1_name'].isin(sorted_groups_df.index) & df_edges['group_2_name'].isin(sorted_groups_df.index)]

        # Create a graph with the edges
        G = nx.from_pandas_edgelist(df_edges, source='group_1_name', target='group_2_name', edge_attr='ochiai')
        
        # add nodes to the graph that does not have any edges
        G.add_nodes_from(sorted_groups_df.index.difference(G.nodes()))
        
        
        G_pos = nx.spring_layout(G)


        # lambdas
        get_all_neighbors = lambda x: list(G.neighbors(x))
        get_all_neighbors_edges = lambda x: [(x, y) for y in get_all_neighbors(x)]
        get_farthest_neighbor = lambda x: max(get_all_neighbors(x), key=lambda y: G[x][y]['ochiai'])
        get_average_ochiai = lambda x: np.mean([G[x][y]['ochiai'] for y in get_all_neighbors(x)])
        get_sum_ochiai = lambda x: np.sum([G[x][y]['ochiai'] for y in get_all_neighbors(x)])
        get_all_neighbor_weights_dict = lambda x: {y: G[x][y]['ochiai'] for y in get_all_neighbors(x)}


        # first update the similarity_bin of the whole table
        for neighbor, weight in get_all_neighbor_weights_dict(first_group_name).items():
            print(f"[DEBUG] neighbor: {neighbor}, weight: {weight}")
            sorted_groups_df.loc[neighbor, 'similarity'] += weight
            sorted_groups_df.loc[neighbor, 'similarity_bin'] += math.ceil(sorted_groups_df.loc[neighbor, 'similarity'] / 10)

        print(f"[DEBUG] graph_remaining_groups_df: {graph_remaining_groups_df.head()}")

        # plot_and_highligh_specific_node(G, first_group_name, f"iteration_{iteration}_first_group.png")

        # remove the first group from sorted_groups_df and from the graph
        
        # TODO REMOVE DEBUG
        # sorted_groups_df.to_csv(f"iteration_ZERO_{iteration}_sorted_groups_df.tsv", sep="\t")
        sorted_groups_df.drop(first_group_name, inplace=True)
        print(f"[REMOVED] {first_group_name} removed |  Graph nodes: {len(G.nodes())}")
        
        print(f"[DEBUG] Graph size: {len(G.nodes())}")
        G.remove_node(first_group_name)
        print(f"[DEBUG] {first_group_name} removed |  Graph size: {len(G.nodes())}")
        # remove first_group_name edges

        # update uncovered_items with groups_to_items[first_group_name]
        intersected_items = self.groups_to_items[first_group_name].intersection(uncovered_items)
        # update uncovered_items
        uncovered_items = uncovered_items.difference(intersected_items)
        items_covered += len(intersected_items)


        # add first group to graph_remaining_groups_df without the simularity_bin column, ignore the index
        self.final_remaining_groups_count += 1
        self.final_remaining_groups.add(first_group_name)
        selected_groups[first_group_name] = {'similarity_bin': len(intersected_items)}

        # if len(uncovered_items) == 0:
        #     return graph_remaining_groups_df.index.tolist()

        while True and sorted_groups_df.shape[0] > 0:
            iteration += 1
            # TODO REMOVE DEBUG
            # sorted_groups_df.to_csv(f"iteration_{iteration}_sorted_groups_df.tsv", sep="\t")

            # sort groups by lowest similarity_bin, lowest modularity, highest average_CSI, highest no_of_items
            print(f"Iteration {iteration} | number of sorted_groups_df = {len(sorted_groups_df)}")
            sorted_groups_df = sorted_groups_df.sort_values(
                by=['similarity_bin', 'modularity', 'average_CSI', 'no_of_items'], ascending=[True, True, False, False]
                )

            # get the first group
            first_group_name = sorted_groups_df.iloc[0].name
            print(f"[CURRENT GROUP] {first_group_name} |  Graph nodes: {len(G.nodes())} | Uncovered items: {len(uncovered_items)} | Items covered: {items_covered} | Iteration: {iteration} | Remaining groups: {len(sorted_groups_df)}")
            # export graph nodes to TSV
            intersected_items = self.groups_to_items[first_group_name].intersection(uncovered_items)
            print(f"[DEBUG] Next group to remove: {first_group_name}")

            group_neighbors = get_all_neighbor_weights_dict(first_group_name)
            for neighbor, weight in group_neighbors.items():
                print(f"[DEBUG] neighbor: {neighbor}, weight: {weight}")
                sorted_groups_df.loc[neighbor, 'similarity'] *= (iteration - 1)
                sorted_groups_df.loc[neighbor, 'similarity'] += weight
                sorted_groups_df.loc[neighbor, 'similarity'] /= iteration
                sorted_groups_df.loc[neighbor, 'similarity_bin'] += math.ceil(sorted_groups_df.loc[neighbor, 'similarity'] / 10)


            # update the graph_remaining_groups_df with the first group
            # add first group to graph_remaining_groups_df without the simularity_bin column, ignore the index
            self.final_remaining_groups_count += 1
            self.final_remaining_groups.add(first_group_name)

            # graph_remaining_groups_df = graph_remaining_groups_df.append(sorted_groups_df.loc[first_group_name])
            # plot_and_highligh_specific_node(G, first_group_name, f"iteration_{iteration}_first_group.png")

            # remove the first group from sorted_groups_df and from the graph
            G.remove_node(first_group_name)
            sorted_groups_df.drop(first_group_name, inplace=True)
            print(f"[REMOVED] {first_group_name} removed |  Graph size: {len(G.nodes())}")
            
            if len(intersected_items) == 0:
                print(f"!!! [DEBUG] iteration: {iteration}, first_group_name: {first_group_name}, intersected_items: {intersected_items}, uncovered_items: {len(uncovered_items)}")
                iteration -= 1
                continue
            
            # update uncovered_items
            uncovered_items = uncovered_items.difference(intersected_items)

            # update items_covered
            items_covered += len(intersected_items)
            selected_groups[first_group_name] = {'similarity_bin': len(intersected_items)}


            print(f"\n[DEBUG] items_coverage until now: {items_covered} / {total_items} ({items_covered / total_items * 100}%)")
            if items_covered / total_items * 100 >= GC:
                break


        print(f"[DEBUG] selected_groups: {selected_groups}")
        self.final_remaining_groups_count = len(selected_groups)
        self.final_remaining_groups.update(selected_groups)

        self.df_groups_metadata.loc[
            ~self.df_groups_metadata['group'].isin(selected_groups) &
            ~self.df_groups_metadata['group'].isin(self.removed_exact_ochiai_groups), 'dedup'] = 'set-cov'
        


    def build_all(self, output_prefix):
        print(f"Building all for {output_prefix}")
        self.process_associations()
        print("Associations processed")

        self.build_group_to_gpi()
        print("Group to GPI built")

        self.build_item_to_CSI()
        print("Item to CSI built")

        self.build_group_to_CSI()
        print("Group to CSI built")

        self.process_pairwise_file(self.main_pairwise_file, self.containment_cutoff)
        print(f"Group modularity computed with containment_cutoff {self.containment_cutoff} on {self.main_pairwise_file}")

        self.build_groups_metadata()
        print("Groups metadata built")

        self.export_item_to_gpi_CSI(f"{output_prefix}_item_to_gpi_CSI.tsv")
        print(f"Item to GPI CSI exported at {output_prefix}_item_to_gpi_CSI.tsv")

        self.remove_exact_ochiai_matches()
        print("Exact Ochiai matches removed")

        self.export_groups_metadata(f"{output_prefix}_PRE-DEDUP_groups_metadata.tsv")
        print(f"PRE Groups metadata exported at {output_prefix}_groups_metadata.tsv")

        print("Deduplication in process ...")
        self.graph_based_setcover(self.GC)
        print(f"Deduplication done. Results exported at {output_prefix}_clusters.tsv")

        self.export_groups_metadata(f"{output_prefix}_groups_metadata.tsv")
        print(f"Groups metadata exported at {output_prefix}_groups_metadata.tsv")

        self.df_logging.to_csv(f"{output_prefix}_set_cover_logging", sep='\t', index=False)
        print(f"Logs exported at {output_prefix}_set_cover_logging.tsv")

        self.export_split_group_metadata(output_prefix)
        print(f"Split groups metadata exported at {output_prefix}_remaining_groups_metadata.tsv and {output_prefix}_removed_groups_metadata.tsv")

        self.write_new_associations_file(f"{output_prefix}_associations.tsv")
        print(f"New associations exported at {output_prefix}_associations.tsv")
        print("-------------------------")      
        print(f"original number of groups {self.total_number_of_groups}")
        number_of_groups_removed_from_set_cover = self.total_number_of_groups - self.filtered_exact_ochiai_groups_count - self.final_remaining_groups_count
        percentage_number_of_groups_removed_from_set_cover = 100 * number_of_groups_removed_from_set_cover / self.total_number_of_groups
        number_of_groups_removed_from_exact_duplicates = self.filtered_exact_ochiai_groups_count
        percentage_number_of_groups_removed_from_exact_duplicates = 100 * number_of_groups_removed_from_exact_duplicates / self.total_number_of_groups
        print(f"groups removed due to exact duplication: {self.filtered_exact_ochiai_groups_count} = {percentage_number_of_groups_removed_from_exact_duplicates}%")
        print(f"number of groups removed from set-cover only: {number_of_groups_removed_from_set_cover} = {percentage_number_of_groups_removed_from_set_cover}%")
        print(f"final number of groups {self.final_remaining_groups_count} = {100 * self.final_remaining_groups_count / self.total_number_of_groups}%")
        total_percentage = percentage_number_of_groups_removed_from_exact_duplicates + percentage_number_of_groups_removed_from_set_cover + (100 * self.final_remaining_groups_count / self.total_number_of_groups)
        print(f"[DEBUG] total percentage of groups: {total_percentage}%")
        
        items_stats = self.calculate_items_stats()
        print(f"original number of items {items_stats['original_no_of_items']}")
        print(f"remaining number of items {items_stats['remaining_no_of_items']}")
        print(f"percentage of items remaining {items_stats['percentage_of_items_remaining']}%")
        print(f"original overlap score {items_stats['original_overlap_score']}")
        print(f"remaining overlap score {items_stats['remaining_overlap_score']}")


def modularity_based_deduplication():
    parser = argparse.ArgumentParser(description='Deduplicate groups')
    parser.add_argument('-c', '--clusters_file', type=str, help='clusters file', required=True)
    parser.add_argument('-a', '--asc', type=str, help='associations file', required=True)
    parser.add_argument('-i', '--index_prefix', type=str, help='index prefix', required=True)
    parser.add_argument('-m', '--containment_cutoff', type=int, help='containment_cutoff', default=80)
    parser.add_argument('-e', '--exact_ochiai_cutoff', type=int, help='exact_ochiai_cutoff', default=80)
    parser.add_argument('-g', '--GC', type=int, help='GC', default=100)
    parser.add_argument('-o', '--output_prefix', type=str, help='output_prefix', required=True)

    dedup = DeduplicateGroups(
        clusters_file=parser.parse_args().clusters_file,
        associations_file=parser.parse_args().asc,
        index_prefix=parser.parse_args().index_prefix,
        containment_cutoff=parser.parse_args().containment_cutoff,
        exact_ochiai_cutoff=parser.parse_args().exact_ochiai_cutoff,
        GC=parser.parse_args().GC,
    )         
    dedup.build_all(parser.parse_args().output_prefix)
    DEDUP_GMT = parser.parse_args().output_prefix + "_DEDUP_GMT.gmt"
    ORIGINAL_GMT = parser.parse_args().output_prefix + "_ORIGINAL_GMT.gmt"
    dedup.export_deduplicated_gmt("DEDUP_GMT")
    dedup.export_original_gmt("NO_ABS_GMT_ORIGINAL_GROUPS")


def graph_based_deduplication():
    parser = argparse.ArgumentParser(description='Deduplicate groups')
    parser.add_argument('-c', '--clusters_file', type=str, help='clusters file', required=True)
    parser.add_argument('-a', '--asc', type=str, help='associations file', required=True)
    parser.add_argument('-i', '--index_prefix', type=str, help='index prefix', required=True)
    
    parser.add_argument('-m', '--containment_cutoff', type=int, help='containment_cutoff', default=80)
    parser.add_argument('-e', '--exact_ochiai_cutoff', type=int, help='exact_ochiai_cutoff', default=80)
    parser.add_argument('-g', '--GC', type=int, help='GC', default=100)
    parser.add_argument('-o', '--output_prefix', type=str, help='output_prefix', required=True)

    dedup = GraphBasedDeduplication(
        clusters_file=parser.parse_args().clusters_file,
        associations_file=parser.parse_args().asc,
        index_prefix=parser.parse_args().index_prefix,
        containment_cutoff=parser.parse_args().containment_cutoff,
        exact_ochiai_cutoff=parser.parse_args().exact_ochiai_cutoff,
        GC=parser.parse_args().GC,
    )    
    dedup.build_all(parser.parse_args().output_prefix)
    DEDUP_GMT = parser.parse_args().output_prefix + "_DEDUP_GMT.gmt"
    ORIGINAL_GMT = parser.parse_args().output_prefix + "_ORIGINAL_GMT.gmt"
    dedup.export_deduplicated_gmt("DEDUP_GMT")
    dedup.export_original_gmt("NO_ABS_GMT_ORIGINAL_GROUPS")

def path_to_absolute_path(ctx, param, value):
    return value if value == "NA" else os.path.abspath(value)

@cli.command(name="setcov", epilog = dbretina_doc.doc_url("setcov"), help_priority=10)
@click.option('-i', '--index-prefix', "index_prefix", required=True, type=click.STRING, help="index file prefix")
@click.option('--modularity', "containment_cutoff", required=False, default=80, show_default = True, type=click.FloatRange(0, 100, clamp=False), help="containment cutoff for modularity calculation")
@click.option('--dedup', "ochiai_cutoff", required=False, default=100, show_default = True, type=click.FloatRange(0, 100, clamp=False), help="deduplication similarity cutoff")
@click.option('--community', "ochiai_community_cutoff", required=False, default=30, show_default = True, type=click.FloatRange(0, 100, clamp=False), help="community detection similarity cutoff")
@click.option('--stop-cov', "GC", required=False, default=100, type=click.FloatRange(0, 100, clamp=False), help="stop when items covered by %")
@click.option('-o', '--output', "output_prefix", required=True, default=None, help="output file prefix")
@click.pass_context
def main(ctx, index_prefix, containment_cutoff, ochiai_cutoff, GC, output_prefix, ochiai_community_cutoff):
    """Apply set-cover algorithm.
    """
    # TODO: make temporary association file (parse the json file directly later)
    # TODO: This is DUMP, but will work for now.
    # TODO: change the `def process_associations(self)` function to accept the json file
    
    raw_associations_json = json.load(open(f"{index_prefix}_raw.json"))['data']
    # write the associations to a two columns tsv file
    _tmp_associations = f".DBRetina_tmp_association_{index_prefix}.tsv"
    with open(_tmp_associations, "w") as f:
        f.write("itemset\titem\n")
        for group, items in raw_associations_json.items():
            for item in items:            
                f.write(f"{group}\t{item}\n")
    
    #### Start main code

    dedup =  DeduplicateGroups(
        associations_file=_tmp_associations,
        index_prefix=index_prefix,
        containment_cutoff=containment_cutoff,
        exact_ochiai_cutoff=ochiai_cutoff,
        ochiai_community_cutoff = ochiai_community_cutoff,
        GC=GC,
        ctx = ctx,
    )
    
    dedup.build_all(output_prefix)
    
    if '.DBRetina' in _tmp_associations:
        os.remove(_tmp_associations)
    