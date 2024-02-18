#!/usr/bin/env python
import pandas as pd
import numpy as np
import csv
from collections import defaultdict
import argparse
import matplotlib.pyplot as plt
import networkx as nx
import math

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


class DeduplicatePathways():
    """
    A class for deduplicating pathways.

    This class provides methods for deduplicating pathways based on various criteria, such as gene coverage, PPI, and PCSI.
    """
    
    communities_clusters_file = ''
    associations_file = ''
    cluster_id_to_pathways = dict()
    pathways_to_genes_no = dict()
    total_number_of_pathways = 0
    main_pairwise_file = ''
    index_prefix = ''
    exact_ochiai_cutoff = 80
    GC = 100
    final_remaining_pathways_count = 0
    final_remaining_pathways = set()
    removed_exact_ochiai_pathways = set()
    
    # Dataframes
    df_genes_metadata = None
    df_pathways_metadata = None
    df_associations = None
    df_gene_to_pcsi = None
    df_pathway_to_avg_pcsi = None
    df_pathway_to_avg_ppi = None
    df_gene_to_ppi = None
    df_pathway_to_modularity = None
    df_logging = None
    df_pathways_per_gene = None
    
    
    
    def __init__(
        self, 
        clusters_file, 
        associations_file, 
        index_prefix,
        GC = 100,
        max_cont_cutoff = 80,
        exact_ochiai_cutoff = 80,
        ):
        self.associations_file = associations_file
        self.communities_clusters_file = clusters_file
        # Initialize a logging DataFrame
        self.df_logging = pd.DataFrame(columns=[
            'pathway', 'no_of_genes', 'coverage %', 'ppi', 'pcsi', 'fragmentation', 'heterogeneity', 'modularity', 'selected', 'cluster_id'
        ])
        self.main_pairwise_file = index_prefix + "_DBRetina_pairwise.tsv"
        self.GC = GC
        self.max_cont_cutoff = max_cont_cutoff
        self.index_prefix = index_prefix
        self.exact_ochiai_cutoff = exact_ochiai_cutoff
    
    # returns a dictionary cluster_id -> pathways
    def cluster_to_pathways(self, clusters_file):
        cluster_id_to_pathways = {}
        header_flag = True
        with open(clusters_file, 'r') as f:
            for line in f:
                if line.startswith("#"): continue
                if header_flag:
                    header_flag = False
                    continue
                cluster_id, cluster_size, cluster_members = line.strip().split("\t")
                cluster_id_to_pathways[int(cluster_id)] = cluster_members.split("|")
        return cluster_id_to_pathways
                
    def process_associations(self):
        self.df_associations = pd.read_csv(self.associations_file, sep="\t")
        # rename 'hgnc_symbol_ids to 'gene'
        self.df_associations.rename(columns={'hgnc_symbol_ids': 'gene'}, inplace=True)
        # all to lower
        self.df_associations['gene'] = self.df_associations['gene'].str.lower()
        self.df_associations['pathway'] = self.df_associations['pathway'].str.lower()
        # calculate pathway to size
        self.pathways_to_genes_no = self.df_associations.groupby('pathway')['gene'].nunique().to_dict()
        # constrcut dictionary of pathway to genes list
        self.pathways_to_genes = self.df_associations.groupby('pathway')['gene'].apply(set).to_dict()


    def build_gene_to_pcsi(self):
        # we have can use df_gene_to_ppi to get the number of pathways and clusters per gene
        
        total_number_of_clusters = self.cluster_id_to_pathways.keys().__len__()
        self.total_number_of_pathways = self.df_associations['pathway'].nunique()
        self.df_pathways_per_gene = self.df_associations.groupby('gene')['pathway'].nunique()
        
        self.df_gene_to_pcsi = pd.DataFrame({
            'gene' : self.df_gene_to_ppi['gene'].values,
            'pathway_count': self.df_gene_to_ppi['pathway_count'].values,
            'cluster_count': self.df_gene_to_ppi['cluster_count'].values,
            'pcsi': 100 * 
                    np.log2(self.df_gene_to_ppi['cluster_count'].values / total_number_of_clusters) / 
                    np.log2(1.0 / total_number_of_clusters)
        })
        
    
    # def old_build_gene_to_pcsi(self):
    #     #1 count number of pathways per gene
    #     self.df_pathways_per_gene = self.df_associations.groupby('gene')['pathway'].nunique()
    #     #2 calculate total number of pathways
    #     self.total_number_of_pathways = self.df_associations['pathway'].nunique()
    #     #3 Calculate the Pathway Specificity Index (PCSI) for each gene
    #     self.df_gene_to_pcsi = pd.DataFrame({
    #         'gene': self.df_pathways_per_gene.index,
    #         'pathways_count': self.df_pathways_per_gene.values,
    #         'pcsi': 100 * np.log2(
    #             self.df_pathways_per_gene.values / self.total_number_of_pathways
    #             ) / np.log2(1.0 / self.total_number_of_pathways)
    #     })
    
    
    def build_pathway_to_pcsi(self):
        #1 merge our original associations with gene_to_pcsi
        merged_df = pd.merge(self.df_associations, self.df_gene_to_pcsi, left_on='gene', right_on='gene', how='left')
        #2 Calculate the sum of PCSI for each pathway and the number of genes per pathway
        pathway_group = merged_df.groupby('pathway')['pcsi'].agg(['sum', 'count'])
        #3 Calculate average PCSI
        self.df_pathway_to_avg_pcsi = pd.DataFrame({
            'average_pcsi': pathway_group['sum'] / pathway_group['count']
        })
        
    def build_pathway_to_ppi(self):
        #1 Get the number of unique clusters and pathways each gene is found in
        self.cluster_id_to_pathways = self.cluster_to_pathways(self.communities_clusters_file)
        # Convert dictionary to dataframe
        cluster_to_pathways_df = pd.DataFrame([(k, v) for k, vv in self.cluster_id_to_pathways.items() for v in vv], 
                                            columns=['cluster_id', 'pathway'])

        #2 map cluster_id to genes
        cluster_to_genes = pd.merge(cluster_to_pathways_df, self.df_associations, on='pathway', how='left')
        cluster_to_genes = pd.merge(cluster_to_pathways_df, self.df_associations, on='pathway', how='left')


        #3. Get the number of unique clusters and pathways each gene is found in
        gene_counts = cluster_to_genes.groupby('gene').agg(
            cluster_count=pd.NamedAgg(column='cluster_id', aggfunc='nunique'),
            pathway_count=pd.NamedAgg(column='pathway', aggfunc='nunique')
        )
        
        #4. Calculate the Pathway Pleiotropy Index (PPI) for each gene (Cd/Nd)
        gene_counts['ppi'] = gene_counts['cluster_count'] / gene_counts['pathway_count']

        #5 All information in one dataframe
        self.df_gene_to_ppi = pd.DataFrame({
            'gene': gene_counts.index,
            'cluster_count': gene_counts['cluster_count'].values,
            'pathway_count': gene_counts['pathway_count'].values,
            'ppi': gene_counts['ppi'].values
        })

        #5. Merge cluster_to_genes with gene_to_ppi to get ppi for each gene in each pathway
        # ['cluster_id', 'pathway', 'gene', 'cluster_count', 'pathway_count', 'ppi']
        merged_df = pd.merge(cluster_to_genes, self.df_gene_to_ppi, on='gene', how='inner')

        # Calculate the average PPI for each pathway
        pathway_group = merged_df.groupby('pathway')['ppi'].mean()

        # Create the DataFrame
        self.df_pathway_to_avg_ppi = pd.DataFrame({
            'pathway': pathway_group.index,
            'average_ppi': 100 * pathway_group.values
        })
             
        
    def export_gene_to_ppi_pcsi(self, file_name):
        # Merge gene_to_ppi with gene_to_pcsi
        merged_df = pd.merge(self.df_gene_to_pcsi, self.df_gene_to_ppi, on='gene', how='left')
        # Export to CSV
        merged_df.to_csv(file_name, index=False, sep='\t')
        return file_name
    
    def export_pathway_to_ppi_pcsi(self, file_name):
        # Merge pathway_to_ppi with pathway_to_pcsi
        merged_df = pd.merge(self.df_pathway_to_avg_pcsi, self.df_pathway_to_avg_ppi, on='pathway', how='left')
        # Export to CSV
        merged_df.to_csv(file_name, index=False, sep='\t')
        return file_name
        
    
    def process_pairwise_file(self, tsv_file, max_cont_threshold=80):
        # Initialize a dictionary to store pathway metrics
        pathway_metrics = defaultdict(lambda: {'fragmentation': 0, 'heterogeneity': 0})
        
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
                # Extract the pathway names from columns 3 and 4
                pathway1 = row[2]
                pathway2 = row[3]
                
                max_containment = float(row[7])
                ochiai_distance = float(row[8])
                
                # instead of re-iterating over the file, do it here
                if ochiai_distance >= self.exact_ochiai_cutoff:
                    self.ochiai_graph.add_edge(pathway1, pathway2)
                
                
                # this for calculating modularity
                if max_containment < max_cont_threshold:
                    continue           

                # Get the lengths of the pathways
                pathway1_len = self.pathways_to_genes_no[pathway1]
                pathway2_len = self.pathways_to_genes_no[pathway2]

                # Compute the metrics
                # Large node fragmented to small nodes
                # Large node heterogenous to small nodes
                if pathway1_len < pathway2_len:
                    pathway_metrics[pathway1]['fragmentation'] -= 1
                    pathway_metrics[pathway2]['heterogeneity'] += 1
                elif pathway1_len > pathway2_len:
                    pathway_metrics[pathway1]['heterogeneity'] += 1
                    pathway_metrics[pathway2]['fragmentation'] -= 1
                # else:  # if they are equal, update both
                #     print("ELSEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE")
                #     pathway_metrics[pathway1]['fragmentation'] -= 1
                #     pathway_metrics[pathway1]['heterogeneity'] += 1
                #     pathway_metrics[pathway2]['fragmentation'] -= 1
                #     pathway_metrics[pathway2]['heterogeneity'] += 1

        # Compute the modularity index for each pathway
        for metrics in pathway_metrics.values():
            fragmentation = metrics['fragmentation']
            heterogeneity = metrics['heterogeneity']
            metrics['modularity'] = abs(fragmentation + heterogeneity)

        # Convert the dictionary to a pandas DataFrame
        self.df_pathway_to_modularity = pd.DataFrame.from_dict(pathway_metrics, orient='index').reset_index()
        self.df_pathway_to_modularity.columns = ['pathway', 'fragmentation', 'heterogeneity', 'modularity']
        
    def build_pathways_metadata(self):
        # This builds the pathways metadata (ppi, pcsi, frag, het, modularity, length)
        #1 merge df_pathway_to_avg_ppi and df_pathway_to_pcsi
        self.df_pathways_metadata = pd.merge(self.df_pathway_to_avg_ppi, self.df_pathway_to_avg_pcsi, on='pathway', how='outer')
        # add deduplication status, default = remain

        #2 merge with df_pathway_indices
        self.df_pathways_metadata = pd.merge(self.df_pathways_metadata, self.df_pathway_to_modularity, on='pathway', how='outer')

        #3 Add no_of_genes column from pathways_to_genes_no dictionary
        self.df_pathways_metadata['no_of_genes'] = self.df_pathways_metadata['pathway'].apply(lambda x: self.pathways_to_genes_no[x])

        #4 Rename the columns
        self.df_pathways_metadata['dedup'] = 'remained'
        self.df_pathways_metadata.columns = ['pathway', 'average_ppi', 'average_pcsi', 'fragmentation', 'heterogeneity', 'modularity', 'no_of_genes', 'dedup']


    def remove_exact_ochiai_matches(self):
        connected_components = self.ochiai_graph.get_connected_components()
        unselected_pathways = set()
        selected_pathways = set()
        for component_pathways in connected_components:
            if len(component_pathways) == 1:
                selected_pathways.add(component_pathways[0])
                continue
            
            # sort the pathways by lowest average_pcsi and highest no_of_genes
            sorted_pathways = self.df_pathways_metadata[self.df_pathways_metadata['pathway'].isin(component_pathways)].sort_values(by=['average_ppi', 'no_of_genes'], ascending=[True, False])
            selected_pathways.add(sorted_pathways.iloc[0]['pathway'])
            unselected_pathways.update(sorted_pathways.iloc[1:]['pathway'].tolist())
                
        # set deduplication status to exact_ochiai
        self.removed_exact_ochiai_pathways = unselected_pathways
        self.filtered_exact_ochiai_pathways_count = len(unselected_pathways)
        
        # set df_pathways_metadata dedup to exact_ochiai if in self.removed_exact_ochiai_pathways
        self.df_pathways_metadata.loc[self.df_pathways_metadata['pathway'].isin(self.removed_exact_ochiai_pathways), 'dedup'] = 'exact_ochiai'

    
    def export_pathways_metadata(self, file_name):
        # Save the dataframe to a CSV file
        # fill NA values with 0
        self.df_pathways_metadata = self.df_pathways_metadata.fillna(0)
        self.df_pathways_metadata.to_csv(file_name, index=False, sep='\t', columns=['pathway', 'no_of_genes', 'average_ppi', 'average_pcsi', 'fragmentation', 'heterogeneity', 'modularity', 'dedup'])
        return file_name

    def cluster_to_universe_set(self, cluster_id):
        # Get pathways for the given cluster ID
        cluster_pathways = self.cluster_id_to_pathways.get(cluster_id, [])

        # Filter df_associations to only include rows with pathways in cluster_pathways
        cluster_df = self.df_associations[self.df_associations['pathway'].isin(cluster_pathways)]

        # Get the unique set of genes associated with these pathways
        return set(cluster_df['gene'].unique())

    def unique_set_of_genes(self):
        return set(self.df_associations['gene'].unique())
    
        

    def modularity_based_set_cover(self, GC = 100):
        
        def _deduplicate_all_pathways_instead_of_communities(GC = 100):
            selected_pathways = {}
            uncovered_genes = self.unique_set_of_genes()
            total_genes = len(uncovered_genes)
            genes_covered = 0
            sorted_pathways_df = self.df_pathways_metadata[
                self.df_pathways_metadata['dedup'] != 'exact_ochiai'].sort_values(
                    by=['modularity', 'average_pcsi', 'no_of_genes'], ascending=[True, False, False]
                    )
            
            # Iterate through sorted pathways
            for _, pathway_row in sorted_pathways_df.iterrows():
                pathway = pathway_row['pathway']
                # Calculate the intersection of the current pathway genes with the uncovered genes
                pathway_genes = self.pathways_to_genes[pathway]
                if common_genes := pathway_genes.intersection(uncovered_genes):
                    # Add the pathway to the set cover
                    coverage = len(common_genes) / total_genes * 100
                    selected_pathways[pathway] = coverage
                    # Remove the covered genes from the uncovered set
                    uncovered_genes -= common_genes
                    # Update the count of covered genes
                    genes_covered += len(common_genes)
                    # If the required gene coverage is met, break the loop
                    if genes_covered / total_genes * 100 >= GC:
                        break

            return selected_pathways
        
        selected_pathways = list(_deduplicate_all_pathways_instead_of_communities(GC).keys())
        self.final_remaining_pathways_count = len(selected_pathways)
        self.final_remaining_pathways.update(selected_pathways)
        
        self.df_pathways_metadata.loc[
            ~self.df_pathways_metadata['pathway'].isin(selected_pathways) &
            ~self.df_pathways_metadata['pathway'].isin(self.removed_exact_ochiai_pathways), 'dedup'] = 'set-cov'

    def write_new_associations_file(self, new_associations_file):
        with open(new_associations_file, 'w') as NEW_ASSOCIATIONS_FILE:
            NEW_ASSOCIATIONS_FILE.write('pathway\tgene\n')
            for pathway in self.final_remaining_pathways:
                for gene in self.pathways_to_genes[pathway]:
                    NEW_ASSOCIATIONS_FILE.write(f'{pathway}\t{gene}\n')    


    def export_split_pathway_metadata(self, file_name):
        remaining_pathways_df = self.df_pathways_metadata[self.df_pathways_metadata['pathway'].isin(self.final_remaining_pathways)]
        removed_pathways_df = self.df_pathways_metadata[~self.df_pathways_metadata['pathway'].isin(self.final_remaining_pathways)]
        remaining_pathways_df.to_csv(f'{file_name}_remaining_pathways_metadata.tsv', index=False, sep='\t', columns=['pathway', 'no_of_genes', 'average_ppi', 'average_pcsi', 'fragmentation', 'heterogeneity', 'modularity'])
        removed_pathways_df.to_csv(f'{file_name}_removed_pathways_metadata.tsv', index=False, sep='\t', columns=['pathway', 'no_of_genes', 'average_ppi', 'average_pcsi', 'fragmentation', 'heterogeneity', 'modularity'])


    def export_deduplicated_gmt(self, file_name):
        file_name = f'{file_name}.gmt'
        pathways_of_interest_df = self.df_associations[
            self.df_associations['pathway'].isin(self.final_remaining_pathways)
        ]
        pathways_to_genes_df = pathways_of_interest_df.groupby('pathway')['gene'].apply(list).reset_index()
        with open(file_name, 'w') as f:
            fake_url = 'http://www.google.com'
            for idx, row in pathways_to_genes_df.iterrows():
                row['gene'] = [gene.upper() for gene in row['gene']]
                f.write('\t'.join([row['pathway'], fake_url] + row['gene']))
                # f.write('\t'.join([row['pathway'], row['pathway']] + row['gene']))
                f.write('\n')

    def export_original_gmt(self, file_name):
        file_name = f'{file_name}.gmt'
        pathways_to_genes_df = self.df_associations.groupby('pathway')['gene'].apply(list).reset_index()
        with open(file_name, 'w') as f:
            for idx, row in pathways_to_genes_df.iterrows():
                # upper case row['gene']
                row['gene'] = [gene.upper() for gene in row['gene']]
                f.write('\t'.join([row['pathway'], row['pathway']] + row['gene']))
                f.write('\n')

    def calculate_genes_stats(self):
        # calculate overlap score for the original data
        original_overlap_score = self.df_pathways_per_gene.mean()
        final_stats = {
            "original_overlap_score": 0,
            "remaining_overlap_score": 0,
            "original_no_of_genes": 0,
            "remaining_no_of_genes": 0,
            "percentage_of_genes_remaining": 0,
            'original_overlap_score': original_overlap_score,
            'original_no_of_genes': len(self.df_pathways_per_gene),
        }
        # create a new dataframe removing pathways not in final_remaining_pathways
        df_remaining_associations = self.df_associations[
            self.df_associations['pathway'].isin(self.final_remaining_pathways)]

        # calculate number of remaining pathways per gene
        df_remaining_pathways_per_gene = df_remaining_associations.groupby('gene')['pathway'].nunique()

        final_stats['remaining_overlap_score'] = df_remaining_pathways_per_gene.mean()
        final_stats['remaining_no_of_genes'] = len(df_remaining_pathways_per_gene)

        final_stats['percentage_of_genes_remaining'] = final_stats['remaining_no_of_genes'] / final_stats['original_no_of_genes'] * 100

        return final_stats
        
    def build_all(self, output_prefix):
        print(f"Building all for {output_prefix}")
        self.process_associations()
        print("Associations processed")

        self.build_pathway_to_ppi()
        print("Pathway to PPI built")

        self.build_gene_to_pcsi()
        print("Gene to PCSI built")

        self.build_pathway_to_pcsi()
        print("Pathway to PCSI built")

        self.process_pairwise_file(self.main_pairwise_file, self.max_cont_cutoff)
        print(f"Pathway modularity computed with max_cont_cutoff {self.max_cont_cutoff} on {self.main_pairwise_file}")

        self.build_pathways_metadata()
        print("Pathways metadata built")

        self.export_gene_to_ppi_pcsi(f"{output_prefix}_gene_to_ppi_pcsi.tsv")
        print(f"Gene to PPI PCSI exported at {output_prefix}_gene_to_ppi_pcsi.tsv")

        self.remove_exact_ochiai_matches()
        print("Exact Ochiai matches removed")

        print("Deduplication in process ...")
        self.modularity_based_set_cover(self.GC)
        print(f"Deduplication done. Results exported at {output_prefix}_clusters.tsv")

        self.export_pathways_metadata(f"{output_prefix}_pathways_metadata.tsv")
        print(f"Pathways metadata exported at {output_prefix}_pathways_metadata.tsv")

        self.df_logging.to_csv(f"{output_prefix}_set_cover_logging", sep='\t', index=False)
        print(f"Logs exported at {output_prefix}_set_cover_logging.tsv")

        self.export_split_pathway_metadata(output_prefix)
        print(f"Split pathways metadata exported at {output_prefix}_remaining_pathways_metadata.tsv and {output_prefix}_removed_pathways_metadata.tsv")

        self.write_new_associations_file(f"{output_prefix}_associations.tsv")
        print(f"New associations exported at {output_prefix}_associations.tsv")
        print("-------------------------")      
        print(f"original number of pathways {self.total_number_of_pathways}")
        number_of_pathways_removed_from_set_cover = self.total_number_of_pathways - self.filtered_exact_ochiai_pathways_count - self.final_remaining_pathways_count
        percentage_number_of_pathways_removed_from_set_cover = 100 * number_of_pathways_removed_from_set_cover / self.total_number_of_pathways
        number_of_pathways_removed_from_exact_duplicates = self.filtered_exact_ochiai_pathways_count
        percentage_number_of_pathways_removed_from_exact_duplicates = 100 * number_of_pathways_removed_from_exact_duplicates / self.total_number_of_pathways
        print(f"pathways removed due to exact duplication: {self.filtered_exact_ochiai_pathways_count} = {percentage_number_of_pathways_removed_from_exact_duplicates}%")
        print(f"number of pathways removed from set-cover only: {number_of_pathways_removed_from_set_cover} = {percentage_number_of_pathways_removed_from_set_cover}%")
        print(f"final number of pathways {self.final_remaining_pathways_count} = {100 * self.final_remaining_pathways_count / self.total_number_of_pathways}%")
        total_percentage = percentage_number_of_pathways_removed_from_exact_duplicates + percentage_number_of_pathways_removed_from_set_cover + (100 * self.final_remaining_pathways_count / self.total_number_of_pathways)
        print(f"[DEBUG] total percentage of pathways: {total_percentage}%")
        
        genes_stats = self.calculate_genes_stats()
        print(f"original number of genes {genes_stats['original_no_of_genes']}")
        print(f"remaining number of genes {genes_stats['remaining_no_of_genes']}")
        print(f"percentage of genes remaining {genes_stats['percentage_of_genes_remaining']}%")
        print(f"original overlap score {genes_stats['original_overlap_score']}")
        print(f"remaining overlap score {genes_stats['remaining_overlap_score']}")


class GraphBasedDeduplication(DeduplicatePathways):
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



        selected_pathways = {}
        uncovered_genes = self.unique_set_of_genes()
        total_genes = len(uncovered_genes)
        genes_covered = 0
        iteration = 1        

        # get initial sorting of pathways
        sorted_pathways_df = self.df_pathways_metadata[
            self.df_pathways_metadata['dedup'] != 'exact_ochiai'].sort_values(
                by=['modularity', 'average_pcsi', 'no_of_genes'], ascending=[True, False, False]
                )

        # drop rows having dedup = 'exact_ochiai'
        sorted_pathways_df = sorted_pathways_df[sorted_pathways_df['dedup'] != 'exact_ochiai']

        # add similarity_bin column to sorted_pathways_df
        sorted_pathways_df['similarity_bin'] = 0
        sorted_pathways_df['similarity'] = 0
        
        sorted_pathways_df.to_csv(f"iteration_{iteration}_loop_sorted_pathways_df.tsv", sep="\t")


        # index sorted_pathways_df by pathway
        sorted_pathways_df = sorted_pathways_df.set_index('pathway')


        graph_remaining_pathways_df = pd.DataFrame(columns = sorted_pathways_df.columns)
        first_pathway_name = sorted_pathways_df.iloc[0].name
        print(f"[DEBUG] first pathway name: {first_pathway_name}")

        # load pairwise tsv file
        df_edges = pd.read_csv(self.main_pairwise_file, sep='\t', comment='#', usecols=['group_1_name', 'group_2_name', 'ochiai'])
        

        # keep edges that have group_1_name and group_2_name in sorted_pathways_df
        df_edges = df_edges[df_edges['group_1_name'].isin(sorted_pathways_df.index) & df_edges['group_2_name'].isin(sorted_pathways_df.index)]

        # Create a graph with the edges
        G = nx.from_pandas_edgelist(df_edges, source='group_1_name', target='group_2_name', edge_attr='ochiai')
        
        # add nodes to the graph that does not have any edges
        G.add_nodes_from(sorted_pathways_df.index.difference(G.nodes()))
        
        
        G_pos = nx.spring_layout(G)


        # lambdas
        get_all_neighbors = lambda x: list(G.neighbors(x))
        get_all_neighbors_edges = lambda x: [(x, y) for y in get_all_neighbors(x)]
        get_farthest_neighbor = lambda x: max(get_all_neighbors(x), key=lambda y: G[x][y]['ochiai'])
        get_average_ochiai = lambda x: np.mean([G[x][y]['ochiai'] for y in get_all_neighbors(x)])
        get_sum_ochiai = lambda x: np.sum([G[x][y]['ochiai'] for y in get_all_neighbors(x)])
        get_all_neighbor_weights_dict = lambda x: {y: G[x][y]['ochiai'] for y in get_all_neighbors(x)}


        # first update the similarity_bin of the whole table
        for neighbor, weight in get_all_neighbor_weights_dict(first_pathway_name).items():
            print(f"[DEBUG] neighbor: {neighbor}, weight: {weight}")
            sorted_pathways_df.loc[neighbor, 'similarity'] += weight
            sorted_pathways_df.loc[neighbor, 'similarity_bin'] += math.ceil(sorted_pathways_df.loc[neighbor, 'similarity'] / 10)

        print(f"[DEBUG] graph_remaining_pathways_df: {graph_remaining_pathways_df.head()}")

        # plot_and_highligh_specific_node(G, first_pathway_name, f"iteration_{iteration}_first_pathway.png")

        # remove the first pathway from sorted_pathways_df and from the graph
        
        # TODO REMOVE DEBUG
        # sorted_pathways_df.to_csv(f"iteration_ZERO_{iteration}_sorted_pathways_df.tsv", sep="\t")
        sorted_pathways_df.drop(first_pathway_name, inplace=True)
        print(f"[REMOVED] {first_pathway_name} removed |  Graph nodes: {len(G.nodes())}")
        
        print(f"[DEBUG] Graph size: {len(G.nodes())}")
        G.remove_node(first_pathway_name)
        print(f"[DEBUG] {first_pathway_name} removed |  Graph size: {len(G.nodes())}")
        # remove first_pathway_name edges

        # update uncovered_genes with pathways_to_genes[first_pathway_name]
        intersected_genes = self.pathways_to_genes[first_pathway_name].intersection(uncovered_genes)
        # update uncovered_genes
        uncovered_genes = uncovered_genes.difference(intersected_genes)
        genes_covered += len(intersected_genes)


        # add first pathway to graph_remaining_pathways_df without the simularity_bin column, ignore the index
        self.final_remaining_pathways_count += 1
        self.final_remaining_pathways.add(first_pathway_name)
        selected_pathways[first_pathway_name] = {'similarity_bin': len(intersected_genes)}

        # if len(uncovered_genes) == 0:
        #     return graph_remaining_pathways_df.index.tolist()

        while True and sorted_pathways_df.shape[0] > 0:
            iteration += 1
            # TODO REMOVE DEBUG
            # sorted_pathways_df.to_csv(f"iteration_{iteration}_sorted_pathways_df.tsv", sep="\t")

            # sort pathways by lowest similarity_bin, lowest modularity, highest average_pcsi, highest no_of_genes
            print(f"Iteration {iteration} | number of sorted_pathways_df = {len(sorted_pathways_df)}")
            sorted_pathways_df = sorted_pathways_df.sort_values(
                by=['similarity_bin', 'modularity', 'average_pcsi', 'no_of_genes'], ascending=[True, True, False, False]
                )

            # get the first pathway
            first_pathway_name = sorted_pathways_df.iloc[0].name
            print(f"[CURRENT PATHWAY] {first_pathway_name} |  Graph nodes: {len(G.nodes())} | Uncovered genes: {len(uncovered_genes)} | Genes covered: {genes_covered} | Iteration: {iteration} | Remaining pathways: {len(sorted_pathways_df)}")
            # export graph nodes to TSV
            intersected_genes = self.pathways_to_genes[first_pathway_name].intersection(uncovered_genes)
            print(f"[DEBUG] Next pathway to remove: {first_pathway_name}")

            pathway_neighbors = get_all_neighbor_weights_dict(first_pathway_name)
            for neighbor, weight in pathway_neighbors.items():
                print(f"[DEBUG] neighbor: {neighbor}, weight: {weight}")
                sorted_pathways_df.loc[neighbor, 'similarity'] *= (iteration - 1)
                sorted_pathways_df.loc[neighbor, 'similarity'] += weight
                sorted_pathways_df.loc[neighbor, 'similarity'] /= iteration
                sorted_pathways_df.loc[neighbor, 'similarity_bin'] += math.ceil(sorted_pathways_df.loc[neighbor, 'similarity'] / 10)


            # update the graph_remaining_pathways_df with the first pathway
            # add first pathway to graph_remaining_pathways_df without the simularity_bin column, ignore the index
            self.final_remaining_pathways_count += 1
            self.final_remaining_pathways.add(first_pathway_name)

            # graph_remaining_pathways_df = graph_remaining_pathways_df.append(sorted_pathways_df.loc[first_pathway_name])
            # plot_and_highligh_specific_node(G, first_pathway_name, f"iteration_{iteration}_first_pathway.png")

            # remove the first pathway from sorted_pathways_df and from the graph
            G.remove_node(first_pathway_name)
            sorted_pathways_df.drop(first_pathway_name, inplace=True)
            print(f"[REMOVED] {first_pathway_name} removed |  Graph size: {len(G.nodes())}")
            
            if len(intersected_genes) == 0:
                print(f"!!! [DEBUG] iteration: {iteration}, first_pathway_name: {first_pathway_name}, intersected_genes: {intersected_genes}, uncovered_genes: {len(uncovered_genes)}")
                iteration -= 1
                continue
            
            # update uncovered_genes
            uncovered_genes = uncovered_genes.difference(intersected_genes)

            # update genes_covered
            genes_covered += len(intersected_genes)
            selected_pathways[first_pathway_name] = {'similarity_bin': len(intersected_genes)}


            print(f"\n[DEBUG] genes_coverage until now: {genes_covered} / {total_genes} ({genes_covered / total_genes * 100}%)")
            if genes_covered / total_genes * 100 >= GC:
                break


        print(f"[DEBUG] selected_pathways: {selected_pathways}")
        self.final_remaining_pathways_count = len(selected_pathways)
        self.final_remaining_pathways.update(selected_pathways)

        self.df_pathways_metadata.loc[
            ~self.df_pathways_metadata['pathway'].isin(selected_pathways) &
            ~self.df_pathways_metadata['pathway'].isin(self.removed_exact_ochiai_pathways), 'dedup'] = 'set-cov'
        


    def build_all(self, output_prefix):
        print(f"Building all for {output_prefix}")
        self.process_associations()
        print("Associations processed")

        self.build_pathway_to_ppi()
        print("Pathway to PPI built")

        self.build_gene_to_pcsi()
        print("Gene to PCSI built")

        self.build_pathway_to_pcsi()
        print("Pathway to PCSI built")

        self.process_pairwise_file(self.main_pairwise_file, self.max_cont_cutoff)
        print(f"Pathway modularity computed with max_cont_cutoff {self.max_cont_cutoff} on {self.main_pairwise_file}")

        self.build_pathways_metadata()
        print("Pathways metadata built")

        self.export_gene_to_ppi_pcsi(f"{output_prefix}_gene_to_ppi_pcsi.tsv")
        print(f"Gene to PPI PCSI exported at {output_prefix}_gene_to_ppi_pcsi.tsv")

        self.remove_exact_ochiai_matches()
        print("Exact Ochiai matches removed")

        self.export_pathways_metadata(f"{output_prefix}_PRE-DEDUP_pathways_metadata.tsv")
        print(f"PRE Pathways metadata exported at {output_prefix}_pathways_metadata.tsv")

        print("Deduplication in process ...")
        self.graph_based_setcover(self.GC)
        print(f"Deduplication done. Results exported at {output_prefix}_clusters.tsv")

        self.export_pathways_metadata(f"{output_prefix}_pathways_metadata.tsv")
        print(f"Pathways metadata exported at {output_prefix}_pathways_metadata.tsv")

        self.df_logging.to_csv(f"{output_prefix}_set_cover_logging", sep='\t', index=False)
        print(f"Logs exported at {output_prefix}_set_cover_logging.tsv")

        self.export_split_pathway_metadata(output_prefix)
        print(f"Split pathways metadata exported at {output_prefix}_remaining_pathways_metadata.tsv and {output_prefix}_removed_pathways_metadata.tsv")

        self.write_new_associations_file(f"{output_prefix}_associations.tsv")
        print(f"New associations exported at {output_prefix}_associations.tsv")
        print("-------------------------")      
        print(f"original number of pathways {self.total_number_of_pathways}")
        number_of_pathways_removed_from_set_cover = self.total_number_of_pathways - self.filtered_exact_ochiai_pathways_count - self.final_remaining_pathways_count
        percentage_number_of_pathways_removed_from_set_cover = 100 * number_of_pathways_removed_from_set_cover / self.total_number_of_pathways
        number_of_pathways_removed_from_exact_duplicates = self.filtered_exact_ochiai_pathways_count
        percentage_number_of_pathways_removed_from_exact_duplicates = 100 * number_of_pathways_removed_from_exact_duplicates / self.total_number_of_pathways
        print(f"pathways removed due to exact duplication: {self.filtered_exact_ochiai_pathways_count} = {percentage_number_of_pathways_removed_from_exact_duplicates}%")
        print(f"number of pathways removed from set-cover only: {number_of_pathways_removed_from_set_cover} = {percentage_number_of_pathways_removed_from_set_cover}%")
        print(f"final number of pathways {self.final_remaining_pathways_count} = {100 * self.final_remaining_pathways_count / self.total_number_of_pathways}%")
        total_percentage = percentage_number_of_pathways_removed_from_exact_duplicates + percentage_number_of_pathways_removed_from_set_cover + (100 * self.final_remaining_pathways_count / self.total_number_of_pathways)
        print(f"[DEBUG] total percentage of pathways: {total_percentage}%")
        
        genes_stats = self.calculate_genes_stats()
        print(f"original number of genes {genes_stats['original_no_of_genes']}")
        print(f"remaining number of genes {genes_stats['remaining_no_of_genes']}")
        print(f"percentage of genes remaining {genes_stats['percentage_of_genes_remaining']}%")
        print(f"original overlap score {genes_stats['original_overlap_score']}")
        print(f"remaining overlap score {genes_stats['remaining_overlap_score']}")


def modularity_based_deduplication():
    parser = argparse.ArgumentParser(description='Deduplicate pathways')
    parser.add_argument('-c', '--clusters_file', type=str, help='clusters file', required=True)
    parser.add_argument('-a', '--asc', type=str, help='associations file', required=True)
    parser.add_argument('-i', '--index_prefix', type=str, help='index prefix', required=True)
    parser.add_argument('-m', '--max_cont_cutoff', type=int, help='max_cont_cutoff', default=80)
    parser.add_argument('-e', '--exact_ochiai_cutoff', type=int, help='exact_ochiai_cutoff', default=80)
    parser.add_argument('-g', '--GC', type=int, help='GC', default=100)
    parser.add_argument('-o', '--output_prefix', type=str, help='output_prefix', required=True)

    dedup = DeduplicatePathways(
        clusters_file=parser.parse_args().clusters_file,
        associations_file=parser.parse_args().asc,
        index_prefix=parser.parse_args().index_prefix,
        max_cont_cutoff=parser.parse_args().max_cont_cutoff,
        exact_ochiai_cutoff=parser.parse_args().exact_ochiai_cutoff,
        GC=parser.parse_args().GC,
    )         
    dedup.build_all(parser.parse_args().output_prefix)
    DEDUP_GMT = parser.parse_args().output_prefix + "_DEDUP_GMT.gmt"
    ORIGINAL_GMT = parser.parse_args().output_prefix + "_ORIGINAL_GMT.gmt"
    dedup.export_deduplicated_gmt("DEDUP_GMT")
    dedup.export_original_gmt("NO_ABS_GMT_ORIGINAL_PATHWAYS")


def graph_based_deduplication():
    parser = argparse.ArgumentParser(description='Deduplicate pathways')
    parser.add_argument('-c', '--clusters_file', type=str, help='clusters file', required=True)
    parser.add_argument('-a', '--asc', type=str, help='associations file', required=True)
    parser.add_argument('-i', '--index_prefix', type=str, help='index prefix', required=True)
    parser.add_argument('-m', '--max_cont_cutoff', type=int, help='max_cont_cutoff', default=80)
    parser.add_argument('-e', '--exact_ochiai_cutoff', type=int, help='exact_ochiai_cutoff', default=80)
    parser.add_argument('-g', '--GC', type=int, help='GC', default=100)
    parser.add_argument('-o', '--output_prefix', type=str, help='output_prefix', required=True)

    dedup = GraphBasedDeduplication(
        clusters_file=parser.parse_args().clusters_file,
        associations_file=parser.parse_args().asc,
        index_prefix=parser.parse_args().index_prefix,
        max_cont_cutoff=parser.parse_args().max_cont_cutoff,
        exact_ochiai_cutoff=parser.parse_args().exact_ochiai_cutoff,
        GC=parser.parse_args().GC,
    )         
    dedup.build_all(parser.parse_args().output_prefix)
    DEDUP_GMT = parser.parse_args().output_prefix + "_DEDUP_GMT.gmt"
    ORIGINAL_GMT = parser.parse_args().output_prefix + "_ORIGINAL_GMT.gmt"
    dedup.export_deduplicated_gmt("DEDUP_GMT")
    dedup.export_original_gmt("NO_ABS_GMT_ORIGINAL_PATHWAYS")


if __name__ == '__main__':
    # modularity_based_deduplication()
    graph_based_deduplication()
    