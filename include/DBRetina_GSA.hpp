#include "kDataFrame.hpp"
#include "parallel_hashmap/phmap.h"
#include <cstdint>
#include <unordered_set>

using namespace std;
using phmap::parallel_flat_hash_map;
using phmap::flat_hash_map;
using phmap::flat_hash_set;
using phmap::parallel_flat_hash_set;


using int_vec_map = parallel_flat_hash_map<uint32_t, vector<uint32_t>, std::hash<uint32_t>, std::equal_to<uint32_t>, std::allocator<std::pair<const uint32_t, vector<uint32_t>>>, 1>;
using str_hashed_set_map = parallel_flat_hash_map<string, flat_hash_set<uint64_t>>;


class GeneSets {
private:
    kDataFrame* gene_to_color;
    int_vec_map color_to_ids;
    flat_hash_map<uint64_t, double> gene_to_PSI;
    // TODO: rename later [pathway-cluster-specificity-index]
    flat_hash_map<uint64_t, double> gene_to_PCSI;
    str_hashed_set_map pathway_to_gene_set;
    uint64_t n_total_pathways;
    double mean_pathway_length;
    string index_prefix;
    flat_hash_map<string, double> pathway_to_average_PSI;
    flat_hash_map<uint32_t, flat_hash_set<string>> clusterToPathways;
    flat_hash_map<uint64_t, double> gene_to_PPI;
    flat_hash_map<uint64_t, string> hashed_gene_to_name;
    flat_hash_map<string, double> pathway_to_average_PPI;
    flat_hash_map<string, double> pathway_to_average_PCSI;
    flat_hash_map<uint64_t, uint64_t> gene_to_no_clusters;
    flat_hash_map<uint32_t, double> cluster_to_average_PPI;
    flat_hash_map<uint32_t, uint32_t> _group_id_to_size;
    flat_hash_map<uint32_t, string> _group_id_to_name;
    unordered_map<string, int> pathway_to_fragmentation;
    unordered_map<string, int> pathway_to_heterogeneity;
    unordered_map<string, int> pathway_to_modularity;
    flat_hash_map<uint64_t, uint64_t> gene_to_associated_pathways;



    flat_hash_set<uint64_t> get_universe_set(flat_hash_set<string>& pathways);
    void build_all_pathways_PSI();
    void calculate_mean_pathway_length();
    void build_gene_to_PSI();
    double get_pathway_PSI(string pathway);
    void build_gene_to_no_clusters();
    void build_gene_to_PPI();
    void build_pathway_to_average_PPI(string pathway);
    double get_cluster_average_ppi(uint32_t cluster_id);
    void build_clusters_to_average_PPI();
    void load_group_sizes(string group_sizes_file);
    void load_namesMap(string namesMapfile);
    vector<string> sort_min_max(const flat_hash_set<string>& pathways);
    vector<string> combine_and_sort(const flat_hash_set<string>& pathways);
    vector<string> combined_sort_frag_ppi_len(flat_hash_set<string>& pathways);
    vector<string> influence_sort_pathway(const flat_hash_set<string>& pathways);

    std::vector<std::pair<std::string, double>> sort_desc_PPI(
        const flat_hash_map<std::string, double>& pathway_scores);

public:
    GeneSets();
    GeneSets(string associations_file);
    void build_from_index(string index_prefix);
    void build_from_clusters_file(string clusters_file);
    unordered_map<string, double> get_pathways_ppi();
    unordered_map<uint32_t, double> get_clusters_ppi();
    unordered_map<string, double> get_pathways_psi();
    // TODO: NEW PCSI
    unordered_map<string, double> get_pathways_pcsi();
    void export_genes_to_ppi_psi_tsv(string filename);
    unordered_map<string, int> get_pathway_lengths();
    void calculate_heterogeneity_and_fragmentation_from_pairwise(string pairwise_file);
    unordered_map<string, int> get_pathway_to_modularity();
    unordered_map<string, int> get_pathway_to_heterogeneity();
    unordered_map<string, int> get_pathway_to_fragmentation();
    void keep_only_these_pathways(string non_redundant_pathways);
    

    flat_hash_set<uint64_t> set_intersection(flat_hash_set<uint64_t>& pathway_genes, flat_hash_set<uint64_t>& universe);
    unordered_map<string, double> greedy_proportional_set_cover(int cluster_id, int GC = 100);
    unordered_map<string, double> proportionalSetCover(int cluster_id, int GC = 100);
    unordered_map<string, double> non_iterative_set_cover(int cluster_id, int GC);

    void clear() {
        this->gene_to_color = nullptr;
        this->cluster_to_average_PPI.clear();
        this->clusterToPathways.clear();
        this->gene_to_no_clusters.clear();
        this->pathway_to_gene_set.clear();
        this->pathway_to_average_PPI.clear();
        this->gene_to_PPI.clear();
        this->gene_to_PSI.clear();
    }
    ~GeneSets();
};
