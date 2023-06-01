#include "kDataFrame.hpp"
#include "parallel_hashmap/phmap.h"
#include <cstdint>

using namespace std;
using phmap::parallel_flat_hash_map;
using phmap::flat_hash_map;
using phmap::flat_hash_set;
using phmap::parallel_flat_hash_set;


using int_vec_map = parallel_flat_hash_map<uint32_t, vector<uint32_t>, std::hash<uint32_t>, std::equal_to<uint32_t>, std::allocator<std::pair<const uint32_t, vector<uint32_t>>>, 1>;
using str_hashed_vec_map = parallel_flat_hash_map<string, parallel_flat_hash_set<uint64_t>>;


class GeneSets {
private:
    kDataFrame* gene_to_color;
    int_vec_map color_to_ids;
    flat_hash_map<uint64_t, double> gene_to_PSI;
    str_hashed_vec_map pathway_to_gene_set;
    uint64_t n_total_pathways;
    flat_hash_map<string, double> pathway_to_average_PSI;
    flat_hash_map<uint32_t, flat_hash_set<string>> clusterToPathways;
    flat_hash_map<uint64_t, double> gene_to_PPI;
    flat_hash_map<string, double> pathway_to_average_PPI;
    flat_hash_map<uint64_t, uint64_t> gene_to_no_clusters;
    flat_hash_map<uint32_t, double> cluster_to_average_PPI;

    flat_hash_set<uint64_t> get_universe_set(flat_hash_set<string> pathways);
    void build_all_pathways_PSI();
    void build_gene_to_PSI();
    double get_pathway_PSI(string pathway);
    void build_gene_to_no_clusters();
    void build_gene_to_PPI();
    void build_pathway_to_average_PPI(string pathway);
    double get_cluster_average_ppi(uint32_t cluster_id);
    void build_clusters_to_average_PPI();

public:
    // empty constructor
    GeneSets();
    GeneSets(string associations_file);
    void build_from_index(string index_prefix);
    void build_from_clusters_file(string clusters_file);
    unordered_map<string, double> get_pathways_ppi();
    unordered_map<uint32_t, double> get_clusters_ppi();
    unordered_map<string, double> get_pathways_psi();
    // Destructor
    ~GeneSets();
    // ~GeneSets() {
    //     this->gene_to_color = nullptr;
    //     this->cluster_to_average_PPI.clear();
    //     this->clusterToPathways.clear();
    //     this->gene_to_no_clusters.clear();
    //     this->pathway_to_gene_set.clear();
    //     this->pathway_to_average_PPI.clear();
    //     this->gene_to_PPI.clear();
    //     this->gene_to_PSI.clear();
    // }
};
