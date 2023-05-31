#include "DBRetina.hpp"


namespace DBRetina {
    class GeneSets {
    private:
        kDataFrame* gene_to_color;
        parallel_flat_hash_map<uint32_t, vector<uint32_t>> color_to_ids;
        parallel_flat_hash_map<uint64_t, double> gene_to_PSI;
        parallel_flat_hash_map<string, parallel_flat_hash_set<uint64_t>> pathway_to_gene_set;
        uint64_t n_total_pathways;
        parallel_flat_hash_map<string, double> pathway_to_average_PSI;
        parallel_flat_hash_map<uint32_t, parallel_flat_hash_set<string>> clusterToPathways;
        parallel_flat_hash_map<uint64_t, double> gene_to_PPI;
        parallel_flat_hash_map<string, double> pathway_to_average_PPI;
        parallel_flat_hash_map<uint64_t, uint64_t> gene_to_no_clusters;
        parallel_flat_hash_map<uint32_t, double> cluster_to_average_PPI;

        GeneSets(string associations_file);
        void build_from_index(string index_prefix);
        void build_from_clusters_file(string clusters_file);
        parallel_flat_hash_set<uint64_t> get_universe_set(parallel_flat_hash_set<string> pathways);
        void build_all_pathways_PSI();
        void build_gene_to_PSI();
        uint64_t get_pathway_PSI(string pathway);
        void build_gene_to_no_clusters();
        void build_gene_to_PPI();
        uint64_t build_pathway_to_average_PPI(string pathway);
        uint64_t get_cluster_average_ppi(uint32_t cluster_id);
        void build_clusters_to_average_PPI();

    public:
        unordered_map<string, double> get_pathways_ppi();
        unordered_map<uint32_t, double> get_clusters_ppi();
        unordered_map<string, double> get_pathways_psi();
    };
}
