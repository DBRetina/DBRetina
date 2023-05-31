#include "DBRetina_GSA.hpp"
#include <kDataFrame.hpp>
#include <colored_kDataFrame.hpp>



inline void set_to_vector(const phmap::flat_hash_set<uint32_t>& set, vector<uint32_t>& vec) {
    vec.clear();
    vec.reserve(set.size());
    for (auto& i : set) {
        vec.push_back(i);
    }
}

inline void load_colors_to_sources(const std::string& filename, int_vec_map* map)
{
    phmap::BinaryInputArchive ar_in(filename.c_str());
    size_t size;
    ar_in.loadBinary(&size);
    map->reserve(size);
    while (size--)
    {
        uint64_t k;
        phmap::flat_hash_set<uint32_t> v;
        vector<uint32_t> vVec;
        ar_in.loadBinary(&k);
        ar_in.loadBinary(&v);
        set_to_vector(v, vVec);
        map->insert_or_assign(std::move(k), std::move(vVec));
    }
}

namespace DBRetina {

    // constructor
    GeneSets::GeneSets(string associations_file) {
        parse_dbretina_json(associations_file, &pathway_to_gene_set);
        n_total_pathways = pathway_to_gene_set.size();
    }

    void GeneSets::build_from_index(string index_prefix) {
        cerr << "[DEBUG] Loading index from: " << index_prefix << endl;
        gene_to_color = kDataFrame::load(index_prefix);
        string colors_map_file = index_prefix + "_color_to_sources.bin";
        cerr << "[DEBUG] Loading colors map from: " << colors_map_file << endl;
        // We need this to find out how many pathways are associated with each gene
        load_colors_to_sources(colors_map_file, &color_to_ids);
        build_all_pathways_PSI();
    }


    void GeneSets::build_from_clusters_file(string clusters_file) {
        std::ifstream file(clusters_file);

        if (file.is_open()) {
            std::string line;

            // Skip header
            std::getline(file, line);

            // Read cluster IDs and pathways
            while (std::getline(file, line)) {
                std::istringstream iss(line);
                std::string field;
                std::getline(iss, field, '\t');  // cluster_id
                int clusterId = std::stoi(field);
                std::getline(iss, field, '\t');  // number_of_pathways
                std::getline(iss, field, '\t');  // pathways_list

                std::istringstream pathwayStream(field);
                std::string pathway;
                while (std::getline(pathwayStream, pathway, '|')) {
                    this->clusterToPathways[clusterId].insert(pathway);
                }
            }
            file.close();
        }

        build_gene_to_no_clusters();
        build_gene_to_PPI();
        for (auto& [cluster_id, pathways] : this->clusterToPathways)
            for (auto& pathway : pathways) build_pathway_to_average_PPI(pathway);

        build_clusters_to_average_PPI();
    }





    flat_hash_set<uint64_t> GeneSets::get_universe_set(flat_hash_set<string> pathways) {
        // TODO optimize pass by reference
        flat_hash_set<uint64_t> universe_set;
        for (auto& pathway : pathways) {
            // check if that pathway does not exist
            if (pathway_to_gene_set.find(pathway) == pathway_to_gene_set.end())
                throw std::invalid_argument("Pathway does not exist");

            auto genes = pathway_to_gene_set[pathway];
            for (auto& gene : genes) universe_set.insert(gene);
        }
        return universe_set;
    }

    // Pathway Specificity Index PSI
    void GeneSets::build_all_pathways_PSI() {
        cerr << "[DEBUG] Building all pathways PSI" << endl;
        for (const auto& [pathway, _] : pathway_to_gene_set)
            pathway_to_average_PSI[pathway] = get_pathway_PSI(pathway);
    }

    void GeneSets::build_gene_to_PSI() {
        cerr << "[DEBUG] Building gene to PSI map" << endl;
        auto it = gene_to_color->begin();
        while (it != gene_to_color->end()) {
            uint64_t gene = it.getHashedKmer();
            uint64_t color = it.getCount();
            uint64_t pathwayCount = color_to_ids[color].size();
            double psi = static_cast<double>(pathwayCount) / n_total_pathways;
            gene_to_PSI[gene] = psi;
            it++;
        }
    }

    uint64_t GeneSets::get_pathway_PSI(string pathway) {
        // check if that pathway does not exist
        if (pathway_to_gene_set.find(pathway) == pathway_to_gene_set.end())
            throw std::invalid_argument("Pathway does not exist");

        auto genes = pathway_to_gene_set[pathway];
        double average_psi = 0.0;
        for (auto& gene : genes) {
            average_psi += gene_to_PSI[gene];
        }
        return average_psi / genes.size();
    }


    // Pathway Pleiotropy Index PPI
    void GeneSets::build_gene_to_no_clusters() {
        // iterate over clusters
        for (const auto& [cluster_id, pathways] : this->clusterToPathways) {
            auto universe_set = get_universe_set(pathways);
            for (auto& gene : universe_set) {
                this->gene_to_no_clusters[gene]++;
            }
        }
    }

    void GeneSets::build_gene_to_PPI() {
        cerr << "[DEBUG] Building gene to PPI map" << endl;
        for (auto& [gene, no_clusters] : gene_to_no_clusters) {
            double ppi = static_cast<double>(no_clusters) / this->clusterToPathways.size();
            gene_to_PPI[gene] = ppi;
        }
    }

    void GeneSets::build_pathway_to_average_PPI(string pathway) {
        // check if that pathway does not exist
        if (pathway_to_gene_set.find(pathway) == pathway_to_gene_set.end())
            throw std::invalid_argument("Pathway does not exist");

        auto & genes = pathway_to_gene_set[pathway];
        double average_ppi = 0.0;
        for (auto& gene : genes) {
            average_ppi += gene_to_PPI[gene];
        }
        pathway_to_average_PPI[pathway] = average_ppi / genes.size();
    }

    uint64_t GeneSets::get_cluster_average_ppi(uint32_t cluster_id) {
        auto pathways = this->clusterToPathways[cluster_id];
        double average_ppi = 0.0;
        for (auto& pathway : pathways) {
            average_ppi += pathway_to_average_PPI[pathway];
        }
        return average_ppi / pathways.size();
    }

    void GeneSets::build_clusters_to_average_PPI() {
        for (auto& [cluster_id, pathways] : this->clusterToPathways) {
            cluster_to_average_PPI[cluster_id] = get_cluster_average_ppi(cluster_id);
        }
    }


    unordered_map<string, double> GeneSets::get_pathways_ppi() {
        auto& before_convert = this->pathway_to_average_PPI;
        unordered_map<string, double> after_convert;
        for (auto& [pathway, ppi] : before_convert) {
            after_convert[pathway] = ppi;
        }
        return after_convert;
    }

    unordered_map<uint32_t, double> GeneSets::get_clusters_ppi() {
        auto& before_convert = this->cluster_to_average_PPI;
        unordered_map<uint32_t, double> after_convert;
        for (auto& [cluster_id, ppi] : before_convert) {
            after_convert[cluster_id] = ppi;
        }
        return after_convert;
    }

    unordered_map<string, double> GeneSets::get_pathways_psi() {
        auto& before_convert = this->pathway_to_average_PSI;
        unordered_map<string, double> after_convert;
        for (auto& [pathway, psi] : before_convert) {
            after_convert[pathway] = psi;
        }
        return after_convert;
    }

};

