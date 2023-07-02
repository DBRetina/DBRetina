#include "DBRetina.hpp"
#include "DBRetina_GSA.hpp"
#include <kDataFrame.hpp>
#include <colored_kDataFrame.hpp>
#include <cstdlib>
#include <filesystem>

inline void set_to_vector(const phmap::flat_hash_set<uint32_t>& set, vector<uint32_t>& vec) {
    vec.clear();
    vec.reserve(set.size());
    for (auto& i : set) {
        vec.push_back(i);
    }
}

inline void create_hash_to_gene_name(const string& json_file, flat_hash_map<uint64_t, string>& hash_to_gene_name) {
    string_hasher hasher = string_hasher();
    string output_prefix = json_file.substr(0, json_file.find_last_of(".")).substr(json_file.find_last_of("/") + 1);

    zstr::ifstream sig_stream(json_file);
    json::value json = json::parse(sig_stream);
    string filetype = string(json["metadata"]["filetype"].as_string());
    auto genes = json["data"].as_object();

    if (filetype == "private") {
        for (auto it = genes.begin(); it != genes.end(); ++it) {
            string parent_name = it->first;
            auto gene = it->second.as_array();
            for (auto it2 = gene.begin(); it2 != gene.end(); ++it2) {
                string str_child = it2->as_string();
                uint64_t hashed_child = hasher(str_child);
                hash_to_gene_name[hashed_child] = str_child;
            }
        }

    }
    else if (filetype == "public") {
        throw "No, it's not expected to use the hashed version!\n";
    }
    else {
        cout << "Unknown filetype: " << filetype << endl;
        throw "Unknown filetype";
    }
}



inline tuple<unordered_map<int, int>, unordered_map<int, vector<int>>, unordered_map<int, vector<int>>> readData(const string& filename) {
    ifstream inFile(filename);
    string line;
    unordered_map<int, int> groupSizes;
    unordered_map<int, vector<int>> outboundEdges;
    unordered_map<int, vector<int>> inboundEdges;

    // skip headers
    getline(inFile, line);
    while (line[0] == '#')
        getline(inFile, line);

    while (getline(inFile, line)) {
        istringstream iss(line);
        int group1ID, group2ID;
        string token;
        // Read group_1_ID
        getline(iss, token, '\t');
        group1ID = stoi(token);
        // Read group_2_ID
        getline(iss, token, '\t');
        group2ID = stoi(token);
        // Assume that you have already filled groupSizes with correct values
        if (groupSizes[group1ID] < groupSizes[group2ID]) {
            outboundEdges[group1ID].push_back(group2ID);
            inboundEdges[group2ID].push_back(group1ID);
        }
        else {
            outboundEdges[group2ID].push_back(group1ID);
            inboundEdges[group1ID].push_back(group2ID);
        }
    }

    return make_tuple(groupSizes, outboundEdges, inboundEdges);
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

GeneSets::GeneSets() {

}

// constructor
GeneSets::GeneSets(string associations_file) {
    parse_dbretina_json(associations_file, &pathway_to_gene_set);
    n_total_pathways = pathway_to_gene_set.size();
    calculate_mean_pathway_length();
    cerr << "[DEBUG] Total number of pathways: " << n_total_pathways << endl;
    cerr << "[DEBUG] Mean pathway length: " << mean_pathway_length << endl;
}

void GeneSets::build_from_index(string index_prefix) {
    this->index_prefix = index_prefix;
    cerr << "[DEBUG] Loading index from: " << index_prefix << endl;
    gene_to_color = kDataFrame::load(index_prefix);
    string colors_map_file = index_prefix + "_color_to_sources.bin";
    cerr << "[DEBUG] Loading colors map from: " << colors_map_file << endl;
    // We need this to find out how many pathways are associated with each gene
    load_colors_to_sources(colors_map_file, &color_to_ids);
    cout << "Loaded colors: " << color_to_ids.size() << endl;
    build_gene_to_PSI();
    build_all_pathways_PSI();
    create_hash_to_gene_name(index_prefix + "_raw.json", this->hashed_gene_to_name);
}


void GeneSets::build_from_clusters_file(string clusters_file) {
    cerr << "[DEBUG] Loading clusters from: " << clusters_file << endl;

    std::ifstream file(clusters_file);

    if (file.is_open()) {
        std::string line;

        // skip comments
        while (std::getline(file, line) && line[0] == '#') {
            // skip
        }


        // Read cluster IDs and pathways
        while (std::getline(file, line)) {
            std::istringstream iss(line);
            std::string field;
            std::getline(iss, field, '\t');  // cluster_id
            uint32_t clusterId = std::stoi(field);
            std::getline(iss, field, '\t');  // number_of_pathways
            std::getline(iss, field, '\t');  // pathways_list

            std::istringstream pathwayStream(field);
            std::string pathway;
            this->clusterToPathways[clusterId] = flat_hash_set<string>();
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


void GeneSets::load_group_sizes(string bin_file) {
    phmap::BinaryInputArchive ar_in_kmer_count(bin_file.c_str());
    this->_group_id_to_size.phmap_load(ar_in_kmer_count);
    assert(this->_group_id_to_size.size());
}

void GeneSets::load_namesMap(string filename) {
    std::ifstream inputFile(filename);

    if (!inputFile.is_open()) {
        throw std::runtime_error("Error opening the file: " + filename);
    }

    std::string line;
    std::getline(inputFile, line); // skip first line
    while (std::getline(inputFile, line)) {
        std::istringstream lineStream(line);
        std::string column1, column2;
        if (std::getline(lineStream, column1, '|') && std::getline(lineStream, column2, '|')) {
            transform(column2.begin(), column2.end(), column2.begin(), ::tolower);
            this->_group_id_to_name.operator[](stoi(column1)) = column2;
        }
        else {
            inputFile.close();
            throw std::runtime_error("Invalid line format: '" + line + "'");
        }
    }

    inputFile.close();
}


unordered_map<string, int> GeneSets::get_pathway_to_heterogeneity() {
    return pathway_to_heterogeneity;
}

unordered_map<string, int> GeneSets::get_pathway_to_fragmentation() {
    return pathway_to_fragmentation;
}

unordered_map<string, int> GeneSets::get_pathway_to_modularity() {
    return pathway_to_modularity;
}


void GeneSets::keep_only_these_pathways(string non_redundant_pathways) {

    flat_hash_map<string, bool> remaining_pathways;
    ifstream inFile(non_redundant_pathways);
    string line;
    while (getline(inFile, line)) {
        string parsed_pathway = line;
        if (pathway_to_gene_set.find(parsed_pathway) == pathway_to_gene_set.end()) {
            throw std::runtime_error("Pathway not found: " + parsed_pathway);
        }
        remaining_pathways[line] = true;
    }
    inFile.close();

    for (auto& [cluster_id, pathways] : this->clusterToPathways) {
        for (auto& pathway : pathways) {
            if (remaining_pathways.find(pathway) == remaining_pathways.end()) {
                this->clusterToPathways[cluster_id].erase(pathway);
            }
        }
    }

}



void GeneSets::calculate_heterogeneity_and_fragmentation_from_pairwise(string pairwise_file) {

    if (!std::filesystem::exists(pairwise_file)) {
        throw std::runtime_error("File not found: " + pairwise_file);
    }

    ifstream inFile(pairwise_file);
    string _group_to_size_file = this->index_prefix + "_groupID_to_featureCount.bin";
    this->load_group_sizes(_group_to_size_file);
    string _namesMap = this->index_prefix + ".namesMap";
    this->load_namesMap(_namesMap);

    string line;
    flat_hash_map<uint32_t, vector<uint32_t>> outboundEdges;
    flat_hash_map<uint32_t, vector<uint32_t>> inboundEdges;

    // skip headers
    getline(inFile, line);
    while (line[0] == '#')
        getline(inFile, line);

    cerr << "[DEBUG] counting in/out bounds from pairwise file: " << pairwise_file << endl;
    while (getline(inFile, line)) {
        istringstream iss(line);
        uint32_t group1ID, group2ID;
        string token;
        // Read group_1_ID
        getline(iss, token, '\t');
        group1ID = stoi(token);
        // Read group_2_ID
        getline(iss, token, '\t');
        group2ID = stoi(token);

        // Create edges based on group size
        if (this->_group_id_to_size.at(group1ID) < _group_id_to_size.at(group2ID)) {
            outboundEdges[group1ID].push_back(group2ID);
            inboundEdges[group2ID].push_back(group1ID);
        }
        else {
            outboundEdges[group2ID].push_back(group1ID);
            inboundEdges[group1ID].push_back(group2ID);
        }
    }

    // Save heterogeneity and fragmentation to this->pathway_to_heterogeneity and this->pathway_to_fragmentation
    // Calculate fragmentation and heterogeneity scores
    cerr << "[DEBUG] calculating fragmentation and heterogeneity scores" << endl;
    for (const auto& pair : this->_group_id_to_name) {
        int id = pair.first;
        string pathway_name = pair.second;
        int fragmentationScore = -outboundEdges[id].size();
        int heterogeneityScore = inboundEdges[id].size();
        this->pathway_to_fragmentation.insert({ pathway_name, fragmentationScore });
        this->pathway_to_heterogeneity.insert({ pathway_name, heterogeneityScore });
        int absolute_total = abs(fragmentationScore + heterogeneityScore);
        this->pathway_to_modularity.insert({ pathway_name, absolute_total });
    }
}


inline double normalize(int value, int min, int max) {
    return (double)(value - min) / (max - min);
}

inline double normalize(double value, double min, double max) {
    return (value - min) / (max - min);
}


vector<string> GeneSets::sort_min_max(const flat_hash_set<string>& pathways) {
    // Convert hashset to vector of tuples
    vector<tuple<string, double, int>> vec;

    // Calculate min and max for normalization
    auto minMaxFrag = minmax_element(pathway_to_fragmentation.begin(), pathway_to_fragmentation.end(),
        [](const auto& a, const auto& b) {
            return a.second < b.second;
        });

    auto minMaxPPI = minmax_element(pathway_to_average_PPI.begin(), pathway_to_average_PPI.end(),
        [](const auto& a, const auto& b) {
            return a.second < b.second;
        });

    // Min and max values
    int minFrag = minMaxFrag.first->second;
    int maxFrag = minMaxFrag.second->second;
    double minPPI = minMaxPPI.first->second;
    double maxPPI = minMaxPPI.second->second;

    for (const auto& pathway : pathways) {
        double normFrag = normalize(pathway_to_fragmentation[pathway], minFrag, maxFrag);
        double normPPI = normalize(pathway_to_average_PPI[pathway], minPPI, maxPPI);
        double compositeScore = (normFrag + normPPI) / 2;
        vec.push_back(make_tuple(pathway, compositeScore, pathway_to_gene_set[pathway].size()));
    }

    // Sort the vector based on the criteria
    sort(vec.begin(), vec.end(), [](const auto& a, const auto& b) {
        auto [path_a, score_a, len_a] = a;
        auto [path_b, score_b, len_b] = b;
        if (tie(score_a, len_a) == tie(score_b, len_b)) {
            // Randomly pick one if there is a complete tie
            return (rand() % 2) == 0;
        }
        return tie(score_a, len_a) < tie(score_b, len_b);
        });

    // Extract sorted keys
    vector<string> sortedKeys;
    transform(vec.begin(), vec.end(), back_inserter(sortedKeys), [](const auto& tuple) {
        return get<0>(tuple);
        });

    return sortedKeys;
}

vector<string> GeneSets::combine_and_sort(const flat_hash_set<string>& pathways) {

    /*
    For each pathway, we fetch the corresponding PPI and fragmentation scores
    from pathway_to_average_PPI and pathway_to_fragmentation, respectively,
    and store these along with the pathway name in the combined vector.
    We then sort the combined vector based on both the PPI and fragmentation scores,
    and finally extract the sorted pathway names from this vector.
    */

    vector<pair<string, pair<double, int>>> combined;

    for (const auto& pathway : pathways) {
        combined.push_back(make_pair(pathway, make_pair(pathway_to_average_PPI[pathway], pathway_to_fragmentation[pathway])));
    }

    sort(combined.begin(), combined.end(), [](const auto& a, const auto& b) {
        return tie(a.second.first, a.second.second) < tie(b.second.first, b.second.second);
        });

    vector<string> sortedKeys;
    for (const auto& item : combined) {
        sortedKeys.push_back(item.first);
    }

    return sortedKeys;
}


vector<string> GeneSets::combined_sort_frag_ppi_len(flat_hash_set<string>& pathways) {

    // Convert to a vector for sorting
    vector<string> pathways_vec(pathways.begin(), pathways.end());

    // Custom sort
    sort(pathways_vec.begin(), pathways_vec.end(), [&](const string& a, const string& b) {
        if (pathway_to_fragmentation[a] != pathway_to_fragmentation[b])
            return pathway_to_fragmentation[a] < pathway_to_fragmentation[b];  // Lower fragmentation first
        if (pathway_to_average_PPI[a] != pathway_to_average_PPI[b])
            return pathway_to_average_PPI[a] < pathway_to_average_PPI[b];  // Lower PPI next
        // Finally, larger pathway size
        return pathway_to_gene_set[a].size() > pathway_to_gene_set[b].size();
        });

    return pathways_vec;
}




flat_hash_set<uint64_t> GeneSets::get_universe_set(flat_hash_set<string>& pathways) {
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

flat_hash_set<uint64_t> GeneSets::set_intersection(flat_hash_set<uint64_t>& pathway_genes, flat_hash_set<uint64_t>& universe) {
    // check if that pathway does not exist

    const auto& smaller = pathway_genes.size() <= universe.size() ? pathway_genes : universe;
    const auto& larger = pathway_genes.size() > universe.size() ? pathway_genes : universe;

    flat_hash_set<uint64_t> intersected_genes;

    for (const auto& element : smaller)
        if (larger.contains(element))
            intersected_genes.insert(element);

    return intersected_genes;
}

unordered_map<string, int> GeneSets::get_pathway_lengths() {
    unordered_map<string, int> pathway_lengths;
    for (auto& [pathway, genes] : pathway_to_gene_set) {
        pathway_lengths[pathway] = genes.size();
    }
    return pathway_lengths;
}



/*
    ------------------------------------------
        Proportional Set-Cover PSC
    ------------------------------------------
*/

// Sort the pathways by their PPI, resolve tie by length, then randomly
std::vector<std::pair<std::string, double>> GeneSets::sort_desc_PPI(
    const flat_hash_map<std::string, double>& pathway_scores) {

    // Copy the data to a vector
    std::vector<std::pair<std::string, double>> score_vec(pathway_scores.begin(),
        pathway_scores.end());

    // Random device
    std::random_device rd;
    std::mt19937 gen(rd());

    // Custom sort
    std::sort(score_vec.begin(), score_vec.end(),
        [this](const auto& a, const auto& b) {
            if (a.second != b.second) {
                return a.second < b.second; // Sort by score
            }
            else {
                // In case of a tie, sort by length
                return this->pathway_to_gene_set.at(a.first).size() > this->pathway_to_gene_set.at(b.first).size();
            }
        }
    );

    // Custom unique operation
    auto it = std::unique(score_vec.begin(), score_vec.end(),
        [this, &gen](auto& a, auto& b) {
            bool same_score = a.second == b.second;
            if (same_score) {
                bool same_length = this->pathway_to_gene_set.at(a.first).size() == this->pathway_to_gene_set.at(b.first).size();
                if (!same_length) {
                    // If scores are equal, keep the one with the longest pathway
                    if (this->pathway_to_gene_set.at(a.first).size() < this->pathway_to_gene_set.at(b.first).size()) {
                        a.first = b.first;
                        a.second = b.second;
                    }
                }
                else {
                    // If scores and lengths are the same, pick one at random
                    std::uniform_int_distribution<> dis(0, 1);
                    if (dis(gen) == 0) {
                        a.first = b.first;
                        a.second = b.second;
                    }
                }
            }
            return same_score;
        }
    );

    score_vec.erase(it, score_vec.end());

    return score_vec;
}

vector<string> GeneSets::influence_sort_pathway(const flat_hash_set<string>& pathways) {
    vector<pair<string, double>> scores;

    // Calculate maximum fragmentation and pathway size for normalization
    int maxFragmentation = 0;
    int maxSize = 0;
    for (const auto& pathway : pathways) {
        maxFragmentation = max(maxFragmentation, pathway_to_fragmentation[pathway]);
        maxSize = max(maxSize, static_cast<int>(pathway_to_gene_set[pathway].size()));
    }

    // Calculate combined scores
    for (const auto& pathway : pathways) {
        double normalizedFragmentation = static_cast<double>(pathway_to_fragmentation[pathway]) / maxFragmentation;
        double normalizedPPI = 1.0 - pathway_to_average_PPI[pathway];  // subtract from 1 because lower PPI is better
        double PSI = pathway_to_average_PSI[pathway];
        double normalizedSize = static_cast<double>(pathway_to_gene_set[pathway].size()) / maxSize;

        // Combine normalized scores with respective weights
        double score = 0.4 * normalizedFragmentation + 0.3 * normalizedPPI + 0.2 * PSI + 0.1 * normalizedSize;

        scores.push_back({ pathway, score });
    }

    // Sort pathways by score in descending order
    sort(scores.begin(), scores.end(), [](const auto& a, const auto& b) {
        return a.second > b.second;
        });

    // Extract and return sorted pathway names
    vector<string> sortedPathways;
    for (const auto& score : scores) {
        sortedPathways.push_back(score.first);
    }

    return sortedPathways;
}


unordered_map<string, double> GeneSets::non_iterative_set_cover(int cluster_id, int GC) {
    cerr << "[DEBUG] Running PSC for cluster: " << cluster_id << endl;
    unordered_map<string, double> selected_pathways;
    cerr << "clusterToPathways size: " << this->clusterToPathways.size() << endl;

    auto cluster_pathways = this->clusterToPathways[cluster_id];
    cerr << "[DEBUG] Cluster " << cluster_id << " has " << cluster_pathways.size() << " pathways" << endl;

    // Get the universe set (all genes)
    flat_hash_set<uint64_t> uncovered_genes = get_universe_set(cluster_pathways);
    // Total genes to cover
    int total_genes = uncovered_genes.size();
    int genes_covered = 0;

    // auto sorted_pathways = this->influence_sort_pathway(cluster_pathways);
    auto sorted_pathways = this->combined_sort_frag_ppi_len(cluster_pathways);

    // Iterate through sorted pathways
    for (const auto& pathway : sorted_pathways) {
        // Calculate the intersection of the current pathway genes with the uncovered genes
        auto& pathway_genes = pathway_to_gene_set[pathway];
        flat_hash_set<uint64_t> common_genes = set_intersection(pathway_genes, uncovered_genes);

        // If the pathway covers any uncovered genes
        if (!common_genes.empty()) {
            // Add the pathway to the set cover
            double coverage = static_cast<double>(common_genes.size()) / total_genes * 100;
            selected_pathways[pathway] = coverage;

            // Remove the covered genes from the uncovered set
            for (const auto& gene : common_genes) {
                uncovered_genes.erase(gene);
            }

            // Update the count of covered genes
            genes_covered += common_genes.size();

            // If the required gene coverage is met, break the loop
            if ((uncovered_genes.size()) && static_cast<double>(genes_covered) / total_genes * 100 >= GC) {
                break;
            }
        }
    }

    return selected_pathways;
}



unordered_map<string, double> GeneSets::proportionalSetCover(int cluster_id, int GC) {
    cerr << "[DEBUG] Running PSC for cluster: " << cluster_id << endl;
    unordered_map<string, double> selected_pathways;
    cerr << "clusterToPathways size: " << this->clusterToPathways.size() << endl;

    auto cluster_pathways = this->clusterToPathways[cluster_id];
    cerr << "[DEBUG] Cluster " << cluster_id << " has " << cluster_pathways.size() << " pathways" << endl;

    for (auto& pathway : cluster_pathways)
        selected_pathways[pathway] = 0.0f;

    // Sort pathways by their PPI, descending.
    flat_hash_map<string, double> cluster_pathways_ppi;
    for (auto& pathway : cluster_pathways) {
        cluster_pathways_ppi[pathway] = this->pathway_to_average_PPI[pathway];
    }

    // ------------------------------------------
    // Sort the pathways by their PPI and fragmentation score
    // ------------------------------------------

    // std::vector<std::pair<std::string, double>> sorted_pathways_by_ppi = this->sort_desc_PPI(cluster_pathways_ppi);
    // auto sorted_pathways_by_ppi = this->sort_min_max(cluster_pathways);
    auto sorted_pathways_by_ppi = this->combine_and_sort(cluster_pathways);


    // Start set cover
    auto uncovered_genes = this->get_universe_set(cluster_pathways);
    int universal_size = uncovered_genes.size();
    int covered_genes = 0;
    int sorted_pathways_index = -1;
    while ((((double)covered_genes / universal_size) * 100 < GC) && !uncovered_genes.empty() && ++sorted_pathways_index < sorted_pathways_by_ppi.size() - 1) {
        cerr << "[DEBUG SETCOV] Covered genes %: " << ((double)covered_genes / universal_size) * 100 << "<" << GC << endl;
        cerr << "[DEBUG SETCOV] Uncovered genes: " << uncovered_genes.size() << endl;
        flat_hash_map<string, double> pathway_scores;
        flat_hash_map<string, int> pathway_to_common_genes;

        auto pathway = sorted_pathways_by_ppi[sorted_pathways_index];
        auto& pathway_genes = this->pathway_to_gene_set[pathway];

        // get genes intersection
        flat_hash_set<uint64_t> common_genes = this->set_intersection(pathway_genes, uncovered_genes);
        pathway_to_common_genes[pathway] = common_genes.size();
        cout << "[DEBUG SETCOV] pathway: " << pathway
            << " | avg_ppi: " << this->pathway_to_average_PPI[pathway]
            << " | fragmentation " << this->pathway_to_fragmentation[pathway]
            << " | common_genes " << common_genes.size()
            << " | Covered genes until now" << covered_genes
            << endl;


        // sometimes there are no common genes due to 
        // duplication full contamination from previous iteration
        if (common_genes.empty())
            continue;

        pathway_scores[pathway] = (double)common_genes.size() / pathway_genes.size();
        cerr << "pathway: " << pathway << " | basic: " << pathway_scores[pathway] << endl;
        // pathway_scores[pathway] += (double)1 / (abs(pathway_genes.size() - this->mean_pathway_length) * 10000);
        pathway_scores[pathway] += 1 - this->pathway_to_average_PSI[pathway]; // small PSI means better pathway.
        cerr << "pathway: " << pathway << " | basic + 1-PPI: " << pathway_scores[pathway] << endl;


        // Select the pathway with the highest overlap
        auto best_pathway = max_element(pathway_scores.begin(), pathway_scores.end(), [](auto& a, auto& b) {
            return a.second < b.second;
            });

        auto& best_pathway_genes = this->pathway_to_gene_set[best_pathway->first];
        covered_genes += pathway_to_common_genes[best_pathway->first];

        cout << "best_pathway: " << best_pathway->first << " score: " << best_pathway->second << endl;
        cout << "best_pathway_genes: " << best_pathway_genes.size() << endl;
        for (const auto& gene : best_pathway_genes)
            uncovered_genes.erase(gene); // TODO this is slow, use a set difference

        selected_pathways[best_pathway->first] = pathway_scores[best_pathway->first];
        cerr << "Now [DEBUG] Covered genes %: " << ((double)covered_genes / universal_size) * 100 << "<" << GC << endl;

    }
    return selected_pathways;
}

unordered_map<string, double> GeneSets::greedy_proportional_set_cover(int cluster_id, int GC) {
    cerr << "[DEBUG] Running PSC for cluster: " << cluster_id << endl;
    unordered_map<string, double> selected_pathways;
    cerr << "clusterToPathways size: " << this->clusterToPathways.size() << endl;

    auto cluster_pathways = this->clusterToPathways[cluster_id];
    cerr << "[DEBUG] Cluster " << cluster_id << " has " << cluster_pathways.size() << " pathways" << endl;

    for (auto& pathway : cluster_pathways)
        selected_pathways[pathway] = 0.0f;

    // Sort pathways by their PPI, descending.
    flat_hash_map<string, double> cluster_pathways_ppi;
    for (auto& pathway : cluster_pathways) {
        cluster_pathways_ppi[pathway] = this->pathway_to_average_PPI[pathway];
    }
    // std::vector<std::pair<std::string, double>> sorted_pathways_by_ppi = this->sort_desc_PPI(cluster_pathways_ppi);
    // auto sorted_pathways_by_ppi = this->sort_min_max(cluster_pathways);
    auto sorted_pathways_by_ppi = this->combine_and_sort(cluster_pathways);


    // Start set cover
    auto uncovered_genes = this->get_universe_set(cluster_pathways);
    int universal_size = uncovered_genes.size();
    int covered_genes = 0;

    while (((double)covered_genes / universal_size) * 100 < GC && !uncovered_genes.empty()) {
        cerr << "[DEBUG SETCOV] Covered genes %: " << ((double)covered_genes / universal_size) * 100 << "<" << GC << endl;
        cerr << "[DEBUG SETCOV] Uncovered genes: " << uncovered_genes.size() << endl;
        flat_hash_map<string, double> pathway_scores;
        flat_hash_map<string, int> pathway_to_common_genes;

        for (auto& pathway : sorted_pathways_by_ppi) {
            auto& pathway_genes = this->pathway_to_gene_set[pathway];

            // get genes intersection
            flat_hash_set<uint64_t> common_genes = this->set_intersection(pathway_genes, uncovered_genes);
            pathway_to_common_genes[pathway] = common_genes.size();
            cout << "[DEBUG SETCOV] pathway: " << pathway
                << " | avg_ppi: " << this->pathway_to_average_PPI[pathway]
                << " | common_genes " << common_genes.size()
                << " | Covered genes until now" << covered_genes
                << endl;


            // sometimes there are no common genes due to 
            // duplication full contamination from previous iteration
            if (common_genes.empty())
                continue;

            pathway_scores[pathway] = (double)common_genes.size() / pathway_genes.size();
            cerr << "pathway: " << pathway << " | basic: " << pathway_scores[pathway] << endl;
            // pathway_scores[pathway] += (double)1 / (abs(pathway_genes.size() - this->mean_pathway_length) * 10000);
            pathway_scores[pathway] += 1 - this->pathway_to_average_PSI[pathway]; // small PSI means better pathway.
            cerr << "pathway: " << pathway << " | basic + 1-PPI: " << pathway_scores[pathway] << endl;
        }

        // Select the pathway with the highest overlap
        auto best_pathway = max_element(pathway_scores.begin(), pathway_scores.end(), [](auto& a, auto& b) {
            return a.second < b.second;
            });

        auto& best_pathway_genes = this->pathway_to_gene_set[best_pathway->first];
        covered_genes += pathway_to_common_genes[best_pathway->first];

        cout << "best_pathway: " << best_pathway->first << " score: " << best_pathway->second << endl;
        cout << "best_pathway_genes: " << best_pathway_genes.size() << endl;
        for (const auto& gene : best_pathway_genes)
            uncovered_genes.erase(gene); // TODO this is slow, use a set difference

        selected_pathways[best_pathway->first] = pathway_scores[best_pathway->first];
        cerr << "Now [DEBUG] Covered genes %: " << ((double)covered_genes / universal_size) * 100 << "<" << GC << endl;

    }
    return selected_pathways;
}


/*
    ------------------------------------------
        Pathway Specificity Index PSI
    ------------------------------------------
*/

void GeneSets::build_all_pathways_PSI() {
    cerr << "[DEBUG] Building all pathways PSI" << endl;
    for (const auto& [pathway, _] : pathway_to_gene_set)
        pathway_to_average_PSI[pathway] = get_pathway_PSI(pathway);
}

void GeneSets::calculate_mean_pathway_length() {
    cerr << "[DEBUG] Calculating mean pathway length" << endl;
    double sum = 0;
    for (const auto& [pathway, genes] : pathway_to_gene_set)
        sum += genes.size();
    mean_pathway_length = sum / pathway_to_gene_set.size();
}

void GeneSets::build_gene_to_PSI() {
    cerr << "[DEBUG] Building gene to PSI map" << endl;
    auto it = gene_to_color->begin();
    while (it != gene_to_color->end()) {
        uint64_t gene = it.getHashedKmer();
        uint64_t color = it.getCount();
        uint64_t pathwayCount = color_to_ids[color].size();
        gene_to_associated_pathways[gene] = pathwayCount;

        // double psi = static_cast<double>(pathwayCount) / n_total_pathways;
        double psi = log2((double)pathwayCount / pathway_to_gene_set.size()) / log2(1.0 / pathway_to_gene_set.size());
        psi *= 100;

        gene_to_PSI[gene] = psi;
        it++;
    }
}

double GeneSets::get_pathway_PSI(string pathway) {
    // check if that pathway does not exist
    if (pathway_to_gene_set.find(pathway) == pathway_to_gene_set.end())
        throw std::invalid_argument("Pathway does not exist");

    auto genes = pathway_to_gene_set[pathway];
    // cout << "Pathway: " << pathway << " genes: " << genes.size() << endl;
    double average_psi = 0.0;
    for (auto& gene : genes) {
        average_psi += gene_to_PSI[gene];
    }
    return average_psi / genes.size();
}


/*
    ------------------------------------------
        Pathway Pleiotropy Index PPI
    ------------------------------------------
*/

void GeneSets::build_gene_to_no_clusters() {
    // iterate over clusters
    for (auto [cluster_id, pathways] : this->clusterToPathways) {
        auto universe_set = get_universe_set(pathways);
        for (auto& gene : universe_set) {
            this->gene_to_no_clusters[gene]++;
        }
    }
}

void GeneSets::build_gene_to_PPI() {
    cerr << "[DEBUG] Building gene to PPI map" << endl;
    for (auto& [gene, no_clusters] : gene_to_no_clusters) {
        
        // TODO: relocate later 
        double pcsi = log2((double)no_clusters / this->clusterToPathways.size()) / log2(1.0 / this->clusterToPathways.size());
        pcsi *= 100;
        gene_to_PCSI[gene] = pcsi;

        double ppi = static_cast<double>(no_clusters) / this->gene_to_associated_pathways[gene];
        ppi *= 100;
        gene_to_PPI[gene] = ppi;
    }
}

void GeneSets::build_pathway_to_average_PPI(string pathway) {
    // check if that pathway does not exist
    if (pathway_to_gene_set.find(pathway) == pathway_to_gene_set.end())
        throw std::invalid_argument("Pathway does not exist");

    auto& genes = pathway_to_gene_set[pathway];
    double average_ppi = 0.0;
    // TODO relocate later
    double average_pcsi = 0.0;
    for (auto& gene : genes) {
        average_ppi += gene_to_PPI[gene];
        average_pcsi += gene_to_PCSI[gene];
    }
    pathway_to_average_PPI[pathway] = average_ppi / genes.size();
    pathway_to_average_PCSI[pathway] = average_pcsi / genes.size();
}

double GeneSets::get_cluster_average_ppi(uint32_t cluster_id) {
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

unordered_map<string, double> GeneSets::get_pathways_pcsi(){
    auto& before_convert = this->pathway_to_average_PCSI;
    unordered_map<string, double> after_convert;
    for (auto& [pathway, pcsi] : before_convert) {
        after_convert[pathway] = pcsi;
    }
    return after_convert;
}

    void GeneSets::export_genes_to_ppi_psi_tsv(string filename){
        // TODO: NEW: adding PCSI
        ofstream myfile;
        myfile.open(filename);
        myfile << "gene\tppi\tpsi\tpcsi\n";
        for (auto& [gene, ppi] : gene_to_PPI) {
            myfile << this->hashed_gene_to_name[gene] << "\t" << ppi << "\t" << gene_to_PSI[gene] << "\t" << gene_to_PCSI[gene] << "\n";
        }
        myfile.close();
    }

// Destructor
GeneSets::~GeneSets() {
}

