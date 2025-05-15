#include <kDataFrame.hpp>
#include <colored_kDataFrame.hpp>
#include <kSpider.hpp>
#include <DBRetina.hpp>


using int_vec_map = parallel_flat_hash_map<uint32_t, vector<uint32_t>, std::hash<uint32_t>, std::equal_to<uint32_t>, std::allocator<std::pair<const uint32_t, vector<uint32_t>>>, 1>;


inline void load_namesMap(string filename, phmap::flat_hash_map<int, std::string>& map) {
    std::ifstream inputFile(filename);
    if (!inputFile.is_open()) {
        std::cerr << "Error opening the file: " << filename << std::endl;
        return;
    }

    // open file for writing

    std::string line;
    std::getline(inputFile, line); // skip first line
    while (std::getline(inputFile, line)) {
        std::istringstream lineStream(line);
        std::string column1, column2;

        if (std::getline(lineStream, column1, '|') && std::getline(lineStream, column2, '|')) {
            cout << column1 << "|" << column2 << std::endl;
            transform(column2.begin(), column2.end(), column2.begin(), ::tolower);
            map.operator[](stoi(column1)) = column2;
        }
        else {
            std::cerr << "Invalid line format: " << line << std::endl;
            inputFile.close();
            return;
        }
    }
    inputFile.close();
}


inline void set_to_vector(const phmap::flat_hash_set<uint32_t>& set, vector<uint32_t>& vec) {
    vec.clear();
    vec.reserve(set.size());
    for (auto& i : set) {
        vec.push_back(i);
    }
}

inline parallel_flat_hash_set<string> load_query_file_single_column(string filename) {
    parallel_flat_hash_set<string> queries;
    ifstream query_file(filename);
    string line;
    while (getline(query_file, line)) {
        // lowecase line
        transform(line.begin(), line.end(), line.begin(), ::tolower);
        // remove double quotes
        line.erase(remove(line.begin(), line.end(), '"'), line.end());
        queries.emplace(line);
    }
    return queries;
}

inline void load_colors_to_matched_supergroups(const std::string& filename, int_vec_map* map)
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
    cerr << "Loaded colors: " << map->size() << endl;
}

inline std::string join(std::vector<std::string>& strings, std::string delim)
{
    return std::accumulate(strings.begin(), strings.end(), std::string(),
        [&delim](std::string& x, std::string& y) {
            return x.empty() ? y : x + delim + y;
        });
}

void query(string index_prefix, string inverted_index_prefix, string query_file, string output_prefix, string commands) {

    // load kDataFrame
    auto* inverted_kf = kDataFramePHMAP::load(inverted_index_prefix);
    int_vec_map inverted_color_to_ids;
    string inverted_colors_map_file = inverted_index_prefix + "_color_to_sources.bin";
    load_colors_to_matched_supergroups(inverted_colors_map_file, &inverted_color_to_ids);
    string inverted_namesmap_file = inverted_index_prefix + ".namesMap";
    phmap::flat_hash_map<int, std::string> inverted_namesmap;
    load_namesMap(inverted_namesmap_file, inverted_namesmap);


    // naming
    string key_val_suffix, val_key_suffix;
    key_val_suffix = "_feature_to_groups.tsv";
    val_key_suffix = "_features_count_per_group.tsv";


    // load query file
    auto queries = load_query_file_single_column(query_file);
    cout << "number of queries: " << queries.size() << endl;

    // print queries
    auto hasher = string_hasher();
    parallel_flat_hash_map<string, vector<string>> results;
    parallel_flat_hash_map<string, uint32_t> inverted_results;
    parallel_flat_hash_map<string, vector<string>> feature_to_groups;


    flat_hash_map<string, bool> super_groups_of_interest;


    ofstream outfile2;
    outfile2.open(output_prefix + key_val_suffix);
    outfile2 << commands << endl;
    outfile2 << "feature\tsupergroups_count\tsupergroups" << endl;

    /*
        input(gene) -> color
        color -> supergroup_ids
        supergroup_ids -> supergroup_names

    */

    for (auto single_query : queries) {
        super_groups_of_interest[single_query] = true;
        uint64_t hashed_q = hasher(single_query);
        auto color = inverted_kf->getCount(hashed_q);
        auto matched_supergroups = inverted_color_to_ids[color];
        assert(matched_supergroups.size() > 0);

        if (matched_supergroups.size() == 0) cerr << "no matched_supergroups found for query: " + single_query;


        for (auto supergroup_id : matched_supergroups) {
            // only if found in namesmap
            if (inverted_namesmap.find(supergroup_id) == inverted_namesmap.end()) {
                throw std::runtime_error("[DEV BUG] supergroup_id not found in namesmap: [" + to_string(supergroup_id) + '\n');
            }
            string supergroup_name = inverted_namesmap[supergroup_id];
            // results[single_query].push_back(supergroup_name);
            inverted_results[supergroup_name]++;
        }
    }

    // free inverted index
    cout << "freeing inverted index" << endl;
    delete inverted_kf;
    inverted_color_to_ids.clear();
    inverted_namesmap.clear();



    // load kDataFrame
    auto* kf = kDataFramePHMAP::load(index_prefix);
    int_vec_map color_to_ids;
    string colors_map_file = index_prefix + "_color_to_sources.bin";
    load_colors_to_matched_supergroups(colors_map_file, &color_to_ids);
    string namesmap_file = index_prefix + ".namesMap";
    phmap::flat_hash_map<int, std::string> namesmap;
    load_namesMap(namesmap_file, namesmap);


    // Export to TSV
    /*
    ofstream outfile;
    outfile.open(output_prefix + key_val_suffix);
    outfile << commands << endl;
    for (auto r : results) {
        outfile << r.first << "\t";
        outfile << join(r.second, "|");
        outfile << endl;
    }
    outfile.close();
    */

    // sort inverted_results by groups_count:
    std::vector<std::pair<string, int>> ordered_inverted_results(inverted_results.begin(), inverted_results.end());

    std::sort(ordered_inverted_results.begin(), ordered_inverted_results.end(),
        [](const std::pair<string, int>& a, const std::pair<string, int>& b) {
            return a.second > b.second;
        });


    for (auto& [feature_name, groups_count] : ordered_inverted_results) {
        outfile2 << feature_name << "\t" << groups_count << '\t';

        uint64_t hashed_q = hasher(feature_name);
        auto color = kf->getCount(hashed_q);
        auto matched_supergroups = color_to_ids[color];
        assert(matched_supergroups.size() > 0);
        string delimiter = "";
        for (int i = 0; i < matched_supergroups.size(); i++) {
            auto supergroup_id = matched_supergroups[i];
            string supergroup_name = namesmap[supergroup_id];
            if (super_groups_of_interest[supergroup_name]) {
                outfile2 << delimiter;
                outfile2 << supergroup_name;
                delimiter = "|";
            }
            // delimiter = "|";
        }
        // auto supergroup_id = matched_supergroups[matched_supergroups.size() - 1];
        // if (super_groups_of_interest[namesmap[supergroup_id]])
        //     outfile2 << namesmap[supergroup_id];
        outfile2 << endl;
    }
    outfile2.close();

}