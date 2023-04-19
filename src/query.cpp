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

    std::string line;
    std::getline(inputFile, line); // skip first line
    while (std::getline(inputFile, line)) {
        std::istringstream lineStream(line);
        std::string column1, column2;

        if (std::getline(lineStream, column1, '|') && std::getline(lineStream, column2, '|')) {
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
    cerr << "Loaded colors: " << map->size() << endl;
}

inline std::string join(std::vector<std::string>& strings, std::string delim)
{
    return std::accumulate(strings.begin(), strings.end(), std::string(),
        [&delim](std::string& x, std::string& y) {
            return x.empty() ? y : x + delim + y;
        });
}

void query(string index_prefix, string query_file, string output_prefix, bool inverted) {
    int_vec_map color_to_ids;
    string colors_map_file = index_prefix + "_color_to_sources.bin";
    load_colors_to_sources(colors_map_file, &color_to_ids);

    string namesmap_file = index_prefix + ".namesMap";
    phmap::flat_hash_map<int, std::string> namesmap;
    load_namesMap(namesmap_file, namesmap);

    // naming
    string key_val_suffix, val_key_suffix;
    if (!inverted) {
        key_val_suffix = "_group_to_genes.tsv";
        val_key_suffix = "_gene_to_groupsCount.tsv";
    }
    else {
        key_val_suffix = "_gene_to_groups.tsv";
        val_key_suffix = "_group_to_genesCount.tsv";
    }


    // load kDataFrame
    auto* kf = kDataFramePHMAP::load(index_prefix);

    // load query file
    auto queries = load_query_file_single_column(query_file);
    cout << "number of queries: " << queries.size() << endl;

    // print queries
    auto hasher = string_hasher();
    parallel_flat_hash_map<string, vector<string>> results;
    parallel_flat_hash_map<string, uint32_t> inverted_results;

    for (auto q : queries) {
        uint64_t hashed_q = hasher(q);
        auto color = kf->getCount(hashed_q);
        auto sources = color_to_ids[color];
        for (auto s : sources) {
            results[q].push_back(namesmap[s]);
            inverted_results[namesmap[s]]++;
        }
    }

    // Export to TSV
    ofstream outfile;
    outfile.open(output_prefix + key_val_suffix + ".tsv");
    for (auto r : results) {
        outfile << r.first << "\t";
        outfile << join(r.second, "|");
        outfile << endl;
    }
    outfile.close();



    ofstream outfile2;
    outfile2.open(output_prefix + val_key_suffix + ".tsv");
    for (auto r : inverted_results) {
        outfile2 << r.first << "\t" << r.second << endl;
    }
    outfile2.close();

}