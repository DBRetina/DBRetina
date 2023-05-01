#include <iostream>
#include <string>
#include "parallel_hashmap/phmap.h"
#include <sstream>
#include <fstream>
#include "parallel_hashmap/phmap_dump.h"
#include <sstream>
#include "cpp-json/json.h"
#include "zstr.hpp"

using phmap::parallel_flat_hash_map;
using phmap::parallel_flat_hash_set;
using namespace std;

using str_vec_map = parallel_flat_hash_map<string, parallel_flat_hash_set<string>>;
using str_hashed_vec_map = parallel_flat_hash_map<string, parallel_flat_hash_set<uint64_t>>;
using str_str_map = parallel_flat_hash_map<string, string>;
using hashed_MAP = parallel_flat_hash_map<string, parallel_flat_hash_set<uint64_t>>;


struct string_hasher
{
    std::hash<string> hasher;

    // initialize the hasher
    string_hasher(): hasher() {}

    // overload the () operator
    size_t operator()(const string& str) const
    {
        return hasher(str);
    }
};

void load_tsv_to_map(string filename, str_vec_map* map);
void load_tsv_to_map_no_names(string filename, str_vec_map* map);
void load_names_tsv_to_map(string filename, str_str_map* map);
void parse_dbretina_json(string json_file, str_hashed_vec_map* map);
void sketch_dbretina(string asc_file, string names_file, string user_prefix = "NA");
void query(string index_prefix, string inverted_index_prefix, string query_file, string output_prefix, string commands);

// class DBRetina_json_parser {

// private:
//     string_hasher hasher;

// public:
//     string json_file_name;
//     zstr::ifstream sig_stream;
//     json::value json;
//     string output_prefix;
//     str_hashed_vec_map* map = new str_hashed_vec_map();

//     void load_json(string json_file);

//     void export_map_to_tsv(string output_file);

//     DBRetina_json_parser(string json_file);

//     ~DBRetina_json_parser(){
//         delete map;
//     };

// };