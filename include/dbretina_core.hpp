#include <iostream>
#include <cstdint>
#include <string>
#include <chrono>
#include "parallel_hashmap/phmap.h"
#include <queue>

using std::string;

namespace dbretina {

    void pairwise(string index_prefix, int user_threads, string cutoff_distance_type, double cutoff_threshold, string full_command, bool calculate_pvalue);
    void dbretina_indexing(string json_file, string user_index_prefix);
    void dbretina_append(string existing_dbri_path, string new_json_file, string new_raw_json_file, string output_dbri_path);
    void dbretina_merge(string index_a_path, string index_b_path, string output_dbri_path);
};
