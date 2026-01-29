#include <iostream>
#include <cstdint>
#include <string>
#include "argh.h"
#include <chrono>
#include "parallel_hashmap/phmap.h"
#include <queue>

using std::string;

namespace kSpider {

    void pairwise(string index_prefix, int user_threads, string cutoff_distance_type, double cutoff_threshold, string full_command, bool calculate_pvalue);
    void dbretina_indexing(string json_file, string user_index_prefix);
};
