#include <iostream>
#include <cstdint>
#include <algorithm>
#include <chrono>
#include <boost/algorithm/string/join.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/math/distributions/hypergeometric.hpp>
#include <boost/math/special_functions/binomial.hpp>
#include <boost/functional/hash.hpp>
#include <ctime>
#include<omp.h>
#include "parallel_hashmap/phmap.h"
#include "parallel_hashmap/phmap_dump.h"
#include <cassert>
#include <math.h>
#include "DBRetinaIndex.hpp"
#include "DBRetinaPairwise.hpp"
#include "ParquetPairwiseWriter.hpp"

using boost::adaptors::transformed;
using boost::algorithm::join;
using namespace std;
using namespace phmap;

// using Map = parallel_flat_hash_map<std::pair<uint32_t, uint32_t>, std::uint64_t, boost::hash<pair<uint32_t, uint32_t>>, std::equal_to<std::pair<uint32_t, uint32_t>>, std::allocator<std::pair<const std::pair<uint32_t, uint32_t>, uint32_t>>, 12, std::mutex>;
using int_int_map = parallel_flat_hash_map<uint32_t, uint32_t, std::hash<uint32_t>, std::equal_to<uint32_t>, std::allocator<std::pair<const uint32_t, uint32_t>>, 1>;
using int_vec_map = parallel_flat_hash_map<uint32_t, vector<uint32_t>, std::hash<uint32_t>, std::equal_to<uint32_t>, std::allocator<std::pair<const uint32_t, vector<uint32_t>>>, 1>;

using PAIRS_COUNTER = phmap::parallel_flat_hash_map<
    std::pair<uint32_t, uint32_t>,
    std::uint64_t,
    boost::hash<pair<uint32_t, uint32_t>>,
    std::equal_to<std::pair<uint32_t, uint32_t>>,
    std::allocator<std::pair<const std::pair<uint32_t, uint32_t>, uint64_t>>, 12, std::mutex>;

using BINS_FEATURE_COUNT = phmap::parallel_flat_hash_map<
    std::string, uint32_t,
    phmap::priv::hash_default_hash<std::string>,
    phmap::priv::hash_default_eq<std::string>,
    std::allocator<std::pair<const std::string, uint32_t>>,
    1,
    std::mutex>;

typedef std::chrono::high_resolution_clock Time;

class Stats {
private:
    flat_hash_map<string, flat_hash_map<string, uint64_t>> stats;

public:
    double min_odds_ratio = 1.0;
    double max_odds_ratio = 1.0;
    Stats() {
        vector<string> distances = { "containment", "ochiai", "jaccard", "csi", "dice" };
        for (string& distance : distances) {
            stats[distance] = flat_hash_map<string, uint64_t>();
            for (int value = 0; value < 100; value += 5) {
                string range = to_string(value) + "-" + to_string(value + 5);
                this->stats[distance][range] = 0;
            }
            this->stats[distance]["100-100"] = 0;
        }
    }

    void update_odds_ratio_stats(double odds_ratio) {
        auto& min = this->min_odds_ratio;
        auto& max = this->max_odds_ratio;
        if (odds_ratio < min) min = odds_ratio;
        if (odds_ratio > max) max = odds_ratio;
    }

    void print_stats_in_json_format() {
        cout << "{" << endl;
        for (auto& [range, stats] : this->stats) {
            cout << "\t\"" << range << "\": {" << endl;
            for (auto& [stat_name, stat_value] : stats) {
                cout << "\t\t\"" << stat_name << "\": " << stat_value << "," << endl;
            }
            cout << "\t}," << endl;
        }
        cout << "}" << endl;
    }

    void stats_to_json_file(string filename) {
        ofstream myfile;
        myfile.open(filename);
        myfile << "{" << endl;
        int json_end_counter = this->stats.size();
        for (auto& [range, stats] : this->stats) {
            json_end_counter--;
            myfile << "\t\"" << range << "\": {" << endl;

            // iterate over stats by index
            for (int i = 0; i < stats.size(); i++) {
                auto& [stat_name, stat_value] = *std::next(stats.begin(), i);
                myfile << "\t\t\"" << stat_name << "\": " << stat_value;
                if (i < stats.size() - 1) {
                    myfile << ",";
                }
                myfile << endl;
            }
            if (json_end_counter > 0) {
                myfile << "\t}," << endl;
            }
            else {
                myfile << "\t}" << endl;
            }
        }
        myfile << "}" << endl;
        myfile.close();
        string new_file_name_without_extension = filename.substr(0, filename.find_last_of("."));
        this->write_odds_ratio_in_file(new_file_name_without_extension + "_odds_ratio.txt");
    }

    void write_odds_ratio_in_file(string filename) {
        ofstream myfile;
        myfile.open(filename);
        myfile << this->min_odds_ratio << endl;
        myfile << this->max_odds_ratio << endl;
        myfile.close();
    }

    string map_value_to_range(double& value) {
        int lower = static_cast<int>(std::floor(value / 5)) * 5;
        int upper = lower + 5;
        if (upper > 100) upper = 100;
        return std::to_string(lower) + "-" + std::to_string(upper);
    }



    void add_stat(string stat_name, double value) {
        string range = this->map_value_to_range(value);
        this->stats[stat_name][range] += 1;
    }


    ~Stats() {
        this->stats.clear();
    }


};

template <typename T>
inline void ascending(T& dFirst, T& dSecond)
{
    if (dFirst > dSecond)
        std::swap(dFirst, dSecond);
}


inline void map_insert(int_int_map& _MAP, uint32_t& key, uint32_t& value) {
    _MAP.insert(make_pair(key, value));
}


inline void load_namesMap(string filename, phmap::flat_hash_map<int, std::string>& map) {
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
            map.operator[](stoi(column1)) = column2;
        }
        else {
            inputFile.close();
            throw std::runtime_error("Invalid line format: '" + line + "'");
        }
    }

    inputFile.close();
}

// NOTE: Population size is now correctly loaded from .dbri metadata
// via dbri.get_population_size() at runtime. See line ~430.


namespace dbretina {

    void set_to_vector(const phmap::flat_hash_set<uint32_t>& set, vector<uint32_t>& vec) {
        vec.clear();
        vec.reserve(set.size());
        for (auto& i : set) {
            vec.push_back(i);
        }
    }

    void load_colors_to_sources(const std::string& filename, int_vec_map* map)
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

    void load_colors_count(const std::string& filename, int_int_map& map) {
        flat_hash_map<uint64_t, uint64_t> tmpMap;
        phmap::BinaryInputArchive ar_in_colorsCount(filename.c_str());
        tmpMap.phmap_load(ar_in_colorsCount);
        assert(tmpMap.size());
        for (auto& i : tmpMap) {
            map.insert_or_assign(i.first, i.second);
        }
    }


    double calc_foldChange(int k, int s, int M, int N) {
        return (double)k * N / (double)(s * M);
    }

    /*
        Guarded hypergeometric CDF evaluation.

        boost::math::hypergeometric_distribution<>(r, n, N) requires the
        distribution parameters to satisfy 0 <= r <= N and 0 <= n <= N, and the
        CDF argument x to lie within the support [max(0, n+r-N), min(n, r)].
        Outside that support boost throws std::domain_error which, under the
        default error policy, calls std::terminate -> SIGABRT (a core dump).
        This happens for small/degenerate gene populations where, e.g., the two
        groups together cover more than the whole population (n + r - N > 0) so
        the observed-overlap minus one falls below the lower support bound.

        We guard both the parameters and the argument here:
          - degenerate distribution parameters (out of the 0..N range, or N<1)
            -> return a non-significant result (cdf treated as 0.0 so the
               over-enrichment pvalue 1 - cdf becomes 1.0);
          - the CDF argument is clamped to the support: below the support the
            cumulative probability is exactly 0.0 and at/above the top it is
            exactly 1.0, which is the mathematically correct continuation of the
            CDF. Inside the support the value is computed by boost unchanged, so
            valid-domain pvalues are bit-for-bit identical to before.
    */
    double safeHyperCDF(int x, int r, int n, int N) {
        // Degenerate / out-of-domain distribution parameters: treat as the
        // empty lower tail (cdf = 0) so callers fall back to a pvalue of 1.0.
        if (N < 1 || r < 0 || r > N || n < 0 || n > N) {
            return 0.0;
        }
        int lower = std::max(0, n + r - N);   // smallest possible overlap
        int upper = std::min(n, r);           // largest possible overlap
        if (x < lower) return 0.0;            // entire mass is above x
        if (x >= upper) return 1.0;           // entire mass is at or below x
        boost::math::hypergeometric_distribution<> hg(r, n, N);
        return boost::math::cdf(hg, x);
    }

    double calcPValue(int k, int s, int M, int N, bool isOverEnrichment) {
        /*
            Here we set the isOverEnriched null hypothesis as boolean set by user.
            If the user set the isOverEnriched to True, then the null hypothesis is that the gene is over-enriched in the source.
            The null hypothesis is a hypothesis that says there is no statistical significance between the two variables in the hypothesis.
            So if it's true, then the alternative hypothesis is that the gene is under-enriched in the source.
            That means with low p-value we can reject the null hypothesis and say that the gene is under-enriched in the source.
            If the user set the isOverEnriched to False, then the null hypothesis is that the gene is under-enriched in the source.
            The alternative hypothesis is that the gene is over-enriched in the source.
            That means with low p-value we can reject the null hypothesis and say that the gene is over-enriched in the source.
        */
        if (isOverEnrichment) {
            return 1 - safeHyperCDF(k - 1, M, s, N);
        }
        return safeHyperCDF(k, M, s, N);
    }

    // Consireding the isOverEnrichment = True
    // (Null hypothesis is that the gene is over-enriched in the source)
    double fastHyperPValue(int k, int s, int M, int N) {
        /*
            The original null hypothesis is that the gene is over-enriched in the source.
            Explanation of the Pvalue calculations between two gene sets with shared genes:
                - The null hypothesis is that the gene is over-enriched in the source.
                - The alternative hypothesis is that the gene is under-enriched in the source.
                - The probability of getting k or more successes in s trials is 1 - the probability of getting k or less successes in s trials.
                - The probability of getting k or less successes in s trials is the cumulative distribution function of the hypergeometric distribution.

            So, if the pvalue is less than the significance level (alpha), we reject the null hypothesis and accept the alternative hypothesis.
            If the pvalue is greater than the significance level (alpha), we accept the null hypothesis and reject the alternative hypothesis.
            smaller pvalues here means that the gene is more over-enriched in the source which means that the gene is more important in the source.

            safeHyperCDF guards boost's domain constraints so small/degenerate
            populations return a non-significant pvalue (1.0) instead of aborting.
        */
        return 1 - safeHyperCDF(k - 1, M, s, N);
    }

    // Disabled for now
    /*
    double fisher_exact(int k_shared_features, int source_1_features, int source_2_features, int population_size) {
        int a = k_shared_features;
        int b = source_1_features - k_shared_features;
        int c = source_2_features - k_shared_features;
        int d = population_size - source_1_features - c;

        int minval = std::min(a + b, a + c);
        int maxval = std::max(0, a - d);
        double p = 0.0;

        boost::math::hypergeometric_distribution<> hgd(a + c, a + b, a + b + c + d);

        for (int i = maxval; i <= minval; ++i) {
            p += boost::math::pdf(hgd, i);
        }
        return p;
    }

    std::tuple<double, double, double> enrichmentAnalysis(int k, int s, int M, int N, bool isOverEnrichment) {
        double fold_change = calc_foldChange(k, s, M, N);
        double expectedSuccesses = calcExpectedSuccesses(s, M, N);
        double pvalue = calcPValue(k, s, M, N, isOverEnrichment);
        return std::make_tuple(pvalue, expectedSuccesses, fold_change);
    }
*/


    double odds_ratio(int k, int s, int M, int N) {
        /*

       Description:
        The odds_ratio function calculates the odds ratio between two gene sets. 
        This measure of association compares the odds of shared genes between two sets. 
        The result provides an idea of whether the overlap between the two gene sets- 
        is due to chance or there is a significant association between them.

        Parameters:
            - `int k`: Represents the number of shared genes between the two sets.
            - `int s`: Total number of genes in the first gene set.
            - `int M`: Total number of genes in the second gene set.
            - `int N`: Total number of genes in the universe.

        The values `a`, `b`, `c`, and `d` correspond to the entries in a 2x2 contingency table as shown below:

        |------------------|-----------------|------------------|
        |                  | Second gene set | Not in second set|
        |------------------|-----------------|------------------|
        | In first set     | a = k           | b = s-k          |
        |------------------|-----------------|------------------|
        | Not in first set | c = M-k         | d = N-(s+M)+k    |
        |------------------|-----------------|------------------|

        The function calculates the odds ratio as `(a * d) / (b * c)`.

        Exceptions:
        The function returns `-1` when `b` or `c` is `0` to prevent division by zero.
        All comparisons with -1 will be removed in any postprocessing step.

        Returns:
        The function returns a `double` that is the calculated odds ratio. 
        If the denominator is zero, it returns `-1`.
        */

        int a, b, c, d;
        a = k;
        b = s - k;
        c = M - k;
        d = N - (s + M) + k;

        // Check for division by zero
        if (b == 0 || c == 0) {
            return -1;
            // throw std::invalid_argument("Odds_ratio Denominator cannot be zero");
        }

        return static_cast<double>(a * d) / (b * c);
    }


    void pairwise(string index_prefix, int user_threads, string cutoff_distance_type, double cutoff_threshold, string full_command, bool calculate_pvalue) {

        vector<string> allowed_distances = { "containment", "ochiai", "jaccard" };
        // cutoff_distance_type must be in allowed_distances
        if (std::find(allowed_distances.begin(), allowed_distances.end(), cutoff_distance_type) == allowed_distances.end()) {
            throw std::invalid_argument("cutoff_distance_type must be in " + string("{containment, ochiai, jaccard}"));
        }

        // Load from .dbri unified index
        std::string dbri_path = index_prefix + ".dbri";
        auto dbri = DBRetinaIndex::open(dbri_path);

        auto begin_time = Time::now();

        // Load population size from metadata (fixes hardcoded 44260 bug)
        uint64_t population_size = dbri.get_population_size();
        cout << "population_size: " << population_size << endl;

        // Load color_to_sources
        int_vec_map color_to_ids;
        dbri.load_color_to_sources(color_to_ids);
        cout << "[dev] mapping colors to groups: " << std::chrono::duration<double, std::milli>(Time::now() - begin_time).count() / 1000 << " secs" << endl;

        // Load namesMap
        flat_hash_map<int, std::string> namesMap;
        dbri.load_names_map(namesMap);
        assert(namesMap.size());

        // Load colorsCount
        begin_time = Time::now();
        int_int_map colorsCount;
        dbri.load_color_count(colorsCount);
        cout << "[dev] parsing index colors: " << std::chrono::duration<double, std::milli>(Time::now() - begin_time).count() / 1000 << " secs" << endl;

        // Load feature counts
        begin_time = Time::now();
        flat_hash_map<uint32_t, uint32_t> groupID_to_featureCount;
        dbri.load_group_feature_count(groupID_to_featureCount);
        assert(groupID_to_featureCount.size());

        std::ofstream fstream_featureCount;
        fstream_featureCount.open(index_prefix + "_DBRetina_featuresNo.tsv");
        fstream_featureCount << "ID\tgroup\tfeatures\n";
        uint64_t counter = 0;
        for (const auto& item : groupID_to_featureCount) {
            fstream_featureCount << ++counter << '\t' << item.first << '\t' << item.second << '\n';
        }
        fstream_featureCount.close();
        cout << "[dev] features counting: " << std::chrono::duration<double, std::milli>(Time::now() - begin_time).count() / 1000 << " secs" << endl;

        // Loading done

        begin_time = Time::now();
        clock_t begin_detailed_pairwise_comb, begin_detailed_pairwise_edges, begin_detailed_pairwise_edges_insertion;
        double detailed_pairwise_comb = 0.0;
        double detailed_pairwise_edges = 0.0;
        double detailed_pairwise_edges_insertion = 0.0;

        PAIRS_COUNTER edges;

        // convert map to vec for parallelization purposes.
        auto vec_color_to_ids = std::vector<std::pair<uint32_t, vector<uint32_t>>>(color_to_ids.begin(), color_to_ids.end());

        cerr << "[dev] number of colors = " << vec_color_to_ids.size() << endl;

        double average_color_size = 0.0;
        for (auto const& item : vec_color_to_ids) {
            average_color_size += item.second.size();
        }
        average_color_size /= vec_color_to_ids.size();

        cerr << "[dev] average color size = " << (int)average_color_size << endl;

        // Longest-processing-time-first: per-color work is O(color_size^2), so
        // dispatch the giant colors before the tail. `edges` is an order-independent
        // reduction, so reordering colors does not change the output. Sort is
        // O(n log n) on the color count -- negligible vs the O(sum k^2) loop.
        std::sort(vec_color_to_ids.begin(), vec_color_to_ids.end(),
            [](const std::pair<uint32_t, vector<uint32_t>>& a,
                const std::pair<uint32_t, vector<uint32_t>>& b) {
                    return a.second.size() > b.second.size();
            });

        size_t n = vec_color_to_ids.size();

        omp_set_num_threads(user_threads);
        begin_time = Time::now();

        // Dynamic scheduling balances the heavy-tailed per-color load: big colors
        // dominate the C(k,2) work and would otherwise pile onto one static block.
#pragma omp parallel for schedule(dynamic)
        for (size_t vec_i = 0; vec_i < n; ++vec_i) {
            const auto& item = vec_color_to_ids[vec_i];
            const auto& ids = item.second;
            // Read ccount once per color (colorsCount is fully loaded before this
            // parallel region and never written here, so concurrent find() is
            // data-race-free; the old mutating operator[] could insert/rehash).
            const auto _cc_it = colorsCount.find(item.first);
            uint32_t ccount = (_cc_it != colorsCount.end()) ? _cc_it->second : 0;

            // Generate the C(k,2) unordered pairs inline straight into `edges`
            // (no Combo / no materialized combs vector). This nested a<b loop
            // enumerates the identical pair set Combo did.
            for (size_t a = 0; a + 1 < ids.size(); ++a) {
                for (size_t b = a + 1; b < ids.size(); ++b) {
                    uint32_t _seq1 = ids[a];
                    uint32_t _seq2 = ids[b];
                    ascending(_seq1, _seq2);

                    auto _p = make_pair(_seq1, _seq2);
                    edges.try_emplace_l(_p,
                        [ccount](PAIRS_COUNTER::value_type& v) { v.second += ccount; }, // called only when key was already present
                        ccount
                    );
                }
            }
        }

        cout << "[dev] pairwise hashmap construction: " << std::chrono::duration<double, std::milli>(Time::now() - begin_time).count() / 1000 << " secs" << endl;
        cout << "[dev] Number of pairwise comparisons: " << edges.size() << endl;
        cout << "[dev] writing pairwise matrix to " << index_prefix << "_DBRetina_pairwise.tsv | Please wait..." << endl;

        Stats distances_stats;

        auto formatDouble = [](double val) {
            char buffer[20];
            std::snprintf(buffer, sizeof(buffer), "%.1f", val);
            return std::string(buffer);
            };

        // --- TSV output (kept for migration) ---
        std::ofstream myfile;
        myfile.open(index_prefix + "_DBRetina_pairwise.tsv");
        myfile << "# DBRetina pairwise output\n";
        myfile << "# population_size: " << population_size << '\n';
        myfile << "# NOTE: This file contains raw (uncorrected) p-values. FDR-corrected results are in *_fdr.tsv if --fdr was used.\n";
        myfile << "# containment = shared/min(s1,s2)\n";
        myfile << "#nodes:" << namesMap.size() << '\n';
        myfile << "#command: " << full_command << '\n';

        myfile
            << "group_1_ID"
            << "\tgroup_2_ID"
            << "\tgroup_1_name"
            << "\tgroup_2_name"
            << "\tshared_features"
            << "\tcontainment"
            << "\tochiai"
            << "\tjaccard"
            << "\tcsi"
            << "\tdice"
            << "\todds_ratio";
        if (calculate_pvalue) { myfile << "\tpvalue"; }

        myfile << '\n';

        // --- .dbrp binary output ---
        DBRetinaPairwise pw;
        uint8_t dbrp_flags = 0x3F;  // all metrics except pvalue
        if (calculate_pvalue) dbrp_flags |= 0x40;

        // Build metadata JSON for .dbrp
        std::string dbrp_metadata = "{\"population_size\":" + std::to_string(population_size)
            + ",\"cutoff_metric\":\"" + cutoff_distance_type + "\""
            + ",\"cutoff_threshold\":" + std::to_string(cutoff_threshold)
            + ",\"command\":\"" + full_command + "\""
            + ",\"num_groups\":" + std::to_string(namesMap.size())
            + "}";

        pw.begin_write(index_prefix + "_DBRetina_pairwise.dbrp", dbrp_flags, namesMap, dbrp_metadata);

        // Stats tracking for .dbrp
        PairwiseStatistics dbrp_stats;
        dbrp_stats.min_odds_ratio = 1.0;
        dbrp_stats.max_odds_ratio = 1.0;
        // Initialize histograms for 5 metrics (containment, ochiai, jaccard, csi, dice)
        for (uint8_t mid = 0; mid < 5; mid++) {
            MetricHistogram hist;
            hist.metric_id = mid;
            hist.bucket_counts.resize(21, 0);
            dbrp_stats.histograms.push_back(hist);
        }

        uint64_t line_count = 0;

        // --- Parquet output (parallel) ---
        std::string parquet_dir = index_prefix + "_DBRetina_pairwise";
        ParquetPairwiseWriter parquet_writer(parquet_dir, calculate_pvalue, user_threads, namesMap);

        // Convert edges to vector for parallel iteration
        auto vec_edges = std::vector<std::pair<std::pair<uint32_t, uint32_t>, uint64_t>>(edges.begin(), edges.end());
        size_t n_edges = vec_edges.size();

        auto parquet_begin = Time::now();

        // Parallel Parquet write
        omp_set_num_threads(user_threads);
#pragma omp parallel
        {
            int tid = omp_get_thread_num();
            int nthreads = omp_get_num_threads();
            // n_edges is size_t, so tid * n_edges promotes to 64-bit before the
            // multiply (no int overflow at >2^31 edges); bounds/index stay size_t.
            size_t edge_start = (size_t)tid * n_edges / nthreads;
            size_t edge_end = (size_t)(tid + 1) * n_edges / nthreads;

            for (size_t ei = edge_start; ei < edge_end; ei++) {
                auto& edge = vec_edges[ei];
                uint64_t shared_features = edge.second;
                uint32_t source_1 = edge.first.first;
                uint32_t source_2 = edge.first.second;
                uint32_t source_1_features = groupID_to_featureCount.at(source_1);
                uint32_t source_2_features = groupID_to_featureCount.at(source_2);
                uint32_t minimum_source_features = min(source_1_features, source_2_features);

                // Similarity metrics with division-by-zero protection
                double containment_val = (minimum_source_features > 0)
                    ? ((double)shared_features / minimum_source_features) * 100 : 0.0;
                double ochiai_val = (source_1_features > 0 && source_2_features > 0)
                    ? 100 * ((double)shared_features / sqrt((double)source_1_features * (double)source_2_features)) : 0.0;
                double jaccard_val = (source_1_features + source_2_features - shared_features > 0)
                    ? 100 * ((double)shared_features / (source_1_features + source_2_features - shared_features)) : 0.0;
                double csi_val = (source_1_features > 0 && source_2_features > 0)
                    ? 100 * ((double)shared_features * shared_features / ((double)source_1_features * source_2_features)) : 0.0;
                double dice_val = (source_1_features + source_2_features > 0)
                    ? 100 * (2.0 * shared_features / (source_1_features + source_2_features)) : 0.0;
                double odds_ratio_val = odds_ratio(shared_features, source_1_features, source_2_features, population_size);

                // Apply cutoff
                double cutoff_val = 0;
                if (cutoff_distance_type == "containment") cutoff_val = containment_val;
                else if (cutoff_distance_type == "ochiai") cutoff_val = ochiai_val;
                else if (cutoff_distance_type == "jaccard") cutoff_val = jaccard_val;
                if (cutoff_val < cutoff_threshold) continue;

                ParquetPairRecord prec;
                prec.group_1_id = source_1;
                prec.group_2_id = source_2;
                prec.shared_features = shared_features;
                prec.containment = static_cast<float>(containment_val);
                prec.ochiai = static_cast<float>(ochiai_val);
                prec.jaccard = static_cast<float>(jaccard_val);
                prec.csi = static_cast<float>(csi_val);
                prec.dice = static_cast<float>(dice_val);
                prec.odds_ratio = static_cast<float>(odds_ratio_val);
                prec.has_pvalue = calculate_pvalue;
                if (calculate_pvalue) {
                    prec.pvalue = fastHyperPValue(shared_features, source_1_features, source_2_features, population_size);
                }

                parquet_writer.write_record(tid, prec);
            }
        }

        parquet_writer.finalize(full_command, population_size, cutoff_distance_type, cutoff_threshold);
        cout << "[dev] parallel parquet write: " << std::chrono::duration<double, std::milli>(Time::now() - parquet_begin).count() / 1000 << " secs" << endl;
        cout << "[dev] wrote parquet pairwise to: " << parquet_dir << "/ (" << parquet_writer.get_stats().num_pairs << " pairs)" << endl;

        // --- Sequential TSV + .dbrp output (legacy, kept for compatibility) ---
        auto legacy_begin = Time::now();

        for (const auto& edge : edges) {

            flat_hash_map<string, double> distance_metrics;

            uint64_t shared_features = edge.second;
            uint32_t source_1 = edge.first.first;
            uint32_t source_2 = edge.first.second;
            uint32_t source_1_features = groupID_to_featureCount.at(source_1);
            uint32_t source_2_features = groupID_to_featureCount.at(source_2);
            uint32_t minimum_source_features = min(source_1_features, source_2_features);
            uint32_t maximum_source_features = max(source_1_features, source_2_features);

            // Warn once if any pair has >1000x size difference
            static bool size_ratio_warned = false;
            if (!size_ratio_warned && minimum_source_features > 0) {
                double size_ratio = (double)maximum_source_features / minimum_source_features;
                if (size_ratio > 1000) {
                    cerr << "NOTE: Some pairs have >1000x size difference. "
                         << "Containment may be more interpretable than Jaccard/Ochiai for such pairs." << endl;
                    size_ratio_warned = true;
                }
            }

            // Similarity metrics with division-by-zero protection
            distance_metrics["containment"] = (minimum_source_features > 0)
                ? ((double)shared_features / minimum_source_features) * 100 : 0.0;
            distance_metrics["ochiai"] = (source_1_features > 0 && source_2_features > 0)
                ? 100 * ((double)shared_features / sqrt((double)source_1_features * (double)source_2_features)) : 0.0;
            distance_metrics["jaccard"] = (source_1_features + source_2_features - shared_features > 0)
                ? 100 * ((double)shared_features / (source_1_features + source_2_features - shared_features)) : 0.0;
            distance_metrics["csi"] = (source_1_features > 0 && source_2_features > 0)
                ? 100 * ((double)shared_features * shared_features / ((double)source_1_features * source_2_features)) : 0.0;
            distance_metrics["dice"] = (source_1_features + source_2_features > 0)
                ? 100 * (2.0 * shared_features / (source_1_features + source_2_features)) : 0.0;

            int k_shared_features = shared_features;
            int s_source_1_features = source_1_features;
            int M_source_2_features = source_2_features;
            int N_population_size = population_size;

            distance_metrics["odds_ratio"] = odds_ratio(k_shared_features, s_source_1_features, M_source_2_features, N_population_size);

            if (distance_metrics[cutoff_distance_type] < cutoff_threshold) continue;

            distances_stats.add_stat("containment", distance_metrics["containment"]);
            distances_stats.add_stat("ochiai", distance_metrics["ochiai"]);
            distances_stats.add_stat("jaccard", distance_metrics["jaccard"]);
            distances_stats.add_stat("csi", distance_metrics["csi"]);
            distances_stats.add_stat("dice", distance_metrics["dice"]);
            distances_stats.update_odds_ratio_stats(distance_metrics["odds_ratio"]);

            // Update .dbrp stats
            auto map_to_bucket = [](double value) -> int {
                int lower = static_cast<int>(std::floor(value / 5)) * 5;
                int idx = lower / 5;
                if (idx >= 20) idx = 20;
                return idx;
            };
            dbrp_stats.histograms[0].bucket_counts[map_to_bucket(distance_metrics["containment"])]++;
            dbrp_stats.histograms[1].bucket_counts[map_to_bucket(distance_metrics["ochiai"])]++;
            dbrp_stats.histograms[2].bucket_counts[map_to_bucket(distance_metrics["jaccard"])]++;
            dbrp_stats.histograms[3].bucket_counts[map_to_bucket(distance_metrics["csi"])]++;
            dbrp_stats.histograms[4].bucket_counts[map_to_bucket(distance_metrics["dice"])]++;

            double or_val = distance_metrics["odds_ratio"];
            if (or_val < dbrp_stats.min_odds_ratio) dbrp_stats.min_odds_ratio = or_val;
            if (or_val > dbrp_stats.max_odds_ratio) dbrp_stats.max_odds_ratio = or_val;

            // Write .dbrp record
            PairRecord rec;
            rec.group_1_id = source_1;
            rec.group_2_id = source_2;
            rec.shared_features = shared_features;
            rec.containment = static_cast<float>(distance_metrics["containment"]);
            rec.ochiai = static_cast<float>(distance_metrics["ochiai"]);
            rec.jaccard = static_cast<float>(distance_metrics["jaccard"]);
            rec.csi = static_cast<float>(distance_metrics["csi"]);
            rec.dice = static_cast<float>(distance_metrics["dice"]);
            rec.odds_ratio = static_cast<float>(distance_metrics["odds_ratio"]);
            // p-value computed once and reused for both the .dbrp record and the TSV
            // column below (it was evaluated twice per edge — the boost hypergeometric
            // CDF is the most expensive scalar op in this loop).
            double pvalue_val = calculate_pvalue
                ? fastHyperPValue(shared_features, source_1_features, source_2_features, population_size)
                : 0.0;
            if (calculate_pvalue) {
                rec.pvalue = pvalue_val;
            }
            pw.write_record(rec);

            // Write TSV line
            myfile << source_1
                << '\t' << source_2
                << '\t' << namesMap[source_1]
                << '\t' << namesMap[source_2]
                << '\t' << shared_features
                << '\t' << formatDouble(distance_metrics["containment"])
                << '\t' << formatDouble(distance_metrics["ochiai"])
                << '\t' << formatDouble(distance_metrics["jaccard"])
                << '\t' << formatDouble(distance_metrics["csi"])
                << '\t' << formatDouble(distance_metrics["dice"])
                << '\t' << formatDouble(distance_metrics["odds_ratio"]);

            if (calculate_pvalue) {
                myfile << '\t' << pvalue_val;
            }

            myfile << '\n';
        }

        myfile.close();
        distances_stats.stats_to_json_file(index_prefix + "_DBRetina_pairwise_stats.json");

        // Finalize .dbrp
        pw.finalize_write(dbrp_stats, dbrp_metadata);
        cout << "[dev] legacy TSV + .dbrp write: " << std::chrono::duration<double, std::milli>(Time::now() - legacy_begin).count() / 1000 << " secs" << endl;
        cout << "[dev] wrote .dbrp binary pairwise file: " << index_prefix << "_DBRetina_pairwise.dbrp" << endl;
    }
}