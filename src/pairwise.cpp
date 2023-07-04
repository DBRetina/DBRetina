#include <iostream>
#include <cstdint>
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

using BINS_KMER_COUNT = phmap::parallel_flat_hash_map<
    std::string, uint32_t,
    phmap::priv::hash_default_hash<std::string>,
    phmap::priv::hash_default_eq<std::string>,
    std::allocator<std::pair<const std::string, uint32_t>>,
    1,
    std::mutex>;

typedef std::chrono::high_resolution_clock Time;

class Combo {

public:
    Combo() = default;

    std::vector<std::pair<uint32_t, uint32_t>> combs;

    void combinations(int n) {
        this->combs.clear();
        this->combs.reserve((n * (n - 1)) / 2);
        this->comb(n, this->r, this->arr);
    }

private:
    int* arr = new int[2];
    int r = 2;

    void comb(int n, int r, int* arr) {
        for (int i = n; i >= r; i--) {
            // choose the first element
            arr[r - 1] = i;
            if (r > 1) { // if still needs to choose
                // recursive into smaller problem
                comb(i - 1, r - 1, arr);

            }
            else {
                this->combs.emplace_back(arr[0] - 1, arr[1] - 1);
            }
        }
    }

};

class Stats {
private:
    flat_hash_map<string, flat_hash_map<string, uint64_t>> stats;

public:
    double min_odds_ratio = 1.0;
    double max_odds_ratio = 1.0;
    Stats() {
        vector<string> distances = { "containment", "ochiai", "jaccard" };
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

inline uint64_t get_population_size(string filename) {
    return 44260; // TODO - remove this

    std::ifstream inputFile(filename);

    if (!inputFile.is_open()) {
        std::cerr << "Error opening the file: " << filename << std::endl;
        return 0;
    }

    std::string line;
    uint64_t counter = 0;
    while (std::getline(inputFile, line)) {
        if (line.find("features:") != std::string::npos) {
            std::string number = line.substr(line.find(":") + 1);
            inputFile.close();
            return stoi(number);
        }
    }
    inputFile.close();
    return 0;

}


namespace kSpider {

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

    // TODO - check how to integrate it in our code
    double calcExpectedSuccesses(int s, int M, int N) {
        // this function calculates the expected successes in the source
        // which means the expected number of genes that are in the source and in the target
        return (double)s * M / N;
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
        boost::math::hypergeometric_distribution<> hg(M, s, N);
        double pvalue;

        if (isOverEnrichment) {
            pvalue = 1 - boost::math::cdf(hg, k - 1);
        }
        else {
            pvalue = boost::math::cdf(hg, k);
        }

        return pvalue;
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
        */
        boost::math::hypergeometric_distribution<> hg(M, s, N);
        return 1 - boost::math::cdf(hg, k - 1);
    }

    // Disabled for now
    /*
    double fisher_exact(int k_shared_kmers, int source_1_kmers, int source_2_kmers, int population_size) {
        int a = k_shared_kmers;
        int b = source_1_kmers - k_shared_kmers;
        int c = source_2_kmers - k_shared_kmers;
        int d = population_size - source_1_kmers - c;

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
            this function calculates the odds ratio of a gene being in the source
            which is the ratio of the odds of the gene being in the source to the odds of the gene being in the population.
            If we have two gene sets with odds ratio equals to X, that means that
            the odds of a gene being in the first set is X times the odds of the gene being in the second set.
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

        // Read colors
        int_vec_map color_to_ids; // = new phmap::flat_hash_map<uint64_t, phmap::flat_hash_set<uint32_t>>;
        string colors_map_file = index_prefix + "_color_to_sources.bin";
        string extra_file = index_prefix + ".extra";
        uint64_t population_size = get_population_size(extra_file);
        cout << "population_size: " << population_size << endl;
        load_colors_to_sources(colors_map_file, &color_to_ids);

        // [DEBUG CODE] dump color to ids in a file
        // ofstream myfile222;
        // myfile222.open(index_prefix + "_color_to_sources.txt");
        // for (auto& color_to_ids_pair : color_to_ids) {
        //     uint64_t color = color_to_ids_pair.first;
        //     vector<uint32_t> ids = color_to_ids_pair.second;
        //     myfile222 << color << ":" << ids.size() << endl;
        // }
        // myfile222.close();

        flat_hash_map<int, std::string> namesMap;
        load_namesMap(index_prefix + ".namesMap", namesMap);
        assert(namesMap.size());

        auto begin_time = Time::now();

        cout << "mapping colors to groups: " << std::chrono::duration<double, std::milli>(Time::now() - begin_time).count() / 1000 << " secs" << endl;

        begin_time = Time::now();
        int_int_map colorsCount;
        load_colors_count(index_prefix + "_color_count.bin", colorsCount);

        // DEBUG
        // dump colors count to file 
        /*
        ofstream myfile3;
        myfile3.open(index_prefix + "_color_count.txt");
        for (auto& color_count_pair : colorsCount) {
            uint64_t color = color_count_pair.first;
            uint64_t count = color_count_pair.second;
            myfile3 << color << ":" << count << endl;
        }
        myfile3.close();


        // TODO: should be csv, rename later.
        std::ifstream data(index_prefix + "_DBRetina_colorCount.tsv");
        if (!data.is_open()) std::exit(EXIT_FAILURE);
        std::string str;
        std::getline(data, str); // skip the first line
        while (std::getline(data, str))
        {
            std::istringstream iss(str);
            std::string token;
            vector<uint32_t> tmp;
            while (std::getline(iss, token, ','))
                tmp.push_back(stoi(token));
            colorsCount.insert(make_pair(tmp[0], tmp[1]));
        }
        */

        cout << "parsing index colors: " << std::chrono::duration<double, std::milli>(Time::now() - begin_time).count() / 1000 << " secs" << endl;
        begin_time = Time::now();

        // for (const auto& record : color_to_ids) {
        //     uint32_t colorCount = colorsCount[record.first];
        //     for (auto group_id : record.second) {
        //         groupID_to_kmerCount[group_id] += colorCount;
        //     }
        // }

        // Loading kmer counts
        flat_hash_map<uint32_t, uint32_t> groupID_to_kmerCount;
        string _file_id_to_kmer_count = index_prefix + "_groupID_to_featureCount.bin";
        phmap::BinaryInputArchive ar_in_kmer_count(_file_id_to_kmer_count.c_str());
        groupID_to_kmerCount.phmap_load(ar_in_kmer_count);
        assert(groupID_to_kmerCount.size());


        std::ofstream fstream_kmerCount;
        fstream_kmerCount.open(index_prefix + "_DBRetina_featuresNo.tsv");
        fstream_kmerCount << "ID\tgroup\tfeatures\n";
        uint64_t counter = 0;
        for (const auto& item : groupID_to_kmerCount) {
            fstream_kmerCount << ++counter << '\t' << item.first << '\t' << item.second << '\n';
        }
        fstream_kmerCount.close();
        cout << "features counting: " << std::chrono::duration<double, std::milli>(Time::now() - begin_time).count() / 1000 << " secs" << endl;

        // Loading done

        begin_time = Time::now();
        clock_t begin_detailed_pairwise_comb, begin_detailed_pairwise_edges, begin_detailed_pairwise_edges_insertion;
        double detailed_pairwise_comb = 0.0;
        double detailed_pairwise_edges = 0.0;
        double detailed_pairwise_edges_insertion = 0.0;

        PAIRS_COUNTER edges;

        // convert map to vec for parallelization purposes.
        auto vec_color_to_ids = std::vector<std::pair<uint32_t, vector<uint32_t>>>(color_to_ids.begin(), color_to_ids.end());

        cerr << "number of colors = " << vec_color_to_ids.size() << endl;

        double average_color_size = 0.0;
        for (auto const& item : vec_color_to_ids) {
            average_color_size += item.second.size();
        }
        average_color_size /= vec_color_to_ids.size();

        cerr << "average color size = " << (int)average_color_size << endl;

        int thread_num, num_threads, start, end, vec_i;
        int n = vec_color_to_ids.size();

        omp_set_num_threads(user_threads);
        begin_time = Time::now();

#pragma omp parallel private(vec_i,thread_num,num_threads,start,end)
        {
            thread_num = omp_get_thread_num();
            num_threads = omp_get_num_threads();
            start = thread_num * n / num_threads;
            end = (thread_num + 1) * n / num_threads;

            for (vec_i = start; vec_i != end; ++vec_i) {
                auto item = vec_color_to_ids[vec_i];
                Combo combo = Combo();
                combo.combinations(item.second.size());
                for (uint32_t i = 0; i < combo.combs.size(); i++) {
                    // for (auto const& seq_pair : combo.combs) {
                    auto const& seq_pair = combo.combs[i];
                    uint32_t _seq1 = item.second[seq_pair.first];
                    uint32_t _seq2 = item.second[seq_pair.second];
                    ascending(_seq1, _seq2);

                    auto _p = make_pair(_seq1, _seq2);
                    uint32_t ccount = colorsCount[item.first];
                    edges.try_emplace_l(_p,
                        [ccount](PAIRS_COUNTER::value_type& v) { v.second += ccount; }, // called only when key was already present
                        ccount
                    );
                }
            }
        }

        cout << "pairwise hashmap construction: " << std::chrono::duration<double, std::milli>(Time::now() - begin_time).count() / 1000 << " secs" << endl;
        cout << "Number of pairwise comparisons: " << edges.size() << endl;
        cout << "writing pairwise matrix to " << index_prefix << "_DBRetina_pairwise.tsv" << endl;

        Stats distances_stats;

        auto formatDouble = [](double val) {
            char buffer[20];
            std::snprintf(buffer, sizeof(buffer), "%.1f", val);
            return std::string(buffer);
            };

        std::ofstream myfile;
        myfile.open(index_prefix + "_DBRetina_pairwise.tsv");
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
            << "\todds_ratio";
        if (calculate_pvalue) { myfile << "\tpvalue"; }
        // << "\texpected_successes"
        // << "\tfold_change"

        myfile << '\n';
        uint64_t line_count = 0;


        for (const auto& edge : edges) {

            flat_hash_map<string, double> distance_metrics;

            uint64_t shared_kmers = edge.second;
            uint32_t source_1 = edge.first.first;
            uint32_t source_2 = edge.first.second;
            uint32_t source_1_kmers = groupID_to_kmerCount[source_1];
            uint32_t source_2_kmers = groupID_to_kmerCount[source_2];
            uint32_t minimum_source_kmers = min(source_1_kmers, source_2_kmers);

            // containment
            distance_metrics["containment"] = ((double)shared_kmers / minimum_source_kmers) * 100;


            // Ochiai distance
            distance_metrics["ochiai"] = 100 * ((double)shared_kmers / sqrt((double)source_1_kmers * (double)source_2_kmers));


            // Jaccard distance (if size of samples is roughly similar)
            // J(A, B) = 1 - |A ∩ B| / (|A| + |B| - |A ∩ B|) <- this is distance not similarity
            distance_metrics["jaccard"] = 100 * ((double)shared_kmers / (source_1_kmers + source_2_kmers - shared_kmers));

            // Kulczynski distance needs abundance of each sample
            // double kulczynski = (double)shared_kmers / (source_1_kmers + source_2_kmers) * 2;

            // p-value using hypergeometric CDF
            int k_shared_kmers = shared_kmers;  // Number of successes
            int s_source_1_kmers = source_1_kmers;  // Sample size
            int M_source_2_kmers = source_2_kmers;  // Number of successes in the population
            int N_population_size = population_size;  // Population size


            // Odds ratio
            distance_metrics["odds_ratio"] = odds_ratio(k_shared_kmers, s_source_1_kmers, M_source_2_kmers, N_population_size);


            // auto [pvalue, expectedSuccesses, fold_change] = enrichmentAnalysis(shared_kmers, source_1_kmers, source_2_kmers, population_size, true);
            // distance_metrics["expected_successes"] = expectedSuccesses;
            // distance_metrics["fold_change"] = fold_change;

            if (distance_metrics[cutoff_distance_type] < cutoff_threshold) continue;


            distances_stats.add_stat("containment", distance_metrics["containment"]);
            distances_stats.add_stat("ochiai", distance_metrics["ochiai"]);
            distances_stats.add_stat("jaccard", distance_metrics["jaccard"]);
            distances_stats.update_odds_ratio_stats(distance_metrics["odds_ratio"]);



            myfile << source_1
                << '\t' << source_2
                << '\t' << namesMap[source_1]
                << '\t' << namesMap[source_2]
                << '\t' << shared_kmers
                << '\t' << formatDouble(distance_metrics["containment"])
                << '\t' << formatDouble(distance_metrics["ochiai"])
                << '\t' << formatDouble(distance_metrics["jaccard"])
                << '\t' << formatDouble(distance_metrics["odds_ratio"]);

            if (calculate_pvalue) {
                myfile << '\t' << fastHyperPValue(shared_kmers, source_1_kmers, source_2_kmers, population_size);
            }

            myfile << '\n';
        }


        myfile.close();
        distances_stats.stats_to_json_file(index_prefix + "_DBRetina_pairwise_stats.json");
    }
}