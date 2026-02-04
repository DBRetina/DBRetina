#pragma once

#include <string>
#include <vector>
#include <cstdint>
#include <mutex>
#include <atomic>
#include <arrow/api.h>
#include <arrow/io/api.h>
#include <parquet/arrow/writer.h>
#include <parquet/properties.h>
#include "parallel_hashmap/phmap.h"

struct ParquetPairRecord {
    uint32_t group_1_id;
    uint32_t group_2_id;
    uint64_t shared_features;
    float containment;
    float ochiai;
    float jaccard;
    float csi;
    float dice;
    float odds_ratio;
    double pvalue;
    bool has_pvalue;
};

struct ParquetWriterStats {
    // Per-metric histogram: 21 buckets (0-5, 5-10, ..., 95-100, 100-100)
    std::vector<std::vector<uint64_t>> histograms;  // [5 metrics][21 buckets]
    double min_odds_ratio;
    double max_odds_ratio;
    uint64_t num_pairs;

    ParquetWriterStats();
    void update(const ParquetPairRecord& rec);
    void merge(const ParquetWriterStats& other);
};

class ParquetPairwiseWriter {
public:
    ParquetPairwiseWriter(const std::string& output_dir,
                          bool has_pvalue,
                          int num_threads,
                          const phmap::flat_hash_map<int, std::string>& names_map);

    ~ParquetPairwiseWriter();

    // Called from each OpenMP thread — thread-safe via thread_id routing
    void write_record(int thread_id, const ParquetPairRecord& rec);

    // Called after all records are written — flushes buffers and writes sidecar files
    void finalize(const std::string& full_command,
                  uint64_t population_size,
                  const std::string& cutoff_metric,
                  double cutoff_threshold);

    // Get merged statistics after finalize
    const ParquetWriterStats& get_stats() const { return merged_stats_; }

private:
    struct ThreadWriter {
        std::shared_ptr<arrow::io::FileOutputStream> file;
        std::unique_ptr<parquet::arrow::FileWriter> writer;

        // Builders for current batch
        std::shared_ptr<arrow::UInt32Builder> group_1_id_builder;
        std::shared_ptr<arrow::UInt32Builder> group_2_id_builder;
        std::shared_ptr<arrow::UInt64Builder> shared_features_builder;
        std::shared_ptr<arrow::FloatBuilder> containment_builder;
        std::shared_ptr<arrow::FloatBuilder> ochiai_builder;
        std::shared_ptr<arrow::FloatBuilder> jaccard_builder;
        std::shared_ptr<arrow::FloatBuilder> csi_builder;
        std::shared_ptr<arrow::FloatBuilder> dice_builder;
        std::shared_ptr<arrow::FloatBuilder> odds_ratio_builder;
        std::shared_ptr<arrow::DoubleBuilder> pvalue_builder;

        ParquetWriterStats stats;
        uint64_t batch_count;

        // Track which group IDs appear in this partition
        phmap::flat_hash_set<uint32_t> group_ids_seen;
    };

    void init_thread_writer(int thread_id);
    void flush_batch(ThreadWriter& tw);
    std::shared_ptr<arrow::Schema> make_schema() const;
    void write_manifest(const std::string& full_command,
                        uint64_t population_size,
                        const std::string& cutoff_metric,
                        double cutoff_threshold);
    void write_names_parquet();
    void write_statistics_json();
    void write_group_index_parquet();

    std::string output_dir_;
    bool has_pvalue_;
    int num_threads_;
    phmap::flat_hash_map<int, std::string> names_map_;
    std::vector<std::unique_ptr<ThreadWriter>> thread_writers_;
    ParquetWriterStats merged_stats_;

    static constexpr size_t BATCH_SIZE = 65536;  // Records per batch before flush
};
