#pragma once

#include <string>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <functional>
#include <cmath>
#include <stdexcept>
#include <cstring>
#include "parallel_hashmap/phmap.h"

// Magic bytes for .dbrp format
static constexpr char DBRP_MAGIC[4] = {'D', 'B', 'R', 'P'};
static constexpr uint32_t DBRP_FORMAT_VERSION = 1;

// Section IDs for .dbrp
enum class DBRPSection : uint32_t {
    NAMES_TABLE   = 1,  // group_id -> group_name
    PAIR_RECORDS  = 2,  // fixed-width pair records (bulk data)
    GROUP_INDEX   = 3,  // per-group record ordinal lists
    STATISTICS    = 4,  // pre-computed histogram buckets
    METADATA      = 5,  // JSON metadata
};

// Metric bit positions in metric_flags
enum class DBRPMetric : uint8_t {
    CONTAINMENT = 0,
    OCHIAI      = 1,
    JACCARD     = 2,
    CSI         = 3,
    DICE        = 4,
    ODDS_RATIO  = 5,
    PVALUE      = 6,
};

struct DBRPSectionEntry {
    uint32_t section_id;
    uint64_t offset;
    uint64_t length;
    uint32_t crc32;
};

struct PairRecord {
    uint32_t group_1_id;
    uint32_t group_2_id;
    uint64_t shared_features;
    float    containment;
    float    ochiai;
    float    jaccard;
    float    csi;
    float    dice;
    float    odds_ratio;
    double   pvalue;   // NaN if not stored

    PairRecord() : group_1_id(0), group_2_id(0), shared_features(0),
                   containment(0), ochiai(0), jaccard(0), csi(0), dice(0),
                   odds_ratio(0), pvalue(std::nan("")) {}
};

// Statistics histogram for a single metric
struct MetricHistogram {
    uint8_t metric_id;
    std::vector<uint64_t> bucket_counts;  // 21 buckets: 0-5, 5-10, ..., 95-100, 100-100
};

struct PairwiseStatistics {
    std::vector<MetricHistogram> histograms;
    double min_odds_ratio;
    double max_odds_ratio;
};

class DBRetinaPairwise {
public:
    DBRetinaPairwise() = default;
    ~DBRetinaPairwise() = default;

    // Non-copyable, movable
    DBRetinaPairwise(const DBRetinaPairwise&) = delete;
    DBRetinaPairwise& operator=(const DBRetinaPairwise&) = delete;
    DBRetinaPairwise(DBRetinaPairwise&&) = default;
    DBRetinaPairwise& operator=(DBRetinaPairwise&&) = default;

    // ----- Writing -----

    // Start writing a new .dbrp file
    // names_map: group_id -> group_name
    void begin_write(
        const std::string& filepath,
        uint8_t metric_flags,
        const phmap::flat_hash_map<int, std::string>& names_map,
        const std::string& metadata_json
    );

    // Write a single pair record
    void write_record(const PairRecord& rec);

    // Write a batch of records
    void write_batch(const PairRecord* records, size_t count);

    // Finalize: writes GROUP_INDEX, STATISTICS, METADATA, TOC; patches header
    void finalize_write(const PairwiseStatistics& stats, const std::string& metadata_json);

    // ----- Reading -----

    static DBRetinaPairwise open(const std::string& filepath);

    // Metadata accessors
    uint64_t get_num_pairs() const { return num_pairs_; }
    uint32_t get_num_groups() const { return num_groups_; }
    uint8_t  get_metric_flags() const { return metric_flags_; }
    bool     has_metric(DBRPMetric m) const { return (metric_flags_ >> static_cast<uint8_t>(m)) & 1; }

    // Name resolution
    std::string get_group_name(uint32_t group_id) const;
    const std::map<uint32_t, std::string>& get_names_map() const { return names_map_; }

    // Load metadata JSON
    std::string get_metadata_json() const;

    // Load statistics
    PairwiseStatistics get_statistics() const;
    std::string get_statistics_json() const;

    // Iterate all pairs with optional metric filter
    // Returns records where metric >= cutoff
    // metric_id: 0=containment, 1=ochiai, ..., 6=pvalue
    // cutoff: minimum value to include (pass 0.0 and metric_id=0 for no filter)
    std::vector<PairRecord> filter_pairs(uint8_t metric_id, double cutoff) const;

    // Get all pairs involving a specific group
    std::vector<PairRecord> query_group(uint32_t group_id, uint8_t metric_id, double cutoff) const;

    // Get all pairs (no filter)
    std::vector<PairRecord> iterate_all() const;

    // Random access by record ordinal (0-based)
    PairRecord read_record(uint64_t record_index) const;

    const std::string& get_filepath() const { return filepath_; }

private:
    std::string filepath_;
    uint64_t num_pairs_ = 0;
    uint32_t num_groups_ = 0;
    uint8_t  metric_flags_ = 0;
    uint8_t  record_size_ = 0;

    // TOC
    std::vector<DBRPSectionEntry> toc_;

    // Loaded on open()
    std::map<uint32_t, std::string> names_map_;

    // Group index (loaded lazily)
    mutable bool group_index_loaded_ = false;
    mutable std::map<uint32_t, std::vector<uint64_t>> group_index_;
    void load_group_index() const;

    // For writing
    std::ofstream write_stream_;
    std::vector<DBRPSectionEntry> write_toc_;
    bool writing_ = false;
    uint64_t write_num_pairs_ = 0;
    uint64_t pair_records_section_offset_ = 0;

    // Group index accumulation during writing
    std::map<uint32_t, std::vector<uint64_t>> write_group_index_;

    // Helpers
    const DBRPSectionEntry* find_section(DBRPSection section_id) const;
    std::vector<char> read_section_raw(DBRPSection section_id) const;
    uint8_t compute_record_size(uint8_t flags) const;
    void serialize_record(const PairRecord& rec, char* buffer) const;
    PairRecord deserialize_record(const char* buffer) const;

    // Write a section (data + TOC entry)
    void write_section(DBRPSection section_id, const char* data, size_t length);
};
