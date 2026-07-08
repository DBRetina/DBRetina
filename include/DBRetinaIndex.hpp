#pragma once

#include <string>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <vector>
#include <queue>
#include <stdexcept>
#include <cstring>
#include "parallel_hashmap/phmap.h"
#include "parallel_hashmap/phmap_dump.h"
#include "DBRetina_kf.hpp"

// Magic bytes for .dbri format
static constexpr char DBRI_MAGIC[4] = {'D', 'B', 'R', 'I'};
static constexpr uint32_t DBRI_FORMAT_VERSION = 1;

// Section IDs
enum class DBRISection : uint32_t {
    PHMAP              = 1,   // feature_hash -> color_id
    COLOR_TO_SOURCES   = 2,   // color_id -> vector<group_id>
    COLOR_COUNT        = 3,   // color_id -> usage_count
    GROUP_FEATURE_COUNT = 4,  // group_id -> feature_count
    NAMES_MAP          = 5,   // group_id -> group_name
    METADATA           = 6,   // JSON metadata string
    RAW_GENE_SETS      = 7,   // original gene set -> genes mapping (optional)
    HASHED_GENE_SETS   = 8,   // gene set -> hashed gene IDs
    TAGS_MAP           = 9,   // color_string -> color_id (for incremental)
    FREE_COLORS        = 10,  // available color IDs for reuse
    GROUP_TO_FEATURE_SET  = 11,  // group_name -> flat_hash_set<uint64_t> (for merge)
};

struct SectionEntry {
    uint32_t section_id;
    uint64_t offset;
    uint64_t length;
    uint32_t crc32;
};

// CRC32 computation (standard polynomial)
uint32_t compute_crc32(const char* data, size_t length);

class DBRetinaIndex {
public:
    // Construction
    DBRetinaIndex() = default;
    ~DBRetinaIndex() = default;

    // Non-copyable (ofstream member), but movable
    DBRetinaIndex(const DBRetinaIndex&) = delete;
    DBRetinaIndex& operator=(const DBRetinaIndex&) = delete;
    DBRetinaIndex(DBRetinaIndex&&) = default;
    DBRetinaIndex& operator=(DBRetinaIndex&&) = default;

    // ----- Writing -----
    // Call begin_write, then write_section for each section, then finalize_write
    void begin_write(const std::string& filepath);

    // Write a raw data section
    void write_section(DBRISection section_id, const char* data, size_t length);

    // Convenience writers for each data structure
    void write_phmap(DBRetina_PHMAP* frame);
    void write_color_to_sources(
        const phmap::flat_hash_map<uint64_t, std::vector<uint32_t>>& legend,
        const phmap::flat_hash_map<uint64_t, uint64_t>& colorsCount);
    void write_color_count(const phmap::flat_hash_map<uint64_t, uint64_t>& colorsCount);
    void write_group_feature_count(const phmap::flat_hash_map<uint32_t, uint32_t>& groupID_to_featureCount);
    void write_names_map(
        const phmap::flat_hash_map<std::string, std::string>& namesMap,
        const phmap::flat_hash_map<std::string, uint64_t>& groupNameMap);
    void write_metadata(const std::string& metadata_json);
    void write_raw_gene_sets(const std::string& json_content);
    void write_hashed_gene_sets(const std::string& json_content);
    void write_tags_map(const phmap::flat_hash_map<std::string, uint64_t>& tagsMap);
    void write_free_colors(const std::priority_queue<uint64_t, std::vector<uint64_t>, std::greater<uint64_t>>& freeColors);
    void write_group_to_feature_set(const phmap::parallel_flat_hash_map<std::string, phmap::flat_hash_set<uint64_t>>& groupName_to_featureSet);

    void finalize_write();

    // ----- Reading -----
    static DBRetinaIndex open(const std::string& filepath);

    // Check if a section exists
    bool has_section(DBRISection section_id) const;

    // Read raw section data
    std::vector<char> read_section_raw(DBRISection section_id) const;

    // Convenience readers for each data structure
    DBRetina_PHMAP* load_phmap() const;

    void load_color_to_sources(
        phmap::parallel_flat_hash_map<uint32_t, std::vector<uint32_t>,
            std::hash<uint32_t>, std::equal_to<uint32_t>,
            std::allocator<std::pair<const uint32_t, std::vector<uint32_t>>>, 1>& color_to_ids) const;

    void load_color_count(
        phmap::parallel_flat_hash_map<uint32_t, uint32_t,
            std::hash<uint32_t>, std::equal_to<uint32_t>,
            std::allocator<std::pair<const uint32_t, uint32_t>>, 1>& colorsCount) const;

    void load_group_feature_count(phmap::flat_hash_map<uint32_t, uint32_t>& groupID_to_featureCount) const;

    void load_names_map(phmap::flat_hash_map<int, std::string>& namesMap) const;

    std::string load_metadata() const;

    std::string load_raw_gene_sets() const;
    std::string load_hashed_gene_sets() const;

    void load_tags_map(phmap::flat_hash_map<std::string, uint64_t>& tagsMap) const;

    void load_free_colors(std::priority_queue<uint64_t, std::vector<uint64_t>, std::greater<uint64_t>>& freeColors) const;

    void load_group_to_feature_set(phmap::parallel_flat_hash_map<std::string, phmap::flat_hash_set<uint64_t>>& groupName_to_featureSet) const;

    // Metadata helpers
    uint64_t get_population_size() const;
    uint32_t get_num_groups() const;
    uint64_t get_max_group_id() const;

    const std::string& get_filepath() const { return filepath_; }

private:
    std::string filepath_;
    std::vector<SectionEntry> toc_;

    // For writing
    std::ofstream write_stream_;
    std::vector<SectionEntry> write_toc_;
    bool writing_ = false;

    // Helper to find section in TOC
    const SectionEntry* find_section(DBRISection section_id) const;
};
