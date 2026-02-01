#include "DBRetinaPairwise.hpp"
#include <sstream>
#include <algorithm>
#include <cstring>

// Use the CRC32 function declared in DBRetinaIndex.hpp
extern uint32_t compute_crc32(const char* data, size_t length);

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

uint8_t DBRetinaPairwise::compute_record_size(uint8_t flags) const {
    // Base: group_1_id(4) + group_2_id(4) + shared_features(8) = 16
    uint8_t size = 16;
    // Bits 0-5: float32 metrics (4 bytes each)
    for (int i = 0; i < 6; i++) {
        if ((flags >> i) & 1) size += 4;
    }
    // Bit 6: pvalue as float64 (8 bytes)
    if ((flags >> 6) & 1) size += 8;
    return size;
}

void DBRetinaPairwise::serialize_record(const PairRecord& rec, char* buffer) const {
    char* ptr = buffer;
    std::memcpy(ptr, &rec.group_1_id, 4); ptr += 4;
    std::memcpy(ptr, &rec.group_2_id, 4); ptr += 4;
    std::memcpy(ptr, &rec.shared_features, 8); ptr += 8;

    if (has_metric(DBRPMetric::CONTAINMENT)) { std::memcpy(ptr, &rec.containment, 4); ptr += 4; }
    if (has_metric(DBRPMetric::OCHIAI))      { std::memcpy(ptr, &rec.ochiai, 4); ptr += 4; }
    if (has_metric(DBRPMetric::JACCARD))     { std::memcpy(ptr, &rec.jaccard, 4); ptr += 4; }
    if (has_metric(DBRPMetric::CSI))         { std::memcpy(ptr, &rec.csi, 4); ptr += 4; }
    if (has_metric(DBRPMetric::DICE))        { std::memcpy(ptr, &rec.dice, 4); ptr += 4; }
    if (has_metric(DBRPMetric::ODDS_RATIO))  { std::memcpy(ptr, &rec.odds_ratio, 4); ptr += 4; }
    if (has_metric(DBRPMetric::PVALUE))      { std::memcpy(ptr, &rec.pvalue, 8); ptr += 8; }
}

PairRecord DBRetinaPairwise::deserialize_record(const char* buffer) const {
    PairRecord rec;
    const char* ptr = buffer;
    std::memcpy(&rec.group_1_id, ptr, 4); ptr += 4;
    std::memcpy(&rec.group_2_id, ptr, 4); ptr += 4;
    std::memcpy(&rec.shared_features, ptr, 8); ptr += 8;

    if (has_metric(DBRPMetric::CONTAINMENT)) { std::memcpy(&rec.containment, ptr, 4); ptr += 4; }
    if (has_metric(DBRPMetric::OCHIAI))      { std::memcpy(&rec.ochiai, ptr, 4); ptr += 4; }
    if (has_metric(DBRPMetric::JACCARD))     { std::memcpy(&rec.jaccard, ptr, 4); ptr += 4; }
    if (has_metric(DBRPMetric::CSI))         { std::memcpy(&rec.csi, ptr, 4); ptr += 4; }
    if (has_metric(DBRPMetric::DICE))        { std::memcpy(&rec.dice, ptr, 4); ptr += 4; }
    if (has_metric(DBRPMetric::ODDS_RATIO))  { std::memcpy(&rec.odds_ratio, ptr, 4); ptr += 4; }
    if (has_metric(DBRPMetric::PVALUE))      { std::memcpy(&rec.pvalue, ptr, 8); ptr += 8; }

    return rec;
}

// Get the float value for a given metric from a record
static float get_metric_value(const PairRecord& rec, uint8_t metric_id) {
    switch (metric_id) {
        case 0: return rec.containment;
        case 1: return rec.ochiai;
        case 2: return rec.jaccard;
        case 3: return rec.csi;
        case 4: return rec.dice;
        case 5: return rec.odds_ratio;
        case 6: return static_cast<float>(rec.pvalue);
        default: return 0.0f;
    }
}

// ---------------------------------------------------------------------------
// Writing
// ---------------------------------------------------------------------------

void DBRetinaPairwise::begin_write(
    const std::string& filepath,
    uint8_t metric_flags,
    const phmap::flat_hash_map<int, std::string>& names_map,
    const std::string& metadata_json)
{
    filepath_ = filepath;
    metric_flags_ = metric_flags;
    record_size_ = compute_record_size(metric_flags);
    num_groups_ = static_cast<uint32_t>(names_map.size());
    writing_ = true;
    write_num_pairs_ = 0;
    write_toc_.clear();
    write_group_index_.clear();

    write_stream_.open(filepath, std::ios::binary);
    if (!write_stream_.is_open()) {
        throw std::runtime_error("Cannot open file for writing: " + filepath);
    }

    // Write file header (48 bytes)
    write_stream_.write(DBRP_MAGIC, 4);
    write_stream_.write(reinterpret_cast<const char*>(&DBRP_FORMAT_VERSION), sizeof(uint32_t));
    uint64_t toc_offset_placeholder = 0;
    write_stream_.write(reinterpret_cast<const char*>(&toc_offset_placeholder), sizeof(uint64_t));
    uint64_t num_pairs_placeholder = 0;
    write_stream_.write(reinterpret_cast<const char*>(&num_pairs_placeholder), sizeof(uint64_t));
    write_stream_.write(reinterpret_cast<const char*>(&num_groups_), sizeof(uint32_t));
    write_stream_.write(reinterpret_cast<const char*>(&metric_flags_), sizeof(uint8_t));
    write_stream_.write(reinterpret_cast<const char*>(&record_size_), sizeof(uint8_t));
    // Reserved: 18 bytes to bring header to 48
    char reserved[18] = {0};
    write_stream_.write(reserved, 18);

    // Write NAMES_TABLE section
    {
        std::ostringstream oss(std::ios::binary);
        uint32_t count = static_cast<uint32_t>(names_map.size());
        oss.write(reinterpret_cast<const char*>(&count), sizeof(uint32_t));
        for (auto& [gid, name] : names_map) {
            uint32_t id = static_cast<uint32_t>(gid);
            uint16_t name_len = static_cast<uint16_t>(name.size());
            oss.write(reinterpret_cast<const char*>(&id), sizeof(uint32_t));
            oss.write(reinterpret_cast<const char*>(&name_len), sizeof(uint16_t));
            oss.write(name.data(), name_len);
        }
        std::string data = oss.str();
        write_section(DBRPSection::NAMES_TABLE, data.data(), data.size());
    }

    // Record the start offset of PAIR_RECORDS section
    pair_records_section_offset_ = static_cast<uint64_t>(write_stream_.tellp());
}

void DBRetinaPairwise::write_record(const PairRecord& rec) {
    if (!writing_) throw std::runtime_error("Not in write mode");

    std::vector<char> buffer(record_size_);
    serialize_record(rec, buffer.data());
    write_stream_.write(buffer.data(), record_size_);

    // Track in group index
    write_group_index_[rec.group_1_id].push_back(write_num_pairs_);
    write_group_index_[rec.group_2_id].push_back(write_num_pairs_);

    write_num_pairs_++;
}

void DBRetinaPairwise::write_batch(const PairRecord* records, size_t count) {
    if (!writing_) throw std::runtime_error("Not in write mode");

    std::vector<char> buffer(record_size_);
    for (size_t i = 0; i < count; i++) {
        serialize_record(records[i], buffer.data());
        write_stream_.write(buffer.data(), record_size_);

        write_group_index_[records[i].group_1_id].push_back(write_num_pairs_);
        write_group_index_[records[i].group_2_id].push_back(write_num_pairs_);

        write_num_pairs_++;
    }
}

void DBRetinaPairwise::write_section(DBRPSection section_id, const char* data, size_t length) {
    if (!writing_) throw std::runtime_error("Not in write mode");

    DBRPSectionEntry entry;
    entry.section_id = static_cast<uint32_t>(section_id);
    entry.offset = static_cast<uint64_t>(write_stream_.tellp());
    entry.length = length;
    entry.crc32 = compute_crc32(data, length);

    write_stream_.write(data, length);
    write_toc_.push_back(entry);
}

void DBRetinaPairwise::finalize_write(const PairwiseStatistics& stats, const std::string& metadata_json) {
    if (!writing_) throw std::runtime_error("Not in write mode");

    // Register the PAIR_RECORDS section in TOC
    {
        uint64_t pair_section_length = static_cast<uint64_t>(write_num_pairs_) * record_size_;

        // Compute CRC32 of the pair records section by re-reading from disk
        // For large files this could be slow, but correctness first
        DBRPSectionEntry entry;
        entry.section_id = static_cast<uint32_t>(DBRPSection::PAIR_RECORDS);
        entry.offset = pair_records_section_offset_;
        entry.length = pair_section_length;

        // For now, compute CRC by re-reading the data
        uint64_t current_pos = static_cast<uint64_t>(write_stream_.tellp());
        write_stream_.flush();

        // Read the pair records data back for CRC
        std::ifstream check_in(filepath_, std::ios::binary);
        check_in.seekg(pair_records_section_offset_);
        uint32_t crc = 0xFFFFFFFF;
        const size_t chunk_size = 1024 * 1024; // 1MB chunks
        std::vector<char> chunk(chunk_size);
        uint64_t remaining = pair_section_length;
        while (remaining > 0) {
            size_t to_read = std::min(static_cast<size_t>(remaining), chunk_size);
            check_in.read(chunk.data(), to_read);
            // Incremental CRC32
            for (size_t i = 0; i < to_read; i++) {
                // Use the same CRC table as compute_crc32
                crc = crc ^ static_cast<uint8_t>(chunk[i]);
                for (int j = 0; j < 8; j++) {
                    if (crc & 1)
                        crc = (crc >> 1) ^ 0xEDB88320;
                    else
                        crc = crc >> 1;
                }
            }
            remaining -= to_read;
        }
        check_in.close();
        entry.crc32 = crc ^ 0xFFFFFFFF;
        write_toc_.push_back(entry);
    }

    // Write GROUP_INDEX section
    {
        std::ostringstream oss(std::ios::binary);
        uint32_t num_entries = static_cast<uint32_t>(write_group_index_.size());
        oss.write(reinterpret_cast<const char*>(&num_entries), sizeof(uint32_t));

        // First pass: compute offsets within the record lists area
        // Each group's list starts right after the index table
        // Index table size = num_entries * (4 + 8 + 4) = num_entries * 16
        uint64_t lists_area_offset = 0;
        struct IndexEntry {
            uint32_t group_id;
            uint64_t list_offset;
            uint32_t list_count;
        };
        std::vector<IndexEntry> index_entries;
        index_entries.reserve(num_entries);

        for (auto& [gid, indices] : write_group_index_) {
            IndexEntry ie;
            ie.group_id = gid;
            ie.list_offset = lists_area_offset;
            ie.list_count = static_cast<uint32_t>(indices.size());
            index_entries.push_back(ie);
            lists_area_offset += indices.size() * sizeof(uint64_t);
        }

        // Write index entries
        for (auto& ie : index_entries) {
            oss.write(reinterpret_cast<const char*>(&ie.group_id), sizeof(uint32_t));
            oss.write(reinterpret_cast<const char*>(&ie.list_offset), sizeof(uint64_t));
            oss.write(reinterpret_cast<const char*>(&ie.list_count), sizeof(uint32_t));
        }

        // Write record lists
        for (auto& [gid, indices] : write_group_index_) {
            for (auto idx : indices) {
                oss.write(reinterpret_cast<const char*>(&idx), sizeof(uint64_t));
            }
        }

        std::string data = oss.str();
        write_section(DBRPSection::GROUP_INDEX, data.data(), data.size());
    }

    // Write STATISTICS section
    {
        std::ostringstream oss(std::ios::binary);
        uint32_t num_metrics = static_cast<uint32_t>(stats.histograms.size());
        oss.write(reinterpret_cast<const char*>(&num_metrics), sizeof(uint32_t));

        for (auto& hist : stats.histograms) {
            oss.write(reinterpret_cast<const char*>(&hist.metric_id), sizeof(uint8_t));
            uint32_t num_buckets = static_cast<uint32_t>(hist.bucket_counts.size());
            oss.write(reinterpret_cast<const char*>(&num_buckets), sizeof(uint32_t));
            for (auto count : hist.bucket_counts) {
                oss.write(reinterpret_cast<const char*>(&count), sizeof(uint64_t));
            }
        }

        oss.write(reinterpret_cast<const char*>(&stats.min_odds_ratio), sizeof(double));
        oss.write(reinterpret_cast<const char*>(&stats.max_odds_ratio), sizeof(double));

        std::string data = oss.str();
        write_section(DBRPSection::STATISTICS, data.data(), data.size());
    }

    // Write METADATA section
    write_section(DBRPSection::METADATA, metadata_json.data(), metadata_json.size());

    // Write TOC at current position
    uint64_t toc_offset = static_cast<uint64_t>(write_stream_.tellp());

    uint32_t num_sections = static_cast<uint32_t>(write_toc_.size());
    write_stream_.write(reinterpret_cast<const char*>(&num_sections), sizeof(uint32_t));

    for (auto& entry : write_toc_) {
        write_stream_.write(reinterpret_cast<const char*>(&entry.section_id), sizeof(uint32_t));
        write_stream_.write(reinterpret_cast<const char*>(&entry.offset), sizeof(uint64_t));
        write_stream_.write(reinterpret_cast<const char*>(&entry.length), sizeof(uint64_t));
        write_stream_.write(reinterpret_cast<const char*>(&entry.crc32), sizeof(uint32_t));
    }

    // Patch header: toc_offset at offset 8, num_pairs at offset 16
    write_stream_.seekp(8);  // after magic(4) + version(4)
    write_stream_.write(reinterpret_cast<const char*>(&toc_offset), sizeof(uint64_t));
    write_stream_.write(reinterpret_cast<const char*>(&write_num_pairs_), sizeof(uint64_t));

    write_stream_.close();
    writing_ = false;
    num_pairs_ = write_num_pairs_;
    toc_ = write_toc_;
}

// ---------------------------------------------------------------------------
// Reading
// ---------------------------------------------------------------------------

DBRetinaPairwise DBRetinaPairwise::open(const std::string& filepath) {
    DBRetinaPairwise pw;
    pw.filepath_ = filepath;

    std::ifstream in(filepath, std::ios::binary);
    if (!in.is_open()) {
        throw std::runtime_error("Cannot open .dbrp file: " + filepath);
    }

    // Read and verify magic
    char magic[4];
    in.read(magic, 4);
    if (std::memcmp(magic, DBRP_MAGIC, 4) != 0) {
        throw std::runtime_error("Invalid .dbrp file (bad magic bytes): " + filepath);
    }

    // Read version
    uint32_t version;
    in.read(reinterpret_cast<char*>(&version), sizeof(uint32_t));
    if (version != DBRP_FORMAT_VERSION) {
        throw std::runtime_error("Unsupported .dbrp format version " + std::to_string(version)
            + " (expected " + std::to_string(DBRP_FORMAT_VERSION) + ")");
    }

    // Read header fields
    in.read(reinterpret_cast<char*>(&pw.num_pairs_), sizeof(uint64_t));  // offset 8: toc_offset
    uint64_t toc_offset = pw.num_pairs_;  // This is actually toc_offset
    in.read(reinterpret_cast<char*>(&pw.num_pairs_), sizeof(uint64_t));  // offset 16: num_pairs
    in.read(reinterpret_cast<char*>(&pw.num_groups_), sizeof(uint32_t));
    in.read(reinterpret_cast<char*>(&pw.metric_flags_), sizeof(uint8_t));
    in.read(reinterpret_cast<char*>(&pw.record_size_), sizeof(uint8_t));

    // Skip reserved (18 bytes)
    in.seekg(48);

    // Seek to TOC and read it
    in.seekg(toc_offset);
    uint32_t num_sections;
    in.read(reinterpret_cast<char*>(&num_sections), sizeof(uint32_t));

    pw.toc_.resize(num_sections);
    for (uint32_t i = 0; i < num_sections; i++) {
        in.read(reinterpret_cast<char*>(&pw.toc_[i].section_id), sizeof(uint32_t));
        in.read(reinterpret_cast<char*>(&pw.toc_[i].offset), sizeof(uint64_t));
        in.read(reinterpret_cast<char*>(&pw.toc_[i].length), sizeof(uint64_t));
        in.read(reinterpret_cast<char*>(&pw.toc_[i].crc32), sizeof(uint32_t));
    }

    // Load names table eagerly
    {
        auto data = pw.read_section_raw(DBRPSection::NAMES_TABLE);
        const char* ptr = data.data();
        uint32_t count;
        std::memcpy(&count, ptr, sizeof(uint32_t));
        ptr += sizeof(uint32_t);

        for (uint32_t i = 0; i < count; i++) {
            uint32_t gid;
            std::memcpy(&gid, ptr, sizeof(uint32_t));
            ptr += sizeof(uint32_t);

            uint16_t name_len;
            std::memcpy(&name_len, ptr, sizeof(uint16_t));
            ptr += sizeof(uint16_t);

            std::string name(ptr, name_len);
            ptr += name_len;

            pw.names_map_[gid] = name;
        }
    }

    in.close();
    return pw;
}

const DBRPSectionEntry* DBRetinaPairwise::find_section(DBRPSection section_id) const {
    uint32_t sid = static_cast<uint32_t>(section_id);
    for (auto& entry : toc_) {
        if (entry.section_id == sid) return &entry;
    }
    return nullptr;
}

std::vector<char> DBRetinaPairwise::read_section_raw(DBRPSection section_id) const {
    const DBRPSectionEntry* entry = find_section(section_id);
    if (!entry) {
        throw std::runtime_error("Section " + std::to_string(static_cast<uint32_t>(section_id))
            + " not found in " + filepath_);
    }

    std::ifstream in(filepath_, std::ios::binary);
    in.seekg(entry->offset);
    std::vector<char> data(entry->length);
    in.read(data.data(), entry->length);
    in.close();

    // Verify CRC32
    uint32_t computed_crc = compute_crc32(data.data(), data.size());
    if (computed_crc != entry->crc32) {
        throw std::runtime_error("CRC32 mismatch in section " + std::to_string(static_cast<uint32_t>(section_id))
            + " of " + filepath_ + " (data corruption)");
    }

    return data;
}

std::string DBRetinaPairwise::get_group_name(uint32_t group_id) const {
    auto it = names_map_.find(group_id);
    if (it == names_map_.end()) {
        return "unknown_" + std::to_string(group_id);
    }
    return it->second;
}

std::string DBRetinaPairwise::get_metadata_json() const {
    auto data = read_section_raw(DBRPSection::METADATA);
    return std::string(data.begin(), data.end());
}

PairwiseStatistics DBRetinaPairwise::get_statistics() const {
    auto data = read_section_raw(DBRPSection::STATISTICS);
    const char* ptr = data.data();

    PairwiseStatistics stats;

    uint32_t num_metrics;
    std::memcpy(&num_metrics, ptr, sizeof(uint32_t));
    ptr += sizeof(uint32_t);

    for (uint32_t i = 0; i < num_metrics; i++) {
        MetricHistogram hist;
        std::memcpy(&hist.metric_id, ptr, sizeof(uint8_t));
        ptr += sizeof(uint8_t);

        uint32_t num_buckets;
        std::memcpy(&num_buckets, ptr, sizeof(uint32_t));
        ptr += sizeof(uint32_t);

        hist.bucket_counts.resize(num_buckets);
        for (uint32_t j = 0; j < num_buckets; j++) {
            std::memcpy(&hist.bucket_counts[j], ptr, sizeof(uint64_t));
            ptr += sizeof(uint64_t);
        }
        stats.histograms.push_back(std::move(hist));
    }

    std::memcpy(&stats.min_odds_ratio, ptr, sizeof(double));
    ptr += sizeof(double);
    std::memcpy(&stats.max_odds_ratio, ptr, sizeof(double));

    return stats;
}

std::string DBRetinaPairwise::get_statistics_json() const {
    auto stats = get_statistics();

    // Convert to JSON format matching the existing _pairwise_stats.json
    std::ostringstream json;
    json << "{\n";

    static const char* metric_names[] = {
        "containment", "ochiai", "jaccard", "csi", "dice"
    };
    static const char* bucket_labels[] = {
        "0-5", "5-10", "10-15", "15-20", "20-25", "25-30", "30-35", "35-40",
        "40-45", "45-50", "50-55", "55-60", "60-65", "65-70", "70-75", "75-80",
        "80-85", "85-90", "90-95", "95-100", "100-100"
    };

    for (size_t h = 0; h < stats.histograms.size(); h++) {
        auto& hist = stats.histograms[h];
        if (hist.metric_id >= 5) continue;  // Skip odds_ratio/pvalue histograms
        json << "\t\"" << metric_names[hist.metric_id] << "\": {\n";
        for (size_t b = 0; b < hist.bucket_counts.size() && b < 21; b++) {
            json << "\t\t\"" << bucket_labels[b] << "\": " << hist.bucket_counts[b];
            if (b < hist.bucket_counts.size() - 1 && b < 20) json << ",";
            json << "\n";
        }
        json << "\t}";
        if (h < stats.histograms.size() - 1) json << ",";
        json << "\n";
    }

    json << "}\n";
    return json.str();
}

void DBRetinaPairwise::load_group_index() const {
    if (group_index_loaded_) return;

    auto data = read_section_raw(DBRPSection::GROUP_INDEX);
    const char* ptr = data.data();

    uint32_t num_entries;
    std::memcpy(&num_entries, ptr, sizeof(uint32_t));
    ptr += sizeof(uint32_t);

    struct IndexEntry {
        uint32_t group_id;
        uint64_t list_offset;
        uint32_t list_count;
    };
    std::vector<IndexEntry> entries(num_entries);
    for (uint32_t i = 0; i < num_entries; i++) {
        std::memcpy(&entries[i].group_id, ptr, sizeof(uint32_t)); ptr += sizeof(uint32_t);
        std::memcpy(&entries[i].list_offset, ptr, sizeof(uint64_t)); ptr += sizeof(uint64_t);
        std::memcpy(&entries[i].list_count, ptr, sizeof(uint32_t)); ptr += sizeof(uint32_t);
    }

    // ptr now points to the start of the record lists area
    const char* lists_base = ptr;

    for (auto& entry : entries) {
        const char* list_ptr = lists_base + entry.list_offset;
        std::vector<uint64_t> indices(entry.list_count);
        for (uint32_t j = 0; j < entry.list_count; j++) {
            std::memcpy(&indices[j], list_ptr, sizeof(uint64_t));
            list_ptr += sizeof(uint64_t);
        }
        group_index_[entry.group_id] = std::move(indices);
    }

    group_index_loaded_ = true;
}

PairRecord DBRetinaPairwise::read_record(uint64_t record_index) const {
    const auto* section = find_section(DBRPSection::PAIR_RECORDS);
    if (!section) {
        throw std::runtime_error("PAIR_RECORDS section not found in " + filepath_);
    }

    uint64_t byte_offset = section->offset + record_index * record_size_;

    std::ifstream in(filepath_, std::ios::binary);
    in.seekg(byte_offset);
    std::vector<char> buffer(record_size_);
    in.read(buffer.data(), record_size_);
    in.close();

    return deserialize_record(buffer.data());
}

std::vector<PairRecord> DBRetinaPairwise::filter_pairs(uint8_t metric_id, double cutoff) const {
    const auto* section = find_section(DBRPSection::PAIR_RECORDS);
    if (!section) {
        throw std::runtime_error("PAIR_RECORDS section not found in " + filepath_);
    }

    std::vector<PairRecord> results;
    std::ifstream in(filepath_, std::ios::binary);
    in.seekg(section->offset);

    const size_t batch_size = 65536;
    std::vector<char> buffer(batch_size * record_size_);

    uint64_t remaining = num_pairs_;
    while (remaining > 0) {
        size_t to_read = std::min(static_cast<size_t>(remaining), batch_size);
        in.read(buffer.data(), to_read * record_size_);

        for (size_t i = 0; i < to_read; i++) {
            PairRecord rec = deserialize_record(buffer.data() + i * record_size_);
            float val = get_metric_value(rec, metric_id);
            if (val >= cutoff) {
                results.push_back(rec);
            }
        }
        remaining -= to_read;
    }

    in.close();
    return results;
}

std::vector<PairRecord> DBRetinaPairwise::query_group(uint32_t group_id, uint8_t metric_id, double cutoff) const {
    load_group_index();

    auto it = group_index_.find(group_id);
    if (it == group_index_.end()) {
        return {};
    }

    std::vector<PairRecord> results;
    for (uint64_t rec_idx : it->second) {
        PairRecord rec = read_record(rec_idx);
        float val = get_metric_value(rec, metric_id);
        if (val >= cutoff) {
            results.push_back(rec);
        }
    }

    return results;
}

std::vector<PairRecord> DBRetinaPairwise::iterate_all() const {
    return filter_pairs(0, 0.0);
}
