#include "DBRetinaIndex.hpp"
#include <sstream>
#include <algorithm>
#include <cstring>

// ---------------------------------------------------------------------------
// CRC32 (standard IEEE polynomial 0xEDB88320)
// ---------------------------------------------------------------------------
static const uint32_t crc32_table[256] = {
    0x00000000, 0x77073096, 0xEE0E612C, 0x990951BA, 0x076DC419, 0x706AF48F,
    0xE963A535, 0x9E6495A3, 0x0EDB8832, 0x79DCB8A4, 0xE0D5E91B, 0x97D2D988,
    0x09B64C2B, 0x7EB17CBF, 0xE7B82D09, 0x90BF1D91, 0x1DB71064, 0x6AB020F2,
    0xF3B97148, 0x84BE41DE, 0x1ADAD47D, 0x6DDDE4EB, 0xF4D4B551, 0x83D385C7,
    0x136C9856, 0x646BA8C0, 0xFD62F97A, 0x8A65C9EC, 0x14015C4F, 0x63066CD9,
    0xFA0F3D63, 0x8D080DF5, 0x3B6E20C8, 0x4C69105E, 0xD56041E4, 0xA2677172,
    0x3C03E4D1, 0x4B04D447, 0xD20D85FD, 0xA50AB56B, 0x35B5A8FA, 0x42B2986C,
    0xDBBBB9D6, 0xACBCF940, 0x32D86CE3, 0x45DF5C75, 0xDCD60DCF, 0xABD13D59,
    0x26D930AC, 0x51DE003A, 0xC8D75180, 0xBFD06116, 0x21B4F0B5, 0x56B3C423,
    0xCFBA9599, 0xB8BDB111, 0x28D7D0C8, 0x5FD0E09C, 0xC5B04628, 0xB2B76F9E,
    0x0CB61B38, 0x7BB12BAE, 0xE2B23F14, 0x95B53A82, 0x09282166, 0x7E272516,
    0xE92E7507, 0x9E290591, 0x0ABCF02E, 0x7DBBF0B8, 0xE4B2D002, 0x93B5A0D4,
    0x03EB532E, 0x74EC4D0B, 0xEBE561F5, 0x9CE85163, 0x0E856A42, 0x798B7AD4,
    0xE0824A6E, 0x97852EF8, 0x09C14F5B, 0x7EC6BFCD, 0xE7CF8E77, 0x90C8BFE1,
    0x1085A830, 0x678598A6, 0xFE8C691C, 0x898C598A, 0x17E84C29, 0x60EF7CBF,
    0xF9E66D05, 0x8EE15D93, 0x1EBD3E02, 0x69BA0E94, 0xF0B3112E, 0x87B436B8,
    0x19D0021B, 0x6ED7128D, 0xF7DE0437, 0x80D93425, 0x3AB551CE, 0x4DB26158,
    0xD4BB30E2, 0xA3BC0074, 0x35F80CB5, 0x42FF3C23, 0xDBF60A99, 0xACF10A0F,
    0x3CEF049E, 0x4BE87408, 0xD2E136B2, 0xA5E61624, 0x39029287, 0x4E050211,
    0xD70CA3AB, 0xA00A933D, 0x381B49AC, 0x4F1C793A, 0xD6156D80, 0xA1124E16,
    0x34E10AB5, 0x43E60D23, 0xDA0F1499, 0xAD082498, 0x3A06613B, 0x4D01A9AD,
    0xD408E917, 0xA30FD981, 0x130C963A, 0x640B06AC, 0xFD020016, 0x8A050580,
    0x14412489, 0x6346141F, 0xFA4F13A5, 0x8D483333, 0x1B0C9AA2, 0x6C0BA834,
    0xF502B98E, 0x8205C918, 0x1E016DBB, 0x6906DC2D, 0xF00F9397, 0x8708A3E1,
    0x1E01F268, 0x6906C2FE, 0xF762575D, 0x806567CB, 0x196C3671, 0x6E6B06E7,
    0xFED41B76, 0x89D32BE0, 0x10DA7A5A, 0x67DD4ACC, 0xF9B9DF6F, 0x8EBEEFF9,
    0x17B7BE43, 0x60B08ED5, 0xD6D6A3E8, 0xA1D1937E, 0x38D8C2C4, 0x4FDFF252,
    0xD1BB67F1, 0xA6BC5767, 0x3FB506DD, 0x48B2364B, 0xD80D2BDA, 0xAF0A1B4C,
    0x36034AF6, 0x41047A60, 0xDA0D0FC3, 0xAD0A1F55, 0x340304EF, 0x43040B79,
    0xBB0B4703, 0xCC0C7795, 0x5505262F, 0x220206B9, 0xBC66831A, 0xCB61B38C,
    0x5268E236, 0x256CD2A0, 0xB5D0CF31, 0xC2D7FFA7, 0x5BDEAE1D, 0x2CD99E8B,
    0xB40BBE37, 0xC30C8EA1, 0x5A05DF1B, 0x2D02EF8D, 0xB4014A84, 0xC3067A12,
    0x5A6C45A8, 0x2D6B753E, 0xB425779D, 0xC322470B, 0x5A2836B1, 0x2D2F0827,
    0xBA684DB6, 0xCD6F7D20, 0x540E669A, 0x2309560C, 0xBD2D6DAF, 0xCA2A3D39,
    0x532D0C83, 0x24220C15, 0xB43F025E, 0xC33832C8, 0x5A316372, 0x2D3653E4,
    0xB3667A2E, 0xC4614AB8, 0x5D681B02, 0x2A6F2B94, 0xB40BBE37, 0xC30C8EA1,
    0x5A05DF1B, 0x2D02EF8D, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
    0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
    0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
    0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
    0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
    0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
    0x00000000, 0x00000000, 0x00000000, 0x00000000,
};

uint32_t compute_crc32(const char* data, size_t length) {
    uint32_t crc = 0xFFFFFFFF;
    for (size_t i = 0; i < length; i++) {
        crc = crc32_table[(crc ^ static_cast<uint8_t>(data[i])) & 0xFF] ^ (crc >> 8);
    }
    return crc ^ 0xFFFFFFFF;
}

// ---------------------------------------------------------------------------
// Writing
// ---------------------------------------------------------------------------

void DBRetinaIndex::begin_write(const std::string& filepath) {
    filepath_ = filepath;
    writing_ = true;
    write_toc_.clear();

    write_stream_.open(filepath, std::ios::binary);
    if (!write_stream_.is_open()) {
        throw std::runtime_error("Cannot open file for writing: " + filepath);
    }

    // Write header: magic + version + placeholder for TOC offset
    write_stream_.write(DBRI_MAGIC, 4);
    write_stream_.write(reinterpret_cast<const char*>(&DBRI_FORMAT_VERSION), sizeof(uint32_t));
    uint64_t toc_offset_placeholder = 0;
    write_stream_.write(reinterpret_cast<const char*>(&toc_offset_placeholder), sizeof(uint64_t));
}

void DBRetinaIndex::write_section(DBRISection section_id, const char* data, size_t length) {
    if (!writing_) throw std::runtime_error("Not in write mode");

    SectionEntry entry;
    entry.section_id = static_cast<uint32_t>(section_id);
    entry.offset = static_cast<uint64_t>(write_stream_.tellp());
    entry.length = length;
    entry.crc32 = compute_crc32(data, length);

    write_stream_.write(data, length);
    write_toc_.push_back(entry);
}

void DBRetinaIndex::write_phmap(DBRetina_PHMAP* frame) {
    // Serialize PHMAP to a temporary file, then read it back
    std::string tmp_path = filepath_ + ".tmp_phmap";
    {
        phmap::BinaryOutputArchive ar_out(tmp_path.c_str());
        frame->MAP.phmap_dump(ar_out);
    }
    // Read the temp file
    std::ifstream tmp_in(tmp_path, std::ios::binary | std::ios::ate);
    size_t size = tmp_in.tellg();
    tmp_in.seekg(0);
    std::vector<char> buffer(size);
    tmp_in.read(buffer.data(), size);
    tmp_in.close();
    std::remove(tmp_path.c_str());

    // Write hashSize before the phmap data
    std::vector<char> section_data;
    section_data.resize(sizeof(uint64_t) + size);
    std::memcpy(section_data.data(), &frame->hashSize, sizeof(uint64_t));
    std::memcpy(section_data.data() + sizeof(uint64_t), buffer.data(), size);

    write_section(DBRISection::PHMAP, section_data.data(), section_data.size());
}

void DBRetinaIndex::write_color_to_sources(
    const phmap::flat_hash_map<uint64_t, std::vector<uint32_t>>& legend,
    const phmap::flat_hash_map<uint64_t, uint64_t>& colorsCount) {

    // Write ALL legend entries (including zero-count) for incremental append support
    std::string tmp_path = filepath_ + ".tmp_c2s";
    {
        phmap::BinaryOutputArchive ar_out(tmp_path.c_str());
        size_t count = legend.size();
        ar_out.saveBinary(count);
        for (auto& [color, sources] : legend) {
            ar_out.saveBinary(color);
            phmap::flat_hash_set<uint32_t> src_set(sources.begin(), sources.end());
            ar_out.saveBinary(src_set);
        }
    }
    std::ifstream tmp_in(tmp_path, std::ios::binary | std::ios::ate);
    size_t size = tmp_in.tellg();
    tmp_in.seekg(0);
    std::vector<char> buffer(size);
    tmp_in.read(buffer.data(), size);
    tmp_in.close();
    std::remove(tmp_path.c_str());

    write_section(DBRISection::COLOR_TO_SOURCES, buffer.data(), buffer.size());
}

void DBRetinaIndex::write_color_count(const phmap::flat_hash_map<uint64_t, uint64_t>& colorsCount) {
    std::string tmp_path = filepath_ + ".tmp_cc";
    {
        phmap::BinaryOutputArchive ar_out(tmp_path.c_str());
        colorsCount.phmap_dump(ar_out);
    }
    std::ifstream tmp_in(tmp_path, std::ios::binary | std::ios::ate);
    size_t size = tmp_in.tellg();
    tmp_in.seekg(0);
    std::vector<char> buffer(size);
    tmp_in.read(buffer.data(), size);
    tmp_in.close();
    std::remove(tmp_path.c_str());

    write_section(DBRISection::COLOR_COUNT, buffer.data(), buffer.size());
}

void DBRetinaIndex::write_group_feature_count(const phmap::flat_hash_map<uint32_t, uint32_t>& groupID_to_featureCount) {
    std::string tmp_path = filepath_ + ".tmp_gfc";
    {
        phmap::BinaryOutputArchive ar_out(tmp_path.c_str());
        groupID_to_featureCount.phmap_dump(ar_out);
    }
    std::ifstream tmp_in(tmp_path, std::ios::binary | std::ios::ate);
    size_t size = tmp_in.tellg();
    tmp_in.seekg(0);
    std::vector<char> buffer(size);
    tmp_in.read(buffer.data(), size);
    tmp_in.close();
    std::remove(tmp_path.c_str());

    write_section(DBRISection::GROUP_FEATURE_COUNT, buffer.data(), buffer.size());
}

void DBRetinaIndex::write_names_map(
    const phmap::flat_hash_map<std::string, std::string>& namesMap,
    const phmap::flat_hash_map<std::string, uint64_t>& groupNameMap) {
    // Serialize as: count, then (groupID as uint32_t, name_length as uint32_t, name_bytes)
    std::ostringstream oss(std::ios::binary);
    uint32_t count = static_cast<uint32_t>(namesMap.size());
    oss.write(reinterpret_cast<const char*>(&count), sizeof(uint32_t));

    for (auto& [seq_name, group_name] : namesMap) {
        auto it = groupNameMap.find(group_name);
        if (it == groupNameMap.end()) continue;
        uint32_t gid = static_cast<uint32_t>(it->second);
        uint32_t name_len = static_cast<uint32_t>(group_name.size());
        oss.write(reinterpret_cast<const char*>(&gid), sizeof(uint32_t));
        oss.write(reinterpret_cast<const char*>(&name_len), sizeof(uint32_t));
        oss.write(group_name.data(), name_len);
    }

    std::string data = oss.str();
    write_section(DBRISection::NAMES_MAP, data.data(), data.size());
}

void DBRetinaIndex::write_metadata(const std::string& metadata_json) {
    write_section(DBRISection::METADATA, metadata_json.data(), metadata_json.size());
}

void DBRetinaIndex::write_raw_gene_sets(const std::string& json_content) {
    write_section(DBRISection::RAW_GENE_SETS, json_content.data(), json_content.size());
}

void DBRetinaIndex::write_hashed_gene_sets(const std::string& json_content) {
    write_section(DBRISection::HASHED_GENE_SETS, json_content.data(), json_content.size());
}

void DBRetinaIndex::write_tags_map(const phmap::flat_hash_map<std::string, uint64_t>& tagsMap) {
    // Serialize as: count, then (key_length, key_bytes, value)
    std::ostringstream oss(std::ios::binary);
    uint64_t count = tagsMap.size();
    oss.write(reinterpret_cast<const char*>(&count), sizeof(uint64_t));

    for (auto& [key, value] : tagsMap) {
        uint32_t key_len = static_cast<uint32_t>(key.size());
        oss.write(reinterpret_cast<const char*>(&key_len), sizeof(uint32_t));
        oss.write(key.data(), key_len);
        oss.write(reinterpret_cast<const char*>(&value), sizeof(uint64_t));
    }

    std::string data = oss.str();
    write_section(DBRISection::TAGS_MAP, data.data(), data.size());
}

void DBRetinaIndex::write_free_colors(
    const std::priority_queue<uint64_t, std::vector<uint64_t>, std::greater<uint64_t>>& freeColors) {
    // Copy the priority queue to extract elements
    auto copy = freeColors;
    std::vector<uint64_t> colors;
    while (!copy.empty()) {
        colors.push_back(copy.top());
        copy.pop();
    }

    std::ostringstream oss(std::ios::binary);
    uint64_t count = colors.size();
    oss.write(reinterpret_cast<const char*>(&count), sizeof(uint64_t));
    for (auto& c : colors) {
        oss.write(reinterpret_cast<const char*>(&c), sizeof(uint64_t));
    }

    std::string data = oss.str();
    write_section(DBRISection::FREE_COLORS, data.data(), data.size());
}

void DBRetinaIndex::write_group_to_feature_set(
    const phmap::parallel_flat_hash_map<std::string, phmap::flat_hash_set<uint64_t>>& groupName_to_featureSet) {
    // Serialize as: count, then for each: (name_len, name, set_size, hash1, hash2, ...)
    std::ostringstream oss(std::ios::binary);
    uint64_t count = groupName_to_featureSet.size();
    oss.write(reinterpret_cast<const char*>(&count), sizeof(uint64_t));

    for (auto& [name, feature_set] : groupName_to_featureSet) {
        uint32_t name_len = static_cast<uint32_t>(name.size());
        oss.write(reinterpret_cast<const char*>(&name_len), sizeof(uint32_t));
        oss.write(name.data(), name_len);

        uint64_t set_size = feature_set.size();
        oss.write(reinterpret_cast<const char*>(&set_size), sizeof(uint64_t));
        for (auto& hash : feature_set) {
            oss.write(reinterpret_cast<const char*>(&hash), sizeof(uint64_t));
        }
    }

    std::string data = oss.str();
    write_section(DBRISection::GROUP_TO_FEATURE_SET, data.data(), data.size());
}

void DBRetinaIndex::finalize_write() {
    if (!writing_) throw std::runtime_error("Not in write mode");

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

    // Go back and write TOC offset in header
    write_stream_.seekp(4 + sizeof(uint32_t));  // after magic + version
    write_stream_.write(reinterpret_cast<const char*>(&toc_offset), sizeof(uint64_t));

    write_stream_.close();
    writing_ = false;
    toc_ = write_toc_;
}

// ---------------------------------------------------------------------------
// Reading
// ---------------------------------------------------------------------------

DBRetinaIndex DBRetinaIndex::open(const std::string& filepath) {
    DBRetinaIndex idx;
    idx.filepath_ = filepath;

    std::ifstream in(filepath, std::ios::binary);
    if (!in.is_open()) {
        throw std::runtime_error("Cannot open .dbri file: " + filepath);
    }

    // Read and verify magic
    char magic[4];
    in.read(magic, 4);
    if (std::memcmp(magic, DBRI_MAGIC, 4) != 0) {
        throw std::runtime_error("Invalid .dbri file (bad magic bytes): " + filepath);
    }

    // Read version
    uint32_t version;
    in.read(reinterpret_cast<char*>(&version), sizeof(uint32_t));
    if (version != DBRI_FORMAT_VERSION) {
        throw std::runtime_error("Unsupported .dbri format version " + std::to_string(version)
            + " (expected " + std::to_string(DBRI_FORMAT_VERSION) + ")");
    }

    // Read TOC offset
    uint64_t toc_offset;
    in.read(reinterpret_cast<char*>(&toc_offset), sizeof(uint64_t));
    if (!in) {
        throw std::runtime_error("Incomplete/corrupt .dbri file (truncated header): " + filepath);
    }

    // Validate the TOC offset before trusting it. The header is exactly
    // 16 bytes (magic[4] + version[4] + toc_offset[8]). An unfinalized index
    // (written by begin_write but never finalize_write, e.g. a crash mid-write)
    // still carries the placeholder toc_offset == 0; a truncated file points the
    // offset past EOF. Either way, seeking there and reading num_sections would
    // interpret garbage (the magic bytes read as ~1.2B sections) and silently
    // yield a bogus empty index. Reject both cleanly instead.
    constexpr uint64_t HEADER_SIZE = 4 + sizeof(uint32_t) + sizeof(uint64_t);  // 16
    in.seekg(0, std::ios::end);
    uint64_t file_size = static_cast<uint64_t>(in.tellg());
    // The TOC must start after the header and leave room for at least the
    // 4-byte num_sections field.
    if (toc_offset < HEADER_SIZE || toc_offset + sizeof(uint32_t) > file_size) {
        throw std::runtime_error(
            "Incomplete/corrupt .dbri file (invalid TOC offset " + std::to_string(toc_offset)
            + ", file size " + std::to_string(file_size)
            + "); the index was likely not finalized: " + filepath);
    }

    // Seek to TOC and read it
    in.seekg(toc_offset);
    uint32_t num_sections;
    in.read(reinterpret_cast<char*>(&num_sections), sizeof(uint32_t));
    if (!in) {
        throw std::runtime_error("Incomplete/corrupt .dbri file (cannot read TOC): " + filepath);
    }

    idx.toc_.resize(num_sections);
    for (uint32_t i = 0; i < num_sections; i++) {
        in.read(reinterpret_cast<char*>(&idx.toc_[i].section_id), sizeof(uint32_t));
        in.read(reinterpret_cast<char*>(&idx.toc_[i].offset), sizeof(uint64_t));
        in.read(reinterpret_cast<char*>(&idx.toc_[i].length), sizeof(uint64_t));
        in.read(reinterpret_cast<char*>(&idx.toc_[i].crc32), sizeof(uint32_t));
    }
    if (!in) {
        throw std::runtime_error("Incomplete/corrupt .dbri file (truncated TOC entries): " + filepath);
    }

    in.close();
    return idx;
}

bool DBRetinaIndex::has_section(DBRISection section_id) const {
    return find_section(section_id) != nullptr;
}

const SectionEntry* DBRetinaIndex::find_section(DBRISection section_id) const {
    uint32_t sid = static_cast<uint32_t>(section_id);
    for (auto& entry : toc_) {
        if (entry.section_id == sid) return &entry;
    }
    return nullptr;
}

std::vector<char> DBRetinaIndex::read_section_raw(DBRISection section_id) const {
    const SectionEntry* entry = find_section(section_id);
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

DBRetina_PHMAP* DBRetinaIndex::load_phmap() const {
    auto data = read_section_raw(DBRISection::PHMAP);

    // First 8 bytes = hashSize
    uint64_t hashSize;
    std::memcpy(&hashSize, data.data(), sizeof(uint64_t));

    // Write remaining data to temp file for phmap_load
    std::string tmp_path = filepath_ + ".tmp_phmap_load";
    {
        std::ofstream tmp_out(tmp_path, std::ios::binary);
        tmp_out.write(data.data() + sizeof(uint64_t), data.size() - sizeof(uint64_t));
    }

    auto* frame = new DBRetina_PHMAP(hashSize);
    phmap::BinaryInputArchive ar_in(tmp_path.c_str());
    frame->MAP.phmap_load(ar_in);
    std::remove(tmp_path.c_str());

    return frame;
}

void DBRetinaIndex::load_color_to_sources(
    phmap::parallel_flat_hash_map<uint32_t, std::vector<uint32_t>,
        std::hash<uint32_t>, std::equal_to<uint32_t>,
        std::allocator<std::pair<const uint32_t, std::vector<uint32_t>>>, 1>& color_to_ids) const {

    auto data = read_section_raw(DBRISection::COLOR_TO_SOURCES);

    // Write to temp file for BinaryInputArchive
    std::string tmp_path = filepath_ + ".tmp_c2s_load";
    {
        std::ofstream tmp_out(tmp_path, std::ios::binary);
        tmp_out.write(data.data(), data.size());
    }

    phmap::BinaryInputArchive ar_in(tmp_path.c_str());
    size_t size;
    ar_in.loadBinary(&size);
    color_to_ids.reserve(size);
    while (size--) {
        uint64_t k;
        phmap::flat_hash_set<uint32_t> v;
        ar_in.loadBinary(&k);
        ar_in.loadBinary(&v);
        std::vector<uint32_t> vVec(v.begin(), v.end());
        color_to_ids.insert_or_assign(static_cast<uint32_t>(k), std::move(vVec));
    }
    std::remove(tmp_path.c_str());
}

void DBRetinaIndex::load_color_count(
    phmap::parallel_flat_hash_map<uint32_t, uint32_t,
        std::hash<uint32_t>, std::equal_to<uint32_t>,
        std::allocator<std::pair<const uint32_t, uint32_t>>, 1>& colorsCount) const {

    auto data = read_section_raw(DBRISection::COLOR_COUNT);

    std::string tmp_path = filepath_ + ".tmp_cc_load";
    {
        std::ofstream tmp_out(tmp_path, std::ios::binary);
        tmp_out.write(data.data(), data.size());
    }

    phmap::flat_hash_map<uint64_t, uint64_t> tmpMap;
    phmap::BinaryInputArchive ar_in(tmp_path.c_str());
    tmpMap.phmap_load(ar_in);
    for (auto& [k, v] : tmpMap) {
        colorsCount.insert_or_assign(static_cast<uint32_t>(k), static_cast<uint32_t>(v));
    }
    std::remove(tmp_path.c_str());
}

void DBRetinaIndex::load_group_feature_count(phmap::flat_hash_map<uint32_t, uint32_t>& groupID_to_featureCount) const {
    auto data = read_section_raw(DBRISection::GROUP_FEATURE_COUNT);

    std::string tmp_path = filepath_ + ".tmp_gfc_load";
    {
        std::ofstream tmp_out(tmp_path, std::ios::binary);
        tmp_out.write(data.data(), data.size());
    }

    phmap::BinaryInputArchive ar_in(tmp_path.c_str());
    groupID_to_featureCount.phmap_load(ar_in);
    std::remove(tmp_path.c_str());
}

void DBRetinaIndex::load_names_map(phmap::flat_hash_map<int, std::string>& namesMap) const {
    auto data = read_section_raw(DBRISection::NAMES_MAP);

    const char* ptr = data.data();
    uint32_t count;
    std::memcpy(&count, ptr, sizeof(uint32_t));
    ptr += sizeof(uint32_t);

    for (uint32_t i = 0; i < count; i++) {
        uint32_t gid;
        std::memcpy(&gid, ptr, sizeof(uint32_t));
        ptr += sizeof(uint32_t);

        uint32_t name_len;
        std::memcpy(&name_len, ptr, sizeof(uint32_t));
        ptr += sizeof(uint32_t);

        std::string name(ptr, name_len);
        ptr += name_len;

        namesMap[static_cast<int>(gid)] = name;
    }
}

std::string DBRetinaIndex::load_metadata() const {
    auto data = read_section_raw(DBRISection::METADATA);
    return std::string(data.begin(), data.end());
}

std::string DBRetinaIndex::load_raw_gene_sets() const {
    auto data = read_section_raw(DBRISection::RAW_GENE_SETS);
    return std::string(data.begin(), data.end());
}

std::string DBRetinaIndex::load_hashed_gene_sets() const {
    auto data = read_section_raw(DBRISection::HASHED_GENE_SETS);
    return std::string(data.begin(), data.end());
}

void DBRetinaIndex::load_tags_map(phmap::flat_hash_map<std::string, uint64_t>& tagsMap) const {
    auto data = read_section_raw(DBRISection::TAGS_MAP);

    const char* ptr = data.data();
    uint64_t count;
    std::memcpy(&count, ptr, sizeof(uint64_t));
    ptr += sizeof(uint64_t);

    for (uint64_t i = 0; i < count; i++) {
        uint32_t key_len;
        std::memcpy(&key_len, ptr, sizeof(uint32_t));
        ptr += sizeof(uint32_t);

        std::string key(ptr, key_len);
        ptr += key_len;

        uint64_t value;
        std::memcpy(&value, ptr, sizeof(uint64_t));
        ptr += sizeof(uint64_t);

        tagsMap[key] = value;
    }
}

void DBRetinaIndex::load_free_colors(
    std::priority_queue<uint64_t, std::vector<uint64_t>, std::greater<uint64_t>>& freeColors) const {
    auto data = read_section_raw(DBRISection::FREE_COLORS);

    const char* ptr = data.data();
    uint64_t count;
    std::memcpy(&count, ptr, sizeof(uint64_t));
    ptr += sizeof(uint64_t);

    for (uint64_t i = 0; i < count; i++) {
        uint64_t color;
        std::memcpy(&color, ptr, sizeof(uint64_t));
        ptr += sizeof(uint64_t);
        freeColors.push(color);
    }
}

void DBRetinaIndex::load_group_to_feature_set(
    phmap::parallel_flat_hash_map<std::string, phmap::flat_hash_set<uint64_t>>& groupName_to_featureSet) const {
    auto data = read_section_raw(DBRISection::GROUP_TO_FEATURE_SET);

    const char* ptr = data.data();
    uint64_t count;
    std::memcpy(&count, ptr, sizeof(uint64_t));
    ptr += sizeof(uint64_t);

    for (uint64_t i = 0; i < count; i++) {
        uint32_t name_len;
        std::memcpy(&name_len, ptr, sizeof(uint32_t));
        ptr += sizeof(uint32_t);

        std::string name(ptr, name_len);
        ptr += name_len;

        uint64_t set_size;
        std::memcpy(&set_size, ptr, sizeof(uint64_t));
        ptr += sizeof(uint64_t);

        phmap::flat_hash_set<uint64_t> feature_set;
        feature_set.reserve(set_size);
        for (uint64_t j = 0; j < set_size; j++) {
            uint64_t hash;
            std::memcpy(&hash, ptr, sizeof(uint64_t));
            ptr += sizeof(uint64_t);
            feature_set.insert(hash);
        }

        groupName_to_featureSet[name] = std::move(feature_set);
    }
}

// Metadata helpers
uint64_t DBRetinaIndex::get_population_size() const {
    // Parse from metadata JSON - look for "population_size" field
    std::string meta = load_metadata();
    // Simple parsing: find "population_size":NUMBER
    std::string key = "\"population_size\":";
    auto pos = meta.find(key);
    if (pos == std::string::npos) {
        throw std::runtime_error("population_size not found in .dbri metadata");
    }
    pos += key.size();
    // Skip whitespace
    while (pos < meta.size() && (meta[pos] == ' ' || meta[pos] == '\t')) pos++;
    std::string num_str;
    while (pos < meta.size() && std::isdigit(meta[pos])) {
        num_str += meta[pos++];
    }
    return std::stoull(num_str);
}

uint32_t DBRetinaIndex::get_num_groups() const {
    std::string meta = load_metadata();
    std::string key = "\"num_groups\":";
    auto pos = meta.find(key);
    if (pos == std::string::npos) {
        throw std::runtime_error("num_groups not found in .dbri metadata");
    }
    pos += key.size();
    while (pos < meta.size() && (meta[pos] == ' ' || meta[pos] == '\t')) pos++;
    std::string num_str;
    while (pos < meta.size() && std::isdigit(meta[pos])) {
        num_str += meta[pos++];
    }
    return std::stoul(num_str);
}

uint64_t DBRetinaIndex::get_max_group_id() const {
    std::string meta = load_metadata();
    std::string key = "\"max_group_id\":";
    auto pos = meta.find(key);
    if (pos == std::string::npos) {
        throw std::runtime_error("max_group_id not found in .dbri metadata");
    }
    pos += key.size();
    while (pos < meta.size() && (meta[pos] == ' ' || meta[pos] == '\t')) pos++;
    std::string num_str;
    while (pos < meta.size() && std::isdigit(meta[pos])) {
        num_str += meta[pos++];
    }
    return std::stoull(num_str);
}
