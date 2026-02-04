#include "ParquetPairwiseWriter.hpp"
#include <arrow/builder.h>
#include <arrow/table.h>
#include <arrow/type.h>
#include <arrow/io/file.h>
#include <parquet/arrow/writer.h>
#include <parquet/properties.h>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <algorithm>

namespace fs = std::filesystem;

// --- ParquetWriterStats ---

ParquetWriterStats::ParquetWriterStats()
    : min_odds_ratio(std::numeric_limits<double>::max()),
      max_odds_ratio(std::numeric_limits<double>::lowest()),
      num_pairs(0)
{
    histograms.resize(5);
    for (auto& h : histograms) {
        h.resize(21, 0);
    }
}

void ParquetWriterStats::update(const ParquetPairRecord& rec) {
    auto to_bucket = [](float value) -> int {
        int lower = static_cast<int>(std::floor(value / 5.0f)) * 5;
        int idx = lower / 5;
        if (idx >= 20) idx = 20;
        if (idx < 0) idx = 0;
        return idx;
    };

    histograms[0][to_bucket(rec.containment)]++;
    histograms[1][to_bucket(rec.ochiai)]++;
    histograms[2][to_bucket(rec.jaccard)]++;
    histograms[3][to_bucket(rec.csi)]++;
    histograms[4][to_bucket(rec.dice)]++;

    if (rec.odds_ratio < min_odds_ratio) min_odds_ratio = rec.odds_ratio;
    if (rec.odds_ratio > max_odds_ratio) max_odds_ratio = rec.odds_ratio;
    num_pairs++;
}

void ParquetWriterStats::merge(const ParquetWriterStats& other) {
    for (int m = 0; m < 5; m++) {
        for (int b = 0; b < 21; b++) {
            histograms[m][b] += other.histograms[m][b];
        }
    }
    if (other.min_odds_ratio < min_odds_ratio) min_odds_ratio = other.min_odds_ratio;
    if (other.max_odds_ratio > max_odds_ratio) max_odds_ratio = other.max_odds_ratio;
    num_pairs += other.num_pairs;
}

// --- ParquetPairwiseWriter ---

ParquetPairwiseWriter::ParquetPairwiseWriter(
    const std::string& output_dir,
    bool has_pvalue,
    int num_threads,
    const phmap::flat_hash_map<int, std::string>& names_map)
    : output_dir_(output_dir),
      has_pvalue_(has_pvalue),
      num_threads_(num_threads),
      names_map_(names_map)
{
    // Create output directory structure
    fs::create_directories(output_dir_ + "/data");

    // Initialize per-thread writers
    thread_writers_.resize(num_threads_);
    for (int i = 0; i < num_threads_; i++) {
        init_thread_writer(i);
    }
}

ParquetPairwiseWriter::~ParquetPairwiseWriter() = default;

std::shared_ptr<arrow::Schema> ParquetPairwiseWriter::make_schema() const {
    std::vector<std::shared_ptr<arrow::Field>> fields = {
        arrow::field("group_1_id", arrow::uint32()),
        arrow::field("group_2_id", arrow::uint32()),
        arrow::field("shared_features", arrow::uint64()),
        arrow::field("containment", arrow::float32()),
        arrow::field("ochiai", arrow::float32()),
        arrow::field("jaccard", arrow::float32()),
        arrow::field("csi", arrow::float32()),
        arrow::field("dice", arrow::float32()),
        arrow::field("odds_ratio", arrow::float32()),
    };
    if (has_pvalue_) {
        fields.push_back(arrow::field("pvalue", arrow::float64()));
    }
    return arrow::schema(fields);
}

void ParquetPairwiseWriter::init_thread_writer(int thread_id) {
    auto tw = std::make_unique<ThreadWriter>();
    tw->batch_count = 0;

    // Create builders
    tw->group_1_id_builder = std::make_shared<arrow::UInt32Builder>();
    tw->group_2_id_builder = std::make_shared<arrow::UInt32Builder>();
    tw->shared_features_builder = std::make_shared<arrow::UInt64Builder>();
    tw->containment_builder = std::make_shared<arrow::FloatBuilder>();
    tw->ochiai_builder = std::make_shared<arrow::FloatBuilder>();
    tw->jaccard_builder = std::make_shared<arrow::FloatBuilder>();
    tw->csi_builder = std::make_shared<arrow::FloatBuilder>();
    tw->dice_builder = std::make_shared<arrow::FloatBuilder>();
    tw->odds_ratio_builder = std::make_shared<arrow::FloatBuilder>();
    if (has_pvalue_) {
        tw->pvalue_builder = std::make_shared<arrow::DoubleBuilder>();
    }

    // Reserve capacity for batch
    tw->group_1_id_builder->Reserve(BATCH_SIZE).ok();
    tw->group_2_id_builder->Reserve(BATCH_SIZE).ok();
    tw->shared_features_builder->Reserve(BATCH_SIZE).ok();
    tw->containment_builder->Reserve(BATCH_SIZE).ok();
    tw->ochiai_builder->Reserve(BATCH_SIZE).ok();
    tw->jaccard_builder->Reserve(BATCH_SIZE).ok();
    tw->csi_builder->Reserve(BATCH_SIZE).ok();
    tw->dice_builder->Reserve(BATCH_SIZE).ok();
    tw->odds_ratio_builder->Reserve(BATCH_SIZE).ok();
    if (has_pvalue_) {
        tw->pvalue_builder->Reserve(BATCH_SIZE).ok();
    }

    // Open partition file
    std::ostringstream filename;
    filename << output_dir_ << "/data/part_" << std::setw(4) << std::setfill('0') << thread_id << ".parquet";

    auto result = arrow::io::FileOutputStream::Open(filename.str());
    if (!result.ok()) {
        throw std::runtime_error("Failed to open parquet file: " + filename.str() + " - " + result.status().ToString());
    }
    tw->file = *result;

    // Create Parquet writer with ZSTD compression
    auto schema = make_schema();

    auto parquet_props = parquet::WriterProperties::Builder()
        .compression(parquet::Compression::ZSTD)
        ->compression_level(3)  // Fast ZSTD level
        ->build();

    auto arrow_props = parquet::ArrowWriterProperties::Builder()
        .store_schema()
        ->build();

    auto writer_result = parquet::arrow::FileWriter::Open(
        *schema, arrow::default_memory_pool(), tw->file, parquet_props, arrow_props);
    if (!writer_result.ok()) {
        throw std::runtime_error("Failed to create parquet writer: " + writer_result.status().ToString());
    }
    tw->writer = std::move(*writer_result);

    thread_writers_[thread_id] = std::move(tw);
}

void ParquetPairwiseWriter::write_record(int thread_id, const ParquetPairRecord& rec) {
    auto& tw = *thread_writers_[thread_id];

    tw.group_1_id_builder->UnsafeAppend(rec.group_1_id);
    tw.group_2_id_builder->UnsafeAppend(rec.group_2_id);
    tw.shared_features_builder->UnsafeAppend(rec.shared_features);
    tw.containment_builder->UnsafeAppend(rec.containment);
    tw.ochiai_builder->UnsafeAppend(rec.ochiai);
    tw.jaccard_builder->UnsafeAppend(rec.jaccard);
    tw.csi_builder->UnsafeAppend(rec.csi);
    tw.dice_builder->UnsafeAppend(rec.dice);
    tw.odds_ratio_builder->UnsafeAppend(rec.odds_ratio);
    if (has_pvalue_ && rec.has_pvalue) {
        tw.pvalue_builder->UnsafeAppend(rec.pvalue);
    }

    tw.group_ids_seen.insert(rec.group_1_id);
    tw.group_ids_seen.insert(rec.group_2_id);

    tw.stats.update(rec);
    tw.batch_count++;

    if (tw.batch_count >= BATCH_SIZE) {
        flush_batch(tw);
    }
}

void ParquetPairwiseWriter::flush_batch(ThreadWriter& tw) {
    if (tw.batch_count == 0) return;

    std::vector<std::shared_ptr<arrow::Array>> arrays;

    auto finish = [](auto& builder) -> std::shared_ptr<arrow::Array> {
        auto result = builder->Finish();
        if (!result.ok()) {
            throw std::runtime_error("Failed to finish array: " + result.status().ToString());
        }
        return *result;
    };

    arrays.push_back(finish(tw.group_1_id_builder));
    arrays.push_back(finish(tw.group_2_id_builder));
    arrays.push_back(finish(tw.shared_features_builder));
    arrays.push_back(finish(tw.containment_builder));
    arrays.push_back(finish(tw.ochiai_builder));
    arrays.push_back(finish(tw.jaccard_builder));
    arrays.push_back(finish(tw.csi_builder));
    arrays.push_back(finish(tw.dice_builder));
    arrays.push_back(finish(tw.odds_ratio_builder));
    if (has_pvalue_) {
        arrays.push_back(finish(tw.pvalue_builder));
    }

    auto batch = arrow::RecordBatch::Make(make_schema(), tw.batch_count, arrays);
    auto table = arrow::Table::Make(make_schema(), batch->columns(), tw.batch_count);

    auto status = tw.writer->WriteTable(*table, tw.batch_count);
    if (!status.ok()) {
        throw std::runtime_error("Failed to write table: " + status.ToString());
    }

    // Re-reserve capacity for next batch
    tw.group_1_id_builder->Reserve(BATCH_SIZE).ok();
    tw.group_2_id_builder->Reserve(BATCH_SIZE).ok();
    tw.shared_features_builder->Reserve(BATCH_SIZE).ok();
    tw.containment_builder->Reserve(BATCH_SIZE).ok();
    tw.ochiai_builder->Reserve(BATCH_SIZE).ok();
    tw.jaccard_builder->Reserve(BATCH_SIZE).ok();
    tw.csi_builder->Reserve(BATCH_SIZE).ok();
    tw.dice_builder->Reserve(BATCH_SIZE).ok();
    tw.odds_ratio_builder->Reserve(BATCH_SIZE).ok();
    if (has_pvalue_) {
        tw.pvalue_builder->Reserve(BATCH_SIZE).ok();
    }

    tw.batch_count = 0;
}

void ParquetPairwiseWriter::finalize(
    const std::string& full_command,
    uint64_t population_size,
    const std::string& cutoff_metric,
    double cutoff_threshold)
{
    // Flush remaining batches and close all thread writers
    for (int i = 0; i < num_threads_; i++) {
        auto& tw = *thread_writers_[i];
        flush_batch(tw);

        auto status = tw.writer->Close();
        if (!status.ok()) {
            throw std::runtime_error("Failed to close parquet writer: " + status.ToString());
        }
        tw.file->Close().ok();

        // Merge stats
        merged_stats_.merge(tw.stats);
    }

    // Remove empty partition files
    for (int i = 0; i < num_threads_; i++) {
        if (thread_writers_[i]->stats.num_pairs == 0) {
            std::ostringstream filename;
            filename << output_dir_ << "/data/part_" << std::setw(4) << std::setfill('0') << i << ".parquet";
            fs::remove(filename.str());
        }
    }

    // Write sidecar metadata files
    write_manifest(full_command, population_size, cutoff_metric, cutoff_threshold);
    write_names_parquet();
    write_statistics_json();
    write_group_index_parquet();
}

void ParquetPairwiseWriter::write_manifest(
    const std::string& full_command,
    uint64_t population_size,
    const std::string& cutoff_metric,
    double cutoff_threshold)
{
    // Escape the command string for JSON
    std::string escaped_command;
    for (char c : full_command) {
        switch (c) {
            case '"':  escaped_command += "\\\""; break;
            case '\\': escaped_command += "\\\\"; break;
            case '\n': escaped_command += "\\n"; break;
            case '\r': escaped_command += "\\r"; break;
            case '\t': escaped_command += "\\t"; break;
            default:
                if (static_cast<unsigned char>(c) < 0x20) {
                    char buf[8];
                    std::snprintf(buf, sizeof(buf), "\\u%04x", static_cast<unsigned char>(c));
                    escaped_command += buf;
                } else {
                    escaped_command += c;
                }
        }
    }

    std::ofstream f(output_dir_ + "/manifest.json");
    f << "{\n"
      << "  \"version\": 1,\n"
      << "  \"format\": \"dbretina_pairwise_parquet\",\n"
      << "  \"num_pairs\": " << merged_stats_.num_pairs << ",\n"
      << "  \"num_groups\": " << names_map_.size() << ",\n"
      << "  \"population_size\": " << population_size << ",\n"
      << "  \"has_pvalue\": " << (has_pvalue_ ? "true" : "false") << ",\n"
      << "  \"cutoff_metric\": \"" << cutoff_metric << "\",\n"
      << "  \"cutoff_threshold\": " << cutoff_threshold << ",\n"
      << "  \"command\": \"" << escaped_command << "\",\n"
      << "  \"num_partitions\": " << num_threads_ << "\n"
      << "}\n";
    f.close();
}

void ParquetPairwiseWriter::write_names_parquet() {
    auto id_builder = std::make_shared<arrow::UInt32Builder>();
    auto name_builder = std::make_shared<arrow::StringBuilder>();

    // Sort by ID for deterministic output
    std::vector<std::pair<int, std::string>> sorted_names(names_map_.begin(), names_map_.end());
    std::sort(sorted_names.begin(), sorted_names.end());

    for (const auto& [id, name] : sorted_names) {
        id_builder->Append(static_cast<uint32_t>(id)).ok();
        name_builder->Append(name).ok();
    }

    auto id_array = *id_builder->Finish();
    auto name_array = *name_builder->Finish();

    auto schema = arrow::schema({
        arrow::field("group_id", arrow::uint32()),
        arrow::field("group_name", arrow::utf8()),
    });

    auto table = arrow::Table::Make(schema, {id_array, name_array});

    auto result = arrow::io::FileOutputStream::Open(output_dir_ + "/names.parquet");
    if (!result.ok()) {
        throw std::runtime_error("Failed to open names.parquet: " + result.status().ToString());
    }

    auto props = parquet::WriterProperties::Builder()
        .compression(parquet::Compression::ZSTD)
        ->build();

    auto status = parquet::arrow::WriteTable(*table, arrow::default_memory_pool(), *result, 65536, props);
    if (!status.ok()) {
        throw std::runtime_error("Failed to write names.parquet: " + status.ToString());
    }
}

void ParquetPairwiseWriter::write_statistics_json() {
    std::vector<std::string> metric_names = {"containment", "ochiai", "jaccard", "csi", "dice"};

    std::ofstream f(output_dir_ + "/statistics.json");
    f << "{\n";
    for (int m = 0; m < 5; m++) {
        f << "  \"" << metric_names[m] << "\": {\n";
        for (int b = 0; b < 21; b++) {
            int lower = b * 5;
            int upper = (b == 20) ? 100 : lower + 5;
            f << "    \"" << lower << "-" << upper << "\": " << merged_stats_.histograms[m][b];
            if (b < 20) f << ",";
            f << "\n";
        }
        f << "  }";
        if (m < 4) f << ",";
        f << "\n";
    }
    f << "}\n";
    f.close();

    // Also write odds ratio stats
    std::ofstream f2(output_dir_ + "/statistics_odds_ratio.txt");
    f2 << merged_stats_.min_odds_ratio << "\n";
    f2 << merged_stats_.max_odds_ratio << "\n";
    f2.close();
}

void ParquetPairwiseWriter::write_group_index_parquet() {
    // Build group_id -> [partition_ids] mapping
    phmap::flat_hash_map<uint32_t, std::vector<int>> group_to_partitions;

    for (int i = 0; i < num_threads_; i++) {
        if (thread_writers_[i]->stats.num_pairs == 0) continue;
        for (uint32_t gid : thread_writers_[i]->group_ids_seen) {
            group_to_partitions[gid].push_back(i);
        }
    }

    // Sort and deduplicate partition lists
    for (auto& [gid, parts] : group_to_partitions) {
        std::sort(parts.begin(), parts.end());
        parts.erase(std::unique(parts.begin(), parts.end()), parts.end());
    }

    // Write as Parquet: group_id (uint32), partition_ids (list<int32>)
    auto gid_builder = std::make_shared<arrow::UInt32Builder>();
    auto list_builder = std::make_shared<arrow::ListBuilder>(
        arrow::default_memory_pool(),
        std::make_shared<arrow::Int32Builder>());
    auto value_builder = static_cast<arrow::Int32Builder*>(list_builder->value_builder());

    // Sort by group ID for deterministic output
    std::vector<std::pair<uint32_t, std::vector<int>>> sorted(
        group_to_partitions.begin(), group_to_partitions.end());
    std::sort(sorted.begin(), sorted.end());

    for (const auto& [gid, parts] : sorted) {
        gid_builder->Append(gid).ok();
        list_builder->Append().ok();
        for (int p : parts) {
            value_builder->Append(p).ok();
        }
    }

    auto gid_array = *gid_builder->Finish();
    auto list_array = *list_builder->Finish();

    auto schema = arrow::schema({
        arrow::field("group_id", arrow::uint32()),
        arrow::field("partition_ids", arrow::list(arrow::int32())),
    });

    auto table = arrow::Table::Make(schema, {gid_array, list_array});

    auto result = arrow::io::FileOutputStream::Open(output_dir_ + "/group_index.parquet");
    if (!result.ok()) {
        throw std::runtime_error("Failed to open group_index.parquet: " + result.status().ToString());
    }

    auto props = parquet::WriterProperties::Builder()
        .compression(parquet::Compression::ZSTD)
        ->build();

    auto status = parquet::arrow::WriteTable(*table, arrow::default_memory_pool(), *result, 65536, props);
    if (!status.ok()) {
        throw std::runtime_error("Failed to write group_index.parquet: " + status.ToString());
    }
}
