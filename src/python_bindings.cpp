#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/unordered_map.h>
#include <nanobind/stl/map.h>
#include <nanobind/stl/vector.h>

#include "dbretina_core.hpp"
#include "DBRetina.hpp"
#include "DBRetina_GSA.hpp"
#include "DBRetinaIndex.hpp"
#include "DBRetinaPairwise.hpp"
#include <nanobind/stl/pair.h>
#include <nanobind/stl/list.h>

namespace nb = nanobind;

NB_MODULE(_dbretina_internal, m) {
    m.doc() = "DBRetina C++ bindings";

    // dbretina namespace functions
    m.def("pairwise", &dbretina::pairwise,
          nb::arg("index_prefix"),
          nb::arg("user_threads"),
          nb::arg("cutoff_distance_type"),
          nb::arg("cutoff_threshold"),
          nb::arg("full_command"),
          nb::arg("calculate_pvalue"));

    m.def("dbretina_indexing", &dbretina::dbretina_indexing,
          nb::arg("json_file"),
          nb::arg("user_index_prefix"));

    m.def("dbretina_append", &dbretina::dbretina_append,
          nb::arg("existing_dbri_path"),
          nb::arg("new_json_file"),
          nb::arg("new_raw_json_file"),
          nb::arg("output_dbri_path"));

    m.def("dbretina_merge", &dbretina::dbretina_merge,
          nb::arg("index_a_path"),
          nb::arg("index_b_path"),
          nb::arg("output_dbri_path"));

    // .dbri extraction helpers for Python consumers
    m.def("dbri_load_raw_gene_sets", [](const std::string& dbri_path) -> std::string {
        auto dbri = DBRetinaIndex::open(dbri_path);
        if (!dbri.has_section(DBRISection::RAW_GENE_SETS))
            return "";
        return dbri.load_raw_gene_sets();
    }, nb::arg("dbri_path"));

    m.def("dbri_load_names_list", [](const std::string& dbri_path) -> std::vector<std::string> {
        auto dbri = DBRetinaIndex::open(dbri_path);
        phmap::flat_hash_map<int, std::string> namesMap;
        dbri.load_names_map(namesMap);
        std::vector<std::string> names;
        names.reserve(namesMap.size());
        for (auto& [id, name] : namesMap) {
            names.push_back(name);
        }
        return names;
    }, nb::arg("dbri_path"));

    m.def("dbri_load_group_feature_counts", [](const std::string& dbri_path) -> std::map<std::string, int> {
        auto dbri = DBRetinaIndex::open(dbri_path);
        phmap::flat_hash_map<int, std::string> namesMap;
        dbri.load_names_map(namesMap);
        phmap::flat_hash_map<uint32_t, uint32_t> counts;
        dbri.load_group_feature_count(counts);
        std::map<std::string, int> result;
        for (auto& [gid, count] : counts) {
            auto it = namesMap.find(gid);
            if (it != namesMap.end()) {
                result[it->second] = count;
            }
        }
        return result;
    }, nb::arg("dbri_path"));

    // Global functions
    m.def("sketch_dbretina", &sketch_dbretina,
          nb::arg("asc_file"),
          nb::arg("names_file"),
          nb::arg("user_prefix") = std::string("NA"));

    m.def("query", &query,
          nb::arg("index_prefix"),
          nb::arg("inverted_index_prefix"),
          nb::arg("query_file"),
          nb::arg("output_prefix"),
          nb::arg("commands"));

    // .dbrp pairwise binary helpers
    m.def("dbrp_get_num_pairs", [](const std::string& dbrp_path) -> uint64_t {
        auto pw = DBRetinaPairwise::open(dbrp_path);
        return pw.get_num_pairs();
    }, nb::arg("dbrp_path"));

    m.def("dbrp_get_num_groups", [](const std::string& dbrp_path) -> uint32_t {
        auto pw = DBRetinaPairwise::open(dbrp_path);
        return pw.get_num_groups();
    }, nb::arg("dbrp_path"));

    m.def("dbrp_load_metadata", [](const std::string& dbrp_path) -> std::string {
        auto pw = DBRetinaPairwise::open(dbrp_path);
        return pw.get_metadata_json();
    }, nb::arg("dbrp_path"));

    m.def("dbrp_load_statistics", [](const std::string& dbrp_path) -> std::string {
        auto pw = DBRetinaPairwise::open(dbrp_path);
        return pw.get_statistics_json();
    }, nb::arg("dbrp_path"));

    m.def("dbrp_load_names", [](const std::string& dbrp_path) -> std::map<uint32_t, std::string> {
        auto pw = DBRetinaPairwise::open(dbrp_path);
        return pw.get_names_map();
    }, nb::arg("dbrp_path"));

    // Filter pairs by metric and cutoff, returns list of dicts
    // metric_id: 0=containment, 1=ochiai, 2=jaccard, 3=csi, 4=dice, 5=odds_ratio, 6=pvalue
    m.def("dbrp_filter_pairs", [](const std::string& dbrp_path, uint8_t metric_id, double cutoff) -> nb::list {
        auto pw = DBRetinaPairwise::open(dbrp_path);
        auto records = pw.filter_pairs(metric_id, cutoff);
        auto names = pw.get_names_map();

        nb::list result;
        for (auto& rec : records) {
            nb::dict d;
            d["group_1_id"] = rec.group_1_id;
            d["group_2_id"] = rec.group_2_id;
            std::string name1 = "unknown", name2 = "unknown";
            auto it1 = names.find(rec.group_1_id);
            if (it1 != names.end()) name1 = it1->second;
            auto it2 = names.find(rec.group_2_id);
            if (it2 != names.end()) name2 = it2->second;
            d["group_1_name"] = name1;
            d["group_2_name"] = name2;
            d["shared_features"] = rec.shared_features;
            d["containment"] = rec.containment;
            d["ochiai"] = rec.ochiai;
            d["jaccard"] = rec.jaccard;
            d["csi"] = rec.csi;
            d["dice"] = rec.dice;
            d["odds_ratio"] = rec.odds_ratio;
            if (pw.has_metric(DBRPMetric::PVALUE)) {
                d["pvalue"] = rec.pvalue;
            }
            result.append(d);
        }
        return result;
    }, nb::arg("dbrp_path"), nb::arg("metric_id"), nb::arg("cutoff"));

    // Query all pairs involving a specific group
    m.def("dbrp_query_group", [](const std::string& dbrp_path, uint32_t group_id, uint8_t metric_id, double cutoff) -> nb::list {
        auto pw = DBRetinaPairwise::open(dbrp_path);
        auto records = pw.query_group(group_id, metric_id, cutoff);
        auto names = pw.get_names_map();

        nb::list result;
        for (auto& rec : records) {
            nb::dict d;
            d["group_1_id"] = rec.group_1_id;
            d["group_2_id"] = rec.group_2_id;
            std::string name1 = "unknown", name2 = "unknown";
            auto it1 = names.find(rec.group_1_id);
            if (it1 != names.end()) name1 = it1->second;
            auto it2 = names.find(rec.group_2_id);
            if (it2 != names.end()) name2 = it2->second;
            d["group_1_name"] = name1;
            d["group_2_name"] = name2;
            d["shared_features"] = rec.shared_features;
            d["containment"] = rec.containment;
            d["ochiai"] = rec.ochiai;
            d["jaccard"] = rec.jaccard;
            d["csi"] = rec.csi;
            d["dice"] = rec.dice;
            d["odds_ratio"] = rec.odds_ratio;
            if (pw.has_metric(DBRPMetric::PVALUE)) {
                d["pvalue"] = rec.pvalue;
            }
            result.append(d);
        }
        return result;
    }, nb::arg("dbrp_path"), nb::arg("group_id"), nb::arg("metric_id"), nb::arg("cutoff"));

    // Iterate all pairs (no filter)
    m.def("dbrp_iterate_all", [](const std::string& dbrp_path) -> nb::list {
        auto pw = DBRetinaPairwise::open(dbrp_path);
        auto records = pw.iterate_all();
        auto names = pw.get_names_map();

        nb::list result;
        for (auto& rec : records) {
            nb::dict d;
            d["group_1_id"] = rec.group_1_id;
            d["group_2_id"] = rec.group_2_id;
            std::string name1 = "unknown", name2 = "unknown";
            auto it1 = names.find(rec.group_1_id);
            if (it1 != names.end()) name1 = it1->second;
            auto it2 = names.find(rec.group_2_id);
            if (it2 != names.end()) name2 = it2->second;
            d["group_1_name"] = name1;
            d["group_2_name"] = name2;
            d["shared_features"] = rec.shared_features;
            d["containment"] = rec.containment;
            d["ochiai"] = rec.ochiai;
            d["jaccard"] = rec.jaccard;
            d["csi"] = rec.csi;
            d["dice"] = rec.dice;
            d["odds_ratio"] = rec.odds_ratio;
            if (pw.has_metric(DBRPMetric::PVALUE)) {
                d["pvalue"] = rec.pvalue;
            }
            result.append(d);
        }
        return result;
    }, nb::arg("dbrp_path"));

    // GeneSets class
    nb::class_<GeneSets>(m, "GeneSets")
        .def(nb::init<>())
        .def(nb::init<std::string>(), nb::arg("associations_file"))
        .def("build_from_index", &GeneSets::build_from_index,
             nb::arg("index_prefix"))
        .def("build_from_clusters_file", &GeneSets::build_from_clusters_file,
             nb::arg("clusters_file"))
        .def("get_pathways_ppi", &GeneSets::get_pathways_ppi)
        .def("get_pathways_psi", &GeneSets::get_pathways_psi)
        .def("get_pathways_pcsi", &GeneSets::get_pathways_pcsi)
        .def("get_pathway_lengths", &GeneSets::get_pathway_lengths)
        .def("export_genes_to_ppi_psi_tsv", &GeneSets::export_genes_to_ppi_psi_tsv,
             nb::arg("filename"))
        .def("calculate_heterogeneity_and_fragmentation_from_pairwise",
             &GeneSets::calculate_heterogeneity_and_fragmentation_from_pairwise,
             nb::arg("pairwise_file"))
        .def("get_pathway_to_modularity", &GeneSets::get_pathway_to_modularity)
        .def("get_pathway_to_heterogeneity", &GeneSets::get_pathway_to_heterogeneity)
        .def("get_pathway_to_fragmentation", &GeneSets::get_pathway_to_fragmentation)
        .def("non_iterative_set_cover", &GeneSets::non_iterative_set_cover,
             nb::arg("cluster_id"), nb::arg("GC"))
        .def("keep_only_these_pathways", &GeneSets::keep_only_these_pathways,
             nb::arg("non_redundant_pathways"));
}
