#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/unordered_map.h>

#include "kSpider.hpp"
#include "DBRetina.hpp"
#include "DBRetina_GSA.hpp"

namespace nb = nanobind;

NB_MODULE(_kSpider_internal, m) {
    m.doc() = "DBRetina C++ bindings";

    // kSpider namespace functions
    m.def("pairwise", &kSpider::pairwise,
          nb::arg("index_prefix"),
          nb::arg("user_threads"),
          nb::arg("cutoff_distance_type"),
          nb::arg("cutoff_threshold"),
          nb::arg("full_command"),
          nb::arg("calculate_pvalue"));

    m.def("dbretina_indexing", &kSpider::dbretina_indexing,
          nb::arg("json_file"),
          nb::arg("user_index_prefix"));

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
