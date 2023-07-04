%module kSpider_internal

%{
#include "DBRetina_GSA.hpp"
#include "kSpider.hpp"
#include "DBRetina.hpp"
%}


using namespace std;
%include std_string.i

namespace kSpider{
    void pairwise(string index_prefix, int user_threads, string cutoff_distance_type, double cutoff_threshold, string full_command, bool calculate_pvalue);
    // void index_kmers(int kSize, string fasta_file, string names_file, int chunk_size, string index_prefix);
    // void index_kmers_nonCanonical(int kSize, string fasta_file, string names_file, int chunk_size, string index_prefix);
    // void index_skipmers(int m, int n, int k, string fasta_file, string names_file, int chunk_size, string index_prefix);
    // void index_protein(int kSize, string fasta_file, string names_file, int chunk_size, string index_prefix);
    // void index_dayhoff(int kSize, string fasta_file, string names_file, int chunk_size, string index_prefix);
    // void index_datasets(string kfs_dir);
    // void sourmash_sigs_indexing(string sigs_dir, int kSize);
    // void paired_end_to_kDataFrame(string r1_file_name, string r2_file_name, int kSize, int chunk_size, int downsampling_ratio, bool remove_singletones);
    // void single_end_to_kDataFrame(string r1_file_name, int kSize, int chunk_size, int downsampling_ration, bool remove_singletones);
    // void protein_to_kDataFrame(string r1_file_name, int kSize, int chunk_size, bool is_dayhoff, string output_prefix, int downsampling_ration = 1);
    void dbretina_indexing(string json_file, string user_index_prefix);
};

void sketch_dbretina(string asc_file, string names_file, string user_prefix);
void parse_dbretina_json(string json_file, str_hashed_set_map* map);
void query(string index_prefix, string inverted_index_prefix, string query_file, string output_prefix, string commands);



%include std_unordered_map.i
%include <stdint.i>


%template(StringDoubleMap) unordered_map<string, double>;
%typemap(out) std::unordered_map<std::string, double> {
    PyObject* obj = PyDict_New();
    if (obj == NULL) {
        SWIG_fail;
    }

    for (auto it = $1.begin(); it != $1.end(); ++it) {
        PyObject* key = SWIG_From_std_string(it->first);
        PyObject* val = PyFloat_FromDouble(it->second);

        if (PyDict_SetItem(obj, key, val) == -1) {
            Py_DECREF(obj);
            SWIG_fail;
        }

        Py_DECREF(key);
        Py_DECREF(val);
    }

    $result = obj;
}

%template(StringBoolMap) unordered_map<string, bool>;
%typemap(out) std::unordered_map<std::string, bool> {
    PyObject* o = PyDict_New();
    if (o == NULL) {
        return NULL;
    }
    
    for (auto it = $1.begin(); it != $1.end(); ++it) {
        PyObject* k = PyUnicode_FromString(it->first.c_str());
        PyObject* v = PyBool_FromLong(it->second);
        if (k == NULL || v == NULL || PyDict_SetItem(o, k, v) == -1) {
            Py_XDECREF(k);
            Py_XDECREF(v);
            Py_DECREF(o);
            return NULL;
        }
        Py_DECREF(k);
        Py_DECREF(v);
    }
    $result = o;
}


%template(StringintMap) unordered_map<string, int>;
%typemap(out) std::unordered_map<std::string, int> {
    PyObject* o = PyDict_New();
    if (o == NULL) {
        return NULL;
    }
    
    for (auto it = $1.begin(); it != $1.end(); ++it) {
        PyObject* k = PyUnicode_FromString(it->first.c_str());
        PyObject* v = PyLong_FromLong(it->second);
        if (k == NULL || v == NULL || PyDict_SetItem(o, k, v) == -1) {
            Py_XDECREF(k);
            Py_XDECREF(v);
            Py_DECREF(o);
            return NULL;
        }
        Py_DECREF(k);
        Py_DECREF(v);
    }
    $result = o;
}


class GeneSets {
public:
    GeneSets();
    GeneSets(string associations_file);
    void build_from_index(string index_prefix);
    void build_from_clusters_file(string clusters_file);
    unordered_map<string, double> get_pathways_ppi();
    unordered_map<string, double> get_pathways_psi();
    unordered_map<string, double> get_pathways_pcsi();
    unordered_map<string, int> get_pathway_lengths();
    void export_genes_to_ppi_psi_tsv(string filename);
    void calculate_heterogeneity_and_fragmentation_from_pairwise(string pairwise_file);
    unordered_map<string, int> get_pathway_to_modularity();
    unordered_map<string, int> get_pathway_to_heterogeneity();
    unordered_map<string, int> get_pathway_to_fragmentation();
    unordered_map<string, double> non_iterative_set_cover(int cluster_id, int GC);
    void keep_only_these_pathways(string non_redundant_pathways);
    ~GeneSets();
};

