%module kSpider_internal

%{
#include "DBRetina_GSA.hpp"
#include "kSpider.hpp"
#include "DBRetina.hpp"
%}


using namespace std;
%include std_string.i

namespace kSpider{
    void pairwise(string index_prefix, int user_threads, string cutoff_distance_type, double cutoff_threshold, string full_command);
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
void parse_dbretina_json(string json_file, str_hashed_vec_map* map);
void query(string index_prefix, string inverted_index_prefix, string query_file, string output_prefix, string commands);



%include std_unordered_map.i

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


class GeneSets {
    // empty constructor
public:
    GeneSets();
    GeneSets(string associations_file);
    void build_from_index(string index_prefix);
    void build_from_clusters_file(string clusters_file);
    unordered_map<string, double> get_pathways_ppi();
    unordered_map<string, double> get_pathways_psi();
    ~GeneSets();
};