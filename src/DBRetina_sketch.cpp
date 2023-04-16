#include <iostream>
#include <string>
#include "parallel_hashmap/phmap.h"
#include <sstream>
#include <fstream>
#include "parallel_hashmap/phmap_dump.h"
#include <sstream>
#include "DBRetina.hpp"

inline uint64_t to_uint64(std::string const& value) {
    uint64_t result = 0;

    size_t const length = value.size();
    switch (length) {
    case 20:    result += (value[length - 20] - '0') * 10000000000000000000ULL;
    case 19:    result += (value[length - 19] - '0') * 1000000000000000000ULL;
    case 18:    result += (value[length - 18] - '0') * 100000000000000000ULL;
    case 17:    result += (value[length - 17] - '0') * 10000000000000000ULL;
    case 16:    result += (value[length - 16] - '0') * 1000000000000000ULL;
    case 15:    result += (value[length - 15] - '0') * 100000000000000ULL;
    case 14:    result += (value[length - 14] - '0') * 10000000000000ULL;
    case 13:    result += (value[length - 13] - '0') * 1000000000000ULL;
    case 12:    result += (value[length - 12] - '0') * 100000000000ULL;
    case 11:    result += (value[length - 11] - '0') * 10000000000ULL;
    case 10:    result += (value[length - 10] - '0') * 1000000000ULL;
    case  9:    result += (value[length - 9] - '0') * 100000000ULL;
    case  8:    result += (value[length - 8] - '0') * 10000000ULL;
    case  7:    result += (value[length - 7] - '0') * 1000000ULL;
    case  6:    result += (value[length - 6] - '0') * 100000ULL;
    case  5:    result += (value[length - 5] - '0') * 10000ULL;
    case  4:    result += (value[length - 4] - '0') * 1000ULL;
    case  3:    result += (value[length - 3] - '0') * 100ULL;
    case  2:    result += (value[length - 2] - '0') * 10ULL;
    case  1:    result += (value[length - 1] - '0');
    }
    return result;
}

struct string_hasher
{
    std::hash<string> hasher;

    // initialize the hasher
    string_hasher(): hasher() {}

    // overload the () operator
    size_t operator()(const string& str) const
    {
        return hasher(str);
    }
};

void load_tsv_to_map(string filename, str_vec_map* map, str_str_map* names_map) {

    std::ifstream inputFile(filename);

    if (!inputFile.is_open()) {
        std::cerr << "Error opening the file: " << filename << std::endl;
        return;
    }

    std::string line;
    std::getline(inputFile, line); // skip first line
    while (std::getline(inputFile, line)) {
        std::istringstream lineStream(line);
        std::string column1, column2;

        if (std::getline(lineStream, column1, '\t') && std::getline(lineStream, column2, '\t')) {
            transform(column1.begin(), column1.end(), column1.begin(), ::tolower);
            transform(column2.begin(), column2.end(), column2.begin(), ::tolower);
            map->operator[](names_map->operator[](column1)).emplace(column2);
        }
        else {
            std::cerr << "Invalid line format: " << line << std::endl;
            inputFile.close();
            return;
        }
    }

    inputFile.close();
}

void load_tsv_to_map_no_names(string filename, str_vec_map* map) {

    std::ifstream inputFile(filename);

    if (!inputFile.is_open()) {
        std::cerr << "Error opening the file: " << filename << std::endl;
        return;
    }

    std::string line;
    std::getline(inputFile, line); // skip first line
    while (std::getline(inputFile, line)) {
        std::istringstream lineStream(line);
        std::string column1, column2;

        if (std::getline(lineStream, column1, '\t') && std::getline(lineStream, column2, '\t')) {
            transform(column1.begin(), column1.end(), column1.begin(), ::tolower);
            transform(column2.begin(), column2.end(), column2.begin(), ::tolower);
            map->operator[](column1).emplace(column2);
        }
        else {
            std::cerr << "Invalid line format: " << line << std::endl;
            inputFile.close();
            return;
        }
    }

    inputFile.close();
}


void load_names_tsv_to_map(string filename, str_str_map* map) {
    std::ifstream inputFile(filename);

    if (!inputFile.is_open()) {
        std::cerr << "Error opening the file: " << filename << std::endl;
        return;
    }

    std::string line;
    std::getline(inputFile, line); // skip first line
    while (std::getline(inputFile, line)) {
        std::istringstream lineStream(line);
        std::string column1, column2;

        if (std::getline(lineStream, column1, '\t') && std::getline(lineStream, column2, '\t')) {
            transform(column1.begin(), column1.end(), column1.begin(), ::tolower);
            transform(column2.begin(), column2.end(), column2.begin(), ::tolower);
            map->operator[](column1) = column2;
        }
        else {
            std::cerr << "Invalid line format: " << line << std::endl;
            inputFile.close();
            return;
        }
    }

    inputFile.close();
}


void sketch_dbretina(string asc_file, string names_file) {
    str_vec_map* asc_map = new str_vec_map();
    str_str_map* names_map = new str_str_map();

    string asc_file_wo_extension = asc_file.substr(0, asc_file.find_last_of("."));
    string asc_file_base_name_without_extension = asc_file_wo_extension.substr(asc_file_wo_extension.find_last_of("/") + 1);
    string private_output_file = asc_file_base_name_without_extension + "_private.json";
    string public_output_file = asc_file_base_name_without_extension + "_public.json";

    if (names_file == "NA") {
        load_tsv_to_map_no_names(asc_file, asc_map);
        cout << "number of loaded groups: " << names_map->size() << endl;
    }
    else {
        load_names_tsv_to_map(names_file, names_map);
        load_tsv_to_map(asc_file, asc_map, names_map);
        cout << "number of loaded groups: " << names_map->size() << endl;
        cout << "number of loaded asc:" << asc_map->size() << endl;
    }

    /* DISABLED NAMESFILE EXPORT
    parallel_flat_hash_map<string, parallel_flat_hash_set<string>> transformed_names;
    for (auto it = names_map->begin(); it != names_map->end(); ++it)
    {
        transformed_names[it->second].insert(it->first);
    }

    // export transformed_names to json
    ofstream file(asc_file_base_name_without_extension + "_names.json");
    file << "{";
    for (auto it = transformed_names.begin(); it != transformed_names.end(); ++it)
    {
        string parent_group = it->first;
        file << "\"" << parent_group << "\": [";
        for (auto it2 = it->second.begin(); it2 != it->second.end(); ++it2)
        {
            file << "\"" << *it2 << "\"";
            if (next(it2) != it->second.end())
            {
                file << ", ";
            }
        }
        file << "]";
        if (next(it) != asc_map->end())
        {
            file << ", ";
        }
    }
    file << "}";
    file.close();
*/


    auto hasher = string_hasher();

    // Public version
    // file.open(private_output_file); .. related to the previous disabling
    ofstream file(private_output_file);
    file << "{";
    file << "\"" << "metadata" << "\": {"
        << "\"" << "filetype" << "\": \"" << "private" << "\""
        << "},";
    file << "\"" << "data" << "\": {";
    for (auto it = asc_map->begin(); it != asc_map->end(); ++it)
    {
        string parent_group = it->first;
        file << "\"" << parent_group << "\": [";
        for (auto it2 = it->second.begin(); it2 != it->second.end(); ++it2)
        {
            file << "\"" << *it2 << "\"";
            if (next(it2) != it->second.end())
            {
                file << ", ";
            }
        }
        file << "]";
        if (next(it) != asc_map->end())
        {
            file << ", ";
        }
    }
    file << "}}";
    file.close();


    // Intermediate step to store reverse hashing
    // parallel_flat_hash_map<uint64_t, string> reverse_hash_map; // disabled reverse hashing


    // dumping json to be indexed
    // redeclare file
    file.open(public_output_file);

    file << "{";
    file << "\"" << "metadata" << "\": {"
        << "\"" << "filetype" << "\": \"" << "public" << "\""
        << "},";
    file << "\"" << "data" << "\": {";
    for (auto it = asc_map->begin(); it != asc_map->end(); ++it)
    {
        string parent_group = it->first;
        file << "\"" << parent_group << "\": [";
        for (auto it2 = it->second.begin(); it2 != it->second.end(); ++it2)
        {
            uint64_t hashed_child = hasher(*it2);
            // reverse_hash_map[hashed_child] = *it2; // disabled reverse hashing
            file << "\"" << hashed_child << "\"";
            if (next(it2) != it->second.end())
            {
                file << ", ";
            }
        }
        file << "]";
        if (next(it) != asc_map->end())
        {
            file << ", ";
        }
    }
    file << "}}";
    file.close();

    // dump reverse hash map to TSV // disabled reverse hashing
    // ofstream reverse_hash_map_file(asc_file_base_name_without_extension + "_reverse.tsv");
    // for (auto it = reverse_hash_map.begin(); it != reverse_hash_map.end(); ++it)
    // {
    //     reverse_hash_map_file << it->first << "\t" << it->second << endl;
    // }

}




void parse_dbretina_json(string json_file, str_hashed_vec_map* map) {
    string_hasher hasher = string_hasher();
    string output_prefix = json_file.substr(0, json_file.find_last_of(".")).substr(json_file.find_last_of("/") + 1);

    zstr::ifstream sig_stream(json_file);
    json::value json = json::parse(sig_stream);
    string filetype = string(json["metadata"]["filetype"].as_string());
    auto genes = json["data"].as_object();

    if (filetype == "private") {
        // Intermediate step to store reverse hashing
        parallel_flat_hash_map<uint64_t, string> reverse_hash_map;

        // iterate over genes:
        for (auto it = genes.begin(); it != genes.end(); ++it) {
            string parent_name = it->first;
            auto gene = it->second.as_array();
            for (auto it2 = gene.begin(); it2 != gene.end(); ++it2) {
                string str_child = it2->as_string();
                uint64_t hashed_child = hasher(str_child);
                reverse_hash_map[hashed_child] = str_child;
                map->operator[](parent_name).emplace(hashed_child);
            }
        }

    }
    else if (filetype == "public") {
        // iterate over genes:
        for (auto it = genes.begin(); it != genes.end(); it++) {
            string parent_name = it->first;
            cout << "parent_name: " << parent_name << endl;
            auto gene = it->second.as_array();
            map->operator[](parent_name) = parallel_flat_hash_set<uint64_t>(gene.size());
            for (auto it2 = gene.begin(); it2 != gene.end(); it2++) {
                map->operator[](parent_name).insert(to_uint64(it2->as_string()));
            }
        }
    }
    else {
        cout << "Unknown filetype: " << filetype << endl;
        exit(0);
    }
}


// void DBRetina_json_parser::load_json(string json_file) {
//     zstr::ifstream sig_stream(json_file);
//     json::value json = json::parse(sig_stream);
//     string filetype = string(json["metadata"]["filetype"].as_string());
//     auto genes = json["data"].as_object();

//     if (filetype == "private") {
//         // Intermediate step to store reverse hashing
//         parallel_flat_hash_map<uint64_t, string> reverse_hash_map;

//         // iterate over genes:
//         for (auto it = genes.begin(); it != genes.end(); ++it) {
//             string parent_name = it->first;
//             auto gene = it->second.as_array();
//             for (auto it2 = gene.begin(); it2 != gene.end(); ++it2) {
//                 string str_child = it2->as_string();
//                 uint64_t hashed_child = hasher(str_child);
//                 reverse_hash_map[hashed_child] = str_child;
//                 map->operator[](parent_name).emplace(hashed_child);
//             }
//         }

//     }
//     else if (filetype == "public") {
//         // iterate over genes:
//         for (auto it = genes.begin(); it != genes.end(); it++) {
//             string parent_name = it->first;
//             cout << "parent_name: " << parent_name << endl;
//             auto gene = it->second.as_array();
//             map->operator[](parent_name) = parallel_flat_hash_set<uint64_t>(gene.size());
//             for (auto it2 = gene.begin(); it2 != gene.end(); it2++) {
//                 map->operator[](parent_name).insert(to_uint64(it2->as_string()));
//             }
//         }
//     }
//     else {
//         cout << "Unknown filetype: " << filetype << endl;
//         exit(0);
//     }

// }

// void DBRetina_json_parser::export_map_to_tsv(string output_file) {
//     ofstream file(output_file);
//     for (auto it = map->begin(); it != map->end(); ++it) {
//         file << it->first << "\t";
//         for (auto it2 = it->second.begin(); it2 != it->second.end(); ++it2) {
//             file << *it2;
//             if (next(it2) != it->second.end()) {
//                 file << ", ";
//             }
//         }
//         file << endl;
//     }
//     file.close();
// }

// DBRetina_json_parser::DBRetina_json_parser(string json_file) {
//     json_file_name = json_file;
//     sig_stream.open(json_file);
//     json = json::parse(sig_stream);
//     hasher = string_hasher();
//     output_prefix = json_file.substr(0, json_file.find_last_of(".")).substr(json_file.find_last_of("/") + 1);
//     load_json(json_file);
//     cout << "exporting map" << endl;
//     export_map_to_tsv(output_prefix + "_map.tsv");
// }




// int main(int argc, char* argv[]) {

//     string asc_file = argv[1];
//     string names_file = argv[2];
//     string asc_file_wo_extension = asc_file.substr(0, asc_file.find_last_of("."));
//     string asc_file_base_name_without_extension = asc_file_wo_extension.substr(asc_file_wo_extension.find_last_of("/") + 1);
//     string private_output_file = asc_file_base_name_without_extension + "_private.json";
//     string public_output_file = asc_file_base_name_without_extension + "_public.json";

//     sketch_dbretine(asc_file, names_file);

//     DBRetina_json_parser parser(private_output_file);

// }