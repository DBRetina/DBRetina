#include "dbretina_core.hpp"
#include <iostream>
#include <cstdint>
#include "DBRetina_kf.hpp"
#include "parallel_hashmap/phmap.h"
#include <vector>
#include <stdexcept>
#include <sstream>
#include <string>
#include <fstream>
#include "parallel_hashmap/phmap_dump.h"
#include "DBRetina.hpp"
#include "DBRetinaIndex.hpp"

using namespace phmap;
using namespace std;

namespace dbretina {

    void dbretina_merge(string index_a_path, string index_b_path, string output_dbri_path) {

        // 1. Load index A fully (base index)
        auto dbriA = DBRetinaIndex::open(index_a_path);

        cerr << "[merge] Loading base index from " << index_a_path << endl;

        DBRetina_PHMAP* frame = dbriA.load_phmap();

        flat_hash_map<int, std::string> namesMap_int;
        dbriA.load_names_map(namesMap_int);

        flat_hash_map<string, string> namesMap;
        flat_hash_map<string, uint64_t> groupNameMap;
        for (auto& [gid, name] : namesMap_int) {
            namesMap[name] = name;
            groupNameMap[name] = gid;
        }

        flat_hash_map<string, uint64_t> tagsMap;
        dbriA.load_tags_map(tagsMap);

        priority_queue<uint64_t, vector<uint64_t>, std::greater<uint64_t>> freeColors;
        dbriA.load_free_colors(freeColors);

        flat_hash_map<uint32_t, uint32_t> groupID_to_featureCount;
        dbriA.load_group_feature_count(groupID_to_featureCount);

        flat_hash_map<string, uint32_t> groupName_to_featureCount;
        for (auto& [gid, count] : groupID_to_featureCount) {
            auto it = namesMap_int.find(gid);
            if (it != namesMap_int.end()) {
                groupName_to_featureCount[it->second] = count;
            }
        }

        phmap::parallel_flat_hash_map<uint32_t, std::vector<uint32_t>,
            std::hash<uint32_t>, std::equal_to<uint32_t>,
            std::allocator<std::pair<const uint32_t, std::vector<uint32_t>>>, 1> color_to_ids_par;
        dbriA.load_color_to_sources(color_to_ids_par);

        auto* legend = new flat_hash_map<uint64_t, std::vector<uint32_t>>();
        for (auto& [color, ids] : color_to_ids_par) {
            (*legend)[color] = ids;
        }

        phmap::parallel_flat_hash_map<uint32_t, uint32_t,
            std::hash<uint32_t>, std::equal_to<uint32_t>,
            std::allocator<std::pair<const uint32_t, uint32_t>>, 1> colorsCount_par;
        dbriA.load_color_count(colorsCount_par);

        flat_hash_map<uint64_t, uint64_t> colorsCount;
        for (auto& [k, v] : colorsCount_par) {
            colorsCount[k] = v;
        }

        // Find max ID across all existing IDs
        uint64_t max_group_id = dbriA.get_max_group_id();
        uint64_t max_color_id = max_group_id;
        for (auto& [color, _] : color_to_ids_par) {
            if (color > max_color_id) max_color_id = color;
        }
        for (auto& [_, id] : tagsMap) {
            if (id > max_color_id) max_color_id = id;
        }
        uint64_t groupID = max_color_id + 1;

        int total_a_groups = namesMap_int.size();
        cerr << "[merge] Base index: " << total_a_groups << " groups" << endl;

        // 2. Load index B's gene sets
        auto dbriB = DBRetinaIndex::open(index_b_path);

        if (!dbriB.has_section(DBRISection::GROUP_TO_FEATURE_SET)) {
            throw std::runtime_error("Index B does not contain GROUP_TO_FEATURE_SET section needed for merge.");
        }

        phmap::parallel_flat_hash_map<std::string, phmap::flat_hash_set<uint64_t>> bFeatureSets;
        dbriB.load_group_to_feature_set(bFeatureSets);

        int total_b_groups = bFeatureSets.size();
        cerr << "[merge] Merging " << total_b_groups << " groups from " << index_b_path << endl;

        // 3. Register B's groups with new IDs
        flat_hash_map<string, uint64_t> groupCounter;
        int new_groups_count = 0;
        for (auto& [group_name, kset] : bFeatureSets) {
            if (groupNameMap.find(group_name) != groupNameMap.end()) {
                throw std::runtime_error("Duplicate gene set name '" + group_name +
                    "' exists in both indexes. Use --prefix to disambiguate.");
            }

            new_groups_count++;
            namesMap[group_name] = group_name;
            groupNameMap[group_name] = groupID;
            groupCounter[group_name] = 1;
            tagsMap[to_string(groupID)] = groupID;
            vector<uint32_t> tmp;
            tmp.push_back(groupID);
            legend->insert(make_pair(groupID, tmp));
            colorsCount.insert(make_pair(groupID, 0));
            groupID++;
        }

        // 4. Build inverse group name map
        flat_hash_map<uint64_t, string> inv_groupNameMap;
        for (auto& [name, gid] : groupNameMap) {
            inv_groupNameMap[gid] = name;
        }

        // 5. Run color update loop for B's gene sets (same as append)
        for (auto& [group_name, featureSet] : bFeatureSets) {
            for (auto it = featureSet.begin(); it != featureSet.end(); ++it) {
                flat_hash_map<uint64_t, uint64_t> convertMap;
                uint64_t groupTag = groupNameMap.find(group_name)->second;

                convertMap.clear();
                convertMap.insert(make_pair(0, groupTag));
                convertMap.insert(make_pair(groupTag, groupTag));

                groupName_to_featureCount[group_name] = featureSet.size();

                auto loaded_features_it = featureSet.begin();
                while (loaded_features_it != featureSet.end()) {
                    uint64_t hashed_feature = *loaded_features_it;
                    uint64_t currentTag = frame->getCount(hashed_feature);
                    auto itc = convertMap.find(currentTag);
                    if (itc == convertMap.end()) {
                        auto legendFind = legend->find(currentTag);
                        if (legendFind == legend->end()) {
                            loaded_features_it++;
                            continue;
                        }
                        vector<uint32_t> colors = legendFind->second;
                        auto tmpiT = find(colors.begin(), colors.end(), groupTag);
                        if (tmpiT == colors.end()) {
                            colors.push_back(groupTag);
                            sort(colors.begin(), colors.end());
                        }

                        std::stringstream ss;
                        if (!colors.empty()) {
                            ss << colors[0];
                            for (size_t k = 1; k < colors.size(); ++k) {
                                ss << ";" << colors[k];
                            }
                        }
                        std::string colorsString = ss.str();

                        auto itTag = tagsMap.find(colorsString);

                        if (itTag == tagsMap.end()) {
                            uint64_t newColor;
                            if (freeColors.empty()) newColor = groupID++;
                            else {
                                newColor = freeColors.top();
                                freeColors.pop();
                            }
                            auto inserted = tagsMap.emplace(colorsString, newColor);
                            (*legend)[newColor] = colors;
                            itTag = inserted.first;
                            colorsCount[newColor] = 0;
                        }
                        uint64_t newColor = itTag->second;

                        convertMap.emplace(currentTag, newColor);
                        itc = convertMap.find(currentTag);
                    }

                    if (itc->second != currentTag) {
                        --colorsCount[currentTag];
                        if (colorsCount[currentTag] == 0 && currentTag != 0) {
                            auto _invGroupNameIT = inv_groupNameMap.find(currentTag);
                            if (_invGroupNameIT == inv_groupNameMap.end()) {
                                freeColors.push(currentTag);
                                auto legendIt = legend->find(currentTag);
                                if (legendIt != legend->end()) {
                                    std::vector<uint32_t> colors = legendIt->second;
                                    std::stringstream ss;
                                    ss << colors[0];
                                    for (unsigned int k = 1; k < colors.size(); k++) {
                                        ss << ";" << colors[k];
                                    }
                                    std::string colorsString = ss.str();
                                    tagsMap.erase(colorsString);
                                }
                                auto convertMapIt = convertMap.find(currentTag);
                                if (convertMapIt != convertMap.end()) {
                                    convertMap.erase(convertMapIt);
                                }
                            }
                        }
                        ++colorsCount[itc->second];
                    }

                    frame->setCount(hashed_feature, itc->second);
                    loaded_features_it++;
                }

                groupCounter[group_name]--;
                if (colorsCount[groupTag] == 0) {
                    if (groupCounter[group_name] == 0) {
                        freeColors.push(groupTag);
                    }
                }
                break;
            }
        }

        cerr << "[merge] Total colors after merge: " << legend->size() << endl;

        // 6. Build updated groupID_to_featureCount
        flat_hash_map<uint32_t, uint32_t> updated_groupID_to_featureCount;
        for (auto& [gname, count] : groupName_to_featureCount) {
            updated_groupID_to_featureCount[groupNameMap[gname]] = count;
        }

        // Find new max group ID
        uint64_t new_max_group_id = 0;
        for (auto& [name, gid] : groupNameMap) {
            if (gid > new_max_group_id) new_max_group_id = gid;
        }

        int total_groups = namesMap.size();

        std::string metadata = "{\"population_size\":" + to_string(frame->size())
            + ",\"num_groups\":" + to_string(total_groups)
            + ",\"max_group_id\":" + to_string(new_max_group_id)
            + ",\"hash_size\":31"
            + "}";

        // 7. Write merged .dbri
        DBRetinaIndex out;
        out.begin_write(output_dbri_path);

        out.write_phmap(frame);
        out.write_color_to_sources(*legend, colorsCount);
        out.write_color_count(colorsCount);
        out.write_group_feature_count(updated_groupID_to_featureCount);
        out.write_names_map(namesMap, groupNameMap);
        out.write_metadata(metadata);
        out.write_tags_map(tagsMap);
        out.write_free_colors(freeColors);

        // Merge feature sets from both indexes
        phmap::parallel_flat_hash_map<std::string, phmap::flat_hash_set<uint64_t>> allFeatureSets;
        if (dbriA.has_section(DBRISection::GROUP_TO_FEATURE_SET)) {
            dbriA.load_group_to_feature_set(allFeatureSets);
        }
        for (auto& [name, kset] : bFeatureSets) {
            allFeatureSets[name] = kset;
        }
        out.write_group_to_feature_set(allFeatureSets);

        // Carry forward embedded JSON from A (B's data would need separate handling)
        if (dbriA.has_section(DBRISection::RAW_GENE_SETS)) {
            std::string existing_raw = dbriA.load_raw_gene_sets();
            out.write_raw_gene_sets(existing_raw);
        }
        if (dbriA.has_section(DBRISection::HASHED_GENE_SETS)) {
            std::string existing_hashed = dbriA.load_hashed_gene_sets();
            out.write_hashed_gene_sets(existing_hashed);
        }

        out.finalize_write();

        cerr << "[merge] Wrote merged index to " << output_dbri_path << endl;
        cerr << "[merge] Total groups: " << total_groups
             << " (" << total_a_groups << " from A + " << new_groups_count << " from B)" << endl;

        delete legend;
        delete frame;
    }

}
