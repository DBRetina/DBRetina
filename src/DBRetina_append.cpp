#include "dbretina_core.hpp"
#include <iostream>
#include <cstdint>
#include <chrono>
#include "DBRetina_kf.hpp"
#include "parallel_hashmap/phmap.h"
#include <vector>
#include <stdexcept>
#include <sstream>
#include <string>
#include <fstream>
#include "cpp-json/json.h"
#include "zstr.hpp"
#include "parallel_hashmap/phmap_dump.h"
#include "DBRetina.hpp"
#include "DBRetinaIndex.hpp"

using namespace phmap;
using namespace std;

namespace dbretina {

    void dbretina_append(string existing_dbri_path, string new_json_file, string output_dbri_path) {

        // 1. Load existing index
        auto dbri = DBRetinaIndex::open(existing_dbri_path);

        cerr << "[append] Loading existing index from " << existing_dbri_path << endl;

        // Load PHMAP
        DBRetina_PHMAP* frame = dbri.load_phmap();

        // Load namesMap
        flat_hash_map<int, std::string> namesMap_int;
        dbri.load_names_map(namesMap_int);

        // Reconstruct namesMap and groupNameMap
        flat_hash_map<string, string> namesMap;
        flat_hash_map<string, uint64_t> groupNameMap;
        for (auto& [gid, name] : namesMap_int) {
            namesMap[name] = name;
            groupNameMap[name] = gid;
        }

        // Load tagsMap
        flat_hash_map<string, uint64_t> tagsMap;
        dbri.load_tags_map(tagsMap);

        // Load freeColors
        priority_queue<uint64_t, vector<uint64_t>, std::greater<uint64_t>> freeColors;
        dbri.load_free_colors(freeColors);

        // Load feature counts
        flat_hash_map<uint32_t, uint32_t> groupID_to_featureCount;
        dbri.load_group_feature_count(groupID_to_featureCount);

        // Reconstruct groupName_to_featureCount
        flat_hash_map<string, uint32_t> groupName_to_featureCount;
        for (auto& [gid, count] : groupID_to_featureCount) {
            auto it = namesMap_int.find(gid);
            if (it != namesMap_int.end()) {
                groupName_to_featureCount[it->second] = count;
            }
        }

        // Load color_to_sources as legend (now persisted with ALL entries)
        phmap::parallel_flat_hash_map<uint32_t, std::vector<uint32_t>,
            std::hash<uint32_t>, std::equal_to<uint32_t>,
            std::allocator<std::pair<const uint32_t, std::vector<uint32_t>>>, 1> color_to_ids_par;
        dbri.load_color_to_sources(color_to_ids_par);

        auto* legend = new flat_hash_map<uint64_t, std::vector<uint32_t>>();
        for (auto& [color, ids] : color_to_ids_par) {
            (*legend)[color] = ids;
        }

        // Load colorsCount
        phmap::parallel_flat_hash_map<uint32_t, uint32_t,
            std::hash<uint32_t>, std::equal_to<uint32_t>,
            std::allocator<std::pair<const uint32_t, uint32_t>>, 1> colorsCount_par;
        dbri.load_color_count(colorsCount_par);

        flat_hash_map<uint64_t, uint64_t> colorsCount;
        for (auto& [k, v] : colorsCount_par) {
            colorsCount[k] = v;
        }

        // Get max group ID and find max color ID in legend
        uint64_t max_group_id = dbri.get_max_group_id();
        uint64_t max_color_id = max_group_id;
        for (auto& [color, _] : color_to_ids_par) {
            if (color > max_color_id) max_color_id = color;
        }
        // Also check tagsMap for any higher IDs
        for (auto& [_, id] : tagsMap) {
            if (id > max_color_id) max_color_id = id;
        }
        uint64_t groupID = max_color_id + 1;

        int total_existing_groups = namesMap_int.size();

        cerr << "[append] Loaded " << total_existing_groups << " existing groups, max_group_id=" << max_group_id << endl;

        // 2. Parse new gene sets
        hashed_MAP newGroupName_to_featureSet;
        parse_dbretina_json(new_json_file, &newGroupName_to_featureSet);

        int new_groups_count = 0;

        // 3. Register new groups
        flat_hash_map<string, uint64_t> groupCounter;
        for (auto& group : newGroupName_to_featureSet) {
            string group_name = group.first;

            // Check for duplicate names
            if (groupNameMap.find(group_name) != groupNameMap.end()) {
                throw std::runtime_error("Duplicate gene set name '" + group_name +
                    "' already exists in index. Use --force to replace.");
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

        cerr << "[append] Adding " << new_groups_count << " new gene sets" << endl;

        // 4. Build inverse group name map
        flat_hash_map<uint64_t, string> inv_groupNameMap;
        for (auto& [name, gid] : groupNameMap) {
            inv_groupNameMap[gid] = name;
        }

        // 5. Run color update loop for new gene sets (same algorithm as indexing)
        for (auto& group : newGroupName_to_featureSet) {
            string group_name = group.first;
            flat_hash_set<uint64_t>& featureSet = group.second;

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
                            // Feature has a stale color not in legend — treat as uncolored
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
                                // Don't erase from legend — keep entries for subsequent groups
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
                        // Don't erase from legend
                    }
                }
                break;
            }
        }

        cerr << "[append] Total colors after append: " << legend->size() << endl;

        // 6. Build updated groupID_to_featureCount
        flat_hash_map<uint32_t, uint32_t> updated_groupID_to_featureCount;
        for (auto& [gname, count] : groupName_to_featureCount) {
            updated_groupID_to_featureCount[groupNameMap[gname]] = count;
        }

        // Compute selected colors
        size_t selected_colors = 0;
        for (auto& [color, sources] : *legend) {
            if (colorsCount[color] > 0) selected_colors++;
        }

        // Find new max group ID
        uint64_t new_max_group_id = 0;
        for (auto& [name, gid] : groupNameMap) {
            if (gid > new_max_group_id) new_max_group_id = gid;
        }

        int total_groups = namesMap.size();

        // Build metadata JSON
        std::string metadata = "{\"population_size\":" + to_string(frame->size())
            + ",\"num_groups\":" + to_string(total_groups)
            + ",\"max_group_id\":" + to_string(new_max_group_id)
            + ",\"hash_size\":31"
            + "}";

        // 7. Write updated .dbri
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

        // Merge existing and new feature sets for GROUP_TO_FEATURE_SET
        if (dbri.has_section(DBRISection::GROUP_TO_FEATURE_SET)) {
            phmap::parallel_flat_hash_map<std::string, phmap::flat_hash_set<uint64_t>> existingFeatureSets;
            dbri.load_group_to_feature_set(existingFeatureSets);
            for (auto& [name, kset] : newGroupName_to_featureSet) {
                existingFeatureSets[name] = kset;
            }
            out.write_group_to_feature_set(existingFeatureSets);
        } else {
            phmap::parallel_flat_hash_map<std::string, phmap::flat_hash_set<uint64_t>> allFeatureSets;
            for (auto& [name, kset] : newGroupName_to_featureSet) {
                allFeatureSets[name] = kset;
            }
            out.write_group_to_feature_set(allFeatureSets);
        }

        // Embed JSON files if they exist
        if (dbri.has_section(DBRISection::RAW_GENE_SETS)) {
            std::string existing_raw = dbri.load_raw_gene_sets();
            out.write_raw_gene_sets(existing_raw);
        }
        if (dbri.has_section(DBRISection::HASHED_GENE_SETS)) {
            std::string existing_hashed = dbri.load_hashed_gene_sets();
            out.write_hashed_gene_sets(existing_hashed);
        }

        out.finalize_write();

        cerr << "[append] Wrote updated index to " << output_dbri_path << endl;
        cerr << "[append] Total groups: " << total_groups
             << " (" << total_existing_groups << " existing + " << new_groups_count << " new)" << endl;

        delete legend;
        delete frame;
    }

}
