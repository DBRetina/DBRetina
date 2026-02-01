#include "dbretina_core.hpp"
#include <iostream>
#include <cstdint>
#include <chrono>
#include "DBRetina_kf.hpp"
#include "parallel_hashmap/phmap.h"
#include <glob.h>
#include <vector>
#include <stdexcept>
#include <sstream>
#include <string>
#include <fstream>
#include "cpp-json/json.h"
#include "zstr.hpp"
#include "parallel_hashmap/phmap.h"
#include "parallel_hashmap/phmap_dump.h"
#include "DBRetina.hpp"
#include "DBRetinaIndex.hpp"



namespace dbretina {

    void dbretina_indexing(string json_file, string user_index_prefix) {

        DBRetina_PHMAP* frame;
        while (json_file.size() > 0 && json_file[json_file.size() - 1] == '/') json_file.erase(json_file.size() - 1, 1);

        std::string json_prefix = json_file.substr(json_file.find_last_of("/\\") + 1);

        while (json_prefix.size() > 0 && json_prefix[json_prefix.size() - 1] == '/') {
            json_prefix.erase(json_prefix.size() - 1, 1);
        }

        // remove json extension
        json_prefix = json_prefix.substr(0, json_prefix.find_last_of("."));
        

        flat_hash_map<string, string> namesMap;
        int selective_hashSize = 31;

        flat_hash_map<string, uint64_t> tagsMap;
        flat_hash_map<string, uint64_t> groupNameMap;
        auto* legend = new flat_hash_map<uint64_t, std::vector<uint32_t>>();
        flat_hash_map<uint64_t, uint64_t> colorsCount;
        uint64_t readID = 0, groupID = 1;
        string seqName, groupName;
        string line;
        priority_queue<uint64_t, vector<uint64_t>, std::greater<uint64_t>> freeColors;
        flat_hash_map<string, uint64_t> groupCounter;
        int detected_hashSize = 0;

        int total_groups_number = 0;
        frame = new DBRetina_PHMAP(selective_hashSize);

        flat_hash_map<string, uint32_t> groupName_to_featureCount;

        // hashed_MAP* groupName_to_featureSet = new(hashed_MAP);
        hashed_MAP groupName_to_featureSet; // = new(hashed_MAP);

        parse_dbretina_json(json_file, &groupName_to_featureSet);


        for (auto& group : groupName_to_featureSet) {
            total_groups_number++;
            string group_name = group.first;

            // Here we can decide
            seqName = group_name;
            groupName = group_name;

            namesMap.insert(make_pair(seqName, groupName));
            auto it = groupNameMap.find(groupName);
            groupCounter[groupName]++;
            if (it == groupNameMap.end()) {
                groupNameMap.insert(make_pair(groupName, groupID));
                tagsMap.insert(make_pair(to_string(groupID), groupID));
                vector<uint32_t> tmp;
                tmp.clear();
                tmp.push_back(groupID);
                legend->insert(make_pair(groupID, tmp));
                colorsCount.insert(make_pair(groupID, 0));
                groupID++;
            }
        }

        // DEBUG
        // cout << "namesmap construction done..." << endl;


        // ----------------------------------------------------------------


        flat_hash_map<uint64_t, string> inv_groupNameMap;
        for (auto& _ : groupNameMap)
            inv_groupNameMap[_.second] = _.first;


        int currIndex = 0;
        string feature;
        uint64_t tagBits = 0;
        uint64_t maxTagValue = (1ULL << tagBits) - 1;

        uint64_t lastTag = 0;
        readID = 0;

        int processed_sigs_count = 0;

        // START
        for (auto& group : groupName_to_featureSet) {
            string group_name = group.first;
            flat_hash_set<uint64_t>& featureSet = group.second;

            //START
            for (auto it = featureSet.begin(); it != featureSet.end(); ++it) {


                // cout << "Processing " << ++processed_sigs_count << "/" << total_groups_number << " | " << group_name << " ... " << endl;


                flat_hash_map<uint64_t, uint64_t> convertMap;

                string readName = group_name;
                string groupName = group_name;

                uint64_t groupTag = groupNameMap.find(groupName)->second;


                convertMap.clear();
                convertMap.insert(make_pair(0, groupTag));
                convertMap.insert(make_pair(groupTag, groupTag));

                auto loaded_features_it = featureSet.begin();
                groupName_to_featureCount[groupName] = featureSet.size();


                while (loaded_features_it != featureSet.end()) {
                    uint64_t hashed_feature = *loaded_features_it;
                    uint64_t currentTag = frame->getCount(hashed_feature);
                    auto itc = convertMap.find(currentTag);
                    if (itc == convertMap.end()) {
                        vector<uint32_t> colors = legend->find(currentTag)->second;
                        auto tmpiT = find(colors.begin(), colors.end(), groupTag);
                        if (tmpiT == colors.end()) {
                            colors.push_back(groupTag);
                            sort(colors.begin(), colors.end());
                        }

                        // [OPTIMIZE] [TODO] Optimize colors concatenation
                        // string colorsString = to_string(colors[0]);
                        // for (int k = 1; k < colors.size(); k++) {
                        //     colorsString += ";" + to_string(colors[k]);
                        // }

                        std::stringstream ss;
                        if (!colors.empty()) {
                            ss << colors[0];
                            for (size_t k = 1; k < colors.size(); ++k) {
                                ss << ";" << colors[k];
                            }
                        }
                        std::string colorsString = ss.str();
                        // END [OPTIMIZE] [TODO] Optimize colors concatenation


                        auto itTag = tagsMap.find(colorsString);

                        // START [OPTIMIZE] [TODO] Optimize
                        // if (itTag == tagsMap.end()) {
                        //     uint64_t newColor;
                        //     if (freeColors.size() == 0) {
                        //         newColor = groupID++;
                        //     }
                        //     else {
                        //         newColor = freeColors.top();
                        //         freeColors.pop();
                        //     }
                        //     tagsMap.insert(make_pair(colorsString, newColor));
                        //     legend->insert(make_pair(newColor, colors));
                        //     itTag = tagsMap.find(colorsString);
                        //     colorsCount[newColor] = 0;
                        // }

                        if (itTag == tagsMap.end()) {
                            uint64_t newColor;
                            if (freeColors.empty()) newColor = groupID++;
                            else {
                                newColor = freeColors.top();
                                freeColors.pop();
                            }
                            auto inserted = tagsMap.emplace(colorsString, newColor);
                            legend->emplace(newColor, colors);
                            itTag = inserted.first;
                            colorsCount[newColor] = 0;
                        }
                        // END [OPTIMIZE] [TODO] Optimize
                        uint64_t newColor = itTag->second;

                        convertMap.emplace(currentTag, newColor);
                        itc = convertMap.find(currentTag);
                    }
                    // START [OPTIMIZE] [TODO] Optimize
                    // if (itc->second != currentTag) {

                    //     colorsCount[currentTag]--;
                    //     if (colorsCount[currentTag] == 0 && currentTag != 0) {

                    //         auto _invGroupNameIT = inv_groupNameMap.find(currentTag);
                    //         if (_invGroupNameIT == inv_groupNameMap.end()) {
                    //             freeColors.push(currentTag);
                    //             vector<uint32_t> colors = legend->find(currentTag)->second;
                    //             string colorsString = to_string(colors[0]);
                    //             for (unsigned int k = 1; k < colors.size(); k++) {
                    //                 colorsString += ";" + to_string(colors[k]);
                    //             }
                    //             tagsMap.erase(colorsString);
                    //             legend->erase(currentTag);
                    //             if (convertMap.find(currentTag) != convertMap.end())
                    //                 convertMap.erase(currentTag);
                    //         }

                    //     }
                    //     colorsCount[itc->second]++;
                    // }

                    if (itc->second != currentTag) {
                        --colorsCount[currentTag];
                        if (colorsCount[currentTag] == 0 && currentTag != 0) {
                            auto _invGroupNameIT = inv_groupNameMap.find(currentTag);
                            if (_invGroupNameIT == inv_groupNameMap.end()) {
                                freeColors.push(currentTag);
                                auto legendIt = legend->find(currentTag);
                                std::vector<uint32_t> colors = legendIt->second;
                                std::stringstream ss;
                                ss << colors[0];
                                for (unsigned int k = 1; k < colors.size(); k++) {
                                    ss << ";" << colors[k];
                                }
                                std::string colorsString = ss.str();
                                tagsMap.erase(colorsString);
                                legend->erase(legendIt);
                                auto convertMapIt = convertMap.find(currentTag);
                                if (convertMapIt != convertMap.end()) {
                                    convertMap.erase(convertMapIt);
                                }
                            }
                        }
                        ++colorsCount[itc->second];
                    }

                    // END [OPTIMIZE] [TODO] Optimize



                    frame->setCount(hashed_feature, itc->second);
                    if (frame->getCount(hashed_feature) != itc->second) {
                        //frame->setC(kmer,itc->second);
                        cout << "Error Founded " << hashed_feature << " from group " << readName << " expected "
                            << itc->second << " found " << frame->getCount(hashed_feature) << endl;
                        exit(1);
                    }
                    loaded_features_it++;
                }
                readID += 1;
                groupCounter[groupName]--;
                if (colorsCount[groupTag] == 0) {
                    if (groupCounter[groupName] == 0) {
                        freeColors.push(groupTag);
                        legend->erase(groupTag);
                    }

                }
                // cout << "   saved_features(~" << frame->size() << ")." << endl;
                // cout << "   colors(~" << legend->size() << ")." << endl << endl;

                break;
            }
            // END

        }

        cerr << "[dev] Total number of colors: " << legend->size() << endl;

        string output_prefix = user_index_prefix;

        // Build groupID_to_featureCount
        flat_hash_map<uint32_t, uint32_t> groupID_to_featureCount;
        for (auto& [groupName, featureCount] : groupName_to_featureCount) {
            groupID_to_featureCount[groupNameMap[groupName]] = featureCount;
        }

        // Compute color statistics for logging
        size_t selected_colors = 0;
        for (auto& [color, sources] : *legend) {
            if (colorsCount[color] > 0) selected_colors++;
        }
        cerr << "[dev] Total selected colors: " << selected_colors << endl;

        // Find max group ID for metadata
        uint64_t max_group_id = 0;
        for (auto& [name, gid] : groupNameMap) {
            if (gid > max_group_id) max_group_id = gid;
        }

        // Build metadata JSON
        std::string metadata = "{\"population_size\":" + to_string(frame->size())
            + ",\"num_groups\":" + to_string(total_groups_number)
            + ",\"max_group_id\":" + to_string(max_group_id)
            + ",\"hash_size\":" + to_string(selective_hashSize)
            + "}";

        // Write unified .dbri index
        DBRetinaIndex dbri;
        std::string dbri_path = output_prefix + ".dbri";
        dbri.begin_write(dbri_path);

        dbri.write_phmap(frame);
        dbri.write_color_to_sources(*legend, colorsCount);
        dbri.write_color_count(colorsCount);
        dbri.write_group_feature_count(groupID_to_featureCount);
        dbri.write_names_map(namesMap, groupNameMap);
        dbri.write_metadata(metadata);
        dbri.write_tags_map(tagsMap);
        dbri.write_free_colors(freeColors);
        dbri.write_group_to_feature_set(groupName_to_featureSet);

        // Embed raw gene sets JSON if it exists
        std::string raw_json_path = output_prefix + "_raw.json";
        {
            std::ifstream raw_in(raw_json_path);
            if (raw_in.is_open()) {
                std::string raw_content((std::istreambuf_iterator<char>(raw_in)),
                                         std::istreambuf_iterator<char>());
                raw_in.close();
                dbri.write_raw_gene_sets(raw_content);
            }
        }

        // Embed hashed gene sets JSON
        {
            std::ifstream hashed_in(json_file);
            if (hashed_in.is_open()) {
                std::string hashed_content((std::istreambuf_iterator<char>(hashed_in)),
                                            std::istreambuf_iterator<char>());
                hashed_in.close();
                dbri.write_hashed_gene_sets(hashed_content);
            }
        }

        dbri.finalize_write();

        cerr << "[dev] Wrote unified index to " << dbri_path << endl;

        delete legend;
        delete frame;
    }

}