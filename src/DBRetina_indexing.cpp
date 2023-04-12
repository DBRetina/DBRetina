#include "kSpider.hpp"
#include <iostream>
#include <cstdint>
#include <chrono>
#include "colored_kDataFrame.hpp"
#include "parallel_hashmap/phmap.h"
#include "kDataFrame.hpp"
#include "algorithms.hpp"
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



namespace kSpider {

    void dbretina_indexing(string json_file) {

        kDataFrame* frame;
        while (json_file.size() > 0 && json_file[json_file.size() - 1] == '/') json_file.erase(json_file.size() - 1, 1);

        std::string json_prefix = json_file.substr(json_file.find_last_of("/\\") + 1);

        while (json_prefix.size() > 0 && json_prefix[json_prefix.size() - 1] == '/') {
            json_prefix.erase(json_prefix.size() - 1, 1);
        }

        // remove json extension
        json_prefix = json_prefix.substr(0, json_prefix.find_last_of("."));

        cout << "json_prefix: " << json_prefix << endl;

        flat_hash_map<string, string> namesMap;
        int selective_kSize = 31;

        flat_hash_map<string, uint64_t> tagsMap;
        flat_hash_map<string, uint64_t> groupNameMap;
        auto* legend = new flat_hash_map<uint64_t, std::vector<uint32_t>>();
        flat_hash_map<uint64_t, uint64_t> colorsCount;
        uint64_t readID = 0, groupID = 1;
        string seqName, groupName;
        string line;
        priority_queue<uint64_t, vector<uint64_t>, std::greater<uint64_t>> freeColors;
        flat_hash_map<string, uint64_t> groupCounter;
        int detected_kSize = 0;

        int total_groups_number = 0;
        frame = new kDataFramePHMAP(selective_kSize, mumur_hasher);

        flat_hash_map<string, uint32_t> groupName_to_kmerCount;

        // hashed_MAP* groupName_to_kmerSet = new(hashed_MAP);
        hashed_MAP groupName_to_kmerSet; // = new(hashed_MAP);

        parse_dbretina_json(json_file, &groupName_to_kmerSet);


        for (auto& group : groupName_to_kmerSet) {
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

        cout << "namesmap construction done..." << endl;


        // ----------------------------------------------------------------


        flat_hash_map<uint64_t, string> inv_groupNameMap;
        for (auto& _ : groupNameMap)
            inv_groupNameMap[_.second] = _.first;


        int currIndex = 0;
        string kmer;
        uint64_t tagBits = 0;
        uint64_t maxTagValue = (1ULL << tagBits) - 1;

        uint64_t lastTag = 0;
        readID = 0;

        int processed_sigs_count = 0;

        // START
        for (auto& group : groupName_to_kmerSet) {
            total_groups_number++;
            string group_name = group.first;
            parallel_flat_hash_set<uint64_t>& kmerSet = group.second;

            //START
            for (auto it = kmerSet.begin(); it != kmerSet.end(); ++it) {


                cout << "Processing " << ++processed_sigs_count << "/" << total_groups_number << " | " << group_name << " k:" << selective_kSize << " ... " << endl;


                flat_hash_map<uint64_t, uint64_t> convertMap;

                string readName = group_name;
                string groupName = group_name;

                uint64_t readTag = groupNameMap.find(groupName)->second;


                convertMap.clear();
                convertMap.insert(make_pair(0, readTag));
                convertMap.insert(make_pair(readTag, readTag));

                auto loaded_kmers_it = kmerSet.begin();
                groupName_to_kmerCount[groupName] = kmerSet.size();


                while (loaded_kmers_it != kmerSet.end()) {
                    uint64_t hashed_kmer = *loaded_kmers_it;
                    uint64_t currentTag = frame->getCount(hashed_kmer);
                    auto itc = convertMap.find(currentTag);
                    if (itc == convertMap.end()) {
                        vector<uint32_t> colors = legend->find(currentTag)->second;
                        auto tmpiT = find(colors.begin(), colors.end(), readTag);
                        if (tmpiT == colors.end()) {
                            colors.push_back(readTag);
                            sort(colors.begin(), colors.end());
                        }

                        string colorsString = to_string(colors[0]);
                        for (int k = 1; k < colors.size(); k++) {
                            colorsString += ";" + to_string(colors[k]);
                        }

                        auto itTag = tagsMap.find(colorsString);
                        if (itTag == tagsMap.end()) {
                            uint64_t newColor;
                            if (freeColors.size() == 0) {
                                newColor = groupID++;
                            }
                            else {
                                newColor = freeColors.top();
                                freeColors.pop();
                            }

                            tagsMap.insert(make_pair(colorsString, newColor));
                            legend->insert(make_pair(newColor, colors));
                            itTag = tagsMap.find(colorsString);
                            colorsCount[newColor] = 0;
                        }
                        uint64_t newColor = itTag->second;

                        convertMap.insert(make_pair(currentTag, newColor));
                        itc = convertMap.find(currentTag);
                    }

                    if (itc->second != currentTag) {

                        colorsCount[currentTag]--;
                        if (colorsCount[currentTag] == 0 && currentTag != 0) {

                            auto _invGroupNameIT = inv_groupNameMap.find(currentTag);
                            if (_invGroupNameIT == inv_groupNameMap.end()) {
                                freeColors.push(currentTag);
                                vector<uint32_t> colors = legend->find(currentTag)->second;
                                string colorsString = to_string(colors[0]);
                                for (unsigned int k = 1; k < colors.size(); k++) {
                                    colorsString += ";" + to_string(colors[k]);
                                }
                                tagsMap.erase(colorsString);
                                legend->erase(currentTag);
                                if (convertMap.find(currentTag) != convertMap.end())
                                    convertMap.erase(currentTag);
                            }

                        }
                        colorsCount[itc->second]++;
                    }

                    frame->setCount(hashed_kmer, itc->second);
                    if (frame->getCount(hashed_kmer) != itc->second) {
                        //frame->setC(kmer,itc->second);
                        cout << "Error Founded " << hashed_kmer << " from sequence " << readName << " expected "
                            << itc->second << " found " << frame->getCount(hashed_kmer) << endl;
                        exit(1);
                    }
                    loaded_kmers_it++;
                }
                readID += 1;
                groupCounter[groupName]--;
                if (colorsCount[readTag] == 0) {
                    if (groupCounter[groupName] == 0) {
                        freeColors.push(readTag);
                        legend->erase(readTag);
                    }

                }
                cout << "   saved_kmers(~" << frame->size() << ")." << endl;
                cout << "   colors(~" << legend->size() << ")." << endl << endl;

                break;
            }
            // END

        }


        string output_prefix = "idx_" + json_prefix;

        // Dump kmer count
        flat_hash_map<uint32_t, uint32_t> groupID_to_kmerCount;
        for (auto& [groupName, kmerCount] : groupName_to_kmerCount) {
            groupID_to_kmerCount[groupNameMap[groupName]] = kmerCount;
        }

        phmap::BinaryOutputArchive ar_out(string(output_prefix + "_groupID_to_kmerCount.bin").c_str());
        groupID_to_kmerCount.phmap_dump(ar_out);


        // Dump color->sources
        auto color_to_sources = new phmap::flat_hash_map<uint64_t, phmap::flat_hash_set<uint32_t>>();
        for (auto it : *legend) {
            phmap::flat_hash_set<uint32_t> tmp(std::make_move_iterator(it.second.begin()), std::make_move_iterator(it.second.end()));
            color_to_sources->operator[](it.first) = tmp;
        }


        phmap::BinaryOutputArchive ar_out_1(string(output_prefix + "_color_to_sources.bin").c_str());
        ar_out_1.saveBinary(color_to_sources->size());
        for (auto& [k, v] : *color_to_sources)
        {
            ar_out_1.saveBinary(k);
            ar_out_1.saveBinary(v);
        }

        // Dump colors count
        phmap::BinaryOutputArchive ar_out_3(string(output_prefix + "_color_count.bin").c_str());
        colorsCount.phmap_dump(ar_out_3);

        // export namesMap
        ofstream namesMapOut(output_prefix + ".namesMap");
        namesMapOut << namesMap.size() << endl;
        for (auto it : namesMap)
        {
            namesMapOut << groupNameMap[it.second] << " " << it.second << endl;
        }
        namesMapOut.close();

        // Write extra info
        ofstream file(output_prefix + ".extra");
        file << selective_kSize << endl;
        file << frame->KD->hash_mode << endl;
        file << frame->KD->slicing_mode << endl;
        file << frame->KD->params_to_string() << endl;
        file.close();

        // ------- Pause serializing index for now.
        /*

        colorTable* colors = new intVectorsTable();
        for (auto it : *legend) {
            colors->setColor(it.first, it.second);
        }

        colored_kDataFrame* res = new colored_kDataFrame();
        res->setColorTable(colors);
        res->setkDataFrame(frame);
        for (auto iit = namesMap.begin(); iit != namesMap.end(); iit++) {
            uint32_t sampleID = groupNameMap[iit->second];
            res->namesMap[sampleID] = iit->second;
            res->namesMapInv[iit->second] = sampleID;
        }
        cout << "saving to " << json_prefix << " ..." << endl;
        res->save(json_prefix);
        */

        // ------ END Pause serializing index for now.

    }

}