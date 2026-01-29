#pragma once

#include <string>
#include <cstdint>
#include <fstream>
#include <iostream>
#include "parallel_hashmap/phmap.h"
#include "parallel_hashmap/phmap_dump.h"

#define PHMAP_SUBMAPS 6

using DBRETINA_PHMAP = phmap::parallel_flat_hash_map<uint64_t, uint64_t,
    std::hash<uint64_t>, std::equal_to<uint64_t>,
    std::allocator<std::pair<const uint64_t, uint64_t>>,
    PHMAP_SUBMAPS>;

struct DBRetina_PHMAP {

    DBRETINA_PHMAP MAP;
    uint64_t kSize;

    DBRetina_PHMAP(uint64_t kSize) : kSize(kSize) {}

    void setCount(uint64_t key, uint64_t value) {
        MAP[key] = value;
    }

    uint64_t getCount(uint64_t key) {
        auto it = MAP.find(key);
        if (it != MAP.end()) return it->second;
        return 0;
    }

    uint64_t getkSize() { return kSize; }

    uint64_t size() { return MAP.size(); }

    auto begin() { return MAP.begin(); }
    auto end() { return MAP.end(); }

    void save(const std::string& prefix) {
        // Write .extra metadata
        std::ofstream extra(prefix + ".extra");
        extra << kSize << "\n";
        extra << 0 << "\n";         // hash_mode  (mumur_hasher)
        extra << 1 << "\n";         // slicing_mode
        extra << kSize << ":0\n";   // params_to_string
        extra << "features:" << kSize << "\n";
        extra.close();

        // Write .phmap binary dump
        phmap::BinaryOutputArchive ar_out((prefix + ".phmap").c_str());
        MAP.phmap_dump(ar_out);
    }

    static DBRetina_PHMAP* load(const std::string& prefix) {
        // Read .extra metadata
        std::ifstream extra(prefix + ".extra");
        if (!extra.is_open()) {
            std::cerr << "Error: cannot open " << prefix << ".extra" << std::endl;
            exit(1);
        }
        std::string line;
        std::getline(extra, line);
        uint64_t kSize = std::stoull(line);
        extra.close();

        auto* frame = new DBRetina_PHMAP(kSize);

        // Read .phmap binary dump
        phmap::BinaryInputArchive ar_in((prefix + ".phmap").c_str());
        frame->MAP.phmap_load(ar_in);

        return frame;
    }
};
