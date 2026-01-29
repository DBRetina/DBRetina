# ----------------------------------------------------------------------------
# DBRetina project
# ----------------------------------------------------------------------------

# Minimum cmake version
cmake_minimum_required (VERSION 3.14...3.20)
# Name of project and that it is C/C++.
project (
        kSpider
        LANGUAGES CXX C
        VERSION 2.0.0
)

set(default_build_type "RELEASE")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -lstdc++fs -fPIC -lgomp -lrt -fopenmp -Ofast")
add_definitions(-w)

# ----------------------------------------------------------------------------
# Includes
# ----------------------------------------------------------------------------

include_directories("${PROJECT_SOURCE_DIR}/include")
if(DEFINED ENV{CONDA_PREFIX})
    include_directories("$ENV{CONDA_PREFIX}/include")
    link_directories("$ENV{CONDA_PREFIX}/lib")
endif()
include_directories("${PROJECT_SOURCE_DIR}/lib/argh")
include_directories("${PROJECT_SOURCE_DIR}/lib/json_parser/lib/include/")
add_subdirectory("lib/zstr")

set(PHMAP_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/lib/parallel-hashmap")
include_directories("${PHMAP_INCLUDE_DIRS}")

add_library(
        kSpider STATIC
        src/pairwise.cpp
        src/DBRetina_sketch.cpp
        src/query.cpp
        src/DBRetina_indexing.cpp
        src/DBRetina_GSA.cpp
        )

set_target_properties(kSpider PROPERTIES POSITION_INDEPENDENT_CODE 1 CXX_STANDARD 17)
target_link_libraries (kSpider z zstr::zstr)
target_include_directories(kSpider INTERFACE ${PHMAP_INCLUDE_DIRS})

add_executable(dump_colors apps/dump_kmer_count.cpp)
target_link_libraries(dump_colors kSpider z zstr::zstr)
