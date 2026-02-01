# ----------------------------------------------------------------------------
# DBRetina project
# ----------------------------------------------------------------------------

cmake_minimum_required(VERSION 3.15...3.27)
project(
        DBRetina
        LANGUAGES CXX C
        VERSION 2.0.0
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# ----------------------------------------------------------------------------
# Python + nanobind
# ----------------------------------------------------------------------------

find_package(Python 3.11 COMPONENTS Interpreter Development.Module REQUIRED)

include(FetchContent)
FetchContent_Declare(
    nanobind
    GIT_REPOSITORY https://github.com/wjakob/nanobind.git
    GIT_TAG        v2.4.0
    GIT_SHALLOW    TRUE
)
FetchContent_MakeAvailable(nanobind)
find_package(nanobind CONFIG REQUIRED PATHS ${nanobind_SOURCE_DIR}/cmake)

# ----------------------------------------------------------------------------
# Includes
# ----------------------------------------------------------------------------

include_directories("${PROJECT_SOURCE_DIR}/include")
if(DEFINED ENV{CONDA_PREFIX})
    include_directories("$ENV{CONDA_PREFIX}/include")
    link_directories("$ENV{CONDA_PREFIX}/lib")
endif()
include_directories("${PROJECT_SOURCE_DIR}/lib/json_parser/lib/include/")
add_subdirectory("lib/zstr")

set(PHMAP_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/lib/parallel-hashmap")
include_directories("${PHMAP_INCLUDE_DIRS}")

# ----------------------------------------------------------------------------
# Static library
# ----------------------------------------------------------------------------

add_library(
        dbretina_core STATIC
        src/pairwise.cpp
        src/DBRetina_sketch.cpp
        src/query.cpp
        src/DBRetina_indexing.cpp
        src/DBRetina_GSA.cpp
        src/DBRetinaIndex.cpp
        src/DBRetina_append.cpp
        src/DBRetina_merge.cpp
        src/DBRetinaPairwise.cpp
        )

set_target_properties(dbretina_core PROPERTIES CXX_STANDARD 17)
target_compile_options(dbretina_core PRIVATE -fopenmp -O3 -Ofast -w)
target_link_libraries(dbretina_core z zstr::zstr)
target_include_directories(dbretina_core INTERFACE ${PHMAP_INCLUDE_DIRS})

# ----------------------------------------------------------------------------
# nanobind Python module
# ----------------------------------------------------------------------------

nanobind_add_module(_dbretina_internal src/python_bindings.cpp)
target_link_libraries(_dbretina_internal PRIVATE dbretina_core z bz2 gomp dl)
target_include_directories(_dbretina_internal PRIVATE
    ${PROJECT_SOURCE_DIR}/include
    ${PHMAP_INCLUDE_DIRS}
    ${PROJECT_SOURCE_DIR}/lib/json_parser/lib/include
    ${PROJECT_SOURCE_DIR}/lib/zstr/src
)
target_compile_options(_dbretina_internal PRIVATE -fopenmp -O3 -Ofast -w)

install(TARGETS _dbretina_internal LIBRARY DESTINATION .)
