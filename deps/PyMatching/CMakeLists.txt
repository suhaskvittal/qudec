cmake_minimum_required(VERSION 3.4...3.18)
project(PyMatching2)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include_directories(src)
set(CMAKE_CXX_STANDARD 20 CACHE STRING "C++ version selection")
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(SIMD_WIDTH 128)
if (NOT(MSVC))
    if (CMAKE_SYSTEM_PROCESSOR MATCHES x86_64)
         set(ARCH_OPT "-O3" "-mno-avx2")
    else ()
         set(ARCH_OPT "-O3")
    endif ()
else ()
    set(ARCH_OPT "-O2")
    # https://stackoverflow.com/a/8591946 (RTC1 incompatible with O2 flag)
    STRING (REGEX REPLACE "/RTC(su|[1su])" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
endif ()

# Note: Stim is provided by the parent project, so we don't fetch it here

set(SOURCE_FILES_NO_MAIN
        src/pymatching/sparse_blossom/driver/namespaced_main.cc
        src/pymatching/sparse_blossom/driver/mwpm_decoding.cc
        src/pymatching/sparse_blossom/flooder/graph.cc
        src/pymatching/sparse_blossom/flooder/detector_node.cc
        src/pymatching/sparse_blossom/flooder_matcher_interop/compressed_edge.cc
        src/pymatching/sparse_blossom/flooder/graph_fill_region.cc
        src/pymatching/sparse_blossom/flooder/match.cc
        src/pymatching/sparse_blossom/flooder/graph_flooder.cc
        src/pymatching/sparse_blossom/matcher/alternating_tree.cc
        src/pymatching/sparse_blossom/matcher/mwpm.cc
        src/pymatching/sparse_blossom/flooder_matcher_interop/region_edge.cc
        src/pymatching/sparse_blossom/flooder_matcher_interop/mwpm_event.cc
        src/pymatching/sparse_blossom/tracker/flood_check_event.cc
        src/pymatching/sparse_blossom/diagram/animation_main.cc
        src/pymatching/sparse_blossom/diagram/mwpm_diagram.cc
        src/pymatching/sparse_blossom/search/search_graph.cc
        src/pymatching/sparse_blossom/search/search_detector_node.cc
        src/pymatching/sparse_blossom/search/search_flooder.cc
        src/pymatching/sparse_blossom/driver/user_graph.cc
        src/pymatching/rand/rand_gen.cc
        )



add_library(libpymatching SHARED ${SOURCE_FILES_NO_MAIN})
set_target_properties(libpymatching PROPERTIES PREFIX "")
target_include_directories(libpymatching PUBLIC src)
target_compile_options(libpymatching PRIVATE ${ARCH_OPT} -fPIC)
if(NOT(MSVC))
    target_link_options(libpymatching PRIVATE -pthread -O3)
endif()
target_link_libraries(libpymatching libstim)
install(TARGETS libpymatching LIBRARY DESTINATION lib)
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/src/" DESTINATION "include" FILES_MATCHING PATTERN "*.h" PATTERN "*.inl")
