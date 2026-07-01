cmake_minimum_required(VERSION 3.20.3)
project(qudec VERSION 1.0.0 LANGUAGES CXX)

# initialize compile options
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(COMPILE_OPTIONS
            -ggdb3
            -Wall
            -Wextra
            -Wpedantic
            -Wunused
            )
else()
    set(COMPILE_OPTIONS
            -Ofast
            -flto=auto
            -fstrict-aliasing
            -march=native
            -Wunused
            )
endif()

# Set C++ standard
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

###################################################
# Dependencies
###################################################

# Add subdirectories for dependencies
add_subdirectory(deps/stim)
add_subdirectory(deps/blossom5)
add_subdirectory(deps/PyMatching)

option(ENABLE_MPI "Enable MPI parallelization" OFF)
if(ENABLE_MPI)
    find_package(MPI REQUIRED)
endif()

###################################################
# Source files
###################################################

set(QUDEC_FILES
    src/decoder/surface_code.cpp
    src/decoder/surface_code/cluster_match.cpp
    src/decoder/surface_code/matching_data.cpp
    src/globals.cpp
    src/sampler.cpp
)


###################################################
# Quantum Decoder Library
###################################################

add_library(qudeclib SHARED ${QUDEC_FILES})
target_include_directories(qudeclib PUBLIC "src" "deps")
target_compile_options(qudeclib PUBLIC ${COMPILE_OPTIONS})
target_link_libraries(qudeclib PUBLIC libstim libpymatching blossom5)
if(ENABLE_MPI)
    target_compile_definitions(qudeclib PUBLIC ENABLE_MPI)
    target_link_libraries(qudeclib PUBLIC MPI::MPI_CXX)
endif()

###################################################
# Verilator (hardware emulation)
###################################################

option(VERILATOR_CLUSTER_MATCH "Enable Verilator hardware emulation for CLUSTER_MATCH" OFF)

if(VERILATOR_CLUSTER_MATCH)
    find_package(verilator REQUIRED HINTS $ENV{VERILATOR_ROOT})

    verilate(qudeclib
        SOURCES systemverilog/cluster_match/filter.sv
        INCLUDE_DIRS systemverilog
        TOP_MODULE FILTER
        VERILATOR_ARGS --cc --assert --x-assign unique --x-initial unique
    )
    verilate(qudeclib
        SOURCES systemverilog/cluster_match/initialize_neighbors.sv
        INCLUDE_DIRS systemverilog
        TOP_MODULE INITIALIZE_NEIGHBORS
        VERILATOR_ARGS --cc --assert --x-assign unique --x-initial unique
    )

    verilate(qudeclib
        SOURCES systemverilog/cluster_match/astrea.sv
        INCLUDE_DIRS systemverilog
        TOP_MODULE ASTREA
        VERILATOR_ARGS --cc --assert --x-assign unique --x-initial unique
    )

    target_compile_definitions(qudeclib PUBLIC VERILATOR_CLUSTER_MATCH)

    add_executable(qudec_verilator main/qudec_verilator.cpp)
    target_compile_options(qudec_verilator PRIVATE ${COMPILE_OPTIONS})
    target_link_libraries(qudec_verilator PRIVATE qudeclib)
endif()

###################################################
# Main executable
###################################################

add_executable(qudec main/qudec.cpp)
target_compile_options(qudec PRIVATE ${COMPILE_OPTIONS})
target_link_libraries(qudec PRIVATE qudeclib)


