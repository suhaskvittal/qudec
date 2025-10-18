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
add_subdirectory(stim)
add_subdirectory(blossom5)
add_subdirectory(PyMatching)

###################################################
# Source files
###################################################

set(QUDEC_FILES
    src/argparse.cpp
    src/decoding_graph.cpp
    src/decoder/surface_code.cpp
    src/io/dem.cpp
    src/gen.cpp
)

###################################################
# Quantum Decoder Library
###################################################

add_library(qudeclib SHARED ${QUDEC_FILES})
target_include_directories(qudeclib PUBLIC "src")
target_compile_options(qudeclib PUBLIC ${COMPILE_OPTIONS})
target_link_libraries(qudeclib PUBLIC libstim blossom5 libpymatching)
#target_compile_definitions(qudeclib PUBLIC -DDEBUG_DECODER)

###################################################
# Main executable
###################################################

add_executable(qudec main/qudec.cpp)
target_compile_options(qudec PRIVATE ${COMPILE_OPTIONS})
target_link_libraries(qudec PRIVATE qudeclib)

