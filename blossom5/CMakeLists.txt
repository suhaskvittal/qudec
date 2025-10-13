# CMakeLists.txt for Blossom5 library
# Minimum cost perfect matching algorithm implementation

cmake_minimum_required(VERSION 3.13)
project(blossom5)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set optimization flags
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Compiler flags
if(NOT MSVC)
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
    set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
else()
    set(CMAKE_CXX_FLAGS_RELEASE "/O2 /DNDEBUG")
    set(CMAKE_CXX_FLAGS_DEBUG "/Od /Zi")
endif()

# Source files from main directory
set(BLOSSOM5_MAIN_SOURCES
    PMduals.cpp
    PMexpand.cpp
    PMinit.cpp
    PMinterface.cpp
    PMmain.cpp
    PMrepair.cpp
    PMshrink.cpp
    misc.cpp
)

# Source files from GEOM subdirectory
set(BLOSSOM5_GEOM_SOURCES
    GEOM/GPMinit.cpp
    GEOM/GPMinterface.cpp
    GEOM/GPMkdtree.cpp
    GEOM/GPMmain.cpp
)

# Source files from MinCost subdirectory
set(BLOSSOM5_MINCOST_SOURCES
    MinCost/MinCost.cpp
)

# Header files
set(BLOSSOM5_HEADERS
    PerfectMatching.h
    PMimplementation.h
    block.h
    LCA.h
    PQ.h
    timer.h
    GEOM/GeomPerfectMatching.h
    GEOM/GPMkdtree.h
    MinCost/MinCost.h
)

# Combine all sources
set(BLOSSOM5_ALL_SOURCES
    ${BLOSSOM5_MAIN_SOURCES}
    ${BLOSSOM5_GEOM_SOURCES}
    ${BLOSSOM5_MINCOST_SOURCES}
)

# Create the library
add_library(blossom5 STATIC ${BLOSSOM5_ALL_SOURCES})

# Set include directories
target_include_directories(blossom5 PUBLIC 
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/GEOM
    ${CMAKE_CURRENT_SOURCE_DIR}/MinCost
)

# Platform-specific libraries
if(NOT APPLE)
    # Add librt for non-Darwin platforms (Linux, etc.)
    find_library(RT_LIBRARY rt)
    if(RT_LIBRARY)
        target_link_libraries(blossom5 ${RT_LIBRARY})
    endif()
endif()

# Set target properties
set_target_properties(blossom5 PROPERTIES
    OUTPUT_NAME "blossom5"
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
)

# Install rules
install(TARGETS blossom5
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
)

install(FILES ${BLOSSOM5_HEADERS}
    DESTINATION include/blossom5
)

# Optional: Create an executable for testing (similar to original Makefile)
option(BUILD_BLOSSOM5_EXECUTABLE "Build blossom5 executable" OFF)

if(BUILD_BLOSSOM5_EXECUTABLE)
    # Check if example.cpp exists for the executable
    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/example.cpp")
        add_executable(blossom5_exe example.cpp)
        target_link_libraries(blossom5_exe blossom5)
        set_target_properties(blossom5_exe PROPERTIES OUTPUT_NAME "blossom5")
    endif()
endif()
