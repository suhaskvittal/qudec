cmake_minimum_required(VERSION 3.16)
project(tesseract_decoder LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(FetchContent)
find_package(Threads REQUIRED)

# === External dependencies ===
# Stim
FetchContent_Declare(
  stim
  GIT_REPOSITORY https://github.com/quantumlib/stim.git
  GIT_TAG bd60b73525fd5a9b30839020eb7554ad369e4337
)
FetchContent_MakeAvailable(stim)

# OpenAI GPT-6: Local packed bitsets replace the former Boost headers.

set(OPT_COPTS -Ofast -fno-fast-math -march=native)

set(TESSERACT_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)

# === Libraries ===
add_library(common ${TESSERACT_SRC_DIR}/common.cc ${TESSERACT_SRC_DIR}/common.h)
target_include_directories(common PUBLIC ${TESSERACT_SRC_DIR})
target_compile_options(common PRIVATE ${OPT_COPTS})
target_link_libraries(common PUBLIC libstim Threads::Threads)

add_library(utils ${TESSERACT_SRC_DIR}/utils.cc ${TESSERACT_SRC_DIR}/utils.h)
target_include_directories(utils PUBLIC ${TESSERACT_SRC_DIR})
target_compile_options(utils PRIVATE ${OPT_COPTS})
target_link_libraries(utils PUBLIC common libstim Threads::Threads)

add_library(visualization ${TESSERACT_SRC_DIR}/visualization.cc ${TESSERACT_SRC_DIR}/visualization.h)
target_include_directories(visualization PUBLIC ${TESSERACT_SRC_DIR})
target_compile_options(visualization PRIVATE ${OPT_COPTS})
target_link_libraries(visualization PUBLIC common)

add_library(tesseract_lib
  ${TESSERACT_SRC_DIR}/tesseract.cc
  ${TESSERACT_SRC_DIR}/tesseract.h
)
target_include_directories(tesseract_lib PUBLIC ${TESSERACT_SRC_DIR})
target_compile_options(tesseract_lib PRIVATE ${OPT_COPTS})
target_link_libraries(tesseract_lib PUBLIC common utils visualization libstim)

# Kept for research purposes (layered beam-search decoder). Not currently
# used by any executable in this tree, but built here to make sure it still
# compiles against the trimmed common/utils/libstim set above.
add_library(tesseract_trellis_lib
  ${TESSERACT_SRC_DIR}/tesseract_trellis.cc
  ${TESSERACT_SRC_DIR}/tesseract_trellis.h
)
target_include_directories(tesseract_trellis_lib PUBLIC ${TESSERACT_SRC_DIR})
target_compile_options(tesseract_trellis_lib PRIVATE ${OPT_COPTS})
target_link_libraries(tesseract_trellis_lib PUBLIC common utils libstim Threads::Threads)
