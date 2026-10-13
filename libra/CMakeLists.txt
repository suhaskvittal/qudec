# Contributor: OpenAI Codex (Libra build entry point).
cmake_minimum_required(VERSION 3.20.3)
project(qudec_libra_compat LANGUAGES CXX)
# Compatibility entry point: the normal library owns the stable implementation.
enable_testing()
get_filename_component(QUDEC_ROOT "${CMAKE_CURRENT_LIST_DIR}/.." ABSOLUTE)
add_subdirectory("${QUDEC_ROOT}" qudec)
