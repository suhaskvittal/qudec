# Copyright 2021 Google LLC
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

cmake_minimum_required(VERSION 3.13)
project(stim)
include_directories(src)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY out)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY out)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY out)

# Convert desired SIMD_WIDTH into machine architecture flags.

if (CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86_64|I386|ARM64)$")
    if(NOT(SIMD_WIDTH))
        set(MACHINE_FLAG "-march=native")
    elseif(SIMD_WIDTH EQUAL 256)
        set(MACHINE_FLAG "-mavx2" "-msse2")
    elseif(SIMD_WIDTH EQUAL 128)
        set(MACHINE_FLAG "-mno-avx2" "-msse2")
    elseif(SIMD_WIDTH EQUAL 64)
        set(MACHINE_FLAG "-mno-avx2" "-mno-sse2")
    endif()
else ()
    set(MACHINE_FLAG "")
endif()

# Define core library source files (excluding CLI commands, diagrams, and main executable)
set(STIM_CORE_SOURCES
    # Core stim functionality
    src/stim.cc

    # Circuit handling
    src/stim/circuit/circuit.cc
    src/stim/circuit/circuit_instruction.cc
    src/stim/circuit/gate_decomposition.cc
    src/stim/circuit/gate_target.cc

    # Detector error model
    src/stim/dem/dem_instruction.cc
    src/stim/dem/detector_error_model.cc

    # Gates
    src/stim/gates/gate_data_annotations.cc
    src/stim/gates/gate_data_blocks.cc
    src/stim/gates/gate_data_collapsing.cc
    src/stim/gates/gate_data_controlled.cc
    src/stim/gates/gate_data_hada.cc
    src/stim/gates/gate_data_heralded.cc
    src/stim/gates/gate_data_noisy.cc
    src/stim/gates/gate_data_pair_measure.cc
    src/stim/gates/gate_data_pauli.cc
    src/stim/gates/gate_data_pauli_product.cc
    src/stim/gates/gate_data_period_3.cc
    src/stim/gates/gate_data_period_4.cc
    src/stim/gates/gate_data_pp.cc
    src/stim/gates/gate_data_swaps.cc
    src/stim/gates/gates.cc

    # Circuit generation
    src/stim/gen/circuit_gen_params.cc
    src/stim/gen/gen_color_code.cc
    src/stim/gen/gen_rep_code.cc
    src/stim/gen/gen_surface_code.cc

    # I/O
    src/stim/io/measure_record.cc
    src/stim/io/measure_record_batch_writer.cc
    src/stim/io/measure_record_writer.cc
    src/stim/io/raii_file.cc
    src/stim/io/sparse_shot.cc
    src/stim/io/stim_data_formats.cc

    # Memory management
    src/stim/mem/bit_ref.cc
    src/stim/mem/simd_util.cc
    src/stim/mem/simd_word.cc
    src/stim/mem/sparse_xor_vec.cc

    # Search algorithms
    src/stim/search/graphlike/algo.cc
    src/stim/search/graphlike/edge.cc
    src/stim/search/graphlike/graph.cc
    src/stim/search/graphlike/node.cc
    src/stim/search/graphlike/search_state.cc
    src/stim/search/hyper/algo.cc
    src/stim/search/hyper/edge.cc
    src/stim/search/hyper/graph.cc
    src/stim/search/hyper/node.cc
    src/stim/search/hyper/search_state.cc
    src/stim/search/sat/wcnf.cc

    # Simulators
    src/stim/simulators/error_analyzer.cc
    src/stim/simulators/error_matcher.cc
    src/stim/simulators/force_streaming.cc
    src/stim/simulators/graph_simulator.cc
    src/stim/simulators/matched_error.cc
    src/stim/simulators/sparse_rev_frame_tracker.cc
    src/stim/simulators/vector_simulator.cc

    # Stabilizers
    src/stim/stabilizers/flex_pauli_string.cc

    # Utilities
    src/stim/util_bot/arg_parse.cc
    src/stim/util_bot/error_decomp.cc
    src/stim/util_bot/probability_util.cc

    # Core utility functions (excluding export utilities)
    src/stim/util_top/circuit_inverse_qec.cc
    src/stim/util_top/circuit_inverse_unitary.cc
    src/stim/util_top/circuit_to_detecting_regions.cc
    src/stim/util_top/circuit_vs_amplitudes.cc
    src/stim/util_top/has_flow.cc
    src/stim/util_top/mbqc_decomposition.cc
    src/stim/util_top/missing_detectors.cc
    src/stim/util_top/reference_sample_tree.cc
    src/stim/util_top/simplified_circuit.cc
    src/stim/util_top/transform_without_feedback.cc
)

# Create the library as a shared/dynamic library
add_library(libstim SHARED ${STIM_CORE_SOURCES})
set_target_properties(libstim PROPERTIES
    PREFIX ""
    POSITION_INDEPENDENT_CODE ON
)
# Set include directories as PUBLIC so they are available to consumers
target_include_directories(libstim PUBLIC src)
if(NOT(MSVC))
    target_compile_options(libstim PRIVATE -O3 -Wall -Wpedantic -fPIC -fno-strict-aliasing ${MACHINE_FLAG})
    target_link_options(libstim PRIVATE -O3)
else()
    target_compile_options(libstim PRIVATE ${MACHINE_FLAG} /O1)
endif()
install(TARGETS libstim LIBRARY DESTINATION)
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/src/" DESTINATION "include" FILES_MATCHING PATTERN "*.h" PATTERN "*.inl")
