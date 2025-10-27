/*
 *  author: Suhas Vittal
 *  date:   23 October 2025
 *
 *  EPR-based surface code circuit generation using sc_epr_generation
 * */

#include "argparse.h"
#include "decoder/surface_code.h"
#include "decoder_eval.h"
#include "gen/epr.h"
#include "qudec_common.h"

#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>

bool GL_DEBUG_DECODER{false};

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

int
main(int argc, char* argv[])
{
    int64_t     code_distance;
    int64_t     num_rounds;
    int64_t     num_trials;
    bool        do_memory_experiment;
    
    // EPR-specific parameters
    double      attenuation_rate;
    double      photonic_link_error;
    int64_t     hw1_round_ns;
    int64_t     hw2_round_ns;
    double      phys_error;
    
    std::string generated_stim_output_file;

    ARGPARSE()
        .optional("-d", "--code-distance", "code distance", code_distance, 3)
        .optional("-r", "--rounds", "number of rounds", num_rounds, 9)
        .optional("-t", "--trials", "number of trials to run", num_trials, 1000000)
        .optional("-m", "--memory-experiment", "do memory experiment (vs stability)", do_memory_experiment, true)
        
        // EPR-specific parameters:
        .optional("-a", "--attenuation-rate", "photon attenuation rate", attenuation_rate, 1e-2)
        .optional("-pl", "--photonic-link-error", "photonic link error rate", photonic_link_error, 1e-2)
        .optional("-hw1", "--hw1-round-ns", "HW1 round time in ns", hw1_round_ns, 1200)
        .optional("-hw2", "--hw2-round-ns", "HW2 round time in ns", hw2_round_ns, 1'200'000)
        .optional("-p", "--phys-error", "physical error rate", phys_error, 1e-3)
        
        // output file:
        .optional("-o", "--output", "output file for generated stim circuit", generated_stim_output_file, "generated.stim.out")
        // decoding:
        .optional("-dd", "--debug-decoder", "set flag debug decoder flag", GL_DEBUG_DECODER, false)
    
        .parse(argc, argv);

    // Create EPR generation configuration
    gen::EPR_GEN_CONFIG config;
    config.attenuation_rate = attenuation_rate;
    config.photonic_link_error = photonic_link_error;
    config.hw1_round_ns = static_cast<uint64_t>(hw1_round_ns);
    config.hw2_round_ns = static_cast<uint64_t>(hw2_round_ns);
    config.phys_error = phys_error;

    // Print configuration
    std::cout << "======================== EPR GENERATION CONFIG ==========================\n";
    print_stat(std::cout, "CODE_DISTANCE", code_distance);
    print_stat(std::cout, "NUM_ROUNDS", num_rounds);
    print_stat(std::cout, "NUM_TRIALS", num_trials);
    print_stat(std::cout, "MEMORY_EXPERIMENT", do_memory_experiment ? "true" : "false");
    print_stat(std::cout, "ATTENUATION_RATE", config.attenuation_rate);
    print_stat(std::cout, "PHOTONIC_LINK_ERROR", config.photonic_link_error);
    print_stat(std::cout, "HW1_ROUND_NS", config.hw1_round_ns);
    print_stat(std::cout, "HW2_ROUND_NS", config.hw2_round_ns);
    print_stat(std::cout, "PHYS_ERROR", config.phys_error);
    print_stat(std::cout, "OUTPUT_FILE", generated_stim_output_file);
    std::cout << "======================================================================\n";

    // Generate the EPR-based surface code circuit
    stim::Circuit circuit = gen::sc_epr_generation(config, num_rounds, code_distance, do_memory_experiment);

    // Print circuit if small enough
    if (code_distance <= 3)
    {
        std::cout << "======================== GENERATED CIRCUIT ==========================\n";
        std::cout << circuit << "\n";
        std::cout << "=====================================================================\n";
    }

    // Write circuit to output file
    std::cout << "Writing circuit to " << generated_stim_output_file << "...\n";
    std::ofstream out(generated_stim_output_file);
    if (!out.is_open())
    {
        std::cerr << "Error: Could not open output file " << generated_stim_output_file << std::endl;
        return 1;
    }
    
    out << circuit.str() << "\n";
    out.close();

    // Print circuit statistics
    std::cout << "EPR circuit generation completed successfully!\n";

    // Run PyMatching decoder on the generated circuit
    std::cout << "\n======================== RUNNING PYMATCHING DECODER ==========================\n";
    print_stat(std::cout, "NUM_TRIALS", num_trials);
    std::cout << "Running decoder...\n";

    DECODER_STATS stats = eval_decoder<PYMATCHING>(circuit, num_trials);

    // Calculate and print results
    double ler = fpdiv(stats.errors, stats.trials);
    double mean_time_us = fpdiv(stats.total_time_us, stats.trials);
    double mean_time_us_nontrivial = fpdiv(stats.total_time_us, stats.trials - stats.trivial_trials);

    std::cout << "======================== DECODER RESULTS ==========================\n";
    print_stat(std::cout, "LOGICAL_ERRORS", stats.errors);
    print_stat(std::cout, "TRIALS", stats.trials);
    print_stat(std::cout, "TRIVIAL_TRIALS", stats.trivial_trials);
    print_stat(std::cout, "LOGICAL_ERROR_RATE", ler);
    print_stat(std::cout, "MEAN_TIME_US", mean_time_us);
    print_stat(std::cout, "MEAN_TIME_US_NONTRIVIAL", mean_time_us_nontrivial);
    std::cout << "===============================================================\n";

    return 0;
}

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
