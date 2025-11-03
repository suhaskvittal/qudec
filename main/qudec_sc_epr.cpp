/*
 *  author: Suhas Vittal
 *  date:   23 October 2025
 *
 *  EPR-based surface code circuit generation using sc_epr_generation
 * */

#include "argparse.h"
#include "decoder/epr_pym.h"
#include "decoder/surface_code.h"
#include "decoder_eval.h"
#include "gen/epr.h"
#include "qudec_common.h"

#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

void
write_stim_circuit_to_file(std::string_view filename, const stim::Circuit& circuit)
{
    std::ofstream out(std::string{filename});
    if (!out.is_open())
    {
        std::cerr << "Error: Could not open output file " << filename << std::endl;
        return;
    }
    out << circuit.str() << "\n";
    out.close();
}

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

int
main(int argc, char* argv[])
{
    int64_t     code_distance;
    int64_t     num_rounds;
    int64_t     num_trials;
    int64_t     num_errors;
    std::string experiment_type;
    
    // EPR-specific parameters
    double      attenuation_rate;
    double      photonic_link_error;
    int64_t     hw1_round_ns;
    int64_t     hw2_round_ns;
    double      phys_error;

    int64_t eval_mode;  // 0 = use pymatching on global,
                        // 1 = use dual pass
                        // -1 = evaluate single hardware EPR only.

    ARGPARSE()
        .optional("-d", "--code-distance", "code distance", code_distance, 3)
        .optional("-r", "--rounds", "number of rounds", num_rounds, 9)
        .optional("-t", "--trials", "number of trials to run", num_trials, 1000000)
        .optional("-k", "--stop-after-errors", "stop after this many errors", num_errors, 25)
        .optional("", "--experiment", "experiment type", experiment_type, "memory")
        
        // EPR-specific parameters:
        .optional("-a", "--attenuation-rate", "photon attenuation rate", attenuation_rate, 1e-2)
        .optional("-pl", "--photonic-link-error", "photonic link error rate", photonic_link_error, 1e-2)
        .optional("-hw1", "--hw1-round-ns", "HW1 round time in ns", hw1_round_ns, 1200)
        .optional("-hw2", "--hw2-round-ns", "HW2 round time in ns", hw2_round_ns, 1'200'000)
        .optional("-p", "--phys-error", "physical error rate", phys_error, 1e-3)
        
        // decoding:
        .optional("-dd", "--debug-decoder", "set flag debug decoder flag", GL_DEBUG_DECODER, false)
        .optional("-v", "--verbose", "set flag for verbose EPR_PYMATCHING", GL_EPR_PYMATCHING_VERBOSE, false)

        // other:
        .optional("-m", "--mode", "0 = global, 1 = dual pass, -1 = single hardware EPR", eval_mode, 0)
    
        .parse(argc, argv);

    bool do_memory_experiment = (experiment_type == "memory");

    // Create EPR generation configuration
    gen::EPR_GEN_CONFIG circuit_config;
    circuit_config.attenuation_rate = attenuation_rate;
    circuit_config.photonic_link_error = photonic_link_error;
    circuit_config.hw1_round_ns = static_cast<uint64_t>(hw1_round_ns);
    circuit_config.hw2_round_ns = static_cast<uint64_t>(hw2_round_ns);
    circuit_config.phys_error = phys_error;

    // Generate the EPR-based surface code circuit
    auto gen_out = gen::sc_epr_generation(circuit_config, num_rounds, code_distance, do_memory_experiment);

    std::cout << "super rounds per round = " << gen_out.num_super_rounds 
            << ", hw1 rounds per super round = " << gen_out.num_hw1_rounds_per_super_round << "\n";

    // Write circuit to output file
    write_stim_circuit_to_file("generated.stim.out", gen_out.circuit);
    write_stim_circuit_to_file("first_pass.stim.out", gen_out.first_pass);
    write_stim_circuit_to_file("second_pass.stim.out", gen_out.second_pass);

    DECODER_EVAL_CONFIG eval_config{.stop_at_k_errors = static_cast<uint64_t>(num_errors)};
    DECODER_STATS stats;
    if (eval_mode == 0)
    {
        stats = eval_decoder<PYMATCHING>(gen_out.circuit, num_trials, eval_config, gen_out.circuit);
    }
    else if (eval_mode == -1)
    {
        stats = eval_decoder<PYMATCHING>(gen_out.second_pass, num_trials, eval_config, gen_out.second_pass);
    }
    else if (eval_mode == 1)
    {
        PYMATCHING reference_decoder(gen_out.circuit);
        EPR_PYMATCHING decoder(gen_out.circuit,
                                gen_out.first_pass,
                                gen_out.second_pass,
                                code_distance,
                                gen_out.num_super_rounds,
                                gen_out.num_hw1_rounds_per_super_round);
                                
        stats = benchmark_decoder(gen_out.circuit, decoder, num_trials,
                                [&reference_decoder] 
                                (syndrome_ref dets, syndrome_ref, syndrome_ref pred, std::ostream& debug_strm)
                                {
                                    std::vector<GRAPH_COMPONENT_ID> detectors;
                                    for (size_t i = 0; i < dets.num_bits_padded(); i++)
                                    {
                                        if (dets[i])
                                            detectors.push_back(static_cast<GRAPH_COMPONENT_ID>(i));
                                    }

                                    auto result = reference_decoder.decode(detectors, debug_strm);

                                    bool mismatch{false};
                                    debug_strm << "reference prediction:";
                                    for (size_t i = 0; i < result.flipped_observables.num_bits_padded(); i++)
                                    {
                                        if (result.flipped_observables[i])
                                            debug_strm << " " << i;
                                        
                                        mismatch |= (result.flipped_observables[i] != pred[i]);
                                    }
                                    debug_strm << "\n";

                                    return mismatch;
                                }, eval_config);
    }

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
