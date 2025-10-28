/*
 *  author: Suhas Vittal
 *  date:   27 October 2025
 *
 *  Sliding window decoder for surface codes
 * */

#include "argparse.h"
#include "decoder_eval.h"
#include "decoder/surface_code.h"
#include "gen.h"
#include "qudec_common.h"

#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>

// Global debug configuration variable definition
bool GL_DEBUG_DECODER{false};

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

int64_t
_safe_compute_detectors_per_round(int64_t total_detectors, int64_t rounds)
{
    if (total_detectors % rounds != 0)
    {
        std::cerr << "total_detectors = " << total_detectors << ", rounds = " << rounds << "\n";
        throw std::runtime_error("total_detectors % rounds != 0");
    }

    return total_detectors / rounds;
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
    int64_t     commit_size;
    
    double phys_error;
    int64_t round_time;
    int64_t t1;
    int64_t t2;
    double e_g1q;
    double e_g2q;
    double e_readout;
    double e_idle;

    std::string experiment;
    std::string generated_stim_output_file;
    std::string decoder_circuit_output_file;

    ARGPARSE()
        .optional("-d", "--code-distance", "code distance", code_distance, 3)
        .optional("-r", "--rounds", "number of rounds", num_rounds, 9)
        .optional("-t", "--trials", "number of trials to run", num_trials, 1'000'000)
        .optional("-c", "--commit-size", "commit size (-1 defaults to `code_distance`)", commit_size, -1)
        .optional("-k", "--stop-after-errors", "stop after this many errors", num_errors, 10)

        // circuit timing:
        .optional("-p", "--phys-error", "physical error rate", phys_error, 1e-3)
        .optional("-rt", "--round-time", "round time in ns", round_time, 1200)
        .optional("-t1", "--t1", "T1 time in us", t1, 1000)
        .optional("-t2", "--t2", "T2 time in us", t2, 500)
        .optional("-e1", "--e-g1q", "gate error rate (1Q)", e_g1q, 1e-4)
        .optional("-e2", "--e-g2q", "gate error rate (2Q)", e_g2q, 1e-3)
        .optional("-em", "--e-readout", "readout error rate", e_readout, 3e-3)
        .optional("-ei", "--e-idle", "idle error rate", e_idle, 1e-4)

        // experiment name:
        .optional("", "--experiment", "experiment name", experiment, "sc_memory_z")

        // decoder:
        .optional("-dd", "--debug-decoder", "enable decoder debug output", GL_DEBUG_DECODER, false)
        .parse(argc, argv);

    if (commit_size < 0)
        commit_size = code_distance;

    // update error and timing based on value of p
    double scale_factor = phys_error / 1e-3;
    t1 *= 1.0/scale_factor;
    t2 *= 1.0/scale_factor;
    e_g1q *= scale_factor;
    e_g2q *= scale_factor;
    e_readout *= scale_factor;
    e_idle *= scale_factor;

    // Generate circuit configuration
    size_t qubit_count;
    if (experiment == "sc_memory_x" || experiment == "sc_memory_z")
        qubit_count = gen::sc_memory_get_qubit_count(code_distance);
    else
        throw std::runtime_error("invalid experiment: " + experiment);

    gen::CIRCUIT_CONFIG circuit_conf = gen::CIRCUIT_CONFIG()
                                .set_qubit_count(qubit_count)
                                .set_round_ns(round_time)
                                .set_t1_ns(t1*1000)
                                .set_t2_ns(t2*1000)
                                .set_e_g1q(e_g1q)
                                .set_e_g2q(e_g2q)
                                .set_e_readout(e_readout)
                                .set_e_idle(e_idle);

    // Generate the full circuit (r rounds) for syndrome generation
    const int64_t window_size{2*commit_size};
    int64_t detectors_per_round;

    stim::Circuit full_circuit, decoder_circuit;
    if (experiment == "sc_memory_x" || experiment == "sc_memory_z")
    {
        full_circuit = gen::sc_memory(circuit_conf, num_rounds, code_distance, experiment == "sc_memory_x");
        decoder_circuit = gen::sc_memory(circuit_conf, window_size+1, code_distance, experiment == "sc_memory_x");

        detectors_per_round = _safe_compute_detectors_per_round(decoder_circuit.count_detectors(), window_size+2);
    }

    PYMATCHING reference_decoder(full_circuit);
    SLIDING_PYMATCHING decoder(decoder_circuit, commit_size, window_size, detectors_per_round, num_rounds);

    DECODER_EVAL_CONFIG eval_conf
    {
        .batch_size = 8192,
        .enable_clock = true,
        .seed = 0,
        .stop_at_k_errors = num_errors
    };
    auto stats = benchmark_decoder(full_circuit, decoder, num_trials,
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
                                }, eval_conf);

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
