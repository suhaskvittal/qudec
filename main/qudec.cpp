/*
 *  author: Suhas Vittal
 *  date:   12 October 2025
 * */

#include "argparse.h"
#include "decoder/surface_code.h"
#include "decoder_eval.h"
#include "gen.h"

#include <stim/gen/gen_surface_code.h>

#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>

//#define USE_STIM_GENERATED_CIRCUITS

// Global debug configuration variable definition
bool GL_DEBUG_DECODER = false;

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

template <class T> void
print_stat(std::ostream& out, std::string name, T value)
{
    out << std::setw(64) << std::left << name;
    if constexpr (std::is_floating_point<T>::value)
    {
        if (value < 1e-3)
            out << std::setw(12) << std::right << std::scientific << std::setprecision(4) << value;
        else
            out << std::setw(12) << std::right << std::fixed << std::setprecision(8) << value;
    }
    else
    {
        out << std::setw(12) << std::right << value;
    }
    out << "\n";
}

template <class INTA, class INTB> double
fpdiv(INTA a, INTB b)
{
    return static_cast<double>(a) / static_cast<double>(b);
}

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

template <class IMPL> void
print_decoder_stats(std::ostream& out, const IMPL& dec)
{
}

template <class IMPL, class... IMPL_ARGS> DECODER_STATS
eval_decoder(const stim::Circuit& circuit, uint64_t num_trials, IMPL_ARGS... args)
{
    IMPL decoder(circuit, std::forward<IMPL_ARGS>(args)...);
    auto out = benchmark_decoder(circuit, decoder, num_trials);

    print_decoder_stats(std::cout, decoder);

    return out;
}

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

int
main(int argc, char* argv[])
{
    std::string stim_file;
    int64_t     code_distance;
    int64_t     num_rounds;
    int64_t     num_trials;

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
    std::string decoder;

    ARGPARSE()
        .optional("-f", "--stim-file", "stim file", stim_file, "")
        .optional("-d", "--code-distance", "code distance", code_distance, 3)
        .optional("-r", "--rounds", "number of rounds", num_rounds, 9)
        .optional("-t", "--trials", "number of trials to run", num_trials, 1'000'000)

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
        .optional("", "--experiment", "experiment name -- do not set if stim-file is used", experiment, "sc_memory_z")
        .optional("", "--generated-stim-output-file", "output file for generated stim circuit", generated_stim_output_file, "generated.stim.out")

        // decoder:
        .optional("", "--decoder", "decoder to use", decoder, "pymatching")
        .optional("-dd", "--debug-decoder", "enable decoder debug output", GL_DEBUG_DECODER, false)
        .parse(argc, argv);

    // update error and timing based on value of p
    double scale_factor = phys_error / 1e-3;
    t1 *= 1.0/scale_factor;
    t2 *= 1.0/scale_factor;
    e_g1q *= scale_factor;
    e_g2q *= scale_factor;
    e_readout *= scale_factor;
    e_idle *= scale_factor;

    stim::Circuit circuit;
    if (stim_file.empty())
    {
        size_t qubit_count;
        if (experiment == "sc_memory_x" || experiment == "sc_memory_z")
            qubit_count = gen::sc_memory_get_qubit_count(code_distance);
        else if (experiment == "sc_stability_x" || experiment == "sc_stability_z")
            qubit_count = gen::sc_stability_get_qubit_count(code_distance);
        else
            throw std::runtime_error("invalid experiment: " + experiment);

        gen::CIRCUIT_CONFIG conf = gen::CIRCUIT_CONFIG()
                                    .set_qubit_count(qubit_count)
                                    .set_round_ns(round_time)
                                    .set_t1_ns(t1*1000)
                                    .set_t2_ns(t2*1000)
                                    .set_e_g1q(e_g1q)
                                    .set_e_g2q(e_g2q)
                                    .set_e_readout(e_readout)
                                    .set_e_idle(e_idle);

        if (experiment == "sc_memory_x" || experiment == "sc_memory_z")
            circuit = gen::sc_memory(conf, num_rounds, code_distance, experiment == "sc_memory_x"); 
        else if (experiment == "sc_stability_x" || experiment == "sc_stability_z")
            circuit = gen::sc_stability(conf, num_rounds, code_distance, experiment == "sc_stability_x");

        if (code_distance <= 3)
        {
            std::cout << "======================== GENERATED CIRCUIT ==========================\n";
            std::cout << circuit << "\n";
            std::cout << "=====================================================================\n";
        }

        // spit out circuit to output file:
        std::ofstream out(generated_stim_output_file);
        out << circuit.str() << "\n";
        out.close();
    }
    else
    {
        FILE* fin = fopen(stim_file.c_str(), "rb");
        circuit = stim::Circuit::from_file(fin);
        fclose(fin);
    }

    FPD_CONFIG fpd_conf;
    fpd_conf.cache_chain_limit = code_distance >> 2;

    DECODER_STATS stats;
    if (decoder == "pymatching")
        stats = eval_decoder<PYMATCHING>(circuit, num_trials);
    else if (decoder == "blossom5")
        stats = eval_decoder<BLOSSOM5>(circuit, num_trials);
    else
        throw std::runtime_error("invalid decoder: " + decoder);

    double ler = fpdiv(stats.errors, stats.trials);
    double mean_time_us = fpdiv(stats.total_time_us, stats.trials);
    double mean_time_us_nontrivial = fpdiv(stats.total_time_us, stats.trials - stats.trivial_trials);

    print_stat(std::cout, "LOGICAL_ERRORS", stats.errors);
    print_stat(std::cout, "TRIALS", stats.trials);
    print_stat(std::cout, "LOGICAL_ERROR_RATE", ler);
    print_stat(std::cout, "MEAN_TIME_US", mean_time_us);
    print_stat(std::cout, "MEAN_TIME_US_NONTRIVIAL", mean_time_us_nontrivial);

    return 0;
}

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
