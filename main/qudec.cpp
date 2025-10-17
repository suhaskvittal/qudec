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
#include <iomanip>
#include <iostream>

#define USE_STIM_GENERATED_CIRCUITS

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

template <class T> void
print_stat(std::ostream& out, std::string name, T value)
{
    out << std::setw(64) << std::left << name;
    if constexpr (std::is_floating_point<T>::value)
        out << std::setw(12) << std::right << std::fixed << std::setprecision(8) << value;
    else
        out << std::setw(12) << std::right << value;
    out << "\n";
}

template <class INTA, class INTB> double
fpdiv(INTA a, INTB b)
{
    return static_cast<double>(a) / static_cast<double>(b);
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

    int64_t round_time;
    int64_t t1;
    int64_t t2;
    double e_g1q;
    double e_g2q;
    double e_readout;
    double e_idle;

    ARGPARSE()
        .optional("-f", "--stim-file", "stim file", stim_file, "")
        .optional("-d", "--code-distance", "code distance", code_distance, 3)
        .optional("-r", "--rounds", "number of rounds", num_rounds, 9)
        .optional("-t", "--trials", "number of trials to run", num_trials, 1'000'000)

        // circuit timing:
        .optional("-rt", "--round-time", "round time in ns", round_time, 1200)
        .optional("-t1", "--t1", "T1 time in ns", t1, 500'000)
        .optional("-t2", "--t2", "T2 time in ns", t2, 250'000)
        .optional("-e1", "--e-g1q", "gate error rate (1Q)", e_g1q, 1e-4)
        .optional("-e2", "--e-g2q", "gate error rate (2Q)", e_g2q, 1e-3)
        .optional("-em", "--e-readout", "readout error rate", e_readout, 1e-3)
        .optional("-ei", "--e-idle", "idle error rate", e_idle, 1e-4)
        .parse(argc, argv);

    stim::Circuit circuit;
    if (stim_file.empty())
    {
        const size_t sc_qubit_count = 2*code_distance*code_distance - 1;

#if defined(USE_STIM_GENERATED_CIRCUITS)
        stim::CircuitGenParameters params(num_rounds, code_distance, "rotated_memory_z");
        params.before_round_data_depolarization = e_g2q;
        params.after_reset_flip_probability = e_g1q;
        params.before_measure_flip_probability = e_readout;
        params.after_clifford_depolarization = e_g2q;
        auto gen_circuit = stim::generate_surface_code_circuit(params);
        circuit = gen_circuit.circuit;
#else
        gen::CIRCUIT_CONFIG conf = gen::CIRCUIT_CONFIG()
                                    .set_qubit_count(sc_qubit_count)
                                    .set_round_ns(round_time)
                                    .set_t1_ns(t1)
                                    .set_t2_ns(t2)
                                    .set_e_g1q(e_g1q)
                                    .set_e_g2q(e_g2q)
                                    .set_e_readout(e_readout)
                                    .set_e_idle(e_idle);
        circuit = gen::sc_memory(conf, num_rounds, code_distance); 

        std::cout << "GENERATED STIM CIRCUIT (d = " << code_distance << ", r = " << num_rounds << ") ===============================\n" 
                    << circuit 
                    << "\n===================================================================\n";
#endif

    }
    else
    {
        FILE* fin = fopen(stim_file.c_str(), "rb");
        circuit = stim::Circuit::from_file(fin);
        fclose(fin);
    }

    // initialize decoder:
//  BLOSSOM5 decoder(circuit);
    PYMATCHING decoder(circuit);

    auto stats = benchmark_decoder(circuit, decoder, num_trials);

    double ler = fpdiv(stats.errors, stats.trials);
    double mean_time_us = fpdiv(stats.total_time_us, stats.trials);
    double mean_time_us_nontrivial = fpdiv(stats.total_time_us, stats.trials - stats.trivial_trials);

    print_stat(std::cout, "LOGICAL_ERROR_RATE", ler);
    print_stat(std::cout, "MEAN_TIME_US", mean_time_us);
    print_stat(std::cout, "MEAN_TIME_US_NONTRIVIAL", mean_time_us_nontrivial);

    return 0;
}

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
