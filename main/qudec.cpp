/*
 *  author: Suhas Vittal
 *  date:   12 October 2025
 * */

#include "argparse.h"
#include "decoder/surface_code.h"
#include "decoder_eval.h"

#include <iomanip>
#include <iostream>

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
    int64_t     num_trials;

    ARGPARSE()
        .required("stim-file", "stim file containing circuit to decode", stim_file)
        .optional("-t", "--trials", "number of trials to run", num_trials, 1'000'000)
        .parse(argc, argv);

    stim::Circuit circuit(stim_file);

    // initialize decoder:
    BLOSSOM5 decoder(circuit);

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
