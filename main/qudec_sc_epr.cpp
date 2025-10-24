/*
 *  author: Suhas Vittal
 *  date:   23 October 2025
 *
 *  EPR-based surface code circuit generation using sc_epr_generation
 * */

#include "argparse.h"
#include "gen/epr.h"

#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>

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

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

int
main(int argc, char* argv[])
{
    int64_t     code_distance;
    int64_t     num_rounds;
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
        .optional("-m", "--memory-experiment", "do memory experiment (vs stability)", do_memory_experiment, true)
        
        // EPR-specific parameters:
        .optional("-a", "--attenuation-rate", "photon attenuation rate", attenuation_rate, 1e-2)
        .optional("-pl", "--photonic-link-error", "photonic link error rate", photonic_link_error, 1e-2)
        .optional("-hw1", "--hw1-round-ns", "HW1 round time in ns", hw1_round_ns, 1200)
        .optional("-hw2", "--hw2-round-ns", "HW2 round time in ns", hw2_round_ns, 1'200'000)
        .optional("-p", "--phys-error", "physical error rate", phys_error, 1e-3)
        
        // output file:
        .optional("-o", "--output", "output file for generated stim circuit", generated_stim_output_file, "generated.stim.out")
        .parse(argc, argv);

    // Create EPR generation configuration
    gen::EPR_GEN_CONFIG config;
    config.attenuation_rate = attenuation_rate;
    config.photonic_link_error = photonic_link_error;
    config.hw1_round_ns = static_cast<uint64_t>(hw1_round_ns);
    config.hw2_round_ns = static_cast<uint64_t>(hw2_round_ns);
    config.phys_error = phys_error;

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
    return 0;
}

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
