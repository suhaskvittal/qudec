/*
 *  author: Suhas Vittal
 * */

#ifndef GEN_h
#define GEN_h

#include <cstdint>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <stim/circuit/circuit.h>

namespace gen
{

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

using stim_qubit_type = uint32_t;
constexpr stim_qubit_type NO_QUBIT{std::numeric_limits<stim_qubit_type>::max()};

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

struct CIRCUIT_CONFIG
{
    struct QUBIT_INFO
    {
        uint64_t t1_ns{1'000'000};
        uint64_t t2_ns{500'000};
        double   e_readout{1e-3};
        double   e_g1q{1e-4};
        double   e_idle{1e-4};
    };

    struct COUPLING_INFO
    {
        double e_g2q{1e-3};
    };

    uint64_t round_ns{1200};

    std::vector<QUBIT_INFO> qubits;
    std::vector<std::vector<COUPLING_INFO>> couplings;

    CIRCUIT_CONFIG() =default;

    /*
     * For ease of use, we allow config setup to be down via the functions below (easy to read):
     * This is for making simple circuits with uniform noise characteristics.
     * */

    CIRCUIT_CONFIG& set_qubit_count(size_t);
    CIRCUIT_CONFIG& set_round_ns(uint64_t);
    CIRCUIT_CONFIG& set_t1_ns(uint64_t);
    CIRCUIT_CONFIG& set_t2_ns(uint64_t);
    CIRCUIT_CONFIG& set_e_readout(double);
    CIRCUIT_CONFIG& set_e_g1q(double);
    CIRCUIT_CONFIG& set_e_g2q(double);
    CIRCUIT_CONFIG& set_e_idle(double);
};

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

// returns X, Y, Z error probabilities given T1, T2, and round latency (all in ns)
std::tuple<double,double,double> pauli_twirling_approx(uint64_t t1_ns, uint64_t t2_ns, uint64_t round_ns);

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

constexpr size_t sc_memory_get_qubit_count(size_t d) { return 2*d*d-1; }
constexpr size_t sc_stability_get_qubit_count(size_t d) { return d*d + (d-1)*(d-1) + 2*d; }


stim::Circuit sc_memory(const CIRCUIT_CONFIG&, size_t rounds, size_t distance, bool is_memory_x=false);
stim::Circuit sc_stability(const CIRCUIT_CONFIG&, size_t rounds, size_t distance, bool is_boundary_x=false);

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

//#include "gen/epr.h"

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

}   // namespace gen

#endif
