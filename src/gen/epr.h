/*
 *  author: Suhas Vittal
 *  date:   12 October 2025
 *
 *  Stim circuits for EPR pair generation
 * */

#ifndef GEN_EPR_h
#define GEN_EPR_h

#include "gen.h"

namespace gen
{

//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

struct EPR_GEN_CONFIG
{
    double attenuation_rate{1e-2};
    double photonic_link_error{1e-2};

    // two different configurations for the two hardware platforms
    // inter-connected via the link.
    CIRCUIT_CONFIG hw1_config;
    CIRCUIT_CONFIG hw2_config;
};

struct PINNED_QUBIT
{
    stim_qubit_type qubit;
    size_t          owner_id;
    stim_qubit_type pinned_check_qubit;
};

//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

stim::Circuit sc_epr_generation(const EPR_GEN_CONFIG&, size_t rounds, size_t distance);

//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

}   // namespace gen

#endif

