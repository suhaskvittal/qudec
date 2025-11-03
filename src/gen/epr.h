/*
 *  author: Suhas Vittal
 *  date:   12 October 2025
 *
 *  Stim circuits for EPR pair generation
 * */

#ifndef GEN_EPR_h
#define GEN_EPR_h

#include "gen.h"
#include "gen/scheduling.h"
#include "gen/utils.h"

namespace gen
{

//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

struct EPR_GEN_CONFIG
{
    double attenuation_rate{1e-2};
    double photonic_link_error{1e-2};

    uint64_t hw1_round_ns{1200};
    uint64_t hw2_round_ns{1'200'000};
    double   phys_error{1e-3};
};

//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

struct SC_EPR_SCHEDULE_INFO
{
    using qubit_coord_type = SC_SCHEDULE_INFO::qubit_coord_type;

    SC_SCHEDULE_INFO sc;

    std::unordered_map<stim_qubit_type, stim_qubit_type>              epr_checks;
    std::unordered_map<stim_qubit_type, std::vector<stim_qubit_type>> epr_cx_order;

    // use this for crumble instead of `sc.qubit_coords`
    std::unordered_map<stim_qubit_type, qubit_coord_type> qubit_coords;

    std::unordered_set<stim_qubit_type> hw1_qubit_set;
    std::unordered_set<stim_qubit_type> hw2_qubit_set;

    SC_EPR_SCHEDULE_INFO(size_t d, bool dual=false);
};

//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

struct SC_EPR_GEN_OUTPUT
{
    stim::Circuit circuit;
    stim::Circuit first_pass;
    stim::Circuit second_pass;

    size_t num_super_rounds;
    size_t num_hw1_rounds_per_super_round;
};

SC_EPR_GEN_OUTPUT sc_epr_generation(const EPR_GEN_CONFIG&, size_t rounds, size_t distance, bool do_memory_experiment);

//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

util::check_meas_map sc_epr_create_super_round(stim::Circuit&, 
                                                const SC_EPR_SCHEDULE_INFO&,
                                                uint64_t hw1_round_ns,
                                                uint64_t hw2_round_ns,
                                                double t1_ns_hw1,
                                                double t2_ns_hw1,
                                                double t1_ns_hw2,
                                                double t2_ns_hw2,
                                                double e_readout,
                                                double e_g1q,
                                                double e_g2q,
                                                double e_idle,
                                                double e_photonic_link,
                                                double hw1_error_scale_factor=1.0);

util::check_meas_map sc_epr_create_hw1_only_circuit(stim::Circuit&, 
                                                        const SC_EPR_SCHEDULE_INFO&,
                                                        uint64_t hw1_round_ns,
                                                        double t1_ns_hw1,
                                                        double t2_ns_hw1,
                                                        double e_readout,
                                                        double e_g1q,
                                                        double e_g2q,
                                                        double e_idle);

//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

/*
 * These are specialized functions for creating detection events for the EPR generation circuit.
 * */

void sc_epr_create_detection_events_super_round(stim::Circuit&,
                                                    const util::stim_qubit_array&,
                                                    const util::check_meas_map& cm_super_round,
                                                    const util::check_meas_map& cm_hw1_only_round,
                                                    size_t hw1_rounds_per_super_round,
                                                    bool is_first_round,
                                                    const SC_EPR_SCHEDULE_INFO&);

void sc_epr_create_detection_events_adjacent_hw1_rounds(stim::Circuit&, 
                                                    const util::stim_qubit_array&,
                                                    const util::check_meas_map& cm_this_round,
                                                    const util::check_meas_map& cm_prev_round,
                                                    const SC_EPR_SCHEDULE_INFO&);

//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

/*
 * `sc_epr_create_detection_events_generic` is a generic function that only creates detection events
 * for checks in `cm` and `checks`. The previous round is assumed to measure the same checks.
 *
 * This is considered "generic" as this is usually how detection events are created (outside of this
 * EPR generation circuit).
 * */

void sc_epr_create_detection_events_generic(stim::Circuit&, 
                                                const util::stim_qubit_array&,
                                                const util::check_meas_map& cm,
                                                bool is_first_round,
                                                const SC_EPR_SCHEDULE_INFO&);

//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

}   // namespace gen

#endif

