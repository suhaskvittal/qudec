/*
 *  author: Suhas Vittal
 *  date:   12 October 2025
 * */

#ifndef GEN_UTILS_h
#define GEN_UTILS_h

#include "gen.h"

#include <stim/circuit/circuit.h>

namespace gen
{
namespace util
{

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

using stim_qubit_array = std::vector<stim_qubit_type>;
using check_set_type = std::unordered_set<stim_qubit_type>;
using check_impl_map = std::unordered_map<stim_qubit_type, std::vector<stim_qubit_type>>;
using check_meas_map = std::unordered_map<stim_qubit_type, size_t>;

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

void init_data_qubits_in_basis(stim::Circuit&, stim_qubit_array, bool init_in_hadamard_basis);
void inject_timing_errors(stim::Circuit&, stim_qubit_array, const CIRCUIT_CONFIG&);
void initialize_parity_qubits(stim::Circuit&, stim_qubit_array z_checks, stim_qubit_array x_checks, const CIRCUIT_CONFIG&, bool inject_errors);

void do_cx_gates(stim::Circuit&, 
                    const check_impl_map&,
                    check_set_type,
                    const CIRCUIT_CONFIG&,
                    bool inject_errors,
                    size_t max_steps=4);

check_meas_map measure_parity_qubits(stim::Circuit&, 
                                        stim_qubit_array z_checks,
                                        stim_qubit_array x_checks,
                                        const CIRCUIT_CONFIG&,
                                        bool inject_errors);

void create_detection_events(stim::Circuit&, stim_qubit_array, const check_meas_map&);
void measure_data_qubits_and_observables(stim::Circuit&, stim_qubit_array, std::vector<stim_qubit_array> observables, bool measure_in_hadamard_basis);

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

void validate_check_cx_order(const check_impl_map&, size_t max_steps=4);

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

}   // namespace util
}   // namespace gen

#endif
