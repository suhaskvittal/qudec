/*
 *  author: Suhas Vittal
 *  date:   12 October 2025
 * */

#include "gen/utils.h"

#include <algorithm>
#include <iomanip>
#include <iostream>

namespace gen
{
namespace util
{

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

void
init_data_qubits_in_basis(stim::Circuit& circuit, std::vector<stim_qubit_type> qubits, bool init_in_hadamard_basis)
{
    circuit.safe_append_u("R", qubits);
    if (init_in_hadamard_basis)
        circuit.safe_append_u("H", qubits);
}

void
inject_timing_errors(stim::Circuit& circuit, std::vector<stim_qubit_type> qubits, const CIRCUIT_CONFIG& config)
{
    for (stim_qubit_type q : qubits)
    {
        const auto& q_info = config.qubits.at(q);
        auto [ex,ey,ez] = pauli_twirling_approx(q_info.t1_ns, q_info.t2_ns, config.round_ns);
        circuit.safe_append_u("PAULI_CHANNEL_1", {q}, {ex,ey,ez});
    }
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

void
initialize_parity_qubits(stim::Circuit& circuit, 
                        stim_qubit_array z_checks,
                        stim_qubit_array x_checks,
                        const CIRCUIT_CONFIG& config,
                        bool inject_errors)
{
    stim_qubit_array all_checks(z_checks);
    all_checks.insert(all_checks.end(), x_checks.begin(), x_checks.end());

    circuit.safe_append_u("R", all_checks);
    if (inject_errors)
    {
        for (auto q : all_checks)
            circuit.safe_append_ua("X_ERROR", {q}, config.qubits.at(q).e_g1q);
    }

    circuit.safe_append_u("H", x_checks);
    if (inject_errors)
    {
        for (auto q : x_checks)
            circuit.safe_append_ua("DEPOLARIZE1", {q}, config.qubits.at(q).e_g1q);
    }
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

void
do_cx_gates(stim::Circuit& circuit, 
            const check_impl_map& check_cx_order,
            check_set_type x_check_set,
            const CIRCUIT_CONFIG& config,
            bool inject_errors,
            size_t max_steps)
{
    for (size_t t = 0; t < max_steps; t++)
    {
        stim_qubit_array targets;
        for (const auto& [c, cx_order] : check_cx_order)
        {
            stim_qubit_type q1{c},
                            q2{cx_order.at(t)};
            if (q2 == NO_QUBIT)
                continue;

            if (!x_check_set.count(c))
                std::swap(q1,q2);
            targets.push_back(q1);
            targets.push_back(q2);
        }

        circuit.safe_append_u("CX", targets);

        if (inject_errors)
        {
            // first do CX errors:
            for (size_t i = 0; i < targets.size(); i += 2)
            {
                stim_qubit_type q1 = targets[i],
                                q2 = targets[i+1];
                const auto& cpl = config.couplings.at(q1).at(q2);
                circuit.safe_append_ua("DEPOLARIZE2", {q1,q2}, cpl.e_g2q);
            }

            // now handle idling errors:
            std::unordered_set<stim_qubit_type> target_set(targets.begin(), targets.end());
            for (stim_qubit_type q = 0; q < config.qubits.size(); q++)
            {
                if (!target_set.count(q))
                    circuit.safe_append_ua("DEPOLARIZE1", {q}, config.qubits.at(q).e_idle);
            }
        }

        circuit.safe_append_u("TICK", {});
    }
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

check_meas_map
measure_parity_qubits(stim::Circuit& circuit,
                        stim_qubit_array z_checks,
                        stim_qubit_array x_checks,
                        const CIRCUIT_CONFIG& config,
                        bool inject_errors)
{
    stim_qubit_array all_checks(z_checks);
    all_checks.insert(all_checks.end(), x_checks.begin(), x_checks.end());

    circuit.safe_append_u("H", x_checks);

    if (inject_errors)
    {
        // inject H errors
        for (auto q : x_checks)
            circuit.safe_append_ua("DEPOLARIZE1", {q}, config.qubits.at(q).e_g1q);
    }

    for (auto q : all_checks)
    {
        if (inject_errors)
            circuit.safe_append_ua("M", {q}, config.qubits.at(q).e_readout);
        else
            circuit.safe_append_u("M", {q});
    }

    check_meas_map check_meas_order;
    for (size_t i = 0; i < all_checks.size(); i++)
        check_meas_order[all_checks[i]] = i;
    return check_meas_order;
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

void
create_detection_events(stim::Circuit& circuit, stim_qubit_array checks, const check_meas_map& check_meas_order)
{
    const uint32_t n_check_meas = check_meas_order.size();
    size_t detector_id{0};
    for (auto q : checks)
    {
        uint32_t meas_idx = check_meas_order.at(q);
        uint32_t base_meas_id = (n_check_meas-meas_idx) | stim::TARGET_RECORD_BIT,
                 prev_meas_id = (2*n_check_meas-meas_idx) | stim::TARGET_RECORD_BIT;

        circuit.safe_append_u("DETECTOR", {prev_meas_id, base_meas_id}, {detector_id, 0});
        detector_id++;
    }
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

void
measure_data_qubits_and_observables(stim::Circuit& circuit, 
                                    stim_qubit_array data_qubits,
                                    std::vector<stim_qubit_array> observables,
                                    bool measure_in_hadamard_basis)
{
    if (measure_in_hadamard_basis)
        circuit.safe_append_u("H", data_qubits);
    circuit.safe_append_u("M", data_qubits);

    std::unordered_map<stim_qubit_type, uint32_t> dq_meas_order;
    for (size_t i = 0; i < data_qubits.size(); i++)
        dq_meas_order[data_qubits.at(i)] = i;

    const uint32_t n_data_meas = data_qubits.size();
    for (size_t i = 0; i < observables.size(); i++)
    {
        const auto& obs = observables.at(i);

        std::vector<uint32_t> obs_meas_id;
        std::transform(obs.begin(), obs.end(), std::back_inserter(obs_meas_id),
                        [&dq_meas_order, n_data_meas] (stim_qubit_type q) 
                        {
                            return (n_data_meas - dq_meas_order.at(q)) | stim::TARGET_RECORD_BIT;
                        });
        circuit.safe_append_ua("OBSERVABLE_INCLUDE", obs_meas_id, i);
    }
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

void
validate_check_cx_order(const check_impl_map& check_cx_order, size_t max_steps)
{
    for (size_t t = 0; t < max_steps; t++)
    {
        std::unordered_set<stim_qubit_type> qubits_in_this_timestep;
        for (const auto& [__unused_check, cx_order] : check_cx_order)
        {
            stim_qubit_type q = cx_order.at(t);
            if (q == NO_QUBIT)
                continue;
            if (qubits_in_this_timestep.count(q))
            {
                // dump schedule:
                for (const auto& [c, cx_order] : check_cx_order)
                {
                    std::cerr << "check " << c << ": ";
                    for (auto q : cx_order)
                    {
                        if (q == NO_QUBIT) std::cerr << std::setw(4) << " ";
                        else               std::cerr << std::setw(4) << q;
                    }
                    std::cerr << "\n";
                }

                throw std::runtime_error("qubit " + std::to_string(q) + " is scheduled at the same timestep twice");
            }
            qubits_in_this_timestep.insert(q);
        }
    }
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

}    // namespace util
}    // namespace gen
