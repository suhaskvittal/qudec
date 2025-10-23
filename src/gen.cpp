/*
 *  author: Suhas Vittal
 * */

#include "gen.h"
#include "gen/scheduling.h"
#include "gen/utils.h"

#include <algorithm>
#include <numeric>

namespace gen
{

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

CIRCUIT_CONFIG&
CIRCUIT_CONFIG::set_qubit_count(size_t n)
{
    qubits = std::vector<QUBIT_INFO>(n, QUBIT_INFO{});
    for (size_t i = 0; i < n; ++i)
        couplings.push_back(std::vector<COUPLING_INFO>(n, COUPLING_INFO{}));
    return *this;
}

CIRCUIT_CONFIG&
CIRCUIT_CONFIG::set_round_ns(uint64_t ns)
{
    this->round_ns = ns;
    return *this;
}

CIRCUIT_CONFIG&
CIRCUIT_CONFIG::set_t1_ns(uint64_t ns)
{
    for (auto& qubit : qubits)
        qubit.t1_ns = ns;
    return *this;
}

CIRCUIT_CONFIG&
CIRCUIT_CONFIG::set_t2_ns(uint64_t ns)
{
    for (auto& qubit : qubits)
        qubit.t2_ns = ns;
    return *this;
}

CIRCUIT_CONFIG&
CIRCUIT_CONFIG::set_e_readout(double error_rate)
{
    for (auto& qubit : qubits)
        qubit.e_readout = error_rate;
    return *this;
}

CIRCUIT_CONFIG&
CIRCUIT_CONFIG::set_e_g1q(double error_rate)
{
    for (auto& qubit : qubits)
        qubit.e_g1q = error_rate;
    return *this;
}

CIRCUIT_CONFIG&
CIRCUIT_CONFIG::set_e_g2q(double error_rate)
{
    for (auto& coupling_row : couplings)
    {
        for (auto& coupling : coupling_row)
            coupling.e_g2q = error_rate;
    }
    return *this;
}

CIRCUIT_CONFIG&
CIRCUIT_CONFIG::set_e_idle(double error_rate)
{
    for (auto& qubit : qubits)
        qubit.e_idle = error_rate;
    return *this;
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

std::tuple<double,double,double>
pauli_twirling_approx(uint64_t t1_ns, uint64_t t2_ns, uint64_t round_ns)
{
    double ft1 = static_cast<double>(t1_ns),
           ft2 = static_cast<double>(t2_ns),
           fr  = static_cast<double>(round_ns);

    double x, y, z;

    x = 0.25 * (1 - exp(-fr/ft1));
    y = x;
    z = 0.25 * (1 - 2.0*exp(-fr/ft2) + exp(-fr/ft1));

    return std::make_tuple(x, y, z);
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

stim::Circuit
sc_memory(const CIRCUIT_CONFIG& config, size_t rounds, size_t distance, bool is_memory_x)
{
    SC_SCHEDULE_INFO sc(distance, distance);

    std::vector<stim_qubit_type> all_check_qubits(sc.x_check_qubits);
    all_check_qubits.insert(all_check_qubits.end(), sc.z_check_qubits.begin(), sc.z_check_qubits.end());
    
    const auto& det_qubits = is_memory_x ? sc.x_check_qubits : sc.z_check_qubits;

    stim::Circuit prolog, error_free_first_round, error_free_last_round, main_circuit, epilog;

    // prolog add coordinates:
    for (const auto& [q, coord] : sc.qubit_coords)
        prolog.safe_append_u("QUBIT_COORDS", {q}, {coord.first, coord.second});

    // initialize prolog -- error free
    util::init_data_qubits_in_basis(prolog, sc.data_qubits, is_memory_x);

    // error free first round:
    util::initialize_parity_qubits(error_free_first_round, sc.z_check_qubits, sc.x_check_qubits, config, false);
    util::do_cx_gates(error_free_first_round, sc.check_cx_order, sc.x_check_set, config, false);
    util::measure_parity_qubits(error_free_first_round, sc.z_check_qubits, sc.x_check_qubits, config, false);
    
    // main circuit:
    util::inject_timing_errors(main_circuit, sc.data_qubits, config);
    util::initialize_parity_qubits(main_circuit, sc.z_check_qubits, sc.x_check_qubits, config, true);
    util::do_cx_gates(main_circuit, sc.check_cx_order, sc.x_check_set, config, true);
    auto check_meas_order = util::measure_parity_qubits(main_circuit, sc.z_check_qubits, sc.x_check_qubits, config, true);

    util::create_detection_events(main_circuit, det_qubits, check_meas_order);
    main_circuit.safe_append_u("SHIFT_COORDS", {}, {0,1});

    // `error_free_last_round` is the same the first round but with detection events:
    error_free_last_round = error_free_first_round;
    util::create_detection_events(error_free_last_round, det_qubits, check_meas_order);

    // epilog
    const auto& obs = is_memory_x ? sc.x_obs : sc.z_obs;
    util::measure_data_qubits_and_observables(epilog, sc.data_qubits, {obs}, is_memory_x);

    // create final circuit:
    stim::Circuit fin;
    fin += prolog;
    fin += error_free_first_round;
    fin += main_circuit * (rounds);
    fin += error_free_last_round;
    fin += epilog;

    return fin;
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

stim::Circuit
sc_stability(const CIRCUIT_CONFIG& config, size_t rounds, size_t distance, bool is_boundary_x)
{
    // assert that `distance` is even
    if (distance % 2 != 0)
        throw std::invalid_argument("sc_stability: distance must be even");

    SC_SCHEDULE_INFO sc(distance, distance, is_boundary_x);
    std::vector<stim_qubit_type> all_check_qubits(sc.x_check_qubits);
    all_check_qubits.insert(all_check_qubits.end(), sc.z_check_qubits.begin(), sc.z_check_qubits.end());
    
    const auto& det_qubits = is_boundary_x ? sc.x_check_qubits : sc.z_check_qubits;

    stim::Circuit prolog, first_round, main_circuit, epilog;

    // prolog add coordinates:
    for (const auto& [q, coord] : sc.qubit_coords)
        prolog.safe_append_u("QUBIT_COORDS", {q}, {coord.first, coord.second});

    // initialize prolog -- error free (note we initialize in opposite basis of boundary)
    util::init_data_qubits_in_basis(prolog, sc.data_qubits, !is_boundary_x);

    // error free first round:
    util::inject_timing_errors(first_round, sc.data_qubits, config);
    util::initialize_parity_qubits(first_round, sc.z_check_qubits, sc.x_check_qubits, config, true);
    util::do_cx_gates(first_round, sc.check_cx_order, sc.x_check_set, config, true);
    auto check_meas_order = util::measure_parity_qubits(first_round, sc.z_check_qubits, sc.x_check_qubits, config, true);

    // main circuit:
    util::inject_timing_errors(main_circuit, sc.data_qubits, config);
    util::initialize_parity_qubits(main_circuit, sc.z_check_qubits, sc.x_check_qubits, config, true);
    util::do_cx_gates(main_circuit, sc.check_cx_order, sc.x_check_set, config, true);
    util::measure_parity_qubits(main_circuit, sc.z_check_qubits, sc.x_check_qubits, config, true);

    util::create_detection_events(main_circuit, det_qubits, check_meas_order);
    main_circuit.safe_append_u("SHIFT_COORDS", {}, {0,1});

    // epilog: observable is all check measurements of the same type as the boundary:
    if (!is_boundary_x)
    {
        epilog.safe_append_u("H", sc.data_qubits);
        for (auto q : sc.data_qubits)
            epilog.safe_append_ua("DEPOLARIZE1", {q}, config.qubits.at(q).e_g1q);
    }

    for (auto q : sc.data_qubits)
        epilog.safe_append_ua("M", {q}, config.qubits.at(q).e_readout);

    const size_t n_data_meas = sc.data_qubits.size();
    const size_t n_check_meas = check_meas_order.size();
    std::vector<uint32_t> obs_meas_id;
    for (auto q : det_qubits)
        obs_meas_id.push_back((n_data_meas+n_check_meas - check_meas_order.at(q)) | stim::TARGET_RECORD_BIT);

    epilog.safe_append_ua("OBSERVABLE_INCLUDE", obs_meas_id, 0);

    // create final circuit:
    stim::Circuit fin;
    fin += prolog;
    fin += first_round;
    fin += main_circuit * (rounds-1);
    fin += epilog;

    return fin;
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

}   // namespace gen
