/*
 *  author: Suhas Vittal
 * */

#include "gen.h"
#include "stim/circuit/circuit.h"

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

std::vector<stim_qubit_type>
get_cx_order(stim_qubit_type nw, stim_qubit_type ne, stim_qubit_type sw, stim_qubit_type se, bool is_x_check)
{
    if (is_x_check) 
        return {nw, ne, sw, se};
    else
        return {nw, sw, ne, se};
}

SC_SCHEDULE_INFO::SC_SCHEDULE_INFO(size_t distance)
    :data_qubits(distance*distance),
    x_obs(distance),
    z_obs(distance)
{
    // initialize data structures describing syndrome extraction circuit:
    std::iota(data_qubits.begin(), data_qubits.end(), 0);

    x_check_qubits.reserve((distance*distance-1)/2);
    z_check_qubits.reserve((distance*distance-1)/2);

    // start with bulk checks first:
    stim_qubit_type check_qubit{distance*distance};
    for (size_t r = 0; r < distance-1; r++)
    {
        for (size_t c = 0; c < distance-1; c++)
        {
            bool is_x_check = ((r+c) & 1) == 0;

            stim_qubit_type nw = r*distance + c;
            stim_qubit_type ne = r*distance + c+1;
            stim_qubit_type sw = (r+1)*distance + c;
            stim_qubit_type se = (r+1)*distance + c+1;

            std::vector<stim_qubit_type> cx_order = get_cx_order(nw, ne, sw, se, is_x_check);
            if (is_x_check)
                x_check_qubits.push_back(check_qubit);
            else
                z_check_qubits.push_back(check_qubit);

            check_cx_order[check_qubit] = cx_order;
            check_qubit++;
        }
    }

    // now handle boundary checks:
    for (size_t i = 0; i < distance-1; i++)
    {
        std::vector<stim_qubit_type> z_cx_order, x_cx_order;
        if (i & 1)  // z check on right boundary, x check on upper boundary
        {
            stim_qubit_type zq1 = (i+1)*distance - 1,
                            zq2 = (i+2)*distance - 1,
                            xq1 = i,
                            xq2 = i+1;
            z_cx_order = get_cx_order(zq1, NO_QUBIT, zq2, NO_QUBIT, false);
            x_cx_order = get_cx_order(NO_QUBIT, NO_QUBIT, xq1, xq2, true);
        }
        else  // z check on left boundary, x check on lower boundary
        {
            stim_qubit_type zq1 = i*distance,
                            zq2 = (i+1)*distance,
                            xq1 = (distance-1)*distance + i,
                            xq2 = (distance-1)*distance + i+1;

            z_cx_order = get_cx_order(NO_QUBIT, zq1, NO_QUBIT, zq2, false);
            x_cx_order = get_cx_order(xq1, xq2, NO_QUBIT, NO_QUBIT, true);
        }

        z_check_qubits.push_back(check_qubit);
        check_cx_order[check_qubit] = z_cx_order;
        check_qubit++;

        x_check_qubits.push_back(check_qubit);
        check_cx_order[check_qubit] = x_cx_order;
        check_qubit++;
    }

    x_check_set = std::unordered_set<stim_qubit_type>(x_check_qubits.begin(), x_check_qubits.end());

    for (size_t i = 0; i < distance; i++)
    {
        x_obs[i] = i*distance;  // leftmost column
        z_obs.at(i) = (distance-1)*distance + i;  // bottom row
    }

#if defined(GEN_VERIFY_CX_SCHEDULE)
    // make sure no qubit is scheduled at the same timestep:
    for (size_t t = 0; t < 4; t++)
    {
        std::unordered_set<stim_qubit_type> qubits_in_this_timestep;
        for (const auto& [__unused_check, cx_order] : check_cx_order)
        {
            stim_qubit_type q = cx_order.at(t);
            if (q == NO_QUBIT)
                continue;
            if (qubits_in_this_timestep.count(q))
                throw std::runtime_error("qubit " + std::to_string(q) + " is scheduled at the same timestep twice");
            qubits_in_this_timestep.insert(q);
        }
    }
#endif
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

stim::Circuit
sc_memory(const CIRCUIT_CONFIG& config, size_t rounds, size_t distance, bool is_memory_x)
{
    SC_SCHEDULE_INFO schedule_info(distance);

    std::vector<stim_qubit_type> all_check_qubits(schedule_info.x_check_qubits);
    all_check_qubits.insert(all_check_qubits.end(), 
                                schedule_info.z_check_qubits.begin(),
                                schedule_info.z_check_qubits.end());

    stim::Circuit prolog, error_free_round_circuit, round_circuit, epilog;

    ///////////////////////
    /// PROLOG          ///
    ///////////////////////

    // initialize prolog -- error free
    prolog.safe_append_u("R", schedule_info.data_qubits);
    if (is_memory_x)
        prolog.safe_append_u("H", schedule_info.data_qubits);

    ///////////////////////
    /// SYNDROME EXT    ///
    ///////////////////////

    // create syndrome extraction round (will also concurrently add instructions to `prolog` for
    // first error-free round that intiializes the code state, but `round_circuit` has errors)
    
    // inject decoherence/dephasing
    for (auto q : schedule_info.data_qubits)
    {
        const auto& q_info = config.qubits.at(q);
        auto [ex,ey,ez] = pauli_twirling_approx(q_info.t1_ns, q_info.t2_ns, config.round_ns);
        if (ex > 0 || ey > 0 || ez > 0)
            round_circuit.safe_append_u("PAULI_CHANNEL_1", {q}, {ex,ey,ez});
    }
    
    // reset parity qubits
    error_free_round_circuit.safe_append_u("R", all_check_qubits);
    round_circuit.safe_append_u("R", all_check_qubits);
    for (auto q : all_check_qubits)
    {
        if (config.qubits.at(q).e_g1q > 0)
            round_circuit.safe_append_ua("X_ERROR", {q}, config.qubits.at(q).e_g1q);
    }

    // initialize parity qubits
    error_free_round_circuit.safe_append_u("H", schedule_info.x_check_qubits);
    round_circuit.safe_append_u("H", schedule_info.x_check_qubits);
    for (auto q : schedule_info.x_check_qubits)
    {
        if (config.qubits.at(q).e_g1q > 0)
            round_circuit.safe_append_ua("DEPOLARIZE1", {q}, config.qubits.at(q).e_g1q);
    }

    // do cnot gates:
    for (size_t t = 0; t < 4; t++)
    {
        std::vector<stim_qubit_type> targets;
        for (const auto& [check_qubit, cx_order] : schedule_info.check_cx_order)
        {
            stim_qubit_type q = cx_order.at(t);
            if (q == SC_SCHEDULE_INFO::NO_QUBIT)
                continue;
            bool is_x_check = schedule_info.x_check_set.count(check_qubit);
            if (is_x_check)
            {
                targets.push_back(check_qubit);
                targets.push_back(q);
            }
            else
            {
                targets.push_back(q);
                targets.push_back(check_qubit);
            }
        }

        error_free_round_circuit.safe_append_u("CX", targets);
        round_circuit.safe_append_u("CX", targets);

        // inject CX error
        for (size_t i = 0; i < targets.size(); i += 2)
        {
            stim_qubit_type q1 = targets[i],
                            q2 = targets[i+1];
            std::vector<stim_qubit_type> args{q1, q2};

            const auto& cpl = config.couplings.at(q1).at(q2);
            if (cpl.e_g2q > 0)
                round_circuit.safe_append_ua("DEPOLARIZE2", {q1,q2}, cpl.e_g2q);
        }

        // inject idling error on any qubits that did not do a CX
        std::unordered_set<stim_qubit_type> cx_targets(targets.begin(), targets.end());
        std::vector<stim_qubit_type> idle_qubits;
        std::copy_if(schedule_info.data_qubits.begin(), schedule_info.data_qubits.end(), std::back_inserter(idle_qubits),
                    [&cx_targets] (stim_qubit_type q) { return !cx_targets.count(q); });
        std::copy_if(all_check_qubits.begin(), all_check_qubits.end(), std::back_inserter(idle_qubits),
                    [&cx_targets] (stim_qubit_type q) { return !cx_targets.count(q); });

        for (auto q : idle_qubits)
        {
            if (config.qubits.at(q).e_idle > 0)
                round_circuit.safe_append_ua("DEPOLARIZE1", {q}, config.qubits.at(q).e_idle);
        }

        error_free_round_circuit.safe_append_u("TICK", {});
        round_circuit.safe_append_u("TICK", {});
    }

    // measure parity qubits
    error_free_round_circuit.safe_append_u("H", schedule_info.x_check_qubits);
    round_circuit.safe_append_u("H", schedule_info.x_check_qubits);
    for (auto q : schedule_info.x_check_qubits)
    {
        if (config.qubits.at(q).e_g1q > 0)
            round_circuit.safe_append_ua("DEPOLARIZE1", {q}, config.qubits.at(q).e_g1q);
    }

    for (auto q : all_check_qubits)
    {
        if (config.qubits.at(q).e_readout > 0)
            round_circuit.safe_append_ua("X_ERROR", {q}, config.qubits.at(q).e_readout);
    }
    error_free_round_circuit.safe_append_u("M", all_check_qubits);
    round_circuit.safe_append_u("M", all_check_qubits);

    // `round_circuit`: initialize detection events
    const size_t n_check_meas = all_check_qubits.size();
    size_t detector_id{0};
    for (size_t i = 0; i < all_check_qubits.size(); i++)
    {
        bool is_x_check = schedule_info.x_check_set.count(all_check_qubits.at(i));
        if (is_x_check != is_memory_x)
            continue;

        uint32_t base_meas_id = static_cast<uint32_t>(n_check_meas-i) | stim::TARGET_RECORD_BIT;
        uint32_t prev_meas_id = static_cast<uint32_t>(2*n_check_meas-i) | stim::TARGET_RECORD_BIT;
        round_circuit.safe_append_u("DETECTOR", {prev_meas_id, base_meas_id}, {detector_id, 0});
        detector_id++;
    }
    // shift detectors:
    round_circuit.safe_append_u("SHIFT_COORDS", {}, {0,1});

    ///////////////////////
    /// EPILOG          ///
    ///////////////////////
    
    if (is_memory_x)
        epilog.safe_append_u("H", schedule_info.data_qubits);
    epilog.safe_append_u("M", schedule_info.data_qubits);

    // compare data qubit measurements with measurement outcomes in last round to handle measurement errors
    std::unordered_map<stim_qubit_type, size_t> dq_meas_order;
    for (size_t i = 0; i < schedule_info.data_qubits.size(); i++)
        dq_meas_order[schedule_info.data_qubits.at(i)] = i;

    const size_t n_data_meas = schedule_info.data_qubits.size();
    detector_id = 0;
    for (size_t i = 0; i < all_check_qubits.size(); i++)
    {
        bool is_x_check = schedule_info.x_check_set.count(all_check_qubits.at(i));
        if (is_x_check != is_memory_x)
            continue;

        auto check_qubit = all_check_qubits.at(i);

        uint32_t check_meas_id = static_cast<uint32_t>(n_data_meas + n_check_meas - i) | stim::TARGET_RECORD_BIT;
        std::vector<uint32_t> meas_list{check_meas_id};
        for (auto q : schedule_info.check_cx_order.at(check_qubit))
        {
            if (q == SC_SCHEDULE_INFO::NO_QUBIT)
                continue;
            meas_list.push_back((n_data_meas - dq_meas_order.at(q)) | stim::TARGET_RECORD_BIT);
        }
        epilog.safe_append_u("DETECTOR", meas_list, {detector_id, 0});
        detector_id++;
    }

    const auto& obs = is_memory_x ? schedule_info.x_obs : schedule_info.z_obs;
    std::vector<uint32_t> obs_meas_id;
    std::transform(obs.begin(), obs.end(), std::back_inserter(obs_meas_id),
                    [&dq_meas_order] (stim_qubit_type q) { return (dq_meas_order[q]+1) | stim::TARGET_RECORD_BIT; });
    epilog.safe_append_ua("OBSERVABLE_INCLUDE", obs_meas_id, 0);

    // create final circuit:
    stim::Circuit fin;
    fin += prolog;
    fin += error_free_round_circuit;
    fin += round_circuit * (rounds);
    fin += epilog;

    return fin;
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

}   // namespace gen
