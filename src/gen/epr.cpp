/*
 *  author: Suhas Vittal
 *  date:   12 October 2025
 * */

#include "gen/epr.h"
#include "gen/scheduling.h"
#include "gen/utils.h"

#include <algorithm>
#include <cmath>
#include <numeric>

namespace gen
{

//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

SC_EPR_SCHEDULE_INFO::SC_EPR_SCHEDULE_INFO(size_t d, bool dual)
    :sc(d, 2*d+1, dual)
{
    // these are the middle rows:
    const size_t mr1 = d,
                 mr2 = d+1;

    // identify all checks in row `mr1` as these will use EPR pairs (also create new CX schedules)
    stim_qubit_type epr_check_qubit{2*(2*d+1)*d - 1};
    for (size_t c = 0; c <= d; c++)
    {
        // assign `c` an epr qubit
        auto q = sc.check_matrix[mr1][c];
        if (q == NO_QUBIT)
            continue;
        epr_checks[q] = epr_check_qubit;

        std::cout << "creating epr pair for check " << q << " --> " << epr_check_qubit << "\n";

        // update CX order
        auto cx_it = sc.check_cx_order.find(q);
        auto orig_cx_order = std::move(cx_it->second);
        sc.check_cx_order.erase(cx_it);

        // split `orig_cx_order` (as there are two parity qubits for the same check now)
        std::vector<stim_qubit_type> cx_order_1, cx_order_2;

        if (sc.x_check_set.count(q))
        {
            // check order is nw, ne, sw, se
            cx_order_1 = {orig_cx_order[0], orig_cx_order[1], NO_QUBIT, NO_QUBIT};
            cx_order_2 = {NO_QUBIT, NO_QUBIT, orig_cx_order[2], orig_cx_order[3]};
        }
        else
        {
            // check order is nw, sw, ne, se
            cx_order_1 = {orig_cx_order[0], NO_QUBIT, orig_cx_order[2], NO_QUBIT};
            cx_order_2 = {NO_QUBIT, orig_cx_order[1], NO_QUBIT, orig_cx_order[3]};
        }

        epr_cx_order[q] = cx_order_1;
        epr_cx_order[epr_check_qubit] = cx_order_2;

        if (sc.x_check_set.count(q))
            sc.x_check_set.insert(epr_check_qubit);

        epr_check_qubit++;
    }

    // re-create qubit coordinates and define `hw1_qubit_set` and `hw2_qubit_set`
    const size_t hw2_y_offset = 2;
    for (size_t r = 0; r < 2*d+1; r++)
    {
        for (size_t c = 0; c < d; c++)
        {
            stim_qubit_type q = r*d + c;
            if (r < d)
            {
                qubit_coords[q] = std::make_pair(c+0.5, r+0.5);
                hw1_qubit_set.insert(q);
            }
            else
            {
                qubit_coords[q] = std::make_pair(c+0.5, r+0.5 + hw2_y_offset);
                hw2_qubit_set.insert(q);
            }
        }
    }

    for (size_t r = 0; r <= 2*d+1; r++)
    {
        for (size_t c = 0; c <= d; c++)
        {
            auto q = sc.check_matrix[r][c];
            if (q == NO_QUBIT)
                continue;
            auto epr_it = epr_checks.find(q);
            if (epr_it != epr_checks.end())
            {
                auto epr_q = epr_it->second;

                // make two coordinates:
                qubit_coords[q] = std::make_pair(c, r);
                qubit_coords[epr_q] = std::make_pair(c, r + hw2_y_offset);

                hw1_qubit_set.insert(q);
                hw2_qubit_set.insert(epr_q);
            }
            else if (r < mr1)
            {
                qubit_coords[q] = std::make_pair(c, r);
                hw1_qubit_set.insert(q);
            }
            else
            {
                qubit_coords[q] = std::make_pair(c, r + hw2_y_offset);
                hw2_qubit_set.insert(q);
            }
        }
    }
}

//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

sc_epr_gen_output_type
sc_epr_generation(const EPR_GEN_CONFIG& config, size_t rounds, size_t distance, bool do_memory_experiment)
{
    // merged surface code has asymmetric distance for ZZ measurement
    SC_EPR_SCHEDULE_INFO epr(distance);

    // identify data and ancilla Z checks in the ancilla region:
    std::vector<stim_qubit_type> anc_data_qubits(distance);
    std::vector<stim_qubit_type> anc_z_checks;

    std::iota(anc_data_qubits.begin(), anc_data_qubits.end(), distance*distance);

    const size_t mr1{distance},
                 mr2{distance+1};
    for (size_t c = 0; c < distance; c += 2)
        anc_z_checks.push_back(epr.sc.check_matrix[mr1][c]);
    for (size_t c = 1; c <= distance; c += 2)
        anc_z_checks.push_back(epr.sc.check_matrix[mr2][c]);

    // create hardware error and timing setup:
    double latency_diff = static_cast<double>(config.hw2_round_ns) / static_cast<double>(config.hw1_round_ns);

    double p = config.phys_error;
    double t1_ns_hw1{500.0/p},
           t2_ns_hw1{250.0/p};
    double t1_ns_hw2{t1_ns_hw1 * latency_diff}, 
           t2_ns_hw2{t2_ns_hw1 * latency_diff};
    double e_readout{3*p},
           e_g1q{p*0.1},
           e_g2q{p},
           e_idle{p*0.1};

    // identify number of super rounds (rounds using EPR pairs) based on attenuation rate:
    double prob_any_attenuation = (epr.epr_checks.size() * config.attenuation_rate);
    double expected_rounds_with_photon_loss = rounds * prob_any_attenuation;
    size_t num_super_rounds = rounds + static_cast<size_t>(std::round(expected_rounds_with_photon_loss));
#if defined(EPR_ONLY_ONE_HW1_ROUND)
    size_t num_hw1_rounds_per_super_round = 1;
#else
    size_t num_hw1_rounds_per_super_round = static_cast<size_t>(std::ceil(latency_diff)) - 1;
#endif
     
    // now that have completed initializing all data structures: create the circuit
    // note stability experiment will not use error free rounds
    stim::Circuit prolog,
                    first_round,
                    hw1_only_first_round,
                    hw1_only_main_round,
                    last_round, 
                    super_round, 
                    epilog;

    // decoder circuits:
    stim::Circuit fp_first_round, fp_main_round, fp_last_round;
    stim::Circuit sp_first_round, sp_main_round, sp_last_round;

    // prolog: add coordinates:
    for (const auto& [q, coord] : epr.qubit_coords)
        prolog.safe_append_u("QUBIT_COORDS", {q}, {coord.first, coord.second});

    // init data qubits: unlike `sc_memory` or `sc_stability` -- initialize these data qubits with noise:
    prolog.safe_append_u("R", epr.sc.data_qubits);
    prolog.safe_append_ua("X_ERROR", epr.sc.data_qubits, e_g1q);

    if (do_memory_experiment)
    {
        prolog.safe_append_u("H", epr.sc.data_qubits);
        prolog.safe_append_ua("DEPOLARIZE1", epr.sc.data_qubits, e_g1q);
    }
    else
    {
        // only `anc_data_qubits` are initialized in X basis for stability experiment
        prolog.safe_append_u("H", anc_data_qubits);
        prolog.safe_append_ua("DEPOLARIZE1", anc_data_qubits, e_g1q);
    }
    prolog.safe_append_u("TICK", {}); 

    // first round (again, erroneous since this is an initialization experiment):
    auto super_check_meas_map = sc_epr_create_super_round(first_round, 
                                                            epr, 
                                                            config.hw1_round_ns,
                                                            config.hw2_round_ns,
                                                            t1_ns_hw1, 
                                                            t2_ns_hw1, 
                                                            t1_ns_hw2, 
                                                            t2_ns_hw2, 
                                                            e_readout, 
                                                            e_g1q, 
                                                            e_g2q, 
                                                            e_idle, 
                                                            config.photonic_link_error);
    // note `super_round` has the same operations as `first_round`
    super_round = first_round;
    last_round = first_round.without_noise();

    // hw1 only round does not measure stabilizers that require EPR pairs:
    auto hw1_only_check_meas_map = sc_epr_create_hw1_only_circuit(hw1_only_first_round, 
                                                                    epr, 
                                                                    config.hw1_round_ns,
                                                                    t1_ns_hw1, 
                                                                    t2_ns_hw1, 
                                                                    e_readout, 
                                                                    e_g1q, 
                                                                    e_g2q, 
                                                                    e_idle);
    hw1_only_main_round = hw1_only_first_round;

    // decoder first pass circuit is like hw1 only
    fp_first_round = hw1_only_first_round;
    fp_main_round = hw1_only_main_round;

    // decoder second pass circuit only has super rounds, but error of hw1 needs to be reduced
    sc_epr_create_super_round(sp_first_round,
                                epr, 
                                config.hw1_round_ns,
                                config.hw2_round_ns,
                                t1_ns_hw1, 
                                t2_ns_hw1, 
                                t1_ns_hw2, 
                                t2_ns_hw2, 
                                e_readout, 
                                e_g1q, 
                                e_g2q, 
                                e_idle, 
                                config.attenuation_rate,
                                1.0 / latency_diff);
    sp_main_round = sp_first_round;

    // create detection events:
    const auto& det_checks = do_memory_experiment ? epr.sc.x_check_qubits : epr.sc.z_check_qubits;
    std::vector<stim_qubit_type> first_round_det_checks(det_checks);
    if (!do_memory_experiment)
    {
        auto it = std::remove_if(first_round_det_checks.begin(), first_round_det_checks.end(),
                                [&anc_z_checks] (auto q) 
                                { 
                                    return std::find(anc_z_checks.begin(), anc_z_checks.end(), q) != anc_z_checks.end();
                                });
        first_round_det_checks.erase(it, first_round_det_checks.end());
    }

    // epr circuit detection events:
    
    sc_epr_create_detection_events_super_round(first_round, 
                                                first_round_det_checks, 
                                                super_check_meas_map, 
                                                super_check_meas_map,
                                                num_hw1_rounds_per_super_round, 
                                                true,
                                                epr);

    sc_epr_create_detection_events_adjacent_hw1_rounds(hw1_only_first_round, 
                                                    det_checks, 
                                                    hw1_only_check_meas_map, 
                                                    super_check_meas_map,
                                                    epr);
    sc_epr_create_detection_events_adjacent_hw1_rounds(hw1_only_main_round, 
                                                    det_checks, 
                                                    hw1_only_check_meas_map, 
                                                    hw1_only_check_meas_map,
                                                    epr);
    sc_epr_create_detection_events_super_round(super_round, 
                                                det_checks, 
                                                super_check_meas_map, 
                                                hw1_only_check_meas_map,
                                                num_hw1_rounds_per_super_round, 
                                                false,
                                                epr);
    sc_epr_create_detection_events_generic(last_round, 
                                                det_checks, 
                                                super_check_meas_map,
                                                false,
                                                epr);

    // decoder pass circuit detection events:
    //
    // first pass is only hw1:
    util::check_meas_map empty_meas_map{};
    sc_epr_create_detection_events_generic(fp_first_round, 
                                                det_checks, 
                                                hw1_only_check_meas_map,
                                                true,
                                                epr);
    sc_epr_create_detection_events_generic(fp_main_round, 
                                                det_checks, 
                                                hw1_only_check_meas_map,
                                                false,
                                                epr);
    fp_last_round = fp_main_round.without_noise();

    // second pass is only super rounds:
    sc_epr_create_detection_events_generic(sp_first_round, 
                                                det_checks, 
                                                super_check_meas_map,
                                                true,
                                                epr);
    sc_epr_create_detection_events_generic(sp_main_round, 
                                                det_checks, 
                                                super_check_meas_map,
                                                false,
                                                epr);
    sp_last_round = sp_main_round.without_noise();

    // create epilog:
    epilog.safe_append_u("TICK", {});

    // measure out data qubits:
    if (do_memory_experiment)
        epilog.safe_append_u("H", epr.sc.data_qubits);
    else
        epilog.safe_append_u("H", anc_data_qubits);
    epilog.safe_append_u("M", epr.sc.data_qubits);

    const size_t n_data_meas = epr.sc.data_qubits.size();
    const size_t n_check_meas = super_check_meas_map.size();
    std::unordered_map<stim_qubit_type, uint32_t> dq_meas_order;
    for (size_t i = 0; i < n_data_meas; i++)
        dq_meas_order[epr.sc.data_qubits[i]] = i;
    std::vector<uint32_t> obs_meas_id;
    if (do_memory_experiment)
    {
        // observable is already defined in `epr.sc.x_obs`
        std::transform(epr.sc.x_obs.begin(), epr.sc.x_obs.end(), std::back_inserter(obs_meas_id),
                        [n_data_meas, &dq_meas_order] (stim_qubit_type q) 
                        {
                            return (n_data_meas - dq_meas_order.at(q)) | stim::TARGET_RECORD_BIT;
                        });
    }
    else
    {
        // observable is the Z measurements in the last round:
        for (auto q : anc_z_checks)
        {
            uint32_t base_meas_id = (n_data_meas+n_check_meas - super_check_meas_map.at(q)) | stim::TARGET_RECORD_BIT;
            obs_meas_id.push_back(base_meas_id);

            if (epr.epr_checks.count(q))
            {
                auto e = epr.epr_checks.at(q);
                uint32_t e_meas_id = (n_data_meas+n_check_meas - super_check_meas_map.at(e)) | stim::TARGET_RECORD_BIT;
                obs_meas_id.push_back(e_meas_id);
            }
        }
    }
    epilog.safe_append_ua("OBSERVABLE_INCLUDE", obs_meas_id, 0);

    stim::Circuit composite_round;
    composite_round += hw1_only_first_round;
#if !defined(EPR_ONLY_ONE_HW1_ROUND)
    composite_round += hw1_only_main_round * (num_hw1_rounds_per_super_round-1);
#endif
    // shift super-round idx and super-round idx back
    composite_round.safe_append_u("SHIFT_COORDS", {}, {0,1,0,1,-static_cast<int>(num_hw1_rounds_per_super_round)});
    composite_round += super_round;

    stim::Circuit fin;
    fin += prolog;
    fin += first_round;
    fin += composite_round * (num_super_rounds-1);
    if (do_memory_experiment)
        fin += last_round;
    fin += epilog;

    stim::Circuit first_pass, second_pass;

    first_pass += prolog;
    first_pass += fp_first_round;
    first_pass += fp_main_round * (2*distance);
    first_pass += fp_last_round;
    first_pass += epilog;

    second_pass += prolog;
    second_pass += sp_first_round;
    second_pass += sp_main_round * (num_super_rounds-1);
    second_pass += sp_last_round;
    second_pass += epilog;

    return sc_epr_gen_output_type{fin, first_pass, second_pass};
}

//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

util::check_meas_map
sc_epr_create_super_round(stim::Circuit& circuit,
                          const SC_EPR_SCHEDULE_INFO& epr,
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
                          double hw1_error_scale_factor)
{
    std::vector<stim_qubit_type> all_qubits(epr.hw1_qubit_set.begin(), epr.hw1_qubit_set.end());
    all_qubits.insert(all_qubits.end(), epr.hw2_qubit_set.begin(), epr.hw2_qubit_set.end());

    std::vector<stim_qubit_type> all_checks(epr.sc.z_check_qubits);
    all_checks.insert(all_checks.end(), epr.sc.x_check_qubits.begin(), epr.sc.x_check_qubits.end());

    auto [hw1_ex, hw1_ey, hw1_ez] = pauli_twirling_approx(t1_ns_hw1, t2_ns_hw1, hw1_round_ns);
    auto [hw2_ex, hw2_ey, hw2_ez] = pauli_twirling_approx(t1_ns_hw2, t2_ns_hw2, hw2_round_ns);

    // start of round error on data qubits:
    std::vector<stim_qubit_type> hw1_data_qubits, hw2_data_qubits;
    std::copy_if(epr.sc.data_qubits.begin(), epr.sc.data_qubits.end(), std::back_inserter(hw1_data_qubits),
                    [&epr] (auto q) { return epr.hw1_qubit_set.count(q); });
    std::copy_if(epr.sc.data_qubits.begin(), epr.sc.data_qubits.end(), std::back_inserter(hw2_data_qubits),
                    [&epr] (auto q) { return epr.hw2_qubit_set.count(q); });

    circuit.safe_append_u("PAULI_CHANNEL_1", hw1_data_qubits, {hw1_ex, hw1_ey, hw1_ez});
    circuit.safe_append_u("PAULI_CHANNEL_1", hw2_data_qubits, {hw2_ex, hw2_ey, hw2_ez});

    // initialization of non-EPR parity check qubits:
    std::vector<stim_qubit_type> all_checks_without_epr, x_checks_without_epr, z_checks_without_epr;
    std::copy_if(all_checks.begin(), all_checks.end(), std::back_inserter(all_checks_without_epr),
                    [&epr] (auto q) { return !epr.epr_checks.count(q); });
    std::copy_if(epr.sc.x_check_qubits.begin(), epr.sc.x_check_qubits.end(), std::back_inserter(x_checks_without_epr),
                    [&epr] (auto q) { return !epr.epr_checks.count(q); });
    std::copy_if(epr.sc.z_check_qubits.begin(), epr.sc.z_check_qubits.end(), std::back_inserter(z_checks_without_epr),
                    [&epr] (auto q) { return !epr.epr_checks.count(q); });

    circuit.safe_append_u("R", all_checks_without_epr);
    circuit.safe_append_ua("X_ERROR", all_checks_without_epr, e_g1q);
    circuit.safe_append_u("TICK", {});

    circuit.safe_append_u("H", x_checks_without_epr);
    circuit.safe_append_ua("DEPOLARIZE1", x_checks_without_epr, e_g1q);
    circuit.safe_append_u("TICK", {});
    
    // initialization of EPR parity checks:
    std::vector<stim_qubit_type> epr_h_pre_targets, epr_cx_targets, epr_h_post_targets;
    for (const auto& [q, epr_q] : epr.epr_checks)
    {
        epr_h_pre_targets.push_back(q);
        epr_cx_targets.push_back(q);
        epr_cx_targets.push_back(epr_q);

        if (!epr.sc.x_check_set.count(q))
        {
            epr_h_post_targets.push_back(q);
            epr_h_post_targets.push_back(epr_q);
        }
    }

    
    circuit.safe_append_u("R", epr_cx_targets);
    circuit.safe_append_u("H", epr_h_pre_targets);
    circuit.safe_append_u("CX", epr_cx_targets);
    circuit.safe_append_ua("DEPOLARIZE2", epr_cx_targets, e_photonic_link);
    circuit.safe_append_u("H", epr_h_post_targets);
    circuit.safe_append_ua("DEPOLARIZE1", epr_h_post_targets, e_g1q);
    circuit.safe_append_u("TICK", {});

    // do CNOTs:
    for (size_t t = 0; t < 4; t++)
    {
        std::vector<stim_qubit_type> cx_targets;
        for (const auto& [c, cx_order] : epr.sc.check_cx_order)
        {
            stim_qubit_type q1{c},
                            q2{cx_order.at(t)};
            if (q2 == NO_QUBIT)
                continue;

            if (!epr.sc.x_check_set.count(c))
                std::swap(q1,q2);
            cx_targets.push_back(q1);
            cx_targets.push_back(q2);
        }

        for (const auto& [c, cx_order] : epr.epr_cx_order)
        {
            stim_qubit_type q1{c},
                            q2{cx_order.at(t)};
            if (q2 == NO_QUBIT)
                continue;

            if (!epr.sc.x_check_set.count(c))
                std::swap(q1,q2);
            cx_targets.push_back(q1);
            cx_targets.push_back(q2);
        }

        circuit.safe_append_u("CX", cx_targets);
        circuit.safe_append_ua("DEPOLARIZE2", cx_targets, e_g2q);

        std::unordered_set<stim_qubit_type> cx_target_set(cx_targets.begin(), cx_targets.end());

        std::vector<stim_qubit_type> idle_qubits;
        std::copy_if(all_qubits.begin(), all_qubits.end(), std::back_inserter(idle_qubits),
                        [&cx_target_set] (auto q) { return !cx_target_set.count(q); });
        circuit.safe_append_ua("DEPOLARIZE1", idle_qubits, e_idle);

        circuit.safe_append_u("TICK", {});
    }

    // measure qubits:
    std::vector<stim_qubit_type> x_meas_qubits(x_checks_without_epr),
                                 z_meas_qubits(z_checks_without_epr);
    for (const auto& [q, epr_q] : epr.epr_checks)
    {
        auto& meas_qubits = epr.sc.x_check_set.count(q) ? x_meas_qubits : z_meas_qubits;
        meas_qubits.push_back(q);
        meas_qubits.push_back(epr_q);
    }

    std::vector<stim_qubit_type> all_meas_qubits(x_meas_qubits);
    all_meas_qubits.insert(all_meas_qubits.end(), z_meas_qubits.begin(), z_meas_qubits.end());

    circuit.safe_append_u("H", x_meas_qubits);
    circuit.safe_append_ua("DEPOLARIZE1", x_meas_qubits, e_g1q);
    circuit.safe_append_ua("M", all_meas_qubits, e_readout);

    util::check_meas_map check_meas_map;
    for (size_t i = 0; i < all_meas_qubits.size(); i++)
        check_meas_map[all_meas_qubits[i]] = i;

    return check_meas_map;
}

//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

util::check_meas_map
sc_epr_create_hw1_only_circuit(stim::Circuit& circuit,
                                const SC_EPR_SCHEDULE_INFO& epr,
                                uint64_t hw1_round_ns,
                                double t1_ns_hw1,
                                double t2_ns_hw1,
                                double e_readout,
                                double e_g1q,
                                double e_g2q,
                                double e_idle)
{
    std::unordered_set<stim_qubit_type> hw1_qubit_set_no_epr(epr.hw1_qubit_set.begin(), epr.hw1_qubit_set.end());
    for (const auto& [q, epr_q] : epr.epr_checks)
        hw1_qubit_set_no_epr.erase(q);

    std::vector<stim_qubit_type> all_qubits(hw1_qubit_set_no_epr.begin(), hw1_qubit_set_no_epr.end());

    std::vector<stim_qubit_type> all_data_qubits;
    std::vector<stim_qubit_type> x_checks, z_checks;
    std::copy_if(epr.sc.data_qubits.begin(), epr.sc.data_qubits.end(), std::back_inserter(all_data_qubits),
                    [&epr, &hw1_qubit_set_no_epr] (auto q) { return hw1_qubit_set_no_epr.count(q); });
    std::copy_if(epr.sc.z_check_qubits.begin(), epr.sc.z_check_qubits.end(), std::back_inserter(z_checks),
                    [&epr, &hw1_qubit_set_no_epr] (auto q) { return hw1_qubit_set_no_epr.count(q); });
    std::copy_if(epr.sc.x_check_qubits.begin(), epr.sc.x_check_qubits.end(), std::back_inserter(x_checks),
                    [&epr, &hw1_qubit_set_no_epr] (auto q) { return hw1_qubit_set_no_epr.count(q); });

    std::vector<stim_qubit_type> all_checks(x_checks);
    all_checks.insert(all_checks.end(), z_checks.begin(), z_checks.end());

    auto [ex, ey, ez] = pauli_twirling_approx(t1_ns_hw1, t2_ns_hw1, hw1_round_ns);

    circuit.safe_append_u("PAULI_CHANNEL_1", all_data_qubits, {ex,ey,ez});
    circuit.safe_append_u("R", all_checks);
    circuit.safe_append_ua("X_ERROR", all_checks, e_g1q);
    circuit.safe_append_u("H", x_checks);
    circuit.safe_append_ua("DEPOLARIZE1", x_checks, e_g1q);

    // do CNOTs:
    for (size_t t = 0; t < 4; t++)
    {
        std::vector<stim_qubit_type> cx_targets;

        for (auto q : all_checks)
        {
            if (epr.epr_checks.count(q))
                continue;
            stim_qubit_type q1{q},
                            q2{epr.sc.check_cx_order.at(q).at(t)};
            if (q2 == NO_QUBIT)
                continue;

            if (!epr.sc.x_check_set.count(q))
                std::swap(q1,q2);
            cx_targets.push_back(q1);
            cx_targets.push_back(q2);
        }

        circuit.safe_append_u("CX", cx_targets);
        circuit.safe_append_ua("DEPOLARIZE2", cx_targets, e_g2q);

        std::unordered_set<stim_qubit_type> cx_target_set(cx_targets.begin(), cx_targets.end());
        std::vector<stim_qubit_type> idle_qubits;
        std::copy_if(all_qubits.begin(), all_qubits.end(), std::back_inserter(idle_qubits),
                        [&cx_target_set] (auto q) { return !cx_target_set.count(q); });
        circuit.safe_append_ua("DEPOLARIZE1", idle_qubits, e_idle);
    }

    // measure checks:
    circuit.safe_append_u("H", x_checks);
    circuit.safe_append_ua("DEPOLARIZE1", x_checks, e_g1q);
    circuit.safe_append_ua("M", all_checks, e_readout);

    util::check_meas_map result;
    for (size_t i = 0; i < all_checks.size(); i++)
        result[all_checks[i]] = i;
    return result;
}

//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

void
sc_epr_create_detection_events_super_round(stim::Circuit& circuit,
                                            const util::stim_qubit_array& checks,
                                            const util::check_meas_map& cm_super_round,
                                            const util::check_meas_map& cm_hw1_only_round,
                                            size_t hw1_rounds_per_super_round,
                                            bool is_first_round,
                                            const SC_EPR_SCHEDULE_INFO& epr)
{
    const uint32_t n_check_meas_s = cm_super_round.size(),
                    n_check_meas_hw1 = cm_hw1_only_round.size();

    const uint32_t n_meas_between_super_rounds = n_check_meas_hw1*hw1_rounds_per_super_round + n_check_meas_s;

    for (size_t i = 0; i < checks.size(); i++)
    {
        auto q = checks.at(i);

        std::vector<uint32_t> targets;

        auto epr_it = epr.epr_checks.find(q);
        if (epr_it != epr.epr_checks.end())
        {
            auto e = epr_it->second;

            // detection event is the XOR of meausrements on `q` and `e`
            for (auto x : {q,e})
            {
                size_t meas_idx = cm_super_round.at(x);
                uint32_t base_meas_id = (n_check_meas_s - meas_idx) | stim::TARGET_RECORD_BIT;
                targets.push_back(base_meas_id);

                if (!is_first_round)
                {
                    uint32_t prev_meas_id = (n_meas_between_super_rounds+n_check_meas_s - meas_idx) 
                                                | stim::TARGET_RECORD_BIT;
                    targets.push_back(prev_meas_id);
                }
            }
        }
        else
        {
            uint32_t base_meas_id = (n_check_meas_s - cm_super_round.at(q)) | stim::TARGET_RECORD_BIT;
            targets.push_back(base_meas_id);

            if (!is_first_round)
            {
                uint32_t prev_meas_id;
                if (epr.hw1_qubit_set.count(q))
                {
                    prev_meas_id = (n_check_meas_s+n_check_meas_hw1 - cm_hw1_only_round.at(q)) 
                                            | stim::TARGET_RECORD_BIT;
                }
                else
                {
                    prev_meas_id = (n_meas_between_super_rounds+n_check_meas_s - cm_super_round.at(q)) 
                                            | stim::TARGET_RECORD_BIT;
                }
                targets.push_back(prev_meas_id);
            }
        }

        circuit.safe_append_u("DETECTOR", targets, {i,0,q});
    }
    circuit.safe_append_u("SHIFT_COORDS", {}, {0,1,0,0,1});  // shift sub-round idx
    circuit.safe_append_u("TICK", {});
}

//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

void
sc_epr_create_detection_events_adjacent_hw1_rounds(stim::Circuit& circuit, 
                                                const util::stim_qubit_array& checks,
                                                const util::check_meas_map& cm_this_round,
                                                const util::check_meas_map& cm_prev_round,
                                                const SC_EPR_SCHEDULE_INFO& epr)
{
    const uint32_t n_check_meas_this_round = cm_this_round.size();
    const uint32_t n_check_meas_prev_round = cm_prev_round.size();

    for (size_t i = 0; i < checks.size(); i++)
    {
        auto q = checks.at(i);
        if (epr.epr_checks.count(q) || !epr.hw1_qubit_set.count(q))
            continue;

        uint32_t meas_idx_this_round = cm_this_round.at(q);
        uint32_t base_meas_id = (n_check_meas_this_round - meas_idx_this_round) 
                                    | stim::TARGET_RECORD_BIT;

        std::vector<uint32_t> targets{base_meas_id};

        if (cm_prev_round.count(q))
        {
            uint32_t meas_idx_prev_round = cm_prev_round.at(q);
            uint32_t prev_meas_id = (n_check_meas_this_round+n_check_meas_prev_round - meas_idx_prev_round) 
                                | stim::TARGET_RECORD_BIT;
            targets.push_back(prev_meas_id);
        }

        circuit.safe_append_u("DETECTOR", targets, {i,0,q});
    }
    circuit.safe_append_u("SHIFT_COORDS", {}, {0,1,0,0,1});  // shift sub-round idx
    circuit.safe_append_u("TICK", {});
}

//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

void
sc_epr_create_detection_events_generic(stim::Circuit& circuit,
                                            const util::stim_qubit_array& checks,
                                            const util::check_meas_map& cm,
                                            bool is_first_round,
                                            const SC_EPR_SCHEDULE_INFO& epr)
{
    const uint32_t n_check_meas = cm.size();
    for (size_t i = 0; i < checks.size(); i++)
    {
        auto q = checks.at(i);

        if (!cm.count(q))
            continue;

        std::vector<uint32_t> targets;

        auto epr_it = epr.epr_checks.find(q);
        if (epr_it != epr.epr_checks.end())
        {
            auto e = epr_it->second;

            if (!cm.count(e))
                continue;

            // detection event is the XOR of meausrements on `q` and `e`
            for (auto x : {q,e})
            {
                size_t meas_idx = cm.at(x);
                uint32_t base_meas_id = (n_check_meas - meas_idx) | stim::TARGET_RECORD_BIT,
                         prev_meas_id = (2*n_check_meas - meas_idx) | stim::TARGET_RECORD_BIT;
                targets.push_back(base_meas_id);
                if (!is_first_round)
                    targets.push_back(prev_meas_id);
            }
        }
        else
        {
            uint32_t base_meas_id = (n_check_meas - cm.at(q)) | stim::TARGET_RECORD_BIT,
                     prev_meas_id = (2*n_check_meas - cm.at(q)) | stim::TARGET_RECORD_BIT;
            targets.push_back(base_meas_id);
            if (!is_first_round)
                targets.push_back(prev_meas_id);
        }

        circuit.safe_append_u("DETECTOR", targets, {i,0,q});
    }
    circuit.safe_append_u("SHIFT_COORDS", {}, {0,1});
    circuit.safe_append_u("TICK", {});
}

//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

}   // namespace gen   
