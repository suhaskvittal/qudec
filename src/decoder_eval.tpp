/*
 *  author: Suhas Vittal
 *  date:   12 October 2025
 */

#include "stim/simulators/frame_simulator.h"
#include "decoder/surface_code.h"

#include <chrono>
#include <type_traits>

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

// Type trait to check if IMPL is PYMATCHING
template<typename T>
constexpr bool is_pymatching_v = std::is_same_v<T, PYMATCHING>;

template <class IMPL> void
decode(IMPL& impl, 
        DECODER_STATS& stats,
        syndrome_type detector_flips,
        syndrome_type observable_flips,
        bool do_not_clock)
{
    // create detector list from `detector_flips`
    std::vector<GRAPH_COMPONENT_ID> detector_list;
    for (size_t i = 0; i < detector_flips.num_bits_padded(); i++)
    {
        if (detector_flips[i])
            detector_list.push_back(static_cast<GRAPH_COMPONENT_ID>(i));
    }

    size_t hw = detector_list.size();

    stats.trials++;
    stats.trials_by_hamming_weight[hw]++;

    // if there are no detector flips, then exit early:
    if (detector_list.empty())
    {
        stats.trivial_trials++;
        return;
    }

    // start clock:
    std::chrono::steady_clock::time_point start_time, end_time;
    if (!do_not_clock)
        start_time = std::chrono::steady_clock::now();

    std::stringstream debug_strm;
    DECODER_RESULT result;

    // Use constexpr type checking for PyMatching (allows for future customization)
    if constexpr (is_pymatching_v<IMPL>) {
        // PyMatching: uses PyMatching's internal MWPM algorithm
        result = impl.decode(detector_list, debug_strm);
    } else {
        // Other decoders: use standard interface
        result = impl.decode(detector_list, debug_strm);
    }

    if (!do_not_clock)
        end_time = std::chrono::steady_clock::now();
    uint64_t time_us = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
    stats.total_time_us += time_us;
    stats.time_us_by_hamming_weight[hw] += time_us;

    // check if result is an error:
    bool any_mismatch{false};
    for (size_t i = 0; i < observable_flips.num_bits_padded(); i++)
    {
        bool truth_is_flipped = observable_flips[i];
        bool prediction_is_flipped = result.flipped_observables.count(static_cast<int64_t>(i)) > 0;
        any_mismatch |= (truth_is_flipped ^ prediction_is_flipped);
    }

    stats.errors += any_mismatch;

#if defined (DEBUG_DECODER)
    if (any_mismatch)
    {
        std::cerr << "TRIAL " << stats.trials << " ==================================== \n";
        std::cerr << "detectors =";
        for (auto d : detector_list)
            std::cerr << " " << d;
        std::cerr << "\ndecoder debug out--------------\n\n";
        std::cerr << debug_strm.str() << "\n";
        std::cerr << "------------------------------\n\n";
        std::cerr << "true flipped observables:";
        for (size_t i = 0; i < observable_flips.num_bits_padded(); i++)
        {
            if (observable_flips[i])
                std::cerr << " " << i;
        }
        std::cerr << "\n";
    }
#endif
}

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

template <class IMPL> DECODER_STATS
benchmark_decoder(const stim::Circuit& circuit,     
                    IMPL& impl, 
                    uint64_t num_trials,
                    uint64_t batch_size,
                    bool do_not_clock,
                    uint64_t seed)
{
    using frame_sim_type = stim::FrameSimulator<stim::MAX_BITWORD_WIDTH>;

    std::mt19937_64 rng(seed);

    size_t num_batches{0};

    DECODER_STATS stats;
    while (num_trials)
    {
        if (num_batches % 20 == 0)
            std::cout << "\n[ trials remaining = " << std::setw(12) << std::right << num_trials << " ]\t";
        (std::cout << ".").flush();

        uint64_t trials_this_batch = std::min(num_trials, batch_size);
        num_trials -= trials_this_batch;

        frame_sim_type sim(circuit.compute_stats(), 
                            stim::FrameSimulatorMode::STORE_DETECTIONS_TO_MEMORY,
                            trials_this_batch,
                            std::move(rng));
        sim.do_circuit(circuit);

        auto detector_table = std::move(sim.det_record.storage);
        auto observable_table = std::move(sim.obs_record);

        // transpose the tables (currently, indices correspond to [detector,shot])
        detector_table = detector_table.transposed();
        observable_table = observable_table.transposed();

        for (uint64_t s = 0; s < trials_this_batch; s++)
            decode(impl, stats, std::move(detector_table[s]), std::move(observable_table[s]), do_not_clock);

        rng = std::move(sim.rng);

        num_batches++;
    }
    std::cout << "\n";

    return stats;
}

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
