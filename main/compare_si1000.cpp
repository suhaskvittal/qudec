/*
 * author: OpenAI GPT-6-Sol
 * date: 24 September 2026
 * purpose: Compare Tesseract and PyMatching on identical SI1000 DEM samples,
 *          recording every shot where only PyMatching has a logical error.
 */

#include "circuit_generator.h"
#include "decoder/surface_code.h"
#include "decoder/tesseract.h"
#include "sampler.h"

#include <stim.h>
#include <stim/util_top/circuit_to_dem.h>

#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>

int main(int argc, char** argv)
{
    if (argc != 4 && argc != 5 && argc != 6)
    {
        std::cerr << "usage: compare_si1000 OUTPUT_PREFIX SHOTS SEED [paper [ORDER_SEED] | UNIFORM_ERROR_COUNT]\n";
        return 2;
    }
    const std::string prefix = argv[1];
    const uint64_t shots = std::stoull(argv[2]);
    const uint64_t seed = std::stoull(argv[3]);
    // OpenAI GPT-6-Sol: Optional fixed-count samples seek rare discordances;
    // their unweighted frequencies are deliberately separate from Monte Carlo.
    // OpenAI GPT-6-Sol: Keep the prior default and conditioned CLI behavior;
    // paper mode selects Appendix B's short beam and coordinate ensemble.
    const bool paper = argc >= 5 && std::string(argv[4]) == "paper";
    const bool conditioned = argc == 5 && !paper;
    const size_t uniform_error_count = conditioned ? std::stoull(argv[4]) : 0;
    const uint64_t order_seed = argc == 6 ? std::stoull(argv[5]) : seed;
    if (argc == 6 && !paper)
        throw std::invalid_argument("six arguments require paper mode");
    constexpr uint32_t d = 9;
    constexpr uint32_t rounds = 9;
    constexpr double p = 1e-3;
    constexpr size_t batch_size = 1024;

    auto circuit = sc_si1000(d, rounds, p, false, true);
    auto dem = stim::circuit_to_dem(circuit, {.decompose_errors = true});
    dec::PyMatching pymatching(dem, false);
    // OpenAI GPT-6-Sol: Match the paper's short-beam search configuration.
    tesseract_decoder::TesseractConfig config{.dem=dem};
    if (paper)
    {
        config.det_beam = 15;
        config.beam_climbing = true;
        config.detector_orders = tesseract_decoder::make_detector_orders(
            16, tesseract_decoder::DetectorOrder::Method::Coordinate, order_seed);
        config.pqlimit = 200000;
        config.no_revisit_dets = true;
        config.det_penalty = 0;
    }
    dec::Tesseract tesseract(std::move(config));

    std::ofstream circuit_out(prefix + ".stim");
    std::ofstream dem_out(prefix + ".dem");
    std::ofstream events(prefix + ".jsonl");
    if (!circuit_out || !dem_out || !events)
        throw std::runtime_error("failed to open comparison output files");
    circuit_out << circuit;
    dem_out << dem;
    circuit_out.close();
    dem_out.close();

    const size_t detector_count = dem.count_detectors();
    const size_t observable_count = dem.count_observables();
    if (observable_count != 1)
        throw std::runtime_error("expected one logical observable");

    events << "{\"type\":\"metadata\",\"schema\":1,\"distance\":" << d
           << ",\"rounds\":" << rounds << ",\"p\":" << p
           << ",\"memory_basis\":\"Z\",\"include_opposite_basis_detectors\":true"
           << ",\"detector_count\":" << detector_count
           << ",\"observable_count\":" << observable_count
           << ",\"requested_shots\":" << shots << ",\"seed\":" << seed
           << ",\"sampling\":\""
           << (conditioned ? "uniform distinct DEM error instructions at fixed count" :
                             "Stim DemSampler from decomposed circuit DEM") << "\""
           << (conditioned ? ",\"fixed_error_count\":" + std::to_string(uniform_error_count) : "")
           << ",\"detector_ids\":\"zero-based IDs in circuit order\""
           << ",\"tesseract_configuration\":\"" << (paper ? "paper_short_beam" : "default") << "\""
           << (paper ? ",\"det_beam\":15,\"beam_climbing\":true,\"detector_orders\":16,\"order_method\":\"coordinate\",\"order_seed\":" + std::to_string(order_seed) + ",\"pqlimit\":200000,\"no_revisit_dets\":true,\"det_penalty\":0" : "")
           << "}\n";
    events.flush();

    stim::DemSampler<64> sampler(dem, std::mt19937_64(seed), batch_size);
    std::mt19937_64 conditioned_rng(seed);
    uint64_t completed = 0;
    uint64_t pymatching_errors = 0;
    uint64_t tesseract_errors = 0;
    uint64_t only_pymatching_errors = 0;
    uint64_t only_tesseract_errors = 0;
    uint64_t both_errors = 0;
    uint64_t tesseract_low_confidence = 0;
    uint64_t captured_low_confidence = 0;
    uint64_t paper_scored_tesseract_errors = 0;
    uint64_t paper_scored_only_pymatching_errors = 0;
    uint64_t uncaptured_raw_only_pymatching = 0;
    uint64_t captured_detector_total = 0;
    uint64_t captured_opposite_detector_total = 0;
    const uint64_t first_body = (d * d - 1) / 2;
    const uint64_t per_body = d * d - 1;
    const uint64_t opposite_per_body = per_body / 2;

    while (completed < shots)
    {
        auto [detectors, observables] = [&]() -> Problem {
            if (conditioned)
                return generate_syndromes_with_k_errors(dem, uniform_error_count, batch_size, conditioned_rng);
            sampler.resample(false);
            return {sampler.det_buffer.transposed(), sampler.obs_buffer.transposed()};
        }();
        for (size_t i = 0; i < batch_size && completed < shots; ++i, ++completed)
        {
            auto pm_result = pymatching.decode(detectors[i], observables[i]);
            auto te_result = tesseract.decode(detectors[i], observables[i]);
            const bool truth = observables[i][0];
            const bool pm_prediction = pm_result.flipped_obs[0];
            const bool te_prediction = te_result.flipped_obs[0];
            const bool pm_error = pm_prediction != truth;
            const bool te_error = te_prediction != truth;
            const bool low_confidence = tesseract.last_low_confidence();
            pymatching_errors += pm_error;
            tesseract_errors += te_error;
            tesseract_low_confidence += low_confidence;
            paper_scored_tesseract_errors += te_error || low_confidence;
            paper_scored_only_pymatching_errors += pm_error && !te_error && !low_confidence;
            only_pymatching_errors += pm_error && !te_error && (!paper || !low_confidence);
            uncaptured_raw_only_pymatching += paper && pm_error && !te_error && low_confidence;
            only_tesseract_errors += !pm_error && te_error;
            both_errors += pm_error && te_error;

            if (pm_error && !te_error && (!paper || !low_confidence))
            {
                events << "{\"type\":\"shot\",\"shot_index\":" << completed
                       << ",\"truth\":" << int(truth)
                       << ",\"pymatching_prediction\":" << int(pm_prediction)
                       << ",\"tesseract_prediction\":" << int(te_prediction)
                       << ",\"tesseract_low_confidence\":" << (low_confidence ? "true" : "false")
                       << ",\"detector_ids\":[";
                bool first = true;
                uint64_t opposite_count = 0;
                for (size_t det = 0; det < detector_count; ++det)
                {
                    if (!detectors[i][det])
                        continue;
                    if (!first) events << ',';
                    events << det;
                    first = false;
                    captured_detector_total++;
                    if (det >= first_body && det < first_body + (rounds - 1) * per_body &&
                        (det - first_body) % per_body >= opposite_per_body)
                        opposite_count++;
                }
                captured_opposite_detector_total += opposite_count;
                captured_low_confidence += low_confidence;
                events << "],\"opposite_basis_detection_count\":" << opposite_count << "}\n";
                // OpenAI GPT-6-Sol: Keep rare captured examples durable during long runs.
                events.flush();
            }
        }
        if (completed % 102400 == 0)
        {
            std::cerr << "completed " << completed << " / " << shots << " shots\n";
            events.flush();
        }
    }

    events << "{\"type\":\"summary\",\"completed_shots\":" << completed
           << ",\"pymatching_errors\":" << pymatching_errors
           << ",\"tesseract_errors\":" << tesseract_errors
           << ",\"only_pymatching_errors\":" << only_pymatching_errors
           << ",\"only_tesseract_errors\":" << only_tesseract_errors
           << ",\"both_errors\":" << both_errors
           << ",\"tesseract_low_confidence\":" << tesseract_low_confidence
           << ",\"captured_low_confidence\":" << captured_low_confidence
           << ",\"captured_detector_total\":" << captured_detector_total
           << ",\"captured_opposite_detector_total\":" << captured_opposite_detector_total
           << ",\"paper_scored_tesseract_errors\":" << paper_scored_tesseract_errors
           << ",\"paper_scored_only_pymatching_errors\":" << paper_scored_only_pymatching_errors
           << ",\"uncaptured_raw_only_pymatching\":" << uncaptured_raw_only_pymatching
           << "}\n";
    std::cout << "shots=" << completed << " pymatching_errors=" << pymatching_errors
              << " tesseract_errors=" << tesseract_errors
              << " only_pymatching=" << only_pymatching_errors
              << " only_tesseract=" << only_tesseract_errors
              << " both=" << both_errors << "\n";
}
