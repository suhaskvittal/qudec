/* Contributor: OpenAI Codex. Pair reporting and configurable PDZ composition. */
#include "circuit_generator.h"
#include "decoder/pdz.h"
#include "decoder/surface_code.h"
#include "decoder/uf.h"
#include <stim/util_top/circuit_to_dem.h>
#include <iostream>
#include <random>
#include <stdexcept>

namespace
{
void check_pairs(const dec::Result& result, dec::SyndromeRef syndrome,
                 size_t detectors, size_t observables)
{
    std::vector<size_t> degree(detectors);
    dec::ObsType parity(observables);
    for (const auto& pair : result.matching_data.assignments)
    {
        for (auto d : {pair.d1, pair.d2})
        {
            if (d < -1 || d >= static_cast<int64_t>(detectors))
                throw std::runtime_error("invalid matched endpoint");
            if (d >= 0) ++degree[d];
        }
        parity ^= pair.frame_flips;
    }
    for (size_t d = 0; d < detectors; ++d)
        if (degree[d] != static_cast<size_t>(syndrome[d]))
            throw std::runtime_error("matching does not cover fired detectors exactly once");
    if (parity != result.flipped_obs)
        throw std::runtime_error("path observable parity differs from UF correction");
}

template <typename Initial>
void check_pipeline(const stim::DetectorErrorModel& dem, Initial& initial,
                    dec::SyndromeRef syndrome, dec::ObsRef truth)
{
    dec::PyMatching residual(dem, {.correlated=true});
    dec::PDZDecoder pipeline(dem, &initial, &residual, {}, 1.0);
    const auto original = dec::SyndromeType(syndrome);
    const auto decoded = pipeline.decode_detailed(syndrome, truth);
    if (original != dec::SyndromeType(syndrome))
        throw std::runtime_error("pipeline modified the input syndrome");
    if (decoded.bypassed)
    {
        if (decoded.remaining.popcnt() || decoded.result.flipped_obs != decoded.initial.flipped_obs)
            throw std::runtime_error("threshold-one bypass is wrong");
    }
    else
    {
        dec::ObsType residual_truth(truth);
        residual_truth ^= decoded.partial.flipped_obs;
        auto expected = residual.decode(decoded.remaining, residual_truth).flipped_obs;
        expected ^= decoded.partial.flipped_obs;
        if (expected != decoded.result.flipped_obs)
            throw std::runtime_error("PDZ correction was not combined with residual correction");
    }
    dec::PDZDecoder bypass(dem, &initial, &residual, {}, 0.0);
    const auto bypassed = bypass.decode_detailed(syndrome, truth);
    if (!bypassed.bypassed || bypassed.result.flipped_obs != bypassed.initial.flipped_obs)
        throw std::runtime_error("threshold-zero bypass changed initial prediction");
}
}

int main()
{
    const stim::DetectorErrorModel toy(
        "error(0.15) D0 D1 L0\nerror(0.2) D1 D2 L1\nerror(0.12) D1 D3\n"
        "error(0.01) D0\nerror(0.02) D2 L0\nerror(0.03) D3 L1\n"
        "error(0.05) D0 D2\n");
    dec::UnionFind ordinary(toy), reporting(toy, {.report_matching_data=true});
    dec::ObsType truth(2);
    for (size_t bits = 0; bits < 16; ++bits)
    {
        dec::SyndromeType syndrome(4);
        for (size_t d = 0; d < 4; ++d) syndrome[d] = (bits >> d) & 1;
        const auto actual = reporting.decode(syndrome, truth);
        if (actual.flipped_obs != ordinary.decode(syndrome, truth).flipped_obs)
            throw std::runtime_error("reporting changed UF prediction");
        check_pairs(actual, syndrome, 4, 2);
    }
    for (int distance : {3, 5})
    {
        const auto dem = stim::circuit_to_dem(sc_si1000(distance, distance, 0.001, false, true),
                                            {.decompose_errors=true});
        dec::UnionFind uf(dem, {.report_matching_data=true}), plain_uf(dem);
        dec::PyMatching matching(dem, {.report_matching_data=true});
        dec::ObsType actual_obs(dem.count_observables());
        std::mt19937_64 rng(42);
        for (size_t shot = 0; shot < 32; ++shot)
        {
            dec::SyndromeType syndrome(dem.count_detectors());
            // Empty, individual DEM mechanisms, and sparse mixed mechanisms.
            dem.iter_flatten_error_instructions([&](const auto& instruction) {
                if (shot == 0 || rng() % 100 >= 2) return;
                for (const auto& target : instruction.target_data)
                    if (target.is_relative_detector_id()) syndrome[target.val()] ^= true;
            });
            const auto out = uf.decode(syndrome, actual_obs);
            check_pairs(out, syndrome, dem.count_detectors(), dem.count_observables());
            if (out.flipped_obs != plain_uf.decode(syndrome, actual_obs).flipped_obs)
                throw std::runtime_error("reporting changed surface-code UF prediction");
            check_pipeline(dem, uf, syndrome, actual_obs);
            check_pipeline(dem, matching, syndrome, actual_obs);
        }
    }
    std::cout << "UF pair coverage/parity and both PDZ pipelines passed\n";
}
