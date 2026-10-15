/* Contributor: OpenAI Codex. Seeded decoder integration regression. */
#include "circuit_generator.h"
#include "decoder/seeded.h"

#include <stim/util_top/circuit_to_dem.h>

#include <stdexcept>
#include <iostream>

int main()
{
    const auto circuit = sc_si1000(3, 3, 1e-3, false, true);
    const auto dem = stim::circuit_to_dem(circuit, {.decompose_errors=true});
    dec::PyMatching reference(dem, {.correlated=true});
    dec::PyMatching seed(dem, {.report_matching_data=true});
    dec::CorrelatedPyMatching second(dem);
    dec::SeededMatching pipeline(dem, seed, second);
    dec::ObsType truth(dem.count_observables());
    dec::SyndromeType syndrome(dem.count_detectors());

    // No filtering means the externally seeded second pass must reproduce
    // PyMatching's own correlated second pass shot by shot.
    for (uint64_t pattern = 0; pattern < 64; ++pattern)
    {
        syndrome.clear();
        for (size_t bit = 0; bit < 6 && bit < dem.count_detectors(); ++bit)
            syndrome[bit] = (pattern >> bit) & 1;
        const auto expected = reference.decode(syndrome, truth);
        const auto actual = pipeline.decode(syndrome, truth);
        if (expected.flipped_obs != actual.flipped_obs)
            throw std::runtime_error("Seeded correlated matching disagrees with the native two-pass decoder");
    }

    // Exercise both gap choices and the optional physical-path filter.
    dec::PyMatching gap_seed(dem, {.report_matching_data=true, .enable_cluster_gap=true});
    dec::Libra libra(dem, {.ensemble_size=2});
    dec::SeededMatching cluster_pipeline(dem, gap_seed, libra,
        {.postprocess=true, .libra_gap_threshold_db=1000});
    dec::PyMatching complementary_seed(dem, {.report_matching_data=true});
    dec::SeededMatching complementary_pipeline(dem, complementary_seed, libra,
        {.libra_gap_mode=dec::SeedGapMode::Complementary, .libra_gap_threshold_db=1000});
    for (uint64_t pattern = 1; pattern < 8; ++pattern)
    {
        syndrome.clear();
        for (size_t bit = 0; bit < 3; ++bit)
            syndrome[bit] = (pattern >> bit) & 1;
        (void)cluster_pipeline.decode(syndrome, truth);
        (void)complementary_pipeline.decode(syndrome, truth);
    }
    std::cout << "Seeded matching tests passed\n";
}
