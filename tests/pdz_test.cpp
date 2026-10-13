/*
 *  contributor: OpenAI Codex
 * */

#include "preprocessor/pdz.h"
#include "decoder/surface_code.h"
#include "circuit_generator.h"
#include <stim/util_top/circuit_to_dem.h>

#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <string>

namespace
{

void
require(bool condition, const char* message)
{
    if (!condition)
        throw std::runtime_error(message);
}

} // anon

int
main()
{
    std::string text = "error(0.01) D0 L0\n";
    for (size_t d = 0; d < 39; d++)
        text += "error(0.01) D" + std::to_string(d) + " D" + std::to_string(d + 1) + "\n";
    const stim::DetectorErrorModel dem(text);
    pp::AffinityTable table(dem);
    for (size_t a = 0; a < 40; a++)
        for (size_t b = 0; b < 40; b++)
        {
            const size_t distance = a > b ? a - b : b - a;
            require(table(a, b) == std::min(distance, size_t(31)),
                    "packed affinity distance or symmetry is incorrect");
        }
    require(table(0, 40) == 1 && table(19, 40) == 20 && table(39, 40) == 31,
            "boundary distances are incorrect");

    // Exercise matching removal and observable accumulation on identical repeated calls.
    pp::PDZ pdz(dem, table, {.aff_keep_threshold=2});
    dec::Result prev;
    dec::ObsType flip(1);
    flip[0] = true;
    prev.matching_data.add({.d1=0, .d2=1, .pr=1, .w_qu=0, .frame_flips=flip});
    prev.matching_data.add({.d1=18, .d2=19, .pr=1, .w_qu=0, .frame_flips=dec::ObsType(1)});
    for (size_t repeat = 0; repeat < 3; repeat++)
    {
        dec::SyndromeType syndrome(40);
        syndrome[0] = syndrome[1] = syndrome[18] = syndrome[19] = true;
        const auto out = pdz.apply(syndrome, prev);
        require(syndrome.popcnt() == 0 && out.matching_data.assignments.size() == 2,
                "accepted matches were not removed consistently");
        require(out.flipped_obs[0] && out.matching_data.frame_flips[0],
                "accepted match observable correction was lost");
    }

    // A negative margin leaves the matching available to the next decoder.
    dec::Result distant_prev;
    distant_prev.matching_data.add({.d1=0, .d2=8, .pr=1, .w_qu=0, .frame_flips=flip});
    distant_prev.matching_data.add({.d1=1, .d2=2, .pr=1, .w_qu=0, .frame_flips=dec::ObsType(1)});
    dec::SyndromeType distant_syndrome(40);
    distant_syndrome[0] = distant_syndrome[1] = distant_syndrome[2] = distant_syndrome[8] = true;
    const auto distant_out = pdz.apply(distant_syndrome, distant_prev);
    require(distant_syndrome.popcnt() == 4 && distant_out.matching_data.assignments.empty() &&
            !distant_out.flipped_obs[0], "negative margin did not preserve the residual syndrome");

    // A distant pair must not hide a nearby competing event.
    dec::Result crowded_prev = prev;
    crowded_prev.matching_data.add({.d1=2, .d2=3, .pr=1, .w_qu=0, .frame_flips=dec::ObsType(1)});
    dec::SyndromeType crowded_syndrome(40);
    for (size_t d : {0, 1, 2, 3, 18, 19})
        crowded_syndrome[d] = true;
    const auto crowded_out = pdz.apply(crowded_syndrome, crowded_prev);
    require(crowded_syndrome.popcnt() == 4 && crowded_syndrome[0] && crowded_syndrome[1] &&
            crowded_out.matching_data.assignments.size() == 1 && !crowded_out.flipped_obs[0],
            "a distant competitor hid a nearby competing event");

    // The infinity sentinel explicitly accepts a pair with no other fired events.
    dec::Result isolated_prev;
    isolated_prev.matching_data.add({.d1=0, .d2=8, .pr=1, .w_qu=0, .frame_flips=flip});
    dec::SyndromeType isolated_syndrome(40);
    isolated_syndrome[0] = isolated_syndrome[8] = true;
    const auto isolated_out = pdz.apply(isolated_syndrome, isolated_prev);
    require(isolated_syndrome.popcnt() == 0 && isolated_out.flipped_obs[0],
            "isolated pair did not pass the no-competitor margin");

    // Multiple assignments may share the boundary, including reversed endpoint order.
    const stim::DetectorErrorModel boundary_dem("error(0.1) D0 L0\nerror(0.1) D1\n");
    pp::PDZ boundary_pdz(boundary_dem, {.aff_keep_threshold=0});
    dec::Result boundary_prev;
    boundary_prev.matching_data.add({.d1=-1, .d2=0, .pr=1, .w_qu=0, .frame_flips=flip});
    boundary_prev.matching_data.add({.d1=1, .d2=-1, .pr=1, .w_qu=0, .frame_flips=dec::ObsType(1)});
    dec::SyndromeType boundary_syndrome(2);
    boundary_syndrome[0] = boundary_syndrome[1] = true;
    const auto boundary_out = boundary_pdz.apply(boundary_syndrome, boundary_prev);
    require(boundary_syndrome.popcnt() == 0 && boundary_out.flipped_obs[0],
            "shared boundary assignments were not handled correctly");

    bool rejected = false;
    try
    {
        dec::SyndromeType invalid(2);
        invalid[0] = true;
        boundary_pdz.apply(invalid, boundary_prev);
    }
    catch (const std::invalid_argument&) { rejected = true; }
    require(rejected, "matching containing an unfired endpoint was accepted");
    boundary_syndrome[0] = boundary_syndrome[1] = true;
    require(boundary_pdz.apply(boundary_syndrome, boundary_prev).flipped_obs[0],
            "failed validation contaminated the next call");

    // Validate GZIP serialization, the file constructor and the delegating constructor.
    const std::string path = "pdz_affinity_test.gz";
    gzFile file = gzopen(path.c_str(), "wb");
    require(file != nullptr, "could not create affinity test file");
    table.write_to_file(file);
    require(gzclose(file) == Z_OK, "could not close affinity test file");
    file = gzopen(path.c_str(), "rb");
    require(file != nullptr, "could not reopen affinity test file");
    pp::AffinityTable loaded(file);
    gzclose(file);
    for (size_t a = 0; a <= 40; a++)
        for (size_t b = 0; b <= 40; b++)
            require(loaded(a, b) == table(a, b), "affinity file round trip changed entries");
    file = gzopen(path.c_str(), "rb");
    pp::PDZ from_file(dem, file);
    gzclose(file);
    pp::PDZ from_dem(dem);

    // Construction and application to real SI1000 opposing-basis detector output.
    for (bool memory_x : {false, true})
    {
        const auto circuit = sc_si1000(3, 3, 0.001, memory_x, true);
        const auto model = stim::circuit_to_dem(circuit, {.decompose_errors=true});
        pp::PDZ predecoder(model, {.aff_keep_threshold=0});
        dec::PyMatching matching(model, {.report_matching_data=true});
        dec::SyndromeType syndrome(model.count_detectors());
        dec::ObsType truth(model.count_observables());
        syndrome[0] = true;
        const auto matched = matching.decode(syndrome, truth);
        const auto partial = predecoder.apply(syndrome, matched);
        require(syndrome.popcnt() == 0 && partial.flipped_obs == matched.flipped_obs,
                "SI1000 partial correction differs from the accepted matching");
    }
    std::cout << "PDZ and affinity tests passed\n";
}
