/*
 *  contributor: OpenAI Codex
 * */

#include "decoder/pymatching.h"

#include <iostream>
#include <stdexcept>

int main()
{
    const stim::DetectorErrorModel dem(
            "error(0.02) D0 L0\nerror(0.11) D1\nerror(0.23) D2 L0\n"
            "error(0.005) D0 D1\nerror(0.4) D0\nerror(0.4) D1\n");
    dec::PyMatching ordinary(dem);
    dec::PyMatching reporting(dem, {.report_matching_data=true});

    dec::SyndromeType syndrome(dem.count_detectors());
    dec::ObsType truth(dem.count_observables());
    for (size_t shot = 0; shot < 4; shot++)
    {
        syndrome[0] = (shot & 1);
        syndrome[1] = ((shot >> 1) & 1);
        syndrome[2] = ((shot + 1) & 1);
        const auto expected = ordinary.decode(syndrome, truth);
        const auto actual = reporting.decode(syndrome, truth);
        if (actual.flipped_obs != expected.flipped_obs)
            throw std::runtime_error("matched-edge observable parity differs from ordinary PyMatching");
        dec::ObsType parity(dem.count_observables());
        for (const auto& edge : actual.matching_data.assignments)
        {
            if (edge.d1 < -1 || edge.d1 >= static_cast<int64_t>(dem.count_detectors()) ||
                edge.d2 < -1 || edge.d2 >= static_cast<int64_t>(dem.count_detectors()))
                throw std::runtime_error("reported matching endpoint is out of range");
            parity ^= edge.frame_flips;
        }
        if (parity != actual.flipped_obs)
            throw std::runtime_error("reported pairwise observable masks do not match prediction");
        for (size_t d = 0; d < dem.count_detectors(); d++)
        {
            size_t degree = 0;
            for (const auto& edge : actual.matching_data.assignments)
                degree += (edge.d1 == static_cast<int64_t>(d)) + (edge.d2 == static_cast<int64_t>(d));
            if (degree != static_cast<size_t>(syndrome[d]))
                throw std::runtime_error("reported matching endpoints do not cover fired detectors exactly once");
        }
    }

    const stim::DetectorErrorModel pair_dem(
            "error(0.2) D0 D1\nerror(0.01) D0\nerror(0.01) D1\n");
    dec::PyMatching pair_decoder(pair_dem, {.report_matching_data=true});
    dec::SyndromeType pair_syndrome(2);
    dec::ObsType pair_truth(0);
    pair_syndrome[0] = pair_syndrome[1] = 1;
    const auto pair = pair_decoder.decode(pair_syndrome, pair_truth);
    bool found_pair = false;
    for (const auto& edge : pair.matching_data.assignments)
        found_pair |= (edge.d1 == 0 && edge.d2 == 1) || (edge.d1 == 1 && edge.d2 == 0);
    if (!found_pair)
        throw std::runtime_error("detector-to-detector match was not reported as a pair");

    const stim::DetectorErrorModel boundary_dem("error(0.1) D0 L0\n");
    dec::PyMatching boundary_decoder(boundary_dem, {.report_matching_data=true});
    dec::SyndromeType boundary_syndrome(1);
    dec::ObsType boundary_truth(1);
    boundary_syndrome[0] = 1;
    const auto boundary = boundary_decoder.decode(boundary_syndrome, boundary_truth);
    if (boundary.matching_data.assignments.size() != 1 ||
        boundary.matching_data.assignments[0].d1 != 0 ||
        boundary.matching_data.assignments[0].d2 != -1)
        throw std::runtime_error("boundary match was not reported with the virtual boundary endpoint");

    const stim::DetectorErrorModel gap_dem("error(0.1) D0 L0\nerror(0.02) D0\n");
    dec::PyMatching gap_ordinary(gap_dem, {.enable_gap_estimation=true});
    dec::PyMatching gap_reporting(gap_dem,
            {.enable_gap_estimation=true, .report_matching_data=true});
    for (size_t fired = 0; fired < 2; fired++)
    {
        boundary_syndrome[0] = fired;
        const auto gap_expected = gap_ordinary.decode(boundary_syndrome, boundary_truth);
        const auto gap_actual = gap_reporting.decode(boundary_syndrome, boundary_truth);
        if (gap_actual.flipped_obs != gap_expected.flipped_obs ||
            gap_actual.matching_data.frame_flips != gap_actual.flipped_obs)
            throw std::runtime_error("gap plus endpoint reporting changed the selected observable parity");
    }
    std::cout << "PyMatching matched-edge reporting passed\n";
}
