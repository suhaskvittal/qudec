/*
 *  contributor: OpenAI Codex
 * */

#include "decoder/pdz_libra.h"

#include <cmath>
#include <iostream>
#include <sstream>
#include <stdexcept>

int
main(int argc, char** argv)
{
    // Two isolated memory-basis pairs are accepted. The four other fired
    // detectors have close competitors, so PDZ rejects their matching pairs.
    const stim::DetectorErrorModel dem(
            "error(0.2) D0 D1\nerror(0.02) D0 L0\nerror(0.02) D1\n"
            "error(0.2) D6 D7\nerror(0.02) D6 L0\nerror(0.02) D7\n"
            // Connect the memory pairs through unfired detectors, satisfying
            // Libra's single memory-component requirement without making the
            // competing fired endpoints close enough for PDZ to reject them.
            "error(0.001) D1 D8\nerror(0.001) D8 D9\n"
            "error(0.001) D9 D10\nerror(0.001) D10 D6\n"
            "error(0.2) D2 D3\nerror(0.2) D4 D5\n"
            "error(0.001) D2 D4\nerror(0.001) D3 D5\n");
    dec::SyndromeType syndrome(dem.count_detectors());
    for (size_t d = 0; d < 8; d++)
        syndrome[d] = true;
    const auto original = syndrome;
    dec::ObsType truth(1);
    dec::PyMatching matching(dem, {.report_matching_data=true});
    const auto expected = matching.decode(syndrome, truth);
    dec::Libra::Config libra;
    libra.ensemble_size = 2;

    dec::PDZLibra at_threshold(dem, {}, libra);
    if (argc == 2)
        at_threshold.record_syndromes(argv[1], 0);
    const auto result = at_threshold.decode(syndrome, truth);
    at_threshold.finish_recording();
    const auto stats = at_threshold.statistics();
    if (stats.input_events != 8 || stats.residual_events != 4 || stats.bypassed_shots != 1)
        throw std::runtime_error("Exactly 50% removal did not bypass Libra at the default threshold");
    if (result.flipped_obs != expected.flipped_obs ||
        result.matching_data.frame_flips != expected.matching_data.frame_flips ||
        result.matching_data.assignments.size() != expected.matching_data.assignments.size())
        throw std::runtime_error("Bypass did not retain the original PyMatching result");
    if (syndrome != original)
        throw std::runtime_error("Bypass changed the input syndrome");
    std::ostringstream output;
    at_threshold.print_stats(output);
    if (output.str().find("Libra shots=0 ensemble_shots=0") == std::string::npos)
        throw std::runtime_error("Bypass called Libra");

    dec::PDZLibra above_threshold(dem, {}, libra, std::nextafter(0.5, 1.0));
    above_threshold.decode(syndrome, truth);
    if (above_threshold.statistics().bypassed_shots)
        throw std::runtime_error("Removal below threshold bypassed Libra");
    output.str("");
    above_threshold.print_stats(output);
    if (output.str().find("Libra shots=1") == std::string::npos)
        throw std::runtime_error("Below-threshold shot did not call Libra");

    dec::PDZLibra disabled(dem, {}, libra, 1);
    disabled.decode(syndrome, truth);
    if (disabled.statistics().bypassed_shots)
        throw std::runtime_error("Threshold 1 bypassed a nonempty residual syndrome");

    syndrome.clear();
    at_threshold.decode(syndrome, truth);
    if (at_threshold.statistics().bypassed_shots != 2)
        throw std::runtime_error("Empty syndrome did not bypass Libra");
    syndrome[0] = syndrome[1] = true;
    disabled.decode(syndrome, truth);
    if (disabled.statistics().bypassed_shots != 1)
        throw std::runtime_error("Threshold 1 did not bypass a fully removed syndrome");
    std::cout << "PDZLibra bypass tests passed\n";
}
